/*
  This is klsupport.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/* Implementation for KLSupport.

  This module provides support code for the Kazhdan-Lusztig computation,
  mostly the management of the list of primitive pairs, and primitivization
  of arbitrary subsets of the block.
*/

#include <cassert>
#include "klsupport.h"

/*
  After some hesitation, I think I _am_ going to assume that the block has
  been sorted by length and root datum involution, and then by R-packet.
  This does make the hard case in the recursion a lot simpler to handle:
  "thickets" of representations are in fact R-packets, so they will be
  consecutively numbered. [Fokko]
*/

namespace atlas {

/*****************************************************************************

        Chapter I -- The KLSupport class

 *****************************************************************************/

namespace klsupport {

KLSupport::KLSupport(const Block_base& b)
  : d_block(b)
  , d_lengthLess()
  , d_descent(size())
  , d_goodAscent(size())
  , d_downset(rank())
  , d_primset(rank())
{
/*
  Make |d_lengthLess| into a vector of size |max(lengths(d_block))+2| such that
  for |0<=l<=max(lengths(d_block))+1|, the |BlockELt| |d_lengthLess[l]| is the
  first one of length at least |l| in |d_block| (or |d_block.size()| if there
  are none, as is the case for $l=\max(lengths(d_block))+1)$. In other words,
  |d_lengthLess[l]| counts the elements in |d_block| of length less than |l|.
*/
  {
    auto l_size = (d_block.size()==0 ? 0 : d_block.length(d_block.size()-1))+2;
    d_lengthLess.resize(l_size);

    d_lengthLess[0]=0; // no elements of length<0
    size_t l=0;
    for (BlockElt z=0; z<d_block.size(); ++z) // invariant |d_block.length(z)>=l|
      while (d_block.length(z)>l)
	d_lengthLess[++l]=z;

    // do not forget the last length!
    d_lengthLess[l+1]=d_block.size(); // here $l=\max(lengths(d_block))$
  }

/*
  Fill in the |downset|, |primset|, |descents| and |goodAscent| bitmap/set
  vectors. Here |downset| and |primset| are vectors indexed by a simple
  reflection |s|, and giving a bitmap over all block elements, while
  |descents| and |goodAscent| are vectors indexed by a block element |z| and
  giving a bitset over all simple reflections. This difference is motivated by
  their use: |downset| and |primset| are used to filter bitmaps over the
  entire block according to some set of simple generators, which is easier if
  the data is grouped by generator. In fact the data computed is stored twice:
  one always has |downset[s].isMember(z) == descents[z].test(s)| and
  |primset[s].isMember(z) != good_ascent[z].test(s)|

  The predicate that |s| is a |descent| for |z| is taken in the weak sense
  that |s| belongs to the "tau-invariant" of |z|, in other words, it is a
  complex descent, real parity (type I or type II), or imaginary compact (the
  final case does not actually allow going down). The |goodAscent| bitset for
  |z| holds the non-decents for |z| that are not imaginary type II, so they
  are either complex ascent, imaginary type I or real nonparity. The |primset|
  bitmap for |s| records the block elements |z| for which |s| is not a
  |goodAscent|, in other words it is either a |descent| or imaginary type II.
*/
  size_t size = d_block.size();

  for (weyl::Generator s=0; s<rank(); ++s)
  {
    d_downset[s].set_capacity(size);
    d_primset[s].set_capacity(size);
    for (BlockElt z = 0; z < size; ++z)
    {
      DescentStatus::Value v = descentValue(s,z);
      if (DescentStatus::isDescent(v))
      {
	d_downset[s].insert(z);
	d_primset[s].insert(z);
	d_descent[z].set(s);
      }
      else // ascents
	if (v == DescentStatus::ImaginaryTypeII)
	  d_primset[s].insert(z);  // s is a "bad" ascent
	else
	  d_goodAscent[z].set(s); // good ascent
    }
  }
} // |KLSupport::KLSupport|

/******** accessors **********************************************************/

/*
  Flag in |b|, which is of size |size()|, those block elements which are
  extremal w.r.t. the simple reflections in |d|, i.e., for which all simple
  generators flagged in |d| are descents. Since |d_downset[s]| flags the
  elements for which |s| is a descent, this amounts to requesting that |z|
  belong to the intersection of all the downsets of generators flagged in |d|.
*/

void KLSupport::filter_extremal(BitMap& b, const RankFlags& d)
  const
{
  for (weyl::Generator s=0; s<rank(); ++s)
    if (d.test(s))
      b &= d_downset[s];
}


/*
  Flag in |b|, which is of size |size()|, those block elements which are
  extremal w.r.t. the simple reflections in |d|, i.e., for which all simple
  generators flagged in |d| are either descents or imaginary type II ascents.
  Since |d_primset[s]| flags the elements for which |s| is a descent or
  imaginary type II ascent, this amounts to requesting that |z| belong to the
  intersection of all the primsets of generators flagged in |d|.
*/
void KLSupport::filter_primitive(BitMap& b, const RankFlags& d) const
{
  for (weyl::Generator s=0; s<rank(); ++s)
    if (d.test(s))
      b &= d_primset[s];
}


/*
  Find for |x| a primitive element for |d| above it, returning that value, or
  return |d_block.size()| if a real nonparity case is hit, or if (in partial
  blocks) ascent through a complex ascent or Cayley transform is attempted,
  which link points outside of the block (it is represented as |UnderBlock|).

  A primitive element for |d| is one for which all elements in |d| are either
  descents or type II imaginary ascents. So if |x| is not primitive, it has an
  ascent in |d| that is either complex, imaginary type I or real nonparity. In
  the first two cases we replace |x| by the (unique) ascended element and
  continue; in the last case, we return |UndefBlock| (for K-L computations, this
  case implies that $P_{x,y}=0$; the value of |UndefBlock| is conveniently
  larger than any valid BlockElt |y|, so this case will be handled effortlessly
  together with triangularity). It is also permissible to pass |x==UndefBlock|,
  for which will |d_block.size()| be returned immediately.
*/
BlockElt
  KLSupport::primitivize(BlockElt x, const RankFlags& d) const
{
  while (x!=UndefBlock)
  {
    RankFlags a = // good ascents for |x| that are descents for |y|
      goodAscentSet(x) & d;
    if (a.none()) // then we have succeeded in making |x| primitive
      return x;
    weyl::Generator s = a.firstBit();
    DescentStatus::Value v = descentValue(s,x);
    x = v == DescentStatus::RealNonparity ? UndefBlock
      : v == DescentStatus::ComplexAscent ? d_block.cross(s,x)
      : d_block.cayley(s,x).first; // imaginary type I
  }
  return d_block.size();
}

#ifndef NDEBUG
void KLSupport::check_sub(const KLSupport& sub, const BlockEltList& embed)
{
  assert(sub.rank()==rank());
  for (unsigned i=1; i<embed.size(); ++i)
    assert(embed[i-1]<embed[i]);
  for (BlockElt z=0; z<sub.block().size(); ++z)
  {
    assert(sub.block().length(z)==d_block.length(embed[z]));
    assert(sub.block().descent(z)==d_block.descent(embed[z]));
    for (weyl::Generator s=0; s<d_block.rank(); ++s)
      if (sub.cross(s,z)!=UndefBlock)
	assert(embed[sub.cross(s,z)]==cross(s,embed[z]));
  }
}
#endif

} // |namespace klsupport|

} // |namespace atlas|
