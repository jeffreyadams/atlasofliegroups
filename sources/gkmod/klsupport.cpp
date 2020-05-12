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
  , d_prim_index(1ul << rank()) // $2^r$ empty slots, with $r$ (semisimple) rank
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


/*
  Fill in the table |d_prim_index|.

  For each |z| in the block, and each descent set |A|, one can from |z| either
  reach a primitive element |z'| for |A| by following 'C+' or 'i1' ascents in
  |A| , or run into an 'rn' ascent in |A| which aborts the search (because when
  this happens, one has $P_{z,y}=0$ for any |y| for which |descentSet(y)==A|).

  After doing this for any |A| and |z|, we store not |z'|, but its index in the
  list of primitives for |A|, recoreded as |prim_index[z]| below. This avoids
  having to use binary search to locate the polynomial $P_{z',y}$, provided it
  is stored at an offset |prim_index[z]| in the appropriate vector. This does
  mean we will have to store null polynomials at primitive elements to ensure
  everything is at its predicted place, while such (fairly common) polynomials
  could be suppressed when using pairs $(x,P_{x,y})$ and binary search on $x$.
*/

void KLSupport::fill_prim_index(prim_index_tp& index_vec,RankFlags A) const
{
  index_vec.resize(d_block.size()); // create slots; we will fill backwards

  BitMap primitives(size()); primitives.fill();
  filter_primitive(primitives,A); // compute all primitive elements for A

  unsigned int prim_count = primitives.size(); // start at high end
  constexpr unsigned int dead_end = -1; // signal no valid primitivization

  for (BlockElt z = d_block.size(); z-->0;)
  {
    auto& dest = index_vec[z]; // the slot to fill during this iteration
    RankFlags a = goodAscentSet(z)&A;
    if (a.none())
    { // then |z| is primitive, record its index
      assert(primitives.isMember(z)); // sanity check
      dest = --prim_count;
    }
    else
    {
      weyl::Generator s = a.firstBit();
      switch (descentValue(s,z))
      {
      case DescentStatus::RealNonparity: dest = dead_end; break;
      case DescentStatus::ComplexAscent:
	{ auto sz=d_block.cross(s,z);
	  dest = sz==UndefBlock ? dead_end : index_vec[sz];
	} break;
      case DescentStatus::ImaginaryTypeI:
	{ auto sz=d_block.cayley(s,z).first;
	  dest = sz==UndefBlock ? dead_end : index_vec[sz];
	} break;
      default: assert(false);
      }
    }
  } // |for(z-->0)|
  assert(prim_count==0); // we've seen all primitives
} // |fill_prim_index|

/******** accessors **********************************************************/

/*
  Flag in |b|, which is of size |size()|, those block elements which are
  extremal w.r.t. the simple reflections in |d|, i.e., for which all simple
  generators flagged in |d| are descents. Since |d_downset[s]| flags the
  elements for which |s| is a descent, this amounts to requesting that |z|
  belong to the intersection of all the downsets of generators flagged in |d|.
*/

void KLSupport::filter_extremal(BitMap& b, const RankFlags& d) const
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


#if 0 // code disabled because replaced by table look-up
/*
  Find for |x| a primitive element for |d| above it, returning that value, or
  return |UndefBlock| if a real nonparity case is hit, or (in partial blocks)
  ascent through an undefined complex ascent or Cayley transform is attempted

  A primitive element for |d| is one for which all elements in |d| are either
  descents or type II imaginary ascents. So if |x| is not primitive, it has an
  ascent in |d| that is either complex, imaginary type I or real nonparity. In
  the first two cases we replace |x| by the (unique) ascended element and
  continue; in the last case, we return |UndefBlock| (for K-L computations, this
  case implies that $P_{x,y}=0$; the value of |UndefBlock| is conveniently
  larger than any valid BlockElt |y|, so this case will be handled effortlessly
  together with triangularity). It is also permissible to pass |x==UndefBlock|,
  which will be returned immediately.
*/
BlockElt
  KLSupport::primitivize(BlockElt x, const RankFlags& d) const
{
  RankFlags a; // good ascents for x that are descents for y

  while (x!=UndefBlock and (a = goodAscentSet(x)&d).any())
  {
    size_t s = a.firstBit();
    DescentStatus::Value v = descentValue(s,x);
    x = v == DescentStatus::RealNonparity ? UndefBlock
      : v == DescentStatus::ComplexAscent ? d_block.cross(s,x)
      : d_block.cayley(s,x).first; // imaginary type I
  }
  return x;
}
#endif // code disabled because replaced by table look-up

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
