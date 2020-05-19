/*
  This is klsupport.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/*
  This module provides support code for the Kazhdan-Lusztig computation,
  mostly the management of the lists of primitive pairs, and primitivization
  of arbitrary subsets of the block.
*/

#include <cassert>
#include "klsupport.h"

/*
  [Original note by Fokko du Cloux (no longer pertinent, thickets are avoided)]
  After some hesitation, I think I _am_ going to assume that the block has
  been sorted by length and root datum involution, and then by R-packet.
  This does make the hard case in the recursion a lot simpler to handle:
  "thickets" of representations are in fact R-packets, so they will be
  consecutively numbered.
*/

namespace atlas {

/*****************************************************************************

				The |KLSupport| class

 *****************************************************************************/

namespace klsupport {

KLSupport::KLSupport(const Block_base& b)
  : d_block(b)
  , info()
  , length_stop()
{
/*
  Make |length_stop| into a vector of size |max(lengths(d_block))+2| such that
  for |0<=l<=max(lengths(d_block))+1|, the |BlockElt| |length_stop[l]| is the
  first one of length at least |l| in |d_block| (or |d_block.size()| if there
  are none, as is the case for $l=1+\max(lengths(d_block))$). In other words,
  |length_stop[l]| counts the elements in |d_block| of length less than |l|.
*/
  {
    length_stop.reserve
      (d_block.size()==0 ? 1 : 2+d_block.length(d_block.size()-1));

    // the following loop could handle length jumps, although that never happens
    for (BlockElt z=0; z<d_block.size(); ++z) // invariant |d_block.length(z)>=l|
      while (length_stop.size()<=d_block.length(z)) // in fact runs at most once
	length_stop.push_back(z);

    // at block size as final |length_stop| (although it appears to be unused)
    length_stop.push_back(d_block.size()); // index is $1+\max(lengths(d_block))$
  }

/*
  Fill in |info|, a vector indexed by a block element |z| and giving two bitsets
  over all simple reflections. The |descents| field of |info[z]| is the
  "tau-invariant" of |z|: bits are set for simple reflections |s| that are
  either a complex descent, real parity (type I or type II), or imaginary
  compact for |z| (the final case does not actually allow descending through
  |s|). The |good_ascents| filed for |z| flags those |s| that are neither
  decents for |z|, nor imaginary type II ascents, so they are either complex
  ascent, imaginary type I or real nonparity.
*/
  info.reserve(d_block.size());
  for (BlockElt z = 0; z < d_block.size(); ++z)
  {
    RankFlags desc, good_asc;
    for (weyl::Generator s=0; s<rank(); ++s)
    {
      DescentStatus::Value v = descent_value(s,z);
      if (DescentStatus::isDescent(v))
	desc.set(s);
      else if (v != DescentStatus::ImaginaryTypeII)
	good_asc.set(s);
    } // |for(s)|
    info.emplace_back(desc,good_asc);
  } // |for(BlockElt z)|
} // |KLSupport::KLSupport|

/******** accessors **********************************************************/

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
  KLSupport::primitivize(BlockElt x, RankFlags desc_y) const
{
  while (x!=UndefBlock)
  {
    RankFlags a = // good ascents for |x| that are descents for |y|
      good_ascent_set(x) & desc_y;
    if (a.none()) // then we have succeeded in making |x| primitive
      return x;
    weyl::Generator s = a.firstBit();
    if (descent_value(s,x) == DescentStatus::RealNonparity)
      break; // and |return d_block.size()|
    x = d_block.unique_ascent(s,x);
  }
  return d_block.size(); // indicate that a dead end was reached
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
