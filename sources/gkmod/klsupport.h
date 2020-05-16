/*
  This is klsupport.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Class definition and function declarations for |KLSupport|.

#ifndef KLSUPPORT_H  /* guard against multiple inclusions */
#define KLSUPPORT_H


#include "bitset.h"	// containment

#include "../Atlas.h"
#include "blocks.h"	// inlining of methods like |cross| and |cayley|

namespace atlas {

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace klsupport {

class KLSupport
{
  const Block_base& d_block;  // non-owned reference

  std::vector<BlockElt> d_lengthLess;
  std::vector<RankFlags> d_descent;
  std::vector<RankFlags> d_goodAscent;
  std::vector<BitMap> d_downset;
  std::vector<BitMap> d_primset;

 public:

// constructors and destructors
  KLSupport(const Block_base&);

// accessors
  const Block_base& block () const { return d_block; }
  size_t rank () const { return d_block.rank(); }
  size_t size () const { return d_block.size(); }

  size_t length (BlockElt z) const { return d_block.length(z); }
  BlockElt length_less (size_t l) const // number of block elements of length < l
    { return d_lengthLess[l]; }

  BlockElt cross (size_t s, BlockElt z) const  { return d_block.cross(s,z); }
  BlockEltPair cayley (size_t s, BlockElt z) const
    { return d_block.cayley(s,z); }

  DescentStatus::Value descent_value (size_t s, BlockElt z) const
    { return d_block.descentValue(s,z); }
  const DescentStatus& descent(BlockElt y) const // combined for all |s|
    { return d_block.descent(y); }

  RankFlags descent_set (BlockElt z) const { return d_descent[z]; }
  RankFlags good_ascent_set (BlockElt z) const { return d_goodAscent[z]; }

  // find ascent for |x| that is descent for |y| if any; |longBits| if none
  unsigned int ascent_descent (BlockElt x,BlockElt y) const
    { return (descent_set(y)-descent_set(x)).firstBit(); }

  BlockElt primitivize (BlockElt x, const RankFlags& A) const;

  // find non-i2 ascent for |x| that is descent for |y| if any; or |longBits|
  unsigned int good_ascent_descent (BlockElt x,BlockElt y) const
    { return (good_ascent_set(x)&descent_set(y)).firstBit(); }

  // the following are filters of the bitmap
  void filter_extremal (BitMap&, const RankFlags&) const;
  void filter_primitive (BitMap&, const RankFlags&) const;

#ifndef NDEBUG
  void check_sub(const KLSupport& sub, const BlockEltList& embed);
#endif
};

} // |namespace klsupport|

} // |namespace atlas|

#endif
