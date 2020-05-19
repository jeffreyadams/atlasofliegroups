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

  struct Elt_info // per block element information
  { RankFlags descents;
    RankFlags good_ascents;
  Elt_info(RankFlags d,RankFlags a): descents(d), good_ascents(a) {}
  };

  std::vector<Elt_info> info;
  std::vector<BlockElt> length_stop; // |length_stop[l]| is first of length |l|

 public:

// constructors and destructors
  KLSupport(const Block_base&);

// accessors
  const Block_base& block () const { return d_block; }
  size_t rank () const { return d_block.rank(); }
  size_t size () const { return d_block.size(); } // also |info.size()|

  size_t length (BlockElt z) const { return d_block.length(z); }
  BlockElt length_less (size_t l) const // number of block elements of length<l
  { return length_stop[l]; }

  BlockElt cross (size_t s, BlockElt z) const  { return d_block.cross(s,z); }
  BlockEltPair cayley (size_t s, BlockElt z) const
    { return d_block.cayley(s,z); }

  DescentStatus::Value descent_value (size_t s, BlockElt z) const
    { return d_block.descentValue(s,z); }
  const DescentStatus& descent(BlockElt y) const // combined for all |s|
    { return d_block.descent(y); }

  RankFlags descent_set (BlockElt z) const { return info[z].descents; }
  RankFlags good_ascent_set (BlockElt z) const { return info[z].good_ascents; }

  // find ascent for |x| that is descent for |y| if any; |longBits| if none
  unsigned int ascent_descent (BlockElt x,BlockElt y) const
    { return (descent_set(y)-descent_set(x)).firstBit(); }

  // find non-i2 ascent for |x| that is descent for |y| if any; or |longBits|
  unsigned int good_ascent_descent (BlockElt x,BlockElt y) const
    { return (good_ascent_set(x)&descent_set(y)).firstBit(); }

  bool is_extremal (BlockElt x, RankFlags descents_y) const
    { return descent_set(x).contains(descents_y); }
  bool is_primitive (BlockElt x, RankFlags descents_y) const
    { return (good_ascent_set(x) & descents_y).none(); } // disjointness
  bool is_extremal (BlockElt x, BlockElt y) const
    { return is_extremal(x,descent_set(y)); }
  bool is_primitive (BlockElt x, BlockElt y) const
    { return is_primitive(x,descent_set(y)); }

  // in practice the main operation for extremals/primitives is reverse traversal

  // set $x$ to last extremal element for $y$ strictly before $x$, or fail
  bool extr_back_up(BlockElt& x, RankFlags desc_y) const
    { while (x-->0) if (is_extremal(x,desc_y)) return true;
      return false; // now |x| has crashed through 0 and should be ignored
    }

  // set $x$ to last primitive element for $y$ strictly before $x$, or fail
  bool prim_back_up(BlockElt& x, RankFlags desc_y) const
    { while (x-->0) if (is_primitive(x,desc_y)) return true;
      return false; // now |x| has crashed through 0 and should be ignored
    }

  // primitive element for |desc_y| above and reachable from |x|, or block size
  BlockElt primitivize (BlockElt x, RankFlags desc_y) const;

#ifndef NDEBUG
  void check_sub(const KLSupport& sub, const BlockEltList& embed);
#endif
};

} // |namespace klsupport|

} // |namespace atlas|

#endif
