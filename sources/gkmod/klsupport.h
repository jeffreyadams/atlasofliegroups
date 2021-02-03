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

#include <cassert>

#include "../Atlas.h"

#include "bitset.h"	// containment
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

  struct prim_index_tp
  {
    std::vector<unsigned int> index; // from |BlockElt| to index of prim'zed
    unsigned int range; // number of primitive elements for this descen set
  prim_index_tp() : index(), range(-1) {}
  };
  std::vector<prim_index_tp> d_prim_index; // indexed by descent set number

 public:

// constructors and destructors
  KLSupport(const Block_base&);

// accessors
  const Block_base& block () const { return d_block; }
  unsigned short rank () const { return d_block.rank(); }
  BlockElt size () const { return d_block.size(); } // also |info.size()|

  unsigned short length (BlockElt z) const { return d_block.length(z); }
  unsigned short l(BlockElt y,BlockElt x) const // length difference
  { assert(length(x)<=length(y)); return length(y)-length(x); }
  BlockElt length_less (size_t l) const // number of block elements of length<l
  { return length_stop[l]; }
  BlockElt length_floor (BlockElt y) const { return length_stop[length(y)]; }

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

  /* computing KL polynomials used to spend a large amount of time evaluating
     primitivize and then to binary-search for the resulting primitive element.
     Since the primitive condition only depends on the descent set, we speed
     things up by tabulating for each descent set the map from block elements
     to the index of their primitivized counterparts.
  */
  void prepare_prim_index (RankFlags A) // call us before any |prim_index(...,A)|
  { prim_index_tp& record=d_prim_index[A.to_ulong()];
    if (record.range==static_cast<unsigned int>(-1))
      fill_prim_index(A);
    assert(record.range!=static_cast<unsigned int>(-1));
  }

  unsigned int prim_index (BlockElt x, RankFlags descent_set) const
  { const prim_index_tp& record=d_prim_index[descent_set.to_ulong()];
    assert(record.range!=static_cast<unsigned int>(-1));
    return x==UndefBlock ? record.range : record.index[x];
  }

  unsigned int nr_of_primitives (RankFlags descent_set) const
  { const prim_index_tp& record=d_prim_index[descent_set.to_ulong()];
    assert(record.range!=static_cast<unsigned int>(-1));
    return record.range;
  }

  // this is where an element |y| occurs in its "own" primitive row
  unsigned int self_index (BlockElt y) const
  { return prim_index(y,descent_set(y)); }

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

  // number of primitive elements for |descent_set(y)| of length less than |y|
  unsigned int col_size (BlockElt y) const
  { BlockElt x=length_floor(y); const RankFlags desc_y=descent_set(y);
    return prim_back_up(x,desc_y) ? prim_index(x,desc_y)+1 : 0;
  }

  // primitive element for |desc_y| above and reachable from |x|, or block size
  BlockElt primitivize (BlockElt x, RankFlags desc_y) const;

  void fill_prim_index(RankFlags A);

#ifndef NDEBUG
  void check_sub(const KLSupport& sub, const BlockEltList& embed);
#endif
};

} // |namespace klsupport|

} // |namespace atlas|

#endif
