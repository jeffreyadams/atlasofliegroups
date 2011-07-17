/*!
\file
\brief Class definition and function declarations for KLSupport.
*/
/*
  This is klsupport.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef KLSUPPORT_H  /* guard against multiple inclusions */
#define KLSUPPORT_H

#include "klsupport_fwd.h"

#include "bitmap.h"
#include "bitset.h"

#include "atlas_types.h"
#include "blocks.h"
#include "descents.h"

namespace atlas {

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace klsupport {

class KLSupport
{
  enum State { DownsetsFilled, LengthLessFilled, Filled, NumStates};

  BitSet<NumStates> d_state;

  Block_base& d_block;  // non-owned reference

  std::vector<RankFlags> d_descent;
  std::vector<RankFlags> d_goodAscent;
  std::vector<BitMap> d_downset;
  std::vector<BitMap> d_primset;
  std::vector<BlockElt> d_lengthLess;

 public:

// constructors and destructors
  KLSupport(Block_base&);

// copy and swap (use automatically generated copy constructor)
  void swap(KLSupport&);

// accessors

  const Block_base& block() const { return d_block; }
  BlockElt cross(size_t s, BlockElt z) const
    { return d_block.cross(s,z); }
  BlockEltPair cayley(size_t s, BlockElt z) const
    { return d_block.cayley(s,z); }
  const RankFlags& descentSet(BlockElt z) const
    { return d_descent[z]; }
  /*!
\brief Descent status of simple root s for block element z. Taken directly from the block.
  */
  DescentStatus::Value descentValue(size_t s, BlockElt z)
    const
    { return d_block.descentValue(s,z); }
  const DescentStatus& descent(BlockElt y) const // full info
    { return d_block.descent(y); }

  size_t rank() const { return d_block.rank(); }
  size_t size() const { return d_block.size(); }

  const RankFlags& goodAscentSet(BlockElt z) const
    { return d_goodAscent[z]; }
  size_t length(BlockElt z) const { return d_block.length(z); }
  /*!
\brief Number of block elements of length strictly less than l.
  */
  BlockElt lengthLess(size_t l) const { return d_lengthLess[l]; }

  BlockElt primitivize
    (BlockElt x, const RankFlags& A) const;

  // the following are filters of the bitmap
  void extremalize(BitMap&, const RankFlags&) const;
  void primitivize(BitMap&, const RankFlags&) const;

// manipulators
  void fill();
  void fillDownsets();

};

} // namespace klsupport

} // namespace atlas

#endif
