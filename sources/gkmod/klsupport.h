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
#include "blocks.h"
#include "descents.h"

namespace atlas {

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace klsupport {

class KLSupport
{
  enum State { DownsetsFilled, LengthLessFilled, Filled, NumStates};

  bitset::BitSet<NumStates> d_state;

  blocks::Block_base& d_block;  // non-owned reference

  std::vector<bitset::RankFlags> d_descent;
  std::vector<bitset::RankFlags> d_goodAscent;
  std::vector<bitmap::BitMap> d_downset;
  std::vector<bitmap::BitMap> d_primset;
  std::vector<blocks::BlockElt> d_lengthLess;

 public:

// constructors and destructors
  KLSupport(blocks::Block_base&);

// copy and swap (use automatically generated copy constructor)
  void swap(KLSupport&);

// accessors

  const blocks::Block_base& block() const { return d_block; }
  blocks::BlockElt cross(size_t s, blocks::BlockElt z) const
    { return d_block.cross(s,z); }
  blocks::BlockEltPair cayley(size_t s, blocks::BlockElt z) const
    { return d_block.cayley(s,z); }
  const bitset::RankFlags& descentSet(blocks::BlockElt z) const
    { return d_descent[z]; }
  /*!
\brief Descent status of simple root s for block element z. Taken directly from the block.
  */
  descents::DescentStatus::Value descentValue(size_t s, blocks::BlockElt z)
    const
    { return d_block.descentValue(s,z); }
  const descents::DescentStatus& descent(blocks::BlockElt y) const // full info
    { return d_block.descent(y); }

  size_t rank() const { return d_block.rank(); }
  size_t size() const { return d_block.size(); }

  const bitset::RankFlags& goodAscentSet(blocks::BlockElt z) const
    { return d_goodAscent[z]; }
  size_t length(blocks::BlockElt z) const { return d_block.length(z); }
  /*!
\brief Number of block elements of length strictly less than l.
  */
  blocks::BlockElt lengthLess(size_t l) const { return d_lengthLess[l]; }

  blocks::BlockElt primitivize
    (blocks::BlockElt x, const bitset::RankFlags& A) const;

  // the following are filters of the bitmap
  void extremalize(bitmap::BitMap&, const bitset::RankFlags&) const;
  void primitivize(bitmap::BitMap&, const bitset::RankFlags&) const;

// manipulators
  void fill();
  void fillDownsets();

};

} // namespace klsupport

} // namespace atlas

#endif
