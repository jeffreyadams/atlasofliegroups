/*!
\file
\brief Class definition and function declarations for KLSupport.
*/
/*
  This is klsupport.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
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

  using blocks::BlockElt;

class KLSupport {

 private:

  enum State {PrimitivizeFilled, DownsetsFilled, LengthLessFilled,
	      Filled, NumStates};

  bitset::BitSet<NumStates> d_state;

  blocks::Block* d_block;  // non-owned pointer
  size_t d_rank;
  BlockElt d_size;

  std::vector<blocks::BlockEltList> d_primitivize;
  std::vector<bitset::RankFlags> d_descent;
  std::vector<bitset::RankFlags> d_goodAscent;
  std::vector<bitmap::BitMap> d_downset;
  std::vector<bitmap::BitMap> d_primset;
  std::vector<BlockElt> d_lengthLess;

 public:

// constructors and destructors
  KLSupport():d_block(0) {}

  KLSupport(blocks::Block&);

  ~KLSupport() {}

// assignment, copy and swap
  void swap(KLSupport&);

// accessors

  const blocks::Block& block() const {
    return *d_block;
  }

  BlockElt cross(size_t s, BlockElt z) const {
    return d_block->cross(s,z);
  }

  blocks::BlockEltPair cayley(size_t s, BlockElt z) const {
    return d_block->cayley(s,z);
  }

  const bitset::RankFlags& descentSet(BlockElt z) const {
    return d_descent[z];
  }

  /*!
\brief Descent status of simple root s for block element z.
  */
  descents::DescentStatus::Value descentValue(size_t s, BlockElt z) const {
    return d_block->descentValue(s,z);
  }

  void extremalize(bitmap::BitMap&, const bitset::RankFlags&) const;

  const bitset::RankFlags& goodAscentSet(BlockElt z) const {
    return d_goodAscent[z];
  }

  size_t length(BlockElt z) const {
    return d_block->length(z);
  }

  /*!
\brief Number of block elements of length strictly less than l.
  */
  BlockElt lengthLess(size_t l) const {
    return d_lengthLess[l];
  }

  void primitivize(bitmap::BitMap&, const bitset::RankFlags&) const;

  BlockElt primitivize(BlockElt x, const bitset::RankFlags& A) const {
    return d_primitivize[A.to_ulong()][x];
  }

  size_t rank() const {
    return d_rank;
  }

  size_t size() const {
    return d_size;
  }

// manipulators
  void fill();

  void fillDownsets();

  void fillPrimitivize();
};

}

}

#endif
