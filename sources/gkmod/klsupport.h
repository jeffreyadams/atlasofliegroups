/*
  This is klsupport.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef KLSUPPORT_H  /* guard against multiple inclusions */
#define KLSUPPORT_H

#include "klsupport_fwd.h"

#include "blocks_fwd.h"

#include "bitmap.h"
#include "bitset.h"
#include "descents.h"

namespace atlas {

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace klsupport {

class KLSupport {

 private:

  enum State { DownsetsFilled, LengthLessFilled, Filled, NumStates };

  bitset::BitSet<NumStates> d_state;

  blocks::Block* d_block;  // non-owned pointer
  size_t d_rank;

  std::vector<ExtremalRow> d_extrPairs;
  std::vector<bitmap::BitMap> d_downset;
  std::vector<size_t> d_lengthLess;

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

  descents::DescentStatus::Value descentValue(size_t, size_t) const;

  void extremalize(bitmap::BitMap&, const descents::DescentStatus&) const;

  const ExtremalRow& extremalRow(size_t z) const {
    return d_extrPairs[z];
  }

  size_t length(size_t) const;

  size_t lengthLess(size_t l) const {
    return d_lengthLess[l];
  }

  size_t numExtremals(size_t z) const {
    return d_extrPairs[z].size();
  }

  size_t rank() const {
    return d_rank;
  }

  size_t size() const {
    return d_extrPairs.size();
  }

// manipulators
  void fill();

  void fillDownsets();

  size_t numExtremals();
};

}

}

#endif
