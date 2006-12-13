/*!
\file
\brief Class definition and function declarations for KLSupport.
*/
/*
  This is klsupport.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef KLSUPPORT_H  /* guard against multiple inclusions */
#define KLSUPPORT_H

#include <map>

#include "klsupport_fwd.h"

#include "bitmap.h"
#include "bitset.h"
#include "blocks.h"
#include "descents.h"

namespace atlas {

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace klsupport {

class KLSupport_new {
};

class KLSupport {

 private:

  enum State {PrimitivizeFilled, DownsetsFilled, LengthLessFilled, 
	      Filled, NumStates};

  bitset::BitSet<NumStates> d_state;

  blocks::Block* d_block;  // non-owned pointer
  size_t d_rank;
  size_t d_size;

  std::vector<blocks::BlockEltList> d_primitivize;
  std::vector<bitset::RankFlags> d_descent;
  std::vector<bitset::RankFlags> d_goodAscent;
  std::vector<bitmap::BitMap> d_downset;
  std::vector<bitmap::BitMap> d_primset;
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

  size_t cross(size_t s, size_t z) const {
    return d_block->cross(s,z);
  }

  blocks::BlockEltPair cayley(size_t s, size_t z) const {
    return d_block->cayley(s,z);
  }

  const bitset::RankFlags& descentSet(size_t z) const {
    return d_descent[z];
  }

  descents::DescentStatus::Value descentValue(size_t s, size_t z) const {
    return d_block->descentValue(s,z);
  }

  void extremalize(bitmap::BitMap&, const bitset::RankFlags&) const;

  const bitset::RankFlags& goodAscentSet(size_t z) const {
    return d_goodAscent[z];
  }

  size_t length(size_t z) const {
    return d_block->length(z);
  }

  size_t lengthLess(size_t l) const {
    return d_lengthLess[l];
  }

  void primitivize(bitmap::BitMap&, const bitset::RankFlags&) const;

  bool primitivize(size_t& x, const bitset::RankFlags& A) const {
    if (d_primitivize[A.to_ulong()][x] == blocks::UndefBlock) return false; 
    x = d_primitivize[A.to_ulong()][x]; 
    return true;
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

  size_t numExtremals();
};

}

}

#endif
