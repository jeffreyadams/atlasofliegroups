/*
  This is blocks.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef BLOCKS_H  /* guard against multiple inclusions */
#define BLOCKS_H

#include "blocks_fwd.h"

#include "bruhat_fwd.h"
#include "descents.h"
#include "kgb_fwd.h"
#include "complexredgp_fwd.h"
#include "weyl_fwd.h"

#include "bitset.h"
#include "realform.h"
#include "weyl.h"

namespace atlas {

/******** type declarations *************************************************/

/******** constant declarations *********************************************/

namespace blocks {

const BlockElt UndefBlock = ~0ul;

}

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace blocks {

class Block {

 protected:

  enum State { BruhatConstructed, NumStates };

  size_t d_rank;
  size_t d_xsize;
  size_t d_ysize;

  kgb::KGBEltList d_x;
  kgb::KGBEltList d_y;
  std::vector<BlockEltList> d_cross;
  std::vector<BlockEltPairList> d_cayley;
  std::vector<BlockEltPairList> d_inverseCayley;
  descents::DescentStatusList d_descent;
  std::vector<size_t> d_length;
  weyl::WeylEltList d_tau;

  realform::RealForm d_realForm;
  realform::RealForm d_dualForm;

  bitset::BitSet<NumStates> d_state;

  bruhat::BruhatOrder* d_bruhat;

  const weyl::WeylGroup* d_weylGroup;

 public:

// constructors and destructors
  Block();

  Block(complexredgp::ComplexReductiveGroup&, realform::RealForm,
	realform::RealForm);

  virtual ~Block();

// copy, assignment and swap
  void swap(Block&);

// accessors
  const bruhat::BruhatOrder& bruhatOrder() const {
    return *d_bruhat;
  }

  BlockEltPair cayley(size_t s, BlockElt z) const {
    return d_cayley[s][z];
  }

  BlockElt cross(size_t s, BlockElt z) const {
    return d_cross[s][z];
  }

  const descents::DescentStatus& descent(BlockElt z) const {
    return d_descent[z];
  }

  descents::DescentStatus::Value descentValue(size_t s, BlockElt z) const {
    return d_descent[z][s];
  }

  realform::RealForm dualForm() const {
    return d_dualForm;
  }

  BlockEltPair inverseCayley(size_t s, BlockElt z) const {
    return d_inverseCayley[s][z];
  }

  bool isStrictAscent(size_t, BlockElt) const;

  bool isStrictDescent(size_t, BlockElt) const;

  size_t length(BlockElt z) const {
    return d_length[z];
  }

  size_t rank() const {
    return d_rank;
  }

  realform::RealForm realForm() const {
    return d_realForm;
  }

  size_t size() const {
    return d_tau.size();
  }

  const weyl::WeylElt& tau(BlockElt z) const;

  const weyl::WeylGroup& weylGroup() const {
    return *d_weylGroup;
  }

  kgb::KGBElt x(BlockElt z) const {
    return d_x[z];
  }

  size_t xsize() const {
    return d_xsize;
  }

  kgb::KGBElt y(BlockElt z) const {
    return d_y[z];
  }

  size_t ysize() const {
    return d_ysize;
  }

  // manipulators
  void fillBruhat();

};

}

}

#endif
