/*
  This is kl.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef KL_H  /* guard against multiple inclusions */
#define KL_H

#include <limits>
#include <set>

#include "kl_fwd.h"

#include "blocks_fwd.h"

#include "bitset.h"
#include "klsupport.h"
#include "polynomials.h"
#include "wgraph.h"

namespace atlas {

/******** constant declarations *********************************************/

namespace kl {

  const KLPol Zero;
  const KLPol One(0);

  const KLCoeff UndefKLCoeff = std::numeric_limits<KLCoeff>::max();
  const KLCoeff UndefMuCoeff = std::numeric_limits<MuCoeff>::max();

}

/******** function declarations *********************************************/

namespace kl {

  void wGraph(wgraph::WGraph&, const KLContext&);

}

/******** type definitions **************************************************/

namespace kl {

class KLContext {

 protected:

  enum State { KLFilled, NumStates };

  bitset::BitSet<NumStates> d_state;

  klsupport::KLSupport* d_support;   // non-owned pointer

  std::vector<klsupport::ExtremalRow> d_extr;

  std::vector<KLRow> d_kl;           // list of polynomial pointers
  std::vector<MuRow> d_mu;           // list of mu-coefficients

  std::set<KLPol> d_store;           // the actual polynomials
  KLPtr d_zero;
  KLPtr d_one;

 public:

// constructors and destructors
  KLContext() {}

  KLContext(klsupport::KLSupport&);

  virtual ~KLContext() {}

// copy, assignment and swap
  KLContext(const KLContext&);

  KLContext& operator= (const KLContext&);

  void swap(KLContext&);

// accessors
  const blocks::Block& block() const {
    return d_support->block();
  }

  const klsupport::ExtremalRow& extremalRow(size_t y) const {
    return d_extr[y];
  }

  const bitset::RankFlags& descentSet(size_t y) const {
    return d_support->descentSet(y);
  }

  bool isZero(const KLPtr p) const {
    return p == d_zero;
  }

  const KLPol& klPol(size_t, size_t) const; // to be defined!

  const KLRow& klRow(size_t y) const {
    return d_kl[y];
  }

  size_t length(size_t y) const {
    return d_support->length(y);
  }

  MuCoeff mu(size_t, size_t) const;

  const MuRow& muRow(size_t y) const {
    return d_mu[y];
  }

  const std::set<KLPol>& polStore() const {
    return d_store;
  }

  const size_t rank() const {
    return d_support->rank();
  }

  const size_t size() const {
    return d_kl.size();
  }

// manipulators
  virtual void fill();

  //  const KLPol& klPol(size_t, size_t);
};

}

}

#endif
