/*
  This is kgb.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef KGB_H  /* guard against multiple inclusions */
#define KGB_H

#include "kgb_fwd.h"

#include "bruhat_fwd.h"
#include "realredgp_fwd.h"
#include "weyl_fwd.h"

#include "gradings.h"

namespace atlas {

/******** function declarations *********************************************/

/******** constant declarations *********************************************/

namespace kgb {

const KGBElt UndefKGB = ~0ul;

}

/******** type definitions **************************************************/

namespace kgb {

class KGB {

 protected:

  enum State { BruhatConstructed, NumStates };

  size_t d_rank;

  std::vector<KGBEltList> d_cross;
  std::vector<KGBEltList> d_cayley;
  std::vector<KGBEltPairList> d_inverseCayley;
  DescentList d_descent;
  std::vector<size_t> d_length;
  weyl::WeylEltList d_tau;
  gradings::StatusList d_status;

  bitset::BitSet<NumStates> d_state;

  bruhat::BruhatOrder* d_bruhat;

  const weyl::WeylGroup* d_weylGroup;

 public:

// constructors and destructors
  explicit KGB(size_t);

  explicit KGB(realredgp::RealReductiveGroup&);

  virtual ~KGB();

// copy, assignment and swap
  void swap(KGB&);

// accessors
  const bruhat::BruhatOrder& bruhatOrder() const {
    return *d_bruhat;
  }

  KGBElt cayley(size_t s, KGBElt x) const {
    return d_cayley[s][x];
  }

  bool compare(KGBElt, KGBElt) const;

  KGBElt cross(size_t s, KGBElt x) const {
    return d_cross[s][x];
  }

  const Descent& descent(KGBElt x) const {
    return d_descent[x];
  }

  KGBEltPair inverseCayley(size_t s, KGBElt x) const {
    return d_inverseCayley[s][x];
  }

  bool isAscent(size_t s, KGBElt x) const;

  bool isDescent(size_t s, KGBElt x) const;

  size_t length(KGBElt x) const {
    return d_length[x];
  }

  size_t rank() const {
    return d_rank;
  }

  size_t size() const {
    return d_length.size();
  }

  const gradings::Status& status(KGBElt x) const {
    return d_status[x];
  }

  gradings::Status::Value status(size_t s, KGBElt x) const {
    return d_status[x][s];
  }

  const weyl::WeylElt& tau(KGBElt) const;

  KGBEltPair tauPacket(const weyl::WeylElt&) const;

  const weyl::WeylGroup& weylGroup() const {
    return *d_weylGroup;
  }

  size_t weylLength(KGBElt) const;

// manipulators
  void fillBruhat(); // might throw

};

}

}

#endif
