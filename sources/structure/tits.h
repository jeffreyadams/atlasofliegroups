/*!
\file
  This is tits.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef TITS_H  /* guard against multiple inclusions */
#define TITS_H

#include "tits_fwd.h"

#include "rootdata_fwd.h"

#include "bitvector.h"
#include "constants.h"
#include "latticetypes.h"
#include "weyl.h"

namespace atlas {

/******** type declarations *************************************************/

namespace tits {

  typedef bitvector::BitVector<constants::RANK_MAX> TorusPart;
  typedef bitvector::BitMatrix<constants::RANK_MAX> TorusMatrix;

}

/******** function declarations *********************************************/

namespace tits {

  void makeTwist(weyl::Twist&, const latticetypes::LatticeMatrix&,
		 const rootdata::RootDatum&);

}

/******** type definitions **************************************************/

namespace tits {

class TitsElt {

 private:
  
  weyl::WeylElt d_w;
  TorusPart d_t;
  
 public:

// constructors and destructors

/*! 
  constructs the identity element in the group
*/
  explicit TitsElt(size_t n):d_t(n) {} 

  TitsElt(const weyl::WeylElt& w, size_t n):d_w(w),d_t(n) {}

// copy and assignment can be left to their defaults

// accessors
  bool operator< (const TitsElt& a) const {
    if (d_w == a.d_w)
      return d_t < a.d_t;
    else
      return d_w < a.d_w;
  }

  const TorusPart& t() const {
    return d_t;
  }

  const weyl::WeylElt& w() const {
    return d_w;
  }

// manipulators

  TitsElt& operator+= (const TorusPart& x) {
    d_t += x;
    return *this;
  }

  TorusPart& t() {
    return d_t;
  }

  weyl::WeylElt& w() {
    return d_w;
  }

};

class TitsGroup {

 private:

  weyl::WeylGroup* d_weyl; // the underlying Weyl group
  
  size_t d_rank;
  std::vector<TorusPart> d_simpleRoot;
  std::vector<TorusPart> d_simpleCoroot;

  TorusMatrix d_involution;
  weyl::Twist d_twist;

// copy and assignment
// reserve and implement when necessary
  TitsGroup(const TitsGroup&);
  TitsGroup& operator= (const TitsGroup&);

 public:

// constructors and destructors
  TitsGroup() {}

  TitsGroup(const rootdata::RootDatum&, const latticetypes::LatticeMatrix&);

  ~TitsGroup();

// accessors
  void conjugate(TitsElt& a, weyl::Generator s) const {
/*
  note: in the Tits group s^{-1} is s.m_s!
*/
    leftProd(a,s);
    prod(a,s);
    a += d_simpleCoroot[s];
  }

  void invConjugate(TorusPart&, const weyl::WeylElt&) const;

  void leftProd(TitsElt&, weyl::Generator) const;

  unsigned long length(const TitsElt& a) const {
    return d_weyl->length(a.w());
  }

  void prod(TitsElt&, weyl::Generator) const;

  void prod(TitsElt&, const TitsElt&) const;

  const size_t rank() const {
    return d_rank;
  }

  void reflection(TorusPart& x, weyl::Generator s) const {
/*
  note: s must be an _outer_ generator
*/
    if (bitvector::scalarProduct(x,d_simpleRoot[s]))
      x += d_simpleCoroot[s];
  }

  const TorusPart& simpleCoroot(size_t j) const {
    return d_simpleCoroot[j];
  }

  const TorusPart& simpleRoot(size_t j) const {
    return d_simpleRoot[j];
  }

  size_t twist(size_t j) const {
    return d_twist[j];
  }

  void twistedConjugate(TitsElt& a, weyl::Generator s) const {
/*
  note: in the Tits group s^{-1} is s.m_s!
*/
    leftProd(a,s);
    prod(a,d_twist[s]);
    a += d_simpleCoroot[d_twist[s]];
  }

  const weyl::WeylGroup& weylGroup() const {
    return *d_weyl;
  }

// manipulators
  void swap(TitsGroup&);
};

}

}

#endif
