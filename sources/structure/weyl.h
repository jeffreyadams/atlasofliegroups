/*
  This is weyl.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef WEYL_H  /* guard against multiple inclusions */
#define WEYL_H

#include "weyl_fwd.h"

#include "constants.h"
#include "error.h"
#include "latticetypes.h"
#include "size.h"
#include "tags.h"

/******** type declarations *************************************************/

namespace atlas {

namespace weyl {

  typedef Generator WeylInterface[constants::RANK_MAX];
  typedef std::pair<EltPiece,Generator> ShiftValue;

  class RowBase;
  typedef RowBase ShiftRow;
  typedef RowBase OutRow;

  class Transducer;

}

/******** constant declarations *********************************************/

namespace weyl {

  const unsigned char UndefValue = constants::ucharMax;
  const EltPiece UndefEltPiece = UndefValue;
  const Generator UndefGenerator = UndefValue;
  const unsigned long UndefOrder = 0; // is never a valid value

  const Twist& identityTwist();

}

/******** function declarations *********************************************/

namespace weyl {

  void copy(Twist&, const Twist&);

  void makeReflections(WeylEltList&, const WeylGroup&);

}

/******** type definitions **************************************************/

namespace weyl {

class WeylElt {

 private:

  EltPiece d_data[constants::RANK_MAX];

 public:

// constructors and destructors
  WeylElt() {
    memset(d_data,0,constants::RANK_MAX);
  }

  explicit WeylElt(EltPiece x) {
    memset(d_data,x,constants::RANK_MAX);
  }

  ~WeylElt() 
    {}

// copy and assignment
  WeylElt(const WeylElt& w) {
    memcpy(d_data,w.d_data,constants::RANK_MAX);
  }

  WeylElt& operator=(const WeylElt& w) {
    memcpy(d_data,w.d_data,constants::RANK_MAX); 
    return *this;
  }

// accessors
  EltPiece operator[] (size_t j) const {
    return d_data[j];
  }

  bool operator< (const WeylElt& w) const {
    return memcmp(d_data,w.d_data,constants::RANK_MAX) < 0;
  }

  bool operator== (const WeylElt& w) const {
    return !memcmp(d_data,w.d_data,constants::RANK_MAX);
  }

  bool operator!= (const WeylElt& w) const {
    return memcmp(d_data,w.d_data,constants::RANK_MAX);
  }

// manipulators
  EltPiece& operator[] (size_t j) {
    return d_data[j];
  }

  void clear() {
    memset(d_data,0,constants::RANK_MAX);
  }
};

const WeylElt Identity;

class RowBase {

 private:

  unsigned char d_data[constants::RANK_MAX];

 public:

// constructors and destructors
  RowBase() {
    memset(d_data,UndefValue,constants::RANK_MAX);
  }

  ~RowBase() {}

// copy and assignment
  RowBase(const RowBase& r) {
    memcpy(d_data,r.d_data,constants::RANK_MAX);
  }

  RowBase& operator=(const RowBase& r) {
    memcpy(d_data,r.d_data,constants::RANK_MAX); return *this;
  }

// accessors
  Generator operator[] (size_t j) const {
    return d_data[j];
  }

// manipulators
  Generator& operator[] (size_t j) {
    return d_data[j];
  }
};

class Transducer {

 private:

  std::vector<ShiftRow> d_shift;
  std::vector<OutRow> d_out;
  std::vector<unsigned long> d_length;
  std::vector<WeylWord> d_piece;

 public:

// constructors and destructors
  Transducer() 
    {}

  Transducer(const latticetypes::LatticeMatrix&, size_t);

  ~Transducer() 
    {}

// accessors
  unsigned long length(EltPiece x) const {
    return d_length[x];
  }

  unsigned long maxlength() const {
    return d_length.back();
  }

  Generator out(EltPiece x, Generator s) const {
    return d_out[x][s];
  }

  EltPiece shift(EltPiece x, Generator s) const {
    return d_shift[x][s];
  }

  unsigned long size() const {
    return d_shift.size();
  }

  const WeylWord& wordPiece(EltPiece x) const {
    return d_piece[x];
  }

// this class should have no manipulators!
};

class WeylGroup {

 private:

  size_t d_rank;
  size::Size d_order;
  unsigned long d_maxlength;
  WeylElt d_longest;
  latticetypes::LatticeMatrix d_coxeterMatrix;
  std::vector<Transducer> d_transducer;
  Twist d_twist;
  WeylInterface d_in;
  WeylInterface d_out;

// private member functions
  void leftProdIn(WeylElt&, Generator) const;
  int prodIn(WeylElt&, Generator) const;
  void prodIn(WeylElt&, const WeylWord&) const;

  const WeylWord& wordPiece(const WeylElt& w, size_t j) const {
    const Transducer& tr = d_transducer[j];
    return tr.wordPiece(w[j]);
  }

 public:

// constructors and destructors
  WeylGroup()
    :d_rank(0UL),d_order(1UL),d_maxlength(0UL) {}

  WeylGroup(const latticetypes::LatticeMatrix&, const Twist* = 0);

  WeylGroup(const WeylGroup&, tags::DualTag);

  ~WeylGroup() {}

// accessors
  void conjugacyClass(WeylEltList&, const WeylElt&, bool twisted = true)
    const;

  void conjugate(WeylElt& w, Generator s) const {
    leftProd(w,s);
    prod(w,s);
  }

  bool hasDescent(Generator, const WeylElt&) const;

  bool hasTwistedCommutation(Generator, const WeylElt&) const;

  void invert(WeylElt&) const;

  unsigned long involutionLength(const weyl::WeylElt&) const;

  void involutionOut(WeylWord&, const WeylElt&) const;

  void invProd(WeylElt&, const WeylWord&) const;

  Generator leftDescent(const WeylElt&) const;

  void leftProd(WeylElt& w, Generator s) const {
    leftProdIn(w,d_in[s]);
  }

  unsigned long length(const WeylElt&) const;

  const WeylElt& longest() const {
    return d_longest;
  }

  unsigned long maxlength() const {
    return d_maxlength;
  }

  const size::Size& order() const {
    return d_order;
  }

  void out(WeylWord&, const WeylElt&) const;

  int prod(WeylElt& w, Generator s) const {
    return prodIn(w,d_in[s]);
  }

  void prod(WeylElt&, const WeylElt&) const;

  void prod(WeylElt&, const WeylWord&) const;

  size_t rank() const {
    return d_rank;
  }

  unsigned long toUlong(const WeylElt&) const;

  WeylElt toWeylElt(unsigned long) const;

  void translate(WeylElt&, const WeylInterface&) const;

  void twist(WeylElt&) const;

  Generator twistGenerator(Generator s) const {
    return d_twist[s];
  }

  void twistedConjugate(WeylElt& w, Generator s) const {
    leftProd(w,s);
    prodIn(w,d_twist[d_in[s]]);
  }

// manipulators
  void swap(WeylGroup&);
};

}

}

#endif
