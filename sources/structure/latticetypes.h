/*
  This is latticetypes.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef LATTICETYPES_H  /* guard against multiple inclusions */
#define LATTICETYPES_H

#include "latticetypes_fwd.h"

// these includes are not required by this file but we put them here
// so that other files don't have to know about matrix and bitvector
// directly.

#include "bitvector.h"
#include "matrix.h"

/******** function declarations **********************************************/

namespace atlas {

namespace latticetypes {

  LatticeElt& operator+= (LatticeElt&, const LatticeElt&);
  LatticeElt& operator-= (LatticeElt&, const LatticeElt&);
  LatticeElt& operator*= (LatticeElt&, const LatticeCoeff&);
  LatticeElt& operator/= (LatticeElt&, const LatticeCoeff&);
  LatticeElt& operator- (LatticeElt&);

  bool isZero(const LatticeElt&);

  LatticeCoeff scalarProduct(const LatticeElt&, const LatticeElt&);
  LatticeCoeff scalarProduct(const RatLatticeElt&, const LatticeElt&);

}

/******** type definitions ***************************************************/

namespace latticetypes {

class RatLatticeElt {

 private:

  LatticeElt d_num;
  LatticeCoeff d_denom;

 public:

// constructors and destructors
  RatLatticeElt() 
    {}

  RatLatticeElt(const LatticeElt& v, LatticeCoeff d)
    :d_num(v), d_denom(d) 
    {}

  RatLatticeElt(size_t n, LatticeCoeff d)
    :d_num(n), d_denom(d) 
    {}

  RatLatticeElt(const RatLatticeElt& v)
    :d_num(v.d_num), d_denom(v.d_denom) 
    {}

  ~RatLatticeElt() 
    {}

// accessors
  LatticeCoeff denominator() const {
    return d_denom;
  }

  const LatticeElt& numerator() const {
    return d_num;
  }

  size_t size() const {
    return d_num.size();
  }

//manipulators
  LatticeCoeff& denominator() {
    return d_denom;
  }

  LatticeElt& numerator() {
    return d_num;
  }

};

}

}

#endif
