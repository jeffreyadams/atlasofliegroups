/*
  This is realweyl.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef REALWEYL_H  /* guard against multiple inclusions */
#define REALWEYL_H


#include "../Atlas.h"

#include "bitvector.h"	// containment
#include "lietype.h"	// containment
#include "weyl.h"	// containment of (vector of) |WeylElt|

/******** function declarations **********************************************/

namespace atlas {

namespace realweyl {

  void blockStabilizerSize(size::Size&, const RealWeyl&);
  void dualRealWeylSize(size::Size&, const RealWeyl&);
  void realWeylSize(size::Size&, const RealWeyl&);

}

/******** type definitions ***************************************************/

namespace realweyl {

class RealWeyl {

  const WeylGroup* d_group;
  RootNbrList d_imaginaryCompact;
  RootNbrList d_imaginaryOrth;
  SmallBitVectorList d_imaginaryR;
  RootNbrList d_imaginary;
  RootNbrList d_realCompact;
  SmallBitVectorList d_realR;
  RootNbrList d_real;
  RootNbrList d_complex;
  LieType d_complexType;
  LieType d_imaginaryType;
  LieType d_imaginaryCompactType;
  RootNbrList d_realOrth;
  LieType d_realType;
  LieType d_realCompactType;

 public:

// constructors and destructors
  RealWeyl():d_group(0) {}

  RealWeyl(const CartanClass&,
	   cartanclass::AdjointFiberElt, cartanclass::AdjointFiberElt,
	   const RootDatum&, const WeylGroup&);

  ~RealWeyl() {}

// accessors
  const RootNbrList& complex() const {
    return d_complex;
  }

  RootNbr complex(size_t j) const {
    return d_complex[j];
  }

  const LieType& complexType() const {
    return d_complexType;
  }

  const RootNbrList& imaginary() const {
    return d_imaginary;
  }

  RootNbr imaginary(size_t j) const {
    return d_imaginary[j];
  }

  const RootNbrList& imaginaryCompact() const {
    return d_imaginaryCompact;
  }

  RootNbr imaginaryCompact(size_t j) const {
    return d_imaginaryCompact[j];
  }

  const LieType& imaginaryCompactType() const {
    return d_imaginaryCompactType;
  }

  RootNbr imaginaryOrth(size_t j) const {
    return d_imaginaryOrth[j];
  }

  const SmallBitVectorList& imaginaryR() const {
    return d_imaginaryR;
  }

  const LieType& imaginaryType() const {
    return d_imaginaryType;
  }

  const SmallBitVector& imaginaryR(size_t j) const {
    return d_imaginaryR[j];
  }

  size_t numComplex() const {
    return d_complex.size();
  }

  size_t numImaginary() const {
    return d_imaginary.size();
  }

  size_t numImaginaryCompact() const {
    return d_imaginaryCompact.size();
  }

  size_t numImaginaryR() const {
    return d_imaginaryR.size();
  }

   size_t numReal() const {
    return d_real.size();
  }

  size_t numRealCompact() const {
    return d_realCompact.size();
  }

  size_t numRealR() const {
    return d_realR.size();
  }

  const RootNbrList& real() const {
    return d_real;
  }

  RootNbr real(size_t j) const {
    return d_real[j];
  }

  const RootNbrList& realCompact() const {
    return d_realCompact;
  }

  RootNbr realCompact(size_t j) const {
    return d_realCompact[j];
  }

  const LieType& realCompactType() const {
    return d_realCompactType;
  }

  RootNbr realOrth(size_t j) const {
    return d_realOrth[j];
  }

  const SmallBitVectorList& realR() const {
    return d_realR;
  }

  const SmallBitVector& realR(size_t j) const {
    return d_realR[j];
  }

  const LieType& realType() const {
    return d_realType;
  }

  const WeylGroup& weylGroup() const {
    return *d_group;
  }
};

class RealWeylGenerators {

 private:
  const WeylGroup* d_group;
  WeylEltList d_imaginaryCompact;
  WeylEltList d_imaginaryR;
  WeylEltList d_imaginary;
  WeylEltList d_realCompact;
  WeylEltList d_realR;
  WeylEltList d_real;
  WeylEltList d_complex;

 public:
// constructors and destructors
  RealWeylGenerators():d_group(0) {}

  RealWeylGenerators(const RealWeyl&, const CartanClass&,
		     const RootDatum&);

  ~RealWeylGenerators() {}

// accessors
  const WeylEltList& complex() const {
    return d_complex;
  }

  const WeylEltList& imaginary() const {
    return d_imaginary;
  }

  const WeylEltList& imaginaryCompact() const {
    return d_imaginaryCompact;
  }

  const WeylEltList& imaginaryR() const {
    return d_imaginaryR;
  }

  const WeylEltList& real() const {
    return d_real;
  }

  const WeylEltList& realCompact() const {
    return d_realCompact;
  }

  const WeylEltList& realR() const {
    return d_realR;
  }

  const WeylGroup& weylGroup() const {
    return *d_group;
  }
};

}

}

#endif
