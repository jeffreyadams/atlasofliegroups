/*
  This is realweyl.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef REALWEYL_H  /* guard against multiple inclusions */
#define REALWEYL_H

#include "realweyl_fwd.h"

#include "cartanclass_fwd.h"
#include "rootdata_fwd.h"

#include "bitvector.h"
#include "lietype.h"
#include "setutils.h"
#include "size.h"
#include "weyl.h"

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

 private:

  const weyl::WeylGroup* d_group;
  rootdata::RootList d_imaginaryCompact;
  rootdata::RootList d_imaginaryOrth;
  latticetypes::ComponentList d_imaginaryR;
  rootdata::RootList d_imaginary;
  rootdata::RootList d_realCompact;
  latticetypes::ComponentList d_realR;
  rootdata::RootList d_real;
  rootdata::RootList d_complex;
  lietype::LieType d_complexType;
  lietype::LieType d_imaginaryType;
  lietype::LieType d_imaginaryCompactType;
  rootdata::RootList d_realOrth;
  lietype::LieType d_realType;
  lietype::LieType d_realCompactType;

 public:

// constructors and destructors
  RealWeyl():d_group(0) {}

  RealWeyl(const cartanclass::CartanClass&, unsigned long, unsigned long, 
	   const rootdata::RootDatum&, const weyl::WeylGroup&);

  ~RealWeyl() {}

// accessors
  const rootdata::RootList& complex() const {
    return d_complex;
  }

  rootdata::RootNbr complex(size_t j) const {
    return d_complex[j];
  }

  const lietype::LieType& complexType() const {
    return d_complexType;
  }

  const rootdata::RootList& imaginary() const {
    return d_imaginary;
  }

  rootdata::RootNbr imaginary(size_t j) const {
    return d_imaginary[j];
  }

  const rootdata::RootList& imaginaryCompact() const {
    return d_imaginaryCompact;
  }

  rootdata::RootNbr imaginaryCompact(size_t j) const {
    return d_imaginaryCompact[j];
  }

  const lietype::LieType& imaginaryCompactType() const {
    return d_imaginaryCompactType;
  }

  rootdata::RootNbr imaginaryOrth(size_t j) const {
    return d_imaginaryOrth[j];
  }

  const latticetypes::ComponentList& imaginaryR() const {
    return d_imaginaryR;
  }

  const lietype::LieType& imaginaryType() const {
    return d_imaginaryType;
  }

  const latticetypes::Component& imaginaryR(size_t j) const {
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

  const rootdata::RootList& real() const {
    return d_real;
  }

  rootdata::RootNbr real(size_t j) const {
    return d_real[j];
  }

  const rootdata::RootList& realCompact() const {
    return d_realCompact;
  }

  rootdata::RootNbr realCompact(size_t j) const {
    return d_realCompact[j];
  }

  const lietype::LieType& realCompactType() const {
    return d_realCompactType;
  }

  rootdata::RootNbr realOrth(size_t j) const {
    return d_realOrth[j];
  }

  const latticetypes::ComponentList& realR() const {
    return d_realR;
  }

  const latticetypes::Component& realR(size_t j) const {
    return d_realR[j];
  }

  const lietype::LieType& realType() const {
    return d_realType;
  }

  const weyl::WeylGroup& weylGroup() const {
    return *d_group;
  }
};

class RealWeylGenerators {

 private:
  const weyl::WeylGroup* d_group;
  weyl::WeylEltList d_imaginaryCompact;
  weyl::WeylEltList d_imaginaryR;
  weyl::WeylEltList d_imaginary;
  weyl::WeylEltList d_realCompact;
  weyl::WeylEltList d_realR;
  weyl::WeylEltList d_real;
  weyl::WeylEltList d_complex;

 public:
// constructors and destructors
  RealWeylGenerators():d_group(0) {}

  RealWeylGenerators(const RealWeyl&, const cartanclass::CartanClass&,
		     const rootdata::RootDatum&);

  ~RealWeylGenerators() {}

// accessors
  const weyl::WeylEltList& complex() const {
    return d_complex;
  }

  const weyl::WeylEltList& imaginary() const {
    return d_imaginary;
  }

  const weyl::WeylEltList& imaginaryCompact() const {
    return d_imaginaryCompact;
  }

  const weyl::WeylEltList& imaginaryR() const {
    return d_imaginaryR;
  }

  const weyl::WeylEltList& real() const {
    return d_real;
  }

  const weyl::WeylEltList& realCompact() const {
    return d_realCompact;
  }

  const weyl::WeylEltList& realR() const {
    return d_realR;
  }

  const weyl::WeylGroup& weylGroup() const {
    return *d_group;
  }
};

}

}

#endif
