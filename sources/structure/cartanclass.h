/*
  This is cartanclass.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef CARTANCLASS_H  /* guard against multiple inclusions */
#define CARTANCLASS_H

#include "cartanclass_fwd.h"

#include "rootdata_fwd.h"
#include "tori_fwd.h"

#include "bitmap.h"
#include "gradings.h"
#include "latticetypes.h"
#include "partition.h"
#include "realform.h"
#include "setutils.h"
#include "size.h"
#include "subquotient.h"

namespace atlas {

/******** function declarations **********************************************/

namespace cartanclass {

  void compactTwoRho(latticetypes::Weight&, unsigned long, const Fiber&, 
		     const rootdata::RootDatum&);

  void restrictGrading(gradings::Grading&, const rootdata::RootSet&,
		       const rootdata::RootList&);

  void specialGrading(gradings::Grading&, const Fiber&, realform::RealForm, 
		      const rootdata::RootDatum&);

  void toMostSplit(rootdata::RootList&, const Fiber&, realform::RealForm,
		   const rootdata::RootDatum&);
}

/******** type definitions ***************************************************/

namespace cartanclass {

class Fiber {

 protected:

  tori::RealTorus* d_torus;
  rootdata::RootSet d_complex;
  rootdata::RootSet d_imaginary;
  rootdata::RootSet d_real;
  rootdata::RootList d_simpleImaginary;
  setutils::Permutation d_rootInvolution;
  latticetypes::ComponentSubquotient d_fiberGroup;
  latticetypes::ComponentSubquotient d_adjointFiberGroup;
  latticetypes::ComponentMap d_toAdjoint;
  latticetypes::ComponentSubspace d_gradingGroup;
  latticetypes::ComponentList d_mAlpha;
  latticetypes::ComponentList d_adjointMAlpha;
  gradings::Grading d_baseGrading;
  gradings::GradingList d_gradingShift;
  rootdata::RootSet d_baseNoncompact;
  rootdata::RootSetList d_noncompactShift;
  partition::Partition d_weakReal;
  partition::Partition d_realFormPartition;
  std::vector<partition::Partition> d_strongReal;
  std::vector<StrongRealFormRep> d_strongRealFormReps;

 public:

// constructors and destructors

  Fiber():d_torus(0)
    {}

  Fiber(const rootdata::RootDatum&, const latticetypes::LatticeMatrix&);

  virtual ~Fiber();

// copy and assignment

  Fiber(const Fiber&);

  Fiber& operator= (const Fiber&);

  void swap(Fiber&);

// accessors
  const latticetypes::ComponentSubquotient& adjointFiberGroup() const {
    return d_adjointFiberGroup;
  }

  size_t adjointFiberRank() const {
    return d_adjointFiberGroup.dimension();
  }

  unsigned long adjointFiberSize() const {
    return d_adjointFiberGroup.size();
  }

  void compactRootSet(rootdata::RootSet&, unsigned long) const;

  const rootdata::RootSet& complexRootSet() const {
    return d_complex;
  }

  const latticetypes::ComponentSubquotient& fiberGroup() const {
    return d_fiberGroup;
  }

  size_t fiberRank() const {
    return d_fiberGroup.dimension();
  }

  unsigned long fiberSize() const {
    return d_fiberGroup.size();
  }

  void grading(gradings::Grading&, unsigned long) const;

  unsigned long gradingRep(const gradings::Grading&) const;

  const rootdata::RootSet& imaginaryRootSet() const {
    return d_imaginary;
  }

  const latticetypes::LatticeMatrix& involution() const;

  void mAlpha(latticetypes::Component&, const rootdata::Root&) const;

  size_t minusRank() const;

  void noncompactRootSet(rootdata::RootSet&, unsigned long) const;

  size_t numRealForms() const {
    return d_weakReal.classCount();
  }

  size_t plusRank() const;

  const partition::Partition& realFormPartition() const {
    return d_realFormPartition;
  }

  const rootdata::RootSet& realRootSet() const {
    return d_real;
  }

  rootdata::RootNbr rootInvolution(rootdata::RootNbr j) const {
    return d_rootInvolution[j];
  }

  const rootdata::RootList& simpleImaginary() const {
    return d_simpleImaginary;
  }

  const partition::Partition& strongReal(size_t j) const {
    return d_strongReal[j];
  }

  const StrongRealFormRep& strongRepresentative(size_t rf) const {
    return d_strongRealFormReps[rf];
  }

  unsigned long toAdjoint(unsigned long) const;

  unsigned long toWeakReal(unsigned long, size_t) const;

  const tori::RealTorus& torus() const {
    return *d_torus;
  }

  const partition::Partition& weakReal() const {
    return d_weakReal;
  }

 };

class CartanClass {

 private:

  Fiber d_fiber;
  Fiber d_dualFiber;
  rootdata::RootList d_simpleComplex;
  size::Size d_orbitSize;

public:

// constructors and destructors

  CartanClass() 
    {}

  CartanClass(const rootdata::RootDatum&, const latticetypes::LatticeMatrix&);

  ~CartanClass() 
    {}

// copy and assignment: defaults are ok for copy and assignment

  void swap(CartanClass&);

// accessors
  const Fiber& dualFiber() const {
    return d_dualFiber;
  }

  const Fiber& fiber() const {
    return d_fiber;
  }

  const rootdata::RootSet& imaginaryRootSet() const {
    return d_fiber.imaginaryRootSet();
  }

  const latticetypes::LatticeMatrix& involution() const {
    return d_fiber.involution();
  }

  bool isMostSplit(unsigned long) const;

  unsigned long numDualRealForms() const {
    return d_dualFiber.numRealForms();
  }

  unsigned long numRealForms() const {
    return d_fiber.numRealForms();
  }

  unsigned long numRealFormClasses() const {
    return d_fiber.realFormPartition().classCount();
  }

   unsigned long orbitSize() const {
    return d_orbitSize.toUlong();
  } 

  const rootdata::RootSet& realRootSet() const {
    return d_fiber.realRootSet();
  }

  rootdata::RootNbr rootInvolution(rootdata::RootNbr j) const {
    return d_fiber.rootInvolution(j);
  }

  const rootdata::RootList& simpleImaginary() const {
    return d_fiber.simpleImaginary();
  }

  const rootdata::RootList& simpleComplex() const {
    return d_simpleComplex;
  }

  const rootdata::RootList& simpleReal() const {
    return d_dualFiber.simpleImaginary();
  }

  const partition::Partition& strongReal(size_t j) const {
    return d_fiber.strongReal(j);
  }

  unsigned long toAdjoint(unsigned long x) const {
    return d_fiber.toAdjoint(x);
  }

  unsigned long toWeakReal(unsigned long c, size_t j) const {
    return d_fiber.toWeakReal(c,j);
  }

  const partition::Partition& weakReal() const {
    return d_fiber.weakReal();
  }

};

}

}

#endif
