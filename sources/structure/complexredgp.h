/*
  This is complexredgp.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef COMPLEXREDGP_H  /* guard against multiple inclusions */
#define COMPLEXREDGP_H

#include "complexredgp_fwd.h"

#include "cartan_fwd.h"
#include "cartanclass_fwd.h"
#include "latticetypes_fwd.h"
#include "poset_fwd.h"
#include "rootdata_fwd.h"
#include "tits_fwd.h"
#include "weyl_fwd.h"

#include "lietype.h"
#include "realform.h"
#include "tags.h"

/******** function declarations **********************************************/

namespace atlas {

namespace complexredgp {

  void lieType(lietype::LieType&, const ComplexReductiveGroup&);

}

/******** type definitions ***************************************************/

namespace complexredgp {

class ComplexReductiveGroup {

 private:

// this class is one of the outer interfaces for the structure library
// hence we use pointers for its data members so that forward declarations
// suffice
  const rootdata::RootDatum* d_rootDatum;
  tits::TitsGroup* d_titsGroup;
  cartan::CartanClasses* d_cartan;

// reserve and implement when necessary
  ComplexReductiveGroup(const ComplexReductiveGroup&);
  ComplexReductiveGroup& operator= (const ComplexReductiveGroup&);

 public:
// constructors and destructors
  ComplexReductiveGroup();

  ComplexReductiveGroup(const rootdata::RootDatum*, 
			const latticetypes::LatticeMatrix&);

  ComplexReductiveGroup(const ComplexReductiveGroup&, tags::DualTag);

  ~ComplexReductiveGroup();

// accessors
  unsigned long blockSize(realform::RealForm, realform::RealForm) const;

  const cartanclass::CartanClass& cartan(size_t) const;

  const poset::Poset& cartanOrdering() const;

  const bitmap::BitMap& cartanSet(realform::RealForm) const;

  const latticetypes::LatticeMatrix& distinguished() const;

  unsigned long dualFiberSize(realform::RealForm, size_t) const;

  const cartanclass::Fiber& dualFundamental() const;

  const realform::RealFormList& dualRealFormLabels(size_t) const;

  unsigned long dualRepresentative(realform::RealForm, size_t) const;

  unsigned long fiberSize(realform::RealForm, size_t) const;

  const cartanclass::Fiber& fundamental() const;

  void grading(rootdata::RootSet&, realform::RealForm) const;

  unsigned long kgbSize(realform::RealForm) const;

  size_t mostSplit(realform::RealForm) const;

  bool noInnerClass() const {
    return d_cartan == 0;
  }

  size_t numCartanClasses() const;

  size_t numDualRealForms() const;

  size_t numRealForms() const;

  size_t rank() const;

  const realform::RealFormList& realFormLabels(size_t) const;

  unsigned long representative(realform::RealForm, size_t) const;

  const rootdata::RootDatum& rootDatum() const {
    return *d_rootDatum;
  }

  size_t semisimpleRank() const;

  const weyl::WeylGroup& weylGroup() const;

  const tits::TitsGroup& titsGroup() const {
    return *d_titsGroup;
  }

// manipulators
  void fillCartan(realform::RealForm rf = 0);

  void swap(ComplexReductiveGroup& G);
};

}

}

#endif
