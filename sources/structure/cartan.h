/*
  This is cartan.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef CARTAN_H  /* guard against multiple inclusions */
#define CARTAN_H

#include "cartan_fwd.h"
#include "rootdata_fwd.h"

#include "bitmap.h"
#include "cartanclass.h"
#include "poset.h"
#include "realform.h"
#include "weyl.h"

/******** function declarations **********************************************/

namespace atlas {

namespace cartan {

  unsigned long blockSize(realform::RealForm, realform::RealForm, 
			  const CartanClasses&);

  unsigned long kgbSize(realform::RealForm, const CartanClasses&);
}

/******** type definitions ***************************************************/

namespace cartan {

class CartanClasses {

 protected:

  cartanclass::Fiber d_fundamental;
  cartanclass::Fiber d_dualFundamental;
  std::vector<cartanclass::CartanClass*> d_cartan;
  weyl::WeylEltList d_twistedInvolution;
  poset::Poset d_ordering;
  std::vector<realform::RealFormList> d_realFormLabels;
  std::vector<realform::RealFormList> d_dualRealFormLabels;
  std::vector<bitmap::BitMap> d_support;
  std::vector<bitmap::BitMap> d_dualSupport;
  bitmap::BitMap d_status;
  std::vector<size_t> d_mostSplit;

 public:

// constructors and destructors
  CartanClasses() {};

  CartanClasses(const rootdata::RootDatum&, 
		const latticetypes::LatticeMatrix&,
		const weyl::WeylGroup&);

  virtual ~CartanClasses();

// copy, assignment and swap
  CartanClasses(const CartanClasses&);

  CartanClasses& operator=(const CartanClasses&);

  void swap(CartanClasses&);

// accessors
  const cartanclass::CartanClass& cartan(size_t j) const {
    return *d_cartan[j];
  }

  const latticetypes::LatticeMatrix& distinguished() const {
    return d_fundamental.involution();
  }

  const latticetypes::LatticeMatrix& dualDistinguished() const {
    return d_dualFundamental.involution();
  }

  unsigned long dualFiberSize(realform::RealForm, size_t) const;

  const cartanclass::Fiber& dualFundamental() const {
    return d_dualFundamental;
  }

  const realform::RealFormList& dualRealFormLabels(size_t j) const {
    return d_dualRealFormLabels[j];
  }

  unsigned long dualRepresentative(realform::RealForm, size_t) const;

  const bitmap::BitMap& dualSupport(realform::RealForm rf) const {
    return d_dualSupport[rf];
  }

  unsigned long fiberSize(realform::RealForm, size_t) const;

  const cartanclass::Fiber& fundamental() const {
    return d_fundamental;
  }

  bool isDefined(realform::RealForm rf, size_t j) const {
    return d_support[rf].isMember(j);
  }

  size_t mostSplit(realform::RealForm rf) const {
    return d_mostSplit[rf];
  }

  void noncompactRootSet(rootdata::RootSet& rs, realform::RealForm rf) const {
    d_fundamental.noncompactRootSet(rs,d_fundamental.weakReal().classRep(rf));
  }

  size_t numCartan() const {
    return d_cartan.size();
  }

  size_t numDualRealForms() const {
    return d_dualFundamental.numRealForms();
  }

  size_t numDualRealForms(size_t j) const {
    return d_cartan[j]->numDualRealForms();
  }

  size_t numRealForms() const {
    return d_fundamental.numRealForms();
  }

  size_t numRealForms(size_t j) const {
    return d_cartan[j]->numRealForms();
  }

  const poset::Poset& ordering() const {
    return d_ordering;
  }

  const realform::RealFormList& realFormLabels(size_t j) const {
    return d_realFormLabels[j];
  }

  unsigned long representative(realform::RealForm, size_t) const;

  const bitmap::BitMap& support(realform::RealForm rf) const {
    return d_support[rf];
  }

// manipulators
  void extend(const weyl::WeylGroup&, const rootdata::RootDatum&, 
	      realform::RealForm);

};

}

}

#endif
