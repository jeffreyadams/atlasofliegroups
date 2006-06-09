/*!
\file
\brief Class definitions and function declarations for CartanClasses.
*/
/*
  This is cartan.h
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

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

  /*!
  \brief Fiber class for the fundamental Cartan subgroup.

  The involution is delta, which preserves the simple roots.
  */ 
  cartanclass::Fiber d_fundamental;

  /*!
  \brief Fiber class for the fundamental Cartan in the dual group.

  The fiber group here is the group of characters of the component
  group of the quasisplit Cartan.
  */ 
  cartanclass::Fiber d_dualFundamental;

  /*!
  \brief List of stable conjugacy classes of Cartan subgroups.

  The list includes only Cartans appearing in real forms considered so far. 
  */ 
  std::vector<cartanclass::CartanClass*> d_cartan;

  /*!
  \brief The twisted involutions for each class of Cartan subgroup.
  */ 
  weyl::WeylEltList d_twistedInvolution;

  /*!
  \brief Partial order of Cartan subgroups.
  */ 
  poset::Poset d_ordering;

  /*!
  \brief Entry \#cn lists the real forms in which Cartan \#cn is defined.
  */ 
  std::vector<realform::RealFormList> d_realFormLabels;

  /*!
  \brief Entry \#cn lists the dual real forms in which dual Cartan
  \#cn is defined.
  */ 
  std::vector<realform::RealFormList> d_dualRealFormLabels;

  /*!
  \brief Entry rf flags the Cartans defined in real form rf. 
  */ 
  std::vector<bitmap::BitMap> d_support;

  /*!
  \brief Entry rf flags the dual Cartans defined in dual real form rf.
  */ 
  std::vector<bitmap::BitMap> d_dualSupport;

  /*!
  \brief  Flags the set of real forms for which the full set of Cartan
  classes is constructed.

  Because of the way the construction proceeds, these are exactly the
  real forms for which the most split Cartan has been reached.
  */ 
  bitmap::BitMap d_status;

  /*!
  \brief Entry rf is the number of the most split Cartan for real form rf.
  */ 
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

  /*!
  \brief 
  */
  const cartanclass::CartanClass& cartan(size_t cn) const {
    return *d_cartan[cn];
  }

  /*!
  \brief 
  */
  const latticetypes::LatticeMatrix& distinguished() const {
    return d_fundamental.involution();
  }

  /*!
  \brief 
  */
  const latticetypes::LatticeMatrix& dualDistinguished() const {
    return d_dualFundamental.involution();
  }

  unsigned long dualFiberSize(realform::RealForm, size_t) const;

  /*!
  \brief 
  */
  const cartanclass::Fiber& dualFundamental() const {
    return d_dualFundamental;
  }

  /*!
  \brief 
  */
  const realform::RealFormList& dualRealFormLabels(size_t cn) const {
    return d_dualRealFormLabels[cn];
  }

  unsigned long dualRepresentative(realform::RealForm, size_t) const;

  /*!
  \brief 
  */
  const bitmap::BitMap& dualSupport(realform::RealForm rf) const {
    return d_dualSupport[rf];
  }

  unsigned long fiberSize(realform::RealForm, size_t) const;

  /*!
  \brief 
  */
  const cartanclass::Fiber& fundamental() const {
    return d_fundamental;
  }

  /*!
  \brief 
  */
  bool isDefined(realform::RealForm rf, size_t cn) const {
    return d_support[rf].isMember(cn);
  }

  /*!
  \brief 
  */
  size_t mostSplit(realform::RealForm rf) const {
    return d_mostSplit[rf];
  }

  /*!
  \brief 
  */
  void noncompactRootSet(rootdata::RootSet& rs, realform::RealForm rf) const {
    d_fundamental.noncompactRootSet(rs,d_fundamental.weakReal().classRep(rf));
  }

  /*!
  \brief 
  */
  size_t numCartan() const {
    return d_cartan.size();
  }

  /*!
  \brief 
  */
  size_t numDualRealForms() const {
    return d_dualFundamental.numRealForms();
  }

  /*!
  \brief 
  */
  size_t numDualRealForms(size_t cn) const {
    return d_cartan[cn]->numDualRealForms();
  }

  size_t numInvolutions() const;

  /*!
  \brief 
  */
  size_t numRealForms() const {
    return d_fundamental.numRealForms();
  }

  /*!
  \brief 
  */
  size_t numRealForms(size_t cn) const {
    return d_cartan[cn]->numRealForms();
  }

  /*!
  \brief 
  */
  const poset::Poset& ordering() const {
    return d_ordering;
  }

  /*!
  \brief 
  */
  realform::RealForm quasisplit() const {
    return 0ul;
  }

  /*!
  \brief 
  */
  const realform::RealFormList& realFormLabels(size_t cn) const {
    return d_realFormLabels[cn];
  }

  unsigned long representative(realform::RealForm, size_t) const;

  /*!
  \brief 
  */
  const bitmap::BitMap& support(realform::RealForm rf) const {
    return d_support[rf];
  }

  /*!
  \brief 
  */
  const weyl::WeylElt& twistedInvolution(size_t cn) const {
    return d_twistedInvolution[cn];
  }

// manipulators
  void extend(const weyl::WeylGroup&, const rootdata::RootDatum&, 
	      realform::RealForm);

};

}

}

#endif
