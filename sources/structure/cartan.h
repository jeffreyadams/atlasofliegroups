/*!
\file
\brief Class definition and function declarations for CartanClasses.
*/

/*
  This is cartan.h
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

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

  /*!
   \brief Stores the set of stable conjugacy classes of Cartan subgroups of G.

Each stable conjugacy classes of Cartan subgroups corresponds to a W-conjugacy
class of involutions in the Gamma-enlarged Weyl group (W semidirect <Gamma>,
where <Gamma>=Z/2Z acts on W), contained in the complement of its subgroup W.
Since such involutions are of the form (w,Gamma), they can be represented by
their element w, which is called a twisted involution. The condition for being
a twisted involution $t$ is $t\Gamma(t)=e$ and "twisted conjugacy" of $t$ by
$w\in W$ is given by $w\cdot t=wt\Gamma(w^{-1})$. The stable conjugacy classes
of Cartan subgroups will each be represented by a canonical representative of
the corresponding twisted conjugacy class of twisted involutions.

In addition to describing the set of Cartan classes, this class provides
access (via the |d_cartan| array) to data for each individual one of them, and
(vie |d_ordering|) to the partial order relation between them. For the latter,
let |tau_i| be involutions acting on the complex torus |H| for various classes
of Cartan subgroups; (H,tau_1) is considered "more compact" than (H,tau_2)
if the identity component of the fixed point set H^tau_2 is W-conjugate to a
subtorus of H^tau_1.

The problem for the dual group of G is identical, the bijection taking
the negative transpose of a twisted involution.  This bijection
reverses the partial order on Cartans.  The class provides also access
to Cartans in the dual group.
  */

using weyl::TwistedInvolution;
using weyl::TwistedInvolutionList;

class CartanClasses {

  /* The following data members are accessible to our Helper class */
 protected:

  /*!
  \brief Fiber class for the fundamental Cartan subgroup.

  The involution is delta, which is stored here. It permutes the simple roots.
  */
  cartanclass::Fiber d_fundamental;

  /*!
  \brief Fiber class for the fundamental Cartan in the dual group.

  The fiber group here is the group of characters (i.e., the dual group)
  of the component group of the quasisplit Cartan.
  */
  cartanclass::Fiber d_dualFundamental;

  /*!
  \brief List of stable conjugacy classes of Cartan subgroups.

  The list includes only Cartans appearing in real forms considered so far.
  */
  std::vector<cartanclass::CartanClass*> d_cartan;

  /*!
  \brief (Representative) twisted involutions for each class of Cartan
  subgroup.
  */
  TwistedInvolutionList d_twistedInvolution;

  /*!
  \brief Partial order of Cartan subgroups.

  This is the ordering by containment of H^theta up to conjugacy: (H,theta_1)
  precedes (H,theta_2) if (H^theta_2)_0 is W-conjugate to a subtorus of
  H^theta_1. Numbering of elements is as in d_twistedInvolution
  */
  poset::Poset d_ordering;

  /*!
  \brief Entry \#n lists the real forms in which Cartan \#n is defined.

  More precisely, d_realFormLabels[n][i] is the inner number of the real form
  that corresponds to part i of the partition cartan(n).fiber().weakReal()

  */
  std::vector<realform::RealFormList> d_realFormLabels;

  /*!
  \brief Entry \#n lists the dual real forms in which dual Cartan \#n
    is defined.

  More precisely, d_dualRealFormLabels[n][i] is the inner number of the dual
  real form that corresponds to part i of the partition
  cartan(n).dualFiber().weakReal()
  */
  std::vector<realform::RealFormList> d_dualRealFormLabels;

  /*!
  \brief Entry \#rf flags the Cartans defined in real form \#rf.
  */
  std::vector<bitmap::BitMap> d_support;

  /*!
  \brief Entry \#rf flags the dual Cartans defined in dual real form \#rf.
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
  \brief Entry \#rf is the number of the most split Cartan for real form \#rf.
  */
  std::vector<size_t> d_mostSplit;

 public:

// constructors and destructors
  CartanClasses() {};

  CartanClasses(const rootdata::RootDatum&,
		const latticetypes::LatticeMatrix&,
		const weyl::WeylGroup&);

  ~CartanClasses();

// copy, assignment and swap
  CartanClasses(const CartanClasses&);

  CartanClasses& operator=(const CartanClasses&);

 private:
  // swap is needed to define Helper class, but should remain private
  void swap(CartanClasses&);

 public:
// accessors

  /*!
  \brief Returns data for stable conjugacy class \#cn of Cartan subgroups.
  */
  const cartanclass::CartanClass& cartan(size_t cn) const {
    return *d_cartan[cn];
  }

  /*!
  \brief Recover the matrix of the involution for the fundamental Cartan.

  This is the one permuting the simple roots, the distinguished one among the
  involutions in this inner class of G.
  */
  const latticetypes::LatticeMatrix& distinguished() const {
    return d_fundamental.involution();
  }

  /*!
  \brief Matrix of involution for the fundamental Cartan in the dual
  group.

  This is -w_0 times the transpose of the fundamental involution.
  */
  const latticetypes::LatticeMatrix& dualDistinguished() const {
    return d_dualFundamental.involution();
  }

  unsigned long dualFiberSize(realform::RealForm, size_t) const;

  /*!
  \brief Fiber class for the fundamental Cartan in the dual group.

  The fiber group here is the group of characters of the component
  group of the quasisplit Cartan.
  */
  const cartanclass::Fiber& dualFundamental() const {
    return d_dualFundamental;
  }

  /*!
  \brief Entry \#cn lists the dual real forms in which dual Cartan
  \#cn is defined.
  */
  const realform::RealFormList& dualRealFormLabels(size_t cn) const {
    return d_dualRealFormLabels[cn];
  }

  unsigned long dualRepresentative(realform::RealForm, size_t) const;

  /*!
  \brief Entry \#rf flags the dual Cartans defined in dual real form \#rf.
  */
  const bitmap::BitMap& dualSupport(realform::RealForm rf) const {
    return d_dualSupport[rf];
  }

  unsigned long fiberSize(realform::RealForm, size_t) const;

  /*!
  \brief Fiber class for the fundamental Cartan subgroup.

  The involution is delta, which preserves the simple roots.
  */
  const cartanclass::Fiber& fundamental() const {
    return d_fundamental;
  }

  /*!
  \brief Tells whether Cartan \#cn is defined in real form \#rf.
  */
  bool isDefined(realform::RealForm rf, size_t cn) const {
    return d_support[rf].isMember(cn);
  }

  /*!
  \brief Entry \#rf is the number of the most split Cartan for real form \#rf.
  */
  size_t mostSplit(realform::RealForm rf) const {
    return d_mostSplit[rf];
  }

  /*!
  \brief Flags in rs the set of noncompact imaginary roots for the
  fundamental Cartan in real form \#rf.
  */
  void noncompactRootSet(rootdata::RootSet& rs, realform::RealForm rf) const {
    d_fundamental.noncompactRootSet(rs,d_fundamental.weakReal().classRep(rf));
  }

  /*!
  \brief Returns the number of stable conjugacy classes of Cartans for G.
  */
  size_t numCartan() const {
    return d_cartan.size();
  }

  /*!
  \brief Returns the number of weak real forms of the dual group of G.
  */
  size_t numDualRealForms() const {
    return d_dualFundamental.numRealForms();
  }

  /*!
  \brief Returns the number of weak real forms of the dual group for
  which the dual of Cartan \#cn is defined.
  */
  size_t numDualRealForms(size_t cn) const {
    return d_cartan[cn]->numDualRealForms();
  }

  size_t numInvolutions() const;

  /*!
  \brief Returns the number of weak real forms of G.
  */
  size_t numRealForms() const {
    return d_fundamental.numRealForms();
  }

  /*!
  \brief Returns the number of weak real forms of G for which Cartan
  \#cn is defined.
  */
  size_t numRealForms(size_t cn) const {
    return d_cartan[cn]->numRealForms();
  }

  /*!
  \brief Returns the partial ordering of the set of Cartans.
  */
  const poset::Poset& ordering() const {
    return d_ordering;
  }

  /*!
  \brief Retruns the (inner) number of the quasisplit real form.
  */
  realform::RealForm quasisplit() const {
    return realform::RealForm(0);
  }

  /*!
  \brief Lists the real forms for which Cartan \#cn is defined.
  */
  const realform::RealFormList& realFormLabels(size_t cn) const {
    return d_realFormLabels[cn];
  }

  unsigned long representative(realform::RealForm, size_t) const;

  /*!
  \brief Entry \#rf flags the Cartans defined in real form \#rf.
  */
  const bitmap::BitMap& support(realform::RealForm rf) const {
    return d_support[rf];
  }

  /*!
  \brief (Representative) twisted involutions for each class of Cartan
  subgroup.
  */
  const TwistedInvolution& twistedInvolution(size_t cn) const {
    return d_twistedInvolution[cn];
  }

// manipulators
  void extend(const weyl::WeylGroup&, const rootdata::RootDatum&,
	      realform::RealForm);

};

}

}

#endif
