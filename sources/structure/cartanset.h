/*!
\file
\brief Class definition and function declarations for CartanClassSet.
*/

/*
  This is cartanset.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef CARTANSET_H  /* guard against multiple inclusions */
#define CARTANSET_H

#include "cartanset_fwd.h"
#include "rootdata_fwd.h"
#include "complexredgp.h"

#include "bitmap.h"
#include "cartanclass.h"
#include "poset.h"
#include "realform.h"
#include "weyl.h"

/******** function declarations **********************************************/

namespace atlas {

namespace cartanset {

  unsigned long blockSize(realform::RealForm, realform::RealForm,
		  	  const CartanClassSet&);

  unsigned long kgbSize(realform::RealForm, const CartanClassSet&);

  void cayley_and_cross_part(rootdata::RootList& cayley,
			     weyl::WeylWord& cross,
			     const weyl::TwistedInvolution& tw,
			     const rootdata::RootDatum& rd,
			     const weyl::WeylGroup& W);
}

/******** type definitions ***************************************************/

namespace cartanset {

using weyl::TwistedInvolution;
using weyl::TwistedInvolutionList;

  /*!
   \brief Stores the set of stable conjugacy classes of Cartan subgroups of G.

Each stable conjugacy classes of Cartan subgroups corresponds to a W-conjugacy
class of involutions in the Gamma-enlarged Weyl group (W semidirect <Gamma>,
where <Gamma>=Z/2Z acts on W), contained in the complement of its subgroup W.
Since such involutions are of the form (w,Gamma), they can be represented by
their element w, which is called a twisted involution. The condition for being
a twisted involution $t$ is $t\Gamma(t)=e$ and "twisted conjugacy" of $t$ by
\f$w\in W\f$ is given by \f$w\cdot t=wt\Gamma(w^{-1})\f$. The stable conjugacy classes
of Cartan subgroups will each be represented by a canonical representative of
the corresponding twisted conjugacy class of twisted involutions.

In addition to describing the set of Cartan classes, this class provides
access (via the |d_cartan| array) to data for each individual one of them, and
(via |d_ordering|) to the partial order relation between them. For the latter,
let |tau_i| be involutions acting on the complex torus |H| for various classes
of Cartan subgroups; (H,tau_1) is considered "more compact" than (H,tau_2)
if the identity component of the fixed point set H^tau_2 is W-conjugate to a
subtorus of H^tau_1.

The problem for the dual group of G is identical, the bijection taking
the negative transpose of a twisted involution.  This bijection
reverses the partial order on Cartans.  The class provides also access
to Cartans in the dual group.
  */

class CartanClassSet {

  /*!
  \brief The inner class to which we are associated (and accessed from)
  */
  const complexredgp::ComplexReductiveGroup& d_parent;

  /*!
  \brief List of stable conjugacy classes of Cartan subgroups.

  The list includes only Cartans appearing in real forms considered so far.
  */
  std::vector<cartanclass::CartanClass*> d_cartan;

  /*!
  \brief (Representative) twisted involutions for each class of Cartan
  subgroup. Choice at |j| must match corresponding |cartan(j).involution()|

  In this version the representative involutions are the canonical ones.
  They satisfy:

  (1) the sum of the positive real roots is dominant (call it $SR$)
  (2) in the subsystem of roots orthogonal to $SR$, which contains all the
      imaginary roots, the sum of the imaginary roots is dominant for the
      subsystem (call it $SI$)
  (3) in the subsystem of roots orthogonal to both $SR$ and $SI$, the
      involution corresponding to twisted involution fixes (globally) the
      dominant chamber of the subsystem (it permutes its simple roots).
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

  /*!
  \brief Entry \#j is the sum of the positive real roots for
  d_twistedInvolution[j].

  [Not yet implemented ]
  */
  // std::vector<latticetypes::LatticeElt> d_PosRealRootSum;

  /*!
  \brief Entry \#j is intended to flag the imaginary roots for
   d_twistedInvolution[j].

  [Not yet implemented.  DV 8/5/06.]
  */
  // rootdata::RootSetList d_imaginary;

  /*!
  \brief Entry \#j is intended to flag the simple imaginary roots for
  d_twistedInvolution[j].

  [Not yet implemented.  DV 8/5/06.]
  */
  // rootdata::RootSetList d_simple_imaginary;

 public:

// constructors and destructors

  // the main and only constructor
  CartanClassSet(const complexredgp::ComplexReductiveGroup& parent,
	         const latticetypes::LatticeMatrix& distinguished);

  ~CartanClassSet();

// copy, assignment and swap

/* the copy constructor and assignment are forbidden, since a copy would point
   to the parent without reciprocal relation. If a need to copy or assign
   |ComplexReductiveGroup| objects should arise, this should be implemented
   using the pseudo copy-constructor declared below (but as yet undefined).
   For similar reasons, a swap operation should not be defined.
 */

 private:
  CartanClassSet(const CartanClassSet&);
  CartanClassSet& operator= (const CartanClassSet&);

 public:
// pseudo copy-constructor that provides a new parent for the copy
  CartanClassSet(const complexredgp::ComplexReductiveGroup& parent,
		 const CartanClassSet&);

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

  /*!
  \brief The size of the fiber orbits corresponding to strong real forms lying
   over weak real form \#rf, in cartan \#cn (all orbits have the same size)

  Precondition: Real form \#rf is defined for cartan \#cn.
*/
unsigned long fiberSize(realform::RealForm rf, size_t cn) const;

  /*!
  \brief Fiber class for the fundamental Cartan subgroup.

  The involution is delta, which preserves the simple roots.
  */
  const cartanclass::Fiber& fundamental() const {
    return d_fundamental;
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
  \brief Returns the set of noncompact imaginary roots for (the representative
  in the adjoint fiber of) the real form \#rf.
  */
  rootdata::RootSet noncompactRoots(realform::RealForm rf) const
  {
    return
      d_fundamental.noncompactRoots(d_fundamental.weakReal().classRep(rf));
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
  size_t numInvolutions(const bitmap::BitMap& Cartan_classes) const;

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
  \brief Returns the (inner) number of the quasisplit real form.
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

 /*!
  \brief Entry \#cn lists the dual real forms in which dual Cartan
  \#cn is defined.
  */
  const realform::RealFormList& dualRealFormLabels(size_t cn) const {
    return d_dualRealFormLabels[cn];
  }

/*!
  \brief An element of the orbit in the adjoint fiber corresponding to |rf|
  in the classification of weak real forms for cartan |\#cn|.

  Precondition: cartan \#cn is defined for rf.
*/
  unsigned long representative(realform::RealForm rf, size_t cn) const;

/*!
  \brief An element of the orbit in the adjoint dual fiber corresponding to
  |drf| in the classification of dual weak real forms for cartan |\#cn|.

  Precondition: cartan \#cn is defined for rf.
*/
  unsigned long dualRepresentative(realform::RealForm, size_t) const;

  const rootdata::RootDatum& rootDatum() const {
    return d_parent.rootDatum();
  }

  /*!
  \brief Entry \#rf flags the Cartans defined in real form \#rf.
  */
  const bitmap::BitMap& support(realform::RealForm rf) const {
    return d_support[rf];
  }

  /*!
  \brief Entry \#rf flags the dual Cartans defined in dual real form \#rf.
  */
  const bitmap::BitMap& dualSupport(realform::RealForm rf) const {
    return d_dualSupport[rf];
  }

  const weyl::WeylGroup& weylGroup() const {
    return d_parent.weylGroup();
  }

  /*!
  \brief (Representative) twisted involutions for each class of Cartan
  subgroup.
  */
  const TwistedInvolution& twistedInvolution(size_t cn) const {
    return d_twistedInvolution[cn];
  }

/*! \brief Make |sigma| canonical and return Weyl group |w| element that
    twisted conjugates the canonical representative back to |sigma|
*/
  const weyl::WeylElt canonicalize(TwistedInvolution&) const;


/*!
\brief matrix giving involution action of |tw| on weight lattice
*/
  latticetypes::LatticeMatrix
    involutionMatrix(const TwistedInvolution& tw) const;


/*!
\brief apply involution action of |tw| on weight lattice to |v|
*/
  void twistedAct
    (const weyl::TwistedInvolution& tw,latticetypes::LatticeElt& v)
  const;

  unsigned long KGB_size(realform::RealForm rf,
			 const bitmap::BitMap& Cartan_classes) const;
  unsigned long block_size(realform::RealForm, realform::RealForm,
			   const bitmap::BitMap& Cartan_classes) const;

  latticetypes::LatticeElt posRealRootSum(const TwistedInvolution&) const;

  latticetypes::LatticeElt posImaginaryRootSum(const TwistedInvolution&) const;


  size_t cayley(size_t, size_t, weyl::WeylElt*) const;

  size_t classNumber(TwistedInvolution) const;

// manipulators
  void extend(realform::RealForm);


 private:
// auxiliary accessors
void leftReflect(weyl::TwistedInvolution&, rootdata::RootNbr) const;

rootdata::RootSet noncompactPosRootSet(realform::RealForm, size_t) const;

std::vector<weyl::WeylEltList> expand() const; // obsolete

// auxiliary manipulators

/*!
  \brief Adds a new cartan, with Cartan involution given by |tw|.
*/
void addCartan(TwistedInvolution tw)
{
  d_cartan.push_back(new cartanclass::CartanClass
		     (rootDatum(),involutionMatrix(tw)));
}

/*!
  \brief Adds a new cartan, obtained from cartan \#j by Cayley transform
  through root \#rn. This unused old version will be phased out.
*/
void addCartan(rootdata::RootNbr rn, size_t j);

size_t addCanonicalRep(TwistedInvolution tw,TwistedInvolutionList& dest);

void correlateForms();

void correlateDualForms(const rootdata::RootDatum& rd,
			const weyl::WeylGroup& W);

void updateStatus(size_t prev_Cartan_size);

void updateSupports(size_t last_Cartan_class_added);

void updateTwistedInvolutions
  (std::vector<weyl::WeylEltList>& known, const TwistedInvolution& tw);

};

}

}

#endif
