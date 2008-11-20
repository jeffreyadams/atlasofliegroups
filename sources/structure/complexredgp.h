/*!
\file
\brief Class definitions and function declarations for the class
ComplexReductiveGroup.
*/
/*
  This is complexredgp.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef COMPLEXREDGP_H  /* guard against multiple inclusions */
#define COMPLEXREDGP_H

#include "complexredgp_fwd.h"
#include "cartanset_fwd.h" // declares |CartanClassSet| sans certain inlines

#include "cartanclass.h"
#include "latticetypes.h"
#include "poset.h"
#include "rootdata.h"
#include "tits.h"
#include "weyl.h"

#include "lietype.h"
#include "realform.h"
#include "tags.h"
#include "setutils.h"

/******** function declarations **********************************************/

namespace atlas {

namespace complexredgp {

  void lieType(lietype::LieType&, const ComplexReductiveGroup&);

}

/******** type definitions ***************************************************/

namespace complexredgp {
  /*!
  \brief Complex reductive group endowed with an inner class of real
  forms.

  This class computes those aspects of the structure theory of (an
  inner class of) real reductive groups G(R) that will be needed to
  describe the Langlands classification of irreducible representations
  of G(R).  Since we look at an inner class of real forms, the first
  problem is to enumerate the different real forms constituting this
  inner class.

  We list in d_cartanSet the conjugacy classes of real Cartan subgroups up to
  stable conjugacy; this classification does not refer to a particular real
  form. However, the enumeration of the real forms actually takes place during
  the construction of the first (fundamental) Cartan subgroup. Each stable
  class corresponds to at most one conjugacy class of real Cartan subgroups in
  each real form; so for each stable class of Cartan subgroups, we enumerate
  the real forms over which it is defined (the fundamental Cartan subgroup is
  defined for all real forms).

  We compute the structure of the real Cartan subgroups (notably the
  groups of connected components); this depends only on the stable
  conjugacy class.  We determine the real Weyl groups of Cartan
  subgroups (which are _almost_ constant across the stable class, but
  not quite).

  Everything is determined by (and computed from) two things: the based root
  datum recorded in the RootDatum class d_rootDatum, and its involutive
  automorphism. Many computations take place inside the Tits group, which is
  an extension of the (complex) Weyl group by the elements of order 2 in the
  torus. (In fact the structure we store in |d_titsGroup| allows computing in
  an even larger group, the semidirect product of the Tits group just
  described by a factor Z/2Z whose action on the other factor is determined by
  the given automorphism of the based root datum, the "twist".)

  The actual structure of this class is subdivided into three parts,
  implemented by three other classes. The field |d_rootDatum| stores the root
  datum, which must have been consructed before. The field |d_titsGroup| holds
  the mentioned (enlarged) Tits group, which is constructed upon entry from
  the root datum and the involution; it also gives access to just the
  (complex) Weyl group when that is necessary. Finally |d_cartanSet| stores all
  the information relative to (stable conjugacy classes of) Cartan subgroupes
  and real forms.

  */
class ComplexReductiveGroup
{
  /*!
  \brief The based root datum.
  */
  const rootdata::RootDatum& d_rootDatum;

  /*!
  \brief The Tits group of the based root datum, extended by an
  involutive automorphism.
  */
  const tits::TitsGroup d_titsGroup;

  //!\brief the permutation of the roots given by the based automorphism
  const setutils::Permutation root_twist;

  /*!
  \brief Storage of data for each stable conjugacy class of Cartan subgroups
  of the inner class of real forms determined by the based root datum with
  involution.
  */
  cartanset::CartanClassSet d_cartanSet;

// copy, assignement and swap are forbidden, and should not be implemented
  ComplexReductiveGroup(const ComplexReductiveGroup&);
  ComplexReductiveGroup& operator= (const ComplexReductiveGroup&);
  void swap(ComplexReductiveGroup& G);

 public:
// constructors and destructors
  ComplexReductiveGroup(const rootdata::RootDatum*,
			const latticetypes::LatticeMatrix&);

  ComplexReductiveGroup(const ComplexReductiveGroup&, tags::DualTag);

  ~ComplexReductiveGroup();

// accessors

/* Most accessors are just forwarding functions to one of the components.

   N.B. Given the immense number of methods that are simply forwarded to
   |d_cartanSet|, one may wonder it it would no have been a better idea to
   derive this class from |CartanClassSet|. To that, one can oppose on one
   hand that an inner class "has" rather than "is" a set of Cartan classes (if
   it is anything, it is a class of real forms rather than Cartan subgroups,
   but even that is not the point of view taken in this software library), and
   on the other hand that this would force us to construct the equivalent of
   |d_cartanSet| _before_ the other data members, which would be problematic.
*/

/*!
  \brief the size of the block defined by the weak real form rf
  and the weak dual real form drf. Requires Cartan for |rf| to be generated
*/
  unsigned long blockSize(realform::RealForm rf, realform::RealForm drf) const;

/*!
  \brief returns Cartan subgroup number \#cn in the group.
*/
  const cartanclass::CartanClass& cartan(size_t cn) const
    { return d_cartanSet.cartan(cn); }

/*!
  \brief returns the ordering of the (currently generated) Cartan classes
*/
  const poset::Poset& cartanOrdering() const{ return d_cartanSet.ordering(); }

/*!
  \brief returns the set of Cartan classes for |rf|
  Requires those Cartan to be already generated
*/
  const bitmap::BitMap& cartanSet(realform::RealForm rf) const
    { return d_cartanSet.support(rf); }

/*!
  \brief returns the set of Cartan classes for the dual real form rf
  Only those Cartan classes already generated will show up
*/
  const bitmap::BitMap& dualCartanSet(realform::RealForm rf) const
    { return d_cartanSet.dualSupport(rf); }


 /*!
  \brief returns the matrix of the distinguished involution.
*/
 const latticetypes::LatticeMatrix& distinguished() const
   { return d_cartanSet.distinguished(); }


/*!
  \brief returns the matrix of the twisted involution |tw|.
*/
  latticetypes::LatticeMatrix involutionMatrix
    (const weyl::TwistedInvolution& tw) const
    { return d_cartanSet.involutionMatrix(tw); }


/*! \brief
  Returns the size of the fiber size corresponding to dual real form \#drf
  and cartan \#cn (which obviously must have been generated before)

  This is a technical function used for size computations.
*/
  unsigned long dualFiberSize(realform::RealForm drf, size_t cn) const
    { return d_cartanSet.dualFiberSize(drf,cn); }

/*!
  \brief returns the dual fundamental fiber.

  This is a technical data structure containing the data for the classification
  of weak and strong dual real forms of G.
*/
  const cartanclass::Fiber& dualFundamental() const
    { return d_cartanSet.dualFundamental(); }


/*!
  \brief returns the dual real form labels for cartan \#cn
*/
  const realform::RealFormList& dualRealFormLabels(size_t cn) const
    { return d_cartanSet.dualRealFormLabels(cn); }


/*!
  \brief returns an element of the orbit corresponding to drf in the
  classification of weak real forms for cartan \#cn.

  Precondition: cartan \#cn is defined for drf.
*/
  unsigned long dualRepresentative(realform::RealForm drf, size_t cn) const
    { return d_cartanSet.dualRepresentative(drf,cn); }

/*!
  \brief returns the size of the fiber size corresponding to real
  form \#rf and cartan \#cn.

  Explanation: this is the size of the orbits, for the shifted action of
  the imaginary Weyl group, that correspond to rf in the classification of
  strong real forms (they all have the same size.)

  This is a technical function used for size computations.
*/
  unsigned long fiberSize(realform::RealForm rf, size_t cn) const
    { return d_cartanSet.fiberSize(rf,cn); }


/*!
  \brief returns the fundamental fiber.

  This is a technical data structure containing the data for the classification
  of weak and strong real forms of G.
*/
  const cartanclass::Fiber& fundamental() const
    { return d_cartanSet.fundamental(); }

/*!
  \brief Returns the set of noncompact imaginary roots for the
  representative of rf
*/
  rootdata::RootSet noncompactRoots(realform::RealForm rf) const
    { return d_cartanSet.noncompactRoots(rf); }

/*!
  \brief returns the number of elements in K\\G/B for real form \#rf.

  Precondition: the Cartan classes for this real form have been generated

  Explanation: this is exactly the number of elements in the one-sided
  parameter set corresponding to any strong real form of G lying over rf.
*/
  unsigned long kgbSize(realform::RealForm rf) const;

/*!
  \brief returns the most split cartan subgroup for real form \#rf.

  Precondition: |fillCartan()| has been called, or at least |fillCartan(rf)|
*/
  size_t mostSplit(realform::RealForm rf) const
    { return d_cartanSet.mostSplit(rf); }

/*!
  \brief returns the number of conjugacy classes of Cartan subgroups
  currently constructed for G. Only after fillCartan() has been called is it
  ensured that this gives the total number of Cartan subgroups for this inner
  class.
*/
  size_t numCartanClasses() const { return d_cartanSet.numCartan(); }

/*!
  \brief returns the number of weak dual real forms for this inner class.
*/
  size_t numDualRealForms() const { return d_cartanSet.numDualRealForms(); }

/*!
  \brief returns the number of involutions for the currently defined Cartans.
*/
  size_t numInvolutions() const { return d_cartanSet.numInvolutions(); }

/*!
  \brief returns the number of involutions for the indicated Cartans.
*/
  size_t numInvolutions(const bitmap::BitMap& Cartan_classes) const
    { return d_cartanSet.numInvolutions(Cartan_classes); }

/*!
  \brief returns the number of weak real forms for this inner class.
*/
  size_t numRealForms() const{ return d_cartanSet.numRealForms(); }

/*!
  \brief returns the quasisplit real form.
*/
  realform::RealForm quasisplit() const { return d_cartanSet.quasisplit(); }

/*!
  \brief returns the rank of the group.
*/
  size_t rank() const { return d_rootDatum.rank(); }

/*!
  \brief returns the real form labels for cartan \#cn

  More precisely, realFormLabels(cn)[i] is the (inner) number of the real form
  that corresponds to part i of the partition cartan(n).fiber().weakReal()
*/
  const realform::RealFormList& realFormLabels(size_t cn) const
    { return d_cartanSet.realFormLabels(cn); }

/*!
  \brief returns an element of the orbit corresponding to rf in the
  classification of weak real forms for cartan \#cn.

  Precondition: cartan \#cn is defined for rf.
*/
  unsigned long representative(realform::RealForm rf, size_t cn) const
    { return d_cartanSet.representative(rf,cn); }

  const rootdata::RootDatum& rootDatum() const { return d_rootDatum; }

//!\brief returns the semisimple rank of the group.
  size_t semisimpleRank() const { return d_rootDatum.semisimpleRank(); }

  const weyl::WeylGroup& weylGroup() const { return d_titsGroup.weylGroup(); }

  //!\brief returns a reference to the Weyl group (owned by the Tits group).
  const tits::TitsGroup& titsGroup() const { return d_titsGroup; }

  const cartanset::CartanClassSet& cartanClasses() const
    { return d_cartanSet; }

  const setutils::Permutation& root_involution() const { return root_twist; }
  rootdata::RootNbr twisted_root(rootdata::RootNbr alpha) const
    { return root_twist[alpha]; }

/*!
  \brief the twisted involution representative for class \#cn of Cartans
*/
  const weyl::TwistedInvolution& twistedInvolution(size_t cn) const
    { return d_cartanSet.twistedInvolution(cn); }

// manipulators

/*!
  \brief fills in the Cartan classes that are defined for the real form |rf|.
 */
  void fillCartan(realform::RealForm rf) { d_cartanSet.extend(rf); }

  /* the following is not done via a default argument to the previous method
     since a default argument cannot refer to a class member (quasisplit) */
  void fillCartan() { fillCartan(quasisplit()); }

 private:
  setutils::Permutation make_root_twist() const;

}; // |class ComplexReductiveGroup|

} // |namespace complexredgp|

} // |namespace atlas|

#endif
