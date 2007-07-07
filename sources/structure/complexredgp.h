/*!
\file
\brief Class definitions and function declarations for the class
ComplexReductiveGroup.
*/
/*
  This is complexredgp.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef COMPLEXREDGP_H  /* guard against multiple inclusions */
#define COMPLEXREDGP_H

#include "complexredgp_fwd.h"

#include "cartanset_fwd.h"
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
  the given automorphism of the based root datum.)

  The actual structure of this class is subdivided into three parts,
  implemented by three other classes. The field |d_rootDatum| stores the root
  datum, which must have been consructed before. The field |d_titsGroup| holds
  the mentioned (enlarged) Tits group, which is constructed upon entry from
  the root datum and the involution; it also gives access to just the
  (complex) Weyl group when that is necessary. Finally |d_cartanSet| stores all
  the information relative to (stable conjugacy classes of) Cartan subgroupes
  and real forms.

  Because this class is one of the outer interfaces for the structure library,
  we prefer to use references for its data members, so that forward
  declarations of their classes suffice for users of this class. The choice to
  use references instead of pointers (as had been done intially) is a
  deliberate one to emphasise the rigid connection between the parts; there
  should be no operation that can dissociate the three parts from the
  top-level structure grouping them together (such as |swap| initially did).
  */
class ComplexReductiveGroup {

 private:

  /*!
  \brief The based root datum.
  */
  const rootdata::RootDatum& d_rootDatum;

  /*!
  \brief The Tits group of the based root datum, extended by an
  involutive automorphism.
  */
  tits::TitsGroup& d_titsGroup;

  /*!
  \brief Storage of data for each stable conjugacy class of Cartan subgroups
  of the inner class of real forms determined by the based root datum with
  involution.
  */
  cartanset::CartanClassSet& d_cartanSet;

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
  unsigned long blockSize(realform::RealForm, realform::RealForm) const;

  const cartanclass::CartanClass& cartan(size_t) const;

  const poset::Poset& cartanOrdering() const;

  const bitmap::BitMap& cartanSet(realform::RealForm) const;

  const bitmap::BitMap& dualCartanSet(realform::RealForm) const;

  const latticetypes::LatticeMatrix& distinguished() const;

  latticetypes::LatticeMatrix involutionMatrix(const weyl::TwistedInvolution&)
    const;

  unsigned long dualFiberSize(realform::RealForm, size_t) const;

  const cartanclass::Fiber& dualFundamental() const;

  const realform::RealFormList& dualRealFormLabels(size_t) const;

  unsigned long dualRepresentative(realform::RealForm, size_t) const;

  unsigned long fiberSize(realform::RealForm, size_t) const;

  const cartanclass::Fiber& fundamental() const;

  rootdata::RootSet noncompactRoots(realform::RealForm) const;

  unsigned long kgbSize(realform::RealForm) const;

  size_t mostSplit(realform::RealForm) const;

  size_t numCartanClasses() const;

  size_t numDualRealForms() const;

  size_t numInvolutions() const;

  size_t numRealForms() const;

  realform::RealForm quasisplit() const;

  size_t rank() const;

  const realform::RealFormList& realFormLabels(size_t) const;

  unsigned long representative(realform::RealForm, size_t) const;

  const rootdata::RootDatum& rootDatum() const {
    return d_rootDatum;
  }

  size_t semisimpleRank() const;

  const weyl::WeylGroup& weylGroup() const;

  const tits::TitsGroup& titsGroup() const {
    return d_titsGroup;
  }

  const cartanset::CartanClassSet& cartanClasses() const {
    return d_cartanSet;
  }

/*!
  \brief the twisted involution representative for class \#cn of Cartans
*/
  const weyl::TwistedInvolution& twistedInvolution(size_t) const;

// manipulators

/* actually |fillCartan| could be declared |const|, since it only modifies the
   value referred to by |d_cartanSet|, not any members of the class itself
   (technically, a const method may use a reference member in a non-const way).
   However, morally it is a manipulator method, not an accessor, whence there
   is no const. In the future a level of indirection could be removed, and in
   that case quilfying fillCartan as const would even become impossible */

  void fillCartan(realform::RealForm rf);

  /* the following is not done via a default argument to the previous method
     since a default argument cannot refer to a class member (quasisplit) */
  void fillCartan() { fillCartan(quasisplit()); }
};

}

}

#endif
