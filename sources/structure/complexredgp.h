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
  /*!
  \brief Complex reductive group endowed with an inner class of real
  forms.

  This class computes those aspects of the structure theory of (an
  inner class of) real reductive groups G(R) that will be needed to
  describe the Langlands classification of irreducible representations
  of G(R).  Since we look at an inner class of real forms, the first
  problem is to enumerate the different real forms constituting this
  inner class.  

  Next we list in d_cartan the conjugacy classes of real Cartan
  subgroups up to stable conjugacy; this classification does not refer
  to a particular real form.  Each stable class corresponds to at most
  one conjugacy class of real Cartan subgroups in each real form; so
  for each stable class of Cartan subgroups, we enumerate the real
  forms over which it is defined.

  We compute the structure of the real Cartan subgroups (notably the
  groups of connected components); this depends only on the stable
  conjugacy class.  We determine the real Weyl groups of Cartan
  subgroups (which are _almost_ constant across the stable class, but
  not quite).

  Everything is determined by (and computed from) two things: the
  based root datum recorded in the RootDatum class d_rootDatum, and
  its involutive automorphism (which is stored inside d_cartan).
  involutive automorphism.  The computations take place mostly inside
  the Tits group d_titsGroup, which is an extension of the (complex)
  Weyl group by the elements of order 2 in the torus.  (More
  precisely, this Tits group is extended by a Z/2Z corresponding to
  the automorphism of the based root datum.)
  
  Because this class is one of the outer interfaces for the structure
  library, we use pointers for its data members, so that forward
  declarations suffice.
  */
class ComplexReductiveGroup {

 private:

// this class is one of the outer interfaces for the structure library
// hence we use pointers for its data members so that forward declarations
// suffice

  /*!
  \brief The based root datum.
  */
  const rootdata::RootDatum* d_rootDatum;
  
  /*!
  \brief The Tits group of the based root datum, extended by an
  involutive automorphism.
  */  
  tits::TitsGroup* d_titsGroup;
  
  /*!
  \brief  List of stable conjugacy classes of Cartan subgroups of the inner
  class of real forms determined by the based root datum with
  involution.  

  (In fact the present constructors provide only those classes defined
  over the real forms that have already been considered by the
  software: if the software has not yet been asked to look at the
  quasisplit real form, then this list will be incomplete.  This
  should probably be regarded as a defect in the software, although it
  causes no mathematical problem.  A related difficulty is that the
  ordering of the list can depend on the order in which real forms
  have been considered.)
  */
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

  size_t numInvolutions() const;

  size_t numRealForms() const;

  realform::RealForm quasisplit() const;

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

  const weyl::WeylElt& twistedInvolution(size_t) const;

// manipulators

  void fillCartan(realform::RealForm rf = 0);
  // default value 0 is d_cartan->quasisplit() but we avoid that reference
  // since it would require us to include cartan.h here
  // the default value adds Cartan subgroups for all possible real forms

  void swap(ComplexReductiveGroup& G);
};

}

}

#endif
