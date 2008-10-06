/*!
\file
  \brief Definitions and declarations for the CartanClass and Fiber classes.
*/

/*
  This is cartanclass.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef CARTANCLASS_H  /* guard against multiple inclusions */
#define CARTANCLASS_H

#include "cartanclass_fwd.h"

#include "tori_fwd.h"

#include "tags.h"
#include "bitmap.h"
#include "gradings.h"
#include "latticetypes.h"
#include "partition.h"
#include "realform.h"
#include "rootdata.h"
#include "setutils.h"
#include "size.h"
#include "subquotient.h"

namespace atlas {

/******** function declarations **********************************************/

namespace cartanclass {

  // quasi-constructor function for |Fiber| class
  Fiber dualFiber
    (const rootdata::RootDatum&, const latticetypes::LatticeMatrix&);

  latticetypes::Weight
    compactTwoRho(AdjointFiberElt, const Fiber&, const rootdata::RootDatum&);

  gradings::Grading
    restrictGrading(const rootdata::RootSet&, const rootdata::RootList&);

  gradings::Grading
    specialGrading(const Fiber&,realform::RealForm,const rootdata::RootDatum&);

  rootdata::RootList
    toMostSplit(const Fiber&,realform::RealForm,const rootdata::RootDatum&);
}

/******** type definitions ***************************************************/

namespace cartanclass {

// this class gathers information associated to a root datum involution
class InvolutionData
{
  setutils::Permutation d_rootInvolution; // permutation of all roots
  rootdata::RootSet d_imaginary, d_real, d_complex;
  rootdata::RootList d_simpleImaginary; // imaginary roots simple wrt subsystem
  rootdata::RootList d_simpleReal; // real roots simple wrt subsystem
 public:
  InvolutionData(const rootdata::RootDatum&,
		 const latticetypes::LatticeMatrix&);
  void swap(InvolutionData&);
  //accessors
  const setutils::Permutation& root_involution() const
    { return d_rootInvolution; }
  rootdata::RootNbr root_involution(rootdata::RootNbr alpha) const
    { return d_rootInvolution[alpha]; }
  const rootdata::RootSet& imaginary_roots() const  { return d_imaginary; }
  const rootdata::RootSet& real_roots() const       { return d_real; }
  const rootdata::RootSet& complex_roots() const    { return d_complex; }
  const rootdata::RootList& imaginary_basis() const
    { return d_simpleImaginary; }
  const rootdata::RootList& real_basis() const
    { return d_simpleReal; }
};

/*!
\brief Describes "the fiber" (over a fixed involution of a torus) in
twisted Tits group.

This is an important but somewhat subtle class. In fact none of the data
stored directly describes the mentioned fiber, nor is it relevant here where
the "elements of the fiber" really live. There is a "fiber group" (which is
a vector space over $Z/2Z$) that acts simply transitively on the fiber, in
other words, the fiber is an affine space over $Z/2Z$ with the fiber group
as space of translations. After a choice of a base point, the fiber itself
can be identified which this fiber group. The way a fiber element manifests
itself is via a grading it defines on the imaginary roots; this grading is
determined for the entire fiber by giving the "base grading" for the base
point, and the "grading shifts" for translations by each of a set of
generators of the grading group.

We fix an involutive automorphism \f$\tau\f$ of the complex torus $H$, and
consider the collection of all strong real forms $x$ of $G$ (which are
elements of "G semidirect delta") that induce tau on H, modulo the conjugation
action of H on the collection. Such a strong real form has square equal to
some element z in Z(G)^delta. Changing the element z by (1+delta)Z is
innocuous, and there are finitely many possibilities for Z^delta/[(1+delta)Z].
For each fixed value of z, this collection of strong real form representatives
in \f$G.\delta\f$ (modulo H conjugation) has a simply transitive action (left
multiplication) of the fiber group H^{-tau}/(1-tau)H, the component group of
the complex group H^{-tau} of fixed points of -tau. The fiber group is a
$Z/2Z$ vector space that can be realized as a subquotient of the group H(2) of
elements of order 2 (or 1) in H.

Equivalence classes of strong real forms (still with a fixed value
of z) give a partition of the fiber group (depending on z).  These
partitions are stored in d_strongReal, and accessed by the function
strongReal.

An equivalence class of strong real forms may be labelled by a pair
of integers: the second labelling the element z (or rather its coset
modulo [(1+delta)Z]), and the first the equivalence class in the
fiber group.
*/
class Fiber {

 private:

  /*!
  \brief Pointer to a torus defined over R.

  Represented as the lattice Z^n endowed with an involutive
  automorphism (represented by its n x n integer matrix).
  */
  tori::RealTorus* d_torus;

  InvolutionData d_involutionData;

  /*!
  \brief Fiber group.

  In terms of the complex torus $H$ and the Cartan involution \f$\tau\f$, it
  is equal to \f$F=H^{-\tau}/(1-\tau)H\f$: the group of fixed points of
  \f$\-tau\f$, modulo its identity component. Here additive notation is used
  for multiplicative composition; e.g., \f$1-\tau\f$ maps \f$z\in H\f$ to
  \f$z/\tau(z)\f$; this notation is natural when \f$\tau\f$ is viewed as
  operating on the weight lattice, in which case \f$1-\tau\f$ acts as the
  identity matrix minus the matrix giving \f$\tau\f$.

  Recall that the group of real points $H(R)$ is a product of $p$ unit circle
  groups, $q$ groups $R^\times$, and $r$ groups $C^\times$. It is easy to see
  that $F=(\Z/2\Z)^p$, with one factor $\Z/2\Z$ coming from each circle.

  The fiber group is represented as a subquotient of the group $H(2)$ of
  elements of order 2 (or 1) in the torus. Write $Y$ for the lattice of
  coweights of $H$; then one has natural isomorphisms $H(2)=(1/2)Y/Y=Y/2Y$.
  The Cartan involution \f$\tau\f$ is always represented by the matrix $q$ of
  its action on the character lattice $X$, which is dual to $Y$; so its action
  on $Y$ is given by the matrix $q^t$. The action of \f$-\tau\f$ on $Y$ is
  given by $-q^t$.

  Because $H$ is the product (not direct) of the connected groups
  \f$(H^\tau)_0\f$ and \f$(H^{-\tau})_0\f$, the group \f$H^\tau\f$ meets every
  component of \f$H^{-tau}\f$. The intersection of \f$H^\tau\f$ and
  \f$H^{-\tau}\f$ consists of \f$H(2)^\tau\f$, the \f$\tau\f$-fixed elements
  of order 2. It follows that the component group $F$ may be identified with
  \f$H(2)^\tau/(H(2)\cap (1-\tau)H)\f$. Thus $F$ can be represented by a
  subquotient of $H(2)=Y/2Y=(Z/2Z)^n$.

  The larger group \f$H(2)^\tau\f$ is computed as the kernel of the action of
  \f$(1+\tau)\f$ on $H(2)$, which via the isomorphism $H(2)=(Z/2Z)^n$ means
  the kernel of the reduction mod 2 of the matrix $I+q^t$. The smaller group
  \f$H(2)\cap (1-\tau)H\f$ consists of the elements of order 2 in the
  connected group \f$(H^{-\tau})_0\f$. It is computed as the reduction mod 2
  of the $-1$ eigenlattice \f$Y^{-\tau}\f$ (the latter is the lattice of
  coweights of \f$(H^{-tau})_0\f$).
  */
  latticetypes::SmallSubquotient d_fiberGroup;

  /*!
  \brief Fiber group for the adjoint group of G.

  Writing H_ad for the complex torus in G_ad, and still writing tau
  for the Cartan involution, this is F_ad=H_ad^{-tau}/(1-tau)H_ad. The
  adjoint covering G-->G_ad defines a natural map F-->F_ad, but this
  map may be neither injective nor surjective.

  The adjoint fiber is computed from the action of tau on the adjoint
  lattice of one-parameter subgroups Y_ad just as the fiber group is
  computed. The lattice Y_ad has as basis the fundamental coweights.
  */
  latticetypes::SmallSubquotient d_adjointFiberGroup;

  /*!
  \brief Flags the noncompact imaginary roots for the base grading
  among all the roots.
  */
  rootdata::RootSet d_baseNoncompact;

  /*!
  \brief Grading with all simple imaginary roots noncompact.

  The grading is a RankFlags; that is, a BitSet<RANK_MAX>.  It flags
  the noncompact imaginary roots among the simple imaginary roots; so
  all the bits up to imaginaryRank are 1.
  */
  gradings::Grading d_baseGrading;

  /*!
  \brief RootSet \#j flags the imaginary roots whose grading is
  changed by canonical basis vector \#j of the adjoint fiber group.
  */
  rootdata::RootSetList d_noncompactShift; // |size()==adjointFiberRank()|

  /*!
  \brief Grading \#j flags the simple imaginary roots whose
  grading is changed by canonical basis vector \#j in the adjoint
  fiber group.
  */
  gradings::GradingList d_gradingShift; // |size()==adjointFiberRank()|

  /*!
  \brief Matrix (over Z/2Z) of the map from the fiber group F for G to the
  fiber group F_ad for G_ad.

  Although the map is induced by Y/2Y-->Y_ad/2Y_ad, whose matrix is the
  reduction mod 2 of the matrix of Y-->Y_ad, that is _not_ the matrix in
  d_toAdjoint, which must take into account the subquotients at both sides.
  Each subquotient F and F_ad is equipped (by the row reduction algorithms in
  Subquotient) with a distinguished basis, and what is recorded in
  d_toAdjoint is the (dim F_ad) x (dim F) matrix with respect to those bases.
  */
  latticetypes::BinaryMap d_toAdjoint;

  /*!
  \brief Partition of the adjoint fiber group according to weak real
  forms.

  The imaginary Weyl group acts on the adjoint fiber group; the
  partition is by orbits of this action.  It is constructed by the
  function weakReal in the Helper class of the unnamed namespace of
  cartanclass.cpp.
  */
  partition::Partition d_weakReal;

  /*!
  \brief Partition of the weak real forms according to the
  corresponding classes in Z(G)^delta/[(1+delta)Z(G)].

  Constructed by the function realFormPartition in the Helper
  class of the unnamed namespace of cartanclass.cpp; details in that
  documentation.
  */
  partition::Partition d_realFormPartition;

  /*!
  \brief Partitions of Fiber group cosets corresponding to the
  possible square classes in Z^delta/[(1+delta)Z].

  The Fiber group acts in a simply transitive way on strong real forms
  (inducing tau on H) with a fixed square in Z^delta.  The number of
  squares that occur (modulo (1+delta)Z) is equal to the number c of
  classes in the partition d_weakReal.  The collection of strong real
  forms is therefore a collection of c cosets of the fiber group F.
  Each of these c cosets is partitioned into W_i-orbits; these partitions
  into orbits are described by the c partitions in d_strongReal.
  */
  std::vector<partition::Partition> d_strongReal;

  /*!
  \brief Representative strong real form for each weak real form.

  A StrongRealFormRep is a pair of numbers.  The second indexes
  the value of the square of the strong real form in
  Z^delta/[(1+delta)Z].  The first indexes a W_im orbit on the
  corresponding coset of the fiber group.
  */
  std::vector<StrongRealFormRep> d_strongRealFormReps;

 public:

// constructors and destructors

  // main and only constructor:
  Fiber(const rootdata::RootDatum&, const latticetypes::LatticeMatrix&);

  ~Fiber();

// copy and assignment

  Fiber(const Fiber&);

  Fiber& operator= (const Fiber&);

  void swap(Fiber&);

// accessors

/*!
  \brief Real torus defined over R.

  Represented as the lattice Z^n endowed with an involutive
  automorphism (represented by its n x n integer matrix).
*/
  const tori::RealTorus& torus() const {
    return *d_torus;
  }

/*!
  \brief RootSet flagging the complex roots.

  That is, a bitmap whose set bits are those corresponding to the
  numbers (within the list of roots in RootDatum) of the complex roots
  (those roots alpha with tau(alpha) equal neither to alpha nor to
  -alpha).
*/
  const rootdata::RootSet& complexRootSet() const {
    return d_involutionData.complex_roots();
  }

/*!
  \brief RootSet flagging the imaginary roots.

  That is, a bitmap whose set bits are those corresponding to the
  numbers (within the list of roots in RootDatum) of the imaginary
  roots (those roots alpha with tau(alpha)=alpha).
*/
  const rootdata::RootSet& imaginaryRootSet() const {
    return d_involutionData.imaginary_roots();
  }

/*!
  \brief RootSet flagging the real roots.

  That is, a bitmap whose set bits are those corresponding to the
  numbers (within the list of roots in RootDatum) of the real roots
  (those roots alpha with tau(alpha)=-alpha).
*/
  const rootdata::RootSet& realRootSet() const {
    return d_involutionData.real_roots();
  }

/*!
  \brief RootList holding the numbers of the simple imaginary roots.

  These are simple for the positive imaginary roots given by the
  (based) RootDatum.  They need not be simple in the entire root
  system.
*/
  const rootdata::RootList& simpleImaginary() const {
    return d_involutionData.imaginary_basis();
  }
  const rootdata::RootNbr simpleImaginary(size_t i) const {
    return d_involutionData.imaginary_basis()[i];
  }

  size_t imaginaryRank() const {
    return simpleImaginary().size();
  }

/*!
  \brief RootList holding the numbers of the simple real roots.

  These are simple for the positive real roots given by the
  (based) RootDatum.  They need not be simple in the entire root
  system.
*/
  const rootdata::RootList& simpleReal() const {
    return d_involutionData.real_basis();
  }
  const rootdata::RootNbr simpleReal(size_t i) const {
    return d_involutionData.real_basis()[i];
  }

  size_t realRank() const {
    return simpleReal().size();
  }

/*!
\brief Action of the Cartan involution on root \#j.
*/
  rootdata::RootNbr involution_image_of_root(rootdata::RootNbr j) const {
    return d_involutionData.root_involution()[j];
  }

/*!
  \brief Fiber group for the adjoint group of G.

  Writing H_ad for the complex torus in G_ad, and still writing tau
  for the Cartan involution, this is F_ad=H_ad^{-tau}/(1-tau)H_ad. The
  adjoint covering G-->G_ad defines a natural map F-->F_ad, but this
  map may be neither injective nor surjective.

  The adjoint fiber is computed from the action of tau on the adjoint
  lattice of one-parameter subgroups Y_ad just as the fiber group is
  computed. The lattice Y_ad has as basis the fundamental coweights.
*/
  const latticetypes::SmallSubquotient& adjointFiberGroup() const {
    return d_adjointFiberGroup;
  }

/*!
\brief Dimension of the adjoint fiber group as a Z/2Z vector space.
*/
  size_t adjointFiberRank() const {
    return d_adjointFiberGroup.dimension();
  }

/*!
\brief Cardinality of the adjoint fiber group.
*/
  size_t adjointFiberSize() const {
    return d_adjointFiberGroup.size();
  }


/*!
  \brief Fiber group.
*/
  const latticetypes::SmallSubquotient& fiberGroup() const {
    return d_fiberGroup;
  }

  /*!
  \brief Dimension of the fiber group as a Z/2Z vector space.
  */
  size_t fiberRank() const {
    return d_fiberGroup.dimension();
  }

  /*!
  \brief Cardinality of the fiber group: 2^dimension.
  */
  size_t fiberSize() const {
    return d_fiberGroup.size();
  }


  rootdata::RootSet compactRoots(AdjointFiberElt x) const;

  rootdata::RootSet noncompactRoots(AdjointFiberElt x) const;

  /*!
  \brief grading associated to an adjoint fiber element
  */
  gradings::Grading grading(AdjointFiberElt x) const;

  AdjointFiberElt gradingRep(const gradings::Grading&) const;

  const latticetypes::LatticeMatrix& involution() const;


  /*!
  \brief Image of a coroot (expressed in weight basis) in the fiber group
  */
  latticetypes::SmallBitVector mAlpha(const rootdata::Root&) const;

  size_t minusRank() const;

/*!
  \brief Number of weak real forms containing this Cartan.

  This is the number of orbits of the imaginary Weyl group on the
  adjoint fiber group.
*/
  size_t numRealForms() const {
    return d_weakReal.classCount();
  }

  size_t plusRank() const;

/*!
  \brief Partition of the weak real forms according to the
  corresponding central square classes in Z(G)/[(1+delta)Z(G)].

  A weak real form (always containing our fixed real torus) is an
  orbit of W_i on the adjoint fiber group.
*/
  const partition::Partition& realFormPartition() const {
    return d_realFormPartition;
  }

/*!
  \brief Returns the central square class to which a weak real form belongs.
*/
  square_class central_square_class (adjoint_fiber_orbit wrf) const
  {
    return square_class(d_realFormPartition(wrf));
  }

/*!
  \brief Returns the base element for a central square class
*/
  AdjointFiberElt class_base(square_class c) const
  {
    return d_weakReal.classRep(d_realFormPartition.classRep(c));
  }

/*!
 \brief Partitions of Fiber group cosets corresponding to the
  possible square classes in Z^delta/[(1+delta)Z].

  The Fiber group acts in a simply transitive way on strong real forms
  (inducing tau on H) with a fixed square in Z^delta.  The number of
  squares that occur (modulo (1+delta)Z) is equal to the number c of
  classes in the partition d_weakReal.  The collection of strong real
  forms is therefore a collection of c cosets of the fiber group F.
  Each of these c cosets is partitioned into W_i orbits; these orbits
  are described by the c partitions in d_strongReal.
*/
  const partition::Partition& strongReal(square_class j) const {
    return d_strongReal[j];
  }

/*!
  \brief Representative strong real form for real form \#rf.

  A StrongRealFormRep is a pair of numbers.  The second indexes
  the value of the square of the strong real form in
  Z^delta/[(1+delta)Z].  The first indexes a W_im orbit on the
  corresponding coset of the fiber group.
*/
  const StrongRealFormRep& strongRepresentative(adjoint_fiber_orbit wrf) const
  {
    return d_strongRealFormReps[wrf];
  }

/*!\brief Natural linear map from fiber group to adjoint fiber group

  This map is induced by the inclusion of the root lattice in the character
  lattice
*/
  AdjointFiberElt toAdjoint(FiberElt) const;

/*!\brief Returns the class number in the weak real form partition of the
  strong real form \#c in real form class \#rfc.
*/
  adjoint_fiber_orbit toWeakReal(fiber_orbit c, square_class csc) const;


/*!
  \brief Partition of the weak real forms according to the
  corresponding classes in Z(G)^delta/[(1+delta)Z(G)].
*/
  const partition::Partition& weakReal() const {
    return d_weakReal;
  }

// private accessors only needed during construction

private:


  latticetypes::SmallSubquotient makeFiberGroup() const;

  latticetypes::SmallSubquotient makeAdjointFiberGroup
    (const rootdata::RootDatum&) const;

  latticetypes::SmallSubspace gradingGroup(const rootdata::RootDatum&) const;

  latticetypes::LatticeMatrix adjointInvolution
    (const rootdata::RootDatum&) const;

  gradings::Grading makeBaseGrading
    (rootdata::RootSet& flagged_roots,const rootdata::RootDatum&) const;

  gradings::GradingList makeGradingShifts
    (rootdata::RootSetList& all_shifts,const rootdata::RootDatum&) const;

  bitset::RankFlagsList adjointMAlphas (const rootdata::RootDatum&) const;

  bitset::RankFlagsList mAlphas(const rootdata::RootDatum&) const;

  latticetypes::BinaryMap makeFiberMap(const rootdata::RootDatum&) const;

  partition::Partition makeWeakReal(const rootdata::RootDatum&) const;

  partition::Partition makeRealFormPartition() const;

  std::vector<partition::Partition> makeStrongReal
    (const rootdata::RootDatum& rd) const;

  std::vector<StrongRealFormRep> makeStrongRepresentatives() const;


}; // class Fiber





/*

			The |CartanClass| class

*/

/*!
\brief Represents a single stable conjugacy class of Cartan subgroups.

Mathematically this means the complex torus in the complex group, together
with an involutive automorphism of this torus (the Cartan involution). (More
precisely, it is a W-conjugacy class of involutive automorphisms.) As the
class is now used in the software, the Cartan involution must be in the inner
class specified by ComplexReductiveGroup. Most of the interesting information
is contained in the two underlying Fiber classes d_fiber and d_dualFiber.
First of all, that is where the Cartan involution lives (since the involution
is needed to define the fibers). But the fiber classes also record for which
real forms this stable conjugacy class of Cartan subgroups is defined; and,
given the real form, which imaginary roots are compact and noncompact.
*/
class CartanClass {

 private:

  /*!
  \brief Class of the fiber group H^{-tau}/[(1-tau)H] for this
  Cartan.

  Elements (very roughly) correspond to possible extensions
  of the real form tau from H to G.
  */
  Fiber d_fiber;

  /*!
  \brief Class of the fiber group for the dual Cartan.

  Elements of the dual fiber group are characters of the group of
  connected components of the real points of H.
  */
  Fiber d_dualFiber;

  /*!
  \brief Roots simple for the "complex factor" of W^tau.

  The subgroup \f$W^\tau\f$ of Weyl group elements commuting with the Cartan
  involution \f$\tau\f$ has two obvious (commuting) normal subgroups: the Weyl
  group $W^R$ of the real (that is, fixed by \f$-\tau\f$) roots, and the Weyl
  group $W^{iR}$ of the imaginary (that is, fixed by \f$\tau\f$) roots. But
  this is not all of \f$W^\tau\f$, as is easily seen in the case of complex
  groups, where both $W^R$ and $W^{iR}$ are trivial (there are no real or
  imaginary roots), yet $W$ is a direct sum of two identical factors
  interchanged by \f$\tau\f$, and the actions of diagonal elements of that sum
  clearly commute with \f$\tau\f$. In general there is a group denoted $W^C
  such that \f$W^\tau\f$ is a semidirect product of $W^R * W^{iR}$ (the normal
  subgroup) with $W^C$. Here is how to describe $W^C$.

  Write $RC$ (standing for "complex roots") for the roots orthogonal to both
  the sum of positive real roots, and the sum of positive imaginary roots
  (this group depends on the choice of positive roots, but all choices lead to
  conjugate subgroups). It turns out that $RC$ as a root system is the direct
  sum of two isomorphic root systems $RC_0$ and $RC_1$ interchanged by $\tau$.
  (There is no canonical choice of this decomposition.) Now $W^C$ is the set
  of $\tau$-fixed elements of $W(RC_0) \times W(RC_1)$, in other words its
  diagonal subgroup.

  In general $W^C$ is not the Weyl group of a root subsystem, but it is
  isomorphic the the Weyl group of (any choice of) $RC_0$ (or of $RC_1$). We
  make a choice for $RC_0$, and |d_simpleComplex| lists its simple roots, so
  that $W^C$ is isomorphic to the Weyl group generated by the reflections
  corresponding to the roots whose numbers are in |d_simpleComplex|.
  */
  rootdata::RootList d_simpleComplex;

  /*!
  \brief Size of the W-conjugacy class of tau.

  The number of distinct involutions defining the same stable
  conjugacy class of Cartan subgroups.
  */
  size::Size d_orbitSize;

public:

// constructors and destructors

  CartanClass(const rootdata::RootDatum&, const latticetypes::LatticeMatrix&);

  ~CartanClass()
    {}

// copy and assignment: defaults are ok for copy and assignment

  void swap(CartanClass&);

// accessors

  /*!
  \brief Returns the matrix of the involution on the weight lattice of
  the Cartan subgroup.
  */
  const latticetypes::LatticeMatrix& involution() const {
    return d_fiber.involution();
  }

  /*!
  \brief Action of the Cartan involution on root \#j.
  */
  rootdata::RootNbr involution_image_of_root(rootdata::RootNbr j) const {
    return d_fiber.involution_image_of_root(j);
  }

  /*!
  \brief RootSet flagging the imaginary roots.

  That is, a bitmap whose set bits are those corresponding to the
  numbers (within the list of roots in RootDatum) of the imaginary
  roots (those roots \f$\alpha\f$ with \f$\tau(\alpha)=\alpha\f$).
  */
  const rootdata::RootSet& imaginaryRootSet() const {
    return d_fiber.imaginaryRootSet();
  }
  /*!
  \brief RootSet flagging the real roots.

  That is, a bitmap whose set bits are those corresponding to the
  numbers (within the list of roots in RootDatum) of the real roots
  (those roots \f$\alpha\f$ with \f$\tau(\alpha)=-\alpha\f$).
  */
  const rootdata::RootSet& realRootSet() const {
    return d_fiber.realRootSet();
  }

  /*!
  \brief RootList holding the numbers of the simple imaginary roots.

  These are simple for the subsystem of imaginary roots. They need not be
  simple in the entire root system.
  */
  const rootdata::RootList& simpleImaginary() const {
    return d_fiber.simpleImaginary();
  }

  /*!
  \brief RootList holding the numbers of the simple real roots.

  These are simple for subsystem of real roots. They need not be simple in the
  entire root system.
  */
  const rootdata::RootList& simpleReal() const {
    return d_fiber.simpleReal();
  }
/* Since only the _numbers_ of simple-real roots are returned, this used to be
   obtained as |d_dualFiber.simpleImaginary()|, before the |InvolutionData|
   contained the basis for the real root system as well. That is indeed the
   same value, given that the constructor for a dual root system preserves the
   numbering of the roots (but exchanging roots and coroots of course). Such a
   root datum differs however from one that could be constructed directly
   (because its coroots rather than its roots are sorted lexicographically);
   as the construction of root data might change in the future, it seems safer
   and more robust to just take the basis of the real root subsystem. MvL.
*/


  /*!
  \brief Class of the fiber group H^{-tau}/[(1-tau)H] for this
  Cartan.

  Elements (very roughly) correspond to possible extensions
  of the real form tau from H to G.
  */
  const Fiber& fiber() const {
    return d_fiber;
  }

  /*!
  \brief Class of the fiber group for the dual Cartan.

  Elements of the dual fiber group are characters of the group of
  connected components of the real points of H.
  */
  const Fiber& dualFiber() const {
    return d_dualFiber;
  }


  bool isMostSplit(adjoint_fiber_orbit wrf) const;

  /*!
  \brief Number of weak real forms for the dual group containing the
  dual Cartan.
  */
  size_t numDualRealForms() const {
    return d_dualFiber.numRealForms();
  }
  /*!
  \brief Number of weak real forms containing this Cartan.

  This is the number of orbits of the imaginary Weyl group on the
  adjoint fiber group.
  */
  size_t numRealForms() const {
    return d_fiber.numRealForms();
  }

  /*!
  \brief Number of possible squares of strong real forms mod
  (1+delta)Z.

  This is the number of classes in the partition of weak real forms
  according to Z^delta/[(1+delta)Z].
  */
  size_t numRealFormClasses() const {
    return d_fiber.realFormPartition().classCount();
  }

  /*!
  \brief Size of the W-conjugacy class of tau.

  The number of distinct involutions defining the same stable conjugacy
  class of Cartan subgroups.
   */
   size_t orbitSize() const {
    return d_orbitSize.toUlong();
  }

  /*!
  \brief Roots simple for the "complex factor" of W^tau.

  The subgroup W^tau of Weyl group elements commuting with the Cartan
  involution tau has two obvious (commuting) normal subgroups: the
  Weyl group W^R of the real (that is, fixed by -tau) roots, and the
  Weyl group W^iR of the imaginary (that is, fixed by tau) roots. But
  this is not all of W^: W^tau is a semidirect product of (W^R x W^iR)
  with a group W^C, the first factor being normal. Here is how to
  describe W^C.

  Write RC (standing for "complex roots") for the roots orthogonal to
  (a) the sum of positive real roots, and also (b) the sum of positive
  imaginary roots. It turns out that RC as a root system is the direct
  sum of two isomorphic root systems RC_0 and RC_1 interchanged by
  tau. (There is no canonical choice of this decomposition.) The group
  W^tau includes W^C, the diagonal subgroup of W(RC_0) x W(RC_1). The
  list d_simpleComplex is the numbers of the simple roots for (a
  choice of) RC_0. That is, the Weyl group W^C is isomorphic to the
  Weyl group generated by the reflections corresponding to the numbers
  in d_simpleComplex.
   */
  const rootdata::RootList& simpleComplex() const {
    return d_simpleComplex;
  }

  /*!
  \brief Partitions of Fiber group cosets corresponding to the
  possible square classes in Z^delta/[(1+delta)Z].

  The Fiber group acts in a simply transitive way on strong real forms
  (inducing tau on H) with a fixed square in Z^delta. The number of
  squares that occur (modulo (1+delta)Z) is equal to the number c of
  classes in the partition d_weakReal. The collection of strong real
  forms is therefore a collection of c cosets of the fiber group
  F. Each of these c cosets is partitioned into W_i orbits; these
  orbits are described by the c partitions in d_strongReal.
  */
  const partition::Partition& strongReal(square_class j) const {
    return d_fiber.strongReal(j);
  }

  /*!
  \brief Returns the image of x in the adjoint fiber group.

  Precondition: x is a valid element in the fiber group.
  */
  AdjointFiberElt toAdjoint(FiberElt x) const {
    return d_fiber.toAdjoint(x);
  }

  /*!
  \brief Returns the class number in the weak real form partition of
  the strong real form \#c in real form class rfc.

  The pair (c,rfc) is the software representation of an equivalence
  class of strong real forms (always assumed to induce tau on H). The
  integer rfc labels an element of Z^delta/[(1+delta)Z], thought of as
  a possible square value for strong real forms. The fiber group acts
  simply transitively on strong real forms with square equal to
  rfc. The integer c labels an orbit of W_i on this fiber group coset;
  this orbit is the equivalence class of strong real forms.
  */
  adjoint_fiber_orbit toWeakReal(fiber_orbit c, square_class j) const {
    return d_fiber.toWeakReal(c,j);
  }

  /*!
  \brief Partition of the weak real forms according to the
  corresponding classes in Z(G)^delta/[(1+delta)Z(G)].
  */
  const partition::Partition& weakReal() const {
    return d_fiber.weakReal();
  }


// private accessors only needed during construction

private:
  rootdata::RootList makeSimpleComplex(const rootdata::RootDatum&) const;

  size::Size makeOrbitSize(const rootdata::RootDatum&) const;

}; // class CartanClass

}

}

#endif
