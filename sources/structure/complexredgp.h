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

#include "bitmap.h"
#include "tags.h"
#include "setutils.h"

#include "cartanclass.h"
#include "latticetypes.h"
#include "poset.h"
#include "rootdata.h"
#include "tits.h"
#include "weyl.h"

#include "lietype.h"
#include "realform.h"

/******** function declarations **********************************************/

namespace atlas {

namespace complexredgp {

  void lieType(lietype::LieType&, const ComplexReductiveGroup&);

/* this function should NOT be made into a method, suppressing the |rd| and |W|
   parameters, as these can be be dual to those in the |ComplexReductiveGroup|!
*/
  void Cayley_and_cross_part(rootdata::RootList& cayley,
			     weyl::WeylWord& cross,
			     const weyl::TwistedInvolution& tw,
			     const rootdata::RootDatum& rd,
			     const weyl::TwistedWeylGroup& W);

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

  The field |d_rootDatum| stores the root datum, which must have been
  constructed before. The field |d_titsGroup| holds the mentioned (enlarged)
  Tits group, which is constructed upon entry from the root datum and the
  involution; it also gives access to just the (complex) Weyl group when that
  is necessary. Finally the other fields store the information relative to
  (stable conjugacy classes of) Cartan subgroups and real forms.

  Each stable conjugacy classes of Cartan subgroups corresponds to a
  W-conjugacy class of involutions in the Gamma-enlarged Weyl group (W
  semidirect <Gamma>, where <Gamma>=Z/2Z acts on W), contained in the
  complement of its subgroup W. Since such involutions are of the form
  (w,Gamma), they can be represented by their element w, which is called a
  twisted involution. The condition for being a twisted involution $t$ is
  $t\Gamma(t)=e$ and "twisted conjugacy" of $t$ by \f$w\in W\f$ is given by
  \f$w\cdot t=wt\Gamma(w^{-1})\f$. The stable conjugacy classes of Cartan
  subgroups will each be represented by a canonical representative of the
  corresponding twisted conjugacy class of twisted involutions.

  In addition to describing the set of Cartan classes, this class provides
  access (via the |d_cartan| array) to data for each individual one of them,
  and (via |Cartan_poset|) to the partial order relation between them. For the
  latter, let |tau_i| be involutions acting on the complex torus |H| for
  various classes of Cartan subgroups; (H,tau_1) is considered "more compact"
  than (H,tau_2) if the identity component of the fixed point set H^tau_2 is
  W-conjugate to a subtorus of H^tau_1.

  The problem for the dual group of G is identical, the bijection taking the
  negative transpose of a twisted involution. This bijection reverses the
  partial order on Cartans. The class also provides access to Cartans in the
  dual group.

  */
class ComplexReductiveGroup
{
  /*!
  \brief The based root datum.
  */
  const rootdata::RootDatum d_rootDatum;

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

  const weyl::WeylGroup W;

  /*!
  \brief The Tits group of the based root datum, extended by an
  involutive automorphism.
  */
  const tits::TitsGroup d_titsGroup;

  //!\brief the permutation of the roots given by the based automorphism
  const setutils::Permutation root_twist;

  /*!\brief
  List of owned pointers of stable conjugacy classes of Cartan subgroups.

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
  weyl::TwistedInvolutionList d_twistedInvolution;


  /*!
  \brief Entry \#rf is the number of the most split Cartan for real form \#rf.
  */
  std::vector<size_t> d_mostSplit;

  /*!
  \brief Partial order of Cartan subgroups.

  This is the ordering by containment of H^theta up to conjugacy: (H,theta_1)
  precedes (H,theta_2) if (H^theta_2)_0 is W-conjugate to a subtorus of
  H^theta_1. Numbering of elements is as in d_twistedInvolution
  */
  poset::Poset Cartan_poset;

  /*!
  \brief Entry \#rf flags the Cartans defined in real form \#rf.
  */
  std::vector<bitmap::BitMap> d_support;

  /*!
  \brief Entry \#rf flags the dual Cartans defined in dual real form \#rf.
  */
  std::vector<bitmap::BitMap> d_dualSupport;

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
  \brief  Flags the set of real forms for which the full set of Cartan
  classes is constructed.

  Because of the way the construction proceeds, these are exactly the
  real forms for which the most split Cartan has been reached.
  */
  bitmap::BitMap d_status;

// copy, assignement and swap are forbidden, and should not be implemented
  ComplexReductiveGroup(const ComplexReductiveGroup&);
  ComplexReductiveGroup& operator= (const ComplexReductiveGroup&);
  void swap(ComplexReductiveGroup& G);

 public:
// constructors and destructors
  ComplexReductiveGroup(const rootdata::RootDatum&,
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
  \brief returns the rank of the group.
*/
  size_t rank() const { return d_rootDatum.rank(); }

//!\brief returns the semisimple rank of the group.
  size_t semisimpleRank() const { return d_rootDatum.semisimpleRank(); }

  const rootdata::RootDatum& rootDatum() const { return d_rootDatum; }

  const weyl::WeylGroup& weylGroup() const { return W; }
  const weyl::TwistedWeylGroup& twistedWeylGroup() const
    { return d_titsGroup; } // in fact its base object

  setutils::Permutation simple_twist() const
    { return setutils::Permutation
	(&twistedWeylGroup().twist()[0],
	 &twistedWeylGroup().twist()[semisimpleRank()]); }

  //!\brief returns a reference to the Weyl group (owned by the Tits group).
  const tits::TitsGroup& titsGroup() const { return d_titsGroup; }

/*!
  \brief Recover the matrix of the involution for the fundamental Cartan.

  This is the one permuting the simple roots, the distinguished one among the
  involutions in this inner class of G.
*/
  const latticetypes::LatticeMatrix& distinguished() const
    { return fundamental().involution(); }

/*!\brief
  Matrix of involution for the fundamental Cartan in the dual group.

  This is -w_0 times the transpose of the fundamental involution.
*/
  const latticetypes::LatticeMatrix& dualDistinguished() const
    { return dualFundamental().involution(); }


//!\brief Returns data for stable conjugacy class \#cn of Cartan subgroups.
  const cartanclass::CartanClass& cartan(size_t cn) const
    { return *d_cartan[cn]; }

/*!
  \brief returns the set of Cartan classes for |rf|
  Requires those Cartan to be already generated
*/
  const bitmap::BitMap& Cartan_set(realform::RealForm rf) const
    { return d_support[rf]; }

/*!
  \brief returns the set of Cartan classes for the dual real form rf
  Only those Cartan classes already generated will show up
*/
  const bitmap::BitMap& dual_Cartan_set(realform::RealForm rf) const
    { return d_dualSupport[rf]; }


//!\brief Returns the partial ordering of (currently generated) Cartan classes
  const poset::Poset& Cartan_ordering() const
    { return Cartan_poset; }

//!\brief matrix giving involution action of |tw| on weight lattice
  latticetypes::LatticeMatrix
    involutionMatrix(const weyl::TwistedInvolution& tw) const;



/*!
  \brief returns the fundamental fiber.

  This is a technical data structure containing the data for the classification
  of weak and strong real forms of G.
*/
  const cartanclass::Fiber& fundamental() const
    { return d_fundamental; }

/*!
  \brief returns the dual fundamental fiber.

  This is a technical data structure containing the data for the classification
  of weak and strong dual real forms of G.
*/
  const cartanclass::Fiber& dualFundamental() const
    { return d_dualFundamental; }

/*!
  \brief returns the number of conjugacy classes of Cartan subgroups
  currently constructed for G. Only after fillCartan() has been called is it
  ensured that this gives the total number of Cartan subgroups for this inner
  class.
*/
  size_t numCartanClasses() const
    { return d_cartan.size(); }

/*!
  \brief returns the number of involutions for the currently defined Cartans.
*/
  size_t numInvolutions() const;

/*!
  \brief returns the number of involutions for the indicated Cartans.
*/
  size_t numInvolutions(const bitmap::BitMap& Cartan_classes) const;

//!\brief Returns the number of weak real forms of G.
  size_t numRealForms() const { return fundamental().numRealForms(); }

/*!\brief
  Returns the number of weak real forms of G for which Cartan \#cn is defined.
*/
  size_t numRealForms(size_t cn) const { return cartan(cn).numRealForms(); }


//!\brief Returns the number of weak real forms of the dual group of G.
  size_t numDualRealForms() const { return dualFundamental().numRealForms(); }

/*!\brief
  Returns the number of weak real forms of the dual group for which the dual
  of Cartan \#cn is defined.
*/
  size_t numDualRealForms(size_t cn) const
    { return cartan(cn).numDualRealForms(); }


/*!\brief
  The size of the fiber orbits corresponding to strong real forms lying
  over weak real form \#rf, in cartan \#cn (all orbits have the same size)

  Explanation: this is the size of the orbits, for the shifted action of
  the imaginary Weyl group, that correspond to rf in the classification of
  strong real forms (they all have the same size.)

  Precondition: Real form \#rf is defined for cartan \#cn, and the latter has
  been generated (obviously, otherwise it's number would not be known)
*/
  unsigned long fiberSize(realform::RealForm rf, size_t cn) const;

  unsigned long dualFiberSize(realform::RealForm, size_t cn) const;

/*!
  \brief returns the number of elements in K\\G/B for real form \#rf.

  Precondition: the Cartan classes for this real form have been generated

  Explanation: this is exactly the number of elements in the one-sided
  parameter set corresponding to any strong real form of G lying over rf.
*/
  unsigned long KGB_size(realform::RealForm rf,
			 const bitmap::BitMap& Cartan_classes) const;

  unsigned long KGB_size(realform::RealForm rf) const
    { return KGB_size(rf,Cartan_set(rf)); }


/*!
  \brief the size of the block defined by the weak real form rf
  and the weak dual real form drf. Requires Cartan for |rf| to be generated
*/
  unsigned long block_size(realform::RealForm, realform::RealForm,
			   const bitmap::BitMap& Cartan_classes) const;

  unsigned long block_size(realform::RealForm rf, realform::RealForm drf)
    const
    { return block_size(rf,drf,Cartan_set(rf)& dual_Cartan_set(drf)); }

//!\brief Returns the (inner) number of the quasisplit real form.
  realform::RealForm quasisplit() const { return realform::RealForm(0); }



/*!
  \brief Entry \#|rf| gives the most split Cartan for real form \#|rf|.

  Precondition: |fillCartan()| has been called, or at least |fillCartan(rf)|
*/
  size_t mostSplit(realform::RealForm rf) const
    { return d_mostSplit[rf]; }

/*!
  \brief returns the real form labels for cartan \#cn

  More precisely, realFormLabels(cn)[i] is the (inner) number of the real form
  that corresponds to part i of the partition cartan(n).fiber().weakReal()
*/
  const realform::RealFormList& realFormLabels(size_t cn) const
    { return d_realFormLabels[cn]; }

/*!
  \brief returns the dual real form labels for cartan \#cn
*/
  const realform::RealFormList& dualRealFormLabels(size_t cn) const
    { return d_dualRealFormLabels[cn]; }


/*! get part in the |weakReal| partition of the fiber in Cartan \#cn
    corresponding to real form |rf|
 */
  cartanclass::adjoint_fiber_orbit
    real_form_part(realform::RealForm rf, size_t cn) const
    { return setutils::find_index(d_realFormLabels[cn],rf); }

/*! get part in the |weakReal| partition of the dual fiber in Cartan \#cn
    corresponding to dual real form |drf|
 */
  cartanclass::adjoint_fiber_orbit
    dual_real_form_part(realform::RealForm drf, size_t cn) const
    { return setutils::find_index(d_dualRealFormLabels[cn],drf); }

/*!\brief
  An element of the orbit in the adjoint fiber corresponding to |rf|
  in the classification of weak real forms for cartan |\#cn|.
  This amounts to searching for |rf| in |d_realFormLabels[cn]|.

  Precondition: cartan \#|cn| is defined for |rf|.
*/
  unsigned long representative(realform::RealForm rf, size_t cn) const
    { return cartan(cn).fiber().weakReal().classRep(real_form_part(rf,cn)); }

/*!\brief
  An element of the orbit in the adjoint dual fiber corresponding to |drf|
  in the classification of dual weak real forms for cartan |\#cn|.
  This amounts to searching for |drf| in |d_dualRealFormLabels[cn]|.

  Precondition: cartan \#|cn| is defined for |drf|.
*/
  unsigned long dualRepresentative(realform::RealForm drf, size_t cn) const
    { return cartan(cn).dualFiber().
	weakReal().classRep(dual_real_form_part(drf,cn));
    }


  const setutils::Permutation& root_involution() const { return root_twist; }

  rootdata::RootNbr twisted_root(rootdata::RootNbr alpha) const
    { return root_twist[alpha]; }

/*! \brief Make |sigma| canonical and return Weyl group |w| element that
    twisted conjugates the canonical representative back to |sigma|
*/
  const weyl::WeylElt
    canonicalize(weyl::TwistedInvolution& sigma, bitset::RankFlags gens) const;

  inline
  const weyl::WeylElt canonicalize(weyl::TwistedInvolution& sigma) const
    { return canonicalize
	(sigma,bitset::RankFlags(constants::lMask[semisimpleRank()]));
    }


/*!\brief
  (Representative) twisted involutions for each class of Cartan subgroup.
*/
  const weyl::TwistedInvolution& twistedInvolution(size_t cn) const
    { return d_twistedInvolution[cn]; }

  size_t class_number(weyl::TwistedInvolution) const;

/*!\brief
  Returns the set of noncompact imaginary roots for (the representative
  in the adjoint fiber of) the real form \#rf.
*/
  rootdata::RootSet noncompactRoots(realform::RealForm rf) const
  {
    return
      fundamental().noncompactRoots(fundamental().weakReal().classRep(rf));
  }

  rootdata::RootSet noncompactPosRootSet(realform::RealForm, size_t) const;

  latticetypes::LatticeElt
    posRealRootSum(const weyl::TwistedInvolution&) const;

  latticetypes::LatticeElt
    posImaginaryRootSum(const weyl::TwistedInvolution&) const;


//!\brief apply involution action of |tw| on weight lattice to |v|
  void twisted_act
    (const weyl::TwistedInvolution& tw,latticetypes::LatticeElt& v) const;

// manipulators

/*!
  \brief fills in the Cartan classes that are defined for the real form |rf|.
*/
void fillCartan(realform::RealForm rf) { extend(rf); }

/* the following is not done via a default argument to the previous method
   since a default argument cannot refer to a class member (quasisplit) */
 void fillCartan() { fillCartan(quasisplit()); }


 private:
// auxiliary accessors
weyl::TwistedInvolution
  reflection(rootdata::RootNbr rn,const weyl::TwistedInvolution& tw) const;

//!\brief Tells whether Cartan \#cn is defined in real form \#rf.
bool is_defined(realform::RealForm rf, size_t cn) const
  { return Cartan_set(rf).isMember(cn); }

// auxiliary manipulators

void extend(realform::RealForm);

void addCartan(weyl::TwistedInvolution tw)
  { d_cartan.push_back
      (new cartanclass::CartanClass(rootDatum(),involutionMatrix(tw)));
  }


void correlateForms();

void correlateDualForms(const rootdata::RootDatum& dual_rd,
			const weyl::TwistedWeylGroup& dual_W);

void updateSupports(size_t last_Cartan_class_added);

void updateStatus(size_t prev_Cartan_size);

}; // |class ComplexReductiveGroup|

} // |namespace complexredgp|

} // |namespace atlas|

#endif
