/*!
\file
\brief Class definitions and function declarations for the class
ComplexReductiveGroup.
*/
/*
  This is complexredgp.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006--2010 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef COMPLEXREDGP_H  /* guard against multiple inclusions */
#define COMPLEXREDGP_H

#include "atlas_types.h"

#include "tags.h"
#include "bitmap.h"	// containment of bitmaps for real forms
#include "permutations.h"// containment of root twist

#include "cartanclass.h"// containment of |Fiber|
#include "involutions.h"// containment of |InvolutionTable|, |Cartan_orbits|
#include "poset.h"	// containment of Cartan poset
#include "rootdata.h"	// containment of root datum and its dual
#include "tits.h"	// containment of Tits group and its dual

/******** function declarations **********************************************/

namespace atlas {

namespace complexredgp {

  WeylElt canonicalize // return value is conjugating element
    (TwistedInvolution& sigma,
     const RootDatum& rd,
     const TwistedWeylGroup& W,
     RankFlags gens);

/* this function should NOT be made into a method, suppressing the |rs| and |W|
   parameters, as these can be be dual to those in the |ComplexReductiveGroup|!
*/
  void Cayley_and_cross_part(RootNbrSet& Cayley,
			     WeylWord& cross,
			     const TwistedInvolution& tw,
			     const RootSystem& rs,
			     const TwistedWeylGroup& W);

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
  // The based root datum. It is stored here (constructed by our constructor)
  const RootDatum d_rootDatum;

  // The dual based root datum. It is also stored and constructed by us
  const RootDatum d_dualRootDatum;

  /*!
  \brief Fiber class for the fundamental Cartan subgroup.

  The involution is delta, which is stored here. It permutes the simple roots.
  */
  Fiber d_fundamental;

  /*!
  \brief Fiber class for the fundamental Cartan in the dual group.

  The fiber group here is the group of characters (i.e., the dual group)
  of the component group of the quasisplit Cartan.
  */
  Fiber d_dualFundamental;

  const WeylGroup* my_W; // pointer to |W| in case we own |W|, or |NULL|
  const WeylGroup& W;    // possibly owned (via |my_W|) reference

  /*!
  \brief The Tits group of the based root datum, extended by an
  involutive automorphism.
  */
  const TitsGroup d_titsGroup;

  /*!
  \brief The Tits group of the dual based root datum, extended by an
  involutive automorphism.
  */
  const TitsGroup d_dualTitsGroup;

  //!\brief the permutation of the roots given by the based automorphism
  const Permutation root_twist;

  typedef std::vector<RankFlags> form_reps; // gradings for real forms

  struct C_info
  { TwistedInvolution tw;
    BitMap real_forms,dual_real_forms; // mark present (dual) real forms
    form_reps rep,dual_rep; // gradings representing those (dual) real forms
    BitMap below;
    CartanClass* class_pt; //!< owned (by parent) pointer, might be NULL
    // remaining fields are set only once |class_pt| has been made non-NULL
    RealFormNbrList real_labels,dual_real_labels;

    C_info(const ComplexReductiveGroup& G, const TwistedInvolution twi,
	   CartanNbr i); // index |i| of this Cartan used to dimension |below|

  };

  std::vector<C_info> Cartan;

  /*!
  \brief Partial order of Cartan subgroups.

  This is the ordering by containment of H^theta up to conjugacy: (H,theta_1)
  precedes (H,theta_2) if (H^theta_2)_0 is W-conjugate to a subtorus of
  H^theta_1. Numbering of elements is as in the vector |Cartan|
  */
  poset::Poset Cartan_poset;

  //! Entry \#rf is the number of the most split Cartan for real form \#rf.
  std::vector<CartanNbr> d_mostSplit;

  Cartan_orbits C_orb;

// copy, assignement and swap are forbidden, and should not be implemented
  ComplexReductiveGroup(const ComplexReductiveGroup&);
  ComplexReductiveGroup& operator= (const ComplexReductiveGroup&);
  void swap(ComplexReductiveGroup& G);

 public:
// constructors and destructors
  ComplexReductiveGroup(const PreRootDatum&, // constructor builds root datum
			const WeightInvolution&);

  ComplexReductiveGroup(const RootDatum&, // alternative that copies root datum
			const WeightInvolution&);

  ComplexReductiveGroup(const ComplexReductiveGroup&, tags::DualTag);

  ~ComplexReductiveGroup();

// accessors


  const RootDatum& rootDatum() const { return d_rootDatum; }
  const RootDatum& dualRootDatum() const { return d_dualRootDatum; }
  const RootSystem& rootSystem() const {return d_rootDatum; } // base object
  const RootSystem& dualRootSystem() const {return d_dualRootDatum; } // base

  const Cartan_orbits& involution_table () const { return C_orb; }


/*!
  \brief returns the rank of the group.
*/
  size_t rank() const { return rootDatum().rank(); }

//!\brief returns the semisimple rank of the group.
  size_t semisimpleRank() const { return rootSystem().rank(); }

  const WeylGroup& weylGroup() const { return W; }
  const TwistedWeylGroup& twistedWeylGroup() const
    { return d_titsGroup; } // in fact its base object
  const TwistedWeylGroup& dualTwistedWeylGroup() const
    { return d_dualTitsGroup; } // in fact its base object

  Permutation simple_twist() const
    { return Permutation
	(&twistedWeylGroup().twist()[0],
	 &twistedWeylGroup().twist()[semisimpleRank()]); }

  //!\brief returns a reference to the Tits group.
  const TitsGroup& titsGroup() const { return d_titsGroup; }

  //!\brief returns a reference to the dual Tits group.
  const TitsGroup& dualTitsGroup() const { return d_dualTitsGroup; }

/*!
  \brief Recover the matrix of the involution for the fundamental Cartan.

  This is the one permuting the simple roots, the distinguished one among the
  involutions in this inner class of G.
*/
  const WeightInvolution& distinguished() const
    { return fundamental().involution(); }

/*!\brief
  Matrix of involution for the fundamental Cartan in the dual group.

  This is -w_0 times the transpose of the fundamental involution.
*/
  const CoweightInvolution& dualDistinguished() const
    { return dualFundamental().involution(); }


//!\brief returns the set of Cartan classes for the real form |rf|
  BitMap Cartan_set(RealFormNbr rf) const;

//!\brief returns the set of Cartan classes for the dual real form |drf|
  BitMap dual_Cartan_set(RealFormNbr drf) const;


//!\brief Returns the partial ordering of Cartan classes
  const poset::Poset& Cartan_ordering() const { return Cartan_poset; }

//!\brief matrix giving involution action of |tw| on weight lattice
  WeightInvolution involutionMatrix(const TwistedInvolution& tw) const;

  InvolutionData involution_data(const TwistedInvolution& tw) const
  { return InvolutionData::build(rootSystem(),twistedWeylGroup(),tw); }

/*!
  \brief returns the fundamental fiber.

  This is a technical data structure containing the data for the classification
  of weak and strong real forms of G.
*/
  const Fiber& fundamental() const { return d_fundamental; }

/*!
  \brief returns the dual fundamental fiber.

  This is a technical data structure containing the data for the classification
  of weak and strong dual real forms of G.
*/
  const Fiber& dualFundamental() const { return d_dualFundamental; }

//! \brief returns the number of conjugacy classes of Cartan subgroups.
  CartanNbr numCartanClasses() const { return Cartan.size(); }


//!\brief Returns the number of weak real forms of G.
  RealFormNbr numRealForms() const { return fundamental().numRealForms(); }

//!\brief Returns the number of weak real forms of the dual group of G.
  RealFormNbr numDualRealForms() const
  { return dualFundamental().numRealForms(); }

/*!\brief
  Returns the number of weak real forms of G for which Cartan \#cn is defined.
*/
  RealFormNbr numRealForms(CartanNbr cn) const
  { return Cartan[cn].real_forms.size(); }

/*!\brief
  Returns the number of weak real forms of the dual group for which the dual
  of Cartan \#cn is defined.
*/
  RealFormNbr numDualRealForms(CartanNbr cn) const
  { return Cartan[cn].dual_real_forms.size(); }

  const BitMap& real_forms(CartanNbr cn) const
  { return Cartan[cn].real_forms; }
  const BitMap& dual_real_forms(CartanNbr cn) const
  { return Cartan[cn].dual_real_forms; }

//!\brief Returns the (inner) number of the quasisplit real form.
  RealFormNbr quasisplit() const { return RealFormNbr(0); }

/*!
  \brief Entry \#|rf| gives the most split Cartan for real form \#|rf|.
*/
  CartanNbr mostSplit(RealFormNbr rf) const
    { return d_mostSplit[rf]; }

  const Permutation& root_involution() const { return root_twist; }

  RootNbr twisted_root(RootNbr alpha) const
    { return root_twist[alpha]; }

/*! \brief Make |sigma| canonical and return Weyl group |w| element that
    twisted conjugates the canonical representative back to |sigma|. Thus
    after the call the involution corresponding to |sigma| satisfies:

  (1) the sum of the positive real roots is dominant (call it $SR$)
  (2) in the subsystem of roots orthogonal to $SR$ (which contains all the
      imaginary roots), the sum of the imaginary roots (call it $SI$) is
      dominant for the subsystem
  (3) in the subsystem of roots orthogonal to both $SR$ and $SI$, the
      involution corresponding to twisted involution fixes (globally) the
      dominant chamber of the subsystem (it permutes its simple roots).
*/
  WeylElt
    canonicalize(TwistedInvolution& sigma, RankFlags gens) const;

  inline
  WeylElt canonicalize(TwistedInvolution& sigma) const
    { return canonicalize
	(sigma,RankFlags(constants::lMask[semisimpleRank()]));
    }


/*!\brief
  (Representative) twisted involutions for each class of Cartan subgroup.
*/
  const TwistedInvolution& twistedInvolution(CartanNbr cn) const
    { return Cartan[cn].tw; }
  TwistedInvolution dualTwistedInvolution(CartanNbr cn) const
    { return W.opposite(Cartan[cn].tw); }

  CartanNbr class_number(TwistedInvolution) const;

/*!\brief
  Returns the set of noncompact imaginary roots at the distinguished involution
  for (the representative in the adjoint fiber of) the real form \#rf.
*/
  RootNbrSet noncompactRoots(RealFormNbr rf) const
  {
    return
      fundamental().noncompactRoots(fundamental().weakReal().classRep(rf));
  }

  Weight
    posRealRootSum(const TwistedInvolution&) const;

  Weight
    posImaginaryRootSum(const TwistedInvolution&) const;

//!\brief apply involution action of |tw| on weight lattice to |v|
  void twisted_act
    (const TwistedInvolution& tw,Weight& v) const;

  TorusPart sample_torus_part(CartanNbr cn, RealFormNbr rf) const
  { return TorusPart(Cartan[cn].rep[rf],semisimpleRank()); }
  TorusPart dual_sample_torus_part(CartanNbr cn, RealFormNbr drf)
    const
  { return TorusPart(Cartan[cn].dual_rep[drf],semisimpleRank()); }


// manipulators

  Cartan_orbits& involution_table () { return C_orb; }

/* The main manipulator is |cartan|, which ensures the |CartanClass| is
   generated; most other manipulators are so because they call |cartan|
 */

//!\brief Returns data for stable conjugacy class \#cn of Cartan subgroups.
  const CartanClass& cartan(CartanNbr cn)
  { if (Cartan[cn].class_pt==NULL)
      add_Cartan(cn);
    return *Cartan[cn].class_pt;
  }

/*!\brief
  An element of the orbit in the adjoint fiber corresponding to |rf|
  in the classification of weak real forms for cartan |\#cn|.
  This amounts to searching for |rf| in |Cartan[cn].real_labels|.
*/
  unsigned long representative(RealFormNbr rf, CartanNbr cn)
    { return cartan(cn).fiber().weakReal().classRep(real_form_part(rf,cn)); }

/*!\brief
  An element of the orbit in the adjoint dual fiber corresponding to |drf|
  in the classification of dual weak real forms for cartan |\#cn|.
  This amounts to searching for |drf| in |Cartan[cn].dual_real_labels|.
*/
  unsigned long dualRepresentative(RealFormNbr drf, CartanNbr cn)
    { return cartan(cn).dualFiber().
	weakReal().classRep(dual_real_form_part(drf,cn));
    }

/*!\brief
  The size of the fiber orbits corresponding to strong real forms lying
  over weak real form \#rf, in cartan \#cn (all orbits have the same size)

  Explanation: this is the size of the orbits, for the shifted action of
  the imaginary Weyl group, that correspond to rf in the classification of
  strong real forms (they all have the same size.)

  Precondition: Real form \#rf is defined for cartan \#cn
*/
  unsigned long fiberSize(RealFormNbr rf, CartanNbr cn);

  unsigned long dualFiberSize(RealFormNbr, CartanNbr cn);

/*!
  \brief returns the number of involutions for the inner class.
*/
  InvolutionNbr numInvolutions();

/*!
  \brief returns the number of involutions for the indicated Cartans.
*/
  InvolutionNbr numInvolutions(const BitMap& Cartan_classes);

  RootNbrSet noncompactPosRootSet(RealFormNbr, size_t);

/*!
  \brief returns the number of elements in K\\G/B for real form \#rf.

  Explanation: this is exactly the number of elements in the one-sided
  parameter set corresponding to any strong real form of G lying over rf.
*/
  unsigned long KGB_size(RealFormNbr rf,
			 const BitMap& Cartan_classes);

  unsigned long KGB_size(RealFormNbr rf)
    { return KGB_size(rf,Cartan_set(rf)); }

/*!
  \brief returns the cardinality of the union of sets \f$K\backslash G/B\f$
  for this inner class.
*/
  unsigned long global_KGB_size();
/*!
  \brief the size of the block defined by the weak real form rf
  and the weak dual real form drf.
*/
  unsigned long block_size(RealFormNbr, RealFormNbr,
			   const BitMap& Cartan_classes);

  unsigned long block_size(RealFormNbr rf, RealFormNbr drf)
    { return block_size(rf,drf,Cartan_set(rf)& dual_Cartan_set(drf)); }

/*!
  \brief returns the real form labels for cartan \#cn

  More precisely, realFormLabels(cn)[i] is the (inner) number of the real form
  that corresponds to part i of the partition cartan(n).fiber().weakReal()
*/
  const RealFormNbrList& realFormLabels(CartanNbr cn)
  {
    cartan(cn); // make sure that labels for Cartan |cn| are generated
    return Cartan[cn].real_labels;
  }

/*!
  \brief returns the dual real form labels for cartan \#cn
*/
  const RealFormNbrList& dualRealFormLabels(CartanNbr cn)
  {
    cartan(cn);// make sure that dual labels for Cartan |cn| are generated
    return Cartan[cn].dual_real_labels;
  }

/*! get part in the |weakReal| partition of the fiber in Cartan \#cn
    corresponding to real form |rf|
 */
  cartanclass::adjoint_fiber_orbit real_form_part(RealFormNbr rf, CartanNbr cn)
  { return permutations::find_index(realFormLabels(cn),rf); }

/*! get part in the |weakReal| partition of the dual fiber in Cartan \#cn
    corresponding to dual real form |drf|
 */
  cartanclass::adjoint_fiber_orbit
    dual_real_form_part(RealFormNbr drf, CartanNbr cn)
  { return permutations::find_index(dualRealFormLabels(cn),drf); }


 private:
// auxiliary accessors
  void construct(); // does essential work, common to two constructors

TwistedInvolution
  reflection(RootNbr rn,const TwistedInvolution& tw) const;

//!\brief Tells whether Cartan \#cn is defined in real form \#rf.
bool is_defined(RealFormNbr rf, CartanNbr cn) const
  { return Cartan[cn].real_forms.isMember(rf); }

// auxiliary manipulators

void add_Cartan(CartanNbr cn);

void map_real_forms(CartanNbr cn);      // set |Cartan[cn].real_labels|
void map_dual_real_forms(CartanNbr cn); // set |Cartan[cn].dual_real_labels|

void correlateForms(CartanNbr cn);      // the old way to do the same things
void correlateDualForms(CartanNbr cn);  // now obsolete

}; // |class ComplexReductiveGroup|

} // |namespace complexredgp|

} // |namespace atlas|

#endif
