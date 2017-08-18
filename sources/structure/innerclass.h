/*
  This is innerclass.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006--2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/*
 Class definitions and function declarations for the class
 InnerClass.
*/

#ifndef INNERCLASS_H  /* guard against multiple inclusions */
#define INNERCLASS_H

#include "../Atlas.h"

#include "tags.h"
#include "bitmap.h"	// containment of bitmaps for real forms
#include "permutations.h"// containment of root twist

#include "cartanclass.h"// containment of |Fiber|
#include "involutions.h"// containment of |InvolutionTable|, |Cartan_orbits|
#include "poset.h"	// containment of Cartan poset
#include "rootdata.h"	// containment of root datum and its dual
#include "tits.h"	// containment of Tits group and its dual

namespace atlas {

namespace innerclass {

/******** function declarations **********************************************/


/* the next function should NOT be made into a method, suppressing the |rs|
   and |W| parameters, as these can be be dual to those in the |InnerClass|!
*/
void Cayley_and_cross_part(RootNbrSet& Cayley,
			   WeylWord& cross,
			   const TwistedInvolution& tw,
			   const RootSystem& rs,
			   const TwistedWeylGroup& W);

RealFormNbr real_form_of // who claims this KGB element?
  (InnerClass& G, TwistedInvolution tw, // by value
   const RatCoweight& torus_factor,
   RatCoweight& cocharacter // additional output
   );

RatCoweight some_coch // some cocharacter whose real form is in class |csc|
  (const InnerClass& G,cartanclass::square_class csc);

Grading grading_of_simples
  (const InnerClass& G, const RatCoweight& coch);

// find compact ones among imaginary simple roots for |G|, as defined by |coch|
 Grading compacts_for(const InnerClass& G, TorusElement coch);

containers::sl_list<TorusPart> preimage
  (const Fiber& fund_f, const cartanclass::square_class csc,
   const cartanclass::FiberElt y, const cartanclass::AdjointFiberElt image);


  // apply involution action of |tw| on weight lattice to |v|
void twisted_act
  (const InnerClass& G, const TwistedInvolution& tw,Weight& v);
void twisted_act
  (const InnerClass& G, const TwistedInvolution& tw,RatWeight& v);

void twisted_act
  (const InnerClass& G,Weight& v, const TwistedInvolution& tw);
void twisted_act
  (const InnerClass& G,RatWeight& v, const TwistedInvolution& tw);



/******** type definitions ***************************************************/


/*   Complex reductive group endowed with an inner class of real forms.

  This class computes those aspects of the structure theory of (an inner class
  of) real reductive groups G(R) that will be needed to describe the Langlands
  classification of irreducible representations of G(R). Since we look at an
  inner class of real forms, the first problem is to enumerate the different
  real forms constituting this inner class.

  We list in d_cartanSet the conjugacy classes of real Cartan subgroups up to
  stable conjugacy; this classification does not refer to a particular real
  form. However, the enumeration of the real forms actually takes place during
  the construction of the first (fundamental) Cartan subgroup. Each stable
  class corresponds to at most one conjugacy class of real Cartan subgroups in
  each real form; so for each stable class of Cartan subgroups, we enumerate
  the real forms over which it is defined (the fundamental Cartan subgroup is
  defined for all real forms).

  We compute the structure of the real Cartan subgroups (notably the groups of
  connected components); this depends only on the stable conjugacy class. We
  determine the real Weyl groups of Cartan subgroups (which are _almost_
  constant across the stable class, but not quite).

  Everything is determined by (and computed from) two things: the based root
  datum recorded in the RootDatum class d_rootDatum, and its involutive
  distinguished automorphism. Many computations take place inside the Tits
  group, which is an extension of the (complex) Weyl group by the elements of
  order 2 in the torus. (In fact the structure we store in |d_titsGroup|
  allows computing in an even larger group, the semidirect product of the Tits
  group just described by a factor Z/2Z whose action on the other factor is
  determined by the given automorphism of the based root datum, the "twist".)

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
  $t\Gamma(t)=e$ and "twisted conjugacy" of $t$ by $w\in W$ is given by
  $w\cdot t=wt\Gamma(w^{-1})$. The stable conjugacy classes of Cartan
  subgroups will each be represented by a canonical representative of the
  corresponding twisted conjugacy class of twisted involutions.

  In addition to describing the set of Cartan classes, this class provides
  access (via the |Cartan| array) to data for each individual one of them, and
  (via |Cartan_poset|) to the partial order relation between them. The
  structure |C_info| contains the imporatant |CartanClass| field |Cc|, which
  holds most of the information about the Cartan class (in two |Fiber|
  structures, an ordinary and dual fiber, see the \.{cartanclass} module);
  this information is generated for all Cartan classes upon construction of
  the |InnerClass|. The remaining fields of |C_info| reflect an
  old and now reverted design decision, in which only a pointer to a
  |CartanClass| was held, which was intially null and "filled" on demand. This
  meant that real forms, which are identified as orbits in the adjoint fiber
  group, were generated without access to the |CartanClass| information. In
  fact, although that information is now available right away, the real forms
  are still genrated without using that information. Instead an "adjoint"
  instance of the Tits group (of which another instance will serve similalrly
  for generation of KGB elements) is used to find elements representing the
  real form at each of its Cartan classes. The other fields of |C_info| store
  minimal information derived from this generation: a single torus part |rep|
  of an adjoint Tits group element for each real form over which the Cartan
  class is defined, a dual counterpart |dual_rep|, and some complementary
  information to help identify real forms across different Cartans.

  For the partial order relation, let |tau_i| be involutions acting on the
  complex torus |H| for various classes of Cartan subgroups; (H,tau_1) is
  considered "more compact" than (H,tau_2) if the identity component of the
  fixed point set H^tau_2 is W-conjugate to a subtorus of H^tau_1.

  The problem for the dual group of G is identical, the bijection taking the
  negative transpose of a twisted involution. This bijection reverses the
  partial order on Cartans. The class also provides access to Cartans in the
  dual group.

  */
class InnerClass
{
  // The based root datum. It is stored here (constructed by our constructor)
  const RootDatum d_rootDatum;

  // The dual based root datum. It is also stored and constructed by us
  const RootDatum d_dualRootDatum;

  const WeylGroup* my_W; // pointer to |W| in case we own |W|, or |NULL|
  const WeylGroup& W;    // possibly owned (via |my_W|) reference

  /*
    Fiber class for the fundamental Cartan subgroup
    The distinguished involution (permuting the simple roots) is stored here
    so for construction it is convient to precede the |d_titsGroup| member
  */
  Fiber d_fundamental;

  /*
    Fiber class for the fundamental Cartan in the dual group.
    The fiber group here is the group of characters (i.e., the dual group)
    of the component group of the quasisplit Cartan.
  */
  Fiber d_dualFundamental;

  // Tits group of the based root datum, extended by an involutive automorphism
  // this member also stores the |TwistedWeylGroup| (which refers to |W|)
  const TitsGroup d_titsGroup;
  // Tits group of the dual based root datum
  const TitsGroup d_dualTitsGroup;
  // the permutation of the roots given by the based automorphism
  const Permutation root_twist;


  struct C_info
  { // gradings of the set of simple (co)roots, for all real forms
    typedef std::vector<RankFlags> form_reps;

    TwistedInvolution tw;
    BitMap real_forms,dual_real_forms; // mark present (dual) real forms
    form_reps rep,dual_rep; // gradings representing those (dual) real forms
    BitMap below; // numbers of Cartan classes below this in partial ordering
    CartanClass Cc; // detailed information, basically |Fiber| and dual |Fiber|
    RealFormNbrList real_labels,dual_real_labels;

    C_info(const InnerClass& G, const TwistedInvolution twi,
	   CartanNbr i); // index |i| of this Cartan used to dimension |below|

  }; // |C_info|
  std::vector<C_info> Cartan;

  // Partial order of Cartan subgroups.
  poset::Poset Cartan_poset;

  // list of the most split Cartan classes for each real form
  std::vector<CartanNbr> d_mostSplit;

  // a general repository for involutions, organised by conjugacy class
  Cartan_orbits C_orb;

 public:
// constructors and destructors
  InnerClass(const PreRootDatum&, // constructor builds root datum
			const WeightInvolution&);

  InnerClass(const RootDatum&, // alternative that copies root datum
			const WeightInvolution&);

  InnerClass(const InnerClass&, tags::DualTag);

// copy, assignment and swap are not needed, and therefore forbidden
  InnerClass(const InnerClass&) = delete;
  InnerClass& operator= (const InnerClass&) = delete;
  void swap(InnerClass& G) = delete;


  ~InnerClass();

// Accessors

// Attributes "inherited" from component objects

  const RootDatum& rootDatum() const { return d_rootDatum; }
  const RootDatum& dualRootDatum() const { return d_dualRootDatum; }
  const RootSystem& rootSystem() const {return d_rootDatum; } // base object
  const RootSystem& dualRootSystem() const {return d_dualRootDatum; } // base
  size_t rank() const { return rootDatum().rank(); }
  size_t semisimpleRank() const { return rootSystem().rank(); }

// Access to certain component objects themselves

  const WeylGroup& weylGroup() const { return W; }
  const TwistedWeylGroup& twistedWeylGroup() const
    { return d_titsGroup; } // in fact its base object
  const TwistedWeylGroup& dualTwistedWeylGroup() const
    { return d_dualTitsGroup; } // in fact its base object
  const TitsGroup& titsGroup() const { return d_titsGroup; }
  const TitsGroup& dualTitsGroup() const { return d_dualTitsGroup; }

  // a general repository for involutions, organised by conjugacy class
  const Cartan_orbits& involution_table () const { return C_orb; }


// Attributes of the inner class as a whole

  Permutation simple_twist() const
    { return Permutation
	(&twistedWeylGroup().twist()[0],
	 &twistedWeylGroup().twist()[semisimpleRank()]); }

  RankFlags simple_roots_imaginary() const; // $\xi$-fixed simple roots
  RankFlags simple_roots_real() const; // $(-w_0\xi^t)$-fixed simple roots

  const Permutation& root_involution() const { return root_twist; }
  RootNbr twisted_root(RootNbr alpha) const { return root_twist[alpha]; }

  /* partial ordering of Cartan classes
     This is the ordering by containment of H^theta up to conjugacy:
     (H,theta_1) precedes (H,theta_2) if (H^theta_2)_0 is W-conjugate to a
     subtorus of H^theta_1. Numbering of elements is as in the vector |Cartan|
  */
  const poset::Poset& Cartan_ordering() const { return Cartan_poset; }

  // the distinguished involution for G and for its dual group
  const WeightInvolution& distinguished() const
    { return d_fundamental.involution(); }
  const CoweightInvolution& dualDistinguished() const
    { return d_dualFundamental.involution(); }

  // number of conjugacy classes of Cartan subgroups.
  CartanNbr numCartanClasses() const { return Cartan.size(); }
  // number of weak real forms of G.
  RealFormNbr numRealForms() const { return d_fundamental.numRealForms(); }
  // number of weak real forms of the dual group of G.
  RealFormNbr numDualRealForms() const
  { return d_dualFundamental.numRealForms(); }
  // total number of involutions for the inner class.
  InvolutionNbr numInvolutions() const;

  // size of a union of sets $K\backslash G/B$ for relevant strong real forms
  unsigned long global_KGB_size() const;

  // the (inner) number of the quasisplit real form.
  RealFormNbr quasisplit() const { return RealFormNbr(0); }


// Information about individual Cartan classes

  // number of weak real forms of G over which this Cartan is defined.
  RealFormNbr numRealForms(CartanNbr cn) const
  { return Cartan[cn].real_forms.size(); }
  // number of dual real forms over which the dual of this Cartan is defined
  RealFormNbr numDualRealForms(CartanNbr cn) const
  { return Cartan[cn].dual_real_forms.size(); }

  // the precise sets over which (the dual of) this Cartan is defined
  const BitMap& real_forms(CartanNbr cn) const
  { return Cartan[cn].real_forms; }
  const BitMap& dual_real_forms(CartanNbr cn) const
  { return Cartan[cn].dual_real_forms; }

  // a twisted involution representative of this Cartan class
  const TwistedInvolution& involution_of_Cartan(CartanNbr cn) const
    { return Cartan[cn].tw; }

  // remaining information is stored in |CartanClass| object
  const CartanClass& cartan(CartanNbr cn) const { return Cartan[cn].Cc; }


// Information selected by subset of the Cartan classes

  // number of involutions for the indicated Cartans
  InvolutionNbr numInvolutions(const BitMap& Cartan_classes) const;


// Information about individual (weak) real forms or dual real forms

  // the subset of Cartan classes defined over the real form |rf|
  BitMap Cartan_set(RealFormNbr rf) const;

  // the set of Cartan classes for the dual real form |drf|
  BitMap dual_Cartan_set(RealFormNbr drf) const;

  // the most split Cartan for this real form
  CartanNbr mostSplit(RealFormNbr rf) const { return d_mostSplit[rf]; }

  /* a set of noncompact imaginary roots at the distinguished involution
     for a representative of this real form */
  RootNbrSet noncompactRoots(RealFormNbr rf) const
  { return d_fundamental.noncompactRoots(d_fundamental.wrf_rep(rf)); }

  /* a set of noncompact imaginary roots at the distinguished involution
     for a representative of this real form */
  RootNbrSet parity_coroots(RealFormNbr drf) const
  { return d_dualFundamental.noncompactRoots(d_dualFundamental.wrf_rep(drf)); }

  // the number of elements in K\\G/B for real form |rf|.
  unsigned long KGB_size(RealFormNbr rf) const
  { return KGB_size(rf,Cartan_set(rf)); }
  // the same limited to indicated Cartan classes only
  unsigned long KGB_size(RealFormNbr rf, const BitMap& Cartan_classes) const;

  // compact simple roots ("grading shift") at first element of KGB for |rf|
  RankFlags simple_roots_x0_compact(RealFormNbr rf) const;

  cartanclass::square_class xi_square(RealFormNbr rf) const;
  RealFormNbr square_class_repr(cartanclass::square_class csc) const;

  TorusPart x0_torus_part(RealFormNbr rf) const;

  // torus parts that remain in the fiber and do not affect any grading
  containers::sl_list<TorusPart> central_fiber(RealFormNbr rf) const;

// Information about a real form or dual real form at a given Cartan class

  // size of the block defined by weak real form |rf| and dual real form| drf|
  arithmetic::big_int block_size(RealFormNbr rf, RealFormNbr drf) const
    { return block_size(rf,drf,Cartan_set(rf)& dual_Cartan_set(drf)); }
  // the same limited to indicated Cartan classes only
  arithmetic::big_int block_size(RealFormNbr rf, RealFormNbr drf,
				 const BitMap& Cartan_classes) const;

  // more functions that interrogate the fundametal fiber
  cartanclass::StrongRealFormRep sample_strong_form (RealFormNbr rf) const
  { return d_fundamental.strongRealForm(rf); }
  unsigned long fundamental_fiber_size() const
  { return d_fundamental.fiberSize(); }
  const Partition&
    fundamental_fiber_partition(cartanclass::square_class csc) const
  { return d_fundamental.fiber_partition(csc); }
  TorusPart lift_from_fundamental_fiber(unsigned long x) const
  { const Fiber& fund=d_fundamental;
    SmallBitVector v (static_cast<RankFlags>(x),fund.fiberRank());
    return fund.fiberGroup().fromBasis(v);
  }
  cartanclass::FiberElt to_fundamental_fiber(TorusPart t) const
  { return d_fundamental.fiberGroup().toBasis(t); }

  // list torus parts mapping to strong form of |y| and further to |image|
  containers::sl_list<TorusPart> torus_parts_for_grading_shift
    (const cartanclass::square_class csc,
     const cartanclass::FiberElt y, const cartanclass::AdjointFiberElt image)
  const;

  const Partition& weak_real_partition() const
  { return d_fundamental.weakReal(); }
  const Partition& dual_weak_real_partition() const
  { return d_dualFundamental.weakReal(); }

  // size of fibers in KGB set for |rf| over any involution in Cartan class |cn|
  unsigned long fiberSize(RealFormNbr rf, CartanNbr cn) const;
  unsigned long dualFiberSize(RealFormNbr drf, CartanNbr cn) const;

  // real forms, vector indexed by parts in |cartan(n).fiber().weakReal()|
  const RealFormNbrList& realFormLabels(CartanNbr cn) const
  { return Cartan[cn].real_labels; }
  // dual real forms, indexed by parts in |cartan(n).dual_fiber().weakReal()|
  const RealFormNbrList& dualRealFormLabels(CartanNbr cn) const
  { return Cartan[cn].dual_real_labels; }

  // part in the |weakReal| partition at Cartan |cn| corresponding to |rf|
  cartanclass::adjoint_fiber_orbit
    real_form_part(RealFormNbr rf, CartanNbr cn) const
  { return permutations::find_index(realFormLabels(cn),rf); }
  cartanclass::adjoint_fiber_orbit
    dual_real_form_part(RealFormNbr drf, CartanNbr cn) const
  { return permutations::find_index(dualRealFormLabels(cn),drf); }

  // adjoint fiber element in fiber of Cartan |cn|, representative of |rf|
  cartanclass::AdjointFiberElt
    representative(RealFormNbr rf, CartanNbr cn) const
  { return cartan(cn).fiber().wrf_rep(real_form_part(rf,cn)); }
  // adjoint fiber element in dual fiber of Cartan |cn|, representative of |drf|
  cartanclass::AdjointFiberElt
    dualRepresentative(RealFormNbr drf, CartanNbr cn) const
  { return cartan(cn).dualFiber().wrf_rep(dual_real_form_part(drf,cn)); }


// Information about individual twisted involutions

  CartanNbr class_number(TwistedInvolution) const; // by value
  Weight posRealRootSum(const TwistedInvolution&) const; // 2rho_{re}
  Weight posImaginaryRootSum(const TwistedInvolution&) const; // 2rho_{im}

  // matrix giving involution action of |tw| on weight lattice (lookup)
  const WeightInvolution& matrix(const TwistedInvolution& tw) const
  { return C_orb.matrix(tw); }

  InvolutionData involution_data(const TwistedInvolution& tw) const
  { return InvolutionData::build(rootSystem(),twistedWeylGroup(),tw); }

/* Make |sigma| canonical and return Weyl group |w| element that twisted
   conjugates the canonical representative back to |sigma|. Thus after the
   call the involution corresponding to |sigma| satisfies:

  (1) the sum of the positive real roots is dominant (call it $SR$)
  (2) in the subsystem of roots orthogonal to $SR$ (which contains all the
      imaginary roots), the sum of the imaginary roots (call it $SI$) is
      dominant \emph{for the subsystem}
  (3) in the subsystem of roots orthogonal to both $SR$ and $SI$, the
      involution corresponding to twisted involution fixes (globally) the
      dominant chamber of the subsystem (it permutes its simple roots).
*/
  WeylWord canonicalize(TwistedInvolution& sigma, RankFlags gens) const;
  WeylWord canonicalize(TwistedInvolution& sigma) const
  { return canonicalize(sigma,RankFlags(constants::lMask[semisimpleRank()])); }


// Manipulators

  void generate_Cartan_orbit (CartanNbr i) { C_orb.add(*this,i); }

// Auxiliary accessors
 private:

  void construct(); // does essential work, common to two constructors

  TorusPart grading_shift_repr(Grading diff) const;

  // whether Cartan \#cn is defined over real form \#rf.
  bool is_defined(RealFormNbr rf, CartanNbr cn) const
  { return Cartan[cn].real_forms.isMember(rf); }

  // adjoint torus parts, in fundamental coweight basis, \emph{are} gradings
  TorusPart sample_torus_part(CartanNbr cn, RealFormNbr rf) const
  { return TorusPart(Cartan[cn].rep[rf],semisimpleRank()); }
  TorusPart dual_sample_torus_part(CartanNbr cn, RealFormNbr drf) const
  { return TorusPart(Cartan[cn].dual_rep[drf],semisimpleRank()); }

// Auxiliary manipulators

  void map_real_forms(CartanNbr cn);      // set |Cartan[cn].real_labels|
  void map_dual_real_forms(CartanNbr cn); // set |Cartan[cn].dual_real_labels|

}; // |class InnerClass|

} // |namespace innerclass|

} // |namespace atlas|

#endif
