/*
  This is complexredgp.cpp.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/


/*
  Implementation for the class ComplexReductiveGroup.

  The ComplexReductiveGroup class will play a central role in the whole
  program. Even though it is entirely defined by its based root datum and an
  involutive automorphism of that datum, it has seemed more natural to use
  this class to collect the wealth of combinatorial data that they give rise
  to, and that will serve as the basis for our description of the
  representation theory of a real group G. Note that the current state of the
  theory, and most notably Vogan duality, makes it natural and necessary to
  consider all the real forms of our complex group (in a given inner class) at
  once; so that is another reason to not choose a real form a priori.
*/

#include "complexredgp.h"
#include "realredgp.h" // for |square_class_choice|

#include <set>

#include "lattice.h"
#include "weyl.h"
#include "kgb.h"


/*****************************************************************************

  This module, together with the cartanclass one, contains code for dealing
  with conjugacy classes of Cartan subgroups, and real forms.

  The enumeration of weak real forms amounts to that of $W^\delta$-orbits in
  the fundamental fiber of the one-sided parameter space for the adjoint group
  (see the "combinatorics" paper on the Atlas website, or the "algorithms"
  paper by Jeff Adams and Fokko du Cloux); this is a very small computation,
  depending only on the Lie algebra (|RootSystem| only).

  The enumeration of conjugacy classes of Cartan subgroups, for the various
  real forms, is part of the enumeration of conjugacy classes of root data
  involutions for the given inner class, or equivalently that of twisted
  involutions in the Weyl group, which is a pure Weyl group computation.
  Fokko's original code proceeded more or less by listing each of these
  conjugacy classes, checking each new twisted involution for membership of
  each class of previously generated ones. Now we find for each conjugacy
  class of twisted involutions a canonical representative, which is relatively
  easy to compute. In this way we avoid enumeration of each conjugacy class.

  The identification of real forms can be done using the list of Cartan
  classes as well, via the unique most split Cartan class for each real form.

  The most delicate part is the "correlation" part: for each Cartan, tell
  which orbit in the corresponding fiber corresponds to which real form, the
  real forms being labelled by the orbits in the fundamental fiber. In the
  current version, the solution to this problem is cleaner then previously,
  and perfectly general: it is obtained by writing out a system of equations
  for the grading defining the real form; this system does not always have a
  unique solution, but all solutions correspond to the same real form.


******************************************************************************/

namespace atlas {

namespace complexredgp {

  void crossTransform(RootNbrList&,
		      const WeylWord&,
		      const RootSystem&);

  unsigned long makeRepresentative(const Grading&,
				   const RootNbrList&,
				   const Fiber&);

  bool checkDecomposition(const TwistedInvolution& ti,
			  const WeylWord& cross,
			  const RootNbrSet& Cayley,
			  const TwistedWeylGroup& W,
			  const RootSystem& rs);

  // make from $(\Q/2\Z)^n$ to $(\Z/2\Z)^n$ taking floors of coefficients
  SmallBitVector floor (Ratvec_Numer_t num, arithmetic::Numer_t d);

  // select |TorusPart| from list with minimal fingerprint for |coch|
  TorusPart minimum(const containers::sl_list<TorusPart>& cf,
		    const ComplexReductiveGroup& G, const TorusElement& coch);

  // replace |coch| by minimum representative modulo image |delta_plus_1|
  void to_minimum_representative(TorusElement& coch,
				 const WeightInvolution& delta);


/*****************************************************************************

        Chapter I -- The ComplexReductiveGroup class

******************************************************************************/

ComplexReductiveGroup::C_info::C_info
  (const ComplexReductiveGroup& G,const TwistedInvolution twi, CartanNbr i)
  : tw(twi)
  , real_forms(G.numRealForms()), dual_real_forms(G.numDualRealForms())
  , rep(G.numRealForms()),        dual_rep(G.numDualRealForms())
  , below(i)
  , Cc(CartanClass(G.rootDatum(),G.dualRootDatum(),
		   G.compute_matrix(tw))) // generate fiber and dual fiber
  , real_labels(), dual_real_labels() // these start out emtpy
  {}


/*
  Main constructor

  Constructs a |ComplexReductiveGroup| from a pre-rootdatum |rd| and a
  distinguished involution |d|, which stabilises the set of simple roots
*/
ComplexReductiveGroup::ComplexReductiveGroup
 (const PreRootDatum& prd, const WeightInvolution& tmp_d)
  : d_rootDatum(prd)
  , d_dualRootDatum(d_rootDatum,tags::DualTag())

  , my_W(new WeylGroup(d_rootDatum.cartanMatrix()))
  , W(*my_W) // owned when this constructor is used

  , d_fundamental(d_rootDatum,tmp_d) // will also be fiber of cartan(0)
  , d_dualFundamental(d_dualRootDatum,dualBasedInvolution(tmp_d,d_rootDatum))
    // dual fundamental fiber is dual fiber of most split Cartan

  , d_titsGroup(d_rootDatum,W,distinguished())
  , d_dualTitsGroup(d_dualRootDatum,W,dualDistinguished())
  , root_twist(d_rootDatum.root_permutation(simple_twist()))

  , Cartan(1,C_info(*this,TwistedInvolution(),0))
  , Cartan_poset() // poset is extended and populated below
  , d_mostSplit(numRealForms(),0) // values 0 may be increased below

  // DON'T use |tmp_d| as involution below: would store a dangling reference!
  , C_orb(d_rootDatum,distinguished(),d_titsGroup) // set up bare table
{
  construct(); // set up |Cartan| data, but only resizes |C_orb|
}

/*
  Variant constructor, differs only by using a constructed root datum

  Constructs a |ComplexReductiveGroup| from a rootdatum |rd| and a
  distinguished involution |d|, which stabilises the set of simple roots
*/
ComplexReductiveGroup::ComplexReductiveGroup
 (const RootDatum& rd, const WeightInvolution& tmp_d)
  : d_rootDatum(rd)
  , d_dualRootDatum(d_rootDatum,tags::DualTag())

  , my_W(new WeylGroup(d_rootDatum.cartanMatrix()))
  , W(*my_W) // owned when this constructor is used

  , d_fundamental(d_rootDatum,tmp_d) // will also be fiber of cartan(0)
  , d_dualFundamental(d_dualRootDatum,dualBasedInvolution(tmp_d,d_rootDatum))
    // dual fundamental fiber is dual fiber of most split Cartan

  , d_titsGroup(d_rootDatum,W,distinguished())
  , d_dualTitsGroup(d_dualRootDatum,W,dualDistinguished())
  , root_twist(d_rootDatum.root_permutation(simple_twist()))

  , Cartan(1,C_info(*this,TwistedInvolution(),0))
  , Cartan_poset() // poset is extended and populated below
  , d_mostSplit(numRealForms(),0) // values 0 may be increased below

  // DON'T use |tmp_d| as involution below: would store a dangling reference!
  , C_orb(d_rootDatum,distinguished(),d_titsGroup) // set up bare table
{
  construct(); // set up |Cartan| data, but only resizes |C_orb|
}

void ComplexReductiveGroup::construct() // common part of two constructors
{
  { // task 1: generate Cartan classes, fill non-dual part of |Cartan|
    { // complete initialisation of |Cartan[0]|
      const Fiber& f=fundamental();
      const Partition& weak_real=f.weakReal();
      // fill initial |form_reps| vector with assignment from |weak_real|
      for (RealFormNbr i=0; i<weak_real.classCount(); ++i)
      {
	Cartan[0].real_forms.insert(i); // Cartan 0 exists at all real forms
	// setting the initial torus part for each real form is subtle.
	Cartan[0].rep[i]=              // in |adj_Tg|, a torus part is a
	  cartanclass::restrictGrading // |Grading| of simple roots, where
	  (f.compactRoots(weak_real.classRep(i)), // compact ones need bit set
	   d_rootDatum.simpleRootList()); // to flip base (noncompact) grading
      }
    }

    const TitsCoset adj_Tg(*this);     // based adjoint Tits group
    const TitsGroup& Tg=adj_Tg.titsGroup(); // same, forgetting "based" stuff

    for (CartanNbr i=0; i<Cartan.size(); ++i) // |Cartan| grows as loop advances
    {
      Cartan_poset.new_max(Cartan[i].below); // include now completed level

#ifndef NDEBUG
      CartanNbr entry_level=Cartan.size(); // is only used in |assert|
#endif
      InvolutionData id =
	InvolutionData::build(d_rootDatum,d_titsGroup,Cartan[i].tw);
      RootNbrSet pos_im = id.imaginary_roots() & d_rootDatum.posRootSet();

      // try to generate new Cartans above current: do Cayleys by |pos_im|
      for (RootNbrSet::iterator it=pos_im.begin(); it(); ++it)
      {
	RootNbr alpha=*it;
	SmallBitVector alpha_bin(d_rootDatum.inSimpleRoots(alpha));

	TitsElt a(Tg,Cartan[i].tw);// test element with null torus part
	WeylWord conjugator; // will right-conjugate from |Cartan[i].tw|

	weyl::Generator s; // outside loop to allow inspection of final value

	// loop while |alpha| not simple, i.e., has a descent |s| not itself
	while (alpha!=
	       d_rootDatum.simpleRootNbr(s=d_rootDatum.find_descent(alpha)))
	{
	  conjugator.push_back(s);
	  adj_Tg.basedTwistedConjugate(a,s);
	  d_rootDatum.simple_reflect_root(s,alpha);
	}

	// fix how test element at |Cartan[i].tw| grades original root $\alpha$
	bool zero_grading = adj_Tg.simple_grading(a,s); // whether noncompact

	TwistedInvolution sigma = // involution after "Cayley transform"
	  W.prod(s,a.tw());       // starting from |a.tw()|
	WeylWord ww = // Weyl word that will conjugate canonical back to current
	  canonicalize(sigma); // and sigma is now canonical elt in new Cartan

	CartanNbr ii;
	for (ii=0; ii<Cartan.size(); ++ii)
	  if (Cartan[ii].tw==sigma) // see if |sigma| matches known canonical
	    break; // found a previously encountered Cartan class

	if (ii==Cartan.size()) // if not seen before, create a new Cartan
	  Cartan.push_back(C_info(*this,sigma,ii)); // with involution |sigma|

	Cartan[ii].below.insert(i); // mark parent Cartan as below |sigma|

	// see which real forms (newly) carry over to the new Cartan class
	for (BitMap::iterator rfi=Cartan[i].real_forms.begin(); rfi(); ++rfi)
	{
	  RealFormNbr rf = *rfi;
	  const RankFlags in_rep = Cartan[i].rep[rf];
	  RankFlags& out_rep = Cartan[ii].rep[rf]; // to be filled in
	  TorusPart tp(in_rep,alpha_bin.size()); // size semisimple rank
	  if (alpha_bin.dot(tp)!=zero_grading // |tp| makes $\alpha$ noncompact
	      and not Cartan[ii].real_forms.isMember(rf)) // and |rf| is new
	  { // this is the first hit of |ii| for |rf|
	    assert(ii>=entry_level); // we may populate only newborn sets
	    Cartan[ii].real_forms.insert(rf); // mark that |rf| has new Cartan

	    TitsElt x(Tg,tp,Cartan[i].tw); // make adjoint Tits group element
	    adj_Tg.basedTwistedConjugate(x,conjugator); // at old Cartan, and
	    adj_Tg.Cayley_transform(x,s);               // move to new Cartan
	    adj_Tg.basedTwistedConjugate(x,ww); // right-conjugate: to $\sigma$
	    assert(x.tw()==sigma); // check we arrived at intended involution

	    out_rep=Tg.left_torus_part(x).data(); // unpack, store torus part
	    d_mostSplit[rf]=ii; // the last such assignment for |rf| sticks
	  }
	} // |for (rf)|
      } // |for (alpha)|; new Cartans above |Cartan[i]| are now found
    } // |for (i<Cartan.size())|

    C_orb.set_size(Cartan.size()); // dimension |C_orb|, but leave it at that
  } // task 1 (Cartan class generation)

  { // task 2: fill remainder of all |Cartan[i]|: |dual_real_forms|, |dual_rep|
    const TitsCoset dual_adj_Tg (*this,tags::DualTag());
    const TitsGroup& dual_Tg = dual_adj_Tg.titsGroup();

    { // first initialise |Cartan.back().dual_rep|
      TwistedInvolution w0(W.longest()); // this is always a twisted inv.
      WeylWord ww = canonicalize(w0); // but not always canonical
      assert(w0==Cartan.back().tw);

      const Fiber& f=dualFundamental();
      const Partition& weak_real=f.weakReal();
      // fill initial |form_reps| vector with assignment from |weak_real|
      for (unsigned long i=0; i<weak_real.classCount(); ++i)
      {
	// most split Cartan exists at all dual real forms
	Cartan.back().dual_real_forms.insert(i);

	RankFlags gr = // as above, torus part is obtained as a grading
	  cartanclass::restrictGrading // values at simple roots give torus part
	  (f.compactRoots(weak_real.classRep(i)), // compact ones need a set bit
	   d_rootDatum.simpleRootList()); // reproduces grading at imag. simples
	TitsElt x(dual_Tg,TorusPart(gr,d_rootDatum.semisimpleRank()));
	dual_adj_Tg.basedTwistedConjugate(x,ww); // transform to canonical
	Cartan.back().dual_rep[i] = dual_Tg.left_torus_part(x).data();
      }
    }

    for (CartanNbr i=Cartan.size(); i-->0; )
    {
      InvolutionData id =
	InvolutionData::build(d_rootDatum,d_titsGroup,Cartan[i].tw);
      RootNbrSet pos_re = id.real_roots() & d_rootDatum.posRootSet();
      for (RootNbrSet::iterator it=pos_re.begin(); it(); ++it)
      {
	RootNbr alpha=*it;
	SmallBitVector alpha_bin(d_rootDatum.inSimpleCoroots(alpha));

	TwistedInvolution tw = Cartan[i].tw; // non-dual
	TwistedInvolution tw_dual = W.opposite(tw);

	// create a test element with null torus part
	TitsElt a(dual_Tg,tw_dual);
	WeylWord conjugator;

	size_t j; // declare outside loop to allow inspection of final value
	while (alpha!=d_rootDatum.simpleRootNbr
				   (j=d_rootDatum.find_descent(alpha)))
	{
	  conjugator.push_back(j);
	  dual_adj_Tg.basedTwistedConjugate(a,j);
	  twistedWeylGroup().twistedConjugate(tw,j);
	  d_rootDatum.simple_reflect_root(j,alpha);
	}

	bool zero_grading = dual_adj_Tg.simple_grading(a,j);

	assert(tw==W.opposite(a.tw())); // coherence with dual group

	W.leftMult(tw,j); // "Cayley transform"
	WeylWord ww=canonicalize(tw); // in non-dual setting

	CartanNbr ii;

	for (ii=i; ii-->0; )
	  if (Cartan[ii].tw==tw)
	    break; // found a previously encountered Cartan class (must happen)
	assert(ii<Cartan.size() and Cartan[i].below.isMember(ii));

	for (BitMap::iterator
	       drfi=Cartan[i].dual_real_forms.begin(); drfi(); ++drfi)
	{
	  RealFormNbr drf = *drfi;
	  const RankFlags in_rep = Cartan[i].dual_rep[drf];
	  RankFlags& out_rep = Cartan[ii].dual_rep[drf];
	  TorusPart tp(in_rep,alpha_bin.size());
	  if (alpha_bin.dot(tp)!=zero_grading)
	  {
	    if (not Cartan[ii].dual_real_forms.isMember(drf))
	    { // this is the first hit of |ii| for |rf|
	      Cartan[ii].dual_real_forms.insert(drf);

	      TitsElt x(dual_Tg,tp,tw_dual);
	      dual_adj_Tg.basedTwistedConjugate(x,conjugator);
	      dual_adj_Tg.Cayley_transform(x,j);
	      dual_adj_Tg.basedTwistedConjugate(x,ww);
	      assert(x.tw()==W.opposite(tw));

	      out_rep=dual_Tg.left_torus_part(x).data();
	    }
	  }
	} // for (rf)
      } // for (alpha)
    } // |for (i=Cartan.size()-->0)|

  } // task 2

  { // task 3: set fields |real_labels|, |dual_real_labels| in each Cartan[i]
    for (CartanNbr cn=0; cn<Cartan.size(); ++cn)
    {
      map_real_forms(cn);      // used to be |correlateForms(cn);|
      map_dual_real_forms(cn); // used to be |correlateDualForms(cn);|
    }
  }
} // |ComplexReductiveGroup::construct|

// Construct the complex reductive group dual to G
ComplexReductiveGroup::ComplexReductiveGroup(const ComplexReductiveGroup& G,
					     tags::DualTag)
  : d_rootDatum(G.d_dualRootDatum)
  , d_dualRootDatum(G.d_rootDatum)

  , my_W(NULL), W(G.W) // not owned here, we depend on existence of |G|

  , d_fundamental(G.d_dualFundamental)
  , d_dualFundamental(G.d_fundamental)

  , d_titsGroup(G.d_dualTitsGroup,W)
  , d_dualTitsGroup(G.d_titsGroup,W)
  , root_twist(d_rootDatum.root_permutation(simple_twist()))

  , Cartan() // filled below
  , Cartan_poset(G.Cartan_poset,tags::DualTag())
  , d_mostSplit(numRealForms(),0) // values 0 may be increased below

  , C_orb(d_rootDatum,d_fundamental.involution(),d_titsGroup)
{
  Cartan.reserve(G.Cartan.size());
  C_orb.set_size(G.Cartan.size());

  for (CartanNbr i=G.Cartan.size(); i-->0; )
  {
    const C_info& src = G.Cartan[i];

    const TwistedInvolution tw_org = W.opposite(src.tw);
    TwistedInvolution canon_tw = tw_org;
    WeylWord conjugator = canonicalize(canon_tw);

    Cartan.push_back(C_info(*this,canon_tw,Cartan.size()));
    C_info& dst = Cartan.back();

    dst.real_forms = src.dual_real_forms;
    dst.dual_real_forms = src.real_forms;
    dst.rep = src.dual_rep; // these are torus parts at |tw_org|, and this
    dst.dual_rep = src.rep; // assignment is mainly to set their size in |dst|

    const TitsCoset adj_Tg(*this);     // based adjoint Tits group
    const TitsGroup& Tg=adj_Tg.titsGroup(); // same, forgetting base
    const TitsCoset dual_adj_Tg (*this,tags::DualTag());
    const TitsGroup& dual_Tg = dual_adj_Tg.titsGroup();

    for (BitMap::iterator it=dst.real_forms.begin(); it(); ++it)
    {
      TitsElt x(Tg,TorusPart(dst.rep[*it],semisimpleRank()),tw_org);
      adj_Tg.basedTwistedConjugate(x,conjugator);
      assert(x.tw()==dst.tw);
      dst.rep[*it] =Tg.left_torus_part(x).data();
      d_mostSplit[*it]=Cartan.size()-1; // last occurrence of |*it| will stick
    }
    for (BitMap::iterator it=dst.dual_real_forms.begin(); it(); ++it)
    {
      TitsElt y(dual_Tg,TorusPart(dst.dual_rep[*it],semisimpleRank()),src.tw);
      dual_adj_Tg.basedTwistedConjugate(y,conjugator);
      assert(y.tw()==W.opposite(dst.tw));
      dst.dual_rep[*it] = dual_Tg.left_torus_part(y).data();
    }

    for (CartanNbr j=i+1; j<G.Cartan.size(); ++j)
      if (G.Cartan[j].below.isMember(i))
	dst.below.insert(G.Cartan.size()-1-j);

    {
      CartanNbr cn=Cartan.size()-1; // ensure |CartanClass object is created|
      map_real_forms(cn);      // used to be |correlateForms(cn);|
      map_dual_real_forms(cn); // used to be |correlateDualForms(cn);|
    }
  }

}


// destruction frees Weyl group if owned
ComplexReductiveGroup::~ComplexReductiveGroup()
{
  delete my_W;
}

/********* copy, assignment and swap *****************************************/

// This should remain empty!

/*****************************************************************************

        Chapter II --- private (auxiliary) methods for ComplexReductiveGroup

******************************************************************************/

/******** private accessors **************************************************/

/* Compose |tw| to the left with the reflection $s_\alpha$

  This is a twisted involution if $s_\alpha$ twisted-commutes with |tw|;
  used in practice with even only with $\alpha$ imaginary for |tw|
*/
TwistedInvolution
ComplexReductiveGroup::reflection(RootNbr alpha,
				  const TwistedInvolution& tw) const
{
  WeylWord rw=rootDatum().reflectionWord(alpha);

  TwistedInvolution result=tw;
  for (size_t i=rw.size(); i-->0; ) // left multiply |tw| by |rw|
    W.leftMult(result,rw[i]);

  return result;
}


/******** manipulators *******************************************************/


/*
  The following function matches real forms between different Cartan classes,
  more precisely it sets the |Cartan[cn].real_labels| for all orbits in the
  adjoint fiber for Cartan |cn| to the corresponding internal real form
  number, identifying an orbit in the fundamental fiber of its inner class.

  This is achieved by first computing a |TorusPart| value |base|, that will,
  when combined with the canonical twisted involution |tw| for the Cartan
  class, give an adjoint Tits element that grades all simple-imaginary roots
  as noncompact. Once this is found, the sample torus parts stored for the
  different real forms over which this Cartan is defined are compared with
  this base torus part, the difference projected to the adjoint fiber group
  (which is a subquotient of the adjoint cocharacter lattice) at this Cartan
  by |toBasis|, and the class of that adjoint fiber group element in the weak
  real form partition of the current fiber selected; the number of that class
  is the index into the |real_labels| of this Cartan that is set to point to
  the real form used.

  Even though there must be an adjoint Tits element associated to the
  quasisplit real form that grades all simple-imaginary roots (at the
  canonical twisted involution for this Cartan) as noncompact, it might not
  have the torus part that was stored as sample in |Cartan[cn].rep[0]| (where
  $0$ is the internal |RealFormNbr| for the quasisplit form). Therefore we
  start by modifying the stored representative, by computing the grading
  |ref_gr| of the imaginary-simple roots the original choice |base| defines,
  looking up a representative adjoint fiber element |rep|, and interpreting it
  (by |fromBasis|) in the adjoint fiber group, and subtracting it from |base|
  so that the resulting |TorusPart| grades all simple-imaginary noncompact.
  The new value replaces the old one (though no other Atlas function uses the
  fact that the representatives for the quasisplit form are thus improved).

 */
void ComplexReductiveGroup::map_real_forms(CartanNbr cn)
{
  TitsCoset adj_Tg(*this);
  const Fiber& f = Cartan[cn].Cc.fiber();
  const Partition& weak_real = f.weakReal();
  const TwistedInvolution& tw = Cartan[cn].tw;
  const RootNbrList& sim = f.simpleImaginary();

  Cartan[cn].real_labels.resize(weak_real.classCount());

  TorusPart base = sample_torus_part(cn,quasisplit());
  TitsElt a(adj_Tg.titsGroup(),base,tw);
  Grading ref_gr; // reference simple-imaginary grading for quasisplit form
  for (size_t i=0; i<sim.size(); ++i)
    ref_gr.set(i,adj_Tg.grading(a,sim[i]));

  cartanclass::AdjointFiberElt rep = f.gradingRep(ref_gr);

  // now lift |rep| to a torus part and subtract from |base|
  SmallBitVector v(RankFlags(rep),f.adjointFiberRank());
  base -= f.adjointFiberGroup().fromBasis(v);
  // now |base| grades all imaginary simple roots noncompact
  Cartan[cn].rep[quasisplit()] = base.data(); // store improved representative

  for (BitMap::iterator rfi=Cartan[cn].real_forms.begin(); rfi(); ++rfi)
  {
    RealFormNbr rf = *rfi;
    TorusPart tp = sample_torus_part(cn,rf);
    cartanclass::AdjointFiberElt rep =
      f.adjointFiberGroup().toBasis(tp-=base).data().to_ulong();
    Cartan[cn].real_labels[weak_real.class_of(rep)]=rf;
  }
  assert(Cartan[cn].real_labels[0]==quasisplit());
} // |ComplexReductiveGroup::map_real_forms|

void ComplexReductiveGroup::map_dual_real_forms(CartanNbr cn)
{
  TitsCoset dual_adj_Tg(*this,tags::DualTag());
  const Fiber& dual_f = Cartan[cn].Cc.dualFiber();
  const Partition& dual_weak_real = dual_f.weakReal();
  const TwistedInvolution dual_tw =W.opposite(Cartan[cn].tw);
  const RootNbrList& sre = dual_f.simpleImaginary(); // simple real

  Cartan[cn].dual_real_labels.resize(dual_weak_real.classCount());
  assert(dual_weak_real.classCount()==Cartan[cn].dual_real_forms.size());

  TorusPart dual_base = dual_sample_torus_part(cn,0);
  TitsElt dual_a(dual_adj_Tg.titsGroup(),dual_base,dual_tw);
  Grading dual_ref_gr; // reference for dual quasisplit form
  for (size_t i=0; i<sre.size(); ++i)
    dual_ref_gr.set(i,dual_adj_Tg.grading(dual_a,sre[i]));

  cartanclass::AdjointFiberElt dual_rep = dual_f.gradingRep(dual_ref_gr);

  // now lift |rep| to a torus part and subtract from |base|
  SmallBitVector v(RankFlags(dual_rep),dual_f.adjointFiberRank());
  dual_base -= dual_f.adjointFiberGroup().fromBasis(v);

  for (BitMap::iterator
	 drfi=Cartan[cn].dual_real_forms.begin(); drfi(); ++drfi)
  {
    RealFormNbr drf = *drfi;
    TorusPart tp = dual_sample_torus_part(cn,drf);
    cartanclass::AdjointFiberElt rep =
      dual_f.adjointFiberGroup().toBasis(tp-=dual_base).data().to_ulong();
    Cartan[cn].dual_real_labels[dual_weak_real.class_of(rep)]=drf;
  }
  assert(Cartan[cn].dual_real_labels[0]==0);
  Cartan[cn].dual_rep[0] = dual_base.data();
} // |ComplexReductiveGroup::map_dual_real_forms|


// The size of the fiber orbits of strong forms for |rf| in Cartan |cn|
unsigned long
ComplexReductiveGroup::fiberSize(RealFormNbr rf, CartanNbr cn) const
{
  cartanclass::adjoint_fiber_orbit wrf = real_form_part(rf,cn);
  // |wrf| indexes a $W_{im}$ orbit on |cartan(cn).fiber().adjointFiberGroup()|

  const Fiber& f = cartan(cn).fiber();
  const cartanclass::StrongRealFormRep& srf = f.strongRealForm(wrf);

  assert(srf.second==f.central_square_class(wrf));

  // get partition of the fiber group according to action for square class
  const Partition& pi =f.fiber_partition(srf.second);
  // |pi| is an (unnormalized) partition of |cartan(cn).fiber().fiberGroup()|

  return pi.classSize(srf.first); // return size of orbit number |srf.first|
}

/*!
  \brief Returns the size of the dual fiber orbits corresponding to
  the dual strong real forms lying over dual real form \#rf, in Cartan
  \#cn.

  Precondition: real form \#rf is defined for cartan \#cn.
*/

unsigned long
ComplexReductiveGroup::dualFiberSize(RealFormNbr rf, CartanNbr cn) const
{
  cartanclass::adjoint_fiber_orbit wrf=dual_real_form_part(rf,cn);

  const Fiber& df = cartan(cn).dualFiber();
  const cartanclass::StrongRealFormRep& srf = df.strongRealForm(wrf);

  assert(srf.second==df.central_square_class(wrf));

  const Partition& pi = df.fiber_partition(srf.second);
  return pi.classSize(srf.first);
}



/*****************************************************************************

        Chapter III --- public accessor methods for ComplexReductiveGroup

******************************************************************************/


/******** accessors **********************************************************/

RankFlags ComplexReductiveGroup::simple_roots_imaginary() const
{
  const Permutation twist = simple_twist();
  RankFlags result;
  for (weyl::Generator s=0; s<semisimpleRank(); ++s)
    result.set(s,s==twist[s]);
  return result;
}

BitMap
ComplexReductiveGroup::Cartan_set(RealFormNbr rf) const
{
  BitMap support(Cartan.size());
  for (CartanNbr i=0; i<Cartan.size(); ++i)
    if (Cartan[i].real_forms.isMember(rf))
      support.insert(i);

  return support;
}

BitMap
ComplexReductiveGroup::dual_Cartan_set(RealFormNbr drf) const
{
  BitMap support(Cartan.size());
  for (CartanNbr i=0; i<Cartan.size(); ++i)
    if (Cartan[i].dual_real_forms.isMember(drf))
      support.insert(i);

  return support;
}

// The total number of involutions (generating Cartans as needed)
InvolutionNbr ComplexReductiveGroup::numInvolutions() const
{
  InvolutionNbr count = 0;

  for (CartanNbr cn=0; cn<numCartanClasses(); ++cn)
    count += cartan(cn).orbitSize();

  return count;
}

// The total number of involutions corresponding to the given set of Cartans
InvolutionNbr ComplexReductiveGroup::numInvolutions
  (const BitMap& Cartan_classes) const
{
  InvolutionNbr count = 0;

  for (BitMap::iterator it=Cartan_classes.begin(); it(); ++it)
    count += cartan(*it).orbitSize();

  return count;
}

WeightInvolution
ComplexReductiveGroup::compute_matrix(const TwistedInvolution& tw)  const
{
  return rootDatum().matrix(weylGroup().word(tw.w())) * distinguished();
}

// Sum of the real roots
Weight
ComplexReductiveGroup::posRealRootSum(const TwistedInvolution& tw) const
{
  return rootDatum().twoRho(involution_data(tw).real_roots());
}

// Sum of the imaginary roots.
Weight
ComplexReductiveGroup::posImaginaryRootSum(const TwistedInvolution& tw) const
{
  return rootDatum().twoRho(involution_data(tw).imaginary_roots());
}

/* Make |sigma| canonical and return Weyl group |w| element that
   twisted conjugates the canonical representative back to original |sigma|.

   We find conjugating generators starting at the original `|sigma|' end, so
   these form the letters of |w| from left (last applied) to right (first).
*/
WeylWord // return value is conjugating element
ComplexReductiveGroup::canonicalize
  (TwistedInvolution &sigma, // element to modify
   RankFlags gens) // subset of generators, "defaults" to all simple generators
  const
{
  return complexredgp::canonicalize(sigma,rootDatum(),twistedWeylGroup(),gens);
}

//!\brief find number of Cartan class containing twisted involution |sigma|
CartanNbr ComplexReductiveGroup::class_number(TwistedInvolution sigma) const
{
  canonicalize(sigma);
  for (CartanNbr i=0; i<Cartan.size(); ++i)
    if (sigma==Cartan[i].tw)
      return i;

  assert(false); // all canonical twisted involutions should occur in |Cartan|
  return ~0;
}


// the size of the subset of KGB for |rf| for involutions in |Cartan_classes|
unsigned long
ComplexReductiveGroup::KGB_size(RealFormNbr rf,
				const BitMap& Cartan_classes) const
{
  unsigned long result=0;
  for (BitMap::iterator it = Cartan_classes.begin(); it(); ++it)
    result +=  cartan(*it).orbitSize() * fiberSize(rf,*it);

  return result;

}

/* cardinality of the union of KGB sets  for this inner class.

  (each real form appears as often as there are strong real forms for it in
   its square class)
*/
unsigned long
ComplexReductiveGroup::global_KGB_size() const
{
  unsigned long result=0;
  for (CartanNbr cn=0; cn<numCartanClasses(); ++cn)
  {
    const CartanClass& cc = cartan(cn);
    result += cc.orbitSize() * cc.numRealFormClasses() * cc.fiber().fiberSize();
  }
  return result;

}

Grading ComplexReductiveGroup::simple_roots_x0_compact(RealFormNbr rf) const
{
#if 1 // actual code
  return
    RankFlags(fundamental().wrf_rep(rf)).unslice(simple_roots_imaginary());
#else // the |unslice| does the change of basis below more efficiently
  const Fiber& fund_f= fundamental();
  RankFlags noncompact(fund_f.wrf_rep(rf)); // among imaginary simples
  SmallBitVector v(noncompact,fund_f.adjointFiberRank());
  return fund_f.adjointFiberGroup().fromBasis(v).data();
#endif
}

cartanclass::square_class
ComplexReductiveGroup::xi_square(RealFormNbr rf) const
{ return fundamental().central_square_class(rf); }

RealFormNbr
ComplexReductiveGroup::square_class_repr(cartanclass::square_class csc) const
{ return fundamental().realFormPartition().classRep(csc); }

// here |diff| grades the simple roots, but only imaginary ones are involved
TorusPart ComplexReductiveGroup::grading_shift_repr (Grading diff) const
{
  const SmallSubquotient& fg = fundamental().fiberGroup();
  const unsigned int f_rank = fg.dimension();
  const unsigned int ssr = semisimpleRank();
  BinaryMap fg2grading (ssr,f_rank);
  for (RankFlags::iterator it=simple_roots_imaginary().begin(); it(); ++it)
  {
    const unsigned int i=*it;
    SmallBitVector ai (rootDatum().simpleRoot(i));
    for (unsigned int j=0; j<f_rank; ++j)
      fg2grading.set(i,j,ai.dot(fg.fromBasis(SmallBitVector(f_rank,j))));
  }
  const BinaryMap grading2fg = fg2grading.section();
  const SmallBitVector v = grading2fg*SmallBitVector(diff,ssr);
  assert((fg2grading*v).data()==diff);
  return fg.fromBasis(v);
}


RatCoweight ComplexReductiveGroup::base_grading_vector(RealFormNbr rf) const
{
  RankFlags cpct = simple_roots_x0_compact(square_class_repr(xi_square(rf)));
  return coch_representative(rootDatum(),cpct);
}

// an inverse operation to |coch_representative| composed with |exp_pi|
Grading compact_simple(const RootDatum& rd, RankFlags imag, TorusElement coch)
{
  Grading result;
  for (auto it=imag.begin(); it(); ++it)
    result.set(*it,coch.negative_at(rd.simpleRoot(*it)));
  return result;
}

/* find all FiberElt values in strong class of |y| for |csc|, mapping to
   |image| under |toAdjoint|, and return a list of |TorusPart| values that
   represent their differences with |y|
*/
containers::sl_list<TorusPart> preimage
(const Fiber& fund_f, const cartanclass::square_class csc,
 const cartanclass::FiberElt y, const cartanclass::AdjointFiberElt image)
{
  const SmallSubquotient& fg = fund_f.fiberGroup();
  const unsigned int f_rk = fund_f.fiberRank();
  const Partition& pi = fund_f.fiber_partition(csc);
  const cartanclass::fiber_orbit srf = pi.class_of(y);
  containers::sl_list<TorusPart> result;
  for (cartanclass::FiberElt i=pi.classRep(srf); i<pi.size(); ++i)
    if (pi.class_of(i)==srf and fund_f.toAdjoint(i)==image)
      result.push_back(fg.fromBasis(SmallBitVector(RankFlags(i^y),f_rk)));
  return result;
} // |preimage|

// torus parts that remain in the fiber and do not affect any grading
containers::sl_list<TorusPart>
ComplexReductiveGroup::central_fiber(RealFormNbr rf) const
{
  const Fiber& fund_f= fundamental();
  const cartanclass::square_class csc = fund_f.central_square_class(rf);
  const Partition& pi = fund_f.fiber_partition(csc);

  // find number of a fiber orbit that maps to our weak real form
  cartanclass::fiber_orbit srf;
  for (srf=0; srf<pi.classCount(); ++srf)
    if (fund_f.toWeakReal(srf,csc)==rf)
      break;

  const cartanclass::FiberElt y = pi.classRep(srf);

  return preimage(fund_f,csc,y,fund_f.toAdjoint(y));
} // |ComplexReductiveGroup::central_fiber|



TorusPart ComplexReductiveGroup::x0_torus_part(RealFormNbr rf) const
{
  Grading base = // find grading that corresponds to square class base
    square_class_grading(*this,xi_square(rf)); // noncompact marked
  base.complement(rootDatum().semisimpleRank()); // now compact is marked
  TorusPart t = grading_shift_repr(base^simple_roots_x0_compact(rf));
  TorusElement coch = // corresponds to (square class) base grading vector
    y_values::exp_pi(coch_representative(rootDatum(),base));
  coch += t; // now |coch| has implied grading for |rf|
  return
    t += minimum(central_fiber(rf),*this,coch); // to preferred orbit repr.tive
}



unsigned long
ComplexReductiveGroup::block_size(RealFormNbr rf,
				  RealFormNbr drf,
				  const BitMap& Cartan_classes) const
{
  unsigned long result=0;
  for (BitMap::iterator it = Cartan_classes.begin(); it(); ++it)
  {
    unsigned long cn=*it;
    result +=
      cartan(cn).orbitSize() * fiberSize(rf,cn) * dualFiberSize(drf,cn);
  }

  return result;

}



/*****************************************************************************

        Chapter IV -- Functions declared in complexredgp.h

******************************************************************************/


WeylWord canonicalize // return value conjugates element new |sigma| to old
  (TwistedInvolution& sigma,
   const RootDatum& rd,
   const TwistedWeylGroup& W,
   RankFlags gens)
{
  const RootNbrList s_image=W.simple_images(rd,sigma);
  InvolutionData id(rd,s_image);

  Weight rrs=rd.twoRho(id.real_roots());
  Weight irs=rd.twoRho(id.imaginary_roots());

/* the code below uses the following fact: if $S$ is a root subsystem of |rd|,
   and $\alpha$ a simple root that does not lie in $S$, then the sum of
   positive roots of $s_\alpha(S)$ is the image by $s_\alpha$ of the sum of
   positive roots of $S$. The reason is that the action of $s_\alpha$ almost
   preseves the notion of positivity; it only fails for the roots $\pm\alpha$,
   which do not occur in $S$ or in $s_\alpha(S)$. The code only applies
   $s_\alpha$ when the sum of positive of roots of $S$ is strictly
   anti-dominant for $\alpha$, and $S$ is either the system of real or
   imaginary roots; then $\alpha$ is a complex root, and in particular
   $\alpha$ does not lie in $S$.
 */

  WeylWord ww; // initialized empty; this will be the result

  { // first phase: make |rrs| dominant for all complex simple roots in |gens|
    // and make |irs| dominant for all such roots that are orthogonal to |rrs|
    RankFlags::iterator it; // allow inspection of final value
    do
      for (it=gens.begin(); it(); ++it)
      {
	weyl::Generator s=*it;
	LatticeCoeff c=rrs.dot(rd.simpleCoroot(s));
	if (c<0 or (c==0 and irs.dot(rd.simpleCoroot(s))<0))
	{
	  rd.simple_reflect(s,rrs);   // apply $s_i$ to real-root sum
	  rd.simple_reflect(s,irs);   // apply $s_i$ to iminary-root sum
	  W.twistedConjugate(sigma,s); // adjust |sigma| accordingly
	  ww.push_back(s);                 // and add generator to |ww|
	  break;     // after this change, continue the |do|-|while| loop
	}
      }
    while (it()); // i.e., until no change occurs any more
  }

/*
  Now that |rrs| and |irs| are dominant vectors, the simple coroots have non
  negative values on them [NOT TRUE: |irs| is only partly assured dominant!].
  Any positive coroot $\alpha$ is the sum of a multiset of simple coroots, and
  if $\alpha$ is orthogonal to |rrs| and |irs|, then its simple constituents
  must be so as well, since there can be no cancellation in the evaluations of
  $\alpha$ on |rrs|, nor on |irs|. Therefore the root subsystem orhogonal to
  |rrs| and |irs| is generated by a subset of the simple roots.
 */

  // clear those simple roots in |gens| not orthogonal to both |rrs| and |irs|
  for (RankFlags::iterator it=gens.begin(); it(); ++it)
    if (rrs.dot(rd.simpleCoroot(*it))>0 or irs.dot(rd.simpleCoroot(*it))>0)
      gens.reset(*it);


/* Now ensure that the involution |theta| associated to the twisted involution
   |sigma| fixes the dominant chamber for the root subsystem now flagged in
   |gens|, which we shall call the complex root subsystem. Since |theta|
   stablises this subsytem globally, this means it must be made to permute its
   positive and negative roots separately. We repeatedly inspect the simple
   roots of this subsystem, searching for some $\alpha_i$ that maps to a
   negative root; each time one is found, we twisted-conjugate |sigma| by $i$,
   which improves the situation (one can think of the twisted-conjugation as
   changing the positivity status of (only) $\alpha_i$, although the new
   statuses actually apply to the $\alpha_i$-reflected images of the roots).
   Eventually all positive roots in the subset will map to positive roots.
*/
  {
    RankFlags::iterator it;
    do
      for (it=gens.begin(); it(); ++it)
      {
	weyl::Generator s=*it;
	RootNbr beta= // image of |rd.simpleRootNbr(s)| by $\theta$
	  rd.permuted_root(W.word(sigma.w()),rd.simpleRootNbr(W.twisted(s)));
	if (rd.is_negroot(beta))
	{
	  W.twistedConjugate(sigma,s); // adjust |sigma|
	  ww.push_back(s);             // and add generator to |ww|
	  break;                       // and continue |do|-|while| loop
	}
      }
    while (it()); // i.e., while |for| loop was interrupted
  }

  return ww; // but the main result is the modfied value left in |sigma|
}

/*!
  \brief Puts into |so| the composite Cayley transform, and into |cross| the
   cross action corresponding to the twisted involution |ti|.

  Explanation: to each root datum involution $q$, we may associate a
  transformation from the fundamental involution to $q$ that factors as the
  composition of a cross action (conjugation by the inverse of an element
  |cross| of |W|), followed by a composite Cayley transform (composition at
  left or right with the product of (commuting) reflections for the roots of a
  strongly orthogonal set |so| of imaginary roots). This function computes
  |cross| and |so|, where $q$ is given by the twisted involution |ti|. Neither
  of the two parts of this decomposition are unique; for instance in the equal
  rank case the initial conjugation is entirely without effect.

  Note that conjugation is by the inverse of |cross| only because conjugation
  uses the letters of a Weyl word successively from right to left, whereas we
  collect the letters of |cross| from left to right as we go from the
  fundamental involution (the trivial twisted involution representing) back to
  $q$. While cross actions and Cayley transforms are found in an interleaved
  fashion, we push the latter systematically to the end (left) in the result;
  this means the corresponding roots are to be reflected by the cross actions
  that were found later, which cross actions themselves remain unchanged.
*/
void Cayley_and_cross_part(RootNbrSet& Cayley,
			   WeylWord& cross,
			   const TwistedInvolution& ti,
			   const RootSystem& rs,
			   const TwistedWeylGroup& W)
{
  weyl::InvolutionWord dec=W.involution_expr(ti);
  TwistedInvolution tw; // to reconstruct |ti| as a check

  RootNbrList so; // values for Cayley; |RootList| is preferable here
  cross.clear(); so.reserve(rs.rank());

  for (size_t j=dec.size(); j-->0; )
    if (dec[j]>=0) // Cayley transform by simple root
    {
      weyl::Generator s=dec[j];
      so.push_back(rs.simpleRootNbr(s));
      W.leftMult(tw,s);
    }
    else // cross action by simple root
    {
      weyl::Generator s=~dec[j];
      cross.push_back(s); // record cross action
      W.twistedConjugate(tw,s); // and twisted-conjugate |tw|
      // and conjugate roots in |so|:
      for (size_t i=0; i<so.size(); ++i)
	rs.simple_reflect_root(s,so[i]); // replace root by reflection image
    }

  assert(tw==ti);

  Cayley.set_capacity(rs.numRoots());
  Cayley.insert(so.begin(),so.end());
  Cayley = rs.long_orthogonalize(Cayley);
}

/* For a very long time real forms in the Atlas software were exclusively
   approached by selecting from a list that is presented after a rather
   elaborate preparatory calculation. The following function marks a different
   possibility, namely by allowing a real form (and in fact a Cartan class and
   a representative adjoint fiber element at that Cartan) to be selected based
   on specifying an involution and a grading of the imaginary roots. For
   practical reasons the grading is specified by its difference with the
   grading that makes noncompact all simple-imaginary roots (which grading
   belongs to the quasisplit real form) and this difference is specified by a
   rational coweight, whose pairing with any imaginary root should be integer,
   and its parity gives the mentioned difference in grading.
 */
RealFormNbr real_form_of // who claims this KGB element?
  (ComplexReductiveGroup& G,
   TwistedInvolution tw, const RatCoweight& grading_shift,
   CartanNbr& cn, WeylWord& conj, cartanclass::AdjointFiberElt& rep // outputs
   )
{
  const TwistedWeylGroup& W = G.twistedWeylGroup();
  const RootDatum& rd=G.rootDatum();

  { // test if |tw| is a proper twisted involution
    WeylElt test = tw;
    W.mult(test,W.twisted(tw));
    if (test!= WeylElt())
      throw std::runtime_error("Not a twisted involution");
  }

  // find the proper Cartan class
  conj = G.canonicalize(tw);
  for (cn=G.numCartanClasses(); cn-->0;)
    if (tw==G.involution_of_Cartan(cn))
      break;
  assert(cn!=~0); // every valid twisted involution should be found here

  // find the grading of the simple-imaginary roots at |tw|
  Grading gr;
  Coweight numer(grading_shift.numerator().begin(),
		 grading_shift.numerator().end()); // copy and convert
  arithmetic::Numer_t denom = grading_shift.denominator();
  rd.dual_act(numer,conj); // convert towards canonical involution
  InvolutionData id = InvolutionData::build(rd,W,tw);
  for (unsigned int i=0; i<id.imaginary_rank(); ++i)
  {
    int eval = rd.root(id.imaginary_basis(i)).dot(numer);
    if (eval%denom!=0)
      throw std::runtime_error
	("Rational coweight nonintegral at imaginary root");
    gr.set(i,eval%(2*denom)==0); // set means noncompact; happens if |eval==0|
  }

  // look up the grading
  const Fiber& f = G.cartan(cn).fiber();
  rep = f.gradingRep(gr);
  cartanclass::adjoint_fiber_orbit orb = f.weakReal().class_of(rep);
  return G.realFormLabels(cn)[orb]; // found our real form!!
} // |real_form_of|


// An attempt to improve the above
/*
   For a very long time real forms in the Atlas software were exclusively
   approached by selecting from a list that is presented after a rather
   elaborate preparatory calculation. The following function marks a different
   possibility, namely by allowing a real form to be selected based on
   specifying an involution and a rational coweight that is to be the
   'torus_factor' of some KGB element in the fiber over that involution. In
   fact the specified torus factor might belong to a different strong real
   form than the one that is implicitly used, so we allow this function to
   export a cocharacter |strong_form_start| that should be used to produce
   the base grading vector for this real form, rather than what the
   |base_grading_vector| method produces based on grading alone. In atlas the
   synthetic real form function is based on this functionality.

   We do insist that |strong_form_start| always produce the same grading (at
   the initial KGB element) as |base_grading_vector| for its weak real form
   does. Also we want any |tw,torus_factor| obtained from a KGB element for
   for a synthetic real form to reproduce the same |strong_form_start| value,
   which implies we must define an objective preference among candidates.
 */
RealFormNbr strong_real_form_of // who claims this KGB element?
  (ComplexReductiveGroup& G,
   TwistedInvolution tw, const RatCoweight& torus_factor,
   TorusElement& strong_form_start // additional output
   )
{
  const GlobalTitsGroup gTg(G);
  const Cartan_orbits& i_tab = G.involution_table();
  const Fiber& fund_f = G.fundamental();
  const RootDatum& rd = G.rootDatum();

  GlobalTitsElement x //  $\exp(\pi\ii torus_factor)$ gives |TorusElement|,
    (TorusElement(torus_factor,false),tw);
  assert(gTg.is_valid(x)); // unless this holds, we cannot hope to succeed

  const unsigned int r = G.semisimpleRank();
  { weyl::Generator s;
    while ((s=gTg.weylGroup().leftDescent(x.tw()))<r)
      if (gTg.hasTwistedCommutation(s,x.tw()))
	gTg.do_inverse_Cayley(s,x);
      else
	gTg.cross_act(s,x);

    x.torus_part().reduce();
    assert(gTg.is_valid(x)); // check that we still have a valid element
    assert (x.tw()==TwistedInvolution()); // and we are at fundamental fiber
  }


  RealFormNbr wrf; // the weak real form to return; find its value now:
  // must do this using the grading of simple-imaginary roots defined by |x|
  Grading gr; // grading of simple-imaginary roots
  G.generate_Cartan_orbit(0); // generate involutions in class of fundamental
  InvolutionNbr inv = i_tab.nr(x.tw()); // without doubt |inv==0|
  for (unsigned i=0; i<i_tab.imaginary_rank(inv); ++i)
    gr.set(i,
      not x.torus_part().negative_at(rd.root(i_tab.imaginary_basis(inv,i))));

  // found the grading, get a corresponding adjoint fiber element, and |wrf|
  cartanclass::AdjointFiberElt here = fund_f.gradingRep(gr);
  wrf = fund_f.weakReal().class_of(here);

  // to adjust the grading is subtle; our freedom depends on the square class
  cartanclass::square_class csc = fund_f.central_square_class(wrf);
  cartanclass::AdjointFiberElt goal = // element that gives x0 grading
    fund_f.wrf_rep(wrf); // get elected (in fact first) representative of |wrf|
  cartanclass::AdjointFiberElt base = fund_f.class_base(csc); // base of coset

  // to lift |here| to a |FiberElt|, must take it relative to coset base
  const SmallBitVector diff(RankFlags(here^base),fund_f.adjointFiberRank());
  cartanclass::FiberElt y = // get |toAdjoint|-preimage |FiberElt| of |diff|
    (fund_f.toAdjoint().section()*diff).data().to_ulong();

  // the orbit of |y| is a strong real form lying over |wrf|, stay within it

  auto pre = // find shifts within fiber orbit of |y| that will move to |goal|
    preimage(fund_f,csc,y,goal^base);

  strong_form_start = x.torus_part() += minimum(pre,G,x.torus_part());
  // now |x.torus_part| is what |torus_factor| should return at first KGB elt

  assert(compact_simple(rd,G.simple_roots_imaginary(),x.torus_part()) ==
	 G.simple_roots_x0_compact(wrf)); // check grading is as expected

  // take into account the torus bits that will be stored at KGB element x0
  strong_form_start += G.x0_torus_part(wrf); // morally subtraction, but mod 2
  to_minimum_representative(strong_form_start,G.distinguished());

  return wrf;

} // |strong_real_form_of|


RatCoweight coch_representative(const RootDatum& rd, Grading compact_simple)
{
  RatWeight result (rd.rank());
  for (RankFlags::iterator it=compact_simple.begin(); it(); ++it)
    result += rd.fundamental_coweight(*it);
  return result;
}

RatCoweight some_square // some value whose $\exp(2i\pi.)$ is in class |csc|
  (const ComplexReductiveGroup& G,cartanclass::square_class csc)
{
  const RootDatum& rd=G.rootDatum();
  auto gr = G.simple_roots_x0_compact(G.square_class_repr(csc));
  // got some grading associated to |csc|; get representative cocharacter for it
  return realredgp::square_class_choice // and standardise that |RatCoweight|
    (G.distinguished(),coch_representative(rd,gr));
}

// mark noncompact and complex simple roots according to |coch|
Grading grading_of_simples
  (const ComplexReductiveGroup& G, const RatCoweight& coch)
{
  const RootDatum& rd=G.rootDatum();
  Grading result; result.fill(rd.semisimpleRank()); // set complex simple roots
  for (auto it=G.simple_roots_imaginary().begin(); it(); ++it)
    result.set(*it,coch.dot(rd.simpleRoot(*it))%2==0); // set noncompact roots
  return result;
} // |grading_of_simples|

Grading square_class_grading(const ComplexReductiveGroup& G,
			     cartanclass::square_class csc)
{
  RatCoweight coch = some_square(G,csc);

  // finally convert that to a grading again
  const RootDatum& rd=G.rootDatum();
  Grading result; result.fill(rd.semisimpleRank()); // set complex simple roots
  for (auto it=G.simple_roots_imaginary().begin(); it(); ++it)
    result.set(*it,coch.dot(rd.simpleRoot(*it))%2==0); // set noncompact roots
  return result;
} // |square_class_grading|

/*****************************************************************************

        Chapter V -- Local Functions

******************************************************************************/

// Cross-transform the roots in |rl| according to |ww| (right to left)
void crossTransform(RootNbrList& rl,
		    const WeylWord& ww,
		    const RootSystem& rs)
{
  for (size_t i = ww.size(); i-->0;)
    rs.simple_root_permutation(ww[i]).renumber(rl);
}

/* Returns an element |x| (interpreted as element of the adjoint fiber
  of |fundf|) such that it grades the elements in |rl| according to |gr|.

  The successive bits of |gr| give the desired grading of the successive roots
  in |rl|. At least one solution for |x| should be known to exist. The roots
  in |rl| should probably be linearly independent, since linearly dependent
  roots would only make the existence of a solution less likely.
*/
unsigned long makeRepresentative(const Grading& gr,
				 const RootNbrList& rl,
				 const Fiber& fundf)
{
  RootNbrSet brs =
    fundf.noncompactRoots(0); // noncompact roots for the base grading
  Grading bgr =
    cartanclass::restrictGrading(brs,rl); // view as grading of roots of |rl|
  SmallBitVector bc(bgr,rl.size()); // transform to binary vector

  // make right hand side
  SmallBitVector rhs(gr,rl.size()); // view |gr| as binary vector (same length)
  rhs += bc;                        // and add the one for the base grading

  // make grading shifts
  SmallBitVectorList cl(fundf.adjointFiberRank(),bc);
  for (unsigned int i = 0; i < cl.size(); ++i)
  {
    Grading gr1 =
      cartanclass::restrictGrading(fundf.noncompactRoots(1 << i),rl);
    cl[i] += // cl[i] is shift for vector e[i]
      SmallBitVector(gr1,rl.size());
  }

  // set up equations
  RankFlags x;
  if (bitvector::combination_exists(cl,rhs,x))
    return x.to_ulong();
  else
    throw std::runtime_error("Representative of impossible grading requested");
} // |makeRepresentative|


// make from $(\Q/2\Z)^n$ to $(\Z/2\Z)^n$ taking floors of coefficients
SmallBitVector floor (Ratvec_Numer_t num, arithmetic::Numer_t d)
{
  SmallBitVector result (num.size()); // will hold the "floor" values
  for (unsigned int i=0; i<num.size(); ++i)
    result.set(i,num[i]%(2*d)/d!=0); // produce bit from floor
  return result;
}


// select |TorusPart| from list with minimal fingerprint when added to |coch|
TorusPart minimum(const containers::sl_list<TorusPart>& cf,
  const ComplexReductiveGroup& G, const TorusElement& coch)
{
  if (cf.size()==1) // trivial (but frequent) case, just return coch
    return cf.front();

  // compute fingerprinting matrix whose kernel annihilates fiber denominator
  const int_Matrix A = G.distinguished().transposed()+1;
  const int_Matrix projector = lattice::row_saturate(A);
  /* |projector| applies to |coch.log_pi(false)|, and the values in $\Q^k$
      produced are to be interpreted in $(\Q/2\Z)^k$. But |coch| is only to be
      shifted by elements of $H(2)$ from |cf|, giving integer vectors under
      |projector|, so for comparison of those we may as well apply the |floor|
      map $\Q/2\Z \to \Z/2\Z$ to the entries of the value so produced.
  */
  const RatCoweight& coch_log = coch.as_Qmod2Z();
  const SmallBitVector ref = // reference bitvector to which shifts are added
    floor(projector*coch_log.numerator(),coch_log.denominator());

  const BinaryMap pr(projector);
  auto it= cf.begin();
  unsigned long v=(ref + pr*(*it)).data().to_ulong(); // current minimal value
  TorusPart min = *it;
  for (++it; it!=cf.end(); ++it)
  {
    const unsigned long val = (ref + pr*(*it)).data().to_ulong();
    if (val<v)
    {
      v = val;
      min = *it;
    }
  }
  return min; // return the offset |TorusPart|, not the modified bitvector
} // |minimum|

// replace |coch| by minimum representative modulo image |delta_plus_1|
void to_minimum_representative(TorusElement& coch,
			       const WeightInvolution& delta)
{
  coch.right_symmetrise(delta); // ensure a delta-fixed representative
  BinaryMap A(delta.transposed()+1);
  const SmallSubspace mod_space(A); // image space of mod 2 (delta^t+1)
  const RatCoweight& coch_log = coch.as_Qmod2Z();
  const SmallBitVector ref = // reference bitvector to which shifts are added
    floor(coch_log.numerator(),coch_log.denominator());

  unsigned long min = ref.data().to_ulong();
  unsigned long min_i = 0;
  const unsigned d = mod_space.dimension();
  const unsigned long N=1ul<<d;
  for (unsigned long i=1; i<N; ++i)
  {
    const SmallBitVector v(RankFlags(i),d);
    const unsigned long val=(ref+mod_space.fromBasis(v)).data().to_ulong();
    if (val<min)
    {
      min = val;
      min_i = i;
    }
  }
  coch += mod_space.fromBasis(SmallBitVector(RankFlags(min_i),d));
}

// Modify |v| through through involution associated to |tw|
void twisted_act
  (const ComplexReductiveGroup& G, const TwistedInvolution& tw,Weight& v)
{
  G.involution_table().matrix(tw).apply_to(v);
}

void twisted_act
  (const ComplexReductiveGroup& G, const TwistedInvolution& tw,RatWeight& v)
{
  G.involution_table().matrix(tw).apply_to(v.numerator());
}

// Modify |v| through through involution associated to |tw|
void twisted_act
  (const ComplexReductiveGroup& G,Weight& v, const TwistedInvolution& tw)
{
  G.involution_table().matrix(tw).right_mult(v);
}

void twisted_act
  (const ComplexReductiveGroup& G,RatWeight& v, const TwistedInvolution& tw)
{
  G.involution_table().matrix(tw).right_mult(v.numerator());
}


#if 0 // functions below are no longer used

/* Add a new real form label list to d_realFormLabels.

  This is called when a new |CartanClass| has just been constructed at index
  |cn|. Then this function finds the labels corresponding to the real forms
  for which this Cartan is defined (the labelling of real forms being defined
  by the adjoint orbit picture in the fundamental fiber.)

  Algorithm: the gradings of the imaginary root system corresponding to the
  various real forms for which the new Cartan is defined are known. We find
  a cross-action followed by a composite Cayley transform, taking the
  fundamental Cartan to the new one. Then for each real form, we take a
  representative grading, and compute a grading for the fundamental Cartan
  transforming to it. This amounts to solving a system of linear equations
  mod 2.
*/
void ComplexReductiveGroup::correlateForms(CartanNbr cn)
{
  const RootSystem& rs = rootDatum();
  const TwistedWeylGroup& tW = twistedWeylGroup();
  const Fiber& fundf = d_fundamental;
  const Fiber& f = cartan(cn).fiber(); // fiber of this Cartan

  const TwistedInvolution& ti = twistedInvolution(cn);

  // find cayley part and cross part
  RootNbrSet so;
  WeylWord ww;
  Cayley_and_cross_part(so,ww,ti,rs,tW);

  assert(checkDecomposition(ti,ww,so,tW,rs));

  const Partition& pi = f.weakReal();
  RealFormNbrList rfl(f.numRealForms());

  // transform gradings and correlate forms
  for (size_t j = 0; j < rfl.size(); ++j)
  {
    unsigned long y = pi.classRep(j); //  orbit |j| has representative |y|
    Grading gr=f.grading(y); // and |gr| is it's grading
    RootNbrList rl = f.simpleImaginary(); // of the roots in |rl|

    gradings::transform_grading(gr,rl,so,rs); // grade those roots at |fundf|
    for (size_t i = 0; i < so.size(); ++i)
      gr.set(rl.size()+i);        // make grading noncompact for roots in |so|
    std::copy(so.begin(),so.end(),back_inserter(rl)); // extend |rl| with |so|
    crossTransform(rl,ww,rs);  // apply cross part of |ti| to roots in |rl|

    /* now |gr| grades the roots in |rl|,
       which are imaginary for the fundamental fiber |fundf| */
    for (size_t i = 0; i < rl.size(); ++i)
      assert(fundf.imaginaryRootSet().isMember(rl[i]));

    unsigned long x = makeRepresentative(gr,rl,fundf);
    RealFormNbr rf = fundf.weakReal()(x); // look up representative
    rfl[j] = rf;
  }

  Cartan[cn].real_labels = rfl;
  assert(rfl[0]==quasisplit()); // adjoint base grading is always quasisplit
} // |ComplexReductiveGroup::correlateForms|

void ComplexReductiveGroup::correlateDualForms(CartanNbr cn)
{
  const RootSystem& rs = dualRootDatum();
  const TwistedWeylGroup& tW = dualTwistedWeylGroup();
  const Fiber& fundf = d_dualFundamental;
  const Fiber& f = cartan(cn).dualFiber();

  TwistedInvolution ti = dualTwistedInvolution(cn);


  // find cayley part and cross part
  RootNbrSet so;
  WeylWord ww;
  Cayley_and_cross_part(so,ww,ti,rs,tW); // computation is in dual setting!

  assert(checkDecomposition(ti,ww,so,tW,rs));

  const Partition& pi = f.weakReal();
  RealFormNbrList rfl(f.numRealForms());

  // transform gradings and correlate forms
  for (size_t j=0; j<rfl.size(); ++j)
  {
    unsigned long y = pi.classRep(j);
    Grading gr=f.grading(y);
    RootNbrList rl = f.simpleImaginary();
    gradings::transform_grading(gr,rl,so,rs);
    for (size_t i=0; i<so.size(); ++i)
      gr.set(rl.size()+i);
    copy(so.begin(),so.end(),back_inserter(rl));
    crossTransform(rl,ww,rs);

// begin testing
    /* now |gr| grades the roots in |rl|,
       which are imaginary for the dual fundamental fiber |fundf| */
    for (size_t i=0; i<rl.size(); ++i)
      assert(fundf.imaginaryRootSet().isMember(rl[i]));
// end testing

    unsigned long x = makeRepresentative(gr,rl,fundf);
    RealFormNbr rf = fundf.weakReal()(x);
    rfl[j] = rf;
  }

  Cartan[cn].dual_real_labels = rfl;
  assert(rfl[0]==0); // adjoint base grading is always quasisplit
} // |ComplexReductiveGroup::correlateDualForms|

/*
  Check whether |ti| decomposes as the composition of the cross-action
  defined by |cross| followed by the Cayley transform defined by |Cayley|.
*/
bool checkDecomposition(const TwistedInvolution& ti,
			const WeylWord& cross,
			const RootNbrSet& Cayley,
			const TwistedWeylGroup& W,
			const RootSystem& rs)
{
  TwistedInvolution tw;

  // cross action part
  for (size_t i=0; i<cross.size(); ++i)
    W.twistedConjugate(tw,cross[i]);

  // cayley transform part
  for (RootNbrSet::iterator it=Cayley.begin(); it(); ++it)
  {
    InvolutionData id(rs,W.simple_images(rs,tw));
    assert(id.imaginary_roots().isMember(*it));
    W.leftMult(tw,rs.reflectionWord(*it));
  }

  return tw == ti;
}

#endif


} // |namespace complexredgp|

} // |namespace atlas|
