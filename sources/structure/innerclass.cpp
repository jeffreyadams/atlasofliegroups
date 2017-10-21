/*
  This is innerclass.cpp.

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/


/*
  Implementation for the class InnerClass.

  The InnerClass class will play a central role in the whole
  program. Even though it is entirely defined by its based root datum and an
  involutive automorphism of that datum, it has seemed more natural to use
  this class to collect the wealth of combinatorial data that they give rise
  to, and that will serve as the basis for our description of the
  representation theory of a real group G. Note that the current state of the
  theory, and most notably Vogan duality, makes it natural and necessary to
  consider all the real forms of our complex group (in a given inner class) at
  once; so that is another reason to not choose a real form a priori.
*/

#include "innerclass.h"
#include "matreduc.h" // |diagonalise| used in |square_class_choice|

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

namespace innerclass {

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


/*****************************************************************************

        Chapter I -- The InnerClass class

******************************************************************************/

// version of later |matrix| method that is usable during construction
inline WeightInvolution compute_matrix
  (const InnerClass& ic,const TwistedInvolution& tw)
{
  return ic.rootDatum().action_matrix(ic.weylGroup().word(tw.w()))
    * ic.distinguished();
}

InnerClass::C_info::C_info
  (const InnerClass& G,const TwistedInvolution twi, CartanNbr i)
  : tw(twi)
  , real_forms(G.numRealForms()), dual_real_forms(G.numDualRealForms())
  , rep(G.numRealForms()),        dual_rep(G.numDualRealForms())
  , below(i)
  , Cc(CartanClass(G.rootDatum(),G.dualRootDatum(),
		   compute_matrix(G,tw))) // generate fiber and dual fiber
  , real_labels(), dual_real_labels() // these start out emtpy
  {}


/*
  Main constructor

  Constructs an |InnerClass| from a pre-rootdatum |rd| and a
  distinguished involution |d|, which stabilises the set of simple roots
*/
InnerClass::InnerClass
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

  Constructs an |InnerClass| from a rootdatum |rd| and a
  distinguished involution |d|, which stabilises the set of simple roots
*/
InnerClass::InnerClass
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

void InnerClass::construct() // common part of two constructors
{
  { // task 1: generate Cartan classes, fill non-dual part of |Cartan|
    { // complete initialisation of |Cartan[0]|
      const Fiber& f=d_fundamental;
      const Partition& weak_real=weak_real_partition();
      // fill initial |form_reps| vector with assignment from |weak_real|
      for (RealFormNbr i=0; i<weak_real.classCount(); ++i)
      {
	Cartan[0].real_forms.insert(i); // Cartan 0 exists at all real forms
	// setting the initial torus part for each real form is subtle.
	Cartan[0].rep[i]=              // in |adj_Tg|, a torus part is a
	  cartanclass::restrictGrading // |Grading| of simple roots, where
	  (f.compactRoots(f.wrf_rep(i)), // compact ones need bit set
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

      const Fiber& f=d_dualFundamental;
      const Partition& weak_real=f.weakReal();
      // fill initial |form_reps| vector with assignment from |weak_real|
      for (unsigned long i=0; i<weak_real.classCount(); ++i)
      {
	// most split Cartan exists at all dual real forms
	Cartan.back().dual_real_forms.insert(i);

	RankFlags gr = // as above, torus part is obtained as a grading
	  cartanclass::restrictGrading // values at simple roots give torus part
	  (f.compactRoots(f.wrf_rep(i)), // compact ones need a set bit
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
} // |InnerClass::construct|

// Construct the complex reductive group dual to G
InnerClass::InnerClass(const InnerClass& G,
					     tags::DualTag)
  : d_rootDatum(G.d_dualRootDatum)
  , d_dualRootDatum(G.d_rootDatum)

  , my_W(nullptr), W(G.W) // not owned here, we depend on existence of |G|

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
InnerClass::~InnerClass()
{
  delete my_W;
}

/********* copy, assignment and swap *****************************************/

// This should remain empty!

/*****************************************************************************

        Chapter II --- private (auxiliary) methods for InnerClass

******************************************************************************/


/******** private manipulators ***************************************/


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
void InnerClass::map_real_forms(CartanNbr cn)
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

  // now subtract transformed coordinates of reference grading from |base|
  base -= f.adjointFiberGroup().fromBasis(f.gradingRep(ref_gr));
  // now |base| grades all imaginary simple roots noncompact
  Cartan[cn].rep[quasisplit()] = base.data(); // store improved representative

  for (BitMap::iterator rfi=Cartan[cn].real_forms.begin(); rfi(); ++rfi)
  {
    RealFormNbr rf = *rfi;
    TorusPart tp = sample_torus_part(cn,rf);
    cartanclass::AdjointFiberElt rep = f.adjointFiberGroup().toBasis(tp-=base);
    Cartan[cn].real_labels[f.adjoint_orbit(rep)]=rf;
  }
  assert(Cartan[cn].real_labels[0]==quasisplit());
} // |InnerClass::map_real_forms|

void InnerClass::map_dual_real_forms(CartanNbr cn)
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

  // now subtract transformed coordinates of reference grading from |dual_base|
  dual_base -=
    dual_f.adjointFiberGroup().fromBasis(dual_f.gradingRep(dual_ref_gr));

  for (BitMap::iterator
	 drfi=Cartan[cn].dual_real_forms.begin(); drfi(); ++drfi)
  {
    RealFormNbr drf = *drfi;
    TorusPart tp = dual_sample_torus_part(cn,drf);
    cartanclass::AdjointFiberElt dual_rep =
      dual_f.adjointFiberGroup().toBasis(tp-=dual_base);
    Cartan[cn].dual_real_labels[dual_f.adjoint_orbit(dual_rep)]=drf;
  }
  assert(Cartan[cn].dual_real_labels[0]==0);
  Cartan[cn].dual_rep[0] = dual_base.data();
} // |InnerClass::map_dual_real_forms|


// The size of the fiber orbits of strong forms for |rf| in Cartan |cn|
unsigned long
InnerClass::fiberSize(RealFormNbr rf, CartanNbr cn) const
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

/*
  Return the size of the dual fiber orbits corresponding to
  the dual strong real forms lying over dual real form \#rf, in Cartan \#cn.

  Precondition: real form \#rf is defined for cartan \#cn.
*/

unsigned long
InnerClass::dualFiberSize(RealFormNbr rf, CartanNbr cn) const
{
  cartanclass::adjoint_fiber_orbit wrf=dual_real_form_part(rf,cn);

  const Fiber& df = cartan(cn).dualFiber();
  const cartanclass::StrongRealFormRep& srf = df.strongRealForm(wrf);

  assert(srf.second==df.central_square_class(wrf));

  const Partition& pi = df.fiber_partition(srf.second);
  return pi.classSize(srf.first);
}



/*****************************************************************************

        Chapter III --- public accessor methods for InnerClass

******************************************************************************/


/******** accessors **********************************************************/

RankFlags InnerClass::simple_roots_imaginary() const
{
  const auto twist = twistedWeylGroup().twist();
  RankFlags result;
  for (weyl::Generator s=0; s<semisimpleRank(); ++s)
    result.set(s,s==twist[s]);
  return result;
}

RankFlags InnerClass::simple_roots_real() const
{
  const auto twist = twistedWeylGroup().dual_twist();
  RankFlags result;
  for (weyl::Generator s=0; s<semisimpleRank(); ++s)
    result.set(s,s==twist[s]);
  return result;
}

BitMap
InnerClass::Cartan_set(RealFormNbr rf) const
{
  BitMap support(Cartan.size());
  for (CartanNbr i=0; i<Cartan.size(); ++i)
    if (Cartan[i].real_forms.isMember(rf))
      support.insert(i);

  return support;
}

BitMap
InnerClass::dual_Cartan_set(RealFormNbr drf) const
{
  BitMap support(Cartan.size());
  for (CartanNbr i=0; i<Cartan.size(); ++i)
    if (Cartan[i].dual_real_forms.isMember(drf))
      support.insert(i);

  return support;
}

// The total number of involutions (generating Cartans as needed)
InvolutionNbr InnerClass::numInvolutions() const
{
  InvolutionNbr count = 0;

  for (CartanNbr cn=0; cn<numCartanClasses(); ++cn)
    count += cartan(cn).orbitSize();

  return count;
}

// The total number of involutions corresponding to the given set of Cartans
InvolutionNbr InnerClass::numInvolutions
  (const BitMap& Cartan_classes) const
{
  InvolutionNbr count = 0;

  for (BitMap::iterator it=Cartan_classes.begin(); it(); ++it)
    count += cartan(*it).orbitSize();

  return count;
}

// Sum of the real roots
Weight
InnerClass::posRealRootSum(const TwistedInvolution& tw) const
{
  return rootDatum().twoRho(involution_data(tw).real_roots());
}

// Sum of the imaginary roots.
Weight
InnerClass::posImaginaryRootSum(const TwistedInvolution& tw) const
{
  return rootDatum().twoRho(involution_data(tw).imaginary_roots());
}

/* Make |sigma| canonical and return Weyl group |w| element that left
   twisted-conjugates the canonical representative back to original |sigma|.

   We find conjugating generators starting at the original `|sigma|' value,
   and from that involution the canonical involution is right-conjugation
*/
WeylWord // letters in result will all be among those specified in |gens|
InnerClass::canonicalize
  (TwistedInvolution &sigma, // element to modify
   RankFlags gens) // subset of generators, "defaults" to all simple generators
  const
{
  const RootDatum& rd=rootDatum();
  const TwistedWeylGroup& tW=twistedWeylGroup();
  const WeylGroup& W=weylGroup();
  InvolutionData id(rd,tW.simple_images(rd,sigma));

  Weight rrs=rd.twoRho(id.real_roots()); // real (pos)root sum
  Weight irs=rd.twoRho(id.imaginary_roots()); // imaginary (pos)root sum

/*
  the code below uses the following fact: if $S$ is a root subsystem of |rd|,
  and $\alpha$ a simple root that does not lie in $S$, then the sum of
  positive roots of $s_\alpha(S)$ is the image by $s_\alpha$ of the sum of
  positive roots of $S$. The reason is that the action of $s_\alpha$ almost
  preseves the notion of positivity; it only fails for the roots $\pm\alpha$,
  which do not occur in $S$ or in $s_\alpha(S)$. The code only applies
  $s_\alpha$ when the sum of positive roots in $S$ is strictly anti-dominant
  for $\alpha$, where $S$ is first the system of real roots, and later the
  system of imaginary roots; in both cases $\alpha$ will always be a complex
  root, and in particular $\alpha$ does never lies in $S$.
*/

  WeylWord ww; // initialized empty; this will be the result

  { // first phase: make |rrs| dominant for all complex simple roots in |gens|
    // and make |irs| dominant for all such roots that are orthogonal to |rrs|
    RankFlags::iterator it; // allow inspection of final value
    do
      for (it=gens.begin(); it(); ++it)
      {
	weyl::Generator s=*it;
	LatticeCoeff c=rd.simpleCoroot(s).dot(rrs);
	if (c<0 or (c==0 and rd.simpleCoroot(s).dot(irs)<0))
	{
	  rd.simple_reflect(s,rrs);   // apply $s_i$ to real-root sum
	  rd.simple_reflect(s,irs);   // apply $s_i$ to imaginary-root sum
	  tW.twistedConjugate(sigma,s); // adjust |sigma| accordingly
	  ww.push_back(s);                 // and add generator to |ww|
	  break;     // after this change, continue the |do|-|while| loop
	}
      }
    while (it()); // i.e., until no change occurs any more
  }

/*
  For any weight $v$ \emph{dominant} for the subsystem |gens|, the further
  subsystem of coroots vanishing on $v$ is generated by its simple coroots. We
  apply this to conceptually reduce |gens|, after making |rrs| dominant for
  |gens|, to its members orthogonal to |rrs|, and then after making |irs|
  dominant for the remaining subsystem to its members that are also orthogonal
  to |irs|. The code below sets |gens| to generate the doubly reduced system.
*/
  // clear those simple roots in |gens| not orthogonal to both |rrs| and |irs|
  for (RankFlags::iterator it=gens.begin(); it(); ++it)
    if (rrs.dot(rd.simpleCoroot(*it))>0 or irs.dot(rd.simpleCoroot(*it))>0)
      gens.reset(*it);


/*
  Now ensure that the involution |theta| associated to the twisted involution
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
	RootNbr beta=rd.simpleRootNbr(tW.twisted(s)); // $\delta(\alpha_s)$
          W.act(rd,sigma.w(),beta); // now $\beta=\theta(\alpha_s)$
	if (rd.is_negroot(beta))
	{
	  tW.twistedConjugate(sigma,s); // adjust |sigma|
	  ww.push_back(s);             // and add generator to |ww|
	  break;                       // and continue |do|-|while| loop
	}
      }
    while (it()); // i.e., while |for| loop was interrupted
  }

  return ww; // but the main result is the modfied value left in |sigma|
} // |canonicalize|

// find the number of the Cartan class containing twisted involution |sigma|
CartanNbr InnerClass::class_number(TwistedInvolution sigma) const
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
InnerClass::KGB_size(RealFormNbr rf,
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
InnerClass::global_KGB_size() const
{
  unsigned long result=0;
  for (CartanNbr cn=0; cn<numCartanClasses(); ++cn)
  {
    const CartanClass& cc = cartan(cn);
    result += cc.orbitSize() * cc.numRealFormClasses() * cc.fiber().fiberSize();
  }
  return result;

}

// weak real forms are grouped by square class, whose index |xi_square| finds
// such forms share a central value for the square of their strong involutions
cartanclass::square_class InnerClass::xi_square // class for square
  (RealFormNbr rf) const
{ return d_fundamental.central_square_class(rf); }

RealFormNbr InnerClass::square_class_repr // class elected real form
  (cartanclass::square_class csc) const
{ return d_fundamental.realFormPartition().classRep(csc); }

// the basic information stored for a real form is a grading of the simple roots
Grading InnerClass::simple_roots_x0_compact(RealFormNbr rf) const
{
#if 1 // exploit that adjoint fiber group structure is rather transparent
  // Weak real forms are constructed as orbits in the adjoint fiber group. The
  // basis used for adjoint X_* is that of the fundamental coweights, on which
  // the distinguished involution acts as a permutation; the adjoint fiber
  // group generators are the fixed points of this (involutive) permutation,
  // which are the fundamental coweights for the imaginary simple coroots. So
  // each bit in a (weak) real form representative directly gives compactness
  // of a imaginary simple coroots; it sufices to move the bits to their place
  RankFlags tmp = d_fundamental.wrf_rep(rf).data(); // |unslice| is non-const
  return tmp.unslice(simple_roots_imaginary());
#else // |unslice| does the following change of basis below more efficiently
  const Fiber& fund_f= d_fundamental;
  RankFlags noncompact(fund_f.wrf_rep(rf)); // among imaginary simples
  SmallBitVector v(noncompact,fund_f.adjointFiberRank());
  return fund_f.adjointFiberGroup().fromBasis(v).data();
#endif
}


/*
  Given a real form cocharacter, reduce to one unique for its square class.

  This basically chooses a representative for the action of the fundamental
  fiber group on strong involutions in the square class, and serves to define
  a base point, to which the bits for the initial $x_0$ will be an offset.

  This means reduce a rational coweight $(X_*)_\Q^{\xi^t}$ modulo the
  sublattice $(X_*)^{\xi^t}$. We cannot just reduce modulo 1, as that may
  produce a non $\xi^t$-fixed coweight and indeed one from which one cannot
  easily recover the desired class of coweights. However once we express the
  coweights in coordinates relative to a basis of $(X_*)^{\xi^t}$, we can just
  reduce those coordinates modulo 1, and then convert back, so we do that.
*/
RatCoweight square_class_choice
  (const WeightInvolution& xi, const RatCoweight& coch)
{
  assert(coch==coch*xi); // assuming $\xi^t$-stable coweights

  int_Matrix fix = xi+1; // projection to $\xi$-fixed weights

  int_Matrix row,col;
  CoeffList diagonal = matreduc::diagonalise(fix,row,col);

  // initial columns of |col| are coordinate weights of $\xi^t$-stable coweights
  int_Matrix col_inv(col.inverse());
  const auto d=diagonal.size(), n=col.numRows();

  // convert to $\xi^t$-stable coordinates, take mod 1, then convert back
  RatCoweight tr = (coch*col.block(0,0,n,d)%=1)*col_inv.block(0,0,d,n);
  return tr.normalize();
}

// For a square class, we'll need to represent the set of compact roots for its
// |square_class_repr| by a square-central |TorusElement|,  via |compacts_for|:


// find compact ones among imaginary simple roots for |G|, as defined by |coch|
Grading compacts_for(const InnerClass& G, TorusElement coch)
{
  const RootDatum& rd=G.rootDatum();
  Grading result;
  for (auto it=G.simple_roots_imaginary().begin(); it(); ++it)
    result.set(*it,coch.negative_at(rd.simpleRoot(*it)));
  return result;
}

// The function |some_coch| defines an elected cocharacter for the square class

// some coweight $t$ with real form of $\exp(i\pi t)\delta_1$ in class |csc|
// this becomes |g_rho_check()| value for all non-synthetic real forms in |csc|
RatCoweight some_coch (const InnerClass& G,cartanclass::square_class csc)
{
  const RootDatum& rd=G.rootDatum();
  auto compacts = G.simple_roots_x0_compact(G.square_class_repr(csc));
  // got some grading associated to |csc|; get representative cocharacter for it
  RatWeight coch_rep (rd.rank());
  for (RankFlags::iterator it=compacts.begin(); it(); ++it)
    coch_rep += rd.fundamental_coweight(*it);
  assert(compacts_for(G,y_values::exp_pi(coch_rep))==compacts);
  return y_values::stable_log // for elected square root of sample square
    (y_values::exp_2pi(coch_rep),G.distinguished().transposed());
} // |some_coch|


// Preparations for finding an initial |TorusPart x0| in each weak real form

// We shall deduce |x0| from difference between set of compact roots (gradings)

// Find a |TorusPart|, assumed to exist, that induces a shift |diff| in grading
// here |diff| grades the simple roots (but only imaginary ones are considered)
TorusPart InnerClass::grading_shift_repr (Grading diff) const
{
  diff.slice(simple_roots_imaginary()); // discard complex roots
  const unsigned int c_rank=simple_roots_imaginary().count();

  const SmallSubquotient& fg = d_fundamental.fiberGroup();
  const SmallBitVectorList basis=fg.basis();
  const unsigned int f_rank = basis.size(); // equals |fg.dimension()|

  BinaryMap fg2grading(c_rank,f_rank);
  auto it=simple_roots_imaginary().begin();
  for (unsigned int i=0; i<c_rank; ++i)
  { SmallBitVector alpha_i(rootDatum().simpleRoot(*it++));
    for (unsigned int j=0; j<f_rank; ++j)
      fg2grading.set(i,j,basis[j].dot(alpha_i));
  }

  const BinaryMap grading2fg = fg2grading.section();
  const SmallBitVector v = grading2fg*SmallBitVector(diff,c_rank);
  assert((fg2grading*v).data()==diff); // check that we solved the equation
  return fg.fromBasis(v);
} // |grading_shift_repr|

/*
  |preimage| is a somewhat technical function to aid |central_fiber|, listing
  what is essentially the automorphism group of a (future) set K\G/B

  Given an element |y| of the fundamental fiber, find all those that (1) map
  under |toAdjoint| to the adjoint fiber element |image| (so that they will
  induce the corresponding shift in grading), and (2) lie in the "same strong
  real form" as |y|, this being an orbit in the fiber group under the action
  of the imaginary Weyl group, the action depending on the square class |csc|;
  return list of |TorusPart| values that represent their differences with |y|.
*/
containers::sl_list<TorusPart>
  InnerClass::torus_parts_for_grading_shift
    (const cartanclass::square_class csc,
     const cartanclass::FiberElt y, const cartanclass::AdjointFiberElt image)
  const
{
  const auto& fund_f = d_fundamental;
  const SmallSubquotient& fg = fund_f.fiberGroup();
  const unsigned int f_rk = fund_f.fiberRank();
  const Partition& pi = fund_f.fiber_partition(csc);
  const cartanclass::fiber_orbit srf = pi.class_of(y.data().to_ulong());
  containers::sl_list<TorusPart> result;
  for (unsigned i=pi.classRep(srf); i<pi.size(); ++i) // only start is optimised
    if (pi.class_of(i)==srf) // beyond that, just test for the right class
    { cartanclass::FiberElt fe(RankFlags(i),f_rk);
      if (fund_f.toAdjoint(fe)==image)
	result.push_back(fg.fromBasis(fe-y));
    }
  return result;
} // |torus_parts_for_grading_shift|

/*
  |preimage| is a somewhat technical function to aid |central_fiber|, listing
  what is essentially the automorphism group of a (future) set K\G/B

  Given an element |y| of the fundamental fiber, find all those that (1) map
  under |toAdjoint| to the adjoint fiber element |image| (so that they will
  induce the corresponding shift in grading), and (2) lie in the "same strong
  real form" as |y|, this being an orbit in the fiber group under the action
  of the imaginary Weyl group, the action depending on the square class |csc|;
  return list of |TorusPart| values that represent their differences with |y|.
*/
containers::sl_list<TorusPart> preimage
(const Fiber& fund_f, const cartanclass::square_class csc,
 const cartanclass::FiberElt y, const cartanclass::AdjointFiberElt image)
{
  const SmallSubquotient& fg = fund_f.fiberGroup();
  const unsigned int f_rk = fund_f.fiberRank();
  const Partition& pi = fund_f.fiber_partition(csc);
  const cartanclass::fiber_orbit srf = pi.class_of(y.data().to_ulong());
  containers::sl_list<TorusPart> result;
  for (unsigned i=pi.classRep(srf); i<pi.size(); ++i) // only start is optimised
    if (pi.class_of(i)==srf) // beyond that, just test for the right class
    { cartanclass::FiberElt fe(RankFlags(i),f_rk);
      if (fund_f.toAdjoint(fe)==image)
	result.push_back(fg.fromBasis(fe+y));
    }
  return result;
} // |preimage|

// now finding all fiber elements for |rf| with identical gradings is easy

// list stabiliser subgroup, in fundamental fiber group, of gradings for |rf|
containers::sl_list<TorusPart>
InnerClass::central_fiber(RealFormNbr rf) const
{
  const Fiber& fund_f= d_fundamental;
  const cartanclass::square_class csc = fund_f.central_square_class(rf);

  // find Fiber element |y| that maps to |fund_f.wrf_rep(rf)| in adjoint fiber
  const cartanclass::AdjointFiberElt target = fund_f.wrf_rep(rf);
  const SmallBitVector diff = target - fund_f.class_base(csc);
  const cartanclass::FiberElt y = fund_f.toAdjoint().section()*diff;
  assert(fund_f.toAdjoint(y)==diff); // check that we found a pre-image

  // and return shifts to other elements with same |toAdjoint| image as |y|
  return preimage(fund_f,csc,y,diff);
} // |InnerClass::central_fiber|


/*
  We need |central_fiber| in |x0_torus_part| only to standardise the choice,
  which the function |minimum| defined below actually accomplishes.

  What preceeds that is straightforward: start with the reference compacts for
  the square class; compare with the desired compacts at |x0| obtained from
  |simple_roots_x0_compact(rf)| to find the grading shift needed, and using
  |grading_shift_repr| find a |TorusPart| that will do this. This is the
  candidate to be standardised; torus element |coch| is needed for |minimum|.
  */
TorusPart InnerClass::x0_torus_part(RealFormNbr rf) const
{
  TorusElement t = // corresponds to (square class) base grading vector
    y_values::exp_pi(some_coch(*this,xi_square(rf)));
  const Grading base = compacts_for(*this,t);
  const Grading rf_cpt = simple_roots_x0_compact(rf);
  // mark compacts for elected cocharacter in square class
  // basic |x0| will be one that relative to base produces the right compacts
  TorusPart bits = grading_shift_repr(base^rf_cpt);

  auto cf = central_fiber(rf);
  assert(not cf.empty());

  auto it= cf.begin();
  auto min = *it;
  while (not (++it).at_end())
    if (bits+*it<bits+min) // offset |bits| does not cancel from this relation!
      min = *it;
  bits+=min;

  assert(tits::compact_simples(rootDatum(),t += bits,simple_roots_imaginary())
	 ==rf_cpt);

  return bits;
} // |x0_torus_part|



arithmetic::big_int
InnerClass::block_size(RealFormNbr rf, RealFormNbr drf,
		       const BitMap& Cartan_classes) const
{
  arithmetic::big_int result(0);
  for (BitMap::iterator it = Cartan_classes.begin(); it(); ++it)
  {
    CartanNbr cn=*it;
    result +=
      arithmetic::big_int::from_unsigned(cartan(cn).orbitSize())
      *arithmetic::big_int::from_unsigned(fiberSize(rf,cn))
      *arithmetic::big_int::from_unsigned(dualFiberSize(drf,cn));
  }

  return result;

}



/*****************************************************************************

        Chapter IV -- Functions declared in innerclass.h

******************************************************************************/


/*
  Put into |so| the composite Cayley transform, and into |cross| the cross
  action corresponding to the twisted involution |ti|.

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
} // |Cayley_and_cross_part|

/* |TitsCoset::TitsCoset| needs the status of all simple roots for |coch|.
  We mark noncompact roots (so this is complementary to |compacts_for|),
  and also mark complex roots that have an even pairing with |coch|.
*/
Grading grading_of_simples
  (const InnerClass& G, const RatCoweight& coch)
{
  const RootDatum& rd=G.rootDatum();
  Grading result;
  for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
    result.set(s,coch.dot(rd.simpleRoot(s))%2==0);
  return result;
}

RealFormNbr real_form_of // who claims this KGB element?
  (InnerClass& G, TwistedInvolution tw, // by value, modified
   const RatCoweight& torus_factor,
   RatCoweight& coch // additional output
  )
{
  const TwistedWeylGroup& W = G.twistedWeylGroup();
  const RootDatum& rd=G.rootDatum();
  const GlobalTitsGroup gTg(G);
  const WeightInvolution& xi = G.distinguished();

  GlobalTitsElement a //  $\exp(\pi\ii torus_factor)$ gives |TorusElement|,
    (TorusElement(torus_factor,false),tw);

  { // compute value of |coch|
    auto square=gTg.square_shifted(a); // yes it is shifted, we know ;-)
    LatticeMatrix simple_roots
      (rd.beginSimpleRoot(),rd.endSimpleRoot(),rd.rank(),tags::IteratorTag());
    assert(is_central(simple_roots,square)); // otherwise we cannot succeed

    coch = stable_log(square,xi.transposed()); // represents elected square root
  }

  auto conj = G.canonicalize(tw);
  gTg.cross_act(a,conj); // move to canonical twisted involution
  assert(a.tw()==tw); // we've moved |a| to the canonical involution

  // find the proper Cartan class
  CartanNbr cn;
  for (cn=G.numCartanClasses(); cn-->0;)
    if (tw==G.involution_of_Cartan(cn))
      break;
  assert(cn!=~0); // every valid twisted involution should be found here


  Grading gr; // will mark noncompact roots among simpl-imaginary ones at |a|
  { InvolutionData id = InvolutionData::build(rd,W,tw);
    for (unsigned int i=0; i<id.imaginary_rank(); ++i)
    {
      bool compact = a.torus_part().negative_at(rd.root(id.imaginary_basis(i)));
      gr.set(i,not compact); // set means noncompact; e.g., if |torus_part()==0|
    }
  }

  // look up the grading at the fiber of Cartan class |cn|
  const Fiber& f = G.cartan(cn).fiber();
  auto rep = f.gradingRep(gr);
  return G.realFormLabels(cn)[f.adjoint_orbit(rep)]; // found real form!
} // |real_form_of|

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
  cartanclass::AdjointFiberElt base(RankFlags(0),fundf.adjointFiberRank());
  RootNbrSet brs =
    fundf.noncompactRoots(base); // noncompact roots for the base grading
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
    cartanclass::AdjointFiberElt afe(fundf.adjointFiberRank(),i); // $e_i$
    Grading gr1 =
      cartanclass::restrictGrading(fundf.noncompactRoots(afe),rl);
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


// Modify |v| through through involution associated to |tw|
void twisted_act
  (const InnerClass& G, const TwistedInvolution& tw,Weight& v)
{
  G.involution_table().matrix(tw).apply_to(v);
}

void twisted_act
  (const InnerClass& G, const TwistedInvolution& tw,RatWeight& v)
{
  G.involution_table().matrix(tw).apply_to(v.numerator());
}

// Modify |v| through through involution associated to |tw|
void twisted_act
  (const InnerClass& G,Weight& v, const TwistedInvolution& tw)
{
  G.involution_table().matrix(tw).right_mult(v);
}

void twisted_act
  (const InnerClass& G,RatWeight& v, const TwistedInvolution& tw)
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
void InnerClass::correlateForms(CartanNbr cn)
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
} // |InnerClass::correlateForms|

void InnerClass::correlateDualForms(CartanNbr cn)
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
} // |InnerClass::correlateDualForms|

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


} // |namespace innerclass|

} // |namespace atlas|
