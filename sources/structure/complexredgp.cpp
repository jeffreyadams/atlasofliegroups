/*!
\file
\brief Implementation for the class ComplexReductiveGroup.

  The ComplexReductiveGroup class will play a central role in the
  whole program.  Even though it is entirely defined by its based root
  datum and an involutive automorphism of that datum, it has seemed
  more natural to use this class to collect the wealth of
  combinatorial data that the root datum gives rise to, and that will
  serve as the basis for our description of the representation theory
  of G. Note that the current state of the theory, and most notably
  Vogan duality, makes it natural and necessary to consider all the
  real forms of our complex group (in a given inner class) at once; so
  that is another reason to not choose a real form a priori.
*/
/*
  This is complexredgp.cpp.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "complexredgp.h"

#include "tags.h"
#include "dynkin.h"
#include "lietype.h"
#include "poset.h"
#include "rootdata.h"
#include "tits.h"
#include "weyl.h"

#include <set>


/*****************************************************************************

  The ComplexReductiveGroup class will play a central role in the whole
  program; even though it is entirely defined by its root datum, it has
  seemed more natural to use this class to collect the wealth of
  combinatorial data that the root datum gives rise to, and that will
  serve as the basis for our description of the representation theory
  of G. Note that the current state of the theory, and most notably
  Vogan duality, makes it natural and necessary to consider all the
  real forms of our complex group (in a given inner class) at once;
  so that is another reason to not choose a real form a priori.

  This module, together with the cartanclass one, contains code for
  dealing with conjugacy classes of Cartan subgroups, and real
  forms. In this version, we have de-emphasized the classification of
  strong real forms, which is perhaps better treated as an
  input-output issue (it can certainly be recovered with a moderate
  effort from the data available here.) It appears that the data
  recorded here is what's most relevant to representation theory.

  The classification of real forms amounts to the classification of
  W^\delta- orbits in the fundamental fiber of the one-sided parameter
  space for the adjoint group (see the "combinatorics" paper on the
  Atlas website, or the forthcoming "algorithms" paper by Jeff Adams
  and Fokko du Cloux.); this is a very small computation, depending
  only on the Lie algebra.

  The classifications of conjugacy classes of Cartan subgroups, for the
  various real forms, all embed into the classification of conjugacy classes
  of root data involutions for the given inner class, equivalent to the
  classification of twisted involutions in the Weyl group, which is a purely
  combinatorial Weyl group computation. Fokko's original code proceeded more
  or less by listing each of these conjugacy classes, checking each new
  twisted involution for membership of each class of previously generated
  ones. Here we find for each conjugacy class of twisted involutions a
  canonical representative, which is relatively easy to compute. In this way
  we avoid enumeration of each conjugacy class.

  The classification of real forms can be recovered in this picture as
  well, by looking at the unique most split Cartan class for each real
  form.

  The most delicate part is the "correlation" part: for each Cartan, tell
  which orbit in the corresponding fiber corresponds to which real form, the
  real forms being labelled by the orbits in the fundamental fiber. In the
  current version, the solution to this problem is cleaner then previously,
  and perfectly general: it is obtained by writing out a system of equations
  for the grading defining the real form; this system does not always have a
  unique solution, but all solutions correspond to the same real form.


******************************************************************************/

namespace atlas {

namespace complexredgp{

  weyl::Twist make_twist(const rootdata::RootDatum& rd,
			 const latticetypes::LatticeMatrix& d);

  void crossTransform(rootdata::RootList&,
		      const weyl::WeylWord&, const rootdata::RootDatum&);

  unsigned long makeRepresentative(const gradings::Grading&,
				   const rootdata::RootList&,
				   const cartanclass::Fiber&);

  void transformGrading(gradings::Grading&,
			const rootdata::RootList&,
			const rootdata::RootList&,
			const rootdata::RootDatum&);

  bool isImaginary(const latticetypes::LatticeElt& v,
		   const weyl::TwistedInvolution& tw,
		   const weyl::WeylGroup&,
		   const rootdata::RootDatum&,
		   const latticetypes::LatticeMatrix& distinguished);

  bool checkDecomposition(const weyl::TwistedInvolution& ti,
			  const weyl::WeylWord& cross,
			  const rootdata::RootList& cayley,
			  const weyl::TwistedWeylGroup&,
			  const rootdata::RootDatum&,
			  const latticetypes::LatticeMatrix& distinguished);

  const size_t UndefMostSplit = ~0ul;
  const realform::RealForm UndefRealForm = ~0ul;

} // |namespace complexredgp|


/*****************************************************************************

        Chapter I -- The ComplexReductiveGroup class

******************************************************************************/

namespace complexredgp {

/*!
  \brief main constructor

  constructs a ComplexReductiveGroup from its root datum and a distinguished
  involution.

  Precondition: d is a based root datum involution: it globally fixes the set
  of positive roots

  NOTE: the ComplexReductiveGroup assumes ownership of the RootDatum pointed
  to by rd; users beware of this.
*/
ComplexReductiveGroup::ComplexReductiveGroup
 (const rootdata::RootDatum& rd, const latticetypes::LatticeMatrix& d)
  : d_rootDatum(rd)
  , d_fundamental(rd,d) // will also be fiber of cartan(0)
  , d_dualFundamental(rootdata::RootDatum(rd,tags::DualTag()),
		      dualBasedInvolution(d,rd)) // dual fiber of most split

  , d_titsGroup(rd,d,make_twist(rd,d))
  , root_twist(rd.root_permutation(simple_twist()))
  // install the fundamental Cartan, with associated data (like numRealForms)
  , d_cartan(1,new cartanclass::CartanClass(rd,d))
  // initalise following structures to correspond to just fundamental Cartan
  , d_twistedInvolution(1,weyl::TwistedInvolution(weyl::WeylElt()))

  // the fundamental fiber can be copied from the fundamental Cartan
  , d_mostSplit(numRealForms(),UndefMostSplit)
  , Cartan_poset(1) // Cartans for a one point poset for now

  /* for now each (dual) real form knows about one Cartan (filled in later) */
  , d_support(numRealForms(),bitmap::BitMap(1))
  , d_dualSupport(numDualRealForms(),bitmap::BitMap(1))
  // for real form labels install single list (filled later)

  , d_realFormLabels(1,realform::RealFormList(numRealForms()))
  // for dual real forms, |correlateDualForms| will add such a list
  , d_dualRealFormLabels()
  , d_status(numRealForms())
{
  /* since the fundamental Cartan is the reference for the numbering of real
     forms, each real form label for it just refers to itself */
  for (size_t i=0; i<d_realFormLabels[0].size(); ++i)
    d_realFormLabels[0][i]=i;

  // for dual real forms, |correlateDualForms| adds the proper singleton list
  correlateDualForms(rootdata::RootDatum(rootDatum(),tags::DualTag())
		    ,weyl::TwistedWeylGroup
		      (twistedWeylGroup(),tags::DualTag()));

  // mark the fundamental Cartan in each real form, and in one dual real form
  updateSupports(0); // take into account Cartan class number |0|

  /* mark the (compact) real form for which the fundamental Cartan is the most
     split one as done in |d_status|, and fill its value in |d_mostSplit| */
  updateStatus(0); // here the |0| means starting from Cartan |0|
}

/*!
  \brief constructs the complex reductive group dual to G.
*/
ComplexReductiveGroup::ComplexReductiveGroup(const ComplexReductiveGroup& G,
					     tags::DualTag)
  : d_rootDatum(G.rootDatum(),tags::DualTag())
  , d_fundamental(G.dualFundamental())
  // for the dual fundamental fiber we similarly copy a fiber, but from |G|
  , d_dualFundamental(G.fundamental())

  , d_titsGroup(d_rootDatum,d_fundamental.involution(),
		make_twist(d_rootDatum,d_fundamental.involution()))

  , root_twist(d_rootDatum.root_permutation(simple_twist()))
  // install the fundamental Cartan, with associated data (like numRealForms)
  , d_cartan(1,new cartanclass::CartanClass
	            (d_rootDatum,d_fundamental.involution()))

  , d_twistedInvolution(1,weyl::TwistedInvolution(weyl::WeylElt()))

  , d_mostSplit(numRealForms(),UndefMostSplit)

  , Cartan_poset(1) // Cartans for a one point poset for now

  /* for now each (dual) real form knows about one Cartan (filled in later) */
  , d_support(numRealForms(),bitmap::BitMap(1))
  , d_dualSupport(numDualRealForms(),bitmap::BitMap(1))
  // for real form labels install single list (filled later)
  , d_realFormLabels(1,realform::RealFormList(numRealForms()))
  // for dual real forms, |correlateDualForms| will add such a list
  , d_dualRealFormLabels()
  , d_status(numRealForms())
{
  /* since the fundamental Cartan is the reference for the numbering of real
     forms, each real form label for it just refers to itself */
  for (size_t i=0; i<d_realFormLabels[0].size(); ++i)
    d_realFormLabels[0][i]=i;

  // for dual real forms, |correlateDualForms| adds the proper singleton list
  correlateDualForms(G.rootDatum(),G.twistedWeylGroup());

  // mark the fundamental Cartan in each real form, and in one dual real form
  updateSupports(0); // take into account Cartan class number |0|

  /* mark the (compact) real form for which the fundamental Cartan is the most
     split one as done in |d_status|, and fill its value in |d_mostSplit| */
  updateStatus(0); // here the |0| means from Cartan |0| on.
}


/*!
  The only data explicitly allocated by the ComplexReductiveGroup is that
  for the CartanClass pointers.
*/
ComplexReductiveGroup::~ComplexReductiveGroup()
{
  for (size_t j = 0; j<d_cartan.size(); ++j)
    delete d_cartan[j];
}

/********* copy, assignment and swap *****************************************/

// This should remain empty!

} // namespace complexredgp

/*****************************************************************************

        Chapter II --- private (auxiliary) methods for ComplexReductiveGroup

******************************************************************************/

namespace complexredgp {

/******** private accessors **************************************************/

/*!
  \brief Returns |tw| composed to the left with the reflection |s_rn|
  corresponding to root \#rn

  This is a twisted involution if |s_rn| twisted-commutes with |tw|;
  in practice root \#rn will in fact be imaginary for |tw|
*/
weyl::TwistedInvolution
ComplexReductiveGroup::reflection(rootdata::RootNbr rn,
				  const weyl::TwistedInvolution& tw) const
{
  const weyl::WeylGroup& W = weylGroup();

  weyl::WeylWord rw=rootDatum().reflectionWord(rn);

  weyl::TwistedInvolution result=tw;
  for (size_t i=rw.size(); i-->0; ) // left multiply |tw| by |rw|
    W.leftMult(result,rw[i]);

  return result;
}


/******** manipulators *******************************************************/

/*!
  \brief Extends the CartanClassSet structure so that it contains all Cartans
  for the real form |rf|.

  This is done by traversing all the known Cartan classes |j| defined for the
  real form |rf| (which includes at least the distinguished Cartan for the
  inner class), and all non-compact positive imaginary roots \f$\alpha\f$ for
  |(rf,j)|; the twisted involution |tw| for the Cartan class |j| is then
  left-multiplied by \f$s_\alpha\f$ to obtain a new twisted involution |ti|,
  and if |ti| is not member of any of the known classes of twisted
  involutions, it starts a new Cartan class defined in |rf|, which will later
  itself be considered in the loop described here.

  While generating the Cartan classes, the ordering is extended by a link from
  the Cartan class |j| to the Cartan class of any twisted involution |ti|
  obtained directly from it.

*/
void ComplexReductiveGroup::extend(realform::RealForm rf)
{
  if (d_status.isMember(rf)) // nothing to be done
    return;

  size_t prev_Cartans=numCartanClasses(); // mark size for roll-back
  bitmap::BitMap prev_status=d_status; // save, hard to roll back

  try
  {
    const rootdata::RootDatum dualRootDatum(rootDatum(),tags::DualTag());
    const weyl::TwistedWeylGroup
      dualWeylGroup(twistedWeylGroup(),tags::DualTag());

    std::set<poset::Link> lks;

    // d_cartan.size() may grow in the course of the loop
    for (size_t j = 0; j<numCartanClasses(); ++j)
    {
      if (not is_defined(rf,j))
	continue;

      rootdata::RootSet rs=noncompactPosRootSet(rf,j);

      for (rootdata::RootSet::iterator it = rs.begin(); it(); ++it)
      {

	// multipy |twistedInvolution(j)| on the left by the reflection |*it|
	weyl::TwistedInvolution tw=reflection(*it,d_twistedInvolution[j]);

	/* check if we have a new twisted involution;
	   if so, extend the cartan structure */

	canonicalize(tw);
	size_t k = setutils::find_index(d_twistedInvolution,tw);
	lks.insert(std::make_pair(j,k)); // insert link in any case
	if (k == d_twistedInvolution.size()) // class |k| is new
	{
	  d_twistedInvolution.push_back(tw);
	     // record new canonical form
	  addCartan(tw); // and add corresponding Cartan class
	  correlateForms();
	  correlateDualForms(dualRootDatum,dualWeylGroup);
	  updateSupports(k); // incorporate Cartan class number |k|
	}
      } // for (it)
    } // for (j)

    // update the order relation
    std::vector<poset::Link> lk(lks.begin(),lks.end());
    Cartan_poset.resize(numCartanClasses());
    Cartan_poset.extend(lk);

    // update |d_status| and |d_mostSplit| values for now completed real forms
    updateStatus(prev_Cartans);

  } // try

  catch (...) // ensure roll-back on (for instance) memory overflow
  { // the restoring code below avoids having to copy everything on entry
    d_cartan.resize(prev_Cartans);
    d_twistedInvolution.resize(prev_Cartans);
    Cartan_poset.resize(prev_Cartans);
    d_realFormLabels.resize(prev_Cartans);
    d_dualRealFormLabels.resize(prev_Cartans);

    for (size_t i=0; i<numRealForms(); ++i)
    {
      d_support[i].set_capacity(prev_Cartans);
      d_dualSupport[i].set_capacity(prev_Cartans);
    }
    d_status.swap(prev_status); // restore set of completed real forms
    prev_status.andnot(d_status); // now |prev_status| flags "new" real forms
    for (bitmap::BitMap::iterator it=prev_status.begin(); it(); ++it)
      d_mostSplit[*it]=UndefMostSplit; // clear unsure most split Cartans
    throw; // rethrow same exception
  }
}

/*!
  \brief Adds a new real form label list to d_realFormLabels.

  This is called when a new Cartan class has just been added to |d_cartan|.
  Then this function finds the labels corresponding to the real forms for
  which this Cartan is defined (the labelling of real forms being defined by
  the adjoint orbit picture in the fundamental fiber.)

  Algorithm: the gradings of the imaginary root system corresponding to the
  various real forms for which the new Cartan is defined are known. We find
  a cross-action followed by a composite Cayley transform, taking the
  fundamental Cartan to the new one. Then for each real form, we take a
  representative grading, and compute a grading for the fundamental Cartan
  transforming to it. This amounts to solving a system of linear equations
  mod 2.
*/
void ComplexReductiveGroup::correlateForms()
{
  const rootdata::RootDatum& rd = rootDatum();

  const cartanclass::Fiber& fundf = d_fundamental;
  const cartanclass::Fiber& f =
    cartan(numCartanClasses()-1).fiber(); // fiber of Cartan just added

  const weyl::TwistedInvolution& ti = d_twistedInvolution.back();

  // find cayley part and cross part
  rootdata::RootList so;
  weyl::WeylWord ww;
  Cayley_and_cross_part(so,ww,ti);

  assert(checkDecomposition
	 (ti,ww,so,twistedWeylGroup(),rd,fundf.involution()));

  const partition::Partition& pi = f.weakReal();
  realform::RealFormList rfl(f.numRealForms());

  // transform gradings and correlate forms
  for (size_t j = 0; j < rfl.size(); ++j) {
    unsigned long y = pi.classRep(j);
    gradings::Grading gr=f.grading(y);
    rootdata::RootList rl =
      f.simpleImaginary(); // the roots of which |gr| are a grading

    transformGrading(gr,rl,so,rd);
    for (size_t i = 0; i < so.size(); ++i)
      gr.set(rl.size()+i);        // make grading noncompact for roots in |so|
    std::copy(so.begin(),so.end(),back_inserter(rl)); // extend |rl| with |so|
    crossTransform(rl,ww,rd);  // apply cross part of |ti| to roots in |rl|

    /* now |gr| grades the roots in |rl|,
       which are imaginary for the fundamental fiber |fundf| */
    for (size_t i = 0; i < rl.size(); ++i)
      assert(fundf.imaginaryRootSet().isMember(rl[i]));

    unsigned long x = makeRepresentative(gr,rl,fundf);
    realform::RealForm rf = fundf.weakReal()(x); // look up representative
    rfl[j] = rf;
  }

  d_realFormLabels.push_back(rfl);
}

void ComplexReductiveGroup::correlateDualForms(const rootdata::RootDatum& rd,
					       const weyl::TwistedWeylGroup& W)
  // arguments are dual root datum and dual group
{
  const cartanclass::Fiber& fundf = d_dualFundamental;
  const cartanclass::Fiber& f = cartan(numCartanClasses()-1).dualFiber();

  /* find dual twisted involution by right-multiplication of |f| by |fundf|.
     However since |word_of_inverse_matrix| is used, we compute |q=fundf*f|. */
  latticetypes::LatticeMatrix q = fundf.involution();
  q *= f.involution();
  weyl::WeylWord tiww=rd.word_of_inverse_matrix(q);
  weyl::TwistedInvolution ti(weyl::WeylElt(tiww,W));

  // find cayley part and cross part
  rootdata::RootList so;
  weyl::WeylWord ww;
  Cayley_and_cross_part(so,ww,ti);
// begin testing
  assert(checkDecomposition(ti,ww,so,W,rd,fundf.involution()));
// end testing

  const partition::Partition& pi = f.weakReal();
  realform::RealFormList rfl(f.numRealForms());

  // transform gradings and correlate forms
  for (size_t j = 0; j < rfl.size(); ++j) {
    unsigned long y = pi.classRep(j);
    gradings::Grading gr=f.grading(y);
    rootdata::RootList rl = f.simpleImaginary();
    transformGrading(gr,rl,so,rd);
    for (size_t i = 0; i < so.size(); ++i)
      gr.set(rl.size()+i);
    copy(so.begin(),so.end(),back_inserter(rl));
    crossTransform(rl,ww,rd);

// begin testing
    /* now |gr| grades the roots in |rl|,
       which are imaginary for the dual fundamental fiber |fundf| */
    for (size_t i = 0; i < rl.size(); ++i)
      assert(fundf.imaginaryRootSet().isMember(rl[i]));
// end testing

    unsigned long x = makeRepresentative(gr,rl,fundf);
    realform::RealForm rf = fundf.weakReal()(x);
    rfl[j] = rf;
  }

  d_dualRealFormLabels.push_back(rfl);
}


/*!
  \brief Updates d_status.

  Precondition: prev is the previous size of d_cartan;

  Explanation: d_status holds the subset of the set of real forms for
  which the full set of Cartan classes is constructed; equivalently, those
  for which the most split Cartan has been reached.
*/
void ComplexReductiveGroup::updateStatus(size_t prev)
{
  for (size_t j=prev; j<numCartanClasses(); ++j) // traverse new Cartan classes
  {

    const cartanclass::Fiber& f = cartan(j).fiber();
    const partition::Partition& pi = f.weakReal();
    const realform::RealFormList& rfl = d_realFormLabels[j];

    for (unsigned long c = 0; c < pi.classCount(); ++c) // traverse weak reals
    {
      if (cartan(j).isMostSplit(c))
      { // flag the form identified by class |c| in fiber of Cartan |j|
	d_status.insert(rfl[c]); // real form |rfl[c]| now fully generated
	d_mostSplit[rfl[c]] = j;
	// might do |break|, since no two real forms share a most split Cartan
      }
    }

  }
}

void ComplexReductiveGroup::updateSupports(size_t last)

/*!
  \brief Updates d_support and d_dualSupport.

  Explanation: for each real form rf, d_support[rf] contains the subset of the
  set of the currently defined cartans which are defined for rf; and
  analogously for d_dualSupport[rf]. This information is in fact a
  "cross-section" of the realFormLabels lists, which for each Cartan list the
  set of real forms for which the Cartan is defined. This function is called
  each time a new Cartan is added, and updates these lists.
*/

{
  for (size_t j = 0; j < d_support.size(); ++j)
    d_support[j].set_capacity(last+1);

  for (size_t j = 0; j < d_dualSupport.size(); ++j)
    d_dualSupport[j].set_capacity(last+1);

  const realform::RealFormList& rfl = d_realFormLabels[last];

  for (size_t j = 0; j < rfl.size(); ++j) {
    d_support[rfl[j]].insert(last);
  }

  const realform::RealFormList& drfl = d_dualRealFormLabels[last];

  for (size_t j = 0; j < drfl.size(); ++j) {
    d_dualSupport[drfl[j]].insert(last);
  }
}

/*****************************************************************************

        Chapter III --- public accessor methods for ComplexReductiveGroup

******************************************************************************/


/******** accessors **********************************************************/

/*!
  \brief Returns the size of the fiber orbits corresponding to the strong
  real forms lying over (weak) real form \#rf, in cartan \#cn.

  Precondition: Real form \#rf is defined for cartan \#cn.
*/
unsigned long
ComplexReductiveGroup::fiberSize(realform::RealForm rf, size_t cn) const
{
  cartanclass::adjoint_fiber_orbit wrf = real_form_part(rf,cn);
  // |wrf| indexes a $W_{im}$ orbit on |cartan(cn).fiber().adjointFiberGroup()|

  const cartanclass::Fiber& f = cartan(cn).fiber();
  const cartanclass::StrongRealFormRep& srf = f.strongRepresentative(wrf);

  assert(srf.second==f.central_square_class(wrf));

  const partition::Partition& pi =f.strongReal(srf.second);
  /* |pi| is an (unnormalized) partition of |cartan(cn).fiber().fiberGroup()|
     whose identifying values label central square classes, and are elements
     of the quotient affine space $adjoint fiber/image of fiber group$
  */

  size_t c = pi(srf.first);
  return pi.classSize(c);
}

/*!
  \brief Returns the size of the dual fiber orbits corresponding to
  the dual strong real forms lying over dual real form \#rf, in Cartan
  \#cn.

  Precondition: real form \#rf is defined for cartan \#cn.
*/

unsigned long
ComplexReductiveGroup::dualFiberSize(realform::RealForm rf, size_t cn) const
{
  cartanclass::adjoint_fiber_orbit wrf=dual_real_form_part(rf,cn);

  const cartanclass::Fiber& df = cartan(cn).dualFiber();
  const cartanclass::StrongRealFormRep& srf = df.strongRepresentative(wrf);

  assert(srf.second==df.central_square_class(wrf));

  const partition::Partition& pi = df.strongReal(srf.second);

  size_t c = pi(srf.first);

  return pi.classSize(c);
}

/*!
  \brief Returns the total number of involutions corresponding to the
  currently defined set of Cartans.
*/
size_t ComplexReductiveGroup::numInvolutions() const
{
  size_t count = 0;

  for (size_t cn = 0; cn<numCartanClasses(); ++cn)
    count += cartan(cn).orbitSize();

  return count;
}

/*!
  \brief Returns the total number of involutions corresponding to the
  indicated set of Cartans.
*/
size_t ComplexReductiveGroup::numInvolutions
  (const bitmap::BitMap& Cartan_classes) const
{
  size_t count = 0;

  for (bitmap::BitMap::iterator it=Cartan_classes.begin(); it(); ++it)
    count += cartan(*it).orbitSize();

  return count;
}


latticetypes::LatticeMatrix
ComplexReductiveGroup::involutionMatrix(const weyl::TwistedInvolution& tw)
  const
{
  latticetypes::LatticeMatrix result;
  rootdata::toMatrix(result,weylGroup().word(tw.w()),rootDatum());
  result *= distinguished();
  return result;
}

bool isImaginary (const latticetypes::LatticeElt& v,
		  const weyl::TwistedInvolution& tw,
		  const weyl::WeylGroup& W,
		  const rootdata::RootDatum& rd,
		  const latticetypes::LatticeMatrix& distinguished)
{
  latticetypes::LatticeElt x=v;
  distinguished.apply(x,x);
  W.act(rd,tw.w(),x);
  return x==v;
}

/*!
  \brief Flags in rs the set of noncompact positive roots for Cartan \#j.
*/
rootdata::RootSet ComplexReductiveGroup::noncompactPosRootSet
  (realform::RealForm rf, size_t j) const
{
  const cartanclass::Fiber& f = cartan(j).fiber();
  unsigned long x = representative(rf,j); // real form orbit-representative

  return f.noncompactRoots(x) & rootDatum().posRootSet();
}

/*!
\brief Sum of the real roots.
*/
latticetypes::LatticeElt
  ComplexReductiveGroup::posRealRootSum(const weyl::TwistedInvolution& tw)
  const
{
  cartanclass::InvolutionData d(*this,tw);
  return rootDatum().twoRho(d.real_roots());
}

  /*!
\brief Sum of the imaginary roots.
  */
latticetypes::LatticeElt
  ComplexReductiveGroup::posImaginaryRootSum(const weyl::TwistedInvolution& tw)
  const
{
  cartanclass::InvolutionData d(*this,tw);
  return rootDatum().twoRho(d.imaginary_roots());
}

/*! \brief Make |sigma| canonical and return Weyl group |w| element that
    twisted conjugates the canonical representative back to original |sigma|.

    We find conjugating generators starting at the original `|sigma|' end, so
    these form the letters of |w| from left (last applied) to right (first).
*/
const weyl::WeylElt // return value is conjugating element
ComplexReductiveGroup::canonicalize
  (weyl::TwistedInvolution &sigma, // element to modify
   bitset::RankFlags gens) // subset of generators
  const
{
  const rootdata::RootDatum& rd=rootDatum();
  const weyl::TwistedWeylGroup& W=twistedWeylGroup();
  weyl::WeylElt w; // this will be the result
  cartanclass::InvolutionData id(*this,sigma);

  latticetypes::LatticeElt rrs=rd.twoRho(id.real_roots());
  latticetypes::LatticeElt irs=rd.twoRho(id.imaginary_roots());

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

  { // first phase: make |rrs| dominant for all complex simple roots in |gens|
    // and make |irs| dominant for all such roots that are orthogonal to |rrs|
    bitset::RankFlags::iterator it; // allow inspection of final value
    do
      for (it=gens.begin(); it(); ++it)
      {
	size_t i=*it;
	latticetypes::LatticeCoeff c=rrs.scalarProduct(rd.simpleCoroot(i));
	if (c<0 or (c==0 and irs.scalarProduct(rd.simpleCoroot(i))<0))
	{
	  rd.reflect(rrs,rd.simpleRootNbr(i));   // apply $s_i$ to re-root sum
	  rd.reflect(irs,rd.simpleRootNbr(i));   // apply $s_i$ to im-root sum
	  W.twistedConjugate(sigma,i); // adjust |sigma| accordingly
	  W.mult(w,i);                 // and add generator to |w|
	  break;     // after this change, continue the |do|-|while| loop
	}
      }
    while (it()); // i.e., until no change occurs any more
  }


  /* Finally deal with simple roots orthogonal to |rrs| and to |irs| */


  // clear those simple roots in |simp_orth| not orthogonal to |irs|
  for (bitset::RankFlags::iterator it=gens.begin(); it(); ++it)
    if (rrs.scalarProduct(rd.simpleCoroot(*it))>0 or
	irs.scalarProduct(rd.simpleCoroot(*it))>0)
      gens.reset(*it);


/* Now ensure that the involution |theta| associated to the twisted involution
   |sigma| fixes the dominant chamber for the root subsystem indicated in
   |simple_orth| (which we shall call the complex root subsystem). The vector
   |x=rd.twoRho()| below is in the interior of the dominant chamber for the
   whole root system, so is a fortiori dominant for the subsystem. If its image
   is not dominant for the complex root subsystem, (twisted) conjugating
   |sigma| by any generator that makes the image more dominant will improve
   the situation (but not in the same way as the action of that generator: the
   image of |rd.twoRho()| under twisted action of |sigma| will get \emph{two}
   steps closer to the dominant chamber!), which we repeat until the image of
   |rd.twoRho()| becomes dominant for the complex root subsystem.
*/
  {
    bitset::RankFlags::iterator it;
    do
      for (it=gens.begin(); it(); ++it)
      {
	size_t i=*it;
	rootdata::RootNbr alpha=rd.simpleRootNbr(i);
	rootdata::RootNbr beta= // image of |alpha| by $\theta$
	  rd.permuted_root(W.word(sigma.w()),twisted_root(alpha));
	if (not rd.isPosRoot(beta))
	{
	  W.twistedConjugate(sigma,i); // adjust |sigma|
	  W.mult(w,i);                 // and add generator to |w|
	  break;                       // and continue |do|-|while| loop
	}
      }
    while (it()); // i.e., while |for| loop was interrupted
  }

  return  w; // but the main result is the modfied value left in |sigma|
}

  /*!
\brief find index of canonical representative of |sigma| in
  |d_twistedInvolution|, under the assumption that it is (already) present
  */
size_t ComplexReductiveGroup::class_number(weyl::TwistedInvolution sigma) const
{
  canonicalize(sigma);
  return setutils::find_index(d_twistedInvolution,sigma);
}


/*!
  \brief Returns the cardinality of the subset of \f$K\backslash G/B\f$
   associated to |rf| whose twisted involutions belong to |Cartan_classes|.

  Precondition: the Cartan classes for this real form have been generated
*/
unsigned long
ComplexReductiveGroup::KGB_size(realform::RealForm rf,
				const bitmap::BitMap& Cartan_classes) const
{
  unsigned long result=0;
  for (bitmap::BitMap::iterator it = Cartan_classes.begin(); it(); ++it)
    result +=  cartan(*it).orbitSize() * fiberSize(rf,*it);

  return result;

}

unsigned long
ComplexReductiveGroup::block_size(realform::RealForm rf,
				  realform::RealForm drf,
				  const bitmap::BitMap& Cartan_classes) const
{
  unsigned long result=0;
  for (bitmap::BitMap::iterator it = Cartan_classes.begin(); it(); ++it)
  {
    unsigned long cn=*it;
    result +=
      cartan(cn).orbitSize() * fiberSize(rf,cn) * dualFiberSize(drf,cn);
  }

  return result;

}

/*!
\brief Modify |v| through through involution associated to |tw|
*/
void ComplexReductiveGroup::twisted_act
  (const weyl::TwistedInvolution& tw,latticetypes::LatticeElt& v) const
{
  distinguished().apply(v,v);
  weylGroup().act(rootDatum(),tw.w(),v);
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
  |cross| and |so|, where $q$ is given by the twisted involution |ti|.

  Note that conjugation is by the inverse of |cross| only because conjugation
  uses the letters of a Weyl word successively from right to left, whereas we
  collect the letters of |cross| from left to right as we go from the (twisted
  involution representing) the fundamental involution back to $q$.
*/
  void ComplexReductiveGroup::Cayley_and_cross_part
    (rootdata::RootList& so,
     weyl::WeylWord& cross,
     const weyl::TwistedInvolution& ti) const
{
  const rootdata::RootDatum& rd=rootDatum();
  const weyl::TwistedWeylGroup& W=twistedWeylGroup();

  std::vector<signed char> dec=W.involution_expr(ti);
  weyl::TwistedInvolution tw; // to reconstruct |ti| as a check

  so.clear(); cross.clear();

  for (size_t j=dec.size(); j-->0; )
    if (dec[j]>=0) // Cayley transform by simple root
    {
      weyl::Generator s=dec[j];
      so.push_back(rd.simpleRootNbr(s));
      W.leftMult(tw,s);
    }
    else // cross action by simple root
    {
      weyl::Generator s=~dec[j];
      cross.push_back(s); // record cross action
      W.twistedConjugate(tw,s); // and twisted-conjugate |tw|
      for (size_t i = 0; i < so.size(); ++i) // and conjugate roots in |so|:
	rd.simple_reflect_root(so[i],s);
    }

  assert(tw==ti);
  rootdata::strongOrthogonalize(so,rd);
}


} // namespace complexredgp

/*****************************************************************************

        Chapter IV -- Functions declared in complexredgp.h

******************************************************************************/

namespace complexredgp {

/*!
  \brief puts in lt the Lie type of G.
*/
void lieType(lietype::LieType& lt, const ComplexReductiveGroup& G)
{
  lt.clear();

  const rootdata::RootDatum& rd = G.rootDatum();
  latticetypes::LatticeMatrix cm;

  rootdata::cartanMatrix(cm,rd);
  dynkin::lieType(lt,cm);

  // add the torus factor

  if (!rd.isSemisimple()) {
    lietype::SimpleLieType slt('T',rd.rank()-rd.semisimpleRank());
    lt.push_back(slt);
  }
}

/*****************************************************************************

        Chapter V -- Local Functions

******************************************************************************/

/*!\brief Returns the twist defined by |d| relative to |rd|.

  Precondition: |d| is an involution of the \emph{based} root datum |rd|.
*/
weyl::Twist make_twist(const rootdata::RootDatum& rd,
		       const latticetypes::LatticeMatrix& d)
{
  latticetypes::WeightList
    simple_roots(rd.beginSimpleRoot(),rd.endSimpleRoot());

  weyl::Twist result;

  for (size_t i = 0; i<simple_roots.size(); ++i)
  {
    result[i] = setutils::find_index(simple_roots,d.apply(simple_roots[i]));
    assert (result[i]<simple_roots.size());
  }

  return result;
}

/*!
  \brief Cross-transforms the roots in |rl| according to |ww|.

  NOTE: the cross-transformations of |ww| are done in right to left order.
*/
void crossTransform(rootdata::RootList& rl,
		    const weyl::WeylWord& ww,
		    const rootdata::RootDatum& rd)
{
  for (size_t j = ww.size(); j-->0;)
    rd.simple_root_permutation(ww[j]).left_mult(rl);
}

/*!
  \brief Returns an element |x| (interpreted as element of the adjoint fiber
  of |fundf|) such that it grades the elements in |rl| according to |gr|.

  The successive bits of |gr| give the desired grading of the successive roots
  in |rl|. At least one solution for |x| should be known to exist. The roots
  in |rl| should probably be linearly independent, since linearly dependent
  roots would only make the existence of a solution less likely.
*/
unsigned long makeRepresentative(const gradings::Grading& gr,
				 const rootdata::RootList& rl,
				 const cartanclass::Fiber& fundf)
{
  rootdata::RootSet brs =
    fundf.noncompactRoots(0); // noncompact roots for the base grading
  gradings::Grading bgr =
    cartanclass::restrictGrading(brs,rl); // view as grading of roots of |rl|
  latticetypes::SmallBitVector bc(bgr,rl.size()); // transform to binary vector

  // make right hand side
  latticetypes::SmallBitVector
                 rhs(gr,rl.size()); // view |gr| as binary vector (same length)
  rhs += bc;                        // and add the one for the base grading

  // make grading shifts
  latticetypes::SmallBitVectorList cl(fundf.adjointFiberRank(),bc);
  for (size_t j = 0; j < cl.size(); ++j) {
    gradings::Grading gr1 =
      cartanclass::restrictGrading(fundf.noncompactRoots(1 << j),rl);
    cl[j] += // cl[j] is shift for vector e[j]
      latticetypes::SmallBitVector(gr1,rl.size());
  }

  // set up equations
  bitset::RankFlags x;
  bitvector::firstSolution(x,cl,rhs);

  return x.to_ulong();
}

void transformGrading(gradings::Grading& gr,
		      const rootdata::RootList& rl,
		      const rootdata::RootList& so,
		      const rootdata::RootDatum& rd)

/*!
  \brief Transforms the grading |gr| of |rl| according to |so|.

  Precondition: |gr| is a grading of the roots in |rl|; |so| is a set of
  strongly orthogonal roots, orthogonal to the roots in |rl|.

  Assume that all roots in |rl| and |so| are imaginary for some involution,
  and that a grading is given that matches |gr| on |rl|, and for which all
  roots of |so| are noncompact (grading 1). Then by Cayley-transforming though
  the roots of |so| (in any order, since they are strongly orthogonal) one
  gets an involution for which the roots in |rl| are still imaginary (those in
  |so| have become real), and a grading of those roots that possibly differs
  from |gr|; this function transforms the grading |gr| into that new grading.

  Formula: the rule for Cayley-transforming a grading through an imaginary
  noncompact root \f$\alpha\f$ in |so| is that the grading of \f$\beta\f$ in
  |rl| is flipped if and only if \f$\alpha+\beta\f$ is a root. So in all the
  grading of \f$\beta\f$ changes iff this condition is met for an odd number
  of \f$\alpha\f$s.
*/

{
  for (size_t i = 0; i < rl.size(); ++i)
    for (size_t j = 0; j < so.size(); ++j)
      if (rd.sumIsRoot(rl[i],so[j]))
	gr.flip(i);
}

/*!
  \brief Checks whether |ti| decomposes as the composition of the cross-action
  defined by |ww| followed by the Cayley transform defined by |so|. Arguments
  |W| and |rd| are needed for the check; giving in addition |q| allows to
  check that Cayley transforms involve imaginary roots. This function is only
  used in an |assert| statement anyway.
*/
bool checkDecomposition(const weyl::TwistedInvolution& ti,
			const weyl::WeylWord& ww,
			const rootdata::RootList& so,
			const weyl::TwistedWeylGroup& W,
			const rootdata::RootDatum& rd,
			const latticetypes::LatticeMatrix& q)
{
  weyl::TwistedInvolution tw;

  // cross action part
  for (size_t j = 0; j < ww.size(); ++j)
    W.twistedConjugate(tw,ww[j]);

  // cayley transform part
  for (size_t j = 0; j < so.size(); ++j) {
    assert(isImaginary(rd.root(so[j]),tw,W,rd,q));
    W.leftMult(tw,rd.reflectionWord(so[j]));
  }

  return tw == ti;
}

} // |namespace complexredgp|

} // |namespace atlas|
