/*! \file \brief Implementation of the class CartanClassSet, storing
information about the set of stable conjugacy classes of Cartan subgroups.

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
*/
/*
  This is cartanset.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "cartanset.h"
#include "complexredgp.h"
#include "setutils.h"
#include <set>

#include <cassert>
#include "rootdata.h"
#include "tags.h"


namespace atlas {

namespace cartanset{

  void crossPart(weyl::WeylWord&, const weyl::WeylWord&,
		   const weyl::WeylGroup&);

  void cayleyPart(rootdata::RootList&, const weyl::WeylWord&,
		    const rootdata::RootDatum&, const weyl::WeylGroup&);

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
			  const weyl::WeylGroup&,
			  const rootdata::RootDatum&,
			  const latticetypes::LatticeMatrix& distinguished);

  const size_t UndefMostSplit = ~0ul;
  const realform::RealForm UndefRealForm = ~0ul;

} // namespace cartanset

/*****************************************************************************

        Chapter I --- The CartanClassSet class, construction

******************************************************************************/

namespace cartanset {

// the unique constructor
CartanClassSet::CartanClassSet
  (const complexredgp::ComplexReductiveGroup& parent,
   const latticetypes::LatticeMatrix& q)
  : d_parent(parent) // fix reference to parent inner class

  // install the fundamental Cartan, with associated data (like numRealForms)
  , d_cartan(1,new cartanclass::CartanClass(parent.rootDatum(),q))

  // initalise following structures to correspond to just fundamental Cartan
  , d_twistedInvolution(1,TwistedInvolution(weyl::WeylElt()))
  , d_ordering(1) // Cartans for a one point poset for now

  // the fundamental fiber can be copied from the fundamental Cartan
  , d_fundamental(d_cartan[0]->fiber())
  /* for the dual fundamental fiber it is tempting to give arguments either
     (d_cartan[0]->dualFiber()) or (parent.rootDatum(),q,tags::DualTag()).
     Both would give the same, wrong, result; one needs |dualBasedInvolution|.
     This fiber will be equal to the dual fiber for the most split Cartan.
  */
  , d_dualFundamental(rootdata::RootDatum(parent.rootDatum(),tags::DualTag())
		     ,dualBasedInvolution(q,parent.rootDatum()))

  // for real form labels install single list (filled later)
  , d_realFormLabels(1,realform::RealFormList(d_fundamental.numRealForms()))
  // for dual real forms, |correlateDualForms| will add such a list
  , d_dualRealFormLabels()

  /* for now each (dual) real form knows about one Cartan (filled in later) */
  , d_support(numRealForms(),bitmap::BitMap(1))
  , d_dualSupport(numDualRealForms(),bitmap::BitMap(1))

  , d_status(numRealForms())
  , d_mostSplit(numRealForms(),UndefMostSplit)
{
  /* since the fundamental Cartan is the reference for the numbering of real
     forms, each real form label for it just refers to itself */
  for (size_t i=0; i<numRealForms(); ++i)
    d_realFormLabels[0][i]=i;

  // for dual real forms, |correlateDualForms| adds the proper singleton list
  correlateDualForms(rootdata::RootDatum(parent.rootDatum(),tags::DualTag())
		    ,weyl::WeylGroup(parent.weylGroup(),tags::DualTag()));

  // mark the fundamental Cartan in each real form, and in one dual real form
  updateSupports(0); // take into account Cartan class number |0|

  /* mark the (compact) real form for which the fundamental Cartan is the most
     split one as done in |d_status|, and fill its value in |d_mostSplit| */
  updateStatus(0); // here the |0| means from Cartan |0| on.
}

/*!
  The only data explicitly allocated by the CartanClass is that for the
  CartanClass pointers.
*/
CartanClassSet::~CartanClassSet()
{
  for (size_t j = 0; j < d_cartan.size(); ++j)
    delete d_cartan[j];
}

/********* copy, assignment and swap *****************************************/

// This should remain empty!

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
void CartanClassSet::extend(realform::RealForm rf)
{
  if (d_status.isMember(rf)) // nothing to be done
    return;

  size_t prev_Cartans=d_cartan.size(); // mark size for possible roll-back
  bitmap::BitMap prev_status=d_status; // save this one, is hard to roll back

  try
  {
    const rootdata::RootDatum dualRootDatum(rootDatum(),tags::DualTag());
    const weyl::WeylGroup dualWeylGroup(weylGroup(),tags::DualTag());

    std::set<poset::Link> lks;

    // d_cartan.size() may grow in the course of the loop
    for (size_t j = 0; j < d_cartan.size(); ++j)
    {
      if (not isDefined(rf,j))
	continue;

      rootdata::RootSet rs=noncompactPosRootSet(rf,j);

      for (rootdata::RootSet::iterator it = rs.begin(); it(); ++it)
      {

	// multipy |twistedInvolution(j)| on the left by the reflection |*it|
	TwistedInvolution tw=reflection(*it,d_twistedInvolution[j]);

	/* check if we have a new twisted involution;
	   if so, extend the cartan structure */

	canonicalize(tw);
	size_t k = setutils::find_index(d_twistedInvolution,tw);
	lks.insert(std::make_pair(j,k)); // insert link in any case
	if (k == d_twistedInvolution.size()) // then class |k| is new
	{
	  d_twistedInvolution.push_back(tw); // record new canonical form
	  addCartan(tw); // and add corresponding Cartan class to |d_cartan|
	  correlateForms();
	  correlateDualForms(dualRootDatum,dualWeylGroup);
	  updateSupports(k); // take into account Cartan class number |k|
	}
      } // for (it)
    } // for (j)

    // update the order relation
    std::vector<poset::Link> lk(lks.begin(),lks.end());
    d_ordering.resize(d_cartan.size());
    d_ordering.extend(lk);

    // update |d_status| and |d_mostSplit| values for now completed real forms
    updateStatus(prev_Cartans);

  } // try

  catch (...) // ensure roll-back on (for instance) memory overflow
  { // the restoring code below avoids having to copy everything on entry
    d_cartan.resize(prev_Cartans);
    d_twistedInvolution.resize(prev_Cartans);
    d_ordering.resize(prev_Cartans);
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
  \brief Returns the various conjugacy classes of twisted involutions
  for the currently known Cartans.

  This method is NO LONGER USED during constructiont of a |CartanClassSet|.

  Note: the lists come out sorted, to allow binary look-up; this is in fact
  dependent on the implementation of |weyl::WeylGroup::twistedConjugacyClass|
*/
std::vector<weyl::WeylEltList> CartanClassSet::expand() const
{
  const weyl::WeylGroup& W = weylGroup();
  std::vector<weyl::WeylEltList> result(d_cartan.size());

  for (size_t j = 0; j < d_cartan.size(); ++j) {
    W.twistedConjugacyClass(result[j],d_twistedInvolution[j]);
  }

  return result;
}

/*!
  \brief Adds a new cartan to |d_cartan|, obtained from cartan \#j by
  Cayley transform through root \#rn.
  NO LONGER USED
*/
void CartanClassSet::addCartan(rootdata::RootNbr rn, size_t j)
{
  const rootdata::RootDatum& rd = rootDatum();
  latticetypes::LatticeMatrix q;
  rd.rootReflection(q,rn);
  q *= cartan(j).involution();
  d_cartan.push_back(new cartanclass::CartanClass(rd,q));
}


} // namespace cartanset

/*****************************************************************************

        Chapter II --- private (auxiliary) methods for CartanClassSet

******************************************************************************/

namespace cartanset {

/******** private accessors **************************************************/



/*!
  \brief Returns |tw| composed to the left with the reflection |s_rn|
  corresponding to root \#rn

  This is a twisted involution if |s_rn| twisted-commutes with |tw|;
  in practice root \#rn will in fact be imaginary for |tw|
*/
TwistedInvolution
CartanClassSet::reflection(rootdata::RootNbr rn,const TwistedInvolution& tw)
  const
{
  const weyl::WeylGroup& W = weylGroup();

  weyl::WeylWord rw=rootDatum().reflectionWord(rn);

  TwistedInvolution result=tw;
  for (size_t i=rw.size(); i-->0; ) // left multiply |tw| by |rw|
    W.leftMult(result,rw[i]);

  return result;
}

rootdata::RootSet CartanClassSet::noncompactPosRootSet
  (realform::RealForm rf, size_t j) const
/*!
  \brief Flags in rs the set of noncompact positive roots for Cartan \#j.
*/

{
  const cartanclass::Fiber& f = cartan(j).fiber();
  unsigned long x = CartanClassSet::representative(rf,j);

  rootdata::RootSet result=f.noncompactRoots(x);
  result &= rootDatum().posRootSet();
  return result;
}

/******** private manipulators **********************************************/

void CartanClassSet::correlateForms()

/*!
  \brief Adds a new real form label list to d_realFormLabels.

  Algorithm: the gradings of the imaginary root system corresponding to the
  various real forms for which the new cartan is defined are known. We find
  a cross-action followed by a composite Cayley transform, taking the
  fundamental cartan to the new one. Then for each real form, we take a
  representative grading, and compute a grading for the fundamental cartan
  transforming to it. This amounts to solving a system of linear equations
  mod 2.

  Explanation: this is called when a new cartan class has just been added
  to d_cartan. Then this function finds the labels corresponding to the
  real forms for which this cartan is defined (the labelling of real forms
  being defined by the adjoint orbit picture in the fundamental fiber.)
*/

{
  using namespace cartanclass;
  using namespace gradings;
  using namespace partition;
  using namespace rootdata;
  using namespace weyl;

  const RootDatum& rd = rootDatum();

  const WeylGroup& W = weylGroup();

  const Fiber& fundf = fundamental();
  const Fiber& f = d_cartan.back()->fiber();

  const TwistedInvolution& ti = d_twistedInvolution.back();

  // find cayley part and cross part
  RootList so;
  WeylWord ww;
  cayley_and_cross_part(so,ww,ti,rd,W);

  assert(checkDecomposition(ti,ww,so,W,rd,distinguished()));

  const Partition& pi = f.weakReal();
  realform::RealFormList rfl(f.numRealForms());

  // transform gradings and correlate forms
  for (size_t j = 0; j < rfl.size(); ++j) {
    unsigned long y = pi.classRep(j);
    Grading gr=f.grading(y);
    RootList rl = f.simpleImaginary(); // the roots of which |gr| are a grading

    transformGrading(gr,rl,so,rd);
    for (size_t i = 0; i < so.size(); ++i)
      gr.set(rl.size()+i);             // make grading set for roots in |so|
    copy(so.begin(),so.end(),back_inserter(rl)); // extend |rl| with |so|
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

void CartanClassSet::correlateDualForms(const rootdata::RootDatum& rd,
					const weyl::WeylGroup& W)
  // arguments are dual root datum and dual group
{
  using namespace cartanclass;
  using namespace gradings;
  using namespace latticetypes;
  using namespace partition;
  using namespace rootdata;
  using namespace weyl;

  const Fiber& fundf = dualFundamental();
  const Fiber& f = d_cartan.back()->dualFiber();

  /* find dual twisted involution by right-multiplication of |f| by |fundf|.
     However since |word_of_inverse_matrix| is used, we compute |q=fundf*f|. */
  LatticeMatrix q = fundf.involution();
  q *= f.involution();
  WeylWord tiww=rd.word_of_inverse_matrix(q);
  TwistedInvolution ti(WeylElt(tiww,W));

  // find cayley part and cross part
  RootList so;
  WeylWord ww;
  cayley_and_cross_part(so,ww,ti,rd,W);
// begin testing
  assert(checkDecomposition(ti,ww,so,W,rd,fundf.involution()));
// end testing

  const Partition& pi = f.weakReal();
  realform::RealFormList rfl(f.numRealForms());

  // transform gradings and correlate forms
  for (size_t j = 0; j < rfl.size(); ++j) {
    unsigned long y = pi.classRep(j);
    Grading gr=f.grading(y);
    RootList rl = f.simpleImaginary();
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

void CartanClassSet::updateStatus(size_t prev)

/*!
  \brief Updates d_status.

  Precondition: prev is the previous size of d_cartan;

  Explanation: d_status holds the subset of the set of real forms for
  which the full set of Cartan classes is constructed; equivalently, those
  for which the most split Cartan has been reached.
*/

{
  for (size_t j = prev; j<d_cartan.size(); ++j) // traverse new Cartan classes
  {

    const cartanclass::Fiber& f = cartan(j).fiber();
    const partition::Partition& pi = f.weakReal();
    const realform::RealFormList& rfl = realFormLabels(j);

    for (unsigned long c = 0; c < pi.classCount(); ++c) // traverse weak reals
    {
      if (cartan(j).isMostSplit(c))
      { // flag the form identified by class |c| in fiber of Cartan |j|
	d_status.insert(rfl[c]);
	d_mostSplit[rfl[c]] = j;
	// might do |break|, since no two real forms share a most split Cartan
      }
    }

  }

  return;
}

void CartanClassSet::updateSupports(size_t last)

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
  using namespace bitmap;

  for (size_t j = 0; j < d_support.size(); ++j)
    d_support[j].set_capacity(last+1);

  for (size_t j = 0; j < d_dualSupport.size(); ++j)
    d_dualSupport[j].set_capacity(last+1);

  const realform::RealFormList& rfl = realFormLabels(last);

  for (size_t j = 0; j < rfl.size(); ++j) {
    d_support[rfl[j]].insert(last);
  }

  const realform::RealFormList& drfl = dualRealFormLabels(last);

  for (size_t j = 0; j < drfl.size(); ++j) {
    d_dualSupport[drfl[j]].insert(last);
  }
}

void CartanClassSet::updateTwistedInvolutions
  (std::vector<weyl::WeylEltList>& known,
   const TwistedInvolution& tw)

/*!
  \brief Updates the known list by adding the twisted class of tw to it.
*/

{

  d_twistedInvolution.push_back(tw);

  weyl::WeylEltList wl; weylGroup().twistedConjugacyClass(wl,tw);

  known.push_back(weyl::WeylEltList());
  known.back().swap(wl);
}

/*****************************************************************************

        Chapter III --- public accessor methods for CartanClassSet

******************************************************************************/


/******** accessors **********************************************************/

/*!
  \brief Returns the size of the fiber orbits corresponding to the strong
  real forms lying over (weak) real form \#rf, in cartan \#cn.

  Precondition: Real form \#rf is defined for cartan \#cn.
*/
unsigned long CartanClassSet::fiberSize(realform::RealForm rf, size_t cn) const
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

unsigned long CartanClassSet::dualFiberSize(realform::RealForm rf, size_t cn)
  const
{
  cartanclass::adjoint_fiber_orbit wrf = dual_real_form_part(rf,cn);

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
size_t CartanClassSet::numInvolutions() const
{
  size_t count = 0;

  for (size_t cn = 0; cn < d_cartan.size(); ++cn)
    count += cartan(cn).orbitSize();

  return count;
}

/*!
  \brief Returns the total number of involutions corresponding to the
  indicated set of Cartans.
*/
size_t CartanClassSet::numInvolutions(const bitmap::BitMap& Cartan_classes)
  const
{
  size_t count = 0;

  for (bitmap::BitMap::iterator it=Cartan_classes.begin(); it(); ++it)
    count += cartan(*it).orbitSize();

  return count;
}

/*!
  \brief Returns a representative for real form \#rf in Cartan \#cn.

  Precondition: cartan \#cn occurs for that real form.

  Explanation: this amounts to searching for rf in d_realFormLabels[cn].
*/
unsigned long CartanClassSet::representative(realform::RealForm rf, size_t cn)
  const
{
  return cartan(cn).fiber().weakReal().classRep(real_form_part(rf,cn));
}

/*!
  \brief Returns a representative for dual real form \#rf in Cartan \#cn.

  Precondition: cartan \#cn occurs for that dual real form.

  Explanation: this amounts to searching for \#rf in d_dualRealFormLabels[cn].
*/
unsigned long CartanClassSet::dualRepresentative(realform::RealForm rf,
						size_t cn) const
{
  return
    cartan(cn).dualFiber().weakReal().classRep(dual_real_form_part(rf,cn));
}


latticetypes::LatticeMatrix
CartanClassSet::involutionMatrix(const weyl::TwistedInvolution& tw)
  const
{
  weyl::WeylWord ww;
  weylGroup().out(ww,tw.w());
  latticetypes::LatticeMatrix result;
  rootdata::toMatrix(result, ww,rootDatum());
  result *= distinguished();
  return result;
}

/*!
\brief Modify |v| through through involution associated to |tw|
*/
void CartanClassSet::twistedAct
  (const weyl::TwistedInvolution& tw,latticetypes::LatticeElt& v) const
{
  distinguished().apply(v,v);
  weylGroup().act(rootDatum(),tw.w(),v);
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
\brief Sum of the real roots.
*/
latticetypes::LatticeElt
  CartanClassSet::posRealRootSum(const TwistedInvolution& tw) const
{
  cartanclass::InvolutionData d(rootDatum(),involutionMatrix(tw));
  return rootDatum().twoRho(d.real_roots());
}

  /*!
\brief Sum of the imaginary roots.
  */
latticetypes::LatticeElt
  CartanClassSet::posImaginaryRootSum(const TwistedInvolution& tw) const
{
  cartanclass::InvolutionData d(rootDatum(),involutionMatrix(tw));
  return rootDatum().twoRho(d.imaginary_roots());
}

/*! \brief Make |sigma| canonical and return Weyl group |w| element that
    twisted conjugates the canonical representative back to original |sigma|.

    We find conjugating generators starting at the original `|sigma|' end, so
    these form the letters of |w| from left (last applied) to right (first).
*/
const weyl::WeylElt
CartanClassSet::canonicalize(TwistedInvolution &sigma) const
{
  const rootdata::RootDatum& rd=rootDatum();

/* the code below uses the following fact: if $S$ is a root subsystem of |rd|,
   and $\alpha$ a simple root that does not lie in $S$, then the sum of
   positive roots of $s_\alpha(S)$ is the image by $s_\alpha$ of the sum of
   positive roots of $S$. The reason is that the action of $s_\alpha$ almost
   preseves the notion of positivity; it only fails for the roots $\pm\alpha$,
   which do not occur in $S$ or in $s_\alpha(S)$. The code only applies
   $s_\alpha$ when the sum of positive of roots of $S$ is anti-dominant for
   $\alpha$, which excludes the case that $\alpha$ lies in $S$.
 */
  weyl::WeylElt w; // this will be the result
  latticetypes::LatticeElt rrs=posRealRootSum(sigma);

  {
    size_t i; // declare outside loop to allow inspection of final value
    do
      for (i=0; i<rd.semisimpleRank(); ++i)
	if (latticetypes::scalarProduct(rd.simpleCoroot(i),rrs) < 0 )
	{
	  rd.reflect(rrs,rd.simpleRootNbr(i));   // apply $s_i$ to root sum
	  weylGroup().twistedConjugate(sigma,i); // adjust |sigma| accordingly
	  weylGroup().mult(w,i);                 // and add generator to |w|
	  break;     // after this change, continue the |do|-|while| loop
	}
    while (i!=rd.semisimpleRank());    // i.e., until no change occurs any more
  }


/* now continue normalization, processing roots orthogonal to |rrs|. Since
   |rrs| is dominant, we can limit our attention to simple roots: any positive
   root orthogonal to |rrs| is a sum of simple roots orthogonal to |rrs|.
*/

  latticetypes::LatticeElt irs=posImaginaryRootSum(sigma);

  bitset::RankFlags simple_orth;

  for (size_t  i=0; i <  rd.semisimpleRank(); ++i)
    if (latticetypes::scalarProduct(rd.simpleCoroot(i),rrs) == 0)
      simple_orth.set(i);

  {
    bitset::RankFlags::iterator it;
    do
      for (it=simple_orth.begin(); it(); ++it)
	if (latticetypes::scalarProduct(rd.simpleCoroot(*it),irs) < 0)
	{
	  rd.reflect(irs,rd.simpleRootNbr(*it));   // apply $s_i$ to root sum
	  weylGroup().twistedConjugate(sigma,*it); // adjust |sigma|
	  weylGroup().mult(w,*it);                 // and add generator to |w|
	  break;           // after this change, continue the |do|-|while| loop
	}
    while (it()); // i.e., until no change occurs any more
  }


  /* Finally deal with simple roots orthogonal to |rrs| and to |irs| */


  // clear those simple roots in |simp_orth| not orthogonal to |irs|
  for (bitset::RankFlags::iterator it=simple_orth.begin(); it(); ++it)
    if (latticetypes::scalarProduct(rd.simpleCoroot(*it),irs) > 0)
      simple_orth.reset(*it);


/* Now ensure that the involution |theta| associated to the twisted involution
   |sigma| fixes the dominant chamber for the root subsystem indicated in
   |simple_orth| (which we shall call the complex root subsystem). The vector
   |x=rd.twoRho()| below is in the interior of the dominant chamber for the
   whole root system, so is a fortiori dominant of the subsystem. If its image
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
    {
      latticetypes::LatticeElt x=rd.twoRho(); // take fresh dominant each time
      twistedAct(sigma,x); // and (re)compute the effect of |sigma| into |x|

      for (it=simple_orth.begin(); it(); ++it)
	if (latticetypes::scalarProduct(rd.simpleCoroot(*it),x) < 0)
	{
	  weylGroup().twistedConjugate(sigma,*it); // adjust |sigma|
	  weylGroup().mult(w,*it);                 // and add generator to |w|
	  break;                             // and continue |do|-|while| loop
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
size_t CartanClassSet::classNumber(TwistedInvolution sigma) const
{
  canonicalize(sigma);
  return setutils::find_index(d_twistedInvolution,sigma);
}

/*!

  \brief returns index of canonical form of the product of twisted involution
  |d_twistedInvolution[j]| with the reflection through its |i|-th imaginary
  simple root. If |conjugator| is non-null, the conjugating element as returnd
  by |canonicalize| is assigned to |conjugator|.

  Note that the mentioned reflection twisted-commutes with
  |d_twistedInvolution[j]|, so that the product is again a twisted involution.
*/
size_t CartanClassSet::cayley(size_t j, size_t i, weyl::WeylElt* conjugator)
  const
{
  cartanclass::InvolutionData d
    (rootDatum(),involutionMatrix(d_twistedInvolution[j]));
  atlas::rootdata::RootNbr rn = d.imaginary_basis()[i];

  weyl::WeylWord rw=rootDatum().reflectionWord(rn);

  TwistedInvolution ti=d_twistedInvolution[j]; weylGroup().leftMult(ti,rw);

  if (conjugator==NULL) canonicalize(ti);
  else *conjugator=canonicalize(ti);

  return setutils::find_index(d_twistedInvolution,ti);
}


/*!
  \brief Returns the cardinality of the subset of \f$K\backslash G/B\f$
   associated to |rf| whose twisted involutions belong to |Cartan_classes|.

  Precondition: the Cartan classes for this real form have been generated
*/
unsigned long
CartanClassSet::KGB_size(realform::RealForm rf,
			 const bitmap::BitMap& Cartan_classes) const
{
  unsigned long result=0;
  for (bitmap::BitMap::iterator it = Cartan_classes.begin(); it(); ++it)
    result +=  cartan(*it).orbitSize() * fiberSize(rf,*it);

  return result;

}

unsigned long
CartanClassSet::block_size(realform::RealForm rf, realform::RealForm drf,
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

} // namespace cartanset




/*****************************************************************************

        Chapter IV --- Functions declared in cartanset.h

******************************************************************************/

namespace cartanset {

unsigned long blockSize(realform::RealForm rf, realform::RealForm drf,
			const CartanClassSet& ccl)

/*!
  \brief Returns the size of the block defined by the real form rf and the
  dual real form drf.

  NOTE: rf and drf are _weak_ real forms; the datum of the underlying weak
  real forms suffices to determine the structure of the block defined by a
  pair of strong real forms.
*/

{
  bitmap::BitMap b = ccl.support(rf);
  b &= ccl.dualSupport(drf);

  return ccl.block_size(rf,drf,b);
}

unsigned long kgbSize(realform::RealForm rf, const CartanClassSet& ccl)

/*!
  \brief Returns the cardinality of K\\G/B for this real form.

  Precondition: the Cartan classes for this real form have been generated

  Explanation: this is the cardinality of the one-sided parameter set
  corresponding to any strong real form over rf.
*/

{
  return ccl.KGB_size(rf,ccl.support(rf));
}

/*!
  \brief Puts into |so| the composite Cayley transform, and into |cross| the
   cross action corresponding to the twisted involution |ti|.

  Explanation: to each root datum involution |q|, we may associate a
  transformation from the fundamental involution to |q| that factors as the
  composition of a cross action (conjugation by the inverse of an element
  |cross| of |W|) and a composite Cayley transform (composition with
  (commuting) reflections for the roots of a strongly orthogonal set |so| of
  imaginary roots). This function computes |cross| and |so|, for the
  involution |q| corresponding to the twisted involution |ti|.

  Note that conjugation is by the inverse of |cross| only because conjugation
  uses the letters of a Weyl word successively from right to left, whereas we
  collect the letters of |cross| from left to right as we go from the (twisted
  ivolution representing) the fundamental involution back to |q|.
*/
void cayley_and_cross_part(rootdata::RootList& so,
			   weyl::WeylWord& cross,
			   const weyl::TwistedInvolution& ti,
			   const rootdata::RootDatum& rd,
			   const weyl::WeylGroup& W)
{
  std::vector<signed char> dec=W.involution_expr(ti);
  TwistedInvolution tw; // to reconstruct |ti| as a check

  so.clear();

  for (size_t j=dec.size(); j-->0; )
    if (dec[j]>=0)
    {
      weyl::Generator s=dec[j];
      so.push_back(rd.simpleRootNbr(s));
      W.leftMult(tw,s);
    }
    else
    {
      weyl::Generator s=~dec[j];
      cross.push_back(s); // record cross action
      W.twistedConjugate(tw,s); // and twisted-conjugate |tw|
      for (size_t i = 0; i < so.size(); ++i) // and conjugate roots in |so|:
	rd.rootReflect(so[i],s);
    }

  assert(tw==ti);
  rootdata::strongOrthogonalize(so,rd);
}

} // namespace cartanset

/*****************************************************************************

        Chapter V --- Auxiliary functions local to this module

******************************************************************************/

namespace cartanset {

void cayleyPart(rootdata::RootList& so,
		const weyl::WeylWord& wi,
		const rootdata::RootDatum& rd,
		const weyl::WeylGroup& W)

/*!
  \brief Puts in so the composite Cayley transform corresponding to wi.

  Precondition: |wi| is a reduced expression for a twisted involution |tw|;

  Explanation: to each root datum involution, we may associate a transformation
  from the fundamental Cartan to the current Cartan, that factors as the
  composition of a cross action and a composite Cayley transform. This function
  puts in |so| a system of strongly orthogonal roots representing the Cayley
  transform.

  The interpretation of the result is that the involution corresponding to
  |tw| can be obtained from a conjugate of the distinguished involution by
  multiplication by the (commuting) reflections for the roots in the strongly
  orthogonal set. The element by which the distinguished involution should be
  conjugated (more precisely its inverse) is given by |crossPart(.,wi,W)|
*/

{
  using namespace rootdata;
  using namespace weyl;

  TwistedInvolution tw;
  so.clear();

  for (size_t j = 0; j < wi.size(); ++j) {
    Generator s = wi[j];
    if (W.hasTwistedCommutation(s,tw)) { // add new root to |so|
      so.push_back(rd.simpleRootNbr(s));
      W.leftMult(tw,s);
    }
    else { // conjugate roots in |so|
      for (size_t i = 0; i < so.size(); ++i)
	so[i] = rd.rootPermutation(s)[so[i]];
      W.twistedConjugate(tw,s); // and twisted-conjugate |tw|
    }
  }

  strongOrthogonalize(so,rd);
}

void crossPart(weyl::WeylWord& ww, const weyl::WeylWord& wi,
		       const weyl::WeylGroup& W)

/*!
  \brief Puts in ww the cross action corresponding to wi.

  Explanation: to each root datum involution, we may associate a transformation
  from the fundamental cartan to the current cartan, that factors as the
  composition of a cross action and a composite Cayley transform. This function
  puts in ww a weyl word representing the cross action part.
*/

{
  using namespace weyl;

  TwistedInvolution tw;
  ww.clear();

  for (size_t j = 0; j < wi.size(); ++j) {
    Generator s = wi[j];
    if (W.hasTwistedCommutation(s,tw)) { // cayley transform
      W.leftMult(tw,s);
    } else { // cross action
      ww.push_back(s);
      W.twistedConjugate(tw,s);
    }
  }
}

void crossTransform(rootdata::RootList& rl,
		    const weyl::WeylWord& ww,
		    const rootdata::RootDatum& rd)

/*!
  \brief Cross-transforms the roots in rl according to ww.

  NOTE: the cross-transformations are done in _reverse_ order.
*/

{

  for (size_t j = ww.size(); j-->0;) {
    setutils::Permutation a = rd.rootPermutation(ww[j]);
    for (size_t i = 0; i < rl.size(); ++i)
      rl[i] = a[rl[i]];
  }
}

unsigned long makeRepresentative(const gradings::Grading& gr,
				 const rootdata::RootList& rl,
				 const cartanclass::Fiber& fundf)
/*!
  \brief Returns an element |x| (interpreted as element of the adjoint fiber
  of |fundf|) such that it grades the elements in |rl| according to |gr|.

  The successive bits of |gr| give the desired grading of the successive roots
  in |rl|. At least one solution for |x| should be known to exist. The roots
  in |rl| should probably be linearly independent, since linearly dependent
  roots would only make the existence of a solution less likely.
*/

{
  using namespace bitset;
  using namespace cartanclass;
  using namespace latticetypes;
  using namespace gradings;
  using namespace rootdata;

  RootSet brs =
    fundf.noncompactRoots(0); // noncompact roots for the base grading
  Grading bgr =
    restrictGrading(brs,rl);    // view it as a grading of the roots of |rl|
  SmallBitVector bc(bgr,rl.size()); // transform into binary vector

  // make right hand side
  SmallBitVector rhs(gr,rl.size()); // view |gr| as binary vector (same length)
  rhs += bc;                        // and add the one for the base grading

  // make grading shifts
  SmallBitVectorList cl(fundf.adjointFiberRank(),bc);
  for (size_t j = 0; j < cl.size(); ++j) {
    Grading gr1 = restrictGrading(fundf.noncompactRoots(1 << j),rl);
    cl[j] += SmallBitVector(gr1,rl.size()); // cl[j] is shift for vector e[j]
  }

  // set up equations
  RankFlags x;
  firstSolution(x,cl,rhs);

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
  noncompact root \f$\alpha\f$ in |so| is that the grading of \f$\beta\f$ in |rl| is
  flipped if and only if \f$\alpha+\beta\f$ is a root. So in all the grading of
  \f$\beta\f$ changes iff this condition is met for an odd number of \f$\alpha\f$s.
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
			const weyl::WeylGroup& W,
			const rootdata::RootDatum& rd,
			const latticetypes::LatticeMatrix& q)
{
  using namespace rootdata;
  using namespace weyl;

  TwistedInvolution tw;

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

} // namespace cartanset
} // namespace atlas
