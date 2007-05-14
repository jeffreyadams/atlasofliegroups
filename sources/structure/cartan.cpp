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
  combinatorial Weyl group computation. This is not always a small
  computation, but can still be managed within a few seconds even for E8. The
  classification of real forms can be recovered in this picture as well, by
  looking at the unique most split Cartan class for each real form.

  The most delicate part is the "correlation" part: for each Cartan, tell
  which orbit in the corresponding fiber corresponds to which real form, the
  real forms being labelled by the orbits in the fundamental fiber. In the
  current version, the solution to this problem is cleaner then previously,
  and perfectly general: it is obtained by writing out a system of equations
  for the grading defining the real form; this system does not always have a
  unique solution, but all solutions correspond to the same real form.
*/
/*
  This is cartan.cpp.
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "cartan.h"

#include <set>

#include <cassert>
#include "rootdata.h"
#include "tags.h"


namespace atlas {

  namespace cartan{

  void crossPart(weyl::WeylWord&, const weyl::WeylWord&,
		   const weyl::WeylGroup&);

  void cayleyPart(rootdata::RootList&, const weyl::WeylWord&,
		    const rootdata::RootDatum&, const weyl::WeylGroup&);

  void crossTransform(rootdata::RootList&, const weyl::WeylWord&,
			const rootdata::RootDatum&);

  unsigned long makeRepresentative(const gradings::Grading&,
				   const rootdata::RootList&,
				   const cartanclass::Fiber&);

  void transformGrading(gradings::Grading&,
			const rootdata::RootList&,
			const rootdata::RootList&,
			const rootdata::RootDatum&);

  bool checkDecomposition(const weyl::WeylWord&, const weyl::WeylWord&,
			  const rootdata::RootList&, const weyl::WeylGroup&,
			  const rootdata::RootDatum&);

  const size_t UndefMostSplit = ~0ul;
  const realform::RealForm UndefRealForm = ~0ul;

} // namespace cartan

/*****************************************************************************

        Chapter I --- The CartanClassSet class, public methods

******************************************************************************/

namespace cartan {

// initial constructor
CartanClassSet::CartanClassSet
  (const complexredgp::ComplexReductiveGroup& parent,
   const latticetypes::LatticeMatrix& d)
  : d_parent(parent) // fix refernce to parent inner class

  /* install the fundamental Cartan with associated data (like numRealForms) */
  , d_cartan(1,new cartanclass::CartanClass(parent.rootDatum(),d))

  // initalise following structures to correspond to just fnudamental Cartan
  , d_twistedInvolution(1,TwistedInvolution(weyl::WeylElt()))
  , d_ordering(1) // Cartans for a one point poset for now

  /* the fundamental fiber can be copied from the fundamental Cartan */
  , d_fundamental(d_cartan[0]->fiber())
  /* for the dual fundamental fiber it is tempting to give arguments either
     (d_cartan[0]->dualFiber()) or (parent.rootDatum(),d,tags::DualTag()).
     Both would give the same, wrong, result; one needs |dualBasedInvolution|.
     This fiber will be equal to the dual fiber for the most split Cartan.
  */
  , d_dualFundamental(rootdata::RootDatum(parent.rootDatum(),tags::DualTag())
		     ,dualBasedInvolution(d,parent.rootDatum()))

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

  /* mark the (compact) real form for which the fundamental Cartan is the most
     split one as done in |d_status|, and fill its value in |d_mostSplit| */
  updateStatus(0); // here the |0| means from Cartan |0|

  // mark the fundamental Cartan in each real form, and in one dual real form
  updateSupports(0); // take into account Cartan class number |0|
}

CartanClassSet::~CartanClassSet()

/*!
  The only data explicitly allocated by the CartanClass is that for the
  CartanClass pointers.
*/

{
  for (size_t j = 0; j < d_cartan.size(); ++j)
    delete d_cartan[j];
}

/********* copy, assignment and swap *****************************************/

// This should remain empty!

/******** accessors **********************************************************/

unsigned long CartanClassSet::dualFiberSize(realform::RealForm rf, size_t cn)
  const
/*!
  \brief Returns the size of the dual fiber orbits corresponding to
  the dual strong real forms lying over dual real form \#rf, in cartan
  \#cn.

  Precondition: real form \#rf is defined for cartan \#cn.
*/

{
  using namespace cartanclass;
  using namespace partition;

  const realform::RealFormList& rfl = d_dualRealFormLabels[cn];
  size_t p = std::find(rfl.begin(),rfl.end(),rf) - rfl.begin();

  const Fiber& df = cartan(cn).dualFiber();
  const StrongRealFormRep& srf = df.strongRepresentative(p);
  const Partition& pi = df.strongReal(srf.second);

  size_t c = pi(srf.first);

  return pi.classSize(c);
}

unsigned long CartanClassSet::dualRepresentative(realform::RealForm rf,
						size_t cn) const

/*!
  \brief Returns a representative for dual real form \#rf in Cartan \#cn.

  Precondition: cartan \#cn occurs for that dual real form.

  Explanation: this amounts to searching for \#rf in d_dualRealFormLabels[cn].
*/

{
  using namespace partition;

  const realform::RealFormList& rfl = d_dualRealFormLabels[cn];
  size_t p = std::find(rfl.begin(),rfl.end(),rf) - rfl.begin();
  const Partition& pi = cartan(cn).dualFiber().weakReal();

  return pi.classRep(p);
}

unsigned long CartanClassSet::fiberSize(realform::RealForm rf, size_t cn) const

/*!
  \brief Returns the size of the fiber orbits corresponding to the strong
  real forms lying over real form \#rf, in cartan \#cn.

  Precondition: Real form \#rf is defined for cartan \#cn.
*/

{
  using namespace cartanclass;
  using namespace partition;

  const realform::RealFormList& rfl = d_realFormLabels[cn];
  size_t p = std::find(rfl.begin(),rfl.end(),rf) - rfl.begin();

  const Fiber& f = cartan(cn).fiber();
  const StrongRealFormRep& srf = f.strongRepresentative(p);
  const Partition& pi =f.strongReal(srf.second);

  size_t c = pi(srf.first);

  return pi.classSize(c);
}

size_t CartanClassSet::numInvolutions() const

/*!
  \brief Returns the total number of involutions corresponding to the
  currently defined set of cartans.
*/

{
  size_t count = 0;

  for (size_t cn = 0; cn < d_cartan.size(); ++cn)
    count += cartan(cn).orbitSize();

  return count;
}

unsigned long CartanClassSet::representative(realform::RealForm rf, size_t cn)
  const

/*!
  \brief Returns a representative for real form \#rf in Cartan \#cn.

  Precondition: cartan \#cn occurs for that real form.

  Explanation: this amounts to searching for rf in d_realFormLabels[cn].
*/

{
  using namespace partition;

  const realform::RealFormList& rfl = d_realFormLabels[cn];
  size_t p = std::find(rfl.begin(),rfl.end(),rf) - rfl.begin();
  const Partition& pi = cartan(cn).fiber().weakReal();

  return pi.classRep(p);
}

latticetypes::LatticeMatrix
CartanClassSet::involutionMatrix(const weyl::TwistedInvolution& tw)
  const
{
  weyl::WeylWord ww;
  weylGroup().out(ww,tw.representative());
  latticetypes::LatticeMatrix result;
  atlas::rootdata::toMatrix(result, ww,rootDatum());
  result *= distinguished();
  return result;
}

void CartanClassSet::twistedAct
  (const weyl::TwistedInvolution& tw,latticetypes::LatticeElt& v) const
{
  distinguished().apply(v,v);
  weylGroup().act(rootDatum(),tw.representative(),v);
}



/******** manipulators *******************************************************/

std::vector<weyl::WeylEltList> CartanClassSet::expand() const

/*!
  \brief Writes in wll the various conjugacy classes of twisted involutions
  for the currently known Cartans.

  Note: the lists come out sorted, to allow binary look-up; this is in fact
  dependent on the implementation of |weyl::WeylGroup::twistedConjugacyClass|
*/

{
  const weyl::WeylGroup& W = weylGroup();
  std::vector<weyl::WeylEltList> result(d_cartan.size());

  for (size_t j = 0; j < d_cartan.size(); ++j) {
    W.twistedConjugacyClass(result[j],d_twistedInvolution[j]);
  }

  return result;
}

/*!
  \brief Returns the index j s.t. w is a member of d_known[j] and j >= first;
         returns d_known.size() if there is no such j.
*/
size_t isMember(const std::vector<weyl::WeylEltList>& known,
		const weyl::WeylElt& w, size_t first)
{
  for (size_t j = first; j < known.size(); ++j)
    if (std::binary_search(known[j].begin(),known[j].end(),w))
      return j;

  return known.size();
}

void CartanClassSet::extend(realform::RealForm rf)

/*!
  \brief Extends the CartanClassSet structure so that it contains all Cartans
  for the real form |rf|.

  This is done by traversing all the known Cartan classes |j| defined for the
  real form |rf| (which includes at least the distinguished Cartan for the
  inner class), and all non-compact positive imaginary roots $\alpha$ for
  |(rf,j)|; the twisted involution |tw| for the Cartan class |j| is then
  left-multiplied by $s_\alpha$ to obtain a new twisted involution |ti|, and
  if |ti| is not member of any of the known classes of twisted involutions, it
  starts a new Cartan class defined in |rf|, which will later itself be
  considered in the loop described here.

  While generating the Cartan classes, the ordering is extended by a link from
  the Cartan class |j| to the Cartan class of any twisted involution |ti|
  obtained directly from it.

*/

{
  if (d_status.isMember(rf)) // nothing to be done
    return;

  size_t prev_Cartans=d_cartan.size(); // mark size for possible roll-back
  bitmap::BitMap prev_status=d_status; // save this one, is hard to roll back

  try
  {
    std::vector<weyl::WeylEltList> known=expand(); // get full twisted classes

    const rootdata::RootDatum dualRootDatum(rootDatum(),tags::DualTag());
    const weyl::WeylGroup dualWeylGroup(weylGroup(),tags::DualTag());

    std::set<poset::Link> lks;

    // d_cartan.size() may grow in the course of the loop
    for (size_t j = 0; j < d_cartan.size(); ++j) {

      if (not isDefined(rf,j))
	continue;

      rootdata::RootSet rs=noncompactPosRootSet(rf,j);

      for (rootdata::RootSet::iterator it = rs.begin(); it(); ++it)
      {

	// multipy |twistedInvolution(j)| on the left by the reflection |*it|
	TwistedInvolution tw=d_twistedInvolution[j]; leftReflect(tw,*it);

	// check if |tw| gives a new class; if so, extend the cartan structure
	size_t k=isMember(known,tw,j+1);// see if |tw| is in some class |k>j|
	lks.insert(std::make_pair(j,k)); // insert a link $j\to k$ in any case

	if (k == known.size()) { // found a new class
	  addCartan(*it,j); // add new |CartanClass| for |tw|
	  updateTwistedInvolutions(known,tw); // add |tw|, and extend |known|
// begin testing that |cartan(k)|
	  assert(cartan(k).orbitSize() == known[k].size());
// end testing
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
      d_support[i].resize(prev_Cartans);
      d_dualSupport[i].resize(prev_Cartans);
    }
    d_status.swap(prev_status); // restore set of completed real forms
    prev_status.andnot(d_status); // now |prev_status| flags "new" real forms
    for (bitmap::BitMap::iterator it=prev_status.begin(); it(); ++it)
      d_mostSplit[*it]=UndefMostSplit; // clear unsure most split Cartans
    throw; // rethrow same exception
  }
}





} // namespace cartan


/*****************************************************************************

        Chapter II --- private (auxiliary) methods for CartanClassSet

******************************************************************************/

namespace cartan {

void CartanClassSet::leftReflect(TwistedInvolution& tw, rootdata::RootNbr rn)
  const
/*!
  \brief Multiplies |tw| to the left by the reflection |s_rn| corresponding to
  root \#rn

  This is a twisted involution if that reflection twisted-commutes with the
  involution corresponding to |tw|; in practice root \#rn will be imaginary
  for |tw|
*/

{
  const weyl::WeylGroup& W = weylGroup();

  weyl::WeylWord rw=rootDatum().reflectionWord(rn);

  for (size_t i=rw.size(); i-->0; ) // left multiply |tw| by |rw|
    W.leftMult(tw,rw[i]);
}

rootdata::RootSet CartanClassSet::noncompactPosRootSet
  (realform::RealForm rf, size_t j) const
/*!
  \brief Flags in rs the set of noncompact positive roots for Cartan \#j.
*/

{
  rootdata::RootSet result;
  const cartanclass::Fiber& f = cartan(j).fiber();
  unsigned long x = CartanClassSet::representative(rf,j);
  f.noncompactRootSet(result,x);
  result &= rootDatum().posRootSet();
  return result;
}

/******** manipulators *******************************************************/

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
  WeylWord wi;
  W.involutionOut(wi,ti);

  // find cayley part and cross part
  RootList so;
  cayleyPart(so,wi,rd,W);
  WeylWord ww;
  crossPart(ww,wi,W);

  const Partition& pi = f.weakReal();
  realform::RealFormList rfl(f.numRealForms());

  // transform gradings and correlate forms
  for (size_t j = 0; j < rfl.size(); ++j) {
    unsigned long y = pi.classRep(j);
    Grading gr; f.grading(gr,y);
    RootList rl = f.simpleImaginary(); // the roots of which |gr| are a grading
    transformGrading(gr,rl,so,rd);
    for (size_t i = 0; i < so.size(); ++i)
      gr.set(rl.size()+i);             // make grading set for roots in |so|
    copy(so.begin(),so.end(),back_inserter(rl)); // extend |rl| with |so|
    crossTransform(rl,ww,rd);  // apply cross part of |ti| to roots in |rl|
    unsigned long x = makeRepresentative(gr,rl,fundf);
    realform::RealForm rf = fundf.weakReal()(x); // look up representative
    rfl[j] = rf;
  }

  d_realFormLabels.push_back(rfl);
}

/*!
  \brief Adds a new cartan, obtained from cartan \#j by Cayley transform
  through root \#rn.
*/
void CartanClassSet::addCartan(rootdata::RootNbr rn, size_t j)
{
  const rootdata::RootDatum& rd = rootDatum();
  latticetypes::LatticeMatrix q;
  rd.rootReflection(q,rn);
  q *= cartan(j).involution();
  d_cartan.push_back(new cartanclass::CartanClass(rd,q));
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

  // find dual twisted involution
  LatticeMatrix q = fundf.involution();
  q *= f.involution();
  WeylWord tiww;
  toWeylWord(tiww,q,rd);
  TwistedInvolution ti(WeylElt(tiww,W));

  WeylWord wi;
  W.involutionOut(wi,ti);

  // find cayley part and cross part
  RootList so;
  cayleyPart(so,wi,rd,W);
  WeylWord ww;
  crossPart(ww,wi,W);
// begin testing
  assert(checkDecomposition(wi,ww,so,W,rd));
// end testing

  const Partition& pi = f.weakReal();
  realform::RealFormList rfl(f.numRealForms());

  // transform gradings and correlate forms
  for (size_t j = 0; j < rfl.size(); ++j) {
    unsigned long y = pi.classRep(j);
    Grading gr;
    f.grading(gr,y);
    RootList rl = f.simpleImaginary();
    transformGrading(gr,rl,so,rd);
    for (size_t i = 0; i < so.size(); ++i)
      gr.set(rl.size()+i);
    copy(so.begin(),so.end(),back_inserter(rl));
    crossTransform(rl,ww,rd);
// begin testing
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
    d_support[j].resize(last+1);

  for (size_t j = 0; j < d_dualSupport.size(); ++j)
    d_dualSupport[j].resize(last+1);

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

} // namespace cartan


/*****************************************************************************

        Chapter III --- Functions declared in cartan.h

******************************************************************************/

namespace cartan {

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
  using namespace bitmap;
  using namespace cartanclass;

  BitMap b = ccl.support(rf);
  b &= ccl.dualSupport(drf);

  unsigned long c = 0;

  BitMap::iterator b_end = b.end();

  for (BitMap::iterator it = b.begin(); it(); ++it) {
    const CartanClass& cc = ccl.cartan(*it);
    unsigned long ci = cc.orbitSize();
    ci *= ccl.fiberSize(rf,*it);
    ci *= ccl.dualFiberSize(drf,*it);
    c += ci;
  }

  return c;
}

unsigned long kgbSize(realform::RealForm rf, const CartanClassSet& ccl)

/*!
  \brief Returns the cardinality of K\\G/B for this real form.

  Precondition: the Cartan classes for this real form have been generated

  Explanation: this is the cardinality of the one-sided parameter set
  corresponding to any strong real form over rf.
*/

{
  using namespace bitmap;
  using namespace cartanclass;

  BitMap b = ccl.support(rf);
  unsigned long c = 0;

  BitMap::iterator b_end = b.end();

  for (BitMap::iterator it = b.begin(); it(); ++it) {
    const CartanClass& cc = ccl.cartan(*it);
    unsigned long ci = cc.orbitSize();
    ci *= ccl.fiberSize(rf,*it);
    c += ci;
  }

  return c;
}

} // namespace cartan

/*****************************************************************************

        Chapter IV --- Auxiliary functions local to this module

******************************************************************************/

namespace cartan {

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
  \brief Returns an element |x| such that the elements in |rl| are graded
  by |x| according to |gr|.

  Precondition: the roots in |rl| are linearly independent, and |gr| is a
  grading of the root system with that basis.
*/

{
  using namespace bitset;
  using namespace cartanclass;
  using namespace latticetypes;
  using namespace gradings;
  using namespace rootdata;

  RootSet brs;
  fundf.noncompactRootSet(brs,0); // the base grading
  Grading bgr;
  restrictGrading(bgr,brs,rl);
  Component bc(bgr,rl.size());

  // make right hand side
  Component rhs(gr,rl.size());
  rhs += bc;

  // make grading shifts
  ComponentList cl(fundf.adjointFiberRank(),bc);
  for (size_t j = 0; j < cl.size(); ++j) {
    RootSet rs;
    fundf.noncompactRootSet(rs, 1 << j);
    Grading gr1;
    restrictGrading(gr1,rs,rl);
    cl[j] += Component(gr1,cl.size());
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
  \brief Transforms the grading |gr| of |rl| according to so.

  Precondition: |gr| is a grading of the roots in |rl|; |so| is a set of
  strongly orthogonal imaginary noncompact roots, orthogonal to the roots in
  |rl|.

  Explanation: the rule for Cayley-transforming a grading through an imaginary
  noncompact root $\alpha$ in |so| is that the grading of $\beta$ in |rl| is
  flipped if and only if $\alpha+\beta$ is a root.
*/

{
  for (size_t i = 0; i < rl.size(); ++i)
    for (size_t j = 0; j < so.size(); ++j)
      if (rootdata::sumIsRoot(rl[i],so[j],rd))
	gr.flip(i);
}

bool checkDecomposition(const weyl::WeylWord& wi,
			const weyl::WeylWord& ww,
			const rootdata::RootList& so,
			const weyl::WeylGroup& W,
			const rootdata::RootDatum& rd)

/*!
  \brief Checks whether wi decomposes as the composition of the
  cross-action defined by ww followed by the Cayley transform defined by so.
*/

{
  using namespace rootdata;
  using namespace weyl;

  TwistedInvolution tw;

  // cross action part
  for (size_t j = 0; j < ww.size(); ++j)
    W.twistedConjugate(tw,ww[j]);

  // cayley transform part
  for (size_t j = 0; j < so.size(); ++j) {
    WeylWord rww;
    toWeylWord(rww,so[j],rd);
    for (size_t i = rww.size(); i-->0;) {
      W.leftMult(tw,rww[i]);
    }
  }

  WeylWord wtest;
  W.involutionOut(wtest,tw);

  return wtest == wi;
}

} // namespace cartan
} // namespace atlas
