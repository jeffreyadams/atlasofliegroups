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

  The classifications of conjugacy classes of Cartan subgroups, for
  the various real forms, all embed into the classification of
  conjugacy classes of root data involutions for the given inner
  class, equivalent to the classification of twisted involutions in
  the Weyl group, which is a purely combinatorial Weyl group
  computation. This is not always a small computation, but can still
  be managed within a few seconds even for E8. The classification of
  real forms can be recovered in this picture as well, by looking at
  the unique most split Cartan class for each real form.

  The most delicate part is the "correlation" part: for each Cartan,
  tell which orbit in the corresponding fiber corresponds to which
  real form, the real forms being labelled by the orbits in the
  fundamental fiber.  In the current version, the solution to this
  problem is cleaner then previously, and perfectly general: it is
  obtained by writing out a system of equations for the grading
  defining the real form; this system does not always have a unique
  solution, but all solutions correspond to the same real form.

  To help with the construction, we introduce a helper class which in
  this case we prefer to be a derived class.
*/
/*
  This is cartan.cpp.
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "cartan.h"

#include <set>

#include "rootdata.h"
#include "tags.h"


namespace atlas {

  namespace cartan{
  namespace helper{

  bool checkDecomposition(const weyl::WeylWord&, const weyl::WeylWord&,
			  const rootdata::RootList&, const weyl::WeylGroup&,
			  const rootdata::RootDatum&);

  const size_t UndefMostSplit = ~0ul;
  const realform::RealForm UndefRealForm = ~0ul;

      /*!
\brief Derived class of CartanClassSet, to carry out the construction
of CartanClassSet.
      */
  class Helper:public CartanClassSet {

  private:

    /* during construction these dual objects will be needed */
    const rootdata::RootDatum d_dualRootDatum;
    const weyl::WeylGroup d_dualWeylGroup;

    /*!
    \brief List of twisted conjugacy classes of twisted
    involutions, one for each non-fundamental Cartan.

    The conjugacy classes can be very large, so this constructor can
    be quite slow.  For a complex group, the twisted conjugacy class
    of the identity has cardinality the order of the Weyl group; so
    for complex An, this is n!.  But the twisted conjugacy class is
    never computed for the fundamental Cartan, so this does not slow
    the program for complex groups.  For the product of a complex
    group and (say) SL(2,R), the full class is computed on the second
    Cartan subgroup; so computing the Cartans for SL(2,R) x (complex
    E7) takes several minutes (while the full Weyl group of E7,
    consisting of 2903040 elements, is traversed).
    */
    std::vector<weyl::WeylEltList> d_known;

  public:
// constructors and destructors

    /* the initial constructor */
    Helper(const complexredgp::ComplexReductiveGroup&,
	   const latticetypes::LatticeMatrix&);

    /* the following constructor is needed when adding new Cartans */
    Helper(const CartanClassSet&); // resurrect Helper object for existing base

    ~Helper();

// accessors
    void cayleyPart(rootdata::RootList&, const weyl::WeylWord&,
		    const rootdata::RootDatum&, const weyl::WeylGroup&) const;

    void crossPart(weyl::WeylWord&, const weyl::WeylWord&,
		   const weyl::WeylGroup&) const;

    void crossTransform(rootdata::RootList&, const weyl::WeylWord&,
			const rootdata::RootDatum&) const;

    const rootdata::RootDatum& dualRootDatum() const {
      return d_dualRootDatum;
    }

    const weyl::WeylGroup& dualWeylGroup() const {
      return d_dualWeylGroup;
    }

    size_t isMember(const weyl::WeylElt&, size_t) const;

    /*! \brief Whether real form #c has Cartan #j as most split Cartan */
    bool isMostSplit(unsigned long c, size_t j) const {
      return cartan(j).isMostSplit(c);
    }

    unsigned long makeRepresentative(const gradings::Grading&,
				     const rootdata::RootList&,
				     const cartanclass::Fiber&) const;

    void makeTwistedInvolution(weyl::TwistedInvolution&,
			       const weyl::TwistedInvolution&,
			       rootdata::RootNbr) const;

    void noncompactPosRootSet(rootdata::RootSet&, realform::RealForm, size_t)
      const;

    void transformGrading(gradings::Grading&,
			  const rootdata::RootList&,
			  const rootdata::RootList&,
			  const rootdata::RootDatum&) const;

// manipulators
    void addCartan(rootdata::RootNbr, size_t);

    void correlateDualForms();

    void correlateForms();

    void expand();

    void extend(realform::RealForm);

    void updateStatus(size_t);

    void updateSupports();

    void updateTwistedInvolutions(const TwistedInvolution&);

  };
}
}

/*****************************************************************************

        Chapter I --- The CartanClassSet class

******************************************************************************/

namespace cartan {

  /* Initial construction.
     Install parent and let Helper construction do the rest
  */

CartanClassSet::CartanClassSet(const complexredgp::ComplexReductiveGroup& p,
			       const latticetypes::LatticeMatrix& d)
  :d_parent(p)
{
  helper::Helper h(p,d);
  swap(h);
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
CartanClassSet::CartanClassSet(const CartanClassSet& source)
  :d_parent(source.d_parent),
   d_fundamental(source.d_fundamental),
   d_dualFundamental(source.d_dualFundamental),
   d_cartan(source.d_cartan),
   d_twistedInvolution(source.d_twistedInvolution),
   d_ordering(source.d_ordering),
   d_realFormLabels(source.d_realFormLabels),
   d_dualRealFormLabels(source.d_dualRealFormLabels),
   d_support(source.d_support),
   d_dualSupport(source.d_dualSupport),
   d_status(source.d_status),
   d_mostSplit(source.d_mostSplit)

{
  // get own copy of each cartan class, since we own them
  for (size_t j = 0; j < d_cartan.size(); ++j)
    if (d_cartan[j]!=NULL)
      d_cartan[j] = new cartanclass::CartanClass(*d_cartan[j]);
}

void CartanClassSet::swap (CartanClassSet& other)
  // Note that parent is not, and cannot be, swapped: they must agree
{
  d_fundamental.swap(other.d_fundamental);
  d_dualFundamental.swap(other.d_dualFundamental);
  d_cartan.swap(other.d_cartan);
  d_twistedInvolution.swap(other.d_twistedInvolution);
  d_ordering.swap(other.d_ordering);
  d_realFormLabels.swap(other.d_realFormLabels);
  d_dualRealFormLabels.swap(other.d_dualRealFormLabels);
  d_support.swap(other.d_support);
  d_dualSupport.swap(other.d_dualSupport);
  d_status.swap(other.d_status);
  d_mostSplit.swap(other.d_mostSplit);
}

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

void CartanClassSet::extend(const weyl::WeylGroup& W,
			   const rootdata::RootDatum& rd,
			   realform::RealForm rf)

/*!
  \brief Extends the CartanClassSet structure so that it contains all Cartans
  for real form \#rf.

  Precondition: W is the Weyl group of the corresponding complex reductive
  group; rd is the root datum;
*/

{
  if (d_status.isMember(rf)) // nothing to be done
    return;

  helper::Helper h(*this);

  h.extend(rf);

  // commit
  swap(h);
}





} // namespace cartan


/*****************************************************************************

        Chapter II --- The Helper class for CartanClassSet

******************************************************************************/

namespace cartan {
  namespace helper {

    /* Initial construction of Helper for CartanClassSet */

Helper::Helper(const complexredgp::ComplexReductiveGroup& parent,
	       const latticetypes::LatticeMatrix& d)
  : CartanClassSet(parent)
  , d_dualRootDatum(parent.rootDatum(),tags::DualTag())
  , d_dualWeylGroup(parent.weylGroup(),tags::DualTag())

{
  // construct fundamental Cartan, serving as seed for all other Cartans
  d_cartan.push_back(new cartanclass::CartanClass(parent.rootDatum(),d));

  // create a singleton poset representing the fundamental Cartan
  d_ordering = poset::Poset(1);

  // the fundamental fiber is (a copy of) the fiber of the fundamental Cartan
  d_fundamental = d_cartan[0]->fiber();

  /* but the dual fundamental fiber is NOT the dual fiber of the fundamental
     Cartan; it is the dual fiber of the most split Cartan that will not be
     constructed here (if constructed at all, it will be the very last Cartan
     for this inner class). Therefore this fiber is constructed separately */

  { // get distinguished involution for the dual group
    latticetypes::LatticeMatrix dd;
    dualBasedInvolution(dd,d,parent.rootDatum());

    // construct the dual fundamental fiber
    d_dualFundamental = cartanclass::Fiber(dualRootDatum(),dd);
  }

  // add corresponding twisted involution
  d_twistedInvolution.push_back(weyl::WeylElt());

  // make the real form labels
  correlateForms();

  // make the dual real form labels
  correlateDualForms();

  // mark off the fundamental forms
  d_status.resize(numRealForms());
  d_mostSplit.assign(numRealForms(),UndefMostSplit);
  updateStatus(0);

  // initialize supports
  d_support.resize(numRealForms());
  d_dualSupport.resize(numDualRealForms());
  updateSupports();
}


/*
 Constructor used when extending a CartanClassSet (via CartanClassSet::extend)
*/
Helper::Helper(const CartanClassSet& source)
  : CartanClassSet(source)
  , d_dualRootDatum(source.rootDatum(),tags::DualTag())
  , d_dualWeylGroup(source.weylGroup(),tags::DualTag())
{}

Helper::~Helper()

{}

/******** copy, assignment and swap ******************************************/

/******** accessors **********************************************************/

void Helper::cayleyPart(rootdata::RootList& so,
			const weyl::WeylWord& wi,
			const rootdata::RootDatum& rd,
			const weyl::WeylGroup& W) const

/*!
  \brief Puts in so the composite Cayley transform corresponding to wi.

  Precondition: |wi| is a reduced expression for a twisted involution |tw|;

  Explanation: to each root datum involution, we may associate a transformation
  from the fundamental cartan to the current cartan, that factors as the
  composition of a cross action and a composite Cayley transform. This function
  puts in |so| a system of strongly orthogonal roots representing the cayley
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

void Helper::crossPart(weyl::WeylWord& ww, const weyl::WeylWord& wi,
		       const weyl::WeylGroup& W) const

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

void Helper::crossTransform(rootdata::RootList& rl,
			    const weyl::WeylWord& ww,
			    const rootdata::RootDatum& rd) const

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

void Helper::makeTwistedInvolution(TwistedInvolution& ti,
				   const TwistedInvolution& tw,
				   rootdata::RootNbr rn) const
/*!
  \brief Puts in |ti| the product of the reflection |s_rn| corresponding to
  root \#rn, and |tw|.

  This is a twisted involution if that reflection commutes with the involution
  corresponding to |tw|; in practice root \#rn will be imaginary for |tw|
*/

{
  using namespace weyl;

  const WeylGroup& W = weylGroup();

  WeylWord rw;
  toWeylWord(rw,rn,rootDatum()); // get reflection as Weyl word

  ti=tw;
  for (size_t i=rw.size(); i-->0; ) // left multiply |ti| by |rw|
    W.leftMult(ti,rw[i]);
}

unsigned long Helper::makeRepresentative(const gradings::Grading& gr,
					 const rootdata::RootList& rl,
					 const cartanclass::Fiber& fundf) const
/*!
  \brief Returns an element x such that the elements in rl are graded
  by x according to gr.

  Precondition: the roots in rl are linearly independent.
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

void Helper::noncompactPosRootSet(rootdata::RootSet& rs, realform::RealForm rf,
				  size_t j) const

/*!
  \brief Flags in rs the set of noncompact positive roots for Cartan \#j.
*/

{
  const cartanclass::Fiber& f = cartan(j).fiber();
  unsigned long x = CartanClassSet::representative(rf,j);
  f.noncompactRootSet(rs,x);
  rs &= rootDatum().posRootSet();
}

void Helper::transformGrading(gradings::Grading& gr,
			      const rootdata::RootList& rl,
			      const rootdata::RootList& so,
			      const rootdata::RootDatum& rd) const

/*!
  \brief Transforms the grading gr of rl according to so.

  Precondition: gr is a grading of the roots in rl; so is a set of strongly
  orthogonal roots, orthogonal to the roots in rl.

  Explanation: the rule for Cayley-tranforming a grading through an
  imaginary noncompact root alpha is that the grading of beta is flipped
  iff alpha+beta is a root.
*/

{
  using namespace rootdata;

  for (size_t i = 0; i < rl.size(); ++i)
    for (size_t j = 0; j < so.size(); ++j)
      if (sumIsRoot(rl[i],so[j],rd))
	gr.flip(i);

  return;
}

/******** manipulators *******************************************************/

void Helper::addCartan(rootdata::RootNbr rn, size_t j)

/*!
  \brief Adds a new cartan, obtained from cartan \#j by Cayley transform
  through root \#rn.
*/

{
  const rootdata::RootDatum& rd = rootDatum();
  latticetypes::LatticeMatrix q;
  rd.rootReflection(q,rn);
  q *= cartan(j).involution();
  d_cartan.push_back(new cartanclass::CartanClass(rd,q));
}

void Helper::correlateDualForms()

{
  using namespace cartanclass;
  using namespace gradings;
  using namespace latticetypes;
  using namespace partition;
  using namespace rootdata;
  using namespace weyl;

  const RootDatum& rd = dualRootDatum();
  const WeylGroup& W = dualWeylGroup();

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

  return;
}

void Helper::correlateForms()

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
    Grading gr;
    f.grading(gr,y);
    RootList rl = f.simpleImaginary();
    transformGrading(gr,rl,so,rd);
    for (size_t i = 0; i < so.size(); ++i)
      gr.set(rl.size()+i);
    copy(so.begin(),so.end(),back_inserter(rl));
    crossTransform(rl,ww,rd);
    unsigned long x = makeRepresentative(gr,rl,fundf);
    realform::RealForm rf = fundf.weakReal()(x);
    rfl[j] = rf;
  }

  d_realFormLabels.push_back(rfl);

  return;
}

void Helper::expand()

/*!
  \brief Writes in wll the various conjugacy classes of twisted involutions
  for the currently known Cartans.

  NOTE: this may forward an Overflow error. Commit-or-rollback is guaranteed.
*/

{
  using namespace weyl;

  const weyl::WeylGroup& W = weylGroup();
  std::vector<weyl::WeylEltList> wll;

  for (size_t j = 0; j < d_cartan.size(); ++j) {
    wll.push_back(WeylEltList());
    W.conjugacyClass(wll.back(),d_twistedInvolution[j],true); // twisted class
  }

  // commit
  d_known.swap(wll);

  return;
}

void Helper::extend(realform::RealForm rf)

/*!
  \brief Extends the CartanClassSet structure so that it contains all Cartans
  for rf.

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
  using namespace bitmap;
  using namespace cartanclass;
  using namespace latticetypes;
  using namespace poset;
  using namespace rootdata;
  using namespace weyl;

  size_t prevSize = d_cartan.size();

  // expand the known classes
  expand(); // this may throw!

  std::set<poset::Link> lks;

  // d_cartan.size() may grow in the course of the loop
  for (size_t j = 0; j < d_cartan.size(); ++j) {

    if (not isDefined(rf,j))
      continue;

    rootdata::RootSet rs;
    noncompactPosRootSet(rs,rf,j);

    TwistedInvolution tw = d_twistedInvolution[j];

    for (rootdata::RootSet::iterator it = rs.begin(); it(); ++it) {

      // multipy |tw| by the reflection |*it| on the left
      TwistedInvolution ti;
      makeTwistedInvolution(ti,tw,*it);

      // check if we have a new twisted involution
      // if yes, extend the cartan structure
      size_t k = isMember(ti,j+1); // see if |ti| is in some class |k>j|
      lks.insert(std::make_pair(j,k)); // insert a link $j\to k$ in any case

      if (k == d_known.size()) { // found a new class
	addCartan(*it,j); // add new |CartanClass| for this twisted involution
	updateTwistedInvolutions(ti);
// begin testing
assert(d_cartan.back()->orbitSize() == d_known.back().size());
// end testing
	correlateForms();
	correlateDualForms();
	updateSupports();
      }
    }
  }

  // update the status flags
  updateStatus(prevSize);

  // update the order relation
  std::vector<poset::Link> lk(lks.begin(),lks.end());
  d_ordering.resize(d_cartan.size());
  d_ordering.extend(lk);
}

size_t Helper::isMember(const weyl::WeylElt& w, size_t first) const

/*!
  \brief Returns the index j s.t. w is a member of d_known[j] and j >=
  first; returns d_known.size() if there is no such j.
*/

{
  for (size_t j = first; j < d_known.size(); ++j)
    if (std::binary_search(d_known[j].begin(),d_known[j].end(),w))
      return j;

  return d_known.size();
}

void Helper::updateStatus(size_t prev)

/*!
  \brief Updates d_status.

  Precondition: prev is the previous size of d_cartan;

  Explanation: d_status holds the subset of the set of real forms for
  which the full set of Cartan classes is constructed; equivalently, those
  for which the most split Cartan has been reached.
*/

{
  using namespace cartanclass;
  using namespace partition;

  for (size_t j = prev; j < d_cartan.size(); ++j) {

    const Fiber& f = cartan(j).fiber();
    const Partition& pi = f.weakReal();
    const realform::RealFormList& rfl = realFormLabels(j);

    for (unsigned long c = 0; c < pi.classCount(); ++c) {
      if (isMostSplit(c,j)) { // flag the form number of x
	d_status.insert(rfl[c]);
	d_mostSplit[rfl[c]] = j;
      }
    }

  }

  return;
}

void Helper::updateSupports()

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

  size_t m = d_cartan.size();

  for (size_t j = 0; j < d_support.size(); ++j)
    d_support[j].resize(m);

  for (size_t j = 0; j < d_dualSupport.size(); ++j)
    d_dualSupport[j].resize(m);

  const realform::RealFormList& rfl = d_realFormLabels.back();

  for (size_t j = 0; j < rfl.size(); ++j) {
    d_support[rfl[j]].insert(m-1);
  }

  const realform::RealFormList& drfl = d_dualRealFormLabels.back();

  for (size_t j = 0; j < drfl.size(); ++j) {
    d_dualSupport[drfl[j]].insert(m-1);
  }

  return;
}

void Helper::updateTwistedInvolutions(const TwistedInvolution& tw)

/*!
  \brief Updates the known list by adding the canonicalRep of tw to it.
*/

{
  weyl::WeylEltList wl;
  weylGroup().twistedConjugacyClass(wl,tw);

  d_twistedInvolution.push_back(tw);

  d_known.push_back(weyl::WeylEltList());
  d_known.back().swap(wl);
}

} // namespace helper
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
  namespace helper {

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

} // namespace helper
} // namespace cartan
} // namespace atlas
