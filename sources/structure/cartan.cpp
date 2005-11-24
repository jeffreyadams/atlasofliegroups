/*
  This is cartan.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "cartan.h"

#include <set>

#include "rootdata.h"
#include "tags.h"

/*****************************************************************************

  This module, together with the cartanclass one, contains code for dealing
  with conjugacy classes of Cartan subgroups, and real forms. In this version,
  we have de-emphasized the classification of strong real forms, which is
  perhaps better treated as an input-output issue (it can certainly be 
  recovered with a moderate effort from the data available here.) It appears
  that the data recorded here is what's most relevant to representation
  theory.

  The classification of real forms amounts to the classification of W^\delta-
  orbits in the fundamental fiber of the one-sided parameter space for the
  adjoint group (see the "combinatorics" paper on the Atlas website, or the 
  forthcoming "algorithms" paper by Jeff Adams and me.); this is a very small
  computation, depending only on the Lie algebra.

  The classifications of conjugacy classes of Cartan subgroups, for the
  various real forms, all embed into the classification of conjugacy classes
  of root data involutions for the given inner class, equivalent to the
  classification of twisted involutions in the Weyl group, which is a purely
  combinatorial Weyl group computation. This is not always a small computation,
  but can still be managed within a few seconds even for E8. The classification
  of real forms can be recovered in this picture as well, by looking at the
  unique most split Cartan class for each real form.

  The most delicate part is the "correlation" part: for each Cartan, tell
  which orbit in the corresponding fiber corresponds to which real form,
  the real forms being labelled by the orbits in the fundamental fiber.
  In the current version, the solution to this problem is cleaner then
  previously, and perfectly general: it is obtained by writing out a system
  of equations for the grading defining the real form; this system does not
  always have a unique solution, but all solutions correspond to the same
  real form.

  To help with the construction, we introduce a helper class which in this
  case we prefer to be a derived class.

******************************************************************************/

namespace atlas {

namespace {

  using namespace cartan;

  bool checkDecomposition(const weyl::WeylWord&, const weyl::WeylWord&, 
			  const rootdata::RootList&, const weyl::WeylGroup&, 
			  const rootdata::RootDatum&);

  const size_t UndefMostSplit = ~0ul;
  const realform::RealForm UndefRealForm = ~0ul;

  class Helper:public CartanClasses {

  private:

    const rootdata::RootDatum* d_rootDatum;
    const rootdata::RootDatum d_dualRootDatum;
    const weyl::WeylGroup* d_weylGroup;
    const weyl::WeylGroup d_dualWeylGroup;

    std::vector<weyl::WeylEltList> d_known;

  public:
// constructors and destructors
    Helper(const CartanClasses&, const rootdata::RootDatum&, 
	   const weyl::WeylGroup&);

    Helper(const rootdata::RootDatum&, const latticetypes::LatticeMatrix&,
	   const weyl::WeylGroup&);

    virtual ~Helper();

// copy, assignment and swap
    void swap(CartanClasses&);

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

    bool isMostSplit(unsigned long c, size_t j) const {
      return cartan(j).isMostSplit(c);
    }

    unsigned long makeRepresentative(const gradings::Grading&, 
				     const rootdata::RootList&,
				     const cartanclass::Fiber&) const;

    void makeTwistedInvolution(weyl::WeylElt&, const weyl::WeylElt&,
			       rootdata::RootNbr) const;

    void noncompactPosRootSet(rootdata::RootSet&, realform::RealForm, size_t) 
      const;

    const rootdata::RootDatum& rootDatum() const {
      return *d_rootDatum;
    }

    void transformGrading(gradings::Grading&, 
			  const rootdata::RootList&, 
			  const rootdata::RootList&, 
			  const rootdata::RootDatum&) const;

    const weyl::WeylGroup& weylGroup() const {
      return *d_weylGroup;
    }

// manipulators
    void addCartan(rootdata::RootNbr, size_t);

    void correlateDualForms();

    void correlateForms();

    void expand();    
    
    void extend(realform::RealForm);
    
    void updateStatus(size_t);

    void updateSupports();

    void updateTwistedInvolutions(const weyl::WeylElt&);

  };
}

/*****************************************************************************

        Chapter I --- The CartanClasses class

  The CartanClasses structure holds the main structural data required for
  the representation theory of our inner class of real forms. The meaning
  of the various fields is as follows:

    - d_fundamental: the fundamental fiber; is also d_cartan[0].fiber().
      Real forms are classified from the weakReal orbit partition afforded
      by this fiber.

    - d_dualFundamental: same on the dual side. This is _not_ d_cartan[0].
      dualFiber()! Rather, it is the dual fiber of the most split Cartan
      of the quasisplit form.

    - d_cartan: a list of pointers to CartanClass structures, one for each
      Cartan currently constructed. Initially, we only construct the
      fundamental one; the structure is grown as required.

    - d_ordering: the order structure on conjugacy classes of Cartans. This
      is a poset of size d_cartan.size(); it is grown as the structure is
      expanded. It is a decreasing subset of the full ordering on root datum
      involutions, which is reached for the quasisplit form.

    - d_realFormLabels: a list of lists, one list for each Cartan class.
      It contains an external labelling of real forms (the internal labelling
      being in terms of the orbits in the adjoint fiber.) This is delicate
      to construct, because it involves correlating the orbit pictures for
      the various Cartans, which requires Cayley transforms.

    - d_dualRealForms: the analogous thing for the dual group.

    - d_support: for each real form, tells which of the currently constructed
      Cartans are defined for that real form; in other words, it is the set
      of Cartans for which that real form label occurs in the corresponding
      realFormLabels list.

    - d_dualSupport: same for the dual group.

    - d_status: a bitmap of size numRealForms. Tells for which real forms
      the Cartan subgroup structure has been fully constructed.

******************************************************************************/

namespace cartan {

CartanClasses::CartanClasses(const rootdata::RootDatum& rd, 
			     const latticetypes::LatticeMatrix& d,
			     const weyl::WeylGroup& W)

{
  Helper h(rd,d,W);
  swap(h);
}

CartanClasses::~CartanClasses()

/*
  The only data explicitly allocated by the CartanClass is that for the
  CartanClass pointers.
*/

{
  for (size_t j = 0; j < d_cartan.size(); ++j)
    delete d_cartan[j];
}

/********* copy, assignment and swap *****************************************/
CartanClasses::CartanClasses(const CartanClasses& source)
  :d_fundamental(source.d_fundamental),
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
  using namespace cartanclass;

  // get own copy of cartan classes
  for (size_t j = 0; j < d_cartan.size(); ++j)
    if (d_cartan[j])
      d_cartan[j] = new CartanClass(*d_cartan[j]);
}

CartanClasses& CartanClasses::operator= (const CartanClasses& source)

{
  // handle self-assignment
  if (&source != this) {
    this->~CartanClasses();
    new(this) CartanClasses(source);
  }

  return *this;
}

void CartanClasses::swap (CartanClasses& other)

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

  return;
}

/******** accessors **********************************************************/

unsigned long CartanClasses::dualFiberSize(realform::RealForm rf, size_t cn) 
  const

/*
  Synopsis: returns the size of the fiber orbits corresponding to the strong
  real forms lying over rf, in cartan #cn.

  Precondition: rf is defined for cartan #cn.
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

unsigned long CartanClasses::dualRepresentative(realform::RealForm rf, 
						size_t cn) const

/*
  Synopsis: returns a representative for dual real form #rf in Cartan #cn.

  Precondition: cartan #cn occurs for that dual real form.

  Explanation: this amounts to searching for rf in d_dualRealFormLabels[cn].
*/

{
  using namespace partition;

  const realform::RealFormList& rfl = d_dualRealFormLabels[cn];
  size_t p = std::find(rfl.begin(),rfl.end(),rf) - rfl.begin();
  const Partition& pi = cartan(cn).dualFiber().weakReal();

  return pi.classRep(p);
}

unsigned long CartanClasses::fiberSize(realform::RealForm rf, size_t cn) const

/*
  Synopsis: returns the size of the fiber orbits corresponding to the strong
  real forms lying over rf, in cartan #cn.

  Precondition: rf is defined for cartan #cn.
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

unsigned long CartanClasses::representative(realform::RealForm rf, size_t cn) 
  const

/*
  Synopsis: returns a representative for real form #rf in Cartan #cn.

  Precondition: cartan #cn occurs for that real form.

  Explanation: this amounts to searching for rf in d_realFormLabels[cn].
*/

{
  using namespace partition;

  const realform::RealFormList& rfl = d_realFormLabels[cn];
  size_t p = std::find(rfl.begin(),rfl.end(),rf) - rfl.begin();
  const Partition& pi = cartan(cn).fiber().weakReal();

  return pi.classRep(p);
}

/******** manipulators *******************************************************/

void CartanClasses::extend(const weyl::WeylGroup& W, 
			   const rootdata::RootDatum& rd,
			   realform::RealForm rf)

/*
  Synopsis: extends the cartanclasses structure so that it contains all Cartans
  for rf.

  Precondition: W is the Weyl group of the corresponding complex reductive 
  group; rd is the root datum;

  NOTE: may forward an Overflow exception from expand.
*/

{
  using namespace bitmap;
  using namespace cartanclass;
  using namespace latticetypes;
  using namespace poset;
  using namespace rootdata;
  using namespace tags;
  using namespace weyl;

  if (d_status.isMember(rf)) // nothing to be done
    return;

  Helper h(*this,rd,W);
  h.extend(rf);

  // commit
  swap(h);

  return;
}

}

/*****************************************************************************

        Chapter II --- The Helper class

  ... explain here when it is stable ...

******************************************************************************/

namespace {

Helper::Helper(const rootdata::RootDatum& rd, 
	       const latticetypes::LatticeMatrix& d,
	       const weyl::WeylGroup& W)
  :d_rootDatum(&rd),
   d_dualRootDatum(rd,tags::DualTag()),
   d_weylGroup(&W),
   d_dualWeylGroup(W,tags::DualTag())

{  
  using namespace bitset;
  using namespace cartanclass;
  using namespace latticetypes;
  using namespace poset;
  using namespace rootdata;
  using namespace weyl;

  d_fundamental = Fiber(rd,d);
  d_ordering = Poset(1);

  LatticeMatrix dd;
  dualBasedInvolution(dd,d,rd);
  d_dualFundamental = Fiber(dualRootDatum(),dd);

  // construct fundamental Cartan
  d_cartan.push_back(new CartanClass(rd,d));

  // add corresponding twisted involution
  d_twistedInvolution.push_back(WeylElt());

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

Helper::Helper(const CartanClasses& source, const rootdata::RootDatum& rd,
	       const weyl::WeylGroup& W)
  :CartanClasses(source),
   d_rootDatum(&rd),
   d_dualRootDatum(rd,tags::DualTag()),
   d_weylGroup(&W),
   d_dualWeylGroup(W,tags::DualTag())

{}

Helper::~Helper() 

{}

/******** copy, assignment and swap ******************************************/

/******** accessors **********************************************************/
void Helper::cayleyPart(rootdata::RootList& so, 
			const weyl::WeylWord& wi,
			const rootdata::RootDatum& rd, 
			const weyl::WeylGroup& W) const

/*
  Synopsis: puts in so the composite Cayley transform corresponding to wi.

  Precondition: wi is a reduced expression for a twisted involution;

  Explanation: to each root datum involution, we may associate a transformation
  from the fundamental cartan to the current cartan, that factors as the
  composition of a cross action and a composite Cayley transform. This function
  puts in so a system of strongly orthogonal roots representing the cayley
  transform.
*/

{
  using namespace rootdata;
  using namespace weyl;

  WeylElt w;
  so.clear();

  for (size_t j = 0; j < wi.size(); ++j) {
    Generator s = wi[j];
    if (W.hasTwistedCommutation(s,w)) { // add new root to so
      so.push_back(rd.simpleRootNbr(s));
      W.leftProd(w,s);
    }
    else { // conjugate roots in so
      for (size_t i = 0; i < so.size(); ++i)
	so[i] = rd.rootPermutation(s)[so[i]];
      W.twistedConjugate(w,s);
    }
  }
    
  strongOrthogonalize(so,rd);

  return;
}

void Helper::crossPart(weyl::WeylWord& ww, const weyl::WeylWord& wi,
		       const weyl::WeylGroup& W) const

/*
  Synopsis: puts in ww the cross action corresponding to ti.

  Explanation: to each root datum involution, we may associate a transformation
  from the fundamental cartan to the current cartan, that factors as the
  composition of a cross action and a composite Cayley transform. This function
  puts in ww a weyl word representing the cross action part.
*/

{
  using namespace weyl;

  WeylElt w;
  ww.clear();

  for (size_t j = 0; j < wi.size(); ++j) {
    Generator s = wi[j];
    if (W.hasTwistedCommutation(s,w)) { // cayley transform
      W.leftProd(w,s);
    } else { // cross action
      ww.push_back(s);
      W.twistedConjugate(w,s);
    }
  }

  return;
}

void Helper::crossTransform(rootdata::RootList& rl, 
			    const weyl::WeylWord& ww,
			    const rootdata::RootDatum& rd) const

/*
  Synopsis: cross-transforms the roots in rl according to ww.

  NOTE: the cross-transformations are done in _reverse_ order.
*/

{
  using namespace rootdata;
  using namespace setutils;

  for (size_t j = ww.size(); j;) {
    --j;
    Permutation a = rd.rootPermutation(ww[j]);
    for (size_t i = 0; i < rl.size(); ++i)
      rl[i] = a[rl[i]];
  }

  return;
}

void Helper::makeTwistedInvolution(weyl::WeylElt& ti, const weyl::WeylElt& w,
				   rootdata::RootNbr rn) const

/*
  Synopsis: puts in ti the product of the reflection s_rn corresponding to
  root #rn, with w.

  Precondition: ti is set to zero (otherwise, it is multiplied on the right
  with the twisted involution constructed here);

  This is a twisted involution if w is, and rn is imaginary w.r.t. the
  involution corresponding to w and the current twist.
*/


{
  using namespace weyl;

  const WeylGroup& W = weylGroup();

  ti.clear();

  WeylWord rw;
  toWeylWord(rw,rn,rootDatum());
  W.prod(ti,rw);
  W.prod(ti,w);
      
  return;
}

unsigned long Helper::makeRepresentative(const gradings::Grading& gr, 
					 const rootdata::RootList& rl,
					 const cartanclass::Fiber& fundf) const

/*
  Synopsis: returns an element x such that the elements in rl are graded
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

/*
  Synopsis: flags in rs the set of noncompact positive roots for Cartan #j.
*/

{    
  using namespace cartanclass;

  const Fiber& f = cartan(j).fiber();
  unsigned long x = representative(rf,j);
  f.noncompactRootSet(rs,x);
  rs &= rootDatum().posRootSet();

  return;
}

void Helper::transformGrading(gradings::Grading& gr, 
			      const rootdata::RootList& rl,
			      const rootdata::RootList& so,
			      const rootdata::RootDatum& rd) const

/*
  Synopsis: transforms the grading gr of rl according to so.

  Precondition: gr is a grading of the roots in rl; so is a set of strongly
  orthogonal roots, orthogonal to the roots in rl.

  Expplanation: the rule for Cayley-tranforming a grading through an
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

/*
  Synopsis: adds a new cartan, obtained from cartan #j by Cayley transform
  through root #j.
*/

{	
  using namespace cartanclass;
  using namespace latticetypes;
  using namespace rootdata;

  const RootDatum& rd = rootDatum();

  LatticeMatrix q;
  rd.rootReflection(q,rn);
  q *= cartan(j).involution();
  d_cartan.push_back(new CartanClass(rd,q));

  return;
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
  WeylElt ti;
  W.prod(ti,tiww);

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

/*
  Synopsis: adds a new real form label list to d_realFormLabels

  Algorithm: the gradings of the imaginary root system corresponding to rhe
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

  const WeylElt& ti = d_twistedInvolution.back();
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

/*
  Synopsis: writes in wll the various conjugacy classes of twisted involutions
  for the currently known Cartans.

  NOTE: this may forward an Overflow error. Commit-or-rollback is guaranteed.
*/

{
  using namespace weyl;

  const weyl::WeylGroup& W = weylGroup();
  std::vector<weyl::WeylEltList> wll;

  for (size_t j = 0; j < d_cartan.size(); ++j) {
    wll.push_back(WeylEltList());
    W.conjugacyClass(wll.back(),d_twistedInvolution[j]);
  }

  // commit
  d_known.swap(wll);

  return;
}

void Helper::extend(realform::RealForm rf)

/*
  Synopsis: extends the cartanclasses structure so that it contains all Cartans
  for rf.

  NOTE: may forward an Overflow exception from expand.
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

  std::set<Link> lks;

  // ccl.d_cartan.size() may grow in the course of the loop
  for (size_t j = 0; j < d_cartan.size(); ++j) {

    if (not isDefined(rf,j))
      continue;

    RootSet rs;
    noncompactPosRootSet(rs,rf,j);

    WeylElt w = d_twistedInvolution[j];

    for (RootSet::iterator i = rs.begin(); i != rs.end(); ++i) {

      // multipy w with the reflection *i on the left
      WeylElt ti;
      makeTwistedInvolution(ti,w,*i);

      // check if we have a new twisted involution
      // if yes, extend the cartan structure
      size_t k = isMember(ti,j+1);
      lks.insert(std::make_pair(j,k));

      if (k == d_known.size()) { // found a new class
	addCartan(*i,j);
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
  std::vector<Link> lk(lks.begin(),lks.end());
  d_ordering.resize(d_cartan.size());
  d_ordering.extend(lk);

  return;
}

size_t Helper::isMember(const weyl::WeylElt& w, size_t first) const

/*
  Synopsis: returns the index j s.t. w is a member of d_known[j] and j >= 
  first; returns d_known.size() if there is no such j.
*/

{
  for (size_t j = first; j < d_known.size(); ++j)
    if (std::binary_search(d_known[j].begin(),d_known[j].end(),w))
      return j;

  return d_known.size();
}

void Helper::updateStatus(size_t prev)

/*
  Synopsis: updates d_status.

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

/*
  Synopsis: updates d_support and d_dualSupport.

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

void Helper::updateTwistedInvolutions(const weyl::WeylElt& w)

/*
  Synopsis: updates the known list by adding the class of w to it.

  NOTE: may forward an Overflow exception. We provide a commit-or-rollback
  guarantee.
*/

{
  using namespace weyl;

  weyl::WeylEltList wl;
  weylGroup().conjugacyClass(wl,w);

  // if we get here, no exception was thrown

  d_twistedInvolution.push_back(w);

  d_known.push_back(WeylEltList());
  d_known.back().swap(wl);

  return;
}

}

/*****************************************************************************

        Chapter III --- Functions declared in cartan.h

  ... explain here when it is stable ...

******************************************************************************/

namespace cartan {

unsigned long blockSize(realform::RealForm rf, realform::RealForm drf, 
			const CartanClasses& ccl)

/*
  Synopsis: returns the size of the block defined by the real form rf and the
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

  for (BitMap::iterator i = b.begin(); i != b.end(); ++i) {
    const CartanClass& cc = ccl.cartan(*i);
    unsigned long ci = cc.orbitSize();
    ci *= ccl.fiberSize(rf,*i);
    ci *= ccl.dualFiberSize(drf,*i);
    c += ci;
  }

  return c;
}

unsigned long kgbSize(realform::RealForm rf, const CartanClasses& ccl)

/*
  Synopsis: returns the cardinality of K\G/B for this real form.

  Explanation: this is the cardinality of the one-sided parameter corresponding
  to any strong real form over rf.
*/

{
  using namespace bitmap;
  using namespace cartanclass;

  BitMap b = ccl.support(rf);
  unsigned long c = 0;

  BitMap::iterator b_end = b.end();

  for (BitMap::iterator i = b.begin(); i != b.end(); ++i) {
    const CartanClass& cc = ccl.cartan(*i);
    unsigned long ci = cc.orbitSize();
    ci *= ccl.fiberSize(rf,*i);
    c += ci;
  }

  return c;
}

}

/*****************************************************************************

        Chapter IV --- Auxiliary functions local to this module

  ... explain here when it is stable ...

******************************************************************************/

namespace {

bool checkDecomposition(const weyl::WeylWord& wi, 
			const weyl::WeylWord& ww, 
			const rootdata::RootList& so,
			const weyl::WeylGroup& W, 
			const rootdata::RootDatum& rd)

/*
  Synopsis: checks whether wi decomposes as the composition of the
  cross-action defined by ww followed by the Cayley transform defined by so.
*/

{
  using namespace rootdata;
  using namespace weyl;

  WeylElt w;

  // cross action part
  for (size_t j = 0; j < ww.size(); ++j)
    W.twistedConjugate(w,ww[j]);

  // cayley transform part
  for (size_t j = 0; j < so.size(); ++j) {
    WeylWord rww;
    toWeylWord(rww,so[j],rd);
    for (size_t i = rww.size(); i;) {
      --i;
      W.leftProd(w,rww[i]);
    }
  }

  WeylWord wtest;
  W.involutionOut(wtest,w);
    
  return wtest == wi;
}

}

}
