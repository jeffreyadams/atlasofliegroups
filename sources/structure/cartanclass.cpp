/*!
\file
  \brief Implementation for the CartanClass and Fiber classes.  
*/
/*
  This is cartanclass.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "cartanclass.h"

#include <cassert>
#include <map>
#include <set>
#include <utility>

#include "dynkin.h"
#include "lietype.h"
#include "rootdata.h"
#include "smithnormal.h"
#include "tori.h"
#include "weyl.h"
#include "weylsize.h"
#include "tags.h"

/*****************************************************************************

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace {

  void pause() {;}

  using namespace cartanclass;

  void makeOrbitSize(size::Size&, const rootdata::RootDatum&, 
		     const CartanClass&);
  void makeSimpleComplex(rootdata::RootList&, const rootdata::RootDatum&,
			 const CartanClass&);

  /*!
  \brief derived class of Fiber, to carry out the construction of Fiber. 
  */
class Helper:public Fiber {

private:

// extra data
  const rootdata::RootDatum* d_rootDatum;

  /*!
  \brief Image of the map from the fiber group to the adjoint fiber group.
  */
  latticetypes::ComponentSubquotient d_fiberImage;
  
public:
// constructors and destructors
  Helper(const rootdata::RootDatum&, const latticetypes::LatticeMatrix&);

  virtual ~Helper();

// accessors
  void adjointInvolution(latticetypes::LatticeMatrix&) const;

  size_t imaginaryRank() const {
    return d_simpleImaginary.size();
  }

  const rootdata::RootDatum& rootDatum() const {
    return *d_rootDatum;
  }

// manipulators
  void adjointFiberGroup();

  void adjointMAlpha();

  void baseGrading();

  void complexRootSet();

  void fiberGroup();

  void gradingGroup();

  void gradingShifts();

  void imaginaryRootSet();

  void mAlpha();

  void makeToAdjoint();

  void realFormPartition();

  void realRootSet();

  void rootInvolution();

  void strongReal();

  void strongRepresentatives();

  void weakReal();
};

  /*!
  \brief Constructs a function object giving the action of simple
  imaginary reflections on a Fiber group coset.

  The function object takes two unsigned long arguments.  The first
  argument s indexes a simple imaginary root.  The second argument x
  indexes an element f_x of the Fiber group.  The returned value y
  indicates that imaginary reflection \#s carries the coset f_x.b to
  the coset f_y.b.  (Here b is the base point of the coset, inducing
  the grading d_baseGrading on the imaginary roots.)
  */  
  class FiberAction {

  private:
    gradings::Grading d_baseGrading;
    gradings::GradingList d_gradingShift;
    bitset::RankFlagsList d_mAlpha;

  public:
  // constructors and destructors
    FiberAction() {}

    FiberAction(const gradings::Grading& gr, const gradings::GradingList& gs,
		const bitset::RankFlagsList& ma)
      :d_baseGrading(gr), d_gradingShift(gs), d_mAlpha(ma) {};

    ~FiberAction() {}

  // accessors
    unsigned long operator() (unsigned long, unsigned long) const;

    bool grading(bitset::RankFlags x, unsigned long s) const {
      return d_baseGrading[s] ^ x.scalarProduct(d_gradingShift[s]);
    }
  };

}

/*****************************************************************************

        Chapter I -- The CartanClass class

  ... explain here when it is stable ...

******************************************************************************/

namespace cartanclass {

CartanClass::CartanClass(const rootdata::RootDatum& rd,
			 const latticetypes::LatticeMatrix& q)
   :d_fiber(rd,q)

   

/*!
  Synopsis: constructs the Cartan class with involution q.

  Precondition: ti is the corresponding twisted involution;
*/

{
  using namespace latticetypes;
  using namespace rootdata;
  using namespace tags;

  // make the dual fiber
  RootDatum rdd(rd,DualTag());
  LatticeMatrix qd = q;
  qd.transpose();
  qd.negate();

  d_dualFiber = Fiber(rdd,qd);

  // make the simple complex roots
  makeSimpleComplex(d_simpleComplex,rd,*this);

  // fill in the orbit size
  makeOrbitSize(d_orbitSize,rd,*this);
}

/******** copy and assignment ************************************************/

void CartanClass::swap(CartanClass& other)

{
  d_fiber.swap(other.d_fiber);
  d_dualFiber.swap(other.d_dualFiber);
  d_simpleComplex.swap(other.d_simpleComplex);
  std::swap(d_orbitSize,other.d_orbitSize);

  return;
}

/******** accessors **********************************************************/

bool CartanClass::isMostSplit(unsigned long c) const

/*!
  \brief Tells whether this cartan class is the most split one for
  weak real form \#c.

  Algorithm: this is the case iff the grading corresponding to c is trivial.
*/

{
  using namespace gradings;

  unsigned long x = fiber().weakReal().classRep(c);

  Grading gr;
  fiber().grading(gr,x);

  return gr.none();
}

}

/*****************************************************************************

        Chapter II -- The Fiber class

  ... explain here when it is stable ...

******************************************************************************/

namespace cartanclass {

Fiber::Fiber(const rootdata::RootDatum& rd, 
	     const latticetypes::LatticeMatrix& q)
  :d_torus(0)

/*!
  \brief Puts in f the fiber for this Cartan class.

  Precondition: q contains the root datum involution for this fiber.
*/

{
  Helper fhelp(rd,q);
  swap(fhelp);
}

Fiber::~Fiber ()

{
  delete d_torus;
}

// copy and assignment

Fiber::Fiber(const Fiber& other)
  :d_torus(other.d_torus),
   d_complex(other.d_complex),
   d_imaginary(other.d_imaginary),
   d_real(other.d_real),
   d_simpleImaginary(other.d_simpleImaginary),
   d_rootInvolution(other.d_rootInvolution),
   d_fiberGroup(other.d_fiberGroup),
   d_adjointFiberGroup(other.d_adjointFiberGroup),
   d_toAdjoint(other.d_toAdjoint),
   d_gradingGroup(other.d_gradingGroup),
   d_mAlpha(other.d_mAlpha),
   d_adjointMAlpha(other.d_adjointMAlpha),
   d_baseGrading(other.d_baseGrading),
   d_gradingShift(other.d_gradingShift),
   d_baseNoncompact(other.d_baseNoncompact),
   d_noncompactShift(other.d_noncompactShift),
   d_weakReal(other.d_weakReal),
   d_realFormPartition(other.d_realFormPartition),
   d_strongReal(other.d_strongReal),
   d_strongRealFormReps(other.d_strongRealFormReps)

{
  using namespace tori;

  if (d_torus) // get own copy
    d_torus = new RealTorus(*d_torus);
}

Fiber& Fiber::operator= (const Fiber& other)

/*!
  Synopsis: assignment operator.

  Use copy constructor. This requires a check for self-assignment, or the
  source would be destroyed!
*/

{
  // handle self-assignment
  if (&other != this) {
    this->~Fiber();
    new(this) Fiber(other);
  }

  return *this;
}

void Fiber::swap(Fiber& other)

{  
  std::swap(d_torus,other.d_torus);
  d_complex.swap(other.d_complex);
  d_imaginary.swap(other.d_imaginary);
  d_real.swap(other.d_real);
  d_simpleImaginary.swap(other.d_simpleImaginary);
  d_rootInvolution.swap(other.d_rootInvolution);
  d_fiberGroup.swap(other.d_fiberGroup);
  d_adjointFiberGroup.swap(other.d_adjointFiberGroup);
  d_toAdjoint.swap(other.d_toAdjoint);
  d_gradingGroup.swap(other.d_gradingGroup);
  d_mAlpha.swap(other.d_mAlpha);
  d_adjointMAlpha.swap(other.d_adjointMAlpha);
  d_baseGrading.swap(other.d_baseGrading);
  d_gradingShift.swap(other.d_gradingShift);
  d_baseNoncompact.swap(other.d_baseNoncompact);
  d_noncompactShift.swap(other.d_noncompactShift);
  d_weakReal.swap(other.d_weakReal);
  d_realFormPartition.swap(other.d_realFormPartition);
  d_strongReal.swap(other.d_strongReal);
  d_strongRealFormReps.swap(other.d_strongRealFormReps);

  return;
}

/******** accessors **********************************************************/

void Fiber::compactRootSet(rootdata::RootSet& rs, unsigned long x) const

/*!
  \brief Flags in rs the compact imaginary roots for elt \#x in the
  adjoint fiber.

  Precondition: x represents an element of the subquotient in 
  adjointFiberGroup.
*/

{
  noncompactRootSet(rs,x);
  ~rs;
  rs &= d_imaginary;

  return;
}

void Fiber::grading(gradings::Grading& gr, unsigned long x) const

/*!
  \brief Flags in gr the noncompact simple imaginary roots for elt \#x in 
  the adjoint fiber.

  Precondition: x represents an element of the subquotient in 
  adjointFiberGroup.
*/

{
  using namespace bitset;
  using namespace latticetypes;

  RankFlags b(x);
  Component v(b,adjointFiberRank());

  gr = d_baseGrading;

  for (size_t j = 0; j < v.size(); ++j)
    if (v.test(j))
      gr ^= d_gradingShift[j];

  return;
}

unsigned long Fiber::gradingRep(const gradings::Grading& gr) const

/*!
  \brief Returns an element of the adjoint fiber group corresponding to gr.

  Algorithm: we know that the grading can be reached from the base grading
  with all ones by adding a linear combination of grading shifts.  [Not
  yet implemented.]
*/

{
  using namespace gradings;

  return 0;
}

const latticetypes::LatticeMatrix& Fiber::involution() const

/*!  
  \brief Returns the matrix of the involution on the weight lattice
  of the Cartan subgroup.

  NOTE: this is not inlined to avoid a dependency upon tori.h in the .h file.
*/

{
  return d_torus->involution();
}

void Fiber::mAlpha(latticetypes::Component& ma, const rootdata::Root& cr) const

/*!
  \brief Puts in ma the m_alpha corresponding to cr.

  Precondition: cr is an imaginary coroot for this Cartan.
*/

{
  using namespace lattice;
  using namespace latticetypes;
  using namespace rootdata;

  ma.resize(d_fiberGroup.dimension());

  Component v;
  mod2(v,cr);
  d_fiberGroup.representative(ma,v);
  d_fiberGroup.toSubquotient(ma);

  return;
}

size_t Fiber::minusRank() const

/*!
  \brief Returns the dimension of the -1 eigenspace of the involution.

  NOTE: this is not inlined to avoid a dependency upon tori.h in the .h file.
*/

{
  return d_torus->minusRank();
}

void Fiber::noncompactRootSet(rootdata::RootSet& rs, unsigned long x) const

/*!
  \brief Flags in rs the noncompact imaginary roots for elt \#x in the 
  adjoint fiber.

  Precondition: x represents an element of the subquotient in 
  adjointFiberGroup.
*/

{
  using namespace bitset;
  using namespace latticetypes;

  RankFlags b(x);
  Component v(b,adjointFiberRank());

  rs = d_baseNoncompact;

  for (size_t j = 0; j < v.size(); ++j)
    if (v.test(j))
      rs ^= d_noncompactShift[j];

  return;
}

size_t Fiber::plusRank() const

/*!
  \brief Returns the dimension of the +1 eigenspace of the involution.

  NOTE: this is not inlined to avoid a dependency upon tori.h in the .h file.
*/

{
  return d_torus->plusRank();
}

unsigned long Fiber::toAdjoint(unsigned long x) const

/*!
  \brief Returns the image of x in the adjoint fiber group.

  Precondition: x is a valid element in the fiber group.
*/

{   
  using namespace bitset;
  using namespace latticetypes;

  RankFlags xf(x);
  Component v(xf,fiberRank());

  Component w(adjointFiberRank());
  d_toAdjoint.apply(w,v);
  
  return w.data().to_ulong();
}

unsigned long Fiber::toWeakReal(unsigned long c, size_t rfc) const

/*!
  \brief Returns the class number in the weak real form partition of the
  strong real form \#c in real form class rfc.

  The pair (c,rfc) is the software representation of an equivalence
  class of strong real forms (always assumed to induce tau on H). The
  integer rfc labels an element of Z^delta/[(1+delta)Z], thought of as
  a possible square value for strong real forms.  The fiber group acts
  simply transitively on strong real forms with square equal to rfc.
  The integer c labels an orbit of W_i on this fiber group coset; this
  orbit is the equivalence class of strong real forms.
  
  This function computes the weak real form (W_i orbit on the adjoint
  fiber group) corresponding to (c,rfc).  

  First, rfc also labels a coset of the fiber group image in the
  adjoint fiber group. The coset rfc is a collection of weak real
  forms (W_i orbits on the adjoint fiber group).  The integer brf
  indexes the base W_i orbit on rfc, and "by" is the base adjoint
  fiber group element in brf.

  The fiber group coset corresponding to rfc is labelled by the fiber
  group itself using a base point with image by.  The integer x is a
  representative in the fiber group of the orbit rfc.  Its image in
  the adjoint fiber group is y.  Translating y by the base point "by"
  gives an adjoint fiber group element representing the weak real form
  we want.
*/

{
  using namespace partition;

  // get the basepoint for this real form class
  unsigned long brf = d_realFormPartition.classRep(rfc);
  unsigned long by = d_weakReal.classRep(brf);

  // get representative of class \#c
  unsigned long x = d_strongReal[rfc].classRep(c);

  // find its image in the adjoint fiber group, and translate
  unsigned long y = toAdjoint(x);
  y ^= by;

  return d_weakReal(y);
}

}

/*****************************************************************************

        Chapter III -- The Helper class

  ... explain here when it is stable ...

******************************************************************************/

namespace {

Helper::Helper(const rootdata::RootDatum& rd, 
	       const latticetypes::LatticeMatrix& q)
  :d_rootDatum(&rd)

/*!
  \brief Does the actual fiber construction.

  Precondition: q contains the root datum involution for this fiber.
*/

{
  using namespace latticetypes;
  using namespace rootdata;
  using namespace tori;

  // construct the torus
  d_torus = new RealTorus(q);

  // make the root involution
  rootInvolution();

  // make root sets
  complexRootSet();
  imaginaryRootSet();
  realRootSet();
  
  // find simple imaginary roots
  RootList rl(d_imaginary.begin(),d_imaginary.end());
  rootBasis(d_simpleImaginary,rl,rd);

  // make fiber groups
  fiberGroup();
  adjointFiberGroup();
  gradingGroup();

  // make m_alpha elements for the imaginary roots
  mAlpha();
  adjointMAlpha();

  // make base grading and grading shifts
  baseGrading();
  gradingShifts();

  // make fiber group map
  makeToAdjoint();

  // make weak real form partition  
  weakReal();

  // make real form classes and strong real forms
  pause();
  realFormPartition();
  strongReal();

  // make strong real form representatives
  strongRepresentatives();
}

Helper::~Helper()

{}

/******** accessors *******************************************************/

void Helper::adjointInvolution(latticetypes::LatticeMatrix& q) const

/*!
  \brief Puts in q the matrix of the negative involution in the simple 
  coweight basis.
*/

{
  using namespace latticetypes;
  using namespace rootdata;

  const RootDatum& rd = rootDatum();

  // write involution in root basis
  RootList rl;

  for (size_t s = 0; s < rd.semisimpleRank(); ++s)
    rl.push_back(d_rootInvolution[rd.simpleRootNbr(s)]);

  WeightList b;

  toRootBasis(rl.begin(),rl.end(),back_inserter(b),rd.simpleRootList(),rd);
  q = LatticeMatrix(b);

  // do negative transpose
  q.transpose();
  q.negate();

  return;
}

/******** manipulators *******************************************************/

void Helper::adjointFiberGroup()

/*!
  \brief Makes the group that acts 1-transitively on the adjoint fiber.

  Algorithm: this is the topology subquotient for the negative transpose
  of the involution induced by tau on the root lattice. We express it in
  the basis of simple coweights (mod 2).
*/

{
  using namespace latticetypes;
  using namespace rootdata;
  using namespace tori;

  // write involution in root basis
  LatticeMatrix qsr;
  adjointInvolution(qsr);

  // construct subquotient
  dualPi0(d_adjointFiberGroup,qsr);

  return;
}

void Helper::adjointMAlpha()

/*!
  \brief Constructs the m_alpha's for the adjoint group, and alpha
  simple imaginary.

  Algorithm: the cocharacter lattice for the adjoint group is the coweight
  lattice. To get the coordinates of an element in the basis of X_* /2X_*
  afforded by the simple coweights, it is enough to pair the corresponding
  coroot with the simple roots. The resulting element is automatically tau-
  invariant. Then what is left to do is to take its image modulo V_- for
  the adjoint group.
*/

{
  using namespace lattice;
  using namespace latticetypes;
  using namespace rootdata;

  const RootDatum& rd = rootDatum();
  size_t n = rd.semisimpleRank();
  d_adjointMAlpha.resize(imaginaryRank(),Component(n));

  for (size_t i = 0; i < imaginaryRank(); ++i) {
    Component v(n);
    // compute pairing with simple roots
    for(size_t j = 0; j < n; ++j) {
      RootNbr r = d_simpleImaginary[i];
      LatticeCoeff c = scalarProduct(rd.coroot(r),rd.simpleRoot(j));
      if (c & 1) // scalar product is odd
	v.set(j);
    }
    // check that v is tau-invariant
    // for testing purposes only
    LatticeMatrix q;
    adjointInvolution(q);
    ComponentMap q2;
    mod2(q2,q);
    Component v1(v.size());
    q2.apply(v1,v);
    assert(v.data() == v1.data());
    // reduce modulo V_-
    d_adjointFiberGroup.representative(d_adjointMAlpha[i],v);
    d_adjointFiberGroup.toSubquotient(d_adjointMAlpha[i]);
  }

  return;
}

void Helper::baseGrading()

/*!
  \brief Fills in d_baseGrading and d_baseNoncompact

  Explanation: the base grading is the one where all the simple imaginary
  roots are noncompact. The baseNoncompact RootSet flags the corresponding
  noncompact roots.

  Algorithm: for the noncompact roots, express imaginary roots in terms of
  the simple ones, and see if the sum of coordinates is even.
*/

{
  using namespace gradings;
  using namespace lattice;
  using namespace latticetypes;
  using namespace rootdata;

  const RootDatum& rd = rootDatum();

  d_baseGrading.set();
  d_baseNoncompact.resize(rd.numRoots());

  // baseGrading is already filled with ones
  d_baseGrading.truncate(imaginaryRank());

  // express imaginary roots in simple imaginary basis
  RootList irl(d_imaginary.begin(),d_imaginary.end());

  WeightList ir;
  toRootBasis(irl.begin(),irl.end(),back_inserter(ir),d_simpleImaginary,rd);
  ComponentList ir2;
  mod2(ir2,ir);

  for (size_t j = 0; j < ir2.size(); ++j)
    if (ir2[j].count() & 1ul)
      d_baseNoncompact.insert(irl[j]);

  return;
}

void Helper::complexRootSet()

/*!
  \brief Flags in d_complex the set of complex roots.
*/

{
  d_complex.resize(rootDatum().numRoots());

  for (size_t j = 0; j < d_rootInvolution.size(); ++j)
    if (d_rootInvolution[j] == j)
      continue;
    else if (d_rootInvolution[j] == rootDatum().rootMinus(j))
      continue;
    else
      d_complex.insert(j);

  return;
}

void Helper::gradingGroup()

/*!
  \brief Makes the stabilizer of the grading in the adjoint fiber group.

  Explanation: each real form defines a grading of the imaginary root system,
  obtained by pairing with the simple imaginary roots. It should be the case
  that the grading entirely defines the real form, i.e., the corresponding
  W_i-orbit. However, this doesn't mean that the map from real form parameters
  to gradings has to be injective. The grading group contains the kernel of
  this map at the level of the adjoint fiber group.
*/

{
  using namespace bitset;
  using namespace lattice;
  using namespace latticetypes;
  using namespace rootdata;

  const RootDatum& rd = rootDatum();

  // define map
  const ComponentList& baf = d_adjointFiberGroup.space().basis();

  WeightList bsi;
  toRootBasis(d_simpleImaginary.begin(),d_simpleImaginary.end(),
	      back_inserter(bsi),rd.simpleRootList(),rd);
  ComponentList bsi2;
  mod2(bsi2,bsi);

  ComponentList b;

  // run through the subquotient basis representatives
  for (RankFlags::iterator j = d_adjointFiberGroup.support().begin(); j(); 
       ++j) {
    Component v(imaginaryRank());
    for (size_t i = 0; i < imaginaryRank(); ++i)
      if (scalarProduct(baf[*j],bsi2[i]))
	v.set(i);
    b.push_back(v);
  }

  ComponentMap q(b);

  // find kernel
  ComponentList ker;
  q.kernel(ker);

  d_gradingGroup = ComponentSubspace(ker,adjointFiberRank());
  assert(d_gradingGroup.dimension() == 0);

  return;
}

void Helper::gradingShifts()

/*!
  \brief Fills in d_gradingShifts and d_noncompactShifts.

  Explanation: d_gradingShift[j] contains the action of the j'th basis
  vector of the canonical basis of the adjoint fiber group, on the
  various simple imaginary roots. Similarly, d_noncompactShift[j] contains
  the action on all imaginary roots. With these data, it is easy to
  compute the grading for an arbitrary parameter.

  Algorithm: the basis for the adjoint fiber group is expressed in the
  simple weight basis. Therefore, in order to apply a root, it is enough
  to express that root in the simple root basis, and do the scalar product
  mod 2.
*/

{
  using namespace bitset;
  using namespace gradings;
  using namespace lattice;
  using namespace latticetypes;
  using namespace rootdata;

  const RootDatum& rd = rootDatum();

  // express imaginary roots in simple root basis
  RootList irl(d_imaginary.begin(),d_imaginary.end());
  WeightList ir;
  toRootBasis(irl.begin(),irl.end(),back_inserter(ir),rd.simpleRootList(),rd);
  ComponentList ir2;
  mod2(ir2,ir);

  RankFlags supp = d_adjointFiberGroup.support();
  const ComponentList b = d_adjointFiberGroup.space().basis();

  for (RankFlags::iterator i = supp.begin(); i(); ++i) {
    RootSet rs(rd.numRoots());
    for (size_t j = 0; j < ir2.size(); ++j)
      if (scalarProduct(b[*i],ir2[j]))
	rs.insert(irl[j]);
    d_noncompactShift.push_back(rs);
  }

  // express simple imaginary roots in simple root basis
  const RootList& sil = d_simpleImaginary;
  WeightList si;
  toRootBasis(sil.begin(),sil.end(),back_inserter(si),rd.simpleRootList(),rd);
  ComponentList si2;
  mod2(si2,si);

  for (RankFlags::iterator i = supp.begin(); i(); ++i) {
    Grading gr;
    for (size_t j = 0; j < si2.size(); ++j)
      if (scalarProduct(b[*i],si2[j]))
	gr.set(j);
    d_gradingShift.push_back(gr);
  }

  return;
}

void Helper::fiberGroup()

/*!
  \brief Makes the group that acts 1-transitively on each subset of the
  fiber with fixed central square.

  Explanation: this is the topology subquotient for the negative transpose
  of the involution.
*/

{ 
  using namespace latticetypes;
  using namespace tori;

  LatticeMatrix q = involution();

  // do negative transpose
  q.transpose();
  q.negate();

  // construct subquotient
  dualPi0(d_fiberGroup,q);

  return;
}

void Helper::imaginaryRootSet()

/*!
  \brief Flags in d_imaginary the set of imaginary roots.
*/

{
  d_imaginary.resize(rootDatum().numRoots());

  for (size_t j = 0; j < d_rootInvolution.size(); ++j)
    if (d_rootInvolution[j] == j)
      d_imaginary.insert(j);

  return;
}

void Helper::mAlpha()

/*!
  \brief Constructs the m_alpha's.

  Explanation: in fact, we simply write down the expression of the coroots
  corresponding to the imaginary simple roots, and compute modulo 2X_*.
*/

{
  using namespace lattice;
  using namespace latticetypes;
  using namespace rootdata;

  const RootDatum& rd = rootDatum();
  size_t n = rd.rank();
  d_mAlpha.resize(imaginaryRank(),Component(n));

  for (size_t j = 0; j < imaginaryRank(); ++j) {
    Component v(n);
    mod2(v,rd.coroot(d_simpleImaginary[j]));
    d_fiberGroup.representative(d_mAlpha[j],v);
    d_fiberGroup.toSubquotient(d_mAlpha[j]);
  }

  return;
}

void Helper::makeToAdjoint()

/*!
  \brief Fills in the toAdjoint matrix.

  Explanation: this is the matrix of the projection from the fiber group
  to the adjoint fiber group.
*/

{
  using namespace bitset;
  using namespace lattice;
  using namespace latticetypes;
  using namespace rootdata;

  const ComponentList b = d_fiberGroup.space().basis();
  const RankFlags& supp = d_fiberGroup.support();
  const RootDatum& rd = rootDatum();
  size_t n = rd.semisimpleRank();

  ComponentList b_sr(n);
  for (size_t j = 0; j < b_sr.size(); ++j)
    mod2(b_sr[j],rd.simpleRoot(j));

  ComponentList b_ad;

  for (RankFlags::iterator i = supp.begin(); i(); ++i) {
    Component v(n);
    for(size_t j = 0; j < n; ++j) {
      if (scalarProduct(b[*i],b_sr[j]))
	v.set(j);
    }
    d_adjointFiberGroup.representative(v,v);
    d_adjointFiberGroup.toSubquotient(v);
    b_ad.push_back(v);
  }

  d_toAdjoint = ComponentMap(b_ad);

  ComponentList b_id;
  initBasis(b_id,d_adjointFiberGroup.dimension());
  d_fiberImage = ComponentSubquotient(b_id,b_ad,
				      d_adjointFiberGroup.dimension());

  return;
}

void Helper::realFormPartition()

/*!
  \brief Makes the partition of the real forms according to central square
  classes.

  Explanation: the classes are indexed by cosets of the image of the
  fiber group in the adjoint fiber group. This also corresponds to the
  various possible values of (1+delta)(x) modulo (1+delta)(Z).

  Algorithm: cl is the partition of the weak real forms according to
  the element of Z(G)^delta/[(1+delta)Z(G)] given by the square of the
  preimage in G^Gamma.  The class of weak real form \#j is given by
  cl[j], which is an integer less than the number of elements of order
  2 in Z(G).

  A weak real form is an orbit of W_im on the adjoint fiber group.
  Because the W_im action lifts to the fiber group, an orbit is
  contained in a single coset of im(Fiber).  Two weak real forms
  define the same central square class if and only if they lie in the
  same coset of im(Fiber).

  For each j, y is the number of a representative in the adjoint fiber
  group for W_im orbit \#j.  Component v is first defined as equal to
  adjoint fiber group element \#y.  Then v is redefined to be the
  canonical coset representative for v in adjointFiber/[image of
  Fiber]. Then cl[j] is the number of element v.  The UnnormalizedTag
  on the Partition constructor call refers to the fact that the
  Partition d_fiberImage uses the values of cl to number the classes. 
*/

{
  using namespace bitset;
  using namespace latticetypes;
  using namespace partition;
  using namespace tags;

  std::vector<unsigned long> cl;
  /* 
   */
  cl.resize(d_weakReal.classCount());

  for (size_t j = 0; j < cl.size(); ++j) {
    unsigned long y = d_weakReal.classRep(j);
    Component v = Component(RankFlags(y),d_adjointFiberGroup.dimension());
    d_fiberImage.representative(v,v);
    d_fiberImage.toSubquotient(v);
    cl[j] = v.data().to_ulong();
  }

  d_realFormPartition = Partition(cl,UnnormalizedTag());

  return;
}

void Helper::realRootSet()

/*!
  \brief Flags in d_real the set of real roots.
*/

{
  d_real.resize(rootDatum().numRoots());

  for (size_t j = 0; j < d_rootInvolution.size(); ++j)
    if (d_rootInvolution[j] == rootDatum().rootMinus(j))
      d_real.insert(j);

  return;
}

void Helper::rootInvolution()

/*!
  \brief Makes the permutation of the roots corresponding to the torus
  involution.
*/

{
  using namespace latticetypes;
  using namespace rootdata;

  const RootDatum& rd = rootDatum();

  d_rootInvolution.resize(rd.numRoots());

  for (size_t j = 0; j < rd.numRoots(); ++j) {
    Weight v(rd.rank());
    involution().apply(v,rd.root(j));
    d_rootInvolution[j] = rd.rootNbr(v);
  }

  return;
}

void Helper::strongReal()

/*!
  \brief Writes the strong real form partitions.

  Algorithm: there is one strong real form partition for each real
  form class in d_weakReal. What we do is walk through the real form
  partition; the first time we encounter a new class, we make the
  corresponding partition, taking a basepoint corresponding to the
  canonical grading representative for our real form, and we make the
  map from strong real forms to weak ones.

  The number of partitions in d_strongReal is
  d_realFormPartition.classCount(), the number of possible values z of
  the square of a strong real form in Z.  (We count z only modulo
  (1+delta)Z.  Changing z by (1+delta)Z changes nothing we want to
  compute.)  Each partition is of the Fiber group; the number n
  computed before makeOrbits below is the cardinality of the Fiber
  group.
*/

{
  using namespace bitset;
  using namespace bitvector;
  using namespace gradings;
  using namespace latticetypes;

  d_strongReal.resize(d_realFormPartition.classCount());

  // get the grading shifts; these are obtained from the images of the
  // canonical basis vectors of the fiber group in the adjoint fiber group,
  // and then transposing
  GradingList gst(fiberRank());

  for (size_t j = 0; j < gst.size(); ++j) {
    Component v(adjointFiberRank());
    d_toAdjoint.column(v,j);
    combination(gst[j],d_gradingShift,v.data());
  }

  GradingList gs(imaginaryRank());

  for (size_t j = 0; j < gs.size(); ++j)
    for (size_t i = 0; i < gst.size(); ++i)
      if (gst[i][j])
	gs[j].set(i);

  // get the m_alpha's
  RankFlagsList ma;

  for (size_t j = 0; j < imaginaryRank(); ++j)
    ma.push_back(d_mAlpha[j].data());

  // make the various orbit pictures
  for (size_t j = 0; j < d_strongReal.size(); ++j) {
    unsigned long rf = d_realFormPartition.classRep(j);
    unsigned long y = d_weakReal.classRep(rf);
    Grading bg;
    grading(bg,y);
    FiberAction fa(bg,gs,ma);
    size_t n = 1ul << fiberRank();
    makeOrbits(d_strongReal[j],fa,imaginaryRank(),n);
  }

  return;
}

void Helper::strongRepresentatives()

/*!
  \brief Fills in the strong real form representatives.

  Algorithm: we run through the various strong real form partitions, and
  for each, we pick as representative the first element lying above a weak real
  form representative.
*/

{
  using namespace bitset;
  using namespace latticetypes;
  using namespace partition;

  ComponentList b(fiberRank());

  for (size_t j = 0; j < b.size(); ++j)
    d_toAdjoint.column(b[j],j);

  d_strongRealFormReps.resize(numRealForms());

  for (size_t rf = 0; rf < numRealForms(); ++rf) {

    // find adjoint representative of rf
    unsigned long y = d_weakReal.classRep(rf);
    RankFlags yf(y);

    // subtract base point
    size_t c = d_realFormPartition(rf);
    size_t rf0 = d_realFormPartition.classRep(c);
    unsigned long y0 = d_weakReal.classRep(rf0);
    RankFlags yf0(y0);

    // find preimage of y in the fiber
    yf ^= yf0;
    Component v(yf,adjointFiberRank());
    RankFlags xf;
    firstSolution(xf,b,v); // there has to be a solution!
    unsigned long x = xf.to_ulong();

    // make representative
    d_strongRealFormReps[rf] = std::make_pair(x,c);
  }

  return;
}

void Helper::weakReal()

/*!
  \brief Writes in d_weakReal the partition of the adjoint fiber 
  corresponding to weak real forms.

  Algorithm: we construct the FiberAction object corresponding to the adjoint
  fiber.  The weak real forms correspond to the orbits of the
  imaginary Weyl group on the adjoint Fiber group.  The group is
  generated by the simple imaginary reflections
*/

{
  using namespace bitset;
  using namespace lattice;
  using namespace latticetypes;
  using namespace gradings;
  using namespace rootdata;

  // grading shifts are given by the "transpose" of the gradingShift data

  GradingList gs(imaginaryRank());

  for (size_t j = 0; j < imaginaryRank(); ++j)
    for (size_t i = 0; i < d_gradingShift.size(); ++i)
      if (d_gradingShift[i][j])
	gs[j].set(i);

  RankFlagsList ma;

  for (size_t j = 0; j < imaginaryRank(); ++j)
    ma.push_back(d_adjointMAlpha[j].data());

  // make orbits
  FiberAction fa(d_baseGrading,gs,ma);
  size_t n = 1ul << adjointFiberRank();
  makeOrbits(d_weakReal,fa,imaginaryRank(),n);

  return;
}

}

/*****************************************************************************

        Chapter IV -- The FiberAction class

  ... explain here when it is stable ...

******************************************************************************/

namespace {

unsigned long FiberAction::operator() (unsigned long s, unsigned long x) const

/*!
  \brief Abstract action operator for the imaginary Weyl group on
  the fiber group.

  Algorithm: we interpret the bits of x as the coordinates of the element in
  the fiber group, in terms of the chosen basepoint. Then the action is
  adding m_alpha[s] if the grading of alpha_s at x is noncompact, and doing
  nothing otherwise.
*/

{
  using namespace bitset;

  RankFlags b(x);

  if (grading(b,s))
    b ^= d_mAlpha[s];

  return b.to_ulong();

}

}

/*****************************************************************************

        Chapter V -- Functions declared in cartanclass.cpp

  ... explain here when it is stable ...

******************************************************************************/

namespace cartanclass {

void compactTwoRho(latticetypes::Weight& tr, unsigned long x,
		   const Fiber& f, const rootdata::RootDatum& rd)

/*!
  \brief Puts in tr the sum of positive compact imaginary roots for x in f.
*/

{
  using namespace rootdata;

  RootSet rs;
  f.compactRootSet(rs,x);
  twoRho(tr,rs,rd);

  return;
}

void restrictGrading(gradings::Grading& gr, const rootdata::RootSet& rs,
		     const rootdata::RootList& rl)

/*!
  \brief Flags in gr the restriction of the grading in rs to rl.
*/

{
  gr.reset();

  for (size_t j = 0; j < rl.size(); ++j)
    gr.set(j,rs.isMember(rl[j]));

  return;
}

void specialGrading(gradings::Grading& gr, const cartanclass::Fiber& f, 
		    realform::RealForm rf, const rootdata::RootDatum& rd)

/*!
  \brief Puts in gr a grading in the orbit corresponding to rf, with the
  smallest possible number of noncompact roots.

  Precondition: f is the fundamental fiber;

  Explanation: for each noncompact noncomplex irreducible real form, there
  is at least one grading with exactly one noncompact simple root. Our choice
  amounts to a grading which induces one of the aforementioned ones on each
  noncompact noncomplex simple factor. The function of this is to enable
  easy type recognition.

  NOTE : the grading is represented in terms of simple roots for the root
  system rd. This is ok; knowledge of the gradings of those roots is enough
  to define the real form.
*/

{
  using namespace bitset;
  using namespace cartanclass;
  using namespace gradings;
  using namespace rootdata;

  std::set<Grading,GradingCompare> grs;
  unsigned long n = f.adjointFiberSize();

  // sort the gradings that occur in this class
  for (unsigned long j = 0; j < n; ++j) {
    if (f.weakReal()(j) != rf)
      continue;
    RootSet rs;
    f.noncompactRootSet(rs,j);
    Grading grx;
    restrictGrading(grx,rs,rd.simpleRootList());
    grs.insert(grx);
  }

  // return the first element
  gr = *(grs.begin());

  return;
}

void toMostSplit(rootdata::RootList& so, const cartanclass::Fiber& fundf, 
		 realform::RealForm rf, const rootdata::RootDatum& rd)

/*!
  \brief Puts in so a set of strongly orthogonal roots for ccl, leading from
  the fundamental Cartan to the most split one for the strong real form rf

  Algorithm: it is very simple: as long as there are noncompact imaginary
  roots, use one of them to get to a less compact Cartan.
*/

{
  using namespace cartanclass;
  using namespace rootdata;

  so.clear();
  RootSet ir = fundf.imaginaryRootSet();
  RootSet rs;
  unsigned long rep = fundf.weakReal().classRep(rf);
  fundf.noncompactRootSet(rs,rep);

  while (not rs.empty()) {
    RootNbr rn = rs.front();
    for (RootSet::iterator i = ir.begin(); i(); ++i) {
      if (not rd.isOrthogonal(rn,*i)) {
	ir.remove(*i);
	rs.remove(*i);
      } else if (sumIsRoot(rn,*i,rd))
	rs.flip(*i);
    }
    so.push_back(rn);
  }

  strongOrthogonalize(so,rd);

  return;
}

}

/*****************************************************************************

        Chapter VI -- Auxiliary functions

  ... explain here when it is stable ...

******************************************************************************/

namespace {

void makeOrbitSize(size::Size& os, const rootdata::RootDatum& rd, 
		   const CartanClass& cc)

/*!
   \brief Puts in os the size of the twisted involution orbit for this class.
*/

{
  using namespace dynkin;
  using namespace latticetypes;
  using namespace lietype;
  using namespace size;
  using namespace weylsize;

  LatticeMatrix cm;
  LieType lt;

  // put size of full Weyl group in os

  cartanMatrix(cm,rd);
  lieType(lt,cm);
  weylSize(os,lt);

  // divide by product of imaginary, real and complex sizes

  Size ws;

  cartanMatrix(cm,cc.simpleImaginary(),rd);
  lieType(lt,cm);
  weylSize(ws,lt);
  os /= ws;

  cartanMatrix(cm,cc.simpleReal(),rd);
  lieType(lt,cm);
  weylSize(ws,lt);
  os /= ws;

  cartanMatrix(cm,cc.simpleComplex(),rd);
  lieType(lt,cm);
  weylSize(ws,lt);
  os /= ws;

  return;
}

void makeSimpleComplex(rootdata::RootList& sc, const rootdata::RootDatum& rd,
		       const CartanClass& cc)

/*!
  \brief Puts in sc the list of the simple roots for a complex factor in
  W^tau.

  Explanation: W^tau is the semidirect product of W^R x W^{iR} (Weyl groups
  of the real and imaginary root systems), with the diagonal subgroup of
  W^C, where W^C is the Weyl group of the root system orthogonal to both the
  half sum of real and of imaginary roots. That root system is complex for
  the involution induced by tau; we put in sc a basis of "half" of it (in
  other words, the complex root system is the disjoint union of the root
  system generated by sc, and its image under tau.)

  NOTE: there was a bad bug here in an earlier version, which amounted to the
  assumption that the standard positive root system for the Phi^C is tau-
  stable; this is very false.
*/

{
  using namespace bitset;
  using namespace dynkin;
  using namespace latticetypes;
  using namespace rootdata;

  RootList rl;

  Weight tri(rd.rank());
  twoRho(tri,cc.imaginaryRootSet(),rd);
  Weight trr(rd.rank());
  twoRho(trr,cc.realRootSet(),rd);

  for (size_t j = 0; j < rd.numRoots(); ++j)
    if (rd.isOrthogonal(tri,j) and 
	rd.isOrthogonal(trr,j))
      rl.push_back(j);

  RootList rb;
  rootBasis(rb,rl,rd);

  LatticeMatrix cm;
  cartanMatrix(cm,rb,rd);

  DynkinDiagram dd(cm);

  RankFlags b;
  b.set();
  b.truncate(dd.rank());

  while (b.any()) {
    size_t s = b.firstBit();
    RankFlags c = dd.component(s);
    b.andnot(c);
    for (RankFlags::iterator i = c.begin(); i(); ++i) {
      sc.push_back(rb[*i]);
      // find image of rb[*i] under the involution
      RootNbr rTau = cc.rootInvolution(rb[*i]);
      // erase all elements of b that are not orthogonal to rTau
      for (RankFlags::iterator j = b.begin(); j(); ++j)
	if (not rd.isOrthogonal(rb[*j],rTau))
	  b.reset(*j);
    }
  }

  return;
}

}

}
