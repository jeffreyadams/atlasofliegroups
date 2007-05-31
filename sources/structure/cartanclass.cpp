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

namespace atlas {

namespace {

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
  \brief Quotient of |d_adjointFiberGroup| by the image of the map
     |d_toAdjoint| from the fiber group to the adjoint fiber group.
  */
  latticetypes::SmallSubquotient d_fiberImage;

public:
// constructors and destructors
  Helper(const rootdata::RootDatum&, const latticetypes::LatticeMatrix&);

  ~Helper();

private:
// accessors (private because only needed during construction)

  size_t imaginaryRank() const {
    return simpleImaginary().size();
  }

  const rootdata::RootDatum& rootDatum() const {
    return *d_rootDatum;
  }

  latticetypes::SmallSubquotient fiberGroup() const;

  latticetypes::SmallSubquotient adjointFiberGroup() const;

  latticetypes::SmallSubspace gradingGroup() const;

  latticetypes::LatticeMatrix adjointInvolution() const;

  gradings::Grading baseGrading(rootdata::RootSet& flagged_roots) const;

  gradings::GradingList gradingShifts(rootdata::RootSetList& all_shifts) const;

  latticetypes::SmallBitVectorList adjointMAlpha() const;

  latticetypes::SmallBitVectorList mAlpha() const;

// manipulators
  void makeToAdjoint();

  void weakReal();

  void realFormPartition();

  void strongReal();

  void strongRepresentatives();
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

        Chapter I -- The InvolutionData and CartanClass classes

******************************************************************************/

namespace cartanclass {

InvolutionData::InvolutionData(const rootdata::RootDatum& rd,
			       const latticetypes::LatticeMatrix& q)
  : d_rootInvolution(rd.rootPermutation(q))
  , d_imaginary(rd.numRoots()), d_real(rd.numRoots()), d_complex(rd.numRoots())
  , d_simpleImaginary()
{

  for (size_t j = 0; j < d_rootInvolution.size(); ++j)
    if (d_rootInvolution[j] == j)
      d_imaginary.insert(j);
    else if (d_rootInvolution[j] == rd.rootMinus(j))
      d_real.insert(j);
    else
      d_complex.insert(j);

  // find simple imaginary roots
  rootBasis(d_simpleImaginary,imaginary_roots(),rd);
}

void InvolutionData::swap(InvolutionData& other)
{
  d_rootInvolution.swap(other.d_rootInvolution);
  d_imaginary.swap(other.d_imaginary);
  d_real.swap(other.d_real);
  d_complex.swap(other.d_complex);
  d_simpleImaginary.swap(other.d_simpleImaginary);
}

/*!
  Synopsis: constructs the Cartan class with involution q.

  Precondition: q is an involution of the root datum rd that can be obtained
  from the distinguished involution by multiplication to the left by the
  action of a Weyl group element w, which is called the twisted involution
  associated to q.
*/
CartanClass::CartanClass(const rootdata::RootDatum& rd,
			 const latticetypes::LatticeMatrix& q)
   : d_fiber(rd,q)
   , d_dualFiber(rd,q,tags::DualTag())
{
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
}

/******** accessors **********************************************************/

bool CartanClass::isMostSplit(unsigned long c) const

/*!
  \brief Tells whether this cartan class is the most split one for
  weak real form corresponding to class \#c in fiber().weakReal().

  Algorithm: this is the case iff the grading corresponding to c is trivial,
  i.e., all imaginary roots are compact.
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

******************************************************************************/

namespace cartanclass {

/*!
  \brief Constructs the fiber for the Cartan class determined by the
  involution |q| of the root datum |rd|.
*/

Fiber::Fiber(const rootdata::RootDatum& rd,
	     const latticetypes::LatticeMatrix& q)
  : d_torus(NULL)
  , d_involutionData(rd,q) // dummy value will disappear after swap
  , d_fiberGroup()
  , d_adjointFiberGroup()
  , d_toAdjoint()
  , d_baseGrading()
  , d_gradingShift()
  , d_baseNoncompact()
  , d_noncompactShift()
  , d_weakReal()
  , d_realFormPartition()
  , d_strongReal()
  , d_strongRealFormReps()
{
  Helper fhelp(rd,q);
  swap(fhelp);
}

/*!
  \brief Constructs the dual fiber for the Cartan class determined by the
  involution q of the root datum rd.
*/
Fiber::Fiber(const rootdata::RootDatum& rd,
	     const latticetypes::LatticeMatrix& q,
	     tags::DualTag)
  : d_torus(NULL)
  , d_involutionData(rd,q) // dummy value will disappear after swap
  , d_fiberGroup()
  , d_adjointFiberGroup()
  , d_toAdjoint()
  , d_baseGrading()
  , d_gradingShift()
  , d_baseNoncompact()
  , d_noncompactShift()
  , d_weakReal()
  , d_realFormPartition()
  , d_strongReal()
  , d_strongRealFormReps()
{
  rootdata::RootDatum rdd(rd,tags::DualTag());
  latticetypes::LatticeMatrix qd = q;
  qd.transpose();
  qd.negate();

  Helper fhelp(rdd,qd);
  swap(fhelp);
}

Fiber::~Fiber ()

{
  delete d_torus;
}

// copy and assignment

Fiber::Fiber(const Fiber& other)
  : d_torus(other.d_torus)
  , d_involutionData(other.d_involutionData)
  , d_fiberGroup(other.d_fiberGroup)
  , d_adjointFiberGroup(other.d_adjointFiberGroup)
  , d_toAdjoint(other.d_toAdjoint)
  , d_baseGrading(other.d_baseGrading)
  , d_gradingShift(other.d_gradingShift)
  , d_baseNoncompact(other.d_baseNoncompact)
  , d_noncompactShift(other.d_noncompactShift)
  , d_weakReal(other.d_weakReal)
  , d_realFormPartition(other.d_realFormPartition)
  , d_strongReal(other.d_strongReal)
  , d_strongRealFormReps(other.d_strongRealFormReps)
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
  d_involutionData.swap(other.d_involutionData);
  d_fiberGroup.swap(other.d_fiberGroup);
  d_adjointFiberGroup.swap(other.d_adjointFiberGroup);
  d_toAdjoint.swap(other.d_toAdjoint);
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
  rs &= imaginaryRootSet();

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
  SmallBitVector v(b,adjointFiberRank());

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

void Fiber::mAlpha(latticetypes::SmallBitVector& ma, const rootdata::Root& cr)
 const

/*!
  \brief Puts in ma the m_alpha corresponding to cr.

  Precondition: cr is an imaginary coroot for this Cartan.
*/

{
  ma.resize(d_fiberGroup.dimension());

  latticetypes::SmallBitVector v(cr);
  ma= d_fiberGroup.toBasis(v);
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
  SmallBitVector v(b,adjointFiberRank());

  rs = d_baseNoncompact;

  for (size_t j = 0; j < v.size(); ++j)
    if (v.test(j))
      rs ^= d_noncompactShift[j];
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
  SmallBitVector v(xf,fiberRank());

  SmallBitVector w(adjointFiberRank());
  d_toAdjoint.apply(w,v);

  return w.data().to_ulong();
}

unsigned long Fiber::toWeakReal(unsigned long c, size_t rfc) const

/*!
  \brief Returns the class number in the weak real form partition of the
  strong real form \#c in real form class \#rfc.

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

} // namespace cartanclass

/*****************************************************************************

        Chapter III -- The Helper class of the Fiber class

******************************************************************************/

namespace {

/*!
  \brief Does the actual fiber construction.

  Precondition: q contains the root datum involution for this fiber.
*/

Helper::Helper(const rootdata::RootDatum& rd,
	       const latticetypes::LatticeMatrix& q)
  : Fiber(new tori::RealTorus(q),rd,q) // auxiliary constructor for base object
  , d_rootDatum(&rd)
  , d_fiberImage()
{
  // make fiber groups
  d_fiberGroup=fiberGroup();
  d_adjointFiberGroup=adjointFiberGroup();
  assert(gradingGroup().dimension()==0);

  // make base grading, its noncompact root set, and grading shifts
  d_baseGrading=baseGrading(d_baseNoncompact);
  d_gradingShift=gradingShifts(d_noncompactShift);

  // make fiber group map
  makeToAdjoint();

  // make weak real form partition
  weakReal();

  // make real form classes and strong real forms
  realFormPartition();
  strongReal();

  // make strong real form representatives
  strongRepresentatives();
}

Helper::~Helper()

{}

/******** accessors *******************************************************/

latticetypes::SmallSubquotient Helper::fiberGroup() const

/*!
  \brief returns the group that acts 1-transitively on each subset of the
  fiber with fixed central square.

  Explanation: the value returned is the subquotient $V_+ + V_-/V_+$ as
  constructed by |tori::dualPi0|, but for the negative transpose of our
  involution.
*/

{
  latticetypes::LatticeMatrix q = involution();

  // do negative transpose
  q.transpose();
  q.negate();

  // construct subquotient
  latticetypes::SmallSubquotient result;
  tori::dualPi0(result,q);

  return result;
}

latticetypes::LatticeMatrix Helper::adjointInvolution() const

/*!
  \brief Returns the matrix of the transformation induced by negative of
  our involution on the coweight lattice, expressed on the dual basis of the
  simple root basis, which is the simple coweight basis.

  So we first transform the involution to one on the root basis, and then take
  the negative transpose as in the |fiberGroup| method.
*/

{
  using namespace latticetypes;
  using namespace rootdata;
  const RootDatum& rd = rootDatum();

  // write involution in root basis
  RootList rl; // root numbers of images of simple roots by our involution

  for (size_t s = 0; s < rd.semisimpleRank(); ++s)
    rl.push_back(involution_image_of_root(rd.simpleRootNbr(s)));

  WeightList b; // will contain members of |rl| expressed in simple roots

  toRootBasis(rl.begin(),rl.end(),back_inserter(b),rd.simpleRootList(),rd);
  latticetypes::LatticeMatrix q = LatticeMatrix(b);

  // do negative transpose
  q.transpose();
  q.negate();

  return q;
}

latticetypes::SmallSubquotient Helper::adjointFiberGroup() const

/*!
  \brief Makes the group that acts 1-transitively on the adjoint fiber.

  Algorithm: this is the subquotient $V_+ + V_-/V_+$ for the negative
  transpose of the involution induced by tau on the root lattice (which is
  computed by |adjointInvolution|).
*/

{
  // write involution in root basis
  latticetypes::LatticeMatrix qsr= adjointInvolution();

  // construct subquotient
  latticetypes::SmallSubquotient result; tori::dualPi0(result,qsr);

  assert(result.rank()==rootDatum().semisimpleRank());

  return result;
}

latticetypes::SmallSubspace Helper::gradingGroup() const

/*!
  \brief Makes the stabilizer of the grading in the adjoint fiber group.

  Explanation: each real form defines a grading of the imaginary root system,
  obtained by pairing with the simple imaginary roots. It should be the case
  that the grading entirely defines the real form, i.e., the corresponding
  W_i-orbit. However, this doesn't mean that the map from real form parameters
  to gradings has to be injective. It could a priori be the case that some
  real form parameters in the same W_i-orbit give identical gradings of the
  imaginary root system; this would mean that the nonzero adjoint fiber
  element that transforms one into the other gives a null pairing with every
  imaginary root. It seems that this never happens, whence the |assert| below.
  Nevertheless, and although never used, we compute the kernel of the map to
  gradings shifts in the adjoint fiber group, calling it the "grading group".
*/

{
  using namespace bitset;
  using namespace lattice;
  using namespace latticetypes;
  using namespace rootdata;

  const rootdata::RootDatum& rd = rootDatum();

  // define map
  const SmallBitVectorList& baf = d_adjointFiberGroup.space().basis();

  WeightList bsi; // express simple imaginary roots on (full) simple roots
  toRootBasis(simpleImaginary().begin(),simpleImaginary().end(),
	      back_inserter(bsi),rd.simpleRootList(),rd);


  SmallBitVectorList bsi2(bsi); // and reduce mod 2

  SmallBitVectorList b;

  /* set b[j][i] = <e_j,bsi2[i]> where e_j=baf[j'], with |baf[j']| the
     subquotient basis representative number |j| */
  for (RankFlags::iterator j = d_adjointFiberGroup.support().begin();
       j(); ++j)
  {
    SmallBitVector v(imaginaryRank());
    for (size_t i = 0; i < imaginaryRank(); ++i)
      v.set(i,scalarProduct(baf[*j],bsi2[i]));
    b.push_back(v);
  }

  BinaryMap q(b);

  // find kernel
  SmallBitVectorList ker; q.kernel(ker);

  return SmallSubspace(ker,adjointFiberRank());
}


gradings::Grading Helper::baseGrading(rootdata::RootSet& flagged_roots) const

/*!
  \brief Returns the base grading (all ones) on simple imaginary roots, while
  flagging all noncompact imaginary roots in |flagged_roots|

  Algorithm: for the noncompact roots, express imaginary roots in terms of
  the simple ones, and see if the sum of coordinates is even.
*/

{
  // express imaginary roots in simple imaginary basis
  rootdata::RootList irl(imaginaryRootSet().begin(),imaginaryRootSet().end());

  latticetypes::WeightList ir;
  rootdata::toRootBasis(irl.begin(),irl.end(),back_inserter(ir)
			,simpleImaginary(),rootDatum());
  latticetypes::SmallBitVectorList ir2(ir); // |ir2.size()==irl.size()()|

  // now flag all roots with noncompact grading
  flagged_roots.resize(rootDatum().numRoots());

  for (size_t j = 0; j < irl.size(); ++j)
    flagged_roots.set_mod2(irl[j],ir2[j].count());

  return gradings::Grading(constants::lMask[imaginaryRank()]); // all ones
}

gradings::GradingList Helper::gradingShifts(rootdata::RootSetList& all_shifts)
  const

/*!
  \brief Computes grading shifts for simple imaginary roots, and also sets
  |all_shifts| to flag the grading of the full set of imaginary roots.

  Explanation: component |j| of the result contains the grading on the set of
  simple imaginary roots given by the canonical basis vector |j| of the
  adjoint fiber group. This is not the grading of a real form, but the amount
  by which the grading changes by the action of that basis vector, whence the
  name grading shitf. Similarly, |all_shifts[j]| is set to contain the action
  on all imaginary roots. With these data and the "base" grading corresponding
  to the real form parameter numbered $0$, it is easy to compute the grading
  for an arbitrary real form parameter (represented by the element of the
  adjoint fiber group moving the base element there).

  Algorithm: the basis for the adjoint fiber group is the reduction modulo 2
  of the simple coweight basis. Therefore, in order to apply to a root, it is
  enough to express that root in the simple root basis, and do the scalar
  product mod 2. Note that in contrast to |baseGrading| above, and in spite
  of the reuse of the same names, we express on the (full) simple root basis
  here, not on the simple imaginary root basis (which by the way is not a
  subset of the former).
*/

{
  using namespace bitset;
  using namespace gradings;
  using namespace lattice;
  using namespace latticetypes;
  using namespace rootdata;

  const rootdata::RootDatum& rd = rootDatum();

  // express imaginary roots in simple root basis
  rootdata::RootList irl(imaginaryRootSet().begin(),imaginaryRootSet().end());

  latticetypes::WeightList ir;
  rootdata::toRootBasis(irl.begin(),irl.end(),back_inserter(ir),
			rd.simpleRootList(),rd);
  latticetypes::SmallBitVectorList ir2(ir); // |ir2.size()==irl.size()|


  // also express simple imaginary roots in simple root basis
  const rootdata::RootList& sil = simpleImaginary();
  latticetypes::WeightList si;
  rootdata::toRootBasis(sil.begin(),sil.end(),back_inserter(si),
			rd.simpleRootList(),rd);
  latticetypes::SmallBitVectorList si2(si); // reduce vectors mod 2

  // now compute all results
  bitset::RankFlags supp = d_adjointFiberGroup.support();
  const SmallBitVectorList& b = d_adjointFiberGroup.space().basis();
  gradings::GradingList result;

  // traverse basis of the subquotient |d_adjointFiberGroup|
  for (RankFlags::iterator i = supp.begin(); i(); ++i) {

    // all imaginary roots part
    RootSet rs(rd.numRoots());
    for (size_t j = 0; j < ir2.size(); ++j)
      rs.set_to(irl[j],scalarProduct(b[*i],ir2[j]));
    all_shifts.push_back(rs);

    // simple imaginary roots part
    Grading gr;
    for (size_t j = 0; j < si2.size(); ++j)
      gr.set(j,scalarProduct(b[*i],si2[j]));
    result.push_back(gr);
  }

  return result;
}

latticetypes::SmallBitVectorList Helper::mAlpha() const

/*!
  \brief Constructs the m_alpha's (images of coroots) in the fiber group,
  for alpha simple imaginary.

  Explanation: in fact, we simply take the coroots corresponding to the
  imaginary simple roots, which in the root datum are already expressed in the
  basis of the coweight lattice $X^*$ dual to the weight lattice, reduce
  modulo 2, and interpret the result in the subquotient |d_fiberGroup| of
  $X^*$.
*/

{
  latticetypes::SmallBitVectorList result;
  const rootdata::RootDatum& rd = rootDatum();
  size_t n = rd.rank();
  result.resize(imaginaryRank(),latticetypes::SmallBitVector(n));

  for (size_t i = 0; i < imaginaryRank(); ++i) {
    latticetypes::SmallBitVector v(rd.coroot(simpleImaginary(i)));
    result[i]= d_fiberGroup.toBasis(v);
  }
  return result;
}

latticetypes::SmallBitVectorList Helper::adjointMAlpha() const

/*!
  \brief Constructs the m_alpha's (images of coroots) in the adjoint fiber
  group, for alpha simple imaginary.

  Algorithm: the cocharacter lattice for the adjoint group is spanned by the
  simple coweights. To get the coordinates of an element in that basis (which
  is dual to that of the simple roots), it is enough to pair it with the
  simple roots. The resulting element for the corot of a simple imaginary
  $\alpha$ is automatically $\tau$-invariant, since $\alpha$ is, so its
  reduction modulo 2 lies in $V_+$ (for the cocharacter lattice). Then what is
  left to do is to convert to the basis of the adjoint fiber group, which
  amounts to reducing modulo $V_-$.
*/

{
  const rootdata::RootDatum& rd = rootDatum();

  latticetypes::SmallBitVectorList result(imaginaryRank());

  for (size_t i = 0; i < imaginaryRank(); ++i) {
    latticetypes::SmallBitVector v(rd.semisimpleRank());
    // compute pairing with simple roots modulo 2
    for(size_t j = 0; j < v.size(); ++j) {
      latticetypes::LatticeCoeff c = latticetypes::scalarProduct
	(rd.coroot(simpleImaginary(i)),rd.simpleRoot(j));
      v.set_mod2(j,c);
    }

    // reduce |v| modulo V_-
    result[i]=d_adjointFiberGroup.toBasis(v);
  }
  return result;
}


/******** manipulators *******************************************************/

void Helper::makeToAdjoint()

/*!
  \brief Fills in the toAdjoint matrix.

  Explanation: this is the matrix of the map from the fiber group to the
  adjoint fiber group. Such a map exists because (1) the roots lattice is a
  sublattice of the weight lattice, so there is a restriction map from the
  coweight lattice $X^*$ (dual to the weight lattice) to the lattice $Y^*$
  spanned by the fundamental coweights (dual to the root lattice), and (2)
  this map sends $X^*_+$ to $Y^*_+$ and $X^*_-$ to $Y^*_-$ and of course
  $2X^*$ to $2Y^*$, so it induces a map on the respective subquotients of the
  form $(V_+ + V_-)/V_-$, where $V_+$ is the modulo $2X^*$ image of $X^*_+$
  (or similarly for $Y^*_-$). While the map $X^* \to Y^*$ is injective (only)
  in the semisimple case, little can be said about the induced map.

  The above argument shows that the vector |v| computed below of scalar
  products of a generator of the fiber group with simple roots can be validly
  incorporated (by calling |d_adjointFiberGroup.toBasis|) into the adjoint
  fiber group; notably, it represents an element of $V_+ + V_-$ in
  $Y^* / 2Y^*$. (without space that formula would have ended our comment!)
*/

{
  using namespace bitset;
  using namespace lattice;
  using namespace latticetypes;
  using namespace rootdata;

  const SmallBitVectorList& b = d_fiberGroup.space().basis();
  const RankFlags& supp = d_fiberGroup.support();
  const RootDatum& rd = rootDatum();
  size_t n = rd.semisimpleRank();

  SmallBitVectorList b_sr(rd.beginSimpleRoot(),rd.endSimpleRoot()); // mod 2

  // images of fiber group basis in adjoint fiber group
  SmallBitVectorList b_ad; b_ad.reserve(d_fiberGroup.dimension());

  for (RankFlags::iterator i = supp.begin(); i(); ++i) {
    SmallBitVector v(n);
    for(size_t j = 0; j < n; ++j) {
      v.set(j,bitvector::scalarProduct(b[*i],b_sr[j]));
    }
    b_ad.push_back(d_adjointFiberGroup.toBasis(v));
  }

  d_toAdjoint = BinaryMap(b_ad);

  SmallBitVectorList b_id;
  bitvector::initBasis(b_id,d_adjointFiberGroup.dimension());

  // construct |d_fiberImage| as quotient of spans of |b_id| (all) and |b_ad|
  // which is a quotient of the adjoint fiber group!
  d_fiberImage = SmallSubquotient(b_id,b_ad,
				      d_adjointFiberGroup.dimension());
  // |d_fiberImage| will be needed in |realFormPartition|
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

  // make m_alpha elements for the imaginary roots
  SmallBitVectorList d_adjointMAlpha=adjointMAlpha();

  RankFlagsList ma;

  for (size_t j = 0; j < imaginaryRank(); ++j)
    ma.push_back(d_adjointMAlpha[j].data());

  // make orbits
  FiberAction fa(d_baseGrading,gs,ma);
  size_t n = 1ul << adjointFiberRank();
  makeOrbits(d_weakReal,fa,imaginaryRank(),n);
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
  std::vector<unsigned long> cl(d_weakReal.classCount());

  for (size_t j = 0; j < cl.size(); ++j) {
    unsigned long y = d_weakReal.classRep(j);
    latticetypes::SmallBitVector v
      (bitset::RankFlags(y),d_adjointFiberGroup.dimension());

    // reduce modulo image of map from fiber group to adjoint fiber group
    cl[j] = d_fiberImage.toBasis(v).data().to_ulong();
  }

  // partition $[0,d_weakReal.classCount()[$ accoring to values of |cl|
  d_realFormPartition = partition::Partition(cl,tags::UnnormalizedTag());
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
    SmallBitVector v(adjointFiberRank());
    d_toAdjoint.column(v,j);
    combination(gst[j],d_gradingShift,v.data());
  }

  GradingList gs(imaginaryRank());

  for (size_t j = 0; j < gs.size(); ++j)
    for (size_t i = 0; i < gst.size(); ++i)
      if (gst[i][j])
	gs[j].set(i);

  SmallBitVectorList d_mAlpha=mAlpha();

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

  SmallBitVectorList b(fiberRank());

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
    SmallBitVector v(yf,adjointFiberRank());
    RankFlags xf;
    firstSolution(xf,b,v); // there has to be a solution!
    unsigned long x = xf.to_ulong();

    // make representative
    d_strongRealFormReps[rf] = std::make_pair(x,c);
  }

  return;
}


} // namespace

/*****************************************************************************

        Chapter IV -- The FiberAction class

******************************************************************************/

namespace {

unsigned long FiberAction::operator() (unsigned long s, unsigned long x) const

/*!
  \brief Abstract action operator for the imaginary Weyl group on
  the fiber group.

  Algorithm: we interpret the bits of x as the coordinates of the element in
  the fiber group, in terms of the chosen basepoint. The choice of the
  basepoint is incorporated in the |d_baseGrading| and |d_gradingShift| fields
  of the |FiberAction| object, and implicitly accessed by the auxiliary methof
  |grading|. Then the action is adding m_alpha[s] if the grading of $\alpha_s$
  given by $x$ is noncompact (i.e., true, odd), and doing nothing otherwise.
*/

{
  bitset::RankFlags b(x);

  if (grading(b,s))
    b ^= d_mAlpha[s];

  return b.to_ulong();

}

}

/*****************************************************************************

        Chapter V -- Functions declared in cartanclass.cpp

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
      RootNbr rTau = cc.involution_image_of_root(rb[*i]);
      // erase all elements of b that are not orthogonal to rTau
      for (RankFlags::iterator j = b.begin(); j(); ++j)
	if (not rd.isOrthogonal(rb[*j],rTau))
	  b.reset(*j);
    }
  }
}

} // namespace

} // namespace atlas
