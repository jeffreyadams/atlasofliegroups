/*!
\file
  \brief Implementation for the CartanClass and Fiber classes.
*/
/*
  This is cartanclass.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "cartanclass.h"

#include <cassert>
#include <map>
#include <set>
#include <utility>
#include <stdexcept>

#include "dynkin.h"
#include "lietype.h"
#include "rootdata.h"
#include "smithnormal.h"
#include "tori.h"
#include "weyl.h"
#include "weylsize.h"
#include "tags.h"
#include "complexredgp.h"

namespace atlas {

namespace {

  using namespace cartanclass;

  /*!
  \brief Constructs a function object defining the action of the simple
  imaginary reflections on a fiber.

  The function object this class provides, and that can be used by
  |partition::makeOrbits| takes the index |s| of an imaginary root and a
  number |x| encoding an element of the fiber; it return a number |y|
  similarly encoding the image of that element under the action of |s|.

  In fact the number |x| describes (in binary form) the element of the fiber
  group that translates the base point to fiber element in question, and the
  interpretation of |y| is the same. The action is determined by (1) the
  grading of the simple imaginary roots associated to the chosen base point of
  the fiber, (2) the grading shifts associated with the generators of the
  fiber group, and (3) the vectors |d_mAlpha[s]| by which each of the simple
  imaginary roots \f$\alpha_s\f$ translates if it acts act non-trivially (this
  is the image in the fiber group of the coroot of \f$\alpha_s\f$). The action
  of \f$\alpha_s\f$ will translate |x| by |d_mAlpha[s]| if the grading at
  \f$\alpha_s\f$ associated to |x| is odd (noncompact). To facilitate the
  determination of that grading, the grading shift information is stored as
  bitsets |d_alpha[s]| for each simple imaginary root. Since in |Fiber|, the
  grading shifts are organised by generator of the fiber group, the neceesary
  "transposition" of the information is done by the constructor below.
  */

  class FiberAction {

  private:
    gradings::Grading d_baseGrading; // indexed by |s|
    bitset::RankFlagsList d_alpha;   // indexed by |s| first
    bitset::RankFlagsList d_mAlpha;  // indexed by |s| first

  public:
  // constructors and destructors
    FiberAction(const gradings::Grading& gr,
		const gradings::GradingList& gs,
		const bitset::RankFlagsList& ma)
      : d_baseGrading(gr)
      , d_alpha(ma.size()) // initialise size only here
      , d_mAlpha(ma)
    {
      for (size_t i = 0; i<ma.size(); ++i)
	for (size_t j = 0; j<gs.size(); ++j)
	    d_alpha[i].set(j,gs[j][i]);
    }

  // accessors
    bool grading(bitset::RankFlags x, unsigned long s) const {
      return d_baseGrading[s] ^ x.scalarProduct(d_alpha[s]);
    }
    unsigned long operator() (unsigned long s, unsigned long x) const
    {
      bitset::RankFlags b(x); if (grading(b,s))  b ^= d_mAlpha[s];
      return b.to_ulong();
    }
  };

} // namespace

/*****************************************************************************

        Chapter I -- The InvolutionData and CartanClass classes

******************************************************************************/

namespace cartanclass {

InvolutionData::InvolutionData(const rootdata::RootDatum& rd,
			       const latticetypes::LatticeMatrix& q)
  : d_rootInvolution(rd.rootPermutation(q))
  , d_imaginary(rd.numRoots()) // we can only dimension root sets for now
  , d_real(rd.numRoots())
  , d_complex(rd.numRoots())
  , d_simpleImaginary()        // here even dimensioning is pointless
  , d_simpleReal()
{
  for (size_t j = 0; j < d_rootInvolution.size(); ++j)
    if (d_rootInvolution[j] == j)
      d_imaginary.insert(j);
    else if (d_rootInvolution[j] == rd.rootMinus(j))
      d_real.insert(j);
    else
      d_complex.insert(j);

  // find simple imaginary roots
  d_simpleImaginary=rd.simpleBasis(imaginary_roots());
  d_simpleReal=rd.simpleBasis(real_roots());
}

InvolutionData::InvolutionData(const complexredgp::ComplexReductiveGroup& G,
			       const weyl::TwistedInvolution& tw)
  : d_rootInvolution()
  , d_imaginary(G.rootDatum().numRoots())
  , d_real(G.rootDatum().numRoots())
  , d_complex(G.rootDatum().numRoots())
  , d_simpleImaginary()
  , d_simpleReal()
{
  const rootdata::RootDatum& rd=G.rootDatum();
  { // follow the inner class twist by the action of |tw| on roots
    std::vector<rootdata::RootNbr> simple_image(rd.semisimpleRank());
    for (size_t i=0; i<rd.semisimpleRank(); ++i)
      simple_image[i]=G.twisted_root(rd.simpleRootNbr(i));

    weyl::WeylWord ww=G.weylGroup().word(tw.w());
    for (size_t i=ww.size(); i-->0;)
      rd.simple_root_permutation(ww[i]).left_mult(simple_image);

    d_rootInvolution=rd.extend_to_roots(simple_image);
  }

  for (size_t j = 0; j < d_rootInvolution.size(); ++j)
    if (d_rootInvolution[j] == j)
      d_imaginary.insert(j);
    else if (d_rootInvolution[j] == rd.rootMinus(j))
      d_real.insert(j);
    else
      d_complex.insert(j);

  // find simple imaginary roots
  d_simpleImaginary=rd.simpleBasis(imaginary_roots());
  d_simpleReal=rd.simpleBasis(real_roots());
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
  , d_dualFiber(cartanclass::dualFiber(rd,q)) // call non-member |dualFiber|
  , d_simpleComplex(makeSimpleComplex(rd))
  , d_orbitSize(makeOrbitSize(rd))
{
  assert(d_dualFiber.simpleReal()==d_fiber.simpleImaginary());
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

/*!
  \brief Tells whether this cartan class is the most split one for
  weak real form corresponding to class \#wrf in fiber().weakReal().

  Algorithm: this is the case iff the grading corresponding to wrf is trivial,
  i.e., all imaginary roots are compact. In this case it does not matter which
  representative fiber element |x| is chosen, since the action of the
  imaginary Weyl group obviously stabilises the trivial grading. [In fact
  there should be only one such representative by injectivity of the map from
  the adjoint fiber to gradings, as given by the assert statement below. MvL]
*/
bool CartanClass::isMostSplit(adjoint_fiber_orbit wrf) const
{
  AdjointFiberElt x = fiber().weakReal().classRep(wrf);

  gradings::Grading gr= fiber().grading(x);

  assert(gr.any() or fiber().weakReal().classSize(wrf)==1);

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
  : d_torus(new tori::RealTorus(q))
  , d_involutionData(rd,q)
  , d_fiberGroup(makeFiberGroup())
  , d_adjointFiberGroup(makeAdjointFiberGroup(rd))
  , d_baseNoncompact()
  , d_baseGrading(makeBaseGrading(d_baseNoncompact,rd))
  , d_noncompactShift()
  , d_gradingShift(makeGradingShifts(d_noncompactShift,rd))
  , d_toAdjoint(makeFiberMap(rd))
  , d_weakReal(makeWeakReal(rd))
  , d_realFormPartition(makeRealFormPartition())
  , d_strongReal(makeStrongReal(rd))
  , d_strongRealFormReps(makeStrongRepresentatives())
{
  assert(gradingGroup(rd).dimension()==0);
}

/*!
  \brief Constructs the dual fiber for the Cartan class determined by the
  involution |q| of the root datum |rd|.

  This is the fiber for the dual root datum and the negative transpose matrix.

  This used to be a separate constructor, but it is more practical to make it
  a function external to the class. This gives us the opportunity to prepare
  the dual root datum and involution locally, and let the main constructor do
  the real work; making this an official constructor would have given us the
  obligation to initialise the member data to irrelevant values, before being
  able to compute the dual root datum and matrix and start the construction.
*/
Fiber dualFiber
  (const rootdata::RootDatum& rd, const latticetypes::LatticeMatrix& q)
{
  latticetypes::LatticeMatrix qd = q; qd.transpose(); qd.negate();

  return Fiber(rootdata::RootDatum(rd,tags::DualTag()),qd);
}

Fiber::~Fiber ()

{
  delete d_torus;
}

// copy and assignment

Fiber::Fiber(const Fiber& other)
  : d_torus(other.d_torus==NULL ? NULL : new tori::RealTorus(*other.d_torus))
  , d_involutionData(other.d_involutionData)
  , d_fiberGroup(other.d_fiberGroup)
  , d_adjointFiberGroup(other.d_adjointFiberGroup)
  , d_baseNoncompact(other.d_baseNoncompact)
  , d_baseGrading(other.d_baseGrading)
  , d_noncompactShift(other.d_noncompactShift)
  , d_gradingShift(other.d_gradingShift)
  , d_toAdjoint(other.d_toAdjoint)
  , d_weakReal(other.d_weakReal)
  , d_realFormPartition(other.d_realFormPartition)
  , d_strongReal(other.d_strongReal)
  , d_strongRealFormReps(other.d_strongRealFormReps)
{}

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
  std::swap(d_torus,other.d_torus); // swap pointers here
  d_involutionData.swap(other.d_involutionData);
  d_fiberGroup.swap(other.d_fiberGroup);
  d_adjointFiberGroup.swap(other.d_adjointFiberGroup);
  d_gradingShift.swap(other.d_gradingShift);
  d_baseGrading.swap(other.d_baseGrading);
  d_noncompactShift.swap(other.d_noncompactShift);
  d_baseNoncompact.swap(other.d_baseNoncompact);
  d_toAdjoint.swap(other.d_toAdjoint);
  d_weakReal.swap(other.d_weakReal);
  d_realFormPartition.swap(other.d_realFormPartition);
  d_strongReal.swap(other.d_strongReal);
  d_strongRealFormReps.swap(other.d_strongRealFormReps);
}

/*       Private accessors of |Fiber| used during construction      */


/*!
  \brief returns the group that acts 1-transitively on each subset of the
  fiber with fixed central square.

  Explanation: the value returned is the subquotient $V_+ + V_-/V_+$ as
  constructed by |tori::dualPi0|, but for the negative transpose of our
  involution. Note that only the involution is used, not the root datum.
*/
latticetypes::SmallSubquotient Fiber::makeFiberGroup() const
{
  // construct subquotient
  latticetypes::SmallSubquotient result;
  tori::dualPi0(result,involution().negative_transposed());

  return result;
}

/*!
  \brief Returns the matrix of the transformation induced by negative of
  our involution on the cocharacter lattice for the adjoint group, expressed
  on the simple coweight basis, which is dual to the simple root basis.

  So we first transform the involution to one on the root basis, and then take
  the negative transpose.
*/
latticetypes::LatticeMatrix Fiber::adjointInvolution
  (const rootdata::RootDatum& rd) const
{
  // write involution in root basis
  rootdata::RootList rl; // root numbers of involution images of simple roots

  for (size_t s = 0; s < rd.semisimpleRank(); ++s)
    rl.push_back(involution_image_of_root(rd.simpleRootNbr(s)));

  latticetypes::WeightList b; // members of |rl| expressed in simple roots

  rd.toRootBasis(rl.begin(),rl.end(),back_inserter(b));

  return latticetypes::LatticeMatrix(b).negative_transposed();
}

/*!
  \brief Makes the group that acts 1-transitively on the adjoint fiber.

  Algorithm: this is the subquotient $V_+ + V_-/V_+$ for the negative
  transpose of the involution induced by \f$\tau\f$ on the root lattice (which
  is computed by |adjointInvolution|).
*/
latticetypes::SmallSubquotient Fiber::makeAdjointFiberGroup
  (const rootdata::RootDatum& rd) const
{
  // construct subquotient
  latticetypes::SmallSubquotient result;
  tori::dualPi0(result,adjointInvolution(rd));

  assert(result.rank()==rd.semisimpleRank());

  return result;
}

/*!
  \brief Makes the stabilizer of the grading in the adjoint fiber group.

  Explanation: each real form parameter defines a grading of the imaginary
  root system, obtained by pairing with the simple imaginary roots. It should
  be the case that the grading entirely defines the real form, i.e., the
  corresponding W_i-orbit. However, this doesn't mean that the map from real
  form parameters to gradings has to be injective. It could a priori be the
  case that some real form parameters in the same W_i-orbit give identical
  gradings of the imaginary root system; this would mean that the nonzero
  adjoint fiber element that transforms one into the other gives a null
  pairing with every imaginary root. It appears that this never happens,
  whence the constructor for |Fiber| just calls this function to |assert| that
  the result is trivial.
*/
latticetypes::SmallSubspace Fiber::gradingGroup
  (const rootdata::RootDatum& rd) const
{
  // define map
  const latticetypes::SmallBitVectorList& baf =
    d_adjointFiberGroup.space().basis();

  // express simple imaginary roots on (full) simple roots
  latticetypes::WeightList bsi;
  rd.toRootBasis(simpleImaginary().begin(),simpleImaginary().end(),
		 back_inserter(bsi));


  latticetypes::SmallBitVectorList bsi2(bsi); // and reduce mod 2

  latticetypes::SmallBitVectorList b;

  /* set b[j][i] = <e_j,bsi2[i]> where e_j=baf[j'], with |baf[j']| the
     subquotient basis representative number |j| */
  for (bitset::RankFlags::iterator j = d_adjointFiberGroup.support().begin();
       j(); ++j)
  {
    latticetypes::SmallBitVector v(imaginaryRank());
    for (size_t i = 0; i < imaginaryRank(); ++i)
      v.set(i,scalarProduct(baf[*j],bsi2[i]));
    b.push_back(v);
  }

  latticetypes::BinaryMap q(b);

  // find kernel
  latticetypes::SmallBitVectorList ker; q.kernel(ker);

  return latticetypes::SmallSubspace(ker,adjointFiberRank());
}


/*!
  \brief Returns the base grading (all ones) on simple imaginary roots, while
  flagging all noncompact imaginary roots in |flagged_roots|

  Algorithm: for the noncompact roots, express imaginary roots in terms of
  the simple ones, and see if the sum of coordinates is even.
*/
gradings::Grading Fiber::makeBaseGrading
  (rootdata::RootSet& flagged_roots,const rootdata::RootDatum& rd) const

{
  // express all imaginary roots in simple imaginary basis
  rootdata::RootList irl(imaginaryRootSet().begin(),imaginaryRootSet().end());
  latticetypes::WeightList ir;
  rd.toRootBasis(irl.begin(),irl.end(),back_inserter(ir),simpleImaginary());

  // now flag all roots with noncompact grading
  flagged_roots.set_capacity(rd.numRoots());

  for (size_t j = 0; j < irl.size(); ++j)
  {
    latticetypes::Weight v=ir[j];
    int count=0;
    for (size_t i=0; i<v.size(); ++i) count+=v[i]; // add coefficients
    flagged_roots.set_mod2(irl[j],count);          // and take result mod 2
  }
  return gradings::Grading(constants::lMask[imaginaryRank()]); // all ones
}

/*!
  \brief Computes, for each basis vector of the adjoint fiber group, the
  grading shifts for simple imaginary roots, and also sets |all_shifts| to
  flag the grading of the full set of imaginary roots.

  Explanation: component |j| of the result contains the grading on the set of
  simple imaginary roots given by the canonical basis vector |j| of the
  adjoint fiber group. This is not the grading of a real form, but the amount
  by which the grading changes by the action of that basis vector, whence the
  name grading shift. Similarly, |all_shifts[j]| is set to contain the action
  on all imaginary roots. With these data and the "base" grading corresponding
  to the real form parameter numbered $0$, it is easy to compute the grading
  for an arbitrary real form parameter (as represented by the element of the
  adjoint fiber group moving the base element there).

  Algorithm: the basis for the adjoint fiber group is expressed in the
  reduction modulo 2 of the simple coweight basis. Therefore, in order to
  apply to a root, it is enough to express that root in the simple root basis,
  and do the scalar product mod 2. Note that in contrast to |makeBaseGrading|
  above, and in spite of the reuse of the same names, we express in the (full)
  simple root basis here, not on the simple imaginary root basis.
*/
gradings::GradingList Fiber::makeGradingShifts
  (rootdata::RootSetList& all_shifts,const rootdata::RootDatum& rd) const
{
  // express imaginary roots in (full) simple root basis
  rootdata::RootList irl(imaginaryRootSet().begin(),imaginaryRootSet().end());

  latticetypes::WeightList ir;
  rd.toRootBasis(irl.begin(),irl.end(),back_inserter(ir));
  latticetypes::SmallBitVectorList ir2(ir); // |ir2.size()==irl.size()|


  // also express simple imaginary roots in (full) simple root basis
  const rootdata::RootList& sil = simpleImaginary();
  latticetypes::WeightList si;
  rd.toRootBasis(sil.begin(),sil.end(),back_inserter(si));
  latticetypes::SmallBitVectorList si2(si); // reduce vectors mod 2

  // now compute all results
  bitset::RankFlags supp = d_adjointFiberGroup.support();
  const latticetypes::SmallBitVectorList& b
    = d_adjointFiberGroup.space().basis();
  gradings::GradingList result;

  // traverse basis of the subquotient |d_adjointFiberGroup|
  for (bitset::RankFlags::iterator i = supp.begin(); i(); ++i) {

    // all imaginary roots part
    rootdata::RootSet rs(rd.numRoots());
    for (size_t j = 0; j < ir2.size(); ++j)
      rs.set_to(irl[j],scalarProduct(b[*i],ir2[j]));
    all_shifts.push_back(rs);

    // simple imaginary roots part
    gradings::Grading gr;
    for (size_t j = 0; j < si2.size(); ++j)
      gr.set(j,scalarProduct(b[*i],si2[j]));
    result.push_back(gr);
  }

  assert(result.size()==adjointFiberRank());

  return result;
}

/*!
  \brief Constructs the \f$m_\alpha\f$s (images of coroots) in the fiber group,
  for \f$\alpha\f$ simple imaginary.

  The effective number of bits of each \f$m_\alpha\f$ is |d_fiberGroup.dimension()|

  We take the coroots corresponding to the imaginary simple roots, which in
  the root datum are already expressed in the basis of the coweight lattice
  $X^*$ dual to the weight lattice; we reduce the coordinates modulo 2 (this
  is hidden in the call of |toBasis|, which converts its argument to a
  |SmallBitVector| first), and interpret the result in the subquotient
  |d_fiberGroup| of $X^* / 2X^*$.
*/
bitset::RankFlagsList Fiber::mAlphas (const rootdata::RootDatum& rd) const
{
  bitset::RankFlagsList result(imaginaryRank());

  for (size_t i = 0; i<result.size(); ++i)
    result[i]= d_fiberGroup.toBasis(rd.coroot(simpleImaginary(i))).data();

  return result;
}

/*!
  \brief Constructs the \f$m_\alpha\f$s (images of coroots) in the adjoint fiber
  group, for \f$\alpha\f$ simple imaginary.

  The number of bits of each \f$m_\alpha\f$ is |d_adjointFiberGroup.dimension()|

  Algorithm: the cocharacter lattice for the adjoint group is spanned by the
  simple coweights. To get the coordinates of an element in that basis (which
  is dual to that of the simple roots), it is enough to pair it with the
  simple roots. The resulting element for the coroot of a simple imaginary
  \f$\alpha\f$ is automatically \f$\tau\f$-invariant, since \f$\alpha\f$ is,
  so its reduction modulo 2 lies in $V_+$ (for the cocharacter lattice). Then
  what is left to do is to convert to the basis of the adjoint fiber group,
  which amounts to reducing modulo $V_-$.
*/
bitset::RankFlagsList
Fiber::adjointMAlphas (const rootdata::RootDatum& rd) const
{
  bitset::RankFlagsList result(imaginaryRank());

  for (size_t i = 0; i<result.size(); ++i) {
    latticetypes::SmallBitVector v(rd.semisimpleRank());
    // compute pairing with simple roots modulo 2
    for(size_t j = 0; j < v.size(); ++j) {
      latticetypes::LatticeCoeff c =
	rd.coroot(simpleImaginary(i)).scalarProduct(rd.simpleRoot(j));
      v.set_mod2(j,c);
    }

    // reduce |v| modulo $V_-$
    result[i]=d_adjointFiberGroup.toBasis(v).data();
  }
  return result;
}

/*!
  \brief Computes the toAdjoint matrix.

  Explanation: this is the matrix of the natural map from the fiber group to
  the adjoint fiber group. Such a map exists because (1) the root lattice is a
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
latticetypes::BinaryMap
Fiber::makeFiberMap(const rootdata::RootDatum& rd) const
{

  latticetypes::SmallBitVectorList b_sr
    (rd.beginSimpleRoot(),rd.endSimpleRoot()); // mod 2

  // images of fiber group basis in adjoint fiber group
  latticetypes::SmallBitVectorList b_ad;
  b_ad.reserve(d_fiberGroup.dimension());

  const latticetypes::SmallBitVectorList& b = d_fiberGroup.space().basis();
  const bitset::RankFlags& supp = d_fiberGroup.support();
  size_t n = rd.semisimpleRank();
  for (bitset::RankFlags::iterator i = supp.begin(); i(); ++i) {
    latticetypes::SmallBitVector v(n);
    for(size_t j = 0; j < n; ++j) {
      v.set(j,bitvector::scalarProduct(b[*i],b_sr[j]));
    }
    b_ad.push_back(d_adjointFiberGroup.toBasis(v));
  }

  return latticetypes::BinaryMap(b_ad); // convert vectors to a matrix
}


/*!
  \brief Computes the partition of the adjoint fiber, whose parts correspond
  to the weak real forms.

  Algorithm: we construct the FiberAction object corresponding to the action
  of the imaginary Weyl group on the adjoint fiber, and then call |makeOrbits|
  to make the partition. For the fiber action we can use the base grading and
  the grading shifts "as is", while the fiber group elements \f$m_\alpha\f$
  are computed by |adjointMAlpha|.
*/
partition::Partition Fiber::makeWeakReal(const rootdata::RootDatum& rd) const
{
  bitset::RankFlagsList ma=adjointMAlphas(rd);

  // make orbits
  partition::Partition result;
  partition::makeOrbits(result,FiberAction(d_baseGrading,d_gradingShift,ma),
			imaginaryRank(),adjointFiberSize());
  return result;
}



/*!
  \brief Computes the partition of the weak real forms according to central
  square classes.

  Explanation: while weak real forms correspond to orbits in the adjoint
  fiber, they can be grouped into even coarser "central square classes",
  defined by the condition that the corresponding orbits belong to the same
  coset in the adjoint fiber (group) by the image (under |toAdjoint|) of the
  fiber group. The imaginary Weyl group $W_{im}$ acts on both the fiber group
  and the adjoint fiber group, and |toAdjoint| intertwines these actions, so
  that an orbit defining a weak real form is always contained in a single
  coset of the image of |toAdjoint|; therefore we really have a partition of
  the set of weak real forms. Its parts are called central square classes,
  since for these weak real forms any \f$x=g.\delta\in G.\delta\f$ whose
  $H$-conjugacy class defines a fiber element in a strong real form lying over
  the weak real form gives the same value of \f$x^2\in Z(G)\f$ modulo
  \f$(1+\delta)(Z(G))\f$. From the above it follows there are $2^m$ central
  square classes where $m=adjointFiberRank-rank(toAdjoint)$.

  Algorithm: we compute the quotient of the adjoint fiber group by the fiber
  group as a |Subquotient| object whose |space()| is the whole adjoint fiber
  group. Its |toBasis| method will map each weak real form to a vector on the
  basis of the (sub)quotient whose |to_ulong()| value will be used to
  characterise the central square class. Since we want to preserve this
  numbering when construting the resulting partition (rather than reordering
  by smallest element), we call the constructor with a |tags::UnnormalizedTag|
  argument.
*/
partition::Partition Fiber::makeRealFormPartition() const
{
  std::vector<unsigned long> cl(numRealForms());

  latticetypes::SmallBitVectorList b_id;
  bitvector::initBasis(b_id,d_adjointFiberGroup.dimension());

  /* construct |modFiberImage| as quotient of spans of |b_id| (i.e., all) and
     image of |d_toAdjoint|; it is a quotient of the adjoint fiber group! */
  latticetypes::SmallSubquotient modFiberImage
    (b_id,d_toAdjoint.image(),d_adjointFiberGroup.dimension());

  // take representatives of weak real forms and reduce modulo fiber image
  for (size_t j = 0; j < cl.size(); ++j) {
    unsigned long y = d_weakReal.classRep(j);
    latticetypes::SmallBitVector v
      (bitset::RankFlags(y),d_adjointFiberGroup.dimension());

    // reduce modulo image of map from fiber group to adjoint fiber group
    cl[j] = modFiberImage.toBasis(v).data().to_ulong();
  }

  /* partition the set $[0,numRealForms()[$ according to the
     |modFiberImage.size()| distinct values in the image of |cl|
   */
  partition::Partition result(cl,tags::UnnormalizedTag());
  assert(result.classCount()==modFiberImage.size());
  return result;
}


/*!
  \brief Computes the strong real form partitions.

  In the software, strong real forms are represented by orbits of "fiber
  elements" under the imaginary Weyl group, where the fiber elements live in a
  union of affine spaces over $Z/2Z$, one for each central square class.
  Each of these affine spaces has the fiber group as associated vector space,
  but the action of the imaginary Weyl group is different for each one.
  Therefore there is a separate partition of the fiber group for each central
  square class of weak real forms (as determined by |makeRealFormPartition|).
  The function |makeStrongReal| computes all of these partitions.

  Method: We traverse the central square classes of |d_realFormPartition|; for
  each class we choose a weak real form |rf| in it, which itself labels an
  orbit on the adjoint fiber, and then choose a point |y| in that orbit.
  Altogether, the values of |y| are just coset representatives of the image of
  the |toAdjoint| map in the adjoint fiber group. The affine space of fiber
  elements for this class of strong real forms will map to that coset by a map
  the sends the base point (the one represented by the number $0$) to |y|, and
  which induces the linear map |toAdjoint| on the associated vector spaces. It
  follows that we must associate to the base point the same grading as to |y|,
  and this determines the action of the imaginary Weyl group on this affine
  space of fiber elements (together with the fact that the grading shift for
  any vector $v$ in the fiber group is to one associated to $toAdjoint(v)$).
  The orbits of this action then define the strong real form partition
  associated to this central square class.

  The grading defined by |y| is \emph{not} independent of the choices of |rf|
  and |y|. However, if $y'$ is another possibility for |y|, then by
  construction there is an element $d$ of the fiber group such that
  $y'=y+toAdjoint[d]$, and the partitions defined by the gradings of $y$ and
  $y'$ will differ by translation in the fiber group over |d|, in other words
  they correspond to another choice of a base point in the affine space.

*/
std::vector<partition::Partition> Fiber::makeStrongReal
  (const rootdata::RootDatum& rd) const
{
  /* get the grading shifts; these are obtained from the images of the
     canonical basis vectors of the fiber group in the adjoint fiber group. */

  gradings::GradingList gs(fiberRank());

  for (size_t j = 0; j < gs.size(); ++j)
    gs[j]=bitvector::combination(d_gradingShift,d_toAdjoint.column(j).data());

  // get the $m_\alpha$s
  bitset::RankFlagsList ma=mAlphas(rd);

  // make the various partitions

  std::vector<partition::Partition> result(d_realFormPartition.classCount());

  for (square_class j = 0; j < result.size(); ++j) {
    gradings::Grading bg=grading(class_base(j));
    size_t n = fiberSize(); // order of fiber group
    partition::makeOrbits(result[j],FiberAction(bg,gs,ma),imaginaryRank(),n);
  }

  return result;
}

/*
  For each weak real form |wrf| we can choose a fiber element |x| (in the
  affine space corresponding to its central square class) that maps to the
  chosen adjoint fiber element representative of |wrf|. The auxiliary method
  |makeStrongRepresentatives| makes a vector of size |numRealForms())| (to be
  stored in |d_strongRealFormReps|) whose element |wrf| is the pair
  $(x,c)$.
*/
std::vector<StrongRealFormRep> Fiber::makeStrongRepresentatives() const
{
  using namespace bitset;
  using namespace latticetypes;
  using namespace partition;

  latticetypes::SmallBitVectorList b(fiberRank());

  for (size_t j = 0; j < b.size(); ++j)
    b[j]=d_toAdjoint.column(j);

  std::vector<StrongRealFormRep> result(numRealForms());

  for (size_t wrf = 0; wrf<result.size(); ++wrf) {

    size_t c = central_square_class(wrf);

    // find representative |yf| of |wrf| in the adjoint fiber (group)
    bitset::RankFlags yf(d_weakReal.classRep(wrf));

    // subtract base point
    yf ^= bitset::RankFlags(class_base(c));

    // find preimage |xf| of |yf| in the fiber
    latticetypes::SmallBitVector v(yf,adjointFiberRank()); // the desired image

    // solve equation |toAdjoint(xf)=v|
    RankFlags xf;
#ifndef NDEBUG
    bool success=bitvector::firstSolution(xf,b,v);
    assert(success);  // there has to be a solution!
#else
    bitvector::firstSolution(xf,b,v);
#endif

    // make representative
    result[wrf] = std::make_pair(FiberElt(xf.to_ulong()),c);
  }

  return result;
}






/*
		     Public accessors of the |Fiber| class
*/



/*!
  \brief Returns the noncompact imaginary roots for elt \#x in the adjoint
  fiber (whose bitset is interpreted as an element of the subquotient
  in |d_adjointFiberGroup|).
*/
rootdata::RootSet Fiber::noncompactRoots(AdjointFiberElt x) const
{
  bitset::RankFlags bit(x);
  rootdata::RootSet result = d_baseNoncompact;

  for (size_t j=0; j<adjointFiberRank() ; ++j)
    if (bit[j])
      result ^= d_noncompactShift[j];
  return result;
}

/*!
  \brief Returns the compact imaginary roots for element \#x in the
   adjoint fiber.
*/
rootdata::RootSet Fiber::compactRoots(AdjointFiberElt x) const
{
  rootdata::RootSet result = imaginaryRootSet();
  result.andnot(noncompactRoots(x));
  return result;
}

/*!
  \brief Flags in gr the noncompact simple imaginary roots for element \#x
  in the adjoint fiber.

  Precondition: |x| represents an element of the subquotient in
  |d_adjointFiberGroup|.
*/
gradings::Grading Fiber::grading(AdjointFiberElt x) const
{
  assert(d_gradingShift.size()==adjointFiberRank()); // length of combination

  gradings::Grading gr = d_baseGrading;
  gr ^= bitvector::combination(d_gradingShift,bitset::RankFlags(x));

  return gr;
}

/*!
  \brief Returns an element |x| of the adjoint fiber group such that
  |grading(x)==gr| (if it exists, it is unique).

  Algorithm: the grading must be reached from the base grading (with all ones)
  by adding a linear combination of grading shifts. So we set up and solve the
  system asking for the appropriate linear combination of grading shifts.
*/
AdjointFiberElt Fiber::gradingRep(const gradings::Grading& gr) const
{
  const size_t ir=imaginaryRank();
  latticetypes::SmallBitVector target(gr,ir);
  target -= latticetypes::SmallBitVector(d_baseGrading,ir);

  latticetypes::SmallBitVectorList shifts(fiberRank());
  for (size_t j = 0; j < shifts.size(); ++j)
    shifts[j]=latticetypes::SmallBitVector(d_gradingShift[j],ir);

  bitset::RankFlags result; bool success=firstSolution(result,shifts,target);
  if (not success)
    throw std::runtime_error("Representative of impossible grading requested");

  return AdjointFiberElt(result.to_ulong());
}

/*!
  \brief Returns the matrix of the involution on the weight lattice
  of the Cartan subgroup.

  NOTE: this is not inlined to avoid a dependency upon tori.h in the .h file.
*/
const latticetypes::LatticeMatrix& Fiber::involution() const
{
  return d_torus->involution();
}

/*!
  \brief Returns the fiber group element \f$m_\alpha\f$ corresponding to |cr|.

  Precondition: |cr| a weight vector for an imaginary coroot \f$\alpha^\vee\f$
  for this Cartan.
*/
latticetypes::SmallBitVector Fiber::mAlpha(const rootdata::Root& cr) const
{
  return d_fiberGroup.toBasis(latticetypes::SmallBitVector(cr));
}

size_t Fiber::plusRank() const

/*!
  \brief Returns the dimension of the +1 eigenspace of the involution.

  NOTE: this is not inlined to avoid a dependency upon tori.h in the .h file.
*/

{
  return d_torus->plusRank();
}

size_t Fiber::minusRank() const

/*!
  \brief Returns the dimension of the -1 eigenspace of the involution.

  NOTE: this is not inlined to avoid a dependency upon tori.h in the .h file.
*/

{
  return d_torus->minusRank();
}


/*!
  \brief Returns the image of x in the adjoint fiber group.

  Precondition: x is a valid element in the fiber group.
*/
AdjointFiberElt Fiber::toAdjoint(FiberElt x) const
{
  bitset::RankFlags xf(x);
  latticetypes::SmallBitVector v(xf,fiberRank());

  latticetypes::SmallBitVector w(adjointFiberRank());
  d_toAdjoint.apply(w,v);

  return AdjointFiberElt(w.data().to_ulong());
}


/*!
  \brief Returns the class number in the weak real form partition of the
  strong real form \#c in central square class \#csc.

  The pair (c,rfc) is the software representation of an equivalence
  class of strong real forms (always assumed to induce \f$\tau\f$ on H). The
  integer |csc| labels an element of Z^delta/[(1+delta)Z], thought of as
  a possible square value for strong real forms.  The fiber group acts
  simply transitively on strong real forms with square equal to |csc|.
  The integer c labels an orbit of W_i on this fiber group coset; this
  orbit is the equivalence class of strong real forms.

  This function computes the weak real form (W_i orbit on the adjoint
  fiber group) corresponding to (c,csc).

  First, |csc| also labels a coset of the fiber group image in the
  adjoint fiber group. The coset |csc| is a union of weak real
  forms (W_i orbits on the adjoint fiber group).  The integer brf
  indexes the base W_i orbit on csc, and "by" is the base adjoint
  fiber group element in brf.

  The fiber group coset corresponding to csc is labelled by the fiber
  group itself using a base point with image by.  The integer x is a
  representative in the fiber group of the orbit csc.  Its image in
  the adjoint fiber group is y.  Translating y by the base point "by"
  gives an adjoint fiber group element representing the weak real form
  we want.
*/
adjoint_fiber_orbit Fiber::toWeakReal(fiber_orbit c, square_class csc) const
{
  using namespace partition;

  // get the basepoint for this real form class
  AdjointFiberElt by = class_base(csc);

  // find image of strong real form element in the adjoint fiber group
  AdjointFiberElt y = toAdjoint(d_strongReal[csc].classRep(c));

  // add, and find orbit number
  return d_weakReal(by^y);
}

} // namespace cartanclass


/*****************************************************************************

        Chapter IV -- Functions declared in cartanclass.cpp

******************************************************************************/

namespace cartanclass {

/*!
  \brief Returns the sum of positive compact imaginary roots for |x| in |f|.
*/
latticetypes::Weight
compactTwoRho(AdjointFiberElt x, const Fiber& f, const rootdata::RootDatum& rd)
{
  return rd.twoRho(f.compactRoots(x));
}

/*!
  \brief Returns the restriction of the grading in |rs| to |rl|.
*/
gradings::Grading
restrictGrading(const rootdata::RootSet& rs, const rootdata::RootList& rl)
{
  gradings::Grading result;
  for (size_t j = 0; j < rl.size(); ++j)
    result.set(j,rs.isMember(rl[j]));

  return result;
}

/*!
  \brief Returns a grading in the orbit corresponding to |rf|, with the
  smallest possible number of noncompact roots.

  Precondition: |f| is the fundamental fiber;

  Explanation: for each noncompact noncomplex irreducible real form, there
  is at least one grading with exactly one noncompact simple root. Our choice
  amounts to a grading which induces one of the aforementioned ones on each
  noncompact noncomplex simple factor. The function of this is to enable
  easy type recognition.

  NOTE : the grading is represented in terms of simple roots for the root
  system |rd|. This is OK; knowledge of the gradings of those roots is enough
  to define the real form.
*/
gradings::Grading
specialGrading(const cartanclass::Fiber& f,
	       realform::RealForm rf, const rootdata::RootDatum& rd)


{
  using namespace bitset;
  using namespace cartanclass;
  using namespace gradings;
  using namespace rootdata;

  std::set<Grading,GradingCompare> grs;
  unsigned long n = f.adjointFiberSize();

  // sort the gradings that occur in this class
  for (unsigned long j = 0; j < n; ++j) {
    if (f.weakReal()(j) == rf)
      grs.insert(restrictGrading(f.noncompactRoots(j),rd.simpleRootList()));
  }

  // return the first element
  return *(grs.begin());
}

/*!
  \brief Returns a set of strongly orthogonal roots, leading from
  the fundamental Cartan to the most split one for the strong real form |rf|

  Algorithm: this is quite simple: as long as there are noncompact imaginary
  roots, use one of them to get to a less compact Cartan.
*/
rootdata::RootList
toMostSplit(const cartanclass::Fiber& fundf,
	    realform::RealForm rf, const rootdata::RootDatum& rd)
{
  using namespace cartanclass;
  using namespace rootdata;

  RootSet ir = fundf.imaginaryRootSet();
  unsigned long rep = fundf.weakReal().classRep(rf);
  RootSet rs=fundf.noncompactRoots(rep);

  rootdata::RootList result;
  while (not rs.empty()) {
    RootNbr rn = rs.front();
    for (RootSet::iterator i = ir.begin(); i(); ++i) {
      if (not rd.isOrthogonal(rn,*i)) {
	ir.remove(*i);
	rs.remove(*i);
      } else if (rd.sumIsRoot(rn,*i))
	rs.flip(*i);
    }
    result.push_back(rn);
  }

  strongOrthogonalize(result,rd);

  return result;
}

}

/*****************************************************************************

        Chapter VI -- Auxiliary methods of |CartanClass|

******************************************************************************/

/*! \brief Computes the list of the simple roots for a complex factor in
  \f$W^\tau\f$, where \f$\tau\f$ is the root datum involution of our Cartan
  class.

  Explanation: \f$W^\tau\f$ is the semidirect product of $W^R x W^{iR}$ (Weyl
  groups of the real and imaginary root systems), with the diagonal subgroup
  of $W^C$, where $W^C$ is the Weyl group of the root system $Phi^C$
  orthogonal to both the sums of positive imaginary and real roots. That root
  system is complex for the involution induced by \f$\tau\f$, i.e., it
  decomposes as orthogonal sum of two subsystems interchanged by \f$\tau\f$;
  we return a basis of one "half" of it.

  NOTE: there was a bad bug here in an earlier version, which amounted to the
  implicit assumption that the standard positive root system for the $Phi^C$
  is \f$\tau\f$-stable; this is very false. It would be true for an involution
  of the based root datum (the distinguished involution of the inner class),
  which would stabilise everyting mentioned here, and in fact it is also true
  for any involution that is canonical in its (twisted conjugation) class; in
  general however although $Phi^C$ is \f$\tau\f$-stable (the sum of positive
  imaginary roots is \f$\tau\f$-fixed, and the sum of positive real roots is
  \f$\tau\f$-negated), its subsets of positive and simple roots is not. As a
  consequence the root |rTau| below need not correspond to any vertex of the
  Dynkin diagram |dd|. The component of |dd| to whose root system the various
  |rTau| found for the component |c| belong (the "other half" that we want to
  exclude) can be characterised as the vertices |j| whose simple roots |rb[j]|
  are non-orthogonal to some of those roots |rTau|. Hence we exclude for each
  |rTau| any nodes that are non-orthogonal to it.
*/
rootdata::RootList
CartanClass::makeSimpleComplex(const rootdata::RootDatum& rd) const
{
  latticetypes::Weight tri=rd.twoRho(imaginaryRootSet());
  latticetypes::Weight trr=rd.twoRho(realRootSet());

  // collect roots orthogonal to both sums of positive roots
  rootdata::RootList rl;
  for (size_t j = 0; j < rd.numRoots(); ++j)
    if (rd.isOrthogonal(tri,j) and
	rd.isOrthogonal(trr,j))
      rl.push_back(j);

  // get a (positive) basis for that root system, and its Cartan matrix
  rootdata::RootList rb=rd.simpleBasis(rl);

  latticetypes::LatticeMatrix cm; cartanMatrix(cm,rb,rd);

  dynkin::DynkinDiagram dd(cm);

  bitset::RankFlags b(constants::lMask[dd.rank()]); // all bits set

  rootdata::RootList result;
  while (b.any())
  {
    size_t s = b.firstBit();

    // get component of chosen vertex, and flag them as done (i.e., reset)
    bitset::RankFlags c = dd.component(s); b.andnot(c);

    for (bitset::RankFlags::iterator i = c.begin(); i(); ++i)
    {
      // include roots corresponding to vertices in this component |c|
      result.push_back(rb[*i]);

      /* exclude matching componont by removing vertices of simple roots
         non-orthogonal to the (non-simple) $\tau$-image of |rb[*i]| */
      rootdata::RootNbr rTau = involution_image_of_root(rb[*i]);
      for (bitset::RankFlags::iterator j = b.begin(); j(); ++j)
	if (not rd.isOrthogonal(rb[*j],rTau))
	  b.reset(*j);
    }
  }
  return result;
}

/*!
   \brief Returns the size of the twisted involution orbit for this class.
*/
size::Size CartanClass::makeOrbitSize(const rootdata::RootDatum& rd) const
{
  lietype::LieType lt=dynkin::lieType(rd.cartanMatrix());

  // initialise |result| to size of full Weyl group
  size::Size result; weylsize::weylSize(result,lt);

  // divide by product of imaginary, real and complex sizes
  size::Size ws;

  weylsize::weylSize(ws,dynkin::lieType(rd.cartanMatrix(simpleImaginary())));
  result /= ws;

  weylsize::weylSize(ws,dynkin::lieType(rd.cartanMatrix(simpleReal())));
  result /= ws;

  weylsize::weylSize(ws,dynkin::lieType(rd.cartanMatrix(simpleComplex())));
  result /= ws;

  return result;
}

} // namespace atlas
