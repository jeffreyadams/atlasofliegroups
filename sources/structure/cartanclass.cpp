/*
  This is cartanclass.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-201g Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/
// Implementation for the CartanClass and Fiber classes.

#include "cartanclass.h"

#include <cassert>
#include <map>
#include <set>
#include <utility>
#include <stdexcept>
#include <string>     // used implicitly in throwing errors

#include "tags.h"
#include "bits.h"
#include "size.h"	// used in orbit size computation

#include "dynkin.h"     // for |makeSimpleComplex| and |orbit_size|
#include "gradings.h"   // |GradingCompare| and |max_orth|
#include "lietype.h"	// used in orbit size computation
#include "weylsize.h"	// used in orbit size computation
#include "rootdata.h"	// used in |InvolutionData|
#include "tori.h"       // |tori::dualPi0| used in |makeFiberGroup|
#include "weyl.h"	// used in one |InvolutionData| pseudo-constructor

// extra defs for windows compilation -spc
#ifdef WIN32
#include <iterator>
#endif

namespace atlas {

  namespace cartanclass {

    namespace {

/*
  A function object defining the action of the simple-imaginary reflections on
  a fiber.

  The function object this class provides, intended for use by
  |partition::orbits|, takes the index |s| of a simple-imaginary root and a
  number |x| encoding an element of the fiber; it returns a number |y| similarly
  encoding the image of that element under the action of |s|.

  In fact the number |x| describes (in binary form) the element of the fiber
  group that translates the base point to fiber element in question, and the
  interpretation of |y| is the same. The action is determined by (1) the
  grading of the simple-imaginary roots associated to the chosen base point of
  the fiber, (2) the grading shifts associated with the generators of the
  fiber group, and (3) the vectors |d_mAlpha[s]| by which each of the simple
  imaginary roots $\alpha_s$ translates if it acts non-trivially (this
  is the image in the fiber group of the coroot of $\alpha_s$). The action
  of $\alpha_s$ will translate |x| by |d_mAlpha[s]| if the grading at
  $\alpha_s$ associated to |x| is odd (noncompact). To facilitate the
  determination of that grading, the grading shift information is stored as
  bitsets |d_alpha[s]| for each simple-imaginary root. Since in |Fiber|, the
  grading shifts are organised by generator of the fiber group, the neceesary
  "transposition" of the information is done by the constructor below.
*/

class FiberAction
{
  Grading d_baseGrading; // indexed by |s|
  RankFlagsList d_alpha;   // indexed by |s| first
  RankFlagsList d_mAlpha;  // indexed by |s| first

public:
// constructors and destructors
  FiberAction(const Grading& gr,
	      const GradingList& gs, // size: fiber rank
	      const RankFlagsList& ma) // size: imaginary rank
    : d_baseGrading(gr)
    , d_alpha(ma.size()) // initialise size only here
    , d_mAlpha(ma)
  {
    for (size_t i = 0; i<ma.size(); ++i)
      for (size_t j = 0; j<gs.size(); ++j)
	  d_alpha[i].set(j,gs[j][i]);
  }

// accessors
  bool grading(RankFlags x, unsigned int s) const
    { return d_baseGrading[s] != x.scalarProduct(d_alpha[s]); }
  unsigned long operator() (unsigned int s, unsigned long x) const
  {
    RankFlags b(x); if (grading(b,s))  b ^= d_mAlpha[s];
    return b.to_ulong();
  }
};

} // |namespace|


/*****************************************************************************

        Chapter I -- The CartanClass class

******************************************************************************/


/*
  Construct the Cartan class with involution theta.

  Precondition: theta is an involution of the root datum rd that can be
  obtained from the distinguished involution by multiplication to the left by
  the action of a Weyl group element w, which is called the twisted involution
  associated to theta.
*/
CartanClass::CartanClass(const RootDatum& rd,
			 const RootDatum& dual_rd,
			 const WeightInvolution& theta)
  : d_fiber(rd,theta)
  , d_dualFiber(dual_rd,theta.negative_transposed())
  , d_simpleComplex(makeSimpleComplex(rd))
  , d_orbitSize(orbit_size(rd))
{
  assert(d_dualFiber.simpleReal()==d_fiber.simpleImaginary());
}

/******** copy and assignment ************************************************/


/******** accessors **********************************************************/

/*
  Tell whether this Cartan class is the most split one for weak real form
  corresponding to class |wrf| in |fiber().weakReal()|.

  Algorithm: this is the case iff the grading corresponding to wrf is trivial,
  i.e., all imaginary roots are compact. Since cross actions by compact
  imaginary roots are trivial, the fiber element |x| below is then unique.
*/
bool CartanClass::isMostSplit(adjoint_fiber_orbit wrf) const
{
  auto& f=fiber(); return f.grading(f.wrf_rep(wrf)).none();
}


/*****************************************************************************

        Chapter II -- The Fiber class

******************************************************************************/


// Fiber constructor for involution |theta| of the root datum |rd|.


Fiber::Fiber(const RootDatum& rd, const WeightInvolution& theta)
  : d_involution(theta)
  , d_involutionData(rd,theta)
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
  assert(grading_kernel(rd).dimension()==0);
}

// Copy and assignment

Fiber::Fiber(const Fiber& other)
  : d_involution(other.d_involution)
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

tori::RealTorus Fiber::torus() const { return tori::RealTorus(d_involution); }

/*       Private accessors of |Fiber| used during construction      */


/*
  The group that acts 1-transitively on each central square class in the fiber.

  This is the subquotient $V_+ + V_-/V_+$ as constructed by |tori::dualPi0|,
  but for $-\theta^t$. So the group is isomorphic to the +1 eigenspace of the
  transpose involution $\theta^t$, modulo the image of $\theta^t+1$, in other
  words it is $\Ker(\theta^t-1)/\Im(\theta^t+1)$, and it also isomorphic to
  the component group of the set of $-\theta$ fixed points of the Cartan $H$.

  Note that only the involution |theta| is used, not the root datum.
*/
SmallSubquotient Fiber::makeFiberGroup() const
{
  return tori::dualPi0(involution().negative_transposed());
}

/*
  The matrix of the transformation induced by |-^theta| on the cocharacter
  lattice for the adjoint group, expressed on the simple coweight basis.

  As this basis is dual to the simple root basis, we first transform the
  involution to one on the root basis, and then take the negative transpose.
*/
CoweightInvolution
adjoint_involution(const RootSystem& rs, const InvolutionData& id)
{
  // write involution in root basis
  WeightList b(rs.rank()); // involution images in simple roots

  for (size_t s=0; s<b.size(); ++s)
    b[s] = rs.root_expr(id.root_involution(rs.simpleRootNbr(s)));

  return WeightInvolution(b,rs.rank()).negative_transposed();
}

/*
  Build the group that acts 1-transitively on the adjoint fiber.

  Algorithm: this is the subquotient $V_+ + V_-/V_+$ for the involution
  induced by |-^theta| on the cocharacter lattice (the |adjoint_involution|).
*/
SmallSubquotient Fiber::makeAdjointFiberGroup(const RootSystem& rs) const
{
  SmallSubquotient result =
    tori::dualPi0(adjoint_involution(rs,d_involutionData));

  assert(result.rank()==rs.rank()); // (rank is that of ambient vector space)

  return result;
}

/*
  Makes the stabilizer of the grading in the adjoint fiber group.

  Explanation: each real form parameter defines a grading of the imaginary
  root system, obtained by pairing with the simple-imaginary roots. It should
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
SmallSubspace Fiber::grading_kernel (const RootSystem& rs) const
{
  // define map
  const SmallBitVectorList& baf =
    d_adjointFiberGroup.space().basis();

  // express simple-imaginary roots on (full) simple roots
  WeightList bsi;
  rs.toRootBasis(simpleImaginary().begin(),simpleImaginary().end(),
		 back_inserter(bsi));


  SmallBitVectorList bsi2(bsi); // and reduce mod 2

  SmallBitVectorList b;

  /* set b[j][i] = <e_j,bsi2[i]> where e_j=baf[j'], with |baf[j']| the
     subquotient basis representative number |j| */
  for (RankFlags::iterator
	 jt = d_adjointFiberGroup.support().begin(); jt(); ++jt)
  {
    SmallBitVector v(imaginaryRank());
    for (size_t i = 0; i < imaginaryRank(); ++i)
      v.set(i,bsi2[i].dot(baf[*jt]));
    b.push_back(v);
  }

  BinaryMap q(b,imaginaryRank());

  return SmallSubspace(q.kernel(),adjointFiberRank());
}


/*
  Flag all noncompact imaginary roots in |flagged_roots|, and return the
  all-ones (quasisplit base) grading of simple-imaginary roots

  Algorithm: for the noncompact roots, express each imaginary root in terms
  of the simple-imaginary ones, and flag it if the sum of coordinates is odd.
*/
Grading Fiber::makeBaseGrading
  (RootNbrSet& flagged_roots,const RootSystem& rs) const

{
  // express all imaginary roots in simple-imaginary basis
  int_VectorList ir;
  rs.toRootBasis(imaginaryRootSet().begin(),
		 imaginaryRootSet().end(),
		 back_inserter(ir),simpleImaginary());

  // now flag, among all roots, those imaginary roots with noncompact grading
  flagged_roots.set_capacity(rs.numRoots());

  for (size_t j = 0; j < ir.size(); ++j)
  {
    const int_Vector& v=ir[j];
    int count=0;
    for (size_t i=0; i<v.size(); ++i)
      count+=v[i];			  // add coefficients (for parity)
    flagged_roots.set_mod2
      (imaginaryRootSet().n_th(j),count); // and flag imaginary root
  }
  return Grading(constants::lMask[imaginaryRank()]); // all ones
}

/*
  Compute, for each basis vector of the adjoint fiber group, the
  grading shifts for simple-imaginary roots, and also sets |all_shifts| to
  flag the grading of the full set of imaginary roots.

  Explanation: component |j| of the result contains the grading on the set of
  simple-imaginary roots given by the canonical basis vector |j| of the
  adjoint fiber group. This is not the grading of a real form, but the amount
  by which the grading changes by the action of that basis vector, whence the
  name grading shift. Similarly, |all_shifts[j]| is set to contain the action
  on all imaginary roots. With these data and the "base" grading corresponding
  to the real form parameter numbered $0$, it is easy to compute the grading
  for an arbitrary real form parameter (as represented by the element of the
  adjoint fiber group moving the base element there).

  Algorithm: the basis for the adjoint fiber group has been expressed in the
  reduction modulo 2 of the simple coweight basis. Therefore, in order to
  apply to a root, it is enough to express that root in the simple root basis,
  and do the scalar product mod 2. Note that in contrast to |makeBaseGrading|
  above, and in spite of the reuse of the same names, we express in the (full)
  simple root basis here, not in the simple-imaginary root basis.
*/
GradingList Fiber::makeGradingShifts
  (std::vector<RootNbrSet>& all_shifts,const RootSystem& rs) const
{
  // express imaginary roots in (full) simple root basis
  const RootNbrList irl
    (imaginaryRootSet().begin(),imaginaryRootSet().end());

  int_VectorList ir;
  rs.toRootBasis(irl.begin(),irl.end(),back_inserter(ir));
  SmallBitVectorList ir2(ir); // |ir2.size()==irl.size()|


  // also express simple-imaginary roots in (full) simple root basis
  const RootNbrList& sil = simpleImaginary();
  int_VectorList si;
  rs.toRootBasis(sil.begin(),sil.end(),back_inserter(si));
  SmallBitVectorList si2(si); // reduce vectors mod 2

  // now compute all results
  const SmallBitVectorList& basis = d_adjointFiberGroup.basis();
  GradingList result;

  // traverse basis of the subquotient |d_adjointFiberGroup|
  for (auto it = basis.begin(); it!=basis.end(); ++it)
  {
    SmallBitVector gen = *it;

    // all imaginary roots part
    RootNbrSet rset(rs.numRoots()); // used only at indices listed in |irl|
    for (size_t j = 0; j < ir2.size(); ++j)
      rset.set_to(irl[j],gen.dot(ir2[j]));
    all_shifts.push_back(rset);

    // simple-imaginary roots part
    Grading gr;
    for (size_t j = 0; j < si2.size(); ++j)
      gr.set(j,gen.dot(si2[j]));
    result.push_back(gr);
  }

  assert(result.size()==adjointFiberRank());

  return result;
}

/*
  Constructs the $m_\alpha$s (images of coroots) in the fiber group,
  for $\alpha$ simple-imaginary.

  The effective number of bits of $m_\alpha$ is |d_fiberGroup.dimension()|

  We take the coroots corresponding to the imaginary simple roots, which in
  the root datum are already expressed in the basis of the coweight lattice
  $X^*$ dual to the weight lattice; we reduce the coordinates modulo 2 (this
  is hidden in the call of |toBasis|, which converts its argument to a
  |SmallBitVector| first), and interpret the result in the subquotient
  |d_fiberGroup| of $X^* / 2X^*$.
*/
RankFlagsList Fiber::mAlphas (const RootDatum& rd) const
{
  RankFlagsList result(imaginaryRank());

  for (size_t i = 0; i<result.size(); ++i)
    result[i]=
      d_fiberGroup.toBasis(SmallBitVector(rd.coroot(simpleImaginary(i))))
      .data();

  return result;
}

/*
  Construct the $m_\alpha$s (images of coroots) in the adjoint fiber group,
  for $\alpha$ simple-imaginary.

  The number of bits of each $m_\alpha$ is
  |d_adjointFiberGroup.dimension()|

  Algorithm: the cocharacter lattice for the adjoint group is spanned by the
  simple coweights. To get the coordinates of an element in that basis (which
  is dual to that of the simple roots), it is enough to pair it with the
  simple roots. The resulting element for the coroot of a simple-imaginary
  $\alpha$ is automatically |theta|-invariant, since $\alpha$ is,
  so its reduction modulo 2 lies in $V_+$ (for the cocharacter lattice). Then
  what is left to do is to convert to the basis of the adjoint fiber group,
  which amounts to reducing modulo $V_-$.
*/
RankFlagsList Fiber::adjointMAlphas (const RootSystem& rs) const
{
  RankFlagsList result(imaginaryRank());

  for (size_t i = 0; i<result.size(); ++i)
  {
    SmallBitVector v(rs.rank());
    // compute pairing with simple roots modulo 2
    for(size_t j = 0; j < v.size(); ++j)
    {
      LatticeCoeff c =
	rs.bracket(rs.simpleRootNbr(j),simpleImaginary(i));
      v.set_mod2(j,c);
    }

    // reduce |v| modulo $V_-$
    result[i]=d_adjointFiberGroup.toBasis(v).data();
  }
  return result;
}

/*
  Compute the toAdjoint matrix.

  Explanation: this is the matrix of the natural map from the fiber group to
  the adjoint fiber group. Such a map exists because (1) the root lattice is a
  sublattice of the weight lattice, so there is a restriction map from the
  cocharacter lattice $X_*$ (dual to the character lattice $X^*$) to the
  lattice $Y_*$ spanned by the fundamental coweights (dual to the root
  lattice), and (2) this map induces a map on the respective subquotients of
  the form $(V_+ + V_-)/V_-$, where $V_+$ is the modulo $2X_*$ image of the
  $\theta^t$ fixed cocharacters, and similarly for $V_-$. While the
  restriction map $X_* \to Y_*$ is injective (only) in the semisimple case,
  little can be said about the induced map.

  The above argument shows that the vector |v| computed below of scalar
  products of a generator of the fiber group with simple roots can be validly
  incorporated (by calling |d_adjointFiberGroup.toBasis|) into the adjoint
  fiber group; notably, it represents an element of $V_+ + V_-$ in
  $Y_* / 2Y_*$. (without space that formula would have ended our comment!)
*/
BinaryMap Fiber::makeFiberMap(const RootDatum& rd) const
{

  SmallBitVectorList b_sr
    (rd.beginSimpleRoot(),rd.endSimpleRoot()); // mod 2

  // images of fiber group basis in adjoint fiber group
  SmallBitVectorList b_ad;
  b_ad.reserve(d_fiberGroup.dimension());

  const SmallBitVectorList& b = d_fiberGroup.space().basis();
  const RankFlags& supp = d_fiberGroup.support();
  size_t n = rd.semisimple_rank();
  for (RankFlags::iterator it = supp.begin(); it(); ++it)
  {
    SmallBitVector v(n);
    for(size_t i = 0; i < n; ++i)
      v.set(i, b_sr[i].dot(b[*it]) );

    b_ad.push_back(d_adjointFiberGroup.toBasis(v));
  }

  return BinaryMap // convert vectors to a matrix
    (b_ad,d_adjointFiberGroup.dimension());
}


/*
  Compute the partition of the adjoint fiber, whose parts correspond to the
  weak real forms.

  Algorithm: we construct the |FiberAction| object corresponding to the action
  of the imaginary Weyl group on the adjoint fiber, and then call the |orbits|
  method to make the partition. For the fiber action we can use the base
  grading and the grading shifts "as is", while the fiber group elements
  $m_\alpha$ are computed by |adjointMAlphas|. */
Partition Fiber::makeWeakReal(const RootSystem& rs) const
{
  RankFlagsList ma=adjointMAlphas(rs);
  return partition::orbits(FiberAction(d_baseGrading,d_gradingShift,ma),
			   imaginaryRank(),adjointFiberSize());
}



/*
  Compute the partition of the weak real forms according to square classes.

  Explanation: while weak real forms correspond to orbits in the adjoint
  fiber, they can be grouped into even coarser "central square classes",
  defined by the condition that the corresponding orbits belong to the same
  coset in the adjoint fiber (group) by the image (under |toAdjoint|) of the
  fiber group. The imaginary Weyl group $W_{im}$ acts on both the fiber group
  and the adjoint fiber group, and |toAdjoint| intertwines these actions, so
  that an orbit defining a weak real form is always contained in a single
  coset of the image of |toAdjoint|; therefore we really have a partition of
  the set of weak real forms. Its parts are called central square classes,
  since for these weak real forms any $x=g.\delta\in G.\delta$ whose
  $H$-conjugacy class defines a fiber element in a strong real form lying over
  the weak real form gives the same value of $x^2\in Z(G)$ modulo
  $(1+\delta)(Z(G))$. From the above it follows there are $2^m$ central
  square classes where $m=adjointFiberRank-rank(toAdjoint)$.

  Algorithm: we compute the quotient of the adjoint fiber group by the fiber
  group as a |Subquotient| object whose |space()| is the whole adjoint fiber
  group. Its |toBasis| method will map each weak real form to a vector expressed
  on the basis of the (sub)quotient whose |to_ulong()| value will be used to
  characterise the central square class. Since we want to preserve this
  numbering when constructing the resulting partition (rather than reordering by
  smallest element), we call the constructor with a |tags::UnnormalizedTag|
  argument.
*/
Partition Fiber::makeRealFormPartition() const
{
  std::vector<unsigned long> cl(numRealForms());

  SmallBitVectorList b_id;
  bitvector::initBasis(b_id,d_adjointFiberGroup.dimension());

  /* construct |modFiberImage| as quotient of spans of |b_id| (i.e., all) and
     image of |d_toAdjoint|; it is a quotient of the adjoint fiber group! */
  SmallSubquotient modFiberImage
    (b_id,d_toAdjoint.image(),d_adjointFiberGroup.dimension());

  // take representatives of weak real forms and reduce modulo fiber image
  for (size_t i = 0; i < cl.size(); ++i)
  {
    unsigned long y = d_weakReal.classRep(i);
    SmallBitVector v
      (RankFlags(y),d_adjointFiberGroup.dimension());

    // reduce modulo image of map from fiber group to adjoint fiber group
    cl[i] = modFiberImage.toBasis(v).data().to_ulong();
  }

  /* partition the set $[0,numRealForms()[$ according to the
     |modFiberImage.size()| distinct values in the image of |cl|
   */
  Partition result(cl,tags::UnnormalizedTag());
  assert(result.classCount()==modFiberImage.size());
  return result;
}


/*
  Compute the strong real form partitions.

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
std::vector<Partition> Fiber::makeStrongReal (const RootDatum& rd) const
{
  /* get the grading shifts; these are obtained from the images of the
     canonical basis vectors of the fiber group in the adjoint fiber group. */

  GradingList gs(fiberRank());

  // basically compute |L*d_toAdjoint| where rows of L are |d_gradingShift|
  for (size_t i = 0; i < gs.size(); ++i)
    gs[i]=bitvector::combination(d_gradingShift,d_toAdjoint.column(i).data());

  // get the $m_\alpha$s
  RankFlagsList ma = mAlphas(rd);

  // make the various partitions

  std::vector<Partition> result;
  result.reserve(d_realFormPartition.classCount());

  size_t order = fiberSize(); // order of fiber group
  for (square_class i=0; i<d_realFormPartition.classCount(); ++i)
  {
    Grading bg = grading(class_base(i));
    result.push_back
      (partition::orbits(FiberAction(bg,gs,ma),imaginaryRank(),order));
  }

  return result;
}

/*
  For each weak real form |wrf| we can choose a fiber element |x| (in the
  affine space corresponding to its central square class |c|) that maps to the
  chosen adjoint fiber element representative of |wrf|. The auxiliary method
  |makeStrongRepresentatives| makes a vector of size |numRealForms())| (to be
  stored in |d_strongRealFormReps|) whose element |wrf| is the pair $(x,c)$.
*/
std::vector<StrongRealFormRep> Fiber::makeStrongRepresentatives() const
{
  SmallBitVectorList b; b.reserve(fiberRank());

  for (size_t j = 0; j < fiberRank(); ++j)
    b.push_back(d_toAdjoint.column(j));

  std::vector<StrongRealFormRep> result(numRealForms());

  for (size_t wrf = 0; wrf<result.size(); ++wrf)
  {
    square_class c = central_square_class(wrf);

    // find representative |yf| of |wrf| in the adjoint fiber (group)
    RankFlags yf(d_weakReal.classRep(wrf));

    // subtract base point
    yf ^= class_base(c).data();

    // find preimage |xf| of |yf| in the fiber
    SmallBitVector v(yf,adjointFiberRank()); // the desired image

    // solve equation |toAdjoint(xf)=v|
    RankFlags xf;
    bool success=bitvector::combination_exists(b,v,xf);
    assert(success);  // there has to be a solution!
    ndebug_use(success); // avoid warning about unused variable

    // make representative
    fiber_orbit x = d_strongReal[c].class_of(xf.to_ulong());
    result[wrf] = std::make_pair(x,c);
  }

  return result;
}






/*
		     Public accessors of the |Fiber| class
*/



/*
  Return the set of noncompact imaginary roots for elt |x| in the adjoint
  fiber (whose bitset is interpreted as an element of the subquotient in
  |d_adjointFiberGroup|).
*/
RootNbrSet Fiber::noncompactRoots(AdjointFiberElt x) const
{
  RootNbrSet result = d_baseNoncompact;
  for (auto it=x.data().begin(); it(); ++it)
    result ^= d_noncompactShift[*it];
  return result;
}

// Return the compact imaginary roots for element |x| in the adjoint fiber.
RootNbrSet Fiber::compactRoots(AdjointFiberElt x) const
{
  RootNbrSet imaginary = imaginaryRootSet(); // take a copy
  return imaginary.andnot(noncompactRoots(x));
}

// Flag the noncompact simple-imaginary roots for adjoint fiber element |x|.
Grading Fiber::grading(AdjointFiberElt x) const
{
  assert(x.size()==d_gradingShift.size()); // length of combination

  Grading gr = d_baseGrading;
  gr ^= bitvector::combination(d_gradingShift,x.data());

  return gr;
}

/*
  Return an element |x| of the adjoint fiber group such that |grading(x)==gr|
  (if it exists, it is unique).

  Algorithm: the grading must be reached from the base grading (with all ones)
  by adding a linear combination of grading shifts. So we set up and solve the
  system asking for the appropriate linear combination of grading shifts.
*/
AdjointFiberElt Fiber::gradingRep(const Grading& gr) const
{
  const size_t ir=imaginaryRank();
  SmallBitVector target(gr,ir);
  target -= SmallBitVector(d_baseGrading,ir); // this complements gradings

  SmallBitVectorList shifts; shifts.reserve(adjointFiberRank());
  for (size_t i = 0; i < adjointFiberRank(); ++i)
    shifts.push_back(SmallBitVector(d_gradingShift[i],ir));

  RankFlags result;
  if (combination_exists(shifts,target,result)) // sets |result|
    return AdjointFiberElt(result,adjointFiberRank());
  else
    throw std::runtime_error("Representative of impossible grading requested");
}

/*
  Return the fiber group element $m_\alpha$ corresponding to |cr|.

  Precondition: |cr| a weight vector for an imaginary coroot $\alpha^\vee$
  for this Cartan.
*/
SmallBitVector Fiber::mAlpha(const rootdata::Root& cr) const
{
  return d_fiberGroup.toBasis(SmallBitVector(cr));
}

/*
  Return the image of |x| in the adjoint fiber group.

  Precondition: |x| is a valid element in the fiber group.
*/
AdjointFiberElt Fiber::toAdjoint(FiberElt x) const
{
  return d_toAdjoint*x;
}


/*
  Return the class number in the weak real form partition of the strong real
  form |c| in central square class |csc|.

  In spite of the name of this method, the value returned is not the number
  globally associated to the weak real form, which is called a real form label
  but which |Fiber| knows nothing about; rather the value is the number of the
  $W_{im}$-orbit in the adjoint fiber (group). However, if the fiber happens
  to be the fundamental fiber, then this value is equal to its real form label.

  The pair |(c,csc)| is the software representation of a strong real form
  (always assumed to induce $\tau$ on $H$). The integer |csc| labels an
  element of $Z^\delta/[(1+\delta)Z]$, thought of as a possible square value
  for strong involutions. The fiber group acts simply transitively on strong
  involutions with square equal to |csc|, so they can be viewed as a coset for
  the fiber group; the imaginary Weyl group $W_{im}$ acts on this coset in a
  manner that depends on |csc|. The integer |c| labels an orbit for this
  action, which is the way we identify a strong real form.

  This function computes the weak real form ($W_{im}$ orbit on the
  \emph{adjoint} fiber group) corresponding to |(c,csc)|.

  First, |csc| also labels a coset of the fiber group image in the adjoint
  fiber group, and |class_base(csc)| is a coset representative (base point).

  The map |toAdjoint| defines a $\Z/2\Z$-linear map from the fiber group to
  the adjoint fiber group, from which the $\Z/2\Z$-affine map for this square
  class is obtained by adding the base point (as a bitvector). It remains to
  (upstream) find a fiber group element in the |W_{im}|-orbit labelled by |c|,
  and (downstream) to find the number of the weak real form of its image under
  the affine map; the former is |d_strongReal[csc].classRep(c)| and the latter
  is obtained by the method |class_of| of the partition |d_weakReal|.
*/
adjoint_fiber_orbit Fiber::toWeakReal(fiber_orbit c, square_class csc) const
{
  FiberElt rep (RankFlags(fiber_partition(csc).classRep(c)),fiberRank());
  // get orbit number of image of strong real form element in the adjoint fiber
  return adjoint_orbit(class_base(csc) + toAdjoint(rep));
}


/*****************************************************************************

        Chapter IV -- Functions declared in cartanclass.cpp

******************************************************************************/

// Return the sum of positive compact imaginary roots for |x| in |f|.
Weight compactTwoRho(AdjointFiberElt x, const Fiber& f, const RootDatum& rd)
{
  return rd.twoRho(f.compactRoots(x));
}

// Returns the restriction of the grading in |rs| to |rl|.
Grading restrictGrading(const RootNbrSet& rs, const RootNbrList& rl)
{
  Grading result;
  for (size_t i=0; i<rl.size(); ++i)
    result.set(i,rs.isMember(rl[i]));

  return result;
}

/*
  Find a grading in the adjoint fiber orbit corresponding to |rf| with among
  such gradings the smallest possible (0 or 1 per simple factor) of simple roots
  marked as noncompact. The grading is represented on the set of all simple
  roots, even though the simple roots that are not imaginary are not really
  being graded, and therefore not candidate for being noncompact.

  Precondition: |f| is the fundamental fiber.

  With the exception of unequal rank A_n, each non-quasicompact non-complex real
  form on a simple type admits at least one grading with exactly one noncompact
  simple root. The position of this noncompact root among the simple roots
  (taking care to consider them in the Bourbaki ordering for the Dynkin diagram)
  though not always unique, easily allows recognising the real form, which will
  ultimately be done in the function |output::print_real_form_name|. In detail,
  for the cases relevant here, skipping the quasicompact forms, and counting
  simple root positions from 0 (for once this turns out slightly awkward):

  A_1: sl(2,R) is unique non-compact form (the simple root is noncompact);
      other non equal rank A_n cases:
       su(p,n+1-p) a single noncompact simple root can be at position p-1 or n-p;
  B_n: so(2p,2n+1-2p) a single noncompact simple root can be at position p-1;
       (forcing the even index first makes this zig-zag correspondence seem easy)
  C_n: sp(p,n-p) a single noncompact simple root can be at position p-1 or n-p-1
       sp(2n,R) a single noncompact simple root can be at position n-1
  D_n equal rank:
       so(2,2n-2) a single noncompact simple root can be at position 0 only
       so(2p,2n-2p), p>1: single noncompact simple root at position p-1 or n-p-1
       n odd:
        so*(2n) a single noncompact simple root can be at position n-2 or n-1
       n even:
        so*(2n)[1,0] single noncompact simple root can be at position n-1-(n/2)%2
        so*(2n)[0,1] single noncompact simple root can be at position n-2+(n/2)%2
      unequal rank:
       so(2p+1,2n-2p-1) single noncompact simple root at position p-1 or n-p-2
  G2: g2(R) is the unique non-compact form (single noncompacts at position 0 or 1)
  F4: f4(so(9)) a single noncompact simple short root (at position 2 or 3)
      f4(R) a single noncompact simple long root (at position 0 or 1)
  E6 equal rank:
      e6(so(10).u(1)) single noncompact simple root can be at position 0 or 5
      e6(su(6).su(2)) single noncompact simple root can be at position 1,2,3,4
     unequal rank:
      e6(R) is the unique non-quasicompact form (single noncompacts at 1 or 3)
  E7: e7(e6.u(1)) a single noncompact simple root can be at position 6
      e7(so(12).su(2)) single noncompact simple root can be at position 0,2,3,5
      e7(R)  a single noncompact simple root can be at position 1 or 4
  E8: e8(e7.su(2))  a single noncompact simple root can be at position 2,3,6,7
      e8(R)  a single noncompact simple root can be at position 0,1,4,5,

  (The unequal rank inner class in type $A_n$ is particular by the rareness of
  imaginary simple roots, candidate for being noncompact: there is at most one,
  and only if $n$ is odd; for this case this still just suffices to distinguish
  sl(n+1,R) where that root is noncompact from sl(n+1/2,H), where it is compact.)

  Our choice of grading for non-simple types will be a grading that has one
  noncompact simple root for each non-quasicompact non-complex simple factor.
  This can be achieved by minimising the number of noncompact simple roots.

  We do not do any effort to pick a particular minimal number grading if there
  should be more than one, as that would involve comparing the position relative
  to the standardised ordering of the Dynkin diagram to make any sense, which is
  not worth the hassle: the remaining ambiguity is dealt with at printing time.
*/
Grading specialGrading(const Fiber& f, RealFormNbr rf, const RootSystem& rs)
{
  unsigned minimum_so_far=rs.rank()+1; // impossibly high value
  Grading result; // should be assigned to at least once below

  unsigned long n = f.adjointFiberSize();

  // sort the gradings that occur in this class
  for (unsigned long i=0; i<n; ++i)
  { AdjointFiberElt x(RankFlags(i),f.adjointFiberRank());
    if (f.adjoint_orbit(x) == rf) // only consider the orbit for |rf|
    { // and only consider the grading of those imaginary roots that are simple
      Grading gr = restrictGrading(f.noncompactRoots(x),rs.simpleRootList());
      if (gr.count()<minimum_so_far) // only minimize the number of noncompacts
      {
	result = gr; // best candidate seen yet
	minimum_so_far = gr.count();
      }
    }
  }
  assert(minimum_so_far<=rs.rank()); // we advanced at least once
  return result;
}

Grading specialGrading
  (const Partition& fg_partition, RealFormNbr rf, RankFlags imaginary_simples)
{
  auto n = fg_partition.size();
  assert(bits::bitCount(n)==1); // we must have a power of $2$ as size

  auto i=fg_partition.classRep(rf); // get first representative element
  RankFlags compacts (i);
  auto max = compacts.count();
  for (++i; i<n; ++i)
    if (fg_partition.class_of(i) == rf and bits::bitCount(i)>=max)
    {
      compacts = RankFlags(i);
      max=compacts.count();
    }

  // convert |compacts| to a grading represented on the simple roots
  compacts^=RankFlags(n-1); // take complement; no need to compute $\log_2(n)$
  return compacts.unslice(imaginary_simples);
}

/*
  a set of strongly orthogonal roots, leading from the fundamental Cartan to
  the most split one for the real form |rf|

  Algorithm: a greedy one: as long as there are noncompact imaginary
  roots, use one of them to get to a less compact Cartan.
*/
RootNbrSet
toMostSplit(const Fiber& fundf, RealFormNbr rf, const RootSystem& rs)
{
  RootNbrSet ir = fundf.imaginaryRootSet();
  AdjointFiberElt rep = fundf.wrf_rep(rf);

  return gradings::max_orth(fundf.noncompactRoots(rep),ir,rs);
}


/*****************************************************************************

        Chapter VI -- Auxiliary methods of |CartanClass|

******************************************************************************/


/*
  The list of the simple roots for a complex factor in $W^\theta$, where
  |theta| is the root datum involution of our Cartan class.

  Explanation: $W^\theta$ is the semidirect product of $W^R x W^{iR}$ (Weyl
  groups of the real and imaginary root systems), with the diagonal subgroup
  of $W^C$, where $W^C$ is the Weyl group of the root system $Phi^C$
  orthogonal to both the sums of positive imaginary and real roots. That root
  system is complex for the involution induced by |theta|, i.e., it
  decomposes as orthogonal sum of two subsystems interchanged by |theta|;
  we return a basis of one "half" of it.

  NOTE: there was a bad bug here in an earlier version, which amounted to the
  implicit assumption that the standard positive root system for the $Phi^C$
  is |theta|-stable; this is very false. It would be true for an involution
  of the based root datum (the distinguished involution of the inner class),
  which would stabilise everything mentioned here, and in fact it is also true
  for any involution that is canonical in its (twisted conjugation) class; in
  general however although $Phi^C$ is |theta|-stable (the sum of positive
  imaginary roots is |theta|-fixed, and the sum of positive real roots is
  |-theta|-fixed), its subsets of positive and simple roots are not. As a
  consequence the root |rTau| below need not correspond to any vertex of the
  Dynkin diagram |dd|. The component of |dd| to whose root system the various
  |rTau| found for the component |c| belong (the "other half" that we want to
  exclude) can be characterised as the vertices |i| whose simple roots |rb[i]|
  are non-orthogonal to some of those roots |rTau|. Hence we exclude for each
  |rTau| any nodes that are non-orthogonal to it.
*/
RootNbrList CartanClass::makeSimpleComplex(const RootDatum& rd) const
{
  Weight tri=rd.twoRho(imaginaryRootSet());
  Weight trr=rd.twoRho(realRootSet());

  // collect roots orthogonal to both sums of positive roots
  RootNbrSet rs(rd.numRoots());
  for (size_t i = 0; i < rd.numRoots(); ++i)
    if (rd.is_orthogonal(tri,i) and
	rd.is_orthogonal(trr,i))
      rs.insert(i);

  // get a (positive) basis for that root system, and its Cartan matrix
  RootNbrList rb=rd.simpleBasis(rs);

  int_Matrix cm =rd.cartanMatrix(rb);

  dynkin::DynkinDiagram dd(cm);

  RankFlags b(constants::lMask[dd.rank()]); // all bits set

  RootNbrList result;
  while (b.any())
  {
    size_t s = b.firstBit();

    // get component of chosen vertex, and flag them as done (i.e., reset)
    RankFlags c = dd.component(s); b.andnot(c);

    for (RankFlags::iterator it = c.begin(); it(); ++it)
    {
      // include roots corresponding to vertices in this component |c|
      result.push_back(rb[*it]);

      /* exclude matching componont by removing vertices of simple roots
         non-orthogonal to the (non-simple) |theta|-image of |rb[*it]| */
      RootNbr rTau = involution_image_of_root(rb[*it]);
      for (RankFlags::iterator jt = b.begin(); jt(); ++jt)
	if (not rd.is_orthogonal(rb[*jt],rTau))
	  b.reset(*jt);
    }
  }
  return result;
}

// The size of the twisted involution orbit for this class.
size_t CartanClass::orbit_size(const RootSystem& rs) const
{
  LieType lt=dynkin::Lie_type(rs.cartanMatrix());

  // initialise |result| to size of full Weyl group
  size::Size result = weylsize::weylSize(lt);

  // divide by product of imaginary, real and complex sizes
  result /=
    weylsize::weylSize(rs.subsystem_type(simpleImaginary()));

  result /=
    weylsize::weylSize(rs.subsystem_type(simpleReal()));

  result /=
    weylsize::weylSize(rs.subsystem_type(simpleComplex()));

  return result.toUlong();
}

}  // |namespace cartanclass|

} // |namespace atlas|
