/*!
\file
\brief Class definitions and function declarations for the RootDatum class.
*/
/*
  This is rootdata.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef ROOTDATA_H  /* guard against multiple inclusions */
#define ROOTDATA_H

#include "rootdata_fwd.h"

#include <algorithm>

#include "prerootdata_fwd.h"
#include "weyl_fwd.h"

#include "arithmetic.h"
#include "bitset.h"
#include "bitmap.h"
#include "latticetypes.h"
#include "lietype.h"
#include "dynkin.h"
#include "setutils.h"
#include "tags.h"

namespace atlas {

namespace rootdata {

  namespace LT = latticetypes;

}

/******** type declarations *************************************************/


/******** function declarations *********************************************/

namespace rootdata {

LT::LatticeMatrix dualBasedInvolution
  (const LT::LatticeMatrix&, const RootDatum&);

RootSet makeOrthogonal(const RootSet& o, const RootSet& subsys,
		       const RootSystem& rs);

void toDistinguished(LT::LatticeMatrix&, const RootDatum&);

LT::LatticeMatrix refl_prod(const RootSet&, const RootDatum&);

weyl::WeylWord toPositive(LT::Weight, const RootDatum&);

RootDatum integrality_datum(const RootDatum& rd,
			    const LT::RatLatticeElt& lambda);

arithmetic::RationalList integrality_points(const RootDatum& rd,
					    LT::RatWeight& nu);

} // namespace rootdata

/******** type definitions **************************************************/

namespace rootdata {


class RootSystem
{
  typedef signed char byte;
  typedef matrix::Vector<byte> Byte_vector;
  struct root_info
  {
    Byte_vector root, dual; // root in root basis, coroot in coroot basis
    bitset::RankFlags descents,ascents; // for reflections by simple roots

    root_info(const Byte_vector& v)
    : root(v), dual(), descents(), ascents() {}
  };
  struct root_compare; // necessary for sorting roots

  size_t rk; // rank of root system

  Byte_vector Cmat;

  std::vector<root_info> ri; //!< List of information about positive roots

//!\brief Root permutations induced by reflections in positive roots.
  std::vector<setutils::Permutation> root_perm;

  // internal access methods
  byte& Cartan_entry(size_t i, size_t j) { return Cmat[i*rk+j]; }
  const byte& Cartan_entry(size_t i, size_t j) const { return Cmat[i*rk+j]; }
  Byte_vector& root(size_t i) { return ri[i].root;}
  Byte_vector& coroot(size_t i) { return ri[i].dual;}
  const Byte_vector& root(size_t i) const { return ri[i].root;}
  const Byte_vector& coroot(size_t i) const { return ri[i].dual;}
  RootNbr abs(RootNbr alpha) const // offset of corresponding positive root
  { return isPosRoot(alpha) ? alpha-numPosRoots() : numPosRoots()-1-alpha; }

  void cons(const latticetypes::LatticeMatrix& Cartan_matrix); // panse bete
 public:

// constructors and destructors

  explicit RootSystem(const latticetypes::LatticeMatrix& Cartan_matrix);

  RootSystem(const RootSystem& rs, tags::DualTag);

// accessors

  size_t rank() const { return rk; } // semisimple rank when part of |RootDatum|
  unsigned long numPosRoots() const { return ri.size(); }
  unsigned long numRoots() const { return 2*numPosRoots(); }

  // Cartan matrix by entry and as a whole
  latticetypes::LatticeCoeff cartan(size_t i, size_t j) const
  { return Cartan_entry(i,j); };
  latticetypes::LatticeMatrix cartanMatrix() const;
  lietype::LieType Lie_type() const
  { return dynkin::Lie_type(cartanMatrix()); }

  // for subsystem
  latticetypes::LatticeMatrix cartanMatrix(const RootList& sub) const;
  lietype::LieType Lie_type(RootList sub) const
  { return dynkin::Lie_type(cartanMatrix(sub)); }



// root list access

  latticetypes::LatticeElt root_expr(RootNbr alpha) const;  // in simple roots
  latticetypes::LatticeElt coroot_expr(RootNbr alpha) const;// in simple coroots

  // convert sequence of root numbers to expressions in the simple roots
  template <typename I, typename O>
    void toRootBasis(const I&, const I&, O) const;
  // convert sequence of root numbers to expressions in subsystem simple weights
  template <typename I, typename O>
    void toSimpleWeights(const I&, const I&, O, const RootList&) const;
  // convert sequence of root numbers to expressions in subsystem simple roots
  template <typename I, typename O>
    void toRootBasis(const I&, const I&, O, const RootList&) const;

  bool isSimpleRoot(RootNbr alpha) const
  { return alpha-numPosRoots()<rk; } // this uses that |RootNbr| is unsigned

  bool isPosRoot(RootNbr alpha) const
  { return alpha>=numPosRoots(); } // second half

  RootNbr simpleRootNbr(size_t i) const
  { assert(i<rk);  return numPosRoots()+i; }

  RootNbr posRootNbr(size_t alpha) const
  { assert(alpha<numPosRoots()); return numPosRoots()+alpha; }

  RootNbr rootMinus(RootNbr alpha) const // roots are ordered symmetrically
  { return numRoots()-1-alpha; }


  RootSet simpleRootSet() const // NOT for iteration over it, seems never used
  {
    RootSet simple_roots(numRoots());
    simple_roots.fill(numPosRoots(),numPosRoots()+rk);
    return simple_roots;
  }

  RootList simpleRootList() const // NOT for iteration over it
  {
    RootList simple_roots(rk);
    for (size_t i=0; i<rk; ++i)
      simple_roots[i]=simpleRootNbr(i);
    return simple_roots;
  }

  RootSet posRootSet() const // NOT for iteration over it, may serve as mask
  {
    RootSet pos_roots(numRoots());
    pos_roots.fill(numPosRoots(),numRoots());
    return pos_roots;
  }

// other accessors

  // the next method only works for _simple_ roots! (whence no RootNbr for |i|)
  const setutils::Permutation& simple_root_permutation(size_t i) const
    { return root_perm[i]; }

  bitset::RankFlags descent_set(RootNbr alpha) const
  {
    RootNbr a = abs(alpha);
    return isPosRoot(alpha) ? ri[a].descents : ri[a].ascents;
  }
  bitset::RankFlags ascent_set(RootNbr alpha) const
  {
    RootNbr a = abs(alpha);
    return isPosRoot(alpha) ? ri[a].ascents : ri[a].descents;
  }

  size_t find_descent(RootNbr alpha) const
  { return descent_set(alpha).firstBit(); }

  bool is_descent(size_t i, RootNbr alpha) const
  { return root_perm[i][alpha]<alpha; } // easier than using |descent_set|

  bool is_ascent(size_t i, RootNbr alpha) const
  { return root_perm[i][alpha]>alpha; }

  RootNbr simple_reflected_root(RootNbr r, size_t i) const
  { return simple_root_permutation(i)[r]; }

  void simple_reflect_root(RootNbr& r, size_t i) const
  { r=simple_root_permutation(i)[r]; }

  RootNbr permuted_root(const weyl::WeylWord& ww, RootNbr r) const
  {
    for (size_t i=ww.size(); i-->0; )
      simple_reflect_root(r,ww[i]);
    return r;
  }


  // for arbitrary roots, reduce root number to positive root offset first
  const setutils::Permutation& root_permutation(RootNbr alpha) const
  { return root_perm[abs(alpha)]; }

  bool isOrthogonal(RootNbr alpha, RootNbr beta) const
  { return root_permutation(alpha)[beta]==beta; }

  // pairing between root |alpha| and coroot |beta|
  latticetypes::LatticeCoeff bracket(RootNbr alpha, RootNbr beta) const;



  // find permutation of roots induced by diagram automorphism
  setutils::Permutation root_permutation(const setutils::Permutation& ) const;

  // extend root datum automorphism given on simple roots to all roots
  setutils::Permutation extend_to_roots(const RootList&) const;


  weyl::WeylWord reflectionWord(RootNbr r) const;

  RootList simpleBasis(RootSet rs) const;

  bool sumIsRoot(RootNbr alpha, RootNbr beta, RootNbr& gamma) const;
  bool sumIsRoot(RootNbr alpha, RootNbr beta) const
  { RootNbr dummy; return sumIsRoot(alpha,beta,dummy); }

  RootSet long_orthogonalize(const RootSet& rest) const;

  RootList high_roots() const;

// manipulators


}; // |class RootSystem|

/*!
\brief Based root datum for a complex reductive group.

What we call a root datum in this program is what is usually called a based
root datum, in other words a fixed choiceof positive roots is always assumed.

The root datum defines the complex reductive group entirely. It consists of a
|RootSystem| that describes the roots and coroots in the lattices they span
themselves, and of additional data that correspond to embeddings of these
lattices into mutually dual free abelian groups (weight and coweight lattices).

The rank |d_rank| is that of the weight and coweight lattices, the root system
itself has rank |semisimpleRank()| which may be smaller. The roots and coroots
are stored in compact form in the |RootSystem|, and again as represented
ointhe weight and coweight lattices, for efficiency of retrieval under this
form. Also constructed are various useful auxiliary things, like d_twoRho (the
sum of the positive roots).

The code is designed to make it preferable always to refer to a root by its
number (index in the root system), for which we use the type name |RootNbr|.
The type |RootList| should in fact (per Fokko) have been called |RootNbrList|.
*/
class RootDatum
: public RootSystem
{

 private:
 /*!
\brief Names describing the  bits of the bitset d_status.

The last enum numFlags is there as a standard programming trick. Its
value - in this case 2 - is one-past-the-last meaningful bit, for
use by accessors.
*/
  enum StatusFlagNames { IsAdjoint, IsSimplyConnected, numFlags };


  typedef bitset::BitSet<numFlags> Status;

/*!
\brief  Rank of the root datum.
*/
  size_t d_rank;

  LT::WeightList d_roots; //!< Full list of roots.
  LT::WeightList d_coroots; //!< Full list of coroots.
  LT::WeightList weight_numer; //!< Simple weight numerators.
  LT::WeightList coweight_numer; //!< Simple coweight numerators.
  LT::WeightList d_radicalBasis; //!< Basis for orthogonal to coroots.
  LT::WeightList d_coradicalBasis; //!< Basis for orthogonal to roots.

/*!
\brief Sum of the positive roots.
*/
  LT::Weight d_2rho;
  LT::Weight d_dual_2rho;

  LT::LatticeCoeff Cartan_denom; //!< Denominator for (co)|weight_numer|

/*!\brief BitSet recording whether the root datum is adjoint/simply connected.

  "Adjoint" here means that the center of the complex group determined by the
  root datum is connected. "Simply connected" means that the derived group of
  that complex group is simply connected. These two properties are exchanged
  by passage to the dual root datum.
*/
  Status d_status;


 public:

// constructors and destructors

 RootDatum()
   : RootSystem(LT::LatticeMatrix(0,0))
   , d_rank(0)
  {}

  explicit RootDatum(const prerootdata::PreRootDatum&);

  RootDatum(const RootDatum&, tags::DualTag);

  RootDatum(LT::LatticeMatrix&, const RootDatum&, tags::DerivedTag);

  RootDatum(LT::LatticeMatrix&, const RootDatum&, tags::SimplyConnectedTag);

  ~RootDatum();

// accessors

  const RootSystem& root_system() const { return *this; } // base object ref
  size_t rank() const { return d_rank; }
  size_t semisimpleRank() const { return RootSystem::rank(); }

// root list access

  LT::WeightList::const_iterator beginRoot() const
    { return d_roots.begin(); }

  LT::WeightList::const_iterator endRoot() const
    { return d_roots.end(); }

  LT::WeightList::const_iterator beginCoroot() const
    { return d_coroots.begin(); }

  LT::WeightList::const_iterator endCoroot() const
    { return d_coroots.end(); }

  LT::WeightList::const_iterator beginRadical() const
    { return d_radicalBasis.begin(); }

  LT::WeightList::const_iterator endRadical() const
    { return d_radicalBasis.end(); }

  LT::WeightList::const_iterator beginCoradical() const
    { return d_coradicalBasis.begin(); }

  LT::WeightList::const_iterator endCoradical() const
    { return d_coradicalBasis.end(); }

  WRootIterator beginSimpleRoot() const // simple roots start halfway
    { return beginRoot()+numPosRoots(); }

  WRootIterator endSimpleRoot() const // and end after |semisimpleRank()|
    { return beginSimpleRoot()+semisimpleRank(); }

  WRootIterator beginPosRoot() const // positive roots start halfway
    { return  beginSimpleRoot(); }

  WRootIterator endPosRoot() const // an continue to the end
    { return endRoot(); }

  WRootIterator beginSimpleCoroot() const
    { return beginCoroot()+numPosRoots(); }

  WRootIterator endSimpleCoroot() const
    { return beginSimpleCoroot()+semisimpleRank(); }

  WRootIterator beginPosCoroot() const // positive coroots start halfway
    { return  beginSimpleCoroot(); }

  WRootIterator endPosCoroot() const
    { return endCoroot(); }


  bool isRoot(const LT::Weight& v) const // ask this of a weight
    { return setutils::find_index(d_roots,v) != d_roots.size(); }

  bool isSemisimple() const { return d_rank == semisimpleRank(); }

  const LT::Weight& root(RootNbr j) const
    { assert(j<numRoots()); return d_roots[j]; }

  const LT::Weight& simpleRoot(size_t j) const
    { assert(j<semisimpleRank()); return *(beginSimpleRoot()+j); }

  const LT::Weight& posRoot(size_t j) const
    { assert(j<numPosRoots()); return *(beginPosRoot()+j); }

  RootNbr rootNbr(const Root& r) const
    { return setutils::find_index(d_roots,r); }


  const LT::Weight& coroot(RootNbr j) const
    { assert(j<numRoots()); return d_coroots[j]; }

  const LT::Weight& simpleCoroot(size_t j) const
    { assert(j<semisimpleRank()); return *(beginSimpleCoroot()+j); }

  const LT::Weight& posCoroot(size_t j) const
    { assert(j<numPosRoots()); return  *(beginPosCoroot()+j); }

// other accessors

/*!
\brief Tells whether the rootdatum is the rootdatum of an adjoint group.

  NOTE: we define a reductive group to be adjoint if its center is
  connected.  An equivalent condition is that the derived group
  of the dual group is simply connected.
*/
  bool isAdjoint() const { return d_status[IsAdjoint]; }

/*!
\brief Tells whether the rootdatum is the rootdatum of a simply connected
  group.

  NOTE: this is the dual condition to being adjoint: it means
  that the derived group is simply connected.  An equivalent condition
  is that the center of the dual group is connected.
*/
  bool isSimplyConnected() const { return d_status[IsSimplyConnected]; }


  const LT::Weight& twoRho() const { return d_2rho; }
  const LT::Weight& dual_twoRho() const { return d_dual_2rho; }


  LT::LatticeCoeff scalarProduct(const LT::Weight& v, RootNbr j) const
  {
    return v.scalarProduct(coroot(j));
  }

  bool isOrthogonal(const latticetypes::Weight& v, RootNbr j) const {
    return v.scalarProduct(coroot(j))==0;
  }

  bool isOrthogonal(RootNbr alpha, RootNbr j) const
  {
    return isOrthogonal(root(alpha),j);
  }

  LT::LatticeCoeff cartan(size_t i, size_t j) const {
    return simpleRoot(i).scalarProduct(simpleCoroot(j));
  }

  //!\brief  Applies to v the reflection about root alpha.
  void reflect(LT::Weight& lambda, RootNbr alpha) const
    { lambda -= d_roots[alpha]*lambda.scalarProduct(d_coroots[alpha]); }
  //!\brief  Applies reflection about coroot |alpha| to a coweight
  void coreflect(LT::Weight& co_lambda, RootNbr alpha) const
    { co_lambda -= d_coroots[alpha]*co_lambda.scalarProduct(d_roots[alpha]); }

  LT::Weight reflection(LT::Weight lambda, RootNbr alpha) const
    { reflect(lambda,alpha); return lambda; }
  LT::Weight coreflection(LT::Weight co_lambda, RootNbr alpha) const
    { coreflect(co_lambda,alpha); return co_lambda; }

  void simpleReflect(LT::Weight& v, size_t i) const
    { reflect(v,simpleRootNbr(i)); }
  void simpleCoreflect(LT::Weight& v, size_t i) const
    { coreflect(v,simpleRootNbr(i)); }

  LT::Weight simpleReflection(LT::Weight lambda, size_t i) const
    { simpleReflect(lambda,i); return lambda; }
  LT::Weight simpleCoreflection(LT::Weight co_lambda, size_t i) const
    { simpleCoreflect(co_lambda,i); return co_lambda; }


  // here any matrix permuting the roots is allowed, e.g., root_reflection(r)
  setutils::Permutation rootPermutation(const LT::LatticeMatrix& q) const;
  // extend diagram automorphism to permutation of all roots

  LT::LatticeMatrix root_reflection(RootNbr r) const;
  LT::LatticeMatrix simple_reflection(size_t i) const
    { return root_reflection(simpleRootNbr(i)); }

  LT::LatticeMatrix matrix(const weyl::WeylWord& ww) const;

  weyl::WeylWord reflectionWord(RootNbr r) const;

  // express root in basis of simple roots
  LT::Weight inSimpleRoots(RootNbr alpha) const { return root_expr(alpha); }
  // express coroot in basis of simple coroots
  LT::Weight inSimpleCoroots(RootNbr alpha) const { return coroot_expr(alpha); }

  LT::Weight twoRho(const RootList&) const;
  LT::Weight twoRho(const RootSet&) const;


  weyl::WeylWord word_of_inverse_matrix(const latticetypes::LatticeMatrix&)
    const;

// manipulators

  void swap(RootDatum&);

// private methods used during construction
 private:

  void fillStatus();


}; // class RootDatum

   /*!
   \brief Iterator for traversing a set of roots.

   This template class has two main instances, with I=bitmap;:BitMap::iterator
   and with I=RootList::const_iterator. Apart from the iterator
   |d_po| of this type the class stores a iterator |d_list|, accessing a set
   of roots. Only |d_pos| changes as the iterator advances, while |d_list|
   always points to the start of some |LT::WeightList|, so it can be |const|
   */
template<typename I>
class RootIterator { // constant Random Access Iterator

 private:

   const LT::WeightList::const_iterator d_list;
   I d_pos;

 public:

// associated types
   typedef LT::Weight value_type;
   typedef ptrdiff_t difference_type;
   typedef const LT::Weight& reference;
   typedef const LT::Weight* pointer;
   typedef std::random_access_iterator_tag iterator_category;

// constructors and destructors
   RootIterator() {}

   RootIterator(const RootIterator& i)
     :d_list(i.d_list),d_pos(i.d_pos) {}

   RootIterator(const LT::WeightList& wl, I i)
     :d_list(wl.begin()), d_pos(i) {}

   RootIterator(LT::WeightList::const_iterator wl, I i)
     :d_list(wl), d_pos(i) {}

   RootIterator(const RootDatum& rd, I i)
     :d_list(rd.beginRoot()), d_pos(i) {}

   ~RootIterator() {}

// assignment
   RootIterator& operator= (const RootIterator& i) {
     d_list = i.d_list;
     d_pos = i.d_pos;
     return *this;
   }

// comparison
   bool operator== (const RootIterator& i) const {
     return d_pos == i.d_pos;
   }

   bool operator!= (const RootIterator& i) const {
     return d_pos != i.d_pos;
   }

// iterator operations
   const LT::Weight& operator* () {
     return d_list[*d_pos];
   }

   RootIterator& operator++ () {
     ++d_pos;
     return *this;
   }

   RootIterator operator++ (int) {
     return RootIterator(d_list,d_pos++);
   }

   RootIterator& operator-- () {
     --d_pos;
     return *this;
   }

   RootIterator operator-- (int) {
     return RootIterator(d_list,d_pos--);
   }

   RootIterator& operator+= (difference_type n) {
     d_pos += n;
     return *this;
   }

   RootIterator operator+ (difference_type n) {
     return RootIterator(d_list,d_pos+n);
   }

   RootIterator& operator-= (difference_type n) {
     d_pos -= n;
     return *this;
   }

   RootIterator operator- (difference_type n) {
     return RootIterator(d_list,d_pos-n);
   }

   difference_type operator- (const RootIterator& i) const {
     return d_pos-i.d_pos;
   }

   const LT::Weight& operator[] (difference_type n) {
     return d_list[d_pos[n]];
   }
}; // |template<typename I> class RootIterator|

// non-member operator + which interchanges operand order
template<typename I>
inline RootIterator<I> operator+ (typename RootIterator<I>::difference_type n,
				  RootIterator<I> i)
  {
    return i+n;
  }


} // |namespace rootdata|

} // |namespace atlas|

#include "rootdata_def.h"

#endif
