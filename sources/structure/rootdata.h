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

#include "bitset.h"
#include "bitmap.h"
#include "latticetypes.h"
#include "lietype.h"
#include "setutils.h"
#include "tags.h"

namespace atlas {

namespace rootdata {

  namespace LT = latticetypes;

}

/******** type declarations *************************************************/

namespace rootdata {

  struct RootBasisTag {};

}

/******** function declarations *********************************************/

namespace rootdata {


void cartanMatrix(LT::LatticeMatrix&, const RootDatum&);

void cartanMatrix(LT::LatticeMatrix&, const RootList&, const RootDatum&);

void dualBasedInvolution(LT::LatticeMatrix&, const LT::LatticeMatrix&,
			 const RootDatum&);

// a functional version of previous one
LT::LatticeMatrix dualBasedInvolution
  (const LT::LatticeMatrix&, const RootDatum&);

void lieType(lietype::LieType&, const RootList&, const RootDatum&);

void longest(LT::LatticeMatrix&, const RootDatum&);

void makeOrthogonal(RootList&, const RootList&, const RootList&,
		    const RootDatum&);

void reflectionMatrix(LT::LatticeMatrix&, RootNbr, const RootDatum&);

/* not implemented (yet?), would give matrix on basis of simple roots
void reflectionMatrix(LT::LatticeMatrix&, RootNbr, const RootDatum&,
                      RootBasisTag);
*/

// compute a basis of a root subsystem
void rootBasis(RootList&, const RootList&, const RootDatum&);

void rootBasis(RootList&, RootSet, const RootDatum&);

void strongOrthogonalize(RootList&, const RootDatum&);

void toDistinguished(LT::LatticeMatrix&, const RootDatum&);

void toMatrix(LT::LatticeMatrix&, const weyl::WeylWord&, const RootDatum&);

void toMatrix(LT::LatticeMatrix&, const RootList&, const RootDatum&);

void toPositive(weyl::WeylWord&, const LT::Weight&, const RootDatum&);

void toWeylWord(weyl::WeylWord&, RootNbr, const RootDatum&);

} // namespace rootdata

/******** type definitions **************************************************/

namespace rootdata {
  /*!
  \brief Based root datum for a complex reductive group.

  What we call a root datum in this program is what is usually called
  a based root datum.

  The root datum defines the complex reductive group entirely.

  The lattices in which the roots and coroots live are both Z^d_rank;
  lists of roots or coroots are lists of vectors of integers, of size
  d_rank. RootDatum begins with the lists of simple roots and coroots,
  then constructs the remaining roots.  Also constructed are various
  useful auxiliary things, like d_twoRho (the sum of the positive
  roots) and d_rootPermutation (a list of permutations of the roots,
  with the jth permutation given by reflection in the jth simple
  root).

  The code is designed to make it preferable always to refer to a root
  by its number (its location in the list d_root).  The type RootList
  should in fact (per Fokko) have been called RootNbrList, since it is
  a list of such numbers.
  */
class RootDatum {

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

/*!
\brief Semisimple rank of the root datum.
*/
  size_t d_semisimpleRank;

  LT::WeightList d_roots; //!< Full list of roots.
  LT::WeightList d_coroots; //!< Full list of coroots.
  LT::WeightList weight_numer; //!< Simple weight numerators.
  LT::WeightList coweight_numer; //!< Simple coweight numerators.
  LT::WeightList d_radicalBasis; //!< Basis for orthogonal to coroots.
  LT::WeightList d_coradicalBasis; //!< Basis for orthogonal to roots.

/*!
\brief Root permutations induced by simple reflections.
*/
  std::vector<setutils::Permutation> d_rootPermutation;

/*!
\brief Sum of the positive roots.
*/
  LT::Weight d_2rho;
  LT::Weight d_dual_2rho;

  LT::LatticeCoeff Cartan_denom;//!< Denominator |weight_numer|,|coweight_numer|
  /*!
\brief BitSet recording whether the root datum is adjoint/simply connected.

  "Adjoint" here means that the center of the complex group determined by the
  root datum is connected. "Simply connected" means that the derived group of
  that complex group is simply connected. These two properties are exchanged
  by passage to the dual root datum.
  */
  Status d_status;


 public:

// constructors and destructors

  RootDatum() : d_rank(0), d_semisimpleRank(0) {}

  explicit RootDatum(const prerootdata::PreRootDatum&);

  RootDatum(const RootDatum&, tags::DualTag);

  RootDatum(LT::LatticeMatrix&, const RootDatum&, tags::DerivedTag);

  RootDatum(LT::LatticeMatrix&, const RootDatum&, tags::SimplyConnectedTag);

  ~RootDatum();

// accessors

  unsigned long numRoots() const { return d_roots.size(); }

  unsigned long numPosRoots() const { return numRoots()/2; }

  size_t rank() const { return d_rank; }

  size_t semisimpleRank() const { return d_semisimpleRank; }

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


  bool isSimpleRoot(RootNbr j) const
    { return j-numPosRoots()<semisimpleRank(); } // uses |RootNbr| is unsigned

  bool isPosRoot(RootNbr j) const { return j>=numPosRoots(); } // second half

  bool isRoot(const LT::Weight& v) const // ask this of a weight
    { return setutils::find_index(d_roots,v) != d_roots.size(); }

  bool isSemisimple() const { return d_rank == d_semisimpleRank; }

  const LT::Weight& root(RootNbr j) const
    { assert(j<numRoots()); return d_roots[j]; }

  RootNbr simpleRootNbr(size_t j) const
    { assert(j<semisimpleRank());  return beginSimpleRoot()-beginRoot()+j; }

  const LT::Weight& simpleRoot(size_t j) const
    { assert(j<semisimpleRank()); return *(beginSimpleRoot()+j); }

  RootNbr posRootNbr(size_t j) const
    { assert(j<numPosRoots()); return beginPosRoot()-beginRoot()+j; }

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


  RootNbr rootMinus(RootNbr j) const // roots are ordered symmetrically
    { return numRoots()-1-j; }


  RootSet simpleRootSet() const // NOT for iteration over it, seems never used
    {
      RootSet simple_roots(numRoots());
      for (size_t i=0; i<semisimpleRank(); ++i)
	simple_roots.insert(simpleRootNbr(i));
      return simple_roots;
    }

  RootList simpleRootList() const // NOT for iteration over it
    {
      RootList simple_roots(semisimpleRank());
      for (size_t i=0; i<semisimpleRank(); ++i)
	simple_roots[i]=simpleRootNbr(i);
      return simple_roots;
    }

  RootSet posRootSet() const // NOT for iteration over it, may serve as mask
    {
      RootSet pos_roots(numRoots());
      for (size_t i=0; i<numPosRoots(); ++i)
	pos_roots.insert(posRootNbr(i));
      return pos_roots;
    }

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



  bool isOrthogonal(const latticetypes::Weight& v, RootNbr j) const {
    return v.scalarProduct(coroot(j))==0;
  }

  bool isOrthogonal(RootNbr i, RootNbr j) const {
    return isOrthogonal(root(i),j);
  }

  const LT::Weight& twoRho() const { return d_2rho; }
  const LT::Weight& dual_twoRho() const { return d_dual_2rho; }


  LT::LatticeCoeff scalarProduct(const LT::Weight& v, RootNbr j) const
  {
    return v.scalarProduct(coroot(j));
  }

  // This method is badly named: the relation is asymmetric
  LT::LatticeCoeff scalarProduct(RootNbr i, RootNbr j) const
  {
    return root(i).scalarProduct(coroot(j));
  }

  LT::LatticeCoeff cartan(size_t i, size_t j) const {
    return simpleRoot(i).scalarProduct(simpleCoroot(j));
  }

  LT::LatticeMatrix cartanMatrix() const;

  LT::LatticeMatrix cartanMatrix(const RootList&) const; // for subsystem

  void reflect(LT::Weight& lambda, RootNbr alpha) const;
  void coreflect(LT::Weight& co_lambda, RootNbr alpha) const;

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


  // the next method only works for _simple_ roots! (whence no RootNbr for |i|)
  const setutils::Permutation& rootPermutation(size_t i) const {
    return d_rootPermutation[i];
  }

  // here any matrix permuting the roots is allowed, e.g., rootReflection(r)
  setutils::Permutation rootPermutation(const LT::LatticeMatrix& q) const;

  void rootReflection(LT::LatticeMatrix& q, RootNbr r) const;

  LT::LatticeMatrix rootReflection(RootNbr r) const;

  void rootReflect(RootNbr& r, size_t s) const  { r=rootPermutation(s)[r]; }

  RootNbr permuted_root(const weyl::WeylWord& ww, RootNbr r) const
    {
      for (size_t i=ww.size(); i-->0; )
	rootReflect(r,ww[i]);
      return r;
    }

  RootNbr reflectedRoot(RootNbr r, size_t s) const {
    return d_rootPermutation[s][r];
  }

  weyl::WeylWord reflectionWord(RootNbr r) const;

  // express root in basis of simple roots
  LT::Weight inSimpleRoots(RootNbr n) const;

  // convert sequence of root numbers to expressions in the simple roots
template <typename I, typename O>
  void toRootBasis(const I&, const I&, O) const;
template <typename I, typename O>
  void toRootBasis(const I&, const I&, O, const RootList&) const;
template <typename I, typename O>
  void toSimpleWeights(const I&, const I&, O, const RootList&) const;

  LT::Weight twoRho(const RootList&) const;
  LT::Weight twoRho(const RootSet&) const;


  weyl::WeylWord word_of_inverse_matrix(const latticetypes::LatticeMatrix&)
    const;

  RootList simpleBasis(const RootList& rl) const;
  RootList simpleBasis(RootSet rs) const;


  bool sumIsRoot(const LT::Weight&, const LT::Weight&) const;
  bool sumIsRoot(RootNbr i, RootNbr j) const
    { return sumIsRoot(root(i),root(j)); }

  RootList high_roots() const;

// manipulators

  void swap(RootDatum&);

// private methods used during construction
 private:

  void fill_positives(latticetypes::WeightList& roots,
		     latticetypes::WeightList& coroots);

  std::vector<size_t>
  sort_roots(const latticetypes::WeightList& roots,
	     const latticetypes::WeightList& simple_coweights) const;

  void fillStatus();


}; // class RootDatum

   /*!
   \brief Iterator for traversing a set of roots.

   This template class has two main instances, with I=bitmap;:BitMap::iterator
   and with I=std::vector<RootNbr>::const_iterator. Apart from the iterator
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
