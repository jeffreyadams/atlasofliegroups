/*!
\file
\brief Class definitions and function declarations for the RootDatum class.
*/
/*
  This is rootdata.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef ROOTDATA_H  /* guard against multiple inclusions */
#define ROOTDATA_H

#include "rootdata_fwd.h"

#include <algorithm>

#include "prerootdata_fwd.h"
#include "weyl_fwd.h"

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

template <typename I, typename O>
  void toRootBasis(const I&, const I&, O, const RootList&, const RootDatum&);

template <typename I, typename O>
  void toSimpleWeights(const I&, const I&, O, const RootList&,
		       const RootDatum&);

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

void reflectionMatrix(LT::LatticeMatrix&, RootNbr, const RootDatum&,
                      RootBasisTag);

void rootBasis(RootList&, const RootList&, const RootDatum&);

void rootBasis(RootList&, RootSet, const RootDatum&);

void strongOrthogonalize(RootList&, const RootDatum&);

bool sumIsRoot(const LT::Weight&, const LT::Weight&, const RootDatum&);

bool sumIsRoot(const RootNbr&, const RootNbr&, const RootDatum&);

void toDistinguished(LT::LatticeMatrix&, const RootDatum&);

void toMatrix(LT::LatticeMatrix&, const weyl::WeylWord&, const RootDatum&);

void toMatrix(LT::LatticeMatrix&, const RootList&, const RootDatum&);

void toPositive(weyl::WeylWord&, const LT::Weight&, const RootDatum&);

void toWeylWord(weyl::WeylWord&, RootNbr, const RootDatum&);

}

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
/*!
\brief  basis for orthogonal to roots.
*/
  LT::WeightList d_coradicalBasis;
/*!
\brief Basis for orthogonal to coroots.
*/
  LT::WeightList d_radicalBasis;
/*!
\brief Full list of roots.
*/
  LT::WeightList d_roots;
/*!
\brief Full list of coroots.
*/
  LT::WeightList d_coroots;
/*!
\brief Lists the negative of each root.
*/
  RootList d_minus;
/*!
\brief Numbers of the positive roots.
*/
  RootList d_posRoots;
/*!
\brief Numbers of the simple roots.
*/
  RootList d_simpleRoots;
/*!
\brief  Simple weights.
*/
  LT::RatWeightList d_weights;
/*!
\brief Simple coweights.
*/
  LT::RatWeightList d_coweights;
/*!
\brief Root permutations induced by simple reflections.
*/
  std::vector<setutils::Permutation> d_rootPermutation;
/*!
\brief BitMap flagging positive roots.
*/
  RootSet d_isPositive;
/*!
\brief BitMap flagging simple roots.
*/
  RootSet d_isSimple;
/*!
\brief Sum of the positive roots.
*/
  LT::Weight d_twoRho;

  /*!
\brief BitSet recording in bit 0 whether the root datum is adjoint, and in
  bit 1 whether the root datum is simply connected.

  "Adjoint" here means that the center of the complex group determined by the
  root datum is connected. "Simply connected" means that the derived group of
  that complex group is simply connected. These two properties are exchanged
  by passage to the dual root datum.
  */
  Status d_status;

  void fillStatus();

 public:

// constructors and destructors

  RootDatum() : d_rank(0), d_semisimpleRank(0) {}

  explicit RootDatum(const prerootdata::PreRootDatum&);

  RootDatum(const RootDatum&, tags::DualTag);

  RootDatum(LT::LatticeMatrix&, const RootDatum&, tags::DerivedTag);

  RootDatum(LT::LatticeMatrix&, const RootDatum&, tags::SimplyConnectedTag);

  ~RootDatum();

// accessors
// root list access

  LT::WeightList::const_iterator beginCoradical() const {
    return d_coradicalBasis.begin();
  }

  LT::WeightList::const_iterator endCoradical() const {
    return d_coradicalBasis.end();
  }

  LT::WeightList::const_iterator beginRadical() const {
    return d_radicalBasis.begin();
  }

  LT::WeightList::const_iterator endRadical() const {
    return d_radicalBasis.end();
  }

  LT::WeightList::const_iterator beginCoroot() const {
    return d_coroots.begin();
  }

  LT::WeightList::const_iterator beginRoot() const {
    return d_roots.begin();
  }

  LT::WeightList::const_iterator endRoot() const {
    return d_roots.end();
  }

  LT::WeightList::const_iterator endCoroot() const {
    return d_coroots.end();
  }

  WRootIterator beginPosCoroot() const;

  WRootIterator endPosCoroot() const;

  WRootIterator beginPosRoot() const;

  WRootIterator endPosRoot() const;

  WRootIterator beginSimpleCoroot() const;

  WRootIterator endSimpleCoroot() const;

  WRootIterator beginSimpleRoot() const;

  WRootIterator endSimpleRoot() const;

  bool isPosRoot(RootNbr j) const {
    return d_isPositive.isMember(j);
  }

  bool isRoot(const LT::Weight& v) const {
    LT::WeightList::const_iterator first = d_roots.begin();
    LT::WeightList::const_iterator last = d_roots.end();
    return std::find(first,last,v) != last;
  }

  bool isSemisimple() const {
    return d_rank == d_semisimpleRank;
  }

  bool isSimpleRoot(RootNbr j) const {
    return d_isSimple.isMember(j);
  }

  const LT::Weight& coroot(RootNbr j) const {
    return d_coroots[j];
  }

  const LT::Weight& posCoroot(size_t j) const {
    return d_coroots[d_posRoots[j]];
  }

  const LT::Weight& posRoot(size_t j) const {
    return d_roots[d_posRoots[j]];
  }

  const LT::Weight& root(RootNbr j) const {
    return d_roots[j];
  }

  const LT::Weight& simpleCoroot(size_t j) const {
    return d_coroots[d_simpleRoots[j]];
  }

  const LT::Weight& simpleRoot(size_t j) const {
    return d_roots[d_simpleRoots[j]];
  }

  RootNbr posRootNbr(size_t j) const {
    return d_posRoots[j];
  }

  RootNbr rootNbr(const Root& r) const {
    LT::WeightList::const_iterator first = d_roots.begin();
    LT::WeightList::const_iterator last = d_roots.end();
    return find(first,last,r) - first;
  }

  RootNbr simpleRootNbr(size_t j) const {
    return d_simpleRoots[j];
  }

  const RootSet& posRootSet() const {
    return d_isPositive;
  }

  const RootList& simpleRootList() const {
    return d_simpleRoots;
  }

  const RootSet& simpleRootSet() const {
    return d_isSimple;
  }

// other accessors

  LT::LatticeCoeff cartan(size_t i, size_t j) const {
    const LT::Weight& arg1 = simpleRoot(i);
    const LT::Weight& arg2 = simpleCoroot(j);
    return LT::scalarProduct(arg1,arg2);
  }

  void coreflection(LT::Weight&, RootNbr) const;

  bool isAdjoint() const;

  bool isOrthogonal(const latticetypes::Weight& v, RootNbr j) const {
    return LT::scalarProduct(v,coroot(j))==0;
  }

  bool isOrthogonal(RootNbr i, RootNbr j) const {
    return isOrthogonal(root(i),j);
  }

  bool isSimplyConnected() const;

  unsigned long numPosRoots() const {
    return d_posRoots.size();
  }

  unsigned long numRoots() const {
    return d_roots.size();
  }

  size_t rank() const {
    return d_rank;
  }

  void reflection(LT::Weight&, RootNbr) const;

  RootNbr rootMinus(RootNbr j) const {
    return d_minus[j];
  }

  // the next method only works for _simple_ roots! (whence no RootNbr for |j|)
  const setutils::Permutation& rootPermutation(size_t j) const {
    return d_rootPermutation[j];
  }

  // here any matrix permuting the roots is allowed, e.g., rootReflection(r)
  setutils::Permutation rootPermutation(const LT::LatticeMatrix& q) const;

  void rootReflection(LT::LatticeMatrix& q, RootNbr r) const;

  LT::LatticeMatrix rootReflection(RootNbr r) const;

  void rootReflect(RootNbr& r, size_t i) const  { r=rootPermutation(i)[r]; }

  weyl::WeylWord reflectionWord(RootNbr r) const;

  LT::LatticeCoeff scalarProduct(const latticetypes::Weight& v, RootNbr j)
    const {
    return LT::scalarProduct(v,coroot(j));
  }

  LT::LatticeCoeff scalarProduct(RootNbr i, RootNbr j) const {
    const LT::Weight& arg1 = root(i);
    const LT::Weight& arg2 = coroot(j);
    return LT::scalarProduct(arg1,arg2);
  }

  size_t semisimpleRank() const {
    return d_semisimpleRank;
  }

  void simpleCoreflection(LT::Weight& v, size_t j) const {
    coreflection(v,d_simpleRoots[j]);
  }

  void simpleReflection(LT::Weight& v, size_t j) const {
    reflection(v,d_simpleRoots[j]);
  }

  LT::LatticeMatrix cartanMatrix() const;

  LT::LatticeMatrix cartanMatrix(const RootList&) const; // for subsystem


  RootNbr rootPermutation(RootNbr i, size_t s) const {
    return d_rootPermutation[s][i];
  }

  template <typename I, typename O>
    void toRootBasis(const I&, const I&, O) const;

  const LT::Weight& twoRho() const {
    return d_twoRho;
  }

  LT::Weight twoRho(const RootList&) const;

  LT::Weight twoRho(const RootSet&) const;


  weyl::WeylWord word_of_inverse_matrix(const latticetypes::LatticeMatrix&)
    const;

  RootList simpleBasis(const RootList& rl) const;
  RootList simpleBasis(RootSet rs) const;

// manipulators

  void swap(RootDatum&);

}; // class RootDatum

// NOTE: this should really be a template, depending on a RandomAccessIterator
// with value_type RootNbr (so that in particular d_pos could be a pointer to
// a RootNbr)

template<typename I>
   /*!
   \brief Iterator for traversing a set of roots.

   This template class is instantiated via WRootIterator, with I set to
   std::vector<RootNbr>::const_iterator.  Such a class stores two
   iterators: d_list, with value type LT::Weight, accessing a set of
   roots, and d_pos, with value type RootNbr, indexing that set of
   roots. Only the second one changes as the iterator advances, the first one
   always points to the start of some |LT::WeightList|, so it can be |const|
   */
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

   difference_type operator- (const RootIterator& j) const {
     return d_pos-j.d_pos;
   }

   const LT::Weight& operator[] (difference_type n) {
     return d_list[d_pos[n]];
   }
};

template<typename I>
inline RootIterator<I> operator+ (typename RootIterator<I>::difference_type n,
				  RootIterator<I> i)
  {
    return i+n;
  }

#if 0
  /*!
  Old non-template version of RootIterator, now excluded by the preprocessor.
  */
class RootIterator { // constant Random Access Iterator

 private:

   LT::WeightList::const_iterator d_list;
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

   RootIterator(const LT::WeightList& wl, RootList::const_iterator i)
     :d_list(wl.begin()), d_pos(i) {}

   RootIterator(LT::WeightList::const_iterator wl, RootList::const_iterator i)
     :d_list(wl), d_pos(i) {}

   RootIterator(const RootDatum& rd, RootList::const_iterator i)
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

   difference_type operator- (const RootIterator& j) const {
     return d_pos-j.d_pos;
   }

   const LT::Weight& operator[] (difference_type n) {
     return d_list[d_pos[n]];
   }
};

inline RootIterator operator+ (RootIterator::difference_type n, RootIterator i)
  {
    return i+n;
  }
#endif

}

}

#include "rootdata_def.h"

#endif
