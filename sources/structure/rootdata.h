/*
  This is rootdata.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

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

void lieType(lietype::LieType&, const RootList&, const RootDatum&);

void longest(LT::LatticeMatrix&, const RootDatum&);

void makeOrthogonal(RootList&, const RootList&, const RootList&,
		    const RootDatum&);

void reflectionMatrix(LT::LatticeMatrix&, RootNbr, const RootDatum&);

void reflectionMatrix(LT::LatticeMatrix&, RootNbr, const RootDatum&,
                      RootBasisTag);

void rootBasis(RootList&, const RootList&, const RootDatum&);

void rootBasis(RootList&, const RootSet&, const RootDatum&);

void rootPermutation(RootList&, const LT::LatticeMatrix&, const RootDatum&);

void strongOrthogonalize(RootList&, const RootDatum&);

bool sumIsRoot(const LT::Weight&, const LT::Weight&, const RootDatum&);

bool sumIsRoot(const RootNbr&, const RootNbr&, const RootDatum&);

void toDistinguished(LT::LatticeMatrix&, const RootDatum&);

void toMatrix(LT::LatticeMatrix&, const weyl::WeylWord&, const RootDatum&);

void toMatrix(LT::LatticeMatrix&, const RootList&, const RootDatum&);

void toPositive(weyl::WeylWord&, const LT::Weight&, const RootDatum&);

void toWeylWord(weyl::WeylWord&, RootNbr, const RootDatum&);

void toWeylWord(weyl::WeylWord&, const latticetypes::LatticeMatrix&, 
		const RootDatum&);

void twoRho(LT::Weight&, const RootList&, const RootDatum&);

void twoRho(LT::Weight&, const RootSet&, const RootDatum&);

}

/******** type definitions **************************************************/

namespace rootdata {

class RootDatum {

 private:

  enum StatusFlagNames { IsAdjoint, IsSimplyConnected, numFlags };
  typedef bitset::BitSet<numFlags> Status;

  size_t d_rank;                    // rank of the group
  size_t d_semisimpleRank;          // rank of the derived group   
  LT::WeightList d_coradicalBasis;  // basis for orthogonal to derived group
  LT::WeightList d_radicalBasis;    // same in the dual group
  LT::WeightList d_roots;           // full list of roots
  LT::WeightList d_coroots;         // full list of coroots
  RootList d_minus;                 // lists the negative of each root
  RootList d_posRoots;              // indices of positive roots
  RootList d_simpleRoots;           // indices of simple roots
  LT::RatWeightList d_weights;      // simple weights
  LT::RatWeightList d_coweights;    // simple coweights
  std::vector<setutils::Permutation> d_rootPermutation; 
                                    // records simple root permutations
  RootSet d_isPositive;             // for easy lookup of positivity
  RootSet d_isSimple;               // for easy lookup of simplicity
  LT::Weight d_twoRho;              // the sum of the positive roots
  Status d_status;

  void fillStatus();

 public:

// constructors and destructors

  RootDatum() 
    {}

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

  RootIterator beginPosCoroot() const;

  RootIterator endPosCoroot() const;

  RootIterator beginPosRoot() const;

  RootIterator endPosRoot() const;

  RootIterator beginSimpleCoroot() const;

  RootIterator endSimpleCoroot() const;

  RootIterator beginSimpleRoot() const;

  RootIterator endSimpleRoot() const;

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
    return !LT::scalarProduct(v,coroot(j));
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

  const setutils::Permutation& rootPermutation(size_t j) const {
    return d_rootPermutation[j];
  }

  void rootReflection(LT::LatticeMatrix& q, RootNbr j) const;

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

  RootNbr rootPermutation(RootNbr i, size_t s) const {
    return d_rootPermutation[s][i];
  }

  template <typename I, typename O> 
    void toRootBasis(const I&, const I&, O) const;

  const LT::Weight& twoRho() const {
    return d_twoRho;
  }

// manipulators

  void swap(RootDatum&);
};

// NOTE: this should really be a template, depending on a RandomAccessIterator
// with value_type RootNbr (so that in particular d_pos could be a pointer to
// a RootNbr)

class RootIterator { // constant Random Access Iterator

 private:

   LT::WeightList::const_iterator d_list;
   RootList::const_iterator d_pos;

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

}

}

#include "rootdata_def.h"

#endif
