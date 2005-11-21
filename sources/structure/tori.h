/*
  This is tori.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef TORI_H  /* guard against multiple inclusions */
#define TORI_H

#include "tori_fwd.h"

#include "bits.h"
#include "latticetypes.h"
#include "subquotient.h"
#include "tags.h"

namespace atlas {

namespace tori {

  namespace LT = latticetypes;

}

/******** function declarations *********************************************/

namespace tori {

  void dualPi0(LT::ComponentSubquotient&, const LT::LatticeMatrix&);

  void minusBasis(LT::WeightList&, const LT::LatticeMatrix&);

  void minusMatrix(LT::LatticeMatrix&, const LT::LatticeMatrix&, 
		   const RealTorus&);

  void plusBasis(LT::WeightList&, const LT::LatticeMatrix&);

  void plusMatrix(LT::LatticeMatrix&, const LT::LatticeMatrix&, 
		  const RealTorus&);
}

/******** type definitions **************************************************/

namespace tori {

class RealTorus {

 private:

  size_t d_rank;                             // rank of torus
  size_t d_complexRank;                      // number of C^x-factors
  LT::LatticeMatrix d_involution;            // matrix of the involution
  LT::WeightList d_plus;                     // basis for +1 eigenlattice
  LT::WeightList d_minus;                    // basis for -1 eigenlattice
  LT::LatticeMatrix d_toPlus;                // projection onto X_+
  LT::LatticeMatrix d_toMinus;               // projection onto X_-

  LT::ComponentSubquotient d_topology;

 public:

// constructors and destructors
  RealTorus() 
    {}

  explicit RealTorus(const LT::LatticeMatrix&);

  RealTorus(const RealTorus&, tags::DualTag);

  ~RealTorus() 
    {}

// accessors
  size_t compactRank() const {
    return d_plus.size()-d_complexRank;
  }

  size_t complexRank() const {
    return d_complexRank;
  }

  void componentMap(LT::ComponentMap&, const LT::LatticeMatrix&, 
		    const RealTorus&) const;

  const LT::LatticeMatrix& involution() const {
    return d_involution;
  }

  bool isSplit() const {
    return d_minus.size() == d_rank;
  }

  const LT::WeightList& minusLattice() const {
    return d_minus;
  }

  size_t minusRank() const {
    return d_minus.size();
  }

  const LT::WeightList& plusLattice() const {
    return d_plus;
  }

  size_t plusRank() const {
    return d_plus.size();
  }

  size_t rank() const {
    return d_rank;
  }

  size_t splitRank() const {
    return d_minus.size()-d_complexRank;
  }

  void toMinus(LT::Weight& dest, const LT::Weight& source) const {
    d_toMinus.apply(dest,source);
  }

  void toPlus(LT::Weight& dest, const LT::Weight& source) const {
    d_toPlus.apply(dest,source);
  }

  const LT::ComponentSubquotient& topology() const {
    return d_topology;
  }

  size_t twoRank() const {
    return d_rank-d_complexRank;
  }
};

}

}

#endif
