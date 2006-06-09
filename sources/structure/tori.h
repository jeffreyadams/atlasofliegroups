/*!
\file
\brief Class definition and function declarations for the class RealTorus.
*/
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

  /*!
  \brief Represents a torus defined over R. 

  This is equivalent to the datum of a lattice with an involution;
  to be consistent with the rest of the program, we use the Cartan
  involution tau which is the negative of the Galois involution on
  the character lattice.

  The fundamental data are the rank d_rank (allowing us to identify the
  character lattice X with d_rank-tuples of integers) and the integer
  matrix d_involution of tau.

  The software constructs two sublattices X_+ and X_-, the +1 and -1
  eigenspaces for the involution. The datum d_plus is a basis of X_+,
  and d_minus is a basis of X_-.  Both X_+ and X_- are supplementable
  in X, but in general X is not equal to the direct sum X_+ + X_-.

... the pi0 part is outdated ...

  The most delicate invariant we shall have to deal with is the
  component group of the group of real points of T. This is an
  elementary abelian 2-group; we shall rather consider its dual
  dpi0(T). To describe its rank is fairly easy : indeed, T may be
  decomposed as a product of compact, split and complex factors;
  denote r_u, r_s and r_c the number of factors of each type, so that
  the rank n of T is r_u + r_s + 2 r_c. Then we have : rk(X_+) = r_u +
  r_c; rk(X_-) = r_s + r_c. Denote V = X/2X, a vector space over the
  two-element field F_2, and denote V_+,V_- the images of X_+, X_- in
  V.  Then it is not hard to show that r_c = dim(V_+\cap V_-); one may
  prove that this is also the image of tau - 1 in V.

  It is a little bit harder to describe dpi0(T) functorially as a
  vector space.  Denote V_+- the intersection V_+ cap V_-. Then one
  may show that the group T(2)(R) of real points of the group of
  elements of order 2 in T is in natural duality with V/V_+-. There is
  a natural surjection from T(2)(R) to pi0(T), so dpi0(T) is a
  sub-vector space of V/V_+-, which may in fact be described as the
  image of V_-. This is how we will consider it.
  */
  class RealTorus {

 private:

  /*!
  rank of torus
  */
  size_t d_rank;                             // rank of torus
  
  /*!
  number of C^x-factors
  */
  size_t d_complexRank;                      // number of C^x-factors
  
  /*!
  matrix of the Cartan involution 
  */
  LT::LatticeMatrix d_involution;            // matrix of the involution
  
  /*!
  basis for +1 eigenlattice of the Cartan involution
  */
  LT::WeightList d_plus;                     // basis for +1 eigenlattice
  
  /*!
  basis for -1 eigenlattice of the Cartan involution
  */
  LT::WeightList d_minus;                    // basis for -1 eigenlattice
  
  /*!
  (some) projection onto +1 eigenlattice X_+; kernel _need not be_ the
  -1 eigenlattice 
  */
  LT::LatticeMatrix d_toPlus;                // projection onto X_+
  
  /*!
  (some) projection onto -1 eigenlattice X_-; kernel _need not be_ the
  +1 eigenlattice
  */
  LT::LatticeMatrix d_toMinus;               // projection onto X_-
  
  /*!
  group of connected components (a vector space over Z/2Z).
  */
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
