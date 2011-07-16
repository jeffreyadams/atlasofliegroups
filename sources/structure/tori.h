/*!
\file
\brief Class definition and function declarations for the class RealTorus.
*/
/*
  This is tori.h
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef TORI_H  /* guard against multiple inclusions */
#define TORI_H

#include "bits.h"
#include "tags.h"

#include "atlas_types.h"
#include "subquotient.h"

namespace atlas {

/******** function declarations *********************************************/

namespace tori {

  void dualPi0(SmallSubquotient&, const WeightInvolution&);

  void plusMatrix(WeightInvolution&, const WeightInvolution&,
		  const RealTorus&);
  void minusMatrix(WeightInvolution&, const WeightInvolution&,
		   const RealTorus&);

  void plusBasis(WeightList&, const WeightInvolution&);

  WeightList plusBasis(const WeightInvolution&);

  void minusBasis(WeightList&, const WeightInvolution&);

  WeightList minusBasis(const WeightInvolution&);

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

  The most delicate invariant we shall have to deal with is the component
  group of the group of real points of T. This is an elementary abelian
  2-group; we shall rather consider its dual dpi0(T). To describe its rank is
  fairly easy. Indeed, T may be decomposed as a product of compact, split and
  complex factors; denote r_u, r_s and r_c the number of factors of each type,
  so that the rank n of T is r_u + r_s + 2 r_c, then the rank of the component
  group is r_s. We have moreover : rk(X_+) = r_u + r_c; rk(X_-) = r_s + r_c.
  Denote V = X/2X, a vector space over the two-element field F_2, and denote
  V_+,V_- the images of X_+, X_- in V. Then it is not hard to show that r_c =
  dim(V_+ cap V_-), which allows computing r_s once rk(X_-) is known.
  One may prove that V_+- := V_+ cap V_- is also the image of tau - 1 in V.

  It is a little bit harder to describe dpi0(T) functorially as a vector
  space. The group T(2)(R) of real points of the group of elements of order 2
  in T is in natural duality with V/V_+-. There is a natural surjection from
  T(2)(R) to pi0(T), so dpi0(T) is a sub-vector space of V/V_+-, which may in
  fact be described as the image of V_-. This is how we will consider it.
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
  WeightInvolution d_involution;            // matrix of the involution

  /*!
  basis for +1 eigenlattice of the Cartan involution
  */
  WeightList d_plus;                     // basis for +1 eigenlattice

  /*!
  basis for -1 eigenlattice of the Cartan involution
  */
  WeightList d_minus;                    // basis for -1 eigenlattice

  /*!
  coordinate transformation from standard basis of $X$ to basis |d_plus| of
  $X_+$; should be applied only to elements of $X_+$
  */
  LatticeMatrix d_toPlus;                // transform coordinates to $X_+$

  /*!
  coordinate transformation from standard basis of $X$ to basis |d_minus| of
  $X_-$; should be applied only to elements of $X_-$
  */
  LatticeMatrix d_toMinus;               // transform coordinates to $X_-$

  /*!
  dual component group of real torus (a vector space over $Z/2Z$), realised
  as the subquotient $(V_+ + V_-)/V_+$ of the $Z/2Z$ vector space $X/2X$
  */
  SmallSubquotient d_topology;

 public:

// constructors and destructors
  RealTorus()
    {}

  explicit RealTorus(const WeightInvolution&);

  RealTorus(const RealTorus&, tags::DualTag);

  ~RealTorus()
    {}

// accessors

  const WeightInvolution& involution() const { return d_involution; }

  size_t rank() const { return d_rank; }
  size_t complexRank() const { return d_complexRank; }
  size_t compactRank() const { return d_plus.size()-d_complexRank; }
  size_t splitRank() const { return d_minus.size()-d_complexRank; }
  size_t plusRank() const { return d_plus.size(); }
  size_t minusRank() const { return d_minus.size(); }
  size_t twoRank() const { return d_rank-d_complexRank; }

  bool isCompact() const { return d_plus.size() == d_rank; }
  bool isSplit() const { return d_minus.size() == d_rank; }

  BinaryMap componentMap(const LatticeMatrix&,
				    const RealTorus&) const;

  const WeightList& plusLattice() const { return d_plus; }
  const WeightList& minusLattice() const { return d_minus; }

  void toPlus(Weight& dest, const Weight& source) const
    { dest=d_toPlus*source; }

  void toMinus(Weight& dest, const Weight& source) const
    { dest=d_toMinus*source; }

  const SmallSubquotient& topology() const { return d_topology; }
};

}

}

#endif
