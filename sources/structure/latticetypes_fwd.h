/*!
\file
\brief Forward declarations of classes and types for namespace latticetypes.

*/
/*
  This is latticetypes_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef LATTICETYPES_FWD_H  /* guard against multiple inclusions */
#define LATTICETYPES_FWD_H

#include <vector>

#include "bitvector_fwd.h"
#include "matrix_fwd.h"
#include "subquotient_fwd.h"

#include "constants.h"

/******** forward type declarations ******************************************/

namespace atlas {

namespace latticetypes {

  typedef int LatticeCoeff; // the instance type of |Vector| and |Matrix|

  typedef matrix::Vector<LatticeCoeff> LatticeElt;

  typedef std::vector<LatticeCoeff> CoeffList; // no vector arithmetic here

  typedef matrix::Matrix<LatticeCoeff> LatticeMatrix;

  class RatLatticeElt;

  typedef LatticeElt Weight;
  typedef RatLatticeElt RatWeight;

  typedef std::vector<Weight> WeightList;
  typedef std::vector<RatWeight> RatWeightList;

  typedef bitvector::BitVector<constants::RANK_MAX> SmallBitVector;

  typedef bitvector::BitVectorList<constants::RANK_MAX> SmallBitVectorList;

  typedef bitvector::BitVector<constants::RANK_MAX+1> BinaryEquation;

  typedef bitvector::BitVectorList<constants::RANK_MAX+1> BinaryEquationList;

  //! \brief Square matrix of size at most RANK_MAX with entries in Z/2Z.
  typedef bitvector::BitMatrix<constants::RANK_MAX> BinaryMap;

  //! \brief Subgroup of (Z/2Z)^RANK_MAX.
  typedef subquotient::Subspace<constants::RANK_MAX> SmallSubspace;

  //! \brief Subquotient of (Z/2Z)^RANK_MAX.
  typedef subquotient::Subquotient<constants::RANK_MAX> SmallSubquotient;

} // |namespace latticetypes|

} // |namespace atlas|

#endif
