/*!
\file
\brief Forward declarations of classes and types for atlas namespace.

*/
/*
  This is atlas_types.h

  Copyright (C) 2011 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef ATLAS_TYPES_H  /* guard against multiple inclusions */
#define ATLAS_TYPES_H

/*

 This module defines types in the (global) atlas namespace that can be
 therefore used in short form by all modules that include this file.

 It is just a forward-declaration file, and is intended to replace several
 such files. One still needs to include the related header files for these
 types to be complete.

 */

#include <vector>

#include "constants.h"

#include "bitmap_fwd.h"
#include "matrix_fwd.h"
#include "ratvec_fwd.h"

/******** forward type declarations ******************************************/

namespace atlas {

  // interpetationless terminology
  typedef matrix::Vector<int> int_Vector;
  typedef matrix::Matrix<int> int_Matrix;
  typedef std::vector<int_Vector> int_VectorList;

  // when related to a root system, these alternatives can be used
  typedef int_Vector Weight;
  typedef int_Vector Coweight;
  typedef ratvec::RationalVector<int> RatWeight;
  typedef ratvec::RationalVector<int> RatCoweight;
  typedef int_Matrix WeightInvolution;
  typedef int_Matrix CoweightInvolution;

  typedef std::vector<Weight> WeightList;
  typedef std::vector<Weight> CoweightList;
  typedef std::vector<RatWeight> RatWeightList;

  // more generic lattice-ralated terms
  typedef int LatticeCoeff; // the instance type of |Vector| and |Matrix|
  typedef std::vector<LatticeCoeff> CoeffList; // no vector arithmetic here
  typedef matrix::Vector<LatticeCoeff> LatticeElt;
  typedef std::vector<LatticeElt> latticeEltList;
  typedef matrix::Matrix<LatticeCoeff> LatticeMatrix;

  namespace bitvector {
    template<size_t> class BitVector;
    template<size_t> class BitVectorList;
    template<size_t> class BitMatrix;
  }
  typedef bitvector::BitVector<constants::RANK_MAX> SmallBitVector;
  typedef bitvector::BitVectorList<constants::RANK_MAX> SmallBitVectorList;
  typedef bitvector::BitVector<constants::RANK_MAX+1> BinaryEquation;
  typedef bitvector::BitVectorList<constants::RANK_MAX+1> BinaryEquationList;
  typedef bitvector::BitMatrix<constants::RANK_MAX> BinaryMap;

  namespace subquotient {
    template<size_t dim> class Subspace;
    template<size_t dim> class Subquotient;
  }
  typedef subquotient::Subspace<constants::RANK_MAX> SmallSubspace;
  typedef subquotient::Subquotient<constants::RANK_MAX> SmallSubquotient;

  namespace lietype {
    class SimpleLieType;
    class LieType;
    class InnerClassType;
    struct Layout;
  }
  using lietype::LieType;

  namespace prerootdata { class PreRootDatum; }
  using prerootdata::PreRootDatum;

  namespace rootdata {
    class RootSystem;
    class RootDatum;
    typedef int_Vector Root;
    typedef WeightList::const_iterator WRootIterator;
  }
  using rootdata::RootSystem;
  using rootdata::RootDatum;
  typedef unsigned short RootNbr;
  typedef std::vector<RootNbr> RootNbrList;
  typedef bitmap::BitMap RootNbrSet;

  namespace weyl {
    typedef unsigned char Generator; // index of simple root / simple reflection
    struct WeylWord : public std::vector<Generator>  { }; // in weyl namesace

    class WeylElt;
    class WeylGroup;
    class TwistedWeylGroup;

    typedef WeylElt TwistedInvolution;
    typedef std::vector<WeylElt> WeylEltList;
    typedef std::vector<TwistedInvolution> TwistedInvolutionList;

    typedef std::vector<signed char> InvolutionWord;
  }
  using weyl::WeylWord;
  using weyl::WeylElt;
  using weyl::WeylGroup;
  using weyl::TwistedWeylGroup;
  using weyl::TwistedInvolution;
  using weyl::WeylEltList;
  using weyl::TwistedInvolutionList;

} // |namespace atlas|

#endif
