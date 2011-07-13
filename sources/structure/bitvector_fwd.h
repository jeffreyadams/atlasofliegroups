/*!
\file
\brief Class declarations for BitVector.
*/
/*
  This is bitvector_fwd.h
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef BITVECTOR_FWD_H  /* guard against multiple inclusions */
#define BITVECTOR_FWD_H

#include <cstddef>
#include <vector>

#include "constants.h"

/******** forward type declarations ******************************************/

namespace atlas {

namespace bitvector {

  template<size_t> class BitVector;
  template<size_t> class BitVectorList;
  template<size_t> class BitMatrix;


  typedef BitVector<constants::RANK_MAX> SmallBitVector;
  typedef BitVectorList<constants::RANK_MAX> SmallBitVectorList;

  typedef BitVector<constants::RANK_MAX+1> BinaryEquation;
  typedef BitVectorList<constants::RANK_MAX+1> BinaryEquationList;

  typedef BitMatrix<constants::RANK_MAX> BinaryMap;

}

}

#endif
