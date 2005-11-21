/*
  This is latticetypes_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
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

  typedef int LatticeCoeff;

  const LatticeCoeff ZeroCoeff = 0;
  const LatticeCoeff OneCoeff = 1;

  typedef std::vector<LatticeCoeff> LatticeElt;
  typedef std::vector<LatticeCoeff> CoeffList;
  class RatLatticeElt;

  typedef LatticeElt Weight;
  typedef RatLatticeElt RatWeight;

  typedef std::vector<Weight> WeightList;
  typedef std::vector<RatWeight> RatWeightList;

  typedef matrix::Matrix<LatticeCoeff> LatticeMatrix;

  typedef bitvector::BitVector<constants::RANK_MAX> Component;
  typedef bitvector::BitVector<2*constants::RANK_MAX> LongComponent;
  typedef bitvector::BitMatrix<constants::RANK_MAX> ComponentMap;

  typedef subquotient::NormalSubspace<constants::RANK_MAX> ComponentSubspace;
  typedef subquotient::NormalSubquotient<constants::RANK_MAX> 
    ComponentSubquotient;

  typedef std::vector<Component> ComponentList;
  typedef std::vector<LongComponent> LongComponentList;

}

}

#endif
