/*
  This is abelian_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef ABELIAN_FWD_H  /* guard against multiple inclusions */
#define ABELIAN_FWD_H

#include <limits>
#include <set>

#include "latticetypes_fwd.h"
#include "matrix_fwd.h"

/******** type declarations **************************************************/

namespace atlas {

namespace abelian {

typedef unsigned long GrpNbr;
typedef std::vector<GrpNbr> GrpNbrList;
typedef std::vector<unsigned long> GrpArr;
typedef std::vector<GrpArr> GrpArrList; 
typedef std::vector<unsigned long> GroupType;
typedef matrix::Matrix<unsigned long> Endomorphism;

class Homomorphism;

class FiniteAbelianGroup;
class GeneralFiniteAbelianGroup;
class ElementaryTwoGroup;

class SubgroupIterator;

}

}

#endif
