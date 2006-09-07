/*
  This is abelian2_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#ifndef ABELIAN2_FWD_H  /* guard against multiple inclusions */
#define ABELIAN2_FWD_H

#include <limits>
#include <set>

#include "latticetypes_fwd.h"
#include "matrix_fwd.h"

/******** type declarations **************************************************/

namespace atlas {

namespace abelian2 {

typedef unsigned long GrpNbr;
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
