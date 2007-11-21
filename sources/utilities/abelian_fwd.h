/*!
\file
  This is abelian_fwd.h
*/
/*  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  For license information see the LICENSE file
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
 class GeneralFiniteAbelianGroup; //not implemented
 class ElementaryTwoGroup; //not implemented

class SubgroupIterator;

}

}

#endif
