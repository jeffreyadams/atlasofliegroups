/*!
\file
  This is abelian_fwd.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef ABELIAN_FWD_H  /* guard against multiple inclusions */
#define ABELIAN_FWD_H

#include <limits>
#include <vector>
#include <set>

#include "matrix_fwd.h"

/******** type declarations **************************************************/

namespace atlas {

namespace abelian {

typedef unsigned long long GrpNbr;  //!< element of Abelian group, compact repr.
typedef std::vector<GrpNbr> GrpNbrList;
typedef std::vector<unsigned long> GrpArr; //!< group element, array repr.
typedef std::vector<GrpArr> GrpArrList;
typedef std::vector<unsigned long> GroupType;
typedef matrix::Matrix_base<unsigned long> Endomorphism;

class Homomorphism;

class FiniteAbelianGroup;
class GeneralFiniteAbelianGroup; //not implemented
class ElementaryTwoGroup; //not implemented

class SubgroupIterator;

} // |namespace abelian|

} // |namespace atlas|

#endif
