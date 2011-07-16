/*!
\file
\brief Forward declaration for the class KLContext and associated types.
*/
/*
  This is kl_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef KL_FWD_H  /* guard against multiple inclusions */
#define KL_FWD_H

#include <set>
#include <vector>

#include "polynomials_fwd.h"
#include "atlas_types.h"

namespace atlas {

/******** type declarations *************************************************/

namespace kl {

class KLContext;

/*!
\brief Coefficient of a KL polynomial.
*/
typedef unsigned int KLCoeff;

/*!
\brief Polynomial with coefficients of type KLCoeff.
*/
typedef polynomials::Safe_Poly<KLCoeff> KLPol;

typedef unsigned int KLIndex; // less than 2^32 distinct polynomials for E8 !


typedef KLCoeff MuCoeff;

// a pair of vectors may be much more compact that a vector of pairs
typedef std::pair<std::vector<BlockElt>,std::vector<MuCoeff> > MuRow;

typedef std::vector<KLPol> KLStore;

typedef KLStore::const_reference KLPolRef;

typedef std::vector<KLIndex> KLRow;

typedef std::vector<BlockElt> PrimitiveRow;

}

}

#endif
