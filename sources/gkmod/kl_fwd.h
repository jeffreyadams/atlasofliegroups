/*!
\file
\brief Forward declaration for the class KLContext and associated types.
*/
/*
  This is kl_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups  

  See file main.cpp for full copyright notice
*/

#ifndef KL_FWD_H  /* guard against multiple inclusions */
#define KL_FWD_H

#include <set>
#include <vector>

#include "polynomials_fwd.h"

namespace atlas {

/******** type declarations *************************************************/

namespace kl {

class KLContext;

/*!
\brief Coefficient of a KL polynomial.

Must be a standard unsigned type; now "unsigned" which [at least on my
Mac; not sure what the standard says - DV 8/14/06] is unsigned long.
*/
typedef unsigned KLCoeff;

/*!
\brief Polynomial with coefficients of type KLCoeff.
*/
typedef polynomials::Polynomial<KLCoeff> KLPol;

/*!
\brief Pointer into a set of KLPol's.
*/
typedef std::set<KLPol>::iterator KLPtr;

typedef std::vector<KLPtr> KLRow;

typedef unsigned MuCoeff;
typedef std::pair<size_t,MuCoeff> MuData;
typedef std::vector<MuData> MuRow;

}

}

#endif
