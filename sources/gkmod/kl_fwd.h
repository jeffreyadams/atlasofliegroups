/*
  This is kl_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

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

typedef unsigned KLCoeff;
typedef polynomials::Polynomial<KLCoeff> KLPol;

typedef std::set<KLPol>::iterator KLPtr;
typedef std::vector<KLPtr> KLRow;

typedef unsigned MuCoeff;
typedef std::pair<size_t,MuCoeff> MuData;
typedef std::vector<MuData> MuRow;

}

}

#endif
