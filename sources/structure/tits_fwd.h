/*!
\file
\brief Class declarations and type definitions for TitsGroup
*/
/*
  This is tits_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#ifndef TITS_FWD_H  /* guard against multiple inclusions */
#define TITS_FWD_H

#include <vector>

#include "constants.h"

/******** forward type declarations ******************************************/

namespace atlas {

namespace tits {

  class TitsElt;
  class TitsGroup;
  class TwistedTitsGroup; //not implemented

  typedef std::vector<TitsElt> TitsEltList;

}

}

#endif
