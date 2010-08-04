/*!
\file
\brief Class declarations and type definitions for TitsGroup
*/
/*
  This is tits_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef TITS_FWD_H  /* guard against multiple inclusions */
#define TITS_FWD_H

#include <vector>

#include "constants.h"

/******** forward type declarations ******************************************/

namespace atlas {

namespace tits {

  class TorusElement;
  class GlobalTitsElement;
  class GlobalTitsGroup;
  class TitsElt;
  class TitsGroup;

  typedef std::vector<TitsElt> TitsEltList;

}

}

#endif
