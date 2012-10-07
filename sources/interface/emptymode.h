/*
  This is emptymode.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#ifndef EMPTYMODE_H  /* guard against multiple inclusions */
#define EMPTYMODE_H

#include "commands_fwd.h"

/******** type declarations *************************************************/

namespace atlas {

namespace emptymode {

  struct EmptymodeTag {};

}

/******** function declarations ********************************************/

namespace emptymode {

commands::CommandMode& emptyMode();

}

}

#endif
