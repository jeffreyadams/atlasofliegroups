/*
  This is emptymode.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
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

  const commands::CommandMode& emptyMode();

}

}

#endif
