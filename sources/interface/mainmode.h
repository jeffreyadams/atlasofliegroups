/*
  This is mainmode.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef MAINMODE_H  /* guard against multiple inclusions */
#define MAINMODE_H

#include "commands_fwd.h"
#include "complexredgp_fwd.h"
#include "complexredgp_io_fwd.h"

/******** type declarations ************************************************/

namespace atlas {

namespace mainmode {

  struct MainmodeTag {};

}

/******** function declarations ********************************************/

namespace mainmode {

  const commands::CommandMode& mainMode();
  complexredgp::ComplexReductiveGroup& currentComplexGroup();
  complexredgp_io::Interface& currentComplexInterface();

}

}

#endif
