/*
  This is helpmode.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef HELPMODE_H  /* guard against multiple inclusions */
#define HELPMODE_H

#include "commands_fwd.h"

/******** function declarations ********************************************/

namespace atlas {

namespace helpmode {

  const commands::CommandMode& helpMode();
  void intro_h();
  void nohelp_h();

}

}

#endif
