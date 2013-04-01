/*
  This is helpmode.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#ifndef HELPMODE_H  /* guard against multiple inclusions */
#define HELPMODE_H

#include "commands_fwd.h"
#include <iostream>

namespace atlas {

namespace commands {

/******** function and variable declarations ********************************/

  commands::CommandNode helpNode();
  extern commands::CommandTree help_mode; // defined in commands.cpp

  void intro_h(); // this is used in emptymode as well

}

}

#endif
