/*
  This is mainhelp.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#ifndef MAINHELP_H  /* guard against multiple inclusions */
#define MAINHELP_H

#include "commands_fwd.h"

/******** function declarations ********************************************/

namespace atlas {

namespace mainhelp {

  void addMainHelp(commands::CommandMode&, commands::TagDict&);

}

}

#endif
