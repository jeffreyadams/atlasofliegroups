/*
  This is realhelp.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  For copyright and license information see the LICENSE file
*/

#ifndef REALHELP_H  /* guard against multiple inclusions */
#define REALHELP_H

#include "commands_fwd.h"

/******** function declarations ********************************************/

namespace atlas {

namespace realhelp {

  void addRealHelp(commands::CommandMode&, commands::TagDict&);

}

}

#endif
