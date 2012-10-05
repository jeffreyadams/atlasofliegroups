/*
  This is mainhelp.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations 

  For copyright and license information see the LICENSE file
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
