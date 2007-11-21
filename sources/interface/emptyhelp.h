/*
  This is emptyhelp.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  For copyright and license information see the LICENSE file
*/

#ifndef EMPTYHELP_H  /* guard against multiple inclusions */
#define EMPTYHELP_H

#include "commands_fwd.h"

/******** function declarations ********************************************/

namespace atlas {

namespace emptyhelp {

  void addEmptyHelp(commands::CommandMode&, commands::TagDict&);

}

}

#endif
