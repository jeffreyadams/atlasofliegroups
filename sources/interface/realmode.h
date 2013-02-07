/*
  This is realmode.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#ifndef REALMODE_H  /* guard against multiple inclusions */
#define REALMODE_H

#include "commands_fwd.h"
#include "atlas_types.h"

namespace atlas {

namespace commands {

/******** type declarations *************************************************/

  struct RealmodeTag {};

/******** function declarations ********************************************/

  commands::CommandNode realNode();
  extern commands::CommandTree& real_mode; // defined in main.cpp
  RealReductiveGroup& currentRealGroup();
  RealFormNbr currentRealForm();
  const Rep_context& currentRepContext();
  Rep_table& currentRepTable();

}

}

#endif
