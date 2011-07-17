/*
  This is realmode.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For copyright and license information see the LICENSE file
*/

#ifndef REALMODE_H  /* guard against multiple inclusions */
#define REALMODE_H

#include "commands_fwd.h"
#include "atlas_types.h"
#include "realredgp_io_fwd.h"

/******** type declarations *************************************************/

namespace atlas {

namespace realmode {

  struct RealmodeTag {};

}

/******** function declarations ********************************************/

namespace realmode {

  commands::CommandMode& realMode();
  RealReductiveGroup& currentRealGroup();
  RealFormNbr currentRealForm();

}

}

#endif
