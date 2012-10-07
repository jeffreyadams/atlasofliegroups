/*
  This is mainmode.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#ifndef MAINMODE_H  /* guard against multiple inclusions */
#define MAINMODE_H

#include "atlas_types.h"
#include "commands_fwd.h"

/******** type declarations ************************************************/

namespace atlas {

namespace mainmode {

  struct MainmodeTag {};

}

/******** function declarations ********************************************/

namespace mainmode {

  commands::CommandMode& mainMode();
  ComplexReductiveGroup& currentComplexGroup();
  ComplexReductiveGroup& current_dual_group();
  complexredgp_io::Interface& currentComplexInterface();
  void replaceComplexGroup(ComplexReductiveGroup*
			   ,complexredgp_io::Interface*);

}

}

#endif
