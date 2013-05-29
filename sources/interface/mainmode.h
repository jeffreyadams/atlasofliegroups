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

namespace atlas {

namespace commands {

/******** type declarations ************************************************/

  struct MainmodeTag {};

/******** function and variable declarations *******************************/

  commands::CommandNode mainNode();
  extern commands::CommandTree& main_mode; // defined in commands.cpp
  ComplexReductiveGroup& currentComplexGroup();
  ComplexReductiveGroup& current_dual_group();
  complexredgp_io::Interface& currentComplexInterface();
  void replaceComplexGroup(ComplexReductiveGroup*
			   ,complexredgp_io::Interface*);

} // |namespace commands|

} // |namespace atlas|

#endif
