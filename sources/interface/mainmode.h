/*
  This is mainmode.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#ifndef MAINMODE_H  /* guard against multiple inclusions */
#define MAINMODE_H

#include "../Atlas.h"
#include "commands_fwd.h"

namespace atlas {

namespace commands {

/******** type declarations ************************************************/

  struct MainmodeTag {};

/******** function and variable declarations *******************************/

  commands::CommandNode mainNode();
  extern commands::CommandTree& main_mode; // defined in commands.cpp
  InnerClass& currentComplexGroup();
  InnerClass& current_dual_group();
  const lietype::Layout& current_layout();
  const WeightList& current_lattice_basis();
  output::Interface& currentComplexInterface();
  void replaceComplexGroup(InnerClass*,output::Interface*);

} // |namespace commands|

} // |namespace atlas|

#endif
