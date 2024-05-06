/*
  This is reprmode.h

  Copyright (C) 2007 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#ifndef REPRMODE_H  /* guard against multiple inclusions */
#define REPRMODE_H

#include "commands_fwd.h"
#include "Atlas.h"
#include "wgraph.h"

namespace atlas {

namespace commands {

/******** type declarations ************************************************/

  struct ReprmodeTag {};
  enum block_type { noblock, full_block, partial_block };

/******** function and variable declarations ********************************/

  commands::CommandNode reprNode();
  extern CommandTree& repr_mode; // defined in commands.cpp
  const SubSystemWithGroup& currentSubSystem();
  blocks::common_block& current_common_block();
  const StandardRepr& currentStandardRepr();
  void ensure_full_block();
  kl::KL_table& current_KL();
  const wgraph::WGraph& current_WGraph();
} // |namespace reprmode|

} // |namespace atlas|

#endif
