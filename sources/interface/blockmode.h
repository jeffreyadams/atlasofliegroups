/*
  This is blockmode.h

  Copyright (C) 2007 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#ifndef BLOCKMODE_H  /* guard against multiple inclusions */
#define BLOCKMODE_H

#include "commands_fwd.h"
#include "atlas_types.h"
#include "wgraph.h"

namespace atlas {

namespace blockmode {

/******** type declarations *************************************************/

  struct BlockmodeTag {};

/******** function declarations ********************************************/

  commands::CommandMode& blockMode();
  ComplexReductiveGroup& currentDualComplexGroup();
  RealReductiveGroup& currentDualRealGroup();
  RealFormNbr currentDualRealForm();
  Block& currentBlock();
  kl::KLContext& currentKL();
  wgraph::WGraph& currentWGraph();

  void addBlockHelp(commands::CommandMode&, commands::TagDict&);


} // namespace blockmode


} // namespace atlas

#endif
