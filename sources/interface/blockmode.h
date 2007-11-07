/*
  This is realmode.h

  Copyright (C) 2007 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef BLOCKMODE_H  /* guard against multiple inclusions */
#define BLOCKMODE_H

#include "commands_fwd.h"
#include "complexredgp_fwd.h"
#include "realredgp_fwd.h"
#include "realform.h"
#include "blocks_fwd.h"
#include "kl_fwd.h"
#include "wgraph.h"

namespace atlas {

namespace blockmode {

/******** type declarations *************************************************/

  struct BlockmodeTag {};

/******** function declarations ********************************************/

  commands::CommandMode& blockMode();
  complexredgp::ComplexReductiveGroup& currentDualComplexGroup();
  realredgp::RealReductiveGroup& currentDualRealGroup();
  realform::RealForm currentDualRealForm();
  blocks::Block& currentBlock();
  kl::KLContext& currentKL();
  wgraph::WGraph& currentWGraph();

  void addBlockHelp(commands::CommandMode&, commands::TagDict&);


} // namespace blockmode


} // namespace atlas

#endif
