/*!
\file
\brief Class declaration and type definitions for class Block.
*/

/*
  This is blocks_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef BLOCKS_FWD_H  /* guard against multiple inclusions */
#define BLOCKS_FWD_H

#include <cstddef>
#include <vector>

namespace atlas {

/******** type declarations *************************************************/

namespace blocks {

  class Block_base;
  class Block;

  typedef unsigned int BlockElt;
  typedef std::vector<BlockElt> BlockEltList;
  typedef std::pair<BlockElt,BlockElt> BlockEltPair;
  typedef std::vector<BlockEltPair> BlockEltPairList;

}

}

#endif
