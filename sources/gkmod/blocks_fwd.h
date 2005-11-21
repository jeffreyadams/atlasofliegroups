/*
  This is blocks_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef BLOCKS_FWD_H  /* guard against multiple inclusions */
#define BLOCKS_FWD_H

#include <cstddef>
#include <vector>

namespace atlas {

/******** type declarations *************************************************/

namespace blocks {

  class Block;

  typedef size_t BlockElt;
  typedef std::vector<BlockElt> BlockEltList;
  typedef std::pair<BlockElt,BlockElt> BlockEltPair;
  typedef std::vector<BlockEltPair> BlockEltPairList;

}

}

#endif
