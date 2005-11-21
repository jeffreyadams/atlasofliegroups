/*
  This is block_io.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef BLOCK_IO_H  /* guard against multiple inclusions */
#define BLOCK_IO_H

#include <iosfwd>

#include "blocks_fwd.h"
#include "descents_fwd.h"

namespace atlas {

/******** function declarations *********************************************/

namespace block_io {

  std::ostream& printBlock(std::ostream&, const blocks::Block&);

  std::ostream& printDescent(std::ostream&, const descents::DescentStatus&,
			     size_t);

}

}

#endif
