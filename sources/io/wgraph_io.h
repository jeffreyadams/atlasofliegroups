/*
  This is wgraph_io.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef WGRAPH_IO_H  /* guard against multiple inclusions */
#define WGRAPH_IO_H

#include <iosfwd>

#include "wgraph.h"

namespace atlas {

/******** function declarations *********************************************/

namespace wgraph_io {

  std::ostream& printCells(std::ostream&, const wgraph::WGraph&);

  std::ostream& printWGraph(std::ostream&, const wgraph::WGraph&);

}

}

#endif
