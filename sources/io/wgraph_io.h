/*
  This is wgraph_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef WGRAPH_IO_H  /* guard against multiple inclusions */
#define WGRAPH_IO_H

#include <iosfwd>

#include "../Atlas.h"

namespace atlas {

/******** function declarations *********************************************/

namespace wgraph_io {

  void printWGraph(std::ostream&, const wgraph::WGraph&);

  void printCells(std::ostream&, const wgraph::WGraph&);

  void printWDecomposition (std::ostream&, const wgraph::DecomposedWGraph&);

}

}

#endif
