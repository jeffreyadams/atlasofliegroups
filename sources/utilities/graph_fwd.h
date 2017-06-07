/*
  This is graph_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef GRAPH_FWD_H  /* guard against multiple inclusions */
#define GRAPH_FWD_H

#include <vector>

/******** type declarations   ******************************************/

namespace atlas {

namespace graph {

  typedef unsigned int Vertex; // assume at most some 4 billion vertices
  typedef std::vector<Vertex> EdgeList; // list of targets of outgoing edges

  class OrientedGraph;

} // |namespace graph|

} // |namespace atlas|

#endif
