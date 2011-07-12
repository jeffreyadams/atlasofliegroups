/*!
\file
  This is graph_fwd.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef GRAPH_FWD_H  /* guard against multiple inclusions */
#define GRAPH_FWD_H

#include "set.h" // that is really a forward-declaration only header file

/******** type declarations   ******************************************/

namespace atlas {

namespace graph {

  typedef set::SetElt Vertex;
  typedef std::vector<Vertex> VertexList;
  typedef Vertex Edge;
  typedef std::vector<Edge> EdgeList;

  class OrientedGraph;

} // |namespace graph|

} // |namespace atlas|

#endif
