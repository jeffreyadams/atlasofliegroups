/*!
\file
  This is graph.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef GRAPH_H  /* guard against multiple inclusions */
#define GRAPH_H

#include "graph_fwd.h"
#include "partition.h"
#include "set.h"

/******** type definitions **************************************************/

namespace atlas {

/******** type declarations *************************************************/

namespace graph {

class OrientedGraph
{

  std::vector<EdgeList> d_edges; // list of edges, outgoing from each vertex

 public:

// constructors and destructors
  OrientedGraph() {}

  explicit OrientedGraph(size_t n):d_edges(n) {}

// copy, construction and swap
  void swap(OrientedGraph& other) { d_edges.swap(other.d_edges); }

// accessors
  void cells(partition::Partition&, OrientedGraph* p = 0) const;

  Vertex edge(Vertex x, size_t j) const { return d_edges[x][j]; } // edge #j

  const EdgeList& edgeList(const Vertex& x) const { return d_edges[x]; }

  size_t size() const { return d_edges.size(); } // number of vertices

// manipulators
  Vertex& edge(Vertex x, size_t j) { return d_edges[x][j]; } // allow clobbering

  EdgeList& edgeList(const Vertex& x) { return d_edges[x]; } // even globally

  Vertex newVertex()
    { Vertex v=size(); d_edges.push_back(EdgeList()); return v; }

  void reset() { d_edges.assign(d_edges.size(),EdgeList()); } // clear edges

  void resize(size_t n) { d_edges.resize(n); } // change number of vertices

  void reverseEdges ();      // make opposite oriented graph
  void reverseNumbering ();  // same graph, but reverse numbering of vertices


  // auxiliary methods
private:
void addLinks
  (const std::vector<const EdgeList*>& out, const partition::Partition& pi);
}; // |class OrientedGraph|

} // |namespace graph|

} // |namespace atlas|

#endif
