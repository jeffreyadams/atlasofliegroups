/*!
\file
  This is graph.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#ifndef GRAPH_H  /* guard against multiple inclusions */
#define GRAPH_H

#include "partition.h"
#include "set.h"

/******** type definitions **************************************************/

namespace atlas {

namespace graph {

  typedef set::SetElt Vertex;
  typedef std::vector<Vertex> VertexList;
  typedef Vertex Edge;
  typedef std::vector<Edge> EdgeList;

  class OrientedGraph;

}

/******** type declarations *************************************************/

namespace graph {

class OrientedGraph {

 private:

  std::vector<EdgeList> d_edges;

 public:

// constructors and destructors
  OrientedGraph() {}

  explicit OrientedGraph(size_t n):d_edges(n) {}

  ~OrientedGraph() {}

// copy, construction and swap
  void swap(OrientedGraph& other) {
    d_edges.swap(other.d_edges);
  }

// accessors
  void cells(partition::Partition&, OrientedGraph* p = 0) const;

  Vertex edge(Vertex x, size_t j) const {
    return d_edges[x][j];
  }

  const EdgeList& edgeList(const Vertex& x) const {
    return d_edges[x];
  }

  size_t size() const {
    return d_edges.size();
  }

// manipulators
  Vertex& edge(Vertex x, size_t j) {
    return d_edges[x][j];
  }

  EdgeList& edgeList(const Vertex& x) {
    return d_edges[x];
  }

  void reset() {
    d_edges.assign(d_edges.size(),EdgeList());
  }

  void resize(size_t n) {
    d_edges.resize(n);
  }
};

}

}

#endif
