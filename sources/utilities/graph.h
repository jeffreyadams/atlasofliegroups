/*
  This is graph.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef GRAPH_H  /* guard against multiple inclusions */
#define GRAPH_H

#include "set.h"

/******** type definitions **************************************************/

namespace atlas {

namespace graph {

  typedef set::SetElt Vertex;
  typedef std::vector<Vertex> VertexList;

  class OrientedGraph;

}

/******** type declarations *************************************************/

namespace graph {

class OrientedGraph {

 private:

  std::vector<VertexList> d_edges;

 public:

/* constructors and destructors */
  explicit OrientedGraph(size_t n):d_edges(n) {}

  ~OrientedGraph() {}

/* accessors */
  Vertex edge(Vertex x, size_t j) const {
    return d_edges[x][j];
  }

  const VertexList& edges(const Vertex& x) const {
    return d_edges[x];
  }

  size_t size() const {
    return d_edges.size();
  }

/* manipulators */
  Vertex& edge(Vertex x, size_t j) {
    return d_edges[x][j];
  }

  VertexList& edges(const Vertex& x) {
    return d_edges[x];
  }

  void resize(size_t n) {
    d_edges.resize(n);
  }
};

}

}

#endif
