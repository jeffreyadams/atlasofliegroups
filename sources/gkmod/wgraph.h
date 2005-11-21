/*
  This is wgraph.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef WGRAPH_H  /* guard against multiple inclusions */
#define WGRAPH_H

#include "bitset.h"
#include "graph.h"
#include "partition.h"

namespace atlas {

/******** type declarations *************************************************/

namespace wgraph {

class WGraph;

 typedef unsigned short Coeff;
 typedef std::vector<Coeff> CoeffList;

}

/******** function declarations *********************************************/

namespace wgraph {

void cells(std::vector<WGraph>&, const WGraph&);

}

/******** type definitions **************************************************/

namespace wgraph {

class WGraph {

 private:

  size_t d_rank;
  graph::OrientedGraph d_graph;
  std::vector<CoeffList> d_coeff;
  std::vector<bitset::RankFlags> d_descent;

 public:

// constructors and destructors
  explicit WGraph(size_t);

  ~WGraph();

// copy, assignment and swap
  void swap(WGraph&);

// accessors
  void cells(partition::Partition& pi, graph::OrientedGraph* p = 0) const {
    d_graph.cells(pi,p);
  }

  const CoeffList& coeffList(graph::Vertex x) const {
    return d_coeff[x];
  }

  const bitset::RankFlags& descent(graph::Vertex x) const {
    return d_descent[x];
  }

  const graph::EdgeList& edgeList(graph::Vertex x) const {
    return d_graph.edgeList(x);
  }

  const graph::OrientedGraph& graph() const {
    return d_graph;
  }

  const size_t rank() const {
    return d_rank;
  }

  size_t size() const {
    return d_graph.size();
  }

// manipulators
  CoeffList& coeffList(graph::Vertex x) {
    return d_coeff[x];
  }

  bitset::RankFlags& descent(graph::Vertex x) {
    return d_descent[x];
  }

  graph::EdgeList& edgeList(graph::Vertex x) {
    return d_graph.edgeList(x);
  }

  void reset();

  void resize(size_t);
};

}

}

#endif
