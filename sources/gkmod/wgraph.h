/*
  This is wgraph.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2007,2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef WGRAPH_H  /* guard against multiple inclusions */
#define WGRAPH_H

#include <iostream>

#include "../Atlas.h"   // must be included before any utility headers are

#include "bitset.h"	// inlines
#include "graph.h"	// containment


namespace atlas {

/******** type declarations  (see ../Atlas.h)  ***************************/


/******** function declarations *********************************************/

namespace wgraph {

std::vector<WGraph> cells(const WGraph&);

// Functions

WGraph wGraph
  ( std::ifstream& block_file
  , std::ifstream& matrix_file
  , std::ifstream& KL_file);

}

/******** type definitions **************************************************/

namespace wgraph {

/*
  The |WGraph| structure stores a graph, edge weights (coefficients), and
  descent sets. The construction of the graph is left to the client.
*/
struct WGraph
{
  typedef unsigned short Coeff;
  typedef std::vector<Coeff> CoeffList;

  size_t d_rank;
  graph::OrientedGraph oriented_graph;
  std::vector<CoeffList> coefficients;
  std::vector<RankFlags> descent_sets;

// constructors and destructors
  explicit WGraph(size_t rank, size_t size);

// copy, assignment and swap: nothing needed beyond defaults

// accessors
  Coeff coefficient(graph::Vertex x,unsigned int index) const
    { return coefficients[x][index]; }

  const RankFlags& descent_set(graph::Vertex x) const
    { return descent_sets[x]; }

  unsigned int degree (graph::Vertex x) const
    { return oriented_graph.edgeList(x).size(); }

  graph::Vertex edge_target (graph::Vertex x,unsigned int index) const
  { return oriented_graph.edgeList(x)[index]; }

  const graph::OrientedGraph& graph() const { return oriented_graph; }

  const size_t rank() const { return d_rank; }

  size_t size() const { return oriented_graph.size(); }

}; // |class WGraph|

class DecomposedWGraph
{
  typedef unsigned int cell_no;

  std::vector<WGraph> d_cell; // the strong components

  std::vector<cell_no> d_part;    // assigns strong component to each BlockElt
  std::vector< std::vector<BlockElt> > d_id; // original vertex numbers

  graph::OrientedGraph d_induced; // induced graph on cells

 public:

// constructors and destructors
  explicit DecomposedWGraph(const WGraph& wg);

// copy, assignment and swap: nothing needed beyond defaults

// accessors
  size_t rank () const { return d_cell[0].rank(); } // all ranks are equal
  size_t cellCount() const { return d_cell.size(); }
  const graph::OrientedGraph& inducedGraph() const { return d_induced; }
  const wgraph::WGraph& cell (size_t c) const { return d_cell[c]; }
  const std::vector<BlockElt>& cellMembers(size_t c) const
    { return d_id[c]; }

}; // |class DecomposedWGraph|

} // |namespace wgraph|

} // |namespace atlas|

#endif
