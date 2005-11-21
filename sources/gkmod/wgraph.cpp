/*
  This is wgraph.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "wgraph.h"

#include <algorithm>

namespace atlas {

/*****************************************************************************

        Chapter I -- The WGraph class

  ... explain here when it is stable ...

******************************************************************************/

namespace wgraph {

/******** constructors and destructors ***************************************/
WGraph::WGraph(size_t r)
  :d_rank(r)

/*
  Synopsis: constructor

  Makes an empty W-graph with rank r.
*/

{}

WGraph::~WGraph()

{}

/******** copy, assignment and swap ******************************************/
void WGraph::swap(WGraph& other)

{
  std::swap(d_rank,other.d_rank);
  d_graph.swap(other.d_graph);
  d_coeff.swap(other.d_coeff);
  d_descent.swap(other.d_descent);

  return;
}

/******** manipulators *******************************************************/
void WGraph::reset()

/*
  Synopsis: resets the W-graph

  Preserves size and rank; resets edge and coefficient lists to empty lists;
  resets descent sets to empty.
*/

{
  d_graph.reset();
  d_coeff.assign(size(),CoeffList());
  d_descent.assign(size(),bitset::RankFlags());

  return;
}

void WGraph::resize(size_t n)

/*
  Synopsis: resizes the lists to size n

  Preserves the existing parts if n > size().
*/

{
  d_graph.resize(n);
  d_coeff.resize(n);
  d_descent.resize(n);

  return;
}

}

/*****************************************************************************

        Chapter II -- Functions declared in wgraph.h

  ... explain here when it is stable ...

******************************************************************************/

namespace wgraph {

void cells(std::vector<WGraph>& wc, const WGraph& wg)

/*
  Synopsis: puts in wc the cells of the W-graph wg.
*/

{
  using namespace graph;
  using namespace partition;

  Partition pi;
  wg.cells(pi);

  for (PartitionIterator i(pi); i(); ++i) {
    wc.push_back(WGraph(wg.rank()));
    WGraph& wci = wc.back();
    wci.resize(i->second - i->first);
    for (size_t z = 0; z < wci.size(); ++z) {
      size_t y = i->first[z];
      wci.descent(z) = wg.descent(y);
      const EdgeList& el = wg.edgeList(y);
      EdgeList& eli = wci.edgeList(z);
      const CoeffList& cl = wg.coeffList(y);
      CoeffList& cli = wci.coeffList(z);
      for (size_t j = 0; j < el.size(); ++j) {
	size_t x = el[j];
	if (pi(x) != pi(y))
	  continue;
	// find relative address of x in this class
	size_t xi = std::lower_bound(i->first,i->second,x) - i->first;
	eli.push_back(xi);
	cli.push_back(cl[j]);
      }
    }
  }

  return;
}

}

}
