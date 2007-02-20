/*
  This is wgraph.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "wgraph.h"
#include "blocks.h"

#include <algorithm>
#include <iostream>

namespace atlas {

/*****************************************************************************

        Chapter I -- The WGraph and DecomposedWGraph classes

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

DecomposedWGraph::DecomposedWGraph(const WGraph& wg)
  : d_cell(), d_part(wg.size()), d_id(), d_induced()
{
  using partition::Partition; using partition::PartitionIterator;

  Partition pi;
  wg.cells(pi,&d_induced);

  d_cell.reserve(pi.classCount()); // there will be this many cells
  d_id.resize(pi.classCount());    // and vectors of identification numbers


  std::vector<unsigned int> relno(wg.size()); // number of element in its cell
  PartitionIterator i(pi);
  for (cell_no n=0; n<d_id.size(); ++n,++i)
  {
    d_cell.push_back(WGraph(wg.rank()));
    WGraph& cur_cell = d_cell.back();
    cur_cell.resize(i->second-i->first);
    std::vector<unsigned int>& idn=d_id[n];
    idn.resize(cur_cell.size());

    for (PartitionIterator::SubIterator j=i->first; j!=i->second; ++j)
    {
      size_t y = *j; size_t z=j-i->first; // |y| gets renamed |z| in cell

      d_part[y]=n; // or equivalently |d_part[y]=pi(y)|
      relno[y]=z; idn[z]=y;

      cur_cell.descent(z) = wg.descent(y);
    }
  }
  // make sure all values |relno[y]| are defined before proceeding

  for (i.rewind(); i(); ++i)
  {
    cell_no n=d_part[*i->first]; // cell number, constant through next loop
    for (PartitionIterator::SubIterator j=i->first; j!=i->second; ++j)
    {
      size_t y = *j; size_t z=relno[y]; // |z==j-i->first|
      const graph::EdgeList& edge = wg.edgeList(y);
      const CoeffList& coeff = wg.coeffList(y);
      graph::EdgeList& cur_el = d_cell[n].edgeList(z);
      CoeffList& cur_cl = d_cell[n].coeffList(z);
      for (size_t k = 0; k < edge.size(); ++k)
	if (d_part[edge[k]]==n) // only look at edges within this cell
	{
	  cur_el.push_back(relno[edge[k]]);
	  cur_cl.push_back(     coeff[k]);
	}
    } // for (j)
  } // for (i)

}

} // namespace wgraph

/*****************************************************************************

        Chapter II -- Functions declared in wgraph.h

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
  wg.cells(pi); // do not collect information about induced graph here

  for (PartitionIterator i(pi); i(); ++i) {
    wc.push_back(WGraph(wg.rank()));
    WGraph& wci = wc.back();
    wci.resize(i->second - i->first);

    /* looping over |z| rather than using |*i| directly implements the
       renumbering of each cell (what was |y=*(i->first+z)| becomes just |z|)
    */
    for (size_t z = 0; z < wci.size(); ++z) {
      size_t y = i->first[z];
      wci.descent(z) = wg.descent(y);
      const EdgeList& el = wg.edgeList(y);
      EdgeList& eli = wci.edgeList(z);
      const CoeffList& cl = wg.coeffList(y);
      CoeffList& cli = wci.coeffList(z);
      for (size_t j = 0; j < el.size(); ++j) {
	size_t x = el[j];
	if (pi(x) != pi(y)) // an edge leading out of the current cell
	  continue;         // ignore these
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
