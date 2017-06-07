/*
  This is wgraph_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2007 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "wgraph_io.h"

#include <iostream>
#include <cassert>

#include "wgraph.h"	// |WGraph|
#include "prettyprint.h" // |printDescentSet|

namespace atlas {

/*****************************************************************************

        Chapter I -- Functions declared in wgraph_io.h

******************************************************************************/

namespace wgraph_io {


/*
  Output the structural data for the W-graph to strm.

  The output format is as follows: one line peer vertex, consisting of three
  colon-separated fields:

    - the vertex number (starting from 0);
    - the descent set (i.e., the tau-invariant), written as {*,*,...}
      (generator numbers starting from 1);
    - the edges originating from that vertex, written as {(*,*),(*,*),...}
      where each (*,*) is a pair (vertex,label); note that the vast majority
      of labels are 1, but not all of them in general;

  This format should be human-readable, but is purposefully kept simple so
  that it can be easily parsed by another program (for instance, a
  prettyprinter?)
*/
void printWGraph(std::ostream& strm, const wgraph::WGraph& wg)
{
  for (size_t z = 0; z < wg.size(); ++z)
  {
    strm << z << ":";
    prettyprint::printDescentSet(strm,wg.descent_set(z),wg.rank());
    assert(wg.coefficients.size()==wg.degree(z));
    strm << ":{";
    for (size_t j = 0; j < wg.degree(z); ++j)
    {
      if (j>0)
	strm << ",";
      strm << "(" << wg.edge_target(z,j) << "," << wg.coefficient(z,j) << ")";
    }
    strm << "}" << std::endl;
  }
} // |printWGraph|

//  Output the various cells in wg, as W-graphs, in the format of printWGraph.
void printCells(std::ostream& strm, const wgraph::WGraph& wg)
{
  auto wc = wgraph::cells(wg);
  for (size_t j = 0; j < wc.size(); ++j)
    printWGraph(strm << "// cell #" << j << std::endl,wc[j]);
}

/*
  We try to improve the output provided by |printCells| above, by providing
  the original |BlockElt| numbers of the vertices, which are lost when
  |wgraph::cells| is called. The |DecomposedWGraph| class was introduced to
  avoid this. At the same time we have occasion to print the induced graph.
*/
void printWDecomposition(std::ostream& strm, const wgraph::DecomposedWGraph& g)
{
  strm << "// Cells and their vertices.\n";
  for (size_t i = 0; i < g.cellCount(); ++i)
    {
      strm << '#' << i << "=";
      std::vector<BlockElt> mem=g.cellMembers(i);
      for (size_t j=0; j<mem.size(); ++j)
	strm << (j==0 ? '{' : ',') << mem[j];
      strm << "}\n";
    }


  strm << "\n// Induced graph on cells.\n";
  for (size_t i = 0; i < g.inducedGraph().size(); ++i)
  {
    strm << '#' << i << ':';
    const graph::EdgeList& el = g.inducedGraph().edgeList(i);
    for (size_t j=0; j<el.size(); ++j)
      strm << (j==0 ? "->#" : ",#") << el[j];
    strm << ".\n";
  }

  strm << "\n// Individual cells.\n";
  for (size_t i = 0; i < g.cellCount(); ++i)
    {
      strm << "// cell #" << i << ":\n";
      wgraph::WGraph cell=g.cell(i);
      std::vector<BlockElt> mem=g.cellMembers(i);
      for (size_t j=0; j<cell.size(); ++j)
      {
	strm << j << '[' << mem[j] << "]: ";
	prettyprint::printDescentSet(strm,cell.descent_set(j),g.rank());
	for (size_t k=0; k<cell.degree(j); ++k)
	{
	  strm << (k==0 ? " --> " : ",");
	  if (cell.coefficient(j,k)==1) strm << cell.edge_target(j,k);
	  else strm << '(' << cell.edge_target(j,k)
		    << ',' << cell.coefficient(j,k) << ')';
	} // for (k)
	strm << "\n"; // end line
      } // for (j)
      strm << "\n"; // empty line between cells
    } // for (i)
} // |printWDecomposition|

}

}
