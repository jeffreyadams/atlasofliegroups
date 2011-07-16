/*
  This is wgraph_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2007 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "wgraph_io.h"

#include <iostream>

#include "prettyprint.h"
#include "blocks_fwd.h"

namespace atlas {

/*****************************************************************************

        Chapter I -- Functions declared in wgraph_io.h

******************************************************************************/

namespace wgraph_io {


void printWGraph(std::ostream& strm, const wgraph::WGraph& wg)

/*
  Synopsis: outputs the structural data for the W-graph to strm.

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

{
  using namespace graph;
  using namespace prettyprint;
  using namespace wgraph;

  for (size_t z = 0; z < wg.size(); ++z) {
    strm << z << ":";
    printDescentSet(strm,wg.descent(z),wg.rank());
    const EdgeList& el = wg.edgeList(z);
    const WCoeffList& cl = wg.coeffList(z);
    // el and cl have the same size
    strm << ":{";
    for (size_t j = 0; j < el.size(); ++j) {
      if (j)
	strm << ",";
      strm << "(" << el[j] << "," << cl[j] << ")";
    }
    strm << "}";
    strm << std::endl;
  }
}

void printCells(std::ostream& strm, const wgraph::WGraph& wg)

/*
  Synopsis: outputs the various cells in wg, as W-graphs, in the format
  of printWGraph.
*/

{
  std::vector<wgraph::WGraph> wc;
  wgraph::cells(wc,wg);

  for (size_t j = 0; j < wc.size(); ++j) {
    strm << "// cell #" << j << std::endl;
    printWGraph(strm,wc[j]);
  }
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
	prettyprint::printDescentSet(strm,cell.descent(j),g.rank());
	const graph::EdgeList& el = cell.edgeList(j);
	const wgraph::WCoeffList& cl = cell.coeffList(j);
	for (size_t k=0; k<el.size(); ++k)
	{
	  strm << (k==0 ? " --> " : ",");
	  if (cl[k]==1) strm << el[k];
	  else strm << '(' << el[k] << ',' << cl[k] << ')';
	} // for (k)
	strm << "\n"; // end line
      } // for (j)
      strm << "\n"; // empty line between cells
    } // for (i)
}

}

}
