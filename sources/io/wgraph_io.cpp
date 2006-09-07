/*
  This is wgraph_io.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "wgraph_io.h"

#include <iostream>

#include "prettyprint.h"

namespace atlas {

/*****************************************************************************

        Chapter I -- Functions declared in wgraph_io.h

******************************************************************************/

namespace wgraph_io {

std::ostream& printCells(std::ostream& strm, const wgraph::WGraph& wg)

/*
  Synopsis: outputs the various cells in wg, as W-graphs, in the format
  of printWGraph.
*/

{
  using namespace wgraph;

  std::vector<WGraph> wc;
  cells(wc,wg);

  for (size_t j = 0; j < wc.size(); ++j) {
    strm << "// cell #" << j << std::endl;
    printWGraph(strm,wc[j]);
  }

  return strm;
}

std::ostream& printWGraph(std::ostream& strm, const wgraph::WGraph& wg)

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
    const CoeffList& cl = wg.coeffList(z);
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

  return strm;
}

}

}
