/*
  This is wgraph.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "wgraph.h"
#include "filekl_in.h"
#include "bitset.h"
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
    } //for (z)
  } // for (i)

} // cells

/* The following function is an alternative to the function |wGraph| defined
   in kl.cpp. Here we do not assume that a KLContext is available, but that
   binary files with information about the block, matrix, and KL polynomials
   are avaialble. As a consequence we must redo the work of
   |kl::Helper::fillMuRow| as well as that to the mentioned |wGraph| function.
*/

WGraph wGraph
  ( std::ifstream& block_file
  , std::ifstream& matrix_file
  , std::ifstream& KL_file)
{
  using blocks::BlockElt;
  typedef std::auto_ptr<filekl::polynomial_info> pol_aptr;

  filekl::matrix_info mi(block_file,matrix_file);
  pol_aptr pol_p(NULL);

  try
  { pol_p=pol_aptr(new filekl::cached_pol_info(KL_file));
  }
  catch (std::exception& e)
  {
    std::cerr << "Failed to use cached polynomials: " << e.what() << std::endl;
    pol_p=pol_aptr(new filekl::polynomial_info(KL_file));
  }

  filekl::polynomial_info& poli=*pol_p;

  size_t max_mu=1;                       // maximal mu found
  std::pair<BlockElt,BlockElt> max_pair; // corresponding (x,y)

  WGraph result(mi.rank()); result.resize(mi.block_size());

  // fill in descent sets
  for (BlockElt y = 0; y < mi.block_size(); ++y)
  {
    bitset::RankFlags d_y = result.descent(y) = mi.descent_set(y);
    size_t ly = mi.length(y);
    if (ly==0) continue; // nothing more to do; avoid negative |d| below

#ifdef VERBOSE
    std::cerr << "reading edges for y=" << y;
    if (max_mu>1) std::cerr << ", maximal mu: " << max_mu;
    std::cerr << '\r';
#endif

    const filekl::strong_prim_list& spy=mi.strongly_primitives(y);

    filekl::strong_prim_list::const_iterator start= spy.begin();
    // traverse lengths |lx| of opposite parity to |ly|, up to |ly-3|
    for (size_t lx=(ly-1)%2,d = (ly-1)/2; d>0; --d,lx+=2) // d=(ly-1-lx)/2
    {
      filekl::strong_prim_list::const_iterator stop =
	std::lower_bound(start,spy.end()-1,mi.first_of_length(lx+1));
      for (start= std::lower_bound(start,stop,mi.first_of_length(lx));
	   start<stop; ++start)
	if (mi.descent_set(*start)!=d_y)
        {
	  BlockElt x = *start;
	  filekl::KLIndex klp = mi.find_pol_nr(x,y);

	  if (poli.degree(klp)==d)
	  {
	    result.edgeList(x).push_back(y);
	    size_t mu=poli.leading_coeff(klp);
	    if (mu>max_mu) { max_mu=mu; max_pair=std::make_pair(x,y); }
	    result.coeffList(x).push_back(mu);
	  }

	} // for (start...) if (descent!=d_y)
    } // for (lx,d)

    // for length |ly-1| we cannot limit ourselves to strongly primitives
    BlockElt end=mi.first_of_length(ly);
    for (BlockElt x=mi.first_of_length(ly-1); x<end; ++x)
    {
      bitset::RankFlags d_x=mi.descent_set(x);
      if (d_x==d_y) continue; // this case would lead nowhere anyway
      filekl::KLIndex klp = mi.find_pol_nr(x,y);
      if (klp!=filekl::KLIndex(0)) // then some edge between |x| and |y| exists
      {
	size_t mu=poli.leading_coeff(klp);
	if (mu>max_mu) { max_mu=mu; max_pair=std::make_pair(x,y); }
	if (not d_y.contains(d_x))
	{
	  result.edgeList(x).push_back(y);
	  result.coeffList(x).push_back(mu);
	}
	if (not d_x.contains(d_y))
	{
	  result.edgeList(y).push_back(x);
	  result.coeffList(y).push_back(mu);
	}
      } // if (klp!=KLIndex(0))
    } // for (x)
  } // for (y)

  size_t n_edges=0;
  for (BlockElt y = 0; y < mi.block_size(); ++y)
    n_edges+=result.edgeList(y).size();

  if (max_mu==1) std::cout << "All edges found are simple.\n";
  else std::cout << "Maximal edge multiplicity " << max_mu << " for ("
		 << max_pair.first << ',' << max_pair.second <<")\n";
  std::cout << "Total number of (directed) edges in graph is " << n_edges
	    << '.' << std::endl;

  return result;
}


}

}
