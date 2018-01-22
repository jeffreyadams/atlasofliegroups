/*
  This is wgraph.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2007,2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "wgraph.h"

#include <memory>
#include <algorithm>
#include <iostream>

#include "filekl_in.h"	// for alternative |wGraph| function

namespace atlas {

/*****************************************************************************

        Chapter I -- The WGraph and DecomposedWGraph classes

******************************************************************************/

namespace wgraph {

/******** constructors and destructors ***************************************/

WGraph::WGraph(size_t r, size_t n)
  : d_rank(r), oriented_graph(n), coefficients(n), descent_sets(n) {}


/******** copy, assignment and swap ******************************************/

/******** manipulators *******************************************************/

DecomposedWGraph::DecomposedWGraph(const WGraph& wg)
  : d_cell(), d_part(wg.size()), d_id(), d_induced()
{
  Partition pi = wg.oriented_graph.cells(&d_induced);

  d_cell.reserve(pi.classCount()); // there will be this many cells
  d_id.resize(pi.classCount());    // and vectors of identification numbers


  std::vector<unsigned int> relno(wg.size()); // number of element in its cell
  Partition::iterator it(pi);
  for (cell_no n=0; n<d_id.size(); ++n,++it)
  {
    d_cell.push_back(WGraph(wg.rank(),it->second-it->first));
    WGraph& cur_cell = d_cell.back();
    std::vector<unsigned int>& idn=d_id[n];
    idn.resize(cur_cell.size());

    for (Partition::iterator::SubIterator j=it->first; j!=it->second; ++j)
    {
      size_t y = *j; size_t z=j-it->first; // |y| gets renamed |z| in cell

      d_part[y]=n; // or equivalently |d_part[y]=pi(y)|
      relno[y]=z; idn[z]=y;

      // transfer descent set unchanged
      cur_cell.descent_sets[z] = wg.descent_sets[y];
    }
  }
  // we have made sure all values |relno[y]| are defined before proceeding

  for (it.rewind(); it(); ++it)
  {
    const cell_no n=d_part[*it->first]; // cell number, fixed for next loop
    auto& cell = d_cell[n];
    auto& graph = cell.oriented_graph;
    for (Partition::iterator::SubIterator j=it->first; j!=it->second; ++j)
    {
      size_t y = *j; size_t z=relno[y]; // |z==j-it->first|
      for (size_t k = 0; k < wg.degree(y); ++k)
	if (d_part[wg.edge_target(y,k)]==n) // ignore edges outside this cell
	{
	  graph.edgeList(z).push_back(relno[wg.edge_target(y,k)]);
	  cell.coefficients[z].push_back(wg.coefficient(y,k));
	}
    } // for (j)
  } // for (it)

} // |DecomposedWGraph::DecomposedWGraph|


/*****************************************************************************

        Chapter II -- Functions declared in wgraph.h

******************************************************************************/


// Return the cells of the W-graph |wg|
std::vector<WGraph> cells(const WGraph& wg)
{
  Partition pi = wg.oriented_graph.cells(nullptr);
  // we do not collect information about induced graph here

  std::vector<WGraph> result; result.reserve(pi.classCount());

  for (Partition::iterator it(pi); it(); ++it) // loop over cells
  {
    result.push_back(WGraph(wg.rank(),it->second - it->first));
    WGraph& cur_cell = result.back();
    auto& graph = cur_cell.oriented_graph;

    /* looping over |z| rather than using |*it| directly implements the
       renumbering of each cell (what was |y=*(it->first+z)| becomes just |z|)
    */
    for (size_t z = 0; z < cur_cell.size(); ++z)
    {
      size_t y = it->first[z];
      cur_cell.descent_sets[z] = wg.descent_sets[y];
      for (size_t j = 0; j < wg.degree(y); ++j)
      {
	graph::Vertex x = wg.edge_target(y,j);
	if (pi.class_of(x) == pi.class_of(y))
	{ // ignore edge leading out of the current cell
	  graph::Vertex xi = // find relative address of x in this class
	    std::lower_bound(it->first,it->second,x) - it->first;
	  graph.edgeList(z).push_back(xi);
	  cur_cell.coefficients[z].push_back(wg.coefficient(y,j));
	}
      }
    } // |for (z)|
  } // |for (it)|
  return result;
} // |cells|

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
  typedef std::unique_ptr<filekl::polynomial_info> pol_uptr;

  filekl::matrix_info mi(block_file,matrix_file);
  pol_uptr pol_p(nullptr);

  try
  { pol_p=pol_uptr(new filekl::cached_pol_info(KL_file));
  }
  catch (std::exception& e)
  {
    std::cerr << "Failed to use cached polynomials: " << e.what() << std::endl;
    pol_p=pol_uptr(new filekl::polynomial_info(KL_file));
  }

  filekl::polynomial_info& poli=*pol_p;

  size_t max_mu=1;                       // maximal mu found
  std::pair<BlockElt,BlockElt> max_pair; // corresponding (x,y)

  WGraph result(mi.rank(),mi.block_size());

  // fill in descent sets
  for (BlockElt y = 0; y < mi.block_size(); ++y)
  {
    RankFlags d_y = result.descent_sets[y] = mi.descent_set(y);
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
	    result.oriented_graph.edgeList(x).push_back(y);
	    size_t mu=poli.leading_coeff(klp);
	    if (mu>max_mu) { max_mu=mu; max_pair=std::make_pair(x,y); }
	    result.coefficients[x].push_back(mu);
	  }

	} // for (start...) if (descent!=d_y)
    } // for (lx,d)

    // for length |ly-1| we cannot limit ourselves to strongly primitives
    BlockElt end=mi.first_of_length(ly);
    for (BlockElt x=mi.first_of_length(ly-1); x<end; ++x)
    {
      RankFlags d_x=mi.descent_set(x);
      if (d_x==d_y) continue; // this case would lead nowhere anyway
      filekl::KLIndex klp = mi.find_pol_nr(x,y);
      if (klp!=filekl::KLIndex(0)) // then some edge between |x| and |y| exists
      {
	size_t mu=poli.leading_coeff(klp);
	if (mu>max_mu) { max_mu=mu; max_pair=std::make_pair(x,y); }
	if (not d_y.contains(d_x))
	{
	  result.oriented_graph.edgeList(x).push_back(y);
	  result.coefficients[x].push_back(mu);
	}
	if (not d_x.contains(d_y))
	{
	  result.oriented_graph.edgeList(y).push_back(x);
	  result.coefficients[y].push_back(mu);
	}
      } // if (klp!=KLIndex(0))
    } // for (x)
  } // for (y)

  size_t n_edges=0;
  for (BlockElt y = 0; y < mi.block_size(); ++y)
    n_edges += result.degree(y);

  if (max_mu==1) std::cout << "All edges found are simple.\n";
  else std::cout << "Maximal edge multiplicity " << max_mu << " for ("
		 << max_pair.first << ',' << max_pair.second <<")\n";
  std::cout << "Total number of (directed) edges in graph is " << n_edges
	    << '.' << std::endl;

  return result;
} // |wGraph| (from file data)


} // |namespace wgraph|

} // |namespace atlas|
