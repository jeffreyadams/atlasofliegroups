/*
  This is wgraph.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2007,2017,2020 Marc van Leeuwen
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
  : d_rank(r), symmetric_graph(n), coefficients(n), descent_sets(n) {}

WGraph::WGraph(unsigned short rank,
	       const containers::sl_list<RankFlags>& tau,
	       const std::vector<containers::sl_list<kl::Mu_pair> >& edge_list)
  : d_rank(rank)
  , symmetric_graph(edge_list.size())
  , coefficients(edge_list.size())
  , descent_sets(tau.begin(),tau.end())
{
  assert(tau.size()==edge_list.size());
  graph::Vertex v = 0;
  for (const auto& list : edge_list)
  {
    auto& coef = coefficients[v];
    auto& edge = symmetric_graph.edgeList(v);
    coef.reserve(list.size()); edge.reserve(list.size());
    for (const auto& pair : list)
    {
      edge.push_back(pair.x);
      coef.push_back(pair.coef);
    }
    ++v;
  }
}

// remove edges $x\to y$ where |descent_set(x)| is contained in |descent_set(y)|
graph::OrientedGraph WGraph::oriented_graph() const
{
  graph::OrientedGraph result(size());
  for (graph::Vertex x=0; x<size(); ++x)
  {
    auto desc_x = descent_set(x);
    containers::sl_list<graph::Vertex> filtered;
    for (graph::Vertex y : edge_list(x))
      if (not descent_set(y).contains(desc_x))
	filtered.push_back(y);
    result.edgeList(x)=filtered.to_vector();
  }
  return result;
}

DecomposedWGraph::DecomposedWGraph(const WGraph& wg)
  : d_cell(), d_part(wg.size()), d_id(), d_induced()
{
  auto oriented = wg.oriented_graph();
  Partition pi = oriented.cells(&d_induced);

  d_cell.reserve(pi.classCount()); // there will be this many cells
  d_id.reserve(pi.classCount());    // and vectors of identification numbers

  std::vector<unsigned long int> relno(wg.size()); // number of element in its cell
  for (Partition::iterator it(pi); it(); ++it) // loop over cells
  {
    d_id.emplace_back(); // create new empty vector
    std::vector<unsigned long int>& idn=d_id.back(); // call it |idn|
    idn.reserve(it->second-it->first);  // then dimension it to the cell size

    containers::sl_list<RankFlags> tau;
    std::vector<containers::sl_list<kl::Mu_pair> > edge_lists;
    edge_lists.reserve(it->second-it->first); // reserve for cell size

    for (Partition::iterator::SubIterator jt=it->first; jt!=it->second; ++jt)
    { // ensure first that the |relno| fields within this cell are set
      idn.push_back(*jt); // vertex number in original graph
      relno[idn.back()]=jt-it->first; // relative number within this cell
    }

    for (Partition::iterator::SubIterator jt=it->first; jt!=it->second; ++jt)
    {
      size_t y = *jt; // vertex number in original graph
      d_part[y]=pi.class_of(y);
      tau.push_back(wg.descent_set(y)); // transfer descent set unchanged
      const auto& ely = wg.edge_list(y);
      edge_lists.emplace_back(); // append empty edge list
      auto& edge_list = edge_lists.back(); // that list, which is to be modified
      for (const graph::Vertex& x : ely)
	if (pi.class_of(x)==pi.class_of(y)) // ignore edges leading out of cell
	  edge_list.emplace_back(relno[x],wg.coefficient(y,&x-&ely[0]));
    }
    d_cell.emplace_back( wg.rank(), tau, edge_lists);
  } // |for (it)| loop over cells


} // |DecomposedWGraph::DecomposedWGraph|


/*****************************************************************************

        Chapter II -- Functions declared in wgraph.h

******************************************************************************/


// Return the cells of the W-graph |wg|
std::vector<WGraph> cells(const WGraph& wg)
{
  auto oriented = wg.oriented_graph();
  Partition pi = oriented.cells(); // no information about induced graph here

  std::vector<WGraph> result; result.reserve(pi.classCount());
  std::vector<unsigned long int> relno(wg.size()); // number of element in its cell

  for (Partition::iterator it(pi); it(); ++it) // loop over cells
  {
    containers::sl_list<RankFlags> tau;
    std::vector<containers::sl_list<kl::Mu_pair> > edge_lists;

    for (Partition::iterator::SubIterator jt=it->first; jt!=it->second; ++jt)
      relno[*jt]=jt-it->first; // relative number within this cell

    for (Partition::iterator::SubIterator jt=it->first; jt!=it->second; ++jt)
    {
      auto y = *jt;
      tau.push_back(wg.descent_set(y)); // transfer descent set unchanged
      const auto& ely = wg.edge_list(y);
      edge_lists.emplace_back(); // append empty edge list
      auto& edge_list = edge_lists.back(); // that list, which is to be modified
      for (const graph::Vertex& x : ely)
	if (pi.class_of(x)==pi.class_of(y)) // ignore edges leading out of cell
	  edge_list.emplace_back(relno[x],wg.coefficient(y,&x-&ely[0]));
    } // |for (z)|
    result.emplace_back( wg.rank(), tau, edge_lists);
  } // |for (it)|
  return result;
} // |cells|

/* The following function is an alternative to the function |wGraph| defined
   in kl.cpp. Here we do not assume that a KL_table is available, but that
   binary files with information about the block, matrix, and KL polynomials
   are available. As a consequence we must redo the work of
   |kl::Helper::fillMuRow| as well as that to the mentioned |wGraph| function.
*/
WGraph wGraph
  ( std::ifstream& block_file
  , std::ifstream& matrix_file
  , std::ifstream& KL_file)
{
  using pol_uptr = std::unique_ptr<filekl::polynomial_info>;

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

  containers::sl_list<RankFlags> tau;
  std::vector<containers::sl_list<kl::Mu_pair> > edge_lists;

  // fill in descent sets
  for (BlockElt y = 0; y < mi.block_size(); ++y)
  {
    auto desc_y = mi.descent_set(y);
    tau.push_back(desc_y);
    edge_lists.emplace_back();
    auto ly = mi.length(y);
    if (ly==0)
      continue; // nothing more to do this time around; avoid negative |d| below

#ifdef VERBOSE
    std::cerr << "reading edges for y=" << y;
    if (max_mu>1) std::cerr << ", maximal mu: " << max_mu;
    std::cerr << '\r';
#endif
    auto& edge_list_y = edge_lists.back();
    const filekl::strong_prim_list& spy=mi.strongly_primitives(y);

    filekl::strong_prim_list::const_iterator start= spy.begin();
    // traverse lengths |lx| of opposite parity to |ly|, up to |ly-3|
    for (unsigned lx=(ly-1)%2,d = (ly-1)/2; d>0; --d,lx+=2) // d=(ly-1-lx)/2
    {
      filekl::strong_prim_list::const_iterator stop =
	std::lower_bound(start,spy.end()-1,mi.first_of_length(lx+1));
      for (start= std::lower_bound(start,stop,mi.first_of_length(lx));
	   start<stop; ++start)
	if (mi.descent_set(*start)!=desc_y)
        {
	  BlockElt x = *start;
	  filekl::KLIndex klp = mi.find_pol_nr(x,y);

	  if (poli.degree(klp)==d)
	  {
	    auto mu=poli.leading_coeff(klp);
	    edge_lists[x].emplace_back(y,mu);
	    if (mu>max_mu) { max_mu=mu; max_pair=std::make_pair(x,y); }
	  }

	} // for (start...) if (descent!=d_y)
    } // for (lx,d)

    // for length |ly-1| we cannot limit ourselves to strongly primitives
    BlockElt end=mi.first_of_length(ly);
    for (BlockElt x=mi.first_of_length(ly-1); x<end; ++x)
    {
      RankFlags desc_x=mi.descent_set(x);
      if (desc_x==desc_y) continue; // this case would lead nowhere anyway
      filekl::KLIndex klp = mi.find_pol_nr(x,y);
      if (klp!=filekl::KLIndex(0)) // then some edge between |x| and |y| exists
      {
	size_t mu=poli.leading_coeff(klp);
	if (mu>max_mu) { max_mu=mu; max_pair=std::make_pair(x,y); }
	if (not desc_y.contains(desc_x))
	  edge_lists[x].emplace_back(y,mu);
	if (not desc_x.contains(desc_y))
	  edge_list_y.emplace_back(x,mu);
      } // if (klp!=KLIndex(0))
    } // for (x)
  } // for (y)

#ifdef VERBOSE
  size_t n_edges=0;
  for (const auto& el : edge_lists)
    n_edges += el.size();

  if (max_mu==1) std::cout << "All edges found are simple.\n";
  else std::cout << "Maximal edge multiplicity " << max_mu << " for ("
		 << max_pair.first << ',' << max_pair.second <<")\n";
  std::cout << "Total number of (directed) edges in graph is " << n_edges
	    << '.' << std::endl;
#endif

  return  {  static_cast<unsigned short>(mi.rank()), tau, edge_lists };
} // |wGraph| (from file data)


} // |namespace wgraph|

} // |namespace atlas|
