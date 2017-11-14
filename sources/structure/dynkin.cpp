/*
  This is dynkin.cpp
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2014,2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Implementation of the class DynkinDiagram.

#include "dynkin.h"

#include <cassert>
#include <algorithm>
#include <string> // used implicitly in throwing |std::runtime_error|

#include "matrix.h"
#include "lietype.h"
#include "error.h"

namespace atlas {
namespace dynkin{

/*****************************************************************************

  Class and functions related to analysis of Dynkin diagrams

******************************************************************************/

namespace {

  Permutation order_by_components
    (const containers::sl_list<RankFlags>& components, unsigned int rank);

  Permutation typeANormalize(const DynkinDiagram& d);
  Permutation typeBNormalize(const DynkinDiagram& d, bool Bourbaki);
  Permutation typeCNormalize(const DynkinDiagram& d, bool Bourbaki);
  Permutation typeDNormalize(const DynkinDiagram& d, bool Bourbaki);
  Permutation typeENormalize(const DynkinDiagram& d);
  Permutation typeFNormalize(const DynkinDiagram& d);
  Permutation typeGNormalize(const DynkinDiagram& d);

} // |namespace|



//        |DynkinDiagram|, constructors

// Construct a presumptive Dynkin diagram from a Cartan matrix
DynkinDiagram::DynkinDiagram(const int_Matrix& c)
  : d_star(c.numColumns())
  , d_downedge()
{
  assert(c.numRows()==c.numColumns());
  for (unsigned int i = 0; i < c.numRows(); ++i)
    for (unsigned int j = 0; j < c.numColumns(); ++j)
      if (i!=j)
      {
	if (c(i,j)<-3 or c(i,j)>0)
	  throw error::CartanError();
	Multiplicity m = -c(i,j);
	if (m!=0)
	{
	  if (c(j,i)==0) // this test ensures a symmetric adjacency matrix
	    throw error::CartanError();
	  d_star[j].set(i); // set |i| as neighbour of |j|
	  if (m>1) // only label multiple edges
	    d_downedge.push_back (std::make_pair(Edge(i,j),m));
	}
      }
      else if (c(i,j)!=2)
	throw error::CartanError();
}


DynkinDiagram DynkinDiagram::subdiagram(RankFlags selection) const
{
  DynkinDiagram result; auto& r_star = result.d_star;
  r_star.reserve(selection.count());
  for (auto it = selection.begin(); it(); ++it)
  {
    r_star.push_back(star(*it));
    r_star.back().slice(selection); // extract bits set in |selection|, repack
  }
  for (auto it=d_downedge.begin(); it!=d_downedge.end(); ++it)
  {
    auto l = it->first.first, s = it->first.second;
    if (selection[l] and selection[s]) // is the edge retained in |selection|
    {
      Edge e(selection.position(l),selection.position(s)); // re-index
      result.d_downedge.emplace_back(e,it->second); // attach multiplicity
    }

  }
  return result;
}


DynkinDiagram DynkinDiagram::folded(const ext_gens& orbits) const
{
  DynkinDiagram result; auto& r_star = result.d_star;
  r_star.resize(orbits.size());

  for (unsigned int i=0; i<orbits.size(); ++i)
  {
    RankFlags neighbours = star(orbits[i].s0); // neighbours of first elt
    if (orbits[i].length()>1)
      neighbours |= star(orbits[i].s1); // add neighbours of second elt
    for (unsigned j=orbits.size(); --j > i; )
      if (neighbours[orbits[j].s0])  // run over neigbours |j| with $i<j<n$
      {
	const unsigned jj = orbits[j].s0;
	unsigned ii = orbits[i].s0;
	if (not are_adjacent(ii,jj)) // then we got wrong orbit elements
	  ii = orbits[i].s1; // fix first one
	assert (are_adjacent(ii,jj));
	r_star[i].set(j);
	r_star[j].set(i); // mark |i| and |j| as neighbours in new diagram

	int diff=orbits[i].length()-orbits[j].length();
	if (diff==0) // equal length
	{ // then copy Cartan entries
	  const Multiplicity m = edge_multiplicity(ii,jj);
	  if (m>1)
	  {
	    if (Cartan_entry(ii,jj)<-1)
	      result.d_downedge.emplace_back(Edge(i,j),m);
	    else
	      result.d_downedge.emplace_back(Edge(j,i),m);
	  }
	}
	else // unequal orbit lengths
	{ // then mark multiplicity 2 edge from longer to shorter orbit
	  if (diff<0)
	    result.d_downedge.emplace_back(Edge(i,j),2);
	  else
	    result.d_downedge.emplace_back(Edge(j,i),2);
	}
      }
  }
  return result;
}


/******** accessors **********************************************************/


// recover Cartan matrix entry from Dynkin diagram
int DynkinDiagram::Cartan_entry(unsigned int i,unsigned int j) const
{
  if (not are_adjacent(i,j))
    return i==j ? 2 : 0;
  for (unsigned int k=0; k<d_downedge.size(); ++k)
    if (d_downedge[k].first.first==i and d_downedge[k].first.second==j)
      return -static_cast<int>(d_downedge[k].second); // -2 or -3

  return -1; // simple edge, or labelled edge in short->long direction
}

// Find the set of terminal nodes (degree one or zero) of the graph.
RankFlags DynkinDiagram::extremities() const
{
  RankFlags e;

  for (unsigned int i = 0; i < d_star.size(); ++i)
    e.set(i,d_star[i].count() <= 1);

  return e;
}


/*
  Find the labelled (multiple) edge in connected diagram; assumed present

  The edge |e| is downwards: |e.first| is longer than |e.second|
*/
Edge DynkinDiagram::labelled_edge() const
{
  assert(d_downedge.size()>0);
  return d_downedge[0].first;
}


/*
  Find the largest multiplicity in the graph.

  Returns 1 in absence of labelled edges, even when there are no edges at all!
*/
Multiplicity DynkinDiagram::edge_label() const
{
  Multiplicity m = 1;

  for (unsigned int i =0; i<d_downedge.size(); ++i)
    if (d_downedge[i].second>m)
      m = d_downedge[i].second;

  return m;
}


/*
  Find a fork node (|degree >= 3|) in a connected non-linear graph.

  Returns exception value |rank()| if the graph does not have a fork node.
*/
unsigned int DynkinDiagram::fork_node() const
{
  for (unsigned int i = 0; i < rank(); ++i)
    if (d_star[i].count() >= 3)
      return i;

  return rank();
}


/*
  Return the connected component containing vertex number|i| in the diagram.

  We use the class invariant that the adjacency matrix is symmetric, so that
  connected components are equivalence classes for reachability

  The algorithm is to start with |i|, and to construct "shells" from there, by
  taking each new shell to be the elements of the union of the stars of the
  old shell, that were not already considered. Since bitset iterators copy
  their bitset at construction, adding to |newElts| in the inner loop will not
  affect that iteration itself (but the logic would not be broken if it did).
 */
RankFlags DynkinDiagram::component(unsigned int i) const
{
  RankFlags result, new_elts;
  new_elts.set(i);

  while (new_elts.any())
  {
    result |= new_elts; // transfer
    for (RankFlags::iterator it = new_elts.begin(); it(); ++it)
      new_elts |= d_star[*it];
    new_elts.andnot(result); // remove any bits that were already previouslt set
  }

  return result;
}

// Return our decomposition into connected components, a list of subsets
containers::sl_list<RankFlags> DynkinDiagram::components() const
{
  containers::sl_list<RankFlags> result;
  RankFlags remainder; remainder.fill(rank());
  while(remainder.any())
  {
    RankFlags c = component(remainder.firstBit());
    remainder.andnot(c); // remove current component from remainder
    result.push_back(c);
  }

  return result;
}

// Return the (semisimple) Lie type of the Cartan matrix |cm|
LieType Lie_type(const int_Matrix& cm)
{
  DynkinDiagram diagram(cm);
  const auto cl = diagram.components();

  LieType result;
  result.reserve(cl.size());
  for (auto it=cl.begin(); it!=cl.end(); ++it)
  {
    auto cd = diagram.subdiagram(*it);
    result.emplace_back(cd.component_kind(),cd.rank());
  }

  return result;
}


/*
  Precondition: any constructed Dynkin diagram is acceptable

  Postcondition: |a| is modified so that the new permutation gives a
  normalized ordering on each component: it gives Bourbaki ordering unless
  |Bourbaki| is false and the type is BCD: then the order is reversed.

  The detected semisimple Lie type is returned.
*/
LieType DynkinDiagram::classify_semisimple(Permutation& a,bool Bourbaki) const
{
  const auto cl = components();
  a = order_by_components(cl,rank());

  unsigned int offset = 0;
  LieType result; result.reserve(cl.size());
  for (auto it=cl.begin(); it!=cl.end(); ++it)
  {
    // make a Dynkin diagram for the connected component
    auto cd = subdiagram(*it);

    // normalize it
    Permutation b;
    result.push_back(cd.classify_simple(b,Bourbaki));

    // piece together the permutation
    permutations::compose(a,b,offset);

    // update offset
    offset += it->count();
  }
  return result;
} // |classify_semisimple|


/*
  Precondition : the current object is a connected Dynkin diagram;

  Postcondition : |pi| holds a permutation which enumerates the vertices of
  |d| in an order that will induce a normal form of |*this|; a CartanError
  will be thrown (either here or in helper) in case it is not a valid diagram

  It is just a dispatching function for the various possible simple types.
*/
SimpleLieType
  DynkinDiagram::classify_simple(Permutation& pi, bool Bourbaki) const
{
  lietype::TypeLetter x = component_kind();

  switch (x)
  {
  case 'A': pi = typeANormalize(*this);
    break;
  case 'B': pi = typeBNormalize(*this,Bourbaki);
    break;
  case 'C': pi = typeCNormalize(*this,Bourbaki);
    break;
  case 'D': pi = typeDNormalize(*this,Bourbaki);
    break;
  case 'E': pi = typeENormalize(*this);
    break;
  case 'F': pi = typeFNormalize(*this);
    break;
  case 'G': pi = typeGNormalize(*this);
    break;
  default:
    pi = Permutation(); // will provoke the error below
  }
  if (pi.size()!=rank())
    throw error::CartanError();

  return SimpleLieType(x,rank());
} // |classify_simple|


/*
  Determines candidate for the (simple) type of a connected Dynkin diagram
  Throws an error if no candidate is found, but if it returns this does not
  guarantee that the Dynkin diagram is correct.

  Precondition : |d| is connected (and therefore not empty)
*/
lietype::TypeLetter DynkinDiagram::component_kind() const
{
  if (rank()<=2) // types A1,A2,B2,C2,G2 are validly possible
    switch (edge_label())
    {
    case 1: return 'A';
    case 3: return 'G';
    case 2: // exceptionally in this case given order in diagram decides type
      return labelled_edge().first==0 ? 'B' : 'C'; // Bourbaki, B starts long
    default: assert(false); // too high label, should have been caught before
   }

  else // |rank>2|
  {
    RankFlags extr = extremities();
    if (extr.count()<2)
      throw error::CartanError();
    unsigned int fork = fork_node();
    if (fork==rank()) // diagram is linear
      switch (edge_label())
      {
      case 1: return 'A';
      case 2: // type is B,C or F
	{
	  Edge e = labelled_edge();
	  return extr.test(e.first) ? 'C' : extr.test(e.second) ? 'B' : 'F';
	}
      default: {} // since G2 was already detected, this canot be right
      }
    else // not a linear diagram
    {
      RankFlags short_arms = star(fork) & extremities();
      if (short_arms.any()) // now |arms| counts short arms
	return short_arms.count() == 1 ? 'E' : 'D';
    }
  }
  throw error::CartanError();
} // |DynkinDiagram::component_kind|


/*
  Return a permutation such that successive intervals of simple roots form
  connected components, numbered as needed for our Weyl group implementation.

  NOTE: the permutation |result| maps new index |i| to old index |result[i]|.
*/
Permutation normalize(const DynkinDiagram& d)
{
  Permutation result;
  d.classify_semisimple(result,false);
  return result;
}



/*
  Return the (semisimple) Lie type of the Cartan matrix cm, also sets |pi|
  to the permutation from the standard ordering of simple roots for that type.

  Standard ordering is taken as Bourbaki ordering if |Bourbaki| holds,
  internal Weyl group implementation ordering otherwise. If |check| holds, a
  complete test is made of all entries in |cm|, throwing a |runtime_error| if
  it fails to be a valid Cartan matrix (throwing an error can also happen when
  |check==false|, but in that case no effort is done to ensure it; therefore
  one should set |check| whenever the validty of |cm| is in doubt).
*/
LieType Lie_type(const int_Matrix& cm,
		 bool Bourbaki, bool check,
		 Permutation& pi)
{
  if (check)
  {
    if (cm.numRows()!=cm.numColumns())
      throw error::CartanError();
    if (cm.numRows()>constants::RANK_MAX) // throw a different error type here
      throw std::runtime_error("Rank of matrix exceeds implementation limit");
  }

  DynkinDiagram d(cm);

  LieType result = d.classify_semisimple(pi,Bourbaki);

  if (check)
  {
   for (unsigned int i=0; i<result.size(); ++i)
    {
      SimpleLieType slt=result[i];
      if ((slt.type()=='E' and slt.rank()>8) or
	  (slt.type()=='F' and slt.rank()>4) or
	  (slt.type()=='G' and slt.rank()>2))
	throw error::CartanError();
    }
    for (unsigned int i=0; i<d.rank(); ++i)
      for (unsigned int j=0; j<d.rank(); ++j)
	if (cm(pi[i],pi[j])!=result.Cartan_entry(i,j))
	  throw error::CartanError();
   }
  return result;

}
/*!
  Synopsis: Returns some permutation that will take |d| to Bourbaki form

  This means that nodes of the diagram |d| taken in the order |a[0],...,a[r-1]|
  traverse each of its connected components consecutively, and in the order
  prescribed by the the Bourbaki conventions for the type of that component
*/
Permutation bourbaki(const DynkinDiagram& d)
{
  // do the normalization as in normalize, but with Bourbaki ordering
  Permutation result;
  d.classify_semisimple(result,true);

  return result;
}



/*****************************************************************************

        Chapter III -- Auxiliary functions for this module

******************************************************************************/

namespace {

/*
  Returns a permutation such that the various components, listed in cl,
  are numbered by successive indices. The result maps these indices back
  to their original positions.
*/
Permutation order_by_components
  (const containers::sl_list<RankFlags>& cl, unsigned int r)
{
  Permutation result; result.reserve(r);

  // traverse each component, write down its elements in sequence
  for (auto it=cl.begin(); it!=cl.end(); ++it)
    std::copy(it->begin(),it->end(),std::back_inserter(result));

  assert (result.size()==r); // check that correct rank was passed
  return result;
} // |order_by_components|

} // |namespace|


namespace {

// an auxiliary function; a first element is already pushed onto |a|
RankFlags linearise(const DynkinDiagram& d, Permutation& a)
{
  RankFlags done,next=d.star(a.back());
  while(done.set(a.back()),(next=d.star(a.back()).andnot(done)).any())
    a.push_back(next.firstBit());

  return done;
}

/*
  Find a permutation that will enumerates |d| along its diagram

  Precondition : |d.edge_label()==1| and |d.fork_node()==d.rank()|. This
  actually ensures this is a valid type An diagram.

  Postcondition : |pi| linearly enumerates the diagram in one of the two
  (except for A1) possible orders
*/
Permutation typeANormalize(const DynkinDiagram& d)
{
  Permutation a;
  a.reserve(d.rank());
  a.push_back(d.extremities().firstBit());

  linearise(d,a);
  return a;
}


/*
  Puts in |a| a permutation that will enumerate |d| in linear order,
  ending with a labelled edge if |Bourbaki| holds, or starting if not.

  Precondition |d.fork_node()==d.rank()| (linear diagram) and
  |d.edge_label()==2|, also this was not classified as type C or F

  There is a unique such ordering
*/
Permutation typeBNormalize(const DynkinDiagram& d, bool Bourbaki)
{
  Permutation a;
  a.reserve(d.rank());
  unsigned int short_node = d.labelled_edge().second; // one of 2 end nodes
  a.push_back // Bourbaki starts with the extramal node that is not short
    (Bourbaki ? d.extremities().reset(short_node).firstBit() : short_node);
  linearise(d,a);
  return a;
}


/*
  Puts in |a| a permutation that will enumerate |d| in linear order,
  ending with a labelled edge if |Bourbaki| holds, or starting if not.

  Precondition : as for type B, but one extremal node was long end of edge

  There is a unique such ordering
*/
Permutation typeCNormalize(const DynkinDiagram& d, bool Bourbaki)
{
  Permutation a;
  a.reserve(d.rank());
  unsigned int long_node = d.labelled_edge().first; // also an end point
  a.push_back // Bourbaki starts with the extramal node that is not long
    (Bourbaki ? d.extremities().reset(long_node).firstBit() : long_node);
  linearise(d,a);
  return a;
}


/*
  Puts in |a| a permutation that will enumerate |d| in alomst linear order
  (only the fork node has one neighbour at index distance 2 from it, which
  index is extremal); the fork node is at index |Bourbaki ? rank-3 : 2|.

  Precondition : there is a fork node (|degree >= 3|) with at least two short
  arms (neighbours that are extremities). Necessarily |d.rank()>=4|.

*/
Permutation typeDNormalize(const DynkinDiagram& d, bool Bourbaki)
{
  unsigned int r = d.rank();
  Permutation a;
  a.reserve(r);

  unsigned int fork = d.fork_node();
  RankFlags short_arms = d.star(fork) & d.extremities();
  RankFlags long_ends = d.extremities().andnot(short_arms);
  RankFlags::iterator it = short_arms.begin();

  if (long_ends.none()) // we either have a D4 diagram or rubbish
  {
    assert(short_arms.count()>=3); // the code below will not run out of |it|
    a.push_back(*it);
    if (Bourbaki)
      a.push_back(fork),a.push_back(*++it);
    else
      a.push_back(*++it),a.push_back(fork);
    a.push_back(*++it);
    // if |r>4| (diagram is rubbish) then |a| incomplete, will throw an error
  }
  else // an extremity not adjacent to the fork node was found
  {
    a.push_back(long_ends.firstBit());
    RankFlags done=linearise(d,a); // will pass through |fork| in valid cases
    if (done[*it])
      ++it; // skip first short arm if already done
    a.push_back(*it); // add an unused short arm vertex (at most one is used)
    if (not Bourbaki)
      std::reverse(a.begin(),a.end());
  }
  return a; // might be incomplete, and completeness does not ensure validity
}


/*
  Put in |a| the permutation enumerating |d| in Bourbaki order.

  Precondition : |d| has a fork node with one short arm

  Postcondition : a holds a permutation for which the node is in position 3
  (counting from 0), position 1 is the extremity of the branch of length 1,
  position 0 is the extremity of a branch of length 2, position 2 is the other
  element of that branch, and the elements of the last branch are enumerated
  from the node. There are two solutions in type E6, one otherwise.
*/
Permutation typeENormalize(const DynkinDiagram& d)
{
  unsigned int r = d.rank();
  unsigned int fork = d.fork_node();
  RankFlags fork_star = d.star(fork);
  RankFlags extr = d.extremities();

  Permutation a;
  if (r<6 or r>8 or fork_star.count()!=3 or extr.count()!=3)
    throw error::CartanError();

  RankFlags short_arms = fork_star & extr;
  assert(short_arms.count()==1); // this was what caused type E classification

  extr.andnot(short_arms);
  RankFlags::iterator it = extr.begin();
  RankFlags inter = d.star(*it) & fork_star;
  if (inter.none()) // skip end point long arm
    inter = d.star(*++it) & fork_star;
  if (inter.none()) // still long arm? then something is wrong
    throw error::CartanError();

  a.push_back(*it);                   // end middle arm
  a.push_back(short_arms.firstBit()); // short arm
  a.push_back(inter.firstBit());      // halfway middle are
  a.push_back(fork);

  inter |= short_arms;     // henceforth |inter| marks done nodes
  fork_star.andnot(inter); // and |fork_star| successors of last node

  while (fork_star.any())
  {
    inter.set(a.back());
    a.push_back(fork_star.firstBit());
    fork_star=d.star(a.back()).andnot(inter);
  }
  return a;
}


/*
  Put in |a| the permutation enumerating |d| in Bourbaki order.

  Precondition : diagram is linear, has multiple edge without extremities

  Postcondition : |a| holds a permutation which enumerates the graph in linear
  order, for which the middle edge is oriented from 1 to 2; (the arrow in F4
  is like in the Bn diagrams); such a permutation is unique.
*/
Permutation typeFNormalize(const DynkinDiagram& d)
{
  Permutation a;
  Edge e = d.labelled_edge(); // there is such an edge
  RankFlags st = d.star(e.first);
  assert(st.count()>1); // this was tested
  st.reset(e.second);

  a.push_back(st.firstBit());
  linearise(d,a);
  return a;
}


// Precondition : |d| has rank 2 and an edge with label 3
Permutation typeGNormalize(const DynkinDiagram& d)
{
  Permutation a(2);

  Edge e = d.labelled_edge();

  a[1] = e.first;
  a[0] = e.second;
  return a;
}

} // |namespace|
} // |namespace dynkin|
} // |namespace atlas|
