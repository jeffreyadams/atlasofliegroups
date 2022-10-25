/*
  This is dynkin.cpp
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2014,2017,2022 Marc van Leeuwen
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

//        |DynkinDiagram|, constructors

// Construct a presumptive Dynkin diagram from a Cartan matrix
DynkinDiagram::DynkinDiagram(const int_Matrix& c)
  : d_star(c.numColumns())
  , down_edges()
{
  assert(c.numRows()==c.numColumns());
  for (unsigned int i = 0; i < c.numRows(); ++i)
    for (unsigned int j = 0; j < c.numColumns(); ++j)
      if (i!=j)
      {
	if (c(i,j)<-3)
	  throw error::Cartan_error("Entry less than -3 present");
	else if (c(i,j)>0)
	  throw error::Cartan_error("Off-diagonal positive entry present");
	Multiplicity m = -c(i,j);
	if (m!=0)
	{
	  if (c(j,i)==0)
	    throw error::Cartan_error("Asymmetric adjacenty relation");
	  d_star[j].set(i); // set |i| as neighbour of |j|
	  if (m>1) // only label multiple edges
	    down_edges.push_back (std::make_pair(Edge(i,j),m));
	}
      }
      else if (c(i,j)!=2)
	throw error::Cartan_error("Diagonal entry unequal to 2");
  classify(c);
}

LieType DynkinDiagram::type() const
{
  LieType result;
  result.reserve(comps.size());
  for (const auto& comp: comps)
    result.emplace_back(comp.type,comp.rank());
  return result;
}

Permutation DynkinDiagram::perm() const
{
  Permutation result;
  result.reserve(rank());
  for (const auto& comp: comps)
    result.insert(result.end(),comp.position.begin(),comp.position.end());
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
	      result.down_edges.emplace_back(Edge(i,j),m);
	    else
	      result.down_edges.emplace_back(Edge(j,i),m);
	  }
	}
	else // unequal orbit lengths
	{ // then mark multiplicity 2 edge from longer to shorter orbit
	  if (diff<0)
	    result.down_edges.emplace_back(Edge(i,j),2);
	  else
	    result.down_edges.emplace_back(Edge(j,i),2);
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
  for (const auto& edge : down_edges)
    if (edge.first.first==i and edge.first.second==j)
      return -static_cast<int>(edge.second); // -2 or -3

  return -1; // simple edge, or labelled edge in short->long direction
}



// Return the (semisimple) Lie type of the Cartan matrix |cm|
LieType Lie_type(const int_Matrix& cm)
{
  return DynkinDiagram(cm).type();
}

DynkinDiagram& DynkinDiagram::classify(const int_Matrix& Cartan)
{
  for (unsigned i=0; i<rank(); ++i)
    if ((d_star[i].to_ulong() & constants::lMask[i])==0)
      comps.emplace_back(i);
    else
    {
      assert(not comps.empty()); // since support intersection is non-empty
      auto it=comps.begin();
      while ((it->support & d_star[i]).none()) // search first match
	assert(not comps.at_end(it)),++it;
      auto& sup = it->support;
      sup.set(i);
      ++it;
      while(not comps.at_end(it))
	if ((it->support & d_star[i]).any()) // search further matches
	{
	  sup |= it->support;
	  comps.erase(it);
	}
        else
	  ++it;
    }

  { // Now compute the offset and root positions for each of the components
    unsigned int i=0;
    for (auto& comp : comps)
    {
      comp.offset=i;
      unsigned size = comp.support.count();
      i += size;
      comp.position.reserve(size);
      classify(Cartan,comp);
    }
  }

  return *this;
}

/* Once a component is isolated in |ci| with |support| (and |offset|) defined,
   complete it by setting |type| and |position| */
void DynkinDiagram::classify(const int_Matrix& Cartan,comp_info& ci) const
{
  unsigned comp_rank = ci.support.count();
  if (comp_rank<=2)
  {
    if (comp_rank==1)
      ci.type='A',ci.position.push_back(ci.support.firstBit());
    else
    { // |comp_rank>2|
      auto it = ci.support.begin();
      unsigned i=*it, j=*++it;
      ci.position.push_back(i);
      ci.position.push_back(j);
      bool increase = Cartan(i,j)==-1;
      switch (Cartan(i,j)*Cartan(j,i))
      {
      case 1: ci.type='A'; break;
      case 2: // exceptionally in this case given order in diagram decides type
	ci.type= increase ? 'C' : 'B'; break;
      case 3: ci.type='G';
	if (not increase) std::swap(ci.position[0],ci.position[1]);
	break;
      default:
	throw error::Cartan_error("Off-diagonal pair both less than -1");
      }
    } // |if (comp_rank==1)|
  } // |if (comp_rank<=2)|
  else
  {
    unsigned fork=rank(), lower=rank(), upper=rank(); // out of bound values
    RankFlags extremities;
    for (unsigned i : ci.support)
      if (d_star[i].count()>2)
      {
	if (d_star[i].count()>3)
	  throw error::Cartan_error("Diagram node with degree more than 3");
	else if (fork<rank())
	  throw error::Cartan_error("Multiple fork nodes in component");
	fork=i;
      }
      else
      { if (d_star[i].count()<2)
	  extremities.set(i);
	for (unsigned j: d_star[i])
	  if (Cartan(i,j)<-1)
	  { if (lower<rank())
	      throw error::Cartan_error("Multiple labelled edges in component");
	    upper=i; lower=j; // this edge goes down from |i| to |j|
	  }
      }
    if (extremities.count() < 2)
      throw error::Cartan_error("Diagram has a loop");
    else if (fork<rank() and lower<rank())
      throw error::Cartan_error("Component with both fork and labelled edge");

    if (lower==rank()) // whether simply laced
      ci.type= fork==rank() ? 'A' :
	(d_star[fork]&extremities).count()==1 ? 'E' : 'D';
    else
      ci.type= extremities.test(lower) ? 'B' :
	extremities.test(upper) ? 'C' : 'F';

    unsigned start; RankFlags remain = ci.support;
    switch (ci.type)
    {
    default: assert(false); // only types ABCDEF just assigned
    case 'A': start=extremities.firstBit(); break; // choose any end
    case 'B': start=extremities.reset(lower).firstBit(); break;
    case 'C': start=extremities.reset(upper).firstBit(); break;
    case 'D':
      if (comp_rank==4)
	start=extremities.firstBit(); // choose any end
      else if (extremities.andnot(d_star[fork]).count()>1)
	throw error::Cartan_error("Fork node without adjacent extremities");
      else
      { assert(extremities.count()==1); // nothing left would be D4; excluded
	start=extremities.firstBit();
      }
      break;
    case 'E':
      { auto short_arm= extremities & d_star[fork];
	assert(short_arm.count()==1); // tested before setting |ci.type='E'|
	start = extremities.andnot(short_arm).firstBit(); // try a longer arm
	if ((d_star[start]&d_star[fork]).none()) // if this is longest arm
	  start = extremities.reset(start).firstBit(); // swap for orth arm
	if ((d_star[start]&d_star[fork]).none()) // if this still a long arm
	  throw error::Cartan_error("Fork node with two too long arms");
	ci.position.push_back(start); remain.reset(start);
	ci.position.push_back(short_arm.firstBit()); remain.andnot(short_arm);
	start=((d_star[start]&d_star[fork]).firstBit());
      }
      break;
    case 'F': start=(extremities.andnot(d_star[lower])).firstBit(); break;
    }

    // now traverse remainder of diagram starting from |start|
    while (ci.position.push_back(start),remain.reset(start).any())
      if (RankFlags cand=d_star[start] & remain)
	start = cand.firstBit();
      else if (ci.type=='D')
      { assert((remain-d_star[fork]).none()); // only a short arm remains
	start = remain.firstBit();
      }
      else
	assert(false); // must have fork, but type 'E' already used short arm
  }
}


/*
  Return the (semisimple) Lie type of the Cartan matrix cm, also sets |pi|
  to the permutation from the standard ordering of simple roots for that type.

  The returned permutation |pi| maps from Bourbaki ordering. If |check| holds,
  a complete test is made of all entries in |cm|, throwing a |runtime_error| if
  it fails to be a valid Cartan matrix (throwing an error can also happen when
  |check==false|, but in that case no effort is done to ensure it; therefore
  one should set |check| whenever the validty of |cm| is in doubt).
*/
LieType Lie_type(const int_Matrix& cm, bool check, Permutation& pi)
{
  if (check)
  {
    if (cm.numRows()!=cm.numColumns())
      throw error::Cartan_error();
    if (cm.numRows()>constants::RANK_MAX) // throw a different error type here
      throw std::runtime_error("Rank of matrix exceeds implementation limit");
  }

  DynkinDiagram d(cm);

  pi=d.perm();
  LieType result = d.type();

  if (check)
  {
   for (unsigned int i=0; i<result.size(); ++i)
    {
      SimpleLieType slt=result[i];
      if ((slt.type()=='E' and slt.rank()>8) or
	  (slt.type()=='F' and slt.rank()>4) or
	  (slt.type()=='G' and slt.rank()>2))
	throw error::Cartan_error("Excessive rank for exceptional type");
    }
    for (unsigned int i=0; i<d.rank(); ++i)
      for (unsigned int j=0; j<d.rank(); ++j)
	if (cm(pi[i],pi[j])!=result.Cartan_entry(i,j))
	  throw error::Cartan_error("");
   }
  return result;

}

} // |namespace dynkin|
} // |namespace atlas|
