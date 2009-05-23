/*!
\file
\brief Implementation of the class DynkinDiagram.
*/
/*
  This is dynkin.cpp
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "dynkin.h"

#include <cassert>

#include "latticetypes.h"
#include "lietype.h"

/*****************************************************************************

  Class and functions related to analysis of Dynkin diagrams

******************************************************************************/

namespace atlas {
namespace dynkin{
namespace {

  setutils::Permutation componentOrder
    (const bitset::RankFlagsList& components, size_t rank);
  lietype::LieType componentNormalize
    (setutils::Permutation&, const bitset::RankFlagsList&,
     const DynkinDiagram&, bool Bourbaki);
  lietype::SimpleLieType
    irreducibleNormalize(setutils::Permutation&,
			 const DynkinDiagram& d, bool Bourbaki);
  lietype::TypeLetter irreducibleType(const DynkinDiagram&);
  void typeANormalize(setutils::Permutation&, const DynkinDiagram&);
  void typeBNormalize(setutils::Permutation&,
		      const DynkinDiagram& d, bool Bourbaki);
  void typeCNormalize(setutils::Permutation&,
		      const DynkinDiagram& d, bool Bourbaki);
  void typeDNormalize(setutils::Permutation&,
		      const DynkinDiagram& d, bool Bourbaki);
  void typeENormalize(setutils::Permutation&, const DynkinDiagram&);
  void typeFNormalize(setutils::Permutation&, const DynkinDiagram&);
  void typeGNormalize(setutils::Permutation&, const DynkinDiagram&);

} // |namespace|
} // |namespace dynkin|

/*****************************************************************************

        Chapter I -- The DynkinDiagram class

  Represents the structure given by a Cartan matrix in graph form

******************************************************************************/

namespace dynkin {

/******** constructors and destructors ***************************************/

/*!
  Constructs a Dynkin diagram from a Cartan matrix. The edge points from i to
  j when the Cartan matrix entry c(i,j) is < -1 (i.e. towards the shorter
  root.)
*/
DynkinDiagram::DynkinDiagram(const latticetypes::LatticeMatrix& c)
  : d_star(c.numColumns())
  , d_downedge()
{
  assert(c.numRows()==c.numColumns());
  for (size_t j = 0; j < c.numColumns(); ++j)
    for (size_t i = 0; i < c.numRows(); ++i)
      if (c(i,j)!=0 and (i != j))
      {
	d_star[j].set(i);
	if (c(i,j) < -1) // only label multiple edges
	  d_downedge.push_back(std::make_pair(Edge(i,j),-c(i,j)));
      }
}

/*!
  Constructs the restriction of |d| to the subset of the vertices flagged
  by |c|.
*/
DynkinDiagram::DynkinDiagram(const bitset::RankFlags& c,
			     const DynkinDiagram& parent)
  : d_star()     // start with empty vector
  , d_downedge() // start without any downward edges
{
  d_star.reserve(c.count()); // vertices selected by |c|, renumbered from 0

  // get the stars of retained vertices by intersection of old star with |c|
  for (bitset::RankFlags::iterator i = c.begin(); i(); ++i)
  {
    bitset::RankFlags st = parent.star(*i);
    st.slice(c); // extract bits set in |c| and repack to size |c.count()|
    d_star.push_back(st); // pack new stars into vector
  }

  // for the downedges we traverse those of |parent|, and see which apply
  for (std::vector<std::pair<Edge,Multiplicity> >::const_iterator
	 it=parent.d_downedge.begin(); it!=parent.d_downedge.end(); ++it)
  {
    size_t l = it->first.first;
    size_t s = it->first.second;
    if (c[l] and c[s]) // both long and short end points retained?
    d_downedge.push_back
      (std::make_pair(std::make_pair(c.position(l),c.position(s)),it->second));

  }
}

/******** accessors **********************************************************/


// recover Cartan matrix entry from Dynkin diagram
int DynkinDiagram::cartanEntry(size_t i,size_t j) const
{
  if (not star(i)[j])
    return i==j ? 2 : 0;
  for (size_t k=0; k<d_downedge.size(); ++k)
    if (d_downedge[k].first.first==i and d_downedge[k].first.second==j)
      return -int(d_downedge[k].second); // -2 or -3

  return -1; // simple edge, or labelled edge in short->long direction
}


/*!
  \brief Returns the connected component containing vertex \#i in the diagram.

  The algorithm is to start with |i|, and to construct "shells" from there, by
  taking each new shell to be the elements of the union of the stars of the
  old shell, that were not already considered. In fact since elements for the
  next shell are added into the current shell |newElts| while it is being
  traversed in the inner loop, we could advance by more than one edge during
  one traversal of the outer loop; if this happens, so much the better.
*/
bitset::RankFlags DynkinDiagram::component(size_t i) const
{
  bitset::RankFlags c;
  bitset::RankFlags newElts;

  for (newElts.set(i); newElts.any(); newElts.andnot(c)) {
    c |= newElts;
    for (bitset::RankFlags::iterator it = newElts.begin(); it(); ++it)
      newElts |= d_star[*it];
  }

  return c;
}


/*!
  \brief returns the set of end nodes of the graph.
*/
bitset::RankFlags DynkinDiagram::extremities() const
{
  bitset::RankFlags e;

  for (size_t i = 0; i < d_star.size(); ++i)
    if (d_star[i].count() <= 1)
      e.set(i);

  return e;
}


/*!
  \brief returns the labelled edge in the graph, assumed to be present.

  In fact it should be unique for valid uses, but the |assert| below should
  not prevent a relevant error to be thrown in erroneous cases
*/
Edge DynkinDiagram::labelEdge() const
{
  assert(d_downedge.size()>=1);
  return d_downedge[0].first;
}


/*!
  \brief returns the largest multiplicity in the graph.

  NOTE : we return 1 when there are no labelled edges, even when there are
  no edges at all!
*/
Multiplicity DynkinDiagram::maxMultiplicity() const
{
  Multiplicity m = 1;

  for (size_t i =0; i<d_downedge.size(); ++i)
    if (d_downedge[i].second>m)
      m = d_downedge[i].second;

  return m;
}


/*!
  Synopsis : returns the fork node of the graph.

  Precondition : the graph is irreducible, of type D or E;

  NOTE : returns exception value ~0 if the graph does not have a fork node.
  If the graph has more than one fork node, it will return the first of them.
*/
size_t DynkinDiagram::node() const
{
  // look for the first element whose star has more than two elements

  for (size_t i = 0; i < d_star.size(); ++i)
    if (d_star[i].count() > 2)
      return i;

  return ~0;
}

}  // |namespace dynkin|

/*****************************************************************************

        Chapter II -- Functions declared in dynkin.h

******************************************************************************/

namespace dynkin {


/*! \brief
  Returns the decomposition of |d| into connected components, a list of subsets

  We assure every vertex is in exactly one component, even if the (unlabelled)
  Dynkin graph should fail to be symmetric
*/
bitset::RankFlagsList components(const DynkinDiagram& d)
{
  bitset::RankFlagsList cl;

  bitset::RankFlags v;
  bitset::set(v,d.rank());

  while(v.any())
  {
    size_t i = v.firstBit();
    bitset::RankFlags c = d.component(i);
    c &= v; // assure no vertex is claimed by two components
    cl.push_back(c);
    v.andnot(c);
  }

  return cl;
}

/*! \brief

  Returns a permutation such that successive intervals of simple roots form
  connected components, numbered as needed for our Weyl group implementation.

  NOTE: the permutation result |pi| maps new index |i| to old index |pi[i]|.
*/
setutils::Permutation normalize(const DynkinDiagram& d)
{
  bitset::RankFlagsList cl = components(d);

  setutils::Permutation result= componentOrder(cl,d.rank());
  componentNormalize(result,cl,d,false);

  return result;
}

/*!
  Returns the (semisimple) Lie type of the Cartan matrix cm.
*/
lietype::LieType Lie_type(const latticetypes::LatticeMatrix& cm)
{

  DynkinDiagram d(cm);
  bitset::RankFlagsList cl = components(d);

  lietype::LieType result;
  result.reserve(cl.size());
  for (size_t i = 0; i < cl.size(); ++i)
  {
    DynkinDiagram cd(cl[i],d);
    result.push_back(lietype::SimpleLieType(irreducibleType(cd),cd.rank()));
  }

  return result;
}


/*! \brief
  Returns the (semisimple) Lie type of the Cartan matrix cm, also sets |pi|
  to the permutation from the standard ordering of simple roots for that type.

  Standard ordering is taken as Bourbaki ordering if |Bourbaki| holds,
  internal Weyl group implementation ordering otherwise. If |check| holds, a
  complete test is made of all entries in |cm|, throwing a |runtime_error| if
  it fails to be a valid Cartan matrix.
*/
lietype::LieType Lie_type(const latticetypes::LatticeMatrix& cm,
			  bool Bourbaki, bool check,
			  setutils::Permutation& pi)
{
  DynkinDiagram d(cm);
  bitset::RankFlagsList cl = components(d);

  // do the normalization as in normalize
  pi = componentOrder(cl,d.rank());
  lietype::LieType result = componentNormalize(pi,cl,d,Bourbaki);
  if (check)
  {
    for (size_t i=0; i<d.rank(); ++i)
      for (size_t j=0; j<d.rank(); ++j)
	if (cm(pi[i],pi[j])!=result.Cartan_entry(i,j))
	  throw std::runtime_error
	    ("Not a Cartan matrix of any semisimple type");
  }
  return result;

}
/*!
  Synopsis: Returns some permutation that will take |d| to Bourbaki form

  This means that nodes of the diagram |d| taken in the order |a[0],...,a[r-1]|
  traverse each of its connected components consecutively, and in the order
  prescribed by the the Bourbaki conventions for the type of that component
*/
setutils::Permutation bourbaki(const DynkinDiagram& d)
{
  bitset::RankFlagsList cl = components(d);

  // do the normalization as in normalize
  setutils::Permutation result = componentOrder(cl,d.rank());
  componentNormalize(result,cl,d,true);

  return result;
}





} // |namespace dynkin|

/*****************************************************************************

        Chapter III -- Auxiliary functions for this module

******************************************************************************/

namespace dynkin {
namespace {

/*!
  Returns a permutation such that the various components, listed in cl,
  are numbered by successive indices. The result maps these indices back
  to their original positions.
*/
setutils::Permutation componentOrder(const bitset::RankFlagsList& cl, size_t r)
{
  setutils::Permutation result; result.reserve(r);

  // traverse each component, write down its elements in sequence
  for (size_t i = 0; i<cl.size(); ++i)
    for (bitset::RankFlags::iterator it = cl[i].begin(); it(); ++it,--r)
      result.push_back(*it);

  assert (r==0); // check that correct rank was passed
  return result;
}


/*!
  Precondition : |a| contains a component ordering of d; cl contains the
  component list.

  Postcondition : |a| is modified so that the new permutation gives a
  normalized ordering on each component: it gives Bourbaki ordering unless
  |Bourbaki| is false and the type is BCD: then the order is reversed.

  The detected semisimple Lie type is returned.
*/
lietype::LieType componentNormalize(setutils::Permutation& a,
				    const bitset::RankFlagsList& cl,
				    const DynkinDiagram& d,
				    bool Bourbaki)
{
  size_t offset = 0;

  lietype::LieType result; result.reserve(cl.size());
  for (size_t i = 0; i < cl.size(); ++i) {

    // make a Dynkin diagram for the component
    DynkinDiagram cd(cl[i],d);

    // normalize it
    setutils::Permutation b;
    result.push_back(irreducibleNormalize(b,cd,Bourbaki));

    // piece together the permutation
    setutils::compose(a,b,offset);

    // update offset
    offset += cl[i].count();
  }
  return result;
}


/*!
  Precondition : d is an _irreducible_ Dynkin diagram;

  Postcondition : a holds a permutation which enumerates the vertices of
  d in an order that will induce a normal form of d;

  It is essentially a dispatching function for the various possible simple
  types.
*/
lietype::SimpleLieType
irreducibleNormalize(setutils::Permutation& a,
		     const DynkinDiagram& d,
		     bool Bourbaki)
{
  lietype::TypeLetter x = irreducibleType(d);

  switch (x) {
  case 'A':
    typeANormalize(a,d);
    break;
  case 'B':
    typeBNormalize(a,d,Bourbaki);
    break;
  case 'C':
    typeCNormalize(a,d,Bourbaki);
    break;
  case 'D':
    typeDNormalize(a,d,Bourbaki);
    break;
  case 'E':
    typeENormalize(a,d);
    break;
  case 'F':
  case 'f':
    typeFNormalize(a,d);
    break;
  case 'G':
  case 'g':
    typeGNormalize(a,d);
    break;
  default: // this should never happen!
    assert(false && "unexpected type in irreducibleNormalize");
    break;
  }
  return lietype::SimpleLieType(x,d.rank());
}


/*!
  Determines the (simple) type of a connected Dynkin diagram

  Precondition : d is connected (and therefore not empty)
*/
lietype::TypeLetter irreducibleType(const DynkinDiagram& d)
{
  const bitset::RankFlagsList& st = d.star();

  bitset::RankFlags extr = d.extremities();

  switch (extr.count())
  {

  case 1: // type is A1
    return 'A';

  case 2: // type is A,B,C,F or G
    switch (d.maxMultiplicity())
    {

    case 1: // type is A
      return 'A';

    case 2: { // type is B,C or F
      if (d.rank() == 2) // type is B
	return 'B';
      Edge e = d.labelEdge();
      if (extr.test(e.first)) // type is C
	return 'C';
      else if (extr.test(e.second)) // type is B
	return 'B';
      else // type is F
	return 'F';
    }

    case 3: // type is G
      return 'G';

    default: // should never happen
      return 0;
    }

  case 3:
    { // type is D or E
      size_t n = d.node();
      bitset::RankFlags sh = st[n]; // sh has three elements
      sh &= extr;           // now sh counts the number of short arms
      if (sh.count() == 1) // type is E
	return 'E';
      else
	return 'D';
    }

  default: // should never happen
    return 0;
  }
}


/*!
  Synopsis : puts in a a permutation that will enumerate d in Bourbaki order.

  Precondition : d is irreducible of type A, of rank >= 1;

  Postcondition : a holds a permutation which linearly enumerates the graph
  (which is a string in this case) --- there are exactly two such except in
  rank one;
*/
void typeANormalize(setutils::Permutation& a, const DynkinDiagram& d)
{
  size_t r = d.rank();
  a.resize(r);

  bitset::RankFlags e = d.extremities(); // e has one or two set elements
  a[0] = e.firstBit();

  bitset::RankFlags done;
  for (size_t i=1; i<r; ++i)
  {
    done.set(a[i-1]);
    a[i]= d.star(a[i-1]).andnot(done).firstBit();
  }
}


/*! \brief
  Puts in |a| a permutation that will enumerate |d| in linear order,
  ending with a labelled edge if |Bourbaki| holds, or starting if not.

  Precondition : d is irreducible of type B, of rank >= 2;

  There is a unique such ordering
*/
void typeBNormalize(setutils::Permutation& a,
		    const DynkinDiagram& d, bool Bourbaki)
{
  size_t r = d.rank();
  a.resize(r);

  size_t short_node = d.labelEdge().second;

  a[0] = Bourbaki ? d.extremities().reset(short_node).firstBit() : short_node;

  bitset::RankFlags done;
  for (size_t i=1; i<r; ++i)
  {
    done.set(a[i-1]);
    a[i]= d.star(a[i-1]).andnot(done).firstBit();
  }
}


/*! \brief
  Puts in |a| a permutation that will enumerate |d| in linear order,
  ending with a labelled edge if |Bourbaki| holds, or starting if not.

  Precondition : d is irreducible of type C, of rank >= 2;

  There is a unique such ordering
*/
void typeCNormalize(setutils::Permutation& a,
		    const DynkinDiagram& d, bool Bourbaki)
{
  size_t r = d.rank();
  a.resize(r);

  size_t long_node = d.labelEdge().first;

  a[0] = Bourbaki ? d.extremities().reset(long_node).firstBit() : long_node;

  bitset::RankFlags done;
  for (size_t i=1; i<r; ++i)
  {
    done.set(a[i-1]);
    a[i]= d.star(a[i-1]).andnot(done).firstBit();
  }
}


/*! \brief
  Puts in |a| a permutation that will enumerate |d| in alomst linear order
  (only the fork node has one neighbour at index distance 2 from it, which
  index is extremal); the fork node is at index |Bourbaki ? rank-3 : 2|.

  Precondition : d is irreducible of type D, with rank >= 4;
*/
void typeDNormalize(setutils::Permutation& a,
		    const DynkinDiagram& d, bool Bourbaki)
{
  size_t r = d.rank();
  a.resize(r);

  size_t fork = d.node();
  bitset::RankFlags st = d.star(fork);
  if (r==4)
  {
    bitset::RankFlags::iterator it = st.begin();
    a[0] = *it;
    a[1] = *++it;
    a[2] = fork;
    a[3] = *++it;
    if (Bourbaki)
      std::swap(a[1],a[2]); // fork node to position 1
    return;
  }

  st &= d.extremities(); // now |st| holds the short arms only

  bitset::RankFlags::iterator it = st.begin();
  bitset::RankFlags done;
  if (Bourbaki)
  {
    a[0] = d.extremities().andnot(st).firstBit();

    bitset::RankFlags done;
    for (size_t i=1; i<r-3; ++i)
    {
      done.set(a[i-1]);
      a[i]= d.star(a[i-1]).andnot(done).firstBit();
    }
    a[r-3] = fork;
    a[r-2] = *it;
    a[r-1] = *++it;
  }
  else // reverse Bourbaki
  {
    a[0] = *it;
    a[1] = *++it;
    a[2] = fork;

    bitset::RankFlags done=st;
    for (size_t i=3; i<r; ++i)
    {
      done.set(a[i-1]);
      a[i]= d.star(a[i-1]).andnot(done).firstBit();
    }
  }
}


/*!
  Synopsis : puts in a a permutation that will enumerate d in Bourbaki order.

  Precondition : d is irreducible of type E;

  Postcondition : a holds a permutation for which the node is in position
  3 (counting from 0), position 1 is the extremity of the branch of length
  1, position 0 is the extremity of a branch of length 2, position 2 is
  the other element of that branch, and the elements of the last branch are
  enumerated from the node. There are two solutions in type E6, one otherwise.
*/
void typeENormalize(setutils::Permutation& a, const DynkinDiagram& d)
{
  size_t r = d.rank();
  a.resize(r);

  size_t fork = d.node();
  a[3] = fork;

  bitset::RankFlags st = d.star(fork);
  bitset::RankFlags extr = d.extremities();

  bitset::RankFlags shortArms = extr;
  // this will make shortArms hold the short arm
  shortArms &= st;
  a[1] = shortArms.firstBit();

  st.andnot(shortArms);
  bitset::RankFlags::iterator it = st.begin();

  size_t x = *it;
  size_t y = *++it;

  bitset::RankFlags st_x = d.star(x);
  bitset::RankFlags st_y = d.star(y);

  bitset::RankFlags e_x = st_x;
  e_x &= extr;

  if (e_x.any()) // x is the origin of a branch of length 2
  {
    a[2] = x;
    a[0] = e_x.firstBit();
    a[4] = y;
  }
  else // y is the origin of a branch of length 2
  {
    a[2] = y;
    bitset::RankFlags e_y = st_y;
    e_y &= extr;
    a[0] = e_y.firstBit();
    a[4] = x;
  }

  // enumerate the last branch

  for (size_t i = 5; i < d.rank(); ++i)
  {
    bitset::RankFlags next = d.star(a[i-1]);
    next.reset(a[i-2]);
    a[i] = next.firstBit();
  }
}


/*!
  Synopsis : puts in a a permutation that will enumerate d in Bourbaki order.

  Precondition : d is irreducible of type F4;

  Postcondition : a holds a permutation which enumerates the graph in linear
  order, for which the middle edge is oriented from 1 to 2; such a permutation
  is unique.
*/
void typeFNormalize(setutils::Permutation& a, const DynkinDiagram& d)
{
  a.resize(4);

  Edge e = d.labelEdge();

  a[1] = e.first;
  a[2] = e.second;

  bitset::RankFlags st = d.star(a[1]);
  st.reset(a[2]);
  a[0] = st.firstBit();

  st = d.star(a[2]);
  st.reset(a[1]);
  a[3] = st.firstBit();
}


/*!
  Synopsis : puts in a a permutation that will enumerate d in Bourbaki order.

  Precondition : d is irreducible of type G2;

  Postcondition : a holds the permutation for which the edge is oriented from
  0 to 1; this is unique.
*/
void typeGNormalize(setutils::Permutation& a, const DynkinDiagram& d)
{
  a.resize(2);

  Edge e = d.labelEdge();

  a[0] = e.first;
  a[1] = e.second;
}

} // |namespace|
} // |namespace dynkin|
} // |namespace atlas|
