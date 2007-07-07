/*!
\file
\brief Implementation of the class DynkinDiagram.
*/
/*
  This is dynkin.cpp
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "dynkin.h"

#include <cassert>

#include "latticetypes.h"
#include "lietype.h"

/*****************************************************************************

  Class and functions related to analysis of Dynkin diagrams

******************************************************************************/

namespace atlas {

namespace {

  using namespace dynkin;

  void componentOrder(setutils::Permutation&, const bitset::RankFlagsList&);
  void componentNormalize(setutils::Permutation&, const bitset::RankFlagsList&,
			  const DynkinDiagram&);
  void irreducibleNormalize(setutils::Permutation&, const DynkinDiagram&);
  lietype::TypeLetter irreducibleType(const DynkinDiagram&);
  void typeANormalize(setutils::Permutation&, const DynkinDiagram&);
  void typeBNormalize(setutils::Permutation&, const DynkinDiagram&);
  void typeCNormalize(setutils::Permutation&, const DynkinDiagram&);
  void typeDNormalize(setutils::Permutation&, const DynkinDiagram&);
  void typeENormalize(setutils::Permutation&, const DynkinDiagram&);
  void typeFNormalize(setutils::Permutation&, const DynkinDiagram&);
  void typeGNormalize(setutils::Permutation&, const DynkinDiagram&);

}

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
  , d_label()
{
  using namespace constants;

  for (size_t j = 0; j < c.numColumns(); ++j)
    for (size_t i = 0; i < c.numRows(); ++i)
      if (c(i,j) and (i != j)) {
	d_star[j].set(i);
	if (c(i,j) < -1) // only label multiple edges
	  d_label.insert(std::make_pair(Edge(i,j),-c(i,j)));
      }
}

/*!
  Constructs the restriction of |d| to the subset of the vertices flagged
  by |c|.
*/
DynkinDiagram::DynkinDiagram(const bitset::RankFlags& c,
			     const DynkinDiagram& d)
  : d_star()  // start with empty vector
  , d_label() // start without any labels
{
  d_star.reserve(c.count()); // vertices selected by |c|, renumbered from 0

  // get the stars of retained vertices by intersection of old star with |c|
  for (bitset::RankFlags::iterator i = c.begin(); i(); ++i) {
    bitset::RankFlags st = d.star(*i);
    st.slice(c); // extract bits set in |c| and repack to size |c.count()|
    d_star.push_back(st); // pack new stars into vector
  }

  // for the labels we just traverse the labels of |d|, and see which apply

  typedef std::map<Edge,Multiplicity>::const_iterator LI;

  LI label_end = d.d_label.end();

  for (LI i = d.d_label.begin(); i != label_end; ++i) {
    Edge e = i->first;
    if (c.test(e.first) and c.test(e.second)) // both end points retained?
    {
      // renumber the label indices according to bit positions in |c|
      e.first = c.position(e.first);
      e.second = c.position(e.second);
      d_label.insert(std::make_pair(e,i->second));
    }
  }
}

/******** accessors **********************************************************/

int DynkinDiagram::cartanEntry(size_t i,size_t j) const
{
  if (star(i).test(j))
  {
    std::map<Edge,Multiplicity>::const_iterator it=d_label.find(Edge(i,j));
    return it==d_label.end() ? -1 : -int(it->second);
  }
  else return i==j ? 2 : 0;
}


/*!
  Returns the component of vertex \#j in the diagram.

  The algorithm is to start with j, and to construct "shells" from there,
  by taking each new shell to be the elements of the union of the stars of
  the old shell, that were not already considered.
*/
bitset::RankFlags DynkinDiagram::component(size_t j) const
{
  using namespace bitset;

  RankFlags c;
  RankFlags newElts;

  for (newElts.set(j); newElts.any(); newElts.andnot(c)) {
    c |= newElts;
    for (RankFlags::iterator i = newElts.begin(); i(); ++i)
      newElts |= d_star[*i];
  }

  return c;
}

bitset::RankFlags DynkinDiagram::extremities() const

/*!
  Synopsis : returns the subset of the vertices which are extremities of the
  graph (i.e. whose star has zero or one element).
*/

{
  using namespace bitset;

  RankFlags e;

  for (size_t j = 0; j < d_star.size(); ++j)
    if (d_star[j].count() <= 1)
      e.set(j);

  return e;
}

Edge DynkinDiagram::labelEdge() const

/*!
  Synopsis : returns the multiple edge of d.

  Precondition : d is irreducible and not simply laced;

  NOTE : if d is not simply laced but not irreducible, will return the first
  multiple edge of d. There is guaranteed to be only one if d is irreducible.
*/

{
  std::pair<Edge,Multiplicity> p = *(d_label.begin());

  return p.first;
}

Multiplicity DynkinDiagram::maxMultiplicity() const

/*!
  Synopsis : returns the largest multiplicity in the graph.

  NOTE : we return 1 when there are no labelled edges, even when there are
  no edges at all!
*/

{
  Multiplicity m = 1;

  typedef std::map<Edge,Multiplicity>::const_iterator LI;

  LI label_end = d_label.end();

  for (LI i = d_label.begin(); i != label_end; ++i) {
    Multiplicity mi = i->second;
    if (mi > m)
      m = mi;
  }

  return m;
}

size_t DynkinDiagram::node() const

/*!
  Synopsis : returns the node of the graph.

  Precondition : the graph is irreducible, of type D or E;

  NOTE : behaviour is undefined if the graph does not have a node (the
  return value is -1 in this case). If the graph has more than one node,
  it will return the first of them.
*/

{
  // look for the first element whose star has more than two elements

  for (size_t j = 0; j < d_star.size(); ++j)
    if (d_star[j].count() > 2)
      return j;

  return ~0;
}

}

/*****************************************************************************

        Chapter II -- Functions declared in dynkin.h

  ... explain here when it is stable ...

******************************************************************************/

namespace dynkin {

void bourbaki(setutils::Permutation& a, const DynkinDiagram& d)

/*!
  Synopsis: writes in |a| some permutation that will take d to Bourbaki form

  This means that nodes of the diagram |d| taken in the order |a[0],...,a[r-1]|
  traverse each of its connected components consecutively, and in the order
  prescribed by the the Bourbaki conventions for the type of that component
*/

{
  using namespace lietype;

  a.resize(d.rank());
  bitset::RankFlagsList cl;

  // do the normalization as in normalize
  components(cl,d);
  componentOrder(a,cl);
  componentNormalize(a,cl,d);

  // examine components

  size_t r = 0;

  for (size_t j = 0; j < cl.size(); ++j) {

    // get type of component
    DynkinDiagram cd(cl[j],d);
    TypeLetter x = irreducibleType(cd);

    // reverse if type is BCD
    if (x == 'B' or x == 'C' or x == 'D') {
      size_t m = cl[j].count();
      for (size_t i = 0; i < m/2; ++i)
	std::swap(a[r+i],a[r+m-1-i]);
    }

    r += cl[j].count();
  }

  return;
}

void components(bitset::RankFlagsList& cl, const DynkinDiagram& d)

/*!
  Writes in cl the list of the connected components of d.
*/

{
  using namespace bitset;
  using namespace constants;

  cl.clear();

  RankFlags v;
  set(v,d.rank());

  for (; v.any();) {
    size_t j = v.firstBit();
    RankFlags c = d.component(j);
    cl.push_back(c);
    v.andnot(c);
  }

  return;
}

/*!
  Synopsis: writes in lt the Lie type of the Cartan matrix cm.
*/
void lieType(lietype::LieType& lt, const latticetypes::LatticeMatrix& cm)
{
  using namespace lietype;

  lt.clear();

  DynkinDiagram d(cm);
  bitset::RankFlagsList cl;

  components(cl,d);
  lt.reserve(cl.size());

  for (size_t j = 0; j < cl.size(); ++j) {
    DynkinDiagram cd(cl[j],d);
    TypeLetter x = irreducibleType(cd);
    SimpleLieType slt(x,cd.rank());
    lt.push_back(slt);
  }
}

lietype::LieType lieType(const latticetypes::LatticeMatrix& cm)
{
  lietype::LieType result; lieType(result,cm); return result;
}

void normalize(setutils::Permutation& a, const DynkinDiagram& d)

/*!
  Synopsis: writes in a a permutation that will take d to normal form

  This means that components of d are numbered consecutively, and that each
  component is normalized according to the Bourbaki conventions for the
  exceptional types, opposite-to-Bourbaki for the classical ones (this is the
  right ordering for the convenient implementation of Weyl groups.)
*/

{
  a.resize(d.rank());
  bitset::RankFlagsList cl;

  components(cl,d);
  componentOrder(a,cl);
  componentNormalize(a,cl,d);

  return;
}

}

/*****************************************************************************

        Chapter III -- Auxiliary functions for this module

  ... explain here when it is stable ...

******************************************************************************/

namespace {

void componentOrder(setutils::Permutation& a, const bitset::RankFlagsList& cl)

/*!
  Returns in a a permutation such that the various components, listed in cl,
  are numbered by successive indices.

  NOTE : it is always very confusing to choose between the permutation and
  its inverse. Our convention is that the _new_ vertex \#i is the _old_ vertex
  \# a[i].
*/

{
  using namespace bitset;

  // traverse each component, write down its elements in sequence

  a.clear();

  for (size_t j = 0; j < cl.size(); ++j)
    for (RankFlags::iterator i = cl[j].begin(); i(); ++i)
      a.push_back(*i);

  return;
}

void componentNormalize(setutils::Permutation& a,
			const bitset::RankFlagsList& cl,
			const DynkinDiagram& d)

/*!
  Precondition : a contains a component ordering of d; cl contains the
  component list.

  Postcondition : a is modified so that the new permutation gives a
  normalized ordering on each component, which is Bourbaki for the
  exceptional types, opposite-to-Bourbaki for the classical ones.
*/

{
  using namespace constants;
  using namespace setutils;

  size_t r = 0;

  for (size_t j = 0; j < cl.size(); ++j) {

    // make a Dynkin diagram for the component
    DynkinDiagram cd(cl[j],d);

    // normalize it
    Permutation b;
    irreducibleNormalize(b,cd);

    // piece together the permutation
    compose(a,b,r);

    // update r
    r += cl[j].count();
  }

  return;
}

void irreducibleNormalize(setutils::Permutation& a, const DynkinDiagram& d)

/*!
  Precondition : d is an _irreducible_ Dynkin diagram;

  Postcondition : a holds a permutation which enumerates the vertices of
  d in an order that will induce a normal form of d;

  It is essentially a dispatching function for the various possible simple
  types.
*/

{
  using namespace lietype;

  lietype::TypeLetter x = irreducibleType(d);

  switch (x) {
  case 'A':
    typeANormalize(a,d);
    break;
  case 'B':
    typeBNormalize(a,d);
    break;
  case 'C':
    typeCNormalize(a,d);
    break;
  case 'D':
    typeDNormalize(a,d);
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

  return;
}

lietype::TypeLetter irreducibleType(const DynkinDiagram& d)

/*!
  Synopsis : returns the type of the Dynkin diagram (one of the letters A-G).

  Precondition : d is irreducible;
*/

{
  using namespace bitset;
  using namespace constants;

  const RankFlagsList& st = d.star();

  RankFlags extr = d.extremities();

  switch (extr.count()) {

  case 1: // type is A1
    return 'A';

  case 2: // type is A,B,C,F or G
    switch (d.maxMultiplicity()) {

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

  case 3: { // type is D or E
    size_t n = d.node();
    RankFlags sh = st[n]; // sh has three elements
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

void typeANormalize(setutils::Permutation& a, const DynkinDiagram& d)

/*!
  Synopsis : puts in a a permutation that will enumerate d in Bourbaki order.

  Precondition : d is irreducible of type A, of rank >= 1;

  Postcondition : a holds a permutation which linearly enumerates the graph
  (which is a string in this case) --- there are exactly two such except in
  rank one;
*/

{
  using namespace bitset;
  using namespace constants;

  size_t r = d.rank();
  a.resize(r);

  RankFlags e = d.extremities(); // e has one or two set elements
  a[0] = e.firstBit();

  if (r == 1)
    return;

  RankFlags st = d.star(a[0]);  // st has a single element
  a[1] = st.firstBit();

  for (size_t j = 2; j < d.rank(); ++j) {
    RankFlags next = d.star(a[j-1]);
    next.reset(a[j-2]);
    a[j] = next.firstBit();
  }

  return;
}

void typeBNormalize(setutils::Permutation& a, const DynkinDiagram& d)

/*!
  Synopsis : puts in a a permutation that will enumerate d in opposite-to-
  Bourbaki order.

  Precondition : d is irreducible of type B, of rank >= 2;

  Postcondition : a holds a permutation which linearly enumerates the graph
  (which is a string in this case), with the label carried by the first edge,
  and the edge pointing toward the first vertex. This is unique.
*/

{
  using namespace bitset;
  using namespace constants;

  size_t r = d.rank();
  a.resize(r);

  Edge e = d.labelEdge();

  a[0] = e.second;
  a[1] = e.first;

  for (size_t j = 2; j < r; ++j) {
    RankFlags next = d.star(a[j-1]);
    next.reset(a[j-2]);
    a[j] = next.firstBit();
  }

  return;
}

void typeCNormalize(setutils::Permutation& a, const DynkinDiagram& d)

/*!
  Synopsis : puts in a a permutation that will enumerate d in opposite-to-
  Bourbaki order.

  Precondition : d is irreducible of type C, of rank >= 2;

  Postcondition : a holds a permutation which linearly enumerates the graph
  (which is a string in this case), with the label carried by the first edge,
  and the edge pointing toward the second vertex. This is unique.
*/

{
  using namespace bitset;
  using namespace constants;

  size_t r = d.rank();
  a.resize(r);

  Edge e = d.labelEdge();

  a[0] = e.first;
  a[1] = e.second;

  for (size_t j = 2; j < r; ++j) {
    RankFlags next = d.star(a[j-1]);
    next.reset(a[j-2]);
    a[j] = next.firstBit();
  }

  return;
}

void typeDNormalize(setutils::Permutation& a, const DynkinDiagram& d)

/*!
  Synopsis : puts in a a permutation that will enumerate d in opposite-to-
  Bourbaki order.

  Precondition : d is irreducible of type D, with rank >= 4;

  Postcondition : a holds a permutation for which the node is at position
  3 from the end, the two short branches are in the two last position, and
  the third branch is enumerated linearly downward from the node. There
  are two such enumerations when the rank is > 4, six when the rank is 4.
*/

{
  using namespace bitset;
  using namespace constants;

  size_t r = d.rank();
  a.resize(r);

  size_t n = d.node();
  a[2] = n;

  RankFlags st = d.star(n);
  RankFlags shortArms = d.extremities();
  // this will make shortArms hold the short arms of the diagram
  shortArms &= st;

  RankFlags::iterator i = shortArms.begin();
  a[0] = *i;
  a[1] = *(++i);

  // a[3] is the last element in st

  st.reset(a[0]);
  st.reset(a[1]);

  a[3] = st.firstBit();

  for (size_t j = 4; j < r; ++j) {
    RankFlags next = d.star(a[j-1]);
    next.reset(a[j-2]);
    a[j] = next.firstBit();
  }

  return;
}

void typeENormalize(setutils::Permutation& a, const DynkinDiagram& d)

/*!
  Synopsis : puts in a a permutation that will enumerate d in Bourbaki order.

  Precondition : d is irreducible of type E;

  Postcondition : a holds a permutation for which the node is in position
  3 (counting from 0), position 1 is the extremity of the branch of length
  1, position 0 is the extremity of a branch of length 2, position 2 is
  the other element of that branch, and the elements of the last branch are
  enumerated from the node. There are two solutions in type E6, one otherwise.
*/

{
  using namespace bitset;
  using namespace constants;

  size_t r = d.rank();
  a.resize(r);

  size_t n = d.node();
  a[3] = n;

  RankFlags st = d.star(n);
  RankFlags extr = d.extremities();

  RankFlags shortArms = extr;
  // this will make shortArms hold the short arm
  shortArms &= st;
  a[1] = shortArms.firstBit();

  st.andnot(shortArms);
  RankFlags::iterator i = st.begin();

  size_t x = *i;
  size_t y = *(++i);

  RankFlags st_x = d.star(x);
  RankFlags st_y = d.star(y);

  RankFlags e_x = st_x;
  e_x &= extr;

  if (e_x.any()) { // x is the origin of a branch of length 2
    a[2] = x;
    a[0] = e_x.firstBit();
    a[4] = y;
  }
  else { // y is the origin of a branch of length 2
    a[2] = y;
    RankFlags e_y = st_y;
    e_y &= extr;
    a[0] = e_y.firstBit();
    a[4] = x;
  }

  // enumerate the last branch

  for (size_t j = 5; j < d.rank(); ++j) {
    RankFlags next = d.star(a[j-1]);
    next.reset(a[j-2]);
    a[j] = next.firstBit();
  }

  return;
}

void typeFNormalize(setutils::Permutation& a, const DynkinDiagram& d)

/*!
  Synopsis : puts in a a permutation that will enumerate d in Bourbaki order.

  Precondition : d is irreducible of type F4;

  Postcondition : a holds a permutation which enumerates the graph in linear
  order, for which the middle edge is oriented from 1 to 2; such a permutation
  is unique.
*/

{
  using namespace bitset;
  using namespace constants;

  a.resize(4);

  Edge e = d.labelEdge();

  a[1] = e.first;
  a[2] = e.second;

  RankFlags st = d.star(a[1]);
  st.reset(a[2]);
  a[0] = st.firstBit();

  st = d.star(a[2]);
  st.reset(a[1]);
  a[3] = st.firstBit();

  return;
}

void typeGNormalize(setutils::Permutation& a, const DynkinDiagram& d)

/*!
  Synopsis : puts in a a permutation that will enumerate d in Bourbaki order.

  Precondition : d is irreducible of type G2;

  Postcondition : a holds the permutation for which the edge is oriented from
  0 to 1; this is unique.
*/

{
  a.resize(2);

  Edge e = d.labelEdge();

  a[0] = e.first;
  a[1] = e.second;

  return;
}

}

}
