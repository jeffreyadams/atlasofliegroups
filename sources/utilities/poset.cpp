/*
  This is poset.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#include "poset.h"

/*****************************************************************************

  This file contains a simple-minded implementation of a poset structure,
  where one explicitly keeps the bit-matrix of the order relation.

******************************************************************************/

/*****************************************************************************

        Chapter I -- The Poset class

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace poset {

/******** constructors and destructors ***************************************/

Poset::Poset(size_t n):d_closure(n,bitmap::BitMap(n))

/*
  Synopsis: constructs the discrete poset of size n.
*/

{
  for (size_t j = 0; j < n; ++j)
    d_closure[j].insert(j);
}

/******** accessors **********************************************************/

void Poset::findMaximals(set::SetEltList& a, const bitmap::BitMap& b) const

/*
  Synopsis: writes in a the elements in b that are maximal for the induced
  order.

  Algorithm: the largest element x in b (if any) is maximal; add that to a,
  remove from b the intersection with the closure of x, and iterate.
*/

{
  using namespace bitmap;

  BitMap bl = b;
  a.clear();

  while (not bl.empty()) {
    unsigned long x = bl.back()-1;
    a.push_back(x);
    bl.andnot(d_closure[x]);
  }

  return;
}

void Poset::hasseDiagram(graph::OrientedGraph& h) const

/*
  Synopsis: puts in h the Hasse diagram of the poset

  Explanation: the Hasse diagram is the oriented graph whose vertices are
  the elements of the poset, with an edge from each vertex to each vertex
  immediately below it.
*/

{
  using namespace bitmap;
  using namespace set;

  h.resize(size());

  for (SetElt x = 0; x < size(); ++x) {
    BitMap b = d_closure[x];
    b.remove(x);
    findMaximals(h.edges(x),b);
  }
}

void Poset::hasseDiagram(graph::OrientedGraph& h, set::SetElt max) const

/*
  Synopsis: puts in h the Hasse diagram of the closure of max.

  Explanation: the Hasse diagram is the oriented graph whose vertices are
  the elements of the poset, with an edge from each vertex to each vertex
  immediately below it.
*/

{
  using namespace bitmap;

  const BitMap& cl = d_closure[max];
  
  h.resize(size());

  for (BitMap::iterator i = cl.begin(); i != cl.end(); ++i) {
    BitMap b = d_closure[*i];
    b.remove(*i);
    findMaximals(h.edges(*i),b);
  }
  
  return;
}

/******** manipulators *******************************************************/
void Poset::resize(unsigned long n)

/*
  Synopsis: resizes the poset to size n, adding only the diagonal for the 
  new rows.
*/

{
  size_t prev = size();
  d_closure.resize(n);

  for (size_t j = 0; j < size(); ++j)
    d_closure[j].resize(n);

  for (size_t j = prev; j < size(); ++j)
    d_closure[j].insert(j);

  return;
}

void Poset::extend(const std::vector<Link>& lk)

/*
  Synopsis: transforms the poset into the weakest ordering stronger than
  the previous one, and for which the relations first < second for the elements
  in lk hold.

  Precondition: lk is sorted in lexicographical order.
*/

{
  for (size_t j = 0; j < lk.size(); ++j)
    d_closure[lk[j].second] |= d_closure[lk[j].first];

  return;
}

}

/*****************************************************************************

        Chapter I -- The SymmetricPoset class

  ... explain here when it is stable ...

******************************************************************************/

namespace poset {

/******** constructors and destructors ***************************************/
SymmetricPoset::SymmetricPoset(size_t n):d_row(n,bitmap::BitMap(n))

/*
  Synopsis: constructs the discrete poset of size n.
*/

{
  for (size_t j = 0; j < n; ++j)
    d_row[j].insert(j);
}

SymmetricPoset::SymmetricPoset(const std::vector<set::SetEltList>& hasse)
  :d_row(hasse.size(),bitmap::BitMap(hasse.size()))

/*
  Synopsis: constructs a symmetric poset from its hasse diagram.

  Precondition: it is assumed that for each x, hasse(x) contains the list
  of elements covered by x, and that the enumeration is such that those all
  have numbers < x.

  Algorithm: the poset part is easy: d_row[x] contains x and the union of
  the d_row[z], where z runs through hasse(x). The dual poset part is obtained
  by symmetrization; it could have been obtained by inverting the covering
  relations.
*/

{
  using namespace bitmap;
  using namespace set;

  // make poset part

  for (size_t x = 0; x < d_row.size(); ++x) {
    BitMap& b = d_row[x];
    b.insert(x);
    const SetEltList& h = hasse[x];
    for (size_t j = 0; j < h.size(); ++j)
      b |= d_row[h[j]];
  }

  // symmetrize; brute force for now

  for (size_t x = 0; x < d_row.size(); ++x) {
    BitMap::iterator b_end = d_row[x].end();
    for (BitMap::iterator i = d_row[x].begin(); i != b_end; ++i) {
      d_row[*i].insert(x);
    }
  }

}

}

}
