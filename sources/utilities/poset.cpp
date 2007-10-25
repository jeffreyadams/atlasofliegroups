/*!
\file
This is poset.cpp.  This file contains a simple-minded
implementation of a poset structure, where one explicitly keeps the
bit-matrix of the order relation.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "poset.h"

/*****************************************************************************

  This file contains a simple-minded implementation of a poset structure,
  where one explicitly keeps the bit-matrix of the order relation.

  However, only strictly sub-diagonal entries of the matrix are stored.

******************************************************************************/

/*****************************************************************************

        Chapter I -- The Poset class

******************************************************************************/

namespace atlas {

namespace poset {

/******** constructors and destructors ***************************************/


/*!
  Synopsis: constructs the discrete poset of size n.
*/
Poset::Poset(size_t n)
: d_below(n)
{
  for (size_t j = 0; j < n; ++j)
  {
    d_below[j].set_capacity(j); // clears all bits up to the diagaonal
  }
}

Poset::Poset(size_t n,const std::vector<Link>& lk)
: d_below(n)
{
  for (size_t j = 0; j < n; ++j)
  {
    d_below[j].set_capacity(j); // clears all bits up to the diagaonal
  }
  extend(lk);
}

/******** accessors **********************************************************/


unsigned long Poset::n_comparable() const
{
  unsigned long n=size(); // account for (absent) diagonal
  for (size_t i=0; i<size(); ++i)
    n+=d_below[i].size();

  return n;
}

/*!
  Synopsis: writes in a the elements in b that are maximal for the induced
  order.

  Algorithm: the largest element x in b (if any) is maximal; add that to a,
  remove from b the intersection with the closure of x, and iterate.
*/
void Poset::findMaximals(set::SetEltList& a, const bitmap::BitMap& b) const
{
  a.clear();
  bitmap::BitMap bl = b; // working copy of |b|
  unsigned long x=bl.capacity();

  while (bl.back_up(x)) {
    a.push_back(x);
    bl.andnot(d_below[x]); // don't care that |bl[x]| remains set
  }
  // in fact we could return |bl| as bitset here if that were asked for
}

set::SetEltList Poset::minima(const bitmap::BitMap& b) const
{
  set::SetEltList result;
  for (bitmap::BitMap::iterator it=b.begin(); it(); ++it)
  {
    size_t n=*it;
    bitmap::BitMap t=b;
    if (not(t&=d_below[n])) // this certainly clears |t[n]|
      result.push_back(n); // if intersection was empty, |n| is minimal
  }
  return result;
}


/*!
\brief Puts in h the Hasse diagram of the poset

  Explanation: the Hasse diagram is the oriented graph whose vertices are
  the elements of the poset, with an edge from each vertex to each vertex
  immediately below it.
*/
void Poset::hasseDiagram(graph::OrientedGraph& h) const
{
  h.resize(size());

  for (set::SetElt x = 0; x < size(); ++x) {
    const bitmap::BitMap& b = d_below[x]; // |x| is already absent from |b|
    findMaximals(h.edgeList(x),b);
  }
}

/*!
\brief Puts in |h| the Hasse diagram of the downward closure of |max|.

  Explanation: the Hasse diagram is the oriented graph whose vertices are
  the elements of the poset, with an edge from each vertex to each vertex
  immediately below it. Closure means all elements less than or equal to max.
*/
void Poset::hasseDiagram(graph::OrientedGraph& h, set::SetElt max) const


{
  bitmap::BitMap cl(max+1);
  cl |= d_below[max]; cl.insert(max); // we must not forget |max| iself.

  h.resize(size());

  for (bitmap::BitMap::iterator i = cl.begin(); i(); ++i) {
    set::SetElt x=*i;
    const bitmap::BitMap& b = d_below[x]; // |x| is already absent from |b|
    findMaximals(h.edgeList(x),b);
  }
}

/******** manipulators *******************************************************/

/*!
\brief Resizes the poset to size n, adding only the diagonal for the
  new rows.
*/
void Poset::resize(unsigned long n)
{
  size_t prev = size();
  d_below.resize(n);

  for (size_t j = prev; j < size(); ++j)
    d_below[j].set_capacity(j);
}

/*!
\brief Transforms the poset into the weakest ordering containing the relations
  it previously contained, plus the relations |first < second| for all elements
  listed in |lk|.

  Precondition: |lk| is sorted in increasing lexicographical order. More
  precisely the following weaker (given the compatibility of the order
  relation with integral ordering) condition is assumed: all occurrences of a
  value |i| as first (smaller) member in a Link must follow all occurrences of
  |i| as second (larger) member in another Link

*/
void Poset::extend(const std::vector<Link>& lk)
{
  for (size_t j = 0; j < lk.size(); ++j)
  {
    bitmap::BitMap& b=d_below[lk[j].second];
    b |= d_below[lk[j].first];
    b.insert(lk[j].first); // to compensate for unrepresented diagonal point
  }
}

}

/*****************************************************************************

        Chapter I -- The SymmetricPoset class

******************************************************************************/

namespace poset {

/******** constructors and destructors ***************************************/

/*!
\brief Constructs the discrete poset (no relations) of size n.
*/
SymmetricPoset::SymmetricPoset(size_t n):d_row(n,bitmap::BitMap(n))
{
  for (size_t j = 0; j < n; ++j)
    d_row[j].insert(j);
}


/*!
\brief Constructs a symmetric poset from its hasse diagram.

  Precondition: it is assumed that for each x, hasse(x) contains the list
  of elements covered by x, and that the enumeration is such that those all
  have numbers < x.

  Algorithm: the poset part is easy: d_row[x] contains x and the union of
  the d_row[z], where z runs through hasse(x). The dual poset part is obtained
  by symmetrization; it could have been obtained by inverting the covering
  relations.
*/
SymmetricPoset::SymmetricPoset(const std::vector<set::SetEltList>& hasse)
  :d_row(hasse.size(),bitmap::BitMap(hasse.size()))
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

  // symmetrize; set bits of row |x| in column |x| (brute force for now)

  for (size_t x = 0; x < d_row.size(); ++x)
    for (BitMap::iterator i = d_row[x].begin(); i(); ++i) {
      d_row[*i].insert(x);
  }

}

unsigned long n_comparable_from_Hasse
  (const std::vector<set::SetEltList>& hasse)
{
  const size_t n=hasse.size();
  set::SetEltList min_after(n+1);
  min_after[n]=n;

  for (size_t i=n; i-->0;)
  {
    set::SetElt min=min_after[i+1];
    for (size_t j=0; j<hasse[i].size(); ++j)
      if (hasse[i][j]<min)
	min=hasse[i][j];
    min_after[i]=min;
  }

  std::vector<bitmap::BitMap> closure(n);
  unsigned long count=0;

  for (size_t i=0; i<n; ++i)
  {
    bitmap::BitMap& cl=closure[i];
    cl.set_capacity(i+1); cl.insert(i);
    for (size_t j=0; j<hasse[i].size(); ++j)
      cl |= closure[hasse[i][j]];

    count+=cl.size();

    // now free memory of bitmaps that will no longer be needed
    for (size_t j=min_after[i]; j<min_after[i+1]; ++j)
      closure[j].set_capacity(0);
  }

  return count;
}


} // namespace poset

} // namespace atlas
