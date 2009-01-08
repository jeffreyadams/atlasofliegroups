/*!
\file
This is poset.cpp.  This file contains a simple-minded
implementation of a poset structure, where one explicitly keeps the
bit-matrix of the order relation.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
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

Poset::Poset(const Poset& p, tags::DualTag)
  : d_below()
{
  size_t n=p.size()-1;
  d_below.reserve(n+1);
  for (size_t i=0; i<=n; ++i)
  {
    d_below.push_back(bitmap::BitMap(i));
    bitmap::BitMap& below_i=d_below.back();
    for (size_t j=0; j<i; ++j)
      if (p.lesseq(n-i,n-j))
	below_i.insert(j); // in this case the new poset makes $j<i$
  }
}


/******** accessors **********************************************************/

bool Poset::operator==(const Poset& other) const
{
  if (size()!=other.size())
    return false;

  for (size_t i=0; i<size(); ++i)
    if (d_below[i]!=other.d_below[i])
    return false;

  return true;
}


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

  while (bl.back_up(x))
  {
    a.push_back(x);  // add currently maximal |x|
    bl.andnot(d_below[x]); // revove below |x| (|bl[x]| may safely remain set)
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

}

/*****************************************************************************

        Chapter I -- The SymmetricPoset class

******************************************************************************/

namespace poset {

/******** constructors and destructors ***************************************/

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
