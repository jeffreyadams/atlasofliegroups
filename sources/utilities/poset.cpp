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
#include "graph.h"

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

bitmap::BitMap Poset::above(set::Elt x) const
{
  bitmap::BitMap result(size());

  for (size_t s=x+1; s<size(); ++s)
    if (d_below[s].isMember(x))
      result.insert(s);
  return result;
}

/*!
  Synopsis: writes in a the elements in b that are maximal for the induced
  order.

  Algorithm: the largest element x in b (if any) is maximal; add that to a,
  remove from b the intersection with the closure of x, and iterate.
*/
set::EltList Poset::maxima(const bitmap::BitMap& b) const
{
  unsigned long x=b.capacity();
  bitmap::BitMap result = b; // working copy; we shall remove covered elements

  while (result.back_up(x))
    result.andnot(d_below[x]); // remove below |x| (|bl[x]| remains set)

  return set::EltList(result.begin(),result.end()); // convert bitmap to vector
}

set::EltList Poset::minima(const bitmap::BitMap& b) const
{
  set::EltList result; // we shall produce a list of elements directly
  for (bitmap::BitMap::iterator it=b.begin(); it(); ++it)
    if (b.disjoint(d_below[*it])) // this certainly clears |t[n]|
      result.push_back(*it); // if intersection was empty, |n| is minimal

  return result;
}

// while dual to |covered_by|, this is harder (and untested)
set::EltList Poset::covers_of(set::Elt x) const
{
  set::EltList result; // we shall produce a list of elements directly
  bitmap::BitMap candidates = above(x);

  // in the next loop |candidates| is sieved while looping; this is allowed
  for (bitmap::BitMap::iterator it=candidates.begin(); it(); ++it)
  {
    result.push_back(*it); // record a covering element
    bitmap::BitMap::iterator jt=it; // copy iterator over |candidates|
    while ((++jt)()) // check remaining candidates
      if (d_below[*jt].isMember(*it)) // then |*jt| above |*it|, no |x|-cover
	candidates.remove(*jt); // so remove |*jt| from |candidates|
  }

  return result;
}


/*!
\brief Puts in h the Hasse diagram of the poset

  Explanation: the Hasse diagram is the oriented graph whose vertices are
  the elements of the poset, with an edge from each vertex to each vertex
  immediately below it.
*/
graph::OrientedGraph Poset::hasseDiagram() const
{
  graph::OrientedGraph h(size()); // empty graph of size of Poset

  for (set::Elt x = 0; x < size(); ++x)
    h.edgeList(x)=maxima(d_below[x]);// |x| is already absent from |d_below[x]|

  return h;
}

/* Puts in |h| the Hasse diagram of the downward closure of |max|. Since we do
   not want to renumber, all elements up to and including |max| will define
   vertices of the graph, but those not comparable to |max| will be isolated
   points. This (in particular the absence of incoming edges) is admittedly
   hard to detect; in practice 0 will be a least element, so that for other
   elements the absence of outgoing edges witnesses incomparability with |max|
 */
graph::OrientedGraph Poset::hasseDiagram(set::Elt max) const
{
  bitmap::BitMap cl(max+1);
  cl |= d_below[max]; cl.insert(max); // we must not forget |max| iself.

  graph::OrientedGraph h(max+1); // empty graph of size of |max+1|

  for (bitmap::BitMap::iterator it = cl.begin(); it(); ++it)
  {
    set::Elt x=*it;
    h.edgeList(x)=maxima(d_below[x]); // contained in |cl| by transitivity
  }
  return h;
}

unsigned long Poset::n_comparable() const
{
  unsigned long n=size(); // account for (absent) diagonal
  for (size_t i=0; i<size(); ++i)
    n+=d_below[i].size();

  return n;
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
  (const std::vector<set::EltList>& hasse)
{
  const size_t n=hasse.size();
  set::EltList min_after(n+1);
  min_after[n]=n;

  for (size_t i=n; i-->0;)
  {
    set::Elt min=min_after[i+1];
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
