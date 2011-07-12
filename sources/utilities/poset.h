/*!
\file
  This is poset.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef POSET_H  /* guard against multiple inclusions */
#define POSET_H

#include <vector>

#include "bitmap.h"
#include "tags.h"

#include "poset_fwd.h"
#include "graph_fwd.h"


namespace atlas {

namespace poset {



/******** type definitions **************************************************/

  /*!
\brief Represents a poset by the matrix of order relations.

  It is required that the ordering be compatible with the natural ordering
  on integers.

  */
class Poset {

 private:

/*!
\brief Matrix of order relations.

Bit i of d_below[j] is set if and only if |i| is less than |j| in the poset.
In other words, viewed as a set of integers, d_below[j]union{j} represents the
downwards closure in the poset of the singleton {j}.

By the assumption on the poset structure, the capacity of |d_below[j]| need
only be |j|.
*/
  std::vector<bitmap::BitMap> d_below;

  //! The basic method to add elementary relations
  void new_cover(unsigned long  x,unsigned long y) // add $x<y$, $y$ maximal
  { (d_below[y] |= d_below[x]).insert(x); }


 public:

// constructors and destructors
  Poset() : d_below() {}

  explicit Poset(size_t n); // poset without any (nontrivial) comparabilities

  //! \brief Build Poset from its Hasse diagram
  template<typename C> // |C| is some container of |set::SetElt|
  explicit Poset(const std::vector<C>& hasse);

  //! \brief Build Poset from arbitrary list of links
  Poset(size_t n,const std::vector<Link>&);

  Poset(const Poset& p, tags::DualTag);

  ~Poset() {}

// swap
  void swap(Poset& other) {
    d_below.swap(other.d_below);
  }

// accessors

  /*!
\brief The order relation itself.
  */

  bool lesseq(set::SetElt i, set::SetElt j) const
  { return i<j ? d_below[j].isMember(i) : i==j; }

  //! \brief Number of comparable pairs (including those on the diagonal)
  unsigned long n_comparable() const;

  void findMaximals(set::SetEltList&, const bitmap::BitMap&) const;
  set::SetEltList minima(const bitmap::BitMap&) const;

  /*!
\brief Size of the poset.
  */
  size_t size() const { return d_below.size(); }

  bool operator==(const Poset& other) const;

  void hasseDiagram(graph::OrientedGraph&) const;

  void hasseDiagram(graph::OrientedGraph&, set::SetElt) const;

// manipulators
  void resize(unsigned long);

/*!
\brief Transforms the poset into the weakest ordering containing the relations
  it previously contained, plus the relations |first < second| for all elements
  listed in |lk|.

  Precondition: |lk| is sorted in increasing lexicographical order, and is
  compatible with relations already present in the poset. More precisely the
  following (weaker, given the compatibility of the order relation with
  integral ordering) condition is assumed: all occurrences of a value |i| as
  first (smaller) member in a Link must follow all occurrences of |i| as
  second (larger) member in another Link, and no element smaller in a
  pre-existing relation should be the larger element of a link. This
  guarantees that the calls of |new_cover| generate the transitive closure.

*/
  void extend(const std::vector<Link>& lks)
  {
    for (size_t i=0; i<lks.size(); ++i)
      new_cover(lks[i].first,lks[i].second);
  }

  template<typename C> void new_max(C container)
  {
    size_t y=d_below.size();
    d_below.push_back(bitmap::BitMap(y));

    for (typename C::const_iterator
	   it=container.begin(); it!=container.end(); ++it)
      new_cover(*it,y);
  }

}; // class Poset


/* ********************** Function(s) *****************************/

// memory-efficient version of |Poset(hasse).n_comparable()|
unsigned long n_comparable_from_Hasse
  (const std::vector<set::SetEltList>& hasse);

/*!
\brief Constructs a Poset from its Hasse diagram.

  Precondition: it is assumed that for each |x|, |hasse(x)| is a container |C|
  listing elements covered by |x|, which elements must be numbers $<x$.

  As a consequence, the closure at |x| can be computed once |hasse(x)|
  is inspected, for increaing |x|.
*/
template<typename C>
Poset::Poset(const std::vector<C>& hasse)
: d_below(hasse.size())
{
  for (size_t i=0; i<hasse.size(); ++i)
  {
    d_below[i].set_capacity(i);
    for (typename C::const_iterator
	   it=hasse[i].begin(); it!=hasse[i].end(); ++it)
      new_cover(*it,i);
  }
}


} // namespace poset

} // namespace atlas

#endif
