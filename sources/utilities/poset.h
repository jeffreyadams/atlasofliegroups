/*
  This is poset.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

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

/*
  A |Poset| represents a poset as the matrix of order relations, but using
  a compact representation of the matrix.

  It is required that the ordering be compatible with the natural ordering
  on integers, so the matrix is upper triangular.
*/
class Poset {

 public:
  typedef unsigned int Elt; // should suffice for any realistic Poset
  typedef std::vector<Elt> EltList;

  // a type used to represent a poset relation in arguments of certain methods
  typedef std::pair<Elt,Elt> Link;

 private:
/*
  List of bitsets representing (upper part of) matrix of order relations.

  Bit i of d_below[j] is set if and only if |i| is less than |j| in the poset.
  In other words, viewed as a set of integers, d_below[j]union{j} represents
  the downwards closure in the poset of the singleton {j}.

  By the assumption on the poset structure, the capacity of |d_below[j]| need
  only be |j|.
*/
  std::vector<bitmap::BitMap> d_below;

 public:

// constructors and destructors
  Poset() : d_below() {}

  explicit Poset(Elt n); // poset without any (nontrivial) comparabilities

  // Construct a |Poset| from its Hasse diagram
  template<typename C> // |C| is some container of |set::Elt|
    explicit Poset(const std::vector<C>& hasse);

  // Build a |Poset| from arbitrary list of (generating) links
  Poset(Elt n,const std::vector<Link>&);

  Poset(const Poset& p, tags::DualTag);

  ~Poset() {}

// swap
  void swap(Poset& other) { d_below.swap(other.d_below); }

// accessors

  // The order relation itself.
  bool lesseq(Elt i, Elt j) const
  { return i<j ? d_below[j].isMember(i) : i==j; }

  bool operator==(const Poset& other) const;

  Elt size() const { return d_below.size(); }

  const bitmap::BitMap& below(Elt y)const { return d_below[y]; }
  bitmap::BitMap above(Elt x) const;

  bitmap::BitMap maxima(const bitmap::BitMap&) const;
  bitmap::BitMap minima(const bitmap::BitMap&) const;

  bitmap::BitMap covered_by(Elt y) const { return maxima(below(y)); }
  bitmap::BitMap covers_of(Elt x) const { return minima(above(x)); }

  graph::OrientedGraph hasseDiagram() const; // full Hasse diagram
  graph::OrientedGraph hasseDiagram(Elt max) const; // part |<=max|

  // Number of comparable pairs (including those on the diagonal)
  unsigned long n_comparable() const;

// manipulators
  void resize(unsigned long);

/*
  Transform the poset into the weakest ordering containing the relations
  it previously contained, plus the relations |first < second| for
  all elements listed in |lks|.

  Precondition: |lks| is sorted in increasing lexicographical order, and is
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
    for (auto it=lks.begin(); it!=lks.end(); ++it)
      new_cover(it->first,it->second);
  }

  // add a new maximal element, comparable with elements in |container|
  template<typename C> void new_max(C container)
  {
    Elt y=d_below.size();
    d_below.push_back(bitmap::BitMap(y));

    for (typename C::const_iterator
	   it=container.begin(); it!=container.end(); ++it)
      new_cover(*it,y);
  }

 private:
  // The basic method to add elementary relations
  void new_cover(unsigned long  x,unsigned long y) // add $x<y$, $y$ maximal
  { (d_below[y] |= d_below[x]).insert(x); }

}; // class Poset


/* ********************** Function(s) *****************************/

// memory-efficient version of |Poset(hasse).n_comparable()|
unsigned long n_comparable_from_Hasse
  (const std::vector<Poset::EltList>& hasse);

/*
  Construct a Poset from its Hasse diagram.

  Precondition: it is assumed that for each |x|, |hasse(x)| is a container |C|
  listing elements covered by |x|, which elements must be numbers $<x$.

  As a consequence, for increasing |x| the closure at |x| can be computed
  once |hasse(x)| is inspected.
*/
template<typename C>
  Poset::Poset(const std::vector<C>& hasse)
: d_below(hasse.size())
{
  for (Elt i=0; i<hasse.size(); ++i)
  {
    d_below[i].set_capacity(i);
    for (typename C::const_iterator
	   it=hasse[i].begin(); it!=hasse[i].end(); ++it)
      new_cover(*it,i);
  }
}


} // |namespace poset|

} // |namespace atlas|

#endif
