/*!
\file
  This is poset.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef POSET_H  /* guard against multiple inclusions */
#define POSET_H

#include <vector>

#include "poset_fwd.h"

#include "bitmap.h"
#include "graph.h"

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

 public:

// constructors and destructors
  Poset() : d_below() {}

  explicit Poset(size_t n); // poset without any (nontrivial) comparabilities

  //! \brief Build Poset from its Hasse diagram
  explicit Poset(const std::vector<set::SetEltList>& hasse);

  //! \brief Build Poset from arbitrary list of links
  Poset(size_t n,const std::vector<Link>&);

  ~Poset() {}

// swap
  void swap(Poset& other) {
    d_below.swap(other.d_below);
  }

// accessors

  /*!
\brief The order relation itself.
  */

  bool lesseq(set::SetElt i, set::SetElt j) const{
    return i<j ? d_below[j].isMember(i) : i==j;
  }

  //! \brief Number of comparable pairs (including those on the diagonal)
  unsigned long Poset::n_comparable() const;

  void findMaximals(set::SetEltList&, const bitmap::BitMap&) const;
  set::SetEltList minima(const bitmap::BitMap&) const;

  /*!
\brief Size of the poset.
  */
  size_t size() const {
    return d_below.size();
  }

  void hasseDiagram(graph::OrientedGraph&) const;

  void hasseDiagram(graph::OrientedGraph&, set::SetElt) const;

// manipulators
  void resize(unsigned long);

  void extend(const std::vector<Link>&);
}; // class Poset

/*!
\brief Represents a poset by the symmetrization of the order relation
matrix.

The poset is {0,1,...,n-1}, and it is assumed that i less than j in
the poset implies i < j as integers.  Entry (j,i) is set if and only
if one of i and j is less than or equal to the other in the poset.
*/
class SymmetricPoset {

 private:

  /*!
\brief Rows of the symmetric poset BitMap.

Row \#j is a BitMap of size n (the size of the poset); bit i is set if
and only i and j are comparable in the poset.
  */
  std::vector<bitmap::BitMap> d_row;

 public:

// constructors and destructors
  SymmetricPoset() {}

  explicit SymmetricPoset(size_t n);

  explicit SymmetricPoset(const std::vector<set::SetEltList>&);

  ~SymmetricPoset() {}

// swap
  void swap(SymmetricPoset& other) {
    d_row.swap(other.d_row);
  }

// accessors

  /*!
\brief Row \#j of the symmetric poset matrix.

For i<j, bit i is set if and only if i is less than j in the poset.
Bit j is always set.  For i>j, bit j is set if and only if i is
greater than j in the poset.
  */
  const bitmap::BitMap& row(size_t j) const {
    return d_row[j];
  }

  /*!
\brief Size of the poset.
  */
  size_t size() const {
    return d_row.size();
  }

// manipulators
}; // class SymmetricPoset


/* ********************** Function(s) *****************************/

// memory-efficient version of |Poset(hasse).n_comparable()|
unsigned long n_comparable_from_Hasse
  (const std::vector<set::SetEltList>& hasse);

}

}

#endif
