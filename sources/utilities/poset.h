/*!
\file
  This is poset.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef POSET_H  /* guard against multiple inclusions */
#define POSET_H

#include <vector>

#include "poset_fwd.h"

#include "bitmap.h"
#include "graph.h"

/******** type definitions **************************************************/

namespace atlas {

namespace poset {

  /*!
\brief Represents a poset by the matrix of order relations.

DV: I do not know whether the poset order is assumed compatible with the
integer order (as it is in SymmetricPoset).
  */
class Poset {

 private:

  /*!
\brief Matrix of order relations.

Bit i of d_closure[j] is set if and only if i is less than or equal to
j in the poset.
  */
std::vector<bitmap::BitMap> d_closure;

 public:

// constructors and destructors
  Poset() {}

  explicit Poset(size_t n);

  explicit Poset(const std::vector<Link>&);

  ~Poset() {}

// swap
  void swap(Poset& other) {
    d_closure.swap(other.d_closure);
  }

// accessors
  void findMaximals(set::SetEltList&, const bitmap::BitMap&) const;

  /*!
\brief Size of the poset.
  */
  size_t size() const {
    return d_closure.size();
  }

  void hasseDiagram(graph::OrientedGraph&) const;

  void hasseDiagram(graph::OrientedGraph&, set::SetElt) const;

// manipulators
  void resize(unsigned long);

  void extend(const std::vector<Link>&);
};

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
};

}

}

#endif
