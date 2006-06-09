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

class Poset {

 private:

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

  size_t size() const {
    return d_closure.size();
  }

  void hasseDiagram(graph::OrientedGraph&) const;

  void hasseDiagram(graph::OrientedGraph&, set::SetElt) const;

// manipulators
  void resize(unsigned long);

  void extend(const std::vector<Link>&);
};

class SymmetricPoset {

 private:

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
  const bitmap::BitMap& row(size_t j) const {
    return d_row[j];
  }

  size_t size() const {
    return d_row.size();
  }

// manipulators
};

}

}

#endif
