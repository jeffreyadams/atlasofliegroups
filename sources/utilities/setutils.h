/*!
\file
  This is setutils.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef SETUTILS_H  /* guard against multiple inclusions */
#define SETUTILS_H

#include <vector>
#include <algorithm>

#include "bitmap_fwd.h"

/******** type declarations *************************************************/

namespace atlas {

namespace setutils {

template<typename T>
  inline size_t find_index(const std::vector<T>& v, const T& x)
  { return std::find(v.begin(),v.end(),x)-v.begin(); }

struct Permutation
  : public std::vector<unsigned long>
  {
    typedef std::vector<unsigned long> Base;
    Permutation() : Base() {}                       // empty
    Permutation(unsigned long n) : Base(n) {}       // dimensioned only
    Permutation(unsigned long n, int unsused);      // identity
    Permutation(const Permutation& pi, int unused); // inverse
    template<typename I> Permutation(I b,I e) : Base(b,e) {} // range copy

  // right-compose with |*this|
    template<typename T> // any assignable type
      std::vector<T> pull_back(const std::vector<T>& v) const;

  // left-compose with |*this|
    template<typename U> // unsigned integral type
      std::vector<U> renumbering(const std::vector<U>& v) const;

    bitmap::BitMap renumbering(const bitmap::BitMap& b) const;

  // left-compose with |*this|, but allowing an exception value
    template<typename U>
      std::vector<U> renumbering(const std::vector<U>& v, U except) const;

  // left-multiply by |*this|; imperative version of |renumbering|
    template<typename U> void left_mult(std::vector<U>& v) const;

  // WARNING: has INVERSE interpretation of |*this| as |pull_back|:
  // right-compose with the inverse of the permutation (defining a left action)
  template<typename T> void permute(std::vector<T>& v) const;

  };

}

/******** function declarations **********************************************/

namespace setutils {

  void compose(Permutation&, const Permutation&, unsigned long n = 0);

  inline void identity(Permutation& p, unsigned long n) { p=Permutation(n,1); }

  inline void invert(Permutation& dst, const Permutation& src)
    { Permutation(src,-1).swap(dst); }

  template<typename U>
    Permutation standardize(const std::vector<U>& a, size_t bound);

}

}

/******** template definitions ***********************************************/

#include "setutils_def.h"

#endif
