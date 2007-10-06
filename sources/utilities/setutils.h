/*!
\file
  This is setutils.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef SETUTILS_H  /* guard against multiple inclusions */
#define SETUTILS_H

#include <vector>
#include <algorithm>

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

  // right-compose with |p|
  template<typename T>
  std::vector<T> pull_back(const std::vector<T>& v) const;

  // left-compose with |p|
  template<typename U>
  std::vector<U> renumber(const std::vector<U>& v) const;

  // left-compose with |p|, but allowing an exception value
  template<typename U>
  std::vector<U> renumber(const std::vector<U>& v, U except) const;

  // WARNING: this one has INVERSE interpretation of the permutation:
  template<typename T> void permute(std::vector<T>& v) const;

  };

}

/******** function declarations **********************************************/

namespace setutils {

  void compose(Permutation&, const Permutation&, unsigned long n = 0);

  inline void identity(Permutation& p, unsigned long n) { p=Permutation(n,1); }

  inline void invert(Permutation& dst, const Permutation& src)
    { Permutation(src,-1).swap(dst); }

}

}

/******** template definitions ***********************************************/

#include "setutils_def.h"

#endif
