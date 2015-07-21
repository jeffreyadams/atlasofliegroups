/*!
\file
  This is permutations.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef PERMUTATIONS_H  /* guard against multiple inclusions */
#define PERMUTATIONS_H

#include "permutations_fwd.h"

#include <vector>
#include <algorithm>

#include "bitset_fwd.h"
#include "bitmap_fwd.h"
#include "matrix_fwd.h"

/******** type declarations *************************************************/

namespace atlas {

namespace permutations {

template<typename T>
  inline size_t find_index(const std::vector<T>& v, const T& x)
  { return std::find(v.begin(),v.end(),x)-v.begin(); }

struct Permutation
  : public std::vector<unsigned long>
  {
    typedef std::vector<unsigned long> Base;
    Permutation() : Base() {}                       // empty
    Permutation(unsigned long n) : Base(n) {}       // dimensioned only
    Permutation(unsigned long n, int unused);       // identity
    Permutation(const Permutation& pi, int unused); // inverse
    template<typename I> Permutation(I b,I e) : Base(b,e) {} // range copy

  // right-compose with |*this|
    template<typename T> // any assignable type
      std::vector<T> pull_back(const std::vector<T>& v) const;

    template<size_t n>
      bitset::BitSet<n> pull_back(const bitset::BitSet<n>& v) const;

  // left-compose with |*this|
    template<typename U> // unsigned integral type
      std::vector<U> renumbering(const std::vector<U>& v) const;

    bitmap::BitMap renumbering(const bitmap::BitMap& b) const;

  // left-multiply by |*this|; imperative version of |renumbering|
    template<typename U> void renumber(std::vector<U>& v) const;

  // left-multiply by |*this|, but allowing an exception value
    template<typename U> void renumber(std::vector<U>& v, U except) const;

  // right-compose with the inverse of the permutation (defining a left action)
    template<typename T> void permute(std::vector<T>& v) const;

  // conjugate by pemutation matrix (to coordinates on permute(standard basis))
    template<typename T> void conjugate(matrix::Matrix_base<T>& M) const;

  // inverse conjugate by basis pemutation (to coordinates on (e_pi[i])_i)
    template<typename T> void inv_conjugate(matrix::Matrix_base<T>& M) const;

  };

}

/******** function declarations **********************************************/

namespace permutations {

  void compose(Permutation&, const Permutation&, unsigned long n = 0);

  template<typename U>
    Permutation standardization(const std::vector<U>& a, size_t bound,
				std::vector<unsigned int>* stops = NULL);

}

}

/******** template definitions ***********************************************/

#include "permutations_def.h"

#endif
