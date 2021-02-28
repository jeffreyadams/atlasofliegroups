/*
  This is permutations.h

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

  template<typename T,typename A>
  inline size_t find_index(const std::vector<T,A>& v, const T& x)
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

    // right-compose with |*this|: |result[i]=v[(*this)[i]]| for all |i|
    template<typename T,typename A> // here |T| is any assignable type
    std::vector<T,A> pull_back(const std::vector<T,A>& v) const;

    template<unsigned int n>
      bitset::BitSet<n> pull_back(const bitset::BitSet<n>& v) const;

    // left-compose with |*this|: |result[i]=(*this)[v[i]]| for all |i|
    template<typename U,typename A> // here |U| is an unsigned integral type
      std::vector<U,A> renumbering(const std::vector<U,A>& v) const;

    // set of values |(*this)[x]| for $x\in b$
    bitmap::BitMap renumbering(const bitmap::BitMap& b) const;

  // left-multiply by |*this|; imperative version of |renumbering|
    template<typename U,typename A> void renumber(std::vector<U,A>& v) const;

  // left-multiply by |*this|, but allowing an exception value
    template<typename U,typename A>
      void renumber(std::vector<U,A>& v, U except) const;

  // right-compose with the inverse of the permutation (defining a left action)
    template<typename T,typename A> void permute(std::vector<T,A>& v) const;

  // conjugate by pemutation matrix (to coordinates on permute(standard basis))
    template<typename T> void conjugate(matrix::Matrix_base<T>& M) const;

  // inverse conjugate by basis pemutation (to coordinates on (e_pi[i])_i)
    template<typename T> void inv_conjugate(matrix::Matrix_base<T>& M) const;

    bool is_negative() const; // whether negative sign
  };


/******** free standing function declarations *****************************/

  inline int sign(const Permutation& pi) { return pi.is_negative() ? -1 : 1; }

  bitmap::BitMap fixed_points(const Permutation& pi); // elements fixed

  void compose(Permutation&, const Permutation&, unsigned long n = 0);

  template<typename InputIt, typename Compare,
	   typename = typename std::enable_if<
  std::is_base_of<std::input_iterator_tag,
		  typename std::iterator_traits<InputIt>::iterator_category
  >::value>::type>
    Permutation standardization(InputIt first, InputIt last, Compare less);
  template<typename U,typename A> // here |U| is an unsigned integral type
  Permutation standardization(const std::vector<U,A>& a, size_t bound,
				std::vector<unsigned int>* stops = NULL);

  // a right action (like |pull_back|) on columns or rows of a matrix
  // produce matrix with columns $(M.column(pi[j]))_{j\in[0..pi.size()[}$
  template<typename T>
    void pull_back_columns(matrix::Matrix_base<T>& M, const Permutation& pi);

  // produce matrix with columns $(M.column(pi[j]))_{j\in[0..pi.size()[}$
  template<typename T>
    void pull_back_rows(matrix::Matrix_base<T>& M, const Permutation& pi);

 } // |namespace permutations|

} // |namespace atlas|

/********* template definitions with implicit instantiations  ******/

#include "permutations_def.h"

#endif
