/*!
\file
  This is setutils.cpp.  This file contains the non-template
  definitions of the functions declared in setutils.h

*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "setutils.h"
#include "bitmap.h"

/*
  This file contains the non-template definitions of the functions declared
  in setutils.h
*/

namespace atlas {

namespace setutils {


Permutation::Permutation(unsigned long n, int) // identity
  : Base(n)
{
  for (size_t i=n; i-->0; ) Base::operator[](i)=i;
}

Permutation::Permutation(const Permutation& pi, int) // inverse
  : Base(pi.size())
{
  for (size_t i=size(); i-->0; ) Base::operator[](pi[i])=i;
}

// Replace each index |i| in |v| by |pi[i]|, where |pi| is our permutation
template<typename U>
std::vector<U> Permutation::renumbering(const std::vector<U>& v) const
{
  const Permutation& pi=*this;
  std::vector<U> result; result.reserve(v.size());
  for (size_t i=0; i<v.size(); ++i)
    result.push_back(pi[v[i]]);
  return result;
}

/*!
  Synopsis: a *= b;

  Precondition : a holds a permutation of [0,N[; b holds a permutation of
  [0,M[; M + n <= N;

  Postcondition : a holds the permutation a_new of [0,N[ where a_new[i+n]
  = a[ b[i] + n ] for i in [0,M[; a is not changed outside the range
  [n,M+n[.

  NOTE : although we know that the permutation can be done in place with
  the aid of just a bitmap, we do the lazy approach here and make a copy
  of the range involved.
*/
void compose(Permutation& a, const Permutation& b, unsigned long n)
{
  std::vector<unsigned long> c(a.begin() + n, a.begin() + n + b.size());

  for (size_t j = 0; j < b.size(); ++j) {
    a[j+n] = c[b[j]];
  }
}

bitmap::BitMap Permutation::renumbering(const bitmap::BitMap& b) const
{
  bitmap::BitMap result(size());
  for (bitmap::BitMap::iterator it=b.begin(); it(); ++it)
    result.insert((*this)[*it]);

  return result;
}

/* Here we are again applying the permutation |p| to each of the entries
   of |v|, but the exceptional value of |except| is passed unchanged */
template<typename U>
std::vector<U>
   Permutation::renumbering(const std::vector<U>& v, U except) const
{
  const Permutation& pi=*this;
  std::vector<U> result; result.reserve(v.size());
  for (size_t i=0; i<v.size(); ++i)
    result.push_back(v[i]==except ? except : pi[v[i]]);
  return result;
}

// Replace each index |i| in |v| by |(*this)[i]|
template<typename U>
void Permutation::left_mult(std::vector<U>& v) const
{
  const Permutation& pi=*this;
  for (size_t i=0; i<v.size(); ++i)
    v[i]=pi[v[i]];
}

/*!
  Applies our permutation |pi| to the vector |v|. In other words, we send each
  entry v[i] to the new position v[pi[i]]; this means that afterwards for all
  |i|: |new_v[pi[i]]==old_v[i]|, or equivalently $new_v[i]=old_v[pi^{-1}[i]]$.

  Note that this is \emph{not} the same notion of permutation of entries used
  in the methods above, nor in the |permute| methods of matrices and bitsets;
  with respect to those we use the inverse permutation

  We are able to perform this permutation essentially in-place, using an
  auxiliary vector<bool>.
*/
template<typename T> void Permutation::permute(std::vector<T>& v) const
{
  assert(v.size()>=size());
  const Permutation& pi=*this;
  std::vector<bool> b(v.size(),false);

  for (unsigned long i = 0; i < v.size(); ++i) {
    if (b[i])
      continue;
    for (unsigned long j = pi[i]; j != i; j=pi[j])
    {
      std::swap(v[i],v[j]);
      b[j] = true;
    }
    b[i] = true;
  }
}

template <typename U>// unsigned type
  Permutation standardize(const std::vector<U>& a, size_t bound)
{
  std::vector<size_t> count(bound,0);
  for (size_t i=a.size(); i-->0; ) // downwards might be faster
  {
    assert(a[i]<bound);
    ++count[a[i]];
  }

  U sum=0;
  for (size_t i=0; i<bound; ++i) // cumulate
  {
    size_t ci=count[i]; count[i]=sum; sum+=ci;
  }
  // now |count[v]| holds number of values less than |v| in |a|

  Permutation result(a.size());
  for (size_t i=0; i<a.size(); ++i )
    result[i] = count[a[i]]++;

  return result;
}



// Instantiation of templates (only these are generated)

template std::vector<unsigned int>
Permutation::renumbering(const std::vector<unsigned int>& v) const; // blocks

template std::vector<unsigned long>
Permutation::renumbering(const std::vector<unsigned long>& v) const; // kgb

template void
Permutation::left_mult(std::vector<unsigned long>& v) const; // rootdata,weyl,.

template Permutation
standardize<unsigned int>(const std::vector<unsigned int>& a, size_t bound);

} // |namespace setutils|

} // |namespace atlas|
