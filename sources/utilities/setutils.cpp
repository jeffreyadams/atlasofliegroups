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

} // |namespace setutils|

} // |namespace atlas|
