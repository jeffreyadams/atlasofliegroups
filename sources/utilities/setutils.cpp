/*!
\file
  This is setutils.cpp.  This file contains the non-template
  definitions of the functions declared in setutils.h

*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "setutils.h"

/*
  This file contains the non-template definitions of the functions declared
  in setutils.h
*/

namespace atlas {

namespace setutils {

void compose(Permutation& a, const Permutation& b, unsigned long n)

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

{
  Permutation c(a.begin() + n, a.begin() + n + b.size());

  for (size_t j = 0; j < b.size(); ++j) {
    a[j+n] = c[b[j]];
  }

  return;
}

void identity(Permutation& a, unsigned long n)

/*!
  Synopsis: Sets a to the identity permutation of size n.
*/

{
  a.resize(n);

  for (unsigned long j = 0; j < n; ++j) {
    a[j] = j;
  }

  return;
}

void invert(Permutation& i, const Permutation& a)

/*!
  Synopsis: writes in i the inverse of a.
*/

{
  Permutation ai(a.size());

  for (unsigned long j = 0; j < a.size(); ++j) {
    ai[a[j]] = j;
  }

  i.swap(ai);

  return;
}

}

}
