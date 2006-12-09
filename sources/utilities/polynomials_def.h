/*!
\file
\brief
  Template definitions for the class Polynomial.

  This module is intended for computations modulo some (small) integer. We
  have unmantled the "safe" constructions from Fokk's original implementation.

  Contrary to the original design, one should think of C as a class (of size
  smaller than that of long or even int), which provides for the arithmetic
  operations. Since modular arithmetic is the intended application, almost any
  operation can potentially make coefficients zero, unlike in the original
  implementation. (MvL)

  Our philsophy is that the degree of the polynomial is always deduced from
  the size of the coefficient vector: it is d_data.size()-1. As a corollary,
  the degree of the zero polynomial is -1 (for the unsigned Degree type);
  this is designated as UndefDegree. A consequence of this is that the
  addition and subtraction operations have to watch out for drop in degree.
  The main advantage of this is that the degree operator, which is used a lot, is trivial.
*/
/*
  This is polynomials_def.h.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include <limits>

#include "error.h"


namespace atlas {

/*****************************************************************************

        Chapter I -- The Polynomial class

 *****************************************************************************/

namespace polynomials {

template<typename C>
Polynomial<C>::Polynomial(Degree d)
  :d_data(d+1,C(0)) // all coefficients are C(0) for now

/*!
  \brief Constructs x^d.

  We construct x^d, and not the zero polynomial, so that our basic assumption
  about the degree is satisfied.
*/

{
  d_data[d] = C(1);
}

/******** accessors **********************************************************/

/******** manipulators *******************************************************/
template<typename C>
void Polynomial<C>::adjustSize()

/*!
  \brief Adjusts the size of d_data so that it corresponds to the degree + 1.

  Steps down through the coefficients of the polynomial, beginning
  with the top coefficient, until it finds a non-zero coefficient (or
  reaches zero).  Then resizes d_data so that this non-zero
  coefficient is the top one (or to have size zero if the polynomial is zero).
*/

{
  size_t j = d_data.size();

  while (j-->0 and d_data[j]==C(0)) { }

  d_data.resize(j+1);

  return;
}

template<typename C>
Polynomial<C>& Polynomial<C>::safeAdd(const Polynomial& q, Degree d, C c)

/*!
  \brief Adds x^d.c.q, to *this

  NOTE: we need to be careful in the case where q = *this, but we can
  avoid making a copy, by doing the addition top-to-bottom.
*/

{
  if (q.isZero()) // do nothing
    return *this;

  // save the degree of q
  Degree dq = q.degree();

  // find maximal possible degree
  if (q.d_data.size()+d > d_data.size())
    d_data.resize(q.d_data.size()+d,C(0));

  for (size_t j = dq+1; j-->0;) {
    C a = q[j];
    a*=c;
    d_data[j+d]+=a;
  }

  // set degree
  adjustSize();

  return *this;
}

template<typename C>
Polynomial<C>& Polynomial<C>::safeSubtract(const Polynomial& q, Degree d, C c)

/*!
  \brief Subtracts x^d.c.q from *this.

  This is now trivially delegated to safeAdd
*/

{ return safeAdd(q,d,-c);
}

}

/*****************************************************************************

        Chapter II -- Functions declared in polynomials.h

 *****************************************************************************/

namespace polynomials {

template<typename C>
bool compare (const Polynomial<C>& p, const Polynomial<C>& q)

/*!
  \brief Polynomial comparison.

  Explanation: p < q if deg(p) < deg(q), or if degrees are equal, and the
  comparison holds for coefficients, in lex-order starting from the top.
*/

{
  if (q.isZero())
    return false;

  // now q is nonzero
  if (p.isZero())
    return true;

  // now both p and q are nonzero
  if (p.degree() < q.degree())
    return true;
  if (q.degree() < p.degree())
    return false;

  // I have left one of Fokko's decreasing loops in the original style, MvL
  // now p and q have same degree
  for (size_t j = p.degree()+1; j;) {
    --j;
    if (p[j] != q[j])
      return p[j] < q[j];
  }

  // now p and q are equal
  return false;
}

template<typename C> void safeAdd(C& a, C b)

/*!
  \brief a += b.

*/

{ a += b;
}

template<typename C> void safeProd(C& a, C b)

/*!
  \brief a *= b.

*/

{
    a *= b;
}

template<typename C> void safeSubtract(C& a, C b)

/*!
  \brief a -= b.
*/

{
    a -= b;
}

}

}
