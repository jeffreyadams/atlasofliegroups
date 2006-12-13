/*!
\file
\brief
  Template definitions for the class Polynomial.

  This module is intended for computations modulo some (small) integer. We
  have disabled the "safe" constructions from Fokko's original implementation.

  Contrary to the original design, one should think of C as a class (of size
  smaller than that of long or even int), which provides for the arithmetic
  operations. Since modular arithmetic is the intended application, almost any
  operation can potentially make coefficients zero, unlike in the original
  implementation. (MvL)

  Our philsophy is that the degree of the polynomial is always deduced from
  the size of the coefficient vector: it is d_data.size()-1. As a corollary,
  the degree of the zero polynomial is -1 (for the unsigned Degree type); this
  is designated as UndefDegree. A consequence of this is that the addition and
  subtraction operations have to watch out for drop in degree. The main
  advantage of this is that the degree operator, which is used a lot, is
  trivial.

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

/*!
  \brief Constructs x^d.

  We construct x^d, and not the zero polynomial, so that our basic assumption
  about the degree is satisfied.
*/

template<typename C>
Polynomial<C>::Polynomial(Degree d)
  :d_data(d+1,C(0)) // all coefficients are C(0) for now
{
  if(d+1!=0) // allow using this constructor with d == -1 to build Zero
    d_data[d] = C(1); // now we really have degree d
  adjustSize();       // small prostration to characteristic 1
}

/******** accessors **********************************************************/

/******** manipulators *******************************************************/

/*!
  \brief Adjusts the size of d_data so that it corresponds to the degree + 1.

  Steps down through the coefficients of the polynomial, beginning
  with the top coefficient, until it finds a non-zero coefficient (or
  reaches zero).  Then resizes d_data so that this non-zero
  coefficient is the top one (or to have size zero if the polynomial is zero).
*/

template<typename C>
void Polynomial<C>::adjustSize()
{
  size_t j = d_data.size();

  while (j-->0 and d_data[j]==C(0)) { }

  d_data.resize(j+1);

  return;
}


/*!
  \brief Adds x^d.c.q, to *this

  NOTE: we need to be careful in the case where q = *this, but we can
  avoid making a copy, by doing the addition top-to-bottom.
*/

template<typename C>
void Polynomial<C>::safeAdd(const Polynomial& q, Degree d, C c)
{
  if (q.isZero()) return; // do nothing

  // make sure d_data has all the coefficients to be modified
  if (q.d_data.size()+d > d_data.size())
    d_data.resize(q.d_data.size()+d,C(0));

  for (size_t j = q.degree()+1; j-->0;) {
    C a = q[j];
    a*=c;
    d_data[j+d]+=a;
  }

  // set degree
  adjustSize();
}

template<typename C>
void Polynomial<C>::safeAdd(PolRef<C> q, Degree d, C c)
{
  if (q.isZero()) return; // do nothing

  // make sure d_data has all the coefficients to be modified
  if (q.end()+d > d_data.size()) d_data.resize(q.end()+d, C(0) );

  for (size_t j = q.end(); j-->0;) d_data[j+d] += q[j]*c;

  // set degree
  adjustSize();
}


} // namespace polynomials

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

template<typename C>
inline void safeAdd(C& a, C b) { a *= b; }

template<typename C>
inline void safeSubtract(C& a, C b) { a -= b; }

}

}
