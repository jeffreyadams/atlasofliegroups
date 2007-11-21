/*!
\file
\brief
  Template definitions for the class Polynomial.

  This module contains the implementation of the few things we need
  about polynomials in this program. Apart from some basic elementary
  operations, the most notable things are the "safe" versions, that
  carefully check for overflow conditions. It is assumed that the
  coefficient type is an _unsigned_ type.

  Our philsophy is that the degree of the polynomial is always deduced from
  the size of the coefficient vector: it is d_data.size()-1. As a corollary,
  the degree of the zero polynomial is -1 (for the unsigned Degree type);
  this is designated as UndefDegree. A consequence of this is that the
  subtraction operations have to watch out for drop in degree (not the addition
  operators, because coefficients are unsigned.) The main advantage of this
  is that the degree operator, which is used a lot, is trivial.
*/
/*
  This is polynomials_def.h.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include <limits>

#include "error.h"

/*
  This module contains the implementation of the few things we need about
  polynomials in this program. Apart from some basic elementary operations,
  the most notable things are the "safe" versions, that carefully check for
  overflow conditions. It is assumed that the coefficient type is an
  _unsigned_ type.

  Our philsophy is that the degree of the polynomial is always deduced from
  the size of the coefficient vector: it is d_data.size()-1. As a corollary,
  the degree of the zero polynomial is -1 (for the unsigned Degree type);
  this is designated as UndefDegree. A consequence of this is that the
  subtraction operations have to watch out for drop in degree (not the addition
  operators, because coefficients are unsigned.) The main advantage of this
  is that the degree operator, which is used a lot, is trivial.
*/

namespace atlas {

/*****************************************************************************

        Chapter I -- The Polynomial class

 *****************************************************************************/

namespace polynomials {

template<typename C>
Polynomial<C>::Polynomial(Degree d)
  :d_data(d+1,0ul)

/*!
  \brief Constructs x^d.

  We construct x^d, and not the zero polynomial, so that our basic assumption
  about the degree is satisfied.
*/

{
  d_data[d] = 1;
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
}

template<typename C>
Polynomial<C>& Polynomial<C>::operator+= (const Polynomial& q)
{
  if (q.isZero()) return *this; // do nothing (not really necessary to check)

  // make sure d_data has all the coefficients to be modified
  if (q.d_data.size() > d_data.size())
    d_data.resize(q.d_data.size(),C(0));

  const C* src=&*q.d_data.end();
  C* dst=&d_data[q.d_data.size()];

  while (src>&q.d_data[0]) *--dst += *--src;
  // set degree
  adjustSize(); return *this;
}

template<typename C>
Polynomial<C>& Polynomial<C>::operator-= (const Polynomial& q)
{
  if (q.isZero()) return *this;
  if (q.d_data.size() > d_data.size()) d_data.resize(q.d_data.size(),C(0));

  const C* src=&*q.d_data.end(); C* dst=&d_data[q.d_data.size()];
  while (src>&q.d_data[0]) *--dst -= *--src;

  adjustSize(); return *this;
}

template<typename C>
Polynomial<C>& Polynomial<C>::subtract_from (const Polynomial& p)
{
  if (p.d_data.size() > d_data.size()) d_data.resize(p.d_data.size(),C(0));

  C* dst=&*d_data.end();
  while (dst>&d_data[p.d_data.size()])
  { --dst; *dst= -*dst; } // negate high coefficients

  // now |dst==&d_data[p.d_data.size()]|
  const C* src=&*p.d_data.end();
  while (dst-->&d_data[0]) { *dst= *--src - *dst; }

  adjustSize(); return *this;
}

template<typename C>
Polynomial<C>& Polynomial<C>::operator*= (C c)
{
  if (c==C(0)) *this=Polynomial();
  else
    for (const C* p=&*d_data.end(); p>&d_data[0]; ) *--p *= c;

  return *this;
}

template<typename C>
Polynomial<C> Polynomial<C>::operator* (const Polynomial& q) const
{
  if (isZero() or q.isZero()) // we must avoid negative degrees!
    return Polynomial();
  Polynomial<C> result(degree()+q.degree());
  result.d_data[result.degree()]=C(0); // clear leading coefficient

  for (size_t j=d_data.size(); j-->0; )
  {
    C c=d_data[j];
    const C* src=&*q.d_data.end();
    C* dst=&result.d_data[q.d_data.size()+j];
    while (src>&q.d_data[0]) *--dst += c * *--src;
  }
  return result;
}


/*!
  \brief Adds x^d.c.q, to *this, watching for overflow.

  NOTE: may forward a NumericOverflow exception.

  NOTE: we need to be careful in the case where q = *this, but we can
  avoid making a copy, by doing the addition top-to-bottom.
*/
template<typename C>
void Polynomial<C>::safeAdd(const Polynomial& q, Degree d, C c)
{
  if (q.isZero()) // do nothing
    return;

  // save the degree of q
  Degree dq = q.degree();

  // find degree
  if (q.d_data.size()+d > d_data.size())
    d_data.resize(q.d_data.size()+d,0ul);

  for (size_t j = dq+1; j;) {
    --j;
    C a = q[j];
    polynomials::safeProd(a,c);          // this may throw
    polynomials::safeAdd(d_data[j+d],a); // this may throw
  }
}

/* A simplified version avoiding multiplication in the common case |c==1| */

template<typename C>
void Polynomial<C>::safeAdd(const Polynomial& q, Degree d)
{
  if (q.isZero()) // do nothing
    return;

  // save the degree of q
  Degree dq = q.degree();

  // find degree
  if (q.d_data.size()+d > d_data.size())
    d_data.resize(q.d_data.size()+d,0ul);

  for (size_t j = dq+1; j;) {
    --j;
    polynomials::safeAdd(d_data[j+d],q[j]); // this may throw
  }
}


/*!
  \brief Subtracts x^d.c.q from *this, watching for underflow.

  NOTE: may forward a NumericUnderflow exception.

  NOTE: q = *this is possible only for d = 0; still, we do the prudent thing
  and subtract backwards.
*/
template<typename C>
void Polynomial<C>::safeSubtract(const Polynomial& q, Degree d, C c)

{
  if (q.isZero()) // do nothing
    return;

  // save the degree of q
  Degree dq = q.degree();

  if (dq+d > degree()) // underflow
    throw error::NumericUnderflow();

  for (size_t j = dq+1; j;) {
    --j;
    C a = q[j];
    polynomials::safeProd(a,c);               // this may throw
    polynomials::safeSubtract(d_data[j+d],a); // this may throw
  }

  // set degree
  adjustSize();
}

/* Again a simplified version deals with the common case |c==1| */

template<typename C>
void Polynomial<C>::safeSubtract(const Polynomial& q, Degree d)

{
  if (q.isZero()) // do nothing
    return;

  // save the degree of q
  Degree dq = q.degree();

  if (dq+d > degree()) // underflow
    throw error::NumericUnderflow();

  for (size_t j = dq+1; j;) {
    --j;
    polynomials::safeSubtract(d_data[j+d],q[j]); // this may throw
  }

  // set degree
  adjustSize();
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

  // now p and q have same degree
  for (size_t j = p.degree()+1; j;) {
    --j;
    if (p[j] < q[j])
      return true;
    if (q[j] < p[j])
      return false;
  }

  // now p and q are equal
  return false;
}

template<typename C> void safeAdd(C& a, C b)

/*!
  \brief a += b.

  Throws a NumericOverflow exception in case of overflow.
*/

{
  using namespace error;

  if (b > std::numeric_limits<C>::max() - a)
    throw NumericOverflow();
  else
    a += b;
}

template<typename C> void safeProd(C& a, C b)

/*!
  \brief a *= b.

  Throws a NumericOverflow exception in case of overflow.
*/

{
  using namespace error;

  if (a == 0) // do nothing
    return;

  if (b > std::numeric_limits<C>::max()/a)
    throw NumericOverflow();
  else
    a *= b;
}

template<typename C> void safeSubtract(C& a, C b)

/*!
  \brief a -= b.

  Throws a NumericUnderflow exception in case of underflow.
*/

{
  using namespace error;

  if (b > a)
    throw NumericUnderflow();
  else
    a -= b;
}

}

}
