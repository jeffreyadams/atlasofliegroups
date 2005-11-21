/*
  This is polynomials_def.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
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

  ... explain here when it is stable ...

 *****************************************************************************/

namespace polynomials {

template<typename C>
Polynomial<C>::Polynomial(Degree d)
  :d_data(d+1,0ul)

/*
  Synopsis: constructs x^d.

  We construct x^d, and not the zero polynomial, so that our basic assumption
  about the degree is satisfied.
*/

{
  d_data[d] = 1;
}

/******** accessors **********************************************************/

/******** manipulators *******************************************************/
template<typename C>
void Polynomial<C>::adjustSize()

/*
  Synopsis: adjust the size of d_data so that it corresponds to the degree + 1.
*/

{
  size_t d = d_data.size();

  for (size_t j = d; j;) {
    --j;
    if (d_data[j])
      break;
    else
      --d;
  }

  d_data.resize(d);

  return;
}

template<typename C>
Polynomial<C>& Polynomial<C>::safeAdd(const Polynomial& q, Degree d, C c)

/*
  Synopsis: adds x^d.c.q, watching for overflow

  NOTE: may forward a NumericOverflow exception.

  NOTE: we need to be careful in the case where q = *this, but we can
  avoid making a copy, by doing the addition top-to-bottom.
*/

{
  if (q.isZero()) // do nothing
    return *this;

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

  return *this;
}

template<typename C>
Polynomial<C>& Polynomial<C>::safeSubtract(const Polynomial& q, Degree d, C c)

/*
  Synopsis: subtracts x^d.c.q, watching for underflow

  NOTE: may forward a NumericUnderflow exception.

  NOTE: q = *this is possible only for d = 0; still, we do the prudent thing
  and subtract backwards.
*/

{
  if (q.isZero()) // do nothing
    return *this;

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

  return *this;
}

}

/*****************************************************************************

        Chapter II -- Functions declared in polynomials.h

  ... explain here when it is stable ...

 *****************************************************************************/

namespace polynomials {

template<typename C>
bool operator< (const Polynomial<C>& p, const Polynomial<C>& q)

/*
  Synopsis: polynomial comparison.

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

template<typename C>
bool operator== (const Polynomial<C>& p, const Polynomial<C>& q)

/*
  Synopsis: equality test.

  Explanation: p and q are equal if they are both zero, or they have same
  degree, and their corresponding coefficients are equal.
*/

{
  if (p.isZero())
    return q.isZero();

  if (p.degree() != q.degree())
    return false;

  for (size_t j = 0; j <= p.degree(); ++j)
    if (p[j] != q[j])
      return false;

  return true;
}

template<typename C> void safeAdd(C& a, C b)

/*
  Synopsis: a += b.

  Throws a NumericOverflow exception in case of overflow.
*/

{
  using namespace error;

  if (b > std::numeric_limits<C>::max() - a)
    throw NumericOverflow();
  else
    a += b;

  return;
}

template<typename C> void safeProd(C& a, C b)

/*
  Synopsis: a *= b.

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

  return;
}

template<typename C> void safeSubtract(C& a, C b)

/*
  Synopsis: a -= b.

  Throws a NumericUnderflow exception in case of underflow.
*/

{
  using namespace error;

  if (b > a)
    throw NumericUnderflow();
  else
    a -= b;

  return;
}

}

}
