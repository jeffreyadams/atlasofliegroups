/*
  This is polynomials_def.h.

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2009 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/*
  Template definitions for the class Polynomial.

  This module contains the implementation of the few things we need about
  polynomials in this program. The elementary operations are implemented in
  the |Polynomial| class template, while the |Safe_Poly| template adds "safe"
  arithmetic operations, that carefully check for overflow conditions,
  assuming that the coefficient type is an _unsigned_ type.

  As a class invariant, the leading coefficient stored is nonzero, so that the
  degree of the polynomial is simply deduced from the size of the coefficient
  vector: it is |size()-1|. As a corollary, the degree of the zero polynomial
  is -1 (for the unsigned |Degree| type). A consequence of this is that the
  subtraction operations have to watch out for drop in degree (not the
  addition operators, because coefficients are unsigned.)
*/


#include <limits>
#include <cassert>
#include <sstream>
#include <algorithm> // |std::fill| and |std::copy|

#include "error.h"

namespace atlas {

/*****************************************************************************

        Chapter I -- The Polynomial class

 *****************************************************************************/

namespace polynomials {

// this contructor assures mostly that |Polynomial<C>(C(0))| is null polynomial
template<typename C>
  Polynomial<C>::Polynomial(C c): d_data()
{
  if (c!=C(0))
    d_data.push_back(c);
}

/*!
  \brief Constructs cX^d.

  We construct cX^d, and not the zero polynomial, so that our basic assumption
  about the degree (leading coefficient is nonzero) is satisfied.
*/
template<typename C>
  Polynomial<C>::Polynomial(Degree d, C c): d_data(d+1,C(0))
{
  assert(c!=C(0));
  d_data[d] = c;
}

// shifted copy of an existing polynomial
template<typename C>
  Polynomial<C>::Polynomial(Degree d,const Polynomial& Q)
  : d_data(Q.isZero() ? 0 : Q.d_data.size()+d)
{
  if (Q.isZero())
    return; // avoid moving zeroes into an empty |d_data| array
  typename std::vector<C>::iterator bottom=d_data.begin()+d;
  std::fill(d_data.begin(),bottom,C(0)); // zero coefficient below bottom
  std::copy(Q.d_data.begin(),Q.d_data.end(),bottom); // remainder copied from Q
}

/******** accessors **********************************************************/

template<typename C>
const C& Polynomial<C>::operator[] (Degree i) const
{ assert(i<d_data.size());
  return d_data[i];
}


/******** manipulators *******************************************************/

template<typename C>
C& Polynomial<C>::operator[] (Degree i)
{ assert(i<d_data.size());
  return d_data[i];
}

/*!
  \brief Adjusts the size of d_data so that it corresponds to the degree + 1.

  Just casts away any leading zero coefficients, possibly all of them
*/

template<typename C>
void Polynomial<C>::adjustSize()
{
  size_t j = size();
  while (j-->0 and d_data[j]==C(0))
    d_data.pop_back(); // this is rare, so less expensive than always |resize|
}

template<typename C>
Polynomial<C>& Polynomial<C>::operator+= (const Polynomial& q)
{
  if (q.isZero())
    return *this; // avoid comparison of |&*q.d_data.end()| and |&q.d_data[0]|

  // make sure our polynomial has all the coefficients to be modified
  if (q.size() > size())
    resize(q.size());

  const C* src=&*q.d_data.end();
  C* dst=&d_data[q.size()];

  while (src>&q.d_data[0]) *--dst += *--src;
  // set degree
  adjustSize(); return *this;
}

template<typename C>
Polynomial<C>& Polynomial<C>::operator-= (const Polynomial& q)
{
  if (q.isZero())
    return *this;
  if (q.size() > size())
    resize(q.size());

  const C* src=&*q.d_data.end();
  C* dst=&d_data[q.size()];

  while (src>&q.d_data[0]) *--dst -= *--src;

  adjustSize(); return *this;
}

template<typename C>
Polynomial<C>& Polynomial<C>::subtract_from (const Polynomial& p)
{
  if (p.size() > size())
    resize(p.size());

  C* dst=&*d_data.end();
  while (dst>&d_data[p.size()])
  { --dst; *dst= -*dst; } // negate high coefficients

  // now |dst==&d_data[p.size()]|
  const C* src=&*p.d_data.end();
  while (dst-->&d_data[0]) { *dst= *--src - *dst; }

  adjustSize(); return *this;
}

template<typename C>
Polynomial<C>& Polynomial<C>::operator*= (C c)
{
  if (c==C(0)) *this=Polynomial();
  else
    for (C* p=&*d_data.end(); p>&d_data[0]; ) *--p *= c;

  return *this;
}

template<typename C>
Polynomial<C> Polynomial<C>::operator* (const Polynomial& q) const
{
  if (isZero() or q.isZero()) // we must avoid negative degrees!
    return Polynomial();
  Polynomial<C> result(degree()+q.degree(),C(1));
  // we have |result.size()==size()+q.size()-1|

  result.d_data[result.degree()]=C(0); // clear leading coefficient

  for (size_t j=size(); j-->0; )
  {
    C c=d_data[j];
    const C* src=&*q.d_data.end();
    C* dst=&result.d_data[q.size()+j]; // beyond end pointer on first iteration
    while (src>&q.d_data[0]) *--dst += c * *--src;
  }
  return result;
}

// rising degree division by $(1+cX)$, return new coefficient of $X^d$
template<typename C>
C Polynomial<C>::factor_by(C c)
{
  if (isZero())
    return C(0);
  C sum = d_data[0];
  for (Degree i=1; i<=degree(); ++i)
    sum = (d_data[i]-=c*sum);
  d_data[degree()]=0; // kill any remaining top coefficient, it is also in |sum|
  adjustSize();
  return sum;
}

template<typename C>
bool Polynomial<C>::multi_term () const
{
  if (isZero())
    return false;
  for (Degree i=degree(); i-->0; ) // skip leading term, which is nonzero
    if ((*this)[i]!=C(0))
      return true;
  return false;
}

/*!
  \brief Polynomial comparison: whether p < q

  Explanation: p < q if deg(p) < deg(q), or if degrees are equal, and the
  comparison holds for coefficients, in lex-order starting from the top.
*/
template<typename C>
bool compare (const Polynomial<C>& p, const Polynomial<C>& q)
{
  if (p.size() != q.size())
    return p.size() < q.size();

  // now p and q have same degree
  for (size_t j = p.size(); j-->0; )
    if (p[j] != q[j])
      return p[j]<q[j];

  // now p and q are equal
  return false;
}

/*
  Synopsis: prints out the monomial c.x^d.

  Preconditions: c is non-zero;

  Explanation: c and d are printed only if non-one, except in degree zero
  where c is always printed; x is the name of the indeterminate.

  The output format for the exponents is tex-like, but "q^13", not "q^{13}".
*/
template<typename C>
void printMonomial
       (std::ostream& strm, C c, polynomials::Degree d, const char* x)
{
  if (d == 0) // output c regardless
    strm << c;
  else
  {
    if (c<C(0) and c == C(-1)) // condition always false for unsigned types
      strm << '-';
    else if (c != C(1))
      strm << c;
    strm << x;
    if (d > 1)
      strm << "^" << d;
  }
}


template<typename C>
std::ostream& Polynomial<C>::print(std::ostream& strm, const char* x) const
{
  const Polynomial<C>& p = *this;
  std::ostringstream o; // accumulate in string for interpretation of width
  if (p.isZero())
    o << "0";
  else
    for (size_t i = p.size(); i-->0; )
      if (p[i]!=C(0)) // guaranteed true the first time
	printMonomial(i<p.degree() and p[i]>C(0) ? o<<'+' : o,p[i],i,x);

  return strm << o.str(); // now |strm.width()| is applied to whole polynomial
}


/*****************************************************************************

        Chapter II -- Safe polynomial arithmetic

 *****************************************************************************/


/*!
  \brief a += b.

  Throws a NumericOverflow exception in case of overflow.
*/
template<typename C> void safeAdd(C& a, C b)
{
  assert(a>=C(0)); // we're try to conserve this; it'd better be true intially
  assert(b>=C(0)); // so the we only need to check for overflow
  if (b > std::numeric_limits<C>::max() - a)
    throw error::NumericOverflow();
  else
    a += b;
}

/*!
  \brief a /= b.

  Throws a NumericOverflow exception in case of nondivisibility.
*/
template<typename C> void safeDivide(C& a, C b)
{
  if (a%b != 0) // safe use of |%|, since test is against |0|
    throw error::NumericOverflow();
  else
    a /= b; // now division is exact, so safe use of |/=|
}


/*!
  \brief a *= b.

  Throws a NumericOverflow exception in case of overflow.
*/
template<typename C> void safeProd(C& a, C b)
{
  assert(a>=C(0)); // we're try to conserve this; it'd better be true intially
  assert(b>=C(0)); // so the we only need to check for overflow
  if (a == 0) // do nothing
    return;

  if (b > std::numeric_limits<C>::max()/a)
    throw error::NumericOverflow();
  else
    a *= b;
}


/*!
  \brief a -= b.

  Throws a NumericUnderflow exception in case of underflow.
*/
template<typename C> void safeSubtract(C& a, C b)
{
  assert(a>=C(0)); // we're try to conserve this; it'd better be true intially
  assert(b>=C(0)); // so the we only need to check for underflow
  if (b > a)
    throw error::NumericUnderflow();
  else
    a -= b;
}


/*!
  \brief Adds x^d.c.q, to *this, watching for overflow, assuming |c>0|.

  NOTE: may forward a NumericOverflow exception.

  NOTE: we need to be careful in the case where q = *this, but we can
  avoid making a copy, by doing the addition top-to-bottom.
*/
template<typename C>
void Safe_Poly<C>::safeAdd(const Safe_Poly& q, Degree d, C c)
{
  if (q.isZero()) // do nothing
    return;

  size_t qs = q.size();   // save the original size of q, it might change

  // find degree of result
  if (q.size()+d > base::size())
    base::resize(q.size()+d);

  for (size_t j = qs; j-->0;)
  {
    C a = q[j];
    polynomials::safeProd(a,c);           // this may throw
    polynomials::safeAdd((*this)[j+d],a); // this may throw
  }
}

/* A simplified version avoiding multiplication in the common case |c==1| */

template<typename C>
void Safe_Poly<C>::safeAdd(const Safe_Poly& q, Degree d)
{
  if (q.isZero()) // do nothing
    return;

  size_t qs = q.size();   // save the original size of q, it might change

  // find degree
  if (q.size()+d > base::size())
    base::resize(q.size()+d);

  for (size_t j = qs; j-->0; )
    polynomials::safeAdd((*this)[j+d],q[j]); // this may throw
}

/*!

\brief Divides polynomial by scalar c, throwing an error is division is inexact
*/
template<typename C>
void Safe_Poly<C>::safeDivide(C c)
{
  for (size_t j = 0; j < base::size(); ++j )
    polynomials::safeDivide((*this)[j],c); //this may throw
}

/* Divides polynomials by $q+1$, imagining if necessary an additional leading
term \mu*q^{d+1} to make division exact, with appropriate scalar \mu and
d=(delta-1)/2 should be whole. We'll have delta==l(y)-l(x), whence the name.

Imagining such a term should be necessary only if the current degree of the
polynomial is precisely $d$, since the quotient must have non-negative
coefficients, and if it has a positive coefficient of $q^d$ then so does its
product by $q+1$.

However it could be that $\mu=0$, in which case the polynomial is already
divisible by $q+1$; then the quotient must have degree strictly less than the
half-integer $(delta-1)/2$. Whence the assertion |2*degree()<delta-1$.
*/
template<typename C>
void Safe_Poly<C>::safeQuotient(Degree delta)
{
  if (base::isZero()) // this avoids problems with |base::degree()|
    return; // need not and cannot invent nonzero \mu*q^{d+1} here
  for (size_t j = 1; j <= base::degree(); ++j)
    polynomials::safeSubtract((*this)[j],(*this)[j-1]); // does c[j] -= c[j-1]
  if ((*this)[base::degree()]==0) // number tested is the candidate for \mu
  { // polynomial was already multiple of q+1
    base::adjustSize(); // decreases degree by exactly 1
    assert(2*base::degree()+1<delta);
  }
  else
    assert(2*base::degree()+1==delta); // a term \mu*q^{d+1} is needed
}

/*!
  \brief Subtracts x^d.c.q from *this, watching for underflow, assuming |c>0|

  NOTE: may forward a NumericUnderflow exception.

  NOTE: q = *this is possible only for d = 0; still, we do the prudent thing
  and subtract backwards.
*/
template<typename C>
void Safe_Poly<C>::safeSubtract(const Safe_Poly& q, Degree d, C c)
{
  if (q.isZero()) // do nothing
    return;

  size_t qs = q.size();   // save the original size of q, it might change

  if (q.size()+d > base::size()) // underflow, leading coef becomes negative
    throw error::NumericUnderflow();

  for (size_t j = qs; j-->0; )
  {
    C a = q[j];
    polynomials::safeProd(a,c);                // this may throw
    polynomials::safeSubtract((*this)[j+d],a); // this may throw
  }

  // set degree
  base::adjustSize();
}

/* Again a simplified version deals with the common case |c==1| */

template<typename C>
void Safe_Poly<C>::safeSubtract(const Safe_Poly& q, Degree d)

{
  if (q.isZero()) // do nothing
    return;

  // save the degree of q
  Degree qs = q.size();

  if (qs+d > base::size()) // underflow
    throw error::NumericUnderflow();

  for (size_t j = qs; j-->0;)
    polynomials::safeSubtract((*this)[j+d],q[j]); // this may throw

  // set degree
  base::adjustSize();
 } // |safeSubtract|

} // namespace polynomials

} // namespace atlas
