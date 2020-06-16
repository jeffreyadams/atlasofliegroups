/*
  This is polynomials_def.h.

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2009,2020 Marc van Leeuwen
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
#include <stdexcept>

#include "error.h" // for ||

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

/*
  Construct $cX^d$.

  The constructor builds a mononomial $cX^d$, and not the zero polynomial, so
  that our invariant about the degree (the leading coefficient is nonzero)
  is satisfied.
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

/*
  Adjust the size of |d_data| so that it corresponds to the degree + 1.

  Just cast away any leading zero coefficients, possibly all of them
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
Polynomial<C>& Polynomial<C>::operator/= (C c)
{
  if (c==C(0))
    throw std::runtime_error("Polynomial division by 0");
  for (C* p=&*d_data.end(); p>&d_data[0]; )
    if (*--p%c==C(0))
      *p /= c;
    else
      throw std::runtime_error("Inexact polynomial integer division");

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

// rising degree division by $(1+cX)$, return remainder in coefficient of $X^d$
template<typename C> C Polynomial<C>::up_remainder(C c, Degree d) const
{
  if (isZero())
    return C(0);
  assert(degree()<=d);
  C remainder = d_data[0];
  for (Degree i=1; i<=d; ++i)
    remainder = coef(i) - c*remainder;
  return remainder; // excess that was found in |coef(d)| for exact division
}

// version of the same, with $c=1$, that also changes |*this| into the quotient
template<typename C> C Polynomial<C>::factor_by_1_plus_q(Degree d)
{
  if (isZero())
    return C(0);
  assert(this->degree()<=d);
  if (d>this->degree())
    d_data.resize(d+1,C(0)); // extend with zero coefficients to make degree $d$
  C remainder = d_data[0];
  for (Degree i=1; i<=d; ++i) // |d_data[i-1]| becomes quotient coefficient
    remainder = (d_data[i] -= remainder); // |d_data[i]| is now remainder
  d_data[d]=0; // kill |remainder| that is now in the top coefficient
  adjustSize(); // and then reduce quotient size to match its actual degree
  return remainder; // excess that was found in |d_data[d]| for exact division
}

// version of the preceding with division by $1+q^k$ rather than by $1+q$
// desired behavior: divide this by 1+q^k, keeping just the terms up to
// degree d-k; last k correspond to remainder. (It might be more useful to have
// the return value the remainder, but that would require returning type
// |Polynomial<C>|, which in the current application in ext_kl is ignored.)
template<typename C>
C Polynomial<C>::factor_by_1_plus_q_to_the(Degree k,Degree d)
{
  if (isZero())
    return C(0);
  assert(this->degree()<=d);
  if (d>this->degree())
    d_data.resize(d+1,C(0)); // extend with zero coefficients to make degree $d$
  for (Degree i=k; i<=d; ++i) // quotient captures |d_data[i-k]|
    d_data[i] -= d_data[i-k]; // and term in degree |i| remains to be treated
  C result = d_data[d];
  for (Degree i=k; i-->0;)
    d_data[d-i]=0; // kill off remainder and lower degree terms
  adjustSize();
  return result;
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

/*
  Polynomial comparison: whether $p < q$

  Explanation: $p < q$ if $\deg(p) < \deg(q)$, or if degrees are equal, and the
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
  Print out the monomial $c.x^d$.

  Preconditions: |c| is non-zero;

  Explanation: |c| and |d| are printed only if non-one, except in degree zero
  where |c| is always printed; |x| is the name of the indeterminate.

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


/*
  Perform |a += b| after testing that it will not overflow capacity of type |C|

  Throw a |NumericOverflow| exception in case that test fails
*/
template<typename C> void safe_add(C& a, C b)
{ // start with some tests that are useful only for signed types |C|
  if (b > std::numeric_limits<C>::max() - a)
    throw error::NumericOverflow();
  a += b;
}

// substraction throws no runtime error; success should be statically ensured
template<typename C> void safe_subtract(C& a, C b)
{
  assert(a>=C(0)); // we're try to conserve this; it'd better be true initially
  assert(b>=C(0)); // so the we only need to check for underflow
  assert(a>=b); // this is what caller should ensure mathematically
  a -= b;
}

/*
  Perform |a *= b| after testing that it will not overflow capacity of type |C|

  Throw a |NumericOverflow| exception in case that test fails
*/
template<typename C> void safe_multiply(C& a, C b)
{
  if (a == C(0))
    return; // do nothing

  if (b > std::numeric_limits<C>::max()/a)
    throw error::NumericOverflow();
  a *= b;
}

template<typename C> void safe_divide(C& a, C b)
{
  assert(b!= C(0));
  assert(a%b == C(0));
  a /= b; // now division is exact, so safe use of |/=|
}



/*
  Add $x^d.c.q$, to |*this|, watching for overflow, assuming |c>0|.

  We need to be careful in the case where |q| aliasses |*this|, but we can
  avoid making a copy, by doing the addition top-to-bottom.
*/
template<typename C>
void Safe_Poly<C>::safeAdd(const Safe_Poly& q, Degree d, C c)
{
  if (q.isZero() or c==C(0))
    return; // do nothing

  const auto qs = q.size();   // save the original size of q, it might change

  // ensure sufficient room for result
  if (qs+d > base::size())
    base::resize(qs+d);

  auto* shift_base = &(*this)[d]; // pointer to lowest coefficient affected
  for (size_t j = qs; j-->0;) // loop must be backwards if |&q == this|
  {
    C a = q[j]; // need a copy for the multiplication
    polynomials::safe_multiply(a,c);        // this may throw
    polynomials::safe_add(shift_base[j],a); // this may throw
  }
} // |safeAdd|

/* A simplified version avoiding multiplication in the common case |c==1| */

template<typename C>
void Safe_Poly<C>::safeAdd(const Safe_Poly& q, Degree d)
{
  if (q.isZero()) // do nothing
    return;

  const auto qs = q.size();   // save the original size of q, it might change

  // find degree
  if (qs+d > base::size())
    base::resize(qs+d);

  auto* shift_base = &(*this)[d]; // pointer to lowest coefficient affected
  for (size_t j = qs; j-->0; ) // loop must be backwards if |&q == this|
    polynomials::safe_add(shift_base[j],q[j]); // this may throw
} // |safeAdd|

/*
  Subtract $x^d.c.q$ from |*this|, asserting no underflow, assuming |c>0|

  We need to be careful in the case where |q| aliasses |*this|, but we can
  avoid making a copy, by doing the addition top-to-bottom.
*/
template<typename C>
void Safe_Poly<C>::safeSubtract(const Safe_Poly& q, Degree d, C c)
{
  if (q.isZero() or c==C(0))
    return; // do nothing

  const auto qs = q.size();   // save the original size of q, it might change

  assert(qs+d<=base::size()); // if not, leading coefficient would get negative

  auto* shift_base = &(*this)[d]; // pointer to lowest coefficient affected
  for (size_t j = qs; j-->0; ) // loop must be backwards if |&q == this|
  {
    C a = q[j]; // need a copy for the multiplication
    polynomials::safe_multiply(a,c);        // this may throw
    polynomials::safe_subtract(shift_base[j],a);
  }

  // set degree
  base::adjustSize();
} // |safeSubtract|

/* Again a simplified version deals with the common case |c==1| */

template<typename C>
void Safe_Poly<C>::safeSubtract(const Safe_Poly& q, Degree d)
{
  if (q.isZero()) // do nothing
    return;

  const auto qs = q.size();   // save the original size of q, it might change

  assert(qs+d<=base::size()); // if not, leading coefficient would get negative

  auto* shift_base = &(*this)[d]; // pointer to lowest coefficient affected
  for (size_t j = qs; j-->0;)
    polynomials::safe_subtract(shift_base[j],q[j]);

  // set degree
  base::adjustSize();
 } // |safeSubtract|

// Divide polynomial by scalar |c|, throwing an error is division is inexact
template<typename C>
void Safe_Poly<C>::safeDivide(C c)
{
  for (C& cur_coef : *this)
    polynomials::safe_divide(cur_coef,c);
} // |safeDivide|

/*
  Divide polynomial by $q+1$, imagining if necessary an additional leading
  term of degree $(delta+1)/2$ (must be integer) to make division exact.

  The following reasoning is applied in order to make the strongest possible
  |assert|, while looping only over the coefficients initially present. Two
  cases are possible: the quotient is of degree strictly less than the maximal
  allowed $d=(delta-1)/2$, allowed, or it has degree $d$. The former holds if
  and only if the original polynomial, which is $q+1$ times the quotient with
  any term of degree $d+1$ suppressed, is divisible by $q+1$; on the other
  hand in the latter case the original polynomial must have had degree $d$, as
  its the coefficient in degree $d$ is the sum of the positive leading
  coefficient and the non-negative preceeding coefficient of the quotient.

  So after doing the upward division up to |this->degree()|, we test the value
  in what used to be the leading coefficient. If it is zero, the division was
  exact, and we get the quotient simply by calling |adjustSize| which will
  precisely drop this one coefficient; then we assert the (quotient) degree is
  less than $d$ (expressed as |2*degree()+1<delta|). Otherwise the nonzero
  value we tested is actually the leading term of the quotient, which we
  assert to be of degree $d$ exactly. (For this case we may imagine the loop
  going on one more step to kill off the fictive coefficient in degree $d+1$.)
  In particular we do not need, in this final case, to call |adjustSize|.
 */
template<typename C>
void Safe_Poly<C>::safe_quotient_by_1_plus_q(Degree delta)
{
  if (base::isZero()) // this avoids problems with |base::degree()|
    return; // need not and cannot invent nonzero \mu*q^{d+1} here
  for (size_t j = 1; j <= base::degree(); ++j)
    polynomials::safe_subtract((*this)[j],(*this)[j-1]); // does c[j] -= c[j-1]
  if ((*this)[base::degree()]==0) // test coefficient in old leading term
  { // then upward division was exact: polynomial was already multiple of q+1
    base::adjustSize(); // decreases degree by exactly 1
    assert(2*base::degree()+1<delta); // quotient had degree less than $d$
  }
  else // we need to imagine a term $\mu*q^{d+1}$ with nonzero $\mu$
    assert(2*base::degree()+1==delta); // and quotient must have degree $d$

} // |safe_quotient_by_1_plus_q|


} // |namespace polynomials|

} // |namespace atlas|
