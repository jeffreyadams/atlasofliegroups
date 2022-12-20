// This is arithmetic.cpp.
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/
/*
  This module contains implementations of some elementary arithmetic operations.
*/

#include "arithmetic.h"

#include <stdexcept> // some cases throw
#include <cassert>   // but most just assert; user should test uncertain cases
#include <limits>
#include <cstdlib>

#include "constants.h"
#include "arithmetic.h"
#include "bits.h"
#include "error.h"

/*****************************************************************************

  This module contains straightforward implementations of some elementary
  arithmetic operations, used in dealing with small abelian groups.

******************************************************************************/

namespace atlas {

namespace arithmetic {


/*
  The classical Euclidean algorithm for positive (indeed unsigned) numbers.
  It is assumed that |b != 0|, but |a| might be zero.
*/
Denom_t unsigned_gcd(Denom_t a, Denom_t b)
{
  do // a double-exit loop
    if ((a %= b)==0)
      return b;
  while ((b %= a)!=0);
  return a;
}


Denom_t dummy_gcd, dummy_mult; // for default arguments to |lcm|

/*
  |lcm(a,b,gcd,mult_a)| returns the least common multiple of |a| and |b|,
  sets |gcd| to their greatest common divisor and |mult_a| to a multiple of
  |a| with |gcd==mult_a%b|. A Bezout coefficient for |gcd| is |mult_a/a|.

  Precondition: a > 0;
*/
Denom_t lcm (Denom_t a, Denom_t b, Denom_t& gcd, Denom_t& mult_a)
{
  Denom_t c=a, d=b; // local variables for the sake of readability
  Denom_t m_c=a, m_d=0; // multiples of |a|, cong. to  |c|,|-d| mod |b|

  do // invariants: |c*m_d+d*m_c==ab|; $m_c\cong c$, $m_d\cong-d$ modulo $b$
  {
    m_d += (d/c)*m_c; d%=c;
    if (d==0)
    { gcd = c; mult_a=m_c; return m_d; }
    m_c += (c/d)*m_d; c%=d;
  } while (c>0);

  gcd = d; mult_a= m_d==0 ? 0 : m_c-m_d; // |m_d<m_c| as we just had |c/d>=1|
  return m_c;
}


/*
  Return a * b mod n.

  Precondition: a < n; b < n.

  NOTE: preliminary implementation. It assumes that |n <= 2^^(longBits/2)|,
  i.e., modular numbers fit in a half-long
  Exit brutally if this is not fulfilled.
*/
Denom_t modProd(Denom_t a, Denom_t b, Denom_t n)
{
  if (n > (1UL << (constants::longBits >> 1)))
    error::FatalError() ("error: overflow in modProd");

  return (a * b) % n;
}

Denom_t power(Denom_t x, unsigned int n)
{ if (n==0)
    return 1;
  if (x<=1)
    return x;
  unsigned int l=bits::lastBit(n); // now $n<2^l$
  Denom_t result=1;
  while (l-->0)
  { result *= result;
    if ((n>>l & 1)!=0)
      result*=x;
  }
  return result;
}

template<typename I>
  Rational<I>& Rational<I>::operator*=(I n)
{ if (n==0)
  { num=0; denom=1; }
  else
  {
    normalize();
    if (n<0)
    { n=-n; num=-num; }
    auto d = unsigned_gcd(denom,n); // converts |n| to unsigned
    denom/=d; // unsigned division (of positive quantities)
    I f = static_cast<I>(n/d); // positive |n| implicitly converted to unsigned
    assert(num<std::numeric_limits<I>::max()/f);
    num *= f;
  }
  return *this; // result will be normalised if |this| was
}

template<typename I>
  Rational<I>& Rational<I>::operator/=(I n)
{ assert(n!=0);
  normalize();
  if (n<0)
  { n=-n; num=-num; }
  I d = static_cast<I>(unsigned_gcd(std::abs(num),n));
  num/=d; // signed divide
  I f = n/d; // signed divide
#ifndef NDEBUG
  using UI = typename std::make_unsigned<I>::type;
  assert(static_cast<UI>(num) < std::numeric_limits<UI>::max()/f);
#endif
  denom *= f;
  return *this; // result will be normalised if |this| was
}

template<typename I>
  Rational<I>& Rational<I>::operator%=(I n)
{ assert(n!=0);
  normalize();
  num = remainder(num,static_cast<I>(denom*std::abs(n)));
  return *this;
}


template<typename I>
  Rational<I> Rational<I>::operator+(Rational<I> q) const
{
  normalize();
  Denom_t sum_denom=lcm(denom,q.denom);
  return Rational<I>(   num*static_cast<I>(sum_denom/denom)
		     +q.num*static_cast<I>(sum_denom/q.denom)
		    , sum_denom);
}
template<typename I>
  Rational<I> Rational<I>::operator-(Rational<I> q) const
{
  normalize();
  Denom_t sum_denom=lcm(denom,q.denom);
  return Rational<I>(  num*static_cast<I>(sum_denom/denom)
		    -q.num*static_cast<I>(sum_denom/q.denom)
		    , sum_denom);
}

template<typename I>
  Rational<I> Rational<I>::operator*(Rational<I> q) const
{
  normalize(); q.normalize();
  Denom_t d = gcd(num,q.denom);
  I new_num = d==1 ? num : (q.denom/=d, num/d);
  auto new_denom = (d = gcd(q.num,denom))==1 ? denom : (q.num/=d, denom/d);
  return Rational<I>(new_num*q.num,new_denom*q.denom);
}

template<typename I>
  Rational<I> Rational<I>::operator/(Rational<I> q) const
{
  assert(q.num!=0);
  normalize(); q.normalize();
  Denom_t d = unsigned_gcd(denom,q.denom);
  auto new_denom = d==1 ? denom : (q.denom/=d, denom/d);
  I sd = gcd(num,std::abs(q.num)); // convert result to signed
  I new_num = sd==1 ? num : (q.num/=sd, num/sd);
  if (q.num>0)
    return Rational<I>(new_num*static_cast<I>(q.denom),new_denom*q.num);
  else
    return Rational<I>(-new_num*static_cast<I>(q.denom),new_denom*-q.num);
}

template<typename I>
  Rational<I> Rational<I>::operator%(Rational<I> q) const
{
  assert(q.num!=0);
  Denom_t d = unsigned_gcd(denom,q.denom);
  if (d==1)
    return Rational<I>
      (remainder
        (num*static_cast<I>(q.denom),static_cast<I>(denom*std::abs(q.num)))
      , denom*q.denom).normalize();
  else
  {
    q.denom/=d;
    return Rational<I>
      (remainder(num*static_cast<I>(q.denom),
		 static_cast<I>((denom/d)*std::abs(q.num)))
      , denom*q.denom).normalize();
  }
}

template<typename I>
  Rational<I>& Rational<I>::power(int n)
{
  normalize();
  Denom_t numer=std::abs(num);
  if (n<0)
  { if (num==0)
      throw std::runtime_error("Negative power of rational zero");
    std::swap(numer,denom); n=-n;
  }
  numer = arithmetic::power(numer,n); denom = arithmetic::power(denom,n);
  num = (num>=0 or n%2==0 ? numer : - static_cast<I>(numer));
  return *this;
}

template<typename I>
  std::ostream& operator<< (std::ostream& out, const Rational<I>& frac)
{
  return out << frac.numerator() << '/' << frac.true_denominator();
}


// Instantiation of templates (only these are generated)

template class Rational<Numer_t>; // the main instance used



} // |namespace arithmetic|

} // |namespace atlas|
