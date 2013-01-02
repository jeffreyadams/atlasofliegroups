/*!
\file
  This is arithmetic.cpp.
  This module contains straightforward implementations of some elementary
  arithmetic operations, used in dealing with small abelian groups.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "arithmetic.h"

#include <stdexcept>

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


/*!
  The classical Euclidian algorithm. It is assumed that b > 0;
*/
Denom_t unsigned_gcd(Denom_t a, Denom_t b)
{
  do
    if ((a %= b)==0)
      return b;
  while ((b %= a)!=0);
  return a;
}


Denom_t dummy_gcd, dummy_mult; // for default arguments to |lcm|

/* |lcm(a,b,gcd,mult_a)| returns the least common multiple of |a| and |b|,
   sets |gcd| to their greatest common divisor and |mult_a| to a multiple of
   |a| with |gcd==mult_a%b|. A Bezout coefficient for |a| is |mult_a/gcd|.

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


/*!
  Synopsis: return a * b mod n.

  Precondition: a < n; b < n.

  NOTE: preliminary implementation. It assumes that |n <= 2^^(longBits/2)|,
  i.e., mudular numbers fit in a half-long
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

Rational& Rational::operator+=(Numer_t n) { num+=n*denom; return *this; }
Rational& Rational::operator-=(Numer_t n) { num-=n*denom; return *this; }

Rational& Rational::operator*=(Numer_t n)
{ if (n==0)
  { num=0; denom=1; }
  else
  {
    if (n<0)
    { n=-n; num=-num; }
    if (denom%n==0)
      denom/=n;
    else
      num*=n;
  }
  return *this;
}

Rational& Rational::operator/=(Numer_t n)
{ if (n<0)
  { n=-n; num=-num; }
  if (num%n==0)
    num/=n;
  else
    denom*=n;
  return *this;
}

Rational Rational::operator+(Rational q) const
{
  Denom_t sum_denom=lcm(denom,q.denom);
  return Rational(num*Numer_t(sum_denom/denom)+q.num*Numer_t(sum_denom/q.denom),
		  sum_denom);
}
Rational Rational::operator-(Rational q) const
{
  Denom_t sum_denom=lcm(denom,q.denom);
  return Rational(num*Numer_t(sum_denom/denom)-q.num*Numer_t(sum_denom/q.denom),
		  sum_denom);
}

Rational Rational::operator*(Rational q) const
{
  return Rational(num*q.num,denom*q.denom).normalize();
}

Rational Rational::operator/(Rational q) const
{
  if (q.num>0)
    return Rational(num*Numer_t(q.denom),denom*q.num).normalize();
  else if (q.num==0)
    throw std::domain_error("Rational division by 0");
  else
    return Rational(-num*Numer_t(q.denom),denom*-q.num).normalize();
}

Rational& Rational::power(int n)
{
  normalize();
  Denom_t numer=arithmetic::abs(num);
  if (n<0)
  { if (num==0)
      throw std::domain_error("Negative power of rational zero");
    std::swap(numer,denom); n=-n;
  }
  numer = arithmetic::power(numer,n); denom = arithmetic::power(denom,n);
  num = (num>=0 ? numer : - Numer_t(numer));
  return *this;
}

std::ostream& operator<< (std::ostream& out, const Rational& frac)
{
  return out << frac.numerator() << '/' << frac.denominator();
}

} // |namespace arithmetic|

} // |namespace atlas|
