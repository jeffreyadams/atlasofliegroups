/*!
\file
  This is arithmetic.cpp.
  This module contains straightforward implementations of some elementary
  arithmetic operations, used in dealing with small abelian groups.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "arithmetic.h"

#include <stdexcept>

#include "constants.h"
#include "intutils.h"
#include "bits.h"
#include "error.h"

/*****************************************************************************

  This module contains straightforward implementations of some elementary
  arithmetic operations, used in dealing with small abelian groups.

******************************************************************************/

namespace atlas {

namespace arithmetic {


/*!
  The classical Euclidian algorithm. It is assumed that

  Precondition: b > 0;
*/
unsigned long unsigned_gcd(unsigned long a, unsigned long b)
{
  unsigned long r;

  while ((r= a%b) != 0)
  {
    a = b;
    b = r;
  }

  return b;
}


unsigned long dummy_gcd, dummy_Bezout; // for default arguments to |lcm|

/*!
  \brief Returns the least common multiple of |a| and |b|, storing their
   greatest common divisor in |gcd| and making |Bezout%a==0|, |Bezout%b==gcd|

  Precondition: a > 0;
*/
unsigned long lcm (unsigned long a, unsigned long b,
		   unsigned long& gcd, unsigned long& Bezout)
{
  unsigned long c=a, d=b; // local variables for the sake of readability
  unsigned long m_c=a, m_d=0; // multiples of |a|, cong. to  |c|,|-d| mod |b|

  do // invariants: |c*m_d+d*m_c==ab|; $m_c\cong c$, $m_d\cong-d$ modulo $b$
  {
    m_d += (d/c)*m_c; d%=c;
    if (d==0)
    { gcd = c; Bezout=m_c; return m_d; }
    m_c += (c/d)*m_d; c%=d;
  } while (c>0);

  gcd = d; Bezout= m_d==0 ? 0 : m_c-m_d; // |m_d<m_c| as we just had |c/d>=1|
  return m_c;
}


/*!
  Synopsis: a *= b mod n.

  Precondition: a < n; b < n.

  NOTE: preliminary implementation. It assumes that |n <= 2^^(longBits/2)|,
  i.e., mudular numbers fit in a half-long
  Exit brutally if this is not fulfilled.
*/
unsigned long& modProd(unsigned long& a, unsigned long b, unsigned long n)
{
  if (n > (1UL << (constants::longBits >> 1)))
    error::FatalError() ("error: overflow in modProd");

  a *= b;
  a %= n;

  return a;
}

unsigned long power(unsigned long x, unsigned int n)
{ if (n==0)
    return 1;
  if (x<=1)
    return x;
  unsigned int l=bits::lastBit(n); // now $n<2^l$
  unsigned long result=1;
  while (l-->0)
  { result *= result;
    if ((n>>l & 1)!=0)
      result*=x;
  }
  return result;
}

Rational Rational::operator+(Rational q) const
{
  unsigned long sum_denom=lcm(denom,q.denom);
  return Rational(num*(sum_denom/denom)+q.num*(sum_denom/q.denom),sum_denom);
}
Rational Rational::operator-(Rational q) const
{
  unsigned long sum_denom=lcm(denom,q.denom);
  return Rational(num*(sum_denom/denom)-q.num*(sum_denom/q.denom),sum_denom);
}

Rational Rational::operator*(Rational q) const
{
  return Rational(num*q.num,denom*q.denom).normalize();
}

Rational Rational::operator/(Rational q) const
{
  if (q.num==0)
    throw std::domain_error("Rational division by 0");
  return Rational(num*q.denom,denom*q.num).normalize();
}

Rational& Rational::power(int n)
{
  normalize();
  unsigned long numer=intutils::abs(num);
  if (n<0)
  { if (num==0)
      throw std::domain_error("Negative power of rational zero");
    std::swap(numer,denom); n=-n;
  }
  numer = arithmetic::power(numer,n); denom = arithmetic::power(denom,n);
  num = (num>0 ? numer : - long(numer));
  return *this;
}

std::ostream& operator<< (std::ostream& out, const Rational& frac)
{
  return out << frac.numerator() << '/' << frac.denominator();
}

} // |namespace arithmetic|

} // |namespace atlas|
