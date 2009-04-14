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


unsigned long dummy; // for default |gcd| argument to |lcm|

/*!
  Synopsis: returns the least common multiple of |a| and |b|,
  while storing their greatest common divisor in |gcd|.

  Precondition: b > 0;
*/
unsigned long lcm (unsigned long a, unsigned long b, unsigned long& gcd)
{
  unsigned long c=a, d=b; // local variables for the sake of readability
  unsigned long m_c=0, m_d=b; // multiples of |b|, cong. to  |-c|,|d| mod |a|

  do // invariant: |c*m_d+d*m_c==ab|
  {
    m_c += (c/d)*m_d; c%=d;
    if (c==0)
    { gcd = d; return m_c; }
    m_d += (d/c)*m_c; d%=c;
  } while (d>0);

  gcd = c; return m_d;
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

Rational Rational::operator+(Rational q) const
{
  unsigned long sum_denom=lcm(denom,q.denom);
  return Rational(num*(sum_denom/denom)+q.num*(sum_denom/q.denom),sum_denom);
}
Rational Rational::operator-(Rational q) const
{
  unsigned long sum_denom=lcm(denom,q.denom);
  return Rational(num*(sum_denom/denom)/q.num*(sum_denom/q.denom),sum_denom);
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

std::ostream& operator<< (std::ostream& out, const Rational& frac)
{
  return out << frac.numerator() << '/' << frac.denominator();
}

} // |namespace arithmetic|

} // |namespace atlas|
