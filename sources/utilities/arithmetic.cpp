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


/*!
  Synopsis: returns the lowest common multiple of a and b.

  Precondition: b > 0;
*/
unsigned long lcm (unsigned long a, unsigned long b)
{
  unsigned long c=a, d=b; // local variables for the sake of readability
  unsigned long m_c=0, m_d=b; // multiples of |b|, cong. to  |-c|,|d| mod |a|

  do // invariant: |c*m_d+d*m_c==ab|
  {
    m_c += (c/d)*m_d; c%=d;
    if (c==0) return m_c;
    m_d += (d/c)*m_c; d%=c;
  } while (d>0);

  return m_d;
}

unsigned long& modProd(unsigned long& a, unsigned long b, unsigned long n)

/*!
  Synopsis: a *= b mod n.

  Precondition: a < n; b < n.

  NOTE: preliminary implementation. It assumes that |n <= 2^^(longBits/2)|,
  i.e., mudular numbers fit in a half-long
  Exit brutally if this is not fulfilled.
*/

{
  if (n > (1UL << (constants::longBits >> 1)))
    error::FatalError() ("error: overflow in modProd");

  a *= b;
  a %= n;

  return a;
}

}

}
