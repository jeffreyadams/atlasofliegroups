/*!
\file
  This is arithmetic.cpp.
  This module contains straightforward implementations of some elementary
  arithmetic operations, used in dealing with small abelian groups.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
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

unsigned int modular_int::modulus=16; // default modulus

unsigned long gcd(unsigned long a, unsigned long b)

/*!
  Synopsis: the classic Euclidian algorithm. It is assumed that

  Precondition: b > 0;
*/

{
  if (a == 0)
    return b;

  if (a < b)  /* exchange a and b */
    return gcd(b,a);

  unsigned long r = a%b;

  while (r != 0)
    {
      a = b;
      b = r;
      r = a%b;
    }

  return b;
}

unsigned long lcm(unsigned long a, unsigned long b)

/*!
  Synopsis: returns the lowest common multiple of a and b.

  Precondition: b > 0;
*/

{
  unsigned long c = gcd(a,b);
  return (a/c)*b;
}

unsigned long& modProd(unsigned long& a, unsigned long b, unsigned long n)

/*!
  Synopsis: a *= b mod n.

  Precondition: a < n; b < n.

  NOTE: preliminary implementation. It assumes that n <= 2^^(longBits/2).
  Will exit brutally if this is not fulfilled.
*/

{
  using namespace constants;
  using namespace error;

  if (n > (1UL << (longBits >> 1)))
    FatalError() ("error: overflow in modProd");

  a *= b;
  a %= n;

  return a;
}

}

}
