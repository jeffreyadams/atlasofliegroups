/*!
\file
  This is size_def.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

/*****************************************************************************

  ... explain here when it is stable ...

******************************************************************************/

#include "bitset.h"
#include "error.h"

/*****************************************************************************

        Chapter I -- The SizeType class

******************************************************************************/

namespace atlas {

namespace size {


/*!
  Synopsis: constructs the SizeType object representing a.

  Triggers a fatal error if a is not representable as a SizeType
*/
template<typename C> SizeType<C>::SizeType(unsigned long a)
{
  if (a == 0)
    error::PrimesError()("error: 0 is not representable as a SizeType");

  memset(d_data,0,PRIMES_MAX);

  for (size_t n = 0; n < PRIMES_MAX; ++n)
    while (a%prime(n) == 0) {
      ++d_data[n];
      a /= prime(n);
    }

  if (a > 1) // a was not representable
    error::PrimesError()("error: number not representable as SizeType");
}

/******** accessors **********************************************************/


/*!
  Synopsis: return the unsigned long value of prime(j)^^d_data[j]

  The algorithm is the classical algorithm with squarings and multiplications,
  logarithmic in j.

  NOTE: in case of overflow, we simply return the value modulo 2^^longBits.
*/
template<typename C>
  unsigned long SizeType<C>::piece(size_t j) const
{
  unsigned long c = 1;
  unsigned long p = prime(j);

  bitset::BitSet<constants::longBits> r(d_data[j]);

  for (size_t j = r.lastBit(); j;) {
    --j;
    c *= c;
    if (r.test(j))
      c *= p;
  }

  return c;
}


/*!
  Synopsis: returns the unsigned long value of *this.

  NOTE: in case of overflow, we simply return the value modulo 2^^longBits.
*/
template<typename C>
  unsigned long SizeType<C>::toUlong() const
{
  unsigned long c = 1;

  for (size_t j = 0; j < PRIMES_MAX; ++j)
    c *= piece(j);

  return c;
}

/******** manipulators ******************************************************/

/*!
  NOTE: multiplication becomes addition in our logarithmic representation.
*/
template<typename C>
  SizeType<C>& SizeType<C>::operator*= (const SizeType& a)
{
  for (size_t j = 0; j < PRIMES_MAX; ++j) {
    C c = a.d_data[j];
    d_data[j] += c;
  }

  return *this;
}


/*!
  Synopsis: current *= SizeType(a).

  Precondition: a is representable as a SizeType.
*/
template<typename C>
  SizeType<C>& SizeType<C>::operator*= (unsigned long a)
{
  for (size_t n = 0; n < PRIMES_MAX; ++n)
    while (a%prime(n) == 0) {
      ++d_data[n];
      a /= prime(n);
    }

  return *this;
}


/*!
  NOTE: division becomes subtraction in our logarithmic representation.
*/
template<typename C>
  SizeType<C>& SizeType<C>::operator/= (const SizeType& a)
{
  for (size_t j = 0; j < PRIMES_MAX; ++j)
    d_data[j] -= a.d_data[j];

  return *this;
}

} // namespace size

/*****************************************************************************

        Chapter I -- The SizeType class

******************************************************************************/

namespace size {


/*!
  Synopsis: puts in a the factorial of n.

  Precondition: n is not greater than the PRIMES_MAX-th prime number.
*/
template<typename C> void factorial (SizeType<C>& a, unsigned long n)
{
  a.reset();

  for (size_t j = 2; j <= n; ++j)
    a *= j;
}

} // namespace size

} // namespace atlas
