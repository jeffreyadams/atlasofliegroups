/*
  This is size_def.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
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

template<typename C> SizeType<C>::SizeType(unsigned long a)

/*
  Synopsis: constructs the SizeType object representing a.

  Triggers a fatal error if a is not representable as a SizeType
*/

{
  using namespace error;

  if (a == 0)
    PrimesError()("error: 0 is not representable as a SizeType");

  memset(d_data,0,PRIMES_MAX);

  for (size_t n = 0; n < PRIMES_MAX; ++n)
    while (a%prime(n) == 0) {
      ++d_data[n];
      a /= prime(n);
    }

  if (a > 1) // a was not representable
    PrimesError()("error: number not representable as SizeType");
}

/******** accessors **********************************************************/

template<typename C>
  unsigned long SizeType<C>::piece(size_t j) const

/*
  Synopsis: return the unsigned long value of prime(j)^^d_data[j]

  The algorithm is the classical algorithm with squarings and multiplications,
  logarithmic in j.

  NOTE: in case of overflow, we simply return the value modulo 2^^longBits.
*/

{
  using namespace bitset;
  using namespace constants;

  unsigned long c = 1;
  unsigned long p = prime(j);

  BitSet<longBits> r(d_data[j]);

  for (size_t j = r.lastBit(); j;) {
    --j;
    c *= c;
    if (r.test(j))
      c *= p;
  }

  return c;
}

template<typename C>
  unsigned long SizeType<C>::toUlong() const

/*
  Synopsis: returns the unsigned long value of *this.

  NOTE: in case of overflow, we simply return the value modulo 2^^longBits.
*/

{
  unsigned long c = 1;

  for (size_t j = 0; j < PRIMES_MAX; ++j)
    c *= piece(j);

  return c;
}

/******** manipulators ******************************************************/

template<typename C>
  SizeType<C>& SizeType<C>::operator*= (const SizeType& a)

/*
  NOTE: multiplication becomes addition in our logarithmic representation.

  NOTE: we explicitly convert to long to avoid compiler warnings.
*/

{
  for (size_t j = 0; j < PRIMES_MAX; ++j) {
    C c = a.d_data[j];
    d_data[j] += c;
  }

  return *this;
}

template<typename C>
  SizeType<C>& SizeType<C>::operator*= (unsigned long a)

/*
  Synopsis: current *= SizeType(a).

  Precondition: a is representable as a SizeType.
*/

{  
  for (size_t n = 0; n < PRIMES_MAX; ++n)
    while (a%prime(n) == 0) {
      ++d_data[n];
      a /= prime(n);
    }

  return *this;
}

template<typename C>
  SizeType<C>& SizeType<C>::operator/= (const SizeType& a)

/*
  NOTE: division becomes subtraction in our logarithmic representation.
*/

{
  for (size_t j = 0; j < PRIMES_MAX; ++j)
    d_data[j] -= a.d_data[j];

  return *this;
}

}

/*****************************************************************************

        Chapter I -- The SizeType class

******************************************************************************/

namespace size {

template<typename C> void factorial (SizeType<C>& a, unsigned long n)

/*
  Synopsis: puts in a the factorial of n.

  Precondition: n is not greater than the PRIMES_MAX-th prime number.
*/

{
  a.reset();

  for (size_t j = 2; j <= n; ++j)
    a *= j;

  return;
}

}

}
