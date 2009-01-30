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

  Implementation of |SizeType| class template

******************************************************************************/

#include "bitset.h"
#include "error.h"

/*****************************************************************************

        Chapter I -- The SizeType class template

******************************************************************************/

namespace atlas {

namespace size {


/*!
  Synopsis: constructs the SizeType object representing value |a|.

  Triggers a fatal error if a is not representable as a SizeType
*/
template<typename C> SizeType<C>::SizeType(unsigned long a)
{
  if (a == 0)
    error::PrimesError()("error: 0 is not representable as a SizeType");

  for (size_t n=0; n<PRIMES_MAX; ++n)
    for (d_exp[n]=0; a%prime(n)==0; a/=prime(n))
      ++d_exp[n];


  if (a > 1) // a was not representable
    error::PrimesError()("error: number not representable as SizeType");
}

/******** accessors **********************************************************/


/*!
  Synopsis: return the unsigned long value of prime(i)^^d_exp[i]

  The algorithm is the classical algorithm with squarings and multiplications,
  logarithmic in i.

  NOTE: in case of overflow, we simply return the value modulo 2^^longBits.
*/
template<typename C>
  unsigned long SizeType<C>::piece(size_t k) const
{
  unsigned long c = 1;
  unsigned long p = prime(k);

  bitset::BitSet<constants::longBits> r(d_exp[k]);

  for (size_t i=r.lastBit(); i-->0;) // |lastBit()| is one beyond last bit
  {
    c *= c;
    if (r.test(i)) // always true on first iteration, if any
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

  for (size_t i=0; i<PRIMES_MAX; ++i)
    c *= piece(i);

  return c;
}

/******** manipulators ******************************************************/

/*!
  NOTE: multiplication becomes addition in our logarithmic representation.
*/
template<typename C>
  SizeType<C>& SizeType<C>::operator*= (const SizeType& a)
{
  for (size_t i=0; i < PRIMES_MAX; ++i)
    d_exp[i] += a.d_exp[i];

  return *this;
}


/*!
  NOTE: division becomes subtraction in our logarithmic representation.
*/
template<typename C>
  SizeType<C>& SizeType<C>::operator/= (const SizeType& a)
{
  for (size_t i=0; i < PRIMES_MAX; ++i)
    d_exp[i] -= a.d_exp[i];

  return *this;
}

} // namespace size

/*****************************************************************************

        Chapter II -- The |factorial| function

******************************************************************************/

namespace size {



} // namespace size

} // namespace atlas
