/*!
\file
  This is size.cpp
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "size.h"
#include "error.h"
#include "bitset.h"

/*****************************************************************************

  Provide a one-time precomputation of the necessary seqence of small primes,
  which are then accessible as |primes(i)|

******************************************************************************/

namespace atlas {

namespace size {

namespace {

// a class instantiated once as a static object; ctor does precomputation
class PrimeHelper
{
  unsigned long d_primes[size::PRIMES_MAX];
public:

  explicit PrimeHelper(unsigned long n); // computes primes up to |n|
  ~PrimeHelper() {}

  unsigned long operator[] (size_t j) { return d_primes[j]; }
}; // |class PrimeHelper|

} // |namespace|


/*****************************************************************************

        Chapter I -- The SizeType class template

******************************************************************************/


/*!
  Synopsis: constructs the SizeType object representing value |a|.

  Triggers a fatal error if a is not representable as a SizeType
*/
template<typename C> SizeType<C>::SizeType(unsigned long a)
{
  if (a == 0)
    error::PrimesError()("error: 0 is not representable as a SizeType");

  for (size_t n=0; n<PRIMES_MAX; ++n) // extract all factors |prime(n)|
    for (d_exp[n]=0; a%prime(n)==0; ++d_exp[n])
      a/=prime(n);


  if (a > 1) // a was not representable
    error::PrimesError()("error: number not representable as SizeType");
}
/*!
  Synopsis: return the unsigned long value of $prime(i)^{d_exp[i]}$

  The algorithm is the classical algorithm with squarings and multiplications,
  logarithmic in i.

  NOTE: in case of overflow, we simply return the value modulo $2^{longBits}$.
*/
template<typename C>
  unsigned long SizeType<C>::piece(size_t k) const
{
  unsigned long c = 1;
  unsigned long p = prime(k);

  bitset::BitSet<constants::longBits> r(d_exp[k]); // get hands on bits

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

  NOTE: in case of overflow, we simply return the value modulo $2^{longBits}$.
*/
template<typename C>
  unsigned long SizeType<C>::toUlong() const
{
  unsigned long c = 1;

  for (size_t i=0; i<PRIMES_MAX; ++i)
    c *= piece(i);

  return c;
}

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

/*****************************************************************************

        Chapter II -- Functions declared in size.h

******************************************************************************/


/*!
  Synopsis: returns the j-th prime number, counting from 0.

  Precondition: j < PRIMES_MAX;
*/
unsigned long prime(size_t j)
{
  static PrimeHelper p(constants::RANK_MAX+1);
  return p[j];
}

//! Return the factorial of |n|.
Size factorial (unsigned long n)
{
  Size result(1);

  for (; n>1; --n)
    result *= Size(n);

  return result;
}


/*****************************************************************************

        Chapter II -- The PrimeHelper class.

******************************************************************************/

namespace {


/*!
  Synopsis: constructs a PrimeHelper object containing the primes <= n.

  Also checks that PRIMES_MAX is large enough to contain the list; triggers
  a (fatal) |error::PrimesError| otherwise.
*/
PrimeHelper::PrimeHelper(unsigned long n)
{
  size_t count=0; // number of primes found

  for (size_t j=2; j <= n; ++j) // |j| is candidate prime
  {
    size_t i; // allow inspection after loop
    for (i=0; i<count; ++i) // test for division by |prime(i)|
    {
      // the next test is optional, only saves time if |n| is fairly large
      if (d_primes[i]*d_primes[i]>j) // no further divisions necessary
      { i=count; break; } // skip useless remainder of loop

      if (j%d_primes[i] == 0) // j is not prime
	break; // exit with |i<count|
    }

    if (i==count) // then |j| is prime
    {
      if (count == size::PRIMES_MAX) // error, cannot add one more prime
	error::PrimesError()("error: value of PRIMES_MAX is too small");
      d_primes[count]=j;
      ++count;
    }
  }
}

} // |namespace|

// instantiate
template class SizeType<BaseType>;

} // |namespace size|

} // |namespace atlas|
