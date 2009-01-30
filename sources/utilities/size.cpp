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

/*****************************************************************************

  Provide a one-time precomputation of the necessary seqence of small primes,
  which are then accessible as |primes(i)|

******************************************************************************/

namespace atlas {

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

        Chapter I -- Functions declared in size.h

******************************************************************************/

namespace size {


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

} // |namespace size|

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

} // |namespace atlas|
