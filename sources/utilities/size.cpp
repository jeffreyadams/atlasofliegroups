/*
  This is size.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "size.h"

#include "error.h"

/*****************************************************************************

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace {

  using namespace size;

  class PrimeHelper {

  private:

    unsigned long d_primes[PRIMES_MAX];

  public:

    // constructors and destructors
    explicit PrimeHelper(unsigned long);
    ~PrimeHelper() {}

    // accessors
    unsigned long operator[] (size_t j) {
      return d_primes[j];
    }
  };

}

/*****************************************************************************

        Chapter I -- Functions declared in size.h

******************************************************************************/

namespace size {

unsigned long prime(size_t j)

/*
  Synopsis: returns the j-th prime number, counting from 0.

  Precondition: j < PRIMES_MAX;
*/

{
  static PrimeHelper p(constants::RANK_MAX+1);

  return p[j];
}

}

/*****************************************************************************

        Chapter II -- The PrimeHelper class.

******************************************************************************/

namespace {

PrimeHelper::PrimeHelper(unsigned long n)

/*
  Synopsis: constructs a PrimeHelper object containing the primes <= n.

  Also checks that PRIMES_MAX is large enough to contain the list; triggers
  a fatal error otherwise.
*/

{
  using namespace error;

  size_t size = 0;

  for (size_t j = 2; j <= n; ++j) {
    for (size_t i = 0; i < size; ++i) {
      if (j%d_primes[i] == 0) // j is not prime
	goto nextj;
      if (j/d_primes[i] < d_primes[i]) // no further divisions necessary
	break;
    }
    if (size == PRIMES_MAX) // error
      PrimesError()("error: value of PRIMES_MAX is too small");
    d_primes[size] = j;
    ++size;
  nextj:
    continue;
  }
}

}

}
