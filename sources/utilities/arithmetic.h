/*!
\file
  This is arithmetic.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef ARITHMETIC_H  /* guard against multiple inclusions */
#define ARITHMETIC_H

/******** function declarations **********************************************/

namespace atlas {

namespace arithmetic {

  unsigned long gcd (long, unsigned long); // signed first argument only!

  unsigned long unsigned_gcd (unsigned long, unsigned long); // avoid overload

  unsigned long lcm (unsigned long, unsigned long);

  // inlined; first two arguments are supposed already reduced modulo third
  unsigned long& modAdd(unsigned long&, unsigned long, unsigned long);

  unsigned long& modProd(unsigned long&, unsigned long, unsigned long);

  template<typename C> unsigned long remainder(C, unsigned long);

} // |namespace arithmetic|


/******** inline function definitions ***************************************/

namespace arithmetic {

  inline unsigned long gcd (long a, unsigned long b) {
    if (a < 0)
      return gcd(static_cast<unsigned long>(-a),b);
    else
      return gcd(static_cast<unsigned long>(a),b);
  }

  inline unsigned long& modAdd(unsigned long& a, unsigned long b,
			       unsigned long n)
    {
      if (a < n-b)
	a += b;
      else
	a -= n-b;
      return a;
    }

} // |namespace arithmetic|

} // |namespace atlas|

#include "arithmetic_def.h"

#endif
