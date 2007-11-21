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

  // inlined
  unsigned long gcd (unsigned long, unsigned long);

  unsigned long lcm (unsigned long, unsigned long);

  // inlined
  unsigned long& modAdd(unsigned long&, unsigned long, unsigned long);

  unsigned long& modProd(unsigned long&, unsigned long, unsigned long);

  template<typename C> unsigned long remainder(C, unsigned long);

}

/******** inline function definitions ***************************************/

namespace arithmetic {

  inline unsigned long gcd (long a, unsigned long b) {
    if (a < 0) 
      return gcd(static_cast<unsigned long>(-a),b); 
    else 
      return gcd(static_cast<unsigned long>(a),b);
  }

  inline unsigned long& modAdd(unsigned long& a, unsigned long b, 
			      unsigned long n) {
    if (a < n-b)
      a += b;
    else
      a -= n-b;
    return a;
  }

}

}

#include "arithmetic_def.h"

#endif
