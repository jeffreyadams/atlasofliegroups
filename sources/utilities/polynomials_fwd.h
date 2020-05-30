/*
  This is polynomials_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef POLYNOMIALS_FWD_H  /* guard against multiple inclusions */
#define POLYNOMIALS_FWD_H

namespace atlas {

namespace polynomials {

template<typename C> class Polynomial; // |C| can be any arithmetic type
template<typename C> class Safe_Poly;  // here |C| must be unsigned integral

// template<typename C> class LaurentPolynomial;

// |Degree| does not need to be short, but using more than 32 bits would be silly
typedef unsigned int Degree; // exponent range; not stored

}

}

#endif
