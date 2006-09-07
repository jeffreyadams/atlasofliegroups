/*!
\file
\brief Forward class declarations for Polynomial and LaurentPolynomial.
*/
/*
  This is polynomials_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#ifndef POLYNOMIALS_FWD_H  /* guard against multiple inclusions */
#define POLYNOMIALS_FWD_H

namespace atlas {

namespace polynomials {

template<typename C> class Polynomial;

template<typename C> class LaurentPolynomial;

typedef size_t Degree;

}

}

#endif
