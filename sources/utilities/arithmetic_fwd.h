/*
  This is arithmetic_fwd.h

  Copyright (C) 2011 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef ARITHMETIC_FWD_H  /* guard against multiple inclusions */
#define ARITHMETIC_FWD_H

#include <vector>

/******** forward type declarations ******************************************/

namespace atlas {

namespace arithmetic {

  using Numer_t = long long int;
  using Denom_t = unsigned long long int;
  template<typename I> class Rational;
  using RatNum = Rational<Numer_t>;
  using RatNumList = std::vector<RatNum>;
  class Split_integer;
  class big_int; // defined in bigint.h but fiorward declared here
  class big_rat; // defined in bigint.h but fiorward declared here

} // |namespace arithmetic|

} // |namespace atlas|

#endif
