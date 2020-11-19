/*
  This is alcoves.h

  Copyright (C) 2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/
// Functions relating to alcoves, for use in deformation computations

#include "../Atlas.h"


#ifndef ALCOVES_H  /* guard against multiple inclusions */
#define ALCOVES_H

namespace atlas {

namespace repr {

RatNum frac_eval(const RootDatum& rd, RootNbr i, const RatWeight& gamma);

// try to change |sr| making |N*gamma| integral weight; report whether changed
bool make_multiple_integral
  (const Rep_context& rc, StandardRepr& sr, long long N);

} // |namespace repr|

} // |namespace atlas|

#endif
