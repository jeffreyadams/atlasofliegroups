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

// walls for alcove containing |gamma|; |integrals| set to those |gamma| lies on
RootNbrSet wall_set
  (const RootDatum& rd, const RatWeight& gamma, RootNbrSet& integrals);

StandardRepr alcove_center(const Rep_context& rc, const StandardRepr& sr);

// try to change |sr| making |N*gamma| integral weight; report whether changed
bool make_multiple_integral
  (const Rep_context& rc, StandardRepr& sr, long long N);

// repeat the above for increasing |N| until |sr*N| has full integral rank
arithmetic::Numer_t  // returns |N|
  simplify(const Rep_context& rc, StandardRepr& sr);

} // |namespace repr|

} // |namespace atlas|

#endif
