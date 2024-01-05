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

namespace weyl {

RatNum frac_eval(const RootDatum& rd, RootNbr i, const RatWeight& gamma);

// walls for alcove containing |gamma|; |integrals| set to those |gamma| lies on
RootNbrSet wall_set
  (const RootDatum& rd, const RatWeight& gamma, RootNbrSet& integrals);

sl_list<RootNbr> sorted_by_label
  (const RootSystem& rs, const RootNbrSet& walls);

// replace |walls| by the walls of the fundamental alcove in some permutation,
// and return Weyl word for element that maps the fundamental alcove to one
// with walls in the same configuration (up to translation) as orignal |walls|
WeylWord from_fundamental_alcove (const RootSystem& rs, RootNbrSet& walls);

StandardRepr alcove_center (const Rep_context& rc, const StandardRepr& sr);

Weight root_vertex_of_alcove (const RootDatum& rd, const RatWeight& gamma);

// try to change |sr| making |N*gamma| integral weight; report whether changed
bool make_multiple_integral
  (const Rep_context& rc, StandardRepr& sr, long long N);

// repeat the above for increasing |N| until |sr*N| has full integral rank
arithmetic::Numer_t  // returns |N|
  simplify(const Rep_context& rc, StandardRepr& sr);


sl_list<WeylElt> finite_subquotient
  (const RootSystem& rs, const WeylGroup& W, RootNbrSet stab, RootNbr alpha);
sl_list<WeylElt> complete_affine_component
  (const RootSystem& rs, const WeylGroup& W, RootNbrSet stab, RootNbr alpha);

sl_list<WeylElt> basic_orbit_ws
(const int_Matrix& Cartan, Generator i,
 const WeylGroup& W, std::vector<WeylElt> gens);
sl_list<WeylElt> affine_orbit_ws
  (const RootDatum& rd, const WeylGroup& W, const RatWeight& gamma);

sl_list<std::pair<WeylElt,sl_list<int_Vector> > > FPP_w_shifts
  (const RootDatum& rd, const WeylGroup& W, const RatWeight& gamma);

sl_list<int_Vector> FPP_orbit_numers
  (const RootDatum& rd, const WeylGroup& W, const RatWeight& gamma);

} // |namespace weyl|

} // |namespace atlas|

#endif
