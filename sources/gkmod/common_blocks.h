/*
  This is common_blocks.h

  Copyright (C) 2019,2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Variant of |param_block| in blocks.h, to be shared when having same KL data

#ifndef COMMON_BLOCKS_H  /* guard against multiple inclusions */
#define COMMON_BLOCKS_H

#include "../Atlas.h"
#include "blocks.h" // we conceptually just extend that module
#include "gradings.h" // type |gradings::Status| must be complete
#include "bruhat.h" // type |BruhatOrder| must be complete for destructor
#include "subsystem.h"
#include "repr.h"
#include "kl.h" // for type |Poly_hash_export|

namespace atlas {

namespace ext_block {

// A variant of |ext_block::param| that avoids fixing |gamma|
// Identical (supplementary) data fields, method absent: |restrict|
struct ext_param // allow public member access; methods ensure no invariants
{
  const repr::Ext_rep_context& ctxt;
  TwistedInvolution tw; // implicitly defines $\theta$

  Coweight l; // with |tw| gives a |GlobalTitsElement|; lifts its |t|
  RatWeight gamma_lambda; // lift of $\gamma-\lambda$ value in a |StandardRepr|
  Weight tau; // a solution to $(1-\theta)*\tau=(1-\delta)gamma_\lambda$
  Coweight t; // a solution to $t(1-theta)=l(\delta-1)$
  bool flipped; // whether tensored with the flipping representation

  ext_param (const repr::Ext_rep_context& ec, const TwistedInvolution& tw,
	   RatWeight gamma_lambda, Weight tau, Coweight l, Coweight t,
	   bool flipped=false);

  // default extension choice:
  ext_param (const repr::Ext_rep_context& ec,
	     KGBElt x, const RatWeight& gamma_lambda, bool flipped=false);
  static ext_param default_extend
  (const repr::Ext_rep_context& ec, const repr::StandardRepr& sr);

  ext_param (const ext_param& p) = default;
  ext_param (ext_param&& p)
  : ctxt(p.ctxt), tw(std::move(p.tw))
  , l(std::move(p.l))
  , gamma_lambda(std::move(p.gamma_lambda))
  , tau(std::move(p.tau))
  , t(std::move(p.t))
  , flipped(p.flipped)
  {}

  ext_param& operator= (const ext_param& p);
  ext_param& operator= (ext_param&& p);

  bool is_flipped() const { return flipped; }

  void flip (bool whether=true) { flipped=(whether!=flipped); }

  const repr::Rep_context& rc() const { return ctxt; } // reference base object
  const WeightInvolution& delta () const { return ctxt.delta(); }
  const WeightInvolution& theta () const;

  KGBElt x() const; // reconstruct |x| component
  // underlying unextended representation
  repr::StandardRepr restrict(const RatWeight& gamma) const;
}; // |ext_param|

} // |namespace ext_block|

} // |namespace atlas|
#endif
