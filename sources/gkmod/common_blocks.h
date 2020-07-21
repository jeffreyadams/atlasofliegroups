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


namespace repr {

/*
  This class is to |paramin| what |ext_block::context| is to |ext_block::param|
  it holds relevant values that remain fixed across extended block
  data fields that are removed: |d_gamma|, |lambda_shifts|
  methods that are absent: |gamma|, |lambda_shift|
*/
class Ext_common_context : public repr::common_context
{
  const WeightInvolution d_delta;
  Permutation pi_delta; // permutation of |delta| on roots of full root datum
  RootNbrSet delta_fixed_roots;
  weyl::Twist twist;
  int_Vector l_shifts; // of size |sub.rank()|; affine center for action on |l|

 public:
  Ext_common_context (RealReductiveGroup& G, const WeightInvolution& delta,
		      const SubSystem& integral_subsystem);

  // accessors
  const WeightInvolution& delta () const { return d_delta; }
  RootNbr delta_of(RootNbr alpha) const { return pi_delta[alpha]; }
  const RootNbrSet& delta_fixed() const { return delta_fixed_roots; }
  weyl::Generator twisted(weyl::Generator s) const { return twist[s]; }
  int l_shift(weyl::Generator s) const { return l_shifts[s]; }

  // whether positive $\alpha$ has $\theta(\alpha)\neq\pm(1|\delta)(\alpha)$
  bool is_very_complex(InvolutionNbr theta, RootNbr alpha) const;
  bool shift_flip(InvolutionNbr theta, InvolutionNbr theta_p,
		  RootNbrSet pos_to_neg) const; // |pos_to_neg| is by value

}; // |Ext_common_context|

} // |namespace repr|

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
