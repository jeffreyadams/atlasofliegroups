/*
  This is common_blocks.h

  Copyright (C) 2019,2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Variant of |param_block| in blocks.h, to be shared when having same KL data

#ifndef COMMON_BLOCKS_H  /* guard against multiple inclusions */
#define COMMON_BLOCKS_H

#include "blocks.h" // we conceptually just extend that module
#include "subsystem.h"

namespace atlas {


namespace repr {
  class Rep_context;
  class StandardReprMod;
}

namespace ext_block {
  class ext_block;
}

namespace blocks {

// a class for blocks of (possibly non integral) parameters
class block_minimal : public Block_base
{
  const Rep_context& rc; // accesses many things, including KGB set for x

  const RatWeight gamma_mod_1;
  const SubSystem integral_datum;

  std::vector<TorusPart> y_bits; // as in |StandardRepr|, indexed by |y|

  y_entry::Pooltype y_pool;
  y_part_hash y_hash;  // hash table allows storing |y| parts by index

  // hash structure to allow rapid lookup of |(x,y)| index pairs
  block_hash xy_hash;

  std::unique_ptr<ext_block::ext_block> extended;

  // group small data members together:
  KGBElt highest_x,highest_y; // maxima over this block

 public:

  // constructor and destructor
  block_minimal
    (const repr::Rep_context& rc,
     const repr::StandardReprMod& srm, // not modified, no "making dominant"
     BlockElt& entry_element	// set to block element matching input
    );
  ~block_minimal(); // cleans up |*extended|, so inline definition impossible

  // accessors that get values via |rc|
  const repr::Rep_context& context() const { return rc; }
  const RootDatum& rootDatum() const;
  const SubSystem& integral_subsystem() const { return integral_datum; }
  const InnerClass& innerClass() const;
  const InvolutionTable& involution_table() const;
  RealReductiveGroup& realGroup() const;

  // with |gamma| unknown, only the difference |gamma-lambda| is meaningful
  RatWeight gamma_lambda(BlockElt z) const;

  BlockElt lookup(const repr::StandardReprMod& srm) const;
  BlockElt lookup(KGBElt x, const RatWeight& gamma_lambda) const;

  ext_gens fold_orbits(const WeightInvolution& delta) const;

  ext_block::ext_block& extended_block(const WeightInvolution& delta);

  // virtual methods
  virtual KGBElt max_x() const { return highest_x; } // might not be final |x|
  virtual KGBElt max_y() const { return highest_y; }

  virtual std::ostream& print // defined in block_io.cpp
    (std::ostream& strm, BlockElt z,bool as_invol_expr) const;


 private:
  void compute_y_bits();

/*
  reverse lengths and order block with them increasing, and by increasing
  |x(z)| among elements of given length; adapt tables accordingly.
*/
  void reverse_length_and_sort();

}; // |class block_minimal|

} // |namespace blocks|

namespace ext_block
{

class paramin_context
{
  const repr::Rep_context& d_rc;
  const WeightInvolution d_delta;

public:
  paramin_context (const repr::Rep_context& rc, const WeightInvolution& delta)
    : d_rc(rc), d_delta(delta) {}

  // accessors
  const repr::Rep_context& rc () const { return d_rc; }
  const RootDatum& root_datum () const;
  const InnerClass& inner_class () const;
  RealReductiveGroup& real_group () const;
  const WeightInvolution& delta () const { return d_delta; }
  const RatCoweight& g_rho_check() const;
  RatCoweight g() const;

}; // |paramin_context|

// this class is to |paramin| what |ext_block::context| is to |ext_block::param|
// it holds relevant values that remain fixed across extended block
// data fields that are removed: |d_gamma|, |lambda_shifts|
// methods that are absent: |gamma|, |lambda_shift|
class context_minimal : public paramin_context
{
  const RootDatum integr_datum; // intgrality datum
  const SubSystem sub; // embeds |integr_datum| into parent root datum
  Permutation pi_delta; // permutation of |delta| on roots of full root datum
  RootNbrSet delta_fixed_roots;
  weyl::Twist twist;
  int_Vector l_shifts; // of size |sub.rank()|; affine center for action on |l|

 public:
  context_minimal (const repr::Rep_context& rc, const WeightInvolution& delta,
		   const SubSystem& integral_subsystem);

  // accessors
  const RootDatum& id() const { return integr_datum; }
  const SubSystem& subsys() const { return sub; }
  RootNbr delta_of(RootNbr alpha) const { return pi_delta[alpha]; }
  const RootNbrSet& delta_fixed() const { return delta_fixed_roots; }
  weyl::Generator twisted(weyl::Generator s) const { return twist[s]; }
  int l_shift(weyl::Generator s) const { return l_shifts[s]; }

  // whether positive $\alpha$ has $\theta(\alpha)\neq\pm(1|\delta)(\alpha)$
  bool is_very_complex (InvolutionNbr theta, RootNbr alpha) const;
  Weight to_simple_shift(InvolutionNbr theta, InvolutionNbr theta_p,
			 RootNbrSet pos_to_neg) const; // |pos_to_neg| by value
  bool shift_flip(InvolutionNbr theta, InvolutionNbr theta_p,
		  RootNbrSet pos_to_neg) const; // |pos_to_neg| is by value

}; // |context_minimal|

// A variant of |ext_block::param| that avoids fixing |gamma|
// Identical (supplementary) data fields, method absent: |restrict|
struct paramin // allow public member access; methods ensure no invariants
{
  const paramin_context& ctxt;
  TwistedInvolution tw; // implicitly defines $\theta$

  Coweight l; // with |tw| gives a |GlobalTitsElement|; lifts its |t|
  RatWeight gamma_lambda; // lift of $\gamma-\lambda$ value in a |StandardRepr|
  Weight tau; // a solution to $(1-\theta)*\tau=(1-\delta)gamma_\lambda$
  Coweight t; // a solution to $t(1-theta)=l(\delta-1)$
  bool flipped; // whether tensored with the flipping representation

  paramin (const paramin_context& ec, const TwistedInvolution& tw,
	   RatWeight gamma_lambda, Weight tau, Coweight l, Coweight t,
	   bool flipped=false);

  // default extension choice:
  paramin (const paramin_context& ec,
	   KGBElt x, const RatWeight& gamma_lambda, bool flipped=false);
  static paramin default_extend
   (const paramin_context& ec, const repr::StandardRepr& sr);

  paramin (const paramin& p) = default;
  paramin (paramin&& p)
  : ctxt(p.ctxt), tw(std::move(p.tw))
  , l(std::move(p.l))
  , gamma_lambda(std::move(p.gamma_lambda))
  , tau(std::move(p.tau))
  , t(std::move(p.t))
  , flipped(p.flipped)
  {}

  paramin& operator= (const paramin& p);
  paramin& operator= (paramin&& p);

  bool is_flipped() const { return flipped; }

  void flip (bool whether=true) { flipped=(whether!=flipped); }

  const repr::Rep_context& rc() const { return ctxt.rc(); }
  const WeightInvolution& delta () const { return ctxt.delta(); }
  const WeightInvolution& theta () const;

  KGBElt x() const; // reconstruct |x| component
  // underlying unextended representation
  repr::StandardRepr restrict(const RatWeight& gamma) const;
}; // |paramin|

} // |namespace ext_block|

} // |namespace atlas|
#endif
