/*
  This is block_minimal.h

  Copyright (C) 2019 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Variant of |param_block| in blocks.h, to be shared when having same KL data

#ifndef BLOCK_MINIMAL_H  /* guard against multiple inclusions */
#define BLOCK_MINIMAL_H

#include "blocks.h" // we conceptually just extend that module
#include "subsystem.h"
#include "repr.h" // for |repr::Rep_context|

namespace atlas {

namespace blocks {

// a class for blocks of (possibly non integral) parameters
class block_minimal : public Block_base
{
  const Rep_context& rc; // accesses many things, including KGB set for x
  const SubSystem integral_datum;

  y_entry::Pooltype y_pool;
  y_part_hash y_hash;  // hash table allows storing |y| parts by index
  std::vector<TorusElement> y_part; // as in the y_values module, indexed by |y|

  // hash structure to allow rapid lookup of |(x,y)| index pairs
  block_hash xy_hash;

  // group small components together:
  KGBElt highest_x,highest_y; // maxima over this block

 public:

  // constructor
  block_minimal
    (const repr::Rep_context& rc,
     StandardRepr sr, // by value,since it will be made dominant before use
     BlockElt& entry_element	// set to block element matching input
    );

 public:
  // accessors that get values via |rc|
  const repr::Rep_context& context() const { return rc; }
  const RootDatum& rootDatum() const;
  const SubSystem& integral_subsystem() const { return integral_datum; }
  const InnerClass& innerClass() const;
  const InvolutionTable& involution_table() const;
  RealReductiveGroup& realGroup() const;

  // with |gamma| unknown, only the difference |gamma-lambda| is meaningful
  RatWeight gamma_lambda(BlockElt z) const;

  BlockElt lookup(const StandardRepr& sr) const;
  BlockElt lookup(KGBElt x, const RatWeight& gamma_lambda) const;

  ext_gens fold_orbits(const WeightInvolution& delta) const;

  // hash code and comparison, to make blocks hashable and rapidly searchable
  bool operator != (const block_minimal& other) const;
  size_t hash_value (KGBElt x, const RatWeight& gamma_lambda) const;
  BlockElt find (const StandardRepr& sr) const; // return index or |UndefBlock|

  // virtual methods
  virtual KGBElt max_x() const { return highest_x; } // might not be final |x|
  virtual KGBElt max_y() const { return highest_y; }

  virtual std::ostream& print // defined in block_io.cpp
    (std::ostream& strm, BlockElt z,bool as_invol_expr) const;


 private:
/*
  reverse lengths and order block with them increasing, and by increasing
  |x(z)| among elements of given length; adapt tables accordingly.
*/
  void reverse_length_and_sort();

}; // |class block_minimal|

size_t hash_value (const repr::Rep_context& rc, const RootNbrSet& ipr,
		   KGBElt x, const RatWeight& gamma_lambda);

// A class that helps look up standard representation agains a list of blocks
class block_hash_table
{
  struct record { const block_minimal* block; BlockElt z; };

  const repr::Rep_context& rc;
  std::vector<record> hash_table;
  size_t count;

  static const record empty_record;
  static const float sparse_factor; // factor by which to oversize |hash_table|

public:

  // static constants
  static constexpr unsigned long empty = ~0ul;

  block_hash_table(const repr::Rep_context& context)
    : rc(context)
    , hash_table(64,empty_record)
    , count(0)
  {}

  unsigned long lookup
    (const RootNbrSet& ipr,KGBElt x, const RatWeight& gamma_lambda) const;

  void add_block(const block_minimal& block);
}; // |block_hash_table|

} // |namespace blocks|

namespace ext_block
{
class context_minimal // holds values that remain fixed across extended block
{
  const repr::Rep_context& d_rc;
  const WeightInvolution d_delta;

  const RootDatum integr_datum; // intgrality datum
  const SubSystem sub; // embeds |integr_datum| into parent root datum
  Permutation pi_delta; // permutation of |delta| on roots of full root datum
  RootNbrSet delta_fixed_roots;
  weyl::Twist twist;
  int_Vector l_shifts; // affine center for complex cross actions on |l|

 public:
  context_minimal (const repr::Rep_context& rc, const WeightInvolution& delta,
		   const SubSystem& integral_subsystem);

  // accessors
  const repr::Rep_context& rc () const { return d_rc; }
  const RootDatum& id() const { return integr_datum; }
  const SubSystem& subsys() const { return sub; }
  const RootDatum& rootDatum() const { return d_rc.rootDatum(); }
  const InnerClass& innerClass () const { return d_rc.innerClass(); }
  RealReductiveGroup& realGroup () const { return d_rc.realGroup(); }
  const WeightInvolution& delta () const { return d_delta; }
  const RatCoweight& g_rho_check() const { return realGroup().g_rho_check(); }
  RatCoweight g() const { return realGroup().g(); }
  RootNbr delta_of(RootNbr alpha) const { return pi_delta[alpha]; }
  const RootNbrSet& delta_fixed() const { return delta_fixed_roots; }
  weyl::Generator twisted(weyl::Generator s) const { return twist[s]; }
  int l_shift(weyl::Generator s) const { return l_shifts[s]; }

  // whether positive $\alpha$ has $\theta(\alpha)\neq\pm(1|\delta)(\alpha)$
  bool is_very_complex (InvolutionNbr theta, RootNbr alpha) const;
  Weight to_simple_shift(InvolutionNbr theta, InvolutionNbr theta_p,
			 RootNbrSet pos_to_neg) const; // |pos_to_neg| is by value
  bool shift_flip(InvolutionNbr theta, InvolutionNbr theta_p,
		  RootNbrSet pos_to_neg) const; // |pos_to_neg| is by value

}; // |context_minimal|

// A variant of |ext_block::param| that avoids fixing |gamma|
struct paramin // allow public member access; methods ensure no invariants
{
  const context_minimal& ctxt;
  TwistedInvolution tw; // implicitly defines $\theta$

  Coweight l; // with |tw| gives a |GlobalTitsElement|; lifts its |t|
  RatWeight gamma_lambda; // lift of $\gamma-\lambda$ value in a |StandardRepr|
  Weight tau; // a solution to $(1-\theta)*\tau=(1-\delta)gamma_\lambda$
  Coweight t; // a solution to $t(1-theta)=l(\delta-1)$
  bool flipped; // whether tensored with the fliiping representation

  paramin (const context_minimal& ec,
	   KGBElt x, const RatWeight& gamma_lambda, bool flipped=false);
  paramin (const context_minimal& ec, const TwistedInvolution& tw,
	   RatWeight gamma_lambda, Weight tau, Coweight l, Coweight t,
	   bool flipped=false);

  paramin (const paramin& p) = default;
  paramin (paramin&& p)
  : ctxt(p.ctxt), tw(std::move(p.tw))
  , l(std::move(p.l))
  , gamma_lambda(std::move(p.gamma_lambda))
  , tau(std::move(p.tau))
  , t(std::move(p.t))
  , flipped(p.flipped)
  {}

  paramin& operator= (const paramin& p)
  { tw=p.tw;
    l=p.l; gamma_lambda=p.gamma_lambda; tau=p.tau; t=p.t;
    flipped=p.flipped;
    return *this;
  }
  paramin& operator= (paramin&& p)
  { tw=std::move(p.tw);
    l=std::move(p.l); gamma_lambda=std::move(p.gamma_lambda);
    tau=std::move(p.tau); t=std::move(p.t);
    flipped=p.flipped;
    return *this;
  }

  bool is_flipped() const { return flipped; }

  void flip (bool whether=true) { flipped=(whether!=flipped); }

  const repr::Rep_context rc() const { return ctxt.rc(); }
  const WeightInvolution& delta () const { return ctxt.delta(); }
  const WeightInvolution& theta () const
    { return ctxt.innerClass().matrix(tw); }

  KGBElt x() const; // reconstruct |x| component
}; // |paramin|

} // |namespace ext_block|

} // |namespace atlas|
#endif
