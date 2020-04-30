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
#include "gradings.h" // type |gradings::Status| must be complete
#include "bruhat.h" // type |BruhatOrder| must be complete for destructor
#include "subsystem.h"
#include "repr.h"

namespace atlas {


namespace repr {
  class Rep_context;
  class common_context;
  class StandardReprMod;
}

namespace ext_block {
  class ext_block;
}

namespace repr {

class Repr_mod_entry
{ KGBElt x; RankFlags y, mask;
public:
  Repr_mod_entry(const Rep_context& rc, const StandardReprMod& srm);

  StandardReprMod srm(const Rep_context& rc,const RatWeight& gamma_mod_1) const;

  // obligatory fields for hashable entry
  using Pooltype =  std::vector<Repr_mod_entry>;
  size_t hashCode(size_t modulus) const
  { return (5*x-11*(y&mask).to_ulong())&(modulus-1); }
  bool operator !=(const Repr_mod_entry& o) const
  { return x!=o.x or (y&mask)!=(o.y&o.mask); }

}; // |Repr_mod_entry|

} // |namespace repr|

namespace blocks {


// a class for blocks of (possibly non integral) parameters
class common_block : public Block_base
{
  const Rep_context& rc; // accesses many things, including KGB set for x

  const RatWeight gamma_mod_1;
  const SubSystem integral_sys;

  std::vector<TorusPart> y_bits; // as in |StandardRepr|, indexed by |y|

  y_entry::Pooltype y_pool;
  y_part_hash y_hash;  // hash table allows storing |y| parts by index

  // hash structure to allow rapid lookup of |(x,y)| index pairs
  block_hash xy_hash;

  // hash structure to facilitate lookup of elements in |StandardReprMod| form
  using repr_hash = HashTable<repr::Repr_mod_entry,BlockElt>;
  repr::Repr_mod_entry::Pooltype z_pool;
  repr_hash srm_hash;

  std::unique_ptr<ext_block::ext_block> extended;

  // group small data members together:
  KGBElt highest_x,highest_y; // maxima over this block
  const bool generated_as_full_block; // tells which constructor was used

 public:

  // constructor and destructor
  common_block // full block
    (const repr::Rep_context& rc,
     const repr::StandardReprMod& srm, // not modified, no "making dominant"
     BlockElt& entry_element	// set to block element matching input
    );
  common_block // partial block
    (const repr::Rep_table& rt,
     const repr::common_context& ctxt,
     containers::sl_list<unsigned long>& elements,
     const RatWeight& gamma_mod_1);
  ~common_block(); // cleans up |*extended|, so inline definition impossible

  // accessors that get values via |rc|
  const repr::Rep_context& context() const { return rc; }
  const RootDatum& root_datum() const;
  const SubSystem& integral_subsystem() const { return integral_sys; }
  const InnerClass& inner_class() const;
  const InvolutionTable& involution_table() const;
  RealReductiveGroup& real_group() const;

  bool is_full () const { return generated_as_full_block; }

  RatWeight gamma_mod1 () const { return gamma_mod_1; }
  // simple coroots of |sub| singular for |gamma|
  RankFlags singular (const RatWeight& gamma) const;

  // with |gamma| unknown, only the difference |gamma-lambda| is meaningful
  RatWeight gamma_lambda(BlockElt z) const;

  BlockElt lookup(const repr::StandardReprMod& srm) const;
  BlockElt lookup(KGBElt x, const RatWeight& gamma_lambda) const;

  repr::StandardReprMod representative (BlockElt z) const
  { return repr::StandardReprMod::build(rc,gamma_mod_1,x(z),gamma_lambda(z)); }
  repr::StandardRepr sr (BlockElt z,const RatWeight& gamma) const;

  ext_gens fold_orbits(const WeightInvolution& delta) const;

  // manipulators
  void swallow // integrate an older partial block, with mapping of elements
    (common_block&& sub, const BlockEltList& embed);
  ext_block::ext_block& extended_block(const WeightInvolution& delta);

  void set_Bruhat
  (containers::sl_list<std::pair<BlockElt,BlockEltList> >&& partial_Hasse);

  // virtual methods
  virtual KGBElt max_x() const { return highest_x; } // might not be final |x|
  virtual KGBElt max_y() const { return highest_y; }

  virtual std::ostream& print // defined in block_io.cpp
    (std::ostream& strm, BlockElt z,bool as_invol_expr) const;


 private:
/*
 sort by increaing length (after reversing if |reverse_length|), within equal
 length groups so by |x(z)|. Permute tables correspondingly
*/
  void sort(unsigned short max_length, bool reverse_length);
  void compute_y_bits();

}; // |class common_block|

} // |namespace blocks|

namespace repr {

// a slight extension of |Rep_context|, fix |delta| for extended representations
class Ext_rep_context : public Rep_context
{
  const WeightInvolution d_delta;
public:
  explicit Ext_rep_context (const repr::Rep_context& rc); // default twisting
  Ext_rep_context (const repr::Rep_context& rc, const WeightInvolution& delta);

  const WeightInvolution& delta () const { return d_delta; }

}; // |class Ext_rep_context|

// another extension of |Rep_context|, fix integral system for common block
class common_context : public Rep_context
{
  const RootDatum integr_datum; // intgrality datum
  const SubSystem sub; // embeds |integr_datum| into parent root datum
public:
  common_context (RealReductiveGroup& G, const SubSystem& integral);

  // accessors
  const RootDatum& id() const { return integr_datum; }
  const SubSystem& subsys() const { return sub; }

  // methods for local common block construction, as in |Rep_context|
  // however, the generator |s| is interpreted for the |integr_datum|
  StandardReprMod cross (weyl::Generator s, const StandardReprMod& z) const;
  StandardReprMod down_Cayley(weyl::Generator s, const StandardReprMod& z) const;
  StandardReprMod up_Cayley(weyl::Generator s, const StandardReprMod& z) const;
  std::pair<gradings::Status::Value,bool> // status and whether a descent/type 1
    status(weyl::Generator s, KGBElt x) const; // with |s| for |integr_datum|
  bool is_parity (weyl::Generator s, const StandardReprMod& z) const;

  Weight to_simple_shift(InvolutionNbr theta, InvolutionNbr theta_p,
			 RootNbrSet pos_to_neg) const; // |pos_to_neg| by value

}; // |class common_context|

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
struct paramin // allow public member access; methods ensure no invariants
{
  const repr::Ext_rep_context& ctxt;
  TwistedInvolution tw; // implicitly defines $\theta$

  Coweight l; // with |tw| gives a |GlobalTitsElement|; lifts its |t|
  RatWeight gamma_lambda; // lift of $\gamma-\lambda$ value in a |StandardRepr|
  Weight tau; // a solution to $(1-\theta)*\tau=(1-\delta)gamma_\lambda$
  Coweight t; // a solution to $t(1-theta)=l(\delta-1)$
  bool flipped; // whether tensored with the flipping representation

  paramin (const repr::Ext_rep_context& ec, const TwistedInvolution& tw,
	   RatWeight gamma_lambda, Weight tau, Coweight l, Coweight t,
	   bool flipped=false);

  // default extension choice:
  paramin (const repr::Ext_rep_context& ec,
	   KGBElt x, const RatWeight& gamma_lambda, bool flipped=false);
  static paramin default_extend
  (const repr::Ext_rep_context& ec, const repr::StandardRepr& sr);

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

  const repr::Rep_context& rc() const { return ctxt; } // reference base object
  const WeightInvolution& delta () const { return ctxt.delta(); }
  const WeightInvolution& theta () const;

  KGBElt x() const; // reconstruct |x| component
  // underlying unextended representation
  repr::StandardRepr restrict(const RatWeight& gamma) const;
}; // |paramin|

} // |namespace ext_block|

} // |namespace atlas|
#endif
