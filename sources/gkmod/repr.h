/*
  This is repr.h

  Copyright (C) 2009-2012 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/
// Classes and utilities for manipulating representations

#ifndef REPR_H  /* guard against multiple inclusions */
#define REPR_H

#include <iostream>

#include "../Atlas.h"

#include "matrix.h"	// containment
#include "ratvec.h"	// containment

#include "rootdata.h" // for |rho|, so |Rep_context::lambda| can be inlined

#include "innerclass.h" // inlines
#include "realredgp.h"	// inlines

#include "hashtable.h"
#include "free_abelian.h"
#include "arithmetic.h" // |SplitInteger|
#include "gradings.h"
#include "subsystem.h"
#include "polynomials.h"

namespace atlas {

namespace repr {

class common_context;

/*
  A parameter of a standard representation is determined by a triplet
  $(x,\tilde\lambda,\gamma)$, where $x\in K\\backslash G/B$ for our fixed real
  form (determining amongs others an involution $\thata$ of $X^*$),
  $\tilde\lambda$ is a genuine character of the $\rho$-cover of $H^{\theta_x}$,
  and $\gamma$ is a character of the complex Lie algebra $h$. The latter two
  values are related; $(1+\theta)\gamma=(1+\theta)\tilde\lambda$, in other words
  $\gamma-\tilde\lambda$ is fixed by $-\theta$; the projection $\lambda_0$ of
  $\tilde\lambda$ or $\gamma$ on the $+1$-eigenspace of $\theta$ is determined
  by this relation. The difference $\gamma-\lambda_0$, i.e., the projection of
  $\gamma$ on the $-1$-eigenspace, is called $\nu$. This $\nu$ is the
  information one is adding here with respect to the values encoded in the
  |standardrepk::StandardRepK| type. The part of $\tilde\lambda$ that is
  independent of $\lambda_0$ is its "torsion part" (due to the disconnectedness
  of $H(R)_c$), which would be represented in the |Block| structure by the
  |TorusPart| component of the |TitsElt| of the dual KGB-element ($y$). In fact
  we also convert it internally to a |TorusPart| here, called |y_bits|. It
  represents an element of the quotient of $(X^*)^{-\theta}$ by the image
  $(1-\theta)X^*$; from it we can recover $\tilde\lambda-\lambda_0$, using
  |involutions::InvolutionTable::unpack|.

  In principle $\gamma$ could take any complex values compatible with
  $\tilde\lambda$. But we are only interested in real values, and in fact
  record a rational value, because we do exact computations, and because all
  interesting phenomena take place at rational infinitesimal character.

  We add a field |height| that is determined by the other ones (in the context
  of a given real reductive group), indeed just by $\lambda_0$. It is useful
  because this is the first statistic used for sorting parameters when they
  are combinined into "parameter polynomials", and this will be faster if we
  don't need the somewhat lengthy computation of the height (twice) for each
  comparison made. On the other hand we now need to compute the height once
  for each parameter constructed, even if it never enters a polynomial.
*/

class StandardRepr
{
  friend class Rep_context;

 protected:
  KGBElt x_part;
  unsigned int hght; // determined by other fields; mainly for (fast) sorting
  TorusPart y_bits; // torsion part of $\lambda$
  RatWeight infinitesimal_char; // $\gamma$ (determines free part of $\lambda$)

  // one should call constructor from |Rep_context| only
  StandardRepr (KGBElt x,TorusPart y,const RatWeight& gamma,unsigned int h)
    : x_part(x), hght(h), y_bits(y), infinitesimal_char(gamma)
  { infinitesimal_char.normalize(); } // to ensure this class invariant

 public:

  const RatWeight& gamma() const { return infinitesimal_char; }
  KGBElt x() const { return x_part; }
  const TorusPart& y() const { return y_bits; }
  unsigned int height() const { return hght; }

  bool operator== (const StandardRepr&) const;

// special members required by HashTable

  typedef std::vector<StandardRepr> Pooltype;
  bool operator!=(const StandardRepr& another) const
    { return not operator==(another); }
  size_t hashCode(size_t modulus) const;
}; // |class StandardRepr|

// a variation that only differs in hashing |infinitesimal_char| modulo $X^*$
class StandardReprMod
{
  friend class Rep_context;

  KGBElt x_part;
  TorusPart y_bits; // torsion part of $\lambda$
  RatWeight inf_char_mod_1; // coset rep. of $\gamma$ in $X^*_\Q / X^*$

  StandardReprMod (StandardRepr&& sr); // private raw constructor

 public:
  // when building, we force integral parts of |gamma_mod1| components to zero
  static StandardReprMod mod_reduce
    (const Rep_context& rc,const StandardRepr& sr);
  static StandardReprMod build
    (const Rep_context& rc, const RatWeight& gamma_mod_1, // must be reduced
     KGBElt x, const RatWeight& gam_lam);

  const RatWeight& gamma_mod1() const { return inf_char_mod_1; }
  KGBElt x() const { return x_part; }
  const TorusPart& y() const { return y_bits; }

  bool operator== (const StandardReprMod& other) const
  { return x_part==other.x_part and y_bits==other.y_bits
    and inf_char_mod_1==other.inf_char_mod_1; }
  typedef std::vector<StandardReprMod> Pooltype;
  bool operator!=(const StandardReprMod& another) const
    { return not operator==(another); }
  size_t hashCode(size_t modulus) const; // this one ignores $X^*$ too
}; // |class StandardReprMod|

class Repr_mod_entry
{ KGBElt x; RankFlags y, mask;
public:
  Repr_mod_entry(const Rep_context& rc, const StandardReprMod& srm);

  StandardReprMod srm(const Rep_context& rc,const RatWeight& gamma_mod_1) const;

  unsigned long y_stripped() const { return(y&mask).to_ulong(); }

  // obligatory fields for hashable entry
  using Pooltype =  std::vector<Repr_mod_entry>;
  size_t hashCode(size_t modulus) const
  { return (5*x-11*y_stripped())&(modulus-1); }
  bool operator !=(const Repr_mod_entry& o) const
    { return x!=o.x or y_stripped()!=o.y_stripped(); }

}; // |Repr_mod_entry|

// This class stores the information necessary to interpret a |StandardRepr|
class Rep_context
{
  RealReductiveGroup& G;
  const KGB& KGB_set;

 public:
  explicit Rep_context(RealReductiveGroup &G);

  // accessors
  RealReductiveGroup& real_group() const { return G; }
  const InnerClass& inner_class() const { return G.innerClass(); }
  const RootDatum& root_datum() const { return G.root_datum(); }
  const TwistedWeylGroup& twisted_Weyl_group() const
    { return G.twistedWeylGroup(); }
  const KGB& kgb() const { return KGB_set; }
  const RatCoweight& g_rho_check() const { return G.g_rho_check(); }
  RatCoweight g() const { return G.g(); }
  size_t rank() const;

  const TwistedInvolution involution_of_Cartan(size_t cn) const;

  RatWeight gamma // compute (representative of) infinitesimal character
    (KGBElt x, const Weight& lambda_rho, const RatWeight& nu) const;
  StandardRepr sr_gamma // use this one when infinitesimal character is known
    (KGBElt x, const Weight& lambda_rho, const RatWeight& gamma) const;
  StandardRepr sr // construct parameter from |(x,\lambda,\nu)| triplet
    (KGBElt x, const Weight& lambda_rho, const RatWeight& nu) const
    { return sr_gamma(x,lambda_rho,gamma(x,lambda_rho,nu)); }
  StandardRepr sr (const StandardReprMod& srm, const RatWeight& gamma) const;

  StandardRepr
    sr(const standardrepk::StandardRepK& srk,
       const standardrepk::SRK_context& srkc,
       const RatWeight& nu) const;

  // component extraction
  const WeightInvolution& theta (const StandardRepr& z) const;

  Weight lambda_rho(const StandardRepr& z) const;
  RatWeight lambda(const StandardRepr& z) const // half-integer
  { return rho(root_datum()).normalize()+lambda_rho(z); }
  RatWeight gamma_lambda
    (InvolutionNbr i_x, const TorusPart& y_bits, const RatWeight& gamma) const;
  RatWeight gamma_lambda(const StandardReprMod& z) const;
  RatWeight gamma_0 // infinitesimal character deformed to $\nu=0$
    (const StandardRepr& z) const;

  RatWeight nu(const StandardRepr& z) const; // rational, $-\theta$-fixed

  // the value of $\exp_{-1}(\gamma-\lambda)$ is $y$ value in a |common_block|
  TorusElement y_as_torus_elt(const StandardRepr& z) const;
  TorusElement y_as_torus_elt(const StandardReprMod& z) const;

  // attributes; they set |witness| only in case they return |false|
  bool is_standard  // whether $I(z)$ is non-virtual: gamma imaginary-dominant
    (const StandardRepr& z, RootNbr& witness) const; // simply-imaginary witness
  bool is_dominant  // whether |gamma| is dominant
    (const StandardRepr& z, RootNbr& witness) const; // simple witness
  bool is_nonzero  // whether $I(z)!=0$: no singular simply-imaginary compact
    (const StandardRepr& z, RootNbr& witness) const; // simply-imaginary witness
  bool is_normal // whether |z==normalise(z)|: has no singular complex descents
    (const StandardRepr& z) const; // complex simple witness
  bool is_semifinal  // whether $I(z)$ unrelated by Hecht-Schmid to more compact
    (const StandardRepr& z, RootNbr& witness) const; // singular real witness
  bool is_final // dominant nonzero without singular descents: all of the above
    (const StandardRepr& z) const;
  bool is_oriented(const StandardRepr& z, RootNbr alpha) const;
  unsigned int orientation_number(const StandardRepr& z) const;

  bool is_twist_fixed(StandardRepr z, const WeightInvolution& delta) const;
  bool is_twist_fixed(const StandardRepr& z) const
  { return is_twist_fixed(z,inner_class().distinguished()); }

  void make_dominant(StandardRepr& z) const; // ensure |z.gamma()| dominant

  // in addition to |make_dominant| apply any singular complex descents
  void normalise(StandardRepr& z) const; // which ensures a normalised form

  bool equivalent(StandardRepr z0, StandardRepr z1) const; // by value

  // deforming the $\nu$ component
  StandardRepr& scale(StandardRepr& sr, const Rational& f) const;
  StandardRepr& scale_0(StandardRepr& sr) const;

  RationalList reducibility_points(const StandardRepr& z) const; // normalised

  // the following take |z| by value, modifying and in some cases returning it
  StandardRepr cross(weyl::Generator s, StandardRepr z) const;
  StandardRepr Cayley(weyl::Generator s, StandardRepr z) const;
  StandardRepr inv_Cayley(weyl::Generator s, StandardRepr z) const;
  StandardRepr inner_twisted(StandardRepr z) const;
  StandardRepr twisted(StandardRepr z, const WeightInvolution& delta) const;

  StandardRepr cross(const Weight& alpha, StandardRepr z) const;
  StandardRepr any_Cayley(const Weight& alpha, StandardRepr z) const;

  class compare
  { Coweight level_vec; // linear form to apply to |gamma| for ordering
  public:
    compare (const Coweight& lv) : level_vec(lv) {}

    bool operator()(const StandardRepr& r,const StandardRepr& s) const;
  }; // |compare|

  compare repr_less() const;

  using poly = Free_Abelian<StandardRepr,Split_integer,compare>;

  poly scale(const poly& P, const Rational& f) const;
  poly scale_0(const poly& P) const;

  containers::sl_list<StandardRepr> finals_for // like |Block_base::finals_for|
    (StandardRepr z) const; // by value
  poly expand_final(StandardRepr z) const; // the same, as |poly| (by value)

  std::ostream& print (std::ostream&,const StandardRepr& z) const;
  std::ostream& print (std::ostream&,const poly& P) const;

 private:
  // make integrally dominant, with precomputed integral subsystem; return path
  WeylWord make_dominant(StandardRepr& z,const SubSystem& subsys) const;
  void singular_cross (weyl::Generator s,StandardRepr& z) const;
  void to_singular_canonical(RankFlags gens, StandardRepr& z) const;
  unsigned int height(Weight theta_plus_1_gamma) const;
}; // |Rep_context|

using SR_poly = Rep_context::poly;

/*
  In addition to providing methods inherited from |Rep_context|, the class
  |Rep_table| provides storage for data that was previously computed for
  various related nonzero final |StandardRepr| values.

  The data stored are: the |blocks::common_block| values encountered for this
  real form, together with their KL data if computed, a table of (twisted) full
  deformation formulae, pools of |kl::KLPol| (positive coeffient) and
  |ext_kl::Pol| (integer coeffient) polynomials that blocks may choose to share
  (they can alternatively choose to maintain their local polynomials themselves).

  This class provides methods for their computation, and handles the data
  storage and retrieval.

  The deformation information is associated to parameters, but is unchanged
  under certain changes of these parameters, so to limit memory usage somewhat
  any parameter is moved to its first reducibility point and expressed in
  nonzezro final ones before looking up or computing its deformation.
*/
class Rep_table : public Rep_context
{
  std::vector<StandardRepr> pool;
  HashTable<StandardRepr,unsigned long> hash;
  std::vector<std::pair<SR_poly,SR_poly> > def_formulae; // ordinary, twisted

  std::vector<StandardReprMod> mod_pool;
  HashTable<StandardReprMod,unsigned long> mod_hash;

  std::vector<kl::KLPol> KL_poly_pool;
  KL_hash_Table KL_poly_hash;

  std::vector<ext_kl::Pol> poly_pool;
  ext_KL_hash_Table poly_hash;

  containers::sl_list<blocks::common_block> block_list;
  using bl_it = containers::sl_list<blocks::common_block>::iterator;
  std::vector<std::pair<bl_it, BlockElt> > place;

 public:
  Rep_table(RealReductiveGroup &G);
  ~Rep_table();
  // both defined out of line because of implicit use |common_block| destructor

  ext_KL_hash_Table* shared_poly_table () { return &poly_hash; }

  const StandardReprMod& srm(unsigned long n) const { return mod_pool[n]; }

  unsigned short length(StandardRepr z); // by value

  unsigned long parameter_number (StandardRepr z) const { return hash.find(z); }
  const SR_poly& deformation_formula(unsigned long h) const
    { assert(h<def_formulae.size()); return def_formulae[h].first; }
  const SR_poly& twisted_deformation_formula(unsigned long h) const
    { assert(h<def_formulae.size()); return def_formulae[h].second; }

  blocks::common_block& lookup_full_block
    (StandardRepr& sr,BlockElt& z); // |sr| is by reference; will be normalised

  blocks::common_block& lookup // constuct only partial block if necessary
    (StandardRepr& sr,BlockElt& z); // |sr| is by reference; will be normalised

  SR_poly KL_column_at_s(StandardRepr z); // by value
  containers::simple_list<std::pair<BlockElt,kl::KLPol> >
    KL_column(StandardRepr z); // by value
  SR_poly twisted_KL_column_at_s(StandardRepr z); // by value

  SR_poly deformation_terms
    (blocks::common_block& block, BlockElt y, const RatWeight& gamma);
#if 0
  SR_poly deformation_terms (unsigned long sr_hash) const;
  // once a parameter has been entered, we can compute this without a block
#endif

  SR_poly deformation(const StandardRepr& z);

  SR_poly twisted_deformation_terms
    (blocks::common_block& block, ext_block::ext_block& eblock,
     BlockElt y, RankFlags singular, const RatWeight& gamma);
#if 0
  SR_poly twisted_deformation_terms (unsigned long sr_hash) const;
  // once a parameter has been entered, we can compute this without a block
#endif

  blocks::common_block& add_block_below // partial; defined in common_blocks.cpp
    (const common_context&, const StandardReprMod& srm, BitMap* subset);

  SR_poly twisted_deformation(StandardRepr z); // by value

 private:
  void block_erase (bl_it pos); // erase from |block_list| in safe manner
  unsigned long formula_index (const StandardRepr&);
  unsigned long add_block(const StandardReprMod&); // full block
  class Bruhat_generator; // helper class: internal |add_block_below| recursion

}; // |Rep_table|


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
  This class is for |paramin| below what |ext_block::context| is for
  |ext_block::param|: it holds relevant values that remain fixed across an
  |extended block|. Data fields that are removed with respect to
  |ext_block::param| are |d_gamma|, |lambda_shifts|. Methods that are absent:
  |gamma|, |lambda_shift|
*/
class Ext_common_context : public common_context
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

// 				Functions

// shift in $\lambda$ component involved in non-simple Cayleys (and crosses)
// gets added to |lambda_rho| in imaginary cases, subtracted in real cases
Weight Cayley_shift (const InnerClass& G,
		     InvolutionNbr theta_upstairs, // at the more split Cartan
		     const WeylWord& to_simple); // acting from the left

SR_poly twisted_KL_column_at_s
  (const Rep_context& rc, StandardRepr z, const WeightInvolution& delta);

} // |namespace repr|

} // |namespace atlas|

#endif
