/*
  This is repr.h

  Copyright (C) 2009-2020 Marc van Leeuwen
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
#include "kgb.h"
#include "subsystem.h"
#include "polynomials.h"

namespace atlas {

namespace repr {

class common_context;

/*
  A parameter of a standard representation is determined by a triplet
  $(x,\tilde\lambda,\gamma)$, where $x\in K\\backslash G/B$ for our fixed real
  form (determining amongst others an involution $\theta$ of $X^*$),
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

  // one should call constructors from |Rep_context| only, so they are |private|
  StandardRepr(KGBElt x,TorusPart y,const RatWeight& gamma,unsigned int h)
    : x_part(x), hght(h), y_bits(y), infinitesimal_char(gamma)
  { infinitesimal_char.normalize(); } // to ensure this class invariant
  StandardRepr(KGBElt x,TorusPart y, RatWeight&& gamma,unsigned int h)
    : x_part(x), hght(h), y_bits(y), infinitesimal_char(std::move(gamma))
  { infinitesimal_char.normalize(); } // to ensure this class invariant

 public:

  const RatWeight& gamma() const { return infinitesimal_char; }
  KGBElt x() const { return x_part; }
  const TorusPart& y() const { return y_bits; }
  unsigned int height() const { return hght; }

  bool operator==(const StandardRepr&) const;

// special members required by HashTable

  using Pooltype = std::vector<StandardRepr>;
  bool operator!=(const StandardRepr& another) const
  { return not operator==(another); }
  size_t hashCode(size_t modulus) const;
}; // |class StandardRepr|

/*
  A variation, representing |StandardRepr| up to parallel shift |gamma|,|lambda|
  which suffices for block construction, and given |gamma| we can easily
  reconstruct a full |StandardRepr|.

  We use a trick to assume by parallel shift that |lambda=rho|, so the
  corresponding value |rgl| of |gamma| now also encodes |gamma-lambda|, and we
  do not need to represent |y_bits| at all. In fact doing without the packing
  and unpacking operations for |y_bits| greatly simplifies our operations.
*/
class StandardReprMod
{
  friend class Rep_context;

  KGBElt x_part;
  RatWeight gamlam; // |real_unique(gamma-lambda)|

  StandardReprMod(KGBElt x_part, RatWeight&& gl) // private raw constructor
    : x_part(x_part), gamlam(std::move(gl.normalize())) {}

 public:
  // the raw constructor is always called through one of two pseudo constructors
  static StandardReprMod mod_reduce
    (const Rep_context& rc,const StandardRepr& sr);
  static StandardReprMod build
    (const Rep_context& rc, KGBElt x, RatWeight gam_lam);

  KGBElt x() const { return x_part; }
  RatWeight gamma_lambda() const { return gamlam; }

  // since pseudo constructors map |rgl| to fundamental domain, equality is easy
  bool operator==(const StandardReprMod& other) const
  { return x_part==other.x_part and gamlam==other.gamlam; }

  typedef std::vector<StandardReprMod> Pooltype;
  bool operator!=(const StandardReprMod& another) const
  { return not operator==(another); }
  size_t hashCode(size_t modulus) const; // this one ignores $X^*$ too
}; // |class StandardReprMod|

// hashable type for |StandardReprMod| upto shift orthogonal to integral system
class Reduced_param
{
  KGBElt x;
  unsigned int int_sys_nr;
  int_Vector evs_reduced; // coroots evaluation (|gamlam| mod $(1-theta)X^*$)

public:
  Reduced_param(const Rep_context& rc, const StandardReprMod& srm);

  Reduced_param(Reduced_param&&) = default;

  using Pooltype = std::vector<Reduced_param>;
  bool operator!=(const Reduced_param& p) const
  { return x!=p.x or int_sys_nr!=p.int_sys_nr or evs_reduced!=p.evs_reduced; }
  bool operator==(const Reduced_param& p) const { return not operator!=(p); }
  size_t hashCode(size_t modulus) const; // this one ignores $X^*$ too
}; // |class Reduced_param|

// This class stores the information necessary to interpret a |StandardRepr|
class Rep_context
{
  RealReductiveGroup& G;
  const KGB& KGB_set;

 public:
  explicit Rep_context(RealReductiveGroup &G);

  // accessors
  RealReductiveGroup& real_group() const { return G; }
  InnerClass& inner_class() const { return G.innerClass(); }
  const InvolutionTable& involution_table() const
  { return inner_class().involution_table(); }
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
  StandardRepr sr_gamma // variant of previous which moves from |gamma|
    (KGBElt x, const Weight& lambda_rho, RatWeight&& gamma) const;
  StandardRepr sr // construct parameter from |(x,\lambda,\nu)| triplet
    (KGBElt x, const Weight& lambda_rho, const RatWeight& nu) const
  { return sr_gamma(x,lambda_rho,gamma(x,lambda_rho,nu)); }

  // reconstruct |StandardRep| from |srm| and difference of |gamma_lambda|s
  StandardRepr sr (const StandardReprMod& srm,const RatWeight& gamma)  const;
  StandardRepr sr
    (const StandardReprMod& srm, const RatWeight& diff, const RatWeight& gamma)
    const;

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
  RatWeight gamma_lambda(const StandardRepr& z) const
  { return gamma_lambda(kgb().inv_nr(z.x()),z.y(),z.gamma()); }
  RatWeight gamma_lambda_rho(const StandardRepr& z) const;
  RatWeight gamma_0 // infinitesimal character deformed to $\nu=0$
    (const StandardRepr& z) const;
  RatWeight nu(const StandardRepr& z) const; // rational, $-\theta$-fixed

  // whether real root with index |i| is parity on |z| deformed to $\nu=0$
  // only depends on |x| and |y| parts of |z|, and is defined for all real roots
  bool is_parity_at_0(RootNbr i,const StandardRepr& z) const;
  // method that assumes an integral real coroot with index |i|, uses |gamma|
  bool is_parity(RootNbr i,const StandardRepr& z) const;

  // |StandardReprMod| handling

  StandardReprMod inner_twisted(const StandardReprMod& z) const;

  // offset in $\gamma-\lambda$ from |srm0| with respect to that of |srm1|
  RatWeight offset
    (const StandardReprMod& sr0, const StandardReprMod& srm1) const;
  RatWeight offset (const StandardRepr& sr, const StandardReprMod& srm) const
  { return offset(StandardReprMod::mod_reduce(*this,sr),srm); }
  // auxiliary for |offset|
  // find element in |(1-theta)X^*| with same evaluation on all integral coroots
  Weight theta_1_preimage
    (const RatWeight& offset, const subsystem::integral_datum_item::codec& codec)
    const;
  StandardReprMod& shift(const RatWeight& diff, StandardReprMod& srm) const;
  StandardReprMod shifted(const RatWeight& diff, StandardReprMod srm) const
  { return shift(diff,srm); } // preform |shift| on a copy and return it

  RatWeight gamma_lambda(const StandardReprMod& z) const
  { return z.gamma_lambda(); }
  RatWeight gamma_lambda_rho(const StandardReprMod& z) const
  { return z.gamma_lambda()+rho(root_datum()); }

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

  // do |make_dominant|, and descend through any singular complex descents
  void deform_readjust(StandardRepr& z) const;
  // do the same, ensuring fixed choice of descent-minimum among equivalent ones
  void normalise(StandardRepr& z) const; // which ensures a normalised form

  bool equivalent(StandardRepr z0, StandardRepr z1) const; // by value

  // deforming the $\nu$ component
  StandardRepr scale(StandardRepr sr, const RatNum& f) const; // |sr| by value
  StandardRepr scale_0(StandardRepr sr) const; // |sr| by value

  RatNumList reducibility_points(const StandardRepr& z) const; // normalised

  // the following take |z| by value, modifying and in some cases returning it
  StandardRepr cross(weyl::Generator s, StandardRepr z) const;
  StandardRepr Cayley(weyl::Generator s, StandardRepr z) const;
  StandardRepr inv_Cayley(weyl::Generator s, StandardRepr z) const;
  StandardRepr inner_twisted(StandardRepr z) const;
  StandardRepr twisted(StandardRepr z, const WeightInvolution& delta) const;

  StandardRepr cross(const Weight& alpha, StandardRepr z) const;
  StandardRepr any_Cayley(const Weight& alpha, StandardRepr z) const;

  struct compare
  {
    bool operator()(const StandardRepr& r,const StandardRepr& s) const;
  }; // |compare|

  using poly = Free_Abelian<StandardRepr,Split_integer,compare>;

  poly scale(const poly& P, const RatNum& f) const;
  poly scale_0(const poly& P) const;

  sl_list<StandardRepr> finals_for // like |Block_base::finals_for|
    (StandardRepr z) const; // by value
  poly expand_final(StandardRepr z) const; // the same, as |poly| (by value)

 private:
  // make integrally dominant, with precomputed integral subsystem; return path
  WeylWord make_dominant(StandardRepr& z,const SubSystem& subsys) const;
  void singular_cross (weyl::Generator s,StandardRepr& z) const;
  void to_singular_canonical(RankFlags gens, StandardRepr& z) const;
  unsigned int height(Weight theta_plus_1_gamma) const;
}; // |Rep_context|

using SR_poly = Rep_context::poly;


class K_type // compact representation of parameters at $\nu=0$
{
  KGBElt d_x;
  Weight lam_rho;

public:
  K_type(const Rep_context& rc, const StandardRepr& sr)
    : d_x(sr.x()), lam_rho(rc.lambda_rho(sr)) {}

  K_type(K_type&&) = default;

  KGBElt x () const { return d_x;  }
  const Weight& lambda_rho () const { return lam_rho; }

  StandardRepr sr (const Rep_context& rc) const // represent as full parameter
  { return rc.sr(d_x,lam_rho,RatWeight(lam_rho.size())); }

  bool operator< (const K_type& another) const
  {
    if (d_x!=another.d_x)
      return d_x<another.d_x;
    assert(lam_rho.size()==another.lam_rho.size()); // this is always assumed
    for (unsigned i=0; i<lam_rho.size(); ++i)
      if (lam_rho[i]!=another.lam_rho[i])
	return lam_rho[i]<another.lam_rho[i];
    return false; // we found equality
  }

  using Pooltype = std::vector<K_type>;
  bool operator!= (const K_type& another) const
  { return d_x!=another.d_x or lam_rho!=another.lam_rho; }
  size_t hashCode (size_t modulus) const
  {
    size_t h = 3*d_x;
    for (auto c : lam_rho)
      h = (17*h&(modulus-1)) + c;
    return h&(modulus-1);
  }
}; // |class K_type|

using K_type_nr = unsigned int; // hashed in |Rep_table| below

using K_term_type = std::pair<K_type_nr,Split_integer>;
using K_type_poly = Free_Abelian_light<K_type_nr,Split_integer>;

/*
  A class to serve as key-value pair for deformation formula lookup.

  A key determines a parameter up to variations in $\nu$ that do not change the
  integral part (always rounding down) of the evaluation on any positive coroot;
  the domain of such changes is called an alcove. (Strictly speaking, set of
  associated infinitesimal characters is the intersection of an alcove with an
  affine subspace parallel to the $-1$ eigenspace of the involution; no
  (discrete) variation is allowed in the direction of the $+1$ eigenspace.)

  For compactness we store a sample |StandardRepr|, but the functions used for
  hashing only take into account aspects that are unchanging within the alcove.
 */
class deformation_unit
{
  friend class Rep_table; // while not essential, allows easier instrumenting

  StandardRepr sample;
  K_type_poly untwisted, twisted;
  const Rep_context& rc; // access coroots etc. necessary for alcove testing
public:
  deformation_unit(const Rep_context& rc, const StandardRepr& sr)
  : sample(sr), untwisted(), twisted(), rc(rc) {}
  deformation_unit(const Rep_context& rc, StandardRepr&& sr)
  : sample(std::move(sr)), untwisted(), twisted(), rc(rc) {}

  deformation_unit(deformation_unit&&) = default; // type is only movable

  bool has_deformation_formula() const { return not untwisted.is_zero(); }
  bool has_twisted_deformation_formula() const { return not twisted.is_zero(); }

  size_t def_form_size () const { return untwisted.size(); }
  size_t twisted_def_form_size () const { return twisted.size(); }

  const K_type_poly& def_formula() const       { return untwisted; }
  const K_type_poly& twisted_def_formula() const { return twisted; }

  const K_type_poly& set_deformation_formula (const K_type_poly& formula)
  { return untwisted = formula.copy(); }
  const K_type_poly& set_deformation_formula (K_type_poly&& formula)
  { return untwisted = std::move(formula); }
  const K_type_poly& set_twisted_deformation_formula (const K_type_poly& formula)
  { return twisted = formula.copy(); }
  const K_type_poly& set_twisted_deformation_formula (K_type_poly&& formula)
  { return twisted = std::move(formula); }

// special members required by HashTable
  using Pooltype = std::vector<deformation_unit>;
  bool operator!=(const deformation_unit& another) const; // distinct alcoves?
  size_t hashCode(size_t modulus) const; // value depending on alcove only
}; // |class deformation_unit|

/*
  In addition to providing methods inherited from |Rep_context|, the class
  |Rep_table| provides storage for data that was previously computed for
  various related nonzero final |StandardRepr| values.

  The data stored are:
  * the |blocks::common_block| values encountered for this real form, together
    with their KL data if computed; these are accessed via |block_list| and
   |place|
  * a table of (twisted) full deformation formulae, associated to alcoves;
    these are stored in |pool| and accessed through |alcove_hash|
  * a table |reduced_hash| of |Reduced_param| values encountered; this table
    guards the creation of blocks, as parameters sharing a |Reduced_param| value
    can share their block: for finding a block, lookup in |reduced_hash| is
    performed, and if found the index is then used via |place| to find the block
  * a table |K_type_hash| of |K_type| values having been found to occur in
    deformation formulas, which allows compact representation of the latter
  * tables |KL_poly_hash| and |poly_hash| of |kl::KLPol| (positive coeffient)
    respectively |ext_kl::Pol| (integer coeffient) polynomials, which can be
    shared among blocks to reduce the size of their tables of such polynomials
    (they may also choose to maintain their local polynomials themselves).
*/
class Rep_table : public Rep_context
{
  deformation_unit::Pooltype pool; // also stores actual deformation formulae
  HashTable<deformation_unit,unsigned long> alcove_hash;

  Reduced_param::Pooltype reduced_pool;
  HashTable<Reduced_param,unsigned long> reduced_hash;

  K_type::Pooltype K_type_pool;
  HashTable<K_type,K_type_nr> K_type_hash;

  PosPolEntry::Pooltype KL_poly_pool;
  KL_hash_Table KL_poly_hash;

  IntPolEntry::Pooltype poly_pool;
  ext_KL_hash_Table poly_hash;

  sl_list<blocks::common_block> block_list;
  using bl_it = sl_list<blocks::common_block>::iterator;
  std::vector<std::pair<bl_it, BlockElt> > place; // parallel to |reduced_pool|

 public:
  Rep_table(RealReductiveGroup &G);
  ~Rep_table();
  // both defined out of line because of implicit use |common_block| destructor

  ext_KL_hash_Table* shared_poly_table () { return &poly_hash; }

  // the |length| method generates a partial block, for best amortised efficiency
  unsigned short length(StandardRepr z); // by value

  blocks::common_block& lookup_full_block
    (StandardRepr& sr,BlockElt& z); // |sr| is by reference; will be normalised

  blocks::common_block& lookup // constuct only partial block if necessary
    (StandardRepr& sr,BlockElt& z); // |sr| is by reference; will be normalised

  SR_poly KL_column_at_s(StandardRepr z); // by value
  simple_list<std::pair<BlockElt,kl::KLPol> >
    KL_column(StandardRepr z); // by value
  SR_poly twisted_KL_column_at_s(StandardRepr z); // by value

  size_t find_reduced_hash(const StandardReprMod& srm) const
  { return reduced_hash.find(Reduced_param(*this,srm)); }
  size_t match_reduced_hash(const StandardReprMod& srm)
  { return reduced_hash.match(Reduced_param(*this,srm)); }

  StandardRepr K_type_sr(K_type_nr i) { return K_type_pool[i].sr(*this); }

  // a signed multiset of final parameters needed to be taken into account
  // (deformations to $\nu=0$ included) when deforming |y| a bit towards $\nu=0$
  sl_list<std::pair<StandardRepr,int> > deformation_terms
    (blocks::common_block& block, BlockElt y,
     const RatWeight& diff, const RatWeight& gamma);

  // full deformation to $\nu=0$ of |z|
  const K_type_poly& deformation(const StandardRepr& z);

  // like |deformation_terms|; caller multiplies returned coefficients by $1-s$
  sl_list<std::pair<StandardRepr,int> > twisted_deformation_terms
    (blocks::common_block& block, ext_block::ext_block& eblock,
     BlockElt y, // in numbering of |block|, not |eblock|
     RankFlags singular, const RatWeight& diff, const RatWeight& gamma);
#if 0
  SR_poly twisted_deformation_terms (unsigned long sr_hash) const;
  // once a parameter has been entered, we can compute this without a block
#endif

  blocks::common_block& add_block_below // partial; defined in common_blocks.cpp
    (const common_context&, const StandardReprMod& srm, BitMap* subset);

  // full twisted deformation, with |flip| telling whether to multiply by |s|
  const K_type_poly& twisted_deformation(StandardRepr z, bool& flip); // by value

 private:
  void block_erase (bl_it pos); // erase from |block_list| in safe manner
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
class common_context
{
  const Rep_context& rep_con;
  const RootDatum integr_datum; // intgrality datum
  const SubSystem sub; // embeds |integr_datum| into parent root datum
public:
  common_context (const Rep_context& rc, const SubSystem& integral);

  // accessors
  const Rep_context& rc() const { return rep_con; }
  const KGB& kgb() const { return rep_con.kgb(); }
  const InvolutionTable& involution_table() const
    { return rep_con.involution_table(); }
  const RootDatum& full_root_datum() const { return rep_con.root_datum(); }
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
  Ext_common_context (const Rep_context& rc, const WeightInvolution& delta,
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

// printing functions for various types are declared in basic_io.h

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
