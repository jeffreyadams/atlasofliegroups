/*
  This is repr.h

  Copyright (C) 2009-2022 Marc van Leeuwen
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
#include "K_repr.h"  // |K_type|, |K_type_pol|

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
  |K_repr::K_type| type. The part of $\tilde\lambda$ that is independent of
  $\lambda_0$ is its "torsion part" (due to the disconnectedness of $H(R)_c$),
  which would be represented in the |Block| structure by the |TorusPart|
  component of the |TitsElt| of the dual KGB-element ($y$). In fact we also
  convert it internally to a |TorusPart| here, called |y_bits|. It represents an
  element of the quotient of $(X^*)^{-\theta}$ by the image $(1-\theta)X^*$;
  from it we can recover $\tilde\lambda-\lambda_0$, using
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

  // precomputed height is first criterion for ordering
  // therefore comparison involves no further |Rep_context| dependency
  bool operator<(const StandardRepr&) const;
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
  RatWeight gamma_lambda() const & { return gamlam; }
  RatWeight&& gamma_lambda() && { return std::move(gamlam); }

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

  K_repr::K_type sr_K(KGBElt x, Weight lambda_rho) const; // 2nd by value
  K_repr::K_type sr_K(const StandardRepr& sr) const // restriction to K
  { return {sr.x(),lambda_rho(sr),sr.height()}; }
  Weight theta_plus_1_lambda (const K_repr::K_type& t) const;

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

  StandardRepr sr (const K_repr::K_type& t) const
  { return sr(t.x(),t.lambda_rho(),RatWeight(t.lambda_rho().size())); }

  // Handling of |K_repr::K_type| values

  bool is_standard  // whether $(1+\theta)*\lambda| is imaginarily dominant
    (const K_repr::K_type& z) const;
  bool is_dominant  // whether $(1+\theta)*\lambda| is dominant
    (const K_repr::K_type& z) const;
  bool is_canonical // whether involution is canonical one for the Cartan class
    (const K_repr::K_type& z) const;
  bool is_theta_stable // absence of complex descents (theta-stable parabolic)
    (const K_repr::K_type& z) const;
  bool is_nonzero  // absence of singular compact simply-imaginary roots
    (const K_repr::K_type& z) const;
  bool is_semifinal  // absence of real parity roots
    (const K_repr::K_type& z) const;
  // call next predicate only for standard, dominant, nonzero, semifinal |z|
  bool is_normal // absence of singular complex descents
    (const K_repr::K_type& z) const; // complex simple witness
  bool is_final // standard, dominant, no singular descents of any kind
    (const K_repr::K_type& z) const;
  bool equivalent(K_repr::K_type z0, K_repr::K_type z1) const; // by value

  // "standard" (|lambda| is dominant for imaginary coroots) is assumed here

  // ensure that |(1+theta)*lambda| is dominant also for complex coroots
  void make_dominant (K_repr::K_type& z) const;
  // exhaust simple complex descents for the involution
  void make_theta_stable (K_repr::K_type& z) const;
  // make the involution the preferred one for the Cartan class
  void to_canonical_involution (K_repr::K_type& z, RankFlags gens) const;
  void to_canonical_involution (K_repr::K_type& z) const
  { return to_canonical_involution
      (z,RankFlags(constants::lMask[root_datum().semisimple_rank()])); }
  void normalise(K_repr::K_type& z) const; // which ensures a normalised form

  simple_list<std::pair<K_repr::K_type,int> >
    finals_for(K_repr::K_type t) const;

  // conservative estimate for lowest height that can be obtained from |lambda|
  level height_bound(RatWeight lambda) const; // "projecting to dominant cone"
  // apart from producing a result, these two methods also make |t| theta-stable
  sl_list<K_repr::K_type> KGP_set (K_repr::K_type& t) const;
  K_repr::K_type_pol K_type_formula (K_repr::K_type& t,level cutoff) const;
  K_repr::K_type_pol K_type_formula (SRK_context& srk, // for pruning
				     K_repr::K_type& t,level cutoff) const;
  K_repr::K_type_pol branch(K_repr::K_type_pol P, level cutoff) const;
  K_repr::K_type_pol branch(SRK_context& srk, // use projections for pruning
			    K_repr::K_type_pol P, level cutoff) const;

  // parameter component extraction
  const WeightInvolution& theta (const StandardRepr& z) const;

  Weight lambda_rho(const StandardRepr& z) const;
  RatWeight lambda(const StandardRepr& z) const // half-integer
  { return rho(root_datum())+lambda_rho(z); }
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
  RatWeight&& gamma_lambda(StandardReprMod&& z) const
  { return std::move(z).gamma_lambda(); }
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
    (const StandardRepr& z) const;
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
  K_repr::K_type scale_0(const StandardRepr& sr) const
  { return sr_K(sr); } // used to return |StandardRepr|, now is just an alias

  RatNumList reducibility_points(const StandardRepr& z) const; // normalised

  // the following take |z| by value, modifying and in some cases returning it
  StandardRepr cross(weyl::Generator s, StandardRepr z) const;
  StandardRepr Cayley(weyl::Generator s, StandardRepr z) const;
  StandardRepr inner_twisted(StandardRepr z) const;
  StandardRepr twisted(StandardRepr z, const WeightInvolution& delta) const;

  StandardRepr cross(const Weight& alpha, StandardRepr z) const;
  StandardRepr any_Cayley(const Weight& alpha, StandardRepr z) const;

  using poly = Free_Abelian<StandardRepr,Split_integer>;

  poly scale(const poly& P, const RatNum& f) const;
  K_repr::K_type_pol scale_0(const poly& P) const;

  simple_list<std::pair<StandardRepr,int> >
    finals_for (StandardRepr z) const; // like |finals_for| K-type (by value)
  poly expand_final(StandardRepr z) const; // the same, as |poly| (by value)

  Weight to_simple_shift(InvolutionNbr theta, InvolutionNbr theta_p,
			 RootNbrSet pos_to_neg) const; // |pos_to_neg| by value

 private:
  // compute $\check\alpha\dot(1+\theta_x)\lambda$, with $(x,\lambda)$ from $t$
  int theta_plus_1_eval (const K_repr::K_type& t, RootNbr alpha) const;
  RankFlags singular_simples (const StandardRepr& z) const;
  WeylWord complex_descent_word (KGBElt x, RankFlags singulars) const;
  // make integrally dominant, with precomputed integral subsystem; return path
  WeylWord make_dominant(StandardRepr& z,const SubSystem& subsys) const;
  void complex_crosses (StandardRepr& z, const WeylWord& ww) const;
  void to_singular_canonical(RankFlags gens, StandardRepr& z) const;
  level height(Weight theta_plus_1_gamma) const;

  K_repr::K_type_pol monomial_product
    (const K_repr::K_type_pol& P, const Weight& e) const;
}; // |Rep_context|

/* In internal computations for deformation, it will be important to have a
   quite compact representation of (twisted) deformation formulas. To this end
   we shall make (inside a |Rep_table|) a table of K-types encountered (the
   number of possible K-types for a given real form is fairly limited), and
   represent linear coefficients as a map from such K-types to the ring of
   coefficients (in practice |Split_integer|). Moreover, for storing formulas as
   collections of key-value pairs, we wish to avoid using the |Free_Abelian|
   container, derived from |std::map| which requires a lot of additional memory
   per node stored, preferring to it our tailor made |Free_Abelian_light| which
   has no per-node memory overhead.
*/

using K_type_nr = unsigned int; // hashed in |Rep_table| below
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

  K_repr::K_type::Pooltype K_type_pool;
  HashTable<K_repr::K_type,K_type_nr> K_type_hash;

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

  K_repr::K_type stored_K_type(K_type_nr i) const
  { return K_type_pool[i].copy(); }

  // a signed multiset of final parameters needed to be considered (i.e., their
  // deformations to $\nu=0$ included) when deforming |y| a bit towards $\nu=0$
  sl_list<std::pair<StandardRepr,int> > deformation_terms
    (blocks::common_block& block, BlockElt y,
     const RatWeight& diff, const RatWeight& gamma);

  // full deformation to $\nu=0$ of |z|
  const K_type_poly& deformation(StandardRepr z); // by value

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


// another extension of |Rep_context|, fix integral system for common block
class common_context
{
  const Rep_context& rep_con;
  unsigned int int_sys_nr;
  const SubSystem& sub; // embeds |integr_datum| into parent root datum
public:
  common_context (const Rep_context& rc, const RatWeight& gamma);
  common_context (const Rep_context& rc, const SubSystem& integral);

  // accessors
  const Rep_context& rc() const { return rep_con; }
  const KGB& kgb() const { return rep_con.kgb(); }
  const InvolutionTable& involution_table() const
    { return rep_con.involution_table(); }
  const RootDatum& full_root_datum() const { return rep_con.root_datum(); }
  const SubSystem& subsys() const { return sub; }

  // methods for local common block construction, as in |Rep_context|
  // however, the generator |s| is interpreted for |subsys()|
  StandardReprMod cross (weyl::Generator s, const StandardReprMod& z) const;
  StandardReprMod down_Cayley(weyl::Generator s, const StandardReprMod& z) const;
  StandardReprMod up_Cayley(weyl::Generator s, const StandardReprMod& z) const;
  std::pair<gradings::Status::Value,bool> // status and whether a descent/type 1
    status(weyl::Generator s, KGBElt x) const;
  bool is_parity (weyl::Generator s, const StandardReprMod& z) const;

}; // |class common_context|

/*
  This class holds relevant values that remain fixed across an extended block,
  and are unchaged under integral changes to |gamma|. Data fields that are
  removed with respect to |ext_block::context| are |d_gamma|, |lambda_shifts|.
  Methods that are absent: |gamma|, |lambda_shift|
*/
class Ext_rep_context
{
  const Rep_context& rep_con;
  const WeightInvolution d_delta;
  Permutation pi_delta; // permutation of |delta| on roots of full root datum
  RootNbrSet delta_fixed_roots; // as subset of full root system
  weyl::Twist twist; // of the full Dynkin diagram

 public:
  explicit Ext_rep_context (const repr::Rep_context& rc); // default twisting
  Ext_rep_context (const Rep_context& rc, const WeightInvolution& delta);

  // accessors
  const Rep_context& rc() const { return rep_con; }
  const WeightInvolution& delta () const { return d_delta; }
  RootNbr delta_of(RootNbr alpha) const { return pi_delta[alpha]; }
  const RootNbrSet& delta_fixed() const { return delta_fixed_roots; }
  weyl::Generator twisted(weyl::Generator s) const { return twist[s]; }

  const RootDatum& root_datum () const { return rep_con.root_datum(); }
  const InnerClass& inner_class() const { return rep_con.inner_class(); }
  RealReductiveGroup& real_group() const { return rep_con.real_group(); }
  const RatCoweight& g_rho_check() const { return rep_con.g_rho_check(); }
  RatWeight gamma_lambda(const StandardReprMod& z) const
  { return rep_con.gamma_lambda(z); }

  Weight to_simple_shift(InvolutionNbr theta, InvolutionNbr theta_p,
			 RootNbrSet pos_to_neg) const // |pos_to_neg| by value
  { return rep_con.to_simple_shift(theta,theta_p,pos_to_neg); }
  // whether positive $\alpha$ has $\theta(\alpha)\neq\pm(1|\delta)(\alpha)$
  bool is_very_complex(InvolutionNbr theta, RootNbr alpha) const;
  bool shift_flip(InvolutionNbr theta, InvolutionNbr theta_p,
		  RootNbrSet pos_to_neg) const; // |pos_to_neg| is by value

}; // |Ext_rep_context|

// 				Functions

// printing functions for various types are declared in basic_io.h

// shift in $\lambda$ component involved in non-simple Cayleys (and crosses)
// gets added to |lambda_rho| in imaginary cases, subtracted in real cases
Weight Cayley_shift (const InnerClass& G,
		     InvolutionNbr theta_upstairs, // at the more split Cartan
		     const WeylWord& to_simple); // acting from the left

SR_poly twisted_KL_column_at_s
  (const Rep_context& rc, StandardRepr z, const WeightInvolution& delta);

K_repr::K_type_pol export_K_type_pol(const Rep_table& rt,const K_type_poly& P);

} // |namespace repr|

} // |namespace atlas|

#endif
