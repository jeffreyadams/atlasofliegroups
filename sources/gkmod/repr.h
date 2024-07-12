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

#include "Atlas.h"

#include "matrix.h"	// containment
#include "ratvec.h"	// containment

#include "rootdata.h" // for |rho|, so |Rep_context::lambda| can be inlined

#include "innerclass.h" // inlines
#include "realredgp.h"	// inlines

#include "hashtable.h"
#include "free_abelian.h"
#include "arithmetic.h" // |Split_integer|
#include "gradings.h"
#include "kgb.h"
#include "subsystem.h"
#include "polynomials.h"
#include "K_repr.h"  // |K_type|, |K_type_pol|

namespace atlas {

namespace repr {

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
  are combined into "parameter polynomials", and this will be faster if we
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
  TorusPart y() const { return y_bits; }
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

  There is always such a representative with |lambda=rho|, and the value that
  |gamma| has for that representative represents |rho+gamma-lambda| for
  arbitrary representatives. That value used to be stored in |StandardReprMod|,
  making separate |y_bits| superfluous, and doing without the packing and
  unpacking operations for |y_bits| greatly simplifies our operations. Currently
  we have decided to not include the |rho| part, so that the stored |gamlam|
  (which is short for |gamma-lambda|) lies in the $-1$ eigenspace for $\theta$.
*/
class StandardReprMod
{
  friend class Rep_context;

  KGBElt x_part;
  RatWeight gamlam; // |real_unique(gamma-lambda)|

  StandardReprMod(KGBElt x_part, RatWeight&& gl) // private raw constructor
    : x_part(x_part), gamlam(std::move(gl.normalize())) {}

 public:
  // the raw constructor is always called through one of these two pseudo constructors:
  static StandardReprMod mod_reduce
    (const Rep_context& rc,const StandardRepr& sr);
  static StandardReprMod build
    (const Rep_context& rc, KGBElt x, RatWeight gam_lam);

  KGBElt x() const { return x_part; }
  const RatWeight& gamma_lambda() const & { return gamlam; }
  RatWeight&& gamma_lambda() && { return std::move(gamlam); }

  // pseudo constructors call |real_unique| on |gamlam|, making equality easy
  bool operator==(const StandardReprMod& other) const
  { return x_part==other.x_part and gamlam==other.gamlam; }

  typedef std::vector<StandardReprMod> Pooltype;
  bool operator!=(const StandardReprMod& another) const
  { return not operator==(another); }
  size_t hashCode(size_t modulus) const; // this one ignores $X^*$ too
}; // |class StandardReprMod|

// hashable type for |StandardReprMod| up to shift orthogonal to integral system
// This is the most reduced for of a parameter for block computation purposes
class Reduced_param
{
  KGBElt x;
  unsigned int int_sys_nr;
  unsigned int evs_reduced; // encodes |codec.internalise(gamlam) mod lattice|

  Reduced_param(KGBElt x, unsigned int i, unsigned int v)
    : x(x), int_sys_nr(i), evs_reduced(v)
  {}

public:
  Reduced_param(Reduced_param&&) = default;
  Reduced_param& operator=(Reduced_param&&) = default;

  static Reduced_param reduce // factory function, also sets |loc| for |gamma|
    (const Rep_context& rc, StandardReprMod srm, const RatWeight& gamma,
     locator& loc);
  static Reduced_param co_reduce // factory function that uses |loc| alcove
    (const Rep_context& rc, StandardReprMod srm, const locator& loc);


  using Pooltype = std::vector<Reduced_param>;
  bool operator!=(const Reduced_param& p) const
  { return x!=p.x or int_sys_nr!=p.int_sys_nr or evs_reduced!=p.evs_reduced; }
  bool operator==(const Reduced_param& p) const { return not operator!=(p); }
  size_t hashCode(size_t modulus) const; // this one ignores $X^*$ too
}; // |class Reduced_param|

/*
 A |codec| structure helps reducing |StandardReprMod| values to |Reduced_param|.
 Instances are constructed once the integral system has been determined, giving
 a |subsystem::integral_datum_item| (that can be stored in the inner class),
 whose |data| method builds a |codec| for the integral system and involution.

 Below, |in| will transform weights in a manner depending only on their coroot
 evaluations (multiplying by |coroots_matrix|) to coordinates on a basis adapted
 to $N=\Im(\theta-1)$, whose multiples by |diagonal| span $N$; this can be
 followed by reductions modulo |diagonal|, then left-multiplication by |out| to
 the lattice $N$, the result being expressed in usual coordinates

 Required properties of |codec|:
 - |internalise| is linear for arguments with the integrality for which the
   |codec| was built, and vanishes on the kernel of |coroots_matrix|
 - |Reduced_param::reduce(gamma)| depends only on remainders of the first
   |diagonal.size()| entries of |internalise(gamma)| (later entries are ignored)
 - For every $\gamma\in(1-\theta)X^*$ (the image space of |out|), one has
   |gamma == out*(internalise(gamma)/diagonal)| (coefficient-wise division)
*/
struct codec
{
  const int_Matrix coroots_matrix;
  std::vector<int> diagonal; // inv.factors for $N$ inside $-1$ eigenlattice
  int_Matrix in, out;     // see above; |in*coroots_matrix*out == diagonal|
  codec
    (const InnerClass& ic,
     InvolutionNbr inv,
     const int_Matrix& int_simp_coroots);

  // main method, compacts |gamma| (with proper integrality) to quotient element
  int_Vector internalise (const RatWeight& gamma) const;
}; // |struct codec|

// This class stores the information necessary to interpret a |StandardRepr|
class Rep_context
{
  const RootDatum& rd;
  InnerClass& ic;
  const TwistedWeylGroup& twisted_W;
  const InvolutionTable& i_tab;
  const KGB& KGB_set;
  RealReductiveGroup& G; // we could have used just this field for |Rep_context|

 public:
  explicit Rep_context(RealReductiveGroup &G);

  // accessors
  unsigned int rank() const { return root_datum().rank(); }
  const RootDatum& root_datum() const { return rd; }
  InnerClass& inner_class() const { return ic; }
  const InvolutionTable& involution_table() const { return i_tab; }
  const WeylGroup& Weyl_group() const
    { return twisted_Weyl_group().Weyl_group(); }
  const TwistedWeylGroup& twisted_Weyl_group() const { return twisted_W; }
  // const WeylGroup& cofolded_W() const { return inner_class().cofolded_W(); }
  const KGB& kgb() const { return KGB_set; }
  RealReductiveGroup& real_group() const { return G; }
  const RatCoweight& g_rho_check() const { return G.g_rho_check(); }
  RatCoweight g() const { return G.g(); }

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

  // reconstruct |StandardRepr| from |srm| and difference of |gamma_lambda|s
  StandardRepr sr (const StandardReprMod& srm,const RatWeight& gamma)  const;
  // modify (stored) |srm| by |bm|, then reconstruct |StandardRepr| at |gamma|
  StandardRepr sr
    (StandardReprMod srm, const block_modifier& bm, const RatWeight& gamma)
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
      (z,RankFlags(constants::lt_mask[root_datum().semisimple_rank()])); }
  void normalise(K_repr::K_type& z) const; // which ensures a normalised form

  simple_list<std::pair<K_repr::K_type,int> > // unsorted list
    finals_for(K_repr::K_type t) const; // defined in K_repr.cpp

  // conservative estimate for lowest height that can be obtained from |lambda|
  level height_bound(RatWeight lambda) const; // "projecting to dominant cone"
  // apart from producing a result, these methods also make |t| theta-stable:
  sl_list<K_repr::K_type> KGP_set (K_repr::K_type& t) const;
  K_repr::KT_pol K_type_formula (K_repr::K_type& t,level cutoff) const;

  K_type_poly branch(K_type_poly P, level cutoff) const;

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

  // shift continuous component of |srm| by |diff| and then renormalise it
  void shift (const RatWeight& diff, StandardReprMod& srm) const;
  StandardReprMod shifted(const RatWeight& diff, StandardReprMod srm) const
  { shift(diff,srm); return srm; } // perform |shift| on a copy and return it

  RatWeight gamma_lambda_rho(const StandardReprMod& z) const
  { return z.gamma_lambda()+rho(root_datum()); }

  bool is_standard  // whether $I(z)$ is non-virtual: gamma imaginary-dominant
    (const StandardRepr& z) const; // simply-imaginary witness
  bool is_dominant (const StandardRepr& z) const
  { return is_dominant_ratweight(root_datum(),z.gamma()); }
  bool is_nonzero  // whether $I(z)!=0$: no singular simply-imaginary compact
    (const StandardRepr& z) const; // simply-imaginary witness
  bool is_normal // whether |z==normalise(z)|: has no singular complex descents
    (const StandardRepr& z) const;
  bool is_semifinal  // whether $I(z)$ unrelated by Hecht-Schmid to more compact
    (const StandardRepr& z) const; // singular real witness
  bool is_final // dominant nonzero without singular descents: all of the above
    (const StandardRepr& z) const;
  unsigned int orientation_number(StandardRepr z) const;

  bool is_fixed(StandardRepr z, const WeightInvolution& delta) const;
  bool is_delta_fixed(const StandardRepr& z) const
  { return is_fixed(z,inner_class().distinguished()); }

  void make_dominant(StandardRepr& z) const; // ensure |z.gamma()| dominant

  // do |make_dominant|, and descend through any singular complex descents
  void deform_readjust(StandardRepr& z) const;
  // do the same, ensuring fixed choice of descent-minimum among equivalent ones
  void normalise(StandardRepr& z) const; // which ensures a normalised form

  bool equivalent(StandardRepr z0, StandardRepr z1) const; // by value

  // transform |srm| by Weyl element |w|, or its inverse
  template<bool left_to_right> void transform
    (const WeylElt& w, StandardReprMod& srm) const;

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

  SR_poly scale(const SR_poly& P, const RatNum& f) const;
  K_type_poly scale_0(const SR_poly& P) const;

  simple_list<std::pair<StandardRepr,int> >  // unsorted list
    finals_for (StandardRepr z) const; // like |finals_for| K-type (by value)
  SR_poly expand_final(StandardRepr z) const; // the same, as |SR_poly| (by value)

  Weight to_simple_shift(InvolutionNbr theta, InvolutionNbr theta_p,
			 RootNbrSet pos_to_neg) const; // |pos_to_neg| by value


  // |gamlam-srm.gamma_lambda()|, adapted by equivalence to be orth to integrals
  RatWeight make_diff_integral_orthogonal
    (const RatWeight& gamlam, const StandardReprMod& srm) const;

  // modifier for block for |srm1, loc1| relative to stored one for |srm0, loc0|
  block_modifier make_relative_to
    (const StandardReprMod& srm0, const locator& loc0,
     StandardReprMod srm1, locator&& loc1) const;

 private:
  // compute $\check\alpha\dot(1+\theta_x)\lambda$, with $(x,\lambda)$ from $t$
  int theta_plus_1_eval (const K_repr::K_type& t, RootNbr alpha) const;
  // auxiliary for |make_diff_integral_orthogonal|
  // find element in |(1-theta)X^*| with same evaluation on all integral coroots
  Weight theta_1_preimage (const RatWeight& offset, const codec& cd) const;
  RankFlags singular_simples (const StandardRepr& z) const;
  WeylElt complex_descent_w (KGBElt x, RankFlags singulars) const; // in |W|
  // make integrally dominant, with precomputed integral subsystem; return path
  WeylWord make_dominant(StandardRepr& z,const SubSystem& subsys) const;
  void complex_crosses (StandardRepr& z, const WeylElt& w) const;
  void to_singular_canonical(RankFlags gens, StandardRepr& z) const;
  level height(Weight theta_plus_1_gamma) const;

  K_repr::KT_pol monomial_product // $P*X^e$; implemented in K_repr.cpp
    (const K_repr::KT_pol& P, const Weight& e) const;
}; // |Rep_context|

/*
  In internal computations for deformation, it will be important to have a quite
  compact representation of (twisted) deformation formulas. To this end we shall
  make (inside a |Rep_table|) a table of K-types encountered (the number of
  possible K-types for a given real form is fairly limited), and represent
  linear coefficients as a map from such K-types to the ring of coefficients (in
  practice |Split_integer|). Moreover, for storing formulas as collections of
  key-value pairs, we wish to avoid using the |Free_Abelian| container, derived
  from |std::map| which requires a lot of additional memory per node stored,
  preferring to it our tailor made |Free_Abelian_light| which, being a short
  list of sorted vectors of key-value pairs, has no per-node memory overhead.
*/

using K_type_nr = unsigned int; // hashed in |Rep_table| below
using K_type_nr_poly = Free_Abelian_light<K_type_nr,Split_integer>;
using KT_nr_pol = Free_Abelian_light<K_type_nr,int>;

/*
  The class |deformation_unit| to serves to store both key and values for
  deformation formula lookup inside the class |Rep_table| below.

  A key determines a parameter up to variations in $\nu$ that do not change the
  integral part (always rounding down) of the evaluation on any positive coroot;
  the domain of such changes is called an alcove here. (Strictly speaking, the
  set of associated infinitesimal characters is the intersection of an alcove
  with an affine subspace parallel to the $-1$ eigenspace of the involution; no
  (discrete) variation is allowed in the direction of the $+1$ eigenspace.)

  Due to the subtle definition of equivalence, the methods |operator!=| and
  |hashCode| for hashing purposes are more elaborate that is typically the case
  for hashable entry types. The fact that associated values are also stored
  inside the entry type (while being ignored by the mentioned methods) is also
  not typical, but valid as a way to use |HashTable| as an association table.

  For compactness we store a sample |StandardRepr|, but the functions used for
  hashing only take into account aspects that are unchanging within the alcove.
 */
class deformation_unit
{
  friend class Rep_table; // while not essential, allows easier instrumenting

  StandardRepr sample; // first parameter found in alcove, its reference point
  // the reference point aspect only serves to determine the default extension

  KT_nr_pol lowest_K_types, def_contrib, twdef_contrib;
  K_type_nr_poly untwisted, twisted;
  const Rep_context& rc; // access coroots etc. necessary for alcove testing

  RankFlags status; // set bits (0-3) when |defcontrib|-|twisted| are defined
public:
  deformation_unit(const Rep_context& rc, const StandardRepr& sr)
    : sample(sr)
    , lowest_K_types(), def_contrib(), twdef_contrib()
    , untwisted(), twisted()
    , rc(rc), status(0)
  {}
  deformation_unit(const Rep_context& rc, StandardRepr&& sr)
  : sample(std::move(sr))
  , lowest_K_types(), def_contrib(), twdef_contrib()
  , untwisted(), twisted()
  , rc(rc), status(0)
  {}

  deformation_unit(deformation_unit&&) = default; // type is only movable


  bool has_def_contrib() const { return status.test(0); }
  bool has_twdef_contrib() const { return status.test(1); }
  bool has_deformation_formula() const  { return status.test(2); }
  bool has_twisted_deformation_formula() const  { return status.test(3); }

  const KT_nr_pol& LKTs() const { return lowest_K_types; }
  const KT_nr_pol& deformation_contribution() const
  { return def_contrib; }
  const KT_nr_pol& twisted_deformation_contribution() const
  { return twdef_contrib; }

  void set_LKTs
    (Rep_table& rt, simple_list<std::pair<K_repr::K_type,int> >&& finals);
  void set_def_contrib(KT_nr_pol&& p)
  { status.set(0); def_contrib=std::move(p); }
  void set_twdef_contrib(KT_nr_pol&& p)
  { status.set(1); twdef_contrib=std::move(p); }

  size_t def_form_size () const { return untwisted.size(); }
  size_t twisted_def_form_size () const { return twisted.size(); }

  K_type_nr_poly def_formula() const       { return untwisted.copy(); }
  K_type_nr_poly twisted_def_formula() const { return twisted.copy(); }


  K_type_nr_poly set_deformation_formula (K_type_nr_poly formula)
  { untwisted = formula.copy(); status.set(2); return std::move(formula); }

  K_type_nr_poly set_twisted_deformation_formula (K_type_nr_poly formula)
  { twisted = formula.copy(); status.set(3); return std::move(formula); }

// special members required by HashTable
  using Pooltype = std::vector<deformation_unit>;
  bool operator!=(const deformation_unit& another) const; // distinct alcoves?
  size_t hashCode(size_t modulus) const; // value depending on alcove only
}; // |class deformation_unit|

// Information about some alcove face relative to one of the fundamental alcove.
struct locator
{
  unsigned int int_syst_nr; // sequence number for integral system in inner class
  WeylElt w; // transformation, applied to a fundamental alcove integral system
  sl_list<RootNbr> simp_int; // images simply integral roots, sorted
  Permutation simple_pi; // to reorder FA integral system simple generators by
};

/*
  A block is stored for only one alcove orientation, so the lookup methods
  return a block reference with |block_modifier|. The stored parameter will
  first be shifted by |shift| then transformed by |w|. The simply integral
  roots for the modified blocks are the |w|-images of those in the stored block,
  but they are reordered increasingly, as they would be if the transformed
  block had been generated freshly; |simple_pi| is the permutation involved.
*/
struct block_modifier
{
  RatWeight shift; // initial change to the stored |gamlam| value
  WeylElt w; // initially from fundamental alcove integral system
  RootNbrList simp_int; // images of simply integral roots of block, sorted
  Permutation simple_pi; // to reorder integral system simple generators by
  block_modifier () = default;
  block_modifier(RatWeight&& shift,locator&& loc)
    : shift(std::move(shift))
    , w(std::move(loc.w))
    , simp_int(loc.simp_int.begin(),loc.simp_int.end())
    , simple_pi(std::move(loc.simple_pi))
  {}
  block_modifier (const common_block& b); // construct trivial modifier
};

struct sub_triple; // implementation specific, needed in auxiliary method type



/*
  In addition to providing methods inherited from |Rep_context|, the class
  |Rep_table| provides storage for data that was previously computed for
  various related nonzero final |StandardRepr| values.

  The data stored are:
  * The |blocks::common_block| values encountered for this real form, together
    with their KL data if computed; these are linked together in |block_list|,
    and accessed via |place|. To this end |StandardRepr| values are mapped to a
    smaller set of |Reduced_param| values that are hashed into |reduced_hash|,
    where the sequence number indexes |place| which provides both an iterator
    into the block list and an element number within the block. The latter
    number need not be unique, in other words the map to |Reduced_param| may
    fail to be injective even restricted to a single block, but block elements
    with the same image are conjugate under some block automorphism. The effect
    is that if a |StandardRepr| value comes along whose |Reduced_param| value is
    already tabled, then there is an isomorphism of its block to a stored block
    for which it corresponds to the element indicated in |place|; when we use
    this, all other block elements will lead to |StandardRepr| values deduced
    via that isomorphism. The isomorphism may involve action of $W^\delta$, and
    a |block_modifier| records the details for deducing |StandardRepr| values.
  * A hash table |alcove_hash| of |deformation_unit|s, each representing a single
    alcove, also holding an ordinary and possibly a twisted deformation formula
  * A table |reduced_hash| of distinct |Reduced_param| values encountered. As
    mentioned the reduction serves to as much as possible share equivalent block
    structures, and each element accesses a block and an element in it.
  * a table |K_type_hash| of |K_type| values having been found to occur in
    deformation formulas, which allows compact representation of the latter
  * tables |KL_poly_hash| and |poly_hash| of |kl::KLPol| (positive coefficient)
    respectively |ext_kl::Pol| (integer coefficient) polynomials, which can be
    shared among blocks to reduce the size of their tables of such polynomials
    (they may choose to instead maintain their local polynomials themselves).
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

  using located_block = std::pair<blocks::common_block,locator>;
  sl_list<located_block> block_list;

  using bl_it = sl_list<located_block>::iterator;
  std::vector<std::pair<bl_it, BlockElt> > place; // parallel to |reduced_pool|

 public:
  Rep_table(RealReductiveGroup &G);
  ~Rep_table();
  // both defined out of line because of implicit use |common_block| destructor

  ext_KL_hash_Table* shared_poly_table () { return &poly_hash; }

  // the |length| method generates a partial block, for best amortised efficiency
  unsigned short length(StandardRepr z); // by value

  // look up block which must be full; if not found construct and swallow crumbs
  blocks::common_block& lookup_full_block
    (StandardRepr& sr, // input argument, but by reference: it will be normalised
     BlockElt& z, block_modifier& bm); // output arguments

  // look up block containing |sr|; if not found constuct one and swallow crumbs
  blocks::common_block& lookup // construct only partial block if necessary
    (StandardRepr& sr, // input argument, but by reference: it will be normalised
     BlockElt& z, block_modifier& bm); // output arguments

  // polynomial representing column of |z| in KL-table, evaluated at $q:=s$
  SR_poly KL_column_at_s(StandardRepr z); // by value
  // column of |z| in KL-table, as (sparse) list of pairs $(x,P_{x,y})$
  simple_list<std::pair<BlockElt,kl::KLPol> >
    KL_column(common_block& block, BlockElt y);
  // length-alternating sum of twisted KL polynomials at $q:=s$ (results tabled)
  SR_poly twisted_KL_column_at_s(StandardRepr z); // by value
  // heights-bounded KL column; computed via dual block, no tabling is done
  SR_poly KL_column_at_s_to_height (StandardRepr p, level height_bound);

#if 0 // this would now require an actual |delta|-fixed |gamma| to be supplied
  size_t find_reduced_hash(const StandardReprMod& srm) const
  { block_modifier bm;
    return reduced_hash.find(Reduced_param::reduce(*this,srm,bm)); }
#endif

  K_type_nr match(K_type&& t) { return K_type_hash.match(std::move(t)); }
  K_type stored_K_type(K_type_nr i) const { return K_type_pool[i].copy(); }

  // a signed multiset of final parameters needed to be considered (i.e., their
  // deformations to $\nu=0$ included) when deforming |y| a bit towards $\nu=0$
  sl_list<std::pair<StandardRepr,int> > deformation_terms
    (blocks::common_block& block, BlockElt y,
     const block_modifier& bm, const RatWeight& gamma);

/*
   A method intended to be an alternative to |deformation|, but ignoring any
   contributions with height above |heign_bound|. In order to do so, the
   recursive "depth first" approach of |deformation| is relaced by a "breadth
   first" approch, that works out all deformation associated with the current
   block at once, so that the resulting terms can all be deformed towards
   $\nu=0$ until (maybe) hitting a reducibility determined by another block.

   Implementation is using dual block, like |KL_column_at_s_to_height|, and
   obtained results are not tabled. Might well be slower than height-unlimited.

   Compute the signed multi-set of final parameters "post deformation"
   (subsequently they will be deformed towards zero without reduction here)
   obtained from the elements of the block of |p|, with their "pre deformation"
   coefficients taken (and removed) from |queue|, where all block elements of
   height strictly above |height_bound| are ignored. The result is a list of
   pairs of a |StandardRepr| with its |Split_integer| multiplicity.
*/
  sl_list<SR_poly::value_type> block_deformation_to_height
    (StandardRepr p, SR_poly& queue, level height_bound); // |p| by value


  // full deformation to $\nu=0$ of |z|
  const deformation_unit& deformation(StandardRepr z); // by value
  K_type_nr_poly full_deformation(StandardRepr z); // by value

  // like |deformation_terms|; caller multiplies returned coefficients by $1-s$
  sl_list<std::pair<StandardRepr,int> > twisted_deformation_terms
    (blocks::common_block& block, ext_block::ext_block& eblock,
     BlockElt y, // in numbering of |block|, not |eblock|
     RankFlags singular, const block_modifier& bm, const RatWeight& gamma);
#if 0
  SR_poly twisted_deformation_terms (unsigned long sr_hash) const;
  // once a parameter has been entered, we can compute this without a block
#endif

  sl_list<StandardReprMod> Bruhat_below
  (const common_context& ctxt, const StandardReprMod& init) const;

  blocks::common_block& add_block_below // partial; defined in common_blocks.cpp
    (const StandardReprMod& srm, BitMap* subset, const locator& loc);

  // full twisted deformation, with |flip| telling whether to multiply by |s|
  const deformation_unit& twisted_deformation(StandardRepr z); // by value
  K_type_nr_poly
  twisted_full_deformation(StandardRepr z, bool& flip); // by value

 private:
  class Bruhat_generator; // helper class: internal |add_block_below| recursion
  void block_erase (bl_it pos); // erase from |block_list| in safe manner
  // add a full block, assuming has trivial attitude
  void add_block (const StandardReprMod& srm, const locator& loc);
  void append_block_containing // maybe appends an item to |sub_blocks|
    (const StandardReprMod& elt, size_t place_limit, locator block_loc,
     sl_list<sub_triple>& sub_blocks);
  void swallow_then_append_singleton // import and erase |subs|, then add |tail|
    (const sl_list<sub_triple>& subs, const locator& loc,
     sl_list<located_block>&& tail);
}; // |Rep_table|


// another extension of |Rep_context|, used by both |common_block| constructors
// given |gamma| or just its facet, pre-compute integral system for common block
class common_context
{
  const Rep_context& rep_con;
  sl_list<RootNbr> simply_int; // simply integral roots for facet, increasing
  const SubSystem sub; // embeds |simply_int| subsystem into full root datum
public:
  // constructor for all normal uses:
  common_context (const Rep_context& rc, const RatWeight& gamma);

  // special constructor for |add_block| and |add_block_below| only:
  common_context (const Rep_context& rc, const locator& loc);

  // accessors
  const Rep_context& rc() const { return rep_con; }
  const KGB& kgb() const { return rep_con.kgb(); }
  const InvolutionTable& involution_table() const
    { return rep_con.involution_table(); }
  const RootDatum& full_root_datum() const { return rep_con.root_datum(); }
  const SubSystem& subsys() const { return sub; }
  RootNbrList simply_integrals() const { return simply_int.to_vector(); }

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
  and are unchanged under integral changes to |gamma|. Data fields that are
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

K_type_poly export_K_type_pol(const Rep_table& rt,const K_type_nr_poly& P);

} // |namespace repr|

} // |namespace atlas|

#endif
