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

#include "innerclass.h" // inlines
#include "realredgp.h"	// inlines

#include "hashtable.h"
#include "free_abelian.h"
#include "arithmetic.h" // |SplitInteger|

#include "rootdata.h" // for |rho|, so |Rep_context::lambda| can be inlined

namespace atlas {

namespace repr {

/*
  A parameter of a standard representation is determined by a triplet
  $(x,\tilde\lambda,gamma)$, where $x\in K\\backslash G/B$ for our fixed real
  form (determining amongs others an involution $\thata$ of $X^*$),
  $\tilde\lambda$ is a genuine character of the $\rho$-cover of
  $H^{\theta_x}$, and $\gamma$ is a character of the complex Lie algebra $h$.
  The latter two values are related; $(1+\theta)\gamma=(1+\theta)\lambda$, in
  other words $\gamma-\tilde\lambda$ is fixed by $-\theta$; the projection of
  $\gamma$ on the $+1$-eigenspace of $\theta$ is determined by this relation
  and is called the free part $\lambda_0$ of $\gamma$. The difference
  $\gamma-\lambda_0$, i.e., the projection of $\gamma$ on the $-1$-eigenspace,
  is called $\nu$. This component is what we are adding with respect to the
  values encoded in the |standardrepk::StandardRepK| type. The part of
  $\tilde\lambda$ that is independent of $\lambda_0$ is its "torsion part"
  (due to the disconnectedness of $H(R)_c$), which would be represented in the
  |Block| structure by the |TorusPart| component of the |TitsElt| of the dual
  KGB-element ($y$). In fact we also convert it internally to a |TorusPart|
  here, called |y_bits|. It represents an element of the quotient of
  $(X^*)^{-\theta}$ by the image $(1-\theta)X^*$; from it we can recover
  $\tilde\lambda-\lambda_0$, using |involutions::InvolutionTable::unpack|.

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


// This class stores the information necessary to interpret a |StandardRepr|
class Rep_context
{
  RealReductiveGroup& G;
  const KGB& KGB_set;

 public:
  explicit Rep_context(RealReductiveGroup &G);

  // accessors
  RealReductiveGroup& realGroup() const { return G; }
  const InnerClass& innerClass() const { return G.innerClass(); }
  const RootDatum& rootDatum() const { return G.rootDatum(); }
  const WeylGroup& weylGroup() const { return G.weylGroup(); }
  const TwistedWeylGroup& twistedWeylGroup() const
    { return G.twistedWeylGroup(); }
  const TitsGroup& titsGroup() const { return G.titsGroup(); }
  const TitsCoset& basedTitsGroup() const { return G.basedTitsGroup(); }
  const KGB& kgb() const { return KGB_set; }
  size_t rank() const;

  const TwistedInvolution involution_of_Cartan(size_t cn) const;

  RatWeight gamma // compute (representative of) infinitesimal character
    (KGBElt x, const Weight& lambda_rho, const RatWeight& nu) const;
  StandardRepr sr_gamma // use this one when infinitesimal character is known
    (KGBElt x, const Weight& lambda_rho, const RatWeight& gamma) const;
  StandardRepr sr // construct parameter from |(x,\lambda,\nu)| triplet
    (KGBElt x, const Weight& lambda_rho, const RatWeight& nu) const
    { return sr_gamma(x,lambda_rho,gamma(x,lambda_rho,nu)); }

  StandardRepr
    sr(const standardrepk::StandardRepK& srk,
       const standardrepk::SRK_context& srkc,
       const RatWeight& nu) const;

  // component extraction
  const WeightInvolution& theta (const StandardRepr& z) const;

  Weight lambda_rho(const StandardRepr& z) const;
  RatWeight lambda(const StandardRepr& z) const // half-integer
  { return rho(rootDatum()).normalize()+lambda_rho(z); }
  RatWeight gamma_0 // infinitesimal character deformed to $\nu=0$
    (const StandardRepr& z) const;

  RatWeight nu(const StandardRepr& z) const; // rational, $-\theta$-fixed

  // the value of $\exp_{-1}(\gamma-\lambda)$ is $y$ value in a |param_block|
  TorusElement y_as_torus_elt(const StandardRepr& z) const;

  // attributes; they set |witness| only in case they return |false|
  bool is_standard  // whether $I(z)$ is non-virtual: gamma imaginary-dominant
    (const StandardRepr& z, RootNbr& witness) const; // simply-imaginary witness
  bool is_dominant  // whether |gamma| is dominant
    (const StandardRepr& z, RootNbr& witness) const; // simple witness
  bool is_nonzero  // whether $I(z)!=0$: no singular simply-imaginary compact
    (const StandardRepr& z, RootNbr& witness) const; // simply-imaginary witness
  bool is_normal // wither |z==normal(z)|; implies no singular complex descents
    (const StandardRepr& z) const; // complex simple witness
  bool is_semifinal  // whether $I(z)$ unrelated by Hecht-Schmid to more compact
    (const StandardRepr& z, RootNbr& witness) const; // singular real witness
  bool is_final // dominant nonzero without singular descents: all of the above
    (const StandardRepr& z) const;
  bool is_oriented(const StandardRepr& z, RootNbr alpha) const;
  unsigned int orientation_number(const StandardRepr& z) const;

  bool is_twist_fixed(StandardRepr z, const WeightInvolution& delta) const;
  bool is_twist_fixed(const StandardRepr& z) const
  { return is_twist_fixed(z,innerClass().distinguished()); }

  void make_dominant(StandardRepr& z) const; // ensure |z.gamma()| dominant

  // in addition to |make_dominant| ensure a normalised form of the parameter
  void normalise(StandardRepr& z) const;

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

  typedef Free_Abelian<StandardRepr,Split_integer,compare> poly;

  poly scale(const poly& P, const Rational& f) const;
  poly scale_0(const poly& P) const;

  containers::sl_list<StandardRepr> finals_for // like |param_block::finals_for|
    (StandardRepr z) const; // by value
  poly expand_final(StandardRepr z) const; // the same, as |poly| (by value)

  std::ostream& print (std::ostream&,const StandardRepr& z) const;
  std::ostream& print (std::ostream&,const poly& P) const;

 private:
  // make integrally dominant, with precomputed integral subsystem; return path
  WeylWord make_dominant(StandardRepr& z,const SubSystem& subsys) const;
  StandardRepr& singular_cross (StandardRepr& z,weyl::Generator s) const;
  void to_singular_canonical(RankFlags gens, StandardRepr& z) const;
  unsigned int height(Weight theta_plus_1_gamma) const;
}; // |Rep_context|

typedef Rep_context::poly SR_poly;

/*
  In addition to providing methods inherited from |Rep_context|, the class
  |Rep_table| provides storage for data that was previously computed for
  various related nonzero final |StandardRepr| values.

  The data stored consists of lengths, (twisted) KL polynomials evaluated as
  $q=s$, and (twisted) full deformation formulae. This class provides methods
  for those computations, and handles the data storage and retrieval.

  The deformation information for a parameter will actually be stored for its
  first reducibility point (thus avoiding some duplication of the same
  information, and avoiding memory usage for parameters without any
  reducibility points), at the normalised form of the parameter there (the
  deformation might have made some complex simple coroots singular, crossing
  their wall). Since no slots are created for parameters that are zero or
  non-final, such parameters must be expanded into finals before attempting
  look-up or storage of deformation formulae; the nonzero and final conditions
  are preserved at intermediate deformation points (though possibly not at the
  terminating point $\nu=0$ of the deformation).
*/
class Rep_table : public Rep_context
{
  std::vector<StandardRepr> pool;
  HashTable<StandardRepr,unsigned long> hash;
  std::vector<unsigned short int> lengths;
  std::vector<SR_poly> KLV_list; // indexed by |hash| values for |StandardRepr|s
  std::vector<SR_poly> def_formula; // idem

  std::vector<SR_poly> twisted_KLV_list; // values at twist-fixed |hash|s only
  std::vector<SR_poly> twisted_def_formula; // idem

 public:
  Rep_table(RealReductiveGroup &G)
    : Rep_context(G), pool(), hash(pool), KLV_list(), def_formula()
  {}

  unsigned int length(StandardRepr z); // by value

  SR_poly KL_column_at_s(StandardRepr z); // by value
  SR_poly twisted_KL_column_at_s(StandardRepr z); // by value

  SR_poly deformation_terms (param_block& block,BlockElt entry_elem);
  // here |block| is non-|const| because it calls |add_block|

  SR_poly deformation(const StandardRepr& z);

  SR_poly twisted_deformation_terms (param_block& block,BlockElt entry_elem);
  // here |block| is non-|const| because it calls |add_block|

  SR_poly twisted_deformation(StandardRepr z); // by value

 private:
  void add_block(param_block& block,
		 containers::sl_list<BlockElt>& extended_finals);
  // here |block| is non-|const| as the method generates KL polynomials in it
  // and |survivors| is non-|const| because the method computes and exports it

  void add_block(ext_block::ext_block& block, param_block& parent,
		 BlockElt top_elt,
		 containers::sl_list<BlockElt>& extended_finals);
  // here |block| is non-|const|; the method generates twisted KLV polys in it

}; // |Rep_table|


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
