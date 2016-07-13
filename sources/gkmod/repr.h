/*! \file \brief Classes and utilities for manipulating representations */
/*
  This is repr.h

  Copyright (C) 2009-2012 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef REPR_H  /* guard against multiple inclusions */
#define REPR_H

#include <iostream>

#include "../Atlas.h"

#include "matrix.h"	// containment
#include "ratvec.h"	// containment

#include "realredgp.h"	// inlines

#include "hashtable.h"
#include "free_abelian.h"
#include "arithmetic.h" // |SplitInteger|

namespace atlas {

namespace repr {

/*
We represent the parameter of a standard representation as a triplet
$(x,\tilde\lambda,gamma)$, where |x| is an element of the set $K\\backslash
G/B$ for our fixed real form (determining amongs others an involution $\thata$
of $X^*$), $\tilde\lambda$ is a genuine character of the $\rho$-cover of
$H^{\theta_x}$, and $\gamma$ is a character of the complex Lie algebra $h$.
The latter two values are related; $(1+\theta)\gamma=(1+\theta)\lambda$, in
other words $\gamma-\tilde\lambda$ is fixed by $-\theta$; the projection of
$\gamma$ on the $+1$-eigenspace of $\theta$ is determined by this relation and
is called the discrete part $\lamda_0$ of $\gamma$. The difference
$\gamma-\lambda_0$, i.e., the projection of $\gamma$ on the $-1$-eigenspace,
is called $\nu$. This component is what we are adding with respect to the
values encoded in the |standardrepk::StandarRepK| type. The part of
$\tilde\lambda$ that is independent of $lambda_0$ is its "torsion part"
(disconnected $H(R)_c$), which would be represented in the |Block| structure
by the |TorusPart| component of the |TitsElt| of the dual KGB-element ($y$).
In fact we convert it intenally to a |TorusPart| here too. It repesents an
element of the quotient of $X^* /2X^*$ by the image of $(X^*)^\theta$, which
can be converted to the difference $\tilde\lambda-\lambda_0$ by the method
|involutions::InvolutionTable::unpack|.

In principle $\gamma$ could take any complex values compatible with
$\tildelambda$, but we shall only be interested in real values, and in fact
record a rational value, because that is all we can do in an exact manner, and
all interesting phenomena take place at rational infinitesimal character.
*/

class StandardRepr
{
  friend class Rep_context;

  KGBElt x_part;
  TorusPart y_bits;
  RatWeight infinitesimal_char; // $\gamma$

  // one should call constructor from |Rep_context| only
  StandardRepr (KGBElt x,TorusPart y,const RatWeight& gamma)
    : x_part(x), y_bits(y), infinitesimal_char(gamma)
  { infinitesimal_char.normalize(); } // ensure this class invariant

 public:

  const RatWeight& gamma() const { return infinitesimal_char; }
  KGBElt x() const { return x_part; }
  const TorusPart& y() const { return y_bits; }

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
  StandardRepr sr(const param_block& b, BlockElt i) const;

  // component extraction
  Weight lambda_rho(const StandardRepr& z) const;

  RatWeight lambda(const StandardRepr& z) const; // half-integer
  RatWeight nu(const StandardRepr& z) const; // rational, $-\theta$-fixed

  // attributes
  bool is_standard  // whether $I(z)$ is non-virtual: gamma imaginary-dominant
    (const StandardRepr& z, RootNbr& witness) const; // simple-imaginary witness
  bool is_zero  // whether $I(z)=0$: exists singular simple-imaginary compact
    (const StandardRepr& z, RootNbr& witness) const; // simple-imaginary witness
  bool is_final  // whether $I(z)$ unrelated by Hecht-Schmid to more compact
    (const StandardRepr& z, RootNbr& witness) const; // singular real witness
  bool is_oriented(const StandardRepr& z, RootNbr alpha) const;
  unsigned int orientation_number(const StandardRepr& z) const;

  // action by equivalence of parameters (not the cross action), changing gamma
  void W_act(const WeylWord& w,StandardRepr& z) const;

  // same, but interpreting the Weyl word in a subsystem
  void W_act(const WeylWord& w,StandardRepr& z,const SubSystem& subsys) const;

  // prepare for |deform|: make |gamma| dominant, and as theta-stable as can be
  // return the sequence of Weyl generators that was applied (to the right)
  WeylWord make_dominant(StandardRepr& z) const;

  // make integrally dominant, with precomputed integral subsystem
  WeylWord make_dominant(StandardRepr& z,const SubSystem& subsys) const;

  RationalList reducibility_points(const StandardRepr& z) const; // normalised

  // the following take |z| by value, modifying and in some cases returning it
  StandardRepr cross(weyl::Generator s, StandardRepr z) const;
  StandardRepr Cayley(weyl::Generator s, StandardRepr z) const;
  StandardRepr inv_Cayley(weyl::Generator s, StandardRepr z) const;
  StandardRepr twist(StandardRepr z) const;

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

  poly expand_final(StandardRepr z) const; // express in final SReprs (by value)

  std::ostream& print (std::ostream&,const StandardRepr& z) const;
  std::ostream& print (std::ostream&,const poly& P) const;

}; // |Rep_context|

typedef Rep_context::poly SR_poly;

class Rep_table : public Rep_context
{
  std::vector<StandardRepr> pool;
  HashTable<StandardRepr,unsigned long> hash;
  std::vector<unsigned short int> lengths;
  std::vector<SR_poly> KL_list; // indexed by |hash| values for |StandardRepr|s
  std::vector<SR_poly> def_formula; // idem

 public:
  Rep_table(RealReductiveGroup &G)
    : Rep_context(G), pool(), hash(pool), KL_list(), def_formula()
  {}

  unsigned int length(StandardRepr z); // by value

  SR_poly KL_column_at_s(StandardRepr z); // by value

  SR_poly deformation_terms (param_block& block,BlockElt entry_elem);
  // here |block| is non-|const| because it calls |add_block|

  SR_poly deformation(const StandardRepr& z);

 private:
  void add_block(param_block& block, BlockEltList& survivors);
  // here |block| is non-|const| as the method generates KL polynomials in it
  // and |survivors| is non-|const| because the method computes and exports it

}; // |Rep_table|


// 				Functions


} // |namespace repr|

} // |namespace atlas|

#endif
