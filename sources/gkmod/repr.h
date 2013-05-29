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

#include "atlas_types.h"

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
$(x,\lambda,gamma)$, where |x| is an element of the set $K\\backslash G/B$ for
our fixed real form, $\lambda$ is a character $\lambda$ of (the $\rho$-cover
of) $H^{\theta_x}$, and $\gamma$ is a character of the complex Lie algebra
$h$. The latter two values are related; $(1+\theta)\gamma=(1+\theta)\lambda$;
the projection of $\gamma$ on the $+1$-eigenspace of $\theta_x$ is determined
by this relation and is called the discrete part of $\gamma$. The difference
with the discrete part, i.e., the projection of $\gamma$ on the
$-1$-eigenspace, is called $\nu$, this is what we are adding with respect to
the values encoded in |standardrepk::StandarRepK| values. The part of
$\lambda$ that is independent of the discrete part of $\gamma$ is its "torsion
part" (disconnected $H(R)_c$), which would be represented in the |Block|
structure by the |TorusPart| component of the |TitsElt| of the dual
KGB-element ($y$). In fact we convert it intenally to a |TorusPart| here too.

Although $\gamma$ could in principle take any complex values compatible with
$\lambda$, we shall only be interested in real values, and in fact record a
rational value because all interesting phenomena take place at rational points.
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
  const ComplexReductiveGroup& complexGroup() const { return G.complexGroup(); }
  const RootDatum& rootDatum() const { return G.rootDatum(); }
  const WeylGroup& weylGroup() const { return G.weylGroup(); }
  const TwistedWeylGroup& twistedWeylGroup() const
    { return G.twistedWeylGroup(); }
  const TitsGroup& titsGroup() const { return G.titsGroup(); }
  const TitsCoset& basedTitsGroup() const { return G.basedTitsGroup(); }
  const KGB& kgb() const { return KGB_set; }
  size_t rank() const;

  const TwistedInvolution twistedInvolution(size_t cn) const;

  StandardRepr
    sr(const standardrepk::StandardRepK& srk,
       const standardrepk::KhatContext& khc,
       const RatWeight& nu) const;

  StandardRepr
    sr(KGBElt x,
       const Weight lambda_rho,
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

  // prepare for |deform|: make |gamma| dominant, and as theta-stable as can be
  void make_dominant(StandardRepr& z) const;

  RationalList reducibility_points(const StandardRepr& z) const;

  class compare
  { Coweight level_vec; // linear form to apply to |gamma| for ordering
  public:
    compare (const Coweight& lv) : level_vec(lv) {}

    bool operator()(const StandardRepr& r,const StandardRepr& s) const;
  };

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
  std::vector<SR_poly> KL_list;
  std::vector<SR_poly> def_formula;

 public:
  Rep_table(RealReductiveGroup &G)
    : Rep_context(G), pool(), hash(pool), KL_list(), def_formula()
  {}

  SR_poly KL_column_at_s(StandardRepr z); // by value

  SR_poly deformation_terms (param_block& block,BlockElt entry_elem);
  // here |block| is non-|const| because it calls |add_block|

  SR_poly deformation(const StandardRepr& z);

 private:
  void add_block(param_block& block, const BlockEltList& survivors);
  // here |block| is non-|const| as the method generates KL polynomials in it

}; // |Rep_table|


// 				Functions


} // |namespace repr|

} // |namespace atlas|

#endif
