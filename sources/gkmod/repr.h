/*!
\file
\brief Classes and utilities for manipulating representations
*/
/*
  This is repr.h

  Copyright (C) 2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef REPR_H  /* guard against multiple inclusions */
#define REPR_H

#include "standardrepk.h"
#include "latticetypes.h"
#include "tits.h"
#include "kgb.h"
#include "realredgp.h"

namespace atlas {

namespace repr {

/*
We represent the parameter of a standard representation as a triplet
$(x,\lambda,gamma)$, where |x| is an element of the set $K\\backslash G/B$ for
our fixed real form, $\lambda$ is a character $\lambda$ of (the $\rho$-cover
of) $H^{\theta_x}$, and $\gamma$ is a character of the complex Lie algebra $h$.
The latter two values are interrelated by $(1+\theta)\gamma=(1+\theta)\lambda$;
the projection of $\gamma$ on the $+1$-eigenspace of $\theta_x$ is determined
by this relation and is called the discrete part of $\gamma$. The difference
with the discrete part, i.e., the projection of $\gamma$ on the
$-1$-eigenspace, is called $\nu$, this is what we are adding with respect to
the values encoded in |standardrepk::StandarRepK| values.

Although $\gamma$ could in principle take any complex values compatible with
$\lambda$, we shall only be interested in real values, and in fact record a
rational value because all interesting phenomena take place at rational points.
*/
class StandardRepr
{
  friend class Rep_context;

  kgb::KGBElt x;
  latticetypes::LatticeElt lambda_rho; // $\lambda-\rho$
  latticetypes::RatWeight infinitesimal_char; // $\gamma$

 public:
  StandardRepr (kgb::KGBElt xx,
		const latticetypes::LatticeElt& lr,
		const latticetypes::RatWeight& gamma)
    : x(xx), lambda_rho(lr), infinitesimal_char(gamma) {}

  const latticetypes::RatWeight& gamma() const { return infinitesimal_char; }

  bool operator== (const StandardRepr&) const;
// special members required by hashtable::HashTable

  typedef std::vector<StandardRepr> Pooltype;
  bool operator!=(const StandardRepr& another) const
    { return not operator==(another); }
  size_t hashCode(size_t modulus) const;



}; // |class StandardRepr|


// This class stores the information necessary to interpret a |StandardRepr|
class Rep_context
{
  realredgp::RealReductiveGroup& G;
  const kgb::KGB_base& KGB_set;

 public:
  Rep_context(realredgp::RealReductiveGroup &G);

  // accessors
  complexredgp::ComplexReductiveGroup& complexGroup() const
  { return G.complexGroup(); }
  const rootdata::RootDatum& rootDatum() const { return G.rootDatum(); }
  const weyl::WeylGroup& weylGroup() const { return G.weylGroup(); }
  const weyl::TwistedWeylGroup& twistedWeylGroup() const
    { return G.twistedWeylGroup(); }
  const tits::TitsGroup& titsGroup() const { return G.titsGroup(); }
  const tits::BasedTitsGroup& basedTitsGroup() const
    { return G.basedTitsGroup(); }

  const weyl::TwistedInvolution twistedInvolution(size_t cn) const
    { return complexGroup().twistedInvolution(cn); }

  StandardRepr
    sr(const standardrepk::StandardRepK& srk,
       const standardrepk::KhatContext& khc,
       const latticetypes::RatWeight& nu) const;

  latticetypes::RatWeight lambda(const StandardRepr& rep) const; // half-integer
  tits::GlobalTitsElement y(const StandardRepr& rep) const;

  kgb::global_KGB dual_KGB(const StandardRepr& rep) const;
}; // |Rep_context|


} // |namespace repr|

} // |namespace atlas|

#endif
