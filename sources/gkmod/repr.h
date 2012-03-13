/*! \file \brief Classes and utilities for manipulating representations */
/*
  This is repr.h

  Copyright (C) 2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef REPR_H  /* guard against multiple inclusions */
#define REPR_H

#include "atlas_types.h"

#include "matrix.h"	// containment
#include "ratvec.h"	// containment

#include "realredgp.h"	// inlines

namespace atlas {

namespace repr {

/*
We represent the parameter of a standard representation as a triplet
$(x,\lambda,gamma)$, where |x| is an element of the set $K\\backslash G/B$ for
our fixed real form, $\lambda$ is a character $\lambda$ of (the $\rho$-cover
of) $H^{\theta_x}$, and $\gamma$ is a character of the complex Lie algebra
$h$. The latter two values are interrelated by
$(1+\theta)\gamma=(1+\theta)\lambda$; the projection of $\gamma$ on the
$+1$-eigenspace of $\theta_x$ is determined by this relation and is called the
discrete part of $\gamma$. The difference with the discrete part, i.e., the
projection of $\gamma$ on the $-1$-eigenspace, is called $\nu$, this is what
we are adding with respect to the values encoded in
|standardrepk::StandarRepK| values. The part of $\lambda$ that is independent
of the discrete part of $\gamma$ is its "torsion part" (disconnected
$H(R)_c$), which would be represented in the |Block| structure by the
|TorusPart| component of the |TitsElt| of the dual KGB-element ($y$).

Although $\gamma$ could in principle take any complex values compatible with
$\lambda$, we shall only be interested in real values, and in fact record a
rational value because all interesting phenomena take place at rational points.
*/
class StandardRepr
{
  friend class Rep_context;

  KGBElt x_part;
  Weight lambda_rho; // $\lambda-\rho$
  RatWeight infinitesimal_char; // $\gamma$

 public:
  StandardRepr (KGBElt x,
		const Weight& lr,
		const RatWeight& gamma)
    : x_part(x), lambda_rho(lr), infinitesimal_char(gamma) {}

  const RatWeight& gamma() const { return infinitesimal_char; }
  KGBElt x() const { return x_part; }

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
  RealReductiveGroup& G;
  const KGB_base& KGB_set;

 public:
  Rep_context(RealReductiveGroup &G);

  // accessors
  ComplexReductiveGroup& complexGroup() const { return G.complexGroup(); }
  const RootDatum& rootDatum() const { return G.rootDatum(); }
  const WeylGroup& weylGroup() const { return G.weylGroup(); }
  const TwistedWeylGroup& twistedWeylGroup() const
    { return G.twistedWeylGroup(); }
  const TitsGroup& titsGroup() const { return G.titsGroup(); }
  const TitsCoset& basedTitsGroup() const { return G.basedTitsGroup(); }
  const KGB_base& kgb() const { return KGB_set; }

  const TwistedInvolution twistedInvolution(size_t cn) const;

  StandardRepr
    sr(const standardrepk::StandardRepK& srk,
       const standardrepk::KhatContext& khc,
       const RatWeight& nu) const;

  StandardRepr
    sr(KGBElt x,
       const Weight lambda_rho,
       const RatWeight& nu) const;

  RatWeight lambda(const StandardRepr& rep) const; // half-integer
  RatWeight nu(const StandardRepr& rep) const; // rational, $-\theta$-fixed
  GlobalTitsElement y(const StandardRepr& rep) const;

  // intended to generate global KGB-set from |y| value; not (yet) implemented
  kgb::global_KGB dual_KGB(const StandardRepr& rep) const;

  bool is_final(const StandardRepr& rep);
  bool is_oriented(const StandardRepr& rep, RootNbr alpha);
  unsigned int orientation_number(const StandardRepr& rep);

}; // |Rep_context|


} // |namespace repr|

} // |namespace atlas|

#endif
