/*
  This is repr.cpp

  Copyright (C) 2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "repr.h"
#include "standardrepk.h"
#include "latticetypes.h"
#include "tits.h"
#include "blocks.h"

namespace atlas {
  namespace repr {

Rep_context::Rep_context(realredgp::RealReductiveGroup &G_R)
  : G(G_R), KGB_set(G_R.kgb())
{}

StandardRepr
  Rep_context::sr
    (const standardrepk::StandardRepK& srk,
     const standardrepk::KhatContext& khc,
     const latticetypes::RatWeight& nu) const
{
  tits::TitsElt a = khc.titsElt(srk); // was reduced during construction |srk|
  kgb::KGBElt x= khc.kgb().lookup(a,titsGroup());
  latticetypes::Weight lambda2 = khc.lift(srk); // doubled coordinates
  latticetypes::RatWeight lambda(lambda2,2);
  latticetypes::LatticeMatrix theta = G.cartan(srk.Cartan()).involution();
  latticetypes::RatWeight diff = lambda - nu;
  latticetypes::RatWeight theta_diff(theta.apply(diff.numerator()),
				     diff.denominator()); // theta(lambda-nu)
  latticetypes::Weight lambda_rho = (lambda2-khc.rootDatum().twoRho())/=2;
  return StandardRepr(x,lambda_rho,((lambda+nu+theta_diff)/=2).normalize());
}

// return $\lambda \in \rho+X^*$ as half-integer rational vector
latticetypes::RatWeight Rep_context::lambda(const StandardRepr& rep) const
{
  latticetypes::LatticeElt num = rep.lambda_rho * 2 + rootDatum().twoRho();
  return latticetypes::RatWeight(num,2);
}

// convert |lambda| and |gamma| into $y$ value
// Formula: $\exp(i\pi(\gamma-\lambda)) \sigma_{tw} \delta_1$
tits::GlobalTitsElement Rep_context::y(const StandardRepr& rep) const
{
  tits::TorusElement t((rep.gamma()-lambda(rep))/=2);
  const weyl::TwistedWeylGroup& W = KGB_set.twistedWeylGroup();
  const weyl::TwistedWeylGroup dual_W (W,tags::DualTag());
  weyl::TwistedInvolution tw =
    blocks::dual_involution(KGB_set.involution(rep.x),W,dual_W);
  return tits::GlobalTitsElement(t,tw);
}

  } // |namespace repr|
} // |namespace atlas|
