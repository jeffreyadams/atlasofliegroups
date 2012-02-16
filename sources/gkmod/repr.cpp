/*
  This is repr.cpp

  Copyright (C) 2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "repr.h"

#include "arithmetic.h"
#include "tits.h"

#include "kgb.h"	// various methods
#include "blocks.h"	// |dual_involution|
#include "standardrepk.h"// |KhatContext| methods

namespace atlas {
  namespace repr {

Rep_context::Rep_context(RealReductiveGroup &G_R)
  : G(G_R), KGB_set(G_R.kgb())
{}

const TwistedInvolution Rep_context::twistedInvolution(size_t cn) const
{ return complexGroup().twistedInvolution(cn); }

StandardRepr
  Rep_context::sr
    (const standardrepk::StandardRepK& srk,
     const standardrepk::KhatContext& khc,
     const RatWeight& nu) const
{
  TitsElt a = khc.titsElt(srk); // was reduced during construction |srk|
  KGBElt x= khc.kgb().lookup(a,titsGroup());
  Weight lambda2 = khc.lift(srk); // doubled coordinates
  RatWeight lambda(lambda2,2);
  WeightInvolution theta = G.cartan(srk.Cartan()).involution();
  RatWeight diff = lambda - nu;
  RatWeight theta_diff(theta*diff.numerator(),
				     diff.denominator()); // theta(lambda-nu)
  Weight lambda_rho = (lambda2-khc.rootDatum().twoRho())/=2;
  return StandardRepr(x,lambda_rho,((lambda+nu+theta_diff)/=2).normalize());
}

StandardRepr Rep_context::sr
  (KGBElt x, const Weight lambda_rho, const RatWeight& nu) const
{
  RatWeight lambda(lambda_rho*2+rootDatum().twoRho(),2);
  WeightInvolution theta =
    complexGroup().involutionMatrix(kgb().involution(x));
  RatWeight diff = lambda - nu;
  RatWeight theta_diff(theta*diff.numerator(),
				     diff.denominator()); // theta(lambda-nu)
  return StandardRepr(x,lambda_rho,((lambda+nu+theta_diff)/=2).normalize());
}

// return $\lambda \in \rho+X^*$ as half-integer rational vector
RatWeight Rep_context::lambda(const StandardRepr& rep) const
{
  Weight num = rep.lambda_rho * 2 + rootDatum().twoRho();
  return RatWeight(num,2);
}

// convert |lambda| and |gamma| into $y$ value
// Formula: $\exp(i\pi(\gamma-\lambda)) \sigma_{tw} \delta_1$
GlobalTitsElement Rep_context::y(const StandardRepr& rep) const
{
  TorusElement t = y_values::exp_pi(rep.gamma()-lambda(rep));
  const TwistedWeylGroup& W = KGB_set.twistedWeylGroup();
  const TwistedWeylGroup dual_W (W,tags::DualTag());
  TwistedInvolution tw =
    blocks::dual_involution(KGB_set.involution(rep.x),W,dual_W);
  return GlobalTitsElement(t,tw);
}

// |rep| final means that no singular real roots satisfy the parity condition
bool Rep_context::is_final(const StandardRepr& rep)
{
  const RootDatum& rd = rootDatum();
  InvolutionNbr inv = kgb().inv_nr(rep.x);
  RootNbrSet pos_real = complexGroup().involution_table().real_roots(inv)
			& rd.posRootSet();
  Weight lambda2_shift=rep.lambda_rho*2 + rd.twoRho()+rd.twoRho(pos_real);

  for (RootNbrSet::iterator it=pos_real.begin(); it(); ++it)
  {
    const Weight& av = rootDatum().coroot(*it);
    if (av.dot(rep.gamma().numerator())==0 and
	av.dot(lambda2_shift)%4 !=0) // singular yet odd on shifted lambda
      return false;
  }
  return true;
}

bool Rep_context::is_oriented(const StandardRepr& rep, RootNbr alpha)
{
  const RootDatum& rd = rootDatum();
  InvolutionNbr inv = kgb().inv_nr(rep.x);
  RootNbrSet real = complexGroup().involution_table().real_roots(inv);

  assert(real.isMember(alpha)); // only real roots should be tested

  const Weight& av = rootDatum().coroot(alpha);
  int numer = av.dot(rep.gamma().numerator());
  int denom = rep.gamma().denominator();
  assert(numer%denom!=0); // and the real root alpha should be non-integral

  Weight lambda2_shift=rep.lambda_rho*2 + rd.twoRho()-rd.twoRho(real);
  int eps = av.dot(lambda2_shift)%4==0 ? 0 : denom;

  return arithmetic::remainder(numer+eps,2*denom)< (unsigned)denom;
}

unsigned int Rep_context::orientation_number(const StandardRepr& rep)
{
  const RootDatum& rd = rootDatum();
  const InvolutionTable& i_tab = complexGroup().involution_table();
  InvolutionNbr inv = kgb().inv_nr(rep.x);
  RootNbrSet real = i_tab.real_roots(inv);
  const Permutation& root_inv = i_tab.root_involution(inv);
  const Weight& numer = rep.gamma().numerator();
  int denom = rep.gamma().denominator();
  Weight lambda2_shift=rep.lambda_rho*2 + rd.twoRho()-rd.twoRho(real);

  unsigned count = 0;

  for (unsigned i=0; i<rd.numPosRoots(); ++i)
  {
    RootNbr alpha = rd.numPosRoots()+i;
    const Weight& av = rootDatum().coroot(alpha);
    int num = av.dot(numer);
    if (num%denom!=0) // skip integral roots
    { if (real.isMember(alpha))
      {
	int eps = av.dot(lambda2_shift)%4==0 ? 0 : denom;
	if ((num>0) == // either positive for gamma and oriente, or neither
	    (arithmetic::remainder(num+eps,2*denom)< (unsigned)denom))
	  ++count;
      }
      else // complex root
      {
	assert(i_tab.complex_roots(inv).isMember(alpha));
	RootNbr beta = root_inv[alpha];
	if (i<rd.rt_abs(beta) // consider only first conjugate "pair"
	    and (num>0)!=(rootDatum().coroot(beta).dot(numer)>0))
	  ++count;
      }
    }
  }
  return count;
}

  } // |namespace repr|
} // |namespace atlas|
