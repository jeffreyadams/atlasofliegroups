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
  InvolutionNbr i_x = kgb().inv_nr(x);
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);

  Weight lambda2 = khc.lift(srk); // doubled coordinates
  RatWeight lambda(lambda2,2);
  RatWeight diff = lambda - nu;
  RatWeight theta_diff(theta*diff.numerator(),
		       diff.denominator()); // theta(lambda-nu)
  Weight lambda_rho = (lambda2-khc.rootDatum().twoRho())/=2;
  return StandardRepr(x,i_tab.pack(i_x,lambda_rho),
		      ((lambda+nu+theta_diff)/=2).normalize());
}

StandardRepr Rep_context::sr
  (KGBElt x, const Weight lambda_rho, const RatWeight& nu) const
{
  InvolutionNbr i_x = kgb().inv_nr(x);
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);
  RatWeight lambda(lambda_rho*2+rootDatum().twoRho(),2);
  RatWeight diff = lambda - nu;
  RatWeight theta_diff(theta*diff.numerator(),
		       diff.denominator()); // theta(lambda-nu)
  return StandardRepr(x,i_tab.pack(i_x,lambda_rho),
		      ((lambda+nu+theta_diff)/=2).normalize());
}

Weight Rep_context::lambda_rho(const StandardRepr& z) const
{
  InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);

  RatWeight gamma_rho = z.gamma() - RatWeight(rootDatum().twoRho(),2);
  Weight im_part2 = gamma_rho.numerator()+theta*gamma_rho.numerator();
  im_part2 /= gamma_rho.denominator(); // exact: $(1+\theta)(\lambda-\rho)$
  return (im_part2 + i_tab.unpack(i_x,z.y_bits))/=2; // division exact again
}

// return $\lambda \in \rho+X^*$ as half-integer rational vector
RatWeight Rep_context::lambda(const StandardRepr& z) const
{
  Weight num = lambda_rho(z) * 2 + rootDatum().twoRho();
  return RatWeight(num,2).normalize();
}

RatWeight Rep_context::nu(const StandardRepr& z) const
{
  InvolutionNbr inv = kgb().inv_nr(z.x());
  const WeightInvolution& theta = complexGroup().involution_table().matrix(inv);
  Weight num = z.gamma().numerator()-theta*z.gamma().numerator();
  return RatWeight(num,2*z.gamma().denominator()).normalize();
}

// convert |lambda| and |gamma| into $y$ value
// Formula: $\exp(i\pi(\gamma-\lambda)) \sigma_{tw} \delta_1$
GlobalTitsElement Rep_context::y(const StandardRepr& z) const
{
  TorusElement t = y_values::exp_pi(z.gamma()-lambda(z));
  const TwistedWeylGroup& W = KGB_set.twistedWeylGroup();
  const TwistedWeylGroup dual_W (W,tags::DualTag()); // VERY EXPENSIVE!
  TwistedInvolution tw =
    blocks::dual_involution(KGB_set.involution(z.x()),W,dual_W);
  return GlobalTitsElement(t,tw);
}

// |z| final means that no singular real roots satisfy the parity condition
bool Rep_context::is_standard(const StandardRepr& z, RootNbr& witness) const
{
  const RootDatum& rd = rootDatum();
  InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = complexGroup().involution_table();
  for (unsigned i=0; i<i_tab.imaginary_rank(i_x); ++i)
  {
    RootNbr alpha = i_tab.imaginary_basis(i_x,i);
    Weight lr = lambda_rho(z);
    int v = lr.scalarProduct(rd.coroot(alpha))+rd.colevel(alpha);
    if (v<0)
      return witness=alpha,false;
  }
  return true;
}

// |z| final means that no singular real roots satisfy the parity condition
bool Rep_context::is_zero(const StandardRepr& z, RootNbr& witness) const
{
  const RootDatum& rd = rootDatum();
  InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = complexGroup().involution_table();
  for (unsigned i=0; i<i_tab.imaginary_rank(i_x); ++i)
  {
    RootNbr alpha = i_tab.imaginary_basis(i_x,i);
    int v = lambda_rho(z).scalarProduct(rd.coroot(alpha))+rd.colevel(alpha);
    bool compact =
      kgb::status(kgb(),z.x(),rd,alpha)==gradings::Status::ImaginaryCompact;
    if (v==0 and compact)
      return witness=alpha,true;
  }
  return false;
}


// |z| final means that no singular real roots satisfy the parity condition
bool Rep_context::is_final(const StandardRepr& z, RootNbr& witness) const
{
  const RootDatum& rd = rootDatum();
  InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = complexGroup().involution_table();
  RootNbrSet pos_real = i_tab.real_roots(i_x) & rd.posRootSet();
  Weight lambda2_shift=lambda_rho(z)*2 + rd.twoRho()-rd.twoRho(pos_real);

  for (RootNbrSet::iterator it=pos_real.begin(); it(); ++it)
  {
    const Weight& av = rootDatum().coroot(*it);
    if (av.dot(z.gamma().numerator())==0 and
	av.dot(lambda2_shift)%4 !=0) // singular yet odd on shifted lambda
      return witness=*it,false;
  }
  return true;
}

bool Rep_context::is_oriented(const StandardRepr& z, RootNbr alpha) const
{
  const RootDatum& rd = rootDatum();
  InvolutionNbr inv = kgb().inv_nr(z.x());
  RootNbrSet real = complexGroup().involution_table().real_roots(inv);

  assert(real.isMember(alpha)); // only real roots should be tested

  const Weight& av = rootDatum().coroot(alpha);
  int numer = av.dot(z.gamma().numerator());
  int denom = z.gamma().denominator();
  assert(numer%denom!=0); // and the real root alpha should be non-integral

  Weight lambda2_shift=lambda_rho(z)*2 + rd.twoRho()-rd.twoRho(real);
  int eps = av.dot(lambda2_shift)%4==0 ? 0 : denom;

  return arithmetic::remainder(numer+eps,2*denom)< (unsigned)denom;
}

unsigned int Rep_context::orientation_number(const StandardRepr& z) const
{
  const RootDatum& rd = rootDatum();
  const InvolutionTable& i_tab = complexGroup().involution_table();
  InvolutionNbr inv = kgb().inv_nr(z.x());
  RootNbrSet real = i_tab.real_roots(inv);
  const Permutation& root_inv = i_tab.root_involution(inv);
  const Weight& numer = z.gamma().numerator();
  int denom = z.gamma().denominator();
  Weight lambda2_shift=lambda_rho(z)*2 + rd.twoRho()-rd.twoRho(real);

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
