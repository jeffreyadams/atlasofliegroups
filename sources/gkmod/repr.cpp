/*
  This is repr.cpp

  Copyright (C) 2009-2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "repr.h"

#include <memory> // for |std::unique_ptr|
#include <map> // used in computing |reducibility_points|
#include <iostream>
#include "error.h"

#include "arithmetic.h"
#include "matreduc.h"

#include "tits.h"
#include "kgb.h"	// various methods
#include "blocks.h"	// |dual_involution|
#include "common_blocks.h" // the |blocks::common_block| class
#include "standardrepk.h"// |KhatContext| methods
#include "subsystem.h" // |SubSystem| methods

#include "kl.h"

#include "ext_block.h"
#include "ext_kl.h"

#include "basic_io.h"

namespace atlas {
  namespace repr {

bool StandardRepr::operator== (const StandardRepr& z) const
{ return x_part==z.x_part and y_bits==z.y_bits
  and infinitesimal_char==z.infinitesimal_char;
}

size_t StandardRepr::hashCode(size_t modulus) const
{ size_t hash=
    x_part + 375*y_bits.data().to_ulong()+83*infinitesimal_char.denominator();
  const Ratvec_Numer_t& num=infinitesimal_char.numerator();
  for (unsigned i=0; i<num.size(); ++i)
    hash= 11*(hash&(modulus-1))+num[i];
  return hash &(modulus-1);
}

StandardReprMod::StandardReprMod (StandardRepr&& sr)
: x_part(sr.x())
, y_bits(sr.y())
, inf_char_mod_1(sr.gamma()) // will be reduced modulo 1 caller; is normalized
{}

StandardReprMod StandardReprMod::mod_reduce
  (const Rep_context& rc, const StandardRepr& sr)
{
  auto gamma_mod1=sr.gamma(); // a normalized rational vector
  auto lam_rho=rc.lambda_rho(sr); // both these valeus are modified below
  auto& num = gamma_mod1.numerator(); assert(num.size()==lam_rho.size());

  const auto d = gamma_mod1.denominator(); // positive value of a signed type
  for (unsigned i=0; i<num.size(); ++i)
  {
    auto q = arithmetic::divide(num[i],d);
    num[i]-= d*q;  // ensure even integral part if |num[i]/d| (0 is even)
    lam_rho[i] -= q; // shift to $\gamma_mod1$ is also applied to $\lambda$ part
  }
  return StandardReprMod(rc.sr_gamma(sr.x(),lam_rho,gamma_mod1));
}

StandardReprMod StandardReprMod::build
  (const Rep_context& rc, const RatWeight& gamma_mod_1, // must be reduced
   KGBElt x, const RatWeight& gam_lam)
{
  const auto gamma_rho = gamma_mod_1 - rho(rc.root_datum());
  const RatWeight lr_rat = (gamma_rho-gam_lam).normalize();
  assert(lr_rat.denominator()==1);
  Weight lam_rho(lr_rat.numerator().begin(),lr_rat.numerator().end());
  return StandardReprMod(rc.sr_gamma(x,lam_rho,gamma_mod_1));
}


size_t StandardReprMod::hashCode(size_t modulus) const
{ size_t hash= x_part +
    243*y_bits.data().to_ulong()+47*inf_char_mod_1.denominator();
  const Ratvec_Numer_t& num=inf_char_mod_1.numerator();
  for (unsigned i=0; i<num.size(); ++i)
    hash= 11*(hash&(modulus-1))+num[i];
  return hash &(modulus-1);
}

Rep_context::Rep_context(RealReductiveGroup &G_R)
  : G(G_R), KGB_set(G_R.kgb())
{}

size_t Rep_context::rank() const { return root_datum().rank(); }

const TwistedInvolution Rep_context::involution_of_Cartan(size_t cn) const
{ return inner_class().involution_of_Cartan(cn); }

StandardRepr Rep_context::sr_gamma
  (KGBElt x, const Weight& lambda_rho, const RatWeight& gamma) const
{ // we use |lambda_rho| only for its real projection |(theta-1)/2*lambda_rho|
  // indeed there is no dependence within its $(1-\theta)(X^*)$-coset either

  int_Matrix theta1 = kgb().involution_matrix(x)+1;
  Weight t1_gamma (gamma.numerator().begin(), gamma.numerator().end());
  // the division in the next computation may throw when |gamma| is bad for |x|
  t1_gamma = theta1*t1_gamma/static_cast<int>(gamma.denominator());
#ifndef NDEBUG // check that constructor below builds a valid |StandardRepr|
  Weight image = // $(\theta+1)(\gamma-\rho)$
    t1_gamma-(theta1*root_datum().twoRho()/2);
  matreduc::find_solution(theta1,image); // a solution must exist
#endif
  const InvolutionTable& i_tab = inner_class().involution_table();
  return StandardRepr(x, i_tab.y_pack(kgb().inv_nr(x),lambda_rho), gamma,
		      height(t1_gamma));
}

StandardRepr Rep_context::sr
  (const StandardReprMod& srm, const RatWeight& gamma) const
{
  const RatWeight gamma_rho = gamma-rho(root_datum());
  const auto lambda_rho = gamma_rho.integer_diff<int>(gamma_lambda(srm));
  return sr_gamma(srm.x(),lambda_rho,gamma);
}

// Height is $\max_{w\in W} \< \rho^v*w , (\theta+1)\gamma >$
unsigned int Rep_context::height(Weight theta_plus_1_gamma) const
{
  const auto& rd=root_datum();
  int result = rd.dual_twoRho().dot(rd.make_dominant(theta_plus_1_gamma));
  assert(result>=0); assert(result%2==0);
  return static_cast<unsigned int>(result/2);
}

RatWeight Rep_context::gamma
  (KGBElt x, const Weight& lambda_rho, const RatWeight& nu) const
{
  const InvolutionTable& i_tab = inner_class().involution_table();
  const RatWeight lambda = rho(root_datum())+lambda_rho;
  const RatWeight diff = lambda - nu;
  const RatWeight theta_diff(i_tab.matrix(kgb().inv_nr(x))*diff.numerator(),
			     diff.denominator()); // theta(lambda-nu)
  return ((lambda+nu+theta_diff)/=2).normalize();
}

StandardRepr
  Rep_context::sr
    (const standardrepk::StandardRepK& srk,
     const standardrepk::SRK_context& srkc,
     const RatWeight& nu) const
{
  const TitsElt a = srkc.titsElt(srk); // was reduced during construction |srk|
  const KGBElt x= kgb().lookup(a);
  Weight lambda_rho = srkc.lift(srk)-root_datum().twoRho();
  lambda_rho/=2; // undo doubled coordinates

  auto result = sr(x,lambda_rho,nu);

  // while standard, final and nonzero, |result| need not be normal
  normalise(result); // prepare |result| for direct inclusion in an |SR_poly|
  return result;
}

const WeightInvolution& Rep_context::theta (const StandardRepr& z) const
{ return inner_class().involution_table().matrix(kgb().inv_nr(z.x())); }

Weight Rep_context::lambda_rho(const StandardRepr& z) const
{
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = inner_class().involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);

  const RatWeight gamma_rho = z.gamma() - rho(root_datum());
  Ratvec_Numer_t im_part2 = gamma_rho.numerator()+theta*gamma_rho.numerator();
  im_part2 /= gamma_rho.denominator(); // exact: $(1+\theta)(\lambda-\rho)$
  Weight i2(im_part2.begin(),im_part2.end()); // convert to |Weight|
  return (i2 + i_tab.y_lift(i_x,z.y()))/2; // division exact again
}

// this function is similar to |common_block::gamma_lambda|
RatWeight Rep_context::gamma_lambda
  (InvolutionNbr i_x,  const TorusPart& y_bits, const RatWeight& gamma) const
{
  const InvolutionTable& i_tab = inner_class().involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);

  const RatWeight gamma_rho = gamma - rho(root_datum());
  return (gamma_rho-theta*gamma_rho - i_tab.y_lift(i_x,y_bits))/2;
}

RatWeight Rep_context::gamma_lambda(const StandardReprMod& z) const
{ return gamma_lambda(kgb().inv_nr(z.x()),z.y(),z.gamma_mod1()); }

RatWeight Rep_context::gamma_0 (const StandardRepr& z) const
{
  const InvolutionTable& i_tab = inner_class().involution_table();
  const auto& theta = i_tab.matrix(kgb().inv_nr(z.x()));
  return ((z.gamma()+theta*z.gamma())/=2).normalize();
}

RatWeight Rep_context::nu(const StandardRepr& z) const
{
  const InvolutionTable& i_tab = inner_class().involution_table();
  const auto& theta = i_tab.matrix(kgb().inv_nr(z.x()));
  return ((z.gamma()-theta*z.gamma())/=2).normalize();
}

TorusElement Rep_context::y_as_torus_elt(const StandardRepr& z) const
{ return y_values::exp_pi(z.gamma()-lambda(z)); }

TorusElement Rep_context::y_as_torus_elt(const StandardReprMod& z) const
{ return y_values::exp_pi(gamma_lambda(z)); }

// |z| standard means (weakly) dominant on the (simple-)imaginary roots
bool Rep_context::is_standard(const StandardRepr& z, RootNbr& witness) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = inner_class().involution_table();
  const Ratvec_Numer_t& numer = z.gamma().numerator();

  for (unsigned i=0; i<i_tab.imaginary_rank(i_x); ++i)
  {
    const RootNbr alpha = i_tab.imaginary_basis(i_x,i);
    if (rd.coroot(alpha).dot(numer)<0)
      return witness=alpha,false;
  }
  return true;
}

// |z| dominant means precisely |gamma| is (weakly) dominant
bool Rep_context::is_dominant(const StandardRepr& z, RootNbr& witness) const
{
  const RootDatum& rd = root_datum();
  const auto& numer = z.gamma().numerator();

  for (auto it=rd.beginSimpleCoroot(); it!=rd.endSimpleCoroot(); ++it)
    if (it->dot(numer)<0)
      return witness=rd.simpleRootNbr(it-rd.beginSimpleCoroot()),false;
  return true;
}

// |z| zero means that no singular simple-imaginary roots are compact; this
// code assumes |is_standard(z)|, namely |gamma| is dominant on imaginary roots
bool Rep_context::is_nonzero(const StandardRepr& z, RootNbr& witness) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = inner_class().involution_table();
  const Ratvec_Numer_t& numer = z.gamma().numerator();

  for (unsigned i=0; i<i_tab.imaginary_rank(i_x); ++i)
  {
    const RootNbr alpha = i_tab.imaginary_basis(i_x,i);
    if (rd.coroot(alpha).dot(numer)==0 and // simple-imaginary, singular
	not kgb().simple_imaginary_grading(z.x(),alpha)) // and compact
      return witness=alpha,false;
  }
  return true;
}

bool Rep_context::is_normal(const StandardRepr& z) const
{
  auto z_normal = z;
  normalise(z_normal);
  return z_normal==z;
}

// |z| final means that no singular real roots satisfy the parity condition
// we do not assume |gamma| to be dominant, so all real roots must be tested
bool Rep_context::is_semifinal(const StandardRepr& z, RootNbr& witness) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = inner_class().involution_table();
  const RootNbrSet pos_real = i_tab.real_roots(i_x) & rd.posRootSet();
  const Weight test_wt = i_tab.y_lift(i_x,z.y()) // $(1-\theta)(\lambda-\rho)$
	   + rd.twoRho()-rd.twoRho(pos_real); // replace $\rho$ by $\rho_R$

  for (RootNbrSet::iterator it=pos_real.begin(); it(); ++it)
  {
    const Weight& av = root_datum().coroot(*it);
    if (av.dot(z.gamma().numerator())==0 and
	av.dot(test_wt)%4 !=0) // singular yet odd on shifted lambda
      return witness=*it,false;
  }
  return true;
}

bool Rep_context::is_final(const StandardRepr& z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionTable& i_tab = inner_class().involution_table();
  const auto& numer = z.gamma().numerator();
  KGBElt x=z.x();
  const InvolutionNbr i_x = kgb().inv_nr(x);

  for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
  {
    auto v=rd.simpleCoroot(s).dot(numer);
    if (v<0)
      return false; // unless |gamma| is dominant, we just say "no"
    else if (v==0)
      switch (kgb().status(s,x))
      {
      case gradings::Status::Complex:
	if (kgb().isDescent(s,x))
	  return false;
	break;
      case gradings::Status::ImaginaryCompact:
	return false; // certainly fails |is_nonzero|
      case gradings::Status::Real:
	if (rd.simpleCoroot(s).dot(i_tab.y_lift(i_x,z.y()))%4!=0)
	  return false;
      default: {} // ImaginaryNoncompact is fine
      }
  }
  RootNbr witness;
  return is_nonzero(z,witness); // check \emph{all} simply-imaginary coroots
}


bool Rep_context::is_oriented(const StandardRepr& z, RootNbr alpha) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = inner_class().involution_table();
  const RootNbrSet real = inner_class().involution_table().real_roots(i_x);

  assert(real.isMember(alpha)); // only real roots should be tested

  const Weight& av = root_datum().coroot(alpha);
  const auto numer = av.dot(z.gamma().numerator());
  const auto denom = z.gamma().denominator();
  assert(numer%denom!=0); // and the real root alpha should be non-integral

  const Weight test_wt =
    i_tab.y_lift(i_x,z.y()) +rd.twoRho() -rd.twoRho(real);
  const auto eps = av.dot(test_wt)%4==0 ? 0 : denom;

  return arithmetic::remainder(numer+eps,2*denom) < denom;
}

unsigned int Rep_context::orientation_number(const StandardRepr& z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionTable& i_tab = inner_class().involution_table();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const RootNbrSet real = i_tab.real_roots(i_x);
  const Permutation& root_inv = i_tab.root_involution(i_x);
  const Ratvec_Numer_t& numer = z.gamma().numerator();
  const arithmetic::Numer_t denom = z.gamma().denominator();
  const Weight test_wt = // representative of a class modulo $2(1-\theta)(X^*)$
    i_tab.y_lift(i_x,z.y()) +rd.twoRho() -rd.twoRho(real);

  unsigned count = 0;

  for (unsigned i=0; i<rd.numPosRoots(); ++i)
  {
    const RootNbr alpha = rd.numPosRoots()+i;
    const Coweight& av = root_datum().coroot(alpha);
    const arithmetic::Numer_t num = av.dot(numer);
    if (num%denom!=0) // skip integral roots
    { if (real.isMember(alpha))
      {
	int eps = av.dot(test_wt)%4==0 ? 0 : denom;
	if ((num>0) == // either positive for gamma and oriented, or neither
	    (arithmetic::remainder(num+eps,2*denom)< (unsigned)denom))
	  ++count;
      }
      else // complex root
      {
	assert(i_tab.complex_roots(i_x).isMember(alpha));
	const RootNbr beta = root_inv[alpha];
	if (i<rd.rt_abs(beta) // consider only first conjugate "pair"
	    and (num>0)!=(root_datum().coroot(beta).dot(numer)>0))
	  ++count;
      }
    }
  }
  return count;
} // |orientation_number|

void Rep_context::make_dominant(StandardRepr& z) const
{
  const RootDatum& rd = root_datum();

  // the following are non-|const|, and modified in the loop below
  Weight lr = lambda_rho(z);
  KGBElt& x = z.x_part; // the |x_part| will be modified in-place
  Ratvec_Numer_t& numer = z.infinitesimal_char.numerator();

  { weyl::Generator s;
    do
      for (s=0; s<rd.semisimpleRank(); ++s)
      {
	int v=rd.simpleCoroot(s).dot(numer);
	if (v<0)
	{
	  rd.simple_reflect(s,numer);
	  int offset; // used to pivot |lr| around $\rho_r-\rho$
	  switch (kgb().status(s,x))
	  {
	  case gradings::Status::Complex: offset = 1; break;
	  case gradings::Status::Real:    offset = 0; break;
	  default: // |s| is an imaginary root; we will not cope with that here
	    throw std::runtime_error("Non standard parameter in make_dominant");
	  }
	  rd.simple_reflect(s,lr,offset);
	  x = kgb().cross(s,x);
	  break; // out of the loop |for(s)|
	} // |if(v<0)|
      } // |for(s)|
    while (s<rd.semisimpleRank()); // wait until inner loop runs to completion
  }
  z.y_bits=inner_class().involution_table().y_pack(kgb().inv_nr(x),lr);
} // |make_dominant|

void Rep_context::singular_cross (weyl::Generator s,StandardRepr& z) const
{
  assert(root_datum().simpleCoroot(s).dot(z.gamma().numerator())==0);
  Weight lr = lambda_rho(z); auto& x=z.x_part;
  root_datum().simple_reflect
    (s, lr, kgb().status(s,x)==gradings::Status::Real ? 0 : 1);
  x = kgb().cross(s,x);
  z.y_bits = // reinsert $y$ bits component
    inner_class().involution_table().y_pack(kgb().inv_nr(x),lr);
}

// auxiliary: move to a canonical for the |gens| (singular) subgroup of $W$
void Rep_context::to_singular_canonical(RankFlags gens, StandardRepr& z) const
{ // simply-singular coroots are simple, so no need to constuct a subsystem
  TwistedInvolution tw = kgb().involution(z.x_part); // copy to be modified
  WeylWord ww = inner_class().canonicalize(tw,gens);
  for (auto it=ww.begin(); it!=ww.end(); ++it) // move to that involution
    singular_cross(*it,z);
  assert(tw == kgb().involution(z.x_part));
}

void Rep_context::normalise(StandardRepr& z) const
{
  make_dominant(z);
  const RootDatum& rd = root_datum();

  RankFlags simple_singulars;
  { const auto& numer = z.infinitesimal_char.numerator();
    for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
      simple_singulars.set(s,rd.simpleCoroot(s).dot(numer)==0);
  }

  to_singular_canonical(simple_singulars,z);

  // the following are non-|const|, and modified in the loop below
  Weight lr = lambda_rho(z);
  KGBElt& x = z.x_part;

  { RankFlags::iterator it;
    do
      for (it=simple_singulars.begin(); it(); ++it)
	if (kgb().isComplexDescent(*it,x))
	{
	  weyl::Generator s=*it;
	  rd.simple_reflect(s,lr,1); // pivot |lr| around $-\rho$
	  x = kgb().cross(s,x);
	  break; // out of the loop |for(s)|
	} // |if(v<0)|
    while (it()); // wait until inner loop runs to completion
  }
  z.y_bits=inner_class().involution_table().y_pack(kgb().inv_nr(x),lr);
} // |normalise|

bool Rep_context::is_twist_fixed
  (StandardRepr z, const WeightInvolution& delta) const
{
  make_dominant(z);
  const RootDatum& rd = root_datum();

  RankFlags simple_singulars;
  { const auto& numer = z.infinitesimal_char.numerator();
    for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
      simple_singulars.set(s,rd.simpleCoroot(s).dot(numer)==0);
  }

  to_singular_canonical(simple_singulars,z);

  return z==twisted(z,delta);
} // |normalise|

// equivalence is equality after |make_dominant| and |to_singular_canonical|
bool Rep_context::equivalent(StandardRepr z0, StandardRepr z1) const
{
  if (kgb().Cartan_class(z0.x_part)!=kgb().Cartan_class(z1.x_part))
    return false; // this non-equivalence can be seen before |make_dominant|

  { // preempt failing of |make_dominant| on non standard parameters
    RootNbr witness;
    if (not (is_standard(z0,witness) and is_standard(z1,witness)))
      return z0==z1; // strict equality unless both are parameters standard
  }

  make_dominant(z0);
  make_dominant(z1);

  if (z0.infinitesimal_char!=z1.infinitesimal_char)
    return false;

  const RootDatum& rd = root_datum();

  RankFlags simple_singulars;
  { const auto& numer = z0.infinitesimal_char.numerator();
    for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
      simple_singulars.set(s,rd.simpleCoroot(s).dot(numer)==0);
  }

  to_singular_canonical(simple_singulars,z0);
  to_singular_canonical(simple_singulars,z1);

  return z0==z1;
}

StandardRepr& Rep_context::scale(StandardRepr& z, const Rational& f) const
{ // we can just replace the |infinitesimal_char|, nothing else changes
  auto image = theta(z)*z.gamma();
  auto diff = z.gamma()-image; // this equals $2\nu(z)$
  z.infinitesimal_char += image;
  z.infinitesimal_char += diff*f; // now we have |(gamma_0(z)+nu(z)*f)*2|
  (z.infinitesimal_char/=2).normalize();
  return z;
}

StandardRepr& Rep_context::scale_0(StandardRepr& z) const
{ z.infinitesimal_char = gamma_0(z); return z; }

RationalList Rep_context::reducibility_points(const StandardRepr& z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = inner_class().involution_table();
  const Permutation& theta = i_tab.root_involution(i_x);

  const RatWeight& gamma = z.gamma();
  const Ratvec_Numer_t& numer = gamma.numerator();
  const arithmetic::Numer_t d = gamma.denominator();
  const Weight lam_rho = lambda_rho(z);

  const RootNbrSet pos_real = i_tab.real_roots(i_x) & rd.posRootSet();
  const Weight two_rho_real = rd.twoRho(pos_real);

  // we shall associate to certain numbers $num>0$ a strict lower bound $lwb$
  // for which we shall then later form fractions $(d/num)*k$ for $k>lwb$
  typedef std::map<long,long> table;

  // because of the parity condition, distinguish cases with even and odd $k$
  table odds,evens; // name indicates the parity that $k$ will have

  for (RootNbrSet::iterator it=pos_real.begin(); it(); ++it)
  {
    arithmetic::Numer_t num =
      rd.coroot(*it).dot(numer); // now $\<\alpha^v,\nu>=num/d$ (real $\alpha$)
    if (num!=0)
    {
      long lam_alpha = lam_rho.dot(rd.coroot(*it))+rd.colevel(*it);
      bool do_odd = (lam_alpha+two_rho_real.dot(rd.coroot(*it))/2)%2 ==0;
      (do_odd ? odds : evens).insert(std::make_pair(std::abs(num),0));
    }
  }

  RootNbrSet pos_complex = i_tab.complex_roots(i_x) & rd.posRootSet();
  for (RootNbrSet::iterator it=pos_complex.begin(); it(); ++it)
  {
    RootNbr alpha=*it, beta=theta[alpha];
    arithmetic::Numer_t vala = rd.coroot(alpha).dot(numer);
    arithmetic::Numer_t valb = rd.coroot(beta).dot(numer);
    arithmetic::Numer_t num = vala - valb; // $2\<\alpha^v,\nu>=num/d$ (complex)
    if (num!=0)
    {
      assert((vala+valb)%d==0); // since $\<a+b,\gamma>=\<a+b,\lambda>$
      long lwb =std::abs(vala+valb)/d;
      std::pair<table::iterator,bool> trial = // try insert |lwb| as |num| value
	(lwb%2==0 ? evens : odds).insert(std::make_pair(std::abs(num),lwb));
      if (not trial.second and lwb<trial.first->second)
	trial.first->second=lwb; // if not new, maybe lower the old bound value
    }
  }

  std::set<Rational> fracs;

  for (table::iterator it= evens.begin(); it!=evens.end(); ++it)
    for (long s= d*(it->second+2); s<=it->first; s+=2*d)
      fracs.insert(Rational(s,it->first));

  for (table::iterator it= odds.begin(); it!=odds.end(); ++it)
    for (long s= it->second==0 ? d : d*(it->second+2); s<=it->first; s+=2*d)
      fracs.insert(Rational(s,it->first));

  return RationalList(fracs.begin(),fracs.end());
} // |reducibility_points|


StandardRepr Rep_context::cross(weyl::Generator s, StandardRepr z) const
{
  make_dominant(z);
  const RatWeight infin_char=z.gamma(); // now get the infinitesimal character
  const RootDatum& rd = root_datum();
  const SubSystem& subsys = SubSystem::integral(rd,infin_char);
  blocks::nblock_help aux(real_group(),subsys);
  blocks::nblock_elt src(z.x(),y_as_torus_elt(z));
  aux.cross_act(s,src);
  const RatWeight& t =	src.y().as_Qmod2Z();
  // InvolutionNbr i_x = kgb().inv_nr(z.x());
  // no need to do |inner_class().involution_table().real_unique(i_x,t)|

  RatWeight lr =  (infin_char - t - rho(rd)).normalize();
  assert(lr.denominator()==1); // we have reconstructed $\lambda-\rho \in X^*$
  return sr_gamma(src.x(),
		  Weight(lr.numerator().begin(),lr.numerator().end()), // mod 2
		  infin_char);
}

StandardRepr Rep_context::cross(const Weight& alpha, StandardRepr z) const
{
  const RootDatum& rd = root_datum();
  KGBElt& x= z.x_part;
  InvolutionNbr i_x = kgb().inv_nr(x);
  const InvolutionTable& i_tab = inner_class().involution_table();

  const RatWeight& gamma=z.infinitesimal_char; // integrally dominant
  RootNbr rt = rd.root_index(alpha);
  if (rt==rd.numRoots())
    throw std::runtime_error("Not a root");
  // the following test ensures that our final |assert| below won't fail
  if (rd.coroot(rt).dot(gamma.numerator())%gamma.denominator()!=0)
    throw std::runtime_error("Not an integral root");

  RatWeight lambda_shifted =
    gamma - lambda(z) + RatWeight(rd.twoRho(i_tab.real_roots(i_x)),2);
  Ratvec_Numer_t& lambda_numer = lambda_shifted.numerator();

  rd.reflect(rt,lambda_numer);
  x = kgb().cross(rd.reflectionWord(rt),x);
  i_x = kgb().inv_nr(x);

  // the addition of $\rho$ below is because |sr_gamma| takes $\lambda-\rho$
  lambda_shifted += RatWeight(rd.twoRho()-rd.twoRho(i_tab.real_roots(i_x)),2);
  lambda_shifted =  (gamma - lambda_shifted).normalize();
  assert(lambda_shifted.denominator()==1);

  return sr_gamma(x,Weight(lambda_numer.begin(),lambda_numer.end()),gamma);
}

StandardRepr Rep_context::Cayley(weyl::Generator s, StandardRepr z) const
{
  make_dominant(z);
  const RatWeight infin_char=z.gamma(); // now get the infinitesimal character
  const RootDatum& rd = root_datum();
  const SubSystem& subsys = SubSystem::integral(rd,infin_char);
  blocks::nblock_help aux(real_group(),subsys);
  blocks::nblock_elt src(z.x(),y_as_torus_elt(z));
  aux.do_up_Cayley(s,src);
  RatWeight t =	 src.y().log_pi(false);
  // InvolutionNbr i_x = kgb().inv_nr(z.x());
  // no need to do |inner_class().involution_table().real_unique(i_x,t)|

  RatWeight lr =  (infin_char - t - rho(rd)).normalize();
  assert(lr.denominator()==1);
  return sr_gamma(src.x(),
		  Weight(lr.numerator().begin(),lr.numerator().end()), // mod 2
		  infin_char);
}

StandardRepr Rep_context::inv_Cayley(weyl::Generator s, StandardRepr z) const
{
  make_dominant(z);
  const RatWeight infin_char=z.gamma(); // now get the infinitesimal character
  const RootDatum& rd = root_datum();
  const SubSystem& subsys = SubSystem::integral(rd,infin_char);
  blocks::nblock_help aux(real_group(),subsys);
  blocks::nblock_elt src(z.x(),y_as_torus_elt(z));
  aux.do_down_Cayley(s,src);
  RatWeight t =	 src.y().log_pi(false);
  // InvolutionNbr i_x = kgb().inv_nr(z.x());
  // no need to do |inner_class().involution_table().real_unique(i_x,t)|

  RatWeight lr =  (infin_char - t - rho(rd)).normalize();
  assert(lr.denominator()==1);
  return sr_gamma(src.x(),
		  Weight(lr.numerator().begin(),lr.numerator().end()), // mod 2
		  infin_char);
}


/*
  Compute shift in |lambda| component of parameter for Cayley transform by a
  non-simple root $\alpha$, from involutions |theta_down| to |theta_up|, where
  |to_simple| left-conjugates root $\alpha$ to some simple root.

  Curiously, this appears to depend only on $\theta$ upstairs and the
  conjugating element |to_simple|; an explanation is needed here. It seems to
  be because \emph{all} upstairs real roots becoming negative by the necessary
  conjugation will be downstairs complex roots (so contribute to the shift).

  Sum of positve real roots becoming negative at $\theta'=^{to\_simple}\theta$
*/
Weight Cayley_shift (const InnerClass& G,
		     InvolutionNbr theta_upstairs, // at the more split Cartan
		     const WeylWord& to_simple)
{ const RootDatum& rd=G.rootDatum();
  const InvolutionTable& i_tab = G.involution_table();
  RootNbrSet S = pos_to_neg(rd,to_simple) & i_tab.real_roots(theta_upstairs);
  return root_sum(rd,S);
}

// a method used to ensure |z| is integrally dominant, used by |any_Cayley|
WeylWord
Rep_context::make_dominant(StandardRepr& z,const SubSystem& subsys) const
{
  const RootDatum& rd = root_datum();
  KGBElt& x= z.x_part;
  InvolutionNbr i_x = kgb().inv_nr(x);
  const InvolutionTable& i_tab = inner_class().involution_table();

  // the following are non-|const|, and modified in the loop below
  Weight lambda2_shifted = (lambda_rho(z)*=2)
    + rd.twoRho() - rd.twoRho(i_tab.real_roots(i_x));
  Ratvec_Numer_t& gamma_num = z.infinitesimal_char.numerator();

  WeylWord result;
  result.reserve(subsys.numPosRoots()); // enough to accommodate the WeylWord

  { weyl::Generator s;
    do
    {
      for (s=0; s<subsys.rank(); ++s)
      {
	RootNbr alpha = subsys.parent_nr_simple(s);
	arithmetic::Numer_t v=rd.coroot(alpha).dot(gamma_num);
	if (v<0)
	{
	  if (i_tab.imaginary_roots(i_x).isMember(alpha))
	    throw std::runtime_error
	      ("Cannot make non-standard parameter integrally dominant");
	  result.push_back(s);

	  // reflect |gamma| by |alpha|
	  gamma_num.subtract(rd.root(alpha).begin(),v);
	  x = kgb().cross(rd.reflectionWord(alpha),x);
	  i_x = kgb().inv_nr(x);
	  rd.reflect(alpha,lambda2_shifted);
	  break; // out of the loop |for(s)|
	} // |if(v<0)|
      } // |for(s)|
    }
    while (s<subsys.rank()); // wait until inner loop runs to completion
  }
  lambda2_shifted -= rd.twoRho() - rd.twoRho(i_tab.real_roots(i_x)); // unshift
  z.y_bits=i_tab.y_pack(i_x,lambda2_shifted/2);
  return result;
} // |make_dominant| (integrally)

StandardRepr Rep_context::any_Cayley(const Weight& alpha, StandardRepr z) const
{
  const RootDatum& rd = root_datum();
  const KGB& kgb = this->kgb();
  const InvolutionTable& i_tab = inner_class().involution_table();
  const SubSystem& subsys = SubSystem::integral(rd,z.infinitesimal_char);

  // prepare: move to a situation with integrally dominant infinitesimal char.
  WeylWord w=make_dominant(z,subsys);
  KGBElt x= z.x_part; // take a working copy; don't disturb |z|
  Weight lr = lambda_rho(z); // use at end to build new parameter
  const RatWeight& infin_char=z.infinitesimal_char; // constant from here on

  // check the root argument, and if OK make the corresponding move to above
  RootNbr rt = subsys.from_parent(rd.root_index(alpha)); // |subsys| numbering
  if (rt == RootNbr(-1)) // either not a root at all or not in subsystem
    throw std::runtime_error("Not an integral root");
  // apply the integrally-dominant-making $W$ element |w| (in |subsys|) to |rt|:
  rt = subsys.permuted_root(rt,w); // now we've got the root to do Cayley by
  const RootNbr n_alpha = subsys.to_parent(rt); // for modified |alpha|

  rt = subsys.rt_abs(rt); // interpret as index of a positive root
  weyl::Generator s=subsys.simple(rt);
  WeylWord ww = subsys.to_simple(rt);
  // neither |alpha| nor |rt| will be used beyond thus point

  // now do the Cayley transform proper
  bool ascent; // whether forward Cayley
  const InvolutionNbr inv0= kgb.inv_nr(x); // initial involution

  x = kgb.cross(ww,x);
  switch (kgb.status(s,x))
  {
  case gradings::Status::ImaginaryNoncompact:
    x = kgb.cayley(s,x); ascent=true; break;
  case gradings::Status::Real: // find out (at inv0) whether root is parity
    { Weight rho2_diff = rd.twoRho() - rd.twoRho(i_tab.real_roots(inv0));
      RatWeight parity_vector = // compute this at the \emph{original} x
	infin_char - lr - RatWeight(std::move(rho2_diff),2);
      if (parity_vector.dot(rd.coroot(n_alpha))%2!=0)
      { // then |alpha| was parity
	x = kgb.inverseCayley(s,x).first; // do inverse Cayley at |inv1|
	ascent=false;
	break;
      }
      // else FALL THROUGH
    }
  case gradings::Status::ImaginaryCompact:
  case gradings::Status::Complex:
    throw error::Cayley_error();
  }
  x = kgb.cross(x,ww); // finally cross back

  lr += // apply shift depending on distance from being simply-real upstairs
    Cayley_shift(inner_class(),ascent ? kgb.inv_nr(x) : inv0,ww);
  z = sr_gamma(x,lr,infin_char);

  return z;
} // |Rep_context::any_Cayley|

StandardRepr Rep_context::inner_twisted(StandardRepr z) const
{
  make_dominant(z);
  RatWeight gamma=z.gamma(); // now get the infinitesimal character
  const RootDatum& rd = root_datum();
  const SubSystem& subsys = SubSystem::integral(rd,gamma);
  blocks::nblock_help aux(real_group(),subsys);
  blocks::nblock_elt src(z.x(),y_as_torus_elt(z));

  inner_class().distinguished().apply_to(gamma.numerator()); // twist |gamma|
  aux.twist(src);
  RatWeight lr = (gamma - src.y().log_pi(false) - rho(rd)).normalize();
  assert(lr.denominator()==1);
  return sr_gamma(src.x(),
		  Weight(lr.numerator().begin(),lr.numerator().end()), // mod 2
		  gamma);
}

StandardRepr Rep_context::twisted
  (StandardRepr z, const WeightInvolution& delta) const
{
  const auto& i_tab = inner_class().involution_table();
  const InvolutionNbr i_x0 = kgb().inv_nr(z.x());
  z.x_part = kgb().twisted(z.x_part,delta);
  const InvolutionNbr i_x1 =  kgb().inv_nr(z.x()); // destination involution
  z.y_bits = i_tab.y_act(i_x0,i_x1,z.y_bits,delta);
  z.infinitesimal_char = delta*z.infinitesimal_char;
  return z;
}

Rep_context::compare Rep_context::repr_less() const
{ return compare(root_datum().dual_twoRho()); }

bool Rep_context::compare::operator()
  (const StandardRepr& r,const StandardRepr& s) const
{
  if (r.height()!=s.height()) // order by increasing height first
    return r.height()<s.height();
  if (r.x()!=s.x()) // then order by decreasing numeric value of |x|
    return r.x()>s.x(); // (height tends to change in opposite sense to |x|)
  if (r.y()!=s.y()) // then order by increasing internal value of |y|
    return r.y()<s.y(); // uses |SmallBitVector::operator<|, internal comparison

  // finally in rare cases individual components of |gamma| need comparison
  auto r_vec = s.gamma().numerator()*r.gamma().denominator(); // cross multiply
  auto s_vec = r.gamma().numerator()*s.gamma().denominator(); // cross multiply

  return r_vec<s_vec;
}


// |Rep_table| methods

Rep_table::Rep_table(RealReductiveGroup &G)
: Rep_context(G)
, pool(), hash(pool), def_formulae()
, mod_pool(), mod_hash(mod_pool), block_list(), place()
{}
Rep_table::~Rep_table() = default;

unsigned short Rep_table::length(StandardRepr sr)
{
  make_dominant(sr); // length should not change in equivalence class
  BlockElt z;
  auto & block = lookup(sr,z); // construct partial block
  return block.length(z);
}

SR_poly Rep_context::scale(const poly& P, const Rational& f) const
{
  poly result(repr_less());
  for (auto it=P.begin(); it!=P.end(); ++it)
  { auto z=it->first; // take a copy for modification
    auto finals = finals_for(scale(z,f));
    for (const StandardRepr& final : finals)
      result.add_term(final,it->second);
  }
  return result;
}

SR_poly Rep_context::scale_0(const poly& P) const
{
  poly result(repr_less());
  for (auto it=P.begin(); it!=P.end(); ++it)
  { auto z=it->first; // take a copy for modification
    auto finals = finals_for(scale_0(z));
    for (const StandardRepr& final : finals)
      result.add_term(final,it->second);
  }
  return result;
}

containers::sl_list<StandardRepr>
  Rep_context::finals_for(StandardRepr z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionTable& i_tab = inner_class().involution_table();

  make_dominant(z); // ensures singular subsystem is generated by simple roots
  const RatWeight& gamma=z.gamma();

  RankFlags singular; // subset of simple roots that is singular at |gamma|
  for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
    singular.set(s,rd.simpleCoroot(s).dot(gamma.numerator())==0);
/* the simple roots flagged in |singular| coincide with the singular
   integrally-simple roots (simple roots of integral subsystem), so we can
   proceed without constructing the integral subsystem (even less the block) */

  containers::sl_list<StandardRepr> result { z };
  auto rit=result.begin();
  do
  { const KGBElt x=rit->x();
    auto it=singular.begin();
    for (; it(); ++it)
    { const auto s=*it;
      if (kgb().status(s,x)==gradings::Status::ImaginaryCompact)
      {
        result.erase(rit); // discard zero parameter
	break;
      }
      else if (kgb().status(s,x)==gradings::Status::Complex)
      { if (kgb().isDescent(s,x))
	{ // replace |*rit| by its complex descent for |s|
	  singular_cross(s,*rit);
	  break; // reconsider all singular roots for the new parameter
	}
      }
      else if (kgb().status(s,x)==gradings::Status::Real)
      { // only do something for parity roots
	const InvolutionNbr i_x = kgb().inv_nr(x);
	if (rd.simpleCoroot(s).dot(i_tab.y_lift(i_x,rit->y()))%4!=0)
	{ // found parity root; |kgb()| can distinguish type 1 and type 2
	  const KGBEltPair p = kgb().inverseCayley(s,x);
	  Weight lr = lambda_rho(*rit);
	  assert(rd.simpleCoroot(s).dot(lr)%2!=0); // parity says this
	  *rit = sr_gamma(p.first,lr,gamma); // replace by first inverse Cayley
	  if (p.second!=UndefKGB) // insert second inverse Cayley after first
	    result.insert(std::next(rit),sr_gamma(p.second,lr,gamma));
	  break;
	}
      }
    }
    if (not it()) // no singular descents found: consolidate |*rit|
      ++rit;
  }
  while (not result.at_end(rit));
  return result;
} // |Rep_context::finals_for|

SR_poly Rep_context::expand_final (StandardRepr z) const
{
  auto finals = finals_for(z);
  poly result (repr_less());
  for (auto it=finals.cbegin(); not finals.at_end(it); ++it)
    result += *it;
  return result;
} // |Rep_context::expand_final|

unsigned long Rep_table::formula_index (const StandardRepr& sr)
{
  const auto prev_size = hash.size();
  const auto h = hash.match(sr);
  if (h>=prev_size)
    def_formulae.push_back
      (std::make_pair(SR_poly(repr_less()),SR_poly(repr_less())));
  return h;
}

unsigned long Rep_table::add_block(const StandardReprMod& srm)
{
  auto first=mod_hash.size(); // future code of first element of this block
  BlockElt srm_in_block; // will hold position of |srm| within that block
  std::unique_ptr<blocks::common_block>
    ptr(new blocks::common_block(*this,srm,srm_in_block));
  auto& block=*ptr;
  block_list.push_back(std::move(ptr));

  const unsigned long result = // future sequence number for our |srm|
    first+srm_in_block;
  const RatWeight gamma_rho = srm.gamma_mod1()-rho(root_datum());
  for (BlockElt z=0; z<block.size(); ++z)
  {
    Weight lambda_rho=gamma_rho.integer_diff<int>(block.gamma_lambda(z));
    auto zm = StandardReprMod::mod_reduce
      (*this, sr_gamma(block.x(z),lambda_rho,srm.gamma_mod1()));
    auto seq = mod_hash.match(zm);
    assert(seq==place.size()); // all block elements should be new
    ndebug_use(seq);
    place.emplace_back(&block,z);
  }
  return result;
}

blocks::common_block& Rep_table::lookup
  (const StandardRepr& sr,BlockElt& z,RankFlags& singular)
{
  auto srm = StandardReprMod::mod_reduce(*this,sr); // modular |z|
  auto h=mod_hash.find(srm); // look up modulo translation in $X^*$
  if (h==mod_hash.empty) // then we are in a new translation family of blocks
    h=add_block(srm); // ensure this block is known
  assert(h<place.size()); // it cannot be |mod_hash.empty| anymore

  auto& block = *place[h].first;
  z = place[h].second;

  singular.reset();
  {
    const SubSystem& sub = block.integral_subsystem();
    for (weyl::Generator s=0; s<block.rank(); ++s)
      singular.set(s,root_datum().coroot(sub.parent_nr_simple(s))
					.dot(sr.gamma().numerator())==0);
  }
  return block;
} // |Rep_table::lookup|

blocks::common_block& Rep_table::lookup (const StandardRepr& sr,BlockElt& which)
{
  auto srm = StandardReprMod::mod_reduce(*this,sr); // modular |z|
  assert(mod_hash.size()==place.size()); // should be in sync at this point
  auto h=mod_hash.find(srm); // look up modulo translation in $X^*$
  if (h!=mod_hash.empty) // then we are in a new translation family of blocks
  {
    assert(h<place.size()); // it cannot be |mod_hash.empty| anymore
    which = place[h].second;
    return *place[h].first;
  }
  common_context ctxt(real_group(),SubSystem::integral(root_datum(),sr.gamma()));
  BitMap subset;
  auto& block= add_block_below(ctxt,srm,&subset); // ensure block is known
  which = last(subset);
  assert(block.representative(which)==srm);
  return block;
} // |Rep_table::lookup|

// in the following type the second component is a mulitplicity so we are in fact
// dealing with a sparse reprensetion of polynomials with |BlockElt| exponents

typedef std::pair<BlockElt,int>  term;
typedef containers::sl_list<term> pair_list;

pair_list combine (pair_list a, pair_list b) // by value, will move from rvalues
{ // |a| and |b| are assumed to be sorted
  a.merge(std::move(b),
	  [](const term& x, const term& y) { return x.first<y.first; });
  // now any like terms are neigbours, combine them whenever this occurs
  for (auto it=a.begin(); not a.at_end(it); ++it)
  {
    auto it1=std::next(it);
    if (not a.at_end(it1) and it->first==it1->first)
    {
      it->second += it1->second;
      a.erase(it1);
    }
  }
  return a;
}

pair_list flip (int sign, pair_list list) // by value
{ if (sign!=1)
    for (auto& p : list)
      p.second *= sign;
  return list;
}

std::vector<pair_list> contributions
  (blocks::common_block& block, RankFlags singular, BlockElt y)
{
  std::vector<pair_list> result(y+1); // initally every |result[z]| is empty
  for (BlockElt z=0; z<=y; ++z) // compute |finals| and |finals_for| in |result|
  {
    const DescentStatus& desc=block.descent(z);
    auto it=singular.begin();
    for ( ; it(); ++it)
      if (DescentStatus::isDescent(desc[*it])) // then |z| is not final
      {
	switch (desc[*it])
	{
	case DescentStatus::ComplexDescent:
	  result[z]=result[block.cross(*it,z)]; // copy from unique descent
	  break;
	case DescentStatus::RealTypeII:
	  result[z]=result[block.inverseCayley(*it,z).first]; // unique descent
	  break;
	case descents::DescentStatus::RealTypeI:
	  {
	    BlockEltPair iC=block.inverseCayley(*it,z);
	    result[z]= combine(result[iC.first],result[iC.second]);
	  }
	default:
	  break; // leave |result[z]| empty in |ImaginaryCompact| case
	}
	// having found one singular descent, we ignore any other ones
	break; // not final, |break| effectively |continue|s outer loop on |z|
      }
    if (not it()) // then previous loop ran to completion
      result[z].emplace_front(z,1); // record singleton contribution to ourselves
    // the fact that |result[z].front()==z| also identifies |z| as "final"
  }
  return result;
} // |Rep_table::contributions|


// compute |extended_finialise| in |BlockElt| form, on initial range of |eblock|
std::vector<pair_list> contributions
  (const ext_block::ext_block& eblock, RankFlags singular_orbits,
   BlockElt limit) // where to stop computing contributions
{
  std::vector<pair_list> result(limit); // each|result[z]| is initially empty
  for (BlockElt z=0; z<limit; ++z)
  {
    auto s = eblock.first_descent_among(singular_orbits,z);
    if (s==eblock.rank()) // then this is an extended final element
      result[z].emplace_front(z,1); // record unit contribution to ourselves
    else
    {
      auto type=eblock.descent_type(s,z);
      if (is_like_compact(type))
	continue; // no descents, |z| represents zero; leave |result[z]| empty
      int sign = eblock.l(z,eblock.some_scent(s,z)) // true link length change
	==2 ? -1 : 1; // due to October surprise
      if (has_double_image(type)) // 1r1f, 2r11
      { auto pair = eblock.Cayleys(s,z);
	result[z] = combine
	  ( flip(sign*eblock.epsilon(s,pair.first,z),result[pair.first]),
	    flip(sign*eblock.epsilon(s,pair.second,z),result[pair.second]));
      }
      else
      { auto x = eblock.some_scent(s,z);
	result[z] = flip(sign*eblock.epsilon(s,x,z),result[x]);
      }
    }
  }
  return result;
} // |contributions|, extended block

SR_poly Rep_table::deformation_terms
  ( blocks::common_block& block, const BlockElt y,
    RankFlags singular, const RatWeight& gamma) const
{ assert(y<block.size());

  SR_poly result(repr_less());
  if (block.length(y)==0)
    return result; // easy case, null result

  std::vector<pair_list> contrib = contributions(block,singular,y);
  containers::sl_list<BlockElt> finals;
  for (BlockElt z=0; z<contrib.size(); ++z)
    if (not contrib[z].empty() and contrib[z].front().first==z)
      finals.push_front(z); // accumulate in reverse order

  const kl::KLContext& klc = block.klc(y,false); // fill silently up to |y|

  std::unique_ptr<unsigned int[]> index // a sparse array, map final to position
    (new unsigned int [block.size()]); // unlike |std::vector| do not initialise

  unsigned pos=0;
  for (auto z : finals)
    index[z]=pos++;

  // since we evaluate at $s=-1$ eventually, we can use integer coefficients
  std::vector<int> acc(finals.size(),0);
  std::vector<int> remainder(finals.size(),0); // coeff.s by |survivor| position
  remainder.front()=1; // we initialised remainder = 1*sr_y
  auto y_parity=block.length(y)%2;

  pos=0;
  // basically |for(BlockElt z:finals)|, but |pos| needs increment on |continue|
  for (auto it=finals.begin(); not finals.at_end(it); ++it,++pos)
  {
    const int c_cur = remainder[pos]; // coefficient of |z| in |remainder|
    if (c_cur==0)
      continue;
    const BlockElt z=*it; // element |pos| of |finals|; value decreases in loop
    const bool contribute = block.length(z)%2!=y_parity;
    for (BlockElt x=z+1; x-->0; ) // for |x| from |z| down to |0| inclusive
    {
      const kl::KLPol& pol = klc.klPol(x,z); // regular KL polynomial
      int eval = 0;
      for (polynomials::Degree d=pol.size(); d-->0; )
	eval = static_cast<int>(pol[d]) - eval; // evaluate at $q = -1$
      if (eval==0)
	continue; // polynomials with $-1$ as root do not contribute; skip
      if ((block.length(z)-block.length(x))%2!=0) // when |l(z)-l(x)| odd
	eval=-eval; // flip sign (do alternating sum of KL column at |-1|)
      for (auto jt=contrib[x].begin(); not contrib[x].at_end(jt); ++jt)
      {
	auto j=index[jt->first]; // position where |P(x,z)| contributes
	assert(j>=pos); // triangularity of KLV polynomials
	int c =c_cur*eval*jt->second;
	remainder[j] -= c;
	if (contribute) // optimisation will apply loop unswitching to this test
	  acc[j] += c; // here we contribute
      }
    }
    assert(remainder[pos]==0); // check relation of being inverse
  }

/* The following could be done inside the previous loop at the end of its body,
   since |acc[pos]| has its definitive value after the iteration for |pos|.
   However we keep loops separate to maybe increase locality of the one above.
   Transform coefficients of |acc| to polynomial |result|, taking into account
   the differences of |orientation_number| values between |y| and (current) |x|.
*/
  {
    const unsigned int orient_y = orientation_number(block.sr(y,gamma));

    auto it=finals.begin();
    for (const int c : acc) // accumulator |acc| runs parallel to |finals|
    {
      const auto z = *it; ++it;
      const auto sr_z = block.sr(z,gamma);

      auto coef = c*arithmetic::exp_i(orient_y-orientation_number(sr_z));
      result.add_term(sr_z,Split_integer(1,-1)*coef);
    }
    assert(it==finals.end());
  }

  return result;
} // |deformation_terms|, common block version

// compute and return sum of KL polynomials at $s$ for final parameter |sr|
SR_poly Rep_table::KL_column_at_s(StandardRepr sr) // |sr| must be final
{
  normalise(sr); // implies that |sr| it will appear at the top of its own block
  assert(is_final(sr));

  BlockElt z;
  RankFlags singular; // subset of integrally-simple roots, singular at |gamma|
  auto& block = lookup(sr,z,singular);

  std::vector<pair_list> contrib = contributions(block,singular,z);
  assert(contrib.size()==z+1 and contrib[z].front().first==z);

  const kl::KLContext& klc = block.klc(z,false); // fill silently up to |z|

  SR_poly result(repr_less());
  const auto& gamma=sr.gamma();
  auto z_length=block.length(z);
  for (BlockElt x=z+1; x-->0; )
  {
    const kl::KLPol& pol = klc.klPol(x,z); // regular KL polynomial
    if (pol.isZero())
      continue;
    Split_integer eval(0);
    for (polynomials::Degree d=pol.size(); d-->0; )
      eval.times_s() += static_cast<int>(pol[d]); // evaluate at $q = s$
    // no test here, nonzero KL polynomials have nonzero evaluation at $q=1$

    if ((z_length-block.length(x))%2!=0) // when |l(z)-l(x)| odd
      eval.negate(); // flip sign (do alternating sum of KL column at |s|)
    for (const auto& pair : contrib[x])
      result.add_term(block.sr(pair.first,gamma),eval*pair.second);
  }

  return result;
} // |Rep_table::KL_column_at_s|

#if 0
SR_poly Rep_table::deformation_terms (unsigned long sr_hash) const
{ // the |StandardRepr| |hash[sr_hash]| is necessarily final (survivor)
  SR_poly result(repr_less());
  SR_poly remainder(hash[sr_hash],repr_less());
  auto y_parity=lengths[sr_hash]%2;

  while(not remainder.empty())
  {
    auto const& leading = *remainder.cbegin(); // least term is leading term
    auto h=hash.find(leading.first); // highest term of |remainder|
    assert(h!=hash.empty); // we remain within the already tabled parameters
    auto c_cur = leading.second;
    const SR_poly& KL_cur = KLV_list[h];
    remainder.add_multiple(KL_cur,-c_cur);
    assert(remainder.empty() or hash.find(remainder.cbegin()->first)!=h);
    if (lengths[h]%2!=y_parity)
      result.add_multiple(KL_cur,c_cur);
  }
  unsigned int orient_y = orientation_number(hash[sr_hash]);
  for (auto& term : result)
  {
    unsigned int orient_express=orient_y-orientation_number(term.first);
    (term.second*= arithmetic::exp_i(orient_express)).times_1_s();
  }

  return result;
} // |deformation_terms|, version without block
#endif

SR_poly Rep_table::deformation(const StandardRepr& z)
// that |z| is dominant and final is a precondition assured in the recursion
// for more general |z|, do the preconditioning outside the recursion
{
  assert(is_final(z));
  StandardRepr z0 = z; scale_0(z0);
  SR_poly result = expand_final(z0); // value without deformation terms

  RationalList rp=reducibility_points(z); // this is OK before |make_dominant|
  if (rp.size()==0) // without deformation terms
    return result; // don't even bother to store the result

  StandardRepr z_near = z; scale(z_near,rp.back());
  normalise(z_near); // so that we may find a stored equivalent parameter
  assert(is_final(z_near));

  { // look up if deformation formula for |z_near| is already known and stored
    unsigned long h=hash.find(z_near);
    if (h!=hash.empty and not def_formulae[h].first.empty())
      return def_formulae[h].first;
  }

  // otherwise compute the deformation terms at all reducibility points
  for (unsigned i=rp.size(); i-->0; )
  {
    auto zi = z; scale(zi,rp[i]);
    normalise(zi); // necessary to ensure the following |assert| will hold
    assert(is_final(zi)); // ensures that |deformation_terms| won't refuse
    BlockElt new_z;
    auto& block = lookup(zi,new_z);

    RankFlags singular;
    const SubSystem& sub = block.integral_subsystem();
    for (weyl::Generator s=0; s<block.rank(); ++s)
      singular.set(s,root_datum().coroot(sub.parent_nr_simple(s))
					.dot(zi.gamma().numerator())==0);

    const SR_poly terms = deformation_terms(block,new_z,singular,zi.gamma());
    for (auto const& term : terms)
      result.add_multiple(deformation(term.first),term.second); // recursion
  }

  const auto h = formula_index(z_near);
  return def_formulae[h].first=result;
} // |Rep_table::deformation|


// basic computation of twisted KL column sum, no tabulation of the result
SR_poly twisted_KL_sum
( ext_block::ext_block& eblock, BlockElt y, const blocks::common_block& parent,
  const RatWeight& gamma) // infinitesimal character, possibly singular
{
  // compute cumulated KL polynomimals $P_{x,y}$ with $x\leq y$ survivors

  // start with computing KL polynomials for the entire block
  std::vector<ext_kl::Pol> pool;
  ext_kl::KL_table twisted_KLV(eblock,&pool);
  twisted_KLV.fill_columns(y+1); // fill table up to |y| inclusive

  // make a copy of |pool| in which polynomials have been evaluated as |s|
  std::vector<Split_integer> pool_at_s; pool_at_s.reserve(pool.size());
  for (unsigned i=0; i<pool.size(); ++i)
    if (pool[i].isZero())
      pool_at_s.push_back(Split_integer(0,0));
    else
    { const auto& P = pool[i];
      auto d=P.degree();
      Split_integer eval(P[d]);
      while (d-->0)
	eval.times_s()+=static_cast<int>(P[d]);
      pool_at_s.push_back(eval);
    }

  RankFlags singular_orbits;
  const auto& ipd = parent.integral_subsystem().pre_root_datum();
  for (weyl::Generator s=0; s<eblock.rank(); ++s)
    singular_orbits.set(s,gamma.dot(ipd.simple_coroot(eblock.orbit(s).s0))==0);

  auto contrib = contributions(eblock,singular_orbits,y+1);

  const auto& rc = parent.context();
  SR_poly result(rc.repr_less());
  unsigned int parity = eblock.length(y)%2;
  for (BlockElt x=0; x<=y; ++x)
  { const auto& p = twisted_KLV.KL_pol_index(x,y);
    Split_integer eval = p.second ? -pool_at_s[p.first] : pool_at_s[p.first];
    if (eblock.length(x)%2!=parity) // flip sign at odd length difference
      eval = -eval;
    for (const auto& pair : contrib[x])
      result.add_term(parent.sr(eblock.z(pair.first),gamma),eval*pair.second);
  }

  return result;
} // |twisted_KL_sum|

// compute and return sum of KL polynomials at $s$ for final parameter |z|
// since |delta| need not be the inner class involution, no storage is done
SR_poly twisted_KL_column_at_s
  (const Rep_context& rc, StandardRepr z, const WeightInvolution& delta)
  // |z| must be delta-fixed, nonzero and final
{
  rc.normalise(z);
  if (not rc.is_final(z))
    throw std::runtime_error("Parameter is not final");
  auto zm = StandardReprMod::mod_reduce(rc,z);
  BlockElt entry; // dummy needed to ensure full block is generated
  blocks::common_block block(rc,zm,entry); // which this constructor does
  ext_block::ext_block eblock(block,delta);

  return twisted_KL_sum(eblock,eblock.element(entry),block,z.gamma());
} // |twisted_KL_column_at_s|

// look up or compute and return the alternating sum of twisted KL polynomials
// at the inner class involution for final parameter |z|, evaluated at $q=s$
SR_poly Rep_table::twisted_KL_column_at_s(StandardRepr sr)
  // |z| must be inner-class-twist-fixed, nonzero and final
{
  normalise(sr);
  const auto& delta = inner_class().distinguished();
  assert(is_final(sr) and sr==inner_twisted(sr));
  BlockElt y0; RankFlags singular;
  auto& block = lookup(sr,y0,singular);
  auto& eblock = block.extended_block(delta);

  RankFlags singular_orbits;
  for (weyl::Generator s=0; s<eblock.rank(); ++s)
    singular_orbits.set(s,singular[eblock.orbit(s).s0]);

  const BlockElt y = eblock.element(y0);
  auto contrib = contributions(eblock,singular_orbits,y+1);

  containers::sl_list<BlockElt> finals; // these have numbering for |eblock|!
  for (BlockElt z=0; z<=y; ++z)
    if (not contrib[z].empty() and contrib[z].front().first==z)
      finals.push_front(z); // accumulate in reverse order

  const auto& kl_tab = eblock.kl_table(y+1);

  SR_poly result(repr_less());
  const auto& gamma=sr.gamma();
  const RatWeight gamma_rho = gamma-rho(block.root_datum());
  auto y_length=block.length(y0);

  for (BlockElt x=y+1; x-->0; )
  {
    const auto& pol = kl_tab.P(x,y); // twisted KL polynomial
    if (pol.isZero())
      continue;
    Split_integer eval(0);
    for (polynomials::Degree d=pol.size(); d-->0; )
      eval.times_s() += pol[d]; // evaluate at $q = s$

    if ((y_length-block.length(eblock.z(x)))%2!=0) // when |l(y)-l(x)| odd
      eval.negate(); // flip sign (do alternating sum of KL column at |s|)
    for (const auto& pair : contrib[x])
      result.add_term(block.sr(eblock.z(pair.first),gamma),eval*pair.second);
  }

  return result;
} // |Rep_table::twisted_KL_column_at_s|

SR_poly Rep_table::twisted_deformation_terms
    (blocks::common_block& block, ext_block::ext_block& eblock,
     BlockElt y, // in numbering of |block|, not |eblock|
     RankFlags singular_orbits, const RatWeight& gamma) const
{
  assert(eblock.is_present(y));
  const BlockElt y_index = eblock.element(y);

  SR_poly result(repr_less());
  if (block.length(y)==0)
    return result; // easy case, null result

  auto contrib = repr::contributions(eblock,singular_orbits,y_index+1);
  containers::sl_list<BlockElt> finals; // these have numbering for |eblock|!
  for (BlockElt z=0; z<contrib.size(); ++z)
    if (not contrib[z].empty() and contrib[z].front().first==z)
      finals.push_front(z); // accumulate in reverse order

  const auto& kl_tab = eblock.kl_table(y_index+1);

  std::vector<int> pool_at_minus_1; // evaluations at $q=-1$ of KL polynomials
  {
    const auto& pool=kl_tab.polys();
    pool_at_minus_1.reserve(pool.size());
    for (const auto& pol: pool)
    {
      int eval=0;
      for (unsigned i=pol.degree()+1; i-->0; )
	eval = pol[i]-eval;
      pool_at_minus_1.push_back(eval);
    }
  }

  std::unique_ptr<unsigned int[]> index // a sparse array, map final to position
    (new unsigned int [eblock.size()]); // unlike |std::vector| do not initialise

  unsigned pos=0;
  for (auto z : finals)
    index[z]=pos++;

  // since we evaluate at $s=-1$ eventually, we can use integer coefficients
  std::vector<int> acc(finals.size(),0);
  std::vector<int> remainder(finals.size(),0); // coeff.s by |survivor| position
  remainder.front()=1; // we initialised remainder = 1*sr_y
  auto y_parity=block.length(y)%2;

  pos=0;
  // basically |for(BlockElt z:finals)|, but |pos| needs increment on |continue|
  for (auto it=finals.begin(); not finals.at_end(it); ++it,++pos)
  {
    const int c_cur = remainder[pos]; // coefficient of |z| in |remainder|
    if (c_cur==0)
      continue;
    const BlockElt z=*it; // element |pos| of |finals|; value decreases in loop
    const bool contribute = block.length(eblock.z(z))%2!=y_parity;
    for (auto x : kl_tab.nonzero_column(z))
    {
      auto p = kl_tab.KL_pol_index(x,z); // pair (index,negate_p)
      if (pool_at_minus_1[p.first]==0)
	continue; // polynomials with $-1$ as root do not contribute; skip
      const int val_xz = p.second!= // XOR stored sign with length diff. parity
	((block.length(eblock.z(x))-block.length(eblock.z(z)))%2!=0)
	? -pool_at_minus_1[p.first] : pool_at_minus_1[p.first];
      for (auto jt=contrib[x].begin(); not contrib[x].at_end(jt); ++jt)
      {
	auto j=index[jt->first]; // position where |P(x,z)| contributes
	assert(j>=pos); // triangularity of KLV polynomials
	int c =c_cur*val_xz*jt->second;
	remainder[j] -= c;
	if (contribute) // optimisation will apply loop unswitching to this test
	  acc[j] += c; // here we contribute
      }
    }
    assert(remainder[pos]==0); // check relation of being inverse
  }
  {
    const unsigned int orient_y = orientation_number(block.sr(y,gamma));

    auto it=acc.begin();
    for (const int f : finals) // accumulator |acc| runs parallel to |finals|
    {
      const int c = *it++;
      if (c==0)
	continue;
      BlockElt z = eblock.z(f); // |block| numbering used to build |StandardRepr|
      const auto sr_z = block.sr(z,gamma);

      auto coef = c*arithmetic::exp_i(orient_y-orientation_number(sr_z));
      result.add_term(sr_z,Split_integer(1,-1)*coef);
    }
    assert(it==acc.end());
  }

  return result;
} // |twisted_deformation_terms(blocks::common_block&,...)|

#if 0
SR_poly Rep_table::twisted_deformation_terms (unsigned long sr_hash) const
{ // the |StandardRepr| |hash[sr_hash]| is necessarily delta-fixed and final
  SR_poly result(repr_less());
  SR_poly remainder(hash[sr_hash],repr_less());
  auto y_parity=lengths[sr_hash]%2;

  while(not remainder.empty())
  {
    auto const& leading = *remainder.cbegin(); // least term is leading term
    auto h=hash.find(leading.first); // highest term of |remainder|
    assert(h!=hash.empty); // we remain within the already tabled parameters
    auto c_cur = leading.second;
    const SR_poly& KL_cur = twisted_KLV_list[h];
    remainder.add_multiple(KL_cur,-c_cur);
    assert(remainder.empty() or hash.find(remainder.cbegin()->first)!=h);
    if (lengths[h]%2!=y_parity)
      result.add_multiple(KL_cur,c_cur);
  }
  unsigned int orient_y = orientation_number(hash[sr_hash]);
  for (auto& term : result)
  {
    unsigned int orient_express=orient_y-orientation_number(term.first);
    (term.second*= arithmetic::exp_i(orient_express)).times_1_s();
  }

  return result;
} // |twisted_deformation_terms|, version without block
#endif

SR_poly Rep_table::twisted_deformation (StandardRepr z)
{
  const auto& delta = inner_class().distinguished();
  RationalList rp=reducibility_points(z);
  bool flip_start=false; // whether a flip in descending to first point
  SR_poly result(repr_less());
  if (rp.empty())
  {
    z = ext_block::scaled_extended_dominant
	  (*this,z,delta,Rational(0,1),flip_start);
    auto L = ext_block::extended_finalise(*this,z,delta);
    for (auto it=L.begin(); it!=L.end(); ++it)
      result.add_term(it->first, it->second==flip_start
				 ? Split_integer(1,0) : Split_integer(0,1) );
    return result;
  }
  else if (rp.back()!=Rational(1,1))
  { // then shrink wrap toward $\nu=0$
    z = ext_block::scaled_extended_dominant(*this,z,delta,rp.back(),flip_start);
    Rational f=rp.back();
    for (auto& a : rp)
      a/=f; // rescale reducibility points to new parameter |z|
    assert(rp.back()==Rational(1,1)); // should make first reduction at |z|
  }

  { // if deformation for |z| was previously stored, return it with |flip_start|
    const auto h=hash.find(z);
    if (h!=hash.empty and not def_formulae[h].second.empty())
      return flip_start // if so we must multiply the stored value by $s$
	? SR_poly(repr_less()).add_multiple
	   (def_formulae[h].second,Split_integer(0,1))
	: def_formulae[h].second;
  }

  { // initialise |result| to fully deformed parameter expanded to finals
    bool flipped; // contrary to |flip_start| this affects value stored for |z|
    auto z0 = ext_block::scaled_extended_dominant
		(*this,z,delta,Rational(0,1),flipped);
    auto L = ext_block::extended_finalise(*this,z0,delta);
    for (const auto& p :L )
      result.add_term(p.first, p.second==flipped // flip means |times_s|
			       ? Split_integer(1,0) : Split_integer(0,1) );
  }

  // compute the deformation terms at all reducibility points
  for (unsigned i=rp.size(); i-->0; )
  {
    Rational r=rp[i]; bool flipped;
    auto zi = ext_block::scaled_extended_dominant(*this,z,delta,r,flipped);
    auto L =
      ext_block::extended_finalise(*this,zi,delta); // rarely a long list

    for (const auto& p : L)
    {
      BlockElt new_z; RankFlags singular;
      auto& block = lookup(p.first,new_z,singular);
      auto& eblock = block.extended_block(delta);

      RankFlags singular_orbits;
      for (weyl::Generator s=0; s<eblock.rank(); ++s)
	singular_orbits.set(s,singular[eblock.orbit(s).s0]);

      const SR_poly terms =
	twisted_deformation_terms(block,eblock,new_z,
				  singular_orbits,zi.gamma());
      const bool flip = flipped!=p.second;
      for (auto const& term : terms)
	result.add_multiple(twisted_deformation(term.first), // recursion
			    flip ? term.second.times_s() : term.second);
    }
  }

  { // store
    const auto h=formula_index(z);
    def_formulae[h].second=result;
  }

  return flip_start // if so we must multiply the stored value by $s$
    ? SR_poly(repr_less()).add_multiple(result,Split_integer(0,1))
    : result;

} // |Rep_table::twisted_deformation (StandardRepr z)|


std::ostream& Rep_context::print (std::ostream& str,const StandardRepr& z)
  const
{
  return
    str << "{x=" << z.x() << ",lambda=" << lambda(z) << ",nu=" << nu(z) << '}';
}

std::ostream& Rep_context::print (std::ostream& str,const SR_poly& P) const
{
  for (SR_poly::const_iterator it=P.begin(); it!=P.end(); ++it)
    print(str << (it==P.begin() ?"":"+") << it->second, it->first)
      << std::endl;
  return str;
}


  } // |namespace repr|
} // |namespace atlas|
