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
#include "block_minimal.h" // the |blocks::block_minimal| class
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

Rep_context::Rep_context(RealReductiveGroup &G_R)
  : G(G_R), KGB_set(G_R.kgb())
{}

size_t Rep_context::rank() const { return rootDatum().rank(); }

const TwistedInvolution Rep_context::involution_of_Cartan(size_t cn) const
{ return innerClass().involution_of_Cartan(cn); }

StandardRepr Rep_context::sr_gamma
  (KGBElt x, const Weight& lambda_rho, const RatWeight& gamma) const
{ // we use |lambda_rho| only for its real projection |(theta-1)/2*lambda_rho|
  int_Matrix theta1 = kgb().involution_matrix(x)+1;
  Weight t1_gamma (gamma.numerator().begin(), gamma.numerator().end());
  // the division in the next computation may throw when |gamma| is bad for |x|
  t1_gamma = theta1*t1_gamma/static_cast<int>(gamma.denominator());
#ifndef NDEBUG // check that constructor below builds a valid |StandardRepr|
  Weight image = // $(\theta-1)(\gamma-rho)$
    t1_gamma-(theta1*rootDatum().twoRho()/2);
  matreduc::find_solution(theta1,image); // a solution must exist
#endif
  const InvolutionTable& i_tab = innerClass().involution_table();
  return StandardRepr(x, i_tab.y_pack(kgb().inv_nr(x),lambda_rho), gamma,
		      height(t1_gamma));
}

// Height is $\max_{w\in W} \< \rho^v*w , (\theta+1)\gamma >$
unsigned int Rep_context::height(Weight theta_plus_1_gamma) const
{
  const auto& rd=rootDatum();
  int result = rd.dual_twoRho().dot(rd.make_dominant(theta_plus_1_gamma));
  assert(result>=0); assert(result%2==0);
  return static_cast<unsigned int>(result/2);
}

RatWeight Rep_context::gamma
  (KGBElt x, const Weight& lambda_rho, const RatWeight& nu) const
{
  const InvolutionTable& i_tab = innerClass().involution_table();
  const RatWeight lambda = rho(rootDatum())+lambda_rho;
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
  Weight lambda_rho = srkc.lift(srk)-rootDatum().twoRho();
  lambda_rho/=2; // undo doubled coordinates

  auto result = sr(x,lambda_rho,nu);

  // while standard, final and nonzero, |result| need not be normal
  normalise(result); // prepare |result| for direct inclusion in an |SR_poly|
  return result;
}

const WeightInvolution& Rep_context::theta (const StandardRepr& z) const
{ return innerClass().involution_table().matrix(kgb().inv_nr(z.x())); }

Weight Rep_context::lambda_rho(const StandardRepr& z) const
{
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = innerClass().involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);

  const RatWeight gamma_rho = z.gamma() - rho(rootDatum());
  Ratvec_Numer_t im_part2 = gamma_rho.numerator()+theta*gamma_rho.numerator();
  im_part2 /= gamma_rho.denominator(); // exact: $(1+\theta)(\lambda-\rho)$
  Weight i2(im_part2.begin(),im_part2.end()); // convert to |Weight|
  return (i2 + i_tab.y_lift(i_x,z.y()))/2; // division exact again
}

RatWeight Rep_context::gamma_0 (const StandardRepr& z) const
{
  const InvolutionTable& i_tab = innerClass().involution_table();
  const auto& theta = i_tab.matrix(kgb().inv_nr(z.x()));
  return ((z.gamma()+theta*z.gamma())/=2).normalize();
}

RatWeight Rep_context::nu(const StandardRepr& z) const
{
  const InvolutionTable& i_tab = innerClass().involution_table();
  const auto& theta = i_tab.matrix(kgb().inv_nr(z.x()));
  return ((z.gamma()-theta*z.gamma())/=2).normalize();
}

TorusElement Rep_context::y_as_torus_elt(const StandardRepr& z) const
{ return y_values::exp_pi(z.gamma()-lambda(z)); }

// |z| standard means (weakly) dominant on the (simple-)imaginary roots
bool Rep_context::is_standard(const StandardRepr& z, RootNbr& witness) const
{
  const RootDatum& rd = rootDatum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = innerClass().involution_table();
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
  const RootDatum& rd = rootDatum();
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
  const RootDatum& rd = rootDatum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = innerClass().involution_table();
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
  const RootDatum& rd = rootDatum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = innerClass().involution_table();
  const RootNbrSet pos_real = i_tab.real_roots(i_x) & rd.posRootSet();
  const Weight test_wt = i_tab.y_lift(i_x,z.y()) // $(1-\theta)(\lambda-\rho)$
	   + rd.twoRho()-rd.twoRho(pos_real); // replace $\rho$ by $\rho_R$

  for (RootNbrSet::iterator it=pos_real.begin(); it(); ++it)
  {
    const Weight& av = rootDatum().coroot(*it);
    if (av.dot(z.gamma().numerator())==0 and
	av.dot(test_wt)%4 !=0) // singular yet odd on shifted lambda
      return witness=*it,false;
  }
  return true;
}

bool Rep_context::is_final(const StandardRepr& z) const
{
  const RootDatum& rd = rootDatum();
  const InvolutionTable& i_tab = innerClass().involution_table();
  const auto& numer = z.gamma().numerator();
  KGBElt x=z.x();
  const InvolutionNbr i_x = kgb().inv_nr(x);

  for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
  {
    auto v=rd.simpleCoroot(s).dot(numer);
    if (v<0)
      return false;
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
  const RootDatum& rd = rootDatum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = innerClass().involution_table();
  const RootNbrSet real = innerClass().involution_table().real_roots(i_x);

  assert(real.isMember(alpha)); // only real roots should be tested

  const Weight& av = rootDatum().coroot(alpha);
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
  const RootDatum& rd = rootDatum();
  const InvolutionTable& i_tab = innerClass().involution_table();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const RootNbrSet real = i_tab.real_roots(i_x);
  const Permutation& root_inv = i_tab.root_involution(i_x);
  const Ratvec_Numer_t& numer = z.gamma().numerator();
  const arithmetic::Numer_t denom = z.gamma().denominator();
  const Weight test_wt =
    i_tab.y_lift(i_x,z.y()) +rd.twoRho() -rd.twoRho(real);

  unsigned count = 0;

  for (unsigned i=0; i<rd.numPosRoots(); ++i)
  {
    const RootNbr alpha = rd.numPosRoots()+i;
    const Weight& av = rootDatum().coroot(alpha);
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
	    and (num>0)!=(rootDatum().coroot(beta).dot(numer)>0))
	  ++count;
      }
    }
  }
  return count;
} // |orientation_number|

void Rep_context::make_dominant(StandardRepr& z) const
{
  const RootDatum& rd = rootDatum();

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
  z.y_bits=innerClass().involution_table().y_pack(kgb().inv_nr(x),lr);
} // |make_dominant|

StandardRepr&
Rep_context::singular_cross (StandardRepr& z,weyl::Generator s) const
{
  assert(rootDatum().simpleCoroot(s).dot(z.gamma().numerator())==0);
  Weight lr = lambda_rho(z); auto& x=z.x_part;
  rootDatum().simple_reflect
    (s, lr, kgb().status(s,x)==gradings::Status::Real ? 0 : 1);
  x = kgb().cross(s,x);
  z.y_bits = // reinsert $y$ bits component
    innerClass().involution_table().y_pack(kgb().inv_nr(x),lr);
  return z;
}

// auxiliary: move to a canonical for the |gens| (singular) subgroup of $W$
void Rep_context::to_singular_canonical(RankFlags gens, StandardRepr& z) const
{ // simply-singular coroots are simple, so no need to constuct a subsystem
  TwistedInvolution tw = kgb().involution(z.x_part); // copy to be modified
  WeylWord ww = innerClass().canonicalize(tw,gens);
  for (auto it=ww.begin(); it!=ww.end(); ++it) // move to that involution
    singular_cross(z,*it);
  assert(tw == kgb().involution(z.x_part));
}

void Rep_context::normalise(StandardRepr& z) const
{
  make_dominant(z);
  const RootDatum& rd = rootDatum();

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
  z.y_bits=innerClass().involution_table().y_pack(kgb().inv_nr(x),lr);
} // |normalise|

bool Rep_context::is_twist_fixed
  (StandardRepr z, const WeightInvolution& delta) const
{
  make_dominant(z);
  const RootDatum& rd = rootDatum();

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

  const RootDatum& rd = rootDatum();

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
  const RootDatum& rd = rootDatum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = innerClass().involution_table();
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
  const RootDatum& rd = rootDatum();
  const SubSystem& subsys = SubSystem::integral(rd,infin_char);
  blocks::nblock_help aux(realGroup(),subsys);
  blocks::nblock_elt src(z.x(),y_as_torus_elt(z));
  aux.cross_act(src,s);
  const RatWeight& t =	src.y().as_Qmod2Z();
  // InvolutionNbr i_x = kgb().inv_nr(z.x());
  // no need to do |innerClass().involution_table().real_unique(i_x,t)|

  RatWeight lr =  (infin_char - t - rho(rd)).normalize();
  assert(lr.denominator()==1); // we have reconstructed $\lambda-\rho \in X^*$
  return sr_gamma(src.x(),
		  Weight(lr.numerator().begin(),lr.numerator().end()), // mod 2
		  infin_char);
}

StandardRepr Rep_context::cross(const Weight& alpha, StandardRepr z) const
{
  const RootDatum& rd = rootDatum();
  KGBElt& x= z.x_part;
  InvolutionNbr i_x = kgb().inv_nr(x);
  const InvolutionTable& i_tab = innerClass().involution_table();

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
  const RootDatum& rd = rootDatum();
  const SubSystem& subsys = SubSystem::integral(rd,infin_char);
  blocks::nblock_help aux(realGroup(),subsys);
  blocks::nblock_elt src(z.x(),y_as_torus_elt(z));
  aux.do_up_Cayley(src,s);
  RatWeight t =	 src.y().log_pi(false);
  // InvolutionNbr i_x = kgb().inv_nr(z.x());
  // no need to do |innerClass().involution_table().real_unique(i_x,t)|

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
  const RootDatum& rd = rootDatum();
  const SubSystem& subsys = SubSystem::integral(rd,infin_char);
  blocks::nblock_help aux(realGroup(),subsys);
  blocks::nblock_elt src(z.x(),y_as_torus_elt(z));
  aux.do_down_Cayley(src,s);
  RatWeight t =	 src.y().log_pi(false);
  // InvolutionNbr i_x = kgb().inv_nr(z.x());
  // no need to do |innerClass().involution_table().real_unique(i_x,t)|

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
  Weight sum(rd.rank(),0); // difference of $\rho_r$ values
  for (auto it=S.begin(); it(); ++it)
    sum += rd.root(*it); // sum real posroots upstairs that |to_simple| negates
  return sum;
}

// a method used to ensure |z| is integrally dominant, used by |any_Cayley|
WeylWord
Rep_context::make_dominant(StandardRepr& z,const SubSystem& subsys) const
{
  const RootDatum& rd = rootDatum();
  KGBElt& x= z.x_part;
  InvolutionNbr i_x = kgb().inv_nr(x);
  const InvolutionTable& i_tab = innerClass().involution_table();

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
  const RootDatum& rd = rootDatum();
  const KGB& kgb = this->kgb();
  const InvolutionTable& i_tab = innerClass().involution_table();
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
    Cayley_shift(innerClass(),ascent ? kgb.inv_nr(x) : inv0,ww);
  z = sr_gamma(x,lr,infin_char);

  return z;
} // |Rep_context::any_Cayley|

StandardRepr Rep_context::inner_twisted(StandardRepr z) const
{
  make_dominant(z);
  RatWeight gamma=z.gamma(); // now get the infinitesimal character
  const RootDatum& rd = rootDatum();
  const SubSystem& subsys = SubSystem::integral(rd,gamma);
  blocks::nblock_help aux(realGroup(),subsys);
  blocks::nblock_elt src(z.x(),y_as_torus_elt(z));

  innerClass().distinguished().apply_to(gamma.numerator()); // twist |gamma|
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
  const auto& i_tab = innerClass().involution_table();
  const InvolutionNbr i_x0 = kgb().inv_nr(z.x());
  z.x_part = kgb().twisted(z.x_part,delta);
  const InvolutionNbr i_x1 =  kgb().inv_nr(z.x()); // destination involution
  z.y_bits = i_tab.y_act(i_x0,i_x1,z.y_bits,delta);
  z.infinitesimal_char = delta*z.infinitesimal_char;
  return z;
}

Rep_context::compare Rep_context::repr_less() const
{ return compare(rootDatum().dual_twoRho()); }

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

unsigned int Rep_table::length(StandardRepr z)
{
  make_dominant(z); // should't hurt, and improves chances of finding |z|
  unsigned long hash_index=hash.find(z);
  if (hash_index!=hash.empty)
    return lengths[hash_index];

  // otherwise do it the hard way, constructing a block up to |z|
  param_block block(*this,z); // compute partial block
  return block.length(block.size()-1);
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
  const RootDatum& rd = rootDatum();
  const InvolutionTable& i_tab = innerClass().involution_table();

  make_dominant(z);
  const RatWeight& gamma=z.gamma();

  RankFlags singular;
  for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
    singular.set(s,rd.simpleCoroot(s).dot(gamma.numerator())==0);

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
	  singular_cross(*rit,s);
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

void Rep_table::add_block
  (param_block& block, containers::sl_list<BlockElt>& survivors)
{
  for (BlockElt x=0; x<block.size(); ++x)
    if (block.survives(x))
      survivors.push_back(x);

  unsigned long old_size = hash.size();
  BlockEltList new_survivors;

  // fill the |hash| table for new surviving parameters in this block
  for (auto x : survivors)
    if (hash.match(block.sr(x))>=old_size)
      new_survivors.push_back(x);

  assert(hash.size()==old_size+new_survivors.size()); // only new surv. added
  if (new_survivors.empty())
    return; // nothing left to do, but we have computed |survivors| for caller

  lengths.resize(hash.size());
  KLV_list.resize(hash.size(),SR_poly(repr_less())); // new slots, init empty
  def_formula.resize(hash.size(),SR_poly(repr_less())); // allocate new slots

  // compute cumulated KL polynomimals $P_{x,y}$ with $x\leq y$ survivors

  // start with computing KL polynomials for the entire block
  const kl::KLContext& klc = block.klc(block.size()-1,false); // silently

  /* get $P(x,z)$ for |x<=z| with |z| among new |survivors|, and contribute
   parameters from |block.finals_for(x)| with coefficient $P(x,z)[q:=s]$
   to the |SR_poly| at |KLV_list[old_size+i], where |z=new_survivors[i]| */

  auto z_start=new_survivors.cbegin();
  for (auto it = z_start; it!=new_survivors.cend(); ++it)
    lengths[old_size+(it-z_start)]=block.length(*it);

  for (BlockElt x=0; x<=new_survivors.back(); ++x)
  {
    auto xs=block.finals_for(x);
    if (xs.empty())
      continue; // no point doing work for |x|'s that don't contribute anywhere

    const unsigned int parity = block.length(x)%2;

    if (*z_start<x)
      ++z_start; // advance so |z| only runs over values with |x<=z|
    assert(z_start!=new_survivors.end() and *z_start>=x);

    for (auto it=z_start; it!=new_survivors.end(); ++it)
    {
      const BlockElt z = *it; // element of |new_survivors| and |x<=z|
      const kl::KLPol& pol = klc.klPol(x,z); // regular KL polynomial
      Split_integer eval(0);
      for (polynomials::Degree d=pol.size(); d-->0; )
	eval.times_s()+=static_cast<int>(pol[d]);
      if (not eval.is_zero())
      {
	auto z_index = old_size+(it-new_survivors.begin());
	assert(hash.find(block.sr(z))==z_index);
	SR_poly& dest = KLV_list[z_index]; // a poly to which |x| contributes
	if (lengths[z_index]%2!=parity)
	  eval.negate(); // incorporate sign for length difference
	for (BlockElt x : xs) // add from |klPol(x,z)| for all |finals_for(x)|
	  dest.add_term(block.sr(x),eval); // contribute term |eval| for each
      }
    } // |for(it)|
  } // |for(x)|
} // |Rep_table::add_block|

// compute and return sum of KL polynomials at $s$ for final parameter |z|
SR_poly Rep_table::KL_column_at_s(StandardRepr z) // |z| must be final
{
  normalise(z); // implies that |z| it will appear at the top of its own block
  assert(is_final(z));
  unsigned long hash_index=hash.find(z);
  if (hash_index==hash.empty) // previously unknown parameter
  { // then we need to compute to find the requested polynomial
    param_block block(*this,z);
    containers::sl_list<BlockElt> survivors;
    add_block(block,survivors);

    hash_index=hash.find(z);
    assert(hash_index!=hash.empty);
  }

  return KLV_list[hash_index];
} // |Rep_table::KL_column_at_s|

SR_poly Rep_table::deformation_terms (param_block& block,BlockElt y)
{
  SR_poly result(repr_less());
  if (not block.survives(y) or block.length(y)==0)
    return result; // easy cases, null result

  containers::sl_list<BlockElt> survivors;
  add_block(block,survivors); // computes survivors, and add anything new

  // map indices of |survivors| to corresponding number in |hash|
  auto sr_y = block.sr(y);
  auto y_parity=block.length(y)%2;
  assert(hash.find(sr_y)!=hash.empty); // |sr_y| should be known now

  { // cut off |survivors| after value |y|, and reverse what is left
    auto it=survivors.cbegin();
    while (*it<y) // skip over everything begore the value |y|
      ++it;
    ++it; // and skip that value itself
    survivors.erase(it,survivors.cend());
    survivors.reverse(); // now |survivors| is decreasing list starting with |y|
  }

  std::vector<unsigned long> remap; // map index into |survivors| to hash nr
  remap.reserve(survivors.size());
  std::unique_ptr<unsigned int[]> index // a sparse array inverting |remap|
    (new unsigned int [hash.size()]); // unlike |std::vector| do not initialise
  for (auto x : survivors)
  {
    auto h=hash.find(block.sr(x));
    assert(h!=hash.empty); // all standard reps for survivors should be found
    index[h]=remap.size(); // point back to entry pushed in next line
    remap.push_back(h);
  }

  // since we evaluate at $s=-1$ eventually, we can use integer coefficients
  std::vector<int> acc(remap.size(),0);
  std::vector<int> remainder(remap.size(),0); // coeff.s by |survivor| position
  remainder.front()=1; // we initialised remainder = 1*sr_y

  for (auto cur : survivors) // value of |cur| is, and must be, decreasing here
  {
    StandardRepr p_cur=block.sr(cur);
    auto h = hash.find(p_cur);
    assert(h!=hash.empty);
    unsigned i=index[h];
    int c_cur = remainder[i];
    const SR_poly& KL_cur = KLV_list[h];
    const bool contribute = block.length(cur)%2!=y_parity;
    for (const auto& pair : KL_cur)
    { auto j = index[hash.find(pair.first)];
      assert(j>=i); // triangularity of |KLV_list|
      int c = c_cur*pair.second.s_to_minus_1();
      remainder[j] -= c;
      if (contribute) // optimisation will apply loop unswitching to this test
	acc[j] += c; // here we contribute
    }
    assert(remainder[i]==0); // check relation of being inverse
  }

  unsigned int orient_ee = orientation_number(sr_y);
  unsigned i=0;
  for (auto x : survivors)
  { auto sr_x=block.sr(x);
    unsigned int orient_express=orient_ee-orientation_number(sr_x);
    auto coef = acc[i++]*arithmetic::exp_i(orient_express);
    result.add_term(sr_x,Split_integer(1,-1)*coef);
  }

  return result;
} // |deformation_terms|

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

  { // look up if closest reducibility point to |z| is already known
    unsigned long h=hash.find(z_near);
    if (h!=hash.empty and not def_formula[h].empty())
      return def_formula[h];
  }

  // otherwise compute the deformation terms at all reducibility points
  for (unsigned i=rp.size(); i-->0; )
  {
    StandardRepr zi = z; scale(zi,rp[i]);
    normalise(zi); // necessary to ensure the following |assert| will hold
    assert(is_final(zi)); // ensures that |deformation_terms| won't refuse
    param_block b(*this,zi); // construct block interval below |zi|
    const SR_poly terms = deformation_terms(b,b.size()-1);
    for (SR_poly::const_iterator it=terms.begin(); it!=terms.end(); ++it)
      result.add_multiple(deformation(it->first),it->second); // recursion
  }

  // now store result for future lookup
  unsigned long h=hash.find(z_near);
  assert(h!=hash.empty); // it should have been added by |deformation_terms|
  def_formula[h]=result;

  return result;
} // |Rep_table::deformation|


// basic computation of twisted KL column sum, no tabulation of the result
SR_poly twisted_KL_sum
( const Rep_context& rc, ext_block::ext_block& eblock, BlockElt y,
  param_block& parent) // its complete unextended block
{
  // compute cumulated KL polynomimals $P_{x,y}$ with $x\leq y$ survivors

  // start with computing KL polynomials for the entire block
  std::vector<ext_kl::Pol> pool;
  ext_kl::KL_table twisted_KLV(eblock,pool);
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
	eval = eval.times_s()+Split_integer(P[d]);
      pool_at_s.push_back(eval);
    }

  // construct a one-column matrix $(P_{x,y}[q:=s])_{x,0}$, range $x$ is block
  matrix::Matrix<Split_integer> P_at_s(y+1,1);
  for (BlockElt x=0; x<=y; ++x)
  { auto pair = twisted_KLV.KL_pol_index(x,y);
    P_at_s(x,0) = // get value from |pool_at_s|, possibly negated
      pair.second ? -pool_at_s[pair.first] : pool_at_s[pair.first];
  }

  // condense |P_at_s| to the extended block elements without singular descents
  containers::simple_list<BlockElt> survivors = eblock.condense(P_at_s,parent);

  // finally transcribe from |P_at_s| result
  SR_poly result(rc.repr_less());
  unsigned int parity = eblock.length(y)%2;
  for (auto it = survivors.begin(); not survivors.at_end(it); ++it)
  {
    BlockElt x = *it;
    auto factor = P_at_s(x,0);
    if (eblock.length(x)%2!=parity) // flip sign at odd length difference
      factor = -factor;
    result.add_term(parent.sr(eblock.z(x)),factor);
  }
  return result;
} // |twisted_KL_sum|

// same computation of twisted KL column sum, but with a |block_minimal|
SR_poly twisted_KL_sum
( ext_block::ext_block& eblock, BlockElt y, const blocks::block_minimal& parent,
  const RatWeight& gamma) // infinitesimal character, possibly singular
{
  // compute cumulated KL polynomimals $P_{x,y}$ with $x\leq y$ survivors

  // start with computing KL polynomials for the entire block
  std::vector<ext_kl::Pol> pool;
  ext_kl::KL_table twisted_KLV(eblock,pool);
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
	eval = eval.times_s()+Split_integer(P[d]);
      pool_at_s.push_back(eval);
    }

  // construct a one-column matrix $(P_{x,y}[q:=s])_{x,0}$, range $x$ is block
  matrix::Matrix<Split_integer> P_at_s(y+1,1);
  for (BlockElt x=0; x<=y; ++x)
  { auto pair = twisted_KLV.KL_pol_index(x,y);
    P_at_s(x,0) = // get value from |pool_at_s|, possibly negated
      pair.second ? -pool_at_s[pair.first] : pool_at_s[pair.first];
  }

  // condense |P_at_s| to the extended block elements without singular descents
  containers::simple_list<BlockElt> survivors =
    eblock.condense(P_at_s,parent,gamma);

  // finally transcribe from |P_at_s| result
  const auto& rc = parent.context();
  const auto gamma_rho = gamma-rho(parent.rootDatum());
  SR_poly result(rc.repr_less());
  unsigned int parity = eblock.length(y)%2;
  for (auto it = survivors.begin(); not survivors.at_end(it); ++it)
  {
    BlockElt elt = *it;
    BlockElt z = eblock.z(elt); // index of |elt| in |parent|
    auto factor = P_at_s(elt,0);
    if (eblock.length(elt)%2!=parity) // flip sign at odd length difference
      factor = -factor;
    const auto lambda_rho = gamma_rho.integer_diff<int>(parent.gamma_lambda(z));
    result.add_term(rc.sr_gamma(parent.x(z),lambda_rho,gamma),factor);
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
  BlockElt entry; // dummy needed to ensure full block is generated
  blocks::block_minimal block(rc,z,entry); // which this constructor does
  ext_block::ext_block eblock(block,delta);

  return twisted_KL_sum(eblock,eblock.element(entry),block,z.gamma());
} // |twisted_KL_column_at_s|

void Rep_table::add_block(ext_block::ext_block& block, // a full extended block
			  param_block& parent, // its complete unextended block
			  BlockElt top_elt, // |block| only used for |<=top_elt|
			  containers::sl_list<BlockElt>& extended_finals)
{
  { containers::sl_list<BlockElt> dummy;  // this exported value will not be used
    add_block(parent,dummy); // but we must ensure parent block is known
  }

  RankFlags singular_orbits = block.singular_orbits(parent);
  ndebug_use(singular_orbits); // it is only used in |assert| statements

  // extend space in twisted tables; zero polynomial means no computed value
  twisted_KLV_list.resize(hash.size(),SR_poly(repr_less())); // init empties
  twisted_def_formula.resize(hash.size(),SR_poly(repr_less())); // init empties

  // compute cumulated KL polynomimals $P_{x,y}$ with $x\leq y$ survivors

  // start with computing KL polynomials for the entire block
  std::vector<ext_kl::Pol> pool;
  ext_kl::KL_table twisted_KLV(block,pool);
  twisted_KLV.fill_columns(); // block is complete, so fill everything

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
	eval = eval.times_s()+Split_integer(P[d]);
      pool_at_s.push_back(eval);
    }

  // construct an upper triangular matrix $(P_{x,y}[q:=s])_{x,y}$ over block
  matrix::Matrix<Split_integer> P_at_s(block.size()); // initialise: identity
  for (BlockElt x=0; x<block.size(); ++x)
    for (BlockElt y=x+1; y<block.size(); ++y)
    { auto pair = twisted_KLV.KL_pol_index(x,y);
      P_at_s(x,y) = // get value from |pool_at_s|, possibly negated
	pair.second ? -pool_at_s[pair.first] : pool_at_s[pair.first];
    }

  // condense |P_mat| to the extended block elements without singular descents
  extended_finals = containers::sl_list<BlockElt> // convert from |simple_list|
    { block.condense(P_at_s,parent) };

  auto top_elt_it = extended_finals.cend();
  // now |at_end(top_elt_it)|, but |top_elt_it| set to "point to" |top_elt| below

  // finally transcribe columns from |P_at_s| into |twisted_KLV_list|
  for (auto it=extended_finals.cbegin(); it!=extended_finals.cend(); ++it)
  { assert(block.first_descent_among(singular_orbits,*it)==block.rank());
    auto y = *it; // block element, index into (extended) |block|
    if (block.z(y)==top_elt)
      top_elt_it=it;
    auto y_index = hash.find(block.sr(y,parent));
    assert (y_index!=hash.empty); // since we looked up everything above
    assert (y_index<twisted_KLV_list.size());
    SR_poly& dest = twisted_KLV_list[y_index];
    if (dest.empty()) // this means the entry was never defined
    { dest = SR_poly(block.sr(y,parent),repr_less()); // coefficient 1
      unsigned int parity = block.length(y)%2;
      for (auto x_it=extended_finals.begin(); x_it!=it; ++x_it) // upper part
      {
	auto x = *x_it;
	auto factor = P_at_s(x,y);
	if (block.length(x)%2!=parity) // flip sign at odd length difference
	  factor = -factor;
	dest.add_term(block.sr(x,parent),factor);
      }
      // since |dest| is a reference, the sum is stored at its destination
    } // |if (y_index>=old_size)|
  } // |for (y : extended_finals)|
  assert(not top_elt_it.at_end()); // we must have set it somewhere on our way
  extended_finals.erase(++top_elt_it,extended_finals.cend()); // keep to |top_elt|
} // |Rep_table::add_block| (extended block)

// look up or compute and return the alternating sum of twisted KL polynomials
// at the inner class involution for final parameter |z|, evaluated at $q=s$
SR_poly Rep_table::twisted_KL_column_at_s(StandardRepr z)
  // |z| must be inner-class-twist-fixed, nonzero and final
{
  normalise(z);
  if (not is_final(z))
    throw std::runtime_error("Parameter not final");
  if (z!=inner_twisted(z))
    throw std::runtime_error("Parameter not twist-fixed");

  unsigned long hash_index=hash.find(z);
  if (hash_index>=twisted_KLV_list.size() // |z| unknown or not extended to, or
      or twisted_KLV_list[hash_index].empty()) // slot created by another block
  {
    BlockElt entry; // dummy needed to ensure full block is generated
    param_block block(*this,z,entry); // which this constructor does
    ext_block::ext_block eblock(block,innerClass().distinguished());

    containers::sl_list<BlockElt> extended_finals;
    add_block(eblock,block,entry,extended_finals);

    if (hash_index==hash.empty) // then reset |hash_index| to now known value
    { hash_index=hash.find(z);
      assert(hash_index!=hash.empty);
    }
  }

  return twisted_KLV_list[hash_index];
} // |Rep_table::twisted_KL_column_at_s|

SR_poly Rep_table::twisted_deformation_terms (param_block& parent,BlockElt y)
{
  const auto& delta = innerClass().distinguished();
  const auto sr_y = parent.sr(y);

  assert(is_twist_fixed(sr_y,delta));

  SR_poly result(repr_less());
  if (not parent.survives(y) or parent.length(y)==0)
    return result; // easy cases, null result

  ext_block::ext_block eblock(parent,delta);
  containers::sl_list<BlockElt> extended_finals;
  add_block(eblock,parent,y,extended_finals);
  assert(eblock.is_present(y)); // since |is_twist_fixed| succeeded
  const BlockElt y_index = eblock.element(y);

  assert(hash.find(sr_y)!=hash.empty); // |sr_y| should be known now
  unsigned int y_parity = eblock.length(y_index)%2;

  extended_finals.reverse();  // now |extended_finals| decreasing from |y|

  std::vector<unsigned long> remap; // from |extended_finals| index to hash nr
  remap.reserve(extended_finals.size());
  std::unique_ptr<unsigned int[]> index // a sparse array inverting |remap|
    (new unsigned int [hash.size()]); // unlike |std::vector| do not initialise
  for (auto x : extended_finals)
  {
    unsigned long h=hash.find(eblock.sr(x,parent));
    assert(h!=hash.empty); // all |extended_finals| standardreps should be found
    index[h]=remap.size(); // point back to entry pushed in next line
    remap.push_back(h);
  }

  // since we evaluate at $s=-1$ eventually, we can use integer coefficients
  std::vector<int> acc(remap.size(),0);
  std::vector<int> remainder(remap.size(),0);
  remainder.front()=1; // we initialised remainder = 1*sr_y

   for (auto cur : extended_finals) // value of |cur| is, and must be, decreasing
  {
    StandardRepr p_cur=eblock.sr(cur,parent);
    auto h = hash.find(p_cur);
    assert(h!=hash.empty);
    unsigned i=index[h];
    int c_cur = remainder[i];
    const SR_poly& KL_cur = twisted_KLV_list[h];
    const bool contribute = eblock.length(cur)%2!=y_parity;
    for (const auto& pair : KL_cur)
    { auto j = index[hash.find(pair.first)];
      assert(j>=i); // triangularity of |twisted_KLV_list|
      int c = c_cur*pair.second.s_to_minus_1();
      remainder[j] -= c;
      if (contribute) // optimisation will apply loop unswitching to this test
	acc[j] += c; // here we contribute
    }
    assert(remainder[i]==0); // check relation of being inverse
  }

  // correct signs in terms of result according to orientation numbers
  unsigned int orient_y = orientation_number(sr_y);
  unsigned i=0;
  for (auto x : extended_finals)
  { auto sr_x=eblock.sr(x,parent);
    unsigned int orient_express=orient_y-orientation_number(sr_x);
    auto coef = acc[i++]*arithmetic::exp_i(orient_express);
    result.add_term(sr_x,Split_integer(1,-1)*coef);
  }

  return result;
} // |twisted_deformation_terms|

// the next function is recursive, so avoid testing properties each time
// notably assure |z| is final and inner-twist |fixed| before calling this
SR_poly Rep_table::twisted_deformation (StandardRepr z)
{
  const auto& delta = innerClass().distinguished();
  RationalList rp=reducibility_points(z);
  bool flip_start=false; // whether a flip in descending to first point
  if (rp.empty()) // then the interesting point is the deformation to $\nu=0$
    z = ext_block::scaled_extended_dominant
	  (*this,z,delta,Rational(0,1),flip_start);
  else
  { if (rp.back()!=Rational(1,1)) // then shrink wrap toward $\nu=0$
    { z =
	ext_block::scaled_extended_dominant(*this,z,delta,rp.back(),flip_start);
      Rational f=rp.back();
      for (auto it=rp.begin(); it!=rp.end(); ++it)
	(*it)/=f; // rescale reducibility points to new parameter |z|
      assert(rp.back()==Rational(1,1)); // should make first reduction at |z|
    }

    // now (still with |not rp.empty()| check if a result was previously stored
    unsigned long h=hash.find(z);
    if (h<twisted_def_formula.size() and not twisted_def_formula[h].empty())
      return flip_start // if so we must multiply the stored value by $s$
	? SR_poly(repr_less()) // need an empty polynomial here
	  .add_multiple(twisted_def_formula[h],Split_integer(0,1))
	: twisted_def_formula[h];
  }

  // this is the first time for |z|, so we must compute
  SR_poly result(repr_less());
  { // initialise |result| to fully deformed parameter expanded to finals
    bool flipped;
    auto z0 = ext_block::scaled_extended_dominant
		(*this,z,delta,Rational(0,1),flipped);
    auto L = ext_block::extended_finalise(*this,z0,delta);
    for (auto it=L.begin(); it!=L.end(); ++it)
      result.add_term(it->first, it->second==flipped
				 ? Split_integer(1,0) : Split_integer(0,1) );
  }

  if (not rp.empty()) // without reducuibilty points, just return |result| now
  {
    BlockElt dummy;
    param_block parent(*this,z,dummy); // full parent block needed for now
    ext_block::ext_block eblock(parent,delta); // full as well
    containers::sl_list<BlockElt> extended_finals;
    add_block(eblock,parent,dummy,extended_finals);
    const unsigned long h=hash.find(z);
    assert(h<twisted_def_formula.size()); // it was just added by |add_block|

    for (unsigned i=rp.size(); i-->0; )
    {
      Rational r=rp[i]; bool flipped;
      auto zi = ext_block::scaled_extended_dominant(*this,z,delta,r,flipped);
      std::unique_ptr<param_block> bp;
      if (i+1<rp.size()) // avoid regenerating same parent block first time
        bp.reset(new param_block(*this,zi,dummy));
      param_block& block = bp.get()==nullptr ? parent : *bp;

      auto L =
	ext_block::extended_finalise(*this,zi,delta); // rarely a long list
      for (auto it=L.begin(); it!=L.end(); ++it)
      { auto zz = block.lookup(it->first);
	SR_poly terms = twisted_deformation_terms(block,zz);
	for (SR_poly::iterator jt=terms.begin(); jt!=terms.end(); ++jt)
	  result.add_multiple(twisted_deformation(jt->first),  // recursion
		      flipped==it->second ? jt->second : jt->second.times_s());
      }
    }
    twisted_def_formula[h]=result; // now store result for future lookup

  }

  return flip_start // if so we must multiply the stored value by $s$
    ? SR_poly(repr_less()).add_multiple(result,Split_integer(0,1))
    : result;
} // |Rep_table::twisted_deformation|

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
