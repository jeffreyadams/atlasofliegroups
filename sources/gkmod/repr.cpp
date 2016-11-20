/*
  This is repr.cpp

  Copyright (C) 2009-2012 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "repr.h"

#include <map> // used in computing |reducibility_points|
#include <iostream>
#include "error.h"

#include "arithmetic.h"
#include "matreduc.h"

#include "tits.h"
#include "kgb.h"	// various methods
#include "blocks.h"	// |dual_involution|
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
{
  const InvolutionTable& i_tab = innerClass().involution_table();
#ifndef NDEBUG // check that constructor below builds a valid StandardRepr
  int_Matrix theta1 = kgb().involution_matrix(x)+1;
  RatWeight g_r = gamma - rho(rootDatum());
  Weight image (g_r.numerator().begin(),g_r.numerator().end()); // convert
  // |gamma| is compatible with |x| if neither of next two lines throws
  image = theta1*image/int(g_r.denominator()); // division must be exact
  matreduc::find_solution(theta1,image); // solution must exist
#endif
  return StandardRepr(x, i_tab.y_pack(kgb().inv_nr(x),lambda_rho), gamma);
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

  return sr(x,lambda_rho,nu);
}

StandardRepr Rep_context::sr(const param_block& b, BlockElt i) const
{
  assert(i<b.size());
  return sr_gamma(b.x(i),b.lambda_rho(i),b.gamma());
}

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

// return $\lambda \in \rho+X^*$ as half-integer rational vector
RatWeight Rep_context::lambda(const StandardRepr& z) const
{
  RatWeight result(rho(rootDatum()));
  return result.normalize()+lambda_rho(z);
}

RatWeight Rep_context::nu(const StandardRepr& z) const
{
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const WeightInvolution& theta = innerClass().involution_table().matrix(i_x);
  const Ratvec_Numer_t num = z.gamma().numerator()-theta*z.gamma().numerator();
  return RatWeight(num,2*z.gamma().denominator()).normalize();
}

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

// |z| zero means that no singular simple-imaginary roots are compact; this
// code assumes |is_standard(z)|, namely |gamma| is dominant on imaginary roots
bool Rep_context::is_zero(const StandardRepr& z, RootNbr& witness) const
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
      return witness=alpha,true;
  }
  return false;
}


// |z| final means that no singular real roots satisfy the parity condition
// we do not assume |gamma| to be dominant, so all real roots must be tested
bool Rep_context::is_final(const StandardRepr& z, RootNbr& witness) const
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

bool Rep_context::is_oriented(const StandardRepr& z, RootNbr alpha) const
{
  const RootDatum& rd = rootDatum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = innerClass().involution_table();
  const RootNbrSet real = innerClass().involution_table().real_roots(i_x);

  assert(real.isMember(alpha)); // only real roots should be tested

  const Weight& av = rootDatum().coroot(alpha);
  const int numer = av.dot(z.gamma().numerator());
  const int denom = z.gamma().denominator();
  assert(numer%denom!=0); // and the real root alpha should be non-integral

  const Weight test_wt =
    i_tab.y_lift(i_x,z.y()) +rd.twoRho() -rd.twoRho(real);
  const int eps = av.dot(test_wt)%4==0 ? 0 : denom;

  return arithmetic::remainder(numer+eps,2*denom)< (unsigned)denom;
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

void Rep_context::W_act(const WeylWord& w,StandardRepr& z) const
{
  const RootDatum& rd = rootDatum();
  Weight lr = lambda_rho(z);
  KGBElt& x = z.x_part;
  Ratvec_Numer_t& numer = z.infinitesimal_char.numerator();

  for (unsigned i=w.size(); i-->0; )
  {
    weyl::Generator s=w[i];
    rd.simple_reflect(s,numer);
    rd.simple_reflect(s,lr);
    if (kgb().status(s,x)!=gradings::Status::Real) // center at $\rho-\rho_r$
      lr -= rd.simpleRoot(s); // so unless |s| is real root compensate
    x = kgb().cross(s,x);
  }
  z.y_bits = // reinsert $y$ bits component
    innerClass().involution_table().y_pack(kgb().inv_nr(x),lr);
}

void
Rep_context::W_act(const WeylWord& w,StandardRepr& z,const SubSystem& subsys)
  const
{
  const RootDatum& rd = rootDatum();
  KGBElt& x= z.x_part;
  InvolutionNbr i_x = kgb().inv_nr(x);
  const InvolutionTable& i_tab = innerClass().involution_table();

  // the following are non-|const|, and modified in the loop below
  Weight lambda2_shifted = (lambda_rho(z)*=2)
    + rd.twoRho() - rd.twoRho(i_tab.real_roots(i_x));
  Ratvec_Numer_t& gamma_num = z.infinitesimal_char.numerator();

  for (unsigned i=w.size(); i-->0; )
  {
    weyl::Generator s=w[i];
    RootNbr alpha = subsys.parent_nr_simple(s);
    rd.reflect(alpha,gamma_num);
    x = kgb().cross(rd.reflectionWord(alpha),x);
    i_x = kgb().inv_nr(x);
    rd.reflect(alpha,lambda2_shifted);
  }
  lambda2_shifted -= rd.twoRho() - rd.twoRho(i_tab.real_roots(i_x)); // unshift
  z.y_bits=i_tab.y_pack(i_x,lambda2_shifted/2);
}

WeylWord Rep_context::make_dominant(StandardRepr& z) const
{
  const RootDatum& rd = rootDatum();

  // the following are non-|const|, and modified in the loop below
  Weight lr = lambda_rho(z);
  KGBElt& x = z.x_part;
  Ratvec_Numer_t& numer = z.infinitesimal_char.numerator();

  WeylWord result;
  result.reserve(rd.numPosRoots()); // enough to accommodate the WeylWord

  { weyl::Generator s;
    do
      for (s=0; s<rd.semisimpleRank(); ++s)
      {
	int v=rd.simpleCoroot(s).dot(numer);
        if (v<0 or (v==0 and kgb().isComplexDescent(s,x)))
        {
	  result.push_back(s);
          rd.simple_reflect(s,numer);
          rd.simple_reflect(s,lr);
	  switch (kgb().status(s,x))
	  {
	  case gradings::Status::ImaginaryCompact:
	  case gradings::Status::ImaginaryNoncompact:
	    throw std::runtime_error("Non standard parameter in make_dominant");
	  case gradings::Status::Complex:
	    lr -= rd.simpleRoot(s); // pivot around $\rho-\rho_r$
	  case gradings::Status::Real: {} // no compensation for real roots
	  }
          x = kgb().cross(s,x);
	  break; // out of the loop |for(s)|
        } // |if(v<0)|
      } // |for(s)|
    while (s<rd.semisimpleRank()); // wait until inner loop runs to completion
  }
  z.y_bits=innerClass().involution_table().y_pack(kgb().inv_nr(x),lr);
  return result;
} // |make_dominant|

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
        if (v<0 or (v==0 and i_tab.is_complex_descent(i_x,alpha)))
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

  // we shall associate to a number $num>0$ a strict lower bound for $k$,
  // for which we shall then later form fractions $(d/num)*k$
  typedef std::map<long,long> table;

  // because of the parity condition, distinguish cases with even and odd $k$
  table odds,evens; // name indicates the parity that $k$ will have

  for (RootNbrSet::iterator it=pos_real.begin(); it(); ++it)
  {
    arithmetic::Numer_t num =
      rd.coroot(*it).dot(numer); // numerator of $\<\alpha^v,\nu>$ in real case
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
    arithmetic::Numer_t num = vala - valb; // numerator of $2\<\alpha^v,\nu>$
    if (num!=0)
    {
      assert((vala+valb)%d==0); // since |\<a+b,\gamma>=\<a+b,\lambda>|
      long lwb =std::abs(vala+valb)/d;
      std::pair<table::iterator,bool> trial =
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
  blocks::nblock_elt src(z.x(),y_values::exp_pi(infin_char-lambda(z)));
  aux.cross_act(src,s);
  const RatWeight& t =  src.y().as_Qmod2Z();
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
  blocks::nblock_elt src(z.x(),y_values::exp_pi(infin_char-lambda(z)));
  aux.do_up_Cayley(src,s);
  RatWeight t =  src.y().log_pi(false);
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
  blocks::nblock_elt src(z.x(),y_values::exp_pi(infin_char-lambda(z)));
  aux.do_down_Cayley(src,s);
  RatWeight t =  src.y().log_pi(false);
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
  if (rt == ~ RootNbr(0)) // either not a root at all or not in subsystem
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

  W_act(w,z); // move back to origingal infinitesimal character representative
  return z;
}

StandardRepr Rep_context::inner_twisted(StandardRepr z) const
{
  make_dominant(z);
  const RatWeight infin_char=z.gamma(); // now get the infinitesimal character
  const RootDatum& rd = rootDatum();
  const SubSystem& subsys = SubSystem::integral(rd,infin_char);
  blocks::nblock_help aux(realGroup(),subsys);
  blocks::nblock_elt src(z.x(),y_values::exp_pi(infin_char-lambda(z)));
  aux.twist(src);
  RatWeight lr =
    (infin_char - src.y().log_pi(false) - rho(rd)).normalize();
  assert(lr.denominator()==1);
  return sr_gamma(src.x(),
		  Weight(lr.numerator().begin(),lr.numerator().end()), // mod 2
		  infin_char);
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
  if (r.x()!=s.x()) // order by |x| component first
    return r.x()<s.x();

  // then compare by scalar product of |gamma()| and |level_vec|
  if (r.gamma()!=s.gamma()) // quick test to avoid work within a same block
  {
    const int rgd=r.gamma().denominator(), sgd=s.gamma().denominator();
    const int lr = sgd*level_vec.dot(r.gamma().numerator()); // cross multiply
    const int ls = rgd*level_vec.dot(s.gamma().numerator());
    if (lr!=ls)
      return lr<ls;

    // next by individual components of |gamma()|
    for (size_t i=0; i<level_vec.size(); ++i)
      if (sgd*r.gamma().numerator()[i]!=rgd*s.gamma().numerator()[i])
	return sgd*r.gamma().numerator()[i]<rgd*s.gamma().numerator()[i];

    assert(false); return false; // cannot happen since |r.gamma()!=s.gamma()|
  }

  // and when neither |x| nor |gamma()| discriminate, use the |y| component
  return r.y()<s.y(); // uses |SmallBitVector::operator<|, internal comparison
}

SR_poly Rep_context::expand_final(StandardRepr z) const // by value
{
  const RootDatum& rd = rootDatum();
  const InvolutionTable& i_tab = innerClass().involution_table();

  make_dominant(z); // this simplifies matters a lot; |z| is unchanged hereafter

  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const RatWeight& gamma=z.gamma();

  RankFlags singular_real_parity;
  for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
    if (rd.simpleCoroot(s).dot(gamma.numerator())==0)
    { if (i_tab.is_real_simple(i_x,s))
	singular_real_parity.set // record whether |s| is a real parity root
	  // |unpack| gives $(1-\theta)(\lambda-\rho)$
	  // real simple coroot odd on $\lambda-\rho$ means it is parity
	  (s,rd.simpleCoroot(s).dot(i_tab.y_lift(i_x,z.y()))%4!=0);
      else if (i_tab.is_imaginary_simple(i_x,s))
      {
	if (kgb().status(s,z.x())==gradings::Status::ImaginaryCompact)
	  return SR_poly(repr_less());; // return a zero result
      }
      else
	assert(not kgb().isComplexDescent(s,z.x())); // because |make_dominant|
    }
  // having made dominant, any non-final is witnessed on a (real) simple root
  if (singular_real_parity.any())
  {
    const weyl::Generator s= singular_real_parity.firstBit();
    const KGBEltPair p = kgb().inverseCayley(s,z.x());
    Weight lr = lambda_rho(z);

    // |lr| may need replacement by an equivalent (at |z.x()|) before it is
    // passed to inverse Cayleys, which will reinterpret is at $x$'s from |p|
    assert(rd.simpleCoroot(s).dot(lr)%2!=0); // we tested this odd above
    lr -= // correct so that $\<\lambda,\alpha^\vee>=0$ (like $\gamma$)
      rd.simpleRoot(s)*((rd.simpleCoroot(s).dot(lr)+1)/2); // project on perp
    assert(rd.simpleCoroot(s).dot(lr)==-1); // because $\<\rho,\alpha^\vee>=1$

    SR_poly result = expand_final(sr(p.first,lr,gamma));
    if (p.second!=UndefKGB)
      result += expand_final(sr(p.second,lr,gamma));
    return result;
  }
  else return SR_poly(z,repr_less());
} // |Rep_context::expand_final|

void Rep_table::add_block(param_block& block, BlockEltList& survivors)
{
  survivors.reserve(block.size());
  for (BlockElt x=0; x<block.size(); ++x)
    if (block.survives(x))
      survivors.push_back(x);

  unsigned long old_size = hash.size();
  BlockEltList new_survivors;

  // fill the |hash| table for new surviving parameters in this block
  for (BlockEltList::const_iterator
	 it=survivors.begin(); it!=survivors.end(); ++it)
    if (hash.match(sr(block,*it))>=old_size)
      new_survivors.push_back(*it);

  assert(hash.size()==old_size+new_survivors.size()); // only new surv. added
  if (new_survivors.empty())
    return; // nothing left to do, but we have computed |survivors| for caller

  lengths.resize(hash.size());
  KL_list.resize(hash.size(),SR_poly(repr_less())); // new slots, init empty
  def_formula.resize(hash.size(),SR_poly(repr_less())); // allocate new slots

  // compute cumulated KL polynomimals $P_{x,y}$ with $x\leq y$ survivors

  // start with computing KL polynomials for the entire block
  const kl::KLContext& klc = block.klc(block.size()-1,false); // silently

  /* get $P(x,z)$ for |x<=z| with |z| among new |survivors|, and contribute
   parameters from |block.survivors_below(x)| with coefficient $P(x,z)[q:=s]$
   to the |SR_poly| at |KL_list[old_size+i], where |z=new_survivors[i]| */

  BlockEltList::const_iterator z_start=new_survivors.begin();
  for (BlockEltList::const_iterator it = z_start; it!=new_survivors.end(); ++it)
    lengths[old_size+(it-z_start)]=block.length(*it);

  for (BlockElt x=0; x<=new_survivors.back(); ++x)
  {
    BlockEltList xs=block.survivors_below(x);
    if (xs.empty())
      continue; // no point doing work for |x|'s that don't contribute anywhere

    const unsigned int parity = block.length(x)%2;

    if (*z_start<x)
      ++z_start; // advance so |z| only runs over values with |x<=z|
    assert(z_start!=new_survivors.end() and *z_start>=x);

    for (BlockEltList::const_iterator it=z_start; it!=new_survivors.end(); ++it)
    {
      const BlockElt z = *it; // element of |new_survivors| and |x<=z|
      const kl::KLPol& pol = klc.klPol(x,z); // regular KL polynomial
      Split_integer eval(0);
      for (polynomials::Degree d=pol.size(); d-->0; )
	eval = eval.times_s()+Split_integer(static_cast<int>(pol[d]));
      if (eval!=Split_integer(0))
      {
	unsigned long z_index = old_size+(it-new_survivors.begin());
	assert(hash.find(sr(block,z))==z_index);
	SR_poly& dest = KL_list[z_index];
	if (lengths[z_index]%2!=parity)
	  eval.negate(); // incorporate sign for length difference
	for (unsigned int i=0; i<xs.size(); ++i)
	  dest.add_term(sr(block,xs[i]),eval);
      }
    } // |for(it)|
  } // |for(x)|
} // |Rep_table::add_block|

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

// compute and return sum of KL polynomials at $s$ for final parameter |z|
SR_poly Rep_table::KL_column_at_s(StandardRepr z) // must be nonzero and final
{
  { RootNbr witness;
    assert(not is_zero(z,witness));
    assert(is_final(z,witness));
    ndebug_use(witness);
  }
  make_dominant(z); // so that |z| it will appear at the top of its own block
  unsigned long hash_index=hash.find(z);
  if (hash_index==hash.empty) // previously unknown parameter
  {
    param_block block(*this,z);
    BlockEltList survivors;
    add_block(block,survivors);

    hash_index=hash.find(z);
    assert(hash_index!=hash.empty);
  }

  return KL_list[hash_index];
} // |Rep_table::KL_column_at_s|

SR_poly Rep_table::deformation_terms (param_block& block,BlockElt entry_elem)
{
  SR_poly result(repr_less());
  if (not block.survives(entry_elem) or block.length(entry_elem)==0)
    return result; // easy cases, null result

  BlockEltList survivors;
  add_block(block,survivors); // computes survivors, and add anything new

  assert(hash.find(sr(block,entry_elem))!=hash.empty); // should be known now

  // count number of survivors of length strictly less than any occurring length
  std::vector<unsigned int> n_surv_length_less
    (block.length(survivors.back())+1); // slots for lengths |<=| largest length
  { // compute |n_surv_length_less| values
    unsigned int l=0;  n_surv_length_less[l]=0;
    for (BlockEltList::const_iterator
	   it=survivors.begin(); it!=survivors.end(); ++it)
      while (l<block.length(*it))
	n_surv_length_less[++l] = it-survivors.begin();
  }

  // map indices of |survivors| to corresponding number in |hash|
  std::vector<unsigned long> remap(survivors.size());
  for (unsigned long i=0; i<survivors.size(); ++i)
  {
    unsigned long h=hash.find(sr(block,survivors[i]));
    assert(h!=hash.empty);
    remap[i]=h;
  }

  SR_poly rem(sr(block,entry_elem),repr_less()); // remainder = 1*entry_elem
  std::vector<Split_integer> acc(survivors.size(),Split_integer(0));

  for (unsigned long i=survivors.size(); i-->0; ) // decreasing essential here
  {
    StandardRepr p_y=sr(block,survivors[i]);
    Split_integer c_y = rem[p_y];
    const SR_poly& KL_y = KL_list[remap[i]];
    rem.add_multiple(KL_y,-c_y);
    assert(rem[p_y]==Split_integer(0)); // check relation of being inverse

    c_y.times_1_s(); // deformation terms are all multiplied by $1-s$
    acc[i]=c_y; // store coefficient at index of survivor
  }
  assert(rem.empty()); // since all terms in |KL_y| should be at most $y$

  // $\sum_{x\leq y<ee}y[l(ee)-l(y) odd] (-1)^{l(x)-l(y)}P_{x,y}*Q(y,ee)$
  unsigned int ll=block.length(entry_elem)-1; // last length of contributing |y|

  for (unsigned int l=ll%2; l<=ll; l+=2) // length of parity opposite |ee|
    for (BlockElt yy=n_surv_length_less[l]; yy<n_surv_length_less[l+1]; ++yy)
      result.add_multiple(KL_list[remap[yy]],acc[yy]);

  // correct signs in terms of result according to orientation numbers
  unsigned int orient_ee = orientation_number(sr(block,entry_elem));
  for (SR_poly::iterator it=result.begin(); it!=result.end(); ++it)
  {
    unsigned int orient_x=orientation_number(it->first);
    assert((orient_ee-orient_x)%2==0);
    int orient_express = (orient_ee-orient_x)/2;
    if (orient_express%2!=0)
      it->second.times_s();
  }

  return result;
} // |deformation_terms|

SR_poly Rep_table::deformation(const StandardRepr& z)
{
  Weight lam_rho = lambda_rho(z);
  RatWeight nu_z =  nu(z);
  StandardRepr z0 = sr(z.x(),lam_rho,RatWeight(rank()));
  SR_poly result = expand_final(z0); // value without deformation terms

  RationalList rp=reducibility_points(z);
  if (rp.size()==0) // without deformation terms
    return result; // don't even bother to store the result

  StandardRepr z_near = sr(z.x(),lam_rho,nu_z*rp.back());
  make_dominant(z_near);

  unsigned long h=hash.find(z_near);
  { // look up if closest reducibility point to |z| is already known
    unsigned long h=hash.find(z_near);
    if (h!=hash.empty and not def_formula[h].empty())
      return def_formula[h];
  }

  for (unsigned i=rp.size(); i-->0; )
  {
    Rational r=rp[i];
    const StandardRepr zi = sr(z.x(),lam_rho,nu_z*r);
    param_block b(*this,zi);
    const SR_poly terms = deformation_terms(b,b.size()-1);
    for (SR_poly::const_iterator it=terms.begin(); it!=terms.end(); ++it)
      result.add_multiple(deformation(it->first),it->second); // recursion
  }

  // now store result for future lookup
  h=hash.find(z_near);
  assert(h!=hash.empty); // it should have been added by |deformation_terms|
  def_formula[h]=result;

  return result;
} // |Rep_table::deformation|

SR_poly Rep_table::twisted_deformation_terms
  (param_block& block,BlockElt entry_elem)
{
  const auto& delta = innerClass().distinguished();
  const auto sr_y = sr(block,entry_elem);

  assert(is_twist_fixed(sr_y,delta));

  SR_poly result(repr_less());
  if (not block.survives(entry_elem) or block.length(entry_elem)==0)
    return result; // easy cases, null result

  ext_block::ext_block eblock(innerClass(),block,kgb(),delta);
  add_block(eblock,block);
  assert(eblock.is_present(entry_elem)); // since |is_twist_fixed| succeeded

  assert(hash.find(sr_y)!=hash.empty); // |sr_y| should be known now
  unsigned int parity = length(sr_y)%2;

  SR_poly rem(sr_y,repr_less()); // remainder = 1*entry_elem

  do
  { const auto& term = *rem.rbegin();
    const StandardRepr& p_x= term.first;
    const Split_integer& c_x = term.second;
    const SR_poly& KL_x = twisted_KLV_list[hash.find(p_x)];
    rem.add_multiple(KL_x,-c_x);
    assert(rem[p_x]==Split_integer(0)); // check relation of being inverse
    if (length(p_x)%2!=parity)
      result.add_multiple(KL_x,-c_x);
  }
  while (not rem.empty());

  // correct signs in terms of result according to orientation numbers
  unsigned int orient_y = orientation_number(sr_y);
  for (SR_poly::iterator it=result.begin(); it!=result.end(); ++it)
  {
    static constexpr Split_integer factor[2] // tabulate $(1-s)*s^i$ for $i=0,1$
      = { Split_integer {1,-1}, Split_integer {-1,1} };
    const unsigned int // conversion to |unsigned| needed for proper div/mod ops
      diff = orient_y-orientation_number(it->first);
    assert(diff%2==0);
    auto ev = it->second.e()-it->second.s(); // evaluation at $s=-1$ matters
    it->second = factor[(diff/2)%2]*ev;
  }

  return result;
} // |twisted_deformation_terms|

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
  for (auto it=survivors.begin(); not survivors.at_end(it); ++it)
  {
    auto x = *it; auto factor = P_at_s(x,0);
    if (eblock.length(x)%2!=parity) // flip sign at odd length difference
      factor = -factor;
    result.add_term(rc.sr(parent,eblock.z(x)),factor);
  }
  return result;
} // |twisted_KL_sum|

void Rep_table::add_block(ext_block::ext_block& block,
			  param_block& parent) // its complete unextended block
{
  const unsigned long old_size = hash.size();
  { BlockEltList survivors;  // this exported value will not be used
    add_block(parent,survivors); // but we must ensure parent block is known
  }
  BlockEltList new_survivors; new_survivors.reserve(block.size());

  // fill the |hash| table for new surviving parameters in this block
  RankFlags singular_orbits = block.singular_orbits(parent);
  for (BlockElt ez = 0; ez<block.size(); ++ez)
    if (hash.find(sr(parent,block.z(ez)))>=old_size and
	block.first_descent_among(singular_orbits,ez)==block.rank())
      new_survivors.push_back(ez);

  // extend space in twisted tables
  twisted_KLV_list.resize(hash.size(),SR_poly(repr_less())); // init empty
  twisted_def_formula.resize(hash.size(),SR_poly(repr_less()));

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
  containers::simple_list<BlockElt> survivors = block.condense(P_at_s,parent);

  // finally transcribe columns from |P_at_s| into |twisted_KLV_list|
  for (auto it=survivors.begin(); not survivors.at_end(it); ++it)
  { assert(block.first_descent_among(singular_orbits,*it)==block.rank());
    auto y = *it;
    auto y_index = hash.find(sr(parent,block.z(y)));
    assert (y_index!=hash.empty); // since we looked up everything above
    assert (y_index<twisted_KLV_list.size());
    SR_poly& dest = twisted_KLV_list[y_index];
    if (dest.empty()) // this means the entry was never defined
    { dest = SR_poly(sr(parent,block.z(y)),repr_less()); // coefficient 1
      unsigned int parity = block.length(y)%2;
      for (auto x_it=survivors.begin(); x_it!=it; ++x_it) // upper part
      {
	auto x = *x_it;
	auto factor = P_at_s(x,y);
	if (block.length(x)%2!=parity) // flip sign at odd length difference
	  factor = -factor;
	dest.add_term(sr(parent,block.z(x)),factor);
      }
      // since |dest| is a reference, the sum is stored at its destination
    } // |if (y_index>=old_size)|
  } // |for (it)|
} // |Rep_table::add_block| (extended block)

// compute and return sum of KL polynomials at $s$ for final parameter |z|
SR_poly twisted_KL_column_at_s
  (const Rep_context& rc, StandardRepr z, const WeightInvolution& delta)
  // |z| must be delta-fixed, nonzero and final
{
  { RootNbr witness; // there should be no imaginary or real singular descents
    if (rc.is_zero(z,witness) or not rc.is_final(z,witness))
      throw std::runtime_error("Representation zero or not final");
  }
  rc.make_dominant(z);
  BlockElt entry; // dummy needed to ensure full block is generated
  param_block block(rc,z,entry); // which this constructor does
  ext_block::ext_block eblock(rc.innerClass(),block,rc.kgb(),delta);

  return twisted_KL_sum(rc,eblock,eblock.element(entry),block);
}

// look up or compute and return the alternating sum of twisted KL polynomials
// at the inner class involution for final parameter |z|, evaluated at $q=s$
SR_poly Rep_table::twisted_KL_column_at_s(StandardRepr z)
  // |z| must be inner-class-twist-fixed, nonzero and final
{
  { RootNbr witness;
    if (is_zero(z,witness) or not is_final(z,witness))
      throw std::runtime_error("Representation zero or not final");
  }
  make_dominant(z);
  unsigned long hash_index=hash.find(z);
  if (hash_index==hash.empty // previously unknown parameter
      or hash_index>=twisted_KLV_list.size() // block known but not extended
      or twisted_KLV_list[hash_index].empty()) // same, but skipped over
  {
    BlockElt entry; // dummy needed to ensure full block is generated
    param_block block(*this,z,entry); // which this constructor does
    const auto &ic = innerClass();
    ext_block::ext_block eblock(ic,block,kgb(),ic.distinguished());

    add_block(eblock,block);

    if (hash_index==hash.empty) // then reset |hash_index| to now known value
    { hash_index=hash.find(z);
      assert(hash_index!=hash.empty);
    }
  }

  return twisted_KLV_list[hash_index];
}


  } // |namespace repr|
} // |namespace atlas|
