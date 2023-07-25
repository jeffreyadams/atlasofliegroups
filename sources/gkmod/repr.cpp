/*
  This is repr.cpp

  Copyright (C) 2009-2022 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "repr.h"

#include <memory> // for |std::unique_ptr|
#include <cstdlib> // for |std::abs|
#include <cassert>
#include <map> // used in computing |reducibility_points|
#include <iostream> // for progress reports and easier debugging
#include "error.h"

#include "arithmetic.h"
#include "matreduc.h"

#include "tits.h"
#include "kgb.h"	// various methods
#include "blocks.h"	// the |blocks::common_block| class, |dual_involution|
#include "subsystem.h" // |SubSystem| methods
#include "alcoves.h"

#include "kl.h"

#include "ext_block.h"
#include "ext_kl.h"

#include "basic_io.h"   // lookup of |operator<<| instances

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

StandardReprMod StandardReprMod::mod_reduce
  (const Rep_context& rc, const StandardRepr& sr)
{
  KGBElt x = sr.x();
  auto gam_lam=sr.gamma()-rho(rc.root_datum())-rc.lambda_rho(sr);
  rc.involution_table().real_unique(rc.kgb().inv_nr(x),gam_lam);
  return StandardReprMod(x,std::move(gam_lam));
}

StandardReprMod StandardReprMod::build
  (const Rep_context& rc, KGBElt x, RatWeight gam_lam)
{
  rc.involution_table().real_unique(rc.kgb().inv_nr(x),gam_lam);
  return StandardReprMod(x,std::move(gam_lam)); // ctor normalises
}

size_t StandardReprMod::hashCode(size_t modulus) const
{ size_t hash = x_part + 47*gamlam.denominator();
  for (auto entry : gamlam.numerator())
    hash= 11*hash+entry;
  return hash &(modulus-1);
}

Reduced_param::Reduced_param
  (const Rep_context& rc, const StandardReprMod& srm)
    : x(srm.x())
    , int_sys_nr()
    , w()
    , evs_reduced() // |int_sys_nr| is set by |integral_eval| below
{
  InnerClass& ic = rc.inner_class();
  const KGB& kgb = rc.kgb();
  const auto& gl = srm.gamma_lambda(); // $\gamma-\lambda$
  auto eval = ic.integral_eval(gl,int_sys_nr,w) * gl.numerator(); // mat * vec
  for (auto& entry : eval)
  {
    assert(entry%gl.denominator()==0);
    entry /= gl.denominator();
  }
  const auto codec =
    ic.int_item(int_sys_nr).data(ic,kgb.inv_nr(x),w);
  evs_reduced = codec.in * // transform coordinates to $1-\theta$-adapted basis
    int_Vector(eval.begin(),eval.end());
  for (unsigned int i=0; i<codec.diagonal.size(); ++i)
    evs_reduced[i] = arithmetic::remainder(evs_reduced[i],codec.diagonal[i]);
}

Reduced_param::Reduced_param
  (const Rep_context& rc, const StandardReprMod& srm,
   unsigned int_sys_nr, const WeylElt& w)
    : x(srm.x())
    , int_sys_nr(int_sys_nr)
    , w(w)
    , evs_reduced() // |int_sys_nr| is set by |integral_eval| below
{
  InnerClass& ic = rc.inner_class();
  const KGB& kgb = rc.kgb();
  const auto& gl = srm.gamma_lambda(); // $\gamma-\lambda$
  int_Matrix M = ic.integral_eval(int_sys_nr,w);
  auto eval = M * gl.numerator(); // mat * vec
  for (auto& entry : eval)
  {
    assert(entry%gl.denominator()==0);
    entry /= gl.denominator();
  }
  const subsystem::integral_datum_item::codec codec { ic, kgb.inv_nr(x), M };
  evs_reduced = codec.in * // transform coordinates to $1-\theta$-adapted basis
    int_Vector(eval.begin(),eval.end());
  for (unsigned int i=0; i<codec.diagonal.size(); ++i)
    evs_reduced[i] = arithmetic::remainder(evs_reduced[i],codec.diagonal[i]);
}

size_t Reduced_param::hashCode(size_t modulus) const
{ size_t hash = 7*x + 83*int_sys_nr;
  for (auto val : evs_reduced)
    hash= 25*hash+val;
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

  const InvolutionTable& i_tab = involution_table();
  auto i_x = kgb().inv_nr(x);
  const auto& theta = i_tab.matrix(i_x);
  auto th1_gamma_num = // numerator of $(1+\theta)*\gamma$ as 64-bits vector
    gamma.numerator()+theta*gamma.numerator();

  // since $(1+\theta)*\gamma = (1+\theta)*\lambda$ it actually lies in $X^*$
  Weight th1_gamma(th1_gamma_num.size());
  for (unsigned i=0; i<th1_gamma_num.size(); ++i)
  {
    assert(th1_gamma_num[i]%gamma.denominator()==0);
    th1_gamma[i] = th1_gamma_num[i]/gamma.denominator();
  }
#ifndef NDEBUG // check that constructor below builds a valid |StandardRepr|
  {
    Weight image = // $(\theta+1)(\gamma-\rho)$
      th1_gamma-i_tab.theta_plus_1_rho(i_x);
    matreduc::find_solution(theta+1,image); // assert that a solution exists
  }
#endif

  return StandardRepr(x, i_tab.y_pack(i_x,lambda_rho), gamma,
		      height(th1_gamma));
} // |sr_gamma|

// the same, but moving from |gamma|
StandardRepr Rep_context::sr_gamma
  (KGBElt x, const Weight& lambda_rho, RatWeight&& gamma) const
{
  const InvolutionTable& i_tab = involution_table();
  auto i_x = kgb().inv_nr(x);
  const auto& theta = i_tab.matrix(i_x);
  auto th1_gamma_num = // numerator of $(1+\theta)*\gamma$ as 64-bits vector
    gamma.numerator()+theta*gamma.numerator();

  // since $(1+\theta)*\gamma = (1+\theta)*\lambda$ it actually lies in $X^*$
  Weight th1_gamma(th1_gamma_num.size());
  for (unsigned i=0; i<th1_gamma_num.size(); ++i)
  {
    assert(th1_gamma_num[i]%gamma.denominator()==0);
    th1_gamma[i] = th1_gamma_num[i]/gamma.denominator();
  }

  return StandardRepr(x, i_tab.y_pack(i_x,lambda_rho),std::move(gamma),
		      height(th1_gamma));
} // |sr_gamma|

StandardRepr Rep_context::sr
  (const StandardReprMod& srm, const RatWeight& gamma) const
{
  const RatWeight gamma_lambda_rho = srm.gamma_lambda()+rho(root_datum());
  const auto lambda_rho = gamma.integer_diff<int>(gamma_lambda_rho);
  return sr_gamma(srm.x(),lambda_rho,gamma);
}

StandardRepr Rep_context::sr
  (const StandardReprMod& srm, const RatWeight& diff, const RatWeight& gamma)
  const
{
  const RatWeight gamma_lambda_rho = srm.gamma_lambda()+rho(root_datum())+diff;
  const auto lambda_rho = gamma.integer_diff<int>(gamma_lambda_rho);
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
  const InvolutionTable& i_tab = involution_table();
  const RatWeight lambda = rho(root_datum())+lambda_rho;
  const RatWeight diff = lambda - nu;
  const RatWeight theta_diff(i_tab.matrix(kgb().inv_nr(x))*diff.numerator(),
			     diff.denominator()); // theta(lambda-nu)
  return ((lambda+nu+theta_diff)/=2).normalize();
}

const WeightInvolution& Rep_context::theta (const StandardRepr& z) const
{ return involution_table().matrix(kgb().inv_nr(z.x())); }

Weight Rep_context::lambda_rho(const StandardRepr& z) const
{
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);


  // recover $\lambda-\rho$ from doubled projections on eigenspaces $\theta$
  const RatWeight gam_rho = z.gamma() - rho(root_datum());
  auto th1_gam_rho_num = // numerator of $(1+\theta)*(\gamma-\rho)$, 64 bits
    gam_rho.numerator() + theta*gam_rho.numerator();

  Weight th1_gam_rho(th1_gam_rho_num.size());
  for (unsigned i=0; i<th1_gam_rho_num.size(); ++i)
  {
    assert(th1_gam_rho_num[i]%gam_rho.denominator()==0);
    th1_gam_rho[i] = th1_gam_rho_num[i]/gam_rho.denominator();
  }

  // the next addition is |Weight::operator+(Weight&&) const &|:
  return ( th1_gam_rho + i_tab.y_lift(i_x,z.y()) ) / 2; // exact division
}

// compute $\gamma-\lambda$ from $(1-\theta)(\gamma-\lambda)=2(\gamma-\lambda)$
RatWeight Rep_context::gamma_lambda
  (InvolutionNbr i_x,  const TorusPart& y_bits, const RatWeight& gamma) const
{
  const InvolutionTable& i_tab = involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);

  // |y_lift(i_x,y_bits==(1-theta)*(lambda-rho)|; get |(1-theta)(gamma-lambda)|
  const RatWeight gamma_rho = gamma - rho(root_datum());
  return (gamma_rho-theta*gamma_rho - i_tab.y_lift(i_x,y_bits))
    /static_cast<arithmetic::Numer_t>(2);
}

// compute $\gamma-\lambda-\rho$ from same information; with respect to above,
// change from subtracting |(1-theta)*rho| to adding |(1+theta)*rho|
RatWeight Rep_context::gamma_lambda_rho (const StandardRepr& sr) const
{
  const InvolutionTable& i_tab = involution_table();
  InvolutionNbr i_x = kgb().inv_nr(sr.x());
  const WeightInvolution& theta = i_tab.matrix(i_x);
  const RatWeight& gamma = sr.gamma();

  return (  gamma - theta*gamma
	 + (i_tab.theta_plus_1_rho(i_x) - i_tab.y_lift(i_x,sr.y()))
	 ) /static_cast<arithmetic::Numer_t>(2);
}

RatWeight Rep_context::gamma_0 (const StandardRepr& z) const
{
  const InvolutionTable& i_tab = involution_table();
  const auto& theta = i_tab.matrix(kgb().inv_nr(z.x()));
  return ((z.gamma()+theta*z.gamma())/=2).normalize();
}

RatWeight Rep_context::nu(const StandardRepr& z) const
{
  const InvolutionTable& i_tab = involution_table();
  const auto& theta = i_tab.matrix(kgb().inv_nr(z.x()));
  return ((z.gamma()-theta*z.gamma())/=2).normalize();
}


bool Rep_context::is_parity_at_0(RootNbr i,const StandardRepr& z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = involution_table();
  const RootNbrSet& real_roots = i_tab.real_roots(i_x);
  assert(real_roots.isMember(i));
  const RootNbrSet non_real_roots = ~real_roots; // the complement

  const Coweight& alpha_hat = rd.coroot(i);

  Weight theta_1_lamrho = i_tab.y_lift(i_x,z.y()); // |(1-theta)*lam_rho|
  int eval = // twice evaluation of |alpha_hat| on |\lambda-\rho_{rea]}|
    alpha_hat.dot(theta_1_lamrho+rd.twoRho(non_real_roots));
  assert (eval%2==0); // because |\lambda-\rho_{rea]}| is integral
  return eval%4!=0;
}

bool Rep_context::is_parity(RootNbr i,const StandardRepr& z) const
{
  assert(involution_table().real_roots(kgb().inv_nr(z.x())).isMember(i));
  RatNum eval = z.gamma().dot_Q(root_datum().coroot(i)).normalize();
  assert(eval.denominator()==1); // must be an integral coroot
  return is_parity_at_0(i,z) == (eval.numerator()%2==0);
}

StandardReprMod Rep_context::inner_twisted(const StandardReprMod& z) const
{
  const auto& delta = inner_class().distinguished();
  return StandardReprMod::build(*this,kgb().twisted(z.x(),delta),
				delta*gamma_lambda(z));
}

/* the |evs_reduced| field of a |Reduced_param| encoode the evaluations of the
   value $\gamma-\lambda$ of a parameter on (simply-) integral coroots for
   $\gamma$. However, since $\lambda$ is defined only up to sifts in the image
   lattice of $1-\theta$, the |ev_reduced| values are effectively reduced modulo
   the image of that lattice under the integral coroot evaluation map. It may
   therefore happen that when subtracting the values of $\gamma-\lambda$ for two
   |StandardReprMod| values producing the same |Reduced_param|, these coroot
   evaluations of the difference are not zero; they must however define a point
   in the latter lattice image. The function |theta_1_preimage| makes this
   explicit: given such a difference |offset| it produces a weight in the image
   of $1-\theta$ that has identical intergral coroot evaluations as |offset|.
 */

// fixed choice in |(1-theta)X^*| of element given integral coroots evaluations
// uses idempotence of: (th_1_image*out*.))o(mod(diagonal))o(in*coroots_mat*.)
Weight Rep_context::theta_1_preimage
  (const RatWeight& offset, const subsystem::integral_datum_item::codec& codec)
  const
{
  RatWeight eval = (codec.coroots_matrix*offset).normalize();
  if (eval.is_zero())
    return Weight(codec.theta_1_image_basis.n_rows(),0);
  assert(eval.denominator()==1);
  auto eval_v = int_Vector(eval.numerator().begin(),eval.numerator().end());
  eval_v = codec.in * eval_v;

  { unsigned int i;
    for (i=0; i<codec.diagonal.size(); ++i)
    {
      assert(eval_v[i]%codec.diagonal[i]==0);
      eval_v[i]/=codec.diagonal[i];
    }
    for (; i<eval_v.size(); ++i)
      assert(eval_v[i]==0);
    eval_v.resize(codec.diagonal.size());
  }

  return codec.theta_1_image_basis * (codec.out * eval_v);
} // |theta_1_preimage|

// difference in $\gamma-\lambda$ from |srm0| with respect to that of |srm1|,
// for representatives of |srm0| and |srm1| with identical integral evaluations
RatWeight Rep_context::offset
  (const StandardReprMod& srm0, const StandardReprMod& srm1) const
{
  const auto& gamlam = srm0.gamma_lambda(); // will also define integral system
  RatWeight result = gamlam - srm1.gamma_lambda();
  if (not result.is_zero()) // optimize out a fairly frequent case
  {
    auto& ic = inner_class();
    const auto codec = ic.integrality_codec(gamlam,kgb().inv_nr(srm0.x()));
    result -= theta_1_preimage(result,codec); // ensure orthogonal to integral sys
    assert((codec.coroots_matrix*result).is_zero()); // check that it was done
  }
  return result;
}

StandardReprMod& Rep_context::shift
  (const RatWeight& shift, StandardReprMod& srm) const
{
  srm.gamlam += shift;
  involution_table().real_unique(kgb().inv_nr(srm.x()),srm.gamlam);
  return srm;
}

// |z| standard means (weakly) dominant on the (simply-)imaginary roots
bool Rep_context::is_standard(const StandardRepr& z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = involution_table();
  const Ratvec_Numer_t& numer = z.gamma().numerator();

  for (unsigned i=0; i<i_tab.imaginary_rank(i_x); ++i)
  {
    const RootNbr alpha = i_tab.imaginary_basis(i_x,i);
    if (rd.coroot(alpha).dot(numer)<0)
      return false;
  }
  return true;
}

// |z| dominant means precisely |gamma| is (weakly) dominant
bool Rep_context::is_dominant(const StandardRepr& z) const
{
  const RootDatum& rd = root_datum();
  const auto& numer = z.gamma().numerator();

  for (auto it=rd.beginSimpleCoroot(); it!=rd.endSimpleCoroot(); ++it)
    if (it->dot(numer)<0)
      return false;
  return true;
}

// |z| zero means that no singular simple-imaginary roots are compact; this
// code assumes |is_standard(z)|, namely |gamma| is dominant on imaginary roots
bool Rep_context::is_nonzero(const StandardRepr& z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = involution_table();
  const Ratvec_Numer_t& numer = z.gamma().numerator();

  for (unsigned i=0; i<i_tab.imaginary_rank(i_x); ++i)
  {
    const RootNbr alpha = i_tab.imaginary_basis(i_x,i);
    if (rd.coroot(alpha).dot(numer)==0 and // simple-imaginary, singular
	not kgb().simple_imaginary_grading(z.x(),alpha)) // and compact
      return false;
  }
  return true;
}

bool Rep_context::is_normal(const StandardRepr& z) const
{
  auto z_normal = z;
  normalise(z_normal);
  return z_normal==z;
}

// |z| semifinal means that no singular real roots satisfy the parity condition
// no assumptions about the real subsystem, so all real roots must be tested
bool Rep_context::is_semifinal(const StandardRepr& z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = involution_table();
  const RootNbrSet pos_real = i_tab.real_roots(i_x) & rd.posroot_set();
  const Weight test_wt = i_tab.y_lift(i_x,z.y()) // $(1-\theta)(\lambda-\rho)$
	   + rd.twoRho()-rd.twoRho(pos_real); // replace $\rho$ by $\rho_R$

  for (RootNbrSet::iterator it=pos_real.begin(); it(); ++it)
  {
    const Weight& av = root_datum().coroot(*it);
    if (av.dot(z.gamma().numerator())==0 and
	av.dot(test_wt)%4 !=0) // singular yet odd on shifted lambda
      return false;
  }
  return true;
}

bool Rep_context::is_final(const StandardRepr& z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionTable& i_tab = involution_table();
  const auto& numer = z.gamma().numerator();
  KGBElt x=z.x();
  const InvolutionNbr i_x = kgb().inv_nr(x);

  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
  {
    auto v = rd.simpleCoroot(s).dot(numer);
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
      } // tests on |v|
  } // |for(s)|
  return is_nonzero(z); // check \emph{all} simply-imaginary coroots
}


bool Rep_context::is_oriented(const StandardRepr& z, RootNbr alpha) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = involution_table();
  const RootNbrSet real = i_tab.real_roots(i_x);

  assert(real.isMember(alpha)); // only real roots should be tested

  const Coweight& av = root_datum().coroot(alpha);
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
  const InvolutionTable& i_tab = involution_table();
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
	auto eps = av.dot(test_wt)%4==0 ? 0 : denom;
	if ((num>0) == // either positive for gamma and oriented, or neither
	    (arithmetic::remainder(num+eps,2*denom) < denom))
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

RankFlags Rep_context::singular_simples (const StandardRepr& z) const
{
  const auto& rd=root_datum();
  assert(is_dominant_ratweight(rd,z.infinitesimal_char));
  const auto& numer = z.infinitesimal_char.numerator();
  RankFlags result;
  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
    result.set(s,rd.simpleCoroot(s).dot(numer)==0);
  return result;
}

WeylWord Rep_context::complex_descent_word (KGBElt x, RankFlags singulars) const
{
  WeylWord result;
  { RankFlags::iterator it;
    do
      for (it=singulars.begin(); it(); ++it)
	if (kgb().isComplexDescent(*it,x))
	{
	  auto s=*it;
	  x = kgb().cross(s,x);
	  result.push_back(s);
	  break; // out of the loop |for(it)|
	} // |if(isComplexDescent)|
    while (it()); // wait until inner loop runs to completion
  }
  return result;
}

void Rep_context::make_dominant(StandardRepr& z) const
{
  const RootDatum& rd = root_datum();

  // the following are non-|const|, and modified in the loop below
  Weight lr = lambda_rho(z);
  KGBElt& x = z.x_part; // the |z.x_part| will be modified in-place
  Ratvec_Numer_t& numer = // and |z.infinitesimal_char| as well
    z.infinitesimal_char.numerator();

  { weyl::Generator s;
    do
      for (s=0; s<rd.semisimple_rank(); ++s)
	if (rd.simpleCoroot(s).dot(numer)<0)
	{
	  int offset; // used to pivot |lr| around $\rho_r-\rho$
	  switch (kgb().status(s,x))
	  {
	  case gradings::Status::Complex: offset = 1; break;
	  case gradings::Status::Real:    offset = 0; break;
	  default: // |s| is an imaginary root; we will not cope with that here
	    throw std::runtime_error("Non standard parameter in make_dominant");
	  }
	  rd.simple_reflect(s,numer); // real or complex reflection of |gamma|
	  rd.simple_reflect(s,lr,offset);
	  x = kgb().cross(s,x);
	  break; // out of the loop |for(s)|
	} // |if(v<0)| and |for(s)|
    while (s<rd.semisimple_rank()); // wait until inner loop runs to completion
  }
  z.y_bits = involution_table().y_pack(kgb().inv_nr(x),lr);
} // |make_dominant|

// apply sequence of cross actions by singular complex simple generators
void Rep_context::complex_crosses (StandardRepr& z, const WeylWord& ww) const
{
  const auto& rd = root_datum();
#ifndef NDEBUG
  const InvolutionTable& i_tab = involution_table();
#endif
  auto& x = z.x_part; // directly operate on |x| component inside |z|
  Weight lr = lambda_rho(z);
  // |z.infinitesimal_char| is unchanged by singular reflections

  for (auto s : ww)
  {
    assert(rd.simpleCoroot(s).dot(z.gamma().numerator())==0);
    assert(i_tab.is_complex_simple(kgb().inv_nr(x),s));
    x = kgb().cross(s,x);
    rd.simple_reflect(s, lr, 1);
  }

  // reinsert $y$ bits component
  z.y_bits = involution_table().y_pack(kgb().inv_nr(x),lr);
}

// auxiliary: move to canonical involution for (singular) |gens| subgroup of $W$
void Rep_context::to_singular_canonical(RankFlags gens, StandardRepr& z) const
{ // simply-singular coroots are simple, so no need to construct a subsystem
  TwistedInvolution tw = kgb().involution(z.x_part); // copy to be modified
  complex_crosses(z,inner_class().canonicalize(tw,gens));
  assert(tw == kgb().involution(z.x_part));
}

// make dominant and descend though any singular complex descents
// this is a version of |make_dominant| also doing singular |complex_crosses|
void Rep_context::deform_readjust(StandardRepr& z) const
{
  const RootDatum& rd = root_datum();

  // the following are non-|const|, and modified in the loop below
  Weight lr = lambda_rho(z);
  KGBElt& x = z.x_part; // the |z.x_part| will be modified in-place
  Ratvec_Numer_t& numer = // and |z.infinitesimal_char| as well
    z.infinitesimal_char.numerator();

  { weyl::Generator s;
    do
      for (s=0; s<rd.semisimple_rank(); ++s)
	if (kgb().status(s,x) == gradings::Status::Complex)
	{
	  auto eval = rd.simpleCoroot(s).dot(numer); // needed for sign only
	  if (eval<0)
	  {
	    rd.simple_reflect(s,numer); // real or complex reflection of |gamma|
	    rd.simple_reflect(s,lr,1);
	    x = kgb().cross(s,x);
	    break; // out of the loop |for(s)|
	  }
	else if (eval==0 and kgb().isDescent(s,x))
	{ // here |numer| will be unchanged
	  rd.simple_reflect(s,lr,1);
	  x = kgb().cross(s,x);
	  break; // out of the loop |for(s)|: other complex descents possible
	}
      } // |for(s)|
    while (s<rd.semisimple_rank()); // wait until inner loop runs to completion
  }
  z.y_bits = involution_table().y_pack(kgb().inv_nr(x),lr);
} // |deform_readjust|

// this also ensures a chosen singular-complex minimum when there are multiple
// but that only arises when singular-real descents exist (not so in deformation)
void Rep_context::normalise(StandardRepr& z) const
{
  make_dominant(z);

  RankFlags singulars = singular_simples(z);
  to_singular_canonical(singulars,z);

  complex_crosses(z,complex_descent_word(z.x(),singulars));
} // |normalise|

bool Rep_context::is_twist_fixed
  (StandardRepr z, const WeightInvolution& delta) const
{
  make_dominant(z);
  to_singular_canonical(singular_simples(z),z);

  return z==twisted(z,delta);
} // |is_twist_fixed|

// equivalence is equality after |make_dominant| and |to_singular_canonical|
bool Rep_context::equivalent(StandardRepr z0, StandardRepr z1) const
{
  if (kgb().Cartan_class(z0.x_part)!=kgb().Cartan_class(z1.x_part))
    return false; // this non-equivalence can be seen before |make_dominant|

  { // preempt failing of |make_dominant| on non standard parameters
    if (not (is_standard(z0) and is_standard(z1)))
      return z0==z1; // strict equality unless both are parameters standard
  }

  make_dominant(z0);
  make_dominant(z1);

  if (z0.infinitesimal_char!=z1.infinitesimal_char)
    return false;

  RankFlags singulars = singular_simples(z0);
  to_singular_canonical(singulars,z0);
  to_singular_canonical(singulars,z1);

  return z0==z1;
} // |Rep_context::equivalent|

StandardRepr Rep_context::scale(StandardRepr z, const RatNum& f) const
{ // we can just replace the |infinitesimal_char|, nothing else changes
  auto image = theta(z)*z.gamma();
  auto diff = z.gamma()-image; // this equals $2\nu(z)$
  z.infinitesimal_char += image;
  z.infinitesimal_char += diff*f; // now we have |(gamma_0(z)+nu(z)*f)*2|
  (z.infinitesimal_char/=2).normalize();
  return z;
}

RatNumList Rep_context::reducibility_points(const StandardRepr& z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = involution_table();
  const Permutation& theta = i_tab.root_involution(i_x);

  const RatWeight& gamma = z.gamma();
  const Ratvec_Numer_t& numer = gamma.numerator();
  const arithmetic::Numer_t d = gamma.denominator();
  const Weight lam_rho = lambda_rho(z);

  const RootNbrSet pos_real = i_tab.real_roots(i_x) & rd.posroot_set();
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

  RootNbrSet pos_complex = i_tab.complex_roots(i_x) & rd.posroot_set();
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

  std::set<RatNum> fracs;

  for (table::iterator it= evens.begin(); it!=evens.end(); ++it)
    for (long s= d*(it->second+2); s<=it->first; s+=2*d)
      fracs.insert(RatNum(s,it->first));

  for (table::iterator it= odds.begin(); it!=odds.end(); ++it)
    for (long s= it->second==0 ? d : d*(it->second+2); s<=it->first; s+=2*d)
      fracs.insert(RatNum(s,it->first));

  return RatNumList(fracs.begin(),fracs.end());
} // |reducibility_points|


StandardRepr Rep_context::cross(weyl::Generator s, StandardRepr z) const
{
  make_dominant(z);
  const RootDatum& rd = root_datum();
  const auto& i_tab = involution_table();
  const KGB& kgb = this->kgb();
  const RatWeight& gamma = z.gamma(); // now get the infinitesimal character
  const SubSystem& subsys = SubSystem::integral(rd,gamma);

  const auto refl = subsys.reflection(s);
  const KGBElt new_x = kgb.cross(refl,z.x());
  RootNbrSet pos_neg = pos_to_neg(rd,refl);
  pos_neg &= i_tab.real_roots(kgb.inv_nr(z.x())); // only real roots for |z|
  RatWeight gamma_lambda = this->gamma_lambda(z);
  gamma_lambda -= root_sum(rd,pos_neg); // correction for $\rho_r$'s
  rd.reflect(subsys.parent_nr_simple(s),gamma_lambda.numerator());

  const Weight lambda_rho = gamma.integer_diff<int>(gamma_lambda+rho(rd));
  return sr_gamma(new_x,lambda_rho,gamma);
}

StandardRepr Rep_context::cross(const Weight& alpha, StandardRepr z) const
{
  // this method does not apply or require any form of making |z| dominant first
  const RootDatum& rd = root_datum();
  KGBElt x= z.x_part;
  InvolutionNbr i_x = kgb().inv_nr(x);
  const InvolutionTable& i_tab = involution_table();

  const RatWeight& gamma=z.infinitesimal_char; // integrally dominant
  const RootNbr rt = rd.root_index(alpha);
  if (rt==rd.numRoots())
    throw std::runtime_error("Not a root");
  // the following test ensures that the |integer_diff| below won't fail
  if (rd.coroot(rt).dot(gamma.numerator())%gamma.denominator()!=0)
    throw std::runtime_error("Not an integral root");

  RatWeight gam_lam_shifted = gamma_lambda(z) +
    RatWeight(rd.twoRho(i_tab.real_roots(i_x)),2); // shift by $\rho_\R$
  Ratvec_Numer_t& lambda_numer = gam_lam_shifted.numerator();

  // transform |x|, |i_x|, and |gam_lam_shifted|, reflecting them by |alpha|
  i_x = kgb().inv_nr( x = kgb().cross(rd.reflection_word(rt),x) );
  rd.reflect(rt,lambda_numer);

  // shift back by $\rho_\R$ at (now) destination |i_x|
  gam_lam_shifted -= RatWeight(rd.twoRho(i_tab.real_roots(i_x)),2);

  const Weight lambda_rho = gamma.integer_diff<int>(gam_lam_shifted+rho(rd));
  return sr_gamma(x,lambda_rho,gamma);
}

StandardRepr Rep_context::Cayley(weyl::Generator s, StandardRepr z) const
{
  make_dominant(z);
  const RootDatum& rd = root_datum();
  const auto& i_tab = involution_table();
  const KGB& kgb = this->kgb();
  const RatWeight& gamma = z.gamma(); // now get the infinitesimal character
  const SubSystem& subsys = SubSystem::integral(rd,gamma);
  const auto parent_s = subsys.parent_nr_simple(s);

  const auto& conj = subsys.to_simple(s); // word in full system
  const KGBElt conj_x = kgb.cross(conj,z.x());
  RootNbrSet pos_neg = pos_to_neg(rd,conj);
  RatWeight gamma_lambda = this->gamma_lambda(z);
  const Coweight& alpha_hat = rd.coroot(parent_s);
  KGBElt new_x;

  if (kgb.status(subsys.simple(s),conj_x)==
      gradings::Status::ImaginaryNoncompact)
  {
    new_x = kgb.cross(kgb.cayley(subsys.simple(s),conj_x),conj);

    const RootNbrSet& upstairs_real_roots = i_tab.real_roots(kgb.inv_nr(new_x));
    RootNbrSet real_flip = upstairs_real_roots;
    real_flip ^= i_tab.real_roots(kgb.inv_nr(z.x())); // remove downstairs reals

    pos_neg &= real_flip; // posroots that change real status and map to negative
    gamma_lambda += root_sum(rd,pos_neg); // correction of $\rho_r$'s

    // correct in case the parity condition fails for our raised |gamma_lambda|
    const int rho_r_corr = // integer since alpha is among |upstairs_real_roots|
      alpha_hat.dot(rd.twoRho(upstairs_real_roots))/2;
    const int eval = gamma_lambda.dot(alpha_hat);
    if ((eval+rho_r_corr)%2==0) // parity condition says it should be 1
      gamma_lambda += RatWeight(rd.root(parent_s),2); // add half-alpha
  }
  else
  {
    const auto& real_roots = i_tab.real_roots(kgb.inv_nr(z.x()));
    if (not real_roots.isMember(parent_s))
      throw error::Cayley_error();
    const int eval = gamma_lambda.dot(alpha_hat);
    const int rho_r_corr = alpha_hat.dot(rd.twoRho(real_roots))/2;
    if ((eval+rho_r_corr)%2==0) // then |s| is real nonparity at |z|
      throw error::Cayley_error();

    new_x = kgb.cross(kgb.inverseCayley(subsys.simple(s),conj_x).first,conj);

    const RootNbrSet real_flip = real_roots^i_tab.real_roots(kgb.inv_nr(new_x));
    pos_neg &= real_flip; // posroots that change real status and map to negative
    gamma_lambda += root_sum(rd,pos_neg); // correction of $\rho_r$'s
    // now |gamma_lambda| is still in the $X^*$-coset of $\gamma-\rho$; it might
    // not be in the $-1$ eigenspace for |new_x|, but |sr_gamma| projects to it
  }

  const Weight lambda_rho = gamma.integer_diff<int>(gamma_lambda+rho(rd));
  return sr_gamma(new_x,lambda_rho,gamma);
}

/*
  Compute shift in |lambda| component of parameter for Cayley transform by a
  non-simple root $\alpha$, from involutions |theta_down| to |theta_up|, where
  |to_simple| left-conjugates root $\alpha$ to some simple root.

  Curiously, this appears to depend only on $\theta$ \emph{upstairs} and the
  conjugating element |to_simple|; an explanation is needed here. It seems to
  be because \emph{all} upstairs real roots becoming negative by the necessary
  conjugation will be downstairs complex roots (so contribute to the shift).

  Sum of positive real roots becoming negative at $\theta'=^{to\_simple}\theta$
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
  KGBElt& x= z.x_part; // thiss component will be modified in place
  InvolutionNbr i_x = kgb().inv_nr(x);
  const InvolutionTable& i_tab = involution_table();

  // the following are non-|const|, and modified in the loop below
  Weight lambda2_shifted = (lambda_rho(z)*=2)
    + rd.twoRho() - rd.twoRho(i_tab.real_roots(i_x));
  Ratvec_Numer_t& gamma_num = z.infinitesimal_char.numerator();

  sl_list<weyl::Generator> result;
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
	  i_x = kgb().inv_nr( x = kgb().cross(rd.reflection_word(alpha),x) );
	  rd.reflect(alpha,lambda2_shifted);
	  break; // out of the loop |for(s)|
	} // |if(v<0)|
      } // |for(s)|
    }
    while (s<subsys.rank()); // wait until inner loop runs to completion
  }
  lambda2_shifted -= rd.twoRho() - rd.twoRho(i_tab.real_roots(i_x)); // unshift
  z.y_bits=i_tab.y_pack(i_x,lambda2_shifted/2); // insert modified bits into |z|
  return { result.to_vector() }; // convert to |WeylWord|
} // |make_dominant| (integrally)

StandardRepr Rep_context::any_Cayley(const Weight& alpha, StandardRepr z) const
{
  const RootDatum& rd = root_datum();
  const KGB& kgb = this->kgb();
  const InvolutionTable& i_tab = involution_table();
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

  // now do the Cayley transform proper, with most work for handling |x|
  bool ascent; // whether forward Cayley, so that we can locate "upstairs"
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
	x = kgb.inverseCayley(s,x).first; // do inverse Cayley from |inv0|
	ascent=false;
	break;
      }
      // else FALL THROUGH
    }
  default: // |ImaginaryCompact| or |Complex|
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
  return twisted(z,inner_class().distinguished());
}

StandardRepr Rep_context::twisted
  (StandardRepr z, const WeightInvolution& delta) const
{
  const auto& i_tab = involution_table();
  const InvolutionNbr i_x0 = kgb().inv_nr(z.x());
  z.x_part = kgb().twisted(z.x_part,delta);
  const InvolutionNbr i_x1 =  kgb().inv_nr(z.x()); // destination involution
  z.y_bits = i_tab.y_act(i_x0,i_x1,z.y_bits,delta);
  z.infinitesimal_char = delta*z.infinitesimal_char;
  return z;
}


bool StandardRepr::operator<(const StandardRepr& s) const
{
  if (height()!=s.height()) // order by increasing height first
    return height()<s.height();
  if (x()!=s.x()) // then order by decreasing numeric value of |x|
    return x()>s.x(); // (height tends to change in opposite sense to |x|)
  if (y()!=s.y()) // then order by increasing internal value of |y|
    return y()<s.y(); // uses |SmallBitVector::operator<|, internal comparison

  // finally in rare cases individual components of |gamma| need comparison
  auto r_vec = s.gamma().numerator()*gamma().denominator(); // cross multiply
  auto s_vec = gamma().numerator()*s.gamma().denominator(); // cross multiply

  return r_vec<s_vec;
}

SR_poly Rep_context::scale(const poly& P, const RatNum& f) const
{
  poly result;
  for (const auto& term : P)
  {
    auto finals_term = finals_for(scale(term.first,f));
    for (auto it = finals_term.begin(); not finals_term.at_end(it); ++it)
      result.add_term(std::move(it->first),term.second*it->second);
  }
  return result;
}

K_repr::K_type_pol Rep_context::scale_0(const poly& P) const
{
  K_repr::K_type_pol result;
  for (auto it=P.begin(); it!=P.end(); ++it)
  { auto z=it->first; // take a copy for modification
    auto finals = finals_for(scale_0(z)); // a |simple_list| of K-type terms
    for (auto it=finals.begin(); not finals.at_end(it); ++it)
      result.add_term(std::move(it->first),Split_integer(it->second));
  }
  return result;
}

using sr_term = std::pair<StandardRepr,int>;
using sr_term_list = simple_list<sr_term>;

// insert (add) a new term into list |L|, assumed sorted decreasingly
void insert_into(sr_term_list& L, StandardRepr&& z, int coef)
{ auto it = L.begin();
  while (not L.at_end(it))
    if (z < it->first)
      ++it; // skip higher terms
    else if (it->first < z) // then we are looking at lower terms
      break; // so break loop and insert before those lower terms
    else  // matching term; operate on coefficient
    { if ((it->second += coef) == 0)
        L.erase(it);
      return; // whether by coefficient update or erasure, we are done
    }
  L.insert(it,std::make_pair(std::move(z),coef));
}

sr_term_list Rep_context::finals_for(StandardRepr z) const
{
  const RootDatum& rd = root_datum();

  sr_term_list result, to_do;
  to_do.emplace_front(std::move(z),1);

  do
  {
    KGBElt x = to_do.front().first.x();
    Weight lr = lambda_rho(to_do.front().first);
#ifndef NDEBUG
    auto height = to_do.front().first.height();
#endif
    RatWeight gamma = std::move(to_do.front().first.infinitesimal_char);
    auto coef = to_do.front().second;
    to_do.pop_front();

  restart:
    for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
      // as |break| from loop is not available within |switch|, use |goto| below
    { auto eval = // morally evaluation coroot at |gamma|, but only sign matters
	rd.simpleCoroot(s).dot(gamma.numerator());
      if (eval>0)
	continue; // when strictly dominant, nothing to do for |s|
      switch (kgb().status(s,x))
      {
      case gradings::Status::ImaginaryCompact:
	if (eval==0)
	  goto drop; // singular imaginary compact |s|, parameter is zero
	rd.simple_reflect(s,lr,1); // $-\rho$-based reflection
	rd.simple_reflect(s,gamma.numerator());
	coef = -coef;
	goto restart;
      case gradings::Status::ImaginaryNoncompact:
	if (eval==0)
	  continue; // nothing to do for singular nci generator |s|
	{ // |eval<0|: reflect, and also add Cayley transform terms
	  KGBElt sx = kgb().cross(s,x);
	  KGBElt Cx = kgb().cayley(s,x);
	  StandardRepr t1 = sr_gamma(Cx,lr,gamma);
	  assert( t1.height() < height );
	  insert_into(to_do,std::move(t1),coef);
	  if (sx==x) // then type 2 Cayley
	  {
	    StandardRepr t2 = sr_gamma(Cx,lr+rd.simpleRoot(s),gamma);
	    assert( t2.height() < height );
	    insert_into(to_do,std::move(t2),coef);
	  }
	  x = sx; // after testing we can update |x| for nci cross action
	  rd.simple_reflect(s,lr,1); // $-\rho$-based reflection
	  rd.simple_reflect(s,gamma.numerator());
	  coef = -coef; // reflect, negate, and continue with modified values
	  goto restart;
	}
      case gradings::Status::Complex:
	if (eval==0 and not kgb().isDescent(s,x))
	  continue; // nothing to do for singular complex ascent
	// now we are either not dominant for |s|, or a complex ascent
	x = kgb().cross(s,x);
	rd.simple_reflect(s,lr,1); // $-\rho$-based reflection
	rd.simple_reflect(s,gamma.numerator());
	// keep |coef| unchanged here
	goto restart;
      case gradings::Status::Real:
	if (eval==0) // singular real root
	{ auto eval_lr = rd.simpleCoroot(s).dot(lr);
	  if (eval_lr%2 == 0) // whether non-parity
	    continue; // nothing to do for a (singular) real nonparity root
	  // now $\alpha_s$ is parity real root: replace by inverse Cayley(s)
	  // |kgb()| can distinguish type 1 and type 2
	  lr -= rd.simpleRoot(s)*((eval_lr+1)/2); // project to wall for |s|
	  assert( rd.simpleCoroot(s).dot(lr) == -1 );
	  const KGBEltPair Cxs = kgb().inverseCayley(s,x);
	  if (Cxs.second!=UndefKGB)
	    insert_into(to_do,sr_gamma(Cxs.second,lr,gamma),coef);
	  insert_into(to_do,sr_gamma(Cxs.first,lr,std::move(gamma)),coef);
	  goto drop; // we have rewritten |current|, don't contribute it
	} // (singular real root)
	// |x = kgb().cross(s,x)|; real roots act trivially on KGB elements
	rd.simple_reflect(s,lr); // $0$-based reflection of $\lambda-\rho$
	rd.simple_reflect(s,gamma.numerator());
	// keep |coef| unchanged here
	goto restart;
      } // |switch|
    } // |for(s)|
    // if loop terminates, then contribute modified, now final, parameter
    result.emplace_front(sr_gamma(x,lr,gamma),coef);
  drop: {} // when jumping here, proceed without contributing
  }
  while (not to_do.empty());
  return result;
} // |Rep_context::finals_for| StandardRepr

SR_poly Rep_context::expand_final (StandardRepr z) const
{
  auto terms = finals_for(std::move(z));
  poly result;
  for (auto it=terms.begin(); not terms.at_end(it); ++it)
    result.add_term(std::move(it->first),Split_integer(it->second));
  return result;
} // |Rep_context::expand_final|


bool deformation_unit::operator!=(const deformation_unit& another) const
{
  if (sample.x()!=another.sample.x() or
      sample.height()!=another.sample.height() or
      sample.y().data()!=another.sample.y().data())
    return true; // easy tests for difference

  auto& i_tab = rc.involution_table();
  const auto& kgb = rc.kgb();
  InvolutionNbr inv_nr = kgb.inv_nr(sample.x());

  {
    const int_Matrix& theta = i_tab.matrix(inv_nr);
    const auto gamma_diff_num =
      (sample.gamma()-another.sample.gamma()).numerator();
    if (not (theta*gamma_diff_num+gamma_diff_num).is_zero())
      return true; // difference in the free part of $\lambda$ spotted
  }

  const auto& rd = rc.root_datum();
  const auto& g0=sample.gamma();
  const auto& g1=another.sample.gamma();
  const auto& num0 = g0.numerator();
  const auto& num1 = g1.numerator();
  const auto d0 = g0.denominator(), d1 = g1.denominator(); // convert to signed

  { RootNbrSet complex_posroots = rd.posroot_set() & i_tab.complex_roots(inv_nr);
    for (auto it=complex_posroots.begin(); it(); ++it)
      if (i_tab.complex_is_descent(inv_nr,*it))
	if (arithmetic::divide(rd.coroot(*it).dot(num0),d0) !=
	    arithmetic::divide(rd.coroot(*it).dot(num1),d1))
	  return true; // distinct integer part of evaluation poscoroot found
  }
  {
    const RootNbrSet real_posroots = rd.posroot_set() & i_tab.real_roots(inv_nr);
    auto lambda_rho_real2 =
      rc.lambda_rho(sample)*2-rd.twoRho(rd.posroot_set()^real_posroots);
    for (auto it=real_posroots.begin(); it(); ++it)
    {
      const auto& alpha_v= rd.coroot(*it);
      if (alpha_v.dot(lambda_rho_real2)%4!=0) // whether parity at |gamma==0|
      { // when parity at 0, compare integer quotients of evaluations by 2
	if (arithmetic::divide(alpha_v.dot(num0),2*d0) !=
	    arithmetic::divide(alpha_v.dot(num1),2*d1))
	  return true; // distinct integer part of evaluation on poscoroot found
      }
      else // nonparity at 0, so shift division cut-off by 1
	if (arithmetic::divide(alpha_v.dot(num0)+d0,2*d0) !=
	    arithmetic::divide(alpha_v.dot(num1)+d1,2*d1))
	  return true; // distinct integer part of evaluation on poscoroot found
    }
  }

  return false; // if no differences detected, consider |another| as equivalent
}

size_t deformation_unit::hashCode(size_t modulus) const
{
  auto& i_tab = rc.involution_table();
  const auto& kgb = rc.kgb();
  InvolutionNbr inv_nr = kgb.inv_nr(sample.x());

  size_t hash = 17*sample.x() + 89*sample.y().data().to_ulong();
  const int_Matrix& theta = i_tab.matrix(inv_nr);
  const auto& g = sample.gamma();
  const auto& num = g.numerator();
  const auto denom = g.denominator();
  // take into account free part of $\lambda$
  for (auto c : (theta*num+num)/denom) // over temporary |arithmetic::Numer_t|
    hash = 21*hash + c;

  const auto& rd = rc.root_datum();
  { RootNbrSet complex_posroots = rd.posroot_set() & i_tab.complex_roots(inv_nr);
    for (auto it=complex_posroots.begin(); it(); ++it)
      if (i_tab.complex_is_descent(inv_nr,*it))
	hash = 5*hash + arithmetic::divide(rd.coroot(*it).dot(num),denom);
  }
  {
    const RootNbrSet real_posroots = rd.posroot_set() & i_tab.real_roots(inv_nr);
    auto lambda_rho_real2 =
      rc.lambda_rho(sample)*2-rd.twoRho(rd.posroot_set()^real_posroots);
    for (auto it=real_posroots.begin(); it(); ++it)
    {
      const auto& alpha_v= rd.coroot(*it);
      const auto shift = alpha_v.dot(lambda_rho_real2)%4==0 ? denom : 0;
      hash = 7*hash + arithmetic::divide(alpha_v.dot(num)+shift,2*denom);
    }
  }

  return hash&(modulus-1);
}

//				|Rep_table| methods


Rep_table::Rep_table(RealReductiveGroup &G)
: Rep_context(G)
, pool(), alcove_hash(pool)
, reduced_pool(), reduced_hash(reduced_pool)
, K_type_pool(), K_type_hash(K_type_pool)
, KL_poly_pool{KLPol(),KLPol(KLCoeff(1))}, KL_poly_hash(KL_poly_pool)
, poly_pool{ext_kl::Pol(0),ext_kl::Pol(1)}, poly_hash(poly_pool)
, block_list(), place()
{}
Rep_table::~Rep_table() = default;

unsigned short Rep_table::length(StandardRepr sr)
{
  make_dominant(sr); // length should not change in equivalence class
  BlockElt z;
  block_modifier bm; // unused, but obligatory as argument
  auto & block = lookup(sr,z,bm); // construct partial block
  return block.length(z);
}


using Mod_hash_tp = HashTable<StandardReprMod,BlockElt>;

// an auxialry structure needed to hash |unsigned long| to shorter |BlockElt|
struct ulong_entry
{ unsigned long val;
  ulong_entry(unsigned long n) : val(n) {}
  typedef std::vector<ulong_entry> Pooltype;
  bool operator != (ulong_entry x) const { return val!=x.val;  }
  size_t hashCode(size_t modulus) const { return (5*val)&(modulus-1); }
};

// |Rep_table| helper class
class Rep_table::Bruhat_generator
{
  Mod_hash_tp& mod_hash;
  const common_context& ctxt;
  std::vector<simple_list<BlockElt> > predecessors;
public:
  Bruhat_generator (Mod_hash_tp& hash, const common_context& ctxt)
    : mod_hash(hash),ctxt(ctxt), predecessors() {}

  bool in_interval (const StandardReprMod& srm) const
  { return mod_hash.find(srm)!=mod_hash.empty; }
  const simple_list<BlockElt>& covered(BlockElt n) const
  { return predecessors.at(n); }
  void block_below(const StandardReprMod& srm);
}; // |class Rep_table::Bruhat_generator|


void Rep_table::Bruhat_generator::block_below (const StandardReprMod& srm)
{
  if (mod_hash.find(srm)!=mod_hash.empty) // then |srm| was seen earlier
    return; // nothing new

  const auto rank = ctxt.subsys().rank();
  sl_list<BlockElt> pred; // list of elements covered by z
  // invariant: |block_below| has been called for every element in |pred|

  weyl::Generator s; // a complex or real type 1 descent to be found, if exists
  for (s=0; s<rank; ++s)
  {
    std::pair<gradings::Status::Value,bool> stat=ctxt.status(s,srm.x());
    if (not stat.second)
      continue; // ignore imaginary, complex ascent or real (potentially) type 2
    if (stat.first==gradings::Status::Complex)
    { // complex descent
      const StandardReprMod sz = ctxt.cross(s,srm);
      block_below(sz); // recursion
      pred.push_back(mod_hash.find(sz)); // |mod_hash.find| after |block_below|
      break; // we shall add $s$-ascents of predecessors of |sz| below
    }
    else if (stat.first==gradings::Status::Real and ctxt.is_parity(s,srm))
    { // |z| has a type 1 real descent at |s|
      const StandardReprMod sz0 = ctxt.down_Cayley(s,srm);
      const StandardReprMod sz1 = ctxt.cross(s,sz0);
      block_below(sz0); // recursion
      block_below(sz1); // recursion
      pred.push_back(mod_hash.find(sz0)); // |mod_hash.find| after |block_below|
      pred.push_back(mod_hash.find(sz1));
      break; // we shall add $s$-ascents of predecessors of |sz_inx| below
    } // |if (real type 1)|

  } // |for (s)|

  // if above loop performed a |break| it found complex or real type 1 descent
  if (s==rank) // otherwise, the only descents are real type 2, if any
  {
    while (s-->0) // we reverse the loop just because |s==rank| already
      if (ctxt.status(s,srm.x()).first==gradings::Status::Real and
	  ctxt.is_parity(s,srm))
      {
	const auto sz = ctxt.down_Cayley(s,srm);
	block_below(sz); // recursion
	pred.push_back(mod_hash.find(sz));
      }
  }
  else // a complex or real type 1 descent |sz==pred.front()| for |s| was found
  { // add |s|-ascents for elements covered by |sz|
    const auto pred_sz = predecessors.at(pred.front());
    for (auto it = pred_sz.begin(); not pred.at_end(it); ++it)
    {
      const BlockElt p = *it; // sequence number of a predecessor of |sz|
      const StandardReprMod& zp = mod_hash[p];
      std::pair<gradings::Status::Value,bool> stat = ctxt.status(s,zp.x());
      switch (stat.first)
      {
      case gradings::Status::Real: case gradings::Status::ImaginaryCompact:
	break; // nothing to do without ascent
      case gradings::Status::Complex:
	if (not stat.second) // complex ascent
	{
	  const StandardReprMod szp = ctxt.cross(s,zp);
	  block_below(szp); // recursion
	  pred.push_back(mod_hash.find(szp));
	} // |if(complex ascent)
	break;
      case gradings::Status::ImaginaryNoncompact:
	{
	  const StandardReprMod szp = ctxt.up_Cayley(s,zp);
	  block_below(szp); // recursion
	  pred.push_back(mod_hash.find(szp));
	  if (not stat.second) // then nci type 2
	  {
	    const StandardReprMod szp1 = ctxt.cross(s,szp);
	    block_below(szp1); // recursion
	    pred.push_back(mod_hash.find(szp1));
	  }
	}
	break;
      } // |switch(status(s,conj_x))|
    } // |for (it)|
  } // |if (s<rank)|

  const auto h = mod_hash.match(srm); // now generate sequence number for |srm|
  assert(h==predecessors.size()); ndebug_use(h);
  predecessors.push_back(pred.undress()); // store |pred| at |h|
} // |Rep_table::Bruhat_generator::block_below|

// a structure used in |Rep_table::add_block_below| and |Rep_table::add_block|
// for each sub_block, record block pointer, entry element, shift
// DO NOT record iterator into |block_list| which might get invalidated,
// access of iterator to |bp| through |place[h]| will be checked and corrected
struct sub_triple {
  common_block* bp; unsigned long h; RatWeight shift;
  sub_triple(common_block* bp, unsigned long h, RatWeight shift)
    : bp(bp),h(h),shift(shift) {}
};

blocks::common_block& Rep_table::add_block_below
  (const common_context& ctxt, const StandardReprMod& init, BitMap* subset)
{
  assert // we are called to add a block for nothing like what is known before
    (reduced_hash.find(Reduced_param
		       (*this,init,ctxt.integral_nr(),ctxt.attitude()))
     ==reduced_hash.empty);

  StandardReprMod::Pooltype pool;
  Mod_hash_tp hash(pool);
  Bruhat_generator gen(hash,ctxt); // object to help generating Bruhat interval
  gen.block_below(init); // generate Bruhat interval below |srm| into |pool|

  const size_t place_limit = place.size();
  sl_list<sub_triple> sub_blocks;
  for (const auto& elt : pool)
  {
    auto h = match_reduced_hash(elt,ctxt); // fingerprint of parameter family
    if (h==place.size()) // block element has new reduced hash value
      place.emplace_back(bl_it(),-1); // create slot; both fields filled later
    else if (h<place_limit) // then a similar parameter was known
    { // record block pointer and offset of |elt| from its buddy in |sub_blocks|
      common_block* sub = &*place[h].first;
      auto hit = [sub] (const sub_triple& tri)->bool { return tri.bp==sub; };
      if (std::none_of(sub_blocks.begin(),sub_blocks.end(),hit))
      {
	StandardReprMod base = sub->representative(place[h].second);
	sub_blocks.emplace_back(sub,h,offset(elt,base));
      }
    }
  }

  // adapt elements for all |sub_blocks|
  size_t limit = pool.size(); // limit of generated Bruhat interval
  const auto rho = rootdata::rho(root_datum());
  for (auto sub : sub_blocks)
  {
    sub.bp->shift(sub.shift); // make representatives match swallowing block
    for (BlockElt z=0; z<sub.bp->size(); ++z)
      hash.match(sub.bp->representative(z)); // if new, to |pool| beyond |limit|
  }

  sl_list<StandardReprMod> elements(pool.begin(),pool.end()); // working copy

  sl_list<blocks::common_block> temp; // must use temporary singleton
  auto& block =
     temp.emplace_back(ctxt,elements); // construct block and get a reference

  *subset=BitMap(block.size()); // this bitmap will be exported via |subset|
  sl_list<std::pair<BlockElt,BlockEltList> > partial_Hasse_diagram;
  for (unsigned int i=0; i<limit; ++i)
  {
    BlockElt i_z = block.lookup(pool[i]); // index of element |i| in new block
    subset->insert(i_z); // mark |z| as element ot the Bruhat interval
    const auto& covered = gen.covered(i);
    BlockEltList row;
    row.reserve(containers::length(covered));
    for (auto it=covered.begin(); not covered.at_end(it); ++it)
    {
      const BlockElt y = block.lookup(pool[*it]); // get relative number
      assert(y!=UndefBlock);
      row.push_back(y); // store covering relation in |Hasse_diagram|
    }
    partial_Hasse_diagram.emplace_back(i_z,row);
  }

  block.set_Bruhat(std::move(partial_Hasse_diagram));
  // remainder of Hasse diagram will be imported from swallowed sub-blocks

  for (const auto& sub : sub_blocks) // swallow sub-blocks
  {
    auto& sub_block = *sub.bp; // already shifted, so ignore |sub.shift|
    BlockEltList embed; embed.reserve(sub_block.size()); // translation array
    for (BlockElt z=0; z<sub_block.size(); ++z)
    {
      const auto& elt = sub_block.representative(z);
      const BlockElt z_rel = block.lookup(elt);
      assert(z_rel!=UndefBlock);
      embed.push_back(z_rel);
    }

    block.swallow(std::move(sub_block),embed,&KL_poly_hash,&poly_hash);
    block_erase(place[sub.h].first); // iterator argument might be corrected
  }


  // only after the |block_erase| upheavals is it safe to link in the new block
  // also, it can only go to the end of the list, to not invalidate any iterators
  const auto new_block_it=block_list.end(); // iterator for new block elements
  block_list.splice(new_block_it,temp,temp.begin()); // link in |block| at end

  for (BlockElt z=block.size(); z-->0; ) // decreasing: least |z| wins below
  { // by using reverse iteration, least elt with same |h| defines |place[h]|
    auto h = find_reduced_hash(block,z);
    assert(h!=reduced_hash.empty);
    place[h] = std::make_pair(new_block_it,z); // extend or replace
  }
  return block;
} // |Rep_table::add_block_below|

// erase node in |block_list| after |pos|, avoiding dangling iterators in |place|
void Rep_table::block_erase (bl_it pos)
{
  assert (not block_list.at_end(pos));
  const auto next_pos = std::next(pos);
  if (not block_list.at_end(next_pos))
  { // then make sure in |place| instances of |next_pos| are replaced by |pos|
    const auto& block = *next_pos;
    for (BlockElt z=0; z<block.size(); ++z)
    {
      auto seq = find_reduced_hash(block,z);
      assert(seq<place.size()); // all elements in |block_list| must have |place|
      if (place[seq].first==next_pos) // could be false if |block| was swallowed
	place[seq].first=pos; // replace iterator that is about to be invalidated
    }
  }
  block_list.erase(pos);
} // |Rep_table::block_erase|

unsigned long Rep_table::add_block
  (const StandardReprMod& srm, const common_context& ctxt)
{
  BlockElt srm_in_block; // will hold position of |srm| within that block
  sl_list<blocks::common_block> temp; // must use temporary singleton
  auto& block = temp.emplace_back(ctxt,srm,srm_in_block); // build full block

  const auto rho = rootdata::rho(root_datum());
  const size_t place_limit = place.size();

  sl_list<sub_triple> sub_blocks;
  for (BlockElt z=0; z<block.size(); ++z)
  {
    auto seq = match_reduced_hash(block,z);
    if (seq==place.size()) // block element has new reduced hash value
      place.emplace_back(bl_it(),z); // create slot; iterator filled later
    else if (seq<place_limit)
    {
      common_block* sub = &*place[seq].first;
      auto hit = [sub] (const sub_triple& tri)->bool { return tri.bp==sub; };
      if (std::none_of(sub_blocks.begin(),sub_blocks.end(),hit))
      {
	StandardReprMod base = sub->representative(place[seq].second);
	sub_blocks.emplace_back(sub,seq,offset(block.representative(z),base));
      }
    }
  }

  // swallow |embeddings|, and remove them from |block_list|
  for (const auto& sub : sub_blocks) // swallow sub-blocks
  {
    sub.bp->shift(sub.shift); // make representatives match swallowing block
    auto& sub_block = *sub.bp;
    BlockEltList embed; embed.reserve(sub_block.size()); // translation array
    for (BlockElt z=0; z<sub_block.size(); ++z)
    {
      const BlockElt z_rel = block.lookup(sub_block.representative(z));
      assert(z_rel!=UndefBlock); // our block is full, lookup should work
      embed.push_back(z_rel);
    }

    block.swallow(std::move(sub_block),embed,&KL_poly_hash,&poly_hash);
    block_erase(place[sub.h].first);
  }

  // only after |block_erase| upheavals is it safe to link in the new block
  // also, it can only go to the end of the list, to not invalidate any iterators
  const auto new_block_it=block_list.end(); // iterator for new block elements
  block_list.splice(new_block_it,temp,temp.begin()); // link in |block| at end

  // now make sure for all |elements| that |place| fields are set for new block
  for (BlockElt z=0; z<block.size(); ++z)
  {
    auto h = find_reduced_hash(block,z);
    place[h].first = new_block_it;
    place[h].second = z;
  }
  return find_reduced_hash(srm,ctxt);
}// |Rep_table::add_block|

blocks::common_block& Rep_table::lookup_full_block
  (StandardRepr& sr,BlockElt& z, block_modifier& bm)
{
  make_dominant(sr); // without this we would not be in any valid block
  auto srm = StandardReprMod::mod_reduce(*this,sr); // modular |z|
  common_context ctxt(*this,srm.gamma_lambda());
  Reduced_param rp(*this,srm,ctxt.integral_nr(),ctxt.attitude());

  auto h = reduced_hash.find(rp); // look up modulo $X^*+integral^\perp$
  if (h==reduced_hash.empty or not place[h].first->is_full()) // then we must
    h=add_block(srm,ctxt); // generate new full block (maybe swallow older ones)
  assert(h<place.size() and place[h].first->is_full());

  auto& block = *place[h].first;
  bm.shift = offset(srm,block.representative(z = place[h].second));
  return block;

} // |Rep_table::lookup_full_block|

blocks::common_block& Rep_table::lookup
  (StandardRepr& sr,BlockElt& which, block_modifier& bm)
{
  normalise(sr); // gives a valid block, and smallest partial block
  auto srm = StandardReprMod::mod_reduce(*this,sr); // modular |z|
  common_context ctxt(*this,srm.gamma_lambda());
  Reduced_param rp(*this,srm,ctxt.integral_nr(),ctxt.attitude());

  assert(reduced_hash.size()==place.size()); // should be in sync at this point
  auto h = reduced_hash.find(rp); // look up modulo $X^*+integral^\perp$
  if (h!=reduced_hash.empty) // then we have found our family of blocks
  {
    assert(h<place.size());
    auto& block = *place[h].first;
    which = place[h].second;
    assert(block.representative(which).x()==srm.x()); // check minimum of sanity
    bm.shift = offset(srm,block.representative(which));
    return block; // use block of related |StandardReprMod| as ours
  }
  BitMap subset;
  auto& block = add_block_below(ctxt,srm,&subset);
  which = last(subset);
  assert(Reduced_param
	 (*this,block.representative(which),ctxt.integral_nr(),ctxt.attitude())
	 ==rp);
  bm.shift = RatWeight(rank());
  return block;
} // |Rep_table::lookup|

/* In the following type the second component is a multiplicity so we are
   dealing with a sparse representation of polynomials with |BlockElt|
   exponents. This is intended for the expansion of non-final elements into
   final ones; we expect few, relative to the block size, terms per polynomial
*/
using BlockElt_term =  std::pair<BlockElt,int>;
using BlockElt_pol = sl_list<BlockElt_term>;

// addition polynomials whose elements are sorted by increasing |BlockElt|
BlockElt_pol combine (BlockElt_pol a, BlockElt_pol b) // by value
{ a.merge(std::move(b),
	  [](const BlockElt_term& x, const BlockElt_term& y)
	    { return x.first<y.first; });
  // now any like terms are neigbours, combine them whenever this occurs
  for (auto it=a.begin(); not a.at_end(it); ++it)
  {
    auto it1=std::next(it);
    if (not a.at_end(it1) and it->first==it1->first)
    {
      it->second += it1->second; // accumulate towards first ot the terms
      a.erase(it1); // then erase the second; next |++it| correctly executes
    }
  }
  return a;
}

BlockElt_pol flip (int sign, BlockElt_pol list) // by value
{ if (sign!=1)
    for (auto& p : list)
      p.second *= sign;
  return list;
}

// Expand block elements up to |y| into final ones for the |singular| system.
// This version is for a |common_block|, so no twist is involved in expansion
std::vector<BlockElt_pol> contributions
  (blocks::common_block& block, RankFlags singular, BlockElt y)
{
  std::vector<BlockElt_pol> result(y+1); // initally every |result[z]| is empty
  for (BlockElt z=0; z<=y; ++z) // compute in |result[z]| the finals for |z|
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
} // |contributions|


// Expand block elements up to |y| into final ones for |singular| system.
// Being for an |ext_block|, the expansion involves signs (twisted case)
std::vector<BlockElt_pol> contributions
  (const ext_block::ext_block& eblock, RankFlags singular_orbits,
   BlockElt limit) // where to stop computing contributions
{
  std::vector<BlockElt_pol> result(limit); // each|result[z]| is initially empty
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

size_t Rep_table::find_reduced_hash(const common_block& block, BlockElt z) const
{ return reduced_hash.find(Reduced_param
			   {*this, block.representative(z),
			    block.integral_nr(),block.attitude()} );
}
size_t Rep_table::find_reduced_hash
  (const StandardReprMod& srm, const common_context& c) const
{ return reduced_hash.find(Reduced_param
			   {*this,srm, c.integral_nr(),c.attitude()} ); }

size_t Rep_table::match_reduced_hash(const common_block& block, BlockElt z)
{ return reduced_hash.match(Reduced_param
			   {*this, block.representative(z),
			    block.integral_nr(),block.attitude()} );
}
size_t Rep_table::match_reduced_hash
  (const StandardReprMod& srm, const common_context& c)
{ return reduced_hash.match(Reduced_param
			    {*this,srm, c.integral_nr(),c.attitude()} ); }

sl_list<std::pair<StandardRepr,int> > Rep_table::deformation_terms
  ( blocks::common_block& block, const BlockElt y,
    const RatWeight& diff, const RatWeight& gamma)
{ assert(y<block.size()); // and |y| is final, see |assert| below

  sl_list<std::pair<StandardRepr,int> > result;
  if (block.length(y)==0)
    return result; // easy case, null result

  std::vector<BlockElt_pol> contrib =
    contributions(block,block.singular(gamma),y);
  sl_list<BlockElt> finals;
  for (BlockElt z=0; z<contrib.size(); ++z)
    if (not contrib[z].empty() and contrib[z].front().first==z)
      finals.push_front(z); // accumulate in reverse order

  assert(not finals.empty() and finals.front()==y); // do not call for non-final
  const kl::KL_table& kl_tab =
    block.kl_tab(&KL_poly_hash,y+1); // fill silently up to |y|

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
      const kl::KLPol& pol = kl_tab.KL_pol(x,z); // regular KL polynomial
      int eval = 0;
      for (polynomials::Degree d=pol.size(); d-->0; )
	eval = static_cast<int>(pol[d]) - eval; // evaluate at $q = -1$
      if (eval==0)
	continue; // polynomials with $-1$ as root do not contribute; skip
      if ((block.length(z)-block.length(x))%2!=0) // when |l(z)-l(x)| odd
	eval=-eval; // flip sign (do alternating sum of KL column at |-1|)
      for (auto jt=contrib[x].wcbegin(); not contrib[x].at_end(jt); ++jt)
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
    const unsigned int orient_y = orientation_number(block.sr(y,diff,gamma));

    auto it=finals.begin();
    for (const int c : acc) // accumulator |acc| runs parallel to |finals|
    {
      const auto z = *it; ++it;
      if (c!=0) // test must follow |++it| !
      {
	const auto sr_z = block.sr(z,diff,gamma);
	auto coef = c*arithmetic::exp_i(orient_y-orientation_number(sr_z));
	result.emplace_back(sr_z,coef);
      }
    }
    assert(it==finals.end());
  }

  return result;
} // |deformation_terms|, common block version

sl_list<SR_poly::value_type> Rep_table::block_deformation_to_height
  (StandardRepr p, SR_poly& queue, level height_bound)
{
  BlockElt start; block_modifier bm;
  auto& block = lookup_full_block(p,start,bm);
  const RatWeight& diff = bm.shift;
  const auto& gamma = p.gamma();
  auto dual_block = blocks::Bare_block::dual(block);
  kl::KL_table& kl_tab = dual_block.kl_tab(nullptr,1);
  // create KL table only minimally filled

  // now fill remainder up to height; prepare for restriction to this subset
  BitMap retained(block.size());
  sl_list<std::pair<const StandardRepr,Split_integer> > result;
  for (BlockElt z=0; z<block.size(); ++z)
  { StandardRepr q = sr(block.representative(z),diff,gamma);
    if (retained.set_to(z,q.height()<=height_bound))
    {
      auto it = queue.find(q);
      if (it==queue.end())
	result.emplace_back(std::move(q),Split_integer(0));
      else
      {
	result.emplace_back(std::move(q),it->second); // push |(q,queue[q])|
	queue.erase(it); // and remove that term it from the |queue|
      }
    }
    else
      kl_tab.plug_hole(block.size()-1-z); // not interested in this parameter
  }
  kl_tab.fill(); // fill whole table; we might go beyond |size()-1-start|

  int_Vector value_at_minus_1;
  value_at_minus_1.reserve(kl_tab.pol_store().size());
  for (auto& entry : kl_tab.pol_store())
  { int val = 0;
    for (unsigned d=entry.size(); d-->0;)
      val = static_cast<int>(entry[d])-val; // Horner evaluate polynomial at -1
    value_at_minus_1.push_back(val);
  }

  const RankFlags singular = block.singular(gamma); // singular simple coroots
  auto it = result.begin();
  for (auto elt: retained) // don't increment |it| here
    if (block.survives(elt,singular))
      ++it;
    else
    { retained.remove(elt);
      assert(it->second.is_zero()); // |queue| cannot have non final parameter
      result.erase(it); // drop term; |it| now points to next term
    }
  result.reverse(); // we shall need to travers elements downwards in |block|

  // viewed from |block|, the |kl_tab| is lower triangular
  // build its transpose, restricted to |retained|, and evaluated at $q=-1$
  int_Matrix Q_mat (result.size()); // initialise to identity matrix
  unsigned int i=0,j;
  unsigned int const top=block.size()-1;
  for (auto it=retained.begin(); it(); ++it,++i)
    for (auto jt=(j=i+1,std::next(it)); jt(); ++jt,++j)
      Q_mat(i,j) = value_at_minus_1[kl_tab.KL_pol_index(top-*jt,top-*it)];

  int_Matrix signed_P = Q_mat.inverse();
  BitMap odd_length(signed_P.n_rows());
  { unsigned int i=0;
    for (const BlockElt z : retained)
      odd_length.set_to(i++,block.length(z)%2!=0);
  }

  matrix::Vector<Split_integer> coef(signed_P.n_rows());
  auto pos = result.size()-1;
  for (auto it=result.begin(); not result.at_end(it); ++it,--pos)
    if (not it->second.is_zero())
    { coef.assign(pos,Split_integer(0));
      // compute into |coef the product of columns of |signed_P| with parity
      // opposite to |pos| with corresponding entries of column |pos| of |Q_mat|
      for (auto j : (odd_length.isMember(pos) ? ~odd_length : odd_length))
      { if (j>=pos)
	  break;
	for (unsigned i=0; i<=j; ++i)
	  coef[i] += signed_P(i,j)*Q_mat(j,pos);
      }
      { auto our_orient = orientation_number(it->first);
	auto j = pos-1;
	for (auto jt = std::next(it); not result.at_end(jt); ++jt,--j)
	{ coef[j] *= it->second;
	  auto diff = our_orient - orientation_number(jt->first);
	  assert(diff%2==0);
	  coef[j].times_1_s();
	  if (diff%4!=0)
	    coef[j].times_s(); // equivalently |coef[j].negate|
	  jt->second += coef[j]; // contribute coefficient to result
	}
      }
  }

  return result;
}

// compute and return sum of KL polynomials at $s$ for final parameter |sr|
SR_poly Rep_table::KL_column_at_s(StandardRepr sr) // |sr| must be final
{
  normalise(sr); // implies that |sr| it will appear at the top of its own block
  assert(is_final(sr));

  BlockElt z; block_modifier bm;
  auto& block = lookup(sr,z,bm);
  const RatWeight& diff = bm.shift;
  assert((involution_table().matrix(kgb().inv_nr(block.x(z)))*diff+diff)
	 .is_zero());

  const auto& gamma=sr.gamma();
  std::vector<BlockElt_pol> contrib =
    contributions(block,block.singular(gamma),z);
  assert(contrib.size()==z+1 and contrib[z].front().first==z);

  const kl::KL_table& kl_tab =
    block.kl_tab(&KL_poly_hash,z+1); // fill silently up to |z|

  SR_poly result;
  auto z_length=block.length(z);
  for (BlockElt x=z+1; x-->0; )
  {
    const kl::KLPol& pol = kl_tab.KL_pol(x,z); // regular KL polynomial
    if (pol.is_zero())
      continue;
    Split_integer eval(0);
    for (polynomials::Degree d=pol.size(); d-->0; )
      eval.times_s() += static_cast<int>(pol[d]); // evaluate at $q = s$
    // no test here, nonzero KL polynomials have nonzero evaluation at $q=1$

    if ((z_length-block.length(x))%2!=0) // when |l(z)-l(x)| odd
      eval.negate(); // flip sign (do alternating sum of KL column at |s|)
    for (const auto& pair : contrib[x])
      result.add_term(block.sr(pair.first,diff,gamma),eval*pair.second);
  }

  return result;
} // |Rep_table::KL_column_at_s|

// compute and return column of KL table for final parameter |sr|
simple_list<std::pair<BlockElt,kl::KLPol> >
  Rep_table::KL_column(StandardRepr sr) // |sr| must be final
{
  assert(is_final(sr));

  BlockElt z; block_modifier bm;
  auto& block = lookup(sr,z,bm);

  const kl::KL_table& kl_tab =
    block.kl_tab(&KL_poly_hash,z+1); // fill silently up to |z|

  simple_list<std::pair<BlockElt,kl::KLPol> > result;
  for (BlockElt x=z+1; x-->0; )
  {
    const kl::KLPol& pol = kl_tab.KL_pol(x,z); // regular KL polynomial
    if (not pol.is_zero())
      result.emplace_front(x,pol);
  }

  return result;
} // |Rep_table::KL_column|


const K_type_poly& Rep_table::deformation(StandardRepr z)
// that |z| is dominant and final is a precondition assured in the recursion
// for more general |z|, do the preconditioning outside the recursion
{
  assert(is_final(z));
  if (z.gamma().denominator() > (1LL<<rank()))
    z = weyl::alcove_center(*this,z);
  RatNumList rp=reducibility_points(z);

  deformation_unit zn(*this,z);
  { // look up if deformation formula for |z| is already known and stored
    unsigned long h=alcove_hash.find(zn);
    if (h!=alcove_hash.empty and pool[h].has_deformation_formula())
      return pool[h].def_formula();
  }

  K_type_poly result {std::less<K_type_nr>()};
  {
    auto base_pol = finals_for(scale_0(z)); // a $\Z$-linear combin. of K-types
    for (auto it=base_pol.begin(); not base_pol.at_end(it); ++it)
    {
      K_type_nr h = K_type_hash.match(std::move(it->first));
      result.add_term(h,Split_integer(it->second)); // purely integer coefficient
    }
  }

  for (unsigned i=rp.size(); i-->0; )
  {
    auto zi = scale(z,rp[i]);
    deform_readjust(zi); // necessary to ensure the following |assert| will hold
    assert(is_final(zi)); // ensures that |deformation_terms| won't refuse
    BlockElt new_z; block_modifier bm;
    auto& block = lookup(zi,new_z,bm);
    const RatWeight& diff = bm.shift;
    assert((involution_table().matrix(kgb().inv_nr(block.x(new_z)))*diff+diff)
	   .is_zero());
    auto dt = deformation_terms(block,new_z,diff,zi.gamma());
    for (auto& term : dt)
    {
      const auto& def = deformation(term.first); // recursion
      result.add_multiple
	(def,Split_integer(term.second,-term.second)); // $(1-s)*c$
    }
  }

  const auto h = alcove_hash.match(std::move(zn)); // allocate a slot in |pool|
  return pool[h].set_deformation_formula(std::move(result).flatten());
} // |Rep_table::deformation|


// basic computation of twisted KL column sum, no tabulation of the result
SR_poly twisted_KL_sum
( ext_block::ext_block& eblock, BlockElt y, const blocks::common_block& parent,
  const RatWeight& diff,
  const RatWeight& gamma) // infinitesimal character, possibly singular
{
  // compute cumulated KL polynomials $P_{x,y}$ with $x\leq y$ survivors

  // start with computing KL polynomials for the entire block
  std::vector<ext_kl::Pol> pool;
  ext_KL_hash_Table hash(pool,4);
  ext_kl::KL_table twisted_KLV(eblock,&hash);
  twisted_KLV.fill_columns(y+1); // fill table up to |y| inclusive

  // make a copy of |pool| in which polynomials have been evaluated as |s|
  std::vector<Split_integer> pool_at_s; pool_at_s.reserve(pool.size());
  for (unsigned i=0; i<pool.size(); ++i)
    if (pool[i].is_zero())
      pool_at_s.push_back(Split_integer(0,0));
    else
    { const auto& P = pool[i];
      auto d=P.degree();
      Split_integer eval(P[d]);
      while (d-->0)
	eval.times_s()+=static_cast<int>(P[d]);
      pool_at_s.push_back(eval);
    }

  const RootDatum& rd = parent.root_datum();
  const RootNbrList int_simples = parent.int_simples();
  RankFlags singular_orbits; // flag singulars among orbits
  for (weyl::Generator s=0; s<eblock.rank(); ++s)
    singular_orbits.set(s,
	     gamma.dot(rd.coroot(int_simples[eblock.orbit(s).s0]))==0);

  auto contrib = contributions(eblock,singular_orbits,y+1);

  SR_poly result;
  unsigned int parity = eblock.length(y)%2;
  for (BlockElt x=0; x<=y; ++x)
  { const auto& p = twisted_KLV.KL_pol_index(x,y);
    Split_integer eval = p.second ? -pool_at_s[p.first] : pool_at_s[p.first];
    if (eblock.length(x)%2!=parity) // flip sign at odd length difference
      eval = -eval;
    for (const auto& pair : contrib[x])
      result.add_term(parent.sr(eblock.z(pair.first),diff,gamma),
		      eval*pair.second);
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
  BlockElt entry;
  common_context ctxt(rc,zm.gamma_lambda());
  blocks::common_block block(ctxt,zm,entry); // build full block
  const RatWeight diff(rc.rank()); // zero: we custom-built our |block| above
  auto eblock = block.extended_block(delta);
  return twisted_KL_sum(eblock,eblock.element(entry),block,diff,z.gamma());
} // |twisted_KL_column_at_s|

// look up or compute and return the alternating sum of twisted KL polynomials
// at the inner class involution for final parameter |z|, evaluated at $q=s$
SR_poly Rep_table::twisted_KL_column_at_s(StandardRepr sr)
  // |z| must be inner-class-twist-fixed, nonzero and final
{
  normalise(sr);
  assert(is_final(sr) and sr==inner_twisted(sr));
  BlockElt y0; block_modifier bm;
  auto& block = lookup(sr,y0,bm);
  const RatWeight& diff = bm.shift;
  block.shift(diff);
  auto& eblock = block.extended_block(&poly_hash);
  block.shift(-diff);

  RankFlags singular=block.singular(sr.gamma());
  RankFlags singular_orbits; // flag singulars among orbits
  for (weyl::Generator s=0; s<eblock.rank(); ++s)
    singular_orbits.set(s,singular[eblock.orbit(s).s0]);

  const BlockElt y = eblock.element(y0);
  auto contrib = contributions(eblock,singular_orbits,y+1);

  sl_list<BlockElt> finals; // these have numbering for |eblock|!
  for (BlockElt z=0; z<=y; ++z)
    if (not contrib[z].empty() and contrib[z].front().first==z)
      finals.push_front(z); // accumulate in reverse order

  const auto& kl_tab = eblock.kl_table(y+1,&poly_hash);

  SR_poly result;
  const auto& gamma=sr.gamma();
  const RatWeight gamma_rho = gamma-rho(block.root_datum());
  auto y_length=block.length(y0);

  for (BlockElt x=y+1; x-->0; )
  {
    const auto& pol = kl_tab.P(x,y); // twisted KL polynomial
    if (pol.is_zero())
      continue;
    Split_integer eval(0);
    for (polynomials::Degree d=pol.size(); d-->0; )
      eval.times_s() += pol[d]; // evaluate at $q = s$

    if ((y_length-block.length(eblock.z(x)))%2!=0) // when |l(y)-l(x)| odd
      eval.negate(); // flip sign (do alternating sum of KL column at |s|)
    for (const auto& pair : contrib[x])
      result.add_term(block.sr(eblock.z(pair.first),diff,gamma),
		      eval*pair.second);
  }

  return result;
} // |Rep_table::twisted_KL_column_at_s|

sl_list<std::pair<StandardRepr,int> >
Rep_table::twisted_deformation_terms
    (blocks::common_block& block, ext_block::ext_block& eblock,
     BlockElt y, // in numbering of |block|, not |eblock|
     RankFlags singular_orbits, const RatWeight& diff, const RatWeight& gamma)
{
  assert(eblock.is_present(y));
  const BlockElt y_index = eblock.element(y);

  sl_list<std::pair<StandardRepr,int> > result;
  if (block.length(y)==0)
    return result; // easy case, null result

  auto contrib = repr::contributions(eblock,singular_orbits,y_index+1);
  sl_list<BlockElt> finals; // these have numbering for |eblock|!
  for (BlockElt z=0; z<contrib.size(); ++z)
    if (not contrib[z].empty() and contrib[z].front().first==z)
      finals.push_front(z); // accumulate in reverse order

  const auto& kl_tab = eblock.kl_table(y_index+1,&poly_hash);

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
      for (auto jt=contrib[x].wcbegin(); not contrib[x].at_end(jt); ++jt)
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
    const unsigned int orient_y = orientation_number(block.sr(y,diff,gamma));

    auto it=acc.begin();
    for (const int f : finals) // accumulator |acc| runs parallel to |finals|
    {
      const int c = *it++;
      if (c==0)
	continue;
      const auto sr_z =
	block.sr(eblock.z(f),diff,gamma); // renumber |f| to |block|

      auto coef = c*arithmetic::exp_i(orient_y-orientation_number(sr_z));
      result.emplace_back(sr_z,coef);
    }
    assert(it==acc.end());
  }

  return result;
} // |twisted_deformation_terms(blocks::common_block&,...)|

#if 0
SR_poly Rep_table::twisted_deformation_terms (unsigned long sr_hash)
{ // the |StandardRepr| |hash[sr_hash]| is necessarily delta-fixed and final
  SR_poly result;
  SR_poly remainder(hash[sr_hash]);
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

const K_type_poly& Rep_table::twisted_deformation(StandardRepr z, bool& flip)
{
  assert(is_final(z));
  if (z.gamma().denominator() > (1LL<<rank()))
    z = weyl::alcove_center(*this,z);
  const auto& delta = inner_class().distinguished();

  RatNumList rp=reducibility_points(z);
  flip = false; // ensure no flip is recorded when shrink wrapping is not done
  if (not rp.empty() and rp.back()!=RatNum(1,1))
  { // shrink wrap toward $\nu=0$ to get useful parameter to store result for
    const RatNum f=rp.back();
    std::pair<StandardRepr,bool> p =
      ext_block::scaled_extended_finalise(*this,z,delta,f);
    z = p.first; flip=p.second;
    for (auto& a : rp)
      a/=f; // rescale reducibility points to new parameter |z|
    assert(rp.back()==RatNum(1,1)); // should make first reduction at |z|
    // here we continue, with |flip| recording whether we already flipped
  }

  deformation_unit zu(*this,z);
  { // if formula for |z| is stored, return it; caller multiplies by |s^flip|
    const auto h=alcove_hash.find(zu);
    if (h!=alcove_hash.empty and pool[h].has_twisted_deformation_formula())
      return pool[h].twisted_def_formula();
  }

  K_type_poly result { std::less<K_type_nr>() };
  { // initialise |result| to restriction of |z| expanded to finals
    auto z_K = ext_block::extended_restrict_to_K(*this,z,delta);
    // the following loop reorders terms by |std::less<K_type_nr>|
    for (auto&& term : z_K) // convert |K_repr::K_type_pol| to |K_type_poly|
      result.add_term(K_type_hash.match(std::move(term.first)),term.second);
  }

  // compute the deformation terms at all reducibility points
  for (unsigned i=rp.size(); i-->0; )
  {
    std::pair<StandardRepr,bool> p =
      ext_block::scaled_extended_finalise(*this,z,delta,rp[i]);
    StandardRepr zi = std::move(p.first);
    const bool flip_p = p.second;
    BlockElt index; block_modifier bm;
    auto& block = lookup(zi,index, bm);
    const RatWeight& diff = bm.shift;
    block.shift(diff); // adapt representatives for extended block construction
    assert(block.representative(index)==
	   StandardReprMod::mod_reduce(*this,zi));
    auto& eblock = block.extended_block(&poly_hash);
    block.shift(-diff);

    RankFlags singular = block.singular(zi.gamma());
    RankFlags singular_orbits; // flag singulars among orbits
    for (weyl::Generator s=0; s<eblock.rank(); ++s)
      singular_orbits.set(s,singular[eblock.orbit(s).s0]);

    auto terms = twisted_deformation_terms(block,eblock,index,
					   singular_orbits,diff,zi.gamma());
    for (auto&& term : terms)
    { bool flip_def;
      const auto& def =
	twisted_deformation(std::move(term.first),flip_def); // recursion
      result.add_multiple
	(def,
	 flip_p!=flip_def ? Split_integer(-term.second,term.second)
			  : Split_integer(term.second,-term.second)
	 );
    }
  }

  const auto h = alcove_hash.match(std::move(zu));  // find or allocate a slot

  return pool[h].set_twisted_deformation_formula(std::move(result).flatten());

} // |Rep_table::twisted_deformation (StandardRepr z)|

K_repr::K_type_pol export_K_type_pol(const Rep_table& rt,const K_type_poly& P)
{
  K_repr::K_type_pol::poly result; // an instance of |std::vector|
  result.reserve(P.size());
  for (const auto& term : P)
    result.emplace_back(rt.stored_K_type(term.first),term.second);
  return { std::move(result), true }; // sort, convert to |K_repr::K_type_poly|
}

//			|common_context| methods

common_context::common_context (const Rep_context& rc, const RatWeight& gamma)
: rep_con(rc)
, int_sys_nr()
, w()
, id_it(rc.inner_class().int_item(gamma,int_sys_nr,w)) // sets |w|
, sub(id_it.int_system(w)) // transform by |w|, and store image subsystem
{} // |common_context::common_context|


std::pair<gradings::Status::Value,bool>
  common_context::status(weyl::Generator s, KGBElt x) const
{
  const auto& conj = sub.to_simple(s); // word in full system
  KGBElt conj_x = kgb().cross(conj,x);
  const auto t=sub.simple(s);
  const auto stat = kgb().status(t,conj_x);
  return std::make_pair(stat,
			stat==gradings::Status::Real
			? kgb().isDoubleCayleyImage(t,conj_x) // real type 1
			: stat==gradings::Status::Complex
			? kgb().isDescent(t,conj_x)
			: conj_x!=kgb().cross(t,conj_x)); // nc imaginary type 1
}

StandardReprMod common_context::cross
    (weyl::Generator s, const StandardReprMod& z) const
{
  const auto& full_datum = full_root_datum();
  const auto& refl = sub.reflection(s); // reflection word in full system
  const KGBElt new_x = kgb().cross(refl,z.x());
  RatWeight gamma_lambda = rc().gamma_lambda(z);

  const auto& i_tab = involution_table();
  RootNbrSet pos_neg = pos_to_neg(full_datum,refl);
  pos_neg &= i_tab.real_roots(kgb().inv_nr(z.x())); // only real roots for |z|
  gamma_lambda -= root_sum(full_datum,pos_neg); // correction for $\rho_r$'s
  subsys().simple_reflect(s,gamma_lambda.numerator()); // integrally simple
  return repr::StandardReprMod::build(rc(),new_x,gamma_lambda);
}

// whether integrally simple root |s|, supposed real, is parity for |z|
bool common_context::is_parity
    (weyl::Generator s, const StandardReprMod& z) const
{
  const auto& full_datum = full_root_datum();
  const auto& i_tab = involution_table();
  const auto& real_roots = i_tab.real_roots(kgb().inv_nr(z.x()));
  assert(real_roots.isMember(sub.parent_nr_simple(s)));
  const Coweight& alpha_hat = subsys().simple_coroot(s);
  const int eval = rc().gamma_lambda(z).dot(alpha_hat);
  const int rho_r_corr = alpha_hat.dot(full_datum.twoRho(real_roots))/2;
  return (eval+rho_r_corr)%2!=0;
}

StandardReprMod common_context::down_Cayley
    (weyl::Generator s, const StandardReprMod& z) const
{
  assert(is_parity(s,z)); // which also asserts that |z| is real for |s|
  const auto& full_datum = full_root_datum();
  const auto& conj = sub.to_simple(s); // word in full system
  const KGBElt conj_x = kgb().cross(conj,z.x());
  const KGBElt new_x =
    kgb().cross(kgb().inverseCayley(sub.simple(s),conj_x).first,conj);
  RatWeight gamma_lambda = rc().gamma_lambda(z); // will be shifted below

  const auto& i_tab = involution_table();
  RootNbrSet pos_neg = pos_to_neg(full_datum,conj);
  RootNbrSet real_flip = i_tab.real_roots(kgb().inv_nr(z.x()));
  real_flip ^= i_tab.real_roots(kgb().inv_nr(new_x));
  pos_neg &= real_flip; // posroots that change real status and map to negative
  gamma_lambda += root_sum(full_datum,pos_neg); // correction of $\rho_r$'s
  return repr::StandardReprMod::build(rc(),new_x,gamma_lambda);
}

StandardReprMod common_context::up_Cayley
    (weyl::Generator s, const StandardReprMod& z) const
{
  const auto& full_datum = full_root_datum();
  const auto& conj = sub.to_simple(s); // word in full system
  const KGBElt conj_x = kgb().cross(conj,z.x());
  assert(kgb().status(sub.simple(s),conj_x)==
	 gradings::Status::ImaginaryNoncompact);
  const auto new_x = kgb().cross(kgb().cayley(sub.simple(s),conj_x),conj);
  RatWeight gamma_lambda = rc().gamma_lambda(z); // will be shifted below

  const auto& i_tab = involution_table();
  const RootNbrSet& upstairs_real_roots = i_tab.real_roots(kgb().inv_nr(new_x));
  RootNbrSet real_flip = upstairs_real_roots;
  real_flip ^= i_tab.real_roots(kgb().inv_nr(z.x())); // remove downstairs reals

  RootNbrSet pos_neg = pos_to_neg(full_datum,conj);
  pos_neg &= real_flip; // posroots that change real status and map to negative
  gamma_lambda += root_sum(full_datum,pos_neg); // correction of $\rho_r$'s

  // correct in case the parity condition fails for our raised |gamma_lambda|
  const Coweight& alpha_hat = subsys().simple_coroot(s);
  const int rho_r_corr = // integer since alpha is among |upstairs_real_roots|
    alpha_hat.dot(full_datum.twoRho(upstairs_real_roots))/2;
  const int eval = gamma_lambda.dot(alpha_hat);
  if ((eval+rho_r_corr)%2==0) // parity condition says it should be 1
    gamma_lambda += RatWeight(subsys().simple_root(s),2); // add half-alpha

  return repr::StandardReprMod::build(rc(),new_x,gamma_lambda);
}


Weight Rep_context::to_simple_shift
  (InvolutionNbr theta, InvolutionNbr theta_p, RootNbrSet S) const
{ const InvolutionTable& i_tab = involution_table();
  S &= (i_tab.real_roots(theta) ^i_tab.real_roots(theta_p));
  return root_sum(root_datum(),S);
}


//			|Ext_rep_context| methods

Ext_rep_context::Ext_rep_context
  (const Rep_context& rc, const WeightInvolution& delta)
    : rep_con(rc)
    , d_delta(delta)
    , pi_delta(rc.root_datum().rootPermutation(delta))
    , delta_fixed_roots(fixed_points(pi_delta))
    , twist()
{
  const RootDatum& rd = rc.root_datum();
  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
    twist[s] = rd.simpleRootIndex(delta_of(rd.simpleRootNbr(s)));

} // |Ext_rep_context::Ext_rep_context|

Ext_rep_context::Ext_rep_context (const repr::Rep_context& rc)
  : Ext_rep_context(rc, rc.inner_class().distinguished()) {}


bool Ext_rep_context::is_very_complex
  (InvolutionNbr theta, RootNbr alpha) const
{ const auto& i_tab = rep_con.involution_table();
  const auto& rd = root_datum();
  assert (rd.is_posroot(alpha)); // this is a precondition
  auto image = i_tab.root_involution(theta,alpha);
  make_positive(rd,image);
  return image!=alpha and image!=delta_of(alpha);
}

/*
  For the conjugation to simple scenario, we compute a set of positive roots
  that become negative under an element of $W^\delta$ that makes the
  integrally-simple root(s) in question simple. From this set |S|, and the
  involutions at both ends of the link in the block, the function |shift_flip|
  computes whether an additional flip is to be added to the link.

  This comes from an action of |delta| on a certain top wedge product of
  root spaces, and the formula below tells whether that action is by $-1$.
*/
bool Ext_rep_context::shift_flip
  (InvolutionNbr theta, InvolutionNbr theta_p, RootNbrSet S) const
{ S.andnot(delta_fixed()); // $\delta$-fixed roots won't contribute

  unsigned count=0; // will count 2-element |delta|-orbit elements
  for (auto it=S.begin(); it(); ++it)
    if (is_very_complex(theta,*it) != is_very_complex(theta_p,*it) and
	not root_datum().sum_is_root(*it,delta_of(*it)))
      ++count;

  assert(count%2==0); // since |pos_to_neg| is supposed to be $\delta$-stable
  return count%4!=0;
}

  } // |namespace repr|
} // |namespace atlas|
