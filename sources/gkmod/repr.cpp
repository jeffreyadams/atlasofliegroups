/*
  This is repr.cpp

  Copyright (C) 2009-2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "repr.h"

#include <memory> // for |std::unique_ptr|
#include <map> // used in computing |reducibility_points|
#include <algorithm> // for |make_heap|
#include <iostream> // for progress reports and easier debugging
#include "error.h"

#include "arithmetic.h"
#include "matreduc.h"

#include "tits.h"
#include "kgb.h"	// various methods
#include "blocks.h"	// the |blocks::common_block| class, |dual_involution|
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

StandardReprMod StandardReprMod::mod_reduce
  (const Rep_context& rc, const StandardRepr& sr)
{
  KGBElt x = sr.x();
  auto rho = rootdata::rho(rc.root_datum());
  auto gam_lam=sr.gamma()-rho-rc.lambda_rho(sr);
  rc.involution_table().real_unique(rc.kgb().inv_nr(x),gam_lam);
  return StandardReprMod(x,rho+gam_lam);
}

StandardReprMod StandardReprMod::build
  (const Rep_context& rc, KGBElt x, RatWeight gam_lam)
{
  rc.involution_table().real_unique(rc.kgb().inv_nr(x),gam_lam);
  return StandardReprMod(x,rho(rc.root_datum())+gam_lam); // ctor normalises
}

size_t StandardReprMod::hashCode(size_t modulus) const
{ size_t hash = x_part + 47*rgl.denominator();
  for (auto entry : rgl.numerator())
    hash= 11*hash+entry;
  return hash &(modulus-1);
}

Reduced_param::Reduced_param
  (const Rep_context& rc, const StandardReprMod& srm)
    : x(srm.x()), evs_reduced() // |int_sys_nr| is set by |integral_eval| below
{
  InnerClass& ic = rc.inner_class();
  const KGB& kgb = rc.kgb();
  const auto& glr = srm.gamma_rep(); // $\gamma-\lambda+\rho$
  auto eval = ic.integral_eval(glr,int_sys_nr) * glr.numerator();
  for (auto& entry : eval)
  {
    assert(entry%glr.denominator()==0);
    entry /= glr.denominator();
  }
  const auto& codec = ic.int_item(int_sys_nr).data(ic,int_sys_nr,kgb.inv_nr(x));
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

  int_Matrix theta1 = kgb().involution_matrix(x)+1;
  Weight t1_gamma (gamma.numerator().begin(), gamma.numerator().end());
  // the division in the next computation may throw when |gamma| is bad for |x|
  t1_gamma = theta1*t1_gamma/static_cast<int>(gamma.denominator());
#ifndef NDEBUG // check that constructor below builds a valid |StandardRepr|
  Weight image = // $(\theta+1)(\gamma-\rho)$
    t1_gamma-(theta1*root_datum().twoRho()/2);
  matreduc::find_solution(theta1,image); // a solution must exist
#endif
  const InvolutionTable& i_tab = involution_table();
  return StandardRepr(x, i_tab.y_pack(kgb().inv_nr(x),lambda_rho), gamma,
		      height(t1_gamma));
}

StandardRepr Rep_context::sr
  (const StandardReprMod& srm, const RatWeight& diff, const RatWeight& gamma)
  const
{
  const RatWeight gamma_lambda_rho = srm.gamma_rep()+diff;
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
{ return involution_table().matrix(kgb().inv_nr(z.x())); }

Weight Rep_context::lambda_rho(const StandardRepr& z) const
{
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);


  // recover $\lambda-\rho$ from doubled projections on eigenspaces $\theta$
  const RatWeight gam_rho = z.gamma() - rho(root_datum());
  auto im_part2 = // $(1+\theta)(\lambda-\rho)$ found as |(1+theta)*(gam_rho)|
    (gam_rho.numerator()+theta*gam_rho.numerator()) / gam_rho.denominator();
  // both divisions, above and below, will be exact
  return (Weight(im_part2.begin(),im_part2.end()) +i_tab.y_lift(i_x,z.y())) / 2;
}

// compute $\gamma-\lambda$ from $(1-\theta)(\gamma-\lambda)=2(\gamma-\lambda)$
RatWeight Rep_context::gamma_lambda
  (InvolutionNbr i_x,  const TorusPart& y_bits, const RatWeight& gamma) const
{
  const InvolutionTable& i_tab = involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);

  const RatWeight gamma_rho = gamma - rho(root_datum());
  return (gamma_rho-theta*gamma_rho - i_tab.y_lift(i_x,y_bits))/2LL;
}

// compute $\gamma-\lambda-\rho$ from same information
RatWeight Rep_context::gamma_lambda_rho (const StandardRepr& sr) const
{
  const InvolutionTable& i_tab = involution_table();
  InvolutionNbr i_x = kgb().inv_nr(sr.x());
  const WeightInvolution& theta = i_tab.matrix(i_x);

  const Weight& rho2 = root_datum().twoRho();
  RatWeight result = sr.gamma()*2LL;
  (result -= theta*(result-rho2)) += rho2;

  return (std::move(result) - i_tab.y_lift(i_x,sr.y())*2)/4LL;
}

Weight Rep_context::theta_1_preimage
  (const RatWeight& offset, const subsystem::integral_datum_item::codec& codec)
  const
{
  auto& ic = inner_class();
  RatWeight eval =
    (ic.int_item(codec.int_sys_nr).coroots_matrix()*offset).normalize();
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

  return ic.involution_table().theta_1_image_basis(codec.inv) *
    (codec.out * eval_v);
}

RatWeight Rep_context::offset
  (const StandardReprMod& srm0, const StandardReprMod& srm1) const
{
  const auto& gam = srm0.gamma_rep(); // will also define integral system
  RatWeight result = gam - srm1.gamma_rep();
  auto& ic = inner_class();
  InvolutionNbr inv = kgb().inv_nr(srm0.x());
  unsigned int int_sys_nr;
  auto codec = ic.integrality_codec(gam,inv,int_sys_nr);
  result -= theta_1_preimage(result,codec);
  assert((ic.integral_eval(gam,int_sys_nr)*result.numerator()).isZero());
  return result;
}

RatWeight Rep_context::gamma_lambda(const StandardReprMod& z) const
{ return z.gamma_lambda(rho(root_datum())); }

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

// |z| standard means (weakly) dominant on the (simple-)imaginary roots
bool Rep_context::is_standard(const StandardRepr& z, RootNbr& witness) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = involution_table();
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
  const InvolutionTable& i_tab = involution_table();
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
  const InvolutionTable& i_tab = involution_table();
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
  const InvolutionTable& i_tab = involution_table();
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
  const InvolutionTable& i_tab = involution_table();
  const RootNbrSet real = i_tab.real_roots(i_x);

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
  z.y_bits = involution_table().y_pack(kgb().inv_nr(x),lr);
} // |make_dominant|

void Rep_context::singular_cross (weyl::Generator s,StandardRepr& z) const
{
  assert(root_datum().simpleCoroot(s).dot(z.gamma().numerator())==0);
  Weight lr = lambda_rho(z); auto& x=z.x_part;
  root_datum().simple_reflect
    (s, lr, kgb().status(s,x)==gradings::Status::Real ? 0 : 1);
  x = kgb().cross(s,x);
  z.y_bits = // reinsert $y$ bits component
    involution_table().y_pack(kgb().inv_nr(x),lr);
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

// make dominant, then descend though any singular complex descents
void Rep_context::deform_readjust(StandardRepr& z) const
{
  make_dominant(z); // typically may have gone negative only for complex coroots
  const RootDatum& rd = root_datum();

  RankFlags simple_singulars;
  { const auto& numer = z.infinitesimal_char.numerator();
    for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
      simple_singulars.set(s,rd.simpleCoroot(s).dot(numer)==0);
  }

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
  z.y_bits = involution_table().y_pack(kgb().inv_nr(x),lr);
} // |deform_readjust|

// this also ensures a chosen singular-complex minumum when there are multiple
// but that only arises when singular-real descents exist (not so in deformation)
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
  z.y_bits = involution_table().y_pack(kgb().inv_nr(x),lr);
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
} // |is_twist_fixed|

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
} // |Rep_context::equivalent|

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
  const InvolutionTable& i_tab = involution_table();
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
  i_x = kgb().inv_nr( x = kgb().cross(rd.reflectionWord(rt),x) );
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
  if (kgb.status(subsys.simple(s),conj_x)!=
      gradings::Status::ImaginaryNoncompact)
    throw error::Cayley_error();
  const auto new_x = kgb.cross(kgb.cayley(subsys.simple(s),conj_x),conj);
  RatWeight gamma_lambda = this->gamma_lambda(z);

  const RootNbrSet& upstairs_real_roots = i_tab.real_roots(kgb.inv_nr(new_x));
  RootNbrSet real_flip = upstairs_real_roots;
  real_flip ^= i_tab.real_roots(kgb.inv_nr(z.x())); // remove downstairs reals

  RootNbrSet pos_neg = pos_to_neg(rd,conj);
  pos_neg &= real_flip; // posroots that change real status and map to negative
  gamma_lambda += root_sum(rd,pos_neg); // correction of $\rho_r$'s

  // correct in case the parity condition fails for our raised |gamma_lambda|
  const Coweight& alpha_hat = rd.coroot(parent_s);
  const int rho_r_corr = // integer since alpha is among |upstairs_real_roots|
    alpha_hat.dot(rd.twoRho(upstairs_real_roots))/2;
  const int eval = gamma_lambda.dot(alpha_hat);
  if ((eval+rho_r_corr)%2==0) // parity condition says it should be 1
    gamma_lambda += RatWeight(rd.root(parent_s),2); // add half-alpha

  const Weight lambda_rho = gamma.integer_diff<int>(gamma_lambda+rho(rd));
  return sr_gamma(new_x,lambda_rho,gamma);
}

StandardRepr Rep_context::inv_Cayley(weyl::Generator s, StandardRepr z) const
{
  make_dominant(z);
  const RootDatum& rd = root_datum();
  const auto& i_tab = involution_table();
  const KGB& kgb = this->kgb();
  const RatWeight& gamma = z.gamma(); // now get the infinitesimal character
  const SubSystem& subsys = SubSystem::integral(rd,gamma);
  const auto parent_s = subsys.parent_nr_simple(s);

  const auto& real_roots = i_tab.real_roots(kgb.inv_nr(z.x()));
  if (not real_roots.isMember(parent_s))
    throw error::Cayley_error();
  RatWeight gamma_lambda = this->gamma_lambda(z);
  const Coweight& alpha_hat = rd.coroot(parent_s);
  const int eval = gamma_lambda.dot(alpha_hat);
  const int rho_r_corr = alpha_hat.dot(rd.twoRho(real_roots))/2;
  if ((eval+rho_r_corr)%2==0) // then |s| is real nonparity at |z|
    throw error::Cayley_error();

  const auto& conj = subsys.to_simple(s); // word in full system
  const KGBElt conj_x = kgb.cross(conj,z.x());
  const KGBElt new_x =
    kgb.cross(kgb.inverseCayley(subsys.simple(s),conj_x).first,conj);

  RootNbrSet pos_neg = pos_to_neg(root_datum(),conj);
  const RootNbrSet real_flip = real_roots^i_tab.real_roots(kgb.inv_nr(new_x));
  pos_neg &= real_flip; // posroots that change real status and map to negative
  gamma_lambda += root_sum(rd,pos_neg); // correction of $\rho_r$'s
  // now |gamma_lambda| is still in the $X^*$-coset of $\gamma-\rho$
  // it might not be in the $-1$ eigenspace for |new_x|, but |sr_gamma| projects

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
	  i_x = kgb().inv_nr( x = kgb().cross(rd.reflectionWord(alpha),x) );
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

SR_poly Rep_context::scale(const poly& P, const Rational& f) const
{
  poly result;
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
  poly result;
  for (auto it=P.begin(); it!=P.end(); ++it)
  { auto z=it->first; // take a copy for modification
    auto finals = finals_for(scale_0(z));
    for (const StandardRepr& final : finals)
      result.add_term(final,it->second);
  }
  return result;
}

sl_list<StandardRepr> Rep_context::finals_for(StandardRepr z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionTable& i_tab = involution_table();

  make_dominant(z); // ensures singular subsystem is generated by simple roots
  const RatWeight& gamma=z.gamma();

  RankFlags singular; // subset of simple roots that is singular at |gamma|
  for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
    singular.set(s,rd.simpleCoroot(s).dot(gamma.numerator())==0);
/* the simple roots flagged in |singular| coincide with the singular
   integrally-simple roots (simple roots of integral subsystem), so we can
   proceed without constructing the integral subsystem (even less the block) */

  sl_list<StandardRepr> result { z };
  auto rit=result.begin();
  do
  { const KGBElt x=rit->x();
    for (const auto s : singular)
      // as |break| from loop is not available within |switch|, use |goto| below
      switch (kgb().status(s,x))
      {
      case gradings::Status::ImaginaryCompact:
	result.erase(rit); // discard zero parameter
	goto repeat;
      case gradings::Status::Complex:
        if (kgb().isDescent(s,x))
	{ // replace |*rit| by its complex descent for |s|
	  singular_cross(s,*rit);
	  goto repeat; // reconsider all singular roots for the new parameter
	}
	break;
      case gradings::Status::Real:
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
	    goto repeat;
	  }
	}
	break;
      case gradings::Status::ImaginaryNoncompact:
        break; // like for complex ascent or real nonparity, do nothing here
      } // end of |switch| and of |for (s:singular)|
    ++rit; // no singular descents found: consolidate |*rit|
  repeat: {} // continue modifications at |*rit| after jump here
  }
  while (not result.at_end(rit));
  return result;
} // |Rep_context::finals_for|

SR_poly Rep_context::expand_final (StandardRepr z) const
{
  poly result;
  for (const auto& sr : finals_for(z))
    result += sr;
  return result;
} // |Rep_context::expand_final|


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

bool deformation_unit::operator!=(const deformation_unit& another) const
{
  if (sample.x()!=another.sample.x() or
      sample.height()!=another.sample.height() or
      sample.y().data()!=another.sample.y().data())
    return true; // easy tests for difference

  {
    const int_Matrix& theta = rc.kgb().involution_matrix(sample.x());
    const auto gamma_diff_num = (sample.gamma()-another.sample.gamma()
				 ).numerator();
    if (not (theta*gamma_diff_num+gamma_diff_num).isZero())
      return true; // difference in the free part of $\lambda$ spotted
  }

  const auto& rd = rc.root_datum();
  const auto& g0=sample.gamma();
  const auto& g1=another.sample.gamma();
  const auto& num0 = g0.numerator();
  const auto& num1 = g1.numerator();
  const int d0=g0.denominator(), d1=g1.denominator(); // convert to |int|
  for (auto it=rd.beginPosCoroot(); it!=rd.endPosCoroot(); ++it)
    if (arithmetic::divide(num0.dot(*it),d0) !=
	arithmetic::divide(num1.dot(*it),d1))
      return true; // different integer part of evaluation on poscoroot found

  return false; // if no differences detected, consider |another| as equivalent
}

size_t deformation_unit::hashCode(size_t modulus) const
{
  size_t hash = 7*sample.x() + 89*sample.y().data().to_ulong();
  const int_Matrix& theta = rc.kgb().involution_matrix(sample.x());
  const auto& num = sample.gamma().numerator();
  const auto free_lambda = (theta*num+num)/sample.gamma().denominator();
  for (int c : free_lambda) // take into account free part of $\lambda$
    hash = 21*(hash&(modulus-1))+c;

  const auto& rd = rc.root_datum();
  const int denom = sample.gamma().denominator(); // convert to |int|
  for (auto it=rd.beginPosCoroot(); it!=rd.endPosCoroot(); ++it)
    hash = 5*(hash&(modulus-1))+arithmetic::divide(num.dot(*it),denom);

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
  auto & block = lookup(sr,z); // construct partial block
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

  const auto rank = ctxt.id().semisimpleRank();
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
    (find_reduced_hash(init)==reduced_hash.empty);

  std::vector<StandardReprMod> pool;
  Mod_hash_tp hash(pool);
  Bruhat_generator gen(hash,ctxt); // object to help generating Bruhat interval
  gen.block_below(init); // generate Bruhat interval below |srm| into |pool|

  const size_t place_limit = place.size();
  sl_list<sub_triple> sub_blocks;
  for (const auto& elt : pool)
  {
    auto h = match_reduced_hash(elt); // fingerprint of parameter family
    if (h==place.size()) // block element has new reduced hash value
      place.emplace_back(bl_it(),-1); // create slot; both fields filled later
    else if (h<place_limit) // then a similar parameter was known
    { // record block pointer and index of |z| in block into |sub_blocks|
      common_block* sub = &*place[h].first;
      auto hit = [sub] (const sub_triple& tri)->bool { return tri.bp==sub; };
      if (std::none_of(sub_blocks.begin(),sub_blocks.end(),hit))
      {
	StandardReprMod base = sub->representative(place[h].second);
	sub_blocks.emplace_back(sub,h,offset(elt,base));
      }
    }
  }

  // reconstruct elements for all |sub_blocks|
  size_t limit = pool.size(); // limit of generated Bruhat interval
  const auto rho = rootdata::rho(root_datum());
  for (auto sub : sub_blocks)
    for (BlockElt z=0; z<sub.bp->size(); ++z)
    {
      auto new_gam_lam =
	sub.bp->representative(z).gamma_lambda(rho)+sub.shift;
      auto shifted =
	StandardReprMod::build(*this,sub.bp->x(z),std::move(new_gam_lam));
      hash.match(shifted); // if new, add |shifted| to |pool| beyong |limit|
    }

  sl_list<StandardReprMod> elements(pool.begin(),pool.end()); // working copy

  sl_list<blocks::common_block> temp; // must use temporary singleton
  auto& block = temp.emplace_back // construct block and get a reference
    (ctxt,elements,init.gamma_rep());

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
    auto& sub_block = *sub.bp;
    BlockEltList embed; embed.reserve(sub_block.size()); // translation array
    for (BlockElt z=0; z<sub_block.size(); ++z)
    {
      const auto& elt = sub_block.representative(z);
      auto new_gam_lam = elt.gamma_lambda(rho)+sub.shift;
      auto elt_shifted =
	StandardReprMod::build(*this,sub_block.x(z),std::move(new_gam_lam));
      assert(find_reduced_hash(elt)==find_reduced_hash(elt_shifted));
      const BlockElt z_rel = block.lookup(elt_shifted);
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
    const StandardReprMod srm =block.representative(z);
    auto h = find_reduced_hash(srm);
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
      auto zm = block.representative(z);
      auto seq = find_reduced_hash(zm);
      assert(seq<place.size()); // all elements in |block_list| must have |place|
      if (place[seq].first==next_pos) // could be false if |block| was swallowed
	place[seq].first=pos; // replace iterator that is about to be invalidated
    }
  }
  block_list.erase(pos);
} // |Rep_table::block_erase|

unsigned long Rep_table::add_block(const StandardReprMod& srm)
{
  BlockElt srm_in_block; // will hold position of |srm| within that block
  sl_list<blocks::common_block> temp; // must use temporary singleton
  auto& block = temp.emplace_back(*this,srm,srm_in_block); // build full block

  const auto rho = rootdata::rho(root_datum());
  const size_t place_limit = place.size();

  sl_list<sub_triple> sub_blocks;
  for (BlockElt z=0; z<block.size(); ++z)
  {
    auto elt = block.representative(z);
    auto seq = match_reduced_hash(elt);
    if (seq==place.size()) // block element has new reduced hash value
      place.emplace_back(bl_it(),z); // create slot; iterator filled later
    else if (seq<place_limit)
    {
      common_block* sub = &*place[seq].first;
      auto hit = [sub] (const sub_triple& tri)->bool { return tri.bp==sub; };
      if (std::none_of(sub_blocks.begin(),sub_blocks.end(),hit))
      {
	StandardReprMod base = sub->representative(place[seq].second);
	sub_blocks.emplace_back(sub,seq,offset(elt,base));
      }
    }
  }

  // swallow |embeddings|, and remove them from |block_list|
  for (const auto& sub : sub_blocks) // swallow sub-blocks
  {
    auto& sub_block = *sub.bp;
    BlockEltList embed; embed.reserve(sub_block.size()); // translation array
    for (BlockElt z=0; z<sub_block.size(); ++z)
    {
      auto new_gam_lam =
	sub_block.representative(z).gamma_lambda(rho)+sub.shift;
      auto shifted =
	StandardReprMod::build(*this,sub_block.x(z),std::move(new_gam_lam));
      const BlockElt z_rel = block.lookup(shifted);
      assert(z_rel!=UndefBlock);
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
    auto h = find_reduced_hash(block.representative(z));
    place[h].first = new_block_it;
    place[h].second = z;
  }
  return find_reduced_hash(srm);
}// |Rep_table::add_block|

blocks::common_block& Rep_table::lookup_full_block (StandardRepr& sr,BlockElt& z)
{
  make_dominant(sr); // without this we would not be in any valid block
  auto srm = StandardReprMod::mod_reduce(*this,sr); // modular |z|
  auto h = find_reduced_hash(srm); // look up modulo $X^*+integral^\perp$
  if (h==reduced_hash.empty or not place[h].first->is_full()) // then we must
    h=add_block(srm); // generate a new full block (possibly swalllow older ones)
  assert(h<place.size() and place[h].first->is_full());

  z = place[h].second;
  return *place[h].first;

} // |Rep_table::lookup_full_block|

blocks::common_block& Rep_table::lookup (StandardRepr& sr,BlockElt& which)
{
  normalise(sr); // gives a valid block, and smallest partial block
  auto srm = StandardReprMod::mod_reduce(*this,sr); // modular |z|
  assert(reduced_hash.size()==place.size()); // should be in sync at this point
  auto h = find_reduced_hash(srm); // look up modulo $X^*+integral^\perp$
  if (h!=reduced_hash.empty) // then we have found our family of blocks
  {
    assert(h<place.size());
    which = place[h].second;
    auto& block = *place[h].first;
    assert(block.representative(which).x()==srm.x()); // check minimum of sanity
    return block; // use block of related |StandardReprMod| as ours
  }
  common_context ctxt(*this,SubSystem::integral(root_datum(),sr.gamma()));
  BitMap subset;
  auto& block= add_block_below(ctxt,srm,&subset); // ensure block is known
  which = last(subset);
  assert(block.representative(which)==srm);
  return block;
} // |Rep_table::lookup|

// in the following type the second component is a mulitplicity so we are in fact
// dealing with a sparse reprensetion of polynomials with |BlockElt| exponents

typedef std::pair<BlockElt,int>  term;
typedef sl_list<term> pair_list;

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

sl_list<std::pair<StandardRepr,int> > Rep_table::deformation_terms
  ( blocks::common_block& block, const BlockElt y,
    const RatWeight& diff, const RatWeight& gamma)
{ assert(y<block.size()); // and |y| is final, see |assert| below

  sl_list<std::pair<StandardRepr,int> > result;
  if (block.length(y)==0)
    return result; // easy case, null result

  std::vector<pair_list> contrib = contributions(block,block.singular(gamma),y);
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

// compute and return sum of KL polynomials at $s$ for final parameter |sr|
SR_poly Rep_table::KL_column_at_s(StandardRepr sr) // |sr| must be final
{
  normalise(sr); // implies that |sr| it will appear at the top of its own block
  assert(is_final(sr));

  BlockElt z;
  auto& block = lookup(sr,z);
  RatWeight diff = offset(sr,block.representative(z));
  assert((involution_table().matrix(kgb().inv_nr(block.x(z)))*diff+diff).isZero());

  const auto& gamma=sr.gamma();
  std::vector<pair_list> contrib = contributions(block,block.singular(gamma),z);
  assert(contrib.size()==z+1 and contrib[z].front().first==z);

  const kl::KL_table& kl_tab =
    block.kl_tab(&KL_poly_hash,z+1); // fill silently up to |z|

  SR_poly result;
  auto z_length=block.length(z);
  for (BlockElt x=z+1; x-->0; )
  {
    const kl::KLPol& pol = kl_tab.KL_pol(x,z); // regular KL polynomial
    if (pol.isZero())
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

  BlockElt z;
  auto& block = lookup(sr,z);

  const kl::KL_table& kl_tab =
    block.kl_tab(&KL_poly_hash,z+1); // fill silently up to |z|

  simple_list<std::pair<BlockElt,kl::KLPol> > result;
  for (BlockElt x=z+1; x-->0; )
  {
    const kl::KLPol& pol = kl_tab.KL_pol(x,z); // regular KL polynomial
    if (not pol.isZero())
      result.emplace_front(x,pol);
  }

  return result;
} // |Rep_table::KL_column|


K_type_poly Rep_table::deformation(const StandardRepr& z)
// that |z| is dominant and final is a precondition assured in the recursion
// for more general |z|, do the preconditioning outside the recursion
{
  assert(is_final(z));
  RationalList rp=reducibility_points(z); // this is OK before |make_dominant|
  StandardRepr z_near = z;
  if (not rp.empty() and rp.back()!=Rational(1,1))
  { // then shrink wrap toward $\nu=0$
    scale(z_near,rp.back()); // snap to nearest reducibility point
    deform_readjust(z_near); // so that we may find a stored equivalent parameter
    assert(is_final(z_near));
  }

  deformation_unit zn(*this,std::move(z_near));
  { // look up if deformation formula for |z_near| is already known and stored
    unsigned long h=alcove_hash.find(zn);
    if (h!=alcove_hash.empty and pool[h].has_deformation_formula())
      return pool[h].def_formula();
  }

  StandardRepr z0 = z; scale_0(z0);
  K_type_poly result {std::less<K_type_nr>()};
  for (const auto& sr : finals_for(z0))
  {
    K_type_nr h = K_type_hash.match(K_type(*this,sr));
    result.add_term(h,Split_integer(1,0));
  }

  if (rp.empty()) // without deformation terms
    return std::move(result).flatten(); // don't even bother to store the result

  // otherwise compute the deformation terms at all reducibility points

  for (unsigned i=rp.size(); i-->0; )
  {
    auto zi = z; scale(zi,rp[i]);
    deform_readjust(zi); // necessary to ensure the following |assert| will hold
    assert(is_final(zi)); // ensures that |deformation_terms| won't refuse
    BlockElt new_z;
    auto& block = lookup(zi,new_z);
    RatWeight diff = offset(zi, block.representative(new_z));
    assert((involution_table().matrix(kgb().inv_nr(block.x(new_z)))*diff+diff)
	   .isZero());
    auto dt = deformation_terms(block,new_z,diff,zi.gamma());
    for (auto const& term : dt)
    { auto def = deformation(term.first);
      result.add_multiple(std::move(def), // recursion
			  Split_integer(term.second,-term.second)); // $(1-s)*c$
    }
  }

  const auto h = alcove_hash.match(zn); // now allocate a slot in |pool|
  return pool[h].set_deformation_formula(std::move(result).flatten());
} // |Rep_table::deformation|


// basic computation of twisted KL column sum, no tabulation of the result
SR_poly twisted_KL_sum
( ext_block::ext_block& eblock, BlockElt y, const blocks::common_block& parent,
  const RatWeight& diff,
  const RatWeight& gamma) // infinitesimal character, possibly singular
{
  // compute cumulated KL polynomimals $P_{x,y}$ with $x\leq y$ survivors

  // start with computing KL polynomials for the entire block
  std::vector<ext_kl::Pol> pool;
  ext_KL_hash_Table hash(pool,4);
  ext_kl::KL_table twisted_KLV(eblock,&hash);
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

  RankFlags singular_orbits; // flag singulars among orbits
  const auto& ipd = parent.integral_subsystem().pre_root_datum();
  for (weyl::Generator s=0; s<eblock.rank(); ++s)
    singular_orbits.set(s,gamma.dot(ipd.simple_coroot(eblock.orbit(s).s0))==0);

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
  blocks::common_block block(rc,zm,entry); // which this constructor does
  RatWeight diff = rc.offset(z, block.representative(entry));
  assert((rc.involution_table().matrix(rc.kgb().inv_nr(z.x()))*diff+diff)
	 .isZero());
  ext_block::ext_block eblock(block,delta,nullptr);

  return twisted_KL_sum(eblock,eblock.element(entry),block,diff,z.gamma());
} // |twisted_KL_column_at_s|

// look up or compute and return the alternating sum of twisted KL polynomials
// at the inner class involution for final parameter |z|, evaluated at $q=s$
SR_poly Rep_table::twisted_KL_column_at_s(StandardRepr sr)
  // |z| must be inner-class-twist-fixed, nonzero and final
{
  normalise(sr);
  assert(is_final(sr) and sr==inner_twisted(sr));
  BlockElt y0;
  auto& block = lookup(sr,y0);
  RatWeight diff = offset(sr,block.representative(y0));

  auto& eblock = block.extended_block(&poly_hash);

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
    if (pol.isZero())
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

K_type_poly Rep_table::twisted_deformation (StandardRepr z)
{
  assert(is_final(z));
  const auto& delta = inner_class().distinguished();
  bool flip_start=false; // whether a flip in descending to first point

  RationalList rp=reducibility_points(z);
  if (not rp.empty() and rp.back()!=Rational(1,1))
  { // then shrink wrap toward $\nu=0$
    const Rational f=rp.back();
    z = ext_block::scaled_extended_dominant(*this,z,delta,f,flip_start);
    for (auto& a : rp)
      a/=f; // rescale reducibility points to new parameter |z|
    assert(rp.back()==Rational(1,1)); // should make first reduction at |z|
    // here we continue, with |flip_start| recording whether we already flipped
  }

  deformation_unit zu(*this,z);
  { // if formula for |z| was previously stored, return it with |s^flip_start|
    const auto h=alcove_hash.find(zu);
    if (h!=alcove_hash.empty and pool[h].has_twisted_deformation_formula())
      return flip_start // if so we must multiply the stored value by $s$
	? K_type_poly().add_multiple
	(pool[h].twisted_def_formula(),Split_integer(0,1))
	: pool[h].twisted_def_formula();
  }

  K_type_poly result { std::less<K_type_nr>() };
  { // initialise |result| to fully deformed parameter expanded to finals
    bool flipped; // contrary to |flip_start| this affects value stored for |z|
    auto z0 = ext_block::scaled_extended_dominant
		(*this,z,delta,Rational(0,1),flipped); // deformation to $\nu=0$
    auto L = ext_block::extended_finalise(*this,z0,delta);
    for (const std::pair<StandardRepr,bool>& p : L)
    {
      auto h = K_type_hash.match(K_type(*this,p.first));
      result.add_term(h,p.second==flipped // flip means |times_s|
			? Split_integer(1,0) : Split_integer(0,1) );
    }
  }

  if (rp.empty())
    return std::move(result).flatten(); // return without storing in easy cases

  // compute the deformation terms at all reducibility points
  for (unsigned i=rp.size(); i-->0; )
  {
    Rational r=rp[i]; bool flipped;
    auto zi = ext_block::scaled_extended_dominant(*this,z,delta,r,flipped);
    auto L =
      ext_block::extended_finalise(*this,zi,delta); // rarely a long list

    for (std::pair<StandardRepr,bool>& p : L)
    {
      BlockElt new_z;
      auto& block = lookup(p.first,new_z);
      RatWeight diff = offset(p.first,block.representative(new_z));
      auto& eblock = block.extended_block(&poly_hash);

      RankFlags singular = block.singular(p.first.gamma());
      RankFlags singular_orbits; // flag singulars among orbits
      for (weyl::Generator s=0; s<eblock.rank(); ++s)
	singular_orbits.set(s,singular[eblock.orbit(s).s0]);

      auto terms = twisted_deformation_terms(block,eblock,new_z,
					     singular_orbits,diff,zi.gamma());
      const bool flip = flipped!=p.second;
      for (auto const& term : terms)
	result.add_multiple(twisted_deformation(term.first), // recursion
			    flip ? Split_integer(-term.second,term.second)
				 : Split_integer(term.second,-term.second));
    }
  }

  const auto h = alcove_hash.match(zu);  // now find or allocate a slot in |pool|
  const auto& res =
    pool[h].set_twisted_deformation_formula(std::move(result).flatten());

  return flip_start // if so we must multiply the stored value by $s$
    ? K_type_poly().add_multiple(res,Split_integer(0,1)) : res;

} // |Rep_table::twisted_deformation (StandardRepr z)|


//			|common_context| methods

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
  integr_datum.simple_reflect(s,gamma_lambda.numerator()); // integrally simple
  return repr::StandardReprMod::build(rc(),new_x,gamma_lambda);
}

bool common_context::is_parity
    (weyl::Generator s, const StandardReprMod& z) const
{
  const auto& full_datum = full_root_datum();
  const auto& i_tab = involution_table();
  const auto& real_roots = i_tab.real_roots(kgb().inv_nr(z.x()));
  assert(real_roots.isMember(sub.parent_nr_simple(s)));
  const Coweight& alpha_hat = integr_datum.simpleCoroot(s);
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
  const Coweight& alpha_hat = integr_datum.simpleCoroot(s);
  const int rho_r_corr = // integer since alpha is among |upstairs_real_roots|
    alpha_hat.dot(full_datum.twoRho(upstairs_real_roots))/2;
  const int eval = gamma_lambda.dot(alpha_hat);
  if ((eval+rho_r_corr)%2==0) // parity condition says it should be 1
    gamma_lambda += RatWeight(integr_datum.root(s),2); // add half-alpha

  return repr::StandardReprMod::build(rc(),new_x,gamma_lambda);
}


Weight common_context::to_simple_shift
  (InvolutionNbr theta, InvolutionNbr theta_p, RootNbrSet S) const
{ const InvolutionTable& i_tab = involution_table();
  S &= (i_tab.real_roots(theta) ^i_tab.real_roots(theta_p));
  return root_sum(full_root_datum(),S);
}


//			|Ext_common_context| methods

Ext_common_context::Ext_common_context
  (const Rep_context& rc, const WeightInvolution& delta, const SubSystem& sub)
    : repr::common_context(rc,sub)
    , d_delta(delta)
    , pi_delta(rc.root_datum().rootPermutation(delta))
    , delta_fixed_roots(fixed_points(pi_delta))
    , twist()
    , l_shifts(id().semisimpleRank())
{
  const RootDatum& rd = rc.root_datum();
  for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
    twist[s] = rd.simpleRootIndex(delta_of(rd.simpleRootNbr(s)));

  // the reflections for |E.l| pivot around |g_rho_check()|
  const RatCoweight minus_g_rho_check = -rc.g_rho_check();
  for (unsigned i=0; i<l_shifts.size(); ++i)
    l_shifts[i] = minus_g_rho_check.dot(id().simpleRoot(i));
} // |Ext_common_context::Ext_common_context|


bool Ext_common_context::is_very_complex
  (InvolutionNbr theta, RootNbr alpha) const
{ const auto& i_tab = involution_table();
  const auto& rd = full_root_datum();
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
bool Ext_common_context::shift_flip
  (InvolutionNbr theta, InvolutionNbr theta_p, RootNbrSet S) const
{ S.andnot(delta_fixed()); // $\delta$-fixed roots won't contribute

  unsigned count=0; // will count 2-element |delta|-orbit elements
  for (auto it=S.begin(); it(); ++it)
    if (is_very_complex(theta,*it) != is_very_complex(theta_p,*it) and
	not full_root_datum().sumIsRoot(*it,delta_of(*it)))
      ++count;

  assert(count%2==0); // since |pos_to_neg| is supposed to be $\delta$-stable
  return count%4!=0;
}

// |Ext_rep_context| methods

Ext_rep_context::Ext_rep_context
  (const repr::Rep_context& rc, const WeightInvolution& delta)
: Rep_context(rc), d_delta(delta) {}

Ext_rep_context::Ext_rep_context (const repr::Rep_context& rc)
: Rep_context(rc), d_delta(rc.inner_class().distinguished()) {}

// |common_context| methods

common_context::common_context (const Rep_context& rc, const SubSystem& sub)
: rep_con(rc)
, integr_datum(sub.pre_root_datum())
, sub(sub)
{} // |common_context::common_context|


  } // |namespace repr|
} // |namespace atlas|
