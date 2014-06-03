/*
  This is repr.cpp

  Copyright (C) 2009-2012 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "repr.h"

#include <map> // used in computing |reducibility_points|
#include <iostream>

#include "arithmetic.h"
#include "tits.h"

#include "kgb.h"	// various methods
#include "blocks.h"	// |dual_involution|
#include "standardrepk.h"// |KhatContext| methods
#include "subsystem.h"

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
{ return complexGroup().involution_of_Cartan(cn); }

StandardRepr
  Rep_context::sr
    (const standardrepk::StandardRepK& srk,
     const standardrepk::KhatContext& khc,
     const RatWeight& nu) const
{
  const KGBElt x= khc.kgb().lookup(khc.titsElt(srk));
  const InvolutionNbr i_x = kgb().inv_nr(x);
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);

  const Weight lambda2 = khc.lift(srk); // doubled coordinates
  const RatWeight lambda(lambda2,2);
  const RatWeight diff = lambda - nu;
  const RatWeight theta_diff(theta*diff.numerator(),
			     diff.denominator()); // theta(lambda-nu)
  const Weight lambda_rho = (lambda2-khc.rootDatum().twoRho())/2;
  return StandardRepr(x,i_tab.pack(i_x,lambda_rho),
		      ((lambda+nu+theta_diff)/=2).normalize());
}

StandardRepr Rep_context::sr_gamma
  (KGBElt x, const Weight& lambda_rho, const RatWeight& gamma) const
{
  const InvolutionTable& i_tab = complexGroup().involution_table();
  return StandardRepr(x, i_tab.pack(kgb().inv_nr(x),lambda_rho), gamma);
}

RatWeight Rep_context::gamma
  (KGBElt x, const Weight& lambda_rho, const RatWeight& nu) const
{
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const RatWeight lambda(lambda_rho*2+rootDatum().twoRho(),2);
  const RatWeight diff = lambda - nu;
  const RatWeight theta_diff(i_tab.matrix(kgb().inv_nr(x))*diff.numerator(),
			     diff.denominator()); // theta(lambda-nu)
  return ((lambda+nu+theta_diff)/=2).normalize();
}

StandardRepr Rep_context::sr(const param_block& b, BlockElt i) const
{
  assert(i<b.size());
  return sr(b.parent_x(i),b.lambda_rho(i),b.gamma());
}

Weight Rep_context::lambda_rho(const StandardRepr& z) const
{
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);

  const RatWeight gamma_rho = z.gamma() - RatWeight(rootDatum().twoRho(),2);
  Ratvec_Numer_t im_part2 = gamma_rho.numerator()+theta*gamma_rho.numerator();
  im_part2 /= gamma_rho.denominator(); // exact: $(1+\theta)(\lambda-\rho)$
  Weight i2(im_part2.begin(),im_part2.end()); // convert to |Weight|
  return (i2 + i_tab.unpack(i_x,z.y()))/2; // division exact again
}

// return $\lambda \in \rho+X^*$ as half-integer rational vector
RatWeight Rep_context::lambda(const StandardRepr& z) const
{
  const Weight num = lambda_rho(z) * 2 + rootDatum().twoRho();
  return RatWeight(num,2).normalize();
}

RatWeight Rep_context::nu(const StandardRepr& z) const
{
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const WeightInvolution& theta = complexGroup().involution_table().matrix(i_x);
  const Ratvec_Numer_t num = z.gamma().numerator()-theta*z.gamma().numerator();
  return RatWeight(num,2*z.gamma().denominator()).normalize();
}

// |z| standard means (weakly) dominant on the (simple-)imaginary roots
bool Rep_context::is_standard(const StandardRepr& z, RootNbr& witness) const
{
  const RootDatum& rd = rootDatum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = complexGroup().involution_table();
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
  const InvolutionTable& i_tab = complexGroup().involution_table();
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
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const RootNbrSet pos_real = i_tab.real_roots(i_x) & rd.posRootSet();
  const Weight test_wt = i_tab.unpack(i_x,z.y()) // $(1-\theta)(\lambda-\rho)$
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
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const RootNbrSet real = complexGroup().involution_table().real_roots(i_x);

  assert(real.isMember(alpha)); // only real roots should be tested

  const Weight& av = rootDatum().coroot(alpha);
  const int numer = av.dot(z.gamma().numerator());
  const int denom = z.gamma().denominator();
  assert(numer%denom!=0); // and the real root alpha should be non-integral

  const Weight test_wt = i_tab.unpack(i_x,z.y()) +rd.twoRho() -rd.twoRho(real);
  const int eps = av.dot(test_wt)%4==0 ? 0 : denom;

  return arithmetic::remainder(numer+eps,2*denom)< (unsigned)denom;
}

unsigned int Rep_context::orientation_number(const StandardRepr& z) const
{
  const RootDatum& rd = rootDatum();
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const RootNbrSet real = i_tab.real_roots(i_x);
  const Permutation& root_inv = i_tab.root_involution(i_x);
  const Ratvec_Numer_t& numer = z.gamma().numerator();
  const arithmetic::Numer_t denom = z.gamma().denominator();
  const Weight test_wt = i_tab.unpack(i_x,z.y()) +rd.twoRho() -rd.twoRho(real);

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
  const InvolutionTable& i_tab = complexGroup().involution_table();

  // the following are non-|const|, and modified in the loop below
  Weight lr = lambda_rho(z);
  KGBElt& x = z.x_part;
  Ratvec_Numer_t& numer = z.infinitesimal_char.numerator();
  InvolutionNbr i_x = kgb().inv_nr(x);

  { weyl::Generator s;
    do
    {
      for (s=0; s<rd.semisimpleRank(); ++s)
      {
	int v=rd.simpleCoroot(s).dot(numer);
        if (v<0 or (v==0 and kgb().isComplexDescent(s,x)))
        {
	  const RootNbr alpha = rd.simpleRootNbr(s);
	  if (i_tab.imaginary_roots(i_x).isMember(alpha))
	    throw std::runtime_error("Non standard parameter in make_dominant");
          rd.simpleReflect(numer,s);
          rd.simpleReflect(lr,s);
	  if (not i_tab.real_roots(i_x).isMember(alpha)) // if |alpha| is real
	    lr -= rd.simpleRoot(s); // then $\rho_r$ cancels $\rho$
          x = kgb().cross(s,x);
	  i_x = kgb().inv_nr(x); // keep up with changing involution
          break;
        }
      }
    }
    while (s<rd.semisimpleRank()); // wait until inner loop runs to completion
  }
  z.y_bits=i_tab.pack(i_x,lr);
} // |make_dominant|

RationalList Rep_context::reducibility_points(const StandardRepr& z) const
{
  const RootDatum& rd = rootDatum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const Permutation& theta = i_tab.root_involution(i_x);

  const RatWeight& gamma = z.gamma();
  const Ratvec_Numer_t& numer = gamma.numerator();
  const arithmetic::Numer_t d = gamma.denominator();
  const Weight lam_rho = lambda_rho(z);

  const RootNbrSet pos_real = i_tab.real_roots(i_x) & rd.posRootSet();
  const Weight two_rho_real = rd.twoRho(pos_real);

  // we shall associate to a first number a strict lower bound for some $k$
  // if first number is $num>0$ we shall later form fractions $(d/num)*k$
  typedef std::map<long,long> table;
  table odds,evens; // name indicates the parity that $k$ will have

  for (RootNbrSet::iterator it=pos_real.begin(); it(); ++it)
  {
    arithmetic::Numer_t num =
      rd.coroot(*it).dot(numer); // numerator of $\<\nu,\alpha^v>$
    if (num!=0)
    {
      long lam_alpha = lam_rho.dot(rd.coroot(*it))+rd.colevel(*it);
      bool do_odd = (lam_alpha+two_rho_real.dot(rd.coroot(*it))/2)%2 ==0;
      (do_odd ? &odds : &evens)->insert(std::make_pair(abs(num),0));
    }
  }

  RootNbrSet pos_complex = i_tab.complex_roots(i_x) & rd.posRootSet();
  for (RootNbrSet::iterator it=pos_complex.begin(); it(); ++it)
  {
    RootNbr alpha=*it, beta=theta[alpha];
    arithmetic::Numer_t vala = rd.coroot(alpha).dot(numer);
    arithmetic::Numer_t valb = rd.coroot(beta).dot(numer);
    arithmetic::Numer_t num = vala - valb; // numerator of $2\<\nu,\alpha^v>$
    if (num!=0)
    {
      assert((vala+valb)%d==0); // since |\<\gamma,a+b>=\<\lambda,a+b>|
      long lwb =abs(vala+valb)/d;
      std::pair<table::iterator,bool> trial =
	(lwb%2==0 ? &evens : &odds)->insert(std::make_pair(abs(num),lwb));
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
  RatWeight t =  src.y().log_pi(false);
  // InvolutionNbr i_x = kgb().inv_nr(z.x());
  // no need to do |complexGroup().involution_table().real_unique(i_x,t)|

  RatWeight lr =  (infin_char - t - RatWeight(rd.twoRho(),2)).normalize();
  assert(lr.denominator()==1);
  return StandardRepr
    (sr_gamma(src.x(),
	      Weight(lr.numerator().begin(),lr.numerator().end()),
	      infin_char));
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
  // no need to do |complexGroup().involution_table().real_unique(i_x,t)|

  RatWeight lr =  (infin_char - t - RatWeight(rd.twoRho(),2)).normalize();
  assert(lr.denominator()==1);
  return StandardRepr
    (sr_gamma(src.x(),
	      Weight(lr.numerator().begin(),lr.numerator().end()),
	      infin_char));
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
  // no need to do |complexGroup().involution_table().real_unique(i_x,t)|

  RatWeight lr =  (infin_char - t - RatWeight(rd.twoRho(),2)).normalize();
  assert(lr.denominator()==1);
  return StandardRepr
    (sr_gamma(src.x(),
	      Weight(lr.numerator().begin(),lr.numerator().end()),
	      infin_char));
}

StandardRepr Rep_context::twist(StandardRepr z) const
{
  make_dominant(z);
  const RatWeight infin_char=z.gamma(); // now get the infinitesimal character
  const RootDatum& rd = rootDatum();
  const SubSystem& subsys = SubSystem::integral(rd,infin_char);
  blocks::nblock_help aux(realGroup(),subsys);
  blocks::nblock_elt src(z.x(),y_values::exp_pi(infin_char-lambda(z)));
  aux.twist(src);
  RatWeight lr =
    (infin_char - src.y().log_pi(false) - RatWeight(rd.twoRho(),2)).normalize();
  assert(lr.denominator()==1);
  return StandardRepr
    (sr_gamma(src.x(),
	      Weight(lr.numerator().begin(),lr.numerator().end()),
	      infin_char));
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
  const InvolutionTable& i_tab = complexGroup().involution_table();

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
	  (s,rd.simpleCoroot(s).dot(i_tab.unpack(i_x,z.y()))%4!=0);
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

  assert(new_survivors.size()>0); // at least top element should be new
  assert(hash.size()==old_size+new_survivors.size()); // only new surv. added

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
  non_integral_block block(*this,z); // compute partial block
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
    non_integral_block block(*this,z);
    BlockEltList survivors;
    add_block(block,survivors);

    hash_index=hash.find(z);
    assert(hash_index!=hash.empty);
  }

  return KL_list[hash_index];
}

SR_poly Rep_table::deformation_terms (param_block& block,BlockElt entry_elem)
{
  SR_poly result(repr_less());
  if (not block.survives(entry_elem) or block.length(entry_elem)==0)
    return result; // easy cases, null result

  BlockEltList survivors;
  if (hash.find(sr(block,entry_elem))==hash.empty) // previously unknown
    add_block(block,survivors);

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

  assert(hash.find(sr(block,entry_elem))!=hash.empty); // should be known now

  // map indices of |survivors| to corresponding number in |hash|
  std::vector<unsigned long> remap(survivors.size());
  for (unsigned long i=0; i<survivors.size(); ++i)
  {
    unsigned long h=hash.find(sr(block,survivors[i]));
    assert(h!=hash.empty);
    remap[i]=h;
  }

  SR_poly Q(sr(block,entry_elem),repr_less()); // remainder, init (1,entry_elem)
  std::vector<Split_integer> acc(survivors.size(),Split_integer(0));

  for (unsigned long i=survivors.size(); i-->0; ) // decreasing essential here
  {
    StandardRepr p_y=sr(block,survivors[i]);
    Split_integer c_y = Q[p_y];
    const SR_poly& KL_y = KL_list[remap[i]];
    Q.add_multiple(KL_y,-c_y);
    assert(Q[p_y]==Split_integer(0)); // check relation of being inverse

    c_y.times_1_s(); // deformation terms are all multiplied by $1-s$
    acc[i]=c_y; // store coefficient at index of survivor
  }
  assert(Q.empty()); // since all terms in |KL_y| should be at most $y$

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
    non_integral_block b(*this,zi);
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


void Rep_table::add_block(ext_block::extended_block& block,
			  param_block& parent) // must be actual parent |block|
{
  unsigned long old_size = hash.size();
  BlockEltList survivors;
  add_block(parent,survivors);

  BlockEltList new_survivors;

  // fill the |hash| table for new surviving parameters in this block
  for (BlockEltList::const_iterator
	   it=survivors.begin(); it!=survivors.end(); ++it)
    if (hash.match(sr(parent,*it))>=old_size and
	parent.Hermitian_dual(*it)==*it)
      new_survivors.push_back(block.element(*it));

  // extend space in twisted tables
  twisted_KLV_list.resize(hash.size(),SR_poly(repr_less())); // init empty
  twisted_def_formula.resize(hash.size(),SR_poly(repr_less()));

  // compute cumulated KL polynomimals $P_{x,y}$ with $x\leq y$ survivors

  // start with computing KL polynomials for the entire block
  std::vector<ext_kl::Pol> pool;
  ext_kl::KL_table twisted_KLV(block,pool);
  twisted_KLV.fill_columns();

  /* get $P(x,y)$ for |x<=y| with |y| among new |survivors|, and contribute
   parameters from |block.survivors_below(x)| with coefficient $P(x,y)[q:=s]$
   to the |SR_poly| at |KL_list[old_size+i], where |y=new_survivors[i]| */
  BlockEltList::const_iterator y_start=new_survivors.begin();

  for (BlockElt x=0; x<=new_survivors.back(); ++x) // elements of twisted block
  {
    BlockElt z = block.z(x); // number of the element for |parent|
    BlockEltList xs=parent.survivors_below(z);
    if (xs.empty())
      continue; // no point doing work for |x|'s that don't contribute anywhere

    const unsigned int parity = block.length(x)%2;

    if (*y_start<x)
      ++y_start; // advance so |y| only runs over values with |x<=y|
    assert(y_start!=new_survivors.end() and *y_start>=x);

    for (BlockEltList::const_iterator it=y_start; it!=new_survivors.end(); ++it)
    {
      const BlockElt y = *it; // element of |new_survivors| and |x<=y|
      const ext_kl::Pol& pol = twisted_KLV.P(x,y); // twisted KLV polynomial
      Split_integer eval(0);
      for (polynomials::Degree d=pol.size(); d-->0; )
	eval = eval.times_s()+Split_integer(static_cast<int>(pol[d]));
      if (eval!=Split_integer(0))
      {
	unsigned long y_index = hash.find(sr(parent,block.z(y)));
	assert (y_index!=hash.empty);
	SR_poly& dest = KL_list[y_index];
	if (lengths[y_index]%2!=parity)
	  eval.negate(); // incorporate sign for length difference
	for (unsigned int i=0; i<xs.size(); ++i)
	  dest.add_term(sr(parent,xs[i]),eval);
      }
    } // |for(it)|
  } // |for(x)|
} // |Rep_table::add_block| (extended block)

// compute and return sum of KL polynomials at $s$ for final parameter |z|
SR_poly Rep_table::twisted_KL_column_at_s(StandardRepr z)
  // |z| must be twist-fixed, nonzero and final
{
  { RootNbr witness;
    if (is_zero(z,witness) or not is_final(z,witness))
      throw std::runtime_error("Representation zero or not final");
  }
  make_dominant(z); // so that |z| it will appear at the top of its own block
  unsigned long hash_index=hash.find(z);
  if (hash_index==hash.empty) // previously unknown parameter
  {
    BlockElt entry; // dummy needed to ensure full block is generated
    non_integral_block block(*this,z,entry); // which this constructor does
    ext_block::extended_block eblock(block,twistedWeylGroup());

    add_block(eblock,block);

    hash_index=hash.find(z);
    assert(hash_index!=hash.empty);
  }

  return twisted_KLV_list[hash_index];
}


  } // |namespace repr|
} // |namespace atlas|
