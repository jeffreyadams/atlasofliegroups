/*
  This is repr.cpp

  Copyright (C) 2009-2012 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "repr.h"

#include <map> // used in computing |reducibility_points|

#include "arithmetic.h"
#include "tits.h"

#include "kgb.h"	// various methods
#include "blocks.h"	// |dual_involution|
#include "standardrepk.h"// |KhatContext| methods

#include "kl.h"

namespace atlas {
  namespace repr {

bool StandardRepr::operator== (const StandardRepr& z) const
{ return x_part==z.x_part and y_bits==z.y_bits
  and infinitesimal_char==z.infinitesimal_char;
}

size_t StandardRepr::hashCode(size_t modulus) const
{ size_t hash=
    x_part + 375*y_bits.data().to_ulong()+83*infinitesimal_char.denominator();
  const Weight& num=infinitesimal_char.numerator();
  for (unsigned i=0; i<num.size(); ++i)
    hash= 11*(hash&(modulus-1))+num[i];
  return hash &(modulus-1);
}

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
  const TitsElt a = khc.titsElt(srk); // was reduced during construction |srk|
  const KGBElt x= khc.kgb().lookup(a,titsGroup());
  const InvolutionNbr i_x = kgb().inv_nr(x);
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);

  const Weight lambda2 = khc.lift(srk); // doubled coordinates
  const RatWeight lambda(lambda2,2);
  const RatWeight diff = lambda - nu;
  const RatWeight theta_diff(theta*diff.numerator(),
			     diff.denominator()); // theta(lambda-nu)
  const Weight lambda_rho = (lambda2-khc.rootDatum().twoRho())/=2;
  return StandardRepr(x,i_tab.pack(i_x,lambda_rho),
		      ((lambda+nu+theta_diff)/=2).normalize());
}

StandardRepr Rep_context::sr
  (KGBElt x, const Weight lambda_rho, const RatWeight& nu) const
{
  const InvolutionNbr i_x = kgb().inv_nr(x);
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);
  const RatWeight lambda(lambda_rho*2+rootDatum().twoRho(),2);
  const RatWeight diff = lambda - nu;
  const RatWeight theta_diff(theta*diff.numerator(),
			     diff.denominator()); // theta(lambda-nu)
  return StandardRepr(x,i_tab.pack(i_x,lambda_rho),
		      ((lambda+nu+theta_diff)/=2).normalize());
}

Weight Rep_context::lambda_rho(const StandardRepr& z) const
{
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const WeightInvolution& theta = i_tab.matrix(i_x);

  const RatWeight gamma_rho = z.gamma() - RatWeight(rootDatum().twoRho(),2);
  Weight im_part2 = gamma_rho.numerator()+theta*gamma_rho.numerator();
  im_part2 /= gamma_rho.denominator(); // exact: $(1+\theta)(\lambda-\rho)$
  return (im_part2 + i_tab.unpack(i_x,z.y()))/=2; // division exact again
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
  const Weight num = z.gamma().numerator()-theta*z.gamma().numerator();
  return RatWeight(num,2*z.gamma().denominator()).normalize();
}

// |z| standard means (weakly) dominant on the (simple-)imaginary roots
bool Rep_context::is_standard(const StandardRepr& z, RootNbr& witness) const
{
  const RootDatum& rd = rootDatum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const Weight& numer = z.gamma().numerator();

  for (unsigned i=0; i<i_tab.imaginary_rank(i_x); ++i)
  {
    const RootNbr alpha = i_tab.imaginary_basis(i_x,i);
    if (numer.dot(rd.coroot(alpha))<0)
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
  const Weight& numer = z.gamma().numerator();

  for (unsigned i=0; i<i_tab.imaginary_rank(i_x); ++i)
  {
    const RootNbr alpha = i_tab.imaginary_basis(i_x,i);
    if (numer.dot(rd.coroot(alpha))==0 and // simple-imaginary, singular
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

void Rep_context::make_dominant(StandardRepr& z) const
{
  const RootDatum& rd = rootDatum();
  const InvolutionTable& i_tab = complexGroup().involution_table();

  // the following are non-|const|, and modified in the loop below
  Weight lr = lambda_rho(z);
  KGBElt& x = z.x_part;
  Weight& numer = z.infinitesimal_char.numerator();
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
}

RationalList Rep_context::reducibility_points(const StandardRepr& z) const
{
  const RootDatum& rd = rootDatum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = complexGroup().involution_table();
  const Permutation& theta = i_tab.root_involution(i_x);

  const RatWeight& gamma = z.gamma();
  const Weight& numer = gamma.numerator();
  const long d = gamma.denominator();
  const Weight lam_rho = lambda_rho(z);

  const RootNbrSet pos_real = i_tab.real_roots(i_x) & rd.posRootSet();
  const Weight two_rho_real = rd.twoRho(pos_real);

  // we shall associate to a first number a strict lower bound for some $k$
  // if first number is $num>0$ we shall later form fractions $(d/num)*k$
  typedef std::map<long,long> table;
  table odds,evens; // name indicates the parity that $k$ will have

  for (RootNbrSet::iterator it=pos_real.begin(); it(); ++it)
  {
    long num = numer.dot(rd.coroot(*it)); // numerator of $\<\nu,\alpha^v>$
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
    long vala = numer.dot(rd.coroot(alpha));
    long valb = numer.dot(rd.coroot(beta));
    long num = vala - valb;   // numerator of $2\<\nu,\alpha^v>$
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
  const Weight& numer = z.gamma().numerator();
  const int denom = z.gamma().denominator();
  const Weight test_wt = i_tab.unpack(i_x,z.y()) +rd.twoRho() -rd.twoRho(real);

  unsigned count = 0;

  for (unsigned i=0; i<rd.numPosRoots(); ++i)
  {
    const RootNbr alpha = rd.numPosRoots()+i;
    const Weight& av = rootDatum().coroot(alpha);
    const int num = av.dot(numer);
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
    if (gamma.numerator().dot(rd.simpleCoroot(s))==0)
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

std::vector<deformation_term_tp>
Rep_table::deformation_terms (non_integral_block& block,BlockElt entry_elem)
{
  if (not block.survives(entry_elem)) // easy case, null result
    return std::vector<deformation_term_tp>();

  std::vector<BlockElt> survivors; survivors.reserve(entry_elem+1);
  for (BlockElt x=0; x<=entry_elem; ++x)
    if (block.survives(x))
      survivors.push_back(x);

  BlockElt n_surv = survivors.size(); // |BlockElt| now indexes |survivors|

  Rep_context rc(block.realGroup());
  std::vector<unsigned int> orient_nr(n_surv);
  for (BlockElt z=0; z<n_surv; ++z)
  {
    const BlockElt zz = survivors[z];
    const StandardRepr r =
      rc.sr(block.parent_x(zz),block.lambda_rho(zz),block.gamma());
    orient_nr[z] = rc.orientation_number(r);
  }

  // compute KL polynomimals $P_{x,y}$ with $x\leq y\leq \\{entry_elem}$
  // storing sign-changed value $P(x,y)=(-1)^{l(y)-l(x)}P_{x,y}$
  const kl::KLContext& klc = block.klc(entry_elem,false); // silently

  // sums of KL polynomials evaluated at $-1$, incorporating length signs
  int_Matrix P(n_surv,n_surv,0);

  // compute $P(x,z)$ for |xx<=zz<=entry_element| with |zz| among |survivors|
  // and contribute to |P(x,z)| where |zz=survivors[z]|, $xx\to survivors[x]$
  for (BlockElt z=n_surv; z-->0; )
  {
    const BlockElt zz = survivors[z];
    const unsigned int parity = block.length(zz)%2;
    for (BlockElt xx=0; xx <= zz; ++xx)
    {
      const kl::KLPol& pol = klc.klPol(xx,zz); // regular KL polynomial
      // evaluate |pol| at $X=-1$ since that is what's needed in the end
      int eval=0;
      for (polynomials::Degree d=pol.size(); d-->0; )
	eval = static_cast<int>(pol[d])-eval;
      if (eval!=0)
      {
	if (block.length(xx)%2!=parity)
	  eval = -eval; // incorporate sign for length difference in evaluation

        // contribute |eval| to all |P(x,z)| for which |xx| descends to
	// |survivors[x]| (expressing the singular $I(xx)$ as a survivor sum)
	BlockEltList nb=block.survivors_below(xx);
	for (size_t i=0; i<nb.size(); ++i)
	{
	  BlockElt x = std::lower_bound // look up |nb[i]| in |survivors|
	    (survivors.begin(),survivors.end(),nb[i])-survivors.begin();
	  assert(survivors[x]==nb[i]); // must be found
	  P(x,z) += eval;
	} // |for (i)| in |nb|
      } // |if(pol!=0)|
    } // |for (x<=z)|
  } // for |z|

  // now compute evaluated polynomials $Q_{x,z}$, for |x <= z==n_surv-1|
  // this is done by inverting the upper unitriangular matrix of the |P(x,z)|
  // only the final column of the inverse is needed

  const BlockElt z=n_surv-1; // $Q$ only needs |entry_elem| as second index
  int_Vector Q_z(n_surv);

  // we solve column $z$ of $Q$ from bottom to top, using the equation $PQ=1$
  Q_z[z]=1; // easy initial case, since $P(z,z)=1$
  for (BlockElt x=z; x-->0; )
  {
    int Q_xz=0;
    for (BlockElt y=x+1; y<=z; ++y) // isolate the term $Q(x,z)$ for $y=x$
      Q_xz -= P(x,y)*Q_z[y];  // set it to minus the sum of other terms
    Q_z[x]=Q_xz;
  }

  assert(survivors[z]==entry_elem);
  unsigned odd = (block.length(entry_elem)+1)%2; // opposite parity

  std::vector<deformation_term_tp> result;
  result.reserve(n_surv-1); // might be quite pessimistic, but not huge anyway

  // $\sum_{x\leq y<z}y[l(z)-l(y) odd] (-1)^{l(x)-l(y)}P_{x,y}*Q(y,z)$
  for (BlockElt x=n_surv-1; x-->0; ) // skip |entry_elem|
  {
    int coef=0;
    for (BlockElt y=x; y<n_surv-1; ++y)
      if (block.length(survivors[y])%2==odd) // length diff with |entry_elem|
	coef += P(x,y)*Q_z[y];
    // here we evaluated  |coef| at $X=-1$: "real" part of $(1-s)*coef[X:=s]$
    // if (orientation_difference(x,z)) coef*=s;
    int orient_express = (orient_nr[n_surv-1]-orient_nr[x])/2;
    if (coef!=0)
      result.push_back(deformation_term_tp
		       (orient_express%2!=0 ? -coef : coef,  survivors[x]));
  }

  return result;
} // |deformation_terms|


size_t Rep_table::rank() const { return context.rootDatum().rank(); }

SR_poly Rep_table::deformation(const StandardRepr& z)
{
  const Rep_context ctxt = rc();
  { unsigned long h=hash.find(z);
    if (h!=hash.empty)
      return def_formula[h];
  }

  Weight lam_rho = ctxt.lambda_rho(z);
  RatWeight nu =  ctxt.nu(z);

  StandardRepr z0 = ctxt.sr(z.x(),lam_rho,RatWeight(rank()));
  SR_poly result = ctxt.expand_final(z0); // value without deformation terms

  RationalList rp=ctxt.reducibility_points(z);
  for (unsigned i=rp.size(); i-->0; )
  {
    Rational r=rp[i];
    const StandardRepr zi = ctxt.sr(z.x(),lam_rho,nu*r);
    non_integral_block b(ctxt,zi);
    const RatWeight& gamma=b.gamma(); // NOT |zi.gamma()|, is made dominant!
    std::vector<deformation_term_tp> def_term = deformation_terms(b,b.size()-1);
    for (unsigned j=0; j<def_term.size(); ++j)
    {
      BlockElt eij=def_term[j].elt;
      StandardRepr zij = ctxt.sr(b.parent_x(eij),b.lambda_rho(eij),gamma);
      Split_integer coef(def_term[j].coef,-def_term[j].coef);
      result.add_multiple(deformation(zij),coef);
    }
  }

  // now store result for future lookup
  hash.match(z);
  def_formula.push_back(result);
  assert(hash.size()==def_formula.size());

  return result;
}

  } // |namespace repr|
} // |namespace atlas|
