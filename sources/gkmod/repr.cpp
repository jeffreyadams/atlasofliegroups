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

size_t Rep_context::rank() const { return rootDatum().rank(); }

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

StandardRepr Rep_context::sr(const non_integral_block& b, BlockElt i) const
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
  if (not block.survives(entry_elem) or block.length(entry_elem)==0)
    return std::vector<deformation_term_tp>(); // easy cases, null result

  ++deformations;

  // count number of survivors of length strictly less than any occurring length
  std::vector<unsigned int> n_surv_length_less(block.length(0),0);
  BlockEltList survivors; survivors.reserve(block.size());
  for (BlockElt x=0; x<block.size(); ++x)
  {
    if (block.length(x)==n_surv_length_less.size())
      n_surv_length_less.push_back(survivors.size());
    if (block.survives(x))
      survivors.push_back(x);
  }

  BlockElt ee = std::lower_bound // look up |entry_elem| in |survivors|
    (survivors.begin(),survivors.end(),entry_elem)-survivors.begin();
  assert(survivors[ee]==entry_elem); // must be found

  unsigned long hash_index=hash.find(sr(block,entry_elem));
  unsigned long cur_block; // set separately in both branches below

  if (hash_index==hash.empty) // previously unknown parameter
  {
    // first add an empty block to the table
    cur_block = block_list.size();
    block_list.push_back(KL_table());
    KL_table& KL=block_list.back();

    // fill the |location| table for new surviving parameters in this block
    for (BlockElt i=0; i<survivors.size(); ++i)
    {
      location here = { cur_block,i };
      assert(hash.size()==loc.size());
      unsigned long h=hash.match(sr(block,survivors[i]));
      if (h==loc.size())
	loc.push_back(loc_list(1,here));
      else
      {
	loc[h].push_back(here);
	++doublures;
      }
    }
    hash_index=hash.find(sr(block,entry_elem));
    assert(hash_index!=hash.empty);
    def_formula.resize(hash.size(),SR_poly(repr_less())); // allocate new slots

    // compute cumulated KL polynomimals $P_{x,y}$ with $x\leq y$ survivors
    // start with computing KL polynomials for the entire block
    const kl::KLContext& klc = block.klc(block.size()-1,false); // silently

    BlockElt n_surv = survivors.size(); // |BlockElt| now indexes |survivors|
    KL.reserve(n_surv);

    // get $P(x,z)$ for |xx<=zz<=entry_element| with |zz| among |survivors|
    // and contribute to |P(x,z)| where |zz=survivors[z]|, $xx\to survivors[x]$
    for (BlockElt zz=0; zz<n_surv; ++zz)
    {
      KL.push_back(std::vector<Split_integer>(zz+1,Split_integer(0)));
      std::vector<Split_integer>& KL_zz=KL.back();

      const BlockElt z = survivors[zz];
      const unsigned int parity = block.length(z)%2;
      for (BlockElt x=0; x <= z; ++x)
      {
	const kl::KLPol& pol = klc.klPol(x,z); // regular KL polynomial
	// evaluate |pol| at $X:=s$ since that is will be stored
	Split_integer eval(0);
	for (polynomials::Degree d=pol.size(); d-->0; )
	  eval = eval.times_s()+Split_integer(static_cast<int>(pol[d]));
	if (eval!=Split_integer(0))
	{
	  if (block.length(x)%2!=parity)
	    eval.negate(); // incorporate sign for length difference

	  // contribute |eval| to all |P(x,z)| for which |x| descends to
	  // |survivors[x]| (expressing the singular $I(x)$ as a survivor sum)
	  const BlockEltList nb=block.survivors_below(x);
	  for (BlockEltList::const_iterator it=nb.begin(); it!=nb.end(); ++it)
	  {
	    BlockElt xx = std::lower_bound // look up |*it| in |survivors|
	      (survivors.begin(),survivors.end(),*it)-survivors.begin();
	    assert(xx<n_surv and survivors[xx]==*it); // must be found
	    KL_zz[xx]+= eval;
	  } // |for (i)| in |nb|
	} // |if(pol!=0)|
      } // |for (x<=z)|
    } // |for(z)|
  } // |if(hash_index==hash.empty)|
  else
    cur_block=loc[hash_index][0].block;

  const KL_table& KL = block_list[cur_block];

  // now extract evaluations at $-1$ from previously filled table
  int_Matrix P(ee+1,ee+1,0);

  { // we must renumber |survivors| to point into block |cur_block|
    BlockEltList remap(ee+1);

    for (BlockElt xx=0; xx<=ee; ++xx)
    {
      unsigned long h=hash.find(sr(block,survivors[xx]));
      assert(h!=hash.empty); // certainly the parameter must be known now
      const loc_list& l=loc[h]; loc_list::const_reverse_iterator it;
      for (it=l.rbegin(); it!=l.rend(); ++it) // lookup, last is most probable
	if (it->block==cur_block)
	{ remap[xx]=it->elt; break; }
      assert(it!=l.rend()); // we should have found it
    }

    for (BlockElt zz=ee+1; zz-->0; )
      for (BlockElt xx=0; xx <= zz; ++xx)
	if (remap[xx]<=remap[zz]) // since only trangular part is stored
	{
	  const Split_integer val = KL[remap[zz]][remap[xx]];
	  P(xx,zz) = val.e()-val.s(); // evaluate stored value at $s:=-1$
	} // |for(xx)|, |for(zz)|
  }
  assert(P(ee,ee)==1); // since there can be no cumulation to the diagonal

  // now compute evaluated polynomials $Q_{xx,zz}$, for |xx <= zz==ee|
  // this is done by inverting the upper unitriangular matrix of the |P(xx,zz)|
  // only the final column of the inverse is needed

  int_Vector Q_ee(ee+1); // $Q$ only needs |ee| as second index

  // we solve column $z$ of $Q$ from bottom to top, using the equation $PQ=1$
  Q_ee[ee]=1; // easy initial case, since $P(ee,ee)=1$
  for (BlockElt xx=ee; xx-->0; )
  {
    int sum=0;
    for (BlockElt yy=xx+1; yy<=ee; ++yy) // skip term $Q(xx,ee)$ (for $yy=xx$)
      sum -= P(xx,yy)*Q_ee[yy];  // set it to minus the sum of other terms
    Q_ee[xx]=sum;
  }


  // $\sum_{x\leq y<ee}y[l(ee)-l(y) odd] (-1)^{l(x)-l(y)}P_{x,y}*Q(y,ee)$
  int_Vector coef(ee,0);
  unsigned int ll=block.length(entry_elem)-1; // last length of contributing |y|

  for (unsigned int l=ll%2; l<=ll; l+=2) // length of parity opposite |ee|
    for (BlockElt yy=n_surv_length_less[l]; yy<n_surv_length_less[l+1]; ++yy)
      for (BlockElt xx=0; xx<=yy; ++xx)
	coef[xx] += P(xx,yy)*Q_ee[yy]; // only contribute those terms

  std::vector<deformation_term_tp> result;
  result.reserve(ee); // might be quite pessimistic, but not huge anyway

  unsigned int orient_ee = orientation_number(sr(block,entry_elem));
  for (BlockElt xx=ee; xx-->0; )
    if (coef[xx]!=0)
    {
      int orient_express =
	(orient_ee-orientation_number(sr(block,survivors[xx])))/2;
      // if |orient_express| odd, we must multiply |coef| by $1-s$; however
      // since only real part of $(1-s)coef$ is stored, this amounts to negation
      if (orient_express%2!=0)
      coef[xx] = -coef[xx]; // factor $(1-s)$ will remain implicit in result

      result.push_back(deformation_term_tp(coef[xx],survivors[xx]));
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

  ++calls; // only count calls that come to this point

  StandardRepr z_near = sr(z.x(),lam_rho,nu_z*rp.back());
  make_dominant(z_near);

  unsigned long h=hash.find(z_near);
  { // look up if closest reducibility point to |z| is already known
    unsigned long h=hash.find(z_near);
    if (h!=hash.empty and not def_formula[h].empty())
    {
      ++ hits;
      return def_formula[h];
    }
  }

  for (unsigned i=rp.size(); i-->0; )
  {
    Rational r=rp[i];
    const StandardRepr zi = sr(z.x(),lam_rho,nu_z*r);
    non_integral_block b(*this,zi);
    std::vector<deformation_term_tp> def_term = deformation_terms(b,b.size()-1);
    for (unsigned j=0; j<def_term.size(); ++j)
    {
      StandardRepr zij = sr(b,def_term[j].elt);
      Split_integer coef(def_term[j].coef,-def_term[j].coef);
      result.add_multiple(deformation(zij),coef);
    }
  }

  // now store result for future lookup
  h=hash.find(z_near);
  assert(h!=hash.empty); // it should have been added by |deformation_terms|
  def_formula[h]=result;

  return result;
} // |Rep_table::deformation|

  } // |namespace repr|
} // |namespace atlas|
