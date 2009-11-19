/*!
\file
\brief Implementation for the classes
StandardRepK and KhatContext.
*/
/*
  This is standardrepk.cpp

  Copyright (C) 2006, 2007 Alfred Noel
  Copyright (C) 2008, 2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "standardrepk.h"

#include "constants.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "realredgp.h"
#include "kgb.h"
#include "tits.h"
#include "descents.h"
#include "lattice.h"
#include "rootdata.h"
#include "matreduc.h"
#include "intutils.h"
#include "ioutils.h"
#include "tags.h"
#include "graph.h"
#include "prettyprint.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <deque>

namespace atlas {

/*
     Chapter 0 -- Local function declarations
*/

namespace {

// an auxiliary function for height computations
latticetypes::LatticeMatrix
orth_projection(const rootdata::RootDatum& rd, bitset::RankFlags gens,
		latticetypes::LatticeCoeff& denom);

} // |namespace|

/*****************************************************************************

        Chapter I -- The HechtSchmid and StandardRepK classes

******************************************************************************/

namespace standardrepk {

Char HechtSchmid::rhs () const
{
  if (rh1==NULL)
    return Char();

  Char result(*rh1);
  if (rh2!=NULL)
    result+= Char(*rh2);

  return result;
}

Char HechtSchmid::equivalent () const
{
  Char result=rhs();
  if (lh2!=NULL)
    result -= Char(*lh2);

  return result;
}

// accessors

bool StandardRepK::operator< (const StandardRepK& rhs) const
{
  if (d_cartan != rhs.d_cartan) return d_cartan < rhs.d_cartan ;
  if (d_fiberElt != rhs.d_fiberElt) return d_fiberElt < rhs.d_fiberElt;
  return d_lambda < rhs.d_lambda;
}

bool StandardRepK::operator== (const StandardRepK& other) const
{
  return d_cartan == other.d_cartan
    and d_fiberElt == other.d_fiberElt
    and d_lambda == other.d_lambda;
  // when these components match, the restrictions to $K$ will be the same
}

size_t StandardRepK::hashCode(size_t modulus) const
{
  size_t hash=13*d_fiberElt.data().to_ulong()+d_cartan;
  for (size_t i=0; i<d_lambda.first.size(); ++i)
    hash+=(hash<<2)+d_lambda.first[i];
  return (hash+(hash<<5)+d_lambda.second.to_ulong())&(modulus-1);
}

} // namespace standardrepk

/*****************************************************************************

        Chapter II -- The SR_rewrites class

******************************************************************************/

namespace standardrepk {

const combination& SR_rewrites::lookup(seq_no n) const
{
  assert(n<system.size());
  return system[n];
}

void SR_rewrites::equate (seq_no n, const combination& rhs)
{
  assert(n==system.size()); // left hand side should be a new |StandardRepK|
  system.push_back(rhs);
}

} // namespace standardrepk



/*****************************************************************************

        Chapter III -- The KhatContext class

******************************************************************************/

namespace standardrepk {

SRK_context::SRK_context(realredgp::RealReductiveGroup &GR)
  : G(GR.complexGroup())
  , Tg(GR.basedTitsGroup())
  , Cartan_set(GR.Cartan_set())
  , C_info(G.numCartanClasses())
  , simple_reflection_mod_2()
{
  const rootdata::RootDatum& rd=rootDatum();
  simple_reflection_mod_2.reserve(G.semisimpleRank());
  for (size_t i=0; i<G.semisimpleRank(); ++i)
    simple_reflection_mod_2.push_back
      (latticetypes::BinaryMap(rd.simple_reflection(i).transposed()));

  size_t n = rootDatum().rank();

  size_t nr=0;
  for (bitmap::BitMap::iterator it=Cartan_set.begin(); it(); ++it,++nr)
  {

    // d_G.cartan[i] is (canonical involution for) (*it)th CartanClass
    Cartan_info& ci=C_info[nr];

    const latticetypes::LatticeMatrix& theta = G.cartan(*it).involution();


    // put in $q$ the matrix of $\theta-1$
    latticetypes::LatticeMatrix q=theta;
    for (size_t i=0; i<n; ++i)
      q(i,i) -= 1;

    // find basis adapted to image of $\theta-1$
    latticetypes::CoeffList factor;
    latticetypes::LatticeMatrix basis = matreduc::adapted_basis(q,factor);
    latticetypes::LatticeMatrix inv_basis=basis.inverse();

    size_t twos=0,l=factor.size();
    bitset::RankFlags torsion;

    for (size_t i=0; i<l; ++i)
      if (factor[i]>1)
      {
	assert(factor[i]==2);
	torsion.set(i);
	++twos;
      }

    // get torsion columns of |basis| to |torsionLift|
    // and rows of |inv_basis|, modulo 2, to |torsionProjector|
    ci.torsionLift.reserve(twos);
    ci.torsionProjector.resize(twos,n);
    size_t i=0;
    for (bitset::RankFlags::iterator it=torsion.begin(); it(); ++it,++i)
    {
      ci.torsionLift.push_back(basis.column(*it));

      for (size_t j=0; j<n; ++j)
	ci.torsionProjector.set_mod2(i,j,inv_basis(*it,j));
    }

    basis.block(0,l,n,n).swap(ci.freeLift); // final |n-l| basis vectors
    inv_basis.block(l,0,n,n).swap(ci.freeProjector); // final |n-l| rows

    latticetypes::SmallSubspace
      (latticetypes::BinaryMap(lattice::eigen_lattice(theta.transposed(),-1)))
      .swap(ci.fiber_modulus);
  } // |for (it)|
}


HCParam SRK_context::project
  (size_t cn, latticetypes::Weight lambda) const
{
  const Cartan_info& ci=info(cn);

  (lambda -= rootDatum().twoRho()) /= 2; // final division is exact

  // now |lambda| actually represents $\lambda-\rho$ in plain coordinates
  return std::make_pair
    (ci.freeProjector.apply(lambda),
     ci.torsionProjector.apply(latticetypes::SmallBitVector(lambda)).data()
     );
}

latticetypes::Weight SRK_context::lift(size_t cn, HCParam p) const
{
  const Cartan_info& ci=info(cn);
  latticetypes::Weight result=ci.freeLift.apply(p.first); // lift free part

  latticetypes::WeightList torsion_lift=ci.torsionLift;
  for (size_t i=0; i<torsion_lift.size(); ++i)
    if (p.second[i])
      result += torsion_lift[i]; // add even vectors representing torsion part
  (result*=2) += rootDatum().twoRho(); // convert $\lambda-\rho$ to $2\lambda$

  return result;
}

StandardRepK SRK_context::std_rep
  (const latticetypes::Weight& two_lambda, tits::TitsElt a) const
{
  const weyl::WeylGroup& W=weylGroup();
  const rootdata::RootDatum& rd=rootDatum();

  weyl::TwistedInvolution sigma=a.tw();
  weyl::WeylElt w = complexGroup().canonicalize(sigma);
  // now |sigma| is canonical and |w| conjugates |sigma| to |a.tw()|


  // conjugate towards canonical element (left-right, for inverse conjugation)
  basedTitsGroup().basedTwistedConjugate(a,W.word(w));

  assert(a.tw()==sigma); // we should now be at canonical twisted involution

  latticetypes::Weight mu=W.imageByInverse(rd,w,two_lambda);

  size_t cn = complexGroup().class_number(sigma);
  StandardRepK result(cn,
		      info(cn).fiber_modulus.mod_image
		        (titsGroup().left_torus_part(a)),
		      project(cn,normalize(mu,cn)));

  return result;
} // |std_rep|

// the following is a variant of |std_rep_rho_plus| intended for |KGB_sum|
// it should only transform the parameters for the Levi factor given by |gens|
// since |lambda| is $\rho$-centered, care should be taken in transforming it
RawRep SRK_context::Levi_rep
    (latticetypes::Weight lambda, tits::TitsElt a, bitset::RankFlags gens)
  const
{
  weyl::TwistedInvolution sigma=a.tw();
  weyl::WeylElt w = complexGroup().canonicalize(sigma,gens);
  // now |sigma| is canonical for |gens|, and |w| conjugates it to |a.tw()|

  const rootdata::RootDatum& rd=rootDatum();

  // conjugate towards canonical element
  {
    weyl::WeylWord ww=weylGroup().word(w);
    for (size_t i=0; i<ww.size(); ++i) // left-to-right for inverse conjugation
    {
      assert(gens.test(ww[i])); // check that we only used elements in $W(L)$
      basedTitsGroup().basedTwistedConjugate(a,ww[i]);
      rd.simpleReflect(lambda,ww[i]);
      lambda-=rd.simpleRoot(ww[i]); // make affine reflection fixing $-\rho$
    }
  }
  assert(a.tw()==sigma); // we should now be at canonical twisted involution

  return RawRep (lambda,a);
} // |Levi_rep|

bool SRK_context::isStandard(const StandardRepK& sr, size_t& witness) const
{
  const rootdata::RootDatum& rd=rootDatum();
  latticetypes::Weight lambda=lift(sr);
  const cartanclass::Fiber& f=fiber(sr);

  for (size_t i=0; i<f.imaginaryRank(); ++i)
    if (lambda.scalarProduct(rd.coroot(f.simpleImaginary(i)))<0)
    {
      witness=i; return false;
    }

  return true;
}

latticetypes::Weight SRK_context::normalize
  (latticetypes::Weight lambda, size_t cn) const
{
  const rootdata::RootDatum& rd = rootDatum();
  const cartanclass::Fiber& f=G.cartan(cn).fiber();

  bitset::RankFlags bi_ortho_simples;
  { // find simple roots orthogonal to |real2rho| and |imaginary2rho|
    latticetypes::Weight real2rho=rd.twoRho(f.realRootSet());
    latticetypes::Weight imaginary2rho=rd.twoRho(f.imaginaryRootSet());
    for (size_t i=0; i<rd.semisimpleRank(); ++i)
      if (rd.isOrthogonal(real2rho,rd.simpleRootNbr(i)) and
	  rd.isOrthogonal(imaginary2rho,rd.simpleRootNbr(i)))
	bi_ortho_simples.set(i);
  }

  const latticetypes::LatticeMatrix& theta = f.involution();
  rootdata::RootSet cplx = f.complexRootSet();

  //  go through the orthogonal list
  //  select the complex roots in the list

  bitset::RankFlags::iterator i;
  do
  {
    latticetypes::Weight mu=theta.apply(lambda);
    mu +=lambda; // $\mu=(1+\theta)\alpha$, simplifies scalar product below

    for (i= bi_ortho_simples.begin(); i(); ++i )
    {
      rootdata::RootNbr alpha = rd.simpleRootNbr(*i);

      if (cplx.isMember(alpha) and rd.scalarProduct(mu,alpha)<0)
      {
	rootdata::RootNbr beta= f.involution_image_of_root(alpha);
	assert (rd.isOrthogonal(alpha,beta));
	rd.reflect(lambda,alpha);
	rd.reflect(lambda,beta);
	break; // and continue do-while loop
      }
    } // for
  }
  while (i()); // while |for|-loop was exited through |break|

  return lambda;
} // normalize

bool SRK_context::isZero(const StandardRepK& sr, size_t& witness) const
{
  const rootdata::RootDatum& rd=rootDatum();
  latticetypes::Weight lambda=lift(sr);
  const cartanclass::Fiber& f=fiber(sr);
  tits::TitsElt a=titsElt(sr);

  for (size_t i=0; i<f.imaginaryRank(); ++i)
    if (not basedTitsGroup().grading(a,f.simpleImaginary(i)) // i.e., compact
	and lambda.scalarProduct(rd.coroot(f.simpleImaginary(i)))==0)
    {
      witness=i; return true;
    }

  return false;
}

bool SRK_context::isFinal(const StandardRepK& sr, size_t& witness) const
{
  const rootdata::RootDatum& rd=rootDatum();
  latticetypes::Weight lambda=lift(sr);
  const cartanclass::Fiber& f=fiber(sr);

  // since coordinates are doubled, the scalar product below is always even
  for (size_t i=0; i<f.realRank(); ++i)
    if (lambda.scalarProduct(rd.coroot(f.simpleReal(i)))%4 == 0)
    {
      witness=i; return false;
    }

  return true;
}


std::ostream& SRK_context::print(std::ostream& strm,const StandardRepK& sr)
  const
{
  prettyprint::printVector(strm,lift(sr)) << '@';
  prettyprint::prettyPrint(strm,sr.d_fiberElt) << '#' << sr.d_cartan;
  return strm;
}

std::ostream& SRK_context::print(std::ostream& strm,const Char& ch) const
{
  if (ch.empty())
    return strm << '0';
  for (Char::const_iterator it=ch.begin(); it!=ch.end(); ++it)
  {
    strm << (it->second>0 ? " + " : " - ");
    long int ac=intutils::abs<long int>(it->second);
    if (ac!=1)
      strm << ac << '*';
    print(strm,it->first);
  }
  return strm;
}

std::ostream& SRK_context::print(std::ostream& strm,const q_Char& ch) const
{
  if (ch.empty())
    return strm << '0';
  for (q_Char::const_iterator it=ch.begin(); it!=ch.end(); ++it)
  {
    if (it->second.degree()==0)
    {
      strm << (it->second[0]>0 ? " + " : " - ");
      long int ac=intutils::abs(it->second[0]);
      if (ac!=1)
	strm << ac << '*';
    }
    else
      prettyprint::printPol(strm<<" + (",it->second,"q")<<")*";
    print(strm,it->first);
  }
  return strm;
}

level
SRK_context::height(const StandardRepK& sr) const
{
  const rootdata::RootDatum& rd=rootDatum();
  const latticetypes::Weight mu=theta_lift(sr);

  level sum=0;
  for (rootdata::WRootIterator
	 it=rd.beginPosCoroot(); it!=rd.endPosCoroot(); ++it)
    sum +=intutils::abs(mu.scalarProduct(*it));

  return sum/2; // each |scalarProduct| above is even (in doubled coordinates)
}


/*
 *
 *
 *           |KhatContext|
 */

KhatContext::KhatContext
  (realredgp::RealReductiveGroup &GR, const kgb::KGB& kgb)
    : SRK_context(GR)
    , d_KGB(kgb)
    , nonfinal_pool(),final_pool()
    , nonfinals(nonfinal_pool), finals(final_pool)
    , height_of()
    , height_graded(height_of)
    , d_rules()
    , proj_pool(), proj_sets(proj_pool), proj_data()
{}

/******** accessors *******************************************************/


std::ostream& KhatContext::print(std::ostream& strm,
				 const combination& ch,
				 bool brief) const
{
  if (ch.empty())
    return strm << '0';
  for (combination::const_iterator it=ch.begin(); it!=ch.end(); ++it)
  {
    strm << (it->second>0 ? " + " : " - ");
    long int ac=intutils::abs<long int>(it->second);
    if (ac!=1)
      strm << ac << '*';
    if (brief)
      strm << 'R' << it->first;
    else print(strm,rep_no(it->first));
  }
  return strm;
}

std::ostream& KhatContext::print
  (std::ostream& strm, const q_combin& ch, bool brief)  const
{
  if (ch.empty())
    return strm << '0';
  for (q_combin::const_iterator it=ch.begin(); it!=ch.end(); ++it)
  {
    if (it->second.degree()==0)
    {
      strm << (it->second[0]>0 ? " + " : " - ");
      long int ac=intutils::abs(it->second[0]);
      if (ac!=1)
	strm << ac << '*';
    }
    else
      prettyprint::printPol(strm<<" + (",it->second,"q")<<")*";
    if (brief)
      strm << 'R' << it->first;
    else print(strm,rep_no(it->first));
  }
  return strm;
}

/* map a character (with terms assumed to be normalized) to one containing
   only Standard terms
   this version ensures the basic |standardize| is recursively called first */
combination KhatContext::standardize(const Char& chi)
{
  combination result(height_graded);

  for (Char::const_iterator i=chi.begin(); i!=chi.end(); ++i)
    result.add_multiple(standardize(i->first),i->second);

  return result;
}

// the basic case, |sr| is assumed normalized here
combination KhatContext::standardize(const StandardRepK& sr)
{
  { // first check if we've already done |sr|
    seq_no n=nonfinals.find(sr);
    if (n!=Hash::empty)
      return d_rules.lookup(n); // in this case an equation should be known
  }

  size_t witness;
  if (isStandard(sr,witness))
  {
    { // now check if we already know |sr| to be Final
      seq_no n=finals.find(sr);
      if (n!=Hash::empty)
	return combination(n,height_graded); // single term known to be final
    }

    if (isZero(sr,witness))
    {
      combination zero(height_graded);
      d_rules.equate(nonfinals.match(sr),zero);
      return zero;
    }

    if (isFinal(sr,witness))
    {
      assert(height_of.size()==final_pool.size());
      height_of.push_back(height(sr));
      return combination(finals.match(sr),height_graded); // single term
    }

    // now |sr| is known to be Standard, but neither Zero nor Final

    HechtSchmid equation= back_HS_id(sr,fiber(sr).simpleReal(witness));
    assert(equation.n_lhs()==1); // |back_HS_id| gives 1-term left hand side
    assert(equation.n_rhs()!=0); // and never a null right hand side

    // now recursively standardize all terms, storing rules
    combination result= standardize(equation.rhs());
    d_rules.equate(nonfinals.match(sr),result); // and add rule for |sr|
    return result;
  } // if (isStandard(sr,witness))

  HechtSchmid equation= HS_id(sr,fiber(sr).simpleImaginary(witness));
  assert(equation.n_lhs()==2); // all cases of |HS_id| produce 2-term lhs

  // recursively standardize all terms of |equation.equivalent()|, storing rules
  combination result= standardize(equation.equivalent());
  d_rules.equate(nonfinals.match(sr),result); // and add rule for |sr|
  return result;
} // standardize



HechtSchmid
KhatContext::HS_id(const StandardRepK& sr, rootdata::RootNbr alpha) const
{
  HechtSchmid id(sr,*this);
  const rootdata::RootDatum& rd=rootDatum();
  tits::TitsElt a=titsElt(sr);
  latticetypes::Weight lambda=lift(sr);
  assert(rd.isPosRoot(alpha)); // indeed |alpha| simple-imaginary for |a.tw()|

  size_t i=0; // simple root index (value will be set in following loop)
  while (true) // we shall exit halfway when $\alpha=\alpha_i$
  {
    while (not rd.is_descent(i,alpha))
    {
      ++i;
      assert(i<rd.semisimpleRank());
    }
    // now $\<\alpha,\alpha_i^\vee> > 0$ where $\alpha$ is simple-imaginary
    // and \alpha_i$ is complex for the involution |a.tw()|

    if (alpha==rd.simpleRootNbr(i)) break; // found it

    // otherwise reflect all data by $s_i$, which decreases level of $\alpha$
    rd.simple_reflect_root(alpha,i);
    rd.simpleReflect(lambda,i);
    basedTitsGroup().basedTwistedConjugate(a,i);
    i=0; // and start over
  }

  latticetypes::Weight mu=rd.simpleReflection(lambda,i);
  if (basedTitsGroup().simple_grading(a,i))
  { // $\alpha_i$ is a non-compact imaginary simple root
    basedTitsGroup().basedTwistedConjugate(a,i); // adds $m_i$ to torus part
    StandardRepK sr0= std_rep(mu, a);
    assert(sr.d_cartan==sr0.d_cartan);
    id.add_lh(sr0,*this);
    // the change to |d_lambda| may involve both components

    /* Now are equivalent:
       = |sr0.d_fiberElt==sr.d_fiberElt|
       = $m_i$ absorbed into quotient ($\alpha_i^vee \in 2X_*+X_*^{-\theta}$)
       = $\alpha^\vee( (X^*)^\theta ) = 2\Z$
       = $\alpha\notin(1-\theta_\alpha)X^*$  ($\theta_\alpha=\theta.s_\alpha$)
       = Cayley transform is type II
    */

    basedTitsGroup().Cayley_transform(a,i);
    StandardRepK sr1= std_rep(lambda, a);
    id.add_rh(sr1,*this);
    if (sr0.d_fiberElt==sr.d_fiberElt) // type II
    {
      lambda += rootDatum().simpleRoot(i); // other possibility
      lambda += rootDatum().simpleRoot(i); // add twice: doubled coordinates
      StandardRepK sr2= std_rep(lambda, a);
      assert(sr1.d_lambda != sr2.d_lambda);
      id.add_rh(sr2,*this);
    }
  }
  else // $\alpha_i$ is a compact root; "easy" Hecht-Schmid identity
    // based twisted conjugation fixed the Tits element; just reflect weight
    id.add_lh(std_rep(mu,a),*this); // and no RHS

  return id;
} // HS_id

/*
  The method |HS_id| is only called when |lambda| is non-dominant for |alpha|
  and therefore gives a second term with a strictly more dominant weight.

  For those simple-imaginary $\alpha$ with $\<\lambda,\alpha^\vee>=0$ the
  corresponding Hecht-Schmid identity has equal weights in the left hand
  terms, and is never generated from |HS_id|, but it is nevertheless valid,
  and (for Standard representations) useful. Its use is as follows:

  - if $\alpha$ is compact, the easy HS identity makes (twice) the
    Standard representation equal to 0.

  - if $\alpha$ is non-compact, the left hand terms are either distinct
    because of their fiber part (type I) or identical (type II).

    - in the type I case the right hand side is a single term whose (modularly
      reduced) weight takes an even value on the (now real) root alpha, and is
      therefore non-final. The identity can be used to express that
      representation as the sum of the (Standard) left hand terms.

    - in the type II case the left hand terms are equal and the two right hand
      terms have distinct non-final parameters (differing in the weight by
      $\alpha$), but designate the same representation (due to a shifted
      action of the real Weyl group). One may therefore equate them and divide
      by 2, so that either of the right hand terms is equated to the original
      representation.

   The purpose of |back_HS_id| is generating the equations for the type I and
   type II cases, for a non-final parameter |sr| and witnessing root |alpha|
 */
HechtSchmid
KhatContext::back_HS_id(const StandardRepK& sr, rootdata::RootNbr alpha) const
{
  HechtSchmid id(sr,*this);
  const rootdata::RootDatum& rd=rootDatum();
  tits::TitsElt a=titsElt(sr);
  latticetypes::Weight lambda=lift(sr);
  assert(rd.isPosRoot(alpha)); // in fact it must be simple-real for |a.tw()|

  // basis used is of $(1/2)X^*$, so scalar product with coroot always even
  assert(rd.scalarProduct(lambda,alpha)%4 == 0); // the non-final condition

  {
    latticetypes::Weight mu=rd.root(alpha);
    mu *= rd.scalarProduct(lambda,alpha)/2; // an even multiple of $\alpha$
    lambda -= mu; // this makes |lambda| invariant under $s_\alpha$
    // due to basis used, |lambda| effectively modified by multiple of $\alpha$
    /* in type I, $\alpha$ is in $(1-\theta)X^*$ and correction is neutral
       in type II, correction need not be in $(1-\theta)X^*$, but adding
       $\alpha$ gives HC parameter designating the same representation
    */
  }

  latticetypes::SmallSubspace mod_space=
    info(sr.d_cartan).fiber_modulus; // make a copy to be modified

  // again, move to situation where $\alpha$ is simple: $\alpha=\alpha_i$
  size_t i=0; // simple root index (value will be set in following loop)

  while (true) // we shall exit halfway when $\alpha=\alpha_i$
  {
    while (not rd.is_descent(i,alpha))
    {
      ++i;
      assert(i<rd.semisimpleRank());
    }
    // now $\<\alpha,\alpha_i^\vee> > 0$ where $\alpha$ is simple-real
    // and \alpha_i$ is complex for the involution |a.tw()|

    if (alpha==rd.simpleRootNbr(i)) break; // found it

    // otherwise reflect all data by $s_i$, which decreases level of $\alpha$
    rd.simple_reflect_root(alpha,i);
    rd.simpleReflect(lambda,i);
    basedTitsGroup().basedTwistedConjugate(a,i);
    mod_space.apply(dual_reflection(i));
    i=0; // and start over
  }

  // one right term is obtained by undoing Cayley for |a|, with lifted |lambda|
  basedTitsGroup().inverse_Cayley_transform(a,i,mod_space);

  StandardRepK sr1 = std_rep(lambda,a);
  id.add_rh(sr1,*this);

  // there will be another term in case of a type I HechtSchmid identity
  // it differs only in the fiber part, by $m_i$; if this vanishes into the
  // quotient, the HechtSchmid identity is type II and nothing is added
  basedTitsGroup().basedTwistedConjugate(a,i);
  StandardRepK sr2=std_rep(lambda,a);
  if (sr1.d_fiberElt!=sr2.d_fiberElt) // type I
    id.add_rh(sr2,*this);

  return id;
} // |back_HS_id|


const proj_info& KhatContext::get_projection(bitset::RankFlags gens)
{
  size_t old_size=proj_data.size();
  size_t h=proj_sets.match(bitset_entry(gens));
  if (h<old_size)
    return proj_data[h];

  proj_data.push_back(proj_info());
  proj_data.back().projection=
    orth_projection(rootDatum(),gens,proj_data.back().denom);
  return proj_data.back();
}

level KhatContext::height_bound(const latticetypes::Weight& lambda)
{
  const rootdata::RootDatum& rd=rootDatum();

  bitset::RankFlags negatives,new_negatives;
  latticetypes::Weight mu;

  do
  {
    new_negatives.reset();
    mu=get_projection(negatives).projection.apply(lambda);
    for (size_t i=0; i<rd.semisimpleRank(); ++i)
      if (not negatives[i] and mu.scalarProduct(rd.simpleCoroot(i))<0)
	new_negatives.set(i);
    negatives |= new_negatives;
  }
  while (new_negatives.any());

  level sp=mu.scalarProduct(rd.dual_twoRho());
  level d=2*get_projection(negatives).denom; // double to match |sum/2| above
  return (sp+d-1)/d; // round upwards, since height is always integer
}

/*!
  The purpose of |theta_stable_parabolic| is to move to a situation in which
  |dom| is dominant, and in which the only positive roots sent to negative
  ones by the involution $\theta$ are the real ones. Then the real roots
  define the Levi factor, and the remaining positive roots the radical of a
  $\theta$-stable parabolic subalgebra. The involution $\theta$ is given by
  the twisted involution in |strong|; we think of keeping it fixed while
  gradually moving around the set of positive roots, each time replacing some
  complex root in the simple basis by its opposite. In realitity the positive
  system is fixed, and the moves conjugate $\theta$ and reflect |dom|. In fact
  we shall need a fiber part as well as an involution once we have obtained
  our goal, so conjugetion is in fact |basedTwistedConjugate| on |strong|.
 */
PSalgebra
KhatContext::theta_stable_parabolic
  (const StandardRepK& sr, weyl::WeylWord& conjugator) const
{
  const rootdata::RootDatum& rd=rootDatum();
  const weyl::TwistedWeylGroup& W=twistedWeylGroup();

  latticetypes::Weight dom=theta_lift(sr);
  tits::TitsElt strong=titsElt(sr);

  conjugator.resize(0); // clear this output parameter

  /* the following loop terminates because we either increase the number of
     positive coroots with strictly positive evaluation on |dom|, or we keep
     that number constant and decrease the number of positive complex roots
     that the involution sends to negative ones.
  */
  while (true) // loop will terminate if inner loop runs to completion
  {
    size_t i;
    for (i=0; i<rd.semisimpleRank(); ++i)
    {
      rootdata::RootNbr alpha=rd.simpleRootNbr(i);
      latticetypes::LatticeCoeff v=dom.scalarProduct(rd.simpleCoroot(i));

      if (v<0) // first priority: |dom| should be made dominant
	break; // found value of |i| to used in conjugation/reflection
      else if (v>0) continue; // don't touch |alpha| in this case

      // now |dom| is on reflection hyperplan for |alpha|

      // second priority give |alpha| and its $\theta$ image the same sign
      rootdata::RootNbr beta= // image of |alpha| by $\theta$
	rd.permuted_root(W.word(strong.w()),rd.simpleRootNbr(W.twisted(i)));
      if (not rd.isPosRoot(beta) and beta!=rd.rootMinus(alpha))
	break; // found |i| in this case as well

    } // for i

    if (i<rd.semisimpleRank()) // then we found a reflection |i| to apply
    {
      basedTitsGroup().basedTwistedConjugate(strong,i);
      rd.simpleReflect(dom,i);
      conjugator.push_back(i);
    }
    else break; // no simple roots give any improvement any more, so stop
  } // |while(true)|

/*
   We have achieved that any real positive root is a sum of real simple roots.
   Here's why. Consider a counterexample that is minimal: a real positive root
   $\alpha$ from which no real simple root can be subtracted while remaining
   positive. Then this must be a sum of simple roots that are either complex
   or imaginary; the $\theta$-images of such simple roots are all positive,
   but their sum is $-\alpha$ which is negative, a contradiction.

   Also |conjugator| is such that |W.twistedConjugate(strong.tw(),conjugator)|
   would make |strong.tw()| equal to its original value again.
*/

  // Build the parabolic subalgebra:

  { // first ensure |strong| is in reduced
    const latticetypes::LatticeMatrix theta =
      complexGroup().involutionMatrix(strong.tw());
    titsGroup().left_torus_reduce(strong,tits::fiber_denom(theta));
  }

  return PSalgebra(strong,d_KGB,complexGroup());

} // theta_stable_parabolic

kgb::KGBEltList KhatContext::sub_KGB(const PSalgebra& q) const
{
  bitmap::BitMap flagged(d_KGB.size());
  tits::TitsElt strong=q.strong_involution();

  kgb::KGBElt root;
  {
    kgb::KGBEltPair p=d_KGB.tauPacket(q.involution());
    kgb::KGBElt x;
    for (x=p.first; x<p.second; ++x)
      if (d_KGB.titsElt(x)==q.strong_involution())
      {
	root=x; break;
      }
    assert(x<p.second); // search should succeed
  }

  flagged.insert(root);
  std::deque<kgb::KGBElt> queue(1,root);
  do
  {
    kgb::KGBElt x=queue.front(); queue.pop_front();
    for (bitset::RankFlags::iterator it=q.Levi_gens().begin(); it(); ++it)
    {
      kgb::KGBElt y=d_KGB.cross(*it,x);
      if (not flagged.isMember(y))
      {
	flagged.insert(y); queue.push_back(y);
      }
      y=d_KGB.inverseCayley(*it,x).first; // second will follow if present
      if (y!=kgb::UndefKGB and not flagged.isMember(y))
      {
	flagged.insert(y); queue.push_back(y);
      }
    }
  }
  while (not queue.empty());

  return kgb::KGBEltList(flagged.begin(),flagged.end());
} // sub_KGB

RawChar KhatContext::KGB_sum(const PSalgebra& q,
			     const latticetypes::Weight& lambda) const
{
  const rootdata::RootDatum& rd=rootDatum();
  kgb::KGBEltList sub=sub_KGB(q); std::reverse(sub.begin(),sub.end());

  std::vector<size_t> sub_inv(d_KGB.size(),~0);

  for (size_t i=0; i<sub.size(); ++i)
    sub_inv[sub[i]]=i; // partially fill array with inverse index

  std::vector<latticetypes::Weight> mu; // list of $\rho$-centered weights,
  mu.reserve(sub.size());               // associated to the elements of |sub|

  mu.push_back(lambda); (mu[0]-=rd.twoRho())/=2; // make $\rho$-centered

  for (size_t i=1; i<sub.size(); ++i)
  {
    kgb::KGBElt x=sub[i];
    bitset::RankFlags::iterator it;
    for (it=q.Levi_gens().begin(); it(); ++it)
    {
      if (d_KGB.cross(*it,x)>x) // then we can use ascending cross action
      {
	size_t k=sub_inv[d_KGB.cross(*it,x)];
	assert(k!=~0ul); // we ought to land in the subset
	mu.push_back(rd.simpleReflection(mu[k],*it)); // $\rho$-centered
	break;
      }
    }
    if (it()) continue; // if we could use a cross action, we're done for |i|

    // now similarly try Cayley transforms
    for (it=q.Levi_gens().begin(); it(); ++it)
    {
      if (d_KGB.cayley(*it,x)!=kgb::UndefKGB) // then we can use this Cayley
      {
	size_t k=sub_inv[d_KGB.cayley(*it,x)];
	assert(k!=~0ul); // we ought to land in the subset
	latticetypes::Weight nu=mu[k]; // $\rho-\lambda$ at split side
	assert(nu.scalarProduct(rd.simpleCoroot(*it))%2 == 0); // finality
	latticetypes::Weight alpha=rd.simpleRoot(*it);
	nu -= (alpha *= nu.scalarProduct(rd.simpleCoroot(*it))/2); // project
	mu.push_back(nu); // use projected weight at compact side of transform
	break;
      }
    }
    assert(it()); // if no cross action worked, some Cayley transform must have
  }
  assert(mu.size()==sub.size());

  size_t max_l=d_KGB.length(sub[0]);

  RawChar result;
  for (size_t i=0; i<sub.size(); ++i)
  {
    kgb::KGBElt x=sub[i];
    RawRep r(mu[i],d_KGB.titsElt(x));
    result += RawChar(r, ((max_l-d_KGB.length(x))%2 == 0 ? 1 : -1));
  }

  return result;
} // KGB_sum

combination KhatContext::truncate(const combination& c, level bound) const
{
  combination result(height_graded);
  for (combination::const_iterator it=c.begin(); it!=c.end(); ++it)
    if (height(it->first)<=bound)
      result.insert(result.end(),*it);

  return result;
}

// Express irreducible K-module as a finite virtual sum of standard ones
CharForm
KhatContext::K_type_formula(const StandardRepK& sr, level bound)
{
  const weyl::WeylGroup& W=weylGroup();
  const rootdata::RootDatum& rd=rootDatum();

  // Get theta stable parabolic subalgebra

  weyl::WeylWord conjugator;
  PSalgebra q = theta_stable_parabolic(sr,conjugator);

  latticetypes::Weight lambda=
    W.imageByInverse(rd,W.element(conjugator),lift(sr));

  RawChar KGB_sum_q= KGB_sum(q,lambda);

  // type of formal linear combination of weights, associated to Tits element

  Char result;
  for (RawChar::const_iterator it=KGB_sum_q.begin(); it!=KGB_sum_q.end(); ++it)
  {
    Char::coef_t c=it->second; // coefficient from |KGB_sum_q|
    const latticetypes::Weight& mu=it->first.first; // weight from |KGB_sum_q|
    const tits::TitsElt& strong=it->first.second; // Tits elt from |KGB_sum_q|
    cartanclass::InvolutionData id(complexGroup(),strong.tw());

    rootdata::RootSet A(rd.numRoots());
    for (bitmap::BitMap::iterator
	   rt=q.radical().begin(); rt!=q.radical().end(); ++rt)
    {
      rootdata::RootNbr alpha=*rt;
      assert(not id.real_roots().isMember(alpha));
      if (id.imaginary_roots().isMember(alpha))
	A.set_to(alpha,basedTitsGroup().grading(strong,alpha)); // add if nc
      else // complex root
      {
	rootdata::RootNbr beta=id.root_involution(alpha);
	assert(rd.isPosRoot(beta));
	A.set_to(alpha,beta>alpha); // add first of two complex roots
      }
    }

//     std::cout << "Sum over subsets of " << A.size() << " roots, giving ";

    typedef free_abelian::Monoid_Ring<latticetypes::Weight> polynomial;
    const latticetypes::LatticeMatrix theta =
      complexGroup().involutionMatrix(strong.tw());

    // compute $X^\mu*\prod_{\alpha\in A}(1-X^\alpha)$ in |pol|
    polynomial pol(mu);
    for (rootdata::RootSet::iterator it=A.begin(); it!=A.end(); ++it)
    {
      polynomial copy=pol; // since |add_multiple| assumes no aliasing
      pol.add_multiple(copy,-1,rd.root(*it));

      // filter out terms that cannot affect anything below |bound|
      for (polynomial::iterator term=pol.begin(); term!=pol.end();)
      {
	latticetypes::Weight lambda=term->first;
	(lambda*=2) += rd.twoRho();
	lambda += theta.apply(lambda);
	if (height_bound(lambda)>bound)
	  pol.erase(term++);
	else
	  term++;
      }
    }
//     std::cout << pol.size() << " terms." << std::endl;

    // iterate over terms in formal sum, taking coef *= |c|
    for (polynomial::const_iterator term=pol.begin(); term!=pol.end(); ++term)
    {
      latticetypes::Weight lambda=term->first;
      polynomial::coef_t coef=term->second;
      result += Char(std_rep_rho_plus(lambda,strong),c*coef); // contribute
    }
  } // for sum over KGB for L
  return std::make_pair(sr, result);
} // K_type_formula

// Apply |K_type_formula| for known Final representation, and |standardize|
equation KhatContext::mu_equation(seq_no n, level bound)
{
  CharForm kf= K_type_formula(rep_no(n),bound);

  equation result(n,combination(height_graded));
  combination& sum=result.second;

  for (Char::const_iterator it=kf.second.begin(); it!=kf.second.end(); ++it)
    sum.add_multiple(truncate(standardize(it->first),bound),it->second);

  return result;
}

Raw_q_Char KhatContext::q_KGB_sum(const PSalgebra& p,
				  const latticetypes::Weight& lambda) const
{
  const rootdata::RootDatum& rd=rootDatum();
  kgb::KGBEltList sub=sub_KGB(p); std::reverse(sub.begin(),sub.end());

  std::vector<size_t> sub_inv(d_KGB.size(),~0);

  for (size_t i=0; i<sub.size(); ++i)
    sub_inv[sub[i]]=i; // partially fill array with inverse index

  std::vector<latticetypes::Weight> mu; // list of $\rho$-centered weights,
  mu.reserve(sub.size());               // associated to the elements of |sub|

  mu.push_back(lambda); (mu[0]-=rd.twoRho())/=2; // make $\rho$-centered

  for (size_t i=1; i<sub.size(); ++i)
  {
    kgb::KGBElt x=sub[i];
    bitset::RankFlags::iterator it;
    for (it=p.Levi_gens().begin(); it(); ++it)
    {
      if (d_KGB.cross(*it,x)>x) // then we can use ascending cross action
      {
	size_t k=sub_inv[d_KGB.cross(*it,x)];
	assert(k!=~0ul); // we ought to land in the subset
	mu.push_back(rd.simpleReflection(mu[k],*it)); // $\rho$-centered
	break;
      }
    }
    if (it()) continue; // if we could use a cross action, we're done for |i|

    // now similarly try Cayley transforms
    for (it=p.Levi_gens().begin(); it(); ++it)
    {
      if (d_KGB.cayley(*it,x)!=kgb::UndefKGB) // then we can use this Cayley
      {
	size_t k=sub_inv[d_KGB.cayley(*it,x)];
	assert(k!=~0ul); // we ought to land in the subset
	latticetypes::Weight nu=mu[k]; // $\rho-\lambda$ at split side
	assert(nu.scalarProduct(rd.simpleCoroot(*it))%2 == 0); // finality
	latticetypes::Weight alpha=rd.simpleRoot(*it);
	nu -= (alpha *= nu.scalarProduct(rd.simpleCoroot(*it))/2); // project
	mu.push_back(nu); // use projected weight at compact side of transform
	break;
      }
    }
    assert(it()); // if no cross action worked, some Cayley transform must have
  }
  assert(mu.size()==sub.size());

  size_t max_l=d_KGB.length(sub[0]);

  Raw_q_Char result;
  for (size_t i=0; i<sub.size(); ++i)
  {
    kgb::KGBElt x=sub[i];
    RawRep r(mu[i],d_KGB.titsElt(x));
    size_t codim=max_l-d_KGB.length(x);
    result += Raw_q_Char(r,q_CharCoeff(codim,codim%2==0 ? 1 : -1)); // $(-q)^c$
  }

  return result;
} // |q_KGB_sum|

// Express irreducible K-module as a finite virtual sum of standard ones
q_CharForm
KhatContext::q_K_type_formula(const StandardRepK& sr, level bound)
{
  const weyl::WeylGroup& W=weylGroup();
  const rootdata::RootDatum& rd=rootDatum();

  // Get theta stable parabolic subalgebra

  weyl::WeylWord conjugator;
  PSalgebra p = theta_stable_parabolic(sr,conjugator);

  latticetypes::Weight lambda=
    W.imageByInverse(rd,W.element(conjugator),lift(sr));

  Raw_q_Char q_KGB_sum_p= q_KGB_sum(p,lambda);

  // type of formal linear combination of weights, associated to Tits element

  q_Char result;
  for (Raw_q_Char::const_iterator
	 it=q_KGB_sum_p.begin(); it!=q_KGB_sum_p.end(); ++it)
  {
    q_CharCoeff c=it->second; // coefficient from |q_KGB_sum|
    const latticetypes::Weight& mu=it->first.first; // weight from |q_KGB_sum|
    const tits::TitsElt& strong=it->first.second; // Tits elt from |q_KGB_sum|
    cartanclass::InvolutionData id(complexGroup(),strong.tw());

    rootdata::RootSet A(rd.numRoots());
    for (bitmap::BitMap::iterator
	   rt=p.radical().begin(); rt!=p.radical().end(); ++rt)
    {
      rootdata::RootNbr alpha=*rt;
      assert(not id.real_roots().isMember(alpha));
      if (id.imaginary_roots().isMember(alpha))
	A.set_to(alpha,basedTitsGroup().grading(strong,alpha)); // add if nc
      else // complex root
      {
	rootdata::RootNbr beta=id.root_involution(alpha);
	assert(rd.isPosRoot(beta));
	A.set_to(alpha,beta>alpha); // add first of two complex roots
      }
    }

//     std::cout << "Sum over subsets of " << A.size() << " roots, giving ";

    typedef free_abelian::Monoid_Ring<latticetypes::Weight,q_CharCoeff>
      polynomial; // with weight exponents and $q$-polynomials as coefficients
    const latticetypes::LatticeMatrix theta =
      complexGroup().involutionMatrix(strong.tw());

    // compute $X^\mu*\prod_{\alpha\in A}(1-X^\alpha)$ in |pol|
    polynomial pol(mu);
    for (rootdata::RootSet::iterator it=A.begin(); it!=A.end(); ++it)
    {
      polynomial copy=pol; // since |add_multiple| assumes no aliasing
      pol.add_multiple(copy,q_CharCoeff(1,-1),rd.root(*it)); // $*(1-qX^\alpha)$

      // filter out terms that cannot affect anything below |bound|
      for (polynomial::iterator term=pol.begin(); term!=pol.end();)
      {
	latticetypes::Weight lambda=term->first;
	(lambda*=2) += rd.twoRho();
	lambda += theta.apply(lambda);
	if (height_bound(lambda)>bound)
	  pol.erase(term++);
	else
	  term++;
      }
    }
//     std::cout << pol.size() << " terms." << std::endl;

    // iterate over terms in formal sum, taking coef *= |c|
    for (polynomial::const_iterator term=pol.begin(); term!=pol.end(); ++term)
    {
      latticetypes::Weight lambda=term->first;
      polynomial::coef_t coef=term->second;
      result += q_Char(std_rep_rho_plus(lambda,strong),c*coef); // contribute
    }
  } // for sum over KGB for L
  return std::make_pair(sr, result);
} // q_K_type_formula

matrix::Matrix<CharCoeff> KhatContext::K_type_matrix
 (std::set<equation>& eq_set,
  level bound,
  std::vector<seq_no>& new_order)
{
  std::vector<equation> system=saturate(eq_set,bound);

  matrix::Matrix<CharCoeff>  m=triangularize(system,new_order);

#ifdef VERBOSE
  std::cout << "Ordering of representations/K-types:\n";
  for (std::vector<seq_no>::const_iterator
	 it=new_order.begin(); it!=new_order.end(); ++it)
    print(std::cout,rep_no(*it)) << ", height " << height(*it)
       << std::endl;

  prettyprint::printMatrix(std::cout<<"Triangular system:\n",m,3);
#endif

  matrix::Matrix<CharCoeff>m_inv=inverse_lower_triangular(m);

  return m_inv;

} // K_type_matrix


combination KhatContext::branch(seq_no s, level bound)
{
  combination result(height_graded); // a linear combination of $K$-types

  if (height(s)>bound)
    return result;

  // a linear combination of Final representations
  combination remainder(s,height_graded); // terms will have |height<=bound|
  do
  {
    combination::iterator head=remainder.begin(); // leading term

    equation eq=mu_equation(head->first,bound);

    result += combination(eq.first,head->second,height_graded);
    remainder.add_multiple(eq.second,-head->second);
  }
  while (not remainder.empty());

  return result;
}



/* convert a system of equations into a list, adding equations for all terms
   recursively (they are generated by |mu_equation|), up to the given
   bound (everything with |height(...)> bound| is pruned away).
   It is assumed that |mu_equation| will ensure all |seq_no|s in
   left and right hand sides are normalized, so that there is no risk of
   trying to add a formula for one term but getting one for another.

   Precondition: the right hand sides contain no terms with |height>bound|; if
   they are obtained from |mu_equation| with the same |bound|, this is assured.
 */
std::vector<equation>
KhatContext::saturate(const std::set<equation>& system, level bound)
{
  bitmap::BitMap lhs(nr_reps()); // left hand sides of all equations seen

  std::deque<equation> queue;

  for (std::set<equation>::iterator
	 it=system.begin(); it!=system.end(); ++it)
    if (height(it->first)<=bound)
    {
      queue.push_back(*it);
      lhs.insert(it->first); // include left hand sides from original system
    }

  std::vector<equation> result;

  while (not queue.empty())
  {
    // ensure bitmap provides space for all current terms
    if (nr_reps()>lhs.capacity())
      lhs.set_capacity( (nr_reps()+constants::posBits)
			& ~constants::posBits); // round up to wordsize

    const equation& cf=queue.front(); // an unprocessed formula
    assert(height(cf.first) <= bound);

    result.push_back(equation(cf.first,combination(height_graded)));
    combination& rhs=result.back().second;

    for (combination::const_iterator
	   term=cf.second.begin(); term!=cf.second.end(); ++term)
    {
      assert(height(term->first) <= bound); // guaranteed by |mu_equation|
      rhs.insert(*term);
      if (not lhs.isMember(term->first)) // no formula for this term seen yet
      {
	lhs.insert(term->first);
	queue.push_back(mu_equation(term->first,bound));
      }
    }

    queue.pop_front(); // we are done with this formula
  }

  return result;
} // saturate

// **************   manipulators **********************

void KhatContext::go(const StandardRepK& initial)
{
  combination chi=standardize(initial);

#ifdef VERBOSE
  if (nonfinal_pool.size()>0)
  {
    const rootdata::RootDatum& rd=rootDatum();
    std::cout << "Intermediate representations:\n";
    for (size_t i=0; i<nonfinal_pool.size(); ++i)
    {
      const StandardRepK& sr=nonfinal_pool[i];
      size_t witness;
      const cartanclass::Fiber& f=fiber(sr);
      print(std::cout << 'N' << i << ": ",sr) << " [" << height(sr) << ']';

      if (not isStandard(sr,witness))
	std::cout << ", non Standard, witness "
		  << rd.coroot(f.simpleImaginary(witness));
      if (isZero(sr,witness))
	std::cout << ", Zero, witness "
		  << rd.coroot(f.simpleImaginary(witness));
      if (not isFinal(sr,witness))
	std::cout << ", non Final, witness "
		  << rd.coroot(f.simpleReal(witness));
      std::cout << std::endl;
    }
  }
#endif

  std::cout << "Standard normal final limit representations:\n";
  for (seq_no i=0; i<nr_reps(); ++i)
  {
    const StandardRepK& sr=rep_no(i);
    print(std::cout << 'R' << i << ": ",sr) << " [" << height(i) << ']'
      << std::endl;
  }

  print(std::cout << "Standardized expression for ",initial) << ":\n";
  {
    std::ostringstream s; print(s,chi,true);
    ioutils::foldLine(std::cout,s.str(),"+\n- ","",1) << std::endl;
  }


} // go

/*****************************************************************************

        Chapter IV -- The PSalgebra class

******************************************************************************/

PSalgebra::PSalgebra (tits::TitsElt base,
		      const kgb::KGB& kgb,
		      const complexredgp::ComplexReductiveGroup& G)
    : strong_inv(base)
    , cn(G.class_number(base.tw()))
    , sub_diagram() // class |RankFlags| needs no dimensioning
    , nilpotents(G.rootDatum().numRoots())
{
  const rootdata::RootDatum& rd=G.rootDatum();
  cartanclass::InvolutionData id(G,base.tw());

  // Put real simple roots into Levi factor
  for (size_t i=0; i<rd.semisimpleRank(); ++i)
    if (id.real_roots().isMember(rd.simpleRootNbr(i)))
      sub_diagram.set(i);


  // put any imaginary or complex positive roots into radical
  for (size_t i=0; i<rd.numPosRoots(); ++i)
  {
    rootdata::RootNbr alpha=rd.posRootNbr(i);
    if (not id.real_roots().isMember(alpha))
      nilpotents.insert(alpha);
  }
}



// ****************** Chapter V -- functions ************************


q_combin to_q(const combination& c)
{
  q_combin result(c.key_comp()); // use same comparison object
  for (combination::const_iterator it=c.begin(); it!=c.end(); ++it)
    result += q_combin(it->first,q_CharCoeff(it->second),c.key_comp());
  return result;
}

matrix::Matrix<CharCoeff>
triangularize (const std::vector<equation>& system,
	       std::vector<seq_no>& new_order)
{
  // order set of equations
  std::vector<equation> equation(system.begin(),system.end());
  size_t n=equation.size();

  matrix::Matrix<CharCoeff> M(n,n,0);
  graph::OrientedGraph incidence(n);

  for (size_t j=0; j<n; ++j) // loop over equations
  {
    size_t n_terms=0;
    for (size_t i=0; i<n; ++i) // loop over left hand sides
      if ((M(i,j)=equation[j].second[equation[i].first])!=0)
      { // |OrientedGraph::cells| puts sinks in front, so record edge $i\to j$.
	incidence.edgeList(i).push_back(j);
	++n_terms;
      }

    if (equation[j].second.size()!=n_terms)
      throw std::runtime_error ("triangularize: system not saturated");
  }

  partition::Partition order; incidence.cells(order,NULL);

  new_order.resize(n);
  for (size_t i=0; i<n; ++i)
  {
    if (order.classSize(i)>1)
      throw std::runtime_error ("triangularize: system has cycles");
    new_order[order(i)]=equation[i].first;
  }

  matrix::Matrix<CharCoeff> result(n,n,0);
  for (size_t i=0; i<n; ++i)
    for (graph::EdgeList::const_iterator p=incidence.edgeList(i).begin();
	 p!=incidence.edgeList(i).end(); ++p) // there is an edge |i| to |*p|
      result(order(i),order(*p))=M(i,*p);     // so |order(i)>=order(*p)|

  return result;
} // triangularize

matrix::Matrix<CharCoeff> inverse_lower_triangular
  (const matrix::Matrix<CharCoeff>& L)
{
  size_t n=L.numColumns();
  if (L.numRows()!=n)
    throw std::runtime_error ("invert triangular: matrix is not square");

  matrix::Matrix<CharCoeff> result(n,n,0);

  for (size_t i=0; i<n; ++i)
  {
    if (L(i,i)!=1)
      throw std::runtime_error ("invert triangular: not unitriangular");
    result(i,i)=1;

    for (size_t j=i; j-->0; )
    {
      CharCoeff sum=0;
      for (size_t k=i; k>j; --k) // $j<k\leq i$
	sum += result(i,k)*L(k,j);
      result(i,j) = -sum;
    }
  }
  return result;
}

} // namespace standardrepk

// ****************** Chapter VI -- local functions ************************

namespace {

// orthogonal projection onto the intersection of kernels of coroots in |gens|
// the projection is parallel to the span of the roots in |gens|
latticetypes::LatticeMatrix
orth_projection(const rootdata::RootDatum& rd, bitset::RankFlags gens,
		latticetypes::LatticeCoeff& denom)
{
  size_t m=gens.count(), r=rd.rank();
  latticetypes::LatticeMatrix root_mat(r,m);
  latticetypes::LatticeMatrix sub_Cartan(m,m); // transposed, later inverted
  latticetypes::LatticeMatrix coroot_mat(m,r);
  for (bitset::RankFlags::iterator i=gens.begin(); i(); ++i)
  {
    size_t ii=gens.position(*i);
    for (bitset::RankFlags::iterator j=gens.begin(); j(); ++j)
      sub_Cartan(ii,gens.position(*j))=rd.cartan(*j,*i);
    for (size_t j=0; j<r; ++j)
    {
      root_mat(j,ii)=rd.simpleRoot(*i)[j];
      coroot_mat(ii,j)=rd.simpleCoroot(*i)[j];
    }
  }
  sub_Cartan.invert(denom); // invert and compute necessary denominator

  latticetypes::LatticeMatrix result(r,r,0); // set to identity scaled |denom|
  for (size_t i=0; i<r; ++i)
    result(i,i)=denom;
  result -= root_mat * sub_Cartan * coroot_mat;
  return result;
}

} // namespace


} // namespace atlas
