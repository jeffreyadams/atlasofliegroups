/*!
\file
\brief Implementation for the classes
StandardRepK and KHatComputations.
*/
/*
  This is standardrepk.cpp

  Copyright (C) 2006, 2007 Alfred Noel
  With comments (C) 2008 by Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups version 0.2.4

  See file main.cpp for full copyright notice
*/

#include "standardrepk.h"

#include "cartanclass.h"
#include "complexredgp.h"
#include "realredgp.h"
#include "kgb.h"
#include "descents.h"
#include "lattice.h"
#include "rootdata.h"
#include "smithnormal.h"
#include "intutils.h"
#include "cartanset.h"
#include "tori.h"
#include "basic_io.h"
#include "ioutils.h"
#include "intutils.h"
#include "tags.h"
#include "graph.h"
#include "prettyprint.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <deque>
#include <algorithm>

namespace atlas {



/*****************************************************************************

        Chapter I -- The StandardRepK class

******************************************************************************/

  namespace standardrepk {


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

const SR_rewrites::combination& SR_rewrites::lookup(seq_no n) const
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

        Chapter III -- The KHatComputations class

******************************************************************************/

namespace standardrepk {

  // for now we only construct a KHatComputations from a KGB structure

KHatComputations::KHatComputations
  (const realredgp::RealReductiveGroup &GR, const kgb::KGB& kgb)
  : d_G(&GR.complexGroup())
  , d_KGB(kgb)
  , nonfinal_pool(),final_pool()
  , nonfinals(nonfinal_pool), finals(final_pool)
  , d_rules()
  , simple_reflection_mod_2()
  , d_realForm(GR.realForm())
  , d_Tg(kgb::EnrichedTitsGroup::for_square_class(GR))
  , d_data(d_G->numCartanClasses())
{
  const rootdata::RootDatum& rd=rootDatum();
  simple_reflection_mod_2.reserve(d_G->semisimpleRank());
  for (size_t i=0; i<d_G->semisimpleRank(); ++i)
    simple_reflection_mod_2.push_back
      (latticetypes::BinaryMap
       (rd.rootReflection(rd.simpleRootNbr(i)).transposed()));

  size_t n = rootDatum().rank();

  // the following loop should be restricted to the Cartan classes for |GR|
  for (size_t r=0; r<d_data.size(); ++r)
  {
    // d_G->cartan[r] is rth CartanClass
    // which by now records the canonical involution for this class
    Cartan_info& ci=d_data[r];

    latticetypes::LatticeMatrix theta = d_G->cartan(r).involution();


    // put in $q$ the matrix of $\theta-1$
    latticetypes::LatticeMatrix q=theta;
    for (size_t j = 0; j < n; ++j)
      q(j,j) -= 1;

    // find Smith basis relative to $q$
    latticetypes::WeightList bs; matrix::initBasis(bs,n);
    latticetypes::CoeffList invf; smithnormal::smithNormal(invf,bs.begin(),q);

    latticetypes::LatticeMatrix basis(bs);
    latticetypes::LatticeMatrix inv_basis=basis.inverse();

//     prettyprint::printMatrix(std::cout,basis);
//     prettyprint::printMatrix(std::cout,inv_basis);


    size_t l = invf.size();
    size_t f=0; while (f<l and invf[f]==1) ++f; // ignore invariant factors 1

//     std::cout << "n=" << n << ", f=" << f << ", l=" << l << std::endl;

    // get columns [f,l) of |basis| to |torsionLift|, doubling the vectors,
    // and corresponding columns of |inv_basis| to |torsionProjector|, mod 2
    ci.torsionProjector.resize(l-f,n);
    for (size_t i=f; i<l; ++i)
    {
      assert(invf[i]==2); // other invariant factors must be 2

      ci.torsionLift.push_back(bs[i]);

      for (size_t j=0; j<n; ++j)
	ci.torsionProjector.set_mod2(i-f,j,inv_basis(i,j));
    }

    if (l==n) // avoid construction from empty list
      ci.freeLift.resize(n,0); // so set matrix to empty rectangle
    else
      ci.freeLift= // or take the final |n-l| basis vectors as columns
	latticetypes::LatticeMatrix(&bs[l],&bs[n],tags::IteratorTag());

    ci.freeProjector.resize(n-l,n);

    for (size_t i =l; i<n; ++i) // copy final rows from inverse
      ci.freeProjector.copyRow(inv_basis,i-l,i); // row |i| to |i-l|

    ci.fiber_modulus=latticetypes::SmallSubspace
      (latticetypes::SmallBitVectorList(tori::minusBasis(theta.transposed())),
       n);
  }
}

/******** accessors *******************************************************/

HCParam KHatComputations::project
  (size_t cn, latticetypes::Weight lambda) const
{
  const Cartan_info& ci=d_data[cn];

  (lambda -= rootDatum().twoRho()) /= 2; // final division is exact

  return std::make_pair
    (ci.freeProjector.apply(lambda),
     ci.torsionProjector.apply(latticetypes::SmallBitVector(lambda)).data()
     );
}

latticetypes::Weight KHatComputations::lift(size_t cn, HCParam p) const
{
  const Cartan_info& ci=d_data[cn];
  latticetypes::Weight result=ci.freeLift.apply(p.first); // lift free part

  latticetypes::WeightList torsion_lift=ci.torsionLift;
  for (size_t i=0; i<torsion_lift.size(); ++i)
    if (p.second[i])
      result += torsion_lift[i]; // add even vectors representing torsion part
  (result *= 2) += rootDatum().twoRho();

  return result;
}

StandardRepK KHatComputations::std_rep
  (const latticetypes::Weight& two_lambda, tits::TitsElt a) const
{
  weyl::TwistedInvolution sigma=a.tw();
  weyl::WeylElt w = d_G->cartanClasses().canonicalize(sigma);
  // now |sigma| is canonical and |w| conjugates |sigma| to |a.tw()|

  const weyl::WeylGroup& W=d_G->weylGroup();

  size_t cn = d_G->cartanClasses().classNumber(sigma);

  // conjugate towards canonical element
  {
    weyl::WeylWord ww=W.word(w);
    for (size_t i=0; i<ww.size(); ++i) // left-to-right for inverse conjugation
      d_Tg.basedTwistedConjugate(a,ww[i]);
  }
  assert(a.tw()==sigma); // we should now be at canonical twisted involution

  const rootdata::RootDatum& rd=rootDatum();
  latticetypes::Weight mu=W.imageByInverse(rd,w,two_lambda);
  StandardRepK result(cn,
		      d_data[cn].fiber_modulus.mod_image
		        (d_Tg.titsGroup().left_torus_part(a)),
		      project(cn,mu));

  const cartanclass::Fiber& cc=d_G->cartan(cn).fiber();
  for (size_t i=0; i<cc.imaginaryRank(); ++i)
    if (rd.scalarProduct(mu,cc.simpleImaginary(i))<0)
      return result; // without setting isStandard

  result.d_status.set(StandardRepK::IsStandard);
  return result;
} // std_rep


std::ostream& KHatComputations::print(std::ostream& strm,StandardRepK sr)
  const
{
  prettyprint::printVector(strm,lift(sr)) << '@';
  prettyprint::prettyPrint(strm,sr.d_fiberElt) << '#' << sr.d_cartan;
  return strm;
}

// map a character to one containing only Standard terms
// this version ensures the basic one is recursively called first
SR_rewrites::combination KHatComputations::standardize(const Char& chi)
{
  SR_rewrites::combination result;
  for (Char::base::const_iterator i=chi.begin(); i!=chi.end(); ++i)
    result.add_multiple(standardize(i->first),i->second);

  return result;
}

// the basic case
SR_rewrites::combination KHatComputations::standardize(StandardRepK sr)
{
  normalize(sr);

  { // first check if we've already done |sr|
    SR_rewrites::seq_no n=nonfinals.find(sr);
    if (n!=Hash::empty)
      return d_rules.lookup(n); // in this case an equation should be known
    if (sr.isStandard() // then also check if we already know |sr| to be Final
	and (n=finals.find(sr)!=Hash::empty))
      return SR_rewrites::combination(n); // single term known to be final
  }

  const cartanclass::Fiber& f=d_G->cartan(sr.d_cartan).fiber();
  const rootdata::RootDatum& rd=rootDatum();
  latticetypes::Weight lambda=lift(sr);

  if (sr.isStandard())
  {
    tits::TitsElt a=titsElt(sr);
    for (size_t i=0; i<f.imaginaryRank(); ++i)
      if (not d_Tg.grading(a,f.simpleImaginary(i))
          and rd.scalarProduct(lambda,f.simpleImaginary(i))==0)
      {
	SR_rewrites::combination zero;
	d_rules.equate(nonfinals.match(sr),zero);
	return zero;
      }
    sr.d_status.set(StandardRepK::IsNonZero);

    // find any simple-real root for which |sr| fails to be Final
    // since coordinates are doubled, the scalar product is 0 or 2 mod 4
    rootdata::RootNbr alpha;
    {
      size_t i;
      for (i=0; i<f.realRank(); ++i)
	if (rd.scalarProduct(lambda,f.simpleReal(i))%4 == 0)
	  break;
      if (i==f.realRank())
      {
	sr.d_status.set(StandardRepK::IsFinal);
	return SR_rewrites::combination(finals.match(sr));
      }
      alpha=f.simpleReal(i);
    }

    HechtSchmid equation= back_HS_id(sr,alpha);
    assert(equation.lh2==NULL); // |back_HS_id| never produces a second member

#ifdef VERBOSE
    print(std::cout << "RHS: ",sr)<<" = ";
#endif

    assert(equation.rh1!=NULL);

    Char rhs(*equation.rh1); // rhs stars as second member negated
#ifdef VERBOSE
    print(std::cout,*equation.rh1);
#endif
    if (equation.rh2!=NULL)
    {
#ifdef VERBOSE
      print(std::cout<<'+',*equation.rh2);
#endif
      rhs+= Char(*equation.rh2);
    }
#ifdef VERBOSE
    std::cout << std::endl;
#endif

    // now recursively standardize all terms, storing rules
    SR_rewrites::combination result= standardize(rhs);
    d_rules.equate(nonfinals.match(sr),result); // and add rule for |sr|
    return result;
  } // if (sr.isStandard())

  // find (again) simple-imaginary root for which |sr| fails to be Standard
  rootdata::RootNbr alpha;
  {
    size_t i;
    for (i=0; i<f.imaginaryRank(); ++i)
      if (rd.scalarProduct(lambda,f.simpleImaginary(i))<0)
	break;
    assert(i<f.imaginaryRank()); // otherwise |sr| would have been Standard
    alpha=f.simpleImaginary(i);
  }

  HechtSchmid equation= HS_id(sr,alpha);
  assert(equation.lh2!=NULL); // all cases of |HS_id| produce a second member

#ifdef VERBOSE
  print(print(std::cout << "HS: ",sr)<<'+',*(equation.lh2)) << " = ";
#endif

  Char rhs(*equation.lh2,-1); // rhs stars as second member negated

  if (equation.rh1!=NULL)
  {
#ifdef VERBOSE
    print(std::cout,*equation.rh1);
#endif
    rhs+= Char(*equation.rh1);
    if (equation.rh2!=NULL)
    {
#ifdef VERBOSE
      print(std::cout<<'+',*equation.rh2);
#endif
      rhs+= Char(*equation.rh2);
    }
  }
#ifdef VERBOSE
  else std::cout << '0';
  std::cout << std::endl;
#endif

  // now recursively standardize all terms, storing rules
  SR_rewrites::combination result= standardize(rhs);
  d_rules.equate(nonfinals.match(sr),result); // and add rule for |sr|
  return result;
} // standardize

void KHatComputations::normalize(StandardRepK& sr) const
{
  const rootdata::RootDatum& rd = rootDatum();
  const cartanclass::Fiber& f=d_G->cartan(sr.d_cartan).fiber();

  latticetypes::Weight lambda=lift(sr);

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
  atlas::rootdata::RootSet cplx = f.complexRootSet();

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


// put the new weakly dominant lambda back in sr

  sr.d_lambda=project(sr.d_cartan,lambda);

} // normalize


HechtSchmid
KHatComputations::HS_id(StandardRepK sr, rootdata::RootNbr alpha) const
{
  HechtSchmid id(sr);
  const rootdata::RootDatum& rd=rootDatum();
  tits::TitsElt a=titsElt(sr);
  latticetypes::Weight lambda=lift(sr);
  assert(rd.isPosRoot(alpha)); // in fact it must be simple-imaginary

  size_t i=0; // simple root index (value will be set in following loop)
  while (true) // we shall exit halfway when $\alpha=\alpha_i$
  {
    while (rd.scalarProduct(alpha,rd.simpleRootNbr(i))<=0)
    {
      ++i;
      assert(i<rd.semisimpleRank());
    }
    // now $\<\alpha,\alpha_i^\vee> > 0$ where $\alpha$ is simple-imaginary
    // and \alpha_i$ is complex for the involution |a.tw()|
    if (alpha==rd.simpleRootNbr(i)) break; // found it
    // otherwise reflect all data by $s_i$, which decreases level of $\alpha$
    rd.rootReflect(alpha,i);
    rd.simpleReflect(lambda,i);
    d_Tg.basedTwistedConjugate(a,i);
    i=0; // and start over
  }

  latticetypes::Weight mu=rd.simpleReflection(lambda,i);
  if (d_Tg.simple_grading(a,i))
  { // |alpha| is a non-compact root
    d_Tg.basedTwistedConjugate(a,i); // adds $m_i$ to torus part
    id.add_lh(std_rep(mu, a));
    assert(sr.d_cartan==id.lh2->d_cartan);
    // the change to |d_lambda| may involve both components

    /* Now are equivalent:
       = |id.lh2->d_fiberElt==sr.d_fiberElt|
       = $m_i$ absorbed into quotient (i.e., $m_i\in(X_*^- mod 2)$)
       = $\alpha^\vee( (X^*)^\theta ) = 2\Z$
       = $\alpha\notin(1-\theta_\alpha)X^*$  ($\theta_\alpha=\theta.s_\alpha$)
       = Cayley transform is type II
    */

    d_Tg.Cayley_transform(a,i);
    id.add_rh(std_rep(lambda, a));
    if (id.lh2->d_fiberElt==sr.d_fiberElt) // type II
    {
      lambda += rootDatum().simpleRoot(i); // other possibility
      lambda += rootDatum().simpleRoot(i); // add twice: doubled coordinates
      id.add_rh(std_rep(lambda, a));
      assert(id.rh1->d_lambda != id.rh2->d_lambda);
    }
  }
  else // $\alpha_i$ is a compact root; "easy" Hecht-Schmid identity
    // based twisted conjugation fixed the Tits element; just reflect weight
    id.add_lh(std_rep(mu, a)); // and no RHS

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
KHatComputations::back_HS_id(StandardRepK sr, rootdata::RootNbr alpha) const
{
  HechtSchmid id(sr);
  const rootdata::RootDatum& rd=rootDatum();
  tits::TitsElt a=titsElt(sr);
  latticetypes::Weight lambda=lift(sr);

  // basis used is of $(1/2)X^*$, so scalar product with coroot always even
  assert(rd.scalarProduct(lambda,alpha)%4 == 0); // the non-final condition

  {
    latticetypes::Weight mu=rd.root(alpha);
    mu *= rd.scalarProduct(lambda,alpha)/2; // an integer multiple of $\alpha$
    lambda -= mu; // this makes lambda invariant under $s_\alpha$
    /* in type I, $\alpha$ is in $(1-\theta)X^*$ and correction is neutral
       in type II, correction need not be in $(1-\theta)X^*$, but adding
       $\alpha$ gives HC parameter designating the same representation
    */
  }

  latticetypes::SmallSubspace mod_space=
    d_data[sr.d_cartan].fiber_modulus; // make a copy to be modified

  // again, move to situation where $\alpha$ is simple: $\alpha=\alpha_i$
  size_t i=0; // simple root index (value will be set in following loop)

  while (true) // we shall exit halfway when $\alpha=\alpha_i$
  {
    while (rd.scalarProduct(alpha,rd.simpleRootNbr(i))<=0)
    {
      ++i;
      assert(i<rd.semisimpleRank());
    }
    // now $\<\alpha,\alpha_i^\vee> > 0$ where $\alpha$ is simple-real
    // and \alpha_i$ is complex for the involution |a.tw()|
    if (alpha==rd.simpleRootNbr(i)) break; // found it
    // otherwise reflect all data by $s_i$, which decreases level of $\alpha$
    rd.rootReflect(alpha,i);
    rd.simpleReflect(lambda,i);
    d_Tg.basedTwistedConjugate(a,i);
    mod_space.apply(simple_reflection_mod_2[i]);
    i=0; // and start over
  }

  // one right term is obtained by undoing Cayley for |a|, with lifted |lambda|
  d_Tg.inverse_Cayley_transform(a,i,mod_space);

  id.add_rh(std_rep(lambda,a));

  // there will be another term in case of a type I HechtSchmid identity
  // it differs only in the fiber part, by $m_i$; if this vanishes into the
  // quotient, the HechtSchmid identity is type II and nothing is added
  d_Tg.basedTwistedConjugate(a,i);
  StandardRepK other=std_rep(lambda,a);
  if (other.d_fiberElt!=id.rh1->d_fiberElt) // type I
    id.add_rh(other);

  return id;
}

atlas::latticetypes::LatticeCoeff
KHatComputations::product_simpleroot
  (const StandardRepK& s, rootdata::RootNbr k) const
{
  latticetypes::Weight lifted=theta_lift(s.d_cartan,s.d_lambda);

  return rootDatum().scalarProduct(lifted,k);
}

atlas::latticetypes::LatticeCoeff
KHatComputations::product_sumposroots(const StandardRepK& s) const
{
  latticetypes::Weight lifted=theta_lift(s.d_cartan,s.d_lambda);

  latticetypes::Weight dual2rho=rootDatum().dual_twoRho();

  return intutils::abs(lifted.scalarProduct(dual2rho));
}

PSalgebra
KHatComputations::theta_stable_parabolic
  (weyl::WeylWord& conjugator,
   const cartanset::CartanClassSet& cs,
   const size_t Cartan_nr) const
{
  const rootdata::RootDatum& rd=cs.rootDatum();
  const weyl::WeylGroup& W=cs.weylGroup();
  weyl::TwistedInvolution twi=cs.twistedInvolution(Cartan_nr);
  // latticetypes::LatticeMatrix theta=cs.involutionMatrix(twi);


/* Conjugate |twi| by simple complex roots to make the positive root system
   more theta stable. There is no hope for real roots, but for complex roots
   we try to achieve that the $\theta$-image of simple roots is positive. As
   our notion of positive roots is fixed, we conjugate $\theta$ itself (in
   fact twisted-conjugating |twi|) which changes the notions of
   imaginary/real/complex roots.
*/
  conjugator.resize(0);

  {
    size_t i;
    do // try to twisted-conjugate by complex simple root mapped to negative
    {
      for (i=0; i<rd.semisimpleRank(); ++i)
      {
	rootdata::RootNbr alpha=rd.simpleRootNbr(i);
	rootdata::RootNbr beta= // image of |alpha| by |twi|
	  rd.permuted_root(W.word(twi.w()),rd.simpleRootNbr(W.twisted(i)));
	if (not rd.isPosRoot(beta) and beta!=rd.rootMinus(alpha))
	{
	  W.twistedConjugate(twi,i);
	  conjugator.push_back(i);
	  break; // and repeat do-while loop
	}
      } // for i
    }
    while (i!=rd.semisimpleRank());
  }
/*
   We have achieved that any real positive root is a sum of real simple roots.
   Here's why. Consider a counterexample that is minimal: a real positive root
   $\alpha$ from which no real simple root can be subtracted while remaining
   positive. Then this must be a sum of simple roots that are either complex
   or imaginary; the $\theta$-images of such simple roots are all positive,
   but their sum is $-\alpha$ which is negative, a contradiction.

   Meanwhile |conjugator| is such that |W.twistedConjugate(twi,conjugator)|
   would make |twi==cs.twistedInvolution(Cartan_nr)| again.
*/

  // Build the parabolic subalgebra:

  return  PSalgebra(cs,twi);

} // theta_stable_parabolic

kgb::KGBEltList KHatComputations::sub_KGB(const PSalgebra& q) const
{
  bitmap::BitMap flagged(d_KGB.size());
  kgb::KGBElt root=d_KGB.tauPacket(q.theta()).first; // choose a torus part
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

free_abelian::Free_Abelian<StandardRepK>
KHatComputations::KGB_sum(const PSalgebra& q,
			  const latticetypes::Weight& lambda) const
{
  const rootdata::RootDatum& rd=rootDatum();
  kgb::KGBEltList sub=sub_KGB(q); std::reverse(sub.begin(),sub.end());

  std::vector<size_t> sub_inv(d_KGB.size(),~0);

  for (size_t i=0; i<sub.size(); ++i)
    sub_inv[sub[i]]=i; // partially fill array with inverse index

  std::vector<latticetypes::Weight> mu; mu.reserve(sub.size());

  mu.push_back(lambda); (mu[0]-=rd.twoRho())/=2; // rho-centered weight

  for (size_t i=1; i<sub.size(); ++i)
  {
    kgb::KGBElt x=sub[i];
    bitset::RankFlags::iterator it;
    for (it=q.Levi_gens().begin(); it(); ++it)
    {
      if (d_KGB.cross(*it,x)>x) // then we can use ascending cross action
      {
	assert(sub_inv[d_KGB.cross(*it,x)]!=~0ul); // we land in the subset
	latticetypes::Weight nu=mu[sub_inv[d_KGB.cross(*it,x)]];
	mu.push_back(rd.simpleReflection(nu,*it)); // rho-centered reflection
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
	assert(k!=~0ul); // we land in the subset
	latticetypes::Weight nu=mu[k];
	assert(nu.scalarProduct(rd.simpleCoroot(*it))%2 == 0); // finality
	latticetypes::Weight alpha=rd.simpleRoot(*it);
	nu -= (alpha *= nu.scalarProduct(rd.simpleCoroot(*it))/2); // project
	mu.push_back(nu); break;
      }
    }
    assert(it()); // if no cross action worked, some Cayley transform must
  }
  assert(mu.size()==sub.size());

  size_t max_l=d_KGB.length(sub[0]);

  free_abelian::Free_Abelian<StandardRepK> result;
  for (size_t i=0; i<sub.size(); ++i)
  {
    kgb::KGBElt x=sub[i];
    standardrepk::StandardRepK sr=std_rep_rho_plus(mu[i],d_KGB.titsElt(x));
    result += free_abelian::Free_Abelian<StandardRepK>
               (sr, ((max_l-d_KGB.length(x))%2 == 0 ? 1 : -1));
  }

  return result;
} // KGB_sum

// Express irreducible K-module as a finite virtual sum of standard ones
CharForm  KHatComputations::character_formula(StandardRepK sr) const
{
  const cartanset::CartanClassSet& cs=d_G->cartanClasses();
  const weyl::WeylGroup& W=cs.weylGroup();
  const rootdata::RootDatum& rd=cs.rootDatum();

  // Get theta stable parabolic subalgebra

  weyl::WeylWord conjugator;
  PSalgebra q = theta_stable_parabolic(conjugator,cs,sr.d_cartan);
  weyl::WeylElt w=W.element(conjugator);

  latticetypes::Weight lambda=W.imageByInverse(rd,w,lift(sr));

  Char multmap(sr);//this handles the empty subset

  latticetypes::LatticeMatrix P=d_data[sr.d_cartan].freeProjector;

  size_t nsubset;

  for (unsigned long i=1; i<nsubset; ++i) // bits of |i| determine the subset
  {
    bitset::RankFlags subset(i);

    latticetypes::Weight lambda = sr.d_lambda.first;
//    for (bitset::RankFlags::iterator j =subset.begin(); j(); j++)
//       if (rpairlist[*j].first == rpairlist[*j].second) // imaginary root
//       { // nothing implemented yet for imaginary roots
//       }
//       else // complex pair of roots
//       {
// 	latticetypes::Weight mu=
// 	  P.apply(cs.rootDatum().root(rpairlist[*j].first));
// 	lambda +=mu; // add the restriction of the first root to lambda
//       }


    StandardRepK new_rep = sr;

    new_rep.d_lambda.first = lambda; // replace |lambda| by modified version

    normalize(new_rep);

    // now add $(-1)^{\#S}$ to the coefficient of |new_rep| in |multmap|
    long sign=subset.count()%2 == 0 ? 1 : -1;
    std::pair<Char::iterator,bool> p=
      multmap.insert(std::make_pair(new_rep,sign));
    if (not p.second) p.first->second += sign;

  } // for (subsets)

  // finally remove items with multiplicity 0

  for (Char::iterator iter = multmap.begin();iter !=multmap.end();)
    if (iter->second == 0)
      multmap.erase(iter++); // must take care to do ++ before erase
    else
      ++iter;

  return std::make_pair(sr, multmap);
} // character_formula

atlas::matrix::Matrix<CharCoeff>
KHatComputations::makeMULTmatrix
 (std::set<CharForm>& column,
  const atlas::latticetypes::LatticeCoeff bound) const
{
  rootdata::RootDatum rd = d_G->rootDatum();
  latticetypes::Weight tr=rd.twoRho();
  latticetypes::Weight cotr = rd.dual_twoRho();

  std::set <atlas::latticetypes::Weight> lookup;

  //it is assumed that the characters in column.first are different

  for (std::set<CharForm>::iterator i = column.begin(); i != column.end();)
    if (product_sumposroots(i->first) > bound)
      column.erase(i++); // increment iterator before actual erase
    else
    {
      lookup.insert(i->first.d_lambda.first);
      ++i;
    }

  assert (column.size()==lookup.size());


  std::vector<CharForm> system=saturate(column,bound);

  std::vector<StandardRepK> new_order;
  atlas::matrix::Matrix<CharCoeff>  m=triangularize(system,new_order);
  for (std::vector<StandardRepK>::const_iterator
	 it=new_order.begin(); it!=new_order.end(); ++it)
  {
    const latticetypes::Weight& l=it->d_lambda.first;
    basic_io::seqPrint(std::cout,l.begin(),l.end(), ", ", "[", "]\n");
  }

  prettyprint::printMatrix(std::cout<<"Triangular system:\n",m,3);

  matrix::Matrix<CharCoeff>m_inv=inverse_upper_triangular(m);

  prettyprint::printMatrix(std::cout<<"Inverse matrix:\n",m_inv,3);

  return m_inv;

} // makeMULTmatrix




/* convert a system of equations into a list, adding equations for all terms
   recursively (they are generated by |character_formula|), up to the given
   bound (everything with |product_sumposroots(...)> bound| is pruned away).
   It is assumed that |character_formula| will ensure all |StandardRepK|s in
   left and right hand sides are normalized, so that there is no risk of
   trying to add a formula for one term but getting one for another.
 */
std::vector<CharForm>
KHatComputations::saturate(std::set<CharForm> system,
			   atlas::latticetypes::LatticeCoeff bound) const
{
  std::set<StandardRepK> lhs; // left hand sides of all formulae seen so far
  for (std::set<CharForm>::iterator it=system.begin(); it!=system.end(); ++it)
    lhs.insert(it->first);

  std::vector<CharForm> result;

  while (not system.empty())
  {
    std::set<CharForm>::iterator current=system.begin(); // choose one
    const CharForm& cf=*current;           // this is an unprocessed formula
    if (product_sumposroots(cf.first) <= bound) // ignore if out of bounds
    {
      result.push_back(CharForm(cf.first,Char())); // start with empty rhs
      Char& rhs=result.back().second;
      for (Char::const_iterator
	     term=cf.second.begin(); term!=cf.second.end(); ++term)
	if (product_sumposroots(term->first) <= bound)
	{
	  rhs.insert(*term);
	  if (lhs.count(term->first)==0) // no formula for this term seen yet
	  {
	    lhs.insert(term->first);
	    system.insert(character_formula(term->first));
	  }
	}

    }
    system.erase(current); // we are done with this formula
  }

  return result;
} // saturate



// **************   manipulators **********************

void KHatComputations::go
  (kgb::KGBElt x, const latticetypes::Weight& lambda)
{
  StandardRepK sr=std_rep_rho_plus(lambda,d_KGB.titsElt(x));

  {
    latticetypes::Weight mu=lambda;
    (mu *= 2) += rootDatum().twoRho();
    prettyprint::printVector(std::cout << "Weight entered: (1/2)",mu);
    prettyprint::printVector(std::cout << " converted to: (1/2)",lift(sr))
      << std::endl;
  }

  SR_rewrites::combination chi=standardize(sr);

#ifdef VERBOSE
  if (nonfinal_pool.size()>0)
  {
    std::cout << "Non-final representations:\n";
    for (size_t i=0; i<nonfinal_pool.size(); ++i)
    {
      const StandardRepK& sr=nonfinal_pool[i];
      assert(not sr.isFinal());
      print(std::cout << 'N' << i << ": ",sr)
		      << ( sr.isNonZero() ? " Non Final."
			   : sr.isStandard() ? " Zero."
			   : " Non Standard."
			   )
		      << std::endl;
    }
  }
#endif

  std::cout << "Standard normal final limit representations:\n";
  for (SR_rewrites::seq_no i=0; i<nr_reps(); ++i)
  {
    const StandardRepK& sr=rep_no(i);
    assert(sr.isFinal());
    print(std::cout << 'R' << i << ": ",sr) << std::endl;
  }

  std::cout << "Standardized expression:\n";
  {
    std::ostringstream s;
    for (SR_rewrites::combination::const_iterator it=chi.begin();
	 it!=chi.end(); ++it)
    {
      s << (it->second>0 ? it==chi.begin() ? "" : " + " : " - ");
      long int ac=intutils::abs<long int>(it->second);
      if (ac!=1)
	s << ac << '*';
      s << 'R' << it->first;
    }
    ioutils::foldLine(std::cout,s.str(),"+-","",1) << std::endl;
  }


} // go

/*****************************************************************************

        Chapter IV -- The PSalgebra class

******************************************************************************/

PSalgebra::PSalgebra (const cartanset::CartanClassSet& cs,
		      const weyl::TwistedInvolution& theta)
    : twi(theta), cartan(cs.classNumber(twi))
    , L(), nilpotents(cs.rootDatum().numRoots())
{
  const rootdata::RootDatum& rd=cs.rootDatum();
  cartanclass::InvolutionData id(rd,cs.involutionMatrix(twi));

  // Put real simple roots into Levi factor
  for (size_t i=0; i<rd.semisimpleRank(); ++i)
    if (id.real_roots().isMember(rd.simpleRootNbr(i)))
      L.set(i);


  // put any imaginary or complex positive roots into radical
  for (rootdata::RootSet::iterator it = rd.posRootSet().begin(); it(); ++it)
    if (not id.real_roots().isMember(*it))
      nilpotents.insert(*it);
}



// ****************** Chapter V -- functions ************************

atlas::matrix::Matrix<CharCoeff>
triangularize (const std::vector<CharForm>& system,
	       std::vector<StandardRepK>& new_order)
{
  std::vector<CharForm> equation(system.begin(),system.end()); // numbering
  size_t n=equation.size();

  atlas::matrix::Matrix<CharCoeff> M(n,n,0);
  graph::OrientedGraph usage(n);
  for (size_t i=0; i<n; ++i) // loop over equations
  {
    size_t n_terms=0;
    for (size_t j=0; j<n; ++j) // loop over left hand sides
    {
      Char::base::const_iterator p= equation[i].second.find(equation[j].first);
      if (p!=equation[i].second.end())
      { // |OrientedGraph::cells| put sinks in front, so record edge $j\to i$.
	usage.edgeList(j).push_back(i); M(i,j)=p->second;
	++n_terms;
      }
    }
    if (equation[i].second.size()!=n_terms)
      throw std::runtime_error ("triangularize: system not saturated");
  }

  partition::Partition order; usage.cells(order,NULL);

  new_order.resize(n);
  for (size_t i=0; i<n; ++i)
  {
    if (order.classSize(i)>1)
      throw std::runtime_error ("triangularize: system has cycles");
    new_order[order(i)]=equation[i].first;
  }

  atlas::matrix::Matrix<CharCoeff> result(n,n,0);
  for (size_t j=0; j<n; ++j)
    for (graph::EdgeList::const_iterator
	   p=usage.edgeList(j).begin(); p != usage.edgeList(j).end(); ++p)
      result(order(*p),order(j))=M(*p,j);
  return result;
} // triangularize

matrix::Matrix<CharCoeff> inverse_upper_triangular
  (const atlas::matrix::Matrix<CharCoeff>& U)
{
  size_t n=U.numColumns();
  if (U.numRows()!=n)
    throw std::runtime_error ("invert triangular: matrix is not square");

  matrix::Matrix<CharCoeff> result(n,n,0);

  for (size_t j=0; j<n; ++j)
  {
    if (U(j,j)!=1)
      throw std::runtime_error ("invert triangular: not unitriangular");
    result(j,j)=1;

    for (size_t i=j; i-->0; )
    {
      CharCoeff sum=0;
      for (size_t k=j; k>i; --k)
	sum += U(i,k)*result(k,j);
      result(i,j) = -sum;
    }
  }
  return result;
}




} // namespace standardrepk
} // namespace atlas
