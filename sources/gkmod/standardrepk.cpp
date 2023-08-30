/*
  Implementation for the classes |StandardRepK| and |KhatContext|.

  This is standardrepk.cpp

  Copyright (C) 2006, 2007 Alfred Noel
  Copyright (C) 2008, 2009 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "standardrepk.h"

#include <cstdlib> // for |std::abs|
#include <cassert>
#include <iostream>
#include <sstream>

#include "tags.h"
#include "matreduc.h"	// |adapted_basis|
#include "polynomials.h"// for $q$-coefficents

#include "lattice.h"	// |eigen_lattice|
#include "rootdata.h"	// various methods
#include "kgb.h"	// associate information with |KGBElt| values
#include "graph.h"	// |OrientedGraph| used in |triangularize|
#include "ioutils.h"	// |foldLine| in |KhatContext::go|
#include "basic_io.h"   // |operator<<| in |SRK_context::print|
#include "prettyprint.h"// |printVector| in |SRK_context::print|

namespace atlas {

/*
     Chapter 0 -- Local function declarations
*/

namespace {

// an auxiliary function for height computations
int_Matrix
orth_projection(const RootDatum& rd, RankFlags gens,
		LatticeCoeff& denom);

} // |namespace|

/*****************************************************************************

        Chapter I -- The HechtSchmid and StandardRepK classes

******************************************************************************/

namespace standardrepk {

Char HechtSchmid::rhs () const
{
  if (rh1==nullptr)
    return Char();

  Char result(*rh1);
  if (rh2!=nullptr)
    result+= Char(*rh2);

  return result;
}

Char HechtSchmid::equivalent () const
{
  Char result=rhs();
  if (lh2!=nullptr)
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
    hash+=(hash<<3)+(hash<<12)+d_lambda.first[i];
  return (hash+(hash<<5)+d_lambda.second.to_ulong())&(modulus-1);
}


/*****************************************************************************

        Chapter II -- The SRK_context class

******************************************************************************/

RootNbr SRK_context::offender = -1; // meaning no error signalled yet

SRK_context::SRK_context(RealReductiveGroup &GR)
  : G(GR)
  , Cartan_set(GR.Cartan_set())
  , C_info(GR.numCartan())
  , simple_reflection_mod_2()
  , proj_pool(), proj_sets(proj_pool), proj_data()
{
  const RootDatum& rd=root_datum();
  simple_reflection_mod_2.reserve(G.semisimple_rank());
  for (size_t i=0; i<G.semisimple_rank(); ++i)
    simple_reflection_mod_2.push_back
      (BinaryMap(rd.simple_reflection(i).transposed()));

  size_t n = root_datum().rank();

  // what remains is filling |C_info|, which is quite a bit of code
  size_t nr=0;
  for (BitMap::iterator it=Cartan_set.begin(); it(); ++it,++nr)
  {

    // d_G.cartan[i] is (canonical involution for) (*it)th CartanClass
    Cartan_info& ci=C_info[nr];

    const Fiber& f=G.cartan(*it).fiber();
    const WeightInvolution& theta = f.involution();

    const WeightInvolution q=theta-1;

    // find basis adapted to image of $\theta-1$
    CoeffList factor;
    int_Matrix basis = matreduc::adapted_basis(q,factor);
    int_Matrix inv_basis=basis.inverse();

    size_t twos=0,l=factor.size();
    RankFlags torsion;

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
    for (RankFlags::iterator it=torsion.begin(); it(); ++it,++i)
    {
      ci.torsionLift.push_back(basis.column(*it));

      for (size_t j=0; j<n; ++j)
	ci.torsionProjector.set_mod2(i,j,inv_basis(*it,j));
    }

    basis.block(0,l,n,n).swap(ci.freeLift); // final |n-l| basis vectors
    inv_basis.block(l,0,n,n).swap(ci.freeProjector); // final |n-l| rows

    SmallSubspace
      (BinaryMap(lattice::eigen_lattice(theta.transposed(),-1)))
      .swap(ci.fiber_modulus);

    { // find simple roots orthogonal to |real2rho| and |imaginary2rho|
      Weight real2rho=rd.twoRho(f.realRootSet());
      Weight imaginary2rho=rd.twoRho(f.imaginaryRootSet());
      for (size_t i=0; i<rd.semisimple_rank(); ++i)
	if (rd.is_orthogonal(real2rho,rd.simpleRootNbr(i)) and
	    rd.is_orthogonal(imaginary2rho,rd.simpleRootNbr(i)))
	{ // test coroot orthogonality
	  RootNbr alpha = rd.simpleRootNbr(i);
	  RootNbr beta= f.involution_image_of_root(alpha);
	  assert (rd.is_simple_root(beta));
	  if (not ci.bi_ortho[rd.simpleRootIndex(beta)]) // skip second of pair
	  {
	    ci.bi_ortho.set(i);
	    ci.sum_coroots.push_back(rd.coroot(alpha)+rd.coroot(beta));
	  }
	}
    } // filling |bi_orth| and |sum_coroots| fields

  } // |for (it)|
} // |SRK_context::SRK_context|


HCParam SRK_context::project (CartanNbr cn, Weight lambda) const
{
  const Cartan_info& ci=info(cn);

  (lambda -= root_datum().twoRho()) /= 2; // final division is exact

  // now |lambda| actually represents $\lambda-\rho$ in plain coordinates
  return std::make_pair
    (ci.freeProjector*lambda,
     (ci.torsionProjector*SmallBitVector(lambda)).data()
     );
}

Weight SRK_context::lift(CartanNbr cn, HCParam p) const
{
  const Cartan_info& ci=info(cn);
  Weight result=ci.freeLift*p.first; // lift free part

  const WeightList& torsion_lift = ci.torsionLift;
  for (size_t i=0; i<torsion_lift.size(); ++i)
    if (p.second[i])
      result += torsion_lift[i]; // add even vectors representing torsion part
  (result*=2) += root_datum().twoRho(); // convert $\lambda-\rho$ to $2\lambda$

  return result;
}

StandardRepK SRK_context::std_rep (const Weight& two_lambda, TitsElt a) const
{
  const RootDatum& rd=root_datum();

  TwistedInvolution sigma=a.tw();
  const WeylWord ww = innerClass().canonicalize(sigma);
  // now |sigma| is canonical and |w| conjugates |sigma| to |a.tw()|


  // conjugate towards canonical element (left-right, for inverse conjugation)
  basedTitsGroup().basedTwistedConjugate(a,ww);

  assert(a.tw()==sigma); // we should now be at canonical twisted involution

  Weight mu=rd.image_by_inverse(two_lambda,ww); // move weight to canonical

  CartanNbr cn = innerClass().class_number(sigma);
  StandardRepK result(cn,
		      info(cn).fiber_modulus.mod_image
		        (Tits_group().left_torus_part(a)),
		      project(cn,mu));

  return result;
} // |std_rep|


#if 0
// the following is a variant of |std_rep_rho_plus| intended for |KGB_sum|
// it should only transform the parameters for the Levi factor given by |gens|
// since |lambda| is $\rho$-centered, care should be taken in transforming it
RawRep SRK_context::Levi_rep (Weight lambda, TitsElt a, RankFlags gens) const
{
  TwistedInvolution sigma=a.tw();
  const WeylWord ww = innerClass().canonicalize(sigma,gens);
  // now |sigma| is canonical for |gens|, and |w| conjugates it to |a.tw()|

  const RootDatum& rd=rootDatum();

  // conjugate towards canonical element
  {
    for (size_t i=0; i<ww.size(); ++i) // left-to-right for inverse conjugation
    {
      assert(gens.test(ww[i])); // check that we only used elements in $W(L)$
      basedTitsGroup().basedTwistedConjugate(a,ww[i]);
      rd.simple_reflect(ww[i],lambda,1); // affine reflection fixing $-\rho$
    }
  }
  assert(a.tw()==sigma); // we should now be at canonical twisted involution

  return RawRep (lambda,a);
} // |Levi_rep|

StandardRepK SRK_context::KGB_elt_rep(KGBElt z) const
{ return std_rep(rootDatum().twoRho(),kgb().titsElt(z)); }
#endif


level SRK_context::height(const StandardRepK& sr) const
{
  const RootDatum& rd=root_datum();
  const Weight mu=theta_lift(sr);

  level sum=0;
  for (auto it=rd.beginPosCoroot(); it!=rd.endPosCoroot(); ++it)
    sum += std::abs(mu.dot(*it));

  return sum/4; // |mu| has doubled coordinates, and use the coroot half-sum
} // |SRK_context::height|


const proj_info& SRK_context::get_projection(RankFlags gens)
{
  size_t old_size=proj_data.size();
  size_t h=proj_sets.match(bitset_entry(gens));
  if (h<old_size)
    return proj_data[h];

  proj_data.push_back(proj_info());
  proj_data.back().projection=
    orth_projection(root_datum(),gens,proj_data.back().denom);
  return proj_data.back();
} // |SRK_context::get_projection|

level SRK_context::height_bound(const Weight& lambda)
{
  const RootDatum& rd=root_datum();

  RankFlags negatives,new_negatives;
  Weight mu;

  do
  {
    new_negatives.reset();
    mu=get_projection(negatives).projection*lambda;
    for (size_t i=0; i<rd.semisimple_rank(); ++i)
      if (not negatives[i] and mu.dot(rd.simpleCoroot(i))<0)
	new_negatives.set(i);
    negatives |= new_negatives;
  }
  while (new_negatives.any());

  level sp=mu.dot(rd.dual_twoRho());
  level d=2*get_projection(negatives).denom; // double to match |height|
  return (sp+d-1)/d; // round upwards, since height is always integer
} // |SRK_context::height_bound|


bool SRK_context::isStandard(const StandardRepK& sr) const
{
  const RootDatum& rd=root_datum();
  Weight lambda=lift(sr);
  const Fiber& f=fiber(sr);

  for (size_t i=0; i<f.imaginaryRank(); ++i)
    if (lambda.dot(rd.coroot(f.simpleImaginary(i)))<0)
    {
      offender=f.simpleImaginary(i); return false;
    }

  return true;
}

bool SRK_context::isNormal(Weight lambda, CartanNbr cn) const
{
  const auto& sc = info(cn).sum_coroots;
  for (unsigned int i=0; i<sc.size(); ++i)
    if (lambda.dot(sc[i])<0)
    {
      offender=i; return false; // |offender| indicates a complex simple root
    }

  return true;
}

bool SRK_context::isZero(const StandardRepK& sr) const
{
  const RootDatum& rd=root_datum();
  Weight lambda=lift(sr);
  const Fiber& f=fiber(sr);
  TitsElt a=titsElt(sr);

  for (size_t i=0; i<f.imaginaryRank(); ++i)
    if (not basedTitsGroup().grading(a,f.simpleImaginary(i)) // i.e., compact
	and lambda.dot(rd.coroot(f.simpleImaginary(i)))==0)
    {
      offender=f.simpleImaginary(i); return true;
    }

  return false;
}

bool SRK_context::isFinal(const StandardRepK& sr) const
{
  const RootDatum& rd=root_datum();
  Weight lambda=lift(sr);
  const Fiber& f=fiber(sr);

  // since coordinates are doubled, the scalar product below is always even
  for (size_t i=0; i<f.realRank(); ++i)
    if (lambda.dot(rd.coroot(f.simpleReal(i)))%4 == 0)
    {
      offender=f.simpleReal(i); return false;
    }

  return true;
}

void SRK_context::normalize(StandardRepK& sr) const
{
  const RootDatum& rd = root_datum();
  CartanNbr cn = sr.Cartan();
  const Cartan_info& ci = info(cn);
  Weight lambda = lift(sr);

  while (not isNormal(lambda,cn))
  { auto ind=witness(); // index into |bi_ortho| bitset
    auto i = ci.bi_ortho.n_th_bit(ind); // simple root for |ind| in |bi_ortho|
    lambda -= rd.simpleRoot(i)*lambda.dot(ci.coroot_sum(ind));
  }

  sr.d_lambda = project(cn,lambda);
} // |normalize|


q_Char SRK_context::q_normalize_eq (const StandardRepK& sr,size_t i) const
{
  return q_reflect_eq(sr,i,lift(sr),info(sr.Cartan()).coroot_sum(i));
} // |SRK_context::q_normalize_eq|


q_Char SRK_context::q_reflect_eq(const StandardRepK& sr,size_t i,
				 Weight lambda, // doubled
				 const Weight& cowt) const
{
  const RootDatum& rd = root_datum();
  const TorusPart& x = sr.d_fiberElt;
  CartanNbr cn = sr.Cartan();

  int n = -lambda.dot(cowt);

  assert (n>0);
  assert (n%2==0);  // because of doubled coordinates
  n/=2;

  Weight a2 = rd.simpleRoot(i)*2; // doubled coordinates

  lambda += a2*n; // $\lambda_0 = s_\alpha(\lambda)$

  q_Char result(StandardRepK(cn,x,project(cn,lambda)),
		q_CharCoeff(1,1)); // start with $q*srep(\lambda_0)$
  int j;
  for (j=1; 2*j<n; ++j)
  {
    q_CharCoeff coef(j+1,1); // $q^{j+1}$
    coef -= q_CharCoeff(j-1,1); // $q^{j-1}$
    lambda -= a2; // $\lambda_0 - j*alpha$
    result+=(q_Char(StandardRepK(cn,x,project(cn,lambda)),coef));
  }

  if (2*j==n) // final term for even length ladders
  {
    q_CharCoeff coef(j,1); // $q^{n/2}$
    coef -= q_CharCoeff(j-1,1); // $q^{n/2-1}$
    lambda -= a2; // $\lambda_0 - n/2*alpha$
    assert (lambda.dot(cowt)==0);
    result+=(q_Char(StandardRepK(cn,x,project(cn,lambda)),coef));
  }
  else
    assert (lambda.dot(cowt)==2); // dot product went from $2n$, steps of $-4$

  return result;
} // |SRK_context::q_reflect_eq|


/*
  The purpose of |theta_stable_parabolic| is to move to a situation in which
  |dom| is dominant, and in which the only positive roots sent to negative
  ones by the involution $\theta$ are the real ones. Then the real roots
  define the Levi factor, and the remaining positive roots the radical of a
  $\theta$-stable parabolic subalgebra. The involution $\theta$ is given by
  the twisted involution in |strong|; we think of keeping it fixed while
  gradually moving around the set of positive roots, each time replacing some
  complex root in the simple basis by its opposite. In reality the positive
  system is fixed, and the moves conjugate $\theta$ and reflect |dom|. In fact
  we shall need a fiber part as well as an involution once we have obtained
  our goal, so conjugation is in fact |basedTwistedConjugate| on |strong|.
 */
PSalgebra
SRK_context::theta_stable_parabolic
  (const StandardRepK& sr, WeylWord& conjugator) const
{
  const RootDatum& rd=root_datum();
  const TwistedWeylGroup& W=twistedWeylGroup();

  Weight dom=theta_lift(sr); // theta-stable weight, to be made dominant
  TitsElt strong=titsElt(sr);

  WeylWord ww; // conjugating element

  /* the following loop terminates because we either increase the number of
     positive coroots with strictly positive evaluation on |dom|, or we keep
     that number constant and decrease the number of positive complex roots
     that the involution sends to negative ones.
  */
  while (true) // loop will terminate if inner loop runs to completion
  {
    weyl::Generator s;
    for (s=0; s<rd.semisimple_rank(); ++s)
    {
      RootNbr alpha=rd.simpleRootNbr(s);
      LatticeCoeff v=dom.dot(rd.simpleCoroot(s));

      if (v<0) // first priority: |dom| should be made dominant
	break; // found value of |s| to use in conjugation/reflection
      else if (v>0) continue; // don't touch |alpha| in this case

      // now |dom| is on reflection hyperplane for |alpha|

      // second priority give |alpha| and its $\theta$ image the same sign
      RootNbr beta= // image of |alpha| by $\theta$
	rd.permuted_root(W.word(strong.w()),rd.simpleRootNbr(W.twisted(s)));
      if (rd.is_negroot(beta) and beta!=rd.rootMinus(alpha))
	break; // found |s| in this case as well

    } // |for(s)|

    if (s<rd.semisimple_rank()) // then we found a reflection |s| to apply
    {
      basedTitsGroup().basedTwistedConjugate(strong,s);
      rd.simple_reflect(s,dom);
      ww.push_back(s);
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

   Also |ww| is such that |W.twistedConjugate(strong.tw(),ww)|
   would make |strong.tw()| equal to its original value again.
*/

  // Build the parabolic subalgebra:

  { // first ensure |strong| is reduced
    const WeightInvolution theta = innerClass().matrix(strong.tw());
    strong.reduce(tits::fiber_denom(theta));
  }

  ww.swap(conjugator); // report the conjugating element we found
  return PSalgebra(strong,innerClass());

} // |theta_stable_parabolic|

KGBEltList SRK_context::sub_KGB(const PSalgebra& q) const
{
  BitMap flagged(kgb().size());

  KGBElt origin=UndefKGB;
  {
    KGBEltPair packet=kgb().tauPacket(q.involution());
    KGBElt x;
    for (x=packet.first; x<packet.second; ++x)
      if (kgb().titsElt(x)==q.strong_involution())
      {
	origin=x; break;
      }
    assert(x<packet.second); // search should succeed
  }

  flagged.insert(origin);
  containers::queue<KGBElt> queue { origin };
  do
  {
    KGBElt x=queue.front(); queue.pop();
    for (RankFlags::iterator it=q.Levi_gens().begin(); it(); ++it)
    {
      KGBElt y=kgb().cross(*it,x);
      if (not flagged.isMember(y))
      {
	flagged.insert(y); queue.push(y);
      }
      y=kgb().inverseCayley(*it,x).first; // second will follow if present
      if (y!=UndefKGB and not flagged.isMember(y))
      {
	flagged.insert(y); queue.push(y);
      }
    }
  }
  while (not queue.empty());

  return KGBEltList(flagged.begin(),flagged.end());
} // |sub_KGB|

RawChar SRK_context::KGB_sum(const PSalgebra& q, const Weight& lambda) const
{
  const RootDatum& rd=root_datum();
  KGBEltList sub=sub_KGB(q); std::reverse(sub.begin(),sub.end());
  // now |sub[0]| is the largest |x| value in |sub|

  std::vector<size_t> sub_inv(kgb().size(),~0ul);

  for (size_t i=0; i<sub.size(); ++i)
    sub_inv[sub[i]]=i; // partially fill array with inverse index

  std::vector<Weight> mu; // list of $-\rho$-centered weights,
  mu.reserve(sub.size()); // associated to the elements of |sub|

  mu.push_back(lambda); (mu[0]-=rd.twoRho())/=2; // |mu| is $-\rho$-centered

  // for remaining |sub| elements, deduce its |mu| from that of some ascent
  for (size_t i=1; i<sub.size(); ++i)
  {
    KGBElt x=sub[i];
    RankFlags::iterator it;
    for (it=q.Levi_gens().begin(); it(); ++it)
    {
      if (kgb().cross(*it,x)>x) // then we can use ascending cross action
      {
	size_t k=sub_inv[kgb().cross(*it,x)];
	assert(k!=~0ul); // we ought to land in the subset
	// apply implicitly $-\rho$-centered refection (sign flip done later)
	mu.push_back(rd.simple_reflection(*it,mu[k]));
	break;
      }
    }
    if (it()) continue; // if we could use a cross action, we're done for |i|

    // now similarly try Cayley transforms
    for (it=q.Levi_gens().begin(); it(); ++it)
    {
      if (kgb().cayley(*it,x)!=UndefKGB) // then we can use this Cayley
      {
	size_t k=sub_inv[kgb().cayley(*it,x)];
	assert(k!=~0ul); // we ought to land in the subset
	Weight nu=mu[k]; // $\rho-\lambda$ upstairs
	assert(nu.dot(rd.simpleCoroot(*it))%2 == 0); // finality (non-parity)
	Weight alpha=rd.simpleRoot(*it);
	nu -= (alpha *= nu.dot(rd.simpleCoroot(*it))/2); // project
	mu.push_back(nu); // use projected weight downstairs
	break;
      }
    } // |for(it)|
    assert(it()); // if no cross ascent worked, some Cayley ascent must have
  } // |for i|
  assert(mu.size()==sub.size());

  size_t max_l=kgb().length(sub[0]);

  RawChar result;
  for (size_t i=0; i<sub.size(); ++i)
  {
    KGBElt x=sub[i];
    RawRep r(mu[i],kgb().titsElt(x));
    result += RawChar(r, ((max_l-kgb().length(x))%2 == 0 ? 1 : -1));
  }

  return result;
} // |KGB_sum|

// Express irreducible K-module as a finite virtual sum of standard ones
CharForm
SRK_context::K_type_formula(const StandardRepK& sr, level bound)
{
  const RootDatum& rd=root_datum();

  // Get from |sr| to some $\theta$-stable parabolic subalgebra

  WeylWord conjugator;
  PSalgebra p = theta_stable_parabolic(sr,conjugator);

  Weight lambda= rd.image_by_inverse(lift(sr),conjugator); // towards |p|

  RawChar KGB_sum_p= KGB_sum(p,lambda);

  typedef free_abelian::Monoid_Ring<Weight> polynomial;
  const InvolutionTable& i_tab = innerClass().involution_table();

  // type of formal linear combination of weights, associated to Tits element

  Char result;
  for (const auto& term : KGB_sum_p)
  {
    Char::coef_t c=term.second; // coefficient from |KGB_sum_p|
    const Weight& mu=term.first.first; // weight from |KGB_sum_p|
    const TitsElt& strong=term.first.second; // Tits elt from |KGB_sum_p|
    const InvolutionNbr i_theta = i_tab.nr(strong.tw());
    const WeightInvolution& theta = i_tab.matrix(i_theta);

    if (height_bound(mu+theta*mu+i_tab.theta_plus_1_rho(i_theta)) > bound)
      continue;

    RootNbrSet A(rd.numRoots()); // oversized: negative roots unused
    // collect in |A| all positive roots that are nci or first of a C+ pair
    for (BitMap::iterator rt=p.radical().begin(); rt(); ++rt)
    {
      RootNbr alpha=*rt;
      assert(not i_tab.real_roots(i_theta).isMember(alpha));
      if (i_tab.imaginary_roots(i_theta).isMember(alpha))
	A.set_to(alpha,basedTitsGroup().grading(strong,alpha)); // add if nc
      else // complex root
      {
	RootNbr beta=i_tab.root_involution(i_theta,alpha);
	assert(rd.is_posroot(beta));
	A.set_to(alpha,beta>alpha); // add first of two paired complex posroots
      }
    }

//     std::cout << "Sum over subsets of " << A.size() << " roots, giving ";


    // compute in |pol| the product $X^\mu*\prod_{\alpha\in A}(1-X^\alpha)$
    polynomial pol(mu);
    for (RootNbrSet::iterator it=A.begin(); it!=A.end(); ++it)
    {
      polynomial copy=pol; // since |add_multiple| assumes no aliasing
      pol.add_multiple(copy,-1,rd.root(*it));

      // filter out terms that cannot affect anything below |bound|
      for (polynomial::iterator term=pol.begin(); term!=pol.end();)
      {
	const Weight& nu=term->first; // $-\rho$-centered
	if (height_bound(nu+theta*nu+i_tab.theta_plus_1_rho(i_theta)) > bound)
	  pol.erase(term++);
	else
	  term++;
      }
    }
//     std::cout << pol.size() << " terms." << std::endl;

    // iterate over terms in formal sum, taking coef *= |c|
    for (polynomial::const_iterator term=pol.begin(); term!=pol.end(); ++term)
    { // contribute |c*term| at |strong| involution to |result|
      Weight lambda_rho=term->first;
      polynomial::coef_t coef=term->second;
      result += Char(std_rep_rho_plus(lambda_rho,strong),c*coef);
    }
  } // |for (term : KGB_sum_p)
  return std::make_pair(sr, result);
} // |K_type_formula|

Raw_q_Char SRK_context::q_KGB_sum(const PSalgebra& p,
				  const Weight& lambda) const
{
  const RootDatum& rd=root_datum();
  KGBEltList sub=sub_KGB(p); std::reverse(sub.begin(),sub.end());

  std::vector<size_t> sub_inv(kgb().size(),~0ul);

  for (size_t i=0; i<sub.size(); ++i)
    sub_inv[sub[i]]=i; // partially fill array with inverse index

  std::vector<Weight> mu; // list of $\rho$-centered weights,
  mu.reserve(sub.size());               // associated to the elements of |sub|

  mu.push_back(lambda); (mu[0]-=rd.twoRho())/=2; // make $\rho$-centered

  for (size_t i=1; i<sub.size(); ++i)
  {
    KGBElt x=sub[i];
    RankFlags::iterator it;
    for (it=p.Levi_gens().begin(); it(); ++it)
    {
      if (kgb().cross(*it,x)>x) // then we can use ascending cross action
      {
	size_t k=sub_inv[kgb().cross(*it,x)];
	assert(k!=~0ul); // we ought to land in the subset
	mu.push_back(rd.simple_reflection(*it,mu[k])); // $\rho$-centered
	break;
      }
    }
    if (it()) continue; // if we could use a cross action, we're done for |i|

    // now similarly try Cayley transforms
    for (it=p.Levi_gens().begin(); it(); ++it)
    {
      if (kgb().cayley(*it,x)!=UndefKGB) // then we can use this Cayley
      {
	size_t k=sub_inv[kgb().cayley(*it,x)];
	assert(k!=~0ul); // we ought to land in the subset
	Weight nu=mu[k]; // $\rho-\lambda$ at split side
	assert(nu.dot(rd.simpleCoroot(*it))%2 == 0); // finality
	Weight alpha=rd.simpleRoot(*it);
	nu -= (alpha *= nu.dot(rd.simpleCoroot(*it))/2); // project
	mu.push_back(nu); // use projected weight at compact side of transform
	break;
      }
    }
    assert(it()); // if no cross action worked, some Cayley transform must have
  }
  assert(mu.size()==sub.size());

  size_t max_l=kgb().length(sub[0]);

  Raw_q_Char result;
  for (size_t i=0; i<sub.size(); ++i)
  {
    KGBElt x=sub[i];
    RawRep r(mu[i],kgb().titsElt(x));
    size_t codim=max_l-kgb().length(x);
    result += Raw_q_Char(r,q_CharCoeff(codim,codim%2==0 ? 1 : -1)); // $(-q)^c$
  }

  return result;
} // |q_KGB_sum|

// Express irreducible K-module as a finite $q$-virtual sum of "standard" ones
q_CharForm
SRK_context::q_K_type_formula(const StandardRepK& sr, level bound)
{
  const RootDatum& rd=root_datum();
  const InvolutionTable& i_tab = innerClass().involution_table();

  // Get theta stable parabolic subalgebra

  WeylWord conjugator;
  PSalgebra p = theta_stable_parabolic(sr,conjugator);

  Weight lambda= rd.image_by_inverse(lift(sr),conjugator); // towords |p|

  Raw_q_Char q_KGB_sum_p= q_KGB_sum(p,lambda);

  // type of formal linear combination of weights, associated to Tits element

  q_Char result;
  for (Raw_q_Char::const_iterator
	 it=q_KGB_sum_p.begin(); it!=q_KGB_sum_p.end(); ++it)
  {
    q_CharCoeff c=it->second; // coefficient from |q_KGB_sum|
    const Weight& mu=it->first.first; // weight from |q_KGB_sum|
    const TitsElt& strong=it->first.second; // Tits elt from |q_KGB_sum|
    InvolutionData id = innerClass().involution_data(strong.tw());

    RootNbrSet A(rd.numRoots());
    for (BitMap::iterator
	   rt=p.radical().begin(); rt!=p.radical().end(); ++rt)
    {
      RootNbr alpha=*rt;
      assert(not id.real_roots().isMember(alpha));
      if (id.imaginary_roots().isMember(alpha))
	A.set_to(alpha,basedTitsGroup().grading(strong,alpha)); // add if nc
      else // complex root
      {
	RootNbr beta=id.root_involution(alpha);
	assert(rd.is_posroot(beta));
	A.set_to(alpha,beta>alpha); // add first of two complex roots
      }
    }

//     std::cout << "Sum over subsets of " << A.size() << " roots, giving ";

    typedef free_abelian::Monoid_Ring<Weight,q_CharCoeff>
      polynomial; // with weight exponents and $q$-polynomials as coefficients
    const InvolutionNbr i_theta = i_tab.nr(strong.tw());
    const WeightInvolution& theta = i_tab.matrix(i_theta);

    // compute $X^\mu*\prod_{\alpha\in A}(1-X^\alpha)$ in |pol|
    polynomial pol(mu);
    for (RootNbrSet::iterator it=A.begin(); it!=A.end(); ++it)
    {
      polynomial copy=pol; // since |add_multiple| assumes no aliasing
      pol.add_multiple(copy,q_CharCoeff(1,-1),rd.root(*it)); // $*(1-qX^\alpha)$

      // filter out terms that cannot affect anything below |bound|
      for (polynomial::iterator term=pol.begin(); term!=pol.end();)
      {
	const Weight& mu=term->first;
	if (height_bound(mu+theta*mu+i_tab.theta_plus_1_rho(i_theta)) > bound)
	  pol.erase(term++);
	else
	  term++;
      }
    }
//     std::cout << pol.size() << " terms." << std::endl;

    // iterate over terms in formal sum, taking coef *= |c|
    for (polynomial::const_iterator term=pol.begin(); term!=pol.end(); ++term)
    {
      Weight lambda=term->first;
      polynomial::coef_t coef=term->second;
      result += q_Char(std_rep_rho_plus(lambda,strong),c*coef); // contribute
    }
  } // for sum over KGB for L
  return std::make_pair(sr, result);
} // |q_K_type_formula|


HechtSchmid
SRK_context::HS_id(const StandardRepK& sr, RootNbr alpha) const
{
  const RootDatum& rd=root_datum();
  TitsElt a=titsElt(sr);
  Weight lambda=lift(sr); // is non-dominant for |alpha|, so |sr| is not normal
  assert(rd.is_posroot(alpha)); // indeed |alpha| simply-imaginary for |a.tw()|

  HechtSchmid id(sr); // start with |sr| on the left hand side
  size_t i=0; // simple root index (value will be set in following loop)
  while (alpha!=rd.simpleRootNbr(i=rd.find_descent(alpha)))
  {
    // now $\alpha_i$ is (not imaginary, so) complex for the involution
    // |a.tw()|; reflect all data by $s_i$, which decreases level of $\alpha$
    rd.simple_reflect_root(i,alpha);
    rd.simple_reflect(i,lambda);
    basedTitsGroup().basedTwistedConjugate(a,i);
  }

  Weight mu=rd.simple_reflection(i,lambda);
  if (basedTitsGroup().simple_grading(a,i))
  { // now $\alpha_i$ is a non-compact imaginary simple root
    basedTitsGroup().basedTwistedConjugate(a,i); // adds $m_i$ to torus part
    StandardRepK sr0= std_rep(mu, a);
    assert(sr.d_cartan==sr0.d_cartan);
    id.add_lh(sr0); // add |sr0| to |sr| as second LHS term; type II: are equal

    /* Now are equivalent:
       = |sr0.d_fiberElt==sr.d_fiberElt|
       = $m_i$ absorbed into quotient ($\alpha_i^vee \in 2X_*+X_*^{-\theta}$)
       = $\alpha^\vee( (X^*)^\theta ) = 2\Z$
       = $\alpha\notin(1-\theta_\alpha)X^*$  ($\theta_\alpha=\theta.s_\alpha$)
       = Cayley transform is type II
    */

    basedTitsGroup().Cayley_transform(a,i);
    StandardRepK sr1= std_rep(lambda, a);
    id.add_rh(sr1);
    if (sr0.d_fiberElt==sr.d_fiberElt) // type II
    {
      lambda += root_datum().simpleRoot(i); // other possibility
      lambda += root_datum().simpleRoot(i); // add twice: doubled coordinates
      StandardRepK sr2= std_rep(lambda, a);
      assert(sr1.d_lambda != sr2.d_lambda);
      id.add_rh(sr2);
    }
  }
  else // $\alpha_i$ is a compact root; "easy" Hecht-Schmid identity
    // based twisted conjugation fixes the Tits element; just reflect weight
    id.add_lh(std_rep(mu,a)); // and no RHS

  return id;
} // |HS_id|

/*
  The method |HS_id| is only called when |lambda| is non-dominant for |alpha|
  and therefore gives a second LHS term with a strictly more dominant weight.

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
SRK_context::back_HS_id(const StandardRepK& sr, RootNbr alpha) const
{
  const RootDatum& rd=root_datum();
  assert(rd.is_posroot(alpha)); // in fact it must be simple-real for |a.tw()|

  TitsElt a=titsElt(sr);

  Weight lambda = lift(sr);
  SmallSubspace mod_space=
    info(sr.d_cartan).fiber_modulus; // make a copy to be modified
  RankFlags orth; // becomes system orthogonal to |tl| below
  { // first assure the theta-lift of sr is dominant
    Weight tl=theta_lift(sr);
    WeylWord ww= rd.to_dominant(tl);
    rd.act(ww,tl); // now it is dominant
    alpha = rd.permuted_root(ww,alpha);
    assert(tl.dot(rd.coroot(alpha))==0); // because $\alpha$ is real
    basedTitsGroup().basedTwistedConjugate(ww,a);
    rd.act(ww,lambda);
    for (size_t i=ww.size(); i-->0; )
      mod_space.apply(dual_reflection(ww[i]));
    for (weyl::Generator i=0; i<rd.semisimple_rank(); ++i)
      orth.set(i,tl.dot(rd.simpleCoroot(i))==0);
  }
  assert(rd.is_posroot(alpha)); // no real reflections; should still be positive
  assert(orth.any()); // since root $\alpha$ is in span of orth. simple roots

  // basis used is of $(1/2)X^*$, so scalar product with coroot always even
  assert(lambda.dot(rd.coroot(alpha))%4 == 0); // the non-final condition

  { // adapt lift $\lambda$ to be orthogonal to $\alpha$
    Weight mu=rd.root(alpha);
    mu *= lambda.dot(rd.coroot(alpha))/2; // an even multiple of $\alpha$
    lambda -= mu; // this makes |lambda| invariant under $s_\alpha$
    assert(lambda.dot(rd.coroot(alpha)) == 0); // check projection
    // due to basis used, |lambda| effectively modified by multiple of $\alpha$
    /* in type I, $\alpha$ is in $(1-\theta)X^*$ and correction is neutral
       in type II, correction need not be in $(1-\theta)X^*$, but adding
       $\alpha$ gives HC parameter designating the same representation
    */
  }

  // the following loop terminates because $\alpha$ is in span of |orth|
  weyl::Generator i; // becomes simple root index of $\alpha$
  while (alpha!=rd.simpleRootNbr(i=(rd.descent_set(alpha)&orth).firstBit()))
  {
    rd.simple_reflect_root(i,alpha);
    basedTitsGroup().basedTwistedConjugate(a,i);
    rd.simple_reflect(i,lambda);
    mod_space.apply(dual_reflection(i));
  }

  // one right term is obtained by undoing Cayley for |a|, with lifted |lambda|
  basedTitsGroup().inverse_Cayley_transform(a,i,mod_space);

  HechtSchmid id(sr); // left hand side is |sr|
  StandardRepK sr1 = std_rep(lambda,a);
  id.add_rh(sr1);

  // there will be another term in case of a type I HechtSchmid identity
  // it differs only in the fiber part, by $m_i$; if this vanishes into the
  // quotient, the HechtSchmid identity is type II and nothing is added
  basedTitsGroup().basedTwistedConjugate(a,i);
  StandardRepK sr2=std_rep(lambda,a);
  if (sr1.d_fiberElt!=sr2.d_fiberElt) // type I
    id.add_rh(sr2);

  return id;
} // |back_HS_id|

q_Char
SRK_context::q_HS_id_eq(const StandardRepK& sr, RootNbr alpha) const
{
  const RootDatum& rd=root_datum();
  TitsElt a=titsElt(sr);
  Weight lambda=lift(sr);
  assert(rd.is_posroot(alpha)); // indeed |alpha| simple-imaginary for |a.tw()|

  // the following test is easiest before we move to |alpha| simple situation
  bool type_II = info(sr.Cartan()).fiber_modulus.contains
    (SmallBitVector(rd.coroot(alpha))); // reduced modulo 2

  size_t i=0; // simple root index (value will be set in following loop)
  while (true) // we shall exit halfway when $\alpha=\alpha_i$
  {
    while (not rd.is_descent(i,alpha))
    {
      ++i;
      assert(i<rd.semisimple_rank());
    }
    // now $\<\alpha,\alpha_i^\vee> > 0$ where $\alpha$ is simple-imaginary
    // and \alpha_i$ is complex for the involution |a.tw()|

    if (alpha==rd.simpleRootNbr(i)) break; // found it

    // otherwise reflect all data by $s_i$, which decreases level of $\alpha$
    rd.simple_reflect_root(i,alpha);
    rd.simple_reflect(i,lambda);
    basedTitsGroup().basedTwistedConjugate(a,i);
    i=0; // and start over
  }

  Weight alpha_vee = rd.simpleCoroot(i);
  q_Char result;
  if (basedTitsGroup().simple_grading(a,i))
  { // $\alpha_i$ is a non-compact imaginary simple root
    if (not type_II) // no point doing this if difference is absorbed anyway
      basedTitsGroup().basedTwistedConjugate(a,i); // adds $m_i$ to torus part

    result -= q_reflect_eq(std_rep(lambda,a),i,lambda,alpha_vee);

    basedTitsGroup().Cayley_transform(a,i);
    result += q_Char(std_rep(lambda,a), // $*q^{floor(n/2)}$, with $n=-\< , >$:
		     q_CharCoeff(-lambda.dot(alpha_vee)/4,1));
    if (type_II)
    {
      lambda += root_datum().simpleRoot(i)*2; // (doubled coordinates)
      result += q_Char(std_rep(lambda,a), // $*q^{same exponent as above}$
		       q_CharCoeff(-lambda.dot(alpha_vee)/4,1));
    }
  }
  else // $\alpha_i$ is a compact root; "easy" Hecht-Schmid identity
    // based twisted conjugation fixes the Tits element; just reflect weight
    result -= q_Char(std_rep(rd.simple_reflection(i,lambda),a),
		     q_CharCoeff(0,1)); // $q^0$

  return result;
} // |q_HS_id_eq|

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
    long int ac=std::abs(it->second);
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
      long int ac=std::abs(it->second[0]);
      if (ac!=1)
	strm << ac << '*';
    }
    else
    {
      bool parens = it->second.multi_term();
      strm<<(parens ? " + (" : " + ") << it->second << (parens ? ")*":"*");
    }
    print(strm,it->first);
  }
  return strm;
}



/*****************************************************************************

        Chapter III -- The |KhatContext| and |qKhetContext| classes

******************************************************************************/


KhatContext::KhatContext
  (RealReductiveGroup &GR)
    : SRK_context(GR)
    , nonfinal_pool(),final_pool()
    , nonfinals(nonfinal_pool), finals(final_pool)
    , height_of()
    , height_graded(height_of)
    , expanded()
{}

qKhatContext::qKhatContext
  (RealReductiveGroup &GR)
    : SRK_context(GR)
    , nonfinal_pool(),final_pool()
    , nonfinals(nonfinal_pool), finals(final_pool)
    , height_of()
    , height_graded(height_of)
    , expanded()
{}

/******** manipulators and accessors ***************************************/

seq_no KhatContext::match_final(const StandardRepK& sr)
{ seq_no n=finals.find(sr);
  if (n!=Hash::empty)
    return n;
  assert (isStandard(sr) and not isZero(sr) and isFinal(sr));
  assert(height_of.size()==final_pool.size());
  height_of.push_back(height(sr)); // store height
  return finals.match(sr); // expand table
}

/* transform a character |chi| to a |combination| of Normal Standard terms
   this version ensures the basic |standardize| is recursively called first
*/
combination KhatContext::standardize(const Char& chi)
{
  combination result(height_graded);

  for (Char::const_iterator it=chi.begin(); it!=chi.end(); ++it)
  {
    StandardRepK sr=it->first; // must have non-|const| value for |normalize|
    normalize(sr);
    result.add_multiple(standardize(sr),it->second);
  }

  return result;
}

// the basic case, |sr| is assumed normalized here
combination KhatContext::standardize(const StandardRepK& sr)
{
  { // first check if we've already done |sr|
    seq_no n=nonfinals.find(sr);
    if (n!=Hash::empty)
      return expanded[n]; // in this case an equation is known
  }

  if (isStandard(sr))
  {
    { // now check if we already know |sr| to be Final
      seq_no n=finals.find(sr);
      if (n!=Hash::empty)
	return combination(n,height_graded); // single term known to be final
    }

    if (isZero(sr))
    {
      witness(); // clear and ignore indication of which root said "zero"
      combination zero(height_graded);
      return equate(nonfinals.match(sr),zero); // add at end of |expanded|
    }

    if (isFinal(sr))
    {
      assert(height_of.size()==final_pool.size());
      height_of.push_back(height(sr));
      return combination(finals.match(sr),height_graded); // single term
    }

    // now |sr| is known to be Standard, but neither Zero nor Final

    HechtSchmid equation= back_HS_id(sr,witness());
    assert(equation.n_lhs()==1); // |back_HS_id| gives 1-term left hand side
    assert(equation.n_rhs()!=0); // and never a null right hand side

    // now recursively standardize all terms, storing rules
    combination result= standardize(equation.rhs());
    return equate(nonfinals.match(sr),result); // and add rule for |sr|
  } // if (isStandard(sr,witness))

  // now |sr| is not Standard; apply a Hecht-Schmid identity
  HechtSchmid equation= HS_id(sr,witness());
  assert(equation.n_lhs()==2); // all cases of |HS_id| produce 2-term lhs

  // recursively standardize all terms of |equation.equivalent()|, storing rules
  combination result= standardize(equation.equivalent());
  return equate(nonfinals.match(sr),result); // and add rule for |sr|
} // |KhatContext::standardize|

/* map a character to one containing only Normal Standard terms
   this version ensures the basic |standardize| is recursively called first */
q_combin qKhatContext::standardize(const q_Char& chi)
{
  q_combin result(height_graded);

  for (q_Char::const_iterator it=chi.begin(); it!=chi.end(); ++it)
    result.add_multiple(standardize(it->first),it->second);

  return result;
}

// the basic case, includes $q$-normalisation as well
q_combin qKhatContext::standardize(const StandardRepK& sr)
{
  { // first check if we've already done |sr|
    seq_no n=nonfinals.find(sr);
    if (n!=Hash::empty) // then seen before
      return expanded[n];
  }

  if (not isNormal(sr))
  {
    q_Char rhs = q_normalize_eq(sr,witness()); // expression to replace |sr| by
    q_combin result= standardize(rhs);         // convert recursively
    return equate(nonfinals.match(sr),result); // and add rule for |sr|
  }

  if (not isStandard(sr))
  {
    q_Char equiv = q_HS_id_eq(sr,witness());

    // recursively standardize all terms of |equiv|, storing rules
    q_combin result= standardize(equiv);
    return equate(nonfinals.match(sr),result); // and add rule for |sr|
  } // if (not isStandard(sr,witness))

  { // now check if we already know |sr| to be Final
    seq_no n=finals.find(sr);
    if (n!=Hash::empty)
      return q_combin(n,height_graded); // single term known to be final
  }

  if (isZero(sr))
  {
    witness(); // to clear |offender|
    q_combin zero(height_graded);
    return equate(nonfinals.match(sr),zero);
  }

  if (isFinal(sr))
  {
    assert(height_of.size()==final_pool.size());
    height_of.push_back(height(sr));
    return q_combin(finals.match(sr),height_graded); // single term
  }

  // now |sr| is known to be Standard, but neither Zero nor Final

  HechtSchmid equation= back_HS_id(sr,witness());
  assert(equation.n_lhs()==1); // |back_HS_id| gives 1-term left hand side
  assert(equation.n_rhs()!=0); // and never a null right hand side

  // now recursively standardize all terms, storing rules
  q_combin result= standardize(to_q(equation.rhs()));
  return equate(nonfinals.match(sr),result); // and add rule for |sr|
} // |qKhatContext::standardize|



combination KhatContext::truncate(const combination& c, level bound) const
{
  combination result(height_graded);
  for (combination::const_iterator it=c.begin(); it!=c.end(); ++it)
    if (height(it->first)<=bound)
      result.insert(result.end(),*it);

  return result;
}


q_combin qKhatContext::truncate(const q_combin& c, level bound) const
{
  q_combin result(height_graded);
  for (q_combin::const_iterator it=c.begin(); it!=c.end(); ++it)
    if (height(it->first)<=bound)
      result.insert(result.end(),*it);

  return result;
}


// Apply |K_type_formula| for known Final representation, and |standardize|
equation KhatContext::mu_equation(seq_no n, level bound)
{
  CharForm kf= K_type_formula(rep_no(n),bound);

  equation result(n,combination(height_graded));
  combination& sum=result.second;

  for (Char::const_iterator it=kf.second.begin(); it!=kf.second.end(); ++it)
  {
    StandardRepK sr=it->first; // must have non-|const| value for |normalize|
    normalize(sr);
    sum.add_multiple(truncate(standardize(sr),bound),it->second);
  }

  return result;
}

// Apply |K_type_formula| for known Final representation, and |standardize|
q_equation qKhatContext::mu_equation(seq_no n, level bound)
{
  q_CharForm kf= q_K_type_formula(rep_no(n),bound);

  q_equation result(n,q_combin(height_graded));
  q_combin& sum=result.second;

  for (q_Char::const_iterator it=kf.second.begin(); it!=kf.second.end(); ++it)
    sum.add_multiple(truncate(standardize(it->first),bound),it->second);

  return result;
}


/*
  Convert a system of equations into a vector, adding equations for all terms
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
  BitMap lhs(nr_reps()); // left hand sides of all equations seen

  containers::queue<equation> queue;

  for (std::set<equation>::iterator
	 it=system.begin(); it!=system.end(); ++it)
    if (height(it->first)<=bound) // filter equations with low enough LHS
    {
      queue.push(*it);
      lhs.insert(it->first); // mark left hand sides from the original system
    }

  std::vector<equation> result;

  while (not queue.empty())
  {
    // ensure bitmap provides space for all current terms
    if (nr_reps()>lhs.capacity())
      lhs.set_capacity( (nr_reps()+constants::posBits)
			& ~constants::posBits); // round up to wordsize

    const equation& cf=queue.front(); // an unprocessed formula
    assert(height(cf.first) <= bound); // all RHS terms are (assumed) low enough

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
	queue.push(mu_equation(term->first,bound));
      }
    }

    queue.pop(); // we are done with this formula
  }

  return result;
} // |KhatContext::saturate|


// this is perfectly similar to the previous method, but for |q_equation|
std::vector<q_equation>
qKhatContext::saturate(const std::set<q_equation>& system, level bound)
{
  BitMap lhs(nr_reps()); // left hand sides of all equations seen

  containers::queue<q_equation> queue;

  for (std::set<q_equation>::iterator
	 it=system.begin(); it!=system.end(); ++it)
    if (height(it->first)<=bound) // filter equations with low enough LHS
    {
      queue.push(*it);
      lhs.insert(it->first); // mark left hand sides from the original system
    }

  std::vector<q_equation> result;

  while (not queue.empty())
  {
    // ensure bitmap provides space for all current terms
    if (nr_reps()>lhs.capacity())
      lhs.set_capacity( (nr_reps()+constants::posBits)
			& ~constants::posBits); // round up to wordsize

    const q_equation& cf=queue.front(); // an unprocessed formula
    assert(height(cf.first) <= bound); // all RHS terms are (assumed) low enough

    result.push_back(q_equation(cf.first,q_combin(height_graded)));
    q_combin& rhs=result.back().second;

    for (q_combin::const_iterator
	   term=cf.second.begin(); term!=cf.second.end(); ++term)
    {
      assert(height(term->first) <= bound); // guaranteed by |mu_equation|
      rhs.insert(*term);
      if (not lhs.isMember(term->first)) // no formula for this term seen yet
      {
	lhs.insert(term->first);
	queue.push(mu_equation(term->first,bound));
      }
    }

    queue.pop(); // we are done with this formula
  }

  return result;
} // |qKhatContext::saturate|


matrix::Matrix<CharCoeff> KhatContext::K_type_matrix
 (std::set<equation>& eq_set,
  level bound,
  std::vector<seq_no>& new_order,
  matrix::Matrix<CharCoeff>* direct_p
  )
{
  std::vector<equation> system=saturate(eq_set,bound);

  matrix::Matrix<CharCoeff> loc; // local matrix, maybe unused
  matrix::Matrix<CharCoeff>& m = direct_p==nullptr ? loc : *direct_p;

  m=triangularize(system,new_order);

  return matrix::inverse_triangular<false>(m);

} // |KhatContext::K_type_matrix|

matrix::Matrix<q_CharCoeff> qKhatContext::K_type_matrix
 (std::set<q_equation>& eq_set,
  level bound,
  std::vector<seq_no>& new_order,
  matrix::Matrix<q_CharCoeff>* direct_p
  )
{
  std::vector<q_equation> system=saturate(eq_set,bound);

  matrix::Matrix<q_CharCoeff> loc; // local matrix, maybe unused
  matrix::Matrix<q_CharCoeff>& m = direct_p==nullptr ? loc : *direct_p;

  m=triangularize(system,new_order);

  return matrix::inverse_triangular<false>(m);

} // |qKhatContext::K_type_matrix|


const combination& KhatContext::equate (seq_no n, const combination& rhs)
{
  assert(n==expanded.size()); // left hand side should be a new |StandardRepK|
  expanded.push_back(rhs);
  return rhs;
}

const q_combin& qKhatContext::equate (seq_no n, const q_combin& rhs)
{
  assert(n==expanded.size()); // left hand side should be a new |StandardRepK|
  expanded.push_back(rhs);
  return rhs;
}

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

    assert(eq.second.begin()->first==head->first); // so our loop will advance

    result += combination(eq.first,head->second,height_graded);
    remainder.add_multiple(eq.second,-head->second);
  }
  while (not remainder.empty());

  return result;
} // |KhatContext::branch|


q_combin qKhatContext::branch(seq_no s, level bound)
{
  q_combin result(height_graded); // a linear q_combin of $K$-types

  if (height(s)>bound)
    return result;

  // a linear q_combin of Final representations
  q_combin remainder(s,height_graded); // terms will have |height<=bound|
  do
  {
    q_combin::iterator head=remainder.begin(); // leading term

    q_equation eq=mu_equation(head->first,bound);

    result += q_combin(eq.first,head->second,height_graded);
    remainder.add_multiple(eq.second,-head->second);
  }
  while (not remainder.empty());

  return result;
} // |qKhatContext::branch|



std::ostream& KhatContext::print(std::ostream& strm,
				 const combination& ch,
				 bool brief) const
{
  if (ch.empty())
    return strm << '0';
  for (combination::const_iterator it=ch.begin(); it!=ch.end(); ++it)
  {
    strm << (it->second>0 ? " + " : " - ");
    long int ac = std::abs(it->second);
    if (ac!=1)
      strm << ac << '*';
    if (brief)
      strm << 'R' << it->first;
    else print(strm,rep_no(it->first));
  }
  return strm;
}

std::ostream& qKhatContext::print
  (std::ostream& strm, const q_combin& ch, bool brief)  const
{
  if (ch.empty())
    return strm << '0';
  for (q_combin::const_iterator it=ch.begin(); it!=ch.end(); ++it)
  {
    if (it->second.degree()==0)
    {
      strm << (it->second[0]>0 ? " + " : " - ");
      long int ac = std::abs(it->second[0]);
      if (ac!=1)
	strm << ac << '*';
    }
    else
    {
      bool parens = it->second.multi_term();
      strm<<(parens ? " + (" : " + ") << it->second << (parens ? ")*":"*");
    }
    if (brief)
      strm << 'R' << it->first;
    else print(strm,rep_no(it->first));
  }
  return strm;
}

void KhatContext::go(const StandardRepK& initial)
{
  StandardRepK sr = initial;
  normalize(sr);
  combination chi=standardize(sr);

#ifdef VERBOSE
  if (nonfinal_pool.size()>0)
  {
    const RootDatum& rd=root_datum();
    std::cout << "Intermediate representations:\n";
    for (size_t i=0; i<nonfinal_pool.size(); ++i)
    {
      const StandardRepK& sr=nonfinal_pool[i];
      print(std::cout << 'N' << i << ": ",sr) << " [" << height(sr) << ']';

      if (not isStandard(sr))
	std::cout << ", non Standard, witness "
		  << rd.coroot(witness());
      if (isZero(sr))
	std::cout << ", Zero, witness "
		  << rd.coroot(witness());
      if (not isFinal(sr))
	std::cout << ", non Final, witness "
		  << rd.coroot(witness());
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


} // |KhatContext::go|


/*****************************************************************************

        Chapter IV -- The PSalgebra class

******************************************************************************/

PSalgebra::PSalgebra(TitsElt base, const InnerClass& G)
    : strong_inv(base)
    , cn(G.class_number(base.tw()))
    , sub_diagram() // class |RankFlags| needs no dimensioning
    , nilpotents(G.root_datum().posroot_set())
{
  const RootDatum& rd=G.root_datum();
  InvolutionData id = G.involution_data(base.tw());

  // get |rd.simple_root_set()&id.real_roots()|, shifted to fit in |RankFlags|
  for (size_t i=0; i<rd.semisimple_rank(); ++i)
    if (id.real_roots().isMember(rd.simpleRootNbr(i)))
      sub_diagram.set(i);

  nilpotents.andnot(id.real_roots());
}



// ****************** Chapter V -- functions ************************


q_combin to_q(const combination& c)
{
  q_combin result(c.key_comp()); // use same comparison object
  for (combination::const_iterator it=c.begin(); it!=c.end(); ++it)
    result += q_combin(it->first,q_CharCoeff(it->second),c.key_comp());
  return result;
}


q_Char to_q(const Char& chi)
{
  q_Char result(chi.key_comp()); // use same comparison object
  for (Char::const_iterator it=chi.begin(); it!=chi.end(); ++it)
    result += q_Char(it->first,q_CharCoeff(it->second),chi.key_comp());
  return result;
}

template <typename C>
  matrix::Matrix_base<C> triangularize
   (const std::vector<
      std::pair<seq_no,
                Free_Abelian<seq_no,C,graded_compare>
               > >& eq,
	       std::vector<seq_no>& new_order)
{
  size_t n=eq.size();

  matrix::Matrix_base<C> M(n,n,C(0));
  graph::OrientedGraph incidence(n);

  for (size_t j=0; j<n; ++j) // loop over equations
  {
    size_t n_terms=0;
    for (size_t i=0; i<n; ++i) // loop over left hand sides
      if ((M(i,j)=eq[j].second[eq[i].first])!=C(0))
      { // |OrientedGraph::cells| puts sinks in front, so record edge $i\to j$.
	incidence.edgeList(i).push_back(j);
	++n_terms;
      }

    if (eq[j].second.size()!=n_terms) // then there are RHS terms outside the system
      throw std::runtime_error ("triangularize: system not saturated");
  }

  Partition order = incidence.cells(nullptr); // run Tarjan for equivalences only

  if (order.classCount()<n) // then there are non-singleton classes: error
  {
    auto sizes = order.class_sizes();
    for (size_t i=0; i<order.classCount(); ++i) // |i| indexes a class
      if (sizes[i]>1)
      {
	std::ostringstream os;bool first=true;
	for (size_t j=0; j<n; ++j)
	  if (order.class_of(j)==i)
	    os << (first?first=false,"triangularize: system has cycle: ":", ") << j;
	throw std::runtime_error (os.str());
      }
  }

  new_order.resize(n);
  for (size_t i=0; i<n; ++i) // |i| indexes class member
    new_order[order.class_of(i)]=eq[i].first; // put LHS equation |i| in its place

  matrix::Matrix_base<C> result(n,n,C(0));
  for (size_t i=0; i<n; ++i)
    for (graph::Vertex v : incidence.edgeList(i))
      result(order.class_of(i),order.class_of(v))=M(i,v);

  return result;
} // |triangularize|

template
  matrix::Matrix_base<Polynomial<int> > triangularize
    (const std::vector<q_equation>& eq, std::vector<seq_no>& new_order);

} // |namespace standardrepk|

// ****************** Chapter VI -- local functions ************************

namespace {

// orthogonal projection onto the intersection of kernels of coroots in |gens|
// the projection is parallel to the span of the roots in |gens|
int_Matrix
orth_projection(const RootDatum& rd, RankFlags gens,
		LatticeCoeff& d)
{
  size_t m=gens.count(), r=rd.rank();
  int_Matrix root_mat(r,m);
  int_Matrix sub_Cartan(m,m); // transposed, later inverted
  int_Matrix coroot_mat(m,r);
  for (RankFlags::iterator i=gens.begin(); i(); ++i)
  {
    size_t ii=gens.position(*i);
    for (RankFlags::iterator j=gens.begin(); j(); ++j)
      sub_Cartan(ii,gens.position(*j))=rd.cartan(*j,*i);
    for (size_t j=0; j<r; ++j)
    {
      root_mat(j,ii)=rd.simpleRoot(*i)[j];
      coroot_mat(ii,j)=rd.simpleCoroot(*i)[j];
    }
  }

  arithmetic::big_int denom;
  auto inv_sub_Cartan = inverse(sub_Cartan,denom);
  d = denom.convert<LatticeCoeff>();

  int_Matrix result(r,r,0); // set to identity scaled |denom|
  result += d;
  result -= root_mat * inv_sub_Cartan * coroot_mat;
  return result;
}

} // |namespace|


} // |namespace atlas|
