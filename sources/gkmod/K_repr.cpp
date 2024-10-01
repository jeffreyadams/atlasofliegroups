/*
  This is K_repr.cpp

  Copyright (C) 2022 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "K_repr.h"
#include "repr.h"

namespace atlas {
  namespace K_repr {

template<typename F>
const K_type_pol& K_type_to_pol_table::put (K_type t, F f)
{
  auto i = hash.find(t);
  if (i!=hash.empty)
    return poly[i];
  K_type_pol value = f(t.copy());
  return poly[hash.match(std::move(t))] = std::move(value);
}

const K_type_pol& K_type_to_pol_table::lookup (const K_type& t) const
{
  auto i = hash.find(t);
  if (i!=hash.empty)
    throw std::runtime_error("Looking up polynomial not stored in table");
  return poly[i];
}

  } // |namespace K_repr|

  namespace repr {

Weight // $(1+\theta_x)\lambda$
Rep_context::theta_plus_1_lambda (const K_repr::K_type& t) const
{
  const InvolutionTable& i_tab = involution_table();
  auto i_x = kgb().inv_nr(t.x());
  const auto& lr = t.lambda_rho();
  const auto& theta = i_tab.matrix(i_x);
  return theta*lr + lr + i_tab.theta_plus_1_rho(i_x);
}

K_repr::K_type Rep_context::sr_K(KGBElt x, Weight lambda_rho) const
{
  const InvolutionTable& i_tab = involution_table();
  auto i_x = kgb().inv_nr(x);
  i_tab.lambda_unique(i_x,lambda_rho); // ensure unique representative
  const auto& theta = i_tab.matrix(i_x);
  auto th1_lambda = lambda_rho+theta*lambda_rho+i_tab.theta_plus_1_rho(i_x);
  return {x,std::move(lambda_rho),height(th1_lambda)};
}

int // $\check\alpha\dot(1+\theta_x)\lambda$, often needed for its sign
Rep_context::theta_plus_1_eval (const K_repr::K_type& t, RootNbr alpha) const
{
  const auto& rd = root_datum();
  RootNbr beta = involution_table().root_involution(kgb().inv_nr(t.x()),alpha);
  const auto& lr = t.lambda_rho();
  return rd.coroot(alpha).dot(lr)+rd.colevel(alpha)
	+rd.coroot(beta).dot(lr)+rd.colevel(beta);
}

// |z| standard means $\lambda (weakly) dominant on the (simply-)imaginary roots
bool Rep_context::is_standard(const K_repr::K_type& z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = involution_table();

  for (RootNbr alpha : i_tab.imaginary_basis(i_x))
    if (rd.coroot(alpha).dot(z.lambda_rho())+rd.colevel(alpha)<0)
      return false;
  return true;
}

// |z| dominant means that $(1+theta_x)\lambda$ fails to be (weakly) dominant
bool Rep_context::is_dominant(const K_repr::K_type& z) const
{
  const RootDatum& rd = root_datum();

  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
    if ( theta_plus_1_eval(z,rd.simpleRootNbr(s)) < 0 )
      return false;

  return true;
}

// |z| zero means that no singular simply-imaginary roots are compact
// we assume |is_standard(z)|, so singular imaginary system is simply generated
bool Rep_context::is_nonzero(const K_repr::K_type& z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = involution_table();

  for (RootNbr alpha : i_tab.imaginary_basis(i_x)) // simple-imaginary
    if (rd.coroot(alpha).dot(z.lambda_rho())+rd.colevel(alpha)==0 // singular
	and not kgb().simple_imaginary_grading(z.x(),alpha)) // and compact
      return false;
  return true;
}

// |z| semifinal means the absence of real parity roots. We do not assume
// |is_theta_stable|, so all really-simple roots must be tested
bool Rep_context::is_semifinal(const K_repr::K_type& z) const
{
  const RootDatum& rd = root_datum();
  const InvolutionNbr i_x = kgb().inv_nr(z.x());
  const InvolutionTable& i_tab = involution_table();
  const RootNbrSet pos_real = i_tab.real_roots(i_x) & rd.posroot_set();
  const Weight test_wt = z.lambda_rho()*2 // $2(\lambda-\rho)$
	   + rd.twoRho()-rd.twoRho(pos_real); // replace $\rho$ by $\rho_R$

  for (RootNbr alpha : i_tab.real_basis(i_x))
    if (rd.coroot(alpha).dot(test_wt)%4 !=0) // doubled odd: parity real root
      return false; // which invalidates the semi-final condition
  return true;
}

// absence of complex singular descents; assumes properties asserted below
bool Rep_context::is_normal(const K_repr::K_type& z) const
{
  assert (is_standard(z)); // otherwise the notion is not defined
  assert (is_dominant(z)); // this is necessary fot the normal form

  // although we could define the predicate without the next two assumptions,
  // this would be harder to implement, while the predicate is less interesting
  assert (is_nonzero(z) and is_semifinal(z));

  const auto& rd = root_datum();

  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
    if (kgb().isComplexDescent(s,z.x()) and
	theta_plus_1_eval(z,rd.simpleRootNbr(s))==0)
      return false;

  return true;
} // |is_normal|

// dominance of $(1+theta_x)\lambda$ and absence of any singular descents
bool Rep_context::is_final(const K_repr::K_type& z) const
{
  const auto& rd = root_datum();
  KGBElt x = z.x();
  const auto& lr = z.lambda_rho();

  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
  {
    auto eval = theta_plus_1_eval(z,rd.simpleRootNbr(s));
    if (eval>0)
      continue;
    else if (eval<0)
      return false;

    // now simple coroot |s| is singular, check for descents
    switch (kgb().status(s,x))
    {
    case gradings::Status::ImaginaryCompact:
      return false; // fails |is_nonzero|
    case gradings::Status::Real:
      if (rd.simpleCoroot(s).dot(lr)%2 != 0) // then $\alpha_s$ is parity
	return false;
      break;
    case gradings::Status::Complex:
      if (kgb().isDescent(s,x))
	return false; // failed |is_normal|
      break;
    default:
      break; // ImaginaryNoncompact is fine
    } // |switch|
  } // |for (s)|

  return true;
}

bool
Rep_context::equivalent(K_repr::K_type z0, K_repr::K_type z1) const // by value
{
  if (kgb().Cartan_class(z0.x())!=kgb().Cartan_class(z1.x()))
    return false; // this non-equivalence requires little work to see
  to_canonical_involution(z0);
  to_canonical_involution(z1);
  // at the canonical involution equivalence means equal |x| and |lambda_rho|
  return z0==z1;
}

void Rep_context::make_dominant (K_repr::K_type& t) const
{
  KGBElt& x=t.d_x; Weight& lr = t.lam_rho; // operate directly on components

  const auto& rd = root_datum();
  const InvolutionTable& i_tab = involution_table();
  Weight im_wt; // imaginary weight to be made dominant for complex coroots

  { // set |im_wt|, and check dominance for imaginary subsystem
    const auto i_x = kgb().inv_nr(x);
    const auto& theta = i_tab.matrix(i_x);
    im_wt = lr + theta*lr + i_tab.theta_plus_1_rho(i_x);
    for (RootNbr alpha : i_tab.imaginary_basis(i_x))
      if (rd.coroot(alpha).dot(lr)+rd.colevel(alpha)<0)
	throw std::runtime_error("Non standard K-type in make_dominant");
  }

  { weyl::Generator s;
    do
      for (s=0; s<rd.semisimple_rank(); ++s)
	if (theta_plus_1_eval(t,rd.simpleRootNbr(s))<0)
	{
	  assert(i_tab.is_complex_simple(kgb().inv_nr(x),s));
	  x = kgb().cross(s,x);
	  rd.simple_reflect(s,im_wt);
	  rd.simple_reflect(s,lr,1); // $-\rho$-based reflection
	  break; // out of the loop |for(s)|
	} // |if(v<0)| and |for(s)|
    while (s<rd.semisimple_rank()); // wait until inner loop runs to completion
  }
  i_tab.lambda_unique(kgb().inv_nr(x),lr); // to preferred coset representative
} // |make_dominant|

// exhaust simple complex descents for the involution
void Rep_context::make_theta_stable (K_repr::K_type& t) const
{
  KGBElt& x=t.d_x; Weight& lr = t.lam_rho; // operate directly on components

  const auto& rd = root_datum();
  const InvolutionTable& i_tab = involution_table();

  { weyl::Generator s;
    do
      for (s=0; s<rd.semisimple_rank(); ++s)
	if (kgb().isComplexDescent(s,x))
	{
	  x = kgb().cross(s,x);
	  rd.simple_reflect(s,lr,1); // $-\rho$-based reflection
	  break; // out of the loop |for(s)|
	} // |if(v<0)| and |for(s)|
    while (s<rd.semisimple_rank()); // wait until inner loop runs to completion
  }
  i_tab.lambda_unique(kgb().inv_nr(x),lr); // to preferred coset representative
} // |make_theta_stable|

void Rep_context::to_canonical_involution
  (K_repr::K_type& t, RankFlags gens) const
{
  KGBElt& x=t.d_x; Weight& lr = t.lam_rho; // operate directly on components

  const auto& rd = root_datum();
  const InvolutionTable& i_tab = involution_table();

  TwistedInvolution sigma = kgb().involution(x);
  WeylWord ww = inner_class().canonicalize(sigma,gens);
  for (auto s : ww)
  {
    assert(i_tab.is_complex_simple(kgb().inv_nr(x),s));
    x = kgb().cross(s,x);
    rd.simple_reflect(s,lr,1); // $-\rho$-based reflection
  }
  i_tab.lambda_unique(kgb().inv_nr(x),lr); // to preferred coset representative
} // |to_canonical_involution|

void Rep_context::normalise(K_repr::K_type& z) const
{
  to_canonical_involution(z);
  // at the canonical involution equivalence means equal |x| and |lambda_rho|
  // so it only remains to ensure moving to a final class member if one exists

  const auto& rd = root_datum();
  KGBElt& x=z.d_x; Weight& lr = z.lam_rho; // operate directly on components

  { weyl::Generator s;
    do // a variation of |make_theta_stable|: exhaust singular complex descents
      for (s=0; s<rd.semisimple_rank(); ++s)
      { auto eval = theta_plus_1_eval(z,rd.simpleRootNbr(s));
	if (kgb().status(x).isComplex(s) and
	    (eval<0 or (eval==0 and kgb().isDescent(s,x))))
	{
	  x = kgb().cross(s,x);
	  rd.simple_reflect(s,lr,1); // $-\rho$-based reflection
	  break; // out of the loop |for(s)|
	} // |if(v<0)| and |for(s)|
      }
    while (s<rd.semisimple_rank()); // wait until inner loop runs to completion
  }
  involution_table().lambda_unique(kgb().inv_nr(x),lr);
}

using term = std::pair<K_repr::K_type, int>;
using term_list = simple_list<term>;

void insert(term&& trm, term_list& L) // assumes |L| is sorted decreasingly
{ auto it = L.begin();
  const auto& t = trm.first;
  while (not L.at_end(it))
    if (t < it->first)
      ++it; // skip higher terms
    else if (it->first < t) // then we are looking at lower terms
      break; // so break loop and insert before those lower terms
    else  // matching term; operate on coefficient
    { if ((it->second += trm.second) == 0)
        L.erase(it);
      return; // whether by coefficient update or erasure, we are done
    }
  L.insert(it,std::move(trm));
}

term_list Rep_context::finals_for(K_repr::K_type t) const
{
  const auto& rd = root_datum();
  const InvolutionTable& i_tab = involution_table();

  term_list result, todo;
  // there is no way to initialise |simple_list| from move-only type values
  todo.push_front(term(std::move(t),1)); // so move |t| to front term manually

  // the following loop is by decreasing height; new terms of lesser height
  // may be inserted into |result| en route, and will be encountered later
  while (not todo.empty())
  {
    auto  current = std::move(todo.front());
    auto& coef = current.second;
    KGBElt& x = current.first.d_x; // access current |K_type|, and
    Weight& lr = current.first.lam_rho; // operate directly on its components
    todo.pop_front();

    Weight im_wt; // imaginary weight to be made dominant
    { // set |im_wt|
      const auto i_x = kgb().inv_nr(x);
      const auto& theta = i_tab.matrix(i_x);
      im_wt = lr + theta*lr + i_tab.theta_plus_1_rho(i_x);
    }

  restart: // go here when |current| and |im_wt| have been modified
    for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
    { const auto eval = rd.simpleCoroot(s).dot(im_wt);
      if (eval<=0)
      { switch(kgb().status(s,x))
	{
	case gradings::Status::ImaginaryCompact:
	  if (eval<0) // then reflect |lambda| and negate sign
	  {
	    rd.simple_reflect(s,im_wt);
	    rd.simple_reflect(s,lr,1); // $-\rho$-based reflection
	    coef = -coef;
	    goto restart;
	  }
	  else goto drop;
	case gradings::Status::ImaginaryNoncompact:
	  if (eval<0) // then also add Cayley transform terms
	  {
	    KGBElt sx = kgb().cross(s,x);
	    KGBElt Cx = kgb().cayley(s,x);
	    K_repr::K_type t1 = sr_K(Cx,lr);
	    assert( t1.height() < current.first.height() );
	    insert(std::make_pair(std::move(t1),coef),todo);
	    if (sx==x) // then type 2 Cayley
	    {
	      K_repr::K_type t2 = sr_K(Cx,lr+rd.simpleRoot(s));
	      assert( t2.height() < current.first.height() );
	      insert(std::make_pair(std::move(t2),coef),todo);
	    }
	    x = sx; // after testing we can update |x| for nci cross action
	    rd.simple_reflect(s,im_wt);
	    rd.simple_reflect(s,lr,1); // $-\rho$-based reflection
	    coef = -coef; // reflect, negate, and continue with |current|
	    goto restart;
	  }
	  else // nothing to do for singular nci generator
	    continue; // continue loop on |s|
	case gradings::Status::Complex:
	  if (eval<0 or kgb().isDescent(s,x))
	  { // do cross action to make |im_wt| more dominant or |x| shorter
	    x = kgb().cross(s,x);
	    rd.simple_reflect(s,im_wt);
	    rd.simple_reflect(s,lr,1); // $-\rho$-based reflection
	    // keep |coef| unchanged
	    goto restart;
	  }
	  else continue; // nothing to do for a singular complex ascent
	case gradings::Status::Real:
	  assert(eval==0);
	  // now evaluation on |lr| counts
	  const auto eval_lr = rd.simpleCoroot(s).dot(lr);
	  if (eval_lr%2 != 0) // then $\alpha_s$ is a parity real root
	  {
	    lr -= rd.simpleRoot(s)*((eval_lr+1)/2);
	    assert( rd.simpleCoroot(s).dot(lr) == -1 );
	    KGBEltPair Cxs = kgb().inverseCayley(s,x);
	    if (Cxs.second!=UndefKGB)
	    {
	      K_repr::K_type t2 = sr_K(Cxs.second,lr);
	      assert( t2.height() == current.first.height() );
	      insert(std::make_pair(std::move(t2),coef),todo);
	    }
	    K_repr::K_type t1 = sr_K(Cxs.first,std::move(lr));
	    assert( t1.height() == current.first.height() );
	    insert(std::make_pair(std::move(t1),coef),todo);
	    goto drop; // we've destroyed |current|, don't contribute it
	  }
	  else continue; // nothing to do for a (singular) real nonparity root
	} // |switch|
      } // |if(eval<=0)|
    } // |for(s)|

    // now we have succeeded in making |current| standard; contriute it
    i_tab.lambda_unique(kgb().inv_nr(x),lr); // preferred coset representative
    result.push_front(std::move(current));

  drop: {} // when jumping here, proceed without contributing |current|
  } // |while(not todo.empty())|

  return result;
} // |Rep_context::finals_for|

sl_list<K_repr::K_type> Rep_context::KGP_set (K_repr::K_type& t) const
{
  const auto& rd = root_datum();
  const InvolutionTable& i_tab = involution_table();
  const auto& kgb = this->kgb();

  make_theta_stable(t);
  const InvolutionNbr i_theta = kgb.inv_nr(t.x());
  std::vector<weyl::Generator> Levi_gens;
  Levi_gens.reserve(rd.semisimple_rank());
  for (weyl::Generator s=0; s<rd.semisimple_rank(); ++s)
    if (i_tab.is_real_simple(i_theta,s))
    {
      assert(rd.simpleCoroot(s).dot(t.lambda_rho())%2==0); // |t| must be final
      Levi_gens.push_back(s);
    }

  BitMap present (kgb.size());
  present.insert(t.x());
  using pre_K_type = std::pair<KGBElt,Weight>;
  containers::queue<pre_K_type> Q { std::make_pair(t.x(), t.lambda_rho()) };
  sl_list<K_repr::K_type> result;
  result.push_back(t.copy()); // keep |t| itself for caller

  do
  {
    pre_K_type current = std::move(Q.front()); Q.pop();
    KGBElt x=current.first;
    const Weight& lam_rho = current.second;
    for (auto s : Levi_gens)
      switch (kgb.status(s,x))
      {
      case gradings::Status::Real:
	{
	  auto pair = kgb.inverseCayley(s,x);
	  KGBElt Csx; auto it = result.end();
	  auto eval = rd.simpleCoroot(s).dot(lam_rho);
	  assert(eval%2==0); // "non-parity"; from final condition
	  Weight new_lr = lam_rho - rd.simpleRoot(s)*(eval/2);
	  // with the first of |pair| more likely to be inserted, try it last
	  if ((Csx=pair.second)!=UndefKGB and not present.isMember(Csx))
	  {
	    present.insert(Csx);
	    result.insert(it,sr_K(Csx,new_lr));
	    Q.push(std::make_pair(Csx,new_lr));
	  }
	  if (not present.isMember(Csx=pair.first))
	  {
	    present.insert(Csx);
	    result.push_back(sr_K(Csx,new_lr));
	    Q.push(std::make_pair(Csx,std::move(new_lr)));
	  }
	}
	break;
      case gradings::Status::Complex:
	{
	  KGBElt sx=kgb.cross(s,x);
	  if (not present.isMember(sx)) // this also excludes complex ascents
	  {
	    present.insert(sx);
	    Weight new_lr = rd.simple_reflection(s,lam_rho);
	    result.push_back(sr_K(sx,new_lr));
	    Q.push(std::make_pair(sx,std::move(new_lr)));
	  }
	}
      default: break;
      } // |switch|
  }
  while (not Q.empty());

  return result;

} // |Rep_context::KGP_sum|

K_repr::K_type_pol Rep_context::monomial_product
  (const K_repr::K_type_pol& P, const Weight& e) const
{
  const InvolutionTable& i_tab = involution_table();
  K_repr::K_type_pol::poly result;
  result.reserve(P.size());
  for (const auto& term : P)
  {
    KGBElt x = term.first.x();
    auto new_exp = term.first.lambda_rho()+e;
    auto i_x = kgb().inv_nr(x);
    i_tab.lambda_unique(i_x,new_exp); // ensure unique representative
    const auto& theta = i_tab.matrix(i_x);
    auto ht = height(new_exp+theta*new_exp+i_tab.theta_plus_1_rho(i_x));
    result.emplace_back(K_repr::K_type{x,std::move(new_exp),ht},term.second);
  } // |for(term)|
  return // convert to |K_repr::K_type_pol|, sorting the shifted terms again
    { std::move(result), true, P.cmp() };
} // |Rep_context::monomial_product|

// compute height of "orthogonal projection to dominant cone" (closest point)
level Rep_context::height_bound (RatWeight lambda) const
/* this projection is dominant, and obtained by orthogonal projection onto the
   intersection of kernels of some set of simple coroots, say indexed by $S$,
   which is moreover such that the projection equals |lambda| plus a positive
   linear combination of the simple roots for $S$
*/
{
  const RootDatum& rd=root_datum();
  assert(lambda.size()==rd.rank());
  struct proj { weyl::Generator s; RatWeight v; };
  sl_list<proj> projectors;

  RankFlags S; weyl::Generator s;
  do
    for (s=0; s<rd.semisimple_rank(); ++s)
      if (not S.test(s) and lambda.dot_Q(rd.simpleCoroot(s)).is_negative())
      {
	RatWeight alpha (rd.simpleRoot(s),1); // to be projected to orthogonal
	for (const auto& p : projectors)
	  alpha -= p.v*alpha.dot_Q(rd.simpleCoroot(p.s));
	// finally ensure that |alpha.dot(rd.simpleCoroot(s)==1|
	(alpha /= alpha.dot_Q(rd.simpleCoroot(s))).normalize();
	// now project |lambda|, already orth to previous, orthogonal to |alpha|
	lambda -= alpha*lambda.dot_Q(rd.simpleCoroot(s));

	// finally extend |projectors| for future adjustments
	projectors.push_back(proj{s,std::move(alpha)});
	S.set(s);
	break; // and repeat outer loop from the beginning
      }
  while (s<rd.semisimple_rank());
  assert(is_dominant_ratweight(rd,lambda));
  auto d = lambda.denominator();
  return // |ceil(height(lambda))|
    (lambda.numerator().dot(rd.dual_twoRho())+d-1)/d;
}

K_repr::K_type_pol Rep_context::K_type_formula
  (K_repr::K_type& t,level max_level) const
{
  const auto& rd = root_datum();
  const InvolutionTable& i_tab = involution_table();
  auto terms = KGP_set(t); // this also moves |t| to a theta-stable parabolic
  // nilpotents of the parabolic subalgebra at |t|:
  auto max_l = kgb().length(t.x());
  RootNbrSet radical_posroots = rd.posroot_set();
  radical_posroots.andnot(i_tab.real_roots(kgb().inv_nr(t.x())));
  K_repr::K_type_pol result;

  for (auto&& term : terms)
  {
    KGBElt x = term.x(); const auto& lr = term.lambda_rho();
    const InvolutionNbr i_x = kgb().inv_nr(x);
    RatWeight lambda_0 // $(1+\theta)/2 * \lambda$
      ( lr + i_tab.matrix(i_x)*lr + i_tab.theta_plus_1_rho(i_x), 2 );
    if (height_bound(lambda_0)>max_level)
      continue;

    RootNbrSet sum_set{rd.numRoots()};
    for (const auto i : radical_posroots)
    { assert(not i_tab.real_roots(i_x).isMember(i)); // no new real roots here
      // though some of the non-radical posroots will no longer be real
       if (i_tab.complex_roots(i_x).isMember(i))
	sum_set.set_to(i,
	  i_tab.root_involution(i_x,i)>i); // first complex of a swapped pair
      else
	sum_set.set_to(i,
	  status(kgb(),x,i)==gradings::Status::ImaginaryNoncompact);
    }
    // |for(i : radical_posroots)|

    int sign = (max_l-kgb().length(x))%2==0 ? 1 : -1;
    K_repr::K_type_pol product (std::move(term),Split_integer(sign));
    for (const auto i : sum_set)
    {
      auto mp = monomial_product(product,rd.root(i));
      for (auto&& t : mp)
      {
	const auto& lr = t.first.lambda_rho();
	RatWeight lambda_0
	  ( lr + i_tab.matrix(i_x)*lr + i_tab.theta_plus_1_rho(i_x), 2);
	if (height_bound(lambda_0)<=max_level)
	  product.add_term(std::move(t.first),-t.second);
      }
    }
    for (auto&& t : product)
    { auto finals = finals_for(std::move(t.first));
      for (auto it=finals.begin(); not finals.at_end(it); ++it)
	if (it->first.height()<=max_level)
	  result.add_term (std::move(it->first),t.second*it->second);
    }
  }
  return result;
} // |Rep_context::K_type_formula|

K_repr::K_type_pol
  Rep_context::branch(K_repr::K_type_pol remainder, repr::level cutoff) const
{
  K_repr::K_type_pol result;
  auto it = remainder.begin();
  if (it==remainder.end())
    return result;
  size_t count=0;
  do
  {
    ++count;
    // periodically flatten |remainder| to remove leading zero terms
    if (count*count>2*remainder.size())
    {
      remainder=std::move(remainder).flatten();
      count=0;
    }

    auto it = remainder.begin();
    if (it->first.height()>cutoff)
    {
      remainder.erase(it); // drop any input terms that are already too high
      continue;
    }
    auto lead = it->first.copy(); // need a modifiable lvalue of K-type
    auto coef = it->second;
    result.add_term(lead.copy(),coef);
    remainder.add_multiple(K_type_formula(lead,cutoff),-coef);
  }
  while (not remainder.is_zero());
  return result;
} // |Rep_context::branch|, basic version


  } // |namespace repr|
} // |namespace atlas|
