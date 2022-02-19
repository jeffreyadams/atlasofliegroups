/*
  This is K_repr.cpp

  Copyright (C) 2022 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "K_repr.h"
#include "repr.h"

namespace atlas {
  namespace repr {

K_repr::K_type Rep_context::sr_K(KGBElt x, Weight lambda_rho) const
{
  const InvolutionTable& i_tab = involution_table();
  auto i_x = kgb().inv_nr(x);
  i_tab.lambda_unique(i_x,lambda_rho); // ensure unique representative
  const auto& theta = i_tab.matrix(i_x);
  auto th1_lambda = lambda_rho+theta*lambda_rho+i_tab.theta_plus_1_rho(i_x);
  return {x,lambda_rho,height(th1_lambda)};
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
    for (RootNbr s : i_tab.imaginary_basis(i_x))
      if (rd.simpleCoroot(s).dot(im_wt)<0)
	throw std::runtime_error("Non standard K-type in make_dominant");
  }

  { weyl::Generator s;
    do
      for (s=0; s<rd.semisimple_rank(); ++s)
	if (rd.simpleCoroot(s).dot(im_wt)<0)
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

using term = std::pair<K_repr::K_type,unsigned int>;
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
  todo.push_front(term(std::move(t),1u)); // so move |t| to front term manually

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
	    insert(std::make_pair(std::move(t1),current.second),todo);
	    if (sx==x) // then type 2 Cayley
	    {
	      K_repr::K_type t2 = sr_K(Cx,lr+rd.simpleRoot(s));
	      assert( t2.height() < current.first.height() );
	      insert(std::make_pair(std::move(t2),current.second),todo);
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
	  auto eval = rd.simpleCoroot(s).dot(lr);
	  if (eval%2 != 0) // then $\alpha_s$ is a parity real root
	  {
	    lr -= rd.simpleRoot(s)*((eval+1)/2);
	    assert( rd.simpleCoroot(s).dot(lr) == -1 );
	    KGBEltPair Cxs = kgb().inverseCayley(s,x);
	    K_repr::K_type t1 = sr_K(Cxs.first,lr);
	    assert( t1.height() == current.first.height() );
	    insert(std::make_pair(std::move(t1),current.second),todo);
	    if (Cxs.second!=UndefKGB)
	    {
	      K_repr::K_type t2 = sr_K(Cxs.second,lr);
	      assert( t2.height() == current.first.height() );
	      insert(std::make_pair(std::move(t2),current.second),todo);
	    }
	    goto drop; // we've destroyed |current|, don't contribute it
	  }
	  else continue; // nothing to do for a (singular) real nonparity root
	} // |switch|
      } // |if(eval<=0)|
    } // |for(s)|
    result.push_front(std::move(current));
  drop: {} // when jumping here, proceed without contributing |current|
  } // |while(not todo.empty())|
  return result;
}
  } // |namespace repr|
} // |namespace atlas|
