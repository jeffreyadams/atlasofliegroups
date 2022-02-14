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
  i_tab.lambda_unique(i_x,lambda_rho);
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
}

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
}

  } // |namespace repr|
} // |namespace atlas|
