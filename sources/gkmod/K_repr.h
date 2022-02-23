/*
  This is repr.h

  Copyright (C) 2022 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/
// Classes and utilities for manipulating $K$-representations

// This compilation unit is meant to ultimately replace the standardrepk unit

#ifndef K_REPR_H  /* guard against multiple inclusions */
#define K_REPR_H

#include "../Atlas.h"

#include "matrix.h"

namespace atlas {

namespace K_repr {

// A |K_type| represents either the $K$-restriction of a standard representation
// or an irreducible representation of $K$ (in other words a $K$-type).

// Either way, we use the representation of parameters, without the $\nu$ part

class K_type // compact representation of parameters at $\nu=0$
{
  friend class repr::Rep_context;
  KGBElt d_x;
  unsigned int hght;
  Weight lam_rho;

public:

  K_type(KGBElt x, const Weight& lam_rho,unsigned int h)
    : d_x(x), hght(h), lam_rho(lam_rho) {}

  K_type(KGBElt x, Weight&& lam_rho,unsigned int h)
    : d_x(x), hght(h), lam_rho(std::move(lam_rho)) {}

  K_type(K_type&&) = default;
  K_type& operator=(K_type&&) = default;

  KGBElt x () const { return d_x;  }
  const Weight& lambda_rho () const { return lam_rho; }
  unsigned int height () const { return hght; }

  K_type copy () const { return {d_x,lam_rho,hght}; }

  bool operator< (const K_type& another) const
  {
    if (hght!=another.hght)
      return hght<another.hght;
    if (d_x!=another.d_x)
      return d_x<another.d_x;
    assert(lam_rho.size()==another.lam_rho.size()); // this is always assumed
    // since |lambda_unique| is always called on |lam_rho|, we can compare them
    for (unsigned i=0; i<lam_rho.size(); ++i)
      if (lam_rho[i]!=another.lam_rho[i])
	return lam_rho[i]<another.lam_rho[i];
    return false; // we found equality
  }

  using Pooltype = std::vector<K_type>;
  bool operator!= (const K_type& another) const
  { return d_x!=another.d_x or lam_rho!=another.lam_rho; }
  size_t hashCode (size_t modulus) const
  {
    size_t h = 3*d_x;
    for (auto c : lam_rho)
      h = (17*h&(modulus-1)) + c;
    return h&(modulus-1);
  }
}; // |class K_type|

using K_type_pol = Free_Abelian_light<K_type,int>;

} // |namespace K_repr|

} // |namespace atlas|

#endif
