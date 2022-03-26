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
#include "hashtable.h"

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

  bool operator== (const K_type& other) const
  {
    return d_x==other.d_x and hght==other.hght and lam_rho==other.lam_rho;
  }
  bool operator != (const K_type& other) const { return not operator==(other); }
  bool operator< (const K_type& other) const
  {
    if (hght!=other.hght)
      return hght<other.hght;
    if (d_x!=other.d_x)
      return d_x<other.d_x;
    assert(lam_rho.size()==other.lam_rho.size()); // this is always assumed
    // since |lambda_unique| is always called on |lam_rho|, we can compare them
    for (unsigned i=0; i<lam_rho.size(); ++i)
      if (lam_rho[i]!=other.lam_rho[i])
	return lam_rho[i]<other.lam_rho[i];
    return false; // we found equality
  }

  using Pooltype = std::vector<K_type>;
  size_t hashCode (size_t modulus) const
  {
    size_t h = 21*d_x;
    for (auto c : lam_rho)
      h = (81*h&(modulus-1)) + c;
    return h&(modulus-1);
  }
}; // |class K_type|

using K_type_pol = Free_Abelian_light<K_type,Split_integer>;

class K_type_to_pol_table
{
  K_type::Pooltype pool;
  HashTable<K_type,unsigned long> hash;
  std::vector<K_type_pol> poly;

public:
  K_type_to_pol_table() : pool(), hash(pool), poly() {}
  template<typename F> const K_type_pol& put (K_type t, F f);
  bool is_present (const K_type& t) const { return hash.find(t)!=hash.empty; }
  const K_type_pol& lookup (const K_type& t) const;

}; // |K_type_to_pol_table|


const K_type_pol&
  branch(K_type t, repr::level cutoff,
	 K_type_to_pol_table& table, const Rep_context& rc);

} // |namespace K_repr|

} // |namespace atlas|

#endif
