/*
  This is ext_kl.h

  Copyright (C) 2013 Marc vna Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef EXT_KL_H  /* guard against multiple inclusions */
#define EXT_KL_H

#include "ext_block.h"
#include "atlas_types.h"
#include "polynomials.h"

namespace atlas {

namespace ext_kl {

class descent_table
{
  std::vector<RankFlags> descents; // sets of descents ($\tau$-invariants)
  std::vector<RankFlags> good_ascents; // ascents usable in primitivization

  std::vector<std::vector<unsigned int> > prim_index;

 public:

// constructors and destructors
  descent_table(const ext_block::extended_block&);

// accessors
  // index of primitive element corresponding to $x$ in row for $y$
  unsigned int x_index(BlockElt x, BlockElt y) const
  { return prim_index[descents[y].to_ulong()][x]; }

  unsigned int self_index(BlockElt y) const { return x_index(y,y); }
  unsigned int col_size(BlockElt y) const;

  // set $x$ to last primitive element for $y$ strictly before $x$, or fail
  bool back_up(BlockElt& x, BlockElt y) const;
}; // |descent_table|

class KL_table
{
  const ext_block::extended_block& block;
  const descent_table& aux;
  kl::KLStore& storage_pool; // the distinct actual polynomials

  std::vector<kl::KLRow> column; // columns are lists of polynomial pointers

 public:
   KL_table(const ext_block::extended_block& b, kl::KLStore& pool)
    : block(b), aux(b), storage_pool(pool), column() {}

  // A constant reference to the Kazhdan-Lusztig-Vogan polynomial P_{x,y}
  kl::KLPolRef P(BlockElt x, BlockElt y) const
  { return storage_pool[KL_pol_index(x,y)]; }
  kl::KLIndex KL_pol_index(BlockElt x, BlockElt y) const;

  // coefficients of P_{x,y} of $q^{(l(y/x)-i)/2}$ (use with i=1,2,3)
  kl::KLCoeff mu(int i,BlockElt x, BlockElt y) const;

  kl::KLPol m(weyl::Generator s,BlockElt x, BlockElt y) const;

  // That polynomial in the form of an index into |storage_pool|

  // manipulator
  void fill_columns(BlockElt y=0)
  { if (y==0)
      y=block.size();
    column.reserve(y);
    while (column.size()<y)
      fill_next_column();
  }
 private:
  void fill_next_column();
  BlockEltList mu1top(weyl::Generator s,BlockElt x, BlockElt y) const;
  BlockEltList mu1bot(weyl::Generator s,BlockElt x, BlockElt y) const;

}; // |KL_table|

} // |namespace ext_kl|

} // |namespace atlas|

#endif
