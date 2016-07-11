/*
  This is ext_kl.h

  Copyright (C) 2013 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef EXT_KL_H  /* guard against multiple inclusions */
#define EXT_KL_H

#include "ext_block.h"
#include "atlas_types.h"
#include "polynomials.h"
#include "kl.h"

namespace atlas {

namespace ext_kl {

typedef Polynomial<int> Pol;
typedef const Pol& PolRef;

Pol qk_plus_1(int k);
inline Pol q_plus_1() { return qk_plus_1(1); }
Pol qk_minus_1(int k);
Pol qk_minus_q(int k);

class descent_table
{
  std::vector<RankFlags> descents; // sets of descents ($\tau$-invariants)
  std::vector<RankFlags> good_ascents; // ascents usable in primitivization

  std::vector<std::vector<unsigned int> > prim_index;
  std::vector<BitMap> prim_flip; // sign for |prim_index|, transposed indexing
 public:
  const ext_block::extended_block& block;

// constructors and destructors
  explicit descent_table(const ext_block::extended_block&);

// accessors

  RankFlags descent_set(BlockElt y) const { return descents[y]; }

  // index of primitive element corresponding to $x$ in row for $y$
  unsigned int x_index(BlockElt x, BlockElt y) const
  { return prim_index[descents[y].to_ulong()][x]; }
  unsigned int self_index(BlockElt y) const { return x_index(y,y); }
  bool flips(BlockElt x,BlockElt y) const
  { return prim_flip[x].isMember(descents[y].to_ulong()); }

  BlockElt length_floor(BlockElt y) const
  { return block.length_first(block.length(y)); }
  unsigned int col_size(BlockElt y) const;

  RankFlags very_easy_set(BlockElt x, BlockElt y) const
  { return good_ascents[x]&descents[y]; }

  RankFlags easy_set(BlockElt x, BlockElt y) const
  { return descents[y]-descents[x]; }

  // set $x$ to last primitive element for $y$ strictly before $x$, or fail
  bool prim_back_up(BlockElt& x, BlockElt y) const;
  bool extr_back_up(BlockElt& x, BlockElt y) const; // same for extremal

}; // |descent_table|
class PolEntry; // class definition will be given in the implementation file

class KL_table
{
  const descent_table aux;
  std::vector<Pol>& storage_pool; // the distinct actual polynomials

  std::vector<kl::KLRow> column; // columns are lists of polynomial pointers

  //TEMPORARY
  kl::KLContext untwisted;

 public:
  KL_table(const ext_block::extended_block& b, std::vector<Pol>& pool);

  size_t rank() const { return aux.block.rank(); }
  size_t size() const { return column.size(); }

  ext_block::DescValue type(weyl::Generator s,BlockElt y) const
  { return aux.block.descent_type(s,y); }

  // The twisted Kazhdan-Lusztig-Vogan polynomial P_{x,y}
  Pol P(BlockElt x, BlockElt y) const;

  bool extremal(BlockElt x, BlockElt y) const
  { return aux.descent_set(x).contains(aux.descent_set(y)); }


  // coefficients of P_{x,y} of $q^{(l(y/x)-i)/2}$ (use with i=1,2,3)
  int mu(int i,BlockElt x, BlockElt y) const;

  // manipulator
  void fill_columns(BlockElt y=0);
 private:
  typedef HashTable<PolEntry,kl::KLIndex> PolHash;
  void fill_next_column(PolHash& hash);

  // component of basis element $a_x$ in product $(T_s+1)C_{sy}$
  Pol product_comp (BlockElt x, weyl::Generator s, BlockElt sy) const;
  Pol extract_M(Pol& Q,unsigned d,unsigned defect) const;

  // compute $M(s,x,sy)$ recursively using vector |Ms| (for direct recursion)
  Pol get_M(weyl::Generator s,BlockElt x, BlockElt sy,
	    const std::vector<Pol>& Ms) const; // previous values $M(s,u,sy)$

  // variant of above for new recursion: omit term if depending on $P_{x_s,y}$
  Pol get_Mp(weyl::Generator s,BlockElt x, BlockElt y,
	     const std::vector<Pol>& Ms) const; // previous values $M(s,u,sy)$

  // look for a direct recursion and return whether possible;
  // if possible also get contributions from $c_s*a_y$ into |out|
  bool has_direct_recursion(BlockElt y,weyl::Generator& s, BlockElt& sy) const;

  void do_new_recursion(BlockElt y,PolHash& hash);

  bool check_polys(BlockElt y) const;

}; // |KL_table|

} // |namespace ext_kl|

} // |namespace atlas|

#endif
