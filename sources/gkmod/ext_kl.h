/*
  This is ext_kl.h

  Copyright (C) 2013-2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef EXT_KL_H  /* guard against multiple inclusions */
#define EXT_KL_H

#include <memory> // for |std::unique_ptr|

#include "../Atlas.h"

#include "ext_block.h"
#include "polynomials.h"
#include "kl.h"

namespace atlas {

namespace ext_kl {

Pol qk_plus_1(int k);
inline Pol q_plus_1() { return qk_plus_1(1); }
Pol qk_minus_1(int k);
Pol qk_minus_q(int k);

class descent_table
{
  struct Elt_info // per block element information
  { RankFlags descents;
    RankFlags good_ascents;
  Elt_info(RankFlags d,RankFlags a): descents(d), good_ascents(a) {}
  };
  std::vector<Elt_info> info;

  std::vector<std::vector<unsigned int> > prim_index;
  std::vector<BitMap> prim_flip; // sign for |prim_index|, transposed indexing
 public:
  const ext_block::ext_block& block;

// constructors and destructors
  explicit descent_table(const ext_block::ext_block&);

// accessors

  RankFlags descent_set(BlockElt y) const { return info[y].descents; }
  RankFlags good_ascent_set(BlockElt y) const { return info[y].good_ascents; }

  bool is_descent(weyl::Generator s, BlockElt y) const
    { return descent_set(y).test(s); }

  RankFlags very_easy_set(BlockElt x, BlockElt y) const
    { return info[x].good_ascents & info[y].descents; }

  RankFlags easy_set(BlockElt x, BlockElt y) const
    { return info[y].descents-info[x].descents; }

  // index of primitive element corresponding to $x$ in row for $y$
  unsigned int x_index(BlockElt x, RankFlags desc) const
    { return prim_index[desc.to_ulong()][x]; }
  unsigned int x_index(BlockElt x, BlockElt y) const
    { return x_index(x,info[y].descents); }
  unsigned int self_index(BlockElt y) const { return x_index(y,y); }
  bool flips(BlockElt x,BlockElt y) const
    { return prim_flip[x].isMember(info[y].descents.to_ulong()); }

  BlockElt length_floor(BlockElt y) const
    { return block.length_first(block.length(y)); }

  // number of primitive elements for |descent_set(y)| of length less than |y|
  unsigned int col_size(BlockElt y) const;

  bool is_extremal (BlockElt x, RankFlags descents_y) const
    { return descent_set(x).contains(descents_y); } // easy set empty
  bool is_primitive (BlockElt x, RankFlags descents_y) const
    { return (good_ascent_set(x) & descents_y).none(); } // very easy set empty

  // set $x$ to last extremal element for $y$ strictly before $x$, or fail
  bool extr_back_up(BlockElt& x, RankFlags desc_y) const
    { while (x-->0) if (is_extremal(x,desc_y)) return true;
      return false; // now |x| has crashed through 0 and should be ignored
    }

  // set $x$ to last primitive element for $y$ strictly before $x$, or fail
  bool prim_back_up(BlockElt& x, RankFlags desc_y) const
    { while (x-->0) if (is_primitive(x,desc_y)) return true;
      return false; // now |x| has crashed through 0 and should be ignored
    }

  // set $x$ to last primitive element for $y$ strictly before $x$, or fail
  bool prim_back_up(BlockElt& x, BlockElt y) const;
  bool extr_back_up(BlockElt& x, BlockElt y) const; // same for extremal

}; // |descent_table|

struct Poly_hash_export // auxiliary to export possibly temporary hash table
{
  std::unique_ptr<ext_KL_hash_Table> own; // maybe own the temporary
  ext_KL_hash_Table& ref;

  Poly_hash_export(ext_KL_hash_Table* hash_ptr)
  : own(nullptr), ref(*hash_ptr) {}
  Poly_hash_export(IntPolEntry::Pooltype& storage)
  : own(new ext_KL_hash_Table(storage,4)), ref(*own) {}
}; // |Poly_hash_export|

class KL_table
{
  const descent_table aux;
  ext_KL_hash_Table* pol_hash; // maybe pointer to shared KL polynomial table
  std::unique_ptr<IntPolEntry::Pooltype> own; // to |storage_pool| if we own
  const IntPolEntry::Pooltype& storage_pool; // the distinct actual polynomials

  using KLColumn = std::vector<ext_kl::KLIndex>;
  std::vector<KLColumn> column; // columns are lists of polynomial pointers

  // the constructors will ensure that |storage_pool| contains 0, 1 at beginning
  enum { zero = 0, one  = 1 }; // indices of polynomials 0,1 in |storage_pool|
  // use |enum| rather than |static constxepr ext_kl::KLIndex|: avoid references

 public:
  KL_table(const ext_block::ext_block& b, ext_KL_hash_Table* poly_hash);

  size_t rank() const { return aux.block.rank(); }
  size_t size() const { return column.size(); }

  RankFlags descent_set (BlockElt y) const
  { return aux.descent_set(y); }

  ext_block::DescValue type(weyl::Generator s,BlockElt y) const
  { return aux.block.descent_type(s,y); }

  unsigned l(BlockElt y, BlockElt x) const { return aux.block.l(y,x); }

  const IntPolEntry::Pooltype& polys() const { return storage_pool;}
  std::pair<ext_kl::KLIndex,bool> KL_pol_index(BlockElt x, BlockElt y) const;

  // The twisted Kazhdan-Lusztig-Vogan polynomial P_{x,y}
  Pol P(BlockElt x, BlockElt y) const;

  bool is_extremal(BlockElt x, BlockElt y) const
  { return aux.easy_set(x,y).none(); }
  bool is_primitive(BlockElt x, BlockElt y) const
  { return aux.very_easy_set(x,y).none(); }

  // list of elements |x| such that $P(x,y)$ is nonzero, decreasing from |y|
  containers::sl_list<BlockElt> nonzero_column(BlockElt y) const;

  // coefficients in $P_{x,y}$ of $q^{(l(y/x)-i)/2}$ (use with i=1,2,3)
  int mu(short unsigned int i,BlockElt x, BlockElt y) const;

  // manipulators
  void fill_columns(BlockElt limit=0); // do all |y<limit|; if |limit==0| do all
  Poly_hash_export polynomial_hash_table ();
  void swallow (KL_table&& sub, const BlockEltList& embed);
 private:
  using PolHash = HashTable<IntPolEntry,ext_kl::KLIndex>;
  void fill_column(BlockElt y,PolHash& hash);

  // component of basis element $a_x$ in product $(T_s+1)C_{sy}$
  Pol product_comp (BlockElt x, weyl::Generator s, BlockElt sy) const;
  Pol extract_M(Pol& Q,unsigned d,unsigned defect) const;

#ifndef NDEBUG
  // compute $M(s,x,sy)$ recursively using vector |Ms| (for direct recursion)
  Pol get_M(weyl::Generator s,BlockElt x, BlockElt sy,
	    const std::vector<Pol>& Ms) const; // previous values $M(s,u,sy)$
#endif

  // variant of above for new recursion: omit term if depending on $P_{x_s,y}$
  Pol get_Mp(weyl::Generator s,BlockElt x, BlockElt y,
	     const std::vector<Pol>& Ms) const; // previous values $M(s,u,sy)$

  // look for a direct recursion and return whether possible
  bool has_direct_recursion(BlockElt y,weyl::Generator& s, BlockElt& sy) const;

  void do_new_recursion(BlockElt y,PolHash& hash);

  bool check_polys(BlockElt y) const;

}; // |ext_kl::KL_table|

// compute matrix of extended KLV polynomials evaluated at $q=-1$
void ext_KL_matrix (const StandardRepr p, const int_Matrix& delta,
		    const Rep_context& rc, // the rest is output
		    std::vector<StandardRepr>& block_list,
		    int_Matrix& P_mat,
		    IntPolEntry::Pooltype& polys);

} // |namespace ext_kl|

} // |namespace atlas|

#endif
