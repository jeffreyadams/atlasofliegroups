/*
  This is kl.h

  Class definitions and function declarations for the class |KL_table|.


  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2012 David Vogan
  Copyright (C) 2006-2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef KL_H  /* guard against multiple inclusions */
#define KL_H

#include <limits>
#include <set>

#include "../Atlas.h"

#include "bitmap.h"
#include "klsupport.h"	// containment
#include "polynomials.h"// containment

namespace atlas {

namespace kl {

/******** function declarations *********************************************/


wgraph::WGraph wGraph(const KL_table&);


/******** type definitions **************************************************/

/* Namely: the definition of KL_table itself */

struct KL_pair
{ BlockElt x; KLIndex P;
  KL_pair (BlockElt x=UndefBlock, KLIndex P=0) : x(x), P(P) {}
  bool operator< (const KL_pair& other) const { return x<other.x; }
};
struct Mu_pair
{ BlockElt x; MuCoeff coef;
  Mu_pair (BlockElt x,MuCoeff coef) : x(x), coef(coef) {}
  bool operator< (const Mu_pair& other) const { return x<other.x; }
};

using KL_column = std::vector<KL_pair>;
using Mu_column = std::vector<Mu_pair>;
using Mu_list = containers::sl_list<Mu_pair>;
using KLStore = PosPolEntry::Pooltype;
using KLPolRef = KLStore::const_reference;

struct Poly_hash_export // auxiliary to export possibly temporary hash table
{
  std::unique_ptr<KL_hash_Table> own; // maybe own the temporary
  KL_hash_Table& ref;

  Poly_hash_export(KL_hash_Table* hash_ptr)
  : own(nullptr), ref(*hash_ptr) {}
  Poly_hash_export(KLStore& storage)
  : own(new KL_hash_Table(storage,4)), ref(*own) {}
}; // |Poly_hash_export|

/*
  |KL_table| is a class that Calculates and stores the
  Kazhdan-Lusztig-Vogan polynomials for a block of representations of $G$.
*/
class KL_table
  : public klsupport::KLSupport // base is needed for full functionality
{

  BitMap d_holes; // columns to fill; its |capacity| limits ambition to do so

// Entry |d_KL[y]| is a sorted vector of pairs |(x,P(x,y))|
  std::vector<KL_column> d_KL;

// Entry |d_mu[y]| is a vector of pairs of an $x$ and corresponding |mu(x,y)|
  std::vector<Mu_column> d_mu;   // lists of $x$'s and their |mu|-coefficients

  KL_hash_Table* const pol_hash; // maybe pointer to shared KL polynomial table
  std::unique_ptr<KLStore> own;  // holds object for |storage_pool| if we own

  const KLStore& storage_pool;

  // the constructors will ensure that |d_store| contains 0, 1 at beginning
  enum { zero = 0, one  = 1 }; // indices of polynomials 0,1 in |d_store|
  // use |enum| rather than |static constxepr KLIndex|: avoid any references

  // copy, assignment, and swap are not needed, and not provided
  KL_table(const KL_table&) = delete;
  KL_table& operator= (const KL_table&) = delete;
  KL_table(KL_table&&) = delete;
  KL_table& operator= (KL_table&&) = delete;
  void swap(KL_table&) = delete;

 public:

// constructors and destructors
  KL_table(const Block_base&, KL_hash_Table* pol_hash=nullptr);

// accessors

  // the first hole, or |block().size()| if none are left
  BlockElt first_hole () const { return d_holes.front(); }

  BlockElt length_floor(BlockElt y) const
  { return length_less(length(y)); }

  // construct lists of primitive elements for |y|
  BitMap primitives (BlockElt y) const;

  bool isZero(const KLIndex p) const { return p == zero; }

  // A constant reference to the Kazhdan-Lusztig-Vogan polynomial P_{x,y}
  KLPolRef KL_pol(BlockElt x, BlockElt y) const;

  // That polynomial in the form of an index into |polStore()==d_store|
  KLIndex KL_pol_index(BlockElt x, BlockElt y) const;

  MuCoeff mu(BlockElt x, BlockElt y) const; // $\mu(x,y)$


  // List of all non-zero KL polynomials for the block, in generation order
  const KLStore& pol_store() const { return storage_pool; }

  const KL_column& KL_data(BlockElt y) const { return d_KL[y]; }

  // List of nonzero $\mu(x,y)$ for |y|, as pairs $(x,\mu(x,y))$
  const Mu_column& mu_column(BlockElt y) const { return d_mu[y]; }

  // get bitmap of primitive elements for column |y| with nonzero KL polynomial
  BitMap prim_map (BlockElt y) const;

// manipulators

  // partial fill, up to column |limit| exclusive; fill all if |limit==0|
  void fill (BlockElt limit=0, bool verbose=false);

  Poly_hash_export polynomial_hash_table ();

  void swallow (KL_table&& sub, const BlockEltList& embed, KL_hash_Table& hash);

  // private methods used during construction
 private:

  //accessors
  weyl::Generator first_direct_recursion(BlockElt y) const;
  weyl::Generator first_nice_and_real(BlockElt x,BlockElt y) const;
  std::pair<weyl::Generator,weyl::Generator>
    first_endgame_pair(BlockElt x, BlockElt y) const;
  BlockEltPair inverse_Cayley(weyl::Generator s, BlockElt y) const;

  // like |KL_pol(x,y)|, but take value from |col|, when |d_KL[y]| not yet ready
  // this assumes either |x| primitive for |y|, or |col| has size one more than
  // that of the block, with entries |One| at |y| and |Zero| beyond
  KLPol lookup(BlockElt x, RankFlags desc_y, const std::vector<KLPol>& col) const
  { return col[primitivize(x,desc_y)]; }

  // manipulators
  void silent_fill(BlockElt limit); // called by public |fill| when not verbose
  void verbose_fill(BlockElt limit); // called by public |fill| when verbose

  // fill column for |y| in the KL-table, all previous ones having been filled
  void fill_KL_column(std::vector<KLPol>& klv, BlockElt y);
  void recursion_column(BlockElt y, weyl::Generator s,
			std::vector<KLPol>& klv);
  void mu_correction(const BlockEltList& extremals,
		     RankFlags desc_y, BlockElt sy, weyl::Generator s,
		     std::vector<KLPol>& klv);
  void complete_primitives(std::vector<KLPol>& klv, BlockElt y);
  void new_recursion_column(std::vector<KLPol>& klv, BlockElt y);
  KLPol mu_new_formula
    (BlockElt x, BlockElt y, weyl::Generator s, const Mu_list& muy);

}; // |class KL_table|


} // |namespace kl|

} // |namespace atlas|

#endif
