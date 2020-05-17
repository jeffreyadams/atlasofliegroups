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

class KLPolEntry; // class definition will given in the implementation file

using KLColumn = std::vector<KLIndex>;
using PrimitiveColumn = std::vector<BlockElt>;
using KLHash = HashTable<KLPolEntry,KLIndex>;

/******** function declarations *********************************************/


wgraph::WGraph wGraph(const KL_table&);


/******** type definitions **************************************************/

/* Namely: the definition of KL_table itself */

struct Mu_pair
{ BlockElt x; MuCoeff coef;
  Mu_pair (BlockElt x,MuCoeff coef) : x(x), coef(coef) {}
  bool operator< (const Mu_pair& other) const { return x<other.x; }
};

using KL_column = std::vector<KLIndex>;
using Mu_column = std::vector<Mu_pair>;
using Mu_list = containers::sl_list<Mu_pair>;


/*
  |KL_table| is a class that Calculates and stores the
  Kazhdan-Lusztig-Vogan polynomials for a block of representations of $G$.
*/
class KL_table
  : public klsupport::KLSupport // base is needed for full functionality
{

  BitMap d_holes; // columns to fill; its |capacity| limits ambition to do so

// Entry |d_KL[y]| is a vector of |KLIndex|es for |P(x,y)| with |x| primitive
  std::vector<KL_column> d_KL;

// Entry |d_mu[y]| is a vector of pairs of an $x$ and corresponding |mu(x,y)|
  std::vector<Mu_column> d_mu;   // lists of $x$'s and their |mu|-coefficients

  KLStore d_store; // the distinct actual polynomials

  // the constructors will ensure that |d_store| contains 0, 1 at beginning
  enum { d_zero = 0, d_one  = 1 }; // indices of polynomials 0,1 in |d_store|
  // using enum rather than |static const int| allows implicit const references

  // copy, assignment, and swap are not needed, and not provided
  KL_table(const KL_table&) = delete;
  KL_table& operator= (const KL_table&) = delete;
  KL_table(KL_table&&) = delete;
  KL_table& operator= (KL_table&&) = delete;
  void swap(KL_table&) = delete;

 public:

// constructors and destructors
  KL_table(const Block_base&); // construct initial base object

// accessors

  BlockElt first_hole () const { return d_holes.front(); }

  // construct lists of extremal respectively primitive elements for |y|
  BitMap extremals (BlockElt y) const;
  BitMap primitives (BlockElt y) const;

  bool isZero(const KLIndex p) const { return p == d_zero; }

  // A constant reference to the Kazhdan-Lusztig-Vogan polynomial P_{x,y}
  KLPolRef KL_pol(BlockElt x, BlockElt y) const;

  // That polynomial in the form of an index into |polStore()==d_store|
  KLIndex KL_pol_index(BlockElt x, BlockElt y) const;

  MuCoeff mu(BlockElt x, BlockElt y) const; // $\mu(x,y)$


  // List of all non-zero KL polynomials for the block, in generation order
  const KLStore& pol_store() const { return d_store; }

  const KL_column& KL_data(BlockElt y) const { return d_KL[y]; }

  // List of nonzero $\mu(x,y)$ for |y|, as pairs $(x,\mu(x,y))$
  const Mu_column& mu_column(BlockElt y) const { return d_mu[y]; }

  // get bitmap of primitive elements for column |y| with nonzero KL polynomial
  BitMap prim_map (BlockElt y) const;

// manipulators

  // partial fill, up to and including the column of |y|
  void fill(BlockElt y, bool verbose=false);

  void fill(bool verbose=false)
  { if (size()>0)
     fill(size()-1,verbose); // simulate forbidden first default argument
  }

  KLHash pol_hash ();

  void swallow (KL_table&& sub, const BlockEltList& embed, KLHash& hash);

  // private methods used during construction
 private:

  //accessors
  weyl::Generator first_direct_recursion(BlockElt y) const;
  weyl::Generator first_nice_and_real(BlockElt x,BlockElt y) const;
  std::pair<weyl::Generator,weyl::Generator>
    first_endgame_pair(BlockElt x, BlockElt y) const;
  BlockEltPair inverse_Cayley(weyl::Generator s, BlockElt y) const;

  // manipulators
  void silent_fill(BlockElt last_y); // called by public |fill| when not verbose
  void verbose_fill(BlockElt last_y); // called by public |fill| when verbose

  void fill_KL_column(BlockElt y, KLHash& hash);
  void recursion_column(std::vector<KLPol> & klv,
			const BitMap& e, BlockElt y, weyl::Generator s);
  void mu_correction(std::vector<KLPol>& klv, const BitMap& e,
		     BlockElt sy, weyl::Generator s);
  void complete_primitives(std::vector<KLPol>& klv, const BitMap& e,
			   BlockElt y, KLHash& hash);
  std::vector<KLIndex>
    new_recursion_column(const BitMap& prims, BlockElt y, KLHash& hash);
  KLPol mu_new_formula
    (BlockElt x, BlockElt y, weyl::Generator s, const Mu_list& muy);

}; // |class KL_table|

// we wrap |KLPol| into a class |KLPolEntry| that can be used in a |HashTable|

/* This associates the type |KLStore| as underlying storage type to |KLPol|,
   and adds the methods |hashCode| (hash function) and |!=| (unequality), for
   use by the |HashTable| template.
 */
class KLPolEntry : public KLPol
{
public:
  // constructors
  KLPolEntry() : KLPol() {} // default constructor builds zero polynomial
  KLPolEntry(const KLPol& p) : KLPol(p) {} // lift polynomial to this class

  // members required for an Entry parameter to the HashTable template
  typedef KLStore Pooltype;		   // associated storage type
  size_t hashCode(size_t modulus) const; // hash function

  // compare polynomial with one from storage
  bool operator!=(Pooltype::const_reference e) const;

}; // |class KLPolEntry|


} // |namespace kl|

} // |namespace atlas|

#endif
