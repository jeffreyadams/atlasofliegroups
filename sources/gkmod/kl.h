/*
  This is kl.h

  Class definitions and function declarations for the class KLContext.


  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright 2012 David Vogan, Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef KL_H  /* guard against multiple inclusions */
#define KL_H

#include <limits>
#include <set>

#include "atlas_types.h"

#include "klsupport.h"	// containment
#include "polynomials.h"// containment

namespace atlas {

/******** function declarations *********************************************/

namespace kl {

  void wGraph(wgraph::WGraph&, const KLContext&);

}

/******** type definitions **************************************************/

/* Namely: the definition of KLContext itself */

namespace kl {

typedef std::vector<std::pair<BlockElt,MuCoeff> > MuRow;

class KLPolEntry;
/*!
\brief Calculates and stores the Kazhdan-Lusztig polynomials for a
  block of representations of G.
  */
class KLContext
  : public klsupport::KLSupport // base is needed for full functionality
{

  BlockElt fill_limit; // all "rows" |y| with |y<fill_limit| have been computed

/*!
  \brief Entry d_prim[y] is a list of the elements x_i that are primitive
  with respect to y and have P_{y,x_i} not zero.
*/
  std::vector<PrimitiveRow> d_prim;

/*
  $d_kl[y]$ is a list of indices into |d_hashtable| of polynomials
  $P_{x_i,y}$ with $x_i=d_prim[i]$
*/
  std::vector<KLRow> d_kl;           // list of polynomial pointers

/*!
  \brief Entry d_mu[y] is a MuRow, which has parallel vectors for x and mu(x,y)
*/
 std::vector<MuRow> d_mu;           // lists of x's and their mu-coefficients

  KLStore d_store; // the distinct actual polynomials

  // the constructors will ensure that |d_store| contains 0, 1 at beginning
  enum { d_zero = 0, d_one  = 1}; // indices of polynomials 0,1 in |d_store|
  // using enum rather than |static const int| allows implicit const references

 private:  // copy and assignment are not needed, and forbidden
  KLContext(const KLContext&);
  KLContext& operator= (const KLContext&);

 public:

// constructors and destructors
  KLContext(const Block_base&); // construct initial base object

// accessors
  // construct lists of extremal respectively primitive elements for |y|
  PrimitiveRow extremalRow(BlockElt y) const;
  PrimitiveRow primitiveRow(BlockElt y) const;

  bool isZero(const KLIndex p) const { return p == d_zero; }

  // A constant reference to the Kazhdan-Lusztig-Vogan polynomial P_{x,y}
  KLPolRef klPol(BlockElt x, BlockElt y) const;

  // That polynomial in the form of an index into |polStore()==d_store|
  KLIndex KL_pol_index(BlockElt x, BlockElt y) const;

/*!
  Returns the list of pointers to the non-zero KL polynomials
  P_{x_i,y} (with x_i = d_prim[i] primitive with respect to y).
*/
  const KLRow& klRow(BlockElt y) const { return d_kl[y]; }

  MuCoeff mu(BlockElt x, BlockElt y) const; // $\mu(x,y)$

  // List of nonzero $\mu(x,y)$ for |y|, as pairs $(x,\mu(x,y))$
  const MuRow& muRow(BlockElt y) const { return d_mu[y]; }

  // List of all non-zero KL polynomials for the block, in generation order
  const KLStore& polStore() const { return d_store; }

  // get bitmap of primitive elements for row |y| with nonzero KL polynomial
  BitMap primMap (BlockElt y) const;

// manipulators

  // partial fill, up to and including the "row" of |y|
  void fill(BlockElt y, bool verbose=true);

  void fill(bool verbose=true)
  { fill(size()-1,verbose); } // simulate forbidden first default argument


  // private methods used during construction
 private:
  typedef HashTable<KLPolEntry,KLIndex> KLHash;

  //accessors
    weyl::Generator firstDirectRecursion(BlockElt y) const;
    weyl::Generator first_nice_and_real(BlockElt x,BlockElt y) const;
    std::pair<weyl::Generator,weyl::Generator>
      first_endgame_pair(BlockElt x, BlockElt y) const;
    BlockEltPair inverseCayley(size_t s, BlockElt y) const;
    std::set<BlockElt> down_set(BlockElt y) const;

    KLPolRef klPol(BlockElt x, BlockElt y,
		   KLRow::const_iterator klv,
		   PrimitiveRow::const_iterator p_begin,
		   PrimitiveRow::const_iterator p_end) const;

    // manipulators
    void silent_fill(BlockElt last_y);
    void verbose_fill(BlockElt last_y);

    // the |size_t| results serve only for statistics; caller may ignore them
    size_t fillKLRow(BlockElt y, KLHash& hash);
    void recursionRow(std::vector<KLPol> & klv,
		      const PrimitiveRow& e, BlockElt y, size_t s);
    void muCorrection(std::vector<KLPol>& klv,
		      const PrimitiveRow& e,
		      BlockElt y, size_t s);
    size_t writeRow(const std::vector<KLPol>& klv,
		    const PrimitiveRow& e, BlockElt y, KLHash& hash);
    size_t remove_zeros(const KLRow& klv,
			const PrimitiveRow& e, BlockElt y);
    void newRecursionRow(KLRow & klv,const PrimitiveRow& pr,
			 BlockElt y, KLHash& hash);
    KLPol muNewFormula(BlockElt x, BlockElt y, size_t s, const MuRow& muy);

}; //class KLContext

} // |namespace kl|

} // |namespace atlas|

#endif
