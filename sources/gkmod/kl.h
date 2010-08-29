/*!
\file
\brief
Class definitions and function declarations for the class KLContext.
*/
/*
  This is kl.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef KL_H  /* guard against multiple inclusions */
#define KL_H

#include <limits>

#include "kl_fwd.h"

#include "blocks_fwd.h"

#include "bitset.h"
#include "klsupport.h"
#include "polynomials.h"
#include "wgraph.h"

namespace atlas {

/******** constant declarations *********************************************/

namespace kl {

extern const KLPol Zero;


} // namespace kl

/******** function declarations *********************************************/

namespace kl {

  void wGraph(wgraph::WGraph&, const KLContext&);

}

/******** type definitions **************************************************/

/* Namely: the definition of KLContext itself */


namespace kl {

  /*!
\brief Calculates and stores the Kazhdan-Lusztig polynomials for a
block of representations of G.
  */
class KLContext
  : public klsupport::KLSupport // base is needed for full functionality
{

 protected:  // permit access of our Helper class to the data members

  /*!
\brief Records whether the KL polynomials for the block have all been computed.
  */
  enum State { KLFilled, NumStates };

  /*!
\brief Bit 0 flags whether the KL polynomials have
all been computed.
  */
  bitset::BitSet<NumStates> d_state;

  /*!
\brief Entry d_prim[y] is a list of the elements x_i that are primitive
with respect to y and have P_{y,x_i} not zero.
  */
  std::vector<PrimitiveRow> d_prim;

  /*!
\brief Entry d_kl[y] is a list of pointers to the polynomials
P_{y,x_i}, numbered as in the list d_prim[y].
  */
  std::vector<KLRow> d_kl;           // list of polynomial pointers

  /*!
\brief Entry d_mu[y] is a MuRow, which has parallel vectors for x and mu(x,y)
  */
  std::vector<MuRow> d_mu;           // lists of x's and their mu-coefficients

  /*!
\brief Set of KL polynomials.
  */
  KLStore d_store;           // the distinct actual polynomials
  /*!
\brief Pointer to the polynomial 0.
  */
  KLIndex d_zero;
  /*!
\brief Pointer to the polynomial 1.
  */
  KLIndex d_one;

// copy and swap are only for helper class
  KLContext(const KLContext&);

  void swap(KLContext&); // needed internally in helper class

 private:
  KLContext& operator= (const KLContext&); // assignment is not needed at all

 public:

// constructors and destructors
  KLContext(blocks::Block_base&); // initial base object

// accessors
  // the following two were moved here from the Helper class
  void makeExtremalRow(PrimitiveRow& e, blocks::BlockElt y) const;

  void makePrimitiveRow(PrimitiveRow& e, blocks::BlockElt y) const;

  /*!
\brief List of the elements x_i that are primitive with respect to y and have
 P_{y,x_i} NOT ZERO. This method is somewhat of a misnomer
  */
  const PrimitiveRow& primitiveRow(blocks::BlockElt y) const
    { return d_prim[y]; }

  bool isZero(const KLIndex p) const { return p == d_zero; }

  /*!
\brief The Kazhdan-Lusztig-Vogan polynomial P_{x,y}
*/
  KLPolRef klPol(blocks::BlockElt x, blocks::BlockElt y) const;

  /*!
\brief Returns the list of pointers to the non-zero KL polynomials
P_{y,x_i} (with x_i primitive with respect to y).
  */
  const KLRow& klRow(blocks::BlockElt y) const {
    return d_kl[y];
  }

  MuCoeff mu(blocks::BlockElt x, blocks::BlockElt y) const;

  /*!
\brief List of MuData, which are pairs (x, top degree coefficient of
P_{y,x}).
  */
  const MuRow& muRow(blocks::BlockElt y) const {
    return d_mu[y];
  }

  /*!
\brief Returns the set of all non-zero KL polynomials for the block.
  */
  const KLStore& polStore() const {
    return d_store;
  }

// get bitmap of primitive elements for row |y| with nonzero KL polynomial
  bitmap::BitMap primMap (blocks::BlockElt y) const;

// manipulators

  // this method used to be virtual, but that seems completely silly. MvL
  void fill();

};

} // |namespace kl|

} // |namespace atlas|

#endif
