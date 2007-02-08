/*!
\file
\brief
Class definitions and function declarations for the class KLContext.
*/
/*
  This is kl.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef KL_H  /* guard against multiple inclusions */
#define KL_H

#include <limits>
#include <set>

#include "kl_fwd.h"

#include "blocks_fwd.h"

#include "bitset.h"
#include "klsupport.h"
#include "polynomials.h"
#include "wgraph.h"

namespace atlas {

/******** constant declarations *********************************************/

namespace kl {

  /*!
\brief Polynomial 0, which is stored as a vector of size 0.
  */
  const KLPol Zero;

  /*!
\brief Polynomial 1.q^0.

The constructor Polynomial(d) gives 1.q^d.
  */

  const KLPol One(0);

  const KLCoeff UndefKLCoeff = std::numeric_limits<KLCoeff>::max();
  const KLCoeff UndefMuCoeff = std::numeric_limits<MuCoeff>::max();


typedef std::vector<KLPtr> KLRow;



} // namespace kl

/******** function declarations *********************************************/

namespace kl {

  void wGraph(wgraph::WGraph&, const KLContext&);

}

/******** type definitions **************************************************/

/* Namely: the definition of KLContext itself */


namespace kl {

 using blocks::BlockElt;

  /*!
\brief Calculates and stores the Kazhdan-Lusztig polynomials for a
block of representations of G.
  */
class KLContext {

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
\brief Pointer to the KLSupport class for this block.
  */
  klsupport::KLSupport* d_support;   // non-owned pointer

  /*!
\brief Entry d_prim[y] is a list of the elements x_i that are primitive
with respect to y and have P_{y,x_i} not zero.
  */
  std::vector<klsupport::PrimitiveRow> d_prim;

  /*!
\brief Entry d_kl[y] is a list of pointers to the polynomials
P_{y,x_i}, numbered as in the list d_prim[y].
  */
  std::vector<KLRow> d_kl;           // list of polynomial pointers

  /*!
\brief Entry d_mu[y] is a list of MuData, which are pairs (x, top
degree coefficient of P_{y,x}).
  */
  std::vector<MuRow> d_mu;           // list of mu-coefficients

  /*!
\brief Set of KL polynomials.
  */
  std::set<KLPol> d_store;           // the actual polynomials
  /*!
\brief Pointer to the polynomial 0.
  */
  KLPtr d_zero;
  /*!
\brief Pointer to the polynomial 1.
  */
  KLPtr d_one;

 public:

// constructors and destructors
  KLContext(klsupport::KLSupport&); // initial base object

  // there is no point in making the destructor virtual
  ~KLContext() {}

// copy, assignment and swap
  KLContext(const KLContext&);

  KLContext& operator= (const KLContext&);

  void swap(KLContext&);

// accessors
  const blocks::Block& block() const {
    return d_support->block();
  }

  /*!
\brief List of the elements x_i that are primitive with respect to y and have
 P_{y,x_i} NOT ZERO. This method is somewhat of a misnomer
  */
  const klsupport::PrimitiveRow& primitiveRow(BlockElt y) const {
    return d_prim[y];
  }

  const bitset::RankFlags& descentSet(BlockElt y) const {
    return d_support->descentSet(y);
  }

  bool isZero(const KLPtr p) const {
    return p == d_zero;
  }

  /*!
\brief The Kazhdan-Lusztig-Vogan polynomial P_{x,y}
*/
  const KLPol& klPol(BlockElt x, BlockElt y) const;

  /*!
\brief Returns the list of pointers to the non-zero KL polynomials
P_{y,x_i} (with x_i primitive with respect to y).
  */
  const KLRow& klRow(BlockElt y) const {
    return d_kl[y];
  }

  /*!
\brief Length of y as a block element.
  */
  size_t length(BlockElt y) const {
    return d_support->length(y);
  }

  MuCoeff mu(BlockElt x, BlockElt y) const;

  /*!
\brief List of MuData, which are pairs (x, top degree coefficient of
P_{y,x}).
  */
  const MuRow& muRow(BlockElt y) const {
    return d_mu[y];
  }

  /*!
\brief Returns the set of all non-zero KL polynomials for the block.
  */
  const std::set<KLPol>& polStore() const {
    return d_store;
  }

  /*!
\brief Rank of the group.
  */
  const size_t rank() const {
    return d_support->rank();
  }

  /*!
\brief Size of the block.
  */
  const size_t size() const {
    return d_kl.size();
  }

// manipulators

  // this method used to be virtual, but that seems completely silly. MvL
  void fill();


};

}

}

#endif
