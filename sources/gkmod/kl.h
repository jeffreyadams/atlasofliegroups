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

#include "hashtable.h"
#include "polynomials.h"
#include "prettyprint.h"
#include "wgraph.h"

namespace atlas {

/******** constant declarations *********************************************/

namespace kl {


std::ostream& operator<<
  (std::ostream& out,  const KLPol p);

 typedef unsigned int KLIndex;

// wrap KLPol into a class that can be used in a HashTable
class KLPolEntry : public KLPol
  {
  public:
    KLPolEntry() : KLPol() {} // default constructor builds zero polynomial
    KLPolEntry(const KLPol& p) : KLPol(p) {} // lift polynomial to this class
    explicit KLPolEntry(polynomials::Degree d) : KLPol(d) {} // represents X^d

    typedef std::vector<KLPolEntry> Pooltype; // associated storage type
    size_t hashCode(KLIndex modulus) const;   // hash function

    // equality is handled by the base type
    bool operator!=(const KLPolEntry& e) const
      { return not (static_cast<const KLPol&>(*this)==
		    static_cast<const KLPol&>(e));
      }
  };

 class HashStore : public hashtable::HashTable<KLPolEntry,KLIndex>
  {
    typedef hashtable::HashTable<KLPolEntry,KLIndex> TableType;

  public:
    HashStore() : TableType() {}

  private:
    KLIndex find(const KLPol& p) const
      {
	KLIndex i=find(p); return i==empty ? end() : i;
      }
    std::pair<KLIndex,bool> insert(const KLPol& p)
      {
	KLIndex s=size(); // will be code of KLIndex if p is new
	KLIndex i=match(p);
	return std::pair<KLIndex,bool>(i,i==s);
      }
  public:

    const KLPol& operator[] (KLIndex i) const // interpret i on our hashTable
      { return pool()[i]; }

  };


typedef HashStore KLStore;

typedef KLIndex KLPtr;

typedef std::vector<KLPtr> KLRow;


  /*!
\brief Polynomial 0, which is stored as a vector of size 0.
  */
  const KLPol Zero;

  /*!
\brief Polynomial 1.q^0.

The constructor Polynomial(d) gives 1.q^d.
  */

  const KLPol One(0);

}

/******** function declarations *********************************************/

namespace kl {

  void wGraph(wgraph::WGraph&, const KLContext&);

}

/******** type definitions **************************************************/

namespace kl {

  /*!
\brief Calculates and stores the Kazhdan-Lusztig polynomials for a
block of representations of G.
  */
class KLContext {

 protected:

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
  KLStore d_store;           // the distinct actual polynomials
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
  KLContext() {}

  KLContext(klsupport::KLSupport&);
  virtual ~KLContext() {}

// copy, assignment and swap
  KLContext(const KLContext&);

  KLContext& operator= (const KLContext&);

  void swap(KLContext&);

// accessors
  const blocks::Block& block() const {
    return d_support->block();
  }

  /*!
\brief List of the elements x_i that are primitive
with respect to y and have P_{y,x_i} not zero.
  */
  const klsupport::PrimitiveRow& primitiveRow(size_t y) const {
    return d_prim[y];
  }

  const bitset::RankFlags& descentSet(size_t y) const {
    return d_support->descentSet(y);
  }

  bool isZero(const KLPtr p) const {
    return p == d_zero;
  }

  /*!
\brief The Kazhdan-Lusztig-Vogan polynomial P_{x,y}

Note that it is returned by value, since we want to allow for the possibility
of compressed storage, in chich case we cannot return an uncompressed
polynomial by reference, but we can return it by value
  */
  const KLPol& klPol(size_t x, size_t y) const;

  /*!
\brief Returns the list of pointers to the non-zero KL polynomials
P_{y,x_i} (with x_i primitive with respect to y).
  */
  const KLRow& klRow(size_t y) const {
    return d_kl[y];
  }

  /*!
\brief Length of y as a block element.
  */
  size_t length(size_t y) const {
    return d_support->length(y);
  }

  MuCoeff mu(size_t x, size_t y) const;

  /*!
\brief List of MuData, which are pairs (x, top degree coefficient of
P_{y,x}).
  */
  const MuRow& muRow(size_t y) const {
    return d_mu[y];
  }

  /*!
\brief Returns the set of all non-zero KL polynomials for the block.
  */
  const KLStore& polStore() const {
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
  virtual void fill();


};

}

}

#endif
