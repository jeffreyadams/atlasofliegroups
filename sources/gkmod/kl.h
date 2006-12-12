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

#include "hashtable-stats.h"
#include "polynomials.h"
#include "prettyprint.h"
#include "wgraph.h"

namespace atlas {


std::ostream& operator<<  (std::ostream& out, const kl::KLPol& p);


namespace kl {


class KLPool; // predeclare the storage class to be used via the hash table


// wrap |KLPol| into a class that can be used in a HashTable

/* This associates the type |KLPool| as underlying storage type to |KLPol|,
   and adds the methods |hashCode| (hash function) and |!=| (unequality), for
   use by the |HashTable| template.
 */

class KLPolEntry : public KLPol
  {
  public:
    // constructors
    KLPolEntry() : KLPol() {} // default constructor builds zero polynomial
    KLPolEntry(const KLPol& p) : KLPol(p) {} // lift polynomial to this class

    // the following constructor is (only) needed during rehashing
    // since the hash function is only defined for KLPolEntry objects
    explicit KLPolEntry(KLPolRef p); // extract polynomial from a PolRef

    // members required for an Entry parameter to the HashTable template
    typedef KLPool Pooltype;		    // associated storage type
    size_t hashCode(KLIndex modulus) const; // hash function

    bool operator!=(KLPolRef e) const;  // compare pol with one from storage

  };


// storage for KL polynomials

/*
   The idea is to store only sequences of coefficients in a huge vector
   |pool|, and to provide minimal overhead (the |index| vector) to be able to
   access individual polynomials. Polynomials are extracted in the form of a
   |KLPolRef| object, which class replaces the type |const KLPol&| that is
   impossible to produce without having an actual |KLPol| value in memory.
   This is why we typedef |const_reference| to be equal to |KLPolRef|.
   Insertion is done by the method |push_back|, and extraction by
   |operator[]|, mimicking the behaviour of |std::vector<KLPol>| except for
   the return type from |operator[]|. Other methods are |size|, needed by the
   hash table code to know the sequence number to assign to a new polynomial
   sent to the storage, and |swap| for the |swap| method of hash tables. The
   method |mem_size| gives the number of bytes used in storage, for
   statistical convenience.

*/

class KLPool
  {
    struct IndexType
    {
      size_t pool_index;
      IndexType(size_t i) : pool_index(i) {}
    };

    std::vector<KLCoeff> pool;
    std::vector<IndexType> index;

  public:
    typedef KLPolRef const_reference; // used by HashTable

    // constructor;
    KLPool() : pool(),index(1,IndexType(0)) { } // index is always 1 ahead

    // accessors
    const_reference operator[] (KLIndex i) const; // select polynomial by nr

    size_t size() const { return index.size()-1; } // number of entries
    size_t mem_size() const                        // memory footprint
      { return sizeof(KLPool)
	  +pool.size()*sizeof(KLCoeff)
	  +index.size()*sizeof(IndexType);
      }

    // manipulators
    void push_back(const KLPol&);

    void swap(KLPool& other)
      { pool.swap(other.pool); index.swap(other.index); }
  }; // class KLPool



typedef hashtable::HashTable<KLPolEntry,KLIndex> KLStore;

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
  KLPolRef klPol(size_t x, size_t y) const;

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
