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
    // a packed structure of size 1
    struct packed_byte
    {
      unsigned char b;

      unsigned int degree() const    { return b&0x1F; } // lower 5 bits
      unsigned int valuation() const { return b>>5; }   // upper 3 bits

      packed_byte(unsigned int d, unsigned int v) : b(d+ (v<<5)) {}
      packed_byte() : b(0) {} // to initialise array members
    };

    // some parameters for storage layout
    static const int int_bits  =std::numeric_limits<unsigned int>::digits;
    static const size_t low_mask; // bit mask for lower order unsigned int

    static const unsigned int deg_limit=32; // (hard) degree limit
    static const unsigned int val_limit=8;  // (soft) valuation limit
    static const int group_size=16; // number of indices packed in a group


/* The following kludge is necessary to repair the broken operation of bit
   shift operations with shift amount exceeding the integer's width.
   Apparently only the final bits of the shift amount are used, so that it is
   interpreted modulo the integer width rather than setting the result to 0
*/

#if ~0ul == ~0u
    static inline unsigned int high_order_int(size_t)   { return 0; }
    static inline size_t set_high_order(unsigned int)   { return 0; }
#else
    static inline unsigned int high_order_int(size_t x)
      { return x>>int_bits; }
    static inline size_t set_high_order(unsigned int x)
      { return size_t(x)<<int_bits; }
#endif


/* The following structure collects information about a group of |group_size|
   consecutive polynomials. The address of the very first coefficient is
   recorded in 32+5=37 bits, and for all but the last polynomial the degree
   |deg<32| and a valuation |val<8| are stored in a |packed_byte|, which
   together detemine the number |1+deg-val| of stored coefficients. The
   valuation of the last polynomial of the group will be stored in the
   |pool_index_high| field of the next |IndexType| structure, and its number
   of coefficients is implicitly determined by the number of coefficients
   remaining between the first coefficient of the current group and that of
   the next one.
*/

    struct IndexType
    {
      // non-data members

      // data members
      unsigned int pool_index_low; // lower order 32 bits of index into pool
      packed_byte pool_index_high; // high order 5 bits of above, +3 bits val
      packed_byte deg_val[group_size-1]; // degree/valuation of polynomials

      IndexType(size_t i, unsigned int v)
	: pool_index_low(i&low_mask)
        , pool_index_high(high_order_int(i),v) {}
    };

    std::vector<KLCoeff> pool;
    std::vector<IndexType> index;
    unsigned int last_index_size; // nr of bytes of last index structe in use

    size_t savings; // gather statistics about savings by using valuations

  public:
    typedef KLPolRef const_reference; // used by HashTable

    // constructor;
    KLPool() : pool(),index(1,IndexType(0,0)), // index is always 1 ahead
      last_index_size(0),
      savings(0) {}

    ~KLPool(); // may print statictics

    // accessors
    const_reference operator[] (KLIndex i) const; // select polynomial by nr

    size_t size() const       // number of entries
      { return group_size*(index.size()-1)+last_index_size; }
    size_t mem_size() const   // memory footprint
      { return sizeof(KLPool)
	  +pool.size()*sizeof(KLCoeff)
	  +index.size()*sizeof(IndexType);
      }

    // manipulators
    void push_back(const KLPol&);

    void swap(KLPool& other)
      {
	pool.swap(other.pool);
	index.swap(other.index);
	std::swap(last_index_size,other.last_index_size);
	std::swap(savings, other.savings);
      }
  }; // class KLPool



typedef hashtable::HashTable<KLPolEntry,KLIndex> KLStore;

typedef KLIndex KLPtr;

typedef std::vector<KLPtr> KLRow;



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
