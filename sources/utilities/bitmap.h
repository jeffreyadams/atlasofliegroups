/*!\file
  \brief Definitions and declarations for the BitMap class.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef BITMAP_H  /* guard against multiple inclusions */
#define BITMAP_H

#include <vector>

#include "bitmap_fwd.h"

#include "constants.h"

/******** type definitions **************************************************/

namespace atlas {

namespace bitmap {
  /*!
  \brief Container of a large (more than twice the machine word size)
  set of bits.

  From the point of view of a user of the class, a BitMap should be
  seen as a container of _unsigned long_, not bits: these unsigned
  longs are the addresses of the set bits.  When the class is used for
  example to flag the noncompact roots in a set of roots, it is most
  convenient to think of it as containing the numbers of the
  noncompact roots (on a fixed list of all roots).  The class obeys
  the semantics of a Forward Container (from the C++ standard
  library).  Its "size" as a container is the number of unsigned long
  that it "contains"; that is, the number of set bits in the bitmap.

  The basic data is in d_map, a vector of unsigned long integers.  Each
  of these integers is a "chunk" (of size longBits, presumably the machine
  word length) of the bit map.  The capacity (in bits) of the BitMap is
  d_capacity; the size of the vector d_map is d_capacity/longBits
  (plus one if there is a remainder in the division).

  We wish to provide bit-address access to this map; for this purpose
  we use the reference trick from vector<bool>. Also we wish to define
  an iterator class, which traverses the _set_ bits of the bitmap; so
  that for instance, b.begin() would be a reference to the first set
  bit. Dereferencing the iterator yields its bit-address.
  */
  class BitMap {

  private:

    std::vector<unsigned long> d_map;
    unsigned long d_capacity;
    static unsigned long posBits;
    static unsigned long baseBits;
    static unsigned long baseShift;

  public:

// type definitions

  /*!
  type for a component of the vector d_map holding the BitMap
  */
    typedef unsigned long value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef ptrdiff_t difference_type;
    typedef unsigned long size_type;

// iterators
    class iterator;
    typedef iterator const_iterator;
    iterator begin() const;
    iterator end() const;
    iterator pos(unsigned long) const;

// constructors and destructors
    BitMap() : d_map(), d_capacity(0) {}

   /*! \brief
       Constructs a zero-initialized bitmap with a capacity of n bits.

      Notice that the size of the vector |d_map| exceeds |n >> baseShift| by
      one unless |longBits| exactly divides |n|.
   */
    explicit BitMap(unsigned long n)
      : d_map((n+posBits)>>baseShift,0), d_capacity(n)
      {}

    BitMap(const BitMap& b) : d_map(b.d_map), d_capacity(b.d_capacity) {}

    template <typename I, typename J>
      BitMap(const I&, const I&, const J&, const J&);

    ~BitMap() {}

// assignment
    BitMap& operator= (const BitMap&);

// accessors
    bool operator< (const BitMap& b) const {
      return d_map < b.d_map;
    }

    bool operator== (const BitMap& b) const {
      return (d_map == b.d_map);
    }

    // decrement |n| (at least once) until it points to a member
    // return value indicates whether this succeeded
    bool back_up(unsigned long& n) const;

    /*!
    Number of bits in use in the bitmap. This is the capacity
    of the BitMap as a standard library container, not d_map.size(), which is
    approximately longBits times smaller.
    */
    unsigned long capacity() const {
      return d_capacity;
    }

    bool contains(const BitMap& b) const;

    bool empty() const;

    unsigned long front() const;

    bool full() const;

    /*!
  Tests whether bit n in the bitmap is set; that is, whether element n
  is a member of the set.
    */
    bool isMember(unsigned long n) const {
      if (n >= d_capacity)
	return false;
      return d_map[n >> baseShift] & constants::bitMask[n & posBits];
    }

    unsigned long n_th(unsigned long) const;

    unsigned long position(unsigned long) const;

    // get a range of bits as unsigned long value; see bitmap.ccp for details
    unsigned long range(unsigned long first, unsigned long number) const;

    size_type size() const; // the number of bits that are set in the bitmap

// manipulators

    // WARNING: the following looks like an accessor, but complements |*this|
    BitMap& operator~ ();

    bool operator&= (const BitMap&);

    BitMap& operator|= (const BitMap&);

    BitMap& operator^= (const BitMap&);

    bool andnot(const BitMap&);

    void fill();

    void fill(unsigned long);

    void flip(unsigned long n) {
      d_map[n >> baseShift] ^= constants::bitMask[n & posBits];
    }
    /*!
    Puts a 1 in position n of the bitmap (that is, inserts an element of
    the set).
    */
    void insert(unsigned long n) {
      d_map[n >> baseShift] |= constants::bitMask[n & posBits];
    }

    /*!
    Puts a 0 in position n of the bitmap (that is, removes an element of
    the set).
    */
    void remove(unsigned long n) {
      d_map[n >> baseShift] &= ~constants::bitMask[n & posBits];
    }

    void set_to(unsigned long n,bool b) {
      if (b) insert(n); else remove(n);
    }

    void set_mod2(unsigned long n, unsigned long v) {
      set_to(n,(v&1)!=0);
    }


    template<typename I>
      void insert(const I&, const I&);

    // the first argument is ignored, just serves to get overloading OK
    iterator insert(iterator, unsigned long n);

    // clear all bits, but do not change the capacity
    void reset() {
      d_map.assign(d_map.size(),0ul);
    }

    // this was called |resize|, but sets |capacity()|, whence the new name
    void set_capacity(unsigned long n);
    void resize(unsigned long n) { set_capacity(n); } // temp. compatibilty

    void setRange(unsigned long, unsigned long, unsigned long);

    void swap(BitMap&);
  };


  /*!
  \brief Traverses the set bits of a BitMap.

  Because of the nature of a bitmap (as a collection of addresses of
  set bits), only constant iterators make sense; one cannot "change
  the value" of an element at a given position because the value _is_
  the position. The most delicate operation is the increment, which
  has to find the position of the next set bit, while avoiding falling
  off the bitmap if there is no such. Because of this we had to pass
  the data for the end of the bitmap into the iterator.
  */
class BitMap::iterator { // is really a const_iterator

 private:

  std::vector<unsigned long>::const_iterator d_chunk;
  unsigned long d_bitAddress;
  unsigned long d_capacity;

 public:

// associated types
  typedef std::forward_iterator_tag iterator_category;
  typedef unsigned long value_type;
  typedef ptrdiff_t difference_type;
  typedef const value_type* pointer;
  typedef const value_type& reference;

// constructors and destructors
  iterator() {}

  iterator(const iterator &j):d_chunk(j.d_chunk), d_bitAddress(j.d_bitAddress),
    d_capacity(j.d_capacity) {}

  iterator(const std::vector<unsigned long>::const_iterator& p,
	   unsigned long n, unsigned long c)
    :d_chunk(p), d_bitAddress(n), d_capacity(c) {}

  ~iterator() {}

// assignment
  iterator& operator= (const iterator& i);

// accessors
  bool operator== (const iterator& i) const {
    return d_bitAddress == i.d_bitAddress; // note that d_chunk is ignored!
  }

  bool operator!= (const iterator& i) const {
    return d_bitAddress != i.d_bitAddress; // note that d_chunk is ignored!
  }

  bool operator() () const {
    return d_bitAddress != d_capacity;
  }

  const value_type& operator* () const {
    return d_bitAddress;
  }

// manipulators
  iterator& operator++ ();

  iterator operator++ (int);
};

}

}

#include "bitmap_def.h"

#endif
