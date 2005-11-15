/*
  This is bitmap.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

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

  class BitMap {

  private:

    std::vector<unsigned long> d_map;
    unsigned long d_capacity;
    static unsigned long posBits;
    static unsigned long baseBits;
    static unsigned long baseShift;

  public:

// type definitions
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
    BitMap() {}

    explicit BitMap(unsigned long n);

    BitMap(const BitMap& b):d_map(b.d_map),d_capacity(b.d_capacity) {}

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

    unsigned long back() const;

    unsigned long capacity() const {
      return d_capacity;
    }

    bool contains(const BitMap& b) const;

    bool empty() const;

    unsigned long front() const;

    bool full() const;

    bool isMember(unsigned long n) const {
      if (n >= d_capacity)
	return false;
      return d_map[n >> baseShift] & constants::bitMask[n & posBits];
    }

    unsigned long n_th(unsigned long) const;

    unsigned long position(unsigned long) const;

    unsigned long range(unsigned long, unsigned long) const;

    size_type size() const;

// manipulators
    BitMap& operator~ ();

    BitMap& operator&= (const BitMap&);

    BitMap& operator|= (const BitMap&);

    BitMap& operator^= (const BitMap&);

    BitMap& andnot(const BitMap&);

    void fill();

    void fill(unsigned long);

    void flip(unsigned long n) {
      d_map[n >> baseShift] ^= constants::bitMask[n & posBits];
    }

    void insert(unsigned long n) {
      d_map[n >> baseShift] |= constants::bitMask[n & posBits];
    }

    template<typename I>
      void insert(const I&, const I&);

    iterator insert(iterator, unsigned long n);

    void remove(unsigned long n) {
      d_map[n >> baseShift] &= ~constants::bitMask[n & posBits];
    }

    void reset() {
      d_map.assign(d_map.size(),0ul);
    }

    void resize(unsigned long n);

    void setRange(unsigned long, unsigned long, unsigned long);

    void swap(BitMap&);
  };

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
    return d_bitAddress == i.d_bitAddress;
  }

  bool operator!= (const iterator& i) const {
    return d_bitAddress != i.d_bitAddress;
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
