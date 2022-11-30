/*
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006,2017,2022 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Definitions and declarations for the BitMap class.

#ifndef BITMAP_H  /* guard against multiple inclusions */
#define BITMAP_H

#include <vector>
#include <cassert>

#include "bitmap_fwd.h"

#include "constants.h"

/******** type definitions **************************************************/

namespace atlas {

namespace bitmap {
/*
  Container of a large (more than twice the machine word size) set of bits.

  From the point of view of a user of the class, a BitMap should be seen as a
  container of |size_t| values, not bits: these unsigned values are the
  addresses of the set bits. When the class is used for example to flag the
  noncompact roots in a set of roots, it is most convenient to think of it as
  containing the numbers of the noncompact roots (on a fixed list of all roots).
  The class obeys the semantics of a Forward Container (from the C++ standard
  library). Its "size" as a container is the number of |size_t| values that it
  "contains"; that is, the number of bits set to 1 in the bitmap.

  The basic data is in |d_map|, a vector of |unsigned long| integers. Each of
  these integers is a "chunk" (of size |longBits|, presumably the machine word
  length) of the bit map. The capacity (in bits) of the |BitMap| is
  |d_capacity|; the size of the vector |d_map| is |d_capacity/longBits| (plus
  one if there is a remainder in the division).

  We provide bit-addressed access to this map, in other words given a value it
  is fast to test or change its membership staturs. Also we wish to define an
  iterator class, which traverses the _set_ bits of the bitmap; so that for
  instance, b.begin() would access the first set bit. Dereferencing the iterator
  yields its bit-address, a |size_t| value.
*/
class BitMap
{
  size_t d_capacity;
  std::vector<unsigned long> d_map;
  static unsigned long posBits;
  static unsigned long baseBits;
  static unsigned long baseShift;

 public:

// type definitions

  using value_type = size_t; // type of values virtually present in a |BitMap|
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef ptrdiff_t difference_type;
  using size_type = size_t; // container capacity same type as |value_type|

// iterators
  class iterator;
  typedef iterator const_iterator;
  iterator begin() const;
  iterator end() const;

// constructors and destructors
 BitMap() : d_capacity(0), d_map() {} // create a bitmap without capacity

  /*
    Construct a zero-initialized bitmap with a capacity of n bits.

    Note that the size of the vector |d_map| exceeds |n >> baseShift| by
    one, unless |longBits| exactly divides |n|.
  */
  explicit BitMap(size_t n)
    : d_capacity(n), d_map((n+posBits)>>baseShift,0)
  {}

  // convert range defined by iterators into a BitMap
  template <typename I>
    BitMap(size_t n, const I& first, const I& last);

  //g an easier version for a full vector
  template <typename U> // unsigned integral type
  BitMap(size_t n,const std::vector<U>& v)
    : d_capacity(n), d_map((n+posBits)>>baseShift)
  {
    for (size_t i=0; i<v.size(); ++i)
      insert(v[i]);
  }

  // Set of offsets into [first,last[ of values (also) found in [fsub,lsub[
  template <typename I, typename J>
    BitMap(const I& first, const I& last, const J& fsub, const J& lsub);


// assignment
  BitMap& operator= (const BitMap&);

// accessors
  /*
   Number of bits in use in the bitmap. This is the capacity
   of the BitMap as a standard library container, not d_map.size(), which is
   approximately longBits times smaller.
  */
  size_t capacity() const { return d_capacity; }
  size_type size() const; // the number of bits that are set in the bitmap

  bool operator< (const BitMap& b) const { return d_map < b.d_map; }
  bool operator== (const BitMap& b) const { return d_map == b.d_map; }
  bool operator!=(const BitMap& b) const { return d_map != b.d_map; }

  bool empty() const; // whether |size()==0|
  bool full() const; // whether |size()==capacity()|
  bool none() const { return empty(); } // naming compatibility with |BitSet|
  bool any() const { return not empty(); }

  /*
    Tests whether bit n in the bitmap is set; that is, whether element n
    is a member of the set.
  */
  bool isMember(size_t n) const
  {
    if (n >= d_capacity)
      return false;
    return (d_map[n >> baseShift] & constants::bitMask[n&posBits])!=0;
  }

  // Whether all elements of |b| satisfy |isMember|
  bool contains(const BitMap& b) const;

  // Whether none of the elements of |b| satisfy |isMember|
  bool disjoint(const BitMap& b) const;

  // Value at index |n| if viewed as list of |size_t| values
  size_t n_th(size_t n) const;

  // Number of values |<n| present (set) in the bitmap
  size_t position(size_t n) const;

  // First value for which |isMember| holds, or |capacity()| if none exist
  size_t front() const;

  // decrement |n| (at least once) until it points to a member
  // return value indicates whether this succeeded
  bool back_up(size_t& n) const;

  // get a range of bits as unsigned long value; see bitmap.ccp for details
  unsigned long range(size_t first, unsigned number) const;

  BitMap operator~ () const { return BitMap(*this).take_complement(); }

  BitMap operator& (const BitMap& other) && { operator&=(other); return *this; }

  BitMap operator| (const BitMap& other) && { return operator|=(other); }

  BitMap operator^ (const BitMap& other) && { return operator^=(other); }

  // bitwise cannot exploit |BitMap&& other|, as |other| allowed to be shorter

  BitMap operator& (const BitMap& other) const &
  { BitMap result(*this); result&=other; return result; }

  BitMap operator| (const BitMap& other) const &
  { BitMap result(*this); return result|=other; }

  BitMap operator^ (const BitMap& other) const &
  { BitMap result(*this); return result^=other; }


// manipulators

  BitMap& take_complement ();

 // these operators allow argument to have less capacity than |*this|
  bool operator&= (const BitMap&); // return whether result is non-empty

  BitMap& operator|= (const BitMap&);

  BitMap& operator^= (const BitMap&);

  BitMap& andnot(const BitMap& b); // remove bits of |b|

  BitMap& operator>>= (size_t delta); // shift right (decrease)
  BitMap& operator<<= (size_t delta); // shift left (increase)

  /*
    Set the bit at position |n| (that is, inserts the value |n| into the set);
    this makes |isMember(n)| hold.
  */
  void insert(size_t n)
  {
    assert(n<d_capacity);
    d_map[n >> baseShift] |= constants::bitMask[n & posBits];
  }

  /*
    Clear the bit at position |n| (that is, removes that element of the set);
    this makes |isMember(n)| false.
  */
  void remove(size_t n)
  {
    assert(n<d_capacity);
    d_map[n >> baseShift] &= ~constants::bitMask[n & posBits];
  }

  bool set_to(size_t n,bool b)
  { if (b) insert(n); else remove(n); return b; }

  bool set_mod2(size_t n, unsigned v) { return set_to(n,(v&1)!=0); }

  void flip(size_t n)
  {
    assert(n<d_capacity);
    d_map[n >> baseShift] ^= constants::bitMask[n & posBits];
  }


 void fill(); // set all bits
 void reset() { d_map.assign(d_map.size(),0ul);  } // clear all bits
 void fill(size_t start,size_t stop); // set consecutive range of bits
 void clear(size_t start,size_t stop); // clear consecutive range of bits
 // insert a range of values (which need not be listed increasingly)
 template<typename I>
   void insert(I, I);

 // this was called |resize|, but sets |capacity()|, whence the new name
 void set_capacity(size_t n); // any new bits will start out cleared
 void extend_capacity(bool b); // extend capacity by |1|, adding member if |b|

 // set an interval of bits from those (least significant ones) of source
 void setRange(size_t start, unsigned amount, unsigned long source);

 void swap(BitMap&);

}; // |class BitMap|

size_t last (const BitMap& set); // final element of non-empty |set|

  /* Iterator to traverse the \emph{set} bits (members present) of a BitMap.

  Because of the nature of a bitmap, only constant iterators make sense (just
  like for iterators into std::set); one cannot "change the value" of an
  element at a given position because the position is determined by the value
  (values are always increasing; in fact a small value-change could be
  realised by swapping a set bit with neighboring unset bits, but that is of
  course not what a non-constant iterator should allow doing). However, these
  iterators do allow changing the underlying bitmap during traversal, even
  though such changes cannot be performed using only the iterator itself (this
  is unlike iterators over a |bitset::BitSet|, which copy the set of bits into
  their own value and therefore will traverse the bits that were set at the
  time of their construction, usually by the |bitset::BitSet::begin| method).

  The most delicate operation here is the increment, which has to find the
  position of the next set bit, while avoiding falling off the bitmap if there
  is no such. Therefore we included the data for the end of the bitmap in the
  iterator. This also allows the |operator()| internal test for exhaustion.
  */
class BitMap::iterator // is really a const_iterator
{
  std::vector<unsigned long>::const_iterator d_chunk; // point to current chunk
  size_t d_bitAddress; // value returned if dereferenced
  size_t d_capacity; // beyond-the-end bit-address, normally constant

 public:

// associated types
  typedef std::forward_iterator_tag iterator_category;
  typedef size_t value_type;
  typedef ptrdiff_t difference_type;
  typedef const value_type* pointer;
  typedef const value_type& reference;

// constructors and destructors
  iterator() {}

  iterator(const iterator &j):d_chunk(j.d_chunk), d_bitAddress(j.d_bitAddress),
    d_capacity(j.d_capacity) {}

  iterator(const std::vector<unsigned long>::const_iterator& p,
	   size_t n, size_t c)
    :d_chunk(p), d_bitAddress(n), d_capacity(c) {}

// assignment
  iterator& operator= (const iterator& i);

// accessors
  bool operator== (const iterator& i) const
    { return d_bitAddress == i.d_bitAddress; } //  |d_chunk| is ignored!
  bool operator!= (const iterator& i) const
    { return d_bitAddress != i.d_bitAddress; } //  |d_chunk| is ignored!
  bool operator() () const { return d_bitAddress != d_capacity; }

  const value_type& operator* () const { return d_bitAddress; }

// manipulators
  iterator& operator++ ();

  iterator operator++ (int);

  void change_owner(const BitMap& b); // adapt |d_chunk| to point into |b|
}; // |class BitMap::iterator|

} // |namespace bitmap|

} // |namespace atlas|

#endif
