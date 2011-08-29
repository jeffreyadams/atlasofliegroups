/*!\file
  \brief Definitions and declarations for the BitMap class.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
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
class BitMap
{
  size_t d_capacity;
  std::vector<unsigned long> d_map;
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

// constructors and destructors
 BitMap() : d_capacity(0), d_map() {} // create a bitmap without capacity

  /*! \brief
    Constructs a zero-initialized bitmap with a capacity of n bits.

    Note that the size of the vector |d_map| exceeds |n >> baseShift| by
    one, unless |longBits| exactly divides |n|.
  */
  explicit BitMap(unsigned long n)
    : d_capacity(n), d_map((d_capacity+posBits)>>baseShift,0)
  {}

  //! \brief Copy constructor
  BitMap(const BitMap& b) : d_capacity(b.d_capacity), d_map(b.d_map) {}

  // convert range defined by iterators into a BitMap
  template <typename I>
    BitMap(unsigned long n, const I& first, const I& last);

  // an easier version for a full vector
  template <typename U> // unsigned integral type
  BitMap(unsigned long n,const std::vector<U>& v)
    : d_capacity(n), d_map((d_capacity+posBits)>>baseShift)
  {
    for (size_t i=0; i<v.size(); ++i)
      insert(v[i]);
  }

  //! Set of offsets into [first,last[ of values (also) found in [fsub,lsub[
  template <typename I, typename J>
    BitMap(const I& first, const I& last, const J& fsub, const J& lsub);


// assignment
 BitMap& operator= (const BitMap&);

// accessors
 /*!
   Number of bits in use in the bitmap. This is the capacity
   of the BitMap as a standard library container, not d_map.size(), which is
   approximately longBits times smaller.
 */
 unsigned long capacity() const { return d_capacity; }
 size_type size() const; // the number of bits that are set in the bitmap

 bool operator< (const BitMap& b) const { return d_map < b.d_map; }
 bool operator== (const BitMap& b) const { return d_map == b.d_map; }
 bool operator!=(const BitMap& b) const { return d_map != b.d_map; }

 bool empty() const; // whether |size()==0|
 bool full() const; // whether |size()==capacity()|

 /*!
   Tests whether bit n in the bitmap is set; that is, whether element n
   is a member of the set.
 */
 bool isMember(unsigned long n) const
 {
   if (n >= d_capacity)
     return false;
   return (d_map[n >> baseShift] & constants::bitMask[n&posBits])!=0;
 }

 // Whether all elements of |b| satisfy |isMember|
 bool contains(const BitMap& b) const;

 // Whether none of the elements of |b| satisfy |isMember|
 bool disjoint(const BitMap& b) const;

 //! \brief Value at index |n| if viewed as list of |unsigned long| values
 unsigned long n_th(unsigned long n) const;

 //! \brief Number of values |<n| present (set) in the bitmap
 unsigned long position(unsigned long n) const;

 unsigned long front() const;

 // decrement |n| (at least once) until it points to a member
 // return value indicates whether this succeeded
 bool back_up(unsigned long& n) const;

 // get a range of bits as unsigned long value; see bitmap.ccp for details
 unsigned long range(unsigned long first, unsigned long number) const;

 BitMap operator& (const BitMap& other) const
 { BitMap result(*this); result&=other; return result; }

 BitMap operator| (const BitMap& other) const
 { BitMap result(*this); result|=other; return result; }

 BitMap operator^ (const BitMap& other) const
 { BitMap result(*this); result^=other; return result; }


// manipulators

 // WARNING: the following looks like an accessor, but complements |*this|
 BitMap& operator~ ();

 // these operators allow argument to have less capacity than |*this|
 bool operator&= (const BitMap&);

 BitMap& operator|= (const BitMap&);

 BitMap& operator^= (const BitMap&);

 bool andnot(const BitMap& b); // remove bits of |b|, return whether any left

 /*!
   Set the bit at position n (that is, inserts the value |n| into the set);
   this makes |isMember(n)| hold.
 */
 void insert(unsigned long n)
 { d_map[n >> baseShift] |= constants::bitMask[n & posBits]; }

 /*!
   Clear the bit at position n (that is, removes an element of the set);
   this makes |isMember(n)| false.
 */
 void remove(unsigned long n)
 { d_map[n >> baseShift] &= ~constants::bitMask[n & posBits]; }

 void set_to(unsigned long n,bool b)
 { if (b) insert(n); else remove(n); }

 void set_mod2(unsigned long n, unsigned long v) { set_to(n,(v&1)!=0); }

 void flip(unsigned long n)
 { d_map[n >> baseShift] ^= constants::bitMask[n & posBits]; }


 void fill(); // set all bits
 void reset() { d_map.assign(d_map.size(),0ul);  } // clear all bits
 void fill(size_t start,size_t stop); // set consecutive range of bits
 void clear(size_t start,size_t stop); // clear consecutive range of bits
 // insert a range of values (which need not be listed increasingly)
 template<typename I>
   void insert(I, I);

 // this was called |resize|, but sets |capacity()|, whence the new name
 void set_capacity(unsigned long n); // any new bits will start out cleared

 // set an interval of bits from those (least significant ones) of source
 void setRange(unsigned long start, unsigned long amount, unsigned long source);

 void swap(BitMap&);

}; // |class BitMap|



  /*!
  \brief Traverses the set bits of a BitMap.

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
  unsigned long d_bitAddress; // value returned if dereferenced
  unsigned long d_capacity; // beyond-the-end bit-address, normally constant

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
