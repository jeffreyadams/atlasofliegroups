/*
  This is bitmap.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "bitmap.h"

#include "bits.h"

/*!****************************************************************************
\file
  \brief Contains the implementation of the BitMap class. 

  A BitMap should be seen as a container of _unsigned long_, not bits; the
  idea is that the unsigned longs it contains are the bit-addresses of
  the set bits. It obeys the semantics of a Forward Container (notion
  from the C++ standard library;

  The basic idea is that a bitmap is a vector of unsigned long, representing
  the "chunks" of the map. We wish to provide bit-address access to this
  map; for this purpose we use the reference trick from vector<bool>. Also
  we wish to define an iterator class, which traverses the _set_ bits of
  the bitmap; so that for instance, b.begin() would be a reference to the
  first set bit. Dereferencing the iterator yields its bit-address.

******************************************************************************/

namespace atlas {

namespace bitmap {

  // constants used to pick a bit-address apart
  // the first one serves as flags for the address within a word.
  // the second one is its logical complement
  // It is assumed that the number of digits in an unsigned long 
  // is a power of two.

  /*!
  Constant used to pick a bit-address apart: serves as flags for the
  address within a word. It is assumed that the number of bits in an
  unsigned long is a power of two.  
  */  
  unsigned long BitMap::posBits = constants::posBits;
  
  /*!
  Constant used to pick a bit-address apart: this is the logical
  complement of posBits. It is assumed that the number of bits in
  an unsigned long is a power of two.
  */
  unsigned long BitMap::baseBits = constants::baseBits;

  /*!
  Constant saying how much we have to shift the BitMap capacity n
  (that is, the power of two by which it much be divided) to get the
  base (that is, the number of unsigned longs in the needed in the
  BitMap).
  */
  unsigned long BitMap::baseShift = constants::baseShift;

}

/*****************************************************************************

        Chapter I -- The BitMap class

******************************************************************************/

namespace bitmap {

/******** constructors and destructors ***************************************/

BitMap::BitMap(unsigned long n)
  :d_map((n >> baseShift) +(bool)(n & posBits),0),d_capacity(n)

/*!
  Constructs a zero-initialized bitmap with a capacity of n bits. Notice
  that the size of the vector d_map is n >> baseShift + 1, except when
  longBits exactly divides n.
*/

{}

/******** assignment *********************************************************/

BitMap& BitMap::operator= (const BitMap& a)

{
  d_map = a.d_map;
  d_capacity = a.d_capacity;

  return *this;
}

/*!******* range bounds *******************************************************/

BitMap::iterator BitMap::begin() const

/*!
  Returns an iterator pointing to the first set bit in the bitmap.
*/

{
  unsigned long n = front();
  unsigned long base = n >> baseShift;

  return iterator(d_map.begin()+base,n,d_capacity);
}

BitMap::iterator BitMap::end() const

/*!
  Synopsis: returns the past-the-end iterator for the bitmap.
*/

{
  return iterator(d_map.end(),d_capacity,d_capacity);
}

BitMap::iterator BitMap::pos(unsigned long n) const

/*!
  Synopsis: returns an iterator with bit-address n.
*/

{  
  unsigned long base = n >> baseShift;
  return iterator(d_map.begin()+base,n,d_capacity);
}

/******** accessors **********************************************************/

unsigned long BitMap::back() const

/*!
  Synopsis: returns the address of the last member of the bitmap plus one,
  0 if there is no such bit (consistent with first, as a reverse iterator.)
*/

{
  using namespace bits;
  using namespace constants;

  if (empty())
    return 0;

  size_t b = d_map.size();

  for (; b;) {
    --b;
    if (d_map[b])
      break;
  }

  unsigned long n = b*longBits;
  n += lastBit(d_map[b]);

  return n;
}

bool BitMap::contains(const BitMap& b) const

/*!
  Tells whether the current bitmap contains b. It is assumed that both have
  the same capacity. It amounts to checking that b andnot this is empty.
*/

{
  for (unsigned long j = 0; j < d_map.size(); ++j)
    if (b.d_map[j] & ~d_map[j])
      return false;

  return true;
}

bool BitMap::empty() const

/*!
  Tells whether the bitmap is empty. Thanks to our convention of zeroing
  unused bits, it is enough to check whether all the components of d_map
  are zero.
*/

{
  for (unsigned long j = 0; j < d_map.size(); ++j)
    if (d_map[j])
      return false;

  return true;
}

unsigned long BitMap::front() const

/*!  
  Synopsis: returns the address of the first member (set bit) of
  the bitmap, past-the-end if there is no such.
*/

{
  using namespace bits;
  using namespace constants;

  if (empty())
    return d_capacity;

  size_t b = 0;

  for (; b < d_map.size(); ++b)
    if (d_map[b])
      break;

  unsigned long n = b*longBits;
  n += firstBit(d_map[b]);

  return n;
}

bool BitMap::full() const

/*!
  Tells whether the bitmap is full. This means that all blocks are full,
  except maybe for the last one where we have to look only at the significant
  part.
*/

{
  if (d_capacity == 0)
    return true;

  if (d_map[d_map.size()-1] != constants::leqMask[(d_capacity-1)&posBits])
    return false;

  if (d_map.size() > 1)
    for (unsigned long j = 0; j < d_map.size()-1; ++j)
      if (~d_map[j])
	return false;

  return true;
}

unsigned long BitMap::n_th(unsigned long n) const

/*!
  Synopsis: returns the position of the n-th set bit; returns d_capacity
  if there is no such.
*/

{
  using namespace bits;
  using namespace constants;

  typedef std::vector<unsigned long>::const_iterator I;

  unsigned long pos = 0;
  unsigned long f;

  I map_end = d_map.end();

  for (I i = d_map.begin(); i != map_end; ++i) {
    unsigned long chunkSize = bits::bitCount(*i);
    if (chunkSize > n) {
      f = *i;
      goto found;
    }
    else {
      n -= chunkSize;
      pos += longBits;
    }
  }

  // if we reach this point the bit is not found
  return d_capacity;

 found:
  // at this point the required bit is going to be in the current chunk

  for (size_t j = 0; j < n; ++j)
    f &= (f-1);

  pos += firstBit(f);

  return pos;
}

unsigned long BitMap::position(unsigned long n) const

/*!
  Synopsis: returns the number of set bits in position < n; if n is set,
  this is also the position of j within the set of set bits.
*/

{
  using namespace bits;
  using namespace constants;

  unsigned long p = 0;
  unsigned long b = n >> baseShift;

  for (size_t j = 0; j < b; ++j)
    p += bitCount(d_map[j]);

  unsigned long f = d_map[b] & lMask[n&posBits];
  p += bitCount(f);

  return p;
}

unsigned long BitMap::range(unsigned long n, unsigned long r) const

/*!
  Synopsis: returns r bits from position n.

  Precondition: r divides longBits, and n is aligned (i.e., n is a multiple
  of r).

  In this way, we are sure that things happen inside a single word in d_map.
*/

{
  using namespace constants;

  unsigned long m = n >> baseShift;

  return (d_map[m] >> (n & posBits)) & lMask[r];

}

unsigned long BitMap::size() const

/*!
  Returns the number of set bits in the bitmap (this is its size as a
  container of unsigned long.)
*/

{
  typedef std::vector<unsigned long>::const_iterator I;

  unsigned long c = 0;
  I map_end = d_map.end();

  for (I i = d_map.begin(); i != map_end; ++i)
    c += bits::bitCount(*i);

  return c;
}

/******** manipulators *******************************************************/

BitMap& BitMap::operator~ ()

/*!
  Synopsis: transforms the complement into its bitwise complement.

  NOTE: one has to be careful about the last chunk, resetting the unused
  bits to zero.
*/

{
  for (unsigned long j = 0; j < d_map.size(); ++j)
    d_map[j] = ~d_map[j];

  d_map[d_map.size()-1] &= constants::leqMask[(d_capacity-1)&posBits];

  return *this;
}

BitMap& BitMap::operator&= (const BitMap& b)

/*!
  Synopsis: intersects the current bitmap with b.
*/

{  
  for (unsigned long j = 0; j < d_map.size(); ++j)
    d_map[j] &= b.d_map[j];

  return *this;
}

BitMap& BitMap::operator|= (const BitMap& b)

/*!
  Synopsis: unites the current bitmap with b
*/

{  
  for (unsigned long j = 0; j < d_map.size(); ++j)
    d_map[j] |= b.d_map[j];

  return *this;
}

BitMap& BitMap::operator^= (const BitMap& b)

/*!
  Synopsis: xor's the current bitmap with b
*/

{  
  for (unsigned long j = 0; j < d_map.size(); ++j)
    d_map[j] ^= b.d_map[j];

  return *this;
}

BitMap& BitMap::andnot(const BitMap& b)

/*!
  Synopsis: takes the current bitmap into its set-difference with b, i.e. we 
  remove from our bitmap its intersection with b.
*/

{  
  for (unsigned long j = 0; j < d_map.size(); ++j)
    d_map[j] &= ~(b.d_map[j]);

  return *this;
}

void BitMap::fill()

/*!
  Synopsis: sets all the bits in the bitmap. 

  As usual we have to be careful to leave the unused bits at the end to zero.
*/

{
  if (d_capacity == 0) // do nothing
    return;

  d_map.assign(d_map.size(),~0ul);
  d_map[d_map.size()-1] &= constants::leqMask[(d_capacity-1)&posBits];

  return;
}

void BitMap::fill(unsigned long n)

/*!
  Synopsis: sets all the bits in position < n.

  Equivalent to fill if n >= d_capacity.
*/

{
  if (n >= d_capacity) { // fill
    fill();
    return;
  }

  size_t map_size = d_map.size();
  d_map.assign(n >> baseShift,~0ul);
  d_map.resize(map_size,0ul);
  d_map[n >> baseShift] = constants::lMask[n&posBits];

  return;
}

BitMap::iterator BitMap::insert(iterator, unsigned long n)

{
  d_map[n >> baseShift] |= constants::bitMask[n & posBits];
  return iterator(d_map.begin()+(n >> baseShift),n,d_capacity);
}
void BitMap::resize(unsigned long n)

/*!
  Synopsis: sets the capacity of the bitmap to n. 

  Does not modify the contents up to the previous size, at least if n is 
  larger. The new elements are initialized to zero.
*/

{
  d_map.resize((n >> baseShift) +(bool)(n & posBits),0);
  d_capacity = n;

  return;
}

void BitMap::setRange(unsigned long n, unsigned long r, unsigned long a)

/*!
  Synopsis: sets r bits from position n to the first r bits in a.

  Precondition: r divides longBits, and n is aligned (i.e., n is a multiple
  of r).

  In this way, we are sure that things happen inside a single word in d_map.
*/

{
  using namespace constants;

  unsigned long m = n >> baseShift;

  d_map[m] &= ~(lMask[r] << (n & posBits));     // set bits to zero
  d_map[m] |= (a & lMask[r]) << (n & posBits);  // insert the bits from a

  return;
}

void BitMap::swap(BitMap& other)

{
  d_map.swap(other.d_map);
  std::swap(d_capacity,other.d_capacity);

  return;
}

}

/*****************************************************************************

        Chapter II -- The BitMap::iterator class

  Because of the nature of a bitmap, only constant iterators make sense;
  one cannot "change the value" of an element at a given position because
  the value _is_ the position. The most delicate operation is the increment,
  which has to find the position of the next set bit, while avoiding falling
  off the bitmap if there is no such. Because of this we had to pass the
  data for the end of the bitmap into the iterator.

******************************************************************************/

namespace bitmap {

BitMap::iterator& BitMap::iterator::operator= (const BitMap::iterator& i)

{
  d_chunk = i.d_chunk;
  d_bitAddress = i.d_bitAddress;
  d_capacity = i.d_capacity;
  
  return *this;
}

BitMap::iterator& BitMap::iterator::operator++ ()

/*!
  The incrementation operator; it has to move the bitAddress to the next
  non-set bit, and move the chunk if necessary.
*/

{
  using namespace bits;
  using namespace constants;

  // put in f the remaining part of the chunk, shifted to the right edge;

  unsigned long f = *d_chunk >> (d_bitAddress & posBits);
  f >>= 1;
  
  if (f) { // we stay in the same chunk
    d_bitAddress += firstBit(f)+1;
    return *this;
  }

  d_bitAddress &= baseBits;  // remove bit-position part

  if ((d_capacity - d_bitAddress) <= longBits) { // we have reached the end
    d_bitAddress = d_capacity;
    return *this;
  }

  // go to beginning of next chunk

  d_bitAddress += longBits;  // first bit in next chunk
  ++d_chunk;
  
  for(; (d_capacity - d_bitAddress) > longBits; ++d_chunk) {
    if (*d_chunk) { // first bit is found
      d_bitAddress += firstBit(*d_chunk);
      return *this;
    }
    else
      d_bitAddress += longBits;
  }

  // when we reach this point, we have reached the last chunk

  if (*d_chunk) 
    d_bitAddress += firstBit(*d_chunk);
  else // past-the-end
    d_bitAddress = d_capacity;

  return *this;
}

BitMap::iterator BitMap::iterator::operator++ (int)

/*!
  Post-increment operator; it should return the value as it was _before_ the
  incrementation.
*/

{
  BitMap::iterator tmp = *this;
  ++(*this);
  return tmp;
}

}

}
