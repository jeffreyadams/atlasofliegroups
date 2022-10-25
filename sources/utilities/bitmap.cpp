/*
  This is bitmap.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006,2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "bitmap.h"

#include <algorithm> // for |lower_bound|, |copy|, |copy_backward|
#include "bits.h"

#include <cassert>

/*****************************************************************************

  The implementation of the BitMap class.

  A BitMap should be seen as a container of some unsigned integer type, with
  values bounded by the capacity of the bitmap, and which will be produced as
  |unsigned long| values (though these can then be narrowed to a smaller type
  that can fit all possible values). The fact that bits are used to signal the
  presence or absence of numbers is an inplementation details that this class
  somewhat hides (though the class name is an obvious giveaway). It obeys the
  semantics of a Forward Container (from the C++ standard library).

  A bitmap is a implemented as a vector of |unsigned long| (no relation with
  the same type used to represent values the container holds), each
  representing a "chunk" of bits in the map. We do not wish to provide
  bit-address access to this map. We do wish to define an iterator class,
  which traverses the _set_ bits of the bitmap; so that dereferencing for
  instance |b.begin()| returns the index of the first set bit in the bitmap.

******************************************************************************/

namespace atlas {

namespace bitmap {

  // constants used to pick a bit-address apart
  // the first one serves as mask for the bit-address within a word.
  // the second one is its logical complement; mask for the word address
  // It is assumed that the number of digits in an unsigned long
  // is a power of two.

  /*
  Constant used to pick a bit-address apart: serves as a bitmask for
  obtaining the bit-address within a word fom a BitMap index. It is assumed
  that the number of bits in an unsigned long is a power of two.
  */
  unsigned long BitMap::posBits = constants::posBits;

  /*
  Constant used to pick a bit-address apart: this is the logical
  complement of posBits, and masks the word-address within a BitMap index
  (which still must be shifted right by baseShift to be interpreted correctly,
  whence this constant is actually little used).

  It is assumed that the number of bits in an unsigned long is a power of two.
  */
  unsigned long BitMap::baseBits = constants::baseBits;

  /*
  Constant saying how much we have to shift the BitMap index n of a bit (that
  is, the power of two by which it much be divided) to get the index of the
  d_map element that contains this bit (it is the number of set bits in
  posBits, typically 5 or 6).
  */
  unsigned long BitMap::baseShift = constants::baseShift;

}

/*****************************************************************************

        Chapter I -- The BitMap class

******************************************************************************/

namespace bitmap {

/******** constructors and destructors ***************************************/

template <typename I> // iterator type
  BitMap::BitMap(unsigned long n, const I& first, const I& last)
  : d_capacity(n)
  , d_map((d_capacity+posBits)>>baseShift,0)
{
  for (I it = first; it!=last; ++it)
    insert(*it);
}

/*
  In this constructor template we assume that I and J are iterator types with
  the same value_type. The idea is that [first,last[ is an ordered range,
  for which we can call lower_bound. Then we construct the bitmap which
  flags the elements from [fsub,lsub[ (not necessarily assumed ordered or
  in range; I should be random-access, but J can basically be any input
  iterator.) It is assumed of course that the elements from [fsub,lsub[
  will be found in [first,last[.
*/
template <typename I, typename J>
  BitMap::BitMap(const I& first, const I& last, const J& fsub, const J& lsub)
  : d_capacity(last-first)
  , d_map((d_capacity+posBits)>>baseShift,0)
{
  for (J j = fsub; j != lsub; ++j)
    insert(lower_bound(first,last,*j)-first);
}



/******** assignment *********************************************************/

BitMap& BitMap::operator= (const BitMap& a)

{
  d_map = a.d_map;
  d_capacity = a.d_capacity;

  return *this;
}

/******** range bounds *******************************************************/

/*
  Return an iterator pointing to the first set bit in the bitmap.

  If the bitset is empty, this is will be equal to |end()|.
*/
BitMap::iterator BitMap::begin() const
{
  unsigned long n=front();
  return iterator(d_map.begin()+(n>>baseShift),n,d_capacity);
}

/*
  Return the past-the-end iterator for the bitmap.

  This is only needed to allow using these iterators in genereric algorithms
  which typically do |for (iterator it=x.begin(), it!=x.end(); ++it)|. In code
  that knows which kind of iterator this is, using |it()| as second clause is
  to be preferred.

  Note that only the middle argument |d_capacity| is of any importance, since
  the only thing one can meaningfully do with end() is test for (in)equality.
  The operator |++| below does not in fact advance to |d_chunk==d_map.end()|!
*/
BitMap::iterator BitMap::end() const
{
  return iterator(d_map.end(),d_capacity,d_capacity);
}

/******** accessors **********************************************************/


/*
  Decrement |n| until it points to a member of the bitset,
  or if none is found returns |false| (in which case |n| is unchanged)
*/
bool BitMap::back_up(unsigned long& n) const
{
  unsigned int i=n>>baseShift;
  auto m=(n&posBits)==0 ? 0 : d_map[i]&constants::lMask[n&posBits];

  while(m==0 and i>0)
    m=d_map[--i];

  if (m==0) return false;
  n = (i<<baseShift)+bits::lastBit(m)-1; // since |lastBit| overshoots by 1

  return true;
}

unsigned long last (const BitMap& set)
{
  unsigned long result = set.capacity();
#ifdef NDEBUG
  set.back_up(result);
#else
  bool success = set.back_up(result);
  assert(success);
#endif
  return result;
}

/*
  Whether the current bitmap contains |b|. It is assumed that
  |b.capacity()<=capacity()|.

  This would amount to |not b.andnot(*this)| if |b| were by-value rather
  than reference, and if capacities were equal.
*/
bool BitMap::contains(const BitMap& b) const
{
  assert(b.capacity()<=capacity());

  for (unsigned long j = 0; j < b.d_map.size(); ++j)
    if ((b.d_map[j] & ~d_map[j])!=0)
      return false;

  return true;
}

// Whether the current bitmap is disjoint from |b|.
bool BitMap::disjoint(const BitMap& b) const
{
  for (unsigned long j = 0; j < d_map.size() and j < b.d_map.size(); ++j)
    if ((d_map[j] & b.d_map[j])!=0)
      return false;

  return true;
}


/*
  Whether the bitmap is empty. Thanks to our convention of zeroing
  unused bits, it is enough to check whether all the components of |d_map|
  are zero.
*/
bool BitMap::empty() const
{
  for (unsigned long j = 0; j < d_map.size(); ++j)
    if (d_map[j]!=0)
      return false;

  return true;
}


/*
  Return the address of the first member (set bit) of the bitmap,
  or past-the-end indicator |d_capacity| if there is no such.
*/
unsigned long BitMap::front() const
{
  size_t b = 0;
  for (; b < d_map.size(); ++b)
    if (d_map[b]!=0) break;

  if (b==d_map.size()) // then all bits were clear
    return d_capacity;
  else
    return (b<<baseShift) + bits::firstBit(d_map[b]);
}

/*
  Tell whether the bitmap is full. This means that all blocks are full,
  except maybe for the last one where we have to look only at the significant
  part.
*/
bool BitMap::full() const
{
  /* Test final block if partial. There should be exactly |n| bits set (in
     other words the value should be |lMask[n]|), where |n>0| is congruent to
     |d_capacity%longBits|. Fokko had used |leqMask[(d_capacity-1)&posBits]|,
     but needed an exception for |d_capacity==0|.
   */
  if ((d_capacity&posBits)!=0
      and d_map.back()!=constants::lMask[d_capacity&posBits])
    return false;

  for (unsigned long j = (d_capacity>>baseShift); j-->0; )
    if (~d_map[j]!=0)
      return false;

  return true;
}


/*
  Return the index of set bit number |i| in the bitset; in other words,
  viewing a bitset |b| as a container of |unsigned long|, |b.n_th(i)| is the
  value of the element |i| of |b|, and the syntax |b[i]| would have been
  logical (as usual, the first element is number 0). This returns |d_capacity|
  if there is no such element, in other words if at most |i| bits are set in
  the bitmap. Whenever |0<=i<=size()| one has |b.position(b.n_th(i))==i|.
*/
unsigned long BitMap::n_th(unsigned long i) const
{
  unsigned long pos = 0;
  unsigned long f;

  for (std::vector<unsigned long>::const_iterator
	 iter = d_map.begin(); iter != d_map.end(); ++iter)
  {
    unsigned long chunkSize = bits::bitCount(*iter);
    if (chunkSize>i)
    {
      f = *iter;
      goto found;
    }
    else
    {
      i -= chunkSize;
      pos += constants::longBits;
    }
  }

  // if we reach this point the bit is not found
  return d_capacity;

 found:
  // at this point the required bit is going to be in the current chunk

  for (size_t j=0; j<i; ++j) // clear |i| initial bits in |f|
    f &= (f-1);

  return pos + bits::firstBit(f); // some bit i
}

/*
  Return the number of set bits in positions |< n|; viewing a bitset |b| as
  a container of |unsigned long|, this is the number of values |< n| that |b|
  contains. If |n| itself is a member of |b|, then |n==b.n_th(b.position(n))|.
*/
unsigned long BitMap::position(unsigned long n) const
{
  unsigned long p = 0;
  unsigned long b = n >> baseShift;

  for (size_t j = 0; j < b; ++j)
    p += bits::bitCount(d_map[j]);

  if ((n&posBits)!=0)
    p += bits::bitCount(d_map[b] & constants::lMask[n&posBits]);

  return p;
}


/*
  Return |r| bits from position |n|.

  Precondition: |r| divides |longBits|, and |n| is a multiple of |r|.

  Thus the bits extracted are found in single element of d_map, and such
  elements define an integral number of disjoint ranges

  It is required that |n<capacity()|, but not that |n+r<=capacity()|; if the
  latter fails, the return value is padded out with (leading) zero bits.
*/
unsigned long BitMap::range(unsigned long n, unsigned long r) const
{
  unsigned long m = n >> constants::baseShift; // index where data is stored

  return (d_map[m] >> (n & constants::posBits)) & constants::lMask[r];

}

/*
  Return the number of set bits in the bitmap (this is its size as a
  container of unsigned long.)

  NOTE: correctness depends on unused bits in the final word being cleared.
*/
unsigned long BitMap::size() const
{
  unsigned long c = 0;

  for (std::vector<unsigned long>::const_iterator i = d_map.begin();
       i != d_map.end(); ++i)
    c += bits::bitCount(*i);

  return c;
}

/******** manipulators *******************************************************/

/*
  Replace a bitmap by its complement

  We are careful about the last chunk, leaving the unused bits equal to zero.
*/
BitMap& BitMap::take_complement ()
{
  for (unsigned long j = 0; j < d_map.size(); ++j)
    d_map[j] = ~d_map[j];

  if ((d_capacity&posBits)!=0) // N.B. also makes |capacity()==0| safe (MvL)
    d_map.back() &= constants::lMask[d_capacity&posBits];

  return *this;
}

// Intersect the current bitmap with |b|, return whether result is non-empty.
bool BitMap::operator&= (const BitMap& b)
{
  assert(b.capacity()<=capacity());
  bool any=false;
  for (unsigned long j = 0; j < b.d_map.size(); ++j)
    if ((d_map[j] &= b.d_map[j])!=0) any=true;

  // don't forget to clear out everything beyond the end of |b|!
  for (unsigned long j = b.d_map.size(); j<d_map.size(); ++j)
    d_map[j] = 0;

  return any;
}

// Unite |b| into the current bitmap.
BitMap& BitMap::operator|= (const BitMap& b)
{
  assert(b.capacity()<=capacity());
  for (unsigned long j = 0; j < b.d_map.size(); ++j)
    d_map[j] |= b.d_map[j];

  return *this;
}

// Perform exclusive or (XOR) of |b| into the current bitmap.
BitMap& BitMap::operator^= (const BitMap& b)
{
  assert(b.capacity()<=capacity());
  for (unsigned long j = 0; j < b.d_map.size(); ++j)
    d_map[j] ^= b.d_map[j];

  return *this;
}

/*
  Take the current bitmap into its set-difference with |b|, i.e.,
  remove from our bitmap any elements appearing in |b|; the latter should not
  exceed the size of the current bitmap (but may be smaller).
  Return whether any bits remain in the result.
*/
BitMap& BitMap::andnot(const BitMap& b)
{
  assert(b.capacity()<=capacity());
  for (unsigned long j = 0; j < b.d_map.size(); ++j)
    d_map[j] &= ~b.d_map[j];

  return *this;
}

BitMap& BitMap::operator<<= (unsigned long delta) // increase values by |delta|
{
  unsigned long delta_rem = delta & posBits;
  delta >>= baseShift; // we must move |delta| words, and then |delta_rem| bits
  if (delta>0) // shifting by |0| is useless, and undefined behavior too
  {
    std::copy_backward(d_map.begin(),d_map.end()-delta,d_map.end());
    std::fill(d_map.begin(),d_map.begin()+delta,0ul);
  }
  if (delta_rem>0 and not d_map.empty())
  {
    std::vector<unsigned long>::iterator it;
    for (it=d_map.end(); --it!=d_map.begin(); ) // reverse loop, omit initial
    {
      *it <<= delta_rem; // shift bits up
      *it |= *(it-1) >> (constants::longBits-delta_rem); // bits shifted in
    }
    *it <<= delta_rem; // shift bits up
  }
  return *this;
}

BitMap& BitMap::operator>>= (unsigned long delta) // decrease values by |delta|
{
  unsigned long delta_rem = delta & posBits;
  delta >>= baseShift; // we must move |delta| words, and then |delta_rem| bits
  if (delta>0) // shifting by |0| is useless, and undefined behavior too
  {
    std::copy(d_map.begin()+delta,d_map.end(),d_map.begin());
    std::fill(d_map.end()-delta,d_map.end(),0ul);
  }
  if (delta_rem>0 and not d_map.empty())
  {
    std::vector<unsigned long>::iterator it;
    for (it=d_map.begin(); it+1!=d_map.end(); ++it)
    {
      *it >>= delta_rem; // shift bits doan
      *it |= *(it+1) << (constants::longBits-delta_rem);
    }
    *it <<= delta_rem; // shift last bits down
  }
  return *this;
}

/*
  Set all the bits in the bitmap.

  As usual we have to be careful to leave the unused bits at the end to zero.
*/
void BitMap::fill()
{
  d_map.assign(d_map.size(),~0ul); // (last word will be overwritten)

  if ((d_capacity&posBits)!=0) // N.B. also makes |capacity()=0| safe (MvL)
    d_map.back() = constants::lMask[d_capacity&posBits];
}

// Set all the bits in positions |i| with |start<=i<stop|.
void BitMap::fill(size_t start, size_t stop)
{
  if (start>=stop) return;
  size_t begin = start >> baseShift;
  size_t end   = stop >> baseShift;
  if (begin==end) // then only one word is affected
    d_map[begin] |= constants::lMask[stop-start] << (start&posBits);
  else
  {
    d_map[begin] |= ~constants::lMask[start&posBits];
    for (size_t i=begin+1; i<end; ++i)
      d_map[i] = ~0ul;
    if ((stop&posBits)!=0) // protect against out-of-bounds setting of 0 bits
      d_map[end] |= constants::lMask[stop&posBits];
  }
}

// Set all the bits in positions |i| with |start<=i<stop|.
void BitMap::clear(size_t start, size_t stop)
{
  if (start>=stop) return;
  size_t begin = start >> baseShift;
  size_t end   = stop >> baseShift;
  if (begin==end) // then only one word is affected
    d_map[begin] &= ~(constants::lMask[stop-start] << (start&posBits));
  else
  {
    d_map[begin] &= constants::lMask[start&posBits];
    for (size_t i=begin+1; i<end; ++i)
      d_map[i] = 0ul;
    if ((stop&posBits)!=0) // protect against out-of-bounds setting of 0 bits
      d_map[end] &= ~constants::lMask[stop&posBits];
  }
}

/*
  Here we assume that I is an iterator whose value_type is unsigned long,
  and we do the sequence of insertions from the range [first,last[.
*/
template<typename I> void BitMap::insert(I first, I last)
{
  while (first!=last)
  {
    insert(*first);
    ++first;
  }
}


/*
  Set the capacity of the bitmap to |n|, shrinking |d_map| if possible

  Does not modify the contents up to the previous size, at least if n is
  larger. The new elements are initialized to zero.
*/
void BitMap::set_capacity(unsigned long n)
{
  if (n<d_capacity and (n&posBits)!=0)
    d_map[n>>baseShift] &= constants::lMask[n&posBits]; // clear partial word
  d_map.resize((n+posBits) >> baseShift,0);
  if (d_map.size()<d_map.capacity()) // |d_map| has become too large
    std::vector<unsigned long>(d_map).swap(d_map); // shrink-wrap vector |d_map|
  d_capacity = n;
}

// Add one more place to bitmap (amortised efficiently), set it to value |b|
void BitMap::extend_capacity(bool b)
{
  if (d_capacity == (d_map.size()<<baseShift))
    d_map.push_back(0); // allocate space; |std::vector| handles efficiency
  set_to(d_capacity++,b);
}

/*
  Set |r| bits from position |n| to the first |r| bits in |a|.

  Precondition: |n| and |n+r-1| have same integer quotient by |longBits| (or
  |r==0| in which case nothing happens), so in particular |r<=longBits|.
  This condition is always satisfied if |r| divides |longBits| and |n| is a
  multiple of |r|. It ensures that only a single word in |d_map| is affected.
*/
void BitMap::setRange(unsigned long n, unsigned long r, unsigned long a)
{
  unsigned long m = n >> baseShift;
  unsigned shift = n & posBits;

  d_map[m] &= ~(constants::lMask[r] << shift);     // set bits to zero
  d_map[m] |= (a & constants::lMask[r]) << shift;  // insert the bits from a
}

void BitMap::swap(BitMap& other)

{
  d_map.swap(other.d_map);
  std::swap(d_capacity,other.d_capacity);
}

}

/*****************************************************************************

        Chapter II -- The BitMap::iterator class

  Because of the nature of a bitmap, only constant iterators make sense (just
  like for iterators into |std::set|); one cannot "change the value" of an
  element at a given position because the position is determined by the value
  (values are always increasing; in fact a small value-change could be
  realised by swapping a set bit with neighboring unset bits, but that is of
  course not what a non-constant iterator should allow doing). On the other
  hand we may modify our bitmap |M| using this iterator |it|, notably to clear
  the set bit |it| points at by |M.remove(*it)|. Unlike the situation for
  |std::set|, such a removal does not invalidate the iterator |it| itself, so
  it is not an invariant of the |BitMap::iterator| class that it always points
  at a set bit.

  The most delicate operation is |operator++|, which has to move to the the
  next set bit, or stop at the end of the bitmap if there is no such.
  Therefore we included the data for the end of the bitmap in the iterator.

******************************************************************************/

namespace bitmap {

BitMap::iterator& BitMap::iterator::operator= (const BitMap::iterator& i)

{
  d_chunk = i.d_chunk;
  d_bitAddress = i.d_bitAddress;
  d_capacity = i.d_capacity;

  return *this;
}

/*
  The incrementation operator; it has to move the bitAddress to the next
  set bit, and move the chunk if necessary.

  This code below assumes that in case of an incomplete last chunk, there are
  no bits set in that chunk beyond the end of the bitmap; if there were,
  |firstBit(f)| below (both instances) could make the iterator advance to such
  a bit when it should have halted at |d_capacity|.
*/

BitMap::iterator& BitMap::iterator::operator++ ()
{
  // re-fetch current chunk, and shift away any bits we already passed over
  unsigned long f = *d_chunk >> (d_bitAddress & constants::posBits);

  //  here |(f&1)!=0|, unless we just did |remove(*it)| for our iterator |it|

  // separating the following shift from the previous one avoids an undefined
  // shift over longBits=posBits+1 when iterator steps into the next chunk
  f >>= 1; // discard the bit we were placed at; we move on!

  if (f!=0) { // if there is still some bit set in this chunk, jump to it
    d_bitAddress += bits::firstBit(f)+1; // +1 to account for the f>>=1
    return *this;
  }

  // if not, we're done with this chunk; we must advance d_chunk at least once
  d_bitAddress &= baseBits;  // prepare to advance by multiples of |longBits|

  do
    if ((d_bitAddress+=constants::longBits)<d_capacity)
      f=*++d_chunk;
    else
    {
      d_bitAddress = d_capacity; // so now we we are out-of-bounds
      return *this;
    }
  while(f==0);

  // now the bit to advance to is in the current chunk, and not out of bounds
  d_bitAddress += bits::firstBit(f);
  return *this;
}


/*
  Post-increment operator; it should return the value as it was _before_ the
  incrementation. This operator can mostly by avoided, as |M.remove(*it++)|
  can safely be replaced by |M.remove(*it),++it|
*/
BitMap::iterator BitMap::iterator::operator++ (int)
{
  BitMap::iterator tmp = *this;
  ++(*this);
  return tmp;
}

void BitMap::iterator::change_owner(const BitMap& b)
{
  assert (b.capacity()==d_capacity);
  d_chunk = b.d_map.begin()+(d_bitAddress>>baseShift);
}

// Instantiations

template void BitMap::insert
 (std::vector<unsigned short>::iterator,
  std::vector<unsigned short>::iterator); // root sets from RootNbrList

using ULVI = std::vector<unsigned long>::iterator;
using VCI = std::vector<int>::const_iterator;

template BitMap::BitMap(unsigned long n, const ULVI& first, const ULVI& last);
template BitMap::BitMap(unsigned long n, const VCI& first, const VCI& last);
// template BitMap::BitMap(const VI& f,const VI& l, const VI& sf,const VI& sl);

} // |namespace bitmap|

} // |namespace atlas|
