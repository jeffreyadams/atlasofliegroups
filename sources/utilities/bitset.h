/*
  This is bitset.h.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Class definitions and function declarations for the BitSet class.

#ifndef BITSET_H  /* guard against multiple inclusions */
#define BITSET_H

#include "bitset_fwd.h"

#include <cstdint> // for |uint_least32_t|
#include <iterator>
#include <cassert>

#include "bits.h"
#include "constants.h"

/******** type definitions ***************************************************/

namespace atlas {

namespace bitset {

/*****************************************************************************

  The |BitSet| class has essentially the functionality of the |stl::bitset
  class| from the STL, but I have redone it because I needed some things that
  bitset doesn't seem to provide. This is a very important data type for the
  program, and it is crucial that it should be efficiently implemented. [Fokko]

  An important difference with |stl::bitset| is that whereas the latter models
  sets of exactly |n| bits, where |n| is the template parameter and therefore
  fixed at compile time, |BitSet| will often be used to model sets of at most
  |n| bits, as in Atlas the exact size of almost nothing is fixed at compile
  time. We leave the available but unused bits simply unset, which allows some
  methods like |any| or |none| to be implemented without knowing the actual
  number of bits in use, while for a few methods like |complement| it means
  they require a second |size| (and in this particular case it cannot be
  called |operator~|).

  Only a few fixed template instances will be defined, with no template code
  handling general values of |n|. The intention is to use |BitSet| for small
  bitmaps with size bounded by a constant (in fact 1 or 2) times the rank, so
  that the corresponding multiple of |constants::RankMax| can be used as
  template argument. This allows for optimal speed, at the cost of some code
  duplication.

  Just as for the implementation of bitset in the STL that I took as a model,
  the template parameter |n| is passed to a base class template |BitSetBase|
  after conversion to number of |unsigned long|s needed rather than bits.
  Only a few instantiations of |BitSetBase| are provided, currently only the
  single-word version is complete and the 2-word version is partially defined.

  The main difference with the STL is the |to_ulong()| method; the one here
  truncates the bitset when the number of bits is too long, instead of
  throwing an exception (however most likely this can never occur).

  No bit-references are implemented, and |operator[]| is defined as accessor,
  equivalent to |test()|. In other respects method naming conventions are not
  yet coherent with those of |BitMap|, a desirable future change.

  The output operator is in io/ioutils.h.

******************************************************************************/

typedef uint_least32_t chunk; // unit of grouping of bits
typedef uint_least64_t wide_uint; // exchange format for all |BitSet| instances

/*
  First a template class |BitSetBase| is defined, which provides the basic
  implementation (a small fixed number of |chunk|s), but not yet the full
  external interface. Because instances of the template class |BitSet| will be
  publicly derived from instances of |BitSetBase| (with a different argument),
  many methods will be directly used by the |BitSet| template instances.
*/

template <unsigned int n> class BitSetBase; // instantiate only a few values |n|

/*
  Base class for a non-empty BitSet that fits in 32 bits.

  If |RANK_MAX<=16|, this template should be the only one instantiated.

  The class |BitSet<n>| for |1<=n<=32| will be derived from this instance.
*/
template<> class BitSetBase<1u>
{
  chunk d_bits; // one chunk of (at least) 32 bits

 protected:

// associated types
  class iterator; // (input) iterator type is an embedded class

// constructors
  BitSetBase<1>() : d_bits(0u) {}
  explicit BitSetBase<1>(wide_uint b): d_bits(static_cast<chunk>(b))  {}

// accessors

 public:
  bool operator==(const BitSetBase<1>& b) const { return d_bits == b.d_bits; }
  bool operator!=(const BitSetBase<1>& b) const { return d_bits != b.d_bits; }
  bool operator< (const BitSetBase<1>& b) const { return d_bits <  b.d_bits; }

  bool any()  const { return d_bits!=0; }
  bool none() const { return d_bits==0; }

  bool any(const BitSetBase<1>& b) const
  { return  (d_bits & b.d_bits)!=0; }
  bool contains (const BitSetBase<1>& b) const
  { return (~d_bits & b.d_bits)==0; }

  unsigned int count() const { return bits::bitCount(d_bits); }
  unsigned int firstBit() const // index of least set bit; |longBits| if none
  { return bits::firstBit(d_bits); }
  unsigned int lastBit() const // index of highest set bit PLUS ONE, 0 if none
  { return bits::lastBit(d_bits); }

  bool test(unsigned int j) const
  { assert(j<32); return (d_bits & constants::bitMask[j])!=0; }

  unsigned int position(unsigned int j) const // rank among set bits of bit |j|
  { assert(j<32); return bits::bitCount(d_bits & constants::lMask[j]); }

  bool scalarProduct(const BitSetBase<1>& b) const
  { return bits::bitCount(d_bits&b.d_bits)%2 != 0; }

  wide_uint to_ulong() const { return d_bits; }

  iterator begin() const;

// manipulators

 protected: // these cannot be called directly: no return value is available
  void operator^= (const BitSetBase<1>& b) { d_bits ^=  b.d_bits; }
  void operator|= (const BitSetBase<1>& b) { d_bits |=  b.d_bits; }
  void operator&= (const BitSetBase<1>& b) { d_bits &=  b.d_bits; }
  void andnot     (const BitSetBase<1>& b) { d_bits &= ~b.d_bits; }

  // next two methods make sure that shift by |32| yields 0.
  void operator<<= (unsigned int c)
  {
    if (c < 32) // bit shifts by more than this are undefined!
      d_bits <<= c;
    else
      d_bits = 0u; // simulate shifting out of all bits
  }
  void operator>>= (unsigned int c)
  {
    if (c < 32) // bit shifts by more than this are undefined!
      d_bits >>= c;
    else
      d_bits = 0u; // simulate shifting out of all bits
  }

  void flip(unsigned int j) { d_bits ^= constants::bitMask[j]; }
  void reset() { d_bits = 0u; }
  void reset(unsigned int j)  { d_bits &= ~constants::bitMask[j]; }
  void set(unsigned int j) { d_bits |= constants::bitMask[j]; }
  void fill(unsigned int limit) { d_bits |= constants::lMask[limit]; }

  void complement(unsigned int limit) { d_bits ^= constants::lMask[limit]; }
  void truncate(unsigned int limit)
    { if (limit<32) d_bits &= constants::lMask[limit]; }

  void slice(const BitSetBase<1>& c); // extract bits set in |c|, compacting
  void unslice(const BitSetBase<1>& c); // expand bits to positions set in |c|
  void swap(BitSetBase<1>& source) { std::swap(d_bits,source.d_bits); }

 }; // |class BitSetBase<1>|

/*
  Base for a non-empty BitSet that fits in two chunks but not one.

  The class BitSet<n>, for |33<=n<=64| will be a derived class of BitSetBase<2>.
*/
template<> class BitSetBase<2>
{
  // use two |chunk|s rather than a |wide_uint| for less strict alignment
  chunk d_bits0,d_bits1;

 protected:

// associated types
  class iterator; // (input) iterator type is an embedded class

// constructors and destructors
  BitSetBase<2>() : d_bits0(0u), d_bits1(0u) {}

/*
  The class |BitSet| will assume that any |BitSetBase| instance has a
  constructor with argument a |wide_uint|, which is guaranteed to contain at
  least 64 bits. We shift them into place into the two |chunk|s here.
*/
  explicit BitSetBase<2>(wide_uint b)
    : d_bits0(b&0xFFFFFFFF), d_bits1((b&0xFFFFFFFF00000000ul)>>32) {}

  BitSetBase<2>(const BitSetBase<2>& b) = default;

// accessors

 public:
  bool operator==(const BitSetBase<2>& b) const
    { return d_bits0==b.d_bits0 and d_bits1==b.d_bits1; }
  bool operator!=(const BitSetBase<2>& b) const
    { return d_bits0!=b.d_bits0 or d_bits1!=b.d_bits1; }
  bool operator< (const BitSetBase<2>& b) const
    { return d_bits1==b.d_bits1 ? d_bits0<b.d_bits0 : d_bits1<b.d_bits1; }

  bool any() const { return d_bits0!=0 or d_bits1!=0; }
  bool none() const { return d_bits0==0 and d_bits1==0; }

  bool any(const BitSetBase<2>& b) const
    { return (d_bits0 & b.d_bits0)!=0 or (d_bits1 & b.d_bits1)!=0; }
  bool contains (const BitSetBase<2>& b) const
    { return (~d_bits0 & b.d_bits0)==0 and (~d_bits1 & b.d_bits1)==0; }

  unsigned int count() const
    { return bits::bitCount(d_bits0) + bits::bitCount(d_bits1); }
  unsigned int firstBit() const;
  unsigned int lastBit() const;
  bool test(unsigned int j) const;
  unsigned int position(unsigned int j) const;
  bool scalarProduct(const BitSetBase<2>& b) const;

  wide_uint to_ulong() const
  { return d_bits0^(static_cast<wide_uint>(d_bits1)<<32u); }

  iterator begin() const;

// manipulators

  void operator^= (const BitSetBase<2>& b);
  void operator|= (const BitSetBase<2>& b);
  void operator&= (const BitSetBase<2>& b);
  void operator<<= (unsigned int c);
  void operator>>= (unsigned int c);
  void andnot(const BitSetBase<2>& b);
  void flip(unsigned int j);

  void reset() { d_bits0 = 0u; d_bits1 = 0u; }
  void reset(unsigned int j) ;
  void set(unsigned int j);
  void fill(unsigned int limit);

  void complement(unsigned int limit);
  void truncate(unsigned int limit);

  void slice(const BitSetBase<2>& c); // extract bits set in |c|, compacting
  void unslice(const BitSetBase<2>& c); // expand bits to positions set in |c|
  void swap(BitSetBase<2>& source);
 }; // |class BitSetBase<2>|


#ifdef constexpr // when this is a macro, dont use it
#define chunks_for(n) (((n) + 31)/32)
#else
 constexpr unsigned int chunks_for(unsigned int n) { return (n+31)/32; }
#endif




// the actual			|BitSet| class


/*
  Bitset of n bits.

  The class is derived from |BitSetBase<m>| for $m=\lceil n/32\rceil$
*/

template<unsigned int n> class BitSet
  : public BitSetBase<chunks_for(n)>
{
  typedef BitSetBase<chunks_for(n)> Base;

 public:

// associated types; there are only constant iterators
  struct iterator; // derived
  typedef iterator const_iterator;

// constructors and destructors
  BitSet() : Base() {}
  explicit BitSet(wide_uint b) : Base(b) {} // set at most |longBits| bits

  template<typename I> // integer type
    explicit BitSet(const std::vector<I>& v); // takes parity bit of each entry

// accessors

#if 0 // these are now publicly inherited methods
  bool operator== (const BitSet& b) const { return Base::operator== (b); }
  bool operator!= (const BitSet& b) const { return Base::operator!= (b); }
  bool operator<  (const BitSet& b) const { return Base::operator<  (b); }

  bool any() const { return Base::any(); }
  bool none() const { return Base::none(); }

  bool any(const BitSet& b) const { return Base::any(b); }
  bool contains(const BitSet& b) const { return Base::contains(b); }

  unsigned int count() const { return Base::count(); }
  unsigned int firstBit() const { return Base::firstBit(); }
  unsigned int lastBit() const { return Base::lastBit(); }

  bool test(unsigned int j) const { return Base::test(j); }

  // rank among the set bits of bit number |j| (assuming it were set)
  unsigned int position(unsigned int j) const { return Base::position(j); }

  bool scalarProduct(const BitSet& b) const { return Base::scalarProduct(b); }

  wide_uint to_ulong() const { return Base::to_ulong(); }
#endif

  iterator begin() const { return iterator(Base::begin()); }
  iterator end() const { return iterator(); } // zero value is end indicator

// manipulators

  BitSet& operator^= (const BitSet& b) { Base::operator^=(b); return *this; }
  BitSet& operator|= (const BitSet& b) { Base::operator|=(b); return *this; }
  BitSet& operator&= (const BitSet& b) { Base::operator&=(b); return *this; }

  BitSet& operator<<= (unsigned int c) { Base::operator<<=(c); return *this; }
  BitSet& operator>>= (unsigned int c) { Base::operator>>=(c); return *this; }

  BitSet& andnot (const BitSet& b) { Base::andnot(b); return *this; }
  BitSet& flip(unsigned int j) { Base::flip(j); return *this; }

  BitSet& reset()         { Base::reset();  return *this; }
  BitSet& reset(unsigned int j) { Base::reset(j); return *this; }
  BitSet& set(unsigned int j) { Base::set(j); return *this; }
  BitSet& set(unsigned int j, bool b)
    { if (b) Base::set(j); else Base::reset(j); return *this; }
  BitSet& fill(unsigned int limit) { Base::fill(limit); return *this; }
  BitSet& complement(unsigned int limit)
    { Base::complement(limit); return *this; }
  BitSet& truncate(unsigned int m) { Base::truncate(m); return *this; }

  BitSet& slice(const BitSet& c) { Base::slice(c); return *this; }
  BitSet& unslice(const BitSet& c) { Base::unslice(c); return *this; }

  void swap(BitSet& source) { Base::swap(source); }

  // accessors (non inherited)

  bool operator[] (unsigned int j) const { return Base::test(j); }

  // non-assignment logical operators added by MvL
  BitSet operator& (const BitSet& b) const // logical AND
    { BitSet t(*this); return t&=b; }
  BitSet operator| (const BitSet& b) const // logical OR
    { BitSet t(*this); return t|=b; }
  BitSet operator^ (const BitSet& b) const // logical XOR
    { BitSet t(*this); return t^=b; }
  BitSet operator- (const BitSet& b) const // logical difference (AND NOT)
    { BitSet t(*this); return t.andnot(b); }


  bool dot(BitSet b) const { return Base::scalarProduct(b); }

}; // |class BitSet|




//			    Iterator classes

// iterators claim just InputIterator status, rather than ForwardIterator:
// while multipass would be possible, no true |reference| type is used

// Iterator through the _set_ bits (like |BitMap::iterator|)
class BitSetBase<1>::iterator
  : public std::iterator<std::input_iterator_tag,unsigned int>
{
  chunk d_bits; // iterator contains a copy of the set iterated over

 public:

// constructors and destructors
  iterator() : d_bits(0u) {}
  explicit iterator(const BitSetBase<1>& b) : d_bits(b.d_bits) {}

// accessors
  bool operator== (iterator i) const { return d_bits==i.d_bits; }
  bool operator!= (iterator i) const { return d_bits!=i.d_bits; }
  unsigned int operator* () const { return bits::firstBit(d_bits); }
  bool operator() () const { return d_bits!=0; }

// manipulators
  iterator& operator++ () { d_bits &= d_bits-1; return *this; }
  iterator operator++ (int) { iterator tmp(*this); ++(*this); return tmp; }
}; // |class BitSetBase<1>::iterator|

class BitSetBase<2>::iterator
  : public std::iterator<std::input_iterator_tag,unsigned int>
{
  uint_least64_t d_bits; // pack bitset data into one unit for efficiency

 public:

// constructors and destructors
  iterator() : d_bits(0ul) {}
  explicit iterator(const BitSetBase<2>& b)
    : d_bits(b.d_bits0^(static_cast<wide_uint>(b.d_bits1)<<32)) {}

// accessors
  bool operator== (iterator i) const { return d_bits==i.d_bits; }
  bool operator!= (iterator i) const { return d_bits!=i.d_bits; }
  unsigned int operator* () const { return bits::firstBit(d_bits); }
  bool operator() () const { return d_bits!=0; }

// manipulators
  iterator& operator++ () { d_bits &= d_bits-1; return *this; }
  iterator operator++ (int) { iterator tmp(*this); ++(*this); return tmp; }
}; // |class BitSetBase<2>::iterator|

// |BitSet| iterator just inherits from the |BitSetBase| iterator type
template<unsigned int n>
  struct BitSet<n>::iterator : public BitSet<n>::Base::iterator
  {
    iterator() : Base::iterator() {} // zero (end) iterator
    explicit iterator(const typename Base::iterator& b) : Base::iterator(b) {}
  }; // |struct Bitset<n>::iterator|

} // |namespace bitset|

} // |namespace atlas|

#endif
