/*!
\file
\brief Class definitions and function declarations for the BitSet class.
*/
/*
  This is bitset.h.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef BITSET_H  /* guard against multiple inclusions */
#define BITSET_H

#include "bitset_fwd.h"

#include "bits.h"
#include "constants.h"
#include "setutils.h"

/******** type definitions ***************************************************/

namespace atlas {

namespace bitset {

/*****************************************************************************

  The |BitSet| class has essentially the functionality of the |stl::bitset
  class| from the STL, but I have redone it because I needed some things that
  bitset doesn't seem to provide. This is a very important data type for the
  program, and it is crucial that it should be efficiently implemented. [Fokko]

  An important difference with |stl::bitset| is that whereas the latter model
  set of exactly |n| bits, where |n| is the template parameter, |BitSet| will
  often be used to model sets of at most |n| bits. The exact size is not
  stored in the |BitSet| for efficiency reasons, and is therefore supplied
  explicitly to various methods; the class guarantees that if the same size is
  systematically used with an object, then no bits beyond that limit will be
  set. Then iterators of the set bits in the |BitSet| can be safely used. One
  consequence of this is that there can be no |operator~|, but rather a method
  |complement| is defined with a second |size| argument.

  Only fixed template instances are defined here, the template for general |n|
  is left undefined. The intention is to use |BitSet| for small bitmaps with
  size bounded at compile time, typically |constants::RankMax|. This allows
  for optimal speed, at the cost of some code duplication.

  Just as for the implementation of bitset in the STL that I took as a model,
  the template paramter |n| is passed to a base class template |BitSetBase|
  after conversion to number of |unsigned long|s needed rather than bits.
  Only a few instantiations of |BitSetBase| are provided, currently only the
  single-word vesrion is complete and the 2-word version is partially defined.

  The main difference with the STL is the |to_ulong()| method; mine will
  truncate the bitset when the number of bits is too long, instead of
  throwing an exception.

  No bit-references are implemented, and |operator[]| is defined as accessor,
  equivalent to |test()|. In other respects method naming conventions are not
  yet coherent with those of |bitmap::BitMap|, a desirable future change.

  The output operator is in io/ioutils.h.

******************************************************************************/

/*
   First a template class BitSetBase is defined, which provides the basic
   implementation (a small fixed number of unsigned long integers), but not
   yet the desired interface. Instances of the template class BitSet will be
   privately derived from instances of BitSetBase (with a different argument).
*/

template <size_t n> class BitSetBase;

  /*!
  \brief Base for a non-empty BitSet that fits in one word.

  With RANK_MAX=16, this template should be the only one instantiated
  on a 32-bit machine.

  The BitSet class BitSet<n>, for n between 1 and the machine word
  length (more precisely, the constant longBits), is a derived class
  of BitSetBase<1>.
  */
template<> class BitSetBase<1>
{
  /*!
  \brief Word that holds the BitSet
  */
  unsigned long d_bits;

 protected:

// associated types
  class iterator;

// constructors
           BitSetBase<1>()               : d_bits(0) {}
  explicit BitSetBase<1>(unsigned long b): d_bits(b) {}

  /*!
\brief Copies from another BitSet size (not BitSetBase) by copying the
  first word (there should not be more than one).
  */
  template<size_t m>
    explicit BitSetBase<1>(const BitSet<m>& b) : d_bits(b.to_ulong(0)) {}

// accessors

 public:
  bool operator==(const BitSetBase<1>& b) const { return d_bits == b.d_bits; }
  bool operator!=(const BitSetBase<1>& b) const { return d_bits != b.d_bits; }
  bool operator< (const BitSetBase<1>& b) const { return d_bits < b.d_bits; }

  bool any() const { return d_bits!=0; }
  bool none() const { return d_bits==0; }

  bool any(const BitSetBase<1>& b) const
  { return (d_bits & b.d_bits)!=0; }
  bool contains (const BitSetBase<1>& b) const
  { return (~d_bits & b.d_bits)==0; }

  size_t count() const { return bits::bitCount(d_bits); }
  size_t firstBit() const // index of least set bit; |longBits| if |none()|
  { return bits::firstBit(d_bits); }
  size_t lastBit() const // index of highest set bit PLUS ONE, 0 if |none()|
  { return bits::lastBit(d_bits); }

  bool test(size_t j) const { return (d_bits & constants::bitMask[j])!=0; }

  size_t position(size_t j) const
  { return bits::bitCount(d_bits & constants::lMask[j]); }

  bool scalarProduct(const BitSetBase<1>& b) const
  { return bits::bitCount(d_bits&b.d_bits)%2 != 0; }

  unsigned long to_ulong() const { return d_bits; }
  unsigned long to_ulong(size_t n) const { return n==0 ? d_bits : 0; }

  iterator begin() const;

// manipulators

 protected: // these cannot be called directly: no return value is available
  void operator^= (const BitSetBase<1>& b) { d_bits ^= b.d_bits; }
  void operator|= (const BitSetBase<1>& b) { d_bits |= b.d_bits; }
  void operator&= (const BitSetBase<1>& b) { d_bits &= b.d_bits; }
  void andnot(const BitSetBase<1>& b)      { d_bits &= ~b.d_bits; }

  // next two methods make sure that shift by |constants::longBits| yields 0.
  void operator<<= (size_t c)
  {
    if (c < constants::longBits) // bit shifts by more than this are undefined!
      d_bits <<= c;
    else
      d_bits = 0ul; // simulate shifting out of all bits
  }
  void operator>>= (size_t c)
  {
    if (c < constants::longBits) // bit shifts by more than this are undefined!
      d_bits >>= c;
    else
      d_bits = 0ul; // simulate shifting out of all bits
  }

  void flip(size_t j) { d_bits ^= constants::bitMask[j]; }
  void permute(const setutils::Permutation& a) { bits::permute(d_bits,a); }
  void reset() { d_bits = 0ul; }
  void reset(size_t j)  { d_bits &= ~constants::bitMask[j]; }
  void set(size_t j) { d_bits |= constants::bitMask[j]; }
  void fill(size_t limit) { d_bits |= constants::lMask[limit]; }

  void complement(size_t limit) { d_bits ^= constants::lMask[limit]; }
  void truncate(size_t limit)
    { if (limit<constants::longBits) d_bits &= constants::lMask[limit]; }

  void slice(const BitSetBase<1>& c); // extract bits set in |c|
  void swap(BitSetBase<1>& source) { std::swap(d_bits,source.d_bits); }

};

  /*!
  \brief Base for a non-empty BitSet that fits in two words but not one.

  The BitSet class BitSet<n>, for n between machine word length + 1
  and twice machine word length, is a derived class of BitSetBase<2>.
  Should not be instantiated on a 32 bit machine with RANK_MAX=16.
  */
template<> class BitSetBase<2>
{
  unsigned long d_bits0,d_bits1;;

 protected:

// associated types
  class iterator;

// constructors and destructors
  BitSetBase<2>() { d_bits0 = 0ul; d_bits1 = 0ul; }

  /*!
  \brief Constructor initializing first word to b and second word to 0.

  [added by DV to let the software compile with RANK_MAX equal to the
  machine word size.]

  The class BitSet assumes that BitSetBase has a constructor with argument an
  unsigned long. This is slightly sloppy coding, since BitSetBase<2> is most
  naturally constructed using two unsigned longs. However such a constructor
  will never be called from |BitSet|, so it would be useless. In fact it seems
  that the current software never needs to construct from an explicit value a
  |BitSet| that needs more than a single |unsigned long|, even if one should
  set |constants::RankMax == constants::longBits|; while for instance
  |gradings::Status::set(size_t,Value)| constructs a |bitset::TwoRankFlags|
  from an |unsigned long|, the latter is in fact only 2 bits wide, which value
  is shifted in place after construction by |bitset::TwoRankFlags::operator<<=|.
  */
  explicit BitSetBase<2>(unsigned long b) { d_bits0 = b; d_bits1 = 0ul; }

// copy constructor, possibly from from other (shorter) size
  template<size_t m>
    explicit BitSetBase<2>(const BitSet<m>& b)
  {
    d_bits0 = b.to_ulong(0);
    d_bits1 = b.to_ulong(1);
  }

// accessors

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
    { return (b.d_bits0 & ~d_bits0)==0 and (b.d_bits1 & ~d_bits1)==0; }

  size_t count() const
    { return bits::bitCount(d_bits0) + bits::bitCount(d_bits1); }
  size_t firstBit() const;
  size_t lastBit() const;
  bool test(size_t j) const;
  size_t position(size_t j) const;
  bool scalarProduct(const BitSetBase<2>& b) const;

  unsigned long to_ulong() const { return d_bits0; }
  unsigned long to_ulong(size_t n) const
    { return n==0 ? d_bits0 : n==1 ? d_bits1 : 0; }

  iterator begin() const;

// manipulators

  void operator^= (const BitSetBase<2>& b);
  void operator|= (const BitSetBase<2>& b);
  void operator&= (const BitSetBase<2>& b);
  void operator<<= (size_t c);
  void operator>>= (size_t c);
  void andnot(const BitSetBase<2>& b);
  void flip(size_t j);
  void permute(const setutils::Permutation& a);

  void reset() { d_bits0 = 0ul; d_bits1 = 0ul; }
  void reset(size_t j) ;
  void set(size_t j);
  void fill(size_t limit);

  void complement(size_t limit);
  void truncate(size_t limit);

  void slice(const BitSetBase<1>& c);
  void swap(BitSetBase<2>& source);
 }; // |class BitSetBase<2>|


/*! \brief The class BaseSize computes (with its member 'value') the base size
  - the number of words needed for a BitSet holding n bits. Since n must be a
  compile time constant, BaseSize<n>::value will be one as well.

  Essentially we must divide n by longBits, but the fractional result
  must be rounded up to the next integer; this is achieved by adding
  constants::posBits=longBits-1 to n before performing the division.
  Code simplified by MvL.

*/
template<size_t n> struct BaseSize
{
  static const size_t value = (n + constants::posBits)/constants::longBits;
};

// the actual BitSet class
/*!
  \brief Set of n bits.

  The first typedef defines Base to be BitSetBase<m>. Consequently function
  references of the form Base::contains refer to BitSetBase<m>::contains.
  [DV doesn't know whether this could have been avoided by replacing
  ":private BitSetBase" by ":public BitSetBase" in the template for BitSet,
  and just invoking the (public) member functions of the base class directly.
  MvL thinks that that would be possible; this depends on the fact (also used
  in the actual implementation below) that the argument of 'contains' is
  implicitly converted to its base class BitSetBase<m> because the method of
  the base class requires this. However, this would expose _all_ public
  methods of BitSetBase to users of BitSet, which is not desired. This could
  however be remedied by making protected rather than public the methods of
  BitSetBase that are only for internal use by BitSet implementations.]
*/
template<size_t n> class BitSet
  : public BitSetBase<BaseSize<n>::value>
{
  typedef BitSetBase<BaseSize<n>::value> Base;

 public:

// associated types; there are only constant iterators
  struct iterator : public Base::iterator
  {
    // associated types

    typedef std::forward_iterator_tag iterator_category;
    typedef size_t value_type;

    iterator() : Base::iterator() {}
    explicit iterator(const typename Base::iterator& b) : Base::iterator(b) {}

  };
  typedef iterator const_iterator;

// constructors and destructors
  BitSet() : Base() {}
  explicit BitSet(unsigned long b) : Base(b) {}

  template<typename I> // integer type
    explicit BitSet(const std::vector<I>& v); // takes parity bit of each entry

  //! \brief Copy from other size BitSets, only to be used with |m<n|
  template<size_t m> BitSet(const BitSet<m>& b) : Base(b) {}

// accessors

#if 0 // these are now publicly inherited methods
  bool operator== (const BitSet& b) const { return Base::operator== (b); }
  bool operator!= (const BitSet& b) const { return Base::operator!= (b); }
  bool operator<  (const BitSet& b) const { return Base::operator<  (b); }

  bool any() const { return Base::any(); }
  bool none() const { return Base::none(); }

  bool any(const BitSet& b) const { return Base::any(b); }
  bool contains(const BitSet& b) const { return Base::contains(b); }

  size_t count() const { return Base::count(); }
  size_t firstBit() const { return Base::firstBit(); }
  size_t lastBit() const { return Base::lastBit(); }

  bool test(size_t j) const { return Base::test(j); }

  size_t position(size_t j) const { return Base::position(j); }

  bool scalarProduct(const BitSet& b) const { return Base::scalarProduct(b); }

  unsigned long to_ulong() const { return Base::to_ulong(); }
  unsigned long to_ulong(size_t i) const { return Base::to_ulong(i); }
#endif

  iterator begin() const { return iterator(Base::begin()); }

// manipulators

  BitSet& operator^= (const BitSet& b) { Base::operator^=(b); return *this; }
  BitSet& operator|= (const BitSet& b) { Base::operator|=(b); return *this; }
  BitSet& operator&= (const BitSet& b) { Base::operator&=(b); return *this; }

  BitSet& operator<<= (size_t c) { Base::operator<<=(c); return *this; }
  BitSet& operator>>= (size_t c) { Base::operator>>=(c); return *this; }

  BitSet& andnot (const BitSet& b) { Base::andnot(b); return *this; }
  BitSet& flip(size_t j) { Base::flip(j); return *this; }

  BitSet& permute(const setutils::Permutation& a)
    { Base::permute(a); return *this; }

  BitSet& reset()         { Base::reset();  return *this; }
  BitSet& reset(size_t j) { Base::reset(j); return *this; }
  BitSet& set(size_t j) { Base::set(j); return *this; }
  BitSet& set(size_t j, bool b)
    { if (b) Base::set(j); else Base::reset(j); return *this; }
  BitSet& fill(size_t limit) { Base::fill(limit); return *this; }
  BitSet& complement(size_t limit) { Base::complement(limit); return *this; }
  BitSet& truncate(size_t m) { Base::truncate(m); return *this; }

  BitSet& slice(const BitSet& c) { Base::slice(c); return *this; }

  void swap(BitSet& source) { Base::swap(source); }

  // accessors (non inherited)

  bool operator[] (size_t j) const { return Base::test(j); }

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

}; // |template<size_t n> class BitSet|



//! \brief Iterator through the _set_ bits (like |BitMap::iterator|)
class BitSetBase<1>::iterator
{
  unsigned long d_bits; // iterator contains a copy of the set iterated over

 public:

// constructors and destructors
  iterator() : d_bits(0ul) {}
  explicit iterator(const BitSetBase<1>& b) : d_bits(b.to_ulong()) {}

// accessors
  bool operator== (const iterator& i) const; // rarely useful
  bool operator!= (const iterator& i) const; // rarely useful
  size_t operator* () const { return bits::firstBit(d_bits); }
  bool operator() () const { return d_bits!=0; }

// manipulators
  iterator& operator++ () { d_bits &= d_bits-1; return *this; }
  iterator operator++ (int) { iterator tmp(*this); ++(*this); return tmp; }
}; // |class BitSetBase<1>::iterator|

class BitSetBase<2>::iterator
{
  unsigned long d_bits0,d_bits1; // copy of bitset data

 public:

// constructors and destructors
  iterator() : d_bits0(0), d_bits1(0) {}
  explicit iterator(const BitSetBase<2>& b)
    : d_bits0(b.to_ulong()), d_bits1(b.to_ulong(1)) {}

// accessors
  bool operator== (const iterator& i) const;
  bool operator!= (const iterator& i) const;
  size_t operator* () const;
  bool operator() () const { return d_bits0!=0 or d_bits1!=0; }

// manipulators
  iterator& operator++ ();
  iterator operator++ (int) { iterator tmp(*this); ++(*this); return tmp; }
}; // |class BitSetBase<2>::iterator|

} // |namespace bitset|

} // |namespace atlas|

#endif
