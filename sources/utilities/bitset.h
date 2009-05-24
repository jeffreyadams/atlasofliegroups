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

/******** function declarations **********************************************/

namespace atlas {

namespace bitset {

  template<size_t n> BitSet<n> operator~ (const BitSet<n>& d_b);
  template<size_t n> void set(BitSet<n>&, size_t);

}

/******** type definitions ***************************************************/

namespace bitset {

/*****************************************************************************

  The BitSet class has essentially the functionality of the stl::bitset class
  from the STL, but I have redone it because I needed some things that
  bitset doesn't seem to provide. This is a very important data type for
  the program, and it is crucial that it should be efficiently implemented.

  The main difference is that _all_ my BitSets are fixed size; they are
  provided for up to 128 bits on 32-bit machines, 256 bits on 64-bit machines.
  (This statement is confusing, since stl::bitset objects are fixed size as
  well; however, they are defined by a general template, while here only some
  concrete template instances are defined. MvL). So the intention is to use
  them for small bitmaps. This allows for optimal speed, at the cost of some
  code duplication (but almost all the functions in the class are easily
  inlined.)

  Just as for the implementation of bitset in the STL that I took as a model,
  the size is rounded of to multiples of the number of bits in an unsigned
  long, so that in effect there are only four instances of the BitSetBase
  to provide (or five, counting the 0-size BitSet.)

  The main difference with the STL is the to_ulong() method; mine will
  truncate the bitset when the number of bits is too long, instead of
  throwing an exception.

  Although I had implemented bit-references, and even though it is an
  interesting C++ trick that shows to what extent you can twist the notation
  to do what you want, I dismissed it now as syntactic sugar. Use
  b.set(i,value) instead of b[i] = value. Operator[] is now defined only as
  an access operator, equivalent to test().

  The output operator is in io/ioutils.h.

******************************************************************************/

/* First a template class BitSetBase is defined, which provides the basic
   implementation (a small fixed number of unsigned long integers), but not
   yet the desired interface. Instances of the template class BitSet will be
   privately derived from instances of BitSetBase (with a different argument).
*/

// this is not supposed to ever be instantiated
template<size_t n> class BitSetBase {
};

  /*!
  \brief  Intended as Base for an empty BitSet.  Fill in if necessary.
  */
template<> class BitSetBase<0> {
};

  /*!
  \brief Base for a non-empty BitSet that fits in one word.

  With RANK_MAX=16, this template should be the only one instantiated
  on a 32-bit machine.

  The BitSet class BitSet<n>, for n between 1 and the machine word
  length (more precisely, the constant longBits), is a derived class
  of BitSetBase<1>.
  */
template<> class BitSetBase<1> {

 private:
  /*!
  \brief Word that holds the BitSet
  */
  unsigned long d_bits;

 protected: // not public, since this class should only be used via BitSet

// associated types
  class iterator;

// constructors and destructors
  BitSetBase<1>():d_bits(0ul) {}

  explicit BitSetBase<1>(unsigned long b)
    :d_bits(b) {}

  ~BitSetBase<1>() {}

// "copy" constructor from BitSet (not Base!) of any size

  template<size_t m>
  /*!
\brief Copies from another bitset size by truncating at the end of the
  first word.

  For each BitSet size, to_ulong is the word holding the first
  longBits (machine word size) of the BitSet.  This copy constructor
  therefore truncates b at the end of its first word.
  */
    explicit BitSetBase<1>(const BitSet<m>& b)
    :d_bits(b.to_ulong()) {}

// assignment from other sizes
  /*!
\brief Assigns from another bitset size by truncating at the end of
  the first word.

  For each BitSet size, to_ulong is the word holding the first
  longBits (machine word size) of the BitSet.  This assignment
  operator therefore truncates b at the end of its first word.
  */
  template<size_t m>
    BitSetBase<1>& operator= (const BitSet<m>& b) {
    d_bits = b.to_ulong();
    return *this;
  }

// accessors

  bool operator==(const BitSetBase<1>& b) const {
    return d_bits == b.d_bits;
  }

  bool operator!=(const BitSetBase<1>& b) const {
    return d_bits != b.d_bits;
  }

  bool operator< (const BitSetBase<1>& b) const {
    return d_bits < b.d_bits;
  }
  /*!
\brief  Returns 1 if any bit of the BitSet is 1, and 0 otherwise.
  */
  bool any() const { return d_bits!=0; }

  bool any(const BitSetBase<1>& b) const { return (d_bits & b.d_bits)!=0; }

  iterator begin() const;

  /*!
\brief Tests whether this bitset contains b.

  Returns whether if every set bit of b is also set in this BitSet.
  */
  bool contains (const BitSetBase<1>& b) const {
    return (b.d_bits & ~d_bits)==0;
  }
  /*!
\brief Number of set bits in BitSet.
  */
  size_t count() const {
    return bits::bitCount(d_bits);
  }

  size_t firstBit() const // index of least set bit; |longBits| if |none()|
  { return bits::firstBit(d_bits); }

  size_t lastBit() const // index of highest set bit PLUS ONE, 0 if |none()|
  { return bits::lastBit(d_bits); }

  bool none() const { return d_bits==0; }

     /*!
\brief Number of set bits in position < j.

If j is itself set, this is the position of bit j among the SET bit
   */
  size_t position(size_t j) const
  { return bits::bitCount(d_bits & constants::lMask[j]); }

  bool scalarProduct(const BitSetBase<1>& b) const;

  //! \brief Returns the value of bit j in the BitSet.
  bool test(size_t j) const { return (d_bits & constants::bitMask[j])!=0; }

  unsigned long to_ulong() const { return d_bits; }

  unsigned long to_ulong1() const { return 0ul; } // second unsigned long slice

// manipulators

  BitSetBase<1>& operator^= (const BitSetBase<1>& b) {
    d_bits ^= b.d_bits;
    return *this;
  }

  BitSetBase<1>& operator|= (const BitSetBase<1>& b) {
    d_bits |= b.d_bits;
    return *this;
  }

  BitSetBase<1>& operator&= (const BitSetBase<1>& b) {
    d_bits &= b.d_bits;
    return *this;
  }

  // We have to make sure that shift by constants::longBits yields 0.
  BitSetBase<1>& operator<<= (size_t c)
  {
    if (c < constants::longBits) // bit shifts by more than this are undefined!
      d_bits <<= c;
    else
      d_bits = 0ul; // simulate shifting out of all bits
    return *this;
  }

  // We have to make sure that shift by constants::longBits yields 0.
  BitSetBase<1>& operator>>= (size_t c)
  {
    if (c < constants::longBits) // bit shifts by more than this are undefined!
      d_bits >>= c;
    else
      d_bits = 0ul; // simulate shifting out of all bits
    return *this;
  }

  BitSetBase<1>& andnot(const BitSetBase<1>& b) {
    d_bits &= ~b.d_bits;
    return *this;
  }

  BitSetBase<1>& flip() {
    d_bits = ~d_bits;
    return *this;
  }

  BitSetBase<1>& flip(size_t j) {
    d_bits ^= constants::bitMask[j]; // lifted from STL !
    return *this;
  }

  BitSetBase<1>& permute(const setutils::Permutation& a) {
    bits::permute(d_bits,a);
    return *this;
  }

  BitSetBase<1>& reset() {
    d_bits = 0ul;
    return *this;
  }

  BitSetBase<1>& reset(size_t j) {
    d_bits &= ~constants::bitMask[j];
    return *this;
  }

  BitSetBase<1>& set() {
    d_bits = ~(0ul);
    return *this;
  }

  BitSetBase<1>& set(size_t j) {
    d_bits |= constants::bitMask[j];
    return *this;
  }

  BitSetBase<1>& set(size_t j, bool b) {
    if (b)
      set(j);
    else
      reset(j);
    return *this;
  }

  BitSetBase<1>& slice(const BitSetBase<1>& c);

  void swap(BitSetBase<1>& source)
  {
    std::swap(d_bits,source.d_bits);
  }

  BitSetBase<1>& truncate(size_t m)
  {
    if (m < constants::longBits)
      d_bits &= constants::lMask[m];
    return *this;
  }

};

  /*!
\brief Iterator through the _set_ bits in the container.

This is the same behaviour as BitMap::iterator.
  */
class BitSetBase<1>::iterator {

 private:

  unsigned long d_bits;

 public:

// associated types

  typedef std::forward_iterator_tag iterator_category;
  typedef size_t value_type;
  typedef ptrdiff_t difference_type;
  typedef const value_type* pointer;
  typedef const value_type& reference;

// constructors and destructors
  iterator()
    :d_bits(0ul) {}

  explicit iterator(unsigned long b)
    :d_bits(b) {}

  ~iterator() {}

// accessors
  bool operator== (const iterator& i) const {
    return bits::firstBit(d_bits) == bits::firstBit(i.d_bits);
  }

  bool operator!= (const iterator& i) const{
    return bits::firstBit(d_bits) != bits::firstBit(i.d_bits);
  }

  size_t operator* () const {
    return bits::firstBit(d_bits);
  }

  bool operator() () const {
    return d_bits!=0;
  }

// manipulators
  iterator& operator++ () {
    d_bits &= d_bits-1; // K-R : erase the first set bit
    return *this;
  }

  iterator operator++ (int) {
    iterator tmp(*this);
    ++(*this);
    return tmp;
  }

}; // |class BitSetBase<1>::iterator|

  /*!
  \brief Base for a non-empty BitSet that fits in two words but not one.

  The BitSet class BitSet<n>, for n between machine word length + 1
  and twice machine word length, is a derived class of BitSetBase<2>.
  Should not be instantiated on a 32 bit machine with RANK_MAX=16.
  */
template<> class BitSetBase<2> {
 private:
  /*!
\brief Array of two words that holds the BitSet.
  */
  unsigned long d_bits[2];

 public:

// associated types
  class iterator;

// constructors and destructors
  BitSetBase<2>() {
    d_bits[0] = 0ul;
    d_bits[1] = 0ul;
  }


  /*!
  \brief Constructor initializing first word to b and second word to 0.

  [added by DV to let the software compile with RANK_MAX equal to the
  machine word size.]

  The class BitSet assumes that BitSetBase has a constructor with
  argument an unsigned long.  This is slightly sloppy coding, since
  BitSetBase<2> is most naturally constructed using an array of two
  unsigned longs (as in the next constructor).  This constructor is
  actually used in the present software only to initialize the first
  one or two bits (for example in the definition of the manipulator
  Status::set(size_t,Value) in gradings.h).  So this crude
  construction seems to work.
  */
  explicit BitSetBase<2>(unsigned long b) {
    d_bits[0] = b;
    d_bits[1] = 0ul;
  }


  explicit BitSetBase<2>(unsigned long b[2]) {
    d_bits[0] = b[0];
    d_bits[1] = b[1];
  }

  ~BitSetBase<2>() {}

// copy constructors from other sizes

  template<size_t m>
    explicit BitSetBase<2>(const BitSet<m>& b) {
    d_bits[0] = b.to_ulong();
    d_bits[1] = b.to_ulong1();
  }

// assignment from other sizes

  template<size_t m>
    BitSetBase<2>& operator= (const BitSet<m>& b) {
    d_bits[0] = b.to_ulong();
    d_bits[1] = b.to_ulong1();
    return *this;
  }

// accessors

  bool operator==(const BitSetBase<2>& b) const {
    return d_bits[0] == b.d_bits[0] and d_bits[1] == b.d_bits[1];
  }

  bool operator!=(const BitSetBase<2>& b) const {
    return d_bits[0] != b.d_bits[0] or d_bits[1] != b.d_bits[1];
  }

  bool operator< (const BitSetBase<2>& b) const
  {
    return d_bits[1] == b.d_bits[1]
      ? d_bits[0] < b.d_bits[0]
      : d_bits[1] < b.d_bits[1];
  }

  bool any() const {
    return d_bits[0]!=0 or d_bits[1]!=0;
  }

  bool any(const BitSetBase<2>& b) const {
    return (d_bits[0] & b.d_bits[0])!=0 or (d_bits[1] & b.d_bits[1])!=0;
  }


  /*!
  \brief Not yet implemented.
  */
  iterator begin() const;

  bool contains (const BitSetBase<2>& b) const {
    return (b.d_bits[0] & ~d_bits[0])==0 and (b.d_bits[1] & ~d_bits[1])==0;
  }

  size_t count() const {
    return bits::bitCount(d_bits[0]) + bits::bitCount(d_bits[1]);
  }

  size_t firstBit() const {
    if (d_bits[0]!=0)
      return bits::firstBit(d_bits[0]);
    else
      return bits::firstBit(d_bits[1]) + constants::longBits;
  }

  size_t lastBit() const {
    if (d_bits[1]!=0)
      return bits::lastBit(d_bits[1]) + constants::longBits;
    else
      return bits::lastBit(d_bits[0]);
  }

  bool none() const { return d_bits[0]==0 and d_bits[1]==0; }

   /*!
\brief Number of set bits in position < j.

If j itself is set, this is the position of bit j among the set bits.
   */
  size_t position(size_t j) const {
    if (j >> constants::baseShift !=0) // two terms
      return bits::bitCount(d_bits[1] &
			    constants::lMask[j & constants::posBits])
	+ bits::bitCount(d_bits[0]);
    else // one term
      return bits::bitCount(d_bits[0] & constants::lMask[j]);
  }

  bool scalarProduct(const BitSetBase<2>& b) const;

  bool test(size_t j) const {
    return (d_bits[j >> constants::baseShift]
           & constants::bitMask[j & constants::posBits])!=0;
  }

  unsigned long to_ulong() const { return d_bits[0]; }

  unsigned long to_ulong1() const { return d_bits[1]; }

// manipulators

  BitSetBase<2>& operator^= (const BitSetBase<2>& b) {
    d_bits[0] ^= b.d_bits[0];
    d_bits[1] ^= b.d_bits[1];
    return *this;
  }

  BitSetBase<2>& operator|= (const BitSetBase<2>& b) {
    d_bits[0] |= b.d_bits[0];
    d_bits[1] |= b.d_bits[1];
    return *this;
  }

  BitSetBase<2>& operator&= (const BitSetBase<2>& b) {
    d_bits[0] &= b.d_bits[0];
    d_bits[1] &= b.d_bits[1];
    return *this;
  }

  BitSetBase<2>& operator<<= (size_t c) {
    if (c == 0) // do nothing
      return *this;
    if (c < constants::longBits) {
      unsigned long f = ~constants::lMask[constants::longBits-c];
      d_bits[1] <<= c;
      // copy top of d_bits[0] onto bottom of d_bits[1]
      d_bits[1] |= ((d_bits[0]&f) >> (constants::longBits - c));
      d_bits[0] <<= c;
    } else if (c == constants::longBits) {
      d_bits[1] = d_bits[0];
      d_bits[0] = 0ul;
    } else if (c < 2*constants::longBits) {
      d_bits[1] = d_bits[0] << (c - constants::longBits);
      d_bits[0] = 0ul;
    } else { // everything is shifted out
      d_bits[1] = 0ul;
      d_bits[0] = 0ul;
    }
    return *this;
  }

  BitSetBase<2>& operator>>= (size_t c) {
   /*!
  We have to make sure that shift by constants::longBits yields 0.
  */
    if (c == 0) // do nothing
      return *this;
    if (c < constants::longBits) {
      d_bits[0] >>= c;
  /*!
  Copy bottom of d_bits[1] onto top of d_bits[0].
  */
      d_bits[0] |= ((d_bits[1] & constants::lMask[c]) <<
		    (constants::longBits - c));
      d_bits[1] >>= c;
    } else if (c == constants::longBits) {
      d_bits[0] = d_bits[1];
      d_bits[1] = 0ul;
    } else if (c < 2*constants::longBits) {
      d_bits[0] = d_bits[1] >> (c - constants::longBits);
      d_bits[1] = 0ul;
    } else { // everything is shifted out
      d_bits[0] = 0ul;
      d_bits[1] = 0ul;
    }
    return *this;
  }

  BitSetBase<2>& andnot(const BitSetBase<2>& b) {
    d_bits[0] &= ~b.d_bits[0];
    d_bits[1] &= ~b.d_bits[1];
    return *this;
  }

  BitSetBase<2>& flip() {
    d_bits[0] = ~d_bits[0];
    d_bits[1] = ~d_bits[1];
    return *this;
  }

  BitSetBase<2>& flip(size_t j) {
    d_bits[j >> constants::baseShift]
      ^= constants::bitMask[j & constants::posBits]; // lifted from STL!
    return *this;
  }

  BitSetBase<2>& permute(const setutils::Permutation& a); // not yet implemented

  BitSetBase<2>& reset() {
    d_bits[0] = 0ul;
    d_bits[1] = 0ul;
    return *this;
  }

  BitSetBase<2>& reset(size_t j) {
    d_bits[j >> constants::baseShift]
      &= ~constants::bitMask[j & constants::posBits];
    return *this;
  }

  BitSetBase<2>& set() {
    d_bits[0] = ~0ul;
    d_bits[1] = ~0ul;
    return *this;
  }

  BitSetBase<2>& set(size_t j) {
    d_bits[j >> constants::baseShift]
      |= constants::bitMask[j & constants::posBits];
    return *this;
  }

  BitSetBase<2>& set(size_t j, bool b) {
    if (b)
      set(j);
    else
      reset(j);
    return *this;
  }

  BitSetBase<2>& slice(const BitSetBase<1>& c);

  void swap(BitSetBase<2>& source) {
    std::swap(d_bits[0],source.d_bits[0]);
    std::swap(d_bits[1],source.d_bits[1]);
  }

  BitSetBase<2>& truncate(size_t m) {
    if (m < constants::longBits)
    {
      d_bits[0] &= constants::lMask[m];
      d_bits[1] = 0ul;
    }
    else if (m < 2*constants::longBits)
      d_bits[1] &= constants::lMask[m & constants::posBits];

    return *this;
  }

};

  /*!
\brief Iterator runs through the _set_ bits in the container.

This is the same behaviour as BitMap::iterator.
  */
class BitSetBase<2>::iterator {
 private:

  // in this case, the most economical way seems to copy the full data
  BitSetBase<2> d_data;

 public:

// associated types

  typedef std::forward_iterator_tag iterator_category;
  typedef size_t value_type;
  typedef ptrdiff_t difference_type;
  typedef const value_type* pointer;
  typedef const value_type& reference;

// constructors and destructors
  iterator() {}

  explicit iterator(const BitSetBase<2>& b):d_data(b) {}

  ~iterator() {}

// accessors
  bool operator== (const iterator& i) const {
    return d_data.firstBit() == i.d_data.firstBit();
  }

  bool operator!= (const iterator& i) const{
    return d_data.firstBit() != i.d_data.firstBit();
  }

  size_t operator* () const {
    return d_data.firstBit();
  }

  bool operator() () const {
    return d_data.any();
  }

// manipulators
  iterator& operator++ () {
    d_data.reset(d_data.firstBit());
    return *this;
  }

  iterator operator++ (int) {
    iterator tmp(*this);
    ++(*this);
    return tmp;
  }
};

// fill in the missing ones when they will be needed

template<> class BitSetBase<3> {
};

template<> class BitSetBase<4> {
};

/*! \brief The class BaseSize computes (with its member 'value') the base size
  - the number of words needed for a BitSet holding n bits. Since n must be a
  compile time constant, BaseSize<n>::value will be one as well.

  Essentially we must divide n by longBits, but the fractional result
  must be rounded up to the next integer; this is achieved by adding
  constants::posBits=longBits-1 to n before performing the division.
  Code simplified by MvL.

*/
template<size_t n> class BaseSize {
  public:
    static const size_t value = (n + constants::posBits)/constants::longBits;
};

// the actual BitSet class
/*!
  \brief Set of n bits.

  The software formally allows for n between 0 and four times
  the machine word length; now up to two times the machine word length
  is implemented.

  The BitSet class has essentially the functionality of the bitset class
  from the STL, but I have redone it because I needed some things that
  bitset doesn't seem to provide. This is a very important data type for
  the program, and it is crucial that it should be efficiently implemented.

  The main difference is that _all_ my BitSets are fixed size; they
  are provided for up to 128 bits on 32-bit machines, 256 bits on
  64-bit machines.  So the intention is to use them for small
  bitmaps. This allows for optimal speed, at the cost of some code
  duplication (but almost all the functions in the class are easily
  inlined.)  [In fact BitSets of size more than twice the word length
  are not fully implemented in the present software.  The gap is that
  BitSet<n> calls on BitSetBase<m>, where m is the smallest integer
  greater than or equal to n/longBits.  (The size m is computed from n
  by the function BaseSize::value.)  The class BitSetBase<m> is
  actually implemented only for m = 1, 2. DV]

  Just as for the implementation of bitset in the STL that I took as a model,
  the size is rounded of to multiples of the number of bits in an unsigned
  long, so that in effect there are only four instances of the BitSetBase
  to provide (or five, counting the 0-size BitSet.)

  The main difference with the STL is the to_ulong() method; mine will
  truncate the bitset when the number of bits is too long, instead of
  throwing an exception.

  Although I had implemented bit-references, and even though it is an
  interesting C++ trick that shows to what extent you can twist the notation
  to do what you want, I dismissed it now as syntactic sugar. Use
  b.set(i,value) instead of b[i] = value. Operator[] is now defined only as
  an access operator, equivalent to test().

  The output operator is in io/ioutils.h.

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
  :private BitSetBase<BaseSize<n>::value>
{
 private:
  /*!
  \brief Class from which BitSet<n> is derived.
   */
    typedef BitSetBase<BaseSize<n>::value> Base;

 public:

// associated types; there are only constant iterators
  typedef typename Base::iterator iterator;
  typedef iterator const_iterator;

// constructors and destructors
  BitSet() {}

  explicit BitSet(unsigned long b):Base(b) {}

  template<typename I> // integer type
    explicit BitSet(const std::vector<I>& v);

  ~BitSet() {}
  /*!
  Copy from other size BitSets is by truncation or extension by zero.
  */
  template<size_t m>
    BitSet(const BitSet<m>& b) : Base(b) { truncate(n); }
  /*!
  Assignment from other size BitSets is by truncation or extension by zero.
  */
  template<size_t m>
    BitSet& operator= (const BitSet<m>& b)
    {
      Base::operator= (b);
      truncate(n);
      return *this;
    }

// accessors

  bool operator== (const BitSet& b) const { return Base::operator== (b); }

  bool operator!= (const BitSet& b) const { return Base::operator!= (b); }

  bool operator< (const BitSet& b) const { return Base::operator< (b); }

  bool operator[] (size_t j) const { return Base::test(j); }

  // non-assignment logical operators added by MvL
  BitSet operator& (const BitSet& b) const // logical AND
    { BitSet t(*this); return t&=b; }
  BitSet operator| (const BitSet& b) const // logical OR
    { BitSet t(*this); return t|=b; }
  BitSet operator^ (const BitSet& b) const // logical XOR
    { BitSet t(*this); return t^=b; }


  /*!
\brief Tests whether BitSet is nonempty.

Returns whether any bit of the BitSet is 1
  */
  bool any() const { return Base::any(); }

  /*!
\brief Tests whether BitSet and b have non-empty intersection.

Returns whether any bit is 1 in both BitSet and b.
  */
  bool any(const BitSet& b) const { return Base::any(b); }

  /*!
  Iterator pointing to the first set bit of BitSet.

Seems not yet to be implemented for BitSetBase<2>.
  */
  iterator begin() const { return Base::begin(); }

  /*!
\brief Tests whether BitSet contains b.
  */
  bool contains(const BitSet& b) const { return Base::contains(b); }

  /*!
\brief Number of set bits in BitSet.

This is the cardinality of the corresponding set.
  */
  size_t count() const { return Base::count(); }

  /*!
\brief Position of the first set bit in BitSet.

Returns the capacity of the BitSet if there is no such bit.
  */
  size_t firstBit() const { return Base::firstBit(); }

  /*!
\brief position of the last set bit PLUS ONE in BitSet, or 0 if |none()|
  */
  size_t lastBit() const { return Base::lastBit(); }

  /*!
\brief Tests whether the corresponding set is empty.

  Returns 1 if no bit of BitSet is set, and 0 otherwise.
  */
  bool none() const { return Base::none(); }

   /*!
\brief Number of set bits in position < j.

If j is set, this  is the position of j in the set of set bits.
   */
  size_t position(size_t j) const { return Base::position(j); }

  bool scalarProduct(const BitSet& b) const { return Base::scalarProduct(b); }
  bool dot(BitSet b) const { return Base::scalarProduct(b); }

  size_t size() const { return n; }

  /*!
\brief Tests bit j of BitSet.
  */
  bool test(size_t j) const { return Base::test(j); }

  /*!
\brief  The first wordlength bits of BitSet, interpreted as an unsigned
  long integer.
  */
  unsigned long to_ulong() const { return Base::to_ulong(); }

  /*!
\brief  The second wordlength bits of BitSet, interpreted as an unsigned
  long integer.
  */
  unsigned long to_ulong1() const { return Base::to_ulong1(); }

// manipulators

  BitSet& operator^= (const BitSet& b) { Base::operator^=(b); return *this; }
  BitSet& operator|= (const BitSet& b) { Base::operator|=(b); return *this; }
  BitSet& operator&= (const BitSet& b) { Base::operator&=(b); return *this; }

  BitSet& operator<<= (size_t c) { Base::operator<<=(c); return *this; }
  BitSet& operator>>= (size_t c) { Base::operator>>=(c); return *this; }

//! Performs a bitwise "and not" of this BitSet with the argument BitSet |b|
  BitSet& andnot (const BitSet& b) { Base::andnot(b); return *this; }
  BitSet& flip() { Base::flip(); truncate(n); return *this; }
  BitSet& flip(size_t j) { Base::flip(j); return *this; }

  /*!
\brief Applies the permutation a to the BitSet.
  */
  BitSet& permute(const setutils::Permutation& a)
  { Base::permute(a); return *this; }

  /*!
\brief Sets every bit of BitSet to zero.  (Empties the set.)
  */
  BitSet& reset()         { Base::reset();  return *this; }
  BitSet& reset(size_t j) { Base::reset(j); return *this; }

  /*!
\brief  Sets the first n bits of BitSet to one.  (Fills the set.)
  */
  BitSet& set()         { Base::set(); truncate(n); return *this; }
  BitSet& set(size_t j) { Base::set(j); return *this; }
  BitSet& set(size_t j, bool b) { Base::set(j,b); return *this; }

  /*!
  \brief Replaces the bitset by the "defragmented" intersection
  with c (i.e., the bits of that intersection are written
  consecutively.)

  More precisely: if c has m set bits, then the BitSet is replaced by
  putting in its first m bits the values previously found in the
  positions flagged by c.  The remaining bits are made 0.

  This is used for linear algebra, as follows.  Suppose a subspace of
  (Z/2Z)^n is given its unique basis in row-reduced form, and that the
  leading bits of the basis vectors are flagged by c.  (This is how
  subspaces are represented in the NormalSubspace class.)  If the
  BitVector belongs to this subspace, then applying slice(c) puts in
  the BitVector its coordinates with respect to the row-reduced basis.
  */
  BitSet& slice(const BitSet& c) { Base::slice(c); return *this; }

  void swap(BitSet& source) { Base::swap(source); }

  /*!
\brief Sets all bits after m to zero.
  */
  BitSet& truncate(size_t m) { Base::truncate(m); return *this; }

}; // |template<size_t n> class BitSet|

} // |namespace bitset|

} // |namespace atlas|

#include "bitset_def.h"

#endif
