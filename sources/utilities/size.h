/*!
\file
  This is size.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef SIZE_H  /* guard against multiple inclusions */
#define SIZE_H

#include <cstring>

#include "constants.h"

/******** type declarations **************************************************/

namespace atlas {

namespace size {

  template<typename C> class SizeType;
  // this may have to be modified if RANK_MAX increases
  typedef signed char BaseType;
  typedef unsigned char UnsignedBaseType;
  typedef SizeType<BaseType> Size;

}

/******** constant declarations **********************************************/

namespace size {

  // useless value
  /*!
  \brief Ordinal (position on the list of primes) of the largest prime
  factor of a Weyl group of rank at most n.

  The ordinal is recorded in PrimesMax<n>::value, and used to make the
  SizeType class for storing Weyl group orders.  This class should
  be instantiated only for n=RANK_MAX, which for other reasons should
  be a power of 2.
  */
  template<unsigned long n> class PrimesMax {
  public:
    static const unsigned long value = 0ul;
  };


  // predefined values for likely instances of RANK_MAX
  /*!
  \brief Position on the list of primes of the
  largest possible prime factor of a Weyl group of rank at most 8.

  This is the fourth prime 7.
  */
  template<> class PrimesMax<8> {
  public:
    static const unsigned long value = 4ul;
  };

  /*!
  \brief Position on the list of primes of the
  largest possible prime factor of a Weyl group of rank at most 16.

  This is the seventh prime 17 (appearing only in the Weyl group S_17
  of type A16).
  */
  template<> class PrimesMax<16> {
  public:
    static const unsigned long value = 7ul;
  };

  /*!
  \brief Position on the list of primes of the
  largest possible prime factor of a Weyl group of rank at most 32.

  This is the eleventh prime 31.
  */

  template<> class PrimesMax<32> {
  public:
    static const unsigned long value = 11ul;
  };

  /*!
  \brief Position on the list of primes of the
  largest possible prime factor of a Weyl group of rank at most 64.

  This is the eighteenth prime 61.
  */

  template<> class PrimesMax<64> {
  public:
    static const unsigned long value = 18ul;
  };

  /*!
  \brief Position on the list of primes of the
  largest possible prime factor of a Weyl group of rank at most RANK_MAX.

  With RANK_MAX=16, this is 7 (for the seventh prime 17).
  */

  const size_t PRIMES_MAX = PrimesMax<constants::RANK_MAX>::value;
}

/******** function declarations **********************************************/

namespace size {

  template<typename C> void factorial (SizeType<C>&, unsigned long);
  unsigned long prime(size_t);

}

/******** type definitions ***************************************************/

namespace size {
  /*!
  \brief Stores a positive integer as product of prime powers, using
  the first PRIMES_MAX primes.

  The exponent of the jth prime is d_data[j].  The reason for using
  this is that the software must occasionally deal with integers too
  large to be unsigned long; mostly these appear as cardinalities of
  Weyl groups.  The constant PRIMES_MAX is chosen so that the
  cardinality of any Weyl group of rank at most RANK_MAX can be
  represented as a SizeType.  (Note from DV: when RANK_MAX is 16, it
  appears to me that the largest possible prime factor of a Weyl group
  order is the sixth prime 13, so one should for elegance choose
  PRIMES_MAX=6.  The constant used above in that case is PRIMES_MAX=7.
  But taking PRIMES_MAX too large is harmless.)
  */
template<typename C> class SizeType {

 private:
  C d_data[PRIMES_MAX];

 public:
// constructors and destructors
  SizeType() {
    memset(d_data,0,PRIMES_MAX);
  }

  explicit SizeType(unsigned long);

  ~SizeType()
    {}

// copy and assignment
  SizeType(const SizeType& a) {
    memcpy(d_data,a.d_data,PRIMES_MAX);
  }

  SizeType& operator=(const SizeType& a) {
    memcpy(d_data,a.d_data,PRIMES_MAX); return *this;
  }

// accessors
  C operator[] (size_t j) const {
    return d_data[j];
  }

  bool operator== (const SizeType& c) const {
    return !memcmp(d_data,c.d_data,PRIMES_MAX);
  }

  bool operator!= (const SizeType& c) const {
    return memcmp(d_data,c.d_data,PRIMES_MAX);
  }

  bool hasOverflow() const;

  bool hasOverflow(size_t) const;

  unsigned long piece(size_t) const;

  unsigned long toUlong() const;

// manipulators
  C& operator[] (size_t j) {
    return d_data[j];
  }

  SizeType& operator*= (const SizeType&);

  SizeType& operator*= (unsigned long);

  SizeType& operator/= (const SizeType&);

  void reset() {
    std::memset(d_data,0,PRIMES_MAX);
  }

  void twoShift(C n) { // multiplication by 2^n; prime #0 is 2
    d_data[0] += n;
  }
};

}

}

#include "size_def.h"

#endif
