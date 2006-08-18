/*!
\file
  This is constants.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef CONSTANTS_H  /* guard against multiple inclusions */
#define CONSTANTS_H

#include <cstddef>
#include <limits>

/******** constants **********************************************************/

namespace atlas {

namespace constants {

  const unsigned char ucharMax = std::numeric_limits<unsigned char>::max();
  /*!
  2**32 - 1 on Mac.
   */
  const unsigned long ulongMax = std::numeric_limits<unsigned long>::max();
  /*!
  8 on Mac.  Fokko says 7 on Debian Linux.
   */
  // I was surprised by the fact that the digits for char returns 7
  const unsigned long charBits = std::numeric_limits<unsigned char>::digits;
  /*!
  32 on Mac.
   */
  const unsigned long longBits = std::numeric_limits<unsigned long>::digits;
  /*!
  32 on Mac.
  */
  const unsigned long sizeBits = std::numeric_limits<size_t>::digits;

  const unsigned long hiBit = 1ul << longBits - 1ul;

  const unsigned long firstChar = (1ul << charBits) - 1ul;
  /*!
  bitMask[j] is 2**j: a single word with a 1 in bit j and 0's elsewhere.
  */
  extern unsigned long bitMask[longBits];
  /*!  
  twoBitMask[j] is 3*(4**j): a single word with a 1 in bits 2j
  and 2j+1, and 0's elsewhere.
  */
  extern unsigned long twoBitMask[longBits >> 1ul];
  extern size_t firstbit[1ul << charBits];
  extern size_t lastbit[1ul << charBits];
  extern unsigned long leqMask[longBits];
  extern unsigned long lMask[longBits];

  
  // check PRIMES_MAX in size.h if you change this!
  /*!
  RANK_MAX is the largest allowed rank for a reductive group.  It
  should not exceed 128, or the data type for Weyl groups will
  overflow (in a way requiring fundamental rewriting of bitSet and
  bitVector).  It should always be a power of two.  The constant
  PRIMES_MAX in size.h has to be set in accordance with RANK_MAX, so
  that Weyl group sizes can be properly stored.
  */
  const size_t RANK_MAX = 16ul; // should not exceed 128 or the data type for
                                // the Weyl groups will overflow
                                // should be a power of two

  // constants used to pick a bit-address apart
  // the first one serves as flags for the address within a word.
  // the second one is its logical complement
  // It is assumed that the number of digits in an unsigned long 
  // is a power of two.
  /*!
  Constant used to pick a bit-address apart: serves as flags for the
  address within a word.  It is assumed that the number of bits in an
  unsigned long is a power of two.
  */
  const unsigned long posBits = longBits - 1ul;
  /*!
  Constant used to pick a bit-address apart: this is the logical
  complement of posBits.  It is assumed that the number of bits in an
  unsigned long is a power of two.
  */
  const unsigned long baseBits = ~posBits;
  
  /*!
  \brief Computes (in the constant BaseShift<n>::value)
  the base two logarithm of n.  
  
  This is used to to compute the base 2 logarithm of longBits at
  compile time.
  */
  template<unsigned long n> class BaseShift {
  public:
    static const unsigned long value = 1ul + BaseShift<n/2ul>::value;
  };

  template<> class BaseShift<1ul> {
  public:
    static const unsigned long value = 0ul;
  };

  /*!
  \brief  This one tells by how much we have to shift a binary number N to get
  the base (number of machine words needed for a bitmap of size N).  

  It is the base 2 logarithm of longBits. The difficulty is to compute
  this at compile time! We use a neat (and well-known) template trick
  for that. The result is correct only when n is a power of two.
  (Should be 5 on 32 bit machine.)
  */
  const unsigned long baseShift = BaseShift<longBits>::value;

  /*!
  \brief Tests whether n is non-zero.  

  The constant Bool<n>::value is 0 if n is zero and 1 if n is not
  zero.
  */
  template<size_t n> class Bool {
  public:
    static const size_t value = 1ul;
  };

  template<> class Bool<(size_t)0ul> {
  public:
    static const size_t value = 0ul;
  };
}

/******** function definitions ***********************************************/

namespace constants {

  void initConstants();

}

}

#endif
