/*
  This is constants.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

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
  const unsigned long ulongMax = std::numeric_limits<unsigned long>::max();

  // I was surprised by the fact that the digits for char returns 7
  const unsigned long charBits = std::numeric_limits<unsigned char>::digits;
  const unsigned long longBits = std::numeric_limits<unsigned long>::digits;
  const unsigned long sizeBits = std::numeric_limits<size_t>::digits;

  const unsigned long hiBit = 1UL << longBits - 1UL;

  const unsigned long firstChar = (1UL << charBits) - 1UL;

  extern unsigned long bitMask[longBits];
  extern unsigned long twoBitMask[longBits >> 1];
  extern size_t firstbit[1 << charBits];
  extern size_t lastbit[1 << charBits];
  extern unsigned long leqMask[longBits];
  extern unsigned long lMask[longBits];

  // check PRIMES_MAX in size.h if you change this!
  const size_t RANK_MAX = 16UL; // should not exceed 128 or the data type for
                                // the Weyl groups will overflow
                                // should be a power of two

  // constants used to pick a bit-address apart
  // the first one serves as flags for the address within a word.
  // the second one is its logical complement
  // It is assumed that the number of digits in an unsigned long 
  // is a power of two.
  
  const unsigned long posBits = longBits - 1UL;
  const unsigned long baseBits = ~posBits;

  // this one tells by how much we have to shift to get the base
  // (it is the base 2 logarithm of longBits). 
  // The difficulty is to compute this at compile time! We use a neat
  // (and well-known) template trick for that. The result is correct
  // only when n is a power of two.

  template<unsigned long n> class BaseShift {
  public:
    static const unsigned long value = 1UL + BaseShift< n >> 1UL>::value;
  };

  template<> class BaseShift<1UL> {
  public:
    static const unsigned long value = 0UL;
  };

  const unsigned long baseShift = BaseShift<longBits>::value;

  template<size_t n> class Bool {
  public:
    static const size_t value = 1UL;
  };

  template<> class Bool<(size_t)0UL> {
  public:
    static const size_t value = 0UL;
  };
}

/******** function definitions ***********************************************/

namespace constants {

  void initConstants();

}

}

#endif
