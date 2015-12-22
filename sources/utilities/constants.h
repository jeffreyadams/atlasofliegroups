/*!
\file
  This is constants.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef CONSTANTS_H  /* guard against multiple inclusions */
#define CONSTANTS_H

// extra defs for windows compilation -spc
#ifdef WIN32
#define and &&
#define or ||
#define not !
#endif

#include <cstddef>
#include <limits>

/******** constants **********************************************************/

namespace atlas {

  namespace aux {
    template<unsigned long n> struct BaseShift // value gives 2log of n
    { static const unsigned long value = 1 + BaseShift<n/2>::value; };
    template<> struct BaseShift<1ul> { static const unsigned long value = 0; };
  }

struct constants
{

  static const unsigned charBits = std::numeric_limits<unsigned char>::digits;
  static const unsigned char ucharMax = ~0; // all bits in byte set
  static const unsigned longBits = std::numeric_limits<unsigned long>::digits;
  static const unsigned sizeBits = std::numeric_limits<size_t>::digits;
  static const unsigned long hiBit = 1ul << (longBits - 1);
  static const unsigned long firstCharMask = (1ul << charBits) - 1;
  static const unsigned long posBits = longBits - 1; // longBits is power of 2
  static const unsigned long baseBits = ~posBits; // mask for non bit-position
  static const unsigned long baseShift = aux::BaseShift<longBits>::value;

  static unsigned long bitMask[longBits];
  static unsigned char firstbit[1ul << charBits];
  static unsigned char lastbit[1ul << charBits]; // last bit PLUS ONE
  static unsigned long leqMask[longBits];
  static unsigned long lMask[longBits+1];

  /*! Implementation bound on rank of reductive group, used to allow working
      with fixed-size bit-arrays or char-arrays. Should (probably) be a power
      of 2, not exceeding 128 (assumption in Weyl group implementation). */
  static const size_t RANK_MAX = 32;

private: // trickery to ensure initialization is performed exactly once
  static constants init();
  static const constants dummy; // a |static| variable can be of own class!
}; // struct constants

} // namespace atlas

#endif
