/*
  This is constants.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2018 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// this module defines various constants, mostly related to bit packing
// it is mostly geared towards the |BitMap| class for large bit-mapped tables

#ifndef CONSTANTS_H  /* guard against multiple inclusions */
#define CONSTANTS_H

#include <cstddef>
#include <limits>

/******** constants **********************************************************/

namespace atlas {

  namespace aux {
    template<unsigned long n> struct BaseShift // value gives 2log of n
    { static constexpr unsigned long value = 1 + BaseShift<n/2>::value; };
    template<> struct BaseShift<1ul>
    { static constexpr unsigned long value = 0; };
  }

struct constants
{
  static constexpr unsigned charBits = // typically (and at least) 8
    std::numeric_limits<unsigned char>::digits;
  static constexpr unsigned char ucharMax = ~0; // |0xFF|, all bits in byte set
  static constexpr unsigned longBits = // typically 64
    std::numeric_limits<unsigned long>::digits;
  static constexpr unsigned long hiBit = // typically |0x80000000ull|
    1ul << (longBits - 1);
  static constexpr unsigned long firstCharMask = // typically |0xFFull|
    (1ul << charBits) - 1;
  static constexpr unsigned posBits = // (often widened), typically |0x3F==63|
    longBits - 1; // part of bit address for position within ull word
  static constexpr unsigned long baseBits = // typically |0xFFFFFFFFFFFFFFC0|
    ~static_cast<unsigned long>(posBits); // mask for word part in a bit address
  static constexpr unsigned long baseShift = // typically 6
    aux::BaseShift<longBits>::value;

  static unsigned char firstbit[1ul << charBits];
  static unsigned char lastbit[1ul << charBits]; // last bit PLUS ONE
  static unsigned long eq_mask[longBits]; // bit set at position |i| only
  static unsigned long lt_mask[longBits+1]; // bits at positions |<i| set
  static const unsigned long* leq_mask; // effectively bits at positions |<=i|

/*
  Implementation bound on rank of reductive group, used to allow working
  with fixed-size bit-arrays or char-arrays. Should (probably) be a power of 2,
  and should not exceed 128 (assumption in Weyl group implementation).
*/
  static const size_t RANK_MAX = 32;

private: // trickery to ensure initialization is performed exactly once
  static constants init();
  static const constants dummy; // a |static| variable can be of own class!
}; // |struct constants|

} // |namespace atlas|

#endif
