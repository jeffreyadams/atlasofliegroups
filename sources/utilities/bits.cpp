/*
  This is bits.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006,2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "bits.h"

#include "constants.h"


/*****************************************************************************

        Chapter I -- Functions declared in bits.h

******************************************************************************/

namespace atlas {

namespace bits {

// Return the sum of the bits (i.e., the number of set bits) in |x|.

unsigned int bitCount(unsigned long long int x)
/*
  This used to be the following code
  {
    unsigned int count = 0;
    while (x!=0) { x &= x-1 ; ++count; }
    return count;
  }

  The straight-line code below, taken from Knuth (pre-fascicle~1a ``Bitwise
  tricks and techniques'' to volume~4 of ``The Art of Computer Programming'',
  p.~11), is attributed to D.~B. Gillies and J.~C.~P. Miller. It is faster on
  the average (assuming non-sparse bitsets), without depending on the number
  of bits in |unsigned long long| (except that it is a multiple of 8 and
  less than 256)
*/
{ static const unsigned long long int b0= ~(~0ull/0x3);
     // |0xAAAA...|; flags odd bit positions
  static const unsigned long long int b1= ~(~0ull/0x5);
     // |0xCCCC...|; flags positions $\cong2,3 \pmod4$
  static const unsigned long int b2= ~(~0ull/0x11);
     // |0xF0F0...|; flags positions $\cong4$--$7\pmod8$
  static const unsigned long int ones= ~0ull/0xFF;
     // |0x0101...|; flags the low bit of each octet
  static const unsigned int high_byte_shift =
    8*(sizeof(unsigned long long int)-1);

  x-=(x&b0)>>1;          // replace pairs of bits $10\to01$ and $11\to10$
  x=(x&~b1)+((x&b1)>>2);
   // sideways add 2 groups of pairs of bits to 4-tuples of bits
  x += x>>4;
   // the sums of octets (bytes) are now in lower 4-tuples of those octets
  return (x&~b2)*ones >> high_byte_shift;
   // add lower 4-tuples of bytes in high octet, and extract
}

/*
  Return the position of the first set bit in |f|.
  Return |constants::longBits| if there is no such bit.
*/
unsigned int firstBit(unsigned long long int f)
{
  if (f == 0)
    return constants::longBits;

  unsigned int fb = 0;

  for (; (f & constants::firstCharMask) == 0; f >>= constants::charBits)
    fb += constants::charBits;

  return fb + constants::firstbit[f & constants::firstCharMask];
}


/*
  Return the position of the last (most significant) set bit in |f|, PLUS ONE.
  Return 0 if there is no such bit.
*/
unsigned int lastBit(unsigned long long int f)
{
  if (f == 0)
    return 0;

  unsigned int lb; // number of skipped low-order bits
  for (lb=0; (f & ~constants::firstCharMask)!=0; f >>= constants::charBits)
    lb += constants::charBits;

  return lb + constants::lastbit[f];
}

} // |namespace bits|

} // |namespace atlas|
