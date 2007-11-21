/*!
\file
  This is bits.cpp
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "bits.h"

#include "constants.h"


/*****************************************************************************

        Chapter I -- Functions declared in bits.h

******************************************************************************/

namespace atlas {

namespace bits {

/*!
  Synopsis: returns the sum of the bits (i.e., the number of set bits) in x.
*/
unsigned int bitCount(unsigned long x)
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
  of bits in |unsigned long| (except that it is a multiple of 8 less than 256)
*/
{ static const unsigned long int b0= ~(~0ul/3);
     // |0xAAAA...|; flags odd bit positions
  static const unsigned long int b1= ~(~0ul/5);
     // |0xCCCC...|; flags positions $\cong2,3 \pmod4$
  static const unsigned long int b2= ~(~0ul/17);
     // |0xF0F0...|; flags positions $\cong4$--$7\pmod8$
  static const unsigned long int ones= ~0ul/255;
     // |0x0101...|; flags the low bit of each octet
  static const unsigned int high_byte_shift=8*(sizeof(unsigned long int)-1);

  x-=(x&b0)>>1;          // replace pairs of bits $10\to01$ and $11\to10$
  x=(x&~b1)+((x&b1)>>2);
   // sideways add 2 groups of pairs of bits to 4-tuples of bits
  x += x>>4;
   // the sums of octets (bytes) are now in lower 4-tuples of those octets
  return (x&~b2)*ones >> high_byte_shift;
   // add lower 4-tuples of bytes in high octet, and extract
}

size_t firstBit(unsigned long f)

/*!
  Synopsis: returns the position of the first set bit in f.

  Returns longBits if there is no such bit.
*/

{
  using namespace constants;

  if (f == 0)
    return longBits;

  size_t fb = 0;

  for (; (f & firstCharMask) == 0; f >>= charBits)
    fb += charBits;

  return fb + firstbit[f & firstCharMask];
}

size_t lastBit(unsigned long f)

/*!
  Synopsis: returns the position of the last (most significant)
  set bit in f, PLUS ONE.  Returns 0 if there is no such bit.
*/

{
  using namespace constants;

  if (f == 0)
    return 0;

  unsigned lb = 0;

  for (; f & ~firstCharMask; f >>= charBits)
    lb += charBits;

  return lb + lastbit[f];
}

void permute(unsigned long& f, const setutils::Permutation& a)

/*!
  Synopsis: permutes the bits of f according to a: if the bits of f are
  interpreted with respect to a sequence of objects x_0,...,x_{n-1}, then the
  bits of the result are interpreted with respect to the permuted sequence
  x_{a[0]},...,x_{a[n-1]}

  Precisely, we transform f = f_old into f_new, where f_new[i] = f_old[a[i]].

  Note that the permutation of bits is not interpreted in the same way as the
  permutation of vector entries in setutils::permute, given in setutils_def.h
*/

{
  using namespace constants;

  unsigned long f_new = 0;

  for (size_t j = 0; j < a.size(); ++j)
    if (f & bitMask[a[j]])
      f_new |= bitMask[j];

  f = f_new;

  return;
}

}

}
