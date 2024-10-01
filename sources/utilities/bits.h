/*
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// This is bits.h

#ifndef BITS_H  /* guard against multiple inclusions */
#define BITS_H

/******** function declarations **********************************************/

namespace atlas {

namespace bits {

  // count number of set bits in the word
  unsigned int bitCount(unsigned long long int);

  // position of first set bit in the word, or |longBits| if none
  unsigned int firstBit(unsigned long long int);

  // position PLUS ONE of last set bit in the word, or |0| if none
  unsigned int lastBit(unsigned long long int);

} // |namespace bits|

} // |namespace atlas|

#endif
