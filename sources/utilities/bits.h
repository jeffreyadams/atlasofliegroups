/*
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// This is bits.h

#ifndef BITS_H  /* guard against multiple inclusions */
#define BITS_H

#include <cstring>
#include "setutils.h"

/******** function declarations **********************************************/

namespace atlas {

namespace bits {

  // count number of set bits in the word
  unsigned int bitCount(unsigned long);

  inline void copy(void* dest, const void* source, size_t d) {
    std::memcpy(dest,source,d);
  }

  size_t firstBit(unsigned long);

  size_t lastBit(unsigned long);

  void permute(unsigned long&, const setutils::Permutation&);

}

}

#endif
