/*
  This is bits.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef BITS_H  /* guard against multiple inclusions */
#define BITS_H

#include <vector>

#include "constants.h"
#include "setutils.h"

/******** function declarations **********************************************/

namespace atlas {

namespace bits {

  unsigned bitCount(unsigned long);

  inline void copy(void* dest, const void* source, size_t d) {
    memcpy(dest,source,d);
  }

  size_t firstBit(unsigned long);

  size_t lastBit(unsigned long);

  void permute(unsigned long&, const setutils::Permutation&);

}

}

#endif
