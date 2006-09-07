/*
  This is constants.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "constants.h"

#include <vector>

#include "error.h"

/*****************************************************************************

  This file provides some constants, particularly for bit-masking.

******************************************************************************/

/******** constant definitions **********************************************/

namespace atlas {

namespace constants {

unsigned long bitMask[longBits];
  unsigned long twoBitMask[longBits >> 1];
size_t firstbit[1 << charBits];
size_t lastbit[1 << charBits];
unsigned long leqMask[longBits];
unsigned long lMask[longBits];

const unsigned long* Primes = 0;

}

/******** function definitions ***********************************************/

namespace constants {

void initConstants()

/*
  This function initializes the following constants :

    - bitMask : bitMask[j] flags the j-th bit;
    - firstbit : firstbit[j] gives position of first set bit in j;
*/

{
  bitMask[0] = 1UL;
  twoBitMask[0] = 3UL;
  leqMask[0] = 1UL;
  lMask[0] = 0;

  for (unsigned long j = 1; j < longBits; ++j) {
    bitMask[j] = bitMask[j-1] << 1;
    leqMask[j] = (leqMask[j-1] << 1) + 1;
    lMask[j] = leqMask[j-1];
  }

  for (unsigned long j = 1; j < (longBits >> 1); ++j) {
    twoBitMask[j] = twoBitMask[j-1] << 2;
  }

  firstbit[0] = charBits;

  for (unsigned long j = 1; j < (1 << charBits-1); ++j)
    firstbit[2*j] = firstbit[j]+1;

  lastbit[0] = 0;

  for (unsigned long j = 1; j < (1 << charBits); ++j)
    lastbit[j] = lastbit[j>>1]+1;

  return;
}

}

}
