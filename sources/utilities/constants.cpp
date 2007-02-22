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
unsigned long twoBitMask[longBits/2];
unsigned char firstbit[1 << charBits];
unsigned char lastbit[1 << charBits];
unsigned long leqMask[longBits];
unsigned long lMask[longBits+1];

const unsigned long* Primes = 0;

}

/******** function definitions ***********************************************/

namespace constants {

void initConstants()

/*
  This function initializes the following constants :

    - bitMask : bitMask[j] flags bit j            : bitMask[j]==1ul<<j
    - leqMask : leqMask[j] flags bits i with i<=j : leqMask[j]==(1ul<<j+1)-1
    - lMask : lMask[j]     flags bits i with i<j  : lMask[j]==(1ul<<j)-1

    for lMask the index j==longBits is useful and used; Fokko forgot this!

    - firstbit : firstbit[j] gives position of first set bit in j (0<j<256),
		 or 8 (charBits) for j==0
    - lastbit : lastbit[j] gives position+1 of last set bit in j (0<j<256),
		 or 0 for j==0;

*/

{
  lMask[0] = 0; // probably unused, but not written in the following loop

  for (unsigned long j = 0; j < longBits; ++j) {
    bitMask[j] = 1ul << j;                    // bit j set
    leqMask[j] = (bitMask[j]-1) | bitMask[j]; // bits 0..j set
    lMask[j+1] = leqMask[j];                  // lMask[j+1] has bits 0..j set
  }

  for (unsigned long j = 0; j < longBits/2; ++j)
    twoBitMask[j] = 3ul << 2*j;  // bits 2*j,2*j+1 set

  firstbit[0] = charBits; // out-of-bounds indication
  firstbit[1] = 0;

  // find least significant set bit by a simple recurrence
  for (unsigned int j = 1; j < (1 << (charBits-1)); ++j) {
    firstbit[2*j]   = firstbit[j]+1; // for even numbers use recursion
    firstbit[2*j+1] = 0;             // for odd numbers, it is 0
  }

  lastbit[0] = 0; // out-of-bounds indication
  lastbit[1] = 1; // this points to bit 0, because of position+1 convention

  // find most significant set bit by a simple recurrence
  for (unsigned int j = 1; j < (1 << charBits-1); ++j) {
    lastbit[2*j]   = lastbit[j]+1; // for even numbers use recursion
    lastbit[2*j+1] = lastbit[j]+1; // for odd numbers use recursion as well
  }
}

}

}
