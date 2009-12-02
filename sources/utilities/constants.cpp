/*
  This is constants.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2008,2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "constants.h"

#include <vector>

#include "error.h"

/*****************************************************************************

  This file provides some constants, particularly for bit-masking.

******************************************************************************/

/******** constant definitions **********************************************/

namespace atlas {

  // definition of static arrays
  unsigned long constants::bitMask[longBits];
  unsigned char constants::firstbit[1 << charBits];
  unsigned char constants::lastbit[1 << charBits];
  unsigned long constants::leqMask[longBits];
  unsigned long constants::lMask[longBits+1];


/******** static member function *******************************************/


/*
  Since we cannot initialize the above static member arrays by initialzer
  lists (we do not even know how long they should be), we provide a static
  member initialisation function, and ensure that it is called (exaclty) once

  This function initializes the following constants :

    - bitMask : bitMask[j] flags bit j            : bitMask[j]==1ul<<j
    - leqMask : leqMask[j] flags bits i with i<=j : leqMask[j]==(1ul<<j+1)-1
    - lMask : lMask[j]     flags bits i with i<j  : lMask[j]==(1ul<<j)-1

    for lMask the index j==longBits is useful and used; Fokko forgot this!

    - firstbit : |firstbit[j]| gives greatest $k$ such that $2^k$ divides $j$
                 (number of trailing bits 0); it is the position of lowest
		 set bit for $j>0$, and $8$ (|charBits|) for $j=0$
    - lastbit : |lastbit[j]| gives least $k$ with $j<2^k$; it is ONE MORE than
                position of highest set bit in $j$ for $j>0$, and $0$ for $j=0$

*/

constants constants::init()
{
  lMask[0] = 0; // probably unused, but not written in the following loop

  for (unsigned long j = 0; j < longBits; ++j)
  {
    bitMask[j] = 1ul << j;                    // bit j set
    leqMask[j] = (bitMask[j]-1) | bitMask[j]; // bits 0..j set
    lMask[j+1] = leqMask[j];                  // lMask[j+1] has bits 0..j set
  }

  firstbit[0] = charBits; // out-of-bounds indication
  firstbit[1] = 0;

  // find least significant set bit by a simple recurrence
  for (unsigned int j = 1; j < (1 << (charBits-1)); ++j)
  {
    firstbit[2*j]   = firstbit[j]+1; // for even numbers use recursion
    firstbit[2*j+1] = 0;             // for odd numbers, it is 0
  }

  lastbit[0] = 0; // out-of-bounds indication
  lastbit[1] = 1; // this points to bit 0, because of position+1 convention

  // find most significant set bit by a simple recurrence
  for (unsigned int j = 1; j < (1 << (charBits-1)); ++j)
  {
    lastbit[2*j]   = lastbit[j]+1; // for even numbers use recursion
    lastbit[2*j+1] = lastbit[j]+1; // for odd numbers use recursion as well
  }

  return constants(); // return an (empty) instance of our own class
}

/* Code that makes sure that |init| is called whenever this module is linked
   into a program
*/

  const constants constants::dummy = constants::init(); // ensure call |init|

} // namespace atlas
