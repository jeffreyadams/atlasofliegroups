/*
  This is constants.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2008,2009 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

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
  unsigned long constants::eq_mask[longBits];
  unsigned char constants::firstbit[1 << charBits];
  unsigned char constants::lastbit[1 << charBits];
  unsigned long constants::lt_mask[longBits+1];
  const unsigned long* constants::leq_mask;


/******** static member function *******************************************/


/*
  Since we cannot initialize the above static member arrays by initializer
  lists (we do not even know how long they should be), we provide a static
  member initialisation function, and ensure that it is called (exactly) once

  This function initializes the following constants :

    - eq_mask : eq_mask[j] flags bit j            : eq_mask[j]==1ul<<j
    - lt_mask : lt_mask[j] flags bits i with i<j  : lt_mask[j]==(1ul<<j)-1
    - leqMask : leqMask[j] flags bits i with i<=j : leqMask[j]==(1ul<<j+1)-1

    for |lt_mask|, the index j==longBits is useful and used; Fokko forgot this!

    - firstbit : |firstbit[j]| gives greatest $k$ such that $2^k$ divides $j$
                 (number of trailing bits 0); it is the position of lowest
		 set bit for $j>0$, and $8$ (|charBits|) for $j=0$
    - lastbit : |lastbit[j]| gives least $k$ with $j<2^k$; it is ONE MORE than
                position of highest set bit in $j$ for $j>0$, and $0$ for $j=0$

*/

constants constants::init()
{
  lt_mask[0] = 0;
  auto* leq_base = &lt_mask[1];
  leq_mask = leq_base; // the same, but as pointer to |const|
  for (unsigned long j = 0; j < longBits; ++j)
  {
    eq_mask[j] = 1ul << j;              // bit j set
    leq_base[j] = lt_mask[j] | eq_mask[j]; // bits 0..j (inclusive) set
  }

  firstbit[0] = charBits; // out-of-bounds indication
  firstbit[1] = 0;
  lastbit[0] = 0; // out-of-bounds indication
  lastbit[1] = 1; // this points to bit 0, because of position+1 convention


  // find least and most significant set bits by a simple recurrence
  for (unsigned int j = 1; j < (1 << (charBits-1)); ++j)
  {
    firstbit[2*j]   = firstbit[j]+1; // for even numbers use recursion
    firstbit[2*j+1] = 0;             // for odd numbers, it is 0
    lastbit[2*j] = lastbit[2*j+1] = lastbit[j]+1; // two same values, recursion
  }

  return constants(); // return an (empty) instance of our own class
} // |constants::init|

/* Code that makes sure that |init| is called whenever this module is linked
   into a program
*/

  const constants constants::dummy = constants::init(); // ensure call |init|

} // |namespace atlas|
