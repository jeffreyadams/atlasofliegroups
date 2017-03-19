/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// This is bits.h

#ifndef BITS_H  /* guard against multiple inclusions */
#define BITS_H

/******** function declarations **********************************************/

namespace atlas {

namespace bits {

  unsigned int bitCount(unsigned long long int);

  unsigned int firstBit(unsigned long long int);

  unsigned int lastBit(unsigned long long int);

} // |namespace bits|

} // |namespace atlas|

#endif
