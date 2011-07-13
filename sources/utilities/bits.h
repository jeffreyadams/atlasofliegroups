/*!
\file
  This is bits.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef BITS_H  /* guard against multiple inclusions */
#define BITS_H

#include <vector>

/******** function declarations **********************************************/

namespace atlas {

namespace bits {

  unsigned bitCount(unsigned long);

  size_t firstBit(unsigned long);

  size_t lastBit(unsigned long);

}

}

#endif
