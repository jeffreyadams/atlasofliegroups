/*!
\file
\brief Class declarations for BitVector.
*/
/*  
  This is bitvector_fwd.h
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef BITVECTOR_FWD_H  /* guard against multiple inclusions */
#define BITVECTOR_FWD_H

#include <cstddef>

/******** forward type declarations ******************************************/

namespace atlas {

namespace bitvector {

  template<size_t> class BitVector;
  template<size_t> class BitMatrix;
  template<size_t> class FirstBit;

}

}

#endif
