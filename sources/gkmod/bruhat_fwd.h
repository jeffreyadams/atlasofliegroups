/*
  This is bruhat_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef BRUHAT_FWD_H  /* guard against multiple inclusions */
#define BRUHAT_FWD_H

#include <cstddef>
#include <vector>

#include "bitmap_fwd.h"

namespace atlas {

/******** type declarations *************************************************/

namespace bruhat {

class BruhatOrder;

typedef size_t BruhatElt;
typedef std::vector<BruhatElt> BruhatEltList;

typedef bitmap::BitMap BruhatRow;

}

}

#endif
