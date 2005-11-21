/*
  This is cartanclass_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef CARTANCLASS_FWD_H  /* guard against multiple inclusions */
#define CARTANCLASS_FWD_H

#include <utility>

/******** forward type definitions *******************************************/

namespace atlas {

namespace cartanclass {

  class CartanClass;
  class Fiber;

  typedef unsigned long RealFormRep;
  typedef std::pair<unsigned long,unsigned long> StrongRealFormRep;
}

}

#endif
