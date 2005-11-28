/*
  This is kltest.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef KLTEST_H  /* guard against multiple inclusions */
#define KLTEST_H

#include "kgb_fwd.h"
#include "kl_fwd.h"

#include "setutils.h"

namespace atlas {

/******** type declarations *************************************************/

namespace kltest {

}

/******** function declarations **********************************************/

namespace kltest {

  bool checkBasePoint(const kgb::KGB&);

  void dualityPermutation(setutils::Permutation&, const kl::KLContext&);
}

/******** type definitions **************************************************/

namespace kltest {

}

}

#endif
