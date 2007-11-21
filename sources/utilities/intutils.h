/*!
\file
  This is intutils.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  For license information see the LICENSE file
*/

#ifndef INTUTILS_H  /* guard against multiple inclusions */
#define INTUTILS_H

/******** function declarations *********************************************/

namespace atlas {

namespace intutils {

  template<typename I> I abs(I);
  template<typename I> I divide(I, I);
  template<typename I> I factorial(I);

}

}

#include "intutils_def.h"

#endif
