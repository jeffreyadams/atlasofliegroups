/*!
\file
  This is setutils.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef SETUTILS_H  /* guard against multiple inclusions */
#define SETUTILS_H

#include <vector>

/******** type declarations *************************************************/

namespace atlas {

namespace setutils {

  typedef std::vector<unsigned long> Permutation;

}

/******** function declarations **********************************************/

namespace setutils {

  void compose(Permutation&, const Permutation&, unsigned long n = 0);

  void identity(Permutation&, unsigned long);

  void invert(Permutation&, const Permutation&);

  template<typename T> void permute(const Permutation&, std::vector<T>&);

}

}

/******** template definitions ***********************************************/

#include "setutils_def.h"

#endif
