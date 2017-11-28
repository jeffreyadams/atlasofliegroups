/*!
\file
  This is set.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#ifndef SET_H  /* guard against multiple inclusions */
#define SET_H

#include <cstddef>
#include <vector>

/******** forward type declarations *****************************************/

namespace atlas {

namespace set {

  typedef size_t SetElt;
  typedef std::vector<SetElt> SetEltList;

}

}

#endif
