/*!
\file
  This is set.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations 

  For license information see the LICENSE file
*/

#ifndef SET_H  /* guard against multiple inclusions */
#define SET_H

#include <cstddef>
#include <vector>

/******** forward type declarations *****************************************/

namespace atlas {

namespace set {

  typedef size_t Elt;
  typedef std::vector<Elt> EltList;

}

}

#endif
