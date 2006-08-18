/*!
\file
  This is realform.h
\brief Type definitions for the namespace realform.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef REALFORM_H  /* guard against multiple inclusions */
#define REALFORM_H

#include <cstddef>
#include <vector>

/******** forward type definitions *******************************************/

namespace atlas {

namespace realform {
/*!
  A variable of this type indicates the position (like 0, 1, 2,...) of a 
  real form on a list of real forms (like the one maintained by
  the class ComplexReductiveGroup in realFormLabels).
*/
  typedef size_t RealForm;
  typedef std::vector<RealForm> RealFormList;
}

}

#endif
