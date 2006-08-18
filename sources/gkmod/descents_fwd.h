/*!
\file
\brief Class declaration and type definitions for class DescentStatus.
*/
/*
  This is descents_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef DESCENTS_FWD_H  /* guard against multiple inclusions */
#define DESCENTS_FWD_H

#include <vector>

namespace atlas {

/******** type declarations *************************************************/

namespace descents {

class DescentStatus;
typedef std::vector<DescentStatus> DescentStatusList;

}

}

#endif
