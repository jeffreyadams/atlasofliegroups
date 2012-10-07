/*
  This is test.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations 

  For license information see the LICENSE file
*/

#ifndef TEST_H  /* guard against multiple inclusions */
#define TEST_H

#include <map>
#include <string>

#include "commands.h"

/******** function declarations ********************************************/

namespace atlas {

namespace test {

template<typename T> void addTestCommands(commands::CommandMode&, T);

template<typename T> 
  void addTestHelp(commands::CommandMode&, commands::TagDict&, T);

}

}

#endif
