/*
  This is test.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
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
