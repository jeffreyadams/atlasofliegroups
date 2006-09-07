/*
  This is commands_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#ifndef COMMANDS_FWD_H  /* guard against multiple inclusions */
#define COMMANDS_FWD_H

#include <map>

/******** forward type declarations *****************************************/

namespace atlas {

namespace commands {

  struct Command;
  class CommandMode;
  class StrCmp;

  typedef std::map<const char*,const char*, StrCmp> TagDict;

  struct EntryError{};

}

}

#endif
