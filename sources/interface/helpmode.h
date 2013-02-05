/*
  This is helpmode.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#ifndef HELPMODE_H  /* guard against multiple inclusions */
#define HELPMODE_H

#include "commands_fwd.h"
#include <iostream>

namespace atlas {

namespace commands {

/******** function declarations ********************************************/

  commands::CommandNode helpNode();
  extern commands::CommandTree help_mode; // defined in main.cpp

  void insertTag(TagDict&, const char*, const char*);
  void printTags(std::ostream&, const TagDict&);

  void intro_h(); // this is used in emptymode as well
  void nohelp_h(); // this may be used in test.cpp and other places

}

}

#endif
