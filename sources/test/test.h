/*
  This is test.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2013 Marc van Leeuwen
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

/* The following function is called to complete the various |CommandNode|
   objects with actions defined in test.cpp, prior to constructing the
   corresponding |CommandTree|. This addition cannot be done on the completed
   |CommandTree|, as the commands would the fail to inherit to descendant
   modes. Therefore the addition cannot be done in one sweep, and the pieces
   for different |CommandNode| values are done in separate calls. Rather than
   use a single function with a switch, we might as well use a template
   function whose instances are distinguished by a dummy type. The node to
   alter is passed as argument.
 */
template<typename T> void addTestCommands(commands::CommandNode&);

}

}

#endif
