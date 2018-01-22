/*
  This is helpmode.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#include "helpmode.h"

#include <cstdio>   // not obviously used, but appears helpful for Windows
#include <iostream>

#include "commands.h"
#include "emptymode.h"
#include "io.h"
#include "mainmode.h"
#include "realmode.h"
#include "blockmode.h"
#include "reprmode.h"
#include "test.h"

/****************************************************************************

  This file contains the definition of the help mode.

*****************************************************************************/

namespace atlas {

namespace commands {

namespace {

  void help_entry();
  void help_exit();

  // help commands

  void questionMark_h();

  // command tags for the help facility

  const char* intro_tag =
    "(in help mode only) prints a message for first time users";
  const char* questionMark_tag =
    "(in help mode only) lists the available commands";

/****************************************************************************

        Chapter I --- The ThisMode class.

  Only one instance of this class will be constructed, on the first
  call of the emptyHelp function.

*****************************************************************************/

void help_entry()

{
  std::cerr << "enter \"?\" for a list of available commands; "
	    << "\"q\" to exit help mode"
	    << std::endl;
  std::cerr << "entering a command name will print an explanatory message"
	    << std::endl;
}

void help_exit()

{}

/*****************************************************************************

        Chapter II --- Functions for the help commands

  This section contains the definitions of the help functions associated to
  the various commands defined in this mode.

******************************************************************************/

void questionMark_h()

{
  std::cerr << std::endl;
  for (TagDict::const_iterator it = tagDict.begin(); it != tagDict.end(); ++it)
    std::cerr << "  - " << it->first << " : " <<  it->second << std::endl;

  std::cerr << std::endl;

  return;
}

} // |namespace|


/*****************************************************************************

        Chapter III --- Functions declared in helpmode.h

******************************************************************************/

/*
  Construct and return the help node.
  For static initialisation, this should be called from *this module* only,
  since it uses |tagDict| that might otherwise be uninitialised
*/
CommandNode helpNode()
{
  CommandNode result ("help: ",help_entry,help_exit);

  result.nohelp_add("?",questionMark_h);
  tagDict.insert(std::make_pair("?",questionMark_tag));

  result.nohelp_add("intro",intro_h);
  tagDict.insert(std::make_pair("intro",intro_tag));

  return result;
}

void intro_h()
{
  io::printFile(std::cerr,"intro_mess",io::MESSAGE_DIR);
}


} // |namespace helpmode|

} // |namespace atlas|
