/*
  This is helpmode.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For copyright and license information see the LICENSE file
*/

#include "helpmode.h"

#include <iostream>

#include "commands.h"
#include "emptyhelp.h"
#include "emptymode.h"
#include "io.h"
#include "mainhelp.h"
#include "mainmode.h"
#include "realhelp.h"
#include "realmode.h"
#include "blockmode.h"
#include "special.h"
#include "test.h"

/****************************************************************************

  This file contains the definition of the help mode.

*****************************************************************************/

namespace atlas {

namespace {
  using namespace helpmode;

  void help_entry();
  void help_exit();

  commands::TagDict tagDict;

  // help commands

  void help_h();
  void qq_h();
  void questionMark_h();

  // command tags for the help facility

  const char* intro_tag =
    "(in help mode only) prints a message for first time users";
  const char* q_tag = "exits the current mode";
  const char* questionMark_tag =
    "(in help mode only) lists the available commands";
}

/****************************************************************************

        Chapter I --- The ThisMode class.

  Only one instance of this class will be constructed, on the first
  call of the emptyHelp function.

*****************************************************************************/

namespace {

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

}

/*****************************************************************************

        Chapter II --- Functions for the help commands

  This section contains the definitions of the help functions associated to
  the various commands defined in this mode.

******************************************************************************/

namespace {

void help_h()

{
  io::printFile(std::cerr,"help.help",io::MESSAGE_DIR);
  return;
}

void qq_h()

{
  io::printFile(std::cerr,"qq.help",io::MESSAGE_DIR);
  return;
}

void questionMark_h()

{
  std::cerr << std::endl;
  commands::printTags(std::cerr,tagDict);
  std::cerr << std::endl;

  return;
}

}

/*****************************************************************************

        Chapter III --- Functions declared in helpmode.h

******************************************************************************/

namespace helpmode {

/*
  Synopsis: returns the help mode.

  It is constructed on the first call.
*/
const commands::CommandMode& helpMode()
{
  static commands::CommandMode help_mode
    ("help: ",help_entry,help_exit);
  if (help_mode.empty()) // true upon first call
  {
    help_mode.add("q",commands::exitMode);
    insertTag(tagDict,"q",q_tag);

    help_mode.add("?",questionMark_h);
    insertTag(tagDict,"?",questionMark_tag);

    help_mode.add("intro",intro_h);
    insertTag(tagDict,"intro",intro_tag);

    // add help functions for the empty mode
    emptyhelp::addEmptyHelp(help_mode,tagDict);

    special::addSpecialHelp(help_mode,tagDict,emptymode::EmptymodeTag());
    test::addTestHelp(help_mode,tagDict,emptymode::EmptymodeTag());

  // add help functions for the main mode
    mainhelp::addMainHelp(help_mode,tagDict);

    special::addSpecialHelp(help_mode,tagDict,mainmode::MainmodeTag());
    test::addTestHelp(help_mode,tagDict,mainmode::MainmodeTag());

  // add help functions for the real mode
    realhelp::addRealHelp(help_mode,tagDict);

    special::addSpecialHelp(help_mode,tagDict,realmode::RealmodeTag());
    test::addTestHelp(help_mode,tagDict,realmode::RealmodeTag());

  // add help functions for the block mode
    blockmode::addBlockHelp(help_mode,tagDict);

    special::addSpecialHelp(help_mode,tagDict,blockmode::BlockmodeTag());
    test::addTestHelp(help_mode,tagDict,blockmode::BlockmodeTag());

  }
  return help_mode;
}



void intro_h()

{
  using namespace io;

  printFile(std::cerr,"intro_mess",io::MESSAGE_DIR);

  return;
}

void nohelp_h()

/*
  Synopsis: help message for when there is no help.
*/

{
  std::cerr << "sorry, no help available for this command" << std::endl;

  return;
}

}

}
