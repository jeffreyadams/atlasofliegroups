/*
  This is helpmode.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
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
#include "special.h"
#include "test.h"

/****************************************************************************

  This file contains the definition of the help mode.

*****************************************************************************/

namespace atlas {

namespace {
  using namespace helpmode;

  class ThisMode:public commands::CommandMode {
  public:
  // constructors and destructors
    ThisMode();
    virtual ~ThisMode();
  };

  void this_entry();
  void this_exit();

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

ThisMode::ThisMode()
  :CommandMode("help: ",this_entry,this_exit)

{
  using namespace commands;
  using namespace emptyhelp;
  using namespace emptymode;
  using namespace mainhelp;
  using namespace mainmode;
  using namespace realhelp;
  using namespace realmode;

  add("q",commands::exitMode);
  insertTag(tagDict,"q",q_tag);
 
  add("?",questionMark_h);
  insertTag(tagDict,"?",questionMark_tag);

  add("intro",intro_h);
  insertTag(tagDict,"intro",intro_tag);

// add help functions for the empty mode
  addEmptyHelp(*this,tagDict);

  special::addSpecialHelp(*this,tagDict,EmptymodeTag());
  test::addTestHelp(*this,tagDict,EmptymodeTag());

// add help functions for the main mode
  addMainHelp(*this,tagDict);

  special::addSpecialHelp(*this,tagDict,MainmodeTag());
  test::addTestHelp(*this,tagDict,MainmodeTag());

// add help functions for the real mode
  addRealHelp(*this,tagDict);

  special::addSpecialHelp(*this,tagDict,RealmodeTag());
  test::addTestHelp(*this,tagDict,RealmodeTag());
}

ThisMode::~ThisMode()

{}

void this_entry()

{
  std::cerr << "enter \"?\" for a list of available commands; "
	    << "\"q\" to exit help mode"
	    << std::endl;
  std::cerr << "entering a command name will print an explanatory message"
	    << std::endl;
}

void this_exit()

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

const commands::CommandMode& helpMode()

/*
  Synopsis: returns the help mode.

  It is constructed on the first call.
*/

{
  static ThisMode mode;
  return mode;
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
