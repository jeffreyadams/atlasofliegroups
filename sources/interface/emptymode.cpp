/*
  This is emptymode.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "emptymode.h"

#include <iostream>

#include "helpmode.h"
#include "io.h"
#include "mainmode.h"
#include "special.h"
#include "test.h"
#include "version.h"

/****************************************************************************

  This file contains the definitions for the "empty" command mode, which
  is the startup mode for the program.

*****************************************************************************/

namespace atlas {

namespace {

  using namespace emptymode;

  class ThisMode:public commands::CommandMode {
  public:
  // constructors and destructors
    ThisMode();
    virtual ~ThisMode();
  // manipulators
    virtual const std::vector<const commands::CommandMode*>& next() const;
  };

  void this_entry();
  void this_exit();

  void printVersion();

  // functions for the predefined commands

  void help_f();
  void q_h();
  void type_f();

}

/****************************************************************************

        Chapter I -- The EmptyMode class

  Only one instance of this class will be constructed, on the first call
  of the function emptyMode().

*****************************************************************************/

namespace {

ThisMode::ThisMode()
  :CommandMode("empty: ",this_entry,this_exit)

{
  using namespace mainmode;

  add("help",help_f);
  add("q",q_h);
  add("qq",commands::exitInteractive);
  // the type command needs to be recognized in the empty mode, or else
  // it will trigger activation of the main mode and _then_ execute, which
  // leads to setting the type twice!
  add("type",type_f);

  special::addSpecialCommands(*this,EmptymodeTag());
  test::addTestCommands(*this,EmptymodeTag());
}

ThisMode::~ThisMode()

{}

const std::vector<const commands::CommandMode*>& ThisMode::next() const

/*
  Synopsis: returns the list of direct descendants of this mode; the list is
  constructed on first call.

  NOTE: apart from adhering to the general principle of lazy construction,
  the construction on first call avoids circular dependencies: in this way,
  the list is constructed with ThisMode object already fully constructed,
  so that the constructors for the descendant modes can call emptyMode() to
  recover their "d_prev" link.
*/

{
  using namespace commands;
  using namespace mainmode;

  static std::vector<const CommandMode*> nextList;

  if (nextList.size() == 0) {
    const CommandMode* modePtr = &mainMode();
    nextList.push_back(modePtr);
  }

  return nextList;
}

void this_entry()

{
  printVersion();
  return;
}

void this_exit()

{}

}

namespace emptymode {

const commands::CommandMode& emptyMode()

/*
  Synopsis: returns the ThisMode object. 

  It is constructed on first call.
*/

{
  static ThisMode mode;
  return mode;
}

}

/*****************************************************************************

        Chapter II --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/

namespace {

void help_f()

{
  helpmode::intro_h();
  activate(helpmode::helpMode());
  return;
}

void q_h()

{  
  io::printFile(std::cerr,"q.help",io::MESSAGE_DIR);
  return;
}

void type_f()

{
  using namespace commands;

  try {
    activate(mainmode::mainMode());
  }
  catch (EntryError) {
    return;
  }

  return;
}

}

/****************************************************************************

        Chapter III --- Utilities.

  This section contains some utility functions used in this module :

    - printVersion() : prints version info on startup;

*****************************************************************************/

namespace {

void printVersion()

/*
  Prints an opening message and the version number.
*/

{      
  std::cout << "This is " << version::NAME << " version " << version::VERSION 
	    << "." << std::endl;
  std::cout << 
    "Enter \"help\" if you need assistance." 
	    << std::endl << std::endl;

  return;
}

}

}
