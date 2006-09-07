/*
  This is special.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "special.h"

#include <iostream>

#include "emptymode.h"
#include "io.h"
#include "mainmode.h"
#include "realmode.h"

/*****************************************************************************

  This module is intended to facilitate the addition os extra commands to
  the program. A "special" command is provided, which can be edited by the
  user, and serves as a model for the addition of other commands.

******************************************************************************/

namespace atlas {

namespace {

  using namespace special;

  enum SpecialMode {EmptyMode, MainMode, RealMode, numSpecialMode};
  SpecialMode specialMode = EmptyMode;

  // functions for the special commands

  void special_f();

  // help functions

  void special_h();

  // tags

  const char* special_tag = "user-definable command";

}

/*****************************************************************************

        Chapter I -- Functions declared by special.h

  This section defines the functions declared in special.h :

    - addSpecialCommands() : adds the special commands to the main command 
      tree;
    - addSpecialHelp() : adds help functionality;

******************************************************************************/

namespace special {

void addSpecialCommands (commands::CommandMode& mode, emptymode::EmptymodeTag)

/*
  Adds to the empty mode the commands defined in special.cpp for that mode.

  NOTE : this is called in emptymode.cpp.
*/

{
  if (specialMode == EmptyMode) {
    mode.add("special",special_f);

    // add additional commands here :

  }

  return;
}

void addSpecialCommands(commands::CommandMode& mode, mainmode::MainmodeTag)

/*
  Adds to the main mode the commands defined in special.cpp for that mode.

  Adds to the real mode the commands defined in special.cpp for that mode.

  NOTE : this is called in realmode.cpp. The special command itself is added
  to the empty mode by default. Change specialMode to MainMode if your 
  definition requires the main command mode (i.e., the knowledge of a
  complex reductive group.)
*/

{
  if (specialMode == MainMode) {
    mode.add("special",special_f);

    // add additional commands here :

  }

  return;
}

void addSpecialCommands(commands::CommandMode& mode, realmode::RealmodeTag)

/*
  Adds to the real mode the commands defined in special.cpp for that mode.

  NOTE : this is called in realmode.cpp. The special command itself is added
  to the empty mode by default. Change specialMode to RealMode if your 
  definition requires the real command mode (i.e., the knowledge of a real 
  form.)
*/

{
  if (specialMode == RealMode) {
    mode.add("special",special_f);

    // add additional commands here :

  }

  return;
}

void addSpecialHelp(commands::CommandMode& mode, commands::TagDict& t, 
		    emptymode::EmptymodeTag)

/*
  Adds the help functions for special commands which require the empty mode 
  (i.e., nothing).

  NOTE: this is called in help.cpp.
*/

{
  using namespace commands;

  if (specialMode == EmptyMode) {
    mode.add("special",special_h);

    // add additional help commands here :

  /******** tags **********************************************************/

    insertTag(t,"special",special_tag);

    // add additional tags here :

  }

  return;
}

void addSpecialHelp(commands::CommandMode& mode, commands::TagDict& t, 
		    mainmode::MainmodeTag)

/*
  Adds to mode the help commands for the commands defined in special.cpp, and 
  which require the main mode; adds the tags to t.

  NOTE : this is called in help.cpp
*/

{
  using namespace commands;

  if (specialMode == MainMode) {
    mode.add("special",special_h);

    // add additional help commands here :

  /******** tags **********************************************************/

    insertTag(t,"special",special_tag);

    // add additional tags here :

  }

  return;
}

void addSpecialHelp(commands::CommandMode& mode, commands::TagDict& t, 
		    realmode::RealmodeTag)

/*
  Adds to mode the help commands for the commands defined in special.cpp, and 
  which require the real mode; adds the tags to t.

  NOTE : this is called in realhelp.cpp
*/

{
  using namespace commands;

  if (specialMode == RealMode) {
    mode.add("special",special_h);

    // add additional help commands here :

  /******** tags **********************************************************/

    insertTag(t,"special",special_tag);

    // add additional tags here :

  }

  return;
}

}

/*****************************************************************************

        Chapter II -- Functions for the special commands

******************************************************************************/

namespace {

void special_f()

/*
  The response to the "special" command is to execute this function. By
  changing the contents (and recompiling the program of course---just say
  "make") you can make it do whatever you like.
*/

{
  std::cout << "not defined for now" << std::endl;
  return;
}

}

/*****************************************************************************

        Chapter III -- Help functions

******************************************************************************/

namespace {

void special_h()

/*
  Help function for the special command. You may want to edit the message
  if you redefine the command.
*/

{  
  io::printFile(std::cerr,"special.help",io::MESSAGE_DIR);
  return;
}

}

}
