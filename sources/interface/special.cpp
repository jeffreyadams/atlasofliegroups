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

  enum SpecialMode {EmptyMode, MainMode, RealMode, BlockMode, numSpecialMode};
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


/*
  Adds to the empty mode the commands defined in special.cpp for that mode.

  NOTE : this is called in emptymode.cpp.
*/
void addSpecialCommands (commands::CommandMode& mode, emptymode::EmptymodeTag)
{
  if (specialMode == EmptyMode) {
    mode.add("special",special_f);

    // add additional commands here :

  }
}


/*
  Adds to the main mode the commands defined in special.cpp for that mode.

  NOTE : this is called in rmainmode.cpp. The special command itself is added
  to the empty mode by default. Change specialMode to MainMode if your
  definition requires the main command mode (i.e., the knowledge of a
  complex reductive group.)
*/
void addSpecialCommands(commands::CommandMode& mode, mainmode::MainmodeTag)
{
  if (specialMode == MainMode) {
    mode.add("special",special_f);

    // add additional commands here :

  }
}


/*
  Adds to the real mode the commands defined in special.cpp for that mode.

  NOTE : this is called in realmode.cpp. The special command itself is added
  to the empty mode by default. Change specialMode to RealMode if your
  definition requires the real command mode (i.e., the knowledge of a real
  form.)
*/
void addSpecialCommands(commands::CommandMode& mode, realmode::RealmodeTag)
{
  if (specialMode == RealMode) {
    mode.add("special",special_f);

    // add additional commands here :

  }
}

/*
  Adds to the block mode the commands defined in special.cpp for that mode.

  NOTE : this is called in blockmode.cpp. The special command itself is added
  to the empty mode by default. Change specialMode to BlockMode if your
  definition requires the block command mode (i.e., the knowledge of a real
  form and a dual real form.)
*/
void addSpecialCommands(commands::CommandMode& mode, blockmode::BlockmodeTag)
{
  if (specialMode == BlockMode) {
    mode.add("special",special_f);

    // add additional commands here :

  }
}


/*
  Adds the help functions for special commands which require the empty mode
  (i.e., nothing).

  NOTE: this is called in help.cpp.
*/
void addSpecialHelp(commands::CommandMode& mode, commands::TagDict& t,
		    emptymode::EmptymodeTag)
{
  using namespace commands;

  if (specialMode == EmptyMode) {
    mode.add("special",special_h);

    // add additional help commands here :

    // since "special" might be assigned to empty mode, define its tag here

    insertTag(t,"special",special_tag);

    // add additional tags here :

  }
}


/*
  Adds to mode the help commands for the commands defined in special.cpp, and
  which require the main mode; adds the tags to t.

  NOTE : this is called in help.cpp
*/
void addSpecialHelp(commands::CommandMode& mode, commands::TagDict& t,
		    mainmode::MainmodeTag)
{
  using namespace commands;

  if (specialMode == MainMode) {
    mode.add("special",special_h);

    // add additional help commands here :

    // since "special" might be assigned to main mode, define its tag here

    insertTag(t,"special",special_tag);

    // add additional tags here :

  }
}

/*
  Adds to mode the help commands for the commands defined in special.cpp, and
  which require the real mode; adds the tags to t.

  NOTE : this is called in realhelp.cpp
*/
void addSpecialHelp(commands::CommandMode& mode, commands::TagDict& t,
		    realmode::RealmodeTag)
{
  using namespace commands;

  if (specialMode == RealMode) {
    mode.add("special",special_h);

    // add additional help commands here :

    // since "special" might be assigned to real mode, define its tag here

    insertTag(t,"special",special_tag);

    // add additional tags here :

  }
}

/*
  Adds to mode the help commands for the commands defined in special.cpp, and
  which require the block mode; adds the tags to t.

  NOTE : this is called in blockhelp.cpp
*/
void addSpecialHelp(commands::CommandMode& mode, commands::TagDict& t,
		    blockmode::BlockmodeTag)
{
  using namespace commands;

  if (specialMode == BlockMode) {
    mode.add("special",special_h);

    // add additional help commands here :

    // since "special" might be assigned to block mode, define its tag here

    insertTag(t,"special",special_tag);

    // add additional tags here :

  }
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
