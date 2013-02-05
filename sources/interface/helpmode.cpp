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
#include "emptyhelp.h"
#include "emptymode.h"
#include "io.h"
#include "mainhelp.h"
#include "mainmode.h"
#include "realhelp.h"
#include "realmode.h"
#include "blockmode.h"
#include "test.h"

/****************************************************************************

  This file contains the definition of the help mode.

*****************************************************************************/

namespace atlas {

namespace commands {

namespace {

  void help_entry();
  void help_exit();

  TagDict tagDict; // static, filled by |helpNode()|: initialisation danger

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
  printTags(std::cerr,tagDict);
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
  insertTag(tagDict,"q",q_tag);

  CommandNode result ("help: ",help_entry,help_exit);
  result.add("q",exitMode);
  insertTag(tagDict,"q",q_tag);

  result.add("?",questionMark_h);
  insertTag(tagDict,"?",questionMark_tag);

  result.add("intro",intro_h);
  insertTag(tagDict,"intro",intro_tag);

  // add help functions for the empty mode
  emptyhelp::addEmptyHelp(result,tagDict);
  test::addTestHelp(result,tagDict,EmptymodeTag());

  // add help functions for the main mode
  mainhelp::addMainHelp(result,tagDict);
  test::addTestHelp(result,tagDict,MainmodeTag());

  // add help functions for the real mode
  realhelp::addRealHelp(result,tagDict);
  test::addTestHelp(result,tagDict,RealmodeTag());

  // add help functions for the block mode
  addBlockHelp(result,tagDict);
  test::addTestHelp(result,tagDict,BlockmodeTag());

  return result;
}

// associate a tag with name in t.
void insertTag(TagDict& t, const char* name, const char* tag)
{
  t.insert(std::make_pair(name,tag));
}


// output the list of commands with their attached tags.
void printTags(std::ostream& strm, const TagDict& t)
{
  for (TagDict::const_iterator it = t.begin(); it != t.end(); ++it)
    strm << "  - " << it->first << " : " <<  it->second << std::endl;
}

void intro_h()
{
  io::printFile(std::cerr,"intro_mess",io::MESSAGE_DIR);
}

/*
  Synopsis: help message for when there is no help.
*/
void nohelp_h()
{
  std::cerr << "sorry, no help available for this command" << std::endl;
}

// variable that must be defined in this module to avoid initialisation fiasco
CommandTree help_mode(helpNode());

} // |namespace helpmode|

} // |namespace atlas|
