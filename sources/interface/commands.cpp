/*
  This is commands.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For copyright and license information see the LICENSE file
*/

#include "commands.h"

#include <cstring>
#include <iostream>
#include <sstream>
#include <stack>

#include "error.h"

/*****************************************************************************

  This module contains the code for the command interface which we use during
  the development stage of the Atlas library. It is directly inspired by the
  interface for Coxeter, but re-written to be (hopefully) C++ conformant. In
  any case, i/o is now done in C++-style, and we also use the containers
  map and stack provided by C++, as well as the C++ string facility. Exceptions
  remain somewhat of a grey area here.

  The basic class is the CommandMode class. This contains a set of recognized
  names, for which it will execute corresponding functions. At each point of
  time, there is a stack of such modes, the local variable modeStack, the top
  of which is the currently active mode. Some commands will lead to pushing
  a new mode on the stack; the entry function of the mode is then executed.
  Similarly, some commands (typically the "q" command) pop the mode stack; this
  results in executing the exit function of the mode.

  Also, to each mode is associated a help mode, which usually provides some
  help information for each command available in that mode. It could do more
  in some cases.

******************************************************************************/

namespace atlas {

namespace { // declarations private to commands.cpp

  using namespace commands;

  std::stack<const CommandMode*> modeStack; // the stack of command modes; the
                                            // active mode is the top of the
                                            // stack

  std::stack<const char*> commandStack;     // lets us postpone some commands;

  input::InputBuffer commandLine;           // the current command line

  bool runFlag;

  // running the command interface

  void ambiguous(const std::vector<const char*>&, const char*);
  void execute(const char* name, const CommandMode* mode);
  const char* getCommand(const CommandMode* mode);
  std::istream& getInteractive(std::istream&, std::string&, const char*);

  // auxiliary functions

  inline bool isEqual(const char* a, const char* b) {
    return std::memcmp(a,b,strlen(a))==0;
  }

  inline bool isInitial(const char* a, const char* b) {
    return std::memcmp(a,b,strlen(a))==0;
  }

}

/****************************************************************************

        Chapter I -- The CommandMode class.

  This is the central class of the command module.

  The following member functions are defined :

  - constructors and destructors :

    - CommandMode(prompt,entry,exit,error) : constructor;
    - ~CommandMode() : destructor;

  - accessors :

    - prompt : prints the prompt;

  - manipulators :

    - add(name,tag,action,rep) : add a new command to the mode;
    - findName : finds a name in the dictionary;
    - setAction(name,a) : sets the action associated to name;
    - setRepeat(name,b) : sets the repeat flag for name;

*****************************************************************************/

namespace commands {


/*
  Constructor for the command mode class.
*/
CommandMode::CommandMode(const char* str,
			 void (*entry)(),
			 void (*exit)(),
			 void (*error)(const char*))
  :d_nextList(), // start out without successors
   d_map(),      // and without commands
   d_prompt(str),
   d_entry(entry),
   d_exit(exit),
   d_error(error)
{} // don't add any commands, so that |empty()| is true initially

/******** accessors *********************************************************/


/*
  Synopsis: adds to |e| the list of command names in mode and its descendants,
  that begin with |name|.
*/

void CommandMode::extensions(std::set<const char*,StrCmp>& e,
			     const char* name) const
{
  for (const_iterator pos = d_map.lower_bound(name); pos != d_map.end();
       ++pos) {
    if (isInitial(name,pos->first))
      e.insert(pos->first);
    else
      break;
  }

  for (size_t j = 0; j < n_desc(); ++j) {
    const CommandMode& mode = nextMode(j);
    mode.extensions(e,name);
  }
}


/*
  Synopsis: puts into |e| the list of command names in mode and its
  descendants, that begin with |name|.

  Forwarded to the set-version, so that repetitions will be automatically
  weeded out.
*/
void CommandMode::extensions(std::vector<const char*>& e,
			     const char* name) const
{
  std::set<const char*,StrCmp> es;

  extensions(es,name);
  e.clear();

  for (std::set<const char*,StrCmp>::const_iterator i = es.begin();
       i != es.end(); ++i)
    e.push_back(*i);

  return;
}


/*
  Synopsis: finds the command in the current mode or one of its submodes.
*/
CommandMode::const_iterator CommandMode::findName(const char* name) const
{
  const_iterator pos = find(name); // search in current mode

  if (pos != end())
    return pos;

  // if not find here, look recursively in descendant modes
  for (size_t j = 0; j < n_desc(); ++j) {
    const CommandMode& next = nextMode(j);
    pos = next.findName(name);
    if (pos != next.end())
      return pos;
  }

  return end();
}

/******** manipulators ******************************************************/


/*
  Synopsis: adds a new command to the mode.

  The parameters have the following meaning :
    - name : name of the command;
    - command: the function to be executed by the command;

  NOTE: if the name was already present, we override it.
*/
void CommandMode::add(const char* const name, const Command& command)
{
  std::pair<const char* const, Command> v(name,command);

  std::pair<CommandDict::iterator,bool> p
    = d_map.insert(v);

  if (not p.second) // then name was already present; override its command!
    p.first->second=command; // and don't touch |p.first->first==name|
}


/*
  Sets the action of the command associated to name to a.

  NOTE : it is assumed that name will be found in mode.
*/
void CommandMode::setAction(const char* name, void (*a)())
{
  CommandMode::iterator pos = find(name);
  pos->second.action = a;
}

}

/****************************************************************************

        Chapter II -- Function definitions.

  This section contains the definitions of the functions declared in
  commands.h :

    - activate(CommandMode*) : makes mode the active mode;
    - defaultError(const char*) : default error handler for command reading;
    - checkName(CommandMode*, const char*) : tries to find name in mode;
    - default_help() : the default help function;
    - exitInteractive() : sets runFlag to false;
    - printTags(stream&,map<std::string,const char*>&) : prints the tag list;
    - quitMode() : exits the current mode;
    - relax_f() : does nothing;
    - run(CommandMode*) :  runs the program;

*****************************************************************************/

namespace commands {


/*
  Attempts to activate the command mode |mode|, by executing its entry
  function, which could throw an |EntryError|. If not, push |mode| onto stack.
*/
void activate(const CommandMode& mode)
{
  mode.entry(); // could throw an EntryError
  modeStack.push(&mode);
}


/*
  Synopsis: inserts the commands from source into dest. This is used when
  going to a "desecendant" mode, to inherit the commands defined for
  the parent mode.

  NOTE: we do the insertion through add, so it will override existing
  commands. Therefore to redefine commands with the same name in the
  descendant mode, add them after making the call to |addCommands|.
*/
void addCommands(CommandMode& dest, const CommandMode& source)
{
  for (CommandMode::const_iterator i = source.begin(); i != source.end(); ++i)
    dest.add(i->first,i->second);
}


/*
  Synopsis: tries to find name in mode, or in one of its descendants.
*/
CheckResult checkName(const CommandMode& mode, const char* name)
{
  CommandMode::const_iterator pos = mode.findName(name);

  if (pos == mode.end()) { // command was not found
    return NotFound;
  }

  // if we get to this point, there is at least one extension of name in mode

  std::vector<const char*> ext;
  mode.extensions(ext,name);

  if (ext.size() > 1) { // ambiguous command
    ambiguous(ext,name);
    return Ambiguous;
  }

  // if we get to this point, the command is found in mode

  return Found;
}

input::InputBuffer& currentLine()

/*
  Synopsis: returns the current command line.

  Explanation: the idea is to enable functions to pick off arguments from the
  command line, so that the user can type forward.
*/

{
  return commandLine;
}


/*
  Synopsis: returns the currently active mode.
  Mostly useful for communication with readline.
*/
const CommandMode* currentMode()
{
  return modeStack.top();
}


/*
  Synopsis: default error handler, which is called when |str| is not a string
  which is recognized by the command tree.

  It prints the name of |str| and an error message.
*/
void defaultError(const char* str)
{
  std::cout << str << ": not found" << std::endl;
}


/*
  Synopsis: quits interactive mode, after exiting all active modes.

  This should be called only if none of the modes in the stack can throw on
  exit, which seems a reasonable assumption.
*/
void exitInteractive()

{
  while (modeStack.size())
    exitMode(); // applies to mode at top of stack, which is then popped

  runFlag = false; // this will cause the |run| loop to terminate
}


/*
  Synopsis: exits the current mode.
*/
void exitMode()  // it is assumed that exit functions don't throw
{
  modeStack.top()->exit();
  modeStack.pop();
}

// The |TagDict| type is addressed below, but no instance is defined here.

/*
  Synopsis: associates tag with name in t.

  NOTE: this is the only way the sun CC compiler will accept it!
*/

void insertTag(TagDict& t, const char* name, const char* tag)
{
  TagDict::value_type v(name,tag); // pair of strings acceptable to |TagDict|
  t.insert(v);
}


/*
  Synopsis: outputs the list of commands with their attached tags.
*/
void printTags(std::ostream& strm, const TagDict& t)
{
  typedef TagDict::const_iterator I;

  for (I i = t.begin(); i != t.end(); ++i) {
    const char* key = i->first;
    const char* tag = i->second;
    strm << "  - " << key << " : " << tag << std::endl;
  }
}

/*
  Synopsis: runs an interactive session of the program.

  Gets commands from the user until it gets the "qq" command, at which time it
  returns control.

  It works as follows : get an input string from the user (leading whitespace
  is chopped off by default in C++); look it up in the current CommandMode, or,
  in case of failure, in its descendants; execute it if it is found, get new
  input otherwise.

  The initMode argument is the startup mode of the interactive session;
  that is, the session starts by executing the entry function of initMode.
*/
void run(const CommandMode& initMode)
{
  try
  {
    activate(initMode); // enter mode, pushing it onto |modeStack|

    const char* name;

    for (runFlag = true; runFlag;)// exit is only through "qq"
    {                             // (or "q" in startup mode)

      const CommandMode* mode = modeStack.top();  // get current active mode
      name = getCommand(mode);

      std::vector<const char*> ext;
      mode->extensions(ext,name);

      switch (ext.size())
      {
      case 0: // command was not found
	mode->error(name);
	break;
      case 1: // command can be unambiguously completed
	execute(ext[0],mode);
	break;
      default: // ambiguous command; execute it if there is an exact match
	if (isEqual(ext[0],name))
	  execute(ext[0],mode);
	else
	  ambiguous(ext,name);
	break;
      }
    } // for(runFlag)
  }
  catch(EntryError) { // something is very wrong
    return;
  }
}

}

/*****************************************************************************

        Chapter IV -- Local function definitions

  This section contains the definition of the private functions for this
  module :

    - ambiguous(ext,name) : handles multiple completions;
    - getIntreractive(strm, name, prompt) : get a command from the user;

******************************************************************************/

namespace {


/*
  Synopsis: outputs an informative message when name has more than one
  completion in the dictionary.

  The message is name: ambiguous ( ... list of possible completions ... )
*/
void ambiguous(const std::vector<const char*>& ext, const char* name)
{
  typedef std::vector<const char*>::const_iterator I;

  std::cout << name << " : ambiguous (";

  for (I i = ext.begin(); i != ext.end();) {
    std::cout << *i;
    ++i;
    if (i != ext.end())
      std::cout << ",";
  }

  std::cout << ")" << std::endl;
}


/*!
  \brief Execute the command "name". Also install default error handling.

  Precondition: name is defined either in the current mode or in one of its
  descendants.

  If name is found in the current mode, we simply execute it. Otherwise,
  we find out in which descendant mode it is found, we attempt the mode
  change, and in case of success, push command on the stack, so that it
  will be executed at the next loop in |run|.
*/
void execute(const char* name, const CommandMode* mode)
{
  CommandMode::const_iterator pos = mode->find(name);

  try
  {
    if (pos != mode->end()) // the command was found in the current mode
    {
      const Command& command = pos->second;
      command();
    }
    else // we have to look in a submode
      for (size_t j = 0; j < mode->n_desc(); ++j)
      {
	const CommandMode& next = mode->nextMode(j);
	pos = next.findName(name);
	if (pos != next.end()) // name is defined in a descendent of next
	{
	  activate(next);
	  commandStack.push(name); // retry command if mode entry successful
	  break; // only attempt to enter the first matching descendant
	}
      }
  }
  catch (commands::EntryError&) { // silently ignore failure to enter mode
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }
}


/*
  Synopsis: gets the name of the next command.

  It is gotten either from the commandStack, if there are commands waiting
  to be processed, or interactively from the user.
*/
const char* getCommand(const CommandMode* mode)
{
  static std::string nameString;
  const char* name;

  if (commandStack.size()) { // there is a command to process
    name = commandStack.top();
    commandStack.pop();
  } else {  // get command from user
    nameString.erase();
    getInteractive(std::cin,nameString,mode->prompt());
    name = nameString.c_str();
  }

  return name;
}


/*
  Synopsis: gets a command interactively from the user.

  The actual input line is gotten through the readline library. For convenience
  we pack it into an InputBuffer (defined in input.h), in order to have a
  C++-like interaction.
*/
std::istream& getInteractive(std::istream& strm, std::string& name,
			     const char* prompt)
{
  using namespace input;

  commandLine.getline(strm,prompt);
  commandLine >> name;

  return strm;
}

} // namespace

} // namespace atlas
