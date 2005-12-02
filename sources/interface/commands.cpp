/*
  This is commands.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "commands.h"

#include <cstring>
#include <iostream>
#include <sstream>
#include <stack>

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

namespace commands {

  std::vector<const CommandMode*> CommandMode::d_empty;

}

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
    return !memcmp(a,b,strlen(a));
  }

  inline bool isInitial(const char* a, const char* b) {
    return !memcmp(a,b,strlen(a));
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

CommandMode::CommandMode(const char* str,
			 void (*entry)(),
			 void (*exit)(), 
			 void (*error)(const char*))
  :d_prompt(str),
   d_entry(entry),
   d_exit(exit),
   d_error(error),
   d_prev(0)

/*
  Constructor for the command mode class. It makes sure that the empty string
  is recognized.
*/

{
  add("",relax_f); // default action is to do nothing
}

CommandMode::~CommandMode()

{}

/******** accessors *********************************************************/

void CommandMode::extensions(std::set<const char*,StrCmp>& e, 
			     const char* name) const

/*
  Synopsis: adds to e the list of commandnames in mode and its descendants,
  that begin with name.
*/

{
  for (const_iterator pos = d_map.lower_bound(name); pos != d_map.end(); 
       ++pos) {
    if (isInitial(name,pos->first))
      e.insert(pos->first);
    else
      break;
  }

  for (size_t j = 0; j < next().size(); ++j) {
    const CommandMode& mode = nextMode(j);
    mode.extensions(e,name);
  }

  return;
}

void CommandMode::extensions(std::vector<const char*>& e, 
			     const char* name) const

/*
  Synopsis: puts in e the list of commandnames in mode and its descendants,
  that begin with name.

  Forwarded to the set-version, so that repetitions will be automatically
  weeded out.
*/

{
  std::set<const char*,StrCmp> es;

  extensions(es,name);
  e.clear();

  for (std::set<const char*,StrCmp>::const_iterator i = es.begin(); 
       i != es.end(); ++i)
    e.push_back(*i);

  return;
}

CommandMode::const_iterator CommandMode::findName(const char* name) const

/*
  Synopsis: finds the command in the current mode or one of its submodes.
*/

{
  const_iterator pos = find(name);

  if (pos != end())
    return pos;

  for (size_t j = 0; j < next().size(); ++j) {
    const CommandMode& next = nextMode(j);
    pos = next.findName(name);
    if (pos != next.end())
      return pos;
  }

  return end();
}

/******** manipulators ******************************************************/

void CommandMode::add(const char* const name, const Command& command)

/*
  Synopsis: adds a new command to the mode.

  The parameters have the following meaning :
    - name : name of the command;
    - command: the function to be executed by the command;

  NOTE: if the name was already present, we override it.
*/

{
  std::pair<const char* const, Command> v(name,command);

  std::pair<CommandDict::iterator,bool> p
    = d_map.insert(v);

  if (!p.second) {
    // name was already present; override!
    d_map.erase(p.first);
    d_map.insert(v);
  }

  return;
}

void CommandMode::setAction(const char* name, void (*a)())

/*
  Sets the action of the command associated to name to a.

  NOTE : it is assumed that name will be found in mode.
*/

{
  CommandMode::iterator pos = find(name);
  pos->second.action = a;

  return;
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
    - pushCommand(const char*) : pushes a command on commandStack;
    - quitMode() : exits the current mode;
    - relax_f() : does nothing;
    - run(CommandMode*) :  runs the program;

*****************************************************************************/

namespace commands {

void activate(const CommandMode& mode)

/*
  Attempts to activate the command mode mode, by executing its entry function.
  Could throw an EntryError.
*/

{
  mode.entry(); // could throw an EntryError

  const CommandMode* modePtr = &mode;
  modeStack.push(modePtr);

  return;
}

void addCommands(CommandMode& dest, const CommandMode& source)

/*
  Synopsis: inserts the commands from source into dest. This is used when
  going to a "derived" mode.

  NOTE: we do the insertion through add, so it will override existing
  commands.
*/

{
  CommandMode::const_iterator source_end = source.end();

  for (CommandMode::const_iterator i = source.begin(); i != source_end; ++i) {
    const char* name = i->first;
    Command command = i->second;
    dest.add(name,command);
  }

  return;
}

CheckResult checkName(const CommandMode& mode, const char* name)

/*
  Synopsis: tries to find name in mode, or in one of its descendants.
*/

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

const CommandMode* currentMode()

/*
  Synopsis: returns the currently active mode.
  Mostly useful for communication with readline.
*/

{
  return modeStack.top();
}

void defaultError(const char* str)

/*
  Synopsis: default error handler when str is not a string which is recognized 
  by the command tree.

  The default behaviour is to print the name of the read string an an error
  message.
*/

{
  std::cout << str << " : not found" << std::endl;
  return;
}

void exitInteractive()

/*
  Synopsis: quits interactive mode. 

  This should be called only when none of the modes in the stack can throw on 
  exit.
*/

{
  while (modeStack.size()) {
    exitMode();
  }

  runFlag = false;
  return;
}

void exitMode()

/*
  Synopsis: exits the current mode. 
*/

{
  const CommandMode* mode = modeStack.top();

  // it is assumed that exit functions don't throw
  mode->exit();

  modeStack.pop();

  return;
}

void insertTag(TagDict& t, const char* name, const char* tag)

/*
  Synopsis: associates tag with name in t.

  NOTE: this is the only way the sun CC compiler will accept it!
*/

{  
  TagDict::value_type v(name,tag);
  t.insert(v);

  return;
}

void printTags(std::ostream& strm, const TagDict& t)

/*
  Synopsis: outputs the list of commands with their attached tags.
*/

{  
  typedef TagDict::const_iterator I;

  for (I i = t.begin(); i != t.end(); ++i) {
    const char* key = i->first;
    const char* tag = i->second;
    strm << "  - " << key << " : " << tag << std::endl;
  }

  return;
}

void pushCommand(const char* name)

/*
  Synopsis: pushes name on commandStack.
*/

{
  commandStack.push(name);
  return;
}

void relax_f()

/*
  Synopsis: does nothing. 

  Useful as a default argument to functions which require a function pointer.
*/

{}

void run(const CommandMode& initMode)

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

{
  try {
    activate(initMode);
  }
  catch(EntryError) { // something is very wrong
    return;
  }

  const char* name;
  
  for (runFlag = true; runFlag;) {   // exit is only through "qq"
                                     // (or "q" in startup mode)

    const CommandMode* mode = modeStack.top();  // get current active mode
    name = getCommand(mode);

    std::vector<const char*> ext;
    mode->extensions(ext,name);

    switch (ext.size()) {
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

void ambiguous(const std::vector<const char*>& ext, const char* name)

/*
  Synopsis: outputs an informative message when name has more than one 
  completion in the dictionary.

  The message is name: ambiguous ( ... list of possible completions ... )
*/

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

  return;
}

void execute(const char* name, const CommandMode* mode)

/*
  Synopsis: executes the command "name".

  Precondition: name is defined either in the current mode or in one of its
  descendants.

  If name is found in the current mode, we simply execute it. Otherwise,
  we find out in which descendant mode it is found, we attempt the mode
  change, and in case of success, push command on the stack, so that it
  will be executed at the next loop in run().
*/

{
  CommandMode::const_iterator pos = mode->find(name);

  if (pos != mode->end()) { // the command was found in the current mode
    const Command& command = pos->second;
    command();
  } else // we have to look in a submode
    for (size_t j = 0; j < mode->next().size(); ++j) {
      const CommandMode& next = mode->nextMode(j);
      pos = next.findName(name);
      if (pos != next.end()) { // name is defined in a descendent of next
	try {
	  activate(next);
	}
	catch (EntryError) {
	  return;
	}
	pushCommand(name);
      }
    }

  return;
}

const char* getCommand(const CommandMode* mode)

/*
  Synopsis: gets the name of the next command.

  It is gotten either from the commandStack, if there are commands waiting
  to be processed, or interactively from the user.
*/

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

std::istream& getInteractive(std::istream& strm, std::string& name, 
			     const char* prompt)

/*
  Synopsis: gets a command interactively from the user.

  The actual input line is gotten through the readline library. For convenience
  we pack it into an InputBuffer (defined in input.h), in order to have a
  C++-like interaction.
*/

{    
  using namespace input;

  commandLine.getline(strm,prompt);
  commandLine >> name;

  return strm;
}

}

}
