/*
  This is commands.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright 2013 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#include "commands.h"

#include <cstdio>   // not obviously used, but appears helpful for Windows
#include <cstring>
#include <iostream>
#include <sstream>
#include <stack>
#include <cassert>

#include "error.h"

// extra defs for windows compilation -spc
#ifdef WIN32
#include "constants.h"
#endif

/*****************************************************************************

  This module contains the code for the command interface which we use during
  the development stage of the Atlas library. It is directly inspired by the
  interface for Coxeter, but re-written to be (hopefully) C++ conformant. In
  any case, i/o is now done in C++-style, and we also use the containers
  map and stack provided by C++, as well as the C++ string facility. Exceptions
  remain somewhat of a grey area here.

  The basic class is the CommandNode class, of which a few instances will be
  explitly constructed at initialisation in separate modules, such as
  mainmode.cpp. Each instance defines a set of recognized names, with
  associated action functions. These instances serve as base class for
  CommandTree, which includes additional links that allow defining a hierarchy
  of "modes", which is constructed in main.cpp. During this construction,
  commands are transitively inherited from ancestor to descendant node, unless
  the child defines a command of the same name as the parent does.

  At each point of time during execution, there is a stack of such modes, held
  in the local variable modeStack, the top of which is the currently active
  mode, with below it its ancestors (the code in the current module does not
  quite enforce the latter relation, but the action functions that may
  changethe state of the stack are defined so as to respect it). Some commands
  will lead to pushing a new mode on the stack; the entry function of the mode
  is then executed. Similarly, some commands (typically the "q" command) pop
  the mode stack; this results in executing the exit function of the mode.

  A command that was originally defined in some CommandNode can be called from
  that mode or from its descendants, and may assume that the entry function
  for the mode has been executed (but not subsequently the exit function).
  This allows using values obtained interactively from the user during the
  entry function, and so limits the amount of dialogue the individual
  mathematical commands need to execute to obtain the arguments they require.
  Ifa command is not defined in the currently active mode, it will be searched
  in descendent modes, and if found there the command will be executed after
  performing in order the entry functions of the modes on the path towards the
  descendant. The variable commandStack serves to retain the name of the
  command to be later executed while performing the entry functions. Some
  commands (like "type" when not called from the empty mode) explitily alter,
  after a user interaction, the values stored in an already active mode. This
  in general makes the values in descendants of that mode invalid; such
  commands therefore should pop any such descendant modes. Since a function
  that was propagated by inheritance does not know from which mode it was
  called, this operation currently requires such functions to be explicitly
  redefined in all decendants of the mode it affects, each redefinition
  executing the proper number of pop operations. This might be improved.

  Also, there is a unique help mode, which provides some help information for
  all commands available. it has no ancestors but can be entered from any
  other mode; since it contains no mode-changing commands except "q" (which
  pops the top mode), and it inherits no commands from other modes, the help
  mode effectively disables executing anything but its own the help commands.

******************************************************************************/

namespace atlas {
namespace commands {

// Local variables to the command.cpp module
namespace {

  std::stack<const CommandTree*> modeStack; // the stack of command modes;
                                             // the active mode is the top
                                             // of the stack

  std::stack<const char*> commandStack;     // pending commands;

  input::InputBuffer commandLine;           // the current command line

  bool runFlag;   // set to false to initiate shutdown of the program

  // running the command interface

  void ambiguous(const std::vector<const char*>&, const char*);
  void execute(const char* name, const CommandNode* mode);
  const char* getCommand(const CommandNode* mode); // get command from |mode|
  void getInteractive(std::string&, const char*);

  // auxiliary functions

  inline bool isEqual(const char* a, const char* b) {
    return std::strcmp(a,b)==0;
  }

  inline bool isInitial(const char* a, const char* b) {
    while (*a!='\0')
      if (*a++ != *b++)
	return false; // found a difference
    return true; // string |a| was a prefix of |b|
  }

} // |namespace|

/****************************************************************************

        Chapter I -- The CommandNode and CommandTree classes.

  These are the central class of the command module.

  - manipulators :

    - add(name,tag,action,rep) : add a new command to the mode;
    - findName : finds a name in the dictionary;
    - setAction(name,a) : sets the action associated to name;
    - setRepeat(name,b) : sets the repeat flag for name;

*****************************************************************************/



/*
  Constructor for the command mode class.
*/
CommandNode::CommandNode(const char* str,
			 void (*entry)(),
			 void (*exit)())
  : d_map()      // start out without commands
  , d_prompt(str)
  , d_entry(entry)
  , d_exit(exit)
{} // don't add any commands, so that |empty()| is true initially

/******** accessors *********************************************************/


/******** manipulators ******************************************************/

/*
  Synopsis: adds a new command during intial construction of a |CommandNode|
  NOTE: names should be unique within the Node, whence the |assert| below
*/
void CommandNode::add(const char* const name, const Command& command)
{
  std::pair<CommandDict::iterator,bool> p
    = d_map.insert(std::make_pair(name,command));

  assert(p.second); // then name should not be already present
}

/*
  Inserts the commands from source into dest. This is used when going to a
  "desecendant" mode, to inherit the commands defined for the parent mode.
  Here existing commands are not overridden!
*/
void CommandNode::addCommands(const CommandNode& source)
{
  for (CommandNode::const_iterator it = source.begin(); it!=source.end(); ++it)
    d_map.insert(*it); // ignore result
}




  //	********		CommandTree			********


/*
  Synopsis: find the command in the current mode or one of its descendant modes.
*/
CommandNode::const_iterator CommandTree::findName(const char* name) const
{
  const_iterator pos = find(name); // search in current mode

  if (pos != end())
    return pos;

  // if not find here, look recursively in descendant modes
  for (size_t j = 0; j < n_desc(); ++j) {
    const CommandTree& next = nextMode(j);
    pos = next.findName(name);
    if (pos != next.end())
      return pos;
  }

  return end();
}

/*
  Synopsis: tries to find name in mode, or in one of its descendants.
*/
CheckResult CommandTree::checkName(const char* name) const
{
  CommandNode::const_iterator pos = findName(name);

  if (pos == end())
    return NotFound; // command was not found

  // if we get to this point, there is at least one extension of name in mode

  std::vector<const char*> ext;
  extensions(ext,name);

  if (ext.size() > 1) { // ambiguous command
    ambiguous(ext,name);
    return Ambiguous;
  }

  // if we get to this point, the command is found in mode

  return Found;
}



/*
  Synopsis: puts into |e| the list of command names in mode and its
  descendants, that begin with |name|.

  Forwarded to the set-version, so that repetitions will be automatically
  weeded out.
*/
void CommandTree::extensions(std::vector<const char*>& e,
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
  Synopsis: adds to |e| the list of command names in mode and its descendants,
  that begin with |name|.
*/
void CommandTree::extensions(std::set<const char*,StrCmp>& e,
			     const char* name) const
{
  for (const_iterator it = find_prefix(name); it != end(); ++it)
    if (isInitial(name,it->first))
      e.insert(it->first);
    else
      break; // any element not starting with |name| ends search

  for (size_t j = 0; j < n_desc(); ++j) {
    const CommandTree& mode = nextMode(j);
    mode.extensions(e,name);
  }
}



/*
  Execute the command "name", which is a complete command name,
  defined either in the current mode or in one of its descendants.

  If name is found in the current mode, we simply execute it. Otherwise,
  we find out in which descendant mode it is found, we attempt the mode
  change, and in case of success, push command on the stack, so that it
  will be executed at the next loop in |run|.

  Failure to enter some mode (which will already have been reported), and all
  other kinds of errors (which we report here) are caught by this method.
*/
void CommandTree::execute(const char* name) const
{
  CommandNode::const_iterator pos = find(name);

  try
  {
    if (pos != end()) // the command was found in the current mode
    {
      const Command& command = pos->second;
      command();
    }
    else // we have to look in a submode
      for (size_t j = 0; j < n_desc(); ++j)
      {
	const CommandTree& next = nextMode(j);
	pos = next.findName(name); // this descends recursively
	if (pos != next.end()) // name is defined in a descendent of next
	{
	  next.activate();
	  commandStack.push(name); // retry command if mode entry successful
	  break; // only attempt to enter the first matching descendant
	}
      }
  }
  catch (commands::EntryError&) { // resume here after failure to enter mode
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e)
  {
    std::cerr << "input for command " << name; e(" aborted");
  }
  catch (std::exception& e)
  {
    std::cerr << "error occurred: " << e.what() << std::endl;
  }
  catch (...)
  {
    std::cerr << std::endl << "unidentified error occurred" << std::endl;
  }
} // |CommandTree::execute|


CommandTree::~CommandTree()
{
  for (unsigned int i=0; i<d_nextList.size(); ++i)
    delete d_nextList[i];
}


/*
  Attempts to activate the current mode tree, executing its entry function
  (which could throw an |EntryError|). Then push the mode tree onto stack.
*/
void CommandTree::activate() const
{
  entry(); // could throw an EntryError
  modeStack.push(this);
}

/*
  Runs an interactive session of the program with this mode as basic mode

  Gets commands from the user until it gets the "qq" command, at which time it
  returns control.

  It works as follows : get an input string from the user (leading whitespace
  is chopped off by default in C++); look it up in the current CommandNode, or,
  in case of failure, in its descendants; execute it if it is found, get new
  input otherwise.
*/
void CommandTree::run() const
{
  try
  {
    activate(); // activate mode, pushing our mode onto |modeStack|

    const char* name;

    for (runFlag = true; runFlag;)// exit is only through "qq"
    {                             // (or "q" in startup mode)

      const CommandTree& mode = *modeStack.top();  // get current active mode
      name = getCommand(&mode); // get user input, as edited by readline

      std::vector<const char*> ext;
      mode.extensions(ext,name); // find all commands with this prefix

      switch (ext.size())
      {
      case 0: // command was not found
	std::cout << name << ": not found" << std::endl;
	break;
      case 1: // command can be unambiguously completed
	mode.execute(ext[0]);
	break;
      default: // ambiguous command; execute it if there is an exact match
	if (isEqual(ext[0],name))
	  mode.execute(ext[0]);
	else
	  ambiguous(ext,name);
	break;
      }
    } // for(runFlag)
  }
  catch(EntryError) { // failed to activate our mode as initial mode
    return;
  }
} // |CommandTree::run|

CommandTree& CommandTree::add_descendant(const CommandNode& c)
{
  d_nextList.push_back(new CommandTree(c));
  CommandTree& child = *d_nextList.back();
  child.addCommands(*this); // will not overwrite commands existing in child
  return child;
}

/****************************************************************************

        Chapter II -- Function definitions.

  This section contains the definitions of the functions declared in
  commands.h :

    - currentLine() : get input buffer holding (remainder of) command line
    - currentMode() : return pointer to mode on top of the stack
    - exitMode() : exits the current mode;
    - drop_to(mode): exit modes until |mode| is left
    - exitInteractive() : sets runFlag to false;
    - relax_f() : does nothing;

*****************************************************************************/


/*
  Returns the InputBuffer in which the current command line was stored

  Explanation: the idea is to enable functions to pick off arguments from the
  command line, so that the user can type forward. Nobody actually does this.
*/
input::InputBuffer& currentLine()
{
  return commandLine;
}


/*
  Synopsis: returns the currently active mode.
  Mostly useful for communication with readline.
*/
const CommandTree* currentMode()
{
  return modeStack.top();
}


/*
  Exit the current mode, and pop it off the stack
*/
void exitMode()  // it is assumed that exit functions don't throw
{
  modeStack.top()->exit();
  modeStack.pop();
}

/* this function helps "mode changing" functions like "type" to be defined
   only twice: once in the parent mode where it just activates the mode, and
   once in the mode itself, where (upon successful obtention of new values for
   the mode) it drops to the current mode and replaces the values

   Avoiding the first definition altogether is not possible, since if only the
   definition in the mode itself would remain, it would have no way to know
   whether that mode had just been entered as a consequence of its own call
   (and in which case it would need to avoid immediately changing anything).
 */
void drop_to(const CommandTree& mode)
{
  while (modeStack.top()!= &mode)
  {
    exitMode();
    assert(not modeStack.empty()); // not finding |mode| is fatal
  }
}

/*
  Quit interactive mode, after exiting all active modes.

  This should be called only if none of the modes in the stack can throw on
  exit, which seems a reasonable assumption.
*/
void exitInteractive()
{
  while (modeStack.size())
    exitMode(); // applies to mode at top of stack, which is then popped

  runFlag = false; // this will cause the |run| loop to terminate
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
  std::cout << name << " : ambiguous (";

  for (std::vector<const char*>::const_iterator
	 it = ext.begin(); it != ext.end(); ++it)
    std::cout << (it == ext.begin() ? "" : ",") << *it;

  std::cout << ")" << std::endl;
}


/*
  Get a command interactively from the user.

  The actual input line is gotten through the readline library. For
  convenience we pack it into a global InputBuffer variable |commandLine|, in
  order to have a C++-like interaction.
*/
void getInteractive(std::string& name, const char* prompt)
{
  commandLine.getline(prompt);
  commandLine >> name;
}


/*
  Get the name of the next command.

  It is gotten either from the commandStack, if there are pending commands,
  or else interactively from the user.
*/
const char* getCommand(const CommandNode* mode)
{
  static std::string nameString; // semi-permanent storage for name
  const char* name;

  if (not commandStack.empty())
  {  // there is a command pending
    name = commandStack.top();
    commandStack.pop();
  } else
  {  // get command from user
    nameString.erase();
    getInteractive(nameString,mode->prompt());
    name = nameString.c_str();
  }

  return name;
}


} // namespace

} // |namespace commands|

} // namespace atlas
