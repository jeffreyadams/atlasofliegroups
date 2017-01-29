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
#include <cassert>

#include "error.h"
#include "io.h"

#include "emptymode.h"
#include "mainmode.h"
#include "realmode.h"
#include "blockmode.h"
#include "reprmode.h"
#include "helpmode.h"

#include "interactive.h" // to clear input buffer

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

  The basic class is the |CommandNode| class, of which a few instances will be
  explitly constructed at initialisation in separate modules, such as
  mainmode.cpp. Each instance defines a set of recognized names, with
  associated action functions. These instances serve as base class for
  |CommandTree|, which includes additional links that allow defining a
  hierarchy of "modes", which is constructed in main.cpp. During this
  construction, commands are transitively inherited from ancestor to
  descendant node, unless the child defines a command of the same name as the
  parent does.

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
  If a command is not defined in the currently active mode, it will be
  searched in descendent modes, and if found there the command will be
  executed after performing in order the entry functions of the modes on the
  path towards the descendant. The variable commandStack serves to retain the
  name of the command to be later executed while performing the entry
  functions. Some commands (like "type" when not called from the empty mode)
  explitily alter, after a user interaction, the values stored in an already
  active mode. This in general makes the values in descendants of that mode
  invalid; such commands therefore should pop any such descendant modes. Since
  a function that was propagated by inheritance does not know from which mode
  it was called, this operation currently requires such functions to be
  explicitly redefined in all decendants of the mode it affects, each
  redefinition executing the proper number of pop operations. This might be
  improved.

  Also, there is a unique help mode, which provides some help information for
  all commands available. it has no ancestors but can be entered from any
  other mode; since it contains no mode-changing commands except "q" (which
  pops the top mode), and it inherits no commands from other modes, the help
  mode effectively disables executing anything but its own the help commands.

******************************************************************************/

namespace atlas {
namespace commands {

  TagDict tagDict; // static, filled by |helpNode()|

  // this mode must be constructed before the others, because the functions
  // called for their construction also populate |help_mode|
  CommandTree help_mode(helpNode());

// static mode variables, declared here to control order of intialisation

  CommandTree empty_mode(emptyNode());
  CommandTree& main_mode = empty_mode.add_descendant(mainNode());
  CommandTree& real_mode = main_mode.add_descendant(realNode());
  CommandTree& block_mode = real_mode.add_descendant(blockNode());
  CommandTree& repr_mode =  real_mode.add_descendant(reprNode());


// Local variables to the command.cpp module
namespace {

  // the stack of command modes; the active mode is at the end
  std::vector<const CommandTree*> modeStack; // non-owned pointers

  input::InputBuffer commandLine;           // the current command line

  bool runFlag;   // set to false to initiate shutdown of the program

  const char* command_name = NULL; // name of currently executing command

  // auxiliaries for running the command interface
  const char* getCommand(const char* prompt); // get command string
  void ambiguous(const std::vector<const char*>&, const char*); // complain

  inline bool isEqual(const char* a, const char* b)
  { return std::strcmp(a,b)==0; }

  inline bool isInitial(const char* a, const char* b)
  {
    while (*a!='\0')
      if (*a++ != *b++)
	return false; // found a difference
    return true; // string |a| was a prefix of |b|
  }

} // |namespace|

/****************************************************************************

        Chapter I -- The CommandNode and CommandTree classes.

  These are the central classes of the command module.

  - manipulators :

    - add(name,action) : add a new command to the mode;
    - addCommands(parent) : inherit commands defined in parent

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
  Synopsis: adds a new command during initial construction of a |CommandNode|
  NOTE: names should be unique within the Node, whence the |assert| below
*/
void CommandNode::add(const char* const name, const Command& command)
{
  std::pair<CommandDict::iterator,bool> p
    = d_map.insert(std::make_pair(name,command));

  assert(p.second); // name should not be doubly defined in one |CommandNode|
  ndebug_use(p);
}

void CommandNode::add(const char* const name, action_pointer f,
		      const char* const tag, action_pointer help_f)
{
  add(name,Command(f));
  std::pair<TagDict::iterator,bool> p =
    tagDict.insert(std::make_pair(name,tag));
  if (p.second) // don't try to add after first definition in help mode
    help_mode.add(name,Command(help_f));
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


/* the following function exists to remedy the fact that |look_up| does not
   record the path to the match it has found (supposing it did); its calling
   interface is already compicated as it is. So we ask our children which one
   knows the mode that was returned. The interrogation is repeated recursively
   down, but in practice very few modes exist at all, and this is very quick.
 */
bool CommandTree::has_descendant (const CommandTree* mode) const
{
  if (mode==this)
    return true; // that was easy
  for (unsigned int i=0; i<n_desc(); ++i)
    if (nextMode(i).has_descendant(mode))
      return true;

  return false;
}


CommandNode::const_iterator CommandTree::look_up
 (const char* name, CheckResult& status, CommandTree const* & where) const
{
  CommandNode::const_iterator result; // might remain undefined
  CommandNode::const_iterator it = find_prefix(name); // in our |mode| only
  if (it != end() and isInitial(name,it->first)) // anything found here?
  {
    if (isEqual(name,it->first)) // got exact match here
      return status=Found, where=this, it;

    // record partial match, then scan ahead for any further partial matches
    if (status!=NotFound) // we found a partial match, and had at least one
      status = Ambiguous;
    else
    {
      status=PartialMatch; where=this; result=it;
      if (++it != end() and isInitial(name,it->first))
	status = Ambiguous;
    }
  }

  bool not_found_before = status==NotFound;
  // now look in descendant modes; if any has exact match, return that
  for (unsigned int i=0; i<n_desc(); ++i)
  {
    it = nextMode(i).look_up(name,status,where);
    if (status==Found) // found exact match, |where| has been set
      return it; // |where| has been set by recursive |lookup|
    if (not_found_before and status!=NotFound) // got a first partial match
    {
      not_found_before = false; // only one can be the first
      result=it; // export result from first successful recursive |look_up|
    }
  }

  // now |status| can be anything except |Found|
  return result; // value should not be used if |status==NotFound|
}



/*
  add to |e| the set of command names in mode and its descendants,
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

  for (size_t j = 0; j < n_desc(); ++j)
  {
    const CommandTree& mode = nextMode(j);
    mode.extensions(e,name);
  }
}

/*
  Put into |e| the list of command names, in our mode and its descendants,
  that begin with |name|.

  Forwarded to the set-version, so that repetitions will be automatically
  weeded out.
*/
std::vector<const char*> CommandTree::extensions(const char* name) const
{
  std::set<const char*,StrCmp> es;
  extensions(es,name);

  return std::vector<const char*>(es.begin(),es.end());
}



CommandTree::~CommandTree()
{
  for (unsigned int i=0; i<d_nextList.size(); ++i)
    delete d_nextList[i];
}

// grow the tree from root to leaves, inheriting down after copying node |c|
CommandTree& CommandTree::add_descendant(const CommandNode& c)
{
  d_nextList.push_back(new CommandTree(c));
  CommandTree& child = *d_nextList.back();
  child.addCommands(*this); // will not overwrite commands existing in child
  return child;
}

/*
  Attempts to activate the current mode tree, executing its entry function
  (which could throw an |EntryError|). Then push the mode tree onto stack.
*/
void CommandTree::activate() const
{
  entry(); // could throw an EntryError
  modeStack.push_back(this);
}

/*
  Runs an interactive session of the program with this mode as basic mode

  Gets commands from the user until it gets the "qq" command, at which time it
  returns control.

  It works as follows : get an input string from the user (leading whitespace
  is chopped off by default in C++); look it up in the current CommandNode or
  in its descendants (a prefix match is allowed if it is unique); execute it
  if it is found, otherwise complain and get new input
*/
void CommandTree::run() const
{
  try { activate(); }
  catch(EntryError) { // we've got off to a very bad start
      std::cout << "Internal error, failed to enter initial mode!" << std::endl;
  }

  runFlag = true;
  const char* name; // needed in catch blocks
  while (runFlag) // exit through |exitInteractive|, i.e., "qq"
  {
    const CommandTree* mode = modeStack.back();  // get current active mode
    name = getCommand(mode->prompt()); // user input, edited by readline
    try
    {
      CheckResult status = NotFound; // must initialise this before the call
      CommandTree const* where;
      CommandTree::const_iterator it = mode->look_up(name,status,where);
      switch (status)
      {
      case NotFound:
	// as a last resort try to locate |name| in descendant of an ancestor
	for (int i=modeStack.size()-1; i-->0; )
	{
	  it = modeStack[i]->look_up(name,status,where);
	  if (status==Found or status==PartialMatch)
	  {
	    mode=modeStack[i];
	    drop_to(*mode); // leave intermediate modes
	    goto found_case;
	  }
	}
	std::cout << name << ": not found" << std::endl;
	break;
      case Ambiguous: // then report all commands with this prefix
	ambiguous(mode->extensions(name),name);
	break;
      case Found: case PartialMatch: // these behave identically now
      found_case:
	while (mode!=where)
	  for (unsigned int i=0; i<mode->n_desc(); ++i)
	    if (mode->nextMode(i).has_descendant(where))
	    {
	      mode->nextMode(i).activate(); // perform entry function
	      mode=&mode->nextMode(i); // change to this descendant mode
	      assert(modeStack.back()==mode); // it should have been pushed
	      break; // from |for| loop, see whether |mode==where| next
	    }
	command_name = it->first; // set pointer to name of command
	it->second(); // finally execute the command in its proper mode
      }
    } // try
    catch (EntryError) {} // resume loop after user abort in |activate|
    catch (error::InputError& e) // user abort in actual command execution
    { std::cerr << "input for command '" << name; e("' aborted");  }
    catch (error::MemoryOverflow& e) { e("error: memory overflow"); }
    catch (std::exception& e)
    {
      std::cerr << "error occurred in command '" << name << "': "
		<< e.what() << std::endl;
    }
    catch (...)
    {
      std::cerr << std::endl << "unidentified error occurred" << std::endl;
    }
  } // |for(runFlag)|
} // |run_from|


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
  return modeStack.back();
}


/*
  Exit the current mode, and pop it off the stack
*/
void exitMode()  // it is assumed that exit functions don't throw
{
  modeStack.back()->exit();
  modeStack.pop_back();
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
  while (modeStack.back()!= &mode)
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

// help message for when there is no help file.
void nohelp_h()
{
  std::cerr << "sorry, no help available for this command" << std::endl;
}

void std_help()
{
  std::string file_name = command_name;
  file_name += ".help";
  io::printFile(std::cerr,file_name.c_str(),io::MESSAGE_DIR);
}

void use_tag() // lacking a help file, we might just print the tag as help
{
  TagDict::const_iterator it = tagDict.find(command_name);
  if (it != tagDict.end())
    std::cerr << it->second << std::endl;
  else
    nohelp_h();
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
  Output an informative message when name has more than one completion in the
  dictionary.

  The message is name: ambiguous ( ... list of possible completions ... )
*/
void ambiguous(const std::vector<const char*>& ext, const char* name)
{
  std::cout << name << ": ambiguous (";

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
const std::string& getInteractive(std::string& name, const char* prompt)
{
  name.erase();  // probably redundant
  interactive::common_input().str(""); // clear any previous type-ahead
  commandLine.getline(prompt);
  commandLine >> name; // get word (up to whitespace) into |name|
  return name;
}


/*
  Get the name of the next command.

  It is gotten either from the commandStack, if there are pending commands,
  or else interactively from the user.
*/
const char* getCommand( const char* prompt)
{
  static std::string nameString; // semi-permanent: most recent command name
  return getInteractive(nameString,prompt).c_str();
}


} // |namespace|

} // |namespace commands|

} // |namespace atlas|
