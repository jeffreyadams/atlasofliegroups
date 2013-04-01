/*
  This is commands.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#ifndef COMMANDS_H  /* guard against multiple inclusions */
#define COMMANDS_H

#include "commands_fwd.h"

#include <map>
#include <set>
#include <vector>
#include <cstring>

#include "input.h"

/******** typedefs **********************************************************/

namespace atlas {

namespace commands {

  typedef void (*action_pointer)(); // pointer to void function, no arguments

  enum CheckResult { NotFound, PartialMatch, Ambiguous, Found };

}

/******** function and variable declarations *******************************/

namespace commands {

  input::InputBuffer& currentLine();
  const CommandTree* currentMode();

  extern TagDict tagDict; // used for helpmode

  void run_from(const CommandTree& initial_mode);
  void drop_to(const CommandTree& mode); // exit modes until |mode| is left
  void exitInteractive(); // used by 'qq' in empty mode
  void exitMode(); // used by 'q' in help mode
  inline void relax_f() {}

  void nohelp_h(); // this is used when no help command is defined (yet)

}

/******* Type definitions **************************************************/

namespace commands {

struct StrCmp // the string "less than" function disguised as function object
{ // a zero-size class with just one accessor method
  bool operator() (const char* a, const char* b) const
  { return strcmp(a,b) < 0; }
};

struct Command // a class wrapping a function pointer into a function object
{
  action_pointer action; // one data field: a function pointer
// constructor
explicit Command(action_pointer a) : action(a) {}; // store |a| as |action|
// accessor
  void operator() () const { action(); }
};

/*
  |CommandNode| instances are ephemerous objects that contain the attributes
  specific to a command mode; they will be passed to the constructor of the
  more long-lived |CommandTree| derived class, of which it (the copy) becomes
  a base object. The contruction of the ephemerous object is done by a
  non-member function, which repeatedly uses the public |add| method to
  populate the node with commands. We have chosen for the |CommandTree|
  objects for the modes to be stored in static rather than automatic
  variables, as this facilitates access by the various (action) functions that
  need to refer to explicitly names modes. This implies however that the
  lifetime of |CommandNode| instances is at static initialisation time, before
  the main program starts; this creates a potential for static initialisation
  fiasco. Therefore no calls of |CommandNode| methods should rely on the prior
  initialisation of static variables, except if these are defined in the same
  compilation unit as the |CommandTree| objects for the modes are (which is
  this unit, commands.cpp). Fortunately the arguments to the |add| method, a
  |C|-style string and a function pointer, can be given by constants that need
  no static initialisation.
 */
class CommandNode
{
  typedef std::map<const char*,Command,StrCmp> CommandDict;

  CommandDict d_map;
  const char* d_prompt;
  action_pointer d_entry, d_exit;

 public:

// Constructors and destructors
  CommandNode(const char*,
	      void (*entry)() = &relax_f,
	      void (*exit)() = &relax_f);

  ~CommandNode() {}

// accessors
  const char* prompt() const { return d_prompt; }

  void addCommands(const CommandNode& source); // inherit commands from |source|

  // bind |name| to |f| in current mode, possibly overriding previous
  // also provide a descriptive tag and possibly a help mode action
  void add(const char* const name, action_pointer f,
	   const char* const tag, action_pointer help_f = nohelp_h);

  void nohelp_add(const char* const name, action_pointer f)
   { add(name,Command(f)); }


  void exit() const { d_exit(); }


 protected: // methods for use by |CommandNode| or |CommandTree| objects only

  typedef CommandDict::const_iterator const_iterator;

  // bind |name| to |action| in current mode, possibly overriding previous
  void add(const char* const name, const Command& action);

  // basic search in our own command map
  const_iterator find(const char* name) const { return d_map.find(name); }

  const_iterator find_prefix(const char* name) const
    { return d_map.lower_bound(name); }

  void entry() const { d_entry(); }

  const_iterator begin() const { return d_map.begin(); }

  const_iterator end() const { return d_map.end(); }

}; // |class CommandNode|

class CommandTree : public CommandNode
{
 // list direct descendant modes
  std::vector<CommandTree*> d_nextList; // owned pointers

 public:
  explicit CommandTree(const CommandNode& root)
    : CommandNode(root) , d_nextList() {}
  ~CommandTree(); // since we have owned pointers

  // extend mode tree, returning reference to added (singleton) subtree
  CommandTree& add_descendant(const CommandNode& c); // make |c| a descendant

  void run() const; // start a command session based on this mode

  void activate() const; // enter (and push) this mode tree

  // the following is for use in guiding readline, and for ambiguity reporting
  std::vector<const char*> extensions(const char*) const;

 private:

  size_t n_desc() const { return d_nextList.size(); }

  const CommandTree& nextMode(unsigned int i) const { return *(d_nextList[i]); }

  bool has_descendant (const CommandTree* mode) const; // search decendants tree

  // look up |name|, return result in |status|, which is an in-out argument.
  // |status| should be anything except |Found| initially, and the result
  // and |where| are set only when it becomes |Found| or |PartialMatch|
  CommandNode::const_iterator look_up(const char* name,
				      CheckResult& status,
				      CommandTree const* & where) const;

  // implementation of the public |extensions| method, eliminates duplicates
  void extensions(std::set<const char*,StrCmp>&, const char*) const;

}; // |class CommandTree|

} // |namespace commands|

} // |namespace atlas|

#endif
