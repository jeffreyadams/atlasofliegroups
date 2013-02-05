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

  enum CheckResult { Ambiguous, Found, NotFound, numCheckResults };

}

/******** function declarations ********************************************/

namespace commands {

  input::InputBuffer& currentLine();
  const CommandTree* currentMode();

  void drop_to(const CommandTree& mode); // exit modes until |mode| is left
  void exitInteractive();
  void exitMode();
  inline void relax_f() {}

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
Command(action_pointer a) : action(a) {}; // store |a| as |action|
// accessor
  void operator() () const { action(); }
};

class CommandNode
{
  typedef std::map<const char*,Command,StrCmp> CommandDict;

  CommandDict d_map;
  const char* d_prompt;
  void (*d_entry)();
  void (*d_exit)();

 public:

  typedef CommandDict::iterator iterator;
  typedef CommandDict::const_iterator const_iterator;

// Constructors and destructors
  CommandNode(const char*,
	      void (*entry)() = &relax_f,
	      void (*exit)() = &relax_f);

  ~CommandNode() {}

// accessors
  const char* prompt() const { return d_prompt; }

  void addCommands(const CommandNode& source); // inherit commands from |source|

  // bind |name| to |f| in current mode, possibly overriding previous
  void add(const char* const name, action_pointer f) { add(name,Command(f)); }

  void exit() const { d_exit(); }


 protected: // methods for use by |CommandNode| or |CommandTree| objects only

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

  // extend mode tree, returning reference to added (singleton) subtree
  CommandTree& add_descendant(const CommandNode& c); // make |c| a descendant
  ~CommandTree(); //

  void run() const; // start a command session based on this mode

  void activate() const; // enter (and push) this mode tree

  // the following is for use in guiding readline
  void extensions(std::vector<const char*>&, const char*) const;

 private:

  size_t n_desc() const { return d_nextList.size(); }

  const CommandTree& nextMode(size_t j) const { return *(d_nextList[j]); }

  // recursive search in this and descendant modes
  CommandNode::const_iterator findName(const char* name) const;

  // recursive search for |name| or unique command stqrting with |name|
  CheckResult checkName(const char* name) const;

  // execute |name| in this or descendant mode, maybe extending mode stack
  void execute(const char* name) const;

  void extensions(std::set<const char*,StrCmp>&, const char*) const;

}; // |class CommandTree|

} // |namespace commands|

} // |namespace atlas|

#endif
