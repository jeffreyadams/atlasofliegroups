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
  const CommandMode* currentMode();
  void defaultError(const char*);
  void exitInteractive();
  void exitMode();
  void insertTag(TagDict&, const char*, const char*);
  void printTags(std::ostream&, const TagDict&);
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

class CommandMode
{
  typedef std::map<const char*,Command,StrCmp> CommandDict;

  std::vector<const CommandMode*> d_nextList; // direct descendant modes

  CommandDict d_map;
  const char* d_prompt;
  void (*d_entry)();
  void (*d_exit)();
  void (*d_error)(const char*);

 public:

  typedef CommandDict::iterator iterator;
  typedef CommandDict::const_iterator const_iterator;

// Constructors and destructors
  CommandMode(const char*,
	      void (*entry)() = &relax_f,
	      void (*exit)() = &relax_f,
	      void (*error)(const char*) = &defaultError);

  ~CommandMode() {}

// accessors
  const char* prompt() const { return d_prompt; }

  bool empty() const { return d_map.empty(); }

  void addCommands(const CommandMode& source); // inherit commands from |source|

  // bind |name| to |f| in current mode, possibly overriding previous
  void add(const char* const name, action_pointer f) { add(name,Command(f)); }

  // extend mode tree
  void add_descendant(CommandMode& c) // make |c| a descendant of |*this|
  { d_nextList.push_back(&c); }

  void run() const; // start a command session based on this mode

  void activate() const;

  void exit() const { d_exit(); }

  // the following is for use in guiding readline
  void extensions(std::vector<const char*>&, const char*) const;

 private: // methods for use by |CommandMode| instances only

  void error(const char* str) const { d_error(str); }

  // bind |name| to |action| in current mode, possibly overriding previous
  void add(const char* const name, const Command& action);

  // basic search in our own command map
  const_iterator find(const char* name) const { return d_map.find(name); }

  iterator find(const char* name) { return d_map.find(name); }

  // recursive search in this and descendant modes
  const_iterator findName(const char* name) const;

  // recursive search for |name| or unique command stqrting with |name|
  CheckResult checkName(const char* name) const;

  // execute |name| in this or descendant mode, maybe extending mode stack
  void execute(const char* name) const;

  void entry() const { d_entry(); }

  const std::vector<const CommandMode*>& next() const { return d_nextList; }

  size_t n_desc() const { return d_nextList.size(); }

  const CommandMode& nextMode(size_t j) const { return *(d_nextList[j]); }

  const_iterator begin() const { return d_map.begin(); }

  const_iterator end() const { return d_map.end(); }

  void extensions(std::set<const char*,StrCmp>&, const char*) const;

}; // class CommandMode

}

}

#endif
