/*
  This is commands.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

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

  void activate(const CommandMode&);
  void addCommands(CommandMode&, const CommandMode&);
  CheckResult checkName(const CommandMode&, const char*);
  input::InputBuffer& currentLine();
  const CommandMode* currentMode();
  void defaultError(const char*);
  void exitInteractive();
  void exitMode();
  void insertTag(TagDict&, const char*, const char*);
  void printTags(std::ostream&, const TagDict&);
  inline void relax_f() {}
  void run(const CommandMode&);

}

/******* Type definitions **************************************************/

namespace commands {

struct StrCmp // the string comparison function disguised as function object
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

class CommandMode {

 private:

  typedef std::map<const char*,Command,StrCmp> CommandDict;

  std::vector<const CommandMode*> d_nextList; // direct successor modes

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

  virtual ~CommandMode() {}

// accessors
  const_iterator begin() const {
    return d_map.begin();
  }

  const_iterator end() const {
    return d_map.end();
  }

  bool empty() const { return d_map.empty(); }

  const_iterator find(const char* name) const {
    return d_map.find(name);
  }

  const char* prompt() const {
    return d_prompt;
  }

  void entry() const {
    d_entry();
  }

  void error(const char* str) const {
    d_error(str);
  }

  void exit() const {
    d_exit();
  }

  void extensions(std::set<const char*,StrCmp>&, const char*) const;

  void extensions(std::vector<const char*>&, const char*) const;

  const_iterator findName(const char* name) const;

  virtual const std::vector<const CommandMode*>& next() const {
    return d_nextList;
  }

  size_t n_desc() const { return d_nextList.size(); }

  const CommandMode& nextMode(size_t j) const {
    return *(d_nextList[j]);
  }

// manipulators
  void add(const char* const name, void (*action)()) {
    add(name,Command(action));
  }

  void add(const char* const, const Command&);

  void add_descendant(CommandMode& c) // make |c| a descendant of |*this|
  { d_nextList.push_back(&c); }

  iterator find(const char* name) {
    return d_map.find(name);
  }

  void setAction(const char*, void (*)());

}; // class CommandMode

}

}

#endif
