/*
  This is commands.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef COMMANDS_H  /* guard against multiple inclusions */
#define COMMANDS_H

#include "commands_fwd.h"

#include <map>
#include <set>
#include <vector>

#include "input.h"

/******** typedefs **********************************************************/

namespace atlas {

namespace commands {

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
  void printTags(std::ostream&, const TagDict&);
  void pushCommand(const char*);
  void relax_f();
  void run(const CommandMode&);

}

/******* Type definitions **************************************************/

namespace commands {

class StrCmp { // a function object for string comparison
 public:
  bool operator() (const char* a, const char* b) const {
    return strcmp(a,b) < 0;
  }
};

struct Command {
  void (*action)();
// Constructors and destructors
  Command(void (*a)()):action(a) {
  };
  ~Command() {
  };
// accessors
  void operator() () const {
    action();
  }
};

class CommandMode {

 private:

  typedef std::map<const char*,Command,StrCmp> CommandDict;

  static std::vector<const CommandMode*> d_empty;

  CommandDict d_map;
  const char* d_prompt;
  void (*d_entry)();
  void (*d_exit)();
  void (*d_error)(const char*);

 protected:

  const CommandMode* d_prev;

 public:

  typedef CommandDict::iterator iterator;
  typedef CommandDict::const_iterator const_iterator;

// Constructors and destructors
  CommandMode(const char*,
	      void (*entry)() = &relax_f,
	      void (*exit)() = &relax_f, 
	      void (*error)(const char*) = &defaultError);

  virtual ~CommandMode();

// accessors
  const_iterator begin() const {
    return d_map.begin();
  }

  const_iterator end() const {
    return d_map.end();
  }

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
    // dummy value for modes without descendants
    return d_empty;
  }

  const CommandMode& nextMode(size_t j) const {
    return next()[j][0];
  }

  const CommandMode& prev() const {
    return *d_prev;
  }

// manipulators
  void add(const char* name, void (*action)()) {
    add(name,Command(action));
  }

  void add(const char*, const Command&);

  iterator find(const char* name) {
    return d_map.find(name);
  }

  void setAction(const char*, void (*)());

};

}

}

#endif
