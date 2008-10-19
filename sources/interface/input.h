/*
  This is input.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For copyright and license information see the LICENSE file
*/

#ifndef INPUT_H  /* guard against multiple inclusions */
#define INPUT_H

#include <iosfwd>
#include <sstream>

#ifndef NREADLINE
#include <readline/history.h>
#endif

/******** type declarations ************************************************/

namespace atlas {

namespace input {

  class InputBuffer;

#ifndef NREADLINE
  class HistoryBuffer;
#else
 typedef InputBuffer HistoryBuffer;
#endif

}

/******** function declarations ********************************************/

namespace input {

  bool hasQuestionMark(InputBuffer&);
  void initReadLine();

}

/******** type definitions ************************************************/

namespace input {

class InputBuffer:public std::istringstream {

 public:

// constructors and destructors
  InputBuffer():std::istringstream()
    {}

  explicit InputBuffer(const std::string& str)
    : std::istringstream(str)
    {}

  virtual ~InputBuffer()
    {}

// manipulators
  virtual std::istream& getline(std::istream&, const char* prompt = "",
				bool toHistory = true);

  virtual void reset();

  virtual void reset(std::streampos);
};

#ifndef NREADLINE

/* A class like |InputBuffer|, but each instance creates a local version of
   the history during its lifetime. This is managed in such a way that it is
   invisible to other instances, or to that static history record maintained
   in the history library.
*/
class HistoryBuffer : public InputBuffer {

 private:

  HISTORY_STATE state; // records our branch of history

 public:

// constructors and destructors
  HistoryBuffer();

  explicit HistoryBuffer(const std::string& str);

  virtual ~HistoryBuffer();

  virtual std::istream& getline(std::istream&, const char* prompt = "",
				bool toHistory = true);

};
#endif

}

}

#endif
