/*
  This is input.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef INPUT_H  /* guard against multiple inclusions */
#define INPUT_H

#include <iosfwd>
#include <sstream>

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

  explicit InputBuffer(const std::string& str):std::istringstream(str) 
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
class HistoryBuffer:public InputBuffer {

 private:

  void* d_history;

 public:

// constructors and destructors
  HistoryBuffer();

  explicit HistoryBuffer(const std::string& str);

  virtual ~HistoryBuffer();
  
};
#endif

}

}

#endif
