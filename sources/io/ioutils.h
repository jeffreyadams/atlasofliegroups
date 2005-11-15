/*
  This is ioutils.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef IOUTILS_H  /* guard against multiple inclusions */
#define IOUTILS_H

#include <iosfwd>
#include <string>

/******** type declarations **************************************************/

namespace atlas {

namespace ioutils {

  class OutputFile;

}

/******** constant declarations **********************************************/

namespace ioutils {

  const size_t LineSize = 79;
  const size_t IndentSize = 4;

}

/******** function declarations *********************************************/

namespace ioutils {

unsigned long digits(unsigned long, unsigned long);

std::ostream& foldLine(std::ostream&, const std::string&, const char*,
		       const char* postHyphens = "", size_t h = IndentSize, 
		       size_t lineSize = LineSize);

std::istream& skipSpaces(std::istream&);

}

/******** type definitions ***************************************************/

namespace ioutils {

class OutputFile {
 private:
  std::ostream* d_stream;
  bool d_foutput;
 public:
  OutputFile();
  ~OutputFile();
  template<typename T> std::ostream& operator<< (const T& arg)
    {return *d_stream << arg;}
  operator std::ostream& () {return *d_stream;}
};

}

}

#endif
