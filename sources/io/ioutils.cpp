/*
  This is ioutils.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "ioutils.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "input.h"

/*****************************************************************************

        Chapter I -- The OutputFile class

  This is a well-known C++ trick : a file which is opened by its constructor,
  and closed by its destructor.

******************************************************************************/

namespace atlas {

namespace ioutils {

OutputFile::OutputFile()

{  
  using namespace input;

  std::string name;
  InputBuffer buf;

  buf.getline(std::cin,"Name an output file (hit return for stdout): ",
	       false);
  buf >> name;

  if (name.empty()) {
    d_foutput = false;
    d_stream = &std::cout;
  }
  else {
    d_foutput = true;
    d_stream = new std::ofstream(name.c_str());
  }
}


OutputFile::~OutputFile()

/*
  Closes *d_stream if it is not std::cout.
*/

{
  if (d_foutput)
    delete d_stream;
}

}

/*****************************************************************************

        Chapter II -- Functions declared in ioutils.h

  ... fill in here when it is stable ...

******************************************************************************/

namespace ioutils {

unsigned long digits(unsigned long a, unsigned long b)

/*
  Synopsis: returns the number of digits of a in base b.

  Precondition: b > 1;
*/

{
  size_t d = 1;

  while (a >= b) {
    ++d;
    a /= b;
  }

  return d;
}

std::ostream& foldLine(std::ostream& strm, const std::string& line, 
		       const char* preHyphens, const char* postHyphens,
		       size_t h, size_t lineSize)

/*
  Synopsis: utility function to fold long lines of output nicely.

  The idea is that when the length of line exceeds lineSize characters, we
  print it over several lines, with well-chosen breakpoints. Here h is the
  indentation after the first line; preHyphens and postHyphens contain
  lists of acceptable breakpoints, either just before for pre, or just
  after for post. If no acceptable breakpoint exists, it breaks off the line 
  brutally at lineSize.
*/

{
  if (line.length() <= lineSize) // line fits on one line
    return strm << line;

  // search for hyphenation point

  size_t bp = line.find_last_of(preHyphens,lineSize);

    if (bp == std::string::npos) { // look for a post-hypenation
      bp = line.find_last_of(postHyphens,lineSize-1);
      if (bp == std::string::npos) // break brutally
	bp = lineSize - h;
      else
	++bp;
    }

  // output first line

  strm << std::string(line,0,bp) << std::endl;

  // print continuation lines

  size_t p = bp;

  for (; p + lineSize < line.length() + h; p += bp) {
    bp = line.find_last_of(preHyphens,p+lineSize-h);
    if (bp == std::string::npos) { // look for a post-hypenation
      bp = line.find_last_of(postHyphens,p+lineSize-h-1);
      if (bp == std::string::npos) // break brutally
	bp = lineSize - h;
      else
	bp -= (p-1);
    } else
      bp -= p;
    strm << std::setw(h) << "";
    strm << std::string(line,p,bp) << std::endl;
  }

  // print last line

  strm << std::setw(h) << "";
  strm << std::string(line,p);

  return strm;
}

std::istream& skipSpaces(std::istream& strm)

/*
  Synopsis: advances stream to the next non-space character.

  Explanation: spaces are the characters recognized by isspace().

  NOTE: this should be a library function, but I couldn't find it!
*/

{
  while (isspace(strm.peek())) // advance strm by one character
    strm.get();

  return strm;
}

}

}
