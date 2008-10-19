/*
  This is ioutils.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "ioutils.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>

namespace atlas {

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

/*
  Synopsis: utility function to fold long lines of output nicely.

  The idea is that when the length of line exceeds |lineSize| characters, we
  print it over several lines, with well-chosen breakpoints. Here |h| is the
  indentation after the first line; |preHyphens| and |postHyphens| contain
  lists of acceptable breakpoints, either just before for |pre|, or just after
  for |post|. If no acceptable breakpoint exists, it breaks off the line
  brutally at lineSize. Breakpoints are strings that must be integrally
  present in the line to allow a break; multiple breakpoints can be specified
  in |pre| or |post|, separated by newline characters
*/
std::ostream& foldLine(std::ostream& strm, const std::string& line,
		       const char* preHyphens, const char* postHyphens,
		       size_t h, size_t lineSize)

{
  if (line.length() <= lineSize) // line fits on one line
    return strm << line;

  assert(h<lineSize); // so continuation lines will have room for something

  // break up breakpoint strings

  std::vector<std::string> pre, post;
  for (const char* p=preHyphens; *p!='\0'; ++p)
  {
    const char* q;
    for (q=p; *q!='\0' and *q!='\n'; ++q) {}
    assert(q!=p); // user should not specify empty separator
    pre.push_back(std::string(p,q));
    p= *q=='\n' ? q : q-1; // skip over string and over '\n'
  }
  for (const char* p=postHyphens; *p!='\0'; ++p)
  {
    const char* q;
    for (q=p; *q!='\0' and *q!='\n'; ++q) {}
    assert(q!=p); // user should not specify empty separator
    post.push_back(std::string(p,q));
    p= *q=='\n' ? q : q-1; // skip over string and over '\n'
  }

  // search for hyphenation point


  size_t point=0; // minimal and useless preakpoint value
  do
  {
    size_t old_break=point;
    size_t indent=old_break==0 ? 0 : h;

    for (size_t i=0; i<pre.size(); ++i)
    {
      size_t bp=line.rfind(pre[i],old_break+lineSize-indent);
      if (bp!=std::string::npos and bp>point) // get maximal non-|npos| value
	point=bp;
    }
    for (size_t i=0; i<post.size(); ++i)
    {
      size_t bp=line.rfind(pre[i],old_break+lineSize-indent);
      if (bp!=std::string::npos and bp+post[i].size()>point)
	point=bp+post[i].size();
    }

    if (point==old_break) // nothing was found
      point += lineSize-indent; // so break brutally

      strm << std::setw(indent) << ""
	   <<line.substr(old_break,point-old_break) << std::endl;
  }
  while (line.length()>point+lineSize-h);

  // print last line

  if (point<line.length()) //
    strm << std::setw(h) << "" << line.substr(point);

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
