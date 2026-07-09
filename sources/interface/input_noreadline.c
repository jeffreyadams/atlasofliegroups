/*
  This is input_noreadline.c

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For copyright and license information see the LICENSE file
*/

/*
  This file is one of the two variants of input.cpp. It is called .c so
  that no object file is created by the makefile.
*/

#include "input.h"

#include <iostream>

#include "commands.h"

namespace atlas {

namespace {

const char* readLine(const char* prompt = "", bool toHistory = true);

}

/*****************************************************************************

        Chapter I -- The InputBuffer class

  An InputBuffer object can hold a line of input, and can be used as in
  |std::istringstream| object (since it is publicly derived from that class),
  in particular one can read from it using the 'source >> variable' syntax.
  In addition, it provides a |getline| method that will fill the string that
  it holds from standard input; this calls readline unless NREADLINE was set.

******************************************************************************/

namespace input {

/*!
  Forget string in input buffer, read a new line, and position at start
*/
void InputBuffer::getline(const char* prompt, bool toHistory)

{
  const char* line = readLine(prompt,toHistory); // non-owned pointer

  str(line==NULL ? std::cout << "qq\n","qq" : line); // 'qq' at end of input
  reset();
}

void InputBuffer::reset()

/*
  Synopsis: rewinds the stream for re-reading.
*/

{
  clear();
  seekg(0,std::ios_base::beg);

  return;
}


/*
  Synopsis: resets the stream to position pos.

  Clears the flags as well; the idea is to undo a peek-forward operation.
*/
void InputBuffer::reset(std::streampos pos)
{
  clear();
  seekg(pos);
}

} // |namespace input|

/*****************************************************************************

        Chapter III -- Functions declared in input.h
                    -- version _without_ readline

******************************************************************************/

namespace input {

bool hasQuestionMark(InputBuffer& buf)

/*
  Synopsis: looks if the next character assignment from buf would return '?'
*/

{
  std::streampos pos = buf.tellg();
  char x = 0;
  buf >> x;
  buf.reset(pos);

  return x == '?';
}

void initReadLine()

{}

}

/*****************************************************************************

        Chapter IV -- Local functions
                   -- version _without_ readline

******************************************************************************/


namespace {

/*
  Synopsis: gets a line of input.
*/
const char* readLine (const char* prompt, bool toHistory)
{
  static std::string line;

  std::cout << prompt;
  std::getline(std::cin,line);

  return line.c_str();
}

} // |namespace|

} // |namespace atlas|
