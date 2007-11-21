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

  ... explain here when it is stable ...

******************************************************************************/

namespace input {

std::istream& InputBuffer::getline(std::istream& is, const char* prompt,
				   bool toHistory)

/*
  Synopsis: reads 
*/

{
  std::string line = readLine(prompt,toHistory);

  str(line);
  reset();

  return is;
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

void InputBuffer::reset(std::streampos pos)

/*
  Synopsis: resets the stream to position pos.

  Clears the flags as well; the idea is to undo a peek-forward operation.
*/

{  
  clear();
  seekg(pos);

  return;
}

}

/*****************************************************************************

        Chapter III -- Functions declared in input.h
                    -- version _without_ readline

  ... explain here when it is stable ...

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

  ... explain here when it is stable ...

******************************************************************************/

namespace {

const char* readLine (const char* prompt, bool toHistory)

/*
  Synopsis: gets a line of input.
*/

{
  static std::string line;

  std::cout << prompt;
  std::getline(std::cin,line);

  return line.c_str();
}

}

}
