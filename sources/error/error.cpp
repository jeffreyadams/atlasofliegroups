/*
  This is error.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "error.h"

#include <iostream>
#include <cstdlib>

/*****************************************************************************

        Chapter I -- The FatalError class

******************************************************************************/

namespace atlas {

namespace error {

// Calling |FatalError("why")| will actually terminate te program!
void FatalError::operator() (const char* mess)
{
  std::cerr << mess << std::endl;
  std::exit(0);
}


/*
  Execute the |InputError|.

  This prints a short message and returns.
*/
void InputError::operator() (const char* mess)
{
  std::cerr << mess << std::endl;
  return;
}


/*
  Executes the |OutputError|.

  This prints a short message and returns.
*/
void OutputError::operator() (const char* mess)
{
  std::cerr << mess << std::endl;
  return;
}


/*
  Execute the |MemoryOverflow| error.

  This prints a short message and returns.
*/
void MemoryOverflow::operator() (const char* mess)
{
  std::cerr << mess << std::endl;
  return;
}

/*
  Execute the |NumericOverflow| error.

  This prints a short message and returns.
*/
void NumericOverflow::operator() (const char* mess)
{
  std::cerr << mess << std::endl;
  return;
}


}

}
