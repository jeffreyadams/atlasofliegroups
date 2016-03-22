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

void InputError::operator() (const char* mess)

/*
  Synopsis: executes the InputError.

  This prints a short message and returns.
*/

{
  std::cerr << mess << std::endl;
  return;
}

void OutputError::operator() (const char* mess)

/*
  Synopsis: executes the OutputError.

  This prints a short message and returns.
*/

{
  std::cerr << mess << std::endl;
  return;
}

void MemoryOverflow::operator() (const char* mess)

/*
  Synopsis: executes the MemoryOverflow error.

  This prints a short message and returns.
*/

{
  std::cerr << mess << std::endl;
  return;
}

void NumericOverflow::operator() (const char* mess)

/*
  Synopsis: executes the NumericOverflow error.

  This prints a short message and returns.
*/

{
  std::cerr << mess << std::endl;
  return;
}

void NumericUnderflow::operator() (const char* mess)

/*
  Synopsis: executes the NumericUnderflow error.

  This prints a short message and returns.
*/

{
  std::cerr << mess << std::endl;
  return;
}

}

}
