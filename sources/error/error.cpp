/*
  This is error.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "error.h"

#include <iostream>

/*****************************************************************************

        Chapter I -- The FatalError class

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace error {

void FatalError::operator() (const char* mess)

/*
  Synopsis: executes the FatalError.

  This prints a short message and exits.
*/

{
  std::cerr << mess << std::endl;
  exit(0);
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

void InnerClassError::operator() (const char* mess)

/*
  Synopsis: executes the InnerClassError.

  This prints a short message and returns.
*/

{
  std::cerr << mess << std::endl;
  return;;
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
