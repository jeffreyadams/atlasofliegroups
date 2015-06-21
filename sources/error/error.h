/*
  This is error.h
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2015 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef ERROR_H  /* guard against multiple inclusions */
#define ERROR_H

#include <stdexcept>

/******** type declarations *************************************************/

namespace atlas {

namespace error {

  struct FatalError;

  // i-o errors
  struct InputError;

  // overflow errors
  struct MemoryOverflow;
  struct NumericOverflow;
  struct NumericUnderflow;

  struct CartanError {};

  typedef FatalError PrimesError;

}

/******** type definitions **************************************************/

namespace error {

struct FatalError {
  void operator() (const char*);
};

// user abandoning input
struct InputError {
  void operator() (const char*);
};

// failure opening or writing output
struct OutputError {
  void operator() (const char*);
};

struct InnerClassError {
  void operator() (const char*);
};

struct MemoryOverflow {
  void operator() (const char*);
};

struct NumericOverflow {
  void operator() (const char*);
};

struct NumericUnderflow {
  void operator() (const char*);
};


// runtime errors for realex

struct Cayley_error : public std::runtime_error
{
  Cayley_error() : std::runtime_error("Undefined Cayley transform") {}
};

  } // |namespace error|

} // |namespace atlas|

#endif
