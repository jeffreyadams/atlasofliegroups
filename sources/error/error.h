/*
  This is error.h
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2015,2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef ERROR_H  /* guard against multiple inclusions */
#define ERROR_H

#include <stdexcept>
#include <string>

namespace atlas {
namespace error {

/******** type declarations *************************************************/


  struct FatalError;

  // i-o errors
  struct InputError;

  // internal overflow errors
  struct NumericOverflow;

  struct Cartan_error;
  struct Cayley_error;

/******** type definitions **************************************************/


// user abandoning input
struct InputError {
  void operator() (const char*);
};

// failure opening or writing output
struct OutputError {
  void operator() (const char*);
};

struct NumericOverflow {
  void operator() (const char*);
};


// runtime errors for atlas
struct Cartan_error : public std::runtime_error
{
  Cartan_error() : std::runtime_error("Illegal Cartan matrix") {}
  Cartan_error(const char* msg)
    : std::runtime_error(std::string("Illegal Cartan matrix: ")+msg) {}
};


struct Cayley_error : public std::runtime_error
{
  Cayley_error() : std::runtime_error("Undefined Cayley transform") {}
};

} // |namespace error|
} // |namespace atlas|

#endif
