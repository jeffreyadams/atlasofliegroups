/*
  This is error.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef ERROR_H  /* guard against multiple inclusions */
#define ERROR_H

/******** type declarations *************************************************/

namespace atlas {

namespace error {

  struct FatalError;

  // i-o errors
  struct InputError;

  // error for interactive constructors
  struct InnerClassError;

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

}

}

#endif
