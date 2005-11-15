/*
  This is allocator_def.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#include "allocator.h"

#include <iostream>
#include <fstream>

#include "typestring.h"

namespace atlas {

namespace {

  std::ostream& memlog();

}

/*****************************************************************************

        Chapter I -- functions declared in allocator.h

*****************************************************************************/

namespace allocator {

void reportAllocation(size_t n, size_t p, const typenumber::TypeData& td)

/*
  Synopsis: allocation report.

  Preconditions: n is the number of objects being allocated;
*/

{
  using namespace typestring;

  memlog() << "allocating " << n*td.size() << " bytes from pool # " << p;
  memlog() << " (" << n << " " << name(td.number()) << ")" << std::endl;

  return;
}

void reportConstruction(size_t n, size_t i, const typenumber::TypeData& td)

/*
  Synopsis: construction report.

  Also makes a convenient breakpoint for debugging.
*/

{
  using namespace typestring;

  std::cerr << "construction of allocator #" << n
	    << " (" << i << " instances)";
  std::cerr << " (" << name(td.number()) << ")" << std::endl;
}

void reportCopyConstruction(size_t n, size_t i, const typenumber::TypeData& td)

/*
  Synopsis: copy construction report.

  Also makes a convenient breakpoint for debugging.
*/

{
  using namespace typestring;

  std::cerr << "copy construction of allocator #" << n
	    << " (" << i << " instances)";
  std::cerr << " (" << name(td.number()) << ")" << std::endl;
}

void reportDeallocation(size_t n, size_t p, const typenumber::TypeData& td)

/*
  Synopsis: deallocation report.

  Preconditions: n is the number of objects being deallocated;
*/

{
  using namespace typestring;

  memlog() << "deallocating " << n*td.size() << " bytes from pool #" << p;
  memlog() << " (" << n << " " << name(td.number()) << ")" << std::endl;

  return;
}

void reportDestruction(size_t n, size_t i, const typenumber::TypeData& td)

{
  using namespace typestring;

  std::cerr << "destruction of allocator #" << n
	    << " (" << i << " instances)";
  std::cerr << " (" << name(td.number()) << ")" << std::endl;
}

void reportDestructionError(size_t n, size_t p)

{
  std::cerr << "allocation count not zero for pool #" << p
	    << " (contains " << n << " objects)" << std::endl;

  return;
}

void reportHeterogeneousCopyConstruction(size_t n)

/*
  Synopsis: hetrogeneous copy construction report.

  Also makes a convenient breakpoint for debugging.
*/

{
  std::cerr << "heterogeneous copy construction of instance #" << n 
	    << std::endl;
}

}

/*****************************************************************************

        Chapter II -- Auxiliary functions

*****************************************************************************/

namespace {

std::ostream& memlog()

/*
  Synopsis: returns the stream on which to output memory management messages.
*/

{
  static std::ofstream strm("memory.out");

  return strm;
}

}

}
