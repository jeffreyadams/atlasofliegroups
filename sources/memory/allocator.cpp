/*!
\file
\brief Implementation for the template class Allocator<T,P>, which
  allocates memory for objects of type T in pools of type P.
*/
/*
  This is allocator.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations 

  For license information see the LICENSE file
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

  // void reportAllocation(size_t n, size_t p, const typenumber::TypeData& td)
  void reportAllocation(size_t n, size_t a, size_t p)

/*!
  \brief Allocation report.

  Preconditions: n is the number of objects being allocated; p is the
  instance number of the pool for objects of this type; and a is the
  size of each object.
*/

{
  using namespace typestring;

  //  memlog() << "allocating " << n*td.size() << " bytes from pool # " << p;
  memlog() << "allocating " << n << "*" << a << " = " << n*a
           << " bytes from pool # " << p << std::endl;
  //  memlog() << " (" << n << " " << name(td.number()) << ")" << std::endl;

  return;
}

void reportConstruction(size_t n, size_t i, const typenumber::TypeData& td)

/*!
  \brief Construction report.

  Also makes a convenient breakpoint for debugging.
*/

{
  using namespace typestring;

  std::cerr << "construction of allocator #" << n
	    << " (" << i << " instances)";
  std::cerr << " (" << name(td.number()) << ")" << std::endl;
}

void reportCopyConstruction(size_t n, size_t i, const typenumber::TypeData& td)

/*!
  \brief Copy construction report.

  Also makes a convenient breakpoint for debugging.
*/

{
  using namespace typestring;

  std::cerr << "copy construction of allocator #" << n
	    << " (" << i << " instances)";
  std::cerr << " (" << name(td.number()) << ")" << std::endl;
}

// void reportDeallocation(size_t n, size_t p, const typenumber::TypeData& td)
void reportDeallocation(size_t n, size_t a, size_t p)
/*!
  \brief Deallocation report.

  Preconditions: n is the number of objects being deallocated; p is
  the instance number of the pool for objects of this type; and a is
  the size of each object.
*/

{
  using namespace typestring;

  // memlog() << "deallocating " << n*td.size() << " bytes from pool #" << p;

  memlog() << "deallocating " << n << "*" << a << " = " << n*a
           << " bytes from pool #" << p <<  std::endl;

  // memlog() << " (" << n << " " << name(td.number()) << ")" << std::endl;

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

/*!
  \brief Heterogeneous copy construction report.

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

/*!
  \brief Returns the stream on which to output memory management messages.
*/

{
  static std::ofstream strm("memory.out");

  return strm;
}

}

}
