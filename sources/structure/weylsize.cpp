/*!
\file
\brief Implementation of  functions for namespace weylsize.

Functions to compute the size of a Weyl group (stored as a prime
factorization, to avoid overflow problems).
  This is weylsize.cpp
*/
/*
  This is weylsize.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "weylsize.h"

#include <cassert>

#include "intutils.h"
#include "lietype.h"

/*****************************************************************************

        Chapter I -- Functions declared in weylsize.h

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace weylsize {

void weylSize(size::Size& c, const lietype::LieType& lt)

/*!
  Synopsis: puts in c the size of the Weyl group with Lie type lt.
*/

{
  using namespace size;

  c.reset();

  for (size_t j = 0; j < lt.size(); ++j) {
    Size a;
    weylSize(a,lt[j]);
    c *= a;
  }
  
  return;
}

void weylSize(size::Size& c, const lietype::SimpleLieType& slt)

/*!
  Synopsis: returns the size of the Weyl group with Lie type slt.

  NOTE: we use some "magic numbers" even though we know how to compute the
  size of the group from first principles; it seems futile to construct the
  Weyl group just to know its size.
*/

{
  using namespace intutils;
  using namespace lietype;

  TypeLetter x = type(slt);
  unsigned long r = rank(slt);

  c.reset();
  
  switch (x) {
  case 'A':
    size::factorial(c,r+1);
    return;
  case 'B':
  case 'C': {
    size::factorial(c,r);
    c[0] += r;
    return;
  }
  case 'D': {
    size::factorial(c,r);
    c[0] += r-1;
    return;
  }
  case 'E':
    switch (r) {
    case 6:
      c[0] = 7;
      c[1] = 4;
      c[2] = 1;
      return;
    case 7:
      c[0] = 10;
      c[1] = 4;
      c[2] = 1;
      c[3] = 1;
      return;
    case 8:
      c[0] = 14;
      c[1] = 5;
      c[2] = 2;
      c[3] = 1;
      return;
    default: // should never happen
      return;
    }
    break;
  case 'F':
  case 'f':
      c[0] = 7;
      c[1] = 2;
      return;
  case 'G':
  case 'g':
      c[0] = 2;
      c[1] = 1;
      return;
  default: // should never happen
    assert(false && "unexpected type in weylSize");
    return;
  }
}

}

}

