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

  For license information see the LICENSE file
*/

#include "weylsize.h"

#include <cassert>

#include "lietype.h"

/*****************************************************************************

        Chapter I -- Functions declared in weylsize.h

******************************************************************************/

namespace atlas {

namespace weylsize {


//! Returns the size of the Weyl group with Lie type lt.
size::Size weylSize(const lietype::LieType& lt)
{
  size::Size c(1);

  for (size_t j=0; j<lt.size(); ++j)
    c *= weylSize(lt[j]);

  return c;
}


/*!
  Synopsis: returns the size of the Weyl group with Lie type slt.

  NOTE: we build in knowledge about the Weyl group here, although we could
  have computed from first principles
*/
size::Size weylSize(const lietype::SimpleLieType& slt)
{
  unsigned long r = slt.rank();

  switch (slt.type())
  {
  case 'A':
    return size::factorial(r+1);
  case 'B':
  case 'C':
    {
      size::Size c = size::factorial(r);
      c[0] += r; // multiply by $2^r$
      return c;
    }
  case 'D':
    {
      size::Size c=size::factorial(r);
      c[0] += r-1; // multiply by $2^{r-1}$
      return c;
    }
  case 'E':
    {
      size::Size c(1); // we'll set components by hand below

      switch (r)
      {
      case 6:
	c[0] = 7;
	c[1] = 4;
	c[2] = 1;
	break; // 51840
      case 7:
	c[0] = 10;
	c[1] = 4;
	c[2] = 1;
	c[3] = 1;
	break; // 2903040
      case 8:
	c[0] = 14;
	c[1] = 5;
	c[2] = 2;
	c[3] = 1;
	break; // 696729600
      }
      return c;
    }
  case 'F':
  case 'f':
    return size::Size(1152); // $2^7 * 3^2$
  case 'G':
  case 'g':
    return size::Size(12); // $2^2 * 3$
  default: // should never happen
    assert(false && "unexpected type in weylSize");
    return size::Size(1);
  }
}

}

}
