/*!
\file
\brief Function declarations for namespace weylsize.

Functions to compute the size of a Weyl group (stored as a prime
factorization, to avoid overflow problems).
*/
/*
  This is weylsize.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef WEYLSIZE_H  /* guard against multiple inclusions */
#define WEYLSIZE_H

#include "atlas_types.h"

#include "size_fwd.h"

/******** function declarations **********************************************/

namespace atlas {

namespace weylsize {

  size::Size weylSize(const LieType&);
  size::Size weylSize(const SimpleLieType&);

}

}

#endif
