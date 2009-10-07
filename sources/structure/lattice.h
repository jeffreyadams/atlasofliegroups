/*!
\file
\brief Function declarations for namespace lattice.
*/
/*
  This is lattice.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef LATTICE_H  /* guard against multiple inclusions */
#define LATTICE_H

#include "bitvector_fwd.h"
#include "latticetypes_fwd.h"

/******** function declarations **********************************************/

namespace atlas {

namespace lattice {

template<typename I, typename O>
  void baseChange(I, I, O, I, I);

template<typename I, typename O>
  void inverseBaseChange(I, I, O, I, I);

latticetypes::WeightList perp(const latticetypes::WeightList&, size_t);

} // |namespace lattice|

} // |namespace atlas|

#endif
