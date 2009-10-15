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

// find basis in dual lattice of perpendicular to given vectors
latticetypes::WeightList perp(const latticetypes::WeightList&, size_t);

// find matrix whose image is the sublattice annihilated by |M|
latticetypes::LatticeMatrix kernel(const latticetypes::LatticeMatrix& M);

latticetypes::LatticeMatrix eigen_lattice
  (latticetypes::LatticeMatrix M, // by value
   latticetypes::LatticeCoeff lambda);

} // |namespace lattice|

} // |namespace atlas|

#endif
