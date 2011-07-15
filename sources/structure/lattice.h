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

#include "atlas_types.h"

/******** function declarations **********************************************/

namespace atlas {

namespace lattice {

template<typename I, typename O>
  void baseChange(I, I, O, I, I);

template<typename I, typename O>
  void inverseBaseChange(I, I, O, I, I);

// find basis in dual lattice of perpendicular to given vectors
CoweightList perp(const WeightList&, size_t);

// find matrix whose image is the sublattice annihilated by |M|
LatticeMatrix kernel(const LatticeMatrix& M);

LatticeMatrix eigen_lattice
  (LatticeMatrix M, // by value
   LatticeCoeff lambda);

// find matrix with same kernel, but with rows spanning saturated sublattice
LatticeMatrix row_saturate(const LatticeMatrix& M);

} // |namespace lattice|

} // |namespace atlas|

#endif
