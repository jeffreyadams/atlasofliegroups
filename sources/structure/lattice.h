/*
  This is lattice.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2009,2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef LATTICE_H  /* guard against multiple inclusions */
#define LATTICE_H

#include "../Atlas.h"

/******** function declarations **********************************************/

namespace atlas {

namespace lattice {

template<typename I, typename O>
  void baseChange (I, I, O, I, I);

template<typename I, typename O>
  void inverseBaseChange (I, I, O, I, I);

// find basis in dual lattice of perpendicular to given vectors
CoweightList perp (const LatticeMatrix&);

// find matrix whose image is the sublattice annihilated by |M|
LatticeMatrix kernel (LatticeMatrix M); // by value

LatticeMatrix eigen_lattice
  (LatticeMatrix M, // by value
   LatticeCoeff lambda);

// find matrix with same kernel, but with rows spanning saturated sublattice
LatticeMatrix row_saturate (const LatticeMatrix& M);

} // |namespace lattice|

} // |namespace atlas|

#endif
