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

template<typename I, typename J, typename O>
  void baseChange(const I&, const I&, O, const J&, const J&);

template<typename I, typename J, typename O>
  void inverseBaseChange(const I&, const I&, O, const J&, const J&);

template<size_t dim>
  void mod2(bitvector::BitVector<dim>&, const latticetypes::LatticeElt&);

template<size_t dim>
  void mod2(std::vector<bitvector::BitVector<dim> >&,
	    const latticetypes::WeightList&);

template<size_t dim> void mod2(bitvector::BitMatrix<dim>&,
			       const latticetypes::LatticeMatrix&);

latticetypes::LatticeMatrix numeratorMatrix(const latticetypes::RatWeightList&);

latticetypes::WeightList perp(const latticetypes::WeightList&, size_t);

latticetypes::RatWeightList toCommonDenominator
(const latticetypes::RatWeightList& l);

}

}

#include "lattice_def.h"

#endif
