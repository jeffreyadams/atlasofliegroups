/*
  This is lattice.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef LATTICE_H  /* guard against multiple inclusions */
#define LATTICE_H

#include "bitvector_fwd.h"
#include "latticetypes_fwd.h"

/******** function declarations **********************************************/

namespace atlas {

namespace lattice {

template<typename J>
  void baseChange(latticetypes::LatticeMatrix&, const J&, const J&);

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

void numeratorMatrix(latticetypes::LatticeMatrix&, 
		     const latticetypes::RatWeightList&);

void perp(latticetypes::WeightList&, const latticetypes::WeightList&, size_t);

void toCommonDenominator(latticetypes::RatWeightList&, 
			 const latticetypes::RatWeightList&);

}

}

#include "lattice_def.h"

#endif
