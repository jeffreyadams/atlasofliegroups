/*
  This is prerootdata.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef PREROOTDATA_H  /* guard against multiple inclusions */
#define PREROOTDATA_H

#include "latticetypes.h"
#include "lietype.h"

namespace atlas {

/******** type declarations *************************************************/

namespace prerootdata {

  class PreRootDatum;

}

/******** function declarations **********************************************/


namespace prerootdata {

  void cartanMatrix(latticetypes::LatticeMatrix&, const lietype::LieType&);
  void cartanMatrix(latticetypes::LatticeMatrix&, 
		    const lietype::SimpleLieType&);

}

/******** type definitions **************************************************/

namespace prerootdata {

class PreRootDatum{

 private:

  latticetypes::WeightList d_roots;
  latticetypes::WeightList d_coroots;
  size_t d_rank;

 public:

// constructors and destructors
  PreRootDatum() {}

  PreRootDatum(const lietype::LieType& lt, const latticetypes::WeightList&);

  ~PreRootDatum() {}

// accessors
  bool operator== (const PreRootDatum& prd) const {
    return (d_roots == prd.d_roots) and (d_coroots == prd.d_coroots) and
       (d_rank == prd.d_rank);
  }

  const latticetypes::WeightList& coroots() const {
    return d_coroots;
  }

  size_t rank() const {
    return d_rank;
  }

  const latticetypes::WeightList& roots() const {
    return d_roots;
  }

// manipulators
  void swap(PreRootDatum&);
};

}

}

#endif
