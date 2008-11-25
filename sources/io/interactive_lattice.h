/*
  This is interactive_lattice.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef INTERACTIVE_LATTICE_H  /* guard against multiple inclusions */
#define INTERACTIVE_LATTICE_H

#include "latticetypes_fwd.h"

#include "error.h"
#include "lietype.h"

/******** function declarations **********************************************/

namespace atlas {

namespace interactive_lattice {

  void adjustBasis(latticetypes::WeightList&, latticetypes::CoeffList&,
		   const latticetypes::WeightList&,
		   const latticetypes::CoeffList&);

  int getGenerators(latticetypes::RatWeightList&,
		    const latticetypes::CoeffList&)
    throw(error::InputError);

  int getLattice(latticetypes::CoeffList&, latticetypes::WeightList&)
    throw(error::InputError);

  void getUniversal(latticetypes::CoeffList&, const latticetypes::CoeffList&);

  void localBasis(latticetypes::WeightList&, const latticetypes::WeightList&,
		  const latticetypes::CoeffList&);

  void smithBasis(latticetypes::CoeffList&, latticetypes::WeightList&,
		  const lietype::LieType&);

}

}

#endif
