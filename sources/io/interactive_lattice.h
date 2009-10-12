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

  int getGenerators(latticetypes::RatWeightList&,
		    const latticetypes::CoeffList&)
    throw(error::InputError);

  int getLattice(const latticetypes::CoeffList&, latticetypes::WeightList&)
    throw(error::InputError);

}

}

#endif
