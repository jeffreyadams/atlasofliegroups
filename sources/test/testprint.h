/*
  This is testprint.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef TESTPRINT_H  /* guard against multiple inclusions */
#define TESTPRINT_H

#include <iostream>

#include "atlas_types.h"

#include "cartanclass_fwd.h"
#include "complexredgp_fwd.h"
#include "complexredgp_io_fwd.h"
#include "realredgp_fwd.h"

#include "gradings.h"

/******** function declarations **********************************************/

namespace atlas {

namespace testprint {

std::ostream& print(std::ostream&, const RootDatum&);

std::ostream& print(std::ostream&, const realredgp::RealReductiveGroup&);

std::ostream& printBlockData(std::ostream&, complexredgp_io::Interface&);

std::ostream& printCartanClasses(std::ostream&, complexredgp_io::Interface&);

std::ostream& printCartanMatrix(std::ostream&, const RootNbrList&,
				const RootSystem&);

std::ostream& printComponents(std::ostream&,
			      const realredgp::RealReductiveGroup&,
			      const char* sep = ",");
}

}

#endif
