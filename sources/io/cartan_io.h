/*
  This is cartan_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef CARTAN_IO_H  /* guard against multiple inclusions */
#define CARTAN_IO_H

#include <iosfwd>

#include "complexredgp_io_fwd.h"
#include "cartanclass_fwd.h"
#include "rootdata_fwd.h"

#include "realform.h"

/******** function declarations *********************************************/

namespace atlas {

namespace cartan_io {

std::ostream&
printCartanClass(std::ostream&, size_t, const complexredgp_io::Interface&);

std::ostream& printFiber(std::ostream&, const cartanclass::Fiber&,
			 const realform::RealFormList&);

std::ostream& printGradings(std::ostream&, const cartanclass::Fiber&,
			    const realform::RealFormList&,
			    const rootdata::RootDatum&);

}

}

#endif
