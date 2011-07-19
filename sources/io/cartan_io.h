/*
  This is cartan_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef CARTAN_IO_H  /* guard against multiple inclusions */
#define CARTAN_IO_H

#include <iosfwd>

#include "atlas_types.h"

/******** function declarations *********************************************/

namespace atlas {

namespace cartan_io {

std::ostream&
printCartanClass(std::ostream&, size_t, complexredgp_io::Interface&);

std::ostream& printFiber(std::ostream&, const Fiber&,
			 const RealFormNbrList&);

std::ostream& printGradings(std::ostream&, const Fiber&,
			    const RealFormNbrList&,
			    const RootSystem&);

}

}

#endif
