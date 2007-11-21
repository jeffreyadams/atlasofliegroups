/*
  This is realweyl_io.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  For license information see the LICENSE file
*/

#ifndef REALWEYL_IO_H  /* guard against multiple inclusions */
#define REALWEYL_IO_H

#include <iosfwd>

#include "realweyl_fwd.h"
#include "rootdata_fwd.h"
#include "weyl_fwd.h"

/******** function declarations *********************************************/

namespace atlas {

namespace realweyl_io {

  std::ostream& printBlockStabilizer(std::ostream&, const realweyl::RealWeyl&,
				     const realweyl::RealWeylGenerators&,
				     const rootdata::RootDatum& rd);

  std::ostream& printDualRealWeyl(std::ostream&, const realweyl::RealWeyl&,
				  const realweyl::RealWeylGenerators&,
				  const rootdata::RootDatum& rd);

  std::ostream& printRealWeyl(std::ostream&, const realweyl::RealWeyl&,
			      const realweyl::RealWeylGenerators&,
			      const rootdata::RootDatum& rd);

}

}

#endif
