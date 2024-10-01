/*
  This is realweyl_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef REALWEYL_IO_H  /* guard against multiple inclusions */
#define REALWEYL_IO_H


#include <iosfwd>
#include "../Atlas.h"


/******** function declarations *********************************************/

namespace atlas {

namespace realweyl_io {

  std::ostream& printBlockStabilizer(std::ostream&, const realweyl::RealWeyl&,
				     const realweyl::RealWeylGenerators&);

  std::ostream& printDualRealWeyl(std::ostream&, const realweyl::RealWeyl&,
				  const realweyl::RealWeylGenerators&);

  std::ostream& printRealWeyl(std::ostream&, const realweyl::RealWeyl&,
			      const realweyl::RealWeylGenerators&);

}

}

#endif
