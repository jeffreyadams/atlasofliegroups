/*
  This is realredgp_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef REALREDGP_IO_H  /* guard against multiple inclusions */
#define REALREDGP_IO_H

#include <iosfwd>

#include "../Atlas.h"


namespace atlas {

/******** function declarations **********************************************/

namespace realredgp_io {

std::ostream& printBlockStabilizer(std::ostream& strm,
				   RealReductiveGroup& G_R,
				   size_t cn,
				   RealFormNbr rf);

std::ostream& printCartanClasses(std::ostream&,
				 RealFormNbr rf,
				 complexredgp_io::Interface&);

std::ostream& printCartanOrder(std::ostream&,
			       const RealReductiveGroup&);

std::ostream& printRealWeyl(std::ostream& strm,
			    RealReductiveGroup& G_R,
			    size_t cn);

std::ostream& printStrongReal(std::ostream& strm,
			      ComplexReductiveGroup& G_C,
			      const realform_io::Interface& rfi,
			      size_t cn);

} // |namespace realredgp_io|

} // |namespace atlas|


#endif
