 /*
  This is kgb_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef KGB_IO_H  /* guard against multiple inclusions */
#define KGB_IO_H

#include <iosfwd>

#include "kgb_fwd.h"
#include "bruhat_fwd.h"
#include "complexredgp_fwd.h"

namespace atlas {

/******** function declarations *********************************************/

namespace kgb_io {

  std::ostream& printKGB(std::ostream&, const kgb::KGB&);
  std::ostream& var_print_KGB(std::ostream&,
			      const complexredgp::ComplexReductiveGroup&,
			      const kgb::KGB&);

  std::ostream& print_sub_KGB(std::ostream& strm,
			      const kgb::KGB& kgb,
			      const kgb::KGBEltList& which);
  std::ostream& printBruhatOrder(std::ostream&, const bruhat::BruhatOrder&);

}

}

#endif
