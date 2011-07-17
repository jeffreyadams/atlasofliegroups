 /*
  This is kgb_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef KGB_IO_H  /* guard against multiple inclusions */
#define KGB_IO_H

#include <iosfwd>

#include "kgb.h"

namespace atlas {

/******** function declarations *********************************************/

namespace kgb_io {

  std::ostream& print(std::ostream& strm,
		      const KGB_base& kgb,
		      bool traditional =false,
		      const ComplexReductiveGroup* G=NULL,
		      const KGBEltList* which = NULL);

  std::ostream& printKGB(std::ostream&, const KGB&);
  std::ostream& var_print_KGB(std::ostream&,
			      const ComplexReductiveGroup&,
			      const KGB&);

  std::ostream& print_sub_KGB(std::ostream& strm,
			      const KGB& kgb,
			      const KGBEltList& which);
  std::ostream& print_X(std::ostream&, const kgb::global_KGB&);

  std::ostream& printBruhatOrder(std::ostream&, const BruhatOrder&);

}

}

#endif
