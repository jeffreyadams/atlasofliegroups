 /*
  This is kgb_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef KGB_IO_H  /* guard against multiple inclusions */
#define KGB_IO_H

#include <iosfwd>

#include "../Atlas.h"

namespace atlas {

/******** function declarations *********************************************/

namespace kgb_io {

  std::ostream& print(std::ostream& strm,
		      const KGB_base& kgb,
		      bool traditional =false,
		      const InnerClass* G=NULL,
		      const KGBEltList* which = NULL);

  std::ostream& printKGB(std::ostream&, const KGB&);
  std::ostream& var_print_KGB(std::ostream&,
			      const InnerClass&,
			      const KGB&);

  std::ostream& print_sub_KGB(std::ostream& strm,
			      const KGB& kgb,
			      const KGBEltList& which);
  std::ostream& print_X(std::ostream&, const global_KGB&);

  std::ostream& printBruhatOrder(std::ostream&, const BruhatOrder&);

  // make a '.dot' file that can be processed by the 'dot' program
  // see www.graphviz.org for more info
  void makeDotFile(std::ostream& strm, const KGB& kgb, const BruhatOrder& bruhat);
}

}

#endif
