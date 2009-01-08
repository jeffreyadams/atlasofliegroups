/*
  This is poset_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef POSET_IO_H  /* guard against multiple inclusions */
#define POSET_IO_H

#include <iosfwd>

#include "poset_fwd.h"

namespace atlas {

/******** function declarations *********************************************/

namespace poset_io {

  std::ostream& printPoset(std::ostream&, const poset::Poset&);

}

}

#endif
