/*
  This is poset_io.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef POSET_IO_H  /* guard against multiple inclusions */
#define POSET_IO_H

#include <iosfwd>

#include "poset_fwd.h"

namespace atlas {

/******** function declarations *********************************************/

namespace poset_io {

  std::ostream& printPoset(std::ostream&, const poset::Poset&);

  std::ostream& printSymmetricPoset(std::ostream&, 
				    const poset::SymmetricPoset&);

}

}

#endif
