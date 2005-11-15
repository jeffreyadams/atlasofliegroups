/*
  This is kgb_io.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef KGB_IO_H  /* guard against multiple inclusions */
#define KGB_IO_H

#include <iosfwd>

#include "kgb_fwd.h"

namespace atlas {

/******** function declarations *********************************************/

namespace kgb_io {

  std::ostream& printKGB(std::ostream&, const kgb::KGB&);

}

}

#endif
