/*
  This is kl_io.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef KL_IO_H  /* guard against multiple inclusions */
#define KL_IO_H

#include <iosfwd>

#include "kl_fwd.h"

namespace atlas {

/******** function declarations *********************************************/

namespace kl_io {

  std::ostream& printAllKL(std::ostream&, kl::KLContext&);

  std::ostream& printExtremalKL(std::ostream&, const kl::KLContext&);

  std::ostream& printKLList(std::ostream&, kl::KLContext&);

  std::ostream& printMu(std::ostream&, const kl::KLContext&);
}

}

#endif
