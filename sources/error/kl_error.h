/*
  This is kl_error.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef KL_ERROR_H  /* guard against multiple inclusions */
#define KL_ERROR_H

#include "kl_fwd.h"

/******** type declarations *************************************************/

namespace atlas {

namespace kl_error {

  struct KLError;

}

/******** type definitions **************************************************/

namespace kl_error {

struct KLError {
  // data
  size_t x;
  size_t y;
  int line;
  const kl::KLContext* klc;
  // constructors
  KLError(size_t x, size_t y, int line, const kl::KLContext& klc)
    :x(x), y(y), line(line), klc(&klc) {}
  // accessors
  void operator() (const char*) const;
};

}

}

#endif
