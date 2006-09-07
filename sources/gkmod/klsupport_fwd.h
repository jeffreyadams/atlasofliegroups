/*!
\file
\brief Forward class declaration and type definitions for KLSupport.
*/
/*
  This is klsupport_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups  

  See file main.cpp for full copyright notice
*/

#ifndef KLSUPPORT_FWD_H  /* guard against multiple inclusions */
#define KLSUPPORT_FWD_H

#include <vector>

#include "polynomials_fwd.h"

namespace atlas {

/******** type declarations *************************************************/

namespace klsupport {

class KLSupport;
class KLSupport_new;

typedef std::vector<size_t> PrimitiveRow;

}

}

#endif
