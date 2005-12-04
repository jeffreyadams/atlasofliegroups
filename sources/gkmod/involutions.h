/*
  This is involutions.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef INVOLUTIONS_H  /* guard against multiple inclusions */
#define INVOLUTIONS_H

#include "involutions_fwd.h"

#include "complexredgp_fwd.h"

namespace atlas {

/******** type declarations *************************************************/

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace involutions {

class InvolutionSet {

 private:

 public:
  // constructors and destructors
  InvolutionSet();

  InvolutionSet(complexredgp::ComplexReductiveGroup&);

};

}

}

#endif
