/*
  This is klcomputations.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef KLCOMPUTATIONS_H  /* guard against multiple inclusions */
#define KLCOMPUTATIONS_H

#include "klcomputations_fwd.h"

#include "blocks_fwd.h"
#include "complexredgp_fwd.h"
#include "kgb_fwd.h"

#include "involutions.h"

namespace atlas {

/******** type declarations *************************************************/

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace klcomputations {

class KLComputations {

 private:

  involutions::InvolutionSet d_involutionSet;
  std::vector<kgb::KGB> d_kgb;
  std::vector<kgb::KGB> d_dualkgb;
  std::vector<blocks::Block> d_block;
  std::vector<blocks::Block> d_dualblock;

  complexredgp::ComplexReductiveGroup* d_G;

 public:

  KLComputations();
  
  KLComputations(complexredgp::ComplexReductiveGroup&);  // un-owned pointer
};

}

}

#endif
