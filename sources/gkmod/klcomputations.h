/*!
\file 
\brief Class definition and function declarations for the class KLComputations.

  [This class has barely been started. It is part of Fokko's plan to
  define an ownership structure for the Kazhdan-Lusztig
  computations. DV 7/23/06.]
*/
/*
  This is klcomputations.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups  

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

  /*!
\brief Holding class for computations related to a block of representations.

  Most prominent among these is the computation of
  Kazhdan-Lusztig-Vogan polynomials. The idea is to mostly apply a
  principle of "lazy evaluation": things are computed when they are
  requested, but always remembered once they are computed.

  We make an exception to the laziness principle for the "involutions" part,
  which is computed right away.

  [This class has barely been started. It is part of Fokko's plan to
  define an ownership structure for the Kazhdan-Lusztig
  computations. DV 7/23/06.]
  */
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
