/*!
\file
\brief Implementation of the KLComputations class.

  [This class has barely been started. It is part of Fokko's plan to
  define an ownership structure for the Kazhdan-Lusztig
  computations. DV 7/23/06.]
*/
/*
  This is klcomputations.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "klcomputations.h"

#include "blocks.h"
#include "kgb.h"

namespace atlas {

/*****************************************************************************

  The KLComputations class is a holding class for the various computations
  related to blocks of representations, and most prominently the computation
  of Kazhdan-Lusztig-Vogan polynomials. The idea is to mostly apply a
  principle of "lazy evaluation": things are computed when they are requested,
  but always remembered once they are computed.

  We make an exception to the lazyness principle for the "involutions" part,
  which is computed right away.

******************************************************************************/

/*****************************************************************************

        Chapter I -- The KLComputations class

  ... explain here when it is stable ...

******************************************************************************/

namespace klcomputations {

KLComputations::KLComputations()
  :d_G(0)

{}

KLComputations::KLComputations(complexredgp::ComplexReductiveGroup& G)
  :d_involutionSet(G),d_G(&G)

/*!
\brief Constructor.

  The involutionSet part is filled right away to make things simpler.
*/

{}

}

}
