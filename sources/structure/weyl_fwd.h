/*!
\file
\brief Class declarations and type definitions for WeylGroup.
*/
/*
  This is weyl_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef WEYL_FWD_H  /* guard against multiple inclusions */
#define WEYL_FWD_H

#include <vector>

#include "constants.h"

/******** forward type declarations ******************************************/

namespace atlas {

namespace weyl {

  class WeylElt;
  typedef WeylElt TwistedInvolution;
  class WeylGroup;
  class TwistedWeylGroup;

  typedef std::vector<WeylElt> WeylEltList;
  typedef std::vector<TwistedInvolution> TwistedInvolutionList;

  /*!
\brief Represents a simple root reflection (one of the standard
generators of W).
  */
  typedef unsigned char Generator;
  // make WeylWord genuine member of our namespace
  struct WeylWord : public std::vector<Generator>
  {
  };


}

}

#endif
