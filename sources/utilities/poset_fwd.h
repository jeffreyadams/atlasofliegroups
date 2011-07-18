/*!
\file
  This is poset_fwd.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef POSET_FWD_H  /* guard against multiple inclusions */
#define POSET_FWD_H

#include "set.h"

/******** forward type declarations *****************************************/

namespace atlas {

namespace poset {

  class Poset;
  typedef std::pair<set::Elt,set::Elt> Link;

}

}

#endif
