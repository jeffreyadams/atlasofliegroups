/*!
\file
\brief  Declarations for the class Status and the type Grading.
*/
/*
  This is gradings_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef GRADINGS_FWD_H  /* guard against multiple inclusions */
#define GRADINGS_FWD_H

#include <vector>

#include "bitset_fwd.h"

namespace atlas {

/******** type declarations *************************************************/

namespace gradings {

  class Status;
  typedef std::vector<Status> StatusList;

  typedef bitset::RankFlags Grading;
  typedef std::vector<Grading> GradingList;

  class GradingCompare;

}

}

#endif
