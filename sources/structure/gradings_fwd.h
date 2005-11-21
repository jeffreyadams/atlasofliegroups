/*
  This is gradings_fwd.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
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
  class FullStatus;
  typedef std::vector<FullStatus> FullStatusList;

  typedef bitset::RankFlags Grading;
  typedef std::vector<Grading> GradingList;

  class GradingCompare;

}

}

#endif
