/*
  This is bruhat.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef BRUHAT_H  /* guard against multiple inclusions */
#define BRUHAT_H

#include "bruhat_fwd.h"

#include "poset.h"

namespace atlas {

/******** type declarations *************************************************/

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace bruhat {

class BruhatOrder {

 private:

  std::vector<set::SetEltList> d_hasse;
  poset::SymmetricPoset d_poset;

 public:

// constructors and destructors
  explicit BruhatOrder(const std::vector<set::SetEltList>& hd);

// accessors
  const set::SetEltList& hasse(size_t x) const {
    return d_hasse[x];
  }

  const poset::SymmetricPoset& poset() const {
    return d_poset;
  }

// manipulators
  set::SetEltList& hasse(size_t x) {
    return d_hasse[x];
  }


};

}

}

#endif
