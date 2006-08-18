/*!
\file
\brief Class definition and function declarations for the class BruhatOrder.
*/

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

  /*!
\brief Intended to represent the Bruhat order on K orbits on G/B, or
on a block of representations.  

[Not seriously used in the present code; I'm not sure whether it is
instantiated.  The classes KGB and Block have manipulators fillBruhat
that would create BruhatOrder classes, but they seem not to be called. 
DV 7/21/06] 
  */ 
  class BruhatOrder {

 private:
  /*!
\brief Hasse diagram for a Bruhat order.

Entry \#j lists the numbers of the immediate predecessors of element
  \#j in the order.
  */
  std::vector<set::SetEltList> d_hasse;
  /*!
\brief Symmetrization of the poset relation.

It is assumed that element \#i can precede element \#j in the poset
only if i < j, so the symmetrized relation determines the relation.
  */  
  poset::SymmetricPoset d_poset;

 public:

// constructors and destructors
  explicit BruhatOrder(const std::vector<set::SetEltList>& hd);

// accessors

/*!
\brief Returns the Hasse diagram for the order.
*/
  const set::SetEltList& hasse(size_t x) const {
    return d_hasse[x];
  }

  /*!
\brief Returns the symmetrization of the poset relation.
   */
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
