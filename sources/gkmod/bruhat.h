/*
  This is bruhat.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/
// Class definition and function declarations for the class |BruhatOrder|

#ifndef BRUHAT_H  /* guard against multiple inclusions */
#define BRUHAT_H

#include "poset.h"	// containment

namespace atlas {

/******** type declarations *************************************************/

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace bruhat {

/*
  Intended to represent the Bruhat order on K orbits on G/B, or
  on a block of representations.

  In fact just stores any given relation as the Hasse diagram, and is capable
  of expanding the full poset by transitivity.
*/
class BruhatOrder
{
/*
  Hasse diagram for a Bruhat order. Entry \#j lists the numbers of the
  immediate predecessors of element \#j in the order.
*/
  // As this is probably sparse; avoid |BitMap|s
  std::vector<poset::Poset::EltList> d_hasse;

/*
   Full poset relation. Only computed on demand, until then has size 0.

   It is assumed (by the class |Poset|) that element |i| can precede
   element |j| in the poset only if |i < j|.
*/
  poset::Poset d_poset;

 public:

// constructors and destructors
  explicit BruhatOrder(const std::vector<poset::Poset::EltList>& Hasse_diagram)
    : d_hasse(Hasse_diagram), d_poset(0) {}

  explicit BruhatOrder(const std::vector<poset::Poset::EltList>&& Hasse_diagram)
    : d_hasse(std::move(Hasse_diagram)), d_poset(0) {}


// accessors

  size_t size() const { return d_hasse.size(); }

  // Return row |x| of the Hasse diagram for the order.
  const std::vector<poset::Poset::Elt>& hasse(size_t x) const
  { return d_hasse[x]; }

  // Return the number of comparable pairs in the order.
  unsigned long n_comparable() const
  { return poset::n_comparable_from_Hasse(d_hasse); }

  // manipulators
  // Return the full poset relation.
  const poset::Poset& poset() { fillPoset(); return d_poset; }

  private:
  void fillPoset();

}; // |class BruhatOrder|

} // |namespace bruhat|

} // |namespace atlas|

#endif
