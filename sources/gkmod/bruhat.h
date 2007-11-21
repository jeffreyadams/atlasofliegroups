/*!
\file
\brief Class definition and function declarations for the class BruhatOrder.
*/

/*
  This is bruhat.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
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

/*!\brief
  Intended to represent the Bruhat order on K orbits on G/B, or
  on a block of representations.

  [Not seriously used in the present code; I'm not sure whether it is
  instantiated.  The classes KGB and Block have manipulators fillBruhat
  that would create BruhatOrder classes, but they seem not to be called.
  DV 7/21/06]

  [In fact the kgb::KGB and blocks::Block classes do contain members of type
  (pointer to) BruhatOrder, and methods fillBruhat that will create instances
  of BruhatOrder accessible from those classes, but these methods are never
  called. MvL 8 Nov 2006]

  [Now instantiated by the new commands blockorder and kgborder.]
*/
class BruhatOrder {

 private:
  /*!
\brief Hasse diagram for a Bruhat order.

Entry \#j lists the numbers of the immediate predecessors of element
  \#j in the order.
  */
  std::vector<set::SetEltList> d_hasse; // probably sparse; avoid |BitMap|s
  /*!
\brief Poset relation.

It is assumed that element \#i can precede element \#j in the poset
only if i < j.
  */
  poset::Poset d_poset;

 public:

// constructors and destructors
  explicit BruhatOrder(const std::vector<set::SetEltList>& Hasse_diagram)
    : d_hasse(Hasse_diagram), d_poset(0) {}


// accessors

  size_t size() const { return d_hasse.size(); }

  //!\brief Returns row |x| of the Hasse diagram for the order.

  const set::SetEltList& hasse(size_t x) const {
    return d_hasse[x];
  }

  /*!
\brief Returns the number of comparable pairs in the order.
  */
  unsigned long n_comparable() const {
    return poset::n_comparable_from_Hasse(d_hasse);
  }

  // manipulators
  /*!
\brief Returns the full poset relation.
   */
  const poset::Poset& poset() {
    fillPoset(); return d_poset;
  }

  private:
  void fillPoset();

}; // class BruhatOrder

} // namespace bruhat

} // namespace atlas

#endif
