/*
  This is mod2_system.h: Solving linear systems over the two-element field

  Copyright (C) 2014, Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef MOD2_SYSTEM_H  /* guard against multiple inclusions */
#define MOD2_SYSTEM_H

#include <vector>
#include <iostream>
#include "bitmap.h"

namespace atlas {

namespace mod2_system {

// We solve systems incrementally, feeding unknowns an equations into a system

class Mod2_System
{
  typedef std::vector<unsigned long> linear_combination;
  struct equation {
    linear_combination lhs; unsigned int rhs; // only for its parity
  equation(unsigned int val) : lhs(), rhs(val&1u) {}
    std::ostream& print(std::ostream& f) const;
  };

  std::vector<equation> eqn; // equations in permuted echelon form
  std::vector<unsigned long> pivot_index; // which |eqn| solves each unknown
  enum { no_pivot=~0ul };

 public:
  // the incremental nature of solving means we only have a default constructor
  Mod2_System() : eqn(), pivot_index() {}

 private:
  Mod2_System(const Mod2_System&); // copy forbidden
  Mod2_System& operator=(const Mod2_System&); // assignment forbidden
 public:

  // manipulators

  // add $n$ unknowns, return \emph{old} size (number of first new unknown)
  unsigned long extend(unsigned int n=1);

  // add an equation, with nonzero LHS terms from begin..end
  template <typename I> // input iterator with unsigned value type
    bool add (I begin, I end, unsigned int rhs); // returns whether consistent

  void reduce(); // simplify system to eliminate all pivot unknowns from |eqn|

  // accessors

  unsigned long size() const; // current number of indeterminates
  bool consistent() const;    // whether at least one solution exists

  // the following two methods return |~0ul| in case of an inconsistent system
  unsigned long rank() const; // effective number of equations
  unsigned long dimension() const; // dimension of solution space

  bitmap::BitMap a_solution() const; // sample solution, set of nonzero unknowns

  // the following calls |reduce| for efficiency, hence is a manipulator
  std::vector<bitmap::BitMap> solution_basis(); // kernel of homogenous system

 private:
  // eliminate set pivot variable of |eqn[inx]| from |dst|, return its |rhs|
  unsigned int eliminate (unsigned long inx, bitmap::BitMap& dst) const;
}; // |Mod2_System|

} // |namespace mod2_system|

} // |namespace atlas|

#endif
