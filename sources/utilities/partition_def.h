/*!
\file
  \brief Template definitions for the class Partition.

  This file contains the definitions of the
  templates declared in partition.h:
    - makeOrbits(Partition&, F a, unsigned long, unsigned long) : constructs
      the orbit partition defined by the "action function" a;

The purpose of the class Partition is to compute the partition of a
finite set given by a group action.  A typical example is a Weyl group
acting on elements of order 2 in a torus.
*/
/*
  This is partition_def.h.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include <stack>

#include "bitmap.h"

/*****************************************************************************

  This file contains the definitions of the templates declared in partition.h:

    - makeOrbits(Partition&, F a, unsigned long, unsigned long) : constructs
      the orbit partition defined by the "action function" a;

******************************************************************************/

namespace atlas {

  /*!
  \brief Functions for working with a partition of a finite set into
  classes.

  Typically classes are orbits of a finite group (like a Weyl group)
  on a vector space over Z/2Z (like elements of order 2 in a torus).
  */
namespace partition {


/*! Constructs the orbit partition defined by the "action function" |a|. The
  type |F| is that of some function object that can be given a pair of
  |unsigned long| argument to produce an |unsigned long| result. The parameter
  |a| is such a function object; its first argument is the index $i$ of a
  generator in the range $[0,c[$, its second argument is the value acted upon,
  which is an integer in the range $[0,n[$; the result is again a vlaue in the
  range $[0,n[$. It is assumed that for any fixed $i$, the function
  $a(i,\cdot)$ defines a permutation of $[0,n[$; thus the action parameter |a|
  defines $c$ generators of a permutation subgroup of $Sym([0,n[)$. The orbits
  for this action are returned as a |Partition| object in |pi|.

  The algorithm is straightforward orbit generation, using a set $B$ of
  elements that have not yet joined any orbit, and a collection $S$ of
  elements that are in the current orbit but whose acted-upon images have yet
  to be generated. The set $B$ is represented as a bitmap on $[0,n[$,
  initially full, while $S$ is realized as a stack (but it could equally well
  have been a queue), initially empty.

  While $B$ is not empty :

    - move an element (the first one left) from $B$ to (the currently empty)
      $S$; start a new orbit with this element

    - while $S$ is not empty :
        - x = S.top(); S.pop();
	- for j in [0,c[
            if x'=a(j,x) is in B:
              remove x' from B, push x' onto S, and add it to the current orbit
*/
template<typename F> // class serving as function object (ul,ul)->ul
  void makeOrbits(Partition& pi,const F& a, unsigned long c, unsigned long n)
{
  pi.resize(n);
  pi.clear();

  bitmap::BitMap b(n);
  b.fill();

  std::stack<unsigned long> toDo;
  const bitmap::BitMap::iterator b_end = b.end();

  for (bitmap::BitMap::iterator i = b.begin(); i != b_end; ++i) {

    pi.newClass(*i);
    unsigned long thisClass = pi(*i);
    toDo.push(*i);
    b.remove(*i);

    while (not toDo.empty()) {

      unsigned long x = toDo.top();
      toDo.pop();

      for (unsigned long j = 0; j < c; ++j) {
	unsigned long y = a(j,x);
	if (b.isMember(y)) { // new element was found
	  b.remove(y);
	  pi.addToClass(thisClass,y);
	  toDo.push(y);
	}
      }
    }
  }
}

} // namespace partition

} // namespace atlas
