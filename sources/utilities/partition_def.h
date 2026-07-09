/*
  Template definitions for the class |Partition|.

  This file contains the definition of the template declared in partition.h:

    - orbits(Partition&, F a, unsigned long, unsigned long) : constructs
      the orbit partition defined by the "action function" a;

  This is partition_def.h.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/


#include "bitmap.h"
#include "sl_list.h"

/*****************************************************************************

  This file contains the definition of the template declared in partition.h:

    - Partition orbits(F a, unsigned long n, unsigned long limit) : constructs
      the orbit partition defined by the "action function" a;

******************************************************************************/

namespace atlas {

namespace partition {

/* Although the function |orbits| currently needs only one instantiation (from
   cartanclass.cpp), it has to be defined as a template, and has to be
   dynamically instantiated: explicit instantiation here is impossible because
   the template argument |F| needs to be specialised to a class local to the
   instantiating code (in the anonymous namespace of the file cartanclass.cpp)
*/

/*! Constructs the orbit partition defined by the "action function" |a|. The
  type |F| is that of some function object that can be given a pair of
  |unsigned long| argument to produce an |unsigned long| result. The parameter
  |a| is such a function object; its first argument is the index $i$ of a
  generator in the range $[0,c[$, its second argument is the value acted upon,
  which is an integer in the range $[0,n[$; the result is again a value in the
  range $[0,n[$. It is assumed that for any fixed $i$, the function
  $a(i,\cdot)$ defines a permutation of $[0,n[$; thus the action parameter |a|
  defines $c$ generators of a permutation subgroup of $Sym([0,n[)$. The orbits
  for this action are returned as a |Partition| object in |pi|.

  The algorithm is straightforward orbit generation, using a set $B$ of
  elements that have not yet joined any orbit, and a collection $S$ of
  elements that are in the current orbit but whose acted-upon images have yet
  to be generated. The set $B$ is represented as a bitmap on $[0,n[$,
  initially full, while $S$ is realized as a queue, initially empty.

  While $B$ is not empty :

    - move an element (the first one left) from $B$ to (the currently empty)
      $S$; start a new orbit with this element

    - while $S$ is not empty :
        - x = S.top(); S.pop();
	- for j in [0,c[
            if x'=a(j,x) is in B:
              remove x' from B, push x' onto S, and add it to the current orbit
*/
template<typename F> // class |F| holds a binary function object: |(ui,ul)->ul|
  Partition orbits(const F& a, unsigned int c, unsigned long n)
{
  Partition result(n);
  std::vector<unsigned long> to_do; to_do.reserve(n);
  bitmap::BitMap b(n); // as yet unclassified elements
  b.fill();            // initially that means everyone

  for (bitmap::BitMap::iterator it=b.begin(); it(); ++it) // |b| will shrink
  {
    unsigned long root = *it; // starting element for a fresh orbit
    unsigned long thisClass = result.new_class(root);
    to_do.clear(); to_do.push_back(root);
    b.remove(root); // avoid looping back to |root| later

    do
    {
      unsigned long x = to_do.back(); to_do.pop_back();

      for (unsigned int i=0; i<c; ++i)
      {
	unsigned long y = a(i,x);
	if (b.isMember(y)) // new element was found
	{
	  b.remove(y);
	  result.addToClass(thisClass,y);
	  to_do.push_back(y);
	}
      }
    }
    while (not to_do.empty());
  }
  return result;
}

} // |namespace partition|

} // |namespace atlas|
