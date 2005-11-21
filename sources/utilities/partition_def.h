/*
  This is partition_def.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include <stack>

#include "bitmap.h"

/*****************************************************************************

  This file contains the definitions of the templates declared in partition.h:

    - makeOrbits(Partition&, F a, unsigned long, unsigned long) : constructs
      the orbit partition defined by the "action function" a;

******************************************************************************/

namespace atlas {

namespace partition {

template<typename F> 
  void makeOrbits(Partition& pi, F& a, unsigned long c, unsigned long n)

/*
  It is assumed that a is a function object which takes two unsigned long
  arguments, one in the range [0,c[, one in the range [0,n[. This is
  interpreted as an action of the set [0,c[ on [0,n[; moreover we assume
  that these actions are given by permutations (so in effect we are giving
  the generators of a permutation subgroup of Sym([0,n[).)

  Then we write in pi the orbit partition for this action. The algorithm is
  like this. Maintain a stack S of elements-under-inspection, initially
  empty, and a bitmap B of not-yet-considered elements, originally [0,n[.

  While B is not empty :

    - delete the first element from B, and push it on S; start new class
      with that element;

    - while S is not empty :
        - x = S.top(); S.pop();
	- for j in [0,c[ : if a(j,x) is in B, delete it from B, add it to 
	  the current orbit, and push it on S
*/

{
  using namespace bitmap;

  pi.resize(n);
  pi.clear();

  BitMap b(n);
  b.fill();

  std::stack<unsigned long> toDo;
  const BitMap::iterator b_end = b.end();

  for (BitMap::iterator i = b.begin(); i != b_end; ++i) {

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

  return;
}

}

}
