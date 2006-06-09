/*!
\file
  This is weyl_fwd.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef WEYL_FWD_H  /* guard against multiple inclusions */
#define WEYL_FWD_H

#include <vector>

#include "constants.h"

/******** forward type declarations ******************************************/

namespace atlas {

namespace weyl {

  class WeylElt;
  /*!
  The Weyl group class is a variation on my favourite implementation in terms
  of transducers.

  I have tried to make a careful choice of datatype for the group elements in
  order to strike the right balance between efficiency and generality. This has
  led me to the choice of _fixed size_ arrays of unsigned chars, representing
  the "parabolic subquotient" representation of a given element; in other
  words, in a group of rank n, the first n elements of the array are used
  (but since it is fixed size, we have to allocate it to some constant value,
  namely RANK_MAX.) This is not so wasteful as it may sound : of course the
  representation as a single number is more compact, but will overflow even
  on 64-bit machines for some groups of rank <= 16, and much sooner of course
  on 32-bit machines; also it imposes some computational overhead for packing
  and unpacking. Any variable-size structure like the STL vector uses already
  three unsigned longs for its control structure (the address of the data, the
  size and the capacity), and then it still has to allocate. This could
  perhaps be simplified to just a pointer (after all the size of the
  allocation is known to the group) but you still have the big overhead of
  allocating and deallocating memory from the heap, and remembering to delete
  the pointers when they go out of scope, or else use autopointers ...

  And if worst comes to worst and one really is worried about a factor 2
  waste for type E8 (which might be significant if one deals with huge
  containers of group elements), then one can still recompile with MAX_RANK=8,
  which will then give a datatype that in 64 bits is the same size as an
  unsigned long.

  Notice that the unsigned char type miraculously suffices for all subquotients
  of all groups up to rank 128 (indeed, the biggest subquotient for B128 is
  of order 256), _provided_ the generators are enumerated in an appropriate
  order. This forces us to do quite a bit of type recognition, which is
  relegated to the dynkin namespace. Because of this reordering, the group
  carries a little interface that will translate back and forth from the
  external ordering and the internal one.
*/
  class WeylGroup;

  typedef std::vector<WeylElt> WeylEltList;

  typedef unsigned char Generator;
  typedef unsigned char EltPiece;

  typedef std::vector<Generator> WeylWord;
  typedef Generator Twist[constants::RANK_MAX];

}

}

#endif
