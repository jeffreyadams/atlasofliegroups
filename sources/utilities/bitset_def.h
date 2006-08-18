/*!
\file
\brief Template definitions for the class BitSet.
*/
/*
  This is bitset_def.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

/*****************************************************************************

        Chapter I -- Template functions declared in bitset.h

******************************************************************************/

namespace atlas {

namespace bitset {

template <size_t n>
BitSet<n> operator~ (const BitSet<n>& d_b) {
  BitSet<n> b = d_b;
  return b.flip();
}

template<size_t n> void set(BitSet<n>& v, size_t d)

/*
  Sets the first d bits of v.
*/

{
  v.set();
  v.truncate(d);

  return;
}

}

}
