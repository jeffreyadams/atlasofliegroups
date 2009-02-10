/*!
\file
\brief Template definitions for the class BitSet.
*/
/*
  This is bitset_def.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

/*****************************************************************************

        Chapter I -- Template functions declared in bitset.h

******************************************************************************/

#include <cassert>

namespace atlas {

namespace bitset {

template <size_t n>
BitSet<n> operator~ (const BitSet<n>& d_b) {
  BitSet<n> b = d_b;
  return b.flip();
}


/*
  Sets the first d bits of v.
*/
template<size_t n> void set(BitSet<n>& v, size_t d)
{
  v.set();
  v.truncate(d);
}

template <size_t n>
  template<typename I>
  BitSet<n>::BitSet(const std::vector<I>& v) : Base()
{
  assert(v.size()<=n);
  for (size_t i=0; i<v.size(); ++i)
    set(i,v[i]%2!=0);
}

} // |namespace bitset|

} // |namespace atlas|
