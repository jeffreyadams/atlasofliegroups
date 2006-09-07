/*!
\file
\brief Implementation of the BitSet class.
*/
/*  
  This is bitset.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "bitset.h"

/*****************************************************************************

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace bitset {

BitSetBase<1>::iterator BitSetBase<1>::begin() const

/*!
  Synopsis : returns an iterator pointing to the first set bit; this is
  essentially just d_bits itself.
*/

{
  return BitSetBase<1>::iterator(d_bits);
}

bool BitSetBase<1>::scalarProduct(const BitSetBase<1>& b) const

/*!
  Synopsis: returns the parity of *this & b.
*/

{
  unsigned long a = d_bits;
  a &= b.d_bits;

  return bits::bitCount(a) & 1ul;
}

BitSetBase<1>& BitSetBase<1>::slice(const BitSetBase<1>& c)

/*
  Synopsis : replaces the bitset by the "defragmented" intersection with
  c (i.e., the bits of that intersection are written consecutively.)
*/

{
  operator&= (c);
  size_t d = 0;

  for (BitSetBase<1>::iterator i = c.begin(); i(); ++i) {
    set(d,test(*i));
    ++d;
  }

  truncate(d);

  return *this;
}

}

}
