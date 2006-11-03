/*!
\file
  This is bits.cpp
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "bits.h"

#include "constants.h"

/*****************************************************************************

  ... explain here when it is stable ...

******************************************************************************/

/*****************************************************************************

        Chapter I -- Functions declared in bits.h

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace bits {

unsigned bitCount(unsigned long f)

/*!
  Synopsis: returns the number of set bits in f.
*/

{
  unsigned count = 0;

  for (; f; f &= f-1)
    ++count;	/* see K&R */

  return count;
}

size_t firstBit(unsigned long f)

/*!
  Synopsis: returns the position of the first set bit in f.

  Returns longBits if there is no such bit.
*/

{
  using namespace constants;

  if (f == 0)
    return longBits;

  size_t fb = 0;

  for (; (f & firstChar) == 0; f >>= charBits)
    fb += charBits;

  return fb + firstbit[f & firstChar];
}

size_t lastBit(unsigned long f)

/*!
  Synopsis: returns the position of the first [should be last - DV]
  set bit in f, plus one.  Returns 0 if there is no such bit.
*/

{
  using namespace constants;

  if (f == 0)
    return 0;

  unsigned lb = 0;

  for (; f & ~firstChar; f >>= charBits)
    lb += charBits;

  return lb + lastbit[f];
}

void permute(unsigned long& f, const setutils::Permutation& a)

/*!
  Synopsis: permutes the bits of f according to a: if the bits of f are
  interpreted with respect to a sequence of objects x_0,...,x_{n-1}, then the
  bits of the result are interpreted with respect to the permuted sequence
  x_{a[0]},...,x_{a[n-1]}

  Precisely, we transform f = f_old into f_new, where f_new[i] = f_old[a[i]].

  Note that the permutation of bits is not interpreted in the same way as the
  permutation of vector entries in setutils::permute, given in setutils_def.h
*/

{
  using namespace constants;

  unsigned long f_new = 0;

  for (size_t j = 0; j < a.size(); ++j)
    if (f & bitMask[a[j]])
      f_new |= bitMask[j];

  f = f_new;

  return;
}

}

}
