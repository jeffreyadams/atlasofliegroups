/*!
\file
  This is intutils_def.h. This file contains some very simple
  templates for integer manipulations. 
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

/*****************************************************************************

  This file contains some very simple templates for integer manipulations.

******************************************************************************/

namespace atlas {

namespace intutils {

template<typename I> I abs(I a)

/*!
  Synopsis: returns the absolute value of a
*/

{
  return (a >= 0) ? a : -a;
}

template<typename I> I divide(I a, I b)

/*!
  The result of divide(a,b) is the unique element q such that a = q.b + r,
  with 0 <= r < b. Here b is assumed to be > 0, the sign of a may be arbitrary.
  For instance, divide(-1,2) should be -1, so that -1 = -1.2 + 1. On my
  machine, -1/2 is 0 (which is the other value accepted by the C standard.)
*/

{
  return (a >= 0) ? a/b : -(-a/b+1);
}

template<typename I> I factorial(I a)

/*!
  Synopsis: returns the factorial of a. It is assumed that a >= 0; we do
  not worry about overflow.
*/

{
  I f = 1;

  for (size_t j = 0; j < a; ++j)
    f *= j+1;

  return f;
}

}

}
