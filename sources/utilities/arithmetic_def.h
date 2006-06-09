/*!
\file
  This is arithmetic_def.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

namespace atlas {

namespace arithmetic {

template<typename C> unsigned long remainder(C a, unsigned long m)

/*!
  Synopsis: returns the remainder of the division of c modulo m.

  The point is to always return the unique number r in [0,m[ such that
  a = q.m + r for a (unique) q. The difficulty arises when C is a signed
  type; if we are not careful the sign conversions will mess things up
  big time.

  NOTE: this looks rather clumsy but I didn't see a better way. It is not
  always true that for a < 0 one should return m - (-a % m); this fails
  iff a is divisible by m.
*/

{
  if (a >= 0)
    return a % m;

  C b = -a % m;

  if (b)
    return m - b;
  else
    return 0;
}

}

}
