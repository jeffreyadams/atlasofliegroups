/*!
\file
  This is arithmetic_def.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

namespace atlas {

namespace arithmetic {


/*!
  Synopsis: returns the remainder of the division of c modulo m.

  The point is to always return the unique number r in [0,m[ such that
  a = q.m + r for a (unique) q. The difficulty arises when C is a signed
  type; if we are not careful the sign conversions will mess things up
  big time.

  NOTE: For $a < 0$ one should not always return |m - (-a % m)|; this fails
  iff $m$ divides $a$. However replacing |-| by |~|, which maps $a\mapsto-1-a$
  and satifies |~(m*q+r)==m*~q+(m+~r)|, the result is always correct.
*/
template<typename C> unsigned long remainder(C a, unsigned long m)
{
  return a>=0 ? a%m : m + ~(~a % m);
}

} // |namespace arithmetic|

} // |namespace atlas|
