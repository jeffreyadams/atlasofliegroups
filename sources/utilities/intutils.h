/*!
\file
  This is intutils.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef INTUTILS_H  /* guard against multiple inclusions */
#define INTUTILS_H

/******** function declarations *********************************************/

namespace atlas {

namespace intutils {

  template<typename I> I abs(I);
  template<typename I> I divide(I, unsigned long);
  template<typename I> unsigned long remainder(I, unsigned long);
  template<typename I> I factorial(I);

  // inline definitions

/*!
  Synopsis: returns the absolute value of a, without changing the size
*/
template<typename I>
  inline I abs(I a) { return a >= 0 ? a : -a; }

/*!
  The result of divide(a,b) is the unique element q such that a = q.b + r,
  with 0 <= r < b. Here the sign of a may be arbitrary, the requirement for r
  assumes b positive, which is why it is passed as unsigned (also this better
  matches the specifiation of |remainder| below). This function would probably
  work correctly even if called with b is given by a negative value of type I,
  since conversion to unsigned long and back to I will be identity, and
  hardware division very probably satisfies x/(-y) = -(x/y) for signed types,
  even if the C++ standard fails to require this unambiguously. Nevertheless,
  callers ought to make sure that b is positive.

  Hardware division probably does _not_ handle negative |a| correctly; for
  instance, divide(-1,2) should be -1, so that -1 = -1.2 + 1, but on my
  machine, -1/2 is 0 (which is the other value accepted by the C standard;
  Fokko.) [Note that the correct symmetry to apply to |a|, one that maps
  classes with the same quotient to each other, is not \f$a\to -a\f$ but
  \f$a\to -1-a\f$, where the latter value can be conveniently written as |~a|
  in C or C++. Amazingly Fokko's incorrect original expresion |-(-a/b -1)|
  never did any harm. MvL]
*/
template<typename I>
  inline I divide(I a, unsigned long b)
  { return a >= 0 ? a/(I)b : ~(~a/(I)b); }

template<> // the following gets instantiated for |Matrix<unsigned long>|
inline unsigned long divide<unsigned long>
  (unsigned long a, unsigned long b) { return a/b; }


/*!
  Synopsis: returns the remainder of the division of a by b.

  The point is to always return the unique number r in [0,m[ such that
  a = q.b + r, in other words with q=divide(a,b) above.

  NOTE: For $a < 0$ one should not always return |m - (-a % m)|; this fails
  iff $m$ divides $a$. However replacing |-| by |~|, which maps $a\mapsto-1-a$
  and satifies |~(m*q+r)==m*~q+(m+~r)|, the result is always correct.
*/
template<typename I>
  inline unsigned long remainder(I a, unsigned long b)
  { return a >= 0 ? a%b : b+~(~a%b); }


} // |namespace intutils|

} // |namespace atlas|

#endif
