/*!
\file
  This is arithmetic.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef ARITHMETIC_H  /* guard against multiple inclusions */
#define ARITHMETIC_H

#include "arithmetic_fwd.h"

#include <iostream>

/******** function declarations **********************************************/

namespace atlas {

namespace arithmetic {

// operations without conversion for all integral types (all are inlined below)
  template<typename I> I abs(I);
  template<typename I> I min(I,I);
  template<typename I> I max(I,I);
  template<typename I> I divide(I, unsigned long);
  template<typename I> unsigned long remainder(I, unsigned long);
  //  template<typename I> I factorial(I); // moved to the size module

  extern unsigned long dummy_gcd, dummy_Bezout;

  unsigned long gcd (long, unsigned long); // signed first argument only!

  unsigned long unsigned_gcd (unsigned long, unsigned long); // avoid overload

  unsigned long div_gcd (unsigned long d, unsigned long a);

  unsigned long lcm (unsigned long, unsigned long,
		     unsigned long& =dummy_gcd, unsigned long& =dummy_Bezout);
  // after |lcm(a,b,d,p)| one has |d=gcd(a,b)| and |p%a==0|, |p%b==d|, |p<lcm|

  unsigned long power(unsigned long base, unsigned int exponent);

  // inlined; first two arguments are supposed already reduced modulo third
  unsigned long& modAdd(unsigned long&, unsigned long, unsigned long);

  unsigned long& modProd(unsigned long&, unsigned long, unsigned long);

} // |namespace arithmetic|

/*        Class definitions        */

namespace arithmetic {

class Rational
{
  long num;
  unsigned long denom;
public:
  explicit Rational (long n=0,unsigned long d=1) : num(n), denom(d)
  { normalize(); }

  long numerator() const            { return num; }
  unsigned long denominator() const { return denom; }

  Rational operator+(Rational q) const;
  Rational operator-(Rational q) const;
  Rational operator*(Rational q) const;
  Rational operator/(Rational q) const;

  Rational& operator=(Rational q)
    { num=q.num; denom=q.denom; return normalize(); }

  Rational& operator+=(Rational q) { return operator=(operator+(q)); }
  Rational& operator-=(Rational q) { return operator=(operator-(q)); }
  Rational& operator*=(Rational q) { return operator=(operator*(q)); }
  Rational& operator/=(Rational q) { return operator=(operator/(q)); }

  Rational& operator+=(long n);
  Rational& operator-=(long n);
  Rational& operator*=(long n);
  Rational& operator/=(long n);

  bool operator==(Rational q) const { return num*q.denom==denom*q.num; }
  bool operator!=(Rational q) const { return num*q.denom!=denom*q.num; }
  bool operator<(Rational q)  const { return num*q.denom<denom*q.num; }
  bool operator<=(Rational q) const { return num*q.denom<=denom*q.num; }
  bool operator>(Rational q)  const { return num*q.denom>denom*q.num; }
  bool operator>=(Rational q) const { return num*q.denom>=denom*q.num; }

  inline Rational& normalize();
  Rational& power(int n); // raise to power |n| and return |*this|

}; // |class Rational|


class Split_integer
{
  int real_part, s_part;
 public:
  explicit Split_integer(int a=0, int b=0) : real_part(a), s_part(b) {}

  int e() const { return real_part; }
  int s() const { return s_part; }

  bool operator== (Split_integer y) const { return e()==y.e() and s()==y.s(); }
  bool operator!= (Split_integer y) const { return not operator==(y); }

  Split_integer& operator +=(Split_integer y)
  { real_part+=y.e(); s_part+=y.s(); return *this; }
  Split_integer& operator -=(Split_integer y)
  { real_part-=y.e(); s_part-=y.s(); return *this; }
  Split_integer operator +(Split_integer y) { return y+= *this; }
  Split_integer operator -(Split_integer y) { return y.negate()+= *this; }
  Split_integer operator -() const { return Split_integer(*this).negate(); }
  Split_integer operator* (Split_integer y) const
  { return Split_integer(e()*y.e()+s()*y.s(),e()*y.s()+s()*y.e()); }

  Split_integer& operator*= (int n) { real_part*=n; s_part*=n; return *this; }
  Split_integer& operator*= (Split_integer y) { return *this=operator*(y); }
  Split_integer& negate() { real_part=-e(); s_part=-s(); return *this; }
  Split_integer& times_s() { std::swap(real_part,s_part); return *this; }
  Split_integer& times_1_s() { real_part= -(s_part-=e()); return *this; }

}; // |class Split_integer|

std::ostream& operator<< (std::ostream& out, const Rational& frac);


/******** inline function definitions ***************************************/

// Return the absolute value/min/max, without changing the size
template<typename I>
  inline I abs(I a) { return a >= 0 ? a : -a; }
template<typename I>
  inline I min(I a,I b) { return a<b ? a : b; }
template<typename I>
  inline I max(I a,I b) { return a<b ? b : a; }

/*! The result of |divide(a,b)| is the unique integer $q$ with $a = q.b + r$,
  and $0 \leq r < b$. Here the sign of |a| may be arbitrary, the requirement
  for |r| assumes |b| positive, which is why it is passed as unsigned (also
  this better matches the specification of |remainder| below). Callers must
  make sure that $b$ is positive, since implicit conversion of a negative
  signed value to unsigned would wreak havoc.

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
  { return a >= 0 ? a/b : ~(~a/b); } // left operand is safely made unsigned

// override for |I=unsigned long| (is instantiated for |Matrix<unsigned long>|)
inline unsigned long divide (unsigned long a, unsigned long b)
  { return a/b; } // unsigned division is OK in this case

/*!
  Synopsis: returns the remainder of the division of a by b.

  The point is to allow |I| to be a signed type, and avoid the catastrophic
  implicit conversion to unsigned when using the '%' operation. Also corrects
  the deficiency of 'signed modulo' by always returning the unique number |r|
  in [0,m[ such that $a = q.b + r$, in other words with |q=divide(a,b)| above.

  NOTE: For $a<0$ one should \emph{not} return |m - (-a % m)|; this fails
  iff $m$ divides $a$. However replacing |-| by |~|, which maps $a\mapsto-1-a$
  and satifies |~(m*q+r)==m*~q+(m+~r)|, the result is always correct.
*/
template<typename I>
  inline unsigned long remainder(I a, unsigned long b)
  { return a >= 0 ? a%b : b+~(~a%b); } // conversions to unsigned are safe here



  inline unsigned long div_gcd (unsigned long d, unsigned long a)
  { return d/unsigned_gcd(a,d); }

  inline unsigned long gcd (long a, unsigned long b)
  {
    if (a < 0)
      return unsigned_gcd(static_cast<unsigned long>(-a),b);
    else
      return unsigned_gcd(static_cast<unsigned long>(a),b);
  }

  inline unsigned long& modAdd(unsigned long& a, unsigned long b,
			       unsigned long n)
  {
    if (a < n-b)
      a += b;
    else
      a -= n-b;
    return a;
  }

  inline Rational& Rational::normalize()
  {
    unsigned long d = gcd(num,denom);
    if (d>1)
      num/=long(d),denom/=d;
    return *this;
  }

} // |namespace arithmetic|

} // |namespace atlas|

#endif
