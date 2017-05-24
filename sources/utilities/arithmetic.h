// This is arithmetic.h
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2016 Marc van Leeuwen
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
  template<typename I> I divide(I, Denom_t);
  template<typename I> Denom_t remainder(I, Denom_t);
  //  template<typename I> I factorial(I); // moved to the size module

  extern Denom_t dummy_gcd, dummy_mult;

  Denom_t gcd (Numer_t, Denom_t); // signed first argument only!

  // the following functions are entirely unsigned
  Denom_t unsigned_gcd (Denom_t, Denom_t); // name choice avoids overloading

  Denom_t div_gcd (Denom_t a, Denom_t b); // $a/\gcd(a,b)$

  Denom_t lcm (Denom_t a, Denom_t b,
	       Denom_t& gcd=dummy_gcd, Denom_t& a_mult=dummy_mult);
  // after |lcm(a,b,d,p)| one has |d=gcd(a,b)| and |p%a==0|, |p%b==d|, |p<lcm|

  Denom_t power(Denom_t base, unsigned int exponent);

  template<typename I> int exp_minus_1 (I n) { return n%2==0 ? 1 : -1; }
  template<typename I> int exp_i (I n)
  { assert(n%2==0); return n%4==0 ? 1 : -1; }

  // inlined; first two arguments are supposed already reduced modulo third
  Denom_t modAdd(Denom_t, Denom_t, Denom_t);

  Denom_t modProd(Denom_t, Denom_t, Denom_t);

} // |namespace arithmetic|

/*        Class definitions        */

namespace arithmetic {

class Rational
{
  Numer_t num;
  Denom_t denom;
public:
  explicit Rational (Numer_t n=0,Denom_t d=1) : num(n), denom(d)
  { normalize(); }

  Numer_t numerator() const   { return num; }

  /* the C++ rule that one unsigned operand silently converts the other
     operand to unsigned as well makes exporting denominator as unsigned too
     error prone; e.g., floor=numerator()/denominator() would wreak havoc */
  Numer_t denominator() const { return Numer_t(denom); }

  // these operators all return normalised results
  Rational operator+(Rational q) const;
  Rational operator-(Rational q) const;
  Rational operator*(Rational q) const;
  Rational operator/(Rational q) const; // assumes $q\neq0$, will not throw
  Rational operator%(Rational q) const; // assumes $q\neq0$, will not throw

  Rational operator-() const { return Rational(-num,denom); } // uniary minus
  Rational inverse() const // multiplicative inverse; nonzero value assumed
  { return num>0 ? Rational(denom,num) : Rational(-denom,-num); }

  Rational& operator=(Rational q)
    { num=q.num; denom=q.denom; return normalize(); }

  Rational& operator+=(Rational q) { return operator=(operator+(q)); }
  Rational& operator-=(Rational q) { return operator=(operator-(q)); }
  Rational& operator*=(Rational q) { return operator=(operator*(q)); }
  Rational& operator/=(Rational q) { return operator=(operator/(q)); }
  Rational& operator%=(Rational q) { return operator=(operator%(q)); }

  // assignment operators with integers have efficient implementations
  Rational& operator+=(Numer_t n);
  Rational& operator-=(Numer_t n);
  Rational& operator*=(Numer_t n);
  Rational& operator/=(Numer_t n); // assumes $n\neq0$, will not throw
  Rational& operator%=(Numer_t n); // assumes $n\neq0$, will not throw

  Numer_t floor () const { return divide(num,denom); }
  Numer_t ceil () const { return -divide(-num,denom); }
  Numer_t quotient (Denom_t n) const { return divide(num,n*denom); }

  // these definitions must use |denominator()| to ensure signed comparison
  bool operator==(Rational q) const
    { return num*q.denominator()==denominator()*q.num; }
  bool operator!=(Rational q) const
    { return num*q.denominator()!=denominator()*q.num; }
  bool operator<(Rational q)  const
    { return num*q.denominator()<denominator()*q.num; }
  bool operator<=(Rational q) const
    { return num*q.denominator()<=denominator()*q.num; }
  bool operator>(Rational q)  const
    { return num*q.denominator()>denominator()*q.num; }
  bool operator>=(Rational q) const
    { return num*q.denominator()>=denominator()*q.num; }

  inline Rational& normalize();
  Rational& power(int n); // raise to power |n| and return |*this|

}; // |class Rational|


class Split_integer
{
  int ev_1, ev_minus_1; // store evaluations at $1$ and $-1$ for efficiency
  // class invariant: |ev_1| and |ev_minus_1| have the same parity
 public:
  explicit constexpr Split_integer(int a=0, int b=0)
  : ev_1(a+b), ev_minus_1(a-b) {}

  int e() const { return (ev_1+ev_minus_1)/2; }
  int s() const { return (ev_1-ev_minus_1)/2; }

  bool operator== (Split_integer y) const
  { return ev_1==y.ev_1 and ev_minus_1==y.ev_minus_1; }
  bool operator!= (Split_integer y) const { return not operator==(y); }
  bool is_zero() const { return ev_1==0 and ev_minus_1==0; }

  Split_integer& operator +=(int n) { ev_1+=n; ev_minus_1+=n; return *this; }
  Split_integer& operator -=(int n) { ev_1-=n; ev_minus_1-=n; return *this; }
  Split_integer& operator +=(Split_integer y)
  { ev_1+=y.ev_1; ev_minus_1+=y.ev_minus_1; return *this; }
  Split_integer& operator -=(Split_integer y)
  { ev_1-=y.ev_1; ev_minus_1-=y.ev_minus_1; return *this; }
  Split_integer operator +(Split_integer y) { return y+= *this; }
  Split_integer operator -(Split_integer y) { return y.negate()+= *this; }
  Split_integer operator -() const { return Split_integer(*this).negate(); }

  Split_integer& operator*= (int n) { ev_1*=n; ev_minus_1*=n; return *this; }
  Split_integer operator* (int n) const { return Split_integer(*this)*=n; }
  // multiplication is where the "split" representation really wins out:
  Split_integer& operator*= (Split_integer y)
  { ev_1*=y.ev_1; ev_minus_1*=y.ev_minus_1; return *this; }
  Split_integer operator* (Split_integer y) const { return y*= *this; }

  Split_integer& negate() { ev_1=-ev_1; ev_minus_1=-ev_minus_1; return *this; }
  Split_integer& times_s() { ev_minus_1=-ev_minus_1; return *this; }
  Split_integer& times_1_s() { ev_1=0; // see "Who Killed the ELectric Car"
    ev_minus_1*=2; return *this; }

  int s_to_1() const { return ev_1; }
  int s_to_minus_1() const { return ev_minus_1; }
}; // |class Split_integer|

std::ostream& operator<< (std::ostream& out, const Rational& frac);


/******** inline function definitions ***************************************/

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
  classes with the same quotient to each other, is not $a\to -a$ but
  $a\to -1-a$, where the latter value can be conveniently written as |~a|
  in C or C++. Amazingly Fokko's incorrect original expresion |-(-a/b -1)|
  never did any harm. MvL]
*/
template<typename I>
  inline I divide(I a, Denom_t b)
  { return a >= 0 ? a/b : ~(~a/b); } // left operand is safely made unsigned

// override for |I=Denom_t| (is instantiated for |Matrix<Denom_t>|)
inline Denom_t divide (Denom_t a, Denom_t b)
  { return a/b; } // unsigned division is OK in this case

/*
  Return the remainder of the division of |a| (signed) by |b| (unsigned).

  The point is to allow |I| to be a signed type, and avoid the catastrophic
  implicit conversion to unsigned when using the '%' operation. Also corrects
  the deficiency of 'signed modulo' by always returning the unique number |r|
  in [0,m[ such that $a = q.b + r$, in other words with |q=divide(a,b)| above.

  NOTE: For $a<0$ one should \emph{not} return |b - (-a % b)|; this fails when
  $b$ divides $a$. However replacing |-| by |~|, which maps $a\mapsto-1-a$
  and satifies |~(q*b+r)==~q*b+(b+~r)|, the result is always correct.
*/
template<typename I>
  inline Denom_t remainder(I a, Denom_t b)
  { return a >= 0 ? a%b // safe implicit conversion to unsigned here
      : b+~(~static_cast<Denom_t>(a)%b); // safe explicit conversion here
  }

  inline Denom_t div_gcd (Denom_t d, Denom_t a) { return d/unsigned_gcd(a,d); }

  inline Denom_t gcd (Numer_t a, Denom_t b)
  {
    if (a < 0)
      return unsigned_gcd(static_cast<Denom_t>(-a),b);
    else
      return unsigned_gcd(static_cast<Denom_t>(a),b);
  }

  // we assume |a| and |b| to be less than |n| here
  inline Denom_t modAdd(Denom_t a, Denom_t b, Denom_t n)
  {
    if (a < n-b)
      return a + b;
    else
      return a - (n-b);
  }

  inline Rational& Rational::normalize()
  {
    Denom_t d = gcd(num,denom);
    if (d>1)      num/=Numer_t(d),denom/=d;
    return *this;
  }

} // |namespace arithmetic|

} // |namespace atlas|

#endif
