// This is arithmetic.h
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef ARITHMETIC_H  /* guard against multiple inclusions */
#define ARITHMETIC_H

#include "arithmetic_fwd.h"

#include <iostream>
#include <cassert>

/******** function declarations **********************************************/

namespace atlas {

namespace arithmetic {

// operations without conversion for signed integral types
// second argument is signed (to simplify interface),
// but it must be non-negative
  template<typename I> I divide(I, I);
  template<typename I> I remainder(I, I);
  //  template<typename I> I factorial(I); // moved to the size module

  extern Denom_t dummy_gcd, dummy_mult;

  Denom_t gcd (Numer_t, Denom_t); // signed first argument only! Defined below.

  // the following functions are entirely unsigned

  // the name of the next function was chosen to avoid any overloading
  Denom_t unsigned_gcd (Denom_t a, Denom_t b); // must have |b!=0|; |a| is free

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


/*        Class definitions        */

template<typename I>
  class Rational
{
  using UI = typename std::make_unsigned<I>::type;
  I num;
  UI denom;
public:
  explicit Rational (Numer_t n=0,Denom_t d=1) : num(n), denom(d)
  { normalize(); }

  I numerator() const   { return num; }

  /*
    The C++ rule that one unsigned operand silently converts the other operand
    to unsigned as well makes exporting denominator as unsigned too error prone;
    e.g., floor=numerator()/denominator() would wreak havoc. Generally speaking,
    casting unsigned to signed is not a good thing, but when this gives a
    negative value here, we risk having (had) overflow anyway, so here it is.
  */
  I denominator() const { return static_cast<I>(denom); }

  // however sometimes (printing, |big_int| conversion, we want the Real Thing
  UI true_denominator() const { return denom; }

  bool is_zero() const { return num==0; }
  bool is_positive() const { return num>0; }
  bool is_negative() const { return num<0; }

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
  Rational& operator+=(I n) { num+=n*denom; return *this; }
  Rational& operator-=(I n) { num-=n*denom; return *this; }
  Rational& operator*=(I n);
  Rational& operator/=(I n); // assumes $n\neq0$, will not throw
  Rational& operator%=(I n); // assumes $n\neq0$, will not throw
  Rational& mod1() { num = remainder(num,denominator()); return *this; };

  I floor () const { return divide(num,static_cast<I>(denom)); }
  I ceil () const { return -divide(-num,static_cast<I>(denom)); }
  I quotient (Denom_t n) const // integer division by positive integer
    { return divide(num,static_cast<I>(n*denom)); }

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

  const Rational& normalize() const; // pseudo-|const|; defined way below
  Rational& power(int n); // raise to power |n| and return |*this|

  Rational& normalize() // non-|const| version delagates to (pseudo) |const| one
  { static_cast<const Rational*>(this)->normalize(); return *this; }
}; // |class Rational|

//				Functions


// arithmetic operations by value
template<typename C> Rational<C> operator+ (Rational<C> q, C n) { return q+=n ;}
template<typename C> Rational<C> operator- (Rational<C> q, C n) { return q-=n ;}
template<typename C> Rational<C> operator* (Rational<C> q, C n) { return q*=n ;}
template<typename C> Rational<C> operator/ (Rational<C> q, C n) { return q/=n ;}
template<typename C> Rational<C> operator% (Rational<C> q, C n) { return q%=n ;}

class Split_integer
{
  int ev_1, ev_minus_1; // store evaluations at $1$ and $-1$ for efficiency
  // class invariant: |ev_1| and |ev_minus_1| have the same parity

  struct raw_tag {}; // to signal use of the private constructor
  explicit constexpr Split_integer(int ev_1, int ev_minus_1, raw_tag)
  : ev_1(ev_1), ev_minus_1(ev_minus_1) {}
 public:
  explicit constexpr Split_integer(int a=0, int b=0)
  : ev_1(a+b), ev_minus_1(a-b) {}

  int e() const { return (static_cast<long long int>(ev_1)+ev_minus_1)/2; }
  int s() const { return (static_cast<long long int>(ev_1)-ev_minus_1)/2; }

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

  Split_integer operator +(Split_integer y) const
  { return Split_integer(ev_1+y.ev_1,ev_minus_1+y.ev_minus_1,raw_tag()); }
  Split_integer operator -(Split_integer y) const
  { return Split_integer(ev_1-y.ev_1,ev_minus_1-y.ev_minus_1,raw_tag()); }
  Split_integer operator -() const
  { return Split_integer(-ev_1,-ev_minus_1,raw_tag()); }

  // multiplication is where the "split" representation really wins out:
  Split_integer& operator*= (int n) { ev_1*=n; ev_minus_1*=n; return *this; }
  Split_integer& operator*= (Split_integer y)
  { ev_1*=y.ev_1; ev_minus_1*=y.ev_minus_1; return *this; }
  Split_integer operator* (int n) const
  { return Split_integer(ev_1*n,ev_minus_1*n,raw_tag()); }
  Split_integer operator* (Split_integer y) const
  { return Split_integer(ev_1*y.ev_1,ev_minus_1*y.ev_minus_1,raw_tag()); }

  Split_integer& negate() { ev_1=-ev_1; ev_minus_1=-ev_minus_1; return *this; }

  Split_integer& times_s() { ev_minus_1=-ev_minus_1; return *this; }
  Split_integer times_s() const
  { return Split_integer(ev_1,-ev_minus_1,raw_tag()); }
  Split_integer& times_1_s() // multiply by |1-s|
  { ev_1=0; /* see "Who Killed the Electric Car" */
    ev_minus_1*=2; return *this; }
  Split_integer times_1_s() const
  { return Split_integer(0,2*ev_minus_1,raw_tag()); }

  Split_integer& times_1_plus_s() // multiply by |1+s|
  { ev_1*=2; ev_minus_1=0; return *this; }
  Split_integer times_1_plus_s() const
  { return Split_integer(2*ev_1,0,raw_tag()); }

  int s_to_1() const { return ev_1; }
  int s_to_minus_1() const { return ev_minus_1; }
}; // |class Split_integer|

template<typename I>
  std::ostream& operator<< (std::ostream& out, const Rational<I>& frac);


/******** inline function definitions ***************************************/

/*
  The result of |divide(a,b)| is the unique integer $q$ with $a = q.b + r$, and
  $0 \leq r < b$. Here the sign of |a| may be arbitrary, the requirement for |r|
  assumes |b| positive. So this function assumes |b>=0|, even though |b| is
  passed in a signed type argument (which also matches the specification of
  |remainder| below). Callers must make sure that $b$ is positive, since
  implicit conversion of a negative signed value to unsigned would wreak havoc.

  Hardware division probably does _not_ handle negative |a| correctly; for
  instance, divide(-1,2) should be -1, so that -1 = -1.2 + 1, but on my
  machine, -1/2 is 0 (which is the other value accepted by the C standard;
  Fokko.) [Note that the correct symmetry to apply to |a|, one that maps
  classes with the same quotient to each other, is not $a\to -a$ but $a\to
  -1-a$. The latter value can be conveniently written as |~a| if |a| is of an
  unsigned type, but that is not the case for the arguments here, and the
  standards refuse to clearly specify complementing on signed values. So for
  negative values |a| we spell out subtraction from $-1$ below, hoping that
  the optimiser will emit complementation for it; however once a signed value
  has been made positive, we do unsigned division (after converting implicitly
  because unsigned operator arguments beat signed), and finally convert back
  to signed, taking care to do any possibly negative computation after the
  conversion to signed (but remainders are never negative).]

  [Amazingly Fokko's incorrect original expression |-(-a/b -1)| never did any
  harm. MvL]
*/

template<typename I>
  I divide(I a, I b_signed) // integer quotient, may be negative
{ using UI = typename std::make_unsigned<I>::type;
  UI b(b_signed); // interpret unsigned, even though it was passed as signed
  // use unsigned division (because |b| is so), then convert back to signed
  return a >= 0 ? static_cast<I>(a/b) : -1-static_cast<I>((-1-a)/b);
}

// override for |I=Denom_t| (is instantiated for |Matrix<Denom_t>|)
template<>
  inline Denom_t divide<Denom_t> (Denom_t a, Denom_t b)
  { return a/b; } // unsigned division is OK in this case

/*
  Return the remainder of the division of |a| (signed) by |b| (unsigned).

  The point is to allow |I| to be a signed type, and avoid the catastrophic
  implicit conversion to unsigned when using the '%' operation. Also corrects
  the deficiency of 'signed modulo' by always returning the unique number |r|
  in [0,m[ such that $a = q.b + r$, in other words with |q=divide(a,b)| above.

  NOTE: For $a<0$ one should \emph{not} return |b - (-a % b)|; this fails when
  $b$ divides $a$. However replacing |-| by |~|, which maps $a\mapsto-1-a$
  and satisfies |~(q*b+r)==~q*b+(b+~r)|, the result is always correct.
*/
template<typename I>
  I remainder(I a, I b_signed)
{ assert(b_signed>0); // should avoid nasty surprises; negatives fail badly
  using UI = typename std::make_unsigned<I>::type;
  UI b(b_signed); // interpret unsigned, even though it was passed as signed
  // use unsigned division (because |b| is so), then convert back to signed
  return static_cast<I>(a >= 0 ? a%b : ~((-1-a)%b) + b); // never negative
}

inline Denom_t div_gcd (Denom_t d, Denom_t a) { return d/unsigned_gcd(a,d); }

inline Denom_t gcd (Numer_t a, Denom_t b) // caller must ensure |b>0|
{
  if (a > 0)
    return unsigned_gcd(static_cast<Denom_t>(a),b);
  else if (a==0) return b; // faster, although |unsigned_gcd| could also cope
  else
    return unsigned_gcd(static_cast<Denom_t>(-a),b);
}

// we assume |a| and |b| to be less than |n| here
inline Denom_t modAdd(Denom_t a, Denom_t b, Denom_t n)
{
  if (a < n-b)
    return a + b;
  else
    return a - (n-b);
}

template<typename I>
  const Rational<I>& Rational<I>::normalize() const
{
  Denom_t d = gcd(num,denom);
  if (d>1)
  { auto& my = const_cast<Rational<I>&>(*this);
    my.num/=static_cast<I>(d);
    my.denom/=d;
  }
  return *this;
}

} // |namespace arithmetic|

} // |namespace atlas|

#endif
