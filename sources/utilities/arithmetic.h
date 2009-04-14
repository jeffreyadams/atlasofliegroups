/*!
\file
  This is arithmetic.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef ARITHMETIC_H  /* guard against multiple inclusions */
#define ARITHMETIC_H

#include <iostream>
#include <vector>

/******** function declarations **********************************************/

namespace atlas {

namespace arithmetic {

  extern unsigned long dummy;

  unsigned long gcd (long, unsigned long); // signed first argument only!

  unsigned long unsigned_gcd (unsigned long, unsigned long); // avoid overload

  unsigned long lcm (unsigned long, unsigned long, unsigned long& = dummy);

  // inlined; first two arguments are supposed already reduced modulo third
  unsigned long& modAdd(unsigned long&, unsigned long, unsigned long);

  unsigned long& modProd(unsigned long&, unsigned long, unsigned long);

  template<typename C> unsigned long remainder(C, unsigned long);

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

  bool operator==(Rational q) const { return num*q.denom==denom*q.num; }
  bool operator!=(Rational q) const { return num*q.denom!=denom*q.num; }
  bool operator<(Rational q)  const { return num*q.denom<denom*q.num; }
  bool operator<=(Rational q) const { return num*q.denom<=denom*q.num; }
  bool operator>(Rational q)  const { return num*q.denom>denom*q.num; }
  bool operator>=(Rational q) const { return num*q.denom>=denom*q.num; }

  inline Rational& normalize();

}; // |class Rational|

std::ostream& operator<< (std::ostream& out, const Rational& frac);

typedef std::vector<Rational> RationalList;

} // |namespace arithmetic|

/******** inline function definitions ***************************************/

namespace arithmetic {

  inline unsigned long gcd (long a, unsigned long b) {
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
      num/=d,denom/=d;
    return *this;
  }

} // |namespace arithmetic|

} // |namespace atlas|

#include "arithmetic_def.h"

#endif
