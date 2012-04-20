/*
  This is y_values.cpp.

  Copyright (C) 2011 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "y_values.h"

#include "tits.h" // |TorusPart|

#include "arithmetic.h"
#include "prerootdata.h"
#include "rootdata.h"

namespace atlas {

namespace y_values {
//              |TorusElement|


// For torus elements we keep numerator entries reduced modulo |2*d_denom|.
// Note this invariant is not enforced in the (private) raw constructor.

// make $\exp(2\pi i r)$, doing affectively *2, then reducing modulo $2X_*$
TorusElement::TorusElement(const RatWeight& r, bool two)
  : repr(r) // but possibly multiplied by 2 below
{ if (two)
    repr*=2;
  unsigned int d=2u*repr.denominator(); // we reduce modulo $2\Z^rank$
  Weight& num=repr.numerator();
  for (size_t i=0; i<num.size(); ++i)
    num[i] = arithmetic::remainder(num[i],d);
}

RatWeight TorusElement::log_pi(bool normalize) const
{
  if (normalize)
    return RatWeight(repr).normalize();
  return repr;
}

RatWeight TorusElement::log_2pi() const
{
  Weight numer = repr.numerator(); // copy
  unsigned int d = 2u*repr.denominator(); // this will make result mod Z, not 2Z
  return RatWeight(numer,d).normalize();
}

// evaluation giving rational number modulo 2
Rational TorusElement::evaluate_at (const Coweight& alpha) const
{
  unsigned int d = repr.denominator();
  int n = arithmetic::remainder(repr.numerator().dot(alpha),d+d);
  return Rational(n,d);
}

TorusElement TorusElement::operator +(const TorusElement& t) const
{
  TorusElement result(repr + t.repr,tags::UnnormalizedTag()); // raw ctor
  int d=2*result.repr.denominator(); // we shall reduce modulo $2\Z^rank$
  Weight& num=result.repr.numerator();
  for (size_t i=0; i<num.size(); ++i)
    if (num[i] >= d) // correct if |result.repr| in interval $[2,4)$
      num[i] -=d;
  return result;
}

TorusElement TorusElement::operator -(const TorusElement& t) const
{
  TorusElement result(repr - t.repr,0); // raw constructor
  int d=2*result.repr.denominator(); // we shall reduce modulo $2\Z^rank$
  Weight& num=result.repr.numerator();
  for (size_t i=0; i<num.size(); ++i)
    if (num[i]<0) // correct if |result.repr| in interval $(-2,0)$
      num[i] +=d;
  return result;
}

TorusElement& TorusElement::operator+=(TorusPart v)
{
  for (size_t i=0; i<v.size(); ++i)
    if (v[i])
    {
      if (repr.numerator()[i]<repr.denominator())
	repr.numerator()[i]+=repr.denominator(); // add 1/2 to coordinate
      else
	repr.numerator()[i]-=repr.denominator(); // subtract 1/2
    }
  return *this;
}


void TorusElement::simple_reflect(const PreRootDatum& prd, weyl::Generator s)
{ prd.simpleReflect(repr.numerator(),s); } // numerator is weight for |prd|

void TorusElement::reflect(const RootDatum& rd, RootNbr alpha)
{ rd.reflect(repr.numerator(),alpha); } // numerator is weight for |rd|


TorusElement TorusElement::simple_imaginary_cross
  (const RootDatum& dual_rd, // dual for pragmatic reasons
   RootNbr alpha) const // any simple-imaginary root
{
  TorusElement t(*this); // make a copy in all cases
  if (not negative_at(dual_rd.coroot(alpha))) // |alpha| is a noncompact root
    t += tits::TorusPart(dual_rd.root(alpha)); // so add $m_\alpha$
  return t;
}


size_t y_entry::hashCode(size_t modulus) const
{
  unsigned long d= fingerprint.denominator()+1;
  const int_Vector& num=fingerprint.numerator();
  size_t h=nr; // start with involution number
  for (size_t i=0; i<num.size(); ++i)
    h=h*d+num[i];
  return h&(modulus-1);
}

bool y_entry::operator !=(const y_entry& y) const
{ return nr!=y.nr or fingerprint!=y.fingerprint; }

} // namsepace

} // namespace atlas
