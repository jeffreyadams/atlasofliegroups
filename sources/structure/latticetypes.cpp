/*!
\file
  \brief This module defines a few simple operations on the types defined in
  latticetypes.h.
*/
/*
  This is latticetypes.cpp.

  Copyright (C) 2009 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "latticetypes.h"
#include "arithmetic.h"
#include <stdexcept>

/*****************************************************************************

  This module implements the non-inline method(s) of the type(s) defined in
  latticetypes.h

******************************************************************************/

namespace atlas {

namespace latticetypes {

bool RatLatticeElt::operator==(const RatLatticeElt& v) const
{ return d_num*(int)v.d_denom == v.d_num*(int)d_denom; } // cross multiply

bool RatLatticeElt::operator<(const RatLatticeElt& v) const
{ // cross multiply component-wise, and compare
  for (size_t i=0; i<d_num.size(); ++i)
  {
    LatticeCoeff d= d_num[i]*(int)v.d_denom - v.d_num[i]*(int)d_denom;
    if (d!=0)
      return d<0;
  }
  return false; // equality if we get here
}

RatLatticeElt RatLatticeElt::operator+(const RatLatticeElt& v) const
{
  unsigned long gcd, m = arithmetic::lcm(d_denom,v.d_denom,gcd);
  unsigned long f = v.d_denom/gcd;
  assert (f==m/d_denom); // if this fails, then there was overflow on m
  RatLatticeElt result(d_num*f,m);
  result.d_num += v.d_num*(d_denom/gcd);
  return result; // don't normalize, better just limit denominator growth
}

// rational vectors are not guaranteed on lowest terms
// however if multiplication can be by cancellation it is done that way
RatLatticeElt& RatLatticeElt::operator*=(int n)
{
  if (n!=0 and (int)d_denom%n==0)
    if (n>0)
      d_denom/=(unsigned int)n; // unsigned arithmetic OK here
    else
    {
      d_num=-d_num;
      d_denom/=(unsigned int)(-n);
    }
  else d_num*=n;
  return *this;
}

// division takes signed argument to avoid catastrophic surprises: if ever a
// negative argument were passed, implicit conversion to unsigned would be fatal
RatLatticeElt& RatLatticeElt::operator/=(int n)
{
  if (n>0)
    d_denom*=(unsigned int)n;
  else
  {
    assert(n!=0);
    d_num=-d_num;
    d_denom*=(unsigned int)(-n);
  }
  return *this;
}

RatLatticeElt& RatLatticeElt::normalize()
{
  unsigned long d=d_denom;

  if (d==0)
    throw std::runtime_error("Denominator 0 in rational vector");

  for (size_t i=0; i<d_num.size(); ++i)
  {
    if (d_num[i]!=0)
      d=arithmetic::gcd(d,arithmetic::abs(d_num[i]));
    if (d==1)
      return *this;
  }

  d_denom/=d;
  d_num/=d;
  return *this;
}

} // namespace latticetypes

} // namespace atlas
