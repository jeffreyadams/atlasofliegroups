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
{ return d_num*v.d_denom == v.d_num*d_denom; } // cross multiply and compare

bool RatLatticeElt::operator<(const RatLatticeElt& v) const
{ // cross multiply component-wise, and compare
  for (size_t i=0; i<d_num.size(); ++i)
  {
    LatticeCoeff d= d_num[i]*v.d_denom - v.d_num[i]*d_denom;
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
  return result.normalize();
}

RatLatticeElt& RatLatticeElt::normalize()
{
  unsigned long d=d_denom;

  if (d==0)
    throw std::runtime_error("Denominator 0 in rational vector");

  for (size_t i=0; i<d_num.size(); ++i)
  {
    if (d_num[i]!=0)
      d=arithmetic::gcd(d,abs(d_num[i]));
    if (d==1)
      return *this;
  }

  d_denom/=d;
  d_num/=d;
  return *this;
}

} // namespace latticetypes

} // namespace atlas
