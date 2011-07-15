/*!
\file
  This is ratvec.cpp. This module contains some simple utilities for matrices.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "ratvec.h"

#include <cassert>
#include <stdexcept>
#include <vector>

#include "arithmetic.h"

namespace atlas {

namespace ratvec {

//		      The RationalVector class template

template<typename C>
RationalVector<C>::RationalVector(const V& v, C d)
 : d_num(v), d_denom(arithmetic::abs(d))
{ if (d<C(0)) d_num*=-C(1); }

// the following implementation assumes |long| can hold cross products
template<typename C>
  bool RationalVector<C>::operator==(const RationalVector<C>& v) const
{ return  // cross multiply
    d_num*(signed long)v.d_denom == v.d_num*(signed long)d_denom;
}

template<typename C>
  bool RationalVector<C>::operator<(const RationalVector<C>& v) const
{ // cross multiply component-wise, and compare
  for (size_t i=0; i<d_num.size(); ++i)
  {
    C d = d_num[i]*(signed long)v.d_denom - v.d_num[i]*(signed long)d_denom;
    if (d!=0)
      return d<0;
  }
  return false; // equality if we get here
}

template<typename C>
RationalVector<C> RationalVector<C>::operator+(const RationalVector<C>& v)
  const
{
  unsigned long gcd, m = arithmetic::lcm(d_denom,v.d_denom,gcd);
  unsigned long f = v.d_denom/gcd;
  assert (f==m/d_denom); // if this fails, then there was overflow on m
  RationalVector<C> result(d_num*f,m);
  result.d_num += v.d_num*(d_denom/gcd);
  return result; // don't normalize, better just limit denominator growth
}

// rational vectors are not guaranteed on lowest terms
// however if multiplication can be by cancellation it is done that way
template<typename C>
RationalVector<C>& RationalVector<C>::operator*=(C n)
{
  if (n!=0 and ((signed long)d_denom)%n==0) // avoid converting |n| to unsigned
    if (n>0)
      d_denom/=(unsigned long)n; // unsigned arithmetic OK here
    else
    {
      d_num=-d_num;
      d_denom/=(unsigned long)(-n);
    }
  else d_num*=n;
  return *this;
}

// division takes signed argument to avoid catastrophic surprises: if ever a
// negative argument were passed, implicit conversion to unsigned would be fatal
template<typename C>
RationalVector<C>& RationalVector<C>::operator/=(C n)
{
  if (n>0)
    d_denom*=(unsigned long)n;
  else
  {
    assert(n!=0);
    d_num=-d_num;
    d_denom*=(unsigned long)(-n);
  }
  return *this;
}

template<typename C>
RationalVector<C>& RationalVector<C>::normalize()
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

  /*

    Instantiation of templates (only these are generated)

  */

template class RationalVector<int>;           // the main instance used


} // |namespace ratvec|

} // |namespace atlas|
