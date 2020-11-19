/*
  This is ratvec.cpp. This module contains some simple utilities for matrices.

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "ratvec.h"

#include <cassert>
#include <cstdlib> // for |std::abs|
#include <stdexcept>
#include <vector>

#include "arithmetic.h"
#include "bigint.h"

namespace atlas {

namespace ratvec {

//		      The RationalVector class template

template<typename C>
template<typename C1>
RationalVector<C>::RationalVector(const matrix::Vector<C1>& v, C d)
  : d_num(v.begin(),v.end()), d_denom(std::abs(d))
{ if (d<C(0))
    d_num.negate();
} // don't try to normalize, caller can do the explicitly if needed

template<typename C>
template<typename C1>
RationalVector<C>::RationalVector
  (const matrix::Vector<C1>& v, const arithmetic::Rational<C>& d)
    : d_num(v.begin(),v.end()), d_denom(d.true_denominator())
{
  d_num *= d.numerator();
} // don't try to normalize, caller can do the explicitly if needed

template<typename C>
RationalVector<C>::RationalVector(V&& v, C d)
  : d_num(std::move(v)), d_denom(std::abs(d))
{ if (d<C(0))
    d_num.negate();
} // don't try to normalize, caller can do the explicitly if needed



// unnormalised comparison without assuming |Numer_t| can hold cross products
template<typename C>
  bool RationalVector<C>::operator==(const RationalVector<C>& v) const
{ typedef arithmetic::big_int bigint;
  if (size()!=v.size())
    return false;
  bigint dv = bigint::from_unsigned(v.d_denom),
    dthis = bigint::from_unsigned(d_denom);
  for (size_t i=0; i<size(); ++i)
    if (bigint::from_signed(d_num[i])*dv !=
	bigint::from_signed(v.d_num[i])*dthis)
      return false;
  return true;
}

template<typename C>
  bool RationalVector<C>::operator<(const RationalVector<C>& v) const
{ // cross multiply component-wise, and compare
  typedef arithmetic::big_int bigint;
  if (size()!=v.size())
    return size()<v.size(); // compare sizes first
  bigint dv = bigint::from_unsigned(v.d_denom),
    dthis = bigint::from_unsigned(d_denom);
  for (size_t i=0; i<d_num.size(); ++i)
  {
    bigint d =
      bigint::from_signed(d_num[i])*dv - bigint::from_signed(v.d_num[i])*dthis;
    if (not d.is_zero())
      return d.is_negative();
  }
  return false; // equality if we get here
}

template<typename C>
  RationalVector<C>& RationalVector<C>::operator+=(const RationalVector<C>& v)
{
  assert(size()==v.size());
  arithmetic::Denom_t gcd, m = arithmetic::lcm(d_denom,v.d_denom,gcd);
  arithmetic::Denom_t f0 = v.d_denom/gcd, f1=d_denom/gcd;
  assert (f0==m/d_denom and f1==m/v.d_denom); // failure indicates overflow on m
  for (unsigned i=d_num.size(); i-->0; )
    d_num[i]=d_num[i]*f0 + v.d_num[i]*f1;
  d_denom=m;
  return normalize();
}

template<typename C>
  RationalVector<C>& RationalVector<C>::operator-=(const RationalVector<C>& v)
{
  assert(size()==v.size());
  arithmetic::Denom_t gcd, m = arithmetic::lcm(d_denom,v.d_denom,gcd);
  arithmetic::Denom_t f0 = v.d_denom/gcd, f1=d_denom/gcd;
  assert (f0==m/d_denom and f1==m/v.d_denom); // failure indicates overflow on m
  for (unsigned i=d_num.size(); i-->0; )
    d_num[i]=d_num[i]*f0 - v.d_num[i]*f1;
  d_denom=m;
  return normalize();
}

template<typename C>
  RationalVector<C>& RationalVector<C>::negate_add(const RationalVector<C>& v)
{
  assert(size()==v.size());
  arithmetic::Denom_t gcd, m = arithmetic::lcm(d_denom,v.d_denom,gcd);
  arithmetic::Denom_t f0 = v.d_denom/gcd, f1=d_denom/gcd;
  assert (f0==m/d_denom and f1==m/v.d_denom); // failure indicates overflow on m
  for (unsigned i=d_num.size(); i-->0; )
    d_num[i]=v.d_num[i]*f1 - d_num[i]*f0;
  d_denom=m;
  return normalize();
}

// rational vectors are not guaranteed on lowest terms
// however if multiplication can be done by cancellation, it is done that way
template<typename C>
RationalVector<C>& RationalVector<C>::operator*=(C n)
{
  if (n!=0 and
      arithmetic::Numer_t(d_denom)%n==0) // avoid converting |n| to unsigned
    if (n>0)
      d_denom/=arithmetic::Denom_t(n); // unsigned arithmetic OK here
    else
    {
      d_num=-d_num;
      d_denom/=arithmetic::Numer_t(-n);
    }
  else d_num*=n;
  return *this;
}

// division takes signed argument to avoid catastrophic surprises: if ever a
// negative argument were passed, implicit conversion to unsigned would be fatal
// like multiplication, this function does not guarantee a normalised result
template<typename C>
RationalVector<C>& RationalVector<C>::operator/=(C n)
{
  assert(n!=0);
  if (n>0)
    d_denom*=arithmetic::Denom_t(n);
  else
  {
    d_num=-d_num;
    d_denom*=-arithmetic::Denom_t(n); // safe to convert to unsigned, negate
  }
  return *this;
}

template<typename C>
RationalVector<C>& RationalVector<C>::operator%=(C n)
{
  assert(n!=0);
  C modulus(d_denom*std::abs(n));
  for (auto it=d_num.begin(); it!=d_num.end(); ++it)
    *it = arithmetic::remainder(*it,modulus);
  return *this;
}

template<typename C>
RationalVector<C>& RationalVector<C>::operator*=(const arithmetic::Rational<C>& r)
{ return (*this /= r.denominator()) *= r.numerator(); }

template<typename C>
RationalVector<C> RationalVector<C>::operator*(const arithmetic::Rational<C>& r)
const
{
  return RationalVector(d_num*r.numerator(),d_denom*r.denominator());
}

template<typename C>
RationalVector<C>& RationalVector<C>::operator/=(const arithmetic::Rational<C>& r)
{ assert (r.numerator()!=0);
  return (*this /= r.numerator()) *= r.denominator();
}

template<typename C>
RationalVector<C> RationalVector<C>::operator/(const arithmetic::Rational<C>& r)
const
{
  return RationalVector(d_num*r.denominator(),d_denom*r.numerator());
}

template<typename C>
const RationalVector<C>& RationalVector<C>::normalize() const
{
  arithmetic::Denom_t d=d_denom;

  if (d==0)
    throw std::runtime_error("Denominator 0 in rational vector");

  for (const auto& entry : d_num)
  {
    d=arithmetic::gcd(entry,d); // |gcd| takes (signed,unsigned) arguments
    if (d==1)
      return *this; // quit early when no common divisors are left
  }

  auto& my = const_cast<RationalVector<C>&>(*this); // now we shall "cheat"
  my.d_denom /= d;    // unsigned exact division, is safe here
  my.d_num   /= C(d); // signed exact division

  return *this; // export as |const| reference
}


template<typename C1, typename C2>
  RationalVector<C2> operator*
  (const matrix::Matrix<C1>& M, const RationalVector<C2>& v)
{
  return RationalVector<C2>(M*v.numerator(),v.denominator());
}

template<typename C1, typename C2>
  RationalVector<C2> operator*
  (const RationalVector<C2>& v,const matrix::Matrix<C1>& M)
{
  return RationalVector<C2>(M.right_prod(v.numerator()),v.denominator());
}

  /*

    Instantiation of templates (only these are generated)

  */

template class RationalVector<arithmetic::Numer_t>; // the main instance used
template RationalVector<arithmetic::Numer_t>::RationalVector
  (const matrix::Vector<int>&, arithmetic::Numer_t);
template RationalVector<arithmetic::Numer_t>::RationalVector
  (const matrix::Vector<int>&, const arithmetic::Rational<arithmetic::Numer_t>&);
template RationalVector<arithmetic::Numer_t>::RationalVector
  (const matrix::Vector<arithmetic::Numer_t>&, arithmetic::Numer_t);
template RationalVector<arithmetic::Numer_t> operator*
  (const matrix::Matrix<int>& M, const RationalVector<arithmetic::Numer_t>& v);
template RationalVector<arithmetic::Numer_t> operator*
  (const RationalVector<arithmetic::Numer_t>& v,const matrix::Matrix<int>& M);

} // |namespace ratvec|

} // |namespace atlas|
