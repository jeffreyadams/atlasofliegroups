/*
  This is ratvec.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2020 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef RATVEC_H  /* guard against multiple inclusions */
#define RATVEC_H

#include "ratvec_fwd.h" // ensure coherence

#include <vector>
#include <cassert>

#include "matrix.h"
#include "arithmetic.h"

// extra defs for windows compilation -spc
#ifdef WIN32
#include "constants.h"
#endif

namespace atlas {

namespace ratvec {

//			       type definitions

/* Element of lattice tensored with rational numbers.

   |LatticeElt d_num| divided by unsigned |LatticeCoeff d_denom|.
*/
template <typename C> // a signed integral type
class RationalVector
{
  typedef matrix::Vector<C> V; // local abbreviation

  V d_num;  // vector of integers, representing the numerators
  arithmetic::Denom_t d_denom; // a positive common denominator

 public:

// constructors and destructors
  explicit RationalVector(size_t r): d_num(r,0), d_denom(1){} // zero vector

  // build the RationalVector with numerator v and denominator d.
  template <typename C1> // possibly convert |v| entries to other type
    RationalVector(const matrix::Vector<C1>& v, C d);

  RationalVector(V&& v, C d);

  RationalVector() : d_num(),d_denom(C(1)) {} // default to empty numerator

// pseudo accessor
  // this method although classified |const| modifies members, but in a manner
  // that is mathematically neutral, so having this done should never harm
  const RationalVector& normalize() const;

// accessors
  // unsigned denominator requires care: plain % or / taboo; so export signed
  C denominator() const { return (C)d_denom; }
  const V& numerator() const { return d_num; }
  size_t size() const { return d_num.size(); }
  bool isZero() const { return d_num.isZero(); }

  bool operator== (const RationalVector& a) const;
  bool operator!= (const RationalVector& a) const
  { return not operator==(a); }
  bool operator< (const RationalVector& a) const; // comparison, for STL use

  RationalVector& operator+=(const RationalVector& v);
  RationalVector& operator-=(const RationalVector& v);
  RationalVector& negate() { d_num.negate(); return *this; }
  RationalVector& negate_add(const RationalVector& v);
  RationalVector& operator*=(C n);
  RationalVector& operator/=(C n);
  RationalVector& operator%=(C n);
  // functional versions are free functions taking first agument by value

  RationalVector operator+(const RationalVector& v) const &
  { RationalVector result(*this); return result+=v; }
  RationalVector operator-(const RationalVector& v) const &
  { RationalVector result(*this); return result-=v; }
  RationalVector operator-() const &
  { RationalVector result(*this); return result.negate(); }

  // when left operand is rvalue reference, use destructive operators
  RationalVector operator+ (const RationalVector& v) && { return *this += v; }
  RationalVector operator- (const RationalVector& v) && { return *this -= v; }
  RationalVector operator- () &&     { return negate(); }

  // when right operand is rvalue reference, use destructive operators too
  RationalVector operator+ (RationalVector&& v) const & { return v+=*this; }
  RationalVector operator- (RationalVector&& v) const &
					     { return v.negate_add(*this); }

  // when both operands are rvalue references, give priority to the left
  RationalVector operator+ (RationalVector&& v) && { return *this += v; }
  RationalVector operator- (RationalVector&& v) && { return *this -= v; }

  template <typename C1>
    RationalVector& operator+=(const matrix::Vector<C1>& v)
    { d_num.add(v.begin(),denominator()); return *this; }
  template <typename C1>
    RationalVector& operator-=(const matrix::Vector<C1>& v)
    { d_num.subtract(v.begin(),denominator()); return *this; }

  template <typename C1>
    RationalVector operator+(const matrix::Vector<C1>& v) const &
    { return RationalVector(*this)+=v; }
  template <typename C1>
    RationalVector operator-(const matrix::Vector<C1>& v) const &
    { return RationalVector(*this)-= v; }
  template <typename C1>
    RationalVector operator+(const matrix::Vector<C1>& v) &&
    { return (*this)+=v; }
  template <typename C1>
    RationalVector operator-(const matrix::Vector<C1>& v) &&
    { return (*this)-= v; }

  RationalVector& operator*=(const arithmetic::Rational<C>& r);
  RationalVector& operator/=(const arithmetic::Rational<C>& r);
  RationalVector operator*(const arithmetic::Rational<C>& r) const;
  RationalVector operator/(const arithmetic::Rational<C>& r) const;

/*
  Returns the scalar product of |*this| and |w|, which are assumed to be of
  same size and such that the scalar product is integral.

  A very long standing bug was to forget to cast |d_denom| to integer before
  the division. With that omission the scalar product is \emph{implicitly}
  cast to |unsigned int| instead, with disastrous consequences for the result.

  NOTE : this implementation does not worry about overflow. It is appropriate
  only for small denominators.
*/
  template <typename C1>
    C dot(const matrix::Vector<C1>& w) const
  {
    auto num = w.dot(d_num);
    assert(num%(C)d_denom==0); // division below must be without remainder
    return num/(C)d_denom; // order is imposed here by return type |C|
  }

// take difference as integer vector (which it is assumed to be here), converting
// entries (without any test) to a possibly different signed integer type |C1|
  template<typename C1>
    matrix::Vector<C1> integer_diff(const RationalVector<C>& v) const
  {
    assert(size()==v.size());
    normalize(); v.normalize();
    assert(d_denom==v.d_denom); // integer difference, so equal normalized denom
    auto d = denominator(); // convert to signed type (avoid unsigned division!)
    matrix::Vector<C1> result(size());
    for (unsigned i=0; i<size(); ++i)
    { assert((d_num[i]-v.d_num[i])%d==0);
      result[i] = (d_num[i]-v.d_num[i])/d;
    }
    return result;
  }

//manipulators
  RationalVector& normalize()
  { static_cast<const RationalVector*>(this) -> normalize(); return *this; }
  V& numerator() { return d_num; } // allow direct manipulation

}; // |template <typename C> class RationalVector|




//				Functions


// right-operate by scalars
template<typename C>
  RationalVector<C> operator* (RationalVector<C> v, C n) { return v*=n ;}
template<typename C>
  RationalVector<C> operator/ (RationalVector<C> v, C n) { return v/=n ;}
template<typename C>
  RationalVector<C> operator% (RationalVector<C> v, C n) { return v%=n ;}

// left-multiply by a matrix
template<typename C1, typename C2>
  RationalVector<C2> operator*
  (const matrix::Matrix<C1>& M, const RationalVector<C2>& v);

// right-multiply by a matrix
template<typename C1, typename C2>
  RationalVector<C2> operator*
  ( const RationalVector<C2>& v,const matrix::Matrix<C1>& M);

// project to fixed points of involution
template<typename C1, typename C2>
  RationalVector<C2>& symmetrise
  (const matrix::Matrix<C1>& M,RationalVector<C2>& v)
{
  v.numerator() += M*v.numerator();
  return (v/=2).normalize();
}

// project to dual fixed points of involution
template<typename C1, typename C2>
  RationalVector<C2>& symmetrise
  (RationalVector<C2>& v,const matrix::Matrix<C1>& M)
{
  v.numerator() += M.right_prod(v.numerator());
  return (v/=2).normalize();
}



} // |namespace ratvec|
} // |namespace atlas|
#endif
