/*!
\file
  This is ratvec.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef RATVEC_H  /* guard against multiple inclusions */
#define RATVEC_H

#include "ratvec_fwd.h" // ensure coherence

#include <vector>
#include "matrix.h"
#include "arithmetic.h"

// extra defs for windows compilation -spc
#ifdef WIN32
#include "constants.h"
#endif

namespace atlas {

namespace ratvec {

//			       type definitions

  /*!
  \brief Element of lattice tensored with rational numbers.

  LatticeElt d_num divided by unsigned LatticeCoeff d_denom.
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

  // Build the RationalVector with numerator v and denominator d.
  template <typename C1>
    RationalVector(const matrix::Vector<C1>& v, C d);

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

  RationalVector operator+(const RationalVector& v) const;
  RationalVector& operator+=(const RationalVector& v)
  { return *this=*this+v; }
  RationalVector operator-() const
  { return RationalVector(-d_num,d_denom); }
  RationalVector operator-(const RationalVector& v) const
  { return *this+-v; }
  RationalVector& operator-=(const RationalVector& v)
  { return *this=*this-v; }
  RationalVector& operator*=(C n);
  RationalVector& operator/=(C n);
  RationalVector& operator%=(C n);

  template <typename C1>
    RationalVector& operator+=(const matrix::Vector<C1>& v)
    { d_num.add(v.begin(),denominator()); return *this; }
  template <typename C1>
    RationalVector& operator-=(const matrix::Vector<C1>& v)
    { d_num.subtract(v.begin(),denominator()); return *this; }

  template <typename C1>
    RationalVector operator+(const matrix::Vector<C1>& v) const
    { return RationalVector(*this)+=v; }
  template <typename C1>
    RationalVector operator-(const matrix::Vector<C1>& v) const
    { return RationalVector(*this)-= v; }

  RationalVector& operator*=(const arithmetic::Rational& r);
  RationalVector& operator/=(const arithmetic::Rational& r);
  RationalVector operator*(const arithmetic::Rational& r) const;
  RationalVector operator/(const arithmetic::Rational& r) const;

/*
  Returns the scalar product of |*this| and |w|, which are assumed to be of
  same size and such that the scalar product is integral.

  A very long standing bug was to forget to cast |d_denom| to integer before
  the division. With that omission the scalar product is \emph{implicitly}
  cast to |unsigned int| instead, with desastrous consequences for the result.

  NOTE : this implementation does not worry about overflow. It is appropriate
  only for small denominators.
*/
  template <typename C1>
    C dot(const matrix::Vector<C1>& w) const
  {
    return w.dot(d_num)/(C)d_denom; // order is imposed here by return type |C|
  }

//manipulators
  RationalVector& normalize();
  V& numerator() { return d_num; } // allow direct manipulation

}; // |template <typename C> class RationalVector|




//				Functions

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
