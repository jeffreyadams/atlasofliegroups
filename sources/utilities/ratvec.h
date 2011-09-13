/*!
\file
  This is ratvec.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef RATVEC_H  /* guard against multiple inclusions */
#define RATVEC_H

#include <vector>
#include <matrix.h>

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
  unsigned long d_denom; // a positive common denominator

 public:

// constructors and destructors
  explicit RationalVector(size_t r): d_num(r,0), d_denom(1){} // zero vector

  /*!
  Builds the RationalVector with numerator v and denominator d.
  */
  RationalVector(const V& v, C d);

// accessors
  // unsigned denominator requires care: plain % or / taboo; so export signed
  C denominator() const { return (C)d_denom; }
  const V& numerator() const { return d_num; }
  size_t size() const { return d_num.size(); }

  bool operator== (const RationalVector<C>& a) const;
  bool operator!= (const RationalVector<C>& a) const { return not operator==(a); }
  bool operator< (const RationalVector<C>& a) const; // comparison, for STL use

  RationalVector<C> operator+(const RationalVector<C>& v) const;
  RationalVector<C>& operator+=(const RationalVector<C>& v)
  { return *this=*this+v; }
  RationalVector<C> operator-() const
  { return RationalVector<C>(-d_num,d_denom); }
  RationalVector<C> operator-(const RationalVector<C>& v) const
  { return *this+-v; }
  RationalVector<C>& operator-=(const RationalVector<C>& v)
  { return *this=*this-v; }
  RationalVector<C>& operator*=(C n);
  RationalVector<C>& operator/=(C n);

/*
  Returns the scalar product of |*this| and |w|, which are assumed to be of
  same size and such that the scalar product is integral.

  A very long standing bug was to forget to cast |d_denom| to integer before
  the division. With that omission the scalar product is \emph{implicitly}
  cast to |unsigned int| instead, with desastrous consequences for the result.

  NOTE : this implementation does not worry about overflow. It is appropriate
  only for small denominators.
*/
  C scalarProduct(const V& w) const
  {
    return d_num.scalarProduct(w)/(C)d_denom;
  }

//manipulators
  RationalVector<C>& normalize();
  V& numerator() { return d_num; } // allow direct manipulation

}; // |template <typename C> class RationalVector|


} // |namespace ratvec|
} // |namespace atlas|
#endif
