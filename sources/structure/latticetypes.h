/*!
\file
\brief Declarations of classes and types for namespace latticetypes.
*/
/*
  This is latticetypes.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef LATTICETYPES_H  /* guard against multiple inclusions */
#define LATTICETYPES_H

#include "latticetypes_fwd.h"
#include "matrix.h" // to make |LatticeElt| a complete type
#include "intutils.h"

/* The following include is not required by this file but we put it here
   so that the types predeclared in "latticetypes_fwd.h" will be complete for
   any module that includes the current file, without it having to separately
   include "bitvector.h" */
#include "bitvector.h"

/******** function declarations **********************************************/

namespace atlas {

namespace latticetypes {

  inline
  BinaryEquation make_equation(const SmallBitVector lhs, bool rhs)
  {
    BinaryEquation eqn(lhs.data(),lhs.size());
    eqn.pushBack(rhs);
    return eqn;
  }

/******** type definitions ***************************************************/

  /*!
  \brief Element of lattice tensored with rational numbers.

  LatticeElt d_num divided by LatticeCoeff d_denom.
  */
class RatLatticeElt {

 private:

  /*!
  Vector of integers, representing the numerators of the RatLatticeElt.
  */
  LatticeElt d_num;

  /*!
  Integer, the common denominator of the RatLatticeElt.
  */
  unsigned int d_denom;

 public:

// constructors and destructors
  explicit RatLatticeElt(size_t r): d_num(r,0), d_denom(1){} // zero vector

  /*!
  Builds the RatLatticeElt with numerator v and denominator d.
  */
  RatLatticeElt(const LatticeElt& v, LatticeCoeff d)
    :d_num(v), d_denom(intutils::abs(d))
  { if (d<0) d_num*=-1; }

// accessors
  unsigned int denominator() const { return d_denom; }
  const LatticeElt& numerator() const { return d_num; }
  size_t size() const { return d_num.size(); }

  bool operator== (const RatLatticeElt& a) const;
  bool operator!= (const RatLatticeElt& a) const { return not operator==(a); }
  bool operator< (const RatLatticeElt& a) const; // comparison, for STL use

  RatLatticeElt operator+(const RatLatticeElt& v) const;
  RatLatticeElt& operator+=(const RatLatticeElt& v) { return *this=*this+v; }
  RatLatticeElt operator-() const { return RatLatticeElt(-d_num,d_denom); }
  RatLatticeElt operator-(const RatLatticeElt& v) const { return *this+-v; }
  RatLatticeElt& operator-=(const RatLatticeElt& v) { return *this=*this-v; }
  RatLatticeElt& operator*=(int n) { d_num*=n; return *this; }
  RatLatticeElt& operator/=(unsigned int n) { d_denom*=n; return *this; }

/*
  Returns the scalar product of |*this| and |w|, which are assumed to be of
  same size and such that the scalar product is integral.

  NOTE : this implementation does not worry about overflow. It is appropriate
  only for small denominators.
*/
  LatticeCoeff  scalarProduct(const LatticeElt& w) const
    {
      return d_num.scalarProduct(w)/d_denom;
    }

//manipulators
  RatLatticeElt& normalize();
  LatticeElt& numerator() { return d_num; } // allow direct manipulation
  unsigned int& denominator() { return d_denom; }

}; // class RatLatticeElt

} // namespace latticetypes

} // namespace atlas

#endif
