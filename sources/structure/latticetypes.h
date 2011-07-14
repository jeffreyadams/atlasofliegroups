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
#include "arithmetic.h"
#include "bitset.h"

/* The following include is not required by this file but we put it here
   so that the types predeclared in "latticetypes_fwd.h" will be complete for
   any module that includes the current file, without it having to separately
   include "bitvector.h" */
#include "bitvector.h"

namespace atlas {

namespace latticetypes {

/******** function declarations **********************************************/


/******** type definitions ***************************************************/

  /*!
  \brief Element of lattice tensored with rational numbers.

  LatticeElt d_num divided by unsigned LatticeCoeff d_denom.
  */
class RatLatticeElt
{
  LatticeElt d_num;  // vector of integers, representing the numerators
  unsigned int d_denom; // a positive common denominator (LatticeCoeff is int)

 public:

// constructors and destructors
  explicit RatLatticeElt(size_t r): d_num(r,0), d_denom(1){} // zero vector

  /*!
  Builds the RatLatticeElt with numerator v and denominator d.
  */
  RatLatticeElt(const LatticeElt& v, LatticeCoeff d)
    :d_num(v), d_denom(arithmetic::abs(d))
  { if (d<0) d_num*=-1; }

// accessors
  // unsigned denominator requires care: plain % or / are taboo; export signed
  LatticeCoeff denominator() const { return (int)d_denom; }
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
  RatLatticeElt& operator*=(int n);
  RatLatticeElt& operator/=(int n);

/*
  Returns the scalar product of |*this| and |w|, which are assumed to be of
  same size and such that the scalar product is integral.

  A very long standing bug was to forget to cast |d_denom| to integer before
  the division. With that omission the scalar product is \emph{implicitly}
  cast to |unsigned int| instead, with desastrous consequences for the result.

  NOTE : this implementation does not worry about overflow. It is appropriate
  only for small denominators.
*/
  LatticeCoeff  scalarProduct(const LatticeElt& w) const
    {
      return d_num.scalarProduct(w)/(int)d_denom;
    }

//manipulators
  RatLatticeElt& normalize();
  LatticeElt& numerator() { return d_num; } // allow direct manipulation

}; // class RatLatticeElt

} // namespace latticetypes

} // namespace atlas

#endif
