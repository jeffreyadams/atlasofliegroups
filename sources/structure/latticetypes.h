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

/* The following includes are not required by this file but we put them here
   so that the types predeclared in "latticetypes_fwd.h" will be complete for
   any module that includes the current file, without it having to separately
   include "bitvector.h" and "matrix.h" */
#include "bitvector.h"
#include "matrix.h"

/******** function declarations **********************************************/

namespace atlas {

namespace latticetypes {

  LatticeElt& operator+= (LatticeElt&, const LatticeElt&);
  LatticeElt& operator-= (LatticeElt&, const LatticeElt&);
  LatticeElt& operator*= (LatticeElt&, LatticeCoeff);
  LatticeElt& operator/= (LatticeElt&, LatticeCoeff);
  LatticeElt& operator- (LatticeElt&); // negates argument in place

  bool isZero(const LatticeElt&);

  LatticeCoeff scalarProduct(const LatticeElt&, const LatticeElt&);
  LatticeCoeff scalarProduct(const RatLatticeElt&, const LatticeElt&);

}

/******** type definitions ***************************************************/

namespace latticetypes {

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
  LatticeCoeff d_denom;

 public:

// constructors and destructors
  RatLatticeElt()
    {}

  /*!
  Builds the RatLatticeElt with numerator v and denominator d.
  */
  RatLatticeElt(const LatticeElt& v, LatticeCoeff d)
    :d_num(v), d_denom(d)
    {}

  /*!
  Builds a RatLatticeElt of in Z^n with denominator d and all entries zero.
  */
  RatLatticeElt(size_t n, LatticeCoeff d)
    :d_num(n), d_denom(d)
    {}

  /*!
  Copies the RatLatticeElt v into  a new RatLatticeElt.
  */
  RatLatticeElt(const RatLatticeElt& v)
    :d_num(v.d_num), d_denom(v.d_denom)
    {}

  ~RatLatticeElt()
    {}

// accessors
  LatticeCoeff denominator() const {
    return d_denom;
  }

  const LatticeElt& numerator() const {
    return d_num;
  }

  size_t size() const {
    return d_num.size();
  }

//manipulators
  LatticeCoeff& denominator() {
    return d_denom;
  }

  LatticeElt& numerator() {
    return d_num;
  }

};

}

}

#endif
