/*
  This is latticetypes.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "latticetypes.h"

/*****************************************************************************

  This module defines a few simple operations on the types defined in 
  latticetypes.h

******************************************************************************/

/*****************************************************************************

        Chapter I : the LatticeElt class

  A LatticeElt is just a vector of some signed integral type, which will
  be fixed in the whole of the program.

******************************************************************************/

namespace atlas {

namespace latticetypes {

LatticeElt& operator+= (LatticeElt& v, const LatticeElt& w)

/*
  Increments v by w.
*/

{  
  for (LatticeElt::size_type j = 0; j < v.size(); ++j)
    v[j] += w[j];

  return v;
}

LatticeElt& operator-= (LatticeElt& v, const LatticeElt& w)

/*
  Decrements v by w.
*/

{
  for (LatticeElt::size_type j = 0; j < v.size(); ++j)
    v[j] -= w[j];

  return v;
}

LatticeElt& operator*= (LatticeElt& v, const LatticeCoeff& c)

/*
  Multiplies each coordinate of v by c.
*/

{  
  for (LatticeElt::size_type j = 0; j < v.size(); ++j)
    v[j] *= c;

  return v;
}

LatticeElt& operator/= (LatticeElt& v, const LatticeCoeff& d)

/*
  Divides each coordinate of v by d. It is the callers responsibility to
  check that v is actually divisible by d; if not, he will just get the
  integral parts of the quotients.
*/

{  
  for (LatticeElt::size_type j = 0; j < v.size(); ++j)
    v[j] /= d;

  return v;
}

LatticeElt& operator- (LatticeElt& v)

/*
  Changes the sign of each coordinate of v.
*/

{  
  for (LatticeElt::size_type j = 0; j < v.size(); ++j)
    v[j] = -v[j];

  return v;
}

}

/*****************************************************************************

        Chapter II -- Functions declared in latticetypes.h

  This section contains the definition of the functions declared in
  latticetypes.h :

    - scalarProduct(const LatticeElt&, const LatticeElt&);

******************************************************************************/

namespace latticetypes {

bool isZero(const LatticeElt& v)

/*
  Returns true if all components of v are zero, false otherwise.
*/

{
  for (size_t j = 0; j < v.size(); ++j)
    if (v[j])
      return false;

  return true;
}

LatticeCoeff scalarProduct(const LatticeElt& v, const LatticeElt& w)

/*
  Returns the scalar product of v and w, which are assumed to be of same
  size.
*/

{
  LatticeCoeff sp = 0;

  for (size_t j = 0; j < v.size(); ++j) {
    sp += v[j]*w[j];
  }

  return sp;
}

LatticeCoeff scalarProduct(const RatLatticeElt& v, const LatticeElt& w)

/*
  Returns the scalar product of v and w, which are assumed to be of same
  size and such that the scalar product is integral.

  NOTE : this implementation does not worry about overflow. It is appropriate
  only for small denominators.
*/

{
  LatticeCoeff sp = 0;
  const LatticeElt& vnum = v.numerator();

  for (size_t j = 0; j < w.size(); ++j) {
    sp += vnum[j]*w[j];
  }

  sp /= v.denominator();

  return sp;
}

}

}
