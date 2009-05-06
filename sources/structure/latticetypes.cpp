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

  RatLatticeElt& RatLatticeElt::normalize()
  {
    unsigned long d=std::abs(d_denom);

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
    for (size_t i=0; i<d_num.size(); ++i)
      d_num/=d;
    return *this;
  }

} // namespace latticetypes

} // namespace atlas
