/*!
\file
  This is intutils.cpp.
  This file contains some very simple templates for integer manipulations.
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

/*****************************************************************************

  This file contains some very simple templates for integer manipulations.

  The template parameter |I| should designate a signed integer type for |abs|
  and |divide| (otherwise these functions are redundant), for |factorial| it
  could be either signed or unsigned.

******************************************************************************/

#include "intutils.h"

namespace atlas {

namespace intutils {



/*!
  Synopsis: returns the factorial of a. It is assumed that a >= 0; we do
  not worry about overflow.
*/
template<typename I> I factorial(I n)
{
  I f = 1;

  for (I i =n; i>1; --i)
    f *= i;

  return f;
}


// There are no instantiations of |factorial| ! See however |size::factorial|


} // namespace intutils

} // namespace atlas
