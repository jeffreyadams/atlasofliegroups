/*!
\file
  \brief This module defines a few simple operations on the types defined in
  latticetypes.h.
*/
/*
  This is latticetypes.cpp.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "latticetypes.h"

/*****************************************************************************

  This module defines a few simple operations on the types defined in
  latticetypes.h

******************************************************************************/

/*****************************************************************************

        Chapter I -- Functions declared in latticetypes.h

  This section contains the definition of the function declared in
  latticetypes.h :

    - isZero(const LatticeElt& v)

******************************************************************************/

namespace atlas {

namespace latticetypes {


/*
  Returns true if all components of v are zero, false otherwise.
*/
bool isZero(const LatticeElt& v)
{
  for (size_t j = 0; j < v.size(); ++j)
    if (v[j])
      return false;

  return true;
}

} // namespace latticetypes

} // namespace atlas
