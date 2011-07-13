/*!
\file
\brief Forward declarations for the classes Subspace and
Subquotient.
*/
/*
  This is subquotient_fwd.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef SUBQUOTIENT_FWD_H  /* guard against multiple inclusions */
#define SUBQUOTIENT_FWD_H

#include <cstddef> // for |size_t|

#include "constants.h" // for |RANK_MAX|

/******** forward type declarations ******************************************/

namespace atlas {

namespace subquotient {

template<size_t dim> class Subspace;
template<size_t dim> class Subquotient;

//! \brief Subgroup of (Z/2Z)^RANK_MAX.
 typedef Subspace<constants::RANK_MAX> SmallSubspace;

//! \brief Subquotient of (Z/2Z)^RANK_MAX.
 typedef Subquotient<constants::RANK_MAX> SmallSubquotient;

}

}

#endif
