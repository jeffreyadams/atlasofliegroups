/*!
\file
\brief Class definition and function declarations for CartanClassSet.
*/

/*
  This is cartanset.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef CARTANSET_H  /* guard against multiple inclusions */
#define CARTANSET_H

#include "cartanset_fwd.h" // get main class declaration and complexredgp_fwd
#include "rootdata_fwd.h"

#include "bitmap.h"
#include "constants.h"
#include "cartanclass.h"
#include "poset.h"
#include "realform.h"
#include "weyl.h"

#include "complexredgp.h" // now |ComplexReductiveGroup| is a complete type

/******** function declarations **********************************************/

namespace atlas {

namespace cartanset {

  unsigned long blockSize(realform::RealForm, realform::RealForm,
		  	  const CartanClassSet&);

  unsigned long kgbSize(realform::RealForm, const CartanClassSet&);

  void cayley_and_cross_part(rootdata::RootList& cayley,
			     weyl::WeylWord& cross,
			     const weyl::TwistedInvolution& tw,
			     const rootdata::RootDatum& rd,
			     const weyl::WeylGroup& W);
}

/******** type definitions ***************************************************/

namespace cartanset {

// inlined members of |CartanClassSet| that call methods of |d_parent|

inline const rootdata::RootDatum& CartanClassSet::rootDatum() const
  { return d_parent.rootDatum(); }
inline const weyl::WeylGroup& CartanClassSet::weylGroup() const
  { return d_parent.weylGroup(); }
inline const weyl::WeylElt
  CartanClassSet::canonicalize(weyl::TwistedInvolution& sigma) const
  { return canonicalize
      (sigma,
       bitset::RankFlags(constants::lMask[d_parent.semisimpleRank()]));
  }

inline void CartanClassSet::addCartan(weyl::TwistedInvolution tw)
{
  d_cartan.push_back(new cartanclass::CartanClass
		     (rootDatum(),involutionMatrix(tw)));
}




} // |namespace cartanset|

} // |namespace atlas|

#endif
