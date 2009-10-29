/*!
\file
\brief Implementation of the class RealReductiveGroup.
*/
/*
  This is realredgp.cpp
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "realredgp.h"

#include "cartanclass.h"
#include "complexredgp.h"
#include "rootdata.h"
#include "tori.h"

#include <cassert>

namespace atlas {

/*****************************************************************************

        Chapter I -- The RealReductiveGroup class

******************************************************************************/

namespace realredgp {


/*!
  Synopsis : constructs a real reductive group from the datum of a complex
  reductive group and a real form.
*/
RealReductiveGroup::RealReductiveGroup
  (complexredgp::ComplexReductiveGroup& G_C, realform::RealForm rf)
  : d_complexGroup(G_C)
  , d_realForm(rf)
  , d_connectivity() // wait for most split torus to be constructed below
  , d_Tg(tits::BasedTitsGroup(G_C,tits::square_class_grading_offset
			(G_C.fundamental(),square_class(),G_C.rootDatum())))
  , d_status()
{
  tori::RealTorus msT = G_C.cartan(G_C.mostSplit(rf)).fiber().torus();
  d_connectivity = topology::Connectivity(msT,G_C.rootDatum());

  d_status.set(IsConnected,d_connectivity.component_rank() == 0);
  d_status.set(IsCompact,msT.isCompact());

  d_status.set(IsQuasisplit,rf == G_C.quasisplit());
  d_status.set(IsSplit,msT.isSplit());
  d_status.set(IsSemisimple,G_C.rank() == G_C.semisimpleRank());

#ifndef NDEBUG
  // construct the torus for the most split Cartan
  const cartanclass::Fiber& fundf = G_C.fundamental();
  rootdata::RootSet so= cartanclass::toMostSplit(fundf,rf,G_C.rootSystem());

  // recompute matrix of most split Cartan
  const rootdata::RootDatum& rd = G_C.rootDatum();
  tori::RealTorus T1
    (refl_prod(so,rd) * G_C.distinguished()); // factors commute in fact

  topology::Connectivity c(T1,rd);
  assert(d_connectivity.component_rank() == c.component_rank());
#endif
}


/******** accessors *********************************************************/



/******** manipulators ******************************************************/


void RealReductiveGroup::swap(RealReductiveGroup& other)
{
  assert(&d_complexGroup==&other.d_complexGroup); // cannot swap references
  std::swap(d_realForm,other.d_realForm);
  d_connectivity.swap(other.d_connectivity);
  std::swap(d_status,other.d_status);
}

} // namespace realredgp

} // namespace atlas
