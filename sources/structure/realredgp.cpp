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

/*****************************************************************************

  The function pause just serves as a convenient debugging break point

******************************************************************************/

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
  , d_connectivity()
  , d_status()
{
  const rootdata::RootDatum& rd = G_C.rootDatum();

  if (rd.rank() == rd.semisimpleRank())
    d_status.set(IsSemisimple);

  if (rf == G_C.quasisplit()) d_status.set(IsQuasisplit);

  // construct the torus for the most split Cartan

  const cartanclass::Fiber& fundf = d_complexGroup.fundamental();
  rootdata::RootList so= cartanclass::toMostSplit(fundf,rf,rd);

  latticetypes::LatticeMatrix q;
  rootdata::toMatrix(q,so,rd);
  q.leftMult(d_complexGroup.distinguished());

  tori::RealTorus T(q);

  if (T.isSplit())
    d_status.set(IsSplit);

  d_connectivity = topology::Connectivity(T,rd);
  if (d_connectivity.dualComponentReps().size() == 0)
    d_status.set(IsConnected);
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
