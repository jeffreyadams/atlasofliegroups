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
RealReductiveGroup::RealReductiveGroup(
		  complexredgp::ComplexReductiveGroup& G_C,
		  realform::RealForm rf)
  : d_complexGroup(&G_C)
  , d_realForm(rf)
  , d_connectivity()
  , d_status()
{
  const rootdata::RootDatum& rd = G_C.rootDatum();

  if (rd.rank() == rd.semisimpleRank())
    d_status.set(IsSemisimple);

  if (rf == G_C.quasisplit()) d_status.set(IsQuasisplit);

  // construct the torus for the most split Cartan

  const cartanclass::Fiber& fundf = d_complexGroup->fundamental();
  rootdata::RootList so= cartanclass::toMostSplit(fundf,rf,rd);

  latticetypes::LatticeMatrix q;
  rootdata::toMatrix(q,so,rd);
  q.leftMult(d_complexGroup->distinguished());

  tori::RealTorus T(q);

  if (T.isSplit())
    d_status.set(IsSplit);

  d_connectivity = topology::Connectivity(T,rd);
  if (d_connectivity.dualComponentReps().size() == 0)
    d_status.set(IsConnected);
}


/******** accessors *********************************************************/


/*!
  Synopsis: returns cartan \#cn in the group.

  Precondition: cn belongs to cartanSet().

  NOTE: this is not inlined to avoid a dependency on complexredegp.h
*/
const cartanclass::CartanClass& RealReductiveGroup::cartan(size_t cn) const
{
  return d_complexGroup->cartan(cn);
}



/*!
  Synopsis: returns the support of the set of Cartan classes for this
  real form.

  NOTE: this is not inlined to avoid a dependency on complexredegp.h
*/
const bitmap::BitMap& RealReductiveGroup::Cartan_set() const
{
  return complexGroup().Cartan_set(d_realForm);
}

size_t RealReductiveGroup::numInvolutions() const {
    return complexGroup().numInvolutions(Cartan_set());
  }


/*!
  Synopsis: returns the distinguished involution of the underlying complex
  group.

  NOTE : this is not inlined to avoid a compiling dependency on complexredgp.h
*/
const latticetypes::LatticeMatrix& RealReductiveGroup::distinguished() const
{
  return d_complexGroup->distinguished();
}

/*!
  Synopsis: Returns the set of noncompact imaginary roots for (the
  representative of) the real form.

  NOTE : this is not inlined to avoid a compiling dependency on complexredgp.h
*/
rootdata::RootSet RealReductiveGroup::noncompactRoots() const
{
  return d_complexGroup->noncompactRoots(d_realForm);
}

/*!
  Synopsis: returns the cardinality of K\\G/B.

  Precondition: fillCartan() has been called.
*/
size_t RealReductiveGroup::KGB_size() const
{
  return d_complexGroup->KGB_size(d_realForm);
}


/*!
  Synopsis: returns the most split cartan subgroup.

  Precondition: fillCartan() has been called.

  NOTE: this is not inlined to avoid a dependency on complexredgp.h
*/
size_t RealReductiveGroup::mostSplit() const
{
  return d_complexGroup->mostSplit(d_realForm);
}

/*!
  Synopsis: returns the number of conjugacy classes of Cartan subgroups.

  NOTE: the value returned is correct only after fullCartan() has been
  called.

  NOTE: this is not inlined in order to avoid a dependency on bitmap.h
*/
size_t RealReductiveGroup::numCartan() const
{
  return Cartan_set().size();
}


/*!
  Synopsis: returns the rank of the group.

  NOTE: this is not inlined to avoid a compiler dependency on rootdata.h
*/
size_t RealReductiveGroup::rank() const
{
  return rootDatum().rank();
}


/*!
  Synopsis: returns the root datum of the underlying complex group.

  NOTE: this is not inlined to avoid a compiler dependency on complexredgp.h
*/
const rootdata::RootDatum& RealReductiveGroup::rootDatum() const
{
  return d_complexGroup->rootDatum();
}


/*!
  Synopsis: returns the semisimple rank of the group.

  NOTE: this is not inlined to avoid a compiler dependency on rootdata.h
*/
size_t RealReductiveGroup::semisimpleRank() const
{
  return rootDatum().semisimpleRank();
}


/*!
  Synopsis: returns a reference to the Tits group.

  NOTE: this is not inlined to avoid a compiler dependency on complexredgp.h
*/
const tits::TitsGroup& RealReductiveGroup::titsGroup() const
{
  return d_complexGroup->titsGroup();
}


/*!
  Synopsis: returns the Weyl group of the associated complex group.

  NOTE: this is not inlined to avoid a compiler dependency on complexredgp.h
*/

const weyl::WeylGroup& RealReductiveGroup::weylGroup() const
{
  return d_complexGroup->weylGroup();
}

/******** manipulators ******************************************************/


/*!
  Synopsis: makes representatives of conjugacy classes of Cartan subgroups
  for this group.

  More precisely, it makes sure that the CartanClasses structure of the
  underlying complex group contains all the cartans for this group.

  This function can be called for any real form, but the atlas code has been
  rewritten to call it only for the quasisplit form, so that the numbering of
  the Cartan classes for an inner class is independent of the order in which
  the various real forms are first considered (and their Cartans generated)
*/
void RealReductiveGroup::fillCartan()
{
  d_complexGroup->fillCartan(d_realForm);
}

void RealReductiveGroup::swap(RealReductiveGroup& other)
{
  std::swap(d_complexGroup,other.d_complexGroup);
  std::swap(d_realForm,other.d_realForm);
  d_connectivity.swap(other.d_connectivity);
  std::swap(d_status,other.d_status);
}

} // namespace realredgp

} // namespace atlas
