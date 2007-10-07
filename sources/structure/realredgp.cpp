/*!
\file
\brief Implementation of the class RealReductiveGroup.
*/
/*
  This is realredgp.cpp
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
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

namespace {

  void pause() {}

}

/*****************************************************************************

        Chapter I -- The RealReductiveGroup class

******************************************************************************/

namespace realredgp {

RealReductiveGroup::RealReductiveGroup(
		  complexredgp::ComplexReductiveGroup& G_C,
		  realform::RealForm rf)
  : d_complexGroup(&G_C)
  , d_realForm(rf)
  , d_connectivity()
  , d_status()

/*!
  Synopsis : constructs a real reductive group from the datum of a complex
  reductive group and a real form.
*/

{
  using namespace cartanclass;
  using namespace latticetypes;
  using namespace rootdata;
  using namespace tori;

  if (rootDatum().rank() == rootDatum().semisimpleRank())
    d_status.set(IsSemisimple);

  if (rf == G_C.quasisplit()) d_status.set(IsQuasisplit);

  // construct the torus for the most split Cartan

  const Fiber& fundf = d_complexGroup->fundamental();
  const RootDatum& rd = d_complexGroup->rootDatum();

  RootList so= toMostSplit(fundf,rf,rd);

  LatticeMatrix q;
  rootdata::toMatrix(q,so,rd);
  q.leftMult(d_complexGroup->distinguished());

  RealTorus T(q);

  if (T.isSplit())
    d_status.set(IsSplit);

  d_connectivity = topology::Connectivity(T,rootDatum());
  if (d_connectivity.dualComponentReps().size() == 0)
    d_status.set(IsConnected);
}

RealReductiveGroup::~RealReductiveGroup()

{}

/******** accessors *********************************************************/

const cartanclass::CartanClass& RealReductiveGroup::cartan(size_t cn) const

/*!
  Synopsis: returns cartan \#cn in the group.

  Precondition: cn belongs to cartanSet().

  NOTE: this is not inlined to avoid a dependency on complexredegp.h
*/

{
  return d_complexGroup->cartan(cn);
}

const poset::Poset& RealReductiveGroup::cartanOrdering() const

/*!
  Synopsis: returns the poset of Cartan classes.

  NOTE: the poset may be only partially constructed. The only guarantee is
  that the part of it for this up to the most split Cartan of this real
  form is valid.

  NOTE: this is not inlined to avoid a dependency on complexredegp.h
*/

{
  return d_complexGroup->cartanOrdering();
}

const bitmap::BitMap& RealReductiveGroup::cartanSet() const

/*!
  Synopsis: returns the support of the set of Cartan classes for this
  real form.

  NOTE: this is not inlined to avoid a dependency on complexredegp.h
*/

{
  return d_complexGroup->cartanSet(d_realForm);
}

const latticetypes::LatticeMatrix& RealReductiveGroup::distinguished() const

/*!
  Synopsis: returns the distinguished involution of the underlying complex
  group.

  NOTE : this is not inlined to avoid a compiling dependency on complexredgp.h
*/

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



size_t RealReductiveGroup::kgbSize() const

/*!
  Synopsis: returns the cardinality of K\\G/B.

  Precondition: fillCartan() has been called.
*/

{
  return d_complexGroup->kgbSize(d_realForm);
}

size_t RealReductiveGroup::mostSplit() const

/*!
  Synopsis: returns the most split cartan subgroup.

  Precondition: fillCartan() has been called.

  NOTE: this is not inlined to avoid a dependency on complexredgp.h
*/

{
  return d_complexGroup->mostSplit(d_realForm);
}

size_t RealReductiveGroup::numCartan() const

/*!
  Synopsis: returns the number of conjugacy classes of Cartan subgroups.

  NOTE: the value returned is correct only after fullCartan() has been
  called.

  NOTE: this is not inlined in order to avoid a dependency on bitmap.h
*/

{
  return cartanSet().size();
}

size_t RealReductiveGroup::rank() const

/*!
  Synopsis: returns the rank of the group.

  NOTE: this is not inlined to avoid a compiler dependency on rootdata.h
*/

{
  return rootDatum().rank();
}

const rootdata::RootDatum& RealReductiveGroup::rootDatum() const

/*!
  Synopsis: returns the root datum of the underlying complex group.

  NOTE: this is not inlined to avoid a compiler dependency on complexredgp.h
*/

{
  return d_complexGroup->rootDatum();
}

size_t RealReductiveGroup::semisimpleRank() const

/*!
  Synopsis: returns the semisimple rank of the group.

  NOTE: this is not inlined to avoid a compiler dependency on rootdata.h
*/

{
  return rootDatum().semisimpleRank();
}

const tits::TitsGroup& RealReductiveGroup::titsGroup() const

/*!
  Synopsis: returns a reference to the Tits group.

  NOTE: this is not inlined to avoid a compiler dependency on complexredgp.h
*/

{
  return d_complexGroup->titsGroup();
}

const weyl::WeylGroup& RealReductiveGroup::weylGroup() const

/*!
  Synopsis: returns the rank of the group.

  NOTE: this is not inlined to avoid a compiler dependency on complexredgp.h
*/

{
  return d_complexGroup->weylGroup();
}

/******** manipulators ******************************************************/

void RealReductiveGroup::fillCartan()

/*!
  Synopsis: makes representatives of conjugacy classes of Cartan subgroups
  for this group.

  More precisely, it makes sure that the CartanClasses structure of the
  underlying complex group contains all the cartans for this group.

  NOTE: may forward an Overflow exception.
*/

{
  using namespace cartanset;

  d_complexGroup->fillCartan(d_realForm);

  return;
}

void RealReductiveGroup::swap(RealReductiveGroup& other)

{
  std::swap(d_complexGroup,other.d_complexGroup);
  std::swap(d_realForm,other.d_realForm);
  d_connectivity.swap(other.d_connectivity);
  std::swap(d_status,other.d_status);

  return;
}

}

}
