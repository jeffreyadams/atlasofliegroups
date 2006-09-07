/*!
\file
\brief Implementation for the class ComplexReductiveGroup.

  The ComplexReductiveGroup class will play a central role in the
  whole program.  Even though it is entirely defined by its based root
  datum and an involutive automorphism of that datum, it has seemed
  more natural to use this class to collect the wealth of
  combinatorial data that the root datum gives rise to, and that will
  serve as the basis for our description of the representation theory
  of G. Note that the current state of the theory, and most notably
  Vogan duality, makes it natural and necessary to consider all the
  real forms of our complex group (in a given inner class) at once; so
  that is another reason to not choose a real form a priori.
*/
/*
  This is complexredgp.cpp.

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "complexredgp.h"

#include <set>

#include "cartan.h"
#include "dynkin.h"
#include "lietype.h"
#include "poset.h"
#include "rootdata.h"
#include "tits.h"
#include "weyl.h"

/*****************************************************************************

  The ComplexReductiveGroup class will play a central role in the whole
  program; even though it is entirely defined by its root datum, it has
  seemed more natural to use this class to collect the wealth of
  combinatorial data that the root datum gives rise to, and that will
  serve as the basis for our description of the representation theory
  of G. Note that the current state of the theory, and most notably
  Vogan duality, makes it natural and necessary to consider all the
  real forms of our complex group (in a given inner class) at once;
  so that is another reason to not choose a real form a priori.

******************************************************************************/

namespace atlas {

namespace {

  void pause() {;}

}

/*****************************************************************************

        Chapter I -- The ComplexReductiveGroup class

  ... explain here when it is stable ...

******************************************************************************/

namespace complexredgp {

ComplexReductiveGroup::ComplexReductiveGroup()
  :d_rootDatum(0),
   d_titsGroup(0),
   d_cartan(0)

/*!
  Synopsis: default constructor.

  This is of course unusable, but it might be created somewhere (for
  instance before an interactive allocation that then fails), and would
  bomb on destruction if the pointers were not set to zero.
*/

{}

ComplexReductiveGroup::
ComplexReductiveGroup(const rootdata::RootDatum* rd,
		      const latticetypes::LatticeMatrix& d)
  :d_rootDatum(rd),
   d_titsGroup(0),
   d_cartan(0)

/*!
  Synopsis: constructs a ComplexReductiveGroup from its root datum and a
  distinguished involution.

  Precondition: d is an involution which fixes a pinning of the group; in
  particular it fixes the Borel.

  NOTE: the ComplexReductiveGroup assumes ownership of the RootDatum pointed
  to by rd. Perhaps this is not entirely safe.
*/

{
  using namespace cartan;
  using namespace latticetypes;
  using namespace rootdata;
  using namespace tits;
  using namespace weyl;

  LatticeMatrix q;
  cartanMatrix(q,rootDatum());
  d_titsGroup = new TitsGroup(rootDatum(),d);
  d_cartan = new CartanClasses(rootDatum(),d,weylGroup());
}

ComplexReductiveGroup::ComplexReductiveGroup(const ComplexReductiveGroup& G, 
					     tags::DualTag)

/*!
  Synopsis: constructs the complex reductive group dual to G.
*/

{
  using namespace cartan;
  using namespace latticetypes;
  using namespace rootdata;
  using namespace tags;
  using namespace tits;

  LatticeMatrix d;
  dualBasedInvolution(d,G.distinguished(),G.rootDatum());

  d_rootDatum = new RootDatum(G.rootDatum(),DualTag());
  d_titsGroup = new TitsGroup(rootDatum(),d);
  d_cartan = new CartanClasses(rootDatum(),d,weylGroup());
}

ComplexReductiveGroup::~ComplexReductiveGroup()

{
  delete d_cartan;
  delete d_titsGroup;
  delete d_rootDatum;
}

/******** accessors **********************************************************/

unsigned long ComplexReductiveGroup::blockSize(realform::RealForm rf, 
					       realform::RealForm drf) const

/*!
  Synopsis: returns the size of the block defined by the weak real form rf
  and the weak dual real form drf.

  NOTE: this is just a forwarding function. It is not inlined to avoid a
  compilation dependency, and because it is expected to be used very little.
*/

{
  return cartan::blockSize(rf,drf,*d_cartan);
}

const cartanclass::CartanClass& ComplexReductiveGroup::cartan(size_t cn) const

/*!
  Synopsis: returns cartan \#cn in the group.

  NOTE: this is not inlined to avoid a dependency on cartan.h
*/

{
  return d_cartan->cartan(cn);
}

const poset::Poset& ComplexReductiveGroup::cartanOrdering() const

/*!
  Synopsis: returns the ordering of the Cartan subgroups

  NOTE: this is not inlined to avoid a dependency on cartan.h
*/

{
  return d_cartan->ordering();
}

const bitmap::BitMap& ComplexReductiveGroup::cartanSet(realform::RealForm rf) 
  const

/*!
  Synopsis: returns the support of the set of Cartan classes for rf

  NOTE : this is not inlined to avoid a compiling dependency on cartan.h
*/

{
  return d_cartan->support(rf);
}

const latticetypes::LatticeMatrix& ComplexReductiveGroup::distinguished() const

/*!
  Synopsis: returns the matrix of the distinguished involution.

  Recall that we always implicitly fix an inner class for our complex reductive
  group.

  NOTE: this is not inlined to avoid a dependency on cartan.h
*/

{
  return d_cartan->distinguished();
}

unsigned long ComplexReductiveGroup::dualFiberSize(realform::RealForm drf, 
						   size_t cn) 
  const

/*!
  Synopsis: returns the size of the fiber size corresponding to dual real
  form \#drf and cartan \#cn.

  Explanation: this is the size of the orbits, for the shifted action of
  the dual imaginary Weyl group (a.k.a. the real Weyl group), that correspond
  to drf in the classification of strong real forms (they all have the same
  size.)

  This is a technical function used for size computations.

  NOTE: this is not inlined to avoid a dependency on cartan.h.
*/

{
  return d_cartan->dualFiberSize(drf,cn);
}

const cartanclass::Fiber& ComplexReductiveGroup::dualFundamental() const

/*!
  Synopsis: returns the dual fundamental fiber.

  This is a technical data structure containing the data for the classification
  of weak and strong dual real forms of G.

  NOTE: this is not inlined to avoid a dependency on cartan.h.
*/

{
  return d_cartan->dualFundamental();
}

const realform::RealFormList& 
ComplexReductiveGroup::dualRealFormLabels(size_t cn) const

/*!
  Synopsis: returns the dual real form labels for cartan \#cn

  NOTE: this is not inlined to avoid a dependency on cartan.h.
*/

{
  return d_cartan->dualRealFormLabels(cn);
}

unsigned long ComplexReductiveGroup::dualRepresentative(realform::RealForm drf, 
							size_t cn)
  const

/*!
  Synopsis: returns an element of the orbit corresponding to drf in the
  classification of weak real forms for cartan \#cn.

  Precondition: cartan \#cn is defined for drf.

  NOTE: this is not inlined to avoid a dependency on cartan.h
*/

{
  return d_cartan->dualRepresentative(drf,cn);
}

unsigned long ComplexReductiveGroup::fiberSize(realform::RealForm rf, size_t cn) 
  const

/*!
  Synopsis: returns the size of the fiber size corresponding to real
  form \#rf and cartan \#cn.

  Explanation: this is the size of the orbits, for the shifted action of
  the imaginary Weyl group, that correspond to rf in the classification of 
  strong real forms (they all have the same size.)

  This is a technical function used for size computations.

  NOTE: this is not inlined to avoid a dependency on cartan.h.
*/

{
  return d_cartan->fiberSize(rf,cn);
}

const cartanclass::Fiber& ComplexReductiveGroup::fundamental() const

/*!
  Synopsis: returns the fundamental fiber.

  This is a technical data structure containing the data for the classification
  of weak and strong real forms of G.

  NOTE: this is not inlined to avoid a dependency on cartan.h.
*/

{
  return d_cartan->fundamental();
}

void ComplexReductiveGroup::grading(rootdata::RootSet& rs, 
				    realform::RealForm rf) 
  const

/*!
  Synopsis: puts in rs the set of noncompact imaginary roots for the
  representative of rf
*/

{
  d_cartan->noncompactRootSet(rs,rf);

  return;
}

unsigned long ComplexReductiveGroup::kgbSize(realform::RealForm rf) const

/*!
  Synopsis: returns the number of elements in K\\G/B for real form \#rf.

  Explanation: this is exactly the number of elements in the one-sided
  parameter set corresponding to any strong real form of G lying over rf.

  NOTE: this is not inlined to avoid a dependency on cartan.h.
*/

{
  return cartan::kgbSize(rf,*d_cartan);
}

size_t ComplexReductiveGroup::numCartanClasses() const

/*!
  Synopsis: returns the number of conjugacy classes of Cartan subgroups
  currently constructed for G.

  NOTE: this is not inlined to avoid a dependency on cartan.h.
*/

{
  return d_cartan->numCartan();
}

size_t ComplexReductiveGroup::mostSplit(realform::RealForm rf) const

/*!
  Synopsis: returns the most split cartan subgroup for real form \#rf.

  Precondition: fillCartan() has been called for rf.

  NOTE: this is not inlined to avoid a dependency on cartan.h
*/

{
  return d_cartan->mostSplit(rf);
}

size_t ComplexReductiveGroup::numDualRealForms() const

/*!
  Synopsis: returns the number of weak dual real forms for this inner class.

  NOTE: this is not inlined to avoid a dependency on cartan.h
*/

{
  return d_cartan->numDualRealForms();
}

size_t ComplexReductiveGroup::numInvolutions() const

/*!
  Synopsis: returns the number of involutions for the currently defined
  cartans.

  NOTE: this is not inlined to avoid a dependency on cartan.h.
*/

{
  return d_cartan->numInvolutions();
}

size_t ComplexReductiveGroup::numRealForms() const

/*!
  Synopsis: returns the number of weak real forms for this inner class.

  NOTE: this is not inlined to avoid a dependency on cartan.h.
*/

{
  return d_cartan->numRealForms();
}

realform::RealForm ComplexReductiveGroup::quasisplit() const

/*!
  Synopsis: returns the quasisplit real form.

  NOTE: this is not inlined to avoid a dependency on cartan.h.
*/

{
  return d_cartan->quasisplit();
}

size_t ComplexReductiveGroup::rank() const

/*!
  Synopsis: returns the rank of the group.

  NOTE: this is not inlined to avoid a dependency on rootdata.h.
*/

{
  return d_rootDatum->rank();
}

const realform::RealFormList& ComplexReductiveGroup::realFormLabels(size_t cn) 
  const

/*!
  Synopsis: returns the real form labels for cartan \#cn

  NOTE: this is not inlined to avoid a dependency on cartan.h.
*/

{
  return d_cartan->realFormLabels(cn);
}

unsigned long ComplexReductiveGroup::representative(realform::RealForm rf, 
						    size_t cn)
  const

/*!
  Synopsis: returns an element of the orbit corresponding to rf in the
  classification of weak real forms for cartan \#cn.

  Precondition: cartan \#cn is defined for rf.

  NOTE: this is not inlined to avoid a dependency on cartan.h.
*/

{
  return d_cartan->representative(rf,cn);
}

size_t ComplexReductiveGroup::semisimpleRank() const

/*!
  Synopsis: returns the semisimple rank of the group.

  NOTE: this is not inlined to avoid a dependency on rootdata.h.
*/

{
  return d_rootDatum->semisimpleRank();
}

const weyl::WeylGroup& ComplexReductiveGroup::weylGroup() const

/*!
  Synopsis: returns a reference to the Weyl group, which is now owned by the
  Tits group.

  NOTE: this is not inlined to avoid a dependency on tits.h
*/

{
  return d_titsGroup->weylGroup();
}

const weyl::WeylElt& ComplexReductiveGroup::twistedInvolution(size_t cn) const

/*!
  Synopsis: returns the twisted involution representative for class \#cn.

  NOTE : this is not inlined to avoid a compiling dependency on cartan.h
*/

{
  return d_cartan->twistedInvolution(cn);
}

/******** manipulators *******************************************************/

void ComplexReductiveGroup::fillCartan(realform::RealForm rf)

/*!
  Synopsis: fills in the Cartan classes that are defined for the real form x.

  NOTE: may forward an Overflow exception.
*/

{
  d_cartan->extend(weylGroup(),rootDatum(),rf);
  return;
}

void ComplexReductiveGroup::swap(ComplexReductiveGroup& other)

/*!
  NOTE: the components being given as pointers, it is enough to swap these!
*/

{
  std::swap(d_cartan,other.d_cartan);
  std::swap(d_rootDatum,other.d_rootDatum);
  std::swap(d_titsGroup,other.d_titsGroup);

  return;
}

}

/*****************************************************************************

        Chapter III -- Functions declared in complexredgp.h

  ... explain here when it is stable ...

******************************************************************************/

namespace complexredgp {

void lieType(lietype::LieType& lt, const ComplexReductiveGroup& G)

/*!
  Synopsis: puts in lt the Lie type of G.
*/

{  
  using namespace latticetypes;
  using namespace lietype;
  using namespace rootdata;

  lt.clear();

  const RootDatum& rd = G.rootDatum();
  LatticeMatrix cm;

  cartanMatrix(cm,rd);
  dynkin::lieType(lt,cm);

  // add the torus factor

  if (!rd.isSemisimple()) {
    SimpleLieType slt('T',rd.rank()-rd.semisimpleRank());
    lt.push_back(slt);
  }

  return;
}

}

}
