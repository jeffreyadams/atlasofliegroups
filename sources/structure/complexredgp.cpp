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

******************************************************************************/

namespace complexredgp {

ComplexReductiveGroup::ComplexReductiveGroup()
  :d_rootDatum(0),
   d_titsGroup(0),
   d_cartan(0)

/*!
  The default constructor. It is called in order to obtain an object
  whose fields will then by replaced by actual pointers to objects. But before
  that is done, we must make sure that there are null pointers, so that if
  the construction is abandoned, proper desctruction will occur.
*/

{}

ComplexReductiveGroup::
ComplexReductiveGroup(const rootdata::RootDatum* rd,
		      const latticetypes::LatticeMatrix& d)
  :d_rootDatum(rd),
   d_titsGroup(0),
   d_cartan(0)

/*!
  The main constructor constructs a ComplexReductiveGroup from its root datum
  and a distinguished involution.

  Precondition: d is an involution of the weight lattice of rd, which fixes rd
  as a _based_ root datum (it permutes the simple roots). In the complex
  group, this corresponds to a an involution fixing a chosen pinning, and in
  particular the Borel subgroup.

  NOTE: the ComplexReductiveGroup assumes ownership of the RootDatum pointed
  to by rd. Perhaps this is not entirely safe.
*/

{
  d_titsGroup = new tits::TitsGroup(rootDatum(),d);
  d_cartan = new cartan::CartanClasses(rootDatum(),d,d_titsGroup->weylGroup());
  /* final call stores d in d_cartan->d_fundamental.d_torus->d_involution,
     and constructs the fundamental fibers for the group and dual group */
}

/*!
  \brief constructs the complex reductive group dual to G.
*/
ComplexReductiveGroup::ComplexReductiveGroup(const ComplexReductiveGroup& G,
					     tags::DualTag)
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
  d_cartan = new CartanClasses(rootDatum(),d,d_titsGroup->weylGroup());
}

ComplexReductiveGroup::~ComplexReductiveGroup()

{
  delete d_cartan;
  delete d_titsGroup;
  delete d_rootDatum;
}

/******** accessors **********************************************************/

/* NOTE: most accessors are just forwarding functions to one of the three
   components of this class.
   They are not inlined to limit the dependencies of te header file.

   N.B. Given the immense number of methods that are simply forwarded to
   |d_cartan|, one may wonder it it would no have been a better idea to derive
   this class from CartanClassSet. To that one can oppose on one hand that an
   inner class "has" rather than "is" a set of Cartan classes (if it is
   anything, it is a class of real forms rather than Cartan subgroups, but
   even that is not the point of view taken in this software library), and on
   the other hand that this would force us to construct the equivalent of
   |d_cartan| _before_ the other data members, which would cause difficulties.
*/

/*!
  \brief the size of the block defined by the weak real form rf
  and the weak dual real form drf.
*/
unsigned long ComplexReductiveGroup::blockSize(realform::RealForm rf,
					       realform::RealForm drf) const
{
  return cartan::blockSize(rf,drf,*d_cartan);
}


/*!
  \brief returns Cartan subgroup number \#cn in the group.
*/
const cartanclass::CartanClass& ComplexReductiveGroup::cartan(size_t cn) const
{
  return d_cartan->cartan(cn);
}

/*!
  \brief returns the ordering of the Cartan subgroups
*/
const poset::Poset& ComplexReductiveGroup::cartanOrdering() const
{
  return d_cartan->ordering();
}

/*!
  \brief returns the support of the set of Cartan classes for rf
*/
const bitmap::BitMap& ComplexReductiveGroup::cartanSet(realform::RealForm rf)
  const
{
  return d_cartan->support(rf);
}

/*!
  \brief returns the support of the set of Cartan classes for the dual real
  form rf
*/
const bitmap::BitMap& ComplexReductiveGroup::dualCartanSet
  (realform::RealForm rf)
  const
{
  return d_cartan->dualSupport(rf);
}

/*!
  \brief returns the matrix of the distinguished involution.
*/
const latticetypes::LatticeMatrix& ComplexReductiveGroup::distinguished() const
{
  return d_cartan->distinguished();
}

/*!
  \brief returns the size of the fiber size corresponding to dual real
  form \#drf and cartan \#cn.

  Explanation: this is the size of the orbits, for the shifted action of
  the dual imaginary Weyl group (a.k.a. the real Weyl group), that correspond
  to drf in the classification of strong real forms (they all have the same
  size.)

  This is a technical function used for size computations.
*/
unsigned long ComplexReductiveGroup::dualFiberSize
  (realform::RealForm drf, size_t cn) const
{
  return d_cartan->dualFiberSize(drf,cn);
}

/*!
  \brief returns the dual fundamental fiber.

  This is a technical data structure containing the data for the classification
  of weak and strong dual real forms of G.
*/
const cartanclass::Fiber& ComplexReductiveGroup::dualFundamental() const
{
  return d_cartan->dualFundamental();
}

/*!
  \brief returns the dual real form labels for cartan \#cn
*/
const realform::RealFormList&
ComplexReductiveGroup::dualRealFormLabels(size_t cn) const
{
  return d_cartan->dualRealFormLabels(cn);
}

/*!
  \brief returns an element of the orbit corresponding to drf in the
  classification of weak real forms for cartan \#cn.

  Precondition: cartan \#cn is defined for drf.
*/
unsigned long ComplexReductiveGroup::dualRepresentative
 (realform::RealForm drf, size_t cn) const
{
  return d_cartan->dualRepresentative(drf,cn);
}

/*!
  \brief returns the size of the fiber size corresponding to real
  form \#rf and cartan \#cn.

  Explanation: this is the size of the orbits, for the shifted action of
  the imaginary Weyl group, that correspond to rf in the classification of
  strong real forms (they all have the same size.)

  This is a technical function used for size computations.
*/
unsigned long ComplexReductiveGroup::fiberSize
  (realform::RealForm rf, size_t cn) const
{
  return d_cartan->fiberSize(rf,cn);
}

/*!
  \brief returns the fundamental fiber.

  This is a technical data structure containing the data for the classification
  of weak and strong real forms of G.
*/
const cartanclass::Fiber& ComplexReductiveGroup::fundamental() const
{
  return d_cartan->fundamental();
}

/*!
  \brief puts in rs the set of noncompact imaginary roots for the
  representative of rf
*/
void ComplexReductiveGroup::grading
  (rootdata::RootSet& rs, realform::RealForm rf) const
{
  d_cartan->noncompactRootSet(rs,rf);
}

/*!
  \brief returns the number of elements in K\\G/B for real form \#rf.

  Precondition: the Cartan classes for this real form have been generated

  Explanation: this is exactly the number of elements in the one-sided
  parameter set corresponding to any strong real form of G lying over rf.
*/
unsigned long ComplexReductiveGroup::kgbSize(realform::RealForm rf) const
{
  return cartan::kgbSize(rf,*d_cartan);
}

/*!
  \brief returns the number of conjugacy classes of Cartan subgroups
  currently constructed for G. Only after fillCartan() has been called is it
  ensured that this gives the total number of Cartan subgroups for this inner
  class.
*/
size_t ComplexReductiveGroup::numCartanClasses() const
{
  return d_cartan->numCartan();
}

size_t ComplexReductiveGroup::mostSplit(realform::RealForm rf) const

/*!
  \brief returns the most split cartan subgroup for real form \#rf.

  Precondition: fillCartan() has been called for rf.
*/

{
  return d_cartan->mostSplit(rf);
}


/*!
  \brief returns the number of weak dual real forms for this inner class.
*/
size_t ComplexReductiveGroup::numDualRealForms() const
{
  return d_cartan->numDualRealForms();
}


/*!
  \brief returns the number of involutions for the currently defined
  cartans.
*/
size_t ComplexReductiveGroup::numInvolutions() const
{
  return d_cartan->numInvolutions();
}


/*!
  \brief returns the number of weak real forms for this inner class.
*/
size_t ComplexReductiveGroup::numRealForms() const
{
  return d_cartan->numRealForms();
}


/*!
  \brief returns the quasisplit real form.
*/
realform::RealForm ComplexReductiveGroup::quasisplit() const
{
  return d_cartan->quasisplit();
}

/*!
  \brief returns the rank of the group.
*/
size_t ComplexReductiveGroup::rank() const
{
  return d_rootDatum->rank();
}


/*!
  \brief returns the real form labels for cartan \#cn

  More precisely, realFormLabels(cn)[i] is the (inner) number of the real form
  that corresponds to part i of the partition cartan(n).fiber().weakReal()
*/
const realform::RealFormList& ComplexReductiveGroup::realFormLabels(size_t cn)
  const
{
  return d_cartan->realFormLabels(cn);
}

/*!
  \brief returns an element of the orbit corresponding to rf in the
  classification of weak real forms for cartan \#cn.

  Precondition: cartan \#cn is defined for rf.
*/
unsigned long ComplexReductiveGroup::representative
  (realform::RealForm rf,size_t cn) const
{
  return d_cartan->representative(rf,cn);
}

/*!
  \brief returns the semisimple rank of the group.
*/
size_t ComplexReductiveGroup::semisimpleRank() const
{
  return d_rootDatum->semisimpleRank();
}

/*!
  \brief returns a reference to the Weyl group, which is now owned by the
  Tits group.
*/
const weyl::WeylGroup& ComplexReductiveGroup::weylGroup() const
{
  return d_titsGroup->weylGroup();
}

/*!
  \brief returns the twisted involution representative for class \#cn.
*/
const weyl::TwistedInvolution&
  ComplexReductiveGroup::twistedInvolution(size_t cn) const
{
  return d_cartan->twistedInvolution(cn);
}

/******** manipulators *******************************************************/

/*!
  \brief fills in the Cartan classes that are defined for the real form x.
*/
void ComplexReductiveGroup::fillCartan(realform::RealForm rf)
{
  d_cartan->extend(d_titsGroup->weylGroup(),rootDatum(),rf);
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

} // namespace complexredgp

/*****************************************************************************

        Chapter III -- Functions declared in complexredgp.h

******************************************************************************/

namespace complexredgp {

/*!
  \brief puts in lt the Lie type of G.
*/
void lieType(lietype::LieType& lt, const ComplexReductiveGroup& G)
{
  lt.clear();

  const rootdata::RootDatum& rd = G.rootDatum();
  latticetypes::LatticeMatrix cm;

  rootdata::cartanMatrix(cm,rd);
  dynkin::lieType(lt,cm);

  // add the torus factor

  if (!rd.isSemisimple()) {
    lietype::SimpleLieType slt('T',rd.rank()-rd.semisimpleRank());
    lt.push_back(slt);
  }
}

}

}
