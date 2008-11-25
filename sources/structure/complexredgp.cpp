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

  For license information see the LICENSE file
*/

#include "complexredgp.h"

#include "tags.h"
#include "cartanset.h"
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


/*****************************************************************************

        Chapter I -- The ComplexReductiveGroup class

******************************************************************************/

namespace complexredgp {

/*!
  \brief main constructor

  constructs a ComplexReductiveGroup from its root datum and a distinguished
  involution.

  Precondition: d is a based root datum involution: it globally fixes the set
  of positive roots

  NOTE: the ComplexReductiveGroup assumes ownership of the RootDatum pointed
  to by rd; users beware of this.
*/
ComplexReductiveGroup::ComplexReductiveGroup
 (const rootdata::RootDatum& rd, const latticetypes::LatticeMatrix& d)
  : d_rootDatum(rd) // we assume ownership
  , d_titsGroup(rootDatum(),d)
  , root_twist(make_root_twist())
  , d_cartanSet(*this,d)
  /* stores |d| in |d_cartanSet.d_fundamental.d_torus->d_involution|,
     and constructs the fundamental fibers for the group and dual group */
{}

/*!
  \brief constructs the complex reductive group dual to G.

  We have to compute the dual based involution twice, since we have no
  variables available during initialisation to store the previous value (and
  d_titsGroup does not save it)
*/
ComplexReductiveGroup::ComplexReductiveGroup(const ComplexReductiveGroup& G,
					     tags::DualTag)
  : d_rootDatum(G.rootDatum(),tags::DualTag())
  , d_titsGroup(rootDatum(),
		dualBasedInvolution(G.distinguished(),G.rootDatum()))
  , root_twist(make_root_twist())
  , d_cartanSet(*this,dualBasedInvolution(G.distinguished(),G.rootDatum()))
{}

setutils::Permutation ComplexReductiveGroup::make_root_twist() const
{
  setutils::Permutation twist(d_rootDatum.semisimpleRank());
  for (size_t i=0; i<twist.size(); ++i)
    twist[i]= d_titsGroup.twisted(i);

  return rootDatum().root_permutation(twist);
}

/******** accessors **********************************************************/

unsigned long ComplexReductiveGroup::blockSize(realform::RealForm rf,
					       realform::RealForm drf) const
{ return cartanset::blockSize(rf,drf,d_cartanSet); }


unsigned long ComplexReductiveGroup::kgbSize(realform::RealForm rf) const
{ return cartanset::kgbSize(rf,d_cartanSet); }


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

} // |namespace complexredgp|

} // |namespace atlas|
