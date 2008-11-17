/*
  This is realredgp_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include <iostream>
#include <sstream>
#include <cassert>

#include "realredgp_io.h"

#include "basic_io.h"
#include "bitmap.h"
#include "cartan_io.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "complexredgp_io.h"
#include "ioutils.h"
#include "partition.h"
#include "poset.h"
#include "realform.h"
#include "realredgp.h"
#include "realweyl.h"
#include "realweyl_io.h"
#include "rootdata.h"
#include "size.h"

namespace atlas {

/*****************************************************************************

        Chapter I -- The Interface class

******************************************************************************/

namespace realredgp_io {

void Interface::swap(Interface& other)

{
  std::swap(d_realGroup,other.d_realGroup);
  std::swap(d_complexInterface,other.d_complexInterface);
}

/******** accessors **********************************************************/

const complexredgp::ComplexReductiveGroup& Interface::complexGroup() const
{
  return d_realGroup->complexGroup();
}

const realform_io::Interface& Interface::realFormInterface() const
{
  return d_complexInterface->realFormInterface();
}

}

/*****************************************************************************

        Chapter II -- Functions declared in realredgp_io.h

******************************************************************************/

namespace realredgp_io {

std::ostream& printBlockStabilizer(std::ostream& strm,
				   const realredgp::RealReductiveGroup& G_R,
				   size_t cn, realform::RealForm drf)

{
  const complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
  const rootdata::RootDatum& rd = G_R.rootDatum();
  const weyl::WeylGroup& W = G_R.weylGroup();

  realform::RealForm rf = G_R.realForm();

  unsigned long x = G_C.representative(rf,cn);
  unsigned long y = G_C.dualRepresentative(drf,cn);
  const cartanclass::CartanClass& cc = G_C.cartan(cn);

  realweyl::RealWeyl rw(cc,x,y,rd,W);
  realweyl::RealWeylGenerators rwg(rw,cc,rd);

  realweyl_io::printBlockStabilizer(strm,rw,rwg);

  // check if the size is correct
  size::Size c;
  realweyl::blockStabilizerSize(c,rw);
  c *= G_C.fiberSize(rf,cn);
  c *= G_C.dualFiberSize(drf,cn);
  c *= cc.orbitSize();
  assert(c == W.order());

  return strm;
}

// Print information about all the Cartan classes for |G_RI.realGroup()|.
std::ostream& printCartanClasses(std::ostream& strm,
				 const realredgp_io::Interface& G_RI)
{
  const realredgp::RealReductiveGroup& G_R = G_RI.realGroup();
  const complexredgp_io::Interface& G_CI = G_RI.complexInterface();

  const bitmap::BitMap& b = G_R.cartanSet();
  bool first = true;

  for (bitmap::BitMap::iterator i = b.begin(); i(); ++i) {
    if (first)
      first = false;
    else
      strm << std::endl << std::endl;
    strm << "Cartan #" << *i << ":" << std::endl;
    cartan_io::printCartanClass(strm,*i,G_CI);
  }

  return strm;
}

/* the function below is much like |kgb_io::printBruhatOrder|, but cannot
   call that function, which needs a different type, nameley |BruhatOrder|.
 */

// Print the Hasse diagram of the Cartan ordering of G_R.
std::ostream& printCartanOrder(std::ostream& strm,
			       const realredgp::RealReductiveGroup& G_R)
{
  const poset::Poset& p = G_R.complexGroup().cartanOrdering();
  size_t cn = G_R.mostSplit();

  graph::OrientedGraph g(0);
  p.hasseDiagram(g,cn); // everything below |cn| is for real form of |G_R|

  strm << "0:" << std::endl; // this is the only minimal element

  // all covering relations are grouped by lower (covered) element
  for (size_t j = 1; j < g.size(); ++j)
  {
    const graph::VertexList& e = g.edgeList(j);
    if (not e.empty()) // suppress non-covered elements; they are not for |G_R|
      basic_io::seqPrint(strm << j << ": ",e.begin(),e.end()) << std::endl;
  }

  return strm;
}


/*
  Print the real Weyl group corresponding to Cartan #cn.

  Precondition: cartan #cn is defined for this real form.
*/
std::ostream& printRealWeyl(std::ostream& strm,
			    const realredgp::RealReductiveGroup& G_R,
			    size_t cn)
{
  const complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();

  realform::RealForm rf = G_R.realForm();

  const rootdata::RootDatum& rd = G_C.rootDatum();
  const weyl::WeylGroup& W = G_C.weylGroup();
  const cartanclass::CartanClass& cc = G_C.cartan(cn);
  cartanclass::AdjointFiberElt x = G_C.representative(rf,cn);

  realweyl::RealWeyl rw(cc,x,0,rd,W);
  realweyl::RealWeylGenerators rwg(rw,cc,rd);

  realweyl_io::printRealWeyl(strm,rw,rwg);

  // check if the size is correct
  size::Size c;
  realweyl::realWeylSize(c,rw);
  c *= G_C.fiberSize(rf,cn);
  c *= cc.orbitSize();
  assert(c == W.order());

  return strm;
}

/*
  Synopsis: outputs information about the strong real forms of G.

  Explanation: the inverse image in \X of a class of weak real forms of G is
  of the form Z.\X(z), where z is an admissible value for x^2 for a strong real
  form of that class; so there is associated to the class of weak real forms
  a coset (1+\delta)(Z).z in Z^\delta. All the various \X(z) for all possible
  choices of z are isomorphic by Z-translation. Therefore it is enough to
  describe the combinatorial structure of one of them, for each class of weak
  real forms. This is a finite problem in all cases.

  We output the orbits of W_im in X(z), which correspond to the various strong
  real forms; we label them with the corresponding weak real form.
*/
std::ostream& printStrongReal(std::ostream& strm,
			      const realredgp::RealReductiveGroup& G_R,
			      const realform_io::Interface& rfi,
			      size_t cn)
{
  const complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
  const cartanclass::CartanClass& cc = G_C.cartan(cn);

  size_t n = cc.numRealFormClasses();

  if (n>1)
    strm << "there are " << n << " real form classes:\n" << std::endl;

  for (size_t j = 0; j < n; ++j)
  {
    if (n>1)
      strm << "class #" << j << ":" << std::endl;

    const partition::Partition& pi = cc.strongReal(j);
    const realform::RealFormList& rfl = G_C.realFormLabels(cn);

    unsigned long c = 0;

    for (partition::PartitionIterator i(pi); i(); ++i,++c) {
      std::ostringstream os;
      realform::RealForm rf = rfl[cc.toWeakReal(c,j)];
      os << "real form #" << rfi.out(rf) << ": ";
      basic_io::seqPrint(os,i->first,i->second,",","[","]")
	<< " (" << i->second - i->first << ")" << std::endl;
      ioutils::foldLine(strm,os.str(),"",",");
    }

    // print information about the center (not implemented)

    if (n>1)
      strm << std::endl;
  }

  return strm;
}

}

}
