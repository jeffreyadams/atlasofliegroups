/*
  This is realredgp_io.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#include <iostream>
#include <sstream>

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

  ... explain here when it is stable ...

******************************************************************************/

namespace realredgp_io {

void Interface::swap(Interface& other)

{
  std::swap(d_realGroup,other.d_realGroup);
  std::swap(d_complexInterface,other.d_complexInterface);

  return;
}

/******** accessors **********************************************************/

const complexredgp::ComplexReductiveGroup& Interface::complexGroup() const

/*
  Synopsis: returns the underyling complex reductive group.

  NOTE: this is not inlined to avoid a dependency on realredgp.h.
*/

{
  return d_realGroup->complexGroup();
}

const realform_io::Interface& Interface::realFormInterface() const

/*
  Synopsis: returns the real form interface of the underlying complex group
  interface.

  NOTE: this is not inlined to avoid a dependency on complexredgp_io.h.
*/

{
  return d_complexInterface->realFormInterface();
}

}

/*****************************************************************************

        Chapter II -- Functions declared in realredgp_io.h

  ... explain here when it is stable ...

******************************************************************************/

namespace realredgp_io {

std::ostream& printBlockStabilizer(std::ostream& strm, const Interface& RI, 
				   size_t cn, realform::RealForm drf)

{
  using namespace cartanclass;
  using namespace complexredgp;
  using namespace realform;
  using namespace realredgp;
  using namespace realweyl;
  using namespace rootdata;
  using namespace size;
  using namespace weyl;

  const RealReductiveGroup& G_R = RI.realGroup();
  const ComplexReductiveGroup& G_C = G_R.complexGroup();
  const RootDatum& rd = G_R.rootDatum();
  const WeylGroup& W = G_R.weylGroup();

  RealForm rf = G_R.realForm();

  unsigned long x = G_C.representative(rf,cn);
  unsigned long y = G_C.dualRepresentative(drf,cn);
  const CartanClass& cc = G_C.cartan(cn);

  RealWeyl rw(cc,x,y,rd,W);
  RealWeylGenerators rwg(rw,cc,rd);

  realweyl_io::printBlockStabilizer(strm,rw,rwg,rd);

  // check if the size is correct
  Size c;
  blockStabilizerSize(c,rw);
  c *= G_C.fiberSize(rf,cn);
  c *= G_C.dualFiberSize(drf,cn);
  c *= cc.orbitSize();
  assert(c == W.order());

  return strm;
}

std::ostream& printCartanClasses(std::ostream& strm, 
				 const realredgp_io::Interface& G_RI)

/*
  Synopsis: outputs information about all the Cartan classes for
  G_RI.realGroup().
*/

{
  using namespace bitmap;
  using namespace complexredgp_io;
  using namespace ioutils;
  using namespace realredgp;

  const RealReductiveGroup& G_R = G_RI.realGroup();
  const complexredgp_io::Interface& G_CI = G_RI.complexInterface();

  const BitMap& b = G_R.cartanSet();
  bool first = true;

  for (BitMap::iterator i = b.begin(); i != b.end(); ++i) {
    if (first)
      first = false;
    else
      strm << std::endl << std::endl;      
    strm << "Cartan #" << *i << ":" << std::endl;
    printCartanClass(strm,*i,G_CI);
  }

  return strm;
}

std::ostream& printCartanOrder(std::ostream& strm, 
			       const realredgp::RealReductiveGroup& G_R)

/*
  Synopsis: prints the Hasse diagram of the Cartan ordering of G_R.
*/

{
  using namespace basic_io;
  using namespace graph;
  using namespace poset;

  const Poset& p = G_R.cartanOrdering();
  size_t cn = G_R.mostSplit();

  OrientedGraph g(0);
  p.hasseDiagram(g,cn);

  strm << "0:" << std::endl;

  for (size_t j = 1; j < g.size(); ++j) {
    const VertexList& e = g.edges(j);
    if (e.empty())
      continue;
    strm << j << ": ";
    VertexList::const_iterator first = e.begin();
    VertexList::const_iterator last = e.end();
    seqPrint(strm,first,last) << std::endl;
  }

  return strm;
}

std::ostream& printRealWeyl(std::ostream& strm, const Interface& RI, size_t cn)

/*
  Synopsis: outputs the real Weyl group corresponding to Cartan #cn.

  Precondition: cartan #cn is defined for this real form.
*/

{
  using namespace cartanclass;
  using namespace complexredgp;
  using namespace realredgp;
  using namespace realweyl;
  using namespace realweyl_io;
  using namespace rootdata;
  using namespace size;
  using namespace weyl;

  const RealReductiveGroup& G_R = RI.realGroup();
  const ComplexReductiveGroup& G_C = G_R.complexGroup();

  unsigned long rf = G_R.realForm();

  const RootDatum& rd = G_C.rootDatum();
  const WeylGroup& W = G_C.weylGroup();
  const CartanClass& cc = G_C.cartan(cn);
  unsigned long x = G_C.representative(rf,cn);

  RealWeyl rw(cc,x,0,rd,W);
  RealWeylGenerators rwg(rw,cc,rd);

  realweyl_io::printRealWeyl(strm,rw,rwg,rd);

  // check if the size is correct
  Size c;
  realWeylSize(c,rw);
  c *= G_C.fiberSize(rf,cn);
  c *= cc.orbitSize();
  assert(c == W.order());

  return strm;
}

std::ostream& printStrongReal(std::ostream& strm, const Interface& RI, 
			      size_t cn)

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

{
  using namespace basic_io;
  using namespace cartanclass;
  using namespace complexredgp;
  using namespace ioutils;
  using namespace partition;
  using namespace realform;
  using namespace rootdata;

  const ComplexReductiveGroup& G_C = RI.complexGroup();
  const CartanClass& cc = G_C.cartan(cn);

  size_t n = cc.numRealFormClasses();

  if (n == 1) { // simplified output

    const Partition& pi = cc.strongReal(0);
    const realform_io::Interface& rfi = RI.realFormInterface();
    const RealFormList& rfl = G_C.realFormLabels(cn);

    unsigned long c = 0;

    for (PartitionIterator i(pi); i(); ++i) {
      std::ostringstream os;
      os << "real form #";
      RealForm rf = rfl[cc.toWeakReal(c,0)];
      os << rfi.out(rf) << ": ";
      std::vector<unsigned long>::const_iterator first = i->first;
      std::vector<unsigned long>::const_iterator last = i->second;
      seqPrint(os,first,last,",","[","]");
      os << " (" << last - first << ")" << std::endl;
      std::string line = os.str();
      foldLine(strm,line,"",",");
      ++c;
    }

    // print information about the center

    return strm;
  }

  strm << "there are " << n << " real form classes:" << std::endl << std::endl;

  for (size_t j = 0; j < n; ++j) {
    strm << "class #" << j << ":" << std::endl;

    const Partition& pi = cc.strongReal(j);
    const realform_io::Interface& rfi = RI.realFormInterface();
    const RealFormList& rfl = G_C.realFormLabels(cn);

    unsigned long c = 0;

    for (PartitionIterator i(pi); i(); ++i) {
      std::ostringstream os;
      os << "real form #";
      RealForm rf = rfl[cc.toWeakReal(c,j)];
      os << rfi.out(rf) << ": ";
      std::vector<unsigned long>::const_iterator first = i->first;
      std::vector<unsigned long>::const_iterator last = i->second;
      seqPrint(os,first,last,",","[","]");
      os << " (" << last - first << ")" << std::endl;
      std::string line = os.str();
      foldLine(strm,line,"",",");
      ++c;
    }

    // print information about the center

    strm << std::endl;
  }

  return strm;
}

}

}
