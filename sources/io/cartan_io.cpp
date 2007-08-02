/*
  This is cartan_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "cartan_io.h"

#include <sstream>

#include "basic_io.h"
#include "bitset.h"
#include "cartanset.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "complexredgp_io.h"
#include "dynkin.h"
#include "gradings.h"
#include "ioutils.h"
#include "lietype.h"
#include "prettyprint.h"
#include "rootdata.h"
#include "setutils.h"
#include "tori.h"


namespace atlas {

namespace {

  typedef abelian::FiniteAbelianGroup AbGrp;
  template<typename T> void ignore(const T&) {}

}

/*****************************************************************************

        Chapter I -- Functions declared in cartan_io.h

******************************************************************************/

namespace cartan_io {

std::ostream& printCartanClass(std::ostream& strm, size_t cn,
			       const complexredgp_io::Interface& CI)

/*
  Synopsis: prints information about the Cartan class #cn.
*/

{
  using namespace basic_io; // for operators <<

  const complexredgp::ComplexReductiveGroup& G = CI.complexGroup();
  const rootdata::RootDatum& rd = G.rootDatum();

  const cartanclass::CartanClass& cc = G.cartan(cn);

  strm << "canonical twisted involution: ";
  prettyprint::printWeylElt
    (strm,G.cartanClasses().twistedInvolution(cn),G.weylGroup())
       << std::endl;

  prettyprint::printTorusType(strm,cc.fiber().torus()) << std::endl;
  strm << "twisted involution orbit size: " << cc.orbitSize() << std::endl;

  // print type of imaginary root system
  lietype::LieType ilt;
  rootdata::lieType(ilt,cc.simpleImaginary(),rd);

  if (ilt.size() == 0)
    strm << "imaginary root system is empty" << std::endl;
  else
    strm << "imaginary root system: " << ilt << std::endl;

  // print type of real root system
  lietype::LieType rlt;
  rootdata::lieType(rlt,cc.simpleReal(),rd);

  if (rlt.size() == 0)
    strm << "real root system is empty" << std::endl;
  else
    strm << "real root system: " << rlt << std::endl;

  // print type of complex root system
  lietype::LieType clt;
  rootdata::lieType(clt,cc.simpleComplex(),rd);

  if (clt.size() == 0)
    strm << "complex factor is empty" << std::endl;
  else
    strm << "complex factor: " << clt << std::endl;

  realform::RealFormList rfl(cc.numRealForms());
  const realform_io::Interface& rfi = CI.realFormInterface();

  for (size_t i = 0; i < rfl.size(); ++i)
    rfl[i] = rfi.out(G.realFormLabels(cn)[i]);

  printFiber(strm,cc.fiber(),rfl);

  return strm;
}

std::ostream& printFiber(std::ostream& strm, const cartanclass::Fiber& f,
			 const realform::RealFormList& rfl)

/*
  Synopsis: prints out the fiber data.
*/

{
  using namespace basic_io;
  using namespace ioutils;
  using namespace partition;

  const Partition& pi = f.weakReal();
  unsigned long c = 0;

  for (PartitionIterator i(pi); i(); ++i,++c) {
    std::ostringstream os;
    os << "real form #";
    os << rfl[c] << ": ";
    std::vector<unsigned long>::const_iterator first = i->first;
    std::vector<unsigned long>::const_iterator last = i->second;
    seqPrint(os,first,last,",","[","]");
    os << " (" << last - first << ")" << std::endl;
    std::string line = os.str();
    foldLine(strm,line,"",",");
  }

  return strm;
}

std::ostream& printGradings(std::ostream& strm, const cartanclass::Fiber& f,
			    const realform::RealFormList& rfl,
			    const rootdata::RootDatum& rd)

/*
  Synopsis: outputs the gradings of the simple imaginary roots corresponding
  to the various real forms defined for f.

  Precondition: rfl contains the outer numbering of the real forms;

  The gradings are output in the same order as the orbit corresponding to the
  real form is output in the "cartan" command.
*/

{
  using namespace dynkin;
  using namespace gradings;
  using namespace ioutils;
  using namespace latticetypes;
  using namespace partition;
  using namespace prettyprint;
  using namespace setutils;
  using namespace rootdata;

  typedef std::vector<unsigned long>::const_iterator VI;

  const RootList& si = f.simpleImaginary();
  LatticeMatrix cm;

  cartanMatrix(cm,si,rd);
  DynkinDiagram d(cm);
  Permutation a;
  bourbaki(a,d);
  cm.permute(a);

  strm << "cartan matrix of imaginary root system is:" << std::endl;
  printMatrix(strm,cm);

  const Partition& pi = f.weakReal();
  unsigned long c = 0;

  for (PartitionIterator i(pi); i(); ++i) {

    std::ostringstream os;
    os << "real form #";
    os << rfl[c] << ": ";

    VI i_first = i->first;
    VI i_last = i->second;

    os << "[";

    for (VI j = i->first; j != i_last; ++j) {
      if (j != i_first)
	os << ",";
      Grading gr= f.grading(*j);
      gr.permute(a);
      prettyPrint(os,gr,f.simpleImaginary().size());
    }

    os << "]" << std::endl;

    std::string line = os.str();
    foldLine(strm,line,"",",");
    ++c;
  }

  return strm;
}

}

}
