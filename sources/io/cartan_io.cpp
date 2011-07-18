/*
  This is cartan_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "cartan_io.h"

#include <sstream>

#include "basic_io.h"
#include "bitset.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "complexredgp_io.h"
#include "dynkin.h"
#include "gradings.h"
#include "ioutils.h"
#include "lietype.h"
#include "prettyprint.h"
#include "rootdata.h"
#include "permutations.h"
#include "tori.h"


namespace atlas {

/*****************************************************************************

        Chapter I -- Functions declared in cartan_io.h

******************************************************************************/

namespace cartan_io {


// Print information about the Cartan class #cn.
std::ostream& printCartanClass(std::ostream& strm, size_t cn,
			       complexredgp_io::Interface& CI)
{
  ComplexReductiveGroup& G = CI.complexGroup();
  const RootSystem& rs = G.rootDatum();

  const CartanClass& cc = G.cartan(cn);
  const Fiber& f = cc.fiber();

  prettyprint::printTorusType(strm,f.torus()) << std::endl;

  {
    std::ostringstream os;
    os << "canonical twisted involution: ";
    prettyprint::printWeylElt(os,G.twistedInvolution(cn),G.weylGroup());
    ioutils::foldLine(strm,os.str(),"",",") << std::endl;
  }

  size_t orbit_size=cc.orbitSize();
  strm << "twisted involution orbit size: " << orbit_size
       << "; fiber size: " << f.fiberSize()
       << "; strong inv: " << orbit_size*f.fiberSize()
       <<std::endl;

  // print type of imaginary root system
  LieType ilt = rs.Lie_type(cc.simpleImaginary());

  if (ilt.size() == 0)
    strm << "imaginary root system is empty" << std::endl;
  else
    strm << "imaginary root system: " << ilt << std::endl;

  // print type of real root system
  LieType rlt = rs.Lie_type(cc.simpleReal());

  if (rlt.size() == 0)
    strm << "real root system is empty" << std::endl;
  else
    strm << "real root system: " << rlt << std::endl;

  // print type of complex root system
  LieType clt = rs.Lie_type(cc.simpleComplex());

  if (clt.size() == 0)
    strm << "complex factor is empty" << std::endl;
  else
    strm << "complex factor: " << clt << std::endl;

  RealFormNbrList rfl(cc.numRealForms());
  const realform_io::Interface& rfi = CI.realFormInterface();

  for (size_t i = 0; i < rfl.size(); ++i)
    rfl[i] = rfi.out(G.realFormLabels(cn)[i]);

  printFiber(strm,f,rfl);

  return strm;
}


// Print the fiber data.
std::ostream& printFiber(std::ostream& strm, const Fiber& f,
			 const RealFormNbrList& rfl)
{
  const partition::Partition& pi = f.weakReal();
  unsigned long c = 0;

  for (partition::PartitionIterator i(pi); i(); ++i,++c)
  {
    std::ostringstream os;
    os << "real form #";
    os << rfl[c] << ": ";
    basic_io::seqPrint(os,i->first,i->second,",","[","]");
    os << " (" << i->second - i->first << ")" << std::endl;
    ioutils::foldLine(strm,os.str(),"",",");
  }

  return strm;
}

/*
  Print the gradings of the simple-imaginary roots corresponding
  to the various real forms defined for |f|.

  Precondition: |rfl| contains the outer numbering of the real forms;

  The gradings are output in the same order as the orbit corresponding to the
  real form is output in the "cartan" command.
*/
std::ostream& printGradings(std::ostream& strm, const Fiber& f,
			    const RealFormNbrList& rfl,
			    const RootSystem& rs)
{
  typedef std::vector<unsigned long>::const_iterator VI;

  const RootNbrList& si = f.simpleImaginary();
  int_Matrix cm = rs.cartanMatrix(si);
  dynkin::DynkinDiagram d(cm);
  permutations::Permutation a = dynkin::bourbaki(d);
  a.inv_conjugate(cm);

  strm << "cartan matrix of imaginary root system is:" << std::endl;
  prettyprint::printMatrix(strm,cm);

  const partition::Partition& pi = f.weakReal();
  unsigned long c = 0;

  for (partition::PartitionIterator i(pi); i(); ++i) {

    std::ostringstream os;
    os << "real form #";
    os << rfl[c] << ": ";

    VI i_first = i->first;
    VI i_last = i->second;

    os << '[';

    for (VI j = i->first; j != i_last; ++j) {
      if (j != i_first)
	os << ',';
      Grading gr= a.pull_back(f.grading(*j));
      prettyprint::prettyPrint(os,gr,f.simpleImaginary().size());
    }

    os << ']' << std::endl;

    std::string line = os.str();
    ioutils::foldLine(strm,line,"",",");
    ++c;
  }

  return strm;
}

} // namespace cartan_io

} // namespace atlas
