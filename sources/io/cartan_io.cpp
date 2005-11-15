/*
  This is cartan_io.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
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
#include "setutils.h"
#include "tori.h"

/*****************************************************************************

  ... explain here when it is stable ...

******************************************************************************/

namespace atlas {

namespace {

  typedef abelian::FiniteAbelianGroup AbGrp;
  template<typename T> void ignore(const T&) {}

}

/*****************************************************************************

        Chapter I -- Functions declared in cartan_io.h

  ... explain here when it is stable ...

******************************************************************************/

namespace cartan_io {

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

  for (PartitionIterator i(pi); i(); ++i) {
    std::ostringstream os;
    os << "real form #";
    os << rfl[c] << ": ";
    std::vector<unsigned long>::const_iterator first = i->first;
    std::vector<unsigned long>::const_iterator last = i->second;
    seqPrint(os,first,last,",","[","]");
    os << " (" << last - first << ")" << std::endl;
    std::string line = os.str();
    foldLine(strm,line,"",",");
    ++c;
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
  Permutation a;
  DynkinDiagram d(cm);
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
      Grading gr;
      f.grading(gr,*j);
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
