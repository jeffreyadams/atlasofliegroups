/*
  This is realweyl_io.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "realweyl_io.h"

#include <iostream>

#include "basic_io.h"
#include "prettyprint.h"
#include "realweyl.h"
#include "rootdata.h"

/*****************************************************************************

        Chapter I -- Functions declared in realweyl_io.h

******************************************************************************/

namespace atlas {

namespace realweyl_io {

std::ostream& printBlockStabilizer(std::ostream& strm,
				   const realweyl::RealWeyl& rw,
				   const realweyl::RealWeylGenerators& rwg,
				   const rootdata::RootDatum& rd)

/*
  Synopsis: prints the stabilizer of a representation in this W-orbit
  for the cross-action.

  Explanation: this is the intersection of the real and dual real Weyl groups.
*/

{
  using namespace basic_io;
  using namespace prettyprint;
  using namespace rootdata;

  strm << "real weyl group is W^C.((A_i.W_ic) x (A_r.W_rc)), where:" 
       << std::endl;

  if (rw.numComplex()) {
    strm << "W^C is isomorphic to a Weyl group of type " 
	 << rw.complexType() << std::endl;
  } else
    strm << "W^C is trivial" << std::endl;

  if (rw.numImaginaryR())
    strm << "A_i is an elementary abelian 2-group of rank "
	 << rw.numImaginaryR() << std::endl;
  else
    strm << "A_i is trivial" << std::endl;

  if (rw.numImaginaryCompact()) {
    strm << "W_ic is a Weyl group of type " 
	 << rw.imaginaryCompactType() << std::endl;
  } else
    strm << "W_ic is trivial" << std::endl;

  if (rw.numRealR())
    strm << "A_r is an elementary abelian 2-group of rank "
	 << rw.numRealR() << std::endl;
  else
    strm << "A_r is trivial" << std::endl;

  if (rw.numRealCompact()) {
    strm << "W_rc is a Weyl group of type " 
	 << rw.realCompactType() << std::endl;
  } else
    strm << "W_rc is trivial" << std::endl;

  strm << std::endl; 

  if (rw.numComplex()) {
    strm << "generators for W^C:" << std::endl;
    printWeylList(strm,rwg.complex(),rwg.weylGroup(),"\n","","\n");
  }

  if (rw.numImaginaryR()) {
    strm << "generators for A_i:" << std::endl;
    printWeylList(strm,rwg.imaginaryR(),rwg.weylGroup(),"\n","","\n");
  }

  if (rw.numImaginaryCompact()) {
    strm << "generators for W_ic:" << std::endl;
    printWeylList(strm,rwg.imaginaryCompact(),rwg.weylGroup(),"\n","","\n");
  }

  if (rw.numRealR()) {
    strm << "generators for A_r:" << std::endl;
    printWeylList(strm,rwg.realR(),rwg.weylGroup(),"\n","","\n");
  }

  if (rw.numRealCompact()) {
    strm << "generators for W_rc:" << std::endl;
    printWeylList(strm,rwg.realCompact(),rwg.weylGroup(),"\n","","\n");
  }

  return strm;
}

std::ostream& printDualRealWeyl(std::ostream& strm, 
				const realweyl::RealWeyl& rw,
				const realweyl::RealWeylGenerators& rwg,
				const rootdata::RootDatum& rd)

/*
  Synopsis: prints the dual real Weyl group.
*/

{
  using namespace basic_io;
  using namespace prettyprint;
  using namespace rootdata;

  strm << "real weyl group is W^C.(W^I x (A.W_rc)), where:" 
       << std::endl;

  if (rw.numComplex()) {
    strm << "W^C is isomorphic to a Weyl group of type " 
	 << rw.complexType() << std::endl;
  } else
    strm << "W^C is trivial" << std::endl;

  if (rw.numImaginary()) {
    strm << "W^I is a Weyl group of type " 
	 << rw.imaginaryType() << std::endl;
  } else
    strm << "W^I is trivial" << std::endl;

  if (rw.numRealR())
    strm << "A is an elementary abelian 2-group of rank "
	 << rw.numRealR() << std::endl;
  else
    strm << "A is trivial" << std::endl;

  if (rw.numRealCompact()) {
    strm << "W_rc is a Weyl group of type " 
	 << rw.realCompactType() << std::endl;
  } else
    strm << "W_rc is trivial" << std::endl;

  strm << std::endl; 

  if (rw.numComplex()) {
    strm << "generators for W^C:" << std::endl;
    printWeylList(strm,rwg.complex(),rwg.weylGroup(),"\n","","\n");
  }

  if (rw.numImaginary()) {
    strm << "generators for W^I:" << std::endl;
    printWeylList(strm,rwg.imaginary(),rwg.weylGroup(),"\n","","\n");
  }

  if (rw.numRealR()) {
    strm << "generators for A:" << std::endl;
    printWeylList(strm,rwg.realR(),rwg.weylGroup(),"\n","","\n");
  }

  if (rw.numRealCompact()) {
    strm << "generators for W_rc:" << std::endl;
    printWeylList(strm,rwg.realCompact(),rwg.weylGroup(),"\n","","\n");
  }

  return strm;
}

std::ostream& printRealWeyl(std::ostream& strm, const realweyl::RealWeyl& rw,
			    const realweyl::RealWeylGenerators& rwg,
			    const rootdata::RootDatum& rd)

/*
  Synopsis: prints the real Weyl group.
*/

{
  using namespace basic_io;
  using namespace prettyprint;
  using namespace rootdata;

  strm << "real weyl group is W^C.((A.W_ic) x W^R), where:" << std::endl;

  if (rw.numComplex()) {
    strm << "W^C is isomorphic to a Weyl group of type " 
	 << rw.complexType() << std::endl;
  } else
    strm << "W^C is trivial" << std::endl;

  if (rw.numImaginaryR())
    strm << "A is an elementary abelian 2-group of rank "
	 << rw.numImaginaryR() << std::endl;
  else
    strm << "A is trivial" << std::endl;

  if (rw.numImaginaryCompact()) {
    strm << "W_ic is a Weyl group of type " 
	 << rw.imaginaryCompactType() << std::endl;
  } else
    strm << "W_ic is trivial" << std::endl;

  if (rw.numReal()) {
    strm << "W^R is a Weyl group of type " 
	 << rw.realType() << std::endl;
  } else
    strm << "W^R is trivial" << std::endl;

  strm << std::endl; 

  if (rw.numComplex()) {
    strm << "generators for W^C:" << std::endl;
    printWeylList(strm,rwg.complex(),rwg.weylGroup(),"\n","","\n");
  }

  if (rw.numImaginaryR()) {
    strm << "generators for A:" << std::endl;
    printWeylList(strm,rwg.imaginaryR(),rwg.weylGroup(),"\n","","\n");
  }

  if (rw.numImaginaryCompact()) {
    strm << "generators for W_ic:" << std::endl;
    printWeylList(strm,rwg.imaginaryCompact(),rwg.weylGroup(),"\n","","\n");
  }

  if (rw.numReal()) {
    strm << "generators for W^R:" << std::endl;
    printWeylList(strm,rwg.real(),rwg.weylGroup(),"\n","","\n");
  }

  return strm;
}

}

}
