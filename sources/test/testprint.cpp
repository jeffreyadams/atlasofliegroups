/*
  This is testprint.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "testprint.h"

#include <iostream>
#include <iomanip>
#include <iterator>
#include <sstream>

#include "basic_io.h"
#include "cartan_io.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "complexredgp_io.h"
#include "dynkin.h"
#include "ioutils.h"
#include "lattice.h"
#include "latticetypes.h"
#include "lietype.h"
#include "prettyprint.h"
#include "realredgp.h"
#include "rootdata.h"
#include "size.h"
#include "tori.h"
#include "weylsize.h"

/*****************************************************************************

  This module contains a number of print commands vor viewing parts of the
  data which would normally be hidden away. It is intended for testing
  purposes.

******************************************************************************/

namespace atlas {

namespace {

struct OrbitData {
  unsigned long size;
  std::vector<unsigned long> fiberSizes;
  std::vector<unsigned long> dualFiberSizes;
};

}

/*****************************************************************************

        Chapter I -- Functions defined in testprint.h

  ... explain here when it is stable ...

******************************************************************************/

namespace testprint {


/*
  This function outputs the following data about the rootdatum :

    - Cartan matrix;
    - set of simple roots and coroots;
    - if not semisimple : bases for radical and coradical;
    - positive roots and coroots;
*/
std::ostream& print(std::ostream& strm, const rootdata::RootDatum& rd)
{
  strm << "cartan matrix :" << std::endl;
  latticetypes::LatticeMatrix q(rd.semisimpleRank());
  for (size_t j = 0; j < rd.semisimpleRank(); ++j)
    for (size_t i = 0; i < rd.semisimpleRank(); ++i)
      q(i,j) = rd.cartan(i,j);
  prettyprint::printMatrix(strm,q) << std::endl;

  strm << "root basis :" << std::endl;
  latticetypes::WeightList r_rb(rd.beginSimpleRoot(),rd.endSimpleRoot());
  prettyprint::printBasis(strm,r_rb) << std::endl;

  strm << "coroot basis :" << std::endl;
  latticetypes::WeightList r_crb(rd.beginSimpleCoroot(),rd.endSimpleCoroot());
  prettyprint::printBasis(strm,r_crb) << std::endl;

  if (not rd.isSemisimple()) { // print radical and coradical bases
    strm << "radical basis :" << std::endl;
    latticetypes::WeightList r_rad(rd.beginRadical(),rd.endRadical());
    prettyprint::printBasis(strm,r_rad) << std::endl;

    strm << "coradical basis :" << std::endl;
    latticetypes::WeightList r_crad(rd.beginCoradical(),rd.endCoradical());
    prettyprint::printBasis(strm,r_crad) << std::endl;
  }

  strm << "positive roots :" << std::endl;
  latticetypes::WeightList r_pr(rd.beginPosRoot(),rd.endPosRoot());
  basic_io::seqPrint(strm,r_pr.begin(),r_pr.end(),"\n","","")
    << std::endl << std::endl;

  strm << "positive coroots :" << std::endl;
  latticetypes::WeightList r_pcr(rd.beginPosCoroot(),rd.endPosCoroot());
  basic_io::seqPrint(strm,r_pcr.begin(),r_pcr.end(),"\n","","")
    << std::endl << std::endl;

  return strm;
}


/*
  Synopsis: prints out data about the real reductive group G.

  NOTE: to be implemented.

std::ostream& print
  (std::ostream& strm, const realredgp::RealReductiveGroup& G)
{
  return strm;
}
*/



/*
  Synopsis: produces the file as posted on the atlas website.

  Precondition: G.fullCartan() has been called successfully.
*/
std::ostream& printBlockData(std::ostream& strm,
			     const complexredgp_io::Interface& CI)
{
  const complexredgp::ComplexReductiveGroup& G = CI.complexGroup();
  const rootdata::RootDatum& rd = G.rootDatum();

  lietype::LieType lt;
  complexredgp::lieType(lt,G);

  strm << std::setw(25) << "" << "Block data"
       << std::endl;

  strm << std::setw(16) << "" << "for the";

  if (rd.isSimplyConnected())
    strm << " simply connected";
  else if (rd.isAdjoint())
    strm << " adjoint";

  strm << " group of type " << lt
       << std::endl;

  strm << std::endl;

  printCartanClasses(strm,CI);

  return strm;
}


/*
  Prints information about the conjugacy classes of Cartan subgroups.
*/
std::ostream& printCartanClasses(std::ostream& strm,
				 const complexredgp_io::Interface& CI)
{
  const complexredgp::ComplexReductiveGroup& G = CI.complexGroup();

  for (unsigned long j = 0; j < G.numCartanClasses(); ++j) {
    strm << "Cartan #" << j << ":" << std::endl;
    cartan_io::printCartanClass(strm,j,CI);
    if (j+1 < G.numCartanClasses())
      strm << std::endl << std::endl;
  }

  return strm;
}


/*
  In this function, rb should contain a basis of some rootsystem; we are
  outputting the matrix c where c(i,j) = root(rb[i]).coroot(rb[j]).
*/
std::ostream& printCartanMatrix(std::ostream& strm,
				const rootdata::RootList& rb,
				const rootdata::RootDatum& rd)
{
  latticetypes::LatticeMatrix q(rb.size());

  for (size_t j = 0; j < rb.size(); ++j)
    for (size_t i = 0; i < rb.size(); ++i)
      q(i,j) = rd.scalarProduct(rb[i],rb[j]);

  prettyprint::printMatrix(strm,q);
  return strm;
}


/*
  Prints the dual component group of G as a sep-separated list.
*/
std::ostream& printComponents(std::ostream& strm,
			      const realredgp::RealReductiveGroup& G,
			      const char* sep)
{
  const latticetypes::SmallBitVectorList& cr = G.dualComponentReps();
  basic_io::seqPrint(strm,cr.begin(),cr.end(),sep);

  return strm;
}

} // namespace testprint

} // namespace atlas
