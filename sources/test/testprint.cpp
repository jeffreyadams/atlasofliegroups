/*
  This is testprint.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "testprint.h"

#include <iostream>
#include <iomanip>
#include <iterator>
#include <sstream>

#include "basic_io.h"
#include "cartanclass.h"
#include "innerclass.h"
#include "output.h"
#include "dynkin.h"
#include "ioutils.h"
#include "lattice.h"
#include "lietype.h"
#include "prettyprint.h"
#include "realredgp.h"
#include "rootdata.h"
#include "size.h"
#include "tori.h"
#include "weylsize.h"


/*****************************************************************************

  This module contains a number of print commands vor viewing parts of the data
  which would normally be hidden away. It is intended for testing purposes.

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

******************************************************************************/

namespace testprint {


/*
  This function outputs the following data about the rootdatum :

    - Cartan matrix;
    - set of simple roots and coroots;
    - if not semisimple : bases for radical and coradical;
    - positive roots and coroots;
*/
std::ostream& print(std::ostream& strm, const RootDatum& rd)
{
  strm << "cartan matrix :" << std::endl;
  int_Matrix q(rd.semisimple_rank(),rd.semisimple_rank());
  for (size_t j = 0; j < rd.semisimple_rank(); ++j)
    for (size_t i = 0; i < rd.semisimple_rank(); ++i)
      q(i,j) = rd.cartan(i,j);
  prettyprint::printMatrix(strm,q) << std::endl;

  strm << "root basis :" << std::endl;
  WeightList r_rb(rd.beginSimpleRoot(),rd.endSimpleRoot());
  prettyprint::printBasis(strm,r_rb) << std::endl;

  strm << "coroot basis :" << std::endl;
  CoweightList r_crb(rd.beginSimpleCoroot(),rd.endSimpleCoroot());
  prettyprint::printBasis(strm,r_crb) << std::endl;

  if (not rd.isSemisimple()) { // print radical and coradical bases
    strm << "radical basis :" << std::endl;
    CoweightList r_rad(rd.beginRadical(),rd.endRadical()); // generate list
    prettyprint::printBasis(strm,r_rad) << std::endl;

    strm << "coradical basis :" << std::endl;
    WeightList r_crad(rd.beginCoradical(),rd.endCoradical());
    prettyprint::printBasis(strm,r_crad) << std::endl;
  }

  strm << "positive roots :" << std::endl;
  WeightList r_pr(rd.beginPosRoot(),rd.endPosRoot()); // generate actual roots
  basic_io::seqPrint(strm,r_pr.cbegin(),r_pr.cend(),"\n","","")
    << std::endl << std::endl;

  strm << "positive coroots :" << std::endl;
  CoweightList r_pcr(rd.beginPosCoroot(),rd.endPosCoroot()); // generate coroots
  basic_io::seqPrint(strm,r_pcr.cbegin(),r_pcr.cend(),"\n","","")
    << std::endl << std::endl;

  return strm;
}


/*
  In this function, |rb| should contain a basis of some root system; we are
  outputting the matrix with M(i,j)=root(rb[i]).scalarProduct(coroot(rb[j]))
*/
std::ostream& printCartanMatrix(std::ostream& strm,
				const RootNbrList& rb,
				const RootSystem& rs)
{
  return prettyprint::printMatrix(strm,rs.Cartan_matrix(rb));
}


// Print the dual component group of G as a sep-separated list.
std::ostream& printComponents(std::ostream& strm,
			      const RealReductiveGroup& G,
			      const char* sep)
{
  const SmallBitVectorList& cr = G.dualComponentReps();
  basic_io::seqPrint(strm,cr.begin(),cr.end(),sep);

  return strm;
}

} // |namespace testprint|

} // |namespace atlas|
