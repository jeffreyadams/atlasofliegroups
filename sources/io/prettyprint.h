/*
  This is prettyprint.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef PRETTYPRINT_H  /* guard against multiple inclusions */
#define PRETTYPRINT_H

#include <iosfwd>
#include <vector>

#include "../Atlas.h"

/******** function declarations *********************************************/

namespace atlas {

namespace prettyprint {

template<unsigned int d>
std::ostream& prettyPrint(std::ostream&, const BitSet<d>&, size_t);

std::ostream& prettyPrint(std::ostream&, const BitMap&, size_t n = 0);

template<unsigned int dim>
std::ostream& prettyPrint(std::ostream&, const BitVector<dim>&);

template<unsigned int dim>
std::ostream& prettyPrint(std::ostream&, const std::vector<BitVector<dim> >&);

template<typename V>
std::ostream& printBasis(std::ostream&, const std::vector<V>&);

std::ostream& printDescentSet(std::ostream&, const RankFlags&, size_t,
			      const char* sep = ",", const char* pre = "{",
			      const char* post = "}");

std::ostream& printInRootBasis(std::ostream&, RootNbr, const RootSystem&);

std::ostream& printInRootBasis(std::ostream&,
			       const RootNbrSet&, const RootSystem&);

std::ostream& printInvolution(std::ostream&,
			      const TwistedInvolution&,
			      const TwistedWeylGroup&);

template<typename C>
std::ostream& printVector(std::ostream&, const std::vector<C>&,
			  unsigned long width = 2);
template<typename C>
std::ostream& printMatrix(std::ostream&, const matrix::Matrix_base<C>&,
			  unsigned long width = 4);

std::ostream& printRootList(std::ostream&,
			    const RootNbrList&,const RootDatum&,
			    const char* sep = "\n");

std::ostream& printCorootList(std::ostream&,
			      const RootNbrList&, const RootDatum&,
			      const char* sep = "\n");

std::ostream& printStatus(std::ostream&, const gradings::Status&, size_t);

std::ostream& printTitsElt(std::ostream&, const TitsElt&, const TitsGroup&);

std::ostream& printTorusType(std::ostream&, const WeightInvolution&);

std::ostream& printWeylElt(std::ostream&, const WeylElt&, const WeylGroup&);

std::ostream& printWeylList(std::ostream&, const WeylEltList&,
			    const WeylGroup&, const char* sep = ",",
			    const char* pre = "", const char* post = "");
}

}

#endif
