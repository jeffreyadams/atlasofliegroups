/*
  This is prettyprint.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef PRETTYPRINT_H  /* guard against multiple inclusions */
#define PRETTYPRINT_H

#include <iosfwd>
#include <vector>

#include "abelian_fwd.h"
#include "bitmap_fwd.h"
#include "bitset_fwd.h"
#include "gradings_fwd.h"
#include "latticetypes_fwd.h"
#include "polynomials_fwd.h"
#include "rootdata_fwd.h"
#include "weyl_fwd.h"
#include "tits_fwd.h"
#include "tori_fwd.h"

#include "lietype.h"

/******** function declarations *********************************************/

namespace atlas {

namespace prettyprint {

std::ostream& prettyPrint(std::ostream&, const bitmap::BitMap&,
			  size_t n = 0);

std::ostream& prettyPrint(std::ostream&, const abelian::GroupType&);

template<size_t d>
std::ostream& prettyPrint(std::ostream&, const bitset::BitSet<d>&, size_t);

template<size_t dim>
std::ostream& prettyPrint(std::ostream&,
			  const bitvector::BitVector<dim>&, size_t);

template<size_t dim>
std::ostream& prettyPrint(std::ostream&,
			  const std::vector<bitvector::BitVector<dim> >&,
			  size_t);

template<size_t dim>
std::ostream& prettyPrint(std::ostream&,
			  const bitvector::BitVector<dim>&, size_t);

template<typename V>
std::ostream& printBasis(std::ostream&, const std::vector<V>&);

std::ostream& printCoroot(std::ostream&, const rootdata::RootNbr&,
			  const rootdata::RootDatum&);

std::ostream& printCorootList(std::ostream&, const rootdata::RootList&,
			      const rootdata::RootDatum&,
			      const char* sep = "\n");

std::ostream& printDescentSet(std::ostream&, const bitset::RankFlags&, size_t,
			      const char* sep = ",", const char* pre = "{",
			      const char* post = "}");

std::ostream& printInRootBasis(std::ostream&, rootdata::RootNbr,
			       const rootdata::RootDatum&);

std::ostream& printInRootBasis(std::ostream&, const rootdata::RootSet&,
			       const rootdata::RootDatum&);

std::ostream& printInvolution(std::ostream&, const weyl::TwistedInvolution&,
			      const weyl::WeylGroup&);

template<typename C>
std::ostream& printVector(std::ostream&, const std::vector<C>&,
			  unsigned long width = 2);
template<typename C>
std::ostream& printMatrix(std::ostream&, const matrix::Matrix<C>&,
			  unsigned long width = 4);

template<typename C>
std::ostream& printMonomial(std::ostream&, C, polynomials::Degree,
			    const char*);

template<typename C>
std::ostream& printPol(std::ostream&, const polynomials::Polynomial<C>&,
		       const char*);

std::ostream& printRoot(std::ostream&, const rootdata::RootNbr&,
			const rootdata::RootDatum&);

std::ostream& printRootList(std::ostream&, const rootdata::RootList&,
			    const rootdata::RootDatum&,
			    const char* sep = "\n");

std::ostream& printStatus(std::ostream&, const gradings::Status&, size_t);

std::ostream& printTitsElt(std::ostream&, const tits::TitsElt&,
			   const tits::TitsGroup&);

std::ostream& printTorusType(std::ostream&, const tori::RealTorus&);

std::ostream& printWeylElt(std::ostream&, const weyl::WeylElt&,
			   const weyl::WeylGroup&);

std::ostream& printWeylList(std::ostream&, const weyl::WeylEltList&,
			    const weyl::WeylGroup&, const char* sep = ",",
			    const char* pre = "", const char* post = "");
}

}

#include "prettyprint_def.h"

#endif
