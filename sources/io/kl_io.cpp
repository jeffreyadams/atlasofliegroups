/*
  This is kl_io.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include <iomanip>
#include <iostream>
#include <set>

#include "kl_io.h"

#include "ioutils.h"
#include "kl.h"
#include "klsupport.h"
#include "prettyprint.h"

/*****************************************************************************

  Input/output functions for the KLContext data structure, defined in 
  sources/kl/kl.h

******************************************************************************/

namespace atlas {

namespace {

/******** constant definitions ***********************************************/

const char* KLIndeterminate = "q";

}

/*****************************************************************************

        Chapter I -- Functions declared in kl_io.h

  ... explain here when it is stable ...

******************************************************************************/

namespace kl_io {

std::ostream& printAllKL(std::ostream& strm, kl::KLContext& klc)

/*
  Synopsis: outputs the non-zero kl polynomials from klc to strm.

*/

{
  using namespace ioutils;
  using namespace kl;
  using namespace klsupport;
  using namespace prettyprint;

  size_t count = 0;
  size_t zeroCount = 0;

  int width = digits(klc.size()-1,10ul);
  int tab = 2;

  for (size_t y = 0; y < klc.size(); ++y) {

    strm << std::setw(width) << y << ": ";
    bool first = true;

    for (size_t x = 0; x <= y; ++x) {
      const KLPol& pol = klc.klPol(x,y);
      if (pol.isZero()) {
	++zeroCount;
	continue;
      }
      if (first) {
	strm << std::setw(width) << x << ": ";
	first = false;
      } else {
	strm << std::setw(width+tab)<< ""
	     << std::setw(width) << x << ": ";
      }
      printPol(strm,pol,KLIndeterminate);
      strm << std::endl;
      ++count;
    }

    strm << std::endl;
  }

  strm << count + zeroCount << " pairs" << std::endl;
  strm << zeroCount << " zero polynomials; "
       << count << " nonzero polynomials" << std::endl;

  return strm;
}

std::ostream& printPrimitiveKL(std::ostream& strm, const kl::KLContext& klc)

/*!
  \brief Outputs the non-zero primitive kl polynomials from klc to strm.
*/

{
  using namespace ioutils;
  using namespace kl;
  using namespace klsupport;
  using namespace prettyprint;

  size_t count = 0;
  //  size_t zeroCount = 0;

  int width = digits(klc.size()-1,10ul);
  int tab = 2;

  for (size_t y = 0; y < klc.size(); ++y) {

    const PrimitiveRow& e = klc.primitiveRow(y);
    const KLRow& klr = klc.klRow(y);
    strm << std::setw(width) << y << ": ";
    bool first = true;
    for (size_t j  = 0; j < e.size(); ++j) {
    //  if (klc.isZero(klr[j])) {
    //    ++zeroCount;
    // 	continue;
    //      }
      if (first) {
	strm << std::setw(width) << e[j] << ": ";
	first = false;
      } else {
	strm << std::setw(width+tab)<< ""
	     << std::setw(width) << e[j] << ": ";
      }
      printPol(strm,*klr[j],KLIndeterminate);
      strm << std::endl;
      ++count;
    }
    strm << std::endl;
  }

  strm << count  << " primitive pairs with nonzero polynomial." << std::endl;
   //  strm << zeroCount << " zero polynomials; "
   //   << count << " nonzero polynomials" << std::endl;

  return strm;
}

std::ostream& printKLList(std::ostream& strm, kl::KLContext& klc)

/*
  Synopsis: outputs the list of all distinct Kazhdan-Lusztig-Vogan
  polynomials for the block.
*/

{
  using namespace kl;
  using namespace prettyprint;

  const std::set<KLPol>& store = klc.polStore();
  std::vector<KLPol> polList;

  std::copy(store.begin(),store.end(),back_inserter(polList));
  std::sort(polList.begin(),polList.end(),polynomials::compare<KLCoeff>);

  for (size_t j = 1; j < polList.size(); ++j) {
    printPol(strm,polList[j],KLIndeterminate);
    strm << std::endl;
  }

  return strm;
}

std::ostream& printMu(std::ostream& strm, const kl::KLContext& klc)

/*
  Synopsis: outputs the mu-coefficients from klc to strm.
*/

{
  using namespace ioutils;
  using namespace kl;

  int width = digits(klc.size()-1,10ul);

  for (size_t y = 0; y < klc.size(); ++y) {
    const MuRow& mrow = klc.muRow(y);
    strm << std::setw(width) << y << ": ";
    for (size_t j = 0; j < mrow.size(); ++j) {
      if (j)
	strm << ",";
      const MuData& md = mrow[j];
      strm << "(" << md.first << "," << md.second << ")";
    }
    strm << std::endl;
  }

  return strm;
}

}

}
