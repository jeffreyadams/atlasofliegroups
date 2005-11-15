/*
  This is kl_io.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#include <iomanip>
#include <iostream>

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

std::ostream& printKL(std::ostream& strm, const kl::KLContext& klc)

/*
  Synopsis: outputs the non-zero extremal kl polynomials from klc to strm.

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
    const ExtremalRow& e = klc.extremalRow(y);
    const KLRow& klr = klc.klRow(y);
    strm << std::setw(width) << y << ": ";
    bool first = true;
    for (size_t j  = 0; j < e.size(); ++j) {
      if (klc.isZero(klr[j])) {
	++zeroCount;
	continue;
      }
      if (first) {
	strm << std::setw(width) << e[j] << ": ";
	first = false;
      } else {
	strm << std::setw(width+2)<< ""
	     << std::setw(width) << e[j] << ": ";
      }
      printPol(strm,*klr[j],KLIndeterminate);
      strm << std::endl;
      ++count;
    }
    strm << std::endl;
  }

  strm << count + zeroCount << " extremal pairs" << std::endl;
  strm << zeroCount << " zero polynomials; "
       << count << " nonzero polynomials" << std::endl;

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
