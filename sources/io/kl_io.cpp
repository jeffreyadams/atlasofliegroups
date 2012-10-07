/*
  This is kl_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "kl_io.h"

#include <iomanip>
#include <iostream>

#include "kl.h"

#include "bruhat.h"	// |BruhatOrder|

#include "ioutils.h"	// |digits|
#include "prettyprint.h" // |printPol|

/*****************************************************************************

  Output functions for the |KLContext|, defined in sources/kl/kl.h

******************************************************************************/

namespace atlas {

namespace {

/******** constant definitions ***********************************************/

const char* KLIndeterminate = "q"; // name used on output for indeterminate

}

/*****************************************************************************

        Chapter I -- Functions declared in kl_io.h

******************************************************************************/

namespace kl_io {


// Print the non-zero Kazhdan-Lusztig-Vogan polynomials from klc to strm.
std::ostream& printAllKL
  (std::ostream& strm, const kl::KLContext& klc, Block& block)
{
  size_t count = 0;

  int width = ioutils::digits(klc.size()-1,10ul);
  int tab = 2;

  for (size_t y = 0; y < klc.size(); ++y)
  {

    strm << std::setw(width) << y << ": ";
    bool first = true;

    for (size_t x = 0; x <= y; ++x) {
      const kl::KLPol& pol = klc.klPol(x,y);
      if (pol.isZero())
	continue;
      if (first)
      {
	strm << std::setw(width) << x << ": ";
	first = false;
      }
      else
      {
	strm << std::setw(width+tab)<< ""
	     << std::setw(width) << x << ": ";
      }
      prettyprint::printPol(strm,pol,KLIndeterminate);
      strm << std::endl;
      ++count;
    }

    strm << std::endl;
  }

  strm << count << " nonzero polynomial" << (count==1 ? "" : "s")
       << std::flush;

  BruhatOrder& Bruhat=block.bruhatOrder(); // non-const!
  size_t comp = Bruhat.n_comparable();
  strm << ", and " << comp-count << " zero polynomial"
       << (comp-count==1 ? "" : "s")
       << ",\n at " << comp << " Bruhat-comparable "
       << (comp==1 ? "pair." : "pairs.") << std::endl;

  return strm;
}


// Print the primitive kl polynomials from klc to strm.
std::ostream& printPrimitiveKL
  (std::ostream& strm, const kl::KLContext& klc, Block& block)
{
  size_t count = 0;
  size_t zero_count = 0;
  size_t incomp_count = 0;

  int width = ioutils::digits(klc.size()-1,10ul);
  int tab = 2;

  BruhatOrder& bo=block.bruhatOrder(); // non-const!
  const poset::Poset& Bruhat=bo.poset(); // full poset is generated here

  for (size_t y = 0; y < klc.size(); ++y)
  {
    kl::PrimitiveRow e = klc.primitiveRow(y); // list of ALL primitive x's

    strm << std::setw(width) << y << ": ";
    bool first = true;
    for (size_t j  = 0; j < e.size(); ++j)
      if (Bruhat.lesseq(e[j],y))
      { // now |x=e[j]| is primitive for |y| and Bruhat-comparable
	++count;
	if ((klc.klPol(e[j],y).isZero()))
	  ++zero_count;
	if (first)
	{
	  strm << std::setw(width) << e[j] << ": ";
	  first = false;
	}
	else
	{
	  strm << std::setw(width+tab)<< ""
	       << std::setw(width) << e[j] << ": ";
	}

	prettyprint::printPol(strm,klc.klPol(e[j],y),KLIndeterminate);
	strm << std::endl;
      }
      else
      {
	assert(klc.klPol(e[j],y).isZero());
	++incomp_count;
      } // |for (j)|

    ++count; // count $P_{y,y}$
    if (not first)
      strm << std::setw(width+tab)<< "";
    strm << std::setw(width) << y << ": 1" << std::endl << std::endl;
  } // |for(y)|

  strm << count  << " Bruhat-comparable primitive "
       << (count==1 ? "pair" : "pairs")
       << ", of which " << zero_count << " ha" << (zero_count==1 ? "s" : "ve")
       << " null polynomial,\n and " << incomp_count
       << " incomparable primitive " << (incomp_count==1 ? "pair" : "pairs")
       << std::endl;

  return strm;
}


// Print the list of all distinct Kazhdan-Lusztig-Vogan polynomials in |klc|
std::ostream& printKLList(std::ostream& strm, const kl::KLContext& klc)
{
  const kl::KLStore& store = klc.polStore();
  std::vector<kl::KLPol> polList;

  // get polynomials, omitting Zero
  for (kl::KLIndex i=0; i<store.size(); ++i)
  {
    const kl::KLPol& r=store[i];
    if (not r.isZero()) polList.push_back(r);
  }

  std::sort(polList.begin(),polList.end(),polynomials::compare<kl::KLCoeff>);

  for (size_t j = 0; j < polList.size(); ++j)
  {
    prettyprint::printPol(strm,polList[j],KLIndeterminate);
    strm << std::endl;
  }

  return strm;
}


// Print the mu-coefficients from |klc| to |strm|.
std::ostream& printMu(std::ostream& strm, const kl::KLContext& klc)
{
  int width = ioutils::digits(klc.size()-1,10ul);

  for (size_t y = 0; y < klc.size(); ++y) {
    const kl::MuRow& mrow = klc.muRow(y);
    strm << std::setw(width) << y << ": ";
    for (size_t j = 0; j < mrow.size(); ++j) {
      if (j>0)
	strm << ",";
      strm << "(" << mrow[j].first << "," << mrow[j].second << ")";
    }
    strm << std::endl;
  }

  return strm;
}

} // namespace kl_io

} // namsespace atlas
