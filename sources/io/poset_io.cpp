/*
  This is poset_io.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  For license information see the LICENSE file
*/

#include <iomanip>
#include <iostream>

#include "poset_io.h"

#include "bitmap.h"
#include "poset.h"

/*****************************************************************************

  Input/output functions for the Poset type and related types, defined in 
  sources/utilities/poset.h

******************************************************************************/

namespace atlas {

/*****************************************************************************

        Chapter I -- Functions declared in poset_io.h

  ... explain here when it is stable ...

******************************************************************************/

namespace poset_io {

std::ostream& printPoset(std::ostream& strm, const poset::Poset& p)

/* not implemented yet */

{
  return strm;
}

std::ostream& printSymmetricPoset(std::ostream& strm, 
				  const poset::SymmetricPoset& p)

/*
  Synopsis: outputs the poset on strm.

  Explanation: prints * for comparable pairs, . for uncomparable ones, and
  a blank on the diagonal.
*/

{
  using namespace bitmap;

  size_t uncomparable = 0;
  size_t comparable = 0;

  for (size_t i = 0; i < p.size(); ++i) {
    const BitMap& row = p.row(i);
    for (size_t j = 0; j < i; ++j) {
      if (row.isMember(j)) {
	strm << '*';
	++comparable;
      } else {
	strm << '.';
	++uncomparable;
      }
    }
    strm << ' ';
    for (size_t j = i+1; j < p.size(); ++j) {
      if (row.isMember(j)) {
	strm << '*';
	++comparable;
      } else {
	strm << '.';
	++uncomparable;
      }
    }
    strm << std::endl;
  }

  std::cerr << "there were " << comparable << " comparable pairs, and "
	    << uncomparable << " uncomparable pairs" << std::endl;

  return strm;
}

}

}
