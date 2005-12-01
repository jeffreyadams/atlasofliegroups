/*
  This is kgb_io.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include <iomanip>
#include <iostream>

#include "kgb_io.h"

#include "ioutils.h"
#include "kgb.h"
#include "prettyprint.h"

/*****************************************************************************

  Input/output functions for the kgb data structure, defined in 
  sources/kl/kgb.h

******************************************************************************/

namespace atlas {

/*****************************************************************************

        Chapter I -- Functions declared in kgb_io.h

  ... explain here when it is stable ...

******************************************************************************/

namespace kgb_io {

std::ostream& printKGB(std::ostream& strm, const kgb::KGB& kgb)

/*
  Synopsis: outputs the data from kgb to strm.

  Explanation: for each parameter, we output the cross-actions and 
  cayley-actions for each generator, the length, and the underlying root
  datum permutation (or rather, the corresponding Weyl group element).
  We use a '*' for undefined cayley actions.

  NOTE: this will print reasonably on 80 columns only for groups that are
  not too large (up to rank 4 or so). We haven't tried to go over to more
  sophisticated formatting for larger groups.
*/

{
  using namespace kgb;
  using namespace prettyprint;

  // compute maximal width of entry
  int width = ioutils::digits(kgb.size()-1,10ul);
  int lwidth = ioutils::digits(kgb.length(kgb.size()-1),10ul);
  const int pad = 2;

  for (size_t j = 0; j < kgb.size(); ++j) {
    strm << std::setw(width) << j << ":  ";

    // print cross actions
    for (size_t s = 0; s < kgb.rank(); ++s) {
      strm << std::setw(width+pad) << kgb.cross(s,j);
    }
    strm << std::setw(pad) << "";

    // print cayley actions
    for (size_t s = 0; s < kgb.rank(); ++s) {
      KGBElt z = kgb.cayley(s,j);
      if (z != UndefKGB)
	strm << std::setw(width+pad) << z;
      else
	strm << std::setw(width+pad) << '*';
    }
    strm << std::setw(pad) << "";

    // print status
    printStatus(strm,kgb.status(j),kgb.rank());
    strm << std::setw(pad) << "";

    // print length
    strm << std::setw(lwidth) << kgb.length(j);
    strm << std::setw(pad) << "";

    // print root datum involution
    printWeylElt(strm,kgb.involution(j),kgb.weylGroup());

    strm << std::endl;
  }

  return strm;
}

}

}
