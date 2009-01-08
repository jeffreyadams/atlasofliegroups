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

******************************************************************************/

namespace poset_io {

std::ostream& printPoset(std::ostream& strm, const poset::Poset& p)

{
  for (size_t i=0; i<p.size(); ++i)
  {
    for (size_t j=i; j<p.size(); ++j)
      strm << (p.lesseq(i,j) ? '*' : '.');
    strm << std::endl;
  }
  return strm;
}

}

}
