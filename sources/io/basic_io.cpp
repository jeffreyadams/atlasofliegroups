/*
  This is basic_io.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "basic_io.h"

#include <iostream>

#include "lattice.h"
#include "lietype.h"

/*****************************************************************************

        Chapter I -- operator<< functions

  NOTE: there are more definitions as templates.

  ... fill in here when it is stable ...

******************************************************************************/

namespace atlas {

namespace basic_io {

/******** from abelian *******************************************************/

std::ostream& operator<< (std::ostream& strm, const abelian::GrpArr& a)

{
  using namespace abelian;

  GrpArr::const_iterator first = a.begin();
  GrpArr::const_iterator last = a.end();

  return seqPrint(strm,first,last,",","[","]");
}

/******** from latticetypes **************************************************/

std::ostream& operator<< (std::ostream& strm, 
			  const latticetypes::LatticeElt& v)

/*
  Synopsis: output of a lattice element.

  Uses the template seqPrint. It is output as a bracket-enclosed, 
  comma-separated list.
*/

{
  using namespace latticetypes;

  LatticeElt::const_iterator first = v.begin();
  LatticeElt::const_iterator last = v.end();

  return seqPrint(strm, first, last, ",", "[", "]");
}

/******** from lietype *******************************************************/

std::ostream& operator<< (std::ostream& strm, 
			  const lietype::SimpleLieType& slt)

/*
  Synopsis: outputs a simple Lie type.
*/

{
  using namespace lietype;

  TypeLetter x = type(slt);
  size_t l = rank(slt);

  return strm << x << l;
}

std::ostream& operator<< (std::ostream& strm, const lietype::LieType& lt)

/*
  Synopsis: outputs the LieType as a dot-separated string of simple Lie types.
*/

{
  using namespace lietype;

  LieType::const_iterator first = lt.begin();
  LieType::const_iterator last = lt.end();

  return seqPrint(strm,first,last,".");
}

/******** from weyl **********************************************************/

std::ostream& operator<< (std::ostream& strm, const weyl::WeylWord& w)

/*
  Synopsis: outputs w as a string of digits

  NOTE: this is satisfactory only if the length is < 10; otherwise, use
  prettyprint::printWeylWord, that will give dot-separated lists.
*/

{
  using namespace weyl;

  for (size_t j = 0; j < w.size(); ++j) {
    unsigned a = w[j]+1;
    strm << a;
  }

  return strm;
}

}

}
