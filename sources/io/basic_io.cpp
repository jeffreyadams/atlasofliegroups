/*
  This is basic_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "basic_io.h"

#include <iostream>
#include <stdexcept>

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

  Originally this was without separators, with the following
  NOTE: this is satisfactory only if the rank is < 10; otherwise, use
  prettyprint::printWeylWord, that will give dot-separated lists.

  Unfortunately this operator is often called automatically, by template
  functions, so there is no place to insert that distinction into the code.
  Therefore we have inserted commas here. MvL
*/

{
  using namespace weyl;

  for (size_t j = 0; j < w.size(); ++j) {
    unsigned a = w[j]+1;
    strm << a;
    if (j+1<w.size()) strm << ',';
  }

  return strm;
}

unsigned long long read_var_bytes(unsigned int n,std::istream& in)
{ switch(n)
  { case 1: return read_bytes<1>(in);
    case 2: return read_bytes<2>(in);
    case 3: return read_bytes<3>(in);
    case 4: return read_bytes<4>(in);
    case 5: return read_bytes<5>(in);
    case 6: return read_bytes<6>(in);
    case 7: return read_bytes<7>(in);
    case 8: return read_bytes<8>(in);
  default: throw std::runtime_error("Illegal read_var_bytes");
  }
}

void put_int (unsigned int val, std::ostream& out) { write_bytes<4>(val,out); }
void write_bytes(unsigned int n, unsigned long long val, std::ostream& out)
{ switch(n)
  { case 1: return write_bytes<1>(val,out);
    case 2: return write_bytes<2>(val,out);
    case 3: return write_bytes<3>(val,out);
    case 4: return write_bytes<4>(val,out);
    case 5: return write_bytes<5>(val,out);
    case 6: return write_bytes<6>(val,out);
    case 7: return write_bytes<7>(val,out);
    case 8: return write_bytes<8>(val,out);
    default: throw std::runtime_error("Illegal write__bytes");
  }
}

} // namespace basic_io

} // namespace atlas
