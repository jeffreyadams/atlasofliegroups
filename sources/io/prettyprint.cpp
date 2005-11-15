/*
  This is prettyprint.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#include "prettyprint.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "abelian.h"
#include "basic_io.h"
#include "gradings.h"
#include "lattice.h"
#include "lietype.h"
#include "rootdata.h"
#include "tits.h"
#include "tori.h"
#include "weyl.h"

/*****************************************************************************

        Chapter I -- Functions declared in prettyprint.h

  ... fill in here when it is stable ...

******************************************************************************/

namespace atlas {

namespace prettyprint {

std::ostream& prettyPrint(std::ostream& strm, const bitmap::BitMap& b,
			  size_t n)

/*
  Synopsis: outputs the first values of the bitmap left-to-right, on a single 
  line
*/

{
  if (n == 0)
    n = b.capacity();

  for (size_t j = 0; j < n; ++j)
    if (b.isMember(j))
      strm << "1";
    else
      strm << "0";

  return strm;
}

std::ostream& prettyPrint(std::ostream& strm, const abelian::GroupType& type)

/*
  Synopsis: outputs the group type as a bracket-enclosed, comma-separated
  list of integers.
*/

{  
  using namespace abelian;
  using namespace basic_io;

  GroupType::const_iterator first = type.begin();
  GroupType::const_iterator last = type.end();
  seqPrint(strm,first,last,",","[","]");

  return strm;
}

std::ostream& printCoroot(std::ostream& strm, const rootdata::RootNbr& j, 
			  const rootdata::RootDatum& rd)

/*
  Synopsis: prints coroot #j in the lattice basis.
*/

{
  using namespace basic_io;

  return strm << rd.coroot(j);
}

std::ostream& printCorootList(std::ostream& strm, const rootdata::RootList& r, 
			      const rootdata::RootDatum& rd, const char* sep)

/*
  Synopsis: prints the coroots in the list in the lattice basis, by default
  as one per line.
*/

{
  for (size_t j = 0; j < r.size(); ++j) {
    printCoroot(strm,r[j],rd);
    if (j+1 < r.size())
      strm << sep;
  }

  return strm;
}

std::ostream& printInRootBasis(std::ostream& strm, const rootdata::RootSet& r, 
			       const rootdata::RootDatum& rd)

/*
  Synopsis: outputs the set of roots contained in r to strm, expressed in root 
  coordinates.
*/

{
  using namespace rootdata;

  RootSet::iterator r_end = r.end();
  RootList rl(r.begin(),r_end);

  RootList::const_iterator first = rl.begin();
  RootList::const_iterator last = rl.end();

  return printInRootBasis(strm,first,last,rd);
}

std::ostream& printInRootBasis(std::ostream& strm, rootdata::RootNbr n, 
			       const rootdata::RootDatum& rd)

/*
  Outputs root #n to strm in the root coordinates.

  NOTE : for now this goes through the general printInRootBasis function.
*/

{
  return printInRootBasis(strm,rd.root(n),rd);
}

std::ostream& printInRootBasis(std::ostream& strm, 
			       const latticetypes::Weight& v, 
			       const rootdata::RootDatum& rd)

/*
  Outputs the element v, assumed to lie in the root lattice, to strm in
  the root coordinates.
*/

{
  using namespace basic_io;
  using namespace latticetypes;

  Weight vrb;
  const Weight* vPtr = &v;
  const Weight* nvPtr = vPtr+1;
  rd.toRootBasis(vPtr,nvPtr,&vrb);

  strm << vrb;

  return strm;
}

std::ostream& printRoot(std::ostream& strm, const rootdata::RootNbr& j, 
			const rootdata::RootDatum& rd)

/*
  Synopsis: prints coroot #j in the lattice basis.
*/

{
  using namespace basic_io;

  return strm << rd.root(j);
}

std::ostream& printRootList(std::ostream& strm, const rootdata::RootList& r, 
			    const rootdata::RootDatum& rd, const char* sep)

/*
  Synopsis: prints the roots in the list in the lattice basis, by default
  as one per line.
*/

{
  for (size_t j = 0; j < r.size(); ++j) {
    printRoot(strm,r[j],rd);
    if (j+1 < r.size())
      strm << sep;
  }

  return strm;
}

std::ostream& printStatus(std::ostream& strm, const gradings::Status& gs, 
			  size_t rank)

/*
  Synopsis: prints the status flags.

  Precondition: there are rank valid fields in gs;

  Explanation: the output is in the format [xxx...] where each entry is
  C for complex, c for (imaginary) compact, n for (imaginary) noncompact,
  and r for real.
*/

{
  using namespace gradings;

  strm << "[";

  for (size_t s = 0; s < rank; ++s) {
    switch (gs[s]) {
    case Status::Complex:
      strm << "C";
      break;
    case Status::ImaginaryCompact:
      strm << "c";
      break;
    case Status::ImaginaryNoncompact:
      strm << "n";
      break;
    case Status::Real:
      strm << "r";
      break;
    default: // should not happen!
      strm << "*";
      break;
    }
  }

  strm << "]";

  return strm;
}

std::ostream& printTitsElt(std::ostream& strm, const tits::TitsElt& a, 
			   const tits::TitsGroup& N)

/*
  Synopsis: outputs a in the format w[...], where the Weyl part w is output
  as a reduced expression, and the torus part as a bitvector.
*/

{
  using namespace basic_io;

  printWeylElt(strm,a.w(),N.weylGroup()) << "[" << a.t() << "]";

  return strm;
}

std::ostream& printTorusType(std::ostream& strm, const tori::RealTorus& T)

/*
  Synopsis: outputs the type of the real torus.

  Explanation: T(R) is of the form (R^x)^p.(U(1))^q.(C^x)^r.
*/

{
  strm << "split: ";
  strm << T.splitRank();

  strm << "; compact: ";
  strm << T.compactRank();

  strm << "; complex: ";
  strm << T.complexRank();

  return strm;
}

std::ostream& printWeylElt(std::ostream& strm, const weyl::WeylElt& w, 
			   const weyl::WeylGroup& W)

/*
  Synopsis: outputs w as a reduced expression.
*/

{
  using namespace basic_io;
  using namespace weyl;

  WeylWord ww;
  W.out(ww,w);
  strm << ww;

  return strm;
}

std::ostream& printWeylList(std::ostream& strm, const weyl::WeylEltList& wl,
			    const weyl::WeylGroup& W, const char* sep, 
			    const char* pre, const char* post)

/*
  Synopsis: outputs the list of WeylElts as words in the outer representation,
  with the given separator, prefix and postfix.
*/

{
  using namespace basic_io;
  using namespace weyl;

  std::vector<WeylWord> wwl(wl.size());

  for (size_t j = 0; j < wl.size(); ++j) {
    WeylWord& ww = wwl[j];
    W.out(ww,wl[j]);
  }

  std::vector<WeylWord>::const_iterator first = wwl.begin();
  std::vector<WeylWord>::const_iterator last = wwl.end();
  seqPrint(strm,first,last,sep,pre,post);

  return strm;
}

}

}
