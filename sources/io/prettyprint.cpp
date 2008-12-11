/*
  This is prettyprint.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "prettyprint.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "basic_io.h"
#include "bitmap.h"
#include "gradings.h"
#include "lattice.h"
#include "lietype.h"
#include "rootdata.h"
#include "tits.h"
#include "tori.h"
#include "weyl.h"

/*****************************************************************************

        Chapter I -- Functions declared in prettyprint.h

******************************************************************************/

namespace atlas {

namespace prettyprint {


/*
  Synopsis: outputs the first values of the bitmap left-to-right, on a single
  line
*/
std::ostream& prettyPrint(std::ostream& strm, const bitmap::BitMap& b,
			  size_t n)
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


/*
  Synopsis: outputs the first values of the bitmap left-to-right, on a single
  line
*/
std::ostream& printCoroot(std::ostream& strm, const rootdata::RootNbr& j,
			  const rootdata::RootDatum& rd)
{
  return strm << rd.coroot(j);
}


/*
  Synopsis: prints the coroots in the list in the lattice basis, by default
  as one per line.
*/
std::ostream& printCorootList(std::ostream& strm, const rootdata::RootList& r,
			      const rootdata::RootDatum& rd, const char* sep)
{
  for (size_t j = 0; j < r.size(); ++j) {
    printCoroot(strm,r[j],rd);
    if (j+1 < r.size())
      strm << sep;
  }

  return strm;
}


/*
  Synopsis: prints the descent set d to strm.

  Here rank is the number of significant bits in d; the output format is
  pre * sep * ... * post, where the * are the bits in d, output as their
  bitposition starting from 1.
*/
std::ostream& printDescentSet(std::ostream& strm, const bitset::RankFlags& d,
			      size_t rank, const char* sep, const char* pre,
			      const char* post)
{
  strm << pre;

  bool first = true;

  for (size_t s = 0; s < rank; ++s)
    if (d.test(s)) {
      if (first)
	first = false;
      else
	strm << sep;
      strm << s+1;
    }

  strm << post;

  return strm;
}


/*
  Outputs root #n to strm in the root coordinates.
*/
std::ostream& printInRootBasis(std::ostream& strm, rootdata::RootNbr n,
			       const rootdata::RootDatum& rd)
{
  return strm << rd.inSimpleRoots(n);
}

/*
  Synopsis: outputs the set of roots contained in r to strm, expressed in root
  coordinates.
*/
std::ostream& printInRootBasis(std::ostream& strm, const rootdata::RootSet& r,
			       const rootdata::RootDatum& rd)
{
  latticetypes::WeightList rl;

  for (rootdata::RootSet::iterator i=r.begin(); i(); ++i)
    rl.push_back(rd.inSimpleRoots(*i));

  basic_io::seqPrint(strm,rl.begin(),rl.end(),"\n","","\n");

  return strm;
}

/*
  Synopsis: outputs an expression for the twisted involution.

  Precondition: w is a (twisted) involution.

  Symbols are to be interpreted from right to left as operations performed on
  an initially empty twisted involution; if the number |s| is followed by a
  period it means left multiplication by a (twisted-commuting) generator |s|,
  if it is followed by an 'x' (for cross action) it means twiseted conjugation
  by |s|.
*/
std::ostream& printInvolution(std::ostream& strm,
			      const weyl::TwistedInvolution& tw,
			      const weyl::TwistedWeylGroup& W)
{
  std::vector<signed char> dec=W.involution_expr(tw);
  for (size_t i=0; i<dec.size(); ++i)
    if (dec[i]>=0) strm << static_cast<char>('1'+dec[i]) << '^';
    else strm << static_cast<char>('1'+~dec[i]) << 'x';

  return strm;
}


/*
  Synopsis: prints coroot #j in the lattice basis.
*/

std::ostream& printRoot(std::ostream& strm, const rootdata::RootNbr& j,
			const rootdata::RootDatum& rd)
{
  return strm << rd.root(j);
}


/*
  Synopsis: prints the roots in the list in the lattice basis, by default
  as one per line.
*/
std::ostream& printRootList(std::ostream& strm, const rootdata::RootList& r,
			    const rootdata::RootDatum& rd, const char* sep)
{
  for (size_t j = 0; j < r.size(); ++j) {
    printRoot(strm,r[j],rd);
    if (j+1 < r.size())
      strm << sep;
  }

  return strm;
}


/*
  Synopsis: prints the status flags.

  Precondition: there are rank valid fields in gs;

  Explanation: the output is in the format [xxx...] where each entry is
  C for complex, c for (imaginary) compact, n for (imaginary) noncompact,
  and r for real.
*/
std::ostream& printStatus(std::ostream& strm, const gradings::Status& gs,
			  size_t rank)
{
  strm << '[';

  for (size_t s = 0; s < rank; ++s)
  {
    if (s>0) strm<<',';
    switch (gs[s])
    {
    case gradings::Status::Complex:
      strm << "C";
      break;
    case gradings::Status::ImaginaryCompact:
      strm << "c";
      break;
    case gradings::Status::ImaginaryNoncompact:
      strm << "n";
      break;
    case gradings::Status::Real:
      strm << "r";
      break;
    }
  }

  strm << ']';

  return strm;
}

std::ostream& printTitsElt(std::ostream& strm, const tits::TitsElt& a,
			   const tits::TitsGroup& Tg)
{
  prettyPrint(strm,Tg.left_torus_part(a));
  printWeylElt(strm,a.w(),Tg.weylGroup());

  return strm;
}


/*
  Synopsis: outputs the type of the real torus.

  Explanation: T(R) is of the form (R^x)^p.(U(1))^q.(C^x)^r.
*/
std::ostream& printTorusType(std::ostream& strm, const tori::RealTorus& T)
{
  strm << "split: ";
  strm << T.splitRank();

  strm << "; compact: ";
  strm << T.compactRank();

  strm << "; complex: ";
  strm << T.complexRank();

  return strm;
}


/*
  Synopsis: outputs w as a reduced expression.
*/
std::ostream& printWeylElt(std::ostream& strm, const weyl::WeylElt& w,
			   const weyl::WeylGroup& W)
{
  strm << W.word(w);
  return strm;
}

/*
  Synopsis: outputs the list of WeylElts as words in the outer representation,
  with the given separator, prefix and postfix.
*/
std::ostream& printWeylList(std::ostream& strm, const weyl::WeylEltList& wl,
			    const weyl::WeylGroup& W, const char* sep,
			    const char* pre, const char* post)
{
  std::vector<weyl::WeylWord> wwl(wl.size());

  for (size_t j = 0; j < wl.size(); ++j)
    wwl[j]=W.word(wl[j]);

  basic_io::seqPrint(strm,wwl.begin(),wwl.end(),sep,pre,post);

  return strm;
}

} // namespace prettyprint

} // namespace atlas
