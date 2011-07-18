/*!
\file
  This is interactive_lietype.cpp
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "interactive_lietype.h"

#include <cassert>

#include <cstdio>
#include <iostream>
#include <sstream>

namespace atlas {

namespace {

  void ignoreSimpleLieType(input::InputBuffer&);
  const char* UnequalRankTypes = "ADET";

}

/*****************************************************************************

        Chapter I -- Functions declared in interactive_lietype.h

******************************************************************************/

namespace interactive_lietype {

/*!
  Synopsis: checks if a valid inner class for lt will be read from buf.
*/

bool checkInnerClass(input::InputBuffer& buf, const LieType& lt,
		     bool output)
{
  std::streampos pos = buf.tellg();
  std::string icl(lietype::innerClassLetters);
  std::string uer(UnequalRankTypes);

  for (size_t j = 0; j < lt.size(); ++j) {

    if (buf.peek() == EOF) {
      if (output)
	std::cerr << "too few inner class symbols" << std::endl;
      buf.reset(pos);
      return false;
    }

    char x = 0;
    buf >> x;
    if (icl.find_first_of(x) == std::string::npos) { // bad type
      if (output)
	std::cerr << "bad inner class symbol " << x
		  << " (should be one of " << lietype::innerClassLetters
		  << ")" << std::endl;
      buf.reset(pos);
      return false;
    }

    if (x == 'C') { // complex case
      SimpleLieType slt = lt[j];
      if (j == lt.size()-1 or lt[j+1] != slt) { // bad type
	if (output)
	  std::cerr << "bad inner class symbol C" << std::endl
		    << "(needs two identical consecutive types)" << std::endl;
	buf.reset(pos);
	return false;
      }
      else // skip next entry of lt
	++j;
    }

    if (x == 'u') { // we must have an unequal-rank inner class
      SimpleLieType slt = lt[j];
      lietype::TypeLetter t = slt.type();
      size_t l = slt.rank();
      if (uer.find_first_of(t) == std::string::npos or
	  (t == 'A' and l == 1) or
	  (t == 'E' and l != 6)) { // bad type
	if (output)
	  std::cerr << "sorry, bad inner class symbol u" << std::endl
		    << "(allowed only for types A_n (n>1), D_n, E6, and T_n)"
		    << std::endl;
	buf.reset(pos);
	return false;
      }

    }
  }

  buf.reset(pos);
  return true;
}


/*!
  Synopsis: checks if buf starts with a valid Lie type.

  A valid Lie type is a dot-separated non-empty string of entities of the form
  where X is a letter in the range [A-G] or T (for torus), and n is a number
  in the appropriate range for X. White between entities is ignored; reading
  terminates when after a pair Xn the next read operation on a character does
  not produce '.' (i.e., when the next non-white character is not a dot.)

  Return value is 0 for correct input, non-zero otherwise.
*/
bool checkLieType(input::InputBuffer& buf)
{
  std::streampos pos = buf.tellg();

  bool notDone = true;

  while (notDone) {
    if (not checkSimpleLieType(buf)) { // bad input
      buf.reset(pos);
      return false;
    }
    ignoreSimpleLieType(buf);
    char x = 0;
    buf >> x;
    if (x != '.')
      break;
  }

  buf.reset(pos);

  if (not checkTotalRank(buf)) { // bad rank
    std::cerr << "sorry, rank should not exceed " << constants::RANK_MAX
	      << std::endl;
    return false;
  }

  return true;
}


/*!
  Synopsis: checks if reading a SimpleLieType from buf will succeed.
*/
bool checkSimpleLieType(input::InputBuffer& buf)
{
  std::streampos pos = buf.tellg();

  lietype::TypeLetter x = 0;
  buf >> x;

  std::string tl(lietype::typeLetters);

  if (tl.find_first_of(x) == std::string::npos) { // bad type
    std::cerr << "sorry, bad type " << x
	      << " (should be one of " << lietype::typeLetters << ")"
	      << std::endl;
    buf.reset(pos);
    return false;
  }

  size_t l = 0;
  buf >> l;

  if (not lietype::checkRank(x,l)) { // bad rank
    printRankMessage(std::cerr,x);
    buf.reset(pos);
    return false;
  }

  buf.reset(pos);
  return true;
}


/*!
  Synopsis: checks that the total rank does not exceed RANK_MAX.

  Precondition: it has already been checked that successive read operations
  on buf will yield a valid simple type, optionally followed by a dot and
  a valid type.
*/
bool checkTotalRank(input::InputBuffer& buf)
{
  std::streampos pos = buf.tellg();

  size_t l_tot = 0;
  bool notDone = true;

  while (notDone) {
    lietype::TypeLetter x;
    size_t l;
    buf >> x;
    buf >> l;
    l_tot += l;
    if (l_tot > constants::RANK_MAX) {
      buf.reset(pos);
      return false;
    }
    buf >> x;
    if (x != '.')
      break;
  }

  buf.reset(pos);
  return true;
}


/*!
  Prints the message appropriate for a bad choice of rank for type x.
*/
std::ostream& printRankMessage(std::ostream& strm, lietype::TypeLetter x)
{
  const unsigned r=constants::RANK_MAX;
  switch (x) {
  case 'A':
    strm << "sorry, in type A the rank must be between 1 and " << r;
    break;
  case 'B':
    strm << "sorry, in type B the rank must be between 2 and " << r;
    break;
  case 'C':
    strm << "sorry, in type C the rank must be between 2 and " << r;
    break;
  case 'D':
    strm << "sorry, in type D the rank must be between 4 and " << r;
    break;
  case 'E':
    strm << "sorry, in type E the rank must be 6, 7 or 8";
    break;
  case 'F':
  case 'f':
    strm << "sorry, in type " << x << " the rank must be 4";
    break;
  case 'G':
  case 'g':
    strm << "sorry, in type " << x << " the rank must be 2";
    break;
  case 'T':
    strm << "sorry, in type T the rank must be between 1 and " << r;
    break;
  default: // cannot happen
    assert(false && "unexpected type in printRankMessage");
  }

  return strm << std::endl;
}


/*!
  Synopsis: reads an inner class type from buf.

  Precondition: checkInnerClass(ict,buf) returns true;

  Maps e ("equal rank") to c; maps "u" (unequal rank) to s except for type
  D_2n.
*/
void readInnerClass(InnerClassType& ict, input::InputBuffer& buf,
		    const LieType& lt)
{
  ict.clear();

  for (size_t j = 0; j < lt.size(); ++j) {
    lietype::TypeLetter x = 0;
    buf >> x;
    if (x == 'e')
      x = 'c';
    if (x == 'u')
      if (lt[j].type() != 'D' or lt[j].rank()%2 != 0)
	x = 's';
    ict.push_back(x);
    if (x == 'C')
      ++j;
  }
}


/*!
  Synopsis: reads the Lie type from buf.

  Precondition: it has been checked that the read operation will succeed:
  successive read operations will yield a valid simple Lie type, optionally
  followed by a dot and a valid Lie type. Reading ends on EOF or when a
  simple Lie type is not imediately followed by a dot.

  To normalize the occurrence of torus factors, this function expands Tn with
  n>1 to n copies of T1
*/
void readLieType(LieType& lt, input::InputBuffer& buf)
{
  bool read = true;

  while (read) {
    lietype::TypeLetter x;
    size_t l;
    buf >> x;
    buf >> l;
    if (x=='T')
      while (l-->0) lt.push_back(SimpleLieType('T',1));
    else
      lt.push_back(SimpleLieType(x,l));
    buf >> x;
    if (x != '.') {
      buf.unget(); // put character back
      break;
    }
  }
}

} // namespace interactive_lietype

/*****************************************************************************

                Chapter III -- Auxiliary functions

******************************************************************************/

namespace {


/*!
  Synopsis: takes a SimpleLieType off the buffer.

  Precondition: a SimpleLieType will be successfully read; otherwise the
  effect is undefined.

  Useful when reading forward for testing purposes
*/
void ignoreSimpleLieType(input::InputBuffer& buf)
{
  lietype::TypeLetter x;
  buf >> x;

  size_t l;
  buf >> l;
}

} // namespace

} // namespace atlas
