/*!
\file
  This is interactive_lietype.cpp
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "interactive_lietype.h"

#include <cassert>

#include <iostream>
#include <sstream>

/*****************************************************************************

  ... explain here when it is stable ....

******************************************************************************/

namespace atlas {

namespace {

  void ignoreSimpleLieType(input::InputBuffer&);
  const char* UnequalRankTypes = "ADE";

}

/*****************************************************************************

        Chapter I -- Functions declared in interactive_lietype.h

  ... explain here when it is stable ...

******************************************************************************/

namespace interactive_lietype {

bool checkInnerClass(input::InputBuffer& buf, const lietype::LieType& lt,
		     bool output)

/*!
  Synopsis: checks if a valid inner class for lt will be read from buf.
*/

{
  using namespace lietype;

  std::streampos pos = buf.tellg();
  std::string icl(innerClassLetters);
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
		  << " (should be one of " << innerClassLetters << ")" 
		  << std::endl;
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
      TypeLetter t = type(slt);
      size_t l = rank(slt);
      if (uer.find_first_of(t) == std::string::npos or
	  (t == 'A' and l == 1) or
	  (t == 'E' and l != 6)) { // bad type
	if (output)
	  std::cerr << "sorry, bad inner class symbol u" << std::endl
		    << "(allowed only for types A_n, n > 1, D_n, and E6)" 
		    << std::endl;
	buf.reset(pos);
	return false;
      }

    }
  }

  buf.reset(pos);
  return true;
}

bool checkLieType(input::InputBuffer& buf)

/*!
  Synopsis: checks if buf starts with a valid Lie type. 

  A valid Lie type is a dot-separated non-empty string of entities of the form 
  where X is a letter in the range [A-G] or T (for torus), and n is a number
  in the appropriate range for X. White between entities is ignored; reading 
  terminates when after a pair Xn the next read operation on a character does
  not produce '.' (i.e., when the next non-white character is not a dot.)

  Return value is 0 for correct input, non-zero otherwise.
*/

{
  using namespace constants;

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
    std::cerr << "sorry, rank should not exceed " << RANK_MAX << std::endl;
    return false;
  }

  return true;
}

bool checkSimpleLieType(input::InputBuffer& buf)

/*!
  Synopsis: checks if reading a SimpleLieType from buf will succeed.
*/

{
  using namespace constants;
  using namespace lietype;

  std::streampos pos = buf.tellg();

  TypeLetter x = 0;
  buf >> x;

  std::string tl(typeLetters);

  if (tl.find_first_of(x) == std::string::npos) { // bad type
    std::cerr << "sorry, bad type " << x 
	      << " (should be one of " << typeLetters << ")" << std::endl;
    buf.reset(pos);
    return false;
  }
 
  size_t l = 0;
  buf >> l;

  if (not checkRank(x,l)) { // bad rank
    printRankMessage(std::cerr,x);
    buf.reset(pos);
    return false;
  }

  buf.reset(pos);
  return true;
}

bool checkTotalRank(input::InputBuffer& buf)

/*!
  Synopsis: checks that the total rank does not exceed RANK_MAX.

  Precondition: it has already been checked that successive read operations
  on buf will yield a valid simple type, optionally followed by a dot and
  a valid type.
*/

{
  using namespace constants;
  using namespace lietype;

  std::streampos pos = buf.tellg();

  size_t l_tot = 0;
  bool notDone = true;

  while (notDone) {
    TypeLetter x;
    size_t l;
    buf >> x;
    buf >> l;
    l_tot += l;
    if (l_tot > RANK_MAX) {
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

std::ostream& printRankMessage(std::ostream& strm, lietype::TypeLetter x)

/*!
  Prints the message appropriate for a bad choice of rank for type x.
*/

{
  using namespace constants;
  using namespace std;

  switch (x) {
  case 'A':
    cerr << "sorry, in type A the rank must be between 1 and " << RANK_MAX
	 << endl;
    break;
  case 'B':
    cerr << "sorry, in type B the rank must be between 2 and " << RANK_MAX
	 << endl;
    break;
  case 'C':
    cerr << "sorry, in type C the rank must be between 2 and " << RANK_MAX
	 << endl;
    break;
  case 'D':
    cerr << "sorry, in type D the rank must be between 4 and " << RANK_MAX
	 << endl;
    break;
  case 'E':
    cerr << "sorry, in type E the rank must be 6, 7 or 8" << endl;
    break;
  case 'F':
  case 'f':
    cerr << "sorry, in type " << x << " the rank must be 4" << endl;
    break;
  case 'G':
  case 'g':
    cerr << "sorry, in type " << x << " the rank must be 2" << endl;
    break;
  case 'T':
    cerr << "sorry, in type T the rank must be between 1 and " << RANK_MAX
	 << endl;
    break;
  default: // cannot happen
    assert(false && "unexpected type in printRankMessage");
  }

  return strm;
}

void readInnerClass(lietype::InnerClassType& ict, input::InputBuffer& buf,
		    const lietype::LieType& lt)

/*!
  Synopsis: reads an inner class type from buf.

  Precondition: checkInnerClass(ict,buf) returns true;

  Maps e ("equal rank") to c; maps "u" (unequal rank) to s except for type
  D_2n.
*/

{
  using namespace lietype;

  ict.clear();

  for (size_t j = 0; j < lt.size(); ++j) {
    TypeLetter x = 0;
    buf >> x;
    if (x == 'e')
      x = 'c';
    if (x == 'u')
      if (type(lt[j]) != 'D' or rank(lt[j])%2 != 0)
	x = 's';
    ict.push_back(x);
    if (x == 'C')
      ++j;
  }

  return;
}

void readLieType(lietype::LieType& lt, input::InputBuffer& buf)

/*!
  Synopsis: reads the Lie type from buf.

  Precondition: it has been checked that the read operation will succeed:
  successive read operations will yield a valid simple Lie type, optionally
  followed by a dot and a valid Lie type. Reading ends on EOF or when a
  simple Lie type is not imediately followed by a dot.
*/

{
  using namespace constants;
  using namespace lietype;

  bool read = true;

  while (read) {
    TypeLetter x;
    size_t l;
    buf >> x;
    buf >> l;
    lt.push_back(SimpleLieType(x,l));
    buf >> x;
    if (x != '.') {
      buf.unget(); // put character back
      break;
    }
  }

  return;
}

}

/*****************************************************************************

                Chapter III -- Auxiliary functions

******************************************************************************/

namespace {

void ignoreSimpleLieType(input::InputBuffer& buf)

/*!
  Synopsis: takes a SimpleLieType off the buffer.

  Precondition: a SimpleLieType will be successfully read; otherwise the
  effect is undefined.

  Useful when reading forward for testing purposes
*/

{
  using namespace constants;
  using namespace lietype;

  TypeLetter x;
  buf >> x;

  size_t l;
  buf >> l;

  return;
}

}

}
