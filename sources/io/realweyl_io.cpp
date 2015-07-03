/*
  This is realweyl_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "realweyl_io.h"

#include <iostream>
#include <sstream>

#include "basic_io.h"
#include "prettyprint.h"
#include "realweyl.h"	// |RealWeyl| class

/*  Chapter 0 -- string constants used for printing */

namespace atlas {

  namespace realweyl_io {

    namespace
{

  /* The following constants are used below in calls to printWeylList, and
     may be set differently to suit your format requirements. For instance
     to get a bracketed list of parenthesized Weyl words separated by
     commas, set pre="[(", sep="), (", and post=")]\n". See also the
     operator << for WeylWord values in basic_io.cpp, which controls
     output of individual Weyl words.
  */

  static const char* const pre="";    // printed before Weyl elements
  static const char* const sep="\n";  // printed between Weyl elements
  static const char* const post="\n"; // printed after Weyl elements

    }
  }
}


namespace atlas {

namespace realweyl_io {

/*  Chapter 0 -- string constants used for printing */

enum which_group { block_stab, real_W, dual_real_W };

namespace {

inline std::string Weyl_type(size_t n,const LieType& t)
{
  if (n==0) return "trivial";
  std::ostringstream os;
  os << "a Weyl group of type " << t;
  return os.str();
}

inline std::string two_type(size_t n)
{
  if (n==0) return "trivial";
  std::ostringstream os;
  os << "an elementary abelian 2-group of rank " << n;
  return os.str();
}

void common_print(which_group kind,
		  std::ostream& strm,
		  const realweyl::RealWeyl& rw,
		  const realweyl::RealWeylGenerators& rwg)
{
  strm << "W^C is " << (rw.numComplex()>0 ? "isomorphic to " : "")
       << Weyl_type(rw.numComplex(),rw.complexType()) << std::endl;

  if (kind==dual_real_W)
    strm << "W^I is "
	 << Weyl_type(rw.numImaginary(),rw.imaginaryType())
	 << std::endl;
  else
    strm << (kind==real_W ? "A" : "A_i") << " is "
	 << two_type(rw.numImaginaryR()) << std::endl
	 << "W_ic is "
	 << Weyl_type(rw.numImaginaryCompact(),rw.imaginaryCompactType())
	 << std::endl;

  if (kind==real_W)
    strm << "W^R is "
	 << Weyl_type(rw.numReal(),rw.realType())
	 << std::endl;
  else
    strm << (kind==dual_real_W ? "A" : "A_r") << " is "
	 << two_type(rw.numRealR()) << std::endl
	 << "W_rc is "
	 << Weyl_type(rw.numRealCompact(),rw.realCompactType())
	 << std::endl;

  strm << std::endl;


  if (rw.numComplex()!=0)
    prettyprint::printWeylList(strm << "generators for W^C:" << std::endl,
			       rwg.complex(),rwg.weylGroup(),sep,pre,post);

  if (kind==dual_real_W)
  {
    if (rw.numImaginary()!=0)
      prettyprint::printWeylList(strm << "generators for W^I:" << std::endl,
				 rwg.imaginary(),rwg.weylGroup(),sep,pre,post);
  }
  else
  {
    if (rw.numImaginaryR()!=0)
      prettyprint::printWeylList
	(strm << "generators for A" << (kind==real_W ? "" : "_i") << std::endl,
	 rwg.imaginaryR(),rwg.weylGroup(),sep,pre,post);

    if (rw.numImaginaryCompact()!=0)
      prettyprint::printWeylList(strm << "generators for W_ic:" << std::endl,
				 rwg.imaginaryCompact(),rwg.weylGroup(),
				 sep,pre,post);
  }

  if (kind==real_W)
  {
    if (rw.numReal()!=0)
      prettyprint::printWeylList(strm << "generators for W^R:" << std::endl,
				 rwg.real(),rwg.weylGroup(),sep,pre,post);
  }
  else
  {
    if (rw.numRealR()!=0)
      prettyprint::printWeylList
	(strm << "generators for A" << (kind==real_W ? "" : "_r") << std::endl,
	 rwg.realR(),rwg.weylGroup(),sep,pre,post);

    if (rw.numRealCompact()!=0)
      prettyprint::printWeylList(strm << "generators for W_rc:" << std::endl,
				 rwg.realCompact(),rwg.weylGroup(),
				 sep,pre,post);
  }
}



} // anonymous namespace

/*****************************************************************************

        Chapter 2 -- Functions declared in realweyl_io.h

******************************************************************************/


/*
  Synopsis: prints the stabilizer of a representation in this W-orbit
  for the cross-action.

  Explanation: this is the intersection of the real and dual real Weyl groups.
*/
std::ostream& printBlockStabilizer(std::ostream& strm,
				   const realweyl::RealWeyl& rw,
				   const realweyl::RealWeylGenerators& rwg)
{
  strm << "block stabilizer is W^C.((A_i.W_ic) x (A_r.W_rc)), where:"
       << std::endl;

  common_print(block_stab,strm,rw,rwg);
  return strm;
}

/*
  Synopsis: prints the real Weyl group.
*/
std::ostream& printRealWeyl(std::ostream& strm, const realweyl::RealWeyl& rw,
			    const realweyl::RealWeylGenerators& rwg)
{
  strm << "real weyl group is W^C.((A.W_ic) x W^R), where:" << std::endl;

  common_print(real_W,strm,rw,rwg);
  return strm;
}

/*
  Synopsis: prints the dual real Weyl group.
*/
std::ostream& printDualRealWeyl(std::ostream& strm,
				const realweyl::RealWeyl& rw,
				const realweyl::RealWeylGenerators& rwg)
{
  strm << "dual real weyl group is W^C.(W^I x (A.W_rc)), where:"
       << std::endl;

  common_print(dual_real_W,strm,rw,rwg);
  return strm;
}


} // |namespace realweyl_io|

} // |namespace atlas|
