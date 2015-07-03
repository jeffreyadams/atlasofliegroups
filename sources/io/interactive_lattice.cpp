/*
  This is interactive_lattice.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "interactive_lattice.h"

#include <iostream>
#include <sstream>
#include <cstdio>
#include <cassert>
#include <cstdlib>

#include "arithmetic.h" // |lcm|
#include "ratvec.h"	// |RatWeight|
#include "interactive.h" // |common_input| variable

#include "matreduc.h"	// |diagonalise|

#include "input.h" // input and history buffers

/*****************************************************************************

 This module deals with functions specific to the user choosing a sublattice.

******************************************************************************/

namespace atlas {

namespace {

  enum GeneratorError { NoError=0, FormatError, BadDenominator,
			NegDenominator, TooFew };

  GeneratorError checkGenerator(input::InputBuffer& buf, size_t& i,
				LatticeCoeff& d,
				const CoeffList& u);

  int_Matrix makeOrthogonal (const RatWeightList& rwl, size_t r);

  std::ostream& printCenter(std::ostream&, const CoeffList&);

  RatWeight readGenerator(size_t n_gen,
			  LatticeCoeff d,
			  input::InputBuffer& buf);

  input::HistoryBuffer kernelgen_input_buffer;
}

/*****************************************************************************

        Chapter I -- Functions declared in interactive_lattice.h

******************************************************************************/

namespace interactive_lattice {

/*
  Gets the generators of X/Q, where Q is the root lattice, from the user.

  It throws an InputError if the interaction with the user does not conclude
  successfully. In that case, d_rwl is not modified.
*/
int getGenerators(RatWeightList& d_rwl, const CoeffList& u)
  throw(error::InputError)
{
  RatWeightList rwl;

  std::string genString;
  interactive::common_input() >> genString;

  if (genString.find("sc") == 0) // match must be a start of string
    return 1; // code for simply connected
  if (genString.find("ad") == 0)
    return 2; // code for adjoint

  std::cout
    << "torsion subgroup of the center of the simply connected group is:"
    << std::endl;

  printCenter(std::cout,u) << std::endl;

  std::cout << "enter kernel generators, one per line" << std::endl;
  std::cout << "(ad for adjoint, sc for simply connected, ? to abort):"
	    << std::endl;

  {
    input::InputBuffer& ib=kernelgen_input_buffer;

    while (true) {
      ib.getline("");
      if (hasQuestionMark(ib))
	throw error::InputError();
      genString.clear();
      ib >> genString;
      if (genString.empty()) // done
	break;
      if (genString.find("sc") == 0)
	return 1; // code for simply connected
      if (genString.find("ad") == 0)
	return 2; // code for adjoint

      ib.reset();
      size_t i;
      LatticeCoeff d;
      switch (checkGenerator(ib,i,d,u))
      {
      case NoError:
	break;
      case FormatError:
	std::cerr << "bad format in entry #" << i+1
		  << " (should be of the form a/";
	if (u[i])
	  std::cerr << u[i] << ")" << std::endl;
	else
	  std::cerr << "b, b > 0)" << std::endl;
	std::cerr << "bad input line --- ignored" << std::endl;
	continue;
      case BadDenominator:
	std::cerr << "denominator in entry #" << i+1
		  << " should be " << u[i] << std::endl
		  << "bad input line --- ignored" << std::endl;
	continue;
      case NegDenominator:
	std::cerr << "denominator in entry #" << i+1
		  << " should be positive" << std::endl
		  << "bad input line --- ignored" << std::endl;
	continue;
      case TooFew:
	std::cerr << "too few valid entries" << std::endl
		  << "bad input line --- ignored" << std::endl;
	continue;
      }
      rwl.push_back(readGenerator(u.size(),d,ib));
    }
  }

  d_rwl.swap(rwl);
  return 0; // "normal" exit
}

/*
  Gets the lattice interactively from the user.

  This works as follows. Initially, |root_lattice_basis| is basis of the
  weight lattice adjusted to the root lattice, namely such that the latter is
  spanned by the multiples in |root_invf| of the corresponding basis elements.
  As a matrix it is in block form according to the simple factors, and for
  each torus factor the corresponding factor in |root_invf| is zero.

  Guided (only) by the non-unit values among those factors, the user has to
  specify "kernel generators" (in |getGenerators|), sequences of rational
  numbers that are transformed into rational coweights; this implicitly
  specifies a full rank (and therefore finite index) sublattice of weights
  taking integral values on all those coweights, and which sublattice contains
  the root lattice (by limitation of the choice of the generators roots
  automatically satisfy the integrality condition).

  Generators for the nontrivial part of this sublattice are computed by
  |makeOrthogonal| into |q|. More precisely, extracting the elements of
  |root_lattice_basis| that are not in the root lattice we obtain a direct
  factor of the weight lattice with the integrality condition trivially
  satisfied on the complementary factor; the columns of the square matrix |q|
  express the linear combinations of the extracted elements that also satisfy
  the integrality condition. To get from this a real sublattice basis, it
  suffices to replace the extracted elements in |root_lattice_basis| by those
  linear combinations.

  This function may throw an |InputError| if the interaction with the user
  fails. In that case the arguments are not modified. Also, if the user
  prefers answering "ad" or "sc" rather than giving any generators, we pass
  this condition as a return code without modifying anything.
*/
int getLattice(const CoeffList& root_invf, WeightList& root_lattice_basis)
  throw(error::InputError)
{
  size_t r = root_lattice_basis.size(); // full rank of root datum
  assert(root_invf.size()==r);

  CoeffList u; // non-unit invariant factors
  WeightList lb; // "local basis" corresponding to non-units
  for (size_t i=0; i<r; ++i)
    if (root_invf[i] != 1)
    {
      u.push_back(root_invf[i]);
      lb.push_back(root_lattice_basis[i]);
    }

  // get generators of character group
  RatWeightList rwl;  // generator list, each of size |u.size()|
  int code=getGenerators(rwl,u);    // input them; might throw an InputError

  if (code>0)
    return code; // bypass computation if user typed "sc" or "ad"

  // make basis elements corresponding to those central elements

  LatticeMatrix q = makeOrthogonal(rwl,u.size()); // local function, see below

  // convert |lb| columns to linear combinations of them according to |q|
  LatticeMatrix lin_comb(lb,r); // matrix with columns |lb|
  lin_comb *= q; // replace them by linear combinations

  // make actual basis, inserting columns |lin_comb| into |root_lattice_basis|
  for (size_t i=0,j=0; i<r; ++i)
    if (root_invf[i] != 1)
    {
      lin_comb.get_column(root_lattice_basis[i],j); // |rlb[i] = lc.column(j)|
      ++j;
    }
  return 0; // normal exit
}

} // |namespace interactive_lattice|

/*****************************************************************************

        Chapter II -- Functions local to interactive_lattice.cpp

******************************************************************************/

namespace {

/*
  Synposis: checks if buf contains data compatible with u.

  Precondition: buf should contain a comma-separated list, with one entry for
  each member of u (extra entries are ignored). The entries should be of the
  form a/b, with b = u[i] if u[i] > 0, or with arbitrary b>0 otherwise.

  Explanation: the element a/b (with b = u[i] if u[i] > 0) represents the
  element exp(2i.pi.a/b) in a one-dimensional torus factor.

  In case of success, the l.c.m. of the various b's and u[i]'s is put in |d|.
  In case of failure, the entry number where the error occurs is put in |r|.

  NOTE: no overflow checking is done on d.
*/
GeneratorError checkGenerator(input::InputBuffer& buf, size_t& r,
			      LatticeCoeff& d,
			      const CoeffList& u)
{
  std::streampos pos = buf.tellg();
  unsigned long lc = 1;

  for (size_t i=0; i<u.size(); ++i) {

    if (buf.peek() == EOF) {
      buf.reset(pos);
      return TooFew;
    }

    LatticeCoeff a;
    buf >> a;
    char x=0;
    buf >> x;
    if (x != '/') {
      r = i;
      buf.reset(pos);
      return FormatError;
    }

    LatticeCoeff b = 0;
    buf >> b;

    // check denominator
    if (u[i] == 0) {
      if (b <= 0) { // negative denominator error
	r = i;
	buf.reset(pos);
	return NegDenominator;
      }
    } else if (b != u[i]) { // bad denominator error
      r = i;
      buf.reset(pos);
      return BadDenominator;
    }

    // update lowest common denominator
    if (a)
      lc = arithmetic::lcm(lc,b);

    if (i+1<u.size()) { // next non-white character should be a comma
      char x = 0;
      buf >> x;
      if (x != ',') {
	buf.reset(pos);
	return TooFew;
      }
    }
  }

  buf.reset(pos);
  d = lc;

  return NoError;
}


/*
  Synopsis: puts in |q| a basis for the lattice "orthogonal" to |rwl|.

  Precondition: each of the elements of |rwl| is a rational weight, size |r|.

  Explanation: the elements of rwl are interpreted as elements of finite
  order in the torus (more precisely, their order divides their denominator.)
  So "orthogonal" means having integral pairing with the rational vector.

  Algorithm: collect in |m| all the denominator vectors after bringing
  everything to a common denominator |d|; then we have the problem of finding
  a basis for the lattice on the dual side that under pairing takes all those
  vectors into $d.\Z$. Write $row.m.col=D$ with |row| and |col| invertible
  such that $D$ is diagonal (as for the Smith normal form, but without
  divisibility condition on the entries of $D$). Then the rows of $row$ form a
  dual basis to a basis $b$ such that the column span of $m$ (or $m*col$) is
  generated by their multiples $D_{i,i} b_i$. It follows that in order for a
  dual vector to send $b_i$ into $d.\Z$, its coefficient of $row[i]$ should be
  a multiple of $d/gcd(d,D_{i,i})$; our "orthogonal" lattice is spanned by the
  multiples $(d/gcd(d,D_{i,i}))*row[i]$. We return a matrix with those columns

  NOTE: this is a sloppy implementation; we don't worry about overflow.
*/
LatticeMatrix makeOrthogonal (const RatWeightList& rwl, size_t r)
{
  // find common denominator
  arithmetic::Denom_t d = 1;
  for (size_t i=0; i<rwl.size(); ++i)
    d = arithmetic::lcm(d,rwl[i].denominator());

  LatticeMatrix m(r,rwl.size()); // matrix of numerators

  for (size_t j=0; j<m.numColumns(); ++j)
  {
    Weight num(rwl[j].numerator().begin(),rwl[j].numerator().end()); // convert
    m.set_column(j,num*LatticeCoeff(d/rwl[j].denominator()));
  }

  // find "orthogonal" basis
  int_Matrix row,col;
  CoeffList factor = matreduc::diagonalise(m,row,col);
  // now |row.row(i)*(m*col).col(j)| is $\delta_{i,j}factor[i]$ (for all $i,j$)

  // initialise |result| to dual basis to |adapted_basis(m)|
  int_Matrix result = row.transposed();

  // integrality requires multiple $d/\gcd(d,factor[j])$ of |result.col(j)|
  for (size_t j=0; j<factor.size(); ++j)
  {
    long mult = arithmetic::div_gcd(d,std::abs(factor[j])); // multiplier
    result.columnMultiply(j,mult);
  }
  return result;
}


/*
  Print the center of the simply connected group.

  Precondition: u contains the necessary data: the orders of a natural set of
  generators for the center of the derived group, and a zero for each torus
  factor.
*/
std::ostream& printCenter(std::ostream& strm, const CoeffList& u)
{
  for (size_t i=0; i<u.size(); ++i)
  {
    if (u[i]==0)
      strm << "Q/Z";
    else
      strm << "Z/" << u[i];
    if (i<u.size()-1)
      strm << ".";
  }

  return strm;
}


/*!
  \brief: returns generator from |buf|, with |n_gen| entries, denominator |d|

*/
RatWeight readGenerator(size_t n_gen,
                                      LatticeCoeff d,
				      input::InputBuffer& buf)
{
  int_Vector numer(n_gen);
  char dummy; // bit bucket for separator characters
  LatticeCoeff b;
  for (size_t i=0; i<n_gen; ++i)
  {
    buf >> numer[i] >> dummy >> b >> dummy;// read fraction, skipping "/", ","
    numer[i] *= d/b; // because denominators are forced to |d|
  }
  return RatWeight(numer,d);
}

} // |namespace|

} // |namespace atlas|
