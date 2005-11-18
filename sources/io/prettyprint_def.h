/*
  This is prettyprint_def.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#include <iomanip>

#include "basic_io.h"
#include "lattice.h"
#include "rootdata.h"

namespace atlas {

namespace prettyprint {

namespace {

// private declarations

template<typename I>
std::ostream& printInRootBasis(std::ostream&, 
			       const I&, 
			       const I&, 
			       const rootdata::RootDatum&,
			       const rootdata::RootNbr&, 
			       const char* sep = ",", 
			       const char* pre = "", 
			       const char* post = "");


template<typename I>
std::ostream& printInRootBasis(std::ostream&, 
			       const I&, 
			       const I&, 
			       const rootdata::RootDatum&,
			       const latticetypes::Weight&, 
			       const char* sep = ",", 
			       const char* pre = "", 
			       const char* post = "");

}

}


/*****************************************************************************

        Chapter I -- Template functions defined in prettyprint.h

  ... explain here when it is stable ...

******************************************************************************/

namespace prettyprint {

template<size_t d> 
std::ostream& prettyPrint(std::ostream& strm, const bitset::BitSet<d>& b,  
			  size_t n)

/*
  Prints the n first bits of v on strm left-to-right.
*/

{
  for (size_t j = 0; j < n; ++j)
    if (b.test(j))
      strm << "1";
    else
      strm << "0";

  return strm;
}

template<size_t dim> 
std::ostream& prettyPrint(std::ostream& strm, 
			  const bitvector::BitVector<dim>& v,  size_t n)

/*
  Prints the n first bits of v on strm in a "vector-like" format.
*/

{
  std::vector<int> vi;

  for (size_t j = 0; j < n; ++j)
    if (v.test(j))
      vi.push_back(1);
    else
      vi.push_back(0);

  seqPrint(strm,vi.begin(),vi.end(),",","(",")");

  return strm;
}

template<size_t dim> 
std::ostream& prettyPrint(std::ostream& strm, 
			  const std::vector<bitvector::BitVector<dim> >& a, 
			  size_t n)

/*
  Pretty-prints a list of bitvectors, one per line.
*/

{
  for (size_t j = 0; j < a.size(); ++j) {
    prettyPrint(strm,a[j],n);
    strm << std::endl;
  }

  return strm;
}

template<typename V>
std::ostream& printBasis(std::ostream& strm, const std::vector<V>& b)

/*
  This is a function to output a basis in matrix form. It is assumed that
  V is a vector type, and that the elements of b all have the same size.
*/

{
  if (b.size() == 0) // do nothing
    return strm;

  size_t dim = b[0].size();

  for (size_t j = 0; j < dim; ++j) { // output one line of the matrix
    for (size_t i = 0; i < b.size(); ++i)
      strm << std::setw(4) << b[i][j];
    strm << std::endl;
  }

  return strm;
}

template<typename I>
std::ostream& printInRootBasis(std::ostream& strm, 
			       const I& first, 
			       const I& last, 
			       const rootdata::RootDatum& rd,
			       const char* sep, 
			       const char* pre, 
			       const char* post)

/*
  Dispatching function for the printInRootBasis template. Basically we
  want to deal with two sorts of iterators : those weight-valued ones
  and RootNbr-valued ones. We call the two specialised functions
  accordingly. C++ magic!
*/

{
  typedef typename std::iterator_traits<I>::value_type VT;

  VT dummy;
  return printInRootBasis(strm,first,last,rd,dummy,sep,pre,post);
}

template<typename C>
std::ostream& printMonomial(std::ostream& strm, C c, polynomials::Degree d,
			    const char* x)

/*
  Synopsis: prints out the monomial c.x^d.

  Preconditions: c is non-zero;

  Explanation: c and d are printed only if non-one, except in degree zero where
  c is always printed; x is the name of the indeterminate.
*/

{
  if (d == 0) // output c
    strm << c;
  else {
    if (c != 1ul)
      strm << c;
    strm << x;
    if (d > 1)
      strm << "^" << d;
  }

  return strm;
}

template<typename C>
std::ostream& printPol(std::ostream& strm, const polynomials::Polynomial<C>& p,
		       const char* x)

/*
  Synopsis: outputs the polynomial p on strm.

  It is assumed that operator<< exists for the coefficient type. The string
  x is the name of the indeterminate. It is assumed that C is an unsigned
  type, so that no minus signs are required.

  The output format is tex-like. Terms are printed only if non-zero, 
  coefficients and exponents are printed only if non-one.
*/

{
  if (p.isZero()) {
    strm << "0";
    return strm;
  }

  bool firstterm = true;

  for (size_t j = p.degree()+1; j;) {
    --j;
    if (p[j]) {
      if (firstterm)
	firstterm = false;
      else
	strm << "+";
      printMonomial(strm,p[j],j,x);
    }
  }

  return strm;
}

/*****************************************************************************

        Chapter II -- Local functions

  ... explain here when it is stable ...

  Of course, the functions are not really local, as this module will be
  included in many others!

******************************************************************************/

namespace {

template<typename I>
std::ostream& printInRootBasis(std::ostream& strm, 
			       const I& first, 
			       const I& last, 
			       const rootdata::RootDatum& rd,
			       const latticetypes::Weight&, 
			       const char* sep, 
			       const char* pre, 
			       const char* post)

/*
  In this function we assume that I is an iterator type whose value_type is
  Weight; furthermore we assume that when i ranges from first to last, *i
  is in the root lattice. This function then outputs the corresponding
  sequence of weights written in the simple root basis, outputting through
  seqPrint with the corresponding separator, prefix and postfix.

  NOTE: the Weight& argument is just a placeholder to make the signature
  unique.
*/

{
  using namespace latticetypes;
  using namespace rootdata;

  WeightList rl;

  rd.toRootBasis(first,last,back_inserter(rl));
  seqPrint(strm,rl.begin(),rl.end(),sep,pre,post);

  return strm;
}

template<typename I>
std::ostream& printInRootBasis(std::ostream& strm, 
			       const I& first, 
			       const I& last, 
			       const rootdata::RootDatum& rd,
			       const rootdata::RootNbr&, 
			       const char* sep, 
			       const char* pre, 
			       const char* post)

/*
  In this function we assume that I is an iterator type whose value_type is
  RootNbr; it is assumed that these numbers refer to the list rd.d_roots. 
  This function then outputs the corresponding sequence of roots written in 
  the simple root basis, outputting through seqPrint with the corresponding 
  separator, prefix and postfix.

  NOTE: the RootNbr& argument is just a placeholder to make the signature 
  unique.
*/

{
  using namespace basic_io;
  using namespace latticetypes;

  WeightList rl;
  
  rootdata::RootIterator firstR(rd.beginRoot(),first);
  rootdata::RootIterator lastR(rd.beginRoot(),last);

  rd.toRootBasis(firstR,lastR,back_inserter(rl));

  WeightList::const_iterator rl_begin = rl.begin();
  WeightList::const_iterator rl_end = rl.end();
  seqPrint(strm,rl_begin,rl_end,sep,pre,post);

  return strm;
}

}

template<typename C>
std::ostream& printMatrix(std::ostream& strm, const matrix::Matrix<C>& m,
			  unsigned long width)

/*
  Outputs the matrix to a stream. It is assumed that operator << is defined
  for C, and that C is "small" (in particular, has no newlines in its output.)
*/

{
  using namespace matrix;

  for (size_t i = 0; i < m.columnSize(); ++i) {
    for (size_t j = 0; j < m.rowSize(); ++j) {
      strm << std::setw(width);
      strm << m(i,j);
    }
    strm << std::endl;
  }

  return strm;
}

}

}
