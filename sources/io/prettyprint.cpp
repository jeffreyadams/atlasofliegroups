/*
  This is prettyprint.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "prettyprint.h"

#include <iostream>
#include <iomanip>
#include <sstream>

#include "bitmap.h"
#include "polynomials.h"

#include "gradings.h"   // |gradings::Status|
#include "rootdata.h"	// |RootSystem|

#include "tits.h"

#include "basic_io.h"	// |operator<<| for vectors, |seqPrint|

/*****************************************************************************

        Chapter I -- Functions declared in prettyprint.h

******************************************************************************/

namespace atlas {

namespace prettyprint {


/*
  Synopsis: outputs the first values of the bitmap left-to-right, on a single
  line
*/
std::ostream& prettyPrint(std::ostream& strm, const BitMap& b,
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


// Prints the n first bits of v on strm left-to-right.
template<size_t d>
std::ostream& prettyPrint(std::ostream& strm, const BitSet<d>& b,
			  size_t n)
{
  for (size_t j = 0; j < n; ++j)
    if (b.test(j))
      strm << "1";
    else
      strm << "0";

  return strm;
}


// Prints the bits of |v| on |strm| in a "vector-like" format.
template<size_t dim>
std::ostream& prettyPrint(std::ostream& strm, const BitVector<dim>& v)
{
  set::EltList vi;

  for (size_t i = 0; i < v.size(); ++i)
    vi.push_back(v[i]?1:0);

  const set::EltList& vir=vi; // this type already used with |seqPrint|
  basic_io::seqPrint(strm,vir.begin(),vir.end(),",","(",")");

  return strm;
}


// Pretty-prints a list of bitvectors, one per line.
template<size_t dim>
std::ostream& prettyPrint(std::ostream& strm,
			  const std::vector<BitVector<dim> >& a)
{
  for (size_t i = 0; i<a.size(); ++i)
  {
    prettyPrint(strm,a[i]);
    strm << std::endl;
  }

  return strm;
}

// This is a function to output a basis as columns in denuded matrix form.
// The component type |V| is indexable, and elements of |b| have same size.
template<typename V>
std::ostream& printBasis(std::ostream& strm, const std::vector<V>& b)
{
  if (b.size() == 0) // do nothing, needed because |b[0].size()| inexistent
    return strm;

  size_t dim = b[0].size();

  for (size_t i = 0; i < dim; ++i) // row index is index into each |b[j]|
  {
    for (size_t j = 0; j < b.size(); ++j)
      strm << std::setw(4) << b[j][i];
    strm << std::endl;
  }

  return strm;
}


/*
  Synopsis: prints the descent set d to strm.

  Here rank is the number of significant bits in d; the output format is
  pre * sep * ... * post, where the * are the bits in d, output as their
  bitposition starting from 1.
*/
std::ostream& printDescentSet(std::ostream& strm, const RankFlags& d,
			      size_t rank, const char* sep, const char* pre,
			      const char* post)
{
  strm << pre;

  bool first = true;

  for (size_t s = 0; s < rank; ++s)
    if (d.test(s))
    {
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
std::ostream& printInRootBasis(std::ostream& strm, RootNbr n,
			       const RootSystem& rs)
{
  return strm << rs.root_expr(n);
}

/*
  Synopsis: outputs the set of roots contained in r to strm, expressed in root
  coordinates.
*/
std::ostream& printInRootBasis(std::ostream& strm, const RootNbrSet& r,
			       const RootSystem& rs)
{
  int_VectorList rl; rl.reserve(r.size());

  for (RootNbrSet::iterator it=r.begin(); it(); ++it)
    rl.push_back(rs.root_expr(*it));

  const int_VectorList& rlr=rl; // share instantiation with |mainmode::roots_f|
  basic_io::seqPrint(strm,rlr.begin(),rlr.end(),"\n","","\n");

  return strm;
}

/*
  Synopsis: prints the roots in the list in the lattice basis, by default
  as one per line.
*/
std::ostream& printRootList(std::ostream& strm, const RootNbrList& r,
			    const RootDatum& rd, const char* sep)
{
  for (size_t i=0; i<r.size(); ++i)
  {
    strm << rd.root(r[i]);
    if (i+1 < r.size())
      strm << sep;
  }

  return strm;
}

/*
  Synopsis: prints the coroots in the list in the lattice basis, by default
  as one per line.
*/
std::ostream& printCorootList(std::ostream& strm, const RootNbrList& r,
			      const RootDatum& rd, const char* sep)
{
  for (size_t j=0; j<r.size(); ++j) {
    strm << rd.coroot(r[j]);
    if (j+1 < r.size())
      strm << sep;
  }

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
			      const TwistedInvolution& tw,
			      const TwistedWeylGroup& W)
{
  weyl::InvolutionWord dec=W.involution_expr(tw);
  for (size_t i=0; i<dec.size(); ++i)
    if (dec[i]>=0) strm << static_cast<char>('1'+dec[i]) << '^';
    else strm << static_cast<char>('1'+~dec[i]) << 'x';

  return strm;
}

template<typename C>
std::ostream& printVector(std::ostream& strm, const std::vector<C>& v,
			  unsigned long width)
{
  for (size_t i = 0; i < v.size(); ++i)
    strm << (i==0 ? '[' : ',') << std::setw(width) << v[i];

  strm << " ]";
  return strm;
}

/*
  Outputs the matrix to a stream. It is assumed that operator << is defined
  for C, and that C is "small" (in particular, has no newlines in its output.)
*/
template<typename C>
std::ostream& printMatrix(std::ostream& strm, const matrix::Matrix_base<C>& m,
			  unsigned long width)
{
  std::vector<unsigned long> widths(m.numColumns(),width);

  { std::ostringstream o;
    for (size_t j=0; j<m.numColumns(); ++j)
      for (size_t i=0; i<m.numRows(); ++i)
      {
        o.str(""); o << m(i,j);
	size_t w=o.str().length()+1;
        if (w>widths[j])
	  widths[j]=w;
      }
  }

  for (size_t i = 0; i < m.numRows(); ++i)
  {
    for (size_t j = 0; j < m.numColumns(); ++j)
      strm << std::setw(widths[j]) << m(i,j);

    strm << std::endl;
  }

  return strm;
}


/*
  Synopsis: prints out the monomial c.x^d.

  Preconditions: c is non-zero;

  Explanation: c and d are printed only if non-one, except in degree zero
  where c is always printed; x is the name of the indeterminate.

  The output format for the exponents is tex-like, but "q^13", not "q^{13}".
*/
template<typename C>
std::ostream& printMonomial(std::ostream& strm, C c, polynomials::Degree d,
			    const char* x)
{
  if (d == 0) // output c regardless
    strm << c;
  else
  {
    if (c<C(0) and c == C(-1)) // condition always false for unsigned types
      strm << '-';
    else if (c != C(1))
      strm << c;
    strm << x;
    if (d > 1)
      strm << "^" << d;
  }

  return strm;
}


/*
  Synopsis: outputs the polynomial  |p| on |strm|.

  It is assumed that operator<< exists for the coefficient type. The string
  |x| is the printed name of the indeterminate. Zero term are suppressed.

*/
template<typename C>
std::ostream& printPol(std::ostream& strm, const Polynomial<C>& p,
		       const char* x)
{
  std::ostringstream o; // accumulate in string for interpretation of width
  if (p.isZero())
    o << "0";
  else
    for (size_t i = p.size(); i-->0; )
      if (p[i]!=C(0)) // guaranteed true the first time
	printMonomial(i<p.degree() and p[i]>C(0) ? o<<'+' : o,p[i],i,x);

  return strm << o.str(); // now |strm.width()| is applied to whole polynomial
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

std::ostream& printTitsElt(std::ostream& strm, const TitsElt& a,
			   const TitsGroup& Tg)
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
std::ostream& printWeylElt(std::ostream& strm, const WeylElt& w,
			   const WeylGroup& W)
{
  strm << W.word(w);
  return strm;
}

/*
  Synopsis: outputs the list of WeylElts as words in the outer representation,
  with the given separator, prefix and postfix.
*/
std::ostream& printWeylList(std::ostream& strm, const WeylEltList& wl,
			    const WeylGroup& W, const char* sep,
			    const char* pre, const char* post)
{
  std::vector<WeylWord> wwl(wl.size());

  for (size_t i = 0; i < wl.size(); ++i)
    wwl[i]=W.word(wl[i]);

  basic_io::seqPrint(strm,wwl.begin(),wwl.end(),sep,pre,post);

  return strm;
}

// Instantiations

template std::ostream& prettyPrint
  (std::ostream&, const RankFlags&, size_t);

template std::ostream& prettyPrint
  (std::ostream&, const BitVector<constants::RANK_MAX>&);

template std::ostream& prettyPrint
  (std::ostream&,
   const std::vector<BitVector<constants::RANK_MAX> >&);

template std::ostream& printBasis
  (std::ostream&, const std::vector<Weight>&);

template std::ostream& printVector
  (std::ostream&, const std::vector<int>&, unsigned long);

template std::ostream& printMatrix
  (std::ostream&, const matrix::Matrix_base<int>&, unsigned long);

template std::ostream& printMatrix // in standardreprk coefficients are |long|
  (std::ostream&, const matrix::Matrix_base<long int>&, unsigned long);

template std::ostream& printMatrix // in |printBlockSizes|: |unsigned long|
  (std::ostream&, const matrix::Matrix_base<unsigned long>&, unsigned long);

template std::ostream& printMatrix
  (std::ostream&,
   const matrix::Matrix_base<Polynomial<int> >&,
   unsigned long);

template std::ostream& printMonomial
  (std::ostream&, int, polynomials::Degree, const char*);

template std::ostream& printPol
  (std::ostream&, const Polynomial<unsigned int>&, const char*);

template std::ostream& printPol
  (std::ostream&, const Polynomial<int>&, const char*);

} // namespace prettyprint

} // namespace atlas
