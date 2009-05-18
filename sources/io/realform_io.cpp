/*
  This is realform_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "realform_io.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

#include "cartanclass.h"
#include "complexredgp.h"
#include "complexredgp_io.h"
#include "gradings.h"
#include "latticetypes.h"
#include "layout.h"
#include "lietype.h"
#include "rootdata.h"
#include "setutils.h"
#include "tori.h"

namespace atlas {

namespace {

  std::ostream& printComplexType(std::ostream&,
				 const lietype::SimpleLieType&);

  std::ostream& printSimpleType(std::ostream&, const gradings::Grading&,
				const lietype::SimpleLieType&,
				lietype::TypeLetter);

  std::ostream& printType(std::ostream& out,
			  const gradings::Grading& gr,
			  const lietype::LieType& lt,
			  const lietype::InnerClassType& ict,
			  const setutils::Permutation& perm);

// a class to temporarily group data necessary for classifying real forms
class RealFormData
{
  realform::RealForm d_realForm;
  gradings::Grading d_grading;
  rootdata::RootSet d_orth;

public:

// constructors and destructors
  RealFormData() {}

  RealFormData(realform::RealForm rf,
	       const gradings::Grading& gr,
	       const rootdata::RootSet& so)
    :d_realForm(rf),d_grading(gr),d_orth(so) {}

// accessors
  const gradings::Grading& grading() const { return d_grading; }
  realform::RealForm realForm() const { return d_realForm; }
  const rootdata::RootSet& orth() const { return d_orth; }
};

// this determines the external numbering of real forms
bool operator< (const RealFormData& first, const RealFormData& second);

}

/*****************************************************************************

        Chapter I -- The Interface class

******************************************************************************/

namespace realform_io {


/*
  Synopsis: builds the standard real form interface for G and lo.

  Explanation: d_in[rf] is the inner numbering (the one used by G) of the
  weak real form rf; d_out[rf] is the outer numbering (the one used by the
  interface) of rf --- so d_in and d_out are inverses of each other. The
  names are names for the corresponding Lie algebras.
*/
Interface::Interface(const complexredgp::ComplexReductiveGroup& G,
		     const layout::Layout& lo)
: d_in(G.numRealForms()), d_out(G.numRealForms()), d_name(G.numRealForms())
{
  const size_t nrf = G.numRealForms();
  const rootdata::RootSystem& rs = G.rootSystem();
  const cartanclass::Fiber& fundf = G.fundamental();

  std::vector<RealFormData> rfd; rfd.reserve(nrf);

  for (realform::RealForm rf = 0; rf<nrf; ++rf)
  {
    rootdata::RootSet so = cartanclass::toMostSplit(fundf,rf,rs);
    gradings::Grading gr = cartanclass::specialGrading(fundf,rf,rs);
    rfd.push_back(RealFormData(rf,gr,so));
  }

  std::sort(rfd.begin(),rfd.end());

  for (size_t j = 0; j<nrf; ++j)
  {
    d_in[j] = rfd[j].realForm();
    d_out[d_in[j]] = j;
  }

  // write names
  std::ostringstream os;

  for (size_t j = 0; j<nrf; ++j)
  {
    os.str("");
    printType(os,rfd[j].grading(),lo.d_type,lo.d_inner,lo.d_perm);
    d_name[j] = os.str();
  }
}


/*
  Synopsis: like the previous one, but for the _dual_ real forms.
*/
Interface::Interface(const complexredgp::ComplexReductiveGroup& G,
		     const layout::Layout& lo, tags::DualTag)
: d_in(G.numDualRealForms())
, d_out(G.numDualRealForms())
, d_name(G.numDualRealForms())
{
  const size_t ndrf = G.numDualRealForms();
  const rootdata::RootSystem& rs = G.dualRootSystem();
  const cartanclass::Fiber& fundf = G.dualFundamental();

  std::vector<RealFormData> rfd; rfd.reserve(ndrf);

  for (realform::RealForm rf = 0; rf<ndrf; ++rf)
  {
    rootdata::RootSet so = cartanclass::toMostSplit(fundf,rf,rs);
    gradings::Grading gr = cartanclass::specialGrading(fundf,rf,rs);
    rfd.push_back(RealFormData(rf,gr,so));
  }

  std::sort(rfd.begin(),rfd.end());

  for (size_t j=0; j<ndrf; ++j)
  {
    d_in[j] = rfd[j].realForm();
    d_out[d_in[j]] = j;
  }

  // write names
  std::ostringstream os;

  lietype::LieType dual_lt = lietype::dual_type(lo.d_type);
  lietype::InnerClassType dual_ict = lietype::dual_type(lo.d_inner,lo.d_type);

  for (size_t j=0; j<ndrf; ++j)
  {
    os.str("");
    printType(os,rfd[j].grading(),dual_lt,dual_ict,lo.d_perm);
    d_name[j] = os.str();
  }
}

/******* copy, assignment and swap *******************************************/
void Interface::swap(Interface& other)

{
  d_in.swap(other.d_in);
  d_out.swap(other.d_out);
  d_name.swap(other.d_name);
}

/******* accessors ***********************************************************/

/*
  Synopsis: returns the name of the real form.

  Precondition: rf is an outer real form number.
*/

const char* Interface::typeName(realform::RealForm rf) const
{
  return d_name[rf].c_str();
}

}

/*****************************************************************************

        Chapter II -- Functions declared in realform_io.h

******************************************************************************/

namespace realform_io {

/*
  Synopsis: outputs the list of real forms in I.
*/

std::ostream& printRealForms(std::ostream& strm, const Interface& I)
{
  for (size_t j = 0; j < I.numRealForms(); ++j) {
    std::cout << j << ": " << I.typeName(j) << std::endl;
  }

  return strm;
}

}

/*****************************************************************************

        Chapter III -- Auxiliary functions local to this module


******************************************************************************/

namespace {

/*
  Synopsis: comparison of RealFormData.

  Sorts the data by compacity (most compact forms first), and for equal
  compacity, by their grading. This ensures that the ordering will be the same
  for all covering groups. However the comparisons of gradings depends on the
  permutation of simple roots relative to the standard ordering of the diagram
  for the Lie type, which could be improved.

*/
bool operator< (const RealFormData& first, const RealFormData& second)
{
  size_t firstSize = first.orth().size();
  size_t secondSize = second.orth().size();

  if (firstSize != secondSize)
    return firstSize < secondSize;

  return first.grading() < second.grading();
}


/*
  Synopsis: prints out the complex form of slt.
*/
std::ostream& printComplexType(std::ostream& strm,
			       const lietype::SimpleLieType& slt)
{
  size_t rk = slt.rank();
  switch (slt.type())
  {
  case 'A':
    strm << "sl(";
    strm << rk+1 << ",C)";
    break;
  case 'B':
    strm << "so(";
    strm << 2*rk+1 << ",C)";
    break;
  case 'C':
    strm << "sp(";
    strm << 2*rk << ",C)";
    break;
  case 'D':
    strm << "so(";
    strm << 2*rk << ",C)";
    break;
  case 'E':
    strm << "e";
    strm << rk << "(C)";
    break;
  case 'F':
  case 'f':
    strm << "f4(C)";
    break;
  case 'G':
  case 'g':
    strm << "g2(C)";
    break;
  case 'T':
    strm << "gl(1,C)";
    if (rk > 1) {
      strm << "^";
      strm << rk;
    }
    break;
  default:
    assert(false && "unexpected type in printComplexType");
    break;
  }

  return strm;
}


/*
  Synopsis: prints out the real form of slt represented by gr.

  Algorithm: it turns out that for each irreducible Dynkin diagram,
  and for each non-compact equal rank real form, there is always at least
  one grading with exactly one noncompact simple imaginary root. In that case,
  gr is the one for which the noncompact root is smallest. In the non-equal
  rank case, this is true for all ireducible types except A_n, where the
  imaginary root system is not irreducible. In that case, of course, the
  only choice is whether gr is trivial or not, so it is easy enough to decide.
*/
std::ostream& printSimpleType(std::ostream& strm, const gradings::Grading& gr,
			      const lietype::SimpleLieType& slt,
			      const lietype::TypeLetter ic)
{
  size_t fb = gr.firstBit();
  size_t rk = slt.rank();

  switch (slt.type())
  {
  case 'A':
    if (rk == 1)
      if (gr.count() == 1)
	strm << "sl(2,R)";
      else
	strm << "su(2)";
    else
      switch (ic)
      {
      case 's':
	if (rk & 1UL) {
	  if (gr.count() == 1) {
	    strm << "sl(";
	    strm << rk+1 << ",R)";
	  } else {
	    strm << "sl(";
	    strm << (rk+1)/2 << ",H)";
	  }
	} else {
	    strm << "sl(";
	    strm << rk+1 << ",R)";
	}
	break;
      case 'c':
	if (gr.count()>0) {
	  size_t p = rk - fb;
	  strm << "su(";
	  strm << p << ",";
	  strm << fb+1 << ")";
	} else {
	  strm << "su(";
	  strm << rk+1 << ")";
	}
	break;
      }
    break;
  case 'B':
    if (gr.count()>0)
    {
      size_t mid = (rk-1)/2;
      if (fb <= mid) {
	strm << "so(";
	strm << 2*rk+1-2*(fb+1) << "," << 2*(fb+1) << ")";
      } else {
	size_t c = rk-fb-1;
	strm << "so(";
	strm << 2*rk-2*c << "," <<  2*c+1 << ")";
      }
    }
    else
    {
      strm << "so(";
      strm << 2*rk+1 << ")";
    }
    break;
  case 'C':
    if (gr.count() == 0) {
      strm << "sp(";
      strm << rk << ")";
    } else if (fb == rk-1) {
      strm << "sp(";
      strm << 2*rk << ",R)";
    } else {
      strm << "sp(";
      strm << rk-fb-1 << "," << fb+1 << ")";
    }
    break;
  case 'D':
    if ((rk & 1UL)!=0)
      switch (ic)
      {
      case 's':
	if (gr.count() == 0) { // distinguished form
	  strm << "so(";
	  strm << 2*rk-1<< "," << 1 << ")";
	} else {
	  strm << "so(";
	  strm << 2*rk-2*fb-3 << "," << 2*fb+3 << ")";
	}
	break;
      case 'c':
	if (gr.count() == 0) {
	  strm << "so(";
	  strm << 2*rk << ")";
	} else if (fb == rk-2) {
	  strm << "so*(";
	  strm << 2*rk << ")";
	} else {
	  strm << "so(";
	  strm << 2*rk-2*fb-2 << "," << 2*fb+2 << ")";
	}
	break;
      case 'u':
        if (gr.count() == 0) { // distinguished form
          strm << "so(";
          strm << 2*rk-1<< "," << 1 << ")";
        } else {
          strm << "so(";
          strm << 2*rk-2*fb-3 << "," << 2*fb+3 << ")";
        }
      }
    else
      switch (ic) {
      case 's':
      case 'c':
	if (gr.count() == 0) { // compact type
	  strm << "so(";
	  strm << 2*rk << ")";
	} else if (fb == rk-2) { // so* type; labelling depends on rank
	  strm << "so*(";
	  if ((rk >> 1) & 1UL) // rank not divisible by 4
	    strm << 2*rk << ")[1,0]";
	  else // rank divisible by 4
	    strm << 2*rk << ")[0,1]";
	} else if (fb == rk-1) { // so* type; labelling depends on rank
	  strm << "so*(";
	  if ((rk >> 1) & 1UL) // rank not divisible by 4
	    strm << 2*rk << ")[0,1]";
	  else // rank divisible by 4
	    strm << 2*rk << ")[1,0]";
	} else { // so(p,q) type
	  strm << "so(";
	  strm << 2*rk-2*fb-2 << "," << 2*fb+2 << ")";
	}
	break;
      case 'u':
	if (gr.count() == 0) { // distinguished form
	  strm << "so(";
	  strm << 2*rk-1<< "," << 1 << ")";
	} else {
	  strm << "so(";
	  strm << 2*rk-2*fb-3 << "," << 2*fb+3 << ")";
	}
	break;
      }
    break;
  case 'E':
    switch (rk)
    {
    case 6:
      switch (ic)
      {
      case 's':
	if (gr.count() == 0)
	  strm << "e6(f4)";
	else
	  strm << "e6(R)";
	break;
      case 'c':
	if (gr.count() == 0)
	  strm << "e6";
	else if (fb == 0)
	  strm << "e6(so(10).u(1))";
	else if (fb == 1)
	  strm << "e6(su(6).su(2))";
	break;
      }
      break;
    case 7:
      if (gr.count() == 0)
	strm << "e7";
      else if (fb == 6)
	strm << "e7(e6.u(1))";
      else if (fb == 0)
	strm << "e7(so(12).su(2))";
      else if (fb == 1)
	strm << "e7(R)";
      break;
    case 8:
      if (gr.count() == 0)
	strm << "e8";
      else if (fb == 2)
	strm << "e8(e7.su(2))";
      else if (fb == 0)
	strm << "e8(R)";
      break;
    }
    break;
  case 'F':
    if (gr.count() == 0)
      strm << "f4";
    else if (fb == 2)
      strm << "f4(so(9))";
    else if (fb == 0)
      strm << "f4(R)";
    break;
  case 'f': // the dual of F; short roots are first here
    if (gr.count() == 0)
      strm << "f4";
    else if (fb == 2)
      strm << "f4(R)";
    else if (fb == 0)
      strm << "f4(so(9))";
    break;
  case 'G':
  case 'g':
    if (gr.count() == 0)
      strm << "g2";
    else
      strm << "g2(R)";
    break;
  case 'T':
    if (ic == 's')
      strm << "gl(1,R)";
    else if (ic == 'c')
      strm << "u(1)";
    if (rk > 1) {
      strm << "^";
      strm << rk;
    }
    break;
  default: // cannot happen
    assert(false && "unexpected type in printSimpleType");
    break;
  }

  return strm;
}


/*
  Synopsis: outputs the real form of lt represented by gr.

  Precondition: |gr| contains a grading of the simple-simaginary roots (simple
  roots of the imaginary roots system for the fundamental involution), which
  represents the given real form, and which contains one noncompact root for
  each noncompact noncomplex simple factor.
*/
std::ostream& printType(std::ostream& strm,
			const gradings::Grading& d_gr,
			const lietype::LieType& lt,
			const lietype::InnerClassType& ict,
			const setutils::Permutation& perm)
{
  gradings::Grading gr = d_gr;
  gr.permute(perm); // adapt grading to standard ordering of type |lt|

  size_t c = 0;

  for (size_t j = 0; j < ict.size(); ++j)
  {
    lietype::SimpleLieType slt = lt[c];
    if (ict[j] == 'C')
    {
      printComplexType(strm,slt);
      gr >>= slt.semisimple_rank();
      ++c;
    }
    else
    {
      gradings::Grading grs = gr;
      grs.truncate(slt.rank());
      printSimpleType(strm,grs,slt,ict[j]);
    }
    gr >>= slt.semisimple_rank();
    ++c;
    if (j < ict.size()-1)
      strm << ".";
  }

  return strm;
}

}

}
