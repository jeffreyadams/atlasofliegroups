/*
  This is realform_io.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
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
#include "tori.h"

namespace atlas {

namespace {

  void pause() {;}

  std::ostream& printComplexType(std::ostream&, 
				 const lietype::SimpleLieType&);

  std::ostream& printSimpleType(std::ostream&, const gradings::Grading&,
				const lietype::SimpleLieType&, 
				lietype::TypeLetter);

  std::ostream& printType(std::ostream&, const gradings::Grading&,
			  const lietype::LieType&, 
			  const lietype::InnerClassType&);


class RealFormData {

  private:

    realform::RealForm d_realForm;
    gradings::Grading d_grading;
    rootdata::RootList d_orth;

  public:

// constructors and destructors
    RealFormData() {}

    RealFormData(realform::RealForm rf, const gradings::Grading& gr,
		 const rootdata::RootList& so)
      :d_realForm(rf),d_grading(gr),d_orth(so) {}

    ~RealFormData() {}

// accessors
    const gradings::Grading& grading() const {
      return d_grading;
    }

    realform::RealForm realForm() const {
      return d_realForm;
    }

    const rootdata::RootList& orth() const {
      return d_orth;
    }
  };

  bool operator< (const RealFormData& first, const RealFormData& second);

}

/*****************************************************************************

        Chapter I -- The Interface class

  ... explain here when it is stable ...

******************************************************************************/

namespace realform_io {

Interface::Interface(const complexredgp::ComplexReductiveGroup& G,
		     const layout::Layout& lo)

/*
  Synopsis: builds the standard real form interface for G and lo.

  Explanation: d_in[rf] is the inner numbering (the one used by G) of the
  weak real form rf; d_out[rf] is the outer numbering (the one used by the
  interface) of rf --- so d_in and d_out are inverses of each other. The
  names are names for the corresponding Lie algebras.
*/

{
  using namespace cartanclass;
  using namespace gradings;
  using namespace realform;
  using namespace rootdata;

  const RootDatum& rd = G.rootDatum();
  const Fiber& fundf = G.fundamental();

  std::vector<RealFormData> rfd(G.numRealForms());

  for (RealForm rf = 0; rf < G.numRealForms(); ++rf) {
    RootList so;
    toMostSplit(so,fundf,rf,rd);
    Grading gr;
    specialGrading(gr,fundf,rf,rd);
    RealFormData data(rf,gr,so);
    rfd[rf] = data;
  }

  std::sort(rfd.begin(),rfd.end());

  // write real form correspondence
  d_in.resize(rfd.size());
  d_out.resize(rfd.size());

  for (size_t j = 0; j < rfd.size(); ++j) {
    d_in[j] = rfd[j].realForm();
    d_out[d_in[j]] = j;
  }

  // write names
  d_name.resize(rfd.size());
  std::ostringstream os;

  for (size_t j = 0; j < rfd.size(); ++j) {
    os.str("");
    printType(os,rfd[j].grading(),lo.d_type,lo.d_inner);
    d_name[j] = os.str();
  }
}

Interface::Interface(const complexredgp::ComplexReductiveGroup& G,
		     const layout::Layout& lo, tags::DualTag)

/*
  Synopsis: like the previous one, but for the _dual_ real forms.
*/

{
  using namespace cartanclass;
  using namespace gradings;
  using namespace lietype;
  using namespace realform;
  using namespace rootdata;
  using namespace tags;

  const RootDatum& rd = RootDatum(G.rootDatum(),DualTag());
  const Fiber& fundf = G.dualFundamental();

  std::vector<RealFormData> rfd(G.numDualRealForms());

  for (RealForm rf = 0; rf < G.numDualRealForms(); ++rf) {
    RootList so;
    toMostSplit(so,fundf,rf,rd);
    Grading gr;
    specialGrading(gr,fundf,rf,rd);
    RealFormData data(rf,gr,so);
    rfd[rf] = data;
  }

  std::sort(rfd.begin(),rfd.end());

  // write real form correspondence
  d_in.resize(rfd.size());
  d_out.resize(rfd.size());

  for (size_t j = 0; j < rfd.size(); ++j) {
    d_in[j] = rfd[j].realForm();
    d_out[d_in[j]] = j;
  }

  // write names
  d_name.resize(rfd.size());
  std::ostringstream os;

  pause();

  LieType dlt;
  dualLieType(dlt,lo.d_type);
  InnerClassType dict;
  dualInnerClassType(dict,lo.d_inner,lo.d_type);

  for (size_t j = 0; j < rfd.size(); ++j) {
    os.str("");
    printType(os,rfd[j].grading(),dlt,dict);
    d_name[j] = os.str();
  }
}

/******* copy, assignment and swap *******************************************/
void Interface::swap(Interface& other)

{  
  d_in.swap(other.d_in);
  d_out.swap(other.d_out);
  d_name.swap(other.d_name);

  return;
}

/******* accessors ***********************************************************/
const char* Interface::typeName(realform::RealForm rf) const

/*
  Synopsis: returns the name of the real form.

  Precondition: rf is an outer real form number.
*/

{
  return d_name[rf].c_str();
}

}

/*****************************************************************************

        Chapter II -- Functions declared in realform_io.h

  ... explain here when it is stable ...

******************************************************************************/

namespace realform_io {

std::ostream& printRealForms(std::ostream& strm, const Interface& I)

/*
  Synopsis: outputs the list of real forms in I.
*/

{  
  for (size_t j = 0; j < I.numRealForms(); ++j) {
    std::cout << j << ": " << I.typeName(j) << std::endl;
  }

  return strm;
}

}

/*****************************************************************************

        Chapter III -- Auxiliary functions local to this module

  ... explain here when it is stable ...

******************************************************************************/

namespace {

bool operator< (const RealFormData& first, const RealFormData& second)

/*
  Synopsis: comparison of RealFormData.

  Sorts the data by compacity (most compact forms first), and for equal
  compacity, by their grading. This ensures that the ordering will be the
  same for all covering groups.
*/

{
  using namespace gradings;

  size_t firstSize = first.orth().size();
  size_t secondSize = second.orth().size();

  if (firstSize < secondSize)
    return true;
  if (firstSize > secondSize)
    return false;

  const Grading& firstGrading = first.grading();
  const Grading& secondGrading = second.grading();

  return firstGrading < secondGrading;
}

std::ostream& printComplexType(std::ostream& strm, 
			       const lietype::SimpleLieType& slt)

/*
  Synopsis: prints out the complex form of slt.
*/

{
  using namespace lietype;

  switch (type(slt)) {
  case 'A':
    strm << "sl(";
    strm << rank(slt)+1 << ",C)";
    break;
  case 'B':
    strm << "so(";
    strm << 2*rank(slt)+1 << ",C)";
    break;
  case 'C':
    strm << "sp(";
    strm << 2*rank(slt) << ",C)";
    break;
  case 'D':
    strm << "so(";
    strm << 2*rank(slt) << ",C)";
    break;
  case 'E':
    strm << "e";
    strm << rank(slt) << "(C)";
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
    if (rank(slt) > 1) {
      strm << "^";
      strm << rank(slt);
    }
    break;
  default:
    assert(false && "unexpected type in printComplexType");
    break;
  }

  return strm;
}

std::ostream& printSimpleType(std::ostream& strm, const gradings::Grading& gr,
			      const lietype::SimpleLieType& slt,
			      const lietype::TypeLetter ic)

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

{
  using namespace lietype;

  size_t fb = gr.firstBit();

  switch (type(slt)) {
  case 'A':
    if (rank(slt) == 1) {
      if (gr.count() == 1)
	strm << "sl(2,R)";
      else
	strm << "su(2)";
    }
    else
      switch (ic) {
      case 's':
	if (rank(slt) & 1UL) {
	  if (gr.count() == 1) {
	    strm << "sl(";
	    strm << rank(slt)+1 << ",R)";
	  } else {
	    strm << "sl(";
	    strm << (rank(slt)+1)/2 << ",H)";
	  }
	} else {
	    strm << "sl(";
	    strm << rank(slt)+1 << ",R)";
	}
	break;
      case 'c':
	if (gr.count()) {
	  size_t p = rank(slt) - fb;
	  strm << "su(";
	  strm << p << ",";
	  strm << fb+1 << ")";
	} else {
	  strm << "su(";
	  strm << rank(slt)+1 << ")";
	}
	break;
      }
    break;
  case 'B':
    if (gr.count()) {
      size_t mid = (rank(slt)-1)/2;
      if (fb < mid) {
	strm << "so(";
	strm << 2*rank(slt)+1-2*(fb+1) << "," << 2*(fb+1) << ")";
      } else {
	size_t c = rank(slt)-fb-1;
	strm << "so(";
	strm << 2*rank(slt)-2*c << "," <<  2*c+1 << ")";
      }
    } else {
      strm << "so(";
      strm << 2*rank(slt)+1 << ")";
    }
    break;
  case 'C':
    if (gr.count() == 0) {
      strm << "sp(";
      strm << rank(slt) << ")";
    } else if (fb == rank(slt)-1) {
      strm << "sp(";
      strm << 2*rank(slt) << ",R)";
    } else {
      strm << "sp(";
      strm << rank(slt)-fb-1 << "," << fb+1 << ")";
    }
    break;
  case 'D':
    if (rank(slt) & 1UL)
      switch (ic) {
      case 's':
	if (gr.count() == 0) { // distinguished form
	  strm << "so(";
	  strm << 2*rank(slt)-1<< "," << 1 << ")";
	} else {
	  strm << "so(";
	  strm << 2*rank(slt)-2*fb-3 << "," << 2*fb+3 << ")";
	}
	break;
      case 'c':
	if (gr.count() == 0) {
	  strm << "so(";
	  strm << 2*rank(slt) << ")";
	} else if (fb == rank(slt)-2) {
	  strm << "so*(";
	  strm << 2*rank(slt) << ")";
	} else {
	  strm << "so(";
	  strm << 2*rank(slt)-2*fb-2 << "," << 2*fb+2 << ")";
	}
	break;
      case 'u':
        if (gr.count() == 0) { // distinguished form
          strm << "so(";
          strm << 2*rank(slt)-1<< "," << 1 << ")";
        } else {
          strm << "so(";
          strm << 2*rank(slt)-2*fb-3 << "," << 2*fb+3 << ")";
        }
      }
    else
      switch (ic) {
      case 's':
      case 'c':
	if (gr.count() == 0) { // compact type
	  strm << "so(";
	  strm << 2*rank(slt) << ")";
	} else if (fb == rank(slt)-2) { // so* type; labelling depends on rank
	  strm << "so*(";
	  if ((rank(slt) >> 1) & 1UL) // rank not divisible by 4
	    strm << 2*rank(slt) << ")[1,0]";
	  else // rank divisible by 4
	    strm << 2*rank(slt) << ")[0,1]";
	} else if (fb == rank(slt)-1) { // so* type; labelling depends on rank
	  strm << "so*(";
	  if ((rank(slt) >> 1) & 1UL) // rank not divisible by 4
	    strm << 2*rank(slt) << ")[0,1]";
	  else // rank divisible by 4
	    strm << 2*rank(slt) << ")[1,0]";
	} else { // so(p,q) type
	  strm << "so(";
	  strm << 2*rank(slt)-2*fb-2 << "," << 2*fb+2 << ")";
	}
	break;
      case 'u':
	if (gr.count() == 0) { // distinguished form
	  strm << "so(";
	  strm << 2*rank(slt)-1<< "," << 1 << ")";
	} else {
	  strm << "so(";
	  strm << 2*rank(slt)-2*fb-3 << "," << 2*fb+3 << ")";
	}
	break;
      }
    break;
  case 'E':
    switch (rank(slt)) {
    case 6:
      switch (ic) {
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
    if (rank(slt) > 1) {
      strm << "^";
      strm << rank(slt);
    }
    break;
  default: // cannot happen
    assert(false && "unexpected type in printSimpleType");
    break;
  }

  return strm;
}

std::ostream& printType(std::ostream& strm, 
			const gradings::Grading& d_gr,
			const lietype::LieType& lt,
			const lietype::InnerClassType& ict)

/*
  Synopsis: outputs the real form of lt represented by gr.

  Precondition: gr contains a grading of the simple roots of the root system
  that are imaginary for the fundamental involution, which represents the
  given inner class, and which contains one noncompact root for each
  noncompact noncomplex simple factor.
*/

{
  using namespace gradings;
  using namespace lietype;

  Grading gr = d_gr;
  size_t c = 0;

  for (size_t j = 0; j < ict.size(); ++j) {
    SimpleLieType slt = lt[c];
    if (ict[j] == 'C') {
      printComplexType(strm,slt);
      gr >>= semisimpleRank(slt);
      ++c;
    } else {
      Grading grs = gr;
      grs.truncate(rank(slt));
      printSimpleType(strm,grs,slt,ict[j]);
    }
    gr >>= semisimpleRank(slt);
    ++c;
    if (j < ict.size()-1)
      strm << ".";
  }

  return strm;
}

}

}

