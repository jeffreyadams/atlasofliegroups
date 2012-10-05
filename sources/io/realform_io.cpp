/*
  This is realform_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "realform_io.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

#include "gradings.h"	// containment
#include "complexredgp.h"// methods
#include "lietype.h"	// methods

namespace atlas {

namespace {

  std::ostream& printComplexType(std::ostream&,
				 const SimpleLieType&);

  std::ostream& printSimpleType(std::ostream&, const Grading&,
				const SimpleLieType&,
				lietype::TypeLetter);

  std::ostream& printType(std::ostream& out,
			  const Grading& gr,
			  const lietype::Layout& lo);

// a class to temporarily group data necessary for classifying real forms
class RealFormData
{
  RealFormNbr d_realForm;
  Grading d_grading;
  RootNbrSet d_orth;

public:

// constructors and destructors
  RealFormData() {}

  RealFormData(RealFormNbr rf,
	       const Grading& gr,
	       const RootNbrSet& so)
    :d_realForm(rf),d_grading(gr),d_orth(so) {}

// accessors
  const Grading& grading() const { return d_grading; }
  RealFormNbr realForm() const { return d_realForm; }
  const RootNbrSet& orth() const { return d_orth; }
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
Interface::Interface(const ComplexReductiveGroup& G,
		     const lietype::Layout& lo)
: d_in(G.numRealForms()), d_out(G.numRealForms()), d_name(G.numRealForms())
{
  const size_t nrf = G.numRealForms();
  const RootSystem& rs = G.rootSystem();
  const Fiber& fundf = G.fundamental();

  std::vector<RealFormData> rf_data; rf_data.reserve(nrf);

  for (RealFormNbr rf = 0; rf<nrf; ++rf)
  {
    RootNbrSet so = cartanclass::toMostSplit(fundf,rf,rs);
    Grading gr = cartanclass::specialGrading(fundf,rf,rs);
    rf_data.push_back(RealFormData(rf,gr,so));
  }

  std::sort(rf_data.begin(),rf_data.end());

  for (size_t i = 0; i<nrf; ++i)
  {
    d_in[i] = rf_data[i].realForm();
    d_out[d_in[i]] = i;
  }

  // write names
  std::ostringstream os;

  for (size_t i = 0; i<nrf; ++i)
  {
    os.str("");
    printType(os,rf_data[i].grading(),lo);
    d_name[i] = os.str();
  }
}


/*
  Synopsis: like the previous one, but for the _dual_ real forms.
*/
Interface::Interface(const ComplexReductiveGroup& G,
		     const lietype::Layout& lo, tags::DualTag)
: d_in(G.numDualRealForms())
, d_out(G.numDualRealForms())
, d_name(G.numDualRealForms())
{
  const size_t ndrf = G.numDualRealForms();
  const RootSystem& drs = G.dualRootSystem();
  const Fiber& dfundf = G.dualFundamental();
  const lietype::Layout dlo = dual(lo);

  std::vector<RealFormData> rf_data; rf_data.reserve(ndrf);

  for (RealFormNbr drf = 0; drf<ndrf; ++drf)
  {
    RootNbrSet so = cartanclass::toMostSplit(dfundf,drf,drs);
    Grading gr = cartanclass::specialGrading(dfundf,drf,drs);
    rf_data.push_back(RealFormData(drf,gr,so));
  }

  std::sort(rf_data.begin(),rf_data.end());

  for (size_t i=0; i<ndrf; ++i)
  {
    d_in[i] = rf_data[i].realForm();
    d_out[d_in[i]] = i;
  }

  // write names
  std::ostringstream os;
  for (size_t i=0; i<ndrf; ++i)
  {
    os.str("");
    printType(os,rf_data[i].grading(),dlo);
    d_name[i] = os.str();
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

const char* Interface::typeName(RealFormNbr rf) const
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
  for (size_t i = 0; i < I.numRealForms(); ++i) {
    std::cout << i << ": " << I.typeName(i) << std::endl;
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
			       const SimpleLieType& slt)
{
  size_t rk = slt.rank();
  switch (slt.type())
  {
  case 'A':
    strm << "sl(" << rk+1 << ",C)";
    break;
  case 'B':
    strm << "so(" << 2*rk+1 << ",C)";
    break;
  case 'C':
    strm << "sp(" << 2*rk << ",C)";
    break;
  case 'D':
    strm << "so(" << 2*rk << ",C)";
    break;
  case 'E':
    strm << "e" << rk << "(C)";
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
    if (rk > 1) // this can no longer occur
    {
      strm << "^" << rk;
    }
    break;
  default:
    assert(false && "unexpected type in printComplexType");
    break;
  }

  return strm;
}

// we often need to print either '(n)' or '(n-m,m)' in printSimpleType
inline std::ostream& split(std::ostream& s,size_t n,size_t m)
{
  s << '(';
  if (m==0)
    s << n;
  else
  {
    if (2*m>n)
      m=n-m; // make sure larger number is listed first
    s << n-m << ',' << m;
  }
  return s << ')';
}

/*
  Synopsis: prints out the real form of slt represented by gr.

  Algorithm: it turns out that for each irreducible Dynkin diagram, and for
  each non-compact equal rank real form, there is always at least one grading
  with exactly one noncompact simple-imaginary root. In that case, this
  function is called with such a grading; we cannot rely on the noncompact
  root being the minimal one possible, as simple roots may have been in an
  unusual order at the point of selecting |gr|, even though such a permutation
  has been straightend out before calling this function. In the non-equal rank
  case, such a grading still exists for all irreducible types except A_n, when
  the imaginary root system is reducible. In that case however the only choice
  is whether |gr| is trivial (sl(n,H)) or not, which is easy enough to decide.
*/
std::ostream& printSimpleType(std::ostream& strm, const Grading& gr,
			      const SimpleLieType& slt,
			      const lietype::TypeLetter ic)
{
  size_t rk = slt.rank();
  bool gr_trivial = gr.count()==0;
  size_t m = gr_trivial ? 0 : gr.firstBit()+1; // usually relevant for printing

  switch (slt.type())
  {
  case 'A':
    if (rk == 1)
      strm << (gr_trivial ? "su(2)" : "sl(2,R)");
    else
    {
      size_t n=rk+1;
      if (ic=='c')
	split(strm << "su",n,m);
      else if (rk%2!=0 and gr_trivial) // both parts of this condition needed
	strm << "sl(" << n/2 << ",H)";
      else
	strm << "sl(" << n << ",R)";
    }
    break;
  case 'B':
    split(strm << "so",2*rk+1,2*m);
    break;
  case 'C':
    if (m == rk)
      strm << "sp(" << 2*rk << ",R)";
    else
      split (strm << "sp",rk,m);
    break;
  case 'D':
    {
      size_t n=2*rk;
      if (ic=='c' or (rk%2==0 and ic=='s')) // equal rank case
	if (m<rk-1) // so(2p,2q) form
	  split (strm << "so",n,2*m);
	else if (rk%2!=0)
	  strm << "so*(" << n << ")";
	else // so* type with label depending on m and parity of |rk/2|
	  strm << "so*(" << n << ((rk%4==0)==(m==rk) ? ")[1,0]" : ")[0,1]");
      else // unequal rank case
	split(strm << "so",n,2*m+1);
    }
    break;
  case 'E':
    strm << 'e' << rk; // it always starts like this
    if (gr_trivial and (ic=='c' or rk>6)) // for |rk>6|, |'s'| means |'c'|
      break; // from |switch|: compact reak form
    strm << '(';
    if (rk==6) // E6, noncompact forms
      if (ic=='c')
	strm << (m==1 or m==6 ? "so(10).u(1)" : "su(6).su(2)");
      else // unequal rank
	strm << (gr_trivial ? "f4" : "R");
    else if (rk==7) // E7, noncompact forms
      strm << (m==7 ? "e6.u(1)" : m==2 or m==5 ? "R" : "so(12).su(2)");
    else // E8, noncompact forms
      strm <<(gr.any(RankFlags(0xCC)) ? "e7.su(2)" : "R");
    strm << ')';
    break;
  case 'f':
    m=5-m; // reverse interpretation of |m| for type 'f'
  case 'F':
    strm << (gr_trivial ? "f4" : m>=3 ? "f4(so(9))" : "f4(R)");
   break;
  case 'g': // for type G root order is unused
  case 'G':
    strm << (gr_trivial ? "g2" : "g2(R)");
    break;
  case 'T':
    strm << (ic=='c' ? "u(1)" : "gl(1,R)");
    if (rk > 1)
      strm << "^" << rk;
    break;
  default: // cannot happen
    assert(false && "unexpected type in printSimpleType");
    break;
  }

  return strm;
}


/*
  Synopsis: outputs the real form of lt represented by gr.

  Precondition: |gr| contains a grading of the simple-imaginary roots (simple
  roots of the imaginary roots system for the fundamental involution), which
  represents the given real form, and which contains one noncompact root for
  each noncompact noncomplex simple factor.
*/
std::ostream& printType(std::ostream& strm,
			const Grading& d_gr,
			const lietype::Layout& lo)
{
  const LieType& lt=lo.d_type;
  const InnerClassType& ict=lo.d_inner;

  // adapt grading to standard ordering of type |lt|
  Grading gr = lo.d_perm.pull_back(d_gr);

  size_t i = 0; // index of simple factor in |lt|

  for (size_t j = 0; j<ict.size(); gr>>=lt[i].semisimple_rank(),++i,++j)
  {
    SimpleLieType slt = lt[i];
    if (ict[j] == 'C')
    {
      printComplexType(strm,slt);
      gr >>= slt.semisimple_rank();  ++i; // skip one more simple Lie type
    }
    else
    {
      Grading grs = gr;
      grs.truncate(slt.rank());
      printSimpleType(strm,grs,slt,ict[j]);
    }
    if (j < ict.size()-1)
      strm << ".";
  }

  return strm;
}

} // |namespace|

} // |namespace atlas|
