/*
  This is basic_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "basic_io.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "matrix.h"  // vectors, matrices
#include "ratvec.h" // rational vectors

#include "partition.h" // instantions

#include "bitvector.h"
#include "lietype.h" // Lie types

#include "prettyprint.h" // |printPol| (why isn't this needed?)

#include "repr.h" // printing of parameters and parameter polynomials

/*****************************************************************************

        Chapter I -- operator<< functions

******************************************************************************/

namespace atlas {

namespace basic_io {

// |seqPrint| was here once, but was moved back to basic_io_def.h

  } // |namespace basic_io|

namespace bitset {

template<size_t d>
  std::ostream& operator<< (std::ostream& strm, const BitSet<d>& b)
{
  for (size_t i = 0; i < d; ++i)
    strm << (b[i]?'1':'0');
  return strm;
}

} // |namespace bitset|

namespace polynomials {

template <typename C>
  std::ostream& operator<< (std::ostream& strm, const Polynomial<C>& P)
{
  return prettyprint::printPol(strm,P,"q");
}

} // |namespace polynomials|

namespace bitvector {

template<size_t dim>
  std::ostream& operator<< (std::ostream& strm,  const BitVector<dim>& v)
{
  for (size_t i=0; i<v.size(); ++i)
    strm << (v[i]?'1':'0');
  return strm;
}

} // |namespace bitvector|


namespace lietype {

std::ostream& operator<< (std::ostream& strm, const SimpleLieType& slt)
{
  return strm << slt.type() << slt.rank();
}

std::ostream& operator<< (std::ostream& strm, const LieType& lt)
{
  return basic_io::seqPrint(strm,lt.begin(),lt.end(),".");
}

std::ostream& operator<< (std::ostream& strm, const InnerClassType& ict)
{
  return basic_io::seqPrint(strm,ict.begin(),ict.end(),"");
}

} // |namespace lietype|


namespace weyl {

/*
  Synopsis: outputs w as a string of digits

  Originally this was without separators, with the following
  NOTE: this is satisfactory only if the rank is < 10; otherwise, use
  prettyprint::printWeylWord, that will give dot-separated lists.

  Unfortunately this operator is often called automatically, by template
  functions, so there is no place to insert that distinction into the code.
  Therefore we have inserted commas here. MvL
*/
  std::ostream& operator<< (std::ostream& strm, const WeylWord& w)
{
  if (w.size()==0)
    return strm << 'e';

  for (size_t j = 0; j < w.size(); ++j) {
    unsigned a = w[j]+1;
    strm << a;
    if (j+1<w.size()) strm << ',';
  }

  return strm;
}

} // |namespace weyl|


namespace matrix { // since |Weight| = |matrix::Vector<int>|

/*
  Synopsis: output of a lattice element.

  Using |seqPrint| it is output as a bracket-enclosed, comma-separated list.
*/

template<typename C>
std::ostream& operator<< (std::ostream& strm, const Vector<C>& v)
{
  std::ostringstream o; // accumulate in string for interpretation of width
  basic_io::seqPrint(o, v.begin(), v.end(), ",", "[", "]");
  return strm << o.str(); // now |strm.width()| is applied to whole vector
}

} // |namespace matrix|

namespace ratvec {

  std::ostream& operator<< (std::ostream& strm,
			    const RatWeight& w)
  {
    std::ostringstream o; // accumulate in string for interpretation of width
    o << w.numerator() << '/' << w.denominator();
    return strm << o.str(); // now |strm.width()| is applied to whole fraction
  }

} // |namespace ratvec|

namespace arithmetic {

  std::ostream& operator<< (std::ostream& strm, const Split_integer& s)
  { return strm << '(' << s.e() << '+' << s.s() << "s)"; }

} // |namespace arithmetic|

namespace repr {

  std::ostream& Rep_context::print (std::ostream& str,const StandardRepr& z)
    const
  {
    return
      str << "{x=" << z.x()
	  << ",lambda=" << lambda(z)
	  << ",nu=" << nu(z) << '}';
  }

  std::ostream& Rep_context::print (std::ostream& str,const SR_poly& P) const
  {
    for (SR_poly::const_iterator it=P.begin(); it!=P.end(); ++it)
      print(str << (it==P.begin() ?"":"+") << it->second, it->first)
	<< std::endl;
    return str;
  }

} // |namespace repr|



// binary input

namespace basic_io {


unsigned long long read_var_bytes(unsigned int n,std::istream& in)
{ switch(n)
  { case 1: return read_bytes<1>(in);
    case 2: return read_bytes<2>(in);
    case 3: return read_bytes<3>(in);
    case 4: return read_bytes<4>(in);
    case 5: return read_bytes<5>(in);
    case 6: return read_bytes<6>(in);
    case 7: return read_bytes<7>(in);
    case 8: return read_bytes<8>(in);
  default: throw std::runtime_error("Illegal read_var_bytes");
  }
}

void put_int (unsigned int val, std::ostream& out) { write_bytes<4>(val,out); }
void write_bytes(unsigned int n, unsigned long long val, std::ostream& out)
{ switch(n)
  { case 1: return write_bytes<1>(val,out);
    case 2: return write_bytes<2>(val,out);
    case 3: return write_bytes<3>(val,out);
    case 4: return write_bytes<4>(val,out);
    case 5: return write_bytes<5>(val,out);
    case 6: return write_bytes<6>(val,out);
    case 7: return write_bytes<7>(val,out);
    case 8: return write_bytes<8>(val,out);
    default: throw std::runtime_error("Illegal write__bytes");
  }
}


} // namespace basic_io

// Instantiations

// The following was an attempt to explicitly instantiate |seqPrint|, but
// giving a complete and non-redundant list proved to be architecture-dependent

// template std::ostream& basic_io::seqPrint // cartan_io
//  (std::ostream& strm,
//   const Partition::iterator::SubIterator& first,
//   const Partition::iterator::SubIterator& last,
//   const char* sep, const char* pre, const char* post);

// template std::ostream& basic_io::seqPrint // interactive
//  (std::ostream& strm,
//   const BitMap::iterator& first,
//   const BitMap::iterator& last,
//   const char* sep, const char* pre, const char* post);

// template std::ostream& basic_io::seqPrint // |kgb_io::printBruhatOrder|
//  (std::ostream& strm,
//   const set::EltList::const_iterator& first,
//   const set::EltList::const_iterator& last,
//   const char* sep, const char* pre, const char* post);

// template std::ostream& basic_io::seqPrint // |prettyprint::printWeylList|
//  (std::ostream& strm,
//   const std::vector<WeylWord>::iterator& first,
//   const std::vector<WeylWord>::iterator& last,
//   const char* sep, const char* pre, const char* post);

// template std::ostream& basic_io::seqPrint // |mainmode::roots_f| and others
//  (std::ostream& strm,
//   const WeightList::const_iterator& first,
//   const WeightList::const_iterator& last,
//   const char* sep, const char* pre, const char* post);

// template std::ostream& basic_io::seqPrint // |testprint::printComponents|
//  (std::ostream& strm,
//   const SmallBitVectorList::const_iterator& first,
//   const SmallBitVectorList::const_iterator& last,
//   const char* sep, const char* pre, const char* post);

// never instantiated
// template std::ostream& bitset::operator<<
//   (std::ostream& strm, const BitSet<constants::RANK_MAX>& b);

template std::ostream& polynomials::operator<<
  (std::ostream& strm, const Polynomial<int>& P);

template std::ostream& bitvector::operator<<
  (std::ostream& strm, const BitVector<constants::RANK_MAX>& b);

namespace matrix {
template std::ostream& operator<< (std::ostream& strm, const Vector<int>& v);
template std::ostream& operator<<
  (std::ostream& strm, const Vector<arithmetic::Numer_t>& v);

} // |namespace matrix|

} // namespace atlas
