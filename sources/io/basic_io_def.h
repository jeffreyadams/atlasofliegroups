/*
  This is basic_io_def.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include <iostream>

/*****************************************************************************

        Chapter I -- Template functions defined in basic_io.h

******************************************************************************/

namespace atlas {

namespace bitvector {

/*
  Prints the bits of v on strm left-to-right
*/
template<size_t dim>
  std::ostream& operator<< (std::ostream& strm,  const BitVector<dim>& v)
{
  for (size_t i=0; i<v.size(); ++i)
    strm << (v[i]?'1':'0');

  return strm;
}

} // |namespace bitvector|

namespace bitset {

/*
  Prints the bits of b on strm left-to-right
*/
template<size_t d>
  std::ostream& operator<< (std::ostream& strm, const BitSet<d>& b)
{
  for (size_t i = 0; i < d; ++i)
    strm << (b[i]?'1':'0');

  return strm;
}

} // |namespace bitset|

namespace basic_io {

template<typename I>
std::ostream& seqPrint(std::ostream& strm, const I& first, const I& last,
		       const char* sep, const char* pre, const char* post)

/*
  This is a function for sequential output. It is assumed that I is an
  iterator type, pointing to a type for which the << operator is defined.
  Then we output the elements of the range [first,last[ as a list,
  with prefix pre, postfix post, and separator sep.
*/

{
  strm << pre;
  bool firstElt = true;

  for (I i = first; i != last; ++i) {
    if (firstElt)
      firstElt = false;
    else
      strm << sep;
    strm << *i;
  }

  strm << post;

  return strm;
}

// binary input
template <unsigned int n>
inline unsigned long long read_bytes(std::istream& in)
{
  return static_cast<unsigned char>(in.get())+(read_bytes<n-1>(in)<<8);
}

template<>
inline unsigned long long read_bytes<1>(std::istream& in)
{
  return static_cast<unsigned char>(in.get());
}

// binary output
template <unsigned int n>
inline void write_bytes(unsigned long long val, std::ostream& out)
{
  out.put(char(val&0xFF)); write_bytes<n-1>(val>>8,out);
}

template <>
inline void write_bytes<1>(unsigned long long val, std::ostream& out)
{
  out.put(char(val));
}


} // namespace basic_io

} // namespace atlas
