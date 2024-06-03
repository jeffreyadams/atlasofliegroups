/*
  This is basic_io.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef BASIC_IO_H  /* guard against multiple inclusions */
#define BASIC_IO_H

#include <iosfwd>
#include <iostream>
#include <vector>

#include "Atlas.h"

#include "K_repr.h"
#include "repr.h"

/******** function declarations *********************************************/

namespace atlas {

namespace basic_io {
  template<typename I>
  std::ostream& seqPrint(std::ostream&, const I&, const I&,
			 const char* sep = ",", const char* pre = "",
			 const char* post = "");
} // |namespace basic_ip|


/* Non-member operators are defined in namespace of an operand, then
   argument-dependent lookup will find the operator if it's known at all.

   However, for |typedef| types, we need the namespace of the definiens
*/

namespace lietype {
  std::ostream& operator<< (std::ostream& strm, const SimpleLieType& slt);
  std::ostream& operator<< (std::ostream&, const LieType&);
  std::ostream& operator<< (std::ostream&, const InnerClassType&);
}

namespace weyl {
  std::ostream& operator<< (std::ostream&, const WeylWord&);
}

namespace matrix {
  template<typename C>
    std::ostream& operator<< (std::ostream& strm, const Vector<C>& v)
  {
    std::ostringstream o; // accumulate in string for interpretation of width
    basic_io::seqPrint(o, v.begin(), v.end(), ",", "[", "]");
    return strm << o.str(); // now |strm.width()| is applied to whole vector
  }
}

namespace polynomials {
template <typename C>
  std::ostream& operator<< (std::ostream& strm, const Polynomial<C>& P);
}

namespace ratvec {
  std::ostream& operator<< (std::ostream&, const RatWeight&);
}

namespace bitset {
template<unsigned int d>
  std::ostream& operator<< (std::ostream&, const BitSet<d>&);
}

namespace bitvector {
template<unsigned int dim>
  std::ostream& operator<< (std::ostream&, const BitVector<dim>&);
}

namespace arithmetic {

  std::ostream& print_split (std::ostream& strm, const Split_integer& s);

} // |namespace arithmetic|

namespace K_repr {

std::ostream& print_K_type
  (std::ostream& out, const K_type& val, const Rep_context& rc);

std::ostream& print_K_type_pol
  (std::ostream& out, const K_type_pol& val, const Rep_context& rc);

} // |namespace K_repr|

namespace repr {

std::ostream& print_stdrep
  (std::ostream& out, const StandardRepr& val, const Rep_context& rc);

std::ostream& print_stdrep_mod
  (std::ostream& out, const StandardReprMod& val, const Rep_context& rc);

std::ostream& print_SR_poly
(std::ostream& out, const repr::SR_poly& val, const Rep_context& rc);

} // |namespace repr|


// other functions

namespace basic_io {

template <unsigned int n>
unsigned long long read_bytes(std::istream& in);

template <unsigned int n>
void write_bytes(unsigned long long val, std::ostream& out);

unsigned long long read_var_bytes(unsigned int n,std::istream& in);

void put_int (unsigned int val, std::ostream& out);
void write_bytes(unsigned int n, unsigned long long val, std::ostream& out);


} // |namespace basic_io|

} // |namespace atlas|

#include "basic_io_def.h"

#endif
