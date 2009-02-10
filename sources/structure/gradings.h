/*!
\file
\brief Class definitions and function declarations for the class Status.
*/
/*
  This is gradings.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef GRADINGS_H  /* guard against multiple inclusions */
#define GRADINGS_H

#include <functional>

#include "gradings_fwd.h"

#include "rootdata_fwd.h"

#include "bitmap.h"
#include "bitset.h"
#include "constants.h"

namespace atlas {

/******** function declarations *********************************************/

namespace gradings {

  rootdata::RootSet max_orth(const rootdata::RootSet& non_compact,
			     const rootdata::RootSet& subsys,
			     const rootdata::RootSystem& rs);

  void transform_grading(gradings::Grading&,
			const rootdata::RootList&,
			const rootdata::RootSet&,
			const rootdata::RootSystem&);

}

/******** class definitions **************************************************/

namespace gradings {

  /*!
  \brief Describes a four-valued root attribute for each simple root: to be
  real, complex, imaginary compact or imaginary noncompact.

  This information fits nicely in two bits, and these bits can then be
  packed in a TwoRankFlags bitset (which for the default bound on the
  rank is 32 bits big.)
  */

class Status {

 private:

  bitset::TwoRankFlags d_flags;

 public:

  enum Value { Complex, ImaginaryCompact, Real, ImaginaryNoncompact };

// constructors and destructors
  Status() {} // all roots are |Complex| by default

  ~Status() {}

// accessors
  Value operator[] (size_t j) const {
    unsigned long f = d_flags.to_ulong();
    f &= constants::twoBitMask[j];
    f >>= 2*j;
    return static_cast<Value>(f);
  }

  bool isComplex(size_t j) const {
    j <<= 1;
    return not d_flags[j] and not d_flags[j+1];
  }

  bool isImaginaryCompact(size_t j) const {
    j <<= 1;
    return d_flags[j] and not d_flags[j+1];
  }

  bool isReal(size_t j) const {
    j <<= 1;
    return not d_flags[j] and d_flags[j+1];
  }

  bool isImaginaryNoncompact(size_t j) const {
    j <<= 1;
    return d_flags[j] and d_flags[j+1];
  }

  bool isImaginary(size_t j) const {
    j <<= 1;
    return d_flags[j];
  }

// manipulators

  void set(size_t j, Value v) {
    unsigned long a = v;
    bitset::TwoRankFlags b(a);
    b <<= 2*j;
    d_flags |= b;
  }

  void set_imaginary(size_t j, bool grading)
  {
    bitset::TwoRankFlags b(grading ? ImaginaryNoncompact : ImaginaryCompact);
    b <<= 2*j;
    d_flags |= b;
  }
}; // class Status


struct GradingCompare
  : public std::binary_function<const Grading& , const Grading& , bool>
{
  bool operator() (const Grading&, const Grading&);
};

}

}

#endif
