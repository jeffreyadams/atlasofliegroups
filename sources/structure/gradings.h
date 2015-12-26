/*!
\file
\brief Class definitions and function declarations for the class Status.
*/
/*
  This is gradings.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef GRADINGS_H  /* guard against multiple inclusions */
#define GRADINGS_H


#include <functional> // base class |std::binary_function| for |GradingCompare|

#include "../Atlas.h"

#include "bitset.h"

 namespace atlas {

/******** function declarations *********************************************/

namespace gradings {

  RootNbrSet max_orth(const RootNbrSet& non_compact,
		      const RootNbrSet& subsys,
		      const RootSystem& rs);

  void transform_grading(Grading&,
			 const RootNbrList&,
			 const RootNbrSet&,
			 const RootSystem&);

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

class Status
{
  TwoRankFlags d_flags;

 public:

  enum Value { Complex, ImaginaryCompact, Real, ImaginaryNoncompact };

// constructors and destructors
  Status() {} // all roots are |Complex| by default

  ~Status() {}

// accessors
  Value operator[] (size_t j) const
  {
    j *= 2;
    unsigned u = d_flags[j] ? 1 : 0;
    if (d_flags[j+1])
      u+=2;
    return static_cast<Value>(u);
  }

  bool isComplex(size_t j) const
  {
    j *= 2;
    return not d_flags[j] and not d_flags[j+1];
  }

  bool isImaginaryCompact(size_t j) const
  {
    j *= 2;
    return d_flags[j] and not d_flags[j+1];
  }

  bool isReal(size_t j) const
  {
    j *= 2;
    return not d_flags[j] and d_flags[j+1];
  }

  bool isImaginaryNoncompact(size_t j) const
  {
    j *= 2;
    return d_flags[j] and d_flags[j+1];
  }

  bool isImaginary(size_t j) const
  {
    return d_flags[2*j];
  }

  bool operator==(const Status& other) const { return d_flags==other.d_flags; }
  bool operator!=(const Status& other) const { return d_flags!=other.d_flags; }

// manipulators

  void set(size_t j, Value v)
  {
    unsigned long a = v;
    TwoRankFlags b(a);
    b <<= 2*j;
    d_flags |= b;
  }

  void set_imaginary(size_t j, bool grading)
  {
    TwoRankFlags b(grading ? ImaginaryNoncompact : ImaginaryCompact);
    b <<= 2*j;
    d_flags |= b;
  }
}; // class Status


struct GradingCompare
  : public std::binary_function<const Grading& , const Grading& , bool>
{
  bool operator() (const Grading&, const Grading&);
};

} // |namespace gradings|

} // |namespace atlas|

#endif
