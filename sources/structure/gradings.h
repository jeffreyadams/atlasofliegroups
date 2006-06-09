/*!
\file
  This is gradings.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef GRADINGS_H  /* guard against multiple inclusions */
#define GRADINGS_H

#include "gradings_fwd.h"

#include "rootdata_fwd.h"

#include "bitmap.h"
#include "bitset.h"
#include "constants.h"

namespace atlas {

/******** function declarations *********************************************/

namespace gradings {

  void findGrading(rootdata::RootSet&, const rootdata::RootList&, 
		      const rootdata::RootList&, const rootdata::RootDatum&);

  void gradingType(rootdata::RootList&, const Grading&, 
		   const rootdata::RootDatum&);

  void compactRoots(rootdata::RootList&, const Grading&, 
		    const rootdata::RootDatum&);

  bool isNonCompact(const rootdata::Root&, const Grading&);

  void makeGradings(GradingList&, const rootdata::RootDatum&);

  void noncompactRoots(rootdata::RootList&, const Grading&, 
		       const rootdata::RootDatum&);

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
  Status() {}

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
};

  /*!
  \brief Early version of Status, no longer instantiated.
  */
class FullStatus {

 private:

  bitmap::BitMap d_map;

 public:

// constructors and destructors

  explicit FullStatus(size_t n):d_map(n << 1) {}

  ~FullStatus() {}

// accessors
  Status::Value operator[] (size_t j) const {
    unsigned long a = d_map.range(j << 1,2);
    return static_cast<Status::Value> (a);
  }

// manipulators
  void set(unsigned long j, Status::Value v) {
    unsigned long a = v;
    d_map.setRange(j << 1,2,a);
  }

};

class GradingCompare {
 public:
  bool operator() (const Grading&, const Grading&);
};

}

}

#endif
