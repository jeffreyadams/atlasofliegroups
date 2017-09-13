/*
  This is descents.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Class definition for the class |DescentStatus|


#ifndef DESCENTS_H  /* guard against multiple inclusions */
#define DESCENTS_H

#include <cstring>

namespace atlas {
namespace descents {

/******** type definitions **************************************************/


/*			     |class DescentStatus|
   The descent status of each simple root for a single representation.

   For each simple root, there are eight possibilities for the corresponding
   descent status of a representation parameter, so this information could be
   packed in three bits per Weyl group generator. However, for memory efficiency
   we have larger fish to fry elsewhere, so we waste five bits per generator,
   and store each status in an octet (|unsigned char|).
*/
class DescentStatus
{
 public:
  // status for a single simple root; first $4$ are ascents, others descents
  enum Value { ComplexAscent, RealNonparity, ImaginaryTypeI, ImaginaryTypeII,
	       ImaginaryCompact, ComplexDescent, RealTypeII, RealTypeI };
 private:
// |d_data[j]| stores status for simple root |j| as a |Value|
  unsigned char d_data[constants::RANK_MAX];

// mask for the bit in |Value| characterising descents
  static constexpr auto DescentMask = 0x4u; // bit 2 set

// mask for bits in |Value| characterising |ComplexDescent| or |RealTypeI|
  static constexpr auto DirectRecursionMask = 0x5u; // bits 0 and 2 set

 public:
  static bool isDescent(Value v) { return (v & DescentMask)!=0; }
  static bool isDirectRecursion(Value v)
    { return (v & DirectRecursionMask) == DirectRecursionMask; }
  static Value dual(Value v) // cooresponding status in dual block
    { static const Value d[] =
	{ ComplexDescent, ImaginaryCompact, RealTypeII, RealTypeI,
	  RealNonparity, ComplexAscent, ImaginaryTypeI, ImaginaryTypeII };
      return d[v];
    }
  DescentStatus dual(unsigned int rank) const
  { DescentStatus result;
    for (unsigned int i=0; i<rank; ++i)
      result.d_data[i]=dual(static_cast<Value>(d_data[i]));
    return result;
  }

// constructors and destructors
  DescentStatus() { // sets statuses of all simple roots to 0 (ComplexAscent)
    std::memset(d_data,0,constants::RANK_MAX);
  }

  ~DescentStatus() {}

// copy and assignment (these copy statuses of all simple roots)
  DescentStatus(const DescentStatus& ds) {
    std::memcpy(d_data,ds.d_data,constants::RANK_MAX);
  }

  DescentStatus& operator=(const DescentStatus& ds) {
    std::memcpy(d_data,ds.d_data,constants::RANK_MAX);
    return *this;
  }

// accessors

/*!
\brief Returns descent status of simple root \#s.
*/
  Value operator[] (size_t s) const {
    return static_cast<Value> (d_data[s]); // cast converts integer to enum
  }

// manipulators

/*!
\brief Sets the descent status of simple root \#s to v.
*/
  void set(size_t s, Value v) {
    d_data[s] = v; // no cast needed here; enum value converts to integral type
  }
}; // |class DescentStatus|

}

}

#endif
