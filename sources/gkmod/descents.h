/*!
\file
\brief Class definition for the class DescentStatus.
*/

/*
  This is descents.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef DESCENTS_H  /* guard against multiple inclusions */
#define DESCENTS_H

#include "descents_fwd.h"

#include "constants.h"

namespace atlas {

/******** constant declarations *********************************************/

namespace descents {}

/******** function declarations *********************************************/

/******** type definitions **************************************************/

namespace descents {


  /*!
\brief Describes the descent status of each simple root for a single
representation.

    There are eight possibilities for the descent status of a representation
    parameter w.r.t. a simple reflection, so this information could be packed
    in three bits. However, this would lead to having some packets lie across
    word boundaries, with the ensuing complications. Four bits is a good
    choice; here we take the lazy way of using even eight bits, as this
    makes the reading even a bit easier. We might come back to  four at
    some later change---this should not require a change in user interface.
  */
class DescentStatus {

 private:

  /*!
\brief Value of byte \#j specifies the descent status of simple root
\#j.

Value should be 0 through 7; values 0 through 3 are ascents, and 4
through 7 are descents.
  */
  unsigned char d_data[constants::RANK_MAX];
  static const unsigned ValMask = constants::charBits - 1;

  /*!
\brief Bitwise "and" of Value with this is non-zero if Value is one of
ImaginaryCompact, ComplexDescent, RealTypeII, or RealTypeI (numbers 4--7)
  */
  static const unsigned DescentMask = 0x4ul;

  /*!
\brief Bitwise "and" of Value with this is equal to this if Value is either
ComplexDescent (5) or RealTypeI (7)
  */
  static const unsigned DirectRecursionMask = 0x5ul;

 public:

  enum Value { ComplexAscent, RealNonparity, ImaginaryTypeI, ImaginaryTypeII,
	       ImaginaryCompact, ComplexDescent, RealTypeII, RealTypeI };

  /*!
\brief Tests whether Value is 4 through 7.  These are the descents.

The simple roots passing this test comprise the tau invariant for the
representation.
  */
  static bool isDescent(Value v) {
    return (v & DescentMask)!=0;
  }

  /*!
\brief Tests whether both bits of DirectRecursionMask are set

In the case of a complex descent or a real type I descent there is a simple
recursion formula for the KL element.
  */
  static bool isDirectRecursion(Value v) {
    return (v & DirectRecursionMask) == DirectRecursionMask;
  }

// constructors and destructors
  DescentStatus() { // sets statuses of all simple roots to 0 (ComplexAscent)
    memset(d_data,0,constants::RANK_MAX);
  }

  ~DescentStatus() {}

// copy and assignment (these copy statuses of all simple roots)
  DescentStatus(const DescentStatus& ds) {
    memcpy(d_data,ds.d_data,constants::RANK_MAX);
  }

  DescentStatus& operator=(const DescentStatus& ds) {
    memcpy(d_data,ds.d_data,constants::RANK_MAX);
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
};

}

}

#endif
