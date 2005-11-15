/*
  This is descents.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

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

class DescentStatus {

 private:

  /*
    There are eight possibilities for the descent status of a representation
    parameter w.r.t. a simple reflection, so this information could be packed
    in three bits. However, this would lead to having some packets lie across
    word boundaries, with the ensuing complications. Four bits is a good
    choice; here we take the lazy way of using even eight bits, as this
    makes the reading even a bit easier. We might come back to  four at
    some later change---this should not require a change in user interface.
  */

  unsigned char d_data[constants::RANK_MAX];
  static const unsigned ValMask = constants::charBits - 1;
  static const unsigned DescentMask = 0x4ul;
  static const unsigned DirectRecursionMask = 0x5ul;

 public:

  enum Value { ComplexAscent, RealNonparity, ImaginaryTypeI, ImaginaryTypeII,
	       ImaginaryCompact, ComplexDescent, RealTypeII, RealTypeI };

  static bool isDescent(Value v) {
    return v & DescentMask;
  }

  static bool isDirectRecursion(Value v) {
    return (v & DirectRecursionMask) == DirectRecursionMask;
  }

// constructors and destructors
  DescentStatus() {
    memset(d_data,0,constants::RANK_MAX);
  }

  ~DescentStatus() {}

// copy and assignment
  DescentStatus(const DescentStatus& ds) {
    memcpy(d_data,ds.d_data,constants::RANK_MAX);
  }

  DescentStatus& operator=(const DescentStatus& ds) {
    memcpy(d_data,ds.d_data,constants::RANK_MAX); 
    return *this;
  }

// accessors
  Value operator[] (size_t s) const {
    return static_cast<Value> (d_data[s]);
  }

// manipulators
  void set(size_t s, Value v) {
    d_data[s] = v;
  }
};

}

}

#endif
