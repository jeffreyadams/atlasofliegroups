/*
  This is lietype.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef LIETYPE_H  /* guard against multiple inclusions */
#define LIETYPE_H

#include "lietype_fwd.h"

#include "latticetypes_fwd.h"

namespace atlas {

/******** constant declarations **********************************************/

namespace lietype {

  const char* const typeLetters = "ABCDEFGT";
  const char* const innerClassLetters = "Ccesu";

}

/******** function declarations **********************************************/

namespace lietype {

  bool checkRank(const TypeLetter&, size_t);

  void dualLieType(LieType&, const LieType&);

  void dualInnerClassType(InnerClassType&, const InnerClassType&,
			  const LieType& lt);

  void involution(latticetypes::LatticeMatrix&, const lietype::LieType&, 
		  const lietype::InnerClassType&);

  size_t rank(const LieType&);

  inline size_t rank(const SimpleLieType& slt) {
    return slt.second;
  }

  size_t semisimpleRank(const LieType&);

  inline size_t semisimpleRank(const SimpleLieType& slt) {
    return slt.first == 'T' ? 0UL : slt.second;
  }

  inline TypeLetter type(const SimpleLieType& slt) {
    return slt.first;
  }

}

}

#endif
