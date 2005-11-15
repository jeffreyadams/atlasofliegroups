/*
  This is layout.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef LAYOUT_H  /* guard against multiple inclusions */
#define LAYOUT_H

#include "latticetypes_fwd.h"
#include "lietype_fwd.h"

namespace atlas {

/******** type declarations **************************************************/

namespace layout {

struct Layout;

}

/******** type definitions **************************************************/

namespace layout {

struct Layout {

  lietype::LieType d_type;
  lietype::InnerClassType d_inner;
  latticetypes::WeightList d_basis;

// constructors and destructors
  Layout() {}

  Layout(const lietype::LieType& lt, 
		  const latticetypes::WeightList& b)
    :d_type(lt),d_basis(b) {}

  ~Layout() {}

// manipulators
  void swap (Layout& source) {
    d_type.swap(source.d_type);
    d_inner.swap(source.d_inner);
    d_basis.swap(source.d_basis);
  }

};

}

}

#endif
