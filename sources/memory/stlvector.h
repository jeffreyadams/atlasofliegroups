/*
  This is stlvector.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef STLVECTOR_H  /* guard against multiple inclusions */
#define STLVECTOR_H

#include <set>
#include <vector>

#include "allocator.h"
#include "pool.h"

namespace atlas {

namespace stlvector {

template<typename T> struct Vector;

template<typename T> struct Vector {
  typedef std::vector<T, allocator::Allocator<T,pool::Pool> > type;
};

}

}

#endif
