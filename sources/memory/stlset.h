/*
  This is stlset.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef STLSET_H  /* guard against multiple inclusions */
#define STLSET_H

#include <set>
#include <vector>

#include "allocator.h"
#include "pool.h"

namespace atlas {

namespace stlset {

template<typename T, typename Compare = std::less<T> > struct Set;

template<typename T, typename Compare> struct Set {
  typedef std::set<T, Compare, allocator::Allocator<T,pool::SimplePool> > type;

};

}

}

#endif
