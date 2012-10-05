/*
  This is stlset.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations 

  For license information see the LICENSE file
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
