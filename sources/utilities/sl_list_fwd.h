/*
  This is sl_list.h, a revolutionary singly linked list container type
*/
/*
  Copyright (C) 2015 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef SL_LIST_FWD_H /* guard against multiple inclusions */
#define SL_LIST_FWD_H

#include <cstddef>
#include <cstdlib>
#include <memory>

#include <stack>
#include <queue>

namespace atlas {

namespace containers {

template<typename T,typename Alloc = std::allocator<T> >
  class simple_list;
template<typename T,typename Alloc = std::allocator<T> >
  class sl_list;

template<typename T, typename Alloc = std::allocator<T> >
  struct sl_list_const_iterator;
template<typename T,typename Alloc = std::allocator<T> >
  class sl_list_iterator;

template<typename T,typename Alloc = std::allocator<T> >
  class mirrored_simple_list; // trivial adapter, to allow use with |std::stack|

template<typename T,typename Alloc = std::allocator<T> >
  class mirrored_sl_list; // trivial adapter, to allow use with |std::stack|

template<typename T,typename Alloc = std::allocator<T> >
  using stack = std::stack<T, mirrored_simple_list<T,Alloc> >;

template<typename T,typename Alloc = std::allocator<T> >
  using queue = std::queue<T, sl_list<T,Alloc> >;

} // |namespace cantainers|

} // |namespace atlas|


#endif

