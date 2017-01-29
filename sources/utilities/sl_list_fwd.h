/*
  This is sl_list_fwd.h, for a revolutionary singly linked list container type

  Copyright (C) 2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef SL_LIST_FWD_H /* guard against multiple inclusions */
#define SL_LIST_FWD_H

#include <cstddef>
#include <cstdlib>
#include <memory>

// include to access the adapter templates, so we can replace default container
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
struct stack : public std::stack<T, mirrored_simple_list<T,Alloc> >
{
#ifndef incompletecpp11
  using std::stack<T, mirrored_simple_list<T,Alloc> >::stack;
#else
  template <typename... Args>
    stack(Args&&... args)
    : std::stack<T, mirrored_simple_list<T,Alloc> >
      (std::forward<Args>(args)...)
  {}
#endif
}; // |struct stack|

template<typename T,typename Alloc = std::allocator<T> >
struct queue : public std::queue<T, sl_list<T,Alloc> >
{
#ifndef incompletecpp11
  using std::queue<T, sl_list<T,Alloc> >::queue;
#else
  template <typename... Args>
    queue(Args&&... args)
    : std::queue<T, sl_list<T,Alloc> > (std::forward<Args>(args)...)
  {}
#endif
  queue() : std::queue<T, sl_list<T,Alloc> > ( sl_list<T,Alloc> {} ) {}
  queue(std::initializer_list<T> l)
    : std::queue<T, sl_list<T,Alloc> > ( sl_list<T,Alloc>(l) ) {}
}; // |struct queue|

} // |namespace cantainers|

} // |namespace atlas|


#endif
