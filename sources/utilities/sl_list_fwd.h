/*
  This is sl_list_fwd.h, for a revolutionary singly linked list container type

  Copyright (C) 2016,2018 Marc van Leeuwen
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

template<typename T,
	 typename Container = mirrored_simple_list<T,std::allocator<T> > >
struct stack : public std::stack<T,Container>
{
#ifndef incompletecpp11
  using std::stack<T,Container>::stack;
#else
  template <typename... Args>
    stack(Args&&... args)
      : std::stack<T,Container> (std::forward<Args>(args)...)
  {}
#endif
}; // |struct stack|

template<typename T, typename Container = sl_list<T,std::allocator<T> > >
struct queue : public std::queue<T,Container>
{
#ifndef incompletecpp11
  using std::queue<T,Container>::queue;
#else
  template <typename... Args>
    queue(Args&&... args)
      : std::queue<T,Container> (std::forward<Args>(args)...)
  {}
#endif
}; // |struct queue|

} // |namespace cantainers|

} // |namespace atlas|


#endif
