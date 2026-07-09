// This is node_allocator.h, defining an allocator for node-based containers

/*
  Copyright (C) 2018 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef NODE_ALLOCATOR_H /* guard against multiple inclusions */
#define NODE_ALLOCATOR_H

#include <cstddef> // for |std::size_t|

namespace atlas {

namespace containers {

template<typename T>
  class node_allocator // stateless allocator only dealing in single |T|s
{
public:
  using value_type = T; // what we allocate for

private:
  static constexpr std::size_t block_size = 0x400; // must be at least 2

  union space // free space for one |T|, linked to next such free space
  {
    T storage; // when given out, it will be so as storage
    space* next; // but otherwise, this link points to next, or is |nullptr|

    space() : next(nullptr) {} // default constructor
    ~space() {} // destructor must be defined, but does nothing
  };

  // data, all |static|
  static space* pool;

public:
  // constructors
  node_allocator () {}
  template<typename U> node_allocator (const node_allocator<U>&) {}

  T* allocate (std::size_t n)
  { // assert(n==1);
    if (pool!=nullptr) // then we can allocate from our free list
    {
      auto result = &pool->storage;
      pool = pool->next; // pop entry from free list
      return result;
    }
    // when our free list runs out, extend it a whole block at a time:
    pool = new space[block_size]; // these blocks will never be |delete[]|d
    size_t i=block_size;
    auto result = &pool[--i].storage;
    pool[--i].next = nullptr;
    while (i-->0) // link together all other entries into an ascending chain
      pool[i].next = &pool[i+1];
    return result;
  }

  void deallocate (T* p,std::size_t n)
  { // assert(n==1);
    auto ps = reinterpret_cast<space*>(p); // formally UB; it would be OK with
    // |reinterpret_cast| on allocation; but then _using_ the storage would be UB
    ps->next = pool;
    pool = ps;
  };
}; // |template<typename T> class node_allocator|

template<typename T>
typename node_allocator<T>::space*
node_allocator<T>::pool = nullptr; // initialise one list for each |T|

template<typename T>
bool operator==(const node_allocator<T>&,const node_allocator<T>&)
  { return true; }
template<typename T>
bool operator!=(const node_allocator<T>&,const node_allocator<T>&)
  { return false; }

} // |namespace containers|

} // |namespace atlas|


#endif
