/*!
  \file
  \brief Class definition for Allocator<T,P>, which allocates memory
  for objects of type T in pools of type P.
*/
/*
  This is allocator_def.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#include "error.h"
#include "pool.h"
#include "typenumber.h"

namespace atlas {

namespace allocator {

template<typename T, typename P> 
typename Allocator<T,P>::pointer 
  Allocator<T,P>::allocate(typename Allocator<T,P>::size_type n, const void*)

/*!
  \brief Allocates raw memory for n instances of T.

  NOTE: very preliminary.
*/

{
  using namespace error;
  using namespace pool;
  using namespace typenumber;

  if (n == 0)
    return 0;

  T* ptr = static_cast<T*> (d_pool->allocate(n));

  if (ptr == 0)  { // allocation failed
    throw MemoryOverflow();
  }

  // reports details of memory allocation for debugging purposes.
  // reportAllocation(n,sizeof(value_type),d_pool->instance());

  return ptr;
}

template<typename T, typename P>
void Allocator<T,P>::deallocate(typename Allocator<T,P>::pointer ptr, 
			      typename Allocator<T,P>::size_type n)

/*
  Synopsis: deallocates raw memory.
*/

{
  using namespace typenumber;

  if (n == 0)
    return;

  d_pool->deallocate(ptr,n);

  // reports details of memory deallocation for debugging purposes.
  // reportDeallocation(n,sizeof(value_type),d_pool->instance());

  return;
}

template<typename T, typename P>
void Allocator<T,P>::construct(typename Allocator<T,P>::pointer p, 
			     typename Allocator<T,P>::const_reference val)

/*
  Synopsis: calls the constructor on raw memory, with placement new.
*/

{           
  new((void*)p) T(val);
}

template<typename T, typename P>
void Allocator<T,P>::destroy(typename Allocator<T,P>::pointer p)

/*
  Synopsis: calls the destructor on *p.
*/

{          
  p->~T();
}

}

}

