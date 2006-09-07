/*!
\file
\brief Function declarations and type definition for Allocator<T,P>,
which allocates memory for objects of type T in pools of type P.
*/

/*
  This is allocator.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

#ifndef ALLOCATOR_H
#define ALLOCATOR_H

#include <cstddef>

#include "typenumber.h"

namespace atlas {

/******** constant declarations **********************************************/

namespace allocator {

static const size_t AllocDigits = 16;

}

/******** type definitions **************************************************/

namespace allocator {

template<typename T, typename P> class Allocator;

}

/******** function definitions ***********************************************/

namespace allocator {

void reportAllocation(size_t, size_t, const typenumber::TypeData&);
void reportConstruction(size_t, size_t, const typenumber::TypeData&);
void reportCopyConstruction(size_t, size_t, const typenumber::TypeData&);
void reportDeallocation(size_t, size_t, const typenumber::TypeData&);
void reportDestruction(size_t, size_t, const typenumber::TypeData&);
void reportDestructionError(size_t, size_t);
void reportHeterogeneousCopyConstruction(size_t);

}

/******** type declarations *************************************************/

namespace allocator {

template<typename T, typename P> class Allocator {

 private:

  static size_t d_instances;
  static size_t d_constructions;

  static P* d_pool;

 public:

// typedefs
  typedef T value_type;

  typedef T* pointer;
  typedef const T* const_pointer;

  typedef T& reference;
  typedef const T& const_reference;

  typedef ptrdiff_t difference_type;
  typedef size_t size_type;

  template<typename U> struct rebind {
    typedef Allocator<U,P> other;
  };

// constructors and destructors
  Allocator() throw(){
    if (d_instances == 0)
      d_pool = new P(sizeof(T),AllocDigits);
    ++d_instances;
  }

  Allocator(const Allocator&) throw() {
    ++d_instances;
  }

  template<typename U> Allocator(const Allocator<U,P>&) throw() {
    if (d_instances == 0)
      d_pool = new P(sizeof(T),AllocDigits);
    ++d_instances;
  }

  ~Allocator() {
    --d_instances;
    if (d_instances == 0) {
      delete d_pool;
    }
  }

// accessors
  size_type max_size() const throw() {
    return P::maxAlloc();
  }

// manipulators
  pointer allocate(size_type, const void* hint = 0);

  void deallocate(pointer, size_type);

  void construct(pointer, const_reference);

  void destroy(pointer);

};

template<typename T, typename P>
  P* Allocator<T,P>::d_pool = 0;

template<typename T, typename P>
  size_t Allocator<T,P>::d_instances = 0;

template<typename T, typename P>
  size_t Allocator<T,P>::d_constructions = 0;

// guarantee that it is instanceless
template<typename T, typename P>
  inline bool operator== (const Allocator<T,P>&, const Allocator<T,P>&) {
  return true;
}

template<typename T, typename P>
  inline bool operator!= (const Allocator<T,P>&, const Allocator<T,P>&) {
  return false;
}

}

}

#include "allocator_def.h"

#endif
