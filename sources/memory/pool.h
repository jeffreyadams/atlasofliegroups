/*!
\file
\brief Class definitions and function declarations for the classes
SimplePool and Pool, which manage memory allocation for the templates
stlset::Set and stlvector::Vector respectively.

The idea is that these substitutes for the standard library templates
std::set and std::vector offer easier access for debugging
memory problems.
*/
/*
  This is pool.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

  See file main.cpp for full copyright notice
*/

/*
  Our memory pool always hands out blocks of T whose sizes are powers of
  2, regardless of the actual request. This seems to mesh well with the way
  push_back is usally implemented.

  In fact, we have one pool 
*/

#ifndef POOL_H
#define POOL_H

#include <cassert>
#include <cstddef>
#include <limits>
#include <vector>

#include "constants.h"

namespace atlas {

/******** type declarations *************************************************/

namespace pool {

class Pool;
class SimplePool;

}

/******** function declarations *********************************************/

namespace pool {

void memoryReport();

}

/******** type definitions **************************************************/

namespace pool {

class Pool {

 private:

  /*!  
    \brief Incremented each time Pool is instantiated, and
    decremented when a Pool is destroyed, so that it counts the number
    of Pool's.
  */
  static size_t instances;

  /*!
    \brief Incremented each time Pool is instantiated, but not
    decremented when a Pool is destroyed.
  */
  static size_t constructions;

  /*!
    \brief Report on the activity of Pool.
  */
  static const char* logfile;
  static bool done;

  /*!
    \brief Size in machine words of each object stored in this Pool.
  */
  const size_t d_atomSize;

  /*!
    \brief Each system memory request is of size
    (2^d_defaultRequest)(d_alignSize) machine words.
  */
  const size_t d_defaultRequest;

  /*
    \brief 
   */
  const size_t d_maxAlloc;
  const size_t d_alignSize;

  /*!
    \brief Which instance of Pool this is.
  */
  const size_t d_instance;

  /*!
    \brief Entry m (for 0 <= m < word size) is a pointer to the first
    free block of 2^m units of memory in the Pool.
  */
  void* d_free[constants::sizeBits];

  /*!
    \brief Entry m (for 0 <= m < word size) tells how many blocks of
    2^m units of memory have been assigned by this Pool.
  */
  size_t d_used[constants::sizeBits];

  /*!
    \brief Entry m (for 0 <= m < word size) tells how many blocks of
    2^m units of memory have been prepared for use by this Pool.
  */
  size_t d_allocated[constants::sizeBits];

  /*!
    \brief Entry \#j is a pointer to the beginning of the jth block
    (of size (d_Size)(2^d_defaultRequest)) assigned by the system
    to this Pool.

    Each assignment is of (d_alignSize)(2^m) bytes, for some m >=
    d_defaultRequest.
  */
  std::vector<void*> d_systemAllocs;

  void* newBlock(size_t m);
  void reportDestruction();

 public:

// constructors and destructors
  Pool(size_t a, size_t d);

  ~Pool();

// accessors
  size_t maxAlloc() const {
    return d_maxAlloc;
  }

  size_t instance() const {
    return d_instance;
  }

// modifiers
  void* allocate(size_t);
  void deallocate(void*, size_t);

// static members
  static void allowReport() {
    done = true;
  }

  static size_t numInstances() {
    return instances;
  }

  static void memoryReport();
};

/*!
\brief Pool for allocating memory in small units.

Each allocation will be of size d_atomSize.  Memory is requested from
the system in blocks of size (d_atomSize)(2^d_systemRequest); pointers
to the beginning of each such request are in d_systemAllocs.
*/
class SimplePool {

 private:

  /*!  
    \brief Incremented each time SimplePool is instantiated, and
    decremented when a SimplePool is destroyed, so that it counts the number
    of SimplePool's.
  */
  static size_t instances;

  /*!
    \brief Incremented each time SimplePool is instantiated, but not
    decremented when a SimplePool is destroyed.
  */
  static size_t constructions;

  /*!
    \brief Report on the activity of SimplePool.
  */
  static const char* logfile;
  static bool done;

  const size_t d_systemRequest;

  /*!
    \brief Size in machine words of each object stored in this SimplePool.
  */
  const size_t d_atomSize;

  /*!
    \brief Which instance of SimplePool this is.
  */
  const size_t d_instance;


  /*!
    \brief Pointer to the first free unit of memory in this SimplePool.
  */
  char* d_free;

  /*!
    \brief Pointer to the first byte beyond the available memory in
    this simple pool.

    Refers always to the most recent system allocation.
   */
  char* d_top;

  /*!
    \brief Tells how many units of memory have been assigned by this
    SimplePool.
  */
  size_t d_used;

  /*!
    \brief Tells how many units of memory have been prepared for use
    by this SimplePool.
  */
  size_t d_allocated;

  /*!
    \brief Entry \#j is a pointer to the beginning of the jth block
    (of size (d_atomSize)(2^d_systemRequest)) assigned by the system
    to this SimplePool. 
  */
  std::vector<void*> d_systemAllocs;

  void* newBlock();
  void reportDestruction();

 public:

// constructors and destructors
  SimplePool(size_t, size_t);

  ~SimplePool();

// accessors
  size_t maxAlloc() const {
    return 1ul;
  }

  size_t instance() const {
    return d_instance;
  }
// modifiers

   /*!
     \brief Returns a pointer to the beginning of n units of memory
     in this SimplePool.

     It is required that n = 1, or there will be serious memory problems.
   */
  void* allocate(size_t n);

  /*!
    \brief Returns n units of memory beginning at ptr to SimplePool.

    It is required that n = 1, or there will be serious memory problems.
  */
  void deallocate(void* ptr, size_t n);

// static members
  static void allowReport() {
    done = true;
  }

  static size_t numInstances() {
    return instances;
  }

  static void memoryReport();
};

}

}

#endif
