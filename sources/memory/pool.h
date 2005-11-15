/*
  This is pool.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

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

  static size_t instances;
  static size_t constructions;
  static const char* logfile;
  static bool done;

  const size_t d_atomSize;
  const size_t d_defaultRequest;
  const size_t d_maxAlloc;
  const size_t d_alignSize;

  const size_t d_instance;

  void* d_free[constants::sizeBits];
  size_t d_used[constants::sizeBits];
  size_t d_allocated[constants::sizeBits];

  std::vector<void*> d_systemAllocs;

  void* newBlock(size_t m);
  void reportDestruction();

 public:

// constructors and destructors
  Pool(size_t, size_t);

  ~Pool();

// accessors
  size_t maxAlloc() const {
    return d_maxAlloc;
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

class SimplePool {

 private:

  static size_t instances;
  static size_t constructions;
  static const char* logfile;
  static bool done;

  const size_t d_systemRequest;
  const size_t d_atomSize;
  const size_t d_instance;

  char* d_free;
  char* d_top;
  size_t d_used;
  size_t d_allocated;

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

}

}

#endif
