/*
  This is pool.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "pool.h"

#include <algorithm>
#include <cstdlib> // for malloc
#include <cstring> // for memset
#include <iostream>
#include <fstream>

namespace atlas {

/******** static member definitions ******************************************/

namespace pool {

  size_t Pool::instances = 0;
  size_t Pool::constructions = 0;
  const char* Pool::logfile = "pool.log";
  bool Pool::done = false;

  size_t SimplePool::instances = 0;
  size_t SimplePool::constructions = 0;
  const char* SimplePool::logfile = "simplepool.log";
  bool SimplePool::done = false;
}

/******** local definitions *************************************************/

namespace {

  using namespace pool;

  size_t lastBit(size_t);
  size_t reducedSize(size_t,size_t);

  struct MemBlock {
    MemBlock* d_next;
  };

  struct PoolDestruct {
    size_t d_instance;
    size_t d_allocated[constants::sizeBits];
    size_t d_systemAllocs;
  };

  inline bool operator< (const PoolDestruct& lhs, const PoolDestruct& rhs) {
    return lhs.d_instance < rhs.d_instance;
  }

  std::vector<PoolDestruct> poolDestructions;

  struct SimplePoolDestruct {
    size_t d_instance;
    size_t d_allocated;
    size_t d_systemAllocs;
  };

  inline bool operator< (const SimplePoolDestruct& lhs, 
			 const SimplePoolDestruct& rhs) {
    return lhs.d_instance < rhs.d_instance;
  }

  std::vector<SimplePoolDestruct> simplePoolDestructions;
}

/*****************************************************************************

        Chapter I -- the Pool class

  ... explain here when it is stable ...

******************************************************************************/

namespace pool {

Pool::Pool(size_t a, size_t d)
  :d_atomSize(a), 
   d_defaultRequest(d),
   d_maxAlloc((1ul << constants::sizeBits - 1ul)/a),
   d_alignSize(a > sizeof(MemBlock*) ? a : sizeof(MemBlock*)),
   d_instance(constructions++)

/*
  Synopsis : returns multiples of r bytes.
*/

{
  using namespace constants;

  ++instances;

  memset(d_free,0,sizeBits*sizeof(void*));
  memset(d_used,0,sizeBits*sizeof(size_t));
  memset(d_allocated,0,sizeBits*sizeof(size_t));
}

Pool::~Pool()

/*
  Synopsis: return the memory to the system.
*/

{
  reportDestruction();

  for (size_t j = 0; j < d_systemAllocs.size(); ++j)
    free(d_systemAllocs[j]);

  --instances;
  if (instances == 0 and done)
    memoryReport();
}

void* Pool::allocate(size_t n)

{
  if (n == 0 or n > d_maxAlloc)
    return 0;

  if (d_atomSize < sizeof(MemBlock*))
    n = reducedSize(n,d_atomSize);

  size_t m = 0;

  if (n > 1)
    m = lastBit(n-1) + 1;

  if (d_free[m] == 0) { // we need a new block
    d_free[m] = newBlock(m);
    if (d_free[m] == 0) // failure in memory allocation
      return 0;
  }

  MemBlock* block = static_cast<MemBlock*> (d_free[m]);
  d_free[m] = block->d_next;
  ++d_used[m];

  return static_cast<void *> (block);
}

void Pool::deallocate(void* ptr, size_t n)

/*
  Synopsis: frees the memory pointed by ptr.

  It is assumed that n is the number of atoms to which ptr is pointing.
*/

{
  if (n == 0)    return;

  // compute the list where ptr belongs
  if (d_atomSize < sizeof(MemBlock*))
    n = reducedSize(n,d_atomSize);

  size_t m = 0;

  if (n > 1)
    m = lastBit(n-1) + 1;

  // relink the list, putting the freed block in front
  MemBlock *block = static_cast<MemBlock*> (ptr);
  block->d_next = static_cast<MemBlock*> (d_free[m]);
  d_free[m] = ptr;
  --d_used[m];

  return;
}

void Pool::memoryReport()

/*
  Synopsis: writes a report on the various pool constructions and destructions
  in the memory log file.
*/

{
  using namespace constants;

  std::ofstream log(logfile);

  log << "total number of pool instances: " << Pool::constructions
	   << std::endl;

  std::sort(poolDestructions.begin(),poolDestructions.end());

  for (size_t j = 0; j < poolDestructions.size(); ++j) {
    const PoolDestruct& pd = poolDestructions[j];
    log << std::endl;
    size_t a = pd.d_systemAllocs;
    log << "pool #" << pd.d_instance << ": "
	     << a << " system allocation" << (a == 1ul ? "" : "s") 
	     << std::endl;
    log << "used: ";
    for (size_t i = 0; i < sizeBits; ++i) {
      log << pd.d_allocated[i];
      if (i < sizeBits-1)
	log << ",";
    }
    log << std::endl;
  }

  return;
}

void* Pool::newBlock(size_t m)

/*
  Synopsis: allocates a new block for a memory request of 2^m units of size
  max(d_atomSize,sizeof(MemBlock*))
*/

{
  for (size_t j = m+1; j <= d_defaultRequest; ++j) {
    if (d_free[j]) /* split this block up */
      {
	void* vptr = d_free[j];
	char* ptr = static_cast<char*> (vptr);
	d_free[j] = static_cast<MemBlock*> (vptr)->d_next;
	--d_allocated[j];
	for (size_t i = m; i < j; ++i) {
	  d_free[i] = ptr + (d_alignSize << i);
	  MemBlock* block = static_cast<MemBlock *> (d_free[i]);
	  block->d_next = 0;
	  ++d_allocated[i];
	}
	MemBlock* block = static_cast<MemBlock*> (d_free[m]);
	block->d_next = static_cast<MemBlock*> (vptr);
	block->d_next->d_next = 0;
	++d_allocated[m];
	return d_free[m];
      }
  }

  /* if we get here we need more memory from the system */

  if (m >= d_defaultRequest) { // get block directly

    d_free[m] = malloc(d_alignSize << m);

    if (d_free[m] == 0) { // memory allocation failed
      return 0;
    }

    d_systemAllocs.push_back(d_free[m]);

    MemBlock* block = static_cast<MemBlock *> (d_free[m]);
    block->d_next = 0;
    ++d_allocated[m]; 
   
    return d_free[m];
  }

  /* if we get here, we allocate and split up a block of size 2^d */

  void* vptr = malloc(d_alignSize << d_defaultRequest);

  if (vptr == 0) { // memory allocation failed
    return 0;
  }

  d_systemAllocs.push_back(vptr);

  char* ptr = static_cast<char*> (vptr);

  for (unsigned j = m; j < d_defaultRequest; ++j) {
    d_free[j] = ptr + (d_alignSize << j);
    MemBlock* block = static_cast<MemBlock *> (d_free[j]);
    block->d_next = 0;
    ++d_allocated[j];
  }

  MemBlock* block = static_cast<MemBlock *> (d_free[m]);
  block->d_next = static_cast<MemBlock *> (vptr);
  ++d_allocated[m];

  block->d_next->d_next = 0;

  return d_free[m];
}

void Pool::reportDestruction()

/*
  Synopsis: writes down data about the memory usage of the pool.
*/

{
  using namespace constants;

  PoolDestruct pd;

  pd.d_instance = d_instance;
  pd.d_systemAllocs = d_systemAllocs.size();
  memcpy(pd.d_allocated,d_allocated,sizeBits*sizeof(size_t));

  poolDestructions.push_back(pd);

  return;
}

}

/*****************************************************************************

        Chapter II -- The SimplePool class

  ... explain here when it is stable ...

******************************************************************************/

namespace pool {

SimplePool::SimplePool(size_t a, size_t d)
  :d_systemRequest(d),
   d_atomSize(a),
   d_instance(constructions++),
   d_free(0),
   d_top(0),
   d_used(0),
   d_allocated(0)

{
  ++instances;
}

SimplePool::~SimplePool()

/*
  Synopsis: return the memory to the system.
*/

{
  reportDestruction();

  for (size_t j = 0; j < d_systemAllocs.size(); ++j)
    free(d_systemAllocs[j]);

  --instances;
  if (instances == 0 and done)
    memoryReport();
}

void* SimplePool::allocate(size_t n)

{
  assert(n == 1);

  if (d_free == d_top) { // we need a new block
    d_free = static_cast<char*> (newBlock());
  }

  void* vptr = d_free;
  d_free += d_atomSize;
  ++d_used;

  return vptr;
}

void SimplePool::deallocate(void* ptr, size_t n)

/*
  Synopsis: does nothing!

  In a SimplePool, memory is never deallocated. It is returned to the system
  on destruction.
*/

{}

void SimplePool::memoryReport()

/*
  Synopsis: writes a report on the various pool constructions and destructions
  in the memory log file.
*/

{
  using namespace constants;

  std::ofstream log(logfile);

  log << "total number of pool instances: " << SimplePool::constructions
	   << std::endl;

  std::sort(simplePoolDestructions.begin(),simplePoolDestructions.end());

  for (size_t j = 0; j < simplePoolDestructions.size(); ++j) {
    const SimplePoolDestruct& spd = simplePoolDestructions[j];
    log << std::endl;
    size_t a = spd.d_systemAllocs;
    log << "pool #" << spd.d_instance << ": "
	<< a << " system allocation" << (a == 1ul ? "" : "s") 
	<< std::endl;
    log << "used: " << spd.d_allocated << std::endl;
  }

  return;
}

void* SimplePool::newBlock()

/*
  Synopsis: allocates a new block 
*/

{
  void* vptr = malloc(d_atomSize << d_systemRequest);
  d_systemAllocs.push_back(vptr);

  if (vptr == 0) // failure in memory allocation
    return 0;

  d_free = static_cast<char*> (vptr);
  d_top = d_free + (d_atomSize << d_systemRequest);
  d_allocated += (1ul << d_systemRequest);

  return vptr;
}

void SimplePool::reportDestruction()

/*
  Synopsis: writes down data about the memory usage of the pool.
*/

{
  using namespace constants;

  SimplePoolDestruct spd;

  spd.d_instance = d_instance;
  spd.d_systemAllocs = d_systemAllocs.size();
  spd.d_allocated = d_used;

  simplePoolDestructions.push_back(spd);

  return;
}

}

/*****************************************************************************

        Chapter III -- Functions declared in pool.h

  ... explain here when it is stable ...

******************************************************************************/

namespace pool {

void memoryReport() 

/*
  Synopsis: prints out a report on the memory usage.

  More precisely, it prints it out now if all pools have already been
  destroyed; it allows the printing after the destruction of the last
  instance if some pools are still in use (even if this function is called
  on exit, there may very well still be static objects alive that are held
  by memory pools; this is even the common situation.)
*/

{
  if (Pool::numInstances() == 0) // all pools hace been destroyed
    Pool::memoryReport();
  else // it will be written when the last instance goes away
    Pool::allowReport();

  if (SimplePool::numInstances() == 0) // all pools hace been destroyed
    SimplePool::memoryReport();
  else // it will be written when the last instance goes away
    SimplePool::allowReport();

  return;
}

void reportConstruction(const char* adj, size_t n)

{
  std::cout << "construction of " << adj << "pool #" << n << std::endl;
}

void reportCopyConstruction(const char* adj, size_t n)
{
  std::cout << "copy construction of " << adj << "pool #" << n << std::endl;
}

void reportDestruction(const char* adj, size_t n)

{
  std::cout << "destruction of " << adj << "pool #" << n << std::endl;
}

}

/*****************************************************************************

        Chapter IV -- Auxiliary functions

  ... explain here when it is stable ...

******************************************************************************/

namespace {

size_t lastBit (size_t d_n)

/*
  Synopsis: straightforward implementation to get the last set bit in n.
*/

{
  if (d_n == 0)
    return std::numeric_limits<size_t>::digits;

  size_t m = 0;

  for (size_t n = d_n >> 1; n; n >>= 1)
    ++m;

  return m;
}

size_t reducedSize(size_t n, size_t a)

/*
  Synopsis: given a < sizeof(MemBlock*), tells how many units of MemBlock*
  have to be allocated to cover n units of a.
*/

{
  return (n*a)/sizeof(MemBlock*) + ((n*a) % sizeof(MemBlock*) ? 1 : 0);
}

}

}
