/*
  This is partition.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#ifndef PARTITION_H  /* guard against multiple inclusions */
#define PARTITION_H

#include <functional>
#include <vector>

#include "tags.h"

/******** type declarations **************************************************/

namespace atlas {

namespace partition {

  class Partition;
  class PartitionIterator;

  struct End{};

}

/******** function declarations **********************************************/

namespace partition {

template<typename F> 
  void makeOrbits(Partition&, F&, unsigned long, unsigned long);

}

/******** type definitions ***************************************************/

namespace partition {

class Partition {

 private:

  std::vector<unsigned long> d_class;
  std::vector<unsigned long> d_classRep;

 public:

// types for unary function object
  typedef unsigned long argument_type;
  typedef unsigned long result_type;

// constructors and destructors
  Partition() {}

  explicit Partition(unsigned long n):d_class(n),d_classRep(0) {}

  explicit Partition(std::vector<unsigned long>&);

  Partition(std::vector<unsigned long>&, tags::UnnormalizedTag);

  ~Partition() {}

// copy and assignment
  void swap(Partition&);

// accessors
  unsigned long operator()(unsigned long j) const {
    return d_class[j];
  }

  bool operator== (const Partition& other) const {
    return d_class == other.d_class;
  }

  unsigned long classCount() const {
    return d_classRep.size();
  }

  unsigned long classRep(unsigned long c) const {
    return d_classRep[c];
  }

  unsigned long classSize(unsigned long) const;

  unsigned long size() const {
    return d_class.size();
  }

// manipulators
  void addToClass(unsigned long c, unsigned long j) {
    // should not be used for the first element in a class!
    // use newClass instead
    d_class[j] = c;
  }

  void clear() {
    d_classRep.clear();
  }

  void newClass(unsigned long c);

  void resize(unsigned long n) {
    d_class.resize(n);
  }
};

class PartitionIterator {

 private:

  std::vector<unsigned long> d_data;
  std::vector<std::vector<unsigned long>::const_iterator> d_stop;
  std::pair<std::vector<unsigned long>::const_iterator,
    std::vector<unsigned long>::const_iterator> d_range;
  std::vector<std::vector<unsigned long>::const_iterator>::const_iterator
    d_currentEnd;

 public:

// associated types
  typedef std::forward_iterator_tag iterator_category;
  typedef std::pair<std::vector<unsigned long>::const_iterator,
    std::vector<unsigned long>::const_iterator> value_type;
  typedef ptrdiff_t difference_type;
  typedef const value_type* pointer;
  typedef const value_type& reference;

// constructors and destructors
  PartitionIterator() {}

  explicit PartitionIterator(const Partition&);

  PartitionIterator(const PartitionIterator& i,End) {
    d_range.first = i.d_data.end();
    d_range.second = i.d_data.end();
  }

  ~PartitionIterator() {}

// accessors
  bool operator== (const PartitionIterator& i) const {
    return d_range.first == i.d_range.first;
  }

  bool operator!= (const PartitionIterator& i) const {
    return d_range.first != i.d_range.first;
  }

  reference operator* () const {
    return d_range;
  }

  pointer operator-> () const {
    return &d_range;
  }

  bool operator() () const {
    return d_range.first != d_data.end();
  }

// manipulators
  PartitionIterator& operator++ ();

  PartitionIterator operator++ (int) {
    PartitionIterator tmp(*this); 
    ++(*this); 
    return tmp;
  }

};

}

}

/******** template definitions ***********************************************/

#include "partition_def.h"

#endif
