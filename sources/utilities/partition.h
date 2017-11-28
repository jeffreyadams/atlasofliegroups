/*!
\file
\brief Class definitions and function declarations for class Partition.

The purpose of the class Partition is to represent the partition of a
finite set given by a group action.  A typical example is a Weyl group
acting on elements of order 2 in a torus.
*/
/*
  This is partition.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef PARTITION_H  /* guard against multiple inclusions */
#define PARTITION_H

#include <functional>
#include <vector>
#include <cstddef>

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

  // the template parameter F is the type of a binary function object
template<typename F>
  void makeOrbits(Partition&, F&, unsigned long, unsigned long);

}

/******** type definitions ***************************************************/

namespace partition {

  /*!
 \brief Partition of some set [0,n[ into classes.

  The partition is represented by a vector d_class of n unsigned longs,
  mapping values to a number in [0,s[ characterizing the class (where s is the
  number of classes), and by a vector d_classRep of s unsigned longs that
  inversely gives a representative element for each class.

  The class is equipped with members allowing it to be used as an unary
  function object ; this objects behaves as the map d_class.

  Main application is to the Fiber class: in that case n=2^m, with m at most
  RANK_MAX; then elements of [0,n[ are interpreted as elements of a vector
  space (Z/2Z)^m. Typical partitions are into the orbits of a Weyl group
  acting on this vector space (such partitions are computed by makeOrbits).
  */
class Partition {

 private:

  /*!
  \brief The number d_class[j] labels the class containing the value j.
  */
  std::vector<unsigned long> d_class;

  /*!
  \brief The value d_classRep[i] is some element of class \#i.
  */
  std::vector<unsigned long> d_classRep;

 public:

// types for unary function object

   /*!
   \brief Required to make Partition an adaptable unary function object

   Alternatively we could have got this by deriving the Parition class from
   std::unary_function<unsigned long,unsigned long>
   */
  typedef unsigned long argument_type;

  /*!
  \brief Required to make Partition an adaptable unary function object

   Alternatively we could have got this by deriving the Parition class from
   std::unary_function<unsigned long,unsigned long>
  */
  typedef unsigned long result_type;

// constructors and destructors
  Partition() {}

  /*!
  \brief The trivial partition of [0,n[ into a single class
  */
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

  /*!
\brief Returns the number of classes in the partition.
  */
  unsigned long classCount() const {
    return d_classRep.size();
  }

  /*!
  \brief Returns the number of an element belong to class \# c.
  */
  unsigned long classRep(unsigned long c) const {
    return d_classRep[c];
  }

  unsigned long classSize(unsigned long) const;

  /*!
  \brief Number of classes in the partition.
   */
  unsigned long size() const {
    return d_class.size();
  }

// manipulators

/*!
  \brief Adds value j to class \#c.

  Assumes that class \#c is already non-empty; that is, that j is not
  the first element of the class.  For the first element of a class,
  use newClass(j) instead.
*/
  void addToClass(unsigned long c, unsigned long j) {
    // should not be used for the first element in a class!
    // use newClass instead
    d_class[j] = c;
  }

  /*!
\brief Clears all entries of d_classRep.
  */
  void clear() {
    d_classRep.clear();
  }

  void newClass(unsigned long c);

  /*!
\brief Resizes the class to be a partition of [0,n[.
  */
  void resize(unsigned long n) {
    d_class.resize(n);
  }
};

class PartitionIterator {
  /*!
  The idea is that a PartitionIterator gives a convenient way of traversing
  the classes of a partition. It is intended that PartitionIterators are
  compared only if they refer to the same partition.

  The vector d_data contains the integers [0,n[, where n is the size of the
  set being partitioned, sorted in the order of the partition values. The
  vector d_stop contains an iterator pointing to the beginning of each
  class, and a final iterator pointing after the end of d_data. We then
  only need to return two consecutive elements in d_stop to bound a class.
  */
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
  typedef std::ptrdiff_t difference_type;
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
