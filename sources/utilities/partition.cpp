/*
  Implementation of Partition.

  The purpose of the class Partition is to compute the partition of a
  finite set given by a group action.  A typical example is a Weyl group
  acting on elements of order 2 in a torus.
*/
/*
  This is partition.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "partition.h"

#include <algorithm>
#include <map>

#include "comparison.h"

/*****************************************************************************

   ... describe here when it is stable ....  (Main application:
computing the partition of a finite set under a group action. Typical
example is a Weyl group acting on elements of order two in a torus.)

******************************************************************************/

/*****************************************************************************

        Chapter I -- The Partition class.

******************************************************************************/

namespace atlas {

namespace partition {

/******** constructors and destructors ***************************************/

Partition::Partition(std::vector<unsigned long>& f)
  : d_class(f.size())
  , d_classRep()
{ /* at this point our object is already in a valid state (although not one
     describing the correct partition), so we feel free to call the methods
     newClass and addToClass for the object under construction */

  unsigned long s = 0; // current number of classes

  std::map<unsigned long,unsigned long> val; // bijective map im(f)->[0,s[

  for (size_t i=0; i<f.size(); ++i)
  {
    // tentatively map f[i] to a new class number |s|
    auto trial = val.insert(std::make_pair(f[i],s));
    if (trial.second) // insertion succeeded, |f[i]| is new
    { // found a new value
      new_class(i); // add a new class containing i to the partition
      ++s;
    }
    else // |f[i]| had already been seen, find its class and record |i| in it
      addToClass(trial.first->second,i);
  }
}

/*
  Count the number of elements in class \#c.

  NOTE: Straightforward implementation. Successively computing |classSize| for
  all classes would cost more time then necessary.
*/
Partition::Partition(std::vector<unsigned long>& f, tags::UnnormalizedTag)
  : d_class(f) // just use (a copy of) |f| as class vector
  , d_classRep()
{
  // find class representatives
  std::map<unsigned long,unsigned long> val; // associates (class number,repr)

  for (size_t j = 0; j < f.size(); ++j)
    val.insert(std::make_pair(f[j],j)); // ignore failures, first value sticks

  d_classRep.resize(val.size()); // resize to size of image of |f|

  std::map<unsigned long,unsigned long>::iterator val_end = val.end();

  // now convert |val| from an associative to an ordinary array
  for (std::map<unsigned long,unsigned long>::iterator i = val.begin();
       i != val_end; ++i)
    d_classRep[i->first] = i->second;
}

/******** copy and assignment ************************************************/

void Partition::swap(Partition& other)
{
  d_class.swap(other.d_class);
  d_classRep.swap(other.d_classRep);
}

/******** manipulators *******************************************************/

// Create a fresh class for element |j|, and return its number
unsigned long Partition::new_class(unsigned long j)
{
  d_class[j] = d_classRep.size();
  d_classRep.push_back(j);
  return d_class[j];
}

/*
  Count the number of elements in class \#c.

  NOTE: Straightforward implementation. Successively computing |classSize| for
  all classes with this method would cost more time then necessary.
*/
unsigned long Partition::classSize(unsigned long c) const
{
  unsigned long n = 0;

  for (size_t j = 0; j < d_class.size(); ++j)
    if (d_class[j] == c)
      ++n;

  return n;
}

}

/*****************************************************************************

        Chapter II -- The PartitionIterator class.

  The idea is that a PartitionIterator gives a convenient way of traversing
  the classes of a partition. It is intended that PartitionIterators are
  compared only if they refer to the same partition.

  The vector d_data contains the integers [0,n[, where n is the size of the
  set being partitioned, sorted in the order of the partition values. The
  vector d_stop contains an iterator pointing to the beginning of each
  class, and a final iterator pointing after the end of d_data. We then
  only need to return two consecutive elements in d_stop to bound a class.

******************************************************************************/

namespace partition {

/*!
  Initializes a PartitionIterator ready to run through the classes of pi.
*/
PartitionIterator::PartitionIterator(const Partition& pi)
  :d_data(pi.size())
{
  d_stop.push_back(d_data.begin());

  if (d_data.size() == 0) // partition of empty set
    goto done;            // we must avoid taking |d_data.front()| below

  for (unsigned long j = 0; j < d_data.size(); ++j)
    d_data[j] = j;

  std::stable_sort(d_data.begin(),d_data.end(),
		   comparison::compare(pi.comparer()));

  {
    SubIterator data_end = d_data.end();
    unsigned long thisClass = pi.class_of(d_data.front());

    for (SubIterator i = d_data.begin(); i != data_end; ++i)
      if (pi.class_of(*i) != thisClass) // then start a new class
      {
	d_stop.push_back(i);
	thisClass = pi.class_of(*i);
      }
  }

  d_stop.push_back(d_data.end());

 done:
  d_stop.push_back(d_data.end());
  rewind();
}

PartitionIterator& PartitionIterator::operator++ ()

{
  ++d_currentEnd;
  d_range.first = d_range.second;
  d_range.second = *d_currentEnd;

  return *this;
}


} // |namespace partition|

} // |namespace atlas|
