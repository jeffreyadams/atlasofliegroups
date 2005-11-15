/*
  This is partition.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.3 

  See file main.cpp for full copyright notice
*/

#include "partition.h"

#include <algorithm>
#include <map>

#include "comparison.h"

/*****************************************************************************

   ... describe here when it is stable ....

******************************************************************************/

/*****************************************************************************

        Chapter I -- The Partition class.

******************************************************************************/

namespace atlas {

namespace partition {

/******** constructors and destructors ***************************************/

Partition::Partition(std::vector<unsigned long>& f)
  :d_class(f.size())

/*
  Synopsis: constructs a partition from the class vector.

  Explanation: the partition is given by the values of f.
*/

{
  std::map<unsigned long,unsigned long> val;

  unsigned long c = 0;

  for (size_t j = 0; j < f.size(); ++j)
    if (val.insert(std::make_pair(f[j],c)).second) { // found a new value
      newClass(j);
      ++c;
    } else {
      d_class[j] = val.find(f[j])->second;
    }
}

Partition::Partition(std::vector<unsigned long>& f, tags::UnnormalizedTag)
  :d_class(f)

/*
  Synopsis: like the previous one, but uses the actual values of f to number
  the classes.

  NOTE: it is recommended that the range of f be of the form [0,a[; I'm not
  sure what might break down if it isn't.
*/

{
  // make class representatives
  std::map<unsigned long,unsigned long> val;

  for (size_t j = 0; j < d_class.size(); ++j)
    val.insert(std::make_pair(d_class[j],j));

  d_classRep.resize(val.size());

  std::map<unsigned long,unsigned long>::iterator val_end = val.end();

  for (std::map<unsigned long,unsigned long>::iterator i = val.begin(); 
       i != val_end; ++i)
    d_classRep[i->first] = i->second;
}

/******** copy and assignment ************************************************/

void Partition::swap(Partition& other)

{
  d_class.swap(other.d_class);
  d_classRep.swap(other.d_classRep);

  return;
}

/******** manipulators *******************************************************/

void Partition::newClass(unsigned long j)

/*
  Synopsis: starts a new class at location j.
*/

{
  d_class[j] = d_classRep.size();
  d_classRep.push_back(j);

  return;
}

unsigned long Partition::classSize(unsigned long c) const

/*
  Synopsis: counts the number of elements in class #c.

  NOTE: Straightforward implementation.
*/

{
  unsigned long n = 0;

  for (size_t j = 0; j < d_class.size(); ++j)
    if (d_class[j] == c)
      ++n;

  return n;
}

}

/*****************************************************************************

        Chapter I -- The PartitionIterator class.

  The idea is that a PartitionIterator gives a convenient way of traversing
  the classes of a partition. It is intended that PartitionIterators are
  compared only if they refer to the same partition.

  The vector d_data contains the integers [0,n[, where n is the size of the
  set being partitioned, sorted in the order of the partition values. The
  vector d_stop contains an iterator pointing to the beginning of each
  class, and a final iterator pointing after the end of d_data. We only
  need then to return two consecutive elements in d_stop to bound a class.

******************************************************************************/

namespace partition {

PartitionIterator::PartitionIterator(const Partition& pi)
  :d_data(pi.size())

/*
  Initializes a PartitionIterator ready to run through the classes of pi.
*/

{
  using namespace comparison;

  d_stop.push_back(d_data.begin());

  if (d_data.size() == 0)
    goto done;

  for (unsigned long j = 0; j < d_data.size(); ++j)
    d_data[j] = j;

  std::stable_sort(d_data.begin(),d_data.end(),Compare<Partition>(pi));

  {
    std::vector<unsigned long>::const_iterator data_end = d_data.end();
    unsigned long thisClass = pi(d_data.front());
    
    for (std::vector<unsigned long>::const_iterator i = d_data.begin(); 
	 i != data_end; ++i)
      if (pi(*i) != thisClass) { // start a new class
	d_stop.push_back(i);
	thisClass = pi(*i);
      }
  }
       
  d_stop.push_back(d_data.end());

 done:
  d_stop.push_back(d_data.end());
  d_range.first = d_stop[0];
  d_range.second = d_stop[1];
  d_currentEnd = d_stop.begin()+1;
}

PartitionIterator& PartitionIterator::operator++ ()

{
  ++d_currentEnd;
  d_range.first = d_range.second;
  d_range.second = *d_currentEnd;

  return *this;
}

}

}
