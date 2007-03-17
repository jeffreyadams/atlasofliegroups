/*!
\file
  This is tally.h
*/
/*
  Copyright (C) 2007 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

/* The purpose of this module it to provide a simple class template |TallyVec|
   that can be used to maintain a histogram in an efficient and flexible way.
   The class behaves as a vector of counters, indexed by unsigned numbers, and
   the basic operation provided is to increment the counter for a given index.
   The resulting vector can be written to are read back from a file, and an
   operation is provided to make a histogram of the counter values themselves.
*/

#ifndef TALLY_H
#define TALLY_H

#include <vector>
#include <map>
#include <iostream>

namespace atlas {
namespace tally {

/* The |TallyVec| class template has one type parameters that should be an
   unsigned integral type: it is the type used for the primary counters,
   determining the range of frequency that can be maintained most efficiently
   (one a counter overflows it will be maintained in a less storage- and
   time-efficient way). Typically it will be a small type to avoid wasting
   space on the average low counter values. The type indexing the counters and
   used for the counter values when they overflow are considered to be of less
   impact, so both will be taken |unsigned long long int| in all cases.
*/

template <typename Count>
class TallyVec
{
 public:
  typedef unsigned long long int Index;
  typedef unsigned long long int ullong; // when used as counter value
 private:
  typedef std::map<Index,ullong> map_type;

  static const Count maxCount; // limit of |Count| type

  std::vector<Count> count; // primary histogram
  map_type overflow; // for multiplicities >=256
  Index max;         // maximal index seen
  ullong total;      // grans total

 public:
  TallyVec(size_t limit) : count(0), overflow(), max(0), total(0)
  { count.reserve(limit); }
  TallyVec (std::istream& file); // recover table dumped to file

  bool tally (Index i); // increase count for i by 1; tell whether new
  bool tally (Index i,ullong multiplicity); // same with multiplicity
  Index size() const { return max+1; } // size of collection now tallied
  inline ullong multiplicity (Index i) const;
  ullong sum() const { return total; }

  template<typename MuCount> TallyVec<MuCount> derived(size_t limit) const;

  void advance(Index& i) const; // like ++ and -- where |i| iterates
  void lower(Index& i) const;   // over indices with nonzero multiplicity

  void write_to(std::ostream& out) const; // dump table to file
}; // class TallyVec

} // namespace tally
} // namespace atlas

#include "tally_def.h"

#endif
