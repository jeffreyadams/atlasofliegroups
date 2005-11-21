/*
  This is klsupport.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "klsupport.h"

#include "blocks.h"
#include "descents.h"

/*
  This module provides support code for the Kazhdan-Lusztig computation,
  mostly the management of the list of extremal pairs, and extremalization
  of arbitrary subsets of the block.

  After some hesitation, I think I _am_ going to assume that the block has
  been sorted by length and root datum involution, and then by R-packet.
  This does make the hard case in the recurdion a lot simpler to handle:
  "thickets" of representations are in fact R-packets, so they will be
  consecutively numbered.
*/

namespace atlas {

namespace {

  void pause() {;}

  void fillLengthLess(std::vector<size_t>&, const blocks::Block&);

}

/*****************************************************************************

        Chapter I -- The KLSupport class

  ... explain here when it is stable ...

 *****************************************************************************/

namespace klsupport {

KLSupport::KLSupport(blocks::Block& b)
  :d_block(&b),
   d_rank(b.rank())

{}

/******** copy, assignment and swap ******************************************/
void KLSupport::swap(KLSupport& other)

{
  d_state.swap(other.d_state);

  std::swap(d_block,other.d_block);
  std::swap(d_rank,other.d_rank);

  d_extrPairs.swap(other.d_extrPairs);
  d_descent.swap(other.d_descent);
  d_downset.swap(other.d_downset);
  d_lengthLess.swap(other.d_lengthLess);

  return;
}

/******** accessors **********************************************************/
descents::DescentStatus::Value KLSupport::descentValue(size_t s, size_t z) 
  const

/*
  Synopsis: returns the descent status of z w.r.t. s.

  NOTE: this is not inlined in order to avoid a compilation dependency on
  blocks.h
*/

{
  return d_block->descentValue(s,z);
}

void KLSupport::extremalize(bitmap::BitMap& b, 
			    const descents::DescentStatus& d) 
  const

/*
  Synopsis: extremalizes b w.r.t. the values in d.

  Preconditions: the capacity of b is the size(); d contains d_rank valid
  flags;

  Explanation: an element z in the block is extremal w.r.t. d, if all the
  descents in d are also descents for z. Since d_downset[s] flags the
  elements for which s is a descent, this amounts to requesting that z
  belong to the intersection of the downsets for the various descents in d.
*/

{
  using namespace descents;

  for (size_t s = 0; s < d_rank; ++s)
    if (DescentStatus::isDescent(d[s]))
      b &= d_downset[s];

  return;
}

size_t KLSupport::length(size_t z) const

/*
  Synopsis: returns the length of z.

  NOTE: this is not inlined in order to avoid a compilation dependency on
  blocks.h
*/

{
  return d_block->length(z);
}

/******** manipulators *******************************************************/
void KLSupport::fill()

/*
  Synopsis: fills the extrPairs list, the downsets, and the lengthLess
  table.

  Explanation: for real reductive Lie groups, we always compute the full
  block; niceties like partial computations are not expected to be very useful
  nor needed in the cases of interest.
*/

{
  using namespace bitmap;

  if (d_state.test(Filled)) // do nothing
    return;

  // fill lengthLess if necessary
  if (d_state.test(LengthLessFilled) == false) {
    fillLengthLess(d_lengthLess,block());
    d_state.set(LengthLessFilled);
  }

  // make the downsets
  fillDownsets();

  // fill in extrPairs
  d_extrPairs.resize(block().size());

  for (size_t z = 0; z < d_extrPairs.size(); ++z) {
    BitMap b(block().size());
    size_t c = d_lengthLess[length(z)];
    b.fill(c);
    b.insert(z);
    // extremalize
    extremalize(b,block().descent(z));
    // copy to list
    std::copy(b.begin(),b.end(),back_inserter(d_extrPairs[z]));
  }

  d_state.set(Filled);

  return;
}

void KLSupport::fillDownsets()

/*
  Synopsis: fills in the downsets bitmaps, and the descent bitsets.

  Explanation: here descent is taken in the weak sense of s belonging to the
  "tau-invariant" of z in b. In other words, s is a complex descent, real
  noncompact (type I or type II), or imaginary compact.

  Sets the DownsetsFilled bit in d_state if successful. Commit-or-rollback
  is guaranteed.
*/

{  
  using namespace bitmap;
  using namespace bitset;
  using namespace descents;

  if (d_state.test(DownsetsFilled))
    return;

  size_t size = d_block->size();
  std::vector<BitMap> ds(d_rank);
  std::vector<RankFlags> descents(size);

  for (size_t s = 0; s < ds.size(); ++s) {
    BitMap& b = ds[s];
    b.resize(size);
    for (size_t z = 0; z < size; ++z)
      if (DescentStatus::isDescent(descentValue(s,z))) {
	b.insert(z);
	descents[z].set(s);
      }
  }

  // commit
  d_downset.swap(ds);
  d_descent.swap(descents);
  d_state.set(DownsetsFilled);

  return;
}

size_t KLSupport::numExtremals()

/*
  Synopsis: returns the total number of extremal pairs.

  The important thing here is that we do it without constructing the extremal
  list. So this function is relatively slow, but does not require memory to
  speak off.
*/

{
  using namespace bitmap;

  if (d_state.test(DownsetsFilled) == false)
    fillDownsets();

  // fill lengthLess if necessary
  if (d_state.test(LengthLessFilled) == false) {
    fillLengthLess(d_lengthLess,block());
    d_state.set(LengthLessFilled);
  }

  size_t size = d_block->size();
  size_t count = 0;

  for (size_t z = 0; z < size; ++z) {
    BitMap b(size);
    b.fill(d_lengthLess[length(z)]);
    b.insert(z);
    // extremalize
    extremalize(b,block().descent(z));
    // copy to list
    count += b.size();
  }

  return count;
}

}

/*****************************************************************************

        Chapter II -- Functions local to this module

  ... explain here when it is stable ...

 *****************************************************************************/

namespace {

void fillLengthLess(std::vector<size_t>& ll, const blocks::Block& b)

/*
  Synopsis: puts in ll[l] the number of elements in b of length < l.

  Precondition: b is sorted by length;
*/

{
  ll.clear();

  size_t l = b.length(0);

  for (size_t j = 0; j <= l; ++j)
    ll.push_back(0);

  for (size_t z = 0; z < b.size(); ++z)
    if (b.length(z) > l) { // new length
      ll.push_back(z);
      l = b.length(z);
    }

  // do not forget the last length!
  ll.push_back(b.size());

  return;
}

}

}
