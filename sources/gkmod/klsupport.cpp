/*!
\file
\brief Implementation for KLSupport.

  This module provides support code for the Kazhdan-Lusztig computation,
  mostly the management of the list of extremal pairs, and extremalization
  of arbitrary subsets of the block.

  After some hesitation, I think I _am_ going to assume that the block has
  been sorted by length and root datum involution, and then by R-packet.
  This does make the hard case in the recursion a lot simpler to handle:
  "thickets" of representations are in fact R-packets, so they will be
  consecutively numbered.
*/
/*
  This is klsupport.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4

  See file main.cpp for full copyright notice
*/

#include "klsupport.h"

#include "blocks.h"
#include "descents.h"

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
   d_rank(b.rank()),
   d_size(b.size())

{}

/******** copy, assignment and swap ******************************************/
void KLSupport::swap(KLSupport& other)

{
  d_state.swap(other.d_state);

  std::swap(d_block,other.d_block);
  std::swap(d_rank,other.d_rank);
  std::swap(d_size,other.d_size);

  d_primitivize.swap(other.d_primitivize);
  d_descent.swap(other.d_descent);
  d_goodAscent.swap(other.d_goodAscent);
  d_downset.swap(other.d_downset);
  d_primset.swap(other.d_primset);
  d_lengthLess.swap(other.d_lengthLess);

  return;
}

/******** accessors **********************************************************/
void KLSupport::extremalize(bitmap::BitMap& b, const bitset::RankFlags& d)
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
  for (size_t s = 0; s < d_rank; ++s)
    if (d.test(s))
      b &= d_downset[s];

  return;
}

void KLSupport::primitivize(bitmap::BitMap& b, const bitset::RankFlags& d)
  const

/*
  Synopsis: primitivizes b w.r.t. the values in d.

  Preconditions: the capacity of b is the size(); d contains d_rank valid
  flags;

  Explanation: an element z in the block is primitive w.r.t. d, if all the
  descents in d are either descents, or imaginary type II ascents for z. Since
  d_primset[s] flags the elements for which s is a descent or imaginary type
  II, this amounts to requesting that z belong to the intersection of the
  primsets for the various descents in d.
*/

{
  for (size_t s = 0; s < d_rank; ++s)
    if (d.test(s))
      b &= d_primset[s];

  return;
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

  //fill the primitivize table
  fillPrimitivize();
  d_state.set(Filled);

  return;
}

void KLSupport::fillDownsets()

/*
  Synopsis: fills in the downsets bitmaps, the primset bitmaps, the descent
  bitsets and the goodAscent bitsets.

  Explanation: here descent is taken in the weak sense of s belonging to the
  "tau-invariant" of z in b. In other words, s is a complex descent, real
  noncompact (type I or type II), or imaginary compact. The primset bitset
  records the x'es for which s is either a descent, or an imaginary type II
  ascent. The goodAscent bitsets hold the ascents that are not imaginary type
  II.

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
  std::vector<BitMap> ps(d_rank);
  std::vector<RankFlags> descents(size);
  std::vector<RankFlags> ga(size);

  for (size_t s = 0; s < ds.size(); ++s) {
    ds[s].resize(size);
    ps[s].resize(size);
    for (size_t z = 0; z < size; ++z) {
      DescentStatus::Value v = descentValue(s,z);
      if (DescentStatus::isDescent(v)) {
	ds[s].insert(z);
	ps[s].insert(z);
	descents[z].set(s);
      } else {
	if (v == DescentStatus::ImaginaryTypeII)
	  ps[s].insert(z);
	else
	  ga[z].set(s);
      }
    }
  }

  // commit
  d_downset.swap(ds);
  d_primset.swap(ps);
  d_descent.swap(descents);
  d_goodAscent.swap(ga);
  d_state.set(DownsetsFilled);

  return;
}

void KLSupport::fillPrimitivize()

/*
  Synopsis: fills in the table d_primitivize.

  Explanation: for each z in the block, and each set A of simple
  roots, d_primitivize[A][z] comes by proper ascents of z through s
  in A (either cross actions increasing length, or Cayley transforms)
  or 0 (if we come to a real non-parity case).

  Sets the PrimitivizeFilled bit in d_state if successful.
*/

{
  using namespace bitset;
  using namespace descents;
  using namespace blocks;

  if (d_state.test(PrimitivizeFilled))
    return;
  d_primitivize.reserve(1ul << rank());
  d_primitivize.resize(1ul << rank());
  size_t blocksize = d_block->size();

for (unsigned long j = 0 ; j >> rank() == 0 ; ++j) {
  BlockEltList prim;
  prim.reserve(blocksize);
  prim.resize(blocksize);
   const bitset::RankFlags A(j);
   for (BlockElt z = blocksize; z != 0;) {
     --z;
  // primitivize
     RankFlags a = goodAscentSet(z);
     a &= A;
     if (a.none()) {
       prim[z] = z;
       continue;
     }
    size_t s = a.firstBit();
    DescentStatus::Value v = descentValue(s,z);
    if (v == DescentStatus::RealNonparity) {
      prim[z] = UndefBlock;
      continue;
    }
    if (v == DescentStatus::ComplexAscent)
      prim[z] = prim[d_block->cross(s,z)];
    else
      prim[z] = prim[d_block->cayley(s,z).first];
  }
   // insert
   d_primitivize[j] = prim;
 }

 d_state.set(PrimitivizeFilled);

 return;
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
