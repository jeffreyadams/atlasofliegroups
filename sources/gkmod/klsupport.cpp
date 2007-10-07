/*!
\file
\brief Implementation for KLSupport.

  This module provides support code for the Kazhdan-Lusztig computation,
  mostly the management of the list of primitive pairs, and primitivization
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
  part of the Atlas of Reductive Lie Groups

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
  This does make the hard case in the recursion a lot simpler to handle:
  "thickets" of representations are in fact R-packets, so they will be
  consecutively numbered. [Fokko]
*/

namespace atlas {

namespace {
  using blocks::BlockElt;

  void pause() {;}

  void fillLengthLess(std::vector<BlockElt>&, const blocks::Block&);

} // namespace

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

  d_descent.swap(other.d_descent);
  d_goodAscent.swap(other.d_goodAscent);
  d_downset.swap(other.d_downset);
  d_primset.swap(other.d_primset);
  d_lengthLess.swap(other.d_lengthLess);

}

/******** accessors **********************************************************/
void KLSupport::extremalize(bitmap::BitMap& b, const bitset::RankFlags& d)
  const

/*!
  \brief: Flags in b those block elements which are extremal w.r.t. the
  simple reflections in d.

  Preconditions: the capacity of b is the size(); d contains d_rank flags;

  Explanation: an element z in the block is extremal w.r.t. d, if all the
  descents in d are also descents for z. Since d_downset[s] flags the
  elements for which s is a descent, this amounts to requesting that z
  belong to the intersection of the downsets for the various descents in d.
*/

{
  for (size_t s = 0; s < d_rank; ++s)
    if (d.test(s))
      b &= d_downset[s];

}

void KLSupport::primitivize(bitmap::BitMap& b, const bitset::RankFlags& d)
  const

/*
  \brief Primitivizes b w.r.t. the values whose descent set is d, i.e.,
  throws out from b the non-primitive elements with respect to d

  Preconditions: the capacity of b is the size(); d contains d_rank flags;

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

}

BlockElt KLSupport::primitivize(BlockElt x, const bitset::RankFlags& d) const

/*!
  \brief Replaces x (the number of a block element) with a primitive
  element above it, and returns true, or returns false, and x is not
  changed.

  Explanation: a primitive element for d is one for which all elements in d
  are descents or type II imaginary ascents. So if x is not primitive, it has
  an ascent that is either complex, imaginary type I or real compact. In the
  first two cases we replace x by the ascended element and continue; in the
  last case, we return UndefBlock (for K-L computations, this will imply that
  $P_{x,y}=0$; the value of UndefBlock is conveniently larger than any valid
  BlockElt y, so this condition can be handled together with triangularity).
*/

{
  using namespace descents;

  bitset::RankFlags a; // good ascents for x that are descents for y

  while ((a = goodAscentSet(x)&d).any()) {
    size_t s = a.firstBit();
    DescentStatus::Value v = descentValue(s,x);
    if (v == DescentStatus::RealNonparity)
      return blocks::UndefBlock; // cop out
    x = v == DescentStatus::ComplexAscent // complex or imaginary type I ?
	? d_block->cross(s,x)
	: d_block->cayley(s,x).first;
  }
  return x;
}

/******** manipulators *******************************************************/
void KLSupport::fill()

/*
  \brief Fills the lengthLess table, and the downsets
*/

{
  using namespace bitmap;

  if (d_state.test(Filled)) // do nothing
    return;

  // fill lengthLess if necessary
  if (not d_state.test(LengthLessFilled)) {
    fillLengthLess(d_lengthLess,block());
    d_state.set(LengthLessFilled);
  }

  // make the downsets
  fillDownsets();

  d_state.set(Filled);

}

void KLSupport::fillDownsets()

/*
  \brief Fills in the downset, primset, descents and goodAscent bitmap/set
  vectors. Here downset and primset are vectors indexed by a simple reflection
  s, and giving a bitmap over all block elements, while descents and
  goodAscent are vectors indexed by a block element z and giving a bitset over
  all simple reflections.

  Explanation: here the predicate that s is a descent for z is taken in the
  weak sense that s belongs to the "tau-invariant" of z in b, in other
  words, it is a complex descent, real noncompact (type I or type II), or
  imaginary compact. The goodAscent bitset for z holds the ascents for z that
  are not imaginary type II. The primset bitmap for s records the block
  elements z for which s is not a goodAscent, in other words it is either
  a descent, or an imaginary type I ascent.

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
    ds[s].set_capacity(size);
    ps[s].set_capacity(size);
    for (BlockElt z = 0; z < size; ++z) {
      DescentStatus::Value v = descentValue(s,z);
      if (DescentStatus::isDescent(v)) {
	ds[s].insert(z);
	ps[s].insert(z);
	descents[z].set(s);
      } else { // ascents
	if (v == DescentStatus::ImaginaryTypeII)
	  ps[s].insert(z);  // s is a "bad" ascent
	else
	  ga[z].set(s); // good ascent
      }
    }
  }

  // commit
  d_downset.swap(ds);
  d_primset.swap(ps);
  d_descent.swap(descents);
  d_goodAscent.swap(ga);
  d_state.set(DownsetsFilled);

}

} // namespace klsupport

/*****************************************************************************

        Chapter II -- Functions local to this module

 *****************************************************************************/

namespace {

void fillLengthLess(std::vector<BlockElt>& ll, const blocks::Block& b)

/*
  \brief Makes ll into a vector of size max(lengths(b))+2 such that for
  0<=l<=max(lengths(b))+1, the element ll[l] is the first BlockElt of length
  at least l in b (or b.size() if there are none, as for l=max(lengths(b))+1)
  In other words, as a number ll[l] counts the elements in b of length < l.

  Precondition: b is sorted by length;
*/

{
  ll.clear(); ll.resize(b.length(b.size()-1)+2);

  ll[0]=0; // no elements of length<0
  size_t l=0;
  for (blocks::BlockElt z=0; z<b.size(); ++z)
    while (b.length(z)>l)
      ll[++l]=z;  // invariant |l<=b.length(z)| after this statement

  // do not forget the last length!
  ll[l+1]=b.size(); // here l=b.length(b.size()-1)=max(lengths(b))
}

}

}
