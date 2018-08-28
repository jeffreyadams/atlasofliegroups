/*
  This is klsupport.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006-2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/* Implementation for KLSupport.

  This module provides support code for the Kazhdan-Lusztig computation,
  mostly the management of the list of primitive pairs, and primitivization
  of arbitrary subsets of the block.
*/

#include "klsupport.h"

/*
  After some hesitation, I think I _am_ going to assume that the block has
  been sorted by length and root datum involution, and then by R-packet.
  This does make the hard case in the recursion a lot simpler to handle:
  "thickets" of representations are in fact R-packets, so they will be
  consecutively numbered. [Fokko]
*/

namespace atlas {

namespace {

  void fillLengthLess
    (std::vector<BlockElt>&, const Block_base&);

} // |namespace|

/*****************************************************************************

        Chapter I -- The KLSupport class

 *****************************************************************************/

namespace klsupport {

KLSupport::KLSupport(const Block_base& b)
  : d_state()
  , d_block(b)
  , d_descent()
  , d_goodAscent()
  , d_downset()
  , d_primset()
  , d_lengthLess()
  , d_prim_index()
{}

/******** copy, assignment and swap ******************************************/

void KLSupport::swap(KLSupport& other)
{
  d_state.swap(other.d_state);

  assert(&d_block==&other.d_block);

  d_descent.swap(other.d_descent);
  d_goodAscent.swap(other.d_goodAscent);
  d_downset.swap(other.d_downset);
  d_primset.swap(other.d_primset);
  d_lengthLess.swap(other.d_lengthLess);
  d_prim_index.swap(other.d_prim_index);

}

/******** accessors **********************************************************/

/*
  Flag in b those block elements which are extremal w.r.t. the simple
  reflections in d.

  Preconditions: the capacity of b is the size(); d contains rank() flags;

  Explanation: an element z in the block is extremal w.r.t. d, if all the
  simple generators flagged in d are descents for z. Since d_downset[s] flags
  the elements for which s is a descent, this amounts to requesting that z
  belong to the intersection of the downsets for the various descents in d.
*/

void KLSupport::filter_extremal(BitMap& b, const RankFlags& d)
  const
{
  for (weyl::Generator s=0; s<rank(); ++s)
    if (d.test(s))
      b &= d_downset[s];
}


/*
  Primitivize b w.r.t. the values whose descent set is d, i.e.,
  throw out from b the non-primitive elements with respect to d

  Preconditions: the capacity of b is the size(); d contains rank() flags;

  Explanation: an element z in the block is primitive w.r.t. d, if all the
  simple generators flagged in d are either descents, or imaginary type II
  ascents for z. Since d_primset[s] flags the elements for which s is a
  descent or imaginary type II, this amounts to requesting that z belong to
  the intersection of the primsets for the various descents in d.
*/
void KLSupport::filter_primitive(BitMap& b, const RankFlags& d) const
{
  for (weyl::Generator s=0; s<rank(); ++s)
    if (d.test(s))
      b &= d_primset[s];
}


/*
  Finds for |x| a primitive element for |d| above it, returning that value, or
  returns |UndefBlock| if a real nonparity case is hit, or (in partial blocks)
  ascent through an undefined complex ascent or Cayley transform is attempted

  Explanation: a primitive element for |d| is one for which all elements in
  |d| are either descents or type II imaginary ascents. So if |x| is not
  primitive, it has an ascent in |d| that is either complex, imaginary type I
  or real nonparity. In the first two cases we replace |x| by the (unique)
  ascended element and continue; in the last case, we return |UndefBlock| (for
  K-L computations, this case implies that $P_{x,y}=0$; the value of
  |UndefBlock| is conveniently larger than any valid BlockElt |y|, so this
  case will be handled effortlessly together with triangularity). It is also
  permissible to pass |x==UndefBlock|, which will be returned immediately.
*/

#if 0 // code disabled because replaced by table look-up
BlockElt
  KLSupport::primitivize(BlockElt x, const RankFlags& d) const
{
  RankFlags a; // good ascents for x that are descents for y

  while (x!=UndefBlock and (a = goodAscentSet(x)&d).any())
  {
    size_t s = a.firstBit();
    DescentStatus::Value v = descentValue(s,x);
    x = v == DescentStatus::RealNonparity ? UndefBlock
      : v == DescentStatus::ComplexAscent ? d_block.cross(s,x)
      : d_block.cayley(s,x).first; // imaginary type I
  }
  return x;
  }
#endif

/******** manipulators *******************************************************/

// Fill the |lengthLess| table, and the downsets
void KLSupport::fill()
{
  if (d_state.test(Filled)) // do nothing
    return;

  // fill lengthLess if necessary
  if (not d_state.test(LengthLessFilled))
  {
    fillLengthLess(d_lengthLess,block());
    d_state.set(LengthLessFilled);
  }

  // make the downsets
  fillDownsets();

  //fill the primitivize table
  fillPrimitivize();

  d_state.set(Filled);

}


/*
  Fill in the |downset|, |primset|, |descents| and |goodAscent| bitmap/set
  vectors. Here |downset| and |primset| are vectors indexed by a simple
  reflection |s|, and giving a bitmap over all block elements, while
  |descents| and |goodAscent| are vectors indexed by a block element |z| and
  giving a bitset over all simple reflections. This difference is motivated by
  their use: |downset| and |primset| are used to filter bitmaps over the
  entire block according to some set of simple generators, which is easier if
  the data is grouped by generator. In fact the data computed is stored twice:
  one always has |downset[s].isMember(z) == descents[z].test(s)| and
  |primset[s].isMember(z) != good_ascent[z].test(s)|

  The predicate that |s| is a |descent| for |z| is taken in the weak sense
  that |s| belongs to the "tau-invariant" of |z|, in other words, it is a
  complex descent, real parity (type I or type II), or imaginary compact (the
  final case does not actually allow going down). The |goodAscent| bitset for
  |z| holds the non-decents for |z| that are not imaginary type II, so they
  are either complex ascent, imaginary type I or real nonparity. The |primset|
  bitmap for |s| records the block elements |z| for which |s| is not a
  |goodAscent|, in other words it is either a |descent| or imaginary type II.
*/
void KLSupport::fillDownsets()
{
  if (d_state.test(DownsetsFilled))
    return;

  size_t size = d_block.size();
  std::vector<BitMap> downset(rank());
  std::vector<BitMap> primset(rank());
  std::vector<RankFlags> descents(size);
  std::vector<RankFlags> good_ascent(size);

  for (weyl::Generator s=0; s<rank(); ++s)
  {
    downset[s].set_capacity(size);
    primset[s].set_capacity(size);
    for (BlockElt z = 0; z < size; ++z)
    {
      DescentStatus::Value v = descentValue(s,z);
      if (DescentStatus::isDescent(v))
      {
	downset[s].insert(z);
	primset[s].insert(z);
	descents[z].set(s);
      }
      else // ascents
	if (v == DescentStatus::ImaginaryTypeII)
	  primset[s].insert(z);  // s is a "bad" ascent
	else
	  good_ascent[z].set(s); // good ascent
    }
  }

  // commit
  d_downset.swap(downset);
  d_primset.swap(primset);
  d_descent.swap(descents);
  d_goodAscent.swap(good_ascent);
  d_state.set(DownsetsFilled);

}

/*
  Synopsis: fills in the table d_prim_index.

  Explanation: for each z in the block, and each descent set |A|, one can
  reach from |z| a primitive element |z'| for |A| by following 'C+' or 'i1'
  ascents, or possibly run into an 'rn' ascent aborting the search (when this
  happens any KL polynomial $P_{z,y}$ with |descentSet(y)==A| will be zero).

  Rather than |z'|, we store the index of |z'| in the list of primitives for
  |A|, this avoiding any form of binary search (but we will have to store
  null polynomials to ensure every polynomial has a predictable place).

  Sets the PrimitivizeFilled bit in d_state if successful.
*/

  void KLSupport::fillPrimitivize()
  {
    using namespace bitset;
    using namespace descents;
    using namespace blocks;

    if (d_state.test(PrimitivizeFilled))
      return;

    size_t two_to_the_r = 1ul << rank();
    d_prim_index.resize(two_to_the_r);

    size_t blocksize = d_block.size();

    for (unsigned long j = 0 ; j<two_to_the_r ; ++j)
    {
      const RankFlags A(j); // the current descent set

      BlockEltList prim;
      std::vector<unsigned int>& prim_index = d_prim_index[j];
      prim_index.resize(blocksize);

      BitMap primitives(size()); primitives.fill();
      filter_primitive(primitives,A); // compute all primitive elements for A

      unsigned int prim_count = primitives.size(); // start at high end

      for (BlockElt z = blocksize; z-->0;)
      {
	// primitivize
	RankFlags a = goodAscentSet(z)&A;
	if (a.none())
	{
	  assert(primitives.isMember(z)); // sanity check
	  prim_index[z] = --prim_count;
	  continue;
	}

        // if \emph{some} ascent is real nonparity, cop out right away
	for (RankFlags::iterator it=a.begin(); it(); ++it)
	  if (descentValue(*it,z) == DescentStatus::RealNonparity)
	    goto finish; // no primitivization

	{ // grouping needed because of the above goto
	  size_t s = a.firstBit();
	  DescentStatus::Value v = descentValue(s,z);
	  assert(v == DescentStatus::ComplexAscent or
		 v==DescentStatus::ImaginaryTypeI);
	  BlockElt sz =  v==DescentStatus::ComplexAscent
	    ? d_block.cross(s,z) : d_block.cayley(s,z).first;
	  if (sz==UndefBlock)
	    goto finish; // link out of partial block: give up primitivization
	  prim_index[z] = prim_index[sz];
	  continue; // avoid executing the |finish| code
	}
      finish: prim_index[z] = ~0; // mark invalid primtivization
      } // |for(z-->0)|
      assert(prim_count==0); // we've seen all primitives
    } // |for(A)|

    d_state.set(PrimitivizeFilled);

    return;
  } // |fillPrimitivize|
} // |namespace klsupport|

/*****************************************************************************

        Chapter II -- Functions local to this module

 *****************************************************************************/

namespace {

/*
  Make ll into a vector of size max(lengths(b))+2 such that for
  0<=l<=max(lengths(b))+1, the element ll[l] is the first BlockElt of length
  at least l in b (or b.size() if there are none, as for l=max(lengths(b))+1)
  In other words, as a number ll[l] counts the elements in b of length < l.

  Precondition: b is sorted by length;
*/
void fillLengthLess
  (std::vector<BlockElt>& ll, const Block_base& b)
{
  ll.clear(); ll.resize(b.size()==0 ? 2 : b.length(b.size()-1)+2);

  ll[0]=0; // no elements of length<0
  size_t l=0;
  for (BlockElt z=0; z<b.size(); ++z)
    while (b.length(z)>l)
      ll[++l]=z;  // invariant |l<=b.length(z)| after this statement

  // do not forget the last length!
  ll[l+1]=b.size(); // here l=b.length(b.size()-1)=max(lengths(b))
}

} // |namsespace|

} // |namespace atlas|
