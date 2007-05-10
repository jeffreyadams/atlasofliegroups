/*!
\file
\brief Implementation of the class Block.

  This module aims to construct the fundamental combinatorial structure
  underlying k-l computations: a block of representations, together with
  cross-actions, Cayley transforms, length function, and descent sets.

  Basically, a block should be constructed from a pair of gradings: one
  for the fundamental torus of G, and one for the fundamental torus of
  G^vee. This suffices to determine the cocycles that allow us to describe
  cosets for the corresponding Tits groups, and therefore the one-sided
  parameter sets for which we then only have to take the restricted
  product. Of course in practice one should be careful to first determine
  the support of the block in the set of root datum involutions, and
  do the construction only for that part, if possible.
*/
/*
  This is blocks.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "blocks.h"

#ifdef VERBOSE
#include <iostream>
#endif

#include <cassert>
#include <memory>
#include <set>

#include "bruhat.h"
#include "complexredgp.h"
#include "descents.h"
#include "kgb.h"
#include "realredgp.h"
#include "tags.h"
#include "weyl.h"

/*
  This module aims to construct the fundamental combinatorial structure
  underlying k-l computations: a block of representations, together with
  cross-actions, Cayley transforms, length function, and descent sets.

  Basically, a block should be constructed from a pair of gradings: one
  for the fundamental torus of G, and one for the fundamental torus of
  G^vee. This suffices to determine the cocycles that allow us to describe
  cosets for the corresponding Tits groups, and therefore the one-sided
  parameter sets for which we then only have to take the restricted
  product. Of course in practice one should be careful to first determine
  the support of the block in the set of root datum involutions, and
  do the construction only for that part, if possible.
*/

namespace atlas {

  //namespace {
  namespace blocks {
    namespace helper {
  void pause() {;}

using namespace blocks;

void makeHasse(std::vector<set::SetEltList>&, const Block&);

void insertAscents(std::set<BlockElt>&, const set::SetEltList&, size_t,
		   const Block&);

      /*!
\brief Derived class of Block, to carry out the construction of Block.
      */
class Helper:public Block {

private:

  size_t d_size; // number of block elements

  kgb::KGB d_kgb;
  kgb::KGB d_dualkgb;

  // d_firstx maps KGB elements x to the first block element z with d_x[z]>=x
  // and moreover d_firstx[d_kgb.size()]==d_size
  std::vector<BlockElt> d_firstx; // of size d_kgb.size()+1

  weyl::WeylInterface d_toDualWeyl;
  const weyl::WeylGroup* d_dualWeylGroup;

public:
  // constructors and destructors
  Helper(realredgp::RealReductiveGroup& G, realredgp::RealReductiveGroup& dG);

  virtual ~Helper() {};

  // accessors
  weyl::TwistedInvolution dualInvolution(const weyl::TwistedInvolution&) const;

  const weyl::WeylGroup& dualWeylGroup() const {
    return *d_dualWeylGroup;
  }

  // manipulators
  void fillCayleyActions();

  void fillCrossActions();

  void fillDescents();

  void fillInvolutions();

  void fillInvolutionSupports();

  void fillLengths();

  void makeWeylCorrelation();

  void orbitPairs();

  void tauCorrelation();

};

    }
  }

/*****************************************************************************

        Chapter I -- The Block class

  ... explain here when it is stable ...

******************************************************************************/

namespace blocks {

  using namespace atlas::blocks::helper;

Block::Block():d_bruhat(0)

{}

Block::Block(complexredgp::ComplexReductiveGroup& G,
	     realform::RealForm rf, realform::RealForm df)
  :d_bruhat(0)

/*!
  \brief Constructor for the block class.

  Constructs a block from the datum of a real form rf for G and a real
  form df for G^vee (_not_ strong real forms: up to isomorphism, the
  result depends only on the underlying real forms!).

  This is a big job; the work is deferred to the Helper constructor.
*/

{
  using namespace complexredgp;
  using namespace realredgp;
  using namespace tags;

#ifdef VERBOSE
  std::cerr << "entering block construction ..." << std::endl;
#endif

  RealReductiveGroup G_R(G,rf);
  ComplexReductiveGroup dG(G,DualTag()); // the dual group
  RealReductiveGroup dG_R(dG,df);
  dG_R.fillCartan();

  Helper help(G_R,dG_R);
  swap(help);

#ifdef VERBOSE
  std::cerr << "done" << std::endl;
#endif
}

Block::~Block()

{
  delete d_bruhat;
}

/******** copy, assignment and swap ******************************************/
void Block::swap(Block& other)

{
  std::swap(d_rank,other.d_rank);
  std::swap(d_xsize,other.d_xsize);
  std::swap(d_ysize,other.d_ysize);

  d_x.swap(other.d_x);
  d_y.swap(other.d_y);
  d_cross.swap(other.d_cross);
  d_cayley.swap(other.d_cayley);
  d_inverseCayley.swap(other.d_inverseCayley);
  d_descent.swap(other.d_descent);
  d_length.swap(other.d_length);
  d_involution.swap(other.d_involution);
  d_involutionSupport.swap(other.d_involutionSupport);

  std::swap(d_realForm,other.d_realForm);
  std::swap(d_dualForm,other.d_dualForm);

  d_state.swap(other.d_state);

  std::swap(d_bruhat,other.d_bruhat);
  std::swap(d_weylGroup,other.d_weylGroup);

  return;
}

/******** accessors **********************************************************/
size_t Block::firstStrictDescent(BlockElt z) const

/*!
  \brief Returns the first descent for z (the number of a simple root) that is
not imaginary compact, or rank() if there is no such descent.
*/

{
  using namespace descents;

  for (size_t s = 0; s < rank(); ++s) {
    DescentStatus::Value v = descentValue(s,z);
    if (DescentStatus::isDescent(v) and v != DescentStatus::ImaginaryCompact) {
      return s;
    }
  }

  return rank();
}

bool Block::isStrictAscent(size_t s, BlockElt z) const

/*!
  \brief Tells if s is a strict ascent generator for z.

  Explanation: this means that d_descent[z][s] is one of ComplexAscent,
  ImaginaryTypeI or ImaginaryTypeII.
*/

{
  using namespace descents;

  DescentStatus::Value v = descentValue(s,z);
  return not DescentStatus::isDescent(v) and v!=DescentStatus::RealNonparity;
}

bool Block::isStrictDescent(size_t s, BlockElt z) const

/*!
  \brief Tells if s is a strict descent generator for z.

  Explanation: this means that d_descent[z][s] is one of ComplexDescent,
  RealTypeI or RealTypeII.
*/

{
  using namespace descents;

  DescentStatus::Value v = descentValue(s,z);
  return DescentStatus::isDescent(v) and v!=DescentStatus::ImaginaryCompact;
}

/*!
  \brief the functor $T_{\alpha,\beta}$

  Precondition: alpha and beta are adjacent roots, of which alpha is a descent
  set for y, while beta is not a sescent for y.

  In fact if this is not satisfied, we return a pair of UndefBlock elements
*/

BlockEltPair Block::link(size_t alpha,size_t beta,BlockElt y) const
{
  const descents::DescentStatus& d=descent(y);

  std::vector<BlockElt> result(2,UndefBlock);
  std::vector<BlockElt>::iterator it=result.begin();

  BlockEltPair p=inverseCayley(alpha,y);
  switch (d[alpha])
  {
  case descents::DescentStatus::ComplexDescent:
    {
      BlockElt y1=cross(alpha,y);
      if (descents::DescentStatus::isDescent(descent(y1)[beta]))
	*it++=y1;
      break;
    }
  case descents::DescentStatus::RealTypeI:
    if (descents::DescentStatus::isDescent(descent(p.second)[beta]))
      *it++=p.second;
    // FALL THROUGH
  case descents::DescentStatus::RealTypeII:
    if (descents::DescentStatus::isDescent(descent(p.first)[beta]))
      *it++=p.first;
    break;
  default: {}
  } // switch(d[alpha])

  p=cayley(beta,y);
  switch (d[beta])
  {
  case descents::DescentStatus::ComplexAscent:
    {
      BlockElt y1=cross(beta,y);
      if (not descents::DescentStatus::isDescent(descent(y1)[alpha]))
	*it++=y1;
      break;
    }
  case descents::DescentStatus::ImaginaryTypeII:
    if (not descents::DescentStatus::isDescent(descent(p.second)[alpha]))
      *it++=p.second;
    // FALL THROUGH
  case descents::DescentStatus::ImaginaryTypeI:
    if (not descents::DescentStatus::isDescent
	(descent(p.first)[alpha]))
      *it++=p.first;
    break;
  default: {}
  } // switch(d[beta])

  assert(&*it<=&result[2]);

  return std::make_pair(result[0],result[1]);
}


/******** manipulators *******************************************************/
void Block::fillBruhat()

/*!
  \brief Constructs the BruhatOrder.

  NOTE: may throw a MemoryOverflow error. Commit-or-rollback is guaranteed.
*/

{
  using namespace bruhat;
  using namespace set;

  if (d_state.test(BruhatConstructed)) // work was already done
    return;

  try {
    std::vector<SetEltList> hd; makeHasse(hd,*this);
    BruhatOrder* bp = new BruhatOrder(hd); // may throw here

    // commit
    delete d_bruhat; // this is overly careful: it must be NULL
    d_bruhat = bp;
    d_state.set(BruhatConstructed);
  }
  catch (std::bad_alloc) { // transform failed allocation into MemoryOverflow
    throw error::MemoryOverflow();
  }
}

}

/*****************************************************************************

        Chapter II -- The Helper class

  ... explain here when it is stable ...

******************************************************************************/

// namespace {
  namespace blocks {
  namespace helper {
Helper::Helper(realredgp::RealReductiveGroup& G,
	       realredgp::RealReductiveGroup& dG)
  :d_kgb(G),
   d_dualkgb(dG)

/*!
  \brief Constructor for the helper class.

  Facilitates the construction of a block through the use of some extra data
  (mostly the kgb orbit structure on both the group side and the dual side.)

  NOTE: simple-minded version. We don't aim for maximal efficiency here. This
  does way too much work for small blocks.
*/

{
  using namespace complexredgp;
  using namespace kgb;
  using namespace weyl;

  pause();

  d_weylGroup = &d_kgb.weylGroup();
  d_dualWeylGroup = &d_dualkgb.weylGroup();
  makeWeylCorrelation();

  d_xsize = d_kgb.size();
  d_ysize = d_dualkgb.size();

  d_realForm = G.realForm();
  d_dualForm = dG.realForm();

  ComplexReductiveGroup& G_C = G.complexGroup();
  d_rank = G_C.semisimpleRank();

  d_size = G_C.blockSize(d_realForm,d_dualForm); // cannot use size() yet!

  // set sizes of main lists
  d_cross.resize(d_rank);
  d_cayley.resize(d_rank);
  d_inverseCayley.resize(d_rank);

  // fill main lists with vectors of d_size default (undefined) values
  for (size_t s = 0; s < d_rank; ++s) {
    d_cross[s].assign(d_size,UndefBlock);
    d_cayley[s].assign(d_size,std::make_pair(UndefBlock,UndefBlock));
    d_inverseCayley[s].assign(d_size,std::make_pair(UndefBlock,UndefBlock));
  }

  d_descent.resize(d_size);
  d_length.assign(d_size,0);
  d_involution.assign(d_size,TwistedInvolution()); // makes size() defined


  d_x.reserve(d_size+1);
  d_y.reserve(d_size+1);
  d_firstx.reserve(d_kgb.size()+1);

  // make list of allowable orbit pairs (block elements) and fill x-offsets
  orbitPairs();

  // fill in cross actions
  fillCrossActions();

  // fill in cayley actions
  fillCayleyActions();

  // fill in involutions
  fillInvolutions();

  // fill in lengths
  fillLengths();

  // fill in descents
  fillDescents();

  // fill in involution supports
  fillInvolutionSupports();
}

/******** accessors **********************************************************/
weyl::TwistedInvolution Helper::dualInvolution
  (const weyl::TwistedInvolution& d_w) const

/*!
  \brief Returns the twisted involution dual to d_w.

  Explanation: we have tau = w.tau_f, with w in the Weyl group, and tau_f
  the fundamental involution of the character lattice. We seek v s.t.
  -{}^t tau = v.tau_f^{vee}, where tau_f^{vee} = -w_0{}^t tau_f. So we
  have:

        (-{}^t tau_f).{}^t w = v.w_0.(-{}^t tau_f)

  which leads to v = ({}^t theta_f(w)).w_0, where theta_f comes from the
  Dynkin diagram involution defining the inner class. The transposition
  takes a root reflection to the corresponding coroot reflection, and reverses
  the order (so in practice, it amounts to going over to the inverse of v)
  and so we get:

        (theta_f(w))^{-1}.w_0 = v.

  NOTE: one extra twist is that we need to move the representation of w from
  the Weyl group of d_kgb to that of d_dualkgb using the d_toDualWeyl table
  (see the comment for makeWeylCorrelation.)
*/

{
  using namespace weyl;

  const WeylGroup& W = weylGroup();

  WeylElt v = W.longest();
  WeylElt w = d_w.representative();
  W.translate(w,d_toDualWeyl);

  W.twist(w);
  W.prod(v,w);
  W.invert(v);

  return TwistedInvolution(v);
}

/******** manipulators *******************************************************/
void Helper::fillCayleyActions()

/*!
  \brief Fills in the cayley actions.

  Explanation: as for cross actions, this is obtained from kgb and dualkgb.
  The cayley transform of (x,y) is obtained by applying the direct Cayley
  transform to x, the inverse Cayley transform to y. The correlation between
  kgb and dualkgb is such that these operations are either both defined or
  both undefined. Just as for cross-actions, the main difficulty is in
  locating the result.
*/

{
  using namespace kgb;

  for (size_t s = 0; s < d_rank; ++s) {
    for (BlockElt z = 0; z < size(); ++z) {
      KGBElt sx = d_kgb.cayley(s,x(z));
      if (sx == UndefKGB) // cayley transform is undefined
	continue;
      KGBEltPair sy = d_dualkgb.inverseCayley(s,y(z));
      assert(sy.first != UndefKGB);
      KGBEltList::iterator i =
	std::find(d_y.begin()+d_firstx[sx],d_y.begin()+d_firstx[sx+1]
		 ,sy.first);
      assert(i  < d_y.begin()+d_firstx[sx+1]);
      BlockElt z0=i - d_y.begin();
      d_cayley[s][z].first = z0;
      {
	// fill in inverse Cayley entry
	BlockEltPair& ic = d_inverseCayley[s][z0];
	if (z < ic.first) {
	  ic.second = ic.first;
	  ic.first = z;
	} else {
	  ic.second = z;
	}
      }
      if (sy.second == UndefKGB) // cayley transform is type I
	continue;
      // if we get here the Cayley transform is two-valued
      i = std::find(d_y.begin()+d_firstx[sx],d_y.begin()+d_firstx[sx+1]
		   ,sy.second);
      assert(  i  < d_y.begin()+d_firstx[sx+1]);
      BlockElt z1=i - d_y.begin();
      d_cayley[s][z].second = z1;
      {
	// fill in inverse Cayley entry
	BlockEltPair& ic = d_inverseCayley[s][z1];
	assert(ic.first == UndefBlock);
	ic.first = z;
      }
    }
  }

  return;
}

void Helper::fillCrossActions()

/*!
  \brief Fills in the cross actions.

  Explanation: in principle this is simple: in terms of orbit pairs, it is
  just the corresponding action on each side. The main difficulty is in
  locating  the result. We do this with the aid off the d_firstx table:
  the result of cross(s,(x,y)) lies in the range beginning at
  firstx(cross(s,x)) and ending at firstx(cross(s,x)+1).

  Note that std::find does a linear (left-to-right) search, so that no other
  occurrence of sy will be found.
*/

{
  using namespace kgb;

  for (size_t s = 0; s < d_rank; ++s) {
    for (BlockElt z = 0; z < size(); ++z) {
      KGBElt sx = d_kgb.cross(s,x(z));
      KGBElt sy = d_dualkgb.cross(s,y(z));
      KGBEltList::iterator i =
	std::find(d_y.begin()+d_firstx[sx],d_y.begin()+d_firstx[sx+1]
		 ,sy);
      assert(i <  d_y.begin()+d_firstx[sx+1]);
      d_cross[s][z] = i - d_y.begin();
    }
  }

  return;
}

void Helper::fillDescents()

/*!
  \brief Fills in the descent table.

  Explanation: a representation can have eight different descent stata:
  descents can be imaginary compact, complex, real type I and real type II;
  ascents can be real non-parity (I think!), complex, imaginary noncompact
  type I and type II. The status of each element of the block w.r.t. each
  generator is easily deduced from the status and descent information for
  the two kgb pieces.
*/

{
  using namespace descents;
  using namespace gradings;
  using namespace kgb;
  using namespace weyl;

  for (BlockElt z = 0; z < d_size; ++z)
    for (size_t s = 0; s < d_rank; ++s) {
      KGBElt x = d_x[z];
      KGBElt y = d_y[z];
      if (d_kgb.status(s,x) == Status::Complex) { // s is complex
	if (d_kgb.isDescent(s,x))
	  d_descent[z].set(s,DescentStatus::ComplexDescent);
	else
	  d_descent[z].set(s,DescentStatus::ComplexAscent);
	continue;
      }
      // now s is  real or imaginary
      if (d_kgb.status(s,x) == Status::ImaginaryNoncompact) {
	if (cross(s,z) != z) // type I
	  d_descent[z].set(s,DescentStatus::ImaginaryTypeI);
	else // type II
	  d_descent[z].set(s,DescentStatus::ImaginaryTypeII);
	continue;
      }
      if (d_dualkgb.status(s,y) == Status::ImaginaryNoncompact) {
	if (cross(s,z) != z) // type II
	  d_descent[z].set(s,DescentStatus::RealTypeII);
	else // type I
	  d_descent[z].set(s,DescentStatus::RealTypeI);
	continue;
      }
      // now s is imaginary compact or real nonparity
      if (d_kgb.status(s,x) == Status::Real)
	d_descent[z].set(s,DescentStatus::RealNonparity);
      else
	d_descent[z].set(s,DescentStatus::ImaginaryCompact);
    }

  return;
}

void Helper::fillInvolutions()

/*!
  \brief Fills in the involution table.

  This is directly deduced from the x-part.
*/

{
  for (BlockElt z = 0; z < d_size; ++z)
    d_involution[z] = d_kgb.involution(d_x[z]);

  return;
}

void Helper::fillInvolutionSupports()

/*!
  \brief Fills in d_involutionSupport.

  Explanation: d_involutionSupport[z] coontains the support of the
  underlying involution, i.e. the generators s that occur either as
  s or as twist(s) in a reduced expression for w as an involution.

  Algorithm: for the minimal length elements, write down the support
  directly from a reduced expression. For the general case, find a
  strict descent s; then if s.z is the descended z, we have supp(z)
  = supp(sz) cup {s} if the descent is a cayley, supp(s.z) cup {s}
  cup{twist(s)} if the descent is a cross.
*/

{
  using namespace descents;
  using namespace weyl;

  const WeylGroup& W = weylGroup();
  d_involutionSupport.resize(size());

  size_t minLength = length(0);
  BlockElt z = 0;

  for (; z < size() and length(z) == minLength; ++z) {
    WeylWord ww;
    W.out(ww,involution(z).representative());
    for (size_t j = 0; j < ww.size(); ++j) {
      d_involutionSupport[z].set(ww[j]);
    }
  }

  for (; z < size(); ++z) {
    size_t s = firstStrictDescent(z);
    DescentStatus::Value v = descentValue(s,z);
    if (v == DescentStatus::ComplexDescent) { // do cross action
      size_t sz = cross(s,z);
      d_involutionSupport[z] = d_involutionSupport[sz];
      d_involutionSupport[z].set(s);
      d_involutionSupport[z].set(W.twistGenerator(s));
    } else { // do inverse cayley transform
      size_t sz = inverseCayley(s,z).first;
      d_involutionSupport[z] = d_involutionSupport[sz];
      d_involutionSupport[z].set(s);
    }
  }

  return;
}

void Helper::fillLengths()

/*!
  \brief Fills in the length table.

  This is directly deduced from the x-part.
*/

{
  for (BlockElt z = 0; z < d_size; ++z)
    d_length[z] = d_kgb.length(d_x[z]);

  return;
}

void Helper::makeWeylCorrelation()

/*!
  \brief Fills d_toDualWeyl.

  Explanation: this is a fairly annoying twist. Because of the normalizations
  that we do for Weyl groups, that cannot distinguish between B2 and C2, F4
  and f4, G2 and g2 (self-dual Dynkin diagrams where the selfduality is not the
  identity), it cannot be guaranteed that in all cases the normal forms we use
  for the Weyl group of the group, and that of the dual group, are compatible.
  The compatibility can only be guaranteed for the _outer_ representations of
  the Weyl group elements. So, d_toDualWeyl is the transition map:

          d_toDualWeyl[s] = d_dualin[d_out[s]]
*/

{
  using namespace weyl;

  const weyl::WeylGroup& W = d_kgb.weylGroup();
  const weyl::WeylGroup& dW = d_dualkgb.weylGroup();

  for (size_t s = 0; s < W.rank(); ++s) {
    WeylElt w;
    W.prod(w,s);
    WeylWord ww;
    dW.out(ww,w);
    d_toDualWeyl[s] = ww[0];
  }

  return;
}

void Helper::orbitPairs()

/*!
  \brief Puts in d_orbitPairs a list of all pairs of orbit parameters that
  belong to the restricted product.

  It also fills in the firstx array. For any given x in d_kgb, this gives
  the offset in d_orbitPairs of the first pair containing x.
*/

{
  using namespace kgb;
  using namespace weyl;

  KGBElt x = 0;
  BlockElt xOffset = 0;

  while (x < d_kgb.size()) {
#ifdef VERBOSE
    std::cerr << x << "\r";
#endif
    const TwistedInvolution& w = d_kgb.involution(x);
    TwistedInvolution dw = dualInvolution(w);
    KGBEltPair yRange = d_dualkgb.tauPacket(dw);
    for (; (x < d_kgb.size()) && (d_kgb.involution(x) == w); ++x) {
      for (size_t y = yRange.first; y < yRange.second; ++y) {
	d_x.push_back(x);
	d_y.push_back(y);
      }
      d_firstx[x] = xOffset;
      xOffset += yRange.second - yRange.first;
    }
  }

  assert(xOffset == d_size); // check that number of pairs is as expected

#ifdef VERBOSE
  std::cerr << std::endl;
#endif

  d_firstx[x] = xOffset; // make d_firstx[d_kgb.size()]==d_size
  // put sentinels at end of d_x and d_y
  d_x.push_back(UndefKGB);
  d_y.push_back(UndefKGB);

}

}
}


/*****************************************************************************

        Chapter III -- Functions local to blocks.cpp

  ... explain here when it is stable ...

******************************************************************************/

// namespace {

  namespace blocks {
  namespace helper {

void insertAscents(std::set<BlockElt>& hs, const set::SetEltList& hr, size_t s,
		   const Block& block)

/*!
  \brief Inserts into hs the ascents from hr through s.

  Explanation: technical function for the hasse construction, that makes the
  part of the coatom list for a given element arising from a given descent.
*/

{
  using namespace descents;

  for (size_t j = 0; j < hr.size(); ++j) {
    BlockElt z = hr[j];
    if (block.isStrictAscent(s,z)) {
      switch (block.descentValue(s,z)) {
      case DescentStatus::ComplexAscent:
	hs.insert(block.cross(s,z));
	break;
      case DescentStatus::ImaginaryTypeI:
	hs.insert(block.cayley(s,z).first);
	break;
      case DescentStatus::ImaginaryTypeII:
	hs.insert(block.cayley(s,z).first);
	hs.insert(block.cayley(s,z).second);
	break;
      default: // this cannot happen!
	break;
      }
    }
  }

 return;
}

void makeHasse(std::vector<set::SetEltList>& hd, const Block& block)

/*!
  \brief Puts in hd the hasse diagram data for the Bruhat
  ordering on the block.

  Explanation: we use the algorithm from Vogan's 1982 Park City notes.
*/

{
  using namespace descents;
  using namespace poset;

  hd.resize(block.size());

  for (BlockElt z = 0; z < block.size(); ++z) {

    std::set<BlockElt> h_z;

    for (size_t s = 0; s < block.rank(); ++s)
      if (block.isStrictDescent(s,z)) {

	DescentStatus::Value v = block.descentValue(s,z);
	BlockEltPair sz;

	switch (v) {
	case DescentStatus::ComplexDescent:
	  sz.first = block.cross(s,z);
	  sz.second = UndefBlock;
	  h_z.insert(sz.first);
	  break;
	case DescentStatus::RealTypeI: // cayley(s,z) is two-valued
	  sz = block.inverseCayley(s,z);
	  h_z.insert(sz.first);
	  h_z.insert(sz.second);
	  break;
	case DescentStatus::RealTypeII: // cayley(s,z) is single-valued
	  sz.first = block.inverseCayley(s,z).first;
	  sz.second = UndefBlock;
	  h_z.insert(sz.first);
	  break;
	default: // this cannot happen!
	  break;
	}

	insertAscents(h_z,hd[sz.first],s,block);

	if (sz.second != UndefBlock) {
	  insertAscents(h_z,hd[sz.second],s,block);
	}
      }

    std::copy(h_z.begin(),h_z.end(),std::back_inserter(hd[z]));

  }

  return;
}

}

}

}
