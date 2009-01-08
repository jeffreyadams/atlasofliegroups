/*
  This is blocks.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Modified by Marc van Leeuwen, 2007
  Part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/
/*!
\file
\brief Implementation of the class Block.

  This module aims to construct the fundamental combinatorial structure
  underlying K-L computations: a block of representations, together with
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

#include "blocks.h"

#ifdef VERBOSE
#include <iostream>
#endif

#include <cassert>
#include <memory>
#include <set>

#include "basic_io.h"
#include "bruhat.h"
#include "complexredgp.h"
#include "descents.h"
#include "kgb.h"
#include "realredgp.h"
#include "tags.h"
#include "weyl.h"
#include "hashtable.h"

/*
  Our task is fairly simple: given the one sided parameter sets for the real
  form and for the dual real form, as provided by the kgb module, which sets
  are fibred over the sets of twisted involutions for the Weyl group and dual
  Weyl group respectively, we must form the fibred product over corresponding
  pairs of twisted involutions, and equip the resulting structure with
  relations inherited from the kgb structures.

  One point to be resolved here is the relation between a twisted involution
  and the corresponding dual twisted involution; it is implemented in the
  function |dualInvolution| below. On the level of involution matrices acting
  on the character and cocharacter lattices, the relation is minus transpose;
  this has to be translated into a relation between Weyl group elements, which
  is slightly complicated by the fact that while the Weyl group and the dual
  Weyl group are the same abstract group, differing only (if at all) by the
  involution of the generating set called |twist|, the groups may differ in
  their internal representation in their numbering of the simple generators.

  Another important point here is matching the local structure of the two KGB
  sets to define the local relations in the block structure. The combinations
  possible of local structures in the KGB sets, as given by |gradings::Status|
  values, are limited by the relationship between involution and dual
  involution: a generator that is complex for one is so as well for the other,
  while the remaining generators are imaginary for one are real for the other.

  It turns out that possibilites are even more restricted then implied by the
  above relation, and every block element $z=(x,y)$ occurs, for each generator
  |s|, in one of the following five configurations. At the left we depict the
  coordinate |x| in the kgb set, to the right of it the coordinte |y| in the
  dual kgb set; at the former case going down increases the length, in the
  latter it increases the length. At the right we describe the local block
  structure; the terminology used there refers mostly to what happend in kgb,
  and in particular going down is considered to be an ascent in the block.

  If |s| is complex for the involution and its dual, one has cross actions

             ( x      ,     y   )    ComplexAscent       z
	       |            |                            |
	       |            |                            |
	     ( x'     ,     y'  )    ComplexDescent      z'

  If |s| is imaginary noncompact for the involution, it will be real, and in
  the image of the Cayley transform, for the dual involution. Either the
  Cayley image of the |x| coordinate is shared with that of another KGB
  element and the |y| coordinate has a single-valued inverse Cayley transform
  (type I situation), or the Cayley image of the |x| coordinate is unshared,
  and the |y| coordinate has a double-valued inverse Cayley transform.

      ( x      x'    , s^y   )    Imaginary Type I (twice)    z     z'
	 \    /         |                                      \   /
	  \  /          |                                       \ /
      (  s^x=s^x'    ,  y    )        Real Type I               s^z



      (    x    ,   s^y=s^y' )    Imaginary Type II              z
           |          / \                                       / \
	   |         /   \                                     /   \
      (   s^x   ,   y     y' )     Real Type II (twice)     s^z_1  s^z_2

  If |s| is imaginary compact for the involution, it will be real and not in
  the image of the Cayley transfor for the dual involution. No Cayley transorm
  will be defined for the |x| coordinate, and no inverse Cayley transform for
  the |y| coordinate, and both are fixed by the cross action; the situation is
  called ImaginaryCompact.

      (    x    ,     y     )     ImaginaryNoncompact            z

  Finally that situation with x and y interchanged is called RealNonparity

      (    x    ,     y     )         RealNonparity              z

  Although the last two cases have no cross action links for |s| to other
  block elements, nor any Cayley or inverse Cayley links for |s|, we consider
  (more in particular |descents::DescentStatus::isDescent| considers) |s| to
  be in the descent set in the ImaginaryNoncompact case, and not in the
  descent set for the RealNonparity case (note that this is opposite to the
  status of imaginary and real generators in the other (parity) cases). These
  cases do not count as strict descent/ascent however, as is indicated in the
  predicates |isStrrictDescent| and |isStrictAscent| below.
*/

namespace atlas {

namespace blocks {
namespace {


using namespace blocks;

weyl::WeylInterface
correlation(const weyl::WeylGroup& W,const weyl::WeylGroup& dW);

descents::DescentStatus descents(kgb::KGBElt x,
				 kgb::KGBElt y,
				 const kgb::KGB& kgb,
				 const kgb::KGB& dual_kgb);

void insertAscents(std::set<BlockElt>&, const set::SetEltList&, size_t,
		   const Block&);
void makeHasse(std::vector<set::SetEltList>&, const Block&);


} // namespace

} // namespace blocks

/*****************************************************************************

        Chapter I -- The Block class

******************************************************************************/

namespace blocks {

/*!
  \brief Constructor for the block class.

  Constructs a block from the datum of a real form rf for G and a real
  form df for G^vee (_not_ strong real forms: up to isomorphism, the
  result depends only on the underlying real forms!).
*/
Block::Block(complexredgp::ComplexReductiveGroup& G,
	     realform::RealForm rf, realform::RealForm df,
	     bool select_Cartans)
  : d_realForm(rf)
  , d_dualForm(df)
  , d_rank(G.semisimpleRank())
  , d_weylGroup(G.twistedWeylGroup())
  , d_xrange(0), d_yrange(0) // set by |generate|
  , d_x(), d_y(), d_first_z_of_x() // filled by |generate|
  , d_cross(d_rank), d_cayley(d_rank) // each entry filled by |generate|
  , d_descent(), d_length(), d_Cartan(), d_involution() // filled by |generate|
  , d_involutionSupport() // filled below
  , d_state()
  , d_bruhat(NULL)
{
  realredgp::RealReductiveGroup G_R(G,rf);
  complexredgp::ComplexReductiveGroup dG(G,tags::DualTag()); // the dual group
  realredgp::RealReductiveGroup dG_R(dG,df);

#ifdef VERBOSE
  std::cerr << "entering block construction... " << std::flush;
#endif

  generate(G_R,dG_R,select_Cartans); // does most of the construction work

  // all that remains is computing the supports in $S$ of twisted involutions
  d_involutionSupport.reserve(size()); // its eventual size
  for (BlockElt z=0; z<size() and d_length[z]==d_length[0]; ++z) // minl length
  {
    if (z==0 or d_involution[z]!=d_involution[z-1])
    { // compute involution support directly from definition
      bitset::RankFlags support;
      weyl::WeylWord ww=weylGroup().word(d_involution[z]);
      for (size_t j=0; j<ww.size(); ++j)
	support.set(ww[j]);
      d_involutionSupport.push_back(support);
    }
    else // unchanged involution
      d_involutionSupport.push_back(d_involutionSupport.back()); // duplicate
  }

  // complete by propagating involution supports
  for (BlockElt z=d_involutionSupport.size(); z<size(); ++z)
  {
    size_t s = firstStrictDescent(z);
    assert (s<d_rank); // must find one, as we are no longer at minimal length
    descents::DescentStatus::Value v = descentValue(s,z);
    if (v == descents::DescentStatus::ComplexDescent) // cross link
    { // use value from shorter cross neighbour, setting |s| and |twist(s)|
      d_involutionSupport[z] = d_involutionSupport[cross(s,z)];
      d_involutionSupport[z].set(s);
      d_involutionSupport[z].set(twistedWeylGroup().twisted(s));
    }
    else // Real Type I or II
    { // use (some) inverse Cayley transform and set |s|
      d_involutionSupport[z] = d_involutionSupport[inverseCayley(s,z).first];
      d_involutionSupport[z].set(s);
    }
  } // |for(z)|


#ifdef VERBOSE
  std::cerr << "done" << std::endl;
#endif
}


Block::~Block()

{
  delete d_bruhat;
}

/******** copy, assignment and swap ******************************************/

/******** accessors **********************************************************/

/*!\brief Look up element by |x|, |y| coordinates

  Precondition: |x| and |y| should be compatible: such a block element exists

  This uses the |d_first_z_of_x| table to locate the range where the |x|
  coordinates are correct; then comparing the given |y| value with the first
  one present for |x| (there must be at least one) we can predict the value
  directly, since for each fixed |x| value the values of |y| are consecutive.
*/
BlockElt Block::element(kgb::KGBElt x,kgb::KGBElt y) const
{
  BlockElt first=d_first_z_of_x[x];
  BlockElt z = first +(y-d_y[first]);
  assert(z<size() and d_x[z]==x and d_y[z]==y); // element should be found
  return z;
}

/*!
  \brief Tells if s is a strict ascent generator for z.

  Explanation: this means that d_descent[z][s] is one of ComplexAscent,
  ImaginaryTypeI or ImaginaryTypeII.
*/
bool Block::isStrictAscent(size_t s, BlockElt z) const
{
  using namespace descents;

  DescentStatus::Value v = descentValue(s,z);
  return not DescentStatus::isDescent(v) and v!=DescentStatus::RealNonparity;
}

/*!
  \brief Tells if s is a strict descent generator for z.

  Explanation: this means that d_descent[z][s] is one of ComplexDescent,
  RealTypeI or RealTypeII.
*/
bool Block::isStrictDescent(size_t s, BlockElt z) const
{
  using descents::DescentStatus;

  DescentStatus::Value v = descentValue(s,z);
  return DescentStatus::isDescent(v) and v!=DescentStatus::ImaginaryCompact;
}

/*!
  \brief Returns the first descent for z (the number of a simple root) that is
not imaginary compact, or rank() if there is no such descent.
*/
size_t Block::firstStrictDescent(BlockElt z) const
{
  for (size_t s = 0; s < rank(); ++s)
    if (isStrictDescent(s,z))
      return s;

  return rank(); // signal nothing was found
}

/*!
  \brief Returns the first descent for z (the number of a simple root) that is
either complex or real type I; if there is no such descent returns |rank()|
*/
size_t Block::firstStrictGoodDescent(BlockElt z) const
{
  for (size_t s = 0; s < rank(); ++s)
    if (isStrictDescent(s,z) and
	descentValue(s,z)!=descents::DescentStatus::RealTypeII)
      return s;

  return rank(); // signal nothing was found
}

/*!
  \brief the functor \f$T_{\alpha,\beta}\f$

  Precondition: alpha and beta are adjacent roots, of which alpha is a (weak)
  descent for y, while beta is not a descent for y.

  In fact if this is not satisfied, we return a pair of UndefBlock elements
*/

BlockEltPair Block::link(size_t alpha,size_t beta,BlockElt y) const
{
  const descents::DescentStatus& desc=descent(y);

  std::vector<BlockElt> result(2,UndefBlock); // overwritten using iterator
  std::vector<BlockElt>::iterator it=result.begin();

  BlockEltPair p=inverseCayley(alpha,y);
  switch (desc[alpha])
  {
  case descents::DescentStatus::ComplexDescent:
    {
      BlockElt y1=cross(alpha,y);
      if (isWeakDescent(beta,y1))
	*it++=y1;
      break;
    }
  case descents::DescentStatus::RealTypeI:
    if (isWeakDescent(beta,p.second))
      *it++=p.second;
    // FALL THROUGH
  case descents::DescentStatus::RealTypeII:
    if (isWeakDescent(beta,p.first))
      *it++=p.first;
    break;
  default: {}
  } // switch(desc[alpha])

  p=cayley(beta,y);
  switch (desc[beta])
  {
  case descents::DescentStatus::ComplexAscent:
    {
      BlockElt y1=cross(beta,y);
      if (not isWeakDescent(alpha,y1))
	*it++=y1;
      break;
    }
  case descents::DescentStatus::ImaginaryTypeII:
    if (not isWeakDescent(alpha,p.second))
      *it++=p.second;
    // FALL THROUGH
  case descents::DescentStatus::ImaginaryTypeI:
    if (not isWeakDescent(alpha,p.first))
      *it++=p.first;
    break;
  default: {}
  } // switch(desc[beta])

  assert(&*it<=&result[2]);

  return std::make_pair(result[0],result[1]);
}


/******** manipulators *******************************************************/

// to be called by the constructor
void Block::generate(realredgp::RealReductiveGroup& G,
		     realredgp::RealReductiveGroup& dG,
		     bool select_Cartans)
{
  kgb::KGB kgb(G,select_Cartans ? common_Cartans(G,dG) : bitmap::BitMap(0));
  kgb::KGB dual_kgb
    (dG,select_Cartans ? common_Cartans(dG,G) : bitmap::BitMap(0));

#ifdef VERBOSE
  std::cerr << "K\\G/B and dual generated... " << std::flush;
#endif

  const weyl::WeylGroup& dual_Weyl_group(dual_kgb.weylGroup());
  weyl::WeylInterface to_dual_Weyl=correlation(d_weylGroup.weylGroup(),
					       dual_Weyl_group);

  // set |d_xrange| and |d_yrange|
  d_xrange = kgb.size();
  d_yrange = dual_kgb.size();

  complexredgp::ComplexReductiveGroup& G_C = G.complexGroup();
  size_t size=G_C.block_size(d_realForm,d_dualForm); // not stored in |Block| !

  d_x.reserve(size);
  d_y.reserve(size);
  d_first_z_of_x.reserve(d_xrange+1);

  d_involution.reserve(size);
  d_length.reserve(size);
  d_Cartan.reserve(size);
  d_descent.reserve(size);

  { // generate the six tables we just dimensioned
    kgb::KGBElt x = 0;
    BlockElt base_z = 0; // block element |z| where generation for |x| starts

    while (x<kgb.size()) // loop over twisted involutions $\tau$
    {
      const weyl::TwistedInvolution& w = kgb.involution(x);

      kgb::KGBEltPair xRange = kgb.tauPacket(w);
      kgb::KGBEltPair yRange =
	dual_kgb.tauPacket(dualInvolution(w,to_dual_Weyl));

      for (assert(x==xRange.first); x < xRange.second; ++x)
      {
	d_first_z_of_x.push_back(base_z); // set even for |x| without any |y|
	for (size_t y = yRange.first; y < yRange.second; ++y)
	{
	  d_x.push_back(x);
	  d_y.push_back(y);
	  d_involution.push_back(w);
	  d_length.push_back(kgb.length(x));
	  d_Cartan.push_back(kgb.Cartan_class(x));
	  d_descent.push_back(descents(x,y,kgb,dual_kgb));
	} // |for(y)|
	base_z += yRange.second - yRange.first; // skip to next |x|
      } // |for(x)|
    } // "for(tau)"

    assert(base_z == size); // check that number of pairs is as expected
    d_first_z_of_x.push_back(base_z);// make |d_first_z_of_x[d_xrange]==d_size|

    assert(d_first_z_of_x.size()==d_xrange+1);
    assert(d_x.size()==size); // similarly for |d_y|, ..., |d_descent|
  } // end of generation of |(x,y)| tables

  // Now |element| can be safely called; install cross and Cayley tables

  for (size_t s = 0; s < d_rank; ++s)
  { // the generation below is completely independent for each |s|
    d_cross[s].resize(size);
    d_cayley[s].resize(size,std::make_pair(UndefBlock,UndefBlock));
    for (BlockElt z = 0; z<size; ++z)
    {
     d_cross[s][z] = element(kgb.cross(s,x(z)),dual_kgb.cross(s,y(z)));
     switch (d_descent[z][s])
     {
     default: break; // leave |d_cayley[s][z]|
     case descents::DescentStatus::ImaginaryTypeII:
       {
	 BlockElt z1=element(kgb.cayley(s,x(z)),
			     dual_kgb.inverseCayley(s,y(z)).second);
	 d_cayley[s][z].second = z1; // double-valued direct Cayley
	 d_cayley[s][z1].first = z; // single-valued inverse Cayley
       }
       // FALL THROUGH
     case descents::DescentStatus::ImaginaryTypeI:
       {
	 BlockElt z0=element(kgb.cayley(s,x(z)),
			     dual_kgb.inverseCayley(s,y(z)).first);
	 d_cayley[s][z].first = z0; // |.second| remains |UndefBlock| in TypeI
	 BlockEltPair& ic = d_cayley[s][z0]; // location for inverse Cayley
	 if (ic.first==UndefBlock) // it may be single or double valued
	   ic.first = z;
	 else
	   ic.second = z;
       }
     } // switch
    } // |for (z)|
  } // |for(s)|


} // generate

/*!
  \brief Constructs the BruhatOrder.

  NOTE: may throw a MemoryOverflow error. Commit-or-rollback is guaranteed.
*/
void Block::fillBruhat()
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


/******** private accessor *************************************************/


//!\brief Returns the twisted involution dual to tw.

/*
  Explanation: we have $\tau = tw.\delta$, with $tw$ in the Weyl group $W$,
  and $\delta$ the fundamental involution of the character lattice. We seek
  the twisted involution $v$ in the dual Weyl group $W^\vee$ such that
  $-\tau^t = v.\delta^\vee$; here $\delta^\vee$ is the fundamental involution
  for $W^\vee$, which acts on the cocharacter lattice and can be expressed as
  $-w_0^t\delta^t$ where $w_0$ is the longest element of $W$. So we have:

        $(-\delta^t). tw^t = v.w_0^t.(-\delta^t)$

  which leads to $v = (w_0.\delta.tw.\delta)^t$. Conjugation by $\delta$ in
  the extended Weyl group defines the automorphism $\theta$ of $W$ coming from
  the Dynkin diagram involution that defines the inner class, so we van write
  $v=(w_0.\theta(tw))^t$. This equation of endomorphisms of the cocharacter
  lattice must be translated into one on the level of abstract Weyl groups.

  Passing from a Weyl group element $w$ acting on the character lattice to the
  induced operation $w^t$ on the cocharacter lattice is an anti-isomorphism
  $D:W\to W^\vee$ that maps each generator of $W$ to the corresponding
  generator of $W^\vee$, so we want to determine $v=D(w_0.\theta(tw))$.
  Computing some $D(w)$ amounts to interpreting $w^{-1}$ in $W^\vee$. This
  reinterpretation is the identity on the level of Weyl words (the same
  numbering of generators is used in the external representation), but since
  we are using the internal representation |WeylElt| here, and the respective
  Weyl groups may use different internal numberings of generators, a call to
  |WeylGroup::translation| is necessary to do the conversion; the |to_daul|
  map should define an automorphism of the Weyl group so that one has

       |dual_W.word(W.translation(weyl::WeylElt(ww,W),to_dual)) == ww|

  for each Weyl word |ww|. Note that this difference of internal numbering
  could have been avoided if the |WeylGroup| constructor had used the Coxeter
  matrix to determine the internal numbering, since this is guaranteed to be
  the same for $W$ and $W^\vee$; however since the Cartan matrix and normal
  renumbering of the corresponding Dynkin diagram is used instead, $W$ and
  $W^\vee$ will differ internally in the cases $B_2$, $G_2$ and $F_4$.
*/
weyl::TwistedInvolution Block::dualInvolution
  (const weyl::TwistedInvolution& tw,weyl::WeylInterface to_dual) const
{
  const weyl::TwistedWeylGroup& tW = twistedWeylGroup();
  const weyl::WeylGroup& W = tW.weylGroup();
  return weyl::TwistedInvolution
    (W.translation(W.inverse(W.prod(W.longest(),tW.twisted(tw.w()))),
		   to_dual));
}

} // namespace blocks


/*****************************************************************************

        Chapter III -- Functions local to blocks.cpp

******************************************************************************/

namespace blocks {
  namespace {

/*!
  \brief Returns mapping of internal numberings from |W| to |dW|

  Explanation: this is a fairly annoying twist. Because of the normalizations
  that we do for Weyl groups based on the Cartan matrix rather than the
  Coxeter matrix, the dual Weyl group may uses a different internal numbering
  than the Weyl group in the (irreducible) cases $B_2$, $G_2$ ad $F_4$. This
  means we cannot reinterpret |WeylElt| values for $W$ as values for $dW$
  without a conversion. The conversion will be performed by
  |WeylGroup::translate| based on the value computed here, which is a mapping
  of external numberings, to be applied by |W| on an element just before
  reinterpreting it in the |dW|. We should have for each |s|, in terms of the
  internal translation arrays:

        |  dW.d_out[W.d_in[d_toDualWeyl[s]]] = s | or equivalently
        |  d_toDualWeyl[s] = W.d_out[dW.d_in[s]] |

  so that the reinterpretation will preserve the outer representation. We do
  not having direct access to those arrays, but we can pass into internal
  representation and back out for another group without acessing them
  explicitly. Note this could not possibly be made to work (in all cases) if
  |WeylGroup::translate| were to use internal numbering, as it originally did.
*/
weyl::WeylInterface
correlation(const weyl::WeylGroup& W,const weyl::WeylGroup& dW)
{
  size_t rank=W.rank();
  weyl::WeylInterface result;
  for (size_t s = 0; s < rank; ++s)
  {
    weyl::WeylElt w=dW.generator(s); // converts |s| to inner numbering |dW|
    weyl::WeylWord ww=W.word(w); // interpret |w| in |dW|; gives singleton
    assert(ww.size()==1);

    /* We want to map |s| to |ww[0]| so that interpreting that internally in
       |W| gives the element that in |dW| will represent |s|
    */
    result[s] = ww[0];
  }
  return result;
}

descents::DescentStatus descents(kgb::KGBElt x,
				 kgb::KGBElt y,
				 const kgb::KGB& kgb,
				 const kgb::KGB& dual_kgb)
{
  using gradings::Status;
  using descents::DescentStatus;

  DescentStatus result;
  for (size_t s = 0; s < kgb.rank(); ++s)
    if (kgb.status(s,x) == Status::Complex) // s is complex
      if (kgb.isDescent(s,x))
	result.set(s,DescentStatus::ComplexDescent);
      else
	result.set(s,DescentStatus::ComplexAscent);
    else if (kgb.status(s,x) == Status::ImaginaryNoncompact)
      if (kgb.cross(s,x)!=x) // type I
	result.set(s,DescentStatus::ImaginaryTypeI);
      else // type II
	result.set(s,DescentStatus::ImaginaryTypeII);
    else if (dual_kgb.status(s,y) == Status::ImaginaryNoncompact)
      if (dual_kgb.cross(s,y)!=y) // type II
	result.set(s,DescentStatus::RealTypeII);
      else // type I
	result.set(s,DescentStatus::RealTypeI);
    // now s is imaginary compact or real nonparity
    else if (kgb.status(s,x) == Status::Real)
      result.set(s,DescentStatus::RealNonparity);
    else
      result.set(s,DescentStatus::ImaginaryCompact);

  return result;
}

/*!
  \brief Inserts into hs the ascents from hr through s.

  Explanation: technical function for the Hasse construction, that makes the
  part of the coatom list for a given element arising from a given descent.
*/
void insertAscents(std::set<BlockElt>& hs, const set::SetEltList& hr, size_t s,
		   const Block& block)
{
  using descents::DescentStatus;

  for (size_t j = 0; j < hr.size(); ++j)
  {
    BlockElt z = hr[j];
    switch (block.descentValue(s,z))
    {
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
    default: // not a strict ascent
      break;
    }
  }
}


/*!
  \brief Puts into |Hasse| the hasse diagram data for the Bruhat
  ordering on |block|.

  Explanation: we used the algorithm from Vogan's 1982 Park City notes...
  which contains a bad definition. Now modified to work like kgb makeHasse:
  seek an ascent s that is complex or type I real. If it exists, use it as in
  kgb. If it doesn't then we're essentially at a split principal series. The
  immediate predecessors of z are just the inverse Cayley transforms.
*/
void makeHasse(std::vector<set::SetEltList>& Hasse, const Block& block)
{
  using descents::DescentStatus;

  Hasse.resize(block.size());

  for (BlockElt z = 0; z < block.size(); ++z)
  {
    std::set<BlockElt> h_z;

    size_t s=block.firstStrictGoodDescent(z);
    if (s<block.rank())
      switch (block.descentValue(s,z))
      {
      default: assert(false); break;
      case DescentStatus::ComplexDescent:
	{
	  BlockElt sz = block.cross(s,z);
	  h_z.insert(sz);
	  insertAscents(h_z,Hasse[sz],s,block);
	}
	break;
      case DescentStatus::RealTypeI: // inverseCayley(s,z) is two-valued
	{
	  BlockEltPair sz = block.inverseCayley(s,z);
	  h_z.insert(sz.first);
	  h_z.insert(sz.second);
	  insertAscents(h_z,Hasse[sz.first],s,block);
	}
      }
    else // now just gather all RealTypeII descents of |z|
      for (size_t s = 0; s < block.rank(); ++s)
	if (block.descentValue(s,z)==DescentStatus::RealTypeII)
	  h_z.insert(block.inverseCayley(s,z).first);

    std::copy(h_z.begin(),h_z.end(),std::back_inserter(Hasse[z])); // set->list
  } // for |z|
}

} // namespace

/*****************************************************************************

      Chapter IV -- Functions exported from blocks.cpp (declared in blocks.h)

******************************************************************************/

std::vector<BlockElt> dual_map(const Block& b, const Block& dual_b)
{

  assert(b.size()==dual_b.size());

  std::vector<BlockElt> result(b.size());
  for (BlockElt i=0; i<b.size(); ++i)
    result[i]=dual_b.element(b.y(i),b.x(i));

  return result;
}

bitmap::BitMap common_Cartans(realredgp::RealReductiveGroup& GR,
			      realredgp::RealReductiveGroup& dGR)
  { return GR.Cartan_set()
      & GR.complexGroup().dual_Cartan_set(dGR.realForm());
  }

} // namespace blocks

} // namespace atlas
