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
#include <vector>
#include <deque>
#include <set> // for |insertAscents|

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
  dual kgb set; in the former case going down increases the length, in the
  latter it decreases the length. At the right we describe the local block
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
      (  s^x=s^x'    ,  y    )         Real Type I              s^z



      (    x    ,   s^y=s^y' )       Imaginary Type II           z
           |          / \                                       / \
	   |         /   \                                     /   \
      (   s^x   ,   y     y' )     Real Type II (twice)     s^z_1  s^z_2

  If |s| is imaginary compact for the involution, it will be real and not in
  the image of the Cayley transform for the dual involution. No Cayley
  transorm will be defined for the |x| coordinate, and no inverse Cayley
  transform for the |y| coordinate, and both are fixed by the cross action;
  the situation is called ImaginaryCompact.

      (    x    ,     y     )       ImaginaryCompact             z

  Finally that situation with x and y interchanged is called RealNonparity

      (    x    ,     y     )         RealNonparity              z

  Although the last two cases have no cross action links for |s| to other
  block elements, nor any Cayley or inverse Cayley links for |s|, we consider
  (more in particular |descents::DescentStatus::isDescent| considers) |s| to
  be in the descent set in the ImaginaryCompact case, and not in the descent
  set for the RealNonparity case (note that this is opposite to the status of
  imaginary and real generators in the other (parity) cases). These cases do
  not count as strict descent/ascent however, as is indicated in the
  predicates |isStrictDescent| and |isStrictAscent| below.
*/

namespace atlas {

namespace blocks {
namespace {

weyl::WeylInterface
correlation(const weyl::WeylGroup& W,const weyl::WeylGroup& dW);

descents::DescentStatus descents(kgb::KGBElt x,
				 kgb::KGBElt y,
				 const kgb::KGB_base& kgb,
				 const kgb::KGB_base& dual_kgb);

void insertAscents(std::set<BlockElt>&, const set::SetEltList&, size_t,
		   const Block&);
void makeHasse(std::vector<set::SetEltList>&, const Block&);


} // namespace

} // namespace blocks

/*****************************************************************************

        Chapter I -- The Block class

******************************************************************************/

namespace blocks {

Block_base::Block_base(const kgb::KGB& kgb,const kgb::KGB& dual_kgb)
  : W(kgb.twistedWeylGroup())
  , xrange(kgb.size()), yrange(dual_kgb.size())
  , d_x(), d_y(), d_first_z_of_x() // filled below
  , d_cross(rank()), d_cayley(rank()) // each entry filled below
  , d_descent(), d_length() // filled below
{
  const weyl::TwistedWeylGroup& dual_W =dual_kgb.twistedWeylGroup();

  std::vector<weyl::TwistedInvolution> dual_w;
  dual_w.reserve(kgb.nr_involutions());
  size_t size=0;
  for (unsigned int i=0; i<kgb.nr_involutions(); ++i)
  {
    const weyl::TwistedInvolution w = kgb.nth_involution(i);
    dual_w.push_back(dual_involution(w,W,dual_W));
    size += kgb.packet_size(w)*dual_kgb.packet_size(dual_w.back());
  }

  d_first_z_of_x.reserve(xrange+1);

  d_x.reserve(size);
  d_y.reserve(size);
  d_length.reserve(size);
  d_descent.reserve(size);

  BlockElt base_z = 0; // block element |z| where generation for |x| starts

  for (unsigned int i=0; i<kgb.nr_involutions(); ++i)
  {
    const weyl::TwistedInvolution w = kgb.nth_involution(i);

    kgb::KGBEltPair x_step = kgb.tauPacket(w);
    kgb::KGBEltPair y_step = dual_kgb.tauPacket(dual_w[i]);

    for (kgb::KGBElt x=x_step.first; x<x_step.second; ++x)
    {
      d_first_z_of_x.push_back(base_z); // set even for |x| without any |y|
      for (size_t y=y_step.first; y<y_step.second; ++y)
      {
	d_x.push_back(x);
	d_y.push_back(y);
	d_descent.push_back(descents(x,y,kgb,dual_kgb));
	d_length.push_back(kgb.length(x));
      } // |for(y)|
      base_z += y_step.second - y_step.first; // skip to next |x|
    } // |for(x)|
  } // |for (i)|
  d_first_z_of_x.push_back(base_z);// make |d_first_z_of_x[d_xrange]==d_size|

  assert(d_first_z_of_x.size()==xrange+1); // all |x|'s have been seen

  // Now |element| can be safely called; install cross and Cayley tables

  for (weyl::Generator s = 0; s<rank(); ++s)
  { // the generation below is completely independent for each |s|
    d_cross[s].resize(size);
    d_cayley[s].resize(size,std::make_pair(UndefBlock,UndefBlock));
    for (BlockElt z=0; z<size; ++z)
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
} // |Block_base::Block_base|


Block_base::Block_base
  (realredgp::RealReductiveGroup& GR,
   const subdatum::SubSystem& sub, // at the dual side
   kgb::KGBElt x,
   const latticetypes::Weight& lambda, // discrete part of parameter
   const latticetypes::RatWeight gamma) // infinitesimal character
: W(GR.twistedWeylGroup())
, xrange(0), yrange(0),d_x(),d_y()
,d_first_z_of_x(),d_cross(sub.rank()),d_cayley(sub.rank())
,d_descent(),d_length()
{
  const kgb::KGB& kgb = GR.kgb();
  const complexredgp::ComplexReductiveGroup& G = GR.complexGroup();
  size_t our_rank = sub.rank(); // this is independent of ranks in |GR|
  weyl::WeylWord dual_involution;

  const tits::GlobalTitsGroup Tg
    (sub,G.involutionMatrix(kgb.involution(x)),dual_involution);

  weyl::TwistedInvolution tw = Tg.weylGroup().element(dual_involution);

  // step 1: get the correct value |y|
  tits::TorusElement t
    (gamma-latticetypes::RatWeight(lambda,1)/=2); // for original group $G$
  tits::GlobalTitsElement y=tits::GlobalTitsElement(t,tw);
  Tg.add(y, kgb.half_rho() - Tg.torus_part_offset () ); // for $G(\gamma)$

  // step 2: move to the minimal fiber
  { // modify |x| and |y|, descending to minimal element for |subsys|
    weyl::Generator s;
    do
    {
      for(s=0; s<our_rank; ++s)
      {
	kgb::KGBElt xx=kgb.cross_act(sub.to_simple(s),x);
	if (kgb.isDescent(sub.simple(s),xx))
	{
	  if (kgb.status(sub.simple(s),xx)==gradings::Status::Complex)
	  {
	    x = kgb.cross_act(kgb.cross(sub.simple(s),xx),sub.to_simple(s));
	    Tg.cross(s,y);
	    break;
	  }
	  else // imaginary
	    if (not Tg.compact(s,y))
	    {
	      xx=kgb.inverseCayley(sub.simple(s),xx).first; // choose one
	      x = kgb.cross_act(xx,sub.to_simple(s));
	      y = Tg.Cayley(s,y);
	      break;
	    }
	} // |if(isDescent)|
      } // |for(s)|
    } while(s<our_rank); // loop until no descents found in |subsys|

    // one might reduce the torus part of |y| here
  } // end of step 2

  weyl::TI_Entry::Pooltype pool;
  hashtable::HashTable<weyl::TI_Entry,unsigned int> hash(pool);
  kgb::GlobalFiberData gfd(sub,hash);
  gfd.add_class(sub,Tg,y.tw());

  std::vector<tits::GlobalTitsElement> y_rep(1,y);
  kgb::KGB_elt_entry::Pooltype y_pool // parallel array with fingerprints
    (1,kgb::KGB_elt_entry(gfd.fingerprint(y),y.tw()));
  hashtable::HashTable<kgb::KGB_elt_entry,kgb::KGBElt> y_hash(y_pool);

  // step 3: generate imaginary fiber-orbit of |x|'s (|y|'s are unaffected)
  std::vector<unsigned int> x_of(kgb.size(),~0);
  std::vector<kgb::KGBElt> kgb_nr_of; // indexed by newly assigned |x| numbers
  kgb_nr_of.reserve(kgb.size()); // maybe wasteful, this avoids reallocation
  {
    cartanclass::InvolutionData id =
      cartanclass::InvolutionData::build(sub,Tg,y.tw());
    // generating reflections are for \emph{real} roots for involution |y.tw()|
    const rootdata::RootList gen_root = id.real_basis();
    const subdatum::SubSystem imaginary(kgb.rootDatum(),gen_root);

    kgb::KGBEltPair p = kgb.packet(x);
    bitmap::BitMap seen(p.second-p.first);
    bitmap::BitMap news = seen;
    news.insert(x-p.first);
    while (news.andnot(seen)) // condition modifies |news| as side effect
    {
      unsigned i=news.front();
      seen.insert(i); // so |i| will be removed from |news| at end of loop
      kgb::KGBElt xx=i+p.first;
      for (weyl::Generator s=0; s<imaginary.rank(); ++s)
	news.insert(kgb.cross_act(imaginary.reflection(s),xx)-p.first);
    }

    // now insert elements from |seen| as first $x$-fiber of block
    size_t fs = seen.size(); // this is lower bound for final size, so reserve
    d_x.reserve(fs); d_y.reserve(fs); d_first_z_of_x.reserve(fs+1);
    d_length.reserve(fs); d_descent.reserve(fs);
    for (weyl::Generator s=0; s<our_rank; ++s)
      { d_cross[s].reserve(fs); d_cayley[s].reserve(fs); }

    for (bitmap::BitMap::iterator it=seen.begin(); it(); ++it)
    {
      kgb::KGBElt kgb_nr = *it+p.first;
      x_of[kgb_nr]=kgb_nr_of.size();
      kgb_nr_of.push_back(kgb_nr);
      d_first_z_of_x.push_back(d_x.size());
      d_x.push_back(x_of[kgb_nr]);
      assert(y_hash.find(kgb::KGB_elt_entry(gfd.fingerprint(y),y.tw()))==0);
      d_y.push_back(0);
      d_descent.push_back(descents::DescentStatus());
      d_length.push_back(0); // we don't know the length of |x| in the sub-KGB
    }

    d_first_z_of_x.push_back(d_x.size()); // ensure one more entry is defined
  } // end of step 3

  // step 4: generate packets for successive involutions

  {
    std::deque<BlockElt> queue(1,d_x.size()); // queue of packet boundaries
    BlockElt next=0;
    while (next<queue.front())
    { // process involution packet of elements from |next| to |queue.front()|
      std::vector<kgb::KGBElt> ys; ys.reserve(0x100);
      for (BlockElt z=next; z<d_x.size() and d_x[z]==d_x[next]; ++z)
	ys.push_back(d_y[z]);
      assert((queue.front()-next)%ys.size()==0);
      unsigned int nr_x= (queue.front()-next)/ys.size();

      for (weyl::Generator s=0; s<our_rank; ++s)
      {
	tits::GlobalTitsElement sample = y_rep[ys[next]];
	int d = Tg.cross(s,sample);
	int l = d_length[next]+d;
	assert(l>=0); // if fails, then starting point not minimal for length

	std::vector<BlockElt>& cross_link = d_cross[s]; // a safe reference

	bool is_new; // whether cross action gives a new twisted involution
	std::vector<kgb::KGBElt> cross_ys(ys.size());
	{
	  size_t old_size =  y_hash.size();
	  for (unsigned int i=0; i<ys.size(); ++i)
	  {
	    tits::GlobalTitsElement y = y_rep[ys[i]];
	    int dd = Tg.cross(s,y);
	    assert(dd==d); // all length changes should be equal
	    cross_ys[i] =
	      y_hash.match(kgb::KGB_elt_entry(gfd.fingerprint(y),y.tw()));
	    if (y_hash.size()>y_rep.size())
	      y_rep.push_back(y);
	    assert(y_hash.size()==y_rep.size());
	  }
	  is_new = y_hash.size()>old_size;
	  assert(not is_new or y_hash.size()==old_size+ys.size()); // all same
	}
	for (unsigned int i=0; i<nr_x; ++i)
	{
	  BlockElt base_z = next+i*ys.size(); // first element of R-packet
	  kgb::KGBElt n = kgb_nr_of[d_x[base_z]];
	  kgb::KGBElt s_x_n = kgb.cross_act(sub.reflection(s),n);
	  assert (is_new == (x_of[s_x_n]==(unsigned int)~0));
	  if (is_new) // add a new R-packet
	  {
	    x_of[s_x_n] = kgb_nr_of.size();
	    kgb_nr_of.push_back(s_x_n);
	    for (unsigned int j=0; j<ys.size(); ++j)
	    {
	      // install link the required this element
	      assert(d_x.size()==d_y.size()); // number of new block element
	      cross_link[base_z+j] = d_x.size();

	      d_x.push_back(x_of[s_x_n]);
	      d_y.push_back(cross_ys[j]);
	      d_descent.push_back(descents::DescentStatus()); // set later
	      d_length.push_back(l);

	    } // |for(j)|
	    d_first_z_of_x.push_back(d_x.size()); // mark end of R-packet
	  }
	  else // install cross links to previously existing element
	    for (unsigned int j=0; j<ys.size(); ++j)
	      cross_link[base_z+j] = element(x_of[s_x_n],cross_ys[j]);
	  // |if(is_new)|
	} // |for(i)|
      } // |for(s)|
      next=queue.front();
      queue.pop_front();
    }
  }
}


/*!\brief Look up element by |x|, |y| coordinates

  Precondition: |x| and |y| should be compatible: such a block element exists

  This uses the |d_first_z_of_x| table to locate the range where the |x|
  coordinates are correct; then comparing the given |y| value with the first
  one present for |x| (there must be at least one) we can predict the value
  directly, since for each fixed |x| value the values of |y| are consecutive.
*/
BlockElt Block_base::element(kgb::KGBElt x,kgb::KGBElt y) const
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
bool Block_base::isStrictAscent(size_t s, BlockElt z) const
{
  descents::DescentStatus::Value v = descentValue(s,z);
  return not descents::DescentStatus::isDescent(v)
    and v!=descents::DescentStatus::RealNonparity;
}

/*!
  \brief Tells if s is a strict descent generator for z.

  Explanation: this means that d_descent[z][s] is one of ComplexDescent,
  RealTypeI or RealTypeII.
*/
bool Block_base::isStrictDescent(size_t s, BlockElt z) const
{
  descents::DescentStatus::Value v = descentValue(s,z);
  return descents::DescentStatus::isDescent(v)
    and v!=descents::DescentStatus::ImaginaryCompact;
}

/*!
  \brief Returns the first descent for z (the number of a simple root) that is
not imaginary compact, or rank() if there is no such descent.
*/
size_t Block_base::firstStrictDescent(BlockElt z) const
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
size_t Block_base::firstStrictGoodDescent(BlockElt z) const
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

BlockEltPair Block_base::link(size_t alpha,size_t beta,BlockElt y) const
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



Block::Block(const kgb::KGB& kgb,const kgb::KGB& dual_kgb)
  : Block_base(kgb,dual_kgb)
  , d_Cartan(), d_involution(), d_involutionSupport() // filled below
  , d_state()
  , d_bruhat(NULL)
{
  d_Cartan.reserve(size());
  d_involution.reserve(size());
  for (BlockElt z=0; z<size(); ++z)
  {
    kgb::KGBElt xx=x(z);
    d_Cartan.push_back(kgb.Cartan_class(xx));
    d_involution.push_back(kgb.involution(xx));
  }

  compute_supports();
}

/*!
  \brief Construction function for the |Block| class.

  Constructs a block from the datum of a real form rf for G and a real
  form df for G^vee (_not_ strong real forms: up to isomorphism, the
  result depends only on the underlying real forms!).
*/
Block // pseudo constructor that ends calling main contructor
  Block::build(complexredgp::ComplexReductiveGroup& G,
	       realform::RealForm rf, realform::RealForm drf,
	       bool select_Cartans)
{
  realredgp::RealReductiveGroup G_R(G,rf);
  complexredgp::ComplexReductiveGroup dG(G,tags::DualTag()); // the dual group
  realredgp::RealReductiveGroup dG_R(dG,drf);

  kgb::KGB kgb
    ( G_R,select_Cartans ? common_Cartans(G_R,dG_R) : bitmap::BitMap(0));
  kgb::KGB dual_kgb
    (dG_R,select_Cartans ? common_Cartans(dG_R,G_R) : bitmap::BitMap(0));

#ifdef VERBOSE
  std::cerr << "K\\G/B and dual generated... " << std::flush;
#endif

  return Block(kgb,dual_kgb); // |kgb| and |dual_kgb| disappear afterwards!
}

// compute the supports in $S$ of twisted involutions
void Block::compute_supports()
{
  d_involutionSupport.reserve(size()); // its eventual size
  for (BlockElt z=0; z<size() and length(z)==length(0); ++z) // minl length
  {
    if (z==0 or involution(z)!=involution(z-1))
    { // compute involution support directly from definition
      bitset::RankFlags support;
      weyl::WeylWord ww=weylGroup().word(involution(z));
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
    assert (s<rank()); // must find one, as we are no longer at minimal length
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

Block::~Block() { delete d_bruhat; } // type of |d_bruhat| is complete here
Block::Block(const Block& b)
  : Block_base(b) // copy
  , d_Cartan(b.d_Cartan)
  , d_involution(b.d_involution)
  , d_involutionSupport(b.d_involutionSupport)
  , d_state(b.d_state)
  , d_bruhat(NULL) // possibly filled below
{
#ifdef VERBOSE // then show that we're called (does not actually happen)
  std::cerr << "copying a block" << std::endl;
#endif

  if (b.d_bruhat!=NULL) // this does not happen when called from |Block::build|
    d_bruhat = new bruhat::BruhatOrder(*b.d_bruhat);
}


/******** copy, assignment and swap ******************************************/

/******** accessors **********************************************************/


/******** manipulators *******************************************************/

/*!
  \brief Constructs the BruhatOrder.

  NOTE: may throw a MemoryOverflow error. Commit-or-rollback is guaranteed.
*/
void Block::fillBruhat()
{
  if (d_state.test(BruhatConstructed)) // work was already done
    return;

  try {
    std::vector<set::SetEltList> hd; makeHasse(hd,*this);
    bruhat::BruhatOrder* bp = new bruhat::BruhatOrder(hd); // may throw here

    // commit
    delete d_bruhat; // this is overly careful: it must be NULL
    d_bruhat = bp;
    d_state.set(BruhatConstructed);
  }
  catch (std::bad_alloc) { // transform failed allocation into MemoryOverflow
    throw error::MemoryOverflow();
  }
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
  Coxeter matrix, the dual Weyl group may use a different internal numbering
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
  not have direct access to those arrays, but we can pass into internal
  representation and back out for another group without acessing them
  explicitly. Note this could not possibly be made to work (in all cases) if
  |WeylGroup::translate| were to use internal numbering, as it originally did.

  This function is left for historical reasons, but is doubly useless, since
  we now share the (untwisted) Weyl group between an inner class and its dual,
  (see the |assert| below), and because we now avoid any reinterpretation of
  elements from the twisted Weylgroup to its dual without passing through the
  external |WeylWord| representation (see |dual_involution| below).
*/
weyl::WeylInterface
correlation(const weyl::WeylGroup& W,const weyl::WeylGroup& dW)
{
  assert(&W==&dW); // so groups are identical, making this function useless!

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
				 const kgb::KGB_base& kgb,
				 const kgb::KGB_base& dual_kgb)
{
  descents::DescentStatus result;
  for (size_t s = 0; s < kgb.rank(); ++s)
    if (kgb.status(s,x) == gradings::Status::Complex) // s is complex
      if (kgb.isDescent(s,x))
	result.set(s,descents::DescentStatus::ComplexDescent);
      else
	result.set(s,descents::DescentStatus::ComplexAscent);
    else if (kgb.status(s,x) == gradings::Status::ImaginaryNoncompact)
      if (kgb.cross(s,x)!=x) // type I
	result.set(s,descents::DescentStatus::ImaginaryTypeI);
      else // type II
	result.set(s,descents::DescentStatus::ImaginaryTypeII);
    else if (dual_kgb.status(s,y) == gradings::Status::ImaginaryNoncompact)
      if (dual_kgb.cross(s,y)!=y) // type II
	result.set(s,descents::DescentStatus::RealTypeII);
      else // type I
	result.set(s,descents::DescentStatus::RealTypeI);
    // now s is imaginary compact or real nonparity
    else if (kgb.status(s,x) == gradings::Status::Real)
      result.set(s,descents::DescentStatus::RealNonparity);
    else
      result.set(s,descents::DescentStatus::ImaginaryCompact);

  return result;
}

/*!
  \brief Inserts into hs the ascents from hr through s.

  Explanation: technical function for the Hasse construction, that makes the
  part of the coatom list for a given element arising from a given descent.
*/
void insertAscents(std::set<BlockElt>& hs,
		   const set::SetEltList& hr,
		   size_t s,
		   const Block& block)
{
  for (size_t j = 0; j < hr.size(); ++j)
  {
    BlockElt z = hr[j];
    switch (block.descentValue(s,z))
    {
    case descents::DescentStatus::ComplexAscent:
      hs.insert(block.cross(s,z));
      break;
    case descents::DescentStatus::ImaginaryTypeI:
      hs.insert(block.cayley(s,z).first);
      break;
    case descents::DescentStatus::ImaginaryTypeII:
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
  Hasse.resize(block.size());

  for (BlockElt z = 0; z < block.size(); ++z)
  {
    std::set<BlockElt> h_z;

    size_t s=block.firstStrictGoodDescent(z);
    if (s<block.rank())
      switch (block.descentValue(s,z))
      {
      default: assert(false); break;
      case descents::DescentStatus::ComplexDescent:
	{
	  BlockElt sz = block.cross(s,z);
	  h_z.insert(sz);
	  insertAscents(h_z,Hasse[sz],s,block);
	}
	break;
      case descents::DescentStatus::RealTypeI: // inverseCayley(s,z) two-valued
	{
	  BlockEltPair sz = block.inverseCayley(s,z);
	  h_z.insert(sz.first);
	  h_z.insert(sz.second);
	  insertAscents(h_z,Hasse[sz.first],s,block);
	}
      }
    else // now just gather all RealTypeII descents of |z|
      for (size_t s = 0; s < block.rank(); ++s)
	if (block.descentValue(s,z)==descents::DescentStatus::RealTypeII)
	  h_z.insert(block.inverseCayley(s,z).first);

    std::copy(h_z.begin(),h_z.end(),std::back_inserter(Hasse[z])); // set->list
  } // for |z|
}

} // namespace

/*****************************************************************************

      Chapter IV -- Functions exported from blocks.cpp (declared in blocks.h)

******************************************************************************/


//!\brief Returns the twisted involution dual to |w|.

/*
  We have $\tau = w.\delta$, with $tw$ in the Weyl group $W$, and $\delta$ the
  fundamental involution of $X^*$. We seek the twisted involution $v$ in the
  dual twisted Weyl group |dual_W| such that $-\tau^t = v.\delta^\vee$. Here
  $\delta^\vee$ is the dual fundamental involution, which acts on $X_*$ as
  minus the transpose of the longest involution $w_0\delta$, where $w_0$ is
  the longest element of $W$. This relation relation can be defined
  independently of being a twisted involution, and leads to a bijection $f$:
  $W \to W$ that is characterised by $f(e)=w_0$ and $f(s.w)=f(w)dwist(s)$ for
  any simple generator $e$ and $w\in W$, where $dwist$ is the twist of the
  dual twisted Weyl group (one also has $f(w.twist(s))=s.f(w)$ so that $f$ is
  intertwines twisted conjugation: $f(s.w.twist(s))=s.f(w).dwist(s)$).

  The implementation below is based directly on the above characterisation, by
  converting |w| into a WeylWord|, and then starting from $w_0$
  right-multiplying successively by the $dwist$ image of the letters of the
  word taken right-to-left. The only thing assumed common to |W| and |W_dual|
  is the external numbering of generators (letters in |ww|), without which the
  notion of dual twisted involution would make no sense; notably the result
  will be correct (when interpreted for |dual_W|) even if the underlying
  (untwisted) Weyl groups of |W| and |dual_W| should differ in their
  external-to-internal mappings. If one assumes that |W| and |dual_W| share
  the same underlying Weyl group (as is currently the case for an inner class
  and its dual) then one could alternatively say simply

    |return W.prod(W.weylGroup().longest(),dual_W.twisted(W.inverse(w)))|

  or

    |return W.prod(W.inverse(W.twisted(w)),W.weylGroup().longest())|
*/
weyl::TwistedInvolution
dual_involution(const weyl::TwistedInvolution& w,
		const weyl::TwistedWeylGroup& W,
		const weyl::TwistedWeylGroup& dual_W)
{
  weyl::WeylWord ww= W.word(w);
  weyl::TwistedInvolution result = dual_W.weylGroup().longest();
  for (size_t i=ww.size(); i-->0; )
    dual_W.mult(result,dual_W.twisted(ww[i]));
  return result;
}


std::vector<BlockElt> dual_map(const Block_base& b, const Block_base& dual_b)
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
