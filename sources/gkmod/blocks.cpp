/*
  This is blocks.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2007--2016 Marc van Leeuwen
  Part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "blocks.h"

#include <cassert>
#include <vector>
#include <set> // for |insertAscents|
#include <algorithm>
#include <iterator>

#include "arithmetic.h"

#include "tags.h"
#include "hashtable.h"

#include "bruhat.h"	// construction
#include "innerclass.h"
#include "realredgp.h"
#include "subsystem.h"
#include "y_values.h"
#include "kgb.h"
#include "weyl.h"
#include "ext_block.h"  // class |ext_block::ext_block| constructor
#include "kl.h"		// destruction

/*
  Our task in the traditional setup (the constructor for |Block|) is fairly
  simple: given the one sided parameter sets (|KGB|) for the real form and for
  the dual real form, which are fibred over the sets of twisted involutions
  for the Weyl group and dual Weyl group respectively, we must form the fibred
  product over corresponding pairs of twisted involutions, and equip the
  resulting structure with relations inherited from the kgb structures.

  One point to be resolved here is the relation between a twisted involution
  and the corresponding dual twisted involution; it is implemented in the
  function |dual_involution| below. On the level of involution matrices acting
  on the character and cocharacter lattices, the relation is minus transpose;
  this has to be translated into a relation between Weyl group elements.

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
  structure; the terminology used there refers mostly to what happened for |x|,
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

  When $\alpha,\alpha^\vee$ are the root and coroot for |s|, the involution of
  $X^*$ at the more compact Cartan (downstairs) is $\theta$ and at the more
  split Cartan (upstairs) $\theta'=s_\alpha*\theta$, then one has equivalences

    \alpha^\vee\notin X_*(1+\theta) \iff \alpha\in(1-\theta')X^* \iff type 1

    \alpha^\vee\in X_*(1+\theta) \iff \alpha\notin(1-\theta')X^* \iff type 2

  Moreover  $\<\alpha^\vee, (X^*)^\theta> = n\Z$  in type $n\in\{1,2\}$ and
  also  $< (X_*)^{-\theta'}, \alpha > = (2/n)\Z$  in type $n\in\{1,2\}$

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
  (more in particular |DescentStatus::isDescent| considers) |s| to
  be in the descent set in the ImaginaryCompact case, and not in the descent
  set for the RealNonparity case (note that this is opposite to the status of
  imaginary and real generators in the other (parity) cases). These cases do
  not count as strict descent/ascent however, as is indicated in the
  predicate methods |isStrictDescent| and |isStrictAscent| below.
*/

namespace atlas {

namespace blocks {




namespace {
  // some auxiliary functions used by methods, but defined near end of file

  // compute descent status of block element, based on its $(x,y)$ parts
DescentStatus descents(KGBElt x, KGBElt y,
		       const KGB_base& kgb, const KGB_base& dual_kgb);

  // compute Hasse diagram of the Bruhat order of a block
std::vector<Poset::EltList> complete_Hasse_diagram
  (const Block_base&,
   std::vector<std::unique_ptr<BlockEltList> >& partial_Hasse_diagram);


} // |namespace|

/*****************************************************************************

        Chapter I -- The |Block_base| class

******************************************************************************/

// an auxiliary function:
// we often need to fill the first empty slot of a |BlockEltPair|
BlockElt& first_free_slot(BlockEltPair& p)
{
  if (p.first==UndefBlock)
    return p.first;
  else
  {
    assert(p.second==UndefBlock); // there should be an empty slot left
    return p.second;
  }
}

Block_base::Block_base(const KGB& kgb)
  : info(), data(kgb.rank()), orbits()
  , dd(kgb.innerClass().rootDatum().cartanMatrix())
  , partial_Hasse_diagram()
  , d_bruhat(nullptr)
  , kl_tab_ptr(nullptr)
{
} // |Block_base::Block_base|

// an almost trivial constructor used for non-integral block derived types
Block_base::Block_base(unsigned int integral_rank)
  : info(), data(integral_rank), orbits()
  , dd()
  , partial_Hasse_diagram()
  , d_bruhat(nullptr)
  , kl_tab_ptr(nullptr)
{}

Block_base::Block_base(const Block_base& b) // copy constructor
  : info(b.info), data(b.data), orbits(b.orbits)
  , dd(b.dd)
  , partial_Hasse_diagram()
  , d_bruhat(nullptr) // don't care to copy; is empty in |Block::build| anyway
  , kl_tab_ptr(nullptr)  // likewise
{
#ifdef VERBOSE // then show that we're called (does not actually happen)
  std::cerr << "copying a block" << std::endl;
#endif
}

Block_base::~Block_base() = default; // but calls deleters implicitly

RankFlags Block_base::descent_generators (BlockElt z) const
{
  RankFlags result;
  for (weyl::Generator s=0; s<rank(); ++s)
    result.set(s,DescentStatus::isDescent(descentValue(s,z)));
  return result;
}

containers::simple_list<BlockElt> down_set(const Block_base& block,BlockElt y)
{
  containers::simple_list<BlockElt> result;

  for (weyl::Generator s : block.descent_generators(y))
    switch (block.descentValue(s,y))
    {
    case DescentStatus::ComplexDescent: result.push_front(block.cross(s,y));
      break;
    case DescentStatus::RealTypeI:
      {
	BlockEltPair sy = block.inverseCayley(s,y);
	result.push_front(sy.first); result.push_front(sy.second);
      }
      break;
    case DescentStatus::RealTypeII:
      result.push_front(block.inverseCayley(s,y).first);
      break;
    default: // |case DescentStatus::ImaginaryCompact| nothing
      break;
    }
  result.sort();
  result.unique();
  return result;

} // |down_set|


/*
  Look up element by |x|, |y| coordinates

  Precondition: |x| and |y| should be compatible: such a block element exists

  This uses the |d_first_z_of_x| table to locate the range where the |x|
  coordinates are correct; then comparing the given |y| value with the first
  one present for |x| (there must be at least one) we can predict the value
  directly, since for each fixed |x| value the values of |y| are consecutive.
*/
BlockElt Block::element(KGBElt xx,KGBElt yy) const
{
  BlockElt first=d_first_z_of_x[xx];
  BlockElt z = first +(yy-y(first));
  assert(z<size() and x(z)==xx and y(z)==yy); // element should be found
  return z;
}

BlockElt Block_base::length_first(size_t l) const
{ // if |length| were an array, we would call |std::lower_bound(begin,end,l)|
  BlockElt min=0, max=size(); // the requested index remains in [min,max]
  while (max>min) // body strictly reduces |max-min| in all cases
  {
    BlockElt z=(min+max)/2;
    if (length(z)>=l)
      max=z; // preserves invariant |length(z)>=l| for all |l>=max|
    else
      min=z+1; // preserves invariant |length(z)<l| for all |l<min|
  }
  assert(min==max);
  return min;
}


/*
  Whether |s| is a strict ascent generator for |z|.

  This means that |descentValue(s,z)| is one of |ComplexAscent|,
  |ImaginaryTypeI| or |ImaginaryTypeII|.
*/
bool Block_base::isStrictAscent(weyl::Generator s, BlockElt z) const
{
  DescentStatus::Value v = descentValue(s,z);
  return not DescentStatus::isDescent(v)
    and v!=DescentStatus::RealNonparity;
}

/*
  Whether |s| is a strict descent generator for |z|.

  This means that |descentValue(s,z)| is one of |ComplexDescent|,
  |RealTypeI| or |RealTypeII|.
*/
bool Block_base::isStrictDescent(weyl::Generator s, BlockElt z) const
{
  DescentStatus::Value v = descentValue(s,z);
  return DescentStatus::isDescent(v)
    and v!=DescentStatus::ImaginaryCompact;
}

/*
  Return the first descent (the number of a simple root) for |z| that is
  not imaginary compact, or |rank()| if there is no such descent.
*/
weyl::Generator Block_base::firstStrictDescent(BlockElt z) const
{
  for (weyl::Generator s = 0; s < rank(); ++s)
    if (isStrictDescent(s,z))
      return s;

  return rank(); // signal nothing was found
}

/*
  Return the first descent (the number of a simple root) for |z| that is either
  complex or real type I; if there is no such descent return |rank()|
*/
weyl::Generator Block_base::firstStrictGoodDescent(BlockElt z) const
{
  for (weyl::Generator s = 0; s < rank(); ++s)
    if (isStrictDescent(s,z) and
	descentValue(s,z)!=DescentStatus::RealTypeII)
      return s;

  return rank(); // signal nothing was found
}

// translation functor from regular to singular $\gamma$ might kill $J_{reg}$
// this depends on the simple coroots for the integral system that vanish on
// the infinitesimal character $\gamma$, namely they make the element zero if
// they define a complex descent, an imaginary compact or a real parity root
bool Block_base::survives(BlockElt z, RankFlags singular) const
{
  const DescentStatus& desc=descent(z);
  for (RankFlags::iterator it=singular.begin(); it(); ++it)
    if (DescentStatus::isDescent(desc[*it]))
      return false;
  return true; // there are no singular simple coroots that are descents
}

// descend through singular simple coroots and return any survivors that were
// reached; they express singular $I(z)$ as sum of 0 or more surviving $I(z')$
containers::sl_list<BlockElt>
  Block_base::finals_for(BlockElt z, RankFlags singular) const
{
  containers::sl_list<BlockElt> result;
  RankFlags::iterator it;
  do
  {
    const descents::DescentStatus& desc=descent(z);
    for (it=singular.begin(); it(); ++it)
      if (DescentStatus::isDescent(desc[*it]))
      {
	switch (desc[*it])
	{
	case DescentStatus::ImaginaryCompact:
	  return result; // 0
	case DescentStatus::ComplexDescent: z = cross(*it,z);
	  break; // follow descent, no branching
	case DescentStatus::RealTypeII:
	  z=inverseCayley(*it,z).first; break; // follow descent, no branching
	case descents::DescentStatus::RealTypeI:
	  {
	    BlockEltPair iC=inverseCayley(*it,z);
	    result.append(finals_for(iC.first,singular));
	    z = iC.second; // continue with right branch, adding its results
	  }
	  break;
	default: assert(false); // should never happen, but compiler wants it
	}
	break; // restart outer loop if a descent was applied
      } // |if(descent(*it,z)|
  }
  while (it()); // terminate on no-break of inner loop
  result.push_back(z);
  return result;
} // |Block_base::finals_for|


// manipulators

void Block_base::set_Bruhat_covered (BlockElt z, BlockEltList&& covered)
{
  assert(z<size());
#ifndef NDEBUG
  for (auto x : covered)
    assert(x<z);
#endif
  partial_Hasse_diagram.resize(size()); // create empty slots for whole block
  if (partial_Hasse_diagram[z].get()==nullptr)
    partial_Hasse_diagram[z].reset(new BlockEltList(std::move(covered)));
}
// Construct the BruhatOrder. Commit-or-rollback is guaranteed.
void Block_base::fill_Bruhat()
{
  if (d_bruhat.get()==nullptr) // if any order is previously stored, just use it
    d_bruhat.reset // otherwise compute it, maybe using |partial_Hasse_diagram|
      (new BruhatOrder(complete_Hasse_diagram(*this,partial_Hasse_diagram)));
}

// computes and stores the KL polynomials
void Block_base::fill_kl_tab(BlockElt limit,
			     KL_hash_Table* pol_hash, bool verbose)
{
  if (kl_tab_ptr.get()==nullptr) // do this only the first time
    kl_tab_ptr.reset(new kl::KL_table(*this,pol_hash));
  kl_tab_ptr->fill(limit,verbose); // extend tables to contain |last_y|
}

// free function

/*
  The functor $T_{\alpha,\beta}$
  List of any descents of |y| by $\alpha$ for which $\beta$ becomes a descent
  and ascents of |y| by |beta$ for which $\alpha$ becomes an ascent

  Here $\alpha$  and $\beta$ shoudl be adjacent roots, of which $\alpha$ is
  a (weak) descent for |y|, while $\beta$ is an ascent for |y|.

  If this is not the case, a pair of |UndefBlock| elements is returned
*/

BlockEltPair link (weyl::Generator alpha,weyl::Generator beta,
		   const Block_base& block, BlockElt y)
{
  BlockElt result[2]; BlockElt* it = &result[0]; // output iterator

  if (block.isStrictDescent(alpha,y))
  {
    if (block.descentValue(alpha,y)==DescentStatus::ComplexDescent)
    {
      BlockElt y1=block.cross(alpha,y);
      if (block.isWeakDescent(beta,y1))
	*it++=y1;
    }
    else // real parity
    {
      BlockEltPair p=block.inverseCayley(alpha,y);
      if (block.isWeakDescent(beta,p.first))
	*it++=p.first;
      if (block.descentValue(alpha,y)==DescentStatus::DescentStatus::RealTypeI
          and block.isWeakDescent(beta,p.second))
	*it++=p.second;
    }
  }

  if (block.isStrictAscent(beta,y))
  {
    if (block.descentValue(beta,y)==DescentStatus::ComplexAscent)
    {
      BlockElt y1=block.cross(beta,y);
      if (not block.isWeakDescent(alpha,y1))
	*it++=y1;
    }
    else // imaginary noncompact
    {
      BlockEltPair p=block.cayley(beta,y);
      if (not block.isWeakDescent(alpha,p.first))
	*it++=p.first;
      if (block.descentValue(beta,y)== DescentStatus::ImaginaryTypeII
	  and not block.isWeakDescent(alpha,p.second))
	*it++=p.second;
    }
  }

  assert(it<=&result[2]);
  while (it<&result[2]) *it++=UndefBlock;

  return std::make_pair(result[0],result[1]);
} // |link|



/*****************************************************************************

        Chapter II -- Derived classes of the Block_base class

******************************************************************************/

/*****			     Bare_block					****/

Bare_block Bare_block::dual(const Block_base& block)
{ auto rank = block.rank();
  auto size = block.size();
  auto max_len = block.length(block.size()-1);

  Bare_block result(rank,block.max_y()+1,block.max_x()+1);
  result.info.reserve(block.size());

  for (BlockElt z=block.size(); z-->0;)
  { result.info.push_back(EltInfo(block.y(z),block.x(z)));
    result.info.back().length = max_len-block.length(z);
    result.info.back().descent = block.descent(z).dual(rank);
  }

  for (unsigned int i=0; i<rank; ++i)
  {
    auto& dst = result.data[i];
    dst.reserve(block.size());
    for (BlockElt z=block.size(); z-->0;)
    {
      dst.emplace_back();
      dst.back().cross_image = size-1-block.cross(i,z);
      const auto& p = block.any_Cayleys(i,z);
      if (p.first!=UndefBlock)
      { dst.back().Cayley_image.first = size-1-p.first;
        if (p.second!=UndefBlock)
	  dst.back().Cayley_image.second = size-1-p.second;
      }
    }
  }

  result.orbits = block.inner_fold_orbits(); // probably not right
  result.dd = block.Dynkin();

  return result;
}


/*****				Block					****/


Block::Block(const Block& b) // obligatory but in practice unused contruction
  : Block_base(b) // copy
  , tW(b.tW) // share
  , d_Cartan(b.d_Cartan)
  , d_involution(b.d_involution)
  , d_first_z_of_x(b.d_first_z_of_x)
  , d_involutionSupport(b.d_involutionSupport)
{}

// Complete the |Block_base| construction, setting |Block|-specific fields
// The real work is done by |Block_base|, |kgb| methods, and |compute_supports|
Block::Block(const KGB& kgb,const KGB& dual_kgb)
  : Block_base(kgb)
  , tW(kgb.twistedWeylGroup())
  , xrange(kgb.size()), yrange(dual_kgb.size())
  , d_Cartan(), d_involution(), d_first_z_of_x(), d_involutionSupport()
    // these fields are filled below
{
  const TwistedWeylGroup& dual_tW =dual_kgb.twistedWeylGroup();

  std::vector<TwistedInvolution> dual_w; // tabulate bijection |tW->dual_tW|
  dual_w.reserve(kgb.nr_involutions());
  size_t size=0;
  for (unsigned int i=0; i<kgb.nr_involutions(); ++i)
  {
    const TwistedInvolution w = kgb.nth_involution(i);
    dual_w.push_back(dual_involution(w,tW,dual_tW));
    size += kgb.packet_size(w)*dual_kgb.packet_size(dual_w.back());
  }

  info.reserve(size);

  // fill |info|
  for (unsigned int i=0; i<kgb.nr_involutions(); ++i)
  {
    // here is where the fibred product via |dual_w| is built
    const TwistedInvolution w = kgb.nth_involution(i);
    const KGBEltPair x_step = kgb.tauPacket(w);
    const KGBEltPair y_step = dual_kgb.tauPacket(dual_w[i]);

    for (KGBElt x=x_step.first; x<x_step.second; ++x)
      for (KGBElt y=y_step.first; y<y_step.second; ++y)
	info.push_back(EltInfo(x,y,descents(x,y,kgb,dual_kgb),kgb.length(x)));
  } // |for (i)|
  compute_first_zs();

  assert(this->size()==size); // check that |info| has exactly |size| elements

  // Now |element| can be safely called; install cross and Cayley tables

  for (weyl::Generator s = 0; s<rank(); ++s)
  { // the generation below is completely independent for each |s|
    data[s].resize(size);
    for (BlockElt z=0; z<size; ++z)
    {
      data[s][z].cross_image
	= element(kgb.cross(s,x(z)),dual_kgb.cross(s,y(z)));
      switch (descentValue(s,z))
      {
      default: break; // most cases leave |data[s][z].Cayley_image| undefined
      case DescentStatus::ImaginaryTypeII:
	{
	  BlockElt z1=element(kgb.cayley(s,x(z)),
			      dual_kgb.inverseCayley(s,y(z)).second);
	  data[s][z].Cayley_image.second = z1; // double-valued direct Cayley
	  data[s][z1].Cayley_image.first = z; // single-valued inverse Cayley
	}
	// FALL THROUGH
      case DescentStatus::ImaginaryTypeI:
	{
	  BlockElt z0=element(kgb.cayley(s,x(z)),
			      dual_kgb.inverseCayley(s,y(z)).first);
	  data[s][z].Cayley_image.first = z0;
	  // in TypeI, |data[s][z].Cayley_image.second| remains |UndefBlock|
	  first_free_slot(data[s][z0].Cayley_image) = z;
	}
      } // switch
    } // |for (z)|
  } // |for(s)|

  // Continue filling the fields of the |Block| derived class proper
  d_Cartan.reserve(size);
  d_involution.reserve(size);
  for (BlockElt z=0; z<size; ++z)
  {
    KGBElt xx=x(z);
    d_Cartan.push_back(kgb.Cartan_class(xx));
    d_involution.push_back(kgb.involution(xx));
  }

  compute_supports();
} // |Block::Block(kgb,dual_kgb)|

// Construction function for the |Block| class.
// It is a pseudo constructor method that ends calling main contructor
Block Block::build(InnerClass& G, RealFormNbr rf, RealFormNbr drf)
{
  RealReductiveGroup G_R(G,rf);
  InnerClass dG(G,tags::DualTag()); // the dual group
  RealReductiveGroup dG_R(dG,drf);

  KGB kgb     (G_R, common_Cartans(G_R,dG_R));
  KGB dual_kgb(dG_R,common_Cartans(dG_R,G_R));
  return Block(kgb,dual_kgb); // |kgb| and |dual_kgb| disappear afterwards!
}

// Given both real group and dual real group, we can just call main contructor
Block Block::build(RealReductiveGroup& G_R, RealReductiveGroup& dG_R)
{
  auto& kgb = G_R.kgb(); auto& dual_kgb = dG_R.kgb(); // temporaries
  return Block(kgb,dual_kgb);
}

// manipulators

void Block::compute_first_zs() // assumes |x| values weakly increase
{
  d_first_z_of_x.resize(xrange+1);
  KGBElt xx=0;
  d_first_z_of_x[xx]=0; // |d_first_z_of_x[xx]| is smallest |z] with |x(z)>=xx|
  for (BlockElt z=0; z<size(); ++z)
    while (xx<x(z)) // no increment in test: often there should be none at all
      d_first_z_of_x[++xx]=z;

  // now |xx==x(size()-1)|; finish off with a sentinel value(s) |size()|
  do // although the largest |x| should be present, so |x==xrange-1| here
    d_first_z_of_x[++xx]=size(); // we don't not depend on that, and fill out
  while (xx<xrange); // stop after setting |d_first_z_of_x[xrange]=size()|
}

// compute the supports in $S$ of twisted involutions
void Block::compute_supports()
{
  d_involutionSupport.reserve(size()); // its eventual size

  // first compute minimal length cases, probably all for the same involution
  for (BlockElt z=0; z<size() and length(z)==length(0); ++z)
  {
    if (z==0 or involution(z)!=involution(z-1))
    { // compute involution support directly from definition
      RankFlags support;
      WeylWord ww=tW.weylGroup().word(involution(z));
      for (size_t j=0; j<ww.size(); ++j)
	support.set(ww[j]);
      d_involutionSupport.push_back(support);
    }
    else // unchanged involution
      d_involutionSupport.push_back(d_involutionSupport.back()); // duplicate
  }

  // complete involution supports at non-minimal lengths, using previous
  for (BlockElt z=d_involutionSupport.size(); z<size(); ++z)
  {
    weyl::Generator s = firstStrictDescent(z);
    assert (s<rank()); // must find one, as we are no longer at minimal length
    DescentStatus::Value v = descentValue(s,z);
    if (v == DescentStatus::ComplexDescent) // cross link
    { // use value from shorter cross neighbour, setting |s| and |twist(s)|
      d_involutionSupport[z] = d_involutionSupport[cross(s,z)];
      d_involutionSupport[z].set(s);
      d_involutionSupport[z].set(tW.twisted(s));
    }
    else // Real Type I or II
    { // use (some) inverse Cayley transform and set |s|
      d_involutionSupport[z] = d_involutionSupport[inverseCayley(s,z).first];
      d_involutionSupport[z].set(s);
    }
  } // |for(z)|
} // |Block::compute_supports|

//		****	     Nothing else for |Block|		****


// auxiliaries to |common_block| constructor not declared in the header file

// add a new block element to |zz_hash| and (therfore) to |info|
void add_z(block_hash& hash,KGBElt x,KGBElt y)
{
  size_t old_size=hash.size();
  BlockElt z=hash.match(block_elt_entry(x,y)); // constructor sets |length==0|
  assert(z==old_size);
  ndebug_use(old_size); ndebug_use(z);
}

// find already constructed element, to be called during construction
BlockElt find_in(const block_hash& hash,KGBElt x,KGBElt y)
{ return hash.find(block_elt_entry(x,y)); }


  // |common_block| methods

RealReductiveGroup& common_block::real_group() const
  { return rc.real_group(); }
const InnerClass& common_block::inner_class() const
  { return rc.inner_class(); }
const InvolutionTable& common_block::involution_table() const
  { return inner_class().involution_table(); }
const RootDatum& common_block::root_datum() const
  { return rc.root_datum(); }

RankFlags common_block::singular (const RatWeight& gamma) const
{
  RankFlags result;
  for (weyl::Generator s=0; s<rank(); ++s)
    result.set(s,root_datum().coroot(integral_sys.parent_nr_simple(s))
				    .dot(gamma.numerator())==0);
  return result;
}

// find value $\gamma-\lambda$ that the parameter for |z| at |gamma%1| would give
RatWeight common_block::gamma_lambda(BlockElt z) const
{
  return rc.gamma_lambda(z_pool[z].srm(rc,gamma_mod_1));
}

common_block::~common_block() = default;


// the full block constructor is only called on explicit user demand
// it is long because of the need to find elements in all corners

common_block::common_block // full block constructor
  (const Rep_context& rc,
   const repr::StandardReprMod& srm, // not modified, |gamma| used mod $X^*$ only
   BlockElt& entry_element	// set to block element matching input
  )
  : Block_base(rootdata::integrality_rank(rc.root_datum(),srm.gamma_mod1()))
  , rc(rc)
  , gamma_mod_1(srm.gamma_mod1()) // already reduced
  , integral_sys(SubSystem::integral(root_datum(),gamma_mod_1))
  , z_pool()
  , srm_hash(z_pool,5)
  , extended(nullptr) // no extended block initially
  , highest_x() // defined below when we have moved to top of block
  , highest_y() // defined below when generation is complete
  , generated_as_full_block(true)
{
  const InnerClass& ic = inner_class();
  const RootDatum& rd = root_datum();

  const InvolutionTable& i_tab = ic.involution_table();
  const KGB& kgb = rc.kgb();

  Block_base::dd = // integral Dynkin diagram, converted from dual side
    DynkinDiagram(integral_sys.cartanMatrix().transposed());

  const unsigned our_rank = integral_sys.rank();

  repr::common_context ctxt(real_group(),integral_sys);

  y_entry::Pooltype y_pool;
  y_part_hash y_hash(y_pool);
  block_hash xy_hash(info);

  // step 1: initialise |z|
  auto z = srm; // get a working copy

  // step 2: move up |z| toward the most split fiber for the current real form
  {
    weyl::Generator s;
    do
      for(s=0; s<our_rank; ++s)
      {
	std::pair<gradings::Status::Value,bool> stat = ctxt.status(s,z.x());
	if (stat.first==gradings::Status::Complex
	    and not stat.second) // complex ascent
	{
	  z = ctxt.cross(s,z);
	  break;
	}
	else if (stat.first==gradings::Status::ImaginaryNoncompact)
	{
	  z = ctxt.up_Cayley(s,z);
	  break;
	}
	// otherwise try next |s|
      } // |for(s)|
    while(s<our_rank); // loop until no ascents found in |integral_sys|
  }
  y_hash.match(i_tab.pack(rc.y_as_torus_elt(z),kgb.inv_nr(highest_x=z.x())));
  // end of step 2

  // step 3: generate imaginary fiber-orbit of |y|'s (|x| is unaffected)
  using LL = containers::sl_list<containers::sl_list<repr::StandardReprMod> >;
  containers::queue<LL>elements; // involution packets, by |x| outer, |y| inner
  {
    const InvolutionNbr theta = kgb.inv_nr(highest_x);
    // generating reflections are by subsystem real roots for |theta0|
    RootNbrSet pos_real = // as subset of full root datum posroots
      integral_sys.positive_roots() & i_tab.real_roots(theta);
    const RootNbrList generator_roots = rd.simpleBasis(pos_real);
    std::vector<WeylWord> reflect(generator_roots.size());
    for (unsigned i=0; i<generator_roots.size(); ++i)
    {
      auto alpha = // generating root expressed as root for |integral_sys|
	integral_sys.from_parent(generator_roots[i]);
      assert(alpha!=RootNbr(-1)); // renaming to subsystem should work
      reflect[i] = integral_sys.reflectionWord(alpha); // word in integral gen's
    }

    containers::sl_list<repr::StandardReprMod> queue { z };
    elements.emplace(); // create empty involution packet at front
    auto& list = // when popping |queue|, move elts here
      elements.front().emplace_back(); // create sublist for the unique |x| value
    do
    {
      auto zz = queue.front();
      list.splice(list.end(),queue,queue.begin()); // move node
      auto prev = srm_hash.size();
      auto seq = srm_hash.match(repr::Repr_mod_entry(rc,zz));
      if (seq==prev and (seq+1)%1000==0)
	std::cout << "Full block step 3:  " << srm_hash.size() << std::endl;
      for (const auto& w : reflect)
      {
	auto new_z = zz;
	for (auto s : w) // order is irrelevant for a reflection word
	  new_z = ctxt.cross(s,new_z);
	assert(new_z.x()==highest_x); // since we have a real reflection

	auto prev = y_pool.size();
	auto new_y = y_hash.match(i_tab.pack(rc.y_as_torus_elt(new_z),theta));
	if (new_y==prev) // then this was really a new y-value
	  queue.push_back(new_z);
      }
    }
    while (not queue.empty());

    // now insert elements from |y_hash| as first R-packet of block
    for (size_t i=0; i<y_hash.size(); ++i)
      add_z(xy_hash,highest_x,i); // adds element to |info|, setting |length==0|

  } // end of step 3

  // step 4: generate packets for successive involutions
  containers::queue<BlockElt> queue { size() }; // involution packet boundaries

  BitMap x_seen(kgb.size()); // for |x| below |highest_x|, record encounters
  x_seen.insert(highest_x);
  BlockElt next=0; // start at beginning of (growing) list of block elements
  containers::sl_list<LL> bundles;

  do // process involution packet of elements from |next| to |queue.front()|
  { // |next| is constant throughout the loop body, popped from |queue| at end
    const KGBElt first_x = x(next);

    // precompute (reversed) length for anything generated this iteration
    const auto next_length = info[next].length+1;

    const LL& bundle = // keep handle on packet being transferred
      elements.pop_splice_to(bundles,bundles.end()); // transfer bundle

    const unsigned int nr_x = bundle.size();
    const unsigned int nr_y = bundle.front().size();
#ifndef NDEBUG // check regularity of the constructed bundle
    assert((queue.front()-next)==nr_x*nr_y);
    const KGBElt first_y=y(next);
    {
      const InvolutionNbr tau = kgb.inv_nr(first_x);
      auto z=next;
      for (const auto& row : bundle)
      {
	assert(kgb.inv_nr(x(z))==tau);
	assert(row.front().x()==x(z));
	auto it = row.wcbegin();
	for (unsigned j=0; j<nr_y; ++j,++z,++it)
	{
	  assert(y_hash[y(z)].nr==tau); // involutions must match
	  assert(y(z)==first_y+j);   // and |y|s are consecutive
	  assert(y_hash.find(i_tab.pack(rc.y_as_torus_elt(*it),tau))==y(z));
	}
      }
    }
#endif

    for (weyl::Generator s=0; s<our_rank; ++s)
    {
      std::vector<block_fields>& tab_s = data[s];
      tab_s.resize(size()); // ensure enough slots for now

      const repr::StandardReprMod& head = bundle.front().front();
      const bool cross_new_involution =
	not x_seen.isMember(ctxt.cross(s,head).x());
      const bool is_real = // whether |s| is real for this involution packet
	ctxt.status(s,first_x).first==gradings::Status::Real;
      const bool is_type1 = is_real and ctxt.status(s,first_x).second;

      containers::sl_list<KGBElt> cross_ys, Cayley_ys;
      const size_t old_y_size = y_pool.size(); // from here up |y|'s are new
      KGBElt sample_x=UndefKGB; // to be set in case of any |Cayley_ys|

      { // compute values into |cross_ys|,
	// check for existence of real parity. in which case compute |Cayley_ys|

	auto it = bundle.front().cbegin();
	for (unsigned int j=0; j<nr_y; ++j,++it)
	{
	  auto sz = ctxt.cross(s,*it);
	  InvolutionNbr theta = kgb.inv_nr(sz.x());
	  TorusElement t = rc.y_as_torus_elt(sz);
	  cross_ys.push_back(y_hash.match(i_tab.pack(t,theta)));
	  if (is_real and ctxt.is_parity(s,*it)) // then do Cayley instead
	  { // looking only for |y| value, so just one Cayley descent suffices
	    sz = ctxt.down_Cayley(s,*it);
	    theta = kgb.inv_nr(sample_x=sz.x()); // same |sample_x| each time
	    t = rc.y_as_torus_elt(sz);
	    Cayley_ys.push_back(y_hash.match(i_tab.pack(t,theta)));
	  }
	}
      } // compute values |cross_ys|

      { // handle cross actions and descent statuses in all cases
	const auto theta = kgb.inv_nr(ctxt.cross(s,head).x());
	BlockElt cur = next; // start of old involution packet

	if (cross_new_involution)
	{ // add a new involution packet
	  LL packet;
	  for (const auto& row : bundle)
	  {
	    auto& packet_list = packet.emplace_back();
	    auto it = cross_ys.cbegin();
	    for (const auto& z : row)
	    {
	      auto &sz = packet_list.emplace_back(ctxt.cross(s,z));
	      assert(*it==y_hash.find(i_tab.pack(rc.y_as_torus_elt(sz),theta)));
	      ndebug_use(theta);
	      tab_s[cur++].cross_image = info.size();
	      add_z(xy_hash,sz.x(),*it);
	      info.back().length=next_length;
	      srm_hash.match(repr::Repr_mod_entry(rc,sz));
	      ++it;
	    }
	    x_seen.insert(packet_list.front().x());
	  } // |for (row)|
	  elements.push(std::move(packet));
	  queue.push(info.size()); // mark end of a new involution packet
	} // |if (cross_new_involution)|
	else
	  for (const auto& row : bundle)
	  {
	    assert(row.front().x()==x(cur));
	    KGBElt s_cross_x =
	      kgb.cross(integral_sys.reflection(s),row.front().x());
	    for (auto y : cross_ys)
	      tab_s[cur++].cross_image = find_in(xy_hash,s_cross_x,y);
	  }

	cur = next; // back up for setting descent status
	if (is_real) // then status depends on |y|
	{ const auto parity = // supposing |y| says "parity", which type is it?
	      is_type1 ? DescentStatus::RealTypeI : DescentStatus::RealTypeII
	    , nonparity = DescentStatus::RealNonparity;
	  for (const auto& row : bundle) // actually |x| is irrelevant
	    for (const auto& z : row)
	      info[cur++].descent.set(s,
				      ctxt.is_parity(s,z) ? parity : nonparity);
	}
	else
	  for (unsigned i=0; i<nr_x; ++i)
	  {
	    auto stat = ctxt.status(s,x(cur));
	    auto block_status =
	      stat.first == gradings::Status::Complex
	      ? stat.second ? DescentStatus::ComplexDescent
	                    : DescentStatus::ComplexAscent
	      : stat.first == gradings::Status::ImaginaryCompact
	                    ? DescentStatus::ImaginaryCompact
	      : stat.second ? DescentStatus::ImaginaryTypeI
	                    : DescentStatus::ImaginaryTypeII;
	    for (unsigned j=0; j<nr_y; ++j)
	      info[cur++].descent.set(s,block_status);
	  }
      } // done for cross action

      if (not Cayley_ys.empty())
      {
	if (not x_seen.isMember(sample_x))
	{ // we must now extend |info| with elements for the new involution
	  // the |x| values of Cayleys of known elements do not suffice; rather
	  // complete |sample_x| to its subsystem fiber over new involution
	  containers::sl_list<containers::sl_list<repr::StandardReprMod> >
	    packet;

	  RootNbrSet pos_imag = // subsystem positive imaginary roots
	    integral_sys.positive_roots() &
	    i_tab.imaginary_roots(kgb.inv_nr(sample_x));
	  RootNbrList imaginary_generators = rd.simpleBasis(pos_imag);

	  auto to_do = containers::queue<KGBElt> { sample_x };
	  do
	  {
	    KGBElt x = to_do.front(); to_do.pop();
	    if (x_seen.isMember(x))
	      continue;

	    x_seen.insert(x);
	    auto& packet_list = packet.emplace_back();

	    for (auto y = old_y_size; y<y_pool.size(); ++y) // distinct new |y|s
	    {
	      add_z(xy_hash,x,y), info.back().length=next_length;
	      auto& new_srm = packet_list.emplace_back
		(repr::StandardReprMod::build
		 (rc,srm.gamma_mod1(), x,y_pool[y].repr().log_pi(false)));
	      auto prev = srm_hash.size();
	      auto seq = srm_hash.match(repr::Repr_mod_entry(rc,new_srm));
	      if (seq==prev and (seq+1)%100==0)
		std::cout << "Full block Cayleys: " << srm_hash.size() << '\n';
	    }

	    // push any new neighbours of |x| onto |to_do|
	    for (RootNbr alpha : imaginary_generators)
	      to_do.push(kgb.cross(rd.reflectionWord(alpha),x));
	  }
	  while (not to_do.empty());

	  elements.push(std::move(packet));
	  queue.push(info.size()); // mark end of a new involution packet
	  tab_s.resize(info.size()); // ensure the Cayley link slots below exist
	} // |if (not x_seen.isMember(sample_x))|: finished extending |info|

	// it remains to set Cayley links in both directions
	BlockElt cur = next; // start of old involution packet
	for (const auto& row : bundle)
	{
	  auto it = Cayley_ys.cbegin();
	  for (const auto& z : row)
	  {
	    if (ctxt.is_parity(s,z))
	    {
	      auto sz = ctxt.down_Cayley(s,z);
	      BlockElt target = find_in(xy_hash,sz.x(),*it);
	      tab_s[cur].Cayley_image.first = target;
	      first_free_slot(tab_s[target].Cayley_image) = cur;
	      if (is_type1)
	      {
		auto other_x = ctxt.cross(s,sz).x();
		assert(x_seen.isMember(other_x));
		target = find_in(xy_hash,other_x,*it);
		tab_s[cur].Cayley_image.second = target;
		first_free_slot(tab_s[target].Cayley_image) = cur;
	      }
	      ++it; // so |it| runs through |Cayley_ys|
	    }
	    ++cur;
	  } // |for(z : row)|
	  assert(Cayley_ys.at_end(it));
	} // |for(row : bundle)|
      } // |if (not Cayley_ys.empty())|
    } // |for(s)|
  }
  while (next=queue.front(), queue.pop(), not queue.empty());
  // end of step 4

  highest_x = last(x_seen); // to be sure; length need not increase with |x|
  highest_y = y_hash.size()-1; // set highest occurring |y| value, for |ysize|

  std::vector<unsigned int> renumber(y_hash.size());
  {
    using pair_tp = std::pair<unsigned long,unsigned int>;
    containers::sl_list<pair_tp> L;
    auto less = [](const pair_tp& a,const pair_tp& b)->bool
      { return a.first<b.first; };
    unsigned int i=0; auto finish = L.end();
    for (const auto& packet : bundles)
    {
      auto start=finish;
      for (const auto& srm : packet.front()) // only traverse first row of packet
      { // gather their |y_stripped| values
	auto y_strip = repr::Repr_mod_entry(rc,srm).y_stripped();
	L.emplace_back(y_strip,i++);
      }
      L.sort(start,finish = L.end(),less);  // and sort interval by those values
      assert(finish==L.end()); // |L.sort| has modified |finish| to achieve this
    }
    i=0;
    for (const auto& pair : L)
      renumber[pair.second]=i++;
  }
  // now we renumber so that for each involution |y| increases with |y_stripped|
  for (auto& entry : info)
    entry.y=renumber[entry.y]; // makes |y_hash| useless; it is dropped anyway

  { // reverse lengths
    auto max_length=info.back().length;
    for (auto& entry : info)
    { assert(entry.length<=max_length);
      entry.length = max_length-entry.length;
    }
  }

  sort();

  entry_element = lookup(srm); // look up element matching the original input

} // |common_block::common_block|, full block version

// the partial block constructor is used for Bruhat intervals below some element
// the precomputed interval is passed as |elements|; this constructor sorts it

common_block::common_block // partial block constructor
    (const repr::Rep_table& rt,
     const repr::common_context& ctxt,
     containers::sl_list<unsigned long>& elements,
     const RatWeight& gamma_mod_1)
  : Block_base(rootdata::integrality_rank(rt.root_datum(),gamma_mod_1))
  , rc(rt)
  , gamma_mod_1(gamma_mod_1) // already reduced
  , integral_sys(SubSystem::integral(root_datum(),gamma_mod_1))
  , z_pool()
  , srm_hash(z_pool)
  , extended(nullptr) // no extended block initially
  , highest_x(0) // it won't be less than this; increased later
  , highest_y(0) // defined when generation is complete
  , generated_as_full_block(false)
{
  info.reserve(elements.size());
  const auto& kgb = rt.kgb();
  const auto& i_tab = inner_class().involution_table();

  Block_base::dd = // integral Dynkin diagram, converted from dual side
    DynkinDiagram(integral_sys.cartanMatrix().transposed());

  y_entry::Pooltype y_pool;
  y_part_hash y_hash(y_pool);

  { // we first fill |y_hash|, carefully ordering those for a same involution
    using y_tab_type = std::pair<unsigned long,TorusPart>;
    std::vector<containers::sl_list<y_tab_type> > y_table
      (inner_class().involution_table().size());
    // every element of |y_table| pairs a |y_stripped| value and a corresponding
    // |TorusPart|; the former is used for sorting, the latter for |gamma_lambda|
    for (unsigned long elt : elements)
    { const auto& srm = rt.srm(elt);
      const KGBElt x = srm.x();
      if (x>highest_x)
	highest_x=x;
      y_tab_type entry(repr::Repr_mod_entry(rc,srm).y_stripped(),srm.y());
      auto& loc = y_table[kgb.inv_nr(x)];
      auto it = std::find_if_not(loc.cbegin(),loc.cend(),
		[&entry](const y_tab_type& t) { return t.first<entry.first; });
      if (it==loc.end() or entry.first<it->first) // only insert |entry| if new
	loc.insert(it,entry);
    }

    for (InvolutionNbr i_x=y_table.size(); i_x-->0; )
    {
      auto old_size = y_hash.size();
      for (const y_tab_type& entry : y_table[i_x])
      {
	RatWeight gamma_lambda = rt.gamma_lambda(i_x,entry.second,gamma_mod_1);
	TorusElement t = y_values::exp_pi(gamma_lambda);
	y_hash.match(i_tab.pack(t,i_x)); // enter this |y_entry| into |y_hash|
      } // we ensured that that |y| increases with |y_stripped| value in packet
      assert(y_hash.size()==old_size+y_table[i_x].size()); // all |y|'s were new
      ndebug_use(old_size);
      highest_y += y_table[i_x].size();
    }
    assert(y_pool.size()==highest_y);
    -- highest_y; // one less than the number of distinct |y| values
  }

  elements.sort // pre-sort by |x| ensures descents precede when setting lengths
    ([&rt](unsigned long a, unsigned long b)
          { return rt.srm(a).x()<rt.srm(b).x(); }
     );

  for (unsigned long elt : elements)
  { const auto& srm=rt.srm(elt);
    const KGBElt x=srm.x();
    auto y = y_hash.find(i_tab.pack(rt.y_as_torus_elt(srm),kgb.inv_nr(x)));
    assert(y!=y_hash.empty);
    info.emplace_back(x,y); // leave descent status unset and |length==0| for now
    auto prev = srm_hash.size();
    auto seq = srm_hash.match(repr::Repr_mod_entry(rc,srm));
    if (seq==prev and (seq+1)%1000==0)
      std::cout << "Partial block:  " << srm_hash.size()
		<< "    " << rt.nr_blocks() << std::endl;
   }

  // allocate link fields with |UndefBlock| entries
  data.assign(integral_sys.rank(),std::vector<block_fields>(elements.size()));

  assert(info.size()==elements.size());
  auto it = elements.cbegin();
  for (BlockElt i=0; i<info.size(); ++i,++it)
  {
    EltInfo& z = info[i];
    const auto srm_z = rt.srm(*it);
    for (weyl::Generator s=0; s<integral_sys.rank(); ++s)
    {
      auto& tab_s = data[s];
      const auto stat = ctxt.status(s,z.x);
      switch (stat.first)
      {
      case gradings::Status::Complex:
	{
	  z.descent.set(s,stat.second
			  ? DescentStatus::ComplexDescent
			  : DescentStatus::ComplexAscent);
	  if (stat.second) // set links both ways when seeing descent
	  {
	    BlockElt sz = lookup(ctxt.cross(s,srm_z));
	    assert(sz!=UndefBlock);
	    if (length(i)==0) // then this is the first descent for |i|
	      info[i].length=length(sz)+1;
	    else
	      assert(length(i)==length(sz)+1);
	    assert(descentValue(s,sz)==DescentStatus::ComplexAscent);
	    tab_s[i].cross_image = sz;
	    tab_s[sz].cross_image = i;
	  }
	}
	break;
      case gradings::Status::Real:
	{
	  if (ctxt.is_parity(s,srm_z))
	  {
	    const auto srm_sz = ctxt.down_Cayley(s,srm_z);
	    BlockElt sz = lookup(srm_sz);
	    assert(sz!=UndefBlock);
	    if (length(i)==0) // then this is the first descent for |i|
	      info[i].length=length(sz)+1;
	    else
	      assert(length(i)==length(sz)+1);
	    tab_s[i].Cayley_image.first = sz; // first Cayley descent
	    if (stat.second)
	    {
	      z.descent.set(s,DescentStatus::RealTypeI);
	      tab_s[i].cross_image = i;
	      assert(descentValue(s,sz)==DescentStatus::ImaginaryTypeI);
	      tab_s[sz].Cayley_image.first = i; // single-valued Cayley ascent
	      sz = lookup(ctxt.cross(s,srm_sz)); // move to other Cayley descent
	      assert(sz!=UndefBlock);
	      assert(descentValue(s,sz)==DescentStatus::ImaginaryTypeI);
	      tab_s[i].Cayley_image.second = sz; // second Cayley descent
	      tab_s[sz].Cayley_image.first = i; // single-valued Cayley ascent
	    }
	    else
	    {
	      z.descent.set(s,DescentStatus::RealTypeII);
	      assert(descentValue(s,sz)==DescentStatus::ImaginaryTypeII);
	      first_free_slot(tab_s[sz].Cayley_image) = i; // one of two ascents
	      tab_s[i].cross_image = lookup(ctxt.cross(s,srm_z)); // maybe Undef
	    }
	  }
	  else
	  {
	    z.descent.set(s,DescentStatus::RealNonparity);
	    tab_s[i].cross_image = i;
	  }
	}
	break;
      case gradings::Status::ImaginaryCompact:
	{
	  z.descent.set(s,DescentStatus::ImaginaryCompact);
	  tab_s[i].cross_image = i;
	}
	break;
      case gradings::Status::ImaginaryNoncompact:
	{
	  if (stat.second)
	  {
	    z.descent.set(s,DescentStatus::ImaginaryTypeI);
	    tab_s[i].cross_image = lookup(ctxt.cross(s,srm_z)); // maybe Undef
	  }
	  else
	  {
	    z.descent.set(s,DescentStatus::ImaginaryTypeII);
	    tab_s[i].cross_image = i;
	  }
	}
	break;
      }
    } // |for (s)|
  } // |for (i)|

  sort();  // finally sort by length, then |x|, then |y|

} // |common_block::common_block|, partial block version

BlockElt common_block::lookup(const repr::StandardReprMod& srm) const
{ // since |srm_hash.empty==UndefBlock|, we can just say:
  return srm_hash.find(repr::Repr_mod_entry(rc,srm));
}

BlockElt common_block::lookup(KGBElt x, const RatWeight& gamma_lambda) const
{
  return lookup(repr::StandardReprMod::build(rc,gamma_mod_1,x,gamma_lambda));
}

repr::StandardRepr common_block::sr (BlockElt z,const RatWeight& gamma) const
{
  const RatWeight gamma_rho = gamma-rho(root_datum());
  const Weight lambda_rho = gamma_rho.integer_diff<int>(gamma_lambda(z));
  return rc.sr_gamma(x(z),lambda_rho,gamma);
}

ext_gens common_block::fold_orbits(const WeightInvolution& delta) const
{
  return rootdata::fold_orbits(integral_sys.pre_root_datum(),delta);
}

ext_block::ext_block& common_block::extended_block(ext_KL_hash_Table* pol_hash)
{
  if (extended.get()==nullptr)
    extended.reset
      (new ext_block::ext_block(*this,inner_class().distinguished(),pol_hash));
  return *extended;
}

kl::Poly_hash_export common_block::KL_hash(KL_hash_Table* KL_pol_hash)
{
  if (kl_tab_ptr.get()==nullptr) // do this only the first time
    kl_tab_ptr.reset(new kl::KL_table(*this,KL_pol_hash));

  return kl_tab_ptr-> polynomial_hash_table();
}

// integrate an older partial block, with mapping of elements
void common_block::swallow
  (common_block&& sub, const BlockEltList& embed,
   KL_hash_Table* KL_pol_hash, ext_KL_hash_Table* ext_KL_pol_hash)
{
  for (BlockElt z=0; z<sub.size(); ++z)
  {
    auto& covered = std::move(sub).Bruhat_order().Hasse(z);
    for (auto& c : covered)
      c=embed[c]; // translate in place
    set_Bruhat_covered(embed[z],std::move(covered));
  }
  if (sub.kl_tab_ptr!=nullptr)
  {
    auto hash_object = KL_hash(KL_pol_hash); // need this for polynomial look-up
    assert (kl_tab_ptr.get()!=nullptr); // because |KL_hash| built |hash|
    kl_tab_ptr->swallow(std::move(*sub.kl_tab_ptr),embed,hash_object.ref);
  }
  if (sub.extended!=nullptr)
  {
    auto& eblock = extended_block(ext_KL_pol_hash); // get/build extended block
    eblock.swallow(std::move(*sub.extended),embed);
  }
}

void common_block::set_Bruhat
  (containers::sl_list<std::pair<BlockElt,BlockEltList> >&& partial_Hasse)
{
  for (auto& pair : partial_Hasse)
    set_Bruhat_covered(pair.first,std::move(pair.second));
}

// comparison of |info| entries: |length|, |x|, |y_stripped()|
bool elt_info_less (const Block_base::EltInfo& a,const Block_base::EltInfo& b)
{ if (a.length!=b.length)
    return a.length<b.length;
  if (a.x!=b.x)
    return a.x<b.x;
  return a.y<b.y; // this implies comparison between |y_stripped| values
}

void common_block::sort()
{
  Permutation ranks = // standardization permutation, to be used for reordering
    permutations::standardization(info.begin(),info.end(),elt_info_less);

  ranks.permute(info); // permute |info|, ordering them by increasing |value|
  ranks.permute(z_pool);
  srm_hash.reconstruct();

  // now adapt |data| tables, assumed to be already computed
  for (weyl::Generator s=0; s<rank(); ++s)
  {
    std::vector<block_fields>& tab_s = data[s];
    ranks.permute(tab_s); // permute fields of |data[s]|
    for (BlockElt z=0; z<size(); ++z) // and update cross and Cayley links
    {
      BlockElt& sz = tab_s[z].cross_image;
      if (sz!=UndefBlock)
	sz = ranks[sz];
      BlockEltPair& p=tab_s[z].Cayley_image;
      if (p.first!=UndefBlock)
      {
	p.first=ranks[p.first];
	if (p.second!=UndefBlock)
	  p.second=ranks[p.second];
      }
    } // |for z|
  } // |for s|

} // |common_block::sort|


// find element in same common block with |x| and |y| parts twisted by |delta|
// may return |UndefBlock|, but success does not mean any representative of the
// $X^*$ coset for this block is actually |delta|-fixed (no such test done)
BlockElt twisted
  (const blocks::common_block& block, BlockElt z,
   const WeightInvolution& delta)
{
  return block.lookup(block.context().kgb().twisted(block.x(z),delta)
		     ,delta*block.gamma_lambda(z));
}


/*****	    |nblock_help| helper class	****/


nblock_help::nblock_help(RealReductiveGroup& GR, const SubSystem& subsys)
  : kgb(GR.kgb()), rd(subsys.parent_datum()), sub(subsys)
  , i_tab(GR.innerClass().involution_table())
  , dual_m_alpha(), half_alpha()
{
  assert(kgb.rank()==rd.semisimpleRank());
  dual_m_alpha.reserve(kgb.rank());
  half_alpha.reserve(kgb.rank());
  for (weyl::Generator s=0; s<kgb.rank(); ++s)
  {
    dual_m_alpha.push_back(TorusPart(rd.simpleRoot(s)));
    half_alpha.push_back(TorusElement(RatWeight(rd.simpleRoot(s),2),false));
  }
}

void nblock_help::check_y (const TorusElement& t, InvolutionNbr i) const
{
  InvolutionData id = sub.involution_data(i_tab.matrix(i));
  const RootNbrList& rb = id.real_basis();
  for (unsigned i=0; i<rb.size(); ++i)
    assert(t.evaluate_at(rd.coroot(rb[i])).normalize().denominator()==1);
}

void nblock_help::parent_cross_act (weyl::Generator s, nblock_elt& z) const
{
  switch (kgb.status(s,z.x()))
  {
  case gradings::Status::Complex:
    z.yy.reflect(rd,rd.simpleRootNbr(s));
    break;
  case gradings::Status::Real:
    z.yy.reflect(rd,rd.simpleRootNbr(s));
    z.yy += dual_m_alpha[s];
    break;
  default: {} // nothing for imaginary (and hence real for |z.yy|) roots
  }
  z.xx=kgb.cross(s,z.xx);
}

void nblock_help::cross_act_parent_word (const WeylWord& ww, nblock_elt& z)
  const
{
  for (size_t i=ww.size(); i-->0; )
    parent_cross_act(ww[i],z);
}

void nblock_help::cross_act (weyl::Generator s, nblock_elt& z) const
{
  cross_act_parent_word(sub.reflection(s),z);
}

void nblock_help::parent_up_Cayley (weyl::Generator s, nblock_elt& z) const
{
  KGBElt cx=kgb.cayley(s,z.xx); // direct Cayley transform on $x$ side
  if (cx == UndefKGB) // undefined Cayley transform: not imaginary noncompact
    return; // silently ignore, done for use from atlas |Cayley| function
  z.xx = cx;

  /* on $y$ side ensure that |z.yy.evaluate_at(rd.simpleCoroot(s))| is even.
   We must adapt by adding a multiple of |simpleRoot(s)|. This may be a
   half-integer multiple even if the initial evaluation is integer, and due to
   that we cannot ensure that the evaluation of |z.yy| on all roots remains
   integer (it will be on real roots, but roots can change their status). In
   the end, adjustment is by a general rational multiple of |simpleRoot(s)|.
  */
  Rational r = z.yy.evaluate_at(rd.simpleCoroot(s))/=2; // in $\Q/\Z$
  int remainder = r.numerator()%r.denominator(); // negative result is OK here
  if (remainder!=0) // odd
    z.yy-=TorusElement(RatWeight(rd.simpleRoot(s)*remainder,r.denominator()),
		       false);
}

void nblock_help::do_up_Cayley (weyl::Generator s, nblock_elt& z) const
{
  const WeylWord& ww=sub.to_simple(s);
  for (size_t i=ww.size(); i-->0; )
    parent_cross_act(ww[i],z);
  parent_up_Cayley(sub.simple(s),z);
  for (size_t i=0; i<ww.size(); ++i)
    parent_cross_act(ww[i],z);
}

bool nblock_help::is_real_nonparity(weyl::Generator s, nblock_elt z) const
{
  cross_act_parent_word(sub.to_simple(s),z);
  assert(kgb.status(sub.simple(s),z.x())==gradings::Status::Real);
  Rational r = z.yy.evaluate_at(rd.simpleCoroot(sub.simple(s))); // modulo $2\Z$
  assert(r.numerator()%r.denominator()==0); // should be integer: real coroot
  return (r.numerator()/r.denominator())%2!=0; // return whether odd
}

void nblock_help::parent_down_Cayley(weyl::Generator s, nblock_elt& z) const
{
  KGBElt cx=kgb.inverseCayley(s,z.xx).first; // inverse Cayley on $x$ side
  if (cx == UndefKGB) // not a real root, so undefined inverse Cayley
    return; // silently ignore, done for use from atlas |inv_Cayley| function

  // on $y$ side just keep the same dual |TorusElement|, so nothing to do
  // however, for non-parity roots, leave $x$ unchanged as well
  Rational r = z.yy.evaluate_at(rd.simpleCoroot(s)); // modulo $2\Z$
  if (r.numerator()%(2*r.denominator())==0) // then it is a parity root
    z.xx = cx; // move $x$ component of |z|
  // for nonparity roots, leave |z| is unchanged for atlas |inv_Cayley|
}

void nblock_help::do_down_Cayley (weyl::Generator s, nblock_elt& z) const
{
  const WeylWord& ww=sub.to_simple(s);
  for (size_t i=ww.size(); i-->0; )
    parent_cross_act(ww[i],z);
  parent_down_Cayley(sub.simple(s),z);
  for (size_t i=0; i<ww.size(); ++i)
    parent_cross_act(ww[i],z);
}


// this essentially modularly reduces the |y| component by taking fingerprint
y_entry nblock_help::pack_y(const nblock_elt& z) const
{
  InvolutionNbr i = kgb.inv_nr(z.x());
#ifndef NDEBUG
  check_y(z.y(),i);
#endif
  return i_tab.pack(z.y(),i);
}



/*****************************************************************************

        Chapter III -- Functions local to blocks.cpp

******************************************************************************/

namespace {

DescentStatus descents(KGBElt x, KGBElt y,
		       const KGB_base& kgb, const KGB_base& dual_kgb)
{
  DescentStatus result;
  for (size_t s = 0; s < kgb.rank(); ++s)
    if (kgb.status(s,x) == gradings::Status::Complex) // s is complex
      if (kgb.isDescent(s,x))
	result.set(s,DescentStatus::ComplexDescent);
      else
	result.set(s,DescentStatus::ComplexAscent);
    else if (kgb.status(s,x) == gradings::Status::ImaginaryNoncompact)
      if (kgb.cross(s,x)!=x) // type I
	result.set(s,DescentStatus::ImaginaryTypeI);
      else // type II
	result.set(s,DescentStatus::ImaginaryTypeII);
    else if (dual_kgb.status(s,y) == gradings::Status::ImaginaryNoncompact)
      if (dual_kgb.cross(s,y)!=y) // type II
	result.set(s,DescentStatus::RealTypeII);
      else // type I
	result.set(s,DescentStatus::RealTypeI);
    // now s is imaginary compact or real nonparity
    else if (kgb.status(s,x) == gradings::Status::Real)
      result.set(s,DescentStatus::RealNonparity);
    else // now |kgb.status(s,x) == gradings::Status::ImaginaryCompact|
      result.set(s,DescentStatus::ImaginaryCompact);

  return result;
}

/*
  Insert into |hs| the ascents through |s| from elements of |hr|.

  This is a technical function for the Hasse construction, that makes the
  part of the co-atom list for a given element arising from a given descent.
*/
  void insert_ascents(const Block_base& block,
		      const Poset::EltList& hr, // ascents of whom
		      size_t s, // ascents by who
		      std::set<BlockElt>& hs) // output
{
  for (BlockElt z : hr)
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


/*
  Put into |Hasse| the Hasse diagram data for the Bruhat ordering on |block|.

  Explanation: we used the algorithm from Vogan's 1982 Park City notes...
  which contains a bad definition. Later modified to work like |kgb::makeHasse|:
  seek a descent |s| that is complex or type I real. If it exists, use it as in
  kgb. If it doesn't then we're essentially at a split principal series. The
  immediate predecessors of |z| are just the inverse Cayley transforms.
*/
std::vector<Poset::EltList> complete_Hasse_diagram
(const Block_base& block,
 std::vector<std::unique_ptr<BlockEltList> >& partial_Hasse_diagram)
{
  std::vector<Poset::EltList> result;
  result.reserve(block.size());

  for (BlockElt z = 0; z < block.size(); ++z)
    if (z<partial_Hasse_diagram.size() and partial_Hasse_diagram[z]!=nullptr)
      result.push_back(std::move(*partial_Hasse_diagram[z])); // accept provided
    else
    {
      std::set<BlockElt> covered;

      auto s = block.firstStrictGoodDescent(z);
      if (s<block.rank())
	switch (block.descentValue(s,z))
	{
	default: assert(false); break;
	case DescentStatus::ComplexDescent:
	  {
	    BlockElt sz = block.cross(s,z);
	    covered.insert(sz);
	    insert_ascents(block,result[sz],s,covered);
	  }
	  break;
	case DescentStatus::RealTypeI: // inverseCayley(s,z) two-valued
	  {
	    BlockEltPair sz = block.inverseCayley(s,z);
	    covered.insert(sz.first);
	    covered.insert(sz.second);
	    insert_ascents(block,result[sz.first],s,covered);
	  }
	}
      else // now just gather all RealTypeII descents of |z|
	for (weyl::Generator s = 0; s < block.rank(); ++s)
	  if (block.descentValue(s,z)==DescentStatus::RealTypeII)
	    covered.insert(block.inverseCayley(s,z).first);

      result.emplace_back(covered.begin(),covered.end()); // convert set->vector
    } // if stored else compute

  partial_Hasse_diagram.clear(); // remove empty shell
  return result;
} // |complete_Hasse_diagram|

} // |namespace|

/*****************************************************************************

      Chapter IV -- Functions exported from blocks.cpp (declared in blocks.h)

******************************************************************************/


// Return the twisted involution dual to |w|.

/*
  We have $\tau = w.\delta$, with $w$ in the Weyl group $W$, and $\delta$ the
  fundamental involution of $X^*$. We seek the twisted involution $v$ in the
  dually twisted Weyl group |dual_tW| such that $-\tau^t = v.\delta^\vee$.
  Here $\delta^\vee$ is the dual fundamental involution, which acts on $X_*$
  as minus the transpose of the longest involution $w_0.\delta$, where $w_0$
  is the longest element of $W$. This relation can be defined independently of
  being a twisted involution, and leads to a bijection $f: W \to W$ that is
  characterised by $f(e)=w_0$ and $f(s.w)=f(w).dwist(s)$ for any simple
  generator $s$ and $w\in W$, where $dwist$ is the twist of the dual twisted
  Weyl group (one also has $f(w.twist(s))=s.f(w)$ so that $f$ intertwines
  twisted conjugation: $f(s.w.twist(s))=s.f(w).dwist(s)$).

  The implementation below is based directly on the above characterisation, by
  converting |w| into a |WeylWord|, and then starting from $w_0$
  right-multiplying successively by the $dwist$ image of the letters of the
  word taken right-to-left. The only thing assumed common to |tW| and
  |dual_tW| is the \emph{external} numbering of generators (letters in |ww|),
  a minimal requirement without which the notion of dual twisted involution
  would make no sense. Notably the result will be correct (when interpreted
  for |dual_tW|) even if the underlying (untwisted) Weyl groups of |W| and
  |dual_W| should differ in their external-to-internal mappings. If one
  assumes that |W| and |dual_W| share the same underlying Weyl group (as is
  currently the case for an inner class and its dual) then one could
  alternatively say simply

    |return tW.prod(tW.weylGroup().longest(),dual_tW.twisted(W.inverse(w)))|

  or

    |return tW.prod(tW.inverse(tW.twisted(w)),tW.weylGroup().longest())|.

  Note that this would involve implicit conversion of an element of |W| as one
  of |dual_W|.
*/
TwistedInvolution
dual_involution(const TwistedInvolution& w,
		const TwistedWeylGroup& tW,
		const TwistedWeylGroup& dual_tW)
{
  WeylWord ww= tW.word(w);
  TwistedInvolution result = dual_tW.weylGroup().longest();
  for (size_t i=ww.size(); i-->0; )
    dual_tW.mult(result,dual_tW.twisted(ww[i]));
  return result;
}

// given dual blocks, map numbers from first block to their partner in second
// this can only work with |Block| values, since the |element| lookup is used
std::vector<BlockElt> dual_map(const Block& b, const Block& dual_b)
{
  assert(b.size()==dual_b.size());

  std::vector<BlockElt> result(b.size());
  for (BlockElt i=0; i<b.size(); ++i)
    result[i]=dual_b.element(b.y(i),b.x(i));

  return result;
}


BitMap common_Cartans(RealReductiveGroup& GR, RealReductiveGroup& dGR)
{ return GR.Cartan_set() & GR.innerClass().dual_Cartan_set(dGR.realForm()); }

} // |namespace blocks|

} // |namespace atlas|
