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
void Block_base::fill_kl_tab(BlockElt last_y,
			     KL_hash_Table* pol_hash, bool verbose)
{
  if (kl_tab_ptr.get()==nullptr) // do this only the first time
    kl_tab_ptr.reset(new kl::KL_table(*this,pol_hash));
  kl_tab_ptr->fill(last_y,verbose); // extend tables to contain |last_y|
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
      // don't set Hermitian dual; they cannot be deduced from |block| at all
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


/*****	    |nblock_help| helper class for |param_block|	****/


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


// auxiliaries to |param_block| constructor not declared in the header file

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


/* |partial_nblock_help| helper class for partrial |param_block| constructor  */

// a derived class whose main additional method computes a Bruhat order ideal
struct partial_nblock_help : public nblock_help
{
  y_part_hash& y_hash; // will reference a field of |param_block|
  block_hash& zz_hash;  // will reference a field of |param_block|

  std::vector<BlockEltList> predecessors; // a list of predecessors for each z

  partial_nblock_help
    (RealReductiveGroup& GR, const SubSystem& subsys,
     y_part_hash& hy, block_hash& hz)
  : nblock_help(GR,subsys)
  , y_hash(hy)
  , zz_hash(hz)
  , predecessors()
{}

  // accessors
  KGBElt conj_in_x(weyl::Generator s,KGBElt x) const
  { return kgb.cross(sub.to_simple(s),x); }
  KGBElt conj_out_x(KGBElt x,weyl::Generator s) const
  { return kgb.cross(x,sub.to_simple(s)); }

  unsigned int length(BlockElt z_inx) const { return zz_hash[z_inx].length; }

  nblock_elt get(BlockElt z_inx) const // convert block index to |nblock_elt|
  {
    const block_elt_entry& e=zz_hash[z_inx];
    return nblock_elt(e.x,y_hash[e.y].repr());
  }
  BlockElt lookup(const block_elt_entry& z_entry) const
  { return zz_hash.find(z_entry); }

  BlockElt lookup(KGBElt x,KGBElt y) const
  { return lookup(block_elt_entry(x,y)); }

  BlockElt lookup(const nblock_elt& z) const
  { KGBElt y = y_hash.find(pack_y(z));
    if (y==y_hash.empty)
      return zz_hash.empty;
    else return lookup(z.x(),y);
  }

  // manipulator

  // the main method, a recursive function
  // construct elements below |z|, at depth |level| from the root call
  BlockElt nblock_below (const nblock_elt& z, unsigned level);

}; // |partial_nblock_help|

// |nblock_below| extends |y_hash|, and |zz_hash| with |z|, having ensured the
// presence of its predecessors. Returns (new, current max) hash index of |z|
BlockElt
  partial_nblock_help::nblock_below (const nblock_elt& z, unsigned level)
{
  { // check if already known, but don't add to |zz_hash| yet if not
    BlockElt res = lookup(z);
    if (res!=zz_hash.empty)
      return res;
  }

  BlockEltList pred; // will hold list of elements covered by z
  // invariant: |nblock_below| has been called for every element in |pred|

  weyl::Generator s; // a complex or real type 1 descent, if such exists

  // |sz_inx| will take a value depending on |s|, but must survive |break|
  BlockElt sz_inx; // will hold number returned by recursive call
  for (s=0; s<sub.rank(); ++s)
  {
    nblock_elt sz = z; // have a copy ready for modification
    KGBElt conj_x= conj_in_x(s,z.x());
    if (kgb.isComplexDescent(sub.simple(s),conj_x))
    {
      cross_act(s,sz);
      sz_inx = nblock_below(sz,level+1); // recursion
      pred.reserve(predecessors[sz_inx].size()+1); // a rough estimate
      pred.push_back(sz_inx); // certainly |sz| is predecessor of |z|
      break; // we shall add $s$-ascents of |predecessors[sz_inx]| below
    }
    else if (kgb.isDoubleCayleyImage(sub.simple(s),conj_x)) // excludes type 2
    {
      if (not is_real_nonparity(s,z)) // excludes real nonparity
      { // so we now know that |z| has a type 1 real descent at |s|
	do_down_Cayley(s,sz);
	sz_inx = nblock_below(sz,level+1); // recursion
	pred.reserve(predecessors[sz_inx].size()+2); // a rough estimate
	pred.push_back(sz_inx);
	cross_act(s,sz); // get other inverse Cayley image of |z|
	pred.push_back(nblock_below(sz,level+1)); // and include it in |pred|
	break; // we shall add $s$-ascents of |predecessors[sz_inx]| below
      } // |if (real_parity)|
    } // |if (doubleCayleyImage)|

  } // |for (s)|

  // if above loop performed a |break| it found complex or real type I descent
  if (s<sub.rank()) // if so, add |s|-ascents for elements covered by |sz|
    for (BlockElt i=0; i<predecessors[sz_inx].size(); ++i)
    {
      nblock_elt c = get(predecessors[sz_inx][i]); // convert from |BlockElt|
      KGBElt conj_x= conj_in_x(s,c.x());
      switch (kgb.status(sub.simple(s),conj_x))
      {
      case gradings::Status::Real: case gradings::Status::ImaginaryCompact:
	break; // nothing to do without ascent
      case gradings::Status::Complex:
	if (not kgb.isDescent(sub.simple(s),conj_x)) // complex ascent
	{
	  cross_act(s,c);
	  pred.push_back(nblock_below(c,level+1)); // recursion
	} // |if(complex ascent)
	break;
      case gradings::Status::ImaginaryNoncompact:
	{
	  bool type_2 = kgb.cross(sub.simple(s),conj_x)==conj_x;
	  do_up_Cayley(s,c);
	  pred.push_back(nblock_below(c,level+1)); // recursion

	  if (type_2)
	  {
	    cross_act(s,c); // this changes |c| since we are in type 2
	    pred.push_back(nblock_below(c,level+1)); // recursion
	  }
	}
	break;
      } // |switch(status(s,conj_x))|
    } // |for (i)|
  else // the loop on |s| above compled, finding only real type II descents
  { // insert and return those descents
    pred.reserve(sub.rank()); // enough, and probably more than that
    while (s-->0) // we reverse the loop just because it looks cute
    {
      nblock_elt sz = z; // have a copy ready for modification
      KGBElt conj_x= conj_in_x(s,z.x());
      if (kgb.status(sub.simple(s),conj_x)==gradings::Status::Real)
      {
	if (not is_real_nonparity(s,z)) // then it was real type II
	{
	  assert (not kgb.isDoubleCayleyImage(sub.simple(s),conj_x));
	  do_down_Cayley(s,sz);
	  pred.push_back(nblock_below(sz,level+1)); // recursion, ignore descents
	}
      }
    } // |while (s-->0)|
  } // |if (s==sub.rank())|

  // finally we can add |z| to |zz_hash|, after all its Bruhat-predecessors
  assert(zz_hash.size()==predecessors.size());
  block_elt_entry e(z.x(),y_hash.match(pack_y(z)),DescentStatus(),level);
  BlockElt res = zz_hash.match(e); // this is where |info| of block grows
  assert(res==predecessors.size()); // |z| must have been added just now
  predecessors.push_back(std::move(pred)); // list of elements covered by |z|
  return res;
} // |partial_nblock_help::nblock_below|



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
