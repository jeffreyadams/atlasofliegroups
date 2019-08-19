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
#include "repr.h" // use of |StandardRepr| in |param_block| constructor

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
std::vector<Poset::EltList> makeHasse(const Block_base&);


} // |namespace|

/*****************************************************************************

        Chapter I -- The Block_base class

******************************************************************************/

// an auxiliary function:
// we often need to fill the first empty slot of a |BlockEltPair|
inline BlockElt& first_free_slot(BlockEltPair& p)
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
  , d_bruhat(nullptr)
  , klc_ptr(nullptr)
{
} // |Block_base::Block_base|

// an almost trivial constructor used for derived non-integral block types
Block_base::Block_base(unsigned int rank)
  : info(), data(rank), orbits()
  , dd()
  , d_bruhat(nullptr)
  , klc_ptr(nullptr)
{}

Block_base::Block_base(const Block_base& b) // copy constructor
  : info(b.info), data(b.data), orbits(b.orbits)
  , dd(b.dd)
  , d_bruhat(nullptr) // don't care to copy; is empty in |Block::build| anyway
  , klc_ptr(nullptr)  // likewise
{
#ifdef VERBOSE // then show that we're called (does not actually happen)
  std::cerr << "copying a block" << std::endl;
#endif
}

Block_base::~Block_base() { delete d_bruhat; delete klc_ptr; }

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
  Tells if s is a strict ascent generator for z.

  Explanation: this means that descentValue(s,z) is one of ComplexAscent,
  ImaginaryTypeI or ImaginaryTypeII.
*/
bool Block_base::isStrictAscent(weyl::Generator s, BlockElt z) const
{
  DescentStatus::Value v = descentValue(s,z);
  return not DescentStatus::isDescent(v)
    and v!=DescentStatus::RealNonparity;
}

/*
  Tells if s is a strict descent generator for z.

  Explanation: this means that descentValue(s,z) is one of ComplexDescent,
  RealTypeI or RealTypeII.
*/
bool Block_base::isStrictDescent(weyl::Generator s, BlockElt z) const
{
  DescentStatus::Value v = descentValue(s,z);
  return DescentStatus::isDescent(v)
    and v!=DescentStatus::ImaginaryCompact;
}

/*
  Returns the first descent for z (the number of a simple root) that is
  not imaginary compact, or rank() if there is no such descent.
*/
weyl::Generator Block_base::firstStrictDescent(BlockElt z) const
{
  for (weyl::Generator s = 0; s < rank(); ++s)
    if (isStrictDescent(s,z))
      return s;

  return rank(); // signal nothing was found
}

/*
 Returns the first descent for z (the number of a simple root) that is either
 complex or real type I; if there is no such descent returns |rank()|
*/
weyl::Generator Block_base::firstStrictGoodDescent(BlockElt z) const
{
  for (weyl::Generator s = 0; s < rank(); ++s)
    if (isStrictDescent(s,z) and
	descentValue(s,z)!=DescentStatus::RealTypeII)
      return s;

  return rank(); // signal nothing was found
}


void param_block::reverse_length_and_sort(bool full_block)
{
  const unsigned max_length =
    // full block ends compact, while partial blcok starts so
    (full_block ? info.back() : info.front()).length;

  const KGBElt x_lim=realGroup().KGB_size(); // limit for |x| values
  std::vector<unsigned> value(size(),0u); // values to be ranked below

  for (BlockElt i=0; i<size(); ++i)
  { assert(length(i)<=max_length);
    auto new_len = info[i].length = max_length-length(i); // reverse length
    value[i]= new_len*x_lim+x(i); // length has priority over value of |x|
  }

  Permutation ranks = // standardization permutation, to be used for reordering
    permutations::standardization(value,(max_length+1)*x_lim,nullptr);

  ranks.permute(info); // permute |info|, ordering them by increasing |value|

  if (not full_block) // data fields are inserted only later for partial block
    return; // so in that case we are done

  // now adapt |data| tables, assumed to be already computed
  for (weyl::Generator s=0; s<rank(); ++s)
  {
    std::vector<block_fields>& tab_s = data[s];
    ranks.permute(tab_s); // permute fields of |data[s]|
    for (BlockElt z=0; z<size(); ++z) // and update cross and Cayley links
    {
      tab_s[z].cross_image = ranks[tab_s[z].cross_image];
      BlockEltPair& p=tab_s[z].Cayley_image;
      if (p.first!=UndefBlock)
      {
	p.first=ranks[p.first];
	if (p.second!=UndefBlock)
	  p.second=ranks[p.second];
      }
    } // |for z|
  } // |for s|

} // |param_block::reverse_length_and_sort|

// Here is one method not related to block construction
/*
  The functor \f$T_{\alpha,\beta}\f$

  Precondition: alpha and beta are adjacent roots, of which alpha is a (weak)
  descent for y, while beta is not a descent for y.

  In fact if this is not satisfied, we return a pair of UndefBlock elements
*/

BlockEltPair Block_base::link
  (weyl::Generator alpha,weyl::Generator beta,BlockElt y) const
{
  const DescentStatus& desc=descent(y);

  BlockElt result[2]; // written using iterator
  BlockElt* it = &result[0];

  BlockEltPair p=inverseCayley(alpha,y); // used only in real parity case
  switch (desc[alpha])
  {
  case DescentStatus::ComplexDescent:
    {
      BlockElt y1=cross(alpha,y);
      if (isWeakDescent(beta,y1))
	*it++=y1;
      break;
    }
  case DescentStatus::RealTypeI:
    if (isWeakDescent(beta,p.second))
      *it++=p.second;
    // FALL THROUGH
  case DescentStatus::RealTypeII:
    if (isWeakDescent(beta,p.first))
      *it++=p.first;
    break;
  default: {}
  } // switch(desc[alpha])

  p=cayley(beta,y);
  switch (desc[beta])
  {
  case DescentStatus::ComplexAscent:
    {
      BlockElt y1=cross(beta,y);
      if (not isWeakDescent(alpha,y1))
	*it++=y1;
      break;
    }
  case DescentStatus::ImaginaryTypeII:
    if (not isWeakDescent(alpha,p.second))
      *it++=p.second;
    // FALL THROUGH
  case DescentStatus::ImaginaryTypeI:
    if (not isWeakDescent(alpha,p.first))
      *it++=p.first;
    break;
  default: {}
  } // switch(desc[beta])

  assert(it<=&result[2]);
  while (it<&result[2]) *it++=UndefBlock;

  return std::make_pair(result[0],result[1]);
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

/*
  Construct the BruhatOrder.
  It could run out of memory, but Commit-or-rollback is guaranteed.
*/
void Block_base::fillBruhat()
{
  if (d_bruhat==nullptr) // do this only the first time
  {
    std::vector<Poset::EltList> hd = makeHasse(*this);
    d_bruhat = new BruhatOrder(hd); // commit iff new completed without throwing
  }
}

// computes and stores the KL polynomials
void Block_base::fill_klc(BlockElt last_y,bool verbose)
{
  if (klc_ptr==nullptr) // do this only the first time
    klc_ptr=new kl::KLContext(*this);

  klc_ptr->fill(last_y,verbose); // extend tables to contain |last_y|
}


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

  // Complete |Block_base| initialisation by installing the dual links.
  if (dual_kgb.Hermitian_dual(0)!=UndefKGB) // then whole block stable
  {
    orbits = tW.twist_orbits(); // orbits on full Dynkin diagram
    // the following is correct since |dual_kgb| was built with dual twist
    for (BlockElt z=0; z<size; ++z)
      info[z].dual =
	element(kgb.Hermitian_dual(x(z)),dual_kgb.Hermitian_dual(y(z)));
  }

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

  KGB kgb     (G_R, common_Cartans(G_R,dG_R),false);
  KGB dual_kgb(dG_R,common_Cartans(dG_R,G_R),true);
  return Block(kgb,dual_kgb); // |kgb| and |dual_kgb| disappear afterwards!
}

// Given both real group and dual real group, we can just call main contructor
Block Block::build(RealReductiveGroup& G_R, RealReductiveGroup& dG_R)
{ return Block(G_R.kgb(),dG_R.kgb_as_dual()); }

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






RealReductiveGroup& param_block::realGroup() const
  { return rc.realGroup(); }
const InnerClass& param_block::innerClass() const
  { return rc.realGroup().innerClass(); }
const InvolutionTable& param_block::involution_table() const
  { return innerClass().involution_table(); }
const RootDatum& param_block::rootDatum() const
  { return innerClass().rootDatum(); }

StandardRepr param_block::sr(BlockElt z) const
  { return rc.sr_gamma(x(z),lambda_rho(z),gamma()); }

RatWeight param_block::nu(BlockElt z) const
{
  InvolutionNbr i_x = rc.kgb().inv_nr(x(z));
  const WeightInvolution& theta = involution_table().matrix(i_x);
  return RatWeight (gamma().numerator()-theta*gamma().numerator()
		    ,2*gamma().denominator()).normalize();
}

Weight param_block::lambda_rho(BlockElt z) const
{
  auto& i_tab = rc.innerClass().involution_table();
  InvolutionNbr i_x = rc.kgb().inv_nr(x(z));
  const WeightInvolution& theta = i_tab.matrix(i_x);

  // do effort to replace |lambda_rho(x)| by what returns from a |Param|
  return (theta*gr_numer + gr_numer // |(1+theta)*gr_numer|, without matrix dup
	  +i_tab.y_lift(i_x,y_bits[y(z)])*gr_denom
	  )/(2*gr_denom);
}

RatWeight param_block::lambda(BlockElt z) const
{
  return rho(rootDatum())+lambda_rho(z);
}


// translation functor from regular to singular $\gamma$ might kill $J_{reg}$
// this depends on the simple coroots for the integral system that vanish on
// the infinitesimal character $\gamma$, namely they make the element zero if
// they define a complex descent, an imaginary compact or a real parity root
bool param_block::survives(BlockElt z) const
{
  const DescentStatus& desc=descent(z);
  for (RankFlags::iterator it=singular.begin(); it(); ++it)
    if (DescentStatus::isDescent(desc[*it]))
      return false;
  return true; // there are no singular simple coroots that are descents
}

// descend through singular simple coroots and return any survivors that were
// reached; they express singular $I(z)$ as sum of 0 or more surviving $I(z')$
containers::sl_list<BlockElt> param_block::finals_for(BlockElt z) const
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
	    result.append(finals_for(iC.first));
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
} // |param_block::finals_for|

// find already constructed element, to be called during construction
BlockElt find_in(const block_hash& hash,KGBElt x,KGBElt y)
{ return hash.find(block_elt_entry(x,y)); } // used during construction



void param_block::compute_duals
  (const y_part_hash& y_hash,const block_hash& hash,
   const InnerClass& G,const SubSystem& rs)
{
  const WeightInvolution& delta = G.distinguished();

  if (delta*infin_char==infin_char) // for stable blocks compute |delta|-action
  {
    WeylWord dummy; // remains empty; the following only serves to get the
    weyl::Twist twist = rs.parent_twist(delta,dummy); // twist induced in |rs|
    {
      unsigned int size=0;
      for (weyl::Generator s=0; s<rs.rank(); ++s)
	if (twist[s]>=s)
	  ++size;
      orbits.reserve(size);
    }

    // analyse the |twist|-orbits on the Dynkin diagram of |rs|
    for (weyl::Generator s=0; s<rs.rank(); ++s)
      if (twist[s]==s)
	orbits.push_back(ext_gen(s));
      else if (twist[s]>s)
	orbits.push_back(ext_gen(rs.cartan(s,twist[s])==0, s,twist[s]));

    for (BlockElt z=0; z<size(); ++z)
    {
      KGBElt dual_x = rc.kgb().Hermitian_dual(x(z));
      if (dual_x!=UndefKGB)
      {
	assert(y_hash[y(z)].nr==rc.kgb().inv_nr(x(z))); // check coherence
	TorusElement t = y_hash[y(z)].t_rep;
	t.act_by(delta); // twist |t| by |delta|
	const KGBElt y =
	  y_hash.find(involution_table().pack(t,rc.kgb().inv_nr(dual_x)));

	info[z].dual = find_in(hash,dual_x,y);
      }
    }
  }
} // |param_block::compute_duals|

BlockElt param_block::lookup(const StandardRepr& sr) const
{ assert(sr.gamma()==gamma());
  return z_hash.find(param_entry{sr.x(),sr.y()});
}

ext_gens param_block::fold_orbits(const WeightInvolution& delta) const
{
  assert(delta*gamma().numerator()==gamma().numerator());
  const auto sub = integrality_datum(innerClass().rootDatum(),infin_char);
  return rootdata::fold_orbits(sub,delta);
}


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

void nblock_help::check_y(const TorusElement& t, InvolutionNbr i) const
{
  InvolutionData id = sub.involution_data(i_tab.matrix(i));
  const RootNbrList& rb = id.real_basis();
  for (unsigned i=0; i<rb.size(); ++i)
    assert(t.evaluate_at(rd.coroot(rb[i])).normalize().denominator()==1);
}

void nblock_help::parent_cross_act(nblock_elt& z, weyl::Generator s) const
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

void nblock_help::cross_act_parent_word(const WeylWord& ww, nblock_elt& z)
  const
{
  for (size_t i=ww.size(); i-->0; )
    parent_cross_act(z,ww[i]);
}

void nblock_help::cross_act (nblock_elt& z, weyl::Generator s) const
{
  cross_act_parent_word(sub.reflection(s),z);
}

void nblock_help::parent_up_Cayley(nblock_elt& z, weyl::Generator s) const
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

void nblock_help::do_up_Cayley (nblock_elt& z, weyl::Generator s) const
{
  const WeylWord& ww=sub.to_simple(s);
  for (size_t i=ww.size(); i-->0; )
    parent_cross_act(z,ww[i]);
  parent_up_Cayley(z,sub.simple(s));
  for (size_t i=0; i<ww.size(); ++i)
    parent_cross_act(z,ww[i]);
}

bool nblock_help::is_real_nonparity(nblock_elt z, weyl::Generator s) const
{
  cross_act_parent_word(sub.to_simple(s),z);
  assert(kgb.status(sub.simple(s),z.x())==gradings::Status::Real);
  Rational r = z.yy.evaluate_at(rd.simpleCoroot(sub.simple(s))); // modulo $2\Z$
  assert(r.numerator()%r.denominator()==0); // should be integer: real coroot
  return (r.numerator()/r.denominator())%2!=0; // return whether odd
}

void nblock_help::parent_down_Cayley(nblock_elt& z, weyl::Generator s) const
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

void nblock_help::do_down_Cayley (nblock_elt& z, weyl::Generator s) const
{
  const WeylWord& ww=sub.to_simple(s);
  for (size_t i=ww.size(); i-->0; )
    parent_cross_act(z,ww[i]);
  parent_down_Cayley(z,sub.simple(s));
  for (size_t i=0; i<ww.size(); ++i)
    parent_cross_act(z,ww[i]);
}

void nblock_help::twist(nblock_elt& z) const
{
  z.xx = kgb.Hermitian_dual(z.xx);
  z.yy.act_by(i_tab.delta); // apply matrix |i_tab.delta| to |RatWeight z.yy|
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


// add a new block element to |zz_hash| and (therfore) to |info|
void add_z(block_hash& hash,KGBElt x,KGBElt y)
{
  size_t old_size=hash.size();
  BlockElt z=hash.match(block_elt_entry(x,y)); // constructor sets |length==0|
  assert(z==old_size);
  assert(z+1==hash.size());
  ndebug_use(old_size);
  ndebug_use(z);
}

param_block::param_block // full block constructor
  (const Rep_context& rc,
   StandardRepr sr,             // by value; made dominant internally
   BlockElt& entry_element	// set to block element matching input
  )
  : Block_base(rootdata::integrality_rank(rc.rootDatum(),sr.gamma()))
  , rc(rc)
  , infin_char(0) // don't set yet
  , gr_numer() // idem
  , y_bits()
  , z_pool()
  , z_hash(z_pool)
  , gr_denom()
  , highest_x(rc.realGroup().KGB_size()-1)
  , highest_y() // defined when generation is complete
  , singular()
{
  y_entry::Pooltype y_pool;
  y_part_hash y_hash(y_pool); // hash table allows storing |y| parts by index
  // A simple structure to pack a pair of already sequenced numbers (indices
  // into the |info| field for some future block) into a hashable value
  block_hash xy_hash(info);

  const InnerClass& G = innerClass();
  const RootDatum& rd = G.rootDatum();

  const InvolutionTable& i_tab = G.involution_table();
  const KGB& kgb = rc.kgb();

  rc.make_dominant(sr); // make dominant before computing subsystem
  infin_char=sr.gamma(); // now we can set the infinitesimal character
  { const RatWeight gamma_rho = (gamma() - rho(rootDatum())).normalize();
    gr_numer=Weight(gamma_rho.numerator().begin(),gamma_rho.numerator().end());
    gr_denom = gamma_rho.denominator();
  }

  const SubSystem sub = SubSystem::integral(rd,infin_char);
  Block_base::dd =
    DynkinDiagram(sub.cartanMatrix().transposed()); // convert from dual side

  size_t our_rank = sub.rank(); // this is independent of ranks in |GR|
  for (weyl::Generator s=0; s<our_rank; ++s)
    singular.set(s,rd.coroot(sub.parent_nr_simple(s))
		            .dot(infin_char.numerator())==0);

  nblock_help aux(realGroup(),sub);

  // step 1: get |y|, which has $y.t=\exp(\pi\ii(\gamma-\lambda))$ (vG based)
  const KGBElt x_org = sr.x();
  const TorusElement y_org = rc.y_as_torus_elt(sr);
  auto z_start = nblock_elt(x_org,y_org); // working copy

  // step 2: move up toward the most split fiber for the current real form
  { // modify |x| and |y|, ascending for |x|, descending for |y|
    weyl::Generator s;
    do
      for(s=0; s<our_rank; ++s)
      {
	KGBElt xx=kgb.cross(sub.to_simple(s),z_start.x());
	weyl::Generator ss=sub.simple(s);
	if (kgb.isAscent(ss,xx))
	{
	  if (kgb.status(ss,xx)==gradings::Status::Complex)
	    aux.cross_act(z_start,s);
	  else // imaginary noncompact
	  {
	    assert(kgb.status(ss,xx) == gradings::Status::ImaginaryNoncompact);
	    aux.do_up_Cayley(z_start,s);
	  }
	  break;
	} // |if(isAscent)|
      } // |for(s)|
    while(s<our_rank); // loop until no ascents found in |subsys|

    // DON'T assert(x+1==kgb.size()); fails e.g. in complex groups, empty |sub|

    y_hash.match(aux.pack_y(z_start)); // save obtained value for |y|
  } // end of step 2

  // step 3: generate imaginary fiber-orbit of |y|'s (|x| is unaffected)
  {
    const KGBElt x0 = z_start.x();
    const InvolutionNbr theta0 = kgb.inv_nr(x0);

    // generating reflections are by subsystem real roots for |theta|
    RootNbrSet pos_real = sub.positive_roots() & i_tab.real_roots(theta0);
    RootNbrList gen_root = rd.simpleBasis(pos_real);
    for (size_t i=0; i<y_hash.size(); ++i) // |y_hash| grows
    {
      const TorusElement y = y_hash[i].repr();
      for (weyl::Generator s=0; s<gen_root.size(); ++s)
      {
	nblock_elt z(x0,y);
	aux.cross_act_parent_word(rd.reflectionWord(gen_root[s]),z);
	assert(z.x()==x0);
	y_hash.match(aux.pack_y(z));
      }
    }

    // now insert elements from |yy_hash| as first R-packet of block
    info.reserve(y_hash.size()); // this is lower bound for final size; reserve
    for (weyl::Generator s=0; s<our_rank; ++s)
      data[s].reserve(y_hash.size());

    for (size_t i=0; i<y_hash.size(); ++i)
      add_z(xy_hash,x0,i); // adds information to |info|; we leave |length==0|

  } // end of step 3

  // step 4: generate packets for successive involutions

  containers::queue<BlockElt> queue { size() }; // invol. packet boundaries
  KGBEltList ys; ys.reserve(0x100); // enough for |1<<RANK_MAX|
  KGBEltList cross_ys(ys.size()); cross_ys.reserve(0x100);
  KGBEltList Cayley_ys(ys.size()); Cayley_ys.reserve(0x100);

  BitMap x_seen(kgb.size());
  x_seen.insert(z_start.x()); // the only value of |x| seen so far

  BlockElt next=0;

  do
  { // process involution packet of elements from |next| to |queue.front()|
    const KGBElt first_x = x(next);
    const InvolutionNbr i_theta = kgb.inv_nr(first_x);

    // precompute (reversed) length for anything generated this iteration
    const auto next_length = info[next].length+1;

    ys.clear();
    size_t nr_y = 0;
    // now traverse R-packet of |first_x|, collecting their |y|'s
    for (BlockElt z=next; z<size() and x(z)==first_x; ++z,++nr_y)
    {
      assert(y_hash[y(z)].nr==i_theta); // involution of |y| must match |x|
      ndebug_use(i_theta);
      assert(y(z)==y(next)+nr_y);   // and |y|s are consecutive
      ys.push_back(y(z));           // so |ys| could have been avoided
    }

    assert((queue.front()-next)%nr_y==0); // |x| values in equal-size R-packets
    unsigned int nr_x= (queue.front()-next)/nr_y; // number of R-packets here

    for (weyl::Generator s=0; s<our_rank; ++s)
    {
      std::vector<block_fields>& tab_s = data[s];
      tab_s.resize(size()); // ensure enough slots for now

      unsigned int y_start=y_hash.size(); // new |y|s numbered from here up

      nblock_elt sample(first_x,y_hash[ys[0]].repr());
      aux.cross_act(sample,s);
      bool new_cross = // whether cross action discovers unseen involution
	y_hash.find(aux.pack_y(sample)) == y_hash.empty;
      assert(x_seen.isMember(sample.x()) == not new_cross);
      // cross action gives a fresh |x| coordinate iff it gives fresh |y|

      cross_ys.clear(); Cayley_ys.clear();
      { // compute values into |cross_ys|
	size_t old_size =  y_hash.size();
	for (unsigned int j=0; j<nr_y; ++j)
	{
	  nblock_elt z(first_x,y_hash[ys[j]].repr());
	  aux.cross_act(z,s);
	  cross_ys.push_back(y_hash.match(aux.pack_y(z)));
	  assert(y_hash.size()== (new_cross ? old_size+j+1 : old_size));
	  ndebug_use(old_size);
	}
      } // compute values |cross_ys|

      for (unsigned int i=0; i<nr_x; ++i) // |i| counts R-packets for involution
      {
	BlockElt base_z = next+i*nr_y; // first element of R-packet
	KGBElt n = x(base_z); // |n| is |x| value for R-packet
	KGBElt s_x_n = kgb.cross(sub.reflection(s),n); // don't need |y| here
	assert (x_seen.isMember(s_x_n) == (not new_cross));
	// cross action on |n| gives a fresh value iff it did for the |ys|

	// set the cross links for this |n| and all corresponding |ys|
	if (new_cross) // add a new R-packet
	{
	  x_seen.insert(s_x_n); // record the new |x| value
	  for (unsigned int j=0; j<nr_y; ++j)
	  {
	    tab_s[base_z+j].cross_image = size(); // link to new element

	    add_z(xy_hash,s_x_n,cross_ys[j]);
	    // same |x| neighbour throughout loop, but |y| neighbour varies

	    info.back().length=next_length;
	  } // |for(j)|
	  // |d_first_z_of_x.push_back(info.size())| would mark end of R-packet
	} // |if(new_cross)|
	else // install cross links to previously existing elements
	  for (unsigned int j=0; j<nr_y; ++j)
	    tab_s[base_z+j].cross_image = find_in(xy_hash,s_x_n,cross_ys[j]);

	// compute component |s| of |info[z].descent|, this |n|, all |y|s
	KGBElt conj_n = kgb.cross(sub.to_simple(s),n); // conjugate
	bool any_Cayley=false; // is made |true| if a real parity case is found
	switch(kgb.status(sub.simple(s),conj_n))
	{
	case gradings::Status::Complex:
	  if (kgb.isDescent(sub.simple(s),conj_n))
	    for (unsigned int j=0; j<nr_y; ++j)
	      info[base_z+j].descent.set(s,DescentStatus::ComplexDescent);
	  else
	    for (unsigned int j=0; j<nr_y; ++j)
	      info[base_z+j].descent.set(s,DescentStatus::ComplexAscent);
	  break;
	case gradings::Status::ImaginaryCompact:
	  for (unsigned int j=0; j<nr_y; ++j)
	    info[base_z+j].descent.set(s,DescentStatus::ImaginaryCompact);
	  break;
	case gradings::Status::ImaginaryNoncompact:
	  if (kgb.cross(sub.simple(s),conj_n)!=conj_n)
	    for (unsigned int j=0; j<nr_y; ++j)
	      info[base_z+j].descent.set(s,DescentStatus::ImaginaryTypeI);
	  else
	    for (unsigned int j=0; j<nr_y; ++j)
	      info[base_z+j].descent.set(s,DescentStatus::ImaginaryTypeII);
	  break;
	case gradings::Status::Real: // now status depends on |y|
	  for (unsigned int j=0; j<nr_y; ++j)
	  {
	    nblock_elt z(n,y_hash[ys[j]].repr()); // element non-conjugated
	    if (aux.is_real_nonparity(z,s))
	      info[base_z+j].descent.set(s,DescentStatus::RealNonparity);
	    else
	    {
	      any_Cayley=true;
	      if (cross_ys[j] != ys[j])
		info[base_z+j].descent.set(s,DescentStatus::RealTypeII);
	      else
		info[base_z+j].descent.set(s,DescentStatus::RealTypeI);
	    }
	  }
	  break;// but this case is re-examined, below
	} // |switch(kgb.status)|


	// now do inverse Cayley transform through |s| if applicable
	if (not any_Cayley)
	  continue; // if there were no valid Cayley links, we're done for |i|

	KGBEltPair Cayleys = kgb.inverseCayley(sub.simple(s),conj_n);
	KGBElt ctx1 = kgb.cross(Cayleys.first,sub.to_simple(s));

	if (i==0)
	{ // |any_Cayley| was independent of |x|; if so, do in first R-packet:
	  assert(y_hash.size()==y_start); // nothing created by cross actions

	  bool new_Cayley = not x_seen.isMember(ctx1);

	  // independent of new creation, fill |Cayley_ys|, extend |y_hash|
	  for (unsigned int j=0; j<nr_y; ++j)
	    if (descentValue(s,base_z+j) != DescentStatus::RealNonparity)
	    {
	      nblock_elt z(n,y_hash[ys[j]].repr());
	      aux.do_down_Cayley(z,s);
	      Cayley_ys.push_back(y_hash.match(aux.pack_y(z)));
	    }

	  if (new_Cayley) // then we need to create new R-packets
	  {
	    // complete fiber of x's over new involution using
	    // subsystem imaginary cross actions
	    RootNbrSet pos_imag = // subsystem positive imaginary roots
	      sub.positive_roots() & i_tab.imaginary_roots(kgb.inv_nr(ctx1));
	    RootNbrList ib = rd.simpleBasis(pos_imag);
	    auto orbit = BitMap{kgb.size()};
	    orbit.insert(ctx1);
	    { auto to_do = containers::queue<KGBElt> { ctx1 };
	      do
	      { KGBElt cur_x = to_do.front(); to_do.pop();
		for (size_t r=0; r<ib.size(); ++r)
		{
		  KGBElt new_x = kgb.cross(rd.reflectionWord(ib[r]),cur_x);
		  if (not orbit.isMember(new_x))
		    orbit.insert(new_x), to_do.push(new_x);
		} // |for(r)|
	      }
	      while (not to_do.empty());
	    }
	    x_seen |= orbit;
	    // then generate corrsponding part of block, combining (x,y)'s
	    for (auto it=orbit.begin(); it(); ++it)
	      for (unsigned int y=y_start; y<y_hash.size(); ++y)
	      {
		add_z(xy_hash,*it,y);
		info.back().length=next_length;
	      }

	    // finallly make sure that Cayley links slots exist for code below
	    tab_s.resize(size());
	  } // |if (new_Cayley)|: finished creating new R-packets
	} // |if (i==0)|: finished work for first |x| when some |y| is parity

	  // now in all cases make links using |element| lookup
	for (unsigned int j=0,p=0; j<nr_y; ++j) // |p| counts parity |y|s only
	  if (descentValue(s,base_z+j)!=DescentStatus::RealNonparity)
	  {
	    KGBElt cty=Cayley_ys[p++]; // unique Cayley transform of |y|
	    BlockElt target = find_in(xy_hash,ctx1,cty);
	    tab_s[base_z+j].Cayley_image.first = target;
	    first_free_slot(tab_s[target].Cayley_image) = base_z+j;
	    if (Cayleys.second!=UndefKGB) // then double valued (type1)
	    {
	      KGBElt ctx2 = kgb.cross(Cayleys.second,sub.to_simple(s));
	      assert (x_seen.isMember(ctx2));
	      target = find_in(xy_hash,ctx2,cty);
	      tab_s[base_z+j].Cayley_image.second = target;
	      first_free_slot(tab_s[target].Cayley_image) = base_z+j;
	    }
	  } // |for(j)|

      } // |for(i)

      if (y_hash.size()>y_start)
	queue.push(size()); // mark end of a new involution packet
    } // |for(s)|
  }
  while (next=queue.front(), queue.pop(), not queue.empty());
  // end of step 4

  highest_y=y_hash.size()-1; // set highest occurring |y| value, for |ysize|
  reverse_length_and_sort(true); // reorder block by increasing value of |x|
  xy_hash.reconstruct(); // adapt to permutation of the block

  compute_duals(y_hash,xy_hash,G,sub); // finally compute Hermitian duals
  compute_y_bits(y_pool);

  // and look up which element matches the original input
  entry_element = lookup(sr);

} // |param_block::param_block|, full block version

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

}; // |class partial_nblock_help|

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
      cross_act(sz,s);
      sz_inx = nblock_below(sz,level+1);
      pred.reserve(predecessors[sz_inx].size()+1); // a rough estimate
      pred.push_back(sz_inx); // certainly |sz| is predecessor of |z|
      break; // we shall add $s$-ascents of |predecessors[sz_inx]| below
    }
    else if (kgb.isDoubleCayleyImage(sub.simple(s),conj_x)) // excludes type 2
    {
      if (not is_real_nonparity(z,s)) // excludes real nonparity
      { // so we now know that |z| has a type 1 real descent at |s|
	do_down_Cayley(sz,s);
	sz_inx = nblock_below(sz,level+1);
	pred.reserve(predecessors[sz_inx].size()+2); // a rough estimate
	pred.push_back(sz_inx);
	cross_act(sz,s); // get other inverse Cayley image of |z|
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
	  cross_act(c,s);
	  pred.push_back(nblock_below(c,level+1));
	} // |if(complex ascent)
	break;
      case gradings::Status::ImaginaryNoncompact:
	{
	  bool type_2 = kgb.cross(sub.simple(s),conj_x)==conj_x;
	  do_up_Cayley(c,s);
	  pred.push_back(nblock_below(c,level+1));

	  if (type_2)
	  {
	    cross_act(c,s); // this changes |c| since we are in type 2
	    pred.push_back(nblock_below(c,level+1));
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
	if (not is_real_nonparity(z,s)) // then it was real type II
	{
	  assert (not kgb.isDoubleCayleyImage(sub.simple(s),conj_x));
	  do_down_Cayley(sz,s);
	  pred.push_back(nblock_below(sz,level+1)); // recurr ignoring descents
	}
      }
    } // |while (s-->0)|
  } // |if (s==sub.rank())|

  // finally we can add |z| to |zz_hash|, after all its Bruhat-predecessors
  assert(zz_hash.size()==predecessors.size());
  block_elt_entry e(z.x(),y_hash.match(pack_y(z)),DescentStatus(),level);
  BlockElt res = zz_hash.match(e);
  assert(res==predecessors.size()); // |z| must have been added just now
  predecessors.push_back(pred); // store list of elements covered by |z|
  return res;
} // |partial_nblock_help::nblock_below|


param_block::param_block // partial block constructor, for interval below |sr|
(const Rep_context& rc, StandardRepr sr) // by value; made dominant internally
  : Block_base(rootdata::integrality_rank(rc.rootDatum(),sr.gamma()))
  , rc(rc)
  , infin_char(0) // don't set yet
  , gr_numer() // idem
  , y_bits()
  , z_pool()
  , z_hash(z_pool)
  , gr_denom()
  , highest_x(0) // it won't be less than this; increased later
  , highest_y() // defined when generation is complete
  , singular() // don't set yet
{
  y_entry::Pooltype y_pool;
  y_part_hash y_hash(y_pool); // hash table allows storing |y| parts by index
  block_hash xy_hash(info);
  const RootDatum& rd = innerClass().rootDatum();

  const KGB& kgb = rc.kgb();

  rc.normalise(sr); // normalise before computing subsystem
  infin_char=sr.gamma(); // now we can set the infinitesimal character
  { const RatWeight gamma_rho = (gamma() - rho(rootDatum())).normalize();
    gr_numer=Weight(gamma_rho.numerator().begin(),gamma_rho.numerator().end());
    gr_denom = gamma_rho.denominator();
  }

  const SubSystem sub = SubSystem::integral(rd,infin_char);
  Block_base::dd =
    DynkinDiagram(sub.cartanMatrix().transposed()); // convert from dual side

  size_t our_rank = sub.rank(); // this is independent of ranks in |GR|
  for (weyl::Generator s=0; s<our_rank; ++s)
    singular.set(s,rd.coroot(sub.parent_nr_simple(s))
			    .dot(infin_char.numerator())==0);

  partial_nblock_help aux(realGroup(),sub,y_hash,xy_hash);
  highest_y=y_hash.size()-1;

  // step 1: get |y|, which has $y.t=\exp(\pi\ii(\gamma-\lambda))$ (vG based)
  const KGBElt x_org = sr.x();
  const nblock_elt org(x_org,rc.y_as_torus_elt(sr));

  // generate partial block in |aux|
  BlockElt last=aux.nblock_below(org,0); // this fills |y_hash| and |xy_hash|

  size_t size= last+1;
  assert(info.size()==size); // |info| should have obtained precisely this size


  reverse_length_and_sort(false); // do reversal operation for partial block
  xy_hash.reconstruct(); // adapt to permutation of the block

  // allocate link fields with |UndefBlock| entries
  data.assign(our_rank,std::vector<block_fields>(size));

  // compute all links for all elements in partial block, by increasing length
  for (BlockElt i=0; i<size; ++i)
  {
    block_elt_entry& z=info[i];
    if (z.x>highest_x) // but first make sure we record highest |x|, for |max_x|
      highest_x=z.x;

    for (weyl::Generator s=0; s<our_rank; ++s)
    {
      auto& tab_s = data[s];
      nblock_elt cur = aux.get(i); // element |z| as |nblock_elt|

      KGBElt conj_x = aux.conj_in_x(s,z.x);
      if (kgb.isDescent(sub.simple(s),conj_x)) // complex descent or real
      {
	if (kgb.status(sub.simple(s),conj_x)==gradings::Status::Complex)
	{
	  aux.cross_act(cur,s);
	  BlockElt sz = aux.lookup(cur);
	  assert(sz!=aux.zz_hash.empty); // should be in generated partial block
	  tab_s[i].cross_image = sz; tab_s[sz].cross_image = i;
	  assert(length(i)==length(sz)+1);
	  z.descent.set(s,DescentStatus::ComplexDescent);
	  assert(descentValue(s,sz)==DescentStatus::ComplexAscent);
	}
	else // |s| is a real root
	{
	  assert(kgb.status(sub.simple(s),conj_x)==gradings::Status::Real);

	  if (aux.is_real_nonparity(cur,s))
	  {
	    tab_s[i].cross_image = i;
	    z.descent.set(s,DescentStatus::RealNonparity);
	  }
	  else // |s| is real parity
	  {
	    aux.do_down_Cayley(cur,s);
	    BlockElt sz = aux.lookup(cur);
	    tab_s[i].Cayley_image.first = sz; // first inverse Cayley
	    assert(length(i)==length(sz)+1);

	    if (kgb.isDoubleCayleyImage(sub.simple(s),conj_x)) // real type 1
	    {
	      z.descent.set(s,DescentStatus::RealTypeI);
	      assert(descentValue(s,sz)==DescentStatus::ImaginaryTypeI);
	      tab_s[i].cross_image = i;
	      tab_s[sz].Cayley_image.first = i; // single-valued Cayley
	      aux.cross_act(cur,s);
	      sz = aux.lookup(cur);
	      assert(descentValue(s,sz)==DescentStatus::ImaginaryTypeI);
	      assert(length(i)==length(sz)+1);
	      tab_s[i].Cayley_image.second = sz; // second inverse Cayley
	      tab_s[sz].Cayley_image.first = i;  // single-valued Cayley
	    }
	    else // real type 2
	    {
	      z.descent.set(s,DescentStatus::RealTypeII);
	      assert(descentValue(s,sz)==DescentStatus::ImaginaryTypeII);
	      first_free_slot(tab_s[sz].Cayley_image) // double-valued Cayley
		= i;
	      cur = aux.get(i); // reset to current element
	      aux.cross_act(cur,s);
	      BlockElt cross_z = aux.lookup(cur);
	      if (cross_z!=aux.zz_hash.empty) // cross neighbour might be absent
		tab_s[i].cross_image = cross_z;
	    } // type 2
	  } // real parity

	} // |s| is real
      } // |if(isDescent)|
      else // ascent, links will be set if and when its ascent appears in block
	if (kgb.status(sub.simple(s),conj_x)==gradings::Status::Complex)
	  z.descent.set(s,DescentStatus::ComplexAscent);
	else // imaginary
	{
	  if (kgb.status(sub.simple(s),conj_x)
	      == gradings::Status::ImaginaryCompact)
	  { z.descent.set(s,DescentStatus::ImaginaryCompact);
	    tab_s[i].cross_image = i;
	  }
	  else if (kgb.cross(sub.simple(s),conj_x)==conj_x)
	  { z.descent.set(s,DescentStatus::ImaginaryTypeII);
	    tab_s[i].cross_image = i;
	  }
	  else // imaginary type 1; here |z| has a nontrivial cross action
	  { z.descent.set(s,DescentStatus::ImaginaryTypeI);
	    KGBElt cross_x = aux.conj_out_x(kgb.cross(sub.simple(s),conj_x),s);
	    BlockElt sz=aux.lookup(cross_x,z.y);
	    if (sz!=aux.zz_hash.empty) // cross neighbour might be absent
	      tab_s[i].cross_image = sz;
	  }
	} // end of imaginary case, and of case distinction

    } // |for(s)|
  } // |for(i)|

  compute_duals(y_hash,xy_hash,innerClass(),sub);
  compute_y_bits(y_pool);

} // |param_block::param_block|, partial block version

void param_block::compute_y_bits(const y_entry::Pooltype& y_pool)
{ y_bits.reserve(y_pool.size());
  const InvolutionTable& i_tab = innerClass().involution_table();
  const RatWeight gamma_rho = gamma() - rho(rootDatum());

  // tabulate some |x| (in fact the first one) for every value |y|
  std::vector<KGBElt> x_of_y(y_pool.size(),UndefKGB);
  for (BlockElt z=0; z<size(); ++z)
    if (x_of_y[y(z)]==UndefKGB)
      x_of_y[y(z)]=x(z);

  for (KGBElt y=0; y<y_pool.size(); ++y)
  { KGBElt x=x_of_y[y];
    assert(x!=UndefKGB); // since every |y| must have at least one matching |x|
    InvolutionNbr i_x = rc.kgb().inv_nr(x);
    RatWeight yr = y_pool[y].repr().log_pi(false);
    yr -= i_tab.matrix(i_x)*yr;
    yr /= 2; // now |yr| has been projected to the $-1$ eigenspace of $\theta$
    const RatWeight lr = (gamma_rho - yr).normalize();
    assert (lr.denominator()==1);
    const Weight lambda_rho(lr.numerator().begin(),lr.numerator().end());
    y_bits.push_back(i_tab.y_pack(i_x,lambda_rho));
  }

  // finally install values into |z_hash|
  z_pool.reserve(size());
  for (BlockElt z=0; z<size(); ++z)
    z_hash.match(param_entry{this->x(z),y_bits[this->y(z)]});
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

  Explanation: technical function for the Hasse construction, that makes the
  part of the coatom list for a given element arising from a given descent.
*/
void insertAscents(std::set<BlockElt>& hs,
		   const Poset::EltList& hr,
		   size_t s,
		   const Block_base& block)
{
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
  Puts into |Hasse| the Hasse diagram data for the Bruhat ordering on |block|.

  Explanation: we used the algorithm from Vogan's 1982 Park City notes...
  which contains a bad definition. Now modified to work like kgb makeHasse:
  seek an ascent s that is complex or type I real. If it exists, use it as in
  kgb. If it doesn't then we're essentially at a split principal series. The
  immediate predecessors of z are just the inverse Cayley transforms.
*/
std::vector<Poset::EltList> makeHasse(const Block_base& block)
{
  std::vector<Poset::EltList> result(block.size());

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
	  insertAscents(h_z,result[sz],s,block);
	}
	break;
      case DescentStatus::RealTypeI: // inverseCayley(s,z) two-valued
	{
	  BlockEltPair sz = block.inverseCayley(s,z);
	  h_z.insert(sz.first);
	  h_z.insert(sz.second);
	  insertAscents(h_z,result[sz.first],s,block);
	}
      }
    else // now just gather all RealTypeII descents of |z|
      for (size_t s = 0; s < block.rank(); ++s)
	if (block.descentValue(s,z)==DescentStatus::RealTypeII)
	  h_z.insert(block.inverseCayley(s,z).first);

    std::copy(h_z.begin(),h_z.end(),std::back_inserter(result[z])); // set->list
  } // for |z|

  return result;
} // |makeHasse|

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
