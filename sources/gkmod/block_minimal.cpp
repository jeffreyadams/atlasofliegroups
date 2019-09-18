/*
  This is block_minimal.cpp

  Copyright (C) 2019 Marc van Leeuwen
  Part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "block_minimal.h"

#include <cassert>
#include <vector>

#include "realredgp.h"
#include "subsystem.h"
#include "kgb.h"
#include "repr.h"

namespace atlas {

namespace blocks {

// Some declarations of undeclared stuff in blocks.cpp that we (continue to) use

void add_z(block_hash& hash,KGBElt x,KGBElt y);
BlockElt find_in(const block_hash& hash,KGBElt x,KGBElt y);
BlockElt& first_free_slot(BlockEltPair& p);

RealReductiveGroup& block_minimal::realGroup() const
  { return rc.realGroup(); }
const InnerClass& block_minimal::innerClass() const
  { return rc.innerClass(); }
const InvolutionTable& block_minimal::involution_table() const
  { return innerClass().involution_table(); }
const RootDatum& block_minimal::rootDatum() const
  { return rc.rootDatum(); }

RatWeight block_minimal::gamma_lambda(BlockElt z) const
{
  return y_pool[y(z)].repr().log_pi(false);
}


block_minimal::block_minimal // full block constructor
  (const Rep_context& rc,
   const StandardRepr sr,       // not modified, |gamma| is mod
   BlockElt& entry_element	// set to block element matching input
  )
  : Block_base(rootdata::integrality_rank(rc.rootDatum(),sr.gamma()))
  , rc(rc)
  , integral_datum(SubSystem::integral(rootDatum(),sr.gamma()))
  , y_pool()
  , y_hash(y_pool)
  , y_part()
  , xy_hash(info)
  , highest_x(rc.realGroup().KGB_size()-1)
  , highest_y() // defined when generation is complete
{
  const InnerClass& ic = innerClass();
  const RootDatum& rd = ic.rootDatum();

  const InvolutionTable& i_tab = ic.involution_table();
  const KGB& kgb = rc.kgb();

  Block_base::dd = // integral Dynkin diagram, converted from dual side
    DynkinDiagram(integral_datum.cartanMatrix().transposed());

  size_t our_rank = integral_datum.rank();

  nblock_help aux(realGroup(),integral_datum);

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
	KGBElt xx=kgb.cross(integral_datum.to_simple(s),z_start.x());
	weyl::Generator ss=integral_datum.simple(s);
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
    while(s<our_rank); // loop until no ascents found in |integral_datum|

    y_hash.match(aux.pack_y(z_start)); // save obtained value for |y|
  } // end of step 2

  // step 3: generate imaginary fiber-orbit of |y|'s (|x| is unaffected)
  {
    const KGBElt x0 = z_start.x();
    const InvolutionNbr theta0 = kgb.inv_nr(x0);

    // generating reflections are by subsystem real roots for |theta0|
    RootNbrSet pos_real =
      integral_datum.positive_roots() & i_tab.real_roots(theta0);
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

  containers::queue<BlockElt> queue { size() }; // involution packet boundaries
  KGBEltList cross_ys; cross_ys.reserve(0x100); // for |1<<RANK_MAX|
  KGBEltList Cayley_ys; Cayley_ys.reserve(0x100);

  BitMap x_seen(kgb.size());
  x_seen.insert(z_start.x()); // the only value of |x| seen so far

  BlockElt next=0; // start at beginning of (growing) list of block elements

  do // process involution packet of elements from |next| to |queue.front()|
  { // |next| is constant throughout the loop body, popped from |queue| at end
    const KGBElt first_x = x(next), first_y=y(next);
    const InvolutionNbr i_theta = kgb.inv_nr(first_x);

    // precompute (reversed) length for anything generated this iteration
    const auto next_length = info[next].length+1;

    unsigned int nr_y = 0; // will be size of the R-packet of |first_x|
    // find it by traversing that R-packet; they have |nr_y| consecutive |y|s
    for (BlockElt z=next; z<size() and x(z)==first_x; ++z)
    {
      assert(y_hash[y(z)].nr==i_theta); // involution of |y| must match |x|
      ndebug_use(i_theta);
      assert(y(z)==first_y+nr_y);   // and |y|s are consecutive
      ++nr_y;
    }

    assert((queue.front()-next)%nr_y==0); // |x| values in equal-size R-packets
    unsigned int nr_x= (queue.front()-next)/nr_y; // number of R-packets here

    for (weyl::Generator s=0; s<our_rank; ++s)
    {
      std::vector<block_fields>& tab_s = data[s];
      tab_s.resize(size()); // ensure enough slots for now

      unsigned int y_start=y_hash.size(); // new |y|s numbered from here up

      nblock_elt sample(first_x,y_hash[first_y].repr());
      aux.cross_act(sample,s);
      bool new_cross = // whether cross action discovers unseen involution
	not x_seen.isMember(sample.x());

      cross_ys.clear(); Cayley_ys.clear();
      { // compute values into |cross_ys|
	size_t old_size =  y_hash.size();
	for (unsigned int j=0; j<nr_y; ++j)
	{
	  nblock_elt z(first_x,y_hash[first_y+j].repr());
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
	KGBElt s_x_n = kgb.cross(integral_datum.reflection(s),n);
	assert (new_cross == not x_seen.isMember(s_x_n)); // check consistency

	// set the cross links for this |n| and all corresponding |y|s
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
	}
	else // if |not new_cross|, install cross links to existing elements
	  for (unsigned int j=0; j<nr_y; ++j)
	    tab_s[base_z+j].cross_image = find_in(xy_hash,s_x_n,cross_ys[j]);

	// compute component |s| of |info[z].descent|, for this |n|, all |y|s
	KGBElt conj_n = kgb.cross(integral_datum.to_simple(s),n); // conjugate
	bool any_Cayley=false; // is made |true| if a real parity case is found
	const auto ids=integral_datum.simple(s);
	switch(kgb.status(ids,conj_n))
	{
	case gradings::Status::Complex:
	  { const auto cplx = kgb.isDescent(ids,conj_n)
	      ? DescentStatus::ComplexDescent
	      : DescentStatus::ComplexDescent;
	    for (unsigned int j=0; j<nr_y; ++j)
	      info[base_z+j].descent.set(s,cplx);
	    break;
	  }
	case gradings::Status::ImaginaryCompact:
	  for (unsigned int j=0; j<nr_y; ++j)
	    info[base_z+j].descent.set(s,DescentStatus::ImaginaryCompact);
	  break;
	case gradings::Status::ImaginaryNoncompact:
	  { const auto imx = kgb.cross(ids,conj_n)!=conj_n
	      ? DescentStatus::ImaginaryTypeI
	      : DescentStatus::ImaginaryTypeII;
	    for (unsigned int j=0; j<nr_y; ++j)
	      info[base_z+j].descent.set(s,imx);
	    break;
	  }
	case gradings::Status::Real: // now status depends on |y|
	  for (unsigned int j=0; j<nr_y; ++j)
	  {
	    nblock_elt z(n,y_hash[first_y+j].repr()); // unconjugated element
	    if (aux.is_real_nonparity(z,s))
	      info[base_z+j].descent.set(s,DescentStatus::RealNonparity);
	    else
	    {
	      any_Cayley=true;
	      if (cross_ys[j] == first_y+j) // whether our own cross neighbour
		info[base_z+j].descent.set(s,DescentStatus::RealTypeI);
	      else
		info[base_z+j].descent.set(s,DescentStatus::RealTypeII);
	    }
	  }
	  break;// but this |Real| case is re-examined, below
	} // |switch(kgb.status(ids,conj_n))|


	// now do inverse Cayley transform through |s| if applicable
	if (not any_Cayley)
	  continue; // if there were no valid Cayley links, we're done for |i|

	KGBEltPair Cayleys = kgb.inverseCayley(ids,conj_n);
	KGBElt ctx1 = kgb.cross(Cayleys.first,integral_datum.to_simple(s));

	if (i==0) // whether in initial R-packet
	{ // |any_Cayley| was independent of |x|; if so, do in first R-packet:
	  assert(y_hash.size()==y_start); // nothing created by cross actions

	  bool new_Cayley = not x_seen.isMember(ctx1);

	  // independent of new creation, fill |Cayley_ys|, extend |y_hash|
	  for (unsigned int j=0; j<nr_y; ++j)
	    if (descentValue(s,base_z+j) != DescentStatus::RealNonparity)
	    {
	      nblock_elt z(n,y_hash[first_y+j].repr());
	      aux.do_down_Cayley(z,s); // anyway there is but a single |y| value
	      Cayley_ys.push_back(y_hash.match(aux.pack_y(z)));
	    }

	  if (new_Cayley) // then we need to create new R-packets
	  {
	    // complete a full (subsystem) fiber of x's over the new involution
	    // using |integral_datum| imaginary cross actions
	    RootNbrSet pos_imag = // subsystem positive imaginary roots
	      integral_datum.positive_roots() &
	      i_tab.imaginary_roots(kgb.inv_nr(ctx1));
	    RootNbrList ib = rd.simpleBasis(pos_imag);

	    {
	      x_seen.insert(ctx1);
	      auto to_do = containers::queue<KGBElt> { ctx1 };
	      do
	      { KGBElt cur_x = to_do.front(); to_do.pop();

		// generate part of block for |cur_x|, and the just added |y|s
		for (unsigned int y=y_start; y<y_hash.size(); ++y)
		{
		  add_z(xy_hash,cur_x,y);
		  info.back().length=next_length;
		}

		// push any new neigghbours of |cur_x| onto |to_do|
		for (size_t r=0; r<ib.size(); ++r)
		{
		  KGBElt new_x = kgb.cross(rd.reflectionWord(ib[r]),cur_x);
		  if (not x_seen.isMember(new_x))
		    x_seen.insert(new_x), to_do.push(new_x);
		} // |for(r)|
	      }
	      while (not to_do.empty());
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
	      KGBElt ctx2 =
		kgb.cross(Cayleys.second,integral_datum.to_simple(s));
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
  reverse_length_and_sort(); // reorder block by increasing value of |x|
  xy_hash.reconstruct(); // adapt to permutation of |info| underlying |xy_hash|


  // and look up which element matches the original input
  entry_element = lookup(sr);

} // |block_minimal::block_minimal|


BlockElt block_minimal::lookup(const StandardRepr& sr) const
{ const auto x = sr.x();
  InvolutionNbr inv = rc.kgb().inv_nr(x);
  const auto y_ent = involution_table().pack(rc.y_as_torus_elt(sr),inv);
  const auto y = y_hash.find(y_ent);
  return y==y_hash.empty ? UndefBlock // the value also known as |xy_hash.empty|
                         : xy_hash.find(EltInfo{x,y});
}
BlockElt block_minimal::lookup(KGBElt x, const RatWeight& gamma_lambda) const
{
  const TorusElement t = y_values::exp_pi(gamma_lambda);
  const auto y_ent = involution_table().pack(t,rc.kgb().inv_nr(x));
  const auto y = y_hash.find(y_ent);
  return y==y_hash.empty ? UndefBlock // the value also known as |xy_hash.empty|
                         : xy_hash.find(EltInfo{x,y});
}

void block_minimal::compute_duals()
{
  const WeightInvolution& delta = innerClass().distinguished();
  weyl::Twist twist;

  {
    WeylWord dummy; // remains empty
    twist = integral_datum.parent_twist(delta,dummy);

    unsigned int size=0;
    for (weyl::Generator s=0; s<integral_datum.rank(); ++s)
      if (twist[s]>=s)
	++size;
    orbits.reserve(size); // dimension |orbits| list in our base object
  }

  // analyse the |twist|-orbits on the Dynkin diagram of |integral_datum|
  for (weyl::Generator s=0; s<integral_datum.rank(); ++s)
    if (twist[s]==s)
      orbits.push_back(ext_gen(s)); // orbit of type |ext_gen::one|
    else if (twist[s]>s) //  orbit of type |ext_gen::two| or |ext_gen::three|
      orbits.push_back
	(ext_gen(integral_datum.cartan(s,twist[s])==0, s,twist[s]));

  const auto& kgb = rc.kgb();
  for (BlockElt z=0; z<size(); ++z)
  {
    KGBElt dual_x = kgb.Hermitian_dual(x(z));
    if (dual_x!=UndefKGB)
    {
      assert(y_hash[y(z)].nr==kgb.inv_nr(x(z))); // check involutions match
      TorusElement t = y_hash[y(z)].t_rep;
      t.act_by(delta); // twist |t| by |delta|
      const KGBElt y =
	y_hash.find(involution_table().pack(t,kgb.inv_nr(dual_x)));

      info[z].dual = find_in(xy_hash,dual_x,y);
    }
  }
} // |block_minimal::compute_duals|

void block_minimal::reverse_length_and_sort()
{
  const unsigned max_length = info.back().length;

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

} // |block_minimal::reverse_length_and_sort|


} // |namespace blocks|

} // |namespace atlas|
