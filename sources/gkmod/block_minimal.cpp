/*
  This is block_minimal.cpp

  Copyright (C) 2019 Marc van Leeuwen
  Part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "block_minimal.h"

#include <cassert>
#include <vector>

#include "matreduc.h"
#include "realredgp.h"
#include "subsystem.h"
#include "kgb.h"
#include "repr.h"
#include "ext_block.h"

namespace atlas {

namespace blocks {

// Some declarations of undeclared stuff in blocks.cpp that we (continue to) use

void add_z(block_hash& hash,KGBElt x,KGBElt y);
BlockElt find_in(const block_hash& hash,KGBElt x,KGBElt y);
BlockElt& first_free_slot(BlockEltPair& p);

  // Function definitions
RealReductiveGroup& block_minimal::realGroup() const
  { return rc.realGroup(); }
const InnerClass& block_minimal::innerClass() const
  { return rc.innerClass(); }
const InvolutionTable& block_minimal::involution_table() const
  { return innerClass().involution_table(); }
const RootDatum& block_minimal::rootDatum() const
  { return rc.rootDatum(); }
const SubSystem& block_minimal::integral_subsystem() const
  { return integral_datum; }

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

    // now insert elements from |y_hash| as first R-packet of block
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
	      : DescentStatus::ComplexAscent;
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

BlockElt twisted
  (const blocks::block_minimal& block, BlockElt z,
   const WeightInvolution& delta)
{
  return block.lookup(block.context().kgb().twisted(block.x(z),delta)
		     ,delta*block.gamma_lambda(z));
}

ext_gens block_minimal::fold_orbits(const WeightInvolution& delta) const
{
  return rootdata::fold_orbits(integral_datum.pre_root_datum(),delta);
}

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

namespace ext_block
{
  // Declarations of some functions re-used from ext_block.cpp

bool in_L_image(Weight beta,WeightInvolution&& A);
bool in_R_image(WeightInvolution&& A,Coweight b);
unsigned int scent_count(DescValue v);
Coweight ell (const KGB& kgb, KGBElt x);

  // Declarations of some local functions
WeylWord fixed_conjugate_simple (const context_minimal& c, RootNbr& alpha);
bool same_standard_reps (const paramin& E, const paramin& F);
bool same_sign (const paramin& E, const paramin& F);
inline bool is_default (const paramin& E)
  { return same_sign(E,paramin(E.ctxt,E.x(),E.gamma_lambda)); }

void z_align (const paramin& E, paramin& F, bool extra_flip);
void z_align (const paramin& E, paramin& F, bool extra_flip, int t_mu);
paramin complex_cross(const ext_gen& p, paramin E);
int level_a (const paramin& E, const Weight& shift, RootNbr alpha);
DescValue star (const paramin& E, const ext_gen& p,
		containers::sl_list<paramin>& links);

bool is_descent (const ext_gen& kappa, const paramin& E);
weyl::Generator first_descent_among
  (RankFlags singular_orbits, const ext_gens& orbits, const paramin& E);


  // |context_minimal| methods and functions

context_minimal::context_minimal
  (const repr::Rep_context& rc, const WeightInvolution& delta,
   const SubSystem& sub)
    : d_rc(rc)
    , d_delta(delta)
    , integr_datum(sub.pre_root_datum())
    , sub(sub)
    , pi_delta(rc.rootDatum().rootPermutation(d_delta))
    , delta_fixed_roots(fixed_points(pi_delta))
    , twist()
    , l_shifts (integr_datum.semisimpleRank())
{
  const RootDatum& rd = rc.rootDatum();
  for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
    twist[s] = rd.simpleRootIndex(delta_of(rd.simpleRootNbr(s)));

  // the reflections for |E.l| pivot around |g_rho_check()|
  const RatCoweight& g_rho_check = this->g_rho_check();
  for (unsigned i=0; i<l_shifts.size(); ++i)
    l_shifts[i] = -g_rho_check.dot(integr_datum.simpleRoot(i));
}

bool context_minimal::is_very_complex (InvolutionNbr theta, RootNbr alpha) const
{ const auto& i_tab = innerClass().involution_table();
  const auto& rd = rootDatum();
  assert (rd.is_posroot(alpha)); // this is a precondition
  auto image = i_tab.root_involution(theta,alpha);
  make_positive(rd,image);
  return image!=alpha and image!=delta_of(alpha);
}

Weight context_minimal::to_simple_shift
  (InvolutionNbr theta, InvolutionNbr theta_p, RootNbrSet S) const
{ const InvolutionTable& i_tab = innerClass().involution_table();
  S &= (i_tab.real_roots(theta) ^i_tab.real_roots(theta_p));
  return root_sum(rootDatum(),S);
}

/*
  For the conjugation to simple scenario, we compute a set of positive roots
  that become negative under an element of $W^\delta$ that makes the
  integrally-simple root(s) in question simple. From this set |S|, and the
  involutions at both ends of the link in the block, the function |shift_flip|
  computes whether an additional flip is to be added to the link.

  This comes from an action of |delta| on a certain top wedge product of
  root spaces, and the formula below tells whether that action is by $-1$.
*/
bool context_minimal::shift_flip
  (InvolutionNbr theta, InvolutionNbr theta_p, RootNbrSet S) const
{ S.andnot(delta_fixed()); // $\delta$-fixed roots won't contribute

  unsigned count=0; // will count 2-element |delta|-orbit elements
  for (auto it=S.begin(); it(); ++it)
    if (is_very_complex(theta,*it) != is_very_complex(theta_p,*it) and
	not rootDatum().sumIsRoot(*it,delta_of(*it)))
      ++count;

  assert(count%2==0); // since |pos_to_neg| is supposed to be $\delta$-stable
  return count%4!=0;
}


/* Try to conjugate |alpha| by product of folded-generators for the (full)
   root system of |c| to a simple root, and return the left-conjugating word
   that was applied. This may fail, if after some conjugation one ends up with
   the long root of a nontrivially folded A2 subsystem (in which case there
   cannot be any solution because |alpha| is fixed by the involution but none
   of the simple roots in its component of the root system is). In this case
   |alpha| is left as that non simple root, and the result conjugates to it.
 */
WeylWord fixed_conjugate_simple (const context_minimal& ctxt, RootNbr& alpha)
{ const RootDatum& rd = ctxt.innerClass().rootDatum();

  WeylWord result;
  while (not rd.is_simple_root(alpha)) // also |break| halfway is possible
  {
    weyl::Generator s = rd.descent_set(alpha)
      .andnot(rd.ascent_set(ctxt.delta_of(alpha))).firstBit();
    assert(s<rd.semisimpleRank()); // exists for positive non-simple roots
    weyl::Generator t = ctxt.twisted(s);
    if (rd.simple_reflected_root(s,alpha)==rd.simpleRootNbr(t))
      break; // |alpha| is sum of (non-commuting) simple roots |s|,|twisted(s)|
    result.push_back(s);
    rd.simple_reflect_root(s,alpha);
    if (s!=t) // second generator for cases of length 2,3
    { result.push_back(t);
      rd.simple_reflect_root(t,alpha);
      if (rd.diagram_linked(s,t)) // then simple reflections |s,t| don't commute
      { // we have $sts=tst$, re-apply |s| to symmetrise the reflection word
	result.push_back(s);
	rd.simple_reflect_root(s,alpha);
      }
    }
  }
  std::reverse(result.begin(),result.end());
  return result;
} // |fixed_conjugate_simple|


  // |paramin| methods and functions


void validate(const paramin& E)
{
#ifndef NDEBUG // make sure this is a no-op when debugging is disabled
  const auto& i_tab = E.rc().innerClass().involution_table();
  const auto& theta = i_tab.matrix(E.tw);
  const auto& delta = E.ctxt.delta();
  assert(delta*theta==theta*delta);
  const Weight symm = E.gamma_lambda.integer_diff<int>(delta*E.gamma_lambda);
  assert(symm == (1-theta)*E.tau);
  assert((delta-1).right_prod(E.l)==(theta+1).right_prod(E.t));
  assert(((E.ctxt.g_rho_check()-E.l)*(1-theta)).numerator().isZero());
  assert(((theta+1)*E.gamma_lambda).numerator().isZero());
#endif
}

paramin::paramin
  (const context_minimal& ec,
   KGBElt x, const RatWeight& gamma_lambda, bool flipped)
  : ctxt(ec)
  , tw(ec.realGroup().kgb().involution(x))
  , l(ell(ec.realGroup().kgb(),x))
  , gamma_lambda(gamma_lambda)
  , tau(matreduc::find_solution
	(1-theta(), gamma_lambda.integer_diff<int>(delta()*gamma_lambda)))
  , t(matreduc::find_solution
	(theta().transposed()+1,(delta()-1).right_prod(l)))
  , flipped(flipped)
{
  validate(*this);
}

paramin::paramin
  (const context_minimal& ec, const TwistedInvolution& tw,
   RatWeight gamma_lambda, Weight tau, Coweight l, Coweight t, bool flipped)
  : ctxt(ec), tw(tw)
  , l(std::move(l))
  , gamma_lambda(std::move(gamma_lambda))
  , tau(std::move(tau))
  , t(std::move(t))
  , flipped(flipped)
{
  validate(*this);
}


// whether |E| and |F| lie over equivalent |StandardRepr| values
bool same_standard_reps (const paramin& E, const paramin& F)
{
  if (&E.ctxt!=&F.ctxt)
  { if (&E.ctxt.innerClass()!=&F.ctxt.innerClass())
      throw std::runtime_error
	("Comparing extended parameters from different inner classes");
    if (E.delta()!=F.delta() or E.ctxt.g_rho_check()!=F.ctxt.g_rho_check())
      return false;
  } // otherwise there might still be a match, so fall through
  return E.theta()==F.theta()
    and in_R_image(E.theta()+1,E.l-F.l)
    and in_L_image(E.gamma_lambda.integer_diff<int>(F.gamma_lambda),E.theta()-1);
}

void z_align (const paramin& E, paramin& F, bool extra_flip)
{ assert(E.t==F.t); // we require preparing |t| upstairs to get this
  int d = E.l.dot((E.delta()-1)*E.tau) - F.l.dot((F.delta()-1)*F.tau);
  assert(d%2==0);
  F.flipped = E.is_flipped()^(d%4!=0)^extra_flip; // XOR 3 Booleans into |F|
}

/*
  In some cases, notably 2i12, one or both of the Cayley transforms may
  involve (in the simple root case) a change |mu| in |lambda_rho| that does
  not necessarily satisfy |t.dot(mu)==0|. In such cases the previous version
  of |z_align| is insufficient, and we should include a contribution from the
  second term of the formula for |z|. But retrieving |mu| from the parameters
  |E| and |F| themselves is complicated by the posssible contribution from
  |Cayley_shift|, which contribution should be ignored; however at the place
  of call the value of |mu| is explicitly available, so we ask here to pass
  |t.dot(mu)| as third argument |t_mu|.
 */
void z_align (const paramin& E, paramin& F, bool extra_flip, int t_mu)
{ z_align(E,F,extra_flip^(t_mu%2!=0)); }


/*
  An auxiliary routine to compute extended parameters across complex links.
  The situation is complicated by the fact that the cross action is by a
  generator of the folded integral system, so we need to exapand it first into
  a product of |length<=3| integral generators, and then have those generators
  act on the components of |E|. For the purpose of changing |E.tw| we further
  develop those generators into reflection words for the full root datum, but
  the reflection action on the other components can be done more directly.

  However, though the integral generators are complex, the action of those
  reflection words need not be purely complex, which implies that the effect
  on the |gamma_lambda| and |l| components are not purely reflections. The
  difference with respect to pure reflection action can be computed comparing
  |rho_r| values (half sums of positive real roots in the full system) with
  the reflected image of that value taken as the starting point, respectively
  (for |l|) the same thing with |rho_check_imaginary|. This is done by the
  "correction" terms below.
 */
paramin complex_cross(const ext_gen& p, paramin E) // by-value for |E|, modified
{ const RootDatum& rd = E.rc().rootDatum();
  const auto& ec = E.ctxt;
  const RootDatum& id = ec.id();
  const InvolutionTable& i_tab = E.rc().innerClass().involution_table();
  auto &tW = E.rc().twistedWeylGroup(); // caution: |p| refers to integr. datum

  InvolutionNbr theta = i_tab.nr(E.tw);
  const RootNbrSet& theta_real_roots = i_tab.real_roots(theta);
  Weight rho_r_shift = rd.twoRho(theta_real_roots);
  Coweight dual_rho_im_shift = rd.dual_twoRho(i_tab.imaginary_roots(theta));

  for (unsigned i=p.w_kappa.size(); i-->0; ) // at most 3 letters, right-to-left
  { weyl::Generator s=p.w_kappa[i]; // generator for integrality datum
    tW.twistedConjugate(ec.subsys().reflection(s),E.tw);
    id.simple_reflect(s,E.gamma_lambda.numerator());
    id.simple_reflect(s,rho_r_shift);
    id.simple_reflect(s,E.tau);
    id.simple_coreflect(E.l,s,ec.l_shift(s));
    id.simple_coreflect(dual_rho_im_shift,s);
    id.simple_coreflect(E.t,s);
  }

  InvolutionNbr new_theta = i_tab.nr(E.tw);
  const RootNbrSet& new_theta_real_roots = i_tab.real_roots(new_theta);
  rho_r_shift -= rd.twoRho(new_theta_real_roots);
  rho_r_shift/=2; // now it is just a sum of (real) roots
  E.gamma_lambda += rho_r_shift;

  assert(ec.delta()*rho_r_shift==rho_r_shift); // diff of $\delta$-fixed

  dual_rho_im_shift -= rd.dual_twoRho(i_tab.imaginary_roots(new_theta));
  dual_rho_im_shift/=2; // now it is just a sum of (imaginary) coroots
  E.l -= dual_rho_im_shift;

  assert(ec.delta().right_prod(dual_rho_im_shift)==dual_rho_im_shift);
  validate(E);

  auto& subs=ec.subsys();
  RootNbr alpha_simple = subs.parent_nr_simple(p.s0);
  const WeylWord to_simple = fixed_conjugate_simple(ec,alpha_simple);
  // by symmetry by $\delta$, |to_simple| conjugates $\delta(\alpha)$ to simple:
  assert(p.length()==1 or rd.is_simple_root(rd.permuted_root(to_simple,
				                subs.parent_nr_simple(p.s1))));
  // apply flip for $\delta$ acting on root set for |to_simple|, as elsewhere
  E.flip(ec.shift_flip(theta,new_theta,pos_to_neg(rd,to_simple)));

  E.flip(p.length()==2); // to parallel the 2i,2r flips

  return E;
} // |complex_cross|


// this implements (comparison using) the formula from Proposition 16 in
// "Parameters for twisted repressentations" (with $\delta-1 = -(1-\delta)$
// the relation is symmetric in |E|, |F|, although not obviously so
bool same_sign (const paramin& E, const paramin& F)
{
  assert(same_standard_reps(E,F));
  const WeightInvolution& delta = E.delta();
  Weight kappa1=E.tau, kappa2=F.tau;
  kappa1 -= delta*kappa1;
  kappa2 -= delta*kappa2;
  int i_exp = E.l.dot(kappa1) - F.l.dot(kappa2);
  assert(i_exp%2==0);
  int n1_exp =
    (F.l-E.l).dot(E.tau) + (E.gamma_lambda-F.gamma_lambda).dot(F.t);
  return ((i_exp/2+n1_exp)%2==0)!=(E.is_flipped()!=F.is_flipped());
}

/*
  for real Cayley transforms, one will add $\rho_r$ to |E.gamma_lambda|
  before projecting it parallel to |alpha| so as to make |alpha_v| vanish on
  |E.gamma_lambda|. Here we compute from |E.gamma_lambda|, corrected by
  that |shift|, the multiple of $\alpha/2$ that such a projection would
  subtract from |E.gamma_lambda|).
*/
int level_a (const paramin& E, const Weight& shift, RootNbr alpha)
{
  const RootDatum& rd = E.rc().rootDatum();
  return (E.gamma_lambda + shift).dot(rd.coroot(alpha));
}


// compute type of |p| for |E|, and export adjacent |paramin| values in |links|
DescValue star (const paramin& E, const ext_gen& p,
		containers::sl_list<paramin>& links)
{
  paramin E0=E; // a copy of |E| that might be modified below to "normalise"
  DescValue result;

  const TwistedWeylGroup& tW = E.rc().twistedWeylGroup();
  const InnerClass& ic = E.rc().innerClass();
  const InvolutionTable& i_tab = ic.involution_table();
  const RootDatum& rd = E.rc().rootDatum();
  const RootDatum& integr_datum = E.ctxt.id();
  const SubSystem& subs = E.ctxt.subsys();
  const InvolutionNbr theta = i_tab.nr(E.tw);
  switch (p.type)
  {
  default: assert(false);
    result=one_complex_ascent; // shut up "maybe uninitialsed" warning
  case ext_gen::one:
    { const Weight& alpha = integr_datum.simpleRoot(p.s0);
      const Coweight& alpha_v = integr_datum.simpleCoroot(p.s0);
      const RootNbr n_alpha = subs.parent_nr_simple(p.s0);
      const RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);

      if (theta_alpha==n_alpha) // length 1 imaginary case
      { // first find out if the simply-integral root $\alpha$ is compact
	int tf_alpha = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	if (tf_alpha%2!=0) // then $\alpha$ is compact
	  return one_imaginary_compact; // quit here, do not collect \$200

	// noncompact case
	const TwistedInvolution new_tw= tW.prod(subs.reflection(p.s0),E.tw);
	const WeightInvolution th_1 = i_tab.matrix(new_tw)-1; // upstairs

	int tau_coef = alpha_v.dot(E.tau); // take $\tau_\alpha$ of table 2

	// try to make $\alpha$ simple by conjugating by $W^\delta$
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const auto theta_p = i_tab.nr(new_tw);
	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = E.ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = E.ctxt.shift_flip(theta,theta_p,S);

	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$
	assert(E.t.dot(alpha)==0); // follows from $\delta*\alpha=\alpha$

	// now separate cases; based on type 1 or 2 first
	if (matreduc::has_solution(th_1,alpha))
	{ // type 1, so extended type is 1i1
	  result = one_imaginary_single;

	  /* if imaginary integrally simple |alpha| were sum of two simple
	     roots these cannot be imaginary (hence integral), so they must be
	     interchanged by $\theta$, but then we are in type 2
	   */
	  assert(rd.is_simple_root(alpha_simple));

	  Weight diff = // called $-\sigma$ in table 2 of [Ptr] (NOTE MINUS)
	      matreduc::find_solution(th_1,alpha); // solutions are equivalent

	  paramin F(E.ctxt,new_tw,
		    E.gamma_lambda- rho_r_shift,
		    E0.tau+diff*tau_coef,
		    E.l+alpha_v*(tf_alpha/2), E.t);

	  E0.l = tf_alpha%4==0 ? F.l+alpha_v : F.l; // for cross
	  assert(not same_standard_reps(E,E0));
	  z_align(E,F,flipped);
	  z_align(F,E0,flipped);

	  links.push_back(std::move(F )); // Cayley link
	  links.push_back(std::move(E0)); // cross link
	} // end of 1i1 case
	else
	{ // imaginary type 2; now we need to distinguish 1i2f and 1i2s

	  auto new_gamma_lambda = E.gamma_lambda; auto new_tau = E.tau;
	  RootNbr first; // maybe a root with |(1-delta)*rd.root(first)==alpha|
	  if (rd.is_simple_root(alpha_simple))
	    first = -1; // invalid value, not used in this case
	  else
	  {
	    --tau_coef; // the parity change and decrease are both relevant
	    weyl::Generator s = // first switched root index
	      rd.find_descent(alpha_simple);
	    first = // corresponding root summand, conjugated back
	      rd.permuted_root(rd.simpleRootNbr(s),ww);
	    assert(alpha == (E.ctxt.delta()+1)*rd.root(first));
	    new_gamma_lambda -= rd.root(first);
	    new_tau -= rd.root(first);
	  }

	  if (tau_coef%2!=0) // was set up so that this means: switched
	  { // no spurious $\tau'$ since $\<\alpha^\vee,(X^*)^\theta>=2\Z$:
	    const auto ratv = (E.ctxt.delta()-1)*(E.gamma_lambda-rho_r_shift);
	    const Weight target
	      { ratv.numerator().begin(),ratv.numerator().end() };
	    assert(not matreduc::has_solution (th_1,target));
	    return one_imaginary_pair_switched; // case 1i2s
	  }
	  result = one_imaginary_pair_fixed;  // what remains is case 1i2f


	  paramin F0(E.ctxt,new_tw,
		     new_gamma_lambda - rho_r_shift,
		     new_tau - alpha*(tau_coef/2),
		     E.l + alpha_v*(tf_alpha/2), E.t);
	  paramin F1(E.ctxt,new_tw,
		     F0.gamma_lambda - alpha, F0.tau, F0.l, E.t);

	  if (not rd.is_simple_root(alpha_simple))
	    flipped = not flipped;

	  z_align(E,F0,flipped);
	  z_align(E,F1,flipped);
	  links.push_back(std::move(F0)); // Cayley link
	  links.push_back(std::move(F1)); // Cayley link
	} // end of type 2 case
      } // end of length 1 imaginary case

      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 1 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const TwistedInvolution new_tw = // downstairs
	  tW.prod(subs.reflection(p.s0),E.tw);

	const auto theta_p=i_tab.nr(new_tw);
	const auto S = pos_to_neg(rd,ww);
	Weight rho_r_shift = E.ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = E.ctxt.shift_flip(theta,theta_p,S);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // as $ww\in W^\delta$

	const int t_alpha = E.t.dot(alpha);
	if (matreduc::has_solution(i_tab.matrix(E.tw)-1,alpha)) // then type 1
	{ // length 1 type 1 real case
	  // for the same reason as |alpha| must be simple in case 1i1, we have:
	  assert(rd.is_simple_root(alpha_simple));
	  const int level = level_a(E,rho_r_shift,n_alpha);
	  if (level%2!=0) // nonparity
	    return one_real_nonparity; // case 1rn, no link added here

	  // now distinguish 1r1f and 1r1s
	  if (t_alpha%2!=0)
	    return one_real_pair_switched; // case 1r1s
	  result = one_real_pair_fixed; // what remains is case 1r1f

	  const RatWeight new_gamma_lambda =
	    E.gamma_lambda + rho_r_shift - alpha*(level/2);
	  assert(new_gamma_lambda.dot(alpha_v)==0); // check effect of |level_a|

	  E0.t -= alpha_v*(t_alpha/2);
	  assert(same_sign(E,E0)); // since only |t| changes

	  paramin F0(E.ctxt,new_tw,new_gamma_lambda, E.tau, E.l          , E0.t);
	  paramin F1(E.ctxt,new_tw,new_gamma_lambda, E.tau, E.l + alpha_v, E0.t);

	  z_align(E0,F0,flipped);
	  z_align(E0,F1,flipped);
	  links.push_back(std::move(F0)); // first Cayley
	  links.push_back(std::move(F1)); // second Cayley

	} // end of 1r1f case
	else // type 2
	{ // length 1 type 2 real

	  int level = level_a(E,rho_r_shift,n_alpha);

	  Weight new_tau=E.tau; // maybe modified below
	  if (not rd.is_simple_root(alpha_simple))
	  { // adapt to integrality based change of lambda
	    RootNbr first = // one of summand roots in |alpha==first+second|
	      rd.permuted_root(rd.simpleRootNbr(rd.find_descent(alpha_simple)),
			       ww);
	    assert(alpha == (E.ctxt.delta()+1)*rd.root(first));
            assert(i_tab.real_roots(theta).isMember(first));

	    rho_r_shift += rd.root(first); // non delta-fixed contribution
	    ++level; // the change in |rho_r_shift| augments its $\alpha$-level

	    // now we must add $d$ to $\tau$ with $(1-\theta')d=(1-\delta)*a0$
	    // where |a0=rd.root(first)|
	    // since $\theta'*a0 = a1 = \delta*a_0$, we can take $d=a0$
	    assert((i_tab.matrix(new_tw))*rd.root(first) ==
		   (E.ctxt.delta()      )*rd.root(first));
	    new_tau += rd.root(first);

	    flipped = not flipped; // flip the Cayley links in this case
	  }
	  if (level%2!=0) // nonparity
	    return one_real_nonparity; // case 1rn, no link added here
	  result = one_real_single; // case 1r2

	  const RatWeight new_gamma_lambda =
	    E.gamma_lambda + rho_r_shift - alpha*(level/2);
	  assert(new_gamma_lambda.dot(alpha_v)==0); // check effect of |level_a|

	  const Coweight diff = // called $s$ in table 2 of [Ptr]
	    matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				    alpha_v);
	  E0.t -= diff*t_alpha;
	  assert(same_sign(E,E0)); // since only |t| changes

	  paramin E1 = E0; // for cross neighbour; share updated value of |t|
	  E1.gamma_lambda -= alpha;
	  assert(not same_standard_reps(E0,E1));

	  paramin F(E.ctxt,new_tw, new_gamma_lambda, new_tau, E.l, E0.t);


	  z_align(E0,F,flipped);
	  z_align(F,E1,flipped);

	  // since |z_align| ignores |gamma_lambda|, we must have equal flips:
	  assert(E0.is_flipped()==E1.is_flipped());

	  links.push_back(std::move(F )); // Cayley link
	  links.push_back(std::move(E1)); // cross link
	} // end of 1r2 case
      }
      else // length 1 complex case
      { result = rd.is_posroot(theta_alpha)
	  ? one_complex_ascent : one_complex_descent ;
	links.push_back(complex_cross(p,E));
      }
    }
    break;

  case ext_gen::two:
    { const Weight& alpha = integr_datum.simpleRoot(p.s0);
      const Coweight& alpha_v = integr_datum.simpleCoroot(p.s0);
      RootNbr n_alpha = subs.parent_nr_simple(p.s0);
      RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);
      const Weight& beta = integr_datum.simpleRoot(p.s1);
      const Coweight& beta_v = integr_datum.simpleCoroot(p.s1);
      RootNbr n_beta = subs.parent_nr_simple(p.s1);
      // RootNbr theta_beta = i_tab.root_involution(theta,n_beta);

      if (theta_alpha==n_alpha) // length 2 imaginary case
      { // first find out if the simply-integral root $\alpha$ is compact
	int tf_alpha = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	int tf_beta = (E.ctxt.g_rho_check() - E.l).dot(beta);
	assert((tf_alpha-tf_beta)%2==0); // same compactness
	if (tf_alpha%2!=0) // then $\alpha$ and $\beta$ are compact
	  return two_imaginary_compact;

	// noncompact case
	const TwistedInvolution new_tw =
	  tW.prod(subs.reflection(p.s1),tW.prod(subs.reflection(p.s0),E.tw));
	// make $\alpha$ simple by conjugating by $W^\delta$
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const auto theta_p = i_tab.nr(new_tw); // upstairs

	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = E.ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = E.ctxt.shift_flip(theta,theta_p,S);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$
	assert(rd.is_simple_root(alpha_simple)); // cannot fail for length 2

	flipped = not flipped; // because of wedge correction for 2i/2r cases

	int at = alpha_v.dot(E.tau); int bt = beta_v.dot(E.tau);
	const WeightInvolution th_1 = i_tab.matrix(new_tw)-1;

	if (matreduc::has_solution(th_1,alpha)) // then type 2i11
	{ result = two_imaginary_single_single;
	  const Weight sigma = matreduc::find_solution(th_1,alpha*at+beta*bt);
	  paramin F (E.ctxt, new_tw,
		     E.gamma_lambda - rho_r_shift,  E.tau + sigma,
		     E.l+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2), E.t);

	  E0.l += alpha_v+beta_v;
	  z_align(E,F,flipped); // no 4th arg, since |E.gamma_lambda| unchanged
	  z_align(F,E0,flipped);
	  links.push_back(std::move(F));  // Cayley link
	  links.push_back(std::move(E0)); // cross link
	}
	else if (matreduc::has_solution(th_1,alpha+beta)) // case 2i12
	{
	  if ((at+bt)%2!=0)
	    return two_imaginary_single_double_switched; // 2i12s
	  result = two_imaginary_single_double_fixed; // 2i12f
	  const int m =  unsigned(at)%2; // safe modular reduction
          const int mm=1-m;

	  // one of the $\tau$ requires upstairs solution for an odd-odd pair:
	  const Weight sigma =
	    matreduc::find_solution(th_1,alpha*(at+mm)+beta*(bt-mm));

	  const Weight new_tau0 = E.tau - alpha*((at+m)/2) - beta*((bt-m)/2);
          const Coweight new_l = E.l+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2);

	  // first Cayley link |F0| will be the one that does not need |sigma|
	  paramin F0(E.ctxt, new_tw,
		     E.gamma_lambda - rho_r_shift - alpha*m, new_tau0,
		     new_l, E.t);
	  paramin F1(E.ctxt, new_tw,
		     E.gamma_lambda - rho_r_shift - alpha*mm, E.tau + sigma,
		     new_l, E.t);

	  int t_alpha=E.t.dot(alpha);
	  z_align(E,F0,flipped,m*t_alpha);
	  z_align(E,F1,flipped,mm*t_alpha);
	  links.push_back(std::move(F0)); // first Cayley
	  links.push_back(std::move(F1)); // second Cayley
	} // end of case 2i12f
	else
	{ // type 2i22
	  result = two_imaginary_double_double;
	  // $\alpha^\vee$ and $\beta^\vee$ are even on $(X^*)^\theta$ and
	  // $(1-\delta)\tau\in(X^*)^\theta+2X^*$ so $<av-bv,\tau>$ is even
	  assert((at-bt)%2==0);
	  int m = static_cast<unsigned int>(at)%2; // safe modular reduction

	  paramin F0(E.ctxt, new_tw,
		     E.gamma_lambda - rho_r_shift + alpha*m,
		     E.tau - alpha*((at+m)/2) - beta*((bt-m)/2),
		     E.l+alpha_v*(tf_alpha/2)+beta_v*(tf_beta/2), E.t);
	  paramin F1(E.ctxt, new_tw,
		     E.gamma_lambda - rho_r_shift - alpha*(1-m) - beta,
		     E.tau - alpha*((at-m)/2) - beta*((bt+m)/2),
		     F0.l,E.t);

	  int ta = E.t.dot(alpha), tb=E.t.dot(beta);
	  z_align(E,F0,flipped,ta*m);
	  z_align(E,F1,flipped,ta*(1-m)+tb);
	  links.push_back(std::move(F0)); // first Cayley
	  links.push_back(std::move(F1)); // second Cayley
	} // end type 2i22 case
      }

      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 2 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	assert(rd.is_simple_root(alpha_simple)); // no complications here
	const TwistedInvolution new_tw = // downstairs
	  tW.prod(subs.reflection(p.s1),tW.prod(subs.reflection(p.s0),E.tw));

	const auto theta_p=i_tab.nr(new_tw);
	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = E.ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = E.ctxt.shift_flip(theta,theta_p,S);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // as $ww\in W^\delta$

	flipped = not flipped; // because of wedge correction for 2i/2r cases

	const int a_level = level_a(E,rho_r_shift,n_alpha);

	if (a_level%2!=0) // nonparity
	   return two_real_nonparity; // no link added here

	const int b_level = level_a(E,rho_r_shift,n_beta);
	assert(b_level%2==0); // since |a_level| and |b_level| have same parity

	WeightInvolution theta_1 = i_tab.matrix(theta)-1; // upstairs

	const RatWeight new_gamma_lambda = E.gamma_lambda + rho_r_shift
	  - alpha*(a_level/2) - beta*(b_level/2);

	int ta = E.t.dot(alpha); int tb = E.t.dot(beta);
	paramin E1=E; // another modifiable copy, like |E0|

	if (matreduc::has_solution(theta_1,alpha))
	{ // type 2r11
	  result = two_real_double_double;
	  // $\alpha$ is even on $(X_*)^{-\theta'}$ (so is $\beta$), and
	  // $t(1-\delta)\in(X_*)^{-\theta'}+2X_*$ so $<t,alpha-beta>$ is even
	  assert((ta-tb)%2==0);
	  int m =  static_cast<unsigned int>(ta)%2;

	  // set two values for |t|; actually the same value in case |m==0|
	  E0.t -= alpha_v*((ta+m)/2) + beta_v*((tb-m)/2);
	  assert(same_sign(E,E0)); // since only |t| changes
	  assert(E0.t.dot(alpha)==-m and E0.t.dot(beta)==m);

	  E1.t -= alpha_v*((ta-m)/2) + beta_v*((tb+m)/2);
	  assert(same_sign(E,E1)); // since only |t| changes
	  assert(E1.t.dot(alpha)==m and E1.t.dot(beta)==-m);

	  paramin F0(E.ctxt, new_tw,
		     new_gamma_lambda,E.tau, E.l+alpha_v*m, E0.t);
	  paramin F1(E.ctxt, new_tw,
		     new_gamma_lambda,E.tau, E.l+alpha_v*(1-m)+beta_v,E1.t);

	  z_align(E0,F0,flipped,m*((b_level-a_level)/2));
	  z_align(E1,F1,flipped,m*((a_level-b_level)/2));

	  // Cayley links
	  links.push_back(std::move(F0));
	  links.push_back(std::move(F1));
	} // end 2r11 case
	else if (matreduc::has_solution(theta_1,alpha+beta))
	{ // type 2r21
	  if ((ta+tb)%2!=0)
	    return two_real_single_double_switched; // 2r21s
	  result = two_real_single_double_fixed; // 2r21f
	  const int m =  static_cast<unsigned int>(ta)%2;
	  const int mm=1-m;

	  // one of the $t$ requires downstairs solution for an odd-odd pair:
	  const Coweight s =
	    matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				    alpha_v*(ta+mm)+beta_v*(tb-mm));

	  // E0 is parameter adapted to Cayley transform that does not need |s|
	  E0.t -= alpha_v*((ta+m)/2) + beta_v*((tb-m)/2);
	  assert(same_sign(E,E0)); // since only |t| changes
	  assert(E0.t.dot(alpha)==-m and E0.t.dot(beta)==m);

	  E1.t -= s;
	  assert(same_sign(E,E1)); // since only |t| changes
	  assert(E1.t.dot(alpha)==-mm and E1.t.dot(beta)==mm);

	  // Cayley links
	  paramin F0(E.ctxt, new_tw,
		     new_gamma_lambda, E.tau, E.l+alpha_v*m, E0.t);
	  paramin F1(E.ctxt, new_tw,
		     new_gamma_lambda, E.tau, E.l+alpha_v*mm, E1.t);

	  z_align(E0,F0,flipped,m *((b_level-a_level)/2));
	  z_align(E1,F1,flipped,mm*((b_level-a_level)/2));
	  links.push_back(std::move(F0));
	  links.push_back(std::move(F1));
	} // end of case 2r21f
	else // case 2r22
	{ result = two_real_single_single;
	  const Coweight s =
	    matreduc::find_solution(i_tab.matrix(new_tw).transposed()+1,
				    alpha_v*ta+beta_v*tb);

	  E0.t -= s; // parameter adapted to Cayley transform |F|
	  assert(same_sign(E,E0)); // since only |t| changes
	  assert(E.t.dot(alpha)==0 and E.t.dot(beta)==0);

	  E1.gamma_lambda -= alpha+beta;
	  E1.t = E0.t; // cross action, keeps adaptation of |t| to |F| below
	  assert(not same_standard_reps(E0,E1));

	  paramin F(E.ctxt, new_tw, new_gamma_lambda, E.tau, E.l, E0.t);

	  z_align(E0,F,flipped); // no 4th arg, as |E.t.dot(alpha)==0| etc.
	  z_align(F,E1,flipped);
	  links.push_back(std::move(F )); // Cayley link
	  links.push_back(std::move(E1)); // cross link
	} // end of case 2r22
      }
      else // length 2 complex case
      { const bool ascent = rd.is_posroot(theta_alpha);
	if (theta_alpha != (ascent ? n_beta : rd.rootMinus(n_beta)))
	{ // twisted non-commutation with |s0.s1|
	  result = ascent ? two_complex_ascent : two_complex_descent;
	  links.push_back(complex_cross(p,E));
	}
	else if (ascent)
	{ // twisted commutation with |s0.s1|: 2Ci
	  result = two_semi_imaginary;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw); // same for |p.s1|

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const auto theta_p = i_tab.nr(new_tw); // upstairs
	  const auto S = pos_to_neg(rd,ww);
	  const Weight rho_r_shift = E.ctxt.to_simple_shift(theta,theta_p,S);
	  bool flipped = E.ctxt.shift_flip(theta,theta_p,S);
	  assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$

	  // downstairs cross by |ww| only has imaginary and complex steps, so
	  // $\alpha_v.(\gamma-\lambda_\rho)$ is unchanged across |ww|
	  const int f = E.gamma_lambda.dot(alpha_v);
	  // number of times $\alpha$ is subtracted from $\gamma-\lambda$


	  const RatWeight new_gamma_lambda =
	    E.gamma_lambda - alpha*f - rho_r_shift;
	  // both $\gamma-\lambda$ and $\tau$ get $f*alpha$ subtracted by
	  // $\alpha$-reflection; adapt $\tau$ for vanishing $1-\delta$ image
	  const Weight new_tau = rd.reflection(n_alpha,E.tau) + alpha*f;

	  // but |dual_v| needs correction by |ell_shift|
	  const int dual_f = (E.ctxt.g_rho_check() - E.l).dot(alpha);

	  const Coweight new_l = E.l + alpha_v*dual_f;
          const Coweight new_t =
	    rd.coreflection(E.t,n_alpha) - alpha_v*dual_f;
	  paramin F (E.ctxt, new_tw, new_gamma_lambda, new_tau, new_l, new_t,
		     E.is_flipped()!=flipped);

	  // do extra conditional flip for 2Ci case
	  int ab_tau = (alpha_v+beta_v).dot(E.tau);
	  assert (ab_tau%2==0);
	  F.flip((ab_tau*dual_f)%4!=0);
	  links.push_back(std::move(F));  // "Cayley" link
	} // end of 2Ci case
	else // twisted commutation with |s0.s1|, and not |ascent|: 2Cr
	{ result = two_semi_real;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw); // same for |p.s1|

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const auto theta_p=i_tab.nr(new_tw);
	  const auto S = pos_to_neg(rd,ww);
	  const Weight rho_r_shift = E.ctxt.to_simple_shift(theta,theta_p,S);
	  bool flipped = E.ctxt.shift_flip(theta,theta_p,S);
	  assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$

	  const int f = level_a(E,rho_r_shift,n_alpha);

	  const RatWeight new_gamma_lambda = // \emph{reflect} parallel to alpha
	    E.gamma_lambda + rho_r_shift - alpha*f;
	  const Weight new_tau = rd.reflection(n_alpha,E.tau) - alpha*f;

	  const int dual_f = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	  const Coweight new_l = E.l + alpha_v*dual_f;
          const Coweight new_t =
	    rd.coreflection(E.t,n_alpha) + alpha_v*dual_f;

	  paramin F (E.ctxt, new_tw, new_gamma_lambda, new_tau, new_l, new_t,
		   E.is_flipped()!=flipped);

	  // do extra conditional flip for 2Cr case
	  int t_ab = E.t.dot(beta-alpha);
	  assert(t_ab%2==0);
	  F.flip((t_ab * (f+alpha_v.dot(E.tau)))%4!=0);
	  links.push_back(std::move(F));  // "Cayley" link
	} // end of 2Cr case
      }
    }
    break;
  case ext_gen::three:
    { const Weight& alpha = integr_datum.simpleRoot(p.s0);
      const Coweight& alpha_v = integr_datum.simpleCoroot(p.s0);
      RootNbr n_alpha = subs.parent_nr_simple(p.s0);
      RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);
      const Weight& beta = integr_datum.simpleRoot(p.s1);
      RootNbr n_beta = subs.parent_nr_simple(p.s1);
      const Coweight& beta_v = integr_datum.simpleCoroot(p.s1);

      RootNbr n_kappa =integr_datum.simple_reflected_root
	 (p.s1, integr_datum.simpleRootNbr(p.s0));
      WeylWord s_kappa = subs.reflection(integr_datum.posRootIndex(n_kappa));

      const Weight& kappa = integr_datum.root(n_kappa);
      assert (kappa==alpha+beta);
      const Coweight& kappa_v = integr_datum.coroot(n_kappa);

      const Weight beta_alpha = beta - alpha;

      const TwistedInvolution new_tw = tW.prod(s_kappa,E.tw); // when applicable

      if (theta_alpha==n_alpha) // length 3 imaginary case
      { // first find out if the simply-integral root $\alpha$ is compact
	int tf_alpha = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	int tf_beta = (E.ctxt.g_rho_check() - E.l).dot(beta);
	assert((tf_alpha-tf_beta)%2==0); // same compactness
	if (tf_alpha%2!=0) // then $\alpha$ and $\beta$ are compact
	  return three_imaginary_compact;

	// noncompact case
	result = three_imaginary_semi; // this is the 3i case

	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	const auto theta_p = i_tab.nr(new_tw); // upstairs

	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = E.ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = E.ctxt.shift_flip(theta,theta_p,S);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$
	assert(rd.is_simple_root(alpha_simple)); // cannot fail for length 3

	E0.tau -= alpha*kappa_v.dot(E.tau); // make |kappa_v.dot(E.tau)==0|
	E0.l += alpha_v*(tf_alpha+tf_beta);
	E0.t += (beta_v-alpha_v)*((tf_alpha+tf_beta)/2);

	paramin F(E.ctxt, new_tw,
		  E0.gamma_lambda-rho_r_shift, E0.tau, E0.l, E0.t);

	flipped = not flipped; // January unsurprise for 3i: delta acts by -1
	z_align(E0,F,flipped^not same_sign(E,E0));

	links.push_back(std::move(F)); // Cayley link
      } // end of 3i case
      else if (theta_alpha==rd.rootMinus(n_alpha)) // length 3 real case
      {
	RootNbr alpha_simple = n_alpha;
	const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	assert(rd.is_simple_root(alpha_simple)); // no complications here

	const auto theta_p=i_tab.nr(new_tw);
	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = E.ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = E.ctxt.shift_flip(theta,theta_p,S);
	assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // as $ww\in W^\delta$

	const int a_level = level_a(E,rho_r_shift,n_alpha);

	if (a_level%2!=0) // nonparity
	   return three_real_nonparity; // no link added here

	// parity case
	result = three_real_semi; // this is the 3r case

	const int b_level = level_a(E,rho_r_shift,n_beta);
	assert(b_level%2==0); // since |a_level| and |b_level| have same parity

	const RatWeight new_gamma_lambda = // make level for |kappa| 0
	  E.gamma_lambda+rho_r_shift - alpha*(a_level+b_level); // even multiple

	E0.t -= alpha_v*kappa.dot(E.t); // makes |E.t.dot(kappa)==0|
	E0.gamma_lambda -= alpha*(a_level+b_level); // even multiple of |alpha|
	E0.tau += beta_alpha*((a_level+b_level)/2);
	assert(same_sign(E,E0)); // neither |t| change nor 2*real_root matter
	assert(E0.gamma_lambda+rho_r_shift==new_gamma_lambda);
	validate(E0);

	flipped = not flipped; // January unsurprise for 3r

	paramin F(E.ctxt, new_tw, new_gamma_lambda,E0.tau,E0.l,E0.t);

	z_align(E0,F,flipped); // no 4th arg since |E.t.dot(kappa)==0|
	links.push_back(std::move(F)); // Cayley link
      } // end of 3r case
      else // length 3 complex case (one of 3Ci or 3Cr or 3C+/-)
      { const bool ascent = rd.is_posroot(theta_alpha);
	if (theta_alpha == (ascent ? n_beta : rd.rootMinus(n_beta)))
	{ // reflection by |alpha+beta| twisted commutes with |E.tw|: 3Ci or 3Cr
	  result = ascent ? three_semi_imaginary : three_semi_real;

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(E.ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const auto theta_p=i_tab.nr(new_tw);
	  const auto S = pos_to_neg(rd,ww);
	  const Weight rho_r_shift = ascent
	    ?  E.ctxt.to_simple_shift(theta,theta_p,S)
	    : -E.ctxt.to_simple_shift(theta,theta_p,S);
	  bool flipped = E.ctxt.shift_flip(theta,theta_p,S);
	  assert(E.ctxt.delta()*rho_r_shift==rho_r_shift); // $ww\in W^\delta$

	  int tf_alpha = (E.ctxt.g_rho_check() - E.l).dot(alpha);
	  int dtf_alpha = E.gamma_lambda.dot(alpha_v);
	  RatWeight new_gamma_lambda = E.gamma_lambda - rho_r_shift; // for now

	  if (ascent) // 3Ci
	  {
	    if (dtf_alpha%2!=0)
	    { new_gamma_lambda -= beta_alpha;
	      E0.gamma_lambda -= beta_alpha;
	      E0.tau -= beta_alpha;
	    }
	    E0.l += kappa_v*tf_alpha;
	    E0.tau -= kappa*(kappa_v.dot(E.tau)/2);
	    validate(E0);
	    assert(E0.t.dot(kappa)==0);

	    paramin F(E.ctxt,new_tw, new_gamma_lambda, E0.tau, E0.l, E0.t);

	    flipped = not flipped; // January unsurprise for 3Ci
	    z_align(E0,F, flipped^(not same_sign(E,E0)));
	    links.push_back(std::move(F)); // Cayley link
	  }
	  else // descent, so 3Cr
	  {
	    E0.gamma_lambda  -= kappa*dtf_alpha;
	    new_gamma_lambda -= kappa*dtf_alpha;

	    E0.t -= kappa_v*(kappa.dot(E.t)/2); // makes |E0.t.dot(kappa)==0|
	    if (tf_alpha%2!=0)
	    {
	      auto b_a = beta_v-alpha_v; // or |kappa_v-alpha_v*2|
	      E0.l += b_a;
	      E0.t -= b_a;
	    }
	    paramin F(E.ctxt, new_tw, new_gamma_lambda, E0.tau, E0.l, E0.t);

	    flipped = not flipped; // January unsurprise for 3Cr
	    z_align(E0,F,flipped^not same_sign(E,E0));
	    // there was no 4th argument there since |E.t.dot(kappa)==0|
	    links.push_back(std::move(F)); // Cayley link
	  }

	} // end of 3ci and 3Cr cases
	else // twisted non-commutation: 3C+ or 3C-
	{
	  result = ascent ? three_complex_ascent : three_complex_descent;
	  links.push_back(complex_cross(p,E));
	}
      }
    }
    break;
  }

  // October surprise: add a flip to links with a length difference of 2
  if (p.length()-(has_defect(result)?1:0)==2)
  { auto it=links.begin(); auto c=scent_count(result);
    for (unsigned i=0; i<c; ++i,++it) // only affect ascent/descent links
      it->flip(); // do the flip
  }

  return result;
} // |star|


// additional |ext_block| methods

ext_block::ext_block
  (const blocks::block_minimal& block, const WeightInvolution& delta)
  : parent(block)
  , orbits(block.fold_orbits(delta))
  , folded(block.Dynkin().folded(orbits))
  , d_delta(delta)
  , info()
  , data(orbits.size()) // create that many empty vectors
  , l_start(parent.length(parent.size()-1)+2,0)
{
  BitMap fixed_points(block.size());

  // compute the delta-fixed points of the block

  // test if twisting some block element lands in the same block
  if (twisted(block,0,delta)==UndefBlock)
    return; // if block not delta-stable, leave |size==0| and quit

  for (BlockElt z=0; z<block.size(); ++z)
    if (twisted(block,z,delta)==z)
      fixed_points.insert(z);

  complete_construction(fixed_points);
  if (not tune_signs(block))
    throw std::runtime_error("Failure detected in extended block construction");

}

template<typename C> // matrix coefficient type (signed)
containers::simple_list<BlockElt> // returns list of elements selected
  ext_block::condense
    (matrix::Matrix<C>& M, const blocks::block_minimal& parent,
     const RatWeight& gamma) const
{
  const auto& integral_pre_datum = parent.integral_subsystem().pre_root_datum();
  RankFlags sing_orbs;
  for (weyl::Generator s=0; s<rank(); ++s)
    if (gamma.dot(integral_pre_datum.simple_coroot(orbits[s].s0))==0)
      sing_orbs.set(s);
  containers::simple_list<BlockElt> result;

  for (BlockElt y=M.numRows(); y-->0; ) // reverse loop is essential here
  { auto s = first_descent_among(sing_orbs,y);
    if (s==rank())
      result.push_front(y); // no singular descents, so a survivor
    else // a singular descent found, so not a survivor
    { // we contribute row |y| to all its descents by |s| with sign |-1|
      // then conceptually we clear row |y|, but don't bother: it gets ignored
      auto type=descent_type(s,y);
      if (is_like_compact(type))
	continue; // no descents, |y| represents zero; nothing to do for |y|

      // length difference and October surprise combine to always give -1
      const C c(-1); // surprise!
      if (has_double_image(type)) // 1r1f, 2r11
      { auto pair = Cayleys(s,y);
	M.rowOperation(pair.first,y,c*epsilon(s,pair.first,y));
	M.rowOperation(pair.second,y,c*epsilon(s,pair.second,y));
      }
      else
      { auto x = some_scent(s,y);
	M.rowOperation(x,y,c*epsilon(s,x,y));
      }
    }
  }
  return result;
} // |ext_block::condense|

bool ext_block::tune_signs(const blocks::block_minimal& block)
{
  context_minimal ctxt (block.context(),delta(),block.integral_subsystem());
  containers::sl_list<paramin> links;
  for (BlockElt n=0; n<size(); ++n)
  { auto z=this->z(n);
    const paramin E(ctxt,block.x(z),block.gamma_lambda(z));
    for (weyl::Generator s=0; s<rank(); ++s)
    { const ext_gen& p=orbit(s); links.clear(); // output arguments for |star|
      auto tp = star(E,p,links);
      if (tp!=descent_type(s,n))
	return false;

      auto it = links.begin();

      switch (tp)
      {
      case one_imaginary_pair_switched: case one_real_pair_switched:
      case one_real_nonparity: case one_imaginary_compact:
      case two_imaginary_single_double_switched:
      case two_real_single_double_switched:
      case two_real_nonparity: case two_imaginary_compact:
      case three_real_nonparity: case three_imaginary_compact:
	assert(links.empty()); break;
      case one_complex_ascent: case one_complex_descent:
      case two_complex_ascent: case two_complex_descent:
      case three_complex_ascent: case three_complex_descent:
	{ assert(links.size()==1);
	  BlockElt m=cross(s,n); // cross neighbour as bare element of |*this|
	  BlockElt cz = this->z(m); // corresponding element of (parent) |block|
	  paramin F(ctxt,block.x(cz),block.gamma_lambda(cz)); // default extn
	  assert(same_standard_reps(*it,F)); // must lie over same parameter
	  if (not same_sign(*it,F))
	    flip_edge(s,n,m);
	} break;
      case one_imaginary_single: case one_real_single:
      case two_imaginary_single_single: case two_real_single_single:
	{ assert(links.size()==2);
	  BlockElt m=some_scent(s,n); // the unique (inverse) Cayley
	  BlockElt Cz = this->z(m); // corresponding element of block
	  paramin F(ctxt,block.x(Cz),block.gamma_lambda(Cz));
	  assert(same_standard_reps(*it,F));
	  if (not same_sign(*it,F))
	    flip_edge(s,n,m);
	  ++it;
	  m=cross(s,n); BlockElt cz = this->z(m);
	  paramin Fc(ctxt,block.x(cz),block.gamma_lambda(cz));
	  assert(same_standard_reps(*it,Fc));
	  if (not same_sign(*it,Fc))
	    flip_edge(s,n,m);
	} break;
      case two_semi_imaginary: case two_semi_real:
      case three_semi_imaginary: case three_real_semi:
      case three_imaginary_semi: case three_semi_real:
	{ assert(links.size()==1);
	  BlockElt m=some_scent(s,n); // the unique (inverse) Cayley
	  BlockElt Cz = this->z(m); // corresponding element of block
	  paramin F(ctxt,block.x(Cz),block.gamma_lambda(Cz));
	  assert(same_standard_reps(*it,F));
	  if (not same_sign(*it,F))
	    flip_edge(s,n,m);
	} break;
      case one_imaginary_pair_fixed: case one_real_pair_fixed:
      case two_imaginary_double_double: case two_real_double_double:
	{ assert(links.size()==2);
	  BlockEltPair m=Cayleys(s,n);
	  BlockElt Cz0 = this->z(m.first); BlockElt Cz1= this->z(m.second);
	  paramin F0(ctxt,block.x(Cz0),block.gamma_lambda(Cz0));
	  paramin F1(ctxt,block.x(Cz1),block.gamma_lambda(Cz1));
	  bool straight=same_standard_reps(*it,F0);
          const auto& node0 = straight ? *it : *std::next(it);
          const auto& node1 = straight ? *std::next(it) : *it;
	  if (not straight)
	    assert(same_standard_reps(node0,F0));
	  assert(same_standard_reps(node1,F1));
	  if (not same_sign(node0,F0))
	    flip_edge(s,n,m.first);
	  if (not same_sign(node1,F1))
	    flip_edge(s,n,m.second);
	} break;
      case two_imaginary_single_double_fixed: case two_real_single_double_fixed:
	{ assert(links.size()==2);
	  BlockEltPair m=Cayleys(s,n);
	  BlockElt Cz0 = this->z(m.first); BlockElt Cz1= this->z(m.second);
	  paramin F0(ctxt,block.x(Cz0),block.gamma_lambda(Cz0));
	  paramin F1(ctxt,block.x(Cz1),block.gamma_lambda(Cz1));
	  bool straight=same_standard_reps(*it,F0);
          const auto& node0 = straight ? *it : *std::next(it);
          const auto& node1 = straight ? *std::next(it) : *it;
	  if (not straight)
	    assert(same_standard_reps(node0,F0));
	  assert(same_standard_reps(node1,F1));
	  if (not same_sign(node0,F0))
	    flip_edge(s,n,m.first);
	  if (not same_sign(node1,F1))
	    flip_edge(s,n,m.second);
	} break;
      } // |switch(tp)|
    } // |for(s)|
  } // |for(n)|
#ifndef NDEBUG // when debugging test braid relations for each extended block
  BitMap dummy(size());
  for (BlockElt x=0; x<size(); ++x)
    for (weyl::Generator s=0; s<rank(); ++s)
    {
      if (not check_quadratic(*this,s,x))
	throw std::runtime_error
	  ("Quadratic relation failure in extended block construction");
      for (weyl::Generator t=s+1; t<rank(); ++t)
	if (not check_braid(*this,s,t,x,dummy))
	  throw std::runtime_error
	    ("Braid relation failure in extended block construction");
    }
#endif
  return true; // report success if we get here
} // |tune_signs|

template containers::simple_list<BlockElt> ext_block::condense
  (matrix::Matrix<Split_integer>& M, const blocks::block_minimal& parent,
   const RatWeight& gamma) const;

} // |namespace ext_block|

} // |namespace atlas|
