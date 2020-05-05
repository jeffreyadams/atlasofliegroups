/*
  This is common_blocks.cpp

  Copyright (C) 2019,2020 Marc van Leeuwen
  Part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/


/*
  Motivation for this compilation unit.

  The computation of KLV polynomials is tightly bound to the notion of block,
  and indeed in takes place in a |kl::KL_table| value that is constructed from
  a block, and owned by a pointer member |kl_tab_ptr| of |blocks::Block_base|.

  In unitarity computations, blocks are computed from |Param| values, which
  internally correspond to the type |StandardRepr|; each one gives a single
  element of a block. This implies each time an entire block, or at least the
  part below the intitial parameter, is generated. To avoid excessive
  wastefulness, values derived from the KLV polynomials are associated in a
  permanent way with the parameters of the block in a |repr::Rep_table|, so that
  for a parameter whose block has already been subject to KLV computations, a
  second computation is avoided. However many parameters can be seen to have
  isomorphic blocks, namely when they share their |KGB| set and the integral
  subdatum (as determined by the infinitesimal character |gamma|, indeed already
  by its class modulo $X^*$); in this unit, blocks of type |common_block| are
  constructed that only use such information, and still allow computation of KLV
  polynomials (ordinary and twisted).

  It is not clear that it is necessary to actually use blocks that are ignorant
  of the infinitesimal character other than through the integral subsystem it
  determines, as it is probably possible to re-use polynomials associated to
  parameters whose block is isomorphic, by applying a shift to their
  infinitesimal character. But in any case showing that the blocks can be
  generated in this way is a good way to get convinced that an isomorphism as
  indicated does exist.
 */
#include "common_blocks.h"

#include <cassert>
#include <vector>
#include <map>
#include <algorithm>

#include "poset.h"
#include "matreduc.h"
#include "realredgp.h"
#include "subsystem.h"
#include "kgb.h"
#include "repr.h"
#include "ext_block.h"
#include "ext_kl.h" // for |ext_block::ext_kl| contructor defined below

namespace atlas {

namespace repr {

Repr_mod_entry::Repr_mod_entry(const Rep_context& rc, const StandardReprMod& srm)
  : x(srm.x())
  , y(srm.y().data())
  , mask(rc.inner_class().involution_table().y_mask(rc.kgb().inv_nr(x)))
{}

StandardReprMod Repr_mod_entry::srm
  (const Rep_context& rc,const RatWeight& gamma_mod_1) const
{ TorusPart yv(y,rc.inner_class().involution_table().tp_sz(rc.kgb().inv_nr(x)));
  return StandardReprMod::build
    (rc,gamma_mod_1, x, rc.gamma_lambda(rc.kgb().inv_nr(x),yv,gamma_mod_1));
}

} // |namespace repr|

namespace blocks {

// Some declarations of undeclared stuff in blocks.cpp that we (continue to) use

void add_z(block_hash& hash,KGBElt x,KGBElt y);
BlockElt find_in(const block_hash& hash,KGBElt x,KGBElt y);
BlockElt& first_free_slot(BlockEltPair& p);

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
  , srm_hash(z_pool)
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
  containers::queue<
    containers::sl_list<containers::sl_list<repr::StandardReprMod> >
    >elements; // involution packets, |x| varies in outer, |y| in inner list
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
      srm_hash.match(repr::Repr_mod_entry(rc,zz));
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
	/*
	  queue.push(repr::StandardReprMod::build
		     (rc,srm.gamma_mod1(),
		      highest_x, y_pool.back().repr().log_pi(false)));
	*/
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

  do // process involution packet of elements from |next| to |queue.front()|
  { // |next| is constant throughout the loop body, popped from |queue| at end
    const KGBElt first_x = x(next), first_y=y(next);

    // precompute (reversed) length for anything generated this iteration
    const auto next_length = info[next].length+1;

    const auto bundle = std::move(elements.front());
    elements.pop();
    const unsigned int nr_x = bundle.size();
    const unsigned int nr_y = bundle.front().size();
#ifndef NDEBUG
    assert((queue.front()-next)==nr_x*nr_y);
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
    ndebug_use(first_y);

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

      { // compute values into |cross_ys|, and check for existence of real parity

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

      { // handle cross actions and descent statusses in all cases
	const auto theta = kgb.inv_nr(ctxt.cross(s,head).x());
	BlockElt cur = next; // start of old involution packet

	if (cross_new_involution)
	{ // add a new involution packet
	  containers::sl_list<containers::sl_list<repr::StandardReprMod> >
	    packet;
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
	      srm_hash.match(repr::Repr_mod_entry(rc,new_srm));
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
	  }
	  assert(Cayley_ys.at_end(it));
	}
      } // |if (not Cayley_ys.empty())|
    } // |for(s)|
  }
  while (next=queue.front(), queue.pop(), not queue.empty());
  // end of step 4

  highest_x = last(x_seen); // to be sure; length need not increase with |x|
  highest_y=y_hash.size()-1; // set highest occurring |y| value, for |ysize|

  // reverse lengths; reorder by increasing length then value of |x|
  sort(info.back().length,true);

  entry_element = lookup(srm); // look up element matching the original input

} // |common_block::common_block|, full

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

  Block_base::dd = // integral Dynkin diagram, converted from dual side
    DynkinDiagram(integral_sys.cartanMatrix().transposed());

  y_entry::Pooltype y_pool;
  y_part_hash y_hash(y_pool);

  std::vector<containers::sl_list<TorusPart> > y_table
    (inner_class().involution_table().size());
  // every element of |y_table| is a list of |TorusPart| values of the same rank
  // therefore, we can (and will) compare list elements and sort the lists
  for (unsigned long elt : elements)
  { const auto& srm = rt.srm(elt);
    const KGBElt x = srm.x();
    const TorusPart& y = srm.y();
    if (x>highest_x)
      highest_x=x;
    auto& loc = y_table[kgb.inv_nr(x)];
    auto it = std::find_if_not(loc.cbegin(),loc.cend(),
			       [&y](const TorusPart& t) { return t<y; });
    if (it==loc.end() or y<*it) // skip if |y| already present in the list
      loc.insert(it,y);
  }

  const auto& i_tab = inner_class().involution_table();
  for (InvolutionNbr i_x=y_table.size(); i_x-->0; )
  {
    highest_y += y_table[i_x].size();
    for (TorusPart& y : y_table[i_x])
    {
      RatWeight gamma_lambda = rt.gamma_lambda(i_x,y,gamma_mod_1);
      TorusElement t = y_values::exp_pi(gamma_lambda);
      y_hash.match(i_tab.pack(t,i_x)); // enter this |y_entry| into |y_hash|
    }
  }
  assert(y_pool.size()==highest_y);
  -- highest_y; // one less than the number of distinct |y| values

  elements.sort // pre-sort by |x| to ensure descents precede in setting lengths
    ([&rt](unsigned long a, unsigned long b)
          { return rt.srm(a).x()<rt.srm(b).x(); }
     );

  for (unsigned long elt : elements)
  { const auto& srm=rt.srm(elt);
    const KGBElt x=srm.x();
    auto y = y_hash.find(i_tab.pack(rt.y_as_torus_elt(srm),kgb.inv_nr(x)));
    assert(y!=y_hash.empty);
    info.emplace_back(x,y); // leave descent status unset and |length==0| for now
    srm_hash.match(repr::Repr_mod_entry(rc,srm));
  }

  // allocate link fields with |UndefBlock| entries
  data.assign(integral_sys.rank(),std::vector<block_fields>(elements.size()));

  assert(info.size()==elements.size());
  unsigned short max_length=0;
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
	    {
	      info[i].length=length(sz)+1;
	      if (info[i].length>max_length)
		max_length=info[i].length;
	    }
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
	    {
	      info[i].length=length(sz)+1;
	      if (info[i].length>max_length)
		max_length=info[i].length;
	    }
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

  sort(max_length,false);  // finally sort by length, then |x|

} // |common_block::common_block|, partial

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

ext_gens common_block::fold_orbits(const WeightInvolution& delta) const
{
  return rootdata::fold_orbits(integral_sys.pre_root_datum(),delta);
}

ext_block::ext_block& common_block::extended_block
  (const WeightInvolution& delta)
{
  if (extended.get()==nullptr)
    extended.reset(new ext_block::ext_block(*this,delta));
  return *extended;
}

// integrate an older partial block, with mapping of elements
void common_block::swallow(common_block&& sub, const BlockEltList& embed)
{
  for (BlockElt z=0; z<sub.size(); ++z)
  {
    auto& covered = std::move(sub).Bruhat_order().Hasse(z);
    for (auto& c : covered)
      c=embed[c]; // translate in place
    set_Bruhat_covered(embed[z],std::move(covered));
  }
}

void common_block::set_Bruhat
  (containers::sl_list<std::pair<BlockElt,BlockEltList> >&& partial_Hasse)
{
  for (auto& pair : partial_Hasse)
    set_Bruhat_covered(pair.first,std::move(pair.second));
}

void common_block::sort(unsigned short max_length, bool reverse_length)
{
  const KGBElt x_lim=highest_x+1; // limit for |x| values
  std::vector<unsigned> value(size(),0u); // values to be ranked below

  for (BlockElt i=0; i<size(); ++i)
  { assert(length(i)<=max_length);
    if (reverse_length)
       info[i].length = max_length-length(i); // reverse length
    value[i]= info[i].length*x_lim+x(i); // length has priority over value of |x|
  }

  Permutation ranks = // standardization permutation, to be used for reordering
    permutations::standardization(value,(max_length+1)*x_lim,nullptr);

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

} // |namespace blocks|

namespace repr {

// |Ext_rep_context| methods

Ext_rep_context::Ext_rep_context
  (const repr::Rep_context& rc, const WeightInvolution& delta)
: Rep_context(rc), d_delta(delta) {}

Ext_rep_context::Ext_rep_context (const repr::Rep_context& rc)
: Rep_context(rc), d_delta(rc.inner_class().distinguished()) {}

// |common_context| methods

common_context::common_context (RealReductiveGroup& G, const SubSystem& sub)
: Rep_context(G)
, integr_datum(sub.pre_root_datum())
, sub(sub)
{} // |common_context::common_context|

struct ulong_entry
{ unsigned long val;
  ulong_entry(unsigned long n) : val(n) {}
  typedef std::vector<ulong_entry> Pooltype;
  bool operator != (ulong_entry x) const { return val!=x.val;  }
  size_t hashCode(size_t modulus) const { return (5*val)&(modulus-1); }
};

// |Rep_table| helper class
class Rep_table::Bruhat_generator
{
  Rep_table& parent;
  const common_context& ctxt;
  std::vector<ulong_entry> pool;
  HashTable<ulong_entry,BlockElt> local_h; // hash table, avoid name |hash|
  std::vector<containers::simple_list<unsigned long> > predecessors;
public:
  Bruhat_generator (Rep_table* caller, const common_context& ctxt)
    : parent(*caller),ctxt(ctxt), pool(), local_h(pool), predecessors() {}

  bool in_interval (unsigned long n) const
  { return local_h.find(n)!=local_h.empty; }
  const containers::simple_list<unsigned long>& covered(unsigned long n) const
  { return predecessors.at(local_h.find(n)); }
  containers::simple_list<unsigned long> block_below(const StandardReprMod& srm);
}; // |class Rep_table::Bruhat_generator|


// |Rep_table| methods
blocks::common_block& Rep_table::add_block_below
  (const common_context& ctxt, const StandardReprMod& srm, BitMap* subset)
{
  assert(mod_hash.find(srm)==mod_hash.empty); // otherwise don't call us
  Bruhat_generator gen(this,ctxt);
  const auto prev_size = mod_pool.size(); // limit of previously known elements
  containers::sl_list<unsigned long> elements(gen.block_below(srm));

  using sub_pair =
    std::pair<blocks::common_block*,containers::sl_list<BlockElt> >;
  containers::sl_list<sub_pair> sub_blocks;
  for (auto z : elements)
    if (z<prev_size) // then |z| was already known as a partial block element
    { // record block pointer and index of |z| in block into |sub_blocks|
      const auto block_p = &*place[z].first;
      const BlockElt z_rel = place[z].second;
      auto it = std::find_if
	(sub_blocks.begin(),sub_blocks.end(),
	 [block_p](const sub_pair& pair) { return &*pair.first==block_p; }
	 );
      if (sub_blocks.at_end(it)) // then we have a fresh sub-block
	sub_blocks.emplace_back(block_p,containers::sl_list<BlockElt>{z_rel});
      else // a sub-block already recorded (at |it|)
	it->second.push_back(z_rel); // append |z_rel| to its list of elements
    }

  for (auto& pair : sub_blocks)
  {
    if (pair.first->size() > pair.second.size()) // then incomplete inclusion
    {
      const auto& block = *pair.first;
      pair.second.sort(); // following loop requires increase
      auto it = pair.second.begin();
      for (BlockElt z=0; z<block.size(); ++z)
	if (not it.at_end() and *it==z)
	  ++it; // skip element already in Bruhat interval
	else // join element outside Bruhat interval to new block
	  elements.push_back(mod_hash.find(block.representative(z)));
    }
    // since all |block| elements are incorporated
    pair.second.clear(); // forget which were in the Bruhat interval
  }

  containers::sl_list<blocks::common_block> temp; // must use temporary singleton
  auto& block = temp.emplace_back // construct block and get a reference
    (*this,ctxt,elements,srm.gamma_mod1());

  *subset=BitMap(block.size());
  containers::sl_list<std::pair<BlockElt,BlockEltList> > partial_Hasse_diagram;
  for (auto z : elements)
    if (gen.in_interval(z)) // these have their covered's in |gen|
    {
      BlockElt i_z = block.lookup(this->srm(z)); // index of |z| in our new block
      subset->insert(i_z); // mark |z| as element ot the Bruhat interval
      const auto& covered = gen.covered(z);
      BlockEltList row;
      row.reserve(atlas::containers::length(covered));
      for (auto it=covered.begin(); not covered.at_end(it); ++it)
      {
	const BlockElt y = block.lookup(this->srm(*it)); // get relative number
	assert(y!=UndefBlock);
	row.push_back(y); // store covering relation in |Hasse_diagram|
      }
      partial_Hasse_diagram.emplace_back(i_z,row);
    }

  block.set_Bruhat(std::move(partial_Hasse_diagram));

  static const std::pair<bl_it, BlockElt> empty(bl_it(),UndefBlock);
  place.resize(mod_pool.size(),empty); // new entries point to |new_block_it|

  for (const auto& pair : sub_blocks) // swallow sub-blocks
  {
    auto& sub_block = *pair.first;
    bl_it block_it; // set belowe
    BlockEltList embed; embed.reserve(sub_block.size());
    for (BlockElt z=0; z<sub_block.size(); ++z)
    {
      auto zm = sub_block.representative(z);
      const BlockElt z_rel = block.lookup(zm);
      assert(z_rel!=UndefBlock);
      embed.push_back(z_rel);
    }
    auto h = mod_hash.find(sub_block.representative(0));
    assert(h!=mod_hash.empty);
    block_it = place[h].first;

    assert(&*block_it==&sub_block); // ensure we erase |sub_block|
    block.swallow(std::move(sub_block),embed);
    block_erase(block_it); // even pilfered, the pointer is still unchanged
  }

  // only after the |block_erase| upheavals is it safe to link in the new block
  // also, it can only go to the end of the list, to not invalidate any iterators
  const auto new_block_it=block_list.end(); // iterator for new block elements
  block_list.splice(new_block_it,temp,temp.begin()); // link in |block| at end
  // now make sure for all |elements| that |place| fields are set for new block
  for (const auto& z : elements)
  {
    const BlockElt z_rel = block.lookup(this->srm(z));
    place[z] = std::make_pair(new_block_it,z_rel); // extend or replace
  }
  return block;
} // |Rep_table::add_block_below|


containers::simple_list<unsigned long> Rep_table::Bruhat_generator::block_below
  (const StandardReprMod& srm)
{
  auto& hash=parent.mod_hash;
  { const auto h=hash.find(srm);
    if (h!=hash.empty) // then |srm| was seen earlier
    {
      unsigned hh = local_h.find(h);
      if (hh!=local_h.empty) // then we visited this element in current recursion
	return containers::simple_list<unsigned long>(); // nothing new
    }
  }
  const auto rank = ctxt.id().semisimpleRank();
  containers::sl_list<unsigned long> pred; // list of elements covered by z
  // invariant: |block_below| has been called for every element in |pred|

  weyl::Generator s; // a complex or real type 1 descent to be found, if exists
  containers::sl_list<containers::simple_list<unsigned long> > results;
  for (s=0; s<rank; ++s)
  {
    std::pair<gradings::Status::Value,bool> stat=ctxt.status(s,srm.x());
    if (not stat.second)
      continue; // ignore imaginary, complex ascent or real (potentially) type 2
    if (stat.first==gradings::Status::Complex)
    { // complex descent
      const StandardReprMod sz = ctxt.cross(s,srm);
      results.push_back(block_below(sz)); // recursion
      pred.push_back(hash.find(sz)); // must call |hash.find| after |block_below|
      break; // we shall add $s$-ascents of predecessors of |sz| below
    }
    else if (stat.first==gradings::Status::Real and ctxt.is_parity(s,srm))
    { // |z| has a type 1 real descent at |s|
      const StandardReprMod sz0 = ctxt.down_Cayley(s,srm);
      const StandardReprMod sz1 = ctxt.cross(s,sz0);
      results.push_back(block_below(sz0)); // recursion
      results.push_back(block_below(sz1)); // recursion
      pred.push_back(hash.find(sz0)); // call |hash.find| after |block_below|
      pred.push_back(hash.find(sz1));
      break; // we shall add $s$-ascents of predecessors of |sz_inx| below
    } // |if (real type 1)|

  } // |for (s)|

  // if above loop performed a |break| it found complex or real type 1 descent
  if (s==rank) // otherwise, the only descents are real type 2, if any
  {
    while (s-->0) // we reverse the loop just because |s==rank| already
      if (ctxt.status(s,srm.x()).first==gradings::Status::Real and
	  ctxt.is_parity(s,srm))
      {
	const auto sz = ctxt.down_Cayley(s,srm);
	results.push_back(block_below(sz)); // recursion
	pred.push_back(hash.find(sz));
      }
  }
  else // a complex or real type 1 descent |sz==pred.front()| for |s| was found
  { // add |s|-ascents for elements covered by |sz|
    const auto pred_sz = predecessors.at(local_h.find(pred.front()));
    for (auto it = pred_sz.begin(); not pred.at_end(it); ++it)
    {
      const auto p = *it; // sequence number of a predecessor of |sz|
      const StandardReprMod zp = parent.mod_pool[p];
      std::pair<gradings::Status::Value,bool> stat=ctxt.status(s,zp.x());
      switch (stat.first)
      {
      case gradings::Status::Real: case gradings::Status::ImaginaryCompact:
	break; // nothing to do without ascent
      case gradings::Status::Complex:
	if (not stat.second) // complex ascent
	{
	  const auto szp = ctxt.cross(s,zp);
	  results.push_back(block_below(szp)); // recursion
	  pred.push_back(hash.find(szp));
	} // |if(complex ascent)
	break;
      case gradings::Status::ImaginaryNoncompact:
	{
	  const auto szp = ctxt.up_Cayley(s,zp);
	  results.push_back(block_below(szp)); // recursion
	  pred.push_back(hash.find(szp));
	  if (not stat.second) // then nci type 2
	  {
	    const auto szp1 = ctxt.cross(s,szp);
	    results.push_back(block_below(szp1)); // recursion
	    pred.push_back(hash.find(szp1));
	  }
	}
	break;
      } // |switch(status(s,conj_x))|
    } // |for (it)|
  } // |if (s<rank)|

  const auto h=hash.match(srm); // finally generate sequence number for |srm|
  { // merge all |results| together and remove duplicates
    results.push_front(containers::simple_list<unsigned long> {h} );
    while (results.size()>1)
    {
      auto it=results.begin();
      auto first = std::move(*it);
      first.merge(std::move(*++it));
      results.erase(results.begin(),++it);
      first.unique();
      results.push_back(std::move(first));
    }
  }
  unsigned hh = local_h.match(h); // local sequence number for |srm|
  assert(hh==predecessors.size()); ndebug_use(hh);
  predecessors.push_back(pred.undress()); // store |pred| at |hh|
  return results.front();
} // |Rep_table::Bruhat_generator::block_below|

std::pair<gradings::Status::Value,bool>
  common_context::status(weyl::Generator s, KGBElt x) const
{
  const auto& conj = sub.to_simple(s); // word in full system
  KGBElt conj_x = kgb().cross(conj,x);
  const auto t=sub.simple(s);
  const auto stat = kgb().status(t,conj_x);
  return std::make_pair(stat,
			stat==gradings::Status::Real
			? kgb().isDoubleCayleyImage(t,conj_x) // real type 1
			: stat==gradings::Status::Complex
			? kgb().isDescent(t,conj_x)
			: conj_x!=kgb().cross(t,conj_x)); // nc imaginary type 1
}

StandardReprMod common_context::cross
    (weyl::Generator s, const StandardReprMod& z) const
{
  const auto& refl = sub.reflection(s); // reflection word in full system
  const KGBElt new_x = kgb().cross(refl,z.x());
  RatWeight gamma_lambda = this->gamma_lambda(z);

  const auto& i_tab = inner_class().involution_table();
  RootNbrSet pos_neg = pos_to_neg(root_datum(),refl);
  pos_neg &= i_tab.real_roots(kgb().inv_nr(z.x())); // only real roots for |z|
  gamma_lambda -= root_sum(root_datum(),pos_neg); // correction for $\rho_r$'s
  integr_datum.simple_reflect(s,gamma_lambda.numerator()); // then reflect
  return repr::StandardReprMod::build(*this,z.gamma_mod1(),new_x,gamma_lambda);
}

StandardReprMod common_context::down_Cayley
    (weyl::Generator s, const StandardReprMod& z) const
{
  assert(is_parity(s,z)); // which also asserts that |z| is real for |s|
  const auto& conj = sub.to_simple(s); // word in full system
  KGBElt conj_x = kgb().cross(conj,z.x());
  conj_x = kgb().inverseCayley(sub.simple(s),conj_x).first;
  const auto new_x = kgb().cross(conj_x,conj);
  RatWeight gamma_lambda = this->gamma_lambda(z);

  const auto& i_tab = inner_class().involution_table();
  RootNbrSet pos_neg = pos_to_neg(root_datum(),conj);
  RootNbrSet real_flip = i_tab.real_roots(kgb().inv_nr(z.x()));
  real_flip ^= i_tab.real_roots(kgb().inv_nr(new_x));
  pos_neg &= real_flip; // posroots that change real status and map to negative
  gamma_lambda += root_sum(root_datum(),pos_neg); // correction of $\rho_r$'s
  return repr::StandardReprMod::build(*this,z.gamma_mod1(),new_x,gamma_lambda);
}

bool common_context::is_parity
    (weyl::Generator s, const StandardReprMod& z) const
{
  const auto& i_tab = inner_class().involution_table();
  const auto& real_roots = i_tab.real_roots(kgb().inv_nr(z.x()));
  assert(real_roots.isMember(sub.parent_nr_simple(s)));
  const Coweight& alpha_hat = integr_datum.simpleCoroot(s);
  const int eval = this->gamma_lambda(z).dot(alpha_hat);
  const int rho_r_corr = alpha_hat.dot(root_datum().twoRho(real_roots))/2;
  return (eval+rho_r_corr)%2!=0;
}

StandardReprMod common_context::up_Cayley
    (weyl::Generator s, const StandardReprMod& z) const
{
  const auto& conj = sub.to_simple(s); // word in full system
  KGBElt conj_x = kgb().cross(conj,z.x());
  conj_x = kgb().cayley(sub.simple(s),conj_x);
  const auto new_x = kgb().cross(conj_x,conj);
  RatWeight gamma_lambda = this->gamma_lambda(z);

  const auto& i_tab = inner_class().involution_table();
  const RootNbrSet& upstairs_real_roots = i_tab.real_roots(kgb().inv_nr(new_x));
  RootNbrSet real_flip = upstairs_real_roots;
  real_flip ^= i_tab.real_roots(kgb().inv_nr(z.x())); // remove downstairs reals

  RootNbrSet pos_neg = pos_to_neg(root_datum(),conj);
  pos_neg &= real_flip; // posroots that change real status and map to negative
  gamma_lambda += root_sum(root_datum(),pos_neg); // correction of $\rho_r$'s

  // correct in case the parity condition fails for our raised |gamma_lambda|
  const Coweight& alpha_hat = integr_datum.simpleCoroot(s);
  const int rho_r_corr = // integer since alpha is among |upstairs_real_roots|
    alpha_hat.dot(root_datum().twoRho(upstairs_real_roots))/2;
  const int eval = gamma_lambda.dot(alpha_hat);
  if ((eval+rho_r_corr)%2==0) // parity condition says it should be 1
    gamma_lambda += RatWeight(integr_datum.root(s),2); // add half-alpha

  return repr::StandardReprMod::build(*this,z.gamma_mod1(),new_x,gamma_lambda);
}


Weight common_context::to_simple_shift
  (InvolutionNbr theta, InvolutionNbr theta_p, RootNbrSet S) const
{ const InvolutionTable& i_tab = inner_class().involution_table();
  S &= (i_tab.real_roots(theta) ^i_tab.real_roots(theta_p));
  return root_sum(root_datum(),S);
}


Ext_common_context::Ext_common_context
  (RealReductiveGroup& G, const WeightInvolution& delta, const SubSystem& sub)
    : repr::common_context(G,sub)
    , d_delta(delta)
    , pi_delta(G.root_datum().rootPermutation(delta))
    , delta_fixed_roots(fixed_points(pi_delta))
    , twist()
    , l_shifts(id().semisimpleRank())
{
  const RootDatum& rd = root_datum();
  for (weyl::Generator s=0; s<rd.semisimpleRank(); ++s)
    twist[s] = rd.simpleRootIndex(delta_of(rd.simpleRootNbr(s)));

  // the reflections for |E.l| pivot around |g_rho_check()|
  const RatCoweight& g_rho_check = this->g_rho_check();
  for (unsigned i=0; i<l_shifts.size(); ++i)
    l_shifts[i] = -g_rho_check.dot(id().simpleRoot(i));
} // |Ext_common_context::Ext_common_context|


bool Ext_common_context::is_very_complex
  (InvolutionNbr theta, RootNbr alpha) const
{ const auto& i_tab = inner_class().involution_table();
  const auto& rd = root_datum();
  assert (rd.is_posroot(alpha)); // this is a precondition
  auto image = i_tab.root_involution(theta,alpha);
  make_positive(rd,image);
  return image!=alpha and image!=delta_of(alpha);
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
bool Ext_common_context::shift_flip
  (InvolutionNbr theta, InvolutionNbr theta_p, RootNbrSet S) const
{ S.andnot(delta_fixed()); // $\delta$-fixed roots won't contribute

  unsigned count=0; // will count 2-element |delta|-orbit elements
  for (auto it=S.begin(); it(); ++it)
    if (is_very_complex(theta,*it) != is_very_complex(theta_p,*it) and
	not root_datum().sumIsRoot(*it,delta_of(*it)))
      ++count;

  assert(count%2==0); // since |pos_to_neg| is supposed to be $\delta$-stable
  return count%4!=0;
}

} // |namespace repr|

namespace ext_block
{
  // Declarations of some functions re-used from ext_block.cpp

bool in_L_image(Weight beta,WeightInvolution&& A);
bool in_R_image(WeightInvolution&& A,Coweight b);
unsigned int scent_count(DescValue v);
Coweight ell (const KGB& kgb, KGBElt x);

  // Declarations of some local functions
  WeylWord fixed_conjugate_simple
    (const repr::Ext_common_context& c, RootNbr& alpha);
bool same_standard_reps (const paramin& E, const paramin& F);
bool same_sign (const paramin& E, const paramin& F);
inline bool is_default (const paramin& E)
  { return same_sign(E,paramin(E.ctxt,E.x(),E.gamma_lambda)); }

void z_align (const paramin& E, paramin& F, bool extra_flip);
void z_align (const paramin& E, paramin& F, bool extra_flip, int t_mu);
paramin complex_cross(const repr::Ext_common_context& ctxt,
		      const ext_gen& p, paramin E);
int level_a (const paramin& E, const Weight& shift, RootNbr alpha);
DescValue star (const repr::Ext_common_context& ctxt,
		const paramin& E, const ext_gen& p,
		containers::sl_list<paramin>& links);

bool is_descent
  (const repr::Ext_common_context& ctxt, const ext_gen& kappa, const paramin& E);
weyl::Generator first_descent_among
  (const repr::Ext_common_context& ctxt, RankFlags singular_orbits,
   const ext_gens& orbits, const paramin& E);


/* Try to conjugate |alpha| by product of folded-generators for the (full)
   root system of |c| to a simple root, and return the left-conjugating word
   that was applied. This may fail, if after some conjugation one ends up with
   the long root of a nontrivially folded A2 subsystem (in which case there
   cannot be any solution because |alpha| is fixed by the involution but none
   of the simple roots in its component of the root system is). In this case
   |alpha| is left as that non simple root, and the result conjugates to it.
 */
  WeylWord fixed_conjugate_simple
    (const repr::Ext_common_context& ctxt, RootNbr& alpha)
{ const RootDatum& rd = ctxt.root_datum();

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

const WeightInvolution& paramin::theta () const
 { return ctxt.inner_class().matrix(tw); }

void validate(const paramin& E)
{
#ifndef NDEBUG // make sure this is a no-op when debugging is disabled
  const auto& i_tab = E.ctxt.inner_class().involution_table();
  const auto& theta = i_tab.matrix(E.tw);
  const auto& delta = E.ctxt.delta();
  assert(delta*theta==theta*delta);
  const Weight diff = E.gamma_lambda.integer_diff<int>(delta*E.gamma_lambda);
  assert(diff == (1-theta)*E.tau);
  assert((delta-1).right_prod(E.l)==(theta+1).right_prod(E.t));
  assert(((E.ctxt.g_rho_check()-E.l)*(1-theta)).numerator().isZero());
  assert(((theta+1)*E.gamma_lambda).numerator().isZero());
#endif
}


paramin::paramin
(const repr::Ext_rep_context& ec, const TwistedInvolution& tw,
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


// contructor used for default extension once |x| and |gamma_lamba| are chosen
paramin::paramin
(const repr::Ext_rep_context& ec,
   KGBElt x, const RatWeight& gamma_lambda, bool flipped)
  : ctxt(ec)
  , tw(ec.real_group().kgb().involution(x)) // now computing |theta()| is valid
  , l(ell(ec.real_group().kgb(),x))
  , gamma_lambda(gamma_lambda)
  , tau(matreduc::find_solution
	(1-theta(), gamma_lambda.integer_diff<int>(delta()*gamma_lambda)))
  , t(matreduc::find_solution
	(theta().transposed()+1,(delta()-1).right_prod(l)))
  , flipped(flipped)
{
  validate(*this);
}

// build a default extended parameter for |sr| in the context |ec|
/*
  Importantly, this does not use |sr.gamma()| otherwise than for asserting its
  $\delta$-stability: the same default is used in |common_block| for an entire
  family of blocks, so dependence on |gamma| must be limited to dependence on
  its reduction modulo 1. Even though |gamma_lambda| is computed at the non
  $\delta$-fixed |srm.gamma_mod1()|, it is also (due to the way |mod_reduce|
  works) a proper value of |gamma_lambda| at |sr|, so |(1-delta)*gamma_lambda|,
  gives a valid value for the equation of which |tau| est une solution. However
  |gamma_lambda| may be a different representative than |rc.gamma_lambda(sr)|,
  so don't use that latter: it would give an undesired dependence on |gamma|.
*/
paramin paramin::default_extend
(const repr::Ext_rep_context& ec, const repr::StandardRepr& sr)
{
  assert(((1-ec.delta())*sr.gamma().numerator()).isZero());

  auto srm =  repr::StandardReprMod::mod_reduce(ec,sr);
  // get default representative at |gamma%1|, normalised
  auto gamma_lambda=ec.gamma_lambda(srm);
  return paramin(ec,sr.x(),gamma_lambda);
}

paramin& paramin::operator= (const paramin& p)
{ assert(&ctxt==&p.ctxt); // assignment should remain in the same context
  tw=p.tw;
  l=p.l; gamma_lambda=p.gamma_lambda; tau=p.tau; t=p.t;
  flipped=p.flipped;
  return *this;
}
paramin& paramin::operator= (paramin&& p)
{ assert(&ctxt==&p.ctxt); // assignment should remain in the same context
  tw=std::move(p.tw);
  l=std::move(p.l); gamma_lambda=std::move(p.gamma_lambda);
  tau=std::move(p.tau); t=std::move(p.t);
  flipped=p.flipped;
  return *this;
}


KGBElt paramin::x() const
{ TitsElt a(ctxt.inner_class().titsGroup(),TorusPart(l),tw);
  return rc().kgb().lookup(a);
}

repr::StandardRepr paramin::restrict(const RatWeight& gamma) const
{
  const RatWeight gamma_rho = gamma-rho(rc().root_datum());
  const auto lambda_rho = gamma_rho.integer_diff<int>(gamma_lambda);
  return rc().sr_gamma(x(),lambda_rho,gamma);
}


// whether |E| and |F| lie over equivalent |StandardRepr| values
bool same_standard_reps (const paramin& E, const paramin& F)
{
  if (&E.ctxt!=&F.ctxt)
  { if (&E.ctxt.inner_class()!=&F.ctxt.inner_class())
      throw std::runtime_error
	("Comparing extended parameters from different inner classes");
    if (E.delta()!=F.delta() or E.ctxt.g_rho_check()!=F.ctxt.g_rho_check())
      return false;
    // otherwise (contexts differ, but agree on vitals), match could still occur
  } // so fall through
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
paramin complex_cross(const repr::Ext_common_context& ctxt,
		      const ext_gen& p, paramin E) // by-value for |E|, modified
{ const RootDatum& rd = E.rc().root_datum();
  const RootDatum& id = ctxt.id();
  const InvolutionTable& i_tab = E.rc().inner_class().involution_table();
  auto &tW = E.rc().twisted_Weyl_group(); // caution: |p| refers to integr. datum
  const SubSystem& subs=ctxt.subsys();

  InvolutionNbr theta = i_tab.nr(E.tw);
  const RootNbrSet& theta_real_roots = i_tab.real_roots(theta);
  Weight rho_r_shift = rd.twoRho(theta_real_roots);
  Coweight dual_rho_im_shift = rd.dual_twoRho(i_tab.imaginary_roots(theta));

  for (unsigned i=p.w_kappa.size(); i-->0; ) // at most 3 letters, right-to-left
  { weyl::Generator s=p.w_kappa[i]; // generator for integrality datum
    tW.twistedConjugate(subs.reflection(s),E.tw);
    id.simple_reflect(s,E.gamma_lambda.numerator());
    id.simple_reflect(s,rho_r_shift);
    id.simple_reflect(s,E.tau);
    id.simple_coreflect(E.l,s,ctxt.l_shift(s));
    id.simple_coreflect(dual_rho_im_shift,s);
    id.simple_coreflect(E.t,s);
  }

  InvolutionNbr new_theta = i_tab.nr(E.tw);
  const RootNbrSet& new_theta_real_roots = i_tab.real_roots(new_theta);
  rho_r_shift -= rd.twoRho(new_theta_real_roots);
  rho_r_shift/=2; // now it is just a sum of (real) roots
  E.gamma_lambda += rho_r_shift;

  assert(ctxt.delta()*rho_r_shift==rho_r_shift); // diff of $\delta$-fixed

  dual_rho_im_shift -= rd.dual_twoRho(i_tab.imaginary_roots(new_theta));
  dual_rho_im_shift/=2; // now it is just a sum of (imaginary) coroots
  E.l -= dual_rho_im_shift;

  assert(ctxt.delta().right_prod(dual_rho_im_shift)==dual_rho_im_shift);
  validate(E);

  RootNbr alpha_simple = subs.parent_nr_simple(p.s0);
  const WeylWord to_simple = fixed_conjugate_simple(ctxt,alpha_simple);
  // by symmetry by $\delta$, |to_simple| conjugates $\delta(\alpha)$ to simple:
  assert(p.length()==1 or rd.is_simple_root(rd.permuted_root(to_simple,
				                subs.parent_nr_simple(p.s1))));
  // apply flip for $\delta$ acting on root set for |to_simple|, as elsewhere
  E.flip(ctxt.shift_flip(theta,new_theta,pos_to_neg(rd,to_simple)));

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
  const RootDatum& rd = E.rc().root_datum();
  return (E.gamma_lambda + shift).dot(rd.coroot(alpha));
}


// compute type of |p| for |E|, and export adjacent |paramin| values in |links|
DescValue star (const repr::Ext_common_context& ctxt,
		const paramin& E, const ext_gen& p,
		containers::sl_list<paramin>& links)
{
  paramin E0=E; // a copy of |E| that might be modified below to "normalise"
  DescValue result;

  const TwistedWeylGroup& tW = E.rc().twisted_Weyl_group();
  const InnerClass& ic = E.rc().inner_class();
  const InvolutionTable& i_tab = ic.involution_table();
  const RootDatum& rd = E.rc().root_datum();
  const RootDatum& integr_datum = ctxt.id();
  const SubSystem& subs = ctxt.subsys();
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
	const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	const auto theta_p = i_tab.nr(new_tw);
	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = ctxt.shift_flip(theta,theta_p,S);

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
	  if (not rd.is_simple_root(alpha_simple))
	  {
	    --tau_coef; // the parity change and decrease are both relevant
	    weyl::Generator s = // first switched root index
	      rd.find_descent(alpha_simple);
	    RootNbr first = // corresponding root summand, conjugated back
	      rd.permuted_root(rd.simpleRootNbr(s),ww);
	    assert(alpha == (E.ctxt.delta()+1)*rd.root(first));
	    new_gamma_lambda -= rd.root(first);
	    new_tau -= rd.root(first);
	  }

	  if (tau_coef%2!=0) // was set up so that this means: switched
	  { // no spurious $\tau'$ since $\<\alpha^\vee,(X^*)^\theta>=2\Z$:
#ifndef NDEBUG
	    auto ratv = (E.ctxt.delta()-1)*(E.gamma_lambda-rho_r_shift);
	    ratv.normalize(); // having a multiple of the weight here won't do!
	    const Weight target
	      { ratv.numerator().begin(),ratv.numerator().end() };
	    assert(not matreduc::has_solution (th_1,target));
#endif
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
	const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	const TwistedInvolution new_tw = // downstairs
	  tW.prod(subs.reflection(p.s0),E.tw);

	const auto theta_p=i_tab.nr(new_tw);
	const auto S = pos_to_neg(rd,ww);
	Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = ctxt.shift_flip(theta,theta_p,S);
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
	links.push_back(complex_cross(ctxt,p,E));
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
	const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	const auto theta_p = i_tab.nr(new_tw); // upstairs

	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = ctxt.shift_flip(theta,theta_p,S);
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
		     E.gamma_lambda - rho_r_shift - alpha*m,
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
	const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	assert(rd.is_simple_root(alpha_simple)); // no complications here
	const TwistedInvolution new_tw = // downstairs
	  tW.prod(subs.reflection(p.s1),tW.prod(subs.reflection(p.s0),E.tw));

	const auto theta_p=i_tab.nr(new_tw);
	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = ctxt.shift_flip(theta,theta_p,S);
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
	  links.push_back(complex_cross(ctxt,p,E));
	}
	else if (ascent)
	{ // twisted commutation with |s0.s1|: 2Ci
	  result = two_semi_imaginary;

	  TwistedInvolution new_tw = E.tw;
	  tW.twistedConjugate(subs.reflection(p.s0),new_tw); // same for |p.s1|

	  RootNbr alpha_simple = n_alpha;
	  const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const auto theta_p = i_tab.nr(new_tw); // upstairs
	  const auto S = pos_to_neg(rd,ww);
	  const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	  bool flipped = ctxt.shift_flip(theta,theta_p,S);
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
	  const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const auto theta_p=i_tab.nr(new_tw);
	  const auto S = pos_to_neg(rd,ww);
	  const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	  bool flipped = ctxt.shift_flip(theta,theta_p,S);
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
	const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	const auto theta_p = i_tab.nr(new_tw); // upstairs

	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = ctxt.shift_flip(theta,theta_p,S);
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
	const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	assert(rd.is_simple_root(alpha_simple)); // no complications here

	const auto theta_p=i_tab.nr(new_tw);
	const auto S = pos_to_neg(rd,ww);
	const Weight rho_r_shift = ctxt.to_simple_shift(theta,theta_p,S);
	bool flipped = ctxt.shift_flip(theta,theta_p,S);
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
	  const WeylWord ww = fixed_conjugate_simple(ctxt,alpha_simple);
	  assert(rd.is_simple_root(alpha_simple)); // no complications here

	  const auto theta_p=i_tab.nr(new_tw);
	  const auto S = pos_to_neg(rd,ww);
	  const Weight rho_r_shift = ascent
	    ?  ctxt.to_simple_shift(theta,theta_p,S)
	    : -ctxt.to_simple_shift(theta,theta_p,S);
	  bool flipped = ctxt.shift_flip(theta,theta_p,S);
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
	  links.push_back(complex_cross(ctxt,p,E));
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
  (const blocks::common_block& block, const WeightInvolution& delta)
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

} // |ext_block::ext_block|, from a |common_block|

bool ext_block::tune_signs(const blocks::common_block& block)
{
  repr::Ext_common_context ctxt
    (block.context().real_group(),delta(),block.integral_subsystem());
  repr::Ext_rep_context param_ctxt(block.context(),delta());
  containers::sl_list<paramin> links;
  for (BlockElt n=0; n<size(); ++n)
  { auto z=this->z(n);
    const paramin E(param_ctxt,block.x(z),block.gamma_lambda(z));
    for (weyl::Generator s=0; s<rank(); ++s)
    { const ext_gen& p=orbit(s); links.clear(); // output arguments for |star|
      auto tp = star(ctxt,E,p,links);
      if (might_be_uncertain(descent_type(s,n)) and
	  data[s][n].links.first==UndefBlock) // then reset the uncertain type
      {
	data[s][n].type=tp; // replace type, leave possible links to |UndefBlock|
	continue; // so in this case there are no link signs to tune
      }
      else if (tp!=descent_type(s,n))
	return false; // something is wrong

      switch (tp)
      {
      case one_imaginary_pair_switched: case one_real_pair_switched:
      case one_real_nonparity: case one_imaginary_compact:
      case two_imaginary_single_double_switched:
      case two_real_single_double_switched:
      case two_real_nonparity: case two_imaginary_compact:
      case three_real_nonparity: case three_imaginary_compact:
	assert(links.empty());
	break;

      case one_complex_ascent: case one_complex_descent:
      case two_complex_ascent: case two_complex_descent:
      case three_complex_ascent: case three_complex_descent:
	{ assert(links.size()==1);
	  const paramin q = *links.begin();
	  BlockElt m=cross(s,n); // cross neighbour as bare element of |*this|
	  if (m==UndefBlock)
	    break; // don't fall off the edge of a partial block
	  BlockElt cz = this->z(m); // corresponding element of (parent) |block|
	  paramin F(param_ctxt,block.x(cz),block.gamma_lambda(cz)); // default
	  assert(same_standard_reps(q,F)); // must lie over same parameter
	  if (not same_sign(q,F))
	    flip_edge(s,n,m);
	}
	break;

      case one_imaginary_single: case one_real_single:
      case two_imaginary_single_single: case two_real_single_single:
	{ assert(links.size()==2);
	  const paramin q0 = *links.begin();
	  const paramin q1 = *std::next(links.begin());
	  BlockElt m=some_scent(s,n); // the unique (inverse) Cayley
	  if (m!=UndefBlock) // don't fall off the edge of a partial block
	  {
	    BlockElt Cz = this->z(m); // corresponding element of block
	    paramin F(param_ctxt,block.x(Cz),block.gamma_lambda(Cz));
	    assert(same_standard_reps(q0,F));
	    if (not same_sign(q0,F))
	      flip_edge(s,n,m);
	  }
	  if ((m=cross(s,n))!=UndefBlock) // cross link, don't fall off the edge
	  {
	    BlockElt cz = this->z(m);
	    paramin Fc(param_ctxt,block.x(cz),block.gamma_lambda(cz));
	    assert(same_standard_reps(q1,Fc));
	    if (not same_sign(q1,Fc))
	      flip_edge(s,n,m);
	  }
	} break;
      case two_semi_imaginary: case two_semi_real:
      case three_semi_imaginary: case three_real_semi:
      case three_imaginary_semi: case three_semi_real:
	{ assert(links.size()==1);
	  const paramin q = *links.begin();
	  BlockElt m=some_scent(s,n); // the unique (inverse) Cayley
	  if (m==UndefBlock)
	    break; // don't fall off the edge of a partial block
	  BlockElt Cz = this->z(m); // corresponding element of block
	  paramin F(param_ctxt,block.x(Cz),block.gamma_lambda(Cz));
	  assert(same_standard_reps(q,F));
	  if (not same_sign(q,F))
	    flip_edge(s,n,m);
	}
	break;

      case one_imaginary_pair_fixed: case one_real_pair_fixed:
      case two_imaginary_double_double: case two_real_double_double:
	{ assert(links.size()==2);
	  const paramin q0 = *links.begin();
	  const paramin q1 = *std::next(links.begin());
	  BlockEltPair m=Cayleys(s,n);

	  if (m.first==UndefBlock)
	    break; // nothing to do if both are undefined

	  BlockElt Cz = this->z(m.first);
	  paramin F0(param_ctxt,block.x(Cz),block.gamma_lambda(Cz));
	  bool straight = same_standard_reps(q0,F0);
	  const auto& node0 = straight ? q0 : q1;
	  assert(same_standard_reps(node0,F0));
	  if (not same_sign(node0,F0))
	    flip_edge(s,n,m.first);

	  if (m.second==UndefBlock)
	    break;

	  Cz = this->z(m.second);
	  paramin F1(param_ctxt,block.x(Cz),block.gamma_lambda(Cz));
	  const auto& node1 = straight ? q1 : q0;
	  assert(same_standard_reps(node1,F1));
	  if (not same_sign(node1,F1))
	    flip_edge(s,n,m.second);
	}
	break;

      case two_imaginary_single_double_fixed: case two_real_single_double_fixed:
	{ assert(links.size()==2);
	  const paramin q0 = *links.begin();
	  const paramin q1 = *std::next(links.begin());
	  BlockEltPair m=Cayleys(s,n);

	  if (m.first==UndefBlock)
	    break; // nothing to do if both are undefined

	  BlockElt Cz = this->z(m.first);
	  paramin F0(param_ctxt,block.x(Cz),block.gamma_lambda(Cz));
	  bool straight=same_standard_reps(q0,F0);
	  const auto& node0 = straight ? q0 : q1;
	  assert(same_standard_reps(node0,F0));
	  if (not same_sign(node0,F0))
	    flip_edge(s,n,m.first);

	  if (m.second==UndefBlock)
	    break;

	  Cz= this->z(m.second);
	  paramin F1(param_ctxt,block.x(Cz),block.gamma_lambda(Cz));
	  const auto& node1 = straight ? q1 : q0;
	  assert(same_standard_reps(node1,F1));
	  if (not same_sign(node1,F1))
	    flip_edge(s,n,m.second);
	}
	break;
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


/*
  This function serves to replace and circumvent |Rep_context::make_dominant|
  applied to a scaled parameter (as occurs in the ordinary deformation function
  by calling |finals_for| or |normalise|, both of which call |make_dominant|,
  after calling |scale|), where |make_dominant| maps any ordinary parameter to
  one with a dominant |gamma| component, and moreover descends through singular
  complex descents in the block to the lowest parameter equivalent to the
  initial parameter. The reason that this is necessary is that scaling only
  affects the |nu| component of the infinitesimal character, so it may make it
  traverse walls of Weyl chambers. Indeed the caller should make sure |sr|
  itself has dominant |gamma|, which moreover is assumed to be fixed by |delta|
  (if not, don't use this function).

  The difference with the functioning of |make_dominant| is that here we keep
  track of all extended parameter components inherited from |sr| (so before
  scaling its |nu| part by |factor|), transforming them from the default choices
  for |sr|, and at the end comparing the transformed values to the default
  choices at the final parameter reached, recording the sign in |flipped|.
 */
StandardRepr scaled_extended_dominant // result will have its |gamma()| dominant
(const Rep_context rc,
 const StandardRepr& sr, const WeightInvolution& delta,
 Rational factor, // |z.nu()| is scaled by |factor| first
 bool& flipped // records whether a net extended flip was computed
 )
{
  const RootDatum& rd=rc.root_datum();
  RealReductiveGroup& G=rc.real_group(); const KGB& kgb = rc.kgb();
  repr::Ext_rep_context ctxt(rc,delta);
  const ext_gens orbits = rootdata::fold_orbits(rd,delta);
  assert(is_dominant_ratweight(rd,sr.gamma())); // dominant
  assert(((delta-1)*sr.gamma().numerator()).isZero()); // $\delta$-fixed

  // first approximation to result is scaled input
  // importantly, $\lambda$ (or equivalently |lambda_rho|) is held fixed here
  auto scaled_sr = rc.sr(sr.x(),rc.lambda_rho(sr),sr.gamma()*factor);
  // it will be convenent to have a working (modifiable) copy of |gamma|
  RatWeight gamma = scaled_sr.gamma(); // a working copy
  KGBElt x = scaled_sr.x(); // another variable, for convenience

  paramin E0 = paramin::default_extend(ctxt,sr);

  E0.gamma_lambda += gamma-sr.gamma(); // shift |E0.gamma_lambda| by $\nu$ change

  int_Vector r_g_eval (rd.semisimpleRank()); // simple root evaluations at |-gr|
  { const RatCoweight& g_r=ctxt.g_rho_check();
    for (unsigned i=0; i<r_g_eval.size(); ++i)
      r_g_eval[i] = -g_r.dot(rd.simpleRoot(i));
  }

  { unsigned i; // index into |orbits|
    do // make |gamma_numer| dominant, uses only complex simple root reflections
      for (i=0; i<orbits.size(); ++i)
	if (kgb.status(x).isComplex(orbits[i].s0))
	{ const auto& s=orbits[i];
	  const auto& alpha_v = rd.simpleCoroot(s.s0);
	  if (alpha_v.dot(gamma.numerator())<0)
	  {
	    rd.act(s.w_kappa,gamma); // change infin.character representative
	    rd.act(s.w_kappa,E0.gamma_lambda);
	    rd.act(s.w_kappa,E0.tau);
	    rd.shifted_dual_act(E0.l,s.w_kappa,r_g_eval);
	    rd.dual_act(E0.t,s.w_kappa);
	    E0.flip(s.length()==2); // record flip for every 2C+/2C- done
	    x = kgb.cross(s.w_kappa,x);
	    break; // indicate we advanced; restart search for |s|
	  }
	} // |for(i)|, if |isComplex|
    while(i<orbits.size()); // continue until above |for| runs to completion
  } // end of transformation of extended parameter components

  // now ensure that |E| gets matching |gamma| and |theta| (for flipped test)
  paramin E1(ctxt,kgb.involution(x),E0.gamma_lambda,E0.tau,E0.l,E0.t,E0.flipped);

  { // descend through complex singular simple descents
    repr::Ext_common_context block_ctxt(G,delta, SubSystem::integral(rd,gamma));
    const auto& int_datum = block_ctxt.id();
    const ext_gens integral_orbits = rootdata::fold_orbits(int_datum,delta);
    const RankFlags singular_orbits = // flag singular among integral orbits
      reduce_to(integral_orbits,singular_generators(int_datum,gamma));
    /*
      |singular_orbits| are in fact orbits of simple roots, because we have
      ensured |gamma| is dominant, so the singular subsystem of |rd|
      is generated by simple (co)roots (evaluting to 0 on |ctxt.gamma()|)
    */
     // record the corresponding simple root indices in |rd|, in order
    containers::sl_list<unsigned> orbit_simple;
    for (auto it=singular_orbits.begin(); it(); ++it)
    { auto alpha=block_ctxt.subsys().parent_nr_simple(integral_orbits[*it].s0);
      orbit_simple.push_back(rd.simpleRootIndex(alpha));
    }

    while(true) // will |break| below if no singular complex descent exists
    {
      auto soit = singular_orbits.begin(); auto it=orbit_simple.begin();
      for (; not orbit_simple.at_end(it); ++it,++soit)
	if (kgb.isComplexDescent(*it,x))
	  break;

      assert((not soit())==orbit_simple.at_end(it));
      if (not soit()) // previous loop ran to completion
	break;

      // find orbit among |integral_orbits| corresponding to that ComplexDescent
      ext_gen p=integral_orbits[*soit];
      assert(block_ctxt.subsys().parent_nr_simple(p.s0)
	     ==rd.simpleRootNbr(*it)); // check that we located it

      containers::sl_list<paramin> links;
      auto type = // compute neighbours in extended block
	star(block_ctxt,E1,p,links);
      assert(is_complex(type) or type==two_semi_real);
      E1 = *links.begin(); // replace |E| by descended parameter
      E1.flip(has_october_surprise(type)); // to undo extra flip |star|
      assert(x>E1.x()); // make sure we advance; we did simple complex descents
      x = E1.x(); // adapt |x| for complex descent test
    } // |while| a singular complex descent exists
  }

  // finally extract |StandardRepr| from |E|
  StandardRepr result = E1.restrict(gamma);

  // but the whole point of this function is to record the relative flip too!
  flipped = // compare |E1| to default
    not same_sign(E1,paramin::default_extend(ctxt,result));
  return result;

} // |scaled_extended_dominant|

/*
  The following function determines whether an extended parameter has a descent
  for generator |kappa|, which is an orbit of singularly-simple roots, and since
  |gamma| is supposed dominant here these are actually simple roots. Rather than
  call |star| here to do the full analysis, we can do a simplified (there is no
  |to_simple_shift|) test of notably the parity condition in the real case.
 */
bool is_descent
  (const repr::Ext_common_context& ctxt, const ext_gen& kappa, const paramin& E)
{ // easy solution would be to |return is_descent(star(E,kappa,dummy))|;
  const InvolutionTable& i_tab = E.rc().inner_class().involution_table();
  const InvolutionNbr theta = i_tab.nr(E.tw); // so use root action of |E.tw|
  const RootNbr n_alpha = ctxt.subsys().parent_nr_simple(kappa.s0);
  const RootNbr theta_alpha = i_tab.root_involution(theta,n_alpha);
  const RootDatum& rd = E.rc().root_datum();
  assert(rd.is_simple_root(n_alpha)); // as explained in the comment above

  // we don't need to inspect |kappa.type|, it does not affect descent status
  if (theta_alpha==n_alpha) // imaginary case, return whether compact
    return (E.ctxt.g_rho_check()-E.l).dot(rd.root(n_alpha)) %2!=0;
  if (theta_alpha==rd.rootMinus(n_alpha)) // real, return whether parity
    // whether $\<\alpha^\vee,\gamma-\lambda>$ even (since $\alpha$ will be
    // singular here, |gamma_lambda| gives the same result as |lambda| would)
    return E.gamma_lambda.dot(rd.coroot(n_alpha)) %2==0;
  else // complex
    return rd.is_negroot(theta_alpha);
} // |is_descent|

weyl::Generator first_descent_among
  (const repr::Ext_common_context& ctxt, RankFlags singular_orbits,
   const ext_gens& orbits, const paramin& E)
{ for (auto it=singular_orbits.begin(); it(); ++it)
    if (is_descent(ctxt,orbits[*it],E))
      return *it;
  return orbits.size(); // no singular descents found
}

/*
 This function is destined to be used after |scaled_extended_dominant|, to
 express the standard representation as linear combination of ones without
 singular descents, while keeping track (unlike |Rep_context::expand_final|)
 of flips that might occur during the process.
 */
containers::sl_list<std::pair<StandardRepr,bool> > extended_finalise
  (const repr::Rep_context& rc,
   const StandardRepr& sr, const WeightInvolution& delta)
{ // in order that |singular_generators| generate the whole singular system:
  assert(is_dominant_ratweight(rc.root_datum(),sr.gamma()));
  // we must assume |gamma| already dominant, DON'T call |make_dominant| here!

  repr::Ext_rep_context param_ctxt(rc,delta);
  repr::Ext_common_context ctxt
    (rc.real_group(),delta,SubSystem::integral(rc.root_datum(),sr.gamma()));

  const ext_gens orbits = rootdata::fold_orbits(ctxt.id(),delta);
  const RankFlags singular_orbits =
    reduce_to(orbits,singular_generators(ctxt.id(),sr.gamma()));

  containers::queue<paramin> to_do { paramin::default_extend(param_ctxt,sr) };
  containers::sl_list<std::pair<StandardRepr,bool> > result;

  do
  { const paramin E= to_do.front();
    to_do.pop(); // we are done with |head|
    auto s = first_descent_among(ctxt,singular_orbits,orbits,E);
    if (s>=orbits.size()) // no singular descents, so append to result
      result.emplace_back
	(std::make_pair(E.restrict(sr.gamma()),not is_default(E)));
    else // |s| is a singular descent orbit
    { containers::sl_list<paramin> links;
      auto type = star(ctxt,E,orbits[s],links);
      if (not is_like_compact(type)) // some descent, push to front of |to_do|
      { bool flip = has_october_surprise(type); // to undo extra flip |star|
	auto l_it=links.begin(); RootNbr witness;
	if (rc.is_nonzero(l_it->restrict(sr.gamma()),witness))
	{
	  l_it->flip(flip);
	  to_do.push(*l_it);
	}
	if (has_double_image(type)  // then maybe add second node after |head|
	    and rc.is_nonzero((++l_it)->restrict(sr.gamma()),witness))
	{ l_it->flip(flip);
	  to_do.push(*l_it);
	}
      }
    }
  }
  while(not to_do.empty());

  return result;
} // |extended_finalise|

} // |namespace ext_block|

} // |namespace atlas|
