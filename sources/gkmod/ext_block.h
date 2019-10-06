/*
  This is ext_block.h

  Copyright (C) 2013-2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Class definition and function declarations for class |ext_block|.

#ifndef EXT_BLOCK_H  /* guard against multiple inclusions */
#define EXT_BLOCK_H

#include "../Atlas.h"

#include <cassert>
#include <iostream>
#include <set>

#include "lietype.h" // for the structure |ext_gen|
#include "dynkin.h"  // for |DynkinDiagram|
#include "innerclass.h"
#include "realredgp.h"
#include "blocks.h" // for some inlined methods (dependency should be removed)
#include "block_minimal.h" // for type |blocks::block_minimal|
#include "subsystem.h" // for inclusion of |SubSystem| field
#include "repr.h" // allows using |Rep_context| methods in this file

namespace atlas {

namespace ext_block {

// for correspondence of enumerations and conventional codes, see |descent_code|
enum DescValue // every even/odd pair is one of associated ascent and descent
{
  one_complex_ascent,
  one_complex_descent,
  one_imaginary_single,         // type 1, single Cayley image
  one_real_pair_fixed,          // type 1, twist-fixed inverse Cayley images
  one_imaginary_pair_fixed,     // type 2, twist-fixed Cayley images
  one_real_single,              // type 2, single inverse Cayley image
  one_imaginary_pair_switched,  // ascent:  type 1i2 with twist-switched images
  one_real_pair_switched,       // descent: type 1r1 with twist-switched images
  one_real_nonparity,		// an ascent
  one_imaginary_compact,	// a descent

  two_complex_ascent,           // distinct commuting complex ascents
  two_complex_descent,          // distinct commuting complex descents
  two_semi_imaginary,           // identical complex ascents
  two_semi_real,                // identical complex descents
  two_imaginary_single_single,  // commuting single-valued Cayleys
  two_real_double_double, // commuting double-valued inverse Cayleys (2-valued)
  two_imaginary_single_double_fixed,  // single-valued Cayleys become double
  two_real_single_double_fixed, // single-valued inverse Cayleys become double
  two_imaginary_double_double,  // commuting double-valued Cayleys (2-valued)
  two_real_single_single,	// commuting single-valued inverse Cayleys
  two_imaginary_single_double_switched,  // ascent: 2i12, twist-switched images
  two_real_single_double_switched,	 // descent:2r12, twist-switched images
  two_real_nonparity,		// an ascent
  two_imaginary_compact,	// a descent

  three_complex_ascent,         // distinct non-commuting complex ascents
  three_complex_descent,        // distinct non-commuting complex descents
  three_semi_imaginary,  // non-commuting complex ascents become imaginary
  three_real_semi, // non-commuting single-valued inverse Cayleys get complex
  three_imaginary_semi,  // non-commuting single-valued Cayleys become complex
  three_semi_real,       // non-commuting complex descents become real
  three_real_nonparity,		// an ascent
  three_imaginary_compact	// a descent
}; // |enum DescValue|

const char* descent_code(DescValue v); // defined in |block_io|
inline bool is_descent(DescValue v) { return v%2!=0; }
bool is_complex(DescValue v); // types 1C+, 1C-, 2C+, 2C-, 3C+, 3C-
bool is_unique_image(DescValue v); // types 1r1f, 1i2f, 2r11, 2i22, and defects
bool has_double_image(DescValue v); // types 1r1f, 1i2f, 2r11, 2i22, and quads
bool is_like_nonparity(DescValue v); // ascents with 0 links
bool is_like_compact(DescValue v);  // descents with 0 links
bool is_like_type_1(DescValue v); // types 1i1, 1r1f, 2i11, 2r11
bool is_like_type_2(DescValue v); // types 1i2f, 1r2, 2i22, 2r22
bool has_defect(DescValue v); // types 2Ci, 2Cr, 3Ci, 3r, 3Cr, 3i
bool has_quadruple(DescValue v); // types 2i12f and 2r21f
bool has_october_surprise(DescValue v); // types 2i**, 2r**, 2C+, 2C-, 3Ci, 3Cr

bool is_proper_ascent(DescValue v); // ascent with at least 1 link

unsigned int generator_length(DescValue v);
unsigned int link_count(DescValue v);

DescValue extended_type(const Block_base& block, BlockElt z, const ext_gen& p,
			BlockElt& first_link);


typedef Polynomial<int> Pol;

class ext_block
{
  struct elt_info // per block element information
  {
    BlockElt z; // index into parent |Block_base| structure
    RankFlags flips[2];
    elt_info(BlockElt zz): z(zz),flips{{},{}} {}

  // methods that will allow building a hashtable with |info| as pool
    typedef std::vector<elt_info> Pooltype;
    size_t hashCode(size_t modulus) const { return (-7*z)&(modulus-1); }
    bool operator != (const elt_info& o) const { return z!=o.z; }

  }; // |struct elt_info|

  struct block_fields // per block element and simple reflection data
  {
    DescValue type;
    BlockEltPair links; // one or two values, depending on |type|

  block_fields(DescValue t) : type(t),links(UndefBlock,UndefBlock) {}
  };

  const Block_base& parent; // block, maybe non-intg. where we are fixed points
  ext_gens orbits; // $\delta$-orbits of generators for the parent block
  const DynkinDiagram folded; // diagram defined on those orbits

  WeightInvolution d_delta; // purely passive: nothing else uses coordinates

  std::vector<elt_info> info; // its size defines the size of the block
  std::vector<std::vector<block_fields> > data;  // size |d_rank| * |size()|
  BlockEltList l_start; // where elements of given length start

 public:

// constructors and destructors
  ext_block(const InnerClass& G,
	    const Block& block,
	    const KGB& kgb, const KGB& dual_kgb, // all are needed
	    const WeightInvolution& delta);
  ext_block(const param_block& block, const WeightInvolution& delta,
	    bool verbose=false);
  ext_block(const blocks::block_minimal& block, const WeightInvolution& delta);

// manipulators
  void flip_edge(weyl::Generator s, BlockElt x, BlockElt y);

// accessors

  size_t rank() const { return orbits.size(); }
  size_t size() const { return info.size(); }

  const Block_base& untwisted() const { return parent; }

  ext_gen orbit(weyl::Generator s) const { return orbits[s]; }
  const DynkinDiagram& Dynkin() const { return folded; }
  const WeightInvolution& delta() const { return d_delta; }

  BlockElt z(BlockElt n) const { assert(n<size()); return info[n].z; }
  StandardRepr sr(BlockElt n, const param_block& parent) const
  { return parent.sr(z(n)); }

  // Look up element by its index in |parent| (if that did define an element)
  // more precisely returns smallest |n| with |z(n)>=z|, or |size()| if none
  BlockElt element(BlockElt z) const; // partial inverse of method |z|
  bool is_present(BlockElt zz) const
  { auto n=element(zz); return n<size() and z(n)==zz; }

  const DescValue descent_type(weyl::Generator s, BlockElt n) const
    { assert(n<size()); assert(s<rank()); return data[s][n].type; }

  unsigned length(BlockElt n) const;
  unsigned l(BlockElt y, BlockElt x) const { return length(y)-length(x); }
  BlockElt length_first(size_t l) const { return l_start[l]; }

  // the following three function return ascents or descents as appropriate
  BlockElt cross(weyl::Generator s, BlockElt n) const;
  BlockElt Cayley(weyl::Generator s, BlockElt n) const; // just one or none
  BlockEltPair Cayleys(weyl::Generator s, BlockElt n) const; // must be two

  // some of the above: an (a/de)scent of |n| in block; assumed to exist
  BlockElt some_scent(weyl::Generator s, BlockElt n) const;

  // whether link for |s| from |x| to |y| has a sign flip attached
  int epsilon(weyl::Generator s, BlockElt x, BlockElt y) const;

  // this works only for blocks generated from parameters, so supply |parent|
  // mark the extended generators (orbits) singular for |parent.gamma()|
  RankFlags singular_orbits(const param_block& parent) const;

  weyl::Generator first_descent_among
    (RankFlags singular_orbits, BlockElt y) const;

  // reduce a matrix to elements without descents among singular generators
  template<typename C> // matrix coefficient type (signed)
  containers::simple_list<BlockElt> // returns list of elements selected
    condense (matrix::Matrix<C>& M, const param_block& parent) const;
  template<typename C> // matrix coefficient type (signed)
  containers::simple_list<BlockElt> // returns list of elements selected
    condense (matrix::Matrix<C>& M, const blocks::block_minimal& parent,
	      const RatWeight& gamma) const;

  // coefficient of neighbour |xx| of |x| in the action $(T_s+1)*a_x$
  Pol T_coef(weyl::Generator s, BlockElt xx, BlockElt x) const;

  BlockEltList down_set(BlockElt y) const;

  // here all elements reached by a link are added to |l|, (a/de)scent first
  void add_neighbours(BlockEltList& dst, weyl::Generator s, BlockElt n) const;

  // print whole block to stream (name chosen to avoid masking by |print|)
  std::ostream& print_to(std::ostream& strm) const; // defined in |block_io|

private:
  void complete_construction(const BitMap& fixed_points);
  bool check(const param_block& block, bool verbose=false);
  bool tune_signs(const blocks::block_minimal& block);

}; // |class ext_block|

// reduce marked set |gen_set| (supposed a union of orbits) to set of |orbits|
RankFlags reduce_to(const ext_gens& orbits, RankFlags gen_set);

// Extended parameters

class context // holds values that remain fixed across extended block
{
  const repr::Rep_context& d_rc;
  const WeightInvolution d_delta;
  RatWeight d_gamma; // dominant representative of infinitesimal character
  // RatCoweight d_g; // we might record |g|, but in fact defer to |realGroup|
  const RootDatum integr_datum; // intgrality datum
  const SubSystem sub; // embeds |integr_datum| into parent root datum
  Permutation pi_delta; // permutation of |delta| on roots of full root datum
  RootNbrSet delta_fixed_roots;
  weyl::Twist twist;
  int_Vector lambda_shifts,l_shifts; // affine centers for complex cross actions

 public:
  context
    (const repr::Rep_context& rc,
     const WeightInvolution& delta,
     const RatWeight& gamma);

  // accessors
  const repr::Rep_context& rc () const { return d_rc; }
  const RootDatum& id() const { return integr_datum; }
  const SubSystem& subsys() const { return sub; }
  const RootDatum& rootDatum() const { return d_rc.rootDatum(); }
  const InnerClass& innerClass () const { return d_rc.innerClass(); }
  RealReductiveGroup& realGroup () const { return d_rc.realGroup(); }
  const WeightInvolution& delta () const { return d_delta; }
  const RatWeight& gamma() const { return d_gamma; }
  const RatCoweight& g_rho_check() const { return realGroup().g_rho_check(); }
  RatCoweight g() const { return realGroup().g(); }
  RootNbr delta_of(RootNbr alpha) const { return pi_delta[alpha]; }
  const RootNbrSet& delta_fixed() const { return delta_fixed_roots; }
  weyl::Generator twisted(weyl::Generator s) const { return twist[s]; }
  int lambda_shift(weyl::Generator s) const { return lambda_shifts[s]; }
  int l_shift(weyl::Generator s) const { return l_shifts[s]; }

  // whether positive $\alpha$ has $\theta(\alpha)\neq\pm(1|\delta)(\alpha)$
  bool is_very_complex (InvolutionNbr theta, RootNbr alpha) const;
  Weight to_simple_shift(InvolutionNbr theta, InvolutionNbr theta_p,
			 RootNbrSet pos_to_neg) const;
  bool shift_flip(InvolutionNbr theta, InvolutionNbr theta_p,
		  RootNbrSet pos_to_neg) const;

  // possible manipulator; |RootDatum|,|SubSystem| need to implement this first
  // void act_on_gamma(const WeylWord& ww); // left-apply |ww| to |d_gamma|

}; // |context|

// detailed parameter data; as defined by Jeff & David
struct param // allow public member access; methods ensure no invariants anyway
{
  const context& ctxt;
  TwistedInvolution tw; // implicitly defines $\theta$

  Coweight l; // with |tw| gives a |GlobalTitsElement|; lifts its |t|
  Weight lambda_rho; // lift of that value in a |StandardRepr|
  Weight tau; // a solution to $(1-\theta)*\tau=(\delta-1)\lambda_\rho$
  Coweight t; // a solution to $t(1-theta)=l(\delta-1)$
  bool flipped; // whether tensored with the fliiping representation

  param (const context& ec, const StandardRepr& sr, bool flipped=false);
  param (const context& ec,
	 KGBElt x, const Weight& lambda_rho, bool flipped=false);
  param (const context& ec, const TwistedInvolution& tw,
	 Weight lambda_rho, Weight tau, Coweight l, Coweight t,
	 bool flipped=false);

  param (const param& p) = default;
  param (param&& p)
  : ctxt(p.ctxt), tw(std::move(p.tw))
  , l(std::move(p.l))
  , lambda_rho(std::move(p.lambda_rho))
  , tau(std::move(p.tau))
  , t(std::move(p.t))
  , flipped(p.flipped)
  {}

  param& operator= (const param& p)
  { tw=p.tw;
    l=p.l; lambda_rho=p.lambda_rho; tau=p.tau; t=p.t;
    flipped=p.flipped;
    return *this;
  }
  param& operator= (param&& p)
  { tw=std::move(p.tw);
    l=std::move(p.l); lambda_rho=std::move(p.lambda_rho);
    tau=std::move(p.tau); t=std::move(p.t);
    flipped=p.flipped;
    return *this;
  }

  bool is_flipped() const { return flipped; }

  void flip (bool whether=true) { flipped=(whether!=flipped); }

  const repr::Rep_context rc() const { return ctxt.rc(); }
  const WeightInvolution& delta () const { return ctxt.delta(); }
  const WeightInvolution& theta () const
    { return ctxt.innerClass().matrix(tw); }

  KGBElt x() const; // reconstruct |x| component
  repr::StandardRepr restrict() const // underlying unexteded representation
    { return ctxt.rc().sr_gamma(x(),lambda_rho,ctxt.gamma()); }
}; // |param|


// a variation of |Rep_context::make_dominant|, used during extended deformation
StandardRepr scaled_extended_dominant // result will have its |gamma()| dominant
(const Rep_context rc, const StandardRepr& sr, const WeightInvolution& delta,
 Rational factor, // |z.nu()| is scaled by |factor| first
 bool& flipped // records whether and extended flip was recorded
 );

// expand parameter into a signed sum of extended nonzero final parameters
containers::sl_list<std::pair<StandardRepr,bool> > extended_finalise
  (const repr::Rep_context& rc,
   const StandardRepr& sr, const WeightInvolution& delta);

// check quadratic relation for |s| at |x|
bool check_quadratic (const ext_block& b, weyl::Generator s, BlockElt x);

// check braid relation at |x|; also mark all involved elements in |cluster|
bool check_braid
  (const ext_block& b, weyl::Generator s, weyl::Generator t, BlockElt x,
   BitMap& cluster);

} // |namespace ext_block|

} // |namespace atlas|
#endif
