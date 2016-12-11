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


KGBElt twisted (const KGB& kgb, KGBElt x, const WeightInvolution& delta);

BlockElt twisted (const Block& block,
		  const KGB& kgb, const KGB& dual_kgb, // all are needed
		  BlockElt z,
		  const WeightInvolution& delta);

typedef Polynomial<int> Pol;

typedef bool (*extended_predicate)(DescValue);

class ext_block
{
  struct elt_info // per block element information
  {
    BlockElt z; // index into parent |Block_base| structure
    RankFlags flips[2];
#ifndef incompletecpp11
    elt_info(BlockElt zz): z(zz),flips{{},{}} {}
#else
  elt_info(BlockElt zz): z(zz),flips() {}
#endif


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
  ext_block(const InnerClass& G,
	    const param_block& block, const WeightInvolution& delta,
	    bool verbose=false);

// manipulators
  void flip_edge(weyl::Generator s, BlockElt x, BlockElt y);

  unsigned int list_edges();  // returns number of toggled pairs

// accessors

  size_t rank() const { return orbits.size(); }
  size_t size() const { return info.size(); }

  const Block_base& untwisted() const { return parent; }

  ext_gen orbit(weyl::Generator s) const { return orbits[s]; }
  const DynkinDiagram& Dynkin() const { return folded; }
  const WeightInvolution& delta() const { return d_delta; }

  BlockElt z(BlockElt n) const { assert(n<size()); return info[n].z; }

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
  void flip_edges(extended_predicate match);

}; // |class ext_block|

RankFlags reduce_to(const ext_gens orbits, RankFlags gen_set);

// Extended parameters

class context // holds values that remain fixed across extended block
{
  const repr::Rep_context& d_rc;
  WeightInvolution d_delta;
  RatWeight d_gamma; // representative of infinitesimal character
  RatCoweight d_g; // chosen lift of the common square for the square class
  RootDatum integr_datum; // intgrality datum
  SubSystem sub; // embeds |integr_datum| into parent root datum

 public:
  context
    (const repr::Rep_context& rc,
     WeightInvolution delta, // by value
     const RatWeight& gamma);

  const repr::Rep_context& rc () const { return d_rc; }
  const RootDatum& id() const { return integr_datum; }
  const SubSystem& subsys() const { return sub; }
  RealReductiveGroup& realGroup () const { return d_rc.realGroup(); }
  const InnerClass& innerClass () const { return realGroup().innerClass(); }
  const WeightInvolution& delta () const { return d_delta; }
  const RatWeight& gamma() const { return d_gamma; }
  const RatCoweight& g() const { return d_g; }
  // the next should match |realForm().g_rho_check()|, but uses stored |d_g|
  RatCoweight g_rho_check() const
  { return (g()-rho_check(rc().rootDatum())).normalize(); }
}; // |context|

// detailed parameter data; as defined by Jeff & David
struct param // prefer |struct| with |const| members for ease of access
{
  const context& ctxt;
  const TwistedInvolution tw; // implicitly defines $\theta$

private:
  Coweight d_l; // with |tw| gives a |GlobalTitsElement|; lifts its |t|
  Weight d_lambda_rho; // lift of that value in a |StandardRepr|
  Weight d_tau; // a solution to $(1-\theta)*\tau=(\delta-1)\lambda_\rho$
  Coweight d_t; // a solution to $t(1-theta)=l(\delta-1)$

public:
  param (const context& ec, const StandardRepr& sr);
  param (const context& ec, KGBElt x, const Weight& lambda_rho);
  param (const context& ec, const TwistedInvolution& tw,
	 Weight lambda_rho, Weight tau, Coweight l, Coweight t);

  param (const param& p) = default;
  param (param&& p)
  : ctxt(p.ctxt), tw(std::move(p.tw))
  , d_l(std::move(p.d_l))
  , d_lambda_rho(std::move(p.d_lambda_rho))
  , d_tau(std::move(p.d_tau))
  , d_t(std::move(p.d_t))
  {}

  param& operator= (const param& p)
  { assert(tw==p.tw); // cannot assign this, so it should match
    d_l=p.d_l; d_lambda_rho=p.d_lambda_rho; d_tau=p.d_tau; d_t=p.d_t;
    return *this;
  }
  param& operator= (param&& p)
  { assert(tw==p.tw); // cannot assign this, so it should match
    d_l=std::move(p.d_l); d_lambda_rho=std::move(p.d_lambda_rho);
    d_tau=std::move(p.d_tau); d_t=std::move(p.d_t);
    return *this;
  }

  const Coweight& l () const { return d_l; }
  const Weight& lambda_rho () const { return d_lambda_rho; }
  const Weight& tau () const { return d_tau; }
  const Coweight& t () const { return d_t; }

  void set_l (Coweight l) { d_l=l; }
  void set_lambda_rho (Weight lambda_rho) { d_lambda_rho=lambda_rho; }
  void set_tau (Weight tau) { d_tau=tau; }
  void set_t (Coweight t) { d_t=t; }

  const repr::Rep_context rc() const { return ctxt.rc(); }
  const WeightInvolution& delta () const { return ctxt.delta(); }
  const WeightInvolution& theta () const
    { return ctxt.innerClass().matrix(tw); }

  KGBElt x() const; // reconstruct |x| component
  repr::StandardRepr restrict() const // underlying unexteded representation
    { return ctxt.rc().sr_gamma(x(),lambda_rho(),ctxt.gamma()); }
}; // |param|

/* Try to conjugate |alpha| by product of folded-generators for the (full)
   root system of |c| to a simple root, and return the left-conjugating word
   that was applied. This may fail, if after some conjugation one ends up with
   the long root of a nontrivially folded A2 subsystem (in which case there
   cannot be any solution because |alpha| is fixed by the involution but none
   of the simple roots in its component of the root system is). In this case
   |alpha| is left as that non simple root, and the result conjugates to it.
 */
WeylWord fixed_conjugate_simple (const context& c, RootNbr& alpha);

// whether |E| and |F| lie over equivalent |StandrdRepr| values
bool same_standard_reps (const param& E, const param& F);
// whether |E| and |F| give same sign, assuming |same_standard_reps(E,F)|
bool same_sign (const param& E, const param& F);
inline int sign_between (const param& E, const param& F)
  { return same_sign(E,F) ? 1 : -1; }

inline bool is_default (const param& E)
{ return same_sign(E,param(E.ctxt,E.x(),E.lambda_rho())); }


// find out type of extended parameters, and push its neighbours onto |links|
DescValue star (const param& E, const ext_gen& p,
		containers::sl_list<std::pair<int,param> >& links);

bool is_descent (const ext_gen& kappa, const param& E);
weyl::Generator first_descent_among
  (RankFlags singular_orbits, const ext_gens& orbits, const param& E);

// a variation of |Rep_context::make_dominant|, used during extended deformation
StandardRepr scaled_extended_dominant // result will have its |gamma()| dominant
(const Rep_context rc, const StandardRepr& sr, const WeightInvolution& delta,
 Rational factor, // |z.nu()| is scaled by |factor| first
 bool& flipped // records whether and extended flip was recorded
 );

// expand parameter into a signed sum of extended nonzero final parameters
containers::sl_list<std::pair<StandardRepr,bool> > extended_finalise
  (const repr::Rep_context& rc,
   StandardRepr sr, // by value: internally |make_dominant| is applied to it
   const WeightInvolution& delta);

// check braid relation at |x|; also mark all involved elements in |cluster|
bool check_braid
  (const ext_block& b, weyl::Generator s, weyl::Generator t, BlockElt x,
   BitMap& cluster);

} // |namespace ext_block|

} // |namespace atlas|
#endif
