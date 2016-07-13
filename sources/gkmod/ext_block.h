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
  one_imaginary_pair_switched,  // type 2, twist-switched Cayley images
  one_real_pair_switched,       // type 1, twist-switched inverse Cayley images
  one_real_nonparity,
  one_imaginary_compact,

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
  two_imaginary_single_double_switched,  // 2i12, twist-switched images
  two_real_single_double_switched,	 // 2r12, twist-switched images
  two_real_nonparity,
  two_imaginary_compact,

  three_complex_ascent,         // distinct non-commuting complex ascents
  three_complex_descent,        // distinct non-commuting complex descents
  three_semi_imaginary,  // non-commuting complex ascents become imaginary
  three_real_semi, // non-commuting single-valued inverse Cayleys get complex
  three_imaginary_semi,  // non-commuting single-valued Cayleys become complex
  three_semi_real,       // non-commuting complex descents become real
  three_real_nonparity,
  three_imaginary_compact
}; // |enum DescValue|

const char* descent_code(DescValue v); // defined in |block_io|
inline bool is_descent(DescValue v) { return v%2!=0; }
bool is_complex(DescValue v);
bool is_unique_image(DescValue v);
bool has_double_image(DescValue v);
bool is_like_nonparity(DescValue v);
bool is_like_compact(DescValue v);
bool is_like_type_1(DescValue v);
bool is_like_type_2(DescValue v);
bool has_defect(DescValue v);
bool has_quadruple(DescValue v); // 2i12/2r21 cases

bool is_proper_ascent(DescValue v);

int generator_length(DescValue v);

DescValue extended_type(const Block_base& block, BlockElt z, const ext_gen& p,
			BlockElt& first_link);


KGBElt twisted (const KGB& kgb, KGBElt x,
		const WeightInvolution& delta, const weyl::Twist& twist);

BlockElt twisted (const Block& block,
		  const KGB& kgb, const KGB& dual_kgb, // all are needed
		  BlockElt z,
		  const WeightInvolution& delta,
		  const weyl::Twist& twist);

typedef Polynomial<int> Pol;

class ext_block
{
  struct elt_info // per block element information
  {
    BlockElt z; // index into parent |Block_base| structure
    unsigned short length; // length for extended group
    RankFlags flips[2];
    elt_info(BlockElt zz): z(zz),length(0),flips{{},{}} {}

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

  WeightInvolution d_delta;

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
	    const param_block& block, const KGB& kgb,
	    const WeightInvolution& delta);

// manipulators
  void flip_edge(weyl::Generator s, BlockElt x, BlockElt y);

  bool // return new value; true means edge was flipped to minus
    toggle_edge(BlockElt x,BlockElt y, bool verbose=true);
  bool set_edge(BlockElt x,BlockElt y);    // always set
  unsigned int list_edges();  // returns number of toggled pairs
  void report_2Ci_toggles() const;

// accessors

  size_t rank() const { return orbits.size(); }
  size_t size() const { return info.size(); }

  ext_gen orbit(weyl::Generator s) const { return orbits[s]; }
  const DynkinDiagram& Dynkin() const { return folded; }
  const WeightInvolution& delta() const { return d_delta; }

  BlockElt z(BlockElt n) const { assert(n<size()); return info[n].z; }

  // Look up element by its index in |parent| (if that did define an element)
  BlockElt element(BlockElt z) const; // partial inverse of method |z|

  const DescValue descent_type(weyl::Generator s, BlockElt n) const
    { assert(n<size()); assert(s<rank()); return data[s][n].type; }

  size_t length(BlockElt n) const { return info[n].length; }
  size_t l(BlockElt y, BlockElt x) const { return length(y)-length(x); }
  BlockElt length_first(size_t l) const { return l_start[l]; }

  // the following three function return ascents or descents as appropriate
  BlockElt cross(weyl::Generator s, BlockElt n) const;
  BlockElt Cayley(weyl::Generator s, BlockElt n) const; // just one or none
  BlockEltPair Cayleys(weyl::Generator s, BlockElt n) const; // must be two

  // some of the above: an (a/de)scent of |n| in block; assumed to exist
  BlockElt some_scent(weyl::Generator s, BlockElt n) const;

  // whether link for |s| from |x| to |y| has a sign flip attached
  int epsilon(weyl::Generator s, BlockElt x, BlockElt y) const;

  // coefficient of neighbour |sx| for $s$ in action $(T_s+1)*a_x$
  Pol T_coef(weyl::Generator s, BlockElt sx, BlockElt x) const;

  BlockEltList down_set(BlockElt y) const;

  // here all elements reached by a link are added to |l|, (a/de)scent first
  void add_neighbours(BlockEltList& dst, weyl::Generator s, BlockElt n) const;

  // print whole block to stream (name chosen to avoid masking by |print|)
  std::ostream& print_to(std::ostream& strm) const; // defined in |block_io|

private:
  void complete_construction(const BitMap& fixed_points);

}; // |class ext_block|



// Extended parameters

class context // holds values that remain fixed across extended block
{
  const repr::Rep_context& d_rc;
  WeightInvolution d_delta;
  RatWeight d_gamma; // representative of infinitesimal character
  RatCoweight d_g; // chosen lift of the common square for the square class
  RootDatum integr_datum; // intgrality datum
  SubSystem sub;

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

KGBElt x(const param& E); // reconstruct |E|

// whether |E| and |F| lie over equivalent |StandrdRepr| values
bool same_standard_reps (const param& E, const param& F);
// whether |E| and |F| give same sign, assuming |same_standard_reps(E,F)|
bool same_sign (const param& E, const param& F);
inline int sign_between (const param& E, const param& F)
  { return same_sign(E,F) ? 1 : -1; }

// find out type of extended parameters, and push its neighbours onto |links|
DescValue type (const param& E, const ext_gen& p,
		containers::sl_list<param>& links);

bool check(ext_block& eb, const param_block& block, bool verbose=false);

// check braid relation at |x|; also mark all involved elements in |cluster|
bool check_braid
  (const ext_block& b, weyl::Generator s, weyl::Generator t, BlockElt x,
   BitMap& cluster);

} // |namespace ext_block|

} // |namespace atlas|
#endif
