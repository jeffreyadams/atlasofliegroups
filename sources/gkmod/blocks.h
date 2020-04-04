/*
  This is blocks.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2007-2016 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Class definition and function declarations for class |Block| and friends.

#ifndef BLOCKS_H  /* guard against multiple inclusions */
#define BLOCKS_H

#include <cassert>
#include <iostream>

#include "ratvec.h"	// containment infinitesimal character

#include "../Atlas.h"
#include "tits.h"	// representative of $y$ in |non_integral_block|
#include "descents.h"	// inline methods
#include "lietype.h"    // |ext_gen|;
#include "dynkin.h"     // |DynkinDiagram|

namespace atlas {

namespace blocks {


/******** function declarations *********************************************/

  // compute the involution in |dual_W| corresponding to |w| in |W|
  TwistedInvolution dual_involution
    (const TwistedInvolution& w,
     const TwistedWeylGroup& W,
     const TwistedWeylGroup& dual_W);

  // map from numbering of |b| to that of |dual_b|, assuming latter is dual
  std::vector<BlockElt> dual_map(const Block& b, const Block& dual_b);

  // flag intersection of |GR.Cartan_set()| with duals of Cartans from |dGR|
  BitMap common_Cartans(RealReductiveGroup& GR,	RealReductiveGroup& dGR);


/******** type definitions **************************************************/

// The class |BlockBase| serves external functionality, not block construction
class Block_base
{
public: // this |struct| must be public, though mainly used in derived classes
  struct EltInfo // per block element information
  {
    KGBElt x,y; // indices into |KGB| sets (which might no longer exist)
    DescentStatus descent;
    unsigned short length;
    BlockElt dual; // number of Hermitian dual of this element, if any
    EltInfo(KGBElt xx,KGBElt yy,DescentStatus dd, unsigned short ll)
      : x(xx),y(yy),descent(dd),length(ll), dual(UndefBlock) {}

  // sometimes leave |descent| and |length| (which |hashCode| ignores) blank
    EltInfo(KGBElt xx,KGBElt yy)
      : x(xx),y(yy),descent(),length(0), dual(UndefBlock) {}

  // methods that will allow building a hashtable with |info| as pool
    typedef std::vector<EltInfo> Pooltype;
    size_t hashCode(size_t modulus) const { return (13*x+21*y)&(modulus-1); }
    bool operator != (const EltInfo& o) const
    { return x!=o.x or y!=o.y; }

  }; // |struct EltInfo|

protected: // all fields may be set in a derived class contructor
  struct block_fields // per block element and simple reflection data
  {
    BlockElt cross_image;
    BlockEltPair Cayley_image;
    block_fields()
      : cross_image(UndefBlock), Cayley_image(UndefBlock,UndefBlock) {}
  };

  std::vector<EltInfo> info; // its size defines the size of the block
  std::vector<std::vector<block_fields> > data; // size |rank| * |size()|
  ext_gens orbits; // orbits of simple generators under distinguished involution

  DynkinDiagram dd; // diagram on simple generators for the block

  // possible tables of Bruhat order and Kazhdan-Lusztig polynomials
  BruhatOrder* d_bruhat;
  kl::KLContext* klc_ptr;

public:
// constructors and destructors
  Block_base(const KGB& kgb); // for |Block|, implicitly at integral inf. char.
  Block_base(unsigned int integral_rank); // only dimensions some vectors

  virtual ~Block_base(); // deletes |d_bruhat| and |klc_ptr| (if non-NULL)

// copy, assignment and swap

  Block_base& operator=(const Block_base& b) = delete;
protected: // derived classes may need to copy-construct their base
  Block_base(const Block_base& b);
public:

// accessors

  unsigned int rank() const { return data.size(); } // integral rank
  unsigned int folded_rank() const { return orbits.size(); }
  BlockElt size() const { return info.size(); }

  virtual KGBElt max_x() const = 0; // used virtually mainly for printing
  virtual KGBElt max_y() const = 0; // used virtually mainly for printing

  const DynkinDiagram& Dynkin() const { return dd; }
  ext_gen orbit(weyl::Generator s) const { return orbits[s]; }
  const ext_gens& inner_fold_orbits() const { return orbits; }

  KGBElt x(BlockElt z) const { assert(z<size()); return info[z].x; }
  KGBElt y(BlockElt z) const { assert(z<size()); return info[z].y; }

  size_t length(BlockElt z) const { return info[z].length; }

  // first element of length (at least) |l|, or |size()| if there are none
  BlockElt length_first(size_t l) const;

  BlockElt cross(weyl::Generator s, BlockElt z) const
  { assert(z<size()); assert(s<rank()); return data[s][z].cross_image; }

  const BlockEltPair& any_Cayleys(weyl::Generator s, BlockElt z) const
  { assert(z<size()); assert(s<rank()); return data[s][z].Cayley_image; }

  BlockEltPair cayley(weyl::Generator s, BlockElt z) const
  { assert(z<size()); assert(s<rank());
    if (not isWeakDescent(s,z))
      return data[s][z].Cayley_image;
    else return BlockEltPair(UndefBlock,UndefBlock);
  }

  BlockEltPair inverseCayley(weyl::Generator s, BlockElt z) const
  { assert(z<size()); assert(s<rank());
    if (isWeakDescent(s,z))
      return data[s][z].Cayley_image;
    else return BlockEltPair(UndefBlock,UndefBlock);
  }

  const DescentStatus& descent(BlockElt z) const
    { assert(z<size()); return info[z].descent; }
  DescentStatus::Value descentValue(weyl::Generator s, BlockElt z) const
    { assert(z<size()); assert(s<rank()); return descent(z)[s]; }

  bool isWeakDescent(weyl::Generator s, BlockElt z) const
    { return DescentStatus::isDescent(descentValue(s,z)); }

  // the following indicate existence of ascending/descending link
  bool isStrictAscent(weyl::Generator, BlockElt) const;
  bool isStrictDescent(weyl::Generator, BlockElt) const;
  weyl::Generator firstStrictDescent(BlockElt z) const;
  weyl::Generator firstStrictGoodDescent(BlockElt z) const;

  BlockElt Hermitian_dual(BlockElt z) const { return info[z].dual; }

  // The functor $T_{\alpha,\beta}$; might have been a non-method function
  BlockEltPair link
    (weyl::Generator alpha,weyl::Generator beta,BlockElt y) const;

  // print whole block to stream (name chosen to avoid masking by |print|)
  std::ostream& print_to
    (std::ostream& strm,bool as_invol_expr) const; // defined in |block_io|

  // print derivated class specific information  for |z| (used in |print_to|)
  virtual std::ostream& print
    (std::ostream& strm, BlockElt z,bool as_invol_expr) const =0;

  // manipulators
  BruhatOrder& bruhatOrder() { fillBruhat(); return *d_bruhat; }
  kl::KLContext& klc(BlockElt last_y, bool verbose)
  { fill_klc(last_y,verbose); return *klc_ptr; }

 private:
  void fillBruhat();
  void fill_klc(BlockElt last_y,bool verbose);

}; // |class Block_base|


// a derived class with minimal implementation to be a concrete class
class Bare_block : public Block_base
{ KGBElt x_size,y_size;
public:
  Bare_block(const Block_base& block)
  : Block_base(block), x_size(block.max_x()+1), y_size(block.max_y()+1) {}
  Bare_block(unsigned int rank, KGBElt x_size, KGBElt y_size)
    : Block_base(rank), x_size(x_size), y_size(y_size) {}
  virtual KGBElt max_x() const { return x_size-1; }
  virtual KGBElt max_y() const { return y_size-1; }
  virtual std::ostream& print
    (std::ostream& strm, BlockElt z,bool as_invol_expr) const { return strm; }

  // pseudo constructors
  static Bare_block dual (const Block_base& block);

}; // |class Bare_block|

/*				|class Block|
  A block constructed from a real form and a dual real form in an inner class.

  [Fokko's original comments; they have little bearing on the class definition]
  For our fixed inner form, orbits of $K$ on $G/B$ are parametrized by classes
  of elements $x$ in $N_G(H).\delta$ (the normalizer in the non-identity
  component $G.\delta$ of the extended group $G^Gamma=G \dot\cup G.\delta$,
  where $\delta$ is (i.e., acts on $G$ as) an involution that itself normalizes
  $H$), modulo the \emph{conjugation} action of $H$. (Dangerous bend: this $H$
  conjugacy class of $x$ is a subset, usually proper, of the coset $xH$. The
  collection of all $x$ is therefore NOT a subset of the extended Weyl group
  $N(H)/H$, but something more subtle.) The requirement on $x$ is that it belong
  to the $G$-conjugacy class of strong involutions defining the inner form.

  Each $x$ therefore defines an involution $\theta_x$ of $H$. Data pertaining to
  the subset of $x$ with a fixed $\theta_x$ is stored in the |Fiber| class.

  A block is characterized by specifying also an inner form of the dual group
  $G^vee$. For this inner form, $K^vee$ orbits on $G^vee/B^vee$ are parametrized
  by elements $y$. The basic theorem is that the block of representations is
  parametrized by pairs $(x,y)$ as above, subject to the requirement that
  $theta_y$ is the negative transpose of $theta_x$.
*/
class Block : public Block_base
{
  const TwistedWeylGroup& tW; // reference is used here only for printing

  size_t xrange;
  size_t yrange;

  // wasteful fields, but we cannot lean on |KGB|, which might no longer exist
  std::vector<size_t> d_Cartan; // of size |size()|
  TwistedInvolutionList d_involution; // of size |size()|

  // map KGB element |x| to the first block element |z| with |this->x(z)>=x|
  BlockEltList d_first_z_of_x; // of size |xrange+1|

  // Flags the generators occurring in reduced expression for |d_involution|.
  std::vector<RankFlags> d_involutionSupport; // of size |size()|


// constructors and destructors
  // the main constructor is private to ensure consistency of twists of KGBs
  Block(const KGB& kgb,const KGB& dual_kgb);

 public:
  // use one of the following two pseudo contructors to build |Block| values
  static Block build // pseudo contructor with small (and forgotten) KGB sets
    (InnerClass&, RealFormNbr rf, RealFormNbr drf);

  static Block build // pseudo contructor with stored KGB sets
    (RealReductiveGroup& G_R, RealReductiveGroup& dG_R);

  ~Block() {}

// copy, assignment and swap
  Block(const Block& b); // copy contructor, must be accessible, but is unused
 private:
  Block& operator=(const Block& b); // we don't however need to assign
 public:

// accessors

  const TwistedWeylGroup& twistedWeylGroup() const { return tW; }
  const WeylGroup& weylGroup() const { return tW.weylGroup(); }

  virtual KGBElt max_x() const { return xrange-1; }
  virtual KGBElt max_y() const { return yrange-1; }

  // Look up element by |x|, |y| coordinates
  BlockElt element(KGBElt x,KGBElt y) const;

  size_t Cartan_class(BlockElt z) const
    { assert(z<size()); return d_Cartan[z]; }

  size_t max_Cartan() const { return Cartan_class(size()-1); } // for printing

/*
  The twisted involution corresponding to the involution $\theta$ for |z|.
  This is the Weyl group element $w$, such that $\theta=w.\delta$
*/
  const TwistedInvolution& involution(BlockElt z) const
    { assert(z<size()); return d_involution[z]; }

  // flag among simple roots those occurring in reduced expr for |involution(z)|
  const RankFlags& involutionSupport(BlockElt z) const
  {
    assert(z<size());
    return d_involutionSupport[z];
  }

  virtual std::ostream& print // defined in block_io.cpp
   (std::ostream& strm, BlockElt z,bool as_invol_expr) const;


  // private accessor and manipulators
private:
  void compute_first_zs(); // set |first_z_of_x| according to |x| values
  void compute_supports(); // used during construction

}; // |class Block|

struct param_entry
{ KGBElt x; TorusPart y;

  // obligatory fields for hashable entry
  typedef std::vector<param_entry> Pooltype;
  size_t hashCode(size_t modulus) const
  { return (5*x-11*y.data().to_ulong())&(modulus-1); }
  bool operator !=(const param_entry& o) const { return x!=o.x or y!=o.y; }

}; // |struct param_entry|

typedef HashTable<y_entry,KGBElt> y_part_hash;
typedef Block_base::EltInfo block_elt_entry;
typedef HashTable<block_elt_entry,BlockElt> block_hash;
typedef HashTable<param_entry,BlockElt> param_hash;

// a class for blocks of (possibly non integral) parameters
class param_block : public Block_base
{
  const Rep_context& rc; // accesses many things, including KGB set for x

  RatWeight infin_char; // infinitesimal character
  Weight gr_numer; // numerator of |gamma-rho|

  std::vector<TorusPart> y_bits; // as in |StandardRepr|, indexed by |y|

  // hash structure to allow rapid lookup of |StandardRepr| values
  param_entry::Pooltype z_pool;
  param_hash z_hash;

  // group small components together:
  int gr_denom;
  KGBElt highest_x,highest_y; // maxima over this (maybe partial) block
  RankFlags singular; // flags, among integrally-simple roots, singular ones

 public:

  param_block // constructor for full block
    (const repr::Rep_context& rc,
     StandardRepr sr, // by value,since it will be made dominant before use
     BlockElt& entry_element	// set to block element matching input
    );

  param_block // alternative constructor, for interval below |sr|
    (const repr::Rep_context& rc,
     StandardRepr sr); // by value,since it will be made dominant before use

 public:
  // accessors that get values via |rc|
  const repr::Rep_context& context() const { return rc; }
  const RootDatum& root_datum() const;
  const InnerClass& inner_class() const;
  const InvolutionTable& involution_table() const;
  RealReductiveGroup& real_group() const;

  const RatWeight& gamma() const { return infin_char; }
  StandardRepr sr(BlockElt z) const; // parameter associated to block element

  RatWeight nu(BlockElt z) const; // "real" projection of |infin_char|
  Weight lambda_rho(BlockElt z) const; // reconstruct from |y_bits| value
  RatWeight lambda(BlockElt z) const; // reconstruct from y value
  RankFlags singular_simple_roots() const { return singular; }
  bool survives(BlockElt z) const; // whether $J(z_{reg})$ survives tr. functor
  containers::sl_list<BlockElt> finals_for(BlockElt z) const; // expr. for $I(z)$

  BlockElt lookup(const StandardRepr& sr) const;

  ext_gens fold_orbits(const WeightInvolution& delta) const;

  // virtual methods
  virtual KGBElt max_x() const { return highest_x; } // might not be final |x|
  virtual KGBElt max_y() const { return highest_y; }

  virtual std::ostream& print // defined in block_io.cpp
    (std::ostream& strm, BlockElt z,bool as_invol_expr) const;


 private:
  void compute_y_bits(const y_entry::Pooltype& y_pool); // set all the |y_bits|
  void compute_duals
  (const y_part_hash& y_hash,const block_hash& hash,
   const InnerClass& G,const SubSystem& rs);

/*
  reverse lengths and order block with them increasing, and by increasing
  |x(z)| among elements of given length; adapt tables accordingly. Argument
  |full_block| tells whether a full (rather than partial) block was generated;
  in this case (only) the |data| fields are adapted to the permutation
 */
  void reverse_length_and_sort(bool full_block);

}; // |class param_block|



class nblock_elt // internal representation during construction
{
  friend class nblock_help;
  KGBElt xx; // identifies element in parent KGB set
  TorusElement yy; // adds "local system" information to |xx|
public:
  nblock_elt (KGBElt x, const TorusElement& t) : xx(x), yy(t) {}

  KGBElt x() const { return xx; }
  const TorusElement& y() const { return yy; }

}; // |class nblock_elt|

class nblock_help // a support class for |nblock_elt|
{
public: // references stored for convenience, no harm in exposing them
  const KGB& kgb;
  const RootDatum& rd;  // the full (parent) root datum
  const SubSystem& sub; // the relevant subsystem
  const InvolutionTable& i_tab; // information about involutions, for |pack|

private:
  std::vector<TorusPart> dual_m_alpha; // the simple roots, reduced modulo 2
  std::vector<TorusElement> half_alpha; // half the simple roots

  void check_y(const TorusElement& t, InvolutionNbr i) const;
  void parent_cross_act(nblock_elt& z, weyl::Generator s) const;
  void parent_up_Cayley(nblock_elt& z, weyl::Generator s) const;
  void parent_down_Cayley(nblock_elt& z, weyl::Generator s) const;

public:
  nblock_help(RealReductiveGroup& GR, const SubSystem& subsys);

  void cross_act(nblock_elt& z, weyl::Generator s) const;
  void cross_act_parent_word(const WeylWord& ww, nblock_elt& z) const;
  void do_up_Cayley (nblock_elt& z, weyl::Generator s) const;
  void do_down_Cayley (nblock_elt& z, weyl::Generator s) const;
  bool is_real_nonparity(nblock_elt z, weyl::Generator s) const; // by value

  void twist(nblock_elt& z) const;

  y_entry pack_y(const nblock_elt& z) const;
}; // |class nblock_help|

} // |namespace blocks|

} // |namespace atlas|
#endif
