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
#include <memory> // for |std::unique_ptr|

#include "../Atlas.h"
#include "ratvec.h"	// containment infinitesimal character
#include "sl_list.h"    // return type |down_set| function

#include "tits.h"	// representative of $y$ in |non_integral_block|
#include "descents.h"	// inline methods
#include "lietype.h"    // |ext_gen|;
#include "dynkin.h"     // |DynkinDiagram|
#include "subsystem.h"  // data memeber in |common_block|
#include "repr.h"       // hash table in |common_block|, |StandardReprMod|

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
    EltInfo(KGBElt xx,KGBElt yy,DescentStatus dd, unsigned short ll)
      : x(xx),y(yy),descent(dd),length(ll) {}

  // sometimes leave |descent| and |length| (which |hashCode| ignores) blank
    EltInfo(KGBElt xx,KGBElt yy)
      : x(xx),y(yy),descent(),length(0) {}

  // currently unused methods that allow building a hashtable around |info|
  // this used to be useful for lookup by |(x,y)| during block construction
    typedef std::vector<EltInfo> Pooltype;
    size_t hashCode(size_t modulus) const { return (13*x+21*y)&(modulus-1); }
    bool operator != (const EltInfo& o) const { return x!=o.x or y!=o.y; }

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
  std::vector<std::unique_ptr<BlockEltList> > partial_Hasse_diagram;
  std::unique_ptr<BruhatOrder> d_bruhat;
  std::unique_ptr<kl::KL_table> kl_tab_ptr;

public:
// constructors and destructors
  Block_base(const KGB& kgb); // for |Block|, implicitly at integral inf. char.
  Block_base(unsigned int integral_rank); // only dimensions some vectors

  virtual ~Block_base(); // out-of-line because of deleters implicitly called

// copy, assignment and swap

  Block_base& operator=(const Block_base& b) = delete;
protected: // derived classes may need to copy-construct their base
  Block_base(const Block_base& b);
public:

// accessors

  unsigned int rank() const { return data.size(); } // integral semisimple rank
  unsigned int folded_rank() const { return orbits.size(); }
  BlockElt size() const { return info.size(); }

  virtual KGBElt max_x() const = 0; // used virtually mainly for printing
  virtual KGBElt max_y() const = 0; // used virtually mainly for printing

  const DynkinDiagram& Dynkin() const { return dd; }
  ext_gen orbit(weyl::Generator s) const { return orbits[s]; }
  const ext_gens& inner_fold_orbits() const { return orbits; }

  KGBElt x(BlockElt z) const { assert(z<size()); return info[z].x; }
  KGBElt y(BlockElt z) const { assert(z<size()); return info[z].y; }

  unsigned short length(BlockElt z) const { return info[z].length; }

  // first element of length (at least) |l|, or |size()| if there are none
  BlockElt length_first(size_t l) const; // does a binary search in the block

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

  BlockElt unique_ascent(weyl::Generator s, BlockElt z) const
  { assert(z<size()); assert(s<rank());
    auto v = descentValue(s,z);
    if (v == DescentStatus::ComplexAscent)
      return data[s][z].cross_image;
    assert(v == DescentStatus::ImaginaryTypeI);
    return data[s][z].Cayley_image.first;
  }

  const DescentStatus& descent(BlockElt z) const
    { assert(z<size()); return info[z].descent; }
  DescentStatus::Value descentValue(weyl::Generator s, BlockElt z) const
    { assert(z<size()); assert(s<rank()); return info[z].descent[s]; }

  RankFlags descent_generators (BlockElt z) const; // all |s| giving weak descent
  bool isWeakDescent(weyl::Generator s, BlockElt z) const
    { return DescentStatus::isDescent(descentValue(s,z)); }

  // the following indicate existence of ascending/descending link
  bool isStrictAscent(weyl::Generator, BlockElt) const;
  bool isStrictDescent(weyl::Generator, BlockElt) const;
  weyl::Generator firstStrictDescent(BlockElt z) const;
  weyl::Generator firstStrictGoodDescent(BlockElt z) const;

  bool survives(BlockElt z, RankFlags singular) const;
    // whether $J(z_{reg})$ survives tr. functor with |singular| singular coroots
  containers::sl_list<BlockElt>
    finals_for(BlockElt z, RankFlags singular) const; // expression for $I(z)$

  // print whole block to stream (name chosen to avoid masking by |print|)
  std::ostream& print_to // defined in |block_io|
    (std::ostream& strm,bool as_invol_expr,RankFlags singular=RankFlags(0))
    const;

  // print derivated class specific information  for |z| (used in |print_to|)
  virtual std::ostream& print
    (std::ostream& strm, BlockElt z,bool as_invol_expr,RankFlags singular)
    const =0;

  // manipulators
  BruhatOrder& bruhatOrder() { fill_Bruhat(); return *d_bruhat; }
  BruhatOrder&& Bruhat_order() && { fill_Bruhat(); return std::move(*d_bruhat); }
  kl::KL_table& kl_tab
    (KL_hash_Table* pol_hash, BlockElt limit=0, bool verbose=false)
  { fill_kl_tab(limit,pol_hash,verbose); return *kl_tab_ptr; }

 protected:
  void set_Bruhat_covered (BlockElt z, BlockEltList&& covered);
 private:
  void fill_Bruhat();
  void fill_kl_tab(BlockElt limit, KL_hash_Table* pol_hash, bool verbose);

}; // |class Block_base|


// The functor $T_{\alpha,\beta}$
BlockEltPair link(weyl::Generator alpha,weyl::Generator beta,
		  const Block_base& block, BlockElt y);

// sorted list of elements reachable by a descent from $y$.
containers::simple_list<BlockElt> down_set(const Block_base& block,BlockElt y);


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
    (std::ostream& strm, BlockElt z,bool as_invol_expr,RankFlags singular) const
    { return strm; }

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
   (std::ostream& strm, BlockElt z,bool as_invol_expr,RankFlags singular) const;


  // private accessor and manipulators
private:
  void compute_first_zs(); // set |first_z_of_x| according to |x| values
  void compute_supports(); // used during construction

}; // |class Block|

// a class for blocks of (possibly non integral) parameters
class common_block : public Block_base
{
  const Rep_context& rc; // accesses many things, including KGB set for x

  const SubSystem integral_sys;

  // hash structure to facilitate lookup of elements in |StandardReprMod| form
  using repr_hash = HashTable<StandardReprMod,BlockElt>;
  StandardReprMod::Pooltype z_pool;
  repr_hash srm_hash;

  std::unique_ptr<ext_block::ext_block> extended;

  // group small data members together:
  KGBElt highest_x,highest_y; // maxima over this block
  const bool generated_as_full_block; // tells which constructor was used

 public:

  // constructor and destructor
  common_block // full block
    (const Rep_context& rc,
     const StandardReprMod& srm, // not modified, no "making dominant"
     BlockElt& entry_element	// set to block element matching input
    );
  common_block // partial block
    (const common_context& ctxt,
     containers::sl_list<StandardReprMod>& elements,
     const RatWeight& gamma_rep);
  ~common_block(); // cleans up |*extended|, so inline definition impossible

  // accessors that get values via |rc|
  const Rep_context& context() const { return rc; }
  const RootDatum& root_datum() const;
  const SubSystem& integral_subsystem() const { return integral_sys; }
  const InnerClass& inner_class() const;
  const InvolutionTable& involution_table() const;
  RealReductiveGroup& real_group() const;

  bool is_full () const { return generated_as_full_block; }

  // simple coroots of |sub| singular for |gamma|
  RankFlags singular (const RatWeight& gamma) const;

  // with |gamma| unknown, only the difference |gamma-lambda| is meaningful
  RatWeight gamma_lambda(BlockElt z) const;
  RatWeight gamma_lambda_rho(BlockElt z) const // that is $\gamma-\lambda+\rho$
  { return z_pool[z].gamma_rep(); } // is actually easier than |gamma_lambda|

  BlockElt lookup(const StandardReprMod& srm) const;
  BlockElt lookup(KGBElt x, RatWeight gamma_lambda) const; // by value

  const StandardReprMod& representative (BlockElt z) const { return z_pool[z]; }

  StandardRepr sr // reconstruct at |gamma| using |diff| of |gamma_rep|s
    (BlockElt z,const RatWeight& diff, const RatWeight& gamma) const;

  ext_gens fold_orbits(const WeightInvolution& delta) const;

  // manipulators
  // obtain KL hash table from |*kl_tab_ptr|, maybe creating it using arugment
  kl::Poly_hash_export KL_hash(KL_hash_Table* KL_pol_hash);
  void swallow // integrate an older partial block, with mapping of elements
    (common_block&& sub, const BlockEltList& embed,
     KL_hash_Table* KL_pol_hash, ext_KL_hash_Table* ext_KL_pol_hash);
  ext_block::ext_block& extended_block(ext_KL_hash_Table* pol_hash);
  // get/build extended block for inner class involution; if built, store it

  void set_Bruhat
  (containers::sl_list<std::pair<BlockElt,BlockEltList> >&& partial_Hasse);

  // virtual methods
  virtual KGBElt max_x() const { return highest_x; } // might not be final |x|
  virtual KGBElt max_y() const { return highest_y; }

  virtual std::ostream& print // defined in block_io.cpp
    (std::ostream& strm, BlockElt z,bool as_invol_expr,RankFlags singular) const;


 private:
// sort by increaing length, then |x|, then |y|; permute tables correspondingly
  void sort();

}; // |class common_block|

BlockElt twisted
  (const blocks::common_block& block, BlockElt z, const WeightInvolution& delta);

} // |namespace blocks|

} // |namespace atlas|
#endif
