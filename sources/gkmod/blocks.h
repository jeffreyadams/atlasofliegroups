/*!
\file
\brief Class definition and function declarations for class Block.
*/
/*
  This is blocks.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef BLOCKS_H  /* guard against multiple inclusions */
#define BLOCKS_H

#include <cassert>
#include <iostream>

#include "ratvec.h"	// containment infinitesimal character

#include "atlas_types.h"
#include "tits.h"	// representative of $y$ in |non_integral_block|
#include "descents.h"	// inline methods


namespace atlas {

namespace blocks {

/******** type declarations *************************************************/

/******** constant declarations *********************************************/

// reserve the last possible unsigned value; it is supposed unused
const BlockElt UndefBlock = ~ BlockElt(0);



/******** function declarations *********************************************/

  // compute the involution in |dual_W| corresponding to |w| in |W|
  TwistedInvolution dual_involution
    (const TwistedInvolution& w,
     const TwistedWeylGroup& W,
     const TwistedWeylGroup& dual_W);

  // map from numbering of |b| to that of |dual_b|, assuming latter is dual
  std::vector<BlockElt> dual_map(const Block_base& b, const Block_base& dual_b);

  BitMap common_Cartans(RealReductiveGroup& GR,	RealReductiveGroup& dGR);



/******** type definitions **************************************************/



// The class |BlockBase| serves external functionality, not block construction
class Block_base
{
 public: // we need this typedef to be public, though used in derived classes
  struct EltInfo // per block element information
  {
    KGBElt x,y; // indices into |KGB| sets (which might no longer exist)
    DescentStatus descent;
    unsigned short length;
    BlockElt dual; // number of Hermitian dual of this element, if any
  EltInfo(KGBElt xx,KGBElt yy,DescentStatus dd, unsigned short ll)
  : x(xx),y(yy),descent(dd),length(ll), dual(UndefBlock) {}
  EltInfo(KGBElt xx,KGBElt yy, unsigned short ll)
  : x(xx),y(yy),descent(),length(ll), dual(UndefBlock) {}

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
  std::vector<std::vector<block_fields> > data;  // size |d_rank| * |size()|

  // map KGB element |x| to the first block element |z| with |this->x(z)>=x|
  // this vector may remain empty if |element| virtual methodis redefined
  std::vector<BlockElt> d_first_z_of_x; // of size |xsize+1|

  // possible tables of Bruhat order and Kazhdan-Lusztig polynomials
  BruhatOrder* d_bruhat;
  kl::KLContext* klc_ptr;

 public:

// constructors and destructors
  Block_base(const KGB& kgb,const KGB& dual_kgb);
  Block_base(unsigned int rank); // only dimensions some vectors

  virtual ~Block_base(); // deletes |d_bruhat| and |klc_ptr| (if non-NULL)

// copy, assignment and swap

  Block_base(const Block_base& b); // implemented but never used (optimized out)
 private:
  Block_base& operator=(const Block_base& b); // not implemented
 public:

// accessors

  size_t rank() const { return data.size(); } // semisimple rank matters
  size_t size() const { return info.size(); }

  virtual KGBElt xsize() const = 0;
  virtual KGBElt ysize() const = 0;

  KGBElt x(BlockElt z) const { assert(z<size()); return info[z].x; }
  KGBElt y(BlockElt z) const { assert(z<size()); return info[z].y; }

  // Look up element by |x|, |y| coordinates
  virtual BlockElt element(KGBElt x,KGBElt y) const;

  size_t length(BlockElt z) const { return info[z].length; }

  BlockElt length_first(size_t l) const; // first element of given length

  BlockElt cross(size_t s, BlockElt z) const //!< cross action
  { assert(z<size()); assert(s<rank()); return data[s][z].cross_image; }

  BlockEltPair cayley(size_t s, BlockElt z) const //!< Cayley transform
  { assert(z<size()); assert(s<rank());
    if (not isWeakDescent(s,z))
      return data[s][z].Cayley_image;
    else return BlockEltPair(UndefBlock,UndefBlock);
  }

  BlockEltPair inverseCayley(size_t s, BlockElt z) const //!< inverse Cayley
  { assert(z<size()); assert(s<rank());
    if (isWeakDescent(s,z))
      return data[s][z].Cayley_image;
    else return BlockEltPair(UndefBlock,UndefBlock);
  }

  const DescentStatus& descent(BlockElt z) const
    { assert(z<size()); return info[z].descent; }
  DescentStatus::Value descentValue(size_t s, BlockElt z) const
    { assert(z<size()); assert(s<rank()); return descent(z)[s]; }

  bool isWeakDescent(size_t s, BlockElt z) const
    { return DescentStatus::isDescent(descentValue(s,z)); }

  bool isStrictAscent(size_t, BlockElt) const;
  bool isStrictDescent(size_t, BlockElt) const;
  size_t firstStrictDescent(BlockElt z) const;
  size_t firstStrictGoodDescent(BlockElt z) const;

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

 protected:
  // a method to straighten out blocks generated in some non standard order
  // renumber |x| through |new_x|, then order increasingly, set |first_z_of_x|
  KGBElt renumber_x(const std::vector<KGBElt>& new_x);
  void compute_first_zs(); // set |first_z_of_x| according to |x| values

 private:
  void fillBruhat();
  void fill_klc(BlockElt last_y,bool verbose);

}; // |class Block_base|

  /*!
\brief Represents a block of representations of an inner form of G.

For our fixed inner form, orbits of $K$ on $G/B$ are parametrized by classes
of elements $x$ in $N_G(H).\delta$ (the normalizer in the non-identity
component $G.\delta$ of the extended group $G^Gamma=G disju G.\delta$, where
$\delta$ is (i.e., acts on $G$ as) an involution that itself normalises $H$),
modulo the \emph{conjugation} action of $H$. (Dangerous bend: this $H$
conjugacy class of $x$ is a subset, usually proper, of the coset $xH$. The
collection of all $x$ is therefore NOT a subset of the extended Weyl group
$N(H)/H$, but something more subtle.) The requirement on $x$ is that it belong
to the $G$-conjugacy class of strong involutions defining the inner form.

Each $x$ therefore defines an involution $\theta_x$ of $H$. Data pertaining to
the subset of $x$ with a fixed $\theta_x$ is stored in the |Fiber| class.

A block is characterized by specifying also an inner form of the dual
group $G^vee$. For this inner form, $K^vee$ orbits on $G^vee/B^vee$ are
parametrized by elements $y$. The basic theorem is that the block of
representations is parametrized by pairs $(x,y)$ as above, subject to
the requirement that $theta_y$ is the negative transpose of $theta_x$.
  */
class Block : public Block_base
{
  const TwistedWeylGroup& tW; // reference is used here only for printing

  size_t xrange;
  size_t yrange;

  // wasteful fields, but we cannot lean on |KGB|, which might no longer exist
  std::vector<size_t> d_Cartan; // of size |size()|
  TwistedInvolutionList d_involution; // of size |size()|

  /*!
\brief Flags the generators occurring in reduced expression for |d_involution|.
  */
  std::vector<RankFlags> d_involutionSupport; // of size |size()|


 public:

// constructors and destructors
  Block(const KGB& kgb,const KGB& dual_kgb);

  static Block build // pseudo contructor with small (and forgotten) KGB sets
    (ComplexReductiveGroup&, RealFormNbr rf, RealFormNbr drf);

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

  virtual KGBElt xsize() const { return xrange; }
  virtual KGBElt ysize() const { return yrange; }

  size_t Cartan_class(BlockElt z) const
    { assert(z<size()); return d_Cartan[z]; }

  size_t max_Cartan() const { return Cartan_class(size()-1); } // for printing

/*!
  \brief Returns the twisted involution corresponding to z.

  This is the corresponding Weyl group element w, such that w.delta is the
  root datum involution tau corresponding to z
*/
  const TwistedInvolution& involution(BlockElt z) const
    { assert(z<size()); return d_involution[z]; }

  //! \brief the simple roots occurring in reduced expression |involution(z)|
  const RankFlags& involutionSupport(size_t z) const
  {
    assert(z<size());
    return d_involutionSupport[z];
  }

  virtual std::ostream& print // defined in block_io.cpp
   (std::ostream& strm, BlockElt z,bool as_invol_expr) const;


  // private accessor and manipulators
private:
  void compute_supports(); // used during construction

}; // |class Block|

typedef HashTable<y_entry,KGBElt> y_part_hash;

// |param_block| is intermediate between |Block_base| and |non_integral_block|
class param_block : public Block_base // blocks of parameters
{
 protected: // everything that is here serves derived classes only
  const KGB& kgb; // on the |x| size we employ a pre-computed KGB structure

  RatWeight infin_char; // infinitesimal character

  std::vector<KGBElt> kgb_nr_of; // maps child |x| numbers to parent |kgb|
  std::vector<KGBElt> x_of;      // inverse mapping, partial

  y_entry::Pooltype y_pool;
  y_part_hash y_hash;

  param_block(const KGB& parent_kgb, unsigned int rank);

  // auxiliary for construction
  void compute_duals(const ComplexReductiveGroup& G);

 public:
  const RatWeight& gamma() const { return infin_char; }
  KGBElt parent_x(BlockElt z) const { return kgb_nr_of[x(z)]; }
  const TorusElement& y_rep(KGBElt y) const { return y_pool[y].repr(); }
  const TwistedInvolution& involution(BlockElt z) const; // obtained from |kgb|

  // virtual methods
  virtual KGBElt xsize() const { return kgb_nr_of.size(); } // child |x| range
  virtual KGBElt ysize() const { return y_hash.size(); }    // child |y| range

}; // |class param_block|

/*
  The class |gamma_block| class was meant to compute blocks at non-integral
  infinitesimal character |gamma| using only the integrality subsystem on the
  dual side. However base-point shifting when converting to and from |y|
  values using this representation is currently not well enough understood,
  and therefore this class does not function reliably. The alternative class
  |nonintegeral_block| defined below avoids this problem and appears correct.
*/
class gamma_block : public param_block
{
  std::vector<TorusElement> y_rep; // representatives, by child |y| numbers

 public:
  gamma_block(RealReductiveGroup& GR,
	      const SubSystemWithGroup& sub,
	      KGBElt x,			// starting |x| value
	      const RatWeight& lambda,	// discrete parameter
	      const RatWeight& gamma,	// infinitesimal character
	      BlockElt& entry_element	// set to block element matching input
	      );


  virtual std::ostream& print // in block_io.cpp
   (std::ostream& strm, BlockElt z,bool as_invol_expr) const;

  // new methods
  RatWeight local_system(BlockElt z) const // reconstruct a |lambda| from |y|
  { assert(z<size()); return y_rep[y(z)].log_2pi(); }

}; // |class gamma_block|


typedef Block_base::EltInfo block_elt_entry;
typedef HashTable<block_elt_entry,BlockElt> block_hash;

class non_integral_block : public param_block
{

  // A simple structure to pack a pair of already sequenced numbers (indices
  // into the |info| field for some (future) block into a hashable value


  RealReductiveGroup& GR; // non-const to allow construction of |Rep_context|
  RankFlags singular; // flags simple roots for which |infin_char| is singular
  block_hash z_hash; //  on |Block_base::info|

 public:
  non_integral_block // rewritten constructor, for full block
    (const repr::Rep_context& rc,
     StandardRepr sr, // by value,since it will be made dominant before use
     BlockElt& entry_element	// set to block element matching input
    );

  non_integral_block // alternative constructor, for interval below |sr|
    (const repr::Rep_context& rc,
     StandardRepr sr); // by value,since it will be made dominant before use

  // "inherited" accessors
  const ComplexReductiveGroup& complexGroup() const;
  const InvolutionTable& involution_table() const;

  virtual BlockElt element(KGBElt x,KGBElt y) const; // redefined using |z_hash|
  virtual std::ostream& print // defined in block_io.cpp
    (std::ostream& strm, BlockElt z,bool as_invol_expr) const;

  // new methods

  RealReductiveGroup& realGroup() const { return GR; }

  RatWeight nu(BlockElt z) const; // "real" projection of |infin_char|
  RatWeight y_part(BlockElt z) const; // raw torus part info, normalized
  Weight lambda_rho(BlockElt z) const; // reconstruct from y value
  RatWeight lambda(BlockElt z) const; // reconstruct from y value
  RankFlags singular_simple_roots() { return singular; }
  bool survives(BlockElt z) const; // whether $J(z_{reg})$ survives tr. functor
  BlockEltList survivors_below(BlockElt z) const; // expression for $I(z)$

 private:
  void add_z(KGBElt x,KGBElt y, unsigned short l);

}; // |class non_integral_block|

} // |namespace blocks|

} // |namespace atlas|
#endif
