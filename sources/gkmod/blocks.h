/*!
\file
\brief Class definition and function declarations for class Block.
*/
/*
  This is blocks.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef BLOCKS_H  /* guard against multiple inclusions */
#define BLOCKS_H

#include "blocks_fwd.h"
#include <cassert>
#include <iostream>

#include "rootdata_fwd.h"
#include "subdatum.h"
#include "bruhat_fwd.h"
#include "descents.h"
#include "kgb.h"
#include "complexredgp_fwd.h"
#include "realredgp_fwd.h"
#include "realform.h"
#include "weyl.h"
#include "tits_fwd.h"

#include "bitset.h"
#include "bitmap.h"
#include "hashtable.h"

namespace atlas {

/******** type declarations *************************************************/

/******** constant declarations *********************************************/

namespace blocks {

// reserve the last possible unsigned value; it is supposed unused
const BlockElt UndefBlock = ~ BlockElt(0);

}

/******** function declarations *********************************************/

namespace blocks {

  weyl::TwistedInvolution
    dual_involution(const weyl::TwistedInvolution& w,
		    const weyl::TwistedWeylGroup& W,
		    const weyl::TwistedWeylGroup& dual_W);

  // map from numbering of |b| to that of |dual_b|, assuming latter is dual
  std::vector<BlockElt>
    dual_map(const Block_base& b, const Block_base& dual_b);

  bitmap::BitMap common_Cartans(realredgp::RealReductiveGroup& GR,
				realredgp::RealReductiveGroup& dGR);

}

/******** type definitions **************************************************/

namespace blocks {

class Block_base {

  const weyl::WeylGroup& W; // used only for |compute_support| and printing

 protected: // other fields may be set in derived class contructor

  kgb::KGBEltList d_x;  // of size |size()|
  kgb::KGBEltList d_y; // of size |size()|

/*!\brief maps KGB element |x| to the first block element |z| with |d_x[z]>=x|.
*/
  std::vector<BlockElt> d_first_z_of_x; // of size |xrange|
  std::vector<BlockEltList> d_cross; // of size |d_rank| * |size()|
  std::vector<BlockEltPairList> d_cayley; // of size |d_rank| * |size()|
  descents::DescentStatusList d_descent; // of size |size()|
  std::vector<size_t> d_length; // of size |size()|

 public:

// constructors and destructors
  Block_base(const kgb::KGB& kgb,const kgb::KGB& dual_kgb);
  Block_base(const subdatum::SubSystem& sub,
	     const weyl::WeylGroup& printing_W);

  virtual ~Block_base() {}

// copy, assignment and swap

// accessors
  const weyl::WeylGroup& weylGroup() const { return W; }

  size_t rank() const { return W.rank(); } // only semisimple rank matters
  size_t size() const { return d_x.size(); }

  virtual size_t xsize() const = 0;
  virtual size_t ysize() const = 0;

  kgb::KGBElt x(BlockElt z) const { assert(z<size()); return d_x[z]; }
  kgb::KGBElt y(BlockElt z) const { assert(z<size()); return d_y[z]; }

  //!\brief Look up element by |x|, |y| coordinates
  BlockElt element(kgb::KGBElt x,kgb::KGBElt y) const;

  size_t length(BlockElt z) const { return d_length[z]; }

  BlockElt length_first(size_t l) const; // first element of given length

  virtual size_t Cartan_class(BlockElt z) const = 0;
  size_t max_Cartan() const // maximal Cartan number, for printing
  { return Cartan_class(size()-1); } // this should be OK in all cases

  virtual const weyl::TwistedInvolution& involution(BlockElt z) const = 0;

  BlockElt cross(size_t s, BlockElt z) const //!< cross action
  { assert(z<size()); assert(s<rank()); return d_cross[s][z]; }

  BlockEltPair cayley(size_t s, BlockElt z) const //!< Cayley transform
  { assert(z<size()); assert(s<rank());
    if (not isWeakDescent(s,z))
      return d_cayley[s][z];
    else return BlockEltPair(UndefBlock,UndefBlock);
  }

  BlockEltPair inverseCayley(size_t s, BlockElt z) const //!< inverse Cayley
  { assert(z<size()); assert(s<rank());
    if (isWeakDescent(s,z))
      return d_cayley[s][z];
    else return BlockEltPair(UndefBlock,UndefBlock);
  }

  const descents::DescentStatus& descent(BlockElt z) const
    { assert(z<size()); return d_descent[z]; }
  descents::DescentStatus::Value descentValue(size_t s, BlockElt z) const
    { assert(z<size()); assert(s<rank()); return d_descent[z][s]; }

  bool isWeakDescent(size_t s, BlockElt z) const
    { return descents::DescentStatus::isDescent(descentValue(s,z)); }

  bool isStrictAscent(size_t, BlockElt) const;
  bool isStrictDescent(size_t, BlockElt) const;
  size_t firstStrictDescent(BlockElt z) const;
  size_t firstStrictGoodDescent(BlockElt z) const;


  /*! \brief the functor \f$T_{\alpha,\beta}\f$ */
  BlockEltPair link(size_t alpha,size_t beta,BlockElt y) const;

  std::pair<BlockElt,BlockElt> R_packet(BlockElt z) const
  {
    assert(z<size());
    BlockElt x=d_x[z];
    return std::make_pair(d_first_z_of_x[x],d_first_z_of_x[x+1]);
  }

  // print derivative class specific per-element information
  virtual std::ostream& print(std::ostream& strm, BlockElt z) const
  { return strm; }

}; // |class Block_base|

  /*!
\brief Represents a block of representations of an inner form of G.

For our fixed inner form, orbits of $K$ on $G/B$ are parametrized by classes
of elements $x$ in $N_G(H).\delta$ (which is the normalizer in the second half
$G.\delta$ of the extended group $G^Gamma=G disju G.\delta$, where $\delta$ is
(i.e., acts on $G$ as) an involution that itself normalises $H$) modulo the
\emph{conjugation} action of $H$. (Dangerous bend: this $H$ conjugacy class of
$x$ is a subset, usually proper, of the coset $xH$. The collection of all $x$
is therefore NOT a subset of the extended Weyl group $N(H)/H$, but something
more subtle.) The requirement on $x$ is that it belong to to the $G$-conjugacy
class of strong involutions defining the inner form.

Each $x$ therefore defines an involution $theta_x$ of $H$.  Describing the
set of $x$ with a fixed involution is accomplished by the Fiber class.

A block is characterized by specifying also an inner form of the dual
group $G^vee$.  For this inner form, $K^vee$ orbits on $G^vee/B^vee$ are
parametrized by elements $y$.  The basic theorem is that the block of
representations is parametrized by pairs $(x,y)$ as above, subject to
the requirement that $theta_y$ is the negative transpose of $theta_x$.
  */
class Block : public Block_base
{

  enum State { BruhatConstructed, NumStates };

  const weyl::TwistedWeylGroup& tW;

  size_t xrange;
  size_t yrange;

  std::vector<size_t> d_Cartan; // of size |size()|
  weyl::TwistedInvolutionList d_involution; // of size |size()|

  /*!
\brief Flags the generators occurring in reduced expression for |d_involution|.
  */
  std::vector<bitset::RankFlags> d_involutionSupport; // of size |size()|

  /*!
\brief Records state bits (in fact one: whether the Bruhat order is computed).
  */
  bitset::BitSet<NumStates> d_state;

  /*!
\brief Bruhat order on the block.

Definition now corrected mathematically from the bad definition of
Vogan's Park City notes to one equivalent to the transitive closure of
non-vanishing KL polynomial.
  */
  bruhat::BruhatOrder* d_bruhat;

 public:

// constructors and destructors
  Block(const kgb::KGB& kgb,const kgb::KGB& dual_kgb);

  static Block build // pseudo contructor
    (complexredgp::ComplexReductiveGroup&,
     realform::RealForm rf, realform::RealForm drf, bool select_Cartans=false);

  ~Block(); // |delete d_bruhat;| but that does not compile here

// copy, assignment and swap
  Block(const Block& b); // copy contructor must handle |d_bruhat|
 private:
  Block& operator=(const Block& b); // we don't however need to assign
 public:

// accessors
  const weyl::TwistedWeylGroup& twistedWeylGroup() const { return tW; }

  virtual size_t xsize() const { return xrange; }
  virtual size_t ysize() const { return yrange; }

  size_t Cartan_class(BlockElt z) const
    { assert(z<size()); return d_Cartan[z]; }

/*!
  \brief Returns the twisted involution corresponding to z.

  This is the corresponding Weyl group element w, such that w.delta is the
  root datum involution tau corresponding to z
*/
  const weyl::TwistedInvolution& involution(BlockElt z) const
    { assert(z<size()); return d_involution[z]; }

  //! \brief the simple roots occurring in reduced expression |involution(z)|
  const bitset::RankFlags& involutionSupport(size_t z) const
  {
    assert(z<size());
    return d_involutionSupport[z];
  }

  // manipulators
  bruhat::BruhatOrder& bruhatOrder()
  {
    fillBruhat(); return *d_bruhat;
  }

  // private accessor and manipulators
private:
  void compute_supports(); // used during construction

  void fillBruhat();
}; // |class Block|


class gamma_block : public Block_base
{
  const kgb::KGB& kgb;

  latticetypes::RatWeight infin_char; // infinitesimal character

  std::vector<kgb::KGBElt> kgb_nr_of; // indexed by child |x| numbers

  struct y_fields
  {
    tits::GlobalTitsElement rep; //representative
    unsigned int Cartan_class;

    y_fields(tits::GlobalTitsElement y, unsigned int cc)
    : rep(y), Cartan_class(cc) {}
  }; // |struct y_fields|

  std::vector<y_fields> y_info; // indexed by child |y| numbers

 public:
  gamma_block(realredgp::RealReductiveGroup& GR,
	      const subdatum::SubSystem& sub,
	      kgb::KGBElt x,
	      const latticetypes::RatWeight& lambda, // discrete parameter
	      const latticetypes::RatWeight& gamma, // infinitesimal character
	      BlockElt& entry_element // set to block element matching the input
	      );

  // virtual methods
  size_t xsize() const { return kgb_nr_of.size(); }
  size_t ysize() const { return y_info.size(); }

  size_t Cartan_class(BlockElt z) const
  { assert(z<size()); return y_info[d_y[z]].Cartan_class; }

  const weyl::TwistedInvolution& involution(BlockElt z) const
  { assert(z<size()); return y_info[d_y[z]].rep.tw(); }

  std::ostream& print(std::ostream& strm, BlockElt z) const;

  // new methods
  latticetypes::RatWeight local_system(BlockElt z) const
  { assert(z<size()); return y_info[d_y[z]].rep.torus_part().as_rational(); }

}; // |class gamma_block|

class non_integral_block : public Block_base
{
  const kgb::KGB& kgb;

  latticetypes::RatWeight infin_char; // infinitesimal character

  std::vector<kgb::KGBElt> kgb_nr_of; // indexed by child |x| numbers

  struct y_fields
  {
    tits::GlobalTitsElement rep; //representative, in ^vG coordinates
    unsigned int Cartan_class; // for ^vG(gamma)

    y_fields(tits::GlobalTitsElement y, unsigned int cc)
    : rep(y), Cartan_class(cc) {}
  }; // |struct y_fields|

  std::vector<y_fields> y_info; // indexed by child |y| numbers

 public:
  non_integral_block
    (realredgp::RealReductiveGroup& GR,
     const subdatum::SubSystem& sub,
     kgb::KGBElt x,
     const latticetypes::RatWeight& lambda, // discrete parameter
     const latticetypes::RatWeight& gamma, // infinitesimal character
     BlockElt& entry_element // set to block element matching the input
    );

  // virtual methods
  size_t xsize() const { return kgb_nr_of.size(); }
  size_t ysize() const { return y_info.size(); }

  size_t Cartan_class(BlockElt z) const
  { assert(z<size()); return y_info[d_y[z]].Cartan_class; }

  const weyl::TwistedInvolution& involution(BlockElt z) const
  { assert(z<size()); return y_info[d_y[z]].rep.tw(); }

  std::ostream& print(std::ostream& strm, BlockElt z) const;

  // new methods
  latticetypes::RatWeight lambda(BlockElt z) const; // reconstruct from y value

}; // |class non_integral_block|

} // |namespace blocks|

} // |namespace atlas|
#endif
