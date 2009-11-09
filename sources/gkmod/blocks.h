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

#include "bruhat_fwd.h"
#include "descents.h"
#include "kgb_fwd.h"
#include "complexredgp_fwd.h"
#include "realredgp_fwd.h"
#include "weyl_fwd.h"

#include "bitset.h"
#include "bitmap.h"
#include "realform.h"
#include "weyl.h"

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
  std::vector<BlockElt>
    dual_map(const Block_base& b, const Block_base& dual_b);

  bitmap::BitMap common_Cartans(realredgp::RealReductiveGroup& GR,
				realredgp::RealReductiveGroup& dGR);

}

/******** type definitions **************************************************/

namespace blocks {

class Block_base {

  const weyl::TwistedWeylGroup& W;
  size_t xrange;
  size_t yrange;

  kgb::KGBEltList d_x;  // of size |size()|
  kgb::KGBEltList d_y; // of size |size()|

/*!\brief maps KGB element |x| to the first block element |z| with |d_x[z]>=x|.
*/
  std::vector<BlockElt> d_first_z_of_x; // of size |xrange|
  std::vector<BlockEltList> d_cross; // of size |d_rank| * |size()|
  std::vector<BlockEltPairList> d_cayley; // of size |d_rank| * |size()|
  descents::DescentStatusList d_descent; // of size |size()|
  std::vector<size_t> d_length; // of size |size()|
  std::vector<size_t> d_Cartan; // of size |size()|
  weyl::TwistedInvolutionList d_involution; // of size |size()|

 public:

// constructors and destructors
  Block_base(const kgb::KGB_base& kgb,const kgb::KGB_base& dual_kgb);

// copy, assignment and swap

// accessors
  const weyl::WeylGroup& weylGroup() const { return W.weylGroup(); }
  const weyl::TwistedWeylGroup& twistedWeylGroup() const { return W; }

  size_t rank() const { return W.rank(); }
  size_t size() const { return d_involution.size(); }

  size_t xsize() const { return xrange; }
  size_t ysize() const { return yrange; }

  kgb::KGBElt x(BlockElt z) const { assert(z<size()); return d_x[z]; }
  kgb::KGBElt y(BlockElt z) const { assert(z<size()); return d_y[z]; }

  //!\brief Look up element by |x|, |y| coordinates
  BlockElt element(kgb::KGBElt x,kgb::KGBElt y) const;

  size_t length(BlockElt z) const { return d_length[z]; }
  size_t Cartan_class(BlockElt z) const
    { assert(z<size()); return d_Cartan[z]; }

  const weyl::TwistedInvolution& involution(BlockElt z) const
    { assert(z<size()); return d_involution[z]; }

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

/*!
  \brief Returns the twisted involution corresponding to z.

This is the corresponding Weyl group element w, such that w.delta is the
root datum involution tau corresponding to z
*/


  std::pair<BlockElt,BlockElt> R_packet(BlockElt z) const
  {
    assert(z<size());
    BlockElt x=d_x[z];
    return std::make_pair(d_first_z_of_x[x],d_first_z_of_x[x+1]);
  }

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
  Block(complexredgp::ComplexReductiveGroup&, realform::RealForm rf,
	realform::RealForm drf, bool select_Cartans=false);

  Block(const kgb::KGB_base& kgb,const kgb::KGB_base& dual_kgb);

  ~Block();

// copy, assignment and swap

// accessors


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
  static Block_base // helper function for main contructor, same arguments
    base_for(complexredgp::ComplexReductiveGroup& G, realform::RealForm rf,
	     realform::RealForm drf, bool select_Cartans);

  void compute_supports(); // used during construction

  void fillBruhat();
}; // class Block

} // namespace blocks

} // namespace atlas
#endif
