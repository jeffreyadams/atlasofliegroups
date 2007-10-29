/*!
\file
\brief Class definition and function declarations for class Block.
*/
/*
  This is blocks.h

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
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

  std::vector<BlockElt> dual_map(const Block& b, const Block& dual_b);
  bitmap::BitMap common_Cartans(realredgp::RealReductiveGroup& GR,
				realredgp::RealReductiveGroup& dGR);

}

/******** type definitions **************************************************/

namespace blocks {

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
class Block {

  enum State { BruhatConstructed, NumStates };

  /*!
\brief Number (in the list maintained by the complex reductive group)
of the real form of G where the block lives.
  */
  realform::RealForm d_realForm;

  /*!
\brief Number of the real form of G^vee defining the block.
  */
  realform::RealForm d_dualForm;

  /*! \brief Semisimple rank of G. */
  size_t d_rank;

  const weyl::WeylGroup& d_weylGroup;

  /*!
\brief Number of K orbits on G/B.
  */
  size_t d_xrange;

  /*!
\brief Number of K^vee orbits on G^vee/B^vee.
  */
  size_t d_yrange;

  /*!
\brief Element d_x[z] identifies the K orbit x on G/B for z.
  */
  kgb::KGBEltList d_x;  // of size |size()|

  /*!
\brief Element d_y[z] identifies the K^vee orbit y on G^vee/B^vee for z.
  */
  kgb::KGBEltList d_y; // of size |size()|

/*!\brief maps KGB element |x| to the first block element |z| with |d_x[z]>=x|.
*/
  std::vector<BlockElt> d_first_z_of_x; // of size |d_xrange|

  /*!
\brief d_cross[s][z] is $s * z$ (for s a simple root, z a BlockElt).
  */
  std::vector<BlockEltList> d_cross; // of size |d_rank| * |size()|

  /*! \brief For $s$ a simple root and $z$ a |BlockElt|, |d_cayley[s][z]| is
the Cayley transform $c_s(z)$ (noncompact imaginary case) or the inverse
Cayley transform $c^s(z)$ (real type 1 or 2) or undefined (otherwise).
  */
  std::vector<BlockEltPairList> d_cayley; // of size |d_rank| * |size()|


  /*!
\brief Entry z flags the descent status of the simple roots for block
element z.
  */
  descents::DescentStatusList d_descent; // of size |size()|

  /*!
\brief Entry z is the length of block element z.
  */
  std::vector<size_t> d_length; // of size |size()|

  /*!
\brief Entry z is the Cartan class of block element z.
  */
  std::vector<size_t> d_Cartan; // of size |size()|


  /*!
\brief Entry z (multiplied by the fixed outer automorphism \delta) is
the involution \f$\tau_z\f$ of H attached to z
(in other words, d_involution[z] is the twisted involution attached to z; MvL)
  */
  weyl::TwistedInvolutionList d_involution; // of size |size()|

  /*!
\brief Entry z flags the simple roots occurring in \f$\theta_z\f$.
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

  ~Block();

// copy, assignment and swap

// accessors
  /*! \brief vector of descent statuses of all simple roots */
  const descents::DescentStatus& descent(BlockElt z) const {
    assert(z<size());
    return d_descent[z];
  }

  /*! \brief descent type of s at z */
  descents::DescentStatus::Value descentValue(size_t s, BlockElt z) const {
    assert(z<size());
    assert(s<d_rank);
    return d_descent[z][s];
  }

  //!\brief whether |s| is a weak descent for |z|
  bool isWeakDescent(size_t s, BlockElt z) const
    { return descents::DescentStatus::isDescent(descentValue(s,z)); }

  realform::RealForm dualForm() const {
    return d_dualForm;
  }

  size_t firstStrictDescent(BlockElt z) const;
  size_t firstStrictGoodDescent(BlockElt z) const;

  /*! \brief cross action */
  BlockElt cross(size_t s, BlockElt z) const {
    assert(z<size());
    assert(s<d_rank);
    return d_cross[s][z];
  }

  /*! \brief Cayley transform */
  BlockEltPair cayley(size_t s, BlockElt z) const {
    assert(z<size());
    assert(s<d_rank);
    if (not isWeakDescent(s,z))
      return d_cayley[s][z];
    else return BlockEltPair(UndefBlock,UndefBlock);
  }

  /*! \brief inverse Cayley transform */
  BlockEltPair inverseCayley(size_t s, BlockElt z) const {
    assert(z<size());
    assert(s<d_rank);
    if (isWeakDescent(s,z))
      return d_cayley[s][z];
    else return BlockEltPair(UndefBlock,UndefBlock);
  }

  /*! \brief the simple roots occurring in \f$\theta_z\f$. */
  const bitset::RankFlags& involutionSupport(size_t z) const {
    assert(z<size());
    return d_involutionSupport[z];
  }

  bool isStrictAscent(size_t, BlockElt) const;

  bool isStrictDescent(size_t, BlockElt) const;

  /*! \brief length of block element */
  size_t length(BlockElt z) const {
    return d_length[z];
  }

  /*! \brief Cartan class of block element */
  size_t Cartan_class(BlockElt z) const {
    assert(z<size());
    return d_Cartan[z];
  }

  /*! \brief the functor \f$T_{\alpha,\beta}\f$ */
  BlockEltPair link(size_t alpha,size_t beta,BlockElt y) const;

  /*! \brief semisimple rank of the group this block is constructed for */
  size_t rank() const {
    return d_rank;
  }

  realform::RealForm realForm() const {
    return d_realForm;
  }

  /*! \brief size of the block */
  size_t size() const {
    return d_involution.size(); // d_involution is one of many possible vectors
  }

/*!
  \brief Returns the twisted involution corresponding to z.

This is the corresponding Weyl group element w, such that w.delta is the
root datum involution tau corresponding to z
*/
  const weyl::TwistedInvolution& involution(BlockElt z) const{
    assert(z<size());
    return d_involution[z];
}


  const weyl::WeylGroup& weylGroup() const {
    return d_weylGroup;
  }

  kgb::KGBElt x(BlockElt z) const {
    assert(z<size());
    return d_x[z];
  }

  kgb::KGBElt y(BlockElt z) const {
    assert(z<size());
    return d_y[z];
  }

  //!\brief Look up element by |x|, |y| coordinates
  BlockElt element(kgb::KGBElt x,kgb::KGBElt y) const;

  std::pair<BlockElt,BlockElt> R_packet(BlockElt z) const
  {
    assert(z<size());
    BlockElt x=d_x[z];
    return std::make_pair(d_first_z_of_x[x],d_first_z_of_x[x+1]);
  }

  size_t xsize() const {
    return d_xrange;
  }
  size_t ysize() const {
    return d_yrange;
  }

  // manipulators
  bruhat::BruhatOrder& bruhatOrder()
  {
    fillBruhat(); return *d_bruhat;
  }

  // private accessor and manipulators
private:
  weyl::TwistedInvolution dualInvolution
    (const weyl::TwistedInvolution& tw,weyl::WeylInterface to_dual) const;

  void generate(realredgp::RealReductiveGroup& G,
		realredgp::RealReductiveGroup& dG,
		bool select_Cartans);

  void fillBruhat();
}; // class Block

} // namespace blocks

} // namespace atlas
#endif
