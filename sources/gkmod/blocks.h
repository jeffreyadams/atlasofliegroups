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

#include "bruhat_fwd.h"
#include "descents.h"
#include "kgb_fwd.h"
#include "complexredgp_fwd.h"
#include "weyl_fwd.h"

#include "bitset.h"
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

/******** type definitions **************************************************/

namespace blocks {

  /*!
\brief Represents a block of representations of an inner form of G.

For our fixed inner form, orbits of K on G/B are parametrized by
certain elements x in N(H) (the normalizer in the extended group
G^Gamma), each defined up to conjugation by H. (Dangerous bend: the H
conjugacy class of x is a subset, usually proper, of the coset xH. The
collection of all x is therefore NOT a subset of the extended Weyl
group N(H)/H, but something more subtle.) The requirement on x is that
it belong to to the G-conjugacy class of strong involutions defining
the inner form.

Each x therefore defines an involution theta_x of H.  Describing the
set of x with a fixed involution is accomplished by the Fiber class.

A block is characterized by specifying also an inner form of the dual
group G^vee.  For this inner form, K^vee orbits on G^vee/B^vee are
parametrized by elements y.  The basic theorem is that the block of
representations is parametrized by pairs (x,y) as above, subject to
the requirement that theta_y is the negative transpose of theta_x.
  */
class Block {

 protected:

  enum State { BruhatConstructed, NumStates };

  /*!
\brief Semisimple rank of G.
  */
  size_t d_rank;

  /*!
\brief Number of K orbits on G/B.
  */
  size_t d_xsize;

  /*!
\brief Number of K^vee orbits on G^vee/B^vee.
  */
  size_t d_ysize;

  /*!
\brief Element d_x[z] identifies the K orbit x on G/B for z.
  */
  kgb::KGBEltList d_x;  // of size size()+1, with UndefKGB sentinel

  /*!
\brief Element d_y[z] identifies the K^vee orbit y on G^vee/B^vee for z.
  */
  kgb::KGBEltList d_y; // of size size()+1, with UndefKGB sentinel

  /*!
\brief d_cross[s][z] is $s \times z$ (for s a simple root, z a BlockElt).
  */
  std::vector<BlockEltList> d_cross; // of size d_rank

  /*!
\brief Element d_cayley[s][z] is the Cayley transform $c_s(z)$
  (for s a simple root, z noncompact imaginary) or undefined (otherwise).
  */
  std::vector<BlockEltPairList> d_cayley; // of size d_rank

  /*!
\brief Element d_inverseCayley[s]][z] is the inverse Cayley transform $c^s(z)$
 (for s a simple root, z is real type 1 or 2) or undefined (otherwise).
  */
  std::vector<BlockEltPairList> d_inverseCayley; // of size d_rank


  /*!
\brief Entry z flags the descent status of the simple roots for block
element z.
  */
  descents::DescentStatusList d_descent; // of size size()

  /*!
\brief Entry z is the length of block element z.
  */
  std::vector<size_t> d_length; // of size size()


  /*!
\brief Entry z (multiplied by the fixed outer automorphism delta) is
the involution $\theta_z$ of H attached to z
(in other words, d_involution[z] is the twisted involution attached to z; MvL)
  */
  weyl::WeylEltList d_involution; // of size size()

  /*!
\brief Entry z flags the simple roots occurring in $\theta_z$.
  */
  std::vector<bitset::RankFlags> d_involutionSupport; // of size size()

  /*!
\brief Number (in the list maintained by the complex reductive group)
of the real form of G where the block lives.
  */
  realform::RealForm d_realForm;

  /*!
\brief Number of the real form of G^vee defining the block.
  */
  realform::RealForm d_dualForm;

  /*!
\brief Records state bits (in fact one: whether the Bruhat order is computed).
  */
  bitset::BitSet<NumStates> d_state;

  /*!
\brief Bruhat order on the block.

Not used in the present code (DV 10/14/06).
  */
  bruhat::BruhatOrder* d_bruhat;

  const weyl::WeylGroup* d_weylGroup;

 public:

// constructors and destructors
  Block();

  Block(complexredgp::ComplexReductiveGroup&, realform::RealForm rf,
	realform::RealForm df);

  virtual ~Block();

// copy, assignment and swap
  void swap(Block&);

// accessors
  const bruhat::BruhatOrder& bruhatOrder() const {
    return *d_bruhat;
  }

  /*! \brief Cayley transform */
  BlockEltPair cayley(size_t s, BlockElt z) const {
    return d_cayley[s][z];
  }

  /*! \brief cross action */
  BlockElt cross(size_t s, BlockElt z) const {
    return d_cross[s][z];
  }

  /*! \brief vector of descent statuses of all simple roots */
  const descents::DescentStatus& descent(BlockElt z) const {
    return d_descent[z];
  }

  /*! \brief descent type of s at z */
  descents::DescentStatus::Value descentValue(size_t s, BlockElt z) const {
    return d_descent[z][s];
  }

  realform::RealForm dualForm() const {
    return d_dualForm;
  }

  size_t firstStrictDescent(BlockElt z) const;

  /*! \brief inverse Cayley transform */
  BlockEltPair inverseCayley(size_t s, BlockElt z) const {
    return d_inverseCayley[s][z];
  }

  /*! \brief the simple roots occurring in $\theta_z$. */
  const bitset::RankFlags& involutionSupport(size_t z) const {
    return d_involutionSupport[z];
  }

  bool isStrictAscent(size_t, BlockElt) const;

  bool isStrictDescent(size_t, BlockElt) const;

  /*! \brief length of block element */
  size_t length(BlockElt z) const {
    return d_length[z];
  }

  /*! \brief the functor $T_{\alpha,\beta}$ */
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
  const weyl::WeylElt& involution(BlockElt z) const {
    return d_involution[z];
}


  const weyl::WeylGroup& weylGroup() const {
    return *d_weylGroup;
  }

  kgb::KGBElt x(BlockElt z) const {
    return d_x[z];
  }

  size_t xsize() const {
    return d_xsize;
  }

  kgb::KGBElt y(BlockElt z) const {
    return d_y[z];
  }

  size_t ysize() const {
    return d_ysize;
  }

  // manipulators
  void fillBruhat();
};

}

}

#endif
