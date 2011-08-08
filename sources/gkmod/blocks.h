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

#include <cassert>
#include <iostream>

#include "ratvec.h"	// containment infinitesimal character

#include "atlas_types.h"
#include "tits.h"	// representative of $y$ in |non_integeral_block|
#include "descents.h"	// inline methods


namespace atlas {

/******** type declarations *************************************************/

/******** constant declarations *********************************************/

namespace blocks {

// reserve the last possible unsigned value; it is supposed unused
const BlockElt UndefBlock = ~ BlockElt(0);

}

/******** function declarations *********************************************/

namespace blocks {

  TwistedInvolution
    dual_involution(const TwistedInvolution& w,
		    const TwistedWeylGroup& W,
		    const TwistedWeylGroup& dual_W);

  // map from numbering of |b| to that of |dual_b|, assuming latter is dual
  std::vector<BlockElt>
    dual_map(const Block_base& b, const Block_base& dual_b);

  BitMap common_Cartans(RealReductiveGroup& GR,
				RealReductiveGroup& dGR);

}

/******** type definitions **************************************************/

namespace blocks {

class Block_base {

  const WeylGroup& W; // used only for |compute_support| and printing

 protected: // other fields may be set in derived class contructor

  KGBEltList d_x;  // of size |size()|
  KGBEltList d_y; // of size |size()|

/*!\brief maps KGB element |x| to the first block element |z| with |d_x[z]>=x|.
*/
  std::vector<BlockElt> d_first_z_of_x; // of size |xsize+1|
  std::vector<BlockEltList> d_cross; // of size |d_rank| * |size()|
  std::vector<BlockEltPairList> d_cayley; // of size |d_rank| * |size()|
  DescentStatusList d_descent; // of size |size()|
  std::vector<size_t> d_length; // of size |size()|

 public:

// constructors and destructors
  Block_base(const KGB& kgb,const KGB& dual_kgb);
  Block_base(const SubSystem& sub,
	     const WeylGroup& printing_W);

  virtual ~Block_base() {}

// copy, assignment and swap

// accessors
  const WeylGroup& weylGroup() const { return W; }

  size_t rank() const { return d_cross.size(); } // semisimple rank matters
  size_t size() const { return d_x.size(); }

  virtual size_t xsize() const = 0;
  virtual size_t ysize() const = 0;

  KGBElt x(BlockElt z) const { assert(z<size()); return d_x[z]; }
  KGBElt y(BlockElt z) const { assert(z<size()); return d_y[z]; }

  //!\brief Look up element by |x|, |y| coordinates
  BlockElt element(KGBElt x,KGBElt y) const;

  size_t length(BlockElt z) const { return d_length[z]; }

  BlockElt length_first(size_t l) const; // first element of given length

  virtual const TwistedInvolution& involution(BlockElt z) const = 0;

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

  const DescentStatus& descent(BlockElt z) const
    { assert(z<size()); return d_descent[z]; }
  DescentStatus::Value descentValue(size_t s, BlockElt z) const
    { assert(z<size()); assert(s<rank()); return d_descent[z][s]; }

  bool isWeakDescent(size_t s, BlockElt z) const
    { return DescentStatus::isDescent(descentValue(s,z)); }

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

  // a method to straighten out blocks generated in some non standard order
  // renumber |x| through |new_x|, order by increasing |x|, set |first_z_of_x|
 protected:
  KGBElt renumber_x(const std::vector<KGBElt>& new_x);

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

  const TwistedWeylGroup& tW; // for interpreting twisted involutions

  size_t xrange;
  size_t yrange;

  std::vector<size_t> d_Cartan; // of size |size()|
  TwistedInvolutionList d_involution; // of size |size()|

  /*!
\brief Flags the generators occurring in reduced expression for |d_involution|.
  */
  std::vector<RankFlags> d_involutionSupport; // of size |size()|

  /*!
\brief Records state bits (in fact one: whether the Bruhat order is computed).
  */
  BitSet<NumStates> d_state;

  /*!
\brief Bruhat order on the block.

Definition now corrected mathematically from the bad definition of
Vogan's Park City notes to one equivalent to the transitive closure of
non-vanishing KL polynomial.
  */
  BruhatOrder* d_bruhat;

 public:

// constructors and destructors
  Block(const KGB& kgb,const KGB& dual_kgb);

  static Block build // pseudo contructor with small (and forgotten) KGB sets
    (ComplexReductiveGroup&, RealFormNbr rf, RealFormNbr drf);

  static Block build // pseudo contructor with stored KGB sets
    (RealReductiveGroup& G_R, RealReductiveGroup& dG_R);

  ~Block(); // |delete d_bruhat;| but that does not compile here

// copy, assignment and swap
  Block(const Block& b); // copy contructor must handle |d_bruhat|
 private:
  Block& operator=(const Block& b); // we don't however need to assign
 public:

// accessors
  const TwistedWeylGroup& twistedWeylGroup() const { return tW; }

  virtual size_t xsize() const { return xrange; }
  virtual size_t ysize() const { return yrange; }

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

  // manipulators
  BruhatOrder& bruhatOrder()
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
  const KGB& kgb;

  RatWeight infin_char; // infinitesimal character

  std::vector<KGBElt> kgb_nr_of; // indexed by child |x| numbers

  struct y_fields
  {
    GlobalTitsElement rep; //representative
    unsigned int Cartan_class;

    y_fields(GlobalTitsElement y, unsigned int cc)
    : rep(y), Cartan_class(cc) {}
  }; // |struct y_fields|

  std::vector<y_fields> y_info; // indexed by child |y| numbers

 public:
  gamma_block(RealReductiveGroup& GR,
	      const SubSystem& sub,
	      KGBElt x,
	      const RatWeight& lambda, // discrete parameter
	      const RatWeight& gamma, // infinitesimal character
	      BlockElt& entry_element // set to block element matching the input
	      );

  // virtual methods
  size_t xsize() const { return kgb_nr_of.size(); }
  size_t ysize() const { return y_info.size(); }

  size_t Cartan_class(BlockElt z) const
  { assert(z<size()); return y_info[d_y[z]].Cartan_class; }

  size_t max_Cartan() const // maximal Cartan number, for printing
  { return Cartan_class(size()-1); } // this should be OK in all cases

  const TwistedInvolution& involution(BlockElt z) const
  { assert(z<size()); return y_info[d_y[z]].rep.tw(); }

  std::ostream& print(std::ostream& strm, BlockElt z) const;

  // new methods
  RatWeight local_system(BlockElt z) const
  { assert(z<size()); return y_info[d_y[z]].rep.torus_part().log_2pi(); }

}; // |class gamma_block|

class non_integral_block : public Block_base
{
  const KGB& kgb;
  const ComplexReductiveGroup& G;
  const SubSystem& sub;

  RankFlags singular;

  const RatWeight infin_char; // infinitesimal character

  std::vector<KGBElt> kgb_nr_of; // indexed by child |x| numbers
  std::vector<GlobalTitsElement> y_info; // indexed by child |y| numbers

 public:
  non_integral_block
    (RealReductiveGroup& GR,
     const SubSystem& subsys,
     KGBElt x,
     const RatWeight& lambda, // discrete parameter
     const RatWeight& gamma, // infinitesimal character
     BlockElt& entry_element // set to block element matching the input
    );

  non_integral_block // alternative constructor, for interval below |x|
    (RealReductiveGroup& GR,
     const SubSystem& subsys,
     KGBElt x,
     const RatWeight& lambda, // discrete parameter
     const RatWeight& gamma // infinitesimal character
    );

  // virtual methods
  size_t xsize() const { return kgb_nr_of.size(); }
  size_t ysize() const { return y_info.size(); }

  const TwistedInvolution& involution(BlockElt z) const
  { assert(z<size()); return y_info[d_y[z]].tw(); }

  std::ostream& print(std::ostream& strm, BlockElt z) const;

  // new methods
  RatWeight lambda(BlockElt z) const; // reconstruct from y value
  RankFlags singular_simple_roots() { return singular; }
  bool is_nonzero(BlockElt z) const; // whether |z| survives singular |gamma|
  BlockEltList nonzeros_below(BlockElt z) const; // reachable nonzeros

}; // |class non_integral_block|

} // |namespace blocks|

} // |namespace atlas|
#endif
