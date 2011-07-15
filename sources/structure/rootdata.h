/*!
\file
\brief Class definitions and function declarations for the RootDatum class.
*/
/*
  This is rootdata.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006--2010 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#ifndef ROOTDATA_H  /* guard against multiple inclusions */
#define ROOTDATA_H

#include "atlas_types.h"

#include <algorithm>

#include "arithmetic_fwd.h"
#include "atlas_types.h"
#include "tags.h"

#include "lietype_fwd.h"
#include "atlas_types.h"

#include "bitset.h" // for ascent and descent sets
#include "matrix.h" // for loads of subobjects
#include "permutations.h" // for storing root permutations

namespace atlas {

/******** type declarations *************************************************/


/******** function declarations *********************************************/

namespace rootdata {

CoweightInvolution dualBasedInvolution
  (const WeightInvolution&, const RootDatum&);

RootNbrSet makeOrthogonal(const RootNbrSet& o, const RootNbrSet& subsys,
		       const RootSystem& rs);

void toDistinguished(WeightInvolution&, const RootDatum&);

// make |Delta| simple (a twist) by Weyl group action; return Weyl word
// whose left-multiplication transforms returned |Delta| into original one
WeylWord wrt_distinguished(const RootSystem& rs, RootNbrList& Delta);

WeightInvolution refl_prod(const RootNbrSet&, const RootDatum&);

RootDatum integrality_datum(const RootDatum& rd,
			    const RatWeight& lambda);

arithmetic::RationalList integrality_points(const RootDatum& rd,
					    RatWeight& nu);

} // namespace rootdata

/******** type definitions **************************************************/

namespace rootdata {


class RootSystem
{
  typedef signed char byte;
  typedef matrix::Vector<byte> Byte_vector;
  struct root_info
  {
    Byte_vector root, dual; // root in root basis, coroot in coroot basis
    bitset::RankFlags descents,ascents; // for reflections by simple roots

    root_info(const Byte_vector& v)
    : root(v), dual(), descents(), ascents() {}
  };
  struct root_compare; // necessary for sorting roots

  size_t rk; // rank of root system

  Byte_vector Cmat; // Cartan matrix in compressed format

  std::vector<root_info> ri; //!< List of information about positive roots

  matrix::Vector<int> two_rho_in_simple_roots;

//!\brief Root permutations induced by reflections in all positive roots.
  std::vector<permutations::Permutation> root_perm;

  // internal access methods
  byte& Cartan_entry(weyl::Generator i, weyl::Generator j)
    { return Cmat[i*rk+j]; }
  const byte& Cartan_entry(weyl::Generator i, weyl::Generator j) const
    { return Cmat[i*rk+j]; }
  Byte_vector& root(RootNbr i) { return ri[i].root;}
  Byte_vector& coroot(RootNbr i) { return ri[i].dual;}
  const Byte_vector& root(RootNbr i) const { return ri[i].root;}
  const Byte_vector& coroot(RootNbr i) const { return ri[i].dual;}
  RootNbr rt_abs(RootNbr alpha) const // offset of corresponding positive root
    { return isPosRoot(alpha) ? alpha-numPosRoots() : numPosRoots()-1-alpha; }

  void cons(const int_Matrix& Cartan_matrix); // panse bete
 public:

// constructors and destructors

  explicit RootSystem(const int_Matrix& Cartan_matrix);

  RootSystem(const RootSystem& rs, tags::DualTag);

// accessors

  size_t rank() const { return rk; } // semisimple rank when part of |RootDatum|
  unsigned long numPosRoots() const { return ri.size(); }
  unsigned long numRoots() const { return 2*numPosRoots(); }

  // Cartan matrix by entry and as a whole
  int cartan(weyl::Generator i, weyl::Generator j) const
  { return Cartan_entry(i,j); };
  int_Matrix cartanMatrix() const;
  lietype::LieType Lie_type() const;

  // for subsystem
  int_Matrix cartanMatrix(const RootNbrList& sub) const;
  lietype::LieType Lie_type(RootNbrList sub) const;



// root list access

  // express root in simple root basis, or coroot in simple coroot basis
  int_Vector root_expr(RootNbr alpha) const;
  int_Vector coroot_expr(RootNbr alpha) const;

  // convert sequence of root numbers to expressions in the simple roots
  template <typename I, typename O>
    void toRootBasis(I, I, O) const;
  // convert sequence of root numbers to expressions in subsystem simple roots
  template <typename I, typename O>
    void toRootBasis(I, I, O, const RootNbrList&) const;
  // convert sequence of root numbers to expressions in subsystem simple weights
  template <typename I, typename O>
    void toSimpleWeights(I, I, O, const RootNbrList&) const;

  bool isSimpleRoot(RootNbr alpha) const
  { return alpha-numPosRoots()<rk; } // this uses that |RootNbr| is unsigned

  bool isPosRoot(RootNbr alpha) const
  { return alpha>=numPosRoots(); } // second half

  RootNbr simpleRootNbr(weyl::Generator i) const
  { assert(i<rk);  return numPosRoots()+i; }

  RootNbr posRootNbr(size_t alpha) const
  { assert(alpha<numPosRoots()); return numPosRoots()+alpha; }

  RootNbr simpleRootIndex(size_t alpha) const
  { assert(isSimpleRoot(alpha));  return alpha-numPosRoots(); }

  RootNbr posRootIndex(size_t alpha) const
    { assert(isPosRoot(alpha)); return alpha-numPosRoots(); }

  RootNbr rootMinus(RootNbr alpha) const // roots are ordered symmetrically
  { return numRoots()-1-alpha; }


  RootNbrSet simpleRootSet() const; // NOT for iteration over it; never used
  RootNbrList simpleRootList() const; // NOT for iteration over it
  RootNbrSet posRootSet() const; // NOT for iteration, may serve as mask

// other accessors

  // the next method only works for _simple_ roots! (whence no RootNbr for |i|)
  const permutations::Permutation& simple_root_permutation(weyl::Generator i)
    const
  { return root_perm[i]; }

  bitset::RankFlags descent_set(RootNbr alpha) const
  {
    RootNbr a = rt_abs(alpha);
    return isPosRoot(alpha) ? ri[a].descents : ri[a].ascents;
  }
  bitset::RankFlags ascent_set(RootNbr alpha) const
  {
    RootNbr a = rt_abs(alpha);
    return isPosRoot(alpha) ? ri[a].ascents : ri[a].descents;
  }

  size_t find_descent(RootNbr alpha) const
  { return descent_set(alpha).firstBit(); }

  bool is_descent(weyl::Generator i, RootNbr alpha) const
  { return root_perm[i][alpha]<alpha; } // easier than using |descent_set|

  bool is_ascent(weyl::Generator i, RootNbr alpha) const
  { return root_perm[i][alpha]>alpha; }

  RootNbr simple_reflected_root(RootNbr r, weyl::Generator i) const
  { return simple_root_permutation(i)[r]; }

  void simple_reflect_root(RootNbr& r, weyl::Generator i) const
  { r=simple_root_permutation(i)[r]; }

  RootNbr permuted_root(const WeylWord& ww, RootNbr r) const
  {
    for (weyl::Generator i=ww.size(); i-->0; )
      simple_reflect_root(r,ww[i]);
    return r;
  }

  RootNbr permuted_root(RootNbr r,const WeylWord& ww) const
  {
    for (weyl::Generator i=0; i<ww.size(); ++i)
      simple_reflect_root(r,ww[i]);
    return r;
  }

  // for arbitrary roots, reduce root number to positive root offset first
  const permutations::Permutation& root_permutation(RootNbr alpha) const
  { return root_perm[rt_abs(alpha)]; }

  bool isOrthogonal(RootNbr alpha, RootNbr beta) const
  { return root_permutation(alpha)[beta]==beta; }

  // pairing between root |alpha| and coroot |beta|
  int bracket(RootNbr alpha, RootNbr beta) const;



  // find permutation of roots induced by diagram automorphism
  permutations::Permutation root_permutation(const permutations::Permutation& ) const;

  // extend root datum automorphism given on simple roots to all roots
  permutations::Permutation extend_to_roots(const RootNbrList&) const;


  WeylWord reflectionWord(RootNbr r) const;

  // express sum of roots positive for |Delta| in fundamental weights
  matrix::Vector<int> pos_system_vec(const RootNbrList& Delta) const;

  RootNbrList simpleBasis(RootNbrSet rs) const; // find simple basis for subsystem

  bool sumIsRoot(RootNbr alpha, RootNbr beta, RootNbr& gamma) const;
  bool sumIsRoot(RootNbr alpha, RootNbr beta) const
  { RootNbr dummy; return sumIsRoot(alpha,beta,dummy); }

  RootNbrSet long_orthogonalize(const RootNbrSet& rest) const;

  RootNbrList high_roots() const;

// manipulators


}; // |class RootSystem|

/*!
\brief Based root datum for a complex reductive group.

What we call a root datum in this program is what is usually called a based
root datum, in other words a fixed choice of positive roots is always assumed.

The root datum defines the complex reductive group entirely. It consists of a
|RootSystem| that describes the roots and coroots in the lattices they span
themselves, and of additional data that correspond to embeddings of these
lattices into mutually dual free abelian groups (weight and coweight lattices).

The rank |d_rank| is that of the weight and coweight lattices, the root system
itself has rank |semisimpleRank()| which may be smaller. The roots and coroots
are stored in compact form in the |RootSystem|, and again as represented in
the weight and coweight lattices, for efficiency of retrieval under this form.
Also constructed are various useful auxiliary things, like d_twoRho (the sum
of the positive roots).

The code is designed to make it preferable always to refer to a root by its
number (index in the root system), for which we use the type name |RootNbr|.
*/

class RootDatum
: public RootSystem
{

 private:
 /*!
\brief Names describing the  bits of the bitset d_status.

The last enum numFlags is there as a standard programming trick. Its
value - in this case 2 - is one-past-the-last meaningful bit, for
use by accessors.
*/
  enum StatusFlagNames { IsAdjoint, IsSimplyConnected, numFlags };


  typedef bitset::BitSet<numFlags> Status;

/*!
\brief  Rank of the root datum.
*/
  size_t d_rank;

  WeightList d_roots; //!< Full list of roots.
  CoweightList d_coroots; //!< Full list of coroots.
  WeightList weight_numer; //!< Fundamental weight numerators.
  CoweightList coweight_numer; //!< Fundamental coweight numerators.
  CoweightList d_radicalBasis; //!< Basis for orthogonal to coroots.
  WeightList d_coradicalBasis; //!< Basis for orthogonal to roots.

/*!
\brief Sum of the positive roots.
*/
  Weight d_2rho;
  Coweight d_dual_2rho;

  int Cartan_denom; //!< Denominator for (co)|weight_numer|

/*!\brief BitSet recording whether the root datum is adjoint/simply connected.

  "Adjoint" here means that the center of the complex group determined by the
  root datum is connected. "Simply connected" means that the derived group of
  that complex group is simply connected. These two properties are exchanged
  by passage to the dual root datum.
*/
  Status d_status;


 public:

// constructors and destructors

 RootDatum()
   : RootSystem(int_Matrix(0,0))
   , d_rank(0)
  {}

  explicit RootDatum(const PreRootDatum&);

  RootDatum(const RootDatum&, tags::DualTag);

  RootDatum(int_Matrix&, const RootDatum&, tags::DerivedTag);

  RootDatum(int_Matrix&, const RootDatum&, tags::SimplyConnectedTag);

  RootDatum sub_datum(const RootNbrList& generators) const; // pseudo-constructor

// accessors

  const RootSystem& root_system() const { return *this; } // base object ref
  size_t rank() const { return d_rank; }
  size_t semisimpleRank() const { return RootSystem::rank(); }

// root list access

  WeightList::const_iterator beginRoot() const
    { return d_roots.begin(); }

  WeightList::const_iterator endRoot() const
    { return d_roots.end(); }

  CoweightList::const_iterator beginCoroot() const
    { return d_coroots.begin(); }

  CoweightList::const_iterator endCoroot() const
    { return d_coroots.end(); }

  CoweightList::const_iterator beginRadical() const
    { return d_radicalBasis.begin(); }

  CoweightList::const_iterator endRadical() const
    { return d_radicalBasis.end(); }

  WeightList::const_iterator beginCoradical() const
    { return d_coradicalBasis.begin(); }

  WeightList::const_iterator endCoradical() const
    { return d_coradicalBasis.end(); }

  // below |WRootIterator| is legacy; it equals |WeightList::const_iterator|
  WRootIterator beginSimpleRoot() const // simple roots start halfway
    { return beginRoot()+numPosRoots(); }

  WRootIterator endSimpleRoot() const // and end after |semisimpleRank()|
    { return beginSimpleRoot()+semisimpleRank(); }

  WRootIterator beginPosRoot() const // positive roots start halfway
    { return  beginSimpleRoot(); }

  WRootIterator endPosRoot() const // an continue to the end
    { return endRoot(); }

  WRootIterator beginSimpleCoroot() const
    { return beginCoroot()+numPosRoots(); }

  WRootIterator endSimpleCoroot() const
    { return beginSimpleCoroot()+semisimpleRank(); }

  WRootIterator beginPosCoroot() const // positive coroots start halfway
    { return  beginSimpleCoroot(); }

  WRootIterator endPosCoroot() const
    { return endCoroot(); }


  bool isRoot(const Weight& v) const // ask this of a weight
    { return permutations::find_index(d_roots,v) != d_roots.size(); }

  bool isSemisimple() const { return d_rank == semisimpleRank(); }

  const Weight& root(RootNbr i) const
    { assert(i<numRoots()); return d_roots[i]; }

  const Weight& simpleRoot(weyl::Generator i) const
    { assert(i<semisimpleRank()); return *(beginSimpleRoot()+i); }

  const Weight& posRoot(size_t i) const
    { assert(i<numPosRoots()); return *(beginPosRoot()+i); }

  RootNbr rootNbr(const Root& r) const
    { return permutations::find_index(d_roots,r); }


  const Coweight& coroot(RootNbr i) const
    { assert(i<numRoots()); return d_coroots[i]; }

  const Coweight& simpleCoroot(weyl::Generator i) const
    { assert(i<semisimpleRank()); return *(beginSimpleCoroot()+i); }

  const Coweight& posCoroot(size_t i) const
    { assert(i<numPosRoots()); return  *(beginPosCoroot()+i); }

  // avoid inlining of the following to not depend on rational vector
  RatWeight fundamental_weight(weyl::Generator i) const;
  RatCoweight fundamental_coweight(weyl::Generator i) const;

// other accessors

/*!
\brief Tells whether the rootdatum is the rootdatum of an adjoint group.

  NOTE: we define a reductive group to be adjoint if its center is
  connected.  An equivalent condition is that the derived group
  of the dual group is simply connected.
*/
  bool isAdjoint() const { return d_status[IsAdjoint]; }

/*!
\brief Tells whether the rootdatum is the rootdatum of a simply connected
  group.

  NOTE: this is the dual condition to being adjoint: it means
  that the derived group is simply connected.  An equivalent condition
  is that the center of the dual group is connected.
*/
  bool isSimplyConnected() const { return d_status[IsSimplyConnected]; }


  const Weight& twoRho() const { return d_2rho; }
  const Coweight& dual_twoRho() const { return d_dual_2rho; }


  int scalarProduct(const Weight& v, RootNbr j) const
    { return v.dot(coroot(j)); }

  bool isOrthogonal(const Weight& v, RootNbr j) const
    { return v.dot(coroot(j))==0; }

  bool isOrthogonal(RootNbr alpha, RootNbr j) const
    { return isOrthogonal(root(alpha),j); }

  int cartan(weyl::Generator i, weyl::Generator j) const
    { return simpleRoot(i).dot(simpleCoroot(j)); }

  //!\brief  Applies to v the reflection about root alpha.
  void reflect(Weight& lambda, RootNbr alpha) const
    { lambda -= root(alpha)*lambda.dot(coroot(alpha)); }
  //!\brief  Applies reflection about coroot |alpha| to a coweight
  void coreflect(Coweight& co_lambda, RootNbr alpha) const
    { co_lambda -= coroot(alpha)*co_lambda.dot(root(alpha)); }

  Weight reflection(Weight lambda, RootNbr alpha) const
    { reflect(lambda,alpha); return lambda; }
  Coweight coreflection(Coweight co_lambda, RootNbr alpha) const
    { coreflect(co_lambda,alpha); return co_lambda; }

  void simpleReflect(Weight& v, weyl::Generator i) const
    { reflect(v,simpleRootNbr(i)); }
  void simpleCoreflect(Coweight& v, weyl::Generator i) const
    { coreflect(v,simpleRootNbr(i)); }

  Weight simpleReflection(Weight lambda, weyl::Generator i) const
    { simpleReflect(lambda,i); return lambda; }
  Coweight simpleCoreflection(Coweight co_lambda, weyl::Generator i) const
    { simpleCoreflect(co_lambda,i); return co_lambda; }

  WeylWord to_dominant(Weight lambda) const; // call by value
  void act(const WeylWord& ww,Weight& lambda) const
    {
      for (weyl::Generator i=ww.size(); i-->0; )
	simpleReflect(lambda,ww[i]);
    }
  Weight image_by(const WeylWord& ww,Weight lambda) const
    { act(ww,lambda); return lambda; }

  // here any matrix permuting the roots is allowed, e.g., root_reflection(r)
  permutations::Permutation rootPermutation(const WeightInvolution& q) const;
  // extend diagram automorphism to permutation of all roots

  WeightInvolution root_reflection(RootNbr r) const;
  WeightInvolution simple_reflection(weyl::Generator i) const
    { return root_reflection(simpleRootNbr(i)); }

  WeightInvolution matrix(const WeylWord& ww) const;

  WeylWord reflectionWord(RootNbr r) const;

  // express root in basis of simple roots
  int_Vector inSimpleRoots(RootNbr alpha) const { return root_expr(alpha); }
  // express coroot in basis of simple coroots
  int_Vector inSimpleCoroots(RootNbr alpha) const { return coroot_expr(alpha); }

  Weight twoRho(const RootNbrList&) const;
  Weight twoRho(const RootNbrSet&) const;
  Coweight dual_twoRho(const RootNbrList&) const;
  Coweight dual_twoRho(const RootNbrSet&) const;


  WeylWord word_of_inverse_matrix(const WeightInvolution&)
    const;

// manipulators

  void swap(RootDatum&);

// private methods used during construction
 private:

  void fillStatus();


}; // |class RootDatum|


} // |namespace rootdata|

} // |namespace atlas|

#endif
