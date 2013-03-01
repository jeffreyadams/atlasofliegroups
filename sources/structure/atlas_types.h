/*!
\file
\brief Forward declarations of classes and types for atlas namespace.

*/
/*
  This is atlas_types.h

  Copyright (C) 2011 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef ATLAS_TYPES_H  /* guard against multiple inclusions */
#define ATLAS_TYPES_H

/*

 This module defines types in the (global) atlas namespace that can be
 therefore used in short form by all modules that include this file.

 It is just a forward-declaration file, and is intended to replace several
 such files. One still needs to include the related header files for these
 types to be complete.

 */

#include <vector>
#include <functional> // for |std::less|

#include "constants.h"

// |ndebug_use| can be used to ensure variables in assert statements are "used"
#ifdef NDEBUG
#define ndebug_use(v) static_cast<void>(v)
#else
#define ndebug_use(v)
#endif

/******** forward type declarations ******************************************/

namespace atlas {

/* We begin with copying the contents of the *_fwd.h files in the utilities
   subdirectory. We justify this abject duplication of information by the fact
   that it avoids opening many very short files many times during
   compitlation. Any *_fwd.h files in the structure and gkmod subirectories
   were simply replaced by this file, avoiding duplication, but the utilities
   modules have the ambition of being reusable independently of the rest of
   the atlas library. In fact this file should always be included in any atlas
   header file other than from the utilitites subdierectory, and alway before
   any header files from that subdirectory, so the definitions here will be
   the only ones seen when compiling the atlas library.
 */
  namespace set {
    typedef size_t Elt;
    typedef std::vector<Elt> EltList;
  }

  namespace bitset {
    template<size_t n> class BitSet;
  }
  using bitset::BitSet;
  typedef BitSet<constants::RANK_MAX> RankFlags;
  typedef BitSet<2*constants::RANK_MAX> TwoRankFlags;
  typedef std::vector<RankFlags> RankFlagsList;

  namespace bitmap { class BitMap; }
  using bitmap::BitMap;

  namespace arithmetic {
    typedef long long int Numer_t;
    typedef unsigned long long int Denom_t;
    class Rational;
    class Split_integer;
  }
  using arithmetic::Rational;
  typedef std::vector<Rational> RationalList;
  using arithmetic::Split_integer;

  namespace matrix {
    template<typename C> class Vector;
    template<typename C> class Matrix_base;
    template<typename C> class Matrix;
  }
  namespace ratvec { template<typename C> class RationalVector; }

  namespace permutations { struct Permutation; }
  using permutations::Permutation;

  namespace partition {
    class Partition;
    class PartitionIterator;
  }
  using partition::Partition;

  namespace poset {
    class Poset;
    typedef std::pair<set::Elt,set::Elt> Link;
  }
  using poset::Poset;

  namespace graph {
    typedef set::Elt Vertex;
    typedef std::vector<Vertex> VertexList;
    typedef Vertex Edge;
    typedef std::vector<Edge> EdgeList;
    class OrientedGraph;
  }
  using graph::OrientedGraph;

  namespace hashtable{
    template <class Entry, typename Number> class HashTable;
  }
  using hashtable::HashTable;

  namespace free_abelian {
    template<typename T, typename C=long int, typename Compare=std::less<T> >
      struct Free_Abelian;
  }
  using free_abelian::Free_Abelian;

  namespace polynomials {
    template<typename C> class Polynomial;
    template<typename C> class Safe_Poly;
    typedef size_t Degree; // exponent range; not stored.
  }
  using polynomials::Polynomial;
  namespace size {
    template<typename C> class SizeType;
    typedef signed char BaseType;
    typedef unsigned char UnsignedBaseType;
    typedef SizeType<BaseType> Size;
  }

// we should now refrain from subsequently reading the original forward files
#define SET_H
#define BITSET_FWD_H
#define BITMAP_FWD_H
#define ARITHMETIC_FWD_H
#define MATRIX_FWD_H
#define RATVEC_FWD_H
#define PERMUTATIONS_FWD_H
#define PARTITIONS_FWD_H
#define POSET_FWD_H
#define HASHTABLE_FWD_H
#define FREE_ABELIAN_FWD_H
#define POLYNOMIALS_FWD_H
#define SIZE_FWD_H

  // interpetationless terminology
  typedef matrix::Vector<int> int_Vector;
  typedef matrix::Matrix<int> int_Matrix;
  typedef std::vector<int_Vector> int_VectorList;

  // when related to a root system, these alternatives can be used
  typedef int_Vector Weight;
  typedef int_Vector Coweight;
  typedef ratvec::RationalVector<arithmetic::Numer_t> RatWeight;
  typedef ratvec::RationalVector<arithmetic::Numer_t> RatCoweight;
  typedef matrix::Vector<arithmetic::Numer_t> Ratvec_Numer_t;
  typedef int_Matrix WeightInvolution;
  typedef int_Matrix CoweightInvolution;

  typedef std::vector<Weight> WeightList;
  typedef std::vector<Weight> CoweightList;
  typedef std::vector<RatWeight> RatWeightList;

  // more generic lattice-ralated terms
  typedef int LatticeCoeff; // the instance type of |Vector| and |Matrix|
  typedef std::vector<LatticeCoeff> CoeffList; // no vector arithmetic here
  typedef matrix::Vector<LatticeCoeff> LatticeElt;
  typedef std::vector<LatticeElt> latticeEltList;
  typedef matrix::Matrix<LatticeCoeff> LatticeMatrix;

  namespace bitvector {
    template<size_t> class BitVector;
    template<size_t> class BitVectorList;
    template<size_t> class BitMatrix;
  }
  using bitvector::BitVector;
  using bitvector::BitVectorList;
  using bitvector::BitMatrix;
  typedef BitVector<constants::RANK_MAX> SmallBitVector;
  typedef BitVectorList<constants::RANK_MAX> SmallBitVectorList;
  typedef BitVector<constants::RANK_MAX+1> BinaryEquation;
  typedef BitVectorList<constants::RANK_MAX+1> BinaryEquationList;
  typedef BitMatrix<constants::RANK_MAX> BinaryMap;

  namespace subquotient {
    template<size_t dim> class Subspace;
    template<size_t dim> class Subquotient;
  }
  typedef subquotient::Subspace<constants::RANK_MAX> SmallSubspace;
  typedef subquotient::Subquotient<constants::RANK_MAX> SmallSubquotient;

  namespace lietype {
    struct SimpleLieType;
    struct LieType;
    struct InnerClassType;
    struct Layout;
  }
  using lietype::SimpleLieType;
  using lietype::LieType;
  using lietype::InnerClassType;

  namespace prerootdata { class PreRootDatum; }
  using prerootdata::PreRootDatum;

  namespace rootdata {
    class RootSystem;
    class RootDatum;
    typedef int_Vector Root;
    typedef WeightList::const_iterator WRootIterator;
  }
  using rootdata::RootSystem;
  using rootdata::RootDatum;
  typedef unsigned short RootNbr;
  typedef std::vector<RootNbr> RootNbrList;
  typedef bitmap::BitMap RootNbrSet;

  namespace weyl {
    class WeylInterface;
    typedef WeylInterface Twist; // use the same implementation
    typedef unsigned char Generator; // index of simple root / simple reflection
    struct WeylWord : public std::vector<Generator>  { }; // in weyl namesace

    class WeylElt;
    class WeylGroup;
    class TwistedWeylGroup;

    typedef WeylElt TwistedInvolution;
    typedef std::vector<WeylElt> WeylEltList;
    typedef std::vector<TwistedInvolution> TwistedInvolutionList;
    struct TI_Entry; // for |TwistedInvolution|s in hash tables

    typedef std::vector<signed char> InvolutionWord;
  }
  using weyl::WeylWord;
  using weyl::WeylElt;
  using weyl::WeylGroup;
  using weyl::TwistedWeylGroup;
  using weyl::TwistedInvolution;
  using weyl::WeylEltList;
  using weyl::TwistedInvolutionList;

  namespace y_values {
    class TorusElement;
    struct y_entry;
  }
  using y_values::TorusElement;
  using y_values::y_entry;

  namespace tits {
    typedef SmallBitVector TorusPart;
    class GlobalTitsElement;
    class GlobalTitsGroup;
    class SubTitsGroup;
    class TitsElt;
    class TitsGroup;
    struct TE_Entry;
    class TitsCoset;
    class EnrichedTitsGroup;
  }
  using tits::TorusPart;
  using tits::GlobalTitsElement;
  using tits::GlobalTitsGroup;
  using tits::TitsElt;
  using tits::TitsGroup;
  using tits::TitsCoset;

  namespace involutions {
    class InvolutionData;
    class InvolutionTable;
    class Cartan_orbit;
    class Cartan_orbits;
  }
  using involutions::InvolutionData;
  using involutions::InvolutionTable;
  typedef unsigned int InvolutionNbr;
  using involutions::Cartan_orbit;
  using involutions::Cartan_orbits;

  namespace gradings {
    class Status;
    struct GradingCompare;
  }
  typedef RankFlags Grading;
  typedef std::vector<Grading> GradingList;

  namespace tori { class RealTorus; }

  namespace topology { class Connectivity; }

  namespace subsystem {class SubSystem; class SubSystemWithGroup; }
  using subsystem::SubSystem;
  using subsystem::SubSystemWithGroup;


  typedef unsigned short CartanNbr; // index of Cartan class
  typedef unsigned short RealFormNbr; // index used in |ComplexReductiveGroup|
  typedef std::vector<RealFormNbr> RealFormNbrList;

  namespace cartanclass {
    class Fiber;
    class CartanClass;
    typedef unsigned int FiberElt;	// element of the fiber group
    typedef unsigned int AdjointFiberElt;  // element of adjoint fiber (-group)
    typedef unsigned short fiber_orbit; // # of W_imag orbit in fiber group
    typedef unsigned short adjoint_fiber_orbit; // same for adjoint fiber group
    typedef unsigned short square_class; // identifies a class of real forms
    typedef std::pair<fiber_orbit,square_class> StrongRealFormRep;
  }
  using cartanclass::Fiber;
  using cartanclass::CartanClass;

  namespace complexredgp { class ComplexReductiveGroup; }
  using complexredgp::ComplexReductiveGroup;

  namespace realredgp { class RealReductiveGroup; }
  using realredgp::RealReductiveGroup;

  namespace bruhat { class BruhatOrder; }
  using bruhat::BruhatOrder;

  namespace kgb {
    class KGB_base;
    struct KGB_elt_entry;
    class GlobalFiberData;
    class global_KGB;
    class KGB;
    class subsys_KGB;
    typedef RankFlags DescentSet;
    class KGP;
  }
  using kgb::KGB_base;
  using kgb::KGB_elt_entry;
  using kgb::global_KGB;
  using kgb::KGB;
  typedef unsigned int KGBElt;
  typedef std::vector<KGBElt> KGBEltList;

  typedef std::pair<KGBElt,KGBElt> KGBEltPair;
  typedef std::vector<KGBEltPair> KGBEltPairList;

  using kgb::DescentSet;
  static const KGBElt UndefKGB = ~0u;

  typedef unsigned int Parabolic;
  typedef unsigned int KGPElt;

  namespace descents { class DescentStatus; }
  using descents::DescentStatus;
  typedef std::vector<DescentStatus> DescentStatusList;

  namespace blocks {
    class Block_base;
    class Block;
    class param_block;
    class gamma_block;
    class non_integral_block;
  }
  using blocks::Block_base;
  using blocks::Block;
  using blocks::param_block;
  using blocks::non_integral_block;
  typedef unsigned int BlockElt;
  typedef std::vector<BlockElt> BlockEltList;
  typedef std::pair<BlockElt,BlockElt> BlockEltPair;
  typedef std::vector<BlockEltPair> BlockEltPairList;
  static const BlockElt UndefBlock = ~0u;

  namespace klsupport { class KLSupport; }
  namespace wgraph {
    class WGraph;
    class DecomposedWGraph;
    typedef std::vector<unsigned short> WCoeffList;
  }
  namespace kl {
    class KLContext;
    typedef unsigned int KLCoeff;
    typedef polynomials::Safe_Poly<KLCoeff> KLPol;
    typedef unsigned int KLIndex; // $<2^{32}$ distinct polynomials for $E_8$!
    typedef KLCoeff MuCoeff;
    typedef std::vector<KLPol> KLStore;
    typedef KLStore::const_reference KLPolRef;
    typedef std::vector<KLIndex> KLRow;
    typedef std::vector<BlockElt> PrimitiveRow;
  }

  namespace standardrepk {
    class StandardRepK;	// standard representation restricted to K
    typedef std::pair <Weight,RankFlags> HCParam; // free part wrt rho, torsion
    typedef Free_Abelian<StandardRepK> Char;// linear combination
    typedef std::pair<StandardRepK,Char> CharForm;
    typedef std::pair<Weight,TitsElt> RawRep;
    typedef Free_Abelian<RawRep> RawChar;
    typedef Free_Abelian<StandardRepK,Polynomial<int> > q_Char;
    typedef std::pair<StandardRepK,q_Char> q_CharForm;// $q$-$K$-type formula
    typedef Free_Abelian<RawRep,Polynomial<int> >Raw_q_Char;
    typedef unsigned int seq_no; // sequence number of stored standard rep|K
    typedef unsigned int level; // unsigned LatticeCoeff
    struct Cartan_info;
    struct bitset_entry;
    class SRK_context;
    class graded_compare;// utility class for comparing by degree first
    class KhatContext;
    class HechtSchmid;	// Hecht-Schmid identity
    class PSalgebra;    // parabolic subalgebra
  }
  using standardrepk::StandardRepK;
  using standardrepk::SRK_context;

  namespace repr {
    class StandardRepr;	// triple $(x,\lambda,\gamma)$
    class Rep_context;	// support class for interpreting |StandardRepr|
    class Rep_table;	// storage class for |StandardRepr| compuations
  }
  using repr::StandardRepr;
  using repr::Rep_context;
  using repr::Rep_table;

  namespace realform_io { class Interface; } // maps internals to names
  namespace complexredgp_io { class Interface; } // a pair of the above
  namespace realweyl {
    class RealWeyl; // for computing real Weyl groups
    class RealWeylGenerators;
    struct PermutationGenerators;
  }

  namespace input {
    class InputBuffer;
#ifndef NREADLINE
    class HistoryBuffer;
#else
    typedef InputBuffer HistoryBuffer;
#endif
  }


} // |namespace atlas|

#endif
