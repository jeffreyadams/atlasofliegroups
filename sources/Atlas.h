/*
  This is Atlas.h

  Copyright (C) 2011 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef ATLAS_H  /* guard against multiple inclusions */
#define ATLAS_H

/*    Forward declarations of classes and types for atlas namespace.

 This module defines types in the (global) atlas namespace that can be
 therefore used in short form by all modules that include this file.

 It is mostly just a forward-declaration file, and is intended to replace
 several such files. We also make system-wide |typedef|s here. One still needs
 to include the related ordinary header files for these types to be complete.

 This file should be included when compiling the Atlas object files everywhere
 BEFORE any of the utility forward files (the ones whose inclusion is excleded
 below) are invoked; this can be accomplished by including it before any other
 Atlas files in each non-utility header file. (Header or other files in the
 utility sub-directory itself should not include or otherwise refer to this
 file, so as to retain their independence of the rest of the Atlas project.) In
 this manner those object files avoid including both such a forward file and its
 information duplicated here. Since the utility object files necessarily do use
 the information from those forward files, their contents and that of Atlas.h
 must nevertheless always remain in sync, lest problems occur at link time.

 */

// exclude forward files from utilities subdirectory whose contents are copied here
#define ARITHMETIC_FWD_H
#define BITMAP_FWD_H
#define BITSET_FWD_H
#define FREE_ABELIAN_FWD_H
#define GRAPH_FWD_H
#define HASHTABLE_FWD_H
#define MATRIX_FWD_H
#define PARTITION_FWD_H
#define PERMUTATIONS_FWD_H
#define POLYNOMIALS_FWD_H
#define POSET_FWD_H
#define RATVEC_FWD_H
#define SIZE_FWD_H
#define SL_LIST_FWD_H


#include <vector>
#include <functional> // for |std::less|
#include <memory> // for |std::allocator|
#include <stack> // for our specialisation below
#include <queue> // for our specialisation below


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
   compilation. Any *_fwd.h files in the structure and gkmod subdirectories
   were simply replaced by this file, avoiding duplication, but the utilities
   modules have the ambition of being reusable independently of the rest of
   the Atlas library. In fact this file should always be included in any Atlas
   header file other than from the utilities subdirectory, and always before
   any header files from that subdirectory, so the definitions here will be
   the only ones seen when compiling the Atlas library.
 */
  namespace bitset {
    template<unsigned int n> class BitSet;
  }
  using bitset::BitSet;
  typedef BitSet<constants::RANK_MAX> RankFlags;
  typedef BitSet<2*constants::RANK_MAX> TwoRankFlags;
  typedef std::vector<RankFlags> RankFlagsList;

  namespace bitmap { class BitMap; }
  using bitmap::BitMap;

  namespace containers {

  template <typename Alloc>
    class allocator_deleter;
  template<typename T,typename Alloc = std::allocator<T> >
    struct sl_node;

  template<typename T, typename Alloc = std::allocator<T> >
    class sl_list_const_iterator;
  template<typename T,typename Alloc = std::allocator<T> >
    class sl_list_iterator;

  template<typename T, typename Alloc = std::allocator<T> >
    class weak_sl_list_const_iterator;
  template<typename T, typename Alloc> class
    weak_sl_list_iterator;

  template<typename T,typename Alloc = std::allocator<T> >
    class simple_list;
  template<typename T,typename Alloc = std::allocator<T> >
    class sl_list;

  template<typename T,typename Alloc = std::allocator<T> >
    class mirrored_simple_list; // trivial adapter, allows use with |std::stack|

  template<typename T,typename Alloc = std::allocator<T> >
    class mirrored_sl_list; // trivial adapter, to allow use with |std::stack|

  template<typename T,typename Alloc = std::allocator<T> > class stack;

  template<typename T,typename Alloc = std::allocator<T> > class queue;

  } // |namespace containers|
  using containers::simple_list;
  using containers::sl_list;
  using containers::sl_list_const_iterator;
  using containers::sl_list_iterator;
  using containers::stack;
  using containers::queue;

  namespace arithmetic {
    using Numer_t = long long int;
    using Denom_t = unsigned long long int;
    template<typename I> class Rational;
    using RatNum = Rational<Numer_t>;
    using RatNumList = std::vector<RatNum>;
    class Split_integer;
    class big_int; // defined in bigint.h
    class big_rat; // defined in bigint.h
  }
  using arithmetic::Rational;
  using arithmetic::RatNum;
  using arithmetic::RatNumList;
  using arithmetic::Split_integer;
  using arithmetic::big_int;
  using arithmetic::big_rat;

  namespace matrix {
    template<typename C> class Vector;
    template<typename C> class Matrix_base;
    template<typename C> class Matrix;
    template<typename C> class PID_Matrix;
    template<typename C> class Vector_cref;
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
  }
  using poset::Poset;

  namespace graph {
    typedef unsigned int Vertex; // assume at most some 4 billion vertices
    typedef std::vector<Vertex> EdgeList; // list of targets of outgoing edges
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
    template<typename T, typename C=long int, typename Compare=std::less<T> >
      struct Monoid_Ring;
    template<typename T, typename C=long int, typename Compare=std::less<T> >
      class Free_Abelian_light;
  }
  using free_abelian::Free_Abelian;
  using free_abelian::Monoid_Ring;
  using free_abelian::Free_Abelian_light;

  namespace polynomials {
    template<typename C> class Polynomial;
    template<typename U> class Safe_Poly;
    template<typename C> class PolEntry;
    template<typename U> class SafePolEntry;
    using Degree = unsigned int; // exponent range; not stored.
  }
  using polynomials::Polynomial;

  namespace size {
    template<typename C> class SizeType;
    using BaseType = signed char;
    using Size= SizeType<signed char>;
  }

// we should now refrain from subsequently reading the original forward files
#define BITSET_FWD_H
#define BITMAP_FWD_H
#define SL_LIST_FWD_H
#define ARITHMETIC_FWD_H
#define MATRIX_FWD_H
#define RATVEC_FWD_H
#define PERMUTATIONS_FWD_H
#define PARTITIONS_FWD_H
#define POSET_FWD_H
#define GRAPH_FWD_H
#define HASHTABLE_FWD_H
#define FREE_ABELIAN_FWD_H
#define POLYNOMIALS_FWD_H
#define SIZE_FWD_H

  // interpretationless terminology
  typedef matrix::Vector<int> int_Vector;
  typedef matrix::PID_Matrix<int> int_Matrix;
  typedef std::vector<int_Vector> int_VectorList;
  typedef matrix::Vector_cref<int> int_Vector_cref;
  typedef ratvec::RationalVector<arithmetic::Numer_t> rat_Vector;
  typedef matrix::Vector<arithmetic::Numer_t> Ratvec_Numer_t;

  // when related to a root system, these alternatives can be used
  typedef int_Vector Weight;
  typedef int_Vector Coweight;
  typedef rat_Vector RatWeight;
  typedef rat_Vector RatCoweight;
  typedef int_Matrix WeightInvolution;
  typedef int_Matrix CoweightInvolution;

  typedef std::vector<Weight> WeightList;
  typedef std::vector<Weight> CoweightList;
  typedef std::vector<RatWeight> RatWeightList;

  // more generic lattice-related terms
  typedef int LatticeCoeff; // the instance type of |Vector| and |Matrix|
  typedef std::vector<LatticeCoeff> CoeffList; // no vector arithmetic here
  typedef matrix::Vector<LatticeCoeff> LatticeElt;
  typedef std::vector<LatticeElt> latticeEltList;
  typedef matrix::PID_Matrix<LatticeCoeff> LatticeMatrix;

  namespace bitvector {
    template<unsigned int> class BitVector;
    template<unsigned int> class BitVectorList;
    template<unsigned int> class BitMatrix;
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
    template<unsigned int dim> class Subspace;
    template<unsigned int dim> class Subquotient;
  }
  typedef subquotient::Subspace<constants::RANK_MAX> SmallSubspace;
  typedef subquotient::Subquotient<constants::RANK_MAX> SmallSubquotient;

  namespace lietype {
    struct SimpleLieType;
    struct LieType;
    struct InnerClassType;
    struct Layout;
    typedef char TypeLetter;
    struct ext_gen;
  }
  using lietype::SimpleLieType;
  using lietype::LieType;
  using lietype::InnerClassType;
  using lietype::ext_gen;
  typedef std::vector<ext_gen> ext_gens;

  namespace prerootdata { class PreRootDatum; }
  using prerootdata::PreRootDatum;

  namespace rootdata {
    class RootSystem;
    class RootDatum;
    typedef int_Vector Root;
  }
  using rootdata::RootSystem;
  using rootdata::RootDatum;
  typedef unsigned short RootNbr;
  typedef std::vector<RootNbr> RootNbrList;
  typedef bitmap::BitMap RootNbrSet;

  namespace dynkin { class DynkinDiagram; }
  using dynkin::DynkinDiagram;

  namespace weyl {
    class Twist; // diagram automorphism (in practice always an involution)
    typedef Twist WeylInterface; // no automorphism, but same implementation
    typedef unsigned char Generator; // index of simple root / simple reflection
    class WeylWord;

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
    struct Cartan_orbit;
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

  namespace subsystem {class SubSystem; class SubSystemWithGroup; }
  using subsystem::SubSystem;
  using subsystem::SubSystemWithGroup;


  typedef unsigned short CartanNbr; // index of Cartan class
  typedef unsigned short RealFormNbr; // index used in |InnerClass|
  typedef std::vector<RealFormNbr> RealFormNbrList;

  namespace cartanclass {
    class Fiber;
    class CartanClass;
    typedef SmallBitVector FiberElt;	// element of the fiber group
    typedef SmallBitVector AdjointFiberElt;  // element of adjoint fiber(-group)
    typedef unsigned short fiber_orbit; // # of W_imag orbit in fiber group
    typedef unsigned short adjoint_fiber_orbit; // same for adjoint fiber group
    typedef unsigned short square_class; // identifies a class of real forms
    typedef std::pair<fiber_orbit,square_class> StrongRealFormRep;
  }
  using cartanclass::Fiber;
  using cartanclass::CartanClass;

  namespace innerclass { class InnerClass; }
  using innerclass::InnerClass;

  namespace realredgp { class RealReductiveGroup; }
  using realredgp::RealReductiveGroup;

  namespace bruhat { class BruhatOrder; }
  using bruhat::BruhatOrder;

  namespace kgb {
    class KGB_base;
    struct KGB_elt_entry;
    class global_KGB;
    class KGB;
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

  typedef RankFlags Parabolic; // set of generators defining parabolic subgroup
  typedef unsigned int KGPElt;

  namespace descents { class DescentStatus; }
  using descents::DescentStatus;
  typedef std::vector<DescentStatus> DescentStatusList;

  namespace blocks {
    class Block_base;
    class Block;
    class common_block;
  }
  using blocks::Block_base;
  using blocks::Block;
  using blocks::common_block;
  typedef unsigned int BlockElt;
  typedef std::vector<BlockElt> BlockEltList;
  typedef std::pair<BlockElt,BlockElt> BlockEltPair;
  typedef std::vector<BlockEltPair> BlockEltPairList;
  static const BlockElt UndefBlock = ~0u;

  namespace klsupport { class KLSupport; }
  namespace wgraph {
    class WGraph;
    class DecomposedWGraph;
  }
  namespace kl {
    struct Poly_hash_export;
    class KL_table;
    using KLCoeff = unsigned int;
    using KLPol = polynomials::Safe_Poly<KLCoeff>;
    using KLIndex = unsigned int; // $<2^{32}$ distinct polynomials for $E_8$!
    using MuCoeff = KLCoeff;
  }
  using kl::KLCoeff;
  using kl::KLPol;
  using PosPolEntry = polynomials::SafePolEntry<KLCoeff>;
  using KL_hash_Table = HashTable<PosPolEntry,kl::KLIndex>;

  namespace ext_kl {
    class KL_table;
    using Coeff = int;
    using Pol = Polynomial<Coeff>;
    using PolRef = const Pol&;
    using KLIndex = unsigned int;
  }
  using IntPolEntry = polynomials::PolEntry<ext_kl::Coeff>;
  using ext_KL_hash_Table = HashTable<IntPolEntry,ext_kl::KLIndex>;

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
  using standardrepk::KhatContext;

  namespace repr {
    using level = unsigned int; // for height statistic of K type and parameters

    class StandardRepr;	// triple $(x,\lambda,\gamma)$
    using SR_poly = Free_Abelian<StandardRepr,Split_integer>;
    class StandardReprMod; // represent an $X^*$-shift class of |StandardRepr|
    class Reduced_param; // even coarser representation, mod integral orthogonal
    class Repr_mod_entry;
    class common_context;
    class Rep_context;	// support class for interpreting |StandardRepr|
    class Rep_table;	// storage class for |StandardRepr| computations
  }
  using repr::StandardRepr;
  using repr::SR_poly;
  using repr::StandardReprMod;
  using repr::Reduced_param;
  using repr::common_context;
  using repr::Rep_context;
  using repr::Rep_table;

  namespace ext_block { class ext_block; }

  namespace output { class FormNumberMap; } // maps internals to names
  namespace output { class Interface; } // a pair of the above
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
