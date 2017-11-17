/*
  This is lietype.h

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2017 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// The type |LieType| and related types, and operations for them


#ifndef LIETYPE_H  /* guard against multiple inclusions */
#define LIETYPE_H


#include <stdexcept>
#include <vector>

#include "../Atlas.h"

#include "permutations.h" // needed in the |Layout| structure

namespace atlas {

namespace lietype {

// type declaration
struct Layout;


/******** constant declarations **********************************************/
  /*
  Used by the interaction in interactive_lietype.cpp to tell the user
  what input is expected for a Lie type. Also used in |SimpleLieType|
  */
const char* const typeLetters = "ABCDEFGT";

  /*!
  Used by the interaction in interactive_lietype.cpp to tell the user
  what input is expected for an inner class type.

  The letter "C" stands for "complex"; it goes with two adjacent equal
  Lie types, to say that the group has the corresponding complex
  simple factor.  (For example, a Lie type of A4.A4 and inner type C
  refers to SL(5,C).)

  The letter "c" stands for "compact"; it designates the inner class
  of forms in which rank(G) = rank(K) (so that there is a compact
  Cartan subgroup.)

  The letter "s" stands for "split," and designates the inner class
  containing the split real form of G.

  The letter "u" stands for "unequal rank," and designates an inner
  class of forms with no compact Cartan and (if possible) no split
  Cartan. This inner class is allowed only in A_n, n > 1, D_n, and
  E_6.  In types A, E6, and D_{2n+1} the unique unequal rank inner
  class is split, and the inner class designator "u" is then mapped to
  "s."  (This is important, because the code implementing the Cartan
  involution in lietype.cpp would not work in types A and E with the
  designator "u.")  In type D_{2n}, the inner class "u" corresponds to
  the Cartan involution exchanging the two simple roots at the
  branched end of the Dynkin diagram.  (In D_4 there are three
  involutions fitting this description.  The software picks just one
  of them.  Since the three differ by an outer automorphism of D_4,
  this loss of generality is harmless.)  The real forms in the inner
  class "u" for D_{2n} are the special orthogonal groups of signature
  (odd,odd).

  The letter "e" stands for "equal rank," and is mapped in every case to
  "c."
  */
const char* const innerClassLetters = "Ccesu";

enum class simple_ict { equal_rank, unequal_rank, complex };

struct SimpleLieType : public std::pair<TypeLetter,size_t>
{ // there are no additional data members
  typedef std::pair<TypeLetter,size_t> base;
  SimpleLieType(TypeLetter t,size_t rank) : base(t,rank) {}
  TypeLetter type() const { return base::first; }
  TypeLetter& type() { return base::first; }
  size_t rank() const { return base::second; }
  size_t& rank() { return base::second; }
  size_t semisimple_rank() const { return type()=='T' ? 0 : rank(); }
  int Cartan_entry(size_t i,size_t j) const;
  int_Matrix Cartan_matrix() const; // square of size rank; null for type $T_n$
  int_Matrix transpose_Cartan_matrix() const;
};

struct LieType : public std::vector<SimpleLieType>
{ // there are no additional data members
  typedef std::vector<SimpleLieType> base;
  LieType() : base() {}
  LieType(const base& b) : base(b) {}

  size_t rank() const;
  size_t semisimple_rank() const;
  int Cartan_entry(size_t i,size_t j) const;
  int_Matrix Cartan_matrix() const; // square of sise |rank()| (maybe zero rows)
  int_Matrix transpose_Cartan_matrix() const;
  int_VectorList Smith_basis(CoeffList& invf) const;

  LieType& append (const LieType& tail)
  { insert(end(),tail.begin(),tail.end()); return *this; }
};

// the follwing rather empty definition serves mainly to make |InnerClassType|
// a genuine member of |namespace lietype| for argument-dependent lookup
struct InnerClassType : public std::vector<TypeLetter>
{ // there are no additional data members
  typedef std::vector<TypeLetter> base;
  InnerClassType() : base() {}
};

/* the |Layout| structure collects a |LieType| and an |InnerClassType|, as
   specified by the user. However dualizing may create Lie types with
   non-standard node labelings (for $G_2$, $F_4$), and in such cases an
   involution of the diagram needs to be indicated, whence the |d_perm| field.
   This is ultimately used (only) to correctly associate a real form name
   (recognised in standard diagram labelling) from a special representative
   grading of the real form (only one bit set), in |output::printType|;
   |d_perm| maps standard (Bourbaki) diagram numbers to simple root indices.
*/
struct Layout
{
  LieType d_type;
  InnerClassType d_inner;
  Permutation d_perm; // of diagram wrt usual order in |d_type|

// constructors and destructors

Layout() : d_type(), d_inner(), d_perm() {} // used in both Fokko and atlas

/* In the old Fokko interface, the Lie type is first provided, and the inner
   class type is later added; this constructor defaults permutation to identity
*/
  Layout(const LieType& lt)
    :d_type(lt),d_inner(),d_perm(lt.rank(),1) {}

// The inner class can also be added right away; permutation defaults to id
  Layout(const LieType& lt, const InnerClassType& ict)
    :d_type(lt),d_inner(ict),d_perm(lt.rank(),1) {}

// Finally we can specify all fieds upon construction
  Layout(const LieType& lt, const InnerClassType& ict, const Permutation& pi)
    :d_type(lt),d_inner(ict),d_perm(pi) {}


}; // |struct Layout|


struct ext_gen // generator of extended Weyl group
{
  enum { one, two, three } type;
  weyl::Generator s0,s1;
  WeylWord w_kappa;

  explicit ext_gen (weyl::Generator s)
    : type(one), s0(s), s1(~0), w_kappa() { w_kappa.push_back(s); }
  ext_gen (bool b) = delete; // defuse implicit conversion
  ext_gen (bool commute, weyl::Generator s, weyl::Generator t)
  : type(commute ? two : three), s0(s), s1(t)
  { w_kappa.reserve(commute ? 2 : 3);
    w_kappa.push_back(s);  w_kappa.push_back(t);
    if (not commute) w_kappa.push_back(s);
  }

  int length() const { return type+1; } // also equals |w_kappa.size()|
};




/******** function declarations **********************************************/

  bool checkRank(const TypeLetter&, size_t);

  WeightInvolution simple_involution(const SimpleLieType& slt, simple_ict tp);

  // involution (permutation) matrix for possibly renumbered Dynkin diagram
  WeightInvolution involution(const Layout& lo);

  // permutation matrix for |ict| in simply connected |lt|, Bourbaki order
  WeightInvolution involution(const LieType& lt,
			      const InnerClassType& ict);
  // return |involution(Layout(lt,ict))|, but |WeightInvolution| incomplete here

  LieType dual_type(LieType lt);

  InnerClassType dual_type(InnerClassType, const LieType& lt);

  Layout dual(const Layout& lo);

}

}

#endif
