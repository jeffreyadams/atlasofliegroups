/*!
\file
  This is abelian.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#ifndef ABELIAN_H  /* guard against multiple inclusions */
#define ABELIAN_H

#include <limits>
#include <set>

#include "abelian_fwd.h"

#include "bitmap.h"
#include "matrix.h"

/******** type declarations **************************************************/

namespace atlas {

namespace abelian {

typedef matrix::Matrix<unsigned long> Endomorphism;

}

/******** function declarations **********************************************/

namespace abelian {

  void basis(latticetypes::WeightList&, const bitmap::BitMap&,
	     const FiniteAbelianGroup&);

  void coset(bitmap::BitMap&, const bitmap::BitMap&, GrpNbr,
	     const FiniteAbelianGroup&);

  const bitmap::BitMap& cycGenerators(const FiniteAbelianGroup&);

  void generateSubgroup(bitmap::BitMap&, GrpNbr, const FiniteAbelianGroup&);

  void generators(GrpNbrList&, const bitmap::BitMap&,
		  const FiniteAbelianGroup&);

  bool isElementaryAbelian(const std::vector<unsigned long>&);

  void quotReps(bitmap::BitMap&, const bitmap::BitMap&,
		const FiniteAbelianGroup&);

  void toArray(GrpArr&, GrpNbr, const GroupType&);

  void toArray(GrpArr&, const latticetypes::Weight&, const GroupType&);

  void toEndomorphism(Endomorphism&, const latticetypes::LatticeMatrix&,
		      const FiniteAbelianGroup&);

  GrpNbr toGrpNbr(const GrpArr&, const GroupType&);

  void transpose(Endomorphism&, const FiniteAbelianGroup&);

}

/******** type definitions ***************************************************/

namespace abelian {

class FiniteAbelianGroup {
  /*!  
  A product of finite cyclic groups of orders d_type[0], d_type[1],
  d_type[2]...  It is assumed that d_type[0] divides d_type[1], which
  divides d_type[2], and so on; these orders are therefore an
  invariant of the group. The vector d_cotype has in its jth
  coordinate the quotient d_type[last]/d_type[j].  This means that the
  jth factor of the group may be written as the quotient of the cyclic
  group of order d_type[last] by the subgroup of order d_cotype[j].
  The integer d_size is the order of the group.
  */

 protected:
  /*!
Cardinality of the group.
  */
  unsigned long d_size;
  /*!
  Sizes of the canonical cyclic factors, arranged so that each divides
  the next.
  */
  GroupType d_type;
  /*!
  If m is the largest order of an element (equal to d_type[last]), then
  the group is a product of quotients of Z/mZ, by the subgroups of
  order d_cotype[j]; each of these orders divides its predecessor, and
  d_cotype[last] is equal to 1.
  */  
  GroupType d_cotype;

 public:

// constructors and destructors
  FiniteAbelianGroup() 
    {}

  explicit FiniteAbelianGroup(const std::vector<unsigned long>&);

  ~FiniteAbelianGroup() 
    {}

// copy and assignment: defaults are ok

// type conversions
  void toArray(GrpArr& a, GrpNbr x) const {
    return abelian::toArray(a,x,d_type);
  }

  void toArray(GrpArr& a, const latticetypes::Weight& v) const {
    return abelian::toArray(a,v,d_type);
  }

  GrpNbr toGrpNbr(const GrpArr& a) const {
    return abelian::toGrpNbr(a,d_type);
  }

  void toWeight(latticetypes::Weight&, const GrpArr&) const;

  void toWeight(latticetypes::Weight&, GrpNbr) const;

// accessors
  bool operator== (const FiniteAbelianGroup& A) const {
    return d_type == A.d_type;
  }

  bool operator!= (const FiniteAbelianGroup& A) const {
    return not operator==(A);
  }

  GrpNbr add(GrpNbr, GrpNbr) const;

  GrpNbr add(GrpNbr, const GrpArr&) const;

  GrpArr& add(GrpArr&, const GrpArr&) const;

  GrpArr& add(GrpArr&, GrpNbr) const;

  unsigned long annihilator() const;

  const GroupType& cotype() const {
    return d_cotype;
  }

  GrpNbr leftApply(GrpNbr, const Endomorphism&) const;

  GrpArr& leftApply(GrpArr&, const Endomorphism&) const;

  GrpNbr minus(GrpNbr) const;

  unsigned long order(GrpNbr) const;

  unsigned long order(const bitmap::BitMap&, GrpNbr) const;

  unsigned long pairing(const GrpArr&, const GrpArr&) const;

  unsigned long pairing(const GrpArr&, const GrpArr&, unsigned long) const;

  GrpNbr prod(GrpNbr, unsigned long) const;

  unsigned long rank() const {
    return d_type.size();
  }

  unsigned long size() const {
    return d_size;
  }

  GrpNbr subtract(GrpNbr, const GrpNbr) const;

  GrpArr& subtract(GrpArr&, const GrpArr&) const;

  const GroupType& type() const {
    return d_type;
  }
};

/*!
  A Homomorphism object represents a homomorphism from a group of type d_source
  to a group of type d_dest, where the elements are represented as arrays.

  As we know, for such a homomorphism to be well defined, the matrix entry
  (i,j) should be a multiple of d_dest[i]/g, with g = 
  gcd(d_dest[i],d_source[j]).

  Here we take a different approach. We set M a common multiple of the
  annihilators of d_source and d_dest. Then
  we allow arbitrary coefficients in Z/M. The corresponding matrix is
  interpreted as follows. It is defined only for those tuples (x_1,...,x_m)
  s.t. for all i, sum_j a_ij q_j x_j is of index t_i in Z/M, where (q_j)
  is the cotype of d_source w.r.t. M; and then the corresponding value y_i is
  (sum_j a_ij q_j x_j)/(M/t_i). It is easy to see that this definition is
  in fact independent of the choice of M (of course we will take the lcm
  of the annihilators.)
*/

class Homomorphism {

 private:

  GroupType d_source;
  GroupType d_dest;
  GroupType d_cosource;
  GroupType d_codest;
  unsigned long d_annihilator;
  matrix::Matrix<unsigned long> d_matrix;

 public:

// constructors and destructors
  Homomorphism(const std::vector<GrpArr>&, const GroupType&, const GroupType&);

// accessors
  void apply(GrpArr&, const GrpArr&) const;

  GrpNbr apply(GrpNbr) const;

  bool defined(const GrpArr&) const;

  bool defined(GrpNbr) const;
};

}

}

#endif
