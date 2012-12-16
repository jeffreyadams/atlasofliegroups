/*!
\file
  This is abelian.h
*/
/*
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef ABELIAN_H  /* guard against multiple inclusions */
#define ABELIAN_H

#include <limits>
#include <set>

#include "abelian_fwd.h"
#include "arithmetic_fwd.h"

#include "bitmap.h"
#include "matrix.h"

/******** type declarations **************************************************/

namespace atlas {


/******** function declarations **********************************************/

namespace abelian {

  void basis(std::vector<matrix::Vector<int> >&, const bitmap::BitMap&,
	     const FiniteAbelianGroup&);

  void coset(bitmap::BitMap&, const bitmap::BitMap&, GrpNbr,
	     const FiniteAbelianGroup&);

  const bitmap::BitMap& cycGenerators(const FiniteAbelianGroup&);

  void generateSubgroup(bitmap::BitMap&, GrpNbr, const FiniteAbelianGroup&);

  void generators(GrpNbrList&, const bitmap::BitMap&,
		  const FiniteAbelianGroup&);

  bool isElementaryAbelian(const std::vector<arithmetic::Denom_t>&);

  void quotReps(bitmap::BitMap&, const bitmap::BitMap&,
		const FiniteAbelianGroup&);

  void to_array(GrpArr&, GrpNbr, const GroupType&);

  void to_array(GrpArr&, const matrix::Vector<int>&, const GroupType&);

  void toEndomorphism(Endomorphism&, const matrix::Matrix<int>&,
		      const FiniteAbelianGroup&);

  GrpNbr to_GrpNbr(const GrpArr&, const GroupType&);

  void transpose(Endomorphism&, const FiniteAbelianGroup&);

}

/******** type definitions ***************************************************/

namespace abelian {

/*!
  A class representing a finite Abelian group as a product of finite cyclic
  groups of orders d_type[0], d_type[1], d_type[2]... It is assumed that
  d_type[0] divides d_type[1], which divides d_type[2], and so on; these
  orders are therefore invariants that characterise the group. The vector
  d_cotype has in its jth coordinate the quotient d_type[last]/d_type[j]. This
  means that the jth factor of the group may be written as the quotient of the
  cyclic group of order d_type[last] by the subgroup of order d_cotype[j]. The
  integer d_size is the order of the group.
*/
class FiniteAbelianGroup {

  /*!
Cardinality of the group.
  */
  unsigned long long d_size; // this limits the order of representable groups
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
  FiniteAbelianGroup() {}

  explicit FiniteAbelianGroup(const GroupType&);

// copy and assignment: defaults are ok

// accessors
  bool operator== (const FiniteAbelianGroup& A) const {
    return d_type == A.d_type;
  }

  bool operator!= (const FiniteAbelianGroup& A) const {
    return not operator==(A);
  }

  unsigned int rank() const { return d_type.size(); }
  const GroupType& type() const {  return d_type; }

  // calling the following member |size| proved to be error prone
  unsigned long long order() const { return d_size; } // order of entire group

// representation conversions
  GrpArr toArray(GrpNbr x) const
    { GrpArr result(rank()); to_array(result,x,d_type); return result; }

  GrpArr toArray(const matrix::Vector<int>& v) const
    { GrpArr result(rank()); to_array(result,v,d_type); return result; }

  GrpNbr toGrpNbr(const GrpArr& a) const { return to_GrpNbr(a,d_type); }

  void toWeight(matrix::Vector<int>&, const GrpArr&) const;

  void toWeight(matrix::Vector<int>&, GrpNbr) const;

  // arithmetic on group elements

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

  GrpNbr subtract(GrpNbr, const GrpNbr) const;

  GrpArr& subtract(GrpArr&, const GrpArr&) const;

}; // |class FiniteAbelianGroup|

/*!
  A |Homomorphism| object represents a homomorphism from a group of type
  |d_source| to a group of type |d_dest|. The essential data is a matrix
  |d_matrix| with integral coefficients. As usual the matrix entry (i,j)
  expresses component $i$ of the image of generator $j$.

  For such a homomorphism to be well defined, the matrix entry (i,j) should be
  a multiple of |d_dest[i]/g|, with |g == gcd(d_dest[i],d_source[j])|.

  Computationally we shall proceed as follows. We set $M$ to a common multiple
  of the annihilators of |d_source| and |d_dest|, and think of all cyclic
  groups as embedded in $Z/M$; this means that before acting upon a source
  element $(x_1,...,x_m)$ we should multiply each component $x_j$ by the
  cotype $q_j=M/s_j$ of that component with respect to $M$, where $s_j$ is
  |d_source[j]|. After forming $y_j=\sum_j a_{i,j} q_j x_j$ for each $i$, the
  component $y_i$ must be in the subgroup of index $t_i$ in $Z/M$, where $t_i$
  is |d_dest[i]|, in other words it must be divisible by $M/t_i$, and the
  compoentn $i$ of the final result will be $y_i/(M/t_i)$. This definition is
  independent of the choice of $M$; we will take the lcm of the annihilators.

  In fact this definition might work for some group elements of the source but
  not for all. In that case we don't really have a homomorphism, but we allow
  applying to those group elements for which it works nontheless. The
  predicates |defined| serve to find out whether a source element is OK.
*/
class Homomorphism {

 private:

  GroupType d_source;
  GroupType d_dest;
  GroupType d_cosource; // cotypes of source components with respect to $M$
  GroupType d_codest;   // cotypes destination components with respect to $M$
  unsigned long d_annihilator; // the common multiple $M$
  matrix::Matrix<unsigned long> d_matrix;

 public:

// constructors and destructors
  Homomorphism(const std::vector<GrpArr>&,
	       const FiniteAbelianGroup&,
	       const FiniteAbelianGroup&);

// accessors
  GrpArr operator*(const GrpArr&) const;

  GrpNbr operator*(GrpNbr) const;
};

} // |namespace abelian|

} // |namespace atlas|

#endif
