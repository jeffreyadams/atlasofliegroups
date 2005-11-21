/*
  This is abelian.h
  
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

 protected:

  unsigned long d_size;
  GroupType d_type;
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
