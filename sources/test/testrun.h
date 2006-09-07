/*
  This is testrun.h
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef TESTRUN_H  /* guard against multiple inclusions */
#define TESTRUN_H

#include "abelian.h"
#include "lietype.h"
#include "prerootdata.h"
#include "realredgp.h"

namespace atlas {

/******** type declarations **************************************************/

namespace testrun {

  enum Category {Simple, Semisimple, Reductive, Adjoint, SimplyConnected,
  Split, Complex};

  typedef std::vector<unsigned long> Shape;

  class GroupIterator;

  class LieTypeIterator;
  class CoveringIterator;
  class RealFormIterator;
  class SubGroupIterator;
  class TorusPartIterator;

}

/******** function declarations **********************************************/

namespace testrun {

}

/******** type definitions ***************************************************/

namespace testrun {

// base class for group iterators

class GroupIterator {

 public:
// associated types
  typedef std::forward_iterator_tag iterator_category;
  typedef realredgp::RealReductiveGroup value_type;
  typedef ptrdiff_t difference_type;
  typedef const value_type* pointer;
  typedef const value_type& reference;

// constructors and destructors
  GroupIterator() {}

  virtual ~GroupIterator() {}

// accessors
  virtual bool operator== (const GroupIterator& i) const = 0;

  bool operator!= (const GroupIterator& i) const {
    return not operator== (i);
  }

  virtual reference operator* () const = 0;

  virtual pointer operator-> () const = 0;

  virtual bool operator() () const = 0;

// manipulators
  virtual GroupIterator& operator++ () = 0;
// note : operator++ (int) cannot be defined here because it must return
// an object.
};

class LieTypeIterator {

 private:

  lietype::LieType d_type;
  size_t d_firstRank;
  size_t d_lastRank;
  Category d_category;
  bool d_done;

 public:

// associated types
  typedef std::forward_iterator_tag iterator_category;
  typedef lietype::LieType value_type;
  typedef ptrdiff_t difference_type;
  typedef const value_type* pointer;
  typedef const value_type& reference;

// constructors and destructors
  LieTypeIterator() {}

  LieTypeIterator(Category, size_t);

  virtual ~LieTypeIterator() {}

// accessors
  bool operator== (const LieTypeIterator& i) const {
    return d_type == i.d_type;
  }

  bool operator!= (const LieTypeIterator& i) const {
    return not operator== (i);
  }

  reference operator* () const {
    return d_type;
  }

  pointer operator-> () const {
    return &d_type;
  }

  bool operator() () const {
    return not d_done;
  }

// manipulators
  LieTypeIterator& operator++ ();

  LieTypeIterator operator++ (int) {
    LieTypeIterator tmp(*this); 
    ++(*this); 
    return tmp;
  }
};

class TorusPartIterator {

 private:

  size_t d_rank;
  std::vector<bitmap::BitMap::iterator> d_data;
  std::vector<abelian::GrpNbr> d_returnValue;
  bitmap::BitMap::iterator d_first;
  bitmap::BitMap::iterator d_last;
  bool d_done;

 public:

// associated types
  typedef std::forward_iterator_tag iterator_category;
  typedef std::vector<abelian::GrpNbr> value_type;
  typedef ptrdiff_t difference_type;
  typedef const value_type* pointer;
  typedef const value_type& reference;

// constructors and destructors
  TorusPartIterator() {}

  TorusPartIterator(size_t, const bitmap::BitMap&);

  TorusPartIterator(size_t, const std::vector<bitmap::BitMap::iterator>&,
		    const bitmap::BitMap::iterator&, 
		    const bitmap::BitMap::iterator&);

  ~TorusPartIterator() {}

// accessors
  bool operator== (const TorusPartIterator& i) const {
    return d_data == i.d_data;
  }

  bool operator!= (const TorusPartIterator& i) const {
    return not operator== (i);
  }

  reference operator* () const {
    return d_returnValue;
  }

  pointer operator-> () const {
    return &d_returnValue;
  }

  bool operator() () const {
    return not d_done;
  }

  size_t rank() const {
    return d_rank;
  }

// manipulators
  TorusPartIterator& operator++ ();

  TorusPartIterator operator++ (int) {
    TorusPartIterator tmp(*this); 
    ++(*this); 
    return tmp;
  }

  void reset(const bitmap::BitMap&);
};

class SubgroupIterator {

 private:

  abelian::FiniteAbelianGroup* d_group;
  std::vector<bitmap::BitMap> d_prevRank;
  std::set<bitmap::BitMap> d_thisRank;
  std::vector<bitmap::BitMap>::const_iterator d_prev;
  bitmap::BitMap d_subgroup;
  bitmap::BitMap d_cycGenerators;
  bitmap::BitMap::iterator d_generator;
  unsigned long d_rank;
  bool d_done;

// private member functions
  void incrementGenerator();
  void resetGenerator();

 public:

// associated types
  typedef std::forward_iterator_tag iterator_category;
  typedef bitmap::BitMap value_type;
  typedef ptrdiff_t difference_type;
  typedef const value_type* pointer;
  typedef const value_type& reference;

// constructors and destructors
  SubgroupIterator() {};

  explicit SubgroupIterator(abelian::FiniteAbelianGroup& A);

  ~SubgroupIterator() {};

// copy and assignment
  SubgroupIterator(const SubgroupIterator&);

  SubgroupIterator& operator= (const SubgroupIterator&);

// accessors
  bool operator== (const SubgroupIterator& i) const {
    return d_subgroup == i.d_subgroup;
  }

  bool operator!= (const SubgroupIterator& i) const {
    return not operator== (i);
  }

  reference operator* () const {
    return d_subgroup;
  }

  bool operator() () const {
    return not d_done;
  }

  const abelian::FiniteAbelianGroup& group() const {
    return *d_group;
  }

  unsigned long rank() const {
    return d_rank;
  }

// manipulators
  SubgroupIterator& operator++ ();

  SubgroupIterator operator++ (int) {
    SubgroupIterator tmp(*this); 
    ++(*this); 
    return tmp;
  }
};

class CoveringIterator {

 private:

  lietype::LieType d_lieType;
  abelian::FiniteAbelianGroup* d_dcenter;
  size_t d_rank;
  size_t d_torusRank;
  size_t d_semisimpleRank;

// iterator management
  bitmap::BitMap d_quotReps;
  SubgroupIterator d_subgroup;
  TorusPartIterator d_torusPart;
  bool d_done;

// data for the PreRootDatum
  latticetypes::WeightList d_smithBasis;
  prerootdata::PreRootDatum d_preRootDatum;

// private member functions
  void adjustBasis(latticetypes::WeightList&, const latticetypes::WeightList&);
  void makeBasis(latticetypes::WeightList&);
  CoveringIterator(const CoveringIterator&);
  CoveringIterator& operator= (const CoveringIterator&);

 public:

// associated types
  typedef std::forward_iterator_tag iterator_category;
  typedef prerootdata::PreRootDatum value_type;
  typedef ptrdiff_t difference_type;
  typedef const value_type* pointer;
  typedef const value_type& reference;

// constructors and destructors
  CoveringIterator() {}

  explicit CoveringIterator(const lietype::LieType&);

  virtual ~CoveringIterator();

// accessors
  bool operator== (const CoveringIterator& i) const {
    return d_preRootDatum == i.d_preRootDatum;
  }

  bool operator!= (const CoveringIterator& i) const {
    return not operator== (i);
  }

  reference operator* () const {
    return d_preRootDatum;
  }

  pointer operator-> () const {
    return &d_preRootDatum;
  }

  bool operator() () const {
    return not d_done;
  }

  const abelian::FiniteAbelianGroup& dcenter() {
    return *d_dcenter;
  }

  const bitmap::BitMap& group() const {
    return *d_subgroup;
  }

// manipulators
  CoveringIterator& operator++ ();
  CoveringIterator operator++ (int) {
    CoveringIterator tmp(*this); 
    ++(*this); 
    return tmp;
  }
};

class RealFormIterator {
};

}

}

#endif
