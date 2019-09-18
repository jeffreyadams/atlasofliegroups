/*
  This is block_minimal.h

  Copyright (C) 2019 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Variant of |param_block| in blocks.h, to be shared when having same KL data

#ifndef BLOCK_MINIMAL_H  /* guard against multiple inclusions */
#define BLOCK_MINIMAL_H

#include "blocks.h" // we conceptually just extend that module
#include "subsystem.h"

namespace atlas {

namespace blocks {

// a class for blocks of (possibly non integral) parameters
class block_minimal : public Block_base
{
  const Rep_context& rc; // accesses many things, including KGB set for x
  const SubSystem integral_datum;

  y_entry::Pooltype y_pool;
  y_part_hash y_hash;  // hash table allows storing |y| parts by index
  std::vector<TorusElement> y_part; // as in the y_values module, indexed by |y|

  // hash structure to allow rapid lookup of |(x,y)| index pairs
  block_hash xy_hash;

  // group small components together:
  KGBElt highest_x,highest_y; // maxima over this (maybe partial) block

 public:

  // constructor
  block_minimal
    (const repr::Rep_context& rc,
     StandardRepr sr, // by value,since it will be made dominant before use
     BlockElt& entry_element	// set to block element matching input
    );

 public:
  // accessors that get values via |rc|
  const repr::Rep_context& context() const { return rc; }
  const RootDatum& rootDatum() const;
  const InnerClass& innerClass() const;
  const InvolutionTable& involution_table() const;
  RealReductiveGroup& realGroup() const;

  // with |gamma| unknown, only the difference |gamma-lambda| is meaningful
  RatWeight gamma_lambda(BlockElt z) const;

  BlockElt lookup(const StandardRepr& sr) const;

  ext_gens fold_orbits(const WeightInvolution& delta) const;

  // virtual methods
  virtual KGBElt max_x() const { return highest_x; } // might not be final |x|
  virtual KGBElt max_y() const { return highest_y; }

  virtual std::ostream& print // defined in block_io.cpp
    (std::ostream& strm, BlockElt z,bool as_invol_expr) const;


 private:
  void compute_duals();

/*
  reverse lengths and order block with them increasing, and by increasing
  |x(z)| among elements of given length; adapt tables accordingly.
*/
  void reverse_length_and_sort();

}; // |class block_minimal|



} // |namespace blocks|

} // |namespace atlas|
#endif
