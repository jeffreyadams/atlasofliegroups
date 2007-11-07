/*!
\file
  This is filekl.h
*/
/*
  Copyright (C) 2007 Marc van Leeuwen
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/


#ifndef FILEKL_H
#define FILEKL_H

#include <iostream>
#include <fstream>

#include "bitset_fwd.h"
#include "blocks.h"
#include "kl.h"

namespace atlas {

  namespace filekl {


const unsigned int magic_code=0x06ABdCF0; // indication of new matrix format


    // The |polynomial_info| class

typedef unsigned long long int ullong; // a 64-bit type even on 32-bit machines

 typedef ullong KLIndex; // NOTE: this differs from the type kl:KLIndex
 // here it must be |unsigned long long| to avoid problems when multiplied

// A class giving access to the polynomial file (which remains open)
class polynomial_info
{
  std::ifstream& file; // non-owned reference to (open) file

  KLIndex n_pols;         // number of polynomials in file
  unsigned int coef_size; // number of bytes per coefficient
  ullong n_coef;          // number of coefficients
  std::streamoff index_begin, coefficients_begin;

public:
  polynomial_info(std::ifstream& coefficient_file);
  virtual ~polynomial_info() { file.close(); }

  KLIndex n_polynomials() const { return n_pols; }
  unsigned int coefficient_size() const { return coef_size; }
  ullong n_coefficients() const { return n_coef; }

  virtual size_t degree(KLIndex i) const;
  std::vector<size_t> coefficients(KLIndex i) const;
  virtual size_t leading_coeff(KLIndex i) const;
  ullong coeff_start(KLIndex i) const; // number of all preceding coefficients
};

// A derived class that caches the degrees, and some leading coefficients
class cached_pol_info
  : public polynomial_info
{
  static const size_t degree_mask = 0x1F; // degree must <32
  mutable std::vector<unsigned char> cache;

public:
  cached_pol_info(std::ifstream& coefficient_file);

  virtual size_t degree(KLIndex i) const;
  virtual size_t leading_coeff(KLIndex i) const;
};

    // The |block_info| class

using blocks::BlockElt;
using bitset::RankFlags;

typedef std::vector<RankFlags> descent_set_vector; // indexed by block element

typedef std::vector<BlockElt> ascent_vector;       // indexed by simple root
typedef std::vector<ascent_vector> ascent_table;   // indexed by block element

typedef std::vector<BlockElt> prim_list;           // list of weak primitives
typedef std::vector<prim_list> prim_table;


struct block_info
{
  unsigned int rank;
  BlockElt size;
  unsigned int max_length; // maximal length of block elements
  std::vector<BlockElt> start_length;
   // array has size |max_length+2|; it defines intervals for each length

  descent_set_vector descent_set;     // descent (bit)set listed per BlockElt

private:
  ascent_table ascents;       // raised element per simple root, per BlockElt
  prim_table primitives_list; // lists of weakly primitives, per descent set

public:
  block_info(std::ifstream& in); // constructor reads, and closes, file

  BlockElt primitivize(BlockElt x, BlockElt y) const;
  const prim_list& prims_for_descents_of(BlockElt y);
private:
  bool is_primitive(BlockElt x, const RankFlags d) const;
};



    // the |matrix_info| class

typedef prim_list strong_prim_list;


class matrix_info
{
  std::ifstream& matrix_file; // non-owned reference to (open) file

  block_info block;

  std::vector<std::streampos> row_pos; // positions where each row starts

// data for currently selected row~|y|
  BlockElt cur_y;		       // row number
  strong_prim_list cur_strong_prims;   // strongly primitives for this row
  std::streampos cur_row_entries; // indices of polynomials for row start here

//private methods
  matrix_info(const matrix_info&); // copying forbidden
  void set_y(BlockElt y);          // install |cur_y| and dependent data

public:
  BlockElt x_prim; // public variable that is set by |find_pol_nr|

// constructor and destructor
  matrix_info(std::ifstream& block_file,std::ifstream& m_file);
  ~matrix_info() {}

// accessors
  size_t rank() const { return block.rank; }
  BlockElt block_size() const { return block.size; }
  size_t length (BlockElt y) const; // length in block
  BlockElt first_of_length (size_t l) const { return block.start_length[l]; }
  RankFlags descent_set (BlockElt y) const { return block.descent_set[y]; }
  std::streamoff row_offset(BlockElt y) const { return row_pos[y]; }
  BlockElt primitivize (BlockElt x,BlockElt y) const
    { return block.primitivize(x,y); }

// manipulators (they are so because they set |cur_y|)
  KLIndex find_pol_nr(BlockElt x,BlockElt y);
  BlockElt prim_nr(unsigned int i,BlockElt y);// find primitive element
  const strong_prim_list& strongly_primitives (BlockElt y)
    { set_y(y); return cur_strong_prims; } // changing |y| invalidates this!
};



    // the progress_info class

class progress_info
{
  std::vector<KLIndex> first_pol; // count distinct polynomials in rows before
public:
  progress_info(std::ifstream& progress_file);

  BlockElt block_size() const { return first_pol.size()-1; }
  KLIndex first_new_in_row(BlockElt y) const // |y==block_size()| is allowed
    { return first_pol[y]; }
  BlockElt first_row_for_pol(KLIndex i) const;
};


//             F u n c t i o n s

void write_block_file(const blocks::Block& block, std::ostream& out);
void write_matrix_file(const kl::KLContext& klc, std::ostream& out);
void write_KL_store(const kl::KLStore& store, std::ostream& out);

} // namespace filekl

} // namespace atlas

#endif
