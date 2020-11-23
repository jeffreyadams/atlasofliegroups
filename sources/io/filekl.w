% Copyright (C) 2007,2008 Marc van Leeuwen
% This file is part of the Atlas of Reductive Lie Groups software (the Atlas)

% This program is made available under the terms stated in the GNU
% General Public License (GPL), see http://www.gnu.org/licences/licence.html

% The Atlas is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

% The Atlas is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with the Atlas; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

@* Introduction.
%
This source file groups functionality for exchanging information between the
\.{Fokko} program and external utility programs that provide information about
Kazhdan-Lusztig polynomials, via some binary files that contain the results of
Fokko's computations. The goal is that such results can be accessed without
having to run \.{Fokko} again, which is quite essential in the case of the
Kazhdan-Lusztig polynomials for the big block of~$E_8$, but which can be
useful for other cases as well.

A variety of simple binary file formats is defined, with in each case an
accompanying class whose constructor will open a file of that type and whose
methods provide the necessary interpretation of the binary format, so that the
utility programs can be written in  a natural style almost as if they had
access to the \.{Fokko} structures that produced the values written to file.

The definitions of those classes is easiest to understand if juxtaposed to the
definitions of the routines that write the binary files. It is important
however to have separate object files for the two kinds of functions, since
the writing functions need to have access to many internal Atlas classes,
whereas the reading functions should not depend on linking together with those
Atlas classes, lest the utility programs become as large as \.{Fokko}
itself. This source file then will write a file \.{filekl.cpp} that compiles
to an object module of \.{Fokko}, and  another file \.{filekl\_in.cpp} whose
object module is to be incorporated into utility programs. Each of modules also
has a corresponding header file, so all in all we produce 4 files for
compilation from the current source text. Our main output file
is \.{filekl.cpp}, with others explicitly added below.

@h "filekl.h"

@c
namespace atlas {
  namespace filekl {
    @< Functions for writing binary files @>@;
  }@;
}@;

@ In \.{filekl.h} we declare everything that we want to export
from \.{filekl.cpp}, and define constant values that should be shared between
that implementation file and its uses. The constant definitions will also be
included in the header for the \.{filekl\_in} module, so that they are ensured
to be the same there.

@( filekl.h @>=

#ifndef FILEKL_H
#define FILEKL_H

@< Includes needed in the header file @>@;
namespace atlas {
  namespace filekl {
    @< Constants common for writing and reading @>@;
    @< Declarations of exported functions @>@;
  }@;
}@;
#endif

@ The header file \.{filekl.h} requires no include files other than the
forward declarations contained in the \.{iosfwd} standard header
and \.{Atlas.h}.

@< Includes needed in the header file @>=
#include <iosfwd>

#include "../Atlas.h"

@ The \.{filekl\_in} implementation does a lot of file reading, and notably uses
the |basic_io::read_bytes| function template to read a fixed number of
consecutive bytes into an |unsigned long long int| value.

@( filekl_in.cpp @>=

#include "filekl_in.h"
#include <stdexcept>

#include "basic_io.h"
@< Includes needed in the input implementation file @>

namespace atlas {
  namespace filekl {

    using basic_io::read_bytes;

    @< Constants common for writing and reading @>@;

    @< Methods for reading binary files @>@;
  }@;
}@;

@ Note that the type |KLIndex| defined below this differs from |kl::KLIndex|.
Here the type must be |unsigned long long| to avoid problems when multiplied.

@< Type definitions for reading files @>=
typedef unsigned long long int ullong; // a 64-bit type even on 32-bit machines

typedef ullong KLIndex;


@
@( filekl_in.h @>=

#ifndef FILEKL_IN_H
#define FILEKL_IN_H


@< Includes needed in the input header file @>@;
namespace atlas {
  namespace filekl {
    @< Type definitions for reading files @>@;
    @< Input class declarations @>@;
  }@;
}@;
#endif

@ In the file \.{filekl\_in.h} we need not only forward class declarations
from \.{iosfwd}, but also typedef symbols like |streampos|, so we include the
more elaborate standard header \.{ios}.

@< Includes needed in the input header file @>=
#include <ios>
@)
#include "bitset.h" // to make |RankFlags| a complete type; used when inlining
#include "../Atlas.h"

@ We need two testable out-of-range values |UndefBlock| and |no_good_ascent| to
record special circumstances in places where |BlockElt| value is expected.

To test whether the special matrix storage file was one that was created by us,
we use a special $32$ bit value |magic_code| that should be present in a
specific place. It's value was chosen in memory of Fokko du Cloux.

@< Constants common for writing and reading @>=

const BlockElt no_good_ascent = UndefBlock-1;
 // value flagging that no good ascent exists
const unsigned int magic_code=0x06ABdCF0; // indication of new matrix format

@* Writing a block file.
Here is how a block file is written.

First the block size in $4$ bytes, then its rank in $1$ byte, the maximal length
for a block element in $1$ byte, then for each non negative value $l$ strictly
less that this maximal length the number of elements of length${}\leq l$, each
in $4$ bytes (this additional information will allow to deduce its length from
any used |BlockElt| value), then for every block element its descent set in $4$
bytes,

@< Declarations of exported functions @>=
void write_block_file(const Block& block, std::ostream& out);

@~@h "blocks.h"
@h "basic_io.h"
@< Functions for writing binary files @>=

void write_block_file(const Block& block, std::ostream& out)
{
  unsigned char rank=block.rank(); // certainly fits in a byte

  basic_io::put_int(block.size(),out);  // block size in 4 bytes
  out.put(rank);                        // rank in 1 byte

  { // output length data
    unsigned char max_length=block.length(block.size()-1);
    out.put(max_length);

    BlockElt z=0;
    for (size_t l=0; l<max_length; ++l)
    {
      while(block.length(z)<=l)
	++z;
      basic_io::put_int(z,out);
      // record that there are |z| elements of |length<=l|
    }
    assert(z<block.size());
    // and don't write |basic_io::put_int(block.size(),out);|
  }

  // write descent sets
  for (BlockElt y=0; y<block.size(); ++y)
  {
    RankFlags d;
    for (size_t s = 0; s < rank; ++s)
      d.set(s,block.isWeakDescent(s,y));
    basic_io::put_int(d.to_ulong(),out); // write |d| as 32-bits value
  }

  // write table of primitivatisation successors
  for (BlockElt x=0; x<block.size(); ++x)
  {
#if VERBOSE
    std::cerr << x << '\r';
#endif
    for (size_t s = 0; s < rank; ++s)
    {
      DescentStatus::Value v = block.descentValue(s,x);
      if (DescentStatus::isDescent(v)
	  or v==DescentStatus::ImaginaryTypeII)
	basic_io::put_int(no_good_ascent,out);
      else if (v == DescentStatus::RealNonparity)
	basic_io::put_int(UndefBlock,out);
      else if (v == DescentStatus::ComplexAscent)
	basic_io::put_int(block.cross(s,x),out);
      else if (v == DescentStatus::ImaginaryTypeI)
	basic_io::put_int(block.cayley(s,x).first,out);
      else assert(false);
    }
  }
}

@* The {\bf block\_info} class.

@< Input class declarations @>=

typedef std::vector<RankFlags>
  descent_set_vector; // indexed by block element

typedef std::vector<BlockElt> ascent_vector; // indexed by simple root
typedef std::vector<ascent_vector> ascent_table;   // indexed by block element

typedef std::vector<BlockElt> prim_list;   // list of weak primitives
typedef std::vector<prim_list> prim_table;


struct block_info
{
  unsigned int rank;
  BlockElt size;
  unsigned int max_length; // maximal length of block elements
  std::vector<BlockElt> start_length;
  @/// array has size |max_length+2|; it defines intervals for each length

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


@*1 Methods of the {\bf block\_info} class.

@< Methods for reading binary files @>=

BlockElt
block_info::primitivize(BlockElt x, BlockElt y) const
{
  RankFlags d=descent_set[y];
start:
  if (x>=y) return x; // possibly with |x==UndefBlock|
  const ascent_vector& ax=ascents[x];
  for (size_t s=0; s<rank; ++s)
    if (d[s] and ax[s]!=no_good_ascent)
    @/{@; x=ax[s]; goto start; } // this should raise $x$, now try another step
  return x; // no raising possible, stop here
}

@
@< Methods for reading binary files @>=

bool
block_info::is_primitive(BlockElt x, const RankFlags d) const
{
  const ascent_vector& ax=ascents[x];
  for (size_t s=0; s<ax.size(); ++s)
    if (d[s] and ax[s]!=no_good_ascent)
      return false;
  return true;
  // now |d[s]| implies |ascents[s]==no_good_ascent| for all simple roots |s|
}

@
@< Methods for reading binary files @>=

const prim_list& block_info::prims_for_descents_of(BlockElt y)
{ RankFlags d=descent_set[y];
  unsigned long s=d.to_ulong();
  prim_list& result=primitives_list[s];
  if (result.empty())
  { for (BlockElt x=0; x<size; ++x)
      if (is_primitive(x,d))
	result.push_back(x);
    prim_list(result).swap(result); // reallocate to fit snugly
  }
  return result;
}

@
@< Includes needed in the input implementation file @>=
#include <fstream>
@~@< Methods for reading binary files @>=

block_info::block_info(std::ifstream& in)
  : rank(), size(), max_length(), start_length()
  , descent_set(), ascents(), primitives_list() // don't initialize yet
{
  in.seekg(0,std::ios_base::beg); // ensure we are reading from the start
  size=read_bytes<4>(in);
  rank=read_bytes<1>(in);
  max_length=read_bytes<1>(in);

  // read intervals of block elements for each length
  start_length.resize(max_length+2);
  start_length[0]=BlockElt(0);
  for (size_t i=1; i<=max_length; ++i) start_length[i]=read_bytes<4>(in);
  start_length[max_length+1]=size;

  // read descent sets
  descent_set.reserve(size);
  for (BlockElt y=0; y<size; ++y)
    descent_set.push_back(RankFlags(read_bytes<4>(in)));

  // read ascent table
  ascents.reserve(size);
  for (BlockElt x=0; x<size; ++x)
    {
      ascents.push_back(ascent_vector());
      ascent_vector& a=ascents.back();
      a.reserve(rank);
      for (size_t s=0; s<rank; ++s)
	a.push_back(read_bytes<4>(in));
    }

  // lists of primitives lazily computed, on demand by |prims_for_descents_of|
  primitives_list.resize(1ul<<rank); // create $2^{rank}$ empty vectors
}


@* Writing a matrix file.
Here is how a matrix file is written.

@h "kl.h"
@h <iostream> // to have |std::streamoff| defined

@< Functions for writing binary files @>=

std::streamoff
write_KL_row(const kl::KL_table& kl_tab, BlockElt y, std::ostream& out)
{
  BitMap prims=kl_tab.prim_map(y); // marks nonzero KL polys among primitives
  const auto& kld=kl_tab.KL_data(y);

  assert(kld.size()+1==prims.capacity()); // check the number of KL polynomials

  // write row number for consistency check on reading
  basic_io::put_int(y,out);

  std::streamoff start_row=out.tellp();

  // write number of primitive elements, plus 1 for |y| itself, for convenience
  basic_io::put_int(prims.capacity(),out);

  // now write the bitmap as a sequence of unsigned int values
  for (size_t i=0; i<prims.capacity(); i+=32)
    basic_io::put_int(prims.range(i,32),out);

  // finally, write the indices of the KL polynomials themselves
  for (size_t i=0; i<kld.size(); ++i)
  {
    assert((kld[i]!=0)==prims.isMember(i));
    if (kld[i]!=0) // only write nonzero indices
      basic_io::put_int(kld[i],out);
  }

  basic_io::put_int(1,out); // write unrecorded final polynomial 1

  // and signal if there was unsufficient space to write the row
  if (not out.good()) throw error::OutputError();

  return start_row;
}

@

@< Declarations of exported functions @>=
void write_matrix_file(const kl::KL_table& kl_tab, std::ostream& out);

@~@< Functions for writing binary files @>=

void write_matrix_file(const kl::KL_table& kl_tab, std::ostream& out)
{
  std::vector<unsigned int> delta(kl_tab.size());
  std::streamoff offset=0;
  for (BlockElt y=0; y<kl_tab.size(); ++y)
  {
    std::streamoff new_offset=write_KL_row(kl_tab,y,out);
    delta[y]=static_cast<unsigned int>((new_offset-offset)/4);
    offset=new_offset;
  }

  // now write the values allowing rapid location of the matrix rows
  for (BlockElt y=0; y<kl_tab.size(); ++y)
    basic_io::put_int(delta[y],out);

  // and finally sign file as being in new format by overwriting 4 bytes
  out.seekp(0,std::ios_base::beg);
  basic_io::put_int(magic_code,out);
}

@*The {\bf matrix\_info} class.

@< Includes needed in the input implementation file @>=
#include "blocks.h"
#include "bitset.h"
@~
@< Input class declarations @>=
typedef prim_list strong_prim_list;

class matrix_info
{
  std::ifstream& matrix_file; // non-owned reference to (open) file

  block_info block;

  std::vector<std::streampos> row_pos; // positions where each row starts

// data for currently selected row~|y|
  BlockElt cur_y;		// row number
  strong_prim_list cur_strong_prims;   // strongly primitives for this row
  std::streampos cur_row_entries; // indices of polynomials for row start here

//private methods
  matrix_info(const matrix_info&); // copying forbidden
  void set_y(BlockElt y);  // install |cur_y| and dependent data

public:
  BlockElt x_prim; // public variable that is set by |find_pol_nr|

// constructor and destructor
  matrix_info(std::ifstream& block_file,std::ifstream& m_file);
  ~matrix_info() {}

// accessors
  size_t rank() const @+{@; return block.rank; }
  BlockElt block_size() const @+{@; return block.size; }
  size_t length (BlockElt y) const; // length in block
  BlockElt first_of_length (size_t l) const
    @+{@; return block.start_length[l]; }
  RankFlags descent_set (BlockElt y) const
    @+{@; return block.descent_set[y]; }
  std::streamoff row_offset(BlockElt y) const @+{@; return row_pos[y]; }
  BlockElt primitivize (BlockElt x,BlockElt y) const
    @+{@; return block.primitivize(x,y); }
@/
// manipulators (they are so because they set |cur_y|)
  KLIndex find_pol_nr(BlockElt x,BlockElt y);
  BlockElt prim_nr(unsigned int i,BlockElt y);
    // find primitive element
  const strong_prim_list& strongly_primitives (BlockElt y)
    {@; set_y(y); return cur_strong_prims; } // changing |y| invalidates this!
};

@*1 Methods of the {\bf matrix\_info} class.


@< Methods for reading binary files @>=

size_t matrix_info::length (BlockElt y) const
{ return
    std::upper_bound(block.start_length.begin(),block.start_length.end(),y)
    -block.start_length.begin() // index of first element of length |l(y)+1|
    -1; // now we have just |l(y)|
}

@
@< Methods for reading binary files @>=

void matrix_info::set_y(BlockElt y)
{
  if (y==cur_y)
    {@; matrix_file.seekg(cur_row_entries); return; }
  cur_y=y;
  const prim_list& weak_prims = block.prims_for_descents_of(y);
  cur_strong_prims.resize(0);
  // restart building from scratch, but don't deallocate storage

  matrix_file.seekg(row_pos[y]);
  size_t n_prim=read_bytes<4>(matrix_file);


  for (size_t i=0; i<n_prim; i+=32)
    {
      unsigned int chunk=read_bytes<4>(matrix_file);
      for (size_t j=0; chunk!=0; ++j,chunk>>=1) // and certainly |j<32|
        if ((chunk&1)!=0) cur_strong_prims.push_back(weak_prims[i+j]);
    }

  {
#ifndef NDEBUG
    const BlockElt* i=
      // point after first block element of length of |y|
      std::upper_bound(&block.start_length[0]
  		    ,&block.start_length[block.max_length+1]
  		    ,y);
    const BlockElt* first=
      // point to first of |weak_prims| of that length
      std::lower_bound(&weak_prims[0],&weak_prims[weak_prims.size()],*(i-1));
    if (cur_strong_prims.back()!=*first)
      {
	std::cerr << "For y=" << y << ", first weak prim found " << *first
		  << " does not match computed value "
		  << cur_strong_prims.back() << ".\n";
	throw std::runtime_error("Panic");
      }
#endif
    cur_strong_prims.back()=y; // replace by |y|
  }

  cur_row_entries=matrix_file.tellg();
}

@
@< Methods for reading binary files @>=

KLIndex matrix_info::find_pol_nr(BlockElt x,BlockElt y)
{
  set_y(y);
  x_prim=block.primitivize(x,y);
  if (x_prim>=y)
    return KLIndex(x_prim==y ? 1 : 0); // primitivisation copped out
  strong_prim_list::const_iterator it=
    std::lower_bound(cur_strong_prims.begin(),cur_strong_prims.end(),x_prim);
  if (it==cur_strong_prims.end() or *it!=x_prim)
    return KLIndex(0); // not strong

  matrix_file.seekg(4*size_t(it-cur_strong_prims.begin()),std::ios_base::cur);
  return KLIndex(read_bytes<4>(matrix_file));
}


@
@< Methods for reading binary files @>=

BlockElt matrix_info::prim_nr(unsigned int i,BlockElt y)
{ const prim_list& weak_prims = block.prims_for_descents_of(y);
  const BlockElt* it=
    // point after first block element of length of |y|
    std::upper_bound(&block.start_length[0]
		    ,&block.start_length[block.max_length+1]
		    ,y);
  const BlockElt* first=
    // point to first of |weak_prims| of that length
    std::lower_bound(&weak_prims[0],&weak_prims[weak_prims.size()],*(it-1));
  size_t limit=first-&weak_prims[0]; // limiting value for |i|
  if (i<limit) return weak_prims[i];
  else if (i==limit) return y;
  std::cerr << "Limit is " << limit << ".\n";
  throw std::runtime_error("Weakly primitive element index too large");
}

@
@< Methods for reading binary files @>=

matrix_info::matrix_info
  (std::ifstream& block_file,std::ifstream& m_file)
: matrix_file(m_file) // store reference to the matrix file
  , block(block_file) // read in block information
  , row_pos(block.size) // dimension these vectors
  , cur_y(UndefBlock), cur_strong_prims(), cur_row_entries(0)
{
  matrix_file.seekg(0,std::ios_base::beg);
  if (read_bytes<4>(matrix_file)==magic_code)
  { matrix_file.seekg(-4*std::streamoff(block_size()),std::ios_base::end);
    std::streamoff cumul=0;
    for (BlockElt y=0; y<block_size(); ++y)
    { cumul+= 4*std::streamoff(read_bytes<4>(matrix_file));
      row_pos[y] = cumul;
    }
  }

  else // read old file format, scanning the whole file
  {
    size_t l=0; // length of y
    matrix_file.seekg(0,std::ios_base::beg);
    for (BlockElt y=0; y<block.size; ++y)
    {
      while (y>=block.start_length[l+1]) ++l;

      if (read_bytes<4>(matrix_file)!=y and y!=0)
      { std::cerr << y << std::endl;
        throw std::runtime_error ("Alignment problem");
      }
      row_pos[y]= matrix_file.tellg(); // record position after row number
      size_t n_prim=read_bytes<4>(matrix_file);

      { size_t n_strong_prim=0;
        { // compute number of entries to skip, while reading bitmap
	  static const unsigned int ulsize=sizeof(unsigned long int);
	  size_t count=n_prim/(8*ulsize); // number of |unsigned long|s to read

	  while (count-->0)
	    n_strong_prim+=bits::bitCount(read_bytes<ulsize>(matrix_file));

	  count=(n_prim%(8*ulsize)+31)/32; // maybe some tetra-byte(s) left
	  while (count-->0)
	    n_strong_prim+=bits::bitCount(read_bytes<4>(matrix_file));
	}

	// now skip over matrix entries
	matrix_file.seekg(4*std::streamoff(n_strong_prim)
			 ,std::ios_base::cur);
      }
      if (not matrix_file.good())
      { std::cerr << y << std::endl;
	throw std::runtime_error ("Premature end of file");
      }
    } // for (BlockElt y...)
  } // |if (...==magic_code)|

  block_file.close(); // success, we no longer need the block file
}


@* Writing a polynomial file.
Here is how the polynomial file is written

@< Declarations of exported functions @>=
void write_KL_store(const kl::KLStore& store, std::ostream& out);

@~This routine prefers a simple format over an extremely space-optimised
representation on disk. After writing the number |N| of polynomials in
$4$~bytes, we write a sequence of |N+1| indices of 5 bytes each, giving for
each polynomial |i| the position of its first (degree 0) coefficient in the
global list and as final 5-byte value (number |N|) the total number of
coefficients. After that, starting from position |9+5*N|, the list of all
coefficients, starting for each polynomial with the constant coefficient and
up to the leading coefficient. The degree of polynomial~$i$ is implicit in the
value of indices $i$~and~$i+1$: their difference, divided by the coefficient
size, is the number of coefficients, the degree is one less.

@< Functions for writing binary files @>=

void write_KL_store(const kl::KLStore& store, std::ostream& out)
{
  const size_t coef_size=4; // dictated (for now) by |basic_io::put_int|

  basic_io::put_int(store.size(),out); // write number of KL poynomials

  // write sequence of 5-byte indices, computed on the fly
  size_t offset=0; // position of first coefficient written
  for (size_t i=0; i<store.size(); ++i)
  {
    kl::KLPolRef p=store[i]; // get reference to polynomial

    // output 5-byte value of offset
    basic_io::put_int(offset&0xFFFFFFFF,out);
    out.put(char(offset>>16>>16)); // >>32 would fail on 32 bits machines

    if (not p.isZero()) // superfluous since polynomials::MinusOne+1==0
      // add number of coefficient bytes to be written
      offset += (p.degree()+1)*coef_size;
  }
  // write final 5-byte value (total size of coefficiant list)
  basic_io::put_int(offset&0xFFFFFFFF,out); out.put(char(offset>>16>>16));

  // now write out coefficients
  for (size_t i=0; i<store.size(); ++i)
  {
    kl::KLPolRef p=store[i]; // get reference to polynomial
    if (not p.isZero())
      for (size_t j=0; j<=p.degree(); ++j)
	basic_io::put_int(p[j],out);
  }
}

@* The {\bf polynomial\_info} class.
%
The class |polynomial_info| gives access to polynomials stored in a file, using
random file access to treat the storage as an indexable repository of
polynomials. The method |coefficients| produces the polynomial at a given index,
and the virtual methods |degree| and |leading_coefficient| give less complete
information, but which is expected to be most frequently accessed. A derived
class might want to devote some memory to speeding up those accesses.

@< Input class declarations @>=
class polynomial_info
{
  std::ifstream& file; // non-owned reference to (open) file

  KLIndex n_pols;         // number of polynomials in file
  unsigned int coef_size; // number of bytes per coefficient
  ullong n_coef;          // number of coefficients
  std::streamoff index_begin, coefficients_begin;

public:
  polynomial_info(std::ifstream& coefficient_file);
  virtual ~polynomial_info();
@)
  KLIndex n_polynomials() const @+{@; return n_pols; }
  unsigned int coefficient_size() const @+{@; return coef_size; }
  ullong n_coefficients() const @+{@; return n_coef; }
@)
  virtual size_t degree(KLIndex i) const;
  std::vector<size_t> coefficients(KLIndex i) const;
  virtual size_t leading_coeff(KLIndex i) const;
};

@*1 Methods of the {\bf polynomial\_info} class.

@< Methods for reading binary files @>=

polynomial_info::polynomial_info (std::ifstream& coefficient_file)
: file(coefficient_file), n_pols(read_bytes<4>(file))
, coef_size(), index_begin(file.tellg()), coefficients_begin()
{ file.seekg(10,std::ios_base::cur); // skip initial 2 indices
  coef_size=read_bytes<5>(file); // size of the |One|
  file.seekg(index_begin+5*n_pols,std::ios_base::beg);
  n_coef=read_bytes<5>(file)/coef_size;
  coefficients_begin=file.tellg();
}

polynomial_info::~polynomial_info() @+{@; file.close(); }

@
@< Methods for reading binary files @>=

size_t polynomial_info::degree(KLIndex i) const
{ if (i<2)
    return i-1; // exit for Zero and One
  file.seekg(index_begin+5*i,std::ios_base::beg);
  ullong index=read_bytes<5>(file);
  ullong next_index=read_bytes<5>(file);
  size_t length=(next_index-index)/coef_size;
  return length-1;
}

@
@< Methods for reading binary files @>=

std::vector<size_t> polynomial_info::coefficients(KLIndex i) const
{ file.seekg(index_begin+5*i,std::ios_base::beg);
  ullong index=read_bytes<5>(file);
  ullong next_index=read_bytes<5>(file);
  size_t length=(next_index-index)/coef_size;

  std::vector<size_t> result(length);
  file.seekg(coefficients_begin+index,std::ios_base::beg);

  for (size_t i=0; i<length; ++i)
    result[i]=basic_io::read_var_bytes(coef_size,file);

  return result;
}

@ The leading coefficient is the last one stored before the next polynomial.

@< Methods for reading binary files @>=

size_t polynomial_info::leading_coeff(KLIndex i) const
{ if (i<2) return i; // this makes "leading coefficient" of Zero return 0
  file.seekg(index_begin+5*(i+1),std::ios_base::beg);
  ullong next_index=read_bytes<5>(file);
  file.seekg(coefficients_begin+next_index-coef_size,std::ios_base::beg);
  return basic_io::read_var_bytes(coef_size,file);
}

@* The {\bf cached\_pol\_info} class.
This is a derived class that caches the degrees, and some leading coefficients.
Since the |const| virtual method |leading_coeff| will use and possibly update the
|cache| member, that member needs to be declared |mutable|.

@< Input class declarations @>=

class cached_pol_info
  : public polynomial_info
{
  static const size_t degree_mask = 0x1F; // degree must be${}<32$
  mutable std::vector<unsigned char> cache;

public:
  cached_pol_info(std::ifstream& coefficient_file);

  virtual size_t degree(KLIndex i) const;
  virtual size_t leading_coeff(KLIndex i) const;
};

@*1 Methods of the {\bf cached\_pol\_info} class.
@< Methods for reading binary files @>=

cached_pol_info::cached_pol_info(std::ifstream& coefficient_file)
  : polynomial_info(coefficient_file)
  , cache(0)
{
  cache.reserve(n_polynomials()-2);
  for (KLIndex i=2; i<n_polynomials(); ++i) // skip first two polynomials: $0,1$
  {
#ifdef VERBOSE
    if ((i&0xFFF)==0) std::cerr << i << '\r';
      // report progress every once in a while
#endif
    auto d=polynomial_info::degree(i);
    if ((d&~degree_mask)!=0)
      throw std::runtime_error("Degree found too large (>=32)");
    cache.push_back(d); // store degrees in lower $5$ bits
  }
#ifdef VERBOSE
    std::cerr << n_polynomials()-1 << '\n';
#endif
}

@ Since recorded nothing for the first two polynomials, the index into |cache|
is shifted by~$2$.
@< Methods for reading binary files @>=

size_t cached_pol_info::degree (KLIndex i) const
@/{@;
  return i<2 ? i-1 : cache[i-2]&degree_mask ;
}

@ Leading coefficients are stored into the $3$ higher order bits of the
|unsigned char| in |cache|, provided they fit. When the bits are$~0$ this either
means the coefficient has not been stored yet, or that they have but did not
fit. in either case |cahced_pol_info::leading_coeff| needs to call the base
method |polynomial_info::leading_coeff| to obtain the leading coefficient by
reading the file.

@< Methods for reading binary files @>=

size_t cached_pol_info::leading_coeff (KLIndex i) const
{
  if (i<2)
    return i;
  if ((cache[i-2]&~degree_mask)!=0) return cache[i-2]/(degree_mask+1);
  size_t lc=polynomial_info::leading_coeff(i); // look up in file
  if (lc<=255/(degree_mask+1)) cache[i-2] |= lc*(degree_mask+1);
  return lc;
}

@* The {\bf progress\_info} class.

@< Input class declarations @>=

class progress_info
{
  std::vector<KLIndex> first_pol; // count distinct polynomials in rows before
public:
  progress_info(std::ifstream& progress_file);

  BlockElt block_size() const { return first_pol.size()-1; }
  KLIndex first_new_in_row(BlockElt y) // |y==block_size()| is allowed
    const
    { return first_pol[y]; }
  BlockElt first_row_for_pol(KLIndex i) const;
};

@  Methods of the |progress_info| class
@< Methods for reading binary files @>=


progress_info::progress_info(std::ifstream& file)
: first_pol()
{ file.seekg(0,std::ios_base::end); // measure |file|
  if (file.tellg()%12!=0)
    throw std::runtime_error("Row file size not a multiple of 12");
  BlockElt size= file.tellg()/12;
  first_pol.reserve(size+1); first_pol.push_back(0);
  file.seekg(0,std::ios_base::beg); // rewind
  for (BlockElt y=0; y<size; ++y)
  { file.seekg(8,std::ios_base::cur); // skip ahead
    first_pol.push_back(first_pol.back()+read_bytes<4>(file));
  }
  file.close();
}

BlockElt progress_info::first_row_for_pol(KLIndex i) const
{ const KLIndex* p=std::upper_bound(&first_pol[1],&*first_pol.end(),i);
  return p-&first_pol[1];
}
