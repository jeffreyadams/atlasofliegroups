% Copyright (C) 2006 Marc van Leeuwen
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


\def\emph#1{{\it #1\/}}

@* Introduction.
%
This little program has as task to read a number of files containing the
matrices with sequence numbers of Kazhdan-Lusztig-Vogan polynomials produced
for the same block, but reducing the coefficients modulo different integers,
and to write a file containing the same matrix information that could have
been produced for that computation, but reducing the coefficients modulo the
least common multiple of those numbers. It also writes auxiliary files that
describe the correspondence between the numbering of the polynomials in the
least common multiple case and those in the base cases. Thus it deals with two
minor complications introduced by doing the modular computations: first that
modular reduction my have given certain pairs of polynomials the same number
using one of the moduli, but not using another, and second that if the
numbering produced in the modular computations should differ from the
canonical one (which could happen in multi-threaded computation) then a
correspondence has to be established between the numberings for the different
moduli.

In fact the values of the moduli play no role at all, apart from determining
the names of the output files. What we shall do is simply to traverse the
matrix files in the canonical order, and for each tuple of polynomial numbers
determine if it has occurred before, if so attribute to it the sequence number
previously associated to that tuple, or else attribute a new sequence number.

@c
@< Type definitions @>
@< Constant definitions @>
@< Function definitions @>
@< Main function @>

@* A simple template class for tuple entries.
%
We shall need to store tuples of polynomial numbers for different moduli in a
hash table; for that purpose we define a template class parametrised by the
number of components of the tuple. Since we are limiting the number of
distinct polynomials to $2^{32}$, all tuple components will be of type
|unsigned int|. Since the only thing we shall do with them is put them in a
hash table, we just define the members necessary for that.

@h <vector>
@h "../utilities/hashtable.h"

@< Type definitions @>=
template<unsigned int n>
  class tuple_entry
  { unsigned int comp[n];
  public:
    tuple_entry() // default constructor builds 0
    {@; for (unsigned int i=0; i<n; ++i) comp[i]=0; }
    tuple_entry(const std::vector<unsigned int>& comps)
    {@; for (unsigned int i=0; i<n; ++i) comp[i]=comps[i]; }
    tuple_entry(const tuple_entry& other) // copy constructor
    {@; for (unsigned int i=0; i<n; ++i) comp[i]=other.comp[i]; }
  @)
    typedef std::vector<tuple_entry> Pooltype;
    size_t hashCode(size_t modulus) const;
    bool operator!= (const tuple_entry& other) const
    { for (unsigned int i=0; i<n; ++i)
        if (comp[i]!=other.comp[i]) return true;
      return false;
    }
  @)
    unsigned int operator[] (unsigned int i) const { return comp[i]; }
  };

@ For the hash function we should consider the fact that in practise most
often all components will be approximately equal, both among each other and to
the eventual sequence number this tuple will get. Since the modulus used when
the tuple is first inserted will therefore usually not be many times larger
than the values of the components, we do not attempt to apply great
multiplicative factors to the components. We shall add up shifted instances of
the components, the last one being unshifted and each previous component being
shifted two bits more than the previous one.

@< Function definitions @>=
template<unsigned int n>
  size_t tuple_entry<n>::hashCode(size_t modulus) const
  { size_t h=0;
    for (unsigned int i=0; i<n; ++i) h=((h<<2)+comp[i])&(modulus-1);
    return h;
  }

@* Principal routines.
%
Our main task will be to read a matrix row from |n| different files, merge the
information about the location of nonzero polynomials among the primitive
pairs, and then enter tuples for all nonzero polynomials into the hash table;
the sequence numbers returned from the hash table can then be written out to
the file for the combined modulus.

@h <iostream>
@h <stdexcept>

@< Function definitions @>=
@< Auxiliary functions @>
template<unsigned int n>
  void combine_rows
    (unsigned int y,
     atlas::hashtable::HashTable<tuple_entry<n>,unsigned int>& hash,
     std::vector<std::istream*>in, std::ostream& out,
     std::vector<unsigned int>& lim)
  { for (unsigned int i=0; i<n; ++i)
      if (read_int(*in[i])!=y) @< Report alignment problem and abort @>
    unsigned int nr_prim=read_int(*in[0]);
    for (unsigned int i=1; i<n; ++i)
      if (read_int(*in[i])!=nr_prim)
        @< Report primitive count problem and abort @>
    write_int(y,out); write_int(nr_prim,out); // reproduce info in output
@)
    @< Read and merge bitmaps from the input file,
       writing the merged bitmap to |out| @>
    @< Traverse primitive elements with nonzero polynomial,
       looking up and writing the corresponding sequence number to |out| @>
  }

@ Since the bitmaps from the various files have the same capacity and
interpretation of bit positions (namely the successive primitive elements of
this row, which is independent of the modulus), merging is simply a matter of
taking the bitwise logical~``or''. Since we shall need the merged sequence and
the individual modular bitmaps in the next phase, we store them in variables.

@h "../utilities/bitmap.h"
@< Read and merge bitmaps from the input file,... @>=

typedef atlas::bitmap::BitMap bit_map;

std::vector<bit_map> in_map(n,bit_map(nr_prim));
bit_map out_map(nr_prim);

for (size_t i=0; i<nr_prim; i+=32)
{ unsigned int b=0;
  for (unsigned int j=0; j<n; ++j)
  @/{@; unsigned int bj=read_int(*in[j]); in_map[j].setRange(i,32,bj);
    b |= bj;
  }
  out_map.setRange(i,32,b); write_int(b,out);
}

@ Here we can use a bitmap-iterator over |out_map|. For each position produced
by this iterator, we read the corresponding number from the file~|in[j]| if
the corresponding bit of |in_map[j]| is set; if not we take the index~|0|
associated to the zero polynomial. Note the use of |tuple_entry<n>| below is
the sole reason that our current function is a template function.

@< Traverse primitive elements with nonzero polynomial,... @>=
for (bit_map::iterator i=out_map.begin(); i(); ++i)
{ std::vector<unsigned int> tuple(n,0);
  for (unsigned int j=0; j<n; ++j)
    if (in_map[j].isMember(*i))
    { tuple[j]=read_int(*in[j]);
      if (tuple[j]>=lim[j]) lim[j]=tuple[j]+1;
        // keep track of limit for modular numbers
    }
    else {} // |tuple[j]| stays 0; and nothing is read from |*in[j]|
  tuple_entry<n> e(tuple);
  write_int(hash.match(e),out);
}

@ When we cannot recognise the start of a row, we say which one it is and
quit.
@< Report alignment problem and abort @>=
{ std::cerr << "y=" << y << ", i=" << i << ":\n";
  throw std::runtime_error("Wrong alignment in source file");
}

@ When we find a mismatch in the number of primitive elements, we say which
one offends and quit.
@< Report primitive count problem and abort @>=
{ std::cerr << "y=" << y << ", i=" << i << ":\n";
  throw std::runtime_error("Primitive count mismatch in source files");
}

@ Here is how integers are read in and written out.

@< Auxiliary functions @>=

unsigned int read_int(std::istream& in)
{ char c; unsigned char uc; unsigned int ui,result;
@/in.get(c); result=uc=c;
@/in.get(c); ui=uc=c; result += ui<<8;
@/in.get(c); ui=uc=c; result += ui<<16;
@/in.get(c); ui=uc=c; return result + (ui<<24);
}

void write_int(unsigned int n, std::ostream& out)
{
  out.put(char(n&0xFF)); n>>=8;
  out.put(char(n&0xFF)); n>>=8;
  out.put(char(n&0xFF)); n>>=8;
  out.put(char(n));
}

@ Here is a function that will set up the hash table and the I/O streams, and
repeatedly call |combine rows|. Again it must be a template function depending
on~|n|.

@h <string>
@h <fstream>
@h <sstream>
@< Function definitions @>=
template<unsigned int n>
void do_work(std::string name_base, std::vector<unsigned int>& modulus)
{ @< Open input and output files @>
@)
  std::vector<tuple_entry<n> > pool;
  atlas::hashtable::HashTable<tuple_entry<n>,unsigned int> hash(pool);
  hash.match(tuple_entry<n>()); // insert index of Zero, it does not occur!
  std::vector<unsigned int> limits(n,1); // limit of modular sequence numbers
@)
  for (unsigned int y=0; in_stream[0]->peek()!=EOF; ++y)
     // something remains to read
  {@; std::cerr << y << '\r';
    combine_rows<n>(y,hash,in_stream,out_file,limits);
  }
  std::cerr << "\ndone!\n";
  for (unsigned int i=0; i<n; ++i) delete in_file[i]; // close files
@)
  @< Report limits of modular numbers and of generated numbers @>
  @< Write files recording the renumbering performed @>
}

@ For opening files in binary modes the following constants are useful.
@s openmode int
@< Constant definitions @>=
const std::ios_base::openmode binary_out=
			    std::ios_base::out
			  | std::ios_base::trunc
			  | std::ios_base::binary;

const std::ios_base::openmode binary_in=
			    std::ios_base::in
			  | std::ios_base::binary;

@ Opening files is easy and a bit repetitive. For input files we need pointers
in order to store them in a vector.

@< Open input and output files @>=
std::vector<std::ifstream*>in_file(n,NULL);
  std::vector<std::istream*>in_stream(n,NULL);
  for (unsigned int i=0; i<n; ++i)
  { std::ostringstream name;
    name << name_base << "-mod" << modulus[i];
    in_file[i]=new std::ifstream(name.str().c_str(),binary_in);
    if (in_file[i]->is_open())
      in_stream[i]=in_file[i]; // get stream underlying file stream
    else
    {@; std::cerr << "Could not open file '" << name.str() << "'.\n";
      exit(1);
    }
  }
@)
  unsigned long out_modulus=modulus[0];
  for (unsigned int i=1; i<n; ++i)
    @< Replace |out_modulus| by its least common multiple with |modulus[i]| @>
@)
  std::ostringstream name;
  name << name_base << "-mod" << out_modulus;
  @< Modify |name| if it coincides with that of one of the input files @>

  std::ofstream out_file(name.str().c_str(),binary_out);
  if (out_file.is_open())
    std::cout << "Output to file: " << name.str() << '\n';
  else
  @/{@; std::cerr << "Could not open output file '" << name.str() << "'.\n";
      exit(1);
    }

@ Here we use the least common multiple function |lcm| from the Atlas.

@h "../utilities/arithmetic.h"
@< Replace |out_modulus| by its least common multiple with |modulus[i]| @>=
out_modulus= atlas::arithmetic::lcm(out_modulus, modulus[i]);

@ It might happen that the output modulus coincides with one of the input
moduli, for instance if one takes twice the same input modulus. In such cases
we add a |'+'| to the file name to avoid that opening it will destroy the
input file.

@< Modify |name| if it coincides with that of one of the input files @>=
{ bool write_protect=false;
  for (ulong i=0; i<n ; ++i)
    if (out_modulus==modulus[i]) write_protect=true;
  if (write_protect) name << '+'; // avoid overwriting file for one modulus
}

@
@h <iomanip>

@< Report limits of modular numbers and of generated numbers @>=
{ std::cout << "Numbers of different polynomials found:\n";
  for (unsigned int i=0; i<n; ++i)
    std::cout << "Mod " << std::setw(10) << modulus[i]
              << ": " << limits[i] << ",\n";
  std::cout << "Mod " << std::setw(10) << out_modulus
            << ": " << hash.size() << ".\n";
}

@
@< Write files recording the renumbering performed @>=
for (unsigned int i=0; i<n; ++i)
{ std::ostringstream name;
  name << name_base << "-renumbering-mod" << modulus[i];
  std::ofstream out_file(name.str().c_str(),binary_out);
  if (out_file.is_open())
    std::cout << "Renumbering output to file: " << name.str() << '\n';
  else
    {@; std::cerr << "Could not open file '" << name.str() << "'.\n";
      exit(1);
    }
@)
  for (unsigned int k=0; k<pool.size(); ++k)
    write_int(pool[k][i],out_file);
}

@ Here is the main routine.
@< Main function @>=
int main(int argc,char** argv)
{ atlas::constants::initConstants(); // don't forget!

  --argc; ++argv; // skip program name
  std::string base;
  if (argc>0) {@; base=*argv++; --argc; }
  else
  {@; std::cout << "File name base (up to '-mod'): " ;
    std::cin >> base;
  }
@)
  std::vector<unsigned int> moduli;
@)
  if (argc==0)
  { std::cout << "Give moduli used, or 0 to terminate.\n";
    while(true)
    { std::cout << "Modulus: ";
      unsigned int m=0; std::cin >> m;
      if (m==0) break;
      moduli.push_back(m);
    }
  }
  else
    while (argc-->0)
    { std::istringstream in(*argv++);
      unsigned int m=0; in >> m;
      if (m!=0) moduli.push_back(m);
      else
      {@; std::cout << "Illegal modulus argument: " << in.str() << "\n";
	exit(1);
      }
    }


  switch (moduli.size())
  { case 1: do_work<1>(base,moduli); break;
    case 2: do_work<2>(base,moduli); break;
    case 3: do_work<3>(base,moduli); break;
    case 4: do_work<4>(base,moduli); break;
    case 5: do_work<5>(base,moduli); break;
    case 6: do_work<6>(base,moduli); break;
    default: std::cout << "I cannot handle " << moduli.size()
		       << " moduli, sorry.\n";
  }
}

@* Index.

% Local IspellDict: british