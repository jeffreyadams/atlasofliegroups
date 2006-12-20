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

@c
template<unsigned int n>
  class tuple_entry
  { unsigned int comp[n];
  public:
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
  };

@ For the hash function we should consider the fact that in practise most
often all components will be approximately equal, both among each other and to
the eventual sequence number this tuple will get. Since the modulus used when
the tuple is first inserted will therefore usually not be many times larger
than the values of the components, we do not attempt to apply great
multiplicative factors to the components. We shall add up shifted instances of
the components, the last one being unshifted and each previous component being
shifted to bits more than the previous one.

@c
template<unsigned int n>
  size_t tuple_entry::hashCode(size_t modulus) const
  { size_t h=0;
    for (unsigned int i=0; i<n; ++i) h=((h<<2)+comp[i])&(modulus-1);
    return h;
  }

@* Program outline.
%


@* Index.

% Local IspellDict: british