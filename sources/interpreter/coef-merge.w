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
\let\cong\equiv
\def\lcm{\mathop{\rm lcm}}

@* Introduction.
%
This little program has as task to read a number of files containing tables of
(Kazhdan-Lusztig-Vogan) polynomials with coefficients reduced modulo certain
distinct moduli, and to compute those polynomials modulo the least common
multiple of those moduli. Since the polynomials may be stored in different
orders, and moreover some might have been identified for certain moduli and
not for others, the program also inputs files for each modulus that defines a
map from a ``canonical'' numbering of the polynomials to the number they have
gotten for that modulus. Apart from I/O matters, this program is therefore
simply an implementation of the Chinese remainder theorem.


@c
@< Type definitions @>
@< Constant definitions @>
@< Function definitions @>
@< Main function @>

@ Honouring the memory of Fokko du Cloux, all calculations in this program
will be done using only |unsigned long| integer arithmetic (even when in some
cases less bits would suffice). Doing the extended Euclidean algorithm and
lifting of modular remainders without using negative numbers is an interesting
challenge that actually leads to a slick solution that is cleaner than would
be obtained using signed integer arithmetic.

@< Type definitions @>=
typedef unsigned long int ulong;

@* Extended Euclidean algorithm.
%
We lift the coefficients from remainders for given moduli to remainders for
their least common multiple, for two moduli $a,b>0$ at the time, repeating if
necessary for more moduli. In order to install this algorithm efficiently, we
first compute the greatest common divisor~$d$ of $a$ and $b$, and a
multiple~$m$ of $a$ such that $d\cong m\pmod{b}$. Then the pair of congruences
$$ \eqalign{x&\cong s\pmod{a}\cr x&\cong t\pmod{b}\cr}$$ can be solved as
follows: it is clearly necessary that $t-s$ be divisible by~$d$; if so write
$t-s=qd$, then $x=s+qm$ is easily seen to be a solution of the pair of
congruences, so that by the Chinese remainder theorem the pair is equivalent
to $$x\cong s+qm\pmod{\lcm(a,b)}.$$

To compute $d$ and $m$, and incidentally the value $\lcm(a,b)$ is done using a
simple double-exit loop. The fact that the value written to |lcm| is indeed
$\lcm(a,b)$ one uses the fact that, apart from the invariant mentioned in the
comments, the determinant $\left\vert {d_0\atop m_0}~{-d_1\atop
m_1}\right\vert=d_om_1+d_1m_0$ is also an invariant of the loop, with value |ab|.

@< Function definitions @>=

ulong extended_gcd(ulong a,ulong b, ulong&lcm, ulong& m)
{ ulong d0=a, m0=a, d1=b,m1=0;
  while(d0!=0)
  // invariant $d_0\cong+m_0\pmod{b}$, \ and $d_1\cong-m_1\pmod{b}$
  { m1+=(d1/d0)*m0; d1%=d0; // i.e., |d1-=(d1/d0)*d0|
    if (d1==0) break;
    // invariant holds, and |d0>d1>0|
    m0+=(d0/d1)*m1; d0%=d1; // i.e., |d0-=(d0/d1)*d1|
  } // invariant holds, and |0<d0<d1|
@)
  if (d1==0) {@; m=m0; lcm=m1; return d0; }
  else {@; m=m0-m1; lcm=m0; return d1; } // |d0==0|
}


@ The formula given above requires computing $q=(t-s)/d$ where $d=\gcd(a,b)$,
but we wish to avoid negative numbers in case $s>t$ (certainly just computing
$(t-s)/d$ with variables of an unsigned type would give incorrect answers).
Therefore we shall use in that case instead as example solution $x=s+q'm'$
where $q=(s-t)/d=-q$ and $m'=\lcm(a,b)-m$; it can easily be seen from the
computation above that $m'\geq0$, so that $x=s+q'm'\geq0$ in this case
(supposing $s$ is chosen non-negative, as it will be). To be able to do both
cases efficiently, we shall pre-compute both $m$ and $m'$ for every pair $a,b$
for which we shall want to lift remainders. The computation will therefore be
performed by a method of the following mini-|class|.

@< Type definitions @>=
class ChineseBox
{ ulong a,b; // the base moduli
  ulong gcd, lcm; // greatest common divisor and least common multiple
  ulong m, mm;
   // multiples of a the are congruent to |gcd| and |-gcd|, respectively
@)
  public:
  ChineseBox(ulong a, ulong b);
@)
  // accessors
  ulong get_gcd() const @+{@; return gcd; }
  ulong get_lcm() const @+{@; return lcm; }
  ulong lift_remainders(ulong s, ulong t) const;
};

@ It should now be obvious how a |ChineseBox| is constructed.
@< Function definitions @>=
ChineseBox::ChineseBox(ulong aa, ulong bb) : a(aa),b(bb)
{@; gcd=extended_gcd(a,b,lcm,m); mm=lcm-m; }

@ Lifting remainders proceeds as described above. We single out the common
case $\gcd(a,b)=1$ (which is actually a hypothesis of the Chinese Remainder
Theorem) because we can then avoid two division operations.

@h <iostream>
@< Function definitions @>=
ulong ChineseBox::lift_remainders(ulong s, ulong t) const
{ if (gcd==1)
    return (s+ (s<=t ? (t-s)*m : (s-t)*mm))%lcm;
  if ((s<=t ? t-s : s-t)%gcd==0)
     return (s+ (s<=t ? (t-s)/gcd*m : (s-t)/gcd*mm))%lcm;
  std::cerr << "Incompatible remainders " @|
              << s << " (mod " << a << ") and "
	      << t << " (mod " << b << ").\n";
  throw false;
}

@* Input and output routines.
%
Now come the more uninteresting but necessary I/O matters. Reading any number
bytes written in little-endian order requires shifting successive bytes more
and more to the left. The simplest solution is a recursive one, which is quite
acceptable given the small number of bytes to be read at once. Writing an
unsigned long value as |n| bytes in little-endian is easier, and can be done
with a simple loop. Since the function will be called quite often, we make a
small optimisation of avoiding a useless shift after writing the final byte.

@< Function definitions @>=
ulong read_bytes(ulong n, std::istream& in)
{
  if (n==0) return 0;
  char c; in.get(c); @+ unsigned char low=c; // make unsigned for arithmetic
  return low+(read_bytes(n-1,in)<<8);
}
@)
inline void write_bytes(ulong val, ulong n, std::ostream& out)
{ while (n-->1) {@; out.put(char(val&0xFF)); val>>=8; } // write |n-1| bytes
  out.put(char(val)); // and one more
}

@ Reading the renumbering table is quite trivial. Since these tables are large
(the number of polynomials times |sizeof(unsigned int)==4| bytes) we take care
to make their capacity equal to their size. This is most efficiently done by
measuring the file size before reading it and resizing the vector. Somewhat
recklessly we do not subsequently test for end-of-file (hopefully nobody is
truncating this file in the mean time).

@h <vector>
@h <fstream>

@< Function definitions @>=
void read_renumbering_table
  (std::ifstream& in, std::vector<unsigned int>& table)
{ in.seekg(0,std::ios_base::end);
  ulong file_size = in.tellg();
  in.seekg(0,std::ios_base::beg); // return to start (do not collect \$200)
@)
  table.resize(file_size/4);
  for (ulong i=0; i<table.size(); ++i) table[i]=read_bytes(4,in);
}

@ The renumbering table will be stored together with other information
pertinent to one modulus in an object of class |modulus_info|.

@s streamoff int
@< Type definitions @>=
class modulus_info
{ ulong modulus;
  ulong nr_polynomials; // number of polynomials for this modulus
  std::streamoff index_begin;
  std::streamoff coefficients_begin;
  ulong nr_coefficients;
  ulong coefficient_size; // 1 for original files, but may be more
  std::vector<unsigned int> renumber;
  std::ifstream* coefficient_file; // owned file pointer
@)
public:
  modulus_info(ulong mod, std::ifstream* ren_file, std::ifstream* coef_file);
  ~modulus_info();
@)
  ulong length(ulong i) const; // length (degree+1) of polynomial |i|
  std::vector<ulong> coefficients(ulong i) const;
    // coefficients of polynomial |i|
  const std::vector<unsigned int>& renumber_vector() const
    @+{@; return renumber; }
};

@ Constructing a |modulus_info| object requires both pertinent files to be
opened and passed to the constructor; the constructor will call
|read_renumbering_table| to set the |renumber| vector. The
|coefficient_size| will be the minimum number of bytes that can hold arbitrary
remainders for |modulus|.

@< Function definitions @>=
modulus_info::modulus_info
  (ulong mod, std::ifstream* ren_file, std::ifstream* coef_file)
  : modulus(mod), coefficient_size(1), coefficient_file(coef_file)
{ coefficient_file->seekg(0,std::ios_base::beg); // begin at the beginning
  nr_polynomials=read_bytes(4,*coefficient_file);
  index_begin=coefficient_file->tellg();
  coefficient_file->seekg(5*nr_polynomials,std::ios_base::cur);
  nr_coefficients =read_bytes(5,*coefficient_file);
  coefficients_begin=coefficient_file->tellg();

  --mod; // make largest remainder
  while ((mod>>=8) != 0) ++coefficient_size;
  read_renumbering_table(*ren_file,renumber);
  delete ren_file; // close file when table is read in
}

@ The destructor for |modulus_info| should morally do |delete
coefficient_file;|, since nobody else is holding a pointer to that
|std::ifstream| object. However, doing so would crash the program: the
|modulus_info| objects are held in a vector, and they cannot be constructed
directly into their final destination; depending on the implementation of
|std::vector| they are either copy-constructed or assigned into that vector,
so the destructor of the temporary original value would already destroy the
pointer. Therefore our destructor is a no-op. Now the file objects never get
destructed, and this could be remedied by having them owned by some global
vector of pointers instead, but since the program will terminate immediately
after our |modulus_info| objects disappear, and this closes the open output
files, we think it is not worth the hassle to do so.

@< Function definitions @>=
modulus_info::~modulus_info() @+{}

@ To get the length of a polynomial, we look up its renumbering to get the
proper index, then compare the index found with the next one.

@< Function definitions @>=
ulong modulus_info::length (ulong i) const
{ coefficient_file->seekg(index_begin+5*renumber[i],std::ios_base::beg);
    // locate index in file
  ulong index=read_bytes(5,*coefficient_file);
  ulong next_index=read_bytes(5,*coefficient_file);
  return (next_index-index)/coefficient_size;
}

@ To get actual coefficients of a polynomial, we find the index in a similar
way, then re-position the coefficients file and read the required number of
blocks of |coefficient_size| bytes.

@< Function definitions @>=
std::vector<ulong> modulus_info::coefficients (ulong i) const
{ coefficient_file->seekg(index_begin+5*renumber[i],std::ios_base::beg);
  ulong index=read_bytes(5,*coefficient_file);
  ulong next_index=read_bytes(5,*coefficient_file);
@)
  coefficient_file->seekg(coefficients_begin+index,std::ios_base::beg);
  std::vector<ulong> result ((next_index-index)/coefficient_size);
  for (ulong i=0; i<result.size(); ++i)
    result[i]=read_bytes(coefficient_size,*coefficient_file);
  return result;
}

@ Next comes the task of writing the index part of the coefficient file. All
new coefficients will have the same size, a small multiple of the size of
$1$~byte used for the original modular coefficients. We only need to inspect
the index part of the modular files to determine the polynomial degrees, and
from that and the size of the new coefficients, we can predict the new indices
to write without actually looking at the coefficients yet. However, to know
the degree of the polynomials, we have to take the maximum of the degrees of
the modular polynomials, since modular reduction may have lowered some of
those degrees.

@< Function definitions @>=
ulong write_indices
 (ulong coefficient_size,
  const std::vector<modulus_info>& mod_info,
  std::ostream& out)
@/// return value is size of (yet unwritten) coefficient part of output file
{ ulong nr_pol=mod_info[0].renumber_vector().size();
   // number of new polynomials
  write_bytes(nr_pol,4,out);
  for (ulong j=1; j<mod_info.size(); ++j)
    if (mod_info[j].renumber_vector().size()!=nr_pol)
    { std::cerr << "Conflicting numbers of polynomials in renumbering files: "
	@|      << nr_pol << "!=" << mod_info[j].renumber_vector().size()
	@|	<< " (modulus nrs O, " << j << ").\n";
      exit(1);
    }
@)
  ulong index=0; // index of current polynomial to be written
  for (ulong i=0; i<nr_pol; ++i)
  { ulong len=0; // maximum of degree+1 of polynomials selected
    for (ulong j=0; j<mod_info.size(); ++j)
    { ulong new_len = mod_info[j].length(i);
      if (new_len>len) len=new_len;
    }
    write_bytes(index,5,out); // write index for polynomial
    index+=len*coefficient_size;
    // and advance by its number of coefficient bytes
  }
  write_bytes(index,5,out); // write final index
  return index;
}

@ Finally we consider writing the coefficient part out the output file, which
is where the application of modular lifting procedure actually takes place.

@< Function definitions @>=
ulong write_coefficients
 (ulong coefficient_size,
  const std::vector<modulus_info>& mod_info,
  const std::vector<ChineseBox>& box,
  std::ostream& out)
@/// return value is maximum of lifted coefficients
{ ulong nr_pol=mod_info[0].renumber_vector().size();
  ulong max=0;
  for (ulong i=0; i<nr_pol; ++i)
  { ulong len=0; // maximum of degree+1 of polynomials selected
    std::vector<std::vector<ulong> > modular_pol;
    for (ulong j=0; j<mod_info.size(); ++j)
    { std::vector<ulong> p=mod_info[j].coefficients(i);
      modular_pol.push_back(p);
      if (modular_pol.back().size()>len) len=modular_pol.back().size();
    }
    std::vector<ulong> lifted_pol(len);
    // the polynomial modulo the lcm of the moduli
@)
    for (ulong d=0; d<len; ++d)
    { ulong c=modular_pol[0][d]; // coefficient for initial modulus
       try
       { for (ulong j=1; j<mod_info.size(); ++j)
         { ulong new_c= d>=modular_pol[j].size() ? 0 : modular_pol[j][d];
           c=box[j-1].lift_remainders(c,new_c); // it happens here!
	 }
	 if (c>max) max=c;
         write_bytes(c, coefficient_size, out);
       }
       catch (bool)
       // incompatibility found during lift; details are already printed
       { std::cerr << "In coefficient " << d
		   << " of polynomial " << i << ".\n";
         exit(1);
       }
    }
  }
  return max;
}

@* The main function.
%
Finally we must put everything together. What is left is mainly getting
arguments, and opening corresponding files, but we must not forget to create
our Chinese boxes! There will be one less of them than there are moduli.

@h <sstream>
@< Main function @>=
int main(int argc, char** argv)
{ --argc; ++argv; // skip program name
  std::string mat_base,coef_base;
  // base names for renumbering and coefficient files

  if (argc>=2) {@; mat_base=*argv++; --argc; coef_base=*argv++; --argc; }
  else test();
  // if two names are not given, go into interactive mode
@)
  if (argc<2) {@; std::cerr<< "Too few moduli"; exit(1); }
  std::vector<ulong> moduli;
  @< Get |moduli| from argument list @>
@)
  std::vector<ChineseBox> box(1,ChineseBox(moduli[0],moduli[1]));
  for (ulong i=2; i<moduli.size(); ++i)
    box.push_back(ChineseBox(box[i-2].get_lcm(),moduli[i]));
    // defines |box[i-1]|
  ulong lcm=box.back().get_lcm();
@)
  ulong coefficient_size=1, rem=lcm-1; // maximal remainder
  while ((rem>>=8)!=0) ++coefficient_size;
@)
  std::vector<modulus_info> mod_info;
   // among other things this will hold the input files
  std::ofstream coefficient_file;
@/@< Open input and output files @>

@)
  ulong nr_c=write_indices(coefficient_size,mod_info,coefficient_file);
  std::cout << "Done writing indices, will now write "
            << nr_c << " coefficient bytes.\n";
  ulong max_coef=
     write_coefficients(coefficient_size,mod_info,box,coefficient_file);
  std::cout << "Done!\nMaximal coefficient found: "
            << max_coef << ".\n";
}

@ When getting moduli, we check that they are indeed numeric and nonzero
(actually it suffices to start with a digit, trailing nonsense is ignored).

@< Get |moduli| from argument list @>=
while (argc-->0)
{ std::istringstream in(*argv++);
  unsigned int m=0; in >> m;
  if (m!=0) moduli.push_back(m);
  else
  {@; std::cerr << "Illegal modulus argument: " << in.str() << "\n";
     exit(1);
  }
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

@ Opening files is easy, but a bit long-winded anyway.

@< Open input and output files @>=
{ for (ulong i=0; i<moduli.size(); ++i)
  { std::ostringstream name0,name1;
    name0 << mat_base << "-renumbering-mod" << moduli[i];
    std::ifstream* renumber_file=
      new std::ifstream(name0.str().c_str(),binary_in);
    if (not renumber_file->is_open())
    @/{@; std::cerr << "Could not open file '" << name0.str() << "'.\n";
      exit(1);
    }
@)
    name1 << coef_base << "-mod" << moduli[i];
    std::ifstream* coef_file=new std::ifstream(name1.str().c_str(),binary_in);
    if (not coef_file->is_open())
    @/{@; std::cerr << "Could not open file '" << name1.str() << "'.\n";
      exit(1);
    }
    mod_info.push_back(modulus_info(moduli[i],renumber_file,coef_file));
  }
@)
  std::ostringstream name;
  name << coef_base << "-mod" << lcm;
  @< Modify |name| if it coincides with that of one of the input files @>
  coefficient_file.open(name.str().c_str(),binary_out);
  if (coefficient_file.is_open())
    std::cout << "Output to file: " << name.str() << '\n';
  else
  @/{@; std::cerr << "Could not open output file '" << name.str() << "'.\n";
      exit(1);
    }
}

@ It might happen that the output modulus coincides with one of the input
moduli, for instance if one takes twice the same input modulus. In Such cases
we add a |'+'| to the file name to avoid that opening it will destroy the
input file.

@< Modify |name| if it coincides with that of one of the input files @>=
{ bool write_protect=false;
  for (ulong i=0; i<moduli.size() ; ++i)
    if (lcm==moduli[i]) write_protect=true;
  if (write_protect) name << '+'; // avoid overwriting file for one modulus
}

@* The test function.
%
This chapter provides an interactive function for playing with modular lifting
interactively.

@< Function definitions @>=
void test()
{ std::vector<unsigned int> moduli;
  std::cout << "Give moduli used, or 0 to terminate.\n";
  while(true)
  { std::cout << "Modulus: ";
    ulong m=0; std::cin >> m;
    if (m==0) break;
    moduli.push_back(m);
  }
@)
  if (moduli.size()<2)
  { std::cerr << "Too few moduli.\n"; exit(1);
  }
  std::vector<ChineseBox> box(1,ChineseBox(moduli[0],moduli[1]));
  for (ulong i=2; i<moduli.size(); ++i)
    box.push_back(ChineseBox(box[i-2].get_lcm(),moduli[i]));
    // defines |box[i-1]|
  ulong lcm=box.back().get_lcm();
@)

  @< Main loop of test function @>

}

@ Here we ask space-separated sequences of remainders, and print out their
lifted values.

@< Main loop of test function @>=
std::vector<ulong> remainder(moduli.size());
while(true)
{ std::cout << "Give remainders mod " << moduli[0];
  for (ulong i=1; i<moduli.size(); ++i) std::cout << ", " << moduli[i];
  std::cout << ": ";
  for (ulong i=0; i<moduli.size(); ++i)
  { remainder[i]=~0ul; std::cin >> remainder[i];
    if (remainder[i]==~0ul)
    {@; std::cout << "Bye.\n"; exit(0); }
  }
@)
  ulong rem=remainder[0];
  try
  { for (ulong i=1; i<moduli.size(); ++i)
      rem=box[i-1].lift_remainders(rem,remainder[i]);
     std::cout << "Solution: " << rem << " (mod " << lcm << ").\n";
  }
  catch (bool)
  { std::cout << "No solution.\n"; }
}

@* Index.

% Local IspellDict: british