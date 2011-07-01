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

The prerequisites for running this program will be to have, for some
collection of moduli~$k$ and fixed prefixes $P$~and~$M$, files $P$\.{-mod}$k$
of polynomials, as written by the |klwrite| command of the modular version of
the Atlas software, and files $M$\.{-renumbering-mod}$k$, as produced by the
\.{matrix-merge} program. We shall also provide an interactive mode in which
the modular remainders can be entered manually by the user.


@c
@< Type definitions @>
@< Constant definitions @>
@< Function definitions @>
@< Main function @>

@ Honouring the memory of Fokko du Cloux, all calculations in this program
will be done using only |unsigned long| integer arithmetic; no negative
numbers will occur anywhere at any time. We shall use the |ulong| type
even when in some cases less bits would probably suffice. Doing the extended
Euclidean algorithm and lifting of modular remainders without using negative
numbers is an interesting challenge that actually leads to a slick solution
that is cleaner than would be obtained using signed integer arithmetic.

However, we want this program to be able to run on $32$-bit architectures even
for processing the files of polynomials for split~$E_8$, whose size is larger
than $2^{32}$ bytes. Since on those architectures the type |unsigned long int|
is usually only $32$~bits wide, we need to use a different type when computing
positions in a file. The appropriate type to handle such computations is
|std::streamoff|, which should be large enough if the underlying operating
system is capable of handling files of more than $2^{32}$ bytes at all; in
practice this means it is almost certainly a $64$-bit unsigned integral type.

As a consequence of these decisions, this program should be capable of merging
arbitrary files of polynomials for up to 4 byte-size moduli on a $32$-bit
architecture, and for up to 8 byte-size moduli on a $64$-bit architecture.


@h <iostream>
@< Type definitions @>=
typedef unsigned long int ulong;
typedef std::streamoff file_pos;

@* Extended Euclidean algorithm.
%
We lift the coefficients from remainders for given moduli to remainders for
their least common multiple, for two moduli $a,b>0$ at the time, repeating if
necessary for more moduli. Our code will work for arbitrary pairs $a,b$,
subject only in some cases to conditions imposed by the limitations of
unsigned long integral arithmetic, but our efficiency considerations will be
based on the assumption that |a<=b|. We shall later try to arrange things so
that this condition is likely to hold, but we do not reorder the moduli
provided by the user, so that one can benefit from the asymmetry in $a$
and~$b$ that will be present in our algorithms to make consistency checks by
permuting the moduli.

In order to install this algorithm efficiently, we
first compute the greatest common divisor~$d$ of $a$ and $b$, and a
multiple~$m$ of $b$ such that $d\cong m\pmod{a}$. Then the pair of congruences
$$ \eqalign{x&\cong s\pmod{a}\cr x&\cong t\pmod{b}\cr}$$ can be solved as
follows: it is clearly necessary that $s-t$ be divisible by~$d$; if so write
$s-t=qd$, then $x=t+qm$ is easily seen to be a solution of the pair of
congruences, so that by the Chinese remainder theorem the pair is equivalent
to $$x\cong t+qm\pmod{\lcm(a,b)}.$$

To compute $d$ and $m$, and incidentally the value $\lcm(a,b)$ is done using a
simple double-exit loop. The fact that the value written to |lcm| is indeed
$\lcm(a,b)$ one uses the fact that, apart from the invariants mentioned in the
comments, the determinant $\left\vert {d_0\atop-m_0}~{d_1\atop m_1}\right\vert
=d_om_1+d_1m_0$ is also an invariant of the loop, with value~|ab|.

@< Function definitions @>=

ulong extended_gcd(ulong a,ulong b, ulong&lcm, ulong& m)
{ ulong d0=a, m0=0, d1=b,m1=b;
  while(d0!=0)
@/ // invariants: $d_0\cong-m_0\pmod{a}$, \ and $d_1\cong+m_1\pmod{a}$
  { m1+=(d1/d0)*m0; d1%=d0; // i.e., |d1-=(d1/d0)*d0|
    if (d1==0) break;
@/   // invariants hold, and |d0>d1>0|
    m0+=(d0/d1)*m1; d0%=d1; // i.e., |d0-=(d0/d1)*d1|
  } // invariants hold, and |0<d0<d1|
@)
  if (d1==0) {@; lcm=m1; m=m1-m0; return d0; }
  else {@; lcm=m0; m=m1; return d1; } // |d0==0|
}


@* Chinese boxes.
%
The formula given above requires computing $q=(s-t)/d$ where $d=\gcd(a,b)$,
but we wish to avoid negative numbers in case $s<t$ (certainly just computing
$(s-t)/d$ with variables of an unsigned type would give incorrect answers).
Therefore we shall use in that case instead as example solution $x=t+q'm'$
where $q'=(t-s)/d=-q$ and $m'=\lcm(a,b)-m$; it can easily be seen from the
computation above that $m'\geq0$, so that $x=t+q'm'\geq0$ in this case
(supposing $t$ is chosen non-negative, as it will be). To be able to do both
cases efficiently, we shall pre-compute both $m$ and $m'$ for every pair $a,b$
for which we shall want to lift remainders. The computation will therefore be
performed by the method~|lift_remainders| of the following mini-|class|, which
performs lifting from modulo~|a| and modulo~|b| to modulo~|lcm|. We make this
method virtual, so that we can easily provide alternative implementations of
the lifting by defining classes derived from |ChineseBox| that redefine this
method. We also make the data members protected rather than private, for the
convenience of derived classes.

@< Type definitions @>=
class ChineseBox
{ protected:
  ulong a,b; // the base moduli
  ulong gcd, lcm; // greatest common divisor and least common multiple
  ulong m, mm;
   // multiples of |b| congruent to |gcd| and |-gcd|, respectively, modulo~|a|
@)
  public:
  ChineseBox(ulong a, ulong b);
  virtual ~ChineseBox() @+ {}
@)
  // accessors
  ulong get_gcd() const @+{@; return gcd; }
  ulong get_lcm() const @+{@; return lcm; }
  virtual ulong lift_remainders(ulong s, ulong t) const;
};

@ It should now be obvious how a |ChineseBox| is constructed.
@< Function definitions @>=
ChineseBox::ChineseBox(ulong aa, ulong bb) : a(aa),b(bb)
{@; gcd=extended_gcd(a,b,lcm,m); mm=lcm-m; }

@ Lifting remainders proceeds as described above. For the moment we only
define a basic case, that works for any pair $a,b$ so that the numbers
|(a/gcd-1)*m| and |b-1+(b/gcd-1)*mm| can be represented by an |ulong| value.

@< Function definitions @>=
ulong ChineseBox::lift_remainders(ulong s, ulong t) const
{ if ((s<=t ? t-s : s-t)%gcd==0)
     return (t+ (s>=t ? (s-t)/gcd*m : (t-s)/gcd*mm))%lcm;
  std::cerr << "Incompatible remainders " @|
              << s << " (mod " << a << ") and "
	      << t << " (mod " << b << ").\n";
  throw false;
}

@ In the common case that $\gcd(a,b)=1$ (which is actually a hypothesis of the
Chinese Remainder Theorem) we can then avoid two division operations in the
above code, so we shall define a derived class |PrimeChineseBox| of
|ChineseBox| that handles this special case more efficiently.

@h <cstdlib>
@< Type definitions @>=
class PrimeChineseBox : public ChineseBox
{ // no additional data
public: // constructor from basic |ChineseBox|
  PrimeChineseBox(const ChineseBox& cb) : ChineseBox(cb)
  { if (gcd!=1)
    { std::cerr << "Non relatively prime numbers, gcd=" << gcd << ".\n";
      std::exit(1); // we should never get here
    }
  }
  virtual ~PrimeChineseBox() @+ {}
@)
  // redefined accessor
  virtual ulong lift_remainders(ulong s, ulong t) const;
};

@ As we said, there is no reason to divide by the |gcd|, nor can there be
incompatible remainders.

@< Function definitions @>=
ulong PrimeChineseBox::lift_remainders(ulong s, ulong t) const
@/{@; return (t+ (s>=t ? (s-t)*m : (t-s)*mm))%lcm; }

@*1 Tabled Chinese boxes.
%
The |lift_remainder| classes of the previous classes depend on computing the
product |(s-t)/gcd*m| or |(t-s)/gcd*mm| as an unsigned long integer, before
reduction modulo the~|lcm|; this will have a serious risk over integral
overflow if |a| and~|b| are even moderately large with respect to the capacity
of such integers: even if both |a| and~|b|, and therefore |t-s| or |s-t|, can
be represented with less than half the number of available bits, the product
might overflow, since |m| and |mm| represent fixed classes modulo~$\lcm(a,b)$
and might require more than half the number of available bits. For this reason
we shall now present some classes in which the multiplications by |m| and~|mm|
modulo the |lcm| are performed via precomputed tables, which allows us to
prevent integral overflow whenever the |lcm| itself can be represented by an
unsigned long integer.

The basic case is handled by the following class. There is a table |m_table|,
that serves to store the products modulo~|lcm| of |m| by numbers up to
|a/gcd|, which is a complete table for this number since |m| is a multiple
of~|b| so that |(a/gcd)*m| give |0| modulo~|lcm|. This table can be used
directly to evaluate |(s-t)/gcd*m| since |s/gcd<a/gcd|; for |(t-s)/gcd*mm| we
can use the same table indirectly via the reverse iterator |mm_table|, since
|mm| represents |-m| modulo~|lcm|.

@< Type definitions @>=
class TabledChineseBox : public ChineseBox
{
protected:
  std::vector<ulong> m_table;
  std::vector<ulong>::const_reverse_iterator mm_table;
public: // constructor from basic |ChineseBox|
  TabledChineseBox(const ChineseBox& cb);
  virtual ~TabledChineseBox() @+ {}  // destroys table
@)
  // redefined accessor
  virtual ulong lift_remainders(ulong s, ulong t) const;
};

@ The multiplication table for~|m| is installed at construction, and for
multiplication by~|mm| we use the same table in reverse. Although not strictly
necessary, we include a final entry |m_table[a/gcd]=0| which will also be
accessible as |mm_table[0]|, so that |mm_table[i]| gives the class
modulo~|lcm| of |i*mm| for all |0<=i<a/gcd|. The construction of the table
must of course be done without actually computing products, since that would
already corrupt the tables by the integral overflow we are trying to avoid. In
the code below we use the fact that |m+mm==lcm|, and arrange to avoid even
negative intermediate results (although that would not have caused any real
problem, as long as only additive operations are used).

@h <cassert>
@< Function definitions @>=
TabledChineseBox::TabledChineseBox(const ChineseBox& cb)
: ChineseBox(cb), m_table(a/gcd+1), mm_table(m_table.rbegin())
{ ulong last = m_table[0]=0;
  for (ulong i=1; i<m_table.size(); ++i)
    m_table[i]= last<mm ? last+=m : last-=mm;
  assert(last==0);
}

@ The lifting routine is now quite simple. There are no longer any modular
reductions by |lcm|. The expression |(s-t)/gcd*m| can be replaced by a direct
look-up in |m_table| since |(s-t)/gcd<a/gcd| always holds. However to replace
the expression |(t-s)/gcd*mm| by a look-up via |mm_table| we need to reduce
the left factor modulo~|a/gcd| to stay inside the table, since |t-s| can be as
large as~|b-1| which could greatly exceed~|a|; this reduction is most easily
obtained by reduction of |t-s| modulo~|a| before division by~|gcd| (the
validity of this reduction can also be seen from the fact that this is the
only place where |s| enters in the formula, whose value is only meaningful
modulo~|a| anyway). We have installed entries~|0| at both ends of |m_table|,
the first of which is only accessed by subscripting |m_table| while the last
is only accessed by subscripting |m_table|. We could have done without the
final~|0|, with |mm_table| pointing one place earlier, if the second variant
defining~|d| below were written |mm_table[(~s+t)%a/gcd]| (where |~s+t| is a
slick way to compute |t-1-s|), but in this case we have preferred to spend one
entry extra for the sake of readability.

@< Function definitions @>=
ulong TabledChineseBox::lift_remainders(ulong s, ulong t) const
{ if ((s>=t ? s-t : t-s)%gcd==0)
  { ulong d=s>=t ? m_table[(s-t)/gcd] : mm_table[(t-s)%a/gcd];
    return t<lcm-d ? t+d : t-(lcm-d);
  }
  std::cerr << "Incompatible remainders " @|
              << s << " (mod " << a << ") and "
	      << t << " (mod " << b << ").\n";
  throw false;
}

@ The case $\gcd(a,b)=1$ can be optimised for |TabledChineseBox| as it could
for |ChineseBox|. The table is constructed by the parent class, so our
constructor only has to add the |gcd| test.

@< Type definitions @>=
class PrimeTabledChineseBox : public TabledChineseBox
{ // no additional data
public: // constructor from basic |ChineseBox|
  PrimeTabledChineseBox(const ChineseBox& cb) : TabledChineseBox(cb)
  { if (gcd!=1)
    { std::cerr << "Non relatively prime numbers, gcd=" << gcd << ".\n";
      exit(1); // we should never get here
    }
  }
  virtual ~PrimeTabledChineseBox() @+ {}
@)
  // redefined accessor
  virtual ulong lift_remainders(ulong s, ulong t) const;
};

@ Again the lifting method is a simplification of the basic case.

@< Function definitions @>=
ulong PrimeTabledChineseBox::lift_remainders(ulong s, ulong t) const
{ ulong d=s>=t ? m_table[s-t] : mm_table[(t-s)%a];
@/return t<lcm-d ? t+d : t-(lcm-d);
}

@*1 Double tabled boxes.
%
As an option of maximal speed at the price of more memory, we define a class
with separate tables for multiplication by |m| and by~|mm|. This avoids all
divisions other than by the |gcd|, essentially by repeating entries in the
larger of the two tables until the largest possible index can be accommodated
without wrapping around. This option is very expensive in memory when the two
moduli have vastly different sizes (supposing that in the case of a single
table the smallest modulus is chosen to be the first one; if not then that
case also requires a lot of memory).

@< Type definitions @>=
class DoubleTabledChineseBox : public ChineseBox
{
protected:
  std::vector<ulong> m_table,mm_table;
public: // constructor from basic |ChineseBox|
  DoubleTabledChineseBox(const ChineseBox& cb);
  virtual ~DoubleTabledChineseBox() @+ {}  // destroys tables
@)
  // redefined accessor
  virtual ulong lift_remainders(ulong s, ulong t) const;
};

@ The multiplication tables are installed at construction. The two cases are
now quite similar, and no extra entry is needed at the far end of either
table, so the sizes are determined by the maximal possible values of
|(s-t)/gcd| and |(t-s)/gcd|, respectively. In the code below we use the fact
that |m+mm=lcm|, and arrange to avoid even negative intermediate results.

@< Function definitions @>=
DoubleTabledChineseBox::DoubleTabledChineseBox(const ChineseBox& cb)
: ChineseBox(cb), m_table(a/gcd), mm_table(b/gcd)
{ ulong last = m_table[0]=0;
  for (ulong i=1; i<m_table.size(); ++i)
    m_table[i]= last<mm ? last+=m : last-=mm;
  assert(last==mm);
  last = mm_table[0]=0;
  for (ulong i=1; i<mm_table.size(); ++i)
    mm_table[i]= last<m ? last+=mm : last-=m;
}

@ The lifting routine is now as simple as it gets for general~$a,b$, without
any modular reductions.

@< Function definitions @>=
ulong DoubleTabledChineseBox::lift_remainders(ulong s, ulong t) const
{ if ((s>=t ? s-t : t-s)%gcd==0)
  { ulong d=s>=t ? m_table[(s-t)/gcd] : mm_table[(t-s)/gcd];
    return t<lcm-d ? t+d : t-(lcm-d);
  }
  std::cerr << "Incompatible remainders " @|
              << s << " (mod " << a << ") and "
	      << t << " (mod " << b << ").\n";
  throw false;
}

@ The case $\gcd(a,b)=1$ can be optimised as before. The table is constructed
by the parent class, so our constructor only has to add the |gcd| test.

@< Type definitions @>=
class PrimeDoubleTabledChineseBox : public DoubleTabledChineseBox
{ // no additional data
public: // constructor from basic |ChineseBox|
  PrimeDoubleTabledChineseBox(const ChineseBox& cb)
  : DoubleTabledChineseBox(cb)
  { if (gcd!=1)
    { std::cerr << "Non relatively prime numbers, gcd=" << gcd << ".\n";
      exit(1); // we should never get here
    }
  }
  virtual ~PrimeDoubleTabledChineseBox() @+ {}
@)
  // redefined accessor
  virtual ulong lift_remainders(ulong s, ulong t) const;
};

@ Again the lifting method is a simplification of the basic case.

@< Function definitions @>=
ulong PrimeDoubleTabledChineseBox::lift_remainders(ulong s, ulong t) const
{ ulong d=s>=t ? m_table[s-t] : mm_table[t-s];
@/return t<lcm-d ? t+d : t-(lcm-d);
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
file_pos read_bytes(ulong n, std::istream& in)
{
  if (n==0) return 0;
  char c; in.get(c); @+ unsigned char low=c; // make unsigned for arithmetic
  return low+(read_bytes(n-1,in)<<8);
}
@)
inline void write_bytes(file_pos val, ulong n, std::ostream& out)
{ while (n-->1) {@; out.put(char(val&0xFF)); val>>=8; } // write |n-1| bytes
  out.put(char(val)); // and one more
}

@ Reading the renumbering table is quite trivial; the function
|read_renumbering_table| will read |nr_pols| $4$-byte numbers from the file
|in| into the vector~|table|. The table will become quite large, namely
|4*nr_pols| bytes. For this reason this function cannot succeed on a $32$-bit
system when processing the more than $1.1*10^9$ polynomials for split~$E_8$;
for that case (and others where there is too little memory) we shall provide a
safer (but slower) alternative solution. The simplest way to find out of there
is sufficient memory to store the table is to try; if there isn't then
|resize| will either throw a |std::length_error| or a |std::bad_alloc|
exception, which may be caught by one of our callers. (More precisely, it will
throw |std::length_error| if the requested size exceeds the maximal possible
size of a vector of this type, as would be the case for split~$E_8$ on a
$32$-bit system, while it will throw |std::bad_alloc| if a vector of the
requested size could exist on the current architecture, but the current amount
of (available) virtual memory does not allow it to be allocated.)

Because of the great size of the tables, we take care to make their
|capacity()| equal to their final |size()| (so no space is uselessly
reserved). The number |nr_pols| used here will in fact have been obtained by
measuring the file size and dividing by~$4$, as we shall see below. Somewhat
recklessly we do not subsequently test for end-of-file (hopefully nobody is
truncating this file in the mean time).

@h <vector>
@h <fstream>

@< Function definitions @>=
void read_renumbering_table
  (ulong nr_pols, std::ifstream& in, std::vector<unsigned int>& table)
  throw (std::length_error,std::bad_alloc)
{ table.resize(nr_pols);
  for (ulong i=0; i<nr_pols; ++i) table[i]=read_bytes(4,in);
}

@ The renumbering table will be stored together with other information
pertinent to one modulus in an object of class |modulus_info_with_table|.
Because of the impossibility to store renumbering tables on certain systems,
we shall however first define a base type |modulus_info| that does without
such tables, and then derive from it a type that does contain and use them,
for use on those systems that do have sufficient virtual memory.

@< Type definitions @>=
class modulus_info
{ protected:
  ulong modulus;
  ulong nr_polynomials;
    // number of polynomials for this modulus, in fact $<2^{32}$
  file_pos index_begin;
  file_pos coefficients_begin;
  file_pos nr_coefficients;
  ulong coefficient_size; // 1 for original files, but may be more
  bool using_renumber,owning_files;
  std::ifstream& renumbering_file; // usually owned file reference
  std::ifstream& coefficient_file; // usually owned file reference
@)
public:
  modulus_info(ulong mod, std::ifstream* ren_file, std::ifstream* coef_file);
  virtual ~modulus_info();
private:
  modulus_info(const modulus_info&); // forbid copy constructor
public:
@) // accessors
  virtual ulong length(file_pos i) const;
  // length (degree+1) of polynomial |i|
  virtual std::vector<ulong> coefficients(file_pos i) const;
    // coefficients of polynomial |i|
  ulong nr_pol() const @+{@; return nr_polynomials; }
@)  // manipulator
  bool set_owning_files(bool b)
  @+{@; bool old=owning_files; owning_files=b; return old; }
};

@ Constructing a |modulus_info| object requires both pertinent files to be
opened and passed to the constructor. The constructor will store references to
both files without for the moment reading from the~|renumbering_file| (in
fact, it might be a null reference). The flag |owning_files| is set to |true|
signalling that normally it is our responsibility to close these files upon
destruction. From the coefficient file the constructor reads some bytes at the
beginning and at the end of the index part, in order to determine the
necessary offsets into the file and the number of coefficients that it
contains. The |coefficient_size| is computed as the minimum number of bytes
that can hold arbitrary remainders for |modulus| (we might check here as well
that this is the number of bytes used to store the unique coefficient of
polynomial number~|1|, which should represent the constant polynomial~|1|).
The number of polynomials for this modulus, which is read from the coefficient
file, is used to determine the end of the index part, but it is not
necessarily the correct number |nr_polynomials| for the final modulus, unless
the coefficient file already uses the canonical numbering. If it is not, it is
the (length of) the renumbering file that determines |nr_polynomials|.

@< Function definitions @>=
modulus_info::modulus_info
  (ulong mod, std::ifstream* ren_file, std::ifstream* coef_file)
  : modulus(mod), nr_polynomials(), index_begin(), coefficients_begin()
  , coefficient_size(1)
@/, using_renumber(ren_file!=NULL), owning_files(true)
@/, renumbering_file(*ren_file)
  , coefficient_file(*coef_file)
{ coefficient_file.seekg(0,std::ios_base::beg); // begin at the beginning
  file_pos nr_mod_pols=read_bytes(4,coefficient_file);
    // number of polynomials for |mod|
  index_begin=coefficient_file.tellg();
  coefficient_file.seekg(5*file_pos(nr_mod_pols),std::ios_base::cur);
  nr_coefficients =read_bytes(5,coefficient_file);
  coefficients_begin=coefficient_file.tellg();
@)
  if (using_renumber)
    @< Determine |nr_polynomials| by measuring the |renumber_file| @>
  else nr_polynomials=nr_mod_pols;
    // in canonical order, assume there were no collisions
@)
  --mod; // make largest remainder
  while ((mod>>=8) != 0)
    ++coefficient_size; // compute number of bytes required
}

@ When we are constructing a |modulus_info| for a modulus for which a
|renumber_file| is present, it is the size of that file rather than the number
of polynomials stored in the |coefficient_file| that determines
|nr_polynomials|; the latter number could be smaller due to collisions caused
by modular reduction, while |nr_polynomials| should match the number of
distinct polynomials \emph{after} lifting of the coefficients. After measuring
the size of the file we rewind it, so that |read_renumbering_table| will read
the file contents from the beginning.

@< Determine |nr_polynomials| by measuring the |renumber_file| @>=
{ renumbering_file.seekg(0,std::ios_base::end);
  nr_polynomials= renumbering_file.tellg()/4;
  renumbering_file.seekg(0,std::ios_base::beg);
    // return to start (do not collect \$200)
}

@ The destructor for |modulus_info| normally deletes the file pointers for
which it holds references, thus closing those files (if it has not already
been done explicitly before). This operation would cause a problem if we were
to use an object of type |std::vector<modulus_info>|, since in that case the
|modulus_info| objects, even if written as the argument to |push_back|, are
not constructed directly into their final destination, but rather
copy-constructed or assigned there (we did not investigate exactly how
|push_back| is implemented). The consequence would be that the the destructor
of the temporary original value is called, and would already destroy the
pointer (and close the associated file) before the copy even gets to use it,
and the destructor of the copy would cause the program to crash.

However, since we have stored references, not pointers, it is impossible to
create an object of type |std::vector<modulus_info>|, since this requires that
an assignment operation be defined, and this cannot be done (at least not in
the standard way) for an object containing a reference. Therefore this problem
cannot arise, and for the same reason we have forbidden the copy constructor,
so that it will not crop up in another setting either. There is however one
more scenario that would cause problems: when we try to build a
|modulus_info_with_table| object (a class derived from |modulus_info|) and the
construction of the base object succeeds, but the allocation of the table
fails for shortage of memory, then the destruction of the already completed
base object would close the files, thwarting an attempt to salvage the
situation by creating and working with just a |modulus_info| object instead.
Therefore we have introduced the flag |owning_flags| that can be made
temporarily false at the point where we may be thrown out of the constructor
for the derived class.

Of course we do need another solution to handle a sequence of |modulus_info|
objects; the solution is to use an object of type |std::vector<modulus_info*>|
instead, and to arrange that at its destruction |delete| is called for the
contained pointers. Then the pointed to |modulus_info| objects are never
copied or assigned, and their destructor can safely call |delete|. An
additional advantage of using pointers is that we will be able to use virtual
methods appropriately.

@< Function definitions @>=
modulus_info::~modulus_info()
{@; if (owning_files)
   {@; delete &renumbering_file; delete &coefficient_file;}}

@ To get the length of a polynomial, we look up its renumbering if appropriate
to get the proper index, then compare the index found with the next one. It is
important that the parameter of this function be declared a |file_pos| rather
than an |ulong|, although the argument value provided will always fit in
$32$-bits; the reason is that otherwise the multiplications below by
$4$~and~$5$ could overflow if |ulong| has only $32$-bits (which is true on
most $32$-bit architectures).

@< Function definitions @>=
ulong modulus_info::length (file_pos i) const
{ if (using_renumber)
  @/{@; renumbering_file.seekg(4*i,std::ios_base::beg);
    i=read_bytes(4,renumbering_file);
  }
  coefficient_file.seekg(index_begin+5*i,std::ios_base::beg);
    // locate index in file
  file_pos index=read_bytes(5,coefficient_file);
  file_pos next_index=read_bytes(5,coefficient_file);
  return (next_index-index)/coefficient_size;
}

@ To get actual coefficients of a polynomial, we find the index in a similar
way, then re-position the coefficients file and read the required number of
blocks of |coefficient_size| bytes.

@< Function definitions @>=
std::vector<ulong> modulus_info::coefficients (file_pos i) const
{ if (using_renumber)
  @/{@; renumbering_file.seekg(4*i,std::ios_base::beg);
    i=read_bytes(4,renumbering_file);
  }
  coefficient_file.seekg(index_begin+5*i,std::ios_base::beg);
  file_pos index=read_bytes(5,coefficient_file);
  file_pos next_index=read_bytes(5,coefficient_file);
@)
  coefficient_file.seekg(coefficients_begin+index,std::ios_base::beg);
  std::vector<ulong> result ((next_index-index)/coefficient_size);
  for (ulong i=0; i<result.size(); ++i)
    result[i]=read_bytes(coefficient_size,coefficient_file);
  return result;
}

@ Next we derive the class |modulus_info_with_table| from |modulus_info|.
Apart from the renumbering table and new versions of the virtual functions,
there is not much to this derivation. Note that the constructor takes the same
arguments as that of |modulus_info|, so that the base object can be
constructed in place; there is no choice since we have forbidden the copy
constructor for~|modulus_info|.

@< Type definitions @>=

class modulus_info_with_table : public modulus_info
{
  std::vector<unsigned int> renumber;
@)
public:
  modulus_info_with_table
   (ulong mod, std::ifstream* ren_file, std::ifstream* coef_file);
  virtual ~modulus_info_with_table() @+{} // nothing extra to destruct
@)
  virtual ulong length(file_pos i) const; // length of polynomial |i|
  virtual std::vector<ulong> coefficients(file_pos i) const;
    // coefficients of polynomial |i|
};


@ Constructing a |modulus_info_with_table| constructs a base |modulus_info|
with the arguments provided. If required, the constructor will then call
|read_renumbering_table| to set the |renumber| vector. It may happen that
|read_renumbering_table| throws a |std::bad_alloc| exception, which we let
escape to be caught higher up; note that |renumbering_file| will not be closed
in that case, so that we may retry later to create just a |modulus_info|
object with the same arguments, which requires considerably less memory.

@< Function definitions @>=
modulus_info_with_table::modulus_info_with_table
  (ulong mod, std::ifstream* ren_file, std::ifstream* coef_file)
  : modulus_info(mod,ren_file,coef_file)
  , renumber()
{ if (using_renumber)
  { set_owning_files(false); // don't close files if next call fails
    read_renumbering_table(nr_polynomials,renumbering_file,renumber);
    renumbering_file.close(); // close file when table is read in
    set_owning_files(true); // upon success, base object takes charge again
  }
}

@ The virtual methods |length| and |coefficients| differ from those in the
base class by replacing a look-up into a file by a subscription of the vector
|renumber|. We believe that this significantly speeds up processing, since it
reduces the number of different places where files are being read from 2~per
modulus to 1~per modulus for |length|, and from 3~per modulus to 2~per modulus
for |coefficients|; one should keep in mind that for reading, unlike for
writing, we have to wait for the disk access to complete, and that reading at
multiple places involves important disk-head motion.

@< Function definitions @>=
ulong modulus_info_with_table::length (file_pos i) const
{ if (using_renumber) i=renumber[i];
  coefficient_file.seekg(index_begin+5*i,std::ios_base::beg);
    // locate index in file
  file_pos index=read_bytes(5,coefficient_file);
  file_pos next_index=read_bytes(5,coefficient_file);
  return (next_index-index)/coefficient_size;
}

@ The |coefficients| method relates to |length| in exactly the same way as it
did for the base class.

@< Function definitions @>=
std::vector<ulong> modulus_info_with_table::coefficients (file_pos i) const
{ if (using_renumber) i=renumber[i];
  coefficient_file.seekg(index_begin+5*i,std::ios_base::beg);
  file_pos index=read_bytes(5,coefficient_file);
  file_pos next_index=read_bytes(5,coefficient_file);
@)
  coefficient_file.seekg(coefficients_begin+index,std::ios_base::beg);
  std::vector<ulong> result ((next_index-index)/coefficient_size);
  for (ulong i=0; i<result.size(); ++i)
    result[i]=read_bytes(coefficient_size,coefficient_file);
  return result;
}


@* Writing the coefficient file.
%
Next comes the task of writing the index part of the coefficient file. All
new coefficients will have the same size, a small multiple of the size of
$1$~byte used for the original modular coefficients. We only need to inspect
the index part of the modular files to determine the polynomial degrees, and
from that and the size of the new coefficients, we can predict the new indices
to write without actually looking at the coefficients yet. However, to know
the degree of the polynomials, we have to take the maximum of the degrees of
the modular polynomials, since modular reduction may have lowered some of
those degrees.

@< Function definitions @>=
file_pos write_indices
 (ulong coefficient_size,
  const std::vector<modulus_info*>& mod_info,
  std::ostream& out,
  bool verbose, bool output)
@/// return value is size of (yet unwritten) coefficient part of output file
{ ulong nr_pols=mod_info[0]->nr_pol(); // number of new polynomials
  write_bytes(nr_pols,4,out);
  for (ulong j=1; j<mod_info.size(); ++j)
    if (mod_info[j]->nr_pol()!=nr_pols)
    @< Complain about wrong number of polynomials, and quit @>
@)
  file_pos index=0; // index of current polynomial to be written
  for (ulong i=0; i<nr_pols; ++i)
  { if (verbose and (i&0xFFF)==0)
      std::cerr << "Polynomial: " << std::setw(10) << i << '\r';
    ulong len=0; // maximum of degree+1 of polynomials selected
    for (ulong j=0; j<mod_info.size(); ++j)
    { ulong new_len = mod_info[j]->length(i);
      @< Check for too large a polynomial; if one is found, quit @>
      if (new_len>len) len=new_len;
    }
    if (output) write_bytes(index,5,out); // write index for polynomial
    index+=len*coefficient_size;
    // and advance by its number of coefficient bytes
  }
  if (verbose) // make final polynomial displayed look right, and stay visible
    std::cerr << "Polynomial: " << std::setw(10) << nr_pols-1 << '\n';
  if (output) write_bytes(index,5,out); // write final index
  return index;
}

@ We must have equal size renumbering files to be confident that we can match
up corresponding polynomials by indexing the renumbering files at the same
index. The size of a renumbering file could be too small if there were
collisions of polynomials even at the combined modulus when one of the
renumbering files was produced by \.{matrix-merge}; if however that program
was run for the same moduli as are used here, then this should not happen.

@< Complain about wrong number of polynomials, and quit @>=
{ std::cout << "Conflicting numbers of polynomials in renumbering files: "
    @|      << nr_pols << "!=" << mod_info[j]->nr_pol()
    @|	<< " (modulus nrs O, " << j << ").\n";
  exit(1);
}

@ The code below is specific for split~$E_8$, and tries to detect anomalies
from misaligned reads. It is curious to note that such problems did occur at
some time, but by the time this code was added to diagnose the situation, the
problem had already been accidentally corrected; therefore this code was never
triggered.

@h <stdexcept>
@< Check for too large a polynomial... @>=
if (new_len>32)
{ std::cout << "Too large length " << new_len << " in polynomial "
	    << i @| << " for modulus nr " << j << ".\n";
@/throw std::logic_error("Oversize polynomial");
}

@ Finally we consider writing the coefficient part out the output file, which
is where the application of modular lifting procedure actually takes place.
The |ChineseBox|es are passed as a vector of pointers, so that the virtual
method |lift_remainders| will be properly called. Lifting is performed in a
bottom-up tree-like pattern, which tries to combine moduli that have the same
distance from the original moduli, and therefore approximately the same size.

@h <iomanip>

@< Function definitions @>=
ulong write_coefficients
 (ulong coefficient_size,
  const std::vector<modulus_info*>& mod_info,
  const std::vector<ChineseBox*>& box,
  std::ostream& out,
  bool verbose, bool output)
@/// return value is maximum of lifted coefficients
{ ulong nr_pols=mod_info[0]->nr_pol();
  ulong n=mod_info.size(); // number of moduli
  ulong max=0;
  std::vector<ulong> rem(2*n-1);
      // remainders for |n| original and |n-1| derived moduli
  std::vector<std::vector<ulong> > modular_pol(n);
   // polynomials from the base moduli
@)
  for (ulong i=0; i<nr_pols; ++i)
  { if (verbose and (i&0xFFF)==0)
      std::cerr << "Polynomial: " << std::setw(10) << i << '\r';
    ulong len=0; // maximum of degree+1 of polynomials selected

@)
    @< Read polynomial coefficients into |modular_pol|,
       and compute their maximal length in |len| @>
    std::vector<ulong> lifted_pol(len);
    // the polynomial modulo the lcm of the moduli
@)
    for (ulong d=0; d<len; ++d)
    { for (ulong j=0; j<n; ++j) // install original remainders
        rem[j]= d>=modular_pol[j].size() ? 0 : modular_pol[j][d];
       try
       { for (ulong j=0; j<n-1; ++j)
           rem[n+j]=box[j]->lift_remainders(rem[2*j],rem[2*j+1]);
           // it happens here!
         ulong c=rem.back();
	 @< Track maximal coefficient and show progress @>
         if (output) write_bytes(c, coefficient_size, out);
       }
       catch (bool)
       // incompatibility found during lift; details are already printed
       { std::cerr << "In coefficient " << d
		   << " of polynomial " << i << ".\n";
         exit(1);
       }
    }
  }
  if (verbose) // make final polynomial displayed look right, and stay visible
    std::cerr << "Polynomial: " << std::setw(10) << nr_pols-1 << '\n';
  return max;
}

@ The modular polynomial coefficients are extracted one polynomial at a time
via the |coefficients| method of the |mod_info| elements, and collected in the
vector |modular_pol| of polynomials (which is re-used repeatedly). Their
length need no be the same for all moduli, so we take for |len| the maximal
length.

@< Read polynomial coefficients... @>=
for (ulong j=0; j<mod_info.size(); ++j)
{ modular_pol[j]=mod_info[j]->coefficients(i);
  if (modular_pol[j].size()>len) len=modular_pol[j].size();
}

@ Keeping track of the maximal coefficient is trivial; each time it increases
we update the current display line. This line is rarely printed and rather
informative, so we print it even in quiet mode.

@< Track maximal coefficient and show progress @>=
if (c>max)
{ max=c;
  std::cout << "Maximal coefficient so far: " @|
            << max << ", in polynomial " << i << std::endl;
}


@* The main function.
%
Finally we must put everything together. What is left is mainly getting
arguments, and opening corresponding files, but we must not forget to create
our Chinese boxes! There will be one less of them than there are moduli, and
they are allocated via |new| to make a vector of pointers.

@h <sstream>
@< Main function @>=
int main(int argc, char** argv)
{ --argc; ++argv; // skip program name
  bool verbose=true, double_tables=false, output=true;
  if (argc>0 and std::string(*argv)=="-q") {@; verbose=false; --argc; ++argv;}
  if (argc>0 and std::string(*argv)=="-d")
    {@; double_tables=true; --argc; ++argv;}
  if (argc>0 and std::string(*argv)=="-nowrite")
    {@; output=false; --argc; ++argv;}
@)
  std::string mat_base,coef_base;
  // base names for renumbering and coefficient files
  std::vector<ulong> moduli;
  bool interactive= argc<2;

  @< Read in base names for files and |moduli| @>
@)
  ulong n=moduli.size();
  std::vector<ChineseBox*> box(n-1,NULL);
  ulong lcm;
  @< Install |n-1| Chinese boxes and set |lcm| @>
@)
  ulong coefficient_size=1, rem=lcm-1; // maximal remainder
  while ((rem>>=8)!=0) ++coefficient_size;
@)
  std::vector<modulus_info*> mod_info;
   // among other things this will hold the input files
  std::ofstream coefficient_file;
@)
  bool caught=false;
  try
  { if (interactive) test(moduli,box); // does not return

    @< Open input and output files @>

    file_pos nr_c=
      write_indices(coefficient_size,mod_info,coefficient_file,verbose,output);
    std::cout << "Done writing indices, will now write "
              << nr_c << " coefficient bytes." << std::endl;
    ulong max_coef=
       write_coefficients
         (coefficient_size,mod_info,box,coefficient_file,verbose,output);
    std::cout << "Maximal coefficient found: "
              << max_coef << "." << std::endl;
  }
  catch (...)
  { caught=true;
    // record that and exception was caught, and proceed to clean-up code
  }

@//* the following code aught to have been done in a destructor, but we found
     it to bothersome to define classes to pack |mod_info| and |box| into
  */
  for (ulong i=0; i<mod_info.size(); ++i) delete mod_info[i];
  for (ulong i=0; i<box.size(); ++i) delete box[i];
  if (caught) throw; // re-throw
}

@ When getting moduli, we check that they are indeed numeric and nonzero
(actually it suffices to start with a digit, trailing nonsense is ignored).

@< Read in base names for files and |moduli| @>=

if (interactive) argc=0; // ignore
else  {@; mat_base=*argv++; --argc; coef_base=*argv++; --argc; }
@)
if (argc==0) @< Get |moduli| from user input @>
else while (argc-->0)
{ std::istringstream in(*argv++);
  unsigned int m=0; in >> m;
  if (m!=0) moduli.push_back(m);
  else
  {@; std::cerr << "Illegal modulus argument: " << in.str() << "\n";
     exit(1);
  }
}
if (moduli.size()<1) {@; std::cerr << "Too few moduli.\n"; exit(1); }

@ We first install simple Chinese boxes, then if their |gcd| field is unequal
to~|1|, we replace it by a |TabledChineseBox|, or else by the more efficient
|PrimeTabledChineseBox|. The order of combining moduli used here matches the
one used in the loop of |write_coefficients|.

@< Install |n-1| Chinese boxes... @>=
{ std::vector<ulong> mod(moduli);
  for (ulong i=0; i<n-1; ++i)
  { ChineseBox b(mod[2*i],mod[2*i+1]);
    mod.push_back(b.get_lcm()); // defines |mod[n+i]|
    if (double_tables)
      box[i]= b.get_gcd()!=1 ? new DoubleTabledChineseBox(b)
			     : new PrimeDoubleTabledChineseBox(b);
    else
      box[i]= b.get_gcd()!=1 ? new TabledChineseBox(b)
			     : new PrimeTabledChineseBox(b);
  }
  lcm=mod.back();
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
    { std::cout << "Assuming canonical order for modulus " << moduli[i]
                << "." << std::endl;
      renumber_file=NULL; // signal absence of renumbering file
    }
@)
    name1 << coef_base << "-mod" << moduli[i];
    std::ifstream* coef_file=new std::ifstream(name1.str().c_str(),binary_in);
    if (not coef_file->is_open())
    @/{@; std::cerr << "Could not open file '" << name1.str() << "'.\n";
      exit(1);
    }
    @< Add a pointer to a new |modulus_info| object for |moduli[i]|,
      |renumber_file|, and |coef_file| to |mod_info|, or if possible a pointer
      to a new |modulus_info_with_table| object @>
  }
@)
  if (output)
  { std::ostringstream name;
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
  else std::cout << "Not writing output files" << std::endl;
}

@ The code below installs the objects that make a dynamic choice of the access
method of the renumbering file possible. If one succeeds in constructing
|modulus_info_with_table| objects, which involves allocating vectors for the
translation, then these will be used; if it fails because of lack of memory,
then we catch the |std::bad_alloc| that was thrown, trust that temporaries
have been cleaned up, and settle for installing a |modulus_info| object
constructed with the same parameters.

@< Add a pointer to a new |modulus_info| object... @>=
try
{ mod_info.push_back
   (new modulus_info_with_table(moduli[i],renumber_file,coef_file));
}
catch (std::length_error)
{ std::cout << "Renumbering table for modulus " << moduli[i] @|
            << " is too large for architecture, doing without" << std::endl;
@/mod_info.push_back(new modulus_info(moduli[i],renumber_file,coef_file));
}
catch (std::bad_alloc)
{ std::cout << "Not enough virtual memory for renumbering table for modulus "
	    << moduli[i] @| << ", doing without." << std::endl;
@/mod_info.push_back(new modulus_info(moduli[i],renumber_file,coef_file));
}


@ It might happen that the output modulus coincides with one of the input
moduli, for instance if one takes twice the same input modulus. In such cases
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
interactively. First the interactive input of moduli, which is also used
if the command line holds file names but no moduli.

@< Get |moduli| from user input @>=
{ std::cout << "Give moduli used, or 0 to terminate.\n";
  while(true)
  { std::cout << "Modulus: ";
    ulong m=0; std::cin >> m;
    if (m==0) break;
    moduli.push_back(m);
  }
}

@ The test function itself. One of the useful aspects of this test is that we
actually check that the lifting is correct, which check we left out for
efficiency reasons in the |write_coefficients| function; this is a quite
efficient detection of lists of moduli that run into integral overflow
problems.

@< Function definitions @>=
void test(std::vector<ulong>& moduli,std::vector<ChineseBox*>& box)
{ ulong n=moduli.size();
  ulong lcm=box.back()->get_lcm();
@)
  std::vector<ulong> remainder(2*n-1);
  // remainders for |n| original and |n-1| derived moduli

  while(true)
  { std::cout << "Give remainders mod " << moduli[0];
    for (ulong i=1; i<n; ++i) std::cout << ", " << moduli[i];
    std::cout << ": ";
    for (ulong i=0; i<n; ++i)
    { remainder[i]=~0ul; std::cin >> remainder[i];
      if (remainder[i]==~0ul)
      {@; std::cout << "Bye.\n"; exit(0); }
      remainder[i]%=moduli[i];
    }
@)
    try
    { for (ulong i=0; i<n-1; ++i)
      { remainder[n+i]=
          box[i]->lift_remainders(remainder[2*i],remainder[2*i+1]);
        std::cout << "Lifted to " << remainder[n+i]
		  << " (mod " << box[i]->get_lcm() << ").\n";
      }
      std::cout << "Solution found: " << remainder.back()
		<< " (mod " << lcm << ").\n";
      for (ulong i=0; i<n; ++i)
        if (remainder.back()%moduli[i]!=remainder[i])
	@/{@; std::cout << "Solution does not pass test.\n"; throw true; }
      std::cout << "Solution checks correctly.\n";
    }
    catch (bool b)
    { if (not b) std::cout << "No solution.\n"; }
  }

}


@* Index.

% Local IspellDict: british
