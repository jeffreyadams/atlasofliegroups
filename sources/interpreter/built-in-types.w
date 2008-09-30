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

\def\emph#1{{\it#1\/}}
\def\Z{{\bf Z}}
\def\lcm{\mathop{\rm lcm}}
\let\eps=\varepsilon

@* Built-in types.
%
This file describes several built-in types related to the Atlas software,
which used by the interpreter as primitive types. It also defines built-in
functions that relate to these types. This part of the program is completely
oriented towards the Atlas library, whereas the \.{evaluator} does things that
are for the most part independent of that library.

@h "built-in-types.h"

@c
namespace atlas { namespace interpreter {
namespace {@; @< Local function definitions @>@; }@;
@< Function definitions @>@;
}@; }@;

@ As usual the external interface is written to the header file associated to
this file.

@( built-in-types.h @>=

#ifndef BUILT_IN_TYPES_H
#define BUILT_IN_TYPES_H

@< Includes needed in the header file @>@;
namespace atlas { namespace interpreter {
@< Type definitions @>@;
@< Declarations of exported functions @>@;
}@; }@;
#endif

@ In execution, the main role of this compilation unit is to install the stuff
it defines into the tables. To that end we export its initialisation function.

@< Declarations of exported functions @>=
void initialise_builtin_types();

@ Since all wrapper functions are defined as local functions, their
definitions will have been seen when the definition below is compiled. This
avoids having to declare each wrapper function.

@< Function definitions @>=
void initialise_builtin_types()
@+{@; @< Install wrapper functions @>
}

@ Before we can define any types we must make sure the types like |value_base|
defined in \.{evaluator.w}, and like |lietype::LieType| defined
in~\.{lietype\_fwd.h}, are known in our header file.

@< Includes needed in the header file @>=
#include "evaluator.h"
#include "lietype_fwd.h"

@*1 Lie types.
Our first chapter concerns Lie types, as indicated by strings like
|"D4.A3.E8.T1"|.

@*2 The primitive type.
%
A first new type corresponds to the type |lietype::LieType| in the Atlas
software. We provide a constructor that incorporates a complete
|lietype::LieType| value, but often one will use the default constructor and
the |add| method. Otherwise most methods are obligatory, except the final two
which are provided for convenience.

@h <stdexcept>
@< Type definitions @>=
struct Lie_type_value : public value_base
{ lietype::LieType val;
@)
  Lie_type_value() : val(0) @+ {}
    // default constructor, produces empty type
  Lie_type_value(lietype::LieType t) : val(t) @+{}
    // constructor from already validated Lie type
@)
  virtual void print(std::ostream& out) const;
  Lie_type_value* clone() const @+{@; return new Lie_type_value(*this); }
  static const char* name() @+{@; return "Lie type"; }
@)
  void add_simple_factor (char,size_t)
    throw(std::bad_alloc,std::runtime_error); // grow
private:
  Lie_type_value(const Lie_type_value& v) : val(v.val) @+{}
    // copy constructor, used by |clone|
public:
  size_t rank() const;
  size_t semisimple_rank() const;
};
@)
typedef std::auto_ptr<Lie_type_value> Lie_type_ptr;

@ The type |lietype::LieType| is defined as |std::vector<SimpleLieType>| where
|SimpleLieType| stands for |std::pair<char,size_t>|. Therefore it
could take arbitrary values, not necessarily sensible ones. To remedy this we
make the method |add_simple_factor|, which is the only proposed way to build
up Lie types, check for the validity.

Since the tests defined in \.{io/interactive\_lietype.cpp} used in the
current interface for the Atlas software are clumsy to use, we perform
our own tests here, emulating |interactive_lietype::checkSimpleLieType|.
Torus factors of rank $r>1$ should be equivalent to $r$ torus factors of
rank~$1$, and it simplifies the software if we rewrite the former form to the
latter on input, so we do that here.

@h "lietype.h"
@h "constants.h"

@< Function definitions @>=
void Lie_type_value::add_simple_factor (char c,size_t rank)
   throw(std::bad_alloc,std::runtime_error)
{ using std::runtime_error;
  static const std::string types=atlas::lietype::typeLetters; // |"ABCDEFGT"|
  size_t t=types.find(c);
  if (t==std::string::npos)
    throw runtime_error(std::string("Invalid type letter '")+c+'\'');
  static const size_t lwb[]={1,2,2,4,6,4,2,0};
  static const size_t r=constants::RANK_MAX;
  static const size_t upb[]={r,r,r,r,8,4,2,r};
  if (rank<lwb[t])
    throw runtime_error("Too small rank "+num(rank)+" for Lie type "+c);
@.Too small rank@>
  if (rank>upb[t])
    if (upb[t]!=r)
      throw runtime_error("Too large rank "+num(rank)+" for Lie type "+c);
@.Too large rank@>
    else
      throw runtime_error
      ("Rank "+num(rank)+" exceeds implementation limit "+num(r));
@.Rank exceeds implementation limit@>
  if (c=='T')
    while (rank-->0) val.push_back(lietype::SimpleLieType('T',1));
  else
    val.push_back(lietype::SimpleLieType(c,rank));
}

@ Now we define a wrapper function that really builds a |Lie_type_value|. We
scan the string looking for sequences of a letter followed by a number. We
allow and ignore sequences of punctuation characters between the two. Since
the correct structure of Lie type strings is so obvious to the human eye, our
error message just cites the entire offending string, rather than trying to
point out the error exactly.

@h <cctype>
@< Local function definitions @>=
inline void skip_punctuation(const char* &p)
{@; while (std::ispunct(*p) || std::isspace(*p)) ++p;}
@)
void Lie_type_wrapper() throw(std::bad_alloc,std::runtime_error)
{ string_ptr s(get<string_value>());
  Lie_type_ptr result(new Lie_type_value);
  size_t total_rank=0;
@/const char* p=s->val.c_str(); skip_punctuation(p);
  while (std::isalpha(*p))
  { char c=*p++;
    if (!std::isdigit(*p)) {@; --p; break; } // and |throw| a |runtime_error|
    size_t rank=*p++-'0'; @+
      while (std::isdigit(*p)) rank=10*rank+(*p++-'0');
    skip_punctuation(p);
    result->add_simple_factor(c,rank);
      // this may |throw| a |runtime_error| as well
    if ((total_rank+=rank)>constants::RANK_MAX)
      throw std::runtime_error
      ("Total rank exceeds implementation limit "+num(constants::RANK_MAX));
@.Total rank exceeds...@>
  }
  if (*p!='\0')
    throw std::runtime_error
      ("Error in type string '"+s->val+"' for Lie type");
  push_value(result);
}

@ We shall call this function \.{Lie\_type}.

@< Install wrapper functions @>=
install_function(Lie_type_wrapper,"Lie_type","(string->LieType)");

@*2 Printing Lie types.
Before we do anything more complicated with this primitive type, we must
ensure that we can print its values. We can use an operator defined in
\.{basic\_io.cpp}.

@h "basic_io.h"
@< Function definitions @>=
void Lie_type_value::print(std::ostream& out) const
{ using atlas::operator<<;
  if (val.empty()) out << "empty Lie type";
  else
    out << "Lie type '" << val << '\'';
}

@*2 Auxiliary functions for Lie types.
As a service for testing routines we provide members that compute the rank
and the semisimple rank.

@< Function definitions @>=
size_t Lie_type_value::rank() const
{ size_t r=0;
  for (size_t i=0; i<val.size(); ++i) r+=val[i].second;
  return r;
}
@)
size_t Lie_type_value::semisimple_rank() const
{ size_t r=0;
  for (size_t i=0; i<val.size(); ++i)
    if (val[i].first!='T') r+=val[i].second;
  return r;
}

@ Here is a function that computes the Cartan matrix for a given Lie type.

@h "prerootdata.h"
@< Local function definitions @>=
void Cartan_matrix_wrapper()
{ Lie_type_ptr t(get<Lie_type_value>());
  matrix_ptr
    result(new matrix_value(latticetypes::LatticeMatrix()));
  prerootdata::cartanMatrix(result->val,t->val);
  push_value(result);
}


@ And here is a function that tries to do the inverse. We do no tests, so one
can play around and see what |dynkin::lieType| does with fairly random
matrices. For instance any null columns in the matrix will be mistaken for
factors of type $A_1$.

@h "dynkin.h"
@< Local function definitions @>=
void type_of_Cartan_matrix_wrapper ()
{ matrix_ptr m(get<matrix_value>());
  lietype::LieType t;
  dynkin::lieType(t,m->val);
  push_value(new Lie_type_value(t));
}
@~Again we install our wrapper functions.
@< Install wrapper functions @>=
install_function(Cartan_matrix_wrapper,"Cartan_matrix","(LieType->mat)");
install_function(type_of_Cartan_matrix_wrapper
		,"type_of_Cartan_matrix","(mat->LieType)");

@*2 Finding lattices for a given Lie type.
%
When a Lie type is fixed, there is still a nontrivial choice to determine the
root datum for a connected Complex reductive group of that type: one has to
choose a sublattice of the weight lattice of the ``simply connected'' group of
that type to become the weight lattice of the chosen Complex reductive group
(a finite quotient of the simply connected group); that sublattice should have
full rank and contain the weight lattice. We shall start with considering some
auxiliary functions to help facilitate this choice.

The function |Smith_basis| below is copied from the Atlas software sources,
where the function is originally defined in the file
\.{io/interactive\_lattice.cpp}. We chose not to link to the corresponding
object file, since that would require including as well a major part of the
Atlas library, some of which we don't need at all, and other parts that at
least are not needed at this point. In fact |Smith_basis| itself requires
nothing but access to the Smith normal form algorithm, which being a template
can be done without linking against any particular object file. Might we
suggest that this function that performs no interaction would have been better
placed elsewhere than in the \.{io} sub-directory?

The purpose of |Smith_basis| is to prepare a basis of the weight lattice of
the simply connected group in terms of which the choice for a sublattice can
be clearly presented; the important property is that there is a subset of
basis vectors (of size the semisimple rank) such that the root lattice is
spanned by pure multiples of these vectors. We quote the remainder of this
comment from its original documentation.

The Lie type is defined as a sequence of simple and torus factors inside |lt|.
This function constructs a ``blockwise Smith normal'' basis for the root
lattice inside the weight lattice (i.e., it does just that for each semisimple
block, and returns the canonical basis for the torus blocks).

The purpose of doing this blockwise instead of globally is to permit a
better reading of the quotient group: this will be presented as a sequence
of factors, corresponding to each simple block.

@h "smithnormal.h"

@< Local function definitions @>=
void smithBasis(latticetypes::CoeffList& invf, latticetypes::WeightList& b,
		const lietype::LieType& lt)

{ matrix::initBasis(b,lietype::rank(lt));
  latticetypes::WeightList::iterator bp = b.begin();

  for (size_t j = 0; j < lt.size(); ++j)
  // Smith-normalize for each simple factor
  { size_t r = lietype::rank(lt[j]);

    if (lietype::type(lt[j]) == 'T')
    // torus type: insert $r$ null ``invariant factors''
      invf.insert(invf.end(),r,latticetypes::ZeroCoeff);
    else {
    latticetypes::LatticeMatrix ms;
    prerootdata::cartanMatrix(ms,lt[j]);
    ms.transpose();
    smithnormal::smithNormal(invf,bp,ms);

    if (lietype::type(lt[j]) == 'D' and r%2 == 0)
      // adjust generators back to canonical basis
      bp[r-2]+=bp[r-1];
    }
    bp += r;
  }
}

@ And here is a wrapper function that applies the previous function, and
returns a block-wise Smith basis for the transposed Cartan matrix for a given
Lie type, and the corresponding invariant factors. In case of torus factors
this description should be interpreted in the sense that the Smith basis for
the corresponding block is the standard basis and the invariant factors are
null.

@< Local function definitions @>=
void Smith_Cartan_wrapper()
{ Lie_type_ptr t(get<Lie_type_value>());
  vector_ptr inv_factors (new vector_value(latticetypes::CoeffList(0)));
  latticetypes::WeightList b; @+
  smithBasis(inv_factors->val,b,t->val);
  matrix_ptr m @| (new matrix_value(latticetypes::LatticeMatrix(b)));
    // convert from |WeightList|
  push_value(m); push_value(inv_factors); wrap_tuple(2);
}

@ The result of |Smith_Cartan| can serve among other things to help specifying
a basis for subgroup of the quotient of the original weight
lattice~$\tilde{X}$ by the sub-lattice spanned by the roots (the columns of
the transposed Cartan matrix). For that purpose only invariant factors unequal
to~$1$ and the corresponding columns of the Smith basis are of interest.
Therefore the following wrapper function, which may operate directly on the
result of the previous one, filters out the invariant factors~$1$ and the
corresponding columns.

This function is particular in that it returns exactly the kind of arguments
that it requires. Therefore rather than calling |push_tuple_components| which
would pop and destroy the tuple, we just extract its fields while leaving it
on the stack; its fields are modified in place, so there is no |push_value| to
return the value either. Therefore we do not own the arguments during the
call, and auto-pointers should not be used; incidentally |erase| and
|eraseColumn| should have no reason whatsoever to throw an exception.

@< Local function definitions @>=
void filter_units_wrapper ()
{ tuple_value* t=force<tuple_value>(execution_stack.back());
  matrix_value* basis=force<matrix_value>(t->val[0]);
  vector_value* inv_f=force<vector_value>(t->val[1]);
  if (inv_f->val.size()!=basis->val.numColumns())
    throw std::runtime_error("Size mismatch "+
      num(inv_f->val.size())+':'+num(basis->val.numColumns())+
      " in filter_units");
@)
  size_t i=0;
  while (i<inv_f->val.size())
    if (inv_f->val[i]!=1) ++i; // keep invariant factor and column
    else
    {@; inv_f->val.erase(inv_f->val.begin()+i);
        basis->val.eraseColumn(i);
    }
}

@ Here is another function, copied from |makeOrthogonal| in
\.{io/interactive\_lattice.cpp} for essentially the same reasons as
|smithBasis| above was. We take the occasion however to change its name and
interface a bit, avoiding the need to introduce rational
vectors and matrices as primitive type, and we also adapt the commentary. In
fact we completely changed this function beyond recognition.

The function |annihilator_modulo| takes as argument an $m\times{n}$
matrix~$M$, and an integer |d|. It returns a $m\times{m}$ matrix~|A|
whose columns span the sub-lattice of $\Z^m$ of vectors $v$ such that
$v^t\cdot{M}\in d\,\Z^n$. The fact that $d$ is called |denominator| below comes
from the alternative interpretation that the transposes of the columns of~$A$
applied to the rational matrix $M/d$ give integral vectors.

The algorithm is quite simple. After finding a matrix~|A| describing the Smith
basis~$(b_j)_{j\in[m]}$ for the lattice spanned by the columns of~|M|, and
invariant factors~$(\lambda_j)_{j\in[r]}$ (where $r$ is its rank, and $[r]$
abbreviates $\{0,\ldots,r-1\}$), we compute the dual basis $(b^*_j)_{j\in[m]}$
as the columns of the transpose inverse matrix of~|A|. Then we know that
$b^*_j$ applied to any column of~$M$ lies in $\lambda_j\Z$ (where we take
$\lambda_j=0$ for $j\notin[r]$), and our result is obtained by multiplying
each column~$j\in[r]$ of the matrix giving the dual basis by
$\lcm(d,\lambda_j)/\lambda_j=d/\gcd(d,\lambda_j)$.

@f lambda NULL

@h "arithmetic.h"
@h "lattice.h"
@h "latticetypes.h"
@h "matrix.h"

@< Local function definitions @>=
latticetypes::LatticeMatrix @|
annihilator_modulo
(const latticetypes::LatticeMatrix& M,
 latticetypes::LatticeCoeff denominator)

{ const size_t m=M.numRows(); // determines dimension of output

  latticetypes::WeightList b; matrix::initBasis(b,m);
    // standard basis of rank $m$

  latticetypes::CoeffList lambda;
  smithnormal::smithNormal(lambda,b.begin(),M);
    // find Smith basis

  latticetypes::LatticeMatrix A
    (latticetypes::LatticeMatrix(b).inverse().transposed());

  for (size_t j = 0; j < lambda.size(); ++j)
  { unsigned long f=(lambda[j]);
    unsigned long c=denominator/arithmetic::gcd(f,denominator);
    for (size_t i=0; i<m; ++i) A(i,j)*=c; // multiply column by |c|
  }
  return A;
}

@ The wrapper function is particularly simple. In fact using auto-pointers
seems overly prudent, but in fact |annihilator_modulo| could cause an
exception during matrix inversion or transposition. We need to put the result
from |annihilator_modulo| into a named variable since |swap| requires a
modifiable reference.

@< Local function definitions @>=
void ann_mod_wrapper()
{ push_tuple_components();
  int_ptr d(get<int_value>());
  matrix_ptr m(get<matrix_value>());
@)
  latticetypes::LatticeMatrix A=
    annihilator_modulo(m->val,d->val);
  m->val.swap(A);
  push_value(m);
}

@ Next a simple administrative routine, needed here because we cannot handle
matrices in our programming language yet. Once one has computed a new lattice
(possibly with the help of |ann_mod|), in the form of vectors to replace those
selected by |filter_units| from the result of |Smith_Cartan|, one needs to
make the replacement. The following function does this, taking its first two
arguments as the result of |Smith_Cartan|, and the third a matrix whose
columns are to be substituted. The second argument serves only to determine,
by the place of its non-unit entries, where the insertion has to take place.
In fact this is so simple that we define the wrapper function directly. And in
fact it will be more practical to take as argument a pair of a pair as
returned by |Smith_Cartan| and a matrix.

@< Local function definitions @>=
void replace_gen_wrapper ()
{ push_tuple_components(); // a pair
  matrix_ptr new_generators(get<matrix_value>());
  push_tuple_components(); // a pair as returned by \.{Smith\_Cartan}
  vector_ptr inv_f(get<vector_value>());
  matrix_ptr old_generators(get<matrix_value>());
@)
  if (new_generators->val.numRows()!=old_generators->val.numRows())
    throw std::runtime_error("Column lengths do not match in replace_gen");
@.Column lengths do not match@>
  if (inv_f->val.size()!=old_generators->val.numColumns())
    throw std::runtime_error("Number of columns mismatch in replace_gen");
@.Size mismatch in replace\_gen@>
@)
  for (size_t j=0,k=0; j<inv_f->val.size(); ++j)
    if (inv_f->val[j]!=1)
       // replace column |j| by column |k| from |new_generators|
    { if (k>=new_generators->val.numColumns())
        throw std::runtime_error
          ("Not enough replacement columns in replace_gen");
@.Not enough replacement columns@>
      for (size_t i=0; i<old_generators->val.numRows(); ++i)
        old_generators->val(i,j)=new_generators->val(i,k);
      ++k;
    }
  push_value(old_generators);
}

@*2 Specifying inner classes. Now we move ahead a bit in the theory, from
functions that help in building root data to functions that help defining
(inner classes of) real forms. The first of such functions is
|lietype::involution|, which takes a Lie type and a |lietype::InnerClassType|
(a vector of characters describing the kind of involution wanted) and produces
a matrix describing the involution, defined on the weight lattice for the
simply connected group of the given type. That function supposes its arguments
have already been checked for validity and undergone some transformation; in
the existing Atlas interface this was done by
|interactive_lietype::checkInnerClass| and
|interactive_lietype::readInnerClass|. This forces us to perform similar
actions before calling |lietype::involution|. We prefer not to use the
functions defined in \.{io/interactive\_lietype}, for the same reason we did
not use |interactive_lietype::checkSimpleLieType| above. Therefore we shall
first define a testing/transformation routine that takes a string describing
an inner class, and transforms it into |lietype::InnerClassType| that is
guaranteed to be valid if returned; the routine throws a |runtime_error| in
case of problems.

@< Local function definitions @>=
lietype::InnerClassType transform_inner_class_type
  (const char* s, const lietype::LieType& lt)
throw (std::bad_alloc, std::runtime_error)
{ static const std::string types=atlas::lietype::innerClassLetters;
    // |"Ccesu"|
  lietype::InnerClassType result(0);
  size_t i=0; // position in simple factors of Lie type |lt|
  for (;skip_punctuation(s),*s!='\0'; ++s)
    @< Test the inner class letter |*s|, and either push a corresponding type
       letter onto |result| while advancing~|i| by the appropriate amount, or
       throw a |runtime_error| @>
  if (i< lt.size()) throw std::runtime_error("Too few inner class symbols");
  return result;
}

@ The type letter |'C'| meaning ``Complex'' is the only one using more that
one Lie type, and both types used must be equal. Types |'s'| ``split'' and
|'c'| ``compact'' with its synonym |'e'| ``equal rank'' require no particular
conditions or treatment; the remaining case is relegated to the next section.

@< Test the inner class letter... @>=
{ if (types.find(*s) == std::string::npos) throw std::runtime_error@|
    (std::string("Unknown inner class symbol `")+*s+"'");
  if (i>= lt.size()) throw std::runtime_error("Too many inner class symbols");
  if (*s=='C') // complex inner class, which requires two equal simple factors
  { if (i+1>=lt.size() or lt[i+1]!=lt[i]) throw std::runtime_error @|
      ("Complex inner class needs two identical consecutive types");
    result.push_back('C');
    i+=2; // advance past both simple factors
  }
  else if (*s=='s') {@; result.push_back('s'); ++i; }
  else if (*s=='c' or *s=='e')
  {@; result.push_back('c'); ++i; } // synonyms
  else if (*s=='u') // unequal rank, the only remaining possibility
  { @< Test and possibly transform the unequal rank inner class symbol @>
    ++i;
  }
}

@ The type letter |'u'| meaning ``unequal rank'' is mathematically only
meaningful for Lie types $A_n$ with $n\geq2$, any legal type $D_n$, or $E_6$.
Moreover in all cases except $D_n$ with $n$ even, this designation is
equivalent to |'s'|, and will be replaced by it. Hence any surviving type
letter |'u'| corresponds to a type of the form~$D_{2n}$.

@< Test and possibly transform... @>=
{ lietype::TypeLetter t = lietype::type(lt[i]); size_t r=lietype::rank(lt[i]);
  if (t=='D') result.push_back(r%2==0 ? 'u' : 's');
  else if (t=='A' and r>=2 or t=='E' and r==6) result.push_back('s');
  else throw std::runtime_error @|
    (std::string("unequal rank class is meaningless for type ")+t+num(r));
}

@ The wrapper function around |lietype::involution| will take a Lie type and a
string of type letters and return a matrix describing the involution
designated by that string, expressed on the fundamental weight basis for the
simply connected group of that type. A variant of this function will later be
defined giving the same matrix on the basis fundamental weights for a complex
group given by a root datum.

@< Local function def... @>=
void basic_involution_wrapper()
{ push_tuple_components();
@/string_ptr str(get<string_value>());
@/Lie_type_ptr t(get<Lie_type_value>());
@)matrix_ptr m
     (new matrix_value(latticetypes::LatticeMatrix()));
@/lietype::involution(m->val,t->val
                     ,transform_inner_class_type(str->val.c_str(),t->val));
  push_value(m);
}

@ The function just defined gives an involution on the basis of fundamental
weights for a simply connected group, or on the basis of the root lattice for
an adjoint group. For a general complex groups however, we need to transform
this involution to an involution of a given sub-lattice, assuming that this is
possible, in other words if that the sub-lattice is globally stable under the
involution specified. Therefore we now provide a function that in addition to
the Lie type takes a matrix specifying a sub-lattice as argument. These two
ingredients will also be used below to construct a root datum; it would be
nice to use the resulting root datum as argument here, but unfortunately the
Lie type and sub-lattice cannot be recovered from it with certainty, and they
are needed separately in the construction of an inner class. Note however that
when constructing inner classes below, we shall nevertheless provide a
function $set\_inner\_class$ that operates on a root datum and a string
specifying an inner class; it will deduce candidates for the Lie type and the
sub-lattice from the root datum and then call the function $based\_involution$
defined here.


@< Local function def... @>=
void based_involution_wrapper()
{ push_tuple_components();
@/string_ptr str(get<string_value>());
@/matrix_ptr basis(get<matrix_value>());
@/Lie_type_ptr type(get<Lie_type_value>());
@)
  size_t r=type->rank();
  if (basis->val.numRows()!=r or basis->val.numRows()!=r)
    throw std::runtime_error @|
    ("lattice matrix should be "+num(r)+'x'+num(r)+ " for this type");
  matrix_ptr m
     (new matrix_value(latticetypes::LatticeMatrix()));
@/lietype::involution
    (m->val,type->val
    ,transform_inner_class_type(str->val.c_str(),type->val));
@)
  latticetypes::LatticeCoeff d;
  m->val = basis->val.inverse(d) * m->val * basis->val;
  if (d==0 or !m->val.divisible(d)) throw std::runtime_error
    ("inner class is not compatible with given lattice");
  m->val/=d; push_value(m);
}


@ All that remains is installing the wrapper functions.

@< Install wrapper functions @>=
install_function(Smith_Cartan_wrapper,"Smith_Cartan","(LieType->mat,vec)");
install_function(filter_units_wrapper,"filter_units","(mat,vec->mat,vec)");
install_function(ann_mod_wrapper,"ann_mod","(mat,int->mat)");
install_function(replace_gen_wrapper,"replace_gen",
		"((mat,vec),mat->mat)");
install_function(basic_involution_wrapper,"basic_involution",
		"(LieType,string->mat)");
install_function(based_involution_wrapper,"based_involution",
		"(LieType,mat,string->mat)");

@*1 Root data.
Now we are ready to introduce a new primitive type for root data.

@< Includes needed in the header file @>=
#include "rootdata.h"

@ The root datum type is laid out just like previous primitive types are.

@< Type definitions @>=
struct root_datum_value : public value_base
{ rootdata::RootDatum val;
@)
  root_datum_value(const rootdata::RootDatum v) : val(v) @+ {}
  virtual void print(std::ostream& out) const;
  root_datum_value* clone() const @+{@; return new root_datum_value(*this); }
  static const char* name() @+{@; return "root datum"; }
private:
  root_datum_value(const root_datum_value& v) : val(v.val) @+{}
};
@)
typedef std::auto_ptr<root_datum_value> root_datum_ptr;

@*2 Printing root data. We shall not print the complete information contained
in the root datum. However we do exercise the routines in \.{dynkin.cpp} to
determine the type of the root datum from its Cartan matrix, as is done in the
function |type_of_Cartan_matrix| above. Contrary to
|prerootdata::cartanMatrix|, the function |rootdata::cartanMatrix| produces
nothing for torus factors, so we add any necessary torus factors at the end.

@< Local fun... @>=
lietype::LieType type_of_datum(const rootdata::RootDatum& rd)
{ latticetypes::LatticeMatrix Cartan; rootdata::cartanMatrix(Cartan,rd);
@/lietype::LieType t; dynkin::lieType(t,Cartan);
  if (!rd.isSemisimple())
    for (size_t i=rd.semisimpleRank(); i<rd.rank(); ++i)
      t.push_back(lietype::SimpleLieType('T',1));
  return t;
}

@ Then this is how we print root data. The mention of |"simply connected"|
and/or |"adjoint"| is in order to make some use of these flags, which are
present anyway.

@< Function definitions @>=
void root_datum_value::print(std::ostream& out) const
{ Lie_type_value type(type_of_datum(val));
  out << (val.isSimplyConnected() ? "simply connected " : "")
   @| << (val.isAdjoint() ? "adjoint " : "")
      << "root datum of " << type;
}

@ We also make the derivation of the type available by a wrapper function.
@< Local fun...@>=
void type_of_root_datum_wrapper()
{ root_datum_ptr rd(get<root_datum_value>());
  push_value(new Lie_type_value(type_of_datum(rd->val)));
}

@*2 Building a root datum.
%
To create a root datum value, the user must specify a Lie type and a square
matrix of the size of the rank of the root datum, which specifies generators
of the desired weight lattice as a sub-lattice of the lattice of weights
associated to the simply connected group of the type given. The given weights
should be independent and span at least the root lattice associated to the
type.

@< Local function definitions @>=
void root_datum_wrapper()
{ push_tuple_components();
  matrix_ptr lattice(get<matrix_value>());
  Lie_type_ptr type(get<Lie_type_value>());
  if (lattice->val.numRows()!=lattice->val.numColumns() @| or
      lattice->val.numRows()!=type->rank())
    throw std::runtime_error
    ("lattice matrix should be " @| +num(type->rank())+'x'+num(type->rank())+
     " for this type");
  @< Test whether the columns of |lattice->val| span the root lattice; if
     not |throw| a |std::runtime_error| @>
  latticetypes::WeightList columns; columnVectors(columns,lattice->val);
  prerootdata::PreRootDatum pre_datum(type->val,columns);
  push_value(new root_datum_value(rootdata::RootDatum(pre_datum)));
}

@ Since the construction of the (simple) roots in a |PreRootDatum| uses
multiplication by an inverse matrix producing an integral matrix, without
testing whether the result was exact, we must laboriously do the same
computation before constructing the |PreRootDatum|, just to assure that the
construction will not fail. If the construction would succeed, even with an
incorrect result, we could test validity more easily be multiplying back the
given lattice basis with the roots in the |PreRootDatum|, and testing for
equality with the transpose Cartan matrix. But this would crash on division be
zero if the matrix were singular (this raises no exception, despite printing
message mentioning a \.{floating point exception}, so it cannot be caught). So
that approach would still require preliminary test for an invertible matrix,
and would end up being as much work as the code below. We note that moreover
in that approach (which we had initially adopted) a subtle error is possible,
because roots that after basis transformation are rounded to zero would be
mistaken for columns coming from torus factors, and the |PreRootDatum| could
therefore have fewer roots than the semisimple rank of the type specified. But
in case our lattice matrix passes the test below, the |PreRootDatum| can
safely be constructed and will be correct in all respects.

@< Test whether the columns of |lattice->val| span the root lattice... @>=
{ latticetypes::LatticeCoeff d;
  latticetypes::LatticeMatrix M=lattice->val.inverse(d);
  if (d==0) throw std::runtime_error
    ("Lattice matrix has dependent columns; in root_datum");
@/latticetypes::LatticeMatrix tC;
  prerootdata::cartanMatrix(tC,type->val); tC.transpose();
  M *= tC;
  // now columns of |M| hold |d| times the simple roots, expressed on the given basis
  if (!M.divisible(d))
    throw std::runtime_error@|
      ("Given lattice does not contain the root lattice; in root_datum");
}

@ To emulate what is done in the Atlas software, we write a function that
integrates some of the previous ones. It is called as $quotient\_basis(t,M)$
where $t$ is a Lie type, and $M$ is a matrix whose columns represent the
denominators of the kernel generators as would be entered in the Atlas
software (the numerators are determined by~$t$ for non-torus factors, and for
the torus factor we take the least common multiple of those other
denominators, which should be large enough for practical purposes).

Let $S=Smith\_Cartan(t)$ and $(C,v)=filter\_units(S)$, then we find a
basis for the sub-lattice needed to build the root datum as follows. The
number of rows of~$M$ should match the length of~$v$. We compute the least
common multiple~$d$ of the nonzero values~$v_j$ (which is~$1$ if $v$ should be
completely null). The entries of~$M$ are brought to a common denominator~$d$
by multiplying each row~$i$ with $v_i\neq0$ by $d/v_i$. If the resulting
matrix is~$M'$, the call $quotient\_basis(t,M)$ amounts to evaluating
$replace\_gen(S,C*ann\_mod(M',d))$.

This wrapper function is particular in that it does most of its work by
calling other wrapper functions; we essentially use the evaluator stack here
many times. In one case we need to duplicate the top of the stack, holding the
value~$S$ (this is somewhat of a coincidence, since the bottom copy happens to
be in place for the call to~$replace\_gen$). We do this by pushing a cloned
version of the stack top.

@< Local function definitions @>=
void quotient_basis_wrapper()
{ push_tuple_components();
  matrix_ptr M(get<matrix_value>());
  // leave Lie type on stack for $Smith\_Cartan$
  Smith_Cartan_wrapper(); // compute $S=Smith\_Cartan(t)$
  value S=execution_stack.back();
    // leave |S| on stack for call of $replace\_gen$
  push_value(S->clone()); // push a copy for call of $filter\_units$
  filter_units_wrapper(); // compute |(C,v)|
  push_tuple_components();
  vector_ptr v(get<vector_value>());
  matrix_ptr C(get<matrix_value>());
@/size_t d=1;
  for (size_t i=0; i<v->val.size(); ++i)
    if (v->val[i]!=0) d=arithmetic::lcm(d,v->val[i]);
  @< Test |M| against |v| and adapt its rows to the common denominator |d| @>
  push_value(C); // for call of $mm\_prod$
  push_value(M); // for call of $ann\_mod$
  push_value(new int_value(d)); // for call of $ann\_mod$
@/wrap_tuple(2); ann_mod_wrapper();
@/wrap_tuple(2); mm_prod_wrapper();
@/wrap_tuple(2); replace_gen_wrapper();
}

@ As said above, |M| must have as many rows as |v| has entries. But as a
service to the user we adapt its size to match that requirement if it has no
entries (either no rows or no columns).

@< Test |M| against |v| and adapt its rows to the common denominator |d| @>=
{ if (M->val.isEmpty()) M->val.resize(v->val.size(),0);
  if (M->val.numRows()!=v->val.size())
    throw std::runtime_error @| ("Number "+num(M->val.numRows())+
      " of rows does not match number "+num(v->val.size())+
      " of kernel generators");
  for (size_t i=0; i<v->val.size(); ++i)
   if (v->val[i]!=0)
     for (size_t j=0; j<M->val.numColumns(); ++j)
       M->val(i,j)*=d/v->val[i];
}

@ The function that integrates all is $quotient\_datum$; the call
$quotient\_datum(t,M)$ is equivalent to $root\_datum(t,quotient\_basis(t,M))$.

@< Local function definitions @>=
void quotient_datum_wrapper()
{ tuple_ptr args(get<tuple_value>());
  push_value(args->val[0]->clone()); // the Lie type, for call of $root\_datum$
  push_value(args); quotient_basis_wrapper();
@/wrap_tuple(2); root_datum_wrapper();
}

@ We define two more wrappers with only a Lie type as argument, for building
the simply connected and the adjoint root data. They are similar to the
previous one in that they mostly call other wrapper functions: the matrices
that specify the sub-lattices are produced by the wrapper functions for
|id_mat| respectively for |Cartan_matrix| and |transpose_matrix|. The only
thing we do ``by hand'' is to make sure, in the case of the adjoint datum,
that all null diagonal entries are replaced by ones.

@< Local function definitions @>=
void simply_connected_datum_wrapper()
{ Lie_type_value* type=get<Lie_type_value>(); push_value(type);
  push_value(new int_value(type->rank()));
  id_mat_wrapper();
@/wrap_tuple(2); root_datum_wrapper();
}
@)
void adjoint_datum_wrapper()
{ Lie_type_value* type=get<Lie_type_value>(); push_value(type);
  push_value(type->clone());
  Cartan_matrix_wrapper(); transpose_mat_wrapper();
  matrix_value* M=get<matrix_value>();
  for (size_t i=0; i<type->rank(); ++i)
    if (M->val(i,i)==0) M->val(i,i)=1;
  push_value(M);
@/wrap_tuple(2); root_datum_wrapper();
}

@ Finally here are two more wrappers to make the root data for the special and
general linear groups, in a form that is the most natural. We supply a matrix
whose columns represent $\eps_1,\ldots,\eps_{n-1}$ for ${\bf SL}_n$ and
$\eps_1,\ldots,\eps_n$ for ${\bf GL}_n$, when the basis of the simply
connected group of type $A_{n-1}$ is given by $\omega_1,\ldots,\omega_{n-1}$
and that of $A_{n-1}T_1$ is given by $\omega_1,\ldots,\omega_{n-1},t$, where
$\omega_k=\sum_{i=1}^k\eps_i-kt$ and $\sum_{i=1}^n\eps_i=nt$. This matrix is
most easily described by an example, for instance for ${\bf GL}_4$ it is
$$
   \pmatrix{1&-1&0&0\cr0&1&-1&0\cr0&0&1&-1\cr1&1&1&1\cr}.
$$
and for ${\bf SL}_n$ it is the $(n-1)\times(n-1)$ top left submatrix of that
for ${\bf GL}_n$.

@< Local function definitions @>=
void SL_wrapper()
{ int_ptr n(get<int_value>());
  if (n->val<1) throw std::runtime_error("Non positive argument for SL");
  const size_t r=n->val-1;
  Lie_type_ptr type(new Lie_type_value());
  if (r>0) type->add_simple_factor('A',r);
  push_value(type);
  matrix_ptr lattice
     (new matrix_value(latticetypes::LatticeMatrix()));
  identityMatrix(lattice->val,r);
  for (size_t i=0; i<r-1; ++i) lattice->val(i,i+1)=-1;
  push_value(lattice);
@/wrap_tuple(2); root_datum_wrapper();
}
@)
void GL_wrapper()
{ int_ptr n(get<int_value>());
  if (n->val<1) throw std::runtime_error("Non positive argument for GL");
  const size_t r=n->val-1;
  Lie_type_ptr type(new Lie_type_value());
  if (r>0) type->add_simple_factor('A',r);
  type->add_simple_factor('T',1);
  push_value(type);
  matrix_ptr lattice
     (new matrix_value(latticetypes::LatticeMatrix()));
  identityMatrix(lattice->val,r+1);
  for (size_t i=0; i<r; ++i)
  @/{@; lattice->val(n->val-1,i)=1; lattice->val(i,i+1)=-1; }
  push_value(lattice);
@/wrap_tuple(2); root_datum_wrapper();
}

@*2 Functions operating on root data.
%
The following functions allow us to look at the roots and co-roots stored in
a root datum value, and the associated Cartan matrix.

@< Local function definitions @>=
void simple_roots_wrapper()
{ root_datum_ptr rd(get<root_datum_value>());
  latticetypes::WeightList l
    (rd->val.beginSimpleRoot(),rd->val.endSimpleRoot());
  push_value(new matrix_value(atlas::latticetypes::LatticeMatrix(l)));
}
@)
void simple_coroots_wrapper()
{ root_datum_ptr rd(get<root_datum_value>());
  latticetypes::WeightList l
    (rd->val.beginSimpleCoroot(),rd->val.endSimpleCoroot());
  push_value(new matrix_value(atlas::latticetypes::LatticeMatrix(l)));
}
@)
void datum_Cartan_wrapper()
{ root_datum_ptr rd(get<root_datum_value>());
  latticetypes::LatticeMatrix M;
  rootdata::cartanMatrix(M,rd->val);
  push_value(new matrix_value(M));
}

@ The following functions allow us to look at the roots and co-roots stored in
a root datum value.

@< Local function definitions @>=
void roots_wrapper()
{ root_datum_ptr rd(get<root_datum_value>());
  latticetypes::WeightList l
    (rd->val.beginRoot(),rd->val.endRoot());
  push_value(new matrix_value(atlas::latticetypes::LatticeMatrix(l)));
}
@)
void coroots_wrapper()
{ root_datum_ptr rd(get<root_datum_value>());
  latticetypes::WeightList l
    (rd->val.beginCoroot(),rd->val.endCoroot());
  push_value(new matrix_value(atlas::latticetypes::LatticeMatrix(l)));
}

@ It is useful to have bases for the sum of the root lattice and the
coradical, and for the sum of the coroot lattice and the radical.

@< Local function definitions @>=
void root_coradical_wrapper()
{ root_datum_ptr rd(get<root_datum_value>());
  latticetypes::WeightList l
    (rd->val.beginSimpleRoot(),rd->val.endSimpleRoot());
  l.insert(l.end(),rd->val.beginCoradical(),rd->val.endCoradical());
  push_value(new matrix_value(atlas::latticetypes::LatticeMatrix(l)));
}
@)
void coroot_radical_wrapper()
{ root_datum_ptr rd(get<root_datum_value>());
  latticetypes::WeightList l
    (rd->val.beginSimpleCoroot(),rd->val.endSimpleCoroot());
  l.insert(l.end(),rd->val.beginRadical(),rd->val.endRadical());
  push_value(new matrix_value(atlas::latticetypes::LatticeMatrix(l)));
}

@ And here is a simple function to dualise a root datum.

@h "tags.h"
@< Local function definitions @>=
void dual_datum_wrapper()
{ root_datum_ptr rd(get<root_datum_value>());
  rootdata::RootDatum dual(rd->val,tags::DualTag());
  rd->val.swap(dual);
  push_value(rd);
}

@ Let us install the above wrapper functions.

@< Install wrapper functions @>=
install_function(type_of_root_datum_wrapper,"type_of_root_datum"
                ,"(RootDatum->LieType)");
install_function(root_datum_wrapper,"root_datum","(LieType,mat->RootDatum)");
install_function(quotient_basis_wrapper
		,"quotient_basis","(LieType,mat->mat)");
install_function(quotient_datum_wrapper
		,"quotient_datum","(LieType,mat->RootDatum)");
install_function(simply_connected_datum_wrapper
		,"simply_connected_datum","(LieType->RootDatum)");
install_function(adjoint_datum_wrapper,"adjoint_datum","(LieType->RootDatum)");
install_function(SL_wrapper,"SL","(int->RootDatum)");
install_function(GL_wrapper,"GL","(int->RootDatum)");
install_function(simple_roots_wrapper,"simple_roots","(RootDatum->mat)");
install_function(simple_coroots_wrapper,"simple_coroots","(RootDatum->mat)");
install_function(datum_Cartan_wrapper,"Cartan_matrix_of_datum"
		,"(RootDatum->mat)");
install_function(roots_wrapper,"roots","(RootDatum->mat)");
install_function(coroots_wrapper,"coroots","(RootDatum->mat)");
install_function(root_coradical_wrapper,"root_coradical","(RootDatum->mat)");
install_function(coroot_radical_wrapper,"coroot_radical","(RootDatum->mat)");
install_function(dual_datum_wrapper,"dual_datum","(RootDatum->RootDatum)");

@*1 A type for complex reductive groups equipped with an involution.
We shall now go ahead to define a primitive type holding a
|complexredgp::ComplexReductiveGroup|, which represents a complex reductive
group equipped with an involution, which defines an ``inner class'' of real
forms for that group. We can construct such an object from a root datum and an
involution. In the Atlas software such involutions are entered indirectly via
a user interaction, and the way it was constructed is stored and used for
certain purposes. For maximal flexibility, we want to be able to provide an
involution produced in any way we like. This means there is some extra work to
do.

@< Includes... @>=
#include "complexredgp.h"

@*2 Analysing involutions. Our constructor for the current built-in type must
do checking to see that a valid involution is entered, and an analysis of the
involution to replace the values otherwise obtained from the user interaction,
in other words we want to find which sequence of inner class letters could
have been used to construct this involution. For simple factors this is
relatively easy, since one just studies the permutation of the simple roots
that the involution defines (if it does not, it is not an automorphism of the
based root datum). Still, the permutation might not correspond to any inner
class letters: it could correspond to a Complex inner class, but on simple
factors that are not adjacent. Also inner class letters are not unique, since
depending on the type at hand certain inner class letters are synonymous. For
torus factors, the problem requires an analysis of the eigen-lattices of the
involution, which is performed in the class |tori::RealTorus|. But in cases
where the group is not a direct product of its derived group and its central
torus, there might be involutions that can arise for several different
combinations of inner class letters, or not at all. Therefore the value we
compute should be taken as an informed guess.

We consider first the case of a pure torus. For any involution given by an
integral matrix, the space can be decomposed into stable rank~$1$ subspaces
where the involution acts as $+1$ or $-1$ (we call these compact and split
factors, respectively) and stable rank~$2$ subspaces for which the involution
exchanges a pair of basis vectors (we call these Complex factors). The number
of each of these subspaces is uniquely determined by the involution, and
determining them is a small part of what the constructor for a
|tori::RealTorus| does. We shall use that class in the function
|classify_involution|, which returns the number of compact and split factors
(the number of Complex factors follows from there by subtraction from the
rank).

@h "tori.h"
@< Local function def...@>=
std::pair<size_t,size_t> classify_involution
  (const latticetypes::LatticeMatrix& M, size_t r)
throw (std::bad_alloc, std::runtime_error)
{ @< Check that |M| is an $r\times{r}$ matrix defining an involution @>
  tori::RealTorus T(M);
  return std::make_pair(T.compactRank(),T.splitRank());
}

@ The test below that |M| is an involution ($M^2=I$) is certainly necessary
when |classify_involution| is called independently. If the code below is
executed in the context of checking the involution for an inner class, it may
seem redundant given the fact that we shall also check that |M| induces an
involutive permutation of the roots; however even there it is not redundant in
the presence of a torus part.

@< Check that |M| is an $r\times{r}$ matrix defining an involution @>=
{ if (M.numRows()!=r or M.numColumns()!=r) throw std::runtime_error
    ("involution should be a "+num(r)+"x"+num(r)+" matrix, got a "
     +num(M.numRows())+"x"+num(M.numColumns())+" matrix");
  latticetypes::LatticeMatrix I,Q(M);
  identityMatrix(I,r); Q*=M; // $Q=M^2$
  if (!(Q==I)) throw std::runtime_error
      ("given transformation is not an involution");
}

@ The function |classify_involution| might be of use to the user, so let us
make a wrapper for it. In fact we shall return the compact, Complex, and split
ranks, in that order.

@< Local function def...@>=
void classify_wrapper()
{ matrix_ptr M(get<matrix_value>());
  size_t r=M->val.numRows();
  std::pair<size_t,size_t> p=classify_involution(M->val,r);
  push_value(new int_value(p.first)); // compact rank
  push_value(new int_value((r-p.first-p.second)/2)); // Complex rank
  push_value(new int_value(p.second)); // split rank
  wrap_tuple(3);
}

@ We now come to the part of the analysis that involves the root datum. Since
the matrix is already expressed on the basis of the weight lattice used by the
root datum, the question of stabilising that lattice is settled, but we must
check that the matrix is indeed an involution, and that it gives an
automorphism of the based root datum. While we are doing that, we can also
determine the Lie type of the root datum and the inner class letters
corresponding to each of its factors. The Lie type will in general be
identical to that of the root datum (and in particular any torus factors will
come at the end), but in case of Complex inner classes we may be forced to
permute the simple factors to make the identical factors associated to such
classes adjacent. Also note that there is no difference between root data of
type $B_2$ and $C_2$, so we shall classify those as $B_2$ regardless of the
type used to create the root datum.

We do not apply the function |classify_involution| defined above to the entire
matrix describing a (purported) involution of the weight lattice, although we
shall use it below for a matrix defined for the central torus part; we can
however reuse the module that tests for being an involution here.

@h "setutils.h"

@< Local function def...@>=
std::pair<lietype::LieType,lietype::InnerClassType> check_involution
 (const latticetypes::LatticeMatrix& M, const rootdata::RootDatum& rd)
 throw (std::bad_alloc, std::runtime_error)
{ size_t r=rd.rank(),s=rd.semisimpleRank();
  @< Check that |M| is an $r\times{r}$ matrix defining an involution @>
@/setutils::Permutation p(s);
  @< Set |p| to the permutation of the simple roots induced by |M|, or throw
     a |runtime_error| if |M| is not an automorphism of |rd| @>
@/std::pair<lietype::LieType,lietype::InnerClassType> result;
@/lietype::LieType& type=result.first;
  lietype::InnerClassType& inner_class=result.second;
  @< Compute the Lie type |type| and the inner class |inner_class| @>
  return result;
}

@ That |M| is an automorphism means that the simple roots are permuted among
each other, and that the Cartan matrix is invariant under that permutation of
its rows and columns.

@< Set |p| to the permutation of the simple roots...@>=
{ rootdata::WRootIterator first=rd.beginSimpleRoot(), last=rd.endSimpleRoot();
  for (size_t i=0; i<s; ++i)
  { rootdata::Root alpha(r); M.apply(alpha,rd.simpleRoot(i));
    size_t p_i=std::find(first,last,alpha)-first;
    if (p_i<s) p[i]=p_i;
    else throw std::runtime_error
      ("given transformation does not permute simple roots");
  }
  for (size_t i=0; i<s; ++i)
    for (size_t j=0; j<s; ++j)
      if (rd.cartan(p[i],p[j])!=rd.cartan(i,j)) throw std::runtime_error@|
      ("given transformation is not a root datum automorphism");
}

@ For each simple factor we look if there are any non-fixed points of the
permutation (if not we have the compact inner class) and if so, whether the
image of that point lies in the same component of the Dynkin diagram. If the
latter is the case we have an unequal rank inner class, which is actually
called split unless the simple factor is of type $D_{2n}$, and if the image
lies in another component of the diagram we have a Complex inner class.

@h "dynkin.h"

@< Compute the Lie type |type| and the inner class |inner_class| @>=
{ latticetypes::LatticeMatrix(C); rootdata::cartanMatrix(C,rd);
@/dynkin::DynkinDiagram diagr(C);
@/lietype::LieType type0; dynkin::lieType(type0,C); // type before rearranging
@/bitset::RankFlagsList component; // list of connected components of diagram
  dynkin::components(component,diagr);
  std::vector<char> letter(component.size());
  size_t nr_Complex_letters=0;
  for (size_t i=0; i<component.size(); ++i)
  { bool equal_rank=true;
    for (bitset::RankFlags::iterator j=component[i].begin(); j(); ++j)
     // traverse component
      if (p[*j]!=*j) {@; equal_rank=false; break; }
    if (equal_rank) letter[i]='c';
      // identity on this component: compact component
    else if(!component[i][p[component[i].firstBit()]])
      // then component not globally fixed
      {@; ++nr_Complex_letters; letter[i]='C'; } // record Complex component
    else letter[i]= type0[i].first=='D' and type0[i].second%2==0 ? 'u' : 's';
    // unequal rank
  }
  @< Make adaptations for any Complex inner classes, and store the final value
  in |type| and |inner_class| @>
}

@ The complex inner classes pose the additional problem of identifying the
corresponding pairs. We build small tables holding for each Complex inner
class symbol the index of its factor in the Lie type, and the index of the
first root of the corresponding component in the Dynkin diagram, and the
number of the matching factor in the Lie type.

@< Make adaptations for any Complex inner classes... @>=
{ std::vector<size_t> pos(nr_Complex_letters);
  std::vector<size_t> first(nr_Complex_letters);
  for (size_t l=0,i=0; l<component.size(); ++l)
    if (letter[l]=='C')
    {@; pos[i]=l; first[i]=component[l].firstBit(); ++i; }
  std::vector<size_t> buddy(nr_Complex_letters,nr_Complex_letters);
   // initialise with sentinels
  for (size_t i=0; i<nr_Complex_letters; ++i)
    if (buddy[i]==nr_Complex_letters) // then value was not yet set by a buddy
    { size_t b=diagr.component(p[component[pos[i]].firstBit()]).firstBit();
        // first bit of companion
      for (size_t j=i+1; j<nr_Complex_letters; ++j)
        if (first[j]==b) {@; buddy[i]=j; buddy[j]=i; break; }
    }
  type.resize(component.size()+(r-s)); // allocate size for type letters
  for (size_t l=0,k=0,i=0; l<component.size(); ++l)
    if (letter[l]!='C')
    {@; type[k++]=type0[l]; inner_class.push_back(letter[l]); }
    else if (buddy[i]<i) ++i; // skip second member of Complex pair
    else // first member of Complex pair
    { type[k++]=type0[l]; type[k++]=type0[pos[buddy[i]]];
       // add equal types for Complex class
      inner_class.push_back('C'); ++i;
    }
  if (r>s)
  @< Add type letters and inner class symbols for the central torus @>
}

@ Finally we have to give inner class symbols for the central torus. While it
is not entirely clear that forming some quotient could not transform an inner
class specified as |"sc"| by the user into one that we will classify as |"C"|
or vice versa, the symbols we generate have a well defined interpretation:
they correspond to the real torus defined by the involution within the central
torus. The weight lattice of the central torus is the quotient of the full
weight lattice by the rational span of the root lattice, so we determine the
action of the involution on this quotient, and classify it using
|classify_involution|. To find the matrix of the action of the involution on
the quotient lattice, we find a Smith basis for the root lattice (of which the
first $r$ vectors span the lattice to be divided out), express the involution
of that basis (which will have zeros in the bottom-left $(r-s)\times{s}$
block), and extract the bottom-right $(r-s)\times(r-s)$ block.

@< Add type letters and inner class symbols for the central torus @>=
{ using latticetypes::WeightList; using latticetypes::LatticeMatrix;
  for (size_t k=component.size(); k<component.size()+(r-s); ++k)
    type[k]=lietype::SimpleLieType('T',1);
  WeightList b; matrix::initBasis(b,r);
  WeightList simple_roots(rd.beginSimpleRoot(),rd.endSimpleRoot());
@/latticetypes::CoeffList ivf;
  smithnormal::smithNormal(ivf,b.begin(),simple_roots);
@/LatticeMatrix inv(LatticeMatrix(M,b),s,s,r,r);
    // involution on quotient by root lattice
  std::pair<size_t,size_t> cl=classify_involution(inv,r-s);
@/size_t& compact_rank=cl.first;
  size_t& split_rank=cl.second;
  size_t Complex_rank=(r-s-compact_rank-split_rank)/2;
@/while (compact_rank-->0) inner_class.push_back('c');
  while (Complex_rank-->0) inner_class.push_back('C');
  while (split_rank-->0) inner_class.push_back('s');
}

@*2 Storing the inner class values.
Although abstractly an inner class value is described completely by an object
of type |complexredgp::ComplexReductiveGroup|, we shall need to record
additional information in order to be able to present meaningful names for the
real forms and dual real forms in this inner class. The above analysis of
involutions was necessary in order to obtain such information; it will be
recorded in values of type |realform_io::Interface|.
@< Includes... @>=
#include "realform_io.h"

@~The class |inner_class_value| will be the first built-in type where we
deviate from the previously used scheme of holding an \.{atlas} object with
the main value in a data member |val|. The reason is that the copy constructor
for |complexredgp::ComplexReductiveGroup| is private (and nowhere defined), so
that the straightforward definition of a copy constructor for such a built-in
type would not work, and the copy constructor is necessary for the |clone|
method. So instead, we shall share the \.{atlas} object when duplicating our
value, and maintain a reference count to allow destruction when the last copy
disappears. The reference count needs to be shared of course, and since the
links between the |inner_class_value| and both the \.{atlas} value and the reference count are indissoluble, we use references for the
members |val| and |ref_count|.

The reference to |val| is not |const|, since at some point we need to apply
the method |fillCartan| to it, which considered to be a manipulator and so is
not defined as a |const| method (even though technically it could have been,
given the fact that it only modifies values at the other side of the
|d_cartan| pointer). The variability of |val| may in fact be noticeable to the
user through effects on the numbering of Cartan classes.

The main constructor takes a auto-pointer to a |ComplexReductiveGroup| as
argument, as a reminder that the caller gives up ownership of this pointer
that should come from a call to~|new|; this pointer will henceforth be owned
by the |inner_class_value| constructed, in shared ownership with any values
later cloned from it: the last one of them to be destroyed will call |delete|
for the pointer. The remaining pair of arguments of the main constructor must
have been computed by |check_involution| above, to ensure their validity (we
prefer not to call that function from inside our constructor, which would have
guaranteed this).

Occasionally we shall need to refer to the dual inner class (for the dual
group); since the construction of an instance even without its complete set of
Cartan classes takes some work, we do not wish to repeat that every time it is
needed, so we create the dual |complexredgp::ComplexReductiveGroup| value upon
construction of the |inner_class_value| and store it in the |dual| field where
it will be available as needed.

Contrary to what was the case for other value types, the copy constructor is
public here. This makes it possible to store an |inner_class_value| inside
another class, when that class depends on the object referenced by our |val|
being alive; doing so does not cost much, and the reference counting mechanism
will then ensure that the inner class remains valid at least as long as the
object of that containing class exists.

@< Type definitions @>=
struct inner_class_value : public value_base
{ complexredgp::ComplexReductiveGroup& val;
  complexredgp::ComplexReductiveGroup& dual;
  size_t& ref_count;
  const realform_io::Interface interface,dual_interface;
@)
  inner_class_value(std::auto_ptr<complexredgp::ComplexReductiveGroup>
   ,@| lietype::LieType,lietype::InnerClassType); // main constructor
  ~inner_class_value();
@)
  virtual void print(std::ostream& out) const;
  inner_class_value* clone() const @+
    {@; return new inner_class_value(*this); }
  static const char* name() @+{@; return "inner class"; }
@)
  inner_class_value(const inner_class_value& v); // copy constructor
  inner_class_value(const inner_class_value& v,tags::DualTag);
   // constructor of dual
};
@)
typedef std::auto_ptr<inner_class_value> inner_class_ptr;

@ Here are the copy constructor and the destructor.
@< Function def...@>=
inner_class_value::inner_class_value(const inner_class_value& v)
: val(v.val), dual(v.dual), ref_count(v.ref_count)
, interface(v.interface), dual_interface(v.dual_interface)
{@; ++ref_count; }

inner_class_value::~inner_class_value()
{@; if (--ref_count==0) {@; delete &val; delete &dual; delete &ref_count;} }

@ The main constructor installs the reference to the Atlas value, creates its
dual value to which a reference is stored, and allocates and initialises the
reference count. To initialise the |interface| and |dual_interface| fields, we
call the appropriate constructors with a |layout::Layout| structure provided
by a constructor that we added to it specifically for the purpose, from values
|lt|, |ict| passed to the main constructor. The field of |layout::Layout|
holding a lattice basis will remain empty, but this field is unused by the
constructors for |realform_io::Interface|.

The current constructor has the defect that the reference to which the |dual|
field is initialised will not be freed if an exception should be thrown during
the subsequent initialisations (which is not likely, but not impossible
either). However, we see no means to correct this defect with the given
structure, since a member of reference type has to be initialised to its
definitive value, and no local variables are available during initialisation
whose destructor could take care of liberating the memory referred to
by~|dual|.

@< Function def...@>=
inner_class_value::inner_class_value
  (std::auto_ptr<complexredgp::ComplexReductiveGroup> g
  ,lietype::LieType lt, lietype::InnerClassType ict)
@/: val(*g)
, dual(*new complexredgp::ComplexReductiveGroup(*g,tags::DualTag()))
@/, ref_count(*new size_t(1))
@/, interface(*g,layout::Layout(lt,ict))
, dual_interface(*g,layout::Layout(lt,ict),tags::DualTag())
 {@; g.release(); } // now that we own |g|, release the auto-pointer

@ We allow construction of a dual |inner_class_value|. Since it can share the
two fields referring to objects of type |complexredgp::ComplexReductiveGroup|
in the opposite order, we can consider it as a member of the same reference
counted family, and share the |ref_count| field. This means this constructor
is more like the copy constructor than like the main constructor, and in
particular the reference count is increased. The dual of the dual inner class
value will behave exactly like a copy of the original inner class.

@< Function def...@>=
inner_class_value::inner_class_value(const inner_class_value& v,tags::DualTag)
@/: val(v.dual), dual(v.val)
, ref_count(v.ref_count)
@/, interface(v.dual_interface), dual_interface(v.interface)
{@; ++ref_count; }


@ One of the most practical informations about a |ComplexReductiveGroup|,
which is available directly after its construction, is the number of real
forms in the inner class defined by it; we print this information when a
|inner_class_value| is printed.

@< Function def...@>=
void inner_class_value::print(std::ostream& out) const
{ out << "Complex reductive group equipped with an involution,\n" @|
         "defining an inner class of "
      << val.numRealForms() @| << " real "
      << (val.numRealForms()==1 ? "form" : "forms") @| << " and "
      << val.numDualRealForms() @| << " dual real "
      << (val.numDualRealForms()==1 ? "form" : "forms");
}

@ So here is our wrapper function for building a complex reductive group with
involution. Apart from the test performed, we have to deal with the question
of giving the complex group ownership of the root datum. Currently the
constructor of the |ComplexReductiveGroup| does not assume ownership
of the root datum until it completes, but this should be corrected in the
constructor, not here; therefore we shall pass a pointer to a copy of the root
datum.

@< Local function def...@>=
void fix_involution_wrapper()
{ push_tuple_components();
  matrix_ptr M(get<matrix_value>());
  root_datum_ptr rd(get<root_datum_value>());
  std::pair<lietype::LieType,lietype::InnerClassType> cl
    =check_involution(M->val,rd->val);
  if (verbosity>0) @< Report the type and inner class found @>
  rootdata::RootDatum* rdp=new rootdata::RootDatum(rd->val);
   // just a shorthand
  std::auto_ptr<complexredgp::ComplexReductiveGroup>@|
    G(new complexredgp::ComplexReductiveGroup(rdp,M->val));
  push_value(new inner_class_value(G,cl.first,cl.second));
}

@ For understanding what is going on, the user may find it useful to look at
the Lie type and inner class that were determined from the root datum and the
involution given.

@< Report the type and inner class found @>=
{ Lie_type_value t(cl.first);
  *output_stream << "Found " << t << ", and inner class '";
  for (size_t i=0; i<cl.second.size(); ++i)
    *output_stream << cl.second[i];
  *output_stream << "'.\n";
}

@ To simulate the functioning of the Atlas software, the function $set\_type$
takes as argument the name of a Lie type, a matrix giving kernel generators,
and a string describing the inner class. The evaluation of the call
$set\_type(lt,gen,ic)$ effectively consists of setting $t={\it
Lie\_type(lt)}$, ${\it basis}={\it quotient\_basis(t,gen)}$, and then
returning the value of ${\it fix\_involution (root\_datum
(t,basis),based\_involution(t,basis,ic))}$. Note that neither the Lie type nor
the inner class type are transmitted directly to the |inner_class_type|
constructor; we are depending on |fix_involution_wrapper()| on reconstructing
them. This can give some surprises such as transforming type $C_2$ into $B_2$,
or inner class letters changing to synonyms (but the latter should not have
externally visible effects).

@< Local function def...@>=
void set_type_wrapper()
{ push_tuple_components();
  string_ptr ict(get<string_value>());
  matrix_ptr gen(get<matrix_value>());
  Lie_type_wrapper(); // convert string to Lie type
  Lie_type_ptr t(get<Lie_type_value>());
@)
  push_value(t->clone()); push_value(gen);
  wrap_tuple(2); quotient_basis_wrapper();
  matrix_ptr basis(get<matrix_value>());
@)
  push_value(t->clone()); push_value(basis->clone());
  wrap_tuple(2); root_datum_wrapper();
@/push_value(t); push_value(basis);
  push_value(ict);
  wrap_tuple(3); based_involution_wrapper();
@/wrap_tuple(2); fix_involution_wrapper();
}

@ It can be awkward to use $set\_type$ which wants to construct the root datum
in its own way, for instance if one already has a root datum as produced by
special functions like $GL$. Unfortunately the precise sub-lattice used to
construct a root datum cannot be recovered from the root datum in all cases,
nor indeed the Lie type, and these are needed to find the proper involution,
as in |based_involution_wrapper|. The problems are only caused by possible
torus factors however, and in most cases one can reconstruct the values with
sufficient accuracy. In fact the function below will always deduce a Lie type
and a lattice basis that \emph{could have} been used to obtain the root datum.

The function call ${\it set\_inner\_class(rd,ict)}$ will deduce the type as
$t={\it type\_of\_root\_datum(rd)}$ and the basis of the sub-lattice as ${\it
basis}= {\it transpose\_mat(coroot\_radical(rd))}$, and then return the value
of the expression ${\it fix\_involution(rd,based\_involution(t,basis,ict))}$.

@< Local function def...@>=
void set_inner_class_wrapper()
{ push_tuple_components();
  string_ptr ict(get<string_value>());
  root_datum_ptr rd(get<root_datum_value>());
@)
  push_value(rd->clone()); type_of_root_datum_wrapper();
  Lie_type_ptr t(get<Lie_type_value>());
@)
  push_value(rd->clone());
  coroot_radical_wrapper(); transpose_mat_wrapper();
  matrix_ptr basis(get<matrix_value>());
@)
  push_value(rd);
@/push_value(t); push_value(basis);
  push_value(ict);
  wrap_tuple(3); based_involution_wrapper();
@/wrap_tuple(2); fix_involution_wrapper();
}

@*2 Functions operating on complex reductive groups.
%
Here are our first functions that access a |inner_class_value|; they
recover the ingredients that were used in the construction, and construct the
dual inner class.

@< Local function def...@>=
void distinguished_involution_wrapper()
{ inner_class_ptr G(get<inner_class_value>());
  push_value(new matrix_value(G->val.distinguished()));
}

void root_datum_of_inner_class_wrapper()
{ inner_class_ptr G(get<inner_class_value>());
  push_value(new root_datum_value(G->val.rootDatum()));
}

void dual_inner_class_wrapper()
{ inner_class_ptr G(get<inner_class_value>());
  push_value(new inner_class_value(*G,tags::DualTag()));
}

@ More interestingly, let us extract the list of names of the real forms.
This uses the interface fields stored in the value. Since they exist for both
the group itself and for the dual group, we define an auxiliary function that
produces the list, and then use it twice.

@< Local function def...@>=
void push_name_list(const realform_io::Interface& interface)
{ row_ptr result
    (new row_value(std::vector<value>()));
  for (size_t i=0; i<interface.numRealForms(); ++i)
    result->val.push_back(new string_value(interface.typeName(i)));
  push_value(result);
}

void form_names_wrapper()
{@; inner_class_ptr G(get<inner_class_value>());
  push_name_list(G->interface);
}

void dual_form_names_wrapper()
{@; inner_class_ptr G(get<inner_class_value>());
  push_name_list(G->dual_interface);
}

@ And now, our first function that really simulates something that can be done
using the Atlas interface, and that uses more than a root datum. This is the
\.{blocksizes} command from \.{mainmode.cpp}, which uses
|innerclass_io::printBlockSizes|, but we have to rewrite it to avoid calling
an output routine. A subtle difference is that we use a matrix of integers
rather than of unsigned long integers to collect the block sizes; this avoids
having to define a new primitive type, and seems to suffice for the cases that
are currently computationally feasible.

@< Local function def...@>=
void block_sizes_wrapper()
{ inner_class_ptr G(get<inner_class_value>());
  G->val.fillCartan();
  matrix_ptr
  M(new matrix_value @|
    (latticetypes::LatticeMatrix(G->val.numRealForms()
                                ,G->val.numDualRealForms())
    ));
  for (size_t i = 0; i < M->val.numRows(); ++i)
    for (size_t j = 0; j < M->val.numColumns(); ++j)
      M->val(i,j) =
      G->val.blockSize(G->interface.in(i),G->dual_interface.in(j));
  push_value(M);
}

@ Next we shall provide a function that displays the occurrence of Cartan
classes for various real forms. The rows will be indexed by real forms, and
the columns by Cartan classes (note the alliteration).

@< Local function def...@>=
void occurrence_matrix_wrapper()
{ inner_class_ptr G(get<inner_class_value>());
  G->val.fillCartan();
  size_t nr=G->val.numRealForms();
  size_t nc=G->val.numCartanClasses();
  matrix_ptr
    M(new matrix_value(latticetypes::LatticeMatrix(nr,nc)));
  for (size_t i=0; i<nr; ++i)
  { bitmap::BitMap b=G->val.cartanSet(G->interface.in(i));
    for (size_t j=0; j<nc; ++j)
      M->val(i,j)= b.isMember(j) ? 1 : 0;
  }
  push_value(M);
}

@ We do the same for dual real forms. Note that we had to introduce the method
|dualCartanSet| for |complexredgp::ComplexReductiveGroup| in order to be able
to write this function.

@< Local function def...@>=
void dual_occurrence_matrix_wrapper()
{ inner_class_ptr G(get<inner_class_value>());
  G->val.fillCartan();
  size_t nr=G->val.numDualRealForms();
  size_t nc=G->val.numCartanClasses();
  matrix_ptr
    M(new matrix_value(latticetypes::LatticeMatrix(nr,nc)));
  for (size_t i=0; i<nr; ++i)
  { bitmap::BitMap b=G->val.dualCartanSet(G->dual_interface.in(i));
    for (size_t j=0; j<nc; ++j)
      M->val(i,j)= b.isMember(j) ? 1 : 0;
  }
  push_value(M);
}

@ Finally we install everything.
@< Install wrapper functions @>=
install_function(classify_wrapper,"classify_involution"
                ,"(mat->int,int,int)");
install_function(fix_involution_wrapper,"fix_involution"
                ,"(RootDatum,mat->InnerClass)");
install_function(set_type_wrapper,"set_type"
                ,"(string,mat,string->InnerClass)");
install_function(set_inner_class_wrapper,"set_inner_class"
                ,"(RootDatum,string->InnerClass)");
install_function(distinguished_involution_wrapper,"distinguished_involution"
                ,"(InnerClass->mat)");
install_function(root_datum_of_inner_class_wrapper,"root_datum_of_inner_class"
                ,"(InnerClass->RootDatum)");
install_function(dual_inner_class_wrapper,"dual_inner_class"
                ,"(InnerClass->InnerClass)");
install_function(form_names_wrapper,"form_names"
                ,"(InnerClass->[string])");
install_function(dual_form_names_wrapper,"dual_form_names"
                ,"(InnerClass->[string])");
install_function(block_sizes_wrapper,"block_sizes"
                ,"(InnerClass->mat)");
install_function(occurrence_matrix_wrapper,"occurrence_matrix"
                ,"(InnerClass->mat)");
install_function(dual_occurrence_matrix_wrapper,"dual_occurrence_matrix"
                ,"(InnerClass->mat)");

@*1 A type for real reductive groups.
A next step in specifying the computational context is choosing a real form in
the inner class of them that was determined by a root datum and an involution.
This determines a, not necessarily connected, real reductive group inside the
connected complex reductive group; the corresponding Atlas class is called
|realredgp::RealReductiveGroup|.

@< Includes... @>=
#include "realredgp.h"

@*2 Class definition.
The layout of this type of value is different from what we have seen before.
An Atlas object of class |realredgp::RealReductiveGroup| is dependent upon
another Atlas object to which it stores a pointer, which is of type
|complexredgp::ComplexReductiveGroup|, so we must make sure that the object
pointed to cannot disappear before it does. The easiest way to do this is to
place a |inner_class_value| object |parent| inside the |real_form_value| class
that we shall now define; the reference-counting scheme introduced above then
guarantees that the data we depend upon will remain in existence sufficiently
long. Since that data can be accessed from inside the
|realredgp::RealReductiveGroup|, we shall mostly mention the |parent| to
access its |interface| and |dual_interface| fields. To remind us that the
|parent| is not there to be changed by us, we declare it |const|. The object
referred to may in fact undergo internal change however, via manipulators of
the |val| field.

We shall later derive from this structure, without adding data members, a
structure to record dual real forms, which are constructed in the context of
the same inner class, but where the |val| field is constructed for the dual
group (and its dual inner class). In order to accommodate for that we provide
an additional protected constructor that will be defined later; this also
explains why the copy constructor is protected rather than private.

@< Type definitions @>=
struct real_form_value : public value_base
{ const inner_class_value parent;
  realredgp::RealReductiveGroup val;
@)
  real_form_value(inner_class_value p,realform::RealForm f)
  : parent(p), val(p.val,f) @+{}
  ~real_form_value() @+{} // everything is handled by destructor of |parent|
@)
  virtual void print(std::ostream& out) const;
  real_form_value* clone() const @+
    {@; return new real_form_value(*this); }
  static const char* name() @+{@; return "real form"; }
protected:
  real_form_value(inner_class_value,realform::RealForm,tags::DualTag);
     // for dual forms
  real_form_value(const real_form_value& v)
  : parent(v.parent), val(v.val) @+{} // copy c'tor
};
@)
typedef std::auto_ptr<real_form_value> real_form_ptr;

@ When printing a real form, we give the name by which it was chosen (where
for once we do use the |parent| field), and provide some information about its
connectedness. Since the names of the real forms are indexed by their outer
number, but the real form itself stores its inner number, we must somewhat
laboriously make the conversion here.

@< Function def...@>=
void real_form_value::print(std::ostream& out) const
{ out << (val.isQuasisplit() ? "quasisplit " : "")
  << "real form '" @|
  << parent.interface.typeName(parent.interface.out(val.realForm())) @|
  << "', defining a"
  << (val.isConnected() ? " connected" : "" ) @|
  << (val.isSplit() ? " split" : "") @|
  << " real group" ;
}

@ To make a real form is easy, one provides an |inner_class_value| and a valid
index into its list of real forms. Since this number coming from the outside
is to be interpreted as an outer index, we must convert it to an inner index
at this point. As a special case we also provide the quasisplit form.

@< Local function def...@>=
void real_form_wrapper()
{ push_tuple_components();
  int_ptr i(get<int_value>());
  inner_class_ptr G(get<inner_class_value>());
  if (size_t(i->val)>=G->val.numRealForms())
    throw std::runtime_error ("illegal real form number: "+num(i->val));
  push_value(new real_form_value(*G,G->interface.in(i->val)));
}
@)
void quasisplit_form_wrapper()
{ inner_class_ptr G(get<inner_class_value>());
  push_value(new real_form_value(*G,G->val.quasisplit()));
}

@*2 Functions operating on real reductive groups.
%
Here is a function that gives information about the dual component group
of a real reductive group: it returns the rank~$r$ such that the component
group is isomorphic to $(\Z/2\Z)^r$.

@< Local function def...@>=
void components_rank_wrapper()
{ real_form_ptr R(get<real_form_value>());
  const latticetypes::ComponentList c=R->val.dualComponentReps();
  push_value(new int_value(c.size()));
}

@ And here is one that counts the number of Cartan classes for the real form.
Note that we must call |fillCartan()| here to ensure that the Cartan classes
are in fact generated.

@< Local function def...@>=
void count_Cartans_wrapper()
{ real_form_ptr rf(get<real_form_value>());
  rf->val.fillCartan();
  push_value(new int_value(rf->val.numCartan()));
}

@ The size of the finite set $K\backslash G/B$ can be determined from the real
form, once the corresponding Cartan classes have been generated.

@< Local function def...@>=
void KGB_size_wrapper()
{ real_form_ptr rf(get<real_form_value>());
  rf->val.fillCartan();
  push_value(new int_value(rf->val.kgbSize()));
}

@ Once the Cartan classes for a real form are constructed, we have a partial
ordering on them. The function |Cartan_order_matrix| displays a matrix for
this partial ordering. This function more or less replaces the \.{corder}
command in atlas.

@h "poset.h"
@< Local function def...@>=
void Cartan_order_matrix_wrapper()
{ real_form_ptr rf(get<real_form_value>());
  rf->val.fillCartan();
  size_t n=rf->val.numCartan();
  matrix_value* M =
     new matrix_value(latticetypes::LatticeMatrix(n,n,0));
  const poset::Poset& p = rf->val.cartanOrdering();
  for (size_t i=0; i<n; ++i)
    for (size_t j=i; j<n; ++j)
      if (p.lesseq(i,j)) M->val(i,j)=1;

  push_value(M);
}

@*2 Dual real forms.
Although they could be considered as real forms for the dual group and dual
inner class, it will be convenient to be able to construct dual real forms in
the context of the |inner_class_value| itself. The resulting objects will have
a different type than ordinary real forms, but the will use the same structure
layout. The interpretation of the fields is as follows for dual real forms:
the |parent| field contains a copy of the original |inner_class_value|, but
the |val| field contains a |realredgp::RealReductiveGroup| object
constructed for the dual inner class, so that its characteristics will be
correct.

@< Function def...@>=

real_form_value::real_form_value(inner_class_value p,realform::RealForm f
				,tags::DualTag)
: parent(p), val(p.dual,f)
{}

@ In order to be able to use the same layout for dual real forms but to
nevertheless make a distinction (for instance in printing), we must derive a
new type from |real_form_value|.

@< Type definitions @>=
struct dual_real_form_value : public real_form_value
{ dual_real_form_value(inner_class_value p,realform::RealForm f)
  : real_form_value(p,f,tags::DualTag()) @+{}
@)
  virtual void print(std::ostream& out) const;
  dual_real_form_value* clone() const @+
    {@; return new dual_real_form_value(*this); }
  static const char* name() @+{@; return "dual real form"; }
private:
  dual_real_form_value(const dual_real_form_value& v)
  : real_form_value(v) @+ {}
};
@)
typedef std::auto_ptr<dual_real_form_value> dual_real_form_ptr;

@ The only real difference with real forms is a slightly modified output
routine.

@< Function def...@>=
void dual_real_form_value::print(std::ostream& out) const
{ out << (val.isQuasisplit() ? "quasisplit " : "")
  << "dual real form '" @|
  << parent.dual_interface.typeName
    (parent.dual_interface.out(val.realForm())) @|
  << "'";
}

@ To make a dual real form, one provides an |inner_class_value| and a valid
index into its list of dual real forms, which will be converted to an inner
index. We also provide the dual quasisplit form.

@< Local function def...@>=
void dual_real_form_wrapper()
{ push_tuple_components();
  int_ptr i(get<int_value>());
  inner_class_ptr G(get<inner_class_value>());
  if (size_t(i->val)>=G->val.numDualRealForms())
    throw std::runtime_error ("illegal dual real form number: "+num(i->val));
  push_value(new dual_real_form_value(*G,G->dual_interface.in(i->val)));
}
@)
void dual_quasisplit_form_wrapper()
{ inner_class_ptr G(get<inner_class_value>());
  push_value(new dual_real_form_value(*G,G->dual.quasisplit()));
}

@ Rather than provide all functions for real forms for dual real forms, we
provide a function that converts it to a real form, which of course will be
associated to the dual inner class.

@< Local function def...@>=
void real_form_from_dual_wrapper()
{ dual_real_form_ptr d(get<dual_real_form_value>());
  push_value(new real_form_value
                 (inner_class_value(d->parent,tags::DualTag())
                 ,d->val.realForm()));
}

@*1 A type for Cartan classes.
Another type of value associated to inner classes are Cartan classes, which
describe equivalence classes (for ``stable conjugacy'') of Cartan subgroups of
real reductive groups in an inner class. The Atlas software associates a fixed
set of Cartan classes to each inner class, and records for each real form the
subset of those Cartan classes that occur for the real form. Contrary to real
forms, the complete list of Cartan classes is not generated directly upon
construction of the |ComplexReductiveGroup|; to ensure that all Cartan classes
for a given real form~|f| of an inner class~|C| exist, one should first call
|C.fillCartan(f)| (or equivalently the method |fillCartan()| for a
|RealReductiveGroup| object constructed using~|f|), or call |C.fillCartan()|
which will generate the full list of Cartan classes for this inner class.

@< Includes... @>=
#include "cartanclass.h"

@*2 Class definition.
The Cartan class is identified by a number within the inner class. This number
is actually its sequence number in the order in which Cartan classes are
generated in this inner class, so it may vary from one run to another;
nevertheless the number is important to record, since some information about
the Cartan class is only stored at the level of the inner class, and the
number is required to obtain this information.

@< Type definitions @>=
struct Cartan_class_value : public value_base
{ const inner_class_value parent;
  size_t number;
  const cartanclass::CartanClass& val;
@)
  Cartan_class_value(const inner_class_value& p,size_t cn);
  ~Cartan_class_value() @+{} // everything is handled by destructor of |parent|
@)
  virtual void print(std::ostream& out) const;
  Cartan_class_value* clone() const @+
    {@; return new Cartan_class_value(*this); }
  static const char* name() @+{@; return "Cartan class"; }
private:
  Cartan_class_value(const Cartan_class_value& v)
  : parent(v.parent), number(v.number), val(v.val) @+{}
};
@)
typedef std::auto_ptr<Cartan_class_value> Cartan_class_ptr;

@ In the constructor we check that the Cartan class with the given number
currently exists. This check \emph{follows} the selection of the pointer to
the |CartanClass| object that initialises the |val| field because we have no
choice (a reference field cannot be set later); the danger in doing this is
limited since the selection itself will probably not crash the program, and
the nonsensical pointer so obtained will probably be stored away without
inspection and then be forgotten when an error is thrown. Moreover, the
function we shall provide to access Cartan classes will ensure that when
constructed, the Cartan class does actually exist.

@< Function def...@>=
Cartan_class_value::Cartan_class_value(const inner_class_value& p,size_t cn)
: parent(p),number(cn),val(p.val.cartan(cn))
{ if (cn>=p.val.numCartanClasses()) throw std::runtime_error
  (std::string("Cartan class number ")+num(cn)+" does not (currently) exist");
}

@ When printing a Cartan class, we show its number, and for how many real
forms and dual real forms it is valid.

@< Function def...@>=
void Cartan_class_value::print(std::ostream& out) const
{ out << "Cartan class #" << number << ", occurring for " @|
  << val.numRealForms() << " real "
  << (val.numRealForms()==1 ? "form" : "forms") << " and for "@|
  << val.numDualRealForms() << " dual real "
  << (val.numDualRealForms()==1 ? "form" : "forms");
}

@ To make a Cartan class, one must provide a |real_form_value| together with a
valid index into its list of Cartan classes.

@< Local function def...@>=
void Cartan_class_wrapper()
{ push_tuple_components();
  int_ptr i(get<int_value>());
  real_form_ptr rf(get<real_form_value>());
  rf->val.fillCartan();
  if (size_t(i->val)>=rf->val.numCartan())
    throw std::runtime_error
    ("illegal Cartan class number: "+num(i->val)
    +", this real form only has "+num(rf->val.numCartan())+" of them");
  bitmap::BitMap cs=rf->val.cartanSet();
  push_value(new Cartan_class_value(rf->parent,cs.n_th(i->val)));
}

@ Like the quasisplit real form for inner classes, there is a particular
Cartan class for each real form, the ``most split Cartan''; it is of
particular importance for knowing which dual real forms can be associated to
this real form. It is (probably) the last Cartan class in the list for the
real form, but we have a direct access to it via the |mostSplit| method for
|realredgp::RealReductiveGroup|.

@< Local function def...@>=
void most_split_Cartan_wrapper()
{ real_form_ptr rf(get<real_form_value>());
  rf->val.fillCartan();
  push_value(new Cartan_class_value(rf->parent,rf->val.mostSplit()));
}


@*2 Functions operating on Cartan classes.
%
This function and the following provide the functionality of the Atlas
command \.{cartan}. They are based on |complexredgp_io::printCartanClass|, but
rewritten to take into account the fact that we have no |Interface| objects
for complex groups. We omit in our first function {\it Cartan\_info} the final
call to |cartan_io::printFiber| that would list all the real forms for which
this Cartan class exists with the corresponding part of the adjoint fiber
group, relegating it instead to a second function {\it fiber\_part} that
operates on a per-real-form basis. This separation seems more natural in a
setup where real forms and Cartan classes are presented as more or less
independent objects, and it also avoids having to include parts of the Atlas
software that we don't really need, like input handling, just because the
compilation unit \.{cartan\_io} uses them.

@h "prettyprint.h"

@< Local function def...@>=
void print_Cartan_info_wrapper()
{ using atlas::operator<<;
  Cartan_class_ptr cc(get<Cartan_class_value>());

  const rootdata::RootDatum& rd=cc->parent.val.rootDatum();

  prettyprint::printTorusType(*output_stream,cc->val.fiber().torus())
  << std::endl;

  *output_stream << "twisted involution orbit size: " << cc->val.orbitSize()
   << std::endl;

@)// print type of imaginary root system
  lietype::LieType ilt,rlt,clt;
  rootdata::lieType(ilt,cc->val.simpleImaginary(),rd);

  if (ilt.size() == 0)
    *output_stream << "imaginary root system is empty" << std::endl;
  else
    *output_stream << "imaginary root system: " << ilt << std::endl;

@)// print type of real root system
  rootdata::lieType(rlt,cc->val.simpleReal(),rd);

  if (rlt.size() == 0)
    *output_stream << "real root system is empty" << std::endl;
  else
    *output_stream << "real root system: " << rlt << std::endl;

@)// print type of complex root system
  rootdata::lieType(clt,cc->val.simpleComplex(),rd);

  if (clt.size() == 0)
    *output_stream << "complex factor is empty" << std::endl;
  else
    *output_stream << "complex factor: " << clt << std::endl;
@)
  wrap_tuple(0);
}

@ A functionality that is implicit in the Atlas command \.{cartan} is the
enumeration of all real forms corresponding to a given Cartan class. While
|cartan_io::printFiber| traverses the real forms in the order corresponding to
the parts of the partition |f.weakReal()| for the fiber~|f| associated to the
Cartan class, it is not necessary to use this order, and we can instead simply
traverse all real forms and check whether the given Cartan class exists for
them. Taking care to convert real form numbers to their inner representation
here, we can in fact return a list of real forms; we also include a version
for dual real forms.

@< Local function def...@>=
void real_forms_of_Cartan_wrapper()
{ Cartan_class_ptr cc(get<Cartan_class_value>());
  const inner_class_value& ic=cc->parent;
@/row_ptr result @|
    (new row_value(std::vector<value>(cc->val.numRealForms())));
  for (size_t i=0,k=0; i<ic.val.numRealForms(); ++i)
  { bitmap::BitMap b(ic.val.cartanSet(ic.interface.in(i)));
    if (b.isMember(cc->number))
      result->val[k++]=new real_form_value(ic,ic.interface.in(i));
  }
  push_value(result);
}
@)
void dual_real_forms_of_Cartan_wrapper()
{ Cartan_class_ptr cc(get<Cartan_class_value>());
  const inner_class_value& ic=cc->parent;
@/row_ptr result @|
    (new row_value(std::vector<value>(cc->val.numDualRealForms())));
  for (size_t i=0,k=0; i<ic.val.numDualRealForms(); ++i)
  { bitmap::BitMap b(ic.val.dualCartanSet(ic.dual_interface.in(i)));
    if (b.isMember(cc->number))
      result->val[k++]=new dual_real_form_value(ic,ic.dual_interface.in(i));
  }
  push_value(result);
}

@ For the fiber group partition information that was not produced by
|print_Cartan_info|, we use a Cartan class and a real form as parameters. This
function returns the part of the |weakReal| partition stored in the fiber of
the Cartan class that corresponds to the given (weak) real form. The numbering
of the parts of that partition are not the numbering of the real forms
themselves, so they must be translated through the |realFormLabels| list for
the Cartan class, which must be obtained from its |parent| inner class. If the
Cartan class does not exist for the given real form, then it will not occur in
that |realFormLabels| list, and the part returned here will be empty. The part
of the partition is returned as a list of integral values.

@f pi NULL

@< Local function def...@>=
void fiber_part_wrapper()
{ push_tuple_components();
  real_form_ptr rf(get<real_form_value>());
  Cartan_class_ptr cc(get<Cartan_class_value>());
  if (&rf->parent.val!=&cc->parent.val)
    throw std::runtime_error
    ("inner class mismatch between real form and Cartan class");
  bitmap::BitMap b(cc->parent.val.cartanSet(rf->val.realForm()));
  if (!b.isMember(cc->number))
    throw std::runtime_error
    ("fiber_part: Cartan class not defined for this real form");
@)
  const partition::Partition& pi = cc->val.fiber().weakReal();
  const realform::RealFormList rf_nr=
     cc->parent.val.realFormLabels(cc->number);
     // translate part number of |pi| to real form
  row_ptr result
    (new row_value(std::vector<value>()));
  for (size_t i=0; i<pi.size(); ++i)
    if ( rf_nr[pi(i)] == rf->val.realForm())
      result->val.push_back(new int_value(i));
  push_value(result);
}

@ The function |print_gradings| gives on a per-real-form basis the
functionality of the Atlas command \.{gradings} that is implemented by
|complexredgp_io::printGradings| and |cartan_io::printGradings|. It therefore
takes, like |fiber_part|, a Cartan class and a real form as parameter. Its
output consist of a list of $\Z/2\Z$-gradings of each of the fiber group
elements in the part corresponding to the real form, where each grading is a
sequence of bits corresponding to the simple imaginary roots.

@f sigma NULL

@< Local function def...@>=
void print_gradings_wrapper()
{ push_tuple_components();
  real_form_ptr rf(get<real_form_value>());
@/Cartan_class_ptr cc(get<Cartan_class_value>());
  if (&rf->parent.val!=&cc->parent.val)
    throw std::runtime_error
    ("inner class mismatch between real form and Cartan class");
  bitmap::BitMap b(cc->parent.val.cartanSet(rf->val.realForm()));
  if (!b.isMember(cc->number))
    throw std::runtime_error
    ("fiber_part: Cartan class not defined for this real form");
@)
  const partition::Partition& pi = cc->val.fiber().weakReal();
  const realform::RealFormList rf_nr=
     cc->parent.val.realFormLabels(cc->number);
     // translate part number of |pi| to real form
@)
  const rootdata::RootList& si = cc->val.fiber().simpleImaginary();
  // simple imaginary roots

  latticetypes::LatticeMatrix cm;
  setutils::Permutation sigma;
@/@< Compute the Cartan matrix |cm| of the root subsystem |si|, and the
     permutation |sigma| giving the Bourbaki numbering of its simple roots @>

  @< Print information about the imaginary root system and its simple roots @>

  @< Print the gradings for the part of |pi| corresponding to the real form @>

  wrap_tuple(0);
}

@ For normalising the ordering of the simple imaginary roots we use two
functions from \.{dynkin.cpp}.

@< Compute the Cartan matrix |cm|... @>=
{ rootdata::cartanMatrix(cm,si,cc->parent.val.rootDatum());
  dynkin::DynkinDiagram d(cm); dynkin::bourbaki(sigma,d);
}

@ The imaginary root system might well be empty, so we make special provisions
for this case.

@< Print information about the imaginary root system... @>=
{ *output_stream << "Imaginary root system is ";
  if (si.size()==0) *output_stream<<"empty.\n";
  else
  { using atlas::operator<<;
    lietype::LieType t; dynkin::lieType(t,cm);

    *output_stream << "of type " << t << ", with simple root"
              << (si.size()==1 ? " " : "s ");
    for (size_t i=0; i<si.size(); ++i)
      *output_stream << si[sigma[i]] << (i<si.size()-1 ? "," : ".\n");
  }
}


@ The gradings are computed by the |grading| method of the fiber of our Cartan
class, and converted to a string of characters |'0'| and |'1'| by
|prettyprint::prettyPrint|. The permutation |sigma| was computed in order to
present these bits in an order adapted to the Dynkin diagram of the imaginary
root system. If that system is empty, the strings giving the gradings are
empty, but the part of the partition |pi| will not be, unless the Cartan class
does not exist for the given real form. Rather than printing directly to
|*output_stream|, we first collect the output in a string contained in a
|std::ostringstream|, and then let |ioutils::foldLine| break the result across
different lines after commas if necessary.

@h <sstream>
@h "ioutils.h"
@< Print the gradings for the part of |pi|... @>=
{ bool first=true; std::ostringstream os;
  for (size_t i=0; i<pi.size(); ++i)
    if ( rf_nr[pi(i)] == rf->val.realForm())
    { os << ( first ? first=false,'[' : ',');
      gradings::Grading gr=cc->val.fiber().grading(i);
      gr.permute(sigma);
      prettyprint::prettyPrint(os,gr,si.size());
    }
  os << "]" << std::endl;
  ioutils::foldLine(*output_stream,os.str(),"",",");
}

@* Test functions.
%
Now we shall make available some commands without actually creating new data
types. This means the values created to perform the computation will be
discarded after the output is produced; our functions will typically return a
$0$-tuple. This is probably not a desirable state of affairs, but the current
Atlas software is functioning like this (the corresponding commands are
defined in \.{realmode.cpp} and in \.{test.cpp} which hardly has any
persistent data) so for the moment it should not be so bad.

@ The \.{realweyl} and \.{strongreal} commands require a real form and a
compatible Cartan class.

@h "realredgp_io.h"
@< Local function def...@>=
void print_realweyl_wrapper()
{ push_tuple_components();
  Cartan_class_ptr cc(get<Cartan_class_value>());
  real_form_ptr rf(get<real_form_value>());
@)
  if (&rf->parent.val!=&cc->parent.val)
    throw std::runtime_error @|
    ("realweyl: inner class mismatch between arguments");
  bitmap::BitMap b(rf->parent.val.cartanSet(rf->val.realForm()));
  if (not b.isMember(cc->number))
    throw std::runtime_error @|
    ("realweyl: Cartan class not defined for real form");
@)
  realredgp_io::printRealWeyl (*output_stream,rf->val,cc->number);
@)
  wrap_tuple(0);
}

@)
void print_strongreal_wrapper()
{ push_tuple_components();
  Cartan_class_ptr cc(get<Cartan_class_value>());
  real_form_ptr rf(get<real_form_value>());
@)
  if (&rf->parent.val!=&cc->parent.val)
    throw std::runtime_error @|
    ("strongreal: inner class mismatch between arguments");
  bitmap::BitMap b(rf->parent.val.cartanSet(rf->val.realForm()));
  if (not b.isMember(cc->number))
    throw std::runtime_error @|
    ("strongreal: Cartan class not defined for real form");
@)
  realredgp_io::printStrongReal
    (*output_stream,rf->val,rf->parent.interface,cc->number);
@)
  wrap_tuple(0);
}


@
We shall next implement the \.{block} command.

@h "blocks.h"
@h "block_io.h"

@< Local function def...@>=
void print_block_wrapper()
{ push_tuple_components();
  dual_real_form_ptr drf(get<dual_real_form_value>());
  real_form_ptr rf(get<real_form_value>());
@)
  rf->val.fillCartan();
  if (&rf->parent.val!=&drf->parent.val)
    throw std::runtime_error @|
    ("inner class mismatch between real form and dual real form");
  bitmap::BitMap b(rf->parent.val.dualCartanSet(drf->val.realForm()));
  if (!b.isMember(rf->val.mostSplit()))
    throw std::runtime_error @|
    ("print_block: real form and dual real form are incompatible");
@)
  blocks::Block block(rf->parent.val
                     ,rf->val.realForm(),drf->val.realForm());
  block_io::printBlock(*output_stream,block);
@)
  wrap_tuple(0);
}

@ We provide functions corresponding to the \.{blockd} and \.{blocku}
variations of \.{block}.

@< Local function def...@>=
void print_blockd_wrapper()
{ push_tuple_components();
  dual_real_form_ptr drf(get<dual_real_form_value>());
  real_form_ptr rf(get<real_form_value>());
@)
  rf->val.fillCartan();
  if (&rf->parent.val!=&drf->parent.val)
    throw std::runtime_error @|
    ("inner class mismatch between real form and dual real form");
  bitmap::BitMap b(rf->parent.val.dualCartanSet(drf->val.realForm()));
  if (!b.isMember(rf->val.mostSplit()))
    throw std::runtime_error @|
    ("print_blockd: real form and dual real form are incompatible");
@)
  blocks::Block block(rf->parent.val
                     ,rf->val.realForm(),drf->val.realForm());
  block_io::printBlockD(*output_stream,block);
@)
  wrap_tuple(0);
}

@)
void print_blocku_wrapper()
{ push_tuple_components();
  dual_real_form_ptr drf(get<dual_real_form_value>());
  real_form_ptr rf(get<real_form_value>());
@)
  rf->val.fillCartan();
  if (&rf->parent.val!=&drf->parent.val)
    throw std::runtime_error @|
    ("inner class mismatch between real form and dual real form");
  bitmap::BitMap b(rf->parent.val.dualCartanSet(drf->val.realForm()));
  if (!b.isMember(rf->val.mostSplit()))
    throw std::runtime_error @|
    ("print_blocku: real form and dual real form are incompatible");
@)
  blocks::Block block(rf->parent.val
                     ,rf->val.realForm(),drf->val.realForm());
  block_io::printBlockU(*output_stream,block);
@)
  wrap_tuple(0);
}

@ The \.{blockstabilizer} command has a slightly different calling scheme than
\.{block} and its friends, in that it requires a real and a dual real form,
and a Cartan class; on the other hand it is simpler in not requiring a block
to be constructed. The signature of |realredgp_io::printBlockStabilizer| is a
bit strange, as it requires a |realredgp::RealReductiveGroup| argument for the
real form, but only numbers for the Cartan class and the dual real form (but
this is understandable, as information about the inner class must be
transmitted in some way). In fact it used to be even a bit stranger, in that
the real form was passed in the form of a |realredgp_io::Interface| value, a
class (not to be confused with |realform_io::Interface|, which does not
specify a particular real form) that we do not use in this program; since only
the |realGroup| field of the |realredgp_io::Interface| was used in
|realredgp_io::printBlockStabilizer|, we have changed its parameter
specification to allow it to be called easily here.

@< Local function def...@>=
void print_blockstabilizer_wrapper()
{ push_tuple_components();
  Cartan_class_ptr cc(get<Cartan_class_value>());
  dual_real_form_ptr drf(get<dual_real_form_value>());
  real_form_ptr rf(get<real_form_value>());
@)
  if (&rf->parent.val!=&drf->parent.val or
      &rf->parent.val!=&cc->parent.val)
    throw std::runtime_error @|
    ("blockstabilizer: inner class mismatch between arguments");
  bitmap::BitMap b(rf->parent.val.cartanSet(rf->val.realForm()));
  b &= bitmap::BitMap(rf->parent.val.dualCartanSet(drf->val.realForm()));
  if (not b.isMember(cc->number))
    throw std::runtime_error @|
    ("blockstabilizer: Cartan class not defined for both real forms");
@)

  realredgp_io::printBlockStabilizer
   (*output_stream,rf->val,cc->number,drf->val.realForm());
@)
  wrap_tuple(0);
}

@ The function |print_KGB| takes only a real form as argument.

@h "kgb.h"
@h "kgb_io.h"

@< Local function def...@>=
void print_KGB_wrapper()
{ real_form_ptr rf(get<real_form_value>());
@)
  rf->val.fillCartan();
@)
  *output_stream
    << "kgbsize: " << rf->val.kgbSize() << std::endl;
  kgb::KGB kgb(rf->val);
  kgb_io::var_print_KGB(*output_stream,rf->val.complexGroup(),kgb);
@)
  wrap_tuple(0);
}

@ The function |print_KL_basis| behaves much like |print_block| as far as
parametrisation is concerned.

@h "kl.h"
@h "klsupport.h"
@h "kl_io.h"
@< Local function def...@>=
void print_KL_basis_wrapper()
{ push_tuple_components();
  dual_real_form_ptr drf(get<dual_real_form_value>());
  real_form_ptr rf(get<real_form_value>());
@)
  rf->val.fillCartan();
  if (&rf->parent.val!=&drf->parent.val)
    throw std::runtime_error @|
    ("inner class mismatch between real form and dual real form");
  bitmap::BitMap b(rf->parent.val.dualCartanSet(drf->val.realForm()));
  if (!b.isMember(rf->val.mostSplit()))
    throw std::runtime_error @|
    ("print_KL_basis: real form and dual real form are incompatible");
@)
  blocks::Block block(rf->parent.val
                     ,rf->val.realForm(),drf->val.realForm());
  klsupport::KLSupport kls(block); kls.fill();

  kl::KLContext klc(kls); klc.fill();
@)
  *output_stream
    << "Full list of non-zero Kazhdan-Lusztig-Vogan polynomials:\n\n";
  kl_io::printAllKL(*output_stream,klc);
@)
  wrap_tuple(0);
}

@ The function |print_prim_KL| is a variation of |print_KL_basis|.

@< Local function def...@>=
void print_prim_KL_wrapper()
{ push_tuple_components();
  dual_real_form_ptr drf(get<dual_real_form_value>());
  real_form_ptr rf(get<real_form_value>());
@)
  rf->val.fillCartan();
  if (&rf->parent.val!=&drf->parent.val)
    throw std::runtime_error @|
    ("inner class mismatch between real form and dual real form");
  bitmap::BitMap b(rf->parent.val.dualCartanSet(drf->val.realForm()));
  if (!b.isMember(rf->val.mostSplit()))
    throw std::runtime_error @|
    ("print_prim_KL: real form and dual real form are incompatible");
@)
  blocks::Block block(rf->parent.val
                     ,rf->val.realForm(),drf->val.realForm());
  klsupport::KLSupport kls(block); kls.fill();

  kl::KLContext klc(kls); klc.fill();
@)
  *output_stream
    << "Non-zero Kazhdan-Lusztig-Vogan polynomials for primitive pairs:\n\n";
  kl_io::printPrimitiveKL(*output_stream,klc);
@)
  wrap_tuple(0);
}

@ The function |print_KL_list| is another variation of |print_KL_basis|, it
outputs just a list of all distinct Kazhdan-Lusztig-Vogan polynomials.

@< Local function def...@>=
void print_KL_list_wrapper()
{ push_tuple_components();
  dual_real_form_ptr drf(get<dual_real_form_value>());
  real_form_ptr rf(get<real_form_value>());
@)
  rf->val.fillCartan();
  if (&rf->parent.val!=&drf->parent.val)
    throw std::runtime_error @|
    ("inner class mismatch between real form and dual real form");
  bitmap::BitMap b(rf->parent.val.dualCartanSet(drf->val.realForm()));
  if (!b.isMember(rf->val.mostSplit()))
    throw std::runtime_error @|
    ("print_KL_list: real form and dual real form are incompatible");
@)
  blocks::Block block(rf->parent.val
                     ,rf->val.realForm(),drf->val.realForm());
  klsupport::KLSupport kls(block); kls.fill();

  kl::KLContext klc(kls); klc.fill();
@)
  kl_io::printKLList(*output_stream,klc);
@)
  wrap_tuple(0);
}

@ We close with two functions for printing the $W$-graph determined by the
polynomials computed. For |print_W_cells| we must construct one more object,
after having built the |klc::KLContext|.

@h "wgraph.h"
@h "wgraph_io.h"

@< Local function def...@>=
void print_W_cells_wrapper()
{ push_tuple_components();
  dual_real_form_ptr drf(get<dual_real_form_value>());
  real_form_ptr rf(get<real_form_value>());
@)
  rf->val.fillCartan();
  if (&rf->parent.val!=&drf->parent.val)
    throw std::runtime_error @|
    ("inner class mismatch between real form and dual real form");
  bitmap::BitMap b(rf->parent.val.dualCartanSet(drf->val.realForm()));
  if (!b.isMember(rf->val.mostSplit()))
    throw std::runtime_error @|
    ("print_W_cells: real form and dual real form are incompatible");
@)
  blocks::Block block(rf->parent.val
                     ,rf->val.realForm(),drf->val.realForm());
  klsupport::KLSupport kls(block); kls.fill();

  kl::KLContext klc(kls); klc.fill();

  wgraph::WGraph wg(klc.rank()); kl::wGraph(wg,klc);
@)
  wgraph_io::printCells(*output_stream,wg);
@)
  wrap_tuple(0);
}

@ And as last function for the moment, |print_W_graph| just gives a variation
of the output routine of |print_W_cells|.

@< Local function def...@>=
void print_W_graph_wrapper()
{ push_tuple_components();
  dual_real_form_ptr drf(get<dual_real_form_value>());
  real_form_ptr rf(get<real_form_value>());
@)
  rf->val.fillCartan();
  if (&rf->parent.val!=&drf->parent.val)
    throw std::runtime_error @|
    ("inner class mismatch between real form and dual real form");
  bitmap::BitMap b(rf->parent.val.dualCartanSet(drf->val.realForm()));
  if (!b.isMember(rf->val.mostSplit()))
    throw std::runtime_error @|
    ("print_W_graph: real form and dual real form are incompatible");
@)
  blocks::Block block(rf->parent.val
                     ,rf->val.realForm(),drf->val.realForm());
  klsupport::KLSupport kls(block); kls.fill();

  kl::KLContext klc(kls); klc.fill();

  wgraph::WGraph wg(klc.rank()); kl::wGraph(wg,klc);
@)
  wgraph_io::printWGraph(*output_stream,wg);
@)
  wrap_tuple(0);
}

@ Finally we install everything (where did we hear that being said before?)

@< Install wrapper functions @>=
install_function(real_form_wrapper,"real_form","(InnerClass,int->RealForm)");
install_function(quasisplit_form_wrapper,"quasisplit_form"
		,"(InnerClass->RealForm)");
install_function(components_rank_wrapper,"components_rank","(RealForm->int)");
install_function(count_Cartans_wrapper,"count_Cartans","(RealForm->int)");
install_function(KGB_size_wrapper,"KGB_size","(RealForm->int)");
install_function(Cartan_order_matrix_wrapper,"Cartan_order_matrix"
					    ,"(RealForm->mat)");
install_function(dual_real_form_wrapper,@|"dual_real_form"
				       ,"(InnerClass,int->DualRealForm)");
install_function(dual_quasisplit_form_wrapper,@|"dual_quasisplit_form"
		,"(InnerClass->DualRealForm)");
install_function(real_form_from_dual_wrapper,@|"real_form_from_dual"
				  ,"(DualRealForm->RealForm)");
install_function(Cartan_class_wrapper,"Cartan_class"
		,"(RealForm,int->CartanClass)");
install_function(most_split_Cartan_wrapper,"most_split_Cartan"
		,"(RealForm->CartanClass)");
install_function(print_Cartan_info_wrapper,"print_Cartan_info"
		,"(CartanClass->)");
install_function(real_forms_of_Cartan_wrapper,"real_forms_of_Cartan"
		,"(CartanClass->[RealForm])");
install_function(dual_real_forms_of_Cartan_wrapper,"dual_real_forms_of_Cartan"
		,"(CartanClass->[DualRealForm])");
install_function(fiber_part_wrapper,"fiber_part"
		,"(CartanClass,RealForm->[int])");
install_function(print_gradings_wrapper,"print_gradings"
		,"(CartanClass,RealForm->)");
install_function(print_realweyl_wrapper,"print_real_Weyl"
		,"(RealForm,CartanClass->)");
install_function(print_strongreal_wrapper,"print_strong_real"
		,"(RealForm,CartanClass->)");
install_function(print_block_wrapper,"print_block"
		,"(RealForm,DualRealForm->)");
install_function(print_blocku_wrapper,"print_blocku"
		,"(RealForm,DualRealForm->)");
install_function(print_blockd_wrapper,"print_blockd"
		,"(RealForm,DualRealForm->)");
install_function(print_blockstabilizer_wrapper,@|"print_blockstabilizer"
		,"(RealForm,DualRealForm,CartanClass->)");
install_function(print_KGB_wrapper,"print_KGB"
		,"(RealForm->)");
install_function(print_KL_basis_wrapper,"print_KL_basis"
		,"(RealForm,DualRealForm->)");
install_function(print_prim_KL_wrapper,"print_prim_KL"
		,"(RealForm,DualRealForm->)");
install_function(print_KL_list_wrapper,"print_KL_list"
		,"(RealForm,DualRealForm->)");
install_function(print_W_cells_wrapper,"print_W_cells"
		,"(RealForm,DualRealForm->)");
install_function(print_W_graph_wrapper,"print_W_graph"
		,"(RealForm,DualRealForm->)");




@* Index.

% Local IspellDict: british
