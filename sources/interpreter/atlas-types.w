% Copyright (C) 2006-2017 Marc van Leeuwen
% This file is part of the Atlas of Lie Groups and Representations (the Atlas)

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
\def\Zee{{\bf Z}} % cwebx uses \Z and \ZZ itself
\def\Qu{{\bf Q}}
\def\ii{{\rm i}} % imaginary unit

@* Atlas types.
%
This file describes several built-in types related to the Atlas software,
which used by the interpreter as primitive types. It also defines built-in
functions that relate to these types. This part of the program is completely
oriented towards the Atlas library, whereas the module \.{axis} implements the
programming language of that name, which is independent of the Atlas library.

@h "atlas-types.h"

@f lambda nullptr
@f pi nullptr
@f alpha nullptr
@f beta nullptr

@c
namespace atlas { namespace interpreter {
@< Global variable definitions @>
namespace {@; @< Local function definitions @>@; }@;
@< Function definitions @>@;
}@; }@;

@ As usual the external interface is written to the header file associated to
this file.

@( atlas-types.h @>=

#ifndef ATLAS_TYPES_H
#define ATLAS_TYPES_H

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
avoids having to declare each wrapper function. We also install all coercions,
and it is important that we do this \emph{before} installing the wrapper
functions, since the presence of coercions should have an effect on the order
in which overloads for the same name are ordered, so the former should be
present before the latter are processed.

@< Function definitions @>=
void initialise_builtin_types()
{ @< Install coercions @>
  @< Install wrapper functions @>
}

@ Before we can define any types we must make sure the types defined
in \.{axis-types.w} from which we shall derive others are known. Including
this as first file from our header file ensures the types are known wherever
they are needed. Some more basic built-in types likes integers, vectors, and
strings are defined in \.{global.w}, and we need their declarations too in our
implementation, but we avoid including its header into out header file.

@h "global.h"

@< Includes needed in the header file @>=
#include "../Atlas.h" // must be very first \.{atlas} include
#include "axis-types.h"

@*1 Lie types.
Our first chapter concerns Lie types, as indicated by strings like
|"D4.A3.E8.T1"|. We base ourselves on the related types defined in the Atlas
library.

@< Includes needed in the header file @>=
#include <stdexcept>
#include "lietype.h"

@*2 The primitive type.
%
A first new type corresponds to the type |LieType| in the Atlas library. We
provide a constructor that incorporates a complete |LieType| value, but often
one will use the default constructor and the |add| method. The remaining methods
are obligatory for a primitive type, though the (private) copy constructor is
only needed because |clone| uses it, and this in turn is only needed because we
use |get_own<Lie_type_value>| below; certain other primitive types do not expose
any functions that are implemented by modifying and existing copy, and can
implement |clone| that simply returns~|nullptr| (wince these calls come in a
context where the actual type is known, using a virtual method |clone| here
probably was a sub-optimal design decision).

@< Type definitions @>=
struct Lie_type_value : public value_base
{ LieType val;
@)
  Lie_type_value() : val() @+ {}
    // default constructor, produces empty type
  Lie_type_value(LieType t) : val(t) @+{}
    // constructor from already validated Lie type
@)
  virtual void print(std::ostream& out) const;
  Lie_type_value* clone() const @+{@; return new Lie_type_value(*this); }
  static const char* name() @+{@; return "Lie type"; }
@)
  void add_simple_factor (char,unsigned int); // grow
private:
  Lie_type_value(const Lie_type_value& v) : val(v.val) @+{}
    // copy constructor, used by |clone|
};
@)
typedef std::shared_ptr<const Lie_type_value> shared_Lie_type;
typedef std::shared_ptr<Lie_type_value> own_Lie_type; // used during construction

@ Before we do anything more complicated with this primitive type, we must
ensure that we can print its values. We can use an operator defined in
\.{basic\_io.cpp}.

@h "basic_io.h"
@< Function definitions @>=
void Lie_type_value::print(std::ostream& out) const
{ if (val.empty()) out << "empty Lie type";
  else
    out << "Lie type '" << val << '\'';
}

@ The type |LieType| is publicly derived from |std::vector<SimpleLieType>|, and
in its turn the type |SimpleLieType| is publicly derived from
|std::pair<char,size_t>|. Therefore these types could in principle take
arbitrary values, not necessarily meaningful ones. To ensure that this cannot
happen to \.{atlas} users, we make the method |add_simple_factor|, which is
invoked to build up Lie types, checks for the validity.

Since the tests defined in \.{io/interactive\_lietype.cpp} used in \.{Fokko} are
clumsy to use, we prefer to perform our own tests here, which emulate
|interactive_lietype::checkSimpleLieType|. One can specify torus factors of
rank~$r>1$, but they are equivalent to $r$ torus factors of rank~$1$, and it
simplifies the software if we rewrite the former form to the latter on input, so
that is what we do here.

@h "constants.h"

@< Function definitions @>=
void Lie_type_value::add_simple_factor (char c,unsigned int rank)
{ static const std::string types=lietype::typeLetters; // |"ABCDEFGT"|
  auto t=types.find(c);
  if (t==std::string::npos)
    throw runtime_error() << "Invalid type letter '" << c << '\'';
@.Invalid type letter@>
  const unsigned int r=constants::RANK_MAX; // for convenience
@/static const unsigned int lwb[]={1,2,2,4,6,4,2,0};
  static const unsigned int upb[]={r,r,r,r,8,4,2,r};
  if (rank<lwb[t])
    throw runtime_error()
    << "Too small rank " << rank << " for Lie type " << c;
@.Too small rank@>
  if (rank>upb[t])
  { if (upb[t]!=r)
      throw runtime_error()
      << "Too large rank " << rank << " for Lie type " << c;
@.Too large rank@>
    else
      throw runtime_error() @|
       << "Rank "+str(rank)+" exceeds implementation limit " << r;
@.Rank exceeds implementation limit@>
  }
@)
  if ((val.rank() + rank)>r)
    throw runtime_error() << "Total rank exceeds implementation limit " << r;
@.Total rank exceeds...@>
@)
  if (c=='T')
    while (rank-->0) val.push_back(SimpleLieType('T',1));
  else
    val.push_back(SimpleLieType(c,rank));
}

@ The function |Lie_type_wrapper| constructs a |Lie_type_value|. We scan the
string looking for sequences of a letter followed by a number, ignoring
punctuation characters between simple factors. Also spaces are allowed anywhere
except inside the number, since formatted input from streams, even of
characters, by default skips spaces. Since the correct structure of Lie type
strings is so obvious to the human eye, our error message just cites the entire
offending string, rather than trying to point out the error exactly.

@h <sstream>
@h <cctype>
@< Local function definitions @>=
inline std::istream& skip_punctuation(std::istream &is)
{@; char c; do is>>c; while (is and std::ispunct(c));
    return is.unget();
}
@)
void Lie_type_wrapper(expression_base::level l)
{ std::istringstream is(get<string_value>()->val);
  own_Lie_type result = std::make_shared<Lie_type_value>();
  char c;
  while (skip_punctuation(is)>>c) // i.e., until |not is.good()|
  { size_t rank;
    if (is.peek()=='-' or not (is>>rank)) // explicitly forbid minus sign
      throw runtime_error ()
        << "Error in type string '" << is.str() << "' for Lie type";
@.Error in type string@>
    result->add_simple_factor(c,rank);
      // this may |throw| a |runtime_error| as well
  }
  if (l!=expression_base::no_value)
    push_value(std::move(result));
}
@)
void Lie_type_coercion()
@+{@; Lie_type_wrapper(expression_base::single_value); }

@ Other useful ways of building Lie types is by combining two of them to a
single one, or in a more incremental fashion by adding a single (type,rank)
pair, which corresponds directly to |add_simple_factor|.

The curious |static_cast<unsigned>| below serves to force creating a temporary
from this constant rather than trying to pass it be reference, as the latter
would require memory to be reserved for it permanently.

@< Local function definitions @>=
void compose_Lie_types_wrapper(expression_base::level l)
{
  shared_Lie_type t2=get<Lie_type_value>();
  own_Lie_type t1=get_own<Lie_type_value>();
  if (t1->val.rank()+t2->val.rank()>constants::RANK_MAX)
      throw runtime_error() @|
        << "Combined rank exceeds implementation limit " @|
        << static_cast<unsigned>(constants::RANK_MAX);
@.Combined rank exceeds...@>
  if (l==expression_base::no_value)
    return;
@)
  t1->val.append(t2->val); // both factors have been tested, so this is safe
  push_value(t1);
}
@)
void extend_Lie_type_wrapper(expression_base::level l)
{
  auto rank = get<int_value>()->int_val();
  auto type_string = get<string_value>()->val;
  char type_letter = type_string.empty() ? 'T' : type_string[0];
  own_Lie_type t=get_own<Lie_type_value>();
  t->add_simple_factor(type_letter,rank);
  if (l!=expression_base::no_value)
    push_value(t);
}

@ We have predicates for testing (in)equality of Lie types. The actual
comparison is done by generic equality operations defined in the \Cpp\ standard
library for |std::vector| values (of equal types).

@< Local function definitions @>=
void Lie_type_eq_wrapper (expression_base::level l)
{ shared_Lie_type lt1 = get<Lie_type_value>();
  shared_Lie_type lt0 = get<Lie_type_value>();
  if (l!=expression_base::no_value)
    push_value(whether(lt0->val==lt1->val));
}
void Lie_type_neq_wrapper (expression_base::level l)
{ shared_Lie_type lt1 = get<Lie_type_value>();
  shared_Lie_type lt0 = get<Lie_type_value>();
  if (l!=expression_base::no_value)
    push_value(whether(lt0->val!=lt1->val));
}


@*2 Auxiliary functions for Lie types.
%
Here is a function that computes the Cartan matrix for a given Lie type. Unlike
Cartan matrices for root data (which just give the pairings of simple roots and
simple coroots), this one is affected by any central torus factors, which will
lead to entirely zero rows and columns.

@h "prerootdata.h"
@< Local function definitions @>=
void Cartan_matrix_wrapper(expression_base::level l)
{ shared_Lie_type t=get<Lie_type_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value>(t->val.Cartan_matrix()));
}


@ And here is a function that tries to do more or less the inverse. However,
since this is intended for Cartan matrices not directly obtained from a Lie type
but rather those computed from a root (sub)system, we use a version that
currently refuses any zero rows and columns such as the |LieType::Cartan_matrix|
method would produce for torus factors (since |dynkin::Lie_type| refuses them),
but on the other hand will recognise ``permuted'' Cartan matrices, not following
the Bourbaki numbering of nodes in a Dynkin diagram. Indeed the permutation
found is exported as a second component of the result, of type \.{[int]}.

We call |dynkin::lieType| in its version that also produces a permutation |pi|
(the one that maps the standard ordering of the diagram to the actual ordering),
and |true| as third argument to request a check that the input was indeed the
Cartan matrix for that situation (if not a runtime error will be thrown).
Without this test invalid Lie types could have been formed, for which root datum
construction would most likely crash.

@h "dynkin.h"
@< Local function definitions @>=
void type_of_Cartan_matrix_wrapper (expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  Permutation pi;
  LieType lt=dynkin::Lie_type(m->val,true,true,pi);
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<Lie_type_value>(lt));
  own_row perm = std::make_shared<row_value>(0);
  perm->val.reserve(pi.size());
  for(auto it=pi.begin(); it!=pi.end(); ++it)
    perm->val.push_back(std::make_shared<int_value>(*it));
  push_value(perm);
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ For programming it is important to be able to analyse a Lie type. To this end
we allow transforming it into a list of $(code,rank)$ pairs, where |code| is a
one-letter string.

@< Local function definitions @>=
void Lie_factors_wrapper(expression_base::level l)
{ shared_Lie_type t=get<Lie_type_value>();
  if (l==expression_base::no_value)
    return;
@)
  own_row result = std::make_shared<row_value>(t->val.size());
  for (unsigned i=0; i<t->val.size(); ++i)
  { const auto& src = t->val[i];
    auto dst=std::make_shared<tuple_value>(2);
    dst->val[0] = std::make_shared<string_value>(std::string(1,src.first));
    dst->val[1] = std::make_shared<int_value>(src.second);
    result->val[i] = std::move(dst);
  }
  push_value(std::move(result));
}

@ We now install all wrapper functions directly associated to Lie types.

@< Install wrapper functions @>=
install_function(Lie_type_wrapper,"Lie_type","(string->LieType)");
install_function(compose_Lie_types_wrapper,"*","(LieType,LieType->LieType)");
install_function(extend_Lie_type_wrapper,@|"extend"
                ,"(LieType,string,int->LieType)");
install_function(Lie_type_eq_wrapper,@|"=","(LieType,LieType->bool)");
install_function(Lie_type_neq_wrapper,@|"!=","(LieType,LieType->bool)");
@)
install_function(Cartan_matrix_wrapper,"Cartan_matrix","(LieType->mat)");
install_function(type_of_Cartan_matrix_wrapper
		,@|"Cartan_matrix_type","(mat->LieType,[int])");
install_function(Lie_factors_wrapper,"%","(LieType->[string,int])");

@*2 Finding lattices for a given Lie type.
%
When a Lie type is fixed, there is still a nontrivial choice to determine the
root datum for a connected complex reductive group of that type: one has to
choose a sub-lattice of the weight lattice of the ``simply connected'' group
of that type to become the weight lattice of the chosen Complex reductive
group (a finite quotient of the simply connected group); that sub-lattice
should have full rank and contain the root lattice. We shall start with
considering some auxiliary functions to help facilitate this choice.

Here is a wrapper function around the |LieType::Smith_basis| method, provided to
match the behaviour of \.{Fokko} in \.{atlas}. It computes a block-wise Smith
basis for the transposed Cartan matrix, and the corresponding ``block invariant
factors''. In case of torus factors this description should be interpreted in
the sense that the Smith basis for those factors is the standard basis and the
invariant factors are null.

@< Local function definitions @>=
void Smith_Cartan_wrapper(expression_base::level l)
{ shared_Lie_type t=get<Lie_type_value>();
  if (l==expression_base::no_value)
    return;
@)
  own_vector inv_factors = std::make_shared<vector_value>(CoeffList());
  push_value(std::make_shared<matrix_value>
    (t->val.Smith_basis(inv_factors->val)));
  push_value(std::move(inv_factors));
  if (l==expression_base::single_value)
     wrap_tuple<2>();
}

@ The result of |Smith_Cartan| can serve among other things to help specifying
a basis for subgroup of the quotient of the original weight
lattice~$\tilde{X}$ by the sub-lattice spanned by the roots (the columns of
the transposed Cartan matrix). For that purpose only invariant factors unequal
to~$1$ and the corresponding columns of the Smith basis are of interest.
Therefore the following wrapper function, which may operate directly on the
result of the previous one, filters out the invariant factors~$1$ and the
corresponding columns.

@< Local function definitions @>=
void filter_units_wrapper (expression_base::level l)
{ own_vector inv_f=get_own<vector_value>();
  own_matrix basis=get_own<matrix_value>();
  unsigned int n=basis->val.numColumns();
  if (inv_f->val.size()>n)
    throw runtime_error @|("Too many factors: "+
@.Too many factors@>
      str(inv_f->val.size()) + " for " +str(n)+ " columns");
  if (l==expression_base::no_value)
    return;
@)
  size_t i=0;
  while (i<n)
    if (i>=inv_f->val.size() or inv_f->val[i]!=1)
      ++i; // keep invariant factor and column
    else
    {@; inv_f->val.erase(inv_f->val.begin()+i);
        basis->val.eraseColumn(i);
    }
  push_value(basis); push_value(inv_f);
  if (l==expression_base::single_value)
      wrap_tuple<2>();
}

@ Here is another function, adapted from the functions |makeOrthogonal| and
|getLattice|, defined locally in the file \.{io/interactive\_lattice.cpp}. We
have taken the occasion to change the name and interface and everything else,
which also avoids the need to introduce rational matrices as primitive type.

The function |annihilator_modulo| takes as argument an $m\times{n}$
matrix~$M$, and an integer |d|. It returns a full rank (so invertible
over~$\Qu$) $m\times{m}$ matrix~|A| such that $A^t\cdot{M}$ is divisible
by~$d$, and whose column span is maximal subject to this constraint: the
columns of~$A$ span the full rank sub-lattice of $\Zee^m$ of vectors $v$ such
that $v^t\cdot{M}\in d\,\Zee^n$. The fact that $d$ is called |denominator|
below comes from the alternative interpretation that the transpose of~$A$
times the rational matrix $M/d$ gives an integral matrix as result.

The algorithm is quite simple. After finding invertible matrices |row|, |col|
such that $row*M*col$ is diagonal with non-zero diagonal entries given in
$\lambda$ (and any zero diagonal entries trailing those), we know that any row
of |row| with factor $\lambda_i$ is a linear form sending the image of $M$ to
$\lambda_i\Zee$, while any remaining rows of |row| (those without
corresponding diagonal entry) annihilate the image altogether. Then all that
is needed it to multiply rows of~|row| of the first kind by
$d/\gcd(d,\lambda_i)$ and transpose the result.

There is a subtlety though, that |matreduc::diagonalise| will not ensure that
the first coefficient in the list is positive, as it prefers to ensure that
the base change matrices |row| and |col| that it sets both have
determinant~$1$ (a detail that does not interest us here, but the function
does not know that); this may necessitate a negative coefficient, which if it
occurs will be the first one. Since the call to |arithmetic::div_gcd|
implicitly converts its second argument to unsigned, it is imperative that we
use the absolute value of this first coefficient.

@h "arithmetic.h"
@h "lattice.h"
@h "matrix.h"
@h "matreduc.h"

@< Local function definitions @>=
LatticeMatrix @|
annihilator_modulo(const LatticeMatrix& M, arithmetic::Denom_t denominator)

{ int_Matrix row,col;
  CoeffList lambda = matreduc::diagonalise(M,row,col);
  if (not lambda.empty())
    lambda[0]=std::abs(lambda[0]); // ensure all entries positive, for |div_gcd|

  for (size_t i=0; i<lambda.size(); ++i)
    row.rowMultiply(i,arithmetic::div_gcd(denominator,lambda[i]));

  return row.transposed();
}

@ The wrapper function is now particularly simple.

@< Local function definitions @>=
void ann_mod_wrapper(expression_base::level l)
{ int d=get<int_value>()->int_val();
  shared_matrix m=get<matrix_value>();
@)
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value>(annihilator_modulo(m->val,d)));
}

@ Next a simple administrative routine, introduced because we could not handle
matrices in our programming language at the time (and still necessary as long
as we want to define |quotient_basis| below as built-in function, since such a
function cannot use the \.{axis} interpreter). Having computed a new
lattice (possibly with the help of |ann_mod|), in the form of vectors to
replace those selected by |filter_units| from the result of |Smith_Cartan|,
one needs to make the replacement. The following function does this, taking as
its first argument the result of |Smith_Cartan|, and as second a matrix whose
columns are to be substituted. The invariant factors in the first argument
serve only to determine, by the place the non-unit ones, where the insertion
has to take place. In fact this is so simple that we define the wrapper
function directly.

@< Local function definitions @>=
void replace_gen_wrapper (expression_base::level l)
{ shared_matrix new_generators=get<matrix_value>();
  push_tuple_components(); // a pair as returned by \.{Smith\_Cartan}
  shared_vector inv_f=get<vector_value>();
  own_matrix generators=get_own<matrix_value>();
  unsigned int n=generators->val.numColumns();
  if (inv_f->val.size()>n)
    throw runtime_error @|("Too many factors: "+
@.Too many factors@>
      str(inv_f->val.size()) + " for " +str(n)+ " columns");
@)
  if (new_generators->val.numRows()!=generators->val.numRows())
    throw runtime_error("Column lengths do not match");
@.Column lengths do not match@>
@)
  size_t k=0; // index to replacement generators
  for (size_t j=0; j<n; ++j)
    if (j>=inv_f->val.size() or inv_f->val[j]!=1)
       // replace column |j| by |k| from |new_generators|
    { if (k>=new_generators->val.numColumns())
        throw runtime_error ("Not enough replacement columns");
@.Not enough replacement columns@>
      generators->val.set_column(j,new_generators->val.column(k));
      ++k;
    }
  if (k<new_generators->val.numColumns())
        throw runtime_error ("Too many replacement columns");
@.Too many replacement columns@>
  if (l!=expression_base::no_value)
    push_value(generators);
}

@ To emulate what is done in \.{Fokko}, we write a function that integrates some
of the previous ones. It is called as |quotient_basis(lt,L)| where $t$ is a Lie
type, and $L$ is a list of rational vectors, each giving a kernel generator as
would be entered in \.{Fokko}. In this interpretation each rational vector
implicitly gives a rational coweight, where each position either corresponds to
a generator of a finite cyclic factor of order $d$ of the centre, with the
rational entry $1/d$ in that position representing that generator, or to a
central $1$-parameter subgroup (torus factor) of the centre with any rational
entry in that position (interpreted modulo~$\Zee$) describing an element (of
finite order) of that subgroup. More precisely one is taking rational linear
combinations of certain coweights in the dual basis to a basis of the weight
lattice adapted to the root lattice; for those weights in the latter basis that
already lie in the root lattice, the corresponding coweight is omitted, since
there are no useful rational multiples of it anyway. This function computes
appropriate integer linear combinations of the relevant subset of the adapted
basis, corresponding to the given choice of kernel generators; it is the
sublattice (of the weight lattice) that they span that is ultimately of
interest.

Let |S==Smith_Cartan(lt)| and $(C,v)=\\{filter\_units}(S)$, then we find a
basis for the sub-lattice needed to build the root datum as follows. The
vector $v$, which describes the orders of the cyclic factors of the centre, is
only used to validate the values in~$L$: each vector in~$L$ should have the
same length as~$v$, and multiplication of corresponding entries should always
give an integer. Then a common denominator~$d$ is found and a matrix $M$ whose
columns form the numerators of the lists of~$L$ brought to the
denominator~$d$. The call $quotient\_basis(lt,L)$ then yields the result of
computing $replace\_gen(S,C*ann\_mod(M,d))$.

This wrapper function does most of its work by calling other wrapper
functions, using the evaluator stack many times. In one case we duplicate the
top of the stack, as~$S$ serves both in the call to $filter\_units$
(immediately) and in the to~$replace\_gen$ (at the end).

@< Local function definitions @>=
void quotient_basis_wrapper(expression_base::level l)
{ shared_row L=get<row_value>();
  // and leave Lie type on stack for $Smith\_Cartan$
  Smith_Cartan_wrapper(expression_base::multi_value);
  shared_value SC_basis = *(execution_stack.end()-2);
  shared_value invf = *(execution_stack.end()-1);
  wrap_tuple<2>(); // |replace_gen| wants as first argument a tuple
@/push_value(SC_basis);
  push_value(invf);
  filter_units_wrapper(expression_base::multi_value);
  shared_vector v=get<vector_value>();
  // and leave $C$ for call to $mm\_prod$
@)
  LatticeMatrix M(v->val.size(),L->length());
  arithmetic::Denom_t d=1;
  @< Compute common denominator |d| of entries in~$L$, and place converted
     denominators into the columns of~$M$; also test validity of entries
     against |v|, and |throw| a runtime error for invalid ones @>
  push_value(std::make_shared<matrix_value>(annihilator_modulo(M,d)));
@/mm_prod_wrapper(expression_base::single_value);
@/replace_gen_wrapper(l); // pass level parameter to final call
}

@ Each vector in |L| must have as many entries as |v|, and multiplying by the
corresponding entry of~$v$ should chase the denominator of each entry. In a
first pass we copy the numerator vectors to columns of~$M$, check divisibility
according to |v| and compute the common denominator |d|; in a second pass the
numerators are reduced modulo their original denominator, and then brought to
the new denominator~|d|.

@< Compute common denominator |d| of entries in~$L$... @>=
{ std::vector<arithmetic::Numer_t> denom(L->length());
  for (unsigned int j=0; j<L->length(); ++j)
  { const RatWeight& gen =
      force<rational_vector_value>(&*L->val[j])->val;
    denom[j] = gen.denominator();
    d=arithmetic::lcm(d,denom[j]);

    if (gen.size()!=v->val.size())
      throw runtime_error () @|
        << "Length mismatch for generator " << j << ": "@|
@.Length mismatch...@>
        << gen.size() << ':' << v->val.size();

    const auto& col=gen.numerator();
    for (unsigned int i=0; i<v->val.size(); ++i)
    { if (v->val[i]*col[i]%gen.denominator()!=0) // must use signed arithmetic!!
	throw runtime_error() << "Improper generator entry: "
@.Improper generator entry@>
         << col[i] << '/' << denom[j] @|
         << " not a multiple of 1/" << v->val[i];
      M(i,j) = @| arithmetic::remainder(col[i],denom[j]);
      // ``mod $\Zee$''; makes |M(i,j)| non-negative
    }
  }
// loop must end here to complete computation of |d|, then restart a new loop
@)
  for (unsigned int j=0; j<L->length(); ++j)
    // convert to common denominator |d|
  { arithmetic::Denom_t f=d/denom[j];
    for (unsigned int i=0; i<v->val.size(); ++i)
      M(i,j) *= f;
  }
}

@*2 Specifying inner classes. Now we move ahead a bit in the theory, from
functions that help in building root data to functions that help defining (inner
classes of) real forms. The first of such functions is |lietype::involution|,
which takes a Lie type and an |InnerClassType| (a vector of characters
describing the kind of involution wanted) and produces a matrix describing the
involution, defined on the weight lattice for the simply connected group of the
given type. That function supposes its arguments have already been checked for
validity and undergone some transformation; in \.{Fokko} this is done by
|checkInnerClass| and |readInnerClass| from the \.{io/interactive\_lietype}
compilation unit. This forces us to perform similar actions before calling
|lietype::involution|. We prefer not to use the functions from that compilation
unit, for the same reason we did not use |checkSimpleLieType| above. Therefore
we shall first define a function |checked_inner_class_type| that transforming a
string describing an inner class into |InnerClassType| that is guaranteed to be
valid if returned; the routine throws a |runtime_error| in case of problems.

@< Local function definitions @>=
InnerClassType checked_inner_class_type
  (const char* s, const LieType& lt)
{ static const std::string types(lietype::innerClassLetters);
    // |"Ccesu"|
  InnerClassType result; // initially empty
  std::istringstream is(s);
  char c;
  unsigned int i=0; // position in simple factors of Lie type |lt|
  while (skip_punctuation(is)>>c)
    @< Test the inner class letter |c|, and either push a corresponding type
       letter onto |result| while advancing~|i| by the appropriate amount, or
       throw a |runtime_error| @>
  if (i< lt.size()) throw runtime_error("Too few inner class symbols");
@.Too few inner class symbols@>
  return result;
}

@ The type letter |'C'| meaning ``Complex'' is the only one using more that
one Lie type, and both types used must be equal. Types |'c'| ``compact'' has
synonym |'e'| ``equal rank'', and means the class containing identity; these
cases are recorded as |'c'|. Type |'s'| ``split'' means the class containing
minus identity, which in some cases equals |'c'|, and similarly |'u'| meaning
``unequal rank'' is often equal to |'s'|; these cases are detailed later.

@< Test the inner class letter... @>=
{ if (types.find(c) == std::string::npos) throw runtime_error@|
    (std::string("Unknown inner class symbol `")+c+"'");
  if (i>= lt.size()) throw runtime_error("Too many inner class symbols");
@.Too many inner class symbols@>
  lietype::TypeLetter t = lt[i].type(); unsigned int r=lt[i].rank();
  if (c=='C') // complex inner class, which requires two equal simple factors
  { if (i+1>=lt.size() or lt[i+1]!=lt[i]) throw runtime_error @|
      ("Complex inner class needs two identical consecutive types");
@.Complex inner class needs...@>
    result.push_back('C');
    i+=2; // advance past both simple factors
  }
  else if (c=='c' or c=='e')
  {@; result.push_back('c'); ++i; } // synonyms
  else if (c=='s')
  {@; @< Contribute |'s'|, or |'c'| in case that means the same thing @>
    ++i;
  }
  else if (c=='u') // unequal rank, the only remaining possibility
  {@; @< Test and possibly transform the unequal rank inner class symbol @>
    ++i;
  }
}

@ The type letter |'s'| means the same thing as |'c'| for those simple types
where the longest Weyl group element acts as $-1$; these types are $A_1$,
$B_n$, $C_n$, $D_{2n}$, $E_7$, $E_8$, $F_4$ and $G_2$. In these cases we
replace the symbol by |'c'|, and in the remaining cases $A_n$ ($n>1$),
$D_{2n+1}$, $E_6$ and $T_1$ we leave the |'s'| as specified.

@< Contribute |'s'|, or |'c'| in case that means the same thing @>=
if (t=='A' and r>=2 or
    t=='D' and r%2!=0 or
    t=='E' and r==6 or
    t=='T')
  result.push_back('s');
else result.push_back('c');

@ The type letter |'u'| is mathematically only
meaningful for Lie types $A_n$ with $n\geq2$, any legal type $D_n$, or $E_6$.
Moreover in all cases except $D_n$ with $n$ even, this designation is
equivalent to |'s'|, and will be replaced by it. Hence any surviving type
letter |'u'| corresponds to a type of the form~$D_{2n}$.

@< Test and possibly transform... @>=
{ if (t=='D') result.push_back(r%2==0 ? 'u' : 's');
  else if (t=='A' and r>=2 or t=='E' and r==6 or t=='T') result.push_back('s');
  else throw runtime_error @|
    (std::string("Unequal rank class is meaningless for type ")+t+str(r));
@.Unequal rank class is meaningless...@>
}

@ Below we shall also need to pass a permutation argument supplied by the user
(to initialise a field of a |lietype::Layout| structure). This requires checking
and conversion from the \.{atlas} type \.{[int]} to a |Permutation| internal
type, which the following function does.

@< Local function def... @>=
Permutation checked_permutation(const std::vector<shared_value>& pi)
{
  auto n=pi.size();
  Permutation result(n); BitMap seen(n);
  auto rit = result.begin();
  for (auto it=pi.begin(); it!=pi.end(); ++it,++rit)
  { auto entry = force<int_value>(it->get())->val.ulong_val();
    if (entry>=pi.size())
      throw runtime_error() << "Permutation entry " << entry << " too big";
    if (seen.isMember(entry))
      throw runtime_error() << "Permutation has repeated entry " << entry;
    seen.insert(entry);
    *rit=entry;
  }
  return result;
}

@ The wrapper function around |lietype::involution| will take a Lie type, a
permutation, and a string of type letters and return a matrix describing the
involution designated by that string, expressed on the fundamental weight basis
for the simply connected group of that type and permutation (under the
correspondence that |type_of_Cartan_matrix| implements). Everything gets packed
into a |lietype::Layout| first.

@< Local function def... @>=
void basic_involution_wrapper(expression_base::level l)
{ shared_string str=get<string_value>();
  shared_row perm = get<row_value>();
  shared_Lie_type t=get<Lie_type_value>();
  if (perm->val.size()!=t->val.rank())
    throw runtime_error() @|
    << "Permutation size " << perm->val.size() << @| " does not match rank "
    << t->val.rank() << " of Lie type";
  if (l==expression_base::no_value)
    return;
@)

  lietype::Layout lo @|
    { t->val
    , checked_inner_class_type(str->val.c_str(),t->val)
    , checked_permutation(perm->val)
    } ;
@/push_value(std::make_shared<matrix_value> (lietype::involution(lo)));
}

@ The function just defined gives an involution on the basis of fundamental
weights for a simply connected group, or on the basis of the root lattice for
an adjoint group. For a general complex groups however, we need to transform
this involution to an involution of a given sub-lattice, assuming that this is
possible, in other words if that the sub-lattice is globally stable under the
involution specified. Therefore we now provide a function that in addition to
the Lie type takes a matrix specifying a sub-lattice as argument, and finally
a string specifying the inner class.

The fist two ingredients are also those used to construct a root datum, and one
might imagine replacing them by a root datum. The function |set_inner_class|
defined later will do that, but it has to accept some ambiguity in recovering
Lie type and sub-lattice from a root datum. Also it may find a permutation of
the simple roots with respect to the standard ordering of the diagram, and has
to deal with that additional generality. So unfortunately it can neither replace
nor call the function |based_involution| defined here; it will have to adapt
some of the work done here to its situation.


@< Local function def... @>=
void based_involution_wrapper(expression_base::level l)
{ shared_string s = get<string_value>();
@/shared_matrix basis = get<matrix_value>();
@/shared_Lie_type type = get<Lie_type_value>();
@)
  unsigned int r=type->val.rank();
  if (basis->val.numRows()!=r or basis->val.numRows()!=r)
    throw runtime_error @|
    ("Basis should be given by "+str(r)+'x'+str(r)+" matrix");
@.Basis should be given...@>
@)
  WeightInvolution inv=lietype::involution
        (type->val,checked_inner_class_type(s->val.c_str(),type->val));
  try
  { push_value(std::make_shared<matrix_value>(inv.on_basis(basis->val)));
    if (l==expression_base::no_value)
      execution_stack.pop_back(); // we needed testing, but not the result
  }
  catch (std::runtime_error&) // relabel |"Inexact integer division"|
  {@; throw runtime_error
      ("Inner class is not compatible with given lattice");
@.Inner class not compatible...@>
  }
}

@ All that remains is installing the wrapper functions.

@< Install wrapper functions @>=
install_function(Smith_Cartan_wrapper,"Smith_Cartan","(LieType->mat,vec)");
install_function(filter_units_wrapper,"filter_units","(mat,vec->mat,vec)");
install_function(ann_mod_wrapper,"ann_mod","(mat,int->mat)");
install_function(replace_gen_wrapper,"replace_gen",
		"((mat,vec),mat->mat)");
install_function(quotient_basis_wrapper
		,@|"quotient_basis","(LieType,[ratvec]->mat)");
install_function(basic_involution_wrapper,"involution",
		"(LieType,[int],string->mat)");
install_function(based_involution_wrapper,"involution",
		"(LieType,mat,string->mat)");

@*1 Root data.
%
We shall now introduce primitive type for root data. There will be a user type
to encapsulate a |RootDatum| object defined in the \.{rootdata} compilation
unit. We also use shared pointers.

@< Includes needed in the header file @>=
#include <memory> // for |std::shared_ptr|
#include "hashtable.h"
#include "rootdata.h"


@ Since root data are the starting point for a hierarchy of more and more
complicated value types, we shall identify identical root data when they are
constructed, so that the values can be shared internally. This sharing will
also help to similarly recognise when identical values are produced further up
the hierarchy. To this end we shall keep a hash table of bare root data, with
just enough information to match identical root data, and associate to each such
value seen a shared pointer to a |RootDatum| value. We can later also attach
more information such as a list of inner classes for this root datum that have
been seen, allowing us to avoid duplication at that level higher up. So here is
a type compatible with the |Hash_table| class template.

@< Type definitions @>=
struct root_datum_entry : public PreRootDatum
{ explicit root_datum_entry(PreRootDatum&& pre)
  : PreRootDatum(std::move(pre)) @+{}

typedef std::vector<root_datum_entry> Pooltype;
  size_t hashCode(size_t modulus) const;
  bool operator!=(const root_datum_entry& x) const
    {@; return not PreRootDatum::operator==(x); }
};

@ The |root_datum_value| has as main purpose to wrap a |RootDatum| object into a
value derived from |value_base|. However, we want to make sure that identical
root data always end up using the \emph{same} |root_datum_value|; this is partly
to avoid wasting storage, but more importantly to be able to also easily
identify identical inner classes and other values based on them. In order to
achieve this, we include static variables with a hash table for all bare root
data seen, and a vector of pointers to corresponding root data values. The test
for existing identical root data should be made prior to calling the
|root_datum_value| constructor, and we should never duplicate such a value, so
the |clone| method does nothing (and there is no reason it should ever actually
get called). A static method |build| takes care of creating a |root_datum_value|
from a |PreRootDatum|; it will either locate an existing one and return a shared
pointer, or create one using |std::make_shared| and the provided constructor. To
ensure that testing for duplicates always takes place, we want to ensure that
clients cannot call the constructor directly. But we cannot make it private,
because it will be called from |std::make_shared|; instead we give it an
additional argument of a local type |token| that only methods of our class (like
|build|) are able to provide. The |val| member is public, but |const|.

@< Type definitions @>=
class root_datum_value;
typedef std::shared_ptr<const root_datum_value> shared_root_datum;
@)
class root_datum_value : public value_base
{ struct token@+{}; // type passed to prove caller has private access;
public:
  const RootDatum val;
  static HashTable<root_datum_entry,unsigned short> hash;
  static std::vector<std::weak_ptr<const root_datum_value> > store;
@)
  root_datum_value(const PreRootDatum& v,token) : val(v) @+ {}
  static shared_root_datum build(PreRootDatum&& pre);
  shared_root_datum dual() const; // get dual datum through |build|
  virtual void print(std::ostream& out) const;
  root_datum_value* clone() const @+{@; return nullptr; }
    // we refuse to create duplicates
  static const char* name() @+{@; return "root datum"; }
};

@ We need to define the static members declared in the class definition.

@< Global variable definitions @>=
root_datum_entry::Pooltype root_data_pool;
HashTable<root_datum_entry,unsigned short> root_datum_value::hash
  (root_data_pool);
std::vector<std::weak_ptr<const root_datum_value> > root_datum_value::store;

@ We have a simple hash function that uses all information in a |PreRootDatum|.
@< Function definitions @>=
size_t root_datum_entry::hashCode(size_t modulus) const
{ size_t h= prefer_coroots() ? 1 : 0;
  for (unsigned int i=0; i<rank(); ++i)
    for (unsigned int j=0; j<semisimple_rank(); ++j)
      h=3*h+simple_root(j)[i];
  for (unsigned int i=0; i<rank(); ++i)
    for (unsigned int j=0; j<semisimple_rank(); ++j)
      h=3*h+simple_coroot(j)[i];
  return h&(modulus-1);
}

@ The static data member |root_datum_entry::store| has weak pointers, since
otherwise its very presence would prevent any |root_datum_value| to ever be
cleaned up when inaccessible. So |build| looks up the |PreRootDatum| passed to
it, and if it finds something (which means and identical |PreRootDatum| was here
before) and the weak pointer can be locked (which means the corresponding
|root_datum_value| is still in use) then it returns that locked pointer; if not
then a shared pointer to a newly created |root_datum_value| is returned, after
either pushing a weak pointer to |store| or storing it in the old slot in case a
the weak pointer there had expired.

The method |dual| facilitates building the dual of a |root_datum_value| while
passing though |build| to ensure that no unnecessary copies are created, notably
if one should later compute the dual of the dual.

@< Function definitions @>=
shared_root_datum root_datum_value::build(PreRootDatum&& pre)
{ auto loc = hash.match(root_datum_entry(std::move(pre)));
  if (loc<store.size())
  { if (auto result=store[loc].lock()) // previous root datum still exists
      return result; // so return it
  }
  auto result =
    std::make_shared<root_datum_value>(hash[loc],token());
    // construct |RootDatum|
  if (loc<store.size()) // happens if identical root datum was cleaned up
    store[loc]=result; // save a weak pointer version in |store|
  else // this must be the first time we ever see a copy of this |PreRootDatum|
@/{@; assert (loc==store.size());
    store.push_back(result);
  } // save a weak pointer version in |store|
  return result;
}
@)
shared_root_datum root_datum_value::dual() const
{@; PreRootDatum pre(val); pre.dualise(); return build(std::move(pre)); }

@*2 Printing root data.
%
We shall not print the complete information contained in the root datum. However
we do exercise some |RootDatum| methods to immediately give some superficial
information. The method |RootDatum::type| produces a complete type, including
central torus factors. We abuse conversion to |Lie_type_value| to ensure the
final part of our output is formatted like the output for a Lie type value.

The mention of |"simply connected"| and/or |"adjoint"| is in order to make some
use of these flags, which are present anyway. However, we only print something
for semisimple root data, since that is mathematically required in order for the
associated complex group to be either simply connected or adjoint. The proper
interpretation for non-semisimple root data would be ``has semisimple derived
group'' for |isSemisimple|, and ``has connected center'' for |isAdjoint|, but
these are quite a mouthful, especially if they need to be expressed in a
grammatically correct way.

@< Function definitions @>=
void root_datum_value::print(std::ostream& out) const
{ if (val.isSemisimple())
    out << (val.isSimplyConnected() ? "simply connected " : "")
     @| << (val.isAdjoint() ? "adjoint " : "");
  out << "root datum of " << Lie_type_value(val.type());
}

@ We also make the Lie type of a root datum (computed by |RootDatum::type|) and
whether the coroots rather than roots were used to determine the root numbering
(from the stored field) available by wrapper functions.
@< Local fun...@>=
void type_of_root_datum_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l!=expression_base::no_value)
    push_value(std::make_shared<Lie_type_value>(rd->val.type()));
}

void coroot_preference_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l!=expression_base::no_value)
    push_value(whether(rd->val.prefer_coroots()));
}

@*2 Building a root datum.
%
The most direct way to create a root datum value is to specify bases of simple
roots and coroots in the form of matrices, which implicitly define a Lie type
and weight lattice. In addition we take a Boolean argument telling whether to
use roots or coroots when generating the full root system from the simple roots
or coroots; this is an argument common to all functions building a fresh root
datum (by contrast those that bases a new root datum on an existing one will
just inherit the attribute). To make sure the root datum construction will
succeed, we must test the ``Cartan'' matrix computed from these data to be a
valid one.

@< Local function definitions @>=
void root_datum_wrapper(expression_base::level l)
{ bool prefer_coroots = get<bool_value>()->val;
@/shared_matrix simple_coroots=get<matrix_value>();
  shared_matrix simple_roots=get<matrix_value>();

  size_t nr = simple_roots->val.numRows(),
         nc = simple_roots->val.numColumns();

  if (simple_coroots->val.numRows()!=nr @| or
      simple_coroots->val.numColumns()!=nc)
    throw runtime_error
    ("Sizes (" +str(nr)+ "," +str(nc) +"),("+
      str(simple_coroots->val.numRows()) +","+
      str(simple_coroots->val.numColumns()) +
      ") of simple (co)root systems differ");
@.Sizes of simple (co)root systems...@>

try @/{
  PreRootDatum prd(simple_roots->val,simple_coroots->val,prefer_coroots);
  prd.test_Cartan_matrix();
  if (l!=expression_base::no_value)
    push_value(root_datum_value::build(std::move(prd)));
}
  catch (error::Cartan_error)
@/{@;
    throw runtime_error("Matrices of (co)roots give invalid Cartan matrix");
}
@.System of (co)roots has invalid...@>
}

@ Alternatively, the user may specify a Lie type and a square matrix of the size
of the rank of the root datum, which specifies generators of the desired weight
lattice as a sub-lattice of the lattice of weights associated to the simply
connected group of the type given. The given weights should be independent and
span at least the root lattice associated to the type. Failure of either
condition will cause |PreRootDatum::quotient| to throw a |std::runtime_error|.

@< Local function definitions @>=
void root_datum_from_type_wrapper(expression_base::level l)
{ bool prefer_coroots = get<bool_value>()->val;
@/shared_matrix lattice=get<matrix_value>();
  shared_Lie_type type=get<Lie_type_value>();
  if (lattice->val.numRows()!=lattice->val.numColumns() @| or
      lattice->val.numRows()!=type->val.rank())
    throw runtime_error
    ("Sub-lattice matrix should have size " @|
@.Sub-lattice matrix should...@>
      +str(type->val.rank())+'x'+str(type->val.rank()));
  PreRootDatum prd(type->val,prefer_coroots);
  prd.quotient(lattice->val);
@.Sub-lattice matrix not square...@>
@.Dependent lattice generators@>
@.Sub-lattice does not contain...@>
  if (l!=expression_base::no_value)
    push_value(root_datum_value::build(std::move(prd)));
}

@ While the previous function takes the sublattice to live in the weight
lattice of the simply connected root datum of the given type, one may more
generally wish in any existing root datum to reduce to a full-rank sublattice
of $X^*$ that contains the root lattice. The following variant of root datum
construction does this. The call to the |quotient| method may throw the same
errors as in the previous function.

@< Local function definitions @>=
void sublattice_root_datum_wrapper(expression_base::level l)
{ shared_matrix lattice=get<matrix_value>();
  shared_root_datum rd=get<root_datum_value>();
  const auto r = rd->val.rank();
  if (lattice->val.numRows()!=r or lattice->val.numColumns()!=r)
    throw runtime_error()
      << "Sub-lattice matrix should have size " @|
@.Sub-lattice matrix should...@>
      << r << 'x' << r;

  PreRootDatum prd = rd->val; // inherits |prefer_coroots| attribute
  prd.quotient(lattice->val); // this may |throw| a |std::runtime_error|
@.Sub-lattice does not contain...@>
@.Dependent lattice generators@>
  if (l!=expression_base::no_value)
    push_value(root_datum_value::build(std::move(prd)));
}

@ We define two more wrappers with only a Lie type as argument, for building the
simply connected and the adjoint root data. They mostly call other wrapper
functions: the matrices that specify the sub-lattices are produced by the
wrapper functions for |id_mat| respectively by those for |Cartan_matrix| and
|transpose_mat|. Thus $simply\_connected\_datum(lt)$ is equivalent to
$root\_datum(lt,id\_mat(Lie\_rank(lt)))$, while $adjoint\_datum(lt)$ is
equivalent to the call $root\_datum(lt,M)$ where $M$ essentially the transpose
matrix of |Cartan_matrix(lt)|; there is one modification that we do here ``by
hand'', namely making sure that all null diagonal entries of~$M$ (which must
come from torus factors) are replaced by ones.

@< Local function definitions @>=
void simply_connected_datum_wrapper(expression_base::level l)
{ bool prefer_coroots = get<bool_value>()->val;
  if (l==expression_base::no_value)
@/{@; execution_stack.pop_back();
    return;
  } // no possibilities of errors, so avoid useless work
@)
  auto rank =
    force<Lie_type_value>(execution_stack.back().get())->val.rank();
  push_value(std::make_shared<int_value>(rank));
  id_mat_wrapper(expression_base::single_value);
  push_value(whether(prefer_coroots));
@/root_datum_from_type_wrapper(expression_base::single_value);
}
@)
void adjoint_datum_wrapper(expression_base::level l)
{ bool prefer_coroots = get<bool_value>()->val;
  if (l==expression_base::no_value)
@/{@; execution_stack.pop_back();
    return;
  } // no possibilities of errors, so avoid useless work
@)
  push_value(execution_stack.back()); // duplicate Lie type argument
  Cartan_matrix_wrapper(expression_base::single_value);
  transpose_mat_wrapper(expression_base::single_value);
  own_matrix M=get_own<matrix_value>();
  for (size_t i=0; i<M->val.numRows(); ++i)
    if (M->val(i,i)==0) M->val(i,i)=1;
  push_value(M);
  push_value(whether(prefer_coroots));
@/root_datum_from_type_wrapper(expression_base::single_value);
}


@*2 Functions operating on root data.
%
The following functions allow us to look at individual simple roots and simple
coroots stored in a root datum value. We adopt a convention that shifts root
indices so that the simple (and positive) roots start at index~$0$, and such
that negative roots have negative indices. This allows the same function to be
used for producing simple, positive, or general roots. Since internally the
full list of roots starts at index~$0$ with the most negative roots, we must
apply a shift by the number of positive roots here.

@< Local function definitions @>=
void root_wrapper(expression_base::level l)
{ int root_index = get<int_value>()->int_val();
  shared_root_datum rd(get<root_datum_value>());
  RootNbr npr = rd->val.numPosRoots();
  RootNbr alpha = npr+root_index;
  if (alpha>=2*npr)
    throw runtime_error("Illegal root index "+str(root_index));
  if (l!=expression_base::no_value)
     push_value(std::make_shared<vector_value>(rd->val.root(alpha)));
}
void coroot_wrapper(expression_base::level l)
{ int root_index = get<int_value>()->int_val();
  shared_root_datum rd(get<root_datum_value>());
  RootNbr npr = rd->val.numPosRoots();
  RootNbr alpha = npr+root_index;
  if (alpha>=2*npr)
    throw runtime_error("Illegal coroot index "+str(root_index));
  if (l!=expression_base::no_value)
     push_value(std::make_shared<vector_value>(rd->val.coroot(alpha)));
}

@ We also allow access to the matrices of all simple or of all positive
(co)roots. For all roots such access is rarely needed, and if so easily
programmed, so we don't provide a built-in function for that.

@< Local function definitions @>=
void simple_roots_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
@)
  WeightList srl
    (rd->val.beginSimpleRoot(),rd->val.endSimpleRoot());
  push_value(std::make_shared<matrix_value>(int_Matrix(srl,rd->val.rank())));
}
@)
void simple_coroots_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
@)
  WeightList scl
    (rd->val.beginSimpleCoroot(),rd->val.endSimpleCoroot());
  push_value(std::make_shared<matrix_value>(int_Matrix(scl,rd->val.rank())));
}
@)
void positive_roots_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  WeightList rl
    (rd->val.beginPosRoot(),rd->val.endPosRoot());
  push_value(std::make_shared<matrix_value>
    (int_Matrix(rl,rd->val.rank())));
}
@)
void positive_coroots_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  WeightList crl
    (rd->val.beginPosCoroot(),rd->val.endPosCoroot());
  push_value(std::make_shared<matrix_value>
    (int_Matrix(crl,rd->val.rank())));
}

@ Here are some important attributes of root data, and look-up
functions for roots and coroots.

@< Local function definitions @>=
void root_datum_eq_wrapper (expression_base::level l)
{ shared_root_datum rd1 = get<root_datum_value>();
  shared_root_datum rd0 = get<root_datum_value>();
  if (l!=expression_base::no_value)
    push_value(whether(rd0.get()==rd1.get())); // compare pointers
}
void root_datum_neq_wrapper (expression_base::level l)
{ shared_root_datum rd1 = get<root_datum_value>();
  shared_root_datum rd0 = get<root_datum_value>();
  if (l!=expression_base::no_value)
    push_value(whether(rd0.get()!=rd1.get())); // compare pointers
}
@)
void datum_Cartan_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  int_Matrix M = rd->val.cartanMatrix();
  push_value(std::make_shared<matrix_value>(M));
}
@)
void rd_rank_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(rd->val.rank()));
}
@)
void rd_semisimple_rank_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(rd->val.semisimpleRank()));
}
@)
void rd_nposroots_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(rd->val.numPosRoots()));
}
@)
void root_index_wrapper(expression_base::level l)
{ shared_vector alpha = get<vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  int index=rd->val.root_index(alpha->val); // ensure signed type here
  index -= static_cast<int>(rd->val.numPosRoots()); // and signed subtract here
  push_value(std::make_shared<int_value>(index));
}
@)
void coroot_index_wrapper(expression_base::level l)
{ shared_vector alpha_v = get<vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  int index=rd->val.coroot_index(alpha_v->val); // ensure signed type here
  index -= static_cast<int>(rd->val.numPosRoots()); // and signed subtract here
  push_value(std::make_shared<int_value>(index));
}


@ It is useful to have bases for the sum of the root lattice and the
coradical, and for the sum of the coroot lattice and the radical; the latter
will in fact be used later in a function to built inner classes.

@< Local function definitions @>=
void root_coradical_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  WeightList srl
    (rd->val.beginSimpleRoot(),rd->val.endSimpleRoot());
  srl.insert(srl.end(),rd->val.beginCoradical(),rd->val.endCoradical());
  push_value(std::make_shared<matrix_value>(int_Matrix(srl,rd->val.rank())));
}
@)
void coroot_radical_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  WeightList scl
    (rd->val.beginSimpleCoroot(),rd->val.endSimpleCoroot());
  scl.insert(scl.end(),rd->val.beginRadical(),rd->val.endRadical());
  push_value(std::make_shared<matrix_value>(int_Matrix(scl,rd->val.rank())));
}

@ We give access to the fundamental weights and coweights on an individual
basis, which is easier since they are rational vectors.

@< Local function definitions @>=
void fundamental_weight_wrapper(expression_base::level l)
{ int i= get<int_value>()->int_val();
  shared_root_datum rd = get<root_datum_value>();
  if (unsigned(i)>=rd->val.semisimpleRank())
    throw runtime_error("Invalid index "+str(i));
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value> @|
      (rd->val.fundamental_weight(i).normalize()));
}
@)
void fundamental_coweight_wrapper(expression_base::level l)
{ int i= get<int_value>()->int_val();
  shared_root_datum rd = get<root_datum_value>();
  if (unsigned(i)>=rd->val.semisimpleRank())
    throw runtime_error("Invalid index "+str(i));
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value> @|
      (rd->val.fundamental_coweight(i).normalize()));
}


@ And here are functions for the dual, derived and quotient by central torus
root data. There is no need for a similar function for the adjoint root datum,
as this is easily synthesised (due to the fact that the simple roots provide a
standard basis for the adjoint character lattice): the simple root matrix is
the identity, the simple coroot matrix the Cartan matrix, and mapping to new
weight coordinates is achieved by pairing with the old coroots.

@h "tags.h"
@< Local function definitions @>=
void dual_datum_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l!=expression_base::no_value)
    push_value(rd->dual());
}
@)
void derived_info_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  int_Matrix projector;
  PreRootDatum pre(projector,rd->val,tags::DerivedTag());
  push_value(root_datum_value::build(std::move(pre)));
  push_value(std::make_shared<matrix_value>(projector));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}
@)
void mod_central_torus_info_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  int_Matrix injector;
  PreRootDatum pre(injector,rd->val,tags::CoderivedTag());
  push_value(root_datum_value::build(std::move(pre)));
  push_value(std::make_shared<matrix_value>(injector));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ This function is more recent; it allows constructing a new root (sub-)datum
by selecting coroots taking integral values on a given rational weight vector.

@< Local function definitions @>=
void integrality_datum_wrapper(expression_base::level l)
{ shared_rational_vector lambda = get<rational_vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (lambda->val.size()!=rd->val.rank())
    throw runtime_error()
      << "Length " << lambda->val.size() @|
      << " of rational vector differs from rank " << rd->val.rank();
@.Length of rational vector...@>
  if (l!=expression_base::no_value)
  @/push_value(root_datum_value::build @|
      (rootdata::integrality_predatum(rd->val,lambda->val)));
}

@ A related function computes a list of fractions of a line segment where the
set of roots with integrality is non-empty.

@< Local function definitions @>=
void integrality_points_wrapper(expression_base::level l)
{ shared_rational_vector lambda = get<rational_vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (lambda->val.size()!=rd->val.rank())
    throw runtime_error()
      << "Length " << lambda->val.size() @|
      << " of rational vector differs from rank " << rd->val.rank();
@.Length of rational vector...@>
  if (l==expression_base::no_value)
    return;
@)
  RationalList ipl = rootdata::integrality_points(rd->val,lambda->val);
    // method normalises rationals
  own_row result = std::make_shared<row_value>(ipl.size());
  for (size_t i=0; i<ipl.size(); ++i)
    result->val[i]=std::make_shared<rat_value>(ipl[i]);
  push_value(std::move(result));
}

@ Let us install the above wrapper functions.

@< Install wrapper functions @>=
install_function(type_of_root_datum_wrapper,@|"Lie_type"
                ,"(RootDatum->LieType)");
install_function(coroot_preference_wrapper,@|"prefers_coroots"
                ,"(RootDatum->bool)");
install_function(root_datum_wrapper,@|"root_datum"
                ,"(mat,mat,bool->RootDatum)");
install_function(root_datum_from_type_wrapper,@|"root_datum"
		,"(LieType,mat,bool->RootDatum)");
install_function(sublattice_root_datum_wrapper,@|"root_datum"
                ,"(RootDatum,mat->RootDatum)");
install_function(simply_connected_datum_wrapper,@|"simply_connected"
                ,"(LieType,bool->RootDatum)");
install_function(adjoint_datum_wrapper,@| "adjoint"
                ,"(LieType,bool->RootDatum)");
install_function(root_wrapper,@|"root","(RootDatum,int->vec)");
install_function(coroot_wrapper,@|"coroot","(RootDatum,int->vec)");
install_function(simple_roots_wrapper,@|"simple_roots","(RootDatum->mat)");
install_function(simple_coroots_wrapper,@|"simple_coroots","(RootDatum->mat)");
install_function(positive_roots_wrapper,@| "posroots","(RootDatum->mat)");
install_function(positive_coroots_wrapper,@| "poscoroots","(RootDatum->mat)");
install_function(root_datum_eq_wrapper,@|"=","(RootDatum,RootDatum->bool)");
install_function(root_datum_neq_wrapper,@|"!=","(RootDatum,RootDatum->bool)");
install_function(datum_Cartan_wrapper,@|"Cartan_matrix","(RootDatum->mat)");
install_function(root_coradical_wrapper,@|"root_coradical","(RootDatum->mat)");
install_function(coroot_radical_wrapper,@|"coroot_radical","(RootDatum->mat)");
install_function(fundamental_weight_wrapper,@|
		 "fundamental_weight","(RootDatum,int->ratvec)");
install_function(fundamental_coweight_wrapper,@|
		 "fundamental_coweight","(RootDatum,int->ratvec)");
install_function(dual_datum_wrapper,@|"dual","(RootDatum->RootDatum)");
install_function(derived_info_wrapper,@|
		 "derived_info","(RootDatum->RootDatum,mat)");
install_function(mod_central_torus_info_wrapper,@|
		 "mod_central_torus_info","(RootDatum->RootDatum,mat)");
install_function(rd_rank_wrapper,@|"rank","(RootDatum->int)");
install_function(rd_semisimple_rank_wrapper@|
		,"semisimple_rank","(RootDatum->int)");
install_function(rd_nposroots_wrapper@|,"nr_of_posroots","(RootDatum->int)");
install_function(root_index_wrapper@|,"root_index","(RootDatum,vec->int)");
install_function(coroot_index_wrapper@|,"coroot_index","(RootDatum,vec->int)");
install_function(integrality_datum_wrapper
                ,@|"integrality_datum","(RootDatum,ratvec->RootDatum)");
install_function(integrality_points_wrapper
                ,@|"integrality_points","(RootDatum,ratvec->[rat])");

@*1 A type for complex reductive groups equipped with an involution.
%
We shall now go ahead to define a primitive type holding an |InnerClass|
object, which represents a complex reductive group equipped with a
distinguished involution defining an ``inner class'' of real forms for that
group. We can construct such an object from a based root datum and an
involution of it. In the \.{Fokko} program, such involutions are entered
indirectly via a user interaction, and the way it was constructed is stored
and used for certain purposes. For maximal flexibility, we want to be able to
provide an involution produced in any way we like. This means there is some
extra work to do.

@< Includes... @>=
#include "innerclass.h"

@*2 Analysing involutions. Our constructor for the current atlas type must
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
exchanges a pair of basis vectors (we call these Complex factors). The number of
each of these subspaces is uniquely determined by the involution, and computed
by |tori::classify|.

@h "tori.h"
@f delta nullptr

@< Local function def...@>=
void classify_wrapper(expression_base::level l)
{ shared_matrix M(get<matrix_value>());
  const auto& delta=M->val;
  const auto r = delta.numRows();
  @< Check that |delta| is an $r\times{r}$ matrix defining an involution @>
  if (l==expression_base::no_value)
    return;
@)
  const auto ranks = tori::classify(delta);
  push_value(std::make_shared<int_value>(std::get<0>(ranks))); // compact rank
  push_value(std::make_shared<int_value>(std::get<1>(ranks))); // C rank
  push_value(std::make_shared<int_value>(std::get<2>(ranks))); // split rank
  if (l==expression_base::single_value)
    wrap_tuple<3>();
}

@ Since |classify_wrapper| can be called with an arbitrary matrix, is certainly
necessary to test that |delta| is an involution ($\delta^2-I=0$) to protect
|tori::classify| against abuse. The code below will also be called in the
context of checking the involution for an inner class, where it may seem
less essential, given the fact that we shall also check that |delta| induces an
involutive permutation of the roots. However even there it is not redundant in
the presence of a central torus part.

@< Check that |delta| is an $r\times{r}$ matrix defining an involution @>=
{ if (delta.numRows()!=r or delta.numColumns()!=r)
     throw runtime_error()
       << "Involution should be a " << r << 'x' << r << " matrix;"@|
@.Involution should be...@>
        " received a "  << delta.numRows() << 'x' << delta.numColumns()
      << " matrix";
  if (not (delta*delta-1).is_zero()) throw runtime_error
      ("Given transformation is not an involution");
@.Given transformation...@>
}

@ We now come to the part of the analysis that involves the root datum. Since
the matrix is already expressed on the basis of the weight lattice used by the
root datum, the question of stabilising that lattice is settled, but we must
check that the matrix is indeed an involution (for which we reuse the above
module), and that it gives an automorphism of the root datum. In fact we need an
automorphism of the \emph{based} root datum, but we will allow any root datum
automorphism to be supplied, and will compute from it a conjugate that sends all
positive roots to positive roots. The following auxiliary function checks
whether |delta|, is an involution that maps simple roots to roots and simple
coroots to coroots (throwing an error if it not), and returns a list of images
of simple roots that can be used by the caller to transform |delta| into a based
root datum involution.

That |delta| is a root datum automorphism means that its action from the left
permutes roots among each other, that its action on the right permutes coroots
among each other (this is not implied by the first condition). If, as we do in
the code below, we only test that left multiplication by~$\delta$
maps \emph{simple} roots to certain roots, and right multiplication
by~$\delta^{-1}$ maps \emph{simple} coroots to certain coroots, then the Cartan
matrix defined for the images will automatically be the same as for the original
simple system. Thus the image of the simple system is an isomorphic simple
system contained within the root datum itself, and this can only be the case if
it is actually obtained from a root datum automorphism; therefore our test is
sufficient.

@f Delta nullptr

@< Local function def...@>=
RootNbrList check_root_datum_involution
  (const RootDatum& rd, const WeightInvolution& delta)
{ const size_t r=rd.rank(), s=rd.semisimpleRank();
  @< Check that |delta| is an $r\times{r}$ matrix defining an involution @>
  RootNbrList Delta(s);
  for (weyl::Generator i=0; i<s; ++i)
  { Delta[i]=rd.root_index(delta*rd.simpleRoot(i));
    if (Delta[i]==rd.numRoots()) // then image not found
      throw runtime_error@|
        ("Matrix maps simple root "+str(i)+" to non-root");
    if (delta.right_prod(rd.simpleCoroot(i))!=rd.coroot(Delta[i]))
      throw runtime_error@|
        ("Matrix does not map simple coroot "+str(i)
        @|+" to coroot "+str(Delta[i]-rd.numPosRoots()));
  }
  return Delta;
}
@)
void check_based_root_datum_involution
  (const RootDatum& rd, const WeightInvolution& delta)
{ const auto Delta=check_root_datum_involution(rd,delta);
  const auto s=rd.semisimpleRank();
  for (weyl::Generator i=0; i<s; ++i)
    if (not rd.is_simple_root(Delta[i]))
      throw runtime_error ("Root datum involution is not distinguished");
}

@ At a more outer level, the function |check_involution| will do all pertinent
checks, and when successful both modify its argument |delta| to a based root
datum involution, and export through the |ww| parameter the Weyl group element
whose conjugation was used; finally the return value is the |Twist| of the
Dynkin diagram corresponding to the modified~|delta| (this allows also using
this function easily to test whether an involution is proper for a given inner
class). After passing the test |check_root_datum_involution|, we are sure that
|wrt_distinguished| will be able to map the images |Delta| to the simple roots
by a Weyl group element.

In addition to this, the function can also determine the Lie type of
the root datum, the permutation possibly needed to map the standard (Bourbaki)
ordering of that Dynkin diagram to the actual one, and the inner class letters
corresponding to each of its factors. This is precisely the data stored in a
|lietype::Layout| structure, so we provide as final argument a pointer |lo| to
such a structure; if non null, those values will be computed and stored there.
Given the amount of effort required below, one might wonder if recording a
|Layout| in an inner class is really necessary, but its presence is essential
to be able to associate real form names to internal data (namely gradings), so
in the current set-up its computation cannot be avoided.

The Lie type will in general be identical to that of the root datum (and in
particular any torus factors will come at the end), but in case of Complex
inner classes we may be forced to permute the simple factors to make the
identical factors associated to such classes adjacent (the permutation will be
adapted to reflect this reordering of simple factors).

We shall not apply the function |classify_involution| defined above to the
entire matrix describing a (purported) involution of the weight lattice,
but rather we shall apply it below to a matrix defined for the central torus
part.

@h "weyl.h"

@< Local function def...@>=
weyl::Twist check_involution
 (WeightInvolution& delta, const RootDatum& rd,
  WeylWord& ww, @| lietype::Layout* lo=nullptr)
{ RootNbrList Delta = check_root_datum_involution(rd,delta);
  ww = wrt_distinguished(rd,Delta);
  const size_t r=rd.rank(), s=rd.semisimpleRank();
  weyl::Twist p; // result
@/ @< Copy the permutation of the simple roots in |Delta| to |p|, and
      left-act on |delta| by the reverse of~|ww| to make it match |Delta| @>
  if (lo==nullptr)
    return p; // if no details are asked for, we are done now
@/LieType& type=lo->d_type;
  InnerClassType& inner_class=lo->d_inner;
  Permutation& pi=lo->d_perm;
  @< Compute the Lie type |type|, the inner class |inner_class|, and the
     permutation |pi| of the simple roots with respect to standard order for
     |type| @>
  if (r>s)
    @< Add type letters to |type| and inner class symbols to |inner_class|
       for the central torus @>
  return p;
}

@ The call to |wrt_distinguished| has left in |Delta| the permuted simple roots,
which information we transfer to~|pi|; then we transform the matrix~|delta| to
match the transformation that was applied to~|Delta|, so that |delta| becomes a
distinguished involution.

@< Copy the permutation...@>=
{ for (weyl::Generator i=0; i<s; ++i)
    p[i]=rd.simpleRootIndex(Delta[i]);
    // this |assert|s that |Delta[i]| is simple
@)
// now adapt |delta| so that it becomes the inner class distinguished involution
  for (unsigned int i=0; i<ww.size(); ++i) // apply elements in generation order
    rd.simple_reflect(ww[i],delta);
}

@ For each simple factor we look if there are any non-fixed points of the
permutation (if not we have the compact inner class) and if so, whether the
image of that point lies in the same component of the Dynkin diagram. If the
latter is the case we have an unequal rank inner class, which is actually
called ``split'' unless the simple factor is of type $D_{2n}$, and if the image
lies in another component of the diagram we have a Complex inner class.

@h "dynkin.h"

@< Compute the Lie type |type|, the inner class... @>=
{ DynkinDiagram diagram(rd.cartanMatrix());
  containers::sl_list<RankFlags> comps =
    diagram.components(); // connected components
  type = diagram.classify_semisimple(pi,true);
    // |pi| normalises to Bourbaki order
  assert(type.size()==comps.size());
@)
  inner_class.reserve(comps.size()+r-s); // certainly enough
  size_t offset=0; // accumulated rank of simple factors seen

  unsigned int i=0; // index into |type|
  for (auto cit=comps.begin(); not comps.at_end(cit); ++cit,++i)
  { bool equal_rank=true;
    size_t comp_rank = cit->count();
    assert (comp_rank==type[i].rank());
       // and |*it| runs through bits in |*cit| in following loop
    for (auto it=&pi[offset]; it!=&pi[offset+comp_rank]; ++it)
      if (p[*it]!=*it) {@; equal_rank=false; break; }
@)  if (equal_rank) inner_class.push_back('c');
      // identity on this component: compact component
    else if (cit->test(p[pi[offset]]))
       // (any) one root stays in component |*cit|: unequal rank
      inner_class.push_back(type[i].first=='D' and comp_rank%2==0 ? 'u' : 's');
    else
    { inner_class.push_back('C'); // record Complex component
      @< Gather elements of Complex inner class component, adapting the values
         of |type|, |comps|, and~|pi| @>
      offset += comp_rank;
      ++cit,++i; // skip over component |i|, loop will skip component |i+1|
    }

    offset += comp_rank;
  } // |for (cit)|
}

@ Complex factors of the inner class involve two simple factors, which requires
some additional care. The corresponding components of the Dynkin diagram might
not be consecutive, in which case we must permute the factors of |type| to make
that true, permute the |comp| subsets to match, and update |pi| as well.
Moreover we wish that, as seen through |pi|, the permutation |p| interchanges
the roots of the two now consecutive factors by a fixed shift in the index by
|comp_rank| (the kind of permutation that |lietype::involution| produces).

Due to this predetermined order (relative to |p|) in which the matching factor
is to be arranged, it is not necessary to record the original values of |pi|
on this factor (which were produced by |dynkin::Lie_type| in increasing order)
nor the type of the factor, before overwriting them by other factors being
shifted up: we can afterwards simply deduce the proper values from the
matching factor |i|.

@< Gather elements of Complex inner class...@>=
{ auto beta = p[pi[offset]];
    // index of simple root, |delta| image of first one in current component
  @< Find the component~|k| after |i| that contains |beta|, rotate entry |k| of
     |comps| to position |i+1|, while shifting values in |type| and |pi|
     upwards by |1| respectively |comp_rank| places @>
@)
  type[i+1]=type[i]; // duplicate factor |i|
  for (size_t j=offset; j<offset+comp_rank; ++j)
    pi[j+comp_rank]=p[pi[j]];
     // reconstruct matching component in matching order, applying~|p|
}


@ When the inner class permutation |p| interchanges a component with another, we
use the image~|beta| of one root in the former component to search for the
latter component. Then, as remarked above, we can simply shift up any
intermediate values of |pi|, |type| and |comp| to their new places; only for the
bitset |comp[k]| it is worth while to save the old value and reinsert it at its
moved-down place.

@< Find the component~|k| after |i| that contains |beta|...@>=
{ auto k=i; auto cit1=cit; // both are to be immediately incremented
  while (++k, not comps.at_end(++cit1))
    if (cit1->test(beta))
      break;
  if (comps.at_end(cit1))
    throw logic_error("Non matching Complex factor");
@.Non matching Complex factor@>

#ifndef NDEBUG
  assert(type[k]==type[i]); // paired simple types for complex factor
  for (unsigned int l=1; l<comp_rank; ++l)
    assert(cit1->test(p[pi[offset+l]]));
        // image by |p| of remainder of |comp[i]| matches |comp[k]|
#endif

  if (k>i+1) // then we need to move component |k| down to |i+1|
  {
    comps.splice(std::next(cit),comps,cit1);
      // rotate node |*cit1| to position after |cit|
    std::copy_backward(&type[i+1],&type[k],&type[k+1]); // shift up |1|
@)
    auto j=offset+comp_rank;
    for (auto it=&type[i+1]; it!=&type[k]; ++it)
      j += it->rank();
    std::copy_backward(&pi[offset+comp_rank],&pi[j],&pi[j+comp_rank]);
    // shift up |comp_rank|
  }

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

@< Add type letters to |type| and inner class symbols to |inner_class|
   for the central torus @>=
{ for (size_t k=0; k<r-s; ++k)
    type.push_back(SimpleLieType('T',1));
  int_Matrix root_lattice
    (rd.beginSimpleRoot(),rd.endSimpleRoot(),r,tags::IteratorTag());
@/CoeffList factor; // values will be unused
  int_Matrix basis =
     matreduc::adapted_basis(root_lattice,factor);
@/WeightInvolution inv =
     basis.inverse().block(s,0,r,r)*delta*basis.block(0,s,r,r);
    // involution on quotient by root lattice
  const auto ranks=tori::classify(inv);
@/auto compact_rank=std::get<0>(ranks);
  auto Complex_rank=std::get<1>(ranks);
  auto split_rank=std::get<2>(ranks);
@/while (compact_rank-->0) inner_class.push_back('c');
  while (Complex_rank-->0) inner_class.push_back('C');
  while (split_rank-->0) inner_class.push_back('s');
}

@*2 Storing the inner class values.
%
Although abstractly an inner class value is described completely by an object of
type |InnerClass|, we shall need to record additional information in order to be
able to present meaningful names for the real forms and dual real forms in this
inner class. The above analysis of involutions was necessary in order to obtain
such information; the class |output::FormNumberMap| will serve to record the
necessary renumbering and naming of real form numbers.

@< Includes... @>=
#include "output.h"

@~The class |inner_class_value| will be the first Atlas type where we deviate
from the previously used scheme of holding an Atlas library object with the main
value in a data member |val|. The reason is that the copy constructor for
|InnerClass| is deleted, so that the straightforward definition of a copy
constructor for such an Atlas type would not work, and the copy constructor is
necessary for the |clone| method. (In fact, now that normal manipulation of
values involves duplicating shared pointers rather than of values, there is
never a need to copy an |inner_class_value|, since |get_own<inner_class_value>|
is never called; however the |clone| method is still defined for possible future
use.) So instead, we shall share the library object when duplicating our value,
and maintain a reference count to allow destruction when the last copy
disappears.

The reference count needs to be shared of course, and since the links between
the |inner_class_value| and both the library value and the reference count
are indissoluble, we use references for the members |val|, |dual| and
|ref_count|. The first two references are not |const|, since some methods will
as a side effect generate |CartanClass| objects in the inner class, whence
they are technically manipulators rather than accessors.

The main constructor takes a unique-pointer to an |InnerClass| as
argument, as a reminder that the caller gives up ownership of this pointer
that should come from a call to~|new|; this pointer will henceforth be owned
by the |inner_class_value| constructed, in shared ownership with any values
later cloned from it: the last one of them to be destroyed will call |delete|
for the pointer. The remaining argument is a |Layout| that must have been
computed by |check_involution| above, in order to ensure its validity.

Occasionally we shall need to refer to the dual inner class (for the dual
group); since the construction of an instance takes some work, we do not wish
to repeat that every time the dual is needed, so we create the dual
|InnerClass| value upon construction of the
|inner_class_value| and store it in the |dual| field where it will be
available as needed.

Unlike for other value types, the copy constructor is public here. Thus a
class depending on our |val| being valid can simply store a copy of our
|inner_class_value|; this does not cost much, and the reference counting
mechanism will then ensure that |val| remains valid while the object of that
containing class exists.

@< Type definitions @>=
struct inner_class_value : public value_base
{ InnerClass& val;
  InnerClass& dual;
  size_t& ref_count;
  shared_root_datum datum, dual_datum; // share if we depend on them
@)
  lietype::LieType rd_type;
  lietype::InnerClassType ic_type;
  const output::FormNumberMap interface,dual_interface;
@)
  inner_class_value(std::unique_ptr<InnerClass> G, const lietype::Layout& lo);
  // main
static inner_class_value build
  (shared_root_datum srd, WeightInvolution& tau, WeylWord* wp=nullptr);
  // to share
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
typedef std::shared_ptr<const inner_class_value> shared_inner_class;

@ Here are the copy constructor and the destructor.
@< Function def...@>=
inner_class_value::inner_class_value(const inner_class_value& v)
: val(v.val), dual(v.dual), ref_count(v.ref_count)
@/, datum(v.datum), dual_datum(v.dual_datum) // possibly double up sharing
@/, rd_type(v.rd_type), ic_type(v.ic_type)
, interface(v.interface), dual_interface(v.dual_interface)
{@; ++ref_count; }

inner_class_value::~inner_class_value()
{@; if (--ref_count==0) {@; delete &val; delete &dual; delete &ref_count;} }

@ The main constructor installs the reference to the Atlas value, creates its
dual value to which a reference is stored, and allocates and initialises the
reference count. The |Layout| argument is stored for convenience, and is used
to initialise the |interface| and |dual_interface| fields.

The current constructor has the defect that the reference to which the |dual|
field is initialised will not be freed if an exception should be thrown during
the subsequent initialisations (which is not likely, but not impossible
either). However, we see no means to correct this defect with the given
structure, since a member of reference type has to be initialised to its
definitive value, and no local variables are available during initialisation
to which we could give temporary ownership of that reference. Even if we
wrapped a function-try-block around the entire constructor, it would not be
allowed to access |dual| to free it, nor could it know whether |dual| had been
successfully initialised. In fact the same problem exists for |ref_count| (if
construction of either of the interfaces should throw, the reference count
variable itself will remain dangling), which means the problem cannot be
solved either by reordering the data members (and therefore their
initialisations). For |val| the problem does not arise because it is owned by
a smart pointer until the construction succeeds.

@< Function def...@>=
inner_class_value::inner_class_value
  (std::unique_ptr<InnerClass> g,
   const lietype::Layout& lo)
@/: val(*g)
, dual(*new InnerClass(*g,tags::DualTag()))
@/, ref_count(*new size_t(1))
@/, datum(nullptr), dual_datum(nullptr)
    // we use independent |InnerClass| objects
@/, rd_type(lo.d_type), ic_type(lo.d_inner)
, interface(*g,lo), dual_interface(*g,lo,tags::DualTag())
 {@; g.release(); } // now that we own |g|, release the unique-pointer

@ We allow construction of a dual |inner_class_value|. Since it can share the
two fields referring to objects of type |InnerClass|
in the opposite order, we can consider it as a member of the same reference
counted family, and share the |ref_count| field. This means this constructor
is more like the copy constructor than like the main constructor, and in
particular the reference count is increased. The dual of the dual inner class
value will behave exactly like a copy of the original inner class.

@< Function def...@>=
inner_class_value::inner_class_value(const inner_class_value& v,tags::DualTag)
@/: val(v.dual), dual(v.val)
, ref_count(v.ref_count)
@/, datum(v.dual_datum), dual_datum(v.datum) // possibly double up sharing
@/, rd_type(lietype::dual_type(v.rd_type))
, ic_type(lietype::dual_type(v.ic_type,rd_type))
, interface(v.dual_interface), dual_interface(v.interface)
{@; ++ref_count; }

@ We often want an |inner_class_value| to use an |InnerClass| that refers to
rather than stores a root datum and its dual, if those can be found in an
existing |root_datum_value| and its dual. To that end we provide a
pseudo-constructor |build| that will build such an |InnerClass| object, and
store shared pointed to the root data it depends on in the associated
|inner_class_value|.

@< Function def...@>=
inner_class_value inner_class_value::build
  (shared_root_datum srd, WeightInvolution& tau, WeylWord* wp)
{ shared_root_datum drd=srd->dual();
  const auto& rd=srd->val;
  WeylWord ww;
  lietype::Layout lo;
  if (wp==nullptr)
    wp=&ww; // use |ww| as dummy output unless |p| points somewhere
  check_involution(tau,rd,*wp,&lo); // may also modify |tau|, and sets |lo|
  std::unique_ptr<InnerClass> p @| (new InnerClass(rd,drd->val,tau));
    // depends on |srd| and |drd|
  inner_class_value result(std::move(p),lo); // use main constructor here
  result.datum=std::move(srd);
  result.dual_datum=std::move(drd); // set dependencies
  return result; // return modified instance of |inner_class_value|
}

@ One of the most practical informations about an |InnerClass|,
which is available directly after its construction, is the number of real
forms in the inner class defined by it; we print this information when a
|inner_class_value| is printed.

@< Function def...@>=
void inner_class_value::print(std::ostream& out) const
{ out << "Complex reductive group of type " << rd_type @|
      << ", with involution defining\n"
         "inner class of type '" << ic_type @|
      << "', with " << val.numRealForms() @| << " real "
      << (val.numRealForms()==1 ? "form" : "forms") @| << " and "
      << val.numDualRealForms() @| << " dual real "
      << (val.numDualRealForms()==1 ? "form" : "forms");
}

@ Here is the wrapper function for constructing an inner class using an
involution, testing its validity in the process. The Weyl word |ww| that was set
by |check_involution| is not used, as it suffices that the call has modified~|M|
to be a distinguished involution; similarly the |weyl::Twist| value returned by
by the call is ignored here. Then the root datum and matrix are passed to a
|InnerClass| constructor that the library provides specifically for this
purpose, and which makes a copy of the root datum; the \.{Fokko} program instead
uses a constructor using a |PreRootDatum| that constructs the |RootDatum|
directly into the |InnerClass|. Using that constructor here would be cumbersome
and even less efficient then copying the existing root datum.

@< Local function def...@>=
void fix_involution_wrapper(expression_base::level l)
{ LatticeMatrix M(get<matrix_value>()->val);
    // copy construct, safely using temporary
  shared_root_datum rd = get<root_datum_value>();
  auto ic_value = inner_class_value::build(rd,M); // build and check
  if (l!=expression_base::no_value)
    push_value(std::make_shared<inner_class_value>(ic_value));
    // light weight copy
}

@ Another wrapper |twisted_involution_wrapper| is similar, but also returns a
second value, which is the twisted involution corresponding in the inner class
to the given involution matrix.

@< Local function def...@>=
void twisted_involution_wrapper(expression_base::level l)
{ LatticeMatrix M(get<matrix_value>()->val);
  shared_root_datum rd = get<root_datum_value>();
  WeylWord ww;
  auto ic_value = inner_class_value::build(rd,M,&ww); // build and check
  if (l==expression_base::no_value)
    return;
@)
 push_value(std::make_shared<inner_class_value>(ic_value));
  // light weight copy |ic_value|
  push_value(std::make_shared<vector_value>
    (std::vector<int>(ww.begin(),ww.end())));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ To simulate the functioning of the \.{Fokko} program, the following overload
of the function |inner_class| in \.{atlas} takes as argument a Lie type, a list
of kernel generators, and a string describing the inner class. The evaluation of
the \.{atlas} call |inner_class(lt,gen,ict)| will successively compute ${\it
basis}={\it quotient\_basis(lt,gen)}$, ${\it rd}={\it root\_datum (lt,basis)}$,
and ${\it M}={\it based\_involution(lt,basis,ict)}$, and then returns the same
value as would |fix_involution(rd,M)|. However we avoid actually calling the
function |fix_involution_wrapper|, which would reconstruct |lt| and |ict| from
|rd| and |M|, while performing tests that are useless given the way $M$ was
computed; rather we store |lt| and |ict| directly in a |Layout| to be stored in
the |inner_class_value|. This follows most closely \.{Fokko} behaviour, and
avoids surprises (however inner class letters do change to synonyms as they
usual do when passing through |checked_inner_class_type|).

@< Local function def...@>=
void inner_class_from_type_wrapper(expression_base::level l)
{ bool prefer_coroots = false;
  shared_string ict = get<string_value>();
    // and leave generators |gen| and type |lt|
  shared_value lt = *(execution_stack.end()-2);
  const LieType& type=force<Lie_type_value>(lt.get())->val;
  lietype::Layout lo(type,checked_inner_class_type(ict->val.c_str(),type));
@)
  quotient_basis_wrapper(expression_base::single_value); @+
  shared_value basis = pop_value();
@)
  push_value(lt); push_value(basis);
  push_value(whether(prefer_coroots));
  root_datum_from_type_wrapper(expression_base::single_value);
  shared_root_datum rd = get<root_datum_value>();
  shared_root_datum dual_rd = rd->dual();
@)
  push_value(lt); push_value(basis);
  push_value(ict);
  based_involution_wrapper(expression_base::single_value);
  shared_matrix M = get<matrix_value>();
  if (l==expression_base::no_value)
    return; // bow out now all possible errors are passed
@)
  std::unique_ptr<InnerClass>@|
    G(new InnerClass(rd->val,dual_rd->val,M->val));
  auto result = std::make_shared<inner_class_value>(std::move(G),lo);
@/result->datum=std::move(rd);
  result->dual_datum=std::move(dual_rd); // set dependencies
  push_value(result);
}

@*2 Functions operating on inner classes.
%
Here are our first functions that access a operate on values of type
|inner_class_value|; they recover the ingredients that were used in the
construction, and the final one constructs the dual inner class.

@< Local function def...@>=
void inner_class_eq_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  shared_inner_class H = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(whether(&G->val==&H->val));
}
void inner_class_neq_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  shared_inner_class H = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(whether(&G->val!=&H->val));
}
@)
void distinguished_involution_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value>(G->val.distinguished()));
}
@)
void root_datum_of_inner_class_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(root_datum_value::build(G->val.rootDatum()));
}
@)
void inner_class_to_root_datum_coercion()
{@; root_datum_of_inner_class_wrapper(expression_base::single_value); }

void dual_inner_class_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<inner_class_value>(*G,tags::DualTag()));
}

@ More interestingly, let us extract the list of names of the real forms.
This uses the interface fields stored in the value. Since they exist for both
the group itself and for the dual group, we define an auxiliary function that
produces the list, and then use it twice.

@< Local function def...@>=
void push_name_list(const output::FormNumberMap& interface)
{ own_row result = std::make_shared<row_value>(0);
  for (size_t i=0; i<interface.numRealForms(); ++i)
    result->val.emplace_back
      (std::make_shared<string_value>(interface.type_name(i)));
  push_value(std::move(result));
}
@)
void form_names_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_name_list(G->interface);
}
@)
void dual_form_names_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_name_list(G->dual_interface);
}

@ We provide functions for counting (dual) real forms and Cartan classes,
although the former two could be computed as the lengths of the list of (dual)
form names.

@< Local function def...@>=
void n_real_forms_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(G->val.numRealForms()));
}
@)
void n_dual_real_forms_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(G->val.numDualRealForms()));
}
@)
void n_Cartan_classes_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(G->val.numCartanClasses()));
}

@ And now, our first function that really simulates something that can be done
with \.{Fokko} using more than just a root datum. This is the \.{blocksizes}
command from \.{mainmode.cpp}, which uses |innerclass_io::printBlockSizes|, but
we have to rewrite it to avoid calling an output routine. A subtle difference is
that we use a matrix of integers rather than of |arithmetic::big_int| values to
collect the block sizes; this avoids having to define a new primitive type, and
probably suffices for the cases where actual computation of the block feasible.
However, we also provide a function the computes a single block size as an
unbounded integer.

@< Local function def...@>=
void block_sizes_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l==expression_base::no_value)
    return;
@)
  own_matrix M = std::make_shared<matrix_value> @|
    (int_Matrix(G->val.numRealForms(),G->val.numDualRealForms()));
  for (size_t i = 0; i < M->val.numRows(); ++i)
    for (size_t j = 0; j < M->val.numColumns(); ++j)
      M->val(i,j) =
      G->val.block_size(G->interface.in(i),G->dual_interface.in(j)).int_val();
  push_value(std::move(M));
}

void block_size_wrapper(expression_base::level l)
{ auto j = get<int_value>()->val.ulong_val();
  auto i = get<int_value>()->val.ulong_val();
  shared_inner_class G = get<inner_class_value>();
  if (i>=G->val.numRealForms())
    throw runtime_error("Real form number "+str(i)+" out of bounds");
  if (i>=G->val.numRealForms())
    throw runtime_error("Real form number "+str(i)+" out of bounds");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value> @|
      (G->val.block_size(G->interface.in(i),G->dual_interface.in(j))));
}



@ Next we shall provide a function that displays the occurrence of Cartan
classes for various real forms. The rows will be indexed by real forms, and
the columns by Cartan classes (note the alliteration).

@< Local function def...@>=
void occurrence_matrix_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l==expression_base::no_value)
    return;
@)
  size_t nr=G->val.numRealForms();
  size_t nc=G->val.numCartanClasses();
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(nr,nc));
  for (size_t i=0; i<nr; ++i)
  { BitMap b=G->val.Cartan_set(G->interface.in(i));
    for (size_t j=0; j<nc; ++j)
      M->val(i,j)= b.isMember(j) ? 1 : 0;
  }
  push_value(std::move(M));
}

@ We do the same for dual real forms. Note that we had to introduce the method
|dualCartanSet| for |InnerClass| in order to be able
to write this function.

@< Local function def...@>=
void dual_occurrence_matrix_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l==expression_base::no_value)
    return;
@)
  size_t nr=G->val.numDualRealForms();
  size_t nc=G->val.numCartanClasses();
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(nr,nc));
  for (size_t i=0; i<nr; ++i)
  { BitMap b=G->val.dual_Cartan_set(G->dual_interface.in(i));
    for (size_t j=0; j<nc; ++j)
      M->val(i,j)= b.isMember(j) ? 1 : 0;
  }
  push_value(std::move(M));
}

@ Finally we install everything related to inner classes.
@< Install wrapper functions @>=

install_function(classify_wrapper,@|"classify_involution"
                ,"(mat->int,int,int)");
install_function(fix_involution_wrapper,@|"inner_class"
                ,"(RootDatum,mat->InnerClass)");
install_function(twisted_involution_wrapper,@|"twisted_involution"
                ,"(RootDatum,mat->InnerClass,vec)");
install_function(inner_class_from_type_wrapper,@|"inner_class"
                ,"(LieType,[ratvec],string->InnerClass)");
install_function(inner_class_eq_wrapper,@|"=","(InnerClass,InnerClass->bool)");
install_function(inner_class_neq_wrapper,@|"!=","(InnerClass,InnerClass->bool)");
install_function(distinguished_involution_wrapper,@|"distinguished_involution"
                ,"(InnerClass->mat)");
install_function(root_datum_of_inner_class_wrapper,@|"root_datum"
                ,"(InnerClass->RootDatum)");
install_function(dual_inner_class_wrapper,@|"dual"
                ,"(InnerClass->InnerClass)");
install_function(form_names_wrapper,@|"form_names"
                ,"(InnerClass->[string])");
install_function(dual_form_names_wrapper,@|"dual_form_names"
                ,"(InnerClass->[string])");
install_function(n_real_forms_wrapper,@|"nr_of_real_forms"
                ,"(InnerClass->int)");
install_function(n_dual_real_forms_wrapper,@|"nr_of_dual_real_forms"
                ,"(InnerClass->int)");
install_function(n_Cartan_classes_wrapper,@|"nr_of_Cartan_classes"
                ,"(InnerClass->int)");
install_function(block_sizes_wrapper,"block_sizes","(InnerClass->mat)");
install_function(block_size_wrapper,"block_size","(InnerClass,int,int->int)");
install_function(occurrence_matrix_wrapper,@|"occurrence_matrix"
                ,"(InnerClass->mat)");
install_function(dual_occurrence_matrix_wrapper,@|"dual_occurrence_matrix"
                ,"(InnerClass->mat)");

@*1 A type for real reductive groups.
A next step in specifying the computational context is choosing a real form in
the inner class of them that was determined by a root datum and an involution.
This determines a, not necessarily connected, real reductive group inside the
connected complex reductive group; the corresponding Atlas class is called
|RealReductiveGroup|.

@< Includes... @>=
#include "realredgp.h"

@*2 Class definition.
The layout of this type of value is different from what we have seen before.
An Atlas object of class |RealReductiveGroup| is dependent upon another Atlas
object to which it stores a pointer, which is of type |InnerClass|,
so we must make sure that the object pointed to cannot disappear before it
does. The easiest way to do this is to place an |inner_class_value| object
|parent| inside the |real_form_value| class that we shall now define; the
reference-counting scheme introduced above then guarantees that the data we
depend upon will remain in existence sufficiently long. Since that data can be
accessed from inside the |RealReductiveGroup|, we shall mostly mention the
|parent| with the purpose of accessing its |interface| and |dual_interface|
fields. To remind us that the |parent| is not there to be changed by us, we
declare it |const|. The object referred to may in fact undergo internal change
however, via manipulators of the |val| field.

This class also serves to store persistent data related to the real form, in
values of type |KhatContext| and |Rep_table|. In order to avoid overhead at
construction, and also to not require including other header files from our
current header file, we store pointer that will only be assigned on first use.

@< Type definitions @>=
struct real_form_value : public value_base
{ const inner_class_value parent;
  RealReductiveGroup val;
@)
  real_form_value(const inner_class_value& p,RealFormNbr f) @/
  : parent(p), val(p.val,f)
  , khc_p(nullptr)
  , rt_p(nullptr) @+{}
  real_form_value
    (const inner_class_value& p,RealFormNbr f
    ,const RatCoweight& coch, TorusPart tp) @/
  : parent(p), val(p.val,f,coch,tp)
  , khc_p(nullptr)
  , rt_p(nullptr) @+{}
  virtual ~real_form_value ();
@)
  virtual void print(std::ostream& out) const;
  real_form_value* clone() const @+
    {@; return new real_form_value(*this); }
  static const char* name() @+{@; return "real form"; }
  const KGB& kgb () @+{@; return val.kgb(); }
   // generate and return $K\backslash G/B$ set
  KhatContext& khc();
  const Rep_context& rc();
  Rep_table& rt();
private:
  KhatContext* khc_p;
  Rep_table* rt_p;
    // owned pointers, initially |nullptr|, assigned at most once
};
@)
typedef std::shared_ptr<const real_form_value> shared_real_form;
typedef std::shared_ptr<real_form_value> own_real_form;

@ The method |khc| ensures a |KhatContext| value is constructed at |*khc_p|,
and similarly |rc| and |rt| ensure a |Rep_table| value is constructed at
|*rt_p|, and then these methods return an appropriate reference. The value so
obtained will serve to manipulate parameters for standard modules, for which
we shall define an Atlas type below. Storing the value here ensures that it
will be shared between different parameters for the same real form, and that
it will live as long as those parameter values do.

@< Function def...@>=
  KhatContext& real_form_value::khc()
    {@; return *(khc_p==nullptr ? khc_p=new KhatContext(val) : khc_p); }
  const Rep_context& real_form_value::rc()
    {@; return *(rt_p==nullptr ? rt_p=new Rep_table(val) : rt_p); }
  Rep_table& real_form_value::rt()
    {@; return *(rt_p==nullptr ? rt_p=new Rep_table(val) : rt_p); }
@)
  real_form_value::~real_form_value () @+{@; delete khc_p; delete rt_p; }

@ When printing a real form, we give the name by which it is known in the
parent inner class, and provide some information about its connectivity.
Since the names of the real forms are indexed by their outer number, but the
real form itself stores its inner number, we must somewhat laboriously make
the conversion here.

@< Function def...@>=
void real_form_value::print(std::ostream& out) const
{ if (val.isCompact()) out << "compact ";
  out << (val.isConnected() ? "connected " : "disconnected " );
  if (val.isQuasisplit())
    out << (val.isSplit() ? "" : "quasi") << "split ";
  out << "real group with Lie algebra '" @|
      << parent.interface.type_name(parent.interface.out(val.realForm())) @|
      << '\'' ;
}

@ To make a real form is easy, one provides an |inner_class_value| and a valid
index into its list of real forms. Since this number coming from the outside
is to be interpreted as an outer index, we must convert it to an inner index
at this point. We also inversely allow finding the index of a real form in its
inner class; of course the opposite conversion applies here.
As a special case we also provide the quasisplit form.

@< Local function def...@>=
void real_form_wrapper(expression_base::level l)
{ int i = get<int_value>()->int_val();
  shared_inner_class G = get<inner_class_value>();
  if (size_t(i)>=G->val.numRealForms())
    throw runtime_error ("Illegal real form number: "+str(i));
@.Illegal real form number@>
  if (l!=expression_base::no_value)
    push_value(std::make_shared<real_form_value>(*G,G->interface.in(i)));
}
@)
void form_number_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>@|
      (rf->parent.interface.out(rf->val.realForm())));
}
@)
void quasisplit_form_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<real_form_value>(*G,G->val.quasisplit()));
}

@*2 Functions operating on real reductive groups.
%
From a real reductive group we can go back to its inner class
@< Local function def...@>=
void inner_class_of_real_form_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<inner_class_value>(rf->parent));
}

void real_form_to_inner_class_coercion()
{@; inner_class_of_real_form_wrapper(expression_base::single_value); }

void real_form_to_root_datum_coercion()
{ shared_real_form rf= get<real_form_value>();
  push_value(root_datum_value::build(rf->parent.val.rootDatum()));
}

@ Here is a function that gives information about the dual component group
of a real reductive group: it returns the rank~$r$ such that the component
group is isomorphic to $(\Zee/2\Zee)^r$.

@< Local function def...@>=
void components_rank_wrapper(expression_base::level l)
{ shared_real_form R= get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(R->val.dualComponentReps().size()));
}

@ And here is one that counts the number of Cartan classes for the real form.

@< Local function def...@>=
void count_Cartans_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(rf->val.numCartan()));
}

@ The size of the finite set $K\backslash G/B$ can be determined from the real
form.

@< Local function def...@>=
void KGB_size_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(rf->val.KGB_size()));
}

@ Here is a somewhat technical function that will facilitate working ``in
coordinates'' with KGB elements for this real form. It returns a rational
coweight that determines the base grading for the real form, and which is an
offset that will be added to |torus_bits| values when computing the
|torus_factor| they represent. It is defined so that a zero value corresponds
to a quasisplit real form, which proves the most useful base point. This does
imply that a standard choice for the ``infinitesimal cocharacter'' for the
real from must differ by ${}^\vee\!\rho$ from the value produced here (and the
|cocharacter| field that stores it does not really have the right name).

We used to ensure that the coweight returned is dominant, while remaining
in the coset of $2X_*$ defined by |rf->cocharacter| so that the torus element
$\exp(\pi\ii t)$ giving the actual base grading is unaffected. There is
however no obvious best way to do this, and the easy way that used to be
employed could give coweights rather far from inside the dominant chamber; in
addition this shift may destroy the ${}^t\xi$-invariance of |rf->cocharacter|
that is hard to reestablish without risk of leaving the coset. Therefore we
just return |rf->cocharacter| as it is stored.

@< Local function def...@>=
void base_grading_vector_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value>(rf->val.g_rho_check()));
}

@ There is a partial ordering on the Cartan classes defined for a real form. A
matrix for this partial ordering is computed by the function
|Cartan_order|, which more or less replaces the \.{corder} command in
\.{Fokko}.

@h "poset.h"
@< Local function def...@>=
void Cartan_order_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l==expression_base::no_value)
    return;
@)
  size_t n=rf->val.numCartan();
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(n,n,0));
  const poset::Poset& p = rf->val.innerClass().Cartan_ordering();
  for (size_t i=0; i<n; ++i)
    for (size_t j=i; j<n; ++j)
      if (p.lesseq(i,j)) M->val(i,j)=1;

  push_value(std::move(M));
}

@ A similar function is |KGB_Hasse| which encodes the Bruhat order on the KGB
set as a matrix.

@h "bruhat.h"

@< Local function def...@>=
void KGB_Hasse_wrapper(expression_base::level l)
{ own_real_form rf= non_const_get<real_form_value>();
  if (l==expression_base::no_value)
    return;
@)
  size_t n=rf->val.KGB_size();
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(n,n,0));
  const auto& Bruhat = rf->val.Bruhat_KGB();
  for (size_t j=0; j<n; ++j)
  { const auto& col = Bruhat.hasse(j);
    for (auto it=col.begin(); it!=col.end(); ++it)
      M->val(*it,j)=1;
  }
  push_value(std::move(M));
}


@ Finally we make available an equality test for real forms. This could easily
be defined in the \.{axis} language itself: two real forms are equal if they
belong to the same inner class, their base grading vectors are equal, and the
torus bits of their initial KGB elements are the same (this should imply their
real form numbers are the same too). However it is useful to define a test
here; not only can we do the tests more efficiently, the same test will later
also be used in other equality tests (for KGB elements, or module parameters);
there one should resist the temptation to test (by pointer equality) for
identical |RealReductiveGroup| objects, which would be too strict. Except if
real form generation above were to be defined in such a way that in the
presence of an equal real form value in its inner class, a newly generated one
would actually return a reference to the old one, which would in fact be a
good idea (but would require additional administration).

@< Local function def...@>=

inline bool operator==
  (const RealReductiveGroup& x, const RealReductiveGroup& y)
{
  return &x.innerClass() == &y.innerClass()
   @| and x.g_rho_check()==y.g_rho_check()
   @| and x.x0_torus_part()==y.x0_torus_part()
   @| and (assert(x.realForm()==y.realForm()),true);
}
inline bool operator!=
  (const RealReductiveGroup& x, const RealReductiveGroup& y)
{@; return not(x==y); }
@)
void real_form_eq_wrapper(expression_base::level l)
{ shared_real_form y = get<real_form_value>();
  shared_real_form x = get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(whether(x->val==y->val));
}

void real_form_neq_wrapper(expression_base::level l)
{ shared_real_form y = get<real_form_value>();
  shared_real_form x = get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(whether(x->val!=y->val));
}

@*2 Dual real forms.
It will be convenient to be able to construct dual real forms in the context
of the |inner_class_value| itself, rather than having to construct the dual
inner class for this purpose. Originally the resulting values were given a
different type than ordinary real forms, but the slight chance of catching
errors (for instance when trying to make a block for a pair of real forms
rather than for a real form and a dual real form) via the type system, rather
than through a runtime check, did not in the end justify the additional
complications this gave to studying the dual situation (and with explicit type
conversions from real form to dual form and vice versa, mentioned type safety
was emptied of it substance).

To make a dual real form, one provides an |inner_class_value| and a valid
index into its list of dual real forms, which will be converted to an inner
index. We also provide the dual quasisplit form.

@< Local function def...@>=
void dual_real_form_wrapper(expression_base::level l)
{ int i =get<int_value>()->int_val();
  shared_inner_class G = get<inner_class_value>();
  if (size_t(i)>=G->val.numDualRealForms())
    throw runtime_error ("Illegal dual real form number: "+str(i));
@.Illegal dual real form number@>
  if (l==expression_base::no_value)
    return;
@)
  inner_class_value G_check(*G,tags::DualTag());
   // tailor make an |inner_class_value|
  push_value(std::make_shared<real_form_value>@|
    (G_check ,G->dual_interface.in(i)));
}
@)
void dual_quasisplit_form_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l==expression_base::no_value)
    return;
@)
  inner_class_value G_check(*G,tags::DualTag());
   // tailor make an |inner_class_value|
  push_value(std::make_shared<real_form_value>(G_check,G->dual.quasisplit()));
}

@*2 Synthetic real forms.
%
It is useful to be able to compute a real form based on other information than
its index within its inner class, namely on a strong involution representative
(involution and torus element). The synthetic \.{atlas} function |real_form|
takes an inner class, a matrix giving an involution~$\theta$, and a rational
coweight~$t$. After projecting $t$ to the $+1$ eigenspace of~$\theta$ to make
it $\theta$-stable, it describes (through $\exp_{-1}:l\mapsto\exp(\pi\ii l)$)
a torus element; the square of this element should be central (meaning in
coordinates that all simple roots have integral evaluation on the projected
rational vector; in the code below the doubled projection is made first, and
evaluations must be even). If this test succeeds, the function then returns
the corresponding real form, but in which a cocharacter value is stored
(representing the square of any strong involution for the real form) deduced
from the given torus element, which may differ from the value for the (weak)
real form selected by the number |rf->val.realForm()| in the inner class. This
difference notably allows the same strong involution representative to be
subsequently used to specify a KGB element for this (strong) real form.

@:synthetic_real_form@>

@< Local function def...@>=
TwistedInvolution twisted_from_involution
  (const InnerClass& G, const WeightInvolution theta0)
{ const RootDatum& rd = G.rootDatum();
  WeightInvolution theta(theta0); // copy allowing modification
  WeylWord ww;
  if (check_involution(theta,rd,ww)!=G.twistedWeylGroup().twist() @| or
      theta!=G.distinguished())
    throw runtime_error("Involution not in this inner class");
  return G.weylGroup().element(ww);
}
@)
void synthetic_real_form_wrapper(expression_base::level l)
{ own_rational_vector torus_factor = get_own<rational_vector_value>();
  shared_matrix theta = get<matrix_value>();
  shared_inner_class G = get<inner_class_value>();
  TwistedInvolution tw = twisted_from_involution(G->val,theta->val);
  if (torus_factor->val.size()!=G->val.rank())
    throw runtime_error ("Torus factor size mismatch");
  {
    Ratvec_Numer_t& num = torus_factor->val.numerator();
    num += theta->val.right_prod(num);
      // make torus factor $\theta$-fixed, temporarily doubled
    TorusElement t(torus_factor->val,false);
      // take a copy as |TorusElement|, using $\exp_{-1}$
    const RootDatum& rd = G->val.rootDatum();
    LatticeMatrix alpha
      (rd.beginSimpleRoot(),rd.endSimpleRoot(),rd.rank(),tags::IteratorTag());
    if (not is_central(alpha,t)) // every root should now have even evaluation
      throw runtime_error
         ("Torus factor does not define a valid strong involution");
@.Torus factor does not...@>
    torus_factor->val /= 2; // now $(1+\theta)/2$ is applied to |torus_factor|
  }

  @< Ensure that the involution table knows about the Cartan class of |tw|
     and those below it in the Cartan ordering @>
  RatCoweight coch(0); // dummy value to be replaced
  RealFormNbr rf = real_form_of(G->val,tw,torus_factor->val,coch);
   // sets |coch|
  TorusPart tp = realredgp::minimal_torus_part @|
   (G->val,rf,coch,std::move(tw),torus_factor->val);
  if (l!=expression_base::no_value)
    push_value(std::make_shared<real_form_value>(*G,rf,coch,tp));
}

@ The call to |minimal_torus_part| uses the involution table in order to be
able to do downward (inverse) Cayley transforms, so we must ensure that any
involutions that can be encountered have been entered into the table.

@< Ensure that the involution table knows about the Cartan class of |tw|... @>=
{ CartanNbr cn = G->val.class_number(tw);
  G->val.generate_Cartan_orbit(cn);
  const BitMap& b = G->val.Cartan_ordering().below(cn);
  for (auto it=b.begin(); it(); ++it)
    G->val.generate_Cartan_orbit(*it);
}

@ The methods |central_fiber| and |x0_torus_part| of |InnerClass|
can be accessed using following functions. The function |central_fiber|
computes those torus parts in the fiber at the distinguished involution that
both remain in the strong real form orbit and are central (do not affect any
gradings).

@< Local function def...@>=
void central_fiber_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l==expression_base::no_value)
    return;
@)
  auto cf = rf->parent.val.central_fiber(rf->val.realForm());
  own_row result = std::make_shared<row_value>(cf.size());
  unsigned int i=0;
  for (auto it=cf.begin(); it!=cf.end(); ++it, ++i)
    result->val[i]= std::make_shared<vector_value>(int_Vector(*it));
  push_value(result);
}

void initial_torus_bits_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value> @|
      (int_Vector(rf->parent.val.x0_torus_part(rf->val.realForm()))));
}


@ Finally we install everything related to real forms.
@< Install wrapper functions @>=
install_function(real_form_wrapper,@|"real_form","(InnerClass,int->RealForm)");
install_function(form_number_wrapper,@|"form_number","(RealForm->int)");
install_function(quasisplit_form_wrapper,@|"quasisplit_form"
		,"(InnerClass->RealForm)");
install_function(inner_class_of_real_form_wrapper
                ,@|"inner_class","(RealForm->InnerClass)");
install_function(components_rank_wrapper,@|"components_rank","(RealForm->int)");
install_function(count_Cartans_wrapper,@|"count_Cartans","(RealForm->int)");
install_function(KGB_size_wrapper,@|"KGB_size","(RealForm->int)");
install_function(base_grading_vector_wrapper
                ,@|"base_grading_vector","(RealForm->ratvec)");
install_function(Cartan_order_wrapper,@|"Cartan_order","(RealForm->mat)");
install_function(KGB_Hasse_wrapper,@|"KGB_Hasse","(RealForm->mat)");
install_function(real_form_eq_wrapper,"=","(RealForm,RealForm->bool)");
install_function(real_form_neq_wrapper,"!=","(RealForm,RealForm->bool)");
install_function(dual_real_form_wrapper,@|"dual_real_form"
				       ,"(InnerClass,int->RealForm)");
install_function(dual_quasisplit_form_wrapper,@|"dual_quasisplit_form"
		,"(InnerClass->RealForm)");
install_function(synthetic_real_form_wrapper,@|"real_form"
		,"(InnerClass,mat,ratvec->RealForm)");
install_function(central_fiber_wrapper,"central_fiber","(RealForm->[vec])");
install_function(initial_torus_bits_wrapper,@|"initial_torus_bits"
                ,"(RealForm->vec)");

@*1 A type for Cartan classes.
%
Another type of value associated to inner classes are Cartan classes, which
describe equivalence classes (for ``stable conjugacy'') of Cartan subgroups of
real reductive groups in an inner class. The Atlas software associates a fixed
set of Cartan classes to each inner class, and records for each real form the
subset of those Cartan classes that occur for the real form. Since versions
0.3.5 of the software, the Cartan classes are identified and numbered upon
construction of an |InnerClass| object, but |CartanClass| objects
are constructed on demand.

@< Includes... @>=
#include "cartanclass.h"

@*2 Class definition.
%
The Cartan class is identified by a number within the inner class. This number
is its sequence number in the order in which Cartan classes are generated in
this inner class.

@< Type definitions @>=
struct Cartan_class_value : public value_base
{ const inner_class_value parent;
  size_t number;
  const CartanClass& val;
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
  : parent(v.parent), number(v.number), val(v.val) @+{} // copy constructor
};
@)
typedef std::shared_ptr<const Cartan_class_value> shared_Cartan_class;

@ In the constructor we used to check that the Cartan class with the given
number currently exists, but now the |InnerClass::cartan| method
assures that one is generated if this should not have been done before. We
therefore call that method in the initialiser; on return it provides a valid
reference.

@< Function def...@>=
Cartan_class_value::Cartan_class_value(const inner_class_value& p,size_t cn)
: parent(p),number(cn),val(p.val.cartan(cn))
@+{}

@ When printing a Cartan class, we show its number, and for how many real
forms and dual real forms it is valid.

@< Function def...@>=
void Cartan_class_value::print(std::ostream& out) const
{ out << "Cartan class #" << number << ", occurring for " @|
  << val.numRealForms() << " real " @|
  << (val.numRealForms()==1 ? "form" : "forms") << " and for "@|
  << val.numDualRealForms() << " dual real "  @|
  << (val.numDualRealForms()==1 ? "form" : "forms");
}

@ To make a Cartan class, one can provide a |inner_class_value| together with
a valid index into its list of Cartan classes.

@< Local function def...@>=
void ic_Cartan_class_wrapper(expression_base::level l)
{ int i=get<int_value>()->int_val();
  shared_inner_class ic = get<inner_class_value>();
  if (size_t(i)>=ic->val.numCartanClasses())
    throw runtime_error ("Illegal Cartan class number: "+str(i)
@.Illegal Cartan class number@>
    +", this inner class only has "+str(ic->val.numCartanClasses())
    +" of them");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<Cartan_class_value>(*ic,i));
}

@ Alternatively (and this used to be the only way) one can provide a
|real_form_value| together with a valid index into \emph{its} list of Cartan
classes. We translate this number into an index into the list for its
containing inner class, and then get the Cartan class from there.

@< Local function def...@>=
void rf_Cartan_class_wrapper(expression_base::level l)
{ int i=get<int_value>()->int_val();
  shared_real_form rf= get<real_form_value>();
  if (size_t(i)>=rf->val.numCartan())
    throw runtime_error ("Illegal Cartan class number: "+str(i)
@.Illegal Cartan class number@>
    +", this real form only has "+str(rf->val.numCartan())+" of them");
  BitMap cs=rf->val.Cartan_set();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<Cartan_class_value>(rf->parent,cs.n_th(i)));
}

@ Like the quasisplit real form for inner classes, there is a particular
Cartan class for each real form, the ``most split Cartan''; it is of
particular importance for knowing which dual real forms can be associated to
this real form. It is (probably) the last Cartan class in the list for the
real form, but we have a direct access to it via the |mostSplit| method for
|RealReductiveGroup|.

@< Local function def...@>=
void most_split_Cartan_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<Cartan_class_value>
      (rf->parent,rf->val.mostSplit()));
}



@*2 Functions operating on Cartan classes.
%
We start with a fundamental attribute of Cartan classes: the associated
(distinguished) involution, in the form of a matrix.

@< Local function def...@>=
void Cartan_involution_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value>(cc->val.involution()));
}


@ This function and the following provide the functionality of the \.{Fokko}
command \.{cartan}. They are based on |output::printCartanClass|, but
rewritten to take into account the fact that we do not know about |Interface|
objects for complex groups, and such that a usable value is returned. We omit
in our function {\it Cartan\_info} the data printed in |printCartanClass| in
the final call to |output::printFiber| (namely all the real forms for which
this Cartan class exists with the corresponding part of the adjoint fiber
group), relegating it instead to a second function {\it fiber\_part} that
operates on a per-real-form basis. This separation seems more natural in a
setup where real forms and Cartan classes are presented as more or less
independent objects.

@h "prettyprint.h"

@< Local function def...@>=
void Cartan_info_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  if (l==expression_base::no_value)
    return;
@)
  auto ranks = tori::classify(cc->val.involution());

  push_value(std::make_shared<int_value>(std::get<0>(ranks)));
  push_value(std::make_shared<int_value>(std::get<1>(ranks)));
  push_value(std::make_shared<int_value>(std::get<2>(ranks)));
  wrap_tuple<3>();

  const weyl::TwistedInvolution& tw =
    cc->parent.val.involution_of_Cartan(cc->number);
  WeylWord ww = cc->parent.val.weylGroup().word(tw);

  std::vector<int> v(ww.begin(),ww.end());
  push_value(std::make_shared<vector_value>(v));

  push_value(std::make_shared<int_value>(cc->val.orbitSize()));
  push_value(std::make_shared<int_value>(cc->val.fiber().fiberSize()));
  wrap_tuple<2>();

  const RootSystem& rs=cc->parent.val.rootDatum();

@)// print types of imaginary and real root systems and of Complex factor
  push_value(std::make_shared<Lie_type_value> @|
    (rs.subsystem_type(cc->val.simpleImaginary())));
  push_value(std::make_shared<Lie_type_value> @|
    (rs.subsystem_type(cc->val.simpleReal())));
  push_value(std::make_shared<Lie_type_value> @|
    (rs.subsystem_type(cc->val.simpleComplex())));
  wrap_tuple<3>();
@)
  if (l==expression_base::single_value)
    wrap_tuple<4>();
}

@ A functionality that is implicit in the \.{Fokko} command \.{cartan} is the
enumeration of all real forms corresponding to a given Cartan class. While
|output::printFiber| traverses the real forms in the order corresponding to
the parts of the partition |f.weakReal()| for the fiber~|f| associated to the
Cartan class, it is not necessary to use this order, and we can instead simply
traverse all real forms and check whether the given Cartan class exists for
them. Taking care to convert real form numbers to their inner representation
here, we can in fact return a list of real forms.

We also include a version for dual real forms. Corresponding Cartan classes in
the dual inner class have a different (opposite) numbering; therefore the test
for |cc->number| should use the |dual_Cartan_set| of the current inner class
rather than the |Cartan_set| of the just constructed dual inner class.

@< Local function def...@>=
void real_forms_of_Cartan_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  if (l==expression_base::no_value)
    return;
@)
  const inner_class_value& ic=cc->parent;
  own_row result = std::make_shared<row_value>(cc->val.numRealForms());
  for (size_t i=0,k=0; i<ic.val.numRealForms(); ++i)
  { RealFormNbr rf = ic.interface.in(i);
    BitMap b(ic.val.Cartan_set(rf));
    if (b.isMember(cc->number))
      result->val[k++] = std::make_shared<real_form_value>(ic,rf);
  }
  push_value(std::move(result));
}
@)
void dual_real_forms_of_Cartan_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  if (l==expression_base::no_value)
    return;
@)
  const inner_class_value dual_ic(cc->parent,tags::DualTag());
  own_row result = std::make_shared<row_value>(cc->val.numDualRealForms());
  for (size_t i=0,k=0; i<dual_ic.val.numRealForms(); ++i)
  { RealFormNbr drf = cc->parent.dual_interface.in(i);
    BitMap b (cc->parent.val.dual_Cartan_set(drf));
    if (b.isMember(cc->number))
      result->val[k++] = std::make_shared<real_form_value>(dual_ic,drf);
  }
  push_value(std::move(result));
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

@< Local function def...@>=
void fiber_partition_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  shared_Cartan_class cc(get<Cartan_class_value>());
  if (&rf->parent.val!=&cc->parent.val)
    throw runtime_error
    ("Inner class mismatch between real form and Cartan class");
@.Inner class mismatch...@>
  BitMap b(cc->parent.val.Cartan_set(rf->val.realForm()));
  if (!b.isMember(cc->number))
    throw runtime_error
    ("Cartan class not defined for this real form");
@.Cartan class not defined@>
  if (l==expression_base::no_value)
    return;
@)
  const Partition& pi = cc->val.fiber().weakReal();
  const RealFormNbrList rf_nr=
     cc->parent.val.realFormLabels(cc->number);
     // translate part number of |pi| to real form
  own_row result =
    std::make_shared<row_value>(0); // cannot predict exact size here
  for (size_t i=0; i<pi.size(); ++i)
    if (rf_nr[pi.class_of(i)] == rf->val.realForm())
      result->val.push_back(std::make_shared<int_value>(i));
  push_value(std::move(result));
}

@ The function |square_classes| returns the set of real forms associated to a
Cartan class, but partitioned according to their square class, determined by
the square of any strong involution representing the real form.
@< Local function def...@>=
void square_classes_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  const output::FormNumberMap rfi = cc->parent.interface;
  const RealFormNbrList& rfl = cc->parent.val.realFormLabels(cc->number);
  if (l==expression_base::no_value)
    return;
@)
  size_t n_sq_classes = cc->val.numRealFormClasses();
  own_row result = std::make_shared<row_value>(n_sq_classes);
  for (cartanclass::square_class csc=0; csc<n_sq_classes; ++csc)
  { const Partition& pi = cc->val.fiber_partition(csc);
    own_row part = std::make_shared<row_value>(pi.classCount());
    for (unsigned long c=0; c<pi.classCount(); ++c)
       part->val[c] =
          std::make_shared<int_value>(rfi.out(rfl[cc->val.toWeakReal(c,csc)]));
    result->val[csc] = std::move(part);
  }
  push_value(std::move(result));
}

@ The function |print_gradings| gives on a per-real-form basis the functionality
of the \.{Fokko} command \.{gradings} that is implemented by
|output::printGradings| and |output::printGradings|. It therefore takes, like
|fiber_partition|, a Cartan class and a real form as parameter. Its output
consist of a list of $\Zee/2\Zee$-gradings of each of the fiber group elements
in the part corresponding to the real form, where each grading is a sequence of
bits corresponding to the simple imaginary roots.

@f sigma nullptr

@< Local function def...@>=
void print_gradings_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
@/shared_Cartan_class cc(get<Cartan_class_value>());
  if (&rf->parent.val!=&cc->parent.val)
    throw runtime_error
    ("Inner class mismatch between real form and Cartan class");
@.Inner class mismatch...@>
  BitMap b(cc->parent.val.Cartan_set(rf->val.realForm()));
  if (!b.isMember(cc->number))
    throw runtime_error ("Cartan class not defined for this real form");
@.Cartan class not defined...@>
@)
  const Partition& pi = cc->val.fiber().weakReal();
  const RealFormNbrList rf_nr=
     cc->parent.val.realFormLabels(cc->number);
     // translate part number of |pi| to real form
@)
  const RootNbrList& si = cc->val.fiber().simpleImaginary();
  // simple imaginary roots

  int_Matrix cm;
  Permutation sigma;
@/@< Compute the Cartan matrix |cm| of the root subsystem |si|, and the
     permutation |sigma| giving the Bourbaki numbering of its simple roots @>

  @< Print information about the imaginary root system and its simple roots @>

  @< Print the gradings for the part of |pi| corresponding to the real form @>

  if (l==expression_base::single_value)
    wrap_tuple<0>(); // |no_value| needs no special care
}

@ For normalising the ordering of the simple imaginary roots we use two
functions from \.{dynkin.cpp}.

@< Compute the Cartan matrix |cm|... @>=
{ cm=cc->parent.val.rootDatum().cartanMatrix(si);
  dynkin::DynkinDiagram d(cm); sigma = dynkin::bourbaki(d);
}

@ The imaginary root system might well be empty, so we make special provisions
for this case.

@< Print information about the imaginary root system... @>=
{ *output_stream << "Imaginary root system is ";
  if (si.size()==0) *output_stream<<"empty.\n";
  else
  { LieType t = dynkin::Lie_type(cm);

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

@h "ioutils.h"
@< Print the gradings for the part of |pi|... @>=
{ bool first=true; std::ostringstream os;
  for (size_t i=0; i<pi.size(); ++i)
    if ( rf_nr[pi.class_of(i)] == rf->val.realForm())
    { const auto& f= cc->val.fiber();
      cartanclass::AdjointFiberElt afe(RankFlags(i),f.adjointFiberRank());
      os << ( first ? first=false,'[' : ',');
      Grading gr=f.grading(afe);
      gr=sigma.pull_back(gr);
      prettyprint::prettyPrint(os,gr,si.size());
    }
  os << "]" << std::endl;
  ioutils::foldLine(*output_stream,os.str(),"",",");
}

@ Finally we install everything related to Cartan classes.
@< Install wrapper functions @>=
install_function(ic_Cartan_class_wrapper,@|"Cartan_class"
		,"(InnerClass,int->CartanClass)");
install_function(rf_Cartan_class_wrapper,@|"Cartan_class"
		,"(RealForm,int->CartanClass)");
install_function(most_split_Cartan_wrapper,@|"most_split_Cartan"
		,"(RealForm->CartanClass)");
install_function(Cartan_involution_wrapper,@|"involution","(CartanClass->mat)");
install_function(Cartan_info_wrapper,@|"Cartan_info"
		,"(CartanClass->(int,int,int),"
                 "vec,(int,int),(LieType,LieType,LieType))");
install_function(real_forms_of_Cartan_wrapper,@|"real_forms"
		,"(CartanClass->[RealForm])");
install_function(dual_real_forms_of_Cartan_wrapper,@|"dual_real_forms"
		,"(CartanClass->[RealForm])");
install_function(fiber_partition_wrapper,@|"fiber_partition"
		,"(CartanClass,RealForm->[int])");
install_function(square_classes_wrapper,@|"square_classes"
                ,"(CartanClass->[[int]])");

@*1 Elements of some $K\backslash G/B$ set.
%
Associated to each real form is a set $K\backslash G/B$, which was already
made available internally through the |real_form_value::kgb| method. We wish
to make available externally a number of functions that operate on (among
other data) elements of this set; for instance each element determines an
involution, corresponding imaginary and real subsystems of the root system,
and a $\Zee/2\Zee$-grading of the imaginary root subsystem. It does not seem
necessary to introduce a type for the set $K\backslash G/B$ (any function that
needs it can take the corresponding real form as argument), but we do provide
one for elements of that set (if desired the user can collect such elements in
an array). These elements contain a pointer |rf | to the |real_form_value|
they were generated from, and thereby to the associated set $K\backslash G/B$,
which allows operations relevant to the KGB element to be defined. Since this
is a shared pointer, the |real_form_value| pointed to is guaranteed to exist
as long as this |KGB_elt_value| does.

@< Type definitions @>=
struct KGB_elt_value : public value_base
{ own_real_form rf;
  KGBElt val;
@)
  KGB_elt_value(const own_real_form& form, KGBElt x) : rf(form), val(x) @+{}
  ~KGB_elt_value() @+{}
@)
  virtual void print(std::ostream& out) const;
  KGB_elt_value* clone() const @+ {@; return new KGB_elt_value(*this); }
  static const char* name() @+{@; return "KGB element"; }
private:
  KGB_elt_value(const KGB_elt_value& v)
  : rf(v.rf), val(v.val) @+{} // copy constructor
};
@)
typedef std::unique_ptr<KGB_elt_value> KGB_elt_ptr;
typedef std::shared_ptr<const KGB_elt_value> shared_KGB_elt;
typedef std::shared_ptr<KGB_elt_value> own_KGB_elt;

@ When printing a KGB element, we print the number. It would be useful to add
some symbolic information, like ``discrete series'', when applicable, but this
is currently not implemented.

@< Function def...@>=
void KGB_elt_value::print(std::ostream& out) const
@+{@; out << "KGB element #" << val;
}

@ To make a KGB element, one provides a |real_form_value| and a valid number.

@< Local function def...@>=
void KGB_elt_wrapper(expression_base::level l)
{ int i = get<int_value>()->int_val();
  own_real_form rf= non_const_get<real_form_value>();
  if (size_t(i)>=rf->val.KGB_size())
    throw runtime_error ("Inexistent KGB element: "+str(i));
@.Inexistent KGB element@>
  if (l!=expression_base::no_value)
    push_value(std::make_shared<KGB_elt_value>(rf,i));
}

@ Working with KGB elements often requires having access to its real form, and
sometimes of its number within its real form.

@< Local function def...@>=
void decompose_KGB_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  if (l==expression_base::no_value)
    return;
@)
  push_value(x->rf);
  push_value(std::make_shared<int_value>(x->val));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ Three important attributes of KGB elements are the associated Cartan class,
the root datum involution and the length.

@< Local function def...@>=
void KGB_Cartan_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<Cartan_class_value>
      (x->rf->parent,kgb.Cartan_class(x->val)));
}

void KGB_involution_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  const InnerClass& G=x->rf->val.innerClass();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<matrix_value>
      (G.matrix(kgb.involution(x->val))));
}

void KGB_length_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(kgb.length(x->val)));
}

@ Cross actions and (inverse) Cayley transforms define the structure of a KGB
set, and we make them available as functions. The inverse Cayley transform may
be double valued of which we only report the first one; the user can easily
test whether it was double valued by applying cross action by the same
generator to the result and testing whether it gives a new KGB element (which
then is the other value of the Cayley transform). Therefore we define
functions that are single-valued in all cases (for undefined cases we just
return the argument KGB element). Given that, there is not much reason to
distinguish forward and inverse Cayley transforms (which have disjoint domains
of definition), so we combine them into one function |Cayley|.

@< Local function def...@>=
inline RootNbr get_reflection_index(int root_index, RootNbr n_posroots)
{ RootNbr alpha= root_index<0 ? -1-root_index : root_index;
  if (alpha>=n_posroots)
    throw runtime_error ("Illegal reflection: "+str(root_index));
  return alpha;
}
@)
void KGB_cross_wrapper(expression_base::level l)
{ own_KGB_elt x = get_own<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  RootNbr npr=kgb.rootDatum().numPosRoots();
  RootNbr alpha = get_reflection_index(get<int_value>()->int_val(),npr);
  if (l==expression_base::no_value)
    return;
@)
  if (alpha<kgb.rank()) // do simple cross action
    x->val= kgb.cross(alpha,x->val);
  else // do non-simple cross action
    x->val = cross(kgb,x->val,npr+alpha);
  push_value(x);
}
@)
void KGB_Cayley_wrapper(expression_base::level l)
{ own_KGB_elt x = get_own<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  RootNbr npr=kgb.rootDatum().numPosRoots();
  RootNbr alpha = get_reflection_index(get<int_value>()->int_val(),npr);
  if (l==expression_base::no_value)
    return;
@)
  if (alpha<kgb.rank())
  {
     KGBElt xv = kgb.any_Cayley(alpha,x->val);
     if (xv!=UndefKGB)
       x->val= xv; // when defined do Cayley transform
  }
  else // do (inverse) Cayley transform through arbitrary root
  { try
    {@; x->val= any_Cayley(kgb,x->val,npr+alpha); }
    catch (std::runtime_error&) {}
      // ignore undefined Cayley error, leave |x| unchanged
  }
  push_value(x);
}

@ One also needs to be able find out the status of roots. Although somewhat
low-level, the simplest thing is to export the |gradings::Status::Value| as an
integer value. It seems however useful to change complex \emph{ascents} from
$0$ to $4$, so that the coding is 0:~Complex descent, 1:~imaginary compact,
2:~real, 3:~imaginary non-compact, 4:Complex ascent. This way a value $v$ is a
descent if |v<3|, imaginary if |v%2==1|, Complex if |v%4==0|, Cayley transform
defined if |v==3|.

@< Local function def...@>=
void KGB_status_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  RootNbr npr=kgb.rootDatum().numPosRoots();
  RootNbr alpha = get_reflection_index(get<int_value>()->int_val(),npr);
  if (l==expression_base::no_value)
    return;
@)
  if (alpha<kgb.rank())
  {
    unsigned stat=kgb.status(alpha,x->val);
    push_value(std::make_shared<int_value>
      (stat==0 and not kgb.isDescent(alpha,x->val) ? 4 : stat));
  }
  else
  {
    alpha += npr; // convert to general root number
    unsigned stat=kgb::status(kgb,x->val,alpha);
    if (stat==0) // $\alpha$ is a complex root, check if it is an ascent
    {
      RootNbr theta_alpha = kgb.innerClass().involution_table().
        root_involution(kgb.inv_nr(x->val),alpha);
      if (kgb.rootDatum().is_posroot(theta_alpha))
       stat = 4; // set status to complex ascent
    }
    push_value(std::make_shared<int_value> (stat));
  }
}

@ In order to ``synthesise'' a KGB element, one may specify a real form, an
involution, and a rational weight that should be the |torus_factor| value. The
latter defines a grading of the corresponding imaginary roots, in the same
manner as for the synthetic |real_form| in section @#synthetic_real_form@>
above, but in fact it even completely describes the KGB element. In order for
this to be possible, the |torus_factor| must be compatible with the
cocharacter |rf->val.g_rho_check()| stored in the real form, but which should
always be right if the real form was itself synthesised from the
|torus_factor| value.

@< Local function def...@>=
void build_KGB_element_wrapper(expression_base::level l)
{ own_rational_vector torus_factor = get_own<rational_vector_value>();
  shared_matrix theta = get<matrix_value>();
  own_real_form rf = non_const_get<real_form_value>();

  if (torus_factor->val.size()!=rf->val.rank())
    throw runtime_error ("Torus factor size mismatch");
@)
  Ratvec_Numer_t& num = torus_factor->val.numerator();
  { // make theta-fixed and remove base grading vector offset
    num += theta->val.right_prod(num);
    ((torus_factor->val /= 2) -=rf->val.g_rho_check()).normalize() ;
    if (torus_factor->val.denominator()!=1)
      throw runtime_error
        ("Torus factor not in cocharacter coset of real form");
@.Torus factor not in cocharacter...@>
  }

  const InnerClass& G = rf->parent.val;
  TitsElt a
   (G.titsGroup(),TorusPart(num),twisted_from_involution(G,theta->val));

  KGBElt x = rf->kgb().lookup(a);
  if (x==UndefKGB)
    throw runtime_error("KGB element not present");

  if (l!=expression_base::no_value)
    push_value(std::make_shared<KGB_elt_value>(rf,x));
}



@ One can conjugate a KGB element by the distinguished involution of the inner
class, or by an explicitly supplied distinguished involution that
should commute with the inner class one; |test_compatible| tests this.

@< Local function def...@>=
void KGB_twist_wrapper(expression_base::level l)
{ own_KGB_elt x = get_own<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  if (l==expression_base::no_value)
    return;
@)
  x->val= kgb.Hermitian_dual(x->val); // do twist
  push_value(x);
}
@)
void test_compatible (const InnerClass& ic, shared_matrix& delta)
{ check_based_root_datum_involution(ic.rootDatum(),delta->val);
  const auto& xi = ic.distinguished();
  if (delta->val*xi!=xi*delta->val)
    throw runtime_error("Non commuting distinguished involution");
}
@)
void KGB_outer_twist_wrapper(expression_base::level l)
{ auto delta = get<matrix_value>();
  own_KGB_elt x = get_own<KGB_elt_value>();
  test_compatible(x->rf->parent.val,delta);
  if (l==expression_base::no_value)
    return;
@)
  x->val= x->rf->kgb().twisted(x->val,delta->val); // do twist
  push_value(x);
}

@ Here is a function that returns the vector of bits that distinguish KGB
elements in the fibre over the same involution.

@< Local function def...@>=
void torus_bits_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  if (l==expression_base::no_value)
    return;
@)
  const KGB& kgb=x->rf->kgb();
  TorusPart t = kgb.torus_part(x->val);
  own_vector result = std::make_shared<vector_value>(lift(t));
  push_value(std::move(result));
}

@ It is often be more useful to obtain an equivalent value, but in a form that
takes into account the cocharacter stored for the real form. The following
function does that, returning the result in the form of a rational vector that
should be interpreted in $(\Qu/2\Zee)^n$ representing
$(X_*)\otimes_\Zee\Qu/2X_*$. The value has two useful properties: it is
invariant under right-multiplication by the involution~$\theta$ associated
to~$x$ (i.e., the ignored part has been projected away), and at the same time it
lies in the $X_*$-coset of |cocharacter| for the real form.

@< Local function def...@>=
void torus_factor_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value> @|
      (x->rf->kgb().torus_factor(x->val)));
}

@ Finally we make available, by popular request, the equality test. This is
straightforward as long as one takes care not to just test the (smart)
pointers |x->rf| and |y->rf| to the |real_form_value| objects for equality,
but to test their |val| fields (of type |RealReductiveGroup|) using the
equality operator defined above, which test the |InnerClass|
objects for identity, and real form numbers for equality. This distinction is
important to make synthesised real forms first class citizens.

@< Local function def...@>=
void KGB_eq_wrapper(expression_base::level l)
{ shared_KGB_elt y = get<KGB_elt_value>();
  shared_KGB_elt x = get<KGB_elt_value>();
  if (l!=expression_base::no_value)
    push_value(whether(x->rf->val==y->rf->val and x->val==y->val));
}

void KGB_neq_wrapper(expression_base::level l)
{ shared_KGB_elt y = get<KGB_elt_value>();
  shared_KGB_elt x = get<KGB_elt_value>();
  if (l!=expression_base::no_value)
    push_value(whether(x->rf->val!=y->rf->val or x->val!=y->val));
}

@ Finally we install everything related to $K\backslash G/B$ elements.
@< Install wrapper functions @>=
install_function(KGB_elt_wrapper,@|"KGB","(RealForm,int->KGBElt)");
install_function(decompose_KGB_wrapper,@|"%","(KGBElt->RealForm,int)");
install_function(KGB_cross_wrapper,@|"cross","(int,KGBElt->KGBElt)");
install_function(KGB_Cayley_wrapper,@|"Cayley","(int,KGBElt->KGBElt)");
install_function(KGB_status_wrapper,@|"status","(int,KGBElt->int)");
install_function(build_KGB_element_wrapper,@|"KGB_elt"
		,"(RealForm,mat,ratvec->KGBElt)");
install_function(KGB_twist_wrapper,@|"twist","(KGBElt->KGBElt)");
install_function(KGB_outer_twist_wrapper,@|"twist","(KGBElt,mat->KGBElt)");
install_function(KGB_Cartan_wrapper,@|"Cartan_class","(KGBElt->CartanClass)");
install_function(KGB_involution_wrapper,@|"involution","(KGBElt->mat)");
install_function(KGB_length_wrapper,@|"length","(KGBElt->int)");
install_function(torus_bits_wrapper,@|"torus_bits","(KGBElt->vec)");
install_function(torus_factor_wrapper,@|"torus_factor","(KGBElt->ratvec)");
install_function(KGB_eq_wrapper,@|"=","(KGBElt,KGBElt->bool)");
install_function(KGB_neq_wrapper,@|"!=","(KGBElt,KGBElt->bool)");


@*1 Blocks associated to a real form and a dual real form.
%
Although blocks as specified by a real form and a dual real form were designed
more for the original \.{Fokko} interface than for use by \.{atlas}, we
provide a data type for such blocks and some simple functionality associated
to them.

@< Includes... @>=

#include "blocks.h"
#include "kl.h"

@ Like other data types we have seen, we include shared pointers to
parent objects to ensure these remain in existence as long as our block does;
in fact we include two such shared pointers, one for each real form. The |val|
field contains an actual |Block| instance, which is constructed van the
|Block_Value| is. We also reserve a field in the structure to store KL
polynomials, though they will only be computed once they are asked for.

@< Type definitions @>=
struct Block_value : public value_base
{ const own_real_form rf; const own_real_form dual_rf;
  Block val; // cannot be |const|, as Bruhat order may be generated implicitly
  kl::KLContext klc;
@)
  Block_value(const own_real_form& form, const own_real_form& dual_form);
  ~Block_value() @+{}
@)
  virtual void print(std::ostream& out) const;
  Block_value* clone() const @+ {@; return new Block_value(*this); }
  static const char* name() @+{@; return "KGB element"; }
private:
  Block_value(const Block_value& v)
  : rf(v.rf), val(v.val), klc(val) @+{} // copy constructor, but |klc| is new
};
@)
typedef std::unique_ptr<Block_value> Block_ptr;
typedef std::shared_ptr<const Block_value> shared_Block;
typedef std::shared_ptr<Block_value> own_Block;

@ The constructor for |Block_value| is relatively elaborate, so we lift it out
of the class declaration. One should avoid calling the version of the |build|
method that takes real form and deal real form numbers (as we originally did
here), as this will generate small versions of the KGB sets, which changes the
numbering and causes inconsistencies when KGB elements are extracted from
block elements. Instead we call the |build| method that accepts a pair of
|RealReductiveGroup| arguments, so that the generated block will be exactly a
classical one. Finally we associate the |klc| field with this block, but this
does not yet do much computation.

@< Function def...@>=
  Block_value::Block_value(const own_real_form& form,
                          const own_real_form& dual_form)
  : rf(form), dual_rf(dual_form)
  , val(Block::build(rf->val,dual_rf->val))
  , klc(val)
  {}


@ When printing a block, we print its size; we shall later provide a separate
print function that tabulates its individual elements.

@< Function def...@>=
void Block_value::print(std::ostream& out) const
@+{@; out << "Block of " << val.size() << " elements";
}

@ To make a block, one provides a |real_form_value| and a |dual_real_form|.
This function is named in a tribute to the creator of the Atlas software, and
of the data structure that is constructed and stored here.

@< Local function def...@>=
void Fokko_block_wrapper(expression_base::level l)
{ own_real_form drf=non_const_get<real_form_value>();
  own_real_form rf=non_const_get<real_form_value>();
@)
  if (&rf->parent.dual!=&drf->parent.val)
    throw runtime_error @|
    ("Inner class mismatch between real form and dual real form");
@.Inner class mismatch...@>
  BitMap b(rf->parent.val.dual_Cartan_set(drf->val.realForm()));
  if (not b.isMember(rf->val.mostSplit()))
    throw runtime_error @|
    ("Real form and dual real form are incompatible");
@.Real form and dual...@>
  if (l!=expression_base::no_value)
    push_value(std::make_shared<Block_value>(rf,drf));
}

@ We provide for now only a couple of basic function on the new type, which
will otherwise serve mainly as argument type for table-generating output
functions (for instance of KLV polynomials). We allow extracting the real
form and dual real forms used to construct the block, find the size of the
block, and allow to created from the number of a block element a pair of KGB
elements for the two forms that uniquely identifies the pair. Since both
components of the pair are KGB elements, but for real forms in mutually dual
inner classes, we must take care to create the dual inner class here, as is
done in |real_form_from_dual_real_form_wrapper|.

@< Local function def...@>=
void decompose_block_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  if (l==expression_base::no_value)
    return;
@)
  push_value(b->rf);
  push_value(b->dual_rf);
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

void size_of_block_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(b->val.size()));
}

void block_element_wrapper(expression_base::level l)
{ int i=get<int_value>()->int_val();
  shared_Block b = get<Block_value>();
  BlockElt z = i; // extract value unsigned
  if (z>=b->val.size())
    throw runtime_error
      ("Block element " +str(z) + " out of range (<" + str(b->val.size())+")");
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<KGB_elt_value>(b->rf,b->val.x(z)));
  inner_class_value dual_ic(b->rf->parent,tags::DualTag());
  own_real_form drf =
    std::make_shared<real_form_value>(dual_ic,b->dual_rf->val.realForm());
  push_value(std::make_shared<KGB_elt_value>(drf,b->val.y(z)));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ The inverse operation of decomposing a block element into a KGB element and
a dual KGB element is also useful. While one could construct the containing
|Block| value from the given values |x|, |y|, this function is probably most
used repeatedly with a fixed known block, so it is more efficient to require
that such a block be passed as a parameter. Much of the effort here goes into
testing that the parameters supplied are coherent with each other; after this,
the method |Block::element| does the actual work of looking up the pair
$(x,y)$ in the block.

@< Local function def...@>=
void block_index_wrapper(expression_base::level l)
{ shared_KGB_elt y = get<KGB_elt_value>();
  shared_KGB_elt x = get<KGB_elt_value>();
  shared_Block b = get<Block_value>();
  if (&b->rf->parent.val!=&x->rf->parent.val)
    throw runtime_error ("Real form not in inner class of block");
  if (&b->rf->parent.val!=&y->rf->parent.dual)
    throw runtime_error ("Dual real form not in inner class of block");
  const KGB& kgb = b->rf->kgb(); const KGB& dual_kgb = b->dual_rf->kgb();
  const TwistedWeylGroup& tw = kgb.twistedWeylGroup();
  const TwistedWeylGroup& dual_tw = dual_kgb.twistedWeylGroup();
  if (blocks::dual_involution(kgb.involution(x->val),tw,dual_tw)
      != dual_kgb.involution(y->val))
    throw runtime_error ("Fiber mismatch KGB and dual KGB elements");

  BlockElt z = b->val.element(x->val,y->val);

  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(z));
}

@ The dual block might be computed from other functions (provided the block
has any elements), but it is very easy to implement directly.

@< Local function def...@>=
void dual_block_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<Block_value>(b->dual_rf,b->rf));
}

@ To make blocks more easily useful we add functions giving a status code for
a Weyl group generator with respect to a block element.

@< Local function def...@>=
void block_status_wrapper(expression_base::level l)
{ BlockElt i = get<int_value>()->int_val();
  shared_Block b = get<Block_value>();
  unsigned int s = get<int_value>()->int_val();
  if (s>=b->rf->val.semisimpleRank())
    throw runtime_error ("Illegal simple reflection: "+str(s));
  if (i>=b->val.size())
    throw runtime_error
      ("Block element " +str(i) + " out of range (<" + str(b->val.size())+")");
  if (l==expression_base::no_value)
    return;
@)
  const DescentStatus::Value dv = b->val.descentValue(s,i);
@/// renumber from |DescentStatus::Value| order to C-,ic,r1,r2,C+,rn,i1,i2
  static const unsigned char tab [] = {4,5,6,7,1,0,3,2};
  push_value(std::make_shared<int_value>(tab[dv]));
}

@ We also allow computing cross actions and (inverse) Cayley transforms.

@< Local function def...@>=

void block_cross_wrapper(expression_base::level l)
{ BlockElt i = get<int_value>()->int_val();
  shared_Block b = get<Block_value>();
  unsigned int s = get<int_value>()->int_val();
  if (s>=b->rf->val.semisimpleRank())
    throw runtime_error ("Illegal simple reflection: "+str(s));
  if (i>=b->val.size())
    throw runtime_error
      ("Block element " +str(i) + " out of range (<" + str(b->val.size())+")");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(b->val.cross(s,i)));
}
@)
void block_Cayley_wrapper(expression_base::level l)
{ shared_int i = get<int_value>();
  unsigned int ii = i->int_val();
  shared_Block b = get<Block_value>();
  unsigned int s = get<int_value>()->int_val();
  if (s>=b->rf->val.semisimpleRank())
    throw runtime_error ("Illegal simple reflection: "+str(s));
  if (ii >= b->val.size())
    throw runtime_error
      ("Block element " +str(ii) + " out of range (<" + str(b->val.size())+")");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt sx = b->val.cayley(s,ii).first;
  if (sx==UndefBlock) // when undefined, return i to indicate so
    push_value(i);
  else
    push_value(std::make_shared<int_value>(sx));
}
@)
void block_inverse_Cayley_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  unsigned int ii = i->int_val();
  shared_Block b = get<Block_value>();
  unsigned int s = get<int_value>()->int_val();
  if (s>=b->rf->val.semisimpleRank())
    throw runtime_error ("Illegal simple reflection: "+str(s));
  if (ii >= b->val.size())
    throw runtime_error
      ("Block element " +str(ii) + " out of range (<" + str(b->val.size())+")");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt sx = b->val.inverseCayley(s,ii).first;
  if (sx==UndefBlock) // when undefined, return i to indicate so
    push_value(i);
  else
    push_value(std::make_shared<int_value>(sx));
}

@ Finally we install everything so far related to blocks (some functions
involving Kazhdan-Lusztig polynomials will be added later).

@< Install wrapper functions @>=
install_function(Fokko_block_wrapper,"block","(RealForm,RealForm->Block)");
install_function(decompose_block_wrapper,@|"%","(Block->RealForm,RealForm)");
install_function(size_of_block_wrapper,"#","(Block->int)");
install_function(block_element_wrapper,"element","(Block,int->KGBElt,KGBElt)");
install_function(block_index_wrapper,"index","(Block,KGBElt,KGBElt->int)");
install_function(dual_block_wrapper,"dual","(Block->Block)");
install_function(block_status_wrapper,"status","(int,Block,int->int)");
install_function(block_cross_wrapper,@|"cross","(int,Block,int->int)");
install_function(block_Cayley_wrapper,@|"Cayley","(int,Block,int->int)");
install_function(block_inverse_Cayley_wrapper,@|"inverse_Cayley"
                ,"(int,Block,int->int)");

@*1 Standard module parameters.
%
We implement a data type for holding parameters that represent standard
modules. Such a parameter is defined by a triple $(x,\lambda,\nu)$ where $x$
is a KGB element, $\lambda$ is a weight in the coset $\rho+X^*$ whose value is
relevant only modulo the sub-lattice $(1-\theta_x)X^*$ where $\theta_x$ is the
involution associated to $x$, and $\nu$ is a rational weight in the kernel of
$1+\theta_x$. Such parameters are stored in instances of the class
|StandardRepr|, which is defined in the file \.{repr.h}.

@<Includes needed in the header file @>=
#include "repr.h"

@*2 Class definition.
Like for KGB elements, we maintain a shared pointer to the real form value, so
that it will be assured to survive as long as parameters for it exist, and we
can access notably the |Rep_context| that it provides.

@< Type definitions @>=
struct module_parameter_value : public value_base
{ own_real_form rf;
  StandardRepr val;
@)
  module_parameter_value(const own_real_form& form, const StandardRepr& v)
  : rf(form), val(v) @+{}
  ~module_parameter_value() @+{}
@)
  virtual void print(std::ostream& out) const;
  module_parameter_value* clone() const
   @+ {@; return new module_parameter_value(*this); }
  static const char* name() @+{@; return "module parameter"; }
@)
  const Rep_context& rc() const @+{@; return rf->rc(); }
  Rep_table& rt() const @+{@; return rf->rt(); }
private:
  module_parameter_value(const module_parameter_value& v)  // copy constructor
  : rf(v.rf),val(v.val) @+ {}
};
@)
typedef std::unique_ptr<module_parameter_value> module_parameter_ptr;
typedef std::shared_ptr<const module_parameter_value> shared_module_parameter;
typedef std::shared_ptr<module_parameter_value> own_module_parameter;

@ When printing a module parameter, we shall indicate a triple
$(x,\lambda,\nu)$ that defines it. Since we shall need to print |StandardRepr|
values in other contexts as well, we shall define an auxiliary output function
|print_stdrep| for such values first, which takes and additional |Rep_context|
argument, and then call that function from |module_parameter_value::print|. By
choosing a name for the auxiliary function different from |print|, we avoid
having that call being mistaken for a recursive call.

@f nu nullptr

@< Local function def...@>=
std::ostream& print_stdrep
  (std::ostream& out,const StandardRepr& val, const Rep_context& rc)
{ return @|
  out << "parameter(x="
      << val.x() << ",lambda="
      << rc.lambda(val) << ",nu="
      << rc.nu(val) << ')';
}

@ Here is virtual method |module_parameter_value::print|, used when printing a
value of type \.{Param} (as opposed to for instance printing a term of
a \.{ParamPol}, which calls |print_stdrep|). Here we prefix the parameter text
proper with additional information about the parameter that may be relevant to
the user.

@< Function definition... @>=
void module_parameter_value::print(std::ostream& out) const
{ RootNbr witness; // dummy needed in call
  out << @< Expression for adjectives that apply to a module parameter @>@;@;;
  print_stdrep(out << ' ',val,rc());
}

@ We provide one of the adjectives ``non-standard'' (when $\gamma$ and
therefore $\lambda$ fails to be imaginary-dominant), ``zero'' (the standard
module vanishes due to the singular infinitesimal character, namely by the
presence of a singular compact simple-imaginary root), ``non-final'' (the
standard module is non-zero, but can be expressed in terms of standard modules
at more compact Cartans using a singular real root satisfying the parity
condition), ``non-normal'' (the parameter differs from its normal form; when
we come to this point it implies there is a complex singular descent), or
finally ``final'' (the good ones that could go into a \.{ParamPol} value; the
condition |is_final| should apply, though it is not tested here).

@< Expression for adjectives... @>=
( not rc().is_standard(val,witness) ? "non-standard"
@|: not rc().is_dominant(val,witness) ? "non-dominant"
@|: not rc().is_nonzero(val,witness) ? "zero"
@|: not rc().is_semifinal(val,witness) ? "non-final"
@|: not rc().is_normal(val) ? "non-normal"
@|: "final")

@ To make a module parameter, one should provide a KGB element~$x$, an
integral weight $\lambda-\rho$, and a rational weight~$\nu$. Since only its
projection on the $-\theta_x$-stable subspace is used, one might specify the
infinitesimal character $\gamma$ in the place of $\nu$.

@< Local function def...@>=
void module_parameter_wrapper(expression_base::level l)
{ shared_rational_vector nu(get<rational_vector_value>());
  shared_vector lam_rho(get<vector_value>());
  shared_KGB_elt x = get<KGB_elt_value>();
  if (nu->val.size()!=lam_rho->val.size()
      or nu->val.size()!=x->rf->val.rank())
    throw runtime_error ("Rank mismatch: ("
        +str(x->rf->val.rank())+","
	+str(lam_rho->val.size())+","+str(nu->val.size())+")");
@.Rank mismatch@>
  if (l!=expression_base::no_value)
    push_value(std::make_shared<module_parameter_value> @| (x->rf,
      x->rf->rc().sr(x->val,lam_rho->val,nu->val)));
}

@*2 Functions operating on module parameters.
%
The following function, which we shall bind to the monadic operator `|%|',
transforms a parameter value into a triple of values $(x,\lambda-\rho,\gamma)$
that defines it, where $\gamma$ is taken to be (a representative of) the
infinitesimal character. (This function used to produce $\nu$ as third
component, but in practice obtaining the infinitesimal character directly
turned out to often be more useful. If needed $\nu$ is easily computed as
$\nu={1+\theta_x\over2}\gamma$; as the opposite conversion requires more work,
we had a separate function returning$~\gamma$ which is now superfluous.) The
triple returned here is not unique, since $\lambda$ is determined only modulo
$(1-\theta_x)X^*$, and this function should make a unique choice. But that
fact (and the choice made) is hidden in the implementation of |StandardRepr|,
of which we just call methods here.

@< Local function def...@>=
void unwrap_parameter_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<KGB_elt_value>(p->rf,p->val.x()));
  push_value(std::make_shared<vector_value>(p->rc().lambda_rho(p->val)));
  push_value(std::make_shared<rational_vector_value>(p->val.gamma()));
  if (l==expression_base::single_value)
    wrap_tuple<3>();
}

@)

void real_form_of_parameter_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(p->rf);
}

@ Here are some more attributes, in the form of predicates.

@< Local function def...@>=
void is_standard_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  RootNbr witness;
  if (l!=expression_base::no_value)
    push_value(whether(p->rc().is_standard(p->val,witness)));
}

void is_zero_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  RootNbr witness;
  if (l!=expression_base::no_value)
    push_value(whether(not p->rc().is_nonzero(p->val,witness)));
}

void is_semifinal_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  RootNbr witness;
  if (l!=expression_base::no_value)
    push_value(whether(p->rc().is_semifinal(p->val,witness)));
}

void is_final_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(whether(p->rc().is_final(p->val)));
}

@ Before constructing (non-integral) blocks, it is essential that the
infinitesimal character is made dominant, and that possible singular (simple)
complex descents are applied to the parameter, as the result of the block
construction will only be mathematically meaningful under these circumstances.
This operation is therefore applied automatically in several places, but it is
useful to give the user an easy way to apply it explicitly. (In fact integral
dominance should be sufficient: only for those coroots with integral evaluation
on~$\gamma$ should it be required that the evaluation be non-negative. But
since making this concrete requires an overhaul of the entire block
construction process, we stick to unqualified dominance for now.)

In accordance with their behaviour when incorporating virtual modules, the
equality operator for parameters will test for \emph{equivalence} when both
are standard; otherwise it tests strict equality. Equivalence of standard
parameters amounts to testing for equality after the parameters are made
dominant (at least that claim was not contested at the time of writing this).
We provide this test, which will be bound to the equality operator. Unlike
earlier equality tests, we \emph{require} the parameters to be associated to
the same real form, giving a runtime error (rather than returning false) if
not; this avoids confusion if there were some subtle difference of real forms
for otherwise similar parameters. If some operation is used to produce
parameters that may of may not be associated to the same real form, then one
should test those forms for equality before testing the parameters.

@< Local function def...@>=
void parameter_dominant_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  if (l!=expression_base::no_value)
  {@; p->rc().make_dominant(p->val);
    push_value(p);
  }
}

void parameter_normal_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  if (l!=expression_base::no_value)
  {@; p->rc().normalise(p->val);
    push_value(p);
  }
}
@)
void parameter_eq_wrapper(expression_base::level l)
{ shared_module_parameter q = get<module_parameter_value>();
  shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(whether(p->rf->val==q->rf->val and p->val==q->val));
}
void parameter_neq_wrapper(expression_base::level l)
{ shared_module_parameter q = get<module_parameter_value>();
  shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(whether(p->rf->val!=q->rf->val or p->val!=q->val));
}

void parameter_equivalent_wrapper(expression_base::level l)
{ shared_module_parameter q = get<module_parameter_value>();
  shared_module_parameter p = get<module_parameter_value>();
  if (p->rf->val!=q->rf->val)
    throw runtime_error @|
      ("Real form mismatch when testing equivalence");
  if (l!=expression_base::no_value)
    push_value(whether(p->rc().equivalent(p->val,q->val)));
}

@ While parameters can be used to compute blocks of (other) parameters, it can
be useful to have available the basic operations of cross actions and Cayley
transforms on individual parameters without going through the construction of
an entire block. This should basically be simple because the construction of
blocks is based on just such operations defined on individual parameters;
however the reality is more complicated because the storage format
|StandardRepr| used in parameter values is not directly suited to the way
(currently) cross actions and Cayley transforms are computed. So the functions
below involve the actual computation sandwiched between unpacking and
repacking operations; this is hidden in the methods |Rep_context::cross| and
friends that are called blow.

Like for KGB elements there is the possibilty of double values, this time both
for the Cayley and inverse Cayley transforms. The ``solution'' to this
difficulty is the same here: the user can find out by herself about a possible
second image by applying a cross action to the result. In the current case
this approach has in fact already been adopted in the methods that are called
here, which present a single-minded interface to these transforms.

@< Local function def...@>=
void parameter_cross_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  int s = get<int_value>()->int_val();
  unsigned int r =
    rootdata::integrality_rank(p->rf->val.rootDatum(),p->val.gamma());
  if (static_cast<unsigned>(s)>=r)
    throw runtime_error
      ("Illegal simple reflection: "+str(s)+ ", should be <"+str(r));
  if (l!=expression_base::no_value)
    push_value(std::make_shared<module_parameter_value>
		(p->rf,p->rc().cross(s,p->val)));
}
@)
void parameter_Cayley_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  int s = get<int_value>()->int_val();
  unsigned int r =
    rootdata::integrality_rank(p->rf->val.rootDatum(),p->val.gamma());
  if (static_cast<unsigned>(s)>=r)
    throw runtime_error
      ("Illegal simple reflection: "+str(s)+ ", should be <"+str(r));
  if (l!=expression_base::no_value)
    push_value(std::make_shared<module_parameter_value>
		(p->rf,p->rc().Cayley(s,p->val)));
}

void parameter_inv_Cayley_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  int s = get<int_value>()->int_val();
  unsigned int r =
    rootdata::integrality_rank(p->rf->val.rootDatum(),p->val.gamma());
  if (static_cast<unsigned>(s)>=r)
    throw runtime_error
      ("Illegal simple reflection: "+str(s)+ ", should be <"+str(r));
  if (l!=expression_base::no_value)
    push_value(std::make_shared<module_parameter_value>
		(p->rf,p->rc().inv_Cayley(s,p->val)));
}

@ The above (old) functions emulate the built-in non-integral block
construction, and in doing so describe the root argument by index \emph{in the
integral subsystem} which is potentially confusing, especially in the case
where the infinitesimal character of the parameter is not dominant, since the
integral subsystem is only fixed after transforming the parameter into an
equivalent one with dominant infinitesimal character. The below function
implement a new approach, less confusing and more general, in which the root
is specified in coordinates. Given the internal implementation, there is in
fact little dealing with the integral subsystem at all, just a test that it
contains the specified root (if not an error is thrown inside the method
called; the error will not happen if |l==expression_base::no_value|, which is
not quite correct).

There is not much difference between integral roots that are simple for the
subsystem and other roots, since what counts in the implementation of Cayley
transforms is being simple for the whole system, and some conjugation is
necessary to achieve that anyway. We also blur the distinction between forward
and inverse Cayley transforms here, as the code will choose what to do
(possibly nothing in undefined cases) depending on the arguments provided.
That code will throw an |error::Cayley_error| in case the Cayley transform is
undefined, but we catch that error here, and return the parameter unchanged;
as long as users have no means to catch errors, this gives them the option to
handle undefined Cayley transforms by testing the returned value.

@h "error.h"
@< Local function def...@>=

void root_parameter_cross_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  shared_vector alpha = get<vector_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<module_parameter_value>
		(p->rf,p->rc().cross(alpha->val,p->val)));
}
@)
void root_parameter_Cayley_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  shared_vector alpha = get<vector_value>();
  if (l==expression_base::no_value)
    return;
@)
  try {
    push_value(std::make_shared<module_parameter_value>
		(p->rf,p->rc().any_Cayley(alpha->val,p->val)));
  }
  catch (error::Cayley_error& e) // ignore undefined Cayley transforms
  {@;
    push_value(p);
  }
}


@ One useful thing to be able to for parameters is to compute their twist by
the distinguished involution.
@< Local function def...@>=
void parameter_twist_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<module_parameter_value>
		(p->rf,p->rc().inner_twisted(p->val)));
}

@ The same can be done for an external twist, specified by a matrix. The
method of |Rep_context| to call is now |twisted| rather than |inner_twisted|,
but we first call |test_compatible|.

@< Local function def...@>=
void parameter_outer_twist_wrapper(expression_base::level l)
{ auto delta = get<matrix_value>();
  auto p = get<module_parameter_value>();
  test_compatible(p->rc().innerClass(),delta);
  if (l!=expression_base::no_value)
    push_value(std::make_shared<module_parameter_value>
		(p->rf,p->rc().twisted(p->val,delta->val)));
}


@ The library can also compute orientation numbers for parameters.

@< Local function def...@>=
void orientation_number_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(p->rc().orientation_number(p->val)));
}

@ Here is a function that computes a list of positive rational values $t\leq1$
such that the parameter obtained by replacing the continuous part~$\nu$ of
by~$t\nu$ is not topmost in its block, so that deformation of the parameter to
$t\nu$ will produce a non-trivial decomposition.

@< Local function def...@>=
void reducibility_points_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l==expression_base::no_value)
    return;
@)
  RationalList rp = p->rc().reducibility_points(p->val);
      // method normalises rationals
  own_row result = std::make_shared<row_value>(rp.size());
  for (size_t i=0; i<rp.size(); ++i)
    result->val[i]=std::make_shared<rat_value>(rp[i]);
  push_value(std::move(result));
}

@ Scaling the continuous component~$\nu$ of a parameter is an important
ingredient for calculating signatures of Hermitian forms. This was for a long
time done with a user defined function, but having this built in is more
efficient. Moreover scaling by~$0$ (called ``deformation to $\nu=0$'') can be
done even more efficiently by specialised code.

@< Local function def...@>=
void scale_parameter_wrapper(expression_base::level l)
{ shared_rat f = get<rat_value>();
  own_module_parameter p = get_own<module_parameter_value>();
  if (l!=expression_base::no_value)
@/{@; p->rc().scale(p->val,f->rat_val());
    push_value(std::move(p));
  }
}

void scale_0_parameter_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  if (l!=expression_base::no_value)
@/{@; p->rc().scale_0(p->val);
    push_value(std::move(p));
  }
}

@ One of the main reasons to introduce module parameter values is that they
allow computing a block, whose elements are again given by module parameters.
Before we define the functions the do that, let us define a common function
|test_standard| they will use to test a parameter for validity. Even though we
shall always have a pointer available when we call |test_standard|, we define
this function to take a reference (requiring us to write a dereferencing at
each call), because the type of pointer (shared or raw) available is not
always the same. The reference is of course not owned by |test_standard|. A
similar test is |test_normal_is_final|, which tests if the parameter will,
after applying |Rep_table::normalise|, satisfy the |is_final| predicate.
When these functions fail, they try to be specific about what condition fails,
in terms of the situation before applying |normalise|

@< Local function def...@>=
void test_standard(const module_parameter_value& p, const char* descr)
{ RootNbr witness;
  if (p.rc().is_standard(p.val,witness))
    return;
  std::ostringstream os; p.print(os << descr << ":\n  ");
  os << "\n  Parameter not standard, negative on coroot #" << witness;
  throw runtime_error(os.str());
}

void test_normal_is_final(const module_parameter_value& p, const char* descr)
{ RootNbr witness; bool nonzero=p.rc().is_nonzero(p.val,witness);
  if (nonzero and p.rc().is_semifinal(p.val,witness))
    return; // nothing to report
  std::ostringstream os; p.print(os << descr << ":\n  ");
@/os << "\n  Parameter is " << (nonzero ? "not semifinal" : "zero")
   @|  <<", as witnessed by coroot #" << witness;
  throw runtime_error(os.str());
}

@ Here is the first block generating function, which just reproduces to output
from the \.{Fokko} program for the \.{nblock} command.

@< Local function def...@>=
void print_n_block_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  BlockElt init_index; // will hold index in the block of the initial element
  param_block block(p->rc(),p->val,init_index);
  *output_stream << "Parameter defines element " << init_index
               @|<< " of the following block:" << std::endl;
  block.print_to(*output_stream,true);
    // print block using involution expressions
  if (l==expression_base::single_value)
    wrap_tuple<0>(); // |no_value| needs no special care
}

@ More interesting than printing the block is to return is to the user as a
list of parameter values. The following function does this, and adds as a
second result the index that the original parameter has in the resulting
block.

@< Local function def...@>=
void block_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt start; // will hold index in the block of the initial element
  param_block block(p->rc(),p->val,start);
  @< Push a list of parameter values for the elements of |block| @>
  push_value(std::make_shared<int_value>(start));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ Construction a list of values is a routine affair. This code must however
also construct a module parameter value for each element of |block|.

@< Push a list of parameter values for the elements of |block| @>=
{ own_row param_list = std::make_shared<row_value>(block.size());
  for (BlockElt z=0; z<block.size(); ++z)
    param_list->val[z] =
	std::make_shared<module_parameter_value>(p->rf,block.sr(z));
  push_value(std::move(param_list));

}

@ There is also a function that computes just a partial block.
@< Local function def...@>=
void partial_block_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  if (l==expression_base::no_value)
    return;
@)
  param_block block(p->rc(),p->val);
  @< Push a list of parameter values for the elements of |block| @>
}

@ Knowing the length in its block of a parameter is of independent interest.
@< Local function def...@>=
void param_length_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot determine block for parameter length");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(p->rt().length(p->val)));
}

@ Here is a version of the |block| command that also exports the table of
Kazhdan-Lusztig polynomials for the block, in the same format as \\{raw\_KL}
that will be defined below.

@< Local function def...@>=
void KL_block_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt start; // will hold index in the block of the initial element
  param_block block(p->rc(),p->val,start);
  @< Push a list of parameter values for the elements of |block| @>
  push_value(std::make_shared<int_value>(start));
  const kl::KLContext& klc = block.klc(block.size()-1,false);
@)
  @< Extract from |klc| an |own_matrix M@;| and |own_row polys@;| @>
@)
  own_vector length_stops = std::make_shared<vector_value>(
     int_Vector(block.length(block.size()-1)+2));
  for (size_t i=0; i<length_stops->val.size(); ++i)
    length_stops->val[i]=block.length_first(i);
@)
  unsigned n_survivors=0;
  for (BlockElt z=0; z<block.size(); ++z)
  @+ if (block.survives(z))
      ++n_survivors;
  own_vector survivor =
    std::make_shared<vector_value>(int_Vector(n_survivors));
  { unsigned i=0;
    for (BlockElt z=0; z<block.size(); ++z)
    @+if (block.survives(z))
        survivor->val[i++]=z;
    assert(i==n_survivors);
  }
  own_matrix contributes_to = std::make_shared<matrix_value>(
    int_Matrix(n_survivors,block.size(),0));
  for (BlockElt z=0; z<block.size(); ++z)
  { BlockEltList sb = block.finals_for(z);
    for (BlockEltList::const_iterator it=sb.begin(); it!=sb.end(); ++it)
    { BlockElt x= permutations::find_index<int>(survivor->val,*it);
        // a row index
      if ((block.length(z)-block.length(*it))%2==0)
        ++contributes_to->val(x,z);
      else
        --contributes_to->val(x,z);
    }
  }
@)
  push_value(std::move(M));
  push_value(std::move(polys));
  push_value(std::move(length_stops));
  push_value(std::move(survivor));
  push_value(std::move(contributes_to));

  if (l==expression_base::single_value)
    wrap_tuple<7>();
}

@ The following module should not be enclosed in braces, as it defines two
variable |M| and |polys|. One reason to extract it is that it can be used
identically in two wrapper functions.

@< Extract from |klc| an |own_matrix M@;| and |own_row polys@;| @>=
own_matrix M = std::make_shared<matrix_value>(int_Matrix(klc.size()));
for (size_t y=1; y<klc.size(); ++y)
  for (size_t x=0; x<y; ++x)
    M->val(x,y)= klc.KL_pol_index(x,y);
@)
own_row polys = std::make_shared<row_value>(0);
polys->val.reserve(klc.polStore().size());
for (size_t i=0; i<klc.polStore().size(); ++i)
{
  const kl::KLPol& pol = klc.polStore()[i];
  std::vector<int> coeffs(pol.begin(),pol.end());
  polys->val.emplace_back(std::make_shared<vector_value>(coeffs));
}

@ Here is a dual variation of the previous function. The main difference is
calling the pseudo constructor |blocks::Bare_block::dual| to transform |block|
into its dual (represented as just a |blocks::Bare_block| which is sufficient)
before invoking the KL computations. The block is reversed with respect to
|block|, so for proper interpretation we reverse the list of
parameters returned, and this that several other result components have to be
transformed as well. On the dual side there should be no condensing of the
polynomial matrix on the ``survivor'' elements, rather just an extraction of a
submatrix of polynomials at the corresponding indices; therefore we leave out
the |contributes_to| matrix altogether.

@< Local function def...@>=
void dual_KL_block_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt start; // will hold index in the block of the initial element
  param_block block(p->rc(),p->val,start);
  auto size1 = block.size()-1;
  @< Push a reversed list of parameter values for the elements of |block| @>
  push_value(std::make_shared<int_value>(size1-start));
  auto dual_block = blocks::Bare_block::dual(block);
  const kl::KLContext& klc = dual_block.klc(block.size()-1,false);
@)
  @< Extract from |klc| an |own_matrix M@;| and |own_row polys@;| @>
@)
  own_vector length_stops = std::make_shared<vector_value>(
     int_Vector(block.length(size1)+2));
  for (size_t i=0; i<length_stops->val.size(); ++i)
    // subtract from |block.size()| here:
    length_stops->val[i] =
      block.size()-block.length_first(length_stops->val.size()-1-i);
@)
  unsigned n_survivors=0;
  for (BlockElt z=0; z<block.size(); ++z)
  @+if (block.survives(z))
      ++n_survivors;
  own_vector survivor =
    std::make_shared<vector_value>(int_Vector(n_survivors));
  { unsigned i=0;
    for (BlockElt z=block.size(); z-->0; )
    @+if (block.survives(z))
        survivor->val[i++]=size1-z;
    assert(i==n_survivors);
  }
@)
  push_value(std::move(M));
  push_value(std::move(polys));
  push_value(std::move(length_stops));
  push_value(std::move(survivor));

  if (l==expression_base::single_value)
    wrap_tuple<6>();
}

@ Reversing a list after it has been pushed to the stack is somewhat
cumbersome, so we prefer to redo a previous module with a slight modification
to reverse the order.

@< Push a reversed list of parameter values for the elements of |block| @>=
{ own_row param_list = std::make_shared<row_value>(0);
  param_list->val.reserve(block.size());
  for (BlockElt z=block.size(); z-->0;)
    param_list->val.push_back @|
	(std::make_shared<module_parameter_value>(p->rf,block.sr(z)));
  push_value(std::move(param_list));

}

@ Here is a version of the |KL_block| that computes just for a partial block
and the Kazhdan-Lusztig polynomials for it. There are six components
in the value returned: the list of parameters forming the partial block (of
which the final one is the initial parameter), a matrix of KL-polynomial
indices, a list of polynomials (as vectors), a vector of length stops (block
element numbers at with the length function increases), a list of block
element numbers for those whose survive the translation-to-singular functor,
and a matrix that indicates which block element contributes to which surviving
element (and with what multiplicity).

@< Local function def...@>=
void partial_KL_block_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  if (l==expression_base::no_value)
    return;
@)
  param_block block(p->rc(),p->val);
  @< Push a list of parameter values for the elements of |block| @>

  const kl::KLContext& klc = block.klc(block.size()-1,false);
  // compute KL polynomials, silently

  own_matrix M = std::make_shared<matrix_value>(int_Matrix(klc.size()));
  for (size_t y=1; y<klc.size(); ++y)
    for (size_t x=0; x<y; ++x)
      M->val(x,y)= klc.KL_pol_index(x,y);
@)
  own_row polys = std::make_shared<row_value>(0);
  polys->val.reserve(klc.polStore().size());
  for (size_t i=0; i<klc.polStore().size(); ++i)
  {
    const kl::KLPol& pol = klc.polStore()[i];
    std::vector<int> coeffs(pol.size());
    for (size_t j=pol.size(); j-->0; )
      coeffs[j]=pol[j];
    polys->val.emplace_back(std::make_shared<vector_value>(coeffs));
  }
@)
  own_vector length_stops = std::make_shared<vector_value>(
     int_Vector(block.length(block.size()-1)+2));
  length_stops->val[0]=0;
  for (size_t i=1; i<length_stops->val.size(); ++i)
    length_stops->val[i]=block.length_first(i);
@)
  unsigned n_survivors=0;
  for (BlockElt z=0; z<block.size(); ++z)
    if (block.survives(z))
      ++n_survivors;
  own_vector survivor =
    std::make_shared<vector_value>(int_Vector(n_survivors));
  { unsigned i=0;
    for (BlockElt z=0; z<block.size(); ++z)
      if (block.survives(z))
        survivor->val[i++]=z;
    assert(i==n_survivors);
  }
  own_matrix contributes_to = std::make_shared<matrix_value>(
    int_Matrix(n_survivors,block.size(),0));
  for (BlockElt z=0; z<block.size(); ++z)
  { BlockEltList sb = block.finals_for(z);
    for (BlockEltList::const_iterator it=sb.begin(); it!=sb.end(); ++it)
    { BlockElt x= permutations::find_index<int>(survivor->val,*it);
        // a row index
      if ((block.length(z)-block.length(*it))%2==0)
        ++contributes_to->val(x,z);
      else
        --contributes_to->val(x,z);
    }
  }
@)
  push_value(std::move(M));
  push_value(std::move(polys));
  push_value(std::move(length_stops));
  push_value(std::move(survivor));
  push_value(std::move(contributes_to));

  if (l==expression_base::single_value)
    wrap_tuple<6>();
}

@ The function |extended_block| intends to make computation of extended
blocks available in \.{atlas}.

@< Local function def...@>=
void extended_block_wrapper(expression_base::level l)
{ auto delta =get<matrix_value>();
  shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  test_compatible(p->rc().innerClass(),delta);
  if (not ((delta->val-1)*p->val.gamma().numerator()).isZero())
    throw runtime_error("Involution does not fix infinitesimal character");
  if (l==expression_base::no_value)
    return;
@)
  const auto& rc = p->rc();
  BlockElt start;
  param_block block(rc,p->val,start);
  @< Construct the extended block, then the return value components,
     calling |push_value| for each of them @>
@)
  if (l==expression_base::single_value)
    wrap_tuple<4>();
}

@ We somewhat laboriously convert internal information from the extended block
into a list of parameters and three tables in the form of matrices.

@< Construct the extended block... @>=
{ ext_block::ext_block eb(rc.innerClass(),block,delta->val);
  own_row params = std::make_shared<row_value>(eb.size());
  int_Matrix types(eb.size(),eb.rank());
@/int_Matrix links0(eb.size(),eb.rank());
  int_Matrix links1(eb.size(),eb.rank());

  for (BlockElt n=0; n<eb.size(); ++n)
  { auto z = eb.z(n); // number of ordinary parameter in |block|
    StandardRepr block_elt_param =
      rc.sr_gamma(block.x(z),block.lambda_rho(z),block.gamma());
    params->val[n] =
      std::make_shared<module_parameter_value>(p->rf,block_elt_param);
    for (weyl::Generator s=0; s<eb.rank(); ++s)
    { auto type = eb.descent_type(s,n);
      types(n,s) = static_cast<int>(type);
      if (is_like_compact(type) or is_like_nonparity(type))
      @/{@; links0(n,s)=eb.size(); links1(n,s)=eb.size(); }
      else
      { links0(n,s)= is_complex(type) ? eb.cross(s,n): eb.Cayley(s,n);
        if (eb.epsilon(s,n,links0(n,s))<0)
	  links0(n,s) = -1-links0(n,s);
        if (link_count(type)==1)
          links1(n,s)=eb.size(); // leave second matrix entry empty
        else
        { links1(n,s)= has_double_image(type)
            ? eb.Cayleys(s,n).second
            : eb.cross(s,n);
          if (eb.epsilon(s,n,links1(n,s))<0)
	    links1(n,s) = -1-links1(n,s);
        }
      }
    }
  }

  push_value(std::move(params));
  push_value(std::make_shared<matrix_value> (std::move(types)));
  push_value(std::make_shared<matrix_value>(std::move(links0)));
  push_value(std::make_shared<matrix_value>(std::move(links1)));
}

@ The following function generates an extended block, computes their extended
KL polynomials and evaluates them at $-1$ (since this turns out to be
sufficient for their use in the deformation algorithm), and then rewrites
elements that have singular descents in terms of those that have not (the
``survivors''), reduces the matrix to be indexed by those elements only, and
finally negates entries at positions with odd length difference for the block
elements corresponding to row and column. All this work is actually performed
inside call to |ext_kl::ext_KL_matrix|.

The function returns the extended block as list of parameters, the matrix just
described, and a list of element lengths.

@< Local function def...@>=
void extended_KL_block_wrapper(expression_base::level l)
{ auto delta = get<matrix_value>();
  auto p = get<module_parameter_value>();
  test_standard(*p,"Cannot generate extended block");
  test_compatible(p->rc().innerClass(),delta);
  if (l==expression_base::no_value)
    return;
@)
  std::vector<StandardRepr> block;
  own_matrix P_mat = std::make_shared<matrix_value>(int_Matrix());
  own_vector lengths = std::make_shared<vector_value>(int_Vector());
  ext_kl::ext_KL_matrix
	(p->val,delta->val,p->rc(),block,P_mat->val,lengths->val);
@)
  own_row param_list = std::make_shared<row_value>(block.size());
  for (BlockElt z=0; z<block.size(); ++z)
    param_list->val[z]=std::make_shared<module_parameter_value>(p->rf,block[z]);
  push_value(std::move(param_list));
  push_value(std::move(P_mat));
  push_value(std::move(lengths));
  if (l==expression_base::single_value)
    wrap_tuple<3>();
}

@ Finally we install everything related to module parameters.
@< Install wrapper functions @>=
install_function(module_parameter_wrapper,@|"param"
                ,"(KGBElt,vec,ratvec->Param)");
install_function(unwrap_parameter_wrapper,@|"%"
                ,"(Param->KGBElt,vec,ratvec)");
install_function(real_form_of_parameter_wrapper,@|"real_form"
		,"(Param->RealForm)");
install_function(is_standard_wrapper,@|"is_standard" ,"(Param->bool)");
install_function(is_zero_wrapper,@|"is_zero" ,"(Param->bool)");
install_function(is_semifinal_wrapper,@|"is_semifinal" ,"(Param->bool)");
install_function(is_final_wrapper,@|"is_final" ,"(Param->bool)");
install_function(parameter_dominant_wrapper,@|"dominant" ,"(Param->Param)");
install_function(parameter_normal_wrapper,@|"normal" ,"(Param->Param)");
install_function(parameter_eq_wrapper,@|"=", "(Param,Param->bool)");
install_function(parameter_neq_wrapper,@|"!=", "(Param,Param->bool)");
install_function(parameter_equivalent_wrapper,@|"equivalent"
                ,"(Param,Param->bool)");
install_function(parameter_cross_wrapper,@|"cross" ,"(int,Param->Param)");
install_function(parameter_Cayley_wrapper,@|"Cayley" ,"(int,Param->Param)");
install_function(parameter_inv_Cayley_wrapper,@|"inv_Cayley"
                ,"(int,Param->Param)");
install_function(root_parameter_cross_wrapper,@|"cross" ,"(vec,Param->Param)");
install_function(root_parameter_Cayley_wrapper,@|"Cayley" ,"(vec,Param->Param)");
install_function(parameter_twist_wrapper,@|"twist" ,"(Param->Param)");
install_function(parameter_outer_twist_wrapper,@|"twist" ,"(Param,mat->Param)");
install_function(orientation_number_wrapper,@|"orientation_nr" ,"(Param->int)");
install_function(reducibility_points_wrapper,@|
		"reducibility_points" ,"(Param->[rat])");
install_function(scale_parameter_wrapper,"*", "(Param,rat->Param)");
install_function(scale_0_parameter_wrapper,"at_nu_0", "(Param->Param)");
install_function(print_n_block_wrapper,@|"print_block","(Param->)");
install_function(block_wrapper,@|"block" ,"(Param->[Param],int)");
install_function(partial_block_wrapper,@|"partial_block","(Param->[Param])");
install_function(param_length_wrapper,@|"length","(Param->int)");
install_function(KL_block_wrapper,@|"KL_block"
                ,"(Param->[Param],int,mat,[vec],vec,vec,mat)");
install_function(dual_KL_block_wrapper,@|"dual_KL_block"
                ,"(Param->[Param],int,mat,[vec],vec,vec)");
install_function(partial_KL_block_wrapper,@|"partial_KL_block"
                ,"(Param->[Param],mat,[vec],vec,vec,mat)");
install_function(extended_block_wrapper,@|"extended_block"
                ,"(Param,mat->[Param],mat,mat,mat)");
install_function(extended_KL_block_wrapper,@|"extended_KL_block"
                ,"(Param,mat->[Param],mat,vec)");

@*1 Polynomials formed from parameters.
%
When working with parameters for standard modules, and notably with the
deformation formulas, the need arises to keep track of formal sums of
standard modules with coefficients of a type that allows keeping track of the
signatures of the modules. These coefficients, which we shall call split
integers and give the type \.{Split} in \.{atlas}, are elements of the group
algebra over $\Zee$ of a (cyclic) group of order~$2$.

@< Includes needed in the header file @>=
#include "arithmetic.h"

@*2 A class for split integers.
%
Although the necessary operations could easily be defined in the \.{axis}
programming language using pairs of integers, it is preferable to make them an
Atlas type, since this allows distinguishing them from pairs of integers used
for other purposes, and to provide special output and conversion facilities.

@< Type definitions @>=

struct split_int_value : public value_base
{ Split_integer val;
@)
  explicit split_int_value(Split_integer v) : val(v) @+ {}
  ~split_int_value()@+ {}
  void print(std::ostream& out) const;
  split_int_value* clone() const @+{@; return new split_int_value(*this); }
  static const char* name() @+{@; return "split integer"; }
private:
  split_int_value(const split_int_value& v) : val(v.val) @+{}
};
@)
typedef std::unique_ptr<split_int_value> split_int_ptr;
typedef std::shared_ptr<const split_int_value> shared_split_int;
typedef std::shared_ptr<split_int_value> own_split_int;

@ Like for parameter values, we first define a printing function on the level
of a bare |Split_integer| value, which can be used in situations where the
method |split_int_value::print| cannot.

@< Local function def...@>=
std::ostream& print (std::ostream& out, const Split_integer& val)
{@;
  return out << '(' << val.e()
             << (val.s()<0?'-':'+') << std::abs(val.s()) << "s)";
}
@ Again the virtual method |print| must not be defined in the anonymous
namespace.

@< Function def... @>=

void split_int_value::print(std::ostream& out) const @+
{@; interpreter::print(out,val); }

@ Here are some basic relations and arithmetic operations.

@< Local function definitions @>=

void split_unary_eq_wrapper(expression_base::level l)
{ Split_integer i=get<split_int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(whether(i.is_zero()));
}
void split_unary_neq_wrapper(expression_base::level l)
{ Split_integer i=get<split_int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(whether(not i.is_zero()));
}
@)
void split_eq_wrapper(expression_base::level l)
{ Split_integer j=get<split_int_value>()->val;
  Split_integer i=get<split_int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(whether(i==j));
}
void split_neq_wrapper(expression_base::level l)
{ Split_integer j=get<split_int_value>()->val;
  Split_integer i=get<split_int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(whether(i!=j));
}
@)
void split_plus_wrapper(expression_base::level l)
{ Split_integer j=get<split_int_value>()->val;
  own_split_int i=get_own<split_int_value>();
  if (l!=expression_base::no_value)
  {@; i->val+=j; push_value(i); }
}

void split_minus_wrapper(expression_base::level l)
{ Split_integer j=get<split_int_value>()->val;
  own_split_int i=get_own<split_int_value>();
  if (l!=expression_base::no_value)
  {@; i->val-=j; push_value(i); }
}
@)
void split_unary_minus_wrapper(expression_base::level l)
{ own_split_int i=get_own<split_int_value>();
  if (l!=expression_base::no_value)
  {@; i->val.negate(); push_value(i); }
}

void split_times_wrapper(expression_base::level l)
{ Split_integer j=get<split_int_value>()->val;
  Split_integer i=get<split_int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(std::make_shared<split_int_value>(i*j));
}

@ We also provide implicit conversions from integers or pairs of integers to
split integers, and an explicit operator for converting back to a pair.

@< Local function definitions @>=

void int_to_split_coercion()
{ int a=get<int_value>()->int_val();
@/push_value(std::make_shared<split_int_value>(Split_integer(a)));
}
@)
void pair_to_split_coercion()
{ push_tuple_components();
  int b=get<int_value>()->int_val();
  int a=get<int_value>()->int_val();
  push_value(std::make_shared<split_int_value>(Split_integer(a,b)));
}
@)
void from_split_wrapper(expression_base::level l)
{ Split_integer si = get<split_int_value>()->val;
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<int_value>(si.e()));
  push_value(std::make_shared<int_value>(si.s()));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@*2 Class definition for virtual modules.
%
The library provides a type |repr::SR_poly| in which such sums can be
efficiently maintained. In order to use it we must have seen the header file
for the module \.{free\_abelian} on which the implementation is based. While
that class itself does not have such an invariant, the handling of these
formal sums in \.{atlas} will be such that all terms are ensured to have the
predicate |is_final| true, which ensures a number of desirable properties,
including having a dominant representative $\gamma$ of the infinitesimal
character. Only under such restriction can it be guaranteed that equivalent
terms (which now must actually be equal) will always be combined, and the test
for the sum being zero therefore mathematically correct.

@< Includes needed in the header file @>=
#include "free_abelian.h" // needed to make |repr::SR_poly| a complete type

@~Like for KGB elements, we maintain a shared pointer to the real form value, so
that it will be assured to survive as long as parameters for it exist.

@< Type definitions @>=
struct virtual_module_value : public value_base
{ own_real_form rf;
  repr::SR_poly val;
@)
  virtual_module_value(const own_real_form& form, const repr::SR_poly& v)
  : rf(form), val(v) @+{}
  ~virtual_module_value() @+{}
@)
  virtual void print(std::ostream& out) const;
  virtual_module_value* clone() const
   @+ {@; return new virtual_module_value(*this); }
  static const char* name() @+{@; return "module parameter"; }
@)
  const Rep_context& rc() const @+{@; return rf->rc(); }
  Rep_table& rt() const @+{@; return rf->rt(); }
private:
  virtual_module_value(const virtual_module_value& v)
  @+ : rf(v.rf),val(v.val) @+{} // copy
};
@)
typedef std::unique_ptr<virtual_module_value> virtual_module_ptr;
typedef std::shared_ptr<const virtual_module_value> shared_virtual_module;
typedef std::shared_ptr<virtual_module_value> own_virtual_module;

@ When printing a virtual module value, we traverse the |std::map| that is
hidden in the |Free_Abelian| class template, and print individual terms using
the auxiliary function that was defined above for printing parameter values.
However when either all coefficients are integers or coefficients are integer
multiples of~$s$, then we suppress the component that is always~$0$; this is
particularly useful if polynomials are used to encode $\Zee$-linear
combinations of parameters.

@h <iomanip> // for |std::setw|
@< Function def...@>=
void virtual_module_value::print(std::ostream& out) const
{ if (val.empty())
    {@; out << "Empty sum of standard modules"; return; }
  bool has_one=false, has_s=false;
  for (repr::SR_poly::const_iterator it=val.begin(); it!=val.end(); ++it)
  { if (it->second.e()!=0)
      has_one=true;
    if (it->second.s()!=0)
      has_s=true;
    if (has_one and has_s)
      break;
  }
  assert (has_one or has_s); // otherwise the module would have been empty
  for (repr::SR_poly::const_iterator it=val.begin(); it!=val.end(); ++it)
  { out << '\n';
    if (has_one and has_s)
      interpreter::print(out,it->second); // print coefficient
    else if (has_one)
      out << it->second.e();
    else
      out << it->second.s() << 's';
    print_stdrep(out << '*',it->first,rc()); // print parameter
    out << " [" << it->first.height() << ']';
  }
}

@*2 Functions for virtual modules.
%
To start off a |virtual_module_value|, one usually takes an empty sum, but
one needs to specify a real form to fill the |rf| field. The information
allows us to extract the real form from a virtual module even if it is empty.
We allow testing the number of terms of the sum, and directly testing the sum
to be empty.

Testing two virtual modules for equality is also implemented. This could be
done by subtracting and then testing the result for being zero (empty), but it
is more efficient to just traverse both in parallel and stop once a difference
is found.

@< Local function def...@>=
void virtual_module_wrapper(expression_base::level l)
{ own_real_form rf = non_const_get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<virtual_module_value> @|
      (rf,repr::SR_poly(rf->rc().repr_less())));
}
@)
void real_form_of_virtual_module_wrapper(expression_base::level l)
{ shared_virtual_module m = get<virtual_module_value>();
  if (l!=expression_base::no_value)
    push_value(m->rf);
}
@)
void virtual_module_unary_eq_wrapper(expression_base::level l)
{ shared_virtual_module m = get<virtual_module_value>();
  if (l!=expression_base::no_value)
    push_value(whether(m->val.empty()));
}

void virtual_module_unary_neq_wrapper(expression_base::level l)
{ shared_virtual_module m = get<virtual_module_value>();
  if (l!=expression_base::no_value)
    push_value(whether(m->val.size()>0));
}

@)
void virtual_module_eq_wrapper(expression_base::level l)
{ shared_virtual_module n = get<virtual_module_value>();
  shared_virtual_module m = get<virtual_module_value>();
  if (l==expression_base::no_value)
    return;
@)
  auto mit=m->val.begin(), nit=n->val.begin(); // keep outside loop
  for (; mit!=m->val.end() and nit!=n->val.end(); ++mit,++nit)
    if (mit->first!=nit->first or mit->second!=nit->second)
    @/{@;  push_value(whether(false));
      return; }
  push_value(whether
      (mit==m->val.end() and nit==n->val.end()));
}
void virtual_module_neq_wrapper(expression_base::level l)
{ shared_virtual_module n = get<virtual_module_value>();
  shared_virtual_module m = get<virtual_module_value>();
  if (l==expression_base::no_value)
    return;
@)
  auto mit=m->val.begin(), nit=n->val.begin(); // keep outside loop
  for (; mit!=m->val.end() and nit!=n->val.end(); ++mit,++nit)
    if (mit->first!=nit->first or mit->second!=nit->second)
    @/{@;  push_value(whether(true));
      return; }
  push_value(whether
    (mit!=m->val.end() or nit!=n->val.end()));
}

@ This function must be global, it is declared in the header
file \.{global.h}.

@< Function definitions @>=
void virtual_module_size_wrapper(expression_base::level l)
{ shared_virtual_module m = get<virtual_module_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(m->val.size()));
}


@ We allow implicitly converting a parameter to a virtual module. This invokes
conversion by the |Rep_context::expand_final| method to \emph{final}
parameters (there can be zero, one, or more of them), to initiate the
invariant that only standard nonzero final parameters with dominant $\gamma$
can be stored in a |virtual_module_value| (the |expand_final| method calls
|normalise| internally, so we don't have to do that here).

@< Local function def...@>=
void param_to_poly()
{ shared_module_parameter p = get<module_parameter_value>();
@/test_standard(*p,"Cannot convert non standard Param to ParamPol");
  const own_real_form& rf=p->rf;
  push_value(std::make_shared<virtual_module_value> @|
    (rf,rf->rc().expand_final(p->val)));
}

@ There also is function to extract the coefficient (multiplicity) of a given
parameter in a virtual module. However, it is bound to the array subscription
syntax, and therefore does not have a wrapper function. Instead, it is
implemented the \.{axis} module, as the |evaluate| method of the
|module_coefficient| class derived from |subscr_base|.

In a subscription of a polynomial by a parameter, the arguments are not
initially on the stack, but come from evaluating the |array| and |index|
fields of the |module_coefficient| expression.

@h "axis.h" // for |module_coefficient|

@< Function def... @>=
void module_coefficient::evaluate(level l) const
{ shared_virtual_module m = (array->eval(),get<virtual_module_value>());
  shared_module_parameter p = (index->eval(),get<module_parameter_value>());
  if (m->rf!=p->rf)
    throw runtime_error @|
      ("Real form mismatch when subscripting ParamPol value");
  test_standard(*p,"In subscription of ParamPol value");
     // it is OK to do this test before |make_dominant|
  StandardRepr sr = p->val;
  p->rc().make_dominant(sr);
  if (l!=expression_base::no_value)
    push_value(std::make_shared<split_int_value>(m->val[sr]));
}

@ The main operations for virtual modules are addition and subtraction of
parameters to or from them.

@< Local function def...@>=
void add_module_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  own_virtual_module accumulator = get_own<virtual_module_value>();
@/test_standard(*p,"Cannot convert non standard Param to term in ParamPol");
  if (accumulator->rf!=p->rf)
    throw runtime_error @|
      ("Real form mismatch when adding standard module to a module");
  if (l!=expression_base::no_value)
  @/{@; accumulator->val+= p->rc().expand_final(p->val);
    push_value(accumulator);
  }
}

void subtract_module_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  own_virtual_module accumulator = get_own<virtual_module_value>();
@/test_standard(*p,"Cannot convert to term in ParamPol");
  if (accumulator->rf!=p->rf)
    throw runtime_error @|
      ("Real form mismatch when subtracting standard module from a module");
  if (l!=expression_base::no_value)
  @/{@; accumulator->val-= p->rc().expand_final(p->val);
    push_value(accumulator);
  }
}

@ More generally than adding or subtracting, we can incorporate a term with
specified coefficient. Here, rather than building a polynomial with the set of
parameters but without the coefficient,as |expand_final| gives us, we directly
iterate over the list produced by |finals_for|, attaching |coef| for each of
its final parameters.

@< Local function def...@>=

void add_module_term_wrapper(expression_base::level l)
{ push_tuple_components(); // second argument is a pair |(coef,p)|
  own_module_parameter p = get_own<module_parameter_value>();
  Split_integer coef=get<split_int_value>()->val;
  own_virtual_module accumulator = get_own<virtual_module_value>();
@/test_standard(*p,"Cannot convert non standard Param to term in ParamPol");
  if (accumulator->rf!=p->rf)
    throw runtime_error @|
      ("Real form mismatch when adding a term to a module");
  if (l==expression_base::no_value)
    return;
@)
  auto finals = p->rc().finals_for(p->val);
  for (auto it=finals.wcbegin(); not finals.at_end(it); ++it)
    accumulator->val.add_term(*it,coef);
  push_value(accumulator);
}

@ Although we initially envisioned allowing conversion from a list of terms to
a virtual module, this could not be defined since it is not possible to know
the real form in case the list of terms is empty (a conversion in the opposite
direction is given below). Therefore we provide instead the addition of an
entire list of terms at once to a virtual module value.

@< Local function... @>=
void add_module_termlist_wrapper(expression_base::level l)
{ shared_row r = get<row_value>();
  own_virtual_module accumulator = get_own<virtual_module_value>();
  if (l==expression_base::no_value)
    return;
@)
  for (auto it=r->val.cbegin(); it!=r->val.cend(); ++it)
  { const tuple_value* t = force<tuple_value>(it->get());
    Split_integer coef=force<split_int_value>(t->val[0].get())->val;
    const module_parameter_value* p =
      force<module_parameter_value>(t->val[1].get());
@/    test_standard(*p,"Cannot convert non standard Param to term in ParamPol");
    if (accumulator->rf!=p->rf)
      throw runtime_error @|
        ("Real form mismatch when adding terms to a module");
     auto finals = p->rc().finals_for(p->val);
     for (auto it=finals.wcbegin(); not finals.at_end(it); ++it)
       accumulator->val.add_term(*it,coef);
   }
  push_value(accumulator);
}

@ Naturally we also want to define addition and scalar multiplication of
virtual modules.
@< Local function... @>=
void add_virtual_modules_wrapper(expression_base::level l)
{
  own_virtual_module accumulator = get_own<virtual_module_value>();
  shared_virtual_module addend = get<virtual_module_value>();
  if (accumulator->rf!=addend->rf)
    throw runtime_error @|("Real form mismatch when adding two modules");
  if (l!=expression_base::no_value)
  @/{@; accumulator->val += addend->val;
    push_value(accumulator);
  }
}
@)
void subtract_virtual_modules_wrapper(expression_base::level l)
{
  shared_virtual_module subtrahend = get<virtual_module_value>();
  own_virtual_module accumulator = get_own<virtual_module_value>();
  if (accumulator->rf!=subtrahend->rf)
    throw runtime_error @|
      ("Real form mismatch when subtracting two modules");
  if (l!=expression_base::no_value)
  @/{@; accumulator->val -= subtrahend->val;
    push_value(accumulator);
  }
}

@ Scalar multiplication potentially makes coefficients zero, in which case the
corresponding terms need to be removed to preserve the invariant that no zero
terms are stored in a virtual module. For integer multiplication we just need
to check for multiplication by $0$, and produce an empty module when this
happens. Because the integer is not on the stack top, this requires a somewhat
unusual manoeuvre. The case of multiplication by zero needs to be handled
separately, since we cannot allow introducing terms with zero coefficients. It
could have been handled more easily though, by testing the factor~|c| just
before the |for| loop, and performing |m->erase()| instead if |c==0|; this is
what we used to do. However that might involve duplicating the virtual module
and then erasing the copy, which is inefficient, and now avoided. This might
seem a rare case, but it is not really: often functions handling
a \.{ParamPol} argument $P$ need to start with an empty module for the same real
form; writing $0*P$ is quite a convenient way to achieve this.

Matters are similar but somewhat subtler for scalar multiplication by split
integers, because these have zero divisors. Therefore we need to
test \emph{each} coefficient produced by multiplication in this case, and
remove the term when the coefficient becomes zero. We must take care to
advance the iterator ``manually'' before doing that, and as a consequence
cannot as usual advance the iterator in the |for| clause. In this case we do
not bother handling the case of an entirely zero split integer
multiplier separately; that case \emph{is} rare, and the given code works
correctly for it (albeit not in the fastest possible way).

@< Local function... @>=

void int_mult_virtual_module_wrapper(expression_base::level l)
{ int c =
    force<int_value>(execution_stack[execution_stack.size()-2].get())->int_val();
  // below top
  if (c==0) // then do multiply by $0$ efficiently:
  { shared_virtual_module m = get<virtual_module_value>();
      // |m| is needed for |m->rc()|
    pop_value();
    if (l!=expression_base::no_value)
    @/push_value@|(std::make_shared<virtual_module_value>
        (m->rf,repr::SR_poly(m->rc().repr_less())));
  }
  else
  { own_virtual_module m = get_own<virtual_module_value>();
     // will modify our copy now
    pop_value();
    assert(c!=0); // we tested that above
    if (l!=expression_base::no_value)
    { for (repr::SR_poly::iterator it=m->val.begin(); it!=m->val.end(); ++it)
        it->second *= c;
      push_value(m);
    }
  }
}
@)
void split_mult_virtual_module_wrapper(expression_base::level l)
{ own_virtual_module m = get_own<virtual_module_value>();
  Split_integer c = get<split_int_value>()->val;
  if (l==expression_base::no_value)
    return;
@)
  for (repr::SR_poly::iterator it=m->val.begin(); it!=m->val.end(); )
    // no |++it| here!
    if ((it->second *= c)==Split_integer(0,0))
      m->val.erase(it++); // advance, then delete the node just abandoned
    else ++it;
  push_value(m);
}

@ For a nonzero virtual module, it is useful to be able to select a component
that is present without looping over all terms. The most useful choice it the
final term, be we allow taking the first term as well.

@< Local function... @>=
void last_term_wrapper (expression_base::level l)
{ shared_virtual_module m = get<virtual_module_value>();
  if (l==expression_base::no_value)
    return;
@)
  if (m->val.empty())
    throw runtime_error("Empty module has no last term");
  const auto& term = *m->val.rbegin();
  push_value(std::make_shared<split_int_value>(term.second));
  push_value(std::make_shared<module_parameter_value>(m->rf,term.first));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}
@)
void first_term_wrapper (expression_base::level l)
{ shared_virtual_module m = get<virtual_module_value>();
  if (l==expression_base::no_value)
    return;
@)
  if (m->val.empty())
    throw runtime_error("Empty module has no first term");
  const auto& term = *m->val.begin();
  push_value(std::make_shared<split_int_value>(term.second));
  push_value(std::make_shared<module_parameter_value>(m->rf,term.first));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ Here are variations of the scaling functions for parameters that operate on
entire virtual modules.

@< Local function def...@>=
void scale_poly_wrapper(expression_base::level l)
{ shared_rat f = get<rat_value>();
  shared_virtual_module P = get<virtual_module_value>();
  if (l!=expression_base::no_value)
    push_value@|(std::make_shared<virtual_module_value>
      (P->rf,P->rc().scale(P->val,f->rat_val())));
}

void scale_0_poly_wrapper(expression_base::level l)
{ shared_virtual_module P = get<virtual_module_value>();
  if (l!=expression_base::no_value)
    push_value@|(std::make_shared<virtual_module_value>
      (P->rf,P->rc().scale_0(P->val)));
}

@*2 Computing with $K$-types.
%
The ``restriction to $K$'' component of the Atlas library has for a long time
had an isolated existence, in part because the ``nonzero final standard
$K$ parameters'' used to designate $K$-types are both hard to decipher and do
not seem to relate easily to values used elsewhere. However, these parameters
do in fact correspond to the subset of module pareters with $\nu=0$.

The function |K_type_formula| converts its parameter to a |StandardRepK|
value~|srk|, for which it calls the method |SRK_context::K_type_formula|. Since
the terms of that formula are not necessarily standard, one needs to pass each
term though |KhatContext::standardize|; finally, to incorporate the resulting
terms into a |virtual_module|, they need to be passed through
|Rep_context::expand_final|. One must not forget to multiply the small integer
coefficients that these two methods produce.

@h "standardrepk.h"

@< Local function def...@>=

void K_type_formula_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  RealReductiveGroup& G = p->rf->val;
  const Rep_context& rc = p->rc();
  KhatContext& khc = p->rf->khc();
  StandardRepK srk =
    khc.std_rep_rho_plus (rc.lambda_rho(p->val),G.kgb().titsElt(p->val.x()));
  @< Check that |srk| is final, and if not |throw| an error @>
  if (l==expression_base::no_value)
    return;
@)
  const standardrepk::Char formula = khc.K_type_formula(srk).second;
   // don't need |first==srk|
  const RatWeight zero_nu(p->rf->val.rank());
@/own_virtual_module acc @|
    (new virtual_module_value(p->rf, repr::SR_poly(p->rc().repr_less())));
  for (auto it=formula.begin(); it!=formula.end(); ++it)
  {
    standardrepk::combination st=khc.standardize(it->first);
    for (auto stit=st.cbegin(); stit!=st.cend(); ++stit)
    {
      StandardRepr term =  rc.sr(khc.rep_no(stit->first),khc,zero_nu);
      Split_integer coef (it->second*stit->second);
      auto finals = p->rc().finals_for(term);
      for (auto jt=finals.wcbegin(); not finals.at_end(jt); ++jt)
         acc->val.add_term(*jt,coef);
    }
  }
  push_value(acc);
}

@ We must test the parameter for being final, or else the method
|K_type_formula| will fail. The error message mentions restriction to $K$,
since the parameter itself reported here might be final.

@< Check that |srk| is final, and if not |throw| an error @>=
{ if (not khc.isFinal(srk))
  { std::ostringstream os;
    print_stdrep(os << "Non final restriction to K: ",p->val,rc) @|
      << "\n  (witness " << khc.rootDatum().coroot(khc.witness()) << ')';
    throw runtime_error(os.str());
  }
}

@ A main function is the actual branching to~$K$: decomposition of a standard
representation into $K$-types, up to a given limit. This used to be limited to
final representations, but it turns out that |KhatContext::standardize| will
expand non-final |StandardRepK| values into sums of final ones, so we decided
to drop that restriction without otherwise changing the code below. Since
|KhatContext::standardize| also does what its name suggests, it seems likely
that we could drop the ``standard'' condition as well, though this will
involve a longer rewriting process (using Hecht-Schmid identities) than the
expansion into finals.

@< Local function def...@>=
void branch_wrapper(expression_base::level l)
{ int bound = get<int_value>()->int_val();
  // not ``branch and bound'' but ``branch up to bound''
  shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Branching of non-standard parameter is not allowed");
  const Rep_context rc = p->rc();
  RealReductiveGroup& G=p->rf->val;
  KhatContext& khc = p->rf->khc();
  StandardRepK srk=
    khc.std_rep_rho_plus (rc.lambda_rho(p->val),G.kgb().titsElt(p->val.x()));
  assert(khc.isStandard(srk)); // should be ensured by |test_standard|
  if (l==expression_base::no_value)
    return;
@)
  khc.normalize(srk);
  standardrepk::combination combo=khc.standardize(srk);
  RatWeight zero_nu(G.rank());
@/own_virtual_module acc @|
    (new virtual_module_value(p->rf, repr::SR_poly(rc.repr_less())));
  for (auto it=combo.begin(); it!=combo.end(); ++it)
    // loop over finals from |srk|
  {
    standardrepk::combination chunk = khc.branch(it->first,bound);
    for (auto jt=chunk.begin(); jt!=chunk.end(); ++jt)
    {
      StandardRepr z = rc.sr(khc.rep_no(jt->first),khc,zero_nu);
      acc->val.add_term(z,Split_integer(it->second*jt->second));
    }
  }
  push_value(acc);
}

@ One can also branch from a polynomial. Since terms of a polynomial are
guaranteed to be final, we can transform its terms to |StandardRepK| values
that can be directly fed to |KhatContext::branch| without passing through
|standardize|. We did need to add a method |KhatContext::match_final| to
transform |srk| into a sequence number inside |khc|, which is what branching
needs.

@< Local function def...@>=
void branch_pol_wrapper(expression_base::level l)
{ int bound = get<int_value>()->int_val();
  // not ``branch and bound'' but ``branch up to bound''
  shared_virtual_module P = get<virtual_module_value>();
  if (l==expression_base::no_value)
    return;
@)
  RealReductiveGroup& G=P->rf->val;
  const Rep_context rc = P->rc();
  KhatContext& khc = P->rf->khc();
  auto P0 = rc.scale_0(P->val);
@/own_virtual_module acc @|
    (new virtual_module_value(P->rf, repr::SR_poly(rc.repr_less())));
  RatWeight zero_nu(G.rank());
  for (auto it=P0.begin(); it!=P0.end(); ++it)
    // loop over terms of |P0|
  {
    StandardRepK srk= khc.std_rep_rho_plus
       (rc.lambda_rho(it->first),G.kgb().titsElt(it->first.x()));
    assert(khc.isNormal(srk));
    standardrepk::combination chunk = khc.branch(khc.match_final(srk),bound);
    for (auto jt=chunk.begin(); jt!=chunk.end(); ++jt)
    {
      StandardRepr z = rc.sr(khc.rep_no(jt->first),khc,zero_nu);
      acc->val.add_term(z,it->second*jt->second);
    }
  }
  push_value(acc);
}

@ In the K-type code, standard representations restricted to $K$ are always
given on the canonical twisted involution of their Cartan class. In order to
be able to understand what a parameter will look like in this representation,
we provide a function that performs a similar transformation of a parameter,
ignoring its $\nu$ component. The operation simply consists of applying
complex cross actions on~$x$ until it is canonical for its Cartan class, and
applying the corresponding simple reflections to~$\lambda$.

@< Local function def...@>=
void to_canonical_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  const InnerClass& G=p->rf->val.innerClass();
  const KGB& kgb = p->rf->kgb();
  const RootDatum& rd=G.rootDatum();
@)
  KGBElt x = p->val.x();
  TwistedInvolution sigma = kgb.involution(x);
  WeylWord w = G.canonicalize(sigma); // $x\times w$ lies over |sigma| canonical
  Weight two_lambda = p->rc().lambda_rho(p->val)*2 + rd.twoRho();
@)
  x = kgb.cross(x,w);
  assert(kgb.involution(x)==sigma);
   // we should now be at canonical twisted involution
  rd.act_inverse(two_lambda,w);
  if (l==expression_base::no_value)
    return;
@)
  RatWeight zero_nu(p->rf->val.rank());
  StandardRepr result = p->rc().sr(x,(two_lambda-rd.twoRho())/2,zero_nu);
  push_value(std::make_shared<module_parameter_value>(p->rf,result));
}

@ Here is one more useful function: computing the height of a parameter
(ignoring the $\nu$ component). This is the same height displayed when
printing \.{ParamPol} values, and it is also used when comparing to the
|bound| argument to |branch| above. While this used to be obtained from the
method |SRK_context::height|, the height is now stored inside |StandardRepr|
values themselves, so we get it from there.

@< Local function def...@>=
void srk_height_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(p->val.height()));
}

@*2 Deformation formulas.
Here is one important application of virtual modules.
%
Using non-integral blocks, we can compute a deformation formula for the given
parameter. This also involves computing Kazhdan-Lusztig polynomials, which
happens inside the method |Rep_table::deformation_terms|, and produces an
|SR_poly| describing the module that is ``split off'' from the standard module
when the continuous part $\mu$ of the parameter is infinitesimally
``deformed'' towards~$0$; these terms involve certain other parameters found
below the parameter in its block. The code below used to apply |expand_final|
to the |deformation_terms|, but that is redundant since that method already
condenses the KL polynomials (and its result) to block elements without
singular descents (so nonzero and final), for which |expand_final| has no
effect.

There is also a variation |twisted_deform| that uses twisted KLV polynomials
instead, for the distinguished involution $\delta$ of the inner class. For the
code here the difference consists mainly of calling the
|Rep_table::twisted_deformation_terms| method instead of
|Rep_table::deformation_terms|. However, that method requires a $\delta$-fixed
involution, so we need to test for that here. If the test fails we report an
error rather than returning for instance a null module, since a twisted
deformation formula for a non-fixed parameter makes little sense; the user
should avoid asking for it. Also, since the construction of an extended block
currently cannot deal with a partial parent block.

@< Local function def...@>=
void deform_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot compute deformation");
  if (l==expression_base::no_value)
    return;
@)
  param_block block(p->rc(),p->val); // partial block construction
  repr::SR_poly terms
     = p->rt().deformation_terms(block,block.size()-1);

  push_value(std::make_shared<virtual_module_value>(p->rf,std::move(terms)));
}
@)
void twisted_deform_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  const auto& rc=p->rc();
  test_standard(*p,"Cannot compute twisted deformation");
  if (not rc.is_twist_fixed(p->val,rc.innerClass().distinguished()))
    throw runtime_error("Parameter not fixed by inner class involution");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt entry_elem;
  param_block block(p->rc(),p->val,entry_elem); // full block
  repr::SR_poly terms
     = p->rt().twisted_deformation_terms(block,entry_elem);

  push_value(std::make_shared<virtual_module_value>(p->rf,std::move(terms)));
}

@ Here is a recursive form of this deformation, which stores intermediate
results for efficiency in the |Rep_table| structure |p->rt()| that is stored
within the |real_form_value|.

@< Local function def...@>=
void full_deform_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot compute full deformation");
  if (l==expression_base::no_value)
    return;
@)
  const auto& rc = p->rc();
  auto finals = rc.finals_for(p->val);
  repr::SR_poly result (rc.repr_less());
  for (auto it=finals.cbegin(); it!=finals.cend(); ++it)
    result += p->rt().deformation(*it);
  push_value(std::make_shared<virtual_module_value>(p->rf,result));
}
@)
void twisted_full_deform_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  const auto& rc=p->rc();
  test_standard(*p,"Cannot compute full twisted deformation");
  auto sr=p->val; // take a copy
  rc.make_dominant(sr); // |is_twist_fixed| and |extended_finalise| like this
  if (not rc.is_twist_fixed(sr))
    throw runtime_error("Parameter not fixed by inner class involution");
  if (l==expression_base::no_value)
    return;
@)
  auto finals =
    ext_block::extended_finalise(rc,sr,rc.innerClass().distinguished());
  repr::SR_poly result (rc.repr_less());
  for (auto it=finals.cbegin(); it!=finals.cend(); ++it)
    result.add_multiple(p->rt().twisted_deformation(it->first) @|
                       ,it->second ? Split_integer(0,1) : Split_integer(1,0));
  push_value(std::make_shared<virtual_module_value>(p->rf,result));
}

@ And here is another way to invoke the Kazhdan-Lusztig computations, which
given a parameter corresponding to $y$ will obtain the formal sum over $x$ in
the block of $y$ (or the Bruhat interval below $y$, where all those giving a
nonzero contribution are located), of the parameter $x$ with as coefficient
the Kazhdan-Lusztig polynomials $P_{x,y}$ multiplied by a sign and evaluated
at the split integer unit~$s$ (since it appears that the information most
frequently needed can be extracted from that evaluation). In formula, this
computes
$$
  \sum_{x\leq y}(-1)^{l(y)-l(x)}P_{x,y}[q:=s] * x
$$
There are in fact two variants, of this function an ordinary one and one using
twisted KLV polynomials, computed for the inner class involution.

@< Local function def...@>=
void KL_sum_at_s_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot compute Kazhdan-Lusztig sum");
  test_normal_is_final(*p,"Cannot compute Kazhdan-Lusztig sum");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<virtual_module_value>@|
      (p->rf,p->rt().KL_column_at_s(p->val)));
}
@)
void twisted_KL_sum_at_s_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot compute Kazhdan-Lusztig sum");
  test_normal_is_final(*p,"Cannot compute Kazhdan-Lusztig sum");
  auto sr=p->val; // take a copy
  p->rc().make_dominant(sr);
    // |is_twist_fixed| and |twisted_KL_column_at_s| like this
  if (not p->rc().is_twist_fixed(sr))
    throw runtime_error("Parameter not fixed by inner class involution");
  if (l!=expression_base::no_value)
    push_value (std::make_shared<virtual_module_value>@|
      (p->rf,p->rt().twisted_KL_column_at_s(sr)));
}

@ We add another function in which the external involution is an argument

@< Local function def...@>=
void external_twisted_KL_sum_at_s_wrapper(expression_base::level l)
{ shared_matrix delta = get<matrix_value>();
  shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot compute Kazhdan-Lusztig sum");
  test_normal_is_final(*p,"Cannot compute Kazhdan-Lusztig sum");
  test_compatible(p->rc().innerClass(),delta);
  if (not p->rc().is_twist_fixed(p->val,delta->val))
    throw runtime_error("Parameter not fixed by given involution");
  if (l!=expression_base::no_value)
    push_value (std::make_shared<virtual_module_value>@|
      (p->rf,twisted_KL_column_at_s(p->rc(),p->val,delta->val)));
}

@ The function |scale_extended| is intended for use with in the deformation
algorithm when interpreting parameters as specifying a representation of he
extended group. One can arrange that deformation starts with a parameter for
which $\gamma$ is dominant, but when deforming this condition may be lost. In
the ordinary deformation algorithm this is taken care of by an implicit
|dominant| conversion of the parameter when it gets used to construct a block
or when it is contributed to a virtual module. However this may potentially
cause the default choice of extended representation associated to the
parameter to flip (at the time of writing this, it appears to be a very rare
event, if it occurs at all). the function below will scale the parameter,
perform the conversion to dominant using extended parameters, and return the
scaled parameter made dominant plus an indication of whether a flip occurred.
The function |scaled_extended_dominant| in the module \\{ext\_block} does the
actual work.

@< Local function def...@>=

void scale_extended_wrapper(expression_base::level l)
{ auto factor = get<rat_value>();
  auto delta = get<matrix_value>();
  auto p = get<module_parameter_value>();
  const StandardRepr sr = p->val;
  const auto& rc = p->rc();
  test_standard(*p,"Cannot scale extended parameter");
  if (not is_dominant_ratweight(rc.rootDatum(),sr.gamma()))
    throw runtime_error("Parameter to be scaled not dominant");
  test_compatible(p->rc().innerClass(),delta);
  if (not rc.is_twist_fixed(sr,delta->val))
    throw runtime_error("Parameter to be scaled not fixed by given involution");
  if (l==expression_base::no_value)
    return;
@)
  bool flipped;
  auto result = @;ext_block::scaled_extended_dominant
    (rc,sr,delta->val,factor->rat_val(),flipped);
  push_value(std::make_shared<module_parameter_value>(p->rf,result));
  push_value(whether(flipped));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ The function |finalize_extended| is useful in expanding a single module
parameter into a linear combination of such parameters. The terms on the list
are paired with a Boolean attribute recording a possible flip of extended
parameters accumulated when the term. This flip is recorded with as
coefficient the split integer unit~$s$, since it should be interpreted as a
signature flip (in the ordinary finalisation procedure flips never occur).
The function |finalise| in the module \\{ext\_block} does the actual work.

@< Local function def...@>=

void finalize_extended_wrapper(expression_base::level l)
{ auto delta = get<matrix_value>();
  auto p = get<module_parameter_value>();
  const auto& rc = p->rc();
  test_standard(*p,"Cannot finalize extended parameter");
  test_compatible(rc.innerClass(),delta);
  if (not p->rc().is_twist_fixed(p->val,delta->val))
    throw runtime_error("Parameter not fixed by given involution");
  if (not is_dominant_ratweight(rc.rootDatum(),p->val.gamma()))
    throw runtime_error("Parameter must have dominant gamma");
  if (l==expression_base::no_value)
    return;
@)
  auto params = @;ext_block::extended_finalise(rc,p->val,delta->val);
  repr::SR_poly result(rc.repr_less());
  for (auto it=params.begin(); it!=params.end(); ++it)
    result.add_term(it->first
                   ,it->second ? Split_integer(0,1) : Split_integer(1,0));
  push_value (std::make_shared<virtual_module_value>(p->rf,std::move(result)));
}


@ Finally we install everything related to polynomials formed from parameters.
@< Install wrapper functions @>=
install_function(split_unary_eq_wrapper,@|"=","(Split->bool)");
install_function(split_unary_neq_wrapper,@|"!=","(Split->bool)");
install_function(split_eq_wrapper,@|"=","(Split,Split->bool)");
install_function(split_neq_wrapper,@|"!=","(Split,Split->bool)");
install_function(split_plus_wrapper,@|"+","(Split,Split->Split)");
install_function(split_minus_wrapper,@|"-","(Split,Split->Split)");
install_function(split_unary_minus_wrapper,@|"-","(Split->Split)");
install_function(split_times_wrapper,@|"*","(Split,Split->Split)");
install_function(from_split_wrapper,@|"%","(Split->int,int)");
install_function(virtual_module_wrapper,@|"null_module","(RealForm->ParamPol)");
install_function(real_form_of_virtual_module_wrapper,@|"real_form"
		,"(ParamPol->RealForm)");
install_function(virtual_module_size_wrapper,@|"#","(ParamPol->int)");
install_function(virtual_module_unary_eq_wrapper,@|"=","(ParamPol->bool)");
install_function(virtual_module_unary_neq_wrapper,@|"!=","(ParamPol->bool)");
install_function(virtual_module_eq_wrapper,@|"=","(ParamPol,ParamPol->bool)");
install_function(virtual_module_neq_wrapper,@|"!=","(ParamPol,ParamPol->bool)");
install_function(add_module_wrapper,@|"+","(ParamPol,Param->ParamPol)");
install_function(subtract_module_wrapper,@|"-","(ParamPol,Param->ParamPol)");
install_function(add_module_term_wrapper,@|"+"
		,"(ParamPol,(Split,Param)->ParamPol)");
install_function(add_module_termlist_wrapper,@|"+"
		,"(ParamPol,[(Split,Param)]->ParamPol)");
install_function(add_virtual_modules_wrapper,@|"+"
		,"(ParamPol,ParamPol->ParamPol)");
install_function(subtract_virtual_modules_wrapper,@|"-"
		,"(ParamPol,ParamPol->ParamPol)");
install_function(int_mult_virtual_module_wrapper,@|"*"
		,"(int,ParamPol->ParamPol)");
install_function(split_mult_virtual_module_wrapper,@|"*"
		,"(Split,ParamPol->ParamPol)");
install_function(last_term_wrapper,"last_term","(ParamPol->Split,Param)");
install_function(first_term_wrapper,"first_term","(ParamPol->Split,Param)");
install_function(scale_poly_wrapper,"*", "(ParamPol,rat->ParamPol)");
install_function(scale_0_poly_wrapper,"at_nu_0", "(ParamPol->ParamPol)");
install_function(K_type_formula_wrapper,@|"K_type_formula" ,"(Param->ParamPol)");
install_function(branch_wrapper,@|"branch" ,"(Param,int->ParamPol)");
install_function(branch_pol_wrapper,@|"branch" ,"(ParamPol,int->ParamPol)");
install_function(to_canonical_wrapper,@|"to_canonical" ,"(Param->Param)");
install_function(srk_height_wrapper,@|"height" ,"(Param->int)");
install_function(deform_wrapper,@|"deform" ,"(Param->ParamPol)");
install_function(twisted_deform_wrapper,@|"twisted_deform" ,"(Param->ParamPol)");
install_function(full_deform_wrapper,@|"full_deform","(Param->ParamPol)");
install_function(twisted_full_deform_wrapper,@|"twisted_full_deform"
                ,"(Param->ParamPol)");
install_function(KL_sum_at_s_wrapper,@|"KL_sum_at_s","(Param->ParamPol)");
install_function(twisted_KL_sum_at_s_wrapper,@|"twisted_KL_sum_at_s"
                ,"(Param->ParamPol)");
install_function(external_twisted_KL_sum_at_s_wrapper,@|"twisted_KL_sum_at_s"
                ,"(Param,mat->ParamPol)");
install_function(scale_extended_wrapper,@|"scale_extended"
                ,"(Param,mat,rat->Param,bool)");
install_function(finalize_extended_wrapper,@|"finalize_extended"
                ,"(Param,mat->ParamPol)");


@*1 Kazhdan-Lusztig tables. We implement a simple function that gives raw
access to the table of Kazhdan-Lusztig polynomials.

@< Local function def...@>=
void raw_KL_wrapper (expression_base::level l)
{ own_Block b = non_const_get<Block_value>();
  const Block& block = b->val;
  if (l==expression_base::no_value)
    return;
@)
  b->klc.fill(false); // this does the actual KL computation
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(b->klc.size()));
  for (size_t y=1; y<b->klc.size(); ++y)
    for (size_t x=0; x<y; ++x)
      M->val(x,y) = b->klc.KL_pol_index(x,y);
@)
  own_row polys = std::make_shared<row_value>(0);
  polys->val.reserve(b->klc.polStore().size());
  for (size_t i=0; i<b->klc.polStore().size(); ++i)
  {
    const kl::KLPol& pol = b->klc.polStore()[i];
    std::vector<int> coeffs(pol.size());
    for (size_t j=pol.size(); j-->0; )
      coeffs[j]=pol[j];
    polys->val.emplace_back(std::make_shared<vector_value>(coeffs));
  }
@)
  std::vector<int> length_stops(block.length(block.size()-1)+2);
  length_stops[0]=0;
  for (size_t i=1; i<length_stops.size(); ++i)
    length_stops[i]=block.length_first(i);
@)
  push_value(std::move(M));
  push_value(std::move(polys));
  push_value(std::make_shared<vector_value>(length_stops));
  if (l==expression_base::single_value)
    wrap_tuple<3>();
}

@ For testing, it is useful to also have the dual Kazhdan-Lusztig tables. In
this case we cannot of course use the field |b->klc| to store the KL
polynomials, so we here us a local |kl::KLContext| variable.

@< Local function def...@>=
void raw_dual_KL_wrapper (expression_base::level l)
{ shared_Block b = get<Block_value>();
  const Block& block = b->val;
  Block dual_block = Block::build(b->dual_rf->val,b->rf->val);

  std::vector<BlockElt> dual=blocks::dual_map(block,dual_block);
  kl::KLContext klc(dual_block); klc.fill(false);
  if (l==expression_base::no_value)
    return;
@)
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(klc.size()));
  for (size_t y=1; y<klc.size(); ++y)
    for (size_t x=0; x<y; ++x)
      M->val(x,y) = klc.KL_pol_index(dual[y],dual[x]);
@)
  own_row polys = std::make_shared<row_value>(0);
  polys->val.reserve(klc.polStore().size());
  for (size_t i=0; i<klc.polStore().size(); ++i)
  {
    const kl::KLPol& pol = klc.polStore()[i];
    std::vector<int> coeffs(pol.size());
    for (size_t j=pol.size(); j-->0; )
      coeffs[j]=pol[j];
    polys->val.emplace_back(std::make_shared<vector_value>(coeffs));
  }
@)
  std::vector<int> length_stops(block.length(block.size()-1)+2);
  length_stops[0]=0;
  for (size_t i=1; i<length_stops.size(); ++i)
    length_stops[i]=block.length_first(i);
@)
  push_value(std::move(M));
  push_value(std::move(polys));
  push_value(std::make_shared<vector_value>(length_stops));
  if (l==expression_base::single_value)
    wrap_tuple<3>();
}

@ In order to have access to extended KLV polynomials, one has the following
wrapper function, which is similar to |raw_KL_wrapper|, but it takes a
parameter and an involution matrix as argument.

@h "ext_kl.h"

@< Local function def...@>=
void raw_ext_KL_wrapper (expression_base::level l)
{ auto delta = get<matrix_value>();
  shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  test_compatible(p->rc().innerClass(),delta);
  if (l==expression_base::no_value)
    return;
@)
  const auto& rc = p->rc();
  BlockElt start;
  param_block block(rc,p->val,start);
  if (not((delta->val-1)*block.gamma().numerator()).isZero())
  { // block not globally stable, so return empty values;
    push_value(std::make_shared<matrix_value>(int_Matrix()));
    push_value(std::make_shared<row_value>(0));
    push_value(std::make_shared<vector_value>(int_Vector()));
  }
  else
  {
    ext_block::ext_block eb(rc.innerClass(),block,delta->val);
    std::vector<Polynomial<int> > pool;
    ext_kl::KL_table klt(eb,pool); klt.fill_columns();
  @)
    own_matrix M = std::make_shared<matrix_value>(int_Matrix(klt.size()));
    for (size_t y=1; y<klt.size(); ++y)
      for (size_t x=0; x<y; ++x)
      @/{@; auto inx = klt.KL_pol_index(x,y);
        M->val(x,y) = inx.second ? -inx.first : inx.first;
      }
  @)
    own_row polys = std::make_shared<row_value>(0);
    polys->val.reserve(pool.size());
    for (size_t i=0; i<pool.size(); ++i)
    {
      std::vector<int> coeffs(pool[i].begin(),pool[i].end());
      polys->val.emplace_back(std::make_shared<vector_value>(coeffs));
    }
  @)
    std::vector<int> length_stops(block.length(block.size()-1)+2);
    length_stops[0]=0;
    for (size_t i=1; i<length_stops.size(); ++i)
      length_stops[i]=eb.length_first(i);
  @)
    push_value(std::move(M));
    push_value(std::move(polys));
    push_value(std::make_shared<vector_value>(length_stops));
  }
  if (l==expression_base::single_value)
    wrap_tuple<3>();
}

@ For old style blocks, the user may want to get information about the
$W$-graph without manually extracting it from the KL table. The following
function provides this information in the form of a value of
type \.{[[int],[int,int]]}, a list of pairs consisting of a $\tau$-invariant
(as list of simple root indices) and a list of outgoing edges, each one a pair
of a destination vertex and a $\mu$-value labelling the edge.

@< Local function def...@>=
void W_graph_wrapper(expression_base::level l)
{ own_Block b = non_const_get<Block_value>();
  if (l=expression_base::no_value)
    return;
@)
  b->klc.fill(false); // this does the actual KL computation
  wgraph::WGraph wg = kl::wGraph(b->klc);
@)
  own_row vertices=std::make_shared<row_value>(0);
  @< Push to |vertices| a list of pairs for each element of |wg|, each
     consisting of a descent set and a list of outgoing labelled edges @>
  push_value(vertices);
}

@ The following code was isolated so that it can be reused below.

@< Push to |vertices| a list of pairs for each element of |wg|, each
   consisting of a descent set and a list of outgoing labelled edges @>=
vertices->val.reserve(wg.size());
for (size_t i = 0; i < wg.size(); ++i)
{ auto ds = wg.descent_set(i);
  own_row descents=std::make_shared<row_value>(0);
  descents->val.reserve(ds.count());
  for (auto it=ds.begin(); it(); ++it)
    descents->val.push_back(std::make_shared<int_value>(*it));
  own_row out_edges = std::make_shared<row_value>(0);
  out_edges->val.reserve(wg.degree(i));
  for (unsigned j=0; j<wg.degree(i); ++j)
  { auto tup = std::make_shared<tuple_value>(2);
    tup->val[0] = std::make_shared<int_value>(wg.edge_target(i,j));
    tup->val[1] = std::make_shared<int_value>(wg.coefficient(i,j));
  @/out_edges->val.push_back(tup);
  }
  auto tup = std::make_shared<tuple_value>(2);
  tup->val[0] = descents;
  tup->val[1] = out_edges;
  vertices->val.push_back(std::move(tup));
}

@ Outputting |W_cells| as row of vectors; (this function was originally
contributed by Jeff Adams).

@< Local function def...@>=
void W_cells_wrapper(expression_base::level l)
{ own_Block b = non_const_get<Block_value>();
  if (l=expression_base::no_value)
    return;
@)
  b->klc.fill(false); // this does the actual KL computation
  wgraph::WGraph wg = kl::wGraph(b->klc);
  wgraph::DecomposedWGraph dg(wg);
@)

  own_row cells=std::make_shared<row_value>(0);
  cells->val.reserve(dg.cellCount());
  for (size_t c = 0; c < dg.cellCount(); ++c)
  { auto& wg=dg.cell(c); // local W-graph of cell
    own_row members =std::make_shared<row_value>(0);
    { const BlockEltList& mem=dg.cellMembers(c);
        // list of members of strong component |c|
      members->val.reserve(mem.size());
      for (auto it=mem.begin(); it!=mem.end(); ++it)
        members->val.push_back(std::make_shared<int_value>(*it));
    }
    own_row vertices=std::make_shared<row_value>(0);
    @< Push to |vertices| a list of pairs for each element of |wg|, each
       consisting of a descent set and a list of outgoing labelled edges @>
    auto tup = std::make_shared<tuple_value>(2);
    tup->val[0] = members;
    tup->val[1] = vertices;
    cells->val.push_back(std::move(tup));
  }
  push_value(cells);
}


@* Test functions.
%
Now we shall make available some commands without actually creating new data
types. This means the values created to perform the computation will be
discarded after the output is produced; our functions will typically return a
$0$-tuple. This is probably not a desirable state of affairs, but in \.{Fokko}
it is no different (the corresponding commands are defined in \.{realmode.cpp}
and in \.{test.cpp} which hardly has any persistent data) so for the moment it
should not be so bad.

@ The \.{realweyl} and \.{strongreal} commands require a real form and a
compatible Cartan class.

@h "output.h"
@< Local function def...@>=
void print_realweyl_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  own_real_form rf= non_const_get<real_form_value>();
@)
  if (&rf->parent.val!=&cc->parent.val)
    throw runtime_error @|
    ("Inner class mismatch between arguments");
@.Inner class mismatch...@>
  BitMap b(rf->parent.val.Cartan_set(rf->val.realForm()));
  if (not b.isMember(cc->number))
    throw runtime_error @|
    ("Cartan class not defined for real form");
@.Cartan class not defined...@>
@)
  output::printRealWeyl (*output_stream,rf->val,cc->number);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@)
void print_strongreal_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
@)
 output::printStrongReal
    (*output_stream,cc->parent.val,cc->parent.interface,cc->number);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}


@
We shall next implement the \.{block} command.

@h "blocks.h"
@h "block_io.h"

@< Local function def...@>=
void print_block_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  const Block& block = b->val;
@)
  block.print_to(*output_stream,false);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@ We provide functions corresponding to the \.{blockd} and \.{blocku}
variations of \.{block}.

@< Local function def...@>=
void print_blockd_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  const Block& block = b->val;
@)
  block.print_to(*output_stream,true);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@)
void print_blocku_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  const Block& block = b->val;
@)
  block_io::printBlockU(*output_stream,block);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@ The \.{blockstabilizer} command has a slightly different calling scheme than
\.{block} and its friends, in that it requires a block and a Cartan class. The
block itself is not actually used, just the real form and dual real form it
holds. The signature of |output::printBlockStabilizer| is a bit strange,
as it requires a |RealReductiveGroup| argument for the real form, but only
numbers for the Cartan class and the dual real form (but this is
understandable, as information about the inner class must be transmitted in
some way). In fact it used to be even a bit stranger, in that the real form
was passed in the form of a |output::Interface| value, a class (no
longer existent, and not to be confused with |output::FormNumberMap|,
which does not specify a particular real form) that we do not use in this
program; since only the |realGroup| field of the |output::Interface| was
used in |output::printBlockStabilizer|, we have changed its parameter
specification to allow it to be called easily here.

@< Local function def...@>=
void print_blockstabilizer_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  shared_Block b = get<Block_value>();
@)
  output::printBlockStabilizer
   (*output_stream, @| b->rf->val,cc->number,b->dual_rf->val.realForm());
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@ The functions |print_KGB|, |print_KGB_order| and |print_KGB_graph| take only
a real form as argument.

@h "kgb.h"
@h "kgb_io.h"

@< Local function def...@>=
void print_KGB_wrapper(expression_base::level l)
{ own_real_form rf= non_const_get<real_form_value>();
@)
  *output_stream
    << "kgbsize: " << rf->val.KGB_size() << std::endl;
  const KGB& kgb=rf->kgb();
  kgb_io::var_print_KGB(*output_stream,rf->val.innerClass(),kgb);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}
@)
void print_KGB_order_wrapper(expression_base::level l)
{ own_real_form rf= non_const_get<real_form_value>();
@)
  *output_stream
    << "kgbsize: " << rf->val.KGB_size() << std::endl;
  kgb_io::printBruhatOrder(*output_stream,rf->val.Bruhat_KGB());
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}
@)
void print_KGB_graph_wrapper(expression_base::level l)
{ own_real_form rf= non_const_get<real_form_value>();
@)
  *output_stream
    << "kgbsize: " << rf->val.KGB_size() << std::endl;
  kgb_io::makeDotFile(*output_stream,rf->kgb(),rf->val.Bruhat_KGB());
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@ The function |print_X| even takes only an inner class as argument.

@< Local function def...@>=
void print_X_wrapper(expression_base::level l)
{ shared_inner_class ic = get<inner_class_value>();
@)
  InnerClass& G=ic->val;
  kgb::global_KGB kgb(G); // build global Tits group, "all" square classes
  kgb_io::print_X(*output_stream,kgb);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@ The function |print_KL_basis| behaves much like |print_block| as far as
parametrisation is concerned.

@h "kl.h"
@h "klsupport.h"
@h "kl_io.h"
@< Local function def...@>=
void print_KL_basis_wrapper(expression_base::level l)
{ own_Block b = non_const_get<Block_value>();
  Block& block = b->val;
@)
  b->klc.fill(false); // this does the actual KL computation
  *output_stream
    << "Full list of non-zero Kazhdan-Lusztig-Vogan polynomials:\n\n";
  kl_io::printAllKL(*output_stream,b->klc,block);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@ The function |print_prim_KL| is a variation of |print_KL_basis|.

@< Local function def...@>=
void print_prim_KL_wrapper(expression_base::level l)
{ own_Block b = non_const_get<Block_value>();
  Block &block = b->val; // this one must be non-|const|
@)
  b->klc.fill(false); // this does the actual KL computation
  *output_stream
    << "Non-zero Kazhdan-Lusztig-Vogan polynomials for primitive pairs:\n\n";
  kl_io::printPrimitiveKL(*output_stream,b->klc,block);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@ The function |print_KL_list| is another variation of |print_KL_basis|, it
outputs just a list of all distinct Kazhdan-Lusztig-Vogan polynomials.

@< Local function def...@>=
void print_KL_list_wrapper(expression_base::level l)
{ own_Block b = non_const_get<Block_value>();
@)
  b->klc.fill(false); // this does the actual KL computation
  kl_io::printKLList(*output_stream,b->klc);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@ We close with two functions for printing the $W$-graph determined by the
polynomials computed. For |print_W_cells| we must construct one more object,
after having built the |klc::KLContext|.

@h "wgraph.h"
@h "wgraph_io.h"

@< Local function def...@>=
void print_W_cells_wrapper(expression_base::level l)
{ own_Block b = non_const_get<Block_value>();
@)
  b->klc.fill(false); // this does the actual KL computation
  wgraph::WGraph wg = kl::wGraph(b->klc);
  wgraph::DecomposedWGraph dg(wg);
@)
  wgraph_io::printWDecomposition(*output_stream,dg);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@ And here is |print_W_graph|, which just gives a variation on the output
routine of |print_W_cells|.

@< Local function def...@>=
void print_W_graph_wrapper(expression_base::level l)
{ own_Block b = non_const_get<Block_value>();
@)
  b->klc.fill(false); // this does the actual KL computation
  wgraph::WGraph wg = kl::wGraph(b->klc);
@)
  wgraph_io::printWGraph(*output_stream,wg);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}


@ Here we install all remaining wrapper functions.

@< Install wrapper functions @>=
install_function(raw_KL_wrapper,@|"raw_KL","(Block->mat,[vec],vec)");
install_function(raw_dual_KL_wrapper,@|"dual_KL","(Block->mat,[vec],vec)");
install_function(raw_ext_KL_wrapper,@|"raw_ext_KL","(Param,mat->mat,[vec],vec)");
install_function(W_graph_wrapper,@|"W_graph","(Block->[[int],[int,int]])");
install_function(W_cells_wrapper,@|"W_cells",
  "(Block->[[int],[[int],[int,int]]])");
@)
install_function(print_gradings_wrapper,@|"print_gradings"
		,"(CartanClass,RealForm->)");
install_function(print_realweyl_wrapper,@|"print_real_Weyl"
		,"(RealForm,CartanClass->)");
install_function(print_strongreal_wrapper,@|"print_strong_real"
		,"(CartanClass->)");
install_function(print_block_wrapper,@|"print_block","(Block->)");
install_function(print_blocku_wrapper,@|"print_blocku","(Block->)");
install_function(print_blockd_wrapper,@|"print_blockd","(Block->)");
install_function(print_blockstabilizer_wrapper,@|"print_blockstabilizer"
		,"(Block,CartanClass->)");
install_function(print_KGB_wrapper,@|"print_KGB","(RealForm->)");
install_function(print_KGB_order_wrapper,@|"print_KGB_order","(RealForm->)");
install_function(print_KGB_graph_wrapper,@|"print_KGB_graph","(RealForm->)");
install_function(print_X_wrapper,@|"print_X","(InnerClass->)");
install_function(print_KL_basis_wrapper,@|"print_KL_basis","(Block->)");
install_function(print_prim_KL_wrapper,@|"print_prim_KL","(Block->)");
install_function(print_KL_list_wrapper,@|"print_KL_list","(Block->)");
install_function(print_W_cells_wrapper,@|"print_W_cells","(Block->)");
install_function(print_W_graph_wrapper,@|"print_W_graph","(Block->)");

@* Installing coercions.
%
Finally we collect here all coercions related to specific Atlas types.

@< Install coercions @>=
{
  coercion(str_type,Lie_type_type,"LT",Lie_type_coercion);
  coercion(ic_type,rd_type,"RdIc",inner_class_to_root_datum_coercion);
  coercion(rf_type,ic_type,"IcRf",real_form_to_inner_class_coercion);
  coercion(rf_type,rd_type,"RdRf",real_form_to_root_datum_coercion);
  coercion(int_type,split_type,"SpI",int_to_split_coercion);
  coercion(int_int_type,split_type,"Sp(I,I)",pair_to_split_coercion);
  coercion(param_type,param_pol_type,"PolP",param_to_poly);
}



@* Index.

% Local IspellDict: british
