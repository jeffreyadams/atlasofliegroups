% Copyright (C) 2006-2021 Marc van Leeuwen
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
namespace {@; @< Local function definitions @>@; }
@< Function definitions @>@;
}}

@ As usual the external interface is written to the header file associated to
this file.

@( atlas-types.h @>=

#ifndef ATLAS_TYPES_H
#define ATLAS_TYPES_H

#include "../Atlas.h" // must be very first \.{atlas} include

@< Includes needed in the header file @>@;
namespace atlas { namespace interpreter {
@< Type definitions @>@;
@< Declarations of exported functions @>@;
}}
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
are obligatory for a primitive type, though the copy constructor is only needed
because we use |get_own<Lie_type_value>| below.

@< Type definitions @>=
struct Lie_type_value : public value_base
{ LieType val;
@)
  Lie_type_value() : val() @+ {}
    // default constructor, produces empty type
  Lie_type_value(const LieType& t) : val(t) @+{}
    // constructor from already validated Lie type
  Lie_type_value(LieType&& t) : val(std::move(t)) @+{} // idem, rvalue reference
@)
  virtual void print(std::ostream& out) const;
  static const char* name() @+{@; return "Lie type"; }
  Lie_type_value (const Lie_type_value& v) = default;
    // we use |get_own<Lie_type_value>|
@)
  void add_simple_factor (char,unsigned int); // grow
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
|std::pair<char,unsigned int>|. Therefore these types could in principle take
arbitrary values, not necessarily meaningful ones. To ensure that this cannot
happen to \.{atlas} users, we make the method |add_simple_factor|, which is
invoked to build up Lie types, checks for the validity.

Since the tests defined in \.{io/interactive\_lietype.cpp} used in \.{Fokko} are
clumsy to use, we prefer to perform our own tests here, which emulate
|interactive_lietype::checkSimpleLieType|. One can specify torus factors of
rank~$r>1$, but they are equivalent to $r$ torus factors of rank~$1$, and it
simplifies the software if we rewrite the former form to the latter on input, so
that is what we do here.

@h <string>
@h "constants.h"

@< Function definitions @>=
void Lie_type_value::add_simple_factor (char c,unsigned int rank)
{ static const std::string types=lietype::typeLetters; // |"ABCDEFGT"|
  auto t=types.find(c);
  if (t==std::string::npos)
  { std::ostringstream o;
    o << "Invalid type letter '" << c << '\'';
    throw runtime_error(o.str());
  }
@.Invalid type letter@>
  const unsigned int r=constants::RANK_MAX; // for convenience
@/static const unsigned int lwb[]={1,2,2,4,6,4,2,0};
  static const unsigned int upb[]={r,r,r,r,8,4,2,r};
  if (rank<lwb[t])
  { std::ostringstream o;
    o << "Too small rank " << rank << " for Lie type " << c;
    throw runtime_error(o.str());
  }
@.Too small rank@>
  if (rank>upb[t])
  { if (upb[t]!=r)
    { std::ostringstream o;
      o << "Too large rank " << rank << " for Lie type " << c;
      throw runtime_error(o.str());
    }
@.Too large rank@>
    else
    { std::ostringstream o;
      o << "Rank " << rank @|
        << " exceeds implementation limit " << r;
      throw runtime_error(o.str());
    }
@.Rank exceeds implementation limit@>
  }
@)
  if (val.rank() + rank > r)
  { std::ostringstream o;
    o << "Total rank " << val.rank() + rank @|
     << " exceeds implementation limit " << r;
    throw runtime_error(o.str());
  }
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
  { unsigned int rank;
    if (is.peek()=='-' or not (is>>rank)) // explicitly forbid minus sign
    { std::ostringstream o;
      o << "Error in string '" << is.str()
        << "' that should specify a Lie type";
@.Error in type string@>
      throw runtime_error(o.str());
    }
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
  { std::ostringstream o;
    o << "Combined rank " << t1->val.rank()+t2->val.rank() @|
        << " exceeds implementation limit " @|
        << static_cast<unsigned>(constants::RANK_MAX);
    throw runtime_error(o.str());
  }
@.Combined rank exceeds...@>
  if (l==expression_base::no_value)
    return;
@)
  t1->val.append(t2->val); // both factors have been tested, so this is safe
  push_value(std::move(t1));
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
    push_value(std::move(t));
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

A user might want to find out whether a call to |Cartan_matrix_type| (maybe
implicit in the construction of a root datum) will succeed, so that an error
stop can be avoided in case it will not. The function |is_Cartan_matrix| serves
this purpose, and is implemented by just calling |dynkin::lieType| as above, but
ignoring the result but catching any |error::Cartan_error| thrown and returning
|false| in that case.

@h "dynkin.h"
@< Local function definitions @>=
void type_of_Cartan_matrix_wrapper (expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  Permutation pi;
  LieType lt=dynkin::Lie_type(m->val,pi);
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<Lie_type_value>(lt));
  own_row perm = std::make_shared<row_value>(0);
  perm->val.reserve(pi.size());
  for(auto it=pi.begin(); it!=pi.end(); ++it)
    perm->val.push_back(std::make_shared<int_value>(*it));
  push_value(std::move(perm));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

void is_Cartan_matrix_wrapper (expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  Permutation pi;
  try
  {
    dynkin::Lie_type(m->val,pi);
    push_value(whether(true));
  }
  catch (error::Cartan_error&)
  {
    push_value(whether(false));
  }
}


@ For programming it is important to be able to analyse a Lie type. To this end
we allow transforming it into a list of $(code,rank)$ pairs, where |code| is a
one-letter string.

@< Local function definitions @>=
void simple_factors_wrapper(expression_base::level l)
{ shared_Lie_type t=get<Lie_type_value>();
  if (l==expression_base::no_value)
    return;
@)
  own_row result = std::make_shared<row_value>(0);
  result->val.reserve(t->val.size());
  for (unsigned i=0; i<t->val.size(); ++i)
  { const auto& src = t->val[i];
    auto dst=std::make_shared<tuple_value>(2);
    if (src.first!='T') // skip torus factors
    { dst->val[0] = std::make_shared<string_value>(std::string(1,src.first));
      dst->val[1] = std::make_shared<int_value>(src.second);
      result->val.push_back(std::move(dst));
    }
  }
  push_value(std::move(result));
}

void rank_of_Lie_type_wrapper
(expression_base::level l)
{ shared_Lie_type t=get<Lie_type_value>();
  if (l==expression_base::no_value)
    return;

  unsigned result=0;
  for (auto it = t->val.begin(); it!=t->val.end(); ++it)
    result += it->second;

  push_value(std::make_shared<int_value>(result));
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
install_function(is_Cartan_matrix_wrapper,@|"is_Cartan_matrix","(mat->bool)");
install_function(simple_factors_wrapper
                ,@|"simple_factors","(LieType->[string,int])");
install_function(rank_of_Lie_type_wrapper,"rank","(LieType->int)");

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
  unsigned int n=basis->val.n_columns();
  if (inv_f->val.size()>n)
  { std::ostringstream o;
    o << "Too many factors: " @|
      << inv_f->val.size() << " for " << n <<" columns";
    throw runtime_error(o.str());
  }
@.Too many factors@>
  if (l==expression_base::no_value)
    return;
@)
  unsigned int i=0;
  while (i<n)
    if (i>=inv_f->val.size() or inv_f->val[i]!=1)
      ++i; // keep invariant factor and column
    else
    {@; inv_f->val.erase(inv_f->val.begin()+i);
        basis->val.eraseColumn(i);
    }
  push_value(std::move(basis)); push_value(std::move(inv_f));
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

@h <cstdlib> // for |std::abs|

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

  for (unsigned int i=0; i<lambda.size(); ++i)
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
  unsigned int n=generators->val.n_columns();
  if (inv_f->val.size()>n)
  { std::ostringstream o;
    o << "Too many factors: " @|
      << inv_f->val.size() << " for " << n <<" columns";
    throw runtime_error(o.str());
  }
@.Too many factors@>
@)
  if (new_generators->val.n_rows()!=generators->val.n_rows())
    throw runtime_error("Column lengths do not match");
@.Column lengths do not match@>
@)
  unsigned int k=0; // index to replacement generators
  for (unsigned int j=0; j<n; ++j)
    if (j>=inv_f->val.size() or inv_f->val[j]!=1)
       // replace column |j| by |k| from |new_generators|
    { if (k>=new_generators->val.n_columns())
        throw runtime_error ("Not enough replacement columns");
@.Not enough replacement columns@>
      generators->val.set_column(j,new_generators->val.column(k));
      ++k;
    }
  if (k<new_generators->val.n_columns())
        throw runtime_error ("Too many replacement columns");
@.Too many replacement columns@>
  if (l!=expression_base::no_value)
    push_value(std::move(generators));
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
@/push_value(std::move(SC_basis));
  push_value(std::move(invf));
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
    { std::ostringstream o;
      o << "Length mismatch for generator " << j << ": "@|
        << gen.size() << ':' << v->val.size();
      throw runtime_error(o.str());
    }
@.Length mismatch...@>

    const auto& col=gen.numerator();
    for (unsigned int i=0; i<v->val.size(); ++i)
    { if (v->val[i]*col[i]%gen.denominator()!=0) // must use signed arithmetic!!
      { std::ostringstream o;
	o << "Improper generator entry: "
          << col[i] << '/' << denom[j] @|
          << " not a multiple of 1/" << v->val[i];
        throw runtime_error(o.str());
      }
@.Improper generator entry@>
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
  else
  { std::ostringstream o;
    o << "Unequal rank class is meaningless for type " @|
      << t << r;
    throw runtime_error(o.str());
  }
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
    { std::ostringstream o;
      o << "Permutation entry " << entry << " too big";
      throw runtime_error(o.str());
    }
    if (seen.isMember(entry))
    { std::ostringstream o;
      o << "Permutation has repeated entry " << entry;
      throw runtime_error(o.str());
    }
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
  { std::ostringstream o;
    o << "Permutation size " << perm->val.size() @|
      << " does not match rank " << t->val.rank() << " of Lie type";
    throw runtime_error(o.str());
  }
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
  if (basis->val.n_rows()!=r or basis->val.n_rows()!=r)
  { std::ostringstream o;
    o << "Basis should be given by " << r << 'x' << r << " matrix";
    throw runtime_error(o.str());
  }
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
  std::size_t hashCode(std::size_t modulus) const;
  bool operator!=(const root_datum_entry& x) const
    {@; return not PreRootDatum::operator==(x); }
};

@ The |root_datum_value| has as main purpose to wrap a |RootDatum| object into a
value derived from |value_base|. However, we want to make sure that identical
root data always end up using the \emph{same} |root_datum_value|; this is partly
to avoid wasting storage, but more importantly to be able to also easily
identify identical inner classes and other values based on them. In order to
achieve this, we include two kinds of special members in |root_datum_value|,
namely |static| ones to provide a look-up table for existing root data, and a
non-|static| member to provide a look-up table for inner classes currently built
upon this |root_datum_value|. In both cases we want to return a copy of an
existing shared pointer instead of creating a new one when this is possible, and
for this purpose we shall store weak pointers: storing a shared pointer would
effectively prevent the objects from ever being reclaimed even if nowhere else
any references remain, and from a raw pointer to the object pointed to by a
shared pointer no copy of the shared pointer can be obtained. Weak pointers are
exactly the right thing for the job here.

@< Type definitions @>=
class root_datum_value;
typedef std::shared_ptr<const root_datum_value> shared_root_datum;
typedef std::weak_ptr<const root_datum_value> root_datum_weak_ptr;
class inner_class_value;
typedef std::shared_ptr<const inner_class_value> shared_inner_class;
typedef std::weak_ptr<const inner_class_value> inner_class_weak_ptr;

@ Root data are immutable mathematical values, so the |val| member of
|root_datum_value| is public, but |const|. The copy constructor is deleted, and
we want all root data construction to take place through the |static| method
|build| that will try look-up first. In case nothing is found, it will need to
call a constructor, so we provide one that takes a |PreRootDatum| as ingredient.
The actual construction is done in the context of |std::make_shared|, not
directly by |build|, so the constructor needs to be public; however to ensure
that clients cannot circumvent |build|, we make the constructor require a
|token| that only methods of our class can supply.

The root datum look-up tables are in the |static| variables |pool|, |hash| and
|store|. Since few inner classes usually exist for a given root datum, we use a
simple linked list |classes| to record them, with the corresponding
distinguished involution. We also include a field for a |WeylGroup|, lazily
generated for clients that need it.

@s mutable const
@h "weyl.h"

@< Type definitions @>=
class root_datum_value : public value_base
{ struct token@+{}; // type passed to prove caller has private access;
  static root_datum_entry::Pooltype pool;
  static HashTable<root_datum_entry,unsigned int> hash;
  static std::vector<root_datum_weak_ptr> store;
  mutable simple_list<std::pair<const WeightInvolution,inner_class_weak_ptr> >
    classes;
  mutable std::shared_ptr<WeylGroup> W_ptr;
public:
  const RootDatum val;
@)
  root_datum_value(const PreRootDatum& v,token)
  : classes(), W_ptr(), val(v) @+ {}
  static shared_root_datum build(PreRootDatum&& pre);
  virtual void print(std::ostream& out) const;
  static const char* name() @+{@; return "root datum"; }
  root_datum_value (const root_datum_value& ) = delete;
@)
  shared_root_datum dual() const; // get dual datum through |build|
  inner_class_weak_ptr& lookup (const WeightInvolution& delta) const;
    // find or create empty
  const WeylGroup& W () const;
};

@ We need to define the static members declared in the class definition. The
|pool| serves as backing storage for |hash| as usual, and will not be addressed
directly; |store| starts out empty.

@< Global variable definitions @>=
root_datum_entry::Pooltype root_datum_value::pool;
HashTable<root_datum_entry,unsigned int> @| root_datum_value::hash
  (root_datum_value::pool);
std::vector<std::weak_ptr<const root_datum_value> > root_datum_value::store;

@ We have a simple hash function that uses all information in a |PreRootDatum|.
@< Function definitions @>=
std::size_t root_datum_entry::hashCode(std::size_t modulus) const
{ std::size_t h= prefer_coroots() ? 1 : 0;
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
@ The method |W| constructs a |WeylGroup| object if this is not already done,
and returns a reference to it.
%
The method |dual| facilitates building the dual of a |root_datum_value| while
passing though |build| to ensure that no unnecessary copies are created, notably
if one should later compute the dual of the dual. We want a root datum and its
dual to share their Weyl group, and the only way to ensure this is to generate
our Weyl group if necessary (by calling |W|) whenever the dual root datum is
formed.

@< Function definitions @>=

const WeylGroup& root_datum_value::W () const
{ if (W_ptr.get()==nullptr)
    W_ptr = std::make_shared<WeylGroup>(val.Cartan_matrix());
  return *W_ptr;
}
@)
shared_root_datum root_datum_value::dual() const
{ PreRootDatum pre(val); pre.dualise();
  auto result = build(std::move(pre));
  if (result->W_ptr.get()==nullptr)
    W(),result->W_ptr = W_ptr;
    // ensure we have a Weyl group, then share pointer
  return result;
}

@ When looking up to see whether an inner class is already known for a given
involution, we use linear search along the linked list. We return a reference to
a weak pointer that was either found on the list or added to it (the result of
applying the |lock| method will tell which case applies, though an expired weak
pointer will behave like a freshly created one, and this is desired); because of
the latter possibility the |classes| field is marked |mutable|, and so |it|
below is a non-|const| iterator.

@< Function definitions @>=
inner_class_weak_ptr& root_datum_value::lookup(const WeightInvolution& delta)
  const
{ auto it = classes.begin();
  while (not classes.at_end(it))
    if (it->first==delta)
      return it->second;
    else ++it;
  classes.insert(it,std::make_pair(delta,inner_class_weak_ptr()));
  return it->second;
}

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

  unsigned int nr = simple_roots->val.n_rows(),
         nc = simple_roots->val.n_columns();

  if (simple_coroots->val.n_rows()!=nr @| or
      simple_coroots->val.n_columns()!=nc)
  { std::ostringstream o;
    o << "Sizes (" << nr << ',' << nc << "),(" @|
      << simple_coroots->val.n_rows() << ','
      << simple_coroots->val.n_columns() @|
      << ") of simple (co)root systems differ";
    throw runtime_error(o.str());
  }
@.Sizes of simple (co)root systems...@>

try @/{
  PreRootDatum prd(simple_roots->val,simple_coroots->val,prefer_coroots);
  prd.test_Cartan_matrix();
  if (l!=expression_base::no_value)
    push_value(root_datum_value::build(std::move(prd)));
}
  catch (error::Cartan_error&)
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
  if (lattice->val.n_rows()!=lattice->val.n_columns() @| or
      lattice->val.n_rows()!=type->val.rank())
  { std::ostringstream o;
    o << "Sub-lattice matrix should have size " @|
      << type->val.rank() << 'x' << type->val.rank();
    throw runtime_error(o.str());
  }
@.Sub-lattice matrix should...@>
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
  if (lattice->val.n_rows()!=r or lattice->val.n_columns()!=r)
  { std::ostringstream o;
    o << "Sub-lattice matrix should have size " @|
      << r << 'x' << r;
    throw runtime_error(o.str());
  }
@.Sub-lattice matrix should...@>

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
  for (unsigned int i=0; i<M->val.n_rows(); ++i)
    if (M->val(i,i)==0) M->val(i,i)=1;
  push_value(std::move(M));
  push_value(whether(prefer_coroots));
@/root_datum_from_type_wrapper(expression_base::single_value);
}


@*2 Functions operating on root data.
%
We start with attributes of root data as a whole.

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
  push_value(std::make_shared<matrix_value>(rd->val.Cartan_matrix()));
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
    push_value(std::make_shared<int_value>(rd->val.semisimple_rank()));
}
@)
void rd_nposroots_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(rd->val.numPosRoots()));
}
@)
void two_rho_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l!=expression_base::no_value)
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value>(rd->val.twoRho()));
}
void two_rho_check_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l!=expression_base::no_value)
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value>(rd->val.dual_twoRho()));
}

@ Here is an auxiliary function that will facilitate the recurring task of
converting a (co)root index from the user range (which has negative values for
negative (co)roots) to the unsigned internal root numbering, while checking that
the index is in the valid range.

@< Local function definitions @>=
RootNbr internal_root_index(const RootDatum& rd, int index, bool is_coroot)
{
  RootNbr npr = rd.numPosRoots();
  RootNbr alpha = npr+index;
  if (alpha>=2*npr)
  { std::ostringstream o;
    o << "Illegal "<< (is_coroot ? "co" : "") << "root index " << index;
    throw runtime_error(o.str());
  }
  return alpha;
}

@ The following functions allow us to look at individual simple roots and simple
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
  RootNbr alpha = internal_root_index(rd->val,root_index,false);
  if (l!=expression_base::no_value)
     push_value(std::make_shared<vector_value>(rd->val.root(alpha)));
}
@)
void coroot_wrapper(expression_base::level l)
{ int coroot_index = get<int_value>()->int_val();
  shared_root_datum rd(get<root_datum_value>());
  RootNbr alpha = internal_root_index(rd->val,coroot_index,true);
  if (l!=expression_base::no_value)
     push_value(std::make_shared<vector_value>(rd->val.coroot(alpha)));
}

@ Also important are look-up functions for roots and coroots.

@< Local function definitions @>=
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

@ The library knows additive decompositions of roots and coroots into simple
roots respectively coroots, and we give the user access to these.

@< Local function definitions @>=
void root_expression_wrapper(expression_base::level l)
{ int root_index = get<int_value>()->int_val();
  shared_root_datum rd = get<root_datum_value>();
  RootNbr alpha = internal_root_index(rd->val,root_index,false);
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value>(rd->val.root_expr(alpha)));
}

void coroot_expression_wrapper(expression_base::level l)
{ int coroot_index = get<int_value>()->int_val();
  shared_root_datum rd = get<root_datum_value>();
  RootNbr alpha = internal_root_index(rd->val,coroot_index,true);
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value>(rd->val.coroot_expr(alpha)));
}

@ Occasionally it is useful to distinguish long and short roots and coroots,
which is a notion relative to the component of the Dynkin diagram the root
belongs to. By convention for simply laced components all roots and coroots are
considered short, while for other components a root is long if and only if the
corresponding coroot is short.

@< Local function def...@>=
void is_long_root_wrapper(expression_base::level l)
{ int root_index = get<int_value>()->int_val();
  shared_root_datum rd = get<root_datum_value>();
  RootNbr alpha = internal_root_index(rd->val,root_index,false);
  if (l!=expression_base::no_value)
    push_value(whether(is_long_root(rd->val,alpha)));
}
@)
void is_long_coroot_wrapper(expression_base::level l)
{ int coroot_index = get<int_value>()->int_val();
  shared_root_datum rd = get<root_datum_value>();
  RootNbr alpha = internal_root_index(rd->val,coroot_index,true);
  if (l!=expression_base::no_value)
    push_value(whether(is_long_coroot(rd->val,alpha)));
}

@ Here we use the built-in method |simple_root_permutation| that, in spite of
its name, works for any positive root index (counting from~$0$).

@< Local function def...@>=
void root_involution_wrapper(expression_base::level l)
{ int root_index = get<int_value>()->int_val();
  shared_root_datum rd = get<root_datum_value>();
  RootNbr alpha = internal_root_index(rd->val,root_index,false);
  if (l==expression_base::no_value)
    return;
  auto perm = rd->val.simple_root_permutation(rd->val.rt_abs(alpha));
  push_value(std::make_shared<vector_value>(perm.begin(),perm.end()));
}

@ An information about roots and coroots that is precomputed in root data and
can be useful for the user tells for each root or coroot $\alpha$ which are the
other roots respectively coroots $\beta$ that are at bottom of a latter
for~$\alpha$, in other words such that $\beta-\alpha$ is not a root. The
following functions extract that information in the form of a |BitMap| (under
the alias |RootNbrSet|), whose coded unsigned values are converted here to the
signed indexing convention of root systems.

@< Local function def...@>=
void root_ladder_bottoms_wrapper(expression_base::level l)
{ int root_index = get<int_value>()->int_val();
  shared_root_datum rd(get<root_datum_value>());
  RootNbr alpha = internal_root_index(rd->val,root_index,false);
  if (l==expression_base::no_value)
    return;

  RootNbr npr = rd->val.numPosRoots();
  const RootNbrSet& bots = rd->val.min_roots_for(alpha);
  own_row result = std::make_shared<row_value>(0);
  result->val.reserve(bots.size());
  for (auto it=bots.begin(); it(); ++it)
    result->val.push_back(std::make_shared<int_value>
      (static_cast<int>(*it-npr)));
  push_value(std::move(result));
}
@)
void coroot_ladder_bottoms_wrapper(expression_base::level l)
{ int coroot_index = get<int_value>()->int_val();
  shared_root_datum rd(get<root_datum_value>());
  RootNbr alpha = internal_root_index(rd->val,coroot_index,true);
  if (l==expression_base::no_value)
    return;

  RootNbr npr = rd->val.numPosRoots();
  const RootNbrSet& bots = rd->val.min_coroots_for(alpha);
  own_row result = std::make_shared<row_value>(0);
  result->val.reserve(bots.size());
  for (auto it=bots.begin(); it(); ++it)
    result->val.push_back(std::make_shared<int_value>
      (static_cast<int>(*it-npr)));
  push_value(std::move(result));
}

@ We provide the fundamental weights and coweights, which are rational
vectors in the span of the roots respectively coroots, such that their pairings
with the simple coroots respectively simple roots are given by the
Kronecker~$\delta$ (in other words, they form dual bases in within those spans).

@< Local function definitions @>=
void fundamental_weight_wrapper(expression_base::level l)
{ int i= get<int_value>()->int_val();
  shared_root_datum rd = get<root_datum_value>();
  if (unsigned(i)>=rd->val.semisimple_rank())
  { std::ostringstream o;
    o << "Invalid index " << i;
    throw runtime_error(o.str());
  }
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value> @|
      (rd->val.fundamental_weight(i)));
}
@)
void fundamental_coweight_wrapper(expression_base::level l)
{ int i= get<int_value>()->int_val();
  shared_root_datum rd = get<root_datum_value>();
  if (unsigned(i)>=rd->val.semisimple_rank())
  { std::ostringstream o;
    o << "Invalid index " << i;
    throw runtime_error(o.str());
  }
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value> @|
      (rd->val.fundamental_coweight(i)));
}

@ We give access to the matrices of all simple or of all positive
(co)roots; for \emph{all} roots we would need index-shifting, so we leave this
to the user program to handle. Our implementation here exploits the by-columns
|matrix_value| constructor, passing iterators to a |int_Matrix| constructor; no
extra copying.

@< Local function definitions @>=
void simple_roots_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<matrix_value> @|
    (rd->val.beginSimpleRoot(),rd->val.endSimpleRoot(),rd->val.rank()));
}
@)
void simple_coroots_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<matrix_value> @|
    (rd->val.beginSimpleCoroot(),rd->val.endSimpleCoroot(),rd->val.rank()));
}
@)
void positive_roots_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<matrix_value> @|
    (rd->val.beginPosRoot(),rd->val.endPosRoot(),rd->val.rank()));
}
@)
void positive_coroots_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<matrix_value> @|
    (rd->val.beginPosCoroot(),rd->val.endPosCoroot(),rd->val.rank()));
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
  std::vector<int_Vector_cref> srl
    (rd->val.beginSimpleRoot(),rd->val.endSimpleRoot());
  srl.insert(srl.end(),rd->val.beginCoradical(),rd->val.endCoradical());
  push_value(std::make_shared<matrix_value> @|
    (srl.begin(),srl.end(),rd->val.rank()));
}
@)
void coroot_radical_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  std::vector<int_Vector_cref> scl
    (rd->val.beginSimpleCoroot(),rd->val.endSimpleCoroot());
  scl.insert(scl.end(),rd->val.beginRadical(),rd->val.endRadical());
  push_value(std::make_shared<matrix_value> @|
    (scl.begin(),scl.end(),rd->val.rank()));
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
  { std::ostringstream o;
    o << "Length " << lambda->val.size() @|
      << " of rational vector differs from rank " << rd->val.rank();
    throw runtime_error(o.str());
  }
@.Length of rational vector...@>
  if (l!=expression_base::no_value)
  @/push_value(root_datum_value::build @|
      (rootdata::integrality_predatum(rd->val,lambda->val)));
}

@ Some function to allow finding the (semisimple) rank and testing dominance for
the integral system without constructing the full integrality datum.

@< Local function definitions @>=
void integrality_rank_wrapper(expression_base::level l)
{ shared_rational_vector lambda = get<rational_vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (lambda->val.size()!=rd->val.rank())
  { std::ostringstream o;
    o << "Length " << lambda->val.size() @|
      << " of rational vector differs from rank " << rd->val.rank();
    throw runtime_error(o.str());
  }
@.Length of rational vector...@>
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>
      (rootdata::integrality_rank(rd->val,lambda->val)));
}
@)
void is_integrally_dominant_wrapper(expression_base::level l)
{ shared_rational_vector lambda = get<rational_vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (lambda->val.size()!=rd->val.rank())
  { std::ostringstream o;
    o << "Length " << lambda->val.size() @|
      << " of rational vector differs from rank " << rd->val.rank();
    throw runtime_error(o.str());
  }
@.Length of rational vector...@>
  if (l==expression_base::no_value)
    return;
  PreRootDatum ipd = rootdata::integrality_predatum(rd->val,lambda->val);
  for (unsigned int j=0; j<ipd.semisimple_rank(); ++j)
    if (ipd.simple_coroot(j).dot(lambda->val.numerator())<0)
    {@;
      push_value(global_false); return;
    }
  push_value(global_true);
}

@ A related function computes a list of fractions of a line segment where the
set of roots with integrality is non-empty.

@< Local function definitions @>=
void integrality_points_wrapper(expression_base::level l)
{ shared_rational_vector lambda = get<rational_vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (lambda->val.size()!=rd->val.rank())
  { std::ostringstream o;
    o << "Length " << lambda->val.size() @|
      << " of rational vector differs from rank " << rd->val.rank();
    throw runtime_error(o.str());
  }
@.Length of rational vector...@>
  if (l==expression_base::no_value)
    return;
@)
  RatNumList ipl = rootdata::integrality_points(rd->val,lambda->val);
    // method normalises rationals
  own_row result = std::make_shared<row_value>(ipl.size());
  for (unsigned int i=0; i<ipl.size(); ++i)
    result->val[i]=std::make_shared<rat_value>(ipl[i]);
  push_value(std::move(result));
}

@ Here are four functions for (hopefully) rapid Weyl group orbit generation,
which can be done using only a root datum. They are implemented by free
functions |Weyl_orbit| and |Weyl_orbit_words| defined in the \.{rootdata}
compilation unit, both in two forms (for weights and coweights) that differ
by the order of their arguments. The first of these pairs return an orbits as a
(usually very wide) |int_Matrix|, the second as a |sl_list<WeylElt>|, the action
of whose elements on the original vector will produce the orbit. In the exported
built-in functions the distinction between the two variants will in each case
again be determined by the order of the arguments.

@< Local function definitions @>=
void Weyl_orbit_wrapper(expression_base::level l)
{
  shared_vector v = get<vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<matrix_value>(Weyl_orbit(rd->val,v->val)));
}
@)
void Weyl_coorbit_wrapper(expression_base::level l)
{
  shared_root_datum rd = get<root_datum_value>();
  shared_vector v = get<vector_value>();
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<matrix_value>(Weyl_orbit(v->val,rd->val)));
}

void Weyl_orbit_ws_wrapper(expression_base::level l)
{
  shared_vector v = get<vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (l==expression_base::no_value)
    return;
@)
  const auto& ws = Weyl_orbit_words(rd->val, rd->W(), v->val);
  own_row result = std::make_shared<row_value>(ws.size());
  size_t i=0;
  for (auto&& w : ws)
    result->val[i++]=std::make_shared<W_elt_value>(rd,w);
  push_value(std::move(result));
}
@)
void Weyl_coorbit_ws_wrapper(expression_base::level l)
{
  shared_root_datum rd = get<root_datum_value>();
  shared_vector v = get<vector_value>();
  if (l==expression_base::no_value)
    return;
@)
  const auto& ws = Weyl_orbit_words(v->val, rd->val, rd->W());
  own_row result = std::make_shared<row_value>(ws.size());
  size_t i=0;
  for (auto&& w : ws)
    result->val[i++]=std::make_shared<W_elt_value>(rd,w);
  push_value(std::move(result));
}

@ Here are two functions making available to the user the |wall_set| and
|alcove_center| functions defined in \.{alcoves.cpp}, which serve to improve the
efficiency and safety against rational number overflow of the deformation
algorithm.

@h "alcoves.h"

@< Local function def...@>=
void walls_wrapper(expression_base::level l)
{
  shared_rational_vector gamma = get<rational_vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (gamma->val.size()!=rd->val.rank())
  { std::ostringstream o;
    o << "Rational weight size mismatch: "
      << gamma->val.size() << ':' << rd->val.rank();
    throw runtime_error(o.str());
  }
  if (l==expression_base::no_value)
    return;
@)
  int npr = rd->val.numPosRoots();
  RootNbrSet integrals,walls = weyl::wall_set(rd->val,gamma->val,integrals);
  own_row roots = std::make_shared<row_value>(0);
  roots->val.reserve(walls.size());
  for (auto it=integrals.begin(); it(); ++it)
  {
    assert(*it-npr < static_cast<unsigned>(npr)); // must be positive root
    roots->val.push_back(std::make_shared<int_value>(*it-npr));
  }

  auto non_integrals = weyl::sorted_by_label(rd->val,walls,integrals);
  for (auto it=non_integrals.begin(); not non_integrals.at_end(it); ++it)
    roots->val.push_back(std::make_shared<int_value>
    // convert to signed root index
       (static_cast<int>(*it)-npr));
  push_value(std::move(roots));
  push_value(std::make_shared<int_value>(integrals.size()));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}
@)
void alcove_center_wrapper(expression_base::level l)
{
  shared_module_parameter p = get<module_parameter_value>();
  if (l==expression_base::no_value)
    return;

  push_value(std::make_shared<module_parameter_value> @|
    (p->rf,weyl::alcove_center(p->rc(),p->val)));
}

@ One interesting property of alcoves is that (projected to the rational span of
the root system) they have exactly one vertex that is in the root lattice.
@< Local function definitions @>=
void alcove_root_vertex_wrapper (expression_base::level l)
{
  shared_rational_vector gamma = get<rational_vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (gamma->val.size()!=rd->val.rank())
  { std::ostringstream o;
    o << "Rational weight size mismatch: "
      << gamma->val.size() << ':' << rd->val.rank();
    throw runtime_error(o.str());
  }
  if (l==expression_base::no_value)
    return;
  push_value(std::make_shared<vector_value> @|
    (weyl::root_vertex_of_alcove(rd->val,gamma->val)));
}

@ Here are more general and more efficient functions for orbit generation, in
fact enumerating quotients of pseudo-Levi subgroups.

@< Local function definitions @>=
void basic_orbit_ws_wrapper(expression_base::level l)
{
  unsigned int stab_rank = get<int_value>()->uint_val();
  shared_row v = get<row_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (v->val.size()<=stab_rank)
    throw runtime_error("Index too large for given list of root numbers");
  RootNbrSet stab; RootNbr final;
  @< Check validity of root indices in |v->val| and the absence acute angles
     among them; store internal number for its entry |stab_rank| in |final|,
     and the set of its earlier entries in |stab| @>
  if (l==expression_base::no_value)
    return;
@)
  bool to_affine_orbit=false;
  @< Set |to_affine_orbit| to whether coroot |final| is linearly dependent
     on the coroots in |stab|; also intersect |stab| with the set of roots in
     the Dynkin diagram component of |final| @>
  const auto ws = to_affine_orbit
    ? weyl::complete_affine_component(rd->val,rd->W(),stab,final)
    : weyl::finite_subquotient(rd->val,rd->W(),stab,final);
  own_row result = std::make_shared<row_value>(ws.size());
  { size_t i=0;
    for (auto&& w : ws)
      result->val[i++]=std::make_shared<W_elt_value>(rd,w);
  }
  push_value(std::move(result));
}
@)
void affine_orbit_ws_wrapper(expression_base::level l)
{
  shared_rational_vector gamma = get<rational_vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (rd->val.rank()!=gamma->val.size())
  { std::ostringstream o;
    o << "Rank and rational weight size mismatch " @|
      << rd->val.rank() << ':' << gamma->val.size();
    throw runtime_error(o.str());
  }
  if (l==expression_base::no_value)
    return;
@)
  auto ws = weyl::affine_orbit_ws(rd->val,rd->W(),gamma->val);
  own_row result = std::make_shared<row_value>(ws.size());
  { size_t i=0;
    for (auto&& w : ws)
      result->val[i++]=std::make_shared<W_elt_value>(rd,std::move(w));
  }
  push_value(std::move(result));
}

@ We temporarily make a list |walls| of root numbers after conversion to
internal root numbering, but in the end only retain in |stab| the set of the
first |stab_rank| entries, and in |final| the one following it.
@< Check validity of root indices in |v->val| and... @>=
{
  RootNbrList walls; walls.reserve(v->val.size());
  const unsigned int npr = rd->val.numPosRoots();
  for (const auto& entry : v->val)
  { int r = force<int_value>(entry.get())->int_val();
      // could be positive or negative
    if (r+npr < rd->val.numRoots()) // unsigned comparison
      walls.push_back(r+npr); // convert to internal index
    else
    { std::ostringstream o;
      o << "Invalid root number " << r;
      throw runtime_error(o.str());
    }
  }

  for (const RootNbr& alpha : walls)
    for (const RootNbr& beta : walls)
      if (&alpha!=&beta and rd->val.bracket(alpha,beta)>0)
      { std::ostringstream o;
        o << "Roots " << int(alpha-npr) << " and " << int(beta-npr)
          << " have acute angle.";
        throw runtime_error(o.str());
      }
  stab = RootNbrSet(rd->val.numRoots(),&walls[0],&walls[stab_rank]);
  final = walls[stab_rank];
}

@ When calling |affine_orbit_ws|, we need to decide whether the internal
function to call is going to be |finite_subquotient| or
|complete_affine_components|; the first one does not use modular vector
arithmetic because the root |final| extends the finite Coxeter group |stab| to a
higher rank finite Coxeter group, but for the second one this extension would be
to an affine (therefore infinite) Coxeter group and this needs to be reduced to
a finite computation by working modulo the root lattice.

@< Set |to_affine_orbit| to whether coroot |final| is linearly dependent
     on the coroots in |stab|; also intersect |stab| with the set of roots in
     the Dynkin diagram component of |final| @>=
{
  RootNbrSet S = stab; S.insert(final);
  for (const auto& comp : rootdata::components(rd->val,S))
    if (comp.isMember(final))
    {
      stab &= comp;
      RootNbrList roots(comp.begin(),comp.end());
      auto dependency = // the number of dependencies among coroots in |comp|
        lattice::kernel(rd->val.Cartan_matrix(roots)).n_columns();
      assert(dependency <= 1);
      to_affine_orbit = dependency>0;
      break;
    }
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

install_function(root_datum_eq_wrapper,@|"=","(RootDatum,RootDatum->bool)");
install_function(root_datum_neq_wrapper,@|"!=","(RootDatum,RootDatum->bool)");
install_function(datum_Cartan_wrapper,@|"Cartan_matrix","(RootDatum->mat)");
install_function(rd_rank_wrapper,@|"rank","(RootDatum->int)");
install_function(rd_semisimple_rank_wrapper,@|
		 "semisimple_rank","(RootDatum->int)");
install_function(rd_nposroots_wrapper@|,"nr_of_posroots","(RootDatum->int)");
install_function(two_rho_wrapper@|,"two_rho","(RootDatum->vec)");
install_function(two_rho_check_wrapper@|,"two_rho_check","(RootDatum->vec)");

install_function(root_wrapper,@|"root","(RootDatum,int->vec)");
install_function(coroot_wrapper,@|"coroot","(RootDatum,int->vec)");
install_function(root_index_wrapper@|,"root_index","(RootDatum,vec->int)");
install_function(coroot_index_wrapper@|,"coroot_index","(RootDatum,vec->int)");
install_function(root_expression_wrapper@|,"root_expression"
		,"(RootDatum,int->vec)");
install_function(coroot_expression_wrapper@|,"coroot_expression"
		,"(RootDatum,int->vec)");
install_function(is_long_root_wrapper@|,"is_long_root" ,"(RootDatum,int->bool)");
install_function(is_long_coroot_wrapper@|,"is_long_coroot"
		,"(RootDatum,int->bool)");
install_function(root_involution_wrapper@|,"root_involution"
		,"(RootDatum,int->vec)");
install_function(root_ladder_bottoms_wrapper,@|"root_ladder_bottoms"
                ,"(RootDatum,int->[int])");
install_function(coroot_ladder_bottoms_wrapper,@|"coroot_ladder_bottoms"
                ,"(RootDatum,int->[int])");
install_function(fundamental_weight_wrapper,@|
		 "fundamental_weight","(RootDatum,int->ratvec)");
install_function(fundamental_coweight_wrapper,@|
		 "fundamental_coweight","(RootDatum,int->ratvec)");

install_function(simple_roots_wrapper,@|"simple_roots","(RootDatum->mat)");
install_function(simple_coroots_wrapper,@|"simple_coroots","(RootDatum->mat)");
install_function(positive_roots_wrapper,@| "posroots","(RootDatum->mat)");
install_function(positive_coroots_wrapper,@| "poscoroots","(RootDatum->mat)");
install_function(root_coradical_wrapper,@|"root_coradical","(RootDatum->mat)");
install_function(coroot_radical_wrapper,@|"coroot_radical","(RootDatum->mat)");

install_function(dual_datum_wrapper,@|"dual","(RootDatum->RootDatum)");
install_function(derived_info_wrapper,@|
		 "derived_info","(RootDatum->RootDatum,mat)");
install_function(mod_central_torus_info_wrapper,@|
		 "mod_central_torus_info","(RootDatum->RootDatum,mat)");
install_function(integrality_datum_wrapper,@|
                 "integrality_datum","(RootDatum,ratvec->RootDatum)");
install_function(integrality_rank_wrapper,@|
                 "integrality_rank","(RootDatum,ratvec->int)");
install_function(is_integrally_dominant_wrapper,@|
                 "is_integrally_dominant","(RootDatum,ratvec->bool)");
install_function(integrality_points_wrapper,@|
                 "integrality_points","(RootDatum,ratvec->[rat])");

install_function(Weyl_orbit_wrapper@|,"Weyl_orbit","(RootDatum,vec->mat)");
install_function(Weyl_coorbit_wrapper@|,"Weyl_orbit","(vec,RootDatum->mat)");
install_function(Weyl_orbit_ws_wrapper@|,"Weyl_orbit_ws",
		"(RootDatum,vec->[WeylElt])");
install_function(Weyl_coorbit_ws_wrapper@|,"Weyl_orbit_ws",
		"(vec,RootDatum->[WeylElt])");
install_function(walls_wrapper,"walls","(RootDatum,ratvec->[int],int)");
install_function(alcove_center_wrapper,"alcove_center","(Param->Param)");
install_function(alcove_root_vertex_wrapper@|,"alcove_root_vertex",
		"(RootDatum,ratvec->vec)");
install_function(basic_orbit_ws_wrapper@|,"basic_orbit_ws",
		"(RootDatum,[int],int->[WeylElt])");
install_function(affine_orbit_ws_wrapper@|,"affine_orbit_ws",
		"(RootDatum,ratvec->[WeylElt])");


@*1 Weyl group elements.
%
Weyl groups were implemented by Fokko du Cloux using the transducer model he
developed, and have been used since the beginnings of Atlas for the
representation of involutions. This only exposes the possibilities of the
implementation in a very limited way, and it is useful for \.{atlas} users to
have direct access to computations with Weyl group elements.

Since the |WeylElt| calls has non remote (outside the object proper) data, as
would be the case had it used |std::vector|, it does not have efficient move
semantics, and there is no point providing a constructor taking |WeylElt| by
rvalue-reference.

@<Type definitions @>=
struct W_elt_value : public value_base
{ shared_root_datum rd;
  const WeylGroup& W;
  WeylElt val;
@)
  W_elt_value(const shared_root_datum& rd, const WeylElt& w)
@/: rd(rd), W(rd->W()), val(w) @+{}
@)
  virtual void print(std::ostream& out) const;
  static const char* name() @+{@; return "Weyl group element"; }
  W_elt_value (const W_elt_value& ) = default;
    // we use |get_own<W_elt_value>|

};
typedef std::shared_ptr<const W_elt_value> shared_W_elt;
typedef std::shared_ptr<W_elt_value> own_W_elt;

@ To distinguish Weyl group elements from other values, we use a format like
\.{<3.2.0.1.2.1.0>}.
@< Function def...@>=
void W_elt_value::print(std::ostream& out) const
{ WeylWord ww = W.word(val); out << '<';
  for (auto it=ww.begin(); it!=ww.end(); ++it)
    out << (it==ww.begin() ? "" : ".") << static_cast<unsigned int>(*it);
  out << '>';
}

@*2 Functions defined for Weyl group elements.
%
A first basic way to make a Weyl group element will be to specify a root datum
and a word in the simple reflections, each represented by the index of a simple
root. We need to test that the given word is valid for the semisimple rank of
the root datum, and since there will be more functions that take such a Weyl
word as an argument, we define a separate function for performing the test.

@< Local function def...@>=
WeylWord check_Weyl_word(const row_value& r, unsigned int rank)
{ WeylWord result; result.reserve(r.val.size());
  for (auto it=r.val.begin(); it!=r.val.end(); ++it)
  { auto v=force<int_value>(it->get())->val.ulong_val();
    if (v<rank)
      result.push_back(v);
    else
  { std::ostringstream o;
    o << "Illegal Weyl word entry " << v @|
         << " (should be <" << rank << ')';
    throw runtime_error(o.str());
  }
@.Illegal Weyl word entry@>
  }
  return result;
}
@)
void W_elt_wrapper(expression_base::level l)
{ shared_row r = get<row_value>();
  shared_root_datum rd = get<root_datum_value>();
  auto ww=check_Weyl_word(*r,rd->val.semisimple_rank());
  if (l!=expression_base::no_value)
    push_value(std::make_shared<W_elt_value>(rd,rd->W().element(ww)));
}

@ We also want to be able to convert a Weyl group element back to a word, to
extract the root datum it was built from, to find its length, to and test for
equality.

@< Local function def... @>=
void W_word_wrapper(expression_base::level l)
{ shared_W_elt w = get<W_elt_value>();
  if (l==expression_base::no_value)
    return;
@)
  auto word = w->W.word(w->val);
  int_Vector wv (word.begin(),word.end()); // convert
  push_value(vector_to_row(wv));
}
@)
void datum_from_W_elt_wrapper(expression_base::level l)
{ shared_W_elt w = get<W_elt_value>();
  if (l!=expression_base::no_value)
    push_value(w->rd);
}
@)
void W_length_wrapper(expression_base::level l)
{ shared_W_elt w = get<W_elt_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(w->W.length(w->val)));
}
@)
void W_elt_unary_eq_wrapper
(expression_base::level l)
{ shared_W_elt w = get<W_elt_value>();
  if (l!=expression_base::no_value)
    push_value(whether(w->val==WeylElt()));
}
void W_elt_unary_neq_wrapper
(expression_base::level l)
{ shared_W_elt w = get<W_elt_value>();
  if (l!=expression_base::no_value)
    push_value(whether(w->val!=WeylElt()));
}
@)
void W_elt_eq_wrapper
(expression_base::level l)
{ shared_W_elt w1 = get<W_elt_value>();
  shared_W_elt w0 = get<W_elt_value>();
  if (&w0->W!=&w1->W)
    throw runtime_error("Weyl group mismatch");
@.Weyl group mismatch@>
  if (l!=expression_base::no_value)
    push_value(whether(w0->val==w1->val));
}
void W_elt_neq_wrapper
(expression_base::level l)
{ shared_W_elt w1 = get<W_elt_value>();
  shared_W_elt w0 = get<W_elt_value>();
  if (&w0->W!=&w1->W)
    throw runtime_error("Weyl group mismatch");
@.Weyl group mismatch@>
  if (l!=expression_base::no_value)
    push_value(whether(w0->val!=w1->val));
}

@ An obvious thing we want to do with Weyl group elements is multiply them, or
take their inverse.

@< Local function def... @>=
void W_elt_prod_wrapper(expression_base::level l)
{ shared_W_elt w1 = get<W_elt_value>();
  own_W_elt w0 = get_own<W_elt_value>();
  if (&w0->W!=&w1->W)
    throw runtime_error("Weyl group mismatch");
@.Weyl group mismatch@>
  if (l==expression_base::no_value)
    return;
@)
  w0->W.mult(w0->val,w1->val);
  push_value(std::move(w0));
}

void W_elt_invert_wrapper(expression_base::level l)
{ shared_W_elt w = get<W_elt_value>();
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<W_elt_value>(w->rd,w->W.inverse(w->val)));
}

@ Since only one argument needs to furnish the Weyl group, we also define
multiplication on either side by a generator represented by its index, or by a
complete Weyl word.

@< Local function def... @>=
void check_Weyl_gen(int s, unsigned int rank)
{ if (static_cast<unsigned int>(s)>=rank)
  { std::ostringstream o;
    o << "Generator " << s @|
      << " out of range for Weyl group (should be <" << rank << ')';
    throw runtime_error(o.str());
  }
}
@)
void W_elt_gen_prod_wrapper(expression_base::level l)
{ auto s = get<int_value>()->int_val();
  own_W_elt w = get_own<W_elt_value>();
  check_Weyl_gen(s,w->W.rank());
  if (l==expression_base::no_value)
    return;
@)
  w->W.mult(w->val,s);
  push_value(std::move(w));
}
@)
void W_gen_elt_prod_wrapper(expression_base::level l)
{ own_W_elt w = get_own<W_elt_value>();
  auto s = get<int_value>()->int_val();
  check_Weyl_gen(s,w->W.rank());
  if (l==expression_base::no_value)
    return;
@)
  w->W.left_multiply(w->val,s);
  push_value(std::move(w));
}

void W_elt_word_prod_wrapper(expression_base::level l)
{ shared_row r = get<row_value>();
  own_W_elt w = get_own<W_elt_value>();
  auto ww=check_Weyl_word(*r,w->W.rank());
  if (l==expression_base::no_value)
    return;
@)
  w->W.mult(w->val,ww);
  push_value(std::move(w));
}
@)
void W_word_elt_prod_wrapper(expression_base::level l)
{ own_W_elt w = get_own<W_elt_value>();
  shared_row r = get<row_value>();
  auto ww=check_Weyl_word(*r,w->W.rank());
  if (l==expression_base::no_value)
    return;
@)
  w->W.left_multiply(w->val,ww);
  push_value(std::move(w));
}

@ Since we have associated the Weyl group to a root datum, we can let Weyl
group elements act on weights and on coweights, and we can convert between them
and (appropriate) square matrices.

@< Local function def... @>=
void W_elt_weight_prod_wrapper(expression_base::level l)
{ shared_vector v = get<vector_value>();
  shared_W_elt w = get<W_elt_value>();
  const auto& rd = w->rd->val;
  if (v->val.size()!=rd.rank())
  { std::ostringstream o;
    o << "Rank and weight size mismatch " @|
      << rd.rank() << ':' << v->val.size();
    throw runtime_error(o.str());
  }
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value>(w->W.image_by(rd,w->val,v->val)));
}
void coweight_W_elt_prod_wrapper(expression_base::level l)
{ shared_W_elt w = get<W_elt_value>();
  shared_vector v = get<vector_value>();
  const auto& rd = w->rd->val;
  if (v->val.size()!=rd.rank())
  { std::ostringstream o;
    o << "Coweight size and rank mismatch " @|
      <<  v->val.size() << ':' << rd.rank();
    throw runtime_error(o.str());
  }
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value>(w->W.image_by(rd,v->val,w->val)));
}

@ A converse operation is to decompose a weight into a Weyl group element and a
dominant weight. There being no value to represent the Weyl group as such, this
function takes a root datum as first argument. We also provide a version for
coweights, which will return a dominant coweight and a Weyl group element to be
applied to its right to obtain the original coweight; to remind of the dual
nature, and that the Weyl word has an opposite interpretation to the first case,
we inverse argument and result order here.

@< Local function def... @>=
void W_decompose_wrapper(expression_base::level l)
{ own_vector v = get_own<vector_value>();
  shared_root_datum rd = get<root_datum_value>();
  if (v->val.size()!=rd->val.rank())
  { std::ostringstream o;
    o << "Rank and weight size mismatch " @|
      << rd->val.rank() << ':' << v->val.size();
    throw runtime_error(o.str());
  }
  if (l==expression_base::no_value)
    return;
@)
  auto ww = rd->val.factor_dominant(v->val);
@/push_value(std::make_shared<W_elt_value>(rd,rd->W().element(ww)));
  push_value(v);
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}
@)
void W_codecompose_wrapper(expression_base::level l)
{ shared_root_datum rd = get<root_datum_value>();
  own_vector v = get_own<vector_value>();
  if (v->val.size()!=rd->val.rank())
  { std::ostringstream o;
    o << "Coweight size and rank mismatch " @|
      <<  v->val.size() << ':' << rd->val.rank();
    throw runtime_error(o.str());
  }
  if (l==expression_base::no_value)
    return;
@)
  auto ww = rd->val.factor_codominant(v->val);
@/push_value(v);
  push_value(std::make_shared<W_elt_value>(rd,rd->W().element(ww)));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ It is easy and efficient to compute the action of Weyl group elements on
roots, thanks to the method |RootSystem::permuted_root|, which taps directly
into the stored tables of root permutations for (simple) root reflections.

@< Local function def...@>=
void root_permutation_wrapper(expression_base::level l)
{ shared_W_elt w = get<W_elt_value>();
  const RootSystem& rs = w->rd->val;
  if (l==expression_base::no_value)
    return;
  const auto ww = w->W.word(w->val);
  int_Vector result(rs.numRoots());
  for (unsigned i=0; i<rs.numPosRoots(); ++i)
  { const auto alpha = rs.posRootNbr(i);
    result[alpha] = rs.permuted_root(ww,alpha);
    result[rs.negRootNbr(i)] = rs.rootMinus(result[alpha]);
  }
  push_value(std::make_shared<vector_value>(std::move(result)));
}


@ Finally we install everything related to Weyl groups elements. Note that since
Weyl words and (co)weights have the same \.{atlas} type \&{vec}, we need to make
a choice whether to associate the wrapper function |W_elt_word_prod_wrapper| or
|W_elt_weight_prod_wrapper| with the \.* operator name, and similarly for
|W_word_elt_prod_wrapper| and |coweight_W_elt_prod_wrapper|; we choose the
latter in both cases, and name the former \.{\#\#} similarly to the
concatenation operation.

@< Install wrapper functions @>=
install_function(W_elt_wrapper,"W_elt","(RootDatum,[int]->WeylElt)");
install_function(W_word_wrapper,"word","(WeylElt->[int])");
install_function(datum_from_W_elt_wrapper,"root_datum","(WeylElt->RootDatum)");
install_function(W_length_wrapper,"length","(WeylElt->int)");
install_function(W_elt_unary_eq_wrapper,"=","(WeylElt->bool)");
install_function(W_elt_unary_neq_wrapper,"!=","(WeylElt->bool)");
install_function(W_elt_eq_wrapper,"=","(WeylElt,WeylElt->bool)");
install_function(W_elt_neq_wrapper,"!=","(WeylElt,WeylElt->bool)");
install_function(W_elt_prod_wrapper,"*","(WeylElt,WeylElt->WeylElt)");
install_function(W_elt_invert_wrapper,"/","(WeylElt->WeylElt)");
install_function(W_elt_gen_prod_wrapper,"#","(WeylElt,int->WeylElt)");
install_function(W_gen_elt_prod_wrapper,"#","(int,WeylElt->WeylElt)");
install_function(W_elt_word_prod_wrapper,"##","(WeylElt,[int]->WeylElt)");
install_function(W_word_elt_prod_wrapper,"##","([int],WeylElt->WeylElt)");
install_function(W_elt_weight_prod_wrapper,"*","(WeylElt,vec->vec)");
install_function(coweight_W_elt_prod_wrapper,"*","(vec,WeylElt->vec)");
install_function(W_decompose_wrapper,@|
                 "from_dominant","(RootDatum,vec->WeylElt,vec)");
install_function(W_codecompose_wrapper,@|
                 "from_dominant","(vec,RootDatum->vec,WeylElt)");
install_function(root_permutation_wrapper,"root_permutation","(WeylElt->vec)");


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
  const auto r = delta.n_rows();
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
{ if (delta.n_rows()!=r or delta.n_columns()!=r)
  { std::ostringstream o;
     o << "Involution should be a " @|
       << r << 'x' << r << " matrix;"
     @| " received a "  << delta.n_rows() << 'x' << delta.n_columns()
      << " matrix";
    throw runtime_error(o.str());
  }
@.Involution should be...@>
  if (not (delta*delta-1).is_zero())
    throw runtime_error("Given transformation is not an involution");
@.Given transformation...@>
}

@ We now come to the part of the analysis that involves the root datum. Since
the matrix is already expressed on the basis of the weight lattice used by the
root datum, the question of stabilising that lattice is settled, but we must
check that the matrix is indeed an involution (for which we reuse the above
module), and that it gives an automorphism of the root datum. Although we need
an automorphism of the \emph{based} root datum, we allow any root datum
automorphism to be supplied; by composing with an element of $W$ we can then
ensure that all positive roots map to positive roots. The following auxiliary
function checks whether |delta| is an involution mapping simple roots to roots
and simple coroots to coroots (throwing an error if it not), and returns a list
of images of simple roots (useful for transforming |delta| into a based root
datum involution).

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
{ const unsigned int r=rd.rank(), s=rd.semisimple_rank();
  @< Check that |delta| is an $r\times{r}$ matrix defining an involution @>
  RootNbrList Delta(s);
  for (weyl::Generator i=0; i<s; ++i)
  { Delta[i]=rd.root_index(delta*rd.simpleRoot(i));
    if (Delta[i]==rd.numRoots()) // then image not found
    { std::ostringstream o;
      o << "Matrix maps simple root " << (unsigned)i @|
						      << " to non-root";
      throw runtime_error(o.str());
    }
    if (delta.right_prod(rd.simpleCoroot(i))!=rd.coroot(Delta[i]))
    { std::ostringstream o;
      o << "Matrix does not map simple coroot "
        << (unsigned)i @| << " to coroot " << Delta[i]-rd.numPosRoots();
      throw runtime_error(o.str());
    }
  }
  return Delta;
}
@)
void check_based_root_datum_involution
  (const RootDatum& rd, const WeightInvolution& delta)
{ const auto Delta=check_root_datum_involution(rd,delta);
  const auto s=rd.semisimple_rank();
  for (weyl::Generator i=0; i<s; ++i)
    if (not rd.is_simple_root(Delta[i]))
      throw runtime_error("Root datum involution is not distinguished");
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

@< Local function def...@>=
weyl::Twist check_involution
 (WeightInvolution& delta, const RootDatum& rd,
  WeylWord& ww, @| lietype::Layout* lo=nullptr)
{ RootNbrList Delta = check_root_datum_involution(rd,delta);
  ww = wrt_distinguished(rd,Delta);
  const unsigned int r=rd.rank(), s=rd.semisimple_rank();
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
  for (unsigned int i=0; i<ww.size(); ++i)
    // apply generators in (original) |theta| towards |delta| order
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
{ DynkinDiagram diagram(rd.Cartan_matrix());
  auto comps = diagram.components(); // connected components
  type = diagram.type();
  pi = diagram.perm(); // |pi| normalises to Bourbaki order
  assert(type.size()==comps.size());
@)
  inner_class.reserve(comps.size()+r-s); // certainly enough
  unsigned int offset=0; // accumulated rank of simple factors seen

  unsigned int i=0; // index into |type|
  for (auto cit=comps.begin(); not comps.at_end(cit); ++cit,++i)
  { bool equal_rank=true;
    unsigned int comp_rank = cit->rank();
    assert (comp_rank==type[i].rank());
       // and |*it| runs through bits in |*cit| in following loop
    for (auto it=&pi[offset]; it!=&pi[offset+comp_rank]; ++it)
      if (p[*it]!=*it) {@; equal_rank=false; break; }
@)  if (equal_rank) inner_class.push_back('c');
      // identity on this component: compact component
    else if (cit->support.test(p[pi[offset]]))
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
  } // |for (cit;;)|
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
  for (unsigned int j=offset; j<offset+comp_rank; ++j)
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
    if (cit1->support.test(beta))
      break;
  if (comps.at_end(cit1))
    throw logic_error("Non matching Complex factor");
@.Non matching Complex factor@>

#ifndef NDEBUG
  assert(type[k]==type[i]); // paired simple types for complex factor
  for (unsigned int l=1; l<comp_rank; ++l)
    assert(cit1->support.test(p[pi[offset+l]]));
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
{ for (unsigned int k=0; k<r-s; ++k)
    type.push_back(SimpleLieType('T',1));
  int_Matrix root_lattice
    (rd.beginSimpleRoot(),rd.endSimpleRoot(),r,tags::IteratorTag());
@/CoeffList factor; // values will be unused
  int_Matrix basis =
     matreduc::adapted_basis(root_lattice,factor);
@/WeightInvolution inv = // involution on quotient by root lattice
     basis.inverse().block(s,0,r,r)*delta*basis.block(0,s,r,r);
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

@~The class |inner_class_value| was the first Atlas type where we deviated from
the previously used scheme of holding an Atlas library object with the main
value in a data member |val|. That was to adapt to the fact that |InnerClass| is
not copy-constructible (its copy constructor is deleted), but since we have
eliminated any need to copy an |inner_class_value| as well, we return to the
situation where the main |InnerClass| object is a member of |inner_class_value|.
Although we handle them using |shared_inner_class|, which is a (smart) pointer
to constant, we occasionally need non-|const| access to the |InnerClass|, namely
to expand the contained table of involutions when required. This is just a form
of lazy initialisation, so does not violate the conceptual constant nature of
the |inner_class_value|; therefore we declare the |val| member to be |mutable|.
When |val| was still a reference, the same effect was produced because \Cpp\
allows a stored non-|const| reference to be exported from a |const| object; the
|mutable| solution is a bit more explicit.

An important property that our approach is set up to ensure, is that we can
compare identity of inner class values by comparing their addresses (checking if
they are the same object internally). Namely, we ensure there cannot coexist two
identical copies of an inner class (similarly to what we did for
|root_datum_value|, and with an implementation depending on that class).

The main constructor takes sufficient data to immediately call the |InnerClass|
constructor, but it will always be invoked indirectly via the |build| static
method; this allows doing some preparations, possibly exporting a |WeylWord| for
the transformation to make the given involution distinguished, and returning a
|shared_inner_class| value ready for use. As in the case of |root_datum_value|,
we need the constructor to be public because it will be called from within
|std::make_shared|, but we include a special |token| argument that only methods
of our |inner_class_value| can supply, so as to make circumventing |build|
impossible.

@< Type definitions @>=
typedef std::shared_ptr<const inner_class_value> shared_inner_class;
class real_form_value;
typedef std::shared_ptr<const real_form_value> shared_real_form;
typedef std::weak_ptr<const real_form_value> real_form_weak_ptr;
@)
class inner_class_value : public value_base
{ struct token@+{}; // serves to control access to the constructor
public:
  mutable InnerClass val; // our |InnerClass| proper is stored here
  shared_root_datum datum, dual_datum;
    // share, our stored |InnerClass| depends on these
  lietype::LieType rd_type;
@+lietype::InnerClassType ic_type; // types used for printing only
  const output::FormNumberMap interface,dual_interface;
@)
  inner_class_value (shared_root_datum rd, shared_root_datum dual_rd, @|
   const WeightInvolution& tau, @|
   const lietype::Layout& lo, token);
static shared_inner_class build
  (shared_root_datum srd, WeightInvolution& tau, WeylWord* wp=nullptr);
@)
  virtual void print(std::ostream& out) const;
  static const char* name() @+{@; return "inner class"; }
  inner_class_value (const inner_class_value& ) = delete;
@)
  shared_inner_class dual () const;
@) // tables for |real_form_value| identification:
  mutable std::vector<real_form_weak_ptr> real_form_wptr;
};

@ The constructor used to be complicated by the storage of a dual |InnerClass|
field and a reference count (mainly to ensure the double dual of an inner class
would test equal to it), but that has been made unnecessary and the constructor
is quite simple. The |Layout| argument is used for convenience (it groups
information of various nature gathered by |check_involution|), and is used to
initialise the |rd_type|, |ic_type|, |interface| and |dual_interface| fields.

@< Function def...@>=
inner_class_value::inner_class_value
  (shared_root_datum rd, shared_root_datum dual_rd, const WeightInvolution& tau,
   const lietype::Layout& lo,
   token)
@/: val(rd->val,dual_rd->val,tau)
   // contain |InnerClass|, referring to |RootDatum| pair
@/, datum(std::move(rd)), dual_datum(std::move(dual_rd))
    // keep the passed shared pointers
@/, rd_type(lo.d_type), ic_type(lo.d_inner)
, interface(val,lo), dual_interface(val,lo,tags::DualTag())
@/,real_form_wptr(val.numRealForms())
{}

@ Since in the \.{atlas} program we always construct inner classes from existing
root data, we want an |inner_class_value| to use an |InnerClass| that refers to
rather than stores a root datum and its dual. To that end the |build| method
prepares two |shared_root_datum| pointers (the logic of that class will ensure
sharing with existing values is done whenever possible), and passes them to the
constructor for storage within in our |inner_class_value|; the |InnerClass|
constructor will capture |RootDatum| references from them.

Matching identical inner classes, namely using the same |root_datum_value| and
equal (distinguished) |WeightInvolution| values, is done within the
|root_datum_value::lookup| method. All we need to do here is see if that returns
a usable pointer to an existing |inner_class_value| (the weak pointer |lock|
method decides this) and if so return a |shared_inner_class| for it; if not,
then we construct an |inner_class_value| here, and before returning store a weak
pointer for it in the location indicated by |lookup|.

@< Function def...@>=
shared_inner_class inner_class_value::build
  (shared_root_datum srd, WeightInvolution& tau, WeylWord* wp)
{ WeylWord ww;
  lietype::Layout lo;
  if (wp==nullptr)
    wp=&ww; // use |ww| as dummy output unless |wp| points somewhere
  check_involution(tau,srd->val,*wp,&lo);
    // may also modify |tau|, and sets |lo|
@)
  auto& w_ptr = srd->lookup(tau);
  if (auto p = w_ptr.lock())
    return p; // reuse existing inner class, if found in |srd|
  auto result =
    std::make_shared<inner_class_value>(srd,srd->dual(),tau,lo,token());
  w_ptr = result; // store a weak pointer version inside |srd| table
  return result; // return pointer to modified instance of |inner_class_value|
}

@ Constructing a dual |inner_class_value| to a given one, is fairly easy, since
|build| will do the necessary to ensure identification of isomorphic inner
classes, notably of the double dual with the inner class itself. The only
subtlety here is that |build| will potentially modify its involution argument
to make it distinguished (though it will not actually here), whence we cannot
pass the constant reference |val.dualDistinguished()| directly here, but need to
store it in a local variable.

@< Function def...@>=
shared_inner_class inner_class_value::dual () const
{ auto delta_for_dual =
    val.dualDistinguished(); // we need a non-|const| reference
  return inner_class_value::build(dual_datum,delta_for_dual);
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
    push_value(std::move(ic_value));
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
  push_value(std::make_shared<W_elt_value>(rd,rd->W().element(ww)));
  push_value(std::move(ic_value));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@*2 Functions operating on inner classes.
%
Here are our first functions that access and operate on values of type
|inner_class_value|. They test for equality (as mentioned above this can be done
by simple pointer comparison), recover the ingredients that were used in the
construction, and the final one constructs the dual inner class.

@< Local function def...@>=
void inner_class_eq_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  shared_inner_class H = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(whether(G.get()==H.get())); // test identical objects
}
void inner_class_neq_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  shared_inner_class H = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(whether(G.get()!=H.get())); // test identical objects
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
    push_value(G->datum);
}
void dual_datum_of_inner_class_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(G->dual_datum);
}
@)
void inner_class_to_root_datum_coercion()
{@; root_datum_of_inner_class_wrapper(expression_base::single_value); }
@)
void dual_inner_class_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(G->dual());
}

@ More interestingly, let us extract the list of names of the real forms.
This uses the interface fields stored in the value. Since they exist for both
the group itself and for the dual group, we define an auxiliary function that
produces the list, and then use it twice.

@< Local function def...@>=
void push_name_list(const output::FormNumberMap& interface)
{ own_row result = std::make_shared<row_value>(0);
  result->val.reserve(interface.numRealForms());
  for (unsigned int i=0; i<interface.numRealForms(); ++i)
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
command from \.{mainmode.cpp}; that commands uses the function
|innerclass_io::printBlockSizes|, but we have to rewrite it to avoid calling an
output routine. A subtle difference is that we use a matrix of integers rather
than of |arithmetic::big_int| values to collect the block sizes; this avoids
having to define a new primitive type, and probably suffices for the cases where
actual computation of the block feasible. However, we also provide a function
the computes a single block size as an unbounded integer.

@< Local function def...@>=
void block_sizes_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l==expression_base::no_value)
    return;
@)
  own_matrix M = std::make_shared<matrix_value> @|
    (int_Matrix(G->val.numRealForms(),G->val.numDualRealForms()));
  for (unsigned int i = 0; i < M->val.n_rows(); ++i)
    for (unsigned int j = 0; j < M->val.n_columns(); ++j)
      M->val(i,j) =
      G->val.block_size(G->interface.in(i),G->dual_interface.in(j)).int_val();
  push_value(std::move(M));
}

void block_size_wrapper(expression_base::level l)
{ auto j = get<int_value>()->val.ulong_val();
  auto i = get<int_value>()->val.ulong_val();
  shared_inner_class G = get<inner_class_value>();
  if (i>=G->val.numRealForms())
  { std::ostringstream o;
    o << "Real form number " << i << " out of bounds";
    throw runtime_error(o.str());
  }
  if (j>=G->val.numDualRealForms())
  { std::ostringstream o;
    o << "Dual real form number " << j << " out of bounds";
    throw runtime_error(o.str());
  }
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
  unsigned int nr=G->val.numRealForms();
  unsigned int nc=G->val.numCartanClasses();
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(nr,nc));
  for (unsigned int i=0; i<nr; ++i)
  { BitMap b=G->val.Cartan_set(G->interface.in(i));
    for (unsigned int j=0; j<nc; ++j)
      M->val(i,j)= b.isMember(j) ? 1 : 0;
  }
  push_value(std::move(M));
}

@ We do the same for dual real forms. Note that we had to introduce the method
|dualCartanSet| for |InnerClass| in order to be able to write this function.

@< Local function def...@>=
void dual_occurrence_matrix_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l==expression_base::no_value)
    return;
@)
  unsigned int nr=G->val.numDualRealForms();
  unsigned int nc=G->val.numCartanClasses();
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(nr,nc));
  for (unsigned int i=0; i<nr; ++i)
  { BitMap b=G->val.dual_Cartan_set(G->dual_interface.in(i));
    for (unsigned int j=0; j<nc; ++j)
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
                ,"(RootDatum,mat->WeylElt,InnerClass)");
install_function(inner_class_eq_wrapper,@|"=","(InnerClass,InnerClass->bool)");
install_function(inner_class_neq_wrapper,@|"!=","(InnerClass,InnerClass->bool)");
install_function(distinguished_involution_wrapper,@|"distinguished_involution"
                ,"(InnerClass->mat)");
install_function(root_datum_of_inner_class_wrapper,@|"root_datum"
                ,"(InnerClass->RootDatum)");
install_function(dual_datum_of_inner_class_wrapper,@|"dual_datum"
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
%
A |RealReductiveGroup| object is dependent upon an |InnerClass|, to which it
stores a reference. In \.{atlas} a |real_form_value| is an independent value, so
we must make sure that |InnerClass| pointed to by its |RealReductiveGroup|
cannot disappear before it does. The easiest way to do this is to place an
|shared_inner_class| pointer |ic_ptr| inside the |real_form_value|.
Since that data can be accessed from inside the |RealReductiveGroup|, we shall
mostly mention the |ic_ptr| with the purpose of accessing its
|interface| and |dual_interface| fields. They however also serve for testing
(by pointer equality) whether two real forms are in the same inner class.

This class also serves to store persistent data related to the real form, in
values of type |KhatContext| and |Rep_table|. In order to avoid overhead at
construction, we lazily construct these tables, storing only pointer, which will
only be assigned on first use. These fields are marked mutable, so that we do
not have to do copy-on-write just to add the table (which would be silly, as
others might benefit from the addition; anyway we forbid copying of any
|real_form_value|) or use |non_const_get| to dissimulate our ``destructive''
intentions. The |val| field itself is |mutable| for similar reasons, notably the
KGB table gets initialised lazily inside the |RealReductiveGroup|. Another
advantage of storing pointers to |KhatContext| and |Rep_table| rather than such
objects, is that it avoids making our header file depend on the header files for
those objects.

@< Type definitions @>=
class real_form_value : public value_base
{ struct token@+{}; // serves to control access to the constructor
public:
  shared_inner_class ic_ptr;
  mutable RealReductiveGroup val;
@)
  real_form_value (shared_inner_class icp,RealFormNbr f,token)
@/: ic_ptr(icp), val(icp->val,f)
  , rt_p(nullptr) @+{}
@)
  real_form_value
    (shared_inner_class icp,RealFormNbr f
    ,const RatCoweight& coch, TorusPart tp) @/
  : ic_ptr(icp), val(icp->val,f,coch,tp)
  , rt_p(nullptr) @+{}
  virtual ~real_form_value ();
@)
  static shared_real_form build(const shared_inner_class& icp,RealFormNbr f);
  static shared_real_form build @|
   (const shared_inner_class& icp,RealFormNbr f,
    const RatCoweight& coch, const TorusPart& tp);
  virtual void print(std::ostream& out) const;
  static const char* name() @+{@; return "real form"; }
  real_form_value (const real_form_value& ) = delete;
@)
  const KGB& kgb () const @+{@; return val.kgb(); }
   // generate and return $K\backslash G/B$ set
  const Rep_context& rc() const;
  Rep_table& rt() const;
private:
  mutable Rep_table* rt_p;
    // owned pointers, initially |nullptr|, assigned at most once
};

@ Like |inner_class_value::build|, the static method |real_form_value::build|
first tests if it can get its result from converting a stored weak pointer to a
shared pointer, and it that fails constructs a fresh |real_form_value|,
returning a shared pointer to it while leaving a weak pointer version of it in
the same location where it had been looked up. In this case looking up is just
indexing the vector |icp->real_form_wptr| with the real form number~|f|.

@< Function def...@>=
shared_real_form real_form_value::build
  (const shared_inner_class& icp,RealFormNbr f)
{
  auto& w_ptr = icp->real_form_wptr[f];
  if (auto p = w_ptr.lock())
    return p; // reuse existing real form, if found in inner class
  auto result =
    std::make_shared<real_form_value>(icp,f,token()); // build new real form
  w_ptr = result; // store a weak pointer version inside |srd| table
  return result; // return pointer to modified instance of |inner_class_value|
}

@ We have a similar method when constructing a real form with an explicitly
given cocharacter |coch| and initial torus part (for the $K\backslash G/B$
construction) |tp|, but in this case we do not guarantee sharing of identically
generated values (which would require fairly extensive tables in the inner
class), except if the special values provided happen to be the same as what
would automatically be chosen if none wore specified, in which case we just drop
them and call the previous |build|~method. The default values can be found in
the relevant |RealReductiveGroup| constructor, and are computed by the function
|some_coch| (defined in the \.{innerclass} compilation unit) and by the
|InnerClass::x0_torus_part| method.

@< Function def...@>=
shared_real_form real_form_value::build @|
   (const shared_inner_class& icp, RealFormNbr f,
    const RatCoweight& coch, const TorusPart& tp)
{
  auto default_coch = some_coch(icp->val,icp->val.xi_square(f));
  if (coch==default_coch and tp==icp->val.x0_torus_part(f))
    return build(icp,f);
  return(std::make_shared<real_form_value>(icp,f,coch,tp));
}


@ The methods |rc| and |rt| ensure a |Rep_table| value is constructed at
|*rt_p|, and then these methods return an appropriate reference. The value so
obtained will serve to manipulate parameters for standard modules, for which
we shall define an Atlas type below. Storing the value here ensures that it
will be shared between different parameters for the same real form, and that
it will live as long as those parameter values do.

@< Function def...@>=
  const Rep_context& real_form_value::rc() const
    {@; return *(rt_p==nullptr ? rt_p=new Rep_table(val) : rt_p); }
  Rep_table& real_form_value::rt() const
    {@; return *(rt_p==nullptr ? rt_p=new Rep_table(val) : rt_p); }
@)
  real_form_value::~real_form_value () @+{@; delete rt_p; }

@ When printing a real form, we give the name by which it is known in the parent
inner class, and provide some information about its topology. The names of the
real forms are indexed by their outer number, but the real form itself stores
its inner number, so we must somewhat laboriously make the conversion here.

@< Function def...@>=
void real_form_value::print(std::ostream& out) const
{ if (val.isCompact()) out << "compact ";
  out << (val.isConnected() ? "connected " : "disconnected " );
  if (val.isQuasisplit())
    out << (val.isSplit() ? "" : "quasi") << "split ";
  out << "real group with Lie algebra '" @|
      << ic_ptr->interface.type_name
          (ic_ptr->interface.out(val.realForm())) @|
      << '\'' ;
}

@ To make a real form is easy: one provides an |inner_class_value| and a valid
index into its list of real forms. Since this number coming from the outside
is to be interpreted as an outer index, we must convert it to an inner index
at this point. We also inversely allow finding the index of a real form in its
inner class; of course the opposite conversion applies here.
As a special case we also provide the quasisplit form.

@< Local function def...@>=
void real_form_wrapper(expression_base::level l)
{ int i = get<int_value>()->int_val();
  shared_inner_class G = get<inner_class_value>();
  if (static_cast<unsigned int>(i)>=G->val.numRealForms())
  { std::ostringstream o;
    o << "Illegal real form number: " << i;
    throw runtime_error(o.str());
  }
@.Illegal real form number@>
  if (l!=expression_base::no_value)
    push_value(real_form_value::build(G,G->interface.in(i)));
}
@)
void form_number_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>@|
      (rf->ic_ptr->interface.out(rf->val.realForm())));
}
@)
void quasisplit_form_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l!=expression_base::no_value)
    push_value(real_form_value::build(G,G->val.quasisplit()));
}

@*2 Functions operating on real reductive groups.
%
From a real reductive group we can go back to its inner class or root datum;
the former is obtained by just copying the stored shared pointer |ic_ptr|.

@< Local function def...@>=
void inner_class_of_real_form_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(rf->ic_ptr);
}

void real_form_to_inner_class_coercion()
{@;push_value(get<real_form_value>()->ic_ptr); }

void real_form_to_root_datum_coercion()
{@;
  push_value(get<real_form_value>()->ic_ptr->datum);
}

@ Here is a function that gives information about the (dual) component group
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

@ Here are two somewhat technical function that will facilitate working ``in
coordinates'' with KGB elements for this real form. The first one, called
|base_grading_vector|, returns a rational coweight that determines the base
grading for the real form, and which is an offset that will be added to
|torus_bits| values when computing the |torus_factor| they represent. It is
defined so that a zero value corresponds to a quasisplit real form, which proves
the most useful base point. This does imply that a standard choice for the
``infinitesimal cocharacter''~$g$ for the real from must differ by
${}^\vee\!\rho$ from the value produced here, which is obtained directly from
the |RealReductiveGroup| class using the |g_rho_check| method (which name should
be read as $g-{}^\vee\!\rho$).

We used to ensure that the coweight returned as |base_grading_vector| is
dominant, while remaining in the coset of $2X_*$ defined by |G.g_rho_check()| so
as to not affect the torus element $\exp(\pi\ii t)$ giving the actual base
grading. There is however no obvious best way to do this, and the easy way that
used to be employed could give coweights unnecessarily far inside the dominant
chamber; in addition this shift may destroy the ${}^t\delta$-invariance of
|G.g_rho_check| that is hard to reestablish without risk of leaving the coset.
Since dominance is not a vital property, we therefore now just return
|G.g_rho_check()| unchanged.

The other function |initial_torus_bits| gives access to the value that is used
internally to initialise the full KGB construction, and which distinguishes this
real from other ones in the same square class (which share the same
|base_grading_vector|).

@< Local function def...@>=
void base_grading_vector_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<rational_vector_value>(rf->val.g_rho_check()));
}

void initial_torus_bits_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<vector_value> @|
      (int_Vector(rf->val.x0_torus_part())));
}


@ There is a partial ordering on the Cartan classes defined for a real form. A
matrix for this partial ordering is computed by the function |Cartan_order|,
which substitutes for the \.{Fokko} command \.{corder}.

@h "poset.h"
@< Local function def...@>=
void Cartan_order_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l==expression_base::no_value)
    return;
@)
  unsigned int n=rf->val.numCartan();
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(n,n,0));
  const poset::Poset& p = rf->val.innerClass().Cartan_ordering();
  for (unsigned int i=0; i<n; ++i)
    for (unsigned int j=i; j<n; ++j)
      if (p.lesseq(i,j)) M->val(i,j)=1;

  push_value(std::move(M));
}

@ A similar function is |KGB_Hasse| which encodes the Bruhat order on the KGB
set as a matrix. This function takes advantage of the fact that
|real_form_value::val| is marked |mutable|, so that the non-|const| method
|RealReductiveGroup::Bruhat_KGB| can be called even though we obtained a
|shared_real_form| (pointer to |const|) from the stack. Originally that field
was not marked |mutable|, and we had to revert |const|-casting in the code
below, via the use of |non_const_get| when pulling the argument from the stack.

@h "bruhat.h"

@< Local function def...@>=
void KGB_Hasse_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l==expression_base::no_value)
    return;
@)
  unsigned int n=rf->val.KGB_size();
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(n,n,0));
  const auto& Bruhat = rf->val.Bruhat_KGB();
  for (unsigned int j=0; j<n; ++j)
  { const auto& col = Bruhat.hasse(j);
    for (auto it=col.begin(); it!=col.end(); ++it)
      M->val(*it,j)=1;
  }
  push_value(std::move(M));
}


@ Finally we make available the equality test for real forms defined in the
library. However, since the |real_form_value::build| will identify real forms
obtained by selecting at equal |RealFormNbr| indices from the same
|inner_class_value|, it is useful to start with a pointer comparison for
efficiency.

@< Local function def...@>=

void real_form_eq_wrapper(expression_base::level l)
{ shared_real_form y = get<real_form_value>();
  shared_real_form x = get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(whether(x.get()==y.get() or x->val==y->val));
}

void real_form_neq_wrapper(expression_base::level l)
{ shared_real_form y = get<real_form_value>();
  shared_real_form x = get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(whether(x.get()!=y.get() and x->val!=y->val));
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
index. We also provide the dual quasisplit form. In both cases, the dual
|inner_calss_value| of~|G| is computed on the fly as |G->dual()|, after which a
|real_form_value| is obtained through |real_form_value|, as in
|real_form_wrapper| above.

@< Local function def...@>=
void dual_real_form_wrapper(expression_base::level l)
{ int i =get<int_value>()->int_val();
  shared_inner_class G = get<inner_class_value>();
  if (static_cast<unsigned int>(i)>=G->val.numDualRealForms())
  { std::ostringstream o;
    o << "Illegal dual real form number: " << i;
    throw runtime_error(o.str());
  }
@.Illegal dual real form number@>
  if (l==expression_base::no_value)
    return;
@)
  push_value(real_form_value::build(G->dual(),G->dual_interface.in(i)));
}
@)
void dual_quasisplit_form_wrapper(expression_base::level l)
{ shared_inner_class G = get<inner_class_value>();
  if (l==expression_base::no_value)
    return;
@)
  auto dual_ic=G->dual();
  push_value(real_form_value::build(dual_ic,dual_ic->val.quasisplit()));
}

@*2 Synthetic real forms.
%
It is useful to be able to compute a real form based on other information than
its index within its inner class, namely on a strong involution representative
(involution and torus element). The synthetic \.{atlas} function |real_form|
takes an inner class, a matrix giving an involution~$\theta$, and a rational
coweight~$t$ that will be projected onto the subspace fixed under (the right
action of)~$\theta$. If $\theta$ is in the inner class and the projected~$t$
describes (through $\exp_{-1}:l\mapsto\exp(\pi\ii l)$) a torus element whose
square is central, the function returns the corresponding real form whose
|RealReductiveGroup::square_class_character| member (which is the value returned
by the |g_check_rho| method) will have been computed from~$t$, inside the
|real_form_of| function (defined in the \.{innerclass} compilation unit). This
value only depends on the square of the torus element, but may differ from the
value stored for this (weak) real form selected by the number
|rf->val.realForm()| in the inner class. This difference notably allows the same
strong involution representative to be subsequently used to specify a KGB
element for this (strong) real form.

@:synthetic_real_form@>
@f theta nullptr

@< Local function def...@>=
TwistedInvolution twisted_from_involution
  (const InnerClass& G, WeightInvolution theta)
{ const RootDatum& rd = G.rootDatum();
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
  if (torus_factor->val.size()!=G->val.rank())
    throw runtime_error ("Torus factor size mismatch");
  TwistedInvolution tw = twisted_from_involution(G->val,theta->val);
  @< Test that |torus_factor| describes torus element whose square is central,
     and project to make it stable under the action of |theta| @>

  @< Make sure the involution table knows about the Cartan class of |tw|
     below in the Cartan ordering @>
  RatCoweight coch(0); // dummy value to be replaced
  RealFormNbr rf = real_form_of(G->val,tw,torus_factor->val,coch);
   // sets |coch|
  TorusPart tp = realredgp::minimal_torus_part @|
   (G->val,rf,coch,tw,torus_factor->val);
   // choose base |TorusPart| for $K\backslash G/B$ construction
  if (l!=expression_base::no_value)
    push_value(real_form_value::build(G,rf,coch,tp));
}

@ The square of the torus element for the projected |torus_factor| is central
if all simple roots have integral evaluation on |torus_factor|. In the code
below we compute the doubled projection first, and use |is_central| to test that
the evaluations are even, before halving to obtain the actual projection.

@< Test that |torus_factor| describes torus element whose square is central,
   and project to make it stable under the action of |theta| @>=
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
  torus_factor->val /= 2;
  // now $(1+\theta)/2$ has been applied to |torus_factor|
}

@ The call to |minimal_torus_part| uses the involution table in order to be
able to do downward (inverse) Cayley transforms, so we must ensure that any
involutions that can be encountered have been entered into the table.

@< Make sure the involution table knows about the Cartan class of |tw|... @>=
{ CartanNbr cn = G->val.class_number(tw);
  G->val.generate_Cartan_orbit(cn);
  const BitMap& b = G->val.Cartan_ordering().below(cn);
  for (auto it=b.begin(); it(); ++it)
    G->val.generate_Cartan_orbit(*it);
}

@ The method |central_fiber| of |InnerClass| can be accessed using following
function. The method computes those torus parts in the fiber at the
distinguished involution that both remain in the strong real form orbit and are
central (do not affect any gradings).

@< Local function def...@>=
void central_fiber_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  if (l==expression_base::no_value)
    return;
@)
  auto cf = rf->val.innerClass().central_fiber(rf->val.realForm());
  own_row result = std::make_shared<row_value>(cf.size());
  unsigned int i=0;
  for (auto it=cf.begin(); it!=cf.end(); ++it, ++i)
    result->val[i]= std::make_shared<vector_value>(int_Vector(*it));
  push_value(std::move(result));
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
install_function(count_Cartans_wrapper
                ,@|"nr_of_Cartan_classes","(RealForm->int)");
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
{ shared_inner_class ic_ptr;
  CartanNbr number;
  const CartanClass& val;
@)
  Cartan_class_value(shared_inner_class icp,CartanNbr cn);
  ~Cartan_class_value() @+{} // everything is handled by destructor of |parent|
@)
  virtual void print(std::ostream& out) const;
  static const char* name() @+{@; return "Cartan class"; }
  Cartan_class_value (const Cartan_class_value& ) = delete;
};
@)
typedef std::shared_ptr<const Cartan_class_value> shared_Cartan_class;

@ In the constructor we used to check that the Cartan class with the given
number currently exists, but now the |InnerClass::cartan| method
assures that one is generated if this should not have been done before. We
therefore call that method in the initialiser; on return it provides a valid
reference.

@< Function def...@>=
Cartan_class_value::Cartan_class_value(shared_inner_class icp,CartanNbr cn)
: ic_ptr(icp), number(cn),val(icp->val.cartan(cn))
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
  if (static_cast<unsigned int>(i)>=ic->val.numCartanClasses())
  { std::ostringstream o;
    o << "Illegal Cartan class number: " << i @|
    << ", this inner class only has " << ic->val.numCartanClasses()
    << " of them";
    throw runtime_error(o.str());
  }
@.Illegal Cartan class number@>
  if (l!=expression_base::no_value)
    push_value(std::make_shared<Cartan_class_value>(ic,i));
}

@ Alternatively (and this used to be the only way) one can provide a
|real_form_value| together with a valid index into \emph{its} list of Cartan
classes. We translate this number into an index into the list for its
containing inner class, and then get the Cartan class from there.

@< Local function def...@>=
void rf_Cartan_class_wrapper(expression_base::level l)
{ int i=get<int_value>()->int_val();
  shared_real_form rf= get<real_form_value>();
  if (static_cast<unsigned int>(i)>=rf->val.numCartan())
  { std::ostringstream o;
    o << "Illegal Cartan class number: " << i @|
      << ", this real form only has " << rf->val.numCartan()
      << " of them";
    throw runtime_error(o.str());
  }
@.Illegal Cartan class number@>
  BitMap cs=rf->val.Cartan_set();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<Cartan_class_value> @|
      (rf->ic_ptr,cs.n_th(i)));
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
    push_value(std::make_shared<Cartan_class_value> @|
      (rf->ic_ptr,rf->val.mostSplit()));
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
    cc->ic_ptr->val.involution_of_Cartan(cc->number);
  WeylWord ww = cc->ic_ptr->val.weylGroup().word(tw);

  push_value(std::make_shared<vector_value>
    (std::vector<int>(ww.begin(),ww.end())));

  push_value(std::make_shared<int_value>(cc->val.orbitSize()));
  push_value(std::make_shared<int_value>(cc->val.fiber().fiberSize()));
  wrap_tuple<2>();

  const RootSystem& rs=cc->ic_ptr->val.rootDatum();

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
  shared_inner_class icp=cc->ic_ptr;
  own_row result = std::make_shared<row_value>(cc->val.numRealForms());
  for (unsigned int i=0,k=0; i<icp->val.numRealForms(); ++i)
  { RealFormNbr rf = icp->interface.in(i);
    BitMap b(icp->val.Cartan_set(rf));
    if (b.isMember(cc->number))
      result->val[k++] = real_form_value::build(icp,rf);
  }
  push_value(std::move(result));
}
@)
void dual_real_forms_of_Cartan_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  if (l==expression_base::no_value)
    return;
@)
  auto dual_ic = cc->ic_ptr->dual();
  own_row result = std::make_shared<row_value>(cc->val.numDualRealForms());
  for (unsigned int i=0,k=0; i<dual_ic->val.numRealForms(); ++i)
  { RealFormNbr drf = cc->ic_ptr->dual_interface.in(i);
    BitMap b (cc->ic_ptr->val.dual_Cartan_set(drf));
    if (b.isMember(cc->number))
      result->val[k++] = real_form_value::build(dual_ic,drf);
  }
  push_value(std::move(result));
}

@ For the fiber group partition information that was not produced by
|print_Cartan_info|, we use a Cartan class and a real form as parameters. This
function returns the part of the |weakReal| partition stored in the fiber of the
Cartan class that corresponds to the given (weak) real form. The numbering of
the parts of that partition are not the numbering of the real forms themselves,
so they must be translated through the |realFormLabels| list for the Cartan
class, which must be obtained from its inner class. If the Cartan class does not
exist for the given real form, then it will not occur in that |realFormLabels|
list, and the part returned here will be empty. The part of the partition is
returned as a list of integral values.

@< Local function def...@>=
void fiber_partition_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
  shared_Cartan_class cc(get<Cartan_class_value>());
  if (rf->ic_ptr.get()!=cc->ic_ptr.get())
    throw runtime_error
    ("Inner class mismatch between real form and Cartan class");
@.Inner class mismatch...@>
  BitMap b(cc->ic_ptr->val.Cartan_set(rf->val.realForm()));
  if (!b.isMember(cc->number))
    throw runtime_error("Cartan class not defined for this real form");
@.Cartan class not defined@>
  if (l==expression_base::no_value)
    return;
@)
  const Partition& pi = cc->val.fiber().weakReal();
  const RealFormNbrList rf_nr=
     cc->ic_ptr->val.realFormLabels(cc->number);
     // translate part number of |pi| to real form
  own_row result =
    std::make_shared<row_value>(0); // cannot predict exact size here
  for (unsigned int i=0; i<pi.size(); ++i)
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
  const output::FormNumberMap rfi = cc->ic_ptr->interface;
  const RealFormNbrList& rfl = cc->ic_ptr->val.realFormLabels(cc->number);
  if (l==expression_base::no_value)
    return;
@)
  unsigned int n_sq_classes = cc->val.numRealFormClasses();
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
  if (rf->ic_ptr.get()!=cc->ic_ptr.get())
    throw runtime_error
    ("Inner class mismatch between real form and Cartan class");
@.Inner class mismatch...@>
  BitMap b(cc->ic_ptr->val.Cartan_set(rf->val.realForm()));
  if (!b.isMember(cc->number))
    throw runtime_error ("Cartan class not defined for this real form");
@.Cartan class not defined...@>
@)
  const Partition& pi = cc->val.fiber().weakReal();
  const RealFormNbrList rf_nr=
     cc->ic_ptr->val.realFormLabels(cc->number);
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
{ cm=cc->ic_ptr->val.rootDatum().Cartan_matrix(si);
  dynkin::DynkinDiagram d(cm); sigma = d.perm();
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
    for (unsigned int i=0; i<si.size(); ++i)
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
  for (unsigned int i=0; i<pi.size(); ++i)
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
{ shared_real_form rf;
  KGBElt val;
@)
  KGB_elt_value(const shared_real_form& form, KGBElt x) : rf(form), val(x) @+{}
@)
  virtual void print(std::ostream& out) const;
  static const char* name() @+{@; return "KGB element"; }
  KGB_elt_value (const KGB_elt_value& ) = default;
    // we use |get_own<KGB_elt_value>|
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
  shared_real_form rf= get<real_form_value>();
  if (static_cast<unsigned int>(i)>=rf->val.KGB_size())
  { std::ostringstream o;
    o << "Inexistent KGB element: " << i;
    throw runtime_error(o.str());
  }
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
    push_value(std::make_shared<Cartan_class_value> @|
      (x->rf->ic_ptr,kgb.Cartan_class(x->val)));
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
  { std::ostringstream o;
    o << "Illegal root index: " << root_index;
    throw runtime_error(o.str());
  }
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
  push_value(std::move(x));
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
  push_value(std::move(x));
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
  shared_real_form rf = get<real_form_value>();

  if (torus_factor->val.size()!=rf->val.rank())
    throw runtime_error("Torus factor size mismatch");
@)
  Ratvec_Numer_t& num = torus_factor->val.numerator();
  { // make theta-fixed and remove base grading vector offset
    num += theta->val.right_prod(num);
    ((torus_factor->val /= 2) -=rf->val.g_rho_check()).normalize() ;
    if (torus_factor->val.denominator()!=1)
      throw runtime_error("Torus factor not in cocharacter coset of real form");
@.Torus factor not in cocharacter...@>
  }

  const InnerClass& G = rf->val.innerClass();
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
  x->val= kgb.twisted(x->val,x->rf->val.innerClass().distinguished());
    // do twist
  push_value(std::move(x));
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
  test_compatible(x->rf->val.innerClass(),delta);
  if (l==expression_base::no_value)
    return;
@)
  x->val= x->rf->kgb().twisted(x->val,delta->val); // do twist
  push_value(std::move(x));
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

@ Like other data types we have seen, we include shared pointers to parent
objects to ensure these remain in existence as long as our block does; in fact
we include two such shared pointers, one for each real form. The |val| field
contains an actual |Block| instance, which is constructed when the |Block_value|
is. We also reserve a field |kl_tab| in the structure to store KL polynomials,
though they will only be computed once they are asked for.

The constructor for |Block_value| should avoid calling the version of the
|Block::build| method that takes real form and deal real form numbers (as we
originally did here), as this will generate small versions of the KGB sets,
which changes the numbering and causes inconsistencies when KGB elements are
extracted from block elements. Instead we call the |build| method that accepts a
pair of |RealReductiveGroup| arguments, so that the generated block will be
exactly a classical one.

@< Type definitions @>=
struct Block_value : public value_base
{ const shared_real_form rf; const shared_real_form dual_rf;
  mutable Block val; // Bruhat order may be generated implicitly
  mutable kl::KL_table kl_tab; // as may KLV polynomials
@)
  Block_value(const shared_real_form& form,
              const shared_real_form& dual_form)
  : rf(form), dual_rf(dual_form)
  , val(Block::build(rf->val,dual_rf->val))
  , kl_tab(val)
  {}
  ~Block_value() @+{}
@)
  virtual void print(std::ostream& out) const;
  static const char* name() @+{@; return "KGB element"; }
  Block_value (const Block_value& ) = delete;
};
@)
typedef std::unique_ptr<Block_value> Block_ptr;
typedef std::shared_ptr<const Block_value> shared_Block;
typedef std::shared_ptr<Block_value> own_Block;

@ When printing a block, we print its size; we shall later provide a separate
print function that tabulates its individual elements.

@< Function def...@>=
void Block_value::print(std::ostream& out) const
{@; out << "Block of " << val.size() << " elements"; }

@ To make a block, one provides a |real_form_value| and a |dual_real_form|.
This function is named in a tribute to the creator of the Atlas software, and
of the data structure that is constructed and stored here.

@< Local function def...@>=
bool is_dual(const shared_inner_class& ic0, const shared_inner_class& ic1)
{ return ic0->dual_datum.get()==ic1->datum.get() and
  ic0->val.dualDistinguished()==ic1->val.distinguished();
}
void Fokko_block_wrapper(expression_base::level l)
{ shared_real_form drf=get<real_form_value>();
  shared_real_form rf=get<real_form_value>();
@)
  if (not is_dual(rf->ic_ptr,drf->ic_ptr))
    throw runtime_error @|
    ("Inner class mismatch between real form and dual real form");
@.Inner class mismatch...@>
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
  { std::ostringstream o;
    o << "Block element " << z
      << " out of range (<" << b->val.size() << ")";
    throw runtime_error(o.str());
  }
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<KGB_elt_value>(b->rf,b->val.x(z)));
  auto dic = b->rf->ic_ptr->dual();
  auto drf = real_form_value::build(dic,b->dual_rf->val.realForm());
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
  if (b->rf->ic_ptr.get()!=x->rf->ic_ptr.get())
    throw runtime_error("Real form not in inner class of block");
  if (not is_dual(b->rf->ic_ptr,y->rf->ic_ptr))
    throw runtime_error("Dual real form not in inner class of block");
  const KGB& kgb = b->rf->kgb(); const KGB& dual_kgb = b->dual_rf->kgb();
  const TwistedWeylGroup& tw = kgb.twistedWeylGroup();
  const TwistedWeylGroup& dual_tw = dual_kgb.twistedWeylGroup();
  if (blocks::dual_involution(kgb.involution(x->val),tw,dual_tw)
      != dual_kgb.involution(y->val))
    throw runtime_error("Fiber mismatch KGB and dual KGB elements");

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
  if (s>=b->rf->val.semisimple_rank())
  { std::ostringstream o;
    o << "Illegal simple reflection: " << s;
    throw runtime_error(o.str());
  }
  if (i>=b->val.size())
  { std::ostringstream o;
    o << "Block element " << i @|
      << " out of range (<" << b->val.size() << ')';
    throw runtime_error(o.str());
  }
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
  if (s>=b->rf->val.semisimple_rank())
  { std::ostringstream o;
    o << "Illegal simple reflection: " << s;
    throw runtime_error(o.str());
  }
  if (i>=b->val.size())
  { std::ostringstream o;
    o << "Block element " << i
      << " out of range (<" << b->val.size() << ")";
    throw runtime_error(o.str());
  }
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(b->val.cross(s,i)));
}
@)
void block_Cayley_wrapper(expression_base::level l)
{ shared_int i = get<int_value>();
  unsigned int ii = i->int_val();
  shared_Block b = get<Block_value>();
  unsigned int s = get<int_value>()->int_val();
  if (s>=b->rf->val.semisimple_rank())
  { std::ostringstream o;
    o << "Illegal simple reflection: " << s;
    throw runtime_error(o.str());
  }
  if (ii >= b->val.size())
  { std::ostringstream o;
    o << "Block element " << ii
      << " out of range (<" << b->val.size() << ")";
    throw runtime_error(o.str());
  }
  if (l==expression_base::no_value)
    return;
@)
  BlockElt sx = b->val.cayley(s,ii).first;
  if (sx==UndefBlock) // when undefined, return i to indicate so
    push_value(std::move(i));
  else
    push_value(std::make_shared<int_value>(sx));
}
@)
void block_inverse_Cayley_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  unsigned int ii = i->int_val();
  shared_Block b = get<Block_value>();
  unsigned int s = get<int_value>()->int_val();
  if (s>=b->rf->val.semisimple_rank())
  { std::ostringstream o;
    o << "Illegal simple reflection: " << s;
    throw runtime_error(o.str());
  }
  if (ii >= b->val.size())
  { std::ostringstream o;
    o << "Block element "  << ii @|
      << " out of range (<" << b->val.size() << ")";
    throw runtime_error(o.str());
  }
  if (l==expression_base::no_value)
    return;
@)
  BlockElt sx = b->val.inverseCayley(s,ii).first;
  if (sx==UndefBlock) // when undefined, return i to indicate so
    push_value(std::move(i));
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

@*1 A class for split integers.
%
We introduce a type extending integers to elements of the group
algebra over $\Zee$ of a (cyclic) group of order~$2$, in other words
$\Zee$-linear combinations of $1$ and another unit, written~$s$ and satisfying
$s^2=1$. We call these numbers split integers, and in the \.{atlas} language
designate the corresponding basic type as \.{Split}. The type will serve notably
to keep track of signatures of Hermitian forms.

@< Includes needed in the header file @>=
#include "arithmetic.h"

@ Although the necessary operations could easily be defined in the \.{axis}
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
  static const char* name() @+{@; return "split integer"; }
  split_int_value (const split_int_value& v) = default;
    // we use |get_own<split_int_value>|
};
@)
typedef std::unique_ptr<split_int_value> split_int_ptr;
typedef std::shared_ptr<const split_int_value> shared_split_int;
typedef std::shared_ptr<split_int_value> own_split_int;

@ Like for parameter values, a printing function |print_split| on the level
of a bare |Split_integer| value is defined in the library, which can be used
also in situations where the method |split_int_value::print| cannot. The latter
method simply calls it.

@< Function def... @>=

void split_int_value::print(std::ostream& out) const @+
{@; print_split(out,val); }

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
  {@; i->val+=j; push_value(std::move(i)); }
}

void split_minus_wrapper(expression_base::level l)
{ Split_integer j=get<split_int_value>()->val;
  own_split_int i=get_own<split_int_value>();
  if (l!=expression_base::no_value)
  {@; i->val-=j; push_value(std::move(i)); }
}
@)
void split_unary_minus_wrapper(expression_base::level l)
{ own_split_int i=get_own<split_int_value>();
  if (l!=expression_base::no_value)
  {@; i->val.negate(); push_value(std::move(i)); }
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

@ Here we install the built-in functions for split integers.

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

@*1 $K$-types.
%
We implement a data type for representing either a standard or irreducible
representation (depending on the context, the latter interpretation being the
unique irreducible quotient of the former) of the complexified maximal compact
subgroup $K$ for a real reductive group. We call these values $K$-types, and in
the \.{atlas} language designate the corresponding basic type as \.{KType}. A
$K$-type is defined by a pair $(x,\lambda-\rho)$ where $x$ is a KGB element,
$\lambda-\rho$ is a weight in $X^*$ that indirectly represents a character
$\lambda$ of the $\rho$-cover of a real Cartan subgroup; the value of the latter
is a class modulo the sub-lattice $(1-\theta_x)X^*$ where $\theta_x$ is the
involution associated to $x$, and $\lambda$ is therefore always reduced to an
elected representative modulo that lattice. The internal type used to represent
$K$-types is called |K_repr::K_type|, defined in the file~\.{K\_repr.h}. Most
actions on values of this type, including their construction, are performed by
methods of the class |Rep_context| defined in the file~\.{repr.h}.

@<Includes needed in the header file @>=
#include "K_repr.h"
#include "repr.h"

@*2 Class definition.
Like for KGB elements, we maintain a shared pointer to the real form value, so
that it will be assured to survive as long as parameters for it exist, and we
can access notably the |Rep_context| that it provides.

@< Type definitions @>=
struct K_type_value : public value_base
{ shared_real_form rf;
  K_repr::K_type val;
@)
  K_type_value(const shared_real_form& form, K_repr::K_type&& v)
  : rf(form), val(std::move(v)) @+{}
  ~K_type_value() @+{}
@)
  virtual void print(std::ostream& out) const;
  static const char* name() @+{@; return "K-type"; }
  K_type_value (const K_type_value& v) : rf(v.rf), val(v.val.copy()) {}
    // we use |get_own<K_type_value>|

@)
  const Rep_context& rc() const @+{@; return rf->rc(); }
  Rep_table& rt() const @+{@; return rf->rt(); }
};
@)
typedef std::unique_ptr<K_type_value> K_type_ptr;
typedef std::shared_ptr<const K_type_value> shared_K_type;
typedef std::shared_ptr<K_type_value> own_K_type;

@ When printing a $K$-type, we shall indicate a pair $(x,\lambda)$ that defines
it. Since we shall need to print |K_repr::K_type| values in other contexts as
well, an auxiliary output function |repr::print_K_type| is defined
in \.{basic\_io.h}, which takes and additional |Rep_context| argument; we shall
call that function from |K_type_value::print|. By choosing a name for the
auxiliary function different from |print|, we avoid having that call being
mistaken for a recursive call.

The virtual method |K_type_value::print|, is used when printing a
value of type \.{KType} (as opposed to for instance printing a term of
a \.{KTypePol}). Here we prefix the $K$-type text proper with additional
information about the $K$-type that may be relevant to the \.{atlas} user.

@< Function definition... @>=
void K_type_value::print(std::ostream& out) const
{
  out << @< Expression for adjectives that apply to a $K$-type @>@;@;;
  print_K_type(out << " K-type",val,rc());
}

@ We call a root singular if the corresponding coroot vanishes on
$(1+\theta_x)\lambda$ (for real roots this is always the case, for imaginary
roots it is equivalent to vanishing on~$\lambda$). We then classify different
levels of ``niceness'' of $K$-type by one of the adjectives ``non-standard''
(when $(1+\theta_x)\lambda$ fails to be imaginary-dominant), ``non-dominant''
(when $(1+\theta_x)\lambda$ fails to be dominant), ``zero'' (there are singular
compact simply-imaginary roots), ``non-final'' (there are real parity roots),
``non-normal'' (there are complex singular descent roots, making the $K$-type
equivalent to one at a lower involution in the Cartan class), or finally
``final'' (the good ones that could go into a \.{KTypePol} value; the condition
|is_final| should apply, though it is not tested here).

@< Expression for adjectives that apply to a $K$-type @>=
( not rc().is_standard(val) ? "non-standard"
@|: not rc().is_dominant(val) ? "non-dominant"
@|: not rc().is_nonzero(val) ? "zero"
@|: not rc().is_semifinal(val) ? "non-final"
@|: not rc().is_normal(val) ? "non-normal"
@|: "final")

@ To make a $K$-type, one should provide a KGB element~$x$, and a weight
$\lambda-\rho\in X^*$.

@< Local function def...@>=
void K_type_wrapper(expression_base::level l)
{ shared_vector lam_rho(get<vector_value>());
  shared_KGB_elt x = get<KGB_elt_value>();
  if (lam_rho->val.size()!=x->rf->val.rank())
  { std::ostringstream o;
    o << "Rank mismatch: (" @|
        << x->rf->val.rank() << ',' << lam_rho->val.size() << ")";
    throw runtime_error(o.str());
  }
@.Rank mismatch@>
  if (l!=expression_base::no_value)
    push_value(std::make_shared<K_type_value> @|
      (x->rf, x->rf->rc().sr_K(x->val,lam_rho->val)));
}

@*2 Functions operating on $K$-types.
%
The following function, which we shall bind to the monadic operator `|%|',
transforms a $K$-type value into a pair of values $(x,\lambda-\rho)$
that defines it. Although $\lambda$ is determined only modulo
$(1-\theta_x)X^*$, a fixed representative is always chosen, and it is impossible
to influence this choice, or to have distinct $K$-types that differ only in the
choice of representative. That functionality is ensured by the implementation of
|K_repr::K_type|, so here we need not worry about it.

@< Local function def...@>=
void unwrap_K_type_wrapper(expression_base::level l)
{ shared_K_type p = get<K_type_value>();
  if (l==expression_base::no_value)
    return;
@)
  push_value(std::make_shared<KGB_elt_value>(p->rf,p->val.x()));
  push_value(std::make_shared<vector_value>(p->val.lambda_rho()));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@)

void real_form_of_K_type_wrapper(expression_base::level l)
{ shared_K_type p = get<K_type_value>();
  if (l!=expression_base::no_value)
    push_value(p->rf);
}

@ Here is one more useful function: computing the height of a $K$-type. This is
the same height when comparing to the |bound| argument in functions like
|K_type_formula| and |branch|. This height is precomputed and stored inside
|K_repr::K_type| values themselves, so we simply get it from there.

@< Local function def...@>=
void K_type_height_wrapper(expression_base::level l)
{ shared_K_type p = get<K_type_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(p->val.height()));
}

@ While the equality and inequality operators for $K$-types test for strict
equality of all components (including of the real forms, which the method
|K_repr::K_type::operator==| simply assumes to be equal, as it has no access to
the real form), a separate function |equivalent| that tests an equivalence
relation generated by action by complex simple reflections (on the $x$ component
by cross action and on the $\lambda$ component) to move between fibres over
involutions for the same Cartan class. For two $K$-types in a fiber over a same
involution, this equivalence amounts to equality. Unlike the equality tests, one
must ensure the arguments are associated to the same real form before calling
|equivalent|, since a runtime error (rather than returning false) will occur
when this is not the case.

@< Local function def...@>=
void K_type_eq_wrapper(expression_base::level l)
{ shared_K_type q = get<K_type_value>();
  shared_K_type p = get<K_type_value>();
  if (l!=expression_base::no_value)
    push_value(whether(p->rf->val==q->rf->val and p->val==q->val));
}
void K_type_neq_wrapper(expression_base::level l)
{ shared_K_type q = get<K_type_value>();
  shared_K_type p = get<K_type_value>();
  if (l!=expression_base::no_value)
    push_value(whether(p->rf->val!=q->rf->val or p->val!=q->val));
}
@)
void K_type_equivalent_wrapper(expression_base::level l)
{ auto q = get_own<K_type_value>();
  auto p = get_own<K_type_value>();
  if (p->rf->val!=q->rf->val)
    throw runtime_error @|
      ("Real form mismatch when testing equivalence");
  if (l!=expression_base::no_value)
    push_value(whether(p->rc().equivalent(std::move(p->val),std::move(q->val))));
}

@ The various predicates that are used above in printing a $K$-type are also
available as functions returning a Boolean value. The function |is_normal| is
excluded here because that method assumes it is only called when |is_nonzero|
and |is_semifinal| both hold, in which case a simple implementation works. We
can ensure that property whenever it is called internally, but not when calling
as a built-in function (without additional code to test for and implement the
more general case). Instead there is a built-in |is_final| that tests for the
absence of \emph{any} singular descent (which is equivalent to |is_normal| in
the case |is_nonzero| and |is_semifinal| both hold), and we define a built-in
function |normalize| below that produces a unique representative in any
equivalence class of $K$-types.

@< Local function def...@>=
void K_type_is_standard_wrapper(expression_base::level l)
{ shared_K_type p = get<K_type_value>();
  if (l!=expression_base::no_value)
    push_value(whether(p->rc().is_standard(p->val)));
}

void K_type_is_dominant_wrapper(expression_base::level l)
{ shared_K_type p = get<K_type_value>();
  if (l!=expression_base::no_value)
    push_value(whether(p->rc().is_dominant(p->val)));
}

void K_type_is_zero_wrapper(expression_base::level l)
{ shared_K_type p = get<K_type_value>();
  if (l!=expression_base::no_value)
    push_value(whether(not p->rc().is_nonzero(p->val)));
}

void K_type_is_semifinal_wrapper(expression_base::level l)
{ shared_K_type p = get<K_type_value>();
  if (l!=expression_base::no_value)
    push_value(whether(p->rc().is_semifinal(p->val)));
}

void K_type_is_final_wrapper(expression_base::level l)
{ shared_K_type p = get<K_type_value>();
  if (l!=expression_base::no_value)
    push_value(whether(p->rc().is_final(p->val)));
}

@ Here are functions to transform a $K$-type to an equivalent one that is
dominant, which is only possible is the parameter is already standard, and to
one that is the elected representative of its equivalence class. The test for
being standard is not present in |K_type_dominant_wrapper|, because the
|make_dominant| method signals failure of being standard on the fly by throwing
an exception. The other three functions below will accept any $K$-type;
|to_canonical_fiber| moves it to the elected fiber for its Cartan class, while
|normal| produces an elected representative of the equivalence class of the
$K$-type, and |theta_stable| ensures the absence of any complex descents. The
$K$-type returned by |normal| has the additional property of satisfying
|is_final| whenever such a $K$-type exists in the equivalence class, and in any
case being complex-dominant and not having any singular (for
$(1+\theta_x)\lambda$) complex descents. Both |theta_stable| and |normal|
greedily perform complex reflections (descents in the former case, ones towards
complex dominance or singular descents in the latter case) in a chosen order,
which leads to a choice of one among possibly multiple candidate results. In
case of |normal| dependence on the initial class representative is removed by
applying |to_canonical_fiber|; also if a final representative exists, it is
unique.

@< Local function def...@>=
void K_type_dominant_wrapper(expression_base::level l)
{ own_K_type p = get_own<K_type_value>();

  if (l!=expression_base::no_value)
  {@; p->rc().make_dominant(p->val);
    push_value(std::move(p));
  }
}
@)
void to_canonical_fiber_wrapper(expression_base::level l)
{ own_K_type p = get_own<K_type_value>();

  if (l!=expression_base::no_value)
  @/{@; p->rc().to_canonical_involution(p->val);
    push_value(std::move(p));
  }
}
@)
void K_type_normal_wrapper(expression_base::level l)
{ own_K_type p = get_own<K_type_value>();
  if (l!=expression_base::no_value)
  {@; p->rc().normalise(p->val);
    push_value(std::move(p));
  }
}
@)
void K_type_theta_stable_wrapper(expression_base::level l)
{ own_K_type p = get_own<K_type_value>();
  if (l!=expression_base::no_value)
  @/{@; p->rc().make_theta_stable(p->val);
    push_value(std::move(p));
  }
}



@*1 Polynomials formed from $K$-types.
%
When working with $K$-types, and notably when produced by the deformation
formulas, the need arises to keep track of formal sums of them with split
integer coefficients. We call these formal sums $K$-type polynomials, and in
the \.{atlas} language designate the corresponding basic type as \.{KTypePol}.

@*2 Class definition.
%
The library provides a type |K_repr::K_type_pol| in which such sums can be
efficiently maintained. In order to use it we must have seen the header file
for the module \.{free\_abelian} on which the implementation is based. While
that class itself does not have such an invariant, the handling of these
formal sums in \.{atlas} will be such that all terms are ensured to have the
predicate |is_final| true, which ensures that all contributions are rewritten
into a unique form under the set of relations given by representation theory.
This will in particular ensure that equivalent terms are always be combined, and
the test for the sum being zero therefore mathematically correct.

@< Includes needed in the header file @>=
#include "free_abelian.h" // needed to make |SR_poly| a complete type

@~Like for KGB elements, we maintain a shared pointer to the real form value, so
that it will be assured to survive as long as $K$-types for it exist. Since
|K_repr::K_type_pol| is a move-only type (it has a move constructor but no copy
constructor), we only allow move construction of a |K_type_pol_value| from a
|K_repr::K_type_pol|. In doing so we have the occasion to simplify the stored
value by calling the |flatten| method, which takes its object by
rvalue-reference and returns the object after modification. Thus many
polynomials produced by built-in functions will start off without zero terms,
but the user may perform operations like adding single terms which can be
realised at less cost if we occasionally allow zero terms to be created in the
internal representation; therefore neither flattened form nor absence of zero
terms will not be a class invariant for |K_type_pol_value|.

@< Type definitions @>=
struct K_type_pol_value : public value_base
{ shared_real_form rf;
  K_repr::K_type_pol val;
@)
  K_type_pol_value(const shared_real_form& form, K_repr::K_type_pol&& v)
  : rf(form), val(std::move(v).flatten()) @+{}
  ~K_type_pol_value() @+{}
@)
  virtual void print(std::ostream& out) const;
  static const char* name() @+{@; return "module $K$-type"; }
  K_type_pol_value (const K_type_pol_value& v);
    // we use |uniquify<K_type_pol_value>|
@)
  const Rep_context& rc() const @+{@; return rf->rc(); }
  Rep_table& rt() const @+{@; return rf->rt(); }
  void assign_coef(const K_type_value& t, const Split_integer& c);
};
@)
typedef std::unique_ptr<K_type_pol_value> K_type_pol_ptr;
typedef std::shared_ptr<const K_type_pol_value> shared_K_type_pol;
typedef std::shared_ptr<K_type_pol_value> own_K_type_pol;

@ Printing a virtual module value calls the free function
|K_repr::print_K_type_pol| to do the actual work. It traverses the |std::map|
that is hidden in the |Free_Abelian| class template, and prints individual terms
by printing the |Split_integer| coefficient, followed by the $K$-type through a
call of |print_stdrep|. When either all coefficients are integers or all
coefficients are (integer) multiples of~$s$, it suppresses the component that is
always~$0$; this is particularly useful if polynomials are used to encode
$\Zee$-linear combinations of $K$-types.

@< Function def...@>=
void K_type_pol_value::print(std::ostream& out) const
{@; print_K_type_pol(out,val,rc()); }

@ For once we need a non-defaulted copy constructor, because |K_repr::K_type|
has no copy constructor, providing instead a method |copy| that must be
explicitly called (this was designed in order to better control the places where
actual copies must be made). Since this is a component type buried inside the
container |K_repr::K_type_pol| we copy the individual terms from the |val| field
first into a vector |accumulator|, and then move-construct a new
|K_repr::K_type_pol| (which has a constructor for that purpose) from the vector.

@< Function def...@>=
K_type_pol_value::K_type_pol_value(const K_type_pol_value& v)
: rf(v.rf), val()
{
  K_repr::K_type_pol::poly accumulator;
  accumulator.reserve(v.val.size());
  for (const auto& term : v.val)
    accumulator.emplace_back(term.first.copy(),term.second);
  val = K_repr::K_type_pol(std::move(accumulator),false,v.val.cmp());
}

@*2 Functions for $K$-type polynomials.
%
To start off a |K_type_pol_value|, one usually takes an empty sum, but one needs
to specify a real form to fill the |rf| field. The information allows us to
extract the real form from a $K$-type polynomial even if it is empty. We allow
testing the number of terms of the sum, and directly testing the sum to be
empty.

Testing two $K$-type polynomials for equality is also implemented. This could be
done by subtracting and then testing the result for being zero (empty), but it
is more efficient to just traverse both in parallel and stop once a difference
is found.

@< Local function def...@>=
void K_type_pol_wrapper(expression_base::level l)
{ shared_real_form rf = get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<K_type_pol_value> @| (rf,K_repr::K_type_pol()));
}
@)
void real_form_of_K_type_pol_wrapper(expression_base::level l)
{ shared_K_type_pol m = get<K_type_pol_value>();
  if (l!=expression_base::no_value)
    push_value(m->rf);
}
@)
void K_type_pol_unary_eq_wrapper(expression_base::level l)
{ shared_K_type_pol m = get<K_type_pol_value>();
  if (l!=expression_base::no_value)
    push_value(whether(m->val.is_zero()));
}

void K_type_pol_unary_neq_wrapper(expression_base::level l)
{ shared_K_type_pol m = get<K_type_pol_value>();
  if (l!=expression_base::no_value)
    push_value(whether(not m->val.is_zero()));
}

@)
void K_type_pol_eq_wrapper(expression_base::level l)
{ shared_K_type_pol n = get<K_type_pol_value>();
  shared_K_type_pol m = get<K_type_pol_value>();
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
void K_type_pol_neq_wrapper(expression_base::level l)
{ shared_K_type_pol n = get<K_type_pol_value>();
  shared_K_type_pol m = get<K_type_pol_value>();
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

@ This function must be global, it is declared in the header file \.{global.h}.

@< Function definitions @>=
void K_type_pol_size_wrapper(expression_base::level l)
{ shared_K_type_pol m = get<K_type_pol_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(m->val.count_terms()));
}


@ We allow implicitly converting a $K$-type to a $K$-type polynomial, which
involve expansion into \emph{final} $K$-types (there can be zero, one, or more
of them, and they can have positive or negative integer coefficients), to
initiate the invariant that only (standard) final $K$-types (with dominant
$(1+\theta_x)\lambda$) can be stored in a |K_type_pol_value|.

@< Local function def...@>=
void K_type_to_poly()
{ shared_K_type p = get<K_type_value>();
  const auto& rf=p->rf;
  auto final_K_types = rf->rc().finals_for(p->val.copy());
  K_repr::K_type_pol result;
  for (auto it=final_K_types.begin(); not final_K_types.at_end(it); ++it)
    result.add_term(std::move(it->first),Split_integer(it->second));
  push_value(std::make_shared<K_type_pol_value>(p->rf,std::move(result)));
}

@ There also is function to extract the coefficient (multiplicity) of a given
$K$-type in a $K$-type polynomial. However, it is bound to the array subscription
syntax, and therefore does not have a wrapper function. Instead, it is
implemented the \.{axis} module, as the |evaluate| method of the
|K_type_pol_coefficient| class derived from |subscr_base|.

In a subscription of a $K$-type polynomial by a $K$-type, the arguments are not
initially on the stack, but come from evaluating the |array| and |index| fields
of the |K_type_coefficient| expression.

@h "axis.h" // for |K_type_pol_coefficient|

@< Function def... @>=
void K_type_pol_coefficient::evaluate(level l) const
{ shared_K_type_pol m = (array->eval(),get<K_type_pol_value>());
  shared_K_type p = (index->eval(),get<K_type_value>());
  if (m->rf!=p->rf and m->rf->val!=p->rf->val)
    // test like |real_form_new_wrapper| does
    throw runtime_error @|
      ("Real form mismatch when subscripting KTypePol value");
  test_final(*p,"In subscription of KTypePol value");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<split_int_value>(m->val[p->val]));
}

@ For modifying an individual coefficient in a |K_type_pol_value| (possibly
adding or removing a term in the process), a method |assign_coef| is provided.
It will be called from the \.{axis} compilation unit (the programming language
interpreter). The implementation is quite simple; the only subtlety is that for
the probably common case of setting a coefficient to~$0$, a call to the method
|clear_coefficient| is made, which is more appropriate for this purpose than
calling |set_coefficient| with a zero value.

@< Function def...@>=
void K_type_pol_value::assign_coef
  (const K_type_value& t, const Split_integer& c)
{
  test_final(t,"In coefficient assignment for KTypePol value");
  if (c.is_zero())
    val.clear_coefficient(t.val);
  else
    val.set_coefficient(t.val.copy(),c);
}

@ The main operations for $K$-type polynomials are addition and subtraction of
$K$-types to or from them. Although only one $K$-type is provided here, we
must apply |finals_for| to expand it into a linear combination of final
$K$-types, and call |accumulator.add_term| with each of them.

@< Local function def...@>=
void add_K_type_wrapper(expression_base::level l)
{ shared_K_type p = get<K_type_value>();
  own_K_type_pol accumulator = get_own<K_type_pol_value>();
  if (accumulator->rf!=p->rf)
    throw runtime_error @|
      ("Real form mismatch when adding a KType to a KTypePol");
  if (l==expression_base::no_value)
    return;
  auto finals = p->rc().finals_for(p->val.copy());
    for (auto it = finals.wbegin(); it!=finals.wend(); ++it)
       accumulator->val.add_term(std::move(it->first),Split_integer(it->second));
  push_value(std::move(accumulator));
}

void subtract_K_type_wrapper(expression_base::level l)
{ shared_K_type p = get<K_type_value>();
  own_K_type_pol accumulator = get_own<K_type_pol_value>();
  if (accumulator->rf!=p->rf)
    throw runtime_error @|
      ("Real form mismatch when subtracting a KType from a KTypePol");
  if (l==expression_base::no_value)
    return;
  auto finals = p->rc().finals_for(p->val.copy());
  for (auto it = finals.wbegin(); it!=finals.wend(); ++it)
    accumulator->val.add_term(std::move(it->first),Split_integer(-it->second));
  push_value(std::move(accumulator));
}

@ More generally than adding or subtracting, we can incorporate a term with
specified coefficient; this is almost the same as the previous two functions.

@< Local function def...@>=

void add_K_type_term_wrapper(expression_base::level l)
{ push_tuple_components(); // second argument is a pair |(coef,p)|
  auto p = get<K_type_value>();
  Split_integer coef=get<split_int_value>()->val;
  auto accumulator = get_own<K_type_pol_value>();

if (accumulator->rf!=p->rf)
    throw runtime_error@|("Real form mismatch when adding a term to a K_type");
  if (l==expression_base::no_value)
    return;
@)
  auto finals = p->rc().finals_for(p->val.copy());
  for (auto it=finals.wbegin(); not finals.at_end(it); ++it)
    accumulator->val.add_term(std::move(it->first),coef*it->second);
  push_value(std::move(accumulator));
}

@ When producing a $K$-type polynomial value, the user will most often want to
provide all terms in a list at once, since repeatedly adding contributions to a
named variable will, currently and without special trickery, lead to above calls
to |get_own| having to duplicate the accumulator before each addition, since the
variable still holds a reference to the value on the stack. The following
built-in function, which will be bound to the operator \.+, does this, but does
require an initial (in practice often empty) $K$-type polynomial as first
operand to provide the real form, in case the list of terms is empty. We try to
modify the initial argument in place (since of course it does not have to start
out empty), but in practice it is more likely that the second argument is
unshared and voluminous (being the result of a |for| loop). For that reason we
do effort here to at least move the individual $K$-types into the returned
polynomial, rather than copying them, if we can. This can only be achieved by
calling |force_own| at several places, so that copies are made if parts of the
term list are shared and cannot be moved from. The most likely cause of the
latter situation is if there is sharing already at the top level of the term
list; if this is detected, we apply different code that does not attempt to get
unique ownership, but instead copies the $K$-types at the point where they get
converted into final $K$-type for the accumulator. Comparing the cases below is
instructive to see what needs to change if one sets out to reuse component
values of a row-of-tuples argument.

@< Local function... @>=
void add_K_type_termlist_wrapper(expression_base::level l)
{ shared_row r = get<row_value>();
  own_K_type_pol accumulator = get_own<K_type_pol_value>();
  if (l==expression_base::no_value)
    return;
@)
  if (r.unique())
  { auto own_r = std::const_pointer_cast<row_value>(r);
    for (auto it=own_r->val.begin(); it!=own_r->val.end(); ++it)
    { auto tup = force_own<tuple_value>(std::move(*it));
      Split_integer coef=force<split_int_value>(tup->val[0].get())->val;
      auto t = force_own<K_type_value>(std::move(tup->val[1]));
      if (accumulator->rf!=t->rf)
        throw runtime_error@|
          ("Real form mismatch when adding terms to a K_type");
       auto finals = t->rc().finals_for(std::move(t->val));
       for (auto it=finals.wbegin(); not finals.at_end(it); ++it)
         accumulator->val.add_term(std::move(it->first),coef*it->second);
     }
  }
  else // there is top level sharing, so avoid duplicating all alon
    for (auto it=r->val.cbegin(); it!=r->val.cend(); ++it)
    { auto tup = force<tuple_value>(it->get());
      Split_integer coef=force<split_int_value>(tup->val[0].get())->val;
      auto t = force<K_type_value>(tup->val[1].get());
      if (accumulator->rf!=t->rf)
        throw runtime_error@|
          ("Real form mismatch when adding terms to a K_type");
       auto finals = t->rc().finals_for(t->val.copy());
       for (auto it=finals.wbegin(); not finals.at_end(it); ++it)
         accumulator->val.add_term(std::move(it->first),coef*it->second);
     }
  push_value(std::move(accumulator));
}

@ Naturally we also want to define addition and subtraction of two $K$-type
polynomials.

@< Local function... @>=
void add_K_type_pols_wrapper(expression_base::level l)
{
  own_K_type_pol addend = get_own<K_type_pol_value>();
  own_K_type_pol accumulator = get_own<K_type_pol_value>();
  if (accumulator->rf!=addend->rf)
    throw runtime_error @|("Real form mismatch when adding two K_types");
  if (l!=expression_base::no_value)
  @/{@; accumulator->val += std::move(addend->val);
    push_value(std::move(accumulator));
  }
}
@)
void subtract_K_type_pols_wrapper(expression_base::level l)
{
  own_K_type_pol subtrahend = get_own<K_type_pol_value>();
  own_K_type_pol accumulator = get_own<K_type_pol_value>();
  if (accumulator->rf!=subtrahend->rf)
    throw runtime_error@|("Real form mismatch when subtracting two K_types");
  if (l!=expression_base::no_value)
  @/{@; accumulator->val -= std::move(subtrahend->val);
    push_value(std::move(accumulator));
  }
}

@ Scalar multiplication potentially makes coefficients zero, in which case the
corresponding terms need to be removed to preserve the invariant that no zero
terms are stored in a $K$-type polynomial. For integer multiplication we just
need to check for multiplication by $0$, and produce an empty module when this
happens. Because the integer is not on the stack top, this requires a somewhat
unusual manoeuvre. The case of multiplication by zero needs to be handled
separately, since we cannot allow introducing terms with zero coefficients. It
could have been handled more easily though, by testing the factor~|c| just
before the |for| loop, and performing |m->erase()| instead if |c==0|; this is
what we used to do. However that might involve duplicating the $K$-type
polynomial and then erasing the copy, which is inefficient, and now avoided.
This might seem a rare case, but it is not really: often functions handling
a \.{KTypePol} argument $P$ need to start with an empty module for the same real
form; writing $0*P$ is quite a convenient way to achieve this.

@< Local function... @>=

void int_mult_K_type_pol_wrapper(expression_base::level l)
{ int c =
    force<int_value>(execution_stack[execution_stack.size()-2].get())
    // value below top
    ->int_val();
  if (c==0) // then do multiply by $0$ efficiently:
  { shared_K_type_pol m = get<K_type_pol_value>();
      // |m| is needed for |m->rc()|
    pop_value();
    if (l!=expression_base::no_value)
    @/push_value@|(std::make_shared<K_type_pol_value>
        (m->rf,K_repr::K_type_pol()));
  }
  else
  { own_K_type_pol m = get_own<K_type_pol_value>();
     // will modify our copy now
    pop_value();
    if (l!=expression_base::no_value)
    { for (auto& term : m->val)
        term.second *= c;
      push_value(std::move(m));
    }
  }
}

@ Matters are similar but somewhat subtler for scalar multiplication by split
integers, because these have zero divisors. We make a $4$-way branch depending
on whether the split integer multiplicand has zero evaluation are $s=1$ and/or
at $s=-1$: if either of these is the case we have a zero divisor, and may expect
some terms to be killed by the multiplication, and if both hold we are
multiplying by $0$. In presence of a nonzero zero divisor, we do not call
|get_own| to true to modify the argument in place, as it seems more efficient to
just reconstruct a potentially small product separately. As in the previous
function we pick up the multiplicand (which is the first argument) when it is
not at the top of the stack; later on, we do not pop it off the stack either,
but rather overwrite it by move-assigning the returned |K_type_pol_value| to
the stack top location |execution_stack.back()|. Since no tests are performed
here, the |no_value| case can be handled right at the beginning, by simply
dropping the arguments from the stack when it applies.

@< Local function... @>=

void split_mult_K_type_pol_wrapper(expression_base::level l)
{  if (l==expression_base::no_value)
   {@; pop_value();
       pop_value();
       return;
   }
   auto c =
    force<split_int_value>(execution_stack[execution_stack.size()-2].get())
    // below top
    ->val;
  if (c.s_to_1()==0)
  // then coefficient is multiple of $1-s$; don't try to modify module in place
  { shared_K_type_pol m = get<K_type_pol_value>();
    K_repr::K_type_pol result;
    if (c.s_to_minus_1()!=0) // otherwise coefficient is $0$ and keep null module
      for (const auto& term: m->val)
        if (term.second.s_to_minus_1()!=0)
          result.add_term(term.first.copy(),term.second*c);
    execution_stack.back()=  // replace stack top
        std::make_shared<K_type_pol_value>(m->rf,std::move(result));
    return;
  }
@)
  else if (c.s_to_minus_1()==0)
  // then coefficient is multiple of $1+s$; don't modify in place
  { shared_K_type_pol m = get<K_type_pol_value>();
    K_repr::K_type_pol result;
    for (const auto& term: m->val)
      if (term.second.s_to_1()!=0)
        result.add_term(term.first.copy(),term.second*c);
    execution_stack.back()= // replace stack top
        std::make_shared<K_type_pol_value>(m->rf,std::move(result));
    return;
  }
  // now we are multiplying by a non zero-divisor
  own_K_type_pol m = get_own<K_type_pol_value>();
     // will modify our copy now
  for (auto& term : m->val)
        term.second *= c; // multiply coefficients by split integer
  execution_stack.back()=std::move(m);  // replace stack top
}

@ For a nonzero $K$-type polynomial, it is useful to be able to select some term
that is present, without looping over all terms. The most useful choice it the
final term, be we allow taking the first term as well.

@< Local function... @>=
void last_K_type_term_wrapper (expression_base::level l)
{ shared_K_type_pol m = get<K_type_pol_value>();
  if (l==expression_base::no_value)
    return;
@)
  if (m->val.is_zero())
    throw runtime_error("Empty KTypePol has no last term");
  const auto& term = *m->val.rbegin();
  push_value(std::make_shared<split_int_value>(term.second));
  push_value(std::make_shared<K_type_value>
	(m->rf,term.first.copy()));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}
@)
void first_K_type_term_wrapper (expression_base::level l)
{ shared_K_type_pol m = get<K_type_pol_value>();
  if (l==expression_base::no_value)
    return;
@)
  if (m->val.is_zero())
    throw runtime_error("Empty KTypePol has no first term");
  const auto& term = *m->val.begin();
  push_value(std::make_shared<split_int_value>(term.second));
  push_value(std::make_shared<K_type_value>
  	(m->rf,term.first.copy()));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ Sometimes we want to ignore $K$-types whose height exceeds a given limit.
Since terms are sorted by height this can be done fairly efficiently, and the
built-in |truncate_above_height| will do this.

@< Local function... @>=
void truncate_K_type_poly_above_wrapper (expression_base::level l)
{ int arg = get<int_value>()->int_val();
  shared_K_type_pol P = get<K_type_pol_value>();
  if (l==expression_base::no_value)
    return;
@)
  sl_list<K_repr::K_type_pol::value_type> L;
  if (arg>=0)
  { unsigned bound = arg; // use |unsigned| value for height comparison
    for (const auto& term : P->val)
      if (term.first.height()<=bound)
        L.emplace_back(term.first.copy(),term.second);
      else
        goto wrap_up;
     push_value(P); // if loop ran to completion, just return |P|
     return; // and skip returning |L|
  }
  wrap_up:
    push_value(std::make_shared<K_type_pol_value>@|
      (P->rf,K_repr::K_type_pol(std::move(L).to_vector(),false)));
      // no need to sort again
}


@*2 Computing with $K$-types.
%
A main application of $K$-types is branching to~$K$: the decomposition of a
standard representation into $K$-types. Because this decomposition is infinite,
we allow computing it up to a given limit in the height of the $K$-types. The
rewriting of (possibly non standard) $K$-types into linear combinations of final
$K$-types used to be a separate functionality provided here, but it is now
incorporated into the general behaviour of $K$-type polynomials. What we are
left with is two auxiliary functions, making these operations that are used
internally by branching also available separately for the user, and the actual
branching function.

The first auxiliary function |KGP_sum|, produces a set of values from which the
second auxiliary function |K_type_formula| builds a finite $K$-type polynomial
that will be used in the long division. The values of the |KGP_sum| are
represented as $K$-types with a sign attached, but since the $K$-types need not
be standard, representing this result as a $K$-type polynomial would not
be right (regardless of whether it would be mathematically correct or not, the
conversion to final $K$-types is not done in the internal implementation of
|K_type_formula|, so it would be unhelpful to do it in the corresponding
built-in function); therefore we return a list of pairs of an integer (sign) and
a $K$-type.

@h <limits> // for |max|

@< Local function def...@>=
void KGP_sum_wrapper(expression_base::level l)
{ shared_K_type p = get<K_type_value>();
@/const Rep_context rc = p->rc();
  auto srk = p->val.copy();
  if (not rc.is_semifinal(srk))
    throw runtime_error@|("K-type has parity real roots (so not semifinal)");
  if (l==expression_base::no_value)
    return;
@)
  auto length = rc.kgb().length(srk.x());
  auto list = rc.KGP_set(srk);
  own_row result = std::make_shared<row_value>(list.size());
  auto res_p=result->val.begin();
  for (auto&& t : list)
  {
    auto tup = std::make_shared<tuple_value>(2);
    auto dl = length - rc.kgb().length(t.x());
    tup->val[0] = std::make_shared<int_value>(dl%2==0 ? 1 : -1);
    tup->val[1] = std::make_shared<K_type_value>(p->rf,std::move(t));
    *res_p++ = std::move(tup);
  }
  push_value(std::move(result));
}

@ While the |KGP_sum| returns a result only depending on the initial $K$-type,
the following function |K_type_formula| allows pruning its result to terms that
will be useful for branching up to a given height, thus in many cases allowing
for much more efficient computation and handling afterwards. Therefore we
provide a height bound argument, which may however be given as a negative number
to indicate that no pruning should be done, an the full (still finite) $K$-type
formula should be computed. Since there is only one method
|Rep_context::K_typ_formula|, we here transform such a negative integer into the
highest possible value of the unsigned type |repr::level|.

@< Local function def...@>=
void K_type_formula_wrapper(expression_base::level l)
{ int bound = get<int_value>()->int_val();
  shared_K_type p = get<K_type_value>();
@/const Rep_context rc = p->rc();
  auto srk = p->val.copy();
  if (not rc.is_semifinal(srk))
    throw runtime_error@|("K-type has parity real roots (so not semifinal)");
  if (l==expression_base::no_value)
    return;
@)
  repr::level h = bound<0 ? std::numeric_limits<repr::level>::max() : bound;
  push_value(std::make_shared<K_type_pol_value>(p->rf,rc.K_type_formula(srk,h)));
}

@ Here is the function that implements branching. Because it is a long division
(of formal power series) type operation, where a remainder expression is
repeatedly modified by subtracting off contributions that can be determined by
the leading (lowest) terms of the current remainder, it is implemented directly
for \.{KTypePol} arguments rather than for individual \.{KType} values; if
needed a \.{KType} can be implicitly converted to \.{KTypePol} first. Also a
bound on the height is an obligatory argument, since only with such a height
bound can the computation terminate.

uses the |K_type_formula| internally.
@< Local function def...@>=
void branch_wrapper(expression_base::level l)
{ int arg = get<int_value>()->int_val();
  own_K_type_pol P = get_own<K_type_pol_value>();
  if (arg<0)
    throw runtime_error("Maximum level in branch cannot be negative");
  if (l==expression_base::no_value)
    return;
@)
  unsigned int bound=arg;
  const Rep_context rc = P->rc();
  auto result = rc.branch(std::move(P->val),bound);
  push_value(std::make_shared<K_type_pol_value>(P->rf,std::move(result)));
}

@ Finally we install everything related to $K$-types.
@< Install wrapper functions @>=
install_function(K_type_wrapper,"K_type","(KGBElt,vec->KType)");
install_function(unwrap_K_type_wrapper,@|"%","(KType->KGBElt,vec)");
install_function(real_form_of_K_type_wrapper,@|"real_form"
		,"(KType->RealForm)");
install_function(K_type_height_wrapper,@|"height" ,"(KType->int)");
install_function(K_type_eq_wrapper,@|"=", "(KType,KType->bool)");
install_function(K_type_neq_wrapper,@|"!=", "(KType,KType->bool)");
install_function(K_type_equivalent_wrapper,@|"equivalent","(KType,KType->bool)");
install_function(K_type_is_standard_wrapper,@|"is_standard" ,"(KType->bool)");
install_function(K_type_is_dominant_wrapper,@|"is_dominant" ,"(KType->bool)");
install_function(K_type_is_zero_wrapper,@|"is_zero" ,"(KType->bool)");
install_function(K_type_is_semifinal_wrapper,@|"is_semifinal" ,"(KType->bool)");
install_function(K_type_is_final_wrapper,@|"is_final" ,"(KType->bool)");
install_function(K_type_dominant_wrapper,@|"dominant" ,"(KType->KType)");
install_function(to_canonical_fiber_wrapper,@|"to_canonical_fiber"
		,"(KType->KType)");
install_function(K_type_normal_wrapper,@|"normal" ,"(KType->KType)");
install_function(K_type_theta_stable_wrapper,@|"theta_stable"
		,"(KType->KType)");
@)
install_function(K_type_pol_wrapper,@|"null_K_module","(RealForm->KTypePol)");
install_function(real_form_of_K_type_pol_wrapper,@|"real_form"
		,"(KTypePol->RealForm)");
install_function(K_type_pol_unary_eq_wrapper,@|"=","(KTypePol->bool)");
install_function(K_type_pol_unary_neq_wrapper,@|"!=","(KTypePol->bool)");
install_function(K_type_pol_eq_wrapper,@|"=","(KTypePol,KTypePol->bool)");
install_function(K_type_pol_neq_wrapper,@|"!=","(KTypePol,KTypePol->bool)");
install_function(K_type_pol_size_wrapper,@|"#","(KTypePol->int)");
install_function(add_K_type_wrapper,@|"+","(KTypePol,KType->KTypePol)");
install_function(subtract_K_type_wrapper,@|"-","(KTypePol,KType->KTypePol)");
install_function(add_K_type_term_wrapper,@|"+"
		,"(KTypePol,(Split,KType)->KTypePol)");
install_function(add_K_type_termlist_wrapper,@|"+"
		,"(KTypePol,[(Split,KType)]->KTypePol)");
install_function(add_K_type_pols_wrapper,@|"+"
		,"(KTypePol,KTypePol->KTypePol)");
install_function(subtract_K_type_pols_wrapper,@|"-"
		,"(KTypePol,KTypePol->KTypePol)");
install_function(int_mult_K_type_pol_wrapper,@|"*"
		,"(int,KTypePol->KTypePol)");
install_function(split_mult_K_type_pol_wrapper,@|"*"
		,"(Split,KTypePol->KTypePol)");
install_function(last_K_type_term_wrapper,@|"last_term"
		,"(KTypePol->Split,KType)");
install_function(first_K_type_term_wrapper,@|"first_term"
		,"(KTypePol->Split,KType)");
install_function(truncate_K_type_poly_above_wrapper,@|"truncate_above_height"
		,"(KTypePol,int->KTypePol)");
@)
install_function(KGP_sum_wrapper,@|"KGP_sum","(KType->[int,KType])");
install_function(K_type_formula_wrapper,@|"K_type_formula"
		,"(KType,int->KTypePol)");
install_function(branch_wrapper,@|"branch" ,"(KTypePol,int->KTypePol)");

@*1 Standard module parameters.
%
We implement a data type for holding parameters that represent a standard module
or, depending on the context, its unique irreducible quotient. We call these
values standard module parameters, often shortened to simply parameters, and in
the \.{atlas} language designate the corresponding basic type as \.{Param}. Such
a parameter is defined by a triple $(x,\lambda,\nu)$ where $x$ is a KGB element,
$\lambda$ is a weight in the coset $\rho+X^*$ whose value is relevant only
modulo the sub-lattice $(1-\theta_x)X^*$ where $\theta_x$ is the involution
associated to $x$, and $\nu$ is a rational weight in the kernel of $1+\theta_x$.
Such parameters are stored in instances of the class |StandardRepr|, which is
defined in the file \.{repr.h}.

@*2 Class definition.
Like for KGB elements, we maintain a shared pointer to the real form value, so
that it will be assured to survive as long as parameters for it exist, and we
can access notably the |Rep_context| that it provides.

@< Type definitions @>=
struct module_parameter_value : public value_base
{ shared_real_form rf;
  StandardRepr val;
@)
  module_parameter_value(const shared_real_form& form, const StandardRepr& v)
  : rf(form), val(v) @+{}
  module_parameter_value(const shared_real_form& form, StandardRepr&& v)
  : rf(form), val(std::move(v)) @+{}
  ~module_parameter_value() @+{}
@)
  virtual void print(std::ostream& out) const;
  static const char* name() @+{@; return "module parameter"; }
  module_parameter_value (const module_parameter_value& ) = default;
    // we use |get_own<module_parameter_value>|

@)
  const Rep_context& rc() const @+{@; return rf->rc(); }
  Rep_table& rt() const @+{@; return rf->rt(); }
};
@)
typedef std::unique_ptr<module_parameter_value> module_parameter_ptr;
typedef std::shared_ptr<const module_parameter_value> shared_module_parameter;
typedef std::shared_ptr<module_parameter_value> own_module_parameter;

@ When printing a module parameter, we shall indicate a triple $(x,\lambda,\nu)$
that defines it. Since we shall need to print |StandardRepr| values in other
contexts as well, an auxiliary output function |repr::print_stdrep| is defined
in \.{basic\_io.h}, which takes and additional |Rep_context| argument; we shall
call that function from |module_parameter_value::print|. By choosing a name for
the auxiliary function different from |print|, we avoid having that call being
mistaken for a recursive call.

The virtual method |module_parameter_value::print|, is used when printing a
value of type \.{Param} (as opposed to for instance printing a term of
a \.{ParamPol}). Here we prefix the parameter text proper with additional
information about the parameter that may be relevant to the \.{atlas} user.

@< Function definition... @>=
void module_parameter_value::print(std::ostream& out) const
{ out << @< Expression for adjectives that apply to a module parameter @>@;@;;
  print_stdrep(out << ' ',val,rc());
}

@ We provide one of the adjectives ``non-standard'' (when $\gamma$ and therefore
$\lambda$ fails to be imaginary-dominant), ``non-dominant'' (when $\gamma$ fails
to be dominant), ``zero'' (the standard module vanishes due to the singular
infinitesimal character, namely by the presence of a singular compact
simply-imaginary root), ``non-final'' (the standard module is non-zero, but can
be expressed in terms of standard modules at more compact Cartans using a
singular real root satisfying the parity condition), ``non-normal'' (the
parameter differs from its normal form; when we come to this point it implies
there is a complex singular descent), or finally ``final'' (the good ones that
could go into a \.{ParamPol} value; the condition |is_final| should apply,
though it is not tested here).

@< Expression for adjectives that apply to a module parameter @>=
( not rc().is_standard(val) ? "non-standard"
@|: not rc().is_dominant(val) ? "non-dominant"
@|: not rc().is_nonzero(val) ? "zero"
@|: not rc().is_semifinal(val) ? "non-final"
@|: not rc().is_normal(val) ? "non-normal"
@|: "final")

@ To make a module parameter, one should provide a KGB element~$x$, an
integral weight $\lambda-\rho$, and a rational weight~$\nu$. Since only its
projection on the $-\theta_x$-stable subspace is used, one might specify the
infinitesimal character $\gamma$ in the place of $\nu$.

@f nu nullptr

@< Local function def...@>=
void module_parameter_wrapper(expression_base::level l)
{ shared_rational_vector nu(get<rational_vector_value>());
  shared_vector lam_rho(get<vector_value>());
  shared_KGB_elt x = get<KGB_elt_value>();
  if (nu->val.size()!=lam_rho->val.size()
      or nu->val.size()!=x->rf->val.rank())
  { std::ostringstream o;
    o << "Rank mismatch: (" @|
        << x->rf->val.rank() << ','
	<< lam_rho->val.size() << ',' << nu->val.size() << ")";
    throw runtime_error(o.str());
  }
@.Rank mismatch@>
  if (l!=expression_base::no_value)
    push_value(std::make_shared<module_parameter_value> @| (x->rf,
      x->rf->rc().sr(x->val,lam_rho->val,nu->val)));
}

@*2 Functions operating on module parameters.
%
The unwrapping function below, which we shall bind to the operator~`|%|',
transforms a parameter value into a triple of values $(x,\lambda-\rho,\gamma)$
that defines it, where $\gamma$ is taken to be (a representative of) the
infinitesimal character. (The function used to produce just the $\nu$ part
of the infinitesimal character as third component, but in practice obtaining the
infinitesimal character directly turns out to often be more useful. If needed,
$\nu$ is easily computed as $\nu={1+\theta_x\over2}\gamma$; the opposite
conversion, which would be needed to be programmed if we returned $\nu$ here,
requires more work.) While triple returned here is not unique, since $\lambda$
is determined only modulo $(1-\theta_x)X^*$, we choose a unique representative
for~|lambda|, which choice is in fact determined by the implementation of
|StandardRepr|.

We also provide a function directly extracting the real form from a module
parameter (which would otherwise require extracting it from the |x| component).

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

@ We provide explicit conversions between a module parameter and a $K$-type,
restricting to $K$ (and ignoring $\nu$) in one direction, and extending with
$\nu=0$ in the opposite direction.

@< Local function def...@>=
void param_to_K_type_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l==expression_base::no_value)
    return;
  const auto& rf=p->rf;
  push_value(std::make_shared<K_type_value> @|
    (rf,rf->rc().sr_K(p->val)));
}
void K_type_to_param_wrapper(expression_base::level l)
{ shared_K_type p = get<K_type_value>();
  if (l==expression_base::no_value)
    return;
  const auto& rf=p->rf;
  push_value(std::make_shared<module_parameter_value> @|
    (rf,rf->rc().sr(p->val)));
}

@ Another simple function is computing the height of a parameter (ignoring the
$\nu$ component). This is the same height as that of the $K$-type to which the
parameter can be restricted, and to the height displayed when
printing \.{ParamPol} values. The height is stored inside each |StandardRepr|
value, so it can be obtained from there there; this is somewhat more efficient
than constructing a $K$-type from the parameter and taking the height of that.

@< Local function def...@>=
void parameter_height_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(p->val.height()));
}


@ While the equality and inequality operators for module parameters test for
strict equality of all components (including of the real forms, which the method
|StandardRepr::operator==| must simply assume to be equal), a separate function
tests for \emph{equivalence}; this is a weaker condition when both parameters
are standard. (When at least one parameter fails to be standard, the equivalence
test reverts to testing strict equality.) Equivalence of standard parameters
amounts to testing for equality after the parameters are made dominant (at least
that claim was not contested at the time of writing this). This test will be
bound to the name |equivalent|. Unlike the equality tests, it \emph{requires}
the parameters to be associated to the same real form, giving a runtime error
(rather than returning false) if not; this avoids confusion if there were some
subtle difference of real forms for otherwise similar parameters. If some
operation is used to produce parameters that may of may not be associated to the
same real form, then one should test those forms for equality before testing
equivalence of the parameters.

@< Local function def...@>=

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

@ Here are some more attributes, in the form of predicates.

@< Local function def...@>=
void is_standard_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(whether(p->rc().is_standard(p->val)));
}

void is_dominant_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(whether(p->rc().is_dominant(p->val)));
}

void is_zero_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(whether(not p->rc().is_nonzero(p->val)));
}

void is_semifinal_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(whether(p->rc().is_semifinal(p->val)));
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


@< Local function def...@>=
void parameter_dominant_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  if (l!=expression_base::no_value)
  {@; p->rc().make_dominant(p->val);
    push_value(std::move(p));
  }
}

void parameter_normal_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  if (l!=expression_base::no_value)
  {@; p->rc().normalise(p->val);
    push_value(std::move(p));
  }
}


@ While parameters can be used to compute blocks of (other) parameters, it can
be useful to have available the basic operations of cross actions and Cayley
transforms on individual parameters without going through the construction of
an entire block. The library provides methods that do this computation directly
in the |StandardRepr| format (calling |make_dominant| on it first). In case they
find that the type of the simple reflection for the integral system is not right
for the Cayley transform demanded (imaginary noncompact respectively real
parity) they throw a |Cayley_error| value, which is caught here and translated
in make the whole function a no-operation (so that the caller gets an occasion
to test the condition).

Like for KGB elements there is the possibility of double values for the Cayley
transforms, in either direction. The ``solution'' to this difficulty is the same
here: the user can find out by herself about a possible second image by applying
a cross action to the result. In the current case this approach has in fact
already been adopted in the method that is called here, which presents a
single-valued interface to the Cayley transform.

@< Local function def...@>=
void parameter_cross_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  int s = get<int_value>()->int_val();
  unsigned int r =
    rootdata::integrality_rank(p->rf->val.root_datum(),p->val.gamma());
  if (static_cast<unsigned>(s)>=r)
  { std::ostringstream o;
    o << "Illegal simple reflection: " << s << ", should be <" << r;
    throw runtime_error(o.str());
  }
  if (l!=expression_base::no_value)
    push_value(std::make_shared<module_parameter_value>
		(p->rf,p->rc().cross(s,p->val)));
}
@)
void parameter_Cayley_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  int s = get<int_value>()->int_val();
  unsigned int r =
    rootdata::integrality_rank(p->rf->val.root_datum(),p->val.gamma());
  if (static_cast<unsigned>(s)>=r)
  { std::ostringstream o;
    o << "Illegal simple reflection: " << s @|
      << ", should be <" << r;
    throw runtime_error(o.str());
  }
  if (l==expression_base::no_value)
    return;

  try {
    push_value(std::make_shared<module_parameter_value>
		(p->rf,p->rc().Cayley(s,p->val)));
  }
  catch (error::Cayley_error& e) // ignore undefined Cayley transforms
  {@;
    push_value(std::move(p));
  }
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
    push_value(std::move(p));
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
  test_compatible(p->rc().inner_class(),delta);
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
  RatNumList rp = p->rc().reducibility_points(p->val);
      // method normalises rationals
  own_row result = std::make_shared<row_value>(rp.size());
  for (unsigned int i=0; i<rp.size(); ++i)
    result->val[i]=std::make_shared<rat_value>(rp[i]);
  push_value(std::move(result));
}

@ Scaling the continuous component~$\nu$ of a parameter is an important
ingredient for calculating signatures of Hermitian forms. This was for a long
time done with a user defined function, but having this built in is more
efficient. For scaling by~$0$ (called ``deformation to $\nu=0$''), the
construction of a $K$-type from a parameter is often a better alternative, as
then the result type makes clear that information about $\nu$ has been erased.

@< Local function def...@>=
void scale_parameter_wrapper(expression_base::level l)
{ shared_rat f = get<rat_value>();
  own_module_parameter p = get_own<module_parameter_value>();
  if (l!=expression_base::no_value)
@/{@; p->val = p->rc().scale(p->val,f->rat_val());
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
{ if (p.rc().is_standard(p.val))
    return;
  std::ostringstream os; p.print(os << descr << ":\n  ");
  os << "\n  Parameter not standard";
  throw runtime_error(os.str());
}
@)
void test_final(const K_type_value& p, const char* descr)
{ if (p.rc().is_final(p.val))
    return;
  std::string reason;
  if (not p.rc().is_standard(p.val))
    reason = "not standard";
  if (not p.rc().is_dominant(p.val))
    reason = "not dominant";
  else if (not p.rc().is_nonzero(p.val))
    reason = "zero";
  else if (not p.rc().is_semifinal(p.val))
    reason = "not semifinal";
  else if (not p.rc().is_normal(p.val)) // this predicate must come last
    reason = "not normal";
  else throw logic_error("Unknown obstruction to K-type finality");
  std::ostringstream os; p.print(os << descr << ":\n  ");
@/os << "\n  K-type is " << reason;
  throw runtime_error(os.str());
}
void test_final(const module_parameter_value& p, const char* descr)
{ std::string reason;
  bool OK = p.rc().is_dominant(p.val);
  if (not OK)
    reason = "not dominant";
  else if (not (OK=p.rc().is_normal(p.val)))
    reason = "not normal";
  else if (not(OK = p.rc().is_nonzero(p.val)))
    reason = "zero";
  else if (not(OK = p.rc().is_semifinal(p.val)))
    reason = "not semifinal";
  else return; // nothing to report
  std::ostringstream os; p.print(os << descr << ":\n  ");
@/os << "\n  Parameter is " << reason;
  throw runtime_error(os.str());
}

@ Here is the first block generating function, which just reproduces to output
of the \.{full\_block} command in the \.{Fokko} program, and a variation for
partial blocks.

@< Local function def...@>=
void print_c_block_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  BlockElt init_index; // will hold index in the block of the initial element
  blocks::common_block& block = p->rt().lookup_full_block(p->val,init_index);
  RatWeight diff = p->rc().offset(p->val, block.representative(init_index));
  *output_stream << "Parameter defines element " << init_index
               @|<< " of the following common block:" << std::endl;
  block.shift(diff);
  block.print_to(*output_stream,block.singular(p->val.gamma()));
    // print block using involution expressions
  block.shift(-diff);
  if (l==expression_base::single_value)
    wrap_tuple<0>(); // |no_value| needs no special care
}

void print_pc_block_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  BlockElt init_index; // will hold index in the block of the initial element
  blocks::common_block& block = p->rt().lookup(p->val,init_index);
  RatWeight diff = p->rc().offset(p->val, block.representative(init_index));
  BitMap less = block.bruhatOrder().poset().below(init_index);
  if (less.full())
  {
    if (init_index+1<block.size())
      *output_stream << "Elements <= " << init_index << " of following block\n";
  }
  else
  {
    *output_stream << "Subset {";
    for (auto n : less)
      *output_stream << n << ',';
    *output_stream << init_index << "} in the following common block:\n";
  }
  block.shift(diff);
  block.print_to(*output_stream,block.singular(p->val.gamma()));
    // print using involution expressions
  block.shift(-diff);
  if (l==expression_base::single_value)
    wrap_tuple<0>(); // |no_value| needs no special care
}

@ More interesting than printing the block is to return is to the user as a
list of parameter values. The following function does this, and adds as a
second result the index that the original parameter has in the resulting
block. If not final, the original parameter will be absent, and the second
value returned~$-1$.

@< Local function def...@>=
void common_block_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt start; // will hold index in the block of the initial element
  auto& block = p->rt().lookup_full_block(p->val,start);
  RatWeight diff = p->rc().offset(p->val,block.representative(start));
  const auto& gamma = p->val.gamma();
  { const RankFlags singular = block.singular(gamma);
    int start_pos = -1;
    own_row param_list = std::make_shared<row_value>(0);
    for (BlockElt z=0; z<block.size(); ++z)
      if (block.survives(z,singular))
      {
        if (z==start)
          start_pos=param_list->val.size();
        param_list->val.push_back @|
          (std::make_shared<module_parameter_value> @|
               (p->rf,p->rc().sr(block.representative(z),diff,gamma)));
      }
    push_value(std::move(param_list));
    push_value(std::make_shared<int_value>(start_pos));
  }
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ There are also a functions that compute just a partial block. We generate
the Bruhat interval inside the block (which might have more elements than just
the requested partial bock because of earlier computations) by completing the
downward closure of the Hasse relation (which is filled anyway by the partial
block |lookup| function); this avoids generating the full
|block.bruhatOrder().poset()| structure.

@< Local function def...@>=
void partial_common_block_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt start;
  blocks::common_block& block = p->rt().lookup(p->val,start);
  RatWeight diff = p->rc().offset(p->val,block.representative(start));
  const auto& gamma = p->val.gamma();
@)
  unsigned long n=block.size();
    // |unsigned long| type is imposed by |BitMap::back_up|
  BitMap subset(n);
  subset.insert(start);
  while (subset.back_up(n)) // compute downward closure
    for (BlockElt y : block.bruhatOrder().hasse(n))
      subset.insert(y);
@)
  { const RankFlags singular = block.singular(gamma);
    own_row param_list = std::make_shared<row_value>(0);
    for (auto z : subset)
      if (block.survives(z,singular))
        param_list->val.push_back @|
          (std::make_shared<module_parameter_value> @|
             (p->rf,p->rc().sr(block.representative(z),diff,gamma)));
    push_value(std::move(param_list));
  }
}

@ Knowing the length in its block of a parameter is of independent interest.
@< Local function def...@>=
void param_length_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot determine block for parameter length");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(p->rt().length(p->val)));
}

@ This function is similar to |KGB_Hasse_wrapper|, but generates the (full)
block on the fly.

@< Local function def...@>=
void block_Hasse_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt init_index; // will hold index in the block of the initial element
  blocks::common_block& block = p->rt().lookup_full_block(p->val,init_index);
  const BruhatOrder& Bruhat = block.bruhatOrder();
  auto n= block.size();
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(n,n,0));
  for (unsigned j=0; j<n; ++j)
    for (unsigned int i : Bruhat.hasse(j))
      M->val(i,j)=1;
  push_value(std::move(M));
}

@ Here is a version of the |block| command that also exports the table of
Kazhdan-Lusztig polynomials for the block. It exports 4 components: the list of
final parameters in the full block, the index of the initial parameter in the
list, a matrix of KL-polynomial indices, and a list of polynomials (as vectors).


@s IntPolEntry BlockElt

@< Local function def...@>=
void KL_block_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  test_standard(*p,"KL_block requires a standard parameter");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt start; // will hold index in the block of the initial element
  auto& block = p->rt().lookup_full_block(p->val,start);
  RatWeight diff = p->rc().offset(p->val,block.representative(start));
  const auto& gamma = p->val.gamma();
  const RankFlags singular = block.singular(gamma);
@)
  sl_list<BlockElt> survivors;
  BlockEltList loc(block.size(),UndefBlock);
  for (BlockElt z=0; z<block.size(); ++z)
    if (block.survives(z,singular))
    @/{@;
      loc[z] = survivors.size();
      survivors.push_back(z);
    }
@)
  const kl::KL_table& kl_tab = block.kl_tab(nullptr);
    // fill full block, not sharing polynomials, silently
  typedef polynomials::Polynomial<int> Pol;
  matrix::Matrix<Pol> M(survivors.size(),survivors.size(),Pol());
@/@< Condense the polynomials from |kl_tab| into the matrix |M| @>
@)
  @< Push list of parameters corresponding to |survivors| in |block|,
     with difference |diff| of $\gamma-\lambda$ values, and
     at infinitesimal character |gamma| @>
  if (loc[start]==UndefBlock)
    push_value(std::make_shared<int_value>(-1));
  else
    push_value(std::make_shared<int_value>(loc[start]));
@)
  @< Group distinct polynomials in |M| into a list, then push a version of |M|
  with polynomials replaced by there indices, and then push a list of the
  distinct polynomials @>
@)
  if (l==expression_base::single_value)
    wrap_tuple<4>();
}

@ Condensing the KL polynomials to the block at possibly singular infinitesimal
character~|gamma|, whose elements are the subset of |block| recorded in
|survivors|, means the following. One use only the columns indexed by
elements~|y| in |survivors|, and for each polynomial $P_{x,y}$ in that column
one computes the set of final elements~|f| to which $x$ contributes at~|gamma|
(which are in |survivors|); then $(-1)^{l(x)-l(f)}P_{x,y}$ is added to
$M_{f',y'}$ where the primes indicate re-indexing of final elements by position
in |survivors|, which renumbering is prepared in the vector~|loc|. Since we are
using unsigned-coefficient KL polynomials~|P| in computations of polynomials
that have signed integer coefficients, we need a conversion that can be obtained
by constructing a new polynomial |Pol(P)|.

@< Condense the polynomials from |kl_tab| into the matrix |M| @>=
{ auto start = survivors.begin();
  for (BlockElt x=0; x<block.size(); ++x)
  {
    while (not start.at_end() and *start<x)
      ++start; // ignore survivors less than |x|
    const auto finals = block.finals_for(x,singular);
    for (const auto f : finals)
    {
      auto i = loc[f];
      if (kl_tab.l(x,f)%2==0)
	for (auto it=start; not survivors.at_end(it); ++it)
	{ BlockElt y = *it;
          const auto& P = kl_tab.KL_pol(x,y);
          if (not P.isZero())
            M(i,loc[y]) += Pol(P);
        }
      else
	for (auto it=start; not survivors.at_end(it); ++it)
	{ BlockElt y = *it;
          const auto& P = kl_tab.KL_pol(x,y);
          if (not P.isZero())
            M(i,loc[y]) -= Pol(P);
        }
    }
  }
}

@ Here is another module that will be shared.
@< Push list of parameters corresponding to |survivors| in |block|,
   with difference |diff| of $\gamma-\lambda$ values, and
   at infinitesimal character |gamma| @>=
{ own_row param_list = std::make_shared<row_value>(0);
  param_list->val.reserve(survivors.size());
  for (BlockElt z : survivors)
    param_list->val.push_back (std::make_shared<module_parameter_value> @|
           (p->rf,p->rc().sr(block.representative(z),diff,gamma)));
  push_value(std::move(param_list));
}

@ And one more such shared module.
@< Group distinct polynomials in |M| into a list,... @>=
{ const auto n_survivors = survivors.size();
  std::vector<Pol> pool = { Pol(), Pol(1) };
  { HashTable<IntPolEntry,unsigned int> hash(pool);
    own_matrix M_ind = std::make_shared<matrix_value>(int_Matrix(n_survivors));
    for (BlockElt i = 0; i<n_survivors; ++i)
      for (BlockElt j = i+1; j<n_survivors; ++j)
         M_ind->val(i,j) = hash.match(M(i,j));
    push_value(std::move(M_ind));
  }
  @< Transfer the coefficient vectors of the polynomials from |pool| to an array,
     and push that array @>
}

@ Since the type |Polynomial<int>| has a method |data| that allows extracting its
coefficient vector, we can move without copying the essential information from
|pool| to the array of vector to be pushed here.

@< Transfer the coefficient vectors of the polynomials from |pool| to an array,
   and push that array @>=
{
  own_row polys = std::make_shared<row_value>(0);
  polys->val.reserve(pool.size());
  for (auto it=pool.begin(); it!=pool.end(); ++it)
    polys->val.emplace_back
      (std::make_shared<vector_value>(std::move(*it).data()));
  push_value(std::move(polys));
}

@ Here is a version of the |KL_block| that computes just for a partial block
and the Kazhdan-Lusztig polynomials for it. There are three components
in the value returned: the list of final parameters in the partial block (of
which the last one is the initial parameter), a matrix of KL-polynomial
indices, and a list of polynomials (as vectors).

@< Local function def...@>=
void partial_KL_block_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  test_standard(*p,"partial_KL_block requires a standard parameter");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt start; // will hold index in the block of the initial element
  auto& block = p->rt().lookup(p->val,start);
  RatWeight diff = p->rc().offset(p->val,block.representative(start));
  const auto& gamma = p->val.gamma();
@)
  unsigned long n=block.size();
    // |unsigned long| type is imposed by |BitMap::back_up|
  BitMap subset(n);
  subset.insert(start);
  while (subset.back_up(n)) // compute downward closure
    for (BlockElt y : block.bruhatOrder().hasse(n))
      subset.insert(y);
  @)
  const RankFlags singular = block.singular(gamma);
  sl_list<BlockElt> survivors;
  BlockEltList loc(block.size(),UndefBlock);
  for (BlockElt z : subset)
    if (block.survives(z,singular))
    @/{@;
      loc[z] = survivors.size();
      survivors.push_back(z);
    }
  @)

  const kl::KL_table& kl_tab = block.kl_tab(nullptr);
    // fill the partial block, not sharing polynomials, silently
  typedef polynomials::Polynomial<int> Pol;
  matrix::Matrix<Pol> M(survivors.size(),survivors.size(),Pol());
@/@< Condense the polynomials from |kl_tab| into the matrix |M| @>
@)
  @< Push list of parameters corresponding to |survivors| in |block|,... @>
  @< Group distinct polynomials in |M| into a list,... @>
@)
  if (l==expression_base::single_value)
    wrap_tuple<3>();
}

@ Here is a dual variation of |KL_block|. The main difference with that function
defined above is that we call the pseudo constructor |blocks::Bare_block::dual|
to transform |block| into its dual (represented as just a |blocks::Bare_block|
which is sufficient) before invoking the KL computations. The block is reversed
with respect to |block|, so for proper interpretation we reverse the list of
parameters returned, and this means that several other result components have to
be transformed as well. On the dual side there should be no ``condensing'' of the
polynomial matrix on the final elements, rather just an extraction of a
submatrix of polynomials at the corresponding indices; therefore we leave out
that part of the computation altogether.

@< Local function def...@>=
void dual_KL_block_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt start; // will hold index into |block| of the initial element
  auto& block = p->rt().lookup_full_block(p->val,start);
@/const auto& gamma = p->val.gamma();
  RatWeight diff = p->rc().offset(p->val,block.representative(start));
  auto dual_block = blocks::Bare_block::dual(block);
  const kl::KL_table& kl_tab = dual_block.kl_tab(nullptr);
  // fill entire KL table, don't share polys
@)
  sl_list<BlockElt> survivors; // indexes into |block|
  BlockEltList loc(block.size(),UndefBlock);
    // map |block| element to index into |survivors|
  @< Fill |survivors| with elements from |block| that survive at |gamma|,
     and for each, put into its slot in |loc| the index at which |survivors|
     contains it @>
@)
  @< Push list of parameters corresponding to |survivors| in |block|,... @>
  push_value(std::make_shared<int_value>@|
    (loc[start]==UndefBlock ? -1 : loc[start]));
@)
  using Pol = polynomials::Polynomial<int>;
  std::vector<Pol> pool { Pol(), Pol(1) };
@/@< Enumerate distinct Kazhdan-Lusztig polynomials from |kl_tab| into |pool|,
     and push lower unitriangular matrix with at position $(x,y)$ the index of
     the polynomial $P_{x',y'}$ returned by |kl_tab.KL_pol| for |dual_block|
     elements $x',y'$ corresponding respectively to $x,y$ @>
@)
  @< Transfer the coefficient vectors of the polynomials from |pool| to an array,
     and push that array @>
@)
  if (l==expression_base::single_value)
    wrap_tuple<4>();
}

@ Whether an element survives as |gamma| is determined using methods |singular|
and |survives| from |blocks::common_block|.

@< Fill |survivors| with elements from |block| that survive at |gamma|,
     and for each, put into its slot in |loc| the index at which |survivors|
     contains it @>=
{
  const RankFlags singular = block.singular(gamma);
  for (BlockElt z=0; z<block.size(); ++z)
    if (block.survives(z,singular))
    @/{@;
      loc[z] = survivors.size();
      survivors.push_back(z);
    }
}

@ The |int_Matrix| constructor with a single |int| argument produces an identity
matrix of the specified size. We proceed to fill just the part strictly below
the diagonal from |kl_tab|. Since that table holds polynomials with |unsigned|
coefficients, we need to convert them (using the |Polynomial<int>| constructor)
to a signed type, so that the vectors can then be moved into the row of vectors
returned to the user. Calling |hash_match| for the converted polynomial both
ensures that it is copied to |pool| if not yet present, and replaces it by its
index into |pool|.

@< Enumerate distinct Kazhdan-Lusztig polynomials from |kl_tab| into |pool|,
   and push lower unitriangular matrix with at position $(x,y)$ the index of
   the polynomial $P_{x',y'}$ returned by |kl_tab.KL_pol| for |dual_block|
   elements $x',y'$ corresponding respectively to $x,y$ @>=
{
  const auto n_survivors = survivors.size();
  HashTable<IntPolEntry,unsigned int> hash(pool);
  own_matrix M_ind = std::make_shared<matrix_value>(int_Matrix(n_survivors));
  const BlockElt last=block.size()-1;
    // mapping |block->dual_block| is subtraction from |last|
  for (auto jt = survivors.begin(); not survivors.at_end(jt); ++jt)
  { BlockElt y = *jt; // index into |block|
    for (auto it = jt; not survivors.at_end(it); ++it)
    { BlockElt x = *it; // index into |block|
       M_ind->val(loc[x],loc[y]) =
         hash.match(Pol(kl_tab.KL_pol(last-x,last-y)));
    }
  }
  push_value(std::move(M_ind));
}

@ Rather than exporting the detailed KL data, the following functions compute
the $W$-graph respectively $W$-cells from the block of the parameter, and
export that.

@< Local function def...@>=
void param_W_graph_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt start; // will hold index in the block of the initial element
  auto& block = p->rt().lookup_full_block(p->val,start);
  push_value(std::make_shared<int_value>(start));
@)
  const kl::KL_table& kl_tab = block.kl_tab(nullptr);
   // this does the actual KL computation
  wgraph::WGraph wg = kl::wGraph(kl_tab);
@)
  own_row vertices=std::make_shared<row_value>(0);
  @< Push to |vertices| a list of pairs for each element of |wg|, each
     consisting of a descent set and a list of outgoing labelled edges @>
  push_value(std::move(vertices));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}
@)
void param_W_cells_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt start; // will hold index in the block of the initial element
  auto& block = p->rt().lookup_full_block(p->val,start);
  push_value(std::make_shared<int_value>(start));
@)
  const kl::KL_table& kl_tab = block.kl_tab(nullptr);
   // this does the actual KL computation
  wgraph::WGraph wg = kl::wGraph(kl_tab);
  wgraph::DecomposedWGraph dg(wg);
@)
  own_row cells=std::make_shared<row_value>(0);
  cells->val.reserve(dg.cellCount());
  for (unsigned int c = 0; c < dg.cellCount(); ++c)
  { auto& wg=dg.cell(c); // local $W$-graph of cell
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
  push_value(std::move(cells));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ The following code was isolated so that it can be reused below.

@< Push to |vertices| a list of pairs for each element of |wg|, each
   consisting of a descent set and a list of outgoing labelled edges @>=
vertices->val.reserve(wg.size());
for (unsigned int i = 0; i < wg.size(); ++i)
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

@ The following wrapper function makes available the Atlas implementation of
Tarjan's strong component algorithm, the one that is used to isolate cells in a
$W$-graph. It actually involves no atlas-specific types at all, as it encodes
its argument graph by a list of lists of integers (of which list~$i$ enumerates
the vertex numbers reachable by an outgoing edge from vertex~$i$), the first
output value encodes the strong components found as lists of vertices, and the
second value the induced graph on strong components, using the same structure as
for the argument graph.

@< Local function def...@>=
void strong_components_wrapper(expression_base::level l)
{
  shared_row graph = get<row_value>();
  const auto size = graph->val.size();
  OrientedGraph G (size);
  for (unsigned i=0; i<size; ++i)
  {
    const row_value* p = force<row_value>(graph->val[i].get());
  @/ auto& edges = G.edgeList(i);
    edges.reserve(p->val.size());
    for (size_t i=0; i<p->val.size(); ++i)
      edges.push_back(force<int_value>(p->val[i].get())->uint_val());
    for (unsigned v : edges)
      if (v>=size)
      { std::ostringstream o;
        o << "Edge target " << v @|
          << " out of bounds (should be <" << size << ")";
        throw runtime_error(o.str());
      }
  }
  if (l==expression_base::no_value)
    return;
@)
  std::unique_ptr<OrientedGraph> induced(new OrientedGraph);
  const auto pi = G.cells(induced.get()); // invoke Tarjan's algorithm
@)
  { std::vector<std::shared_ptr<row_value> > part(pi.classCount());
    auto sizes = pi.class_sizes();
    for (unsigned i=0; i<part.size(); ++i)
    @/{@;
      part[i]=std::make_shared<row_value>(0);
      part[i]->val.reserve(sizes[i]);
    }
    for (unsigned long n=0; n<size; ++n)
      part[pi.class_of(n)]->val.push_back(std::make_shared<int_value>(n));
    own_row partition_list = std::make_shared<row_value>(pi.classCount());
    for (unsigned i=0; i<part.size(); ++i)
      partition_list->val[i] = std::move(part[i]); // widen shared pointer
    push_value(std::move(partition_list));
  }
@)
  {
    own_row edge_list_list = std::make_shared<row_value>(induced->size());
    for (unsigned i=0; i<induced->size(); ++i)
    {
      const auto& edges = induced->edgeList(i);
      auto dest = std::make_shared<row_value>(edges.size());
      for (unsigned j=0; j<edges.size(); ++j)
        dest->val[j] = std::make_shared<int_value>(edges[j]);
      edge_list_list->val[i] = dest; // widen shared pointer

    }
    push_value(std::move(edge_list_list));
  }
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}


@ The function |extended_block| makes computation of extended blocks available
directly in \.{atlas}.

@< Local function def...@>=
void extended_block_wrapper(expression_base::level l)
{ auto delta =get<matrix_value>();
  shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  test_compatible(p->rc().inner_class(),delta);
  if (not ((delta->val-1)*p->val.gamma().numerator()).isZero())
    throw runtime_error@|("Involution does not fix infinitesimal character");
  if (l==expression_base::no_value)
    return;
@)
  const auto& rc = p->rc();
  BlockElt start;
  auto zm = repr::StandardReprMod::mod_reduce(rc,p->val);
  common_context ctxt(rc,zm.gamma_lambda());
  blocks::common_block block(ctxt,zm,start); // build full block
  @< Construct the extended block, then the return value components,
     calling |push_value| for each of them @>
@)
  if (l==expression_base::single_value)
    wrap_tuple<4>();
}

@ We somewhat laboriously convert internal information from the extended block
into a list of parameters and three tables in the form of matrices.

@< Construct the extended block... @>=
{ auto eb = block.extended_block(rc.inner_class().distinguished());
  own_row params = std::make_shared<row_value>(eb.size());
  int_Matrix types(eb.size(),eb.rank());
@/int_Matrix links0(eb.size(),eb.rank());
  int_Matrix links1(eb.size(),eb.rank());

  const auto& gamma=p->val.gamma();
  const RatWeight gamma_rho = gamma-rho(block.root_datum());
  for (BlockElt n=0; n<eb.size(); ++n)
  { auto z = eb.z(n); // number of ordinary parameter in |block|
    const Weight lambda_rho=gamma_rho.integer_diff<int>(block.gamma_lambda(z));
    params->val[n] = std::make_shared<module_parameter_value> @|
      (p->rf,rc.sr_gamma(block.x(z),lambda_rho,gamma));
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
KL polynomials, and then rewrites elements that have singular descents in terms
of those that have not (the ``survivors''), reduces the matrix to be indexed by
those elements only, and finally negates entries at positions with odd length
difference for the block elements corresponding to row and column. All this work
is actually performed inside call to |ext_kl::ext_KL_matrix|.

The function returns the extended block as list of parameters, a matrix, and
finally a list of vectors, to be interpreted as polynomials, and into which list
the matrix entries are indices.

@< Local function def...@>=
void extended_KL_block_wrapper(expression_base::level l)
{ auto delta = get<matrix_value>();
  auto p = get<module_parameter_value>();
  test_standard(*p,"Cannot generate extended block");
  test_compatible(p->rc().inner_class(),delta);
  if (l==expression_base::no_value)
    return;
@)
  std::vector<StandardRepr> block;
  own_matrix P_mat = std::make_shared<matrix_value>(int_Matrix());
  std::vector<ext_kl::Pol> pool;
  ext_kl::ext_KL_matrix(p->val,delta->val,p->rc(),block,P_mat->val,pool);
@)
  own_row param_list = std::make_shared<row_value>(block.size());
  for (BlockElt z=0; z<block.size(); ++z)
    param_list->val[z]=std::make_shared<module_parameter_value>
      (p->rf,std::move(block[z]));
  push_value(std::move(param_list));
  push_value(std::move(P_mat));
  @< Transfer the coefficient vectors of the polynomials from |pool|... @>

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
install_function(parameter_height_wrapper,@|"height" ,"(Param->int)");
install_function(param_to_K_type_wrapper,@|"K_type", "(Param->KType)");
install_function(K_type_to_param_wrapper,@|"param", "(KType->Param)");
install_function(parameter_eq_wrapper,@|"=", "(Param,Param->bool)");
install_function(parameter_neq_wrapper,@|"!=", "(Param,Param->bool)");
install_function(parameter_equivalent_wrapper,@|"equivalent"
                ,"(Param,Param->bool)");
install_function(is_standard_wrapper,@|"is_standard" ,"(Param->bool)");
install_function(is_dominant_wrapper,@|"is_dominant" ,"(Param->bool)");
install_function(is_zero_wrapper,@|"is_zero" ,"(Param->bool)");
install_function(is_semifinal_wrapper,@|"is_semifinal" ,"(Param->bool)");
install_function(is_final_wrapper,@|"is_final" ,"(Param->bool)");
install_function(parameter_dominant_wrapper,@|"dominant" ,"(Param->Param)");
install_function(parameter_normal_wrapper,@|"normal" ,"(Param->Param)");
install_function(parameter_cross_wrapper,@|"cross" ,"(int,Param->Param)");
install_function(parameter_Cayley_wrapper,@|"Cayley" ,"(int,Param->Param)");
install_function(root_parameter_cross_wrapper,@|"cross" ,"(vec,Param->Param)");
install_function(root_parameter_Cayley_wrapper,@|"Cayley" ,"(vec,Param->Param)");
install_function(parameter_twist_wrapper,@|"twist" ,"(Param->Param)");
install_function(parameter_outer_twist_wrapper,@|"twist" ,"(Param,mat->Param)");
install_function(orientation_number_wrapper,@|"orientation_nr" ,"(Param->int)");
install_function(reducibility_points_wrapper,@|
		"reducibility_points" ,"(Param->[rat])");
install_function(scale_parameter_wrapper,"*", "(Param,rat->Param)");
@)
install_function(print_c_block_wrapper,@|"print_block","(Param->)");
install_function(print_pc_block_wrapper,@|"print_partial_block","(Param->)");
install_function(common_block_wrapper,@|"block" ,"(Param->[Param],int)");
install_function(partial_common_block_wrapper,@|"partial_block"
                ,"(Param->[Param])");
install_function(param_length_wrapper,@|"length","(Param->int)");
install_function(block_Hasse_wrapper,@|"block_Hasse","(Param->mat)");
install_function(KL_block_wrapper,@|"KL_block"
                ,"(Param->[Param],int,mat,[vec])");
install_function(dual_KL_block_wrapper,@|"dual_KL_block"
                ,"(Param->[Param],int,mat,[vec])");
install_function(partial_KL_block_wrapper,@|"partial_KL_block"
                ,"(Param->[Param],mat,[vec])");
install_function(param_W_graph_wrapper,@|"W_graph"
		,"(Param->int,[[int],[int,int]])");
install_function(param_W_cells_wrapper,@|"W_cells"
                ,"(Param->int,[[int],[[int],[int,int]]])");
install_function(strong_components_wrapper,@|"strong_components"
                ,"([[int]]->[[int]],[[int]])");
install_function(extended_block_wrapper,@|"extended_block"
                ,"(Param,mat->[Param],mat,mat,mat)");
install_function(extended_KL_block_wrapper,@|"partial_extended_KL_block"
                ,"(Param,mat->[Param],mat,[vec])");

@*1 Polynomials formed from parameters.
%
When working with parameters for standard modules, and notably with the
deformation formulas, the need arises to keep track of formal sums of standard
modules, or their irreducible quotients (the interpretation will depend on the
context) with split integer coefficients. We call these formal sums virtual
modules, and in the \.{atlas} language designate the corresponding basic type
as \.{ParamPol}.

@*2 Class definition.
%
The library provides a type |SR_poly| in which such sums can be
efficiently maintained. In order to use it we must have seen the header file
for the module \.{free\_abelian} on which the implementation is based. While
that class itself does not have such an invariant, the handling of these
formal sums in \.{atlas} will be such that all terms are ensured to have the
predicate |is_final| true, which ensures a number of desirable properties,
including having a dominant representative $\gamma$ of the infinitesimal
character. Only under such restriction can it be guaranteed that equivalent
terms (which now must actually be equal) will always be combined, and the test
for the sum being zero therefore mathematically correct.

@~Like for KGB elements, we maintain a shared pointer to the real form value, so
that it will be assured to survive as long as parameters for it exist.

@< Type definitions @>=
struct virtual_module_value : public value_base
{ shared_real_form rf;
  SR_poly val;
@)
  virtual_module_value(const shared_real_form& form, const SR_poly& v)
  : rf(form), val(v) @+{}
  virtual_module_value(const shared_real_form& form, SR_poly&& v)
  : rf(form), val(std::move(v)) @+{}
  ~virtual_module_value() @+{}
@)
  virtual void print(std::ostream& out) const;
  static const char* name() @+{@; return "module parameter"; }
  virtual_module_value (const virtual_module_value& v) = default;
    // we use |get_own<virtual_module_value>|
@)
  const Rep_context& rc() const @+{@; return rf->rc(); }
  Rep_table& rt() const @+{@; return rf->rt(); }
  void assign_coef(const module_parameter_value& t, const Split_integer& c);
};
@)
typedef std::unique_ptr<virtual_module_value> virtual_module_ptr;
typedef std::shared_ptr<const virtual_module_value> shared_virtual_module;
typedef std::shared_ptr<virtual_module_value> own_virtual_module;

@ Printing a virtual module value calls the free function |repr::print_SR_poly|
to do the actual work. It traverses the |std::map| that is hidden in the
|Free_Abelian| class template, and prints individual terms by printing the
|Split_integer| coefficient, followed by the parameter through a call of
|print_stdrep|. When either all coefficients are integers or all coefficients
are (integer) multiples of~$s$, it suppresses the component that is always~$0$;
this is particularly useful if polynomials are used to encode $\Zee$-linear
combinations of parameters.

@< Function def...@>=
void virtual_module_value::print(std::ostream& out) const
{@; print_SR_poly(out,val,rc()); }

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
{ shared_real_form rf = get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<virtual_module_value> @| (rf,SR_poly()));
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


@ We provide a conversion from a module parameter to a $K$-type (restricting to
$K$, and therefore forgetting the $\nu$ component of the infinitesimal
character). We also allow implicit conversion of a module parameter to a virtual
module, which involves expansion by the |Rep_context::expand_final| method
to \emph{final} parameters (there can be zero, one, or more of them, and they
can have positive or negative integer coefficients), to initiate the invariant
that only (dominant, standard, nonzero) final parameters can be stored in a
|virtual_module_value|. Then there is a conversion from module parameters
directly to $K$-type polynomials, the conversion doing first the restriction to
$K$ and then the expansion into final $K$-types. Finally, the conversion from
module parameters to $K$-types can be extended linearly to a map from virtual
modules to $K$-type polynomials, still defined mathematically by restriction
to$~K$.

We don't want the restriction maps to be implicit conversions. For the
restriction on the level of polynomials, such conversion would make it
impossible to define certain operations with the same name both for
types \.{KTypePol} and \.{ParamPol}~: for instance we have (built-in) instances
of the operator \.* with arguments types \.{(int,KTypePol)}
and \.{(Split,ParamPol)}, and in the presence of an implicit conversion
from \.{ParamPol} to \.{KTypePol} this combination would be forbidden, due to
possible ambiguity for \.* with argument types \.{int} and \.{ParamPol}. An
implicit conversion from module parameters to $K$-types would be problematic for
similar reasons: it would create divergent implicit conversions from
the \.{Param} type, leading to ambiguity when and argument could match two
different overloads using one or the other of these conversions. Although
currently the system does not test for such divergences, its logic assumes that
when there are implicit conversions from a same type to two different types,
then there is also at least one conversion between those two types so that
there is either a hierarchy or an equivalence of inter-convertible types.

@< Local function def...@>=
void param_to_poly()
{ shared_module_parameter p = get<module_parameter_value>();
  const auto& rf=p->rf;
  push_value(std::make_shared<virtual_module_value> @|
    (rf,rf->rc().expand_final(p->val)));
}
@)
void param_poly_to_K_type_poly_wrapper(expression_base::level l)
{ shared_virtual_module p = get<virtual_module_value>();
  if (l==expression_base::no_value)
    return;
  const auto& rf=p->rf;
  K_repr::K_type_pol result;
  for (const auto& term : p->val)
  { auto finals = rf->rc().finals_for(rf->rc().sr_K(term.first));
    for (auto it = finals.begin(); not finals.at_end(it); ++it)
       result.add_term(std::move(it->first),term.second*it->second);
  }
  push_value(std::make_shared<K_type_pol_value>(p->rf,std::move(result)));
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
  if (m->rf!=p->rf and m->rf->val!=p->rf->val)
    // test like |real_form_new_wrapper| does
    throw runtime_error @|
      ("Real form mismatch when subscripting ParamPol value");
  test_standard(*p,"In subscription of ParamPol value");
     // it is OK to do this test before |make_dominant|
  StandardRepr sr = p->val;
  p->rc().make_dominant(sr);
  if (l!=expression_base::no_value)
    push_value(std::make_shared<split_int_value>(m->val[sr]));
}

@ For modifying an individual coefficient in a |virtual_module_value| (possibly
adding or removing a term in the process), a method |assign_coef| is provided.
It will be called from the \.{axis} compilation unit (the programming language
interpreter). The implementation is similar to that of |Free_Abelian::add_term|,
using the |std:map| interface from which |Free_Abelian| was derived; as is the
case there a term can get created, modified, or deleted, of nothing can happen
at all,

@< Function def...@>=
void virtual_module_value::assign_coef
  (const module_parameter_value& t, const Split_integer& c)
{
  test_final(t,"In coefficient assignment for ParamPol value");
  auto interval = val.equal_range(t.val);
  if (interval.first==interval.second) // no term present
  { if (not c.is_zero())
    val.insert(interval.first,std::make_pair(t.val,c));
  }
  else if (c.is_zero())
    val.erase(interval.first);
  else interval.first->second = c;
}

@ The main operations for virtual modules are addition and subtraction of
parameters to or from them.

@< Local function def...@>=
void add_module_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  own_virtual_module accumulator = get_own<virtual_module_value>();
  if (accumulator->rf!=p->rf)
    throw runtime_error @|
      ("Real form mismatch when adding a Param to a ParamPol");
  if (l!=expression_base::no_value)
  @/{@; accumulator->val+= p->rc().expand_final(p->val);
    push_value(std::move(accumulator));
  }
}

void subtract_module_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  own_virtual_module accumulator = get_own<virtual_module_value>();
  if (accumulator->rf!=p->rf)
    throw runtime_error @|
      ("Real form mismatch when subtracting a Param from a ParamPol");
  if (l!=expression_base::no_value)
  @/{@; accumulator->val-= p->rc().expand_final(p->val);
    push_value(std::move(accumulator));
  }
}

@ More generally than adding or subtracting, we can incorporate a term with
specified coefficient. Here, rather than building a polynomial with the set of
parameters but without the coefficient, as |expand_final| gives us, we directly
iterate over the list produced by |finals_for|, attaching |coef| for each of
its final parameters.

@< Local function def...@>=

void add_module_term_wrapper(expression_base::level l)
{ push_tuple_components(); // second argument is a pair |(coef,p)|
  shared_module_parameter p = get<module_parameter_value>();
  Split_integer coef=get<split_int_value>()->val;
  own_virtual_module accumulator = get_own<virtual_module_value>();

if (accumulator->rf!=p->rf)
    throw runtime_error@|("Real form mismatch when adding a term to a module");
  if (l==expression_base::no_value)
    return;
@)
  auto finals = p->rc().finals_for(p->val);
  for (auto it=finals.wcbegin(); not finals.at_end(it); ++it)
    accumulator->val.add_term(it->first,coef*it->second);
  push_value(std::move(accumulator));
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
    if (accumulator->rf!=p->rf)
      throw runtime_error@|("Real form mismatch when adding terms to a module");
     auto finals = p->rc().finals_for(p->val);
     for (auto it=finals.wcbegin(); not finals.at_end(it); ++it)
       accumulator->val.add_term(it->first,coef*it->second);
   }
  push_value(std::move(accumulator));
}

@ Naturally we also want to define addition and subtraction of two virtual
modules.

@< Local function... @>=
void add_virtual_modules_wrapper(expression_base::level l)
{
  own_virtual_module accumulator = get_own<virtual_module_value>();
  shared_virtual_module addend = get<virtual_module_value>();
  if (accumulator->rf!=addend->rf)
    throw runtime_error @|("Real form mismatch when adding two modules");
  if (l!=expression_base::no_value)
  @/{@; accumulator->val += addend->val;
    push_value(std::move(accumulator));
  }
}
@)
void subtract_virtual_modules_wrapper(expression_base::level l)
{
  shared_virtual_module subtrahend = get<virtual_module_value>();
  own_virtual_module accumulator = get_own<virtual_module_value>();
  if (accumulator->rf!=subtrahend->rf)
    throw runtime_error@|("Real form mismatch when subtracting two modules");
  if (l!=expression_base::no_value)
  @/{@; accumulator->val -= subtrahend->val;
    push_value(std::move(accumulator));
  }
}

@ Scalar multiplication potentially makes coefficients zero, in which case the
corresponding terms need to be removed to preserve the invariant that no zero
terms are stored in a virtual module. For integer multiplication we just need to
check for multiplication by $0$, and produce an empty module when this happens.
Because the integer is not on the stack top, this requires a somewhat unusual
manoeuvre. The case of multiplication by zero needs to be handled separately,
since we cannot allow introducing terms with zero coefficients. It could have
been handled more easily though, by testing the factor~|c| just before the |for|
loop, and performing |m->erase()| instead if |c==0|; this is what we used to do.
However that might involve duplicating the virtual module and then erasing the
copy, which is inefficient, and now avoided. This might seem a rare case, but it
is not really: often functions handling a \.{ParamPol} argument $P$ need to
start with an empty module for the same real form; writing $0*P$ is quite a
convenient way to achieve this.

@< Local function... @>=

void int_mult_virtual_module_wrapper(expression_base::level l)
{ int c =
    force<int_value>(execution_stack[execution_stack.size()-2].get())
    // value below top
    ->int_val();
  if (c==0) // then do multiply by $0$ efficiently:
  { shared_virtual_module m = get<virtual_module_value>();
      // |m| is needed for |m->rc()|
    pop_value();
    if (l!=expression_base::no_value)
    @/push_value@|(std::make_shared<virtual_module_value>
        (m->rf,SR_poly()));
  }
  else
  { own_virtual_module m = get_own<virtual_module_value>();
     // will modify our copy now
    pop_value();
    if (l!=expression_base::no_value)
    { for (auto& term : m->val)
        term.second *= c;
      push_value(std::move(m));
    }
  }
}

@ Matters are similar but somewhat subtler for scalar multiplication by split
integers, because these have zero divisors. We make a $4$-way branch depending
on whether the split integer multiplicand has zero evaluation are $s=1$ and/or
at $s=-1$: if either of these is the case we have a zero divisor, and may expect
some terms to be killed by the multiplication, and if both hold we are
multiplying by $0$. In presence of a nonzero zero divisor, we do not call
|get_own| to true to modify the argument in place, as it seems more efficient to
just reconstruct a potentially small product separately. As in the previous
function we pick up the multiplicand (which is the first argument) when it is
not at the top of the stack; later on, we do not pop it off the stack either,
but rather overwrite it by move-assigning the returned |virtual_module_value| to
the stack top location |execution_stack.back()|. Since no tests are performed
here, the |no_value| case can be handled right at the beginning, by simply
dropping the arguments from the stack when it applies.

@< Local function... @>=

void split_mult_virtual_module_wrapper(expression_base::level l)
{  if (l==expression_base::no_value)
   {@; pop_value();
       pop_value();
       return;
   }
   auto c =
    force<split_int_value>(execution_stack[execution_stack.size()-2].get())
    // below top
    ->val;
  if (c.s_to_1()==0)
  // then coefficient is multiple of $1-s$; don't try to modify module in place
  { shared_virtual_module m = get<virtual_module_value>();
    SR_poly result;
    if (c.s_to_minus_1()!=0) // otherwise coefficient is $0$ and keep null module
      for (const auto& term: m->val)
        if (term.second.s_to_minus_1()!=0)
          result.add_term(term.first,term.second*c);
    execution_stack.back()=  // replace stack top
        std::make_shared<virtual_module_value>(m->rf,result);
    return;
  }
@)
  else if (c.s_to_minus_1()==0)
  // then coefficient is multiple of $1+s$; don't modify in place
  { shared_virtual_module m = get<virtual_module_value>();
    SR_poly result;
    for (const auto& term: m->val)
      if (term.second.s_to_1()!=0)
        result.add_term(term.first,term.second*c);
    execution_stack.back()= // replace stack top
        std::make_shared<virtual_module_value>(m->rf,result);
    return;
  }
  // now we are multiplying by a non zero-divisor
  own_virtual_module m = get_own<virtual_module_value>();
     // will modify our copy now
  for (auto& term : m->val)
        term.second *= c; // multiply coefficients by split integer
  execution_stack.back()=std::move(m);  // replace stack top
}

@ For a nonzero virtual module, it is useful to be able to select some term that
is present, without looping over all terms. The most useful choice it the final
term, be we allow taking the first term as well.

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
  push_value(std::make_shared<module_parameter_value>
	(m->rf,std::move(term.first)));
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
  push_value(std::make_shared<module_parameter_value>
  	(m->rf,std::move(term.first)));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ Sometimes we want to ignore parameters whose height exceeds a given limit.
Since are sorted by height this can be done fairly efficiently, and the built-in
|truncate_above_height| will do this.

@h <algorithm>
@< Local function... @>=

void truncate_param_poly_above_wrapper (expression_base::level l)
{ int arg = get<int_value>()->int_val();
  shared_virtual_module P = get<virtual_module_value>();
  if (l==expression_base::no_value)
    return;
@)
  auto point = P->val.begin(); // default to no terms
  if (arg>=0)
  {
    unsigned bound = arg; // use |unsigned| value for height comparison
    auto predicate = @/ @[ [bound](const SR_poly::value_type& term)
      {@; return term.first.height()<=bound; } @];
    point = std::partition_point(point,P->val.end(),predicate);
  }
  if (point== P->val.end())
    push_value(P); // if all terms are selected, just return |P|
  else
    push_value(std::make_shared<virtual_module_value>@|
      (P->rf,SR_poly(P->val.begin(),point)));
}

@ Here is a variation of the scaling function for parameters that operates on
entire virtual module.

@< Local function def...@>=
void scale_poly_wrapper(expression_base::level l)
{ shared_rat f = get<rat_value>();
  shared_virtual_module P = get<virtual_module_value>();
  if (l!=expression_base::no_value)
    push_value@|(std::make_shared<virtual_module_value>
      (P->rf,P->rc().scale(P->val,f->rat_val())));
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
effect. On the other hand, since the |deformation_terms| method assumes its
arguments to be final, we apply |finals_for| to |p->val| and sum over any
(final) parameters this might produce.

@< Local function def...@>=
void deform_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  if (l==expression_base::no_value)
    return;
@)
  auto& rt=p->rt();
  const auto& rc = p->rc();
@/SR_poly result;
  auto finals = rc.finals_for(p->val);
  for (auto it=finals.begin(); not finals.at_end(it); ++it)
  {
    auto& q = it->first;
    BlockElt q_index; // will hold index of |q| in the block
    auto& block = rt.lookup(q,q_index); // generate partial common block
    RatWeight diff = rc.offset(q,block.representative(q_index));
    for (auto&& term : rt.deformation_terms(block,q_index,diff,q.gamma()))
    result.add_term(std::move(term.first),
                    Split_integer(term.second,-term.second)*it->second);
  }

  push_value(std::make_shared<virtual_module_value>(p->rf,std::move(result)));
}

@ There is also a variation |twisted_deform| that uses twisted KLV polynomials
instead, for the distinguished involution $\delta$ of the inner class. For the
code here the difference consists mainly of calling the
|Rep_table::twisted_deformation_terms| method instead of
|Rep_table::deformation_terms|. However, that method requires a $\delta$-fixed
involution, so we need to test for that here. If the test fails we report an
error rather than returning for instance a null module, since a twisted
deformation formula for a non-fixed parameter makes little sense; the user
should avoid asking for it. Similarly the twisted variant cannot allow non
dominant parameters, as this would internally produce an |SR_poly| value with
non-dominant terms, which should never happen.

@< Local function def...@>=
void twisted_deform_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  auto& rt=p->rt();
  const auto& delta=rt.inner_class().distinguished();
  test_standard(*p,"Cannot compute twisted deformation terms");
  if (not rt.is_twist_fixed(p->val,delta))
    throw runtime_error@|("Parameter not fixed by inner class involution");
  test_final(*p,"Twisted deformation requires final parameter");
  if (l==expression_base::no_value)
    return;
@)
  BlockElt entry_elem;
  auto& block = rt.lookup(p->val,entry_elem);
    // though by reference, does not change |p->val|
  RatWeight diff = rt.offset(p->val, block.representative(entry_elem));
  block.shift(diff);
  auto& eblock = block.extended_block(rt.shared_poly_table());
  block.shift(-diff);
@)
  RankFlags singular = block.singular(p->val.gamma());
  RankFlags singular_orbits;
  for (weyl::Generator s=0; s<eblock.rank(); ++s)
    singular_orbits.set(s,singular[eblock.orbit(s).s0]);
@)
  auto terms = rt.twisted_deformation_terms@|(block,eblock,entry_elem,
					     singular_orbits,
                                             diff,p->val.gamma());
  SR_poly result;
  for (auto&& term : terms)
    result.add_term(std::move(term.first),
                    Split_integer(term.second,-term.second));

  push_value(std::make_shared<virtual_module_value>(p->rf,std::move(result)));
}

@ An intermediate between a single deformation at a parameter and a full
recursive deformation to tempered parameters of all terms produced by such a
deformation (plus the original parameter itself), here is a function that
continues deformations within a given block until no more terms are produced,
and then slides down all the contributions from this block, either to the next
reducibility point or all the way to a tempered parameter. Due to the way this
is set up, it is suited for a truncated computation where all terms above a
given height bound are ignored, so we provide such a bound argument. Also, since
we are dealing with a whole block, it will be more efficient if we collectively
treat a set of terms at hand that all live in a same block. This is achieved by
providing a polynomial |accumulator| and an initial parameter (presumably
corresponding to a term in |accumulator|); all terms in the block of |p| are
extracted from |accumulator| and deformed to give the first component of the
result, with the second component being the remainder of |accumulator|.
Returning two parts can be helpful in understanding the details of the
deformation, but in practice the deformed terms are probably to be added back to
the accumulator after which another block is deformed.

@s SR_poly vector

@< Local function def...@>=
void block_deform_wrapper(expression_base::level l)
{ int bound = get<int_value>()->int_val();
  own_virtual_module accumulator = get_own<virtual_module_value>();
  own_module_parameter p = get_own<module_parameter_value>();
  if (l==expression_base::no_value)
    return;
@)
  SR_poly result;
  if (not p->rc().nu(p->val).is_zero())
  {
    auto deformed = p->rt().block_deformation_to_height @|
      (p->val,accumulator->val
      ,bound>=0 ? static_cast<repr::level>(bound) : repr::level(-1));
    for (const auto& term : deformed)
    { auto rps = p->rc().reducibility_points(term.first);
      auto i =
        rps.size()>0 and rps.back()==RatNum(1) ? rps.size()-1 : rps.size();
      RatNum f = i>0 ? rps[i-1] : RatNum(0);
      result.add_term(p->rc().scale(term.first,f),term.second);
    }
  }
  push_value(std::make_shared<virtual_module_value>(p->rf,std::move(result)));
  push_value(accumulator);
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ Here is a recursive form of the deformation, which stores intermediate
results for efficiency in the |Rep_table| structure |p->rt()| that is stored
within the |real_form_value|.

@s K_type_poly vector

@< Local function def...@>=
void full_deform_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  if (l==expression_base::no_value)
    return;
@)
  repr::K_type_poly result;
    // this is the data type used by |Rep_table::deformation|
  auto finals = p->rc().finals_for(p->val);
  for (auto it=finals.begin(); not finals.at_end(it); ++it)
    for (auto&& term : p->rt().deformation(std::move(it->first)))
      result.add_term(std::move(term.first),term.second*it->second);
@) // now convert from (tabled) |repr::K_type_poly| to |K_repr::K_type_pol|
  push_value(std::make_shared<K_type_pol_value>@|
    (p->rf,export_K_type_pol(p->rt(),result)));
}
@)
void twisted_full_deform_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  const auto& rc=p->rc(); auto& rt=p->rt();
  test_standard(*p,"Cannot compute full twisted deformation");
  if (not rc.is_twist_fixed(p->val))
    throw runtime_error@|("Parameter not fixed by inner class involution");
  if (l==expression_base::no_value)
    return;
@)
  auto finals = @;ext_block::
    extended_finalise(rc,p->val,rc.inner_class().distinguished());
  repr::K_type_poly result;
    // this is the data type used by |Rep_table::deformation|
  for (auto it=finals.cbegin(); it!=finals.cend(); ++it)
  { bool flip;
    const auto& def = rt.twisted_deformation(std::move(it->first),flip);
    result.add_multiple(def,
	flip!=it->second ? Split_integer(0,1) : Split_integer(1,0));
  }
@) // now convert from (tabled) |repr::K_type_poly| to |K_repr::K_type_pol|
  push_value(std::make_shared<K_type_pol_value>@|
    (p->rf,export_K_type_pol(p->rt(),result)));
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
  test_final(*p,"Cannot compute Kazhdan-Lusztig sum");
  if (l!=expression_base::no_value)
    push_value(std::make_shared<virtual_module_value>@|
      (p->rf,p->rt().KL_column_at_s(p->val)));
}
@)
void twisted_KL_sum_at_s_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot compute Kazhdan-Lusztig sum");
  test_final(*p,"Cannot compute Kazhdan-Lusztig sum");
  auto sr=p->val; // take a copy
  p->rc().make_dominant(sr);
    // |is_twist_fixed| and |twisted_KL_column_at_s| like this
  if (not p->rc().is_twist_fixed(sr))
    throw runtime_error@|("Parameter not fixed by inner class involution");
  if (l!=expression_base::no_value)
    push_value (std::make_shared<virtual_module_value>@|
      (p->rf,p->rt().twisted_KL_column_at_s(sr)));
}

@ Here is a function to directly access a stored Kazhdan-Lusztig polynomial

@< Local function def...@>=
void KL_column_wrapper(expression_base::level l)
{ own_module_parameter p = get_own<module_parameter_value>();
  test_standard(*p,"Cannot compute Kazhdan-Lusztig column");
  test_final(*p,"Cannot compute Kazhdan-Lusztig column");
  if (l==expression_base::no_value)
    return;
@)
  auto col = p->rt().KL_column(p->val);
  BlockElt z;
  const blocks::common_block& block = p->rt().lookup(p->val,z);
  RatWeight diff = p->rc().offset(p->val, block.representative(z));
  own_row column = std::make_shared<row_value>(0);
  column->val.reserve(length(col));
  for (auto it=col.wcbegin(); not col.at_end(it); ++it)
  {
    StandardRepr sr = block.sr(it->first,diff,p->val.gamma());
    auto tup = std::make_shared<tuple_value>(3);
    tup->val[0] = std::make_shared<int_value>(it->first);
    tup->val[1] = std::make_shared<module_parameter_value>(p->rf,std::move(sr));
    tup->val[2] = std::make_shared<vector_value>@|(
      std::vector<int>(it->second.begin(),it->second.end()));
    column->val.push_back(std::move(tup));
  }
  push_value(std::move(column));
}
@)

@ We add another function in which the external involution is an argument

@< Local function def...@>=
void external_twisted_KL_sum_at_s_wrapper(expression_base::level l)
{ shared_matrix delta = get<matrix_value>();
  shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p,"Cannot compute Kazhdan-Lusztig sum");
  test_final(*p,"Cannot compute Kazhdan-Lusztig sum");
  test_compatible(p->rc().inner_class(),delta);
  if (not p->rc().is_twist_fixed(p->val,delta->val))
    throw runtime_error("Parameter not fixed by given involution");
  if (l!=expression_base::no_value)
    push_value (std::make_shared<virtual_module_value>@|
      (p->rf,twisted_KL_column_at_s(p->rc(),p->val,delta->val)));
}

@ The function |scale_extended| is intended for use with in the deformation
algorithm when interpreting parameters as specifying a representation of the
extended group. One can arrange that deformation starts with a final parameter,
for which $\gamma$ is dominant and which has no singular descent (simple)
reflections, but when deforming this condition may be lost. In the ordinary
deformation algorithm this is taken care of by an implicit conversion of the
parameter when it gets used to construct a block or when it is contributed to a
virtual module. However this may potentially cause the default choice of
extended representation associated to the parameter to flip. The function below
will scale the parameter, perform the conversion to a final parameter using
extended parameters, and return the resulting parameter plus an indication of
whether a flip occurred. The function |scaled_extended_finalise| in the
module \\{ext\_block} does the actual work.

@< Local function def...@>=

void scale_extended_wrapper(expression_base::level l)
{ auto factor = get<rat_value>();
  auto delta = get<matrix_value>();
  auto p = get<module_parameter_value>();
  const StandardRepr sr = p->val;
  const auto& rc = p->rc();
  test_final(*p,"Cannot scale extended parameter");
  if (not factor->val.is_positive())
    throw runtime_error("Factor in scale_extended must be positive");
  test_compatible(p->rc().inner_class(),delta);
  if (not rc.is_twist_fixed(sr,delta->val))
    throw runtime_error@|
      ("Parameter to be scaled not fixed by given involution");
  if (l==expression_base::no_value)
    return;
@)
  auto result = @;ext_block::scaled_extended_finalise
    (rc,sr,delta->val,factor->rat_val());
  push_value(std::make_shared<module_parameter_value>
    (p->rf,std::move(result.first)));
  push_value(whether(result.second));
  if (l==expression_base::single_value)
    wrap_tuple<2>();
}

@ The function |K_type_pol_extended| is like |scale_extended|, but is used for
scaling by a factor~$0$. As a consequence the result can be represented in the
format of $K$-types rather than parameters, and there need not a single $K$-type
produced from a given parameter, due to the possibility of singular real and
imaginary reflections being applied (the former might produce one or two
contributions from a Hecht-Schmid identity, the latter might remove a
contribution). Since we handle failure of dominance and of finality for
$K$-types anyway, we need only require the argument parameter to satisfy
|is_standard| here. The function |extended_restrict_to_K| in the
module \\{ext\_block} does the actual work.

@< Local function def...@>=

void K_type_pol_extended_wrapper(expression_base::level l)
{ auto delta = get<matrix_value>();
  auto p = get<module_parameter_value>();
  const auto& rc = p->rc();
  test_standard(*p,"Parameter in K_type_pol_extended| must be standard");
  test_compatible(rc.inner_class(),delta);
  if (not p->rc().is_twist_fixed(p->val,delta->val))
    throw runtime_error("Parameter not fixed by given involution");
  if (l==expression_base::no_value)
    return;
@)
  auto result = @;ext_block::extended_restrict_to_K(rc,p->val,delta->val);
  push_value (std::make_shared<K_type_pol_value>(p->rf,std::move(result)));
}


@ The function |finalize_extended| is useful in expanding a single module
parameter into a linear combination of such parameters. The terms on the list
are paired with a Boolean attribute recording a possible flip of extended
parameters accumulated when the term. This flip is recorded with as coefficient
the split integer unit~$s$, since it should be interpreted as a signature flip
(in the ordinary finalisation procedure, coefficients are purely integer, i.e.,
split integers without any $s$ component). The function |extended_finalise| in
the module \\{ext\_block} does the actual work.

@< Local function def...@>=

void finalize_extended_wrapper(expression_base::level l)
{ auto delta = get<matrix_value>();
  auto p = get<module_parameter_value>();
  const auto& rc = p->rc();
  test_standard(*p,"Cannot finalize extended parameter");
  test_compatible(rc.inner_class(),delta);
  if (not p->rc().is_twist_fixed(p->val,delta->val))
    throw runtime_error("Parameter not fixed by given involution");
  if (l==expression_base::no_value)
    return;
@)
  auto params = @;ext_block::extended_finalise(rc,p->val,delta->val);
  SR_poly result;
  for (auto it=params.begin(); it!=params.end(); ++it)
    result.add_term(it->first
                   ,it->second ? Split_integer(0,1) : Split_integer(1,0));
  push_value (std::make_shared<virtual_module_value>(p->rf,std::move(result)));
}


@ Finally we install everything related to polynomials formed from parameters.
@< Install wrapper functions @>=
install_function(virtual_module_wrapper,@|"null_module","(RealForm->ParamPol)");
install_function(real_form_of_virtual_module_wrapper,@|"real_form"
		,"(ParamPol->RealForm)");
install_function(virtual_module_size_wrapper,@|"#","(ParamPol->int)");
install_function(param_poly_to_K_type_poly_wrapper,@|"K_type_pol"
		,"(ParamPol->KTypePol)");
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
install_function(truncate_param_poly_above_wrapper,@|"truncate_above_height"
		,"(ParamPol,int->ParamPol)");
install_function(scale_poly_wrapper,"*", "(ParamPol,rat->ParamPol)");

install_function(deform_wrapper,@|"deform" ,"(Param->ParamPol)");
install_function(twisted_deform_wrapper,@|"twisted_deform" ,"(Param->ParamPol)");
install_function(block_deform_wrapper,@|"block_deform"
                ,"(Param,ParamPol,int->ParamPol,ParamPol)");
install_function(full_deform_wrapper,@|"full_deform","(Param->KTypePol)");
install_function(twisted_full_deform_wrapper,@|"twisted_full_deform"
                ,"(Param->KTypePol)");
install_function(KL_sum_at_s_wrapper,@|"KL_sum_at_s","(Param->ParamPol)");
install_function(twisted_KL_sum_at_s_wrapper,@|"twisted_KL_sum_at_s"
                ,"(Param->ParamPol)");
install_function(KL_column_wrapper,@|"KL_column","(Param->[int,Param,vec])");
install_function(external_twisted_KL_sum_at_s_wrapper,@|"twisted_KL_sum_at_s"
                ,"(Param,mat->ParamPol)");
install_function(scale_extended_wrapper,@|"scale_extended"
                ,"(Param,mat,rat->Param,bool)");
install_function(K_type_pol_extended_wrapper,@|"K_type_pol_extended"
                ,"(Param,mat->KTypePol)");
install_function(finalize_extended_wrapper,@|"finalize_extended"
                ,"(Param,mat->ParamPol)");


@*1 Kazhdan-Lusztig tables. We implement a simple function that gives raw
access to the table of Kazhdan-Lusztig polynomials.

@< Local function def...@>=
void raw_KL_wrapper (expression_base::level l)
{ shared_Block b = get<Block_value>();
  const Block& block = b->val;
  if (l==expression_base::no_value)
    return;
@)
  b->kl_tab.fill(); // this does the actual KL computation
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(b->kl_tab.size()));
  for (unsigned int y=1; y<b->kl_tab.size(); ++y)
    for (unsigned int x=0; x<y; ++x)
      M->val(x,y) = b->kl_tab.KL_pol_index(x,y);
@)
  own_row polys = std::make_shared<row_value>(0);
  const auto& store = b->kl_tab.pol_store();
  polys->val.reserve(store.size());
  for (auto it=store.begin(); it!=store.end(); ++it)
    polys->val.emplace_back(std::make_shared<vector_value> @|
       (std::vector<int>(it->begin(),it->end())));
@)
  std::vector<int> length_stops
    (block.size()==0 ? 2 :block.length(block.size()-1)+2);
  length_stops[0]=0;
  for (unsigned int i=1; i<length_stops.size(); ++i)
    length_stops[i]=block.length_first(i);
@)
  push_value(std::move(M));
  push_value(std::move(polys));
  push_value(std::make_shared<vector_value>(length_stops));
  if (l==expression_base::single_value)
    wrap_tuple<3>();
}

@ For testing, it is useful to also have the dual Kazhdan-Lusztig tables. In
this case we cannot of course use the field |b->kl_tab| to store the KL
polynomials, so we here us a local |kl::KL_table| variable.

@< Local function def...@>=
void raw_dual_KL_wrapper (expression_base::level l)
{ shared_Block b = get<Block_value>();
  const Block& block = b->val;
  Block dual_block = Block::build(b->dual_rf->val,b->rf->val);

  std::vector<BlockElt> dual=blocks::dual_map(block,dual_block);
  kl::KL_table kl_tab(dual_block); kl_tab.fill();
  if (l==expression_base::no_value)
    return;
@)
  own_matrix M = std::make_shared<matrix_value>(int_Matrix(kl_tab.size()));
  for (unsigned int y=1; y<kl_tab.size(); ++y)
    for (unsigned int x=0; x<y; ++x)
      M->val(x,y) = kl_tab.KL_pol_index(dual[y],dual[x]);
@)
  own_row polys = std::make_shared<row_value>(0);
  const auto& store = kl_tab.pol_store();
  polys->val.reserve(store.size());
  for (auto it=store.begin(); it!=store.end(); ++it)
    polys->val.emplace_back(std::make_shared<vector_value> @|
       (std::vector<int>(it->begin(),it->end())));
@)
  std::vector<int> length_stops
    (block.size()==0 ? 2 :block.length(block.size()-1)+2);
  length_stops[0]=0;
  for (unsigned int i=1; i<length_stops.size(); ++i)
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
  own_module_parameter p = get_own<module_parameter_value>();
  test_standard(*p,"Cannot generate block");
  test_compatible(p->rc().inner_class(),delta);
  if (l==expression_base::no_value)
    return;
@)
  const auto& rc = p->rc();
  rc.make_dominant(p->val);
  const auto gamma = p->val.gamma();
  const auto srm = repr::StandardReprMod::mod_reduce(rc,p->val);
  BlockElt start;
  common_context ctxt(rc,srm.gamma_lambda());
  blocks::common_block block(ctxt,srm,start); // build full block
  if (not((delta->val-1)*gamma.numerator()).isZero())
  { // block not globally stable, so return empty values;
    push_value(std::make_shared<matrix_value>(int_Matrix()));
    push_value(std::make_shared<row_value>(0));
    push_value(std::make_shared<vector_value>(int_Vector()));
  }
  else
  {
    @;ext_block::ext_block eb = block.extended_block(delta->val);
    std::vector<ext_kl::Pol> pool;
    ext_KL_hash_Table hash(pool,4);
    ext_kl::KL_table klt(eb,&hash); klt.fill_columns();
  @)
    own_matrix M = std::make_shared<matrix_value>(int_Matrix(klt.size()));
    for (unsigned int y=1; y<klt.size(); ++y)
      for (unsigned int x=0; x<y; ++x)
      @/{@; auto inx = klt.KL_pol_index(x,y);
        M->val(x,y) = inx.second ? -inx.first : inx.first;
      }
  @)
    std::vector<int> length_stops(block.length(block.size()-1)+2);
    length_stops[0]=0;
    for (unsigned int i=1; i<length_stops.size(); ++i)
      length_stops[i]=eb.length_first(i);
  @)
    push_value(std::move(M));
    @< Transfer the coefficient vectors of the polynomials from |pool|... @>
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
{ shared_Block b = get<Block_value>();
  if (l=expression_base::no_value)
    return;
@)
  b->kl_tab.fill(); // this does the actual KL computation
  wgraph::WGraph wg = kl::wGraph(b->kl_tab);
@)
  own_row vertices=std::make_shared<row_value>(0);
  @< Push to |vertices| a list of pairs for each element of |wg|, each
     consisting of a descent set and a list of outgoing labelled edges @>
  push_value(std::move(vertices));
}

@ This function computes |W_cells| for a block, as list of nested integer
structures.

@< Local function def...@>=
void W_cells_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  if (l=expression_base::no_value)
    return;
@)
  b->kl_tab.fill(); // this does the actual KL computation
  wgraph::WGraph wg = kl::wGraph(b->kl_tab);
  wgraph::DecomposedWGraph dg(wg);
@)

  own_row cells=std::make_shared<row_value>(0);
  cells->val.reserve(dg.cellCount());
  for (unsigned int c = 0; c < dg.cellCount(); ++c)
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
  push_value(std::move(cells));
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
  shared_real_form rf= get<real_form_value>();
@)
  if (rf->ic_ptr.get()!=cc->ic_ptr.get())
    throw runtime_error("Inner class mismatch between arguments");
@.Inner class mismatch...@>
  BitMap b(rf->val.innerClass().Cartan_set(rf->val.realForm()));
  if (not b.isMember(cc->number))
    throw runtime_error("Cartan class not defined for real form");
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
    (*output_stream,
     cc->ic_ptr->val,cc->ic_ptr->interface,cc->number);
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
a real form as argument. Since full KGB listings can be very long, we provide a
version of |print_KGB| that in addition takes a list of KGB elements, and limits
the output to those elements; the internal function |kgb_io::print| that is
(also) used to implement the full listing, already provides for such a limitation
through a final argument. Since that argument requires a |KGBEltList| (a vector
of |KGBElt|), we do the conversion here, which also ensures there are no
out-of-range values in the list.

@h "kgb.h"
@h "kgb_io.h"

@< Local function def...@>=
void print_KGB_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
@)
  *output_stream
    << "kgbsize: " << rf->val.KGB_size() << std::endl;
  const KGB& kgb=rf->kgb();
  kgb_io::var_print_KGB(*output_stream,rf->val.innerClass(),kgb);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}
void print_KGB_selection_wrapper(expression_base::level l)
{ shared_row selection = get<row_value>();
  shared_real_form rf= get<real_form_value>();
@)
  KGBEltList sel; sel.reserve(selection->val.size());
  for (auto p : selection->val)
  { const auto* x = force<KGB_elt_value>(p.get());
    if (x->rf!=rf)
      throw runtime_error("Real form mismatch when printing KGB element");
    sel.push_back(x->val);
  }
@)
  const KGB& kgb=rf->kgb();
  kgb_io::print(*output_stream,kgb,false,&rf->val.innerClass(),&sel);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}
@)
void print_KGB_order_wrapper(expression_base::level l)
{ shared_real_form rf= get<real_form_value>();
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
{ shared_real_form rf= get<real_form_value>();
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
{ shared_Block b = get<Block_value>();
  Block& block = b->val;
@)
  b->kl_tab.fill(); // this does the actual KL computation
  *output_stream
    << "Full list of non-zero Kazhdan-Lusztig-Vogan polynomials:\n\n";
  kl_io::printAllKL(*output_stream,b->kl_tab,block);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@ The function |print_prim_KL| is a variation of |print_KL_basis|.

@< Local function def...@>=
void print_prim_KL_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  Block &block = b->val; // this one must be non-|const|
@)
  b->kl_tab.fill(); // this does the actual KL computation
  *output_stream
    << "Non-zero Kazhdan-Lusztig-Vogan polynomials for primitive pairs:\n\n";
  kl_io::printPrimitiveKL(*output_stream,b->kl_tab,block);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@ The function |print_KL_list| is another variation of |print_KL_basis|, it
outputs just a list of all distinct Kazhdan-Lusztig-Vogan polynomials.

@< Local function def...@>=
void print_KL_list_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
@)
  b->kl_tab.fill(); // this does the actual KL computation
  kl_io::printKLList(*output_stream,b->kl_tab);
@)
  if (l==expression_base::single_value)
    wrap_tuple<0>();
}

@ We close with two functions for printing the $W$-graph determined by the
polynomials computed. For |print_W_cells| we must construct one more object,
after having built the |kl_tab::KL_table|.

@h "wgraph.h"
@h "wgraph_io.h"

@< Local function def...@>=
void print_W_cells_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
@)
  b->kl_tab.fill(); // this does the actual KL computation
  wgraph::WGraph wg = kl::wGraph(b->kl_tab);
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
{ shared_Block b = get<Block_value>();
@)
  b->kl_tab.fill(); // this does the actual KL computation
  wgraph::WGraph wg = kl::wGraph(b->kl_tab);
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
install_function(print_KGB_selection_wrapper,@|"print_KGB"
		,"(RealForm,[KGBElt]->)");
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
  coercion(KType_type,KTypePol_type,"KpolK",K_type_to_poly);
  coercion(param_type,param_pol_type,"PolP",param_to_poly);
}



@* Index.

% Local IspellDict: british
