% Copyright (C) 2006-2012 Marc van Leeuwen
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
\def\Z{{\bf Z}}
\def\lcm{\mathop{\rm lcm}}
\let\eps=\varepsilon
\def\rk{\mathop{\rm rk}}

@* Built-in types.
%
This file describes several built-in types related to the Atlas software,
which used by the interpreter as primitive types. It also defines built-in
functions that relate to these types. This part of the program is completely
oriented towards the Atlas library, whereas the \.{evaluator} does things that
are for the most part independent of that library.

@h "built-in-types.h"

@f lambda NULL
@f pi NULL
@f alpha NULL
@f beta NULL

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
{ @< Install wrapper functions @>
  @< Install coercions @>
}

@ Before we can define any types we must make sure the types defined
in \.{types.w} from which we shall derive others are known. Including this as
first file from our header file ensures the types are known wherever they are
needed. Some more basic built-in types likes integers, vectors, and strings
are defined in \.{global.w}, and we need their declarations too in our
implementation, but we avoid including its header into out header file.

@h "global.h"

@< Includes needed in the header file @>=
#include "types.h"

@*1 Lie types.
Our first chapter concerns Lie types, as indicated by strings like
|"D4.A3.E8.T1"|. We base ourselves on the related types defined in the Atlas
library.

@< Includes needed in the header file @>=
#include "atlas_types.h"
#include <stdexcept>
#include "lietype.h"

@*2 The primitive type.
%
A first new type corresponds to the type |LieType| in the Atlas
library. We provide a constructor that incorporates a complete
|LieType| value, but often one will use the default constructor and
the |add| method. The remaining methods are obligatory for a primitive type.

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
  void add_simple_factor (char,size_t)
    throw(std::bad_alloc,std::runtime_error); // grow
private:
  Lie_type_value(const Lie_type_value& v) : val(v.val) @+{}
    // copy constructor, used by |clone|
};
@)
typedef std::auto_ptr<Lie_type_value> Lie_type_ptr;
typedef std::tr1::shared_ptr<Lie_type_value> shared_Lie_type;

@ The type |LieType| is publicly derived from |std::vector<SimpleLieType>|,
and in its turn the type |SimpleLieType| is publicly derived from
|std::pair<char,size_t>|. Therefore these types could take arbitrary values,
not necessarily meaningful ones. To remedy this we make the method
|add_simple_factor|, which is the only proposed way to build up Lie types,
check for the validity.

Since the tests defined in \.{io/interactive\_lietype.cpp} used in the
current interface for the Atlas software are clumsy to use, we perform
our own tests here, emulating |interactive_lietype::checkSimpleLieType|.
Torus factors of rank $r>1$ should be equivalent to $r$ torus factors of
rank~$1$, and it simplifies the software if we rewrite the former form to the
latter on input, so we do that here.

@h "constants.h"

@< Function definitions @>=
void Lie_type_value::add_simple_factor (char c,size_t rank)
   throw(std::bad_alloc,std::runtime_error)
{ static const std::string types=lietype::typeLetters; // |"ABCDEFGT"|
  size_t t=types.find(c);
  if (t==std::string::npos)
    throw std::runtime_error(std::string("Invalid type letter '")+c+'\'');
@.Invalid type letter@>
  const size_t r=constants::RANK_MAX;
@/static const size_t lwb[]={1,2,2,4,6,4,2,0};
  static const size_t upb[]={r,r,r,r,8,4,2,r};
  if (rank<lwb[t])
    throw std::runtime_error("Too small rank "+str(rank)+" for Lie type "+c);
@.Too small rank@>
  if (rank>upb[t])
    if (upb[t]!=r)
      throw std::runtime_error("Too large rank "+str(rank)+" for Lie type "+c);
@.Too large rank@>
    else
      throw std::runtime_error @|
      ("Rank "+str(rank)+" exceeds implementation limit "+str(r));
@.Rank exceeds implementation limit@>
  if (c=='T')
    while (rank-->0) val.push_back(SimpleLieType('T',1));
  else
    val.push_back(SimpleLieType(c,rank));
}

@ Now we define a wrapper function that really builds a |Lie_type_value|. We
scan the string looking for sequences of a letter followed by a number. We
allow and ignore punctuation characters between simple factors, also spaces
are allowed anywhere except inside the number, since formatted input from
streams, even of characters, by default skips spaces. Since the correct
structure of Lie type strings is so obvious to the human eye, our error
message just cites the entire offending string, rather than trying to point
out the error exactly.

@h <sstream>
@h <cctype>
@< Local function definitions @>=
inline std::istream& skip_punctuation(std::istream &is)
{@; char c; do is>>c; while (is and std::ispunct(c));
    return is.unget();
}
@)
void Lie_type_wrapper(expression_base::level l)
   throw(std::bad_alloc,std::runtime_error)
{ std::istringstream is(get<string_value>()->val);
  Lie_type_ptr result(new Lie_type_value);
  size_t total_rank=0; char c;
  while (skip_punctuation(is)>>c) // i.e., until |not is.good()|
  { size_t rank;
    if (is.peek()=='-' or not (is>>rank)) // explicitly forbid minus sign
      throw std::runtime_error
       ("Error in type string '"+is.str()+"' for Lie type");
@.Error in type string@>
    result->add_simple_factor(c,rank);
      // this may |throw| a |runtime_error| as well
    if ((total_rank+=rank)>constants::RANK_MAX)
      throw std::runtime_error
      ("Total rank exceeds implementation limit "+str(constants::RANK_MAX));
@.Total rank exceeds...@>
  }
  if (l!=expression_base::no_value)
    push_value(result);
}
@)
void Lie_type_coercion()
@+{@; Lie_type_wrapper(expression_base::single_value); }

@ We shall call this function \.{Lie\_type}.

@< Install wrapper functions @>=
install_function(Lie_type_wrapper,"Lie_type","(string->LieType)");

@ In fact we shall also make it available as an implicit conversion. This is
very easy to implement.

@< Install coercions @>=



@*2 Auxiliary functions for Lie types.
Before we do anything more complicated with this primitive type, we must
ensure that we can print its values. We can use an operator defined in
\.{basic\_io.cpp}.

@h "basic_io.h"
@< Function definitions @>=
void Lie_type_value::print(std::ostream& out) const
{ if (val.empty()) out << "empty Lie type";
  else
    out << "Lie type '" << val << '\'';
}

@ Here is a function that computes the Cartan matrix for a given Lie type.

@h "prerootdata.h"
@< Local function definitions @>=
void Cartan_matrix_wrapper(expression_base::level l)
{ shared_Lie_type t=get<Lie_type_value>();
  if (l==expression_base::no_value)
    return;
  matrix_ptr result(new matrix_value(t->val.Cartan_matrix()));
  push_value(result);
}


@ And here is a function that tries to do the inverse. We call
|dynkin::lieType| in its version that also produces a permutation |pi| needed
to map the standard ordering for the type to the actual ordering, and |true|
as third argument to request a check that the input was indeed the Cartan
matrix for that situation (if not a runtime error will be thrown). Without
this test invalid Lie types could have been formed, for which root datum
construction would most likely crash.

@h "dynkin.h"
@< Local function definitions @>=
void type_of_Cartan_matrix_wrapper (expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  Permutation pi;
  LieType lt=dynkin::Lie_type(m->val,true,true,pi);
  if (l==expression_base::no_value)
    return;
  push_value(new Lie_type_value(lt));
  push_value(new vector_value(CoeffList(pi.begin(),pi.end())));
  if (l==expression_base::single_value)
    wrap_tuple(2);
}

@ And some small utilities for finding the (Lie) rank and semisimple rank of
the Lie type, and the naked type string.

@< Local function definitions @>=
void Lie_rank_wrapper(expression_base::level l)
{ shared_Lie_type t=get<Lie_type_value>();
  if (l!=expression_base::no_value)
    push_value(new int_value(t->val.rank()));
}
void semisimple_rank_wrapper(expression_base::level l)
{ shared_Lie_type t=get<Lie_type_value>();
  if (l!=expression_base::no_value)
    push_value(new int_value(t->val.semisimple_rank()));
}
void Lie_type_string_wrapper(expression_base::level l)
{ std::ostringstream s; s << get<Lie_type_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new string_value(s.str()));
}

@ For programming it is important to know the number of factors in a
Lie type, so that for instance the correct number of inner class letters can
be prepared (for this purpose type $T_n$ counts as $n$ factors). Since we
expanded any $T_n$ into factors $T_1$, we can simply call the |size| method of
the stored |LieType| value. To allow per-factor treatment of Lie type, we also
define a function that converts a Lie type into a row of Lie types, one for
every simple factor or torus factor.

@< Local function definitions @>=
void nr_factors_wrapper(expression_base::level l)
{ shared_Lie_type t=get<Lie_type_value>();
  if (l!=expression_base::no_value)
    push_value(new int_value(t->val.size()));
}
@)
void Lie_factors_wrapper(expression_base::level l)
{ shared_Lie_type t=get<Lie_type_value>();
  if (l==expression_base::no_value)
    return;
  row_ptr result(new row_value(t->val.size()));
  for (unsigned i=0; i<t->val.size(); ++i)
  { std::vector<SimpleLieType> factor(1,t->val[i]);
    result->val[i]=shared_Lie_type(new Lie_type_value(LieType(factor)));
  }
  push_value(result);
}

@ Again we install our wrapper functions.
@< Install wrapper functions @>=
install_function(Cartan_matrix_wrapper,"Cartan_matrix","(LieType->mat)");
install_function(type_of_Cartan_matrix_wrapper
		,@|"Cartan_matrix_type","(mat->LieType,vec)");
install_function(Lie_rank_wrapper,"rank","(LieType->int)");
install_function(semisimple_rank_wrapper,"semisimple_rank","(LieType->int)");
install_function(Lie_type_string_wrapper,"str","(LieType->string)");
install_function(nr_factors_wrapper,"#","(LieType->int)");
install_function(Lie_factors_wrapper,"%","(LieType->[LieType])");

@*2 Finding lattices for a given Lie type.
%
When a Lie type is fixed, there is still a nontrivial choice to determine the
root datum for a connected complex reductive group of that type: one has to
choose a sub-lattice of the weight lattice of the ``simply connected'' group
of that type to become the weight lattice of the chosen Complex reductive
group (a finite quotient of the simply connected group); that sub-lattice
should have full rank and contain the root lattice. We shall start with
considering some auxiliary functions to help facilitate this choice.

Here is a wrapper function around the |LieType::Smith_basis| method, which
computes a block-wise Smith basis for the transposed Cartan matrix, and the
corresponding invariant factors. In case of torus factors this description
should be interpreted in the sense that the Smith basis for those factors is
the standard basis and the invariant factors are null.

@< Local function definitions @>=
void Smith_Cartan_wrapper(expression_base::level l)
{ shared_Lie_type t=get<Lie_type_value>();
  if (l==expression_base::no_value)
    return;
  vector_ptr inv_factors (new vector_value(CoeffList()));
  WeightList b = t->val.Smith_basis(inv_factors->val);
  push_value(new matrix_value(LatticeMatrix(b,b.size())));
  push_value(inv_factors);
  if (l==expression_base::single_value)
     wrap_tuple(2);
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
{ shared_vector inv_f=get_own<vector_value>();
  shared_matrix basis=get_own<matrix_value>();
  if (inv_f->val.size()!=basis->val.numColumns())
    throw std::runtime_error @|("Size mismatch "+
@.Size mismatch@>
      str(inv_f->val.size())+':'+str(basis->val.numColumns()));
  if (l==expression_base::no_value)
    return;
@)
  size_t i=0;
  while (i<inv_f->val.size())
    if (inv_f->val[i]!=1) ++i; // keep invariant factor and column
    else
    {@; inv_f->val.erase(inv_f->val.begin()+i);
        basis->val.eraseColumn(i);
    }
  push_value(basis); push_value(inv_f);
  if (l==expression_base::single_value)
      wrap_tuple(2);
}

@ Here is another function, adapted from the functions |makeOrthogonal| and
|getLattice|, defined locally in the file \.{io/interactive\_lattice.cpp}. We
have taken the occasion to change the name and interface and everything else,
which also avoids the need to introduce rational matrices as primitive type.

The function |annihilator_modulo| takes as argument an $m\times{n}$
matrix~$M$, and an integer |d|. It returns a $m\times{m}$ matrix~|A| whose
columns span the full rank sub-lattice of $\Z^m$ of vectors $v$ such that
$v^t\cdot{M}\in d\,\Z^n$. The fact that $d$ is called |denominator| below
comes from the alternative interpretation that the transpose of~$A$ times the
rational matrix $M/d$ gives an integral matrix as result.

The algorithm is quite simple. After finding invertible matrices |row|, |col|
such that $row*M*col$ is diagonal with non-zero diagonal entries given in
$\lambda$ (and any zero diagonal entries trailing those), we know that any row
of |row| with factor $\lambda_i$ is a linear form sending the image of $M$ to
$\lambda_i\Z$, while any remaining row of |row| annihilates the image
altogether. Then all that is needed it to multiply rows of~|row| of the first
kind by $d/\gcd(d,\lambda_i)$ and transpose the result.

@h "arithmetic.h"
@h "lattice.h"
@h "matrix.h"
@h "matreduc.h"

@< Local function definitions @>=
LatticeMatrix @|
annihilator_modulo(const LatticeMatrix& M, LatticeCoeff denominator)

{ int_Matrix row,col;
  CoeffList lambda = matreduc::diagonalise(M,row,col);

  for (size_t i=0; i<lambda.size(); ++i)
  @/{@; unsigned long c = arithmetic::div_gcd(denominator,lambda[i]);
    row.rowMultiply(i,c);
  }
  return row.transposed();
}

@ The wrapper function is particularly simple.

@< Local function definitions @>=
void ann_mod_wrapper(expression_base::level l)
{ shared_int d=get<int_value>();
  shared_matrix m=get<matrix_value>();
@)
  if (l!=expression_base::no_value)
    push_value(new matrix_value(annihilator_modulo(m->val,d->val)));
}

@ Next a simple administrative routine, needed here because we cannot handle
matrices in our programming language yet. Once one has computed a new lattice
(possibly with the help of |ann_mod|), in the form of vectors to replace those
selected by |filter_units| from the result of |Smith_Cartan|, one needs to
make the replacement. The following function does this, taking as its first
argument the result of |Smith_Cartan|, and as second a matrix whose columns
are to be substituted. The invariant factors in the first argument serve only
to determine, by the place the non-unit ones, where the insertion has to take
place. In fact this is so simple that we define the wrapper function directly.

@< Local function definitions @>=
void replace_gen_wrapper (expression_base::level l)
{ shared_matrix new_generators=get<matrix_value>();
  push_tuple_components(); // a pair as returned by \.{Smith\_Cartan}
  shared_vector inv_f=get<vector_value>();
  shared_matrix generators=get_own<matrix_value>();
   // start with old generators
@)
  if (new_generators->val.numRows()!=generators->val.numRows())
    throw std::runtime_error("Column lengths do not match");
@.Column lengths do not match@>
  if (inv_f->val.size()!=generators->val.numColumns())
    throw std::runtime_error("Number of columns mismatch");
@.Number of columns mismatch@>
@)
  size_t k=0; // index to replacement generators
  for (size_t j=0; j<inv_f->val.size(); ++j)
    if (inv_f->val[j]!=1)
       // replace column |j| by column |k| from |new_generators|
    { if (k>=new_generators->val.numColumns())
        throw std::runtime_error ("Not enough replacement columns");
@.Not enough replacement columns@>
      generators->val.set_column(j,new_generators->val.column(k));
      ++k;
    }
  if (k<new_generators->val.numColumns())
        throw std::runtime_error ("Too many replacement columns");
@.Too many replacement columns@>
  if (l!=expression_base::no_value)
    push_value(generators);
}

@*2 Specifying inner classes. Now we move ahead a bit in the theory, from
functions that help in building root data to functions that help defining
(inner classes of) real forms. The first of such functions is
|lietype::involution|, which takes a Lie type and a |InnerClassType|
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
an inner class, and transforms it into |InnerClassType| that is
guaranteed to be valid if returned; the routine throws a |runtime_error| in
case of problems.

@< Local function definitions @>=
InnerClassType transform_inner_class_type
  (const char* s, const LieType& lt)
throw (std::bad_alloc, std::runtime_error)
{ static const std::string types(lietype::innerClassLetters);
    // |"Ccesu"|
  InnerClassType result; // initially empty
  std::istringstream is(s);
  char c;
  size_t i=0; // position in simple factors of Lie type |lt|
  while (skip_punctuation(is)>>c)
    @< Test the inner class letter |c|, and either push a corresponding type
       letter onto |result| while advancing~|i| by the appropriate amount, or
       throw a |runtime_error| @>
  if (i< lt.size()) throw std::runtime_error("Too few inner class symbols");
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
{ if (types.find(c) == std::string::npos) throw std::runtime_error@|
    (std::string("Unknown inner class symbol `")+c+"'");
  if (i>= lt.size()) throw std::runtime_error("Too many inner class symbols");
@.Too many inner class symbols@>
  lietype::TypeLetter t = lt[i].type(); size_t r=lt[i].rank();
  if (c=='C') // complex inner class, which requires two equal simple factors
  { if (i+1>=lt.size() or lt[i+1]!=lt[i]) throw std::runtime_error @|
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
  else throw std::runtime_error @|
    (std::string("Unequal rank class is meaningless for type ")+t+str(r));
@.Unequal rank class is meaningless...@>
}

@ The wrapper function around |lietype::involution| will take a Lie type and a
string of type letters and return a matrix describing the involution
designated by that string, expressed on the fundamental weight basis for the
simply connected group of that type.

@< Local function def... @>=
void basic_involution_wrapper(expression_base::level l)
{ shared_string str=get<string_value>();
  shared_Lie_type t=get<Lie_type_value>();
@/push_value(new matrix_value @| (lietype::involution
           (t->val,transform_inner_class_type(str->val.c_str(),t->val))));
}

@ The function just defined gives an involution on the basis of fundamental
weights for a simply connected group, or on the basis of the root lattice for
an adjoint group. For a general complex groups however, we need to transform
this involution to an involution of a given sub-lattice, assuming that this is
possible, in other words if that the sub-lattice is globally stable under the
involution specified. Therefore we now provide a function that in addition to
the Lie type takes a matrix specifying a sub-lattice as argument, and finally
a string specifying the inner class.

The fist two ingredients are also those used to construct a root datum, and
one might imagine replacing them by a root datum. The function
$set\_inner\_class$ defined later will do that, but it has to accept some
ambiguity in recovering Lie type and sub-lattice from a root datum. Also it
may find a permutation of the simple roots with respect to the standard
ordering of the diagram, and has to deal with that additional generality. So
unfortunately it can neither replace nor call the function $based\_involution$
defined here; it will have to adapt some of the work done here to its
situation.


@< Local function def... @>=
void based_involution_wrapper(expression_base::level l)
{ shared_string s = get<string_value>();
@/shared_matrix basis = get<matrix_value>();
@/shared_Lie_type type = get<Lie_type_value>();
@)
  size_t r=type->val.rank();
  if (basis->val.numRows()!=r or basis->val.numRows()!=r)
    throw std::runtime_error @|
    ("Basis should be given by "+str(r)+'x'+str(r)+" matrix");
@.Basis should be given...@>
@)
  WeightInvolution inv=lietype::involution
        (type->val,transform_inner_class_type(s->val.c_str(),type->val));
  try
  {@; push_value(new matrix_value(inv.on_basis(basis->val.columns()))); }
  catch (std::runtime_error&) // relabel |"Inexact integer division"|
  {@; throw std::runtime_error
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
install_function(basic_involution_wrapper,"involution",
		"(LieType,string->mat)");
install_function(based_involution_wrapper,"involution",
		"(LieType,mat,string->mat)");

@*1 Root data.
Now we are ready to introduce a new primitive type for root data.

@< Includes needed in the header file @>=
#include "rootdata.h"

@ The root datum type is laid out just like previous primitive types are.

@< Type definitions @>=
struct root_datum_value : public value_base
{ RootDatum val;
@)
  root_datum_value(const RootDatum& v) : val(v) @+ {}
  virtual void print(std::ostream& out) const;
  root_datum_value* clone() const @+{@; return new root_datum_value(*this); }
  static const char* name() @+{@; return "root datum"; }
private:
  root_datum_value(const root_datum_value& v) : val(v.val) @+{}
};
@)
typedef std::auto_ptr<root_datum_value> root_datum_ptr;
typedef std::tr1::shared_ptr<root_datum_value> shared_root_datum;

@*2 Printing root data. We shall not print the complete information contained
in the root datum. However we do exercise the routines in \.{dynkin.cpp} to
determine the type of the root datum from its Cartan matrix, as is done in the
function |type_of_Cartan_matrix| above. Contrary to |LieType::Cartan_matrix|,
the function |RootDatum::cartanMatrix| produces nothing for torus factors (in
fact only their number is well defined); hence we add the necessary number of
torus factors at the end.

@< Local fun... @>=
LieType type_of_datum(const RootDatum& rd)
{ int_Matrix Cartan = rd.cartanMatrix();
@/LieType t = dynkin::Lie_type(Cartan);
  if (!rd.isSemisimple())
    for (size_t i=rd.semisimpleRank(); i<rd.rank(); ++i)
      t.push_back(SimpleLieType('T',1));
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
void type_of_root_datum_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  push_value(new Lie_type_value(type_of_datum(rd->val)));
}

@*2 Building a root datum.
%
To create a root datum value, the user may specify a Lie type and a square
matrix of the size of the rank of the root datum, which specifies generators
of the desired weight lattice as a sub-lattice of the lattice of weights
associated to the simply connected group of the type given. The given weights
should be independent and span at least the root lattice associated to the
type. Failure of either condition will cause the |PreRootDatum| constructor to
throw a |std::runtime_error|, but in case a root should fail to be expressed
in the sub-lattice basis this happens in a call buried fairly deeply, so we
prefer to replace the |"Inexact integer division"| message by a somewhat more
descriptive one.

@< Local function definitions @>=
void root_datum_wrapper(expression_base::level l)
{ shared_matrix lattice=get<matrix_value>();
  shared_Lie_type type=get<Lie_type_value>();
  if (lattice->val.numRows()!=lattice->val.numColumns() @| or
      lattice->val.numRows()!=type->val.rank())
    throw std::runtime_error
    ("Sub-lattice matrix should have size " @|
@.Sub-lattice matrix should...@>
      +str(type->val.rank())+'x'+str(type->val.rank()));
  try
  {
    PreRootDatum prd(type->val,lattice->val.columns());
    if (l!=expression_base::no_value)
      push_value(new root_datum_value @| (RootDatum(prd)));
  }
  catch (std::runtime_error& e)
  { if (e.what()[0]=='I') // |"Inexact integer division"|
      throw std::runtime_error
        ("Sub-lattice does not contain the root lattice");
@.Sub-lattice does not contain...@>
    throw; // |"Dependent lattice generators"|, no need to rephrase
@.Dependent lattice generators@>
  }
}

@ Alternatively, a user may just specify bases of simple roots and coroots in
the form of matrices, which implicitly define a Lie type and weight lattice.
To make sure the root datum construction will succeed, we must test the
``Cartan'' matrix computed from these data to be a valid one.

@< Local function definitions @>=
void raw_root_datum_wrapper(expression_base::level l)
{ shared_matrix simple_coroots=get<matrix_value>();
  shared_matrix simple_roots=get<matrix_value>();

  size_t nr = simple_roots->val.numRows(),
         nc = simple_roots->val.numColumns();

  if (simple_coroots->val.numRows()!=nr @| or
      simple_coroots->val.numColumns()!=nc)
    throw std::runtime_error
    ("Sizes (" +str(nr)+ "," +str(nc) +"),("+
      str(simple_coroots->val.numRows()) +","+
      str(simple_coroots->val.numColumns()) +
      ") of simple (co)root systems differ");
@.Sizes of simple (co)root systems...@>

  WeightList s; CoweightList c;
  s.reserve(nc);
  c.reserve(nc);

  for (size_t j=0; j<nc; ++j)
@/{@; s.push_back(simple_roots->val.column(j));
      c.push_back(simple_coroots->val.column(j));
  }

  PreRootDatum prd(s,c,nr);
  try @/{@; Permutation dummy;
    dynkin::Lie_type(prd.Cartan_matrix(),true,true,dummy);
  }
  catch (error::CartanError)
@/{@;
    throw std::runtime_error("System of (co)roots has invalid Cartan matrix");
}
@.System of (co)roots has invalid...@>
  if (l!=expression_base::no_value)
    push_value(new root_datum_value @| (RootDatum(prd)));
}

@ To emulate what is done in the Atlas software, we write a function that
integrates some of the previous ones. It is called as |quotient_basis(lt,L)|
where $t$ is a Lie type, and $L$ is a list of rational vectors (interpreted as
coweights), each giving a kernel generator as would be entered in the Atlas
software.

Let |S==Smith_Cartan(lt)| and |(C,v)=filter\_units(S)|, then we find a basis
for the sub-lattice needed to build the root datum as follows. Each vector
in~$L$ should have the same length as~$v$, and multiplication of corresponding
entries should always give an integer. Then a common denominator~$d$ is found
and a matrix $M$ whose columns form the numerators of the lists of~$L$ brought
to the denominator~$d$. The call $quotient\_basis(lt,L)$ then yields the result
of computing $replace\_gen(S,C*ann\_mod(M,d))$.

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
  wrap_tuple(2); // |replace_gen| wants as first argument a tuple
@/push_value(SC_basis);
  push_value(invf);
  filter_units_wrapper(expression_base::multi_value);
  shared_vector v=get<vector_value>();
  // and leave |C| for call to $mm\_prod$
@)
  LatticeMatrix M(v->val.size(),L->length());
  size_t d=1;
  @< Compute common denominator |d| of entries in~$L$, and place converted
     denominators into the columns of $M$; also test validity of entries and
     |throw| a runtime error for invalid ones @>
  push_value(new matrix_value(annihilator_modulo(M,d)));
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
{ std::vector<unsigned long> denom(L->length());
  for (size_t j=0; j<L->length(); ++j)
  { const RatWeight& gen =
      force<rational_vector_value>(&*L->val[j])->val;
    denom[j] = gen.denominator();
    d=arithmetic::lcm(d,denom[j]);

    if (gen.size()!=v->val.size())
      throw std::runtime_error@|
        ("Length mismatch for generator "+str(j) +": "@|
@.Length mismatch...@>
        +str(gen.size()) + ':' + str(v->val.size()));

    Weight col(gen.numerator().begin(),gen.numerator().end()); // convert
    M.set_column(j,col);
    for (size_t i=0; i<v->val.size(); ++i)
      if (v->val[i]*M(i,j)%long(denom[j])!=0) // must use signed arithmetic!!
	throw std::runtime_error("Improper generator entry: "
@.Improper generator entry@>
         +str(M(i,j))+'/'+str(denom[j])+" not a multiple of 1/"
         +str(v->val[i]));
  }
@)
  for (size_t j=0; j<L->length(); ++j)
    // convert modulo $\Z$ and to common denominator |d|
  { size_t f=d/denom[j];
    for (size_t i=0; i<v->val.size(); ++i)
      M(i,j) = arithmetic::remainder(M(i,j),denom[j])*f;
  }
}

@ The function that integrates all is $quotient\_datum$; the call
$quotient\_datum(lt,L)$ is equivalent to
$root\_datum(lt,quotient\_basis(lt,L))$.

@< Local function definitions @>=
void quotient_datum_wrapper(expression_base::level l)
{ shared_value L = pop_value();
  shared_value lt= pop_value();
  push_value(lt); // the Lie type, for call of $root\_datum$
@/push_value(lt);
  push_value(L); quotient_basis_wrapper(expression_base::single_value);
@/root_datum_wrapper(l); // pass level parameter to final call
}

@ We define two more wrappers with only a Lie type as argument, for building
the simply connected and the adjoint root data. They are similar to the
previous one in that they mostly call other wrapper functions: the matrices
that specify the sub-lattices are produced by the wrapper functions for
|id_mat| respectively for |Cartan_matrix| and |transpose_mat|. In fact the
call $simply\_connected\_datum(lt)$ is equivalent to
$root\_datum(lt,id\_mat(Lie\_rank(lt)))$, and $adjoint\_datum(lt)$ is almost
equivalent to $root\_datum(lt,M)$ where
$M=transpose\_mat(Cartan\_matrix(lt)))$; the only thing we do ``by hand'' is
to make sure that all null diagonal entries of~$M$ (which must come from torus
factors) are replaced by ones.

@< Local function definitions @>=
void simply_connected_datum_wrapper(expression_base::level l)
{ if (l==expression_base::no_value)
  { execution_stack.pop_back();
    return; // no possibilities of errors, so avoid useless work
  }
  size_t rank =
    force<Lie_type_value>(execution_stack.back().get())->val.rank();
  push_value(new int_value(rank));
  id_mat_wrapper(expression_base::single_value);
@/root_datum_wrapper(expression_base::single_value);
}
@)
void adjoint_datum_wrapper(expression_base::level l)
{ if (l==expression_base::no_value)
  { execution_stack.pop_back();
    return; // no possibilities of errors, so avoid useless work
  }
  push_value(execution_stack.back()); // duplicate Lie type argument
  Cartan_matrix_wrapper(expression_base::single_value);
  transpose_mat_wrapper(expression_base::single_value);
  shared_matrix M=get_own<matrix_value>();
  for (size_t i=0; i<M->val.numRows(); ++i)
    if (M->val(i,i)==0) M->val(i,i)=1;
  push_value(M);
@/root_datum_wrapper(expression_base::single_value);
}


@*2 Functions operating on root data.
%
The following functions allow us to look at the simple roots and simple
coroots stored in a root datum value, and at the associated Cartan matrix.

@< Local function definitions @>=
void simple_roots_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
  WeightList srl
    (rd->val.beginSimpleRoot(),rd->val.endSimpleRoot());
  push_value(new matrix_value
    (int_Matrix(srl,rd->val.rank())));
}
@)
void simple_coroots_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
  WeightList scl
    (rd->val.beginSimpleCoroot(),rd->val.endSimpleCoroot());
  push_value(new matrix_value
    (int_Matrix(scl,rd->val.rank())));
}
@)
void datum_Cartan_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
  int_Matrix M = rd->val.cartanMatrix();
  push_value(new matrix_value(M));
}

@ We also allow access to all positive (co)roots, or all of them (positive or
negative).

@< Local function definitions @>=
void positive_roots_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
  WeightList rl
    (rd->val.beginPosRoot(),rd->val.endPosRoot());
  push_value(new matrix_value
    (int_Matrix(rl,rd->val.rank())));
}
@)
void positive_coroots_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
  WeightList crl
    (rd->val.beginPosCoroot(),rd->val.endPosCoroot());
  push_value(new matrix_value
    (int_Matrix(crl,rd->val.rank())));
}
void roots_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
  WeightList rl
    (rd->val.beginRoot(),rd->val.endRoot());
  push_value(new matrix_value(int_Matrix(rl,rd->val.rank())));
}
@)
void coroots_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
  WeightList crl
    (rd->val.beginCoroot(),rd->val.endCoroot());
  push_value(new matrix_value(int_Matrix(crl,rd->val.rank())));
}

@ Here are three important numeric attributes of root data.

@< Local function definitions @>=
void rd_rank_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l!=expression_base::no_value)
    push_value(new int_value(rd->val.rank()));
}
@)
void rd_semisimple_rank_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l!=expression_base::no_value)
    push_value(new int_value(rd->val.semisimpleRank()));
}
@)
void rd_nposroots_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l!=expression_base::no_value)
    push_value(new int_value(rd->val.numPosRoots()));
}

@ It is useful to have bases for the sum of the root lattice and the
coradical, and for the sum of the coroot lattice and the radical; the latter
will in fact be used later in a function to built inner classes.

@< Local function definitions @>=
void root_coradical_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
  WeightList srl
    (rd->val.beginSimpleRoot(),rd->val.endSimpleRoot());
  srl.insert(srl.end(),rd->val.beginCoradical(),rd->val.endCoradical());
  push_value(new matrix_value(int_Matrix(srl,rd->val.rank())));
}
@)
void coroot_radical_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l==expression_base::no_value)
    return;
  WeightList scl
    (rd->val.beginSimpleCoroot(),rd->val.endSimpleCoroot());
  scl.insert(scl.end(),rd->val.beginRadical(),rd->val.endRadical());
  push_value(new matrix_value(int_Matrix(scl,rd->val.rank())));
}

@ We give access to the fundamental weights and coweights on an individual
basis, which is easier since they are rational vectors.

@< Local function definitions @>=
void fundamental_weight_wrapper(expression_base::level l)
{ int i= get<int_value>()->val;
  shared_root_datum rd(get<root_datum_value>());
  if (unsigned(i)>=rd->val.semisimpleRank())
    throw std::runtime_error("Invalid index "+str(i));
  if (l!=expression_base::no_value)
    push_value(new rational_vector_value
      (rd->val.fundamental_weight(i).normalize()));
}
@)
void fundamental_coweight_wrapper(expression_base::level l)
{ int i= get<int_value>()->val;
  shared_root_datum rd(get<root_datum_value>());
  if (unsigned(i)>=rd->val.semisimpleRank())
    throw std::runtime_error("Invalid index "+str(i));
  if (l!=expression_base::no_value)
    push_value(new rational_vector_value
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
{ shared_root_datum rd(get<root_datum_value>());
  if (l!=expression_base::no_value)
    push_value(new root_datum_value@|
      (RootDatum(rd->val,tags::DualTag())));
}
@)
void derived_info_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l!=expression_base::no_value)
  { int_Matrix projector;
    push_value(new root_datum_value@|
      (RootDatum(projector,rd->val,tags::DerivedTag())));
    push_value(new matrix_value(projector));
    if (l==expression_base::single_value)
      wrap_tuple(2);
  }
}
@)
void mod_central_torus_info_wrapper(expression_base::level l)
{ shared_root_datum rd(get<root_datum_value>());
  if (l!=expression_base::no_value)
  { int_Matrix injector;
    push_value(new root_datum_value@|
      (RootDatum(injector,rd->val,tags::AdjointTag())));
    push_value(new matrix_value(injector));
    if (l==expression_base::single_value)
      wrap_tuple(2);
  }
}

@ This function is more recent; it allows constructing a new root (sub-)datum
by selecting coroots taking integral values on a given rational weight vector.

@< Local function definitions @>=
void integrality_datum_wrapper(expression_base::level l)
{ shared_rational_vector lambda = get<rational_vector_value>();
  shared_root_datum rd(get<root_datum_value>());
  if (lambda->val.size()!=rd->val.rank())
  { std::ostringstream o;
    o << "Length " << lambda->val.size()
      << " of rational vector differs from rank " << rd->val.rank();
    throw std::runtime_error(o.str());
@.Length of rational vector...@>
  }
  if (l!=expression_base::no_value)
    push_value(new root_datum_value @|
      (rootdata::integrality_datum(rd->val,lambda->val)));
}

@ A related function computes a list of fractions of a line segment where the
set of roots with integrality is non-empty.

@< Local function definitions @>=
void integrality_points_wrapper(expression_base::level l)
{ shared_rational_vector lambda = get<rational_vector_value>();
  shared_root_datum rd(get<root_datum_value>());
  if (lambda->val.size()!=rd->val.rank())
  { std::ostringstream o;
    o << "Length " << lambda->val.size()
      << " of rational vector differs from rank " << rd->val.rank();
    throw std::runtime_error(o.str());
@.Length of rational vector...@>
  }
  if (l==expression_base::no_value)
    return;
@)
  RationalList ipl = rootdata::integrality_points(rd->val,lambda->val);
    // method normalises rationals
  row_ptr result (new row_value(ipl.size()));
  for (size_t i=0; i<ipl.size(); ++i)
    result->val[i]=shared_value(new rat_value(ipl[i]));
  push_value(result);
}

@ Let us install the above wrapper functions.

@< Install wrapper functions @>=
install_function(type_of_root_datum_wrapper,@|"Lie_type"
                ,"(RootDatum->LieType)");
install_function(root_datum_wrapper,@|"root_datum","(LieType,mat->RootDatum)");
install_function(raw_root_datum_wrapper,
                 @|"root_datum","(mat,mat->RootDatum)");
install_function(quotient_basis_wrapper
		,@|"quotient_basis","(LieType,[ratvec]->mat)");
install_function(quotient_datum_wrapper
		,@|"root_datum","(LieType,[ratvec]->RootDatum)");
install_function(simply_connected_datum_wrapper
		,@|"simply_connected","(LieType->RootDatum)");
install_function(adjoint_datum_wrapper,@| "adjoint","(LieType->RootDatum)");
install_function(simple_roots_wrapper,@|"simple_roots","(RootDatum->mat)");
install_function(simple_coroots_wrapper,@|"simple_coroots","(RootDatum->mat)");
install_function(positive_roots_wrapper,@| "posroots","(RootDatum->mat)");
install_function(positive_coroots_wrapper,@| "poscoroots","(RootDatum->mat)");
install_function(datum_Cartan_wrapper,@|"Cartan_matrix","(RootDatum->mat)");
install_function(roots_wrapper,@|"roots","(RootDatum->mat)");
install_function(coroots_wrapper,@|"coroots","(RootDatum->mat)");
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
install_function(rd_nposroots_wrapper@|
		,"nr_of_posroots","(RootDatum->int)");
install_function(integrality_datum_wrapper
                ,@|"integrality_datum","(RootDatum,ratvec->RootDatum)");
install_function(integrality_points_wrapper
                ,@|"integrality_points","(RootDatum,ratvec->[rat])");

@*1 A type for complex reductive groups equipped with an involution.
We shall now go ahead to define a primitive type holding a
|ComplexReductiveGroup|, which represents a complex reductive
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
  (const WeightInvolution& M)
throw (std::bad_alloc, std::runtime_error)
{ size_t r=M.numRows();
  @< Check that |M| is an $r\times{r}$ matrix defining an involution @>
  tori::RealTorus Tor(M);
  return std::make_pair(Tor.compactRank(),Tor.splitRank());
}

@ The test below that |M| is an involution ($M^2=\\{Id}$) is certainly
necessary when |classify_involution| is called independently. If the code
below is executed in the context of checking the involution for an inner
class, it may seem redundant given the fact that we shall also check that |M|
induces an involutive permutation of the roots; however even there it is not
redundant in the presence of a torus part.

@< Check that |M| is an $r\times{r}$ matrix defining an involution @>=
{ if (M.numRows()!=r or M.numColumns()!=r) throw std::runtime_error
    ("Involution should be a "+str(r)+"x"+str(r)+" matrix, got a "
@.Involution should be...@>
     +str(M.numRows())+"x"+str(M.numColumns())+" matrix");
  WeightInvolution Id(r),Q(M*M);
  if (Q!=Id) throw std::runtime_error
      ("Given transformation is not an involution");
@.Given transformation...@>
}

@ The function |classify_involution| might be of use to the user, so let us
make a wrapper for it. In fact we shall return the compact, Complex, and split
ranks, in that order.

@< Local function def...@>=
void classify_wrapper(expression_base::level l)
{ shared_matrix M(get<matrix_value>());
  std::pair<size_t,size_t> p=classify_involution(M->val);
  if (l==expression_base::no_value)
    return;
  push_value(new int_value(p.first)); // compact rank
  push_value(new int_value((M->val.numRows()-p.first-p.second)/2));
    // Complex rank
  push_value(new int_value(p.second)); // split rank
  if (l==expression_base::single_value)
    wrap_tuple(3);
}

@ We now come to the part of the analysis that involves the root datum. Since
the matrix is already expressed on the basis of the weight lattice used by the
root datum, the question of stabilising that lattice is settled, but we must
check that the matrix is indeed an involution, and that it gives an
automorphism of the root datum. In fact we need an automorphism of the based
root datum, but we will allow a conjugate of such an automorphism to be
supplied as~$M$, and export thought the |ww| parameter a Weyl group element
conjugates the based root datum automorphism to~$M$. While we are doing all
this, we can also determine the Lie type of the root datum, the permutation
possibly needed to map the standard (Bourbaki) ordering of that Dynkin diagram
to the actual one, and the inner class letters corresponding to each of its
factors. This is precisely the data stored in a |lietype::Layout| structure
(which is later also used to associate real form names to internal data
(gradings)).

The Lie type will in general be identical to that of the root datum (and in
particular any torus factors will come at the end), but in case of Complex
inner classes we may be forced to permute the simple factors to make the
identical factors associated to such classes adjacent (the permutation will be
adapted to reflect this reordering of simple factors).

We do not apply the function |classify_involution| defined above to the entire
matrix describing a (purported) involution of the weight lattice, although we
shall use it below for a matrix defined for the central torus part; we can
however reuse the module that tests for being an involution here.

@h "weyl.h"

@< Local function def...@>=
lietype::Layout check_involution
 (const WeightInvolution& M, const RootDatum& rd,
  WeylWord& ww)
 throw (std::bad_alloc, std::runtime_error)
{ size_t r=rd.rank(),s=rd.semisimpleRank();
  @< Check that |M| is an $r\times{r}$ matrix defining an involution @>
@/Permutation p(s);
  @< Set |ww| to the reversed Weyl group element needed to be applied after
  the action of |M| in order to map positive roots to positive roots, and |p|
  to the permutation of the simple roots so obtained; throw a |runtime_error|
  if |M| is not an automorphism of |rd| @>
@/lietype::Layout result;
@/LieType& type=result.d_type;
  InnerClassType& inner_class=result.d_inner;
  Permutation& pi=result.d_perm;
  @< Compute the Lie type |type|, the inner class |inner_class|, and the
     permutation |pi| of the simple roots with respect to standard order for
     |type| @>
  if (r>s)
    @< Add type letters and inner class symbols for the central torus @>
  return result;
}

@ That |M| is an automorphism means that the roots are permuted among
each other, and that after applying the Weyl group action to map simple roots
to simple roots, the result is a diagram automorphism; this is tested by
checking that the Cartan matrix is invariant under the corresponding
permutation of its rows and columns.

@f Delta NULL

@< Set |ww| to the reversed Weyl group element...@>=
{ RootNbrList Delta(s);
  for (weyl::Generator i=0; i<s; ++i)
  { Delta[i]=rd.rootNbr(M*rd.simpleRoot(i));
    if (Delta[i]==rd.numRoots()) // then image not found
      throw std::runtime_error@|
        ("Matrix maps simple root "+str(i)+" to non-root");
  }
  ww = wrt_distinguished(rd,Delta);
  for (weyl::Generator i=0; i<s; ++i)
    if (rd.isSimpleRoot(Delta[i])) // should have made every root simple
      p[i]=rd.simpleRootIndex(Delta[i]);
    else throw std::runtime_error
      ("Matrix does not define a root datum automorphism");
@.Matrix does not define...@>
  for (weyl::Generator i=0; i<s; ++i)
    for (weyl::Generator j=0; j<s; ++j)
      if (rd.cartan(p[i],p[j])!=rd.cartan(i,j)) throw std::runtime_error@|
      ("Matrix does not define a root datum automorphism");
}

@ For each simple factor we look if there are any non-fixed points of the
permutation (if not we have the compact inner class) and if so, whether the
image of that point lies in the same component of the Dynkin diagram. If the
latter is the case we have an unequal rank inner class, which is actually
called split unless the simple factor is of type $D_{2n}$, and if the image
lies in another component of the diagram we have a Complex inner class.

@h "dynkin.h"

@< Compute the Lie type |type|, the inner class... @>=
{ int_Matrix Cartan = rd.cartanMatrix();
  RankFlagsList comp =
    dynkin::components(dynkin::DynkinDiagram(Cartan));

  type.reserve(comp.size()+r-s);
  inner_class.reserve(comp.size()+r-s); // certainly enough
  type = dynkin::Lie_type(Cartan,true,false,pi);
    // no need to check validity of |Cartan|
  assert(type.size()==comp.size());
  size_t offset=0; // accumulated rank of simple factors seen

  for (size_t i=0; i<type.size(); ++i)
  { bool equal_rank=true;
    size_t comp_rank = type[i].rank();
    assert (comp_rank==comp[i].count());
       // and |pi[j]| runs through |comp[i]| in following loop
    for (size_t j=offset; j<offset+comp_rank; ++j) // traverse component
      if (p[pi[j]]!=pi[j]) {@; equal_rank=false; break; }
@)  if (equal_rank) inner_class.push_back('c');
      // identity on this component: compact component
    else if(comp[i].test(p[pi[offset]]))
       // (any) one root stays in |comp[i]|; unequal rank
      inner_class.push_back(type[i].first=='D' and comp_rank%2==0 ? 'u' : 's');
    else
    { inner_class.push_back('C'); // record Complex component
      @< Gather elements of Complex inner class component, adapting the values
         of |type|, |comp| and |pi| @>
      offset += comp_rank;
      ++i; // skip over component |i|, loop will skip component |i+1|
    }

    offset += comp_rank;
  } // |for (i)|
}

@ Complex factors of the inner class involve two simple factors, which
requires some additional care. The corresponding components of the Dynkin
diagram might not be consecutive, in which case we must permute the factors of
|type| to make that true, permute the |comp| subsets to match, and update |pi|
as well. Moreover we wish that, as seen through |pi|, the permutation |p|
interchanges the roots of the two now consecutive factors by a fixed shift in
the index by |comp_rank| (the kind that |lietype::involution| produces).

Due to this predetermined order (relative to |p|) in which the matching factor
is to be arranged, it is not necessary to record the original values of |pi|
on this factor (which were produced by |dynkin::Lie_type| in increasing order)
nor the type of the factor, before overwriting then by other factors being
shifted up: we can afterwards simply deduce the proper values from the
matching factor |i|.

@< Gather elements of Complex inner class...@>=
{ size_t beta = p[pi[offset]];
  @< Find the component~|k| after |i| that contains |beta|, move |comp[k]| to
     |comp[i+1]| while shifting any intermediate components correspondingly
     upwards in |pi|, |type| and |comp| @>
@)
  type[i+1]=type[i]; // duplicate factor |i|
  for (size_t j=offset; j<offset+comp_rank; ++j)
    pi[j+comp_rank]=p[pi[j]];
     // reconstruct matching component in matching order
}


@ When the inner class permutation |p| interchanges a component with another,
we search for that component. Then, as remarked above, we can simply shift up
any intermediate values of |pi|, |type| and |comp| to their new places; only
for the bitset |comp[k]| it is worth while to save the old value and reinsert
it at its moved-down place.

@< Find the component~|k| after |i| that contains |beta|...@>=
{ size_t j, k;
  for (j=offset+comp_rank,k=i+1; k<type.size(); ++k)
    if (comp[k].test(beta))
      break;
    else
      j+=type[k].second;
  if (k==type.size())
    throw std::logic_error("Non matching Complex factor");
@.Non matching Complex factor@>

#ifndef NDEBUG
  assert(type[k]==type[i]); // paired simple types for complex factor
  for (size_t l=1; l<comp_rank; ++l)
    assert(comp[k].test(p[pi[offset+l]]));
        // image by |p| of remainder of |comp[i]| matches |comp[k]|
#endif

  if (k>i+1) // then we need to move component |k| down to |i+1|
  {
    while (j-->offset+comp_rank) // shift up intermediate components in |pi|
      pi[j+comp_rank]=pi[j];

    RankFlags match_comp=comp[k];
      // save component matching |comp[i]| under |p|

    while(k-->i+1) // shift intermediate simple factors
    @/{@; type[k+1]=type[k]; comp[k+1]=comp[k]; }

    comp[i+1]=match_comp; // reinsert mathcing component after |comp[i]|
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

@< Add type letters and inner class symbols for the central torus @>=
{ for (size_t k=0; k<r-s; ++k)
    type.push_back(SimpleLieType('T',1));
  int_Matrix root_lattice
    (rd.beginSimpleRoot(),rd.endSimpleRoot(),r,tags::IteratorTag());
@/CoeffList factor; // values will be unused
  int_Matrix basis =
     matreduc::adapted_basis(root_lattice,factor);
@/WeightInvolution inv =
     basis.inverse().block(s,0,r,r)*M*basis.block(0,s,r,r);
    // involution on quotient by root lattice
  std::pair<size_t,size_t> cl=classify_involution(inv);
@/size_t& compact_rank=cl.first;
  size_t& split_rank=cl.second;
  size_t Complex_rank=(r-s-compact_rank-split_rank)/2;
@/while (compact_rank-->0) inner_class.push_back('c');
  while (Complex_rank-->0) inner_class.push_back('C');
  while (split_rank-->0) inner_class.push_back('s');
}

@*2 Storing the inner class values.
Although abstractly an inner class value is described completely by an object
of type |ComplexReductiveGroup|, we shall need to record additional
information in order to be able to present meaningful names for the real forms
and dual real forms in this inner class. The above analysis of involutions was
necessary in order to obtain such information; it will be recorded in values
of type |realform_io::Interface|.

@< Includes... @>=
#include "realform_io.h"

@~The class |inner_class_value| will be the first built-in type where we
deviate from the previously used scheme of holding an \.{atlas} object with
the main value in a data member |val|. The reason is that the copy constructor
for |ComplexReductiveGroup| is private (and nowhere defined), so
that the straightforward definition of a copy constructor for such a built-in
type would not work, and the copy constructor is necessary for the |clone|
method. (In fact, now that normal manipulation of values involves duplicating
shared pointers rather than of values, there is never a need to copy an
|inner_class_value|, since |get_own<inner_class_value>| is never called;
however the |clone| method is still defined for possible future use.) So
instead, we shall share the \.{atlas} object when duplicating our value, and
maintain a reference count to allow destruction when the last copy disappears.

The reference count needs to be shared of course, and since the links between
the |inner_class_value| and both the \.{atlas} value and the reference count
are indissoluble, we use references for the members |val|, |dual| and
|ref_count|. The first two references are not |const|, since some methods will
as a side effect generate |CartanClass| objects in the inner class, whence
they are technically manipulators rather than accessors.

The main constructor takes a auto-pointer to a |ComplexReductiveGroup| as
argument, as a reminder that the caller gives up ownership of this pointer
that should come from a call to~|new|; this pointer will henceforth be owned
by the |inner_class_value| constructed, in shared ownership with any values
later cloned from it: the last one of them to be destroyed will call |delete|
for the pointer. The remaining argument is a |Layout| that must have been
computed by |check_involution| above, in order to ensure its validity.

Occasionally we shall need to refer to the dual inner class (for the dual
group); since the construction of an instance takes some work, we do not wish
to repeat that every time the dual is needed, so we create the dual
|ComplexReductiveGroup| value upon construction of the
|inner_class_value| and store it in the |dual| field where it will be
available as needed.

Unlike for other value types, the copy constructor is public here. Thus a
class depending on our |val| being valid can simply store a copy of our
|inner_class_value|; this does not cost much, and the reference counting
mechanism will then ensure that |val| remains valid while the object of that
containing class exists.

@< Type definitions @>=
struct inner_class_value : public value_base
{ ComplexReductiveGroup& val;
  ComplexReductiveGroup& dual;
  size_t& ref_count;
@)
  lietype::Layout type_info;
  const realform_io::Interface interface,dual_interface;
@)
  inner_class_value(std::auto_ptr<ComplexReductiveGroup> G
   ,@| const lietype::Layout& lo); // main constructor
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
typedef std::tr1::shared_ptr<inner_class_value> shared_inner_class;

@ Here are the copy constructor and the destructor.
@< Function def...@>=
inner_class_value::inner_class_value(const inner_class_value& v)
: val(v.val), dual(v.dual), ref_count(v.ref_count)
, type_info(v.type_info)
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
an auto-pointer until the construction succeeds.

@< Function def...@>=
inner_class_value::inner_class_value
  (std::auto_ptr<ComplexReductiveGroup> g,
   const lietype::Layout& lo)
@/: val(*g)
, dual(*new ComplexReductiveGroup(*g,tags::DualTag()))
@/, ref_count(*new size_t(1))
@/, type_info(lo), interface(*g,lo), dual_interface(*g,lo,tags::DualTag())
 {@; g.release(); } // now that we own |g|, release the auto-pointer

@ We allow construction of a dual |inner_class_value|. Since it can share the
two fields referring to objects of type |ComplexReductiveGroup|
in the opposite order, we can consider it as a member of the same reference
counted family, and share the |ref_count| field. This means this constructor
is more like the copy constructor than like the main constructor, and in
particular the reference count is increased. The dual of the dual inner class
value will behave exactly like a copy of the original inner class.

@< Function def...@>=
inner_class_value::inner_class_value(const inner_class_value& v,tags::DualTag)
@/: val(v.dual), dual(v.val)
, ref_count(v.ref_count)
@/,type_info(v.type_info)
, interface(v.dual_interface), dual_interface(v.interface)
{ ++ref_count;
@/type_info.d_type=lietype::dual_type(type_info.d_type);
  type_info.d_inner=lietype::dual_type(type_info.d_inner,type_info.d_type);
@/// |d_perm| remains identical to |v.d_perm|
}


@ One of the most practical informations about a |ComplexReductiveGroup|,
which is available directly after its construction, is the number of real
forms in the inner class defined by it; we print this information when a
|inner_class_value| is printed.

@< Function def...@>=
void inner_class_value::print(std::ostream& out) const
{ out << "Complex reductive group of type " << type_info.d_type @|
      << ", with involution defining\n"
         "inner class of type '" << type_info.d_inner @|
      << "', with " << val.numRealForms() @| << " real "
      << (val.numRealForms()==1 ? "form" : "forms") @| << " and "
      << val.numDualRealForms() @| << " dual real "
      << (val.numDualRealForms()==1 ? "form" : "forms");
}

@ Our wrapper function builds a complex reductive group with an involution,
testing its validity. The Weyl word |ww| that was needed in the test to make
the involution into one of the based root datum must be applied to the matrix,
since the test does not actually modify its matrix argument. Then the root
datum is passed to a |ComplexReductiveGroup| constructor that the library
provides specifically for this purpose, and which makes a copy of the root
datum; the \.{atlas} program instead uses a constructor using a |PreRootDatum|
that constructs the |RootDatum| directly into the |ComplexReductiveGroup|, but
using that constructor here would be cumbersome and even less efficient then
copying the existing root datum. Another wrapper |twisted_involution_wrapper|
is similar, but also returns a second value, which is the twisted involution
corresponding in the inner class to the given involution matrix.

@< Local function def...@>=
void fix_involution_wrapper(expression_base::level l)
{ LatticeMatrix M(get<matrix_value>()->val); // safe use of temporary
  shared_root_datum rd(get<root_datum_value>());
  WeylWord ww;
  lietype::Layout lo = check_involution(M,rd->val,ww);
  if (l==expression_base::no_value)
    return;
@)
  for (unsigned int i=0; i<ww.size(); ++i) // apply elements in generation order
    rd->val.simple_reflect(ww[i],M);
  std::auto_ptr<ComplexReductiveGroup> G(new ComplexReductiveGroup(rd->val,M));
  push_value(new inner_class_value(G,lo));
}

void twisted_involution_wrapper(expression_base::level l)
{ LatticeMatrix M(get<matrix_value>()->val);
  shared_root_datum rd(get<root_datum_value>());
  WeylWord ww;
  lietype::Layout lo = check_involution(M,rd->val,ww);
  if (l==expression_base::no_value)
    return;
@)
  for (unsigned int i=ww.size(); i-->0;)
    rd->val.simple_reflect(ww[i],M);
  std::auto_ptr<ComplexReductiveGroup> G(new ComplexReductiveGroup(rd->val,M));
  push_value(new inner_class_value(G,lo));
  push_value(new vector_value(std::vector<int>(ww.begin(),ww.end())));
  if (l==expression_base::single_value)
    wrap_tuple(2);
}

@ To simulate the functioning of the \.{atlas} program, the function
|set_type| takes as argument a Lie type, a list of kernel generators, and a
string describing the inner class. The evaluation of the call
|set_type(lt,gen,ict)| computes ${\it basis}={\it quotient\_basis(lt,gen)}$,
${\it rd}={\it root\_datum (lt,basis)}$, ${\it M}={\it
based\_involution(lt,basis,ict)}$, and then returns the same value as would
|fix_involution(rd,M)|. However we avoid actually calling the function
|fix_involution_wrapper|, which would reconstruct |lt| and |ict| from |rd| and
|M|, while performing tests that are useless given the way $M$ was computed;
rather we store |lt| and |ict| directly in a |Layout| to be stored in the
|inner_class_value|. This follows most closely \.{altas} behaviour, and avoids
surprises (however inner class letters do change as usual to synonyms when
passing through |transform_inner_class_type|).

@< Local function def...@>=
void set_type_wrapper(expression_base::level l)
{ shared_string ict = get<string_value>();
    // and leave generators |gen| and type |lt|
  shared_value lt = *(execution_stack.end()-2);
  LieType& type=force<Lie_type_value>(lt.get())->val;
  lietype::Layout lo(type,transform_inner_class_type(ict->val.c_str(),type));
@)
  quotient_basis_wrapper(expression_base::single_value); @+
  shared_value basis = pop_value();
@)
  push_value(lt); push_value(basis);
  root_datum_wrapper(expression_base::single_value);
  shared_root_datum rd = get<root_datum_value>();
@)
  push_value(lt); push_value(basis);
  push_value(ict);
  based_involution_wrapper(expression_base::single_value);
  shared_matrix M = get<matrix_value>();
  if (l==expression_base::no_value)
    return; // bow out now all possible errors are passed
@)
  std::auto_ptr<ComplexReductiveGroup>@|
    G(new ComplexReductiveGroup(rd->val,M->val));
  push_value(new inner_class_value(G,lo));
}

@ It can be awkward to use |set_type| which wants to construct the root datum
in its own way, for instance if one already has a root datum as produced by
special functions like $GL$. The meaning of an inner class string, in the
sense of defining an involution of the based root datum, is almost
unambiguous, so we introduce a function |set_inner_class| that from a root
datum and an inner class specification produces an inner class for that datum.
A small ambiguity can be present for torus factors, in that there is no
obvious basis of the part of the root lattice orthogonal to the roots on which
to interpret the inner class specification for torus factors. We shall use the
basis implicitly defined by the radical generators in the root datum; if that
should fail to give the desired involution, the use can always take recourse
to |fix_involution| to explicitly specify one.

The wrapper for |set_inner_class| below expresses the permutation matrix~$P$
produced by |lietype::involution| on the basis of the lattice of the root
datum~|rd|. That basis is already the canonical basis of $\Z^n$ used in~|rd|,
so the difficulty is to find the basis~$b$ on which $P$ is expressed. If~|rd|
was entered by the user as in the \.{atlas} program or in
|root_datum_wrapper|, this was the canonical basis used until the
|PreRootDatum| constructor expressed everything on the provided basis of a
sub-lattice; we wish to recover the (integral) matrix~$M$ that specified the
sub-lattice, since it expresses our root datum lattice on~$b$. Unfortunately
the |PreRootDatum| constructor only expresses roots and coroots on this basis
and its dual, so we cannot recover~$M$ reliably from the root datum (for
instance is there are no roots). In fact we shall use the information stored
in the coroots, whose coordinates give the entries of~$M$, except for its rows
of coordinates on the torus factors (cf.\ |prerootdata::corootBasis|). We
complement this with coordinates of the radical generators stored in the root
datum to find a matrix~$M$ that \emph{could have} been used to obtain the root
datum: using it with the Lie type of~|rd| certainly would give the same
coroots, while the roots, being determined by their values on coroots and the
radical, would be identical as well; all other data in |rd| is deduced from
these.

The |Layout| structure records a Lie type~|lt|, a permutation~|pi| with
respect to the standard diagram ordering, and an inner class string. The call
to |dynkin::Lie_type| deduces |lt| and~|pi|, and we record the given inner
class type~|ict|. Then we find |M=transpose_mat(coroot_radical(rd))|, call
|lietype::involution| to produce an involution for |lt| and~|pi|, express it
as~|inv| on the basis~$M$, and finally call |fix_involution(rd,inv)| to
produce the result. The |on_basis| method may throw a |std::runtime_error|
when the involution does not stabilise the sub-lattice; we catch this error
and re-throw with a more explicit error indication.

@< Local function def...@>=
void set_inner_class_wrapper(expression_base::level l)
{ shared_string ict(get<string_value>());
  shared_root_datum rdv(get<root_datum_value>());
  const RootDatum& rd=rdv->val;
@)
  lietype::Layout lo;
  lo.d_type = dynkin::Lie_type(rd.cartanMatrix(),true,false,lo.d_perm);
   // get (permuted) type
  for (size_t i=rd.semisimpleRank(); i<rd.rank(); ++i) // if not semisimple
  { lo.d_type.push_back(SimpleLieType('T',1)); // add a torus factor
    lo.d_perm.push_back(i);
      // and a fixed point of permutation, needed by |lietype::involution|
  }
  lo.d_inner=transform_inner_class_type(ict->val.c_str(),lo.d_type);
@)
  push_value(rdv);
  coroot_radical_wrapper(expression_base::single_value);
  transpose_mat_wrapper(expression_base::single_value);
  shared_matrix M(get<matrix_value>());
  size_t r=lo.d_type.rank();
  assert(M->val.numRows()==r and M->val.numRows()==r);
  ndebug_use(r);
@)
  push_value(rdv);
  try
  {@; push_value(new matrix_value
       (lietype::involution(lo).on_basis(M->val.columns())));
  }
  catch (std::runtime_error&) // relabel inexact division error
  {@; throw std::runtime_error @|
    ("Inner class is not compatible with root datum");
@.Inner class is not compatible...@>
  }
  fix_involution_wrapper(l);
}

@*2 Functions operating on inner classes.
%
Here are our first functions that access a operate on values of type
|inner_class_value|; they recover the ingredients that were used in the
construction, and the final one constructs the dual inner class.

@< Local function def...@>=
void inner_class_equality_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  shared_inner_class H(get<inner_class_value>());
  if (l!=expression_base::no_value)
    push_value(new bool_value(&G->val==&H->val));
}
@)
void distinguished_involution_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l!=expression_base::no_value)
    push_value(new matrix_value(G->val.distinguished()));
}
@)
void root_datum_of_inner_class_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l!=expression_base::no_value)
    push_value(new root_datum_value(G->val.rootDatum()));
}
@)
void inner_class_to_root_datum_coercion()
{@; root_datum_of_inner_class_wrapper(expression_base::single_value); }

void dual_inner_class_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l!=expression_base::no_value)
    push_value(new inner_class_value(*G,tags::DualTag()));
}

@ More interestingly, let us extract the list of names of the real forms.
This uses the interface fields stored in the value. Since they exist for both
the group itself and for the dual group, we define an auxiliary function that
produces the list, and then use it twice.

@< Local function def...@>=
void push_name_list(const realform_io::Interface& interface)
{ row_ptr result(new row_value(0));
  for (size_t i=0; i<interface.numRealForms(); ++i)
    result->val.push_back
      (shared_value(new string_value(interface.typeName(i))));
  push_value(result);
}
@)
void form_names_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l!=expression_base::no_value)
    push_name_list(G->interface);
}
@)
void dual_form_names_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l!=expression_base::no_value)
    push_name_list(G->dual_interface);
}

@ We provide functions for counting (dual) real forms and Cartan classes,
although the former two could be computed as the lengths of the list of (dual)
form names.

@< Local function def...@>=
void n_real_forms_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l!=expression_base::no_value)
    push_value(new int_value(G->val.numRealForms()));
}
@)
void n_dual_real_forms_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l!=expression_base::no_value)
    push_value(new int_value(G->val.numDualRealForms()));
}
@)
void n_Cartan_classes_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l!=expression_base::no_value)
    push_value(new int_value(G->val.numCartanClasses()));
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
void block_sizes_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l==expression_base::no_value)
    return;
  matrix_ptr M(new matrix_value @|
    (int_Matrix(G->val.numRealForms()
                                ,G->val.numDualRealForms())
    ));
  for (size_t i = 0; i < M->val.numRows(); ++i)
    for (size_t j = 0; j < M->val.numColumns(); ++j)
      M->val(i,j) =
      G->val.block_size(G->interface.in(i),G->dual_interface.in(j));
  push_value(M);
}

@ Next we shall provide a function that displays the occurrence of Cartan
classes for various real forms. The rows will be indexed by real forms, and
the columns by Cartan classes (note the alliteration).

@< Local function def...@>=
void occurrence_matrix_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l==expression_base::no_value)
    return;
  size_t nr=G->val.numRealForms();
  size_t nc=G->val.numCartanClasses();
  matrix_ptr M(new matrix_value(int_Matrix(nr,nc)));
  for (size_t i=0; i<nr; ++i)
  { BitMap b=G->val.Cartan_set(G->interface.in(i));
    for (size_t j=0; j<nc; ++j)
      M->val(i,j)= b.isMember(j) ? 1 : 0;
  }
  push_value(M);
}

@ We do the same for dual real forms. Note that we had to introduce the method
|dualCartanSet| for |ComplexReductiveGroup| in order to be able
to write this function.

@< Local function def...@>=
void dual_occurrence_matrix_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l==expression_base::no_value)
    return;
  size_t nr=G->val.numDualRealForms();
  size_t nc=G->val.numCartanClasses();
  matrix_ptr M(new matrix_value(int_Matrix(nr,nc)));
  for (size_t i=0; i<nr; ++i)
  { BitMap b=G->val.dual_Cartan_set(G->dual_interface.in(i));
    for (size_t j=0; j<nc; ++j)
      M->val(i,j)= b.isMember(j) ? 1 : 0;
  }
  push_value(M);
}

@ Finally we install everything related to inner classes.
@< Install wrapper functions @>=

install_function(classify_wrapper,@|"classify_involution"
                ,"(mat->int,int,int)");
install_function(fix_involution_wrapper,@|"inner_class"
                ,"(RootDatum,mat->InnerClass)");
install_function(twisted_involution_wrapper,@|"twisted_involution"
                ,"(RootDatum,mat->InnerClass,vec)");
install_function(set_type_wrapper,@|"inner_class"
                ,"(LieType,[ratvec],string->InnerClass)");
install_function(set_inner_class_wrapper,@|"inner_class"
                ,"(RootDatum,string->InnerClass)");
install_function(inner_class_equality_wrapper,@|"="
                ,"(InnerClass,InnerClass->bool)");
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
install_function(block_sizes_wrapper,@|"block_sizes"
                ,"(InnerClass->mat)");
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
#include "kgb.h"

@*2 Class definition.
The layout of this type of value is different from what we have seen before.
An Atlas object of class |RealReductiveGroup| is dependent upon another Atlas
object to which it stores a pointer, which is of type |ComplexReductiveGroup|,
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

We shall later derive from this structure, without adding data members, a
structure to record dual real forms, which are constructed in the context of
the same inner class, but where the |val| field is constructed for the dual
group (and its dual inner class). In order to accommodate for that we provide
an additional protected constructor that will be defined later; this also
explains why the copy constructor is protected rather than private.

@< Type definitions @>=
struct real_form_value : public value_base
{ const inner_class_value parent;
  RealReductiveGroup val;
  RealFormNbr form_index; // internal number associated to real form in |parent|
@)
  real_form_value(const inner_class_value& p,RealFormNbr f) @/
  : parent(p), val(p.val,f),form_index( f), rt_p(NULL) @+{}
@)
  virtual void print(std::ostream& out) const;
  real_form_value* clone() const @+
    {@; return new real_form_value(*this); }
  static const char* name() @+{@; return "real form"; }
  const KGB& kgb () @+{@; return val.kgb(); }
   // generate and return $K\backslash G/B$ set
  const Rep_context& rc();
  Rep_table& rt();
  ~real_form_value() @+{@; delete rt_p; }
protected:
  real_form_value(const inner_class_value&,RealFormNbr,tags::DualTag);
  real_form_value(const real_form_value& v)
  : parent(v.parent), val(v.val), rt_p(v.rt_p) @+{}
private:
  Rep_table* rt_p;
    // owned pointer, initially |NULL|, assigned at most once
};
@)
typedef std::auto_ptr<real_form_value> real_form_ptr;
typedef std::tr1::shared_ptr<real_form_value> shared_real_form;

@ The methods |rc| and |rt| ensure a |Rep_table| value is constructed at
|*rt_p|, and returns a reference. The value so obtained will serve to
manipulate parameters for standard modules, for which we shall define a
built-in type below. Storing the value here ensures that it will be shared
between different parameters, and that it will live as long as those parameter
values do. The value itself does not take much space, but constructing it
implicitly calls the |val.kgb| method, so we avoid doing this until there is a
concrete need.

@< Function def...@>=
  const Rep_context& real_form_value::rc()
    {@; return *(rt_p==NULL ? rt_p=new Rep_table(val) : rt_p); }
  Rep_table& real_form_value::rt()
    {@; return *(rt_p==NULL ? rt_p=new Rep_table(val) : rt_p); }

@ When printing a real form, we give the name by which it is known in the
parent inner class, and provide some information about its connectedness.
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
      << parent.interface.typeName(parent.interface.out(val.realForm())) @|
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
{ shared_int i(get<int_value>());
  shared_inner_class G(get<inner_class_value>());
  if (size_t(i->val)>=G->val.numRealForms())
    throw std::runtime_error ("Illegal real form number: "+str(i->val));
@.Illegal real form number@>
  if (l!=expression_base::no_value)
    push_value(new real_form_value(*G,G->interface.in(i->val)));
}
@)
void form_number_wrapper(expression_base::level l)
{ shared_real_form rf(get<real_form_value>());
  if (l!=expression_base::no_value)
    push_value(new int_value(rf->parent.interface.out(rf->val.realForm())));
}
@)
void quasisplit_form_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l!=expression_base::no_value)
    push_value(new real_form_value(*G,G->val.quasisplit()));
}

@*2 Functions operating on real reductive groups.
%
From a real reductive group we can go back to its inner class
@< Local function def...@>=
void inner_class_of_real_form_wrapper(expression_base::level l)
{ shared_real_form rf(get<real_form_value>());
  if (l!=expression_base::no_value)
    push_value(new inner_class_value(rf->parent));
}

void real_form_to_inner_class_coercion()
{@; inner_class_of_real_form_wrapper(expression_base::single_value); }

void real_form_to_root_datum_coercion()
{ shared_real_form rf(get<real_form_value>());
  push_value(new root_datum_value(rf->parent.val.rootDatum()));
}

@ Here is a function that gives information about the dual component group
of a real reductive group: it returns the rank~$r$ such that the component
group is isomorphic to $(\Z/2\Z)^r$.

@< Local function def...@>=
void components_rank_wrapper(expression_base::level l)
{ shared_real_form R(get<real_form_value>());
  if (l!=expression_base::no_value)
    push_value(new int_value(R->val.dualComponentReps().size()));
}

@ And here is one that counts the number of Cartan classes for the real form.

@< Local function def...@>=
void count_Cartans_wrapper(expression_base::level l)
{ shared_real_form rf(get<real_form_value>());
  if (l!=expression_base::no_value)
    push_value(new int_value(rf->val.numCartan()));
}

@ The size of the finite set $K\backslash G/B$ can be determined from the real
form, (the necessary |CartanClass| objects will be automatically generated
when doing so).

@< Local function def...@>=
void KGB_size_wrapper(expression_base::level l)
{ shared_real_form rf(get<real_form_value>());
  if (l!=expression_base::no_value)
    push_value(new int_value(rf->val.KGB_size()));
}

@ Here is a somewhat technical function that will facilitate working ``in
coordinates'' with KGB elements for this real form. It returns a rational
vector that determines the base grading for the real form, which is an offset
that should be added to |torus_bits| values. It is defined so as to be zero
for quasisplit real forms, so in order to compute gradings, the
imaginary-$\rho$ value (dependent on the involution) should be added to it.

@< Local function def...@>=
void base_grading_vector_wrapper(expression_base::level l)
{ shared_real_form rf(get<real_form_value>());
  if (l!=expression_base::no_value)
  { const KGB& kgb=rf->kgb();
    RatCoweight t = kgb.base_grading_vector();
    push_value(new rational_vector_value(t));
  }
}

@ There is a partial ordering on the Cartan classes defined for a real form. A
matrix for this partial ordering is computed by the function
|Cartan_order_matrix|, which more or less replaces the \.{corder} command in
\.{atlas}.

@h "poset.h"
@< Local function def...@>=
void Cartan_order_matrix_wrapper(expression_base::level l)
{ shared_real_form rf(get<real_form_value>());
  if (l==expression_base::no_value)
    return;
  size_t n=rf->val.numCartan();
  matrix_ptr M(new matrix_value(int_Matrix(n,n,0)));
  const poset::Poset& p = rf->val.complexGroup().Cartan_ordering();
  for (size_t i=0; i<n; ++i)
    for (size_t j=i; j<n; ++j)
      if (p.lesseq(i,j)) M->val(i,j)=1;

  push_value(M);
}

@*2 Dual real forms.
Although they could be considered as real forms for the dual group and dual
inner class, it will be convenient to be able to construct dual real forms in
the context of the |inner_class_value| itself. The resulting objects will have
a different type than ordinary real forms, but they will use the same structure
layout. The interpretation of the fields is as follows for dual real forms:
the |parent| field contains a copy of the original |inner_class_value|, but
the |val| field contains a |RealReductiveGroup| object
constructed for the dual inner class, so that its characteristics will be
correct.

@< Function def...@>=

real_form_value::real_form_value(const inner_class_value& p,RealFormNbr f
				,tags::DualTag)
: parent(p), val(p.dual,f), rt_p(NULL)
@+{}

@ In order to be able to use the same layout for dual real forms but to
nevertheless make a distinction (for instance in printing), we must derive a
new type from |real_form_value|.

@< Type definitions @>=
struct dual_real_form_value : public real_form_value
{ dual_real_form_value(const inner_class_value& p,RealFormNbr f)
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
typedef std::tr1::shared_ptr<dual_real_form_value> shared_dual_real_form;

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
void dual_real_form_wrapper(expression_base::level l)
{ shared_int i(get<int_value>());
  shared_inner_class G(get<inner_class_value>());
  if (size_t(i->val)>=G->val.numDualRealForms())
    throw std::runtime_error ("Illegal dual real form number: "+str(i->val));
@.Illegal dual real form number@>
  if (l!=expression_base::no_value)
    push_value(new dual_real_form_value(*G,G->dual_interface.in(i->val)));
}
@)
void dual_quasisplit_form_wrapper(expression_base::level l)
{ shared_inner_class G(get<inner_class_value>());
  if (l!=expression_base::no_value)
    push_value(new dual_real_form_value(*G,G->dual.quasisplit()));
}

@ Rather than provide all functions for real forms for dual real forms, we
provide functions that convert a dual real form to a real form, which of
course will be associated to the dual inner class, and back.

@< Local function def...@>=
void real_form_from_dual_wrapper(expression_base::level l)
{ shared_dual_real_form d(get<dual_real_form_value>());
  if (l!=expression_base::no_value)
    push_value(new real_form_value
                 (inner_class_value(d->parent,tags::DualTag())
                 ,d->val.realForm()));
}

void dual_real_form_from_wrapper(expression_base::level l)
{ shared_real_form rf = get<real_form_value>();
  if (l!=expression_base::no_value)
    push_value(new dual_real_form_value
                 (inner_class_value(rf->parent,tags::DualTag())
                 ,rf->val.realForm()));
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
install_function(Cartan_order_matrix_wrapper,@|"Cartan_order"
					    ,"(RealForm->mat)");
install_function(dual_real_form_wrapper,@|"dual_real_form"
				       ,"(InnerClass,int->DualRealForm)");
install_function(dual_quasisplit_form_wrapper,@|"dual_quasisplit_form"
		,"(InnerClass->DualRealForm)");
install_function(real_form_from_dual_wrapper,@|"real_form"
				  ,"(DualRealForm->RealForm)");
install_function(dual_real_form_from_wrapper,@|"dual_real_form"
				  ,"(RealForm->DualRealForm)");

@*1 A type for Cartan classes.
%
Another type of value associated to inner classes are Cartan classes, which
describe equivalence classes (for ``stable conjugacy'') of Cartan subgroups of
real reductive groups in an inner class. The Atlas software associates a fixed
set of Cartan classes to each inner class, and records for each real form the
subset of those Cartan classes that occur for the real form. Since versions
0.3.5 of the software, the Cartan classes are identified and numbered upon
construction of a |ComplexReductiveGroup| object, but |CartanClass| objects
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
typedef std::auto_ptr<Cartan_class_value> Cartan_class_ptr;
typedef std::tr1::shared_ptr<Cartan_class_value> shared_Cartan_class;

@ In the constructor we used to check that the Cartan class with the given
number currently exists, but now the |ComplexReductiveGroup::cartan| method
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
{ shared_int i(get<int_value>());
  shared_inner_class ic(get<inner_class_value>());
  if (size_t(i->val)>=ic->val.numCartanClasses())
    throw std::runtime_error ("Illegal Cartan class number: "+str(i->val)
@.Illegal Cartan class number@>
    +", this inner class only has "+str(ic->val.numCartanClasses())
    +" of them");
  if (l!=expression_base::no_value)
    push_value(new Cartan_class_value(*ic,i->val));
}

@ Alternatively (and this used to be the only way) one can provide a
|real_form_value| together with a valid index into \emph{its} list of Cartan
classes. We translate this number into an index into the list for its
containing inner class, and then get the Cartan class from there.

@< Local function def...@>=
void rf_Cartan_class_wrapper(expression_base::level l)
{ shared_int i(get<int_value>());
  shared_real_form rf(get<real_form_value>());
  if (size_t(i->val)>=rf->val.numCartan())
    throw std::runtime_error ("Illegal Cartan class number: "+str(i->val)
@.Illegal Cartan class number@>
    +", this real form only has "+str(rf->val.numCartan())+" of them");
  BitMap cs=rf->val.Cartan_set();
  if (l!=expression_base::no_value)
    push_value(new Cartan_class_value(rf->parent,cs.n_th(i->val)));
}

@ Like the quasisplit real form for inner classes, there is a particular
Cartan class for each real form, the ``most split Cartan''; it is of
particular importance for knowing which dual real forms can be associated to
this real form. It is (probably) the last Cartan class in the list for the
real form, but we have a direct access to it via the |mostSplit| method for
|RealReductiveGroup|.

@< Local function def...@>=
void most_split_Cartan_wrapper(expression_base::level l)
{ shared_real_form rf(get<real_form_value>());
  if (l!=expression_base::no_value)
    push_value(new Cartan_class_value(rf->parent,rf->val.mostSplit()));
}

@ It is useful to be able to compute a real form or Cartan class based on
other information than their enumeration within an inner class. The \.{realex}
function |Cartan_class_real_form| takes an inner class, a matrix giving an
involution, and a rational co-weight describing the grading of the imaginary
roots, or more precisely the offset with respect to the grading that makes all
simple-imaginary roots non-compact; it returns a pair consisting of the Cartan
class for the involution and the real form associated there to the grading.

@:Cartan_class_real_form@>

@< Local function def...@>=
TwistedInvolution twisted_from_involution
  (const ComplexReductiveGroup& G, const WeightInvolution theta)
{ const RootDatum& rd = G.rootDatum();
  unsigned int ssr =  rd.semisimpleRank();
  RootNbrList Delta(ssr);
  for (weyl::Generator i=0; i<ssr; ++i)
  { Delta[i]=rd.rootNbr(theta*rd.simpleRoot(i));
    if (Delta[i]==rd.numRoots()) // then image not found
      throw std::runtime_error@|
        ("Matrix maps simple root "+str(i)+" to non-root");
  }
  WeylWord ww = wrt_distinguished(rd,Delta);
  for (weyl::Generator i=0; i<ssr; ++i)
    if (Delta[i]!=G.twisted_root(rd.simpleRootNbr(i)))
        // |Delta| should match distinguished involution
      throw std::runtime_error@|
        ("Matrix does not define a root datum automorphism");
  return G.weylGroup().element(ww);
}

void Cartan_class_real_form_wrapper(expression_base::level l)
{ shared_rational_vector grading_shift = get<rational_vector_value>();
  shared_matrix theta = get<matrix_value>();
  shared_inner_class G(get<inner_class_value>());
  { Coweight num(grading_shift->val.numerator().begin(),
                 grading_shift->val.numerator().end());
    if (theta->val.right_mult(num)!=num)
      throw std::runtime_error ("Grading coweight not fixed by involution");
@.Grading coweight not fixed@>
  }

  TwistedInvolution tw = twisted_from_involution(G->val,theta->val);

  CartanNbr cn; WeylWord ww; cartanclass::AdjointFiberElt rep;
  RealFormNbr rf = real_form_of(G->val,tw,grading_shift->val,cn,ww,rep);
  if (l==expression_base::no_value)
     return;
  push_value(new Cartan_class_value(*G,cn));
  push_value(new real_form_value(*G,rf));
  if (l==expression_base::single_value)
    wrap_tuple(2);
}



@*2 Functions operating on Cartan classes.
%
We start with a fundamental attribute of Cartan classes: the associated
(distinguished) involution, in the form of a matrix.

@< Local function def...@>=
void Cartan_involution_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  if (l==expression_base::no_value)
    return;
  push_value(new matrix_value(cc->val.involution()));
}


@ This function and the following provide the functionality of the \.{atlas}
command \.{cartan}. They are based on |cartan_io::printCartanClass|, but
rewritten to take into account the fact that we do not know about |Interface|
objects for complex groups, and such that a usable value is returned. We omit
in our function {\it Cartan\_info} the data printed in |printCartanClass| in
the final call to |cartan_io::printFiber| (namely all the real forms for which
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

  push_value(new int_value(cc->val.fiber().torus().compactRank()));
  push_value(new int_value(cc->val.fiber().torus().complexRank()));
  push_value(new int_value(cc->val.fiber().torus().splitRank()));
  wrap_tuple(3);

  const weyl::TwistedInvolution& tw =
    cc->parent.val.involution_of_Cartan(cc->number);
  WeylWord ww = cc->parent.val.weylGroup().word(tw);

  std::vector<int> v(ww.begin(),ww.end());
  push_value(new vector_value(v));

  push_value(new int_value(cc->val.orbitSize()));
  push_value(new int_value(cc->val.fiber().fiberSize()));
  wrap_tuple(2);

  const RootSystem& rs=cc->parent.val.rootDatum();

@)// print types of imaginary and real root systems and of Complex factor
  push_value(new Lie_type_value(rs.Lie_type(cc->val.simpleImaginary())));
  push_value(new Lie_type_value(rs.Lie_type(cc->val.simpleReal())));
  push_value(new Lie_type_value(rs.Lie_type(cc->val.simpleComplex())));
  wrap_tuple(3);
@)
  if (l==expression_base::single_value)
    wrap_tuple(4);
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
void real_forms_of_Cartan_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  const inner_class_value& ic=cc->parent;
  if (l==expression_base::no_value)
    return;
  row_ptr result @| (new row_value(cc->val.numRealForms()));
  for (size_t i=0,k=0; i<ic.val.numRealForms(); ++i)
  { BitMap b(ic.val.Cartan_set(ic.interface.in(i)));
    if (b.isMember(cc->number))
      result->val[k++] =
	shared_value(new real_form_value(ic,ic.interface.in(i)));
  }
  push_value(result);
}
@)
void dual_real_forms_of_Cartan_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  const inner_class_value& ic=cc->parent;
  if (l==expression_base::no_value)
    return;
  row_ptr result @| (new row_value(cc->val.numDualRealForms()));
  for (size_t i=0,k=0; i<ic.val.numDualRealForms(); ++i)
  { BitMap b(ic.val.dual_Cartan_set(ic.dual_interface.in(i)));
    if (b.isMember(cc->number))
      result->val[k++] =
	shared_value(new dual_real_form_value(ic,ic.dual_interface.in(i)));
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

@< Local function def...@>=
void fiber_partition_wrapper(expression_base::level l)
{ shared_real_form rf(get<real_form_value>());
  shared_Cartan_class cc(get<Cartan_class_value>());
  if (&rf->parent.val!=&cc->parent.val)
    throw std::runtime_error
    ("Inner class mismatch between real form and Cartan class");
@.Inner class mismatch...@>
  BitMap b(cc->parent.val.Cartan_set(rf->val.realForm()));
  if (!b.isMember(cc->number))
    throw std::runtime_error
    ("Cartan class not defined for this real form");
@.Cartan class not defined@>
  if (l==expression_base::no_value)
    return;
@)
  const Partition& pi = cc->val.fiber().weakReal();
  const RealFormNbrList rf_nr=
     cc->parent.val.realFormLabels(cc->number);
     // translate part number of |pi| to real form
  row_ptr result (new row_value(0)); // cannot predict exact size here
  for (size_t i=0; i<pi.size(); ++i)
    if (rf_nr[pi.class_of(i)] == rf->val.realForm())
      result->val.push_back(shared_value(new int_value(i)));
  push_value(result);
}

@ The function |square_classes| returns the set of real forms associated to a
Cartan class, but partitioned according to their square class, determined by
the square of any strong involution representing the real form.
@< Local function def...@>=
void square_classes_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  const realform_io::Interface rfi = cc->parent.interface;
  const RealFormNbrList& rfl = cc->parent.val.realFormLabels(cc->number);
  if (l==expression_base::no_value)
    return;
  size_t n_sq_classes = cc->val.numRealFormClasses();
  row_ptr result @| (new row_value(n_sq_classes));
  for (cartanclass::square_class csc=0; csc<n_sq_classes; ++csc)
  { const Partition& pi = cc->val.fiber_partition(csc);
    row_ptr part @| (new row_value(pi.classCount()));
    for (unsigned long c=0; c<pi.classCount(); ++c)
       part->val[c] =
          shared_value(new int_value(rfi.out(rfl[cc->val.toWeakReal(c,csc)])));
    result->val[csc] = part;
  }
  push_value(result);
}

@ The function |print_gradings| gives on a per-real-form basis the
functionality of the Atlas command \.{gradings} that is implemented by
|complexredgp_io::printGradings| and |cartan_io::printGradings|. It therefore
takes, like |fiber_partition|, a Cartan class and a real form as parameter. Its
output consist of a list of $\Z/2\Z$-gradings of each of the fiber group
elements in the part corresponding to the real form, where each grading is a
sequence of bits corresponding to the simple imaginary roots.

@f sigma NULL

@< Local function def...@>=
void print_gradings_wrapper(expression_base::level l)
{ shared_real_form rf(get<real_form_value>());
@/shared_Cartan_class cc(get<Cartan_class_value>());
  if (&rf->parent.val!=&cc->parent.val)
    throw std::runtime_error
    ("Inner class mismatch between real form and Cartan class");
@.Inner class mismatch...@>
  BitMap b(cc->parent.val.Cartan_set(rf->val.realForm()));
  if (!b.isMember(cc->number))
    throw std::runtime_error ("Cartan class not defined for this real form");
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
    wrap_tuple(0);
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
    { os << ( first ? first=false,'[' : ',');
      Grading gr=cc->val.fiber().grading(i);
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
install_function(Cartan_class_real_form_wrapper,@|"Cartan_class_real_form"
		,"(InnerClass,mat,ratvec->CartanClass,RealForm)");
install_function(Cartan_involution_wrapper,@|"involution","(CartanClass->mat)");
install_function(Cartan_info_wrapper,@|"Cartan_info"
		,"(CartanClass->(int,int,int),"
                 "vec,(int,int),(LieType,LieType,LieType))");
install_function(real_forms_of_Cartan_wrapper,@|"real_forms"
		,"(CartanClass->[RealForm])");
install_function(dual_real_forms_of_Cartan_wrapper,@|"dual_real_forms"
		,"(CartanClass->[DualRealForm])");
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
and a $\Z/2\Z$-grading of the imaginary root subsystem. It does not seem
necessary to introduce a type for the set $K\backslash G/B$ (any function that
needs it can take the corresponding real form as argument), but we do provide
one for elements of that set (if desired the user can collect such elements in
an array). These elements contain a pointer |rf | to the |real_form_value| they
were generated from, and thereby to the associated set $K\backslash G/B$,
which allows operations relevant to the KGB element to be defined. Since this
is a shared pointer, the |real_form_value| pointed to is guaranteed to exist
as long as this |KGB_elt_value| does.

@< Type definitions @>=
struct KGB_elt_value : public value_base
{ shared_real_form rf;
  KGBElt val;
@)
  KGB_elt_value(const shared_real_form& form, KGBElt x) : rf(form), val(x) @+{}
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
typedef std::auto_ptr<KGB_elt_value> KGB_elt_ptr;
typedef std::tr1::shared_ptr<KGB_elt_value> shared_KGB_elt;

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
{ int i = get<int_value>()->val;
  shared_real_form rf(get<real_form_value>());
  if (size_t(i)>=rf->val.KGB_size())
    throw std::runtime_error ("Inexistent KGB element: "+str(i));
@.Inexistent KGB element@>
  if (l!=expression_base::no_value)
    push_value(new KGB_elt_value(rf,i));
}

@ Working with KGB elements often requires having access to its real form, and
sometimes of its number within its real form.

@< Local function def...@>=
void decompose_KGB_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  if (l==expression_base::no_value)
    return;
  push_value(x->rf);
  push_value(new int_value(x->val));
  if (l==expression_base::single_value)
    wrap_tuple(2);
}

@ Three important attributes of KGB elements are the associated Cartan class,
the root datum involution and the length.

@< Local function def...@>=
void KGB_Cartan_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  if (l!=expression_base::no_value)
    push_value(new Cartan_class_value(x->rf->parent,kgb.Cartan_class(x->val)));
}

void KGB_involution_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  const ComplexReductiveGroup& G=x->rf->val.complexGroup();
  if (l!=expression_base::no_value)
    push_value(new matrix_value(G.involutionMatrix(kgb.involution(x->val))));
}

void KGB_length_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  if (l!=expression_base::no_value)
    push_value(new int_value(kgb.length(x->val)));
}

@ Cross actions and (inverse) Cayley transforms define the structure of a KGB
set, and we make them available as functions. The inverse Cayley transform may
be double valued of which we only report the first one; the user can easily
test whether it was double valued by applying cross action by the same
generator to the result and testing whether it gives a new KGB element (which
then is the other value of the Cayley transform).

@< Local function def...@>=
void KGB_cross_wrapper(expression_base::level l)
{ shared_KGB_elt x = get_own<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  int s = get<int_value>()->val;
  if (static_cast<unsigned>(s)>=kgb.rank())
    throw std::runtime_error ("Illegal simple reflection: "+str(s));
  if (l==expression_base::no_value)
    return;
  x->val= kgb.cross(s,x->val); // do cross action
  push_value(x);
}
@)
void KGB_Cayley_wrapper(expression_base::level l)
{ shared_KGB_elt x = get_own<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  int s = get<int_value>()->val;
  if (static_cast<unsigned>(s)>=kgb.rank())
    throw std::runtime_error ("Illegal simple reflection: "+str(s));
  if (l==expression_base::no_value)
    return;
  if (kgb.cayley(s,x->val)!=UndefKGB) // when defined
    x->val= kgb.cayley(s,x->val); // do Cayley transform
  push_value(x);
}
@)
void KGB_inv_Cayley_wrapper(expression_base::level l)
{ shared_KGB_elt x = get_own<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  int s = get<int_value>()->val;
  if (static_cast<unsigned>(s)>=kgb.rank())
    throw std::runtime_error ("Illegal simple reflection: "+str(s));
  if (l==expression_base::no_value)
    return;
  if (kgb.inverseCayley(s,x->val).first!=UndefKGB) // when defined, do first
    x->val= kgb.inverseCayley(s,x->val).first; // inverse Cayley transform
  push_value(x);
}

@ One also needs to be able find out the status of simple roots. Although
somewhat low-level, the simplest thing is to export the
|gradings::Status::Value| as an integer value. It seems however useful to
change complex ascents from $0$ to $4$, so that the coding is 0:~Complex
descent, 1:~imaginary compact, 2:~real, 3:~imaginary non-compact, 4:Complex
ascent. This way a value $v$ is a descent if |v<3|, imaginary if |v%2==1|,
Complex if |v%4==0|, Cayley transform defined if |v==3|.

@< Local function def...@>=
void KGB_status_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  int s = get<int_value>()->val;
  if (static_cast<unsigned>(s)>=kgb.rank())
    throw std::runtime_error ("Illegal simple reflection: "+str(s));
  if (l==expression_base::no_value)
    return;
  unsigned stat=kgb.status(s,x->val);
  push_value(new int_value
    (stat==0 and not kgb.isDescent(s,x->val) ? 4 : stat));
}

@ For a given KGB element, all imaginary roots are classified into compact and
non-compact roots. While this could be deduced from other attributes of the
element, this would be laborious and error-prone, so we supply this
information as built-in function. The root is transmitted in coordinates
rather than as index into the list of positive roots, as this avoids possible
confusion about the interpretation of the index, and ambiguity with the
previous instance of |status|; it is also more convenient in those cases where
the root results from a computation (as opposed to selection from the list of
roots).

@< Local function def...@>=
void root_status_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  const ComplexReductiveGroup& G=x->rf->parent.val;
  const RootDatum& rd = G.rootDatum();
  shared_vector alpha_vec = get<vector_value>();
  RootNbr alpha = rd.rootNbr(alpha_vec->val);
  if (alpha>=rd.numRoots())
    throw std::runtime_error ("Vector is not a root");
@.Vector is not a root@>
  if (l==expression_base::no_value)
    return;
  unsigned stat=kgb::status(kgb,x->val, rd,alpha);
  if (stat==0) // $\alpha$ is a complex root, check if it is an ascent
  {
    RootNbr theta_alpha =
      G.involution_table().root_involution(kgb.inv_nr(x->val),alpha);
    if (rd.isPosRoot(alpha)==rd.isPosRoot(theta_alpha))
      stat = 4; // set status to complex ascent
  }
  push_value(new int_value (stat));
}

@ In order to ``synthesise'' a KGB element, one may specify a real form and a
rational vector that defines a grading of the corresponding imaginary roots,
in the same manner as for |Cartan_class_real_form| in section
@#Cartan_class_real_form@> above.

@< Local function def...@>=
void build_KGB_element_wrapper(expression_base::level l)
{ shared_rational_vector grading_shift = get<rational_vector_value>();
  shared_matrix theta = get<matrix_value>();
  shared_real_form rf = get<real_form_value>();

  { Coweight num(grading_shift->val.numerator().begin(),
                 grading_shift->val.numerator().end());
    if (theta->val.right_mult(num)!=num)
      throw std::runtime_error ("Grading coweight not fixed by involution");
@.Grading coweight not fixed@>
  }
  RatCoweight tv = grading_shift->val - rf->kgb().base_grading_vector();
  if (tv.normalize().denominator()!=1)
    throw std::runtime_error("Grading coweight not appropriate for real form");

  const ComplexReductiveGroup& G = rf->parent.val;
  TorusPart t(tv.numerator());
  TwistedInvolution tw = twisted_from_involution(G,theta->val);
  TitsElt a (G.titsGroup(),t,tw);

  KGBElt x = rf->kgb().lookup(a);
  if (x== rf->kgb().size())
    throw std::runtime_error("KGB element not present");

  if (l==expression_base::no_value)
    return;
  push_value(new KGB_elt_value (rf,x));
}



@ One can conjugate a KGB element by the distinguished involution of the inner
class.

@< Local function def...@>=
void KGB_twist_wrapper(expression_base::level l)
{ shared_KGB_elt x = get_own<KGB_elt_value>();
  const KGB& kgb=x->rf->kgb();
  if (l==expression_base::no_value)
    return;
  x->val= kgb.Hermitian_dual(x->val); // do twist
  push_value(x);
}

@ Here is a function that returns the vector of bits that distinguish KGB
elements in the fibre over the same involution.

@< Local function def...@>=
void torus_bits_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  if (l!=expression_base::no_value)
  { const KGB& kgb=x->rf->kgb();
    TorusPart t = kgb.torus_part(x->val);
    vector_value* p;
    vector_ptr result(p=new vector_value(int_Vector(kgb.torus_rank(),0)));
    for (unsigned int i=0; i<kgb.torus_rank(); ++i)
      p->val[i]=t[i];
    push_value(result);
  }
}

@ It might be more useful to export the same value in a form that takes into
account the base grading of the KGB set. The following function does that,
returning the result in the form of a rational vector that should be
interpreted in $({\bf Q}/2{\bf Z})^n$.

@< Local function def...@>=
void torus_factor_wrapper(expression_base::level l)
{ shared_KGB_elt x = get<KGB_elt_value>();
  if (l!=expression_base::no_value)
  { const KGB& kgb=x->rf->kgb();
    RatCoweight t = kgb.torus_part_global(x->val);
    push_value(new rational_vector_value(t));
  }
}

@ We haven't defined much useful other things to do with KGB elements yet, but
by popular request we make available the straightforward equality test.

@< Local function def...@>=
void KGB_equals_wrapper(expression_base::level l)
{ shared_KGB_elt y = get<KGB_elt_value>();
  shared_KGB_elt x = get<KGB_elt_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(x->rf==y->rf and x->val==y->val));
}

@ Finally we install everything related to $K\backslash G/B$ elements.
@< Install wrapper functions @>=
install_function(KGB_elt_wrapper,@|"KGB","(RealForm,int->KGBElt)");
install_function(decompose_KGB_wrapper,@|"%","(KGBElt->RealForm,int)");
install_function(KGB_cross_wrapper,@|"cross","(int,KGBElt->KGBElt)");
install_function(KGB_Cayley_wrapper,@|"Cayley","(int,KGBElt->KGBElt)");
install_function(KGB_inv_Cayley_wrapper,@|"inv_Cayley","(int,KGBElt->KGBElt)");
install_function(KGB_status_wrapper,@|"status","(int,KGBElt->int)");
install_function(root_status_wrapper,@|"status","(vec,KGBElt->int)");
install_function(build_KGB_element_wrapper,@|"KGB_elt"
		,"(RealForm,mat,ratvec->KGBElt)");
install_function(KGB_twist_wrapper,@|"twist","(KGBElt->KGBElt)");
install_function(KGB_Cartan_wrapper,@|"Cartan_class","(KGBElt->CartanClass)");
install_function(KGB_involution_wrapper,@|"involution","(KGBElt->mat)");
install_function(KGB_length_wrapper,@|"length","(KGBElt->int)");
install_function(torus_bits_wrapper,@|"torus_bits","(KGBElt->vec)");
install_function(torus_factor_wrapper,@|"torus_factor","(KGBElt->ratvec)");
install_function(KGB_equals_wrapper,@|"=","(KGBElt,KGBElt->bool)");


@*1 Blocks associated to a real form and a dual real form.
%
Although blocks as specified by a real form and a dual real form were designed
more for the original \.{atlas} interface than for use by \.{realex}, we
provide a data type for such blocks and some simple functionality associated
to them.

@< Includes... @>=

#include "blocks.h"

@ Like other data types we have seen, we include shared pointers to
parent objects to ensure these remain in existence as long as our block does;
in fact we include two such shared pointers, one for each real form. The |val|
field contains an actual |Block| instance, which is constructed van the
|Block_Value| is.

@< Type definitions @>=
struct Block_value : public value_base
{ const shared_real_form rf; const shared_dual_real_form dual_rf;
  Block val; // cannot be |const|, as Bruhat order may be generated implicitly
@)
  Block_value(const shared_real_form& form,
         const shared_dual_real_form& dual_form)
  : rf(form), dual_rf(dual_form)
  ,
  val(Block::build(rf->val.complexGroup(),
                   rf->form_index,dual_rf->form_index)) {}
  ~Block_value() @+{}
@)
  virtual void print(std::ostream& out) const;
  Block_value* clone() const @+ {@; return new Block_value(*this); }
  static const char* name() @+{@; return "KGB element"; }
private:
  Block_value(const Block_value& v)
  : rf(v.rf), val(v.val) @+{} // copy constructor
};
@)
typedef std::auto_ptr<Block_value> Block_ptr;
typedef std::tr1::shared_ptr<Block_value> shared_Block;

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
{ shared_dual_real_form drf=get<dual_real_form_value>();
  shared_real_form rf=get<real_form_value>();
@)
  if (&rf->parent.val!=&drf->parent.val)
    throw std::runtime_error @|
    ("Inner class mismatch between real form and dual real form");
@.Inner class mismatch...@>
  BitMap b(rf->parent.val.dual_Cartan_set(drf->val.realForm()));
  if (!b.isMember(rf->val.mostSplit()))
    throw std::runtime_error @|
    ("Real form and dual real form are incompatible");
@.Real form and dual...@>
  if (l!=expression_base::no_value)
    push_value(new Block_value(rf,drf));
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
  push_value(b->rf);
  push_value(b->dual_rf);
  if (l==expression_base::single_value)
    wrap_tuple(2);
}

void block_size_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  if (l!=expression_base::no_value)
    push_value(new int_value(b->val.size()));
}

void block_element_wrapper(expression_base::level l)
{ shared_int i(get<int_value>());
  shared_Block b = get<Block_value>();
  BlockElt z = i->val; // extract value unsigned
  if (z>=b->val.size())
    throw std::runtime_error("Block element out of range");
  if (l==expression_base::no_value)
    return;
  push_value(new KGB_elt_value(b->rf,b->val.x(z)));
  inner_class_value dual_ic(b->rf->parent,tags::DualTag());
  shared_real_form drf (new real_form_value(dual_ic,b->dual_rf->form_index));
  push_value(new KGB_elt_value(drf,b->val.y(z)));
  if (l==expression_base::single_value)
    wrap_tuple(2);
}


@ Finally we install everything related to blocks.

@< Install wrapper functions @>=
install_function(Fokko_block_wrapper,"block","(RealForm,DualRealForm->Block)");
install_function(decompose_block_wrapper,@|"%"
                ,"(Block->RealForm,DualRealForm)");
install_function(block_size_wrapper,"#","(Block->int)");
install_function(block_element_wrapper,"element","(Block,int->KGBElt,KGBElt)");

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
{ shared_real_form rf;
  StandardRepr val;
@)
  module_parameter_value(const shared_real_form& form, const StandardRepr& v)
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
typedef std::auto_ptr<module_parameter_value> module_parameter_ptr;
typedef std::tr1::shared_ptr<module_parameter_value> shared_module_parameter;

@ When printing a module parameter, we shall indicate a triple
$(x,\lambda,\nu)$ that defines it. Since we shall need to print |StandardRepr|
values in other contexts as well, we shall define an auxiliary output function
of such values first and then use that. The auxiliary function needs the
|Rep_context|, so we pass that explicitly. Using the same name |print| for the
auxiliary seems natural, but forces us to qualify upon calling.

@f nu NULL

@< Function def...@>=
std::ostream& print
  (std::ostream& out,const StandardRepr& val, const Rep_context& rc)
{ RootNbr witness; // dummy needed in call
  return
  out << @< Expression for adjectives that apply to a module parameter @>@;@;
@/    << " parameter (x="
      << val.x() << ",lambda="
      << rc.lambda(val) << ",nu="
      << rc.nu(val) << ')';
}
@)
void module_parameter_value::print(std::ostream& out) const
{@; interpreter::print(out,val,rc()); }

@ We provide one of the adjectives ``non-standard'' (when $\lambda$ fails to
be imaginary-dominant; in this case little can be done with the parameter),
``zero'' (the standard module vanishes due to the singular infinitesimal
character, namely by the presence of a singular compact simple-imaginary
root), ``non-final'' (the standard module is non-zero, but can be expressed in
terms of standard modules at more compact Cartans using a singular real root
satisfying the parity condition) or ``final'' (the good ones; the condition
implies ``standard'' an ``non-zero'').

@< Expression for adjectives... @>=
( rc.is_standard(val,witness) ?
    rc.is_zero(val,witness) ? "zero" :
      rc.is_final(val,witness) ? "final" : "non-final"
  : "non-standard" )

@ To make a module parameter, one should provide a KGB element~$x$, an
integral weight $\lambda-\rho$, and a rational weight~$\nu$. Since only its
projection on the $-\theta_x$-stable subspace is used, one might specify the
infinitesimal character $\gamma$ in the place of $\nu$.

@< Local function def...@>=
void module_parameter_wrapper(expression_base::level l)
{ shared_rational_vector nu(get<rational_vector_value>());
  shared_vector lam_rho(get<vector_value>());
  shared_KGB_elt x(get<KGB_elt_value>());
  if (nu->val.size()!=lam_rho->val.size()
      or nu->val.size()!=x->rf->val.rank())
    throw std::runtime_error ("Rank mismatch: ("
        +str(x->rf->val.rank())+","
	+str(lam_rho->val.size())+","+str(nu->val.size())+")");
@.Rank mismatch@>
  if (l!=expression_base::no_value)
    push_value(new@| module_parameter_value(x->rf,
      x->rf->rc().sr(x->val,lam_rho->val,nu->val)));
}

@ The following function, which we shall bind to the monadic operator `|%|',
transforms a parameter value into a tuple of values that defines it. This
tuple is not unique (since $\lambda$ is determined only modulo
$(1-\theta_x)X^*$) and this function should make a unique choice. Whether that
is really the case depends on the implementation of |StandardRepr|
though; the current code, like the printing routine, just uses the methods to
extract the components.

@< Local function def...@>=
void unwrap_parameter_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
  { push_value(new KGB_elt_value(p->rf,p->val.x()));
    push_value(new vector_value(p->rc().lambda_rho(p->val)));
    push_value(new rational_vector_value(p->rc().nu(p->val)));
    if (l==expression_base::single_value)
      wrap_tuple(3);
  }
}

@*2 Functions operating on module parameters.
A crucial attribute of module parameters is their infinitesimal character.

@< Local function def...@>=
void infinitesimal_character_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(new rational_vector_value(p->val.gamma()));
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
    push_value(new bool_value(p->rc().is_standard(p->val,witness)));
}

void is_zero_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  RootNbr witness;
  if (l!=expression_base::no_value)
    push_value(new bool_value(p->rc().is_zero(p->val,witness)));
}

void is_final_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  RootNbr witness;
  if (l!=expression_base::no_value)
    push_value(new bool_value(p->rc().is_final(p->val,witness)));
}

@ Before constructing (non-integral) blocks, it is essential that the
infinitesimal character is made dominant, and that possible singular (simple)
complex descents are applied to the parameter, as the result of the block
construction will only be mathematically meaningful under these circumstances.
This operation is therefore applied automatically in several places, but it is
useful to give the user an easy way to apply it explicitly.

Testing for equivalence of parameters amounts to testing for equality after
the parameters are made dominant (at least that claim was not contested at the
time of writing this). We provide this test, which will be bound to the
equality operator.

@< Local function def...@>=
void parameter_dominant_wrapper(expression_base::level l)
{ shared_module_parameter p = get_own<module_parameter_value>();
  if (l!=expression_base::no_value)
  {@; p->rc().make_dominant(p->val);
    push_value(p);
  }
}

void parameter_equivalent_wrapper(expression_base::level l)
{ shared_module_parameter q = get<module_parameter_value>();
  shared_module_parameter p = get<module_parameter_value>();
  if (p->rf!=q->rf)
    throw std::runtime_error @|
      ("Real form mismatch when testing equivalence");
  if (l!=expression_base::no_value)
  { StandardRepr z0=p->val, z1=q->val; // copy
    p->rc().make_dominant(z0); q->rc().make_dominant(z1);
    push_value(new bool_value(z0==z1));
  }
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
  int s = get<int_value>()->val;
  unsigned int r =
    rootdata::integrality_rank(p->rf->val.rootDatum(),p->val.gamma());
  if (static_cast<unsigned>(s)>=r)
    throw std::runtime_error
      ("Illegal simple reflection: "+str(s)+ ", should be <"+str(r));
  if (l!=expression_base::no_value)
    push_value(new module_parameter_value(p->rf,p->rc().cross(s,p->val)));
}

void parameter_Cayley_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  int s = get<int_value>()->val;
  unsigned int r =
    rootdata::integrality_rank(p->rf->val.rootDatum(),p->val.gamma());
  if (static_cast<unsigned>(s)>=r)
    throw std::runtime_error
      ("Illegal simple reflection: "+str(s)+ ", should be <"+str(r));
  if (l!=expression_base::no_value)
    push_value(new module_parameter_value(p->rf,p->rc().Cayley(s,p->val)));
}
void parameter_inv_Cayley_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  int s = get<int_value>()->val;
  unsigned int r =
    rootdata::integrality_rank(p->rf->val.rootDatum(),p->val.gamma());
  if (static_cast<unsigned>(s)>=r)
    throw std::runtime_error
      ("Illegal simple reflection: "+str(s)+ ", should be <"+str(r));
  if (l!=expression_base::no_value)
    push_value(new module_parameter_value(p->rf,p->rc().inv_Cayley(s,p->val)));
}


@ The library can also compute orientation numbers for parameters.

@< Local function def...@>=
void orientation_number_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
    push_value(new int_value(p->rc().orientation_number(p->val)));
}

@ Here is a function that computes a list of positive rational values $t\leq1$
such that the parameter obtained by replacing the continuous part~$\nu$ of
by~$t\nu$ is not topmost in its block, so that deformation of the parameter to
$t\nu$ will produce a non-trivial decomposition.

@< Local function def...@>=
void reducibility_points_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  if (l!=expression_base::no_value)
  {
    RationalList rp = p->rc().reducibility_points(p->val);
      // method normalises rationals
    row_ptr result(new row_value(rp.size()));
    for (size_t i=0; i<rp.size(); ++i)
      result->val[i]=shared_value(new rat_value(rp[i]));
    push_value(result);
  }
}

@ One of the main reasons to introduce module parameter values is that they
allow computing a block, whose elements are again given by module parameters.
Before we define the functions the do that, let us define a common function
they will use to test a parameter for validity. Even though we shall always
have a pointer available when we call |test_standard|, we define this function
to take a reference (requiring us to write a dereferencing at each call),
because the type of pointer (shared or raw) available is not always the same.
The reference is of course not owned by |test_standard|. A similar test is
|is_nonzero_final|.

@< Local function def...@>=
void test_standard(const module_parameter_value& p)
{ RootNbr witness;
  if (p.rc().is_standard(p.val,witness))
    return;
  std::ostringstream os; p.print(os);
  os << "\nParameter not standard, negative on coroot #" << witness;
  throw std::runtime_error(os.str());
}

void test_nonzero_final(const module_parameter_value& p)
{ RootNbr witness; bool zero=p.rc().is_zero(p.val,witness);
  if (not zero and p.rc().is_final(p.val,witness))
    return; // nothing to report
  std::ostringstream os; p.print(os);
@/os << "\nParameter is " << (zero ? "zero" : "not final")
   @|  <<", as witnessed by coroot #" << witness;
  throw std::runtime_error(os.str());
}

@ Here is the first block generating function, which just reproduces to output
from the \.{atlas} program for the \.{nblock} command.

@< Local function def...@>=
void print_n_block_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p);
  BlockElt init_index; // will hold index in the block of the initial element
  non_integral_block block(p->rc(),p->val,init_index);
  *output_stream << "Parameter defines element " << init_index
               @|<< " of the following block:" << std::endl;
  block.print_to(*output_stream,true);
    // print block using involution expressions
  if (l==expression_base::single_value)
    wrap_tuple(0);
}

@ More interesting than printing the block is to return is to the user as a
list of parameter values. The following function does this, and adds as a
second result the index that the original parameter has in the resulting
block.

@< Local function def...@>=
void block_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p);
  if (l!=expression_base::no_value)
  {
    BlockElt start; // will hold index in the block of the initial element
    non_integral_block block(p->rc(),p->val,start);
    @< Push a list of parameter values for the elements of |block| @>
    push_value(new int_value(start));
    if (l==expression_base::single_value)
      wrap_tuple(2);
  }
}

@ Construction a list of values is a routine affair. This code must however
also construct a module parameter value for each element of |block|.

@< Push a list of parameter values for the elements of |block| @>=
{ row_ptr param_list (new row_value(block.size()));
  const RatWeight& gamma=block.gamma();
  for (BlockElt z=0; z<block.size(); ++z)
  { StandardRepr block_elt_param =
      p->rc().sr(block.parent_x(z),block.lambda_rho(z),gamma);
    param_list->val[z] =
	shared_value(new module_parameter_value(p->rf,block_elt_param));
  }
  push_value(param_list);

}

@ There is also a function that computes just a partial block.
@< Local function def...@>=
void partial_block_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p);
  if (l!=expression_base::no_value)
  {
    non_integral_block block(p->rc(),p->val);
    @< Push a list of parameter values for the elements of |block| @>
  }
}

@ Knowing the length in its block of a parameter is of independent interest.
@< Local function def...@>=
void param_length_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p);
  if (l!=expression_base::no_value)
    push_value(new int_value(p->rt().length(p->val)));
}

@ Here is a version of the |block| command that also exports the table of
Kazhdan-Lusztig polynomials for the block, in the same format as \\{raw\_KL}
that will be defined below.

@< Local function def...@>=
void KL_block_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p);
  if (l!=expression_base::no_value)
  {
    BlockElt start; // will hold index in the block of the initial element
    non_integral_block block(p->rc(),p->val,start);
    @< Push a list of parameter values for the elements of |block| @>
    push_value(new int_value(start));
    const kl::KLContext& klc = block.klc(block.size()-1,false);

    matrix_ptr M(new matrix_value(int_Matrix(klc.size())));
    for (size_t y=1; y<klc.size(); ++y)
      for (size_t x=0; x<y; ++x)
        M->val(x,y)= klc.KL_pol_index(x,y);
@)
    row_ptr polys(new row_value(0)); polys->val.reserve(klc.polStore().size());
    for (size_t i=0; i<klc.polStore().size(); ++i)
    {
      const kl::KLPol& pol = klc.polStore()[i];
      std::vector<int> coeffs(pol.size());
      for (size_t j=pol.size(); j-->0; )
        coeffs[j]=pol[j];
      polys->val.push_back(shared_value(new vector_value(coeffs)));
    }
@)
    vector_ptr length_stops(new vector_value(
       int_Vector(block.length(block.size()-1)+1)));
    length_stops->val[0]=0;
    for (size_t i=1; i<length_stops->val.size(); ++i)
      length_stops->val[i]=block.length_first(i);
@)
    unsigned n_survivors=0;
    for (BlockElt z=0; z<block.size(); ++z)
      if (block.survives(z))
        ++n_survivors;
    vector_ptr survivor(new vector_value(int_Vector(n_survivors)));
    { unsigned i=0;
      for (BlockElt z=0; z<block.size(); ++z)
        if (block.survives(z))
          survivor->val[i++]=z;
      assert(i==n_survivors);
    }
    matrix_ptr contributes_to(new matrix_value(
      int_Matrix(n_survivors,block.size(),0)));
    for (BlockElt z=0; z<block.size(); ++z)
    { BlockEltList sb = block.survivors_below(z);
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
    push_value(M);
    push_value(polys);
    push_value(length_stops);
    push_value(survivor);
    push_value(contributes_to);

    if (l==expression_base::single_value)
      wrap_tuple(7);
  }
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
  test_standard(*p);
  if (l!=expression_base::no_value)
  {
    non_integral_block block(p->rc(),p->val);
    @< Push a list of parameter values for the elements of |block| @>

    const kl::KLContext& klc = block.klc(block.size()-1,false);
    // compute KL polynomials, silently

    matrix_ptr M(new matrix_value(int_Matrix(klc.size())));
    for (size_t y=1; y<klc.size(); ++y)
      for (size_t x=0; x<y; ++x)
        M->val(x,y)= klc.KL_pol_index(x,y);
@)
    row_ptr polys(new row_value(0)); polys->val.reserve(klc.polStore().size());
    for (size_t i=0; i<klc.polStore().size(); ++i)
    {
      const kl::KLPol& pol = klc.polStore()[i];
      std::vector<int> coeffs(pol.size());
      for (size_t j=pol.size(); j-->0; )
        coeffs[j]=pol[j];
      polys->val.push_back(shared_value(new vector_value(coeffs)));
    }
@)
    vector_ptr length_stops(new vector_value(
       int_Vector(block.length(block.size()-1)+1)));
    length_stops->val[0]=0;
    for (size_t i=1; i<length_stops->val.size(); ++i)
      length_stops->val[i]=block.length_first(i);
@)
    unsigned n_survivors=0;
    for (BlockElt z=0; z<block.size(); ++z)
      if (block.survives(z))
        ++n_survivors;
    vector_ptr survivor(new vector_value(int_Vector(n_survivors)));
    { unsigned i=0;
      for (BlockElt z=0; z<block.size(); ++z)
        if (block.survives(z))
          survivor->val[i++]=z;
      assert(i==n_survivors);
    }
    matrix_ptr contributes_to(new matrix_value(
      int_Matrix(n_survivors,block.size(),0)));
    for (BlockElt z=0; z<block.size(); ++z)
    { BlockEltList sb = block.survivors_below(z);
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
    push_value(M);
    push_value(polys);
    push_value(length_stops);
    push_value(survivor);
    push_value(contributes_to);

    if (l==expression_base::single_value)
      wrap_tuple(6);
  }
}

@ Finally we install everything related to module parameters.
@< Install wrapper functions @>=
install_function(module_parameter_wrapper,@|"param"
                ,"(KGBElt,vec,ratvec->Param)");
install_function(unwrap_parameter_wrapper,@|"%"
                ,"(Param->KGBElt,vec,ratvec)");
install_function(infinitesimal_character_wrapper,@|"infinitesimal_character"
                ,"(Param->ratvec)");
install_function(real_form_of_parameter_wrapper,@|"real_form"
		,"(Param->RealForm)");
install_function(is_standard_wrapper,@|"is_standard" ,"(Param->bool)");
install_function(is_zero_wrapper,@|"is_zero" ,"(Param->bool)");
install_function(is_final_wrapper,@|"is_final" ,"(Param->bool)");
install_function(parameter_dominant_wrapper,@|"dominant" ,"(Param->Param)");
install_function(parameter_equivalent_wrapper,@|"=" ,"(Param,Param->bool)");
install_function(parameter_cross_wrapper,@|"cross" ,"(int,Param->Param)");
install_function(parameter_Cayley_wrapper,@|"Cayley" ,"(int,Param->Param)");
install_function(parameter_inv_Cayley_wrapper,@|"inv_Cayley"
                ,"(int,Param->Param)");
install_function(orientation_number_wrapper,@|"orientation_nr" ,"(Param->int)");
install_function(reducibility_points_wrapper,@|
		"reducibility_points" ,"(Param->[rat])");
install_function(print_n_block_wrapper,@|"print_block","(Param->)");
install_function(block_wrapper,@|"block" ,"(Param->[Param],int)");
install_function(partial_block_wrapper,@|"partial_block","(Param->[Param])");
install_function(param_length_wrapper,@|"length","(Param->int)");
install_function(KL_block_wrapper,@|"KL_block"
                ,"(Param->[Param],int,mat,[vec],vec,vec,mat)");
install_function(partial_KL_block_wrapper,@|"partial_KL_block"
                ,"(Param->[Param],mat,[vec],vec,vec,mat)");

@*1 Polynomials formed from parameters.
%
When working with parameters for standard modules, and notably with the
deformation formulas, the need arises to keep track of formal sums of
standard modules with coefficients of a type that allows keeping track of the
signatures of the modules. These coefficients, which we shall call split
integers and give the type \.{Split} in \.{realex}, are elements of the group
algebra over $\Z$ of a (cyclic) group of order~$2$.

@< Includes needed in the header file @>=
#include "arithmetic.h"

@*2 A class for split integers.
Although the necessary operations could easily be defined in the \.{realex}
programming language using pairs of integers, it is preferable to make them a
built-in type, since this allows distinguishing them from pairs of integers
used for other purposes, and to provide special output and conversion
facilities.

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
typedef std::auto_ptr<split_int_value> split_int_ptr;
typedef std::tr1::shared_ptr<split_int_value> shared_split_int;

@ Like for parameter values, we first define a printing function on the level
of a bare |Split_integer| value, which can be used in situations where the
method |split_int_value::print| cannot.

@< Function def...@>=
std::ostream& print (std::ostream& out, const Split_integer& val)
{@;
  return out << '(' << val.e() << (val.s()<0?'-':'+') << abs(val.s()) << "s)";
}
@)
void split_int_value::print(std::ostream& out) const @+
{@; interpreter::print(out,val); }

@ Here are the basic arithmetic operations.

@< Local function definitions @>=

void split_plus_wrapper(expression_base::level l)
{ Split_integer j=get<split_int_value>()->val;
  shared_split_int i=get_own<split_int_value>();
  if (l!=expression_base::no_value)
  {@; i->val+=j; push_value(i); }
}
@)
void split_minus_wrapper(expression_base::level l)
{ Split_integer j=get<split_int_value>()->val;
  shared_split_int i=get_own<split_int_value>();
  if (l!=expression_base::no_value)
  {@; i->val-=j; push_value(i); }
}
@)
void split_unary_minus_wrapper(expression_base::level l)
{ shared_split_int i=get_own<split_int_value>();
  if (l!=expression_base::no_value)
  {@; i->val.negate(); push_value(i); }
}
@)
void split_times_wrapper(expression_base::level l)
{ Split_integer j=get<split_int_value>()->val;
  Split_integer i=get<split_int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new split_int_value(i*j));
}

@ We also provide implicit conversions from integers or pairs of integers to
split integers, and an explicit operator for converting back to a pair.

@< Local function definitions @>=

void int_to_split_coercion()
{ int a=get<int_value>()->val;
  push_value(new split_int_value(Split_integer(a)));
}
@)
void pair_to_split_coercion()
{ push_tuple_components();
  int b=get<int_value>()->val;
  int a=get<int_value>()->val;
  push_value(new split_int_value(Split_integer(a,b)));
}
@)
void from_split_wrapper(expression_base::level l)
{ Split_integer si = get<split_int_value>()->val;
  if (l!=expression_base::no_value)
  {  push_value(new int_value(si.e()));
     push_value(new int_value(si.s()));
     if (l==expression_base::single_value)
       wrap_tuple(2);
  }
}

@*2 Class definition for virtual modules.
%
The library provides a type |repr::SR_poly| in which such sums can be
efficiently maintained. In order to use it we must have seen the header file
for the module \.{free\_abelian} on which the implementation is based.

@< Includes needed in the header file @>=
#include "free_abelian.h" // needed to make |repr::SR_poly| a complete type

@ Like for KGB elements, we maintain a shared pointer to the real form value, so
that it will be assured to survive as long as parameters for it exist.

@< Type definitions @>=
struct virtual_module_value : public value_base
{ shared_real_form rf;
  repr::SR_poly val;
@)
  virtual_module_value(const shared_real_form& form, const repr::SR_poly& v)
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
typedef std::auto_ptr<virtual_module_value> virtual_module_ptr;
typedef std::tr1::shared_ptr<virtual_module_value> shared_virtual_module;

@ When printing a virtual module value, we traverse the |std::map| that is
hidden in the |Free_Abelian| class template, and print individual terms using
the auxiliary function that was defined above for printing parameter values.
However when either all coefficients are integers or  coefficients are integer
multiples of~$s$, then we suppress the component that is always~$0$; this is
particularly useful if polynomials are used to encode $\Z$-linear combinations
of parameters.

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
    interpreter::print(out << '*',it->first,rc()); // print parameter
  }
}

@ To start off a |virtual_module_value|, one usually takes an empty sum, but
one needs to specify a real form to fill the |rf| field. The information
allows us to extract the real form from a virtual module even if it is empty.
We allow testing the number of terms of the sum, notably for testing the sum
to be empty.

@< Local function def...@>=
void virtual_module_wrapper(expression_base::level l)
{ shared_real_form rf (get<real_form_value>());
  if (l!=expression_base::no_value)
    push_value(new@|
       virtual_module_value(rf,repr::SR_poly(rf->rc().repr_less())));
}
@)
void real_form_of_virtual_module_wrapper(expression_base::level l)
{ shared_virtual_module m = get<virtual_module_value>();
  if (l!=expression_base::no_value)
    push_value(m->rf);
}
@)
void virtual_module_size_wrapper(expression_base::level l)
{ shared_virtual_module m = get<virtual_module_value>();
  if (l!=expression_base::no_value)
    push_value(new int_value(m->val.size()));
}

@ Here is function to extract the coefficient (multiplicity) of a given
parameter in a virtual module. It is bound to the array subscription syntax,
and therefor implemented as the |evaluate| method of the appropriate class
derived from |subscr_base|.

@ In a subscription of a polynomial by a parameter, the arguments are not
initially on the stack, but come from evaluating the |array| and |index|
fields of the |module_coefficient| expression.

@h "evaluator.h" // for |module_coefficient|

@< Function def... @>=
void module_coefficient::evaluate(level l) const
{ shared_virtual_module m = (array->eval(),get<virtual_module_value>());
  shared_module_parameter p = (index->eval(),get<module_parameter_value>());
  if (l!=expression_base::no_value)
    push_value(new split_int_value(m->val[p->val]));
}

@ We also allow implicitly converting a parameter to a virtual module.

@< Local function def...@>=
void param_to_poly()
{ shared_module_parameter p = get<module_parameter_value>();
@/test_standard(*p);
  const shared_real_form& rf=p->rf;
  push_value(new@|
    virtual_module_value(rf,rf->rc().expand_final(p->val)));
}

@ The main operations for virtual modules are addition and subtraction of
parameters, or of other virtual modules.

@< Local function def...@>=
void add_module_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  shared_virtual_module accumulator = get_own<virtual_module_value>();
@/test_standard(*p);
  if (accumulator->rf!=p->rf)
    throw std::runtime_error @|
      ("Real form mismatch when adding standard module to a module");
  if (l!=expression_base::no_value)
  @/{@; accumulator->val+= p->rc().expand_final(p->val);
    push_value(accumulator);
  }
}

void add_module_term_wrapper(expression_base::level l)
{ push_tuple_components(); // second argument is a pair |(coef,p)|
  shared_module_parameter p = get_own<module_parameter_value>();
  Split_integer coef=get<split_int_value>()->val;
  shared_virtual_module accumulator = get_own<virtual_module_value>();
@/test_standard(*p);
  if (accumulator->rf!=p->rf)
    throw std::runtime_error @|
      ("Real form mismatch when adding a term to a module");
  if (l!=expression_base::no_value)
  @/{@; accumulator->val.add_multiple(p->rc().expand_final(p->val),coef);
    push_value(accumulator);
  }
}

@ Although we initially envisioned allowing conversion from a list of terms to
a virtual module, this could not be defined since it is not possible to not
know the real form in case the list of terms is empty (a conversion in the
opposite direction is given below). Therefore we provide instead the addition
of an entire list of terms at once to a virtual module value.

@< Local function... @>=
void add_module_termlist_wrapper(expression_base::level l)
{ shared_row r = get<row_value>();
  shared_virtual_module accumulator = get_own<virtual_module_value>();
  if (l!=expression_base::no_value)
  { for (std::vector<shared_value>::const_iterator
           it=r->val.begin(); it!=r->val.end(); ++it)
    { const tuple_value* t = force<tuple_value>(it->get());
      Split_integer coef=force<split_int_value>(t->val[0].get())->val;
      const module_parameter_value* p =
        force<module_parameter_value>(t->val[1].get());
@/    test_standard(*p);
      if (accumulator->rf!=p->rf)
        throw std::runtime_error @|
          ("Real form mismatch when adding terms to a module");
      accumulator->val.add_multiple(p->rc().expand_final(p->val),coef);
     }
    push_value(accumulator);
  }
}

@ Naturally we also want to define addition and scalar multiplication of
virtual modules. Scalar multiplication potentially makes coefficients zero, in
which case the corresponding terms need to be removed to preserve the
invariant that no zero terms are stored in a virtual module. For integer
multiplication we just need to check for multiplication by $0$ and destroy the
whole module when this happens. However it is somewhat subtler for scalar
multiplication by split integers, because these have zero divisors. Therefore
we test each coefficient produced by multiplication in this case, and remove
the term when the coefficient becomes zero. We must take care to advance the
iterator ``manually'' before doing that, and as a consequence cannot as usual
advance the iterator in the |for| clause.

@< Local function... @>=
void add_virtual_modules_wrapper(expression_base::level l)
{ shared_virtual_module accumulator = get_own<virtual_module_value>();
  shared_virtual_module addend = get<virtual_module_value>();
  if (accumulator->rf!=addend->rf)
    throw std::runtime_error @|("Real form mismatch when adding two modules");
  if (l!=expression_base::no_value)
  { for (repr::SR_poly::const_iterator
         it=addend->val.begin(); it!=addend->val.end(); ++it)
      accumulator->val.add_term(it->first,it->second);
    push_value(accumulator);
  }
}
@)
void int_mult_virtual_module_wrapper(expression_base::level l)
{ shared_virtual_module m = get_own<virtual_module_value>();
  int c = get<int_value>()->val;
  if (l!=expression_base::no_value)
  { if (c==0)
      m->val.clear(); // avoid creating null terms into an |SR_poly|
    else
      for (repr::SR_poly::iterator it=m->val.begin(); it!=m->val.end(); ++it)
        it->second *= c;
    push_value(m);
  }
}
@)
void split_mult_virtual_module_wrapper(expression_base::level l)
{ shared_virtual_module m = get_own<virtual_module_value>();
  Split_integer c = get<split_int_value>()->val;
  if (l!=expression_base::no_value)
  { for (repr::SR_poly::iterator it=m->val.begin(); it!=m->val.end(); )
      // no |++it| here!
      if ((it->second *= c)==Split_integer(0,0))
	m->val.erase(it++); // advance, then delete the node just abandoned
      else ++it;
    push_value(m);
  }
}

@ Here is our principal application of virtual modules.
%
Using the computation of non-integral blocks, we can compute a deformation
formula for the given parameter. This also involves computing Kazhdan-Lusztig
polynomials, which happens inside the method |deformation_terms|, and produces
an expression for the ``deformed'' parameter (meaning $\nu$ is infinitesimally
decreased towards~$0$) in terms of certain other parameters found in the
block.

@< Local function def...@>=
void deform_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p);
  if (l!=expression_base::no_value)
  {
    non_integral_block block(p->rc(),p->val); // partial block construction
    repr::SR_poly terms
       = p->rt().deformation_terms(block,block.size()-1);

    virtual_module_ptr acc
      (new virtual_module_value(p->rf, repr::SR_poly(p->rc().repr_less())));
    for (repr::SR_poly::const_iterator it=terms.begin(); it!=terms.end(); ++it)
      acc->val.add_multiple(p->rc().expand_final(it->first),it->second);

    push_value(acc);
  }
}

@ Here is a recursive form of this deformation, which stores intermediate
results for efficiency in a |Rep_table| structure. Though this structure
should really be associated to and saved with the real form, we currently hold
it in a local variable of this wrapper function.

@< Local function def...@>=
void full_deform_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p);
  if (l!=expression_base::no_value)
  {
    repr::SR_poly result = p->rt().deformation(p->val);
    push_value (new virtual_module_value(p->rf,result));
  }
}

@ And here is another way to invoke the Kazhdan-Lusztig computations, which
given a parameter corresponding to $y$ will obtain the formal sum of
Kazhdan-Lusztig polynomials $P_{x,y}$ where $x$ ranges over the values in the
block of $y$ (or the Bruhat interval below $y$, where all those giving a
nonzero contribution are located), multiplied by a sign and evaluated at the
split integer unit~$s$ (since it appears that the information most frequently
needed can be extracted from that evaluation). In formula, this computes
$$
  \sum_{x\leq y}(-1)^{l(y)-l(x)}P_{x,y}[q:=s]
$$
@< Local function def...@>=
void KL_sum_at_s_wrapper(expression_base::level l)
{ shared_module_parameter p = get<module_parameter_value>();
  test_standard(*p);
  test_nonzero_final(*p);
  if (l!=expression_base::no_value)
  {
    repr::SR_poly result = p->rt().KL_column_at_s(p->val);
    push_value (new virtual_module_value(p->rf,result));
  }
}



@ Finally we install everything related to polynomials formed from parameters.
@< Install wrapper functions @>=
install_function(split_plus_wrapper,@|"+","(Split,Split->Split)");
install_function(split_minus_wrapper,@|"-","(Split,Split->Split)");
install_function(split_unary_minus_wrapper,@|"-","(Split->Split)");
install_function(split_times_wrapper,@|"*","(Split,Split->Split)");
install_function(from_split_wrapper,@|"%","(Split->int,int)");
install_function(virtual_module_wrapper,@|"null_module","(RealForm->ParamPol)");
install_function(real_form_of_virtual_module_wrapper,@|"real_form"
		,"(ParamPol->RealForm)");
install_function(virtual_module_size_wrapper,@|"#","(ParamPol->int)");
install_function(add_module_wrapper,@|"+","(ParamPol,Param->ParamPol)");
install_function(add_module_term_wrapper,@|"+"
		,"(ParamPol,(Split,Param)->ParamPol)");
install_function(add_module_termlist_wrapper,@|"+"
		,"(ParamPol,[(Split,Param)]->ParamPol)");
install_function(add_virtual_modules_wrapper,@|"+"
		,"(ParamPol,ParamPol->ParamPol)");
install_function(int_mult_virtual_module_wrapper,@|"*"
		,"(int,ParamPol->ParamPol)");
install_function(split_mult_virtual_module_wrapper,@|"*"
		,"(Split,ParamPol->ParamPol)");
install_function(deform_wrapper,@|"deform" ,"(Param->ParamPol)");
install_function(full_deform_wrapper,@|"full_deform","(Param->ParamPol)");
install_function(KL_sum_at_s_wrapper,@|"KL_sum_at_s","(Param->ParamPol)");


@*1 Kazhdan-Lusztig tables. We implement a simple function that gives raw
access to the table of Kazhdan-Lusztig polynomials.

@< Local function def...@>=
void raw_KL_wrapper (expression_base::level l)
{ shared_Block b = get<Block_value>();
  const Block& block = b->val;
  kl::KLContext klc(block); klc.fill();
@)
  if (l==expression_base::no_value)
    return;
  matrix_ptr M(new matrix_value(int_Matrix(klc.size())));
  for (size_t y=1; y<klc.size(); ++y)
    for (size_t x=0; x<y; ++x)
      M->val(x,y) = klc.KL_pol_index(x,y);
@)
  row_ptr polys(new row_value(0)); polys->val.reserve(klc.polStore().size());
  for (size_t i=0; i<klc.polStore().size(); ++i)
  {
    const kl::KLPol& pol = klc.polStore()[i];
    std::vector<int> coeffs(pol.size());
    for (size_t j=pol.size(); j-->0; )
      coeffs[j]=pol[j];
    polys->val.push_back(shared_value(new vector_value(coeffs)));
  }
@)
  std::vector<int> length_stops(block.length(block.size()-1)+1);
  length_stops[0]=0;
  for (size_t i=1; i<length_stops.size(); ++i)
    length_stops[i]=block.length_first(i);
@)
  push_value(M);
  push_value(polys);
  push_value(new vector_value(length_stops));
  if (l==expression_base::single_value)
    wrap_tuple(3);
}

@ For testing, it is useful to also have the dual Kazhdan-Lusztig tables.

@< Local function def...@>=
void raw_dual_KL_wrapper (expression_base::level l)
{ shared_Block b = get<Block_value>();
  const Block& block = b->val;
  Block dual_block = Block::build(b->dual_rf->val,b->rf->val);

  std::vector<BlockElt> dual=blocks::dual_map(block,dual_block);
  kl::KLContext klc(dual_block); klc.fill();
@)
  if (l==expression_base::no_value)
    return;
  matrix_ptr M(new matrix_value(int_Matrix(klc.size())));
  for (size_t y=1; y<klc.size(); ++y)
    for (size_t x=0; x<y; ++x)
      M->val(x,y) = klc.KL_pol_index(dual[y],dual[x]);
@)
  row_ptr polys(new row_value(0)); polys->val.reserve(klc.polStore().size());
  for (size_t i=0; i<klc.polStore().size(); ++i)
  {
    const kl::KLPol& pol = klc.polStore()[i];
    std::vector<int> coeffs(pol.size());
    for (size_t j=pol.size(); j-->0; )
      coeffs[j]=pol[j];
    polys->val.push_back(shared_value(new vector_value(coeffs)));
  }
@)
  std::vector<int> length_stops(block.length(block.size()-1)+1);
  length_stops[0]=0;
  for (size_t i=1; i<length_stops.size(); ++i)
    length_stops[i]=block.length_first(i);
@)
  push_value(M);
  push_value(polys);
  push_value(new vector_value(length_stops));
  if (l==expression_base::single_value)
    wrap_tuple(3);
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
void print_realweyl_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  shared_real_form rf(get<real_form_value>());
@)
  if (&rf->parent.val!=&cc->parent.val)
    throw std::runtime_error @|
    ("Inner class mismatch between arguments");
@.Inner class mismatch...@>
  BitMap b(rf->parent.val.Cartan_set(rf->val.realForm()));
  if (not b.isMember(cc->number))
    throw std::runtime_error @|
    ("Cartan class not defined for real form");
@.Cartan class not defined...@>
@)
  realredgp_io::printRealWeyl (*output_stream,rf->val,cc->number);
@)
  if (l==expression_base::single_value)
    wrap_tuple(0);
}

@)
void print_strongreal_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
@)
 realredgp_io::printStrongReal
    (*output_stream,cc->parent.val,cc->parent.interface,cc->number);
@)
  if (l==expression_base::single_value)
    wrap_tuple(0);
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
    wrap_tuple(0);
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
    wrap_tuple(0);
}

@)
void print_blocku_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  const Block& block = b->val;
@)
  block_io::printBlockU(*output_stream,block);
@)
  if (l==expression_base::single_value)
    wrap_tuple(0);
}

@ The \.{blockstabilizer} command has a slightly different calling scheme than
\.{block} and its friends, in that it requires a block and a Cartan class. The
lock itself is not actually used, just the real form and dial real form it
holds. The signature of |realredgp_io::printBlockStabilizer| is a bit strange,
as it requires a |RealReductiveGroup| argument for the real form, but only
numbers for the Cartan class and the dual real form (but this is
understandable, as information about the inner class must be transmitted in
some way). In fact it used to be even a bit stranger, in that the real form
was passed in the form of a |realredgp_io::Interface| value, a class (not to
be confused with |realform_io::Interface|, which does not specify a particular
real form) that we do not use in this program; since only the |realGroup|
field of the |realredgp_io::Interface| was used in
|realredgp_io::printBlockStabilizer|, we have changed its parameter
specification to allow it to be called easily here.

@< Local function def...@>=
void print_blockstabilizer_wrapper(expression_base::level l)
{ shared_Cartan_class cc(get<Cartan_class_value>());
  shared_Block b = get<Block_value>();
@)
  realredgp_io::printBlockStabilizer
   (*output_stream,b->rf->val,cc->number,b->dual_rf->val.realForm());
@)
  if (l==expression_base::single_value)
    wrap_tuple(0);
}

@ The function |print_KGB| takes only a real form as argument.

@h "kgb.h"
@h "kgb_io.h"

@< Local function def...@>=
void print_KGB_wrapper(expression_base::level l)
{ shared_real_form rf(get<real_form_value>());
@)
  *output_stream
    << "kgbsize: " << rf->val.KGB_size() << std::endl;
  const KGB& kgb=rf->kgb();
  kgb_io::var_print_KGB(*output_stream,rf->val.complexGroup(),kgb);
@)
  if (l==expression_base::single_value)
    wrap_tuple(0);
}

@ The function |print_X| even takes only an inner class as argument.

@< Local function def...@>=
void print_X_wrapper(expression_base::level l)
{ shared_inner_class ic(get<inner_class_value>());
@)
  ComplexReductiveGroup& G=ic->val;
  kgb::global_KGB kgb(G); // build global Tits group, "all" square classes
  kgb_io::print_X(*output_stream,kgb);
@)
  if (l==expression_base::single_value)
    wrap_tuple(0);
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
  kl::KLContext klc(block); klc.fill(false);
@)
  *output_stream
    << "Full list of non-zero Kazhdan-Lusztig-Vogan polynomials:\n\n";
  kl_io::printAllKL(*output_stream,klc,block);
@)
  if (l==expression_base::single_value)
    wrap_tuple(0);
}

@ The function |print_prim_KL| is a variation of |print_KL_basis|.

@< Local function def...@>=
void print_prim_KL_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  Block &block = b->val; // this one must be non-|const|
  kl::KLContext klc(block); klc.fill(false);
@)
  *output_stream
    << "Non-zero Kazhdan-Lusztig-Vogan polynomials for primitive pairs:\n\n";
  kl_io::printPrimitiveKL(*output_stream,klc,block);
@)
  if (l==expression_base::single_value)
    wrap_tuple(0);
}

@ The function |print_KL_list| is another variation of |print_KL_basis|, it
outputs just a list of all distinct Kazhdan-Lusztig-Vogan polynomials.

@< Local function def...@>=
void print_KL_list_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  const Block &block = b->val;
  kl::KLContext klc(block); klc.fill(false);
@)
  kl_io::printKLList(*output_stream,klc);
@)
  if (l==expression_base::single_value)
    wrap_tuple(0);
}

@ We close with two functions for printing the $W$-graph determined by the
polynomials computed. For |print_W_cells| we must construct one more object,
after having built the |klc::KLContext|.

@h "wgraph.h"
@h "wgraph_io.h"

@< Local function def...@>=
void print_W_cells_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  const Block &block = b->val;
  kl::KLContext klc(block); klc.fill(false);
@)

  wgraph::WGraph wg(klc.rank()); kl::wGraph(wg,klc);
  wgraph::DecomposedWGraph dg(wg);
@)
  wgraph_io::printWDecomposition(*output_stream,dg);
@)
  if (l==expression_base::single_value)
    wrap_tuple(0);
}

@ And here is |print_W_graph|, which just gives a variation on the output
routine of |print_W_cells|.

@< Local function def...@>=
void print_W_graph_wrapper(expression_base::level l)
{ shared_Block b = get<Block_value>();
  const Block &block = b->val;
  kl::KLContext klc(block); klc.fill(false);
@)

  wgraph::WGraph wg(klc.rank()); kl::wGraph(wg,klc);
@)
  wgraph_io::printWGraph(*output_stream,wg);
@)
  if (l==expression_base::single_value)
    wrap_tuple(0);
}


@ Here we install all remaining wrapper functions.

@< Install wrapper functions @>=
install_function(raw_KL_wrapper,@|"raw_KL","(Block->mat,[vec],vec)");
install_function(raw_dual_KL_wrapper,@|"dual_KL","(Block->mat,[vec],vec)");
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
