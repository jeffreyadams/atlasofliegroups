\def\emph#1{{\it#1\/}}
\def\Z{{\bf Z}}
\def\lcm{\mathop{\rm lcm}}
\let\eps=\varepsilon

@* Built-in types.
This file describes several built-in types related to the Atlas software,
which used by the interpreter as primitive types. It also defines built-in
functions that relate to these types.

@h "built-in-types.h"

@c
namespace atlas { namespace interpreter {
@< Global variable definitions @>@;
@< Declarations of local functions @>@;
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
@< Declarations of global variables @>@;
@< Declarations of exported functions @>@;
}@; }@;
#endif

@*1 Lie types.
@ Before we can define any types we must make sure the types like |value_base|
defined in \.{evaluator.w}, and like |lietype::LieType| defined
in~\.{lietype\_fwd.h}, are known in our header file.

@< Includes needed in the header file @>=
#include "evaluator.h"
#include "lietype_fwd.h"

@*2 The primitive type.
A first new type will be the Lie type for complex groups. This corresponds
to the type |lietype::LieType| in the Atlas software.

@h <stdexcept>
@< Type definitions @>=
struct Lie_type_value : public value_base
{ lietype::LieType value;
  Lie_type_value() : value(0) @+ {}
    // default constructor, produces empty type
  Lie_type_value(lietype::LieType t) : value(t) @+{}
  virtual void print(std::ostream& out) const;
  Lie_type_value* clone() const @+{@; return new Lie_type_value(*this); }
  void add_simple_factor (char,size_t)
    throw(std::bad_alloc,std::runtime_error);
private:
  Lie_type_value(const Lie_type_value& v) : value(v.value) @+{}
public:
  size_t rank() const;
  size_t semisimple_rank() const;
};

@ The type |lietype::LieType| is defined as |std::vector<SimpleLieType>| where
|SimpleLieType| stands for |std::pair<char,size_t>|. Therefore it
could take arbitrary values, not necessarily sensible ones. To remedy this we
make the method |add_simple_factor|, which is the only proposed way to build
up Lie types, check for the validity. We first define a small function to help
giving sensible error messages.

@h <sstream>
@< Function definitions @>=
std::string num(size_t n)
@+{@; std::ostringstream s; s<<n; return s.str(); }

@ Since the testing in the current (old) interface for the Atlas software is
clumsy to use, we perform our own tests here, emulating the ones in
|lietype::checkSimpleLieType|.

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
  value.push_back(lietype::SimpleLieType(c,rank));
}

@ Now we define a wrapper function that really builds a |Lie_type_value|.
We scan the string looking for sequences of a letter followed by a number.
We allow and ignore sequences of punctuation characters between the two.

@h <cctype>
@< Function definitions @>=
void scan_lie_type_wrapper() throw(std::bad_alloc,std::runtime_error)
{ string_value* s=get_string();
  std::auto_ptr<Lie_type_value> result(new Lie_type_value);
  size_t total_rank=0;
@/const char* p=s->value.c_str();
  while (std::ispunct(*p) || std::isspace(*p)) ++p;
  // skip any punctuation or space
  while (std::isalpha(*p))
  { char c=*p++;
    if (!std::isdigit(*p)) {@; --p; break; } // and |throw| a |runtime_error|
    size_t rank=*p++-'0'; @+
      while (std::isdigit(*p)) rank=10*rank+(*p++-'0');
    while (std::ispunct(*p) || std::isspace(*p)) ++p;
    // skip any punctuation or space
    result->add_simple_factor(c,rank);
      // this may |throw| a |runtime_error| as well
    if ((total_rank+=rank)>constants::RANK_MAX)
      throw std::runtime_error
      ("Total rank exceeds implementation limit "+num(constants::RANK_MAX));
@.Total rank exceeds...@>
  }
  if (*p!='\0')
    throw std::runtime_error
      ("Error in type string '"+s->value+"' for complex Lie type");
  push_value(result.release());
}

@ We shall call this function \.{Ctype} (for ``complex type'' as opposed to
``real type'').

@< Install wrapper functions @>=
install_function(scan_lie_type_wrapper,"Ctype","(string->Ctype)");

@*2 Printing Lie types.
Before we do anything more complicated with this primitive type, we must
ensure that we can print its values.

@< Function definitions @>=
void Lie_type_value::print(std::ostream& out) const
{ if (value.empty()) out << "empty Complex Lie type";
  else
  { out << "complex Lie type '" << value[0].first << value[0].second;
    for (size_t i=1; i<value.size(); ++i)
      out<< '.' << value[i].first << value[i].second;
    out << '\'';
  }
}

@*2 Auxiliary functions for Lie types.
As a service for testing routines we provide members that compute the rank
and the semisimple rank.

@< Function definitions @>=
size_t Lie_type_value::rank() const
{ size_t r=0;
  for (size_t i=0; i<value.size(); ++i) r+=value[i].second;
  return r;
}
@)
size_t Lie_type_value::semisimple_rank() const
{ size_t r=0;
  for (size_t i=0; i<value.size(); ++i)
    if (value[i].first!='T') r+=value[i].second;
  return r;
}

@ We shall need a new variation of |get_int| to unpack complex Lie types.

@< Declarations of exported functions @>=
Lie_type_value* get_Ctype() throw(std::logic_error);

@~There is nothing surprising here.
@< Function definitions @>=
Lie_type_value* get_Ctype() throw(std::logic_error)
{ value_ptr p=pop_arg();
  Lie_type_value* result=dynamic_cast<Lie_type_value*>(p);
  if (result==0)
    {@; delete p; throw std::logic_error("Argument is not a Lie_type"); }
@.Argument is not a Lie\_type@>
  return result;
}

@ Here is a function that computes the Cartan matrix for a given complex type.

@h "prerootdata.h"
@< Function definitions @>=
void Cartan_matrix_wrapper()
{ std::auto_ptr<Lie_type_value> t(get_Ctype());
  std::auto_ptr<latmat_value>
    result(new latmat_value(latticetypes::LatticeMatrix()));
  prerootdata::cartanMatrix(result->value,t->value);
  push_value(result.release());
}


@ And here is a function that tries to do the inverse. We do not tests, so one
can play around and see what |dynkin::lieType| does with fairly random
matrices. For instance any null columns in the matrix will be mistaken for
factors of type $A_1$.

@h "dynkin.h"
@< Function definitions @>=
void type_of_wrapper ()
{ std::auto_ptr<latmat_value> m(get_mat());
  lietype::LieType t;
  dynkin::lieType(t,m->value);
  push_value(new Lie_type_value(t));
}
@~
@< Install wrapper functions @>=
install_function(Cartan_matrix_wrapper,"Cartan_matrix","(Ctype->mat)");
install_function(type_of_wrapper,"type_of","(mat->Ctype)");

@*2 Finding lattices for a given Lie type.
The following function is copied from the atlas software sources, from the file
\.{io/interactive\_lattice.cpp}, since linking to the corresponding object
file seemed to pull in a major part of the Atlas software, while the function
we want to use requires nothing but access to the Smith normal form algorithm,
which being a template can be done by looking at the header file. Might we
suggest that this function that performs no interaction would have been better
placed elsewhere? We quote from its original documentation:

The Lie type is defined as a sequence of simple and torus factors inside |lt|.
This function constructs a ``blockwise Smith normal'' basis for the root
lattice inside the weight lattice (i.e., it does just that for each semisimple
block, and returns the canonical basis for the torus blocks.)

The purpose of doing this blockwise instead of globally is to permit a
better reading of the quotient group: this will be presented as a sequence
of factors, corresponding to each simple block.

% Of course, it is still assumed that the caller is a ``savvy'' user, and
% knows that some blocks will not contribute (because the adjoint group is
% simply connected), and that blocks $D_n$ with $n$ even, contribute
% \emph{two} factors $\Z/2\Z$. Torus blocks $T_n$ contribute $n$ factors $\Z$
% (this will be reflected in the ``missing'' invariant factors.)

@h "smithnormal_def.h"

@< Function definitions @>=
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
      latticetypes::operator+=(bp[r-2],bp[r-1]);
    }
    bp += r;
  }
}

@ And here is a function that applies the previous function, and returns a
block-wise Smith basis for the transposed Cartan matrix for a given complex
type, and the corresponding invariant factors. In case of torus factors this
description should be interpreted in the sense that the Smith basis for the
corresponding block is the standard basis and the invariant factors are null.

@< Function definitions @>=
void Smith_Cartan_wrapper()
{ std::auto_ptr<Lie_type_value> t(get_Ctype());
  std::auto_ptr<latmat_value>
    m(new latmat_value(latticetypes::LatticeMatrix()));
  std::auto_ptr<weight_value> inv_factors
     @| (new weight_value(latticetypes::CoeffList(0)));
  latticetypes::WeightList b;
  smithBasis(inv_factors->value,b,t->value);
  latticetypes::LatticeMatrix basis(b); // convert from list of vectors
  m->value.swap(basis);
  push_value(m.release()); push_value(inv_factors.release()); wrap_tuple(2);
}

@ The result of that function can serve among other things to help specifying
a basis for subgroup of the quotient of the original weight
lattice~$\tilde{X}$ by the sub-lattice spanned by the roots (the columns of
the transposed Cartan matrix). For that purpose only invariant factors unequal
to~$1$ and the corresponding columns of the Smith basis are of interest.
Therefore the following wrapper function, which may operate directly on the
result of the previous one, filters out the invariant factors~$1$ and the
corresponding columns.

This function is particular in that it returns exactly the kind of arguments
that it requires. Therefore rather that calling |push_tuple_components| which
would pop and destroy the tuple, we just extract its fields while leaving it
on the stack; its fields are modified in place, so there is no |push_value| to
return the value either. Therefore we do not own the arguments during the
call, and auto-pointers should not be used; incidentally |erase| and
|eraseColumn| should have no reason whatsoever to throw an exception. We do
take care to check that the argument is a $2$-tuple, but to avoid even more
|throw| expressions, we use |get_mat| and |get_vec| after temporarily pushing
the components onto the stack.

@< Function definitions @>=
void filter_units_wrapper ()
{ tuple_value* t=get_tuple(); push_value(t); // get a non-owned copy
  if (t->value.size()!=2) throw std::logic_error("Argument is not a pair");
  execution_stack.push_back(t->value[0]);
  latmat_value* basis=get_mat();
  execution_stack.push_back(t->value[1]);
  weight_value* inv_f=get_vec();
  if (inv_f->value.size()!=basis->value.numColumns())
    throw std::runtime_error("Size mismatch "+
      num(inv_f->value.size())+':'+num(basis->value.numColumns())+
      " in filter_units");
@)
  size_t i=0;
  while (i<inv_f->value.size())
    if (inv_f->value[i]!=1) ++i; // keep invariant factor and column
    else
    {@; inv_f->value.erase(inv_f->value.begin()+i);
        basis->value.eraseColumn(i);
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

The algorithm is quite simple. After finding a Smith basis~$(b_j)_{j\in[m]}$
and invariant factors~$(\lambda_j)_{j\in[r]}$ for the lattice spanned by the
columns of~$M$ (where $r$ is its rank), we compute the dual basis
$(b^*_j)_{j\in[m]}$ as the columns of transpose inverse matrix. Then we know
that $b^*_j$ applied to any column of~$M$ lies in $\lambda_j\Z$ (where we take
$\lambda_j=0$ for $j\notin[r]$), and we can find $A$ from the matrix giving
the dual basis by multiplying column~$j$ by
$\lcm(d,\lambda_j)/\lambda_j=d/\gcd(d,\lambda_j)$, for $j\in[r]$.

@f lambda NULL

@h "arithmetic.h"
@h "lattice.h"
@h "latticetypes.h"
@h "matrix.h"
@h "smithnormal.h"

@< Function definitions @>=
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

  latticetypes::LatticeMatrix A(b); A.invert(); A.transpose();

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

@< Function definitions @>=
void ann_mod_wrapper()
{ push_tuple_components();
@/std::auto_ptr<int_value> d(get_int());
@/std::auto_ptr<latmat_value> m(get_mat());
@)
  latticetypes::LatticeMatrix A=
    annihilator_modulo(m->value,d->value);
  m->value.swap(A);
  push_value(m.release());
}

@ Next a simple administrative routine, needed here because we cannot handle
matrices in our programming language yet. Once one has computed a new lattice,
in the form of vectors to replace those selected by \.{filter\_units} from the
result of \.{Smith\_Cartan}, possibly with the help of \.{ann\_mod}, one needs
to make the replacement. The following function does this, taking its first
two arguments as the result of \.{Smith\_Cartan}, and the third a matrix whose
columns are to be substituted. The second argument serves only to determine,
by the place of its non-unit entries, where the insertion has to take place.
In fact this is so simple that we define the wrapper function directly. And in
fact it will be more practical to take as argument a pair of a pair as
returned by \.{Smith\_Cartan} and a matrix.

@< Function definitions @>=
void replace_gen_wrapper ()
{ push_tuple_components(); // a pair
  std::auto_ptr<latmat_value> new_generators(get_mat());
  push_tuple_components(); // a pair as returned by \.{Smith\_Cartan}
  std::auto_ptr<weight_value> inv_f(get_vec());
  std::auto_ptr<latmat_value> old_generators(get_mat());
@)
  if (new_generators->value.numRows()!=old_generators->value.numRows())
    throw std::runtime_error("Column lengths do not match in replace_gen");
@.Column lengths do not match@>
  if (inv_f->value.size()!=old_generators->value.numColumns())
    throw std::runtime_error("Number of columns mismatch in replace_gen");
@.Size mismatch in replace\_gen@>
@)
  for (size_t j=0,k=0; j<inv_f->value.size(); ++j)
    if (inv_f->value[j]!=1)
       // replace column |j| by column |k| from |new_generators|
    { if (k>=new_generators->value.numColumns())
        throw std::runtime_error
          ("Not enough replacement columns in replace_gen");
@.Not enough replacement columns@>
      for (size_t i=0; i<old_generators->value.numRows(); ++i)
        old_generators->value(i,j)=new_generators->value(i,k);
      ++k;
    }
  push_value(old_generators.release());
}

@ All that remains is installing the wrapper functions.

@< Install wrapper functions @>=
install_function(Smith_Cartan_wrapper,"Smith_Cartan","(Ctype->mat,vec)");
install_function(filter_units_wrapper,"filter_units","(mat,vec->mat,vec)");
install_function(ann_mod_wrapper,"ann_mod","(mat,int->mat)");
install_function(replace_gen_wrapper,"replace_gen",
                 "((mat,vec),mat->mat)");

@*1 Root data.
Now we are ready to introduce a new primitive type for root data.

@< Includes needed in the header file @>=
#include "rootdata.h"

@ The root datum type is laid out just like other primitive types.

@< Type definitions @>=
struct root_datum_value : public value_base
{ rootdata::RootDatum value;
  root_datum_value(const rootdata::RootDatum v) : value(v) @+ {}
  virtual void print(std::ostream& out) const;
  root_datum_value* clone() const @+{@; return new root_datum_value(*this); }
private:
  root_datum_value(const root_datum_value& v) : value(v.value) @+{}
};

@*2 Printing root data.
We shall not print the complete information contained in the root datum.
However we do exercise the routines in \.{dynkin.cpp} to determine the type of
the root datum from its Cartan matrix, as is done in the function
|type_of_wrapper| above. Contrary to |prerootdata::cartanMatrix|, the function
|rootdata::cartanMatrix| produces nothing for torus factors, so we add one
cumulative torus factor at the end.

@< Function definitions @>=
void root_datum_value::print(std::ostream& out) const
{ lietype::LieType t; latticetypes::LatticeMatrix Cartan;
  rootdata::cartanMatrix(Cartan,value);
  dynkin::lieType(t,Cartan);
  Lie_type_value type(t);
  if (!value.isSemisimple())
    type.add_simple_factor('T',value.rank()-value.semisimpleRank());
  out << (value.isSimplyConnected() ? "simply connected " : "")
      << (value.isAdjoint() ? "adjoint " : "")
      << "root datum of " << type;
}

@*2 Building a root datum.
Strangely, there is no function in the atlas to convert a value of type
|latticetypes::LatticeMatrix| into a |latticetypes::WeightList| value
(the sequence of its columns). We provide a conversion function here.

@< Function definitions @>=
latticetypes::WeightList columns (const latticetypes::LatticeMatrix M)
{ latticetypes::WeightList result(M.numColumns());
  for (size_t i=0; i<M.numColumns(); ++i)
    M.column(result[i],i);
  return result;
}

@ To create a root datum value, the user must specify a Lie type and a square
matrix of the size of the rank of the root datum, which specifies generators
of the desired weight lattice as a sub-lattice of the lattice of weights
associated to the simply connected group of the type given. The given weights
should be independent and span at least the root lattice associated to the
type; for the moment we do not perform these tests however.

@< Function definitions @>=
void root_datum_wrapper()
{ push_tuple_components();
  std::auto_ptr<latmat_value> lattice(get_mat());
  std::auto_ptr<Lie_type_value> type(get_Ctype());
  if (lattice->value.numRows()!=lattice->value.numColumns() or
      lattice->value.numRows()!=type->rank())
    throw std::runtime_error
    ("Lattice matrix should be "+num(type->rank())+'x'+num(type->rank())+
     " for this type");
  @< Test whether the columns of |lattice->value| span the root lattice; if
     not |throw| a |std::runtime_error| @>
  prerootdata::PreRootDatum pre(type->value,columns(lattice->value));
  push_value(new root_datum_value(rootdata::RootDatum(pre)));
}

@ Since the construction of the (simple) roots in a |PreRootDatum| uses
multiplication by an inverse matrix producing an integral matrix without
testing whether the result was exact, we must laboriously do the same
computation before constructing the |PreRootDatum|, just to assure that the
construction will not fail. If the construction would succeed, even with an
incorrect result, we could test validity more easily be multiplying back the
given lattice basis with the roots in the |PreRootDatum|, and testing for
equality with the transpose Cartan matrix. But this would crash on division be
zero if the matrix were singular (this raises no exception, despite printing
message mentioning a \.{floating point exception}, so it cannot be caught). So
that approach would still require preliminary test for an invertible matrix,
and would end up being as much work as the code below. We note that in that
approach, which we had initially adopted, a subtle error is possible because
of roots that after basis transformation are rounded to zero would be mistaken
for columns coming from torus factors, and the |PreRootDatum| could therefore
have fewer roots than the semisimple rank of the type specified. But in case
the lattice matrix passes the test below, the |PreRootDatum| can safely be
constructed and will be correct in all respects.

@< Test whether the columns of |lattice->value| span the root lattice... @>=
{ latticetypes::LatticeMatrix M(lattice->value);
  latticetypes::LatticeCoeff d;
  M.invert(d);
  if (d==0) throw std::runtime_error
    ("Lattice matrix has dependent columns; in root_datum");
@/latticetypes::LatticeMatrix tC;
  prerootdata::cartanMatrix(tC,type->value); tC.transpose();
  matrix::leftProd(tC,M);
  if (!tC.divisible(d))
    throw std::runtime_error@|
      ("Given lattice does not contain the root lattice; in root_datum");
}

@ To emulate what is done in the atlas software, we write a function that
integrates some of the previous ones. Suppose $S=Smith\_Cartan(t)$ and
$(C,v)=filter\_units(S)$, then it builds a root datum from two inputs: a Lie
type~$t$ and a matrix~$M$ whose number of rows should equal the length of~$v$,
as follows. The entries in row~$i$ of~$M$ may be considered as numerators,
with corresponding denominator~$v_i$ (the order of the generator~$i$ of the
center), unless $v_i=0$ in which case the denominator is take to be the least
common multiple~$d$ of the nonzero values~$v_j$ (which is~$1$ if $v$ should be
completely null). If $M'$ is the result of bringing $M$ to a common
denominator~$d$ by multiplying each row~$i$ of~$M$ by $d/v_i$ whenever
$v_i\neq0$, then the root datum produced will be the result of
$root\_datum(t,replace\_gen(S,C*ann\_mod(M',d)))$.

This wrapper function is particular in that it does most of its work by
calling other wrapper functions; we essentially use the evaluator stack here
many times. In some cases we need to duplicate the top of the stack (this is
somewhat of a coincidence, since the bottom copy happens to be in place for a
future call). We do this by popping the stack into a named pointer variable,
pushing it straight back on and then pushing a cloned version; we could have
avoided popping and naming by calling |push_value(execution_stack.back())| at
the price of readability and some type safety in case we made some coding
errors here (since we are bypassing the type checker).

@< Function definitions @>=
void quotient_datum_wrapper()
{ push_tuple_components();
  std::auto_ptr<latmat_value> M(get_mat());
  Lie_type_value* type=get_Ctype(); // this copy is not owned
  push_value(type); // push back on stack for call of $root\_datum$
  push_value(type->clone()); // push a cloned copy for $Smith\_Cartan$
  Smith_Cartan_wrapper(); // compute $S=Smith\_Cartan(t)$
  tuple_value* S=get_tuple(); push_value(S); // for call of $replace\_gen$
  push_value(S->clone()); // push a cloned copy for call of $filter\_units$
  filter_units_wrapper(); // compute |(C,v)|
  push_tuple_components();
  std::auto_ptr<weight_value> v(get_vec());
  std::auto_ptr<latmat_value> C(get_mat());
  size_t d=1;
  for (size_t i=0; i<v->value.size(); ++i)
    if (v->value[i]!=0) d=arithmetic::lcm(d,v->value[i]);
  @< Test |M| against |v| and adapt its rows to the common denominator |d| @>
  push_value(C.release()); // for call of $mm\_prod$
  push_value(M.release()); // for call of $ann\_mod$
  push_value(new int_value(d)); // for call of $ann\_mod$
@/wrap_tuple(2); ann_mod_wrapper();
@/wrap_tuple(2); mm_prod_wrapper();
@/wrap_tuple(2); replace_gen_wrapper();
@/wrap_tuple(2); root_datum_wrapper();
}

@ As said above, |M| must have as many rows as |v| has entries.

@< Test |M| against |v| and adapt its rows to the common denominator |d| @>=
{ if (M->value.numRows()!=v->value.size())
    throw std::runtime_error ("Number "+num(M->value.numRows())+
      " of rows does not match number "+num(v->value.size())+
      " of kernel generators");
  for (size_t i=0; i<v->value.size(); ++i)
   if (v->value[i]!=0)
     for (size_t j=0; j<M->value.numColumns(); ++j)
       M->value(i,j)*=d/v->value[i];
}

@ We define two more wrappers with only a Lie type as argument, for building
the simply connected and the adjoint root data. They are similar to the
previous one in that they mostly call other wrapper functions. Only to make
the adjoint root datum we make sure to replace any null diagonal entries by
ones.

@< Function definitions @>=
void simply_connected_datum_wrapper()
{ Lie_type_value* type=get_Ctype(); push_value(type);
  push_value(new int_value(type->rank()));
  id_mat_wrapper();
@/wrap_tuple(2); root_datum_wrapper();
}
@)
void adjoint_datum_wrapper()
{ Lie_type_value* type=get_Ctype(); push_value(type);
  push_value(type->clone());
  Cartan_matrix_wrapper(); transmat_wrapper();
  latmat_value* M=get_mat();
  for (size_t i=0; i<type->rank(); ++i)
    if (M->value(i,i)==0) M->value(i,i)=1;
  push_value(M);
@/wrap_tuple(2); root_datum_wrapper();
}

@ Finally here is one more wrapper to make the root data for the special and
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

@< Function definitions @>=
void SL_wrapper()
{ std::auto_ptr<int_value> n(get_int());
  if (n->value<1) throw std::runtime_error("Non positive argument for GL");
  const size_t r=n->value-1;
  std::auto_ptr<Lie_type_value>type(new Lie_type_value());
  if (r>0) type->add_simple_factor('A',r);
  push_value(type.release());
  std::auto_ptr<latmat_value> lattice
     (new latmat_value(latticetypes::LatticeMatrix()));
  identityMatrix(lattice->value,r);
  for (size_t i=0; i<r-1; ++i) lattice->value(i,i+1)=-1;
  push_value(lattice.release());
@/wrap_tuple(2); root_datum_wrapper();
}
@)
void GL_wrapper()
{ std::auto_ptr<int_value> n(get_int());
  if (n->value<1) throw std::runtime_error("Non positive argument for GL");
  const size_t r=n->value-1;
  std::auto_ptr<Lie_type_value>type(new Lie_type_value());
  if (r>0) type->add_simple_factor('A',r);
  type->add_simple_factor('T',1);
  push_value(type.release());
  std::auto_ptr<latmat_value> lattice
     (new latmat_value(latticetypes::LatticeMatrix()));
  identityMatrix(lattice->value,r+1);
  for (size_t i=0; i<r; ++i)
  @/{@; lattice->value(n->value-1,i)=1; lattice->value(i,i+1)=-1; }
  push_value(lattice.release());
@/wrap_tuple(2); root_datum_wrapper();
}

@*2 Functions operating on root data.
We shall need a another variation of |get_int| to unpack root datum values.

@< Declarations of exported functions @>=
root_datum_value* get_root_datum() throw(std::logic_error);

@~
@< Function definitions @>=
root_datum_value* get_root_datum() throw(std::logic_error)
{ value_ptr p=pop_arg();
  root_datum_value* result=dynamic_cast<root_datum_value*>(p);
  if (result==0)
    {@; delete p; throw std::logic_error("Argument is not a root_datum"); }
@.Argument is not a Lie\_type@>
  return result;
}

@ The following functions allow us to look at the roots and co-roots stored in
a root datum value, and the associated Cartan matrix.

@< Function definitions @>=
void simple_roots_wrapper()
{ std::auto_ptr<root_datum_value> rd(get_root_datum());
  latticetypes::WeightList l
    (rd->value.beginSimpleRoot(),rd->value.endSimpleRoot());
  latmat_value* result=new latmat_value(atlas::latticetypes::LatticeMatrix(l));
  push_value(result);
}
@)
void simple_coroots_wrapper()
{ std::auto_ptr<root_datum_value> rd(get_root_datum());
  latticetypes::WeightList l
    (rd->value.beginSimpleCoroot(),rd->value.endSimpleCoroot());
  latmat_value* result=new latmat_value(atlas::latticetypes::LatticeMatrix(l));
  push_value(result);
}
@)
void Cartan_wrapper()
{ std::auto_ptr<root_datum_value> rd(get_root_datum());
  latticetypes::LatticeMatrix M;
  rootdata::cartanMatrix(M,rd->value);
  latmat_value* result=new latmat_value(M);
  push_value(result);
}

@ The following functions allow us to look at the roots and co-roots stored in
a root datum value.

@< Function definitions @>=
void roots_wrapper()
{ std::auto_ptr<root_datum_value> rd(get_root_datum());
  latticetypes::WeightList l
    (rd->value.beginRoot(),rd->value.endRoot());
  latmat_value* result=new latmat_value(atlas::latticetypes::LatticeMatrix(l));
  push_value(result);
}
@)
void coroots_wrapper()
{ std::auto_ptr<root_datum_value> rd(get_root_datum());
  latticetypes::WeightList l
    (rd->value.beginCoroot(),rd->value.endCoroot());
  latmat_value* result=new latmat_value(atlas::latticetypes::LatticeMatrix(l));
  push_value(result);
}

@ It is useful to have bases for the sum of the root lattice and the
coradical, and for the sum of the coroot lattice and the radical.

@< Function definitions @>=
void root_coradical_wrapper()
{ std::auto_ptr<root_datum_value> rd(get_root_datum());
  latticetypes::WeightList l
    (rd->value.beginSimpleRoot(),rd->value.endSimpleRoot());
  l.insert(l.end(),rd->value.beginCoradical(),rd->value.endCoradical());
  latmat_value* result=new latmat_value(atlas::latticetypes::LatticeMatrix(l));
  push_value(result);
}
@)
void coroot_radical_wrapper()
{ std::auto_ptr<root_datum_value> rd(get_root_datum());
  latticetypes::WeightList l
    (rd->value.beginSimpleCoroot(),rd->value.endSimpleCoroot());
  l.insert(l.end(),rd->value.beginRadical(),rd->value.endRadical());
  latmat_value* result=new latmat_value(atlas::latticetypes::LatticeMatrix(l));
  push_value(result);
}

@ And here is a simple function to dualise a root datum.

@h "tags.h"
@< Function definitions @>=
void dual_datum_wrapper()
{ std::auto_ptr<root_datum_value> rd(get_root_datum());
  rootdata::RootDatum dual(rd->value,tags::DualTag());
  rd->value.swap(dual);
  push_value(rd.release());
}


@ Let us install the above wrapper functions.

@< Install wrapper functions @>=
install_function(root_datum_wrapper,"root_datum","(Ctype,mat->RootDatum)");
install_function(quotient_datum_wrapper
		,"quotient_datum","(Ctype,mat->RootDatum)");
install_function(simply_connected_datum_wrapper
		,"simply_connected_datum","(Ctype->RootDatum)");
install_function(adjoint_datum_wrapper,"adjoint_datum","(Ctype->RootDatum)");
install_function(SL_wrapper,"SL","(int->RootDatum)");
install_function(GL_wrapper,"GL","(int->RootDatum)");
install_function(simple_roots_wrapper,"simple_roots","(RootDatum->mat)");
install_function(simple_coroots_wrapper,"simple_coroots","(RootDatum->mat)");
install_function(Cartan_wrapper,"Cartan","(RootDatum->mat)");
install_function(roots_wrapper,"roots","(RootDatum->mat)");
install_function(coroots_wrapper,"coroots","(RootDatum->mat)");
install_function(root_coradical_wrapper,"root_coradical","(RootDatum->mat)");
install_function(coroot_radical_wrapper,"coroot_radical","(RootDatum->mat)");
install_function(dual_datum_wrapper,"dual_datum","(RootDatum->RootDatum)");

@*1 A bitset type.
Here is a type designed to study what happens inside |dynkin|.

@< Includes... @>=
#include "bitset.h"
#include "setutils.h"

@
@< Type definitions @>=
struct bitset_value : public value_base
{ bitset::RankFlags value;
  bitset_value(unsigned long n) : value(bitset::RankFlags(n)) @+ {}
  virtual void print(std::ostream& out) const
   {@; out << "Bit set " << value.to_ulong(); }
  bitset_value* clone() const @+{@; return new bitset_value(*this); }
  void add_simple_factor (char,size_t)
    throw(std::bad_alloc,std::runtime_error);
private:
  bitset_value(const bitset_value& v) : value(v.value) @+{}
};

@ Here are wrapper functions that build, and set a bit set. It took long to
figure out that for a bitset~|b| and number~|n|, both |set(b,n)| and
|b.set(n)| are defined but do not means the same thing.

@< Function definitions @>=
void bitset_wrapper()
{ std::auto_ptr<int_value> n(get_int());
  push_value(new bitset_value(n->value));
}
@)
bitset_value* get_bitset() throw(std::logic_error)
{ value_ptr p=pop_arg();
  bitset_value* result=dynamic_cast<bitset_value*>(p);
  if (result==0)
    {@; delete p; throw std::logic_error("Argument is not a bitset"); }
@.Argument is not a bitset@>
  return result;
}
@)
void set_bits_wrapper()
{ push_tuple_components();
  std::auto_ptr<int_value> n(get_int());
  bitset_value* b=get_bitset();
  set(b->value,n->value);
  push_value(b);
}

@
@< Install wrapper functions @>=
install_function(bitset_wrapper,"bitset","(int->bitset)");
install_function(set_bits_wrapper,"set_bits","(bitset,int->bitset)");

@* Epilogue.
Here are some empty modules, which are place-holders for if anything in these
categories should be necessary (if not, the modules can be dropped).

@< Declarations of local functions @>=
@
@< Global variable definitions @>=
@
@< Declarations of global variables @>=

@ In execution, the main role of this compilation unit is to install the stuff
it defines into the tables. To that end we export an initialisation function

@< Declarations of exported functions @>=
void initialise_builtin_types();

@ This function is defined last, so that all wrapper functions are seen when
it is compiled; this avoids having to declare each wrapper function.

@< Function definitions @>=
void initialise_builtin_types()
@+{@; @< Install wrapper functions @>
}



@* Index.

% Local IspellDict: default
