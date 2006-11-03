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
its definitions into the tables. To that end we export its initialisation
function.

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
A first new type corresponds to the type |lietype::LieType| in the Atlas
software.

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

@ We first define a small function to help giving sensible error messages.

@h <sstream>
@< Local function definitions @>=
std::string num(size_t n)
@+{@; std::ostringstream s; s<<n; return s.str(); }

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
latter on input.

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
    while (rank-->0) value.push_back(lietype::SimpleLieType('T',1));
  else
    value.push_back(lietype::SimpleLieType(c,rank));
}

@ Now we define a wrapper function that really builds a |Lie_type_value|.
We scan the string looking for sequences of a letter followed by a number.
We allow and ignore sequences of punctuation characters between the two.

@h <cctype>
@< Local function definitions @>=
inline void skip_punctuation(const char* &p)
{@; while (std::ispunct(*p) || std::isspace(*p)) ++p;}
@)
void Lie_type_wrapper() throw(std::bad_alloc,std::runtime_error)
{ string_value* s=get_string();
  std::auto_ptr<Lie_type_value> result(new Lie_type_value);
  size_t total_rank=0;
@/const char* p=s->value.c_str(); skip_punctuation(p);
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
      ("Error in type string '"+s->value+"' for Lie type");
  push_value(result.release());
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
{ if (value.empty()) out << "empty Lie type";
  else
  {@; using basic_io::operator<<;
    out << "Lie type '" << value << '\'';
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

@ We shall need a new variation of |get_int| to unpack Lie types.

@< Declarations of exported functions @>=
Lie_type_value* get_Lie_type() throw(std::logic_error);

@~There is nothing surprising here.
@< Function definitions @>=
Lie_type_value* get_Lie_type() throw(std::logic_error)
{ value_ptr p=pop_arg();
  Lie_type_value* result=dynamic_cast<Lie_type_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not a Lie type"); }
@.Argument is not a Lie type@>
  return result;
}

@ Here is a function that computes the Cartan matrix for a given Lie type.

@h "prerootdata.h"
@< Local function definitions @>=
void Cartan_matrix_wrapper()
{ std::auto_ptr<Lie_type_value> t(get_Lie_type());
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
@< Local function definitions @>=
void type_of_Cartan_matrix_wrapper ()
{ std::auto_ptr<latmat_value> m(get_mat());
  lietype::LieType t;
  dynkin::lieType(t,m->value);
  push_value(new Lie_type_value(t));
}
@~
@< Install wrapper functions @>=
install_function(Cartan_matrix_wrapper,"Cartan_matrix","(LieType->mat)");
install_function(type_of_Cartan_matrix_wrapper
		,"type_of_Cartan_matrix","(mat->LieType)");

@*2 Finding lattices for a given Lie type.
The following function is copied from the Atlas software sources, from the file
\.{io/interactive\_lattice.cpp}, since linking to the corresponding object
file seemed to pull in a major part of the Atlas software, while the function
we want to use requires nothing but access to the Smith normal form algorithm,
which being a template can be done by looking at the header file. Might we
suggest that this function that performs no interaction would have been better
placed elsewhere? We quote from its original documentation:

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
      latticetypes::operator+=(bp[r-2],bp[r-1]);
    }
    bp += r;
  }
}

@ And here is a function that applies the previous function, and returns a
block-wise Smith basis for the transposed Cartan matrix for a given Lie
type, and the corresponding invariant factors. In case of torus factors this
description should be interpreted in the sense that the Smith basis for the
corresponding block is the standard basis and the invariant factors are null.

@< Local function definitions @>=
void Smith_Cartan_wrapper()
{ std::auto_ptr<Lie_type_value> t(get_Lie_type());
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

@ The result of |Smith_Cartan| can serve among other things to help specifying
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

@< Local function definitions @>=
void filter_units_wrapper ()
{ tuple_value* t=get_tuple(); push_value(t); // get a non-owned copy
  if (t->value.size()!=2) throw std::logic_error("Argument is not a pair");
  execution_stack.push_back(t->value[0]);@+
  latmat_value* basis=get_mat();
  execution_stack.push_back(t->value[1]);@+
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

@< Local function definitions @>=
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

@< Local function definitions @>=
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
meaningful for Lie types $A_n$ with $n\geq2$, any legal type $D_n$, or $E_n$
with $n\neq6$. Moreover in all cases except $D_n$ with $n$ even, this
designation is equivalent to |'s'|, and will be replaced by it. Therefore any
surviving type letter |'u'| corresponds to a type of the form $D_{2n}$.

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
@/std::auto_ptr<string_value> str(get_string());
@/std::auto_ptr<Lie_type_value> t(get_Lie_type());
@)std::auto_ptr<latmat_value> m
     (new latmat_value(latticetypes::LatticeMatrix()));
@/lietype::involution(m->value,t->value
                     ,transform_inner_class_type(str->value.c_str(),t->value));
  push_value(m.release());
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
@/std::auto_ptr<string_value> str(get_string());
@/std::auto_ptr<latmat_value> basis(get_mat());
@/std::auto_ptr<Lie_type_value> type(get_Lie_type());
@)
  size_t r=type->rank();
  if (basis->value.numRows()!=r or basis->value.numRows()!=r)
    throw std::runtime_error @|
    ("lattice matrix should be "+num(r)+'x'+num(r)+ " for this type");
  std::auto_ptr<latmat_value> m
     (new latmat_value(latticetypes::LatticeMatrix()));
@/lietype::involution
    (m->value,type->value
    ,transform_inner_class_type(str->value.c_str(),type->value));
@)m->value *= basis->value;
  latticetypes::LatticeCoeff d; basis->value.invert(d);
  matrix::leftProd(m->value,basis->value);
  if (d==0 or !m->value.divisible(d)) throw std::runtime_error
    ("inner class is not compatible with given lattice");
  m->value/=d; push_value(m.release());
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
|type_of_Cartan_matrix_wrapper| above. Contrary to
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
{ Lie_type_value type(type_of_datum(value));
  out << (value.isSimplyConnected() ? "simply connected " : "")
   @| << (value.isAdjoint() ? "adjoint " : "")
      << "root datum of " << type;
}

@ We also make the derivation of the type available by a wrapper function.
@< Local fun...@>=
void type_of_root_datum_wrapper()
{ std::auto_ptr<root_datum_value> rd(get_root_datum());
  push_value(new Lie_type_value(type_of_datum(rd->value)));
}

@*2 Building a root datum.
Strangely, there is no function in the Atlas to convert a value of type
|latticetypes::LatticeMatrix| into a |latticetypes::WeightList| value
(the sequence of its columns). We provide a conversion function here.

@< Local function definitions @>=
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
type.

@< Local function definitions @>=
void root_datum_wrapper()
{ push_tuple_components();
  std::auto_ptr<latmat_value> lattice(get_mat());
  std::auto_ptr<Lie_type_value> type(get_Lie_type());
  if (lattice->value.numRows()!=lattice->value.numColumns() or
      lattice->value.numRows()!=type->rank())
    throw std::runtime_error
    ("lattice matrix should be "+num(type->rank())+'x'+num(type->rank())+
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
be in place for the call to~$replace\_gen$). We do this by popping
the stack into a named pointer variable, pushing it straight back on and then
pushing a cloned version; we could have avoided popping and naming by calling
|push_value(execution_stack.back())| at the price of readability and some type
safety in case we made some coding errors here (since we are bypassing the
type checker).

@< Local function definitions @>=
void quotient_basis_wrapper()
{ push_tuple_components();
  std::auto_ptr<latmat_value> M(get_mat());
   // leave Lie type one stack for $Smith\_Cartan$
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
}

@ As said above, |M| must have as many rows as |v| has entries. But as a
service to the user we adapt its size to match that requirement if it has no
entries (either no rows or no columns).

@< Test |M| against |v| and adapt its rows to the common denominator |d| @>=
{ if (M->value.isEmpty()) M->value.resize(v->value.size(),0);
  if (M->value.numRows()!=v->value.size())
    throw std::runtime_error @| ("Number "+num(M->value.numRows())+
      " of rows does not match number "+num(v->value.size())+
      " of kernel generators");
  for (size_t i=0; i<v->value.size(); ++i)
   if (v->value[i]!=0)
     for (size_t j=0; j<M->value.numColumns(); ++j)
       M->value(i,j)*=d/v->value[i];
}

@ The function that integrates all is $quotient\_datum$; the call
$quotient\_datum(t,M)$ is equivalent to $root\_datum(t,quotient\_basis(t,M))$.

@< Local function definitions @>=
void quotient_datum_wrapper()
{ std::auto_ptr<tuple_value> args(get_tuple());
  push_value(args->value[0]->clone()); // the Lie type, for call of $root\_datum$
  push_value(args.release()); quotient_basis_wrapper();
@/wrap_tuple(2); root_datum_wrapper();
}

@ We define two more wrappers with only a Lie type as argument, for building
the simply connected and the adjoint root data. They are similar to the
previous one in that they mostly call other wrapper functions. Only to make
the adjoint root datum we make sure to replace any null diagonal entries by
ones.

@< Local function definitions @>=
void simply_connected_datum_wrapper()
{ Lie_type_value* type=get_Lie_type(); push_value(type);
  push_value(new int_value(type->rank()));
  id_mat_wrapper();
@/wrap_tuple(2); root_datum_wrapper();
}
@)
void adjoint_datum_wrapper()
{ Lie_type_value* type=get_Lie_type(); push_value(type);
  push_value(type->clone());
  Cartan_matrix_wrapper(); transpose_mat_wrapper();
  latmat_value* M=get_mat();
  for (size_t i=0; i<type->rank(); ++i)
    if (M->value(i,i)==0) M->value(i,i)=1;
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

@~Again, there are no surprises.
@< Function definitions @>=
root_datum_value* get_root_datum() throw(std::logic_error)
{ value_ptr p=pop_arg();
  root_datum_value* result=dynamic_cast<root_datum_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not a root_datum"); }
@.Argument is not a root datum@>
  return result;
}

@ The following functions allow us to look at the roots and co-roots stored in
a root datum value, and the associated Cartan matrix.

@< Local function definitions @>=
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
void datum_Cartan_wrapper()
{ std::auto_ptr<root_datum_value> rd(get_root_datum());
  latticetypes::LatticeMatrix M;
  rootdata::cartanMatrix(M,rd->value);
  latmat_value* result=new latmat_value(M);
  push_value(result);
}

@ The following functions allow us to look at the roots and co-roots stored in
a root datum value.

@< Local function definitions @>=
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

@< Local function definitions @>=
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
@< Local function definitions @>=
void dual_datum_wrapper()
{ std::auto_ptr<root_datum_value> rd(get_root_datum());
  rootdata::RootDatum dual(rd->value,tags::DualTag());
  rd->value.swap(dual);
  push_value(rd.release());
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

@*2 Analysing involutions.
Our constructor for the current built-in type must do checking to see that a
valid involution is entered, and an analysis of the involution to replace the
values otherwise obtained from the user interaction, in other words we want to
find which sequence of inner class letters could have been used to construct
this involution. For simple factors this is relatively easy, since one just
studies the permutation of the simple roots that the involution defines (if it
does not, it is not an automorphism of the based root datum). For torus
factors, this requires an analysis of the eigen-lattices of the involution,
which is performed in the class |tori::RealTorus|. The problem of finding the
inner class letters for the torus factors might not have a unique solution:
certainly involutions that correspond to a Complex inner class but on simple
factors that are not adjacent cannot be specified by inner class letters, and
on the other hand certain inner class letters are synonymous; even apart from
this there might, in cases where the group is not a direct product of its
derived group and its central torus, be involutions that can arise in several
different ways or not at all. Therefore the value we compute should be taken
as an informed guess.

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

@ The test that |M| is an involution ($M^2=I$) is certainly necessary when
|classify_involution(M)| is called independently; if called in the context of
checking the involution for an inner class it may seem redundant given the
fact that we shall also check that |M| induces an involutive permutation of
the roots, but it is not so in case there is a torus part.

@< Check that |M| is an $r\times{r}$ matrix defining an involution @>=
{ if (M.numRows()!=r or M.numColumns()!=r) throw std::runtime_error
    ("involution should be a "+num(r)+"x"+num(r)+" matrix, got a "
     +num(M.numRows())+"x"+num(M.numColumns())+" matrix");
  latticetypes::LatticeMatrix I,Q(M);
  identityMatrix(I,r); Q*=M; // $Q=M^2$
  if (!(Q==I)) throw std::runtime_error
      ("given transformation is not an involution");
}

@ That function might be of use to the user, so let us make a wrapper for it.
In fact we shall return the compact, Complex, and split ranks, in that order.

@< Local function def...@>=
void classify_wrapper()
{ std::auto_ptr<latmat_value> M(get_mat());
  size_t r=M->value.numRows();
  std::pair<size_t,size_t> p=classify_involution(M->value,r);
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
corresponds to each of its factors. The Lie type will in general be identical
to that of the root datum (and in particular any torus factors will come at
the end), but in case of Complex inner classes we may be forced to permute the
simple factors to make the identical factors associated to such classes
adjacent.

We do not apply the function |classify_involution| to the entire matrix,
although we shall use it below for a matrix defined for the central torus
part; we can however reuse the module that tests for being an involution here.

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
its rows and columns. We keep track of the number of fixed points of the
permutation, which in the semisimple case each account for a compact factor of
the real torus, at least in the simply connected or adjoint case.

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
@/bitset::RankFlagsList comps; dynkin::components(comps,diagr);
  std::vector<char> letter(comps.size());
  size_t nr_Complex_letters=0;
  for (size_t i=0; i<comps.size(); ++i)
  { bool equal_rank=true;
    for (bitset::RankFlags::iterator j=comps[i].begin(); j(); ++j)
      if (p[*j]!=*j) {@; equal_rank=false; break; }
    if (equal_rank) letter[i]='c';
      // identity on this component: compact
    else if(!comps[i][p[comps[i].firstBit()]])
      // exchange with other component: Complex
      {@; ++nr_Complex_letters; letter[i]='C'; }
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
  for (size_t l=0,i=0; l<comps.size(); ++l)
    if (letter[l]=='C')
    {@; pos[i]=l; first[i]=comps[l].firstBit(); ++i; }
  std::vector<size_t> buddy(nr_Complex_letters,nr_Complex_letters);
  for (size_t i=0; i<nr_Complex_letters; ++i)
    if (buddy[i]==nr_Complex_letters) // value was not set by buddy
    { size_t b=diagr.component(p[comps[pos[i]].firstBit()]).firstBit();
      for (size_t j=i+1; j<nr_Complex_letters; ++j)
        if (first[j]==b) {@; buddy[i]=j; buddy[j]=i; break; }
    }
  type.resize(comps.size()+(r-s));
  for (size_t l=0,k=0,i=0; l<comps.size(); ++l)
    if (letter[l]!='C')
    {@; type[k++]=type0[l]; inner_class.push_back(letter[l]); }
    else if (buddy[i]<i) ++i; // skip second member of Complex pair
    else // first member of Complex pair
    { type[k++]=type0[l]; type[k++]=type0[pos[buddy[i]]]; // should be equal
      inner_class.push_back('C'); ++i;
    }
  if (r>s)
  @< Add type letters and inner class symbols for the central torus @>
}

@ Finally we have to give inner class symbols for the central torus. While it
is not entirely clear that forming some quotient could not transform an inner
class specified as |"sc"| by the user into one that we will classify as |"C"|
or vice versa, the symbols we generate have a well defined interpretation: the
correspond to the real torus in the central torus defined by the involution.
The weight lattice of the central torus is the quotient of the full weight
lattice by the rational span of the root lattice, so we determine the action
of the involution on this quotient, and classify it using
|classify_involution|. To find the matrix of the action of the involution on
the quotient lattice, we find a Smith basis for the root lattice (of which the
first $r$ vectors span the lattice to be divided out), express the involution
of that basis (which will have zeros in the bottom-left $(r-s)\times{s}$
block), and extract the bottom-right $(r-s)\times(r-s)$ block.

@< Add type letters and inner class symbols for the central torus @>=
{ using latticetypes::WeightList; using latticetypes::LatticeMatrix;
  for (size_t k=comps.size(); k<comps.size()+(r-s); ++k)
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
This is the first built-in type where we deviate from the previously used
scheme of just having one data field |value| holding a value of a class
defined in the Atlas. The reason is that the copy constructor for
|complexredgp::ComplexReductiveGroup| is private (and nowhere defined), so
that the straightforward definition of a copy constructor for such a built-in
type would not work, and the copy constructor is necessary for the |clone|
method. So instead, we shall share the Atlas object when duplicating our
value, and maintain a reference count to allow destruction when the last copy
disappears. The reference count needs to be shared of course, and since the
links between the built in value and both the Atlas value it represents and
the reference count are indissoluble, we use references for the data fields.

The reference to the Atlas value is not constant, since at some point we need
to apply the method |fillCartan| to it; this is not defined as a |const|
method (in other words it is considered a manipulator), even though it could
have been, given the fact that it only modifies values at the other side of
the |d_cartan| pointer. The Atlas value so referred to then may change in the
course of the computation, to store information about the group as it is
computed, but apart from efficiency considerations the sharing should be
transparent (not noticeable by the user). Comments suggest however that the
numbering of Cartan subgroups may be affected by the order of operations
requested, so this transparency may not be complete.

@< Includes... @>=
#include "realform_io.h"

@~The main constructor takes a pointer to a |ComplexReductiveGroup| as
argument, which should come from a call to~|new|; this pointer will be owned
by the |inner_class_value| constructed and all values cloned from it, with
the last one to be destroyed calling |delete| for the pointer. This argument
is followed by the values that are computed by |check_involution| above.

Occasionally we shall need to refer to the dual inner class (for the dual
group); since the construction of an instance even without its complete set of
Cartan classes takes some work, we do not wish to repeat that every time it is
needed, so we create the dual |complexredgp::ComplexReductiveGroup| value upon
construction of the |inner_class_value| and store it in the |dual| field where
it will be available as needed.

Contrary to what was the case for other value types, the copy constructor is
public here, which will make it possible to store an |inner_class_value|
inside another class that depends on the object referenced by our |value|
field staying alive; the reference counting mechanism will then ensure that
this is the case as long as the object of that class exists.

@< Type definitions @>=
struct inner_class_value : public value_base
{ complexredgp::ComplexReductiveGroup& value;
  complexredgp::ComplexReductiveGroup& dual;
  size_t& ref_count;
  const realform_io::Interface interface,dual_interface;
@)
  inner_class_value(complexredgp::ComplexReductiveGroup*
                    ,lietype::LieType, lietype::InnerClassType);
  ~inner_class_value();
@)
  virtual void print(std::ostream& out) const;
  inner_class_value* clone() const @+
    {@; return new inner_class_value(*this); }
  inner_class_value(const inner_class_value& v);
  inner_class_value(const inner_class_value& v,tags::DualTag);
};

@ Here are the copy constructor and the destructor.
@< Function def...@>=
inner_class_value::inner_class_value(const inner_class_value& v)
: value(v.value), dual(v.dual), ref_count(v.ref_count)
, interface(v.interface), dual_interface(v.dual_interface)
{@; ++ref_count; }

inner_class_value::~inner_class_value()
{@; if (--ref_count==0) {@; delete &value; delete &dual; delete &ref_count;} }

@ The constructor installs the reference to the Atlas value, creates its dual
value to which a reference is stored, and allocates and initialises the
reference count. To initialise the |interface| and |dual_interface| fields, we
call the appropriate constructors with a |layout::Layout| structure provided
by a constructor that we added to it specifically for the purpose. This
constructor leaves the field holding a lattice basis empty, but this field is
unused by the constructors for |realform_io::Interface|.

@< Function def...@>=
inner_class_value::inner_class_value
  (complexredgp::ComplexReductiveGroup* g
  ,lietype::LieType lt, lietype::InnerClassType ict)
@/: value(*g)
, dual(*new complexredgp::ComplexReductiveGroup(*g,tags::DualTag()))
@/, ref_count(*new size_t(1))
@/, interface(*g,layout::Layout(lt,ict))
, dual_interface(*g,layout::Layout(lt,ict),tags::DualTag())
 @+ {}

@ We allow construction of a dual |inner_class_value|. Since it can share the
two fields referring to objects of type |complexredgp::ComplexReductiveGroup|
in the opposite order, we can consider it as a member of the same reference
counted family, and share the |ref_count| field. This means this constructor
is more like the copy constructor than like the main constructor, and in
particular the reference count is increased. The dual of the dual inner class
value will behave exactly like a copy of the original inner class.

@< Function def...@>=
inner_class_value::inner_class_value(const inner_class_value& v,tags::DualTag)
@/: value(v.dual), dual(v.value)
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
      << value.numRealForms() @| << " real "
      << (value.numRealForms()==1 ? "form" : "forms") @| << " and "
      << value.numDualRealForms() @| << " dual real "
      << (value.numDualRealForms()==1 ? "form" : "forms");
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
  std::auto_ptr<latmat_value> M(get_mat());
  std::auto_ptr<root_datum_value> rd(get_root_datum());
  std::pair<lietype::LieType,lietype::InnerClassType> cl
    =check_involution(M->value,rd->value);
  if (verbosity>0) @< Report the type and inner class found @>
  rootdata::RootDatum* rdp=new rootdata::RootDatum(rd->value);
  std::auto_ptr<complexredgp::ComplexReductiveGroup>@|
    G(new complexredgp::ComplexReductiveGroup(rdp,M->value));
  inner_class_value* result=new
    inner_class_value(G.get(),cl.first,cl.second);
  G.release(); push_value(result);
}

@ For understanding what is going on, the user may find it useful to look at
the Lie type and inner class that were determined from the root datum and the
involution given.

@< Report the type and inner class found @>=
{ Lie_type_value t(cl.first);
  std::cout << "Found " << t << ", and inner class '";
  for (size_t i=0; i<cl.second.size(); ++i)
    std::cout << cl.second[i];
  std::cout << "'.\n";
}

@ To simulate the functioning of the Atlas software, the function $set\_type$
takes as argument the name of a Lie type, a matrix giving kernel generators,
and a string describing the inner class. The evaluation of the call
$set\_type(lt,gen,ic)$ effectively consists of setting $t={\it Lie\_type(lt)}$,
${\it basis}={\it quotient\_basis(t,gen)}$, and then returning the value of
${\it fix\_involution (root\_datum (t,basis),based\_involution(t,basis,ic))}$.

@< Local function def...@>=
void set_type_wrapper()
{ push_tuple_components();
  std::auto_ptr<string_value> ict(get_string());
  std::auto_ptr<latmat_value> gen(get_mat());
  Lie_type_wrapper(); // convert string to Lie type
  std::auto_ptr<Lie_type_value> t(get_Lie_type());
@)
  push_value(t->clone()); push_value(gen.release());
  wrap_tuple(2); quotient_basis_wrapper();
  std::auto_ptr<latmat_value> basis(get_mat());
@)
  push_value(t->clone()); push_value(basis->clone());
  wrap_tuple(2); root_datum_wrapper();
@/push_value(t.release()); push_value(basis.release());
  push_value(ict.release());
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
  std::auto_ptr<string_value> ict(get_string());
  std::auto_ptr<root_datum_value> rd(get_root_datum());
@)
  push_value(rd->clone()); type_of_root_datum_wrapper();
  std::auto_ptr<Lie_type_value> t(get_Lie_type());
@)
  push_value(rd->clone());
  coroot_radical_wrapper(); transpose_mat_wrapper();
  std::auto_ptr<latmat_value> basis(get_mat());
@)
  push_value(rd.release());
@/push_value(t.release()); push_value(basis.release());
  push_value(ict.release());
  wrap_tuple(3); based_involution_wrapper();
@/wrap_tuple(2); fix_involution_wrapper();
}

@*2 Functions operating on complex reductive groups.
We shall need a another variation of |get_int| to unpack our new values.

@< Declarations of exported functions @>=
inner_class_value* get_complexredgp() throw(std::logic_error);

@~Again, there are no surprises.
@< Function definitions @>=
inner_class_value* get_complexredgp() throw(std::logic_error)
{ value_ptr p=pop_arg();
  inner_class_value* result=dynamic_cast<inner_class_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not a complex group"); }
@.Argument is not a complex group@>
  return result;
}

@ Here are our first functions that access a |inner_class_value|; they
recover the ingredients that were used in the construction, and construct the
dual inner class.

@< Local function def...@>=
void distinguished_involution_wrapper()
{ std::auto_ptr<inner_class_value> G(get_complexredgp());
  push_value(new latmat_value(G->value.distinguished()));
}

void root_datum_of_inner_class_wrapper()
{ std::auto_ptr<inner_class_value> G(get_complexredgp());
  push_value(new root_datum_value(G->value.rootDatum()));
}

void dual_inner_class_wrapper()
{ std::auto_ptr<inner_class_value> G(get_complexredgp());
  push_value(new inner_class_value(*G,tags::DualTag()));
}

@ More interestingly, let us extract the list of names of the real forms.
This uses the interface fields stored in the value. Since they exist for both
the group itself and for the dual group, we define an auxiliary function that
produces the list, and then use it twice.

@< Local function def...@>=
void push_name_list(const realform_io::Interface& interface)
{ std::auto_ptr<row_value> result
    (new row_value(std::vector<value_ptr>()));
  for (size_t i=0; i<interface.numRealForms(); ++i)
    result->value.push_back(new string_value(interface.typeName(i)));
  push_value(result.release());
}

void form_names_wrapper()
{@; std::auto_ptr<inner_class_value> G(get_complexredgp());
  push_name_list(G->interface);
}

void dual_form_names_wrapper()
{@; std::auto_ptr<inner_class_value> G(get_complexredgp());
  push_name_list(G->dual_interface);
}

@ And now, our first function that really simulates something that can be done
using the atlas interface, and that uses more than a root datum. This is the
\.{blocksizes} command from \.{mainmode.cpp}, which uses
|innerclass_io::printBlockSizes|, but we have to rewrite it to avoid calling
an output routine. A subtle difference is that we use a matrix of integers
rather than of unsigned long integers to collect the block sizes; this avoids
having to define a new primitive type.

@< Local function def...@>=
void block_sizes_wrapper()
{ std::auto_ptr<inner_class_value> G(get_complexredgp());
  G->value.fillCartan();
  std::auto_ptr<latmat_value>
  M(new latmat_value @|
    (latticetypes::LatticeMatrix(G->value.numRealForms()
                                ,G->value.numDualRealForms())
    ));
  for (size_t i = 0; i < M->value.numRows(); ++i)
    for (size_t j = 0; j < M->value.numColumns(); ++j)
      M->value(i,j) =
      G->value.blockSize(G->interface.in(i),G->dual_interface.in(j));
  push_value(M.release());
}

@ Next we shall provide a function that displays the occurrence of Cartan
matrices for various real forms. The rows will be indexed by real forms, and
the columns by Cartan matrices (note the alliteration).

@< Local function def...@>=
void occurrence_matrix_wrapper()
{ std::auto_ptr<inner_class_value> G(get_complexredgp());
  G->value.fillCartan();
  size_t nr=G->value.numRealForms();
  size_t nc=G->value.numCartanClasses();
  std::auto_ptr<latmat_value>
    M(new latmat_value(latticetypes::LatticeMatrix(nr,nc)));
  for (size_t i=0; i<nr; ++i)
  { bitmap::BitMap b=G->value.cartanSet(G->interface.in(i));
    for (size_t j=0; j<nc; ++j)
      M->value(i,j)= b.isMember(j) ? 1 : 0;
  }
  push_value(M.release());
}

@ We do the same for dual real forms. Note that we had to introduce the method
|dualCartanSet| for |complexredgp::ComplexReductiveGroup| in order to be able
to write this function.

@< Local function def...@>=
void dual_occurrence_matrix_wrapper()
{ std::auto_ptr<inner_class_value> G(get_complexredgp());
  G->value.fillCartan();
  size_t nr=G->value.numDualRealForms();
  size_t nc=G->value.numCartanClasses();
  std::auto_ptr<latmat_value>
    M(new latmat_value(latticetypes::LatticeMatrix(nr,nc)));
  for (size_t i=0; i<nr; ++i)
  { bitmap::BitMap b=G->value.dualCartanSet(G->dual_interface.in(i));
    for (size_t j=0; j<nc; ++j)
      M->value(i,j)= b.isMember(j) ? 1 : 0;
  }
  push_value(M.release());
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
the |value| field.

We shall later derive from this structure, without adding data members, a
structure to record dual real forms, which are constructed in the context of
the same inner class, but where the |value| field is constructed for the dual
group (and its dual inner class). In order to accommodate for that we provide
an additional protected constructor that will be defined later; this also
explains why the copy constructor is protected rather than private.

@< Type definitions @>=
struct real_form_value : public value_base
{ const inner_class_value parent;
  realredgp::RealReductiveGroup value;
@)
  real_form_value(inner_class_value p,realform::RealForm f)
  : parent(p), value(p.value,f) @+{}
  ~real_form_value() @+{} // everything is handled by destructor of |parent|
@)
  virtual void print(std::ostream& out) const;
  real_form_value* clone() const @+
    {@; return new real_form_value(*this); }
protected:
  real_form_value(inner_class_value,realform::RealForm,tags::DualTag);
  real_form_value(const real_form_value& v)
  : parent(v.parent), value(v.value) @+{}
};

@ When printing a real form, we give the name by which it was chosen (where
for once we do use the |parent| field), and provide some information about its
connectedness. Since the names of the real forms are indexed by their outer
number, but the real form itself stores its inner number, we must somewhat
laboriously make the conversion here.

@< Function def...@>=
void real_form_value::print(std::ostream& out) const
{ out << (value.isQuasisplit() ? "quasisplit " : "")
  << "real form '" @|
  << parent.interface.typeName(parent.interface.out(value.realForm())) @|
  << "', defining a"
  << (value.isConnected() ? " connected" : "" ) @|
  << (value.isSplit() ? " split" : "") @|
  << " real group" ;
}

@ To make a real form is easy, one provides a |inner_class_value| and a valid
index into its list of real forms. Since this number coming from the outside
is to be interpreted as an outer index, we must convert it to an inner index
at this point. As a special case we also provide the quasisplit form.

@< Local function def...@>=
void real_form_wrapper()
{ push_tuple_components();
  std::auto_ptr<int_value> i(get_int());
  std::auto_ptr<inner_class_value> G(get_complexredgp());
  if (i->value<0 || size_t(i->value)>=G->value.numRealForms())
    throw std::runtime_error ("illegal real form number: "+num(i->value));
  push_value(new real_form_value(*G,G->interface.in(i->value)));
}
@)
void quasisplit_form_wrapper()
{ std::auto_ptr<inner_class_value> G(get_complexredgp());
  push_value(new real_form_value(*G,G->value.quasisplit()));
}

@*2 Functions operating on real reductive groups.
We shall need a another variation of |get_int| to unpack our new values.

@< Declarations of exported functions @>=
real_form_value* get_real_form() throw(std::logic_error);

@~Again, there are no surprises.
@< Function definitions @>=
real_form_value* get_real_form() throw(std::logic_error)
{ value_ptr p=pop_arg();
  real_form_value* result=dynamic_cast<real_form_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not a real form"); }
@.Argument is not a real form@>
  return result;
}

@ Here is a function that gives full information about the component group of
a reals reductive group: it returns the rank~$r$ such that the component group
is isomorphic to $(\Z/2\Z)^r$.

@< Local function def...@>=
void components_rank_wrapper()
{ std::auto_ptr<real_form_value> R(get_real_form());
  const latticetypes::ComponentList c=R->value.componentReps();
  push_value(new int_value(c.size()));
}

@ And here is one that counts the number of Cartan classes for the real form.
Note that we must call |fillCartan()| here to ensure that the Cartan classes
are in fact generated.

@< Local function def...@>=
void count_Cartans_wrapper()
{ std::auto_ptr<real_form_value> R(get_real_form());
  R->value.fillCartan();
  push_value(new int_value(R->value.numCartan()));
}

@*2 Dual real forms.
Although they could be considered as real forms for the dual group and dual
inner class, it will be convenient to be able to construct dual real forms in
the context of the |inner_class_value| itself. The resulting objects will have
a different type than ordinary real forms, but the will use the same structure
layout. The interpretation of the fields is as follows for dual real forms:
the |parent| field contains a copy of the original |inner_class_value|, but
the |value| field contains a |realredgp::RealReductiveGroup| object
constructed for the dual inner class, so that its characteristics will be
correct.

@< Function def...@>=

real_form_value::real_form_value(inner_class_value p,realform::RealForm f
				,tags::DualTag)
: parent(p), value(p.dual,f)
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
private:
  dual_real_form_value(const dual_real_form_value& v)
  : real_form_value(v) @+ {}
};

@ The only real difference with real forms is a slightly modified output
routine.

@< Function def...@>=
void dual_real_form_value::print(std::ostream& out) const
{ out << (value.isQuasisplit() ? "quasisplit " : "")
  << "dual real form '" @|
  << parent.dual_interface.typeName
    (parent.dual_interface.out(value.realForm())) @|
  << "'";
}

@ To make a dual real form, one provides a |inner_class_value| and a valid
index into its list of dual real forms, which will be converted to an inner
index. We also provide the dual quasisplit form.

@< Local function def...@>=
void dual_real_form_wrapper()
{ push_tuple_components();
  std::auto_ptr<int_value> i(get_int());
  std::auto_ptr<inner_class_value> G(get_complexredgp());
  if (i->value<0 || size_t(i->value)>=G->value.numDualRealForms())
    throw std::runtime_error ("illegal dual real form number: "+num(i->value));
  push_value(new dual_real_form_value(*G,G->dual_interface.in(i->value)));
}
@)
void dual_quasisplit_form_wrapper()
{ std::auto_ptr<inner_class_value> G(get_complexredgp());
  push_value(new dual_real_form_value(*G,G->dual.quasisplit()));
}

@ And here is another obligatory part

@< Declarations of exported functions @>=
dual_real_form_value* get_dual_real_form() throw(std::logic_error);

@~Again, there are no surprises.
@< Function definitions @>=
dual_real_form_value* get_dual_real_form() throw(std::logic_error)
{ value_ptr p=pop_arg();
  dual_real_form_value* result=dynamic_cast<dual_real_form_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not a dual real form"); }
@.Argument is not a dual real form@>
  return result;
}

@ Rather than provide all functions for real forms for dual real forms, we
provide a function that converts it to a real form, which of course will be
associated to the dual inner class.

@< Local function def...@>=
void real_form_from_dual_wrapper()
{ std::auto_ptr<dual_real_form_value> d(get_dual_real_form());
  push_value(new real_form_value
                 (inner_class_value(d->parent,tags::DualTag())
                 ,d->value.realForm()));
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

@ The Cartan class is identified by a number within the inner class. This
number is actually its sequence number in the order in which Cartan classes
are generated in this inner class, so it may vary from one run to another;
nevertheless the number is important to record, since some information about
the Cartan class is only stored at the level of the inner class, and the
number is required to obtain this information.

@< Type definitions @>=
struct Cartan_class_value : public value_base
{ const inner_class_value parent;
  size_t number;
  const cartanclass::CartanClass& value;
@)
  Cartan_class_value(const inner_class_value& p,size_t cn);
  ~Cartan_class_value() @+{} // everything is handled by destructor of |parent|
@)
  virtual void print(std::ostream& out) const;
  Cartan_class_value* clone() const @+
    {@; return new Cartan_class_value(*this); }
private:
  Cartan_class_value(const Cartan_class_value& v)
  : parent(v.parent), number(v.number), value(v.value) @+{}
};

@ In the constructor we check that the Cartan class with the given number
currently exists. This check \emph{follows} the selection of the pointer to
the |CartanClass| object that initialises the |value| field because we have no
choice (a reference field cannot be set later); the danger is doing this is
limited since the selection itself will probably not crash the program, and
the nonsensical pointer so obtained will probably be stored away without
inspection and then be forgotten when an error is thrown. Moreover, the
function we shall provide to access Cartan classes will ensure that when
constructed, the Cartan class does actually exist.

@< Function def...@>=
Cartan_class_value::Cartan_class_value(const inner_class_value& p,size_t cn)
: parent(p),number(cn),value(p.value.cartan(cn))
{ if (cn>=p.value.numCartanClasses()) throw std::runtime_error
  (std::string("Cartan class number ")+num(cn)+" does not (currently) exist");
}

@ When printing a Cartan class, we show its number, whether it is the most
split one, and for how many real forms and dual real forms it is valid.

@< Function def...@>=
void Cartan_class_value::print(std::ostream& out) const
{ out << "Cartan class #" << number << ", occurring for " @|
  << value.numRealForms() << " real "
  << (value.numRealForms()==1 ? "form" : "forms") << " and for "@|
  << value.numDualRealForms() << " dual real "
  << (value.numDualRealForms()==1 ? "form" : "forms");
}

@ To make a Cartan class, one must provide a |real_form_value| together with a
valid index into its list of Cartan classes.

@< Local function def...@>=
void Cartan_class_wrapper()
{ push_tuple_components();
  std::auto_ptr<int_value> i(get_int());
  std::auto_ptr<real_form_value> rf(get_real_form());
  rf->value.fillCartan();
  if (i->value<0 || size_t(i->value)>=rf->value.numCartan())
    throw std::runtime_error
    ("illegal Cartan class number: "+num(i->value)
    +", this real form only has "+num(rf->value.numCartan())+" of them");
  bitmap::BitMap cs=rf->value.cartanSet();
  push_value(new Cartan_class_value(rf->parent,cs.n_th(i->value)));
}

@ We need yet another function for getting a value from the stack.

@< Declarations of exported functions @>=
Cartan_class_value* get_Cartan_class() throw(std::logic_error);

@~Again, there are no surprises.
@< Function definitions @>=
Cartan_class_value* get_Cartan_class() throw(std::logic_error)
{ value_ptr p=pop_arg();
  Cartan_class_value* result=dynamic_cast<Cartan_class_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not a Cartan class"); }
@.Argument is not a Cartan class@>
  return result;
}

@ This function and the following provide the functionality of the Atlas
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
{ using basic_io::operator<<;

@/std::auto_ptr<Cartan_class_value> cc(get_Cartan_class());

  const rootdata::RootDatum& rd=cc->parent.value.rootDatum();

  prettyprint::printTorusType(std::cout,cc->value.fiber().torus())
  << std::endl;

  std::cout << "twisted involution orbit size: " << cc->value.orbitSize()
   << std::endl;

@)// print type of imaginary root system
  lietype::LieType ilt,rlt,clt;
  rootdata::lieType(ilt,cc->value.simpleImaginary(),rd);

  if (ilt.size() == 0)
    std::cout << "imaginary root system is empty" << std::endl;
  else
    std::cout << "imaginary root system: " << ilt << std::endl;

@)// print type of real root system
  rootdata::lieType(rlt,cc->value.simpleReal(),rd);

  if (rlt.size() == 0)
    std::cout << "real root system is empty" << std::endl;
  else
    std::cout << "real root system: " << rlt << std::endl;

@)// print type of complex root system
  rootdata::lieType(clt,cc->value.simpleComplex(),rd);

  if (clt.size() == 0)
    std::cout << "complex factor is empty" << std::endl;
  else
    std::cout << "complex factor: " << clt << std::endl;
@)
  wrap_tuple(0);
}

@ For the fiber group partition information that was not produced by the
previous function, we use a Cartan class and a real form as parameters. This
function returns the part of the |weakReal| partition stored in the fiber of
the Cartan class that corresponds to the given (weak) real form. The numbering
of the parts of that partition are not the numbering of the real forms
themselves, so they must be translated through the |realFormLabels| list for
the Cartan class, which must be obtained from its |parent| inner class. If the
Cartan class does not exist for the given real form, then it will not occur in
that |realFormLabels| list, and the part returned here will be empty.

Currently the part of the partition is returned as a list of integral values.

@f pi NULL

@< Local function def...@>=
void fiber_part_wrapper()
{ push_tuple_components();
  std::auto_ptr<real_form_value> rf(get_real_form());
@/std::auto_ptr<Cartan_class_value> cc(get_Cartan_class());
  if (&rf->parent.value!=&cc->parent.value)
    throw std::runtime_error
    ("inner class mismatch between real form and Cartan class");
  bitmap::BitMap b(cc->parent.value.cartanSet(rf->value.realForm()));
  if (!b.isMember(cc->number))
    throw std::runtime_error
    ("fiber_part: inner class not defined for this real form");
@)
  const partition::Partition& pi = cc->value.fiber().weakReal();
  const realform::RealFormList rf_nr=
     cc->parent.value.realFormLabels(cc->number);
     // translate part number of |pi| to real form
  std::auto_ptr<row_value> result
    (new row_value(std::vector<value_ptr>()));
  for (size_t i=0; i<pi.size(); ++i)
    if ( rf_nr[pi(i)] == rf->value.realForm())
      result->value.push_back(new int_value(i));
  push_value(result.release());
}

@ The function |print_grading| gives on a per-real-form basis the
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
  std::auto_ptr<real_form_value> rf(get_real_form());
@/std::auto_ptr<Cartan_class_value> cc(get_Cartan_class());
  if (&rf->parent.value!=&cc->parent.value)
    throw std::runtime_error
    ("inner class mismatch between real form and Cartan class");
  bitmap::BitMap b(cc->parent.value.cartanSet(rf->value.realForm()));
  if (!b.isMember(cc->number))
    throw std::runtime_error
    ("fiber_part: inner class not defined for this real form");
@)
  const partition::Partition& pi = cc->value.fiber().weakReal();
  const realform::RealFormList rf_nr=
     cc->parent.value.realFormLabels(cc->number);
     // translate part number of |pi| to real form
@)
  const rootdata::RootList& si = cc->value.fiber().simpleImaginary();
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
{ rootdata::cartanMatrix(cm,si,cc->parent.value.rootDatum());
  dynkin::DynkinDiagram d(cm); dynkin::bourbaki(sigma,d);
}

@ The imaginary root system might well be empty, so we make special provisions
for this case.

@< Print information about the imaginary root system... @>=
{ std::cout << "Imaginary root system is ";
  if (si.size()==0) std::cout<<"empty.\n";
  else
  { using basic_io::operator<<;
    lietype::LieType t; dynkin::lieType(t,cm);

    std::cout << "of type " << t << ", with simple root"
              << (si.size()==1 ? " " : "s ");
    for (size_t i=0; i<si.size(); ++i)
      std::cout << si[sigma[i]] << (i<si.size()-1 ? "," : ".\n");
  }
}


@ The gradings are computed by the |grading| method of the fiber of our Cartan
class, and converted to a string of characters |'0'| and |'1'| by
|prettyprint::prettyPrint|. The permutation |sigma| was computed in order to
present these bits in an order adapted to the Dynkin diagram of the imaginary
root system. If that system is empty, the strings giving the gradings are
empty, but the part of the partition |pi| will not be, unless the Cartan class
does not exist for the given real form.

@< Print the gradings for the part of |pi|... @>=
{ bool first=true;
  for (size_t i=0; i<pi.size(); ++i)
    if ( rf_nr[pi(i)] == rf->value.realForm())
    { std::cout << ( first ? first=false,'[' : ',');
      gradings::Grading gr; cc->value.fiber().grading(gr,i);
      gr.permute(sigma);
      prettyprint::prettyPrint(std::cout,gr,si.size());
    }
  std::cout << "]" << std::endl;
}

@ Finally we install everything (where did we hear that being said before?)

@< Install wrapper functions @>=
install_function(real_form_wrapper,"real_form","(InnerClass,int->RealForm)");
install_function(quasisplit_form_wrapper,"quasisplit_form"
		,"(InnerClass->RealForm)");
install_function(components_rank_wrapper,"components_rank","(RealForm->int)");
install_function(count_Cartans_wrapper,"count_Cartans","(RealForm->int)");
install_function(dual_real_form_wrapper,"dual_real_form"
				       ,"(InnerClass,int->DualRealForm)");
install_function(dual_quasisplit_form_wrapper,"dual_quasisplit_form"
		,"(InnerClass->DualRealForm)");
install_function(real_form_from_dual_wrapper,"real_form_from_dual"
				  ,"(DualRealForm->RealForm)");
install_function(Cartan_class_wrapper,"Cartan_class"
		,"(RealForm,int->CartanClass)");
install_function(print_Cartan_info_wrapper,"print_Cartan_info"
		,"(CartanClass->)");
install_function(fiber_part_wrapper,"fiber_part"
		,"(CartanClass,RealForm->[int])");
install_function(print_gradings_wrapper,"print_gradings"
		,"(CartanClass,RealForm->)");





@* Index.

% Local IspellDict: default
