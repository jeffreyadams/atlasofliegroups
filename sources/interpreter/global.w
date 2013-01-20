% Copyright (C) 2012 Marc van Leeuwen
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
\def\foreign#1{{\sl#1\/}}
\def\id{\mathop{\rm id}}

@* Outline.
%
This file originated from splitting off a part of the module \.{evaluator.w}
that was getting too large. It collects functions that are important to the
evaluation process, but which are not part of the recursive machinery of
type-checking and evaluating all different expression forms.

This module has three major, largely unrelated, parts. The first part groups
some peripheral operations to the main evaluation process: the calling
interface to the type checking and conversion process, and operations that
implement global changes such as the introduction of new global identifiers.
The second part is dedicated to some fundamental types, like integers,
Booleans, strings, without which the programming language would be an empty
shell (but types more specialised to the Atlas software are defined in another
module, \.{built-in-types}. Finally there is a large section with basic
functions related to these types.

@( global.h @>=

#ifndef GLOBAL_H
#define GLOBAL_H

@< Includes needed in the header file @>@;
namespace atlas { namespace interpreter {
@< Type definitions @>@;
@< Declarations of global variables @>@;
@< Declarations of exported functions @>@;
@< Template and inline function definitions @>@;
}@; }@;
#endif

@ The implementation unit follows the usual pattern, although not all possible
sections are present.

@h "global.h"
@h <cstdlib>
@c
namespace atlas { namespace interpreter {
@< Global variable definitions @>@;
namespace {@;
@< Local function definitions @>@;
}@;
@< Global function definitions @>@;
}@; }@;

@* Global operations.
%
We start with several operations at the outer level of the interpreter such as
initialisations and the initial invocation of the type checker.

Before executing anything the evaluator needs some initialisation, called
from the main program.

@< Declarations of exported functions @>=
void initialise_evaluator();

@ The details of this initialisation will be given when the variables involved
are introduced.

@< Global function definitions @>=
void initialise_evaluator()
@+{@; @< Initialise evaluator @> }

@~Although not necessary, the following will avoid some early reallocations of
|execution_stack|, a vector variable defined in \.{types.w}.

@< Initialise evaluator @>=
execution_stack.reserve(16); // avoid some early reallocations

@*1 Invoking the type checker.
%
Let us recapitulate the general organisation of the evaluator, explained in
more detail in the introduction of \.{evaluator.w}. The parser reads what the
user types, and returns an |expr| value representing the abstract syntax tree.
Then the highly recursive function |convert_expr| is called for this value,
which will either produce (a pointer to) an executable object of a type
derived from |expression_base|, or throw an exception in case a type error or
other problem is detected. In the former case we are ready to call the
|evaluate| method of the value returned by |convert_expr|; after this main
program will print the result. The initial call to |convert_expr| however is
done via |analyse_types|, which takes care of catching any exceptions thrown,
and printing error messages.

@< Declarations of exported functions @>=
type_ptr analyse_types(const expr& e,expression_ptr& p)
   throw(std::bad_alloc,std::runtime_error);

@~The function |analyse_types| switches the roles of the output parameter
|type| of |convert_expr| and its return value: the former becomes the return
value and the latter is assigned to the output parameter~|p|. The initial
value of |type| passed to |convert_expr| is a completely unknown type. Since
we cannot return any type from |analyse_types| in the presence of errors, we
map these errors to |std::runtime_error| after printing their error message;
that error is an exception for which the code that calls us will have to
provide a handler anyway, and which handler will serve as a more practical
point to really resume after an error.

@h "evaluator.h"

@< Global function definitions @>=
type_ptr analyse_types(const expr& e,expression_ptr& p)
  throw(std::bad_alloc,std::runtime_error)
{ try
  {@; type_ptr type=copy(unknown_type);
    p.reset(convert_expr(e,*type));
    return type;
  }
  catch (type_error& err)
  { std::cerr << err.what() << ":\n  Subexpression " << err.offender << @|
    " has wrong type: found " << *err.actual << @|
    " while " << *err.required << " was needed.\n";
@.Subexpression has wrong type@>
  }
  catch (expr_error& err)
  { std::cerr << "Error in expression " << err.offender << "\n  "
              << err.what() << std::endl;
  }
  catch (program_error& err)
  {@; std::cerr << "Error during analysis of expression " << e << "\n  "
                << err.what() << std::endl;
  }
  throw std::runtime_error("Type check failed");
@.Type check failed@>
}


@*1 Operations other than evaluation of expressions.
%
This section will be devoted to some interactions between user and program
that do not consist just of evaluating expressions.

The function |global_set_identifier| handles introducing identifiers, either
normal ones or overloaded instances of functions, using the \&{set} syntax.
The function |global_declare_identifier| just introduces an identifier into
the global (non-overloaded) table with a definite type, but does not provide a
value (it will then have to be assigned to before it can be validly used).
Conversely |global_forget_identifier| removes an identifier from the global
table. The last three functions serve to provide the user with information
about the state of the global tables (but invoking |type_of_expr| will
actually go through the full type analysis and conversion process).

@< Declarations of exported functions @>=
void global_set_identifier(struct id_pat id, expr e, int overload);
void global_declare_identifier(id_type id, type_p type);
void global_forget_identifier(id_type id);
void global_forget_overload(id_type id, type_p type);
void show_ids();
void type_of_expr(expr e);
void show_overloads(id_type id);

@ Global identifiers can be introduced (or modified) by the function
|global_set_identifier|. We allow (in the \&{set} syntax) the same
possibilities in a global identifier definition as in a local one (the \&{let}
syntax), so we take an |id_pat| as argument. We also handle definitions of
overloaded function instances in case the parameter |overload| is nonzero. In
a change from our initial implementation, that parameter set by the parser
allows overloading but does not force it. Allowing the parameter to be cleared
here then actually serves to allow more cases to be handled using the overload
table, since the parser will now set |overload>0| more freely. Indeed the
parser currently passes |overload==0| only when the
``\\{identifier}\.:\\{value}'' syntax is used to introduce a new identifier.

However, the code below sets |overload=0| also whenever the defining
expression has anything other that a function type (which must in addition
take at least one argument); in particular this makes it impossible to add
multiple items to the overload table with a single \&{set} command. This
restriction is mostly motivated by the complications that allowing mixing of
overloaded and non-overloaded definitions would entail (for instance when
reporting the resulting updates to the user), and by the fact that operator
definitions would in any case be excluded from multiple-overload situations
for syntactic reasons.

We follow the logic for type-analysis of a let-expression, and for evaluation
we follow the logic of binding identifiers in a user-defined function (these
are defined in \.{evaluator.w}). However we use |analyse_types| here (which
catches and reports errors) rather than calling |convert_expr| directly. To
provide some feedback to the user we report any types assigned, but not the
values.

@< Global function definitions @>=
void global_set_identifier(id_pat pat, expr rhs, int overload)
{ size_t n_id=count_identifiers(pat);
  static const char* phase_name[3] = {"type_check","evaluation","definition"};
  int phase=0; // needs to be declared outside the |try|, is used in |catch|
  try
  { expression_ptr e;
    type_ptr t=analyse_types(rhs,e);
    if (not pattern_type(pat)->specialise(*t))
      @< Report that type of |rhs| does not have required structure,
         and |throw| @>
    if (overload!=0)
      @< Set |overload=0| if type |t| is not an appropriate function type @>
@)
    phase=1;
    bindings b(n_id);
    thread_bindings(pat,*t,b); // match identifiers and their future types

    std::vector<shared_value> v;
    v.reserve(n_id);
@/  e->eval();
    thread_components(pat,pop_value(),v); // associate values with identifiers
@)
    phase=2;
    @< Emit indentation corresponding to the input level to |std::cout| @>
    if (overload==0)
      @< Add instance of identifiers in |b| with values in |v| to
         |global_id_table| @>
    else
      @< Add instance of identifier in |b[0]| with value in |v[0]| to
         |global_overload_table| @>

    std::cout << std::endl;
  }
  @< Catch block for errors thrown during a global identifier definition @>
}

@ When |overload>0|, choosing whether the definition enters into the overload
table or into the global identifier table is determined by the type of the
defining expression (in particular this allows operators to be defined by an
arbitrary expression). However, this creates the possibility (if the defining
expression should have non-function type) of causing an operator to be added
the global identifier table, which is pointless (since the syntax does not
allow such a value to be retrieved). Therefore this case needs some attention:
the parser will pass |overload==2| in this case, signalling that it must not
be cleared to~$0$, but rather result in an error message in cases where
setting it to~$0$ is attempted.

@< Set |overload=0| if type |t| is not an appropriate function type @>=
{ bool clear = t->kind!=function_type;
    // cannot overload with a non-function value
  if (not clear)
  { type_expr& arg=t->func->arg_type;
  @/clear = arg.kind==tuple_type and arg.tuple==NULL;
     // nor parameterless functions
  }
  if (clear and overload==2) // inappropriate function type with operator
  { std::string which(t->kind==function_type ? "parameterless " : "non-");
    throw std::runtime_error
      ("Cannot set operator to a "+which+"function value");
  }
  if (clear)
    overload=0;
}

@ For identifier definitions we print their names and types (paying attention
to the very common singular case), before calling |global_id_table->add|.
@< Add instance of identifiers in |b| with values in |v| to
   |global_id_table| @>=
{ if (n_id>0)
    std::cout << "Identifier";
  for (size_t i=0; i<n_id; ++i)
  { std::cout << (i==0 ? n_id==1 ? " " : "s " : ", ") @|
              << main_hash_table->name_of(b[i].first);
    if (global_id_table->type_of(b[i].first)!=NULL)
      std::cout << " (overriding previous)";
    std::cout << ": " << *b[i].second;
    global_id_table->add(b[i].first,v[i],copy(*b[i].second));
  }
}

@ For overloaded definitions the main difference is calling the |add| method
of |global_overload_table| instead of that of |global_id_table|, and the
different wording of the report to the user. However another difference is
that here the |add| method may throw because of a conflict of a new definition
with an existing one; we therefore do not print anything before the |add|
method has successfully completed. Multiple overloaded definitions in a
single \&{set} statement are excluded by the fact that a tuple type for the
defining expression will cause setting |overload=0|, so that we are prevented
from coming here either by that assignment or by a pattern mismatch error
being thrown during type analysis; the |assert| statement below checks this
exclusion. If the logic above were relaxed to allow such multiple definitions,
one could introduce a loop below as in the ordinary definition case; then
however error handling would also need adaptation, since a failed definition
need no be the first one, and the previous ones would need to be either undone
or not reported as failed.

@< Add instance of identifier in |b[0]| with value in |v[0]| to
   |global_overload_table| @>=
{ assert(n_id=1);
  size_t old_n=global_overload_table->variants(b[0].first).size();
  global_overload_table->add(b[0].first,v[0],copy(*b[0].second));
    // insert or replace table entry
  size_t n=global_overload_table->variants(b[0].first).size();
  if (n==old_n)
    std::cout << "Redefined ";
  else if (n==1)
    std::cout << "Defined ";
  else
    std::cout << "Added definition [" << n << "] of ";
  std::cout << main_hash_table->name_of(b[0].first) << ": " << *b[0].second;
}

@ For readability of the output produced during input from auxiliary files, we
emit two spaces for every current input level. The required information is
available from the |main_input_buffer|.

@< Emit indentation corresponding to the input level to |std::cout| @>=
{ unsigned int input_level = main_input_buffer->include_depth();
  std::cout << std::setw(2*input_level) << "";
}

@ When the right hand side type does not match the requested pattern, we throw
a |runtime_error| signalling this fact; we have to re-generate the required
pattern using |pattern_type| to do this.

@< Report that type of |rhs| does not have required structure, and |throw| @>=
{ std::ostringstream o;
  o << "Type " << *t @|
    << " of right hand side does not match required pattern "
    << *pattern_type(pat);
  throw std::runtime_error(o.str());
}

@ A |std::runtime_error| may be thrown either during type check, matching with
the identifier pattern, or evaluation; we catch all those cases here. Whether
or not an error is caught, the pattern |pat| and the expression |rhs| should
not be destroyed here, since the parser which aborts after calling this
function should do that while clearing its parsing stack.

@< Catch block for errors thrown during a global identifier definition @>=
catch (std::runtime_error& err)
{ std::cerr << err.what() << '\n';
  if (n_id>0)
  { std::vector<Hash_table::id_type> names; names.reserve(n_id);
    list_identifiers(pat,names);
    std::cerr << "  Identifier" << (n_id==1 ? "" : "s");
    for (size_t i=0; i<n_id; ++i)
      std::cerr << (i==0 ? " " : ", ") << main_hash_table->name_of(names[i]);
    std::cerr << " not " << (overload==0 ? "created." : "overloaded.")
              << std::endl;
  }
  reset_evaluator(); main_input_buffer->close_includes();
}
catch (std::logic_error& err)
{ std::cerr << "Unexpected error: " << err.what() << ", " @|
            << phase_name[phase]
            << " aborted.\n";
@/reset_evaluator(); main_input_buffer->close_includes();
}
catch (std::exception& err)
{ std::cerr << err.what() << ", "
            << phase_name[phase]
            << " aborted.\n";
@/reset_evaluator(); main_input_buffer->close_includes();
}

@ The following function is called when an identifier is declared with type
but undefined value.

@< Global function definitions @>=
void global_declare_identifier(Hash_table::id_type id, type_p t)
{ value undef=NULL;
  const type_expr& type=*t;
  global_id_table->add(id,shared_value(undef),copy(type));
  std::cout << "Identifier " << main_hash_table->name_of(id)
            << " : " << type << std::endl;
}

@ Finally the user may wish to forget the value of an identifier, which the
following function achieves.

@< Global function definitions @>=
void global_forget_identifier(Hash_table::id_type id)
{ std::cout << "Identifier " << main_hash_table->name_of(id)
            << (global_id_table->remove(id) ? " forgotten" : " not known")
            << std::endl;
}

@ Forgetting the binding of an overloaded identifier at a given type is
similar.

@< Global function definitions @>=
void global_forget_overload(Hash_table::id_type id, type_p t)
{ const type_expr& type=*t;
  std::cout << "Definition of " << main_hash_table->name_of(id)
            << '@@' << type @|
            << (global_overload_table->remove(id,type)
               ? " forgotten"
               : " not known")
            << std::endl;
}

@ It is useful to print type information, either for a single expression or
for all identifiers in the table. We here define and export the pointer
variable that is used for all normal output from the evaluator.

@< Declarations of global variables @>=
extern std::ostream* output_stream;

@ The |output_stream| will normally point to |std::cout|, but the pointer
may be assigned to in the main program, which causes output redirection.

@< Global variable definitions @>=
std::ostream* output_stream= &std::cout;

@ The function |type_of_expr| prints the type of a single expression, without
evaluating it. Since we allow arbitrary expressions, we must cater for the
possibility of failing type analysis, in which case |analyse_types|, after
catching it, will re-throw a |std::runtime_error|. By in fact catching and
reporting any |std::exception| that may be thrown, we also ensure ourselves
against unlikely events like |bad_alloc|.

@< Global function definitions @>=
void type_of_expr(expr e)
{ try
  {@; expression_ptr p;
    *output_stream << "type: " << *analyse_types(e,p) << std::endl;
  }
  catch (std::exception& err) {@; std::cerr<<err.what()<<std::endl; }
}

@ The function |show_overloads| has a similar purpose to |type_of_expr|,
namely to find out the types of overloaded symbols. It does not however need
to call |analyse_types|, as it just has to look into the overload table and
extract the types stored there.

@< Global function definitions @>=
void show_overloads(id_type id)
{ const overload_table::variant_list& variants =
   global_overload_table->variants(id);
  *output_stream
   << (variants.empty() ? "No overloads for " : "Overloaded instances of ") @|
   << main_hash_table->name_of(id) << std::endl;
 for (size_t i=0; i<variants.size(); ++i)
   *output_stream << "  "
    << variants[i].type->arg_type << "->" << variants[i].type->result_type @|
    << std::endl;
}

@ The function |show_ids| prints a table of all known identifiers, their
types, and the stored values; it does this both for overloaded symbols and for
global identifiers.

@< Global function definitions @>=
void show_ids()
{ *output_stream << "Overloaded operators and functions:\n"
                 << *global_overload_table @|
                 << "Global variables:\n" << *global_id_table;
}

@ Here is a tiny bit of global state that can be set from the main program and
inspected by any module that cares to (and that reads \.{global.h}).

@< Declarations of global variables @>=
extern int verbosity;

@~By raising the value of |verbosity|, some trace of internal operations can
be activated.

@< Global variable definitions @>=
int verbosity=0;

@ We shall define a small template function |str| to help giving sensible
error messages, by converting integers into strings that can be incorporated
into thrown values. The reason for using a template is that without them it is
hard to do a decent job for both signed and unsigned types.

@< Includes needed in the header file @>=
#include <string>
#include <sstream>

@~Using string streams, the definition of~|str| is trivial; in fact the
overloads of the output operator ``|<<|'' determine the exact conversion. As a
consequence of this implementation, the template function will in fact turn
anything printable into a string.

@< Template and inline function definitions @>=
template <typename T>
  std::string str(T n) @+{@; std::ostringstream s; s<<n; return s.str(); }

@* Basic types.
%
This section is devoted to primitive types that are not not very
Atlas-specific, ranging from integers to matrices, and which often have some
related functionality in the programming language (like conditional clauses
for Boolean values), which functionality is defined in \.{evaluator.w}. There
are also implicit conversions related to these types, and these will be
defined in the current module.

This section can be seen as in introduction to the large
module \.{built-in-types.w}, in which many more types and functions are
defined that provide Atlas-specific functionality. In fact the type for
rational numbers defined here is based on the |RatWeight| class defined in the
Atlas library, so we must include a header file (which defines the necessary
class template) into ours.

@<Includes needed in the header file @>=
#include "arithmetic.h"

@*1 First primitive types: integer, rational, string and Boolean values.
%
We derive the first ``primitive'' value types. For each type we define a
corresponding auto-pointer type (whose name has the \&{\_ptr} suffix),
since we shall often need to hold such values by pointers, and the risk of
exceptions is ever present. Whenever the values get stored in a more permanent
place (which includes being bound to an identifier, but also being places on
the |execution_stack|), we shall in fact use smart pointers of the
|shared_ptr| kind, and for the corresponding types we shall used
the \&{shared\_} prefix.

@< Type definitions @>=

struct int_value : public value_base
{ int val;
@)
  explicit int_value(int v) : val(v) @+ {}
  ~int_value()@+ {}
  void print(std::ostream& out) const @+{@; out << val; }
  int_value* clone() const @+{@; return new int_value(*this); }
  static const char* name() @+{@; return "integer"; }
private:
  int_value(const int_value& v) : val(v.val) @+{}
};
@)
typedef std::auto_ptr<int_value> int_ptr;
typedef std::tr1::shared_ptr<int_value> shared_int;
@)
struct rat_value : public value_base
{ Rational val;
@)
  explicit rat_value(Rational v) : val(v) @+ {}
  ~rat_value()@+ {}
  void print(std::ostream& out) const @+{@; out << val; }
  rat_value* clone() const @+{@; return new rat_value(*this); }
  static const char* name() @+{@; return "integer"; }
private:
  rat_value(const rat_value& v) : val(v.val) @+{}
};
@)
typedef std::auto_ptr<rat_value> rat_ptr;
typedef std::tr1::shared_ptr<rat_value> shared_rat;

@ Here are two more; this is quite repetitive.

@< Type definitions @>=

struct string_value : public value_base
{ std::string val;
@)
  explicit string_value(const std::string& s) : val(s) @+ {}
  ~string_value()@+ {}
  void print(std::ostream& out) const @+{@; out << '"' << val << '"'; }
  string_value* clone() const @+{@; return new string_value(*this); }
  static const char* name() @+{@; return "string"; }
private:
  string_value(const string_value& v) : val(v.val) @+{}
};
@)
typedef std::auto_ptr<string_value> string_ptr;
typedef std::tr1::shared_ptr<string_value> shared_string;
@)

struct bool_value : public value_base
{ bool val;
@)
  explicit bool_value(bool v) : val(v) @+ {}
  ~bool_value()@+ {}
  void print(std::ostream& out) const @+{@; out << std::boolalpha << val; }
  bool_value* clone() const @+{@; return new bool_value(*this); }
  static const char* name() @+{@; return "Boolean"; }
private:
  bool_value(const bool_value& v) : val(v.val) @+{}
};
@)
typedef std::auto_ptr<bool_value> bool_ptr;
typedef std::tr1::shared_ptr<bool_value> shared_bool;

@*1 Primitive types for vectors and matrices.
%
The interpreter distinguishes its own types like \.{[int]} ``row of integer''
from similar built-in types of the library, like \.{vec} ``vector'', which it
will consider to be primitive types. In fact a value of type ``vector''
represents an object of the Atlas type |int_Vector|, and similarly other
primitive types will stand for other Atlas types. We prefer using a basic type
name rather than something resembling the more mathematically charged
equivalents like |Weight| for |int_Vector|, as that might be more confusing
that helpful to users. In any case, the interpretation of the values is not at
all fixed (vectors are used for coweights and (co)roots as well as for
weights, and matrices could denote either a basis or an automorphism of a
lattice). The header \.{atlas\_types.h} makes sure all types are pre-declared,
but we need to see the actual type definitions in order to incorporated these
values in ours.

@< Includes needed in the header file @>=
#include "atlas_types.h" // type declarations that are ``common knowledge''
#include "types.h" // base types we derive from
#include "matrix.h" // to make |int_Vector| and |int_Matrix| complete types
#include "ratvec.h" // to make |RatWeight| a complete type

@ The definition of |vector_value| is much like the preceding ones. In its
constructor, the argument is a reference to |std::vector<int>|, from which
|int_Vector| is derived (without adding data members); since a constructor for
the latter from the former is defined, we can do with just one constructor for
|vector_value|.

@< Type definitions @>=

struct vector_value : public value_base
{ int_Vector val;
@)
  explicit vector_value(const std::vector<int>& v) : val(v) @+ {}
  ~vector_value()@+ {}
  virtual void print(std::ostream& out) const;
  vector_value* clone() const @+{@; return new vector_value(*this); }
  static const char* name() @+{@; return "vector"; }
private:
  vector_value(const vector_value& v) : val(v.val) @+{}
};
@)
typedef std::auto_ptr<vector_value> vector_ptr;
typedef std::tr1::shared_ptr<vector_value> shared_vector;

@ Matrices and rational vectors follow the same pattern, but in this case the
constructors take a constant reference to a type identical to the one that
will be stored.

@< Type definitions @>=
struct matrix_value : public value_base
{ matrix::Matrix<int> val;
@)
  explicit matrix_value(const matrix::Matrix<int>& v) : val(v) @+ {}
  ~matrix_value()@+ {}
  virtual void print(std::ostream& out) const;
  matrix_value* clone() const @+{@; return new matrix_value(*this); }
  static const char* name() @+{@; return "matrix"; }
private:
  matrix_value(const matrix_value& v) : val(v.val) @+{}
};
@)
typedef std::auto_ptr<matrix_value> matrix_ptr;
typedef std::tr1::shared_ptr<matrix_value> shared_matrix;
@)
struct rational_vector_value : public value_base
{ RatWeight val;
@)
  explicit rational_vector_value(const RatWeight& v):val(v)@+{}
  rational_vector_value(const int_Vector& v,int d)
   : val(v,d) @+ {@; val.normalize(); }
  ~rational_vector_value()@+ {}
  virtual void print(std::ostream& out) const;
  rational_vector_value* clone() const
   @+{@; return new rational_vector_value(*this); }
  static const char* name() @+{@; return "rational vector"; }
private:
  rational_vector_value(const rational_vector_value& v) : val(v.val) @+{}
};
@)
typedef std::auto_ptr<rational_vector_value> rational_vector_ptr;
typedef std::tr1::shared_ptr<rational_vector_value> shared_rational_vector;
@)

@ To make a small but visible difference in printing between vectors and lists
of integers, weights will be printed in equal width fields one longer than the
minimum necessary. Rational vectors are a small veriation.

@h<sstream>
@h<iomanip>

@< Global function def... @>=
void vector_value::print(std::ostream& out) const
{ size_t l=val.size(),w=0; std::vector<std::string> tmp(l);
  for (size_t i=0; i<l; ++i)
  { std::ostringstream s; s<<val[i]; tmp[i]=s.str();
    if (tmp[i].length()>w) w=tmp[i].length();
  }
  if (l==0) out << "[ ]";
  else
  { w+=1; out << std::right << '[';
    for (size_t i=0; i<l; ++i)
      out << std::setw(w) << tmp[i] << (i<l-1 ? "," : " ]");
  }
}
@)
void rational_vector_value::print(std::ostream& out) const
{ size_t l=val.size(),w=0; std::vector<std::string> tmp(l);
  for (size_t i=0; i<l; ++i)
  { std::ostringstream s; s<<val.numerator()[i]; tmp[i]=s.str();
    if (tmp[i].length()>w) w=tmp[i].length();
  }
  if (l==0) out << "[ ]";
  else
  { w+=1; out << std::right << '[';
    for (size_t i=0; i<l; ++i)
      out << std::setw(w) << tmp[i] << (i<l-1 ? "," : " ]");
  }
  out << '/' << val.denominator();
}

@ For matrices we align columns, and print vertical bars along the sides.
However if there are no entries, we print the dimensions of the matrix.

@< Global function def... @>=
void matrix_value::print(std::ostream& out) const
{ size_t k=val.numRows(),l=val.numColumns();
  if (k==0 or l==0)
  {@;  out << "The " << k << 'x' << l << " matrix"; return; }
  std::vector<size_t> w(l,0);
  for (size_t i=0; i<k; ++i)
    for (size_t j=0; j<l; ++j)
    { std::ostringstream s; s<<val(i,j); size_t len=s.str().length();
      if (len>w[j]) w[j]=len;
    }
  out << std::endl << std::right;
  for (size_t i=0; i<k; ++i)
  { out << '|';
    for (size_t j=0; j<l; ++j)
      out << std::setw(w[j]+1) << val(i,j) << (j<l-1 ? ',' : ' ');
    out << '|' << std::endl;
  }
}

@*1 Implementing some conversion functions.
%
Here we define the set of implicit conversions that apply to types in the base
language; there will also be implicit conversions added that concern
Atlas-specific types, such as the conversion from (certain) strings to Lie
types. The conversions defined here are from integers to rationals,
from lists of rationals to rational vectors and back, from vectors to rational
vectors, and from (nested) lists of integers to vectors and matrices and back.

@ Here are the first conversions involving rational numbers and rational
vectors.

@< Local function def... @>=
void rational_convert() // convert integer to rational (with denominator~1)
{@; shared_int i = get<int_value>();
    push_value(new rat_value(Rational(i->val)));
}
@)
void ratlist_ratvec_convert() // convert list of rationals to rational vector
{ shared_row r = get <row_value>();
  int_Vector numer(r->val.size()),denom(r->val.size());
  unsigned int d=1;
  for (size_t i=0; i<r->val.size(); ++i)
  { Rational frac = force<rat_value>(r->val[i].get())->val;
    numer[i]=frac.numerator();
    denom[i]=frac.denominator();
    d=arithmetic::lcm(d,denom[i]);
  }
  for (size_t i=0; i<r->val.size(); ++i)
    numer[i]*= d/denom[i]; // adjust numerators to common denominator

  push_value(new rational_vector_value(numer,d)); // normalises
}
@)
void ratvec_ratlist_convert() // convert rational vector to list of rationals
{ shared_rational_vector rv = get<rational_vector_value>();
  row_ptr result(new row_value(rv->val.size()));
  for (size_t i=0; i<rv->val.size(); ++i)
  { Rational q(rv->val.numerator()[i],rv->val.denominator());
    result->val[i] = shared_value(new rat_value(q.normalize()));
  }
  push_value(result);
}
@)
void vec_ratvec_convert() // convert vector to rational vector
{ shared_vector v = get<vector_value>();
  push_value(new rational_vector_value(RatWeight(v->val,1)));
}

@ The conversions into vectors or matrices use an auxiliary function
|row_to_weight|, which constructs a new |Weight| from a row of integers,
leaving the task to clean up that row to their caller.

Note that |row_to_weight| returns its result by value (rather than by
assignment to a reference parameter); in principle this involves making a copy
of the vector |result| (which includes duplicating its entries). However since
there is only one |return| statement, the code generated by a decent compiler
like \.{g++} creates the vector object directly at its destination in this
case (whose location is passed as a hidden pointer argument), and does not
call the copy constructor at all. More generally this is avoided if all
|return| statements return the same variable, or if all of them return
expressions instead. So we shall not hesitate to use this idiom henceforth; in
fact, while it was exceptional when this code was first written, it is now
being used throughout the Atlas library. In fact this is so much true that
this comment can be considered to be a mere relic to remind where this change
of idiom was first applied within the Atlas software.

@< Local function def... @>=
int_Vector row_to_weight(const row_value& r)
{ int_Vector result(r.val.size());
  for(size_t i=0; i<r.val.size(); ++i)
    result[i]=force<int_value>(r.val[i].get())->val;
  return result;
}
@)
void intlist_vector_convert()
{@; shared_row r(get<row_value>());
  push_value(new vector_value(row_to_weight(*r)));
}
@)
void intlist_ratvec_convert()
{@; shared_row r(get<row_value>());
  push_value(new rational_vector_value(RatWeight(row_to_weight(*r),1)));
}

@ The conversion |veclist_matrix_convert| is longer, but still
straightforward. Any ragged columns are silently extended with null entries to
make a rectangular shape for the matrix. This has as inevitable consequence
that for an empty list of vectors the number of rows of the resulting empty
matrix is based on an uninformed and usually wrong guess (namely that it is
$0$, like the number of columns. Functions that need a set of vectors as
argument can better take matrix as argument, so that the size of a vector can
be deduced even if there aren't any, but if the user nevertheless supplies a
list of vectors using this coercion, then the vector size will be wrongly
deduced for an empty list. For this reason functions that really need the
vector size do better to include it as a separate, though often redundant,
argument.

@< Local function def... @>=
void veclist_matrix_convert()
{ shared_row r(get<row_value>());
@/std::vector<int_Vector > column_list;
  column_list.reserve(r->val.size());
  size_t depth=0; // maximal length of vectors
  for(size_t i=0; i<r->val.size(); ++i)
  { column_list.push_back(force<vector_value>(r->val[i].get())->val);
    if (column_list[i].size()>depth) depth=column_list[i].size();
  }
  for(size_t i=0; i<column_list.size(); ++i)
    // extend weights with 0's if necessary
    if (column_list[i].size()<depth)
    { size_t j=column_list[i].size();
      column_list[i].resize(depth);
      for (;j<depth; ++j) column_list[i][j]=0;
    }
  push_value(new matrix_value(matrix::Matrix<int>(column_list,depth)));
}

@ There remains one ``internalising'' conversion function, from row of row of
integer to matrix.

@< Local function def... @>=
void intlistlist_matrix_convert()
{ shared_row r(get<row_value>());
@/std::vector<int_Vector > column_list;
  column_list.reserve(r->val.size());
  size_t depth=0; // maximal length of vectors
  for(size_t i=0; i<r->val.size(); ++i)
  { column_list.push_back(row_to_weight(*force<row_value>(r->val[i].get())));
    if (column_list[i].size()>depth) depth=column_list[i].size();
  }
  for(size_t i=0; i<column_list.size(); ++i)
    // extend weights with 0's if necessary
    if (column_list[i].size()<depth)
    { size_t j=column_list[i].size();
      column_list[i].resize(depth);
      for (;j<depth; ++j) column_list[i][j]=0;
    }
  push_value(new matrix_value(matrix::Matrix<int>(column_list,depth)));

}

@ There remain the ``externalising'' conversions (towards lists of values) of
the vector and matrix conversions. It will be handy to have a basic function
|weight_to_row| that performs more or less the inverse transformation of
|row_to_weight|, but rather than returning a |row_value| it returns a
|row_ptr| pointing to it.

@< Local function def... @>=
row_ptr weight_to_row(const int_Vector& v)
{ row_ptr result (new row_value(v.size()));
  for(size_t i=0; i<v.size(); ++i)
    result->val[i]=shared_value(new int_value(v[i]));
  return result;
}
@)
void vector_intlist_convert()
{@; shared_vector v(get<vector_value>());
  push_value(weight_to_row(v->val));
}
@)
void matrix_veclist_convert()
{ shared_matrix m=get<matrix_value>();
  row_ptr result(new row_value(m->val.numColumns()));
  for(size_t i=0; i<m->val.numColumns(); ++i)
    result->val[i]=shared_value(new vector_value(m->val.column(i)));
  push_value(result);
}
@)
void matrix_intlistlist_convert()
{ shared_matrix m=get<matrix_value>();
  row_ptr result(new row_value(m->val.numColumns()));
  for(size_t i=0; i<m->val.numColumns(); ++i)
    result->val[i]=shared_value(weight_to_row(m->val.column(i)).release());

  push_value(result);
}

@ All that remains is to initialise the |coerce_table|.
@< Initialise evaluator @>=
coercion(int_type,rat_type, "QI", rational_convert); @/
coercion(row_of_rat_type,ratvec_type, "Qv[Q]", ratlist_ratvec_convert); @/
coercion(ratvec_type,row_of_rat_type, "[Q]Qv", ratvec_ratlist_convert); @/
coercion(vec_type,ratvec_type,"QvV", vec_ratvec_convert); @/
coercion(row_of_int_type,ratvec_type,"Rv[I]", intlist_ratvec_convert);
@)
coercion(row_of_int_type, vec_type, "V[I]", intlist_vector_convert); @/
coercion(row_of_vec_type,mat_type, "M[V]", veclist_matrix_convert);
coercion(row_row_of_int_type,mat_type, "M[[I]]", intlistlist_matrix_convert); @/
coercion(vec_type,row_of_int_type, "[I]V", vector_intlist_convert); @/
coercion(mat_type,row_of_vec_type, "[V]M", matrix_veclist_convert); @/
coercion(mat_type,row_row_of_int_type, "[[I]]M", matrix_intlistlist_convert); @/

@* Wrapper functions.
%
We now come to defining wrapper functions. The arguments and results of
wrapper functions will be transferred from and to stack as a |shared_value|,
so a wrapper function has neither arguments nor a result type. Variables that
refer to a wrapper function have the type |wrapper_function| defined below;
the |level| parameter serves the same function as for |evaluate| methods of
classes derived from |expression_base|, as described in \.{evaluator.w}: to
inform whether a result value should be produced at all, and if so whether it
should be expanded on the |execution_stack| in case it is a tuple.

@< Type definitions @>=
typedef void (* wrapper_function)(expression_base::level);

@ The following function will greatly
facilitate the later repetitive task of installing wrapper functions.

@< Declarations of exported functions @>=
void install_function
 (wrapper_function f,const char*name, const char* type_string);

@ We start by determining the specified type, and building a print-name for
the function that appends the argument type (since there will potentially be
many instances with the same name). Then we construct a |builtin_value| object
and finally add it to |global_overload_table|. Although currently there are no
built-in functions with void argument type, we make a provision for them in
case they would be needed later; notably they should not be overloaded and are
added to |global_id_table| instead.

@< Global function def... @>=
void install_function
 (wrapper_function f,const char*name, const char* type_string)
{ type_ptr type = make_type(type_string);
  std::ostringstream print_name; print_name<<name;
  if (type->kind!=function_type)
    throw std::logic_error
     ("Built-in with non-function type: "+print_name.str());
  if (type->func->arg_type==void_type)
  { shared_value val(new builtin_value(f,print_name.str()));
    global_id_table->add (main_hash_table->match_literal(name),val,type);
  }
  else
  { print_name << '@@' << type->func->arg_type;
    shared_value val(new builtin_value(f,print_name.str()));
    global_overload_table->add(main_hash_table->match_literal(name),val,type);
  }
}

@*1 Integer functions.
%
Our first built-in functions implement integer arithmetic. Arithmetic
operators are implemented by wrapper functions with two integer arguments.
Since arguments to built-in functions are evaluated with |level| parameter
|multi_value|, two separate values will be present on the stack. Note that
these are pulled from the stack in reverse order, which is important for the
non-commutative operations like `|-|' and `|/|'. Since values are shared, we
must allocate new value objects for the results.

@< Local function definitions @>=

void plus_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new int_value(i+j));
}
@)
void minus_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new int_value(i-j));
}
@)
void times_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new int_value(i*j));
}

@ We take the occasion of defining a division operation to repair the integer
division operation |operator/| built into \Cpp, which is traditionally broken
for negative dividends. This is done by using |arithmetic::divide| that
handles such cases correctly (rounding the quotient systematically downwards).
Since |arithmetic::divide| takes an unsigned second argument, we handle the
case of a negative divisor ourselves. Incidentally, this Euclidean division
operation will be bound to the operator ``$\backslash$'', because ``$/$'' is
used to form rational numbers.

@< Local function definitions @>=
void divide_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (j==0) throw std::runtime_error("Division by zero");
  if (l!=expression_base::no_value)
    push_value(new int_value
     (j>0 ? arithmetic::divide(i,j) : -arithmetic::divide(i,-j)));
}

@ We also define a remainder operation |modulo|, a combined
quotient-and-remainder operation |divmod|, unary subtraction, and an integer
power operation (defined whenever the result is integer).

@< Local function definitions @>=
void modulo_wrapper(expression_base::level l)
{ int  j=get<int_value>()->val; int i=get<int_value>()->val;
  if (j==0) throw std::runtime_error("Modulo zero");
  if (l!=expression_base::no_value)
    push_value(new int_value(arithmetic::remainder(i,arithmetic::abs(j))));
}
@)
void divmod_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (j==0) throw std::runtime_error("DivMod by zero");
  if (l!=expression_base::no_value)
  { push_value(new int_value
     (j>0 ? arithmetic::divide(i,j) : -arithmetic::divide(i,-j)));
    push_value(new int_value(arithmetic::remainder(i,arithmetic::abs(j))));
    if (l==expression_base::single_value)
      wrap_tuple(2);
  }
}
@)
void unary_minus_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new int_value(-i));
}
@)
void power_wrapper(expression_base::level l)
{ static shared_int one(new int_value(1));
@/int n=get<int_value>()->val; shared_int i=get<int_value>();
  if (arithmetic::abs(i->val)!=1 and n<0)
    throw std::runtime_error("Negative power of integer");
  if (l==expression_base::no_value)
    return;
@)
  if (i->val==1)
  {@; push_value(one);
      return;
  }
  if (i->val==-1)
  {@; push_value(n%2==0 ? one : i);
      return;
  }
@)
  push_value(new int_value(arithmetic::power(i->val,n)));
}

@*1 Rationals.
%
As mentioned above the operator `/' applied to integers will not denote
integer division, but rather formation of fractions (rational numbers). Since
the |Rational| constructor requires an unsigned denominator, we must make sure
the integer passed to it is positive. The opposite operation of separating a
rational number into numerator and denominator is also provided; this
operation is essential in order to be able to get from rationals back into the
world of integers. Currently this splitting operation has an intrinsic danger,
since the \Cpp-type |int| used in |int_value| is smaller than the types used
for numerator and denominator of |Rational| values (which are |long long int|
respectively |unsigned long long int|).

@< Local function definitions @>=

void fraction_wrapper(expression_base::level l)
{ int d=get<int_value>()->val; int n=get<int_value>()->val;
  if (d==0) throw std::runtime_error("fraction with zero denominator");
  if (d<0) {@; d=-d; n=-n; } // ensure denominator is positive
  if (l!=expression_base::no_value)
    push_value(new rat_value(Rational(n,d)));
}
@)

void unfraction_wrapper(expression_base::level l)
{ Rational q=get<rat_value>()->val;
  if (l!=expression_base::no_value)
  { push_value(new int_value(q.numerator()));
    push_value(new int_value(q.denominator()));
    if (l==expression_base::single_value)
      wrap_tuple(2);
  }
}

@ We define arithmetic operations for rational numbers, made possible thanks to
operator overloading.

@< Local function definitions @>=

void rat_plus_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val;
  Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new rat_value(i+j));
}
@)
void rat_minus_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val;
  Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new rat_value(i-j));
}
@)
void rat_times_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val;
  Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new rat_value(i*j));
}
@)
void rat_divide_wrapper(expression_base::level l)
{ Rational j=get<rat_value>()->val;
  Rational i=get<rat_value>()->val;
  if (j.numerator()==0)
    throw std::runtime_error("Rational division by zero");
  if (l!=expression_base::no_value)
    push_value(new rat_value(i/j));
}
@)
void rat_unary_minus_wrapper(expression_base::level l)
{@; Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new rat_value(Rational(0)-i)); }
@)
void rat_inverse_wrapper(expression_base::level l)
{@; Rational i=get<rat_value>()->val;
  if (i.numerator()==0)
    throw std::runtime_error("Inverse of zero");
  if (l!=expression_base::no_value)
    push_value(new rat_value(Rational(1)/i)); }
@)
void rat_power_wrapper(expression_base::level l)
{ int n=get<int_value>()->val; Rational b=get<rat_value>()->val;
  if (b.numerator()==0 and n<0)
    throw std::runtime_error("Negative power of zero");
  if (l!=expression_base::no_value)
    push_value(new rat_value(b.power(n)));
}

@*1 Booleans.
%
Relational operators are of the same flavour.
@< Local function definitions @>=

void int_eq_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val==j->val));
}
@)
void int_neq_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val!=j->val));
}
@)
void int_less_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val<j->val));
}
@)
void int_lesseq_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val<=j->val));
}
@)
void int_greater_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val>j->val));
}
@)
void int_greatereq_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val>=j->val));
}

@ We do that again for rational numbers

@< Local function definitions @>=

void rat_eq_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>(); shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val==j->val));
}
@)
void rat_neq_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>(); shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val!=j->val));
}
@)
void rat_less_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>(); shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val<j->val));
}
@)
void rat_lesseq_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>(); shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val<=j->val));
}
@)
void rat_greater_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>(); shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val>j->val));
}
@)
void rat_greatereq_wrapper(expression_base::level l)
{ shared_rat j=get<rat_value>(); shared_rat i=get<rat_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val>=j->val));
}

@ For booleans we also have equality and ineqality.
@< Local function definitions @>=

void equiv_wrapper(expression_base::level l)
{ bool a=get<bool_value>()->val; bool b=get<bool_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new bool_value(a==b));
}
@)
void inequiv_wrapper(expression_base::level l)
{ bool a=get<bool_value>()->val; bool b=get<bool_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new bool_value(a!=b));
}

@*1 Strings.
%
The string type is intended mostly for preparing output to be printed, so few
operations are defined for it. We define functions for comparing and
for concatenating them, and one for converting integers to their string
representation.

@< Local function definitions @>=

void string_eq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val==j->val));
}
@)
void string_leq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val<=j->val));
}
@)
void concatenate_wrapper(expression_base::level l)
{ shared_string b=get<string_value>(); shared_string a=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(new string_value(a->val+b->val));
}
@)
void int_format_wrapper(expression_base::level l)
{ shared_int n=get<int_value>();
  std::ostringstream o; o<<n->val;
  if (l!=expression_base::no_value)
    push_value(new string_value(o.str()));
}

@ To give a rudimentary capability of analysing strings, we provide, in
addition to the subscripting operation, a function to convert the first
character of a string into a numeric value.

@< Local function definitions @>=

void string_to_ascii_wrapper(expression_base::level l)
{ shared_string c=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(new int_value
      (c->val.size()==0 ? -1 : (unsigned char)c->val[0]));
}
@)
void ascii_char_wrapper(expression_base::level l)
{ int c=get<int_value>()->val;
  if (c<' ' or c>'~')
    throw std::runtime_error("Value "+str(c)+" out of range");
  if (l!=expression_base::no_value)
    push_value(new string_value(std::string(1,c)));
}


@*1 Size-of and other generic operators.
%
For the size-of operator we provide several specific bindings: for strings,
vectors and matrices.

@< Local function definitions @>=
void sizeof_string_wrapper(expression_base::level l)
{ size_t s=get<string_value>()->val.size();
  if (l!=expression_base::no_value)
    push_value(new int_value(s));
}
@)
void sizeof_vector_wrapper(expression_base::level l)
{ size_t s=get<vector_value>()->val.size();
  if (l!=expression_base::no_value)
    push_value(new int_value(s));
}

@)
void matrix_bounds_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l==expression_base::no_value)
    return;
  push_value(new int_value(m->val.numRows()));
  push_value(new int_value(m->val.numColumns()));
  if (l==expression_base::single_value)
    wrap_tuple(2);
}

@ Here are functions for extending vectors one or many elements at a time.

@< Local function definitions @>=
void vector_suffix_wrapper(expression_base::level l)
{ int e=get<int_value>()->val;
  shared_vector r=get_own<vector_value>();
  if (l!=expression_base::no_value)
  {@; r->val.push_back(e);
    push_value(r);
  }
}
@)
void vector_prefix_wrapper(expression_base::level l)
{ shared_vector r=get_own<vector_value>();
  int e=get<int_value>()->val;
  if (l!=expression_base::no_value)
  {@; r->val.insert(r->val.begin(),e);
    push_value(r);
  }
}
@)
void join_vectors_wrapper(expression_base::level l)
{ shared_vector y=get<vector_value>();
  shared_vector x=get<vector_value>();
  if (l!=expression_base::no_value)
  { vector_ptr result(new vector_value(std::vector<int>()));
    result->val.reserve(x->val.size()+y->val.size());
    result->val.insert(result->val.end(),x->val.begin(),x->val.end());
    result->val.insert(result->val.end(),y->val.begin(),y->val.end());
    push_value(result);
  }

}

@ Finally, as last function of general utility, one that breaks off
computation with an error message.

@< Local function definitions @>=
void error_wrapper(expression_base::level l)
{@; throw std::runtime_error(get<string_value>()->val); }

@*1 Vectors and matrices.
%
We now define a few functions, to really exercise something, even if it is
modest, from the Atlas library. These wrapper function are not really to be
considered part of the interpreter, but a first step to its interface with the
Atlas library, which is developed in much more detail in the compilation
unit \.{built-in-types}. In fact we shall make some of these wrapper functions
externally callable, so they can be directly used from that compilation unit.

We start with vector and matrix equality comparisons, which are quite similar
to what we saw for rationals, for instance.

@< Local function definitions @>=
void vec_eq_wrapper(expression_base::level l)
{ shared_vector j=get<vector_value>(); shared_vector i=get<vector_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val==j->val));
}
@)
void vec_neq_wrapper(expression_base::level l)
{ shared_vector j=get<vector_value>(); shared_vector i=get<vector_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val!=j->val));
}
@)
void mat_eq_wrapper(expression_base::level l)
{ shared_matrix j=get<matrix_value>(); shared_matrix i=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val==j->val));
}
@)
void mat_neq_wrapper(expression_base::level l)
{ shared_matrix j=get<matrix_value>(); shared_matrix i=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val!=j->val));
}

@ Null vectors and matrices are particularly useful as starting values. In
addition, the latter can produce empty matrices without any (null) entries,
when either the number of rows or column is zero but the other is not; such
matrices (which are hard to obtain by other means) are good starting points
for iterations that consist of adding a number of rows or columns of equal
size, and they determine this size even if none turn out to be contributed.
Since vectors are treated as column vectors, their transpose is a one-line
matrix; such matrices (like null matrices) cannot be obtained from the special
matrix-building expressions.

Since in general built-in functions may throw exceptions, we hold the pointers
to local values in smart pointers; for values popped from the stack this would
in fact be hard to avoid.

@< Local function definitions @>=
void null_vec_wrapper(expression_base::level lev)
{ int l=get<int_value>()->val;
  if (lev!=expression_base::no_value)
    push_value(new vector_value(int_Vector(std::abs(l),0)));
}
@) void null_mat_wrapper(expression_base::level lev)
{ int l=get<int_value>()->val;
  int k=get<int_value>()->val;
  if (lev!=expression_base::no_value)
    push_value(new matrix_value
      (matrix::Matrix<int>(std::abs(k),std::abs(l),0)));
}
void transpose_vec_wrapper(expression_base::level l)
{ shared_vector v=get<vector_value>();
  if (l!=expression_base::no_value)
  { matrix_ptr m (new matrix_value(matrix::Matrix<int>(1,v->val.size())));
    for (size_t j=0; j<v->val.size(); ++j)
      m->val(0,j)=v->val[j];
    push_value(m);
  }
}

@ The wrapper functions for matrix transposition and identity matrix are
called from \.{built-in-types.w}.

@< Declarations of exported functions @>=
void transpose_mat_wrapper (expression_base::level);
void id_mat_wrapper(expression_base::level l);

@ Their definitions are particularly simple, as they just call a matrix method
to do the work.

@< Global function definitions @>=
@) void transpose_mat_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(new matrix_value(m->val.transposed()));
}
@)
void id_mat_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new matrix_value(matrix::Matrix<int>(std::abs(i)))); // identity
}

@ We also define |diagonal_wrapper|, a slight generalisation of
|id_mat_wrapper| that produces a diagonal matrix from a vector. The function
|vector_div_wrapper| produces a rational vector, for which we also provide
addition and subtraction.

@< Local function def... @>=
void diagonal_wrapper(expression_base::level l)
{ shared_vector d=get<vector_value>();
  if (l==expression_base::no_value)
    return;
  size_t n=d->val.size();
  matrix_ptr m (new matrix_value(matrix::Matrix<int>(n)));
  for (size_t i=0; i<n; ++i)
    m->val(i,i)=d->val[i];
  push_value(m);
}
@)
void vector_div_wrapper(expression_base::level l)
{ int n=get<int_value>()->val;
  shared_vector v=get<vector_value>();
  if (l!=expression_base::no_value)
    push_value(new rational_vector_value(v->val,n));
}
@)
void ratvec_unfraction_wrapper(expression_base::level l)
{ shared_rational_vector v = get<rational_vector_value>();
  if (l!=expression_base::no_value)
  { Weight num(v->val.numerator().begin(),v->val.numerator().end()); // convert
    push_value(new vector_value(num));
    push_value(new int_value(v->val.denominator()));
    if (l==expression_base::single_value)
      wrap_tuple(2);
  }
}
@)
void ratvec_plus_wrapper(expression_base::level l)
{ shared_rational_vector v1= get<rational_vector_value>();
  shared_rational_vector v0= get<rational_vector_value>();
  if (l!=expression_base::no_value)
    push_value(new rational_vector_value((v0->val+v1->val).normalize()));
}
@)
void ratvec_minus_wrapper(expression_base::level l)
{ shared_rational_vector v1= get<rational_vector_value>();
  shared_rational_vector v0= get<rational_vector_value>();
  if (l!=expression_base::no_value)
    push_value(new rational_vector_value((v0->val-v1->val).normalize()));
}


@ Now the products between vector and/or matrices. We make the wrapper
|mm_prod_wrapper| around matrix multiplication callable from other compilation
units; for the other wrappers this is not necessary and the will be kept
local.
@< Declarations of exported functions @>=
void mm_prod_wrapper (expression_base::level);

@ For wrapper functions with multiple arguments, we must always remember that
they are to be popped from the stack in reverse order.

@< Global function definitions @>=
void mm_prod_wrapper(expression_base::level l)
{ shared_matrix rf=get<matrix_value>(); // right factor
  shared_matrix lf=get<matrix_value>(); // left factor
  if (lf->val.numColumns()!=rf->val.numRows())
  { std::ostringstream s;
    s<< "Size mismatch " << lf->val.numColumns() << ":" << rf->val.numRows();
    throw std::runtime_error(s.str());
  }
  if (l!=expression_base::no_value)
    push_value(new matrix_value(lf->val*rf->val));
}

@ The other product operations are very similar. As a historic note, the
function corresponding to |mv_prod_wrapper| was in fact our first function
with more than one argument (arithmetic on integer constants was done inside
the parser at that time).

@< Local function definitions @>=
void mv_prod_wrapper(expression_base::level l)
{ shared_vector v=get<vector_value>();
  shared_matrix m=get<matrix_value>();
  if (m->val.numColumns()!=v->val.size())
    throw std::runtime_error(std::string("Size mismatch ")@|
     + str(m->val.numColumns()) + ":" + str(v->val.size()));
  if (l!=expression_base::no_value)
    push_value(new vector_value(m->val*v->val));
}
@)
void mrv_prod_wrapper(expression_base::level l)
{ shared_rational_vector v=get<rational_vector_value>();
  shared_matrix m=get<matrix_value>();
  if (m->val.numColumns()!=v->val.size())
    throw std::runtime_error(std::string("Size mismatch ")@|
     + str(m->val.numColumns()) + ":" + str(v->val.size()));
  if (l!=expression_base::no_value)
    push_value(new rational_vector_value(m->val*v->val));
}
@)
void vv_prod_wrapper(expression_base::level l)
{ shared_vector w=get<vector_value>();
  shared_vector v=get<vector_value>();
  if (v->val.size()!=w->val.size())
    throw std::runtime_error(std::string("Size mismatch ")@|
     + str(v->val.size()) + ":" + str(w->val.size()));
  if (l!=expression_base::no_value)
    push_value(new int_value(v->val.dot(w->val)));
}
@)
void vm_prod_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>(); // right factor
  shared_vector v=get<vector_value>(); // left factor
  if (v->val.size()!=m->val.numRows())
  { std::ostringstream s;
    s<< "Size mismatch " << v->val.size() << ":" << m->val.numRows();
    throw std::runtime_error(s.str());
  }
  if (l!=expression_base::no_value)
    push_value(new vector_value(m->val.right_mult(v->val)));
}

@ Here is the column echelon function.
@h "matreduc.h"
@h "bitmap.h"

@<Local function definitions @>=
void echelon_wrapper(expression_base::level l)
{ shared_matrix M=get_own<matrix_value>();
  if (l!=expression_base::no_value)
  { BitMap pivots=matreduc::column_echelon(M->val);
    push_value(M);
    row_ptr p_list (new row_value(0)); p_list->val.reserve(pivots.size());
    for (BitMap::iterator it=pivots.begin(); it(); ++it)
      p_list->val.push_back(shared_value(new int_value(*it)));
    push_value(p_list);
    if (l==expression_base::single_value)
      wrap_tuple(2);
  }
}

@ And here are general functions |diagonalize| and |adapted_basis|, rather
similar to Smith normal form, but without divisibility guarantee on diagonal
entries. While |diagonalize| provides the matrices applied on the left and
right to obtain diagonal form, |adapted_basis| gives only the left factor (row
operations applied) and gives it inverted, so that this matrix
right-multiplied by the diagonal matrix has the same image as the original
matrix.

@<Local function definitions @>=
void diagonalize_wrapper(expression_base::level l)
{ shared_matrix M=get<matrix_value>();
  if (l!=expression_base::no_value)
  { matrix_ptr row(new matrix_value(matrix::Matrix<int>())),
            column(new matrix_value(matrix::Matrix<int>()));
    vector_ptr diagonal(
       new vector_value(matreduc::diagonalise(M->val,row->val,column->val)));
    push_value(diagonal);
    push_value(row);
    push_value(column);
    if (l==expression_base::single_value)
      wrap_tuple(3);
  }
}
@)
void adapted_basis_wrapper(expression_base::level l)
{ shared_matrix M=get<matrix_value>();
  if (l!=expression_base::no_value)
  { vector_ptr diagonal(new vector_value(std::vector<int>()));
    matrix_ptr basis
      (new matrix_value(matreduc::adapted_basis(M->val,diagonal->val)));
    push_value(basis);
    push_value(diagonal);
    if (l==expression_base::single_value)
      wrap_tuple(2);
  }
}

@ There are two particular applications of the |diagonalise| function defined
in the Atlas library: |kernel| which for any matrix~$M$ will find another
whose image is precisely the kernel of~$M$, and |eigen_lattice| which is a
special case for square matrices~$A$, where the kernel of $A-\lambda\id$ is
computed. We include them here to enable testing these functions.

@h "lattice.h"

@<Local function definitions @>=
void kernel_wrapper(expression_base::level l)
{ shared_matrix M=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(new matrix_value(lattice::kernel(M->val)));
}
@)
void eigen_lattice_wrapper(expression_base::level l)
{ int eigen_value = get<int_value>()->val;
  shared_matrix M=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value (new matrix_value(lattice::eigen_lattice(M->val,eigen_value)));
}
@)
void row_saturate_wrapper(expression_base::level l)
{ shared_matrix M=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(new matrix_value(lattice::row_saturate(M->val)));
}

@ As a last example, here is the Smith normal form algorithm. We provide both
the invariant factors and the rewritten basis on which the normal for is
assumed, as separate functions, and the two combined into a single function.

@< Local function definitions @>=
void invfact_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l==expression_base::no_value)
    return;
  vector_ptr inv_factors (new vector_value(std::vector<int>()));
@/matreduc::Smith_basis(m->val,inv_factors->val);
  push_value(inv_factors);
}
@)
void Smith_basis_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l==expression_base::no_value)
    return;
  vector_ptr inv_factors (new vector_value(std::vector<int>()));
@/push_value(new matrix_value(matreduc::Smith_basis(m->val,inv_factors->val)));
}
@)
void Smith_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l==expression_base::no_value)
    return;
  vector_ptr inv_factors (new vector_value(std::vector<int>()));
@/push_value(new matrix_value(matreduc::Smith_basis(m->val,inv_factors->val)));
  push_value(inv_factors);
  if (l==expression_base::single_value)
    wrap_tuple(2);
}

@ Here is one more wrapper function that uses the Smith normal form algorithm,
but behind the scenes, namely to invert a matrix. Since this cannot be done in
general over the integers, we return an integral matrix and a common
denominator to be applied to all coefficients.
@< Local function definitions @>=
void invert_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (m->val.numRows()!=m->val.numColumns())
  { std::ostringstream s;
    s<< "Cannot invert a " @|
     << m->val.numRows() << "x" << m->val.numColumns() << " matrix";
    throw std::runtime_error(s.str());
  }
  if (l==expression_base::no_value)
    return;
  int_ptr denom(new int_value(0));
@/push_value(new matrix_value(m->val.inverse(denom->val)));
  push_value(denom);
  if (l==expression_base::single_value)
    wrap_tuple(2);
}

@ We must not forget to install what we have defined. The names of the
arithmetic operators correspond to the ones used in the parser definition
file \.{parser.y}.

@< Initialise... @>=
install_function(plus_wrapper,"+","(int,int->int)");
install_function(minus_wrapper,"-","(int,int->int)");
install_function(times_wrapper,"*","(int,int->int)");
install_function(divide_wrapper,"\\","(int,int->int)");
install_function(modulo_wrapper,"%","(int,int->int)");
install_function(divmod_wrapper,"\\%","(int,int->int,int)");
install_function(unary_minus_wrapper,"-","(int->int)");
install_function(power_wrapper,"^","(int,int->int)");
install_function(fraction_wrapper,"/","(int,int->rat)");
install_function(unfraction_wrapper,"%","(rat->int,int)");
   // unary \% means ``break open''
install_function(rat_plus_wrapper,"+","(rat,rat->rat)");
install_function(rat_minus_wrapper,"-","(rat,rat->rat)");
install_function(rat_times_wrapper,"*","(rat,rat->rat)");
install_function(rat_divide_wrapper,"/","(rat,rat->rat)");
install_function(rat_unary_minus_wrapper,"-","(rat->rat)");
install_function(rat_inverse_wrapper,"/","(rat->rat)");
install_function(rat_power_wrapper,"^","(rat,int->rat)");
install_function(int_eq_wrapper,"=","(int,int->bool)");
install_function(int_neq_wrapper,"!=","(int,int->bool)");
install_function(int_less_wrapper,"<","(int,int->bool)");
install_function(int_lesseq_wrapper,"<=","(int,int->bool)");
install_function(int_greater_wrapper,">","(int,int->bool)");
install_function(int_greatereq_wrapper,">=","(int,int->bool)");
install_function(rat_eq_wrapper,"=","(rat,rat->bool)");
install_function(rat_neq_wrapper,"!=","(rat,rat->bool)");
install_function(rat_less_wrapper,"<","(rat,rat->bool)");
install_function(rat_lesseq_wrapper,"<=","(rat,rat->bool)");
install_function(rat_greater_wrapper,">","(rat,rat->bool)");
install_function(rat_greatereq_wrapper,">=","(rat,rat->bool)");
install_function(equiv_wrapper,"=","(bool,bool->bool)");
install_function(inequiv_wrapper,"!=","(bool,bool->bool)");
install_function(string_eq_wrapper,"=","(string,string->bool)");
install_function(string_leq_wrapper,"<=","(string,string->bool)");
install_function(concatenate_wrapper,"#","(string,string->string)");
install_function(int_format_wrapper,"int_format","(int->string)");
install_function(string_to_ascii_wrapper,"ascii","(string->int)");
install_function(ascii_char_wrapper,"ascii","(int->string)");
install_function(sizeof_string_wrapper,"#","(string->int)");
install_function(sizeof_vector_wrapper,"#","(vec->int)");
install_function(matrix_bounds_wrapper,"#","(mat->int,int)");
install_function(vector_div_wrapper,"/","(vec,int->ratvec)");
install_function(ratvec_unfraction_wrapper,"%","(ratvec->vec,int)");
install_function(ratvec_plus_wrapper,"+","(ratvec,ratvec->ratvec)");
install_function(ratvec_minus_wrapper,"-","(ratvec,ratvec->ratvec)");
install_function(null_vec_wrapper,"null","(int->vec)");
install_function(null_mat_wrapper,"null","(int,int->mat)");
install_function(id_mat_wrapper,"id_mat","(int->mat)");
install_function(error_wrapper,"error","(string->*)");
install_function(vector_suffix_wrapper,"#","(vec,int->vec)");
install_function(vector_prefix_wrapper,"#","(int,vec->vec)");
install_function(join_vectors_wrapper,"#","(vec,vec->vec)");
install_function(vec_eq_wrapper,"=","(vec,vec->bool)");
install_function(vec_neq_wrapper,"!=","(vec,vec->bool)");
install_function(mat_eq_wrapper,"=","(mat,mat->bool)");
install_function(mat_neq_wrapper,"!=","(mat,mat->bool)");
install_function(transpose_mat_wrapper,"^","(mat->mat)");
install_function(transpose_vec_wrapper,"^","(vec->mat)");
install_function(diagonal_wrapper,"diagonal","(vec->mat)");
install_function(vv_prod_wrapper,"*","(vec,vec->int)");
install_function(mrv_prod_wrapper,"*","(mat,ratvec->ratvec)");
install_function(mv_prod_wrapper,"*","(mat,vec->vec)");
install_function(mm_prod_wrapper,"*","(mat,mat->mat)");
install_function(vm_prod_wrapper,"*","(vec,mat->vec)");
install_function(echelon_wrapper,"echelon","(mat->mat,[int])");
install_function(diagonalize_wrapper,"diagonalize","(mat->vec,mat,mat)");
install_function(adapted_basis_wrapper,"adapted_basis","(mat->mat,vec)");
install_function(kernel_wrapper,"kernel","(mat->mat)");
install_function(eigen_lattice_wrapper,"eigen_lattice","(mat,int->mat)");
install_function(row_saturate_wrapper,"row_saturate","(mat->mat)");
install_function(invfact_wrapper,"inv_fact","(mat->vec)");
install_function(Smith_basis_wrapper,"Smith_basis","(mat->mat)");
install_function(Smith_wrapper,"Smith","(mat->mat,vec)");
install_function(invert_wrapper,"invert","(mat->mat,int)");

@* Index.

% Local IspellDict: british
