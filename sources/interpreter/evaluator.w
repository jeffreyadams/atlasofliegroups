% Copyright (C) 2006-2009 Marc van Leeuwen
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
\def\foreign#1{{\sl#1\/}}
\def\id{\mathop{\rm id}}

@* Outline.
%
This file describes the central part of the interpreter for the (new) command
language of the Atlas of Lie Groups and Representation software. This part is
concerned with the analysis and execution of expressions that have already
been processed by the parser. This requires substantial work; we may decide to
relegate the implementation of some basic user types like vectors and
matrices, which does not interact much with the general evaluation procedure,
to another compilation unit, but this separation has not yet been made.
@( evaluator.h @>=

#ifndef EVALUATOR_H
#define EVALUATOR_H

#include "types.h"

@< Includes needed in the header file @>@;
namespace atlas { namespace interpreter {
@< Type definitions @>@;
@< Declarations of global variables @>@;
@< Declarations of exported functions @>@;
@< Template and inline function definitions @>@;
}@; }@;
#endif

@ The implementation unit follows a somewhat similar pattern.

@h "evaluator.h"
@h <cstdlib>
@c
namespace atlas { namespace interpreter {
@< Global variable definitions @>@;
namespace {
@< Local variable definitions @>@;
@< Local function definitions @>@;
}@; // |namespace|
@< Function definitions @>@;
}@; }@;

@* Preliminaries.
%
We start with a small template function to help giving sensible error messages.

@< Includes needed in the header file @>=
#include <string>
#include <sstream>

@~The template function |str| just returns an integer (or any printable value)
represented as a string (in decimal). The reason for using a template is that
without them it is hard to do a decent job for both signed and unsigned types.
Using string streams, the definition of~|str| is trivial.

@< Template and inline function definitions @>=
template <typename T>
  std::string str(T n) @+{@; std::ostringstream s; s<<n; return s.str(); }

@ Before executing anything the evaluator needs some initialisation, called
from the main program.

@< Declarations of exported functions @>=
void initialise_evaluator();

@ The details of this initialisation will be given when the variables involved
are introduced.

@< Function definitions @>=
void initialise_evaluator()
@+{@; @< Initialise evaluator @> }

@~Although not necessary, the following will avoid some early reallocations of
|execution_stack|.

@< Initialise evaluator @>=
execution_stack.reserve(16); // avoid some early reallocations

@ Here is one more function related to the execution stack. The stack owns the
values it contains, but there was no reason to wrap it into a class with a
destructor, since we never intend to destroy the stack entirely: if our
program exits either peacefully or by an uncaught exception we don't care
about some values that are not destroyed. We must remember however to |delete|
the values whenever we empty the stack after catching a runtime error. We
provide a function to reset the evaluator after catching an exception that
does this, as well as other actions needed to have the evaluator start with a
clear slate.

@< Declarations of exported functions @>=
void reset_evaluator ();

@~We shall monitor disappearing values when clearing the stack. This provides
some output that may not be easy to interpret for an unsuspecting user, so we
limit this to the case that the |verbosity| parameter (defined below) is
nonzero. Other actions necessary for a clean restart will be added later.

@< Function definitions @>=
void reset_evaluator()
{ if (!execution_stack.empty())
  { if (verbosity==0)
      execution_stack.clear();
    else
    {
      std::cerr << "Discarding from execution stack:" << std::endl;
      do
        std::cerr << *pop_value() << std::endl;
      while (!execution_stack.empty());
    }
  }
  @< Actions to reset the evaluator @>
}

@ Here are some more predefined types, less fundamental than those declared in
\.{types.h} (since not actually used in any of the core language constructs),
but useful for instance for specifying various coercions.

@< Declarations of global variables @>=
extern const type_expr rat_type; // \.{rat}
extern const type_expr str_type; // \.{string}
extern const type_expr vec_type; // \.{vec}
extern const type_expr ratvec_type; // \.{ratvec}
extern const type_expr mat_type; // \.{mat}
extern const type_expr row_of_int_type; // \.{[int]}
extern const type_expr row_of_rat_type; // \.{[rat]}
extern const type_expr row_of_vec_type; // \.{[vec]}
extern const type_expr row_row_of_int_type; // \.{[[int]]}
extern const type_expr pair_type; // \.{(*,*)}
extern const type_expr int_int_type; // \.{(int,int)}

@ The construction of type constants follows the same pattern as before,
calling |copy| in the case of composite types.

@< Global variable definitions @>=
const type_expr rat_type(rational_type);
const type_expr str_type(string_type);
const type_expr vec_type(vector_type);
const type_expr ratvec_type(rational_vector_type);
const type_expr mat_type(matrix_type);
const type_expr row_of_int_type(copy(int_type));
const type_expr row_of_rat_type(copy(rat_type));
const type_expr row_of_vec_type(copy(vec_type));
const type_expr row_row_of_int_type(copy(row_of_int_type));
const type_expr pair_type(*unknown_tuple(2));  // copy and destroy original
const type_expr int_int_type(*make_type("(int,int)")); // idem


@ Now we derive the first ``primitive'' value types. The type for rational
numbers is in fact implemented in the atlas library, so we must include a
header file into ours.

@<Includes needed in the header file @>=
#include "arithmetic.h"

@ For each type we define a corresponding auto-pointer type (whose name ends
with \&{\_ptr}), since we shall often need to hold such values by pointers,
and the risk of exceptions is ever present.

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
{ arithmetic::Rational val;
@)
  explicit rat_value(arithmetic::Rational v) : val(v) @+ {}
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

@ Now let us define our first class derived from |expression_base|, which is
|denotation|; it simply stores a |value|, which it returns upon evaluation.
The constructor is passed a smart pointer for the usual reason: if the
constructor should be part of a |new| expression (as is practice it will) then
the allocation for that |new| happens between the moment of passing the
constructor arguments and the invocation of the constructor, and if that
allocation throws it would produce a memory leak if raw a pointer were passed.
This does force us to cast the arguments to |shared_value| below.

@< Type definitions @>=
struct denotation : public expression_base
{ shared_value denoted_value;
@)
  explicit denotation(shared_value v) : denoted_value(v) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; denoted_value->print(out); }
};

@ The following function is not defined inside the class definition, because
that definition precedes the one of the |inline| function |push_value| in the
header file. It is not even made |inline| itself, since there is little point
in doing so for virtual methods (calls via the vtable cannot be inlined). For
the first time we see the |level| argument in action; whenever |l==no_value|
we should avoid producing a result value.

@< Function def... @>=
void denotation::evaluate(level l) const
@+{@; if (l!=no_value) push_value(denoted_value); }


@* Type-checking and conversion to executable form  of expressions.
%
The expression returned by the parser is type-checked before execution starts.
During type checking, there are subexpressions for which a definite type is
required, and others where none is (like the complete expression). The
difference is important, as in the former case conversions can be inserted to
make types match, for instance between a list of integers and a vector value;
this is in fact the only way the user can produce the latter kind of values.
However, both cases are handled by a single function |convert_expr|, which in
addition builds (upon success) an |expression| value, i.e., one in terms of
which evaluation is defined. Apart from the |expr e@;| produced by the parser,
|convert_expr| takes a type as argument, in the form of a non-constant
reference |type_expr& type@;|. If |type| is undefined initially, then it
will be set to the type derived for the expression; if it is defined then it
will just be tested. It could also be partially defined initially (such as
`\.{(int,*)}', meaning ``pair on an integer and something''), in which case
derivation and testing functionality are combined; this gives flexibility to
|convert_expr|. Upon successful completion, |type| will usually have become a
completely defined. The object |type| should be owned by the caller, who will
automatically gain ownership of any new nodes added, which will be so due to
calling the |specialise| method for |type| or for its descendants. In some
cases |type| will remain partly undefined, like for an emtpy list display
which gets type~`\.{[*]}'; however if |type| remains completely undefined
`\.*' (as would happen for the selection of a value from an empty list, or for
a function that is unconditionally recursive) then it can be seen that
evaluation cannot possibly complete without error, so we might treat this case
as a type error (we do so occasionally if it simplifies our code).

With the translation of expressions into |expression| values, there is no need
for a separate evaluator function: we shall just call the |evaluate| method of
the expression object, which will do its work with the help of the similar
methods of its subexpressions, resulting in an indirectly recursive evaluation
framework.

@< Declarations of exported functions @>=
expression convert_expr(const expr& e, type_expr& type)
  throw(std::bad_alloc,program_error);

@ In the function |convert_expr| we shall need a type for storing bindings
between identifiers and types. We use a vector, and since these cannot hold
auto-pointers, we need to define a destructor to clean up the types. Different
such binding vectors will be stacked, for nested scopes, but using an STL
container for that would necessitate defining a copy constructor, which would
be a painful operation: (1)~it must call |copy| on the types held in the
bindings, in order to avoid double destruction, (2)~these calls to |copy|
could throw an exception, and (3)~the constructor won't be complete, and the
destructor therefore not activated, until the final entry is copied, so we
would need a |try|\dots|catch| in the constructor to avoid a memory leak.

So instead we chain the different bindings into a linked list, and forbid any
copy or assignment. In fact the nesting of the various bindings will embed in
the nested calls of |convert_expr|, so we can afford to allocate each
|bindings| as local variable to (some instace of) |convert_expr|, and there is
no ownership of them down the linked list. We do provide methods |push| and
|pop| to prepend and detach a |bindings| from a list.

@< Type def... @>=
class bindings
: public std::vector<std::pair<Hash_table::id_type,type_p> >
{ typedef std::vector<std::pair<Hash_table::id_type,type_p> >
    base;
  bindings* next; // non-owned pointer
@)bindings(const bindings&); // copy constructor
  void operator= (const bindings&); // assignment operator
public:
  bindings(size_t n=0) : @[ base() @], next(NULL) @+{@; base::reserve(n); }
    // predict size |n|, informative
  ~bindings () @+
  {@; for (base::iterator it=begin(); it!=end(); ++it)
      delete it->second;
  }
  void add(Hash_table::id_type id,type_ptr t);
  type_p lookup
    (Hash_table::id_type id, size_t& depth, size_t& offset) const;
  void push (bindings*& sp) @+{@; next=sp; sp=this; }
  void pop (bindings*& sp) const @+{@; sp=next; }
};

@ The method |add| adds a pair to the vector of bindings, taking care not to
release the type pointer until the pair is successfully allocated. This is
also a good place to check for the presence of identical identifiers.

@< Function def... @>=
void bindings::add(Hash_table::id_type id,type_ptr t)
{ for (base::const_iterator it=begin(); it!=end(); ++it)
    if (it->first==id)
      throw program_error(std::string("Multiple binding of ")
                          +main_hash_table->name_of(id)
                          +" in same scope");
  push_back(std::make_pair(id,(type_p)NULL));
  back().second=t.release();
}

@ The method |lookup| runs through the linked list of bindings and returns a
pointer to the type if a match was found, also assigning the coordinates to
output arguments. If no match is found a null pointer is returned and the
output parameters are unchanged.

@< Function def... @>=
type_p bindings::lookup
  (Hash_table::id_type id, size_t& depth, size_t& offset) const
{ size_t i=0;
  for (const bindings* p=this; p!=NULL; p=p->next,++i)
    for (size_t j=0; j<p->size(); ++j)
      if ((*p)[j].first==id)
      {@; depth=i; offset=j; return (*p)[j].second; }
  return NULL;
}

@ During conversion of expressions, we keep a stack |id_context| of identifier
bindings in order to detemine their (lexical) bindind and type. We could have
passed this context as an argument to |convert_expr|, but it changes only on
very few occasions, so it is more convenient to use a static variable to hold
the context. This variable could have been |static| local to the function
|convert_expr|, if it were not for the fact that conversion could be aborted by
throwing an exception during type analysis, in which case all active
invocations of |convert_expr| will be terminated, and unless we were willing
to add a |catch| clause at every change to |id_context|, we would have no
occasion to clear the context for the next type analysis. So we decide to
declare this variable local to this compilation unit, and provide a separate
action for clearing it.

@< Local var... @>=
bindings* id_context;

@~The cleaning up action just clears the pointer; in fact the pointed to
objects have already disappeared from the \Cpp\ runtime stack at the time this
code gets executed.

@< Actions to reset the evaluator @>=
id_context=NULL;

@ The function |convert_expr| returns a pointer to the conversion of the
|expr| to |expression|, of which the caller should take ownership. The reason
we don't return an |expression_ptr| is entirely pragmatic. The code below
takes into account the possibility that a denotation is converted immediately
to some other type, for instance integer denotations can be used where a
rational number is expected. The function |coerce| tests for this possibility,
and may modify its final argument correspondingly.

Altogether this is a quite extensive function, with as many cases in the
switch as there are variants of |expr_union|, and for many of those branches a
considerable amount of work to be done. It is therefore convenient to postpone
most of these cases and treat them one syntactic construction at the time.

@< Function definitions @>=
expression convert_expr(const expr& e, type_expr& type)
  throw(std::bad_alloc,program_error)
{

  switch(e.kind)
  { case integer_denotation:
    { expression_ptr d@|(new denotation
        (shared_value(new int_value(e.e.int_denotation_variant))));
      return conform_types(int_type,type,d,e);
    }
   case string_denotation:
    { expression_ptr d@|(new denotation
        (shared_value(new string_value(e.e.str_denotation_variant))));
      return conform_types(str_type,type,d,e);
    }
   case boolean_denotation:
    { expression_ptr d@|(new denotation
          (shared_value(new bool_value(e.e.int_denotation_variant))));
      return conform_types(bool_type,type,d,e);
    }
   @\@< Other cases for type-checking and converting expression~|e| against
   |type|, all of which either |return| or |throw| a |type_error| @>
 }
 return NULL; // keep compiler happy
}

@* List displays.
%
Now we shall consider the handling of list displays, expressions that build a
list of values by constructing each component by an explicit expression. We
start by defining the structure used to represent list displays after
conversion in the type check. The logical name |list_display| for this
structure is already taken, so we call it a |list_expression| instead.

@< Type definitions @>=
struct list_expression : public expression_base
{ std::vector<expression> component;
@)
  explicit list_expression(size_t n) : component(n,NULL) @+{}
   // always start out with null pointers
  virtual ~list_expression();
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ A |list_expression| owns its component expressions, so the destructor must
take care to delete them.

@< Function def... @>=
list_expression::~list_expression()
{@; for (size_t i=0; i<component.size(); ++i) delete component[i];
}

@ When we print a list display, we just print the component expressions,
enclosed in brackets and separated by commas, just as in their input syntax.

@< Function def... @>=
void list_expression::print(std::ostream& out) const
{ out << '[';
  if (component.size()==0) out << ']';
  else
    for (size_t i=0; i<component.size(); ++i)
    out << *component[i] << (i<component.size()-1 ? ',' : ']');
}


@ If a list display has multiple components, they must all have the same
type~$t$, and the result type will be ``row of''~$t$. If a type of that form
is required, the components will be required to have type $t$, if no
particular type is require then the components will just be required to have
equal types; in addition we want to allow additional cases which can be made
to conform by inserting conversion routines.

The function |convert_expr| handles all this in an elegant way. As we proceed
along the components, our type may get specialised, and in case of definitely
unequal types we may look up if a conversion between them is defined. It may
even happen that the type remains partly undetermined, in case of an empty
list. The current code does have a slight asymmetry, in that types from
previous components may guide conversions applied to later components but not
vice versa: if some components are vectors and other are lists of integers,
then everything will be converted to the first kind that occurs (unless the
context requires a particular type, in which case all components are converted
accordingly).

We also have the obligation here to build a converted |list_expression|
representing the list display; we know its size beforehand, and upon
allocation we fill it with null pointers that will progressively be replaced,
in order to ensure proper destruction in case a |type_error| is thrown by a
component expression. When a type other than ``row of'' is expected, we must
of course take explicit action to see whether some type conversion can resolve
the conflict.

@< Other cases for type-checking and converting... @>=
case list_display:
  if (type.specialise(row_of_type))
  { std::auto_ptr<list_expression> result (new list_expression(0));
    result->component.reserve(length(e.e.sublist));
    for (expr_list l=e.e.sublist; l!=NULL; l=l->next)
      result->component.push_back(convert_expr(l->e,*type.component_type));
    return result.release(); // and convert (derived|->|base) to |expression|
  }
  else
  @< If |type| can be converted from some row-of type, check the components of
     |e.e.sublist| against the required type, and apply a conversion function
     to the converted expression; otherwise |throw| a |type_error| @>

@ When in |convert_expr| we encounter a list display when a non-row is
expected, we single out the cases that a conversion from a row type to the
required type is available; in that case we continue to convert the component
expressions with as expected type the corresponding component type (if
multiple coercions to the required type are known, the first one in the table
gets preference; this occurs for required type \.{mat}, and means that the
component type will then be \.{vec} rather than \.{[int]}).

@< If |type| can be converted from some row-of type, check the components of
   |e.e.sublist|... @>=
{ type_expr comp_type;
  conversion_record* conv = row_coercion(type,comp_type);
  if (conv==NULL)
    throw type_error(e,copy(row_of_type),copy(type));
@)
  std::auto_ptr<list_expression> display@|
      (new list_expression@|(length(e.e.sublist)));
  @/size_t i=0;
    for (expr_list l=e.e.sublist; l!=NULL; l=l->next,++i)
      display->component[i]=convert_expr(l->e,comp_type);
    return new conversion(*conv,expression_ptr(display));
}

@ The evaluation of a |list_expression| evaluates the components in a simple
loop. If |l==no_value|, we only evaluate for side effects, otherwise we wrap
the result in a single |row_value|. Since evaluation of component expressions
pushes the resulting value onto the execution stack, we pop each one off
immediately t integrate it into the result. We take care to hold the partial
result via an auto-pointer |result|, so that in case of a runtime error during
the evaluation of one of the component expressions the values already computed
are cleaned up.

@< Function def... @>=
void list_expression::evaluate(level l) const
{ if (l==no_value)
    for (size_t i=0; i<component.size(); ++i)
      component[i]->void_eval();
  else
  { row_ptr result(new row_value(component.size()));
    for (size_t i=0; i<component.size(); ++i)
      component[i]->eval(),result->val[i]=pop_value();
    push_value(result); // result will be shared from here on
  }
}

@* Tuples.
%
Besides types for uniform sequences we shall need types for non-uniform
sequences of values, usually short and with a given sequence of types for
their components; their most obvious and simple use is for argument lists of
functions. Without using these tuple types one might define function types
taking multiple arguments, but then a function could still only yield a single
value. With tuple types multiple results can be produced, and there will be no
need to explicitly cater for functions with multiple arguments.

@*1 Tuple displays.
%
Tuple values can be produced by tuple displays. Syntactically these are very
much like list displays, listing explicitly their component expressions. After
type-checking, they are represented by a |tuple_expression| object (the
name |tuple_display| was already taken); we derive this from
|list_expression|.


@< Type definitions @>=
struct tuple_expression : public list_expression
{ explicit tuple_expression(size_t n) : list_expression(n)@+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ When we print a tuple display, we just print the component expressions,
enclosed in parentheses and separated by commas, to match their input syntax.


@< Function def... @>=
void tuple_expression::print(std::ostream& out) const
{ out << '(';
  if (component.size()==0) out << ')';
  else
    for (size_t i=0; i<component.size(); ++i)
    out << *component[i] << (i<component.size()-1 ? ',' : ')');
}


@ When converting a tuple expression, we first try to specialise |type| to a
tuple type with the right number of unknown components; unless |type| was
completely undetermined, this just amounts to a test that it is a tuple type
with the right number of components. We report a wrong number of components
via a type pattern, which is probably as clear as mentioning too few or too
many components.

@< Other cases for type-checking and converting... @>=
case tuple_display:
{ type_ptr tup=unknown_tuple(length(e.e.sublist));
  if (type.specialise(*tup))
  { std::auto_ptr<tuple_expression> result(new tuple_expression(0));
    result->component.reserve(length(e.e.sublist));
  @/type_list tl=type.tuple;
    for (expr_list el=e.e.sublist; el!=NULL; el=el->next,tl=tl->next)
      result->component.push_back(convert_expr(el->e,tl->t));
    return result.release();  // and convert (derived|->|base) to |expression|
  }
  else throw type_error(e,tup,copy(type));
}

@*1 Evaluating tuple displays.
%
Evaluating a tuple display is much like evaluating a list display: we evaluate
the components in a simple loop. If |l==no_value| this is done for side
effects only, otherwise each component produces (via the |eval| method) a
single value on the stack. Afterwards the result needs to be grouped into a
single value only if |l==single_value|, which is accomplished by |wrap_tuple|.

@< Function def... @>=
void tuple_expression::evaluate(level l) const
{ if (l==no_value)
    for (size_t i=0; i<component.size(); ++i)
      component[i]->void_eval();
  else
  { for (size_t i=0; i<component.size(); ++i)
      component[i]->eval();
    if (l==single_value)
      wrap_tuple(component.size());
  }
}

@* Array subscription.
%
While we have seen expressions to build list, and vectors and matrices out of
them, we so far are not able to access their components once they are
constructed. To that end we shall now introduce operations to index such
values. We allow subscription of rows, but also of vectors, rational vectors
and matrices. Since after type analysis we know which of the cases applies, we
define several classes. These differ mostly by their |evaluate| method, so we
first derive an intermediate class from |expression_base|, and derive the
others from it. This class also serves to host an enumeration type that will
serve later.

@< Type definitions @>=
struct subscr_base : public expression_base
{ enum sub_type @+
  { row_entry, vector_entry, ratvec_entry, matrix_entry, matrix_column };
  expression array, index;
@)
  subscr_base(expression_ptr a, expression_ptr i)
  : array(a.release()),index(i.release()) @+{}
  ~subscr_base() @+ {@; delete array; delete index; }
  void print(std::ostream& out) const;
  static bool indexable
  (const type_expr& aggr,
   const type_expr& index,
         type_expr& subscr,
         sub_type& kind);
};

@)
struct row_subscription : public subscr_base
{ row_subscription(expression_ptr a, expression_ptr i)
  : subscr_base(a,i) @+{}
  virtual void evaluate(level l) const;
};

@)
struct vector_subscription : public subscr_base
{ vector_subscription(expression_ptr a, expression_ptr i)
  : subscr_base(a,i) @+{}
  virtual void evaluate(level l) const;
};
@)
struct ratvec_subscription : public subscr_base
{ ratvec_subscription(expression_ptr a, expression_ptr i)
  : subscr_base(a,i) @+{}
  virtual void evaluate(level l) const;
};
@)
struct matrix_subscription : public subscr_base
{ matrix_subscription(expression_ptr a, expression_ptr ij)
  : subscr_base(a,ij) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};
@)
struct matrix_slice : public subscr_base
{ matrix_slice(expression_ptr a, expression_ptr j)
  : subscr_base(a,j) @+{}
  virtual void evaluate(level l) const;
};

@ These subscriptions are printed in the usual subscription syntax. For matrix
subscriptions, where the index type is \.{(int,int)}, the index expression is
quite likely to be a tuple display, in which case we suppress parentheses.
Since we have passed the type check here, we know that any tuple display is
necessarily a pair.

@< Function definitions @>=
void subscr_base::print(std::ostream& out) const
{@; out << *array << '[' << *index << ']';
}
@)
void matrix_subscription::print(std::ostream& out) const
{ tuple_expression* p=dynamic_cast<tuple_expression*>(index);
  if (p==NULL) out << *array << '[' << *index << ']';
  else
    out << *array << '[' << *p->component[0] << ',' << *p->component[1] << ']';
}

@ It shall be useful to have a function recognising valid aggregate-index
combinations.

@< Function def... @>=
bool subscr_base::indexable
  (const type_expr& aggr,
   const type_expr& index,
         type_expr& subscr,
         sub_type& kind)
{ if (aggr.kind==row_type and index==int_type)
  @/{@; kind=row_entry;
        return subscr.specialise(*aggr.component_type);
  }
  if (aggr==vec_type and index==int_type)
  @/{@; kind=vector_entry;
        return subscr.specialise(int_type);
  }
  if (aggr==ratvec_type and index==int_type)
  @/{@; kind=ratvec_entry;
        return subscr.specialise(rat_type);
  }
  if (aggr!=mat_type)
    return false;
  if (index==int_int_type)
  @/{@; kind=matrix_entry;
        return subscr.specialise(int_type);
  }
  kind=matrix_column;
  return index==int_type and subscr.specialise(vec_type);
}

@ When encountering a subscription in |convert_expr|, we determine the types
of array and of the indexing expression separately, ignoring so far any type
required by the context. Then we look if the types agree with any of the four
types of subscription expressions that we can convert to, throwing an error if
it does not. Finally we check is the a priori type |subscr_type| of the
subscripted expression equals or specialises to the required |type|, or can be
converted to it by |coerce|, again throwing an error if nothing works. For the
indexing expression only equality of types is admitted, since nothing can be
coerced to \.{int} or to \.{(int,int)} anyway, and there is little point in
catering for (indexing) expressions having completely undetermined type (which
can only happen for expressions that cannot be evaluated without error).

@< Other cases for type-checking and converting... @>=
case subscription:
{ type_expr array_type, index_type, subscr_type;
    // initialised to |undetermined_type|
  expression_ptr array
    (convert_expr(e.e.subscription_variant->array,array_type));
  expression_ptr index
    (convert_expr(e.e.subscription_variant->index,index_type));
  subscr_base::sub_type kind;
  expression_ptr subscr;
  if (subscr_base::indexable(array_type,index_type,subscr_type,kind))
    switch (kind)
    { case subscr_base::row_entry:
      subscr.reset(new row_subscription(array,index));
    break;
    case subscr_base::vector_entry:
      subscr.reset(new vector_subscription(array,index));
    break;
    case subscr_base::ratvec_entry:
      subscr.reset(new ratvec_subscription(array,index));
    break;
    case subscr_base::matrix_entry:
      subscr.reset(new matrix_subscription(array,index));
    break;
    case subscr_base::matrix_column:
      subscr.reset(new matrix_slice(array,index));
    break;
    }
  else
  { std::ostringstream o;
    o << "Cannot subscript " << array_type << " value with index of type "
      << index_type;
    throw expr_error(e,o.str());
  }
@)
  return conform_types(subscr_type,type,subscr,e);
}


@ Here are the |evaluate| methods for the various subscription expressions.
For |matrix_subscription|, note that |push_tuple_components()| takes care of
giving access to the individual index values without ownership conflict (by
the time |j| and |i| are accessed, the tuple is already destroyed, but its
components survive).

@< Function definitions @>=
inline std::string range_mess(int i,size_t n,const expression_base* e)
{ std::ostringstream o;
  e->print(o << "index " << i << " out of range (<" << n
             << ") in subscription ");
  return o.str();
}
@)
void row_subscription::evaluate(level l) const
{ shared_int i=((index->eval(),get<int_value>()));
  shared_row r=((array->eval(),get<row_value>()));
  if (static_cast<unsigned int>(i->val)>=r->val.size())
    throw std::runtime_error(range_mess(i->val,r->val.size(),this));
  push_expanded(l,r->val[i->val]);
}
@)
void vector_subscription::evaluate(level l) const
{ shared_int i=((index->eval(),get<int_value>()));
  shared_vector v=((array->eval(),get<vector_value>()));
  if (static_cast<unsigned int>(i->val)>=v->val.size())
    throw std::runtime_error(range_mess(i->val,v->val.size(),this));
  if (l!=no_value)
    push_value(new int_value(v->val[i->val]));
}
@)
void ratvec_subscription::evaluate(level l) const
{ shared_int i=((index->eval(),get<int_value>()));
  shared_rational_vector v=((array->eval(),get<rational_vector_value>()));
  if (static_cast<unsigned int>(i->val)>=v->val.size())
    throw std::runtime_error(range_mess(i->val,v->val.size(),this));
  if (l!=no_value)
    push_value(new rat_value(arithmetic::Rational @|
       (v->val.numerator()[i->val],v->val.denominator())));
}
@)
void matrix_subscription::evaluate(level l) const
{ index->multi_eval(); @+
  shared_int j=get<int_value>();
  shared_int i=get<int_value>();
  shared_matrix m=((array->eval(),get<matrix_value>()));
  if (static_cast<unsigned int>(i->val)>=m->val.numRows())
    throw std::runtime_error
     ("initial "+range_mess(i->val,m->val.numRows(),this));
  if (static_cast<unsigned int>(j->val)>=m->val.numColumns())
    throw std::runtime_error
     ("final "+range_mess(j->val,m->val.numColumns(),this));
  if (l!=no_value)
    push_value(new int_value(m->val(i->val,j->val)));
}
@)
void matrix_slice::evaluate(level l) const
{ shared_int j=((index->eval(),get<int_value>()));
  shared_matrix m=((array->eval(),get<matrix_value>()));
  if (static_cast<unsigned int>(j->val)>=m->val.numColumns())
    throw std::runtime_error(range_mess(j->val,m->val.numColumns(),this));
  if (l!=no_value)
    push_value(new vector_value(m->val.column(j->val)));
}


@* Identifiers.
%
Calling functions, whether built-in or user-defined, requires using
identifiers, but they occur in other contexts as well, and their treatment is
no different between these cases. They come in two kinds, local and global
identifiers, which are treated in fairly different way, and during type
analysis the two are converted into different kinds of |expression|. The most
fundamental difference is that for global identifiers a value is already known
at the time the identifier expression is type-checked, and the type of this
value can be used; for local identifiers only a type is associated to the
identifier at that time, and indeed during different evaluations the same
local identifier may find itself bound to different values.

Global identifiers values will be stored in a global identifier table holding
values and their types (initially it does so for the built-in functions). The
values of local identifiers will be stored at runtime in a stack of variable
bindings, but not their types (these are held elsewhere during type analysis,
and have disappeared at evaluation time).

In spite of these differences there is some common ground for global and local
identifiers: they have a name that can be printed. For this reason we derive
an intermediate structure from |expression_base| that will serve as base for
both kinds of applied identifier expressions.

@< Type definitions @>=
struct identifier : public expression_base
{ Hash_table::id_type code;
@)
  explicit identifier(Hash_table::id_type id) : code(id) @+{}
  virtual ~identifier() @+ {}
  const char* name() const;
  virtual void print(std::ostream& out) const;
};

@ To print an identifier, we get its name from the main hash table.

@< Function definitions @>=
const char* identifier::name() const
@+{@; return main_hash_table->name_of(code); }
void identifier::print(std::ostream& out) const
@+{@; out<< name(); }

@*1 The global identifier table.
%
We need an identifier table to record the values of globally bound identifiers
(such as those for built-in functions) and their types. The values are held in
shared pointers, so that we can evaluate a global identifier without
duplicating the value in the table itself. Modifying the value of such an
identifier by an assignment will produce a new pointer, so that
``shareholders'' of the old value will not see any change. There is another
level of sharing, which affects applied occurrences of the identifier as
converted during type analysis. The value accessed by such identifiers (which
could be contained in user-defined function bodies and therefore have long
lifetime) are expected to undergo change when a new value is assigned to the
global variable; they will therefore access the location of the shared value
pointer rather than the value pointed to. However, if a new identifier of the
same name should be introduced, a new value pointer stored in a different
location will be created, while existing applied occurrences of the identifier
will continue to access the old value, avoiding the possibility of accessing a
value of unexpected type. In such a circumstance, the old shared pointer
location itself will no longer be owned by the identifier table, so we should
arrange for shared ownership of that location. This explains that the
|id_data| structure used for entries in the table has a shared pointer to a
shared pointer.

For the type component on the other hand, the identifier table will assume
strict ownership. Giving ownership of the type directly to |id_data| would
complicate its duplication, and therefore its insertion into the table. We
then only allow construction of |id_data| objects in an empty state; the
pointers should be set only after insertion of the |id_data| into the table,
which then assumes ownership of the type.

@< Type definitions @>=

typedef std::tr1::shared_ptr<shared_value> shared_share;
struct id_data
{ shared_share val; @+ type_p type;
  id_data() : val(),type(NULL)@+ {}
};

@ We cannot store auto pointers in a table, so upon entering into the table we
convert type pointers to ordinary pointers, and this is what |type_of| lookup
will return (the table retains ownership); destruction of the type expressions
referred to only takes place when the table itself is destructed.

@< Includes needed in the header file @>=
#include <map>
#include "lexer.h" // for the identifier hash table

@~We currently do not do overloading, so a simple associative table with the
identifier as key is used.

@< Type definitions @>=
class Id_table
{ typedef std::map<Hash_table::id_type,id_data> map_type;
  map_type table;
  Id_table(const Id_table&); // copying forbidden
  Id_table& operator=(const Id_table&); // assignment forbidden
public:
  Id_table() : table() @+{} // the default and only accessible constructor
  ~Id_table(); // destructor of all values referenced in the table
@)
  void add(Hash_table::id_type id, shared_value v, type_ptr t); // insertion
  bool remove(Hash_table::id_type id); // deletion
  shared_share address_of(Hash_table::id_type id); // locate
@)
  type_p type_of(Hash_table::id_type id) const; // lookup
  shared_value value_of(Hash_table::id_type id) const; // lookup
@)
  size_t size() const @+{@; return table.size(); }
  void print(std::ostream&) const;
};

@ As the table has strict ownership of the contained types, the destructor
must explicitly delete them.

@< Function def... @>=
Id_table::~Id_table()
{ for (map_type::const_iterator p=table.begin(); p!=table.end(); ++p)
  @/{@; delete p->second.type; }
}

@ The method |add| tries to insert the new mapping from the key |id| to a new
value-type pair. Doing so, it must distinguish two cases. It tentatively
inserts a new empty entry for the identifier into the table; this returns a
pair with a second boolean component telling whether a new entry was added
(the identifier was unknown as global identifier). If this is the case, we
continue to fill the slot with a newly allocated shared pointer and the
provided type. If the identifier was already present, we abandon the old
shared pointer location, resetting the pointer to it to point to a newly
allocated one, destroy the old type and insert the new type. All in all, the
only difference in the code of the two branches is the deletion of the old
type.

@< Function def... @>=
void Id_table::add(Hash_table::id_type id, shared_value val, type_ptr type)
{ std::pair<map_type::iterator,bool> trial=
    table.insert(std::make_pair(id,id_data()));
  id_data& slot=trial.first->second;
  if (trial.second) // a fresh global identifier
  {@; slot.val.reset(new shared_value(val));
    slot.type=type.release();
  }
  else
  {@; slot.val.reset(new shared_value(val));
    delete slot.type;
    slot.type=type.release();
  }
}

@ The |remove| method removes an identifier if present, and returns whether
this was the case. We should not forget to clean up the owned type that
disappears.

@< Function def... @>=
bool Id_table::remove(Hash_table::id_type id)
{ map_type::iterator p = table.find(id);
  if (p==table.end())
    return false;
  delete p->second.type; table.erase(p); return true;
}

@ In order to have |const| lookup methods, we must refrain from inserting into
the table if the key is not found; we return a null pointer in that case. The
pointer returned by |type_of| remains owned by the table. Although it does not
immediately modify the table, |address_of| is classified as a manipulator,
since the pointer returned may be used to modify a stored value.

@< Function def... @>=
type_p Id_table::type_of(Hash_table::id_type id) const
{@; map_type::const_iterator p=table.find(id);
  return p==table.end() ? NULL : p->second.type;
}
shared_value Id_table::value_of(Hash_table::id_type id) const
{ map_type::const_iterator p=table.find(id);
  return p==table.end() ? shared_value(value(NULL)) : *p->second.val;
}
shared_share Id_table::address_of(Hash_table::id_type id)
{ map_type::iterator p=table.find(id);
  if (p==table.end())
    throw std::logic_error @|
    (std::string("Identifier without value:")+main_hash_table->name_of(id));
@.Identifier without value@>
  return p->second.val;
}

@ We provide a |print| member that shows the contents of the entire table.
@< Function def... @>=

void Id_table::print(std::ostream& out) const
{ for (map_type::const_iterator p=table.begin(); p!=table.end(); ++p)
    out << main_hash_table->name_of(p->first) << ": " @|
        << *p->second.type << ": " << **p->second.val << std::endl;
}

std::ostream& operator<< (std::ostream& out, const Id_table& p)
@+{@; p.print(out); return out; }

@~We shouldn't forget to declare that operator, or it won't be found, giving
kilometres of error message.

@< Declarations of exported functions @>=
std::ostream& operator<< (std::ostream& out, const Id_table& p);

@ We declare just a pointer to the global identifier table here.
@< Declarations of global variables @>=
extern Id_table* global_id_table;

@~Here we set the pointer to a null value; the main program will actually
create the table.

@< Global variable definitions @>=
Id_table* global_id_table=NULL;

@*1 Global identifiers.
%
When during type checking an identifiers binds to a value in the global
identifier table, it will be converted into a |global_identifier| object.
Since a value is already available at this time, we can record the location of
the (pointer to the shared) value in the |global_identifier| object. Apart
from avoiding look-up at evaluation time, this measure also allows detaching
an applied global identifier from a newly assigned value if the latter should
have a different type, as was done in the |add| method of the identifier
table. This precaution ensures that evaluation will always result in a value
of the expected type.

@< Type definitions @>=
class global_identifier : public identifier
{ const shared_share address;
public:
  explicit global_identifier(Hash_table::id_type id);
  virtual ~global_identifier() @+ {}
  virtual void evaluate(level l) const;
};

@ The constructor for |global_identifier::evaluate| locates the value
associated to the identifier in the global identifier table.

@< Function definitions @>=
global_identifier::global_identifier(Hash_table::id_type id)
: identifier(id), address(global_id_table->address_of(id))
@+{}


@ Evaluating a global identifier returns the value stored in the location
|address|, possibly expanded if |l==multi_value|, or nothing at all if
|l==no_value|. However, since initially undefined global variables were added
to the language, we have to watch out for a (shared) null pointer at
|*address|.

@< Function definitions @>=
void global_identifier::evaluate(level l) const
{ if (address->get()==NULL)
  { std::ostringstream o;
    o << "Taking value of uninitialized variable " << name();
    throw std::runtime_error(o.str());
  }
  push_expanded(l,*address);
}

@*1 Local identifiers.
%
Local identifiers will be accessed from the current execution context, which
is a stack of variable bindings independent of the |execution_stack| (the
latter bieng used for anonymous components of expressions being evaluated).
This stack is implemented as a singly linked list, and accessed through a
(smart) pointer. That pointer is declared local to the evaluator; compared
with making it a parameter to the |evaluate| methods, this has the advantage
of not encumbering the numerous such methods that neither use nor modify the
context in any way (those not involving identifiers or user defined
functions).

@< Local var... @>=
context_ptr execution_context;

@ A disadvantege of using a static variable is that in case of exceptions it
retains the value current before throwing. Therefore we need to explicitly
reset the execution in such cases. Since it is a smart pointer, resetting
automatically takes care of adjusting reference counts and maybe deleting
values that are part of the discarded context.

@< Actions... @>=
execution_context.reset();

@ We derive the class of local identifiers from that of global ones, which
takes care of its |print| method.

@< Type definitions @>=
class local_identifier : public identifier
{ size_t depth, offset;
public:
  explicit local_identifier(Hash_table::id_type id, size_t i, size_t j)
     : identifier(id), depth(i), offset(j) @+{}
  virtual void evaluate(level l) const; // only this method is redefined
};

@ The method |local_identifier::evaluate| looks up a value in the
|execution_context|.

@< Function definitions @>=
void local_identifier::evaluate(level l) const
{@; push_expanded(l,execution_context->elem(depth,offset)); }

@ For an applied identifier, we first look in |id_context| for a binding of
the identifier, and if found it will be a local identifier, and otherwise we
look in |global_id_table|. If found in either way, the associated type must
equal the expected type (if any), or be convertible to it using |coerce|.
There is a subtlety in that the identifier may have a more general type than
|type| required by the context (for instance if it was \&{let} equal to an
empty list, and a concrete type of list is required). In this case
|specialise| succeeds without modifying |type| and we specialise the
identifier to |type| instead, so that it cannot be subsequently used with an
incompatible specialisation (notably any further assignments to the variable
must respect the more specific type).

@< Other cases for type-checking and converting... @>=
case applied_identifier:
{ type_p it; expression_ptr id; size_t i,j;
  if ((it=id_context->lookup(e.e.identifier_variant,i,j))!=NULL)
    id.reset(new local_identifier(e.e.identifier_variant,i,j));
  else if ((it=global_id_table->type_of(e.e.identifier_variant))!=NULL)
    id.reset(new global_identifier(e.e.identifier_variant));
  else throw program_error  @|
       (std::string("Undefined identifier ")
	+main_hash_table->name_of(e.e.identifier_variant));
@.Undefined identifier@>
  if (type.specialise(*it))
    {@; it->specialise(type); return id.release(); }
  else if (coerce(*it,type,id))
    return id.release();
  else throw type_error(e,copy(*it),copy(type));
}

@*1 Operator and function overloading.
%
While the simple mechanism of identifier identification given above suffices
for many purposes, it would be very restrictive in case of operators (since it
allows only one function, with fixes argument types, to be bound globally to
an operator symbol identifier), and to a somewhat lesser measure (since one
could vary the identifier name according to the argument types) for functions.
So we definitely want to allow operator overloading (defining the same
operator for different combinations of argument types), and with such a
mechanism in place, it is easy to allow function overloading as well, which
will for instance allow the intuitive convention of simply naming built-in
functions after the mathematical meaning of the result they compute, even if
such a result can be computed from different sets of input data.

To implement overloading we use a similar structure as the ordinary global
identifier table. However the basic table entry for an overloading needs a
level of sharing less, since the function bound for given argument types
cannot be changed by assignment, so the call will refer directly to the value
stored rather than to its location. We also take into account that the stored
types are always function types (this saves space, while for normal use we do
not need to access the full function type as a |type_expr|, something this
representation makes rather difficult). Remarks about ownership of the type
apply without change however.

@< Type definitions @>=

struct overload_data
{ shared_value val; @+ func_type_p type;
  overload_data() : val(),type(NULL)@+ {}
};

@ Looking up an overloaded identifier leads to a vector of value-type pairs.
This is preferable to using a |std::multimap| multi-mapping identifiers to
individual value-type pairs, as that would give us no control over the order
in which the pairs for the same identifier are ordered. This is important
since we want to try matching more specific (harder to convert to) argument
types before trying less specific ones. We envisage as only method for making
an identifier overloaded actually giving an overloaded definition, so the
table will not normally associate an empty vector to an identifier (though it
could do so after an entry is explicitly removed). The |variants| method will
signal absence of an identifier by returning an empty list of variants, and no
separate test for this condition is provided.

@< Type definitions @>=

class overload_table
{
public:
  typedef std::vector<overload_data> variant_list;
  typedef std::map<Hash_table::id_type,variant_list> map_type;
private:
  map_type table;
  overload_table(const Id_table&); // copying forbidden
  overload_table& operator=(const Id_table&); // assignment forbidden
public:
  overload_table() : table() @+{} // the default and only accessible constructor
  ~overload_table(); // destructor of all values referenced in the table
@) // accessors
  const variant_list& variants(Hash_table::id_type id) const;
  size_t size() const @+{@; return table.size(); }
  void print(std::ostream&) const;
@) // manipulators
  void add(Hash_table::id_type id, shared_value v, type_ptr t);
   // insertion
  bool remove(Hash_table::id_type id, const type_expr& arg_t); //deletion
};

@ As for the ordinary identifier table, the table owns the types, so the
destructor must clean then up.

@< Function definitions @>=

overload_table::~overload_table()
{ for (map_type::const_iterator p=table.begin(); p!=table.end(); ++p)
    for (size_t i=0; i<p->second.size(); ++i)
    {@; delete p->second[i].type; }
}

@ The |variants| method just returns a reference to the stored vector of
overload instances, or to a static empty vector if nothing is found.

@< Function definitions @>=
const overload_table::variant_list& overload_table::variants
  (Hash_table::id_type id) const
{ static const variant_list empty;
  map_type::const_iterator p=table.find(id);
  return p==table.end() ? empty : p->second;
}

@ The |add| method is what introduces and controls overloading. We first try
to associate an singleton vector with the identifier, and upon success insert
the given value-type pair. In the contrary case (the identifier already had a
vector associate to it), we must test the new pair against existing elements,
reject it is there is a conflicting entry present, and otherwise make sure it
is inserted before any strictly less specific overloaded instances.

@< Function def... @>=
void overload_table::add
  (Hash_table::id_type id, shared_value val, type_ptr tp)
{ if (tp->kind!=function_type)
    throw std::logic_error("Wrapper function has non-function type");
  func_type_ptr type(tp->func);
  tp->kind=undetermined_type; // release
  if (type->arg_type==void_type)
    throw program_error("Cannot define function overload without arguments");
  std::pair<map_type::iterator,bool> trial=
    table.insert(std::make_pair(id,variant_list(1,overload_data())));
  variant_list& slot=trial.first->second;
  if (trial.second) // a fresh overloaded identifier
  {@; slot[0].val=val;
    slot[0].type=type.release();
  }
  else
  @< Compare |type| against entries of |slot|, if none are close then add
  |val| and |type| at the end, if any is close without being one-way
  convertible to or from it throw an error, and in the remaining case make
  sure |type| is added after any types that convert to it and before any types
  it converts to @>
}

@ We call |is_close| for each existing argument type; if it returns a nonzero
value it must be either |0x6|, in which case insertion must be after that
entry, or |0x5|, in which case insertion must be no later than at this
position, so that the entry in question (after shifting forward) stays ahead.
The last of the former cases and the first of the latter are recorded, and
their requirements should be compatible.

Although the module name does not mention it, we allow one case of close and
mutually convertible types, namely identical types; in this case we simple
replace the old definition for this type by the new one. This could still
change the result type, but that does not matter because if any calls that
were type-checked against the old definition should survive (in a closure),
they have been also bound to the (function) \emph{value} that was previously
accessed by that definition, and will continue to use it; their operation is
in no way altered by the replacement of the definition.

@< Compare |type| against entries of |slot|... @>=
{ size_t lwb=0; size_t upb=slot.size();
  for (size_t i=0; i<slot.size(); ++i)
  { unsigned int cmp= is_close(type->arg_type,slot[i].type->arg_type);
    switch (cmp)
    {
      case 0x6: lwb=i+1; break; // existent type |i| converts to |type|
      case 0x5: @+ if (upb>i) upb=i; @+ break; // |type| converts to type |i|
      case 0x7:
        if (slot[i].type->arg_type==type->arg_type)
        @/{@; slot[i].val=val;
            delete slot[i].type;
            slot[i].type=type.release();
            return;
          }
      @/// |else| {\bf fall through}
      case 0x4: // conflicting cases
        { std::ostringstream o;
          o << "Cannot overload `" << main_hash_table->name_of(id) << "', " @|
               "previous type " << slot[i].type->arg_type
            << " is too close to "@| << type->arg_type
            << ",\nmaking overloading potentially ambiguous." @|
               " Priority cannot\ndisambiguate, as "
            << (cmp==0x4 ? "neither" :"either") @|
            << " type converts to the other";
          throw program_error(o.str());
        }
      default: @+{} // nothing for unrelated argument types
    }
  }
  if(lwb>upb)
    throw std::logic_error("Conflicting order related types");
  @< Insert |val| and |type| at entry |upb| after shifting remainder up @>
}

@ Shifting entries during a call of |insert| will probably create duplicated
type pointers at some points; however these do not risk double deletion since
no errors can be thrown at such points. Indeed, the only point where throwing
may occur is when extending the vector which happens at the beginning, and in
case of successful reallocation entries are copied (and temporarily
duplicated) anyway; the |insert| method involves nothing more dangerous than
this.

@< Insert |val| and |type| at entry |upb| after shifting remainder up @>=
{ slot.insert(slot.begin()+upb,overload_data());
@/slot[upb].val=val; slot[upb].type=type.release();
}


@ The |remove| method allows removing an entry from the overload table, for
instance to make place for another one. It returns a boolean telling whether
any such binding was found (and removed). The |variants| array might become
empty, but remains present and will be reused upon future additions.

@< Function def... @>=
bool overload_table::remove(Hash_table::id_type id, const type_expr& arg_t)
{ map_type::iterator p=table.find(id);
  if (p==table.end()) return false;
  variant_list& variants=p->second;
  for (size_t i=0; i<variants.size(); ++i)
    if (variants[i].type->arg_type==arg_t)
    @/{@; delete variants[i].type;
      variants.erase(variants.begin()+i);
      return true;
    }
  return false;
}

@ We provide a |print| member that shows the contents of the entire table,
just like for identifier tables. Only this one prints multiple entries per
identifier.

@< Function def... @>=

void overload_table::print(std::ostream& out) const
{ type_expr type; type.kind=function_type;
  for (map_type::const_iterator p=table.begin(); p!=table.end(); ++p)
    for (size_t i=0; i<p->second.size(); ++i)
    { type.func = p->second[i].type;
      out << main_hash_table->name_of(p->first) << ": " @|
        << type << ": " << *p->second[i].val << std::endl;
    }
  type.kind=undetermined_type; // avoid destruction of |type.func|
}

std::ostream& operator<< (std::ostream& out, const overload_table& p)
{@; p.print(out); return out; }

@~We shouldn't forget to declare that operator, if we want to use it.

@< Declarations of exported functions @>=
std::ostream& operator<< (std::ostream& out, const overload_table& p);

@ We introduce a single overload table in the same way as the global
identifier table.

@< Declarations of global variables @>=
extern overload_table* global_overload_table;

@~Here we set the pointer to a null value; the main program will actually
create the table.

@< Global variable definitions @>=
overload_table* global_overload_table=NULL;

@ Overloading resolution will be called from the case in |convert_expr| for
function applications, after testing that overloads exist, so we can transmit
the relevant |variants| together with the other parameters.

@< Declarations of exported functions @>=
expression resolve_overload
  (const expr& e,
   type_expr& type,
   const overload_table::variant_list& variants);

@~To resolve overloading, we used to plunge into each variant, catching and
ignoring errors, until one succeeded without error. For long formulae this
traverses a search tree of exponential size, giving unacceptably inefficient
expression analysis. So instead we now first try to find a matching variant,
using an \foreign{a priori} type of the operand(s). This implies that the
operand must be correctly typed without the benefit of a known result type,
but the type found need not be an exact match with the one specified in the
variant. Our matching condition is that |is_close| should hold, with the bit
set that indicates a possible conversion from |a_priori_type| to the operand
type for which the variant is defined. Coercions may need to be inserted in
the operand expression, and since that expression could be arbitrarily
complex, inserting coercions explicitly after the fact would be very hard to
program. We choose the easier approach is to cast away the converted
expression once the matching variant is found, and redo the analysis with the
now known result type so that coercions get inserted during conversion. But
then we again risk exponential time (although with powers of~2 rather than
powers of the number of variants); since probably coercions will not be needed
at all levels, we mitigate this risk by not redoing any work in case of an
exact match.

Apart from those in |variants|, we also test for certain argument types that
will match without being in any table; for instance the size-of operator~`\#'
can be applied to any row type to give its number of components. Being more
generic bindings, we test for them after the more specific ones fail. The
details of these cases, like those of the actual construction of a call for a
matching overloaded function, will be given later.

@< Function definitions @>=
expression resolve_overload
  (const expr& e,
   type_expr& type,
   const overload_table::variant_list& variants)
{ const expr& args = e.e.call_variant->arg;
  type_expr a_priori_type;
  expression_ptr arg(convert_expr(args,a_priori_type));
  Hash_table::id_type id =  e.e.call_variant->fun.e.identifier_variant;
  @< If a special operator like size-of matches |id|, |return| a call of it
     with argument |args| @>
  for (size_t i=0; i<variants.size(); ++i)
  { const overload_data& v=variants[i];
    if ((is_close(a_priori_type,v.type->arg_type)&0x1)!=0)
      // could first convert to second?
    { if (a_priori_type!=v.type->arg_type)
        arg=expression_ptr(convert_expr(args,v.type->arg_type));
          // redo conversion
      @< Return a call of variant |v| with argument |arg|, or |throw| if
         result type mismatches |type| @>
    }
  }

  @< Complain about failing overload resolution @>
}

@ Failing overload resolution causes a |program_error| explaining the
is matching identifier and type. Since most function definitions will be in
the overload table even when only one definition is present, we produce a
|type_error| in that case, so that the message will mention the unique
expected argument type.

@< Complain about failing overload resolution @>=
if (variants.size()==1)
  throw type_error(args,copy(a_priori_type),copy(variants[0].type->arg_type));
else
{ std::ostringstream o;
  o << "Failed to match `"
    << main_hash_table->name_of(id) @|
    << "' with argument type "
    << a_priori_type;
  throw expr_error(e,o.str());
}

@* Function calls.
%
One of the most basic tasks of the evaluator is to allow function calls, which
may involve either buit-in or user-defined functions. We start with
introducing a type for representing general function calls after type
checking; the function is not assumes to be given by an identifier, and if it
should anyway, it is a non-overloaded call that dynamically takes the value
bound to the (possibly local) identifier.

@< Type def... @>=
struct call_expression : public expression_base
{ expression function, argument;
@)
  call_expression(expression_ptr f,expression_ptr a)
   : function(f.release()),argument(a.release()) @+{}
  virtual ~call_expression() @+ {@; delete function; delete argument; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ To print a function call we print the function expression in parentheses,
and the argument, latter enclosed in parentheses unless it is a tuple
expression (which already has parentheses), which condition is tested by a
dynamic cast.

@< Function definitions @>=
void call_expression::print(std::ostream& out) const
{ if (dynamic_cast<identifier*>(function)!=NULL) out << *function;
  else out << '(' << *function << ')';
  if (dynamic_cast<tuple_expression*>(argument)!=NULL) out << *argument;
  else out << '(' << *argument << ')';
}

@ The evaluation of the call of a built-in function executes a ``wrapper
function'', that usually consists of a call to a library function sandwiched
between unpacking and repacking statements; in some simple cases a wrapper
function may decide to do the entire job itself.

The arguments and results of wrapper functions will be transferred from and to
stack as a |shared_value|, so a wrapper function has neither arguments nor a
result type. Thus variables that refer to a wrapper function have the type
|wrapper_function| defined below; the |level| parameter serves the same
function as for |evaluate| methods, to inform whether a result value should be
produced at all, and if so whether it should be expanded on the
|execution_stack| in case it is a tuple. We shall need to bind values of this
type to identifiers representing built-in functions, so we derive an
associated ``primitive type'' from |value_base|.

@< Type definitions @>=
typedef void (* wrapper_function)(expression_base::level);
@)
struct builtin_value : public value_base
{ wrapper_function val;
  std::string print_name;
@)
  builtin_value(wrapper_function v,const std::string& n)
  : val(v), print_name(n) @+ {}
  virtual void print(std::ostream& out) const
  @+{@; out << '{' << print_name << '}'; }
  builtin_value* clone() const @+{@; return new builtin_value(*this); }
  static const char* name() @+{@; return "built-in function"; }
private:
  builtin_value(const builtin_value& v)
  : val(v.val), print_name(v.print_name)
  @+{}
};

@ While syntactically more complicated than ordinary function calls, the call
of overloaded functions is actually simpler at run time, because the function
is necessarily referred to by an identifier instead of by an arbitrary
expression, and overloading resolution results in a function \emph{value}
rather than in the description of a location where the function can be found
at run time. If that value happens to be a built-in function, the call will be
translated into an |overloaded_builtin_call| rather than into a
|call_expression|.

@< Type definitions @>=
struct overloaded_builtin_call : public expression_base
{ wrapper_function f;
  std::string print_name;
  expression argument;
@)
  overloaded_builtin_call(wrapper_function v,const char* n,expression_ptr a)
  : f(v), print_name(n), argument(a.release())@+ {}
  virtual ~overloaded_builtin_call() @+ {@; delete argument; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ When printing, we ignore the stored wrapper function (which does not record
its name) and use the overloaded function name; otherwise we proceed as for
general function calls with an identifier as function.

@< Function definitions @>=
void overloaded_builtin_call::print(std::ostream& out) const
{ out << print_name;
  if (dynamic_cast<tuple_expression*>(argument)!=NULL) out << *argument;
  else out << '(' << *argument << ')';
}

@ Some built-in functions like |print| accept arguments of any types, and in
particular tuples of any length. For such functions we cannot adopt the method
used for other built-in functions of expanding argument tuples on the stack,
since there would then be no way to recover their number. Fortunately such
functions are necessarily accessed through overloading, so we detect their use
at analysis time, and record it in the type off call expression generated.
Therefore we derive a type from |overloaded_builtin_call| that will override
only the |evaluate| method.

@< Type definitions @>=
struct generic_builtin_call : public overloaded_builtin_call
{ typedef overloaded_builtin_call base;
@)
  generic_builtin_call(wrapper_function v,const char* n,expression_ptr a)
  : base(v,n,a)@+ {}
  virtual void evaluate(level l) const;
};

@*1 Type-checking function calls.
%
The function in a call can be any type of expression; in case of a non-local
identifier for which overloads are defined we attempt overload resolution (and
ignore any value possibly present in the global identifier table). Otherwise
the function expression determines its own type, and once this is known, its
argument and result types can be used to help converting the argument
expression and the call expression itself. Thus we first get the type of the
expression in the function position, requiring only that it be a function
type, then type-check and convert the argument expression using the obtained
result type, and build a converted function call~|call|. Finally we test if
the required type matches the return type (in which case we simply
return~|call|), or if the return type can be coerced to it (in which case we
return |call| as transformed by |coerce|); if neither is possible we throw
a~|type_error|.

@< Other cases for type-checking and converting... @>=
case function_call:
{ if (e.e.call_variant->fun.kind==applied_identifier)
    @< Convert and |return| an overloaded function call if
    |e.e.call_variant->fun| is not a local identifier and is known in
    |global_overload_table| @>
  type_ptr f_type=copy(gen_func_type); // start with generic function type
  expression_ptr fun(convert_expr(e.e.call_variant->fun,*f_type));
  expression_ptr arg
    (convert_expr(e.e.call_variant->arg,f_type->func->arg_type));
  expression_ptr call (new call_expression(fun,arg));
  return conform_types(f_type->func->result_type,type,call,e);
}

@ The main work here has been relegated to |resolve_overload|; otherwise we
just need to take care of the things mentioned in the module name. In fact
there is one more case where overload resolution is not invoked than that name
mentions, namely when the argument expression is an empty tuple display, since
overloading with void argument type is forbidden. But this special treatment
of empty arguments makes a global value with type function-without-arguments
almost behave like an overloaded instance; the user should just avoid
providing a nonempty argument of void type, which is bad practice anyway.

The cases relegated to |resolve_overload| include calls of special operators
like the size-of operator~`\#', even if such an operator should not occur in
the overload table.

@< Convert and |return| an overloaded function call... @>=
{ const Hash_table::id_type id =e.e.call_variant->fun.e.identifier_variant;
  const expr arg=e.e.call_variant->arg;
  size_t i,j;
  if (not is_empty(arg) and id_context->lookup(id,i,j)==NULL)
  { const overload_table::variant_list& variants
      = global_overload_table->variants(id);
    if (variants.size()>0 or is_special_operator(id))
      return resolve_overload(e,type,variants);
  }
}

@ For overloaded function calls, once the overloading is resolved, we proceed
in a similar fashion to non-overloaded calls, except that there is no function
expression to convert (the overload table contains an already evaluated
function value, either built-in or user-defined). We deal with the built-in
case here, and will give the user-defined case later when we have discussed
the necessary value types.

@< Return a call of variant |v|... @>=
{ expression_ptr call;
  builtin_value* f = dynamic_cast<builtin_value*>(v.val.get());
  if (f!=NULL)
    call = expression_ptr
      (new overloaded_builtin_call(f->val,f->print_name.c_str(),arg));
  else @< Set |call| to the call of the user-defined function |v| @>
  return conform_types(v.type->result_type,type,call,e);
}

@ The names of special operators are tested for each time analyse an
overloaded call; to avoid having to look them up in |main_hash_table| each
time, we store each one in a static variable inside a local function.

@< Local function definitions @>=
Hash_table::id_type size_of_name()
{@; static Hash_table::id_type name=main_hash_table->match_literal("#");
  return name;
}
Hash_table::id_type print_name()
{@; static Hash_table::id_type name=main_hash_table->match_literal("print");
  return name;
}
Hash_table::id_type prints_name()
{@; static Hash_table::id_type name=main_hash_table->match_literal("prints");
  return name;
}
@)
inline bool is_special_operator(Hash_table::id_type id)
{@; return id==size_of_name()
        or id==print_name()
        or id==prints_name(); }

@ For operator symbols that satisfy |is_special_operator(id)|, we test generic
argument type patterns before we test instances in the overload table, because
the latter could otherwise mask some generic ones due to coercion. Therefore
if we fail to find a match, we simply fall through; however if we match an
argument type but fail to match the returned type, we throw a |type_error|.


@: sizeof section @>

@< If a special operator like size-of... @>=
{ if (id==size_of_name())
  { if (a_priori_type.kind==row_type)
    { expression_ptr call(new overloaded_builtin_call(sizeof_wrapper,"#",arg));
      return conform_types(int_type,type,call,e);
    }
    else if (a_priori_type.kind!=undetermined_type and
             a_priori_type.specialise(pair_type))
    @< Recognise and return 2-argument versions of `\#', or fall through in
       case of failure @>
  }
  else if (id==print_name() or id==prints_name()) // these always match
  { expression c = id==print_name() @|
      ? new generic_builtin_call(print_wrapper,"print",arg) @|
      : new generic_builtin_call(prints_wrapper,"prints",arg);
    expression_ptr call(c); // get ownership
    if (type.specialise(void_type))
      return call.release();
    throw type_error(e,copy(type),copy(void_type));
  }

}

@ For dyadic use of the operator `\#' we shall encounter the somewhat unusual
situation that we have already converted the entire argument expression at the
point where we discover, based on the type found, that one of the arguments
might need to be coerced. This means that such a coercion must be inserted
into an already constructed expression, whereas usually it is applied on the
outside of an expression under construction (notably in ordinary overloading
if coercion is needed we simply convert the operands a second time, inserting
the coercions while doing so). Concretely this means that although we can use
the |coerce| function, it wants a reference to an auto-pointer for the
expression needing modification (so that it can manage ownership during its
operation), but the expression here is held in an ordinary pointer inside a
|tuple_expression|. Therefore we must copy the ordinary |expression| pointer
temporarily to an |expression_ptr|, simulating the auto-pointer actions
manually: after construction of the |expression_ptr| we set the |expression|
(temporarily) to~|NULL| to avoid potential double destruction, and after the
call to |coerce| we release the (possibly modified) |expression_ptr| back into
the |expression|.

Another complication is that we decided having a dyadic use of `\#' based on
finding a 2-tuple type, but this is no guarantee there are actually two
operand subexpressions (in a |tuple_expression|); we do a dynamic cast to find
that out, and in the case the user was so contrived as to use monadic `\#' on
a non-tuple expression of 2-tuple type, we just report that no coercion of
operands is done (after all we need a subexpression to be able to insert any
conversion). These complications warrant defining a separate function to
handle them.

@< Local function definitions @>=
bool can_coerce_arg
  (expression e,size_t i,const type_expr& from,const type_expr& to)
{ tuple_expression* tup= dynamic_cast<tuple_expression*>(e);
  if (tup==NULL) return false;
  expression_ptr comp(tup->component[i]); tup->component[i]=NULL;
  bool result = coerce(from,to,comp); tup->component[i]=comp.release();
  return result;
}

@ The operator `\#' can be used also as infix operator, to join (concatenate)
two row values of the same type or to extend one on either end by a single
element. In the former case we require that both arguments have identical row
type, in the latter case we allow the single element to be converted to the
component type of the row value.

There is a subtlety in the order here, due to the fact that one of the
arguments could be the empty list, with undetermined row type. Then there is
actual ambiguity: with `\#' interpreted as join, the result is just the other
operand, while suffixing or prefixing to the empty list gives a singleton of
the other operand, and to some operands the empty list itself can also be
suffixed or prefixed. The simplest resolution of this ambiguity is to say that
join never applies with one of the arguments of undetermined list type, and
that suffixing is preferred over prefixing (for instance in $[\,]\#[[2]]$).
This can be obtained by testing for suffixing before testing for join: after
the former test we know the first argument does not have type \.{[*]}, so it
will match the second argument type only if both are determined row types.

@< Recognise and return 2-argument versions of `\#'... @>=
{ type_expr& arg_tp0 = a_priori_type.tuple->t;
  type_expr& arg_tp1 = a_priori_type.tuple->next->t;
  if (arg_tp0.kind==row_type)
  { if (arg_tp0.component_type->specialise(arg_tp1) or @|
        can_coerce_arg(arg.get(),1,arg_tp1,*arg_tp0.component_type)) // suffix
    { expression_ptr call(new overloaded_builtin_call
        (suffix_element_wrapper,"#",arg));
      return conform_types(arg_tp0,type,call,e);
    }
    if (arg_tp0==arg_tp1) // join
    { expression_ptr call(new overloaded_builtin_call
        (join_rows_wrapper,"#",arg));
      return conform_types(arg_tp0,type,call,e);
    }
  }
  if (arg_tp1.kind==row_type and @|
         (arg_tp1.component_type->specialise(arg_tp0) or @|
          can_coerce_arg(arg.get(),0,arg_tp0,*arg_tp1.component_type)))
          // prefix
  { expression_ptr call(new overloaded_builtin_call
      (prefix_element_wrapper,"#",arg));
    return conform_types(arg_tp1,type,call,e);
  }
}

@*1 Evaluating built-in function calls.
%
To evaluate a |call_expression| object we evaluate the function, and then
test whether it is a built-in function. In the former case we evaluate the
arguments expanded on th stack and call the built-in function, passing the
|level| parameter so that if necessary the call can in its turn return and
expanded result (or no result at all). If not a built-in function, it must be
a user-defined function, whose execution will be detailed later, but in this
case it will be more useful to have the argument as a single value. As a
general mechanism to aid locating errors, we signal if an error was produced
during the evaluation of a function call, but we make sure the evaluation of
the arguments(s) is done outside this |try| block, since reporting functions
that have not yet started executing would be confusing.

@< Function definitions @>=
void call_expression::evaluate(level l) const
{ function->eval(); @+ shared_value fun=pop_value();
@/builtin_value* f=dynamic_cast<builtin_value*>(fun.get());
  argument->evaluate(f==NULL ? single_value : multi_value);
  try
  { if (f==NULL)
      @< Call user-defined function |fun| with argument on |execution_stack| @>
    else // built-in functions
      (*f->val)(l); // call the wrapper function, handling |l| appropriately
  }
  @< Catch-block for exceptions thrown within function calls @>
}

@ Although we catch all |std::exception| errors thrown during the execution of
a function call, we only report it if the function was referred to by an
identifier, for otherwise the message would become too messy. In the mentioned
case we tack a line with the function name to the error string, and re-throw
the error. This will result in a traceback, inner to outer, of interrupted
(non anonymous) function calls. Different types of error could be concerned:
|std::runtime_error| thrown in some |evaluate| method or built-in function is
the most common case, but there could be some |std::logic_error| as well, and
|std::bad_alloc| is also always a possibility. Since the base class
|std::exception| provides no means to influence the message produced by the
|what| method, nor does |std::bad_alloc|, we necessarily have to relabel the
latter as |std::runtime_error| in order to extend the error message. However
we can maintain the distinction between a |logic_error| and a |runtime_error|
using a dynamic cast.

@< Catch-block for exceptions thrown within function calls @>=
catch (const std::exception& e)
{ identifier* p=dynamic_cast<identifier*>(function);
  if (p!=NULL) // named function
  {
    const std::logic_error* l_err= dynamic_cast<const std::logic_error*>(&e);
    if (l_err!=NULL)
      throw std::logic_error
        (std::string(e.what())+"\n(in call of "+p->name()+')');
    throw std::runtime_error
      (std::string(e.what())+"\n(in call of "+p->name()+')');
  }
  throw; // for anonymous function calls, just rethrow the error unchanged
}

@*1 Evaluating overloaded built-in function calls.
%
Calling an overloaded built-in function calls the wrapper function after
evaluating the argument(s). We provide the same trace of interrupted functions
by temporarily catching errors as in the case of non-overloaded function
calls, but here we need not test that the call was one of a named function.

@< Function definitions @>=
void overloaded_builtin_call::evaluate(level l) const
{ argument->multi_eval();
  try
  {@; (*f)(l); }
  catch (const std::exception& e)
  { const std::logic_error* l_err= dynamic_cast<const std::logic_error*>(&e);
    if (l_err!=NULL)
      throw std::logic_error
        (std::string(e.what())+"\n(in call of "+print_name+')');
    throw std::runtime_error
      (std::string(e.what())+"\n(in call of "+print_name+')');
  }
}

@ For generic built-in functions like |print| we only change the fact that
arguments are evaluated using |eval| to a single value on the stack.

@< Function definitions @>=
void generic_builtin_call::evaluate(level l) const
{ argument->eval();
  try
  {@; (*f)(l); }
  catch (const std::exception& e)
  { const std::logic_error* l_err= dynamic_cast<const std::logic_error*>(&e);
    if (l_err!=NULL)
      throw std::logic_error
        (std::string(e.what())+"\n(in call of "+print_name+')');
    throw std::runtime_error
      (std::string(e.what())+"\n(in call of "+print_name+')');
  }
}

@* Let-expressions.
%
We shall now consider a simple type of expression in which local variables
occur, the let-expression. It is equivalent to an anonymous function
($\lambda$-expression) applied to the expression(s) in the let-declarations,
but no types need to be declared for the function parameters, since the types
of the corresponding expressions can be used for this. Nevertheless, we shall
in converting expressions to internal form forget the syntactic origin of the
expression, and translate to an application of an anonymous function, giving
us an occasion to introduce such functions as a new form of~|expression|.

@ We prepare the definition of $\lambda$-expression with the introduction of
auxiliary types, needed to deal with the general patterns by which formal
function parameters can be given. The parser produces values accessed by
pointers to |struct id_pat@;| to describe such patterns. The structure can be
used directly in a $\lambda$-expression, but we must handle some questions of
ownership. The value produced by the parser will be destroyed after the
command is processed, at which time the $\lambda$-expression could still
exist, so we cannot simply use a non-owned pointer. Moreover a
$\lambda$-expression can be evaluated one or more times, yielding ``closure''
values the might outlive the $\lambda$-expression, so appropriate duplication
or sharing must be organised. We opt for sharing, which will be done by a
|shared_ptr| supplied with a specialised deleter, since the sharing is done at
the outermost level only, but complete destruction requires calling
|destroy_id_pat| (defined in the file \.{parsetree.w}). To get this deleter
without runtime storage, we wrap it into a zero-size structure.

@< Type def... @>=
struct id_pat_deleter
{@; void operator()(id_pat* p) @+{@; destroy_id_pat(p); delete p; }};
typedef std::tr1::shared_ptr<id_pat> shared_pattern;
  // deleter type does not enter into this

@ We must also treat the question of obtaining ownership of an |id_pat|
structure. An initial node to be shared must certainly be allocated (the
parser never excutes |new idpat|), but we might steal a possible |sublist|
field, replacing it in the original by a null pointer to avoid double
destruction. However doing a deep copy is a cleaner solution, which avoids
modifying the parsed expression and isolates us from the allocation policy
used in the parser, which might change. To provide exception safety during the
copy, it seems for once easier to explicitly catch and clean up than to
introduce an intermediate class only for exception safety. At each point where
an exception might be thrown the argument to |destroy_id_pat| is properly
|NULL|-terminated.

@< Local function def... @>=
id_pat copy_id_pat(const id_pat& p)
{
  id_pat result=p; // shallow copy
  try
  { if ((p.kind&0x2)!=0)
    { result.sublist=NULL;
      for (patlist s=p.sublist,*d=&result.sublist;
           s!=NULL; s=s->next,d=&(*d)->next)
      @/{@;
        *d=new pattern_node; (*d)->next=NULL; (*d)->body=copy_id_pat(s->body);
      }
    }
  }
  catch (...) @+{@; destroy_id_pat(&result); throw; }
  return result;
}

@ Now we can define our $\lambda$-expression to use a |shared_pattern|; the
body is also shared.

@< Type def... @>=
struct lambda_expression : public expression_base
{ const shared_pattern param;
  shared_expression body;
@)
  lambda_expression(const id_pat& p, expression_ptr b);
  virtual ~lambda_expression() @+{} // subobjects do all the work
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ The main constructor cannot be inside the class definition, as it requires
the local function |copy_id_pat|. It creates a new node at the head of
|param|, which will henceforth be shared, and fills it with a deep copy. For
the body we create sharing as well, which is simpler since the passed
auto-pointer already gives us ownership.

@< Function def... @>=
inline
lambda_expression::lambda_expression(const id_pat& p, expression_ptr b)
: param(new id_pat(copy_id_pat(p)),id_pat_deleter()) , body(b.release())
@+{}

@ To print an anonymous function, we print the parameter, enclosed in
parentheses if the full parameter is named, followed by a colon and by the
function body. The parameter list cannot include types with the current setup,
as they are not explicitly stored after type analysis.

@< Function definitions @>=
void lambda_expression::print(std::ostream& out) const
{ if (param.get()==NULL)
    out << "()";
  else if ((param->kind&0x1)!=0)
    out << '(' << *param << ')';
  else out << *param;
  out << ": " << *body;
}

@ For handling declarations with patterns as left hand side, we need a
corresponding type pattern; for instance \\{whole}:$(x,,z:(f,))$ requires the
type \.{(*,*,(*,*))}. These recursive functions construct such types.

@< Local function def... @>=
type_list_ptr pattern_list(const patlist p);
type_ptr pattern_type(const id_pat& pat)
{@; return (pat.kind&0x2)==0
  ? copy(unknown_type)
  : make_tuple_type(pattern_list(pat.sublist));
}
@)
type_list_ptr pattern_list(const patlist p)
{@; return p==NULL ? type_list_ptr(NULL)
  : make_type_list(pattern_type(p->body),pattern_list(p->next));
}

@ We shall need some other functions to deal with patterns, all with a
similar structure. Here we count the number or list the identifiers in a
pattern.

@< Local function def... @>=
size_t count_identifiers(const id_pat& pat)
{ size_t result= pat.kind & 0x1; // 1 if |pat.name| is defined, 0 otherwise
  if ((pat.kind & 0x2)!=0) // then a list of subpatterns is present
    for (patlist p=pat.sublist; p!=NULL; p=p->next)
      result+=count_identifiers(p->body);
  return result;
}

void list_identifiers(const id_pat& pat, std::vector<Hash_table::id_type>& d)
{ if ((pat.kind & 0x1)!=0)
    d.push_back(pat.name);
  if ((pat.kind & 0x2)!=0) // then a list of subpatterns is present
    for (patlist p=pat.sublist; p!=NULL; p=p->next)
      list_identifiers(p->body,d);
}

@ Here we do a similar traversal, using a type of the proper structure,
pushing pairs onto a |bindings|.

@< Local function def... @>=
void thread_bindings
(const id_pat& pat,const type_expr& type, bindings& dst)
{ if ((pat.kind & 0x1)!=0) dst.add(pat.name,copy(type));
  if ((pat.kind & 0x2)!=0)
  { assert(type.kind==tuple_type);
    type_list l=type.tuple;
    for (patlist p=pat.sublist; p!=NULL; p=p->next,l=l->next)
      thread_bindings(p->body,l->t,dst);
  }
}

@ Finally, at runtime we shall perform a similar manipulation with an
appropriate |shared_value|.

@< Local function def... @>=
void thread_components
(const id_pat& pat,const shared_value& val, std::vector<shared_value>& dst)
{ if ((pat.kind & 0x1)!=0) dst.push_back(val);
  if ((pat.kind & 0x2)!=0)
  { tuple_value* t=force<tuple_value>(val.get());
    size_t i=0;
    for (patlist p=pat.sublist; p!=NULL; p=p->next,++i)
      thread_components(p->body,t->val[i],dst);
  }
}

@ To convert a let-expression, we first deduce the type of the declared
identifiers from the right hand side of its declaration, then set up new
bindings for those identifiers with the type found, and finally convert the
body to the required type in the extended context. Note that the constructed
|bindings| remains a local variable, but it is temporarily preprended to the
context by calling |push| and |pop| with the pointer variable |id_context|.
The |expression| obtained from converting the body is first turned into a
|lambda_expression|, and then an application of that expression to the
argument is produced and returned.

@< Other cases for type-checking and converting... @>=
case let_expr:
{ let lexp=e.e.let_variant;
  id_pat& pat=lexp->pattern;
  type_ptr decl_type=pattern_type(pat);
  expression_ptr arg(convert_expr(lexp->val,*decl_type));
  size_t n_id=count_identifiers(pat);
@/bindings new_bindings(n_id);
  thread_bindings(pat,*decl_type,new_bindings);
  new_bindings.push(id_context);
  expression_ptr body(convert_expr(lexp->body,type));
  new_bindings.pop(id_context);
  expression_ptr func(new lambda_expression(pat,body));
  return new call_expression(func,arg);
}

@ Before we can dicuss the evaluation of user-defined functions, we need to
introduce a type for the intermediate value produced by the anonymous
function, before it is actually called. Such a value is traditionally called a
closure, and it contains (a reference to) the expression body, as well as the
evaluation context current at the point the anonymous function is produced.

@< Type def... @>=
struct closure_value : public value_base
{ context_ptr cont;
  shared_pattern param;
   // used in evaluation only to count arguments
  shared_expression body;
@)
  closure_value@|(context_ptr c,
                  const shared_pattern& p,
                  const shared_expression& b) : cont(c), param(p), body(b) @+{}
  void print(std::ostream& out) const;
  closure_value* clone() const @+
  {@; return new closure_value(cont,param,body); }
  static const char* name() @+{@; return "closure"; }
};
typedef std::auto_ptr<closure_value> closure_ptr;
typedef std::tr1::shared_ptr<closure_value> shared_closure;

@ For now a closure prints just like the |lambda_expression| from which it was
obtained. One could imagine printing after this body ``where'' followed by the
bindings held in the |context| field. Even better only the bindings for
relevant (because referenced) identifiers could be printed.

@< Function def... @>=
void closure_value::print(std::ostream& out) const
{ if (param.get()==NULL)
    out << "()";
  else if ((param->kind&0x1)!=0)
    out << '(' << *param << ')';
  else out << *param;
  out << ": " << *body;
}

@ Evaluating a $\lambda$-expression just forms a closure and returns that.

@< Function def... @>=
void lambda_expression::evaluate(level l) const
{ if (l!=no_value)
  @/{@;closure_ptr result(new closure_value(execution_context,param,body));
     push_value(result);
  }
}

@ A call of a user-defined function passes through the same code as that of a
builtin function; after all, the type check does not make a difference between
the two kinds, so the distinction can only be made by a dynamic test during
evaluation (which test was already presented). After the test we come to the
code below, which in essence a call-by-value $\lambda$-calculus
evaluator. We must now have a closure as function value, and its evaluation
just temporarily replaces the current execution context from the one stored in
the closure, pushes a new frame defined by the argument and the evaluates the
function body.

@< Call user-defined function |fun| with argument on |execution_stack| @>=
{ closure_value* f=force<closure_value>(fun.get());
@)std::vector<shared_value> new_frame;
  new_frame.reserve(count_identifiers(*f->param));
  thread_components(*f->param,pop_value(),new_frame);
@)
  context_ptr saved_context(execution_context);
  execution_context.reset(new context(f->cont,new_frame));
  f->body->evaluate(l); // pass evaluation level |l| to function body
  execution_context = saved_context;
}

@ For function overloads given by a user-defined function, we need a new
expression type which is capable of storing a closure value.

@< Type definitions @>=
struct overloaded_closure_call : public expression_base
{ shared_closure fun;
  std::string print_name;
  expression argument;
@)
  overloaded_closure_call
   (shared_closure f,const std::string& n,expression_ptr a)
  : fun(f), print_name(n), argument(a.release())@+ {}
  virtual ~overloaded_closure_call() @+ {@; delete argument; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ When printing it, we ignore the closure and use the overloaded function
name.

@< Function definitions @>=
void overloaded_closure_call::print(std::ostream& out) const
{ out << print_name;
  if (dynamic_cast<tuple_expression*>(argument)!=NULL) out << *argument;
  else out << '(' << *argument << ')';
}

@ The identification of this case is done inside |resolve_overload|, after
testing that the value bound is not a built-in function. Since an |overload|
table should hold only values of function type, we must have a |closure_value|
if it was not a |builtin_value|.

@< Set |call| to the call of the user-defined function |v| @>=
{ shared_closure fun =
   std::tr1::dynamic_pointer_cast<closure_value>(v.val);
  if (fun==NULL)
    throw std::logic_error("Overloaded value is not a function");
  std::ostringstream name;
  name << main_hash_table->name_of(id) << '@@' << v.type->arg_type;
  call = expression_ptr (new overloaded_closure_call(fun,name.str(),arg));
}

@ Evaluation of an overloaded function call bound to a closure is a simplified
version the part of |call_expression::evaluate| dedicated to a closures,
including the temporary catching of errors in order to produce a trace of
interrupted function calls in the error message. The simplification consists
of the fact that the closure is already evaluated and stored, and that in
particular we don't have to distinguish dynamically between built-in functions
and closures, nor between calls of anonymous or named functions (we are always
in the latter case) for producing the error trace.

@< Function definitions @>=
void overloaded_closure_call::evaluate(level l) const
{ argument->eval();
  try
  { std::vector<shared_value> new_frame;
    new_frame.reserve(count_identifiers(*fun->param));
    thread_components(*fun->param,pop_value(),new_frame);
@)
    context_ptr saved_context(execution_context);
    execution_context.reset(new context(fun->cont,new_frame));
    fun->body->evaluate(l); // pass evaluation level |l| to function body
    execution_context = saved_context;
  }
  catch (const std::exception& e)
  { const std::logic_error* l_err= dynamic_cast<const std::logic_error*>(&e);
    if (l_err!=NULL)
      throw std::logic_error
        (std::string(e.what())+"\n(in call of "+print_name+')');
    throw std::runtime_error
      (std::string(e.what())+"\n(in call of "+print_name+')');
  }
}

@* User-defined functions.
%
Now we shall consider the general case of a user-defined function. In fact all
that needs to be done is type-check and convert the case |lambda_expr| of an
|expr| constructed by the parser; the necessary types derived from
|expression| that provide their implementation were already introduced.

We first test if the required |type| specialises to a function type, i.e.,
either it was some function type or undefined. Then we get the argument type
|arg_type| from the function expression the parser provided; we need to
statically cast from a void pointer that was used to hide from the parser the
class |type_expr| that a \Cee-compiler does not understand. We further
specialise the argument type of |type| to the argument type of the function
(signalling a type error in the rare cases that a different type was
expected). Then a call to |thread_bindings| extracts from the specified
pattern |id_pat| the identifiers it contains, and couples them to the
associated types extracted from |arg_type|, storing the result in
|new_bindings|. These bindings are activated during the type-check and
conversion of the function body, and if all this went well, we check that the

@< Other cases for type-checking and converting... @>=
case lambda_expr:
{ if (not type.specialise(gen_func_type))
    throw type_error(e,copy(gen_func_type),copy(type));
  lambda fun=e.e.lambda_variant;
  id_pat& pat=fun->pattern;
  type_expr& arg_type=*static_cast<type_p>(fun->arg_type);
  if (not type.func->arg_type.specialise(arg_type))
  @/throw type_error(e,
                     make_function_type(copy(arg_type),copy(unknown_type)),
                     copy(type));
  size_t n_id=count_identifiers(pat);
@/bindings new_bindings(n_id);
  thread_bindings(pat,arg_type,new_bindings);
  new_bindings.push(id_context);
  expression_ptr body(convert_expr(fun->body,type.func->result_type));
  new_bindings.pop(id_context);
  return new lambda_expression(pat,body);
}

@* Control structures.
%
We shall now introduce conventional control structures, which must of course
be part of any serious programming language; yet they were implemented only
after plenty of other language elements were in place, such as
let-expressions, functions, rows and selection form them, implicit
conversions.

@*1 Conditional expressions.
%
A first control structure it the conditional expression.

@< Type def... @>=
struct conditional_expression : public expression_base
{ expression condition, then_branch, else_branch;
@)
  conditional_expression(expression_ptr c,expression_ptr t, expression_ptr e)
   : condition(c.release()),then_branch(t.release()), else_branch(e.release())
  @+{}
  virtual ~conditional_expression()
  {@; delete condition; delete then_branch; delete else_branch; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ To print a conditional expression, we reconstruct \&{elif} constructions
that were eliminated in the parser (and even those that the user did not
employ, but could have).

@< Function definitions @>=
void conditional_expression::print(std::ostream& out) const
{ out << " if "; const conditional_expression* cur=this;
  do
  { out << *cur->condition << " then " << *cur->then_branch;
    conditional_expression* p =
      dynamic_cast<conditional_expression*>(cur->else_branch);
    if (p==NULL)
      break;
    out << " elif "; cur=p;
  }
  while(true); // \Cpp\ does not allow using |p| in final condition
  out << " else " << *cur->else_branch << " fi ";
}

@ For type-checking conditional expressions we are in a somewhat similar
situation as for list displays: both branches need to be of the same type, but
we might not know which. After checking that the |condition| yields a Boolean
value, we first convert the else-branch and then the then-branch. This unusual
order is explained by the possibility of an absent else-branch, which the
parser will replace by an empty tuple. If this happens in a context that does
not impose a result type, testing the else-branch first will set |type| to
|void_type|, after which the voiding coercion is available when type-checking
the then-branch; the opposite order might result in an error in converting the
void else-branch to the type of the then-branch.

This is the first place where we use that |bool_type| in not |const|; if it
had been, a possible solution would be to take a copy of it, or to start with
an undefined type and test afterwards that it has become equal to |bool_type|.
The latter option would exclude any coercions to \&{bool} in the condition; at
the time of writing no such coercions exist

@< Other cases for type-checking and converting... @>=
case conditional_expr:
{ expression_ptr c (convert_expr(e.e.if_variant->condition,bool_type));
  expression_ptr el (convert_expr(e.e.if_variant->else_branch,type));
  expression_ptr th (convert_expr(e.e.if_variant->then_branch,type));
  return new conditional_expression(c,th,el);
}

@ Evaluating a conditional expression ends up evaluating either the
then-branch or the else-branch.

@< Function definitions @>=
void conditional_expression::evaluate(level l) const
{ condition->eval();
  if (get<bool_value>()->val)
   then_branch->evaluate(l);
  else else_branch->evaluate(l);
}

@*1 While loops.
%
Next we consider |while| loops, which have three parts (the final one is
optional; if absent it will be a null pointer).

@< Type def... @>=
struct while_expression : public expression_base
{ expression condition, body;
@)
  while_expression(expression_ptr c,expression_ptr b)
   : condition(c.release()),body(b.release())
  @+{}
  virtual ~while_expression() @+ {@; delete condition; delete body; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ Printing a |while| expression is straightforward.

@< Function definitions @>=
void while_expression::print(std::ostream& out) const
{ out << " while " <<  *condition << " do " << *body << " od ";
}

@ Type checking is a bit more complicated for |while| loops that for
conditional expressions, because a row result must be produced. If the context
requires void type, we shall require the same for the body, knowing that
generation of a row value will be unconditionally suppressed in these cases
anyway. In all other cases we proceed for the body expression as for the
components of a row display (except that there is only one expression in this
case).

@< Other cases for type-checking and converting... @>=
case while_expr:
{ w_loop w=e.e.while_variant;
  expression_ptr c (convert_expr(w->condition,bool_type));
  if (type==void_type or type.specialise(row_of_type))
  { expression_ptr b
     (convert_expr(w->body, @|
                   type==void_type ? void_type :*type.component_type));
    @/return new while_expression(c,b);
  }
  else
  @< If |type| can be converted from some row-of type, check |w->body|
     against its component type, construct the |while_expression|, and apply
     the appropriate conversion function to it; otherwise |throw| a
     |type_error| @>
}

@ For |while| loops we follow the same logic for finding an appropriate
component type as for list displays.

@< If |type| can be converted from some row-of type, check |w->body| against
   its component type, construct the |while_expression|, and apply the
   appropriate conversion function to it; otherwise |throw| a |type_error| @>=
{ type_expr comp_type;
  conversion_record* conv = row_coercion(type,comp_type);
  if (conv==NULL)
    throw type_error(e,copy(row_of_type),copy(type));
@)
  expression_ptr b(convert_expr(w->body,comp_type));
  expression_ptr loop(new while_expression(c,b));
  return new conversion(*conv,loop);
}


@ Of course evaluating is what most distinguishes loops from conditionals.

@< Function definitions @>=
void while_expression::evaluate(level l) const
{ if (l==no_value)
    while (condition->eval(),get<bool_value>()->val)
       body->void_eval();
  else
  { row_ptr result (new row_value(0));
    while (condition->eval(),get<bool_value>()->val)
    @/{@; body->eval();
      result->val.push_back(pop_value());
    }
    push_value(result);
  }
}

@*1 For loops.
%
Next we consider |for| loops, which also have three parts, an identifier
pattern defining the loop variable(s), an in-part giving the object whose
components are iterated over, and a body that may produce a new value for each
component of the in-part. We allow iteration over vectors and matrices,
non-row types which are indexable by integers (so in the matrix case iteration
is over its columns), in which case the implicit subscription of the in-part
must be handled appropriately; therefore

@< Type def... @>=
struct for_expression : public expression_base
{ id_pat pattern; expression in_part, body; subscr_base::sub_type kind;
@)
  for_expression
   (id_pat p, expression_ptr i, expression_ptr b, subscr_base::sub_type k);
  virtual ~for_expression() @+
  {@; destroy_id_pat(&pattern); delete in_part; delete body; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ We needed to lift this out of the header file so that the local function
|copy_id_pat| could be used.
@< Function definitions @>=
inline
for_expression::for_expression@|
 (id_pat p, expression_ptr i, expression_ptr b,subscr_base::sub_type k)
   : pattern(copy_id_pat(p)), in_part(i.release()), body(b.release()), kind(k)
  @+{}


@ Printing a |for| expression is straightforward, taking care that
|next_part| could be null.

@< Function definitions @>=
void for_expression::print(std::ostream& out) const
{ out << " for " << pattern.sublist->next->body;
    if (pattern.sublist->body.kind==0x1)
      out << '@@' << pattern.sublist->body;
    out << " in " << *in_part << " do " << *body << " od ";
}

@ Type checking is more complicated for |for| loops than for |while| loops,
since more types and potential coercions are involved. we start by processing
the in-part in a neutral type context, which will on success set |in_type| to
its a priori type. This type must be indexable by integers (so it is either a
row-type or vector or matrix), and the call to |subscr_base::indexable| will
set |comp_type| to the component type resulting from integer subscription.

@< Other cases for type-checking and converting... @>=
case for_expr:
{ f_loop f=e.e.for_variant; type_expr in_type;
  expression_ptr in_expr(convert_expr(f->in_part,in_type));
@/type_expr comp_type; subscr_base::sub_type which;
  if (not subscr_base::indexable(in_type,int_type,comp_type,which))
  { std::ostringstream o;
    o << "Cannot iterate over value of type " << in_type;
    throw expr_error(e,o.str());
  }
  type_ptr pt = pattern_type(f->id);
  type_ptr it_type=make_tuple_type(make_type_list@|
    (copy(int_type),make_type_singleton(copy(comp_type))));
  if (not pt->specialise(*it_type))
    throw expr_error(e,"Improper structure of loop variable pattern");
  size_t n_id = count_identifiers(f->id);
  bindings bind(n_id); thread_bindings(f->id,*it_type,bind);
  type_expr body_type, *btp; conversion_record* conv=NULL;
  if (type==void_type)
    btp=&void_type;
  else if (type.specialise(row_of_type))
    btp=type.component_type;
  else if ((conv=row_coercion(type,body_type))!=NULL)
    btp=&body_type;
  else throw type_error(e,copy(row_of_type),copy(type));
  bind.push(id_context);
  expression_ptr body(convert_expr (f->body,*btp));
  bind.pop(id_context);
  expression_ptr loop(new for_expression(f->id,in_expr,body,which));
  return conv==NULL ? loop.release() : new conversion(*conv,loop);
}

@ We can start evaluating the |in_part| regardless of |kind|, but for deducing
the number of iterations we must already distinguish on |kind| to predict the
type of the in-part. A |loop_frame| is constructed for the new variable(s),
but the shared pointers it contains must be different for each iteration,
because the body might get hold, through a closure that incorporates the
|context| constructed below, of a copy of those pointers, and this closure
should not get to see subsequent values of variables. On the other hand, the
syntax does not allow users to name the entire tuple |loop_var| formed of the
loop index and the in-part component, so this pointer cannot end up in the
|loop_frame| and it may remain the same between iterations.

@< Function definitions @>=
void for_expression::evaluate(level l) const
{ size_t n_id = count_identifiers(pattern);
  context_ptr saved_context=execution_context;
  in_part->eval();
  shared_tuple loop_var(new tuple_value(2));
       // this is safe to re-use between iterations
  std::vector<shared_value> loop_frame; loop_frame.reserve(n_id);

  row_ptr result(NULL);
  @< Evaluate the loop, dispatching the various possibilities for |kind|, and
  setting |result| @>

  execution_context = saved_context;
  if (l!=no_value)
    push_value(result);
}

@ For evaluating |for| loops we must take care to interpret the |kind| field
when selecting a component from the in-part. Because of differences in the
type of |in_val|, some code must be duplicated, which we do as much as
possible by sharing a module between the various loop bodies.

@< Evaluate the loop, dispatching the various possibilities for |kind|... @>=
switch (kind)
{ case subscr_base::row_entry:
  { shared_row in_val = get<row_value>();
    size_t n=in_val->val.size();
    if (l!=no_value)
      result = row_ptr(new row_value(n));
    for (size_t i=0; unsigned(i)<n; ++i,loop_frame.clear())
    { loop_var->val[1]=in_val->val[i];
      @< Set |loop_var->val[0]| to |i|, fill |loop_frame| according to
      |pattern| with values from |loop_var|, create a new |context| and
      evaluate the |loop_body| in it, and maybe assign |result->val[i]|
      from it @>
    }
  }
  break;
  case subscr_base::vector_entry:
  { shared_vector in_val = get<vector_value>();
    size_t n=in_val->val.size();
    if (l!=no_value)
      result = row_ptr(new row_value(n));
    for (size_t i=0; unsigned(i)<n; ++i,loop_frame.clear())
    { loop_var->val[1].reset(new int_value(in_val->val[i]));
      @< Set |loop_var->val[0]| to |i|,... @>
    }
  }
  break;
  case subscr_base::ratvec_entry:
  { shared_rational_vector in_val = get<rational_vector_value>();
    size_t n=in_val->val.size();
    if (l!=no_value)
      result = row_ptr(new row_value(n));
    for (size_t i=0; unsigned(i)<n; ++i,loop_frame.clear())
    { loop_var->val[1].reset(new rat_value(arithmetic::Rational @|
        (in_val->val.numerator()[i],in_val->val.denominator())));
      @< Set |loop_var->val[0]| to |i|,... @>
    }
  }
  break;
  case subscr_base::matrix_column:
  { shared_matrix in_val = get<matrix_value>();
    size_t n=in_val->val.numColumns();
    if (l!=no_value)
      result = row_ptr(new row_value(n));
    for (size_t i=0; unsigned(i)<n; ++i,loop_frame.clear())
    { loop_var->val[1].reset(new vector_value(in_val->val.column(i)));
      @< Set |loop_var->val[0]| to |i|,... @>
    }
  }
  break;
  case subscr_base::matrix_entry: break; // excluded in type analysis
}


@ Since the module below exists only for the sake of source-code sharing, we
don't bother to put braces around its expansion, as they are not needed in the
uses above.

We set the in-part component stored in |loop_var->val[1]| separately for the
various values of |kind|, but |loop_var->val[0]| is always the (integral) loop
index. Once initialised, |loop_var| is passed through the function
|thread_components| to set up |loop_frame|, whose pointers are copied into a
new |context| that extends the initial |saved_context| to form the new
|execution_context|. Like for |loop_var->val[0]|, it is important that
|execution_context| be set to point to a newly created node at each iteration,
since any closure values in the loop body will incorporate its current value;
there would be no point in supplying fresh pointers in |loop_var| if they were
subsequently copied to overwrite the pointers in the same |context| object
each time. Once these things have been handled, the evaluation of the loop
body is standard.

@< Set |loop_var->val[0]| to |i|,... @>=
loop_var->val[0].reset(new int_value(i)); // must be newly created each time
thread_components(pattern,loop_var,loop_frame);
execution_context.reset(new context(saved_context,loop_frame));
  // this one too
if (l==no_value)
  body->void_eval();
else
{@; body->eval(); result->val[i]=pop_value(); }

@*1 Counted loops.
%
Next we consider counted |for| loops. Increasing and decreasing loops give
distinct types.

@< Type def... @>=
struct inc_for_expression : public expression_base
{ expression count, bound, body; Hash_table::id_type id;
@)
  inc_for_expression@/
   (Hash_table::id_type i, expression_ptr cnt, expression_ptr bnd,
    expression_ptr b)
  : count(cnt.release()),bound(bnd.release()),body(b.release()),id(i)
  @+{}
  virtual ~inc_for_expression() @+
  {@; delete count; delete bound; delete body; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

struct dec_for_expression : public expression_base
{ expression count, bound, body; Hash_table::id_type id;
@)
  dec_for_expression@/
   (Hash_table::id_type i, expression_ptr cnt, expression_ptr bnd,
    expression_ptr b)
  : count(cnt.release()),bound(bnd.release()),body(b.release()),id(i)
  @+{}
  virtual ~dec_for_expression() @+
  {@; delete count; delete bound; delete body; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ Printing a counted |for| expression is straightforward; we don't bother to
suppress ``\&{from}~0''.

@< Function definitions @>=
void inc_for_expression::print(std::ostream& out) const
{ out << " for " << main_hash_table->name_of(id)  << " = " << *count @|
      << " from " << *bound << " do " << *body << " od ";
}
void dec_for_expression::print(std::ostream& out) const
{ out << " for " << main_hash_table->name_of(id)  << " = " << *count @|
      << " downto " << *bound << " do " << *body << " od ";
}

@ Type-checking counted |for| loops is rather like that of |while| loops, but
we must extend the context with the loop variable while processing the loop
body.

@< Other cases for type-checking and converting... @>=
case cfor_expr:
{ c_loop c=e.e.cfor_variant;
  expression_ptr count_expr(convert_expr(c->count,int_type));
  static shared_value zero=shared_value(new int_value(0));
    // avoid repeated allocation
  expression_ptr bound_expr
    (is_empty(c->bound) ? new denotation(zero)
    : convert_expr(c->bound,int_type)
    );
@)
  bindings bind(1); bind.add(c->id,copy(int_type));
  type_expr body_type, *btp; conversion_record* conv=NULL;
  if (type==void_type)
    btp=&void_type;
  else if (type.specialise(row_of_type))
    btp=type.component_type;
  else if ((conv=row_coercion(type,body_type))!=NULL)
    btp=&body_type;
  else throw type_error(e,copy(row_of_type),copy(type));
  bind.push(id_context);
  expression_ptr body(convert_expr (c->body,*btp));
  bind.pop(id_context);
  expression e;
  if (c->up!=0)
    e=new inc_for_expression(c->id,count_expr,bound_expr,body);
  else
    e=new dec_for_expression(c->id,count_expr,bound_expr,body);
  expression_ptr loop(e);
  return conv==NULL ? loop.release() : new conversion(*conv,loop);

}

@ Executing a loop is a simple variation of what we have seen before for
|while| and |for| loops.

@< Function definitions @>=
void inc_for_expression::evaluate(level l) const
{ int b=(bound->eval(),get<int_value>()->val);
  int c=(count->eval(),get<int_value>()->val);
  if (c<0)
    c=0; // no negative size result

  context_ptr saved_context=execution_context;
  std::vector<shared_value>loop_frame(1);
  shared_value& loop_var=loop_frame[0];
  if (l==no_value)
  { c+=b;
    for (int i=b; i<c; ++i)
    { loop_var.reset(new int_value(i));
      execution_context.reset(new context(saved_context,loop_frame));
      body->void_eval();
    }
  }
  else
  { row_ptr result (new row_value(0)); result->val.reserve(c);
    c+=b;
    for (int i=b; i<c; ++i)
    { loop_var.reset(new int_value(i));
      execution_context.reset(new context(saved_context,loop_frame));
      body->eval(); result->val.push_back(pop_value());
    }
    push_value(result);
  }
  execution_context=saved_context;
}

@ Downward loops are not much different, but they actually use a |while| loop.

@< Function definitions @>=
void dec_for_expression::evaluate(level l) const
{ int b=(bound->eval(),get<int_value>()->val);
  int i=(count->eval(),get<int_value>()->val);
  if (i<0)
    i=0; // no negative size result

  context_ptr saved_context=execution_context;
  std::vector<shared_value>loop_frame(1);
  shared_value& loop_var=loop_frame[0];
  if (l==no_value)
  { i+=b;
    while (i-->b)
    { loop_var.reset(new int_value(i));
      execution_context.reset(new context(saved_context,loop_frame));
      body->void_eval();
    }
  }
  else
  { row_ptr result (new row_value(0)); result->val.reserve(i);
    i+=b;
    while (i-->b)
    { loop_var.reset(new int_value(i));
      execution_context.reset(new context(saved_context,loop_frame));
      body->eval(); result->val.push_back(pop_value());
    }
    push_value(result);
  }
  execution_context=saved_context;
}

@* Casts.
%
Casts are very simple to process; they do not need any |expression| type to
represent them.

@< Other cases for type-checking and converting... @>=
case cast_expr:
{ cast c=e.e.cast_variant;
  type_expr& ctype=*static_cast<type_p>(c->type);
  expression_ptr p(convert_expr(c->exp,ctype));
  return conform_types(ctype,type,p,e);
}

@ The overload table stores type information in a |func_type| value, which
cannot be handed directly to the |specialise| method. The following function
simulates specialisation to a function type |from|$\to$|to|.

@< Local function definitions @>=
inline bool spec_func(type_expr& t, const type_expr& from, const type_expr& to)
{ return t.specialise(gen_func_type) @|
  and t.func->arg_type.specialise(from) @|
  and t.func->result_type.specialise(to);
}

@ Operation casts similarly only access existing kinds of expression. We must
however access the global overload table to find the value. Since upon success
we find a bare function value, we must abuse the |denotation| class a bit to
serve as wrapper that upon evaluation will return the value again.

@< Other cases for type-checking and converting... @>=
case op_cast_expr:
{ op_cast c=e.e.op_cast_variant;
  const overload_table::variant_list& variants =
   global_overload_table->variants(c->oper);
  type_expr& ctype=*static_cast<type_p>(c->type);
  if (is_special_operator(c->oper))
    @< Test special argument patterns, and on match |return| an appropriate
       denotation @>
  for (size_t i=0; i<variants.size(); ++i)
    if (variants[i].type->arg_type==ctype)
    {
      expression_ptr p(new denotation(variants[i].val));
      if (spec_func(type,ctype,variants[i].type->result_type))
        return p.release();
      type_ptr ftype=make_function_type
	(copy(ctype),copy(variants[i].type->result_type));
      throw type_error(e,ftype,copy(type));
    }
  std::ostringstream o;
  o << "Cannot resolve " << main_hash_table->name_of(c->oper) @|
    << " at argument type " << ctype;
  throw program_error(o.str());
}
break;

@ For now size-of is the only identifier with special operand patterns.

@< Test special argument patterns... @>=
{ if (c->oper==size_of_name())
  { if (ctype.kind==row_type)
    { if (spec_func(type,ctype,int_type))
      return new denotation(shared_value
        (new builtin_value(sizeof_wrapper,"#@@[T]")));
      throw type_error(e,copy(ctype),copy(type));
    }
    else if (ctype.specialise(pair_type))
    { type_expr& arg_tp0 = ctype.tuple->t;
      type_expr& arg_tp1 = ctype.tuple->next->t;
      if (arg_tp0.kind==row_type)
      { if (arg_tp0==arg_tp1)
        { if (spec_func(type,ctype,arg_tp0))
          return new denotation(shared_value @|
            (new builtin_value(join_rows_wrapper,"#@@([T],[T]->[T])")));
          throw type_error(e,copy(ctype),copy(type));
        }
	else if (*arg_tp0.component_type==arg_tp1)
        { if (spec_func(type,ctype,arg_tp0))
          return new denotation(shared_value @|
            (new builtin_value(suffix_element_wrapper,"#@@([T],T->[T])")));
          throw type_error(e,copy(ctype),copy(type));
        }
      }
      if (arg_tp1.kind==row_type and *arg_tp1.component_type==arg_tp0)
      { if (spec_func(type,ctype,arg_tp1))
        return new denotation(shared_value @|
          (new builtin_value(prefix_element_wrapper,"#@@(T,[T]->[T])")));
        throw type_error(e,copy(ctype),copy(type));
      }
    }
  }
}

@* Assignments.
%
Syntactically there is hardly anything simpler than simple assignment
statements. However, semantically we distinguish assignments to local and to
global variables; the two kinds will be converted into objects of classes
derived from |assignment_expr|.

@< Type definitions @>=
struct assignment_expr : public expression_base
{ Hash_table::id_type lhs;
  expression rhs;
@)
  assignment_expr(Hash_table::id_type l,expression_ptr r)
   : lhs(l),rhs(r.release()) @+{}
  virtual ~assignment_expr() @+{@; delete rhs; }
  virtual void print(std::ostream& out) const;
};

@ Printing proceeds as usual for identifiers and converted expressions.
@< Function def... @>=
void assignment_expr::print(std::ostream& out) const
{@; out << main_hash_table->name_of(lhs) << ":=" << *rhs; }

@ For global assignments we need to access the location where the identifier
is stored, as in the case of applied identifiers, but we need a non-|const|
access to it.

@< Type definitions @>=
class global_assignment : public assignment_expr
{ shared_share address;
public:
  global_assignment(Hash_table::id_type l,expression_ptr r);
  virtual ~global_assignment() @+{}
  virtual void evaluate(level l) const;
};

@ The constructor for |global_assignment| is rather similar to that for
|global_identifier|.

@< Function def... @>=
global_assignment::global_assignment(Hash_table::id_type l,expression_ptr r)
: assignment_expr(l,r), address(global_id_table->address_of(l)) @+{}

@ Evaluating a global assignment evaluates the left hand side, and replaces
the old value stored at |*address| by the new (shared pointer) value.

@< Function def... @>=
void global_assignment::evaluate(level l) const
{@; rhs->eval();
  *address = pop_value();
  push_expanded(l,*address);
}

@ For local assignments we also need to access the location where the
identifier is stored, but as in the case of applied identifiers this location
cannot be known when the statement is type-checked and converted (and in fact
it may vary between evaluations of the same assignment) so we need to store
coordinates of the identifier in the evaluation context.

@< Type definitions @>=
class local_assignment : public assignment_expr
{ size_t depth, offset;
public:
  local_assignment(Hash_table::id_type l, size_t i,size_t j, expression_ptr r);
  virtual ~local_assignment() @+{}
  virtual void evaluate(level l) const;
};

@ The constructor for |local_assignment| is straightforward.

@< Function def... @>=
local_assignment::local_assignment
 (Hash_table::id_type l, size_t i,size_t j, expression_ptr r)
: assignment_expr(l,r), depth(i), offset(j) @+{}

@ Evaluating a local assignment evaluates the left hand side, and replaces the
old value stored at |execution_context->elem(depth,offset)| by the new (shared
pointer) value.

@< Function def... @>=
void local_assignment::evaluate(level l) const
{ rhs->eval();
  shared_value& dest =  execution_context->elem(depth,offset);
  dest= pop_value();
  push_expanded(l,dest);
}


@ Type-checking and converting assignment statements follows the same lines as
that of applied identifiers as far as discriminating between local and global
is concerned. We first look in |id_context| for a local binding of the
identifier, and if not found we look in |global_id_table|. If found in either
way,

@< Other cases for type-checking and converting... @>=
case ass_stat:
{ Hash_table::id_type lhs=e.e.assign_variant->lhs;
  const expr& rhs=e.e.assign_variant->rhs;
  type_p it; expression_ptr assign; size_t i,j;
  if ((it=id_context->lookup(lhs,i,j))!=NULL)
  @/{@; expression_ptr r(convert_expr(rhs,*it));
    assign.reset(new local_assignment(lhs,i,j,r));
  }
  else if ((it=global_id_table->type_of(lhs))!=NULL)
  @/{@; expression_ptr r(convert_expr(rhs,*it));
    assign.reset(new global_assignment(lhs,r));
  }
  else throw program_error @|
    (std::string("Undefined identifier in assignment: ")
     +main_hash_table->name_of(lhs));
@.Undefined identifier in assignment@>
  return conform_types(*it,type,assign,e);
}

@*1 Component assignments.
%
The language we are implementing does not employ the notion of sub-object; in
other words if one sets $b=a[i]$ for some list, vector or matrix $a$, then $s$
will behave as a copy of the entry $a[i]$ rather than as an alias, so
subsequent assignment to $b$ will not affect~$a$ or vice versa. (This does no
prevent us to share storage between $b$ and $a$ initially, it just means the
sharing should be broken if $b$ or $a$ are modified.) This simplifies the
semantic model considerably, but if we want to allow creating composite values
by subsequently setting their components, we need to allow assignments of the
form $a[i]:=c$. The meaning of this is the same as assigning a new value to
all of $a$ that differs from the original value only at index~$i$; it may
however be expected to be implemented more efficiently if the storage of $a$
is not currently shared, as would usually be the case at least from the second
such assignment to~$a$ on. The interpreter will have to treat such component
assignments as a whole (with three components $a,i,c$), which also means that
it will not be able to handle something like $a[i][j]:=c$ even when that would
seem to make sense (however $m[i,j]:=c$ for matrix values $m$ will be
supported).

In fact we need to implement a whole range of component assignments: there are
assignments to general row-value components, to vector and matrix components
and to matrix columns, and all this for local variables as well as for global
ones. Like for general assignments we can start with a base class that
implements common methods, and derive the specialised versions from it later.
In fact by deriving from |assignment_expr| we only need to add the index as
data member. We also provide a method |assign| that will do the real work for
the |evaluate| methods of the derived classes, after those have located
address of the aggregate to be modified and the type of component assignment
to apply.

@< Type definitions @>=
struct component_assignment : public assignment_expr
{ expression index;
  component_assignment
   (Hash_table::id_type a,expression_ptr i,expression_ptr r)
   : assignment_expr(a,r), index(i.release()) @+{}

  virtual ~component_assignment() @+{@; delete index; }
  virtual void print (std::ostream& out) const;
@)
  void assign(level l,shared_value& aggregate,subscr_base::sub_type kind) const;
};

@ Printing reassembles the subexpressions according to the input syntax.
@< Function def...@>=
void component_assignment::print(std::ostream& out) const
{@; out << main_hash_table->name_of(lhs) << '[' << *index << "]:=" << *rhs; }

@ For global assignments we need to non-|const| access the location where the
identifier is stored.

@< Type definitions @>=
class global_component_assignment : public component_assignment
{ subscr_base::sub_type kind;
  shared_share address;
public:
  global_component_assignment
    (Hash_table::id_type a,expression_ptr i,expression_ptr r,
     subscr_base::sub_type k);
  virtual void evaluate(level l) const;
};

@ The constructor for |global_component_assignment| stores the address of the
aggregate object and the component kind.

@< Function def... @>=
global_component_assignment::global_component_assignment @|
  (Hash_table::id_type a,expression_ptr i,expression_ptr r,
   subscr_base::sub_type k)
: component_assignment(a,i,r)
, kind(k),address(global_id_table->address_of(a)) @+{}

@ It is in evaluation that component assignments differ most from ordinary
ones. The work is delegated to the |assign| method of the base class, which is
given a reference to the |shared_value| pointer holding the current value of
the aggregate; it is this pointer that is in principle modified. Like when
fetching the value of a global variable, we must be aware of a possible
undefined value in the variable.

@< Function def... @>=
void global_component_assignment::evaluate(level l) const
{ if (address->get()==NULL)
  { std::ostringstream o;
    o << "Assigning to component of uninitialized variable "
      << main_hash_table->name_of(lhs);
    throw std::runtime_error(o.str());
  }
  assign(l,*address,kind);
}

@ The |assign| method, which will also be called for local component
assignments, starts by the common work of evaluating the index and the value
to be assigned, and of making sure the aggregate variable is made to point to
a unique copy of its current value, which copy can then be modified in place.
For actually changing the aggregate, we must distinguish cases according to
the kind o component assignment at hand.

@< Function def... @>=
void component_assignment::assign
  (level l,shared_value& aggregate, subscr_base::sub_type kind) const
{ rhs->eval();
  uniquify(aggregate);
  value loc=aggregate.get(); // simple reference from shared pointer
  switch (kind)
  { case subscr_base::row_entry:
  @/@< Replace component at |index| in row |loc| by value on stack @>
  @+break;
    case subscr_base::vector_entry:
  @/@< Replace entry at |index| in vector |loc| by value on stack @>
  @+break;
    case subscr_base::matrix_entry:
  @/@< Replace entry at |index| in matrix |loc| by value on stack @>
  @+break;
    case subscr_base::matrix_column:
  @/@< Replace columns at |index| in matrix |loc| by value on stack @>
  @+break;
  case subscr_base::ratvec_entry: {} // case is eliminated in type analysis
  }
}

@ A |row_value| component assignment is the simplest kind. The variable |loc|
holds a generic pointer, known to refer to a |row_value|. Since we need their
value, but not for incorporation into a result, we use |force| to get ordinary
pointers of the right kind. Then we do a bound check, and on success replace a
component of the value held in |a| by the stack-top value. Afterwards,
depending on |l|, we may put back the stack-top value as result of the
component assignment, possibly expanding a tuple in the process.

@< Replace component at |index| in row |loc|... @>=
{ unsigned int i=(index->eval(),get<int_value>()->val);
  row_value& a=*force<row_value>(loc);
  if (i>=a.val.size())
    throw std::runtime_error(range_mess(i,a.val.size(),this));
  a.val[i]= pop_value();
  push_expanded(l,a.val[i]);
}

@ For |vec_value| entry assignments the type of the aggregate object is
vector, and the value assigned always an integer. The latter certainly needs
no expansion, so we either leave it on the stack, or remove it if the value of
the component assignment expression is not used.

@< Replace entry at |index| in vector |loc|... @>=
{ unsigned int i=(index->eval(),get<int_value>()->val);
  vector_value& v=*force<vector_value>(loc);
  if (i>=v.val.size())
    throw std::runtime_error(range_mess(i,v.val.size(),this));
  v.val[i]= force<int_value>(execution_stack.back().get())->val;
  if (l==no_value)
    execution_stack.pop_back();
}

@ For |mat_value| entry assignments at |index| must be split into a pair of
indices, and there are two bound checks.

@< Replace entry at |index| in matrix |loc|... @>=
{ index->multi_eval();
  unsigned int j=get<int_value>()->val;
  unsigned int i=get<int_value>()->val;
@/
  matrix_value& m=*force<matrix_value>(loc);
  if (i>=m.val.numRows())
    throw std::runtime_error(range_mess(i,m.val.numRows(),this));
  if (j>=m.val.numColumns())
    throw std::runtime_error(range_mess(j,m.val.numColumns(),this));
  m.val(i,j)= force<int_value>(execution_stack.back().get())->val;
  if (l==no_value)
    execution_stack.pop_back();
}

@ A |matrix_value| column assignment is like that of a vector entry, but we
add a test for matching column length.

@< Replace columns at |index| in matrix |loc|... @>=
{ unsigned int j=(index->eval(),get<int_value>()->val);
  matrix_value& m=*force<matrix_value>(loc);
  const vector_value& v=*force<vector_value>(execution_stack.back().get());
  if (j>=m.val.numColumns())
    throw std::runtime_error(range_mess(j,m.val.numColumns(),this));
  if (v.val.size()!=m.val.numRows())
    throw std::runtime_error
      (std::string("Cannot replace column of size ")+str(m.val.numRows())+
       " by one of size "+str(v.val.size()));
  m.val.set_column(j,v.val);
  if (l==no_value)
    execution_stack.pop_back();
}

@ For local assignments we also need to access the location where the
identifier is stored, which as before is done by storing coordinates of the
identifier in the evaluation context.

@< Type definitions @>=
class local_component_assignment : public component_assignment
{ subscr_base::sub_type kind;
  size_t depth, offset;
public:
  local_component_assignment @|
   (Hash_table::id_type l, expression_ptr i,size_t d, size_t o,
    expression_ptr r, subscr_base::sub_type k);
  virtual void evaluate(level l) const;
};

@ The constructor for |local_component_assignment| is straightforward, in
spite of the number of arguments.

@< Function def... @>=
local_component_assignment::local_component_assignment
 (Hash_table::id_type l, expression_ptr i,size_t d, size_t o, expression_ptr r,
  subscr_base::sub_type k)
: component_assignment(l,i,r), kind(k), depth(d), offset(o) @+{}

@ The |evaluate| method locates the |shared_value| pointer of the aggregate,
calls |assign| to do the work.

@< Function def... @>=
void local_component_assignment::evaluate(level l) const
{@; assign(l,execution_context->elem(depth,offset),kind); }

@ Type-checking and converting component assignment statements follows the
same lines as that of ordinary assignment statements, but must also
distinguish different aggregate types.

@< Other cases for type-checking and converting... @>=
case comp_ass_stat:
{ Hash_table::id_type aggr=e.e.comp_assign_variant->aggr;
  const expr& index=e.e.comp_assign_variant->index;
  const expr& rhs=e.e.comp_assign_variant->rhs;
@/type_p aggr_t; type_expr ind_t; type_expr comp_t;
  expression_ptr assign; size_t d,o; bool is_local;
  if ((aggr_t=id_context->lookup(aggr,d,o))!=NULL)
    is_local=true;
  else if ((aggr_t=global_id_table->type_of(aggr))!=NULL)
    is_local=false;
  else throw program_error @|
    (std::string("Undefined identifier in component assignment: ")
     +main_hash_table->name_of(aggr));
@.Undefined identifier in assignment@>
@)
  expression_ptr i(convert_expr(index,ind_t));
  subscr_base::sub_type kind;
  if (not subscr_base::indexable(*aggr_t,ind_t,comp_t,kind)
      or kind==subscr_base::ratvec_entry)
  { std::ostringstream o;
    o << "Cannot subscript " << *aggr_t << @| " value with index of type "
      << ind_t << " in assignment";
    throw expr_error(e,o.str());
  }
  expression_ptr r(convert_expr(rhs,comp_t));
  if (is_local)
    assign.reset(new local_component_assignment(aggr,i,d,o,r,kind));
  else
    assign.reset(new global_component_assignment(aggr,i,r,kind));

  return conform_types(comp_t,type,assign,e);
}

@* Sequence expressions.
%
Since sequences are probably short on average, we use a chained
representation after type analysis, just like before. The forward and reverse
variants are implemented by similar but distinct types derived from
|expression_base|.

@< Type def... @>=
struct seq_expression : public expression_base
{ expression first,last;
@)
  seq_expression(expression_ptr f,expression_ptr l)
   : first(f.release()),last(l.release()) @+{}
  virtual ~seq_expression() @+ {@; delete first; delete last; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};
@)
struct next_expression : public expression_base
{ expression first,last;
@)
  next_expression(expression_ptr f,expression_ptr l)
   : first(f.release()),last(l.release()) @+{}
  virtual ~next_expression() @+ {@; delete first; delete last; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ To print a sequence, we print the two expressions separated by a semicolon
or \.{next}.

@< Function definitions @>=
void seq_expression::print(std::ostream& out) const
{@; out << *first << ';' << *last; }
void next_expression::print(std::ostream& out) const
{@; out << *first << " next " << *last; }

@ Evaluating a sequence expression evaluates the |first| for side effects
only, and then the |last| expression.

@< Function def... @>=
void seq_expression::evaluate(level l) const
{@; first->void_eval(); last->evaluate(l); }

@ For a next-expression, it is the value of the first expression that is
retained as result.

@< Function def... @>=
void next_expression::evaluate(level l) const
{@; first->evaluate(l); last->void_eval(); }

@ It remains to type-check and convert sequence expressions, which is easy.

@< Other cases for type-checking and converting... @>=
case seq_expr:
{ sequence seq=e.e.sequence_variant;
  if (seq->forward!=0)
  { expression_ptr first(convert_expr(seq->first,void_type));
    expression_ptr last(convert_expr(seq->last,type));
    return new seq_expression(first,last);
  }
  else
  { expression_ptr first(convert_expr(seq->first,type));
    expression_ptr last(convert_expr(seq->last,void_type));
    return new next_expression(first,last);
  }
}

@* Invoking the type checker.
%
Let us recapitulate what will happen. The parser will read what the user
types, and returns an |expr| value. Then we shall call |convert_expr| for this
value, which will either produce (a pointer to) an executable object of a type
derived from |expression_base|, or throw an exception in case a type error or
other problem is detected. In the former case we are ready to call the
|evaluate| method of the value returned by |convert_expr|; after this main
program will print the result. The call to |convert_expr| is done via
|analyse_types|, takes care of printing error messages for exceptions thrown.

@< Declarations of exported functions @>=
type_ptr analyse_types(const expr& e,expression_ptr& p)
   throw(std::bad_alloc,std::runtime_error);

@~An important reason for having this function is that it will give us an
occasion to catch any thrown |type_error| and |program_error| exceptions,
something we did not want to do inside the recursive function |convert_expr|;
since we cannot return normally from |analyse_types| in the presence of these
errors, we map these errors to |std::runtime_error|, an exception for which
the code that calls us will have to provide a handler anyway.

@< Function definitions @>=
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

@* Primitive types for vectors and matrices.
%
The interpreter distinguishes its own types like \.{[int]} ``row of integer''
from similar built-in types of the library, like \.{vec} ``vector'', which it
will consider to be primitive types. In fact a value of type ``vector''
represents an object of the Atlas type |atlas::matrix::Vector<int>|, and
similarly other primitive types will stand for other Atlas types. We prefer
using a basic type name rather than something resembling the more
mathematically charged equivalents like |latticetypes::Weight| for
|atlas::matrix::Vector<int>|, as that might be more confusing that helpful to
users. In any case, the interpretation of the values is not at all fixed
(vectors are used for coweights and (co)roots as well as for weights, and
matrices could denote either a basis or an automorphism of a lattice).

There is one type that is genuinely defined in the \.{latticetypes} unit,
namely |latticetypes::RatWeight|, and we shall use that, but call the
resulting type just rational vector (\.{ratvec} for the user).

@< Includes needed in the header file @>=
#include "latticetypes.h"

@ We start with deriving |vector_value| from |value_base|. In its constructor,
the argument is a reference to |std::vector<int>|, from which
|matrix::Vector<int>| is derived (without adding data members); since a
constructor for the latter from the former is defined, we can do with just one
constructor for |vector_value|.

@< Type definitions @>=

struct vector_value : public value_base
{ matrix::Vector<int> val;
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

@ Matrices and rational vectors follow the same pattern, but have constructors
from identical types to the one stored.

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
{ latticetypes::RatWeight val;
@)
  explicit rational_vector_value(const latticetypes::RatWeight& v):val(v)@+{}
  rational_vector_value(const matrix::Vector<int>& v,int d)
   : val(v,d) @+ { val.normalize(); }
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

@< Function def... @>=
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

@< Function def... @>=
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
The conversions into vectors or matrices these use an auxiliary function
|row_to_weight|, which constructs a new |Weight| from a row of integers,
leaving the task to clean up that row to their caller.

Note that |row_to_weight| returns its result by value (rather than by
assignment to a reference parameter); in principle this involves making a copy
of the vector |result| (which includes duplicating its entries). However since
there is only one |return| statement, the code generated by a decent compiler
like \.{g++} creates the vector object directly at its destination in this
case (whose location is passed as a hidden pointer argument), and does not
call the copy constructor at all More generally this is avoided if all
|return| statements return the same variable, or if all of them return
expressions instead. So we shall not hesitate to use this idiom henceforth; in
fact, while it was exceptional when this code was first written, it is now
being used throughout the Atlas library.

@< Local function def... @>=
matrix::Vector<int> row_to_weight(const row_value& r)
{ matrix::Vector<int> result(r.val.size());
  for(size_t i=0; i<r.val.size(); ++i)
    result[i]=force<int_value>(r.val[i].get())->val;
  return result;
}
@)
void vector_convert()
{@; shared_row r(get<row_value>());
  push_value(new vector_value(row_to_weight(*r)));
}

@ For |matrix_conversion|, the |evaluate| method is longer, but still
straightforward. Since the type of the argument was checked to be \.{[vec]},
we can safely cast the argument pointer to |(row_value*)|, and each of its
component pointers to |(vector_value*)|. Any ragged columns are silently
extended with null entries to make a rectangular shape for the matrix.

@< Local function def... @>=
void matrix_convert()
{ shared_row r(get<row_value>());
@/std::vector<matrix::Vector<int> > column_list;
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

@ All that remains is to initialise the |coerce_table|.
@< Initialise evaluator @>=
coercion(row_of_int_type, vec_type, "V", vector_convert); @/
coercion(row_of_vec_type,mat_type, "M", matrix_convert);

@ Here are conversions involving rational numbers and vectors.

@< Local function def... @>=
void rational_convert() // convert integer to rational (with denominator~1)
{@; shared_int i = get<int_value>();
    push_value(new rat_value(arithmetic::Rational(i->val)));
}
@)
void ratvec_convert() // convert list of rationals to rational vector
{ shared_row r = get <row_value>();
  matrix::Vector<int> numer(r->val.size()),denom(r->val.size());
  unsigned int d=1;
  for (size_t i=0; i<r->val.size(); ++i)
  { arithmetic::Rational frac = force<rat_value>(r->val[i].get())->val;
    numer[i]=frac.numerator();
    denom[i]=frac.denominator();
    d=arithmetic::lcm(d,denom[i]);
  }
  for (size_t i=0; i<r->val.size(); ++i)
    numer[i]*= d/denom[i]; // adjust numerators to common denominator

  push_value(new rational_vector_value(latticetypes::RatWeight(numer,d)));
}
@)
void rat_list_convert() // convert rational vector to list of rationals
{ shared_rational_vector rv = get<rational_vector_value>();
  row_ptr result(new row_value(rv->val.size()));
  for (size_t i=0; i<rv->val.size(); ++i)
  { arithmetic::Rational q(rv->val.numerator()[i],rv->val.denominator());
    result->val[i] = shared_value(new rat_value(q.normalize()));
  }
  push_value(result);
}


@ There remains one ``internalising'' conversion function, from row of row of
integer to matrix, and the ``externalising'' counterparts (towards lists) of
the vector and matrix conversions. The former is similar to |matrix_convert|.

@< Local function def... @>=
void matrix2_convert()
{ shared_row r(get<row_value>());
@/std::vector<matrix::Vector<int> > column_list;
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

@ For the ``externalising'' conversions, it will be handy to have a basic
function |weight_to_row| that performs more or less the inverse transformation
of |row_to_weight|, but rather than returning a |row_value| it returns a
|row_ptr| pointing to it.

@< Local function def... @>=
row_ptr weight_to_row(const matrix::Vector<int>& v)
{ row_ptr result (new row_value(v.size()));
  for(size_t i=0; i<v.size(); ++i)
    result->val[i]=shared_value(new int_value(v[i]));
  return result;
}
@)
void int_list_convert()
{@; shared_vector v(get<vector_value>());
  push_value(weight_to_row(v->val));
}
@)
void vec_list_convert()
{ shared_matrix m=get<matrix_value>();
  row_ptr result(new row_value(m->val.numColumns()));
  for(size_t i=0; i<m->val.numColumns(); ++i)
    result->val[i]=shared_value(new vector_value(m->val.column(i)));
  push_value(result);
}
@)
void int_list_list_convert()
{ shared_matrix m=get<matrix_value>();
  row_ptr result(new row_value(m->val.numColumns()));
  for(size_t i=0; i<m->val.numColumns(); ++i)
    result->val[i]=shared_value(weight_to_row(m->val.column(i)).release());

  push_value(result);
}

@ All that remains is to initialise the |coerce_table|.
@< Initialise evaluator @>=
coercion(int_type,rat_type, "Q", rational_convert); @/
coercion(row_of_rat_type,ratvec_type, "QV", ratvec_convert); @/
coercion(ratvec_type,row_of_rat_type, "[Q]", rat_list_convert);
@)
coercion(row_row_of_int_type,mat_type, "M2", matrix2_convert); @/
coercion(vec_type,row_of_int_type, "[I]", int_list_convert); @/
coercion(mat_type,row_of_vec_type, "[V]", vec_list_convert); @/
coercion(mat_type,row_row_of_int_type, "[[I]]", int_list_list_convert); @/

@* Wrapper functions.
%
We now come to defining wrapper functions. The following function will greatly
facilitate the latter repetitive task of installing them.

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

@< Function def... @>=
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

@ Our first built-in functions implement with integer arithmetic. Arithmetic
operators are implemented by wrapper functions with two integer arguments.
Since arguments top built-in functions are evaluated with |level| parameter
|multi_value|, two separate value will be produced on the stack. Note that
these are pulled from the stack in reverse order, which is important for the
non-commutative operations like `|-|' and `|/|'. Since values are shared, we
must allocate new value objects for the results.

@h "intutils.h"

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

@ We take the occasion to repair the integer division operation for negative
arguments, by using |intutils::divide| rather than |operator/|. Since that
takes an unsigned second argument, we handle the case |j<0| ourselves.

@< Local function definitions @>=
void divide_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (j==0) throw std::runtime_error("Division by zero");
  if (l!=expression_base::no_value)
    push_value(new int_value
     (j>0 ? intutils::divide(i,j) : -intutils::divide(i,-j)));
}

@ We also define a remainder operation |modulo|, a combined
quotient-and-remainder operation |divmod|, unary subtraction, exact division
of integers producing a rational number, and an integer power operation
(defined whenever the result is integer).

@< Local function definitions @>=
void modulo_wrapper(expression_base::level l)
{ int  j=get<int_value>()->val; int i=get<int_value>()->val;
  if (j==0) throw std::runtime_error("Modulo zero");
  if (l!=expression_base::no_value)
    push_value(new int_value(intutils::remainder(i,intutils::abs(j))));
}
@)
void divmod_wrapper(expression_base::level l)
{ int j=get<int_value>()->val; int i=get<int_value>()->val;
  if (j==0) throw std::runtime_error("DivMod by zero");
  if (l!=expression_base::no_value)
  { push_value(new int_value
     (j>0 ? intutils::divide(i,j) : -intutils::divide(i,-j)));
    push_value(new int_value(intutils::remainder(i,intutils::abs(j))));
    if (l==expression_base::single_value)
      wrap_tuple(2);
  }
}
@)
void unary_minus_wrapper(expression_base::level l)
{@; int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new int_value(-i)); }
@)
void power_wrapper(expression_base::level l)
{ static shared_int one(new int_value(1));
@/int n=get<int_value>()->val; shared_int i=get<int_value>();
  if (intutils::abs(i->val)!=1 and n<0)
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

@ The operator `/' will not denote integer division, but rather formation of
fractions (rational numbers). The opposite operation of separating a rational
number into numerator and denominator is also provided; it is essential in
order to be able to get from rationals back into the world of integers.

@< Local function definitions @>=

void fraction_wrapper(expression_base::level l)
{ int d=get<int_value>()->val; int n=get<int_value>()->val;
  if (d==0) throw std::runtime_error("fraction with zero denominator");
  if (l!=expression_base::no_value)
    push_value(new rat_value(arithmetic::Rational(n,d)));
}
@)

void unfraction_wrapper(expression_base::level l)
{ arithmetic::Rational q=get<rat_value>()->val;
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
{ arithmetic::Rational j=get<rat_value>()->val;
  arithmetic::Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new rat_value(i+j));
}
@)
void rat_minus_wrapper(expression_base::level l)
{ arithmetic::Rational j=get<rat_value>()->val;
  arithmetic::Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new rat_value(i-j));
}
@)
void rat_times_wrapper(expression_base::level l)
{ arithmetic::Rational j=get<rat_value>()->val;
  arithmetic::Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new rat_value(i*j));
}
@)
void rat_divide_wrapper(expression_base::level l)
{ arithmetic::Rational j=get<rat_value>()->val;
  arithmetic::Rational i=get<rat_value>()->val;
  if (j.numerator()==0)
    throw std::runtime_error("Rational division by zero");
  if (l!=expression_base::no_value)
    push_value(new rat_value(i/j));
}
@)
void rat_unary_minus_wrapper(expression_base::level l)
{@; arithmetic::Rational i=get<rat_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new rat_value(arithmetic::Rational(0)-i)); }
@)
void rat_inverse_wrapper(expression_base::level l)
{@; arithmetic::Rational i=get<rat_value>()->val;
  if (i.numerator()==0)
    throw std::runtime_error("Inverse of zero");
  if (l!=expression_base::no_value)
    push_value(new rat_value(arithmetic::Rational(1)/i)); }
@)
void rat_power_wrapper(expression_base::level l)
{ int n=get<int_value>()->val; arithmetic::Rational b=get<rat_value>()->val;
  if (b.numerator()==0 and n<0)
    throw std::runtime_error("Negative power of zero");
  if (l!=expression_base::no_value)
    push_value(new rat_value(b.power(n)));
}

@ Relational operators are of the same flavour.
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

@ To have some operations on strings, we define a function for comparing and
for concatenating them, and one for converting integers to their string
representation (of course this remains a very limited repertoire).

@< Local function definitions @>=

void string_eq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val==j->val));
}
@)
void string_neq_wrapper(expression_base::level l)
{ shared_string j=get<string_value>(); shared_string i=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val!=j->val));
}
@)
void concatenate_wrapper(expression_base::level l)
{ shared_string b=get<string_value>(); shared_string a=get<string_value>();
  if (l!=expression_base::no_value)
    push_value(new string_value(a->val+b->val));
}
void int_format_wrapper(expression_base::level l)
{ shared_int n=get<int_value>();
  std::ostringstream o; o<<n->val;
  if (l!=expression_base::no_value)
    push_value(new string_value(o.str()));
}

@ Here is a simple function that outputs any value, in the format used by the
interpreter itself. This function has an argument of unknown type; we just
pass the popped value to the |operator<<|.

This is the first place in this file where we produce user
output to a file. In general, rather than writing directly to |std::cout|, we
shall pass via a pointer whose |output_stream| value is maintained in the main
program, so that redirecting output to a different stream can be easily
implemented. Since this is a wrapper function there is no other way to convey
the output stream to be used than via a dedicated global variable.

@< Local function definitions @>=
void print_wrapper(expression_base::level l)
{ *output_stream << *pop_value() << std::endl;
  if (l==expression_base::single_value)
    wrap_tuple(0); // don't forget to return a value if asked for
}

@ Sometimes the user may want to use a stripped version of the |print| output:
no quotes in case of a string value, or no parentheses or commas in case of a
tuple value (so that a single statement can chain several texts on the same
line). The |prints_wrapper| does this down to the level of omitting quotes in
individual argument strings, using dynamic casts to determine the case that
applies.

@< Local function definitions @>=
void prints_wrapper(expression_base::level l)
{ shared_value v=pop_value();
  string_value* s=dynamic_cast<string_value*>(v.get());
  if (s!=NULL)
    *output_stream << s->val << std::endl;
  else
  { tuple_value* t=dynamic_cast<tuple_value*>(v.get());
    if (t!=NULL)
    { for (size_t i=0; i<t->val.size(); ++i)
      { s=dynamic_cast<string_value*>(t->val[i].get());
        if (s!=NULL)
	  *output_stream << s->val;
        else
           *output_stream << *t->val[i];
      }
      *output_stream << std::endl;
    }
    else
      *output_stream << *v << std::endl; // just like |print| in other cases
  }
  if (l==expression_base::single_value)
    wrap_tuple(0); // don't forget to return a value if asked for
}

@ For the size-of operator we provide several specific bindings: for strings,
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

@ The generic size-of wrapper was explicitly referenced in
section@# sizeof section @>, it is not entered into any table.

@< Local function definitions @>=
void sizeof_wrapper(expression_base::level l)
{ size_t s=get<row_value>()->val.size();
  if (l!=expression_base::no_value)
    push_value(new int_value(s));
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

@ Here are functions for adding individual elements to a row value, and for
joining two such values.

@< Local function definitions @>=
void suffix_element_wrapper(expression_base::level l)
{ shared_value e=pop_value();
  shared_row r=get_own<row_value>();
  if (l!=expression_base::no_value)
  {@; r->val.push_back(e);
    push_value(r);
  }
}
@)
void prefix_element_wrapper(expression_base::level l)
{ shared_row r=get_own<row_value>();
  shared_value e=pop_value();
  if (l!=expression_base::no_value)
  {@; r->val.insert(r->val.begin(),e);
    push_value(r);
  }
}
@)
void join_rows_wrapper(expression_base::level l)
{ shared_row y=get<row_value>();
  shared_row x=get<row_value>();
  if (l!=expression_base::no_value)
  { row_ptr result(new row_value(0));
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

@ We now define a few functions, to really exercise something, even if it is
modest, from the Atlas library. These wrapper function are not really to be
considered part of the interpreter, but a first step to its interface with the
Atlas library, which is developed in much more detail in the compilation
unit \.{built-in-types}. In fact we shall make some of these wrapper functions
externally callable, so they can be directly used from that compilation unit.

We start with vector and matrix equality comparisons.
@< Declarations of exported functions @>=
void vec_eq_wrapper (expression_base::level);
void vec_neq_wrapper (expression_base::level);
void mat_eq_wrapper (expression_base::level);
void mat_neq_wrapper (expression_base::level);

@ This is of course quite similar to what we saw for rationals, for instance.

@< Function definitions @>=
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

@ Next we have the null vector and matrices, the identity matrix and matrix
transposition.

@< Declarations of exported functions @>=
void null_vec_wrapper (expression_base::level);
void null_mat_wrapper (expression_base::level);
void id_mat_wrapper (expression_base::level);
void transpose_mat_wrapper (expression_base::level);
void transpose_vec_wrapper (expression_base::level);

@ Null vectors and matrices are particularly useful as starting values. In
addition, the latter can produce empty matrices without any (null) entries,
when either the number of rows or column is zero but the other is not; such
matrices (which are hard to obtain by other means) are good starting points
for iterations that consist of adding a number of rows or columns of equal
size, and they determine this size even if none turn out to be contributed.

Since in general built-in functions may throw exceptions (even for such simple
operations as |transposed|) we hold the pointers to local values in smart
pointers; for values popped from the stack this would in fact be hard to avoid.

@< Function definitions @>=
void null_vec_wrapper(expression_base::level lev)
{ int l=get<int_value>()->val;
  if (lev!=expression_base::no_value)
    push_value(new vector_value(matrix::Vector<int>(std::abs(l),0)));
}
@) void null_mat_wrapper(expression_base::level lev)
{ int l=get<int_value>()->val;
  int k=get<int_value>()->val;
  if (lev!=expression_base::no_value)
    push_value(new matrix_value
      (matrix::Matrix<int>(std::abs(k),std::abs(l),0)));
}
@) void id_mat_wrapper(expression_base::level l)
{ int i=get<int_value>()->val;
  if (l!=expression_base::no_value)
    push_value(new matrix_value(matrix::Matrix<int>(std::abs(i)))); // identity
}
@) void transpose_mat_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(new matrix_value(m->val.transposed()));
}
@) void transpose_vec_wrapper(expression_base::level l)
{ shared_vector v=get<vector_value>();
  if (l!=expression_base::no_value)
  { matrix_ptr m (new matrix_value(matrix::Matrix<int>(1,v->val.size())));
    for (size_t j=0; j<v->val.size(); ++j)
      m->val(0,j)=v->val[j];
    push_value(m);
  }
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
  { push_value(new vector_value(v->val.numerator()));
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
    push_value(new rational_vector_value(v0->val+v1->val));
}
@)
void ratvec_minus_wrapper(expression_base::level l)
{ shared_rational_vector v1= get<rational_vector_value>();
  shared_rational_vector v0= get<rational_vector_value>();
  if (l!=expression_base::no_value)
    push_value(new rational_vector_value(v0->val-v1->val));
}


@ Now the products between vector and/or matrices. The function |mv_prod| was
in fact our first function with more than one argument (arithmetic on integer
constants was done inside the parser at that time). We make them callable from
other compilation units.

@< Declarations of exported functions @>=
void vv_prod_wrapper (expression_base::level);
void mv_prod_wrapper (expression_base::level);
void mm_prod_wrapper (expression_base::level);
void vm_prod_wrapper(expression_base::level l);

@ For wrapper functions with multiple arguments, we must always remember that
they are to be popped from the stack in reverse order; here in fact this only
matters for error reporting.

@< Function definitions @>=
void vv_prod_wrapper(expression_base::level l)
{ shared_vector w=get<vector_value>();
  shared_vector v=get<vector_value>();
  if (v->val.size()!=w->val.size())
    throw std::runtime_error(std::string("Size mismatch ")@|
     + str(v->val.size()) + ":" + str(w->val.size()));
  if (l!=expression_base::no_value)
    push_value(new int_value(v->val.dot(w->val)));
}

@ The other product operations are very similar.

@< Function definitions @>=
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
  { bitmap::BitMap pivots=matreduc::column_echelon(M->val);
    push_value(M);
    row_ptr p_list (new row_value(0)); p_list->val.reserve(pivots.size());
    for (bitmap::BitMap::iterator it=pivots.begin(); it(); ++it)
      p_list->val.push_back(shared_value(new int_value(*it)));
    push_value(p_list);
    if (l==expression_base::single_value)
      wrap_tuple(2);
  }
}

@ And here are general functions |diagonalise| and |adapted_basis|, rather
similar to Smith normal form, but without divisibility guarantee on diagonal
entries. While |diagonalise| provides the matrices applied on the left and
right to obtain diagonal form, |adapted_basis| gives only the left factor (row
operations applied) and gives it inverted, so that this matrix
right-multiplied by the diagonal matrix has the same image as the original
matrix.

@<Local function definitions @>=
void diagonalise_wrapper(expression_base::level l)
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
install_function(int_format_wrapper,"int_format","(int->string)");
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
install_function(error_wrapper,"error","(string->)");
install_function(string_eq_wrapper,"=","(string,string->bool)");
install_function(string_neq_wrapper,"!=","(string,string->bool)");
install_function(concatenate_wrapper,"#","(string,string->string)");
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
install_function(mv_prod_wrapper,"*","(mat,vec->vec)");
install_function(mm_prod_wrapper,"*","(mat,mat->mat)");
install_function(vm_prod_wrapper,"*","(vec,mat->vec)");
install_function(echelon_wrapper,"echelon","(mat->mat,[int])");
install_function(diagonalise_wrapper,"diagonalise","(mat->vec,mat,mat)");
install_function(adapted_basis_wrapper,"adapted_basis","(mat->mat,vec)");
install_function(kernel_wrapper,"kernel","(mat->mat)");
install_function(eigen_lattice_wrapper,"eigen_lattice","(mat,int->mat)");
install_function(row_saturate_wrapper,"row_saturate","(mat->mat)");
install_function(invfact_wrapper,"inv_fact","(mat->vec)");
install_function(Smith_basis_wrapper,"Smith_basis","(mat->mat)");
install_function(Smith_wrapper,"Smith","(mat->mat,vec)");
install_function(invert_wrapper,"invert","(mat->mat,int)");

@* Operations other than evaluation of expressions.
This section will be devoted to some other interactions between user and
program that do not consist just of evaluating expressions. What will
presented is not particularly related to the evaluator, and is present here
for somewhat opportunistic purposes.

@ Applied identifiers can be introduced (or modified) by the function
|global_set_identifier|; it was declared in \.{parsetree.h} since it is called
directly by the parser, and therefore it has \Cee-linkage. We define it here
since it uses the services of the evaluator.

We allow the same possibilities in a global identifier definition as in a
local one, so we take an |id_pat| as argument. In addition we handle
definitions of overloaded function instances if |overload| is true (nonzero).
We follow the logic for type-analysis of a let-expression, and that of binding
identifiers in a user-defined function for evaluation. However we use
|analyse_types| (which reports errors) rather than calling |convert_expr|
directly. To provide some feedback to the user we report any types assigned,
but not the values.

@< Function definitions @>=
extern "C"
void global_set_identifier(id_pat pat, expr rhs, int overload)
{ using namespace atlas::interpreter;
  size_t n_id=count_identifiers(pat);
  int phase=0;
  static const char* phase_name[3] = {"type_check","evaluation","definition"};
  try
  { expression_ptr e;
    type_ptr t=analyse_types(rhs,e);
    if (not pattern_type(pat)->specialise(*t))
      @< Report that type of |rhs| does not have required structure,
         and |throw| @>
@)
    phase=1;
    bindings b(n_id);
    thread_bindings(pat,*t,b); // match identifiers and their future types

    std::vector<shared_value> v;
    v.reserve(n_id);
@/  e->eval();
    thread_components(pat,pop_value(),v);
@)
    phase=2;
    if (overload==0)
      @< Add instance of identifiers in |b| with values in |v| to
         |global_id_table| @>
    else
      @< Add instance of identifiers in |b| with values in |v| to
         |global_overload_table| @>

    std::cout << std::endl;
  }
  @< Catch block for errors thrown during a global identifier definition @>
}

@ For identifier definitions we print their names and types (paying attention
to the very common singular case), before calling |global_id_table->add|.
@< Add instance of identifiers in |b| with values in |v| to
   |global_id_table| @>=
{ if (n_id>0)
    std::cout << "Identifier";
  for (size_t i=0; i<n_id; ++i)
  { std::cout << (i==0 ? n_id==1 ? " " : "s " : ", ") @|
              << main_hash_table->name_of(b[i].first) << ": "
              << *b[i].second;
    global_id_table->add(b[i].first,v[i],copy(*b[i].second));
  }
}

@ For overloaded definitions the main difference is calling the |add| method
of |global_overload_table| instead of that of |global_id_table|. However
another difference is that overloaded definitions may be rejected because of a
conflict with an existing one, so we do not print anything before the |add|
method has successfully completed. Multiple overloaded definitions in a
single \&{set} statement are non currently allowed syntactically, which the
|assert| below tests. If such multiple definitions should be made possible
syntactically, one could introduce a loop below as in the ordinary definition
case; then however error handling would also need adaptation since a failed
definition need no be the first one, and the previous ones would need to be
either undone or not reported as failed.

@< Add instance of identifiers in |b| with values in |v| to
   |global_overload_table| @>=
{ assert(n_id=1);
  size_t old_n=global_overload_table->variants(b[0].first).size();
  global_overload_table->add(b[0].first,v[0],copy(*b[0].second));
  size_t n=global_overload_table->variants(b[0].first).size();
  if (n==old_n)
    std::cout << "Redefined ";
  else if (n==1)
    std::cout << "Defined ";
  else
    std::cout << "Added definition [" << n << "] of ";
  std::cout << main_hash_table->name_of(b[0].first) << ": " << *b[0].second;
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

@< Function definitions @>=
extern "C"
void global_declare_identifier(Hash_table::id_type id, ptr t)
{ value undef=NULL;
  const type_expr& type=*static_cast<type_p>(t);
  global_id_table->add(id,shared_value(undef),copy(type));
  std::cout << "Identifier " << main_hash_table->name_of(id)
            << " : " << type << std::endl;
}

@ Finally the user may wish to forget the value of an identifier, which the
following function achieves.

@< Function definitions @>=
extern "C"
void global_forget_identifier(Hash_table::id_type id)
{ std::cout << "Identifier " << main_hash_table->name_of(id)
            << (global_id_table->remove(id) ? " forgotten" : " not known")
            << std::endl;
}

@ Forgetting the binding of an overloaded identifier at a given type is
similar.

@< Function definitions @>=
extern "C"
void global_forget_overload(Hash_table::id_type id, ptr t)
{ const type_expr& type=*static_cast<type_p>(t);
  std::cout << "Definition of " << main_hash_table->name_of(id)
            << '@@' << type @|
            << (global_overload_table->remove(id,type)
               ? " forgotten"
               : " not known")
            << std::endl;
}

@ It is useful to print type information, either for a single expression or for
all identifiers in the table. We declare the pointer that was already used in
|print_wrapper|.

@< Declarations of global variables @>=
extern std::ostream* output_stream;
@ The |output_stream| will normally point to |std::cout|.

@< Global variable definitions @>=
std::ostream* output_stream= &std::cout;

@ The function |type_of_expr| prints the type of a single expression, without
evaluating it. Since we allows arbitrary expressions, we must cater for the
possibility of failing type analysis, in which case |analyse types| after
catching it re-throws a |std::runtime_error|. By catching and |std::exception|
we ensure ourselves against unlikely events like |bad_alloc|.

@< Function definitions @>=
extern "C"
void type_of_expr(expr e)
{ try
  {@; expression_ptr p;
    *output_stream << "type: " << *analyse_types(e,p) << std::endl;
  }
  catch (std::exception& err) {@; std::cerr<<err.what()<<std::endl; }
}

@ The function |show_overloads| has a similar purpose, namely to find out the
types of overloaded symbols. It is however much simpler, since it just has to
look into the overload table and extract the types stored there.

@< Function definitions @>=
extern "C"
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

@ The function |show_ids| prints a table of all known identifiers and their
types.

@< Function definitions @>=
extern "C"
void show_ids()
{ *output_stream << "Overloaded operators and functions:\n"
                 << *global_overload_table @|
                 << "Global values:\n" << *global_id_table;
}

@ Here is a tiny bit of global state that can be set from the main program and
inspected by any module that cares to (and that reads \.{evaluator.h}).

@< Declarations of global variables @>=
extern int verbosity;

@~By raising the value of |verbosity|, some trace of internal operations can
be activated.

@< Global variable definitions @>=
int verbosity=0;

@* Index.

% Local IspellDict: british
