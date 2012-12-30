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
\def\foreign#1{{\sl#1\/}}
\def\id{\mathop{\rm id}}

@* Outline.
%
This file describes the central part of the interpreter for the command
language of the Atlas of Lie Groups and Representation software
called \.{realex} (which could be taken to stand for Redesigned
Expression-based Atlas of Lie groups EXecutable). This part is concerned with
the analysis and execution of expressions that have already been processed by
the parser. These are highly recursive processes, and this rather large module
has been limited to those functions that play a part in this recursion. Other
more one-time matters like initialisation and setting global variables that
were originally done in this module have been relegated to a separate
module \.{global.w}.

@( evaluator.h @>=

#ifndef EVALUATOR_H
#define EVALUATOR_H

#include "types.h"

@< Includes needed in the header file @>@;
namespace atlas { namespace interpreter {
@< Type definitions @>@;
@< Declarations of global variables @>@;
@< Declarations of exported functions @>@;
}@; }@;
#endif

@ The implementation unit follows a somewhat similar pattern.

@h "evaluator.h"
@h <cstdlib>
@c
namespace atlas { namespace interpreter {
@< Global variable definitions @>@;
namespace {@;
@< Local variable definitions @>@;
@< Local function definitions @>@;
}@;
@< Function definitions @>@;
}@; }@;

@ Although initialising the evaluator will be handled in \.{global.w}, we
define a function that resets the evaluator here. The stack owns the values it
contains, but there was no reason to wrap it into a class with a destructor,
since we never intend to destroy the stack entirely: if our program exits
either peacefully or by an uncaught exception we don't care about some values
that are not destroyed. We must remember however to |delete| the values
whenever we empty the stack after catching a runtime error. We provide a
function to reset the evaluator after catching an exception that does this, as
well as other actions needed to have the evaluator start with a clear slate.

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


@* Outline of the evaluation process.
%
This module is concerned with the transformation of an abstract syntax tree as
produced by the parser into ultimately actions and computed values. This
transformation consist of two separate parts: type analysis, which also
transforms the syntax tree into a more directly executable form, and execution
of those transformed expressions.

The expression returned by the parser, of type |expr|, and the conversion to
the executable format |expression| (a type defined in \.{types.w} as a pointer
to the base class |expression_base|, from which many more specialised classes
will be derived) is performed by the function |convert_expr|. This is a large
and highly recursive function, and a large part of the current module is
dedicated to its definition. The execution of the converted value is performed
by calling the (purely) virtual method |expression_base::evaluate|, so that
the code describing the actual execution of expressions is distributed among
the many definitions of that method in derived classes, and this definition is
only implicitly (mutually) recursive through calls to the
|expression_base::evaluate| method.

@ During type checking, it may happen for certain subexpressions that a
definite type is required for them, while for others nothing is known
beforehand about their type (for instance this is the case for the complete
expression entered by the user for evaluation). The difference is important,
as in the former case conversions can be inserted to make types match, for
instance between a list of integers and a vector value; this is in fact the
only way the user can initially produce vector values. However, both cases are
handled by a single function |convert_expr|, which in addition builds (upon
success) an |expression| value. As arguments |convert_expr| takes an |expr
e@;| value produced by the parser, and a type in the form of a non-constant
reference |type_expr& type@;|. If |type| is undefined initially, then it will
be set to the type derived for the expression; if it is defined then it may
guide the conversion process, and the type eventually found will have to match
it. It could also be that |type| is initially partially defined (such as
`\.{(int,*)}', meaning ``pair on an integer and something''), in which case
derivation and testing functionality are combined; this gives flexibility to
|convert_expr|.

Upon successful completion, |type| will usually have become a completely
defined. The object |type| should be owned by the caller, who will
automatically gain ownership of any new nodes added. The latter, if it
happens, will be so due to calling the |specialise| method for |type| or for
its descendants. In some cases |type| will remain partly undefined, like for
an emtpy list display which gets type~`\.{[*]}'; however if |type| remains
completely undefined `\.*' (as would happen for the selection of a value from
an empty list, or for a function that is unconditionally recursive) then it
can be seen that evaluation cannot possibly complete without error, so we
might treat this case as a type error (and we shall do so occasionally if this
allows us to simplify our code).

@< Declarations of exported functions @>=
expression convert_expr(const expr& e, type_expr& type)
 @/ throw(std::bad_alloc,program_error);

@ In the function |convert_expr| we shall need a type for storing bindings
between identifiers and types, and this will be the |bindings| class. It uses
a vector of individual bindings, and since a vector cannot hold auto-pointers,
we need to define a destructor to clean up the types. Different such binding
vectors will be stacked, for nested scopes, but using an STL container for
that would necessitate defining a copy constructor, which would be a painful
operation: (1)~it must call |copy| on the types held in the bindings, in order
to avoid double destruction, (2)~these calls to |copy| could throw an
exception, and (3)~the constructor won't be complete, and the destructor
therefore not activated, until the final entry is copied, so we would need a
|try|\dots|catch| in the constructor to avoid a memory leak.

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
bindings in order to detemine their (lexical) binding and type. We could have
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
we don't return an |expression_ptr| is entirely pragmatic: there are many
return statements in this function, and ensuring conversion to an auto-pointer
before return would be somewhat laborious, whereas the conversion can be
easily and safely achieved by the caller; also sometimes the resulting pointer
is stored in a data structure that cannot accommodate auto-pointers.

The code below takes into account the possibility that a denotation is
converted immediately to some other type, for instance integer denotations can
be used where a rational number is expected. The function |coerce| tests for
this possibility, and may modify its final argument correspondingly.

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

@*1 Executable expression objects.
%
Let us define a first class derived from |expression_base|, which is
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
particular type is required then the components will just be required to have
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

@ When in |convert_expr| we encounter a list display when a non-row type is
expected, we single out the cases that a conversion from a row type to the
required type is available; in that case we continue to convert the component
expressions with as expected type the corresponding component type (if
multiple coercions to the required type are known, the first one in the table
gets preference; this occurs for required type \.{mat}, and means that the
component type will then be \.{vec} rather than \.{[int]}).

@< If |type| can be converted from some row-of type, check the components of
   |e.e.sublist|... @>=
{ type_expr comp_type;
  const conversion_record* conv = row_coercion(type,comp_type);
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
  bool tuple_expected=type.specialise(*tup); // whether |type| is a tuple
  tuple_expression* tup_exp;
  expression_ptr result(tup_exp=new tuple_expression(0));
  tup_exp->component.reserve(length(e.e.sublist));
  type_list tl= tuple_expected ? type.tuple : tup->tuple;
  for (expr_list el=e.e.sublist; el!=NULL; el=el->next,tl=tl->next)
    tup_exp->component.push_back(convert_expr(el->e,tl->t));
  if (tuple_expected)
    return result.release();  // and convert (derived|->|base) to |expression|
  else if (coerce(*tup,type,result))
    return result.release();
  else throw type_error(e,copy(*tup),copy(type));
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
While we have seen expressions to build lists, and to make vectors and
matrices out of them, we so far are not able to access their components once
they are constructed. To that end we shall now introduce operations to index
such values. We allow subscription of rows, but also of vectors, rational
vectors, matrices, and strings. Since after type analysis we know which of the
cases applies, we define several classes. These differ mostly by their
|evaluate| method, so we first derive an intermediate class from
|expression_base|, and derive the others from it. This class also serves to
host an enumeration type that will serve later. We include cases here that are
related to types define in \.{built-in-types}.

@< Type definitions @>=
struct subscr_base : public expression_base
{ enum sub_type
  { row_entry, vector_entry, ratvec_entry, string_char
  , matrix_entry, matrix_column, mod_poly_term };
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
  static bool assignable(sub_type);
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
struct string_subscription : public subscr_base
{ string_subscription(expression_ptr a, expression_ptr i)
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
@)
struct module_coefficient : public subscr_base
{ module_coefficient(expression_ptr pol, expression_ptr param)
  : subscr_base(pol,param) @+{}
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
{ if (aggr.kind==primitive_type)
  switch (aggr.prim)
  {  default:
  break; case vector_type: if (index==int_type)
    {@; kind=vector_entry; return subscr.specialise(int_type);  }
  break; case rational_vector_type: if(index==int_type)
    {@; kind=ratvec_entry; return subscr.specialise(rat_type); }
  break; case string_type: if(index==int_type)
    {@; kind=string_char; return subscr.specialise(str_type); }
  break; case matrix_type:
    if (index==int_int_type)
    {@; kind=matrix_entry; return subscr.specialise(int_type); }
    else if (index==int_type)
    {@; kind=matrix_column; return subscr.specialise(vec_type); }
  break; case virtual_module_type: if (index==param_type)
    {@; kind=mod_poly_term; return subscr.specialise(split_type); }
  }
  else if (aggr.kind==row_type and index==int_type)
  @/{@; kind=row_entry;
        return subscr.specialise(*aggr.component_type);
  }
  return false;
}

@ Some cases, although valid as subscriptions, do not allow a new value to be
assigned to the component value (this holds for instance selecting a character
from a string).

@< Function def... @>=
bool subscr_base::assignable(subscr_base::sub_type t)
{ switch (t)
  {@; case ratvec_entry: case string_char: case mod_poly_term: return false;
    default: return true;
  }
}


@ When encountering a subscription in |convert_expr|, we determine the types
of array and of the indexing expression separately, ignoring so far any type
required by the context. Then we look if the types agree with any of the four
types of subscription expressions that we can convert to, throwing an error if
it does not. Finally we check is the a priori type |subscr_type| of the
subscripted expression equals or specialises to the required |type|, or can be
converted to it by |coerce|, again throwing an error if nothing works. For the
indexing expression only equality of types is admitted, since the basic
language has no conversions that could apply, and if an extension does provide
some (indeed we shall later a add conversion from split integer to pair of
integers), we would not want to apply them implicitly in index positions. Also
there is little point in catering for (indexing) expressions having completely
undetermined type, as such a type can only apply to an expression that can
never return (it cannot be evaluated without error).

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
    case subscr_base::string_char:
      subscr.reset(new string_subscription(array,index));
    break;
    case subscr_base::matrix_entry:
      subscr.reset(new matrix_subscription(array,index));
    break;
    case subscr_base::matrix_column:
      subscr.reset(new matrix_slice(array,index));
    break;
    case subscr_base::mod_poly_term:
      subscr.reset(new module_coefficient(array,index));
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
{ shared_int i=(index->eval(),get<int_value>());
  shared_row r=(array->eval(),get<row_value>());
  if (static_cast<unsigned int>(i->val)>=r->val.size())
    throw std::runtime_error(range_mess(i->val,r->val.size(),this));
  push_expanded(l,r->val[i->val]);
}
@)
void vector_subscription::evaluate(level l) const
{ shared_int i=(index->eval(),get<int_value>());
  shared_vector v=(array->eval(),get<vector_value>());
  if (static_cast<unsigned int>(i->val)>=v->val.size())
    throw std::runtime_error(range_mess(i->val,v->val.size(),this));
  if (l!=no_value)
    push_value(new int_value(v->val[i->val]));
}
@)
void ratvec_subscription::evaluate(level l) const
{ shared_int i=(index->eval(),get<int_value>());
  shared_rational_vector v=(array->eval(),get<rational_vector_value>());
  if (static_cast<unsigned int>(i->val)>=v->val.size())
    throw std::runtime_error(range_mess(i->val,v->val.size(),this));
  if (l!=no_value)
    push_value(new rat_value(Rational @|
       (v->val.numerator()[i->val],v->val.denominator())));
}
@)
void string_subscription::evaluate(level l) const
{ shared_int i=(index->eval(),get<int_value>());
  shared_string s=(array->eval(),get<string_value>());
  if (static_cast<unsigned int>(i->val)>=s->val.size())
    throw std::runtime_error(range_mess(i->val,s->val.size(),this));
  if (l!=no_value)
    push_value(new string_value(s->val.substr(i->val,1)));
}

@ And here are the cases for matrix indexing.

@< Function definitions @>=
void matrix_subscription::evaluate(level l) const
{ index->multi_eval(); @+
  shared_int j=get<int_value>();
  shared_int i=get<int_value>();
  shared_matrix m=(array->eval(),get<matrix_value>());
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
{ shared_int j=(index->eval(),get<int_value>());
  shared_matrix m=(array->eval(),get<matrix_value>());
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
strict ownership. If one would instead give ownership of the |type| field
directly to individual |id_data| entries, this would greatly complicate their
duplication, and therefore their insertion into the table. We therefore only
allow construction of |id_data| objects in an empty state; the pointers they
contain should be set only \emph{after} the insertion of the object into the
table, which then immediately assumes ownership of the type.

@< Type definitions @>=

typedef std::tr1::shared_ptr<shared_value> shared_share;
struct id_data
{ shared_share val; @+ type_p type;
  id_data() : val(),type(NULL)@+ {}
};

@ We cannot store auto pointers in a table, so upon entering into the table we
convert type pointers to ordinary pointers, and this is what |type_of| lookup
will return (the table retains ownership); destruction of the type expressions
referred to will be handled explicitly when the entry, or the table itself, is
destructed.

@< Includes needed in the header file @>=
#include <map>
#include "lexer.h" // for the identifier hash table

@~Overloading is not done in this table, so a simple associative table with
the identifier as key is used.

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
    (std::string("Identifier without table entry:")
     +main_hash_table->name_of(id));
@.Identifier without value@>
  return p->second.val;
}

@ We provide a |print| member that shows the contents of the entire table.
Since identifiers might have undefined values, we must test for that condition
and print dummy output in that case. This case is distinct from the one where
a user explicitly asks for the value of an uninitialised variable; the latter
will be caught by the evaluator at the point where the variable is evaluated.

@< Function def... @>=

void Id_table::print(std::ostream& out) const
{ for (map_type::const_iterator p=table.begin(); p!=table.end(); ++p)
  { out << main_hash_table->name_of(p->first) << ": " @|
        << *p->second.type << ": ";
    if (*p->second.val==NULL)
      out << '*';
    else
      out << **p->second.val;
    out << std::endl;
   }
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
Id_table* global_id_table=NULL; // will never be |NULL| at run time

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

@~A disadvantege of using a static variable is that in case of exceptions it
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
program. We choose the easier approach to cast away the converted expression
once the matching variant is found, and to redo the analysis with the now
known result type so that coercions get inserted during conversion. But then
we again risk exponential time (although with powers of~2 rather than powers
of the number of variants); since probably coercions will not be needed at all
levels, we mitigate this risk by not redoing any work in case of an exact
match.

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
function'', which is defined in \.{global.h}.

@< Includes needed in the header file @>=

#include "global.h" // for |wrapper_function|

@ A wrapper function usually consists of a call to a library function
sandwiched between unpacking and repacking statements; in some simple cases a
wrapper function may decide to do the entire job itself.

@< Type definitions @>=

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
like the size-of operator~`\#', even in case such an operator should not occur
in the overload table.

@< Convert and |return| an overloaded function call... @>=
{ const Hash_table::id_type id =e.e.call_variant->fun.e.identifier_variant;
  const expr arg=e.e.call_variant->arg;
  size_t i,j; // dummies; local binding not used here
  if (not is_empty(arg) and id_context->lookup(id,i,j)==NULL)
  { const overload_table::variant_list& variants
      = global_overload_table->variants(id);
    if (variants.size()>0 or is_special_operator(id))
      return resolve_overload(e,type,variants);
  }
}

@ The names of some special operators are tested for each time analyse an
overloaded call; to avoid having to look them up in |main_hash_table| each
time, we store each one in a static variable inside a local function.

@< Local function definitions @>=
Hash_table::id_type equals_name()
{@; static Hash_table::id_type name=main_hash_table->match_literal("=");
  return name;
}

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

@ For overloaded function calls, once the overloading is resolved, we proceed
in a similar fashion to non-overloaded calls, except that there is no function
expression to convert (the overload table contains an already evaluated
function value, either built-in or user-defined). We deal with the built-in
case here, and will give the user-defined case later when we have discussed
the necessary value types.

As a special safety measure against the easily made error of writing `\.='
instead of an assignment operator~`\.{:=}', we forbid converting the result of
an (always overloaded) call to the equality operator to the void type,
treating this case as a type error instead. In the unlikely case that the user
defines an overloaded instance of `\.=' with void result type, calls to this
operator will still be accepted.

@< Return a call of variant |v|... @>=
{ expression_ptr call;
  builtin_value* f = dynamic_cast<builtin_value*>(v.val.get());
  if (f!=NULL)
    call = expression_ptr
      (new overloaded_builtin_call(f->val,f->print_name.c_str(),arg));
  else @< Set |call| to the call of the user-defined function |v| @>

  if (type==void_type and
      id==equals_name() and
      v.type->result_type!=void_type)
    throw type_error(e,copy(v.type->result_type),copy(type));
  return conform_types(v.type->result_type,type,call,e);
}

@ For operator symbols that satisfy |is_special_operator(id)|, we test generic
argument type patterns before we test instances in the overload table, because
the latter could otherwise mask some generic ones due to coercion. Therefore
if we fail to find a match, we simply fall through; however if we match an
argument type but fail to match the returned type, we throw a |type_error|.

The function |print| (but not |prints|) will return the value printed if
required, so it has the type of a generic identity function. This is done so
that inserting |print| around subexpressions for debugging purposes can be
done without other modifications of the user program. To ensure that coercions
that would apply in the absence of |print| still get their chance, we do a
recursive call to |coerce| after accepting without check the call to |print|.

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
  else if (id==print_name()) // this one always matches
  { expression c = new generic_builtin_call(print_wrapper,"print",arg);
    expression_ptr call(c); // get ownership
    if (type.specialise(a_priori_type) or coerce(a_priori_type,type,call))
      return call.release();
    else throw type_error(e,copy(a_priori_type),copy(type));
 }
  else if(id==prints_name()) // this always matches as well
  { expression c = new generic_builtin_call(prints_wrapper,"prints",arg);
    expression_ptr call(c); // get ownership
    if (type.specialise(void_type))
      return call.release();
    throw type_error(e,copy(void_type),copy(type));
  }
}

@*1 Wrappers that are not accessed through tables.
%
Here is a simple function that outputs any value, in the format used by the
interpreter itself. This function has an argument of unknown type; we just
pass the popped value to the |operator<<|. The function returns its argument
as result.

This is the first place in this file where we produce user output to a file.
In general, rather than writing directly to |std::cout|, we shall pass via a
pointer whose |output_stream| value is maintained in the main program, so that
redirecting output to a different stream can be easily implemented. Since this
is a wrapper function there is no other way to convey the output stream to be
used than via a dedicated global variable.

@< Local function definitions @>=
void print_wrapper(expression_base::level l)
{
  *output_stream << *execution_stack.back() << std::endl;
  if (l!=expression_base::single_value)
    push_expanded(l,pop_value()); // remove and possibly expand value
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

@ The generic size-of wrapper was explicitly referenced in
section@# sizeof section @> and never accessed in any other way. Therefore we
do not enter it into any table like wrapper functions usually are, but just
define it here.

@< Local function definitions @>=
void sizeof_wrapper(expression_base::level l)
{ size_t s=get<row_value>()->val.size();
  if (l!=expression_base::no_value)
    push_value(new int_value(s));
}

@ The operator `\#' can be used also as infix operator, to join (concatenate)
two row values of the same type or to extend one on either end by a single
element. In the former case we require that both arguments have identical row
type, in the latter case we allow the single element to be, or to be converted
to, the component type of the row value.

There is a subtlety in the order here, due to the fact that one of the
arguments could be the empty list, with undetermined row type. Then there is
actual ambiguity: with `\#' interpreted as join, the result is just the other
operand, while suffixing or prefixing to the empty list gives a singleton of
the other operand, and to some operands the empty list itself can also be
suffixed or prefixed. The simplest resolution of this ambiguity is to say that
join never applies with one of the arguments of undetermined list type, and
that suffixing is preferred over prefixing (for instance $[[2]]\#[\,]$ will
give $[[2],[\,]]$, but $[\,]\#[[2]]$ will give $[[[2]]]$). This can be
obtained by testing for suffixing before testing for join: after the former
test we know the first argument does not have type \.{[*]}, so it will match
the second argument type only if both are determined row types.

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

@ We called |can_coerce_arg| when type-checking the dyadic use of the operator
`\#', in order to see if one of its arguments can be coerced from its type
|from| to |to|. The situation there is somewhat unusual, in that we have
already converted the entire argument expression at the point where we
discover, based on the type found, that one of the arguments might need to be
coerced. This means that such a coercion must be inserted into an already
constructed expression, whereas usually it is applied on the outside of an
expression under construction (notably in ordinary overloading if coercion is
needed we simply convert the operands a second time, inserting the coercions
while doing so). Concretely this means that although we can use the |coerce|
function, it wants a reference to an auto-pointer for the expression needing
modification (so that it can manage ownership during its operation), but the
expression here is held in an ordinary pointer inside a |tuple_expression|.
Therefore we must copy the ordinary |expression| pointer temporarily to an
|expression_ptr|, simulating the auto-pointer actions manually: after
construction of the |expression_ptr| we set the |expression| (temporarily)
to~|NULL| to avoid potential double destruction, and after the call to
|coerce| we release the (possibly modified) |expression_ptr| back into the
|expression|.

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


@*1 Evaluating built-in function calls.
%
To evaluate a |call_expression| object we evaluate the function, and then
test whether it is a built-in function. In the former case we evaluate the
arguments expanded on the stack and call the built-in function, passing the
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

@ We shall need some other functions to deal with patterns, all with a
similar structure.

@< Declarations of exported functions @>=
type_ptr pattern_type(const id_pat& pat);
size_t count_identifiers(const id_pat& pat);
void list_identifiers(const id_pat& pat, std::vector<Hash_table::id_type>& d);
void thread_bindings(const id_pat& pat,const type_expr& type, bindings& dst);
void thread_components
  (const id_pat& pat,const shared_value& val, std::vector<shared_value>& dst);

@ For handling declarations with patterns as left hand side, we need a
corresponding type pattern; for instance \\{whole}:$(x,,z:(f,))$ requires the
type \.{(*,*,(*,*))}. These recursive functions construct such types.

@< Function definitions @>=
type_list_ptr pattern_list(const patlist p)
{@; return p==NULL ? type_list_ptr(NULL)
  : make_type_list(pattern_type(p->body),pattern_list(p->next));
}
@)
type_ptr pattern_type(const id_pat& pat)
{@; return (pat.kind&0x2)==0
  ? copy(unknown_type)
  : make_tuple_type(pattern_list(pat.sublist));
}

@ Here we count the number or list the identifiers in a pattern. The latter
uses an output parameter rather than a return value since this avoid doing any
concatenation of vectors.

@< Function definitions @>=
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

@< Function definitions @>=
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

@< Function definitions @>=
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
had been, possible alternative solutions would be to take a copy of it, or
else to start with an undefined type and test afterwards that it has become
equal to |bool_type|. The latter option would be slightly different by
excluding any coercions to \&{bool} in the condition. However this would not
be noticeable since (currently) no such coercions exist anyway.

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
{@; out << " while " <<  *condition << " do " << *body << " od ";
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
    expression_ptr result(new while_expression(c,b));
    return type==void_type ? new voiding(result) : result.release();
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
  const conversion_record* conv = row_coercion(type,comp_type);
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
Next we consider |for| loops over components of a value. They also have three
parts, an identifier pattern defining the loop variable(s), an in-part giving
the object whose components are iterated over, and a body that may produce a
new value for each component of the in-part. We allow iteration over vectors
and matrices, non-row types which are indexable by integers (for now this
means strings), and iteration over the terms of a parameter polynomial
(representing isotypical components of a virtual module). The syntactic
structure of the for loop is the same for all these cases.

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

@ We could not inline the following constructor definition in the class
declaration, as it uses the local function |copy_id_pat| that is not known in
the header file.

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
its a priori type.

@< Other cases for type-checking and converting... @>=
case for_expr:
{ f_loop f=e.e.for_variant;
  bindings bind(count_identifiers(f->id));
   // for identifier(s) introduced in this loop
  type_expr in_type;
  expression_ptr in_expr(convert_expr(f->in_part,in_type));  // \&{in} part
  subscr_base::sub_type which; // the kind of aggregate iterated over
  @< Set |which| according to |in_type|, and set |bind| according to the
     identifiers contained in |f->id| @>
  type_expr body_type, *btp; const conversion_record* conv=NULL;
  if (type==void_type)
    btp=&void_type; // void context is more permissive for body
  else if (type.specialise(row_of_type))
    btp=type.component_type;
  else if ((conv=row_coercion(type,body_type))!=NULL)
    btp=&body_type;
  else throw type_error(e,copy(row_of_type),copy(type));
  bind.push(id_context);
  expression_ptr body(convert_expr (f->body,*btp));
@/bind.pop(id_context);
  expression_ptr loop(new for_expression(f->id,in_expr,body,which));
@/return type==void_type ? new voiding(loop) :
      @| conv!=NULL ? new conversion(*conv,loop) : @| loop.release() ;
}

@ This type must be indexable by integers (so it is either a row-type or
vector, matrix or string), or it must be loop over the coefficients of a
polynomial. The call to |subscr_base::indexable| will set |comp_type| to the
component type resulting from such a subscription.

@< Set |which| according to |in_type|, and set |bind| according to the
   identifiers contained in |f->id| @>=
{ type_expr comp_type; const type_expr* tp;
  if (subscr_base::indexable(in_type,*(tp=&int_type),comp_type,which) @|
   or subscr_base::indexable(in_type,*(tp=&param_type),comp_type,which))
  { type_ptr pt = pattern_type(f->id), @|
      it_type=make_tuple_type(make_type_list@|
      (copy(*tp),make_type_singleton(copy(comp_type))));
    if (not pt->specialise(*it_type))
      throw expr_error(e,"Improper structure of loop variable pattern");
    thread_bindings(f->id,*it_type,bind);
  }
  else
  { std::ostringstream o;
    o << "Cannot iterate over value of type " << in_type;
    throw expr_error(e,o.str());
  }
}

@ We can start evaluating the |in_part| regardless of |kind|, but for deducing
the number of iterations we must already distinguish on |kind| to predict the
type of the in-part. A |loop_frame| is constructed that will contain the
bindings for the new variables whose scope is the loop body: the loop variable
(with as value a component of the |in_part|), any named sub-components of it,
and the optional name given to the loop index. The shared pointers held in
this frame must be different on each iteration, because the body might get
hold, through a closure that incorporates the |context| constructed below, of
a copy of those pointers; this closure should not get to see changing values
of variables during subsequent iterations. On the other hand, the syntax does
not allow users to name the entire tuple |loop_var| formed of the loop index
and the in-part component, so this pointer cannot end up in the |loop_frame|
and it may remain the same between iterations.

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
    { loop_var->val[1]=in_val->val[i]; // the row current component
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
    { loop_var->val[1].reset(new rat_value(Rational @|
        (in_val->val.numerator()[i],in_val->val.denominator())));
      @< Set |loop_var->val[0]| to |i|,... @>
    }
  }
  break;
  case subscr_base::string_char:
  { shared_string in_val = get<string_value>();
    size_t n=in_val->val.size();
    if (l!=no_value)
      result = row_ptr(new row_value(n));
    for (size_t i=0; unsigned(i)<n; ++i,loop_frame.clear())
    { loop_var->val[1].reset(new string_value(in_val->val.substr(i,1)));
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
  case subscr_base::mod_poly_term:
  @< Perform a loop over the terms of a virtual module @>
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
since any closure values in the loop body will incorporate its current
instance; there would be no point in supplying fresh pointers in |loop_var| if
they were subsequently copied to overwrite the pointers in the same |context|
object each time. Once these things have been handled, the evaluation of the
loop body is standard.

@< Set |loop_var->val[0]| to |i|,... @>=
loop_var->val[0].reset(new int_value(i)); // index; newly created each time
thread_components(pattern,loop_var,loop_frame);
execution_context.reset(new context(saved_context,loop_frame));
  // this one too
if (l==no_value)
  body->void_eval();
else
{@; body->eval(); result->val[i]=pop_value(); }

@ The loop over terms of a virtual module is slightly different, and since it
handles values defined in the modules \.{built-in-types.w} we shall include
its header file.

@h "built-in-types.h"
@< Perform a loop over the terms of a virtual module @>=
{ shared_virtual_module pol_val = get<virtual_module_value>();
  size_t n=pol_val->val.size(),i=0;
  if (l!=no_value)
    result = row_ptr(new row_value(n));
  for (repr::SR_poly::base::const_iterator it=pol_val->val.begin();
       it!=pol_val->val.end(); ++it,++i,loop_frame.clear())
  { loop_var->val[0].reset(new module_parameter_value(pol_val->rf,it->first));
    loop_var->val[1].reset(new split_int_value(it->second));
    thread_components(pattern,loop_var,loop_frame);
    execution_context.reset(new context(saved_context,loop_frame));
    if (l==no_value)
      body->void_eval();
    else
      {@; body->eval(); result->val[i]=pop_value(); }
  }
}

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
{ out << " for " << main_hash_table->name_of(id)  << ": " << *count @|
      << " from " << *bound << " do " << *body << " od ";
}
void dec_for_expression::print(std::ostream& out) const
{ out << " for " << main_hash_table->name_of(id)  << " : " << *count @|
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
  type_expr body_type, *btp; const conversion_record* conv=NULL;
  if (type==void_type)
    btp=&void_type;
  else if (type.specialise(row_of_type))
    btp=type.component_type;
  else if ((conv=row_coercion(type,body_type))!=NULL)
    btp=&body_type;
  else throw type_error(e,copy(row_of_type),copy(type));
  bind.push(id_context);
  expression_ptr body(convert_expr (c->body,*btp));
@/bind.pop(id_context);
  expression_ptr loop;
  if (c->up!=0)
    loop.reset(new inc_for_expression(c->id,count_expr,bound_expr,body));
  else
    loop.reset(new dec_for_expression(c->id,count_expr,bound_expr,body));
  return type==void_type ? new voiding(loop) : @|
         conv!=NULL ? new conversion(*conv,loop) : @|  loop.release();

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
the kind of component assignment at hand. Assignments to components of
rational vectors and of strings will be forbidden, see module
@#comp_ass_type_check@>.

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
  default: {} // remaining cases are eliminated in type analysis
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
@:comp_ass_type_check@>

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
  if (subscr_base::indexable(*aggr_t,ind_t,comp_t,kind)
      and subscr_base::assignable(kind))
  { expression_ptr r(convert_expr(rhs,comp_t));
    if (is_local)
      assign.reset(new local_component_assignment(aggr,i,d,o,r,kind));
    else
      assign.reset(new global_component_assignment(aggr,i,r,kind));

    return conform_types(comp_t,type,assign,e);
  }
  else
  { std::ostringstream o;
    o << "Cannot subscript " << *aggr_t << @| " value with index of type "
      << ind_t << " in assignment";
    throw expr_error(e,o.str());
  }
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


@* Index.

% Local IspellDict: british
