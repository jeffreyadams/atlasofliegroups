% Copyright (C) 2006-2015 Marc van Leeuwen
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
called \.{axis} (after the second cervical vertebra, the one below atlas).
This part is concerned with the analysis and execution of expressions that
have already been processed by the parser. These are highly recursive
processes, and this rather large module has been limited to those functions
that play a part in this recursion. Other more one-time matters like
initialisation and setting global variables that were originally done in this
module have been relegated to a separate module \.{global.w}.

@( axis.h @>=

#ifndef EVALUATOR_H
#define EVALUATOR_H

#include "axis-types.h"

@< Includes needed in the header file @>@;
namespace atlas { namespace interpreter {
@< Type definitions @>@;
@< Declarations of global variables @>@;
@< Declarations of exported functions @>@;
}@; }@;
#endif

@ The implementation unit follows a somewhat similar pattern.

@h "axis.h"
@h <cstdlib>
@c
namespace atlas { namespace interpreter {
@< Global variable definitions @>@;
namespace {
@< Local class definitions @>@;
@< Local variable definitions @>@;
@< Local function definitions @>@;
@< Static variable definitions that refer to local functions @>@;
}@;
@< Function definitions @>@;
}@; }@;

@ Although initialising the evaluator will be handled in \.{global.w}, we
define a function that resets the evaluator here (since it effectively
functions as destructor, in case of exceptions only, for some values that are
used in this module).

@< Declarations of exported functions @>=
void reset_evaluator ();

@~The |execution_stack|, held in a global variable defined in \.{axis-types.w},
owns the values it contains, but there was no reason to wrap it into a class
with a destructor, since we never intend to destroy the stack entirely: if our
program exits either peacefully or by an uncaught exception we don't care
about some values that are not destroyed. We must however remember to empty
the stack after catching a runtime error; since the stack holds |shared_value|
smart pointers, this will automatically clean up the objects remaining on the
stack. We provide a function to reset the evaluator after catching an
exception that does this, as well as other actions (defined later) needed to
have the evaluator start with a clear slate. We allow monitoring disappearing
values when they are cleared off the stack. As this output may not be easy to
interpret for an unsuspecting user, we limit this to the case that the
|verbosity| variable (defined in \.{global.h}) is nonzero.

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
}


@* Outline of the analysis and evaluation processes.
%
This module is concerned with the processing of the abstract syntax tree as
produced by the parser, ultimately producing actions and computed values. This
processing consist of two separate stages: type analysis, which also
transforms the syntax tree into a more directly executable form and therefore
might be called compilation, and execution of those transformed expressions.

The expression returned by the parser, of type |expr|, and the conversion to
the executable format |expression| (a type defined in \.{axis-types.w} as a
pointer to the base class |expression_base|, from which many more specialised
classes will be derived) is performed by the function |convert_expr|. This is
a large and highly recursive function, and a large part of the current module
is dedicated to its definition. The execution of the converted value is
performed by calling the (purely) virtual method |expression_base::evaluate|,
so that the code describing the actual execution of expressions is distributed
among the many definitions of that method in derived classes, and this
definition is only implicitly (mutually) recursive through calls to the
|expression_base::evaluate| method.

@ During type checking, it may happen for certain subexpressions that a
definite type is required for them, while for others nothing is known
beforehand about their type (for instance this is the case for the complete
expression entered by the user for evaluation). The difference is important,
as in the former case it is possible to insert implicit conversions to make
the types match, for instance between a list of integers (built using the
facilities of the interpreter) and a vector value (one that can be directly
used by the Atlas library); this is in fact the only way the user can
construct such vector values. However, both cases (with known or unknown
result type), and some intermediate cases (where the result type is partially
known) are handled by a single function |convert_expr|. In addition to doing
types analysis, it builds (upon success) an |expression| value. As arguments
|convert_expr| takes an |const expr& e@;| referring to a value produced by the
parser, and a type in the form of a modifiable reference |type_expr& type@;|.
If |type| is undefined initially, then it will be set to the type derived for
the expression; if it is defined then it may guide the conversion process, and
the type eventually found will have to match it. It could also be that |type|
is initially partially defined (such as `\.{(int,*)}', meaning ``pair on an
integer and something''), in which case derivation and testing functionality
are combined; this gives flexibility to |convert_expr|.

Upon successful completion, |type| will usually have become completely
defined. The object |type| should be owned by the caller, who will
automatically gain ownership of any new nodes added in the process, and
accessible from |type|. The latter, if it happens, will be caused by calls of
the |specialise| method for |type|, or for its descendants, within
|convert_expr|. Although it is a modifiable lvalue reference argument in which
normally the result of the type analysis is returned, callers should ensure
that changes to the |type| argument remain local until successful completion,
so that when the function instead terminates by throwing an exception no
modifications to |type| will leave permanent traces (any updating of tables
should be explicitly done using the value of |type| after return from
|convert_expr|). Consequently |convert_expr| is free to move from |type| in
case it throws: the stealing will not leave any noticeable effects.

In some cases |type| will remain partly undefined, like for an emtpy list
display with an unknown type, which gets specialised only to~`\.{[*]}'.
However if |type| remains completely undefined `\.*' (as will happen for an
expression that selects a value from an empty list, or that calls |error|),
then this means that evaluation cannot possibly complete without error (since
no resulting value could have all possible types at once). This situation
would also happen for a function defined recursively without any terminating
case, if it were possible to specify such a function in \.{axis} (in reality
the somewhat tedious method that currently is the only way to define recursive
functions in \.{axis} requires the type of the function to be fixed
beforehand, so this case does not occur). Therefore we might treat the case
where |convert_expr| leaves |type| completely undetermined as a type error;
currently this is not signalled as such, but occasionally we do choose to
ignore certain scenarios in which the type derived for a subexpression is
`\.*', if this allows us to simplify our code.

@< Declarations of exported functions @>=
expression_ptr convert_expr(const expr& e, type_expr& type);

@*1 Layers of lexical context.
%
In the function |convert_expr| we shall need a type for storing bindings
between identifiers and types, and this will be the |layer| class. It stores a
vector of type |layer::vec| of individual bindings, while it automatically
pushes (a pointer to) itself on a stack |layer::lexical_context| of such
vectors. Instances of |layer| shall always be stored in ordinary
(stack-allocated) variables; since their unique constructor pushes the new
object, and their destructor pops it, |layer::lexical_context| is guaranteed
to hold at all times a list of references to all current |layer| objects, from
newest to oldest. This invariant is maintained even in the presence of
exceptions that may be thrown during type analysis. Since the space on the
runtime stack occupied by a |vector| is hardly larger than that of the
reference to it, this is mostly a proof-of-concept, namely that lists can be
handled with storage for the elements in automatic variables residing in the
runtime stack. The fact that pushing and popping is automatically managed in
an exception-safe manner is also cute; however explicitly clearing the list in
|reset_evaluator| in case of an exception, as used to be done, would also
work.

We also use this structure to maintain an additional attributes telling
whether we are inside a function and how deeply nested inside loops we are.
This is to facilitate checking the legality of |break| and |return|
expressions. Correct use of these could have been left to the parser at the
price of increasing the size of the grammar. These attribute may change only
when pushing/popping a layer. In most cases such a layer was needed anyway,
but in some cases (|while| loops) a layer is added specifically to mark this
change.

@< Type def... @>=
class layer
{
  typedef containers::simple_list<layer*> list;
public:
  typedef std::vector<std::pair<id_type,type_expr> > vec;
  static list lexical_context; // the unique |layer::list| in existence
private:
  vec variable;
  BitMap constness;
  const unsigned loop_depth; // number of nested loops we are in
  const type_p return_type;
    // return type of current function, if any (non owned)
public:
  layer(const layer&) = @[delete@]; // no ordinary copy constructor
  layer& operator= (const layer&) = @[delete@]; // nor assignment operator
  layer(size_t n); // non-function non-loop layer
  layer(size_t n,type_p return_type);
    // function or (with |return_type==nullptr|) loop layer
  ~layer () @+{@; lexical_context.pop_front(); }
@)
  void add(id_type id,type_expr&& t, bool is_const);
  static const_type_p lookup
    (id_type id, size_t& depth, size_t& offset, bool& is_const);
  static void specialise (size_t depth, size_t offset,const type_expr& t);
@)
  bool empty() const @+{@; return variable.empty(); }
  std::pair<id_type,type_expr>& operator[] (size_t i)
  @+{@; return variable[i]; }
  vec::iterator begin() @+{@; return variable.begin(); }
  vec::iterator end() @+{@; return variable.end(); }
  vec::const_iterator cbegin() const @+{@; return variable.begin(); }
  vec::const_iterator cend() const @+{@; return variable.end(); }
  bool is_const (vec::const_iterator it) const
  @+{@; return constness.isMember(it-cbegin()); }
  static bool may_break(unsigned depth)
  {@; return not lexical_context.empty()
      and lexical_context.front()->loop_depth > depth; }
  static bool may_return()
  {@; return not lexical_context.empty()
      and lexical_context.front()->return_type!=nullptr; }
  static type_expr& current_return_type()
  {@; return *lexical_context.front()->return_type; }
};

@ Here are the constructors, which are used on three kinds of occasions: the
first one for \&{let} expressions, and the second one for loops (in which case
|return_type==nullptr|) and for user-defined functions (in which case
|return_type!=nullptr|).

@< Function def... @>=
layer::layer(size_t n) // non-function non-loop layer
: variable(), constness(n)
, loop_depth(lexical_context.empty() ? 0 : lexical_context.front()->loop_depth)
, return_type(lexical_context.empty() ? nullptr
             : lexical_context.front()->return_type)
{@; variable.reserve(n); lexical_context.push_front(this); }
layer::layer(size_t n,type_p return_type) // function or loop layer
: variable(), constness(n)
,@/ loop_depth(return_type!=nullptr ? 0
            : lexical_context.empty() ? 1
            : lexical_context.front()->loop_depth+1)
,@/ return_type(return_type!=nullptr ? return_type @|
             : lexical_context.empty() ? nullptr
             : lexical_context.front()->return_type)
{@; variable.reserve(n); lexical_context.push_front(this); }
@q layer (size_t n,type_p rt) : layer(n) // delegating constructor @>
@q @@+{@@; if (rt==nullptr) ++loop_depth; else loop_depth=0,return_type=rt; } @>


@ The method |add| adds a pair to the vector of bindings; the type is moved
into the |layer| object. This is also a good place to check for the presence
of identical identifiers.

@< Function def... @>=
void layer::add(id_type id,type_expr&& t,bool is_const)
{ for (auto it=variable.begin(); it!=variable.end(); ++it)
  // traverse |variable| vector
    if (it->first==id)
      throw program_error @/
       (std::string("Multiple binding of '")
                    +main_hash_table->name_of(id)
                    +"' in same scope");
  constness.set_to(variable.size(),is_const);
  variable.emplace_back( id, std::move(t) );
}

@ During conversion of expressions, we keep a stack |layer::lexical_context|
of identifier bindings in order to determine their (lexical) binding and type.
Making this a static member if the |layer| class means we cannot start an
independent call of |convert_expr| (one that does not build upon the current
lexical context) while some instance of |convert_expr| is still active. Since
we do not intend to do that, this is not a problem.

@< Global var... @>=
layer::list layer::lexical_context;

@ The method |lookup| runs through the linked list of layers and returns a
pointer to the type if a match for the identifier |id| was found, also
assigning its static binding coordinates to output arguments |depth| and
|offset|. If no match is found a null pointer is returned and the output
parameters are unchanged. We allow having a |layer| with no variables
at all for which no stack frame will correspond at all (as a frame would incur
a runtime cost for no good at all). The method |lookup| will skip such layers
without increasing the |depth| it reports.

@< Function def... @>=
const_type_p layer::lookup
  (id_type id, size_t& depth, size_t& offset, bool& is_const)
{ size_t i=0;
  for (auto range=lexical_context.cbegin(); not lexical_context.at_end(range);
       ++range)
    if (not (*range)->variable.empty())
    { for (auto it=(*range)->cbegin(); it!=(*range)->cend(); ++it)
        if (it->first==id) // then found; now set output values
          { depth=i;
            offset=it-(*range)->begin();
            is_const = (*range)->is_const(it);
            return &it->second;
          }
      ++i; // increment depth for non-empty layers only
    }
  return nullptr;
}

@ The method |specialise| serves a specific detail that could have been (and
used to be) handled by having |lookup| return |type_p| rather than
|const_type_p|. This is the possibility that a caller will afterwards need to
specialise the type found for an identifier from its uses, if the type
initially had an unknown component as in `\.{[*]}'. Rather than modifying the
look-up type, one achieves this by calling |specialise| with the |depth| and
|offset| returned by |lookup|. This method must skip without counting the same
layers that |lookup| skipped.

@< Function def... @>=

void layer::specialise (size_t depth, size_t offset,const type_expr& t)
{ auto range=lexical_context.cbegin();
  while (depth-->0)
    do ++range;
    while ((*range)->variable.empty()); // advance layer, skipping empty ones
  (**range)[offset].second.specialise(t);
}

@ The function |convert_expr| returns a owning pointer |expression_ptr| to the
result of converting the |expr| to an |expression|.

Altogether this is a quite extensive function, with as many cases in the
switch as there are variants of |expr|, and for many of those branches a
considerable amount of work to be done. It is therefore convenient to postpone
these cases and treat them one syntactic construction at the time.

@< Function definitions @>=
expression_ptr convert_expr(const expr& e, type_expr& type)
{ switch(e.kind)
  {
   @\@< Cases for type-checking and converting expression~|e| against
   |type|, all of which either |return| or |throw| a |type_error| @>
   case no_expr: assert(false);
  }
  return expression_ptr(nullptr); // keep compiler happy
}

@* Denotations, capturing values in a program.
%
A first class derived from |expression_base|, called |denotation|, simply
stores a known |value|, which it returns upon evaluation. The value may be
passed by constant reference or by rvalue reference to |shared_ptr|; in the
former case it creates an additional sharing, and in the latter case ownership
is transferred upon construction of the |denoted_value| field. The latter case
will notably apply in case the argument expression is the result of calling
|std::make_shared|.

Since the |denotation| object stores a (constant, shared) value inside it, the
evaluation simply consists of copying the pointer to the |execution_stack|.
For the first time we see the |level| argument in action, but it is actually
handled inside the call to |push_expanded|. It may be surprising that
denotations may need to be expanded as a tuple (as long as there are no such
things as denotations for complex numbers), and indeed this was not done until
the possibility to refer to the last value computed was added to the language
(see just below); that value is wrapped into a denotation, and it could well
be a tuple.

@< Type definitions @>=
struct denotation : public expression_base
{ shared_value denoted_value;
@)
  explicit denotation(const shared_value& v) : denoted_value(v) @+{}
  explicit denotation(shared_value&& v) : denoted_value(std::move(v)) @+{}
  virtual void evaluate(level l) const @+{@; push_expanded(l,denoted_value); }
  virtual void print(std::ostream& out) const
  @+{@; denoted_value->print(out); }
};

@ Here are the first examples of the conversions done in |convert_expr|. Each
time we extract a \Cpp\ value from the |expr| produced by the parser,
construct and |new|-allocate a \.{axis} value (for instance |int_value|)
from it making the pointer shared using |std::make_shared|, pass that
pointer to the |denotation| constructor, and convert the resulting pointer to
a unique pointer.

The code below takes into account the possibility that a denotation is
converted immediately to some other type, for instance integer denotations can
be used where a rational number is expected. The function |conform_types|
(defined in \.{axis-types.w}) will test whether the denotation provides or can be
converted to the required type, and may modify its final argument in the
latter case.

@< Cases for type-checking and converting... @>=
case integer_denotation:
  { expression_ptr d@|(new denotation
      (std::make_shared<int_value>(e.int_denotation_variant)));
    return conform_types(int_type,type,std::move(d),e);
  }
case string_denotation:
  { expression_ptr d@|(new denotation
      (std::make_shared<string_value>(e.str_denotation_variant)));
    return conform_types(str_type,type,std::move(d),e);
  }
case boolean_denotation:
  { expression_ptr d@|(new denotation
        (whether(e.bool_denotation_variant)));
    return conform_types(bool_type,type,std::move(d),e);
  }

@ We allow using the last value computed to be used in an expression, using
the symbol `\.\$'.

@< Declarations of global variables @>=
extern type_expr last_type;
extern shared_value last_value;

@~We set the pointers to |nullptr| here, but the |main| function will give them
more appropriate starting values.

@< Global variable definitions @>=
type_expr last_type;
shared_value last_value;

@ In some occasions, a previously computed value can be captured in an
expression (currently this applies to `\.\$', which captures the last computed
value, and of operator casts, which capture a function value from the overload
table). In those cases we shall use an expression type that is like
|denotation| so that evaluation will give back the captured value; however for
the purpose of printing (if this expression occurs inside a function body) it
is undesirable to embark on printing the whole captured value, so we derive
and override the |print| method.

@<Type definitions @>=
class capture_expression : public denotation
{ std::string print_name;
public:
  capture_expression(const shared_value& v, const std::string& name)
  : denotation(v), print_name(name) @+{}
  virtual void print(std::ostream& out) const @+{@; out << print_name; }
};

@ Upon parsing `\.\$', an |expr| value with |kind==last_value_computed| is
transmitted. Upon type-checking we capture the value in a |denotation|
structure, which may or may not be evaluated soon after; even if the value
gets captured in a function value, it will remain immutable. For printing the
expression so formed, we suppress the actual value, but record the type as
stored in |last_type| at the time of type checking.

@< Cases for type-checking and converting... @>=
case last_value_computed:
{ std::ostringstream o;
  o << '(' << last_type << ":$)"; @q$@>
@/return conform_types
    (last_type
    ,type
    ,expression_ptr(new capture_expression(last_value,o.str()))
    ,e);
}

@*1 The suicidal expression.
%
In some cases, notably for facilitating recursive definitions, it is useful to
have a placeholder expression \&{die} that will assume any type the context
requires (we even allow an undetermined type, although in its application to
recursion the type will always be defined). The evaluation of this placeholder
is not intended, and trying to evaluate is throws a runtime error; for this
reason the expression is written \&{die}. One could require that a function
call returns \&{true} by writing ``$f(...)$~\&{or die}''.

These expressions must be representable at runtime, so we define an empty
shell for them.

@< Type definitions @>=
struct shell : public expression_base
{
virtual void evaluate (level l) const;
virtual void print(std::ostream& out) const @+{@; out << " die "; }
};

@ As said above, attempting to evaluate a |shell| is suicidal.

@< Function definitions @>=
void shell::evaluate (level l) const
{@; throw runtime_error("I die"); } // our |shell| explodes


@ The main point of \&{die} is not trying to evaluate, but allowing it to
pass type checking successfully. It does so trivially.

@< Cases for type-checking and converting... @>=
case die_expr:
{@; return expression_ptr(new shell); }

@*1 The break expression.
%
The expression \&{break} allows a premature exit from any kind of loop.

@< Type definitions @>=
struct breaker: public expression_base
{ unsigned depth;
  breaker(unsigned depth) : depth(depth) @+{}
virtual void evaluate (level l) const;
virtual void print(std::ostream& out) const
{@; out << " break ";
   if (depth>0)
     out << depth << ' ';
}
};

@ The break is realised by the \Cpp\ exception mechanism. We shall make sure
it can only by used in places where it will be caught.

@< Function definitions @>=
void breaker::evaluate (level l) const
{@; throw loop_break(depth); }


@ The only check we do for \&{break} is that it occurs in a loop.

@< Cases for type-checking and converting... @>=
case break_expr:
{ if (layer::may_break(e.break_variant))
    return expression_ptr(new breaker(e.break_variant));
  throw expr_error(e,"not in the reach of a loop");
}

@*1 The return expression.
%
The expression |return exp| allows a premature exit from a user-defined
function, returning the value of the expressions |exp|.

@< Type definitions @>=
struct returner: public expression_base
{ expression_ptr exp;
  returner(expression_ptr&& exp) : exp(std::move(exp)) @+{}
virtual void evaluate (level l) const;
virtual void print(std::ostream& out) const @+
{@; out << " return " << *exp;}
};

@ The return is realised by the \Cpp\ exception mechanism. We shall make sure
it can only by used in places where it will be caught. As usual the value of
the enclosed expression is evaluated to the stack, from where in case of
success we move it into a |function_return| object and throw that.

@< Function definitions @>=
void returner::evaluate (level l) const
{@; exp->eval(); throw function_return(pop_value()); }

@ For a \&{return} expression we check that it occurs in a function body. From
the |layer| structure we get a (modifiable) reference
|layer::current_return_type()| to the return type of the current function,
which will provide the type context for the expression after~|return|.

@< Cases for type-checking and converting... @>=
case return_expr:
{ if (layer::may_return())
    return expression_ptr(new @|
      returner(convert_expr(*e.return_variant,layer::current_return_type())));
  throw expr_error(e,"not inside a function body");
}

@* Tuple displays.
%
Tuples are sequences of values of non-uniform type, usually short and with a
given sequence of types for their components; their most obvious and simple
use is for argument lists of functions. Without using these tuple types one
might define function types taking multiple arguments, but then a function
could still only yield a single value. With tuple types multiple results can
be produced, and there will be no need to explicitly cater for functions with
multiple arguments.

Tuple values can be produced by tuple displays, which list explicitly their
component expressions. After type-checking, they are given by a
|tuple_expression| object (the name |tuple_display| was already taken).


@< Type definitions @>=
struct tuple_expression : public expression_base
{ std::vector<expression_ptr> component;
@)
#ifdef incompletecpp11
  explicit tuple_expression(size_t n) : component()
    // avoid copying result of single |expression_ptr()|
  { component.reserve(n);
    while (n-->0)
      component.push_back(expression_ptr());
  }
#else
  explicit tuple_expression(size_t n) : component(n) @+{}
#endif
   // always start out with null pointers
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ When we print a tuple display, we just print the component expressions,
enclosed in parentheses and separated by commas, to match their input syntax.


@< Function def... @>=
void tuple_expression::print(std::ostream& out) const
{ out << '(';
  for (auto it=component.begin(); it!=component.end(); ++it)
    out << (it==component.begin() ? "" : ",") << **it;
  out << ')';
}

@ When converting a tuple expression, we first try to specialise |type| to a
tuple type with the right number of unknown components; unless |type| was
completely undetermined, this just amounts to a test that it is a tuple type
with the right number of components. Even if the specialisation fails, it need
not be a type error: the tuple might be convertible to the proper type.
(Currently the conversion from a pair of integers to a split integer is the
unique conversion of this kind that is defined.) For this reason we continue
to type-check components even when |tuple_expected| is false, using the a
priori (tuple) type found as starting point for a possible coercion. In this
case we already know that |type| cannot be specialised to |*tup|, so we call
|coerce| directly rather than using |conform_types| as happens elsewhere. If
we do find a type error due to a wrong number of components, this is reported
via a type pattern; this is probably as clear as mentioning too few or too
many components explicitly.

@< Cases for type-checking and converting... @>=
case tuple_display:
{ type_expr tup=unknown_tuple(length(e.sublist));
  std::unique_ptr<tuple_expression> tup_exp(new tuple_expression(0));
  std::vector<expression_ptr>& comp = tup_exp->component;
  comp.reserve(length(e.sublist));
  bool tuple_expected = type.specialise(tup);
  // whether |type| is a tuple of correct size
  wtl_iterator tl_it (tuple_expected ? type.tupple : tup.tupple);
  for (wel_const_iterator it(e.sublist); not it.at_end(); ++it,++tl_it)
  { comp.push_back(convert_expr(*it,*tl_it));
    if (*tl_it==void_type and not is_empty(*it))
      comp.back().reset(new voiding(std::move(comp.back())));
  }
  expression_ptr result(std::move(tup_exp));
  if (tuple_expected or coerce(tup,type,result))
    return result;
  throw type_error(e,std::move(tup),std::move(type));
}

@*1 Evaluating tuple displays.
%
Evaluating a tuple display evaluates its components in a simple loop. If
|l==no_value| this is done for side effects only (a rare case), which is
achieved by calling the |void_eval| method of the component expressions.
If a single value is to be produced we prepare the result tuple at the right
size (but filled with null pointers), then in a loop evaluate the component
values using the |eval| method that leaves a single value on the stack, which
we move (by |pop_value|) into to corresponding slot of |result|; when done,
the |result| is pushed onto the stack. Finally when |l=multi_eval| we proceed
similarly but without preparing a result; we simple leave the component values
on the execution stack.

The final case is one of the rare occasions where we leave values sitting on
the |execution_stack| for some time, and the only one where a |break| might
occur at such a point (this should be quite rare, but it is possible). It is
therefore easier to clean the stack up in the code below, than to mark the
execution stack before entering any loop that might terminate with |break|,
and resetting it to the marked position when the |loop_break| is caught (in
particular since there are many kinds of loops that would need consideration).
However, while we thus localise the clean-up in a single place, this code does
get executed very often (every time a built-in function is called, except if
it has exactly $1$ argument). We therefore take care (by not marking anything
at loop entry) that in the non-throwing case (the vast majority) no cycles are
wasted at all (assuming, as seems reasonable, that simply entering a
|try|-block does not involve any work).

@< Function def... @>=
void tuple_expression::evaluate(level l) const
{ switch(l)
  {
  case no_value:
    for (auto it=component.begin(); it!=component.end(); ++it)
      (*it)->void_eval();
    break;
  case single_value:
    { auto result = std::make_shared<tuple_value>(component.size());
      auto dst_it = result->val.begin();
      for (auto it=component.cbegin(); it!=component.cend(); ++it,++dst_it)
        {@; (*it)->eval(); *dst_it=pop_value(); }
      push_value(result);
    } break;
  case multi_value:
    { auto it=component.begin();
      try
      {@; for (; it!=component.end(); ++it)
          (*it)->eval();
      }
      catch (const loop_break&) // clean up execution stack locally
      { for (; it!=component.begin(); --it)
          execution_stack.pop_back();
        throw; // propagate the |break|
      }
      catch (const function_return&) // clean up execution stack locally
      { for (; it!=component.begin(); --it)
          execution_stack.pop_back();
        throw; // propagate the |break|
      }
    }
  } // |switch(l)|
}

@* List displays.
%
Besides tuples we have rows, which are sequences of values of a uniform type.
Row values can be given by list displays, expressions that build a list of
values by constructing each component by an explicit expression. We start by
defining the structure used to represent list displays after conversion in the
type check. The logical name |list_display| for this structure is already
taken, so we call it a |list_expression| instead. Syntactically these are very
much like tuple displays, so we derive this type from |tuple_expression|.

@< Type definitions @>=
struct list_expression : public tuple_expression
{ explicit list_expression(size_t n) : tuple_expression(n)@+{}
@)virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ When we print a list display, we just print the component expressions,
enclosed in brackets and separated by commas, just as in their input syntax.

@< Function def... @>=
void list_expression::print(std::ostream& out) const
{ out << '[';
  for (auto it=component.begin(); it!=component.end(); ++it)
    out << (it==component.begin() ? "" : ",") << **it;
  out << ']';
}


@ The evaluation of a |list_expression| evaluates the components in a simple
loop. If |l==no_value|, we only evaluate for side effects, otherwise we wrap
the result in a single |row_value|. Since evaluation of component expressions
pushes the resulting value onto the execution stack, we pop each one off
immediately to integrate it into the result. We take care to hold the partial
result via a unique-pointer |result|, so that in case of a runtime error during
the evaluation of one of the component expressions the values already computed
are cleaned up.

@< Function def... @>=
void list_expression::evaluate(level l) const
{ if (l==no_value)
    for (auto it=component.begin(); it!=component.end(); ++it)
      (*it)->void_eval();
  else
  { own_row result= std::make_shared<row_value>(0);
    result->val.reserve(component.size());
    for (auto it=component.begin(); it!=component.end(); ++it)
  @/{@;
      (*it)->eval();
      result->val.push_back(pop_value());
    }
    push_value(std::move(result)); // result will be shared from here on
  }
}

@*1 Type balancing.
%
Type-checking of list displays involves a difficulty that also exists for
conditional expressions and possibly other cases: all component expressions
need to be of the same type, but we might not know which. We adopt the rule
(borrowed from Algol~68) that at least one of the components gets the common
type when converted in the context of the original type (pattern), while the
others can be converted in the context of that common type (possibly using
coercions). This process is called type balancing (the mental image is
comparing the ``weights'' of the component expression types, to see which one
is the firmest in determining the result type).

@ Here is the general set-up for balancing. We try to find in |common| a type
to which all branches can conform. Branch types incomparable with the current
value of |common| are put aside in |conflicts|. At the end, |common| may have
become broad enough to accommodate them (for instance |void| can accommodate
every possible type). So we prune |conflict| before possibly reporting an
error.

In case of success, any branches that were originally found to have a
different (narrower) type, or whose conversion threw a balancing error inside
the branch, are converted again, replacing a possible previous result in
|components|. Branches that threw a |balance_error| are certain to satisfy the
test |comp_type[i]!=target| below (which in their case means that |target| has
changed since a copy was taken to initialise |comp_type[i]|), since such an
error will have contributed to |conflicts| type that were strict
specialisations of the original |target|, and which cannot have been pruned
unless |target| was changed to a broader type since. The new conversion may
insert coercions that were absent in the original conversion (it may also
throw some other type error if conversion in that context turns out to be
impossible after all).

@< Local function definitions @>=

void balance
   ( type_expr& target // component type required from context
   , raw_expr_list elist // list of expression to be balanced
   , const expr& e // containing expression, for error reporting
   , const char* description // what kind of components, for error message
   , std::vector<expression_ptr>& components // output, converted expressions
   )
{
  unsigned n = length(elist); components.reserve(n);
  std::vector<type_expr> comp_type; comp_type.reserve(n);
  type_expr common;
    // greatest common denominator that branch types convert to\dots
  containers::sl_list<type_expr> conflicts;
    // except those branch types that are put aside here
  @< Convert each expression in |elist| in the context of a copy of |target|,
     pushing the results to |components|; maintain |common| as balancing type,
     record in |conflicts| non conforming component types @>
  @< Prune from |conflicts| any types that now test narrower than |common|
     and if nothing is left specialise |target| to |common|;
     otherwise |throw| a |balance_error| mentioning |common| and |conflicts| @>

  wel_const_iterator it(elist);
  for (unsigned i=0; i<n; ++i,++it)
    if (comp_type[i]!=target)
      components[i] = convert_expr(*it,target);
      // redo conversion with broader |common| type
}

@ We try to maintain |common| as the maximal type between the branches for the
|broader_eq| relation. If this fails due to incomparable types we move the
non-conforming type to |conflicts|. That list also collects types from
|balance_error| if thrown directly by one of the calls to |convert_expr|.
Since such type collections have internal incompatibilities, they never
provide the common type; no comparison with |common| is needed.

When catching a |balance_error|, we re-|throw| if the error was produced in a
subexpression of the branch, as this indicates an error independent of out
balancing. If that branch was a list display, the reported types are somewhat
laboriously wrapped in a ``row-of'' to produce the component type that
interests us here.

@< Convert each expression in |elist| in the context... @>=
for (wel_const_iterator it(elist); not it.at_end(); ++it)
{ try
  { comp_type.push_back(target.copy());
        // start each with a copy of original |target| type
    type_expr& ctype =comp_type.back(); // call that copy |ctype|
    components.push_back(expression_ptr());
       // push, whether or not |convert_expr| succeeds
    components.back()=convert_expr(*it,ctype);
    if (not broader_eq(common,ctype))
      { if (broader_eq(ctype,common))
          common = ctype.copy();
        else
          conflicts.push_back(ctype.copy());
          // record type not convertible to |common|
      }
  }
  catch (balance_error& err)
  { if (&err.offender!=&*it) // only incorporate top-level balancing errors
      throw; // any deeper error is propagated to be reported
    else if (err.offender.kind==list_display)
      // then wrap variants in row-of
      for (auto it=err.variants.wbegin(); not err.variants.at_end(it); ++it)
        it->set_from(type_expr(type_ptr(new type_expr(std::move(*it)))));
        // row-of |*it|
    conflicts.append(std::move(err.variants)); // then join to our |conflicts|
  }
}

@ Pruning is quite simple, and gives us an occasion to exercise the |erase|
method of |containers::sl_list|. In such loops one should not forget
to \emph{not advance} the iterator in case a node is erased in front of it.

Only if if at least one conflicting type remains do we report an error; if so,
the type |common| is added as first type to the error object, so that one has
a complete list of types that caused to balancing to fail.


@< Prune from |conflicts| any types... @>=
{ for (auto it=conflicts.cbegin(); not conflicts.at_end(it); )
    if (broader_eq(common,*it))
      conflicts.erase(it);
    else
      ++it;
  if (conflicts.empty())
     // then balancing succeeded, so set |target| to |common|
  { bool success = target.specialise(common); ndebug_use(success);
    assert(success);
    // since |common| was obtained by |convert_expr| from some copy of |target|
  }
  else
  { balance_error err(e);
    (err.message += " between ") += description ;
    err.variants.push_back(std::move(common));
    err.variants.append(std::move(conflicts));
    throw std::move(err);
  }
}

@ With balancing implemented, converting a list display become fairly easy.
The simplest case is one where a |type| a row type (or |undefined_type| that
can be specialised to such). In that case we prepare an initially empty
|list_expression|, then call |balance| with the component type of |type|,
which if successful will have converted to component types into our
|list_expression|, and it remains to |return| that object.

The case where a list display occurs in a void context is rare but valid.
(Since no list will be created, even though all component expression will be
evaluated, the choice of writing a list display is rather curious.) For it we
perform balancing with an undetermined component type (as if the display were
in undetermined type context).

In the remaining case we call |row_coercion| (defined in \.{axis-types.w}) to
see if a coercion to |type| from some row of |comp_type| exists. If this is
successful, we continue balancing with |comp_type| as expected component type.
As before balancing performs the conversion of the component expressions into
|result|, which here we wrap up in the appropriate |conversion| found by
|row_coercion|.

Finally if nothing works we report a type error with ``found type'' \.{[*]}
(we cannot be more specific; indeed none of the component expressions has been
analysed at this point).

@:list display conversion@>

@< Cases for type-checking and converting... @>=
case list_display:
{ std::unique_ptr<list_expression> result (new list_expression(0));
  auto& comps = result->component;
  result->component.reserve(length(e.sublist));
@/static const char* const str = "components of list expression";
  if (type.specialise(row_of_type))
  {
    balance(*type.component_type,e.sublist,e,str,comps);
    if (*type.component_type==void_type)
      @< Insert voiding coercions into members of |comps| that need it @>
    return std::move(result);
  }
@)
  type_expr comp_type;
  if (type==void_type) // in void context leave undetermined target type
  { balance(comp_type,e.sublist,e,str,result->component);
    return std::move(result);
    // and forget |comp_type|
  }
@)
  const conversion_record* conv = row_coercion(type,comp_type);
  if (conv!=nullptr)
  { balance(comp_type,e.sublist,e,str,result->component);
    return expression_ptr(new
      conversion(*conv,expression_ptr(std::move(result))));
  }
@)
  throw type_error(e,row_of_type.copy(),std::move(type));
  // |type| incompatible with any list
}

@ It is rare but legal to have a list display of type row-of-void (though it
could happen unintentionally if balancing found one of the components to have
void type). In any case, we must ensure that in such cases the component
expressions get a voiding coercion inserted, to ensure the invariant that
subexpressions in void context always get evaluated with |no_value| set.
However, as before for tuple display, we exempt the empty display ``()'', which
does not need, nor would it benefit from, being evaluated with |l=no_value|.

@< Insert voiding coercions into members of |comps| that need it @>=
{ wel_const_iterator it(e.sublist);
  for (auto cit=comps.begin(); cit!=comps.end(); ++it,++cit)
    if (not is_empty(*it))
      cit->reset(new voiding(std::move(*cit)));
}

@* Identifiers.
%
Identifiers are used to access values of all types, and also for designating
overloaded functions. In the latter usage a single identifier can be used to
access (depending on its immediate context) one of many values, and it
requires a certain amount of work to determine which one; this will be
discussed later, so for the moment we stick to simple applied identifiers that
identify the closest defining occurrence of that identifier in the current
lexical context. This identification can result in two outcomes: it may be
bound to a local or to a global name, which two cases are treated in fairly
different way. In particular after type analysis the two cases are converted
into different kinds of |expression|. The most fundamental difference is that
for global identifiers, the value (object) it refers to is already known at
the time the identifier expression is type-checked; for local identifiers the
|lookup| just find a position in some layer, but that will correspond to
different objects at different executions of the expression. However in either
case the type of the identifier is known at analysis time, and will not change.

Global identifiers values will be stored in a global identifier table holding
values and their types. The values of local identifiers will be stored at
runtime in a stack |frame::current| of variable bindings, but not their types
(these are held in |layer| values during type analysis, and have disappeared
at evaluation time).

In spite of these differences there is some common ground for global and local
identifiers: they have a name that can be printed. For this reason we derive
an intermediate structure from |expression_base| that will serve as base for
both kinds of applied identifier expressions.

@< Type definitions @>=
#ifdef incompletecpp11
#define nothing_new_here {} // patch for gcc 4.6
#else
#define nothing_new_here @[@[@]=@[default@]
#endif

struct identifier : public expression_base
{ id_type code;
@)
  explicit identifier(id_type id) : code(id) @+{}
  virtual ~@[identifier() nothing_new_here@];
  const char* name() const;
  virtual void print(std::ostream& out) const;
};

@ To print an identifier, we get its name from the main hash table.

@< Function definitions @>=
const char* identifier::name() const
@+{@; return main_hash_table->name_of(code); }
void identifier::print(std::ostream& out) const
@+{@; out<< name(); }

@*1 Global identifiers.
%
When during type checking an identifiers binds to a value in the global
identifier table, it will be converted into a |global_identifier| object.
Since a value is already available at this time, we can record the location of
the (pointer to the shared) value in the |global_identifier| object. Apart
from avoiding look-up at evaluation time, this binding can remain intact in
case the global identifier should be defined anew using the |add| method of
the identifier table, and will continue to access the old value. This
precaution is necessary because (in contrast to assignment) a new definition
may change the type of the identifier, but the applied identifier expression
has already been type-checked and should not be allowed to return a value of a
different type than it did originally.

@< Type definitions @>=
class global_identifier : public identifier
{ const shared_share address;
public:
  explicit global_identifier(id_type id);
  virtual ~@[global_identifier() nothing_new_here@];
  virtual void evaluate(level l) const;
};

@ The constructor for |global_identifier::evaluate| locates the value
associated to the identifier in the global identifier table.

@< Function definitions @>=
global_identifier::global_identifier(id_type id)
: identifier(id), address(global_id_table->address_of(id))
@+{}


@ Evaluating a global identifier returns the value currently stored in the
location |address|, possibly expanded if |l==multi_value|, or nothing at all
if |l==no_value|. However, since undefined global variables have been made
possible in the language, we have to watch out for a (shared) null pointer at
|*address|.

@< Function definitions @>=
void global_identifier::evaluate(level l) const
{ if (address->get()==nullptr)
  { std::ostringstream o;
    o << "Taking value of uninitialized variable '" << name() << '\'';
    throw runtime_error(o.str());
  }
  push_expanded(l,*address);
}

@*1 Local identifiers.
%
Local identifiers will be accessed from the current execution context, which
is a stack of variable bindings independent of the |execution_stack| (the
latter being used for anonymous components of expressions being evaluated).
This stack is implemented as a singly linked list, and accessed through a
shared pointer. It is not an instance of |containers::simple_list| mainly
because of sharing of parts between different contexts, which arises when
closures are formed as will be described later. The structure pointed to by
|shared_context| is described in \.{axis-types.w}; essentially, each node of the
list is a vector of values associated with identifiers introduced in the same
lexical |layer|. Although each value is associated with an identifier, they
are stored anonymously; the proper location of an applied identifier is
determined by its position in the list of lexical layers at the time of type
checking, and recorded as a pair of a relative depth (of the defining
occurrence with respect to the applied occurrence) and an offset within the
layer.

Thus using applied identifiers requires no looking up at run time, although
traversing of the linked list of corresponding |frame|s, up to the specified
depth, is necessary. One might imagine keeping a stack of |layer| pointers
cached to speed up the evaluation of applied identifiers of large depth, but
such a cache would have to be renewed at each context switch, such as those
that occur when calling or returning from a user-defined function; it is
doubtful whether this would actually result in more rapid evaluation.

The pointer holding the current execution context is declared a as static
variable of a local class |frame| to be detailed in section @#frame class@>.
Having a static variable has the advantage, compared with making it a
parameter to the |evaluate| methods, of not encumbering the numerous such
methods that neither use nor modify the context in any way (those not
involving identifiers or user defined functions).

@< Local var... @>=
shared_context frame::current; // points to topmost current frame

@ We derive the class of local identifiers from that of global ones, which
takes care of its |print| method. The data stored are |depth| identifying a
layer, and |offset| locating the proper values within the layer.

@< Type definitions @>=
class local_identifier : public identifier
{ size_t depth, offset;
public:
  explicit local_identifier(id_type id, size_t i, size_t j)
     : identifier(id), depth(i), offset(j) @+{}
  virtual void evaluate(level l) const; // only this method is redefined
};

@ The method |local_identifier::evaluate| looks up a value in the evaluation
context |frame::current|, by calling the method |evaluation_context::elem|.

@< Function definitions @>=
void local_identifier::evaluate(level l) const
{@; push_expanded(l,frame::current->elem(depth,offset)); }

@ When type-checking an applied identifier, we first look in
|layer::lexical_context| for a binding of the identifier; if found it will be
a local identifier, and otherwise we look in |global_id_table|. If found in
either way, the associated type must equal the expected type (if any), or be
convertible to it using |coerce|.

There is a subtlety in that the identifier may have a more general type than
|type| required by the context (for instance if it was \&{let} equal to an
empty list, and a concrete type of list is required). In this case the first
call of |specialise| below succeeds without making |type| equal to the
identifier type |*id_t|, and if this happens we specialise the latter instead
to |type|, using the |specialise| method either of the |layer| class (a static
method) or of |global_id_table|. This ensures that the same local identifier
cannot be subsequently used with an incompatible specialisation (notably any
further assignments to the variable must respect the more specific type). It
remains a rare circumstance that an applied occurrence (rather than an
assignment) of a local identifier specialises its type; it could happen if the
identifier is used in a cast. However type safety requires that we always
record the type to which the identifier value was specialised, since if one
allows different specialisations of the same identifier type to be made in
different subexpressions, then a devious program can manage to exploit this to
get false type predictions.

@< Cases for type-checking and converting... @>=
case applied_identifier:
{ const id_type id=e.identifier_variant;
  const_type_p id_t; size_t i,j; bool is_const;
  const bool is_local=(id_t=layer::lookup(id,i,j,is_const))!=nullptr;
  if (not is_local and (id_t=global_id_table->type_of(id,is_const))==nullptr)
  {
    std::ostringstream o;
    o << "Undefined identifier '" << main_hash_table->name_of(id) << '\'';
    if (e.loc.file!=Hash_table::empty)
      o << ' ' << e.loc;
    throw program_error (o.str());
  }
@.Undefined identifier@>
  expression_ptr id_expr = @| is_local
  ? expression_ptr(new local_identifier(id,i,j))
  : expression_ptr(new global_identifier(id));
  if (type.specialise(*id_t)) // then required type admits known identifier type
    { if (type!=*id_t)
      // usage has made type of identifier more specialised
      { if (is_local)
          layer::specialise(i,j,type);
          // then refine its type in the local context
        else
          global_id_table->specialise(id,type); // or in |global_id_table|
      }
      return id_expr;
    }
  else if (coerce(*id_t,type,id_expr))
    return id_expr;
  throw type_error(e,id_t->copy(),std::move(type));
}

@*1 Resolution of operator and function overloading.
%
While the simple mechanism of identifier identification given above suffices
for values, it would be very restrictive in case of operators (since it
allows only one function, with fixed argument types, to be bound globally to
an operator symbol identifier), and to a somewhat lesser measure for functions.
We definitely want to allow operator overloading (defining the same
operator for different combinations of argument types), and with such a
mechanism in place, it is easy to allow function overloading as well, which
turns out to be very convenient. In our discussion below we shall talk about
operators and operands, but everything applies to functions and arguments as
well.

We shall call |resolve_overload| from the case for function applications in
|convert_expr|, after testing that a non empty set of overloads exists.
Therefore the caller can pass the relevant list of |variants| as a parameter.

@< Declarations of exported functions @>=
expression_ptr resolve_overload
  (const expr& e,
   type_expr& type,
   const overload_table::variant_list& variants);

@ To resolve overloading, we used to try each variant, recursively calling
|convert_expr| to convert the arguments under the hypothesis that they could
be so in the strong type context of the type expected by that variant, but
catching (type) errors, interpreting them simply as an indication that this
variant is not the right one. However in long formulae operator applications
are nested, and |convert_expr| would call |resolve_overload| one level down,
leading to a recursive repetition of the same try-all-variants scenario. This
led to a search tree growing exponentially with the size of the formula, and
unacceptably long analysis times. So instead the rules were made slightly
stricter: every operand combination must be such that is can be analysed in
isolation, giving rise to an \foreign{a priori} type, which might however not
quite be the actual type for the intended variant. Overload resolution will
then try to find a variant, whose expected type is close enough to the
found \foreign{a priori} type to expect that coercion will allow making the
match. The necessary coercions might however have to ``creep inside'' the
argument expression (most obviously so in the common cases of an operand pair
or argument tuple) on order to apply. Since matching was based not on
the \emph{structure} of the argument expression but on type only, there
remains a possibility that coercion is not possible after all; such a case
will be considered an error rather than just a failed match (additional
language restrictions are imposed to ensure that such an error will not
mask a possible successful match).

To implement the new rules we shall, after having found a variant where the
expected operand type is unequal to the found operand type, but where
|is_close| reports that a coercion may be possible, redo the conversion of the
operand expression in the context of the expected type of the variant. Doing
so we throw away the previously converted argument expression. In principle
this again gives a potential exponential growth of the search tree, albeit
with base$~2$ rather than the number of variants: both the original and
the new call of |convert_expr| can recursively call |resolve_overload|. In
practice this is not a problem, as it requires a deeply nested formula with
coercions required at all levels, which would be a highly artificial
situation.

@ The matching condition is that the call |is_close(a_priori_type,arg_type)|
sets the bit indicating a possible conversion from the |a_priori_type| to the
expected |arg_type|. (Since |is_close| treats |void| as just the $0$-tuple
type, this means that an instance with |arg_type==void_type| (no arguments)
will not match any call with an argument (unless it has void type), even
though any type can be voided to |void|; this is intended behaviour.)

Apart from the cases listed in |variants|, there are also cases that match a
``generic'' operation (in fact an operation that has a second order type, but
our language does not have such types yet); for instance the operator~`\#'
represents several generic operations related to general row types, like
taking the size or adding an element. Contrary to ordinary overloading, we
require (with one exception) such generic operations to have |a_priori_type|
exactly matching their argument type pattern in order to apply, in other words
we allow no coercion in their arguments. The exception is |print| which should
behave transparently: its argument is matched as if the |print| was absent, so
can contain coercions if they would be there with |print| absent.

We test for generic function after the (more specifically typed) operations
with the same name from the overload table, so that a user definition can
override a generic operator for specific cases. This however only applies when
the table produced an exact match, since we consider a (necessarily exact)
generic match better than a specific variant that requires a coercion. For
this reason our table matching makes a distinction between exact matches, for
which a result is immediately returned, and inexact matches, for which we
temporarily store the call expression for the operator and its resulting type,
and which will only become the final result if no generic match is found. It
is important that the final call to |conform_types| is postponed in case of an
inexact match, because it might irreversibly specialise |type| according to
the tentative match, and this should be avoided if a generic match turns out
to override this match.

The parts of this function that actually construct a function calls are
postponed to be detailed later.

@:resolve_overload@>

@< Function definitions @>=
expression_ptr resolve_overload
  (const expr& e,
   type_expr& type,
   const overload_table::variant_list& variants)
{ const expr& args = e.call_variant->arg;
  type_expr apt; // modifiable temporary
  expression_ptr arg = convert_expr(args,apt);
    // get \foreign{a priori} type for argument
  const type_expr& a_priori_type = apt;
  id_type id =  e.call_variant->fun.identifier_variant;
  expression_ptr call(nullptr);
    // will be assigned in case of match requiring coercion
  const_type_p result_type = nullptr; // corresponding result type from table
  for (auto it=variants.begin(); it!=variants.end(); ++it)
  { const overload_data& v=*it;
    const auto& arg_type=v.type().arg_type;
    if (a_priori_type==arg_type) // exact match
    { @< Set |call| to a call of the function value |*v.val| with argument
         moved from |arg| @>
      return conform_types(v.type().result_type,type,std::move(call),e);
      // exact match, return
    }
    else
    if ((is_close(a_priori_type,arg_type)&0x1)!=0)
    {
// inexact match, so tentatively convert again using |arg_type|, hoping coercion helps
      try
      { expression_ptr arg=convert_expr(args,as_lvalue(arg_type.copy()));
        // redo conversion
        @< Set |call| to a call of the function value |*v.val|... @>
        result_type = &v.type().result_type;
        // just store the |call| and the its |result_type| now
      }
      catch (const type_error&) @+{}
       // if coercion fails ignore the match, but quit search anyway
      break;
    }
  }
  @< If |id| is a special operator like \# and it matches
     |a_priori_type|, |return| a call |id(arg)| @>
  if (result_type!=nullptr)
    // return inexact match if special operators did not do exact match
    return conform_types(*result_type,type,std::move(call),e);
  @< Complain about failing overload resolution @>
}

@ Here is the final part of |resolve_overload|, reached when no valid match
could be found. In that case we |throw| a |expr_error| explaining the is
matching identifier and type. As most function definitions will be in the
overload table even if only one definition is present, we produce a more
specific |type_error| in that case, whose message will mention the unique
expected argument type.

@< Complain about failing overload resolution @>=
if (variants.size()==1)
  throw type_error(args,a_priori_type.copy(),variants[0].type().arg_type.copy());
else
{ std::ostringstream o;
  o << "Failed to match '"
    << main_hash_table->name_of(id) @|
    << "' with argument type "
    << a_priori_type;
  throw expr_error(e,o.str());
}

@ Here we define |is_special_operator|, a function called when a function or
operator name is absent from the overload, so see if |resolve_overload| should
be called anyway. It must compare against the numeric codes of these
identifiers, which are not known explicitly at compile time, but which will
not change once the tables are initialised. To avoid having to look up these
codes in |main_hash_table| each time, we store each one in a static variable
inside a dedicated local function. While the equality operator is not a
generic one, it is tested for elsewhere, so we also give it its local function
with static variable.

@h "lexer.h" // for |main_hash_table|

@< Local function definitions @>=
id_type size_of_name()
{@; static id_type name=main_hash_table->match_literal("#");
  return name;
}
id_type concatenate_name()
{@; static id_type name=main_hash_table->match_literal("##");
  return name;
}
id_type print_name()
{@; static id_type name=main_hash_table->match_literal("print");
  return name;
}
id_type to_string_name()
{@; static id_type name=main_hash_table->match_literal("to_string");
  return name;
}
id_type prints_name()
{@; static id_type name=main_hash_table->match_literal("prints");
  return name;
}
id_type error_name()
{@; static id_type name=main_hash_table->match_literal("error");
  return name;
}
@)
inline bool is_special_operator(id_type id)
{@; return id==size_of_name()
    @|  or id==concatenate_name()
    @|  or id==print_name()
    @|  or id==to_string_name()
    @|  or id==prints_name()
    @|  or id==error_name(); }
@)
id_type equals_name()
{@; static id_type name=main_hash_table->match_literal("=");
  return name;
}

@* Function calls.
%
One of the most basic and important tasks of the evaluator is to allow
function calls, which may involve either built-in or user-defined functions.
This central part of the evaluator will be presented with an initial focus on
built-in functions, while leaving the particulars of user defined functions
(also known as $\lambda$-expressions) and their evaluation aside until
somewhat later. This corresponds more or less to the development history of
the interpreter, in which initially only built-in functions were catered for;
however many of the aspects that we deal with right away, notably function
overloading, are in fact much more recent additions than used-defined
functions were.

There will be several classes of expressions to represent function calls,
differing in the degree to which the called function has been identified
during type analysis. An intermediate class |call_base| between
|expression_base| and these classes is derived, to group some functionality
common to them. All call expressions take a general argument expression, and
location information is stored to allow the calling expression to be
identified during an error trace-back. Apart from the location information, an
error trace will also provide a name of the called function (which is more
readable than trying to reproduce the whole function call expression), which
will be obtained from the virtual method |function_name|.

@< Type def... @>=
struct call_base : public expression_base
{ expression_ptr argument;
  source_location loc;
@)
  call_base(expression_ptr&& arg, const source_location& loc)
  : argument(arg.release()), loc(loc) @+{}
  virtual ~@[call_base() nothing_new_here@];
  virtual std::string function_name() const=0;
};

@ We start with introducing a type for representing general function calls
after type checking. This is the general form where function can be given by
any kind of expression, not necessarily an applied identifier; indeed most
cases where a named function is called will handled by another kind of
expression, the overloaded call. In contrast with that, this type of call will
dynamically evaluate the function part (as opposed to argument) of the call,
possibly resulting in different functions between evaluations.

@< Type def... @>=
struct call_expression : public call_base
{ expression_ptr function;
@)
  call_expression
    (expression_ptr&& f,expression_ptr&& a, const source_location& loc)
   : call_base(std::move(a),loc), function(f.release()) @+{}
  virtual ~@[call_expression() nothing_new_here@];
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
  virtual std::string function_name() const
    // here we just print-wrap the function expression
    {@; std::ostringstream o; o << *function; return o.str(); }

};

@ To print a function call we print the function expression, enclosed in
parentheses unless it is an identifier, and the argument, enclosed in
parentheses unless it is a tuple expression (which already has parentheses).
The conditions for suppressing parentheses are tested using dynamic casts.

@< Function definitions @>=
void call_expression::print(std::ostream& out) const
{ if (dynamic_cast<identifier*>(function.get())!=nullptr)
    out << *function;
  else out << '(' << *function << ')';
  if (dynamic_cast<tuple_expression*>(argument.get())!=nullptr)
    out << *argument;
  else out << '(' << *argument << ')';
}

@ When a call involves a built-in function, what is executed is a value of
type |wrapper_function|, which is a |typedef| for a specific kind of function
pointer, defined in \.{global.h}.

@< Includes needed in the header file @>=

#include "global.h" // for |wrapper_function|

@ The class of dynamic values holding a wrapper function is called
|builtin_value|. Besides the function pointer it also stores a print name,
which is used when the wrapper function, rather than being called, gets
printed as (part of) a value in its own right; it is also used when reporting
an error during the execution of the built-in function. Most |builtin_value|
instances are constructed at start-up time when functions are entered into the
global overload table; their |print_name| will stick, even if the user binds
it to a new name.

@< Type definitions @>=

struct builtin_value : public value_base
{ wrapper_function val;
  std::string print_name;
@)
  builtin_value(wrapper_function v,const std::string& n)
  : val(v), print_name(n) @+ {}
  virtual ~ @[builtin_value() nothing_new_here@];
  virtual void print(std::ostream& out) const
  @+{@; out << '{' << print_name << '}'; }
  virtual builtin_value* clone() const @+{@; return new builtin_value(*this); }
  static const char* name() @+{@; return "built-in function"; }
private:
  builtin_value@[(const builtin_value& v) = default@];
};
typedef std::shared_ptr<const builtin_value> shared_builtin;

@ While syntactically more complicated than ordinary function calls, the call
of overloaded functions is actually simpler at run time, because the function
is necessarily referred to by an identifier (or operator) instead of by an
arbitrary expression, and overloading resolution results in a
function \emph{value} that has been identified at analysis time. If that value
happens to be a built-in function, the call will be translated into an
|overloaded_builtin_call| rather than into a |call_expression| (otherwise the
call will become an |overloaded_closure_call| that will be defined below).
Here we store a shared pointer to the |builtin_value|, which has the advantage
of not duplicating the |print_name| string for every call expression. To avoid
that this costs an extra pointer dereference at each call, we copy the
function pointer directly into |overloaded_builtin_call| as its field~|f|.

@< Type definitions @>=
struct overloaded_builtin_call : public call_base
{ wrapper_function f; // shortcut to implementing function
  shared_builtin fun; // points to the full |builtin_value|
@)
  overloaded_builtin_call
    (const shared_builtin& fun,expression_ptr&& a,const source_location& loc)
  : call_base(std::move(a),loc), f(fun->val), fun(fun) @+ {}
  virtual ~@[overloaded_builtin_call() nothing_new_here@];
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
  virtual std::string function_name() const @+{@; return fun->print_name; }
};

@ When printing, we use the |fun| field for its |print_name|, which the method
|function_name| achieves; we ensure it is called non-virtually to avoid the
overhead (in fact no derived class redefines the method, but the compiler
cannot know that). For the argument list we proceed as for general function
calls.

@< Function definitions @>=
void overloaded_builtin_call::print(std::ostream& out) const
{ out << overloaded_builtin_call::function_name();
  if (dynamic_cast<tuple_expression*>(argument.get())!=nullptr)
    out << *argument;
  else out << '(' << *argument << ')';
}

@ Some built-in functions like |print| accept arguments of any types, and in
particular tuples of any length. For such functions there is no use in
adopting the approach used for other built-in functions of expanding argument
tuples on the stack; instead the argument is always considered as one value.
Fortunately such functions are necessarily accessed through overloading, so we
detect the fact that they are being used at analysis time. This fact is then
recorded it in the type of call expression generated, and the |evaluate|
method will ask for an unexpanded argument on the execution stack. Therefore
we derive a type from |overloaded_builtin_call| that will override only the
|evaluate| method.

@< Type definitions @>=
struct variadic_builtin_call : public overloaded_builtin_call
{ typedef overloaded_builtin_call base;
@)
  variadic_builtin_call @| (const shared_builtin& fun,expression_ptr&& a,
    const source_location& loc)
  : base(fun,std::move(a),loc)@+ {}
  virtual void evaluate(level l) const;
};

@ The following function builds a |variadic_builtin_call|, ensuring that the
argument is wrapped up in a |voiding| if |needs_voiding| holds. The Boolean
will be set to hold in the (rare) case that the argument type is |void| but
the argument expression is not just ``()''.

@< Local fun... @>=
expression_ptr make_variadic_call
  (const shared_builtin& fun
  ,expression_ptr&& a, bool needs_voiding
  ,const source_location& loc)
{ if (needs_voiding)
    a.reset(new voiding(std::move(a))); // wrap argument in voiding
  return expression_ptr(new variadic_builtin_call(fun,std::move(a),loc));
}

@*1 Evaluating calls of built-in functions.
%
We now discuss how at run time built-in functions are called. Basically the
task consists of evaluating the arguments, placing them on the
|execution_stack|, and then calling the wrapper function through its pointer.

Some complication is added to this in order to be able to provide a back-trace
in case an error occurs during the function call. In debugging mode, that is
when |verbosity>0|, the arguments with which the function was called will be
printed in the back trace; since these arguments no longer need to exist at
the time the error occurs, we need to anticipate this possibility, and we do
so by recording the argument(s) as a string in a local variable of the
evaluation method.

Of the three classes derived from |call_base|, the one with the simplest
|evaluate| method is the class |variadic_builtin_call|, as it treats whatever
arguments it receives as a single value, that can be converted to a string
simply by performing output to an |ostringstream|.

In order to provide a trace of interrupted functions in case of an error,
function calls are executed in a |try| block (this will be true as well for
cases to be given later). We have made sure that the evaluation of the
arguments(s) of the function were done outside this |try| block, since
reporting functions that have not yet started executing would be confusing. We
detach the code for the |catch| block, so that it can be textually shared
with another |evaluate| method.

@< Function definitions @>=
void variadic_builtin_call::evaluate(level l) const
{ std::string arg_string;
  argument->eval(); ; // always evaluate to single value
  if (verbosity>0) // then record argument(s) as string
  {@; std::ostringstream o;
    o << *execution_stack.back();
    arg_string = o.str();
  }
@)
  try
  {@; (*f)(l); } // call the built-in function
  @< Catch-block for exceptions thrown within function calls @>
}

@ To provide back-trace, we catch and re-throw an error after extending the
stored error string. The result is a list of interrupted named function calls,
from inner to outer.

The work of modifying the error string is common to several such |catch|
blocks, and relegated to a function |extend_message| to be defined presently.
The error string is modified withing the existing error object; this is a
possibility that error objects derived from our |error_base| provide, contrary
to standard error objects like |std::runtime_error|. Nonetheless, we need to
deal with some errors derived from |std::exception| but not from our
|error_base|; notably the Atlas library may throw |std::runtime_error| rather
than our (\.{axis}) |runtime_error|, and |std::bad_alloc| can be thrown from
many places. Those errors do not have a modifiable message field (and
|std::bad_alloc| cannot even be raised with a provided error string at all),
so we re-brand those exceptions as |runtime_error|, by throwing the latter
after initialising its message field from |e.what()| and extending it through
a call of~|extend_message|.

@:Catch to trace back calls@>

@< Catch-block for exceptions thrown within function calls @>=
catch (error_base& e)
{@; extend_message(e,this,fun,arg_string);
  throw;
}
catch (const std::exception& e)
{ runtime_error new_error(e.what());
  extend_message(new_error,this,fun,arg_string);
  throw new_error;
}

@ The function |extend_message| facilitates appending information to error
messages in |catch| blocks. It is called with, apart from the error~|e| whose
message is to be modified, the expression~|call| whose evaluation was
interrupted by the error, the value |fun| of the function called (either a
|builtin_value| or a |closure_value|), and a string~|arg| that in debug mode
describes the arguments (when not in debug mode the string will be empty and
is ignored).

We report the source location of the call expression and the name of the
function called (both obtained from |call|), and a source location for the
definition of the called function (obtained from |fun|) in case it is
user-defined; when the called function was built in we just report that.

@< Local fun... @>=
void extend_message
  (error_base& e,const call_base* call, const shared_value& fun,
   const std::string& arg)
{ std::ostringstream o;
  o << "\n(in call " << call->loc << " of " << call->function_name() << ", ";
  auto f=dynamic_cast<const closure_value*>(fun.get());
  if (f==nullptr)
     o << "built-in";
  else o << "defined " << f->p->loc;
  o << ')';
  if (verbosity>0)
    o << "\n  argument" << (arg[0]=='(' ? "s: " : ": ") << arg;
  e.message.append(o.str());
}

@ The |evaluate| method for ordinary built-in functions is similar to that of
generic functions, but is somewhat complicated by the fact that for efficiency
reasons arguments are directly evaluated onto the |execution_stack|, without
ever constructing a tuple for them (this is achieved by using the |multi_eval|
method). This somewhat complicates the code for recording the argument as a
string in debugging mode, since it needs to find out how many separate
arguments have been evaluated. This is done by recording the stack pointer
before arguments are evaluated, and comparing with its value afterwards.

Currently there are no built-in function that take no arguments, but the code
below caters for the possibility anyway. Recording the arguments as strings
happens for all function calls (we cannot predict which ones will throw an
error) and obviously has a performance penalty; this is the reason why this
work is only done in debug mode. By using the same variable name
|arg_string|, we can reuse the |catch| block defined before.

@< Function definitions @>=
void overloaded_builtin_call::evaluate(level l) const
{ std::string arg_string;
  if (verbosity==0)
    argument->multi_eval();
  else // record argument(s) as string
  { auto sp = execution_stack.size();
    argument->multi_eval(); // mark stack before evaluation
    std::ostringstream o;
    if (execution_stack.size()>sp+1) // multiple arguments
      for (o << '(';
           sp<execution_stack.size();
           o << (sp<execution_stack.size() ? ',' : ')')
          )
        o << *execution_stack[sp++];
    else if (execution_stack.size()==sp) o << "()"; // no arguments
    else
      o << *execution_stack.back(); // single argument case
    arg_string = o.str();
  }
@)
  try
  {@; (*f)(l); } // call the built-in function
  @< Catch-block for exceptions thrown within function calls @>
}


@ Finally we consider the case where evaluating a |call_expression| results in
calling a built-in function. Since the function to be called is here produced
by evaluating an expression (maybe as simple as an identifier), the fact that
it is a built-in rather than user-defined function can here only be
determined at run time. The part of this method that deals with the case of a
user defined function is split off, and will be presented later once we have
discussed the representation of user defined functions.

To evaluate a |call_expression| object, in which the function part can be any
expression, we must evaluate this function part, and then dynamically test
whether it is a built-in or a user-defined function. In the former case we
evaluate the arguments, expanding them on the |execution_stack|, and then call
the built-in function. In that case we pass on the |level| parameter that was
passed to |call_expression::evaluate| method to the built-in function, so that
if necessary it can in its turn return and expanded result (or no result at
all). The evaluation of user-defined functions will be detailed later, but we
can already say that in this case it will be more useful to receive the
argument on the stack as a single value.

We reuse the previous |catch| block literally a third time; this time not only
do we judiciously choose the name |arg_string| to match what we did before,
but also the local variable name |fun| to math the field name
|overloaded_builtin_call::fun| that the cited module referred to in previous
instances.

@< Function definitions @>=
void call_expression::evaluate(level l) const
{ function->eval(); @+ shared_value fun=pop_value();
  auto f = dynamic_cast<const builtin_value*>(fun.get());
  const bool user_defined = f==nullptr;
  std::string arg_string;
  if (verbosity==0)
    argument->evaluate(user_defined ? single_value : multi_value);
  else
  { auto sp = execution_stack.size();
    argument->evaluate(user_defined ? single_value : multi_value);
    std::ostringstream o;
    if (execution_stack.size()>sp+1) // multiple arguments
      for (o << '(';
           sp<execution_stack.size();
           o << (sp<execution_stack.size() ? ',' : ')')
          )
        o << *execution_stack[sp++];
    else if (execution_stack.size()==sp) o << "()"; // no arguments
    else
      o << *execution_stack.back(); // single argument case
    arg_string = o.str();
  }
@)
  try
  { if (user_defined)
      @< Call user-defined function |fun| with argument on |execution_stack| @>
    else // built-in functions
      (*f->val)(l); // call the wrapper function, handling |l| appropriately
  }
  @< Catch-block for exceptions thrown within function calls @>
}


@*1 Type-checking function calls.
%
We now discuss the treatment of function calls at the time of type analysis,
and how the instances of classes derived from |call_base| come to be.

When we type-check a function call, we must expect the function part to be any
type of expression. But when it is a single identifier (possibly an operator
symbol) for which one or more overloads are defined then we attempt overload
resolution, unless the same identifier is locally bound with function type as
such bindings take precedence (however we ignore a possible binding for the
identifier in the global identifier table, even if it should have function
type). In all other cases (including that of a local function identifier), the
known type of the function expression gives the argument and result
types, and can be used to help converting the argument expression and the call
expression itself. Thus in such cases we first get the type of the expression
in the function position, requiring only that it be a function type, then
type-check and convert the argument expression using the obtained result type,
and build a converted function call~|call|. Finally (and this is done by
|conform_types|) we test if the required type matches the return type (in
which case we simply return~|call|), or if the return type can be coerced to
it (in which case we return |call| as transformed by |coerce|); if neither is
possible |conform_types| will throw a~|type_error|.

@:type-check calls@>

@< Cases for type-checking and converting... @>=
case function_call:
{ const application_node& call = *e.call_variant;
  if (call.fun.kind==applied_identifier)
    @< Convert and |return| an overloaded function call
    if |call.fun| is known in |global_overload_table|,
    unless it is a local function identifier @>
  type_expr f_type=gen_func_type.copy(); // start with generic function type
  expression_ptr fun = convert_expr(call.fun,f_type);
  expression_ptr arg = convert_expr(call.arg,f_type.func->arg_type);
  if (f_type.func->arg_type==void_type and not is_empty(call.arg))
    arg.reset(new voiding(std::move(arg)));
  expression_ptr result(new @|
     call_expression(std::move(fun),std::move(arg),e.loc));
  return conform_types(f_type.func->result_type,type,std::move(result),e);
}

@ When a call expression has an identifier in the place of the function (as is
often the case; this also includes all operators applied in formulae), we
prepare a call to the function |resolve_overload| defined above in
section@#resolve_overload@>. Before doing that, we check if the identifier has
a local binding with function type, in which we fall through the code below to
make a |call_expression| as in the general case, without any overloading. The
call is also omitted when the identifier is absent from the overload table
altogether; in that case it might still be a global identifier with function
type.

The cases relegated to |resolve_overload| include, due to call to
|is_special_operator| below, calls of special operators like the size-of
operator~`\#', even if such an operator should be absent from the overload
table.

@< Convert and |return| an overloaded function call... @>=
{ const id_type id =call.fun.identifier_variant;
  size_t i,j; bool b; // dummies; local binding not used here
  auto local_type_p=layer::lookup(id,i,j,b);
  if (local_type_p==nullptr or local_type_p->kind!=function_type)
 // not calling by local identifier
  { const overload_table::variant_list& variants
      = global_overload_table->variants(id);
    if (variants.size()>0 or is_special_operator(id))
    @/return resolve_overload(e,type,variants);
  }
}

@ Here is how we match special operators with generic argument type patterns;
they have an identifier~|id| that satisfies the predicate
|is_special_operator| to be defined below, and which is tested to ensure
such function calls get to the current code even if the operator should be
absent from the overload table. The built-in function objects that are
inserted into the calls here do not come from the |global_overload_table|, but
from a collection of static variables whose name ends with |_builtin|, and
which are initialised in a module given later, using calls to
|std::make_shared| so that they refer to unique shared instances.

The function |print| (but not |prints|) will return the value printed if
required, so it has the type of a generic identity function. This is done so
that inserting |print| around subexpressions for debugging purposes can be
done without other modifications of the user program. For this to work in all
cases, we treat the argument of |print| as if it we directly in the context of
the call to |print| (unless that is a void context, lest the argument get
voided before |print| sees it): if necessary, we re-convert the argument in a
|type| context (then coercion will be done inside the argument, and no
coercion applies directly to the |print| call itself). It is still
theoretically possible that inserting a call to print into valid code results
in an error, namely for argument expressions that fail to produce
an \foreign{a priori} type at all in the initial conversion (done without the
|type|context); in such cases an error will have been reported even before we
even get to the code below. These cases are quite rare though, and can be
overcome by inserting a cast inside the |print|.

In the case of |prints|, the context must either expect or accept a |void|
type, which is the condition that the call |type.specialise(void_type)| below
tests. For |to_string| the context must similarly either expect or accept a
|string| type. The case of |error| is like |prints| for its arguments, but
will not return, so nothing at all is demanded of the context type, like in
the case of \&{die}.

@< If |id| is a special operator like \#... @>=
{ if (id==size_of_name())
  { if (a_priori_type.kind==row_type)
    { expression_ptr call(new @|
        overloaded_builtin_call(sizeof_row_builtin,std::move(arg),e.loc));
      return conform_types(int_type,type,std::move(call),e);
    }
    else if (is_pair_type(a_priori_type))
    @< Recognise and return 2-argument versions of `\#', or fall through @>
  }
  else if (id==concatenate_name())
    @< Recognise and return instances of `\#\#', or fall through @>
  else // remaining cases always match
  { const bool needs_voiding = a_priori_type==void_type and not is_empty(args);
    if (id==print_name())
    { if (type!=void_type and not type.specialise(a_priori_type))
        arg = convert_expr(args,type); // redo conversion, now in |type| context
      return
        make_variadic_call(print_builtin,std::move(arg),needs_voiding,e.loc);
    }
    else if(id==to_string_name()) // this always matches as well
    {
      expression_ptr call =
        make_variadic_call(to_string_builtin,std::move(arg),needs_voiding,e.loc);
      return conform_types(str_type,type,std::move(call),e);
    }
    else if(id==prints_name()) // this always matches as well
    { expression_ptr call =
        make_variadic_call(prints_builtin,std::move(arg),needs_voiding,e.loc);
      return conform_types(void_type,type,std::move(call),e);
      // check that |type==void_type|
    }
    else if(id==error_name()) // this always matches as well
      return
        make_variadic_call(error_builtin,std::move(arg),needs_voiding,e.loc);
  }
}

@ We used the following simple type predicate above. Contrary to specialising
to |pair_type|, this function cannot alter its argument.

@< Local function definitions @>=
bool is_pair_type(const type_expr& t)
@+{@; return t.kind==tuple_type and length(t.tupple)==2; }

@ The operator `\#' can also be used as infix operator, to extend row value on
either end by a single element. The corresponding wrapper functions will be
defined in
%
section@#hash wrappers@>.
%
We require that one operand has row-of type, with component type equal to the
type of the other operand (we allow no coercion of operands, as for almost all
generic operators). We do however make an exception to the rule that operands
with type \.{[*]} (from the expression $[\,]$, or more often from identifiers
that have been initialised with that expression) will not match in overloading
unless that exact type is specified: for the operator `\#', if one of the
operands has this type, then its unknown component type will be specialised to
the type of the other operand.

These rules allow for a few ambiguous cases, when both operands have row type,
but these are quite rare. We resolve the ambiguity by stipulating that suffix
is preferred over prefix when both can apply. This implies that for instance
$[[2]]\#[\,]$ will give as result $[[2],[\,]]$, while $[\,]\#[[2]]$ will give
$[[[2]]]$.

@< Recognise and return 2-argument versions of `\#'... @>=
{
  const type_expr& ap_tp0 = a_priori_type.tupple->contents;
  const type_expr& ap_tp1 = a_priori_type.tupple->next->contents;
  if (ap_tp0.kind==row_type and
      ap_tp0.component_type->specialise(ap_tp1)) // suffix case
  { expression_ptr call(new @| overloaded_builtin_call
      (suffix_elt_builtin,std::move(arg),e.loc));
    return conform_types(ap_tp0,type,std::move(call),e);
  }
  if (ap_tp1.kind==row_type and
      ap_tp1.component_type->specialise(ap_tp0)) // prefix case
  { expression_ptr call(new @| overloaded_builtin_call
      (prefix_elt_builtin,std::move(arg),e.loc));
    return conform_types(ap_tp1,type,std::move(call),e);
  }
}

@ For `\#\#' we have instances concatenating all elements of a row of rows,
and another concatenating two rows of the same type.

@< Recognise and return instances of `\#\#'... @>=
{ if (a_priori_type.kind==row_type and
      a_priori_type.component_type->kind==row_type)
  { expression_ptr call(new @|
      overloaded_builtin_call(join_rows_row_builtin,std::move(arg),e.loc));
    return conform_types(*a_priori_type.component_type,type,std::move(call),e);
  }
  if (a_priori_type.kind==tuple_type and
    @|length(a_priori_type.tupple)==2 and
    @|a_priori_type.tupple->contents==a_priori_type.tupple->next->contents and
    @|a_priori_type.tupple->contents.kind==row_type)
  { expression_ptr call(new @| overloaded_builtin_call
        (join_rows_builtin,std::move(arg),e.loc));
    return conform_types
         (a_priori_type.tupple->contents,type,std::move(call),e);
  }
}

@* Let-expressions, and identifier patterns.
%
We shall now consider a simple type of expression in which local variables
occur, the let-expression. It is equivalent to an anonymous function
($\lambda$-expression) applied to the initialiser expression(s) in the
let-declarations, but no types need to be declared for the variable(s)
introduced, since the types of the corresponding expressions can be used for
this. Nevertheless, let-expressions could be represented internally, and
indeed were for a very long time, as a call in which the function is given by
a $\lambda$-expression, thus hiding the syntactic origin of the expression.
Currently they instead use a separate internal representation as
|let_expression|, whose evaluation is somewhat more efficient than that of the
equivalent call form. In our presentation we shall discuss let-expressions
before $\lambda$-expressions, which given an opportunity to compare them.

@ We prepare the definition of let-expressions with the introduction of
auxiliary types, needed to deal with the general patterns by which new
variables (and later formal function parameters) can be introduced. The parser
produces values of type |id_pat| to describe such patterns, and they are
needed at runtime to guide the decomposition of the corresponding values (in
fact the names used in the pattern are not needed for this, but they are
useful for printing the let-expression). The |id_pat| structure can be
used directly in a let-expression, but ownership must be properly
handled. The value produced by the parser will be destroyed after the command
is processed, at which time the let-expression could still exist (if held in a
function body), so we cannot simply use a non-owned pointer.

The top-level |id_pat| structure for the bound variable(s) will be stored in
the |let_expression| itself. We might move it there from the |id_pat| produced
by the parser and thereby steal a possible |sublist| field and further nodes
accessible from it. However this would amputate that expression (which would
be weird if the expression were printed in an error message), so instead doing
a deep copy is a cleaner solution, and patterns are rarely very deep. The
function |copy_id_pat| accomplishes making the copy. The implicit recursion of
|copy_id_pat| is achieved by passing itself as final argument to
|std::transform|; recursion terminates when |sublist.empty()| holds. This is
certainly a relative high on the coolness scale in implicit \Cpp\ programming.

@h <algorithm>
@< Local function def... @>=
id_pat copy_id_pat(const id_pat& p)
{
  containers::sl_list<id_pat> rl;
  std::transform(p.sublist.begin(),end(p.sublist)
                , std::back_inserter(rl), copy_id_pat);
  return id_pat (p.name, p.kind, rl.undress());
}

@ Now we can define our let-expression to use a |shared_pattern|. Initialiser
and body are owned, and held in unique-pointers.

@< Type def... @>=
struct let_expression : public expression_base
{ const id_pat variable;
  expression_ptr initialiser, body;
@)
  let_expression(const id_pat& v, expression_ptr&& ini, expression_ptr&& b);
  virtual ~@[let_expression() nothing_new_here@]; // subobjects do all the work
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ The main constructor cannot be inside the class definition, as it requires
the local function |copy_id_pat|. The variable (pattern) is copied into the
|let_expression| node, while for the two expression components ownership is
transferred from the passed unique-pointer values.

@< Function def... @>=
inline
let_expression::let_expression @|
  (const id_pat& v, expression_ptr&& ini, expression_ptr&& b)
: variable(copy_id_pat(v))
, initialiser(ini.release())
, body(b.release())
@+{}

@ To print a let expression, we reproduce the input syntax, as far as
possible; restructuring done during parsing makes it impossible to know the
exact form used at input.

@< Function definitions @>=
void let_expression::print(std::ostream& out) const
{@; out << " let " << variable << "=" << *initialiser << " in " << *body;
}

@ We shall need some other functions to deal with patterns, all with a
similar structure.

@s back_insert_iterator vector

@< Declarations of exported functions @>=
type_expr pattern_type(const id_pat& pat);
size_t count_identifiers(const id_pat& pat);
void list_identifiers(const id_pat& pat, std::vector<id_type>& d);
void thread_bindings
  (const id_pat& pat,const type_expr& type, layer& dst, bool is_const);
void thread_components
  (const id_pat& pat,const shared_value& val,
   std::back_insert_iterator<std::vector<shared_value> > dst);

@ For handling declarations with patterns as left hand side, we need a
corresponding type pattern; for instance $(x,,(f,):z)$:\\{whole} requires the
type \.{(*,*,(*,*))}. These recursive functions construct such types.

@< Function definitions @>=
type_list pattern_list_types(const patlist& p)
{ dressed_type_list result;
  for (auto it = p.begin(); not p.at_end(it); ++it)
    result.push_back(pattern_type(*it));
  return result.undress();
}
@)
type_expr pattern_type(const id_pat& pat)
{@; return (pat.kind&0x2)==0
  ? unknown_type.copy()
  : type_expr(pattern_list_types(pat.sublist));
}

@ Here we count the number or list the identifiers in a pattern. The latter
function uses an output parameter rather than a return value since this avoids
doing any concatenation of vectors. Instead of a modifiable reference~|d| to a
vector it could have used an output iterator, but in practice we always
collect the results in a vector, and this avoids having to call
|std::back_inserter| all the time.

@< Function definitions @>=
size_t count_identifiers(const id_pat& pat)
{ size_t result= pat.kind & 0x1; // 1 if |pat.name| is defined, 0 otherwise
  if ((pat.kind & 0x2)!=0) // then a list of subpatterns is present
    for (auto it=pat.sublist.begin(); not pat.sublist.at_end(it); ++it)
      result+=count_identifiers(*it);
  return result;
}
@)
void list_identifiers(const id_pat& pat, std::vector<id_type>& d)
{ if ((pat.kind & 0x1)!=0)
    d.push_back(pat.name);
  if ((pat.kind & 0x2)!=0) // then a list of subpatterns is present
    for (auto it=pat.sublist.begin(); not pat.sublist.at_end(it); ++it)
      list_identifiers(*it,d);
}

@ Here we do a similar traversal, using a type with structure matching |pat|;
we push pairs onto a |layer|. If |is_const| holds, all identifiers will be
const; otherwise a bit from |pat.kind| determines constness of
|pat.name|.

@< Function definitions @>=
void thread_bindings
(const id_pat& pat,const type_expr& type, layer& dst, bool is_const)
{ if ((pat.kind & 0x1)!=0)
    dst.add(pat.name,type.copy(),is_const or (pat.kind & 0x4)!=0);
  if ((pat.kind & 0x2)!=0)
  { assert(type.kind==tuple_type);
    wtl_const_iterator t_it(type.tupple);
    for (auto p_it=pat.sublist.begin(); not pat.sublist.at_end(p_it);
         ++p_it,++t_it)
      thread_bindings(*p_it,*t_it,dst,is_const);
  }
}

@ Finally, there is similar processing of an appropriate |shared_value| at
runtime. This time we do use an output iterator. We could have written |*dst++
= val| below simply as |dst=val|, omitting operators that just return their
arguments; however the given expression is more in the spirit of general
iterator handling.

@< Function definitions @>=
void thread_components
  (const id_pat& pat,const shared_value& val
  , std::back_insert_iterator<std::vector<shared_value> > dst)
{ if ((pat.kind & 0x1)!=0)
     *dst++ = val; // copy |shared_value| pointer |val|, creating sharing

  if ((pat.kind & 0x2)!=0)
  { const tuple_value* t=force<tuple_value>(val.get());
    assert(t->val.size()==length(pat.sublist));
    auto src = t->val.begin();
    for (auto it=pat.sublist.begin(); not pat.sublist.at_end(it); ++it,++src)
      thread_components(*it,*src,dst);
  }
}

@ To convert a let-expression, we first deduce the type of the declared
identifiers from the right hand side of its declaration, then set up new
bindings for those identifiers with the type found, and finally convert the
body to the required type in the extended context. Note that the constructed
|layer| is a local variable whose constructor pushes it onto the
|layer::lexical_context| list, and whose destructor will pop it off.

If there are no identifiers at all, we should avoid that the execution of the
|let_expression| push an empty frame on the evaluation context, as this would
wreak havoc due to the fact that we made applied identifiers not count empty
layers. Rather than have the |let_expression::evaluate| method handle
this dynamically, we avoid generating such a |let_expression| altogether, and
instead generate a |seq_expression| whose |evaluate| method does exactly what
is needed in this case.

@< Cases for type-checking and converting... @>=
case let_expr:
{ const auto& lexp=*e.let_variant;
  const id_pat& pat=lexp.pattern;
  type_expr decl_type=pattern_type(pat);
  expression_ptr arg = convert_expr(lexp.val,decl_type);
@/auto n=count_identifiers(pat);
  if (n==0)
  // then avoid frame without identifiers, so compile as sequence expression
    return expression_ptr(new @|
      seq_expression(std::move(arg),convert_expr(lexp.body,type)));
  if (decl_type==void_type and not is_empty(lexp.val))
    // rare case, introducing void identifier
    arg.reset(new voiding(std::move(arg)));
  layer new_layer(n);
  thread_bindings(pat,decl_type,new_layer,false);
  return expression_ptr(new @|
    let_expression(pat,std::move(arg),convert_expr(lexp.body,type)));
}

@ Here is a class whose main purpose, like that of |layer| before, is to have
a constructor-destructor pair, in this case one that temporarily suspends the
current execution context, replacing it by a new one determined by an
identifier pattern and an execution context for the enclosing lexical layers.
All instances of this class should be automatic (local) variables, to ensure
that they have nested lifetimes.

We take care when pushing a new |evaluation_context| to avoid changing the
reference count of the pointer in |frame::current| as would happen when
copying it, by \emph{moving} the pointer to the tail of the new node. When
popping on destruction however we need to copy, since the node being popped
could have become accessed independently (through a |closure|, to be discussed
later), so we are not free to move from the tail of this node.

@:frame class@>

@< Local class definitions @>=
class frame
{
  const id_pat& pattern;
public:
  static shared_context current;
@)
  frame (const id_pat& pattern)
  : pattern(pattern)
  { assert(count_identifiers(pattern)>0); // avoid frames without identifiers
    current = std::make_shared<evaluation_context>(std::move(current));
  }
  ~frame() @+{@; current = current->tail(); } // don't use |std::move| here!
@)
  void bind (const shared_value& val)
    { current->reserve(count_identifiers(pattern));
      thread_components(pattern,val,current->back_inserter());
    }
  std::vector<id_type> id_list() const; // list identifiers, for back-tracing
};

@ This method is only called during exception handling, so a simple access to
the identifiers is more important than an efficient one. Therefore we convert
the pattern to a vector, using a  call to |list_identifiers|.

@< Local function definitions @>=
std::vector<id_type> frame::id_list() const
{ std::vector<id_type> names; names.reserve(count_identifiers(pattern));
  list_identifiers(pattern,names);
  return names;
}

@ Evaluating a let expression is now straightforward: evaluate the initialiser
to produce a value on the stack; then create a new |frame| in which this value
is bound to the pattern |variable|, and in this extended context evaluate
the~|body|.

The stack will be automatically popped when the lifetime of the |frame| ends.
Usually this happens when we return from the |let_expression::evaluate|, but
it also happens in case an error is thrown. The latter also includes instances
of |break| or |return| premature exits from a loop or a user-defined function
in an axis program, which are implemented by throwing and catching specific
exceptions. In all cases, the static variable |frame::current| will
automatically be set to the proper value by the |~frame| destructor.

@< Function def... @>=
void let_expression::evaluate(level l) const
{ initialiser->eval(); // evaluate on stack as single value
@)
  frame fr(variable); // save context, create new one for |f|
  fr.bind(pop_value()); // decompose arguments(s) and bind values in |fr|
  try {@; body->evaluate(l); }
    // call, passing evaluation level |l| to function body
  @< Catch block for providing a trace-back of local variables @>
} // restore context upon destruction of |fr|


@ Providing the current values of local variables at an error stop is useful
and quite easy. We obtain the names of the identifiers from the frame |fr| to
which they are copied rather than directly from our |variable|, with the main
purpose of reusing the module identically as |catch| block for user defined
functions, where the current (call) expression does not have an identifier
pattern available, but there is a frame |fr| from which is can be obtained.

@< Catch block for providing a trace-back of local variables @>=
catch (error_base& e)
{ std::vector<id_type> names = fr.id_list();
  auto id_it = names.cbegin(); std::ostringstream o; o << "\n  [";
  for (auto it = frame::current->begin(); it!=frame::current->end();
       ++it,++id_it)
  { if (it!=frame::current->begin())
      o << ", ";
    o << main_hash_table->name_of(*id_it) << '=' << **it;
  }
  o << ']';
  e.message.append(o.str());
  throw;
}


@*1 Lambda-expressions (user-defined functions).
%
In contrast to let-expressions, a $\lambda$-expression can be evaluated one
or more times, yielding ``closure'' values that need to refer to the pattern,
and which might outlive the $\lambda$-expression, so appropriate duplication
or sharing must be organised. We opt for sharing between the
$\lambda$-expression and any closures obtained from it. Since the evaluator
handles expressions by reference to |expression_base|, one cannot achieve
sharing directly to |lambda_expression|; instead we store a shared pointer to
a structure with the necessary components.

Since this is the kind
of runtime value that will hold the result of a user function definition, we
provide a field |loc| to record the source location.

@< Type def... @>=
struct lambda_struct
{ id_pat param; @+ expression_ptr body; @+ source_location loc;
  lambda_struct
      (id_pat&& param, expression_ptr&& body, const source_location& loc)
  : param(std::move(param)), body(std::move(body)), loc(loc) @+{}
};
typedef std::shared_ptr<lambda_struct> shared_lambda;
@)
struct lambda_expression : public expression_base
{ shared_lambda p;
  @)
  lambda_expression
    (const id_pat& p, expression_ptr&& b, const source_location& loc);
  virtual ~@[lambda_expression() nothing_new_here@];
    // subobjects do all the work
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};


@ The main constructor cannot be inside the class definition, as it requires
the local function |copy_id_pat|. It copies the pattern and creates a new
shared reference to the copy (further sharing will occur when the
$\lambda$-expression is evaluated). The |loc| field is copy-constructed from
the one passed, which resides in a |const|-qualified |expr| object produced by
the parser; therefore moving is not an option here, and since this is plain
data it wouldn't be more efficient anyway.

@< Function def... @>=
inline
lambda_expression::lambda_expression @|
  (const id_pat& p, expression_ptr&& b, const source_location& loc)
: p(std::make_shared<lambda_struct>(copy_id_pat(p),std::move(b),loc))
@+{}

@ To print an anonymous function, we print the parameter list, followed by a
colon and by the function body. If the parameter list contains a name for the
whole, as happens in particular when there is just a single parameter, then it
must be enclosed in parentheses to resemble in the input syntax, but if it is
an unnamed nonempty tuple then it will supply its own parentheses; this leaves
the case of an empty parameter list, for which we reconstitute the input
syntax~`\.@@'. The printed parameter list cannot include types with the
current setup, as they are not explicitly stored after type analysis. It could
be made possible to print types if a type were explicitly stored in the
|lambda_expression| structure; at the time of writing this would seem possible
because each function has to have a definite type, but if the type system were
extended with second order types (which would be quite useful), then this
might no longer be true. For now leaving out the types indicates to the
(knowledgeable) user that a runtime value is being printed rather than just a
syntax tree representing the user input (as happens in messages from the type
checker).

@< Function definitions @>=
std::ostream& operator<<(std::ostream& out, const lambda_struct& l)
{ if ((l.param.kind&0x1)!=0)
    if ((l.param.kind&0x4)!=0)
      out << "(!" << l.param << ')';
    else out << '(' << l.param << ')';
  else if ((l.param.kind&0x2)!=0 and not l.param.sublist.empty())
    out << l.param;
  else out << '@@';
  return out << ": " << *l.body;
}
void lambda_expression::print(std::ostream& out) const @+{@; out << *p; }

@ Handling of user-defined functions in type analysis is usually uneventful
(apart from the analysis of its body, which can of course be arbitrarily
complex), but still has to be done with some care to correctly handle rare
circumstances.

The main work here is to call |thread_bindings| to extract from the specified
argument pattern the identifiers it contains, couple them to their specified
types, and store the result in |new_layer|, to be used during the type-check
and conversion of the function body. It is rare for a $\lambda$-expression to
appear in a context that requires a specific type, but it can happen (the
right hand side of an assignment for instance) and necessitates various tests
in the code below. After conversion of the function body, no implicit
conversions are ever needed (as there are none that ever apply to expressions
of function type), so there is no point in calling |conform_types| at the end.
However it is possible that the context requires a void type; like for
denotations, this is silly (evaluation will have no side effect) but legal.
So we handle both cases here, sharing code where possible.

In non-void context we specialise the required |type| (often undetermined
initially) to a function type with argument type the one given in the
$\lambda$-expression (signalling a type error if a different type was
expected), then we convert the function body in the new context, specialising
the return type. In void context we do only the conversion (just for error
checking), and ignore the return type.

@< Cases for type-checking and converting... @>=
case lambda_expr:
{ const lambda_node& fun=*e.lambda_variant;
  const id_pat& pat=fun.pattern;
  const type_expr& arg_type=fun.parameter_type;
    // argument type specified in |fun|
  if (not pattern_type(pat).specialise(arg_type))
    // do |pat| structure and |arg_type| conflict?
    throw expr_error(e,"Function argument pattern does not match its type");
  type_expr* rt; type_expr dummy;
  if (type.specialise(gen_func_type)
             and type.func->arg_type.specialise(arg_type))
    rt = &type.func->result_type; // we can now safely access this
  else if (type==void_type)
    rt=&dummy; // in void context there is no return type to set
  else
    @/throw type_error(e,
                       type_expr(arg_type.copy(),unknown_type.copy()),
                       std::move(type));
@/layer new_layer(count_identifiers(pat),rt);
  thread_bindings(pat,arg_type,new_layer,false);
@/return expression_ptr(new @|
      lambda_expression(pat, convert_expr(fun.body,*rt), std::move(e.loc)));
}

@* Closures.
%
In first approximation $\lambda$-expression are like denotations of
user-defined functions: their evaluation just returns the stored function
body. However, this evaluation also captures the current execution context:
the bindings of the local variables that may occur as free identifiers in the
function body (any used global variables can be bound at compile time, so they
do not need any special consideration). Therefore the evaluation of a
$\lambda$~expressions actually yields an intermediate value that is
traditionally called a closure. It contains a shared pointer~|p| to the
|lambda_struct| holding the function body, as well as the execution context
|context| that was current at the point in time the $\lambda$~expression was
encountered.

Sharing the |lambda_struct| among different closures obtained from the same
$\lambda$-expression is efficient in terms of space, but would require double
dereference upon evaluation. Since the latter occurs frequently, we speed up
evaluation by also using a reference |body| directly to the function body.
Note that closures are formed when the \emph{definition} of a user-defined
function is processed, so this optimisation should make evaluation of globally
defined functions a bit faster. For local functions, the closure is formed
during the execution of the outer function, so the optimisation only helps if
the closure formed will be called more than once; this is still probable,
though there are usage patters (for instance simulating a case statement by
selecting a closure from an array, and then calling it) for which local
closures are actually executed less than once on average; in such cases we are
actually wasting effort here.

@< Type def... @>=
struct closure_value : public value_base
{ shared_context context;
  shared_lambda p;
  const expression_base& body; // shortcut to function body
@)
  closure_value@|(const shared_context& c, const shared_lambda& l)
  : context(c), p(l), body(*p->body) @+{}
  virtual ~ @[closure_value() nothing_new_here@];
  virtual void print(std::ostream& out) const;
  virtual closure_value* clone() const @+
  {@; return new closure_value(context,p); }
  static const char* name() @+{@; return "closure"; }
};
typedef std::unique_ptr<closure_value> closure_ptr;
typedef std::shared_ptr<const closure_value> shared_closure;

@ A closure prints the |lambda_expression| from which it was obtained, but we
also print an indication of where the function was defined (this was not
useful for |lambda_expression|, since these never get printed directly, only
as part of printing a |closure_value|). One could imagine printing after this
body ``where'' followed by the bindings held in the |c| field. Even better
only the bindings for relevant (because referenced) identifiers could be
printed. But it's not done yet.

@< Function def... @>=
void closure_value::print(std::ostream& out) const
{@; out << "Function defined " << p->loc << std::endl << *p;
}

@ Evaluating a $\lambda$-expression just forms a closure using the current
execution context, and returns that.

While this code looks rather innocent, note that the sharing of
|frame::current| created here may survive after one or more frames on the list
|frame::current| get removed after evaluation of the expression that returned
the closure; these frames then get an extended lifetime through the closure
formed here. This implies that the execution context cannot be embedded in any
kind of stack, in particular it cannot be embedded in the \Cpp\ runtime stack
(while the layers of the lexical context could).

@< Function def... @>=
void lambda_expression::evaluate(level l) const
{@;if (l!=no_value)
     push_value(std::make_shared<closure_value>(frame::current,p));
}

@*1 Calling user-defined functions.
%
In order to implement calling of user-defined functions, we define a variation
of the class |frame|. Again the purpose is to have a constructor-destructor
pair that temporarily suspends the current execution context, replacing it by
a new one determined by the parameter(s) of the $\lambda$-expression, on top
of the execution context stored in the closure. Again instances of this class
should be automatic variables, to ensure that they have nested lifetimes.

This context switching is a crucial and recurrent step in the evaluation
process, so we take care to not uselessly change the reference count of
|frame::current|. It is \emph{moved} into |saved| upon construction, and upon
destruction moved back again to |frame::current|. In contrast to |frame|, the
constructor here needs a try block for exception safety, as the call to
|std::make_shared<evaluation_context>| may throw an exception after
|frame::current| has been moved from, but before out constructor completes;
since the destructor would in this scenario \emph{not} be called, we then need
to move the pointer back explicitly in the |catch| block.

Since a |lambda_frame::bind| is called every time a user defined function is
called, this is a convenient point to check whether the signal handler has set
the interrupt flag, and to bail out if it did. (Choosing the points to do this
test is a somewhat delicate matter. One wants to do it at points that are
regularly encountered during evaluation, but not so frequently that the
checks incur a serious performance penalty. Just the test here does not
provide an absolute guarantee of rapid detection of a signalled interrupt.)

If one tried to derive this class from |frame|, one would have to construct
the base (which modifies |frame::current|) before doing anything else; this
would make saving the value of |frame::current| problematic. For that reason
it is better to just repeat some of the things done for |frame| independently,
similarly but with a few important changes. We do of course have to make use
of the static member |frame::current| of that class, which is the whole point
of defining the |lambda_frame| class.

@< Local class definitions @>=
class lambda_frame
{
  const id_pat& pattern;
  const shared_context saved;
  const bool empty;
public:
  lambda_frame (const id_pat& pattern, const shared_context& outer)
  : pattern(pattern)
  , saved(std::move(frame::current))
  , empty(count_identifiers(pattern)==0)
  { assert(&outer!=&frame::current); // for excluded case use |frame| instead
    try {@;
      frame::current =
        empty ? outer : std::make_shared<evaluation_context>(outer);
    }
    catch(...)
    {@; frame::current = std::move(saved); throw; }
    // restore as destructor would do
  }
  ~lambda_frame() @+{@; frame::current = std::move(saved); }
@)
  bool is_empty() const @+{@; return empty; }
  void bind (const shared_value& val)
    { assert(not empty);
      // one should not call |bind| for an |empty| lambda frame
      if (interrupt_flag!=0)
        throw user_interrupt();
      frame::current->reserve(count_identifiers(pattern));
      thread_components(pattern,val,frame::current->back_inserter());
    }
  std::vector<id_type> id_list() const; // list identifiers, for back-tracing
};

@ This method is identical to the one in |frame|, but as said, we cannot use
inheritance.

@< Local function definitions @>=
std::vector<id_type> lambda_frame::id_list() const
{ std::vector<id_type> names; names.reserve(count_identifiers(pattern));
  list_identifiers(pattern,names);
  return names;
}

@ In general a closure formed from a $\lambda$ expression can be handled in
various ways (like being passed as argument, returned, stored) before being
applied as a function, in which case the call is performed by
|call_expression::evaluate| described above; the actual code this executes is
given in section@# lambda evaluation @> below. However, in most cases the path
from definition to call is more direct: the closure from a user-defined
function is bound to an identifier (or operator) in the global overload table,
and located during type-checking of a call expression. As this special but
frequent case can be handled more efficiently than by building a
|call_expression|, we introduce a new |expression| type that is capable of
directly storing a closure value.

Closures themselves are anonymous, so the |print_name| reflects the overloaded
name that was used to identify this function; it can vary separately from the
closure |fun| if the latter is entered more than once in the the tables. This
is in contrast to |overloaded_builtin_call| where the name is taken from the
stored |builtin_value|, and cannot be dissociated from the wrapper function.

@< Type definitions @>=
struct overloaded_closure_call : public call_base
{ shared_closure fun;
  std::string print_name;
@)
  overloaded_closure_call @|
   (shared_closure f,const std::string& n,expression_ptr&& a
   ,const source_location& loc)
  : call_base(std::move(a),loc), fun(f), print_name(n) @+ {}
  virtual ~@[overloaded_closure_call() nothing_new_here@];
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
  virtual std::string function_name() const @+{@; return print_name; }
};

@ When printing it, we ignore the closure and use the overloaded function
name.

@< Function definitions @>=
void overloaded_closure_call::print(std::ostream& out) const
{ out << print_name;
  if (dynamic_cast<tuple_expression*>(argument.get())!=nullptr)
     out << *argument;
  else out << '(' << *argument << ')';
}

@ The definition of the function |resolve_overload| left one module to be
specified, namely building a function call once a match is found
in the overload table; we will give this part now. This code is in fact
included in the function twice, once for the case the \foreign{a priori} type
of the argument exactly matches the expected argument type, and once for an
inexact match where the argument had to be reconverted for the expected type;
however the difference between these scenarios has no effect on what is to be
done here.

Having found a match and converted the argument expression to |arg|, we need
to construct a function call object. We proceed in a similar fashion as in
section@#type-check calls@> when building a |call_expression|, but with the
argument already converted. There is no function part to convert either: the
overload table contains an already evaluated function value~|v|, which is
either a |builtin_value| representing a built-in function, or a
|closure_value| representing a user-defined function and its evaluation
context. Here we build an |overloaded_builtin_call| respectively an
|overloaded_closure_call| for those cases.

This is one of the places where we might have to insert a |voiding|, in the
rare case that a function with void argument type is called with a nonempty
argument expression. An alternative would be to replace the call by a sequence
expression, evaluating the argument expression separately and then the
function call with an empty argument expression. In any case the
test \emph{can} be made here, since we have type argument type in the variable
|arg_type|, and the argument expression in |args|.

As a special safety measure against the easily made error of writing `\.='
instead of an assignment operator~`\.{:=}', we forbid converting to void the
result of an (always overloaded) call to the equality operator, treating this
case as a type error instead. In the unlikely case that the user defines an
overloaded instance of `\.=' with void result type, calls to this operator
will still be accepted.

@< Set |call| to a call of the function value |*v.val|... @>=
{ if (arg_type==void_type and not is_empty(args))
    arg.reset(new voiding(std::move(arg)));
  if (@[auto f = std::dynamic_pointer_cast<const builtin_value>(v.val)@;@])
    call = expression_ptr (new @| overloaded_builtin_call
      (f,std::move(arg),e.loc));
  else
  if (@[auto fun = std::dynamic_pointer_cast<const closure_value>(v.val)@;@])
  { std::ostringstream name;
    name << main_hash_table->name_of(id) << '@@' << arg_type;
    call = expression_ptr(new
             overloaded_closure_call(fun,name.str(),std::move(arg),e.loc));
  }
  else
    throw logic_error("Overloaded value is not a function");
@)
  if (type==void_type and
      id==equals_name() and
      v.type().result_type!=void_type)
  { std::ostringstream o;
    o << "Use op equality operator '=' in void context; " @|
      << "did you mean ':=' instead?\n  If you really want " @|
      << "the result of '=' to be voided, use a cast to " @|
      << v.type().result_type << '.';
    throw expr_error(e,o.str());
      // importantly this is not a |type_error|, which might be caught
  }
}

@ When calling a function in a non overloading manner, we come to the code
below if it turns out not to be a built-in function. This code basically
implements a call-by-value $\lambda$-calculus evaluator, in concert with the
|evaluation_context| data structure and the translation of applied local
identifiers into references by relative layer number and offset. The names of
the identifiers are not used at all at runtime (they are present in the
|id_pat| structure for printing purposes, but ignored by the |bind| method
used here); our implementation can be classified as one using ``nameless
dummies'' (also known as ``de Bruijn indices'').

When we come here |f| must be a |closure_value|, the argument has already been
evaluated, and is available as a single value on the |execution_stack|. The
evaluation of the call temporarily replaces the current execution context
|frame::current| by one composed of |f->context| stored in the closure and a
new frame defined by the parameter list |f->param| and the argument obtained
as |pop_value()|; the function body is evaluated in this extended context.
Afterwards the original context is restored by the destructor of~|fr|, whether
the call completes normally or is terminated by a runtime error. The most
important advantage of this approach is in case of abnormal exits from loops
and functions, which are implemented by throwing and catching an exception at
run-time and will therefore unwind the \Cpp\ stack. (Actually, performing
|break| from a loop should never lead to destructing any |lambda_frame|,
though it might destruct some |frame|s.)

By naming our frame |fr|, we can textually reuse a |catch| block, as mentioned
at its definition.

@: lambda evaluation @>

@< Call user-defined function |fun| with argument on |execution_stack| @>=
try
{ const closure_value& f=*force<closure_value>(fun.get());
@)
  lambda_frame fr(f.p->param,f.context);
    // save context, create new one for |f|
  if (fr.is_empty()) // we must test for functions without named arguments
  {@; execution_stack.pop_back();
    f.body.evaluate(l);
  } //  drop argument, evaluate avoiding |bind|
  else
  { fr.bind(pop_value()); // decompose arguments(s) and bind values in |fr|
    try {@; f.body.evaluate(l); }
      // call, passing evaluation level |l| to function body
    @< Catch block for providing a trace-back of local variables @>
  }
} // restore context upon destruction of |fr|
@< Catch block for explicit |return| from functions @>

@ When the call |f.body.evaluate(l)| ends up executing an explicit |return|
expression, it won't have put anything on the |execution_stack|, but the value
to return will be stored inside the error object. (This is the right thing to
do in the unlikely case that intermediate values are removed from
|execution_stack| during the \Cpp\ stack unwinding.) It suffices to place the
value on the stack, which is just what |push_expanded| does.

@< Catch block for explicit |return| from functions @>=
catch (function_return& err) @+
{@; push_expanded(l,std::move(err.val)); }

@ Evaluation of an overloaded function call bound to a closure consists of a
simplified version of the part of |call_expression::evaluate| dedicated to
closures (including the temporary catching of errors as defined in
section@#Catch to trace back calls@> in order to produce a trace of
interrupted function calls in the error message) together with the evaluation
part of section@# lambda evaluation @> just given. The simplification
consists of the fact that the closure is already evaluated and stored, and
that in particular we don't have to distinguish dynamically between built-in
functions and closures.

Not having varied our naming conventions (|arg_string| and |fun|), we can make
a fourth textual reuse of the |catch| block for function calls, as well as
(|fr|) a third reuse of the |catch| block for local variables.

@< Function definitions @>=
void overloaded_closure_call::evaluate(level l) const
{ argument->eval(); // evaluate arguments as a single value
  std::string arg_string;
  if (verbosity!=0) // then record argument(s) as string
  {@; std::ostringstream o;
    o << *execution_stack.back();
    arg_string = o.str();
  }
  try
  { lambda_frame fr(fun->p->param,fun->context);
      // save context, create new one for |fun|
@)
    if (fr.is_empty()) // we must test for functions without named arguments
      {@; execution_stack.pop_back(); fun->body.evaluate(l); }
      //  avoid |bind|, evaluate
    else
    { fr.bind(pop_value()); // decompose arguments(s) and bind values in |fr|
      try {@; fun->body.evaluate(l); }
      // call, passing evaluation level |l| to function body
      @< Catch block for providing a trace-back of local variables @>
    }
  } // restore context upon destruction of |fr|
  @< Catch block for explicit |return| from functions @>
  @< Catch-block for exceptions thrown within function calls @>
}

@* Sequence expressions.
%
Sequence expressions are used to evaluate two expressions one after the other,
discarding any value from on of them (although it is possible that the other
value also gets discarded by voiding). In most cases it is the first
expression whose value is discarded (as for the comma-operator in \Cee/\Cpp),
and the semicolon is used to indicate this; occasionally however it is useful
to retain the first value and evaluate a second expression for its side
effects \emph{afterwards}, and the \.{next} keyword is used instead of a
semicolon for this purpose.

In practice it turns out that long sequences of expressions chained by
semicolons are rare, since often the need to introduce new local variables
will interrupt the chain. Therefore we see no need to use a |std::vector|
representation for such sequences of expressions, and prefer to use a chained
representation instead, just like before type analysis. In other words, we
use an expression node with two descendents each time. The forward and
reverse variants are implemented by similar but distinct types derived from
|expression_base|.

@< Type def... @>=
struct seq_expression : public expression_base
{ expression_ptr first,last;
@)
  seq_expression(expression_ptr&& f,expression_ptr&& l)
   : first(f.release()),last(l.release()) @+{}
  virtual ~@[seq_expression() nothing_new_here@];
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};
@)
struct next_expression : public expression_base
{ expression_ptr first,last;
@)
  next_expression(expression_ptr&& f,expression_ptr&& l)
   : first(f.release()),last(l.release()) @+{}
  virtual ~@[next_expression() nothing_new_here@];
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
retained as result. This is one occasion where a value pushed to the
|execution_stack| may remain there for some time (during all of
|last->void_eval();|); the other one is when |tuple_display::evaluate| is
called with a |l==multi_value| argument).

@< Function def... @>=
void next_expression::evaluate(level l) const
{@; first->evaluate(l); last->void_eval(); }

@ It remains to type-check and convert sequence expressions, which is easy.

@< Cases for type-checking and converting... @>=
case seq_expr:
{ const sequence_node& seq=*e.sequence_variant;
  expression_ptr first = convert_expr(seq.first,as_lvalue(void_type.copy()));
  expression_ptr last  = convert_expr(seq.last,type);
  return expression_ptr(new seq_expression(std::move(first),std::move(last)));
}
break;
case next_expr:
{ const sequence_node& seq=*e.sequence_variant;
  expression_ptr first = convert_expr(seq.first,type);
  expression_ptr last  = convert_expr(seq.last,as_lvalue(void_type.copy()));
  return expression_ptr(new next_expression(std::move(first),std::move(last)));
}
break;

@* Array subscription and slicing.
%
We have seen expressions to build lists, and although vectors and matrices can
be made out of them using coercions, we so far are not able to access their
components once they are constructed. To that end we shall now introduce
operations to index such values. We allow subscription of rows, but also of
vectors, rational vectors, matrices, strings, and of the Atlas \.{ParamPol}
values. Since after type analysis we know which of the cases applies for a
given expression, we define several classes among which type analysis will
choose. These classes differ mostly by their |evaluate| method, so we first
derive an intermediate class from |expression_base|, and derive the others
from it. This class also serves to host an enumeration type and some static
methods that will serve later. We include a case here, |mod_poly_term|, that
is related to a type defined in \.{atlas-types}.

At the same time we shall define ``slices'', which differ from subscriptions
in selecting a whole range of index values. The name is not really
appropriate, as it more evokes subscripting at certain index positions while
leaving other positions vary (an example would be selecting a row from
a matrix); we don't (yet?) support that as a form of slicing, though selecting
an entire column from a matrix is possible as a special form of subscripting.

@< Type definitions @>=
struct subscr_base : public expression_base
{ enum sub_type
  { row_entry, vector_entry, ratvec_entry, string_char
  , matrix_entry, matrix_column, mod_poly_term, not_so };
  expression_ptr array, index; // the two parts of the subscription expression
@)
  subscr_base(expression_ptr&& a, expression_ptr&& i)
@/: array(a.release())
  , index(i.release())
  @+{}
  virtual ~@[subscr_base() nothing_new_here@] ;
@)
  void print (std::ostream& out, bool reversed) const; // non |virtual|
  static sub_type index_kind
  (const type_expr& aggr,
   const type_expr& index,
         type_expr& subscr);
  static sub_type slice_kind(const type_expr& aggr);
  static bool assignable(sub_type);
};
@)
struct slice_base : public expression_base
{ expression_ptr array, lower,upper; // the three parts of the slice expression
@)
  slice_base(expression_ptr&& a, expression_ptr&& l, expression_ptr&& u)
@/: array(a.release())
  , lower(l.release())
  , upper(u.release())
  @+{}
  virtual ~@[slice_base() nothing_new_here@] ;
@)
  void print(std::ostream& out, unsigned flags) const; // non |virtual|
};

@ We derive a number of types from |subscr_base| which only differ by their
|evaluate| methods.

@< Type definitions @>=

template <bool reversed>
struct row_subscription : public subscr_base
{ row_subscription(expression_ptr&& a, expression_ptr&& i)
@/: subscr_base(std::move(a),std::move(i)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; subscr_base::print(out,reversed); }
};

@)
template <bool reversed>
struct vector_subscription : public subscr_base
{ vector_subscription(expression_ptr&& a, expression_ptr&& i)
@/: subscr_base(std::move(a),std::move(i)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; subscr_base::print(out,reversed); }
};
@)
template <bool reversed>
struct ratvec_subscription : public subscr_base
{ ratvec_subscription(expression_ptr&& a, expression_ptr&& i)
@/: subscr_base(std::move(a),std::move(i)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; subscr_base::print(out,reversed); }
};
@)
template <bool reversed>
struct string_subscription : public subscr_base
{ string_subscription(expression_ptr&& a, expression_ptr&& i)
@/: subscr_base(std::move(a),std::move(i)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; subscr_base::print(out,reversed); }
};
@)
template <bool reversed>
struct matrix_subscription : public subscr_base
{ matrix_subscription(expression_ptr&& a, expression_ptr&& ij)
@/: subscr_base(std::move(a),std::move(ij)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const; // slightly more complicated
};
@)
template <bool reversed>
struct matrix_get_column : public subscr_base
{ matrix_get_column(expression_ptr&& a, expression_ptr&& j)
@/: subscr_base(std::move(a),std::move(j)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; subscr_base::print(out,reversed); }
};
@)
struct module_coefficient : public subscr_base
{ module_coefficient(expression_ptr&& pol, expression_ptr&& param)
@/: subscr_base(std::move(pol),std::move(param)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; out << *array << '[' << ']'; }
};


@ We derive a number of types, templated over |unsigned|, from |slice_base|.
The template value represents the optional reversals ($8$ possibilities). All
of these differ only by their |evaluate| method.

@< Type definitions @>=

template <unsigned flags>
struct row_slice : public slice_base
{ row_slice(expression_ptr&& a, expression_ptr&& l,  expression_ptr&& u)
@/: slice_base(std::move(a),std::move(l),std::move(u)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; slice_base::print(out,flags); }
};

@)
template <unsigned flags>
struct vector_slice : public slice_base
{ vector_slice(expression_ptr&& a, expression_ptr&& l,  expression_ptr&& u)
@/: slice_base(std::move(a),std::move(l),std::move(u)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; slice_base::print(out,flags); }
};
@)
template <unsigned flags>
struct ratvec_slice : public slice_base
{ ratvec_slice(expression_ptr&& a, expression_ptr&& l,  expression_ptr&& u)
@/: slice_base(std::move(a),std::move(l),std::move(u)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; slice_base::print(out,flags); }
};
@)
template <unsigned flags>
struct string_slice : public slice_base
{ string_slice(expression_ptr&& a, expression_ptr&& l,  expression_ptr&& u)
@/: slice_base(std::move(a),std::move(l),std::move(u)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; slice_base::print(out,flags); }
};
@)
template <unsigned flags>
struct matrix_slice : public slice_base
{ matrix_slice(expression_ptr&& a, expression_ptr&& l,  expression_ptr&& u)
@/: slice_base(std::move(a),std::move(l),std::move(u)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; slice_base::print(out,flags); }
};


@ These subscriptions and slices are printed in the usual syntax. The
templated virtual |print| methods transform their template argument to a
runtime value, and call the non |virtual| member |slice_base::print| with that
value as second argument. For matrix subscriptions, where the index type
is \.{(int,int)}, the index expression is quite likely to be a tuple display,
in which case we suppress parentheses. Since we have passed the type check
here, we know that any tuple display is necessarily a pair.

@< Function definitions @>=
void subscr_base::print(std::ostream& out,bool reversed) const
{@; out << *array << (reversed ? "~[" : "[") << *index << ']';
}
@)
template <bool reversed>
void matrix_subscription<reversed>::print(std::ostream& out) const
{ out  << *array << (reversed ? "~[" : "[");
  tuple_expression* p=dynamic_cast<tuple_expression*>(index.get());
  if (p==nullptr) out << *index << ']';
  else
    out << *p->component[0] << ',' << *p->component[1] << ']';
}
@)
void slice_base::print(std::ostream& out, unsigned flags) const
{  out << *array << ((flags&0x1)==0 ? "[" : "~[")
@|     << *lower << ((flags&0x2)==0 ? ":" : "~:")
@|     << *upper << ((flags&0x4)==0 ? "]" : "~]");
}

@ It shall be useful to have a function recognising valid aggregate-index
combinations. Upon success, the last two parameters serve to store the type
the subscription will result in, and an element of the |sub_type| enumeration
that indicates the kind of subscription that was found. The |mod_poly_term|
case indicates that a ``parameter polynomial'' can be subscripted with a
parameter to return a split integer result.

@< Function def... @>=
subscr_base::sub_type subscr_base::index_kind
  (const type_expr& aggr,
   const type_expr& index,
         type_expr& subscr)
{ if (aggr.kind==primitive_type)
    switch (aggr.prim)
    {  default:
    break; case vector_type:
      if (index==int_type and subscr.specialise(int_type))
        return vector_entry; @+
    break; case rational_vector_type:
      if (index==int_type and subscr.specialise(rat_type))
        return ratvec_entry; @+
    break; case string_type:
      if (index==int_type and subscr.specialise(str_type))
        return string_char; @+
    break; case matrix_type:
      if (index==int_int_type and subscr.specialise(int_type))
        return matrix_entry;
      else if (index==int_type and subscr.specialise(vec_type))
        return matrix_column;
    break; case virtual_module_type:
      if (index==param_type and subscr.specialise(split_type))
        return mod_poly_term;
    }
  else if (aggr.kind==row_type and index==int_type and
           subscr.specialise(*aggr.component_type))
         return row_entry;
  return not_so;
}
@)
subscr_base::sub_type subscr_base::slice_kind (const type_expr& aggr)
{ if (aggr.kind==primitive_type)
    switch (aggr.prim)
    {
    case vector_type: return vector_entry;
    case rational_vector_type: return ratvec_entry;
    case string_type: return string_char;
    case matrix_type: return matrix_column;
    default: return not_so;
    }
  else if (aggr.kind==row_type)
    return row_entry;
  else return not_so;
}

@ Some cases, although valid as subscriptions, do not allow a new value to be
assigned to the component value (this holds for instance selecting a character
from a string).

@< Function def... @>=
bool subscr_base::assignable(subscr_base::sub_type t)
{ switch (t)
  { case ratvec_entry: case string_char: case mod_poly_term:
    case not_so: return false;
    default: return true;
  }
}


@ When encountering a subscription in |convert_expr|, we determine the types
of array and of the indexing expression separately, ignoring so far any type
required by the context. Then we look if the types agree with any of the types
of subscription expressions that we can convert to, throwing an error if it
does not. Finally we check if the a priori type |subscr_type| of the
subscripted expression equals or specialises to the required |type|, or can be
converted to it by |coerce|, again throwing an error if nothing works. For the
indexing expression only equality of types is admitted, since the basic
language has no conversions that could apply, and if an extension does provide
some (like the implicit conversion from split integer to pair of integers that
at some time existed), we would not want to apply them implicitly in index
positions. Also there is little point in catering for (indexing) expressions
having completely undetermined type, as such a type can only apply to an
expression that can never return (it cannot be evaluated without error).

@< Cases for type-checking and converting... @>=
case subscription:
{ type_expr array_type, index_type, subscr_type;
    // all initialised to |undetermined_type|
  auto& subsn = *e.subscription_variant;
    // alias |struct| pointed to by the |unique_ptr|
  expression_ptr array = convert_expr(subsn.array,array_type);
  expression_ptr index = convert_expr(subsn.index,index_type);
  expression_ptr subscr;
  @< Set |subscr| to a pointer to a subscription of a kind determined by
     |array_type|, |index_type| while setting |subscr_type|, and holding
     pointers moved from |array| and |index|, or |throw| an error @>
  return conform_types(subscr_type,type,std::move(subscr),e);
}

@ This is a large |switch| statement (the first of several) that is required to
separately and explicitly specify each class template instance that our
program uses (there are $13$ of them here).

The decision whether the subscription is allowed, and what will be the
resulting |subscr_type| are made by the static method |subscr_base::index_kind|.

@< Set |subscr| to a pointer to a subscription of a kind determined by...@>=
switch (subscr_base::index_kind(array_type,index_type,subscr_type))
{ case subscr_base::row_entry:
  if (subsn.reversed)
    subscr.reset(new
      row_subscription<true>(std::move(array),std::move(index)));
  else
    subscr.reset(new
      row_subscription<false>(std::move(array),std::move(index)));
break;
case subscr_base::vector_entry:
  if (subsn.reversed)
    subscr.reset(new
      vector_subscription<true>(std::move(array),std::move(index)));
  else
    subscr.reset(new
      vector_subscription<false>(std::move(array),std::move(index)));
break;
case subscr_base::ratvec_entry:
  if (subsn.reversed)
    subscr.reset(new
      ratvec_subscription<true>(std::move(array),std::move(index)));
  else
    subscr.reset(new
      ratvec_subscription<false>(std::move(array),std::move(index)));
break;
case subscr_base::string_char:
  if (subsn.reversed)
    subscr.reset(new
      string_subscription<true>(std::move(array),std::move(index)));
  else
    subscr.reset(new
      string_subscription<false>(std::move(array),std::move(index)));
break;
case subscr_base::matrix_entry:
  if (subsn.reversed)
    subscr.reset(new
      matrix_subscription<true>(std::move(array),std::move(index)));
  else
    subscr.reset(new
      matrix_subscription<false>(std::move(array),std::move(index)));
break;
case subscr_base::matrix_column:
  if (subsn.reversed)
    subscr.reset(new
      matrix_get_column<true>(std::move(array),std::move(index)));
  else
    subscr.reset(new
      matrix_get_column<false>(std::move(array),std::move(index)));
break;
case subscr_base::mod_poly_term:
  if (subsn.reversed)
    throw expr_error(e,"Cannot do reversed subscription of a ParamPol");
  subscr.reset(new module_coefficient(std::move(array),std::move(index)));
break;
case subscr_base::not_so:
  std::ostringstream o;
  o << "Cannot subscript value of type " << array_type @|
    << " with index of type " << index_type;
  throw expr_error(e,o.str());
}

@ For slices we shall similarly need $5$ kinds of slice each with $8$ values
of the template parameter |flags| for $40$ classes in all. Convert a runtime
values |flags| a template argument (which must be a compile time constant) can
basically only be done by listing all applicable values. To avoid extreme
repetitiveness, we use for this a function that is itself templated over the
class template that takes |flags| as template argument; there will be $5$
different such class templates used in calls of |make_slice|.

@< Local function definitions @>=
template < @[ template < unsigned > class @+ slice @] >
expression make_slice(unsigned flags
  ,expression_ptr&& a, expression_ptr&& b, expression_ptr&& c)
{ switch (flags)
  {
  case 0: return new slice<0>(std::move(a),std::move(b),std::move(c));
  case 1: return new slice<1>(std::move(a),std::move(b),std::move(c));
  case 2: return new slice<2>(std::move(a),std::move(b),std::move(c));
  case 3: return new slice<3>(std::move(a),std::move(b),std::move(c));
  case 4: return new slice<4>(std::move(a),std::move(b),std::move(c));
  case 5: return new slice<5>(std::move(a),std::move(b),std::move(c));
  case 6: return new slice<6>(std::move(a),std::move(b),std::move(c));
  case 7: return new slice<7>(std::move(a),std::move(b),std::move(c));
  default: assert(false); return(nullptr);
  }
}

@ When encountering a slice in |convert_expr|, we convert the array expression
with the same type required by the context (since the slice will not change
the type), and the bound expressions with integer type required. Then we look
if array type is one that can be sliced at all, throwing an error if it cannot.

@< Cases for type-checking and converting... @>=
case slice:
{ type_expr array_type; // initialised to |undetermined_type|
  expression_ptr array = convert_expr(e.slice_variant->array,array_type);
  expression_ptr lower =
    convert_expr(e.slice_variant->lower,as_lvalue(int_type.copy()));
  expression_ptr upper =
    convert_expr(e.slice_variant->upper,as_lvalue(int_type.copy()));
  expression_ptr subscr; const unsigned fl = e.slice_variant->flags.to_ulong();
  switch (subscr_base::slice_kind(array_type))
    { case subscr_base::row_entry: subscr.reset(
    @| make_slice<row_slice>
        (fl,std::move(array),std::move(lower),std::move(upper)));
    break;
    case subscr_base::vector_entry: subscr.reset(
    @| make_slice<vector_slice>
        (fl,std::move(array),std::move(lower),std::move(upper)));
    break;
    case subscr_base::ratvec_entry: subscr.reset(
    @| make_slice<ratvec_slice>
        (fl,std::move(array),std::move(lower),std::move(upper)));
    break;
    case subscr_base::string_char: subscr.reset(
    @| make_slice<string_slice>
        (fl,std::move(array),std::move(lower),std::move(upper)));
    break;
    case subscr_base::matrix_column: subscr.reset(
    @| make_slice<matrix_slice>
        (fl,std::move(array),std::move(lower),std::move(upper)));
    break;
    default: std::ostringstream o;
      o << "Cannot slice value of type " << array_type;
      throw expr_error(e,o.str());
    }
@)
  return conform_types(array_type,type,std::move(subscr),e);
}


@ Here are the |evaluate| methods for the various subscription expressions.
They all follow the same straightforward pattern, and differ only in the way
the result value push on the stack is constructed. The |static_cast<unsigned
int>| allows a range check of the (signed) integer index with a single
comparison against the unsigned array size; in the error message the signed
quantity is transmitted however.

@< Function definitions @>=
inline std::string range_mess
  (int i,size_t n,const expression_base* e,const char* where)
{ std::ostringstream o;
  e->print(o << "index " << i << " out of range (0<= . <" << n
             << ") in " << where << ' ');
  return o.str();
}
@)
template <bool reversed>
void row_subscription<reversed>::evaluate(level l) const
{ int i=(index->eval(),get<int_value>()->val);
  shared_row r=(array->eval(),get<row_value>());
  size_t n = r->val.size();
  if (reversed)
    i=n-1-i;
  if (static_cast<unsigned int>(i)>=n)
    throw runtime_error(range_mess(i,n,this,"subscription"));
  push_expanded(l,r->val[i]);
}
@)
template <bool reversed>
void vector_subscription<reversed>::evaluate(level l) const
{ int i=(index->eval(),get<int_value>()->val);
  shared_vector v=(array->eval(),get<vector_value>());
  size_t n = v->val.size();
  if (reversed)
    i=n-1-i;
  if (static_cast<unsigned int>(i)>=n)
    throw runtime_error(range_mess(i,n,this,"subscription"));
  if (l!=no_value)
    push_value(std::make_shared<int_value>(v->val[i]));
}
@)
template <bool reversed>
void ratvec_subscription<reversed>::evaluate(level l) const
{ int i=(index->eval(),get<int_value>()->val);
  shared_rational_vector v=(array->eval(),get<rational_vector_value>());
  size_t n = v->val.size();
  if (reversed)
    i=n-1-i;
  if (static_cast<unsigned int>(i)>=n)
    throw runtime_error(range_mess(i,n,this,"subscription"));
  if (l!=no_value)
    push_value(std::make_shared<rat_value>(Rational @|
       (v->val.numerator()[i],v->val.denominator())));
}
@)
template <bool reversed>
void string_subscription<reversed>::evaluate(level l) const
{ int i=(index->eval(),get<int_value>()->val);
  shared_string s=(array->eval(),get<string_value>());
  size_t n = s->val.size();
  if (reversed)
    i=n-1-i;
  if (static_cast<unsigned int>(i)>=n)
    throw runtime_error(range_mess(i,n,this,"subscription"));
  if (l!=no_value)
    push_value(std::make_shared<string_value>(s->val.substr(i,1)));
}

@ And here are the cases for matrix indexing and column selection, which are
just slightly more complicated.

@< Function definitions @>=
template <bool reversed>
void matrix_subscription<reversed>::evaluate(level l) const
{ index->multi_eval(); @+
  int j=get<int_value>()->val;
  int i=get<int_value>()->val;
  shared_matrix m=(array->eval(),get<matrix_value>());
  size_t r = m->val.numRows();
  size_t c = m->val.numColumns();
  if (reversed)
  {@;  i=r-1-i; j=c-1-j; }
  if (static_cast<unsigned int>(i)>=r)
    throw runtime_error
     ("initial "+range_mess(i,r,this,"matrix subscription"));
  if (static_cast<unsigned int>(j)>=c)
    throw runtime_error
     ("final "+range_mess(j,c,this,"matrix subscription"));
  if (l!=no_value)
    push_value(std::make_shared<int_value>(m->val(i,j)));
}
@)
template <bool reversed>
void matrix_get_column<reversed>::evaluate(level l) const
{ int j=(index->eval(),get<int_value>()->val);
  shared_matrix m=(array->eval(),get<matrix_value>());
  size_t c = m->val.numColumns();
  if (reversed)
    j=c-1-j;
  if (static_cast<unsigned int>(j)>=c)
    throw runtime_error(range_mess(j,c,this,"matrix column selection"));
  if (l!=no_value)
    push_value(std::make_shared<vector_value>(m->val.column(j)));
}

@ For slice these are template functions. This is where the actual reversals
happen.
@< Function definitions @>=
void slice_range_error
  (int lwb, int upb,int n,unsigned flags,const expression_base* e)
{ std::ostringstream o;
  if (upb>n)
    if (lwb<0)
      o << "both bounds " << lwb << ':' << upb;
    else
      o << "upper bound " << upb;
  else
    o << "lower bound " << lwb;
  o << " out of range (should be ";
  if (lwb<0)
    o << ">=0" << (upb>n ? " respectively " : ")");
  if (upb>n)
    o << "<=" << n << ')';
  e->print(o << " in slice ");
  throw runtime_error(o.str());
}
@)
template <unsigned flags>
void row_slice<flags>::evaluate(level l) const
{ int upb=(upper->eval(),get<int_value>()->val);
  int lwb=(lower->eval(),get<int_value>()->val);
  shared_row arr=(array->eval(),get<row_value>());
  const auto& r = arr->val;
  int n = r.size();
  if ((flags&0x2)!=0)
    lwb = n - lwb;
  if ((flags&0x4)!=0)
    upb = n - upb;
  if (lwb<0 or upb>n)
    slice_range_error(lwb,upb,n,flags,this);
  if (lwb>=upb)
  {@; push_value(std::make_shared<row_value>(0)); return; }
  push_value( (flags&0x1)==0
@|  ? std::make_shared<row_value>(r.begin()+lwb,r.begin()+upb)
@|  : std::make_shared<row_value>(r.rbegin()+lwb,r.rbegin()+upb)
    );
}

@ For vectors this is quite similar.
@< Function definitions @>=

template <unsigned flags>
void vector_slice<flags>::evaluate(level l) const
{ int upb=(upper->eval(),get<int_value>()->val);
  int lwb=(lower->eval(),get<int_value>()->val);
  shared_vector arr=(array->eval(),get<vector_value>());
  const auto& r = arr->val;
  int n = r.size();
  if ((flags&0x2)!=0)
    lwb = n - lwb;
  if ((flags&0x4)!=0)
    upb = n - upb;
  if (lwb<0 or upb>n)
    slice_range_error(lwb,upb,n,flags,this);
  if (lwb>=upb)
  {@; push_value(std::make_shared<vector_value>(int_Vector(0))); return; }
  push_value( (flags&0x1)==0
@|  ? std::make_shared<vector_value>(r.begin()+lwb,r.begin()+upb)
@|  : std::make_shared<vector_value>(r.rbegin()+lwb,r.rbegin()+upb)
    );
}

@ For rational vectors this is only slightly more complicated. We select from
the numerator, and provide a common denominator; the |rational_vector_value|
constructor takes care of normalising the result.

@< Function definitions @>=

template <unsigned flags>
void ratvec_slice<flags>::evaluate(level l) const
{ int upb=(upper->eval(),get<int_value>()->val);
  int lwb=(lower->eval(),get<int_value>()->val);
  shared_rational_vector arr=(array->eval(),get<rational_vector_value>());
  const auto& r = arr->val.numerator();
  int n = r.size();
  if ((flags&0x2)!=0)
    lwb = n - lwb;
  if ((flags&0x4)!=0)
    upb = n - upb;
  if (lwb<0 or upb>n)
    slice_range_error(lwb,upb,n,flags,this);
  if (lwb>=upb)
  {@; push_value(std::make_shared<rational_vector_value>(RatWeight(0)));
      return; }
  const int d=arr->val.denominator();
  push_value( (flags&0x1)==0
@|  ? std::make_shared<rational_vector_value>(r.begin()+lwb,r.begin()+upb,d)
@|  : std::make_shared<rational_vector_value>(r.rbegin()+lwb,r.rbegin()+upb,d)
    );
}


@ For strings it is quite easy too.
@< Function definitions @>=

template <unsigned flags>
void string_slice<flags>::evaluate(level l) const
{ int upb=(upper->eval(),get<int_value>()->val);
  int lwb=(lower->eval(),get<int_value>()->val);
  shared_string arr=(array->eval(),get<string_value>());
  const auto& r = arr->val;
  int n = r.size();
  if ((flags&0x2)!=0)
    lwb = n - lwb;
  if ((flags&0x4)!=0)
    upb = n - upb;
  if (lwb<0 or upb>n)
    slice_range_error(lwb,upb,n,flags,this);
  if (lwb>=upb)
  {@; push_value(std::make_shared<string_value>(std::string())); return; }
  push_value
    ((flags&0x1)==0
@|  ? std::make_shared<string_value>(std::string(&r[lwb],&r[upb]))
@|  : std::make_shared<string_value>(@|
       std::string(r.rbegin()+lwb,r.rbegin()+upb)));
}

@ And here is the version for matrices. Here we cannot use reverse iterators,
but it is not hard to do the right thing using the integer indices |lwb|,
|upb| directly.

@< Function definitions @>=

template <unsigned flags>
void matrix_slice<flags>::evaluate(level l) const
{ int upb=(upper->eval(),get<int_value>()->val);
  int lwb=(lower->eval(),get<int_value>()->val);
  shared_matrix mat=(array->eval(),get<matrix_value>());
  const auto& m = mat->val;
  int n = m.numColumns();
  if ((flags&0x2)!=0)
    lwb = n - lwb;
  if ((flags&0x4)!=0)
    upb = n - upb;
  if (lwb<0 or upb>n)
    slice_range_error(lwb,upb,n,flags,this);
  if (lwb>=upb)
@/{@; push_value(std::make_shared<matrix_value>(int_Matrix(m.numRows(),0)));
      return; }
  own_matrix result =
    std::make_shared<matrix_value>(int_Matrix(m.numRows(),upb-lwb));
  auto& r = result->val;
  for (unsigned int j=0; lwb<upb; ++j)
    r.set_column(j,m.column((flags&0x1)==0 ? lwb++ : n - ++lwb));
  push_value(std::move(result));
}

@* Control structures.
%
We shall now introduce conventional control structures, which must of course
be part of any serious programming language; yet they were implemented only
after plenty of other language elements were in place, such as
let-expressions, functions, rows and selection form them, implicit
conversions.

In fact, the power of functions is such that certain control structures can be
simulated with functions. For instance, there are currently no multi-way
choice expressions (also known as |case| or |switch| statements), but one can
define an row of anonymous functions without parameters, use row selection to
pick out one of them based on an integer value, and then call it. Recursion is
not directly supported either, because anonymous $\lambda$-expressions have no
way of referring to themselves; however recursion can be realised assigning to
a local function variable as new value a $\lambda$-expression referring to
that same variable.

Nonetheless a good repertoire of control structures is essential to easy
programming. We provide the conditional expression as unique selection
statement, and a large variety of iterative statements. Somewhat unusual is
the fact that all control structures are expressions that may yield a value;
in the case of loop statements, a value of ``row-of'' type is returned.

@*1 The Boolean negation.
%
Though not really a control structure, Boolean negation is mostly used in
conjunction with control structures. When it is, it can most often be handled
by translation into a variation of the control structure, such as swapping the
branches of a conditional expression. However some uses cannot be so absorbed
and do need some runtime action; we translate these into calls of a
built-in function |boolean_negate_builtin| that will be defined below, along
with the generic functions.

The type check is easy; the argument must be of |bool| type (or convertible
to it, but nothing else is so) and so must (become) the required |type|.

@< Cases for type-checking and converting... @>=
case negation_expr:
{ type_expr b=bool_type.copy();
   expression_ptr arg = convert_expr(*e.negation_variant,b);
  if (not type.specialise(b)) // |not| preserves the |bool| type
    throw type_error(e,std::move(b),std::move(type));
  return expression_ptr(new @| overloaded_builtin_call
     (boolean_negate_builtin,std::move(arg),e.loc));
}

@ Often however a Boolean negation can be eliminated entirely from the
translated program by changing the syntax tree. For instance if present as
outermost operation in the condition of a conditional expression, then it
suffices to interchange the two branches to achieve the effect of the
negation. To help make such simplifications, here is a function that detects
and removes a top level negation from an expression, and returns whether it
found one. (The caller cannot afford ignoring the return value, since the
expression has already been (potentially) modified when this function
returns.) Although it is unlikely to find multiple directly nested
negations, it cost us very little to handle that case as well, so we do.

The memory management is somewhat subtle here, but is handled inside the
move-assignment operator for |expr|, which we had to rewrite when this code
was added (previously it did not handle this specific case correctly). In case
a negation was found, we shall copy a subexpression of~|e| to~|e|. This
requires moving the (top-level) value to be copied to a temporary variable
(while marking the original as having become |no_expr|) before destructing the
old value of~|e|, lest the value to be copied be thrown away in the process of
recursive destruction. This is similar to how move-assignment of
|std::unique_ptr| instances ensures the proper timing (releasing the pointer
to be assigned to a temporary raw pointer before destructing the old value of
the pointer being assigned to) in case of recursive structures linked by such
smart pointers. But we could not use that here, since in |expr| the links are
raw pointers to structures containing |expr| fields: our caller does not hold
any pointer pointing directly to the argument~|e| it passes, so we cannot
operate merely by changing some pointers.

The negation hunting works through expression types that return one of their
component expressions: |let| expressions, sequence expressions, and \&{next}
expressions; this is easily achieved by recursive calls. More importantly it
also handles the |do|-expressions that combine condition and loop body inside
a |while| loop, in which case the recursion goes into the (first) condition
part, since that is the part where the negations are to be removed.

@< Local function definitions @>=

bool was_negated (expr& e, bool b)
{ switch (e.kind)
  { default: return b;
  case negation_expr:
    e=std::move(*e.negation_variant); // remove top level |not|
    return was_negated(e,not b); // recursively continue
  case let_expr: return was_negated(e.let_variant->body,b);
  case seq_expr: return was_negated(e.sequence_variant->last,b);
  case next_expr: return was_negated(e.sequence_variant->first,b);
  }
}

bool was_negated(expr& e) @+{@; return was_negated(e,false); }

@*1 Conditional expressions.
%
A first control structure it the conditional expression.

@< Type def... @>=
struct conditional_expression : public expression_base
{ expression_ptr condition, then_branch, else_branch;
@)
  conditional_expression
   (expression_ptr&& c,expression_ptr&& t, expression_ptr&& e)
   : condition(c.release()),then_branch(t.release()), else_branch(e.release())
  @+{}
  virtual ~@[conditional_expression() nothing_new_here@];
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ To print a conditional expression, we reconstruct \&{elif} constructions
that were eliminated in the parser (and even those that the user did not
employ, but could have). To this end, we remain in a loop as long as the
|else|-part is itself a conditional expression. This makes a loop with exit in
the middle a natural solution, and in any case \Cpp\ does not allow using a
variable introduced in the body, like |p| below, to be used in the condition
of a |while| or |for| controlling the loop.

@< Function definitions @>=
void conditional_expression::print(std::ostream& out) const
{ out << " if "; const conditional_expression* cur=this;
  do
  { out << *cur->condition << " then " << *cur->then_branch;
    conditional_expression* p =
      dynamic_cast<conditional_expression*>(cur->else_branch.get());
    if (p==nullptr)
      break;
    out << " elif "; cur=p;
  }
  while(true);
  out << " else " << *cur->else_branch << " fi ";
}

@ Evaluating a conditional expression ends up evaluating either the
then-branch or the else-branch.

@< Function definitions @>=
void conditional_expression::evaluate(level l) const
{ condition->eval();
  if (get<bool_value>()->val)
    then_branch->evaluate(l);
  @+ else else_branch->evaluate(l);
}

@ With the function |balance| defined above, and thanks to the fact that we
somewhat artificially arranged the then and else branches to be accessible as a
|raw_expr_list|, conversion of conditionals has become easy.

@< Cases for type-checking and converting... @>=
case conditional_expr:
{ auto& exp = *e.if_variant;
  if (was_negated(exp.condition)) // eliminate negated conditions
    exp.branches.contents.swap(exp.branches.next->contents);
  expression_ptr c  =  convert_expr(exp.condition,as_lvalue(bool_type.copy()));
  std::vector<expression_ptr> conv;
  balance(type,&exp.branches,e,"branches of conditional",conv);
@/return expression_ptr(new @|
    conditional_expression(std::move(c),std::move(conv[0]),std::move(conv[1])));
}

@*1 Integer controlled case expressions.
%
The integer case expression (multi-way branch controlled by an integer value)
is quite similar to a conditional, but contains a list of branches rather than
two of them.

@< Type def... @>=
struct int_case_expression : public expression_base
{ expression_ptr condition; std::vector<expression_ptr> branches;
@)
  int_case_expression
   (expression_ptr&& c,std::vector<expression_ptr>&& b)
   : condition(c.release()),branches(std::move(b))
  @+{}
  virtual ~@[int_case_expression() nothing_new_here@];
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ To print a case expression is straightforward.

@< Function definitions @>=
void int_case_expression::print(std::ostream& out) const
{ auto it = branches.cbegin();
  assert(it!=branches.cend());
  out << " case " << *condition << " in " << **it;
  while (++it!=branches.cend())
    out << ", " << **it;
  out << " esac ";
}

@ Evaluating a case expression ends up evaluating one of the |branches|.

@< Function definitions @>=
void int_case_expression::evaluate(level l) const
{ condition->eval();
  int i = get<int_value>()->val;
  if (static_cast<unsigned>(i)>=branches.size())
    throw runtime_error(range_mess(i,branches.size(),this,"case expression"));
  branches[i]->evaluate(l);
}

@ With the function |balance| defined above, conversion of case expressions
has become easy.

@< Cases for type-checking and converting... @>=
case int_case_expr:
{ auto& exp = *e.if_variant;
  expression_ptr c  =  convert_expr(exp.condition,as_lvalue(int_type.copy()));
  std::vector<expression_ptr> conv;
  balance(type,&exp.branches,e,"branches of case",conv);
@/return expression_ptr(new @|
    int_case_expression(std::move(c),std::move(conv)));
}

@*1 While loops.
%
Next we consider different kinds of loops. Apart from the distinction between
loops controlled by a condition (|while|), by the components of some value
(|for|) or by a simple counter (counted |for|), each kind may have several
variants, determined at compile time and transformed into a template argument
to the classes implementing the loops. The variant is determined by a
combination of bits, like one indicating reverse traversal of the ``input''
range (for |for| loops only) and one indicating reverse accumulation of the
loop body values into a row (for any kind of loop). For readability these bits
will be tested by the following predicates; since |flags| will be a template
parameter at all uses, the predicates will have a constant value in all cases,
whence the |constexpr| declaration.

@s constexpr const

@< Local function def... @>=
constexpr bool in_reversed(unsigned flags) @+{@; return (flags&0x1)!=0; }
constexpr bool in_forward(unsigned flags) @+{@; return (flags&0x1)==0; }
constexpr bool out_reversed(unsigned flags) @+{@; return (flags&0x2)!=0; }
constexpr bool out_forward(unsigned flags) @+{@; return (flags&0x2)==0; }
constexpr bool no_frame(unsigned flags) @+{@; return (flags&0x4)!=0; }
constexpr bool has_frame(unsigned flags) @+{@; return (flags&0x4)==0; }
constexpr bool yields_count(unsigned flags) @+{@; return (flags&0x8)!=0; }

@ Let us start with considering |while| loops.
%
Although they contain two parts, a condition before |do| and a body after it,
the two are present in the following structure as a single |body|. The reason
for this is that we want to allow declarations in the condition to remain
valid in the body, which can be realised by having both contained in a
|do_expression| structure that does not have to be a direct descendent of the
|while_expression|, but can involve a number of intermediate |let_expression|
and |seq_expression| nodes. The grammar guarantees that eventually such a
|do_expression| will be reached.

@< Type def... @>=
template<unsigned flags>
struct while_expression : public expression_base
{ expression_ptr body;
  while_expression(expression_ptr&& b): body(std::move(b)) @+{}
  virtual ~@[while_expression() nothing_new_here@];
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ In printing of while loops, the symbol \.{do} or \.{\char\~do} with appear
in the output for |*body|.

@< Function definitions @>=
template<unsigned flags>
void while_expression<flags>::print(std::ostream& out) const
{@; out << " while " << *body << (out_reversed(flags) ? " ~od " : " od "); }

@ The following function helps instantiating the class template, selecting a
constant template parameter equal to |flags|. The only flag bits used in while
loops are |0x2| for reverse accumulation of values from the loop body, and
|0x8| indicating that only the number of iterations done should be returned as
value for the loop. The two are mutually exclusive, so we have three cases
here.

@< Local function def... @>=
expression make_while_loop (unsigned flags,expression_ptr&& body)
{ switch(flags)
  {
  case 0x0: return new while_expression<0>(std::move(body));
  case 0x2: return new while_expression<2>(std::move(body));
  case 0x8: return new while_expression<8>(std::move(body));
  default: assert(false); return nullptr;
  }
}

@ Type-checking for |while| loops has a few complications because possibly a
row result must be produced from the loop body expression. If the context
requires void type, we require the same for the body type, knowing that
generation of a row value will be suppressed in these cases anyway. (This is
more flexible than leaving it undetermined as we did in the somewhat similar
case of a list display; in contrast to that case, a loop in void context is
quite normal, and making their body void allows for instance conditionals
contained in them to have branches of not necessarily compatible types.) In
all other cases we proceed for the body expression as for the components of a
list display (except that there is only one expression in this case).

@< Cases for type-checking and converting... @>=
case while_expr:
{ while_node& w=*e.while_variant;
  layer bind(0,nullptr); // no local variables for loop, but allow |break|
@)
  if (type==void_type)
  { expression_ptr result(make_while_loop @| (w.flags.to_ulong(),
       convert_expr(w.body, as_lvalue(void_type.copy()))));
    return expression_ptr(new voiding(std::move(result)));
  }
  else if (type==int_type)
    return expression_ptr(make_while_loop @| (0x8,
       convert_expr(w.body, as_lvalue(void_type.copy()))));
  else if (type.specialise(row_of_type))
  { auto& comp_type = *type.component_type;
    expression_ptr b = convert_expr(w.body,comp_type);
    if (comp_type==void_type and not is_empty(w.body))
      b.reset(new voiding(std::move(b)));
    return expression_ptr (make_while_loop @|
       (w.flags.to_ulong(),std::move(b)));
  }
  else
  @< If |type| can be converted from some row-of type, check |w.body|
     against its component type, construct the |while_expression|, and apply
     the appropriate conversion function to it; otherwise |throw| a
     |type_error| @>
}

@ For |while| loops we follow the same logic for finding an appropriate
component type as for list displays, in section@#list display conversion@>.

@< If |type| can be converted from some row-of type, check |w.body| against
   its component type, construct the |while_expression|, and apply the
   appropriate conversion function to it; otherwise |throw| a |type_error| @>=
{ type_expr comp_type;
  const conversion_record* conv = row_coercion(type,comp_type);
  if (conv==nullptr)
    throw type_error(e,row_of_type.copy(),std::move(type));
@)
  return expression_ptr(new conversion(*conv, expression_ptr(make_while_loop @|
       (w.flags.to_ulong(),convert_expr(w.body,comp_type)))));
}

@ Of course evaluating is what most distinguishes loops from conditionals. The
work is partly delegated to the evaluation of the special expression in the
body that determines both the loop guard condition and if positive the loop
body. Afterwards it needs to transmit what it has done back to the enclosing
loop. It could do this on the |execution_stack|, but since this condition is
quite short-lived (it just has to maybe traverse the destruction of a number
of stack |frame|s for |let| declarations that terminate), we can avoid
wrapping into a |shared_value| and unwrapping by reserving a static |bool|
variable for this purpose.

@< Local variable definitions @>=
bool while_condition_result;

@ The mentioned special expression is a |do_expresision|, which is a slight
variation of |seq_expression|.

@< Type def... @>=
template <bool negated>
struct do_expression : public expression_base
{ expression_ptr condition,body;
@)
  do_expression(expression_ptr&& c,expression_ptr&& b)
   : condition(c.release()),body(b.release()) @+{}
  virtual ~@[do_expression() nothing_new_here@];
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ For efficiency we also define two variants without condition, for in case
the condition is given simply by |true| or |false|. In the former case we just
have a |body| that will always be evaluated and (if no interruption occurs)
set |while_condition_result=true|, in the latter case not even a body is
necessary, and evaluation just sets |while_condition_result=false|.

@< Type def... @>=
struct forever_expression : public expression_base
{ expression_ptr body;
@)
  forever_expression(expression_ptr&& b) : body(b.release()) @+{}
  virtual ~@[forever_expression() nothing_new_here@];
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};
struct dont_expression : public expression_base
{ dont_expression() @+{}
  virtual ~@[dont_expression() nothing_new_here@];
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};


@ Printing a |do_expression| differs from that of |next_expression| only in
the keyword, where we print a tilde to indicate an inverted condition (giving
an until-loop rather than a while-loop so to speak), even though this is
(currently) not an allowed input syntax (writing an explicit |not| that gives
the inverted condition is deemed more readable).

@< Function definitions @>=
template <bool negated>
  void do_expression<negated>::print(std::ostream& out) const
{@; out << *condition << (negated ? " ~do " : " do ") << *body; }
void forever_expression::print(std::ostream& out) const @+
{@; out << " do " << *body; }
void dont_expression::print(std::ostream& out) const @+
{@; out << " dont "; }

@ The analysis of a |do| expression is similar to the of a sequence
expression, notably it is the type of the final (body) subexpression that is
deemed to be the type of the whole expression; however its initial
subexpression has a |bool| rather than |void| context. If the condition was
negated, we generate the negated template instance of |do_expression|.

@< Cases for type-checking and converting... @>=
case do_expr:
{ sequence_node& seq=*e.sequence_variant;
  bool neg = was_negated(seq.first);
  expression_ptr body = convert_expr(seq.last,type);
    // body needs type-checking in all cases
  if (seq.first.kind==boolean_denotation)
  { if (seq.first.bool_denotation_variant==neg)
      return expression_ptr(new dont_expression()); // and drop |body|
    else
      return expression_ptr(new forever_expression(std::move(body)));
  }
  expression_ptr condition = convert_expr(seq.first,as_lvalue(bool_type.copy()));
  if (neg)
    return expression_ptr(new @|
      do_expression<true>(std::move(condition),std::move(body)));
  else
    return expression_ptr(new @|
       do_expression<false>(std::move(condition),std::move(body)));
}
break;

@ The evaluation of a |do_expression| is somewhat special in that whether it
places a value on the runtime stack (assuming that |l!=no_value|) depends on
the outcome of evaluation the condition. Nonetheless the definition is very
simple. The only subtle point is to postpone setting |while_condition_result|
to the end of the evaluation, since the call |body->evaluate(l)| may clobber
this variable. The two explicit instances below could easily be combined into
one template definition, but since the whole point of inverted conditions is
to make this (marginally) more efficient (than evaluating a call to
|boolean_negate_builtin|), we prefer for once to be explicit rather than to
depend on the \Cpp\ compiler's constant folding doing the simplification.

@< Function def... @>=
template <>
  void do_expression<false>::evaluate(level l) const
{ condition->eval();
  if (get<bool_value>()->val)
    {@; body->evaluate(l); while_condition_result=true; }
  else while_condition_result=false;
}
template <>
  void do_expression<true>::evaluate(level l) const
{ condition->eval();
  if (get<bool_value>()->val)
     while_condition_result=false;
  else {@; body->evaluate(l); while_condition_result=true; }
}
@)
void forever_expression::evaluate(level l) const
{@; body->evaluate(l); while_condition_result=true; }
void dont_expression::evaluate(level l) const
{@; while_condition_result=false; }


@ Now we can consider the evaluation of |while| loops themselves. If no value
is asked for, we simply perform a |while|-loop at the \Cpp~level (applying
|void_eval| to the body expression). When the |yields_count| predicate holds,
we do similarly but keep a count of the loop iterations and return that value.
In the remaining cases a row value is produced; it will be detailed later.

Curiously, due to the combined evaluation of condition and loop body, of the
|void| case corresponds most naturally to a |do|-|while| loop in \Cpp, since
nothing remains to be done after testing the termination condition. This is
not true for the value-producing while loop, since the vale produced by the
loop body needs to be popped from the |execution_stack|; therefore we have to
use the comma-operator in the condition of the |while| loop here.

Since |while|-loops can run indefinitely and don't necessarily involve calls
to user-defined functions (where a test is done too), we put a test of the
user interrupt flag in the body of the evaluation of any |while|-loop. This
does of course carry a slight performance penalty, but that will hardly be
noticeable unless evaluating |body| does very little, and it is precisely in
such situations (for instance where the user forgot some updating that would
alter the value of the loop condition) where the possibility of a user
interrupt may be vital.

@< Function definitions @>=
template<unsigned flags>
void while_expression<flags>::evaluate(level l) const
{ if (l==no_value)
  { try
    { do
      { body->void_eval();
        if (interrupt_flag!=0)
          throw user_interrupt();
      }
      while (while_condition_result);
    }
    catch (loop_break& err) @+
    {@; if (err.depth-- > 0)
          throw;
    }
  }
  else if (yields_count(flags))
  { int count=0;
    try
    { while(body->void_eval(), while_condition_result)
      { ++count;
        if (interrupt_flag!=0)
          throw user_interrupt();
      }
    }
    catch (loop_break& err) @+
    {@; if (err.depth-- > 0)
          throw;
    }
    push_value(std::make_shared<int_value>(count));
  }
  else
  @< Perform a |while| loop, accumulating values from the loop bodies into a
     row that is pushed onto |execution_stack| @>
}

@ When accumulating values into a row we also use a \Cpp\ |while|-loop, but
use |eval| to produce a value on |execution_stack| each time around, popping
it off and pushing it to the front of a |simple_list| (which therefore
accumulates values in reverse order), which we move-convert into a vector, in
the appropriate order, when the loop terminates.

@< Perform a |while| loop, accumulating values... @>=
{ containers::simple_list<shared_value> result;
  size_t s=0;
  try
  { while (body->eval(),while_condition_result)
    { if (interrupt_flag!=0)
         throw user_interrupt();
      result.push_front(pop_value());
      ++s;
    }
  }
  catch (loop_break& err) @+
  {@; if (err.depth-- > 0)
        throw;
  }
  own_row r = std::make_shared<row_value>(s);
  if (out_forward(flags)) // forward accumulating while loop
  { auto dst = r->val.rbegin();
      // use reverse iterator, since we reverse accumulated
    for (auto it = result.wbegin(); not it.at_end(); ++it,++dst)
      *dst = std::move(*it);
  }
  else // reverse accumulating while loop
  { auto dst = r->val.begin(); // use ordinary iterator
    for (auto it = result.wbegin(); not it.at_end(); ++it,++dst)
      *dst = std::move(*it);
  }
  push_value(std::move(r));
}

@*1 For loops.
%
We now consider |for| loops over components of a value. They have three
parts, an identifier pattern defining the loop variable(s), an in-part giving
the object whose components are iterated over, and a body that may produce a
new value for each component of the in-part. Due to the way |for_expressions|
are constructed, |pattern| is always an unnamed $2$-tuple, the first component
of which is the name (if given) of the index of the component selected, and
the second component is the pattern for the component itself. This
organisation of the pattern is already present in the structure \&{for\_loop}
of the relevant variant of |expr|, and ultimately it is determined by parser
actions building the pattern list. (Looking at those rules, the left-to-right
reversal of component and index may not be obvious, but the function
|make_pattern_list| used takes the \emph{tail} of the list as first argument,
and the head as second; this is adapted to the subsequent reversal that
usually takes place when a pattern list is completed, but this reversal does
not happen for the $2$-element list used for the patterns in for-loops.)

Apart from iterating over row value of any type, we allow iteration over
vectors, rational vectors, strings, and matrix columns (types that are
indexable by integers), and also over the terms of a parameter polynomial
(representing isotypical components of a virtual module). The syntax of the
for loop is the same for all these cases. However the constructed |expression|
will be of a class templated on |kind| describing the kind of value iterated
over, so that their |evaluate| methods will be specialised to that kind. We
also template over the |flags| that indicate the reversal options (determined
by by the precise syntax used), so that these get built into the evaluate
method as well. Bit |0x1| of |flags| controls reverse traversal at input,
while bit |0x2| indicates reverse accumulation of values for the result.
Altogether we use $4\times6=24$ instances of |for_expression|.

@< Type def... @>=
template <unsigned flags, subscr_base::sub_type kind>
struct for_expression : public expression_base
{ id_pat pattern; expression_ptr in_part;
  expression_ptr body;
@)
  for_expression (const id_pat& p, expression_ptr&& i, expression_ptr&& b);
  virtual ~@[for_expression() nothing_new_here@];
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ We could not inline the following constructor definition in the class
declaration, as it uses the local function |copy_id_pat| that is not known in
the header file, but otherwise it is quite straightforward.

@< Function definitions @>=
template <unsigned flags, subscr_base::sub_type kind>
for_expression<flags,kind>::for_expression@|
 (const id_pat& p, expression_ptr&& i, expression_ptr&& b)
   : pattern(copy_id_pat(p)), in_part(i.release()), body(b.release())
  @+{}


@ Printing a |for| expression is straightforward, taking note that the index
part is first in the pattern though written \emph{after} \.@@, and that it
could be absent. We split off a (non templated) part for reuse in other loops,
which shows how the bits of |flags| correspond to optional tildes before the
keywords |do| and~\&{od}.

@< Function definitions @>=
void print_body(std::ostream& out,const expression_ptr& body,unsigned flags)
{@; out << (in_reversed(flags) ? " ~do " : " do ")  << *body
        << (out_reversed(flags) ? " ~od " : " od ");
}
@)
template <unsigned flags, subscr_base::sub_type kind>
void for_expression<flags,kind>::print(std::ostream& out) const
{ out << " for " << *std::next(pattern.sublist.begin());
    if (pattern.sublist.front().kind==0x1)
      out << '@@' << pattern.sublist.front();
  print_body(out << " in " << *in_part,body,flags);
}

@ As in |make_slice| above, we need to convert runtime values for |flags| and
the |sub_type t@;| to template arguments. This is quite boring, but is the
price we have to pay to get selection at type-check time of an |evaluate|
method for which the dynamic choice has been constant-folded away by the \Cpp\
compiler.

@< Local function definitions @>=
expression make_for_loop
  (unsigned flags, const id_pat& id, expression_ptr&& i, expression_ptr&& b
  ,subscr_base::sub_type t)
{ switch(t)
  {
  case subscr_base::row_entry:
    switch (flags)
    {
    case 0: return new for_expression<@[0,subscr_base::row_entry@]>
      (id,std::move(i),std::move(b));
    case 1: return new for_expression<@[1,subscr_base::row_entry@]>
      (id,std::move(i),std::move(b));
    case 2: return new for_expression<@[2,subscr_base::row_entry@]>
      (id,std::move(i),std::move(b));
    case 3: return new for_expression<@[3,subscr_base::row_entry@]>
      (id,std::move(i),std::move(b));
    default: assert(false); return(nullptr);
    }
  case subscr_base::vector_entry:
    switch (flags)
    {
    case 0: return new for_expression<@[0,subscr_base::vector_entry@]>
      (id,std::move(i),std::move(b));
    case 1: return new for_expression<@[1,subscr_base::vector_entry@]>
      (id,std::move(i),std::move(b));
    case 2: return new for_expression<@[2,subscr_base::vector_entry@]>
      (id,std::move(i),std::move(b));
    case 3: return new for_expression<@[3,subscr_base::vector_entry@]>
      (id,std::move(i),std::move(b));
    default: assert(false); return(nullptr);
    }
  case subscr_base::ratvec_entry:
    switch (flags)
    {
    case 0: return new for_expression<@[0,subscr_base::ratvec_entry@]>
      (id,std::move(i),std::move(b));
    case 1: return new for_expression<@[1,subscr_base::ratvec_entry@]>
      (id,std::move(i),std::move(b));
    case 2: return new for_expression<@[2,subscr_base::ratvec_entry@]>
      (id,std::move(i),std::move(b));
    case 3: return new for_expression<@[3,subscr_base::ratvec_entry@]>
      (id,std::move(i),std::move(b));
    default: assert(false); return(nullptr);
    }
  case subscr_base::string_char:
    switch (flags)
    {
    case 0: return new for_expression<@[0,subscr_base::string_char@]>
      (id,std::move(i),std::move(b));
    case 1: return new for_expression<@[1,subscr_base::string_char@]>
      (id,std::move(i),std::move(b));
    case 2: return new for_expression<@[2,subscr_base::string_char@]>
      (id,std::move(i),std::move(b));
    case 3: return new for_expression<@[3,subscr_base::string_char@]>
      (id,std::move(i),std::move(b));
    default: assert(false); return(nullptr);
    }
  case subscr_base::matrix_column:
    switch (flags)
    {
    case 0: return new for_expression<@[0,subscr_base::matrix_column@]>
      (id,std::move(i),std::move(b));
    case 1: return new for_expression<@[1,subscr_base::matrix_column@]>
      (id,std::move(i),std::move(b));
    case 2: return new for_expression<@[2,subscr_base::matrix_column@]>
      (id,std::move(i),std::move(b));
    case 3: return new for_expression<@[3,subscr_base::matrix_column@]>
      (id,std::move(i),std::move(b));
    default: assert(false); return(nullptr);
    }
  case subscr_base::mod_poly_term:
    switch (flags)
    {
    case 0: return new for_expression<@[0,subscr_base::mod_poly_term@]>
      (id,std::move(i),std::move(b));
    case 1: return new for_expression<@[1,subscr_base::mod_poly_term@]>
      (id,std::move(i),std::move(b));
    case 2: return new for_expression<@[2,subscr_base::mod_poly_term@]>
      (id,std::move(i),std::move(b));
    case 3: return new for_expression<@[3,subscr_base::mod_poly_term@]>
      (id,std::move(i),std::move(b));
    default: assert(false); return(nullptr);
    }
  default: assert(false); return(nullptr);
  }
}

@ Type-checking is more complicated for |for| loops than for |while| loops,
since more types and potential coercions are involved. We start by processing
the in-part in a neutral type context, which will on success set |in_type| to
its a priori type. Then after binding the loop variable(s) in a new |layer|,
we process the loop body either in neutral type context from the fresh and
subsequently ignored type |body_type| (if the loop occurs in void context), or
passing |*type.component_type| if |type| is a row type, or else if
|row_coercion| finds an applicable coercion, passing again |body_type| but now
set to the required type (if none of these apply a |type_error| is thrown).
After converting the loop, we must not forget to maybe apply voiding or a
coercion.

There is a slight twist for the rare occasion that a loop over components
introduced no identifiers at all, since |for_expression::evaluate| should not
push a forbidden empty |frame| on the evaluation context. The case is marginal
since it just allows repeated evaluation of the loop body as many times as the
number of components looped over, but it is no excluded syntactically (and
would be hard to). To avoid having to worry about this at runtime, we shall
substitute a counted loop for such loops.

@< Cases for type-checking and converting... @>=
case for_expr:
{ const for_node& f=*e.for_variant;
  type_expr in_type;
  expression_ptr in_expr = convert_expr(f.in_part,in_type);  // \&{in} part
  subscr_base::sub_type which; // the kind of aggregate iterated over
  layer bind(count_identifiers(f.id),nullptr);
   // for identifier(s) introduced in this loop
  @< Set |which| according to |in_type|, and set |bind| according to the
     identifiers contained in |f.id| @>
  type_expr body_type;
  type_p btp; // will point to the place to record body type
  const conversion_record* conv=nullptr;
    // initialise in case we don't reach the assignment below
  if (type==void_type)
    btp=&type; // we can reuse this type; no risk of specialisation
  else if (type.specialise(row_of_type))
    btp=type.component_type;
  else if ((conv=row_coercion(type,body_type))!=nullptr)
    btp=&body_type;
  else throw type_error(e,row_of_type.copy(),std::move(type));
@)
  expression_ptr body(convert_expr (f.body,*btp));
  if (type!=void_type and *btp==void_type and not is_empty(f.body))
    body.reset(new voiding(std::move(body)));
@/expression_ptr loop;
  if (bind.empty()) // we must avoid having an empty |frame| produced at runtime
    @< Set |loop| to a index-less counted |for| loop controlled by the
       size of |in_expr|, and containing |body| @>
  else loop.reset(make_for_loop@|
    (f.flags.to_ulong(),f.id,std::move(in_expr),std::move(body),which));
@/return type==void_type ? expression_ptr(new voiding(std::move(loop)))
  : @| conv!=nullptr ? expression_ptr(new conversion(*conv,std::move(loop)))
  : @| std::move(loop) ;
}

@ The |in_type| must be indexable by integers (so it is either a row-type or
vector, rational vector, matrix or string), or it must be the polynomial type,
indexable by |param_type|. If so, the first or second call to
|subscr_base::index_kind| will set |comp_type| to the component type resulting
from such a subscription. We also make |tp| point to the index type used.

@< Set |which| according to |in_type|, and set |bind| according to the
   identifiers contained in |f.id| @>=
{ type_expr comp_type; const_type_p inx_type;
  which = subscr_base::index_kind(in_type,*(inx_type=&int_type),comp_type);
  if (which==subscr_base::not_so)
    // if not integer-indexable, try parameter-indexable
    which = subscr_base::index_kind(in_type,*(inx_type=&param_type),comp_type);
  if (which==subscr_base::not_so) // if its not that either, it is wrong
  { std::ostringstream o;
    o << "Cannot iterate over value of type " << in_type;
    throw expr_error(e,o.str());
  }
  type_list it_comps; // type of ``iterator'' value (pattern) named in the loop
  it_comps.push_front(std::move(comp_type));
  it_comps.push_front(type_expr(inx_type->copy()));
  type_expr it_type(std::move(it_comps));
    // build tuple type from index and component types
  if (not pattern_type(f.id).specialise(it_type))
    throw expr_error(e,"Improper structure of loop variable pattern");
  thread_bindings(f.id,it_type,bind,true); // force all identifiers constant
}

@ Now follows the code that actually implements various kinds of loops. It is
templated both by |flags| indicating one of the $4$ reversal variants, and by
|kind| indicating the kind of value looped over (row, vector, rational
vector, string, matrix, parameter polynomial). This code is therefore compiled
in many variations, and we hope that constant folding and elimination of dead
code will reduce each instance to its specialised essence, which avoids
making the choices such as the |switch| below at evaluation time. We shall try
to reuse as much common source code among the variants as possible. This
presentation seems to be preferable and easier to maintain than providing $24$
explicit hand-specialised template instances.

We can start evaluating the |in_part| regardless of |kind|, but for deducing
the number of iterations we must already distinguish on |kind| to predict the
type of the in-part. A |loop_var| pair of values is constructed that will
temporarily contain the values of the loop index and loop component before
they are transferred to the execution context. The latter will happen by
calling the method |frame::bind| which explains why a shared pointer is formed
here. However because the corresponding pattern is limited, we know the shared
pointer itself will not be copied to the context, which explains why this pair
has to be allocated just once for all iterations of the loop. However, the
shared pointers that might get copied to the frame must be different on each
iteration, because the body might get hold, through a closure evaluated in the
loop body, of a copy of those pointers. In this case it is undesirable
that the closure should get to ``see'' changing values of variables during
subsequent iterations: in practice the closures would upon execution (after the
loop has run) each just see the final value of the variables rather than that
of the iteration in which the closure was formed.

@< Function definitions @>=
template <unsigned flags, subscr_base::sub_type kind>
void for_expression<flags,kind>::evaluate(level l) const
{ in_part->eval();
  own_tuple loop_var = std::make_shared<tuple_value>(2);
       // safe to re-use among iterations
  own_row result;
  std::vector<shared_value>::iterator dst;
  // set to |result->val.begin()| or |result->val.end()|
  try
  { switch (kind)
    @/{@;
      @< Cases for evaluating a loop over components of a value,
         each setting |result| @>
    }
  }
  catch (loop_break& err)
  { if (err.depth-- > 0)
      throw;
    if (l!=no_value)
    { if (out_reversed(flags))
        // doing |break| in reverse-gathering loop requires a shift
        dst=std::move(dst,result->val.end(),result->val.begin());
          // left-align |result|
      result->val.resize(dst-result->val.begin());
      // resize needed in all cases
    }
  }

  if (l!=no_value)
    push_value(std::move(result));
}

@ For evaluating |for| loops we must take care to interpret the |kind| field
when selecting a component from the in-part. Because of differences in the
type of |in_val|, some code must be duplicated, which we do as much as
possible by sharing modules between the various loop bodies textually.

@< Cases for evaluating a loop over components of a value... @>=
case subscr_base::row_entry:
  { shared_row in_val = get<row_value>();
    size_t n=in_val->val.size();
    @< Define loop index |i|, allocate |result| and initialise iterator |dst| @>
    while (i!=(in_forward(flags) ? n : 0))
    { loop_var->val[1]=in_val->val[in_forward(flags) ? i : i-1];
        // move index into |loop_var| pair
      @< Set |loop_var->val[0]| to |i++| or to |--i|, create a new |frame| for
      |pattern| binding |loop_var|, and evaluate the |loop_body| in it;
      maybe assign |*dst++| or |*--dst| from it @>
    }
  }
  @+break;
case subscr_base::vector_entry:
  { shared_vector in_val = get<vector_value>();
    size_t n=in_val->val.size();
    @< Define loop index |i|, allocate |result| and initialise iterator |dst| @>
    while (i!=(in_forward(flags) ? n : 0))
    { loop_var->val[1] = std::make_shared<int_value>
        (in_val->val[in_forward(flags) ? i : i-1]);
      @< Set |loop_var->val[0]| to... @>
    }
  }
  @+break;
case subscr_base::ratvec_entry:
  { shared_rational_vector in_val = get<rational_vector_value>();
    size_t n=in_val->val.size();
    @< Define loop index |i|, allocate |result| and initialise iterator |dst| @>
    while (i!=(in_forward(flags) ? n : 0))
    { loop_var->val[1] = std::make_shared<rat_value> @|
      (Rational
        (in_val->val.numerator()[in_forward(flags) ? i : i-1]
        ,in_val->val.denominator()));
      @< Set |loop_var->val[0]| to... @>
    }
  }
  @+break;
case subscr_base::string_char:
  { shared_string in_val = get<string_value>();
    size_t n=in_val->val.size();
    @< Define loop index |i|, allocate |result| and initialise iterator |dst| @>
    while (i!=(in_forward(flags) ? n : 0))
    { loop_var->val[1] = std::make_shared<string_value>
            (in_val->val.substr(in_forward(flags) ? i : i-1,1));
      @< Set |loop_var->val[0]| to... @>
    }
  }
  @+break;

@ Here are the remaining cases. The case |matrix_column| is ever so slightly
different because the iteration count~|n| is given by the number of columns,
rather than the size, of |inv_val->val|. The case |mod_poly_term| has more
important differences, to be detailed later. The other two cases should never
arise.

@< Cases for evaluating a loop over components of a value... @>=
case subscr_base::matrix_column:
  { shared_matrix in_val = get<matrix_value>();
    size_t n=in_val->val.numColumns();
    @< Define loop index |i|, allocate |result| and initialise iterator |dst| @>
    while (i!=(in_forward(flags) ? n : 0))
    { loop_var->val[1] = std::make_shared<vector_value>
        (in_val->val.column(in_forward(flags) ? i : i-1));
      @< Set |loop_var->val[0]| to... @>
    }
  }
  @+break;
case subscr_base::mod_poly_term:
  @< Perform a loop over the terms of a virtual module @>
  break;
case subscr_base::matrix_entry:; // excluded in type analysis
case subscr_base::not_so: assert(false);

@ The following code, which occurs five times, used both the input and output
direction attributes.

@< Define loop index |i|, allocate |result| and initialise iterator |dst| @>=
size_t i= in_forward(flags) ? 0 : n;
if (l!=no_value)
{ result = std::make_shared<row_value>(n);
  dst = out_forward(flags) ? result->val.begin() : result->val.end();
}

@ This code too occurs identically five times. We set the in-part component
stored in |loop_var->val[1]| separately for the various values of |kind|, but
|loop_var->val[0]| is always the (integral) loop index. Once initialised,
|loop_var| is passed by the method |frame::bind| through the function
|thread_components| to set up |loop_frame|, whose constructor has pushed it
onto |frame::current| to form the new evaluation context. Like for
|loop_var->val[0]|, it is important that |frame::current| be set to point to a
newly created frame at each iteration, since any closure values in the loop
body will incorporate its current instance by reference; there would be no
point in supplying fresh pointers in |loop_var| if they were subsequently
copied to overwrite the pointers in the same |evaluation_context| object each
time. Once these things have been handled, the evaluation of the loop body is
standard.

@< Set |loop_var->val[0]| to... @>=
{ loop_var->val[0] = std::make_shared<int_value>(in_forward(flags) ? i++ : --i);
    // create a fresh index each time
  frame loop_frame (pattern);
  loop_frame.bind(loop_var);
  if (l==no_value)
    body->void_eval();
  else
  {@; body->eval();
     *(out_forward(flags) ? dst++ : --dst) = pop_value();
  }
} // restore context upon destruction of |loop_frame|

@ The loop over terms of a virtual module is slightly different, and since it
handles values defined in the modules \.{atlas-types.w} we shall include its
header file. We implement the $4$ reversal variants, even though reversal at
the source makes little sense unless the internal order of the terms in the
polynomial (over which the user has no control) are meaningful to the user.
This reversal is implemented by using reverse iterators to control the loop;
the loop body itself is textually identical, though the type of |it| differs
between them.

@h "atlas-types.h"
@< Perform a loop over the terms of a virtual module @>=
{ shared_virtual_module pol_val = get<virtual_module_value>();
  size_t n=pol_val->val.size();
  if (l!=no_value)
  { result = std::make_shared<row_value>(n);
    dst = out_forward(flags) ? result->val.begin() : result->val.end();
  }
  if (in_forward(flags))
    for (auto it=pol_val->val.cbegin(); it!=pol_val->val.cend(); ++it)
      @< Loop body for iterating over terms of a virtual module @>
  else
    for (auto it=pol_val->val.crbegin(); it!=pol_val->val.crend(); ++it)
      @< Loop body for iterating over terms of a virtual module @>
}

@~And here is that loop body, included twice identically.

@< Loop body for iterating over terms of a virtual module @>=
{ loop_var->val[0] =
    std::make_shared<module_parameter_value>(pol_val->rf,it->first);
  loop_var->val[1] = std::make_shared<split_int_value>(it->second);
  frame loop_frame(pattern);
  loop_frame.bind(loop_var);
  if (l==no_value)
    body->void_eval();
  else
  {@; body->eval();
    *(out_forward(flags) ? dst++ : --dst) = pop_value();
  }
} // restore context upon destruction of |loop_frame|

@*1 Counted loops.
%
Next we consider counted |for| loops. Like with other loops there are quite a
few variations. For efficiency we shall handle cases of omitted identifier and
(lower) |bound=0| (most likely due to an omitted \&{from} clause) specially.
We could do the all distinctions detectable at compile time through the
template argument |flags|, which would avoid any runtime tests, but give a lot
of cases. We choose to represent in |flags| all distinction \emph{except} that
of an absent |bound| expression. The latter will be tested for presence when
initialising the lower bound; this gives a minute runtime cost when the bound
is present, but halves the number of template instances used.

@< Type def... @>=
template <unsigned flags>
struct counted_for_expression : public expression_base
{ id_type id; // may be $-1$, if |no_frame(flags)|
  expression_ptr count, bound, body;
  // we allow |bound| (but not |count|) to hold |nullptr|
@)
  counted_for_expression
@| (id_type i, expression_ptr&& cnt, expression_ptr&& bnd,
    expression_ptr&& b)
@/: id(i), count(cnt.release()),bound(bnd.release()), body(b.release())
  @+{}
  virtual ~@[counted_for_expression() nothing_new_here@];
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ Printing a counted |for| expression is straightforward, omitting optional
parts if absent.

@< Function definitions @>=
template <unsigned flags>
void counted_for_expression<flags>::print(std::ostream& out) const
{ if (has_frame(flags))
    out << " for " << main_hash_table->name_of(id) << ": " << *count;
  else out << " for : " << *count;  // omit nonexistent identifier
  if (bound.get()!=nullptr)
    out << " from " << *bound;
  print_body(out,body,flags);
}

@ As in |make_slice| and |make_for_loop| above, we need to convert runtime
values for |flags|. This time there are $6$ template instances, since a loop
index may be omitted (bit position $2$) but in that case it makes no sense to
reverse the order of traversal of the ``input sequence'' (and the syntax will
not allow it).

@< Local function definitions @>=
expression make_counted_loop (unsigned flags, id_type id, @|
   expression_ptr&& count, expression_ptr&& bound, expression_ptr&& body)
{ switch(flags)
  {
  case 0: return new @| counted_for_expression<0>
    (id,std::move(count),std::move(bound),std::move(body));
  case 1: return new @| counted_for_expression<1>
    (id,std::move(count),std::move(bound),std::move(body));
  case 2: return new @| counted_for_expression<2>
    (id,std::move(count),std::move(bound),std::move(body));
  case 3: return new @| counted_for_expression<3>
    (id,std::move(count),std::move(bound),std::move(body));
   case 4: return new @| counted_for_expression<4>
    (id,std::move(count),std::move(bound),std::move(body));
  case 6: return new @| counted_for_expression<6>
    (id,std::move(count),std::move(bound),std::move(body));
  default: assert(false); return(nullptr);
  }
}

@ Here is another way we can obtain a counted |for|-loop, from a
loop-over-components that does not bind any identifiers at all; a case that we
set aside as explained earlier.

In order to substitute a counted loop without index for such a loop, we need
to assemble a call to the appropriate size-computing built-in function. The
needed |shared_builtin| values will be assembled later, and some of the actual
built-in functions needed to be declared as global functions (and in one case
even defined in the first place) specifically for this purpose.

@< Set |loop| to a index-less counted |for| loop... @>=
{ expression_ptr call; const source_location &loc = f.in_part.loc;
  switch(which)
  {
  case subscr_base::row_entry: call.reset(new @|
      overloaded_builtin_call(sizeof_row_builtin,std::move(in_expr),loc));
  break; case subscr_base::vector_entry: call.reset(new @|
      overloaded_builtin_call(sizeof_vector_builtin,std::move(in_expr),loc));
  break; case subscr_base::ratvec_entry: call.reset(new @|
      overloaded_builtin_call(sizeof_ratvec_builtin,std::move(in_expr),loc));
  break; case subscr_base::string_char: call.reset(new @|
      overloaded_builtin_call(sizeof_string_builtin,std::move(in_expr),loc));
  break; case subscr_base::matrix_column: call.reset(new @|
      overloaded_builtin_call(matrix_columns_builtin,std::move(in_expr),loc));
  break; case subscr_base::mod_poly_term: call.reset(new @|
      overloaded_builtin_call(sizeof_parampol_builtin,std::move(in_expr),loc));
  break; default: assert(false);
  }

  if (f.flags.test(1)) // whether reversed assembly of return value
    loop.reset(new @| counted_for_expression<6>
      (-1,std::move(call),nullptr,std::move(body)));
  else
    loop.reset(new @| counted_for_expression<4>
      (-1,std::move(call),nullptr,std::move(body)));
}

@ Type-checking counted |for| loops is rather like that of other |for| loops,
but we must extend the context with the loop variable while processing the loop
body.

@< Cases for type-checking and converting... @>=
case cfor_expr:
{ const cfor_node& c=*e.cfor_variant;
  expression_ptr count_expr = convert_expr(c.count,as_lvalue(int_type.copy()));
  expression_ptr bound_expr = is_empty(c.bound)
    ? nullptr
    : convert_expr(c.bound,as_lvalue(int_type.copy())) ;
@)
  type_expr body_type;
  type_expr *btp=&body_type; // point to place to record body type
  const conversion_record* conv=nullptr;
  if (type==void_type)
    btp=&type; // we can reuse this type; no risk of specialisation
  else if (type.specialise(row_of_type))
    btp=type.component_type;
  else if ((conv=row_coercion(type,body_type))==nullptr)
    throw type_error(e,row_of_type.copy(),std::move(type));
@)
  if (c.flags[2]) // case of absent loop variable
  { layer bind(0,nullptr);  // no local variables for loop, but allow |break|
    expression_ptr body(convert_expr (c.body,*btp));
    if (type!=void_type and *btp==void_type and not is_empty(c.body))
      body.reset(new voiding(std::move(body)));
  @/expression_ptr loop(make_counted_loop(c.flags.to_ulong(), @|
      c.id,std::move(count_expr),std::move(bound_expr),std::move(body)));
  return type==void_type ? expression_ptr(new voiding(std::move(loop))) : @|
         conv!=nullptr ? expression_ptr(new conversion(*conv,std::move(loop)))
                       : @| std::move(loop);
  }
  else // case of a present loop variable
  { layer bind(1,nullptr);
    bind.add(c.id,int_type.copy(),true); // add |id| as constant
    expression_ptr body(convert_expr (c.body,*btp));
    if (type!=void_type and *btp==void_type and not is_empty(c.body))
      body.reset(new voiding(std::move(body)));
  @/expression_ptr loop(make_counted_loop(c.flags.to_ulong(), @|
      c.id,std::move(count_expr),std::move(bound_expr),std::move(body)));
    return type==void_type ? expression_ptr(new voiding(std::move(loop))) : @|
         conv!=nullptr ? expression_ptr(new conversion(*conv,std::move(loop)))
                       : @| std::move(loop);
  }
}

@ Executing a loop is a simple variation of what we have seen before for
|while| and |for| loops over value components. We distinguish a number of
cases, some by template argument and others dynamically, but always before
entering the loop, in order to allow an optimal adaptation to the task. We en
up always using a \Cpp\ |while| loop to implement the |for| loop.

A special case that is optimised for is when the user gave no name to loop
index, so bit |flags&0x4| is set: we can then omit introducing a |frame| for
the loop altogether. Also, since the syntax ensures that absence of a name for
the loop index implies absence of a lower bound expression (which would be
unused anyway) we can omit trying to evaluate |bound| here. We choose to
always use a decreasing loop counter internally in such cases, as this
simplifies the termination condition and might therefore be marginally faster.

@< Function definitions @>=
template <unsigned flags>
void counted_for_expression<flags>::evaluate(level l) const
{ int c=(count->eval(),get<int_value>()->val);
  if (c<0)
    c=0; // no negative size result

  if (has_frame(flags)) // then loop uses index
  { int b=(bound.get()==nullptr ? 0 : (bound->eval(),get<int_value>()->val));
    c+=b; // set to upper bound, exclusive
    id_pat pattern(id);
    if (l==no_value)
      @< Perform counted loop that uses an index, without storing result,
         between lower bound |b| and exclusive upper bound |c| @>
    else
      @< Perform counted loop that uses an index, pushing result to
         |execution_stack|,
         between lower bound |b| and exclusive upper bound |c| @>
  }
  else if (l==no_value) // counted loop without index and no value
  { try {@; while (c-->0)
      body->void_eval();
    }
    catch (loop_break& err) @+
    {@; if (err.depth-- > 0)
          throw;
    }
  }
  else // counted loop without index producing a value
  { own_row result = std::make_shared<row_value>(c);
    auto dst = out_forward(flags) ? result->val.begin() : result->val.end();
    try @/{@;
      while (c-->0)
      {@; body->eval();
        *(out_forward(flags)? dst++:--dst) = pop_value();
      }
    }
    catch (loop_break& err)
    { if (err.depth-- > 0)
        throw;
      if (out_reversed(flags))
        dst=std::move(dst,result->val.end(),result->val.begin());
        // after break, left-align |result|
      result->val.resize(dst-result->val.begin());
    }
    push_value(std::move(result));
  }

}

@ For a counted loop using its index, we remain at the \Cpp\ level close to
the \.{axis} loop being implemented, but we transform the bound |b| or|c| into
our loop index. For decreasing loops we distinguish the case where the lower
bound is~$0$, the default value, since a test against a constant~$0$ is a bit
more efficient.

@< Perform counted loop that uses an index, without storing result,
   between lower bound |b| and exclusive upper bound |c| @>=
{ try
  { if (in_forward(flags)) // increasing loop
      while (b<c)
      @/{@; frame fr(pattern);
        fr.bind(std::make_shared<int_value>(b++));
        body->void_eval();
      }
    else if (b!=0)
      while (c-->b)
      @/{@; frame fr(pattern);
        fr.bind(std::make_shared<int_value>(c));
        body->void_eval();
      }
    else // same with |b==0|, but this is marginally faster
       while (c-->0)
      @/{@; frame fr(pattern);
        fr.bind(std::make_shared<int_value>(c));
        body->void_eval();
      }
  }
  catch (loop_break& err) @+
  {@; if (err.depth-- > 0)
        throw;
  }
}

@ The case where a result is accumulated just differs by adding the necessary
bits of stuff.

@< Perform counted loop that uses an index, pushing result to |execution_stack|,
   between lower bound |b| and exclusive upper bound |c| @>=
{ own_row result = std::make_shared<row_value>(c);
  auto dst = out_forward(flags) ? result->val.begin() : result->val.end();
  try
  { if (in_forward(flags)) // increasing loop
      while (b<c)
      { frame fr(pattern);
        fr.bind(std::make_shared<int_value>(b++));
        body->eval();
        *(out_forward(flags)? dst++:--dst) = pop_value();
      }
    else if (b!=0)
      while (c-->b)
      { frame fr(pattern);
        fr.bind(std::make_shared<int_value>(c));
        body->eval();
        *(out_forward(flags)? dst++:--dst) = pop_value();
      }
    else // same with |b==0|, but this is marginally faster
      while (c-->0)
      { frame fr(pattern);
        fr.bind(std::make_shared<int_value>(c));
        body->eval();
        *(out_forward(flags)? dst++:--dst) = pop_value();
     }
  }
  catch (loop_break& err)
  { if (err.depth-- > 0)
      throw;
    if (out_reversed(flags))
      dst=std::move(dst,result->val.end(),result->val.begin());
      // after break, left-align |result|
    result->val.resize(dst-result->val.begin());
  }
  push_value(std::move(result));
}

@* Casts and operator casts.
%
Casts are very simple to process. They do not need any |expression| type to
represent them, so type-checking is all there is to it. Nonetheless there is a
subtlety in the code below, which starts with |type.specialise(c.type)| even
though the final |conform_types| will attempt the same specialisation. One
consequence is that if the cast should (rather sillily) specify a less
specific type than |type| from the context, as the second cast
in \.{[rat]:[*]:[0]}, then the stronger context will be used in the
conversion. More importantly though, the specialisation now takes
place \emph{before} the conversion, which in case our cast is a function body
means that the result type of that function gets specialised before the rest
of the body is analysed, and any |return| expressions in the body will profit
from the cast's strong type context.

@< Cases for type-checking and converting... @>=
case cast_expr:
{ cast_node& c=*e.cast_variant;
  if (type.specialise(c.type)) // see if we can do without conversion
    return convert_expr(c.exp,type); // in which case use now specialised |type|
  expression_ptr p = convert_expr(c.exp,c.type); // otherwise use |c.type|
  return conform_types(c.type,type,std::move(p),e);
}

@ Another kind of cast is the operator cast, which selects an operator or
overloaded function instance as it would for arguments of specified types, but
without giving actual arguments, so that the selected function itself can be
handled as a value. In order to do so we shall need the following function.

The overload table stores type information in a |func_type| value, which
cannot be handed directly to the |specialise| method. The following function
simulates specialisation to a function type |from|$\to$|to|.

@< Local function definitions @>=
inline bool functype_specialise
  (type_expr& t, const type_expr& from, const type_expr& to)
{ return t.specialise(gen_func_type) @|
  and t.func->arg_type.specialise(from) @|
  and t.func->result_type.specialise(to);
}

@ Operator casts only access already existing values. In most cases we must
access the global overload table to find the value. Since upon success we find
a bare |shared_function|, we must (as we did for~`\.\$') use the
|capture_expression| class to serve as wrapper that upon evaluation will
return the value again.

@< Cases for type-checking and converting... @>=
case op_cast_expr:
{ const op_cast& c=e.op_cast_variant;
  const overload_table::variant_list& variants =
   global_overload_table->variants(c->oper);
  const type_expr& ctype=c->type;
  std::ostringstream o;
  o << main_hash_table->name_of(c->oper) << '@@' << ctype;
@)
  size_t i;
  for (i=0; i<variants.size(); ++i)
    if (variants[i].type().arg_type==ctype)
      break;
  if (i<variants.size()) // something was found
  {
    expression_ptr p(new capture_expression(variants[i].val,o.str()));
    const type_expr& res_t = variants[i].type().result_type;
    if (functype_specialise(type,ctype,res_t) or type==void_type)
      return p;
    throw type_error(e,type_expr(ctype.copy(),res_t.copy()),std::move(type));
  }
@)// now we have no match from the overload table, try generic operations
  if (is_special_operator(c->oper))
    @< Test special argument patterns, and on match |return| an appropriate
       denotation @>
@/throw program_error("No instance for "+o.str()+" found");
}
break;

@ For our special operators, |print|, |prints|, |to_string|, |error|, we
select their wrapper function here always, since they accept any argument
type. We signal an error only if the context requires a type that cannot be
specialised to the type of operator found. For $\#$ the situation will be
slightly more complicated.

@< Test special argument patterns... @>=
{ if (c->oper==print_name())
  { if (functype_specialise(type,ctype,ctype))
    return expression_ptr(new @| capture_expression (print_builtin,o.str()));
  }
  else if (c->oper==prints_name())
  { if (functype_specialise(type,ctype,void_type))
    return expression_ptr(new @| capture_expression (prints_builtin,o.str()));
  }
  else if (c->oper==to_string_name())
  { if (functype_specialise(type,ctype,str_type))
    return expression_ptr(new @| capture_expression (to_string_builtin,o.str()));
  }
  else if (c->oper==error_name())
  { if (functype_specialise(type,ctype,unknown_type))
    return expression_ptr(new @| capture_expression (error_builtin,o.str()));
  }
  else if (c->oper==size_of_name())
    @< Select the proper instance of the \.\# operator,
       or fall through if none applies @>
  else if (c->oper==concatenate_name())
    @< Select the proper instance of the \.{\#\#} operator,
       or fall through if none applies @>
}

@ For the \.\# operator, we select from four possible variants that deliver
different wrapper functions. We signal an error if we found a match but the
type of the resulting operator does not match the type required by the
context. If no match is found here, there can still be one in the overload
table.

@< Select the proper instance of the \.\# operator,... @>=
{ if (ctype.kind==row_type)
  { if (functype_specialise(type,ctype,int_type))
  @/return expression_ptr(new @|
      capture_expression (sizeof_row_builtin,o.str()));
    throw type_error(e,ctype.copy(),std::move(type));
  }
  else if (is_pair_type(ctype))
  {
    type_expr& arg_tp0 = ctype.tupple->contents;
    type_expr& arg_tp1 = ctype.tupple->next->contents;
    if (arg_tp0.kind==row_type and *arg_tp0.component_type==arg_tp1)
    { if (functype_specialise(type,ctype,arg_tp0))
        return expression_ptr(new @|
          capture_expression(suffix_elt_builtin,o.str()));
      throw type_error(e,ctype.copy(),std::move(type));
    }
    if (arg_tp1.kind==row_type and *arg_tp1.component_type==arg_tp0)
    { if (functype_specialise(type,ctype,arg_tp1))
      return expression_ptr(new @|
        capture_expression (prefix_elt_builtin,o.str()));
      throw type_error(e,ctype.copy(),std::move(type));
    }
  }
}

@ Like for \.\#, we select for the \.{\#\#} operator the variants that deliver
different wrapper functions.

@< Select the proper instance of the \.{\#\#} operator,... @>=
{ if (ctype.kind==row_type and ctype.component_type->kind==row_type)
  { if (functype_specialise(type,ctype,*ctype.component_type))
  @/return expression_ptr(new @|
      capture_expression (join_rows_row_builtin,o.str()));
    throw type_error(e,ctype.copy(),std::move(type));
  }
  else if (is_pair_type(ctype))
  {
    type_expr& arg_tp0 = ctype.tupple->contents;
    type_expr& arg_tp1 = ctype.tupple->next->contents;
    if (arg_tp0.kind==row_type and arg_tp1==arg_tp0)
    { if (functype_specialise(type,ctype,arg_tp0))
        return expression_ptr(new @|
          capture_expression (join_rows_builtin,o.str()));
      throw type_error(e,ctype.copy(),std::move(type));
    }
  }
}

@* Assignments.
%
Syntactically there is hardly anything simpler than simple assignment
statements. However, semantically we distinguish assignments to local and to
global variables. Then there are ``component assignments'' which modify
repetitive values like row values by changing one component; these too will
distinguish local and global versions. Finally, while not present in the
initial language design, a multiple assignment statement was added to the
language that can take apart tuple components, just as can be done in
definitions of new (global or local) variables. As these can mix local and
global destinations, we will introduce a separate expression type for these
somewhat more expensive assignments. We retain the following intermediate
class for all assignment statements except multiple assignments (the latter
will derive directly from |expression_base|), as it allows to avoid a bit of
code duplication.

@< Type definitions @>=
struct assignment_expr : public expression_base
{ id_type lhs;
  expression_ptr rhs;
@)
  assignment_expr(id_type l,expression_ptr&& r)
   : lhs(l),rhs(r.release()) @+{}
  virtual ~@[assignment_expr() nothing_new_here@];
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
  global_assignment(id_type l,expression_ptr&& r);
  virtual ~@[global_assignment() nothing_new_here@];
  virtual void evaluate(level l) const;
};

@ The constructor for |global_assignment| is rather similar to that for
|global_identifier|.

@< Function def... @>=
global_assignment::global_assignment(id_type l,expression_ptr&& r)
: assignment_expr(l,std::move(r)), address(global_id_table->address_of(l)) @+{}

@ Evaluating a global assignment evaluates the left hand side, and replaces
the old value stored at |*address| by the new (shared pointer) value. The
value is then also pushed on the evaluation stack according to the level |l|
(which will however often be |no_value| in which case nothing happens).

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
coordinates of the identifier in the execution context.

@< Type definitions @>=
class local_assignment : public assignment_expr
{ size_t depth, offset;
public:
  local_assignment(id_type l, size_t i,size_t j, expression_ptr&& r);
  virtual ~@[local_assignment() nothing_new_here@];
  virtual void evaluate(level l) const;
};

@ The constructor for |local_assignment| is straightforward.

@< Function def... @>=
local_assignment::local_assignment
 (id_type l, size_t i,size_t j, expression_ptr&& r)
: assignment_expr(l,std::move(r)), depth(i), offset(j) @+{}

@ Evaluating a local assignment evaluates the left hand side, and replaces the
old value stored at |frame::current->elem(depth,offset)| by the new (shared
pointer) value.

@< Function def... @>=
void local_assignment::evaluate(level l) const
{ rhs->eval();
  shared_value& dest =  frame::current->elem(depth,offset);
  dest= pop_value();
  push_expanded(l,dest);
}

@ The type for multiple assignments has to cater for a mixture of global and
local names present in the destination pattern. This is don by having
(possibly empty) vectors for both types of destination, and a |Bitmap| telling
for each name in left-to-right order whether it is global.

The constructor here is more elaborate than for simple assignments, because
the various identifiers in the pattern |lhs| must be looked up, the fields of
the class initialised to reflect the results (or an error could be thrown),
and a corresponding type needs to be exported to serve as a cast for the right
hand side expression.

@h "bitmap.h"

@< Type definitions @>=
class multiple_assignment : public expression_base
{
 public:
  struct local_dest {@; size_t depth, offset; };
  typedef containers::simple_list<local_dest> loc_list;
  typedef containers::simple_list<shared_share> glob_list;
 private:
  id_pat lhs;
  expression_ptr rhs;
  loc_list locals;
  glob_list globals;
  BitMap is_global;
public:
  multiple_assignment
    (const id_pat& lhs,expression_ptr&& r
    ,loc_list&& ll, glob_list&& gl, BitMap&& bm);
  virtual ~@[multiple_assignment() nothing_new_here@];
  virtual void print(std::ostream& out) const;
  virtual void evaluate(level l) const;
};

@ The constructor needs to be defined outside the class definition because it
sues |copy_id_pat|.
@< Function definitions @>=
multiple_assignment::multiple_assignment @|
    (const id_pat& lhs,expression_ptr&& r
    ,loc_list&& ll, glob_list&& gl, BitMap&& bm)
: lhs(copy_id_pat(lhs)), rhs(std::move(r))
, locals(std::move(ll)),globals(std::move(gl)), is_global(std::move(bm))@+{}

@ Printing reflects the special syntax for multiple assignments.
@< Function def... @>=
void multiple_assignment::print(std::ostream& out) const
{@; out << "set " << lhs << ":=" << *rhs; }

@ Evaluating a multiple assignment requires in general a recursive traversal
of the left hand side pattern, with a corresponding traversal of the right
hand side value, to assign to each destination identifier the corresponding
component of that value. Traversal is in post-order, which is motivated by
considerations of type-checking as explained later. Preparing the required
destinations in order is sufficiently complicated to merit the definition of a
special class |dest_iterator| dedicated to this. It acts as an output
iterator, whose sole purpose is to accept the stream of values produced. But
rather than using the composite syntax |*dst++=v| we write |dst.receive(v)|,
the method |receive| doing the required assignment and iterator advancing.

@< Local function definitions @>=
void thread_assign
  (const id_pat& pat, const shared_value& val, dest_iterator& dst)
{ if ((pat.kind&0x2)!=0)
  { const tuple_value* t=force<tuple_value>(val.get());
    assert(t->val.size()==length(pat.sublist));
    auto src = t->val.begin();
    for (auto it=pat.sublist.begin(); not pat.sublist.at_end(it); ++it,++src)
      thread_assign(*it,*src,dst);
  }
  if ((pat.kind&0x1)!=0)
    dst.receive(val);
}

@ It is now easy to write the |multiple_assignment::evaluate| method. The
|dest_iterator| will need to know about the internal lists |locals,globals|
and bitmap |is_global| to be able to do its job.

@< Function definitions @>=
void multiple_assignment::evaluate (level l) const
{ rhs->eval();
  shared_value v = pop_value();
  dest_iterator dst(locals,globals,is_global);
  thread_assign(lhs,v,dst);
  push_expanded(l,v);
}

@ So now we are left with describing the output iterator. Internally it
suffices to store weak iterators into the lists of local and global
destinations of the assignment, and a current position in the bitmap.

@< Local class definitions @>=
class dest_iterator
{
  multiple_assignment::loc_list::weak_const_iterator l_it;
  multiple_assignment::glob_list::weak_const_iterator g_it;
  const BitMap& is_global; @+
  unsigned long n; // index into |is_global|
 public:
  dest_iterator @|
   (const multiple_assignment::loc_list& locs
   ,const multiple_assignment::glob_list& globs
   ,const BitMap& is_global)
  : l_it(locs.wcbegin())
  , g_it(globs.wcbegin())
  , is_global(is_global)
  , n(0)
  @+{}
  void receive (const shared_value& val);
};

@ And finally here is the |receive| method. With the proper preparations, it
is not really difficult to implement, and easily understood by comparing to
statements in the |evaluate| methods for simple assignments.

@< Local function definitions @>=
void dest_iterator::receive (const shared_value& val)
{ if (is_global.isMember(n++))
  @/{@; assert (not g_it.at_end());
    *(*g_it++) = val;
  } // send to |shared_share| stored in |globs|
  else
  @/{@;
    assert (not l_it.at_end());
    frame::current->elem(l_it->depth,l_it->offset)=val;
    ++l_it;
  }
}

@ Here are some simple functions that will be called for errors both in simple
assignments and component assignments.

@< Local function definitions @>=
void report_undefined (id_type id,const expr& e,const char* where)
{ std::ostringstream o;
  o << "Undefined identifier '" << main_hash_table->name_of(id)
  @|<< "' in " << where << ' ' << e;
  if (e.loc.file!=Hash_table::empty)
    o << ' ' << e.loc;
  throw program_error (o.str());
}
@)
void report_constant_modified (id_type id,const expr& e,const char* where)
{ std::ostringstream o;
  o << "Name '" << main_hash_table->name_of(id)
    << "' is constant in " << where << ' ' << e;
  if (e.loc.file!=Hash_table::empty)
    o << ' ' << e.loc;
  throw program_error (o.str());
}

@ Converting assignment statements follows the same lines as for applied
identifiers, as far as discriminating between local and global is concerned.
We first look in |id_context| for a local binding of the identifier, and then
maybe in |global_id_table|. If found in either way, the right hand side is
converted in a type context given by the type of the variable found. After
forming the proper kind of assignment expression, we must as usual allow for a
coercion to be applied to the result of the assignment, if the required |type|
demands this.

While in most successful cases the type of the variable may direct the
conversion of the right hand side, the type of the right hand side may
occasionally be more specific than the previously known type of the variable
(only if it is a specialisation of the latter will the conversion succeed).
For instance this happens when assigning a row of concrete type to a variable
initialised with an empty row. In those cases we call the
|specialise| method of |frame| or of |global_id_table| to make sure the type
assumed by the variable is recorded.

Since variables of |void| type are allowed and can be assigned to (and due to
the voiding coercion a right hand side of any type will be accepted), we must
take care to insert a |voiding| in such rare cases, to ensure no value will be
computed and assigned.

@< Cases for type-checking and converting... @>=
case ass_stat:
if ( e.assign_variant->lhs.kind==0x1) // single identifier, do simple assign
{
  id_type lhs=e.assign_variant->lhs.name;
  const_type_p id_t; size_t i,j; bool is_const;
  const bool is_local = (id_t=layer::lookup(lhs,i,j,is_const))!=nullptr;
  if (not is_local and (id_t=global_id_table->type_of(lhs,is_const))==nullptr)
    report_undefined(lhs,e,"assignment");
@.Undefined identifier@>
  if (is_const)
    report_constant_modified(lhs,e,"assignment");
@.Name is constant @>
@)
  type_expr rhs_type = id_t->copy(); // provide a modifiable copy
  expression_ptr r(convert_expr(e.assign_variant->rhs,rhs_type));
  if (rhs_type!=*id_t)
    // assignment will specialise identifier, record to which type it does
  {@; if (is_local)
      layer::specialise(i,j,rhs_type);
      else
      global_id_table->specialise(lhs,rhs_type);
  }
  if (rhs_type==void_type and not is_empty(e.assign_variant->rhs))
    r.reset(new voiding(std::move(r)));
@)
  expression_ptr assign = is_local
  ? expression_ptr(new local_assignment(lhs,i,j,std::move(r)))
@/: expression_ptr(new global_assignment(lhs,std::move(r)));
  return conform_types(rhs_type,type,std::move(assign),e);
}
else @< Generate and |return| a |multiple_assignment| @>

@ For traversing the left hand side pattern in a multiple definition, we need
some semi-local variables, to be accessible from within the recursive function
but not renewed for each recursive call. The solution of passing around a
reference to a structure containing those variables is elegantly realised by
definition the traversal function as a recursive method of that structure (the
implicit reference |*this| is passed around unchanged).

@< Local class definitions @>=
struct threader
{ typedef containers::sl_list<multiple_assignment::local_dest> loc_list;
  typedef containers::sl_list<shared_share> glob_list;
@)
  const expr& e;
  loc_list locs;
  glob_list globs;
  BitMap is_global;
  containers::sl_list<std::pair<id_type,const_type_p> > assoc;
@)
  threader (const expr& e) : e(e), locs(), globs(), is_global(), assoc() @+{}
  void thread (const id_pat& pat,type_expr& type); // recursively analyse |pat|
  void refine () const; // maybe specialise stored identifiers
};

@ The left hand side pattern is traversed in post-order: when there is both an
identifier for the whole and a sub-list, the former is handled after the
latter. This simplifies testing of type compatibility in the destination
pattern, where the only possible error now is that a type for a ``parent''
identifier does not match the (possibly partly specified) tuple type
established by its children.

@< Function definitions @>=
void threader::thread(const id_pat& pat,type_expr& type)
{ if ((pat.kind&0x4)!=0)
    @< Throw an error to signal forbidden qualifier \.! before |pat.name| @>
  if ((pat.kind&0x2)!=0) // first treat any sublist
  { type.specialise(unknown_tuple(length(pat.sublist)));
    assert(type.kind==tuple_type); // this should succeed
    wtl_iterator t_it(type.tupple);
    for (auto it=pat.sublist.begin(); not pat.sublist.at_end(it); ++it,++t_it)
      thread(*it,*t_it);
  }
  if ((pat.kind&0x1)!=0)
  { id_type id = pat.name;
    const_type_p id_t; // will point to type of local or global |id|
    @< Check that |id| did not occur previously in this left hand side @>
    size_t i,j; bool is_const;
    const bool is_local = (id_t=layer::lookup(id,i,j,is_const))!=nullptr;
    if (not is_local and (id_t = global_id_table->type_of(id,is_const))==nullptr)
      report_undefined(id,e,"multiple assignment");
    if (is_const)
      report_constant_modified(id,e,"multiple assignment");
    is_global.extend_capacity(not is_local);
@)
    if (not type.specialise(*id_t)) // incorporate found type into |type|
      @< Throw an error to signal type incompatibility for |id| @>
    assoc.push_back(std::make_pair(id,&type));
      // record pointer to |type| for later refinement of |id|
    if (is_local)
      locs.push_back(@[multiple_assignment::local_dest{i,j}@]);
    else
      globs.push_back(global_id_table->address_of(id));
  }
}

@ The error signalled here should really be a syntax error, but the fact that
a generator without conflicts can be generated for our grammar depends on the
pattern after \.{set} being independent of whether \.= or \.{:=} follows it;
this is why we allowed these qualifiers to arrive up to this point.

@< Throw an error to signal forbidden qualifier \.! before |pat.name| @>=
{ std::ostringstream o;
  o << "Cannot constant-qualify '!' identifier '"
    << main_hash_table->name_of(pat.name) @| << "' in multi-assignment";
  throw expr_error(e,o.str());
}

@ Multiple assignment is only well defined if all target variables are
disjoint. Component assignments are not possible here, so we simply check that
all identifiers used are distinct.

@< Check that |id| did not occur previously in this left hand side @>=
for (auto it=assoc.wcbegin(); not it.at_end(); ++it)
  if (it->first==id)
   { std::ostringstream o;
     o << "Multiple assignments to same identifier '"
       << main_hash_table->name_of(pat.name) @| << "' in multi-assignment";
     throw expr_error(e,o.str());
   }


@ The error of incompatible left hand side patterns does not match the format
of a |type_error|, so we throw a more general |expr_error| instead, after
assembling the data to identify the error.

@< Throw an error to signal type incompatibility for |id| @>=
{ std::ostringstream o;
  o << "Incompatible type for '" << main_hash_table->name_of(id)
  @|<< "' in multi-assignment: type " << *id_t
  @|<< " does no match pattern " << type;
  throw expr_error(e,o.str());
}

@ Just like for simple assignments, there is a remote possibility that the
type of the right hand side specialises types of one or more variables used in
the left hand side pattern. This means that after converting, the type that
was found for that identifier will have been specialised. The method |refine|
traverses all identifiers and specialises their type, as stored either in
|global_id_table| or in |layer::lexical_context|; most of the time this will
do nothing, but when there is need, this will do what is required.

The implementation traverses the |assoc| list, and in parallel the bits in
|is_global| to determine whether a global or local identifier type has to be
specialised. In the case of global identifiers the identifier stored in
|assoc| is used, but for local identifiers the |depth| and |offset| stored in
|locs| are used instead, which requires a separate iterator |loc_it|.

@< Function definitions @>=
void threader::refine() const
{ unsigned long n=0;
  auto loc_it = locs.cbegin();
  for (auto it = assoc.cbegin(); not assoc.at_end(it); ++it)
    if (is_global.isMember(n++))
      global_id_table->specialise(it->first,*it->second);
    else
    {@;
      layer::specialise(loc_it->depth,loc_it->offset,*it->second);
      ++loc_it;
    }
}

@ For a multiple assignment, we first get information about type |rhs_type| of
the assignment from a call to |threader::thread| with the pattern~|pat| of the
left hand side (this may leave some slots undefined for components of the
value without destination variable in the assignment), then we convert the
right hand side in this context (which may fill in the open slots, though that
information will not be used), and finally call |refine| to adapt the types of
the left hand side variables to possible specialisations a variable type made
during of that conversion.

As before, we need to take care to insert a |voiding| in the rare case that
the type assigned is void (and the RHS is nonempty); this really can only
arise when |pat| specifies a $0$-tuple, which is a completely useless case.
Finally our |multiple_assignment| is constructed, mostly from information
stored in |thr|.

@< Generate and |return| a |multiple_assignment| @>=
{ const id_pat& pat=e.assign_variant->lhs;
  type_expr rhs_type;
  threader thr(e);
  thr.thread(pat,rhs_type);
  expression_ptr r = convert_expr(e.assign_variant->rhs,rhs_type);
  thr.refine();
  if (rhs_type==void_type and not is_empty(e.assign_variant->rhs))
    r.reset(new voiding(std::move(r)));
  expression_ptr m_ass (
    new @| multiple_assignment
      (pat,std::move(r)
      ,thr.locs.undress(),thr.globs.undress(),std::move(thr.is_global)));
  return conform_types(rhs_type,type,std::move(m_ass),e);
}

@*1 Component assignments.
%
The language we are implementing does not employ the notion of sub-object; in
other words if one sets $b=a[i]$ for some list, vector or matrix $a$, then $s$
will behave as a copy of the entry $a[i]$ rather than as an alias, so
subsequent assignment to $b$ will not affect~$a$ or vice versa. (This does no
prevent us to share storage between $b$ and $a$ initially, it just means the
sharing should be broken if $b$ or $a$ are modified; we practice
copy-on-write.) This simplifies the semantic model considerably, but if we
want to allow creating composite values by subsequently setting their
components, we need to allow assignments of the form $a[i]:=c$. The meaning of
this is the same as assigning a new value to all of $a$ that differs from the
original value only at index~$i$; it may however be expected to be implemented
more efficiently if the storage of $a$ is not currently shared, as would
usually be the case at least from the second such assignment to~$a$ on. The
interpreter will have to treat such component assignments as a whole (with
three components $a,i,c$, in which $a$ must be an identifier), which also
means that it will not be able to handle something like $a[i][j]:=c$ even when
that would seem to make sense (however $m[i,j]:=c$ for matrix values $m$ will
be supported).

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
template <bool reversed>
struct component_assignment : public assignment_expr
{ expression_ptr index;
@)
  component_assignment
   (id_type a,expression_ptr&& i,expression_ptr&& r)
   : assignment_expr(a,std::move(r)), index(i.release()) @+{}
  virtual ~@[component_assignment() nothing_new_here@];

  virtual void print (std::ostream& out) const;
@)
  void assign(level l,shared_value& aggregate,subscr_base::sub_type kind) const;
};

@ Printing reassembles the subexpressions according to the input syntax.
@< Function def...@>=
template <bool reversed>
void component_assignment<reversed>::print(std::ostream& out) const
{@; out << main_hash_table->name_of(lhs) << (reversed ? "~[" : "[")
        << *index << "]:=" << *rhs;
}

@ For global assignments, we need to have non-|const| access the location
where the identifier is stored.

@< Type definitions @>=
template <bool reversed>
class global_component_assignment : public component_assignment<reversed>
{ typedef component_assignment<reversed> base;
@)
  subscr_base::sub_type kind;
  shared_share address;
public:
  global_component_assignment
    (id_type a,expression_ptr&& i,expression_ptr&& r,
     subscr_base::sub_type k);
  virtual void evaluate(expression_base::level l) const;
};

@ The constructor for |global_component_assignment| stores the address of the
aggregate object and the component kind.

@< Function def... @>=
template <bool reversed>
global_component_assignment<reversed>::global_component_assignment @|
  (id_type a,expression_ptr&& i,expression_ptr&& r, subscr_base::sub_type k)
: base(a,std::move(i),std::move(r))
, kind(k),address(global_id_table->address_of(a)) @+{}

@ It is in evaluation that component assignments differ most from ordinary
ones. The work is delegated to the |assign| method of the base class, which is
given a reference to the |shared_value| pointer holding the current value of
the aggregate; it is this pointer that is in principle modified. Like when
fetching the value of a global variable, we must be aware of a possible
undefined value in the variable.

@< Function def... @>=
template <bool reversed>
void global_component_assignment<reversed>::evaluate(expression_base::level l)
  const
{ if (address->get()==nullptr)
  { std::ostringstream o;
    o << "Assigning to component of uninitialized variable "
      << main_hash_table->name_of(this->lhs);
    throw runtime_error(o.str());
  }
  base::assign(l,*address,kind);
}

@ The |assign| method, which will also be called for local component
assignments, starts by the common work of evaluating the (component) value to
be assigned, and of then making sure the aggregate variable is made to point
to a unique copy of its current value, which copy can then be modified in
place. The index is not yet evaluated at this point, but this will be done
inside the |switch| statement; this is because possible expansion of a tuple
index value depends on~|kind|. For actually changing the aggregate, we must
distinguish cases according to the kind of component assignment at hand.
Assignments to components of rational vectors and of strings will be
forbidden, see module @#comp_ass_type_check@>.

@< Function def... @>=
template <bool reversed>
void component_assignment<reversed>::assign
  (level lev,shared_value& aggregate, subscr_base::sub_type kind) const
{ rhs->eval();
  value loc=uniquify(aggregate);
    // simple pointer to modifiable value from shared pointer
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
  @/@< Replace column at |index| in matrix |loc| by value on stack @>
  @+break;
  default: {} // remaining cases are eliminated in type analysis
  }
}

@ A |row_value| component assignment is the simplest kind. The variable |loc|
holds a generic pointer, known to refer to a |row_value|. Since we need to
access the vector of shared pointers, we use |force| to get an ordinary
pointer, and then select the |val| field. Then we do a bound check, and on
success replace a component of the value held in |a| by the stack-top value.
Afterwards, depending on |l|, we may put back the stack-top value as result of
the component assignment, possibly expanding a tuple in the process.

@< Replace component at |index| in row |loc|... @>=
{ unsigned int i=(index->eval(),get<int_value>()->val);
  std::vector<shared_value>& a=force<row_value>(loc)->val;
  size_t n=a.size();
  if (i>=n)
    throw runtime_error(range_mess(i,a.size(),this,"component assignment"));
  auto& ai = a[reversed ? n-1-i : i];
  ai = pop_value(); // assign non-expanded value
  push_expanded(lev,ai); // return value may need expansion, or be omitted
}

@ For |vec_value| entry assignments the type of the aggregate object is
vector, and the value assigned always an integer. The latter certainly needs
no expansion, so we either leave it on the stack, or remove it if the value of
the component assignment expression is not used.

@< Replace entry at |index| in vector |loc|... @>=
{ unsigned int i=(index->eval(),get<int_value>()->val);
  std::vector<int>& v=force<vector_value>(loc)->val;
  size_t n=v.size();
  if (i>=n)
    throw runtime_error(range_mess(i,v.size(),this,"component assignment"));
  v[reversed ? n-1-i : i]= force<int_value>(execution_stack.back().get())->val;
    // assign |int| from un-popped top
  if (lev==no_value)
    execution_stack.pop_back(); // pop it anyway if result not needed
}

@ For matrix entry assignments at |index| must be split into a pair of
indices, and there are two bound checks.

@< Replace entry at |index| in matrix |loc|... @>=
{ index->multi_eval();
  unsigned int j=get<int_value>()->val;
  unsigned int i=get<int_value>()->val;
@/
  int_Matrix& m=force<matrix_value>(loc)->val;
  size_t k=m.numRows(),l=m.numColumns();
  if (i>=k)
    throw runtime_error(range_mess(i,m.numRows(),this,"matrix entry assignment"));
  if (j>=l)
    throw runtime_error(
      range_mess(j,m.numColumns(),this,"matrix entry assignment"));
  m(reversed ? k-1-i : i,reversed ? l-1-j : j)=
    force<int_value>(execution_stack.back().get())->val;
    // assign |int| from un-popped top
  if (lev==no_value)
    execution_stack.pop_back(); // pop it anyway if result not needed
}

@ A matrix column assignment is like that of a vector entry, but with a test
for matching column length.

@< Replace column at |index| in matrix |loc|... @>=
{ unsigned int j=(index->eval(),get<int_value>()->val);
  int_Matrix& m=force<matrix_value>(loc)->val;
@/const int_Vector& v=force<vector_value>(execution_stack.back().get())->val;
    // don't pop
  size_t l=m.numColumns();
  if (j>=l)
    throw runtime_error(
      range_mess(j,m.numColumns(),this,"matrix column assignment"));
  if (v.size()!=m.numRows())
    throw runtime_error
      (std::string("Cannot replace column of size ")+str(m.numRows())+
       " by one of size "+str(v.size()));
  m.set_column(reversed ? l-j-1 : j,v);
    // copy value of |int_Vector| into the matrix
  if (lev==no_value)
    execution_stack.pop_back(); // pop the vector if result not needed
}

@ For local assignments we also need to access the location where the
identifier is stored, which as before is done by storing coordinates of the
identifier in the execution context.

@< Type definitions @>=
template <bool reversed>
class local_component_assignment : public component_assignment<reversed>
{ typedef component_assignment<reversed> base;
@)
  subscr_base::sub_type kind;
  size_t depth, offset;
public:
  local_component_assignment @|
   (id_type arr, expression_ptr&& i,size_t d, size_t o,
    expression_ptr&& r, subscr_base::sub_type k);
  virtual void evaluate(expression_base::level l) const;
};

@ The constructor for |local_component_assignment| is straightforward, in
spite of the number of arguments.

@< Function def... @>=
template <bool reversed>
local_component_assignment<reversed>::local_component_assignment
 (id_type arr, expression_ptr&& i,size_t d, size_t o, expression_ptr&& r,
  subscr_base::sub_type k)
: base(arr,std::move(i),std::move(r)), kind(k), depth(d), offset(o) @+{}

@ The |evaluate| method locates the |shared_value| pointer of the aggregate,
calls |assign| to do the work.

@< Function def... @>=
template <bool reversed>
void local_component_assignment<reversed>::evaluate(expression_base::level l)
  const
{@; base::assign (l,frame::current->elem(depth,offset),kind); }

@ Type-checking and converting component assignment statements follows the
same lines as that of ordinary assignment statements, but must also
distinguish different aggregate types. Most of the code is straightforward,
but there is a subtle point that in case of a component assignment to a
variable of previously undetermined row type, the component type must be
recorded with the variable. This is achieved by the call to the |specialise|
for the pointer found at |aggr_t->component_type|, which is part of the type
for the variable in the local or global table. Here for one time we abuse of
the fact that, although |aggr_t| is a pointer-to-constant, we are still
allowed to call a non-|const| method for the |type_expr| that
|aggr_t->component_type| points to; otherwise we would have to call a
|specialise| method for the local or global variable table that holds the
identifier |aggr|, as was done in the case of ordinary assignments.

@:comp_ass_type_check@>

@< Cases for type-checking and converting... @>=
case comp_ass_stat:
{ id_type aggr=e.comp_assign_variant->aggr;
  const expr& index=e.comp_assign_variant->index;
  const expr& rhs=e.comp_assign_variant->rhs;
@/const_type_p aggr_t; size_t d,o; bool is_const;
  bool is_local = (aggr_t=layer::lookup(aggr,d,o,is_const))!=nullptr;
  if (not is_local and (aggr_t=global_id_table->type_of(aggr,is_const))==nullptr)
    report_undefined(aggr,e,"component assignment");
@.Undefined identifier@>
  if (is_const)
    report_constant_modified(aggr,e,"component assignment");
@.Name is constant @>
@)
  type_expr ind_t;
  expression_ptr i = convert_expr(index,ind_t);
@/type_expr comp_t;
  subscr_base::sub_type kind=subscr_base::index_kind(*aggr_t,ind_t,comp_t);
  if (not subscr_base::assignable(kind))
  { std::ostringstream o;
    o << "Cannot subscript value of type " << *aggr_t @|
      << " with index of type " << ind_t << " in assignment";
    throw expr_error(e,o.str());
  }
  expression_ptr r = convert_expr(rhs,comp_t);
  if (aggr_t->kind==row_type)
    aggr_t->component_type->specialise(comp_t); // record type
  if (comp_t==void_type and not is_empty(rhs))
    r.reset(new voiding(std::move(r)));
  expression_ptr p;
  if (is_local)
    if (e.comp_assign_variant->reversed)
      p.reset(new local_component_assignment<true>
        (aggr,std::move(i),d,o,std::move(r),kind));
    else
      p.reset(new local_component_assignment<false>
        (aggr,std::move(i),d,o,std::move(r),kind));
  else
    if (e.comp_assign_variant->reversed)
      p.reset(new global_component_assignment<true>
        (aggr,std::move(i),std::move(r),kind));
    else
      p.reset(new global_component_assignment<false>
        (aggr,std::move(i),std::move(r),kind));
  return conform_types(comp_t,type,std::move(p),e);
}

@* Some special wrapper functions.
%
In this chapter we define some wrapper functions that are not accessed through
the overload table; they must be directly visible to the type-checking code
that inserts them, which is why they are defined as local functions to the
current \.{axis.w} module.

@< Static variable definitions that refer to local functions @>=
static shared_builtin sizeof_row_builtin =
    std::make_shared<const builtin_value>(sizeof_wrapper,"#@@[T]");
static shared_builtin sizeof_vector_builtin =
    std::make_shared<const builtin_value>(sizeof_vector_wrapper,"#@@vec");
static shared_builtin sizeof_ratvec_builtin =
    std::make_shared<const builtin_value>(sizeof_ratvec_wrapper,"#@@ratvec");
static shared_builtin sizeof_string_builtin =
    std::make_shared<const builtin_value>(sizeof_string_wrapper,"#@@string");
static shared_builtin matrix_columns_builtin =
    std::make_shared<const builtin_value>(matrix_ncols_wrapper,"ncols@@mat");
static shared_builtin sizeof_parampol_builtin =
    std::make_shared<const builtin_value>(virtual_module_size_wrapper,
    "#@@ParamPol");
static shared_builtin print_builtin =
  std::make_shared<const builtin_value>(print_wrapper,"print@@T");
static shared_builtin to_string_builtin =
  std::make_shared<const builtin_value>(to_string_wrapper,"to_string@@T");
static shared_builtin prints_builtin =
  std::make_shared<const builtin_value>(prints_wrapper,"prints@@T");
static shared_builtin error_builtin =
  std::make_shared<const builtin_value>(error_wrapper,"error@@T");
static shared_builtin prefix_elt_builtin =
  std::make_shared<const builtin_value>(prefix_element_wrapper,"#@@(T,[T])");
static shared_builtin suffix_elt_builtin =
  std::make_shared<const builtin_value>(suffix_element_wrapper,"#@@([T],T)");
static shared_builtin join_rows_builtin =
  std::make_shared<const builtin_value>(join_rows_wrapper,"##@@([T],[T])");
static shared_builtin join_rows_row_builtin =
  std::make_shared<const builtin_value>(join_rows_row_wrapper,"##@@([[T]])");
static shared_builtin boolean_negate_builtin =
  std::make_shared<const builtin_value>(bool_not_wrapper,"not@@bool");

@ The function |print| outputs any value in the format used by the interpreter
itself. This function has an argument of unknown type; we just pass the popped
value to the |operator<<|. The function returns its argument unchanged as
result, which facilitates inserting |print| statements for debugging purposes.

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
  if (l!=expression_base::single_value) // in |single_value| case we are done
    push_expanded(l,pop_value()); // otherwise remove and possibly expand value
}

@ Sometimes the user may want to use a stripped version of the |print| output:
no quotes in case of a string value, or no parentheses or commas in case of a
tuple value (so that a single statement can chain several texts on the same
line). The |prints_wrapper| does this down to the level of omitting quotes in
individual argument strings, using dynamic casts to determine the case that
applies. The function |error| does the same, but collects the output into a
string which it then throws as |runtime_error|.

@< Local function definitions @>=
std::ostream& to_string_aux(std::ostream& o, expression_base::level l)
{ shared_value v=pop_value();
@)
  const string_value* s=dynamic_cast<const string_value*>(v.get());
  if (s!=nullptr)
    o << s->val; // single string without quotes
  else
  { const tuple_value* t=dynamic_cast<const tuple_value*>(v.get());
    if (t!=nullptr)
    { for (auto it=t->val.begin(); it!=t->val.end(); ++it)
      { s=dynamic_cast<const string_value*>(it->get());
        if (s!=nullptr)
	  o << s->val; // string components without quotes
        else
           o << *it->get(); // treat non-string tuple components as |print|
      }
    }
    else
      o << *v; // output like |print| unless string or tuple
  }
  return o;
}

void to_string_wrapper(expression_base::level l)
{ std::ostringstream o;
  to_string_aux(o,l);
  if (l!=expression_base::no_value)
    push_value(std::make_shared<string_value>(o.str()));
}

void prints_wrapper(expression_base::level l)
{ to_string_aux(*output_stream,l) << std::endl;
  if (l==expression_base::single_value)
    wrap_tuple<0>(); // don't forget to return a value if asked for
}

void error_wrapper(expression_base::level l)
@/{@; std::ostringstream o;
  to_string_aux(o,l);
  throw runtime_error(o.str());
}


@ The generic size-of wrapper is used to find the length of any ``row-of''
value. Finding sizes of other objects like vectors, matrices, polynomials,
will require more specialised unary overloads of the `\#' operator.

@< Local function definitions @>=
void sizeof_wrapper(expression_base::level l)
{ size_t s=get<row_value>()->val.size();
  if (l!=expression_base::no_value)
    push_value(std::make_shared<int_value>(s));
}


@ Here are functions for adding individual elements to a row value, and for
joining two such values.

@:hash wrappers@>

@< Local function definitions @>=
void suffix_element_wrapper(expression_base::level l)
{ shared_value e=pop_value();
  own_row r=get_own<row_value>();
  if (l!=expression_base::no_value)
  {@; r->val.push_back(std::move(e));
    push_value(std::move(r));
  }
}
@)
void prefix_element_wrapper(expression_base::level l)
{ own_row r=get_own<row_value>();
  shared_value e=pop_value();
  if (l!=expression_base::no_value)
  {@; r->val.insert(r->val.begin(),e);
    push_value(std::move(r));
  }
}
@)
void join_rows_wrapper(expression_base::level l)
{ shared_row second=get<row_value>();
  shared_row first=get<row_value>();
  if (l==expression_base::no_value)
    return;
  const auto& x=first->val;
  const auto& y=second->val;
  own_row result = std::make_shared<row_value>(x.size()+y.size());
@/std::copy(y.begin(),y.end(),
      @+ std::copy(x.begin(),x.end(),result->val.begin()) @+ );
@/push_value(std::move(result));
}

@)
void join_rows_row_wrapper(expression_base::level l)
{ shared_row arg=get<row_value>();
  if (l==expression_base::no_value)
    return;
  const std::vector<shared_value>& x=arg->val;
  std::vector< const std::vector<shared_value>*> p; p.reserve(x.size());
  for (auto it=x.cbegin(); it!=x.cend(); ++it)
    p.push_back(&force<row_value>(it->get())->val);
  size_t s=0;
  for (auto it=p.cbegin(); it!=p.cend(); ++it)
    s+=(*it)->size();
  own_row result = std::make_shared<row_value>(s);
  auto dst=result->val.begin();
  for (auto it=p.cbegin(); it!=p.cend(); ++it)
    dst=std::copy((*it)->cbegin(),(*it)->cend(),dst);
  assert(dst==result->val.end());
@/push_value(std::move(result));
}

@ Finally we define the Boolean negation wrapper function.
@< Local function definitions @>=
void bool_not_wrapper(expression_base::level l)
{ bool b=get<bool_value>()->val;
  if (l!=expression_base::no_value)
    push_value(whether(not b));
}
@)


@* Index.

% Local IspellDict: british
