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

#ifndef AXIS_H
#define AXIS_H

#include "Atlas.h" // must be very first \.{atlas} include

#include "axis-types.h"

@< Includes needed in the header file @>@;
namespace @;atlas { namespace interpreter {
@< Type definitions @>@;
@< Declarations of global variables @>@;
@< Declarations of exported functions @>@;
} }
#endif

@ The implementation unit follows a somewhat similar pattern.

@h "axis.h"
@h <cstdlib>
@c
namespace @;atlas { namespace interpreter {
@< Global variable definitions @>@;
namespace {
@< Local class definitions @>@;
@< Local variable definitions @>@;
@< Local function definitions @>@;
@< Static variable definitions that refer to local functions @>@;
}
@< Function definitions @>@;
} }

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

@ When a user interrupts the computation, we wish to return to the main
interpreter loop. To that end we shall at certain points in the program check
the status of a flag that the signal handler sets, and if it is raised call an
error. We also want to provide a mechanism to have functions time-out when their
computation is taking too long. These things require some system header files to
be included before our header file is read.

@< Includes needed... @>=
#include <csignal>
#include <chrono>

@~The flag must be a global variable, and we choose the traditional type
for it (even though our program currently does not use threads).

@< Declarations of global variables @>=

extern volatile std::sig_atomic_t interrupt_flag;

@~We shall define the variable right away, although our compilation unit does
not use it. It is cleared initially, and will be cleared every time an
interrupt is caught.

@<Global variable definitions @>=

volatile std::sig_atomic_t interrupt_flag=0;

@~When raised we shall call the following simple error value, with which no
data is associated.

@< Type definitions @>=

struct user_interrupt : public error_base {
  user_interrupt() : error_base("\nUser interrupt") @+{}
};
struct time_out : public error_base {
  time_out() : error_base("\nTimed out") @+{}
};

@ At several points in the interpreter we shall check whether a user interrupt
has raised |interrupt_flag|. Since we want to also be able to interrupt certain
built-in functions, we need to export this function.

@< Declarations of exported functions @>=
void set_timer(long long int milliseconds);
void clear_timer();
void check_interrupt();

@ We want to potentially use the user interrupt check also to test for time-out
of some functions that were called with a restriction on the time they are
allowed to spend on their computation.

@< Declarations of global variables @>=

static constexpr auto eternity = std::chrono::steady_clock::time_point::max();
static auto timer_out = eternity;

@ We provide simple functions to set and reset a time-out interval.

@< Function definitions @>=

void set_timer(long long int period)
{
  using namespace std::literals::chrono_literals;
  timer_out = std::chrono::steady_clock::now() + period * 1ms;
}
void clear_timer()
{@; timer_out = eternity;
}

@ When a user interrupt is handled, we clear the flag and throw an exception.
For time-out purposes the value |eternity| will supposedly never be surpassed,
so we could test |std::chrono::steady_clock::now()>timer_out| directly
regardless of what value was set in |timer_out|. Nonetheless we still test that
value against |eternity| first, since this avoids making a system call to fetch
the current time each time in the common cases where no time-out period has been
set in the first place.

@< Function definitions @>=
void check_interrupt()
{@;
  if (interrupt_flag!=0)
  {@; interrupt_flag=0;
    throw user_interrupt();
  }
  if (timer_out!=eternity and std::chrono::steady_clock::now()>timer_out)
  {@; clear_timer();
    throw time_out();
  }
}


@* Outline of the analysis and evaluation processes.
%
This module is concerned with the processing of the abstract syntax tree as
produced by the parser, ultimately producing actions and computed values. This
processing consist of two separate stages: type analysis, which also
transforms the syntax tree into a more directly executable form and therefore
might be called compilation, and execution of those transformed expressions.

The expression returned by the parser, of type |expr|, and the conversion to the
executable format |expression| (a type defined in \.{axis-types.w} as a pointer
to the base class |expression_base|, from which many more specialised classes
will be derived) is performed by the function |convert_expr|. This is a large
and highly recursive function, and a large part of the current compilation unit
is dedicated to its definition. The execution of the converted value is
performed by calling the (purely) virtual method |evaluate|, so that the code
describing the actual execution of expressions is distributed among the many
definitions of that method in derived classes, and this definition is only
implicitly (mutually) recursive through calls to the |expression_base::evaluate|
method.

@ During type checking, it may happen for certain subexpressions that a definite
type is required for them, while for others nothing is known beforehand about
their type (for instance this is the case for the complete expression entered by
the user for evaluation). The difference is important, as in the former case it
is possible to insert implicit conversions to make the types match, for instance
between a list of integers (built using the facilities of the interpreter) and a
vector value (one that can be directly used by the Atlas library); this is in
fact the only way the user can construct such vector values. However, both cases
(with known or unknown result type), and some intermediate cases (where the
result type is partially known) are handled by a single function |convert_expr|.
There are even cases where the expected type is any instance of a given
polymorphic type, which is why |convert_expr| takes a |type| as second argument
(where originally it took a |type_expr|). The first argument is an a |const
expr& e@;| referring to a value produced by the parser, which represents the
expression to be analysed.

If there are no restrictions on the final type, then |tp| will be undetermined
initially; if a precise type is required then it is set to that. It could also
be that there are partial requirements; then |tp| is a partially determined type
expression reflecting those. Then during the call to |convert_expr|, the value
of |tp| can be changed by calling |specialise| on it or on its subexpressions.
As polymorphic types were introduced, we discovered (fairly late) that we cannot
limit ourselves to partially defined context types determined by patterns with
undetermined parts, but we need to allow the requirement of being an instance of
some polymorphic type: while overloaded function calls must evaluate their
argument(s) in an undetermined type context, even if some function instances
are polymorphic, non overloaded calls do impose a type on their argument, which
may be a polymorphic type. For this reason we pass (by reference) a |type|
argument rather than a |type_expr|, and expect its type assignment to pick up
information during the type checking process.

The argument |tp| can be changed to a more specific one in various ways: by
calling |unify|, or calling |specialise| for (a subexpressions of) the contained
|type_expr|, or possibly by assigning a new more specialised value obtained in
another way. We also occasionally call the |expand| method to replace a tabled
sub-type by an equivalent one that makes the top level structure apparent. One
should not modify |tp| in other ways, and in particular not |std::move| from it.
Also |tp| should not refer to a table entry like the type of a variable, or to
any other permanent value, to avoid changes to the |tp| argument causing
permanent changes to them. (But after returning from |convert_expr|, |tp| can be
used in a permanent manner.)

@< Declarations of exported functions @>=
expression_ptr convert_expr (const expr& e, type& tp);

@ In some cases |tp| used to remain partly undefined, like for an empty list
display with an unknown type, which got specialised only to~`\.{[*]}'. Now that
|tp| is a |type|, that result will be given as the polymorphic type `\.{[A]}'
(to which that old result could have been converted by calling |type::wrap|). If
|tp| remains completely undefined (as will happen for an expression that selects
a value from an empty list, or that calls |error|) this will produce the
uninhabited type `\.A'; an expression having that type means that its evaluation
cannot possibly complete without error. (If one would try to deduce the return
type of a recursive function from its body, this situation would also happen if
such a function were to have no terminating case; in reality, the syntax insists
that recursive functions declare their return type explicitly, avoiding this
situation.)

@*1 Layers of lexical context.
%
In the function |convert_expr| we shall need a type for storing bindings between
identifiers and types, and this will be the |layer| class. It stores a
|layer::vec| of identifier bindings, while it automatically pushes (a pointer
to) itself onto a stack |layer::lexical_context| of such vectors. Instances of
|layer| shall always be stored in ordinary (stack-allocated) variables; since
their unique constructor pushes the new object, and their destructor pops it,
|layer::lexical_context| is guaranteed to hold at all times a list of pointers
to all current |layer| objects, from newest to oldest. Besides being cute, this
approach is exception safe: as long as all |layer| objects are held in automatic
variables their destructors will be called as the \Cpp\ stack unwinds, popping
them from the context. We do commit to not using threads during type analysis,
as this would break the strict stack discipline.

Ordinary methods are used to prepare the |layer| object inside the block where
it is constructed, while generally accessible lookup methods are implemented as
|static| methods, which use |lexical_context| to access the stack of layers.

We also use this structure to maintain additional attributes telling
whether we are inside a function and how deeply nested inside loops we are.
This is to facilitate checking the legality of \&{break} and \&{return}
expressions. Correct use of these could have been left to the parser at the
price of increasing the size of the grammar. These attributes may change only
when pushing/popping a layer. In most cases such a layer was needed anyway,
but in some cases (|while| loops) a layer is added specifically to mark this
change.

@< Type def... @>=
class layer
{ static simple_list<layer*> lexical_context;
public:
  struct id_data
  { id_type id; @+ type tp;
    id_data(id_type id, type&& tp) : id(id), tp(std::move(tp)) @+{}
  };
  using vec = std::vector<id_data>;
private:
  vec variable;
  BitMap constness;
  const unsigned loop_depth; // number of nested loops we are in
  type* const return_type;
    // address of return type of current function, if any (non owned)
public:
  layer(const layer&) = delete; // no ordinary copy constructor
  layer& operator= (const layer&) = delete; // nor assignment operator
  explicit layer(size_t n); // non-function non-loop layer
  layer(size_t n,type* return_type);
    // function or (with |return_type==nullptr|) loop layer
  ~layer () @+{@; lexical_context.pop_front(); }
@)
  void add(id_type id,type&& t, unsigned char flags);
  bool empty() const @+{@; return variable.empty(); }
  id_data& operator[] (size_t i) @+{@; return variable[i]; }
  vec::iterator begin() @+{@; return variable.begin(); }
  vec::iterator end() @+{@; return variable.end(); }
  vec::const_iterator cbegin() const @+{@; return variable.begin(); }
  vec::const_iterator cend() const @+{@; return variable.end(); }
  bool is_const (vec::const_iterator it) const
  @+{@; return constness.isMember(it-cbegin()); }
@)
  static const type* lookup
    (id_type id, size_t& depth, size_t& offset, bool& is_const);
  static const type* lookup (id_type id, size_t& depth, size_t& offset);
  static bool may_break(unsigned depth);
    // whether \&{break}~|depth| is legal here
  static bool may_return(); // whether \&{return} is legal here
  static type& current_return_type() @+
  {@; return *lexical_context.front()->return_type; }
};

@ Here are the constructors, which are used on three kinds of occasions: the
first one for \&{let} expressions, and the second one for loops (in which case
|return_type==nullptr|) and for user-defined functions (in which case
|return_type!=nullptr|). Effectively we have a loop counter starting at~$0$,
incremented for every loop and reset to~$0$ for every function body entered, and
a pointer to a return type that set for every function body, and otherwise
remains constant in newer layers.

@< Function def... @>=
layer::layer(size_t n) // non-function non-loop layer
: variable(), constness(n)
, loop_depth(lexical_context.empty() ? 0 : lexical_context.front()->loop_depth)
, return_type(lexical_context.empty() ? nullptr
             : lexical_context.front()->return_type)
{@; variable.reserve(n); lexical_context.push_front(this); }
layer::layer(size_t n,type* return_type) // function or loop layer
: variable(), constness(n)
,@/ loop_depth(return_type!=nullptr ? 0
            : lexical_context.empty() ? 1
            : lexical_context.front()->loop_depth+1)
,@/ return_type(return_type!=nullptr ? return_type @|
             : lexical_context.empty() ? nullptr
             : lexical_context.front()->return_type)
{@; variable.reserve(n); lexical_context.push_front(this); }

@ A statement of $n$ successive \&{break} keywords breaks out of $n$ levels of
loops at once, and will call |layer::may_break| with |depth==n-1| to see if it
is used legally. This means that a current |layer| must be present with
|loop_depth>depth|. For \&{return} to be legal it suffices to be anywhere inside
a function body, and the above constructors ensure that |return_type!=nullptr|
in the current layer gives this condition.

@< Function def... @>=
bool layer::may_break(unsigned depth)
{@; return not lexical_context.empty()
      and lexical_context.front()->loop_depth > depth; }
bool layer::may_return()
  {@; return not lexical_context.empty()
      and lexical_context.front()->return_type!=nullptr; }


@ The method |add| adds a pair to the vector of bindings; the type is moved
into the |layer| object.

@h <string>

@< Function def... @>=
void layer::add(id_type id,type&& tp, unsigned char flags)
{ @< Check that |id| is not already bound in our |layer| @>
  @< Check that we are not binding an operator to a non-function value @>
   constness.set_to(variable.size(),(flags&0x4)!=0);
  variable.emplace_back( id, std::move(tp) );
}

@ This is a good place to check for the presence of identical identifiers.

@< Check that |id| is not already bound in our |layer| @>=
{
  for (auto it=variable.begin(); it!=variable.end(); ++it)
  // traverse |variable| vector
    if (it->id==id)
      throw program_error @/
       (std::string("Multiple binding of '")
                    +main_hash_table->name_of(id)
                    +"' in same scope");
}

@ Just like for global definitions, we forbid locally binding a operator symbol
to an expression of non-function type.

@< Check that we are not binding an operator to a non-function value @>=
{
  if ((flags&0x8)!=0 and tp.top_kind()!=function_type)
    { std::ostringstream o;
      o << "Cannot bind operator '" << main_hash_table->name_of(id) @|
        << "' to an expression of non-function type " << tp;
      throw program_error(o.str());
    }
}

@ We need to define the |static| variable |layer::lexical_context| outside the
class definition; the default |sl_list| constructor ensures that it is initially
empty. The actual |layer|s are all local to the recursive function
|convert_expr|, and when the root call of that function terminates, the list
will be empty again.

@< Global var... @>=
simple_list<layer*> layer::lexical_context;

@ The method |layer::lookup| runs through the linked list of layers, and if a
match for the identifier |id| was found it returns a pointer to its type, while
also assigning its static binding coordinates to output arguments |depth| and
|offset|. If no match is found, a null pointer is returned and the output
parameters are unchanged. We allow having a |layer| with no variables, for which
no run time stack frame will correspond at all (as such a frame would incur a
runtime cost for no good at all). The method |lookup| skips such layers without
increasing the |depth| it reports.

The loop variable |range| below is a const-iterator over a |simple_list|, where
|*range| is a list item which is a pointer to |layer|. We therefore access
|layer| members using a double dereference of |range| (the second one hidden in
the arrow symbol).

@< Function def... @>=
const type* layer::lookup (id_type id, size_t& depth, size_t& offset)
{@; bool dummy; return lookup(id,depth,offset,dummy); }
const type* layer::lookup
  (id_type id, size_t& depth, size_t& offset, bool& is_const)
{ size_t i=0;
  for (auto range=lexical_context.cbegin();
       not lexical_context.at_end(range);
       ++range)
    if (not (*range)->variable.empty())
    { for (auto it=(*range)->cbegin(); it!=(*range)->cend(); ++it)
        if (it->id==id) // then found; now set output values
          { depth=i;
            offset=it-(*range)->begin();
            is_const = (*range)->is_const(it);
            return &it->tp;
          }
      ++i; // increment depth for non-empty layers only
    }
  return nullptr;
}

@ We come to our highly recursive main type analysis function |convert_expr|. It
returns a owning pointer |expression_ptr| to the result of converting the |expr|
to an |expression|.

Altogether this is a quite extensive function, with as many cases in the
switch as there are variants of |expr|, and for many of those branches a
considerable amount of work to be done. It is therefore convenient to postpone
these cases and treat them one syntactic construction at the time.

@< Function definitions @>=
expression_ptr convert_expr(const expr& e, type& tp)
{ const auto fc = tp.floor(); // short for ``fixed count''
  switch(e.kind)
  {
   @\@< Cases for type-checking and converting expression~|e| against
   |tp|, all of which either |return| or |throw| a |type_error| @>
   case no_expr: {} // should not occur
  }
  assert(false);
  return expression_ptr(nullptr); // keep compiler happy
}

@ Sometimes we need to convert an expression in the context of a fixed fully
determined type, which will therefore not change in the call to |convert_expr|.
Nonetheless it requires a non-|const| argument that cannot be provided by an
rvalue as would result from calling the |copy| method in the argument
expression. So to cater for that, we define a local function
|convert_expr_strongly| that will provide the required variable, and |assert|
that it is being used correctly.

@< Local function definitions @>=
expression_ptr convert_expr_strongly
  (const expr& e, unsigned int fc, const type_expr& te)
{
  type tp = type::wrap(te,fc);
  assert(te.raw_kind()!=undetermined_type and tp.degree()==0);
  return convert_expr(e,tp);
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
be a tuple. Since then we also added constant folding to the interpreter, so
that denotations for values of any type can also come into existence in other
ways.

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
construct and |new|-allocate an \.{axis} value (for instance |int_value|)
from it, making the pointer shared using |std::make_shared|, pass that
pointer to the |denotation| constructor, and convert the resulting pointer to
a unique pointer.

The code below takes into account the possibility that a denotation is
converted immediately to some other type, for instance integer denotations can
be used where a rational number is expected. The function |conform_types|
(defined in \.{axis-types.w}) will test whether the denotation provides or can
be converted to the required type, and may modify its final argument in the
latter case.

@< Cases for type-checking and converting... @>=
case integer_denotation:
  { expression_ptr d(new denotation @|
      (std::make_shared<int_value>(
         big_int(e.str_denotation_variant->c_str(),10))));
    return conform_types(int_type,tp,std::move(d),e);
  }
case string_denotation:
  { expression_ptr d(new denotation @|
      (std::make_shared<string_value>(*e.str_denotation_variant)));
    return conform_types(str_type,tp,std::move(d),e);
  }
case boolean_denotation:
  { expression_ptr d(new denotation @|
        (whether(e.bool_denotation_variant)));
    return conform_types(bool_type,tp,std::move(d),e);
  }

@ We allow using the last value computed to be used in an expression, using
the symbol `\.\$'.

@< Declarations of global variables @>=
extern type last_type;
extern shared_value last_value;

@~We set the pointers to |nullptr| here, but the |main| function will give them
more appropriate starting values.

@< Global variable definitions @>=
type last_type = type::bottom(0);
shared_value last_value;

@ In some occasions, a previously computed value can be captured in an
expression (currently this applies to `\.\$', which captures the last computed
value, and of operator casts, which capture a function value from the overload
table). In those cases we shall use an expression type that is like
|denotation| so that evaluation will give back the captured value; however for
the purpose of printing (if this expression occurs inside a function body) it
is undesirable to embark on printing the whole captured value, so we derive
and override the |print| method.

@< Includes needed in the header file @>=

#include <string>

@~@<Type definitions @>=
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
    ,tp
    ,expression_ptr(new capture_expression(last_value,o.str()))
    ,e);
}

@*1 The suicidal expression.
%
In some cases, notably for facilitating recursive definitions, it is useful to
have a placeholder expression \&{die} that will assume any type the context
requires (we even allow an undetermined type, although in its application to
recursion the type will always be defined). The evaluation of this placeholder
is not intended, and trying to evaluate it throws a runtime error; for this
reason the expression is written \&{die}. One could require that a call to
function~$f$ returns |true| by writing ``$f(...)$~\&{or die}''.

These expressions must be representable at run time, so we define an empty
|shell| for them.

@< Type definitions @>=
struct shell : public expression_base
{
virtual void evaluate (level l) const;
virtual void print(std::ostream& out) const @+{@; out << " die "; }
};

@ As said above, attempting to evaluate a |shell| is lethal.

@< Function definitions @>=
void shell::evaluate (level l) const
{@; throw runtime_error("I die"); } // our |shell| explodes


@ The main point of \&{die} is not trying to evaluate, but allowing it to pass
type checking successfully. It does so trivially. Leaving |tp| an undetermined
type will cause, due to how |type::wrap| functions, the type deduced for this
expression to be the impossible ``bottom'' (or ``any type'') type.

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
{@; unsigned d=depth;
  do out << " break";
  while (d-->0);
  out << ' ';
}
};

@ The \&{break} is realised by the \Cpp\ exception mechanism, which allows its
evaluation to terminate that of intermediate clauses, without burdening the
latter with the need to test for this possibility. An error type |loop_break| is
defined for this purpose. During type checking we only allow \&{break} in places
where it is certain to be caught: inside the required number of nested loops,
with no function definitions (lambda expressions) intervening between
the \&{break} and the outermost of those loops.

@< Function definitions @>=
void breaker::evaluate (level l) const
@+{@; throw loop_break(depth); }


@ The only check we do for \&{break} is that the lexical context is suitable for
breaking out of the required number of loops, which number $n+1$ where $n$ is
the number stored in |e.break_variant|. The test for this condition is performed in
the method |layer::may_break|.


@< Cases for type-checking and converting... @>=
case break_expr:
{ if (layer::may_break(e.break_variant))
    return expression_ptr(new breaker(e.break_variant));
  if (e.break_variant==0)
    throw expr_error(e,"Using 'break' not in the reach of any loop");
  else
  { std::ostringstream o;
    o << "Using 'break";
    for (unsigned i=0 ; i<e.break_variant; ++i)
      o << " break";
    o << "' requires " << e.break_variant+1 << " nested levels of loops";
    throw expr_error(e,o.str());
  }
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
|layer::current_return_type()| to the return type of the current function. It
provides both the type context for the expression after~\&{return}, and in case
that was undetermined initially, a place where the type determined here by
|convert_expr| can be transmitted to later places where a return value for the
same function is produced (either by another |return| expression, or by one
directly producing the value of the function body).

When a function has multiple expressions that (in different circumstances) will
produce its value, the |return| expressions do not participate in any
``balancing'' (a notion to be described later) to determine the result type of
the function. Indeed, it would be difficult to do that (with our implementation
of balancing) since it may require converting the expressions a second time in a
more strict type context, but at the time this is discovered the
|layer::lexical_context| may have changed from the initial conversion, and would
be hard to reestablish. Instead it will be the first explicit type found (in the
order of type checking, which is generally from left to right) that determines
the result type. If the function body is a cast, that gives the type, otherwise
the first |return| expression gives the type, and if there are none, any further
expressions occurring in branches of conditionals or integer case expressions
are balanced to determine the result type. This rule results implicitly from the
fact that the function body and all |return| expressions share a same |type|
variable to record their type, which is the one |layer::current_return_type()|
refers to.

@< Cases for type-checking and converting... @>=
case return_expr:
{ if (layer::may_return())
  { type& rt = layer::current_return_type(); // shared return type
    return expression_ptr(new returner(convert_expr(*e.return_variant,rt)));
  }
  throw expr_error(e,"One can only use 'return' inside a function body");
}

@* Tuple displays.
%
Tuples are sequences of values of non-uniform type, usually short and with a
given sequence of types for their components; their most obvious and simple
use is for argument lists of functions. Without using these tuple types, one
might have a mechanism explicitly allowing functions to take multiple arguments,
as most programming languages do, but that would not allow for functions to
\emph{return} more than one value. With tuple types, we cater at the same time
for functions that need multiple arguments and for those that need to return
multiple results.

Tuple values can be produced by tuple displays, which list explicitly their
component expressions. After type-checking, they are given by a
|tuple_expression| object (the name |tuple_display| was already taken).

Sometimes we wish to indicate that the components of a tuple expression should
be evaluated right-to-left; we make this distinction possible through a Boolean
template parameter |r_to_l|.

@< Type definitions @>=
template<bool r_to_l> struct tuple_expression_tmpl : public expression_base
{ std::vector<expression_ptr> component;
@)
  explicit tuple_expression_tmpl(size_t n) : component(n) @+{}
   // always start out with null pointers
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

using tuple_expression = tuple_expression_tmpl<false>;

@ When we print a tuple display, we just print the component expressions,
enclosed in parentheses and separated by commas, to match their input syntax.


@< Function def... @>=
template<bool r_to_l>
  void tuple_expression_tmpl<r_to_l>::print(std::ostream& out) const
{ out << '(';
  for (auto it=component.begin(); it!=component.end(); ++it)
    out << (it==component.begin() ? "" : ",") << **it;
  out << (r_to_l ? "~)" : ")");
}

@ When converting a tuple expression, a difficulty to keep in mind is that there
may be components that turn out to have a polymorphic type, in which case we
must ensure that they are converted at an unchanged level of fixed type
variables (because any new type variables introduced inside the component
expression will be numbers by the lexical analyser at that level), but that in
the resulting type the polymorphic variable may then need renumbering to make
then disjoint form those of previous components. In a case where this might
happen, we do not want to pass a component of |tp| to the recursive call of
|convert_expr| for the component, as it will get specialised to the wrong type
and do not want to then later assign a renumbered component type (maybe it would
be acceptable since we are still assigning a specialised version of the original
type pattern, but we do not want to break a convention this is adhered to
everywhere else). So whenever this might happen we pass a |type_expr| disjoint
from |tp| to |convert_expr|, and specialise |tp| only after possibly creating a
renumbered polymorphic type.

One property that we do exploit here is that any determined components of |tp|
will be so to a monomorphic type. (We used to think that |tp| is either
completely undetermined or completely determined and monomorphic, after tweaking
the processing of \&{let} expressions to make it hold there, but we forgot about
multiple assignments: when the left hand side fails to specify a destination
for some part of the value, the right hand side has no type for that
part, while it has for other parts.) So we choose to provide a disjoint type
when calling |convert_expr| for a component whenever the corresponding
component of |tp| is unstable (as long as |gap==0| we could avoid the
renumbering of a polymorphic type but we still would need to update |gap|, so
we choose not to treat this case separately).

We also need to keep in mind the case where the required type can only be
obtained using |coerce| (currently the only possible coercion is the
Atlas-specific conversion from a pair of integers to a |Split_integer|). This
sets |n_tuple_expected==false| at the beginning, and uses a separate $n$-tuple
type |ntt| in place of |tp| to record the tuple type before coercion.

If any type error is reported here (and not from a nested call of
|convert_expr|), then it must be because |coerce| was attempted and failed. In
that case |ntt| now holds the a priori $n$-tuple type of the expression, and
|tp| is unmodified since entering this case. Therefore the error message given
by |type_error| describes faithfully the actual and expected type patterns.

@h <memory> // for |std::unique_ptr|
@< Cases for type-checking and converting... @>=
case tuple_display:
{
  auto n=length(e.sublist);
  type_expr ntt=unknown_tuple(n); // ensure an $n$-tuple type
  tp.unify_specialise(ntt);
    // maybe fill components of |ntt| according to |tp|
  bool is_constant = true; // whether all components are denotations
@)
  std::unique_ptr<tuple_expression> tup_exp(new tuple_expression(0));
  std::vector<expression_ptr>& comp = tup_exp->component;
  comp.reserve(n);
  sl_list<type> comp_types;
@)
  wtl_const_iterator tl_it (ntt.tuple());
  for (wel_const_iterator it(e.sublist); not it.at_end(); ++it,++tl_it)
  {
    type& comp_type = comp_types.push_back(type::wrap(*tl_it,fc));
    if (tl_it->is_unstable())
      @< Call |convert_expr| for |*it| with |comp_type| @>
    else
      comp.push_back(convert_expr(*it,comp_type));
    @< Handle |is_constant| maintenance and possible need for voiding @>
  }
@)
  expression_ptr result(std::move(tup_exp));
    // convert |tup_exp| to a more generic pointer
  if (is_constant)
     make_row_denotation<false>(result);
     // wrap tuple inside a denotation, see below
  auto tup_tp = type::wrap_tuple(std::move(comp_types));
  if (tp.unify_to(tup_tp) or coerce(tup_tp.bake(),tp.bake(),result,e.loc))
    return result;
  throw type_error(e,tup_tp.bake_off(),tp.bake());
}

@ This module was greatly simplified when allowing |tp| to be a polymorphic
type. We already prepared |comp_type| from the |type_expr| value |*tl_it| to
import any type requirement from the context, and it aliases the current element
of |comp_types|, which will later be wrapped into a type to be unified with
|tp|. So all the remains to be done here is call |convert_expr|, and push the
converted expression onto the vector~|comp|.

@< Call |convert_expr| for |*it| with |comp_type| @>=
{
  comp.push_back(convert_expr(*it,comp_type));
}

@ After processing a tuple component in whichever way, we perform two small
auxiliary tasks. We record in |is_constant| whether we remain on track for
constant-folding the tuple expression into a constant tuple value, and we test
tot see whether a |voiding| needs to be inserted. The latter applies only in the
rare case that a void type is found (possibly because to was required by |tp|)
but the component expression is not an empty tuple display; when it does, the
inserted |voiding| invalidates any |in_constant| that we might still be hoping
for.

@< Handle |is_constant| maintenance and possible need for voiding @>=
{
  is_constant = is_constant and
    dynamic_cast<const denotation*>(comp.back().get()) != nullptr;
  if (*tl_it==void_type and not is_empty(*it))
    is_constant=false,comp.back().reset(new voiding(std::move(comp.back())));
}

@ When all components of a tuple expression are denotations, we make the tuple
expression into a a denotation with a tuple value. Since a very similar
operation will be needed for row displays, we implement this by a function
template with a Boolean argument, which when |true| will produce a |denotation|
holding a |row_value|, rather than a |tuple_value| as it does in the call from
the code above.

@< Local function definitions @> =
template<bool is_row>
  void make_row_denotation(expression_ptr& p)
{
  const std::vector<expression_ptr>& comp = is_row @|
  ? static_cast<const list_expression*>(p.get())->component @|
  : static_cast<const tuple_expression*>(p.get())->component;
  tuple_value tup(comp.size());
  for (size_t i=0; i<comp.size(); ++i)
    // transfer denotation contents into |tup|
    tup.val[i] = static_cast<const denotation*>(comp[i].get())->denoted_value;
  shared_value val;
  if (is_row)
    val = std::make_shared<row_value>(std::move(tup));
  else
    val = std::make_shared<tuple_value>(std::move(tup));
  p.reset(new denotation(std::move(val)));
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
the |execution_stack| for some time, and the only one where a \&{break} might
occur at such a point (this should be quite rare, but it is possible). It is
therefore easier to clean the stack up in the code below, than to mark the
execution stack before entering any loop that might terminate with \&{break},
and resetting it to the marked position when the |loop_break| is caught (in
particular since there are many kinds of loops that would need consideration).
However, while we thus localise the clean-up in a single place, this code does
get executed very often (every time a built-in function is called, except if
it has exactly $1$ argument). We therefore take care (by not marking anything
at loop entry) that in the non-throwing case (the vast majority) no cycles are
wasted at all (assuming, as seems reasonable, that simply entering a
|try|-block does not involve any work).

@< Function def... @>=
template<>
  void tuple_expression::evaluate(level l) const
{ switch(l)
  {
  case level::no_value:
    for (auto it=component.begin(); it!=component.end(); ++it)
      (*it)->void_eval();
    break;
  case level::single_value:
    { auto result = std::make_shared<tuple_value>(component.size());
      auto dst_it = result->val.begin();
      for (auto it=component.cbegin(); it!=component.cend(); ++it,++dst_it)
      @/{@; (*it)->eval(); *dst_it=pop_value(); }
      push_value(result);
    } break;
  case level::multi_value:
    { auto it=component.begin();
      // this variable needs to survive |try| into the |catch| blocks
@)
      try
      {@; for (; it!=component.end(); ++it)
          (*it)->eval();
      }
@)
      catch (const loop_break&) // clean up execution stack locally
      { for (; it!=component.begin(); --it)
          execution_stack.pop_back();
        throw; // propagate the \&{break}
      }
      catch (const function_return&) // clean up execution stack locally
      { for (; it!=component.begin(); --it)
          execution_stack.pop_back();
        throw; // propagate the \&{return}
      }
    }
  } // |switch(l)|
}

@ And here is the right-to-left evaluation function. For the |no_value| case we
simply |void_eval| the components in reverse order; for the |single_value| case
the results of |eval| are placed as |component| value of a pre-allocated tuple
from back to front. In the |multi_value| case the values are needed separately
on the |execution_stack|, but in reverse order to what would be produced by the
successive |eval| calls. We therefore pop the values off after each evaluation,
storing them temporarily in a local |stack<shared_value>| variable, from which
they are then popped again afterwards to be placed on the |execution_stack| in
reverse order. Since the values are now held in a local variable, the
|execution_stack| remains unchanged in between these evaluations, and the is no
need for a |try|--|catch| construction to correctly handle possible exceptions
thrown during the evaluation of a component.

@< Function def... @>=
template<>
  void tuple_expression_tmpl<true>::evaluate(level l) const
{ switch(l)
  {
  case level::no_value:
    for (auto it=component.crbegin(); it!=component.crend(); ++it)
      (*it)->void_eval();
    break;
  case level::single_value:
    { auto result = std::make_shared<tuple_value>(component.size());
      auto dst_it = result->val.rbegin(); // fill |result| tuple from rear to front
      for (auto it=component.crbegin(); it!=component.crend(); ++it,++dst_it)
      @/{@; (*it)->eval(); *dst_it=pop_value(); }
      push_value(std::move(result));
    } break;
  case level::multi_value:
    { stack<shared_value> args;
      for (auto it=component.crbegin(); it!=component.crend(); ++it)
      @/{@; (*it)->eval(); args.push(pop_value()); }

@)    while(not args.empty()) // push components onto stack again in reverse order
        push_value(std::move(args.top())),args.pop();
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
{ if (l==level::no_value)
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
conditional expressions and possibly other cases: all component expressions need
to be of the same type, but we might not know which. We adopt the rule (borrowed
from Algol~68) that at least one of the components already gets the common type
when converted in the context of the original type (pattern), while the others
can be accepted in a context requiring that common type, possibly either by
inserting coercions or in case they have a polymorphic type by certain
substitutions to certain type variables. For instance a list display containing
expressions of \foreign{a priori} types \.{int} and \.{rat} will compile as a
list of \.{rat} expressions, or for a polymorphic example the list display
$[([\,],[3/4,5]),([2],[\,])]$ will compile as type \.{[([int],[rat])]}. (The
latter example is polymorphic due to the presence of empty lists. This somewhat
contrived example incidentally did not compile without a cast before the
introduction of polymorphic types into the language. It should be said that
since our algorithm chooses between either using coercions, for monomorphic
types, or using substitutions, for polymorphic types, it still cannot handle
cases where these could operate in parallel, as in
$[([3/4,5/7],[\.{true}]),([2],[\,])]$; this is because of the way the auxiliary
function |join_to| used operates.) This process is called type balancing (the
mental image is comparing the ``weights'' of the component expression types, to
see which one is the most accepting in determining the result type).

The general situation where balancing can be applied is that we have a list of
expressions that can all be type checked successively, producing a common type
and a list of converted expressions. This means that for nested conditional
expressions, we only balance two subexpressions each time; if one or both of
these subexpression are themselves balanced internally, that involves separate
calls of |balance|. Also, if for instance we have a row display some of whose
component expressions are also row displays, the latter are handled by separate
calls of |balance| (but there may be interaction with the outer call of
|balance| in case the inner balance fails, through the catching of a
|balance_error|).

@ Here is the set-up for balancing. We try to find in |common| a type that can
accept all component expressions. Branch types incomparable with the current
value of |common| are put aside in |conflicts|. At the end, |common| may have
become accepting enough to accommodate them (for instance |void| accepts every
possible type; we will not however make |common| be |void| unless at least one
expression \emph{requires} a void type). So at the end we prune |conflict|
before possibly reporting an error.

In case of success, the context type |target| that was used as initial goal for
the conversion of the balanced expressions is set to |common|. The change here
can only involve specialisation; although |target| initially cannot be
polymorphic (just undetermined), the specialisation can make it polymorphic if
|common| is. Then any branches that were originally found to have a different
type than |common|, but apparently were saved by pruning (so the final value of
|common| accepts them), are converted again in the context of the |target| type.
The second conversion can only differ from the first by insertion of coercions,
so we skip this phase altogether if |common| is polymorphic; in that case the
fact that pruning succeeded shows that |common| accepts all types of components,
and this suffices.

@< Local function definitions @>=

void balance
   ( type& target // component type required from context
   , raw_expr_list elist // list of expression to be balanced
   , const expr& e // containing expression, for error reporting
   , const char* description // what kind of components, for error message
   , std::vector<expression_ptr>& components // output, converted expressions
   )
{
  unsigned n = length(elist); components.reserve(n);
  std::vector<type> comp_type; comp_type.reserve(n);
  type common // will become the common type accepting all branch types
    = target.copy(); // every expression initially sees the |target| type
  containers::sl_list<type> conflicts;
    // except those branch types that are put aside here
  @< Convert each expression in |elist| in the context of a copy of |target|,
     pushing the results to |components|; maintain |common| as unified type,
     record in |conflicts| non conforming component types @>
  @< Prune from |conflicts| any types that |common| now accepts, and if any
     conflicts remain, |throw| a |balance_error| mentioning |common| and all
     |conflicts|; otherwise specialise |target| to |common| @>

  if (not common.is_polymorphic())
    @< Redo conversion with context type |target| for components that do not
       already have that type @>
}

@ We try to maintain |common| as unifying type between the branches, as computed
by repeated application of |join_for|. If this fails due to incomparable types
we move the non-conforming type to |conflicts|. Changing |common| to incorporate
the presence of |ctype| as a priori type of a branch is what the call
|join_to(common,ctype)| accomplishes; it decides between changing according to a
coercion, unifying a polymorphic type, doing nothing if |common| already accepts
|ctype|, or returning |false| if nothing can be done.

The list |conflicts| also collects types from any |balance_error| thrown
directly by one of the calls to |convert_expr|. This is because it might happen
that the conflicting types causing the |balance_error| can nonetheless all be
converted to a final |common| type that ends up getting chosen due to other
branches. So we |catch| such an error and just record the types involved; if
there ends up being a |balance_error| for the parent expression anyway, they
will be reported among the offending types. However, when catching a
|balance_error| that was produced in a lower level subexpression, we do
re-|throw|, as a more accepting type |common| cannot salvage in this case: the
thrown |balance_error| is in that case independent of our current balancing and
to be reported directly. For the case where we do move the types from the
|balance_error| to |conflicts|, the fact that they apparently have internal
incompatibilities means that they can never provide the |common| type here,
which justifies making no comparison with that type.

If a branch that produces a maybe salvageable |balance_error| is a list display,
the type that it will give for our current balancing has an additional
``row-of'' with respect to the types balanced in the subexpression. The
|wrap_row| call make the required adjustment. These wrapped types will not be
the |common| type, but at least they will compare correctly to it for the
purpose of pruning.

@< Convert each expression in |elist| in the context... @>=
for (wel_const_iterator it(elist); not it.at_end(); ++it)
{ try
  { components.push_back(expression_ptr());
       // push, whether or not |convert_expr| succeeds
    type ctype = target.copy();
    components.back() = convert_expr(*it,ctype);
    comp_type.push_back(ctype.copy());
    if (not join_to(common,std::move(ctype)))
      conflicts.push_back(std::move(ctype));
      // record type not convertible to |common|
  }
  catch (balance_error& err)
  { if (&err.offender!=&*it) // only incorporate top-level balancing errors
      throw; // any deeper error is propagated to be reported
    else if (err.offender.kind==list_display)
      // then wrap variants in row-of
      for (auto jt=err.variants.wbegin(); not err.variants.at_end(jt); ++jt)
        jt->wrap_row();
        // replace |*jt| by its ``row-of'' type
    conflicts.append(std::move(err.variants)); // then join to our |conflicts|
  }
}

@ Pruning is quite simple, and gives us an occasion to exercise the |erase|
method of |containers::sl_list|. In such loops one should not forget
to \emph{not advance} the iterator in case a node is erased in front of it.

Only if at least one conflicting type remains do we report an error; if so, the
type |common| is added as first type to the error object, unless it is unchanged
from the |type::bottom(fc)| value it was initialised to (which may happen if
every branch threw a |balance_error| that was caught), so that one has a
complete list of types that caused balancing to fail. Upon success, we
specialise |target| to a type accepting |common| (which might be polymorphic, or
even |type::bottom(fc)|, as happens in an empty list display which has $0$
branches); the method |type::unify_specialise| accomplishes this.


@< Prune from |conflicts| any types... @>=
{ for (auto it=conflicts.begin(); not conflicts.at_end(it); )
      // no increment here!
    if (join_to(common,std::move(*it))) // only actually moves if it succeeds
      conflicts.erase(it);
    else
      ++it;
  if (not conflicts.empty())
  { balance_error err(e,description);
    if (not common.is_bottom())
      err.variants.push_back(std::move(common));
    err.variants.append(std::move(conflicts));
    throw std::move(err);
  }
  if (not common.unify(target))
    throw type_error(e,common.bake_off(),target.bake());
}

@ When converting a component expression a second time in the hope it can adapt
to the |target| type that results from balancing, the result of a successful
conversion replaces the previous result in the |components| vector. There is no
guarantee the new conversion succeeds, but we won't catch any exceptions this
second time.

The same treatment is given to branches whose conversion threw a balancing error
when initially converted (meaning they had internal balancing that failed in
their initial type context), but for which all types involved in the balancing
were later pruned in the outer balancing. Such branches leave their |comp_type|
at the value |target| had for the attempted initial conversion, which must
differ from its final value if pruning removed the offending types. Therefore
the test that |comp_type[i].unwrap()!=target| holds for them, and
|components[i]|, which was left unset by the failed initial conversion, gets set
by the new conversion (if it is successful)).

@< Redo conversion with context type |target| for components that do not
   already have that type @>=
{ wel_const_iterator it(elist);
  type_expr te = target.bake();
  for (unsigned i=0; i<n; ++i,++it)
    if (not comp_type[i].is_polymorphic() and
        comp_type[i].bake_off()!=te)
      components[i] = convert_expr(*it,target);
      // redo conversion with unifying |target| type
}

@ With balancing implemented, converting a list display becomes fairly easy.
The simplest case is one where |tp| is a row type (or |undefined_type|, which
can be specialised to such). In that case we prepare an initially empty
|list_expression|, then call |balance| with the component type of |tp|,
which if successful will have converted to component types into our
|list_expression|, and it remains to |return| that object.

The case where a list display occurs in a void context is rare but valid.
(Since no list will be created, even though all component expressions will be
evaluated, the choice of writing a list display is rather curious.) For it we
perform balancing with an undetermined component type (as if the display were
in undetermined type context).

In the remaining case we call |row_coercion| (defined in \.{axis-types.w}) to
see if a coercion to |tp| from some row of |comp_type| exists. If this is
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
{ auto* const lp = new list_expression(0);
@/auto& comps = lp->component;
  comps.reserve(length(e.sublist));
  expression_ptr result(lp); // convert to |std::unique_ptr<expression_base>|
  auto is_const = @[ [] (const expression_ptr& p)
    {@; return dynamic_cast<const denotation*>(p.get()) != nullptr; } @];
@)
  static const char* const str = "components of list expression";
  if (tp.is_void()) // in void context leave undetermined target type
  { type comp_type = type::bottom(fc);
    balance(comp_type,e.sublist,e,str,comps);
      // no need to export back to |comp_type|
    return result; // and forget |comp_type| and result of |balance|
  }
@)
  type_expr pattern = row_of_type.copy();
  if (tp.unify_specialise(pattern))
  {
    type comp_type = type::wrap(pattern.component_type(),fc);
    balance(comp_type,e.sublist,e,str,comps);
    if (comp_type.is_void())
      @< Insert voiding coercions into members of |comps| that need it @>
    if (std::all_of(comps.begin(),comps.end(),is_const))
      make_row_denotation<true>(result);
    tp.unify(comp_type.wrap_row());
    return result;
  }
@)
  type_expr comp_type;
  const conversion_record* conv = row_coercion(tp.bake(),comp_type);
  if (conv!=nullptr)
  { type target = type::wrap(comp_type,fc);
    balance(target,e.sublist,e,str,comps);
      // no need to export |target| back to |comp_type|
    if (std::all_of(comps.begin(),comps.end(),is_const))
    @/{@;
      make_row_denotation<true>(result);
      do_conversion(result,*conv,e.loc);
    }
    else
      result.reset(new conversion(*conv,std::move(result)));
    return result;
  }
@)
  throw type_error(e,std::move(pattern),tp.bake());
  // |tp| incompatible with any list
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

@* Identifiers, global and local.
%
Identifiers are used to access values of all types, and also for designating
overloaded functions. In the latter usage a single identifier can be used to
access (depending on its immediate context) one of many values, and it
requires a certain amount of work to determine which one; this will be
discussed later, so for the moment we stick to simple applied identifiers that
identify the closest defining occurrence of that identifier in the current
lexical context. This identification can result in two outcomes: it may be
bound to a local or to a global name, which two cases are treated in fairly
different ways. In particular after type analysis the two cases are converted
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
struct identifier : public expression_base
{ id_type code;
@)
  explicit identifier(id_type id) : code(id) @+{}
  virtual ~identifier() = default;
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

We add a Boolean template parameter to this class, in order to have a variant
where the evaluation empties the variable itself, which the compiler can employ
for efficiency purposes.

@< Type definitions @>=
template<bool pilfer=false>
  class global_identifier : public identifier
{ const shared_share address;
public:
  explicit global_identifier(id_type id);
  virtual ~global_identifier() = default;
  virtual void print(std::ostream& out) const;
  virtual void evaluate(level l) const;
};

@ The constructor for |global_identifier| locates the value
associated to the identifier in the global identifier table. Because of the
reference to |global_id_table| we made it definition out-of-line.

@< Function definitions @>=
template<bool pilfer>
  global_identifier<pilfer>::global_identifier(id_type id)
: identifier(id), address(global_id_table->address_of(id))
@+{}

@)
template<>
  void global_identifier<false>::print(std::ostream& out) const
@+{@; out << name(); }
template<>
  void global_identifier<true>::print(std::ostream& out) const
@+{@; out << '$' << name(); }

@ Evaluating a global identifier returns the value currently stored in the
location |address|, possibly expanded if |l==multi_value|, or nothing at all
if |l==level::no_value|. However, since undefined global variables have been made
possible in the language, we have to watch out for a (shared) null pointer at
|*address|.

@< Function definitions @>=
template<bool pilfer>
  void global_identifier<pilfer>::evaluate(level l) const
{ if (address->get()==nullptr)
  { std::ostringstream o;
    o << "Taking value of uninitialized variable '" << name() << '\'';
    throw runtime_error(o.str());
  }
  if (pilfer)
    push_expanded(l,std::move(*address));
  else
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

Like for global identifiers we add a Boolean template parameter to this class,
which indicates whether the evaluation empties the variable itself, which again
the compiler can employ for efficiency purposes.

@< Type definitions @>=
template<bool pilfer=false>
  class local_identifier : public identifier
{ size_t depth, offset;
public:
  explicit local_identifier(id_type id, size_t i, size_t j)
     : identifier(id), depth(i), offset(j) @+{}
  virtual void print(std::ostream& out) const;
  virtual void evaluate(level l) const; // only this method is redefined
};

@ The method |local_identifier::evaluate| looks up a value in the evaluation
context |frame::current|, by calling the method |evaluation_context::elem|.

@< Function definitions @>=
template<>
  void local_identifier<false>::print(std::ostream& out) const
@+{@; out << name(); }
template<>
  void local_identifier<true>::print(std::ostream& out) const
@+{@; out << '$' << name(); }

template<bool pilfer>
  void local_identifier<pilfer>::evaluate(level l) const
{
  if (pilfer)
    push_expanded(l,std::move(frame::current->elem(depth,offset)));
  else
    push_expanded(l,frame::current->elem(depth,offset));
}

@*1 Type-checking applied identifiers.
%
When type-checking an applied identifier, we first look in
|layer::lexical_context| for a binding of the identifier; if found it will be
a local identifier, and otherwise we look in |global_id_table|. If found in
either way, the associated type must equal the expected type (if any), or be
convertible to it using |coerce|.

Our answer to the question what to do if the identifier has a type |*id_t| that
is more general than the type~|tp| required by the context has changed over
time. Originally it was allowed, and had as side effect that the identifier type
was specialised to |tp| for future use. That is incorrect since the identifier
may have been used before in a setting where the specialisation would not be
allowed. Then this situation was made illegal (a type error), but finally we
accept it (due to unification inside |conform_types|) but without side effects;
indeed identifier types never change after the introduction of the identifier.

@< Cases for type-checking and converting... @>=
case applied_identifier:
{ const id_type id=e.identifier_variant;
  const type* id_t; size_t i,j;
  std::ostringstream o;
@/
  const bool is_local=(id_t=layer::lookup(id,i,j))!=nullptr;
  if (is_local or (id_t=global_id_table->type_of(id))!=nullptr)
  { expression_ptr id_expr = @| is_local
    ? expression_ptr(new local_identifier<false>(id,i,j))
    : expression_ptr(new global_identifier<false>(id));
    return conform_types(*id_t,tp,std::move(id_expr),e);
  }
  else if (@[auto* vars=global_overload_table->variants(id)@;@])
  @< See if a unique member of |*vars| matches |tp|, and if so |return| a
     |capture_expression| holding the value of that variant @>
     o << "Undefined identifier '" << main_hash_table->name_of(id) << '\'';
  if (e.loc.file!=Hash_table::empty)
    o << ' ' << e.loc;
  throw expr_error(e,o.str());
@.Undefined identifier@>
}

@ When an identifier is found neither in the current |lexical_context| nor in
the |global_id_table|, we used to flag an error, but as a service to the user we
now try instead to find something from the |global_overload_table| first,
provided there is a unique instance of the identifier there that matches the
type requirement |tp| from the context (which may be no requirement at all).
If the context type does not accept any function type we skip this code (the
overload table can hold only functions), and if it does we distinguish between
accepting functions of any argument type, or making some restriction on the
argument type. This is mainly to be able to give more meaningful error messages
when our attempt fails, but the case of no restrictions also can be handled with
somewhat simpler logic. if the type requirements should rule out all available
variants, then we fall through this code as if there were no overloads at all.

@< See if a unique member of |*vars| matches |tp|, and if so |return| a
   |capture_expression| holding the value of that variant @>=
{ type_expr f_type = gen_func_type.copy();
  if (tp.unify_specialise(f_type))
  { if (tp.assign().is_polymorphic(f_type.func()->arg_type))
    @< If |*vars| has a unique variant, |return| its value wrapped in a
       |capture_expression|, otherwise |throw| an error signalling ambiguity @>
    else
    @< If a unique variant among |*vars| unifies to the type |tp|,
       |return| its value wrapped in a |capture_expression|; if more than one
       does, |throw| an error signalling ambiguity, and if none does
       fall through @>
  }
}

@ When the context poses no requirements, we basically look for a unique variant
for the identifier and return it. However we must specialise |tp| to the
possibly polymorphic function type of |variant|, which |functype_absorb| should
accomplish. If there are multiple variants, we warn the user that one cannot
simply use an overloaded symbol as an identifier with nothing to disambiguate
the usage.

@< If |*vars| has a unique variant, |return| its value wrapped...@>=
{ if (vars->size()==1)
  { const auto& variant = vars->front();
    if (functype_absorb(tp,variant))
    { o << main_hash_table->name_of(id) << '@@' << tp.func()->arg_type;
      return
        expression_ptr(new capture_expression(variant.value(),o.str()));
    }
    throw logic_error
      ("Partially determined function type failed to specialise");
  }
  else
  { o << "Use of overloaded '" << main_hash_table->name_of(id)
    @| << "' is ambiguous, specify argument type to disambiguate";
    throw expr_error(e,o.str());
  }
}

@ In the case where the context poses a restriction, we filter the variants for
one that matches the requirement, and use the result if it is unique. We proceed
basically as in the case of no restriction, but here the call of
|functype_absorb| can meaningfully fail, indicating that the variant does not
satisfy the restriction. In case of such a failure, we must undo any type
assignments that the failed attempt to match may have made. Since this must be
done repeatedly, it is simplest to just apply any assignment that might be
pending initially, by calling |tp.wring_out|, after which any type assignments
we have must be new, and can be undone by calling |tp.clear|. The fact that we
call |tp.wring_out| does mean that callers to |convert_expr| cannot expect that
changes to |tp| are limited to acquiring type assignments; we do not believe any
caller make that assumption, but if that should be the case, the code below must
be more careful, and instead record |tp.polymorphics()| initially, and then use
|restore_polymorphics| instead of |clear|.

@< If a unique variant among |*vars| unifies to the type |tp|... @>=
{ tp.expunge(); // ensure no assignments are pending
  const overload_data* prev_match=nullptr;
  expression_ptr result;
  for (const auto& variant : *vars)
  {
    if (functype_absorb(tp,variant))
    { if (result!=nullptr)
        @< Throw an error reporting ambiguous overloaded symbol usage @>
      o << main_hash_table->name_of(id) << '@@' << tp.func()->arg_type;
      result.reset(new capture_expression(variant.value(),o.str()));
      prev_match = &variant; // record for purpose of possible error message
    }
    tp.clear(); // reset initial type
  }
  if (result!=nullptr)
    return result;
}

@ The code below resembles what we will later do to report an ambiguous
overloaded function call, but here we report the full function types.

@< Throw an error reporting ambiguous overloaded symbol usage @>=
{
  tp.clear(); // forget the type assignment matching current variant
  o << "Ambiguous overloaded symbol " << main_hash_table->name_of(id)
    <<  ", context type " << tp
  @|<< " matches both (" << prev_match->f_tp().arg_type
  @|<< "->" << prev_match->f_tp().result_type
  @|<< ") and (" << variant.f_tp().arg_type
  @|<< "->" << variant.f_tp().result_type << ')';
  throw expr_error(e,o.str());
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

To resolve overloading, we used to try each variant, each time recursively
calling |convert_expr| to convert the arguments under the hypothesis that they
were in the strong type context of the type expected by that variant; this while
catching (type) errors, interpreting them simply as an indication that this
variant fails to match. However in long formulae, operator applications are
nested, and each |convert_expr| would call |resolve_overload| one level down,
leading to a recursive repetition of the same try-all-variants scenario. This
led to a search tree growing exponentially with the size of the formula, and
unacceptably long analysis times. So instead the rules were made slightly
stricter: every operand combination must be such that is can be analysed in
isolation, giving rise to an \foreign{a priori} type, which might however not
quite be the actual type for the intended variant. Overload resolution will then
try to find a variant whose expected type is close enough to the
found \foreign{a priori} type to expect that zero or more coercions will allow
making the match. If any coercions are necessary to do so, they might need to
``creep inside'' the argument expression in order to apply (most obviously so in
the common cases of an operand pair or argument tuple), so that it might not be
possible to use the expression as initially converted; if this appears to be the
case, the conversion will be attempted as second time, this time in the strong
context of the type expected by the matched variant. The second attempt is only
made once all variants fail to make an exact match, and if may still fail (if
the expression structure of the argument does not allow for inserting the
coercion), which will then be considered to be an error rather than just a
failed match. Additional language restrictions are imposed to ensure that such
an error cannot mask a possible successful match of a later variant.

If we do decide to redo the initially converted argument expression, this again
gives a potential exponential growth of the search tree, albeit with base$~2$
rather than the number of variants: both the original and the new call of
|convert_expr| can recursively call |resolve_overload|. With a few specific
adaptations, we avoid useless re-conversions in most cases, ensuring that such
growth is exceedingly unlikely to happen in practice, although it can still be
provoked in artificial examples.

@ The function |resolve_overload| must decide which of the |variants| best
matches the argument expression. It starts by converting the latter without any
imposed type, finding an |a_priori_type|. Then when matching a variant, we can
get an exact match with |a_priori_type|, in which case the variant is selected
and an application of its function is returned. When, lacking an exact match, a
possible match after coercions is suspected, a second conversion is attempted.

We perform (at most) two passes over the variants, first looking for exact
match, then for an inexact match. The first pass includes matching generic
operators, their matching against |a_priori_type| being done by unification. If
this succeeds, it is considered to be an exact match; we will skip over any
polymorphic variants during a second pass if we need to do one. During the
second pass we may detected sufficient type proximity is for a potential inexact
match, and if so, a second conversion of the argument with the expected type for
the variant is done, in the hope that inserted implicit conversions can help
obtaining that type. This may however still fail (the structure of the argument
expression may not allow for insertion of the coercions), and if this happens
for some variant, we give up without testing any more variants. If the second
pass completes without finding any potential match, we report a ``failed to
match'' error.

@:resolve_overload@>

@< Local function definitions @>=

@< Define the function |equals_name| @>
expression_ptr resolve_overload
  ( const expr& e
  , type& tp
  , const overload_table::variant_list& variants
  )
{
  const expr& args = e.call_variant->arg;
  const auto n_args = args.kind==tuple_display ? length(args.sublist) : 1;
  auto tup_exp = std::make_unique<tuple_expression>(n_args);
  std::vector<expression_ptr>& arg_vector = tup_exp->component;
@)
  type a_priori_type = type::bottom(tp.floor()); bool is_constant=false;
@/@< Fill |arg_vector| with expressions converted from the components of |args|,
     and set |a_priori_type| to the type found; also set |is_constant=true|
     if there are multiple arguments and they are all constant expressions @>
@)
  id_type id =  e.call_variant->fun.identifier_variant;
  a_priori_type.expunge(); // make sure any type assignments are baked in
  @< Try to find an element of |variants| whose |arg_type| exactly matches
     |a_priori_type|, and for it |return| a corresponding call with argument
     from |tup_exp| @>
@) // no exact match; retry all cases for an inexact match
  @< Try to find an element of |variants| whose |arg_type| matches
     |a_priori_type| after inserting implicit type conversions, and for
     it |return| a corresponding call with argument from |tup_exp| @>
@)
  @< Complain about failing overload resolution @>

}

@ The following code rather uncharacteristically treats all argument lists as
tuples, even when |n_args==1| indicates that no tuple expression is actually
present. The reason that we do not simply call |convert_expr| for the entire
argument tuple, is that we later want to be able to easily reconvert one
argument expression at a time, and it would be hard to get access to these
individual arguments if the argument tuple had been converted to a single
|expression_ptr|. Since we are thus bypassing the normal conversion of tuple
expressions, we also get the responsibility for doing constant folding across
them (which is otherwise done as a part of converting tuple expressions). This
in fact allows us to postpone constant folding until the correct types have been
found, and any required coercions are inserted.

The variable |is_constant| is set for later use, namely to decide whether before
building the call we should convert an argument that is a tuple of constants
(represented as denotations) into a constant tuple; since this question does not
arise in case of a single argument (and the conversion should not be attempted)
we set the variable to |false| in that case.

@< Fill |arg_vector| with expressions converted from... @>=
{ is_constant = n_args>1;
  // will record whether all |n_args>1| arguments are denotations
  if (n_args==1)
    arg_vector[0] = convert_expr(args,a_priori_type);
    // get \foreign{a priori} type for argument
  else
  {
    sl_list<type> comp_types;
    auto arg_vector_it = arg_vector.begin();
    for (wel_const_iterator it(args.sublist); not it.at_end();
         ++it,++arg_vector_it)
    {
      type& comp_tp = comp_types.push_back(type::bottom(tp.floor()));
      *arg_vector_it = convert_expr(*it,comp_tp);
      is_constant = is_constant and
         dynamic_cast<const denotation*>(arg_vector_it->get()) != nullptr;
    }
    a_priori_type = type::wrap_tuple(std::move(comp_types));
      // combine argument types
  }
}

@ It may be that more than one variant can produce a match, in which case we
always prefer an exact match if there is one. Also, an exact match should be the
unique such match, which the code below detects and refuses at the second exact
match. Before going on to look for a possible second match, we convert and set
aside a call for the first match, and leave a pointer |prev_match| to the
current variant, which will be used for error reporting in case an ambiguity is
found.

Whenever an argument matches, be it exact of inexact, |conform_types| will be
called for the variant result type and the expected |tp|; this may throw an
error, also aborting the matching process.

@< Try to find an element of |variants| whose |arg_type| exactly matches...@>=
{
  const unsigned int apt_deg = a_priori_type.degree();
  expression_ptr result; // buffer for storage of result
  const overload_data* prev_match=nullptr;
  for (const auto& variant : variants)
  {
    unsigned int shift_amount;
    if (a_priori_type.matches
         (variant.f_tp().arg_type,variant.poly_degree(),shift_amount))
    { // exact match
      if (prev_match!=nullptr)
        @< Throw an error reporting an ambiguous exact match @>
      expression_ptr call;
      expression_ptr arg = n_args==1 ? std::move(arg_vector[0])
			 : expression_ptr(std::move(tup_exp));
      const type_expr& arg_type = a_priori_type.unwrap();
        // actual argument type
@/    @< Assign to |call| a converted call expression of the function value
        |variant.value()| with argument |arg|, which is of type |arg_type| @>
      const type res_type = type::wrap @|
        ( a_priori_type.assign().substitution
            (variant.f_tp().result_type ,shift_amount)
        , a_priori_type.floor()
        );
      result = conform_types(res_type,tp,std::move(call),e);
@/    prev_match = &variant;
@/ // |res_type| is recorded in |tp|, and |arg_type| has served and is forgotten
    }
    a_priori_type.clear(apt_deg);
      // remove type variable introduced by |match|
   }
@)
   if (result!=nullptr) // then a unique match was found
     return result;
}

@ Here we list the matching types, which is easy due to the |prev_match|
pointer.
@< Throw an error reporting an ambiguous exact match @>=
{
  a_priori_type.clear(); // forget the type assignment matching current variant
  std::ostringstream o;
  o << "Ambiguous argument in function call, argument type " << a_priori_type
  @|<< " matches both " << prev_match->f_tp().arg_type
  @|<< " and " << variant.f_tp().arg_type;
  throw expr_error(e,o.str());
}

@ We shall compare frequently for the `\.=' operator name, so it pays to look up
its identifier code once and for all.

@h "lexer.h" // for |main_hash_table|

@< Define the function |equals_name| @>=
id_type equals_name()
{@; static id_type name=main_hash_table->match_literal("=");
  return name;
}

@ Having found a match (exact or inexact; the code below is included twice), and
having converted the argument expression to |arg|, we need to construct a
function call object. The overload table contains an already evaluated function
value~|variant.value()|, which is has type |shared_function|, more specific than
|shared_value| as it points to a value that represents a function object. Such
values can either hold a built-in function, or a user-defined function together
with its evaluation context, and we shall add to the list of possibilities
later. Each type has a corresponding specialised function call type, and we
shall provide a virtual method |function_base::build_call| that binds an
argument expression to the function object to form a complete call; this renders
the expression-building part of code below particularly simple. One noteworthy
point is that |build_call| may need to store a shared pointer back to the
|function_base| derived object that created it, and we therefore pass the shared
pointer |variant.value()| that gave us access to the |function_base| as first
argument to its |build_call| virtual method. The underlying raw pointer equals
|this| of that object, but we need a |std::shared_ptr| so that we also transfer
shared ownership.

As a special safety measure against the easily made error of writing `\.='
instead of an assignment operator~`\.{:=}', we forbid converting to void the
result of an (always overloaded) call to the equality operator, treating this
case as a type error instead. In the unlikely case that the user defines an
overloaded instance of `\.=' with void result type, calls to this operator
will still be accepted.

This is one of the places where we might have to insert a |voiding|, in the rare
case that a function with void argument type is called with a nonempty argument
expression of void type. An alternative solution would be to replace such a call
by a sequence expression, evaluating the argument expression separately and then
the function call with an empty argument expression. In any case the
test \emph{can} be made here, since we have the argument type in the variable
|arg_type|, and the argument expression in |args|.

@< Assign to |call| a converted call expression of the function value
   |variant.value()|... @>=
{ if (tp.is_void() and
      id==equals_name() and
      variant.f_tp().result_type!=void_type)
  { std::ostringstream o;
    o << "Use of equality operator '=' in void context; " @|
      << "did you mean ':=' instead?\n  If you really want " @|
      << "the result of '=' to be voided, use a cast to " @|
      << variant.f_tp().result_type << '.';
    throw expr_error(e,o.str());
  }
@)
  // now select unique component or consolidate tuple
  std::ostringstream name;
  name << main_hash_table->name_of(id) << '@@' << arg_type;
@)
  if (arg_type==void_type)
@/{@; if (not is_empty(args))
      arg.reset(new voiding(std::move(arg)));
  }
  else if (is_constant)
    make_row_denotation<false>(arg); // wrap tuple inside a denotation
  call = variant.value()->build_call
           (variant.value(),name.str(),std::move(arg),e.loc);
}

@ Inexact matches are only considered for variants with a completely specific
(i.e., monomorphic) type; moreover, when left with a choice among several
inexact matches, we prefer one that invokes the fewest coercions. These priority
rules are implemented by ordering our tests so that the first match is also the
best; for this purpose the variants have been partially sorted to ensure this.
To test for an inexact type match, we call |is_close| with |a_priori_type| and
the expected type of a variant.

Using |is_close| here, rather than |coercible|, means that the possibility of
voiding coercions is not taken into account: when |arg_type==void_type|, the
only match accepted is when no arguments are present (though pedantically
speaking, an argument with a priori type void also matches); this is intended
behaviour.


@< Try to find an element of |variants| whose |arg_type| matches... @>=
for (const auto& variant : variants)
  if (not variant.is_polymorphic() and
      (is_close(a_priori_type.unwrap(),variant.f_tp().arg_type)&0x1)!=0)
    // inexact match
  { const type_expr arg_type = variant.f_tp().arg_type.expanded();
    expression_ptr arg; // will hold the final converted argument expression
    @< Apply an implicit conversion for every argument whose type in
       |a_priori_type| differs from that required in |arg_type|, then assign
       converted expression to |arg| @>
    expression_ptr call;
@/  @< Assign to |call| a converted call expression of the function value
       |variant.value()| with argument |arg|, which is of type |arg_type| @>
    const type res_type = type::wrap(variant.f_tp().result_type,tp.floor());
    return conform_types(res_type,tp,std::move(call),e);
  }


@ When we come here, an inexact match is found; we traverse the arguments (maybe
a single one) of the call, and modify only those whose a priori type does not
match the type expected for that argument of the match. By treating arguments
independently, we mostly avoid a potential inefficiency of our approach that had
previously been present for years. When multiple arguments were treated as a
single tuple-display argument, it would require being entirely converted anew if
\emph{any} component did not produce the exact type required: this could lead to
repeated conversion of another argument whose type does not need conversion. In
nested formulae, like an explicit repeated vector addition where all argument
vectors need conversion from integer lists, this could lead to processing time
increasing exponentially with the formula size.

The code below avoids a new call of |convert_expr| in two circumstances, namely
when the type found for this argument is already an exact match for the
corresponding position of |arg_type| (which can only happen if |n_args>1|), and
when the form of the argument is such (for instance a formula or function call)
that implicit conversion is only possible on the outside of the expression
(contrary to for instance list display or conditional expressions where the
conversion can ``creep into'' the expression itself). The condition on the form
of the expression that allows avoiding a call of |convert_expr| is detailed in a
separate module called |@< The form of |cur_arg| shields... @>|. The fact that
we do these tests makes a scenario with exponential increase of processing time,
thought theoretically still possible, exceedingly unlikely.

In the case of a single argument there is no tuple expression for the arguments.
Nonetheless |arg_vector| will be a vector of length~$1$ in this case, but
(probably) none of |a_priori_type|, |arg_type| and |args| hold vectors or lists
at all, so we could not define the iterators used below in the case of more than
one argument. We therefore separate the two cases, even though the actions are
the same.

@< Apply an implicit conversion for every argument... @>=
{
  if (n_args==1)
  { const expr& cur_arg = args; // rename to allow sharing the following module:
    if (@< The form of |cur_arg| shields out the context type @>@;@;)
    { arg = std::move(arg_vector[0]);
      if (not coerce(a_priori_type.unwrap(),arg_type,arg,e.loc))
        throw type_error(e,a_priori_type.unwrap().copy(),arg_type.copy());
      }
    else arg =
      convert_expr_strongly(cur_arg,tp.floor(),arg_type); // redo conversion
  }
  else // |n_args!=1|
  { assert(a_priori_type.kind()==tuple_type
       and arg_type.raw_kind()==tuple_type);
    wtl_const_iterator apt_it(a_priori_type.tuple())
    , argt_it(arg_type.tuple());
    wel_const_iterator exp_it(args.sublist);
    for (auto res_it = arg_vector.begin(); res_it!=arg_vector.end();
       ++res_it,++apt_it,++argt_it,++exp_it)
      if (*apt_it!=*argt_it) // only act if a priori component type is not exact
      { const expr& cur_arg = *exp_it;
        if (@< The form of |cur_arg| shields out the context type @>@;@;)
        { if (not coerce(*apt_it,*argt_it,*res_it,e.loc))
            throw type_error(e,apt_it->copy(),argt_it->copy());
        }
        else
          *res_it = convert_expr_strongly(cur_arg,tp.floor(),*argt_it);
          // redo conversion
      }
    arg = expression_ptr(std::move(tup_exp));
    // wrap up modified tuple expression
  }
}

@ Many expression types listed in |enum expr_kind@;| have an interpretation that
is insensible to the type expected by the context, so that any type mismatch can
only be accommodated by a conversion applied externally. For these expression
type we therefore avoid calling |convert_expr| again. This is a long logical
disjunction, where the main cases that might be quite large subexpressions are
mentioned first: function applications (which includes formulae), assignments
statements and casts. Some cases mentioned here are only for conceptual
completeness, as they can only apply in erroneous situations: for instance no
implicit conversion can apply to expressions of type |bool|, so no
|boolean_denotation| or |negation_expression| can ever come here (with an
inexact match). A clever compiler could transform the test below into testing a
bit at position |cur_arg.kind| in a fixed bitset.

@< The form of |cur_arg| shields out the context type @>=
cur_arg.kind==function_call @| or
cur_arg.kind==ass_stat @| or
cur_arg.kind==comp_ass_stat @| or
cur_arg.kind==field_ass_stat @| or
cur_arg.kind==cast_expr @| or
cur_arg.kind==op_cast_expr @| or
cur_arg.kind==subscription @| or
cur_arg.kind==slice @| or
cur_arg.kind==applied_identifier @| or
cur_arg.kind==integer_denotation @| or
cur_arg.kind==string_denotation @| or
cur_arg.kind==boolean_denotation @| or
cur_arg.kind==lambda_expr @| or
cur_arg.kind==rec_lambda_expr @| or
cur_arg.kind==last_value_computed @| or
cur_arg.kind==negation_expr

@ Here is the final part of |resolve_overload|, reached when no valid match
could be found. In that case we |throw| a |expr_error| explaining the operator
or function identifier, and the a priori type of the operand. Most function
definitions will be in the overload table even if just one definition is
present; in the latter case the ``Failed to match'' error might seem
unnecessarily vague, so we produce instead a more specific |type_error|, whose
message will also mention that the type that was expected by the unique
instance.

@< Complain about failing overload resolution @>=
if (variants.singleton())
  throw type_error(args,a_priori_type.unwrap().copy(),
                   variants.front().f_tp().arg_type.copy());
else
{ std::ostringstream o;
  o << "Failed to match '"
    << main_hash_table->name_of(id) @|
    << "' with argument type "
    << a_priori_type;
  throw expr_error(e,o.str());
}

@* Function calls.
%
One of the most basic and important tasks of the evaluator is to allow
function calls, which may involve either built-in or user-defined functions;
to this list may be added other types of runtime values that behave like
functions, such as selectors from a tuple.

This central part of the evaluator will be presented with an initial focus on
built-in functions, while leaving the particulars of user defined functions
(also known as $\lambda$-expressions) and other types aside until somewhat
later. This corresponds more or less to the development history of the
interpreter, in which initially only built-in functions were catered for;
however many of the aspects that we deal with right away, notably function
overloading, are in fact much more recent additions than user-defined
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
  : expression_base(), argument(arg.release()), loc(loc) @+{}
  virtual ~call_base() = default;
  virtual std::string function_name() const=0;
};
@)

@ We similarly define an intermediate class |function_base| between |value_base|
and the concrete classes that will define function objects, like built-in
functions. The main (virtual) methods introduced here are |apply|, which will
serve in |call_expression::evaluate| below to implement a call of the function
object once arguments have been evaluated to the stack, and |build_call| that is
instead used to build a specialised call expression when a function value is
identified at analysis time (in overloaded calls). In addition |argument_policy|
tells how the function object wants its arguments prepared, |maybe_push| is a
hook that does nothing except for recursive functions that use it for their
implementation, and |report_origin| which serves in forming an back-trace in
case of errors during execution of the function.

@< Type def... @>=
// \.{global.h} predeclares |function_base|, and defines:
// |typedef std::shared_ptr<const function_base> shared_function;|
@)
struct function_base : public value_base
{
  function_base() : value_base() @+{}
  virtual ~function_base() = default;
  static const char* name() @+{@; return "function value"; }
@)
  virtual void apply(eval_level l) const=0;
    // go; arg.\ values on |execution_stack|
  virtual eval_level argument_policy() const=0;
    // form to prepare arguments in
  virtual void maybe_push(const std::shared_ptr<const function_base>& p) const
    @+{}
  virtual void report_origin(std::ostream& o) const=0;
    // tell where we are from
  virtual expression_ptr build_call
    (const shared_function& master,const std::string& name,
     expression_ptr&& arg, const source_location& loc) const=0;
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
  virtual ~call_expression() = default;
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
{ if (dynamic_cast<const identifier*>(function.get())!=nullptr)
    out << *function;
  else out << '(' << *function << ')';
  if (dynamic_cast<const tuple_expression*>(argument.get())!=nullptr)
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

Some built-in functions like |print| accept arguments of any types, and in
particular tuples of any length. For such functions there is no use in
adopting the approach used for other built-in functions of expanding argument
tuples on the stack; instead the argument is always considered as one value.
We cater for the distinction between the two variants at run-time by
making this a class template with a Boolean template argument |variadic|. This
is necessary because operator casts make it possible to use specialisations
of variadic functions as function values (of the correspondingly specialised
type), so they can occur not only in overloaded calls, but in any place that
other built-in or user-defined functions can.

As another supplementary information, we store an indication of whether this
built-in function wants to produce its result by modifying one of its arguments;
for instance, |suffix_element_wrapper| likes to modify the row value that is its
first element in-place by adding a suffix. That operation can be done without
creating a copy of the argument only if it is not shared in any way, and in
certain cases we want to be able to optimise performance by arranging for that
argument to not be shared if it can be avoided. The |hunger| field tells whether
this is such a function: a value $0$ means no desire to eat anything, a value
$1$ means it wishes to modify the first of two arguments in-place, a value of
$2$ that it wishes to do so for the second argument, and a value of $3$ means
that there is a unique argument, which it wishes to modify.

@< Type definitions @>=
template <bool variadic>
  struct builtin_value : public function_base
{ wrapper_function val;
  std::string print_name;
  unsigned char hunger; // (always |0| when |variadic| holds)
@)
  builtin_value(wrapper_function v,const std::string& n, unsigned char hunger)
  : function_base(), val(v), print_name(n), hunger(hunger) @+ {}
  virtual ~builtin_value() = default;
  virtual void print(std::ostream& out) const
  @+{@; out << '{' << print_name << '}'; }
  virtual void apply(eval_level l) const @+
    {@; (*val)(l); } // apply function pointer
  virtual eval_level argument_policy() const
  {@; return variadic
      ? eval_level::single_value : eval_level::multi_value; }
  virtual void report_origin(std::ostream& o) const @+
  {@; o << "built-in"; }
  virtual expression_ptr build_call
    (const shared_function& master,const std::string& name,
     expression_ptr&& arg, const source_location& loc) const;
@)
  static const char* name() @+{@; return "built-in function"; }
  builtin_value(const builtin_value& v) = delete;
};

@ While syntactically more complicated than ordinary function calls, the call
of overloaded functions is actually more direct at run time, because the
function is necessarily referred to by an identifier (or operator) instead of
by an arbitrary expression, and overloading resolution results in that
identifier being replaced by a function \emph{value}, known at analysis time.
Different kinds of function values, derived from |function_base|, can then
give rise to different kind of overloaded calls. We introduce
|overloaded_call| as another abstract class between |call_base| and those
kinds of call; its main purpose is to store a name that reflects how
overloading was resolved.
@< Type definitions @>=
struct overloaded_call : public call_base
{ std::string name;
@)
  overloaded_call
    (const std::string& name,
     expression_ptr&& arg,
     const source_location& loc)
@/: call_base(std::move(arg),loc), name(name) @+ {}
  virtual ~overloaded_call() = default;
  virtual std::string function_name() const @+{@; return name; }
  virtual void print(std::ostream& out) const;
};

@ This base class already has the elements necessary for printing the call.

@< Function definitions @>=
void overloaded_call::print(std::ostream& out) const
{ out << name;
  if (dynamic_cast<const tuple_expression*>(argument.get())!=nullptr)
     out << *argument;
  else out << '(' << *argument << ')';
}

@ If its function value is a |builtin_value|, an overloaded call will become an
|overloaded_builtin_call| rather than a |call_expression|. The
|overloaded_builtin_call| stores a pointer |f| back to |builtin_value| that
created it, which is useful to identify the actual function when an error is
thrown during its execution. To ensure that the |builtin_value| remains in
existence as long as the |overloaded_builtin_call| is, this is in fact a shared
pointer. (This back reference might have been avoided by copying relevant
information from the |builtin_value| into the |overloaded_builtin_call|, but a
given built-in function might generate many calls, so it would be wasteful to do
so.) What is most frequently needed about concerning the |builtin_value| is its
stored function pointer, which is needed every time the
|overloaded_builtin_call| is executed, so we do copy this pointer into a member
|f_ptr| of the latter, which avoid an extra indirection when using it.

When accessed through overloading, the condition whether a built-in function
is variadic or not is known at compile time, so we can make this a class
template with a |variadic| template argument, just like |builtin_value|;
indeed the template argument serves exclusively to adapt the type of the
member~|f|. An additional constructor without |name| argument is provided for
convenience to places where our interpreter directly produces calls to
certain built-in functions, without going through the overload table; the
function name is then deduced from the one stored in the |builtin_value|.
The first constructor take the shared pointer |fun| is by value since it is most
often, but not always, held in a local variable of the caller from which it can
be moved; for the second constructor we pass by constant reference since the
argument will not come from a local variable.

@< Type definitions @>=
template <bool variadic>
  struct overloaded_builtin_call : public overloaded_call
{ typedef std::shared_ptr<const builtin_value<variadic> > ptr_to_builtin;
@)
  ptr_to_builtin f;
   // points to the full |builtin_value|, for back-tracing interrupted calls
  wrapper_function f_ptr; // shortcut to implementing function
@)
  overloaded_builtin_call
    (ptr_to_builtin fun,
     const std::string& name,
     expression_ptr&& arg,
     const source_location& loc)
@/: overloaded_call(name,std::move(arg),loc)
  , f(std::move(fun)), f_ptr(f->val) @+ {}
  overloaded_builtin_call
    (const ptr_to_builtin& fun,
     expression_ptr&& arg,
     const source_location& loc)
@/: overloaded_call(fun->print_name,std::move(arg),loc)
  , f(fun), f_ptr(fun->val) @+ {}
  virtual ~overloaded_builtin_call() = default;
  virtual void evaluate(level l) const;
};
@)
typedef overloaded_builtin_call<false> builtin_call;
typedef overloaded_builtin_call<true> variadic_builtin_call;

@ A |builtin_value| can turn itself into an |overloaded_builtin_call| when
provided with a shared pointer |master| to itself, an argument expression, a
|name| to call itself and a |source_location| for the call. The
|static_pointer_cast| reverts |master| back to the type it had before it was
up-cast in the caller to the pointer-to-base type |shared_function|; one cannot
make a virtual method argument of covariant pointer type.

@< Function def... @>=
template <bool variadic>
expression_ptr builtin_value<variadic>::build_call
    (const shared_function& master,const std::string& name,
     expression_ptr&& arg, const source_location& loc) const
{ assert(master.get()==this);
  auto f=std::static_pointer_cast<const builtin_value<variadic> >(master);
  return expression_ptr(new @|
    overloaded_builtin_call<variadic>(f,name,std::move(arg),loc));
}

@ Sometimes we want a builtin operator to have special behaviour at compile
time, in the sense that when an application is formed with the global overload
table matching the |builtin_value| for this operator, the virtual |build_call|
method inspects the arguments and in some cases replaces the operator by another
builtin function with a part of the arguments; the typical example is an
application $E+1$ that transforms itself into |succ(E)|. To accommodate this we
define a class derived from |builtin_value| that provides an alternative |apply|
method that tries this substitution before reverting to
|builtin_value::build_call| if special argument values were not found.

@< Type definitions @>=
struct special_builtin : public builtin
{
  typedef expression_ptr (*tester)
    (expression_ptr&,const shared_builtin&,const source_location&);
  using test_data = std::pair<tester,shared_builtin>;
  sl_list<test_data> tests;
@)
  special_builtin
    (wrapper_function v,const std::string& n, unsigned char hunger)
  : builtin(v,n,hunger), tests() @+{}
  virtual ~special_builtin() = default;

  virtual expression_ptr build_call
    (const shared_function& master,const std::string& name,
     expression_ptr&& arg, const source_location& loc) const
  {  for (const auto& test : tests)
     { auto p = (*test.first)(arg,test.second,loc);
       if (p!=nullptr)
         return p;
     }
     return builtin::build_call(master,name,std::move(arg),loc);
  }
};

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
template<>
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
  {@; (*f_ptr)(l); } // call the built-in function
  @< Catch block for exceptions thrown within call of |f| with |arg_string| @>
}

@ To provide back-trace, we catch and re-throw an error after extending the
stored error string. The result is a list of interrupted named function calls,
from inner to outer.

The work of modifying the error string is common to several such |catch|
blocks, and relegated to a function |extend_message| to be defined presently.
The error string is modified within the existing error object; this is a
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

@< Catch block for exceptions thrown within call of |f| with |arg_string| @>=
catch (error_base& e)
{@; extend_message(e,this->function_name(),loc,f,arg_string);
  throw;
}
catch (const std::exception& e)
{ runtime_error new_error(e.what());
  extend_message(new_error,this->function_name(),loc,f,arg_string);
  throw new_error;
}

@ The function |extend_message| facilitates appending information to error
messages in |catch| blocks. It is called with, apart from the error~|e| whose
message is to be modified, the expression~|call| whose evaluation was
interrupted by the error, the value |f| of the function called (either a
|builtin_value| or a |closure_value|), and a string~|arg| that in debug mode
describes the arguments (when not in debug mode the string will be empty and
is ignored).

We report the source location of the call expression and the name of the
function called (both obtained from |call|), and a source location for the
definition of the called function (obtained from |f|) in case it is
user-defined; when the called function was built in we just report that.

@h "parsetree.h" // for output of |source_location| value

@< Local fun... @>=
void extend_message
  (error_base& e,
   const std::string& name, const source_location& loc,
   shared_function f,
   const std::string& arg)
{ std::ostringstream o;
  o << "In call of " << name << ' ' << loc << ", ";
  f->report_origin(o);
  o << '.';
  if (verbosity>0)
    o << "\n  argument" << (arg[0]=='(' ? "s: " : ": ") << arg;
  e.trace(o.str());
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
template<>
void builtin_call::evaluate(level l) const
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
  {@; (*f_ptr)(l); } // call the built-in function
  @< Catch block for exceptions thrown within call of |f|... @>
}


@ Finally we consider the case where evaluation is done though a
|call_expression|, as happens when the function object cannot be identified at
compile time. Since the function to be called is here produced by evaluating an
expression (maybe as simple as an identifier), the question whether a built-in
rather than user-defined (or other kind of) function is being called can only be
determined at run time. The distinction is handled though the virtual methods
|argument_policy| and |apply| of |functions_base|, so the code below handles all
kinds of function objects at once.

To evaluate a |call_expression| object, we first evaluate this function part,
and then dynamically cast it to a |function_base f@;| (which should always
succeed, assuming type checking is sound). Then we evaluate the arguments
according to the |f->argument_policy()|, and then call |f->apply|, which should
find the arguments in the form it expects.
In that case we pass on the |level| parameter that was passed to
|call_expression::evaluate| method to the built-in function, so that if
necessary it can in its turn return and expanded result (or no result at all).
The evaluation of user-defined functions will be detailed later, but we can
already say that in this case it will be more useful to receive the argument on
the stack as a single value.

We reuse the previous |catch| block literally a third time; this time not only
do we judiciously choose the name |arg_string| to match what we did before,
but also the local variable name |f| to math the field name
|builtin_call::f| that the cited module referred to in previous
instances.

@: general call evaluate @>
@< Function definitions @>=
void call_expression::evaluate(level l) const
{ function->eval();
  auto f = std::dynamic_pointer_cast<const function_base>(pop_value());
  if (f.get()==nullptr)
    throw logic_error
      ("Non-function value found for function in a call expression");
  f->maybe_push(f); // provide shared pointer to self for recursive functions
  std::string arg_string;
  if (verbosity==0)
    argument->evaluate(f->argument_policy());
  else
  { auto sp = execution_stack.size();
    argument->evaluate(f->argument_policy());
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
  try {@; f->apply(l); } // apply the function, handling |l| appropriately
  @< Catch block for exceptions thrown within call of |f|... @>
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
    @< Convert and |return| an overloaded function call if |call.fun| is known
       in |global_overload_table|,
       unless it is a local identifier with function type @>
  type f_type=type::bottom(fc); // start with completely generic type
  expression_ptr fun = convert_expr(call.fun,f_type);
    // consolidate to a (maybe polymorphic) |type|
  type_expr f_tp = gen_func_type.copy();
  if (not f_type.unify_specialise(f_tp))
    @< Complain that function part of |e| has non-function type |f_type| @>
  type arg_type = type::wrap(f_tp.func()->arg_type,fc);
  expression_ptr arg = convert_expr(call.arg,arg_type);
  if (arg_type.is_void() and not is_empty(call.arg))
    arg.reset(new voiding(std::move(arg)));
  if (not f_type.matches_argument(arg_type))
    throw type_error(e,arg_type.bake_off(),std::move(f_type.func()->arg_type));
  expression_ptr re(new @|
     call_expression(std::move(fun),std::move(arg),e.loc));
  const type result_type = type::wrap
    (f_type.assign().substitution(f_type.func()->result_type) ,fc);
  return conform_types(result_type,tp,std::move(re),e);
}

@ Rather than imposing a generic function type when converting the function part
of a call, we let the conversion freely determine the type, and only detect the
problem when that type fails to unify with generic function type. This makes it
possible to give a more specific error message, which we do here.

@< Complain that function part of |e| has non-function type |f_type| @>=
{ std::ostringstream o;
  o << "Expression in function position has non-function type " << f_type;
  throw expr_error(call.fun,o.str());
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

@< Convert and |return| an overloaded function call... @>=
{ const id_type id =call.fun.identifier_variant;
  size_t i,j; // dummies; local binding not used here
  auto local_type_p=layer::lookup(id,i,j);
  if (local_type_p==nullptr or local_type_p->top_kind()!=function_type)
     // not calling by local identifier
  { if (@[const auto* variants = global_overload_table->variants(id)@;@])
      return resolve_overload(e,tp,*variants);
  } // |else| fall through to try non-overloaded call
}

@*1 Support for constant folding.
%
Our interpreter is based on the conversion by |convert_expr| of an |expr|
representing the abstract syntax tree into an executable form as an
|expression_ptr|, which is then executed. But there are some points where we
want to perform parts of the evaluation during the conversion process itself, so
that constant subexpressions can be produced from a |denotation| even if not
written as one in the input (which is possible only for natural numbers, Boolean
values, and strings).

As a small intermezzo, we define a function that implements constant folding
across implicit conversions, in other words that ensures that if an implicit
conversion is applied to a constant expression (represented as a |denotation|),
then the conversion is done during analysis of the program (compile time) rather
then being postponed to run time, as it must be if the argument is non constant.
A similar function |do_builtin| is called to test for certain built-in functions
whether their arguments are constant; unlike |do_conversion| which builds an
expression in all cases, |do_builtin| does nothing if the arguments are not
constant, so it returns a Boolean telling whether it did transform a denotation.

The function |frozen_error| is an auxiliary for these transformations, to
transform errors that might occur at compile time into expressions that postpone
the error until run time, reproducing the error message.

@< Declarations of exported functions @> =
void do_conversion
  (expression_ptr& e, const conversion_info& ci, const source_location& loc);
bool do_builtin
  (expression_ptr& e, wrapper_function f, const source_location& loc);
expression_ptr frozen_error(std::string message, const source_location& loc);

@ Implementing |do_conversion| is not hard, but does imply some mixing of stages
in the evaluation process. First we need to look into the expression that the
|e| points to using a |dynamic_cast|, to see whether its actual type is
|denotation|. If so, we apply the conversion to the value held inside, for which
that value is temporarily placed on the execution stack. Placing the converted
value back into the |denotation| requires casting away the |const| in the type
of |den_ptr|, inherited from that of |e|. This could have been avoided by
wrapping a fresh |denotation| around the new value and assigning that to |e|; it
is a non-|const| reference, so that would have been possible without trickery.

However, the current solution is both simpler and more efficient, and is really
an indication that the old decision to make |expression_ptr| a pointer-to-const
results in resistance against the evolution of the language. That choice
reflects the fact that the \emph{evaluation process} never needs to modify the
tree of executable expressions, but as the implementation involves, it is now
becoming more and more common to restructure the expression tree after its
initial construction, at compile time. Apart from needing non-|const| access to
the |denotation| structure, the code below also uses the fact that its
|denoted_value| field is not declared |const|; if it were not for the presence
of the functions here, that would have been a reasonable thing to do.

Since we are performing some evaluation during compile time here, we must
consider the possibility that this produces a ``runtime'' error. Our solution is
(for the moment) to not throw a |runtime_error| during the type analysis, but to
compile in, using |frozen_error|, a call to |error_builtin| reproducing the
error upon evaluation. Maybe at some future point we shall decide that it is
actually preferable to already emit a warning at compile time when this
happens, rather than silently compiling a time bomb.

@< Function def... @> =
void do_conversion(expression_ptr& e, const conversion_info& ci,
  const source_location& loc)
{
  if (@[auto* den_ptr = dynamic_cast<const denotation*>(e.get())@;@])
  { try
    {
      auto& den_val = const_cast<denotation*>(den_ptr)->denoted_value;
      push_value(std::move(den_val));
      ci.convert();
      den_val=pop_value();
    }
    catch(std::exception& err)
    {@;
      e = frozen_error(err.what(),loc);
    }
  }
  else
    e.reset(new conversion(ci,std::move(e)));
}

@ The implementation of |do_builtin| is very similar, including the fact that it
replaces |e| by a call to |error_builtin| if a runtime error if one was thrown
by calling~|f|. A difference is the if |e| was not a constant expression, we
simply return |false|. We also need to take into account the fact that builtin
functions with multiple arguments expect them to be present as separate values
on the execution stack, while the |denotation| holds them as a single tuple;
therefore we need to apply |push_expanded| with |multi_value| to it to prepare
to calling the wrapper function~|*f| of the built-in. The latter call should in
all cases produce a single value, since it needs to be reinserted at the
location |den_val| where the arguments used to be.

@< Function def... @> =
bool do_builtin(expression_ptr& e, wrapper_function f,
  const source_location& loc)
{
  if (@[auto* den_ptr = dynamic_cast<const denotation*>(e.get())@;@])
  { try
    {
      auto& den_val = const_cast<denotation*>(den_ptr)->denoted_value;
      push_expanded(eval_level::multi_value,std::move(den_val));
@/    (*f)(eval_level::single_value);
@/    den_val=pop_value();
    }
    catch(std::exception& err)
    {@;
      e = frozen_error(err.what(),loc);
    }
    return true;
  }
  return false;
}

@ The implementation of |frozen_error| is straightforward: we make a
|shared_value| for the string, wrap in into a |denotation|, and then build and
return a call to |error_builtin| with this string as argument.

@< Function def... @> =
@)
expression_ptr frozen_error(std::string message, const source_location& loc)
{
  auto mess_val = std::make_shared<string_value>(std::move(message));
  expression_ptr arg(new denotation(std::move(mess_val)));
  return expression_ptr (new variadic_builtin_call
@|  (error_builtin,"error@@string",std::move(arg),loc) );
}



@* Let-expressions, and identifier patterns.
%
We shall now consider a simple type of expression in which local variables
occur, the \&{let}-expression. It is equivalent to an anonymous function (or
$\lambda$-expression) applied to the explicitly given initialiser expression(s).
But the presence of these expressions makes it unnecessary to declare types
explicitly for the variable(s) introduced, since they can be taken to be the
types of the corresponding initialiser expressions. Nevertheless,
let-expressions could be represented internally as a call in which the function
is a $\lambda$-expression, thus hiding the syntactic origin of the expression;
indeed this was our implementation for a long time. Currently we instead use a
separate internal representation as |let_expression|, whose evaluation is
somewhat more efficient than that of the equivalent $\lambda$-expression call.
In our presentation here we discuss \&{let}-expressions before
$\lambda$-expressions, doing the simpler things first.

@ The parser produces values of type |id_pat| to describe such patterns, and
they are needed at run time to guide the decomposition of the corresponding
values (in fact the names used in the pattern are not needed for this, but they
are useful for printing the \&{let}-expression). The |id_pat| structure can be
used directly in a \&{let}-expression, and will handle ownership properly if we
place it inside the |let_expression| structure (rather than using |raw_id_pat|
as the parser often does, which contains a raw pointer at the top level).
The initialiser expression and body are owned, and held in unique-pointers.

@< Type def... @>=
struct let_expression : public expression_base
{ const id_pat variable;
  expression_ptr initialiser, body;
@)
  let_expression(const id_pat& v, expression_ptr&& ini, expression_ptr&& b);
  virtual ~let_expression() = default; // subobjects do all the work
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ When constructing a |let_expression|, we might place |id_pat| produced by the
parser there using move semantics. Thereby we would steal a possible |sublist|
field and further nodes accessible from it, amputating the expression; that
would produce a weird effect if the expression were then printed in an error
message. So instead we make a deep copy (and patterns are rarely very deep)
using the function~|copy_id_pat|. The implicit recursion of |copy_id_pat| is
achieved by passing (a pointer to) the function itself as final argument to
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

@ The main |let_expression| constructor could not be inside its class
definition, as it uses the local function |copy_id_pat|. The variable (pattern)
is copied into the |let_expression| node, while for the two expression
components ownership is transferred from the passed unique-pointer values.

@< Function def... @>=
inline
let_expression::let_expression @|
  (const id_pat& v, expression_ptr&& ini, expression_ptr&& b)
: variable(copy_id_pat(v))
, initialiser(std::move(ini))
, body(std::move(b))
@+{}

@ To print a \&{let}-expression, we reproduce the input syntax, as far as
possible: any restructuring that was done during parsing cannot be undone.
So we just print a basic form that could have been the input.

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
  (const id_pat& pat,const type_expr& type, unsigned int lvl,
   layer& dst, bool is_const);
void thread_components
  (const id_pat& pat,const shared_value& val,
   std::back_insert_iterator<std::vector<shared_value> > dst);

@ For handling declarations with patterns as left hand side, we need a
corresponding type pattern; for instance $(x,,(f,):z)$:\\{whole} requires the
type pattern \.{(*,*,(*,*))}. The following two mutually recursive functions
construct such type patterns.

@< Function definitions @>=
type_list pattern_list_types(const patlist& p)
{ dressed_type_list result;
  for (auto it = p.begin(); not p.at_end(it); ++it)
    result.push_back(pattern_type(*it));
  return result.undress();
}
@)
type_expr pattern_type(const id_pat& pat)
{@;
  return (pat.kind&0x2)==0
  ? type_expr()
  : type_expr::tuple(pattern_list_types(pat.sublist));
}

@ Here we count the number or list the identifiers in a pattern. The latter
function uses an output parameter rather than a return value since this avoids
doing any concatenation of vectors. Instead of a modifiable reference~|d| to a
vector it could have used an output iterator, but in practice we always
collect the results in a vector, and this avoids having to call
|std::back_inserter| all the time.

@< Function definitions @>=
size_t count_identifiers(const id_pat& pat)
{ size_t result= (pat.kind & 0x1)==0 ? 0 : 1;
                 // or we might assign |pat.kind & 0x1| directly
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

@ Here we do a similar traversal, using a type with structure matching |pat|; we
push a pair onto a |layer| for every identifier in |pat|. All types created will
have the same number |lvl| of fixed type variables, which is passed
unchanged in recursive calls.

If |is_const| holds, all identifiers will be immutable; otherwise any identifier
that is either flagged as such by the user (which is recorded in a bit from
|pat.kind|) or has polymorphic type will be made immutable. The reason that
having polymorphic type implies immutability is explained in the introduction to
polymorphic types in \.{axis-types.w}; essentially, we wish to avoid
ever \emph{requiring} values to have (sufficiently) polymorphic type, as would
be the case when assigning a new value to such an identifier.

@< Function definitions @>=
void thread_bindings
  (const id_pat& pat,const type_expr& te, unsigned int lvl,
   layer& dst, bool is_const)
{ if ((pat.kind & 0x1)!=0)
  {
    type tp = type::wrap(te,lvl);
    unsigned char flags = pat.kind;
    if (tp.is_polymorphic())
      flags |= 0x4; // polymorphic type implies constant
    dst.add(pat.name,std::move(tp),flags);
  }
  if ((pat.kind & 0x2)!=0)
    // recursively traverse sub-list for a tuple of identifiers
  { auto tex = te.expanded(); // ensure substitution into any argument types
    assert(tex.raw_kind()==tuple_type);
    wtl_const_iterator t_it(tex.tuple());
    for (auto p_it=pat.sublist.begin(); not pat.sublist.at_end(p_it);
         ++p_it,++t_it)
      thread_bindings(*p_it,*t_it,lvl,dst,is_const);
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

@ To convert a \&{let}-expression, we first deduce the type of the declared
identifiers from the right hand side of its declaration, then set up new
bindings for those identifiers with the type found, and finally convert the
body to the required type in the extended context. Note that the constructed
|layer| is a local variable whose constructor pushes it onto the
|layer::lexical_context| list, and whose destructor will pop it off.

If there are no identifiers at all, we should avoid that the execution of the
|let_expression| push an empty frame on the evaluation context, as this would
wreak havoc due to the fact that we made applied identifiers not count empty
layers. Rather than have the |let_expression::evaluate| method handle this
dynamically, we avoid generating such a |let_expression| altogether, and instead
generate a |seq_expression| whose |evaluate| method does exactly what is needed
in this case. Note that the most likely scenario where this happens is when the
identifier pattern is $(\,)$ specifying a $0$-tuple (something that happens
quite frequently in actual user code, to effectively insert a statement into a
sequence of \&{let}-declarations), in which case it is particularly important
that we initialised |decl_type| to |pattern_type(pat)| before converting the
right hand side, as the latter will then have |decl_type==void_type|, allowing
expressions of arbitrary type to be accepted due to the voiding coercion. If
instead a specialisation to |pattern_type(pat)| were attempted after the
conversion, this would fail in cases where $(\,)={}$ is followed by an
expression not spontaneously having void type.

@< Cases for type-checking and converting... @>=
case let_expr:
{ const auto& lexp=*e.let_variant;
  const id_pat& pat=lexp.pattern;
  const auto& rhs = lexp.val;
  type decl_type=type::bottom(fc);
  expression_ptr arg = convert_expr(rhs,decl_type);
  type_expr tuple_tp = pattern_type(pat);
  if (not decl_type.unify_specialise(tuple_tp))
  @< Complain that the type |decl_type| of |rhs| does not match the identifier
      pattern |pat| in a let expression @>
@/auto n=count_identifiers(pat);
  if (n==0)
  // then avoid frame without identifiers, so compile as sequence expression
    return expression_ptr(new @|
      seq_expression(std::move(arg),convert_expr(lexp.body,tp)));
  if (decl_type.is_void() and not is_empty(rhs))
    // rare case, introducing void identifier
    arg.reset(new voiding(std::move(arg)));
  layer new_layer(n);
  thread_bindings(pat,tuple_tp,fc,new_layer,false);
  return expression_ptr(new @|
    let_expression(pat,std::move(arg),convert_expr(lexp.body,tp)));
}

@ When reporting, in a \&{let} binding, a mismatch between the right hand side
type and the left hand side pattern, we throw an |expr_error| that mentions just
the |rhs| expression, as the containing \&{let} expression could be much larger.

@< Complain that the type |decl_type| of |rhs| does not match the identifier
   pattern |pat| in a let expression @>=
{ std::ostringstream o;
  o << "Expression has type '" << decl_type
    << "' which is not suited to bind identifiers " << pat;
  throw expr_error(rhs,o.str());
}

@ Here is a class whose main purpose, like that of |layer| before, is to have
a constructor-destructor pair, in this case one that temporarily suspends the
current execution context, replacing it by a new one determined by an
identifier pattern and an execution context for the enclosing lexical layers.
All instances of this class should be automatic (local) variables, to ensure
that they have nested lifetimes. The actual execution context controlled by this
class is a linked structure entirely stored on the heap, with the exception of
the initial |static| pointer |frame::current| that gives access to it. The
|evaluation_context| class is defined in \.{axis-types.w}; each node contains a
|next| link to the next node, and a vector of values for all identifiers
introduced in the pattern corresponding to this node (but that pattern is not
itself recorded in the code).

We take care when pushing a new |evaluation_context| to avoid changing the
reference count of the pointer in |frame::current| as would happen when copying
it, by \emph{moving} the pointer to the tail of the new node. When popping on
destruction however we need to copy, since the node being popped could have
become accessed independently (through a |closure|, to be discussed later), so
we are not free to move from the tail (or any other part) of this node.

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
  { assert(count_identifiers(pattern)>0); // we avoid frames without identifiers
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

@ The method |id_list| is only called during exception handling, so a simple
access to the identifiers is more important than an efficient one. Therefore we
convert the pattern to a vector, following the same order as used in a
corresponding |shared_context| node, using a call to |list_identifiers|.

@< Local function definitions @>=
std::vector<id_type> frame::id_list() const
{ std::vector<id_type> names; names.reserve(count_identifiers(pattern));
  list_identifiers(pattern,names);
  return names;
}

@ Evaluating a \&{let}-expression is now straightforward: evaluate the
initialiser to produce a value on the stack; then create a new |frame| in which
this value is bound to the pattern |variable|, and in this extended context
evaluate the~|body|.

The stack will be automatically popped when the lifetime of the |frame| ends.
Usually this happens when we return from the |let_expression::evaluate|, but
it also happens in case an error is thrown. The latter also includes instances
of \&{break} or \&{return} premature exits from a loop or a user-defined function
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


@ Providing the current values of local variables at an error stop is useful and
quite easy. By calling |fr.id_list| we obtain the names of the identifiers from
the frame |fr| to which they are copied, rather than directly from our
|variable|; the main purpose of this is that it allows reusing the module
identically as |catch| block for user defined functions, where the current
(call) expression does not have an identifier pattern available, but there is a
frame |fr| from which is can be obtained. This will also be reused on multiple
other occasions where new identifiers, such as loop variables, are locally
introduced.

While users cannot leave an undefined value in a local name, certain
optimisations will temporarily move a value from a local variable, so that there
is no sharing that prevents its modification in-place; this happens just before
a new value (quite likely that modified value) gets assigned to the variable. It
can happen however that the function computing the new value throws a runtime
error, and if this happens the code below will find the variable with its pants
down. Therefore we are careful to handle the case where |*it==nullptr|, and
print an indication of absence of any value at this point, rather than crash
the \.{atlas} program.

@< Catch block for providing a trace-back of local variables @>=
catch (error_base& e)
{ std::ostringstream o;
  std::vector<id_type> names = fr.id_list();
  auto id_it = names.cbegin();
  for (auto it = frame::current->begin(); it!=frame::current->end();
       ++it,++id_it)
  { o << (it==frame::current->begin() ? "{ " :", ")
   @| << main_hash_table->name_of(*id_it) << '=';
    if (*it==nullptr)
      o << "(.)"; // indicate a value that is gone
    else
      o << **it;
  }
  o << " }";
  e.trace(o.str());
  throw;
}


@*1 Lambda-expressions (user-defined functions).
%
When a user defines a function, either globally with \&{set} or inside an
expression (possibly another user defined function) using \&{let}, this creates,
just like for a variable definition, a binding between a name and a value, the
two being (as usual) treated as separate entities. The value created is a called
a $\lambda$-expression (in honour of the $\lambda$ calculus; however using the
Greek letter $\lambda$ would be a dreadful choice of notation, which we avoid)
and it represents the function, with its argument pattern and body (the
definiens), but without the name given (the definiendum). In fact one can write
a $\lambda$-expression without giving it a name at all, creating an anonymous
user defined function; this can be quite practical, for instance to pass the
function to another function. A $\lambda$-expression is basically a denotation
for a piece of code, whose evaluation involves no action and returns that code.
However a $\lambda$-expression may refer to identifiers known in the context in
which it is written, so evaluating it captures the bindings of those identifiers
at the time of evaluation into a runtime value know as a closure.

So in contrast to \&{let}-expressions which just extend the context with new
bindings, a closure captures the current context inside a value, which after
possibly being passed around, can later be combined with argument values to
provide an extended context in which its body is evaluated. A closure might
outlive the expression that contained its $\lambda$-expression, and evaluations
of the same textual $\lambda$-expression can create multiple closures whose
lifetimes may overlap.

We wish to create closures without copying the $\lambda$-expression into them,
rather by storing a reference, yet memory for the $\lambda$-expression should be
freed when the last reference to it disappears, so we must use a
|std::shared_ptr| based mechanism for sharing the $\lambda$-expression among its
closures. Now the $\lambda$-expression is already referred to by an
|expression_ptr| smart pointer from its containing expression, but that is a
unique pointer that cannot be used for sharing among closures. So instead we
define |shared_lambda| as shared pointer to a |lambda_struct| that contains the
data of a $\lambda$-expression, and |lambda_expression| as a class derived from
|expression_base| that contains a |shared_lambda| as unique data member.

We thus make a distinction between the $\lambda$-expression object (class
|lambda_expression|) that is executable, yielding a closure upon evaluation, and
the |lambda_struct| object that it shares with its closures. The former is of
type |lambda_expression|, accessed through unique pointers of type
|expression_ptr|, while |lambda_struct| will be accessed via shared pointers of
type |shared_lambda|. Upon construction of a |lambda_expression|, it dynamically
creates a subsidiary object of class |lambda_struct| where it stores all its
information (argument pattern, body, source location), and which will be shared
with the closures it spawns.

@< Type def... @>=
struct lambda_struct
{ id_pat param; @+ expression_ptr body; @+ source_location loc;
@)
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
    (id_pat&& pat, expression_ptr&& b, const source_location& loc)
    : p(std::make_shared<lambda_struct>(std::move(pat),std::move(b),loc)) @+{}

  virtual ~lambda_expression() = default;
    // subobjects do all the work
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ In a function definition the $\lambda$-expression itself does not have access
to the name the user binds it to (it is included in the context
only \emph{after} the definition is processed), so an ordinary
$\lambda$-expression cannot call itself recursively. In order give explicit
support for recursive functions, we define a second structure |recursive_lambda|
for a $\lambda$-expression whose body can also refer to the $\lambda$-expression
itself. In addition to a |lambda_expression| it will store a recursive
identifier the user gave to it (called |self_id| below). In spite of this
addition, we derive |recursive_lambda| from |lambda_expression| without adding
any data members: |self_id| will be combined into the identifier pattern already
used for the ordinary argument names in |lambda_expression|. The distinctive
interpretation of this pattern will be handled by the virtual method
|recursive_lambda::evaluate|, and by similar methods of the recursive closures
it generates.

@< Type def... @>=
struct recursive_lambda : public lambda_expression
{ recursive_lambda@|(id_type self_id,
     id_pat&& p, expression_ptr&& b, const source_location& loc);
  virtual ~recursive_lambda() = default;
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ The |recursive_lambda| constructor forms as |id_pat| a pair with the recursive
identifier as first component and the usual argument pattern as second field,
and passes that to the |lambda_expression| constructor, which in fact places it
into its |lambda_struct| object, where it can be seen by closures created by the
|recursive_lambda|. Although the recursive identifier cannot be assigned to,
there is no need to set the constant bit in its pattern here, as that property
will be enforced when setting up the context for type checking the body of a
recursive function.

@< Function def... @>=
id_pat rec_pair(id_type s,id_pat && p)
{@;
  patlist pl;
  pl.push_front(std::move(p));
  pl.push_front(id_pat(s));
  return id_pat(std::move(pl));
}
@)
inline recursive_lambda::recursive_lambda@|(id_type self_id,
     id_pat&& p, expression_ptr&& b, const source_location& loc)
  : lambda_expression(rec_pair(self_id,std::move(p)),std::move(b),loc) @+{}


@ To print an anonymous function, we print the parameter list, followed by a
colon and by the function body. If the parameter list contains a name for the
whole, as happens in particular when there is just a single parameter, then it
must be enclosed in parentheses to resemble in the input syntax, but if it is an
unnamed nonempty tuple then it will supply its own parentheses; this leaves the
case of an empty parameter list, for which we reconstitute the input
syntax~`\.@@'. The printed parameter list cannot include types with the current
implementation, as they are not explicitly stored after type analysis. It could
be made possible to print types if a type were explicitly stored in the
|lambda_expression| structure; at the time of writing this would seem possible
because each function has to have a definite type, but if the type system were
extended with second order types (which would be quite useful), then this might
no longer be true. For now leaving out the types indicates to the
(knowledgeable) user that a runtime value is being printed rather than just a
syntax tree representing the user input (as happens in messages from the type
checker).

@< Function definitions @>=
void print_lambda(std::ostream& out,
  const id_pat& param, const expression_ptr& body)
{ if ((param.kind&0x1)!=0) // single argument
    out << '(' << param << ')';
  else if ((param.kind&0x2)!=0 and not param.sublist.empty())
    out << param;
  else out << '@@';
  out << ": " << *body;
}
@)
void lambda_expression::print(std::ostream& out) const
{@; print_lambda(out,p->param,p->body); }
@)
void recursive_lambda::print(std::ostream& out) const
{ auto it=p->param.sublist.begin();
  id_type self_id=it->name;
  const id_pat& param=*++it;
  out << "(recfun " << main_hash_table->name_of(self_id) << ' ';
  print_lambda(out,param,p->body);
  out << ')';
}

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

In non-void context we specialise the required |tp| (often undetermined
initially) to a function type with argument type the one given in the
$\lambda$-expression (signalling a type error if a different type was
expected), then we convert the function body in the new context, specialising
the return type. In void context we do only the conversion (just for error
checking), and ignore the return type.

@< Cases for type-checking and converting... @>=
case lambda_expr:
{ const lambda_node& fun=*e.lambda_variant;
  const id_pat& pat=fun.pattern;
  type_expr fun_type = gen_func_type.copy();
  type_expr& par_tp = fun_type.func()->arg_type;
  par_tp.specialise(global_id_table->expand(fun.parameter_type));
      // argument type specified in |fun|
  if (not par_tp.can_specialise(pattern_type(pat)))
    throw expr_error
      (e,"Specified parameter type does not match identifier pattern");
  if (not(tp.unify_specialise(fun_type) or tp.is_void()))
    throw type_error(e,std::move(fun_type), tp.bake());
@)
  type body_type = type::wrap(fun_type.func()->result_type,fc);
    // probably undetermined here
@/layer new_layer(count_identifiers(pat),&body_type);
    // create layer for function body
  thread_bindings(pat,par_tp,fc,new_layer,false);
@/expression_ptr body_ptr = convert_expr(fun.body,body_type);
  if (not tp.is_void()) // then we must integrate |body_type| into |tp|
  { const type_expr& sub_tp = tp.top_expr().func()->result_type;
    if (not tp.unify_to(sub_tp,body_type))
      throw logic_error("Failed to export result type of lambda-expression");
  }
@/return expression_ptr(new @| lambda_expression
   (copy_id_pat(pat), std::move(body_ptr), std::move(e.loc)));
}

@ Type-checking recursive lambda expressions is slightly different. We must, in
addition to binding the argument pattern to the argument type, bind the
recursive identifier to the function type given by |fun.parameter_type| and
|fun.result_type|; we build that type as |f_type| rather than specialise and
then use |tp|, just to allow cases with |type==void_type| initially, which
even if silly is legal. The type |fun.return_type| will be determined (due to
syntactic constraints on recursive functions) so |f_type| will be full
determined, and |tp| will be specialised to it during the condition
|tp.specialise(f_type)|. The call to |thread_bindings| for the recursive
identifier has final argument |true| indicating the identifier should be treated
as constant during the call to |convert_expr| for the body. That call needs an
lvalue for its expected type, and the |layer| constructor needs the same in the
form of a pointer to non-|const| (in both cases because the type might be
specialised, from the type of the body respectively from |return| clauses inside
it), For these uses we supply |f_type.func()->result_type| (aliased |r_type|)
rather than using |tp.func_type()->return_type|, which since |f_type| is a
local variable means any specialisation done to it would be lost; however there
should be no such changes since |fun.result_type| is already fully specialised.

@< Cases for type-checking and converting... @>=
case rec_lambda_expr:
{ const rec_lambda_node& fun=*e.rec_lambda_variant;
  const id_pat& pat=fun.pattern;
  type_expr fun_type = type_expr::function
    ( global_id_table->expand(fun.parameter_type)
    , global_id_table->expand(fun.result_type)
    );
  type_expr& par_tp = fun_type.func()->arg_type;
  if (not par_tp.can_specialise(pattern_type(pat)))
    throw expr_error
      (e,"Specified parameter type does not match identifier pattern");
  if (not (tp.unify_specialise(fun_type) or tp.is_void()))
      @/throw type_error(e, std::move(fun_type), tp.bake());
@)
  type body_type = type::wrap(fun_type.func()->result_type,fc);
  layer new_layer(1+count_identifiers(pat),&body_type);
@/thread_bindings(id_pat(fun.self_id),fun_type,fc,new_layer,true);
  thread_bindings(pat,par_tp,fc,new_layer,false);
@/expression_ptr body_ptr = convert_expr(fun.body,body_type);
  if (not tp.is_void()) // then we must integrate |body_type| into |tp|
  { const type_expr& sub_tp = tp.top_expr().func()->result_type;
    if (not tp.unify_to(sub_tp,body_type))
      throw logic_error("Failed to export result type of lambda-expression");
  }
@/return expression_ptr(new @| recursive_lambda
   (fun.self_id,copy_id_pat(pat), std::move(body_ptr), std::move(e.loc)));
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
dereference upon evaluation. Since the latter occurs frequently at run time, we
speed up evaluation by also using a reference |body| directly to the function
body. Note that closures are formed when the \emph{definition} of a user-defined
function is processed, so this optimisation should make evaluation of globally
defined functions a bit faster. For local functions, the closure is formed
during the execution of the outer function, so the optimisation only helps if
the closure formed will be called more than once; this is still probable, though
there are usage patterns (for instance finding the first of a sequence of
conditions that is satisfied, the conditions being produced as a list of
parameterless functions in a loop) for which local closures are actually
executed less than once on average; in such cases we are actually wasting effort
here. It is however impossible to know here whether the closure will be globally
of locally bound.

Because of slight differences in evaluation that we do not want to implement
using runtime tests, we define three types of closure. The case of non recursive
$\lambda$-expressions being split into those that introduce no parameters at all
(usually but not necessarily because the argument type is |void|), and those
that introduce at least one name; this distinction reflects our implementation
choice to omit empty layers on the stack of local bindings. For recursive
$\lambda$-expressions there is (only) a third type of closures; calling them
will always push at least the recursive identifier.

We defined a virtual method |function_base::maybe_push| that will push a shared
pointer (which in fact will be one to the closure itself) on the stack, as is
needed for the implementation of recursive functions. Since the default
implementation of this method is to do nothing, we need to implement it only for
|kind==recursive_closure|, but it is most convenient to define it here
regardless of |kind|, and in its implementation use a test of the
template argument, which the compiler should hopefully optimise away.

@< Type def... @>=

enum Closure_kind @+{ parameterless, with_parameters, recursive_closure };
@)
template<Closure_kind kind>
struct closure_value : public function_base
{ shared_context context;
  shared_lambda p;
  const expression_base& body; // shortcut to function body
@)
  closure_value@|(const shared_context& c, const shared_lambda& l)
  : function_base(), context(c), p(l), body(*p->body) @+{}
  virtual ~closure_value() = default;
  virtual void print(std::ostream& out) const;
  virtual void apply(eval_level l) const;
  virtual eval_level argument_policy() const
  {@; return eval_level::single_value; }
  virtual void maybe_push(const std::shared_ptr<const function_base>& p) const
  {@; if (kind==recursive_closure)
      push_value(p);
  }
  virtual void report_origin(std::ostream& o) const;
  virtual expression_ptr build_call
    (const shared_function& master,const std::string& name,
     expression_ptr&& arg, const source_location& loc) const;
@)
  static const char* name() @+
   {@; return kind==recursive_closure ? "recursive closure" : "closure"; }
  closure_value (const closure_value& ) = delete;
};

@ For readability, we define |shared_closure| as a type template.

@s shared_closure vector
@< Type def... @>=
template<Closure_kind kind>
using shared_closure = std::shared_ptr<const closure_value<kind> >;

@ A closure prints the |lambda_expression| from which it was obtained, but we
also print an indication of where the function was defined (this was not
useful for |lambda_expression|, since these never get printed directly, only
as part of printing a |closure_value|). One could imagine printing after this
body ``where'' followed by the bindings held in the |context| field. Even
better only the bindings for relevant (because referenced) identifiers could
be printed. But it's not done yet.

@< Function def... @>=
template<>
void closure_value<parameterless>::print(std::ostream& out) const
{@; print_lambda(out << "Function defined " << p->loc << std::endl,
                 p->param,p->body);
}
template<>
void closure_value<with_parameters>::print(std::ostream& out) const
{@; print_lambda(out << "Function defined " << p->loc << std::endl,
                 p->param,p->body);
}
template<>
void closure_value<recursive_closure>::print(std::ostream& out) const
{ auto it=p->param.sublist.begin();
  id_type self_id=it->name;
  const id_pat& param=*++it;
  out << "Recursive function defined " << p->loc << std::endl;
  print_lambda(out<< main_hash_table->name_of(self_id) << " = ",param,p->body);
}

template<Closure_kind kind>
void closure_value<kind>::report_origin(std::ostream& o) const
{@; o << "defined " << p->loc; }

@ Evaluating a $\lambda$-expression just forms a closure using the current
execution context, and returns that. Since the |recursive_lambda| already has
done the work of wrapping the recursive identifier into a pattern with that of
the arguments, the recursive case only differs from the ordinary one here by the
setting of the template argument of the |closure_value| constructor. For non
recursive $\lambda$-expression, this is the place where we test for the absence
of any introduced identifiers. The actual difference in implementation between
the three types of closure will be mostly evident in the
|closure_value<kind>::apply| methods.

While this code looks rather innocent, note that the sharing of
|frame::current| created here may survive after one or more frames on the list
|frame::current| get removed after evaluation of the expression that returned
the closure; these frames then get an extended lifetime through the closure
formed here. This implies that the execution context cannot be embedded in any
kind of stack, in particular it cannot be embedded in the \Cpp\ runtime stack
(while the layers of the lexical context could).

@< Function def... @>=
void lambda_expression::evaluate(level l) const
{ if (l!=level::no_value)
  { if (count_identifiers(p->param)==0)
      push_value(std::make_shared<closure_value<parameterless> >
       (frame::current,p));
    else
      push_value(std::make_shared<closure_value<with_parameters> >
       (frame::current,p));
  }
}
void recursive_lambda::evaluate(level l) const
{ if (l!=level::no_value)
    push_value(std::make_shared<closure_value<recursive_closure> >
      (frame::current,p));
}

@*1 Calling user-defined functions.
%
In order to implement calling of user-defined functions, we define a variation
of the class |frame|. Again the purpose is to have a constructor-destructor pair
that temporarily suspends the current execution context, replacing it by a new
one determined by the parameter(s) of the $\lambda$-expression, on top of the
execution context stored in the closure. Here too, instances of this class
should be automatic variables, to ensure that they have nested lifetimes. More
precisely, |lambda_frame| variables are local to functions performing evaluation
of user defined functions, and their construction modifies |frame::current| by
swapping it out with the context of the closure, which means replacement of the
run time stack for local user variables by another that is partially or wholly
disjoint from the previous stack. The new stack remains in effect until the end
of the function containing the variable (when the destructor reinstates the old
stack). The stack records themselves are in dynamic memory, not on the \Cpp\
runtime stack, and will remain in existence as long as anybody might access
them.

This context switching is a crucial and recurrent step in the evaluation
process, so we take care to not uselessly change the reference count of
|frame::current|. It is \emph{moved} into |saved| upon construction, and upon
destruction moved back again to |frame::current|. We also use a template
parameter to indicate whether we are called for a |parameterless| closure or
not, which avoids any run time test for this condition. The distinction between
recursive and non recursive closures is not relevant here.

It might seem that we could have derived this class from |frame|, which already
provides a similar constructor-destructor pair. However, we would then have to
construct the base (which modifies the value |frame::current|) before doing
anything else; this would make saving the value of |frame::current| problematic.
For that reason it is better to keep the two independent, and just repeat some
of the things done for |frame|, similarly but with a few important changes. We
do of course have to make use of the static member |frame::current| of that
class, which is the whole point of defining the |lambda_frame| class.

@< Local class definitions @>=
template<bool no_names>
class lambda_frame
{
  const id_pat& pattern;
  const shared_context saved;
public:
  lambda_frame (const id_pat& pattern, const shared_context& outer);
  ~lambda_frame() @+{@; frame::current = std::move(saved); }
@)
  void bind (const shared_value& val);
  std::vector<id_type> id_list() const; // list identifiers, for back-tracing
};

@ In contrast to |frame|, the constructor here needs a try block for exception
safety, as an exception may be thrown during construction (in the call to
|std::make_shared<evaluation_context>|), after |frame::current| has been moved
from, but before out constructor completes; since the destructor would in this
scenario \emph{not} be called, we then need to move the pointer back explicitly
in the |catch| block.

This function requires that the argument |outer| is not an alias of
|frame_current|, as its logic would then fail. This condition is satisfied
whenever |outer| is a |context| member of a closure, as it should. When
|no_names| holds, the addition of a new |evaluation_context| stack frame (for
the function parameters) is omitted here.

@< Local function definitions @>=
template<bool no_names>
  lambda_frame<no_names>::lambda_frame
    (const id_pat& pattern, const shared_context& outer)
  : pattern(pattern)
  , saved(std::move(frame::current))
  { assert(&outer!=&frame::current); // for excluded case, use |frame| instead
    try {@;
      frame::current =
        no_names ? outer : std::make_shared<evaluation_context>(outer);
    }
    catch(...)
    {@; frame::current = std::move(saved); throw; }
    // restore as destructor would do
  }

@ Since a |lambda_frame::bind| will be called every time a user defined function
with arguments is called, it provides a convenient point to check whether the
signal handler has set the interrupt flag, and to bail out if it did. (Choosing
the points to do this test is a somewhat delicate matter. One wants to do it at
points that are regularly encountered during evaluation, but not so frequently
that the checks incur a serious performance penalty. Just having the test here
does not provide an absolute guarantee of rapid detection of a signalled
interrupt.) This call to |check_interrupt| is the reason this method definition
was lifted out of the class definition. Since this method is never called when
the template argument |no_name| holds, we take the opportunity to only define it
for the other case (in practice it makes no difference to define a declared
method that is never called, but we want to avoid the impression that
|check_interrupt| is called for parameterless functions).

@< Local function definitions @>=
template<>
void lambda_frame<false>::bind (const shared_value& val)
{ check_interrupt();
  frame::current->reserve(count_identifiers(pattern));
  thread_components(pattern,val,frame::current->back_inserter());
}

@ This method is identical to the one in |frame|, but as said, we cannot use
inheritance. Again, it is never called if |no_names| holds.

@< Local function definitions @>=
template<>
std::vector<id_type> lambda_frame<false>::id_list () const
{ std::vector<id_type> names; names.reserve(count_identifiers(pattern));
  list_identifiers(pattern,names);
  return names;
}

@ In general a closure formed from a $\lambda$ expression can be handled in
various ways (like being passed as argument, returned, stored) before being
applied as a function, in which case the call is performed by
|call_expression::evaluate| described above in
%
section@# general call evaluate @>; the actual code this executes (through the
virtual method |apply|) is given in section@# lambda evaluation @> below.
However, in many cases the path from definition to call is more direct: the
closure from a user-defined function is bound to an identifier (or operator) in
the global overload table, and located during type-checking of a call
expression. As this special but frequent case can be handled more efficiently
than by building a |call_expression|, we introduce a new |expression| type
|closure_call| that is capable of directly storing a closure value, evaluating
only its argument sub-expression at run time.

Closures themselves are anonymous, so the name |n| that becomes the |name| field
of the |overloaded_call| base object reflects the overloaded name that was used
to identify this function; it can vary separately from the closure |f| if the
latter is entered more than once in the the tables. This is in contrast to
|builtin_call| where the name is taken from the stored |builtin_value|, and
cannot be dissociated from the wrapper function.

@< Type definitions @>=
template<Closure_kind kind>
class closure_call : public overloaded_call
{ shared_closure<kind> f; // remaining fields are shortcuts into |f|
  const id_pat& param;
  const shared_context& context;
  const expression_base& body;
public:
  closure_call @|
   (shared_closure<kind >&& f_ref,const std::string& n,expression_ptr&& a
   ,const source_location& loc)
@/: overloaded_call(n,std::move(a),loc)
  , f(std::move(f_ref))
@/, param(f->p->param)
  , context(f->context)
  , body(f->body) @+ {}
  virtual ~closure_call() = default;
  virtual void evaluate(level l) const;
};

@ Here is how a |closure_value| can turn itself into a |closure_call| when
provided with an argument expression, as well as a |name| to call itself and a
|source_location| for the call. Every call of |build_call| will provide the
shared pointer to the |functions_base| derived object it is called for as first
argument |master|; for lack of covariance this pointer has been up-cast to
|shared_function|. In order to provide the constructed |closure_call| with a
shared pointer |f| to our |closure_value|, we perform a (static) down-cast of
|master|. The provided argument(s) |arg|, and |name| are also passed into the
|closure_call|.

@< Function def... @>=
template<Closure_kind kind>
expression_ptr closure_value<kind>::build_call
    (const shared_function& master,const std::string& name,
     expression_ptr&& arg, const source_location& loc) const
{ assert(master.get()==this);
  auto p=std::static_pointer_cast<const closure_value<kind> >(master);
@/return expression_ptr(new @| closure_call<kind>
    (std::move(p),name,std::move(arg),loc));
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

When we come here, the argument of our |closure_value| has already been
evaluated, and is available as a single value on the |execution_stack|. The
evaluation of the call temporarily replaces the current execution context
|frame::current| by one composed of |f->context| stored in the closure and a new
frame defined by the parameter list |f->param| and the argument obtained as
|pop_value()|; the function body is evaluated in this extended context.
Afterwards the original context is restored by the destructor of~|fr|, whether
the call completes normally or is terminated by a runtime error. The most
important advantage of this approach is in case of abnormal exits from loops and
functions, which are implemented by throwing and catching an exception at
run-time and will therefore unwind the \Cpp\ stack. (Actually, performing
\&{break} from a loop should never lead to destructing any |lambda_frame|,
though it might destruct some |frame|s.)

By naming our frame |fr|, and because |lambda_frame| has a method |id_list|
with the same signature as |frame::id_list|, we can textually reuse a |catch|
block, as was mentioned when that block was defined earlier.

@: lambda evaluation @>

@< Function def... @>=
template<Closure_kind kind>
void closure_value<kind>::apply(eval_level l) const
{
  static const bool no_names=kind==parameterless;
  try
  { lambda_frame<no_names> fr(p->param,context);
      // save context, create new one for |*this|
    if (no_names) // functions without named arguments are different
    {@; execution_stack.pop_back();
      body.evaluate(l);
    } //  drop arg, evaluate avoiding |bind|
    else
    { fr.bind(pop_value()); // decompose arguments(s) and bind values in |fr|
      try {@; body.evaluate(l); }
        // call, passing evaluation level |l| to function body
      @< Catch block for providing a trace-back of local variables @>
    }
  } // restore context upon destruction of |fr|
  @< Catch block for explicit \&{return} from functions @>
}

@ Here is where recursive functions differ from non-recursive ones. Instead of
just binding the arguments popped from the stack, we build a pair consisting of
our current (recursive) |closure_value| itself and the argument, and bind that
to the pattern of the |lambda_struct| access from the closure. Since the
recursive name is always present, we do not have to worry about the possibility
of a |lambda_frame| without any identifiers.

@< Function def... @>=
template<>
void closure_value<recursive_closure>::apply(eval_level l) const
{
  try
  { lambda_frame<false> fr(p->param,context);
      // save context, create new one for |*this|
    wrap_tuple<2>(); // combine pre-pushed self object with pushed argument(s)
    fr.bind(pop_value()); // bind self value and arguments in |fr|
    try {@; body.evaluate(l); }
      // call, passing evaluation level |l| to function body
      @< Catch block for providing a trace-back of local variables @>
  } // restore context upon destruction of |fr|
  @< Catch block for explicit \&{return} from functions @>
}

@ When the call |f.body.evaluate(l)| ends up executing an explicit \&{return}
expression, it won't have put anything on the |execution_stack|, but the value
to return will be stored inside the error object. (This is the right thing to
do in the unlikely case that intermediate values are removed from
|execution_stack| during the \Cpp\ stack unwinding.) It suffices to place the
value on the stack, which is just what |push_expanded| does.

@< Catch block for explicit \&{return} from functions @>=
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

Not having varied our naming conventions (|f| and |arg_string|), we can make
a fourth textual reuse of the |catch| block for function calls, as well as
(|fr|) a third reuse of the |catch| block for local variables.

@< Function definitions @>=
template<Closure_kind kind>
void closure_call<kind>::evaluate(level l) const
{ static const bool no_names=kind==parameterless;
  if (kind==recursive_closure)
    push_value(f); // duplicate our closure to the execution stack
  argument->eval(); // evaluate arguments as a single value
  std::string arg_string;
  if (verbosity!=0) // then record argument(s) as string
  {@; std::ostringstream o;
    o << *execution_stack.back();
    arg_string = o.str();
  }
  try
  { lambda_frame<no_names> fr(param,context);
      // save context, create new one for |f|
@)
    if (no_names)
       // we must test for functions without named arguments
      {@; execution_stack.pop_back(); body.evaluate(l); }
      //  avoid |bind|, evaluate
    else
    { if (kind==recursive_closure)
        wrap_tuple<2>(); // combine our closure with actual arguments
      fr.bind(pop_value()); // decompose arguments(s) and bind values in |fr|
      try {@; body.evaluate(l); }
      // call, passing evaluation level |l| to function body
      @< Catch block for providing a trace-back of local variables @>
    }
  } // restore context upon destruction of |fr|
  @< Catch block for explicit \&{return} from functions @>
  @< Catch block for exceptions thrown within call of |f|... @>
}

@* Type abstraction expressions.
%
The extension of the axis language with second order types was made remarkably
late, in view of the fact the execution mechanism always has been able to handle
expressions that can evaluate to a value of unknown type. The main difficulty
was adapting the originally monomorphic type system, and related code like
resolving overloaded function calls, to allowing polymorphic types throughout.

Once this is done, a number of built-in polymorphic functions (like those doing
row concatenation) that used to be given special treatment during type checking,
can be handled like any other functions. Similarly the empty list display and
expressions that never return like calls of |error| naturally get a polymorphic
type. But a mechanism is still needed to allow users to define polymorphic
functions. In many (pure or not quite pure) functional languages
like \.{Haskell} and \.{OCaml} this is done by in general \emph{inferring} types
of identifiers and expressions from the way they are used in an equation-solving
style. But this method is incompatible with our style of overloading in function
calls, since the same symbol can represent functions of many unrelated types;
this makes it impossible to associate an equation (which would need to be in
terms of unknown argument types) with such calls. Instead we rely on identifiers
receiving a type directly at their definition, which in particular requires
explicitly declaring types for function parameters.

To allow defining polymorphic functions, we then use the simple mechanism of
allowing local type variables to be introduced whose scope is a given
expression. In that scope, these type variables designate types that are treated
as primitive types, but without any names or values that relate to them. This
makes it impossible to in any way ``look inside'' values of such types (not even
to test equality), but they can be moved around and stored in data structures
and so on. If the expression that forms the scope of the type variables is
properly type checked, its type may involve those type variables; when
interpreted in its context, the type variable changes status to become
universally quantified in a polymorphic type.

The expression type that introduces type variables is |type_abstraction_expr|.
The only information it receives about the type variables is their number
|a.count|, since the lexical analyser has been set up to recognise the
introduced identifiers within their scope as type variables, and to number them
successively from $0$ upwards. All that remains to be done in type checking is
then to raise the level |fc| of fixed type variables by this number when going
in. When coming back out, the type variables are implicitly ``unfixed'',
becoming free variables of a polymorphic type, which is exported to the context
type~|tp| as usual by calling |tp.specialise|. It should be  noted that |tp| was
not passed inwards in the call of |convert_expr|, which instead passes an
initially undetermined type |result_type|. This reflects our decision to not
allow context to \emph{require} polymorphic types, so nothing would be gained
from passing |tp| in; also we don't have to worry about any change of
interpretation of type variables when going in.

@< Cases for type-checking and converting... @>=
case type_abstraction_expr:
{ const abstr_node& a = *e.abstr_variant;
  tp.assign().set_floor(fc+a.count);
    // increase abstraction level in the scope |a.exp|
  expression_ptr result  = convert_expr(a.exp,tp);
  tp.assign().set_floor(fc); // restore outer abstraction level
  return result;
}


@* Sequence expressions.
%
Sequence expressions are used to evaluate two expressions one after the other,
discarding any value from one of them (although it is possible that the other
value also gets discarded by voiding). In most cases it is the first
expression whose value is discarded (as for the comma-operator in \Cee/\Cpp),
and the semicolon is used to indicate this; occasionally however it is useful
to retain the first value and evaluate a second expression for its side
effects \emph{afterwards}, and the \.{next} keyword is used instead of a
semicolon for this purpose.

Long sequences of expressions chained by semicolons are rare, since often the
need to introduce new local variables interrupts the chain. Therefore we see no
need for a single expression representing an arbitrarily long sequence (using a
|std::vector| internally), and prefer to use an expression node with just two
descendent expressions. The semicolon and \&{next} variants are implemented by
similar but distinct types derived from |expression_base|.

@< Type def... @>=
struct seq_expression : public expression_base
{ expression_ptr first,last;
@)
  seq_expression(expression_ptr&& f,expression_ptr&& l)
   : first(f.release()),last(l.release()) @+{}
  virtual ~seq_expression() = default;
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};
@)
struct next_expression : public expression_base
{ expression_ptr first,last;
@)
  next_expression(expression_ptr&& f,expression_ptr&& l)
   : first(f.release()),last(l.release()) @+{}
  virtual ~next_expression() = default;
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
  expression_ptr first = convert_expr_strongly(seq.first,fc,void_type);
  expression_ptr last  = convert_expr(seq.last,tp);
  return expression_ptr(new seq_expression(std::move(first),std::move(last)));
}
break;
case next_expr:
{ const sequence_node& seq=*e.sequence_variant;
  type void_tp = type::wrap(void_type,fc);
  expression_ptr first = convert_expr(seq.first,tp);
  expression_ptr last  = convert_expr_strongly(seq.last,fc,void_type);
  return expression_ptr(new next_expression(std::move(first),std::move(last)));
}
break;

@* Array subscription and slicing.
%
We have seen expressions to build lists, and although vectors and matrices can
be made out of them using coercions, we so far are not able to access their
components once they are constructed. To that end we shall now introduce
operations to index such values. We allow subscription of rows, but also of
vectors, rational vectors, matrices, strings, and of value of the
Atlas \.{KTypePol} and \.{ParamPol} types. Since after type analysis we know
which of the cases applies for a given expression, we define several classes
among which type analysis will choose. These classes differ mostly by their
|evaluate| method, so we first derive an intermediate class from
|expression_base|, and derive the others from it. This class also serves to host
an enumeration type and some static methods that will serve later. We include
two cases here, |K_type_poly_term| and |mod_poly_term|, that are related to
types defined in \.{atlas-types}.

At the same time we shall define ``slices'', which differ from subscriptions in
selecting a whole range of index values, producing a value of the same type as
the original, though with fewer elements. The name is not really appropriate, as
it more evokes subscripting at certain index positions while leaving other
positions vary (an example would be selecting a row from a matrix); we don't
(yet?) support that as a form of slicing, though selecting an entire column from
a matrix is possible as a special form of subscripting.

@< Type definitions @>=
struct subscr_base : public expression_base
{ enum sub_type
  { row_entry, vector_entry, ratvec_entry, string_char
  , matrix_entry, matrix_column, K_type_poly_term, mod_poly_term, not_so };
  expression_ptr array, index; // the two parts of the subscription expression
@)
  subscr_base(expression_ptr&& a, expression_ptr&& i)
@/: array(a.release())
  , index(i.release())
  @+{}
  virtual ~subscr_base() = default;
@)
  virtual void print(std::ostream&out) const
  // used only for cases without reversal option
  {@; out << *array << '[' << *index << ']'; }
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
  slice_base(expression_ptr&& a, expression_ptr&& lwb, expression_ptr&& upb)
@/: array(a.release())
  , lower(lwb.release())
  , upper(upb.release())
  @+{}
  virtual ~slice_base() = default;
@)
  virtual void print(std::ostream&out) const@+{} // never used; avoid warning
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

@ The final two classes derived from |subscr_base| are slightly different, in
that the are not templated over |reversed| (since these are indexing associative
arrays, there is no notion of reverse indexing). As a consequence they need not
redefine the virtual |subscr_base::print| method.

@< Type definitions @>=

struct K_type_pol_coefficient : public subscr_base
{ K_type_pol_coefficient(expression_ptr&& pol, expression_ptr&& K_type)
@/: subscr_base(std::move(pol),std::move(K_type)) @+{}
  virtual void evaluate(level l) const;
};
@)
struct module_coefficient : public subscr_base
{ module_coefficient(expression_ptr&& pol, expression_ptr&& param)
@/: subscr_base(std::move(pol),std::move(param)) @+{}
  virtual void evaluate(level l) const;
};


@ We derive a number of types, templated over |unsigned|, from |slice_base|.
The template value represents the optional reversals ($8$ possibilities). All
of these differ only by their |evaluate| method.

@< Type definitions @>=

template <unsigned flags>
struct row_slice : public slice_base
{ row_slice(expression_ptr&& a, expression_ptr&& lwb,  expression_ptr&& upb)
@/: slice_base(std::move(a),std::move(lwb),std::move(upb)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; slice_base::print(out,flags); }
};

@)
template <unsigned flags>
struct vector_slice : public slice_base
{ vector_slice(expression_ptr&& a, expression_ptr&& lwb,  expression_ptr&& upb)
@/: slice_base(std::move(a),std::move(lwb),std::move(upb)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; slice_base::print(out,flags); }
};
@)
template <unsigned flags>
struct ratvec_slice : public slice_base
{ ratvec_slice(expression_ptr&& a, expression_ptr&& lwb,  expression_ptr&& upb)
@/: slice_base(std::move(a),std::move(lwb),std::move(upb)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; slice_base::print(out,flags); }
};
@)
template <unsigned flags>
struct string_slice : public slice_base
{ string_slice(expression_ptr&& a, expression_ptr&& lwb,  expression_ptr&& upb)
@/: slice_base(std::move(a),std::move(lwb),std::move(upb)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; slice_base::print(out,flags); }
};
@)
template <unsigned flags>
struct matrix_slice : public slice_base
{ matrix_slice(expression_ptr&& a, expression_ptr&& lwb,  expression_ptr&& upb)
@/: slice_base(std::move(a),std::move(lwb),std::move(upb)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const
  @+{@; slice_base::print(out,flags); }
};


@ These subscriptions and slices are printed in the usual source syntax, with a
tilde after the array name if it is conceptually reversed before the selection
or slice, and in case of slices possibly after each bound, in case that bound
counts displacement from the rear rather than from the front end. The templated
virtual |print| methods transform their template argument to a runtime value,
and call the non |virtual| member |slice_base::print| with that value as second
argument. For matrix subscriptions, where the index type is \.{(int,int)}, the
index expression is quite likely to be a tuple display; therefore we do some
effort to suppress parentheses for that case, in accordance with allowed source
syntax. Since we have passed the type check here, we know that any tuple display
in the index position is necessarily a pair.

@< Function definitions @>=
void subscr_base::print(std::ostream& out,bool reversed) const
{@; out << *array << (reversed ? "~[" : "[") << *index << ']';
}
@)
template <bool reversed>
void matrix_subscription<reversed>::print(std::ostream& out) const
{ out  << *array << (reversed ? "~[" : "[");
  auto* p=dynamic_cast<const tuple_expression*>(index.get());
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
combinations. Upon success, the last two parameters serve to store the type the
subscription will result in, and an element of the |sub_type| enumeration that
indicates the kind of subscription that was found. The |K_type_poly_term| case
indicates that a ``$K$-type polynomial'' can be subscripted with a $K$-type to
return a split integer result, and the |mod_poly_term| case similarly indicates
that a ``parameter polynomial'' can be subscripted with a parameter to return a
split integer result.

@< Function def... @>=
subscr_base::sub_type subscr_base::index_kind
  (const type_expr& aggr,
   const type_expr& index,
         type_expr& subscr)
{ type_expr ag_tp = aggr.expanded();
  if (ag_tp.raw_kind()==primitive_type)
    switch (ag_tp.prim())
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
    break; case K_type_pol_type:
      if (index==KType_type and subscr.specialise(split_type))
        return K_type_poly_term;
    break; case virtual_module_type:
      if (index==param_type and subscr.specialise(split_type))
        return mod_poly_term;
    }
  else if (ag_tp.raw_kind()==row_type and index==int_type and
           subscr.specialise(ag_tp.component_type()))
         return row_entry;
  return not_so;
}
@)
subscr_base::sub_type subscr_base::slice_kind (const type_expr& aggr)
{ type_expr ag_tp = aggr.expanded();
  if (ag_tp.raw_kind()==primitive_type)
    switch (ag_tp.prim())
    {
    case vector_type: return vector_entry;
    case rational_vector_type: return ratvec_entry;
    case string_type: return string_char;
    case matrix_type: return matrix_column;
    default: return not_so;
    }
  else if (ag_tp.raw_kind()==row_type)
    return row_entry;
  else return not_so;
}

@ Some cases, although valid as subscriptions, do not allow a new value to be
assigned to the component value (this holds for instance selecting a character
from a string).

@< Function def... @>=
bool subscr_base::assignable(subscr_base::sub_type t)
{ switch (t)
  { case ratvec_entry: case string_char:
    case not_so: return false;
    default: return true;
  }
}


@ When encountering a subscription in |convert_expr|, we determine the types
of array and of the indexing expression separately, ignoring so far any type
required by the context. Then we look if the types agree with any of the types
of subscription expressions that we can convert to, throwing an error if it
does not. Finally we check if the \foreign{a priori} type |subscr_type| of the
subscripted expression equals or specialises to the required |tp|, or can be
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
{ type array_type = type::bottom(fc),
       index_type = type::bottom(fc);
  auto& subsn = *e.subscription_variant;
    // alias the |struct| pointed to by the |unique_ptr|
  expression_ptr array = convert_expr(subsn.array,array_type);
  expression_ptr index = convert_expr(subsn.index,index_type);
  expression_ptr subscr; type_expr subscr_type;
  @< Set |subscr| to a pointer to a subscription of a kind determined by
     |array_type|, |index_type| while setting |subscr_type|, and holding
     pointers moved from |array| and |index|, or |throw| an error @>
  type result_type = type::wrap(subscr_type,fc);
  return conform_types(result_type,tp,std::move(subscr),e);
}

@ This is a large |switch| statement (the first of several) that is required to
separately and explicitly specify each class template instance that our
program uses (there are $14$ of them here).

The decision whether the subscription is allowed, and what will be the
resulting |subscr_type| are made by the static method |subscr_base::index_kind|.

@< Set |subscr| to a pointer to a subscription of a kind determined by...@>=
switch (subscr_base::index_kind(array_type.bake(),index_type.bake(),subscr_type))
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
case subscr_base::K_type_poly_term:
  if (subsn.reversed)
    throw expr_error(e,"Cannot do reversed subscription of a KTypePol");
  subscr.reset(new K_type_pol_coefficient(std::move(array),std::move(index)));
break;
case subscr_base::mod_poly_term:
  if (subsn.reversed)
    throw expr_error(e,"Cannot do reversed subscription of a ParamPol");
  subscr.reset(new module_coefficient(std::move(array),std::move(index)));
break;
case subscr_base::not_so:
  std::ostringstream o;
  o << "Cannot subscript value of type " << array_type.bake_off() @|
    << " with index of type " << index_type.bake_off();
  throw expr_error(e,o.str());
}

@ For slices we shall similarly need $5$ kinds of slice each with $8$ values
of the template parameter |flags| for $40$ classes in all. Converting a runtime
value |flags| into a template argument (which must be a compile time constant)
can basically only be done by listing all applicable values. To avoid extreme
repetitiveness, we use for this a function that is itself templated over the
class template |slice| that takes |flags| as template argument. In calls to
|make_slice|, its template argument will be specified to one of the $5$
different class templates |row_slice|, |vector_slice|, |ratvec_slice|,
|string_slice|, and |matrix_slice|.

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
{ type array_type = type::bottom(fc),
       index_type = type::wrap(int_type,fc);
  expression_ptr array = convert_expr(e.slice_variant->array,array_type);
  expression_ptr lower =
    convert_expr_strongly(e.slice_variant->lower,fc,int_type);
  expression_ptr upper =
    convert_expr_strongly(e.slice_variant->upper,fc,int_type);
  expression_ptr subscr; const unsigned fl = e.slice_variant->flags.to_ulong();
  switch (subscr_base::slice_kind(array_type.bake()))
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
  return conform_types(array_type,tp,std::move(subscr),e);
    // slicing does not change array type

}


@ Here are the |evaluate| methods for the simpler subscription expressions.
They all follow the same straightforward pattern, and differ only in the way
the result value push on the stack is constructed. The |static_cast<unsigned
int>| allows a range check of the (signed) integer index with a single
comparison against the unsigned array size; in the error message the signed
quantity is transmitted however.

@< Function definitions @>=
inline std::string range_mess
  (arithmetic::Numer_t i,size_t n,const expression_base* e,const char* where)
{ std::ostringstream o;
  e->print(o << "index " << i << " out of range (0<= . <" << n
             << ") in " << where << ' ');
  return o.str();
}
@)
template <bool reversed>
void row_subscription<reversed>::evaluate(level l) const
{ auto i=(index->eval(),get<int_value>()->long_val());
  shared_row r=(array->eval(),get<row_value>());
  size_t n = r->val.size();
  if (static_cast<size_t>(i)>=n)
    throw runtime_error(range_mess(i,n,this,"subscription"));
  push_expanded(l,r->val[reversed ? n-1-i : i]);
}
@)
template <bool reversed>
void vector_subscription<reversed>::evaluate(level l) const
{ auto i=(index->eval(),get<int_value>()->long_val());
  shared_vector v=(array->eval(),get<vector_value>());
  size_t n = v->val.size();
  if (static_cast<size_t>(i)>=n)
    throw runtime_error(range_mess(i,n,this,"subscription"));
  if (l!=level::no_value)
    push_value(std::make_shared<int_value>(v->val[reversed ? n-1-i : i]));
}
@)
template <bool reversed>
void ratvec_subscription<reversed>::evaluate(level l) const
{ auto i=(index->eval(),get<int_value>()->long_val());
  shared_rational_vector v=(array->eval(),get<rational_vector_value>());
  size_t n = v->val.size();
  if (static_cast<size_t>(i)>=n)
    throw runtime_error(range_mess(i,n,this,"subscription"));
  if (l!=level::no_value)
    push_value(std::make_shared<rat_value>(RatNum @|
       (v->val.numerator()[reversed ? n-1-i : i],v->val.denominator())));
}
@)
template <bool reversed>
void string_subscription<reversed>::evaluate(level l) const
{ auto i=(index->eval(),get<int_value>()->long_val());
  shared_string s=(array->eval(),get<string_value>());
  size_t n = s->val.size();
  if (static_cast<size_t>(i)>=n)
    throw runtime_error(range_mess(i,n,this,"subscription"));
  if (l!=level::no_value)
    push_value(std::make_shared<string_value>
      (s->val.substr(reversed ? n-1-i : i,1)));
}

@ And here are the cases for matrix indexing and column selection, which are
just slightly more complicated. The remaining Atlas-related classes
|K_type_coefficient| and |module_coefficient| have their |evaluate| methods
defined in \.{atlas-types.w}.

@< Function definitions @>=
template <bool reversed>
void matrix_subscription<reversed>::evaluate(level l) const
{ index->multi_eval(); @+
  auto j=(get<int_value>()->long_val());
  auto i=(get<int_value>()->long_val());
  shared_matrix m=(array->eval(),get<matrix_value>());
  size_t r = m->val.n_rows();
  size_t c = m->val.n_columns();
  if (static_cast<size_t>(i)>=r)
    throw runtime_error
     ("initial "+range_mess(i,r,this,"matrix subscription"));
  if (static_cast<size_t>(j)>=c)
    throw runtime_error
     ("final "+range_mess(j,c,this,"matrix subscription"));
  if (l!=level::no_value)
    push_value(std::make_shared<int_value>
      (reversed ? m->val(r-1-i,c-1-j): m->val(i,j)));
}
@)
template <bool reversed>
void matrix_get_column<reversed>::evaluate(level l) const
{ auto j=(index->eval(),get<int_value>()->long_val());
  shared_matrix m=(array->eval(),get<matrix_value>());
  size_t c = m->val.n_columns();
  if (static_cast<size_t>(j)>=c)
    throw runtime_error(range_mess(j,c,this,"matrix column selection"));
  if (l!=level::no_value)
    push_value(std::make_shared<vector_value>
      (m->val.column(reversed ? c-1-j : j)));
}

@ For slices these are template functions. This is where the actual reversals
happen.
@< Function definitions @>=
void slice_range_error
  (arithmetic::Numer_t lwb
  ,arithmetic::Numer_t upb
  ,size_t n
  ,unsigned flags,const expression_base* e)
{ std::ostringstream o; bool u_high = upb>=0 and static_cast<size_t>(upb)>n;
  if (u_high)
    if (lwb<0)
      o << "both bounds " << lwb << ':' << upb;
    else
      o << "upper bound " << upb;
  else
    o << "lower bound " << lwb;
  o << " out of range (should be ";
  if (lwb<0)
    o << ">=0" << (u_high ? " respectively " : ")");
  if (u_high)
    o << "<=" << n << ')';
  e->print(o << " in slice ");
  throw runtime_error(o.str());
}
@)
template <unsigned flags>
void row_slice<flags>::evaluate(level l) const
{ auto upb=(upper->eval(),get<int_value>()->long_val());
  auto lwb=(lower->eval(),get<int_value>()->long_val());
  shared_row arr=(array->eval(),get<row_value>());
  const auto& r = arr->val;
  size_t n = r.size();
  if ((flags&0x2)!=0)
    lwb = n - lwb;
  if ((flags&0x4)!=0)
    upb = n - upb;
  if (lwb<0 or upb>=0 and static_cast<size_t>(upb)>n)
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
{ auto upb=(upper->eval(),get<int_value>()->long_val());
  auto lwb=(lower->eval(),get<int_value>()->long_val());
  shared_vector arr=(array->eval(),get<vector_value>());
  const auto& r = arr->val;
  size_t n = r.size();
  if ((flags&0x2)!=0)
    lwb = n - lwb;
  if ((flags&0x4)!=0)
    upb = n - upb;
  if (lwb<0 or (upb>=0 and static_cast<size_t>(upb)>n))
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
{ auto upb=(upper->eval(),get<int_value>()->long_val());
  auto lwb=(lower->eval(),get<int_value>()->long_val());
  shared_rational_vector arr=(array->eval(),get<rational_vector_value>());
  const auto& r = arr->val.numerator();
  size_t n = r.size();
  if ((flags&0x2)!=0)
    lwb = n - lwb;
  if ((flags&0x4)!=0)
    upb = n - upb;
  if (lwb<0 or (upb>=0 and static_cast<size_t>(upb)>n))
    slice_range_error(lwb,upb,n,flags,this);
  if (lwb>=upb)
  @/{@; push_value(std::make_shared<rational_vector_value>(RatWeight(0)));
      return;
  }
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
{ auto upb=(upper->eval(),get<int_value>()->long_val());
  auto lwb=(lower->eval(),get<int_value>()->long_val());
  shared_string arr=(array->eval(),get<string_value>());
  const auto& r = arr->val;
  size_t n = r.size();
  if ((flags&0x2)!=0)
    lwb = n - lwb;
  if ((flags&0x4)!=0)
    upb = n - upb;
  if (lwb<0 or (upb>=0 and static_cast<size_t>(upb)>n))
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
{ auto upb=(upper->eval(),get<int_value>()->long_val());
  auto lwb=(lower->eval(),get<int_value>()->long_val());
  shared_matrix mat=(array->eval(),get<matrix_value>());
  const auto& m = mat->val;
  size_t n = m.n_columns();
  if ((flags&0x2)!=0)
    lwb = n - lwb;
  if ((flags&0x4)!=0)
    upb = n - upb;
  if (lwb<0 or (upb>=0 and static_cast<size_t>(upb)>n))
    slice_range_error(lwb,upb,n,flags,this);
  if (lwb>=upb)
@/{@; push_value(std::make_shared<matrix_value>(int_Matrix(m.n_rows(),0)));
      return; }
  own_matrix result =
    std::make_shared<matrix_value>(int_Matrix(m.n_rows(),upb-lwb));
  auto& r = result->val;
  for (unsigned int j=0; lwb<upb; ++j)
    r.set_column(j,m.column((flags&0x1)==0 ? lwb++ : n - ++lwb));
  push_value(std::move(result));
}

@* Projection functions.
%
The axis language does not have an absolute need for operations of selection
from a tuple, since binding of patterns can achieve the same effect; indeed for
a long time no such operation existed. We introduce them here nonetheless,
because they can be more practical in certain situations, and also because
discriminated unions will have similar injection operations that would be more
cumbersome to do without. These values are also essential for component
assignments.

@< Type def... @>=
struct projector_value : public function_base
{ type_expr type; // used for printing purposes only
  unsigned position;
  id_type id;
  source_location loc;
@)
  projector_value
     (const type_expr& t,unsigned i,id_type id,const source_location& loc)
  : type(t.copy()),position(i),id(id),loc(loc) @+ {}
  virtual ~projector_value() = default;
  virtual void print(std::ostream& out) const;
  virtual void apply(eval_level l) const;
  virtual eval_level argument_policy() const;
  virtual void report_origin(std::ostream& o) const;
  virtual expression_ptr build_call
    (const shared_function& master,const std::string& name,
     expression_ptr&& arg, const source_location& loc) const;
@)
  static const char* name() @+{@; return "built-in function"; }
  projector_value (const projector_value& ) = delete;
};

@ Here are two virtual methods. We print the position selected and the type
selected from. The (tuple) operand of a selection is to be computed as a
|single_value|, since it is unlikely to be given as a tuple display anyway,
and this makes it easy to replace the value on the stack by one of its
components.

@< Function def... @>=
void projector_value::print(std::ostream& out) const
  { out << "{." << main_hash_table->name_of(id) << ": projector_" << position
     @| << '('  << type << ") }"; }
eval_level projector_value::argument_policy() const
  {@; return eval_level::single_value; }
void projector_value::report_origin(std::ostream& o) const
 {@; o << "projector defined " << loc; }

@ Applying a selector is quite easy: we pop the tuple value from the stack,
selecting our component from it, then hand it to |push_expanded| to place the
component back on the stack, according to the policy |l| requested from us.

@< Function def... @>=
void projector_value::apply(eval_level l) const
{@; push_expanded(l,get<tuple_value>()->val[position]); }

@ When a projector is found in the global overload table (as is almost always
the case) its application produces an executable object of class
|projector_call|. It records the projector position and its name.

@< Type def... @>=
struct projector_call : public overloaded_call
{ unsigned position; id_type id;
@)
  projector_call @|
   (const projector_value& f,const std::string& n,expression_ptr&& a
   ,const source_location& loc)
  : overloaded_call(n,std::move(a),loc), position(f.position), id(f.id) @+ {}
  virtual ~projector_call() = default;
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ Here is how a |projector_value| can turn itself into an
|projector_call| when provided with an argument expression, as well
as a |name| to call itself and a |source_location| for the call. We don't store
a (shared) pointer to the |projector_value|, so we ignore the initial
argument here.

@< Function def... @>=
expression_ptr projector_value::build_call
    (const shared_function&,const std::string& name,
     expression_ptr&& arg, const source_location& loc) const
{@; return expression_ptr(new projector_call(*this,name,std::move(arg),loc));
}

@ The virtual methods for |projector_call| are easy. Since nothing can go
wrong with selection, |projector_call::evaluate| does no effort to contribute
to an error back-trace.

@< Function def... @>=
void projector_call::evaluate(level l) const
{@;
  argument->eval(); // evaluate arguments as a single value
  push_expanded(l,get<tuple_value>()->val[position]);
}
@)
void projector_call::print(std::ostream& out) const
{@; out << *argument << '.' << main_hash_table->name_of(id); }

@* Injection functions.
%
For union types there is notion corresponding to projection functions for
tuple types, but dual to it: the injection functions map each of the
constituent types of the union into the union type. Contrary to projection
functions, the use of these is necessary for union types, as the language does
not provide any other means of forming values of these types.

@< Type def... @>=
struct injector_value : public function_base
{ type_expr type; // used for printing purposes only
  unsigned position;
  id_type id;
  source_location loc;
@)
  injector_value
     (const type_expr& t,unsigned i,id_type id,const source_location& loc)
  : type(t.copy()),position(i),id(id),loc(loc) @+ {}
  virtual ~injector_value() = default;
  virtual void print(std::ostream& out) const;
  virtual void apply(eval_level l) const;
  virtual eval_level argument_policy() const;
  virtual void report_origin(std::ostream& o) const;
  virtual expression_ptr build_call
    (const shared_function& master,const std::string& name,
     expression_ptr&& arg, const source_location& loc) const;
@)
  static const char* name() @+{@; return "built-in function"; }
  injector_value (const injector_value& ) = delete;
};

@ Here are two virtual methods. We print the position selected and the type
selected from. The (tuple) operand of a selection is to be computed as a
|single_value|, since it is unlikely to be given as a tuple display anyway,
and this makes it easy to replace the value on the stack by one of its
components.

@< Function def... @>=
void injector_value::print(std::ostream& out) const
  { out << "{."<< main_hash_table->name_of(id) << ": injector_" << position
     @| << '(' << type << ") }"; }
eval_level injector_value::argument_policy() const
  {@; return eval_level::single_value; }
void injector_value::report_origin(std::ostream& o) const
 {@; o << "injector defined " << loc; }

@ Applying an injector is quite easy: we pop the value from the stack, attach
the proper tag (number), and then push that to the stack again.

@< Function def... @>=
void injector_value::apply(eval_level l) const
{ shared_value component = pop_value();
  push_value(std::make_shared<union_value>(position,std::move(component),id));
}

@ Like for projectors, injectors are usually applied in special call
expressions built at ``compile time''.

@< Type def... @>=
struct injector_call : public overloaded_call
{ unsigned position; id_type id;
@)
  injector_call @|
   (const injector_value& f,const std::string& n,expression_ptr&& a
   ,const source_location& loc)
  : overloaded_call(n,std::move(a),loc), position(f.position), id(f.id) @+ {}
  virtual ~injector_call() = default;
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ Here is how a |injector_value| can turn itself into an |injector_call| when
provided with an argument expression, as well as a |name| to call itself and a
|source_location| for the call. We don't store a (shared) pointer to the
|injector_value|, so we ignore the initial argument here.

@< Function def... @>=
expression_ptr injector_value::build_call
    (const shared_function&,const std::string& name,
     expression_ptr&& arg, const source_location& loc) const
{@; return expression_ptr(new injector_call(*this,name,std::move(arg),loc));
}

@ The virtual methods for |injector_call| are easy. Since nothing can go
wrong with injection, the method |injector_call::evaluate| does no effort to
contribute to an error back-trace.

@< Function def... @>=
void injector_call::evaluate(level l) const
{
  argument->eval(); // evaluate arguments as a single value
  push_value(std::make_shared<union_value>(position,pop_value(),id));
}
@)
void injector_call::print(std::ostream& out) const
{@; out << *argument << '.' << main_hash_table->name_of(id); }

@* Control structures.
%
We shall now introduce conventional control structures, which must of course
be part of any serious programming language; yet they were implemented only
after plenty of other language elements were in place, such as
\&{let}-expressions, functions, rows and selection form them, implicit
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
to it, but nothing else is so) and so must (become) the required |tp|.

@< Cases for type-checking and converting... @>=
case negation_expr:
{ expression_ptr arg = convert_expr_strongly(*e.negation_variant,fc,bool_type);
  if (not tp.unify_to(type::wrap(bool_type,fc)))
    // |not| preserves the |bool| type
    throw type_error(e,bool_type.copy(),tp.bake());
  return expression_ptr(new @| builtin_call
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
  virtual ~conditional_expression() = default;
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
    auto* p =
      dynamic_cast<const conditional_expression*>(cur->else_branch.get());
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
  expression_ptr c = convert_expr_strongly(exp.condition,fc,bool_type);
  std::vector<expression_ptr> conv;
  balance(tp,&exp.branches,e,"branches of conditional",conv);
@/return expression_ptr(new @|
    conditional_expression(std::move(c),std::move(conv[0]),std::move(conv[1])));
}

@*1 Integer controlled case expressions.
%
The integer case expression (multi-way branch controlled by an integer value)
is quite similar to a conditional, but contains a list of branches rather than
two of them. It used to have just that, but it turns out to be useful to give
the user control of how to handle out of range selection values; for that
purpose $0$, $1$, or $2$ additional expressions can be provided as an
alternative for the listed branches.

@< Type def... @>=
struct int_case_expression : public expression_base
{ expression_ptr condition; std::vector<expression_ptr> branches;
@)
  int_case_expression
   (expression_ptr&& c,std::vector<expression_ptr>&& b)
   : condition(c.release()),branches(std::move(b))
  @+{}
  virtual ~int_case_expression() = default;
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};
struct int_case_else_expression : public expression_base
{ expression_ptr condition;
  expression_ptr out_branch;
  std::vector<expression_ptr> branches;
@)
  int_case_else_expression
   (expression_ptr&& c,expression_ptr&& o,std::vector<expression_ptr>&& b)
   : condition(c.release()),out_branch(std::move(o)),branches(std::move(b))
  @+{}
  virtual ~int_case_else_expression() = default;
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};
struct int_case_then_else_expression : public expression_base
{ expression_ptr condition;
  expression_ptr pre_branch, post_branch;
  std::vector<expression_ptr> branches;
@)
  int_case_then_else_expression
   (expression_ptr&& c,expression_ptr&& t,expression_ptr&& e,
    std::vector<expression_ptr>&& b)
   : condition(c.release())
   , pre_branch(std::move(t))
   , post_branch(std::move(e))
   , branches(std::move(b))
  @+{}
  virtual ~int_case_then_else_expression() = default;
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
void int_case_else_expression::print(std::ostream& out) const
{ auto it = branches.cbegin();
  assert(it!=branches.cend());
  out << " case " << *condition << " in " << **it;
  while (++it!=branches.cend())
    out << ", " << **it;
  out << " else " << *out_branch << " esac ";
}
void int_case_then_else_expression::print(std::ostream& out) const
{ auto it = branches.cbegin();
  assert(it!=branches.cend());
  out << " case " << *condition << " then " << *pre_branch << " in " << **it;
  while (++it!=branches.cend())
    out << ", " << **it;
  out << " else " << *post_branch << " esac ";
}

@ Evaluating a case expression ends up evaluating one of the |branches|. In all
cases a condition that indexes one of the \&{in} branches selects that branch.
After conversion to an unsigned type, a single comparison suffices to detect
this case. For an out-of-bounds condition, if no clause is provided for that
case, we apply |arithmetic::remainder| and select the corresponding \&{in}
branch. Otherwise we select the \&{else} branch, or in case of a negative
condition the \&{then} branch if one was provided.

@< Function definitions @>=
void int_case_expression::evaluate(level l) const
{ condition->eval();
  auto i = get<int_value>();
  if (i->val.size()==1) // in normal cases avoid using |shift_modulo|
  { auto ii = i->int_val();
    if (ii >= 0 and static_cast<unsigned>(ii) < branches.size())
      branches[ii]->evaluate(l);
    else
      branches
       [arithmetic::remainder(i->int_val(),static_cast<int>(branches.size()))]
      ->evaluate(l);
  }
  else
    branches[big_int(i->val).shift_modulo(branches.size())]->evaluate(l);
}
void int_case_else_expression::evaluate(level l) const
{ condition->eval();
  auto i = get<int_value>();
  if (i->val.size()==1) // necessary condition for selecting an \&{in} branch
  { auto ii = i->int_val();
    if (ii >= 0 and static_cast<unsigned>(ii) < branches.size())
      return branches[ii]->evaluate(l);
  }
  out_branch->evaluate(l);
}
void int_case_then_else_expression::evaluate(level l) const
{ condition->eval();
  auto i = get<int_value>();
  if (i->val.is_negative())
    return pre_branch->evaluate(l);
  if (i->val.size()==1) // necessary condition for selecting an \&{in} branch
  { auto ii = i->int_val();
    if (ii >= 0 and static_cast<unsigned>(ii) < branches.size())
      return branches[ii]->evaluate(l);
  }
  post_branch->evaluate(l);
}

@ With the function |balance| defined above, conversion of case expressions
has become easy.

@< Cases for type-checking and converting... @>=
case int_case_expr0:
case int_case_expr1:
case int_case_expr2:
{ auto& exp = *e.if_variant;
  expression_ptr c = convert_expr_strongly (exp.condition,fc,int_type);
  std::vector<expression_ptr> conv;
  balance(tp,&exp.branches,e,"branches of case",conv);
  if (e.kind==int_case_expr0)
    return expression_ptr(new @|
      int_case_expression(std::move(c),std::move(conv)));
  expression_ptr p0(std::move(conv[0]));
  conv.erase(conv.begin());
  if (e.kind==int_case_expr1)
    return expression_ptr(new @|
      int_case_else_expression(std::move(c),std::move(p0),std::move(conv)));
  expression_ptr p1(std::move(conv[0]));
  conv.erase(conv.begin());
@/return expression_ptr(new @|
    int_case_then_else_expression(std::move(c),std::move(p0),std::move(p1),
                                 std::move(conv)));
}

@*1 Discrimination clauses.
%
The discrimination clause (a multi-way branch controlled by a union value) is
syntactically similar to the integer case expression, but semantically a bit
more complicated. The main difference is that each branch can bind, depending on
its actual variant, to the value of union type being discriminated upon. We
define two syntactic forms of discrimination clauses, a simple one where each
branch is like a \&{let}-expression and they bind to the variants in order, and
a more elaborated form where the branches identify their corresponding variant
of the union type by a tag associated to it. The latter form allows stating the
branches in any order, and using a default branch to replace one or more
branches; its usage presupposes the user having given a (general) type
definition for the union type, which associates tags to each of the variants.

The difference between the two forms is merely syntactic; after conversion, both
forms will give rise to a |discrimination_expression|. Unless a default branch
is present (which cannot be used in the first form), the converted form in no
way witnesses the form that was used to produce it. Note that exceptionally this
structure uses shared pointers for the branches. This is because a default
branch expression can stand for multiple branches, in which case several copies
of the pointer to it will be present in |branches|.

@< Type def... @>=
using choice_part = std::pair<id_pat,shared_expression>;

struct discrimination_expression : public expression_base
{ expression_ptr subject; std::vector<choice_part> branches;
@)
  discrimination_expression
   (expression_ptr&& subj,std::vector<choice_part>&& br)
   : subject(subj.release()),branches(std::move(br))
  @+{}
  virtual ~discrimination_expression() = default;
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ To print a discrimination clause is straightforward.

@< Function definitions @>=
void discrimination_expression::print(std::ostream& out) const
{ auto it = branches.cbegin();
  assert(it!=branches.cend());
  out << " case " << *subject;
  do out << " | (" << it->first << "):" << *it->second;
  while (++it!=branches.cend());
  out << " esac ";
}

@ Evaluating a discrimination clause begins by evaluating the |subject| (the
expression upon which one discriminates), pops that value from the stack into
|discriminant|, and then selects a |branch| depending on its (numeric)
|variant()| attribute. Then that branch is evaluated like a \&{let}-expression
introducing an identifier pattern |branch.first|, which is bound to the value of
|subject| within its variant, in which context the body |branch.second| is
evaluated to produce the result of the discrimination clause. However we must be
aware that some branches bind no identifiers (for instance this is always the
case for a defaulted branch); we test for that before creating a |frame| for
|branch.first|.

@< Function definitions @>=
void discrimination_expression::evaluate(level l) const
{ subject->eval();
  shared_union discriminant = get<union_value>();
  const auto& branch = branches[discriminant->variant()]; // make selection
  if (count_identifiers(branch.first)==0)
    branch.second->evaluate(l); // avoid any empty |frame|
  else
  { frame fr(branch.first);
    fr.bind(discriminant->contents());
    try @+{@; branch.second->evaluate(l); }
    @< Catch block for providing a trace-back of local variables @>
  }
}


@ When type checking discrimination clauses, the case where |has_tags| holds
is the more complicated one: we use the tags associated to the variants of the
union type to identify and possibly reorder branches, and allow for some
branches to be handled together in a default branch (in which case no
information from the evaluated subject expression, other than the fact that it
selects none of the explicitly treated branches, is passed to the selected
default branch). So we insist that the union type used be present
in |type_expr::type_map| (presumably specifying tags) in that case.

@< Cases for type-checking and converting... @>=
case discrimination_expr:
{ auto& exp = *e.disc_variant;
  size_t n_branches=length(&exp.branches);
  std::ostringstream o; // for use in various error messages
@)
  type subject_type = type::bottom(fc);
  expression_ptr c  =  convert_expr(exp.subject,subject_type);
  if (subject_type.expunge().expand().kind()!=union_type)
    @< Report that a discrimination clause needs to be a union type,
       |exp.subject| has non-union type |subject_type| @>
  const auto* variants = subject_type.tuple();
  size_t n_variants=length(variants);
  std::vector<choice_part> choices(n_variants);
@)
  if (exp.has_tags())
  {
    sl_list<id_type> tags = exp.tag_set();
    auto candidates = type_expr::matching_bindings(subject_type);
    if (candidates.empty())
      @< Report that union type |subject_type| cannot be used in a discrimination
         clause without having been defined @>
    bool unique_candidate = candidates.singleton();
    simple_list<unsigned int> positions; // to record indexes of the |tags|
    @< Filter from |candidates| those that do not know all |tags|, and report an
       error if this does not reduce it to a singleton;
       otherwise set |positions| according to that remaining candidate @>
    const wtl_const_iterator candidate_var_tps(candidates.front()->tp.tuple());
  @)
    @< Process |branches| and assign them to |choices|, possibly reordering them
      according to the specified variant names and taking into account a
      possible default branch @>
  }
  else
    @< Process |branches| in order @>
@/return expression_ptr(new @|
    discrimination_expression(std::move(c),std::move(choices)));
}

@ We start with the easiest case, namely where the user does not use any tags,
but simply lists branches in the order of the union type. Processing then
proceeds by sequentially following variants of |subject_type| which should match
the branches. We do set up this code (notably using an index $k$ rather than yet
another iterator) in such a way that the actual branch processing can be shared
with the case with tags.

@< Process |branches| in order @>=
{ if (n_branches!=n_variants)
    @< Report mismatching number of branches @>
  case_list::weak_const_iterator branch_it(&exp.branches);
  wtl_const_iterator type_it(variants);
  for (size_t k=0; k<n_branches; ++k,++branch_it,++type_it)
  { const auto& branch = *branch_it;
    const type variant_type =
      type::wrap(*type_it,fc); // type of current variant of the union
    @/@< Type-check branch |branch.branch|, with the names of |branch.pattern|
         bound to |variant_type|, requiring result type |tp|; insert the
         resulting |choice_part| value into |choices[k]| @>
  }
}

@ When tags were specified, their set is given in |tags| and the candidate union
types that were retrieved from |type_expr::type_map| are collected in the list
|candidates|. Before insisting that there is a unique candidate that governs the
correspondence between tags and positions in the union, we remove and candidates
whose set of tags do not provide for all tags that were used in the clause.

@< Filter from |candidates| those that do not know all |tags|,...@>=
{ for (auto it = candidates.begin(); not candidates.at_end(it); )
  // no |++it| here!
  { sl_list<unsigned int> pos;
    const std::vector<id_type>& fields = (*it)->fields;
    assert(fields.size()==n_variants);
      // because |(*it)->tp| unifies with |subject_type|
    bool failure = false;
    for (id_type tag : tags)
    { unsigned int i;
      for (i=0; i<n_variants; ++i)
        if (fields[i]==tag)
        {@;
          pos.push_back(i);
          break;
        }
      if (i==n_variants) // then |tag| was not found
      { if (unique_candidate)
          @< Report that |tag| does not match any of the tags of our unique
             candidate type @>
        failure = true;
        break;
      }
    } // |for(tag)|
    if (failure)
      candidates.erase(it);
    else
    { if (positions.empty())
        positions=pos.undress();
      ++it;
    }
  } // |for (it)| looping over |candidates|
  if (not candidates.singleton())
    @< Report that the number of |candidates| accommodating all |tags|
       is not exactly one @>
}

@ In order to be able to share the default branch expression among multiple
branches, we define a shared expression |default_choice| that a default branch
if present will set.

The branches of a discrimination clause should have bodies of equal types. It
would then seem natural to balance these types like for an integer case
expression, but our current implementation of balancing is not up to the task.
This is because it will convert, and occasionally re-convert, all expressions to
be balanced in the same (lexical) context, while the branches of a
discrimination clause set up different local bindings for each of them. So we
implement conversion of discrimination clauses without balancing, simply using
for each branch the type provided by the context as possibly modified by the
conversion of previous branches.

@< Process |branches| and assign them to |choices|... @>=
{ shared_expression default_choice;
  auto pos_it = positions.cbegin();
  for (case_list::weak_const_iterator br_it(&exp.branches);
       not br_it.at_end(); ++br_it)
  { const auto& branch = *br_it;
    if (branch.is_default())
      @< Use |branch| to set |default_choice| @>
    else
    { size_t k = *pos_it;
      ++pos_it;
      @< Check that |choices[k]| has not been filled before @>
      const wtl_const_iterator types_start(variants); // base for ``indexing''
      type variant_type =
        type::wrap(*std::next(types_start,k),fc);
      // type of variant for this branch
      type tag_type = type::wrap(*std::next(candidate_var_tps,k),0,fc);
      bool success = variant_type.unify_to(tag_type);
      assert(success); ndebug_use(success);
    @/@< Type-check branch |branch.branch|, with the names of |branch.pattern|
         bound to |variant_type|, requiring result type |tp|; insert the
         resulting |choice_part| value into |choices[k]| @>
    }
  }
  @< Check valid use of default branch, and insert it into the empty slots
      of |choices| @>
}

@ Type-checking a branch is easy once everything is put in place. We have the
identifier pattern for the branch in |branch.pattern|, and its type in
|variant_type|, so the |layer| data type and its |thread_bindings| method will
take care of setting up the evaluation context, and then |convert_expr| will
to the actual type checking. We just need to take care to ensure that an
entirely empty |id_pat| is recorded whenever |branch.pattern| does not bind
any identifiers at all, which the |evaluate| method then will recognise if
this branch is chosen, and suppress creating a |frame| for the branch.

@< Type-check branch |branch.branch|, with the names of |branch.pattern|... @>=
{ if (not variant_type.unwrap().can_specialise(pattern_type(branch.pattern)))
  { o << "Pattern " << branch.pattern @|
      << " does not match type " << variant_type @|
      << " for variant " <<  main_hash_table->name_of(branch.label);
    throw expr_error(e,o.str());
  }
@)
  expression_ptr result;
  layer branch_layer(count_identifiers(branch.pattern));
  thread_bindings(branch.pattern,variant_type.unwrap(),fc,branch_layer,false);
  result=convert_expr(branch.branch,tp);
@/choices[k] = choice_part (copy_id_pat(branch.pattern),std::move(result));
}

@ When a pair of identical branch labels is found, we print the whole
discrimination expression.

@< Check that |choices[k]| has not been filled before... @>=
{ if (choices[k].second.get()!=nullptr)
  { o << "Multiple branches with label "
      << main_hash_table->name_of(branch.label);
    throw expr_error(e,o.str());
  }
}

@ Reporting a non-union subject is done by producing a union type with unknown
components, and throwing a |type_error| with that as required type. Since we
don't know which union type the user intended, we produce their number by
counting branches od the discrimination clause, but ensuring that this number is
at least~$2$ so that a union type can be recognised in the error message.

@< Report that a discrimination clause needs to be a union type,
       |exp.subject| has non-union type |subject_type| @>=

{ type_list l; auto n = std::max<size_t>(2,n_branches);
  while (n-->0)
    l.push_front(unknown_type.copy());
  throw type_error(exp.subject,subject_type.bake_off(),
                   type_expr::tuple_or_union(union_type,std::move(l)));
}

@ The error message for a wrong number of branches just uses the branch count,
but does report the full union type for clarity rather than just its number of
variants.

@< Report mismatching number of branches @>=
{ std::ostringstream o;
  o << "Discrimination clause has " << n_branches @|
    << "branches, which does not match\nthe number of variants ("
    << n_variants @| << ") of the type " << subject_type
    << " discriminated upon";
  throw expr_error(e,o.str());
}

@ Here our error message tries to suggest the two solutions to attempting to use
a union type without declared tags: either declare them or use the tag-less
version of a discrimination clause.

@< Report that union type |subject_type| cannot be used in a discrimination
   clause without having been defined @>=

{ std::ostringstream o;
  o << "Discrimination on expression of type " << subject_type @|
    << " with clause using tags, but none are known.\n" @|
     "  Either use 'set_type' with tag names first, " @|
     "or use discrimination clause without tags.";
  throw expr_error(e,o.str());
}

@ When there is an unambiguous defined union type matching |subject_type|, we
report the first tag that fails to match it as the offender.

@< Report that |tag| does not match any of the tags... @>=
{ std::ostringstream o;
  o << "Identifier " << main_hash_table->name_of(tag) @|
    << " is not a tag associated with the union type" @|
    << (candidates.front()->arity==0 ? " " : " constructor " ) ;
  if (candidates.front()->name==type_binding::no_id)
    o << candidates.front()->tp;
  else
    o << main_hash_table->name_of(candidates.front()->name);
  throw expr_error(e,o.str());
}

@ We already tested there were defined union types, or type constructors, with
tags, that match the |subject_type| of the discrimination clause, so if we find
there are none left, the error message should focus on the tags that failed to
all match for one of the candidates.

@< Report that the number of |candidates| accommodating all |tags| is not
   exactly one @>=
{ std::ostringstream o;
  if (candidates.empty())
  {
    o << "No union definition found to accommodate the tag" @|
      << (tags.singleton() ? ": '" : "s: '") @|
      << main_hash_table->name_of(tags.front());
    for (auto it=std::next(tags.begin()); not tags.at_end(it); ++it)
      o << "', '" << main_hash_table->name_of(*it);
    o << "',\n  used in discrimination clause";
  }
  else
  { o << "Ambiguity in discrimination clause, possible types are:\n";
    for (const auto& cand : candidates)
    { o << "    ";
      if (cand->name==type_binding::no_id)
         o << cand->tp << '\n';
      else
        o << main_hash_table->name_of(cand->name) << '\n';
    }
  }
  throw expr_error(e,o.str());
}

@ A default branch is just an expression, which we store unless a default
branch was already defined.

@< Use |branch| to set |default_choice| @>=
{ if (default_choice.get()!=nullptr)
    throw expr_error(e,"Multiple default branches present");
  default_choice = convert_expr(branch.branch,tp);
}

@ A default branch should be present if and only if the number of other
branches specified is strictly less than the number of variants, which means
that the number |n_branches| (which includes a possible default) should never
exceed |n_branches|. We also test for possible a missing default branch.

When no errors are signalled, we need to implement the default branch
mechanism, which is easy thanks to our use of |shared_expression|. By having
the identifier pattern of defaulted branches be empty, the value from the
|discriminant| expression will not be bound to anything for such branches.

@< Check valid use of default branch, and insert it into the empty slots
   of |choices| @>=
{
  if (n_branches>n_variants)
    throw expr_error(e,"Spurious default branch present");
  if (default_choice.get()!=nullptr)
  // if a default was given, insert for omitted branches
  { for (auto it=choices.begin(); it!=choices.end(); ++it)
      if (it->second.get()==nullptr)
        *it = choice_part(id_pat(),default_choice);
              // share |default_choice| here
  }
  else if (n_branches<n_variants)
    @< Report a missing branch @>
}

@ If |n_branches<n_variants| and no default branch was seen, we give an error
message that reports a variant of the union for which a branch is actually
missing. We mention its tag if it has one, and otherwise mention its position
among the variants.

@< Report a missing branch @>=
{
  auto it=choices.begin();
  while(it->second.get()!=nullptr) ++it;
  auto tag = candidates.front()->fields[it-choices.begin()];
  if (tag==type_binding::no_id)
    o << "Missing branch for anonymous variant " << it-choices.begin();
  else
    o << "Missing branch for variant " << main_hash_table->name_of(tag);
  o << " of type " << subject_type << " in discrimination clause";
  throw expr_error(e,o.str());
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
constexpr bool has_frame(unsigned flags) @+{@; return (flags&0x4)==0; }
constexpr bool yields_count(unsigned flags) @+{@; return (flags&0x8)!=0; }

@ Let us start with considering |while| loops.
%
Although they contain two parts, a condition before |do| and a body after it,
the two are present in the following structure as a single |body|. The reason
for this is that we want to allow declarations in the condition to remain valid
in the body, which can be realised by having both contained in a |do_expression|
structure. That expression does not have to be a direct descendent of the
|while_expression|, but its relation can involve descending to the right child
of a |let_expression| or |seq_expression| nodes any number of times, or passing
into a branch of a choice expression. The grammar guarantees that eventually
such a |do_expression| will be reached, and the same for each branch in the case
of choice expressions.

@< Type def... @>=
template<unsigned flags>
struct while_expression : public expression_base
{ expression_ptr body;
  while_expression(expression_ptr&& b): body(std::move(b)) @+{}
  virtual ~while_expression() = default;
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ In printing of while loops, the symbol \.{do} will appear in the output for
|*body| (possibly more than once, if a choice expression is involved).

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
  layer bind(0,nullptr); // no local variables for loop, but allow \&{break}
@)
  if (tp.is_void())
    return expression_ptr@|(make_while_loop (w.flags.to_ulong(),
       convert_expr_strongly(w.body,fc,void_type)));
  else if (tp.bake()==int_type)
    return expression_ptr(make_while_loop @| (0x8,
       convert_expr_strongly(w.body,fc,void_type)));
  type_expr loop_type = row_of_type.copy();
  if (tp.unify_specialise(loop_type))
  { type comp_type = type::wrap(loop_type.component_type(),fc);
    expression_ptr b = convert_expr(w.body,comp_type);
    if (comp_type.is_void() and not is_empty(w.body))
      b.reset(new voiding(std::move(b)));
    return expression_ptr (make_while_loop @|
       (w.flags.to_ulong(),std::move(b)));
  }
  else
  @< If |tp| can be converted from some row-of type, check |w.body|
     against its component type, construct the |while_expression|, and apply
     the appropriate conversion function to it; otherwise |throw| a
     |type_error| @>
}

@ For |while| loops we follow the same logic for finding an appropriate
component type as for list displays, in section@#list display conversion@>.

@< If |tp| can be converted from some row-of type, check |w.body| against
   its component type, construct the |while_expression|, and apply the
   appropriate conversion function to it; otherwise |throw| a |type_error| @>=
{ type_expr comp_type;
  const conversion_record* conv = row_coercion(tp.bake(),comp_type);
  if (conv==nullptr)
    throw type_error(e,row_of_type.copy(),tp.bake());
@)
  return expression_ptr(new conversion(*conv, expression_ptr(make_while_loop @|
       (w.flags.to_ulong(),convert_expr_strongly(w.body,fc,comp_type)))));
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
  virtual ~do_expression() = default;
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
  virtual ~forever_expression() = default;
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};
struct dont_expression : public expression_base
{ dont_expression() @+{}
  virtual ~dont_expression() = default;
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
  expression_ptr body = convert_expr(seq.last,tp);
    // body needs type-checking in all cases
  if (seq.first.kind==boolean_denotation)
  { if (seq.first.bool_denotation_variant==neg)
      return expression_ptr(new dont_expression()); // and drop |body|
    else
      return expression_ptr(new forever_expression(std::move(body)));
  }
  expression_ptr condition = convert_expr_strongly(seq.first,fc,bool_type);
  if (neg)
    return expression_ptr(new @|
      do_expression<true>(std::move(condition),std::move(body)));
  else
    return expression_ptr(new @|
       do_expression<false>(std::move(condition),std::move(body)));
}
break;

@ The evaluation of a |do_expression| is somewhat special in that whether it
places a value on the runtime stack (assuming that |l!=level::no_value|) depends on
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
{ if (l==level::no_value)
  { try
    { do
      { body->void_eval();
        check_interrupt();
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
        check_interrupt();
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
    { check_interrupt();
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
|make_pattern_list| that is used takes the \emph{tail} of the list as first
argument, and the head as second; this is adapted to the subsequent reversal that
usually takes place when a pattern list is completed, but this reversal does
not happen for the $2$-element list used for the patterns in for-loops.)

Apart from iterating over row value of any type, we allow iteration over
vectors, rational vectors, strings, and matrix columns (types that are indexable
by integers), and also over the terms of a $K$-type polynomial or a parameter
polynomial (also called virtual module, of which each term represents an
isotypical component). The only thing we can select by subscription but not loop
over is individual matrix entries (the corresponding |subscr_base::sub_type|
value is |matrix_entry|). The syntax of the for loop is the same for all
these cases. However the constructed |expression| will be of a class templated
on |kind| describing the kind of value iterated over, so that their |evaluate|
methods will be specialised to that kind. We also template over the |flags| that
indicate the reversal options (determined by by the precise syntax used), so
that these get built into the evaluate method as well. Bit |0x1| of |flags|
controls reverse traversal at input, while bit |0x2| indicates reverse
accumulation of values for the result. Altogether we use $4\times7=28$ instances
of |for_expression|.

@< Type def... @>=
template <unsigned flags, subscr_base::sub_type kind>
struct for_expression : public expression_base
{ id_pat pattern; expression_ptr in_part;
  expression_ptr body;
@)
  for_expression (const id_pat& p, expression_ptr&& i, expression_ptr&& b);
  virtual ~for_expression() = default;
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
{ out << (in_reversed(flags) ? " ~do " : " do ")  << *body
  @|  << (out_reversed(flags) ? " ~od " : " od ");
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
  case subscr_base::K_type_poly_term:
    switch (flags)
    {
    case 0: return new for_expression<@[0,subscr_base::K_type_poly_term@]>
      (id,std::move(i),std::move(b));
    case 1: return new for_expression<@[1,subscr_base::K_type_poly_term@]>
      (id,std::move(i),std::move(b));
    case 2: return new for_expression<@[2,subscr_base::K_type_poly_term@]>
      (id,std::move(i),std::move(b));
    case 3: return new for_expression<@[3,subscr_base::K_type_poly_term@]>
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
its a priori type. Then after binding the loop variable(s) in a new |layer|, we
process the loop body requiring either the type |*tp.component_type()| if |tp|
is a row type, of void if the loop occurs in a void context, or a type that
|row_coercion| finds will allow a subsequent coercion of the row type to~|tp|.
If none of these apply a |type_error| is thrown, which indicates the row type of
the first option as expected type. After converting the loop, we must not forget
to maybe apply voiding (if the body has void type that is not due to a void
context for the loop itself) or a coercion.

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
  type in_type = type::bottom(fc);
  expression_ptr in_expr = convert_expr(f.in_part,in_type);  // \&{in} part
  subscr_base::sub_type which; // the kind of aggregate iterated over
  layer bind(count_identifiers(f.id),nullptr);
   // for identifier(s) introduced in this loop
  @< Set |which| according to |in_type|, and set |bind| according to the
     identifiers contained in |f.id| @>

  type body_type = type::bottom(fc);
  const conversion_record* conv=nullptr;
    // initialise in case we don't reach the assignment below
  type_expr tmp_tp = row_of_type.copy();
  bool do_export = tp.unify_specialise(tmp_tp);
  if (do_export)
    body_type = type::wrap(tmp_tp.component_type(),fc);
  else if (tp.is_void())
    body_type =tp.copy();
  else if ((conv=row_coercion(tp.bake(),tmp_tp))!=nullptr)
    body_type = type::wrap(tmp_tp,fc);
  else throw type_error(e,row_of_type.copy(),tp.bake());
@)
  expression_ptr body(convert_expr(f.body,body_type));
  if (do_export)
    tp.unify_to(body_type.wrap_row());
  if (body_type.is_void() and not tp.is_void() and not is_empty(f.body))
    body.reset(new voiding(std::move(body)));
@/expression_ptr loop;
  if (bind.empty()) // we must avoid having an empty |frame| produced at runtime
    @< Set |loop| to a index-less counted |for| loop controlled by the
       size of |in_expr|, and containing |body| @>
  else loop.reset(make_for_loop@|
    (f.flags.to_ulong(),f.id,std::move(in_expr),std::move(body),which));
@/return conv==nullptr ? std::move(loop)
  : @| expression_ptr(new conversion(*conv,std::move(loop))) ;
}

@ The |in_type| must be indexable by integers (so it is either a row-type or
vector, rational vector, matrix or string), or it must be the polynomial type,
indexable by |param_type|. If so, the first or second call to
|subscr_base::index_kind| will set |comp_type| to the component type resulting
from such a subscription. We also make |tp| point to the index type used.

@< Set |which| according to |in_type|, and set |bind| according to the
   identifiers contained in |f.id| @>=
{ type_expr in_tp = in_type.bake(), comp_type;
  const_type_p inx_type;
  which = subscr_base::index_kind(in_tp,*(inx_type=&int_type),comp_type);
  if (which==subscr_base::not_so)
    // if not integer-indexable, try $K$-type indexable
    which = subscr_base::index_kind(in_tp,*(inx_type=&KType_type),comp_type);
  if (which==subscr_base::not_so)
    // if not integer-indexable, try parameter indexable
    which = subscr_base::index_kind(in_tp,*(inx_type=&param_type),comp_type);
  if (which==subscr_base::not_so) // if its not that either, it is wrong
  { std::ostringstream o;
    o << "Cannot iterate over value of type " << in_type;
    throw expr_error(e,o.str());
  }
  type_list it_comps; // type of "iterator" value (pattern) named in the loop
  it_comps.push_front(std::move(comp_type));
  it_comps.push_front(type_expr(inx_type->copy()));
  const type_expr it_type = type_expr::tuple(std::move(it_comps));
    // build tuple type from index and component types
  if (not it_type.can_specialise(pattern_type(f.id)))
    throw expr_error(e,"Improper structure of loop variable pattern");
  thread_bindings(f.id,it_type,fc,bind,true); // force all identifiers constant
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
    if (l!=level::no_value)
    { if (out_reversed(flags))
        // doing \&{break} in reverse-gathering loop requires a shift
        dst=std::move(dst,result->val.end(),result->val.begin());
          // left-align |result|
      result->val.resize(dst-result->val.begin());
      // resize needed in all cases
    }
  }

  if (l!=level::no_value)
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
    try {
      while (i!=(in_forward(flags) ? n : 0))
      { loop_var->val[1]=in_val->val[in_forward(flags) ? i : i-1];
          // move index into |loop_var| pair
        @< Set |loop_var->val[0]| to |i++| or to |--i|, create a new |frame| for
        |pattern| binding |loop_var|, and evaluate the |loop_body| in it;
        maybe assign |*dst++| or |*--dst| from it @>
      }
    }
    @< Catch block for reporting iteration number within loop that threw @>
  }
  @+break;
case subscr_base::vector_entry:
  { shared_vector in_val = get<vector_value>();
    size_t n=in_val->val.size();
    @< Define loop index |i|, allocate |result| and initialise iterator |dst| @>
    try {
      while (i!=(in_forward(flags) ? n : 0))
      { loop_var->val[1] = std::make_shared<int_value>
          (in_val->val[in_forward(flags) ? i : i-1]);
        @< Set |loop_var->val[0]| to... @>
      }
    }
    @< Catch block for reporting iteration number within loop that threw @>
  }
  @+break;
case subscr_base::ratvec_entry:
  { shared_rational_vector in_val = get<rational_vector_value>();
    size_t n=in_val->val.size();
    @< Define loop index |i|, allocate |result| and initialise iterator |dst| @>
    try {
      while (i!=(in_forward(flags) ? n : 0))
      { loop_var->val[1] = std::make_shared<rat_value> @|
        (RatNum
          (in_val->val.numerator()[in_forward(flags) ? i : i-1]
          ,in_val->val.denominator()));
        @< Set |loop_var->val[0]| to... @>
      }
    }
    @< Catch block for reporting iteration number within loop that threw @>
  }
  @+break;
case subscr_base::string_char:
  { shared_string in_val = get<string_value>();
    size_t n=in_val->val.size();
    @< Define loop index |i|, allocate |result| and initialise iterator |dst| @>
    try {
      while (i!=(in_forward(flags) ? n : 0))
      { loop_var->val[1] = std::make_shared<string_value>
              (in_val->val.substr(in_forward(flags) ? i : i-1,1));
        @< Set |loop_var->val[0]| to... @>
      }
    }
    @< Catch block for reporting iteration number within loop that threw @>
  }
  @+break;

@ Here are the remaining cases. The case |matrix_column| is ever so slightly
different because the iteration count~|n| is given by the number of columns,
rather than the size, of |inv_val->val|. The cases |K_type_poly_term| and
|mod_poly_term| have more important differences, to be detailed later. The other
two cases should never arise.

@< Cases for evaluating a loop over components of a value... @>=
case subscr_base::matrix_column:
  { shared_matrix in_val = get<matrix_value>();
    size_t n=in_val->val.n_columns();
    @< Define loop index |i|, allocate |result| and initialise iterator |dst| @>
    try {
      while (i!=(in_forward(flags) ? n : 0))
      { loop_var->val[1] = std::make_shared<vector_value>
          (in_val->val.column(in_forward(flags) ? i : i-1));
        @< Set |loop_var->val[0]| to... @>
      }
    }
    @< Catch block for reporting iteration number within loop that threw @>
  }
  @+break;
case subscr_base::K_type_poly_term:
  @< Perform a loop over the terms of a $K$-type polynomial @>
  break;
case subscr_base::mod_poly_term:
  @< Perform a loop over the terms of a virtual module @>
  break;
case subscr_base::matrix_entry:; // excluded in type analysis
case subscr_base::not_so: assert(false);

@ The following code, which occurs five times, uses both the input and output
direction attributes.

@< Define loop index |i|, allocate |result| and initialise iterator |dst| @>=
size_t i= in_forward(flags) ? 0 : n;
if (l!=level::no_value)
{ result = std::make_shared<row_value>(n);
  dst = out_forward(flags) ? result->val.begin() : result->val.end();
}

@ Correspondingly, when an error is thrown from the loop body, we indicate the
iteration number in the back trace. This code too is shared among the initial
$5$ kinds of for-loop evaluators. We are using |while| loops in which the
index~|i| has been increased or decreased \emph{before} the loop body is
executed, so the code below compensates for that to print the actual count
(from~$0$) of the iteration.

@< Catch block for reporting iteration number within loop that threw @>=
catch (error_base& e)
{
  std::ostringstream o;
  o << "During iteration " << (in_forward(flags) ? i-1 : n-1-i) @|
    << " of the " << (in_forward(flags) ? "" : "reversed ") << "for-loop";
  e.trace(o.str());
  throw;
}

@ This code too occurs identically five times. We have separately set the
in-part component stored in |loop_var->val[1]| for the various values of |kind|,
but |loop_var->val[0]| is always the (integral) loop index. Once |loop_var|
initialised, we set up a new |frame| (named |fr| so as to be able to reuse a
catch block textually) according to |pattern|, which pushes a new
|evaluation_context| onto |frame::current|. Like it was important for
|loop_var->val[0]| to build a fresh |shared_value|, it is important that a fresh
|evaluation_context| be created at each iteration, rather than using a single
one for the entire loop. Any closure values formed inside the loop body will
incorporate the |evaluation_context| by reference; these will only be able to
access distinct instances of |loop_var| in the absences of aliasing either at
the level |evaluation_context| of the |shared_value| pointers stored in it. With
these things properly handled, the evaluation of the loop body is standard.

@< Set |loop_var->val[0]| to... @>=
{ loop_var->val[0] = std::make_shared<int_value>(in_forward(flags) ? i++ : --i);
    // create a fresh index each time
  frame fr (pattern);
  fr.bind(loop_var);
  try {
    if (l==level::no_value)
      body->void_eval();
    else
    {@; body->eval();
       *(out_forward(flags) ? dst++ : --dst) = pop_value();
    }
   }
   @< Catch block for providing a trace-back of local variables @>
} // restore context upon destruction of |fr|

@ The loop over terms of a $K$-type polynomial is slightly different, and since
it handles values defined in the compilation unit \.{atlas-types}, we shall
include its header file. We implement the $4$ reversal variants; since terms are
sorted by ``height'' of a $K$-type, reversal at the source may make sense in
some cases. This reversal is implemented by using reverse iterators to control
the loop; the loop body itself is textually identical, though the type of |it|
differs between them. Something that is specific for the |Free_Abelian_light|
container class template used to implement |K_type_pol_value| is that its |size|
method only produces an upper bound for the actual number of (nonzero) terms
encountered during an iteration; therefore we must check after iteration whether
the expected number of items was copied into |result|. These checks are
independent of the traversal direction, so it is outside the conditional on
|in_forward(flags)| (the checks have to test |out_forward(flags)| instead).
These checks are conditional on |l!=level::no_value| however, since otherwise |result|
still holds a null (shared) pointer and the test would crash.

@h "atlas-types.h"
@< Perform a loop over the terms of a $K$-type polynomial @>=
{ shared_K_type_pol pol_val = get<K_type_pol_value>();
  size_t n=pol_val->val.size(); // an upper bound for the number of nonzero terms
  if (l!=level::no_value)
  { result = std::make_shared<row_value>(n);
    dst = out_forward(flags) ? result->val.begin() : result->val.end();
  }
  if (in_forward(flags))
  { const auto start=pol_val->val.begin();
    auto it = start; // need these in |catch| clause
    try
    {
      for ( ; it!=pol_val->val.end(); ++it)
        @< Loop body for iterating over terms of a $K$-type polynomial @>
    }
    @< Catch block for reporting iteration number within loop over
       terms in a $K$-type polynomial @>
  }
  else // not |in_forward(flags)|
  { const auto start=pol_val->val.rbegin();
    auto it = start; // need these in |catch| clause
    try
    {
      for (; it!=pol_val->val.rend(); ++it)
        @< Loop body for iterating over terms of a $K$-type polynomial @>
    }
    @< Catch block for reporting iteration number within loop over
       terms in a $K$-type polynomial @>
  }
  if (l!=level::no_value)
    @< Adjust |result| in case loop produced fewer items that predicted @>
}

@ The catch clause is similar to those for other types of loops,
but we use |std::distance| to compute the offset of |it| from its initial value
|start|.
@< Catch block for reporting iteration number within loop over terms in
   a $K$-type polynomial @>=
catch (error_base& e)
{
  std::ostringstream o;
  o << "During iteration " << std::distance(start,it) @|
    << " of the for-loop over KTypePol";
  e.trace(o.str());
  throw;
}


@~And here is the loop body  that is included twice identically.

@< Loop body for iterating over terms of a $K$-type polynomial @>=
{ loop_var->val[0] =
    std::make_shared<K_type_value>(pol_val->rf,it->first.copy());
  loop_var->val[1] = std::make_shared<split_int_value>(it->second);
  frame fr(pattern);
  fr.bind(loop_var);
  try {
    if (l==level::no_value)
      body->void_eval();
    else
    {@; body->eval();
      *(out_forward(flags) ? dst++ : --dst) = pop_value();
    }
   }
   @< Catch block for providing a trace-back of local variables @>
} // restore context upon destruction of |fr|

@ The kind of adjustments to |result| that are needed when fewer items than
expected were contributed are the same as when any kind of loop is interrupted
by an explicit |break_expr|: we must drop the slots in |result| that were not
filled, and in case of an output-reversed loop, we must first shift the actual
contributions to the beginning of the vector.

@< Adjust |result| in case loop produced fewer items that predicted @>=
{
  if (out_forward(flags))
  { if (dst!=result->val.end())
      result->val.resize(dst-result->val.begin());
  }
  else
  { if (dst!=result->val.begin())
  {
    dst = std::move(dst,result->val.end(), result->val.begin());
    result->val.resize(dst-result->val.begin());
  }}
}

@ The loop over terms of a virtual module is similar to that over those of a
$K$-type polynomial; however if the loop completes normally, there is no need to
check for an incomplete |result| since the exact number of iterations was
predicted (from |pol_val->val.size()|) before the loop started.

@< Perform a loop over the terms of a virtual module @>=
{ shared_virtual_module pol_val = get<virtual_module_value>();
  size_t n=pol_val->val.size();
  if (l!=level::no_value)
  { result = std::make_shared<row_value>(n);
    dst = out_forward(flags) ? result->val.begin() : result->val.end();
  }
  if (in_forward(flags))
  { const auto start=pol_val->val.cbegin();
    auto it = start; // need these in |catch| clause
    try {
    for ( ; it!=pol_val->val.cend(); ++it)
      @< Loop body for iterating over terms of a virtual module @>
    }
    @< Catch block for reporting iteration number within loop over
       terms in a virtual module @>
  }
  else
  { const auto start=pol_val->val.crbegin();
    auto it = start; // need these in |catch| clause
    try {
    for (; it!=pol_val->val.crend(); ++it)
      @< Loop body for iterating over terms of a virtual module @>
    }
    @< Catch block for reporting iteration number within loop over
       terms in a virtual module @>
  }
}

@ The catch clause is similar to those for other types of loops,
but we use |std::distance| to compute the offset of |it| from its initial value
|start|.
@< Catch block for reporting iteration number within loop over terms in
   a virtual module @>=
catch (error_base& e)
{
  std::ostringstream o;
  o << "During iteration " << std::distance(start,it) @|
    << " of the for-loop over ParamPol";
  e.trace(o.str());
  throw;
}


@~And here is the loop body  that is included twice identically.

@< Loop body for iterating over terms of a virtual module @>=
{ loop_var->val[0] =
    std::make_shared<module_parameter_value>(pol_val->rf,it->first);
  loop_var->val[1] = std::make_shared<split_int_value>(it->second);
  frame fr(pattern);
  fr.bind(loop_var);
  try {
    if (l==level::no_value)
      body->void_eval();
    else
    {@; body->eval();
      *(out_forward(flags) ? dst++ : --dst) = pop_value();
    }
   }
   @< Catch block for providing a trace-back of local variables @>
} // restore context upon destruction of |fr|



@*1 Counted loops.
%
Next we consider counted |for| loops. Like with other loops there are quite a
few variations. For efficiency we shall handle cases of omitted identifier and
(lower) |bound=0| (most likely due to an omitted \&{from} clause) specially.
We could do the all distinctions detectable at compile time through the
template argument |flags|, which would avoid any runtime tests, but give a lot
of cases. We choose to represent in |flags| all distinctions \emph{except} that
of an absent |bound| expression. The latter will be tested for presence when
initialising the lower bound; this gives a minute runtime cost when the bound
is present, but halves the number of template instances used.

@< Type def... @>=
template <unsigned flags>
struct counted_for_expression : public expression_base
{ id_type id; // may be $-1$, if |not has_frame(flags)|
  expression_ptr count, bound, body;
  // we allow |bound| (but not |count|) to hold |nullptr|
@)
  counted_for_expression
@| (id_type i, expression_ptr&& cnt, expression_ptr&& bnd,
    expression_ptr&& b)
@/: id(i), count(cnt.release()),bound(bnd.release()), body(b.release())
  @+{}
  virtual ~counted_for_expression() = default;
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

In order to substitute a counted loop without index for such a loop, we need to
assemble a call to the appropriate size-computing built-in function. The needed
|shared_builtin| values are held in variables whose definition will be given
later. It is specifically for this purpose that some of the actual built-in
functions needed here are declared (and in one case even defined in the first
place) as global rather than as local functions.

@< Set |loop| to a index-less counted |for| loop... @>=
{ expression_ptr call; const source_location &loc = f.in_part.loc;
  switch(which)
  {
  case subscr_base::row_entry: call.reset(new @|
      builtin_call(sizeof_row_builtin,std::move(in_expr),loc));
  break; case subscr_base::vector_entry: call.reset(new @|
      builtin_call(sizeof_vector_builtin,std::move(in_expr),loc));
  break; case subscr_base::ratvec_entry: call.reset(new @|
      builtin_call(sizeof_ratvec_builtin,std::move(in_expr),loc));
  break; case subscr_base::string_char: call.reset(new @|
      builtin_call(sizeof_string_builtin,std::move(in_expr),loc));
  break; case subscr_base::matrix_column: call.reset(new @|
      builtin_call(matrix_columns_builtin,std::move(in_expr),loc));
  break; case subscr_base::mod_poly_term: call.reset(new @|
      builtin_call(sizeof_parampol_builtin,std::move(in_expr),loc));
  break; default: assert(false);
  }

  if (f.flags.test(1)) // whether reversed assembly of return value
    loop.reset(new @| counted_for_expression<6>
      (-1,std::move(call),expression_ptr(nullptr),std::move(body)));
  else
    loop.reset(new @| counted_for_expression<4>
      (-1,std::move(call),expression_ptr(nullptr),std::move(body)));
}

@ Type-checking counted |for| loops is rather like that of other |for| loops,
but we must extend the context with the loop variable while processing the loop
body.

@< Cases for type-checking and converting... @>=
case cfor_expr:
{ const cfor_node& c=*e.cfor_variant;
  expression_ptr count_expr = convert_expr_strongly(c.count,fc,int_type);
  expression_ptr bound_expr = is_empty(c.bound)
    ? expression_ptr(nullptr)
    : convert_expr_strongly(c.bound,fc,int_type) ;
@)
  type body_type = type::bottom(fc);
  const conversion_record* conv=nullptr;
  type_expr tmp_tp = row_of_type.copy();
  bool do_export = tp.unify_specialise(tmp_tp);
  if (do_export)
    body_type = type::wrap(tmp_tp.component_type(),fc);
  else if (tp.is_void())
    body_type =tp.copy();
  else if ((conv=row_coercion(tp.bake(),tmp_tp))!=nullptr)
    body_type = type::wrap(tmp_tp,fc);
  else throw type_error(e,row_of_type.copy(),tp.bake());
@)
  if (c.flags[2]) // case of absent loop variable
  { layer bind(0,nullptr);  // no local variables for loop, but allow \&{break}
    expression_ptr body(convert_expr(c.body,body_type));
    if (do_export)
      tp.unify_to(body_type.wrap_row());
    if (not tp.is_void() and body_type.is_void() and not is_empty(c.body))
      body.reset(new voiding(std::move(body)));
  @/expression_ptr loop(make_counted_loop(c.flags.to_ulong(), @|
      c.id,std::move(count_expr),std::move(bound_expr),std::move(body)));
    return conv==nullptr
           ? std::move(loop) @|
           : expression_ptr(new conversion(*conv,std::move(loop)));
  }
  else // case of a present loop variable
  { layer bind(1,nullptr);
    bind.add(c.id,type::wrap(int_type,fc),true); // add |id| as constant
    expression_ptr body(convert_expr(c.body,body_type));
    if (do_export)
      tp.unify_to(body_type.wrap_row());
    if (not tp.is_void() and body_type.is_void() and not is_empty(c.body))
      body.reset(new voiding(std::move(body)));
  @/expression_ptr loop(make_counted_loop(c.flags.to_ulong(), @|
      c.id,std::move(count_expr),std::move(bound_expr),std::move(body)));
    return conv==nullptr
           ? std::move(loop) @|
           : expression_ptr(new conversion(*conv,std::move(loop)));
  }
}

@ Executing a loop is a simple variation of what we have seen before for
|while| and |for| loops over value components. Just like for the index while
looping over a row or vector object, there is a \Cpp\ counter controlling the
loop directly, and which the user loop variable shadows. However here the
counter is initialised and bounded by user program values which are |big_int|,
so the principled way to implement this would be to use a |big_int| counter;
however this probably would make all loops a bit slower, just in order to be
able to handle some very exotic cases, do we prefer to refuse values that will
not fit in a signed $64$-bit value, and use |long long int| (which is what
|long_val| returns) for the loop counter and bounds.

We distinguish a number of cases, some by template argument and others
dynamically, but always before entering the loop, in order to allow an optimal
adaptation to the task. We shall end up always using a \Cpp\ |while| loop to
implement the |for| loop. A special case that is optimised for is when the user
gave no name to the loop index: we can then omit introducing a |frame| for the
loop altogether. Also, since the syntax ensures that absence of a name for the
loop index implies absence of a lower bound expression (which would be unused
anyway) we can omit trying to evaluate |bound| here.

@< Function definitions @>=
template <unsigned flags>
void counted_for_expression<flags>::evaluate(level l) const
{ const auto n=(count->eval(),get<int_value>()->long_val());
  const auto lwb=(bound.get()==nullptr
          ? 0 : (bound->eval(),get<int_value>()->long_val()));
  long long int c = n<0 ? 0 : n; // no negative size result

  if (has_frame(flags)) // then loop uses index
  { long long int b=lwb;
    id_pat pattern(id);
    if (l==level::no_value)
      @< Perform counted loop that uses an index, without storing result,
         doing |c| iterations with lower bound |b| @>
    else
      @< Perform counted loop that uses an index, pushing result to
         |execution_stack|,
         doing |c| iterations with lower bound |b| @>
  }
  else if (l==level::no_value)
    @< Perform counted loop without index and no value @>
  else
    @< Perform counted loop without index producing a value @>
}

@ In cases of an anonymously counted loop, we choose to use a decreasing loop
counter internally. This simplifies the termination condition and might
therefore be marginally faster.

@< Perform counted loop without index and no value @>=
{ try {@;
    while (c-->0)
      body->void_eval();
  }
  catch (loop_break& err) @+
    {@; if (err.depth-- > 0)
          throw;
    }
  @< Catch block for reporting iteration number within
     anonymously counted loop @>

}

@ If a value is produced from each loop body evaluation, it is stored using the
|dst| iterator. In this case one must also take care, as before, that upon
catching a |break| executed inside the loop body, the returned row is shorter
than initially computed, and it may need to be shifted in the case of an
output-reversed loop.

@< Perform counted loop without index producing a value @>=
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
  @< Catch block for reporting iteration number within
     anonymously counted loop @>
@)
  push_value(std::move(result));
}

@ In case the loop body throws an error, no loop counter will be shown in the
back trace, but we can record an iteration count as was done for |while| loops.
@< Catch block for reporting iteration number within anonymously counted loop @>=
catch (error_base& e)
{
  std::ostringstream o;
  o << "During iteration " <<  n-1-c << " of the  counted for-loop";
  e.trace(o.str());
  throw;
}

@ For a counted loop using its index, we remain at the \Cpp\ level close to
the \.{axis} loop being implemented, but we transform the count |c| into an
upper bound, and then use |b| (for increasing loops) or |c| (for decreasing
loops) as our loop index. For decreasing loops we distinguish treat separately
the case where the lower bound is~$0$, the default value, since a test against a
constant~$0$ is probably a bit more efficient.

@< Perform counted loop that uses an index, without storing result,
   doing |c| iterations with lower bound |b| @>=
{ c+=b; // set to upper bound, exclusive
  if (c<b) // then additive overflow occurred
    throw runtime_error("Loop range would overflow signed long integer counter");
  try
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
  @< Catch block for reporting iteration number within indexed loop @>
}

@ The case where a result is accumulated just differs by adding the necessary
bits of stuff.

@< Perform counted loop that uses an index, pushing result to |execution_stack|,
   doing |c| iterations with lower bound |b| @>=
{ own_row result = std::make_shared<row_value>(c);
  c+=b; // set to upper bound, exclusive
  if (c<b) // then additive overflow occurred
    throw runtime_error("Loop range would overflow signed long integer counter");
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
  @< Catch block for reporting iteration number within indexed loop @>
  push_value(std::move(result));
}

@ Although we will not generate a local variable trace line for indexed loops,
we can report the name and value of the loop counter when reporting the loop
iteration number.
@< Catch block for reporting iteration number within indexed loop @>=
catch (error_base& e)
{
  std::ostringstream o;
  o << "During iteration ";
  if (in_forward(flags))
  { --b; // back up to actual loop index
    o << b-lwb << " (" << main_hash_table->name_of(id) << '=' << b @|
      << ") of the counted for-loop";
  }
  else
  {
    o << (b+n-1)-c << " (" << main_hash_table->name_of(id) << '=' << c @|
      << ") of the counted reversed for-loop";
  }
  e.trace(o.str());
  throw;
}

@* Casts and operator casts.
%
Casts are very simple to process. They do not need any |expression| type to
represent them, so type-checking is all there is to it. Nonetheless there is a
subtlety in the code below, in that |tp.specialise(c.dst_tp)| is called and used
to make (just) a separate call of |convert_expr| for the common case that |tp|
allows itself to be specialised to the cast type |c.dst_pt|, and therefore no
coercion will be need to be applied to the result of the cast. It would be
possible to leave this separation of paths out entirely, since the call to
|convert_expr_strongly| followed by |conform_types| (which starts trying to
|tp.specialise(c.dst_tp)|) can handle this no-coercion case as well. There is
however an advantage to our current approach, namely that if our cast is a
function body (as will often be the case), the specialisation of |tp| will
immediately set the return type variable of the function, which not only affects
ordinary type expected for the function body (which is our cast), but also of
any \&{return} expressions in that body.

Instead of a coercion, one might imagine that a discrepancy between the required
|tp| and the type |c.dst_pt| of the cast could be resolved by unification,
namely by assigning certain concrete types for certain polymorphic type
variables in the latter to obtain the former; that could be achieved by wrapping
|c.dst_tp| in a |type| variable, and using the alternative form of
|conform_types|, as we do for instance in the case of applied identifiers. The
reason we do not do so here is that we deliberately no allow a cast to specify a
polymorphic type (since allowing contexts that \emph{require} a polymorphic type
causes headaches we wish to avoid, and there is no obvious advantage of allowing
such casts); the situation allowing such unification following a cast therefore
simply cannot arise. A cast can of course require a type using abstract type
variables that are in scope, and the resulting type will then become polymorphic
if it is exported from that scope.

@< Cases for type-checking and converting... @>=
case cast_expr:
{ cast_node& c=*e.cast_variant;
  type_expr cast_tp = global_id_table->expand(c.dst_tp);
  if (tp.unify_specialise(cast_tp)) // see if we can do without coercion
    return convert_expr(c.exp,tp); // in which case use now specialised |tp|
  expression_ptr p = convert_expr_strongly(c.exp,fc,cast_tp);
    // otherwise use |cast_tp|
  return conform_types(c.dst_tp,tp,std::move(p),e);
}

@ Another kind of cast is the operator cast, which selects an operator or
overloaded function instance as it would for arguments of specified types, but
without giving actual arguments, so that the selected function itself can be
handled as a value. In order to do so we shall need the following function.

The overload table stores type information in a |func_type| value, whose two
components may contain polymorphic values numbered from~$0$. We want to
integrate these into a |type tp@;| that was initialised to a generic function
type with a possibly nonzero |floor()| value. We could use the subexpression
form of |unify_specialise| here with properly shifted copies of the |func_type|
components, but it seems a bit more straightforward to use |unify|.

@< Local function definitions @>=
inline bool functype_absorb (type& tp, const overload_data& entry)
{ auto te = type_expr::function @|
   (entry.f_tp().arg_type.copy(),
    entry.f_tp().result_type.copy());
  type model = type::wrap(te,0,tp.floor()); // proper shift facilitates |unify|
  return tp.unify(model);
}

@ Operator casts only access already existing values, which are looked up in the
global overload table. Since upon success we find a bare |shared_function|, we
must (as we did for~`\.\$') use the |capture_expression| class to serve as
wrapper that upon evaluation will return the value again.

We do two attempts to match the specified type |c_type| to an entry in the
|global_overload_table|. In the first we deduce from |c_type| a type as it would
be specified in a table entry for an exact match. Then as a second try, we try
to find a unique polymorphic overload that accepts |c_type| as it would in a
function call, and if one is found we use that. Although an exact match would
also be found in the second try, the first one allows specifying any overload,
even if the exact type of that overload also matches another, more general,
polymorphic overload. (In that case it is not possible to call the intended
overload through ordinary overload resolution, since no match to it would ever
be unique; we do want to be able to specify it in an operator cast, if only to
be able to apply \.{forget}, and thereby get rid of the useless overload.)

For the first try we use the method |overload_table::entry|, which in the call
below finds an entry at symbol |c->oper|, if there exists one, whose argument
type with any polymorphic variables shifted up by |fc| is |c_type|. That shift
is what we want, because at the point where the operator cast is written, any
polymorphic types introduced start at |fc|, while in the overload table there
are no fixed type variables and they start at~$0$. When |overload_table::entry|
returns a null pointer no exact match was found, and we pas to the second try,
which uses |overload_table::variants|.

@< Cases for type-checking and converting... @>=
case op_cast_expr:
{ const op_cast& c=e.op_cast_variant;
  const type_expr c_type = global_id_table->expand(c->arg_type);
  std::ostringstream o; // prepare name for value, and for error message
  o << main_hash_table->name_of(c->oper) << '@@' << c_type;
  expression_ptr result;
  type deduced_type = type::wrap(gen_func_type,fc);
  if (const auto* entry = global_overload_table->entry(c->oper,c_type,fc))
  {
    result.reset(new capture_expression(entry->value(),o.str()));
    bool success = functype_absorb(deduced_type,*entry);
    assert(success); ndebug_use(success);
  }
  else
  {
    type target = type::wrap(c_type,fc);
    const unsigned int target_deg=target.degree();
    const overload_data* prev_match=nullptr;
    if (@[auto* vars = global_overload_table->variants(c->oper)@;@])
      for (const auto& variant : *vars)
    @< See if |target| matches the argument type of a unique variant; if so,
       assign to |result| a |capture_expression| around its value, and to
       |deduced_type| the type with any substitutions to polymorphic type
       variables to match |target| applied to it @>
  }
  if (result==nullptr)
    throw expr_error(e,"No instance for "+o.str()+" found");
  return conform_types(deduced_type,tp,std::move(result),e);
}
break;

@ Contrary to ordinary casts, operator casts can specify a polymorphic
(argument) type, and this is indeed a natural way to select a polymorphic
variant; of course overloaded entries can have a polymorphic type as well.
Finding the right variant then resembles overloading resolution in function
calls, and can be done using the |type::matches| method; here too we insist of
having a unique variant match the specified parameter type. Our logic follows
that overload resolution closely, including the fact that a match is stored away
temporarily to see if a second matching variant exists, in which case we throw
an error for ambiguity instead of returning the initial match. Here too the call
of |matches| sets a possibly nonzero |shift_amount| by which the polymorphic
type variables from the overloaded binding are to be shifted upwards to steer
clear of any type variables (fixed or not) in |target|; the need to do so arises
even when |target| is monomorphic but involves type variables fixed in the
context. A difference is that here we need construct the entire function type
|deduced_type| obtained by unification, rather than just the return type. The
factory function |type_expr::function| will move from its argument types, so we
perform two separate substitutions to provide those. If nothing is found here,
we fall through this code, leading to a ``no instance found'' error.

@< See if |target| matches the argument type of a unique variant... @>=
{
  unsigned int op_deg = variant.poly_degree(), shift_amount;
  if (target.matches(variant.f_tp().arg_type,op_deg,shift_amount))
  {
    if (prev_match!=nullptr)
      @< Throw an error reporting an ambiguous match in operator cast @>
    @< Write to |o| the name of operator |c->oper| with the argument type of
       |variant|, to which the substitutions in |target.assign()| have been
       applied @>
    result.reset(new capture_expression(variant.value(),o.str()));
    auto te = type_expr::function @|
      (target.assign().substitution(variant.f_tp().arg_type,shift_amount)
      ,target.assign().substitution(variant.f_tp().result_type,shift_amount)
      );
    type model = type::wrap(te,target.floor());
    deduced_type.unify(model);
    prev_match = &variant;
  }
  target.clear(target_deg);
}

@ Similarly to what we do for ambiguous exact overload matches, we use the
|prev_match| pointer to build an error report.

@< Throw an error reporting an ambiguous match in operator cast @>=
{
  o.str(std::string()); // clear previous string prepared in |o|
  o << "Ambiguous argument in function call, specified type " << c_type
  @|<< " matches both " << prev_match->f_tp().arg_type
  @|<< " and " << variant.f_tp().arg_type;
  throw expr_error(e,o.str());
}

@ We store with the value returned a string that specifies both the polymorphic
argument type of the identified variant, and the substitutions that we made to
match the specified type.

@< Write to |o| the name of operator |c->oper| with the argument type of
   |variant|, to which the substitutions in |target.assign()| have been
   applied @>=
{
  o.str(std::string()); // clear previous string prepared in |o|
  o << main_hash_table->name_of(c->oper) << '@@' << variant.f_tp().arg_type;
  bool first=true;
  for (unsigned int i=0; i<op_deg; ++i)
    if (@[auto* p=target.assign().equivalent(i+shift_amount)@;@])
    {
      o << (first ? (first=false, '[') : ',')
        << *mk_type_variable(i) << '='
        << target.assign().substitution(*p);
    }
  if (not first)
    o << ']';
}


@* Assignments.
%
Syntactically there is hardly anything simpler than simple assignment
statements. However, semantically we distinguish assignments to local and to
global variables. Then there are ``component assignments'' which modify
composite values like row values by changing just one component; these too will
distinguish local and global versions. Finally, while not present in the initial
language design, a multiple assignment statement was added to the language that
can take apart tuple components, just as can be done in definitions of new
(global or local) variables. As these can mix local and global destinations, we
will introduce a separate expression type for these somewhat more expensive
assignments. We retain the following intermediate class for all assignment
statements except multiple assignments (the latter will derive directly from
|expression_base|), as it allows to avoid a bit of code duplication.

@< Type definitions @>=
struct assignment_expr : public expression_base
{ id_type lhs;
  expression_ptr rhs;
@)
  assignment_expr(id_type l,expression_ptr&& r)
   : lhs(l),rhs(r.release()) @+{}
  virtual ~assignment_expr() = default;
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
  virtual ~global_assignment() = default;
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
  virtual ~local_assignment() = default;
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
local names present in the destination pattern. This is done by having
(possibly empty) vectors for both types of destination, and a |BitMap| telling
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
    // type to describe a local binding
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
  virtual ~multiple_assignment() = default;
  virtual void print(std::ostream& out) const;
  virtual void evaluate(level l) const;
};

@ The constructor needs to be defined outside the class definition because it
uses |copy_id_pat|.
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

@ Assignment statements follows the same rules for locating the binding of their
left hand side as are used for applied identifiers. We first look in
|id_context| for a local binding of the identifier, and then maybe in
|global_id_table|. If found in either way, the right hand side is converted in a
type context given by the type of the variable found. After forming the proper
kind of assignment expression, we must as usual allow for a coercion to be
applied to the result of the assignment, if the externally required |tp|
demands this.

While in most successful cases the type of the variable governs the conversion
of the right hand side, the type of the right hand side may occasionally be more
specific than the previously known type of the variable (but any type that is
not a specialisation of the known type will cause the conversion to fail). For
instance, this specialisation happens when assigning a row of concrete type to a
variable initialised with an empty row. In those cases we call the |specialise|
method of |frame| or of |global_id_table| to make sure the type assumed by the
variable is recorded.

Since variables of |void| type are allowed (even if pretty useless), and can be
assigned to, we must take care to test for the necessity of a |voiding|. Indeed,
due to the voiding coercion a right hand side of any type will be accepted in an
assignment to a variable of void type, and since in this case the voiding is not
implied by the structure of the expression, we need to insert an explicit
|voiding| in such rare cases, to ensure that no actual (non void) value will be
computed and assigned.

After the call to |convert_expr|, we insert some code that tries to apply an
optimisation for certain built-in operations, to be detailed in the following
sections. Since identification of operations depends on types of their
arguments, this code needs to come after types have been checked, and must
operate on the converted expression |r| rather than on~|e|, even though this is
more difficult.

We could have written |conform_types(*id_t,tp,fc,std::move(assign),e)| as final
statement of the |if| branch, suggesting that the type |*id_t| of the variable
might be unified (reduced in polymorphic generality) to lead to the type |tp|
required in the context. But we leave out this possibility, for a similar reason
to why we omitted it for casts: the language will not allow assigning to
variables of polymorphic type (since this would give a strong polymorphic
context for the left hand side of the assignment). So |*id_t| is necessarily
monomorphic, and the case where unification would be useful cannot arise. By the
same token, we need not bother with unification of the results of other forms of
assignment, like multiple assignments or component assignments.

@< Cases for type-checking and converting... @>=
case ass_stat:
if (e.assign_variant->lhs.kind==0x1) // single identifier, do simple assign
{
  id_type lhs=e.assign_variant->lhs.name;
  const type* id_t; size_t depth,offset; bool is_const;
  const bool is_local = (id_t=layer::lookup(lhs,depth,offset,is_const))!=nullptr;
  if (not is_local and (id_t=global_id_table->type_of(lhs,is_const))==nullptr)
    report_undefined(lhs,e,"assignment");
@.Undefined identifier@>
  if (is_const)
    report_constant_modified(lhs,e,"assignment");
@.Name is constant @>
  assert(not id_t->is_polymorphic()); // polymorphic variables are made constant
@)
  const type_expr rhs_type = id_t->bake();
  expression_ptr r =
    convert_expr_strongly(e.assign_variant->rhs,fc,rhs_type);
  @< Check whether |r| refers to an |builtin_call| of a function with nonzero
     |hunger|, and with the identifier |lhs| as its corresponding argument;
     if so modify that argument, and possibly the application, accordingly @>

  if (id_t->is_void() and not is_empty(e.assign_variant->rhs))
    r.reset(new voiding(std::move(r)));

@)
  expression_ptr assign = is_local
  ? expression_ptr(new local_assignment(lhs,depth,offset,std::move(r)))
@/: expression_ptr(new global_assignment(lhs,std::move(r)));
  return conform_types(rhs_type,tp,std::move(assign),e);
}
else @< Generate and |return| a |multiple_assignment| @>

@ Here we seek top optimise assignments of the form $x:=\star x$ where $\star$
is a built-in operation that wants to change its argument in place. Normally the
potential access via the variable~$x$ makes the access to the operand of~$\star$
non-exclusive, so that a copy must be made for $\star$ to modify. However, if as
it is in the example the result of the operation is immediately assigned to~$x$,
one can defeat this unwanted and unneeded sharing by detaching the operand value
from~$x$ before invoking~$\star$, so that during a short time $x$ is not bound
to any value. In cases where this can be done safely, this is realised by
replacing the applied identifier expression for $x$ by a version that moves the
value out of~$x$, leaving it empty until the assignment.

This optimisation only applies to certain built-in operations, and none of them
are variadic, so we test whether |rhs| refers to a |builtin_call|, namely the
|overloaded_builtin_call| template instance with temple parameter |variadic| set
to |false|. The cases catered for are either single-argument operations, or
two-argument operations where the variable to be modified is expected to be a
specific one of the two operands; the attribute |hunger| of the |builtin_value|
template informs about this, with value $0$ stands for no desire to eat any
argument, values $1$ and $2$ for wanting to gobble up the left respectively
right argument out of two, and a value $3$ indicates a wish to transform a
unique argument.

@< Check whether |r| refers to an |builtin_call|... @>=
{
  auto* rhs = dynamic_cast<const builtin_call*>(r.get());
  if (rhs!=nullptr and rhs->f->hunger!=0)
  {
    const unsigned char h = rhs->f->hunger;
    if (h==3)
      @< See if |rhs->argument| is the variable |lhs|, and if so assign to~|r|
         a new call in which this variable gets emptied before the call @>
    else
      @< See if the argument of |rhs| indicated by $h\in\{1,2\}$ is the
         variable |lhs|, and if so change the argument list |rhs->argument| so
         that that argument gets evaluated last, emptying the variable while
         doing so @>
  }
}

@ Although we only need to change the |argument| field of the call pointed to by
|rhs|, the type of |rhs| is pointer to |const builtin_call|, which only gives
|const| access to |rhs->argument| (and in initialising |rhs| above, this
|const|ness was inherited from the |expression_ptr| type of |r|, which is a
unique pointer to |const expression_base|; the |dynamic_cast| used would fail to
work if we had omitted the |const|). So rather than assigning to |rhs->argument|
we build a new call, using mostly the values found at |rhs|, and assign it
to~|r|. In contrast to the original call, the new one sets the |pilfer| template
argument to |local_identifier| or |global_identifier| to |true|; at run time
this will cause evaluation to empty the variable, as is our goal here. An
alternative would be to defeat the |const|ness by assigning to
|const_cast<expression_ptr&>(c)|. Such casts are frowned upon, so we shall not
use one here. But in the next module the reconstruction efforts needed
would be even more tedious, and there we shall choose to cheat.

@< See if |rhs->argument| is the variable |lhs|, and if so...@>=
{ const expression_ptr& c = rhs->argument;
  const identifier* a = dynamic_cast<const identifier *>(c.get());
  if (a!=nullptr and a->code == lhs)
  { auto new_arg = is_local
      ? expression_ptr(new local_identifier<true>(a->code,depth,offset))
      : expression_ptr(new global_identifier<true>(a->code));
    r = expression_ptr(new builtin_call @|
        (rhs->f,rhs->name,std::move(new_arg),rhs->loc));
  }
}

@ When the built-in function that might benefit from changing a value in-place
takes two arguments, the value of |h| indicates which argument could so be
modified: left for $h=1$ and right for $h=2$. The change can only be applied if
there is a $2$-tuple of arguments, and the appropriate argument expression is an
applied identifier the coincides with the destination variable (since both
occurrences of the identifier arise in the same context, having the same name
ensures they will identify the same variable). When it applies, we want the
variable to be evaluated after the other argument, so that in the event where
that other argument also references the same variable, it will not find that
variable already detached from its value. Therefore the case $h=2$ is simpler
here, as the normal left-to-right evaluation order can be used.

In either case as before we first build a ``pilfering applied identifier''
expression |new_arg| (the pilfering being signalled by the |true| template
argument) with the otherwise same characteristics as the old applied identifier.
We still face the difficulty signalled above that we want to (move-)assign
|new_arg| to |c|, but that is a reference to |const expression_ptr|. Rather than
to construct anew all subexpressions that are to contain |new_arg|, we choose to
remove the |const| from the reference using a |const_cast|; doing so is safe
here, as no one else shares the expression under construction yet. For $h=2$,
that trick directly gives us what we want, but for $h=1$ we must also rebuild the
tuple expression as one with right-to-left evaluation semantics. In this case we
must not only defeat the |const|-ness when assigning to |rhs->argument|, but we
must also defeat the |const|-ness of the \emph{other} component of the old
$2$-tuple: we want to move from that $2$-tuple that is not going to get used,
to avoid having to make a deep copy, but moving from a reference to constant
will not work, and would instead try to copy construct a temporary, which is not
possible for |std::unique_ptr| instances. These uses of |const_cast| amount
to paying the price for wanting to rebuild the tree structure accessed by
|expressions_ptr| that was not designed to allow alteration after construction.

@< See if the argument of |rhs| indicated by $h\in\{1,2\}$ is...@>=
{
  const tuple_expression* p =
    dynamic_cast<const tuple_expression*>(rhs->argument.get());
  if (p!=nullptr and p->component.size()==2)
  {
    const expression_ptr& c = p->component[h-1];
    const identifier* a = dynamic_cast<const identifier *>(c.get());
    if (a!=nullptr and a->code == lhs)
    { auto new_arg = is_local
      ? expression_ptr(new local_identifier<true>(a->code,depth,offset))
      : expression_ptr(new global_identifier<true>(a->code));
      if (h==1) // hungry for first argument: evaluate arguments right-to-left
      {
        std::unique_ptr<tuple_expression_tmpl<true> > args @|
          (new tuple_expression_tmpl<true>(2));
        args->component[0] = std::move(new_arg);
        args->component[1] =
          // move other argument into new pair, after removing |const|-ness
          std::move(const_cast<expression_ptr&>(p->component[1]));
        const_cast<expression_ptr&>(rhs->argument) = std::move(args);
      }
      else
        const_cast<expression_ptr&>(c)= std::move(new_arg);
    }
  }
}

@ For traversing the left hand side pattern in a multiple assignment, we need
some semi-local variables, to be accessible from within the recursive function
but not renewed for each recursive call. The solution of passing around a
reference to a structure containing those variables is elegantly realised by
defining the traversal function as a recursive method |thread| of that structure
(the implicit reference |*this| is passed around unchanged). The goal of
|thread| is on one hand to determine the type expected for the right hand side
of the assignment, and to collect the target variables over which the components
of the right hand side value will be distributed; for the former an output
parameter |tp| is used, while the latter is stored in the fields of the
|threader| structure itself, which fields are |public|, so no methods need to be
declared to be able to recover them. After |thread| has completed, we can use
another method |refine| to specialise, if necessary, the types associated to the
variables in question. Those variables have been stored in the |locs| and
|globs| fields, and may get modified there. All in all, and somewhat
surprisingly, both methods return |void|.

@< Local class definitions @>=
struct threader
{ typedef containers::sl_list<multiple_assignment::local_dest> loc_list;
  typedef containers::sl_list<shared_share> glob_list;
@)
  const expr& e; // the multiple assignment expression we are working on
  loc_list locs; // local variables occurring, in order
  glob_list globs; // global variables occurring, in order
  BitMap is_global; // tells how locals and globals are interspersed
  containers::sl_list<std::pair<id_type,const_type_p> > assoc;
  // types found for them
@)
  threader (const expr& e) : e(e), locs(), globs(), is_global(), assoc() @+{}
  void thread (const id_pat& pat,type_expr& type); // recursively analyse |pat|
};

@ The left hand side pattern is traversed in post-order: when there is both an
identifier for the whole and a sub-list, the former is handled after the
latter. This simplifies testing of type compatibility in the destination
pattern, where the only possible error now is that a type for a ``parent''
identifier does not match the (possibly partly specified) tuple type
established by its children. The main information obtained by the traversal is
recorded in the output parameter |te|, which is preferred here over a return
value because its value is obtained by multiple calls of the |specialise|
method, and in addition possibly expansion of tabled type subexpressions.

@< Function definitions @>=
void threader::thread(const id_pat& pat,type_expr& te)
{ if ((pat.kind&0x4)!=0)
    @< Throw an error to signal forbidden qualifier \.! before |pat.name| @>
  if ((pat.kind&0x2)!=0) // first treat any sublist
  { te.expand().specialise(unknown_tuple(length(pat.sublist)));
    assert(te.raw_kind()==tuple_type); // succeeds due to successful |specialise|
    wtl_iterator t_it(te.tuple());
    for (auto it=pat.sublist.begin(); not pat.sublist.at_end(it); ++it,++t_it)
      thread(*it,*t_it);
  }
  if ((pat.kind&0x1)!=0)
  @< Look up type associated to |pat.name|, and after making some checks,
     record it in |te|, updating our fields |locs|, |globs|, |is_global| and
     |assoc| @>
}

@ While there are several things to do when processing each target, everything
is quite straightforward here. We need to check for the absence of repeated
identifiers, look up each identifier locally and maybe globally, refuse
assigning to identifiers that were marked as being constant, transfer the type
information from that lookup into |te| using a |specialise| call, and finally
recording the localisation of the identifiers in our various fields.

@< Look up type associated to |pat.name|, and after making some checks... @>=
{ id_type id = pat.name;
  const type* id_t; // will point to type of local or global |id|
  @< Check that |id| did not occur previously in this left hand side @>
  size_t i,j; bool is_const;
  const bool is_local = (id_t=layer::lookup(id,i,j,is_const))!=nullptr;
  if (not is_local and (id_t = global_id_table->type_of(id,is_const))==nullptr)
    report_undefined(id,e,"multiple assignment");
  if (is_const)
    report_constant_modified(id,e,"multiple assignment");
  assert(not id_t->is_polymorphic()); // polymorphic variables are made constant

  is_global.extend_capacity(not is_local); // push one bit onto the |BitMap|
@)
  if (not te.specialise(id_t->unwrap()))
  // incorporate type found for |id| into |te|
    @< Throw an error to signal type incompatibility for |id| @>
  assoc.push_back(std::make_pair(id,&te));
    // record pointer to |te| for later refinement of |id|
  if (is_local)
    locs.push_back(multiple_assignment::local_dest{i,j});
  else
    globs.push_back(global_id_table->address_of(id));
}

@ The error signalled here is quite silly: the user has qualified a target
identifier in the multiple assignment with ``\.!''. This should really be a
syntax error, as it is indeed in the case of simple assignments. But the fact
that a parser without conflicts can be generated for our grammar depends on the
fact that the pattern allowed after \.{set} is independent of whether \.=
or \.{:=} follows it; this is why we allowed these qualifiers to sneak and only
be detected here during context sensitive analysis. We are in fact paying here
for the use of the same keyword for two different purposes,

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
  @|<< " does no match pattern " << te;
  throw expr_error(e,o.str());
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
  type_expr lhs_te;
  threader thr(e);
  thr.thread(pat,lhs_te);
  type rhs_type = type::wrap(lhs_te,fc); // this can be polymorphic
  expression_ptr r = convert_expr(e.assign_variant->rhs,rhs_type);
  if (rhs_type.is_void() and not is_empty(e.assign_variant->rhs))
    r.reset(new voiding(std::move(r)));
  expression_ptr m_ass (
    new @| multiple_assignment
      (pat,std::move(r)
      ,thr.locs.undress(),thr.globs.undress(),std::move(thr.is_global)));
  return conform_types(rhs_type,tp,std::move(m_ass),e);
}

@*1 Component assignments.
%
The language we are implementing does not employ the notion of sub-object; in
other words if one sets $b=a[i]$ for some list, vector or matrix $a$, then $b$
will behave as a copy of the entry $a[i]$ rather than as an alias, so subsequent
assignment to $b$ will not affect~$a$ or vice versa. (This does no prevent us to
share storage between $b$ and $a$ initially, it just means the sharing should be
broken if $b$ or $a$ are modified; we practice copy-on-write.) This simplifies
the semantic model considerably; it makes no distinction between primitive and
composite values, avoids the distinction necessary for instance in Python or
Java between a (compound) value and the object that holds it. In \.{axis},
values that share the same memory behave exactly like values in separate memory
locations that happen to be equal.

However, if we want to allow creating composite values by sequentially setting
their components, we need to allow assignments of the form $a[i]:=c$ to achieve
the modifications, and we cannot view this as an operation on the ``subobject''
$a[i]$ of~$a$ (not having such a notion). The meaning of this is assignment will
be taken to be that of assigning a new value to all of $a$, which differs from
the original value only at index~$i$ (it will however be implemented more
efficiently if the storage of $a$ is not currently shared, as would usually be
the case, at least from the second such assignment to~$a$ onward). The
interpreter will treat such component assignments as a whole, using an
expression type with three components $a,i,c$, in which $a$ must be an
identifier. The type of this identifier may be one of several cases that allow
component assignments: any row type, vector, matrix, or an Atlas-specific
polynomial type. The class |component_assignment| below is a base class from
which specific classes for local and global assignments will be derived; its
only new data member is an expression |index| which at run time determines the
component that is to be changed (the aggregate name |lhs| and expression |rhs|
for the value to be assigned are members of its |assignment_expr| base class).
This class provides a method |assign| that will do the real work for the
|evaluate| methods of the derived classes, after those have located address of
the aggregate to be modified and the type of component assignment to apply.
Also, the class itself is templated over a Boolean |reversed| to allow for
reversed indexing. We similarly define a base class |field_assignment| for
assignments to a field of a value of some tuple type (whose usage requires
declaring named field selectors for that type).


@< Type definitions @>=

template <bool reversed>
struct component_assignment : public assignment_expr
{ expression_ptr index;
@)
  component_assignment
   (id_type a,expression_ptr&& i,expression_ptr&& r)
   : assignment_expr(a,std::move(r)), index(i.release()) @+{}
  virtual ~component_assignment() = default;

  virtual void print (std::ostream& out) const;
@)
  void assign(level l,shared_value& aggregate,subscr_base::sub_type kind) const;
};
@)
struct field_assignment : public assignment_expr
{ const unsigned position;
  id_type id;
@)
  field_assignment
   (id_type a,unsigned pos,id_type id,expression_ptr&& r)
   : assignment_expr(a,std::move(r)), position(pos), id(id) @+{}
  virtual ~field_assignment() = default;

  virtual void print (std::ostream& out) const;
@)
  void assign(level l,shared_value& tupple) const;
};

@ We define a variant of |component_assignment| called |component_transform|, in
which the new value of the component is computed with aid of its previous value;
this both avoids computing the indexing expression twice, and allows to try to
optimise the case where the component can be modified without duplication when
the transformation allows for modification in-place. The second point is
relevant only in the case of a component of a row-of type aggregate (rather than
a vector or matrix, whose components cannot be shared anyway), so we limit this
to cases corresponding to |kind==subscr_base::row_entry|, and there is no need
for a |kind| argument to |transform|.

@< Type definitions @>=

template <bool reversed>
struct component_transform : public assignment_expr
{ source_location loc; // as in |call_base|
  std::string name; // as in |overloaded_call|
  shared_builtin f; // as in |builtin_call|
  wrapper_function f_ptr; // shortcut, as in |builtin_call|
  expression_ptr index; // as in |component_assignment|
@)
  component_transform
   (id_type a, expression_ptr&& i,expression_ptr&& r,
    const shared_builtin& fun,
    const std::string& name, const source_location& loc)
   : assignment_expr(a,std::move(r))
   , loc(loc), name(fun->print_name), f(fun), f_ptr(fun->val)
   , index(i.release()) @+{}
  virtual ~component_transform() = default;

  virtual void print (std::ostream& out) const;
@)
  void transform(level l,shared_value& aggregate) const;
};
@)
struct field_transform : public assignment_expr
{ source_location loc; // as in |call_base|
  std::string name; // as in |overloaded_call|
  shared_builtin f; // as in |builtin_call|
  wrapper_function f_ptr; // shortcut, as in |builtin_call|
  const unsigned position;
  id_type id;
@)
  field_transform
   (id_type a,unsigned pos,id_type id,expression_ptr&& r,
    const shared_builtin& fun,
    const std::string& name, const source_location& loc)
   : assignment_expr(a,std::move(r))
   , loc(loc), name(fun->print_name), f(fun), f_ptr(fun->val)
   , position(pos), id(id) @+{}
  virtual ~field_transform() = default;

  virtual void print (std::ostream& out) const;
@)
  void transform(level l,shared_value& tupple) const;
};

@ Printing reassembles the subexpressions according to the input syntax, except
for field assignments which just print the position to be modified. As we shall
see below, the |right hand side| field can be null for the \&{transform}
structures, so we take care not to crash the program when this is the case.

@< Function def...@>=
template <bool reversed>
void component_assignment<reversed>::print (std::ostream& out) const
{ out << main_hash_table->name_of(lhs) << (reversed ? "~[" : "[")
      << *index << "]:=" << *rhs;
}
@)
template <bool reversed>
void component_transform<reversed>::print (std::ostream& out) const
{ out << main_hash_table->name_of(lhs) << (reversed ? "~[" : "[")
      << *index << "] " @| << name << ":= ";
  if (rhs==nullptr) out << "()"; @+ else out << *rhs;
}
@)
void field_assignment::print (std::ostream& out) const
{ out << main_hash_table->name_of(lhs) << '.' @|
      << main_hash_table->name_of(id)
      << '(' << this->position << ") := " @|
      << *rhs;
}
@)
void field_transform::print (std::ostream& out) const
{ out << main_hash_table->name_of(lhs) << '.' @|
      << main_hash_table->name_of(id)
      << '(' << this->position << ") " @|
      << name << ":= ";
  if (rhs==nullptr) out << "()"; @+ else out << *rhs;
}

@ For global assignments or transforms, we need to have non-|const| access the
location where the identifier is stored, whence the |shared_share| fields in the
definitions below.

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
  virtual void evaluate(eval_level l) const;
};
@)
template <bool reversed>
class global_component_transform
: public component_transform<reversed>
{ using base = component_transform<reversed>;
@)
  shared_share address;
public:
  global_component_transform
    (id_type a, expression_ptr&& i,expression_ptr&& r,
    const shared_builtin& fun,
    const std::string& name, const source_location& loc);
  virtual void evaluate(eval_level l) const;
};
@)
class global_field_assignment : public field_assignment
{ shared_share address;
public:
  global_field_assignment
    (id_type a, unsigned pos,id_type id,expression_ptr&& r);
  virtual void evaluate(eval_level l) const;
};
@)
class global_field_transform : public field_transform
{ shared_share address;
public:
  global_field_transform (id_type a, unsigned pos,id_type id,expression_ptr&& r,
    const shared_builtin& fun,
    const std::string& name, const source_location& loc);
  virtual void evaluate(eval_level l) const;
};

@ The constructor for |global_component_assignment| stores the address of the
aggregate object, the expression to be assigned, and the component kind. The
other cases don not store a component kind.

@< Function def... @>=
template <bool reversed>
global_component_assignment<reversed>::global_component_assignment @|
  (id_type a,expression_ptr&& i,expression_ptr&& r, subscr_base::sub_type k)
: base(a,std::move(i),std::move(r))
, kind(k),address(global_id_table->address_of(a)) @+{}
@)
template <bool reversed>
global_component_transform<reversed>::global_component_transform @|
    (id_type a, expression_ptr&& i,expression_ptr&& r,
    const shared_builtin& fun,
    const std::string& name, const source_location& loc)
: base(a,std::move(i),std::move(r),fun,name,loc)
, address(global_id_table->address_of(a)) @+{}
@)
global_field_assignment::global_field_assignment @|
  (id_type a, unsigned pos,id_type id,expression_ptr&& r)
: field_assignment(a,pos,id,std::move(r))
, address(global_id_table->address_of(a)) @+{}
@)
global_field_transform::global_field_transform @|
  (id_type a, unsigned pos,id_type id,expression_ptr&& r,
    const shared_builtin& fun,
    const std::string& name, const source_location& loc)
: field_transform(a,pos,id,std::move(r),fun,name,loc)
, address(global_id_table->address_of(a)) @+{}

@ It is in evaluation that component assignments differ most from ordinary ones.
The work is delegated to the |assign| or |transform| method of the base class
(to be defined below), which is given a reference to the |shared_value| pointer
holding the current value of the aggregate; it is this pointer that is in
principle modified. In the templated context of |global_component_assignment|
and |global_component_transform|, the base class must be explicitly mentioned
using the local type name |base| when calling its |assign| method. On the other
hand, while |global_field_assignment::evaluate| is also calling a method of that
name from its base class, no local type name is needed (nor is it defined) in
this case. Like when fetching the value of a global variable, we must be aware
of a possible undefined value in the variable.

@< Function def... @>=
template <bool reversed>
void global_component_assignment<reversed>::evaluate(eval_level l)
  const
{ if (address->get()==nullptr)
  { std::ostringstream o;
    o << "Assigning to component of uninitialized variable " @|
      << main_hash_table->name_of(this->lhs);
    throw runtime_error(o.str());
  }
  base::assign(l,*address,kind);
}
@)
template <bool reversed>
void global_component_transform<reversed>::evaluate(eval_level l)
  const
{ if (address->get()==nullptr)
  { std::ostringstream o;
    o << "Transforming component of uninitialized variable " @|
      << main_hash_table->name_of(this->lhs);
    throw runtime_error(o.str());
  }
  base::transform(l,*address);
}
@)
void global_field_assignment::evaluate(eval_level l) const
{ if (address->get()==nullptr)
  { std::ostringstream o;
    o << "Assigning to field of uninitialized variable " @|
      << main_hash_table->name_of(this->lhs);
    throw runtime_error(o.str());
  }
  assign(l,*address); // call method from base class (not called |base| here)
}
@)
void global_field_transform::evaluate(eval_level l) const
{ if (address->get()==nullptr)
  { std::ostringstream o;
    o << "Transforming field of uninitialized variable " @|
      << main_hash_table->name_of(this->lhs);
    throw runtime_error(o.str());
  }
  transform(l,*address); // call method from base class (not called |base| here)
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
  virtual void evaluate(eval_level l) const;
};
@)
template <bool reversed>
class local_component_transform : public component_transform<reversed>
{ using base = component_transform<reversed>;
@)
  size_t depth, offset;
public:
  local_component_transform @|
   (id_type arr, expression_ptr&& i,size_t d, size_t o, expression_ptr&& r,
    const shared_builtin& fun,
    const std::string& name, const source_location& loc);
  virtual void evaluate(eval_level l) const;
};
@)
class local_field_assignment : public field_assignment
{ size_t depth, offset;
public:
  local_field_assignment
    (id_type a, unsigned pos,id_type id,size_t d, size_t o, expression_ptr&& r);
  virtual void evaluate(eval_level l) const;
};
@)
class local_field_transform : public field_transform
{ size_t depth, offset;
public:
  local_field_transform
    (id_type a, unsigned pos,id_type id,size_t d, size_t o, expression_ptr&& r,
    const shared_builtin& fun,
    const std::string& name, const source_location& loc);
  virtual void evaluate(eval_level l) const;
};

@ The constructors for these structures are all quite straightforward, in
spite of their number of arguments.

@< Function def... @>=
template <bool reversed>
local_component_assignment<reversed>::local_component_assignment
 (id_type arr, expression_ptr&& i,size_t d, size_t o, expression_ptr&& r,
  subscr_base::sub_type k)
: base(arr,std::move(i),std::move(r)), kind(k), depth(d), offset(o) @+{}
@)
template <bool reversed>
local_component_transform<reversed>::local_component_transform @|
    (id_type a, expression_ptr&& i,size_t d, size_t o,expression_ptr&& r,
    const shared_builtin& fun,
    const std::string& name, const source_location& loc)
: base(a,std::move(i),std::move(r),fun,name,loc), depth(d), offset(o) @+{}
@)
local_field_assignment::local_field_assignment @|
  (id_type a, unsigned pos,id_type id,size_t d, size_t o, expression_ptr&& r)
: field_assignment(a,pos,id,std::move(r)), depth(d), offset(o) @+{}
@)
local_field_transform::local_field_transform @|
  (id_type a, unsigned pos,id_type id,size_t d, size_t o, expression_ptr&& r,
    const shared_builtin& fun,
    const std::string& name, const source_location& loc)
: field_transform(a,pos,id,std::move(r),fun,name,loc), depth(d), offset(o) @+{}

@ The |evaluate| methods locate the |shared_value| pointer of the aggregate,
then |assign| or |transform| does its job.

@< Function def... @>=
template <bool reversed>
void local_component_assignment<reversed>::evaluate(eval_level l)
  const
{@; base::assign (l,frame::current->elem(depth,offset),kind); }
@)
template <bool reversed>
void local_component_transform<reversed>::evaluate(eval_level l)
  const
{@; base::transform (l,frame::current->elem(depth,offset)); }
@)
void local_field_assignment::evaluate(eval_level l) const
{@; assign (l,frame::current->elem(depth,offset)); }
@)
void local_field_transform::evaluate(eval_level l) const
{@; transform (l,frame::current->elem(depth,offset)); }

@ The |assign| method, which will also be called for local component
assignments, starts by the common work of evaluating the (component) value to be
assigned. For actually changing the aggregate, we must distinguish cases
according to the kind of component assignment at hand. Assignments to components
of rational vectors and of strings will be forbidden, see
module@#comp_ass_type_check@>. The evaluation of the aggregate index is done
inside this case distinction, because possible expansion of a tuple index value
depends on~|kind|. Finally we shall make sure we hold a unique copy of the
aggregate; since |uniquify|, which does this operation, needs to know the type
of the aggregate in a template argument, its call has to be done inside each
branch.

For assignments to a field of a tuple value, no case distinction is necessary,
and the code is quite simple.

@< Function def... @>=
template <bool reversed>
void component_assignment<reversed>::assign
  (level lev,shared_value& aggregate, subscr_base::sub_type kind) const
{ rhs->eval();
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
    case subscr_base::K_type_poly_term:
  @/@< Replace coefficient at |index| in $K$-type polynomial |loc|
       by value on stack @>
  @+break;
    case subscr_base::mod_poly_term:
  @/@< Replace coefficient at |index| in virtual module |loc|
       by value on stack @>
  @+break;
  default: {} // remaining cases are eliminated in type analysis
  }
}
@)
void field_assignment::assign (level lev,shared_value& tupple) const
{ shared_value& field=uniquify<tuple_value>(tupple)->val[position];
  rhs->eval();
  push_expanded(lev,field=pop_value());
}

@ A |row_value| component assignment is the simplest kind. The variable
|aggregate| holds a generic |shared_value|, known to refer to a |row_value|.
Since we need to access the vector of shared pointers, we use |uniquify| to
ensure unshared and unconstrained access to its |val| field. After a bound check
we then replace a component by the stack-top value. Afterwards, depending on
|lev|, we may put back the stack-top value as result of the component assignment,
possibly expanding a tuple in the process.

@< Replace component at |index| in row |loc|... @>=
{ auto i = (index->eval(),get<int_value>()->long_val());
  auto& a = uniquify<row_value>(aggregate)->val;
  size_t n=a.size();
  if (static_cast<size_t>(i) >= n)
    throw runtime_error(range_mess(i,a.size(),this,"component assignment"));
  auto& ai = a[reversed ? n-1-i : i];
  push_expanded(lev,ai = pop_value()); // assign component and yield that value
}

@ The |transform| method of |component_transform| is similar to the |assign|
methods above. However, instead of assigning the evaluation of |rhs| into the
aggregate, it combines it with the previous value of the destination component
using |f_ptr| (which is a |shared_builtin|). Since the first operand of |f_ptr|,
which is the old value of the row component |ai|, is moved out of the row to the
stack, we make sure to evaluate the second operand before it, so that during its
evaluation the row |a| is still intact. This requires temporarily moving that
second argument off the stack before moving it back on. (The expression for that
operand is called |rhs|, in the |assignment| class we inherit from.) The code
takes into account the possibility of an absent second argument, indicated by
the condition |rhs=nullptr|; this is because an optimisation may have replaced a
call of an operator with two arguments by a call of a function with only one
argument, as in $v[i]\mathrel+:=1$ where the addition gets replaced by a call of
|succ|. In that case |f_ptr| will point to the replacement function and a null
pointer is substituted for~|rhs|.

The method |field_transform::transform| is similar but simpler, and shares the
part calling |f_ptr|.

@< Function def... @>=
template <bool reversed>
void component_transform<reversed>::transform
  (level lev,shared_value& aggregate) const
{ auto op2 = (rhs==nullptr ? nullptr : (rhs->eval(),pop_value()));
   // put aside additional operand
  auto i = (index->eval(),get<int_value>()->long_val());
  auto& a = uniquify<row_value>(aggregate)->val;
  size_t n=a.size();
  if (static_cast<size_t>(i) >= n)
    throw runtime_error(range_mess(i,a.size(),this,"component assignment"));
  auto& ai = a[reversed ? n-1-i : i];
  push_value(std::move(ai)); // move-push component before transformation
  if (op2!=nullptr)
    push_value(std::move(op2)); // and possibly additional argument
  @< Call |*f_ptr| to produce a single value, taking measures for back tracing @>
  push_expanded(lev,ai = pop_value()); // assign component and yield that value
}
@)
void field_transform::transform (level lev,shared_value& tupple) const
{ auto op2 = (rhs==nullptr ? nullptr : (rhs->eval(),pop_value()));
   // put aside additional operand
  shared_value& field=uniquify<tuple_value>(tupple)->val[position];
  push_value(std::move(field)); // move-push field before transformation
  if (op2!=nullptr)
    push_value(std::move(op2)); // and possibly additional argument
  @< Call |*f_ptr| to produce a single value, taking measures for back tracing @>
  push_expanded(lev,field=pop_value());
}

@ This code is similar to that of |built_in::evaluate|, except that we know we
have exactly two arguments, and that they are already placed on the
|execution_stack|.

@< Call |*f_ptr| to produce a single value, taking measures for back tracing @>=
{ std::string arg_string;
  if (verbosity!=0) // record argument(s) as string
  { std::ostringstream o;
    if (rhs==nullptr)
      o << '(' << *execution_stack[execution_stack.size()-1] << ')';
    else
    { const auto* p = &execution_stack[execution_stack.size()-2];
      o << '(' << *p[0] << ',' << *p[1] << ')';
    }
    arg_string = o.str();
  }
@)
  try
  {@; (*f_ptr)(eval_level::single_value); } // call the built-in function
  catch (error_base& e)
  {@; extend_message(e,name,loc,f,arg_string);
    throw;
  }
  catch (const std::exception& e)
  { runtime_error new_error(e.what());
    extend_message(new_error,name,loc,f,arg_string);
    throw new_error;
  }
}

@ We complete our definition with the non-row component assignments, starting
with the case of |vec_value| entry assignments. Here the type of the aggregate
object is vector, and the value assigned always an integer. The latter certainly
needs no expansion, so we either leave it on the stack, or remove it if the
value of the component assignment expression is not used.

@< Replace entry at |index| in vector |loc|... @>=
{ auto i=(index->eval(),get<int_value>()->long_val());
  auto& v = uniquify<vector_value>(aggregate)->val;
  size_t n=v.size();
  if (static_cast<size_t>(i) >= n)
    throw runtime_error(range_mess(i,v.size(),this,"component assignment"));
  v[reversed ? n-1-i : i] =
  // assign |int| extracted from stack top without popping
    force<int_value>(execution_stack.back().get())->int_val();
  if (lev==level::no_value)
    execution_stack.pop_back(); // pop it anyway if result not needed
}

@ For matrix entry assignments, the value |index| must be split into a pair of
indices, and there are two bound checks.

@< Replace entry at |index| in matrix |loc|... @>=
{ index->multi_eval();
  auto j=get<int_value>()->long_val();
  auto i=get<int_value>()->long_val();
@/
  auto& m = uniquify<matrix_value>(aggregate)->val;
  size_t k=m.n_rows(),l=m.n_columns();
  if (static_cast<size_t>(i) >= k)
    throw runtime_error@|
      ("initial "+range_mess(i,m.n_rows(),this,"matrix entry assignment"));
  if (static_cast<size_t>(j) >= l)
    throw runtime_error@|
      ("final "+range_mess(j,m.n_columns(),this,"matrix entry assignment"));
  m(reversed ? k-1-i : i,reversed ? l-1-j : j)=
    // assign |int| from un-popped top
    force<int_value>(execution_stack.back().get())->int_val();
  if (lev==level::no_value)
    execution_stack.pop_back(); // pop it anyway if result not needed
}

@ A matrix column assignment is like that of a vector entry, but with in
addition a necessary test for matching column length of the vector assigned.

@< Replace column at |index| in matrix |loc|... @>=
{ auto j=(index->eval(),get<int_value>()->long_val());
  auto& m = uniquify<matrix_value>(aggregate)->val;
@/const int_Vector& v=force<vector_value>(execution_stack.back().get())->val;
    // don't pop
  size_t l=m.n_columns();
  if (static_cast<size_t>(j) >= l)
    throw runtime_error(
      range_mess(j,m.n_columns(),this,"matrix column assignment"));
  if (v.size()!=m.n_rows())
    throw runtime_error
      (std::string("Cannot replace column of size ")+str(m.n_rows())+
       " by one of size "+str(v.size()));
  m.set_column(reversed ? l-j-1 : j,v);
    // copy value of |int_Vector| into the matrix
  if (lev==level::no_value)
    execution_stack.pop_back(); // pop the vector if result not needed
}

@ For |K_type_pol_value| coefficient assignments the type of the aggregate
object is $K$-type polynomial, and the value assigned a split integer. The
latter certainly needs no expansion, so we either leave it on the stack, or
remove it if the value of the component assignment expression is not used.

@< Replace coefficient at |index| in $K$-type polynomial |loc|... @>=
{ index->eval();
  auto t = get<K_type_value>();
  auto* pol = uniquify<K_type_pol_value>(aggregate);
  const auto& top = force<split_int_value>(execution_stack.back().get());
  pol->assign_coef(*t,top->val);
  if (lev==level::no_value)
    execution_stack.pop_back(); // pop the vector if result not needed
}

@ For |virtual_module_value| coefficient assignments the type of the aggregate
object is ``virtual module'', and the value assigned a split integer, again
needing no expansion. The main differences between this module and the previous
one are hidden in the respective |assign_coef| methods.

@< Replace coefficient at |index| in virtual module |loc|... @>=
{ index->eval();
  auto t = get<module_parameter_value>();
  auto* pol = uniquify<virtual_module_value>(aggregate);
  const auto& top = force<split_int_value>(execution_stack.back().get());
  pol->assign_coef(*t,top->val);
  if (lev==level::no_value)
    execution_stack.pop_back(); // pop the vector if result not needed
}

@*2 Type-checking component assignments.
%
Type-checking and converting component assignment statements follows the
same lines as that of ordinary assignment statements, but must also
distinguish different aggregate types.

@:comp_ass_type_check@>

@< Cases for type-checking and converting... @>=
case comp_ass_stat:
{ id_type aggr=e.comp_assign_variant->aggr;
  const expr& index=e.comp_assign_variant->index;
  const expr& rhs=e.comp_assign_variant->rhs;
@/const type* aggr_tp; size_t d,o; bool is_const;
  bool is_local = (aggr_tp=layer::lookup(aggr,d,o,is_const))!=nullptr;
  if (not is_local and (aggr_tp=global_id_table->type_of(aggr,is_const))==nullptr)
    report_undefined(aggr,e,"component assignment");
@.Undefined identifier@>
  if (is_const)
    report_constant_modified(aggr,e,"component assignment");
@.Name is constant @>
  assert(not aggr_tp->is_polymorphic());
  // polymorphic variables are made constant
@)
  type ind_tp = type::bottom(fc);
  expression_ptr i = convert_expr(index,ind_tp);
@/type_expr comp_te;
  subscr_base::sub_type kind =
    subscr_base::index_kind(aggr_tp->expanded(),ind_tp.bake(),comp_te);
  if (not subscr_base::assignable(kind))
  { std::ostringstream o;
    o << "Cannot subscript value of type " << *aggr_tp @|
      << " with index of type " << ind_tp << " in assignment";
    throw expr_error(e,o.str());
  }
  type comp_tp = type::wrap(comp_te,fc);
  expression_ptr r = convert_expr(rhs,comp_tp);
  if (comp_tp.is_void() and not is_empty(rhs))
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
  return conform_types(comp_tp,tp,std::move(p),e);
}

@ And here's the analogous expression analysis for field assignments. We must
find the field selector in the overload table with the exact (tuple) type, and
it must be bound to a |projector_value| (so some user defined function that
selects the proper field will not work here).

@< Cases for type-checking and converting... @>=
case field_ass_stat:
{ id_type tuple=e.field_assign_variant->aggr;
  id_type selector =e.field_assign_variant->selector;
  const expr& rhs=e.field_assign_variant->rhs;
@/const type* tuple_tp; size_t d,o; bool is_const;
  bool is_local = (tuple_tp=layer::lookup(tuple,d,o,is_const))!=nullptr;
  if (not is_local
      and (tuple_tp=global_id_table->type_of(tuple,is_const))==nullptr)
    report_undefined(tuple,e,"field assignment");
@.Undefined identifier@>
  if (is_const)
    report_constant_modified(tuple,e,"field assignment");
@.Name is constant @>
  assert(not tuple_tp->is_polymorphic());
  // polymorphic variables are made constant
@)
  unsigned pos; type_expr component;
  @< Look up a field of |*tuple_tp| named |selector|, and if found assign
     its position to |pos| and set |component| to the corresponding
     component of |*tuple_tp|; on failure |throw expr_error| @>
  expression_ptr r = convert_expr_strongly(rhs,fc,component);
  expression_ptr p;
  if (is_local)
    p.reset(new local_field_assignment(tuple,pos,selector,d,o,std::move(r)));
  else
    p.reset(new global_field_assignment(tuple,pos,selector,std::move(r)));
  return conform_types(component,tp,std::move(p),e);
}

@ If either no function at all doing the requested projection is found, or if
the function found is not a projector, then we signal failure. Moving |pos|
places forward in the linked list can be done by calling |std::next| after
converting the raw node pointer |tuple_tp->tuple()| to a weak type list iterator.

@< Look up a field of |*tuple_tp| named |selector|... @>=
{ std::ostringstream o;
  if (tuple_tp->top_kind()!=tuple_type)
    throw expr_error
      (e,"Field assignment with variable of non tuple type");
  auto candidates = type_expr::matching_bindings(*tuple_tp);
  if (candidates.empty())
    @< Report that field assignments for tuple type |*tuple_tp| require a
       type definition with field names @>
  for (auto it=candidates.begin(); not candidates.at_end(it); ) // no |++it|
  { unsigned int i;
    for (i=0; i<(*it)->fields.size(); ++i)
      if ((*it)->fields[i]==selector)
      {@; pos=i;
        break;
      }
    if (i==(*it)->fields.size())
      candidates.erase(it);
    else
      ++it;
  }
  if (candidates.empty())
    @< Report that |*tuple_tp| has no defined field named |selector| @>
  else if (not candidates.singleton())
    @< Report selecting |selector| from |*tuple_tp| is ambiguous @>
@)
  type_expr tup_exp = tuple_tp->expanded(); // ensure tabled type is expanded
  component = std::move(*std::next(wtl_iterator(tup_exp.tuple()),pos));
}

@ We try to give error messages that identify clearly what has gone wrong.

@< Report that field assignments for tuple type |*tuple_tp| require a
   type definition with field names @>=
{ o << "Type " << *tuple_tp
    << " of variable in field assignment has no associated field names";
  throw expr_error(e,o.str());
}

@ When there are field selectors for the tuple type, but the selector matches
none of the candidates, we blame the selector rather than the type.

@< Report that |*tuple_tp| has no defined field named |selector| @>=
{ o << "Type " << *tuple_tp @|
    << " of variable in field assignment has no field '" @|
    << main_hash_table->name_of(selector) << '\'';
  throw expr_error(e,o.str());
}

@ In the rare case that there is more than one tuple type (constructor) that
both matches the type of the variable being selected from and have a field
called |selector|, we just say the situation is ambiguous.

@< Report selecting |selector| from |*tuple_tp| is ambiguous @>=
{ o << "Type " << *tuple_tp @|
    << " of variable matches more than one definition with field name '" @|
    << main_hash_table->name_of(selector) << '\'';
  throw expr_error(e,o.str());
}

@ Although type-checking component and field transform statements is similar to
their assignment counterparts, converting them is quite complicated. Since we
come here before type checking is done, we have to handle all expressions of the
syntactic form in question (operation-assigning to a field selection from an
identifier expression), whether or not that gives any occasion to invoke
in-place modification. The conditions for that optimisation to be possible are
that the transformation is to be performed by a built-in operator, and that its
first operand and result types coincide with that of the selected field. The
identification of the operator in fact uses that type for the first operand, so
that part of the condition is likely to be satisfied if the operator can be
found at all, but the type condition still requires the absence of implicit
conversions; if there were any, that would frustrate any in-place operation
anyway. The condition of being built-in (which implies being found in the global
overload table, as we cannot determine the actual value in local bindings) is
rather restrictive, and we would have like to not impose it, but the in-place
field transformation semantics we want to apply mean that the tuple has a hole
in it at the moment the operation is applied, and for user defined functions we
cannot ensure that they cannot notice this circumstance.

Concretely, we type check the call of the operator ignoring field-transforming
context (but supplying the required return type) and then test whether we can
use a |field_transform|; if we can, we build it using pieces of the converted
expression, which includes the identity of the built-in operation that was
found. In the contrary case, we reassemble the pieces together differently, to
form an ordinary |field_assignment| instead, with a value produced by an
ordinary function call of the operator.

@< Cases for type-checking and converting... @>=
case field_trans_stat:
{ expr& lhs = e.comp_trans_variant->dest;
  id_type op = e.comp_trans_variant->op;
  expr& rhs=e.comp_trans_variant->arg;
@/
  assert(lhs.kind == function_call);
  app dot = lhs.call_variant; // call of the field selector
  assert(dot->arg.kind==applied_identifier);
  id_type tuple = dot->arg.identifier_variant; // name of tuple selected from
@/assert(dot->fun.kind==applied_identifier);
  id_type selector = dot->fun.identifier_variant; // field name
@)const type* tuple_tp; size_t d,o; bool is_const;
  bool is_local = (tuple_tp=layer::lookup(tuple,d,o,is_const))!=nullptr;
  if (not is_local and
      (tuple_tp=global_id_table->type_of(tuple,is_const))==nullptr)
    report_undefined(tuple,e,"field transform");
@.Undefined identifier@>
  if (is_const)
    report_constant_modified(tuple,e,"field transform");
@.Name is constant @>
  assert(not tuple_tp->is_polymorphic());
  // polymorphic variables are made constant
@)
  unsigned pos; type_expr component;
  @< Look up a field of |*tuple_tp| named |selector|... @>
  expression_ptr call;
  @< Assign to |call| the |convert_expr| of the application of |op| to
     an argument pair formed of |lhs|... @>
  @< Construct, from |tuple|, |pos| and |*call|... @>
}

@ Here we build the application |appl| of the symbol |op| mentioned in the title
as an |expr| structure, and pass it to |convert_expr| that most probably will
call |resolve_overload|. This mainly serves to find the relevant instance of
|op|, but the converted expression |call| or part of it will also be used. The
reason |appl| is |static| is for correct error reporting, as explained below.

@< Assign to |call| the |convert_expr| of the application of |op| to an
   argument pair formed of |lhs| and |rhs|, converted to type |component| @>=
{
  static expr_ptr appl; expr_ptr saved_appl;
  @< Set |appl| to the application of |op| to |lhs| and |rhs|... @>
  call = convert_expr_strongly(*appl,fc,component);
  @< Restore initial state of |*e.comp_trans_variant| and of |appl| @>
}

@ The elements of the argument pair we construct are stolen (i.e., moved) from
the (component transformation) expression |e| we are processing. This is
necessary because we cannot copy |expr| values; after converting |appl| we shall
move these parts back into |lhs| and |rhs| so that |e| is intact on successful
return, and can be used for later error messages. For similar reasons we must
move aside the previous contents of the |static| variable |appl| so that it
won't be destructed in the assignment; in fact it could well hold an expression
that contains |e| as subexpression, in which case its destruction would have
catastrophic consequences.

Moving the two arguments is done in two steps, since in a |tuple_expression| the
components are accessed by an |expr_ptr|, which cannot point to an |expr| value
that is contained in a larger structure, as is the case for |lhs| and |rhs|.

@< Set |appl| to the application of |op| to |lhs| and |rhs|... @>=
{ saved_appl = std::move(appl);
  expr_ptr arg1(new expr(std::move(lhs)));
    // move top level data into isolated |expr|
  expr_ptr arg2(new expr(std::move(rhs))); // likewise
@/appl =
    internal_binary_call(op,std::move(arg1),std::move(arg2),
                         e.loc,e.comp_trans_variant->op_loc);
}

@ We have moved parts from the expression |*e.comp_trans_variant| into the one
accessed by |appl|, but our caller must see the whole expression |e| intact, in
case it is a subexpression of a larger one for which an |expr_error| (or derived
instance) will be thrown. So upon successful conversion, we must dig into the
parts of |appl| and move them back where they came from. Then (and only then) we
must also restore the value of the |static| variable |appl| itself, which was
saved in a local variable |save|; this will also clean up the parts of the
|appl| expression that we built ourselves (rather than moved). These manoeuvres
could have been avoided if |expr| were copy constructible, but writing a
(recursive, deep) copy constructor would be even more work.

The reason |appl| is a static variable is that in case the |convert_expr|
call above should throw an error, a reference to |*appl| will be stored in the
|expr_error| object, and it will be caught only after stack unwinding has
destroyed all local variables of our (recursive) function |convert_expr|; if
|appl| were such a local variable the mentioned reference would become a
dangling one. As it is, the static variable will keep alive the |expr| it points
to, even after reporting it, which is useless. However this memory wastage is a
one-off (a new error thrown from the same place will replace and clean up the
|expr| value), so we don't care.

@< Restore initial state of |*e.comp_trans_variant| and of |appl| @>=
{
  auto* args = appl->call_variant->arg.sublist;
@/lhs = std::move(args->contents);
  rhs = std::move(args->next->contents);
@/appl = std::move(saved_appl); // restore, so our caller will not notice
}

@ The code below illustrates one way to solve the coding problem of avoiding
multiple identical |else| clauses in a situation where the positive option (here
that of using a |field_transform| rather than a |field_assignment|) is dependent
on the conjunction of several conditions. Here the positive option requires that
the |call| is an application of a built-in operator (necessarily found in the
global overload table), and that its argument is either a pair with a projector
call (i.e., a field selection) as first element, or just a projector call.
Testing this involves introducing several intermediate values |c|, |arg|, and
|tup| that depend on previous ones, and may be used (only) in the branch for
this option. The straightforward approach of using nested |if| expressions to
successively test the conditions would introduce several identical |else|
clauses. Those clauses could be fused into a single one using |goto|, but that
solution is distinctly ugly. The solution adopted here is to have a sequence of
tests conditionally setting the intermediate variables, which are pointers left
null otherwise, and then leave the initial |if| expression and follow it by a
second one that tests the final value. (In fact the test is a dynamic cast of
|arg|, which can fail either if |arg==nullptr| or if the dynamic cast finds a
wrong pointer; the cast pointer itself is not needed in the sequel.)

The complication that |arg| can be defined in two was is due to the optimisation
that, although we used |build_binary_call| which involves two arguments, on
optimisation during conversion may have replaced it by a call with a single
argument, like replacing |x.a+1| by |succ(x.a)|. If this happens and the
|field_transform| branch applies, it passes a null pointer in place of the |rhe|
right hand expression to the |field_transform| constructor, which the evaluation
functions will detect to avoid actually evaluating a second argument.

The actual work to be done is quite straightforward, since all the pieces from
which we want to construct either a |field_transform| or a |field_assignment|
have already been converted. In the former case we use just the converted second
argument |rhe| of the |call| of |op|, while in the latter case we use |call| as
a whole.

@< Construct, from |tuple|, |pos| and |*call|, either a structure derived from
   |field_transform| or one derived from |field_assignment|, and |return|
   the result of passing it through |conform_types| @>=
{ const tuple_expression* tup = nullptr;
  const expression_base* arg = nullptr;
  auto* c = dynamic_cast<const builtin_call*>(call.get());
  if (c!=nullptr)
  {
    tup = dynamic_cast<const tuple_expression *>(arg = c->argument.get());
    if (tup != nullptr)
      arg = tup->component[0].get();
  }
@)
  expression_ptr re;
  if (dynamic_cast<const projector_call*>(arg)!=nullptr)
  {
    expression_ptr nil(nullptr);
    auto& rhe =
      tup==nullptr ? nil : const_cast<expression_ptr&>(tup->component[1]);
    if (is_local)
      re.reset (new local_field_transform@|
        (tuple,pos,selector,d,o,std::move(rhe),c->f,c->name,e.loc));
    else
      re.reset (new global_field_transform@|
        (tuple,pos,selector,std::move(rhe),c->f,c->name,e.loc));
  }
  else
  { if (is_local)
      re.reset(new local_field_assignment
        (tuple,pos,selector,d,o,std::move(call)));
    else
      re.reset(new global_field_assignment(tuple,selector,pos,std::move(call)));
  }
  return conform_types(component,tp,std::move(re),e);
}

@ All that was done for the case of field transformations in a tuple must also
be done for transformations of a component in a row, with some extra
complications. Again we handle all expressions of the syntactic form
(operation-assigning to a subscripted name), and decide only after the
identification of the operation whether we can actually generate a
|component_transform| or whether we expand to an ordinary
|component_assignment|. The additional complications here are that there is an
index expression whose double evaluation must be avoided, which may require
additional rewriting of the syntax tree when we revert to a
|component_assignment|, and the presence of a |reversed| attribute which makes
the code for generating |component_transform| or |component_assignment| more
repetitive.

The conditions that need to be satisfied are those of a field transformation,
plus the fact that name we are subscripting must have row-of type: the
|component_transform| class was not designed to handle any other |kind| of
subscription, which indeed do not appear to be able to benefit from in-place
transformation. As before we shall type check a call of the operator as if
|component_transform| is involved, and then from the result try to find out
whether our conditions are satisfied. If they are, then we proceed much like in
the field transformation case. In the other case we reassemble the pieces into a
|component_assignment|.

By consistent choice of variable names, we can reuse here the pieces of code
that assemble the pieces of our original |comp_transform_node| into a function
application, and then later restore everything to its initial state. The two
pieces are further apart here, because as we shall see the function call
expression |appl| may need to be converted a second time in a different context,
so it is left intact until that is behind us. This also means that we cannot
have and |return| expressions before the restoring is done (as that would skip
their execution), and instead we just have a variable |result| that is set
differently in different cases.

@< Cases for type-checking and converting... @>=
case comp_trans_stat:
{ expr& lhs=e.comp_trans_variant->dest;
  expr& rhs=e.comp_trans_variant->arg;
  assert(lhs.kind == subscription);
  sub s = lhs.subscription_variant;
  assert(s->array.kind==applied_identifier); // grammar ensures this
  id_type aggr=s->array.identifier_variant;
  expr& index=s->index;
  bool reversed=s->reversed;
  id_type op = e.comp_trans_variant->op;
@/const type* aggr_tp; size_t d,o; bool is_const;
  bool is_local = (aggr_tp=layer::lookup(aggr,d,o,is_const))!=nullptr;
  if (not is_local and (aggr_tp=global_id_table->type_of(aggr,is_const))==nullptr)
    report_undefined(aggr,e,"component transform");
@.Undefined identifier@>
  if (is_const)
    report_constant_modified(aggr,e,"component transform");
@.Name is constant @>
  assert(not aggr_tp->is_polymorphic());
  // polymorphic variables are made constant
@)
  static expr_ptr appl; expr_ptr saved_appl;
  @< Set |appl| to the application of |op| to |lhs| and |rhs| while saving the
     previous value in |saved_appl| @>
  expression_ptr ind;
  type ind_tp = type::bottom(fc);
  type_expr comp_te;
  subscr_base::sub_type kind;
  @< Convert |index| to |ind|... @>
  expression_ptr call = convert_expr_strongly(*appl,fc,comp_te);

  expression_ptr result;
  @< If the conditions for an optimised in-place component transformation... @>
  @< Restore initial state of |*e.comp_trans_variant| and of |appl| @>
  return result;
}

@ All variables in the title were declared before, and the values set here may
be used in subsequent modules.

@< Convert |index| to |ind|, set |ind_t| to its type, and |comp_t| to the
   type of subscription of |*aggr_tp| by it; throw an exception if the |kind| of
   subscription does not allow assignment @>=
{
  ind = convert_expr(index,ind_tp);
@/kind=subscr_base::index_kind(aggr_tp->expanded(),ind_tp.bake(),comp_te);
  if (not subscr_base::assignable(kind))
  { std::ostringstream o;
    o << "Cannot assign to component of value of type " << *aggr_tp @|
      << " selected by index of type " << ind_tp
      << " in transforming assignment";
    throw expr_error(e,o.str());
  }
}

@ As said above, we can only construct a |component_transform| under certain
conditions. To sum up, we only consider row-of aggregates, we must find a
|builtin_call| after overload resolution, and we refuse implicit conversions.
The absence of the latter is tested by seeing if the converted expression has
the precise structure of the |expr| from which it was converted, as witnessed by
succeeding |dynamic_cast| invocations. It is because these conditions can be
determined only after type checking that we had to convert |appl|, even though
this |call| will not be used when producing a |component_transform|. In that
case we shall need |aggr|, |ind|, the resolved |c->f| and |c->name|, and the
second argument |tup->component[1]| of the call.

In the case that we must use a |component_assignment|, there is an additional
complication that we must avoid to use the index expression being evaluated
twice. That will be handled by evaluating the index in a |let| expression
generated on the spot, but since that requires some effort and has a slight run
time penalty, we avoid this complication for simple enough index expressions:
integer denotations and applied identifiers.

@< If the conditions for an optimised in-place component transformation are met,
   construct a structure derived from |component_transform| from pieces
   of |*call| and set |result| by passing it through |conform_types|,
   otherwise build a |comp_assignment| @>=
{ const builtin_call* c;
  const tuple_expression* tup = nullptr;
  const expression_base* arg = nullptr;
  if (kind==subscr_base::row_entry)
  { c = dynamic_cast<const builtin_call*>(call.get());
    if (c!=nullptr)
    {
      tup = dynamic_cast<const tuple_expression *>(arg=c->argument.get());
      if (tup!=nullptr)
        arg=tup->component[0].get();
    }
  }
  if (dynamic_cast<const subscr_base*>(arg)!=nullptr)
  @< Set |result| to a |component_transform| assembled from |aggr|, |ind|,
     |c->f|, |c->name|, and maybe |tup->component[1]|, passed through
   |conform_types| from |comp_t| to |tp| @>
  else
    if (@< |index| is a constant expression @>@;@;)
@/@< Set |result| to a |component_assignment| assembled from |aggr|, |ind|,
     |call|, and |kind|,
     passed through |conform_types| from |comp_t| to |tp| @>
   else
  @< Set |result| to a \&{let} expression that binds the index expression
     to a hidden identifier, within which |convert_expr| is applied (again) to
     a |component_assignment| containing a modified version of |appl| @>

}

@ This is straightforward; the $4$ classes derived from
|component_transform| all need their own line.

@< Set |result| to a |component_transform| assembled... @>=
{ expression_ptr nil(nullptr);
  auto& rhe =
    tup==nullptr ? nil : const_cast<expression_ptr&>(tup->component[1]);
  if (is_local)
  { if (reversed)
    @/ result.reset
      (new local_component_transform<true>@|
        (aggr,std::move(ind),d,o,std::move(rhe),c->f,c->name,e.loc));
    else
    @/ result.reset
      (new local_component_transform<false>@|
        (aggr,std::move(ind),d,o,std::move(rhe),c->f,c->name,e.loc));
  }
  else
  { if (reversed)
    @/ result.reset
      (new global_component_transform<true>@|
        (aggr,std::move(ind),std::move(rhe),c->f,c->name,e.loc));
    else
    @/ result.reset
      (new global_component_transform<false>@|
        (aggr,std::move(ind),std::move(rhe),c->f,c->name,e.loc));
  }
  result = conform_types(comp_te,tp,std::move(result),e);
}

@ When we cannot build a |component_transform|, we build a
|component_assignment|, translating $v[i]\mathrel\star:=E$ as
``$v[i]:=v[i]\star{E}$''. In this case we do use the full converted |call|, as
well as |kind|, which is not restricted (apart from being assignable as was
already tested) here.

@< Set |result| to a |component_assignment|... @>=
{ if (is_local)
  { if (reversed)
    @/ result.reset
      (new local_component_assignment<true>@|
        (aggr,std::move(ind),d,o,std::move(call),kind));
    else
    @/ result.reset
      (new local_component_assignment<false>@|
        (aggr,std::move(ind),d,o,std::move(call),kind));
  }
  else
  { if (reversed)
    @/ result.reset
      (new global_component_assignment<true>@|
        (aggr,std::move(ind),std::move(call),kind));
    else
    @/ result.reset
      (new global_component_assignment<false>@|
        (aggr,std::move(ind),std::move(call),kind));
  }
  result = conform_types(comp_te,tp,std::move(result),e);
}

@ The most frequent case of an expression type that we shall recognise as
guaranteed without side effects is that of |applied_identifier| expressions, but
we shall also recognise |integer_denotation| expressions (for various kinds of
subscription) and, only for matrix subscriptions, $2$-tuples of two identifier
or denotation expressions. The last part a bit cumbersome to test, and inside an
expression we cannot introduce variables to help us, but since we already
successfully type-checked the subscription expression, we are at least sure that
any tuple expression used as index is a $2$-tuple, so we don't test for that.

It is somewhat questionable whether using a |let| expression to bind a constant
index pair for a matrix subscription, and using it twice, would actually have
been less efficient than evaluating the pair of constant indices twice.
However, the |let| expression would force the actual formation of a $2$-tuple to
be bound to the hidden variable, and then use |push_expanded| to get the
components back on the stack (where they were before the $2$-tuple was formed),
while the direct approach never forms a $2$-tuple; therefore our guess would be
that the |let| solution is indeed less efficient, so that trying to avoid it (as
is done here) is justified.

@< |index| is a constant expression @>=
index.kind==applied_identifier or
index.kind==integer_denotation or @|
(index.kind==tuple_display and @|
 (index.sublist->contents.kind==applied_identifier or
  index.sublist->contents.kind==integer_denotation
 )
 and @|
 (index.sublist->next->contents.kind==applied_identifier or
  index.sublist->next->contents.kind==integer_denotation
 )
)

@ In this final case $v[I]\mathrel\star:= E$ is translated into the equivalent
of ``\&{let}~$\$=I$~\&{in}~$v[\$]:=v[\$]\star E$'', where $\$$ is a local
variable that can't conflict with $v$ or any names used in the expression~$E$.
(By the way, the original implementation of operation-assign to an aggregate
component was to always do this rewriting; it was achieved during syntax tree
construction in \.{parsetree.w} with relative ease.) Since here the lexical
level of $E$ is deeper here than in the earlier conversion, we cannot extract
and use a part of the converted |call| as we did above, and rather convert the
expression again in a modified setting.

Here is the plan: get the |id_type| value for the hidden identifier $\$$,
construct an applied identifier expression for this identifier, swap it out with
the index expression $I$ that is held, at location |index|, as subexpression of
the |call| that was set to represent $v[I]\mathrel\star{E}$ (and which already
has its converted form in |ind|), then build a |layer| as when
processing \&{let}~$\$=I$, then convert the modified |appl| in this new context
and as a part of a new |let_expression| using |ind| that wraps everything up,
and finally return the |let_expression|.

Since we now pass our locally built |comp_assignment_node| through another call
of |convert_expr|, we again need a |static| variable to hold it for correct
error reporting, and the old contents needs to be saved in a local variable as
was the case for |appl|. We must also make sure that the |index| expression,
which does not need a new conversion but must be swapped out for a variable, is
put back in place so that our expressions |e| will be intact on successful
completion,

@< Set |result| to a \&{let} expression... @>=
{
  id_type hidden = lookup_identifier("$");@q$@>
  id_pat dollar(hidden);
  expr temp(hidden,index.loc,expr::identifier_tag()); // applied $\$$
  index.swap(temp); // modify |call| to use $\$$
  static expr_ptr ca; expr_ptr saved_ca = std::move(ca);
  ca.reset(new expr(new comp_assignment_node @|
     {aggr
     ,expr(hidden,index.loc,expr::identifier_tag())
     ,std::move(*appl)
     ,reversed},e.loc));
  layer let_layer(1);
  thread_bindings(dollar,ind_tp.bake_off(),fc,let_layer,true);
  result.reset(new let_expression @|
    (dollar,std::move(ind),convert_expr(*ca,tp)));
  *appl = std::move(ca->comp_assign_variant->rhs);
    // restore |appl| for further restoration
  index.swap(temp); // restore our original expression for outer error reporting
  ca = std::move(saved_ca);
    // restore state of static variable now that no error was thrown
}

@* Some special wrapper functions.
%
In this chapter we define some wrapper functions that are not accessed through
the overload table; they must be directly visible to the type-checking code
that inserts them, which is why they are defined as local functions to the
current \.{axis.w} module.

@< Static variable definitions that refer to local functions @>=
static shared_builtin sizeof_row_builtin =
    std::make_shared<const builtin>(sizeof_wrapper,"#@@[T]",0);
static shared_builtin sizeof_vector_builtin =
    std::make_shared<const builtin>
      (sizeof_vector_wrapper,"#@@vec",0);
static shared_builtin sizeof_ratvec_builtin =
    std::make_shared<const builtin>
      (sizeof_ratvec_wrapper,"#@@ratvec",0);
static shared_builtin sizeof_string_builtin =
    std::make_shared<const builtin>
       (sizeof_string_wrapper,"#@@string",0);
static shared_builtin matrix_columns_builtin =
    std::make_shared<const builtin>
     (matrix_ncols_wrapper,"#@@mat",0);
static shared_builtin sizeof_parampol_builtin =
    std::make_shared<const builtin>
      (virtual_module_size_wrapper, "#@@ParamPol",0);
static shared_variadic_builtin print_builtin =
  std::make_shared<const builtin_value<true> >(print_wrapper,"print@@T",0);
static shared_variadic_builtin to_string_builtin =
  std::make_shared<const builtin_value<true> >
  (to_string_wrapper,"to_string@@T",0);
static shared_variadic_builtin prints_builtin =
  std::make_shared<const builtin_value<true> >(prints_wrapper,"prints@@T",0);
static shared_variadic_builtin error_builtin =
  std::make_shared<const builtin_value<true> >(error_wrapper,"error@@T",0);
static shared_builtin prefix_elt_builtin =
  std::make_shared<const builtin>
    (prefix_element_wrapper,"#@@(T,[T])",2);
static shared_builtin suffix_elt_builtin =
  std::make_shared<const builtin>
    (suffix_element_wrapper,"#@@([T],T)",1);
static shared_builtin join_rows_builtin =
  std::make_shared<const builtin>
    (join_rows_wrapper,"##@@([T],[T])",0);
static shared_builtin join_rows_row_builtin =
  std::make_shared<const builtin>
    (join_rows_row_wrapper,"##@@([[T]])",0);
static shared_builtin boolean_negate_builtin =
  std::make_shared<const builtin>(bool_not_wrapper,"not@@bool",0);

@ Finally we define the Boolean negation wrapper function.
@< Local function definitions @>=
void bool_not_wrapper(eval_level l)
{ bool b=get<bool_value>()->val;
  if (l!=eval_level::no_value)
    push_value(whether(not b));
}
@)


@* Index.

% Local IspellDict: british
