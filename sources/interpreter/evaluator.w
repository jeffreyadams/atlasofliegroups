% Copyright (C) 2006 Marc van Leeuwen
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

@* Outline.
%
This file describes the central part of the interpreter for the (new) command
language of the Atlas of Lie Groups and Representation software. This part is
concerned with the analysis and execution of expressions that have already
been processed by the parser. While for a long time during development this
part consisted of a single compilation unit, the increasing size of that unit
has led us to break it up into several parts that can be compiled
individually. The initial compilation unit describes the fundamental type
declarations used throughout the interpreter, as well as some functions that
operate on types, for instance to see if one con be converted into another.
This compilation unit is called \.{types}, which refers both to user types as
to meta-types of the interpreter used to represent user types and other
fundamental notions.

For the global organisation of these compilation unit we need many similar
modules. In order to distinguish the units, we have systematically used module
names starting with ``Initial'' for the outer level modules of the \.{types}
compilation unit; by choosing otherwise similar names it is easy to move code
between the compilation units if necessary, by changing a module name in the
defining section for the code.

@( types.h @>=

#ifndef TYPES_H
#define TYPES_H

@< Includes needed in \.{types.h} @>@;
namespace atlas { namespace interpreter {
@< Initial type definitions @>@;
@< Initial declarations of global variables @>@;
@< Initial declarations of exported functions @>@;
@< Initial template and inline function definitions @>@;
}@; }@;
#endif

@ The second compilation unit \.{evaluator} implements type analysis and
evaluation of expressions, and the current source file takes its name from
this compilation unit. We may decide to relegate the implementation of some
basic user types like vectors and matrices, which does not interact much with
the general evaluation procedure, to yet another compilation unit, but this
separation has not yet been made. As for the textual organisation of this
source file, we shall try to group contributions according to the compilation
units, from most fundamental to most applied; this will facilitate breaking up
the source file itself if its size should get unwieldy, and of more immediate
importance, will cause the compilable files of the earlier units to remain
unchanged (in particular in its source file line references) and therefore
without need of recompilation when only later units are being edited.

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

@ Each compilation unit follows a similar global pattern. While some modules
are absent in each case, the order in the implementation files is local type
definitions, global variables, local variables, local functions, global
functions (the last point being is the main goal of the implementation unit).

@( types.cpp @>=

#include "types.h"
#include <cstdlib>

namespace atlas { namespace interpreter {
namespace {
@< Initial local type definitions @>@;
}@; // |namespace|
@< Initial global variable definitions @>@;
namespace {
@< Initial local function definitions @>@;
}@; // |namespace|
@< Initial function definitions @>@;
}@; }@;

@ For the \.{evaluator} compilation unit we omit the ``Initial'' designation
in module names.

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


@* Types used throughout the evaluator.
%
The parser produces a parse tree, in the form of a value of type |expr|
defined in the unit \.{parsetree}. Our task is to take such an expression,
analyse it and then evaluate it to obtain a value. At various points errors
can occur, which we shall have to handle gracefully: during the analysis
static ``type'' errors may prevent us from undertaking any meaningful action,
and in absence of such errors there still can be dynamic ``runtime'' errors.

We first define some data types used by the evaluator. First of all we shall
consider ``type'' values, represented by the structure |type_expr|, in
terms of which the static analysis is performed. Then we shall define the
structure of values for user data, the kind manipulated by the evaluator at
runtime; these are of classes derived from~|value_base|. Then we shall define
values that represent the user expression after type analysis, and which serve
to control the runtime actions; these are of classes derived from
|expression_base|. Originally the |expr| value itself (possibly slightly
modified during type analysis) was used for this purpose, but it turns out to
be useful to rebuild the expression tree; as a side benefit we can liberate
ourselves from the constraint on data types built by the parser, that they can
be understood by a \Cee-program. Finally we need some types for errors that we
might have to throw; these are the classes |program_error|, and |type_error|
which is derived from it.

Most of these types are recursive in a manner more complicated that simple
linked lists or binary trees: at each level they represent various alternative
possibilities, each of which might recursively refer to the original type one
or more times. Such recursion is represented by linked structures that cannot
easily be hidden in a (template) class with a simple interface, such as
happens in the contained classes of the Standard Template Library. For this
reason the code using these types is to some extent exposed to the linked
nature of the data, and some concern for proper memory management will be
necessary. We shall however provide functions to facilitate safe and
consistent handling of storage of these data, while avoiding excessive
duplication during manipulation. Altogether this part of the program will be
of a quite different nature than that of the main mathematical library.


@< Includes needed in \.{types.h} @>=
#include <memory> // for |std::auto_ptr|

@*1 Types to represent types.
%
We make a fundamental choice to check types before attempting to execute any
expression entered by the use; thus type errors can be signalled earlier and
expressed more understandably than if would wait until an accident is about to
happen at runtime. This means that types are to be represented independently
of runtime values. Of course the interpreter executes the same code regardless
of the types of the values it manipulates; those values are then accessed via
generic pointers. The runtime values do also contain an indication of their
actual type, which will allow us to do some dynamic testing to trap possible
errors in the logic of our interpreter, but the user should never notice this.

Types are represented by small tree structures, with nodes of types
|type_expr| that will be detailed later. These expose the structure of the
represented type rather than attempting to hide it behind a type identifier
(as \Cpp~classes do). Note however that classes defined in the Atlas software
library itself will be presented to the user as primitive types, so we do
provide abstraction there. A simple but important example is a function type
that specifies the types of the individual function arguments and returned
values.

Different types will not share any sub-trees among each other, which simplifies
memory management. This choice means some recursive copying of tree structures
is sometimes required, but we avoid doing so more than absolutely necessary.

Usually types are built from the bottom up (from leaves to the root), although
during type checking the reverse also occurs. In bottom-up handling of types,
they should not reside in local |type_expr| variables, since these would
require (recursive) copying in order to become a descendant type of another
type. Instead they will be referred to by local pointer values, and in order
to ensure exception-safety, these should be auto-pointers; for this reason we
define |type_ptr| to be an auto-pointer type. For the links between a node and
its descendant types however we use ordinary pointers of type~|type_p|. We
also define a |type_list| type that will be used inside tuple types, and a
corresponding auto-pointer type |type_list_ptr|.

@< Initial type definitions @>=
struct type_expr;
typedef type_expr* type_p;
typedef std::auto_ptr<type_expr> type_ptr;
typedef struct type_node* type_list;
typedef std::auto_ptr<type_node> type_list_ptr;

@ Since types and type lists own their trees, their copy constructors must
make a deep copy. Most invocations of the copy constructor for types, from
outside the implementation of the |type_expr| class itself, will take
place by an explicit call of a function |copy| for improved visibility. This
function converts a reference to a |type_expr| into an auto-pointer to a
copied instance of that type expression.

@< Initial declarations of exported functions @>=
type_ptr copy(const type_expr& t);

@~The function |copy| simply creates an auto-pointer from the call of |new|
for |type_expr(t)|, where the latter invokes the copy constructor of
|type_expr|. The latter, has a recursive definition (given in section
@#type expression copy@>) that performs a deep copy, but the details don't
concern us here.

@< Initial function definitions @>=
type_ptr copy(const type_expr& t) @+
{@; return type_ptr(new type_expr(t)); }

@*2 Type lists.
%
The auxiliary type |type_list| gives a preview of matters related to
ownership, that will also apply, in a more complex setting, to
|type_expr|. It is a simply linked list, whose nodes contain a
|type_expr| as a sub-object. We might instead have used a pointer there,
but that would be less compact and require more memory (de)allocations. On the
other hand it forced us to introduce the |type_expr::set_from| method
(see below) to make a shallow copy, avoiding the deep copy that would result
from initialising |t(*head)|. Note how the pointer on the new node is passed
as an auto-pointer that is immediately released when inserted into the node;
this is the unique way to handle this that assures there is exactly one owner
at all times (if an ordinary pointer had been passed, the pointed-to structure
might fail to be deleted if the allocation of this node should throw an
exception).

Because a |type_node| contains a
|type_expr| we must make sure its definition given here \emph{follows}
that of |type_expr|, at least in the tangled code.

@< Definition of |struct type_node@;| @>=
struct type_node
{ type_expr t; @+ type_list next;
@)
  type_node(type_ptr head, type_list_ptr tail) : t(),next(tail.release())
  @+{@; t.set_from(head); }
  type_node(const type_node& n);
  ~type_node();
};

@ Copy-constructing a type list is easy: the recursive copy of the tail list
is copy-constructed into the fresh node. This is safe because the node
contains only one pointer: if |new| succeeds, the construction is immediately
terminated, without any occasion for an exception causing loss of the pointer
before our destructor is in place. If there had been two pointers in the node,
then some protection against a memory leak would be needed (since |new| might
throw); another advantage of not using a pointer for the |t| field.

@< Initial function definitions @>=
type_node::type_node(const type_node& n)
 : t(n.t),next(n.next==NULL ? NULL : new type_node(*n.next))
@+{}

@ The destructor for |type_node| will possibly be called implicitly by the
destruction of a |type_list_ptr|. The destructor for variants containing a
|type_list| should call |delete| for that pointer, which will also call the
destructor below, which will recursively clean up the whole list.

@< Initial function def... @>=
type_node::~type_node()
@+{@; delete next; }

@ Type lists are usually built by calling |make_type_singleton| (to get
started) or |make_type_list| (to prefix a type to an existing list), rather
than using a constructor. These functions are efficient, due to the use of
auto-pointers. A caller holding auto-pointers to the descendants transfers
ownership to the new node at the call (the copied auto-pointers are released
at precisely the right moment), and regains ownership of the new structure
upon return. Using ordinary pointers one would run into ownership conflicts
that can only be solved by calling the (recursive) copy constructor, which
would also be necessary when using local |type_node| variables.

@< Initial declarations of exported functions @>=
type_list_ptr make_type_singleton(type_ptr t);
type_list_ptr make_type_list(type_ptr t,type_list_ptr l);
size_t length(type_list l);

@ We pass the auto-pointers by value, rather than by (necessarily non-constant)
reference. This implies a small overhead due to multiple ownership transfer,
but allows (due to a clever definition of auto-pointers) passing anonymous
values obtained from a constructor or construction function as argument, which
would be forbidden if using non-constant reference parameters.

@< Initial function def... @>=

type_list_ptr make_type_singleton(type_ptr t)
{@; return type_list_ptr(new type_node(t,type_list_ptr(NULL))); }
@)
type_list_ptr make_type_list(type_ptr t,type_list_ptr l)
{@; return type_list_ptr(new type_node(t,l)); }

@ Here is one more function that is convenient to have around.

@< Initial function definitions @>=
size_t length(type_list l)
@+{@; size_t len=0;
  for (; l!=NULL; l=l->next) ++len;
  return len;
}

@*2 Primitive types.
%
Types can be formed in a variety of ways (primitive types, tuple types,
function types, etc.), and correspondingly will involve a |union| of various
alternatives. So before we can define |type_expr|, we need to enumerate
its variants and the possibilities for primitive types. Some of them will be
introduced later.

@< Initial type definitions @>=
enum type_tag @+
{ undetermined_type, primitive_type, row_type, tuple_type, function_type };

enum primitive_tag
{ integral_type, rational_type, string_type, boolean_type
  , @< Other primitive type tags @>@;@; @|
  @+nr_of_primitive_types };

@ For printing types (and later for parsing them) we shall need names for the
primitive ones.

@< Initial declarations of global variables @>=

extern const char* prim_names[];

@~Here is the list of names of the primitive types (some are given later),
terminated by a null pointer. The last name |"void"| is not a primitive name,
corresponds to |nr_of_primitive_types|, and will be treated exceptionally in
|make_prim_type| to make an empty tuple type instead.

@< Initial global variable definitions @>=
const char* prim_names[]=@/
{"int","rat","string","bool",@< Other primitive type names@>@;@;@+
 "void", NULL };

@*2 Type expressions.
%
We have a simple but flexible type model. There is a finite number of
``primitive'' types, many of which are abstractions of complicated classes
defined in the Atlas library, such as root data or reductive groups.
Furthermore one has types that are ``row of'' some other type, tuples
(Cartesian products) of some given sequence of types, and function types. We
also allow for an undetermined type, which can serve as a wild-card to specify
type patterns. Here are enumerations of tags for the basic kinds of types, and
for the sub-cases of |primitive_type|.

Type expressions are defined by a tagged union. Code accessing the union will
in general test the tag and access the corresponding variant directly, so we
give public access to the data fields (it is not obvious how coherent variant
selection could be enforced by private data and public methods without greatly
complicating usage; the \Cpp~model of abstraction does not seem particularly
suited for tagged unions). By using an anonymous union, the field selectors
like |prim| of the variants in the union can be used directly on the level of
the |type_expr|, thus avoiding an additional level of selection to
access them.

There is one restriction on types that is not visible in the definition below,
namely that the list of types referred to by the |tuple| field cannot have
length~$1$ (but length~$0$ is allowed). This is because anything that would
suggest a $1$-tuple (for instance a parenthesised expression) is identified
with its unique component.

We forbid the assignment operator, since the default definition would not do
the right thing, the proper definition would be somewhat complicated, and
without much use since types do not change dynamically. Two special cases are
provided for however, for replacing an undefined (or partially undefined) type
by a more specific type (note that in these cases nothing disappears, so no
clean-up is necessary). The first case, the |set_from| method, makes a shallow
copy of its argument into the object for which is was called, which is
required to be undetermined initially; it is intended for use in constructors
to avoid the deep copy of the copy constructor. Since the shallow copy
incorporates any descendants without copying, the caller of |set_from| should
have and relinquish ownership of the argument type passed; therefore the
argument is passed as an auto-pointer by value. The second case, the
|specialise| method, is used during type analysis, to see if our type matches
a given pattern, or in case it was (partially) undefined whether it can be
made to match the pattern by if necessary replacing some undetermined
descendants by more specific ones. The call returns a value indicating whether
this was possible, and if so makes the necessary specialisations to our type.
That is done by copying, so the caller does not require or lose ownership of
the pattern for this method.

@< Initial type definitions @>=

struct func_type; // must be predeclared for |type_expr|
@)
struct type_expr
{ type_tag kind;
  union
  { primitive_tag prim; // when |kind==primitive|
    type_p component_type; // when |kind==row_type|
    type_list tuple; // when |kind==tuple_type|
    func_type* func; // when |kind==function_type|
  };
@)
  type_expr() : kind(undetermined_type) @+{}
  explicit type_expr(primitive_tag p)
    : kind(primitive_type) @+{@; prim=p; }
  explicit type_expr(type_ptr c)
    : kind(row_type) @+{@; component_type=c.release(); }
  explicit type_expr(type_list_ptr l) : kind(tuple_type)
  @+{@; tuple=l.release(); }
  inline type_expr(type_ptr arg, type_ptr result); // for function types
@)
  type_expr(const type_expr& t); // copy constructor
  ~type_expr();
@)
  void set_from(type_ptr p);
  bool specialise(const type_expr& pattern);
private:
  type_expr& operator=(const type_expr& t);
   // assignment is forbidden
};

@< Definition of |struct type_node@;| @>@; // this must \emph{follow}

@ The constructor for function types cannot be defined inside the structure
definition, since |func_type| is not (and cannot be) a complete type there.
The constructor for |func_type| used will be defined below.

@< Initial template and inline function definitions @>=
type_expr::type_expr(type_ptr arg, type_ptr result)
   : kind(function_type) {@; func=new func_type(arg,result); }

@ The variant for function types needs both an argument type and a result
type. We cannot define the appropriate structure directly within the |union|
of |type_expr| where it is needed, since the scope of that definition
would then be too limited to perform the appropriate |new| in the constructor.
Therefore we must declare and name the structure type before using it in
|type_expr|. Afterwards we define the structure, and while we are doing
that, we also define a constructor to fill the structure (it was used above)
and a copy constructor.

The |func_type| structure has two sub-objects of type |type_expr|,
rather than pointers to them. This is more compact and practical for most
purposes, but we should take care that the basic constructor only makes a
shallow copy of the topmost nodes pointed to. This is achieved by using the
|set_from| method to fill the initially undetermined slots with a shallow copy
of the types, in the same way as was done in the constructor for~|type_node|.
The copy constructor needs no special care, as it makes a deep copy.

@s result_type normal

@< Initial type definitions @>=
struct func_type
{ type_expr arg_type, result_type;
@)
  func_type(type_ptr a, type_ptr r)
   : arg_type(), result_type()
  {@; arg_type.set_from(a); result_type.set_from(r); }
  func_type(const func_type& f)
   : arg_type(f.arg_type),result_type(f.result_type) @+{}
};


@ As we remarked above, the copy constructor for a |type_expr|
recursively copies the descendant types; this is necessary because these
descendants are owned by the object. As usual the recursion is implicit in the
use of a copy constructor within |new|. In all cases there is only one object
directly pointed to, so there is no need for intermediate auto-pointers: no
exception can be thrown between the moment that the call to |new| returns (if
it does), and the completion of our constructor, which incorporates the
returned pointer.
@:type expression copy@>

@< Initial function definitions @>=
type_expr::type_expr(const type_expr& t) : kind(t.kind)
{ switch (kind)
  { case undetermined_type: break;
    case primitive_type: prim=t.prim; break;
    case row_type:
      component_type=new type_expr(*t.component_type);
    break;
    case tuple_type:
      tuple= t.tuple==NULL ? NULL : new type_node(*t.tuple);
    break;
    case function_type: func=new func_type(*t.func); break;
  }
}

@ The destructor must similarly clean up afterwards, with the recursion again
being implicit.

@< Initial function definitions @>=
type_expr::~type_expr()
{ switch (kind)
  { case undetermined_type: case primitive_type: break;
    case row_type: delete component_type; break;
    case tuple_type: delete tuple; break;
    case function_type: delete func; break;
  }
}

@ The method |set_from| is like an assignment operator, but it avoids making a
deep copy. Its argument is an auto-pointer passed by value, to remind users
that this argument must be |new|-allocated and that they are giving up
ownership of this argument. The contents of its top-level structure will be
copied to the current |type_expr|, but without invoking the copy
constructor; afterwards the original node will be destroyed after detaching
any possible descendants from it. This operation is only safe if the
|type_expr| previously had no descendants, and in fact we insist that it
had |kind==undetermined_type|; if this condition fails we signal a
|std::logic_error|. In a sense this is like a |swap| method, but only defined
if the type was undetermined to begin with; in modern parlance, it implements
move semantics.

@h <stdexcept>
@< Initial function definitions @>=
void type_expr::set_from(type_ptr p)
{ if (kind!=undetermined_type) throw std::logic_error("Illegal set_from");
  kind=p->kind;
  switch(kind) // copy top node
  { case undetermined_type: break;
    case primitive_type: prim=p->prim; break;
    case row_type: component_type=p->component_type; break;
    case function_type: func=p->func; break;
    case tuple_type: tuple=p->tuple; break;
  }
  p->kind=undetermined_type;
  // detach descendants, so auto-pointer destroys top-level only
}


@ The |specialise| method is mostly used to either set a completely
undetermined type to a given pattern or to test if it already matches it;
however, we do not exclude the possibility that a partly determined type is
modified by specialisation of one of its descendants to match the given
pattern. Matching a pattern means being at least as specific, and specialising
means replacing by something more specific, so we are (upon success) replacing
the type for which the method is called by the most general unifier (i.e.,
common specialisation) of its previous value and |pattern|. Therefore the name
of this method is somewhat misleading, in that the specialisation is not
necessarily to |pattern|, but to something matching |pattern|.

In the case of an |undetermined_type|, |specialise| calls the copy constructor
of the specified type via a placement-|new| into the undetermined type
indication. This overwrites the old value, which is harmless since no
destruction is necessary for the |undetermined_type| case. In the other cases
we only continue if the top levels of both type declarers match, in which case
we try to recursively specialise all descendants. We do not guarantee
commit-or-roll-back, in other words, when the specialisation fails, some
modifications to our type may still have been made. This is not very serious,
since types are not shared and the initial call of specialisation always
starts with a fresh type, so at worst it might lead to a slightly misleading
type appearing in an error message (but even this is not easy to produce).

@< Initial function definitions @>=
bool type_expr::specialise(const type_expr& pattern)
{ if (pattern.kind==undetermined_type)
    return true; // specialisation to \.* trivially succeeds.
  if (kind==undetermined_type)
    {@;new(this) type_expr(pattern); return true; }
  if (pattern.kind!=kind) return false; // impossible to refine
  switch(kind)
  { case primitive_type: return prim==pattern.prim;
    case row_type: return component_type->specialise(*pattern.component_type);
    case function_type:
      return func->arg_type.specialise(pattern.func->arg_type) @|
         and func->result_type.specialise(pattern.func->result_type);
    case tuple_type:
     @< Try to specialise types in |tuple| to those in |pattern.tuple|,
        and |return| whether this succeeded @>
    default: return true; // to keep the compiler happy, cannot be reached
  }
}

@ For tuples, specialisation is done component by component. We do not check
beforehand that the lengths of the lists match, so we must be prepared for
either one of the lists running out before the other does.

@< Try to specialise types in |tuple| to those in |pattern.tuple|... @>=
{ type_list l0=tuple, l1=pattern.tuple;
  while (l0!=NULL and l1!=NULL and l0->t.specialise(l1->t))
  {@; l0=l0->next; l1=l1->next; }
  return l0==NULL and l1==NULL;
  // we succeeded only if both lists terminate simultaneously
}

@ For printing types, we shall pass |type_expr| values to the
operator~`|<<|' by constant reference, which seems more decent than doing
so by pointer (which would override the definition that simply prints the
hexadecimal address of a pointer); we shall not define instances of~`|<<|' for
other pointer types either. Since we often hold types in |type_ptr| values,
this does mean the we must dereference explicitly in printing.

@< Initial declarations of exported functions @>=
std::ostream& operator<<(std::ostream& out, const type_expr& t);

@~The cases for printing the types are fairly straightforward. Only
function types are somewhat more involved, since we  want to suppress
additional parentheses around argument and result types in case these are
tuple types; defining a separate operator for a |type_list| facilitates our
task a bit.

@< Initial function definitions @>=

std::ostream& operator<<(std::ostream& out, type_list l)
{ for (; l!=NULL; l=l->next) out << l->t << ( l->next!=NULL ? "," : "" );
  return out;
}
@)
std::ostream& operator<<(std::ostream& out, const type_expr& t)
{ switch(t.kind)
  { case undetermined_type: out << '*'; break;
    case primitive_type: out << prim_names[t.prim]; break;
    case row_type: out << '[' << *t.component_type << ']'; break;
    case tuple_type:
      out << '(' << t.tuple << ')' ;
    break;
    case function_type:
      out << '(';
      if (t.func->arg_type.kind==tuple_type)
         out << t.func->arg_type.tuple; // naked tuple
      else out << t.func->arg_type; // other component type
      out << "->";
      if (t.func->result_type.kind==tuple_type)
         out << t.func->result_type.tuple; // naked tuple
      else out << t.func->result_type; // other component type
      out << ')'; break;
  }
  return out;
}

@ Finally we need a comparison for structural equality of type
declarators.

@< Initial declarations of exported functions @>=
bool operator== (const type_expr& x,const type_expr& y);
inline bool operator!= (const type_expr& x,const type_expr& y)
{@; return !(x==y); }

@~This code is quite similar to the |specialise| method; in fact one could
often use that method instead of the equality operator, but here we want both
operands to be |const|.

@< Initial function definitions @>=
bool operator== (const type_expr& x,const type_expr& y)
{ if (x.kind!=y.kind) return false;
  switch (x.kind)
  { @+ default:
// all cases are listed below, but somehow this keeps \.{g++} from complaining
  @\case undetermined_type: return true;
    case primitive_type: return x.prim==y.prim;
    case row_type: return *x.component_type==*y.component_type;
    case tuple_type:
    { type_list l0=x.tuple, l1=y.tuple;
      while (l0!=NULL and l1!=NULL and l0->t==l1->t)
	{@; l0=l0->next; l1=l1->next; }
      return l0==NULL and l1==NULL;
      // lists must end simultaneously for success
    }
    case function_type:
      return  x.func->arg_type==y.func->arg_type
	 and   x.func->result_type==y.func->result_type;
  }
}

@ Instead of using the constructors directly, we usually use the constructing
functions below, which take and return |type_ptr| or |type_list_ptr| values,
which are auto-pointers.
@< Initial declarations of exported functions @>=
type_ptr make_undetermined_type();
type_ptr make_prim_type(primitive_tag p);
type_ptr make_row_type(type_ptr c);
type_ptr make_tuple_type (type_list_ptr l);
type_ptr make_function_type (type_ptr a, type_ptr r);

@ The reasons for, and the way of, using them are the
same as explained above for |make_type_list|. Note that |make_prim_type|
contrary to what its name suggests returns and empty tuple type for the type
name |"void"|.

@< Initial function definitions @>=
type_ptr make_undetermined_type()
@+{@; return type_ptr(new type_expr); }
@)
type_ptr make_prim_type(primitive_tag p)
{ return p<nr_of_primitive_types ?
    type_ptr(new type_expr(p)) :
    type_ptr(new type_expr(type_list_ptr(NULL)));
}
@)
type_ptr make_row_type(type_ptr c)
{@; return type_ptr(new type_expr(c)); }
@)
type_ptr make_tuple_type (type_list_ptr l)
{@; return type_ptr(new type_expr(l)); }
@)
type_ptr make_function_type (type_ptr a, type_ptr r)
{@; return type_ptr(new type_expr(a,r));
}

@*2 A parser interface to constructing types.
%
In order for the parser, which is (currently) compiled as \Cee\ code, to be
able to call the type-constructing functions, we provide functions
with \Cee-linkage and types understandable from \Cee, which call the above
functions. Their prototypes, visible to the parser, were defined in the
module \.{parsetree.w}.

@< Includes needed in \.{types.h} @>=
#include "parsetree.h"

@ All that happens here is a painful conversion from
void pointer (the type |ptr|) to auto-pointer (or from integer to enumeration
type), and a somewhat less painful conversion back. Nothing can throw during
these conversions, so passing bare pointers is exception-safe.

@< Initial function def... @>=
extern "C"
ptr mk_type_singleton(ptr t)
{@;
  return make_type_singleton(type_ptr(static_cast<type_p>(t)))
  .release();
}
extern "C"
ptr mk_type_list(ptr t,ptr l)
{ return make_type_list(type_ptr(static_cast<type_p>(t)),@|
                        type_list_ptr(static_cast<type_list>(l))).release(); }
extern "C"
ptr mk_prim_type(int p)
{@; return make_prim_type(static_cast<primitive_tag>(p)).release(); }
extern "C"
ptr mk_row_type(ptr c)
{@; return make_row_type(type_ptr(static_cast<type_p>(c)))
  .release(); }
extern "C"
ptr mk_tuple_type(ptr l)
{@; return make_tuple_type(type_list_ptr(static_cast<type_list>(l)))
  .release(); }
extern "C"
ptr mk_function_type(ptr a,ptr r)
{ return make_function_type(type_ptr(static_cast<type_p>(a)),@|
    type_ptr(static_cast<type_p>(r))).release(); }
@)
extern "C"
void destroy_type(ptr t)@+
{@; delete static_cast<type_p>(t); }
extern "C"
void destroy_type_list(ptr t)@+
{@; delete static_cast<type_list>(t); }

ptr first_type(ptr typel)@+
{@; return &static_cast<type_list>(typel)->t; }

std::ostream& print_type(std::ostream& out, ptr type)
{@; return out << *static_cast<type_p>(type); }

@*2 Specifying types by strings.
%

In practice we shall rarely call functions like |make_prim_type| and
|make_row_type| directly to make explicit types, since this is rather
laborious. Instead, such explicit types will be constructed by the function
|make_type| that parses a string, and correspondingly calls the appropriate
type constructing functions.

@< Initial declarations of exported functions @>=
type_ptr make_type(const char* s);

@ The task of converting a properly formatted string into a type is one of
parsing a simple kind of expressions. We are not going to write incorrect
strings (we hope) so we don't care if the error handling is crude. The
simplest way of parsing ``by hand'' is recursive descent, so that is what we
shall use. By passing a character pointer by reference, we allow the recursive
calls to advance the index within the string read.

@< Initial function definitions @>=
type_ptr scan_type(const char*& s);
type_list_ptr scan_type_list(const char*& s);
@)
type_ptr make_type(const char* s)
{ const char* orig=s;
  try {@; return scan_type(s); }   // provide an lvalue
  catch (std::logic_error e)
  { std::cerr << e.what() << "; original string: '" << orig @|
              << "' text remaining: '" << s << "'\n";
    throw@[@];
  // make the error hard to ignore; if thrown probably aborts the program
  }
}

type_ptr scan_type(const char*& s)
{ if (*s=='*') return ++s,make_undetermined_type();
  else if (*s=='[')
    @< Scan and |return| a row type, or |throw| a |logic_error| @>
  else if (*s=='(')
    @< Scan and |return| a tuple or function type,
       or |throw| a |logic_error| @>
  else @< Scan and |return| a primitive type, or |throw| a |logic_error| @>
}


@ Since we did not advance the pointer when testing for |'['|, we must start
with that.

@< Scan and |return| a row type, or |throw| a |logic_error| @>=
{ type_ptr c=scan_type(++s);
  if (*s++!=']') throw std::logic_error("Missing ']' in type");
    return make_row_type(c);
}

@ Now we do tuple and function types. Here again we start by advancing the
pointer. Otherwise the descent is still straightforward, thanks to
|scan_type_list|. After scanning a first list we decide whether this will be a
tuple type (if a right parenthesis follows) and record it in a boolean
variable~|is_tuple| to be able to share the code for constructing the tuple 
type.

The only complication is that single parenthesised types, and single argument
or return types should not be converted into tuple types with one component,
but just into the constituent type. This is done by assigning to the
|type_ptr| variables |a| and |r| either the extracted singleton type or the
tuple type constructed from the non-singleton list. Note that this is one of
the few places where we really use the ownership-tracking semantics of
auto-pointers, in the sense that their destruction behaviour at a certain
point is variable: the list pointed to by |l0| and |l1| will only be deleted
in the case of a singleton list, where the auto-pointer was not passed on in a
call to |make_tuple_type|.

@< Scan and |return| a tuple or function type, or |throw| a |logic_error| @>=
{ type_list_ptr l0=scan_type_list(++s), l1(NULL);
  bool is_tuple=*s==')';
  if (*s=='-' and *++s=='>') l1=scan_type_list(++s);
  if (*s++!=')') throw std::logic_error("Missing ')' in type");
  type_ptr a =
    l0.get()!=NULL and l0->next==NULL ? copy(l0->t) : make_tuple_type(l0);
  if (is_tuple) return a;
  type_ptr r =
    l1.get()!=NULL and l1->next==NULL ? copy(l1->t) : make_tuple_type(l1);
  return make_function_type(a,r);
}

@ A comma-separated list of types is handled by a straightforward recursion.
We must not forget that the list could be empty, which happens only of the
very first character we see is |')'| or |'-'|.

@< Initial function definitions @>=
@)
type_list_ptr scan_type_list(const char*& s)
{ if (*s==')' or *s=='-') return type_list_ptr(NULL);
  type_ptr head=scan_type(s);
  if (*s!=',') return make_type_singleton(head);
  return make_type_list(head,scan_type_list(++s));
}

@ For primitive types we use the same strings as for printing them. We test as
many characters as the type name has, and the fact that no alphanumeric
character follows, so that a longer type name will not match a prefix of it.

@h <cctype>
@< Scan and |return| a primitive type, or |throw| a |logic_error| @>=
{ for (size_t i=0; i<nr_of_primitive_types; ++i)
  { std::string name=prim_names[i];
    if (name.compare(0,name.length(),s,name.length())==0 and
        not isalpha(s[name.length()]))
    @/{@; s+=name.length();
      return make_prim_type(static_cast<primitive_tag>(i));
    }
  }
  throw std::logic_error("Type unrecognised");
}

@*2 Predefined type expressions.
%
We shall often need to refer to certain types for comparison or for providing
a required type context. Instead of generating them on the fly each time using
|make_type|, we define constant values that can be used everywhere. Three of
these types, for \.{int}, \.{bool} and \.{void}, need to be non-|const|, since
in the role of required type (as argument to |convert_expr| defined later)
they could potentially be specialised; if they were const, we would be obliged
to call |copy| for every such use. However since these particular types cannot
possibly be specialised, they are assured to remain constant even though they
are not declared as such.

@< Initial declarations of global variables @>=
extern const type_expr unknown_type; // \.{*}
extern type_expr void_type; // \.{()}
extern type_expr int_type; // \.{int}
extern type_expr bool_type; // \.{bool}
extern const type_expr row_of_type; // \.{[*]}
extern const type_expr gen_func_type; // \.{(*->*)}

@ The definition of the variables uses the constructors we have seen above,
rather than functions like |make_primitive_type| and |make_row_type|, so that
no dynamic allocation is required for the top level structure. For generic row
and function types we construct the |type_expr| from auto-pointers (of which
the constructor takes possession) pointing to other |type_expr|s produced by
calling |copy| for previous type constants.

@< Initial global variable definitions @>=

const type_expr unknown_type; // uses default constructor
 type_expr void_type(type_list_ptr(NULL));
 type_expr int_type(integral_type);
 type_expr bool_type(boolean_type);
const type_expr row_of_type(copy(unknown_type));
const type_expr gen_func_type(copy(unknown_type),copy(unknown_type));

@ We shall also need a tuple pattern with any number of unknown components.

@< Initial declarations of exported functions @>=
type_ptr unknown_tuple(size_t n);

@~This pattern is built up in a simple loop.

@< Initial function definitions @>=
type_ptr unknown_tuple(size_t n)
{ type_list_ptr tl(NULL);
  while (n-->0) tl=make_type_list(copy(unknown_type),tl);
  return make_tuple_type(tl);
}

@*1 Dynamically typed values.
%
Now we shall consider runtime values. As we mentioned before, the interpreter
must access values via generic pointers in order to be able to manipulate them
regardless of their types, which could be arbitrarily complicated. We could
either use void pointers to represent generic values and cast them when
necessary, or use inheritance and the dynamic cast feature of \Cpp. We choose
the second option, which is quite convenient to use, although this means that
in reality we have dynamic type information stored within the values, even
though that information had already been determined during type analysis. We
shall in fact use this information to double-check our type analysis at
runtime.

@< Includes needed in \.{types.h} @>=
#include <iostream> // needed for specification of |print| method below
#include <tr1/memory>

@~We start with a base class for values. There must be at least one virtual
function in the class, which can conveniently be a function for printing. This
allows the base class to be defined abstract (one cannot declare a destructor
purely virtual since it will always be called, after the destructor for a
derived class). The printing function will demonstrate the ease of using
dynamic typing via inheritance. It does not even require any dynamic casting,
but other operations on values will. Apart from |print| we define another
(purely) virtual method, |clone|, which allows making a copy of a runtime
value of any type derived from |value_base|.

The method |name| is useful in reporting logic errors from template functions,
notably the failure of a value to be of the predicted type. Since the template
function may know the type (via a template argument) but need not have any
object of the type at hand, we define |name| as a |static| rather than
|virtual| method. We forbid assignment of |value_base| objects, since they
should always be handled by reference; the base class is abstract anyway, but
this ensures us that for no derived class an implicitly defined assignment
operator is accidentally invoked. Copy constructors will in fact be defined
for all derived types, as they are needed to implement the |clone| method;
these will be |private| or |protected| as well, so as to forbid accidental use
elsewhere, but they don't copy-construct the |value_base| base object (rather
they default-construct it). We can then declare the |value_base| copy
constructor |private| so that in case of accidental omission the use of a
synthesised constructor will be caught here as well.

As mentioned values are always handled via pointers. We a raw pointer type
|value|, an auto-pointer |owned_value| (which cannot be stored in STL
containers), and a shared smart pointer |shared_value| (which by contrast can
be stored in STL containers).

@< Initial type definitions @>=
struct value_base
{ value_base() @+ {};
  virtual ~value_base() @+ {};
  virtual void print(std::ostream& out) const =0;
  virtual value_base* clone() const =0;
  static const char* name(); // just a model; this instance remains undefined
private: //copying and assignment forbidden
  value_base& operator=(const value_base& x);
  value_base(const value_base& x);
};
@)
typedef value_base* value;
typedef std::auto_ptr<value_base> owned_value;
typedef std::tr1::shared_ptr<value_base> shared_value;

@ We can already make sure that the operator~`|<<|' will do the right thing
for any of our values.

@< Initial declarations of exported functions @>=
std::ostream& operator<< (std::ostream& out, const value_base& v);

@~The operator~`|<<|' handles calling the virtual |print| method of the actual
value as usual.

@< Initial function definitions @>=
std::ostream& operator<< (std::ostream& out, const value_base& v)
{@; v.print(out); return out; }

@ Often we know what variant a |value| object takes, based on the type
analysis. We can convert to that type using a |dynamic_cast|, but at such
moments we wish to throw a |logic_error| in case our type prediction was
wrong. To avoid having such casts and |throw| statements all over the place,
we define a template function to do the casting and throwing. It is defined at
the level of ordinary pointers, and it is not intended for use where the
caller assumes ownership of the result; original pointer is assumed to retain
ownership as long as the result of this call survives, and in particular that
result should not be converted to a smart pointer, lest double deletion would
ensue.

@< Initial template and inline function definitions @>=
template <typename D> // |D| is a type derived from |value_base|
 D* force(value v) throw(std::logic_error)
{ D* p=dynamic_cast<D*>(v);
  if (p==NULL) throw
    std::logic_error(std::string("forced value is no ")+D::name());
  return p;
}

@ Here we define a first type derived from |value_base|, namely the type for
``row of'' types. They are implemented using vectors from the standard
template library.

@< Includes needed in \.{types.h} @>=
#include <vector>
#include <cassert>

@~Since the actual values accessed will be of types derived
from |value_base|, we must pass through a level of indirection, so we have a
vector of pointers. We define these pointers to be |shared_value| pointers, so
that the row takes (shared) ownership of its components without needing a
explicit destructor. This as the additional advantage over explicit ownership
management that the copy constructor, needed for the |clone| method, can
safely just copy-construct the vector of pointers: a possible exception thrown
during the copy is guaranteed to clean up any pointers present in the vector.

Of course ownership of pointers to |row_value| objects also needs to be
managed, which could be either by an auto-pointer |row_ptr| (if the pointer is
known to be unshared) or by a shared pointer |shared_row|.

@< Initial type definitions @>=
struct row_value : public value_base
{ std::vector<shared_value> val;
@)
  explicit row_value(size_t n) : val(n) @+{} // start with |n| empty pointers
  void print(std::ostream& out) const;
  size_t length() const @+{@; return val.size(); }
  row_value* clone() const @+{@; return new row_value(*this); }
  static const char* name() @+{@; return "row value"; }
protected:
  row_value(const row_value& v) : val(v.val) @+{}
};
@)
typedef std::auto_ptr<row_value> row_ptr;
typedef std::tr1::shared_ptr<row_value> shared_row;

@ So here is the first occasion where we shall use virtual functions. For the
moment the output routine performs an immediate recursion; later we shall try
to make this more elegant by computing the width needed to output component
values, and adapt the formatting to that.

@< Initial function definitions @>=
void row_value::print(std::ostream& out) const
{ if (val.empty()) out << "[]";
  else
  { out << '[';
    std::vector<shared_value>::const_iterator p=val.begin();
    do {@; (*p)->print(out); ++p; out << (p==val.end() ? ']' : ','); }
    while (p!=val.end());
  }
}


@ Since we use dynamically typed values internally, we can collect the
components of a tuple in a vector without problem. In fact we could reuse the
type |row_value| to hold the components of a tuple, if it weren't for the fact
that it would then print with brackets. Therefore we trivially derive a new
class from |row_value|.

@< Initial type definitions @>=
struct tuple_value : public row_value
{ tuple_value(size_t n) : row_value(n) @+{}
  tuple_value* clone() const @+{@; return new tuple_value(*this); }
  void print(std::ostream& out) const;
  static const char* name() @+{@; return "tuple value"; }
private:
  tuple_value(const tuple_value& v) : row_value(v)@+{}
 // copy constructor; used by |clone|
};
@)
typedef std::auto_ptr<tuple_value> tuple_ptr;
typedef std::tr1::shared_ptr<tuple_value> shared_tuple;

@ We just need to redefine the |print| method.
@< Initial function definitions @>=
void tuple_value::print(std::ostream& out) const
{ if (val.empty()) out << "()";
  else
  { out << '(';
    std::vector<shared_value>::const_iterator p=val.begin();
    do {@; (*p)->print(out); ++p; out << (p==val.end() ? ')' : ','); }
    while (p!=val.end());
  }
}

@ Here are functions that pack and unpack tuples from values on the stack;
they will be used by wrapper functions around functions from the Atlas
library, and in the case of |wrap_tuple| for the evaluation of tuple
expressions. The function |push_tuple_components| will be called by wrapper
functions that need the tuple components on the stack; the function call
|wrap_tuple(n)| inversely builds a tuple from $n$ components on the stack.


@< Initial declarations of exported functions @>=
void push_tuple_components();
void wrap_tuple(size_t n);

@~These functions use the same convention for stack order: the last tuple
component is on top of the stack. In |push_tuple_components| we start by
popping a value and checking it to be a tuple, which is done by
|get<tuple_value>| that will be defined later. The |shared_value| type takes
care of ownership; there is (at least) double shared ownership of the
components as the stack expands, but this is normal, and afterwards this
sharing disappears with |tuple|.

@< Initial function definitions @>=
void push_tuple_components()
{ shared_tuple tuple(get<tuple_value>());
  for (size_t i=0; i<tuple->length(); ++i)
    push_value(tuple->val[i]); // push component
}

@ We need no auto-pointer in |wrap_tuple|, as shrinking the stack will not
throw any exceptions.

@< Initial function definitions @>=
void wrap_tuple(size_t n)
{ shared_tuple result(new tuple_value(n));
  while (n-->0) // standard idiom; not |(--n>=0)|, since |n| is unsigned!
    result->val[n]=pop_value();
  push_value(result);
}

@*2 Representation of an evaluation context.
%
While evaluating user programs, values will be given to local identifiers such
as arguments of functions being called. The identification of identifiers is
determined during type analysis (static binding); it results for local
identifiers in a method to locate the associated value in the evaluation
context, which is formed by a stack of frames, each holding a vector of
values.
Frames are actually allocated on the heap, and their lifetimes do not follow
a stack regime unless a very limited use is made of user-defined functions
(never passing such a function as value out of the expression in which it was
defined), so it is better to just say they are linked lists of frames. A
singly linked list suffices, and by using shared pointers as links,
destruction of frames once inaccessible is automatic.

@< Initial type definitions @>=
class context;
typedef std::tr1::shared_ptr<context> context_ptr;
class context
{ const context_ptr next;
  std::vector<shared_value> frame;
  context(const context&); // contexts should not be copied, just shared
public:
  context(const context_ptr& n,
          const std::vector<shared_value>& f) : next(n), frame(f) @+{}
  shared_value& elem(size_t i,size_t j);
};

@ The method |context::elem| descends the stack and then selects a value from
the proper frame.

@< Initial function def... @>=
shared_value& context::elem(size_t i, size_t j)
{
  context* p=this;
  for (; i-->0; p=p->next.get())
    assert(p->next.use_count()>0 and p->next.get()!=NULL);
  assert(j<p->frame.size());
  return p->frame[j];
}

@*1 Values representing type-checked expressions.
%
The parser is (currently) a \Cee-program that upon success returns a value of
type |expr| representing this parse tree. While analysing this expression for
type-correctness, it will be convenient to transform it into a value that can
be efficiently evaluated. This value will be a pointer to an object of one of
a number of classes directly derived from an empty base class (in the same way
as runtime values are pointers to an object of a class derived from
|value_base|), which classes have a virtual method |evaluate| that performs
the operation described by the expression. We shall now define these classes.

A fundamental choice is whether to make the result type of the |evaluate| type
equal to |value|. Although this would seem the natural choice, we prefer
to make its result |void|, and to handle all value-passing via an execution
stack; thus we hope to be able to avoid packing and unpacking of tuple values
in most cases when calling functions. In order to do that effectively, we
dynamically pass a parameter to the |evaluate| method telling whether the
result is expected to be ``expanded'' on the runtime stack in case it is of a
tuple type.

@< Initial type definitions @>=
struct expression_base
{ enum level @+{ no_value, single_value, multi_value };
@)
  expression_base() @+ {}
  virtual ~expression_base() @+ {}
  virtual void evaluate(level l) const =0;
  virtual void print(std::ostream& out) const =0;
@)
  void void_eval() const @+{@; evaluate(no_value); }
  void eval() const @+{@; evaluate(single_value); }
  void multi_eval() const @+{@; evaluate(multi_value); }
};
@)
typedef expression_base* expression;
typedef std::auto_ptr<expression_base> expression_ptr;
typedef std::tr1::shared_ptr<expression_base> shared_expression;

@ Like for values, we can assure right away that printing converted
expressions will work.

@< Initial declarations of exported functions @>=
inline std::ostream& operator<< (std::ostream& out, const expression_base& e)
{@; e.print(out); return out; }


@ Since our |evaluate| methods will put values on the execution stack, let us
declare it right away. We decide that values on the execution stack can be
shared with other values (for instance when the user subscripts a row, vector
or matrix bound to an identifier, it would be wasteful to duplicate that
entire structure just so that it can briefly reside on the execution stack),
whence we use |shared_value| smart pointers in the stack. This choice will
have consequences in many places in the evaluator, since once a value is
referred to by such a smart pointer, its ownership cannot be transferred to
any other regime; when strict ownership should be needed, the only option will
be to make a copy by calling |clone|.

@< Initial declarations of global variables @>=
extern std::vector<shared_value> execution_stack;

@~We define the stack as a static variable of this compilation unit; it is
initially empty. All usable built-in functions will be provided with a small
wrapper function that takes it values from the stack and places its results
there again. Parameters are placed on the stack in order, and should therefore
be popped from the stack in reverse order.

@< Initial global variable definitions @>=
std::vector<shared_value> execution_stack;

@ We shall define some inline functions to facilitate manipulating the stack.
The function |push_value| does what its name suggests. For exception safety it
takes either an auto-pointer or a shared pointer as argument; the former is
converted into the latter, in which case the |use_count| will become~$1$. For
convenience we make these template functions that accept a smart pointer to
any type derived from |value_base| (since a conversion of such pointers from
derived to base is not possible without a cast in a function argument
position). For even more convenience we also provide a variant taking an
ordinary pointer, so that expressions using |new| can be written without cast
in the argument of |push_value|. Since |push_value| has only one argument,
such use of does not compromise exception safety: nothing can throw between
the return of |new| and the conversion of its result into an auto-pointer.

@< Initial template and inline function definitions @>=
template<typename D> // |D| is a type derived from |value_base|
  inline void push_value(std::auto_ptr<D> v)
  @+{@; execution_stack.push_back(std::tr1::shared_ptr<D>(v)); }

template<typename D> // |D| is a type derived from |value_base|
  inline void push_value(std::tr1::shared_ptr<D> v)
  @+{@; execution_stack.push_back(v); }

inline void push_value(value_base* p) @+{@; push_value(owned_value(p)); }

@ There is a counterpart |pop_value| to |push_value|. Most often the result
pust be dynamically cast to the type they are known to have because we passed
the type checker; should the cast fail we shall throw a |std::logic_error|.
The template function |get| with explicitly provided type serves for this
purpose; it is very much like the template function |force|, but returns a
shared pointer (because the value on the stack might be shared).

@< Initial template and inline function definitions @>=

inline shared_value pop_value()
{@; shared_value arg=execution_stack.back();
  execution_stack.pop_back();
  return arg;
}
@)
template <typename D> // |D| is a type derived from |value_base|
 inline std::tr1::shared_ptr<D> get() throw(std::logic_error)
{ std::tr1::shared_ptr<D> p=std::tr1::dynamic_pointer_cast<D>(pop_value());
  if (p.get()==NULL)
    throw std::logic_error(std::string("Argument is no ")+D::name());
  return p;
}
@.Argument is no ...@>

@ The function |push_expanded| will be called when a stored value is retrieved
in a context where a tuple value might need expansion, depending on the |level
l@;|, but the value is not known to be of tuple type (since type information
is not retained in compiled expression values).

@< Function definitions @>=
void push_expanded(expression_base::level l, const shared_value& v)
{ if(l==expression_base::single_value)
    push_value(v);
  else if (l==expression_base::multi_value)
  { shared_tuple p = std::tr1::dynamic_pointer_cast<tuple_value>(v);
    if (p==NULL)
      push_value(v);
    else
      for (size_t i=0; i<p->length(); ++i)
        push_value(p->val[i]); // push components
  }
}

@ In many cases we will want to get unique access to an object, duplicating it
if necessary. The operation |uniquify| implements this; it was originally
introduced in order to implement component assignments efficiently (the code
for this will be given later) but is useful much more generally. In fact if we
just computed the value in question from a function call it is virtually
guaranteed to be unshared, but we shall call |uniquify| anyway to make clear
our (destructive) intentions. Here we also use it right away to provide a
variant template function |get_own| of |get|, which returns a privately owned
copy of the value from the stack, so that modifications can be made to it
without danger of altering shared instances. In spite of the uniqueness
guarantee, |get_own| must be declared to return a |shared_ptr| in order to
avoid having to call |clone|: there is no way to persuade a |shared_ptr| to
release its ownership, even if it happens to be the unique owner.


@< Initial template and inline function def... @>=
inline void uniquify(shared_value& v)
{@; if (not v.unique()) v=shared_value(v->clone()); }
@)
template <typename D> // |D| is a type derived from |value_base|
 std::tr1::shared_ptr<D> get_own() throw(std::logic_error)
{@; uniquify(execution_stack.back());
    return get<D>();
}

@*1 Implicit conversion of values between types.
%
When interfacing this generic interpreter with a concrete library such as that
of the Atlas of Lie Groups and Representations, a mechanism must be provided
to convert data in the interpreter (represented essentially as nested lists)
into the internal format of the library. For transparency of this mechanism we
have chosen to provide the conversion through implicit operations that are
accompanied by type changes; thus when the user enters a list of lists of
integers in a position where an integral matrix is required, the necessary
conversions are automatically inserted during type analysis. In fact we shall
put in place a general mechanism of automatic type conversions, which will for
instance also provide the inverse conversions where appropriate, and on some
occasions merely provides convenience to the user, for instance by allowing
integers in positions where a rational number is required.

@ We shall derive a single class |conversion| from |expression_base| to
represent any expression that is to be converted using one of the various
conversions. Which conversion is to be applied is determined by |type|, a
reference to a |conversion_info| structure containing a function pointer
|convert| that provides the actual conversion routine; this is easier than
deriving a plethora of classes differing only in their virtual method
|evaluate|, although it is slightly less efficient. The type |conv_f| of
conversion function pointer specifies no argument or return value, as we know
beforehand that the |shared_value| objects serving for this are to be found
and left on the runtime stack.

@< Initial type definitions @>=
struct conversion_info
{ typedef void (*conv_f)();
  conv_f convert;
  std::string name;
  conversion_info(const char* s,conv_f c) : convert(c),name(s) @+{}
};
@)
class conversion : public expression_base
{ const conversion_info& type;
  expression exp;
public:
  conversion(const conversion_info& t,expression_ptr e)
   :type(t),exp(e.release()) @+{}
  virtual ~conversion()@;{@; delete exp; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ The |evaluate| method for conversions dispatches to the |convert| member,
after evaluating |exp|. Although automatic conversions are only inserted when
the type analysis requires a non-empty result type, it is still in principle
possible that at run time this method is called with |l==no_value|, so we
cater for that. The case is so rare that we don't mind the inefficiency of
performing the conversion and then discarding the result; this will allow a
failing conversion to be signalled as an error in such cases.

@< Initial function def...@>=
void conversion::evaluate(level l) const
{@; exp->eval();
  (*type.convert)();
  if (l==no_value)
    execution_stack.pop_back();
}
@)
void conversion::print(std::ostream& out) const
@+{@; out << type.name << ':' << *exp; }

@*2 Coercion of types.
%
An important aspect of automatic conversions is that they can be applied in
situations where the result type and only part of the source type is known:
for instance one can be inserted when a list display occurs in a context
requiring a vector, because no row type equal the (primitive) vector type.
While this use requires particular consideration according to the syntactic
form of the expression (list display), there is also a simpler form of
automatic conversion that can be applied to a large variety of expressions
(identifiers, function calls, \dots) whenever they are found to have a
different type from what the context requires. The a function |coerce| will
try to insert an automatic conversion in such situations, if this can resolve
the type conflict. We present this mechanism first, since the table it employs
can then be re-used to handle the more subtle cases of automatic conversions.

The function |coerce| requires to fully determined types |from_type| and
|to_type|, and its final argument~|e| is a reference to the previously
converted expression. If a conversion of value of |from_type| to |to_type| is
available, then |coerce| will modify |e| by insertion of a |conversion| around
it; the return value of |coerce| indicates whether an applicable conversion
was found.

@< Initial declarations of exported functions @>=

bool coerce(const type_expr& from_type, const type_expr& to_type,
            expression_ptr& e);

@ The implementation of |coerce| will be determined by a simple table lookup.
The records in this table contain a |conversion_info| structure (in fact they
are derived from it) and in addition indications of the types converted from
and to. The table entries store these types by pointer, so that table entries
are assignable, and the table can be allowed to grow (which is a technical
necessity in order to allow other compilation units to contribute coercions).

@< Initial type definitions @>=
struct conversion_record : public conversion_info
{ const type_expr* from,* to; // non-owned pointers
  conversion_record (const type_expr& from_type,
                     const type_expr& to_type,
                     const char* s, conv_f c)
   : conversion_info(s,c), from(&from_type),to(&to_type) @+{}
};

@ Here is the lookup table. It is defined as vector so that other compilation
units can extend it; however it should not be extended once evaluation starts,
since |conversion| objects will store references to table entries, which would
become invalid in case of reallocation.

@< Initial global variable definitions @>=

std::vector<conversion_record> coerce_table;

@ The following function simplifies filling the coercion table; it is
externally callable.

@< Initial declarations of exported functions @>=
void coercion(const type_expr& from,
              const type_expr& to,
              const char* s, conversion_info::conv_f f);
@~The action is simply extending |coerce_table| with a new |conversion_record|.
@< Initial function def... @>=
void coercion(const type_expr& from,
              const type_expr& to,
              const char* s, conversion_info::conv_f f)
{@; coerce_table.push_back(conversion_record(from,to,s,f)); }

@ There is one coercion that is not stored in the lookup table, since it can
operate on any input type: the voiding coercion. It is necessary for instance
to allow a conditional expression that is intended for its side effects only
to have branches that evaluate to different types (including the possibility
of absent branches which will be taken to deliver an empty tuple): the
conditional expression will get \.{void} type to which the result of each of
its branches can be coerced.

In fact voiding is not dealt with using a |conversion| expression, since as we
have seen its |evaluate| method evaluates its argument using the |eval|
method, so with |level| parameter set to |single_value|. However for a voiding
coercion this level should be |no_value|. So we introduce a type derived from
|expression| whose main virtue is that its |evaluate| method sets the level to
|no_value| by calling |void_eval|.

@< Initial local type definitions @>=
class voiding : public expression_base
{ expression exp;
public:
  voiding(expression_ptr e) : exp(e.release()) @+{}
  virtual ~voiding()@;{@; delete exp; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ The |evaluate| method should not ignore its |level| argument completely:
when |l==single_value| an actual empty tuple should be produced, which
|wrap_tuple(0)| does.

@< Initial local function definitions @>=
void voiding::evaluate(level l) const
{@; exp->void_eval();
  if (l==single_value)
    wrap_tuple(0);
}
@)
void voiding::print(std::ostream& out) const
@+{@; out << "voided:" << *exp; }


@ The function |coerce| simply traverses the |coerce_table| looking for an
appropriate entry, and wraps |e| into a corresponding |conversion| it finds
one. Ownership of the expression pointed to by |e| is handled implicitly: it
is released during the construction of the |conversion|, and immediately
afterwards |reset| reclaims ownership of the pointer to that |conversion|.

@< Initial function definitions @>=
bool coerce(const type_expr& from_type, const type_expr& to_type,
	    expression_ptr& e)
{ for (conversion_record*
       it=&coerce_table[0]; it<&*coerce_table.end(); ++it)
    if (from_type==*it->from and to_type==*it->to)
    @/{@; e.reset(new conversion(*it,e));
      return true;
    }
  if (to_type==void_type)    {@; e.reset(new voiding(e));  return true; }
  return false;
}

@ List displays and loops produce a row of values of arbitrary (but identical)
type; when they occur in a context requiring a non-row type, we may be able to
find a coercion that reconciles the requirements. The following function finds
whether this is possible, and sets |components_type| to the type required for
the components; by using this mechanism the components themselves obtain a
context that may generate further conversions to obtain this type.

@< Initial declarations of exported functions @>=
conversion_record* row_coercion(const type_expr& final_type,
                                      type_expr& component_type);
@~The implementation is simply to look in |coerce_table| for conversions from
some row type.

@< Initial function def... @>=
conversion_record* row_coercion(const type_expr& final_type,
                                      type_expr& component_type)
{ for (conversion_record* it=&coerce_table[0]; it<&*coerce_table.end(); ++it)
    if (final_type==*it->to and it->from->kind==row_type)
    @/{@; component_type.specialise(*it->from->component_type); return it; }
  return NULL;
}

@*2 Proximity of types.
%
The above implicit conversions of types pose a limitation to the possibilities
of overloading operator symbols and function identifiers. If a symbol should
be overloaded for too closely related operand types, situations could occur in
which given operand expressions can be converted to either of the operand
types. This would either produce unpredictable behaviour, or necessitate a
complicated set of rules to determine which of the definitions of the symbol
is to be used (overloading resolution in \Cpp\ is a good example of such
complications). There are two ways to avoid the occurrence of difficulties by
restricting the rules of the language: either forbid type conversions in
arguments of overloaded symbols, or forbidding simultaneous definitions of
such symbols for too closely related types. (A mixture of both is also
conceivable, allowing only certain conversions and forbidding overloading
between types related by them; the language Algol~68 is a good example of an
approach along these lines). Forbidding all automatic type conversions in case
of overloading would defeat to a large extent the purpose of overloading,
namely as a convenience to the user; therefore we choose the latter solution
of forbidding overloading in certain cases.

Two types are too closely related to allow overloading between them if any
conceivable expression could be converted to either of them if the context so
specified (this implies that \.{void} should be forbidden as operand type,
which is not a great limitation for overloading). The relation we shall
implement is intended to be an equivalence relation not finer than the one
just described, which will actually be the case under the assumption of some
decency on the set of conversions in |coerce_table|; notably if some type can
be directly converted to two different types, then one of those two types
can also be converted (possibly indirectly) to the other. In practice we have
for instance that \.{mat} can be converted to either \.{[vec]} or
to \.{[[int]]}, which is all right since each \.{vec} entry of a  \.{[vec]}
value can be converted into a \.{[int]} (and vice versa).

@< Initial declarations of exported functions @>=
bool is_close (const type_expr& x, const type_expr& y);

@ So here is the (recursive) definition of the relation |is_close|. A
primitive type is in the relation |is_close| to another type only is it is
identical or if there is a direct conversion between the types. Two distinct
primitive types are in the relation |is_close| if and only if they are either
both row types or both tuple types with the same number of component types,
and if each pair of corresponding component types is in the relation
|is_close|.

The expression |dummy| may be prepended to by |coerce|, but is then abandoned
(and cleaned up).

@< Initial function definitions @>=
bool is_close (const type_expr& x, const type_expr& y)
{ expression_ptr dummy(NULL);
  if (x.kind==primitive_type or y.kind==primitive_type)
    return x.prim==y.prim or coerce(x,y,dummy) or coerce(y,x,dummy);
  if (x.kind!=y.kind)
    return false;
  if (x.kind==row_type)
    return is_close(*x.component_type,*y.component_type);
  if (x.kind!=tuple_type)
    return x==y; // non-aggregate types are only close if equal
  type_list l0=x.tuple, l1=y.tuple; // recurse for tuples, as in |operator==|
  while (l0!=NULL and l1!=NULL and is_close(l0->t,l1->t))
    {@; l0=l0->next; l1=l1->next; }
  return l0==NULL and l1==NULL; // lists must end simultaneously for success
}

@*1 Error values.
Before we describe evaluation of expressions we must realise that evaluation
can cause runtime errors. The evaluator may throw exceptions due to
inconsistency of our (rather than the user's) program, which are classified as
|std::logic_error|. It may also throw exceptions due to errors not caught by
the type checker in the user input, such as size mismatch in matrix
operations; in such cases it will throw |std::runtime_error|. These standard
exception types will be used without any type derivation.

@< Includes needed in \.{types.h} @>=
#include <stdexcept>

@ For errors detected before execution starts, we first derive a general
exception class |program_error| from |std::exception|; it represents any kind
of error of the user input determined by static analysis (for instance use of
undefined variables).

Although it does nothing explicitly (the string will be destructed anyway), we
must explicitly define a destructor for |program_error|: the automatically
generated one would lack the |throw()| specifier.

@< Initial type definitions @>=
class program_error : public std::exception
{ std::string message;
public:
  explicit program_error(const std::string& s) : message(s) @+{}
  virtual ~program_error() throw() @+{@;} // obligatory definition
  const char* what() const throw() @+{@; return message.c_str(); }
};

@ We derive from |program_error| an exception type |expr_error| that stores in
addition to the error message an expression to which the message applies. It
is declared a |struct|, as we leave it up to the |catch| code to incorporate
the offending expression in a message in addition to the one produced by
|what()|. For this type no virtual methods are defined at all, in particular
we do not need a virtual destructor: the (automatically generated) destruction
code that takes care of destructing |offender| is not called through a vtable
(which in fact is absent). It is not quite clear whether this explains that
the compiler does not complain here that the undeclared destructor has a too
loose (because absent) |throw| specification, while it would for
|program_error|.

@< Initial type definitions @>=
struct expr_error : public program_error
{ expr offender; // the subexpression causing a problem
@)
  expr_error (const expr& e,const std::string& s) throw()
    : program_error(s),offender(e) @+{}
};

@ We derive from |expr_error| an even more specific exception type
|type_error| that we shall throw when an expression fails to have the right
type. In addition to the offending subexpression it stores pointers to the
types that failed to match. For the latter, the error value owns the types
pointed to, so the caller should relinquish ownership of, or copy (in case it
did not own), the types passed when throwing the exception.

@< Initial type definitions @>=
struct type_error : public expr_error
{ type_p actual; @+
  type_p required; // the types that conflicted
@)
  type_error (const expr& e, type_ptr a, type_ptr r) throw() @/
    : expr_error(e,"Type error") @|
      ,actual(a.release()),required(r.release()) @+{}
  type_error(const type_error& e);
  ~type_error() throw() @+{@; delete actual; delete required; }
};

@ We do not plan to use the copy constructor for |type_error| since its use on
throwing can be avoided by good compilers (and in addition we shall always
catch by reference, although the compiler cannot know this at the point where
it throws). However the language requires that a copy constructor be
accessible for objects thrown, so we cannot make the copy constructor private.
This being the case, we cannot let it be synthesised by the compiler either,
because a shallow copy would lead to the types pointed to being destroyed
twice.

So we must laboriously define the proper way to do it, hoping that it will
never be used. This is not easy either, since two copies of type expressions
must be made, which can throw exceptions (let us hope the hypothetical call to
the copy constructor is made \emph{before} the |throw| of |type_error| is
executed, or else we will be dead due to interwoven exceptions\dots), and
during the second copy the pointer to the first copy must be held in an
auto-pointer to avoid a memory leak. This excludes any attempt to initialise
|actual| and |required| directly to their proper values, sigh! In the code
below, we could have said instead |actual=new type_expr(*e.actual)|,
thus avoiding the use of the anonymous auto-pointer produced by the second
call of |copy|, at the price of some loss of symmetry and readability.

@< Initial function definitions @>=
type_error::type_error(const type_error& e)
 : expr_error(e)
{@; type_ptr p=copy(*e.required);
  actual=copy(*e.actual).release(); required=p.release();
}

@*1 Enumeration of primitive types.
%
The precise list of primitive types is of minor importance for the design of
the evaluator, but they must be given so that the type |type_expr| is
complete. So we collect here the list of primitive type tags, to be completed
whenever one is added in a different compilation unit. Most of these
enumeration values are never directly used, but they must be present so that
the evaluator can be aware of the number of primitive types (via the final
enumeration value |nr_of_primitive_types|).

@< Other primitive type tags @>=
vector_type, matrix_type, rational_vector_type,
complex_lie_type_type , root_datum_type, inner_class_type, real_form_type,
dual_real_form_type, Cartan_class_type, @[@]


@~The following list must match that of the previous module, for proper
functioning of I/O.

@< Other primitive type names @>=
"vec", "mat", "ratvec",
"LieType","RootDatum", "InnerClass", "RealForm", "DualRealForm",
"CartanClass", @[@]

@* The evaluator.
%
We now come to the actual routines involved in analysis and evaluation.

@*1 Preliminaries.
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
const type_expr int_int_type(*make_type("(int,int)"));
  // copy and destroy original


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


@*1 Type-checking and conversion to executable form  of expressions.
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
rational number is expected. The function |coerce|, defined later, tests for
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
      if(type.specialise(int_type) or coerce(int_type,type,d))
        return d.release();
      else throw type_error(e,copy(int_type),copy(type));
    }
   case string_denotation:
    { expression_ptr d@|(new denotation
        (shared_value(new string_value(e.e.str_denotation_variant))));
      if (type.specialise(str_type) or coerce(str_type,type,d))
        return d.release();
      else throw type_error(e,copy(str_type),copy(type));
    }
   case boolean_denotation:
    { expression_ptr d@|(new denotation
          (shared_value(new bool_value(e.e.int_denotation_variant))));
      if (type.specialise(bool_type) or coerce(bool_type,type,d))
        return d.release();
      else throw type_error(e,copy(bool_type),copy(type));
    }
   @\@< Other cases for type-checking and converting expression~|e| against
   |type|, all of which either |return| or |throw| a |type_error| @>
 }
 return NULL; // keep compiler happy
}

@*1 List displays.
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

@*1 Tuples.
%
Besides types for uniform sequences we shall need types for non-uniform
sequences of values, usually short and with a given sequence of types for
their components; their most obvious and simple use is for argument lists of
functions. Without using these tuple types one might define function types
taking multiple arguments, but then a function could still only yield a single
value. With tuple types multiple results can be produced, and there will be no
need to explicitly cater for functions with multiple arguments.

@*2 Tuple displays.
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
break;

@*2 Evaluating tuple displays.
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

@*1 Array subscription.
%
While we have seen expressions to build list, and vectors and matrices out of
them, we so far are not able to access their components once they are
constructed. To that end we shall now introduce operations to index such
values. We allow subscription of rows, but also of vectors and matrices. Since
after type analysis we know which of the cases applies, we define several
classes. These differ mostly by their |evaluate| method, so we first derive an
intermediate class from |expression_base|, and derive the others from it. This
class also serves to host an enumeration type that will serve later.

@< Type definitions @>=
struct subscr_base : public expression_base
{ enum sub_type @+{ row_entry, vector_entry, matrix_entry, matrix_column };
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
  if (type.specialise(subscr_type) or coerce(subscr_type,type,subscr))
    return subscr.release(); // and convert (derived|->|base) to |expression|
  else throw type_error(e,copy(subscr_type),copy(type));
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


@*1 Identifiers.
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

@*2 The global identifier table.
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
  type_p type_of(Hash_table::id_type id) const; // lookup
  shared_value value_of(Hash_table::id_type id) const; // lookup
  shared_share address_of(Hash_table::id_type id); // locate
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

@*2 Global identifiers.
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

@*2 Local identifiers.
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

@*1 Function calls.
%
One of the most basic tasks of the evaluator is to allow function calls, which
may involve either buit-in or user-defined functions. We start with
introducing a type for representing function calls after type checking.

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

@*2 Type-checking function calls.
%
The function in a call, whether bound to an identifier or given by some other
type of expression, determines its own type; once this is known, its argument
and result types can be used to help converting the argument expression and
the call expression itself. Therefore when encountering a call in
|convert_expr|, we first get the type of the expression in the function
position, requiring only that it be a function type, then type-check and
convert the argument expression using the obtained result type, and build a
converted function call~|call|. Finally we test if the required type matches
the return type (in which case we simply return~|call|), or if the return type
can be coerced to it (in which case we return |call| as transformed by
|coerce|); if neither is possible we throw a~|type_error|.

@< Other cases for type-checking and converting... @>=
case function_call:
{ type_ptr f_type=copy(gen_func_type); // start with generic function type
  expression_ptr fun(convert_expr(e.e.call_variant->fun,*f_type));
  expression_ptr arg
    (convert_expr(e.e.call_variant->arg,f_type->func->arg_type));
  expression_ptr call (new call_expression(fun,arg));
  if (type.specialise(f_type->func->result_type) or
      coerce(f_type->func->result_type,type,call))
    return call.release();
  else throw type_error(e,copy(f_type->func->result_type),copy(type));
}

@*2 Evaluating function calls.
%
The evaluation of the call of a built-in function executes a ``wrapper
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
  builtin_value(wrapper_function v,const char* n)
  : val(v), print_name(n) @+ {}
  virtual void print(std::ostream& out) const
  @+{@; out << ':' << print_name << ':'; }
  builtin_value* clone() const @+{@; return new builtin_value(*this); }
  static const char* name() @+{@; return "built-in function"; }
private:
  builtin_value(const builtin_value& v)
  : val(v.val), print_name(v.print_name)
  @+{}
};
@)
typedef std::auto_ptr<builtin_value> builtin_ptr;
typedef std::tr1::shared_ptr<builtin_value> shared_builtin;

@ To evaluate a |call_expression| object we evaluate the function, and then
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
  argument->evaluate(f==0 ? single_value : multi_value);
  try
  { if (f==NULL)
      @< Call user-defined function |fun| with argument on |execution_stack| @>
    else // built-in functions
      f->val(l); // call the wrapper function, leaving result on the stack
  }
  @< Catch-block for exceptions thrown withing function calls @>
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

@< Catch-block for exceptions thrown withing function calls @>=
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

@*1 Let-expressions.
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

@*1 User-defined functions.
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

@*1 Control structures.
%
We shall now introduce conventional control structures, which must of course
be part of any serious programming language; yet they were implemented only
after plenty of other language elements were in place, such as
let-expressions, functions, rows and selection form them, implicit
conversions.

@*2 Conditional expressions.
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
break;

@ Evaluating a conditional expression ends up evaluating either the
then-branch or the else-branch.

@< Function definitions @>=
void conditional_expression::evaluate(level l) const
{ condition->eval();
  if (get<bool_value>()->val)
   then_branch->evaluate(l);
  else else_branch->evaluate(l);
}

@*2 While loops.
%
Next we consider |while| loops, which have three parts (the final one is
optional; if absent it will be a null pointer).

@< Type def... @>=
struct while_expression : public expression_base
{ expression condition, body, next_part;
@)
  while_expression(expression_ptr c,expression_ptr b, expression_ptr n)
   : condition(c.release()),body(b.release()), next_part(n.release())
  @+{}
  virtual ~while_expression() @+
  {@; delete condition; delete body; delete next_part; }
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ Printing a |while| expression is straightforward, taking care that
|next_part| could be null.

@< Function definitions @>=
void while_expression::print(std::ostream& out) const
{ out << " while " <<  *condition << " do " << *body;
  if (next_part!=NULL)
    out << " next " << *next_part;
  out << " od ";
}

@ Type checking is a bit more complicated for |while| loops that for
conditional expressions, because a row result must be produced. If the context
requires void type, we shall require the same for the body, knowing that
generation of a row value will be unconditionally suppressed in these cases
anyway. In all other cases we proceed for the body expression as for the
components of a row display (except that there is only one expression in this
case). The next-part, if present, is always converted requiring void type.

@< Other cases for type-checking and converting... @>=
case while_expr:
  { w_loop w=e.e.while_variant;
    expression_ptr c (convert_expr(w->condition,bool_type));
    expression_ptr n
     (w->next_part.kind==tuple_display and w->next_part.e.sublist==NULL
     ? NULL @|
     : convert_expr(e.e.while_variant->next_part,void_type));
    if (type==void_type or type.specialise(row_of_type))
    { expression_ptr b
       (convert_expr(w->body, @|
                     type==void_type ? void_type :*type.component_type));
      @/return new while_expression(c,b,n);
    }
    else
    @< If |type| can be converted from some row-of type, check |w->body|
       against its component type, construct the |while_expression|, and apply
       the appropriate conversion function to it; otherwise |throw| a
       |type_error| @>
  }
break;

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
  expression_ptr loop(new while_expression(c,b,n));
  return new conversion(*conv,loop);
}


@ Of course evaluating is what most distinguishes loops from conditionals.

@< Function definitions @>=
void while_expression::evaluate(level l) const
{ if (l==no_value)
    while (condition->eval(),get<bool_value>()->val)
    { body->void_eval();
      if (next_part!=NULL)
        next_part->void_eval();
    }
  else
  { row_ptr result (new row_value(0));
    while (condition->eval(),get<bool_value>()->val)
    { body->eval(); result->val.push_back(pop_value());
      if (next_part!=NULL)
        next_part->void_eval();
    }
    push_value(result);
  }
}

@*2 For loops.
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
break;

@ For evaluating |for| loops we must take care to interpret the |kind| field
when selecting a component from the in-part. Because of differences in the thy
of |in_val|, some code must be duplicated, which we do as mush as possible by
sharing a module between the various loop bodies.

We can start evaluating the |in_part| regardless of |kind|, but for deducing
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
  if (kind==subscr_base::row_entry)
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
  else if (kind==subscr_base::vector_entry)
  { shared_vector in_val = get<vector_value>();
    size_t n=in_val->val.size();
    if (l!=no_value)
      result = row_ptr(new row_value(n));
    for (size_t i=0; unsigned(i)<n; ++i,loop_frame.clear())
    { loop_var->val[1].reset(new int_value(in_val->val[i]));
      @< Set |loop_var->val[0]| to |i|,... @>
    }
  }
  else if (kind==subscr_base::matrix_column)
  { shared_matrix in_val = get<matrix_value>();
    size_t n=in_val->val.numColumns();
    if (l!=no_value)
      result = row_ptr(new row_value(n));
    for (size_t i=0; unsigned(i)<n; ++i,loop_frame.clear())
    { loop_var->val[1].reset(new vector_value(in_val->val.column(i)));
      @< Set |loop_var->val[0]| to |i|,... @>
    }
  }
  execution_context = saved_context;
  if (l!=no_value)
    push_value(result);
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

@*2 Counted loops.
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
    (c->bound.kind!=tuple_display or c->bound.e.sublist!=NULL
    ? convert_expr(c->bound,int_type)
    : new denotation(zero)
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
break;

@ Executing a loop is a simple variation of what we have seen before for
|while| and |for| loops.

@< Function definitions @>=
void inc_for_expression::evaluate(level l) const
{ int b=(bound->eval(),get<int_value>()->val);
  int c=(count->eval(),get<int_value>()->val);

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

@*1 Casts.
%
Casts are very simple to process; they do not need any |expression| type to
represent them.

@< Other cases for type-checking and converting... @>=
case cast_expr:
{ cast c=e.e.cast_variant;
  type_expr& ctype=*static_cast<type_p>(c->type);
  expression_ptr p(convert_expr(c->exp,ctype));
  if (type.specialise(ctype) or coerce(ctype,type,p))
    return p.release();
  else throw type_error(e,copy(ctype),copy(type));
}
break;


@*1 Assignments.
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
  if (type.specialise(*it) or coerce(*it,type,assign))
    return assign.release();
  else throw type_error(e,copy(*it),copy(type));
}
break;

@*2 Component assignments.
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
  if (subscr_base::indexable(*aggr_t,ind_t,comp_t,kind))
  { expression_ptr r(convert_expr(rhs,comp_t));
    if (is_local)
      assign.reset(new local_component_assignment(aggr,i,d,o,r,kind));
    else
      assign.reset(new global_component_assignment(aggr,i,r,kind));
  }
  else
  { std::ostringstream o;
    o << "Cannot subscript " << *aggr_t << @| " value with index of type "
      << ind_t << " in assignment";
    throw expr_error(e,o.str());
  }

  if (type.specialise(*aggr_t) or coerce(*aggr_t,type,assign))
    return assign.release();
  else throw type_error(e,copy(*aggr_t),copy(type));
}
break;

@*1 Sequence expressions.
%
Since sequences are probably short on average, we use a chained
representation after type analysis, just like before.

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

@ To print a sequence, we print the two expressions separated by a semicolon.

@< Function definitions @>=
void seq_expression::print(std::ostream& out) const
{@; out << *first << ';' << *last; }

@ Evaluating a sequence expression evaluates the |first| for side effects
only, and then the |last| expression.

@< Function def... @>=
void seq_expression::evaluate(level l) const
{@; first->void_eval(); last->evaluate(l); }

@ It remains to type-check and convert sequence expressions, which is easy.

@< Other cases for type-checking and converting... @>=
case seq_expr:
{ sequence seq=e.e.sequence_variant;
  type_expr voided;
  expression_ptr first(convert_expr(seq->first,voided));
  expression_ptr last(convert_expr(seq->last,type));
  return new seq_expression(first,last);
}
break;

@*1 Invoking the type checker.
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

@*1 Primitive types for vectors and matrices.
%
The interpreter distinguishes its own types like \.{[int]} ``row of integer''
from similar built-in types of the library, like \.{vec} ``vector'', which it
will consider to be primitive types. In fact a value of type ``vector''
represents an object of the Atlas type |atlas::latticetypes::Weight|, which
stands for |atlas::matrix::Vector<int>|, and similarly other primitive types
will stand for other Atlas types. In choosing the short names, we choose to
hide the specific mathematical meaning that was implied in the names (like
|Weight|) these types have in the Atlas software. We believe that a more
extensive name might be more confusing that helpful to users; besides, the
interpretation of the values is not entirely fixed (vectors are used for
coweights and (co)roots as well as for weights, and matrices could denote
either a basis or an automorphism of a lattice).

By including the file~\.{latticetypes.h} we know about the basic vector and
matrix types.

@< Includes needed in the header file @>=
#include "latticetypes.h"

@ We start with deriving |vector_value| from |value_base|. In its constructor,
the argument is a reference to |latticetypes::CoeffList|, which stands for
|std::vector<latticetypes::LatticeCoeff>|, from which |latticetypes::Weight|
is derived (without adding data members), and since a constructor for the
latter from the former is defined, we can do with just one constructor for
|vector_value|.

@< Type definitions @>=

struct vector_value : public value_base
{ latticetypes::Weight val;
@)
  explicit vector_value(const latticetypes::CoeffList& v) : val(v) @+ {}
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
{ latticetypes::LatticeMatrix val;
@)
  explicit matrix_value(const latticetypes::LatticeMatrix& v) : val(v) @+ {}
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
  rational_vector_value
    (const latticetypes::LatticeElt& v,latticetypes::LatticeCoeff d)
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
  {@;  out << "The unique " << k << 'x' << l << " matrix"; return; }
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

@*2 Implementing some conversion functions.
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
latticetypes::Weight row_to_weight(const row_value& r)
{ latticetypes::Weight result(r.val.size());
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
extended will null entries to make a rectangular shape for the matrix.

@< Local function def... @>=
void matrix_convert()
{ shared_row r(get<row_value>());
@/latticetypes::WeightList column_list;
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
  push_value(new matrix_value(latticetypes::LatticeMatrix(column_list)));
}

@ All that remains is to initialise the |coerce_table|.
@< Initialise evaluator @>=
coercion(row_of_int_type, vec_type, "V", vector_convert); @/
coercion(row_of_vec_type,mat_type, "M", matrix_convert);

@ Here are conversions involving rational numbers and vectors.

@< Local function def... @>=
void rational_convert()
{@; shared_int i = get<int_value>();
    push_value(new rat_value(arithmetic::Rational(i->val)));
}
@)
void ratvec_convert()
{ shared_row r = get <row_value>();
  latticetypes::LatticeElt numer(r->val.size()),denom(r->val.size());
  unsigned int d=1;
  for (size_t i=0; i<r->val.size(); ++i)
  { arithmetic::Rational frac = force<rat_value>(r->val[i].get())->val;
    numer[i]=frac.numerator();
    denom[i]=frac.denominator();
    d=arithmetic::lcm(d,denom[i]);
  }
  for (size_t i=0; i<r->val.size(); ++i)
    numer[i]*= d/denom[i];

  push_value(new rational_vector_value(latticetypes::RatWeight(numer,d)));
}
@)
void rat_list_convert()
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
@/latticetypes::WeightList column_list;
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
  push_value(new matrix_value(latticetypes::LatticeMatrix(column_list)));

}

@ For the ``externalising'' conversions, it will be handy to have a basic
function |weight_to_row| that performs more or less the inverse transformation
of |row_to_weight|, but rather than returning a |row_value| it returns a
|row_ptr| pointing to it.

@< Local function def... @>=
row_ptr weight_to_row(const latticetypes::Weight& v)
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
{ shared_matrix m(get<matrix_value>());
  row_ptr result(new row_value(m->val.numColumns()));
  for(size_t i=0; i<m->val.numColumns(); ++i)
    result->val[i]=shared_value(new vector_value(m->val.column(i)));
  push_value(result);
}
@)
void int_list_list_convert()
{ shared_matrix m(get<matrix_value>());
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

@*1 Wrapper functions.
%
We have not defined any wrapper functions yet, and therefore have nothing in
the |global_id_table|. The following function will greatly facilitate the
repetitive task of installing wrapper functions.

@< Template and inline... @>=
inline void install_function
 (wrapper_function f,const char*name, const char* type)
{ shared_value val(new builtin_value(f,name));
  global_id_table->add
    (main_hash_table->match_literal(name),val,make_type(type));
}

@ Our first built-in functions implement with integer arithmetic. Arithmetic
operators are implemented by wrapper functions with two integer arguments.
Since arguments top built-in functions are evaluated with |level| parameter
|multi_value|, two separate value will be produced on the stack. Note that
these are pulled from the stack in reverse order, which is important for the
non-commutative operations like `|-|' and `|/|'. Since values are shared, we
must allocate new value objects for the results.

@< Local function definitions @>=

void plus_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new int_value(i->val+j->val));
}
@)
void minus_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new int_value(i->val-j->val));
}
@)
void times_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new int_value(i->val*j->val));
}
@)
void divide_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (j->val==0) throw std::runtime_error("Division by zero");
  if (l!=expression_base::no_value)
    push_value(new int_value(i->val/j->val));
}

@ We also define a remainder operation |modulo|, a combined
quotient-and-remainder operation |divmod|, and unary subtraction. Finally we
define an exact devision operator that constructs a rational number.

@< Local function definitions @>=
void modulo_wrapper(expression_base::level l)
{ shared_int  j=get<int_value>(); shared_int i=get<int_value>();
  if (j->val==0) throw std::runtime_error("Modulo zero");
  if (l!=expression_base::no_value)
    push_value(new int_value(i->val%j->val));
}
@)
void divmod_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (j->val==0) throw std::runtime_error("DivMod by zero");
  if (l!=expression_base::no_value)
  { push_value(new int_value(i->val/j->val));
    push_value(new int_value(i->val%j->val));
    if (l==expression_base::single_value)
      wrap_tuple(2);
  }
}
@)
void unary_minus_wrapper(expression_base::level l)
{@; shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new int_value(-i->val)); }
@)
void fraction_wrapper(expression_base::level l)
{ shared_int d=get<int_value>(); shared_int n=get<int_value>();
  if (d->val==0) throw std::runtime_error("fraction with zero denominator");
  if (l!=expression_base::no_value)
    push_value(new rat_value(arithmetic::Rational(n->val,d->val)));
}

@ Relational operators are of the same flavour.
@< Local function definitions @>=

void eq_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val==j->val));
}
@)
void neq_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val!=j->val));
}
@)
void less_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val<j->val));
}
@)
void lesseq_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val<=j->val));
}
@)
void greater_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val>j->val));
}
@)
void greatereq_wrapper(expression_base::level l)
{ shared_int j=get<int_value>(); shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
    push_value(new bool_value(i->val>=j->val));
}

@ To have some operations on strings, we define a function for concatenating
them, and one for converting integers to their string representation (of
course this remains a very limited repertoire).
@< Local function definitions @>=
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

@ Here is a simple function that outputs a string without its quotes, but with
a terminating newline. This is the first place in this file where we produce
output to a file. In general, rather than writing directly to |std::cout|, we
shall pass via a pointer whose |output_stream| value is maintained in the main
program, so that redirecting output to a different stream can be easily
implemented. Since this is a wrapper function (hence without arguments) there
is no other way to convey the output stream to be used than via a dedicated
global variable.

@< Local function definitions @>=
void print_wrapper(expression_base::level l)
{ shared_string s=get<string_value>();
  *output_stream << s->val << std::endl;
  if (l==expression_base::single_value)
    wrap_tuple(0); // don't forget to return a value if asked for
}


@ We now define a few functions, to really exercise something, even if it is
modest, from the Atlas library. These wrapper function are not really to be
considered part of the interpreter, but a first step to its interface with the
Atlas library, which is developed in much more detail in the compilation
unit \.{built-in-types}. In fact we shall make some of these wrapper functions
externally callable, so they can be directly used from that compilation unit.
First of all we have the identity matrix and matrix transposition.

@< Declarations of exported functions @>=
void id_mat_wrapper (expression_base::level);
void transpose_mat_wrapper (expression_base::level);

@ In |id_mat_wrapper| we create a |matrix_value| around an empty
|LatticeMatrix| rather than build a filled matrix object first. If a
constructor or function producing an identity matrix were available we could
have used that, but there is (currently) no such constructor or function; the
main reason seems to be that |matrix::Matrix| is a class template, and
specifying the desired entry type would be somewhat awkward.

Since in general built-in functions may throw exceptions (even for such simple
operations as |transposed|) we hold the pointers to local values in smart
pointers; for values popped from the stack this would in fact be hard to avoid.

@< Function definitions @>=
void id_mat_wrapper(expression_base::level l)
{ shared_int i=get<int_value>();
  if (l!=expression_base::no_value)
  { matrix_ptr m (new matrix_value(latticetypes::LatticeMatrix()));
    matrix::identityMatrix(m->val,std::abs(i->val)); push_value(m);
  }
}
@) void
transpose_mat_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l!=expression_base::no_value)
    push_value(new matrix_value(m->val.transposed()));
}

@ We also define |diagonal_wrapper|, a slight generalisation of
|id_mat_wrapper| that produces a diagonal matrix from a vector. The function
|vector_div_wrapper| produces a rational vector.

@< Local function def... @>=
void diagonal_wrapper(expression_base::level l)
{ shared_vector d=get<vector_value>();
  if (l==expression_base::no_value)
    return;
  size_t n=d->val.size();
  matrix_ptr m (new matrix_value(latticetypes::LatticeMatrix(n,n,0)));
  for (size_t i=0; i<n; ++i)
    m->val(i,i)=d->val[i];
  push_value(m);
}

void vector_div_wrapper(expression_base::level l)
{ shared_int n=get<int_value>();
  shared_vector v=get<vector_value>();
  if (l!=expression_base::no_value)
    push_value(new rational_vector_value(v->val,n->val));
}

@ Now the product of a matrix and a vector or matrix. The first of these was
in fact our first function with more than one argument (arithmetic on integer
constants was done inside the parser at that time). We make them callable from
other compilation units.

@< Declarations of exported functions @>=
void mv_prod_wrapper (expression_base::level);
void mm_prod_wrapper (expression_base::level);

@ In |mv_prod_wrapper| we use the matrix method |apply| which requires its
output vector to be already of the proper size, but which nevertheless could
throw an exception (since it copies its input argument, for in case the output
argument should coincide with it). Therefore we use an auto-pointer for the
output vector~|w|. In |mm_prod_wrapper|, the method |operator*=| does its own
resizing (and may also throw an exception).

For wrapper functions with multiple arguments, we must always remember that
they are to be popped from the stack in reverse order, but in the case of
|mm_prod_wrapper| this is particularly important.

@< Function definitions @>=
void mv_prod_wrapper(expression_base::level l)
{ shared_vector v=get<vector_value>();
  shared_matrix m=get<matrix_value>();
  if (m->val.numColumns()!=v->val.size())
    throw std::runtime_error(std::string("Size mismatch ")@|
     + str(m->val.numColumns()) + ":" + str(v->val.size()) + " in mv_prod");
  if (l!=expression_base::no_value)
    push_value(new vector_value(m->val.apply(v->val)));
}
@)
void mm_prod_wrapper(expression_base::level l)
{ shared_matrix rf=get<matrix_value>(); // right factor
  shared_matrix lf=get<matrix_value>(); // left factor
  if (lf->val.numColumns()!=rf->val.numRows())
  { std::ostringstream s;
    s<< "Size mismatch " << lf->val.numColumns() << ":" << rf->val.numRows()
    @| << " in mm_prod";
    throw std::runtime_error(s.str());
  }
  if (l!=expression_base::no_value)
    push_value(new matrix_value(lf->val*rf->val));
}

@ As a last example, here is the Smith normal form algorithm. We provide both
the invariant factors and the rewritten basis on which the normal for is
assumed, as separate functions, and to illustrate the possibilities of tuples,
the two combined into a single function.

@h "smithnormal.h"

@< Local function definitions @>=
void invfact_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l==expression_base::no_value)
    return;
  size_t nr=m->val.numRows();
  latticetypes::WeightList b; @+ matrix::initBasis(b,nr);
  vector_ptr inv_factors (new vector_value(latticetypes::Weight(0)));
  smithnormal::smithNormal(inv_factors->val,b.begin(),m->val);
  push_value(inv_factors);
}
@)
void Smith_basis_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l==expression_base::no_value)
    return;
  size_t nr=m->val.numRows();
  latticetypes::WeightList b; @+ matrix::initBasis(b,nr);
  latticetypes::Weight inv_factors(0);
  smithnormal::smithNormal(inv_factors,b.begin(),m->val);
  push_value(new matrix_value(latticetypes::LatticeMatrix(b)));
}
@)

void Smith_wrapper(expression_base::level l)
{ shared_matrix m=get<matrix_value>();
  if (l==expression_base::no_value)
    return;
  size_t nr=m->val.numRows();
  latticetypes::WeightList b; @+ matrix::initBasis(b,nr);
  vector_ptr inv_factors (new vector_value(latticetypes::Weight(0)));
  smithnormal::smithNormal(inv_factors->val,b.begin(),m->val);
@/push_value(new matrix_value(latticetypes::LatticeMatrix(b)));
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
install_function(unary_minus_wrapper,"-u","(int->int)");
install_function(fraction_wrapper,"/","(int,int->rat)");
install_function(eq_wrapper,"=","(int,int->bool)");
install_function(neq_wrapper,"!=","(int,int->bool)");
install_function(less_wrapper,"<","(int,int->bool)");
install_function(lesseq_wrapper,"<=","(int,int->bool)");
install_function(greater_wrapper,">","(int,int->bool)");
install_function(greatereq_wrapper,">=","(int,int->bool)");
install_function(concatenate_wrapper,"concatenate","(string,string->string)");
install_function(int_format_wrapper,"int_format","(int->string)");
install_function(print_wrapper,"print","(string->)");
install_function(vector_div_wrapper,"vec_div","(vec,int->ratvec)");
install_function(id_mat_wrapper,"id_mat","(int->mat)");
install_function(transpose_mat_wrapper,"transpose_mat","(mat->mat)");
install_function(diagonal_wrapper,"diagonal_mat","(vec->mat)");
install_function(mv_prod_wrapper,"mv_prod","(mat,vec->vec)");
install_function(mm_prod_wrapper,"mm_prod","(mat,mat->mat)");
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
local one, so we take an |id_pat| as argument, and follow the logic for
type-analysis of a let-expression, and for evaluation that of binding
identifiers in a user-define function. However we use |analyse_types| (which
reports errors) rather than calling |convert_expr| directly. To provide some
feedback to the user we report any types assigned, but not the values.

@< Function definitions @>=
extern "C"
void global_set_identifier(id_pat pat, expr rhs)
{ using namespace atlas::interpreter;
  size_t n_id=count_identifiers(pat);
  try
  { expression_ptr e;
    type_ptr t=analyse_types(rhs,e);
    if (not pattern_type(pat)->specialise(*t))
      @< Report that type of |rhs| does not have required structure,
         and |throw| @>

    bindings b(n_id);
    thread_bindings(pat,*t,b); // match identifiers and their future types

    std::vector<shared_value> v;
    v.reserve(n_id);
@/  e->eval();
    thread_components(pat,pop_value(),v);

    if (n_id>0)
      std::cout << "Identifier";
    for (size_t i=0; i<n_id; ++i)
    { std::cout << (i==0 ? n_id==1 ? " " : "s " : ", ") @|
                << main_hash_table->name_of(b[i].first) << ": "
                << *b[i].second;
      global_id_table->add(b[i].first,v[i],copy(*b[i].second));
    }
    std::cout << std::endl;
  }
  @< Catch block for errors thrown during a global identifier definition @>
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
    std::cerr << " not created." << std::endl;
  }
  reset_evaluator(); main_input_buffer->close_includes();
}
catch (std::logic_error& err)
{ std::cerr << "Unexpected error: " << err.what()
            << ", evaluation aborted.\n";
  reset_evaluator(); main_input_buffer->close_includes();
}
catch (std::exception& err)
{ std::cerr << err.what() << ", evaluation aborted.\n";
  reset_evaluator(); main_input_buffer->close_includes();
}

@ The following function is called when an identifier is declared with type
but undefined value.

@< Function definitions @>=
extern "C"
void global_declare_identifier(Hash_table::id_type id, ptr t)
{@; value undef=NULL;
   global_id_table->add(id,shared_value(undef),copy(*static_cast<type_p>(t)));
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

@ The function |show_ids| prints a table of all known identifiers and their
types.

@< Function definitions @>=
extern "C"
void show_ids()
@+{@; *output_stream << *global_id_table; }

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
