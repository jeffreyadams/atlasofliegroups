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

@* Outline.
%
This file describes a fundamental central part of the interpreter for the
(new) command language of the Atlas of Lie Groups and Representation software.
This part describes the fundamental type declarations used throughout the
interpreter, as well as some functions that operate on types, for instance to
see if one con be converted into another. This compilation unit is
called \.{types}, which refers both to user types as to meta-types of the
interpreter used to represent user types and other fundamental notions.

@( types.h @>=

#ifndef TYPES_H
#define TYPES_H

@< Includes needed in \.{types.h} @>@;
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

@h "types.h"
@h <cstdlib>

@c

namespace atlas { namespace interpreter {
namespace {
@< Local type definitions @>@;
}@; // |namespace|
@< Global variable definitions @>@;
@< Function definitions @>@;
}@; }@;


@ The parser produces a parse tree, in the form of a value of type |expr|
defined in the unit \.{parsetree}. Our task is to take such an expression,
analyse it and then evaluate it to obtain a value. At various points errors
can occur, which we shall have to handle gracefully: during the analysis
static ``type'' errors may prevent us from undertaking any meaningful action,
and in absence of such errors there still can be dynamic ``runtime'' errors.

We first define some data types used by the evaluator. First of all we shall
consider ``type'' values, represented by the structure |type_expr|, in terms
of which the static analysis is performed. Then we shall define the structure
of values for user data, the kind manipulated by the evaluator at runtime;
these are of classes derived from~|value_base|. Then we shall define values
that represent the user expression after type analysis, and which serve to
control the runtime actions; these are of classes derived from
|expression_base|. Finally we need some types for errors that we might have
to throw; these are the classes |program_error|, and |type_error| which is
derived from it.

Originally the |expr| value itself (possibly slightly modified during type
analysis) was used for evaluation, but it turns out to be useful to rebuild
the expression tree with a modified structure. It used to be a side benefit
that we could liberate ourselves from the constraint on data types built by
the parser, that they could be understood by a \Cee-program; currently
however \Cpp~language types can be, and are, part of the type |expr|.
Nonetheless the structure built on |expression_base| is different in many
respects from |expr|, and the transformation between these structures involves
certain processes like overloading resolution that more resemble compilation
than interpretation.

Most of these types are recursive in a manner more complicated than simple
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

@* Types to represent types.
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

@< Type definitions @>=
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

@< Declarations of exported functions @>=
type_ptr copy(const type_expr& t);

@~The function |copy| simply creates an auto-pointer from the call of |new|
for |type_expr(t)|, where the latter invokes the copy constructor of
|type_expr|. The latter, has a recursive definition (given in section
@#type expression copy@>) that performs a deep copy, but the details don't
concern us here.

@< Function definitions @>=
type_ptr copy(const type_expr& t) @+
{@; return type_ptr(new type_expr(t)); }

@*1 Type lists.
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

@< Function definitions @>=
type_node::type_node(const type_node& n)
 : t(n.t),next(n.next==NULL ? NULL : new type_node(*n.next))
@+{}

@ The destructor for |type_node| will possibly be called implicitly by the
destruction of a |type_list_ptr|. The destructor for variants containing a
|type_list| should call |delete| for that pointer, which will also call the
destructor below, which will recursively clean up the whole list.

@< Function def... @>=
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

@< Declarations of exported functions @>=
type_list_ptr make_type_singleton(type_ptr t);
type_list_ptr make_type_list(type_ptr t,type_list_ptr l);
size_t length(type_list l);

@ We pass the auto-pointers by value, rather than by (necessarily non-constant)
reference. This implies a small overhead due to multiple ownership transfer,
but allows (due to a clever definition of auto-pointers) passing anonymous
values obtained from a constructor or construction function as argument, which
would be forbidden if using non-constant reference parameters.

@< Function def... @>=

type_list_ptr make_type_singleton(type_ptr t)
{@; return type_list_ptr(new type_node(t,type_list_ptr(NULL))); }
@)
type_list_ptr make_type_list(type_ptr t,type_list_ptr l)
{@; return type_list_ptr(new type_node(t,l)); }

@ Here is one more function that is convenient to have around.

@< Function definitions @>=
size_t length(type_list l)
@+{@; size_t len=0;
  for (; l!=NULL; l=l->next) ++len;
  return len;
}

@*1 Primitive types.
%
Types can be formed in a variety of ways (primitive types, tuple types,
function types, etc.), and correspondingly will involve a |union| of various
alternatives. So before we can define |type_expr|, we need to enumerate
its variants and the possibilities for primitive types. Some of them will be
introduced later.

@< Type definitions @>=
enum type_tag @+
{ undetermined_type, primitive_type, row_type, tuple_type, function_type };

enum primitive_tag
{ integral_type, rational_type, string_type, boolean_type
  , @< Other primitive type tags @>@;@; @|
  @+nr_of_primitive_types };

@ For printing types (and later for parsing them) we shall need names for the
primitive ones.

@< Declarations of global variables @>=

extern const char* prim_names[];

@~Here is the list of names of the primitive types (some are given later),
terminated by a null pointer. The last name |"void"| is not a primitive name,
corresponds to |nr_of_primitive_types|, and will be treated exceptionally in
|make_prim_type| to make an empty tuple type instead.

@< Global variable definitions @>=
const char* prim_names[]=@/
{"int","rat","string","bool",@< Other primitive type names@>@;@;@+
 "void", NULL };

@*1 Type expressions.
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

@< Definition of |type_expr| @>=
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
  bool can_specialise(const type_expr& pattern) const;
private:
  type_expr& operator=(const type_expr& t);
   // assignment is forbidden
};


@ For that definition to be processed properly, we must pay some attention to
ordering of type definitions, because of the recursions present. The structure
|func_type| will contain instances of |type_expr|, so its definition must
follow, but since we used a pointer to such a structure, it must be declared
as a structure before the definition of |type_expr| is seen. Similarly the
structure |type_node| used in type lists must be declared before the
definition of |type_expr|, but we already did that above; on the other hand
the definition of that structure, which is also already presented, should be
seen by the compiler \emph{after} that of |type_expr|, so we tend to that when
inserting the module reference for that definition.

@< Type definitions @>=

struct func_type; // must be predeclared for |type_expr|
@< Definition of |type_expr| @>

@< Definition of |struct type_node@;| @>@; // this must \emph{follow}

@ The constructor for function types cannot be defined inside the structure
definition, since |func_type| is not (and cannot be) a complete type there.
The constructor for |func_type| used will be defined below.

@< Template and inline function definitions @>=
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

@< Type definitions @>=
struct func_type
{ type_expr arg_type, result_type;
@)
  func_type(type_ptr a, type_ptr r)
   : arg_type(), result_type()
  {@; arg_type.set_from(a); result_type.set_from(r); }
  func_type(const func_type& f)
   : arg_type(f.arg_type),result_type(f.result_type) @+{}
};
typedef func_type* func_type_p;
typedef std::auto_ptr<func_type> func_type_ptr;

@ As we remarked above, the copy constructor for a |type_expr|
recursively copies the descendant types; this is necessary because these
descendants are owned by the object. As usual the recursion is implicit in the
use of a copy constructor within |new|. In all cases there is only one object
directly pointed to, so there is no need for intermediate auto-pointers: no
exception can be thrown between the moment that the call to |new| returns (if
it does), and the completion of our constructor, which incorporates the
returned pointer.
@:type expression copy@>

@< Function definitions @>=
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

@< Function definitions @>=
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
@< Function definitions @>=
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
modifications to our type may still have been made. This is no problem in most
situations, since failure to specialise $t_1$ to $t_2$ will usually be
followed by an attempt to coerce $t_2$ to $t_1$, or by throwing of an error;
here any specialisation that brings $t_1$ closer to $t_2$ cannot be harmful
(it probably makes no difference at all). And we provide an accessor method
|can_specialise| that could be tested before calling |specialise| in cases
where having commit-or-roll-back is important.

@< Function definitions @>=
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


@ Here is the definition of the accessor |can_specialise|, which is quite
similar.

@< Function definitions @>=
bool type_expr::can_specialise(const type_expr& pattern) const
{ if (pattern.kind==undetermined_type)
    return true; // specialisation to \.* trivially succeeds.
  if (kind==undetermined_type)
    return true;
  if (pattern.kind!=kind) return false; // impossible to refine
  switch(kind)
  { case primitive_type: return prim==pattern.prim;
    case row_type:
      return component_type->can_specialise(*pattern.component_type);
    case function_type:
      return func->arg_type.can_specialise(pattern.func->arg_type) @|
         and func->result_type.can_specialise(pattern.func->result_type);
    case tuple_type:
    { type_list l0=tuple, l1=pattern.tuple;
      while (l0!=NULL and l1!=NULL and l0->t.can_specialise(l1->t))
        {@; l0=l0->next; l1=l1->next; }
      return l0==NULL and l1==NULL;
    }
    default: return true; // to keep the compiler happy, cannot be reached
  }
}

@ For printing types, we shall pass |type_expr| values to the
operator~`|<<|' by constant reference, which seems more decent than doing
so by pointer (which would override the definition that simply prints the
hexadecimal address of a pointer); we shall not define instances of~`|<<|' for
other pointer types either. Since we often hold types in |type_ptr| values,
this does mean the we must dereference explicitly in printing.

@< Declarations of exported functions @>=
std::ostream& operator<<(std::ostream& out, const type_expr& t);

@~The cases for printing the types are fairly straightforward. Only
function types are somewhat more involved, since we  want to suppress
additional parentheses around argument and result types in case these are
tuple types; defining a separate operator for a |type_list| facilitates our
task a bit.

@< Function definitions @>=

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

@< Declarations of exported functions @>=
bool operator== (const type_expr& x,const type_expr& y);
inline bool operator!= (const type_expr& x,const type_expr& y)
{@; return !(x==y); }

@~This code is quite similar to the |specialise| method; in fact one could
often use that method instead of the equality operator, but here we want both
operands to be |const|.

@< Function definitions @>=
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
@< Declarations of exported functions @>=
type_ptr make_undetermined_type();
type_ptr make_prim_type(primitive_tag p);
type_ptr make_row_type(type_ptr c);
type_ptr make_tuple_type (type_list_ptr l);
type_ptr make_function_type (type_ptr a, type_ptr r);

@ The reasons for, and the way of, using them are the
same as explained above for |make_type_list|. Note that |make_prim_type|
contrary to what its name suggests returns and empty tuple type for the type
name |"void"|.

@< Function definitions @>=
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

@*1 Specifying types by strings.
%

In practice we shall rarely call functions like |make_prim_type| and
|make_row_type| directly to make explicit types, since this is rather
laborious. Instead, such explicit types will be constructed by the function
|make_type| that parses a string, and correspondingly calls the appropriate
type constructing functions.

@< Declarations of exported functions @>=
type_ptr make_type(const char* s);

@ The task of converting a properly formatted string into a type is one of
parsing a simple kind of expressions. We are not going to write incorrect
strings (we hope) so we don't care if the error handling is crude. The
simplest way of parsing ``by hand'' is recursive descent, so that is what we
shall use. By passing a character pointer by reference, we allow the recursive
calls to advance the index within the string read.

@< Function definitions @>=
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

@< Function definitions @>=
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

@*1 Predefined type expressions.
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

@< Declarations of global variables @>=
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

@< Global variable definitions @>=

const type_expr unknown_type; // uses default constructor
 type_expr void_type(type_list_ptr(NULL));
 type_expr int_type(integral_type);
 type_expr bool_type(boolean_type);
const type_expr row_of_type(copy(unknown_type));
const type_expr gen_func_type(copy(unknown_type),copy(unknown_type));

@ There are more such statically allocated type expressions, which are used in
the evaluator. They are less fundamental, as they are not actually used in any
of the core language constructs, but useful for instance for specifying
various coercions. The reason they should be initialised in this compilation
unit is explained in the next section.

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
extern const type_expr Lie_type_type; // \.{LieType}
extern const type_expr rd_type; // \.{RootDatum}
extern const type_expr ic_type; // \.{InnerClass}
extern const type_expr rf_type; // \.{RealForm}
extern const type_expr split_type; // \.{Split}
extern const type_expr param_type; // \.{RealForm}
extern const type_expr param_pol_type; // \.{RealForm}

@ Since some of these types are built from earlier defined ones, it is vital
that they are initialised in order, and since we cannot control the relative
order of static initialisation between compilation units, we must initialise
them here (this used no not be the case, and led to a subtle bug whose
appearance depended on the precise compiler version used!). The construction
of type constants follows the same pattern as before, calling |copy| in the
case of composite types. In the final case we choose the simplest solution of
calling |make_type| and copying the resulting nested structure into the static
variable before destroying the function result.

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
const type_expr Lie_type_type(complex_lie_type_type);
const type_expr rd_type(root_datum_type);
const type_expr ic_type(inner_class_type);
const type_expr rf_type(real_form_type);
const type_expr split_type(split_integer_type);
const type_expr param_type(module_parameter_type);
const type_expr param_pol_type(virtual_module_type);

@ We shall also need a tuple pattern with any number of unknown components;
such a type is built by calling |unknown_tuple|.

@< Declarations of exported functions @>=
type_ptr unknown_tuple(size_t n);

@~This pattern is built up in a simple loop.

@< Function definitions @>=
type_ptr unknown_tuple(size_t n)
{ type_list_ptr tl(NULL);
  while (n-->0) tl=make_type_list(copy(unknown_type),tl);
  return make_tuple_type(tl);
}

@* Dynamically typed values.
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

As mentioned values are always handled via pointers. We define a raw pointer
type |value|, an auto-pointer |owned_value| (which cannot be stored in STL
containers), and a shared smart pointer |shared_value| (which by contrast can
be stored in STL containers).

@< Type definitions @>=
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

@< Declarations of exported functions @>=
std::ostream& operator<< (std::ostream& out, const value_base& v);

@~The operator~`|<<|' handles calling the virtual |print| method of the actual
value as usual.

@< Function definitions @>=
std::ostream& operator<< (std::ostream& out, const value_base& v)
{@; v.print(out); return out; }

@ Often we know what variant a |value| object takes, based on the type
analysis. We can convert to that type using a |dynamic_cast|, but at such
moments we wish to throw a |logic_error| in case our type prediction was
wrong. To avoid having such casts and |throw| statements all over the place,
we define a template function to do the casting and throwing. It is defined at
the level of ordinary pointers, and it is not intended for use where the
caller assumes ownership of the result; the original pointer is assumed to
retain ownership as long as the result of this call survives, and in
particular that result should not be converted to a smart pointer, lest double
deletion would ensue.

@< Template and inline function definitions @>=
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

@~Since the actual values accessed will be of types derived from |value_base|,
we must pass through a level of indirection, so we have a vector of pointers.
We define these pointers to be |shared_value| pointers, so that the row takes
(shared) ownership of its components without needing a explicit destructor.
This has the additional advantage over explicit ownership management that the
copy constructor, needed for the |clone| method, can safely just
copy-construct the vector of pointers: a possible exception thrown during the
copy is guaranteed to clean up any pointers present in the vector. Note also
that default-constructed shared pointers are set to null pointers, so the
constructor below, which already reserves space for |n| shared pointers, has
set them to exception-safe values while waiting for the slots to be filled.

Of course ownership of pointers to |row_value| objects also needs to be
managed, which could be either by an auto-pointer |row_ptr| (if the pointer is
known to be unshared) or by a shared pointer |shared_row|.

@< Type definitions @>=
struct row_value : public value_base
{ std::vector<shared_value> val;
@)
  explicit row_value(size_t n) : val(n) @+{} // start with |n| null pointers
  void print(std::ostream& out) const;
  size_t length() const @+{@; return val.size(); }
  row_value* clone() const @+{@; return new row_value(*this); }
    // copy the outer level vector
  static const char* name() @+{@; return "row value"; }
protected:
  row_value(const row_value& v) : val(v.val) @+{}
    // copy still shares the individual entries
};
@)
typedef std::auto_ptr<row_value> row_ptr;
typedef std::tr1::shared_ptr<row_value> shared_row;

@ So here is the first occasion where we shall use virtual functions. For the
moment the output routine performs an immediate recursion; later we shall try
to make this more elegant by computing the width needed to output component
values, and adapt the formatting to that.

@< Function definitions @>=
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

@< Type definitions @>=
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
@< Function definitions @>=
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


@< Declarations of exported functions @>=
void push_tuple_components();
void wrap_tuple(size_t n);

@~These functions use the same convention for stack order: the last tuple
component is on top of the stack. In |push_tuple_components| we start by
popping a value and checking it to be a tuple, which is done by
|get<tuple_value>| that will be defined later. The |shared_value| type takes
care of ownership; there is (at least) double shared ownership of the
components as the stack expands, but this is normal, and afterwards this
sharing disappears with |tuple|.

@< Function definitions @>=
void push_tuple_components()
{ shared_tuple tuple=get<tuple_value>();
  for (size_t i=0; i<tuple->length(); ++i)
    push_value(tuple->val[i]); // push component
}

@ We need no auto-pointer in |wrap_tuple|, as shrinking the stack will not
throw any exceptions.

@< Function definitions @>=
void wrap_tuple(size_t n)
{ shared_tuple result(new tuple_value(n));
  while (n-->0) // standard idiom; not |(--n>=0)|, since |n| is unsigned!
    result->val[n]=pop_value();
  push_value(result);
}

@*1 Representation of an evaluation context.
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

@< Type definitions @>=
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

@< Function def... @>=
shared_value& context::elem(size_t i, size_t j)
{
  context* p=this;
  for (; i-->0; p=p->next.get())
    assert(p->next.use_count()>0 and p->next.get()!=NULL);
  assert(j<p->frame.size());
  return p->frame[j];
}

@* Values representing type-checked expressions.
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

@< Type definitions @>=
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

@< Declarations of exported functions @>=
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

@< Declarations of global variables @>=
extern std::vector<shared_value> execution_stack;

@~We define the stack as a static variable of this compilation unit; it is
initially empty. All usable built-in functions will be provided with a small
wrapper function that takes it values from the stack and places its results
there again. Parameters are placed on the stack in order, and should therefore
be popped from the stack in reverse order.

@< Global variable definitions @>=
std::vector<shared_value> execution_stack;

@ We shall define some inline functions to facilitate manipulating the stack.
The function |push_value| does what its name suggests. For exception safety it
takes either an auto-pointer or a shared pointer as argument; the former is
converted into the latter, in which case the |use_count| will become~$1$. The
former form cannot be made to take a reference argument: if the reference were
constant it would be impossible to transfer ownership to the shared pointer,
and if it were non-constant it would be impossible to bind the argument to
expressions other than variables (i.e., to rvalues), such as |owned_value(p)|
below. The shared pointer version does not have this restriction, as a copy
construction and assignment from constant references are possible here.

For convenience we make these template functions that accept a smart pointer
to any type derived from |value_base| (since a conversion of such pointers
from derived to base is not possible without a cast in a function argument
position). For even more convenience we also provide a variant taking an
ordinary pointer, so that expressions using |new| can be written without cast
in the argument of |push_value|. Since |push_value| has only one argument,
such use of does not compromise exception safety: nothing can throw between
the return of |new| and the conversion of its result into an auto-pointer.

@< Template and inline function definitions @>=
template<typename D> // |D| is a type derived from |value_base|
  inline void push_value(std::auto_ptr<D> v)
  {@; execution_stack.push_back(std::tr1::shared_ptr<D>(v)); }

template<typename D> // |D| is a type derived from |value_base|
  inline void push_value(const std::tr1::shared_ptr<D>& v)
  @+{@; execution_stack.push_back(v); }

inline void push_value(value_base* p) @+{@; push_value(owned_value(p)); }

@ There is a counterpart |pop_value| to |push_value|. Most often the result
must be dynamically cast to the type they are known to have because we passed
the type checker; should the cast fail we shall throw a |std::logic_error|.
The template function |get| with explicitly provided type serves for this
purpose; it is very much like the template function |force|, but returns a
shared pointer (because the value on the stack might be shared).

@< Template and inline function definitions @>=

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

@ Sometimes we may need to expand a value into tuple components separately
pushed onto the stack, but only if the |level l@;| so indicates and the value
is indeed of tuple type; the function |push_expanded| will help doing this.

@< Declarations of exported functions @>=
void push_expanded(expression_base::level l, const shared_value& v);

@~Type information is not retained in compiled expression values, so
|push_expanded| cannot know which type had been found for |v|, but it can use
a dynamic cast do determine whether it actually is a tuple value or not.

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


@< Template and inline function def... @>=
inline void uniquify(shared_value& v)
{@; if (not v.unique()) v=shared_value(v->clone()); }
@)
template <typename D> // |D| is a type derived from |value_base|
 std::tr1::shared_ptr<D> get_own() throw(std::logic_error)
{@; uniquify(execution_stack.back());
    return get<D>();
}

@* Implicit conversion of values between types.
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

The function |coerce| requires two fully determined types |from_type| and
|to_type|, and its final argument~|e| is a reference to the previously
converted expression. If a conversion of value of |from_type| to |to_type| is
available, then |coerce| will modify |e| by insertion of a conversion around
it; the return value of |coerce| indicates whether an applicable conversion
was found. The function |conform_types| first tries to specialise the type
|required| to the one |found|, and if this fails tries to coerce |found| to
|required|, in the latter case enveloping the translated expression |d| in the
applied conversion function; if both fail an error mentioning the
expression~|e| is thrown. The function |row_coercion| specialises if possible
|component_type| in such a way that the corresponding row type can be coerced
to |final_type|, and returns a pointer to the |conversion_record| for the
coercion in question. The function |coercion| serves for filling the coercion
table.

@< Declarations of exported functions @>=

bool coerce(const type_expr& from_type, const type_expr& to_type,
            expression_ptr& e);
expression conform_types
  (const type_expr& found, type_expr& required, expression_ptr d, expr e);
const conversion_record* row_coercion(const type_expr& final_type,
                                            type_expr& component_type);
void coercion(const type_expr& from,
              const type_expr& to,
              const char* s, conversion_info::conv_f f);

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

@< Type definitions @>=
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
after evaluating |exp|. The |level| argument is not passed to the |convert|
function, which will always replace a one or more values on the stack by a
single value. There is a language design decision implicit in this
implementation: there are no implicit conversions that return a tuple type. In
fact we tried some such conversions, for instance from a rational number (and
later from a split integer) to a pair of integers, and this was unsatisfactory
even when correctly implemented. One reason is that it disturbs operator and
function overloading: one can no longer define operators for the
converted-from type if the operator or function already exists for the tuple
type converted to, for instance one could not define unary minus for rational
numbers because binary minus for integers was already defined. Another reason
is that using decomposition of tuples in a let-expression to disassemble the
converted-from type will not work without a cast-to-a-tuple, since in this
context the mere desire to have some unspecified tuple does not suffice to
activate the implicit conversion. For these reasons it is preferable to always
make the conversion to a tuple explicit.

Although automatic conversions are only inserted when the type analysis
requires a non-empty result type, it is still possible that at run time this
method is called with |l==no_value|, so we do cater for that here. The case is
fairly rare, so we don't mind the inefficiency of performing the conversion
and then discarding the result; this will allow a failing conversion to be
signalled as an error in such cases.

@< Function def...@>=
void conversion::evaluate(level l) const
{@; exp->eval();
  (*type.convert)();
  if (l==no_value)
    execution_stack.pop_back();
}
@)
void conversion::print(std::ostream& out) const
@+{@; out << type.name << ':' << *exp; }

@*1 Coercion of types.
%
An important aspect of automatic conversions is that they can be applied in
situations where the result type and only part of the source type is known:
for instance one can be inserted when a list display occurs in a context
requiring a vector, because no row type equals the (primitive) vector type.
While this use requires particular consideration according to the syntactic
form of the expression (list display), there is also a simpler form of
automatic conversion that can be applied to a large variety of expressions
(identifiers, function calls, \dots) whenever they are found to have a
different type from what the context requires. The function |coerce| will
try to insert an automatic conversion in such situations, if this can resolve
the type conflict. We present this mechanism first, since the table it employs
can then be re-used to handle the more subtle cases of automatic conversions.

@ The implementation of |coerce| will be determined by a simple table lookup.
The records in this table contain a |conversion_info| structure (in fact they
are derived from it) and in addition indications of the types converted from
and to. The table entries store these types by pointer, so that table entries
are assignable, and the table can be allowed to grow (which is a technical
necessity in order to allow other compilation units to contribute coercions).

@< Type definitions @>=
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

@< Global variable definitions @>=

std::vector<conversion_record> coerce_table;

@ We will iterate over |coerce_table| several times; the following |typedef|
makes this easier.

@< Local type def... @>=

typedef std::vector<conversion_record>::const_iterator coerce_iter;

@ The function |coercion| simplifies filling the coercion table; it is
externally callable. Its action is simply extending |coerce_table| with a new
|conversion_record|.
@< Function def... @>=
void coercion(const type_expr& from,
              const type_expr& to,
              const char* s, conversion_info::conv_f f)
{@; coerce_table.push_back(conversion_record(from,to,s,f)); }

@ There is one coercion that is not stored in the lookup table, since it can
operate on any input type: the voiding coercion. It is necessary for instance
to allow a conditional expression that is intended for its side effects only,
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

@< Type definitions @>=
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

@< Function definitions @>=
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
Note that in the |conversion| constructor, the first argument |*it| is used
only for its |conversion_info| base type; the |from| and |to| fields are not
accessible to the newly built expression. As last resort we build a |voiding|
if |to_type==void_type|.

@< Function definitions @>=
bool coerce(const type_expr& from_type, const type_expr& to_type,
	    expression_ptr& e)
{ for (coerce_iter
       it=coerce_table.begin(); it!=coerce_table.end(); ++it)
    if (from_type==*it->from and to_type==*it->to)
    @/{@; e.reset(new conversion(*it,e));
      return true;
    }
  if (to_type==void_type)
  {@; e.reset(new voiding(e));
     return true;
  }
  return false;
}

@ Often we first try to specialise a required type to the available type of a
subexpression, or else (if the first fails) coerce the available type to the
one required. The function |conform_types| will facilitate this. The argument
|d| is a possibly already partially converted expression, which should be
further wrapped in a conversion call if appropriate, while |e| is the original
expression that should be mentioned in an error message if both attempts fail.
A call to |conform_types| will be invariably followed (upon success) by
returning the expression~|d| (possibly modified) from the calling function; we
can save some work by returning the required value already from
|conform_types|. Given that we do, we might as well take ownership of~|d| (by
having it passed by value) so that the caller knows it should afterwards use
our return value, rather than~|d|.

@~The copied auto-pointer |d| provides the modifiable reference that |coerce|
needs. If both attempts to conform the types fail, we must take a copy of both
type expressions, since the originals are not owned by us, and will probably
be destructed before the error is caught.

@< Function def... @>=
expression conform_types
(const type_expr& found, type_expr& required, expression_ptr d, expr e)
{ if (not required.specialise(found) and not coerce(found,required,d))
    throw type_error(e,copy(found),copy(required));
  return d.release();
}


@ List displays and loops produce a row of values of arbitrary (but identical)
type; when they occur in a context requiring a non-row type, we may be able to
find a coercion that reconciles the requirements. The following function finds
whether this is possible, and if so sets |components_type| to the type
required for the components; if not it will return~|NULL|. By using this
mechanism the components themselves obtain a context that may generate further
conversions to obtain this type. The implementation is simply to look in
|coerce_table| for conversions from some row type.

@< Function def... @>=
const conversion_record* row_coercion(const type_expr& final_type,
                                            type_expr& component_type)
{ for (coerce_iter it=coerce_table.begin(); it!=coerce_table.end(); ++it)
    if (final_type==*it->to and it->from->kind==row_type)
    @/{@; component_type.specialise(*it->from->component_type); return &*it; }
  return NULL;
}

@*1 Proximity of types.
%
The above implicit conversions of types pose a limitation to the possibilities
of overloading operator symbols and function identifiers. If a symbol should
be overloaded for too closely related operand types, situations could occur in
which given operand expressions can be converted to either of the operand
types. This would either produce unpredictable behaviour, or necessitate a
complicated set of rules to determine which of the definitions of the symbol
is to be used (overloading resolution in \Cpp\ is a good example of such
complications). There are two ways to avoid the occurrence of complications by
restricting the rules of the language: either forbid type conversions in
arguments of overloaded symbols, or forbid simultaneous definitions of such
symbols for too closely related types. (A mixture of both is also conceivable,
allowing only certain conversions and forbidding overloading between types
related by them; the language Algol~68 is a good example of an approach along
these lines). Forbidding all automatic type conversions in case of overloading
would defeat to a large extent the purpose of overloading, namely as a
convenience to the user; therefore we choose the latter solution of forbidding
overloading in certain cases. The predicate |is_close| will be used to
characterise pairs of argument types that are mutually exclusive for
overloading purposes.

@< Declarations of exported functions @>=
unsigned int is_close (const type_expr& x, const type_expr& y);

@ We do allow simultaneous overloading between closely related types in some
cases, namely if they can be ordered so that if $t_1$ precedes $t_2$ then some
expression of type~$t_1$ can be converted to type~$t_2$ but no expression of
type~$t_2$ can be converted to type~$t_1$; in such cases reasonable behaviour
can be obtained by trying a match for~$t_1$ before trying one for~$t_2$, and
this allows for instance arithmetic operators to be defined for
type \.{(int,int)} as well as for type \.{(rat,rat}). Therefore |is_close|
returns a value composed of 3~bits: one indicating whether the types are close
at all, and two others for indicating the existence of conversions in one
direction or the other. Thus |is_close| returns an integer rather than a
boolean value, with the following interpretations: |0x0| means the types are
unrelated, |0x4| means the types are mutually exclusive but neither can be
converted to the other (as for instance \.{(int,rat)} and \.{(rat,int)}),
|0x5| means the types are close and (only) the first can be converted to the
second (example, \.{int} and \.{rat}), |0x6| is the opposite relation, and
|0x7| means both types can be converted to each other (like \.{vec}
and \.{[int]}, or any case of equal types).

For types $t_1$ and~$t_2$ which do not admit a relative priority, we want to
disallow simultaneous overloading with arguments types $t_1$ and~$t_2$ if any
expression given as argument could be converted to either of them. Deciding
the existence of such an expression would require study of all available
language constructs, but the situation is somewhat simplified by the fact
that, for efficiency reasons, overloading resolution is not done using the
argument expression, but only its type. In fact matching will be done using
calls to the very function |is_close| we are discussing here, testing the bit
for conversion towards the required argument type; this provides us with an
opportunity to adjust rules for possible type conversions of arguments at the
same time as defining the exclusion rules.

We must forbid \.{void} altogether as operand type of overloaded functions,
since anything can be converted to that type; this is not a great limitation.
Apart from that, the function |coerce| provides us with the basic information
about possible conversions; however we do not limit ourselves to these, because
conversions might be possible inside row or tuple displays, and we want the
consider for instance \.{[[rat]]} as matching a required type \.{[ratvec]}.
Empty row displays, or more precisely arguments of type~\.{[*]}, pose a
difficulty: either we forbid using such arguments in overloaded calls, or we
must accept that for any pair of row types, no matter how different their
components, there exists arguments that can be converted to either type, so
that they are mutually exclusive for overloading.

Since empty rows as arguments are probably quite common and forcing a specific
type on them relatively tedious, we opt for the latter solution. Thus when
comparing two row types, |is_close| will always set the bit for closeness, but
the other two bits will be set to indicate the convertibility of the component
types. Also, in accordance with the choice to allow~\.{[*]} as operand type,
the (component) type \.{*} will be considered to convert to any type. On the
other hand we disallow an empty row where a primitive type with conversion
from some row type (like \.{vec} or \.{mat}) is required, so that these types
can coexist with an unrelated row type for overloading purposes.

@ So here is the (recursive) definition of the relation |is_close|. Equal
types are always close, while undetermined types behave as convertible to any
type. A primitive type is in the relation |is_close| to another type only if
it is identical or if there is a direct conversion between the types, as
decided by~|coerce|. Two row types are always in the relation |is_close|, and
two tuple types are so if they have the same number of component types, and if
each pair of corresponding component types is (recursively) in the relation
|is_close|.

The above applies to the most significant of the three bits used in the result
of the function |is_close|. The two other bits indicate whether, by applying
conversions from~|coerce| to corresponding (nested) component types of row and
tuple types, one of the types can be converted to the other. If the |is_close|
relation is false all bits will be zero, but in the contrary case any of the 4
possible combinations of the remaining bits could be set.

The expression |dummy| may be prepended to by |coerce|, but is then abandoned
(and cleaned up).

@< Function definitions @>=
unsigned int is_close (const type_expr& x, const type_expr& y)
{ expression_ptr dummy(NULL);
  if (x==y)
    return 0x7;
  if (x.kind==undetermined_type)
    return 0x5; // |x| matches when |y| is required
  if (y.kind==undetermined_type)
    return 0x6; // |y| matches when |x| is required
  if (x.kind==primitive_type or y.kind==primitive_type)
  { unsigned int flags=0x0;
    if (coerce(x,y,dummy)) flags |= 0x1;
    if (coerce(y,x,dummy)) flags |= 0x2;
    return flags==0 ? flags : flags|0x4;
  }
  if (x.kind!=y.kind)
    return 0x0;
  if (x.kind==row_type)
    return 0x4 | is_close(*x.component_type,*y.component_type); // always close
  if (x.kind!=tuple_type)
    return 0x0; // non-aggregate types are only close if equal
  type_list l0=x.tuple, l1=y.tuple; // recurse for tuples, as in |operator==|
  unsigned int flags=0x7;
  while (l0!=NULL and l1!=NULL and (flags&=is_close(l0->t,l1->t))!=0)
    @/{@; l0=l0->next; l1=l1->next; }
  return l0==NULL and l1==NULL ? flags : 0x0;
    // lists must end simultaneously for success
}

@* Error values.
Before we describe evaluation of expressions we must realise that evaluation
can cause runtime errors. The evaluator may throw exceptions due to
inconsistency of our (rather than the user's) program, which are classified as
|std::logic_error|. It may also throw exceptions due to errors not caught by
the type checker in the user input, such as size mismatch in matrix
operations; in such cases it will throw |std::runtime_error|. These standard
exception types will be used without any type derivation.

@< Includes needed in \.{types.h} @>=
#include <stdexcept>
#include "parsetree.h" // type |expr| will be used in error classes

@ For errors detected before execution starts, we first derive a general
exception class |program_error| from |std::exception|; it represents any kind
of error of the user input determined by static analysis (for instance use of
undefined variables).

Although it does nothing explicitly (the string will be destructed anyway), we
must explicitly define a destructor for |program_error|: the automatically
generated one would lack the |throw()| specifier.

@< Type definitions @>=
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

@< Type definitions @>=
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

@< Type definitions @>=
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

@< Function definitions @>=
type_error::type_error(const type_error& e)
 : expr_error(e)
{@; type_ptr p=copy(*e.required);
  actual=copy(*e.actual).release(); required=p.release();
}

@* Enumeration of primitive types.
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
dual_real_form_type, Cartan_class_type, KGB_element_type,
module_parameter_type, split_integer_type, virtual_module_type, @[@]


@~The following list must match that of the previous module, for proper
functioning of I/O.

@< Other primitive type names @>=
"vec", "mat", "ratvec",
"LieType","RootDatum", "InnerClass", "RealForm", "DualRealForm",
"CartanClass", "KGBElt", "Param", "Split", "ParamPol", @[@]

@* Index.

% Local IspellDict: british
