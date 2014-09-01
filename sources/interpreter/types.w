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
#include <memory> // for |std::unique_ptr|, |std::shared_ptr|
#include "sl_list.h" // for singly linked lists

@* Types to represent types.
%
We make a fundamental choice to check types before attempting to execute any
expression entered by the user; thus type errors can be signalled earlier and
expressed more understandably than if one would wait until an accident is
about to happen at runtime. This means that types are to be represented
independently of runtime values. Of course the interpreter executes its same
internal code regardless of the types of the values it manipulates; those
values are then accessed via generic pointers. The runtime values also contain
an indication of their actual type, which will allow us to do some dynamic
testing to trap possible errors in the logic of our interpreter, but the user
should never notice this.

Types are represented by small tree structures, with nodes of types
|type_expr| that will be detailed later. For now we provide no user-defined
encapsulation (as \Cpp~classes do), so a type expression exposes the
accessible operations: a type printed as \.{(int->[(bool,mat)])} for instance
will specify a function mapping integers to lists of pairs of a Boolean value
and a matrix, which implies the kind of operation that can be performed with
such a value. Note however that classes defined in the Atlas software library
itself will be presented to the user as primitive types, so we do provide
abstraction there. One kind of type is a function type, which specifies the
types of the individual function arguments and of the values it returns.

Different types will not share any sub-trees among each other, as is the case
for STL container types; this simplifies memory management. This choice means
some recursive copying of tree structures is sometimes required, but (although
runtime cost of type handling is negligible with respect to other factors in
the interpreter) we avoid doing so more than absolutely necessary.

Usually types are built from the bottom up (from leaves to the root), although
during type checking the reverse also occurs. This means that nodes of the
structure are held in local variables before being moved to dynamically
managed storage; it is important to be able to do this while transferring
ownership of data linked to, which requires resetting the original node so
that it can be cleaned up without destroying dependent data. This is called
move semantics, and was originally realised using auto-pointers for individual
pointers, and by a special method on the level of complete nodes. With the
advent of \Cpp11, special syntactic support was added, so that one can
distinguish calls that should be realised with shallow copies and ownership
transfer from the occasionally needed deep copy. At the same time
auto-pointers were replaced by unique-pointers, with essentially the same
functionality, but better adapted to the syntactic facilities. The structure
also uses ordinary pointers of type~|type_p| whose ownership is managed by the
containing structure.

@< Type definitions @>=
struct type_expr;
typedef type_expr* type_p;
typedef const type_expr* const_type_p;
typedef std::unique_ptr<type_expr> type_ptr;

@*2 Type lists.
%
We also need type lists as building block for types (for instance for the
arguments of a function). These used to be defined in a similar manner to
types themselves, but since an STL-compatible singly-linked list container
class templates |simple_list| and |sl_list| was added to the Atlas utilities
library, these were used to replace the implementation of type lists. The
former is built into the structure itself, the latter which is mode flexible
but requires more space for the class instance itself is occasionally used for
temporary variables. This provide a good test for the usability of the new
container type.

@< Type definitions @>=
typedef containers::simple_list<type_expr> type_list;
typedef atlas::containers::sl_node<type_expr>* raw_type_list;
typedef containers::sl_list<type_expr> dressed_type_list;

@ Since types and type lists own their trees, their copy constructors must
make a deep copy. The class |type_expr| will provide no copy constructor but
instead a more visible |copy| method to do a deep copy. On some occasions
however one rather needs a copy at the pointer level: one has access to a
pointer to some |type_expr| (possibly by simply taking its address), but one
needs an owning |type_ptr| instead. The function |acquire| will achieve this.

@< Declarations of exported functions @>=
type_ptr acquire(const type_expr* t);

@~The function |acquire| simply creates a unique-pointer from the call of
|new| whose constructor invokes the |copy| method of |type_expr| pointed to.
The latter, has a recursive definition that performs a deep copy, but the
details (given in section @#type expression copy@>) don't concern us here.

@< Function definitions @>=
type_ptr acquire(const type_expr* t) @+
{@; return type_ptr(new type_expr(t->copy())); }


@ Type lists are usually built by starting with a default-constructed
|type_list| (which is equivalent to initialising it with |empty_tuple()|
defined below, and which function can help avoid a Most Vexing Parse
situation), and repeatedly calling |prefix| to add nodes in front. This
function is efficient, due its use of move-semantics, both for the node, and
for the list itself, although the latter is a modifiable lvalue reference
because it will continue to hold the extended list.

@< Declarations of exported functions @>=
type_list empty_tuple();
type_list& prefix(type_expr&& t, type_list& dst);
dressed_type_list& prefix(type_expr&& t, dressed_type_list& dst);

@ Move semantics, introduced in \Cpp11, solves any difficulty with ownership
management that this code previously had to deal with.

@< Function def... @>=
type_list empty_tuple() @+{@; return type_list(); }
@)
type_list& prefix(type_expr&& t, type_list& dst)
{@; dst.push_front(std::move(t));
   return dst;
}
dressed_type_list& prefix(type_expr&& t, dressed_type_list& dst)
{@; dst.push_front(std::move(t));
  return dst;
}

@*2 Primitive types.
%
We have a simple but flexible type model. There is a finite number of
``primitive'' types, many of which are abstractions of complicated classes
defined in the Atlas library, such as root data or reductive groups.
Furthermore one has types that are ``row of'' some other type, tuples
(Cartesian products) of some given sequence of types, and function types. We
also allow for an undetermined type, which can serve as a wild-card to specify
type patterns. Before we can define |type_expr|, we need to enumerate its
variants and the possibilities for primitive types. Here are enumerations of
tags for the basic kinds of types, and for the sub-cases of |primitive_type|
(some of them will be introduced later).

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
|mk_prim_type| to make an empty tuple type instead.

@< Global variable definitions @>=
const char* prim_names[]=@/
{"int","rat","string","bool",@< Other primitive type names@>@;@;@+
 "void", nullptr };

@*1 Type expressions.
%
Type expressions are defined by a tagged union. Code accessing the union will
in general test the tag and access the corresponding variant directly, so we
give public access to the data fields (it is not obvious how coherent variant
selection could be enforced by private data and public methods without greatly
complicating usage; the \Cpp~model of abstraction does not seem particularly
suited for tagged unions). By using an anonymous union, the field selectors
like |prim| of the variants in the union can be used directly on the level of
the |type_expr|, thus avoiding an additional level of selection to
access them.

The field name |tupple| was chose because when working with a program
involving some |unique_pre| instance, the \.{gdb} debugger cannot access any
variable or field named |tuple|; the was reported as bug~17098
on \.{sourceware.org/bugzilla}. Awaiting resolution of the bug, circumvent it
by using a voluntary misspelling of the name.

There is one restriction on types that is not visible in the definition below,
namely that the list of types referred to by the |tupple| field cannot have
length~$1$ (but length~$0$ is allowed). This is because anything that would
suggest a $1$-tuple (for instance a parenthesised expression) is identified
with its unique component.

@< Definition of |type_expr| @>=
struct type_expr
{ type_tag kind;
  union
  { primitive_tag prim; // when |kind==primitive|
    type_p component_type; // when |kind==row_type|
    type_list tupple; // when |kind==tuple_type|
    func_type* func; // when |kind==function_type|
  };
@)
  @< Methods of the |type_expr| structure @>@;
 };

@ Every variant of the union gets its own constructor that directly constructs
the |type_expr| structure under this variant of the structure, and in addition
there is a default constructor that sets |kind==undetermined_type| and
constructs no variant at all (in \Cpp\ \emph{at most} one variant of a union
is active). For cases where it will be necessary (or simply more convenient)
to first default-construct a |type_expr| and then change its variant, one will
have to take care to create the variant first by placement-|new|, unless the
variant is a POD type in which case one can simply assign to the field
instead. The destructor for |type_expr| must similarly take care to explicitly
call the destructor for the active field.

A move constructor for |type_expr| is provided, but no copy constructor;
instead of the latter we provide a |copy| method explicitly creating a deep
copy value, to which one may then apply move-construction. We similarly
provide a move assignment operator, which implicitly makes the copy assignment
operator deleted (rather then implicitly provided), and a |swap| method (as
proof-of-concept, it is not currently used). Two special cases of (partial)
move semantics are provided for however, for replacing an undefined (or
partially undefined) type by a more specific type (note that in these cases
nothing disappears, so no clean-up is necessary). The first case, the
|set_from| method, makes a shallow copy of its argument into the object for
which is was called, which is required to be undetermined initially. This
method that was present long before \Cpp11 allowed proper move semantics is in
fact used to implement the move constructor and move assignment operator. The
second case, the |specialise| method, is used during type analysis, to see if
our type matches a given pattern, or in case it was (partially) undefined
whether it can be made to match the pattern by if necessary replacing some
undetermined descendants by more specific ones. The call returns a value
indicating whether this was possible, and if so makes the necessary
specialisations to our type. That is done by copying, so the caller does not
require or lose ownership of the pattern for this method.

The constructors for the row and tuple types receive pointers that will be
directly inserted into our |struct| and become owned by it; the ownership
management is simplest when such pointers are ``passed by check'', i.e., as
rvalue references to smart pointers.

@< Methods of the |type_expr| structure @>=

type_expr() noexcept : kind(undetermined_type) @+{}
explicit type_expr(primitive_tag p) noexcept
  : kind(primitive_type), prim(p) @+{}
explicit type_expr(type_ptr&& c) noexcept
  : kind(row_type), component_type(c.release()) @+{}
explicit type_expr(type_list&& l) noexcept; // tuple types
explicit type_expr(dressed_type_list&& l) noexcept; // tuple types
inline type_expr(type_expr&& arg, type_expr&& result) noexcept;
 // for function types
@)
type_expr(type_expr&& t) noexcept; // move constructor
void clear() noexcept;
  // resets to undefined state, cleaning up owned pointers
~type_expr() noexcept @+{@; clear(); }
  // that is all explcitly needed for destruction
type_expr& operator=(type_expr&& t) noexcept; // do move assignment only
void swap(type_expr& t) noexcept;
@)
type_expr copy() const; // in lieu of deep copy constructor
void set_from(type_expr&& p) noexcept; // shallow copy
bool specialise(const type_expr& pattern);
  // try to match pattern, possibly modifying |*this|
bool can_specialise(const type_expr& pattern) const;
  // tell whether |specialse| would succeed

@ For that definition to be processed properly, we must pay some attention to
ordering of type definitions, because of the recursions present. The structure
|func_type| will contain instances of |type_expr|, so its definition must
follow, but since we used a pointer to such a structure, it must be declared
as a structure before the definition of |type_expr| is seen. For type lists
this kind of difficulty is taken care of because |simple_list<type_expr>|
could be used in a |typedef| before |type_expr| was a complete type.

@< Type definitions @>=

struct func_type; // must be predeclared for |type_expr|
@< Definition of |type_expr| @>

@ The constructor for the |type_list| variant with |tupple| field is a move
constructor, and is easily implemented since |type_expr| has a move
constructor.
@< Function definitions @>=
type_expr::type_expr(type_list&& l) noexcept
  : kind(tuple_type), tupple(std::move(l))
  @+{}

@ Instead of a regular copy constructor, which would have to make a deep copy
(because these descendants are owned by the object), |type_expr| provides a
method |copy| that recursively copies the descendant types; this way the
inadvertent making of deep copies is avoided. If necessary the |type_expr|
move constructor can be applied to the temporary for the result from |copy|,
but this can probably always be avoided by return value optimisation (and in
any case move construction costs little). Using a named function here results
in the recursion being visible in the argument the calls of |new|. Since this
is not a constructor, the pointers returned from |new| are immediately owned
when stored in a component of |type_expr| (and even if this had been a
constructor, no exception can happen between the storing the pointer and
termination).

@:type expression copy@>

@< Function definitions @>=
type_expr type_expr::copy() const
{ type_expr result;
  switch (result.kind=kind)
  { case undetermined_type: break;
    case primitive_type: result.prim=prim; break;
    case row_type:
      result.component_type=new type_expr(component_type->copy());
    break;
    case tuple_type:
      @< Placement-construct a deep copy of |tupple| into |result.tupple| @>
    break;
    case function_type: result.func=new func_type(func->copy()); break;
  }
  return result;
}

@ First off, as the module name says we must make sure a valid object is
constructed into the field |result.tupple|, because we default constructed
|result| with no variant active. We cannot use a copy constructor of
|type_list| to do this, because since |type_expr| has a deleted copy
constructor, |type_list| has no copy constructor either (more precisely, the
class template declares such a constructor, but an attempt to use it will not
compile because of the missing |type_expr| has no copy constructor). We must
instead construct an empty shell first and then fill it, applying the |copy|
method for all members of the list |tupple|. This is achieved by first
default constructing |result.tupple|, and then adding fields, keeping an
(output) iterator |oit| pointing at the end of the list, where |insert| can be
used to add node (note that the class template |simple_list| defines no method
|end|, as this would be expensive is used like |end| methods usually are).

@h <algorithm>

@< Placement-construct a deep copy of |tupple| into |result.tupple| @>=
{ new (&result.tupple) @[type_list@]; // construct empty object
  auto oit = result.tupple.begin();
  for (auto it = tupple.begin(); not tupple.at_end(it); ++it,++oit)
    result.tupple.insert(oit,it->copy());
}

@ The method |clear|, doing the work for the destructor, must similarly clean
up afterwards, with the recursion again being implicit. We set |kind =
undetermined_type| to indicate that no variant is active; this condition is
tested in the |set_from| method (and also in by the destructor, though if
|clear| is \emph{called by} the destructor, that test is already behind us).

There is no reason to use a smart pointers as variant of a union, since one
must still explicitly propagate destruction to the active variant, as happens
here for the |tupple| field. If it were done here, then all that would change
is that the code below would have to call their destructors rather than call
|delete| directly.

@< Function definitions @>=
void type_expr::clear() noexcept
{ switch (kind)
  { case undetermined_type: case primitive_type: break;
    case row_type: delete component_type; break;
    case tuple_type: tupple.~type_list(); break;
    case function_type: delete func; break;
  }
  kind = undetermined_type;
}

@ The method |set_from| makes a shallow copy of the structure; it implements
 move semantics. Correspondingly it takes as argument another |type_expr| by
 modifiable rvalue reference. The contents of its top-level structure will be
 moved to the current |type_expr| using move assignment where applicable, and
 then set to a empty |undetermined_type| value, which effectively detaches any
 possible descendants from it. This operation requires that |type_expr|
 previously had |kind==undetermined_type|. We test this condition using an
 |assert| statement (rather than throwing |std::logic_error|) to honour the
 |noexcept| specification.

@h <stdexcept>
@< Function definitions @>=
void type_expr::set_from(type_expr&& p) noexcept
{ assert (kind==undetermined_type);
   // logic should ensure this; we promised not to |throw|
  switch(kind=p.kind) // copy top node
  { case undetermined_type: break;
    case primitive_type: prim=p.prim; break;
    case row_type: component_type=p.component_type; break;
    case function_type: func=p.func; break;
    case tuple_type: new (&tupple) type_list(std::move(p.tupple)); break;
  }
  p.kind=undetermined_type;
  // detach descendants, so |p.clear()| will destroy top-level only
}
@)
type_expr::type_expr(type_expr&& x) noexcept // move constructor
: kind(x.kind)
{ switch(kind) // move top node
  { case undetermined_type: break;
    case primitive_type: prim=x.prim; break;
    case row_type: component_type=x.component_type; break;
    case function_type: func=x.func; break;
    case tuple_type: new(&tupple) type_list(std::move(x.tupple)); break;
  }
  x.kind=undetermined_type;
  // detach descendants, so destructor of |x| will do nothing
}

@ For move assignment we reuse the |set_from| method, and for |swap| we do a
move construction and two |set_from| calls, unless the |kind| fields match, in
which case we can call |std::swap| directly on the matching variant fields.

@< Function definitions @>=
type_expr& type_expr::operator=(type_expr&& x) noexcept // move assignment
{ if (this!=&x)
  { clear(); // detach anything previously linked
    set_from(std::move(x)); // move top level structure
  }
  return *this;
}
@)
void type_expr::swap(type_expr& other) noexcept
{
  if (kind==other.kind) // an easy case
    switch(kind)
    { case undetermined_type: break; // no need to swap |nothing| fields
      case primitive_type: std::swap(prim,other.prim); break;
      case row_type: std::swap(component_type,other.component_type); break;
      case function_type: std::swap(func,other.func); break;
      case tuple_type: std::swap(tupple,other.tupple); break;
    }
  else
  {
    type_expr t(std::move(other));
    other.set_from(std::move(*this));
    this->set_from(std::move(t));
  }
}

@ The |specialise| method is mostly used to either set a completely
undetermined type to a given pattern, or to test if it already matches it;
however, we do not exclude the possibility that a partly determined type is
modified by specialisation of one of its descendants to match the given
pattern. Matching a pattern means being at least as specific, and specialising
means replacing by something more specific, so we are (upon success) replacing
the type for which the method is called by the most general unifier (i.e.,
common specialisation) of its previous value and |pattern|. Therefore the name
of this method is somewhat misleading, in that the specialisation is not
necessarily to |pattern|, but to something matching |pattern|.

In the case of an |undetermined_type|, |specialise| uses the |set_from| method
to make |*this| a copy of |pattern|. In the other cases we only continue if
the top levels of both type declarers match, in which case we try to
recursively specialise all descendants. We do not guarantee
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
  if (kind==undetermined_type) // specialising \.* also always succeeds,
    {@; set_from(pattern.copy()); return true; }
     // by setting |*this| to |pattern|
  if (pattern.kind!=kind) return false; // impossible to refine
  switch(kind)
  { case primitive_type: return prim==pattern.prim;
    case row_type: return component_type->specialise(*pattern.component_type);
    case function_type:
      return func->arg_type.specialise(pattern.func->arg_type) @|
         and func->result_type.specialise(pattern.func->result_type);
    case tuple_type:
     @< Try to specialise types in |tupple| to those in |pattern.tupple|,
        and |return| whether this succeeded @>
    default: return true; // to keep the compiler happy, cannot be reached
  }
}

@ For tuples, specialisation is done component by component. We check
beforehand that the lengths of the lists match, without which there is no hope
of finding a proper specialisation.

@< Try to specialise types in |tupple| to those in |pattern.tupple|... @>=
{ auto l0=tupple.begin();
  auto l1=pattern.tupple.begin();
  while (not tupple.at_end(l0) and not pattern.tupple.at_end(l1)
         and l0->specialise(*l1))
    @/{@; ++l0; ++l1; }
  return tupple.at_end(l0) and pattern.tupple.at_end(l1);
  // whether both lists terminated
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
    { auto l0=tupple.begin(), l1=pattern.tupple.begin();
      while (not tupple.at_end(l0) and not pattern.tupple.at_end(l1)
             and l0->can_specialise(*l1))
      @/{@; ++l0; ++l1; }
      return tupple.at_end(l0) and pattern.tupple.at_end(l1);
      // whether both lists terminated
    }
    default: return true; // to keep the compiler happy, cannot be reached
  }
}

@ The constructor for function types cannot be defined inside the structure
definition, since |func_type| is not (and cannot be) a complete type there.
The constructor for |func_type| used will be defined below.

@< Template and inline function definitions @>=
type_expr::type_expr(type_expr&& arg, type_expr&& result) noexcept
@/ : kind(function_type)
   , func(new func_type(std::move(arg),std::move(result)))
   @+{}

@ The variant for function types needs both an argument type and a result
type. We cannot define the appropriate structure directly within the |union|
of |type_expr| where it is needed, since the scope of that definition
would then be too limited to perform the appropriate |new| in the constructor.
Therefore we must declare and name the structure type before using it in
|type_expr|. Afterwards we define the structure, and while we are doing
that, we also define a constructor to fill the structure (it was used above)
and a move constructor.

The |func_type| structure has two sub-objects of type |type_expr|,
rather than pointers to them. A move constructor only makes a
shallow copy of the topmost nodes in the same way as was done in the move
constructor for~|type_node|. Like for |type_expr| we do not provide an
ordinary copy constructor, but a value-returning |copy| accessor method.

@s result_type normal

@< Type definitions @>=
struct func_type
{ type_expr arg_type, result_type;
@)
  func_type(type_expr&& a, type_expr&& r)
@/ : arg_type(std::move(a)), result_type(std::move(r)) @+{}
  func_type(func_type&& f) = @[default@]; // move constructor
  func_type& operator=(func_type&& f) = @[default@]; // move assignment
  func_type copy() const // in lieu of a copy contructor
  {@; return func_type(arg_type.copy(),result_type.copy()); }
};
typedef func_type* func_type_p;
typedef std::unique_ptr<func_type> func_type_ptr;

@ For printing types, we shall pass |type_expr| values to the
operator~`|<<|' by constant reference, which seems more decent than doing
so by pointer (which would override the definition that simply prints the
hexadecimal address of a pointer); we shall not define instances of~`|<<|' for
other pointer types either. Since we often hold types in |type_ptr| values,
this does mean the we must then dereference explicitly in printing.

@< Declarations of exported functions @>=
std::ostream& operator<<(std::ostream& out, const type_expr& t);
std::ostream& operator<<(std::ostream& out, const func_type& f);

@~The cases for printing the types are fairly straightforward. Only
function types are somewhat more involved, since we  want to suppress
additional parentheses around argument and result types in case these are
tuple types; defining a separate operator for a |type_list| facilitates our
task a bit.

@< Function definitions @>=

std::ostream& operator<<(std::ostream& out, const type_list& l)
{ for (auto it=l.begin(); not l.at_end(it); ++it)
    out << *it << ( l.at_end(std::next(it)) ? "" : "," );
  return out;
}
std::ostream& operator<<(std::ostream& out, const func_type& f)
{
  out << '(';
  if (f.arg_type.kind==tuple_type)
     out << f.arg_type.tupple; // naked tuple
  else out << f.arg_type; // other component type
  out << "->";
  if (f.result_type.kind==tuple_type)
     out << f.result_type.tupple; // naked tuple
  else out << f.result_type; // other component type
  out << ')';
  return out;
}
@)
std::ostream& operator<<(std::ostream& out, const type_expr& t)
{ switch(t.kind)
  { case undetermined_type: out << '*'; break;
    case primitive_type: out << prim_names[t.prim]; break;
    case row_type: out << '[' << *t.component_type << ']'; break;
    case tuple_type:
      out << '(' << t.tupple << ')' ;
    break;
    case function_type: out << *t.func;
    break;
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
// all cases are listed below, but compilers want a |default| to |return|
  @\case undetermined_type: return true;
    case primitive_type: return x.prim==y.prim;
    case row_type: return *x.component_type==*y.component_type;
    case tuple_type:
    { auto it0=x.tupple.begin(), it1=y.tupple.begin();
      while (not x.tupple.at_end(it0) and not y.tupple.at_end(it1)
             and *it0==*it1)
      @/{@; ++it0; ++it1; }
      return x.tupple.at_end(it0) and y.tupple.at_end(it1);
    }
    case function_type:
      return x.func->arg_type==y.func->arg_type
	 and x.func->result_type==y.func->result_type;
  }
}

@ Instead of using the constructors directly, we often use the constructing
functions below. They all and return |type_ptr| values owning the constructed
expression. They also take such smart pointers, or |type_list| values, as
argument whenever the underlying pointer is to be directly inserted into the
structure built. However for |mk_function_type| this is not the case, so
instead of insisting that the caller hold a unique-pointer to the argument
types it suffices to hold a |type_expr| whose contents can be moved into the
type to be constructed. If |t| is a |type_ptr| held by a client, it can pass
|std::move(*t)| as argument to |mk_function_type|, which will move the
contents of the node |t| points to into the function type, and after return
the destructor of |t| will eventually delete the now empty node.

@< Declarations of exported functions @>=
type_ptr mk_prim_type(primitive_tag p);
type_ptr mk_row_type(type_ptr&& c);
type_ptr mk_tuple_type (type_list&& l);
type_ptr mk_function_type(type_expr&& a, type_expr&& r);
@)
type_p make_prim_type(int p);
type_p make_row_type(type_p c);
type_p make_tuple_type(raw_type_list l);
type_p make_function_type(type_p a,type_p r);
@)
raw_type_list make_type_singleton(type_p raw);
raw_type_list make_type_list(type_p t,raw_type_list l);

@ The functions like |mk_prim| below simply wrap the call to the constructor
into one of the operator |new|, and then capture of the resulting pointer into
a |type_ptr| result. The functions like |make_prim| wraps these functions into
a parser interface, which converts it arguments to unique pointers, and
converts such pointers to raw pointers for the return value; these functions
should be used exclusively in the parser. Using this setup it is ensured that
all pointers are considered owning their target, implicitly so while being
manipulated by the parser (it guarantees that every pointer placed on the
parsing stack will be argument of an interface function exactly once, possibly
some |destroy| function in case it pops symbols during error recovery.

Note that we make a special provision that |mk_prim_type| will return an
empty tuple type when called with the type name for |"void"|, although this is
contrary to what the name of the function suggests.

@< Function definitions @>=
type_ptr mk_prim_type(primitive_tag p)
{ return p<nr_of_primitive_types ?
    type_ptr(new type_expr(p)) :
    type_ptr(mk_tuple_type(empty_tuple()));
}

type_ptr mk_row_type(type_ptr&& c)
{@; return type_ptr (new type_expr(std::move(c))); }

type_ptr mk_tuple_type (type_list&& l)
{@; return type_ptr(new type_expr(std::move(l))); }

type_ptr mk_function_type (type_expr&& a, type_expr&& r)
{@; return type_ptr(new type_expr(std::move(a),std::move(r)));
}
@)
type_p make_prim_type(int p)
{@; return mk_prim_type(static_cast<primitive_tag>(p)).release(); }

type_p make_row_type(type_p c)
{@; return mk_row_type(type_ptr(c)).release(); }

type_p make_tuple_type(raw_type_list l)
{@; return mk_tuple_type(type_list(l)).release(); }

type_p make_function_type(type_p a,type_p r)
{@; return
    mk_function_type(std::move(*type_ptr(a)),std::move(*type_ptr(r))).release();
}
@)

raw_type_list make_type_singleton(type_p t)
{@; type_list result;
  result.push_front(std::move(*type_ptr(t)));
  return result.release();
}

raw_type_list make_type_list(type_p t,raw_type_list l)
{ type_list tmp(l); // since |prefix| needs second argument an lvalue reference
  return prefix(std::move(*type_ptr(t)),tmp).release();
}

@*1 Specifying types by strings.
%
In practice we shall rarely call functions like |mk_prim_type| and
|mk_row_type| directly to make explicit types, since this is rather laborious.
Instead, such explicit types will be constructed by the function |mk_type|
that parses a (\Cee~type) string, and correspondingly calls the appropriate
type constructing functions.

@< Declarations of exported functions @>=
type_ptr mk_type(const char* s);
type_expr mk_type_expr(const char* s);
  // ``exported'' for our global variable initialisation

@ The task of converting a properly formatted string into a type is one of
parsing a simple kind of expressions. The strings used here come from string
denotations in the source code (mostly in calls installing built-in \.{realex}
function) rather than from user input, and we are not going to write incorrect
strings (we hope). Therefore we don't care if the error handling is crude
here. The simplest way of parsing ``by hand'' is recursive descent, so that is
what we shall use. By passing a character pointer by reference, we allow the
recursive calls to advance the index within the string read.

The function |scan_type| does the real parsing, |mk_type_expr| calls it,
providing a local modifiable pointer to bind to its reference parameter (which
is important because |scan_type| cannot directly accept a \Cee-string constant
as argument) while also doing error reporting, and |mk_type| is just a
wrapper around |mk_type_expr| that converts the result from a |type_expr| to
a smart pointer to (a freshly allocated instance of) such. Currently the
function |mk_type_expr| is called only during the start-up phase
of \.{realex}, and if an error encountered (of type |std::logic_error|, since
it indicates an error in the \.{realex} program itself), printing of the error
message will be followed by termination of the program.

@< Function definitions @>=
dressed_type_list scan_type_list(const char*& s);
@)
type_expr scan_type(const char*& s)
{ if (*s=='*') return ++s,@[type_expr()@]; // undetermined type
  else if (*s=='[')
    @< Scan and |return| a row type, or |throw| a |logic_error| @>
  else if (*s=='(')
    @< Scan and |return| a tuple or function type,
       or |throw| a |logic_error| @>
  else @< Scan and |return| a primitive type, or |throw| a |logic_error| @>
}
@)
type_expr mk_type_expr(const char* s)
{ const char* orig=s;
  try
  {@; return type_expr(scan_type(s)); }
  catch (std::logic_error e)
  { std::cerr << e.what() << "; original string: '" << orig @|
              << "' text remaining: '" << s << "'\n";
    throw@[@];
  // make the error hard to ignore; if thrown probably aborts the program
  }
}
@)
type_ptr mk_type(const char* s)
{@; return type_ptr(new type_expr(mk_type_expr(s))); }
  // wrap up |type_expr| in |type_ptr|


@ The following code demonstrates how simple recursive descent parsing can be.
The only subtle point is that did not advance the pointer |s| when testing for
|'['| above, so we must start with doing that.

@< Scan and |return| a row type, or |throw| a |logic_error| @>=
{ type_expr comp=scan_type(++s);
  if (*s++!=']') throw std::logic_error("Missing ']' in type");
    return type_expr(type_ptr(new type_expr(std::move(comp))));
}

@ Now we do tuple and function types. Here again we start by advancing the
pointer. Otherwise the descent is still straightforward, thanks to
|scan_type_list|. After scanning a first list we decide whether this will be a
tuple type (if a right parenthesis follows) and record it in a boolean
variable~|is_tuple| to be able to share the code for constructing the tuple
type.

The only complication is that single parenthesised types, and single argument
or return types should not be converted into tuple types with one component,
but just into the constituent type. This is done by assigning to the local
|type_ptr| variables |a| and |r| either the extracted singleton type or the
tuple type constructed from the non-singleton list. When extracting a
singleton type, the |type_expr| reference |l0->t| or |l1->t| used is one to
part of a |type_node| structure, so there is no other option here than calling
|copy| to turn this reference into a |type_ptr| owning a fresh copy of that
singleton type.

Note that this is one of the few places where we really use the
ownership-tracking semantics of unique-pointers, in the sense that their
destruction behaviour at a certain point is runtime-dependent: the list nodes
pointed to by |l0| and |l1| will only be deleted in the case of a singleton
list, where |type_list_ptr| was not captured and set to null in the
|type_expr| constructor for a tuple type; in that case the node will be
deleted, but its contents has been moved out by then.

@< Scan and |return| a tuple or function type, or |throw| a |logic_error| @>=
{ dressed_type_list l0=scan_type_list(++s), l1;
  bool is_tuple=*s==')';
  if (*s=='-' and *++s=='>') l1=scan_type_list(++s);
  if (*s!=')') throw std::logic_error("Missing ')' in type");
   @+ else ++s;
  type_expr a =
    l0.size()==1
    ? std::move(l0.front())
    : type_expr(l0.undress()); // construct tuple type
  if (is_tuple)
    return a;
  type_expr r =
    l1.size()==1
    ? std::move(l1.front())
    : type_expr(l1.undress()); // construct tuple type
  return type_expr(std::move(a),std::move(r)); // construct function type
}

@ A comma-separated list of types is handled by a straightforward loop.
We must not forget that the list could be empty, which happens only of the
very first character we see is |')'| or |'-'|. Since we advance our pointer in
the test for the presence of a comma extending the list, we must back it up
when the test fails, i.e., after loop exit.

@< Function definitions @>=
@)
dressed_type_list scan_type_list(const char*& s)
{ dressed_type_list result; // start with empty tuple
  if (*s==')' or *s=='-')
    return result;
  do result.push_back(scan_type(s));
  while (*s++==',');
  --s; // back up to character that was not a comma
  return result;
}

@ For primitive types we use the same strings as for printing them. We test as
many characters as the type name has, and the fact that no alphanumeric
character follows, so that a longer type name will not match a prefix of it.

In this module we use the fact lather the order in the list |prim_names|
matches that in the enumeration type |primitive_tag|, by casting the integer
index into the former list to an element of that enumeration.

@h <cstring> // |strlen|, |strncmp|
@h <cctype> // |isalpha|
@< Scan and |return| a primitive type, or |throw| a |logic_error| @>=
{ for (size_t i=0; i<nr_of_primitive_types; ++i)
  { const char* name=prim_names[i];
    size_t l= std::strlen(name);
    if (std::strncmp(s,name,l)==0 and not isalpha(s[l]))
    { s+=l;
      primitive_tag tag = static_cast<primitive_tag>(i);
      if (tag == nr_of_primitive_types) // then name scanned was |"void"|
        return type_expr(empty_tuple()); // which is equivalent to |"()"|
      return type_expr(tag); // otherwise this is a true primitive type
    }
  }
  throw std::logic_error("Type unrecognised");
}

@*1 Predefined type expressions.
%
We shall often need to refer to certain types for comparison or for providing
a required type context. Instead of generating them on the fly each time using
|mk_type|, we define constant values that can be used everywhere.

@< Declarations of global variables @>=
extern const type_expr unknown_type; // \.{*}
extern const type_expr void_type; // \.{()}
extern const type_expr int_type; // \.{int}
extern const type_expr bool_type; // \.{bool}
extern const type_expr row_of_type; // \.{[*]}
extern const type_expr gen_func_type; // \.{(*->*)}

@ In some cases we need temporary copies of |void_type|, |int_type| and
|bool_type| to be used in the position of a modifiable lvalue argument. In
order to provide these temporary copies as arguments without having to bind
them to named variables, we define a function template that will produce a
modifiable lvalue from the modifiable rvalue, such as the result of calling
the |expr::copy| method.

@< Template and inline... @>=
template<typename T> T& as_lvalue(T&& rvalue) @+{@; return rvalue; }

@ The definition of the variables uses the constructors we have seen above, or
calls to |mk_type_expr|, rather than functions like |mk_prim_type| and
|mk_row_type|, so that no dynamic allocation is required for the top level
structure.

@< Global variable definitions @>=

@: first types section @>

const type_expr unknown_type; // uses default constructor
const type_expr void_type(empty_tuple());
const type_expr int_type(integral_type);
const type_expr bool_type(boolean_type);
const type_expr row_of_type(mk_type_expr("[*]"));
const type_expr gen_func_type(mk_type_expr("(*->*)"));

@ There are more such statically allocated type expressions, which are used in
the evaluator. They are less fundamental, as they are not actually used in any
of the core language constructs, but useful for instance for specifying
various coercions. The definition of these constants therefore might be moved
to another compilation unit, but they are defined here so that in case their
initialisation should use other such constants, the order of initialisation
will be controlled (this is not the case between initialisations in different
compilation unit, which can lead to the so-called static initialisation
fiasco; indeed at some point we had for these constants a subtle bug whose
appearance depended on the precise compiler version used).

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

@ The definitions below have all become self-contained, due to the use of
|mk_type_expr| for non-primitive types. Indeed, since we cannot have sharing
between type (sub-)expressions, the economy of using the |copy| method for
previously constructed type constants would be truly marginal. So in their
current form, some of these definitions could now (again) be moved to other
compilation units, where they might even be just local constants.

@: second types section @>

@< Global variable definitions @>=
const type_expr rat_type(rational_type);
const type_expr str_type(string_type);
const type_expr vec_type(vector_type);
const type_expr ratvec_type(rational_vector_type);
const type_expr mat_type(matrix_type);
const type_expr row_of_int_type(mk_type_expr("[int]"));
const type_expr row_of_rat_type(mk_type_expr("[rat]"));
const type_expr row_of_vec_type(mk_type_expr("[vec]"));
const type_expr row_row_of_int_type(mk_type_expr("[[int]]"));
const type_expr pair_type(mk_type_expr("(*,*)"));
const type_expr int_int_type(mk_type_expr("(int,int)"));
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
type_expr unknown_tuple(size_t n);

@~This pattern is built up in a simple loop.

@< Function definitions @>=
type_expr unknown_tuple(size_t n)
{ type_list tl;
  while (n-->0)
    prefix(type_expr(),tl);
  return type_expr(std::move(tl));
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

@~We start with a base class for values. For it to be an abstract class, there
must be at least one pure virtual function in the class; the destructor having
to be virtual anyway, we make it pure virtual (this does not mean it is
unimplemented, in fact it must be implemented as it will always be called,
after the destructor for a derived class; it just means derived
classes \emph{must} override the default). The printing function does not have
a useful default, so we make it pure virtual as well, without providing an
implementation (in the base class). This |print| method will demonstrate the
ease of using dynamic typing via inheritance; it will not do any dynamic
casting, but other operations on values will. Apart from |print| we define
another virtual method, |clone|, which allows making a copy of a runtime value
of any type derived from |value_base|.

The method |name| is useful in reporting logic errors from function templates,
notably the failure of a value to be of the predicted type. Since the template
function may know the type (via a template argument) but need not have any
object of the type at hand, we define |name| as a |static| rather than
|virtual| method. We disable assignment of |value_base| objects, since they
should always be handled by reference; the base class is abstract anyway, but
this ensures us that for no derived class an implicitly defined assignment
operator is accidentally invoked. Copy constructors will in fact be defined
for most derived types, as they are needed to implement the |clone| method;
these will be |private| or |protected| as well, so as to forbid accidental use
elsewhere, but they do copy-construct their |value_base| base object (which is
constructor is therefore made |protected|). Those derived classes might have
value-constructed (i.e., use the no-arguments constructor) the base object
instead (and indeed this used to be the case), but copy-construct is what the
default copy-constructor for the derived class does, so not deleting the
copy-constructor here allows those defaults to be used. Cloning only happens
in specific contexts (namely the function templates |uniquify| and |get_own|
defined below) where the actual type will in fact be known, so some derived
types may choose to never use this and not implement |clone| at all, which is
why it is not defined pure virtual. In fact this state of affairs suggests one
could do duplication with a non-virtual (but systematically named) method
instead.

Values are always handled via pointers. The raw pointer type is |value|, and a
shared smart pointer-to-constant is |shared_value|. The const-ness of the
latter reflects a copy-on-write policy: we rarely need to modify values
in-place, but when we do, we ensure our shared pointer is actually unique, and
then |const_cast| it to |own_value| for modification (this will be hidden in a
function template defined later).

@< Type definitions @>=
struct value_base
{ value_base() @+ {};
  virtual ~value_base() = 0;
  virtual void print(std::ostream& out) const =0;
  virtual value_base* clone() const @+{@; assert(false); return nullptr; }
  static const char* name(); // just a model; this instance remains undefined
protected:
  value_base(const value_base& x) = @[default@];
public:
  value_base& operator=(const value_base& x) = @[delete@];
};
inline value_base::~value_base() @+{} // necessary but empty implementation
@)
typedef value_base* value;
typedef std::shared_ptr<const value_base> shared_value;
typedef std::shared_ptr<value_base> own_value;

@ We can already make sure that the operator~`|<<|' will do the right thing
for any of our values.

@< Declarations of exported functions @>=
std::ostream& operator<< (std::ostream& out, const value_base& v);

@~The operator~`|<<|' handles calling the virtual |print| method of the actual
value as usual.

@< Function definitions @>=
std::ostream& operator<< (std::ostream& out, const value_base& v)
{@; v.print(out); return out; }

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
copy is guaranteed to clean up any pointers present in the vector (resetting
their reference counts to their original values). Note also that
default-constructed shared pointers are set to null pointers, so the
constructor below, which already reserves space for |n| shared pointers, has
set them to exception-safe values while waiting for the slots to be filled.

Of course ownership of pointers to |row_value| objects also needs to be
managed. The type |own_row| will be used after constructing the row while
filling in the contents, or after ensuring unique ownership of the row in
order to perform destructive operations (so although a
|shared_ptr<value_base>|, it is known to actually be unique); at all other
times the pointer converted to |shared_row| (and possibly down-cast to
|shared_value|) will be used.

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
typedef std::shared_ptr<const row_value> shared_row;
typedef std::shared_ptr<row_value> own_row;

@ So here is the first occasion where we shall use virtual functions. For the
moment the output routine performs an immediate recursion; later we shall try
to make this more elegant by computing the width needed to output component
values, and adapt the formatting to that.

@< Function definitions @>=
void row_value::print(std::ostream& out) const
{ out << '[';
  for (auto it=val.begin(); it!=val.end(); ++it)
    out << (it==val.begin() ? "" : ",") << **it;
   out << ']';
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
typedef std::unique_ptr<tuple_value> tuple_ptr;
typedef std::shared_ptr<const tuple_value> shared_tuple;
typedef std::shared_ptr<tuple_value> own_tuple;

@ We just need to redefine the |print| method.
@< Function definitions @>=
void tuple_value::print(std::ostream& out) const
{ out << '(';
  for (auto it=val.begin(); it!=val.end(); ++it)
    out << (it==val.begin() ? "" : ",") << **it;
   out << ')';
}

@ Here are functions that pack and unpack tuples from values on the stack, in
terms of the operations |push_value| and |pop_value| defined in sections
@# Push execution stack @> and @# Pop execution stack @> below.
They will be used by wrapper functions around functions from the Atlas
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
sharing disappears with the destruction of |tuple|.

@< Function definitions @>=
void push_tuple_components()
{ shared_tuple tuple=get<tuple_value>();
  for (size_t i=0; i<tuple->length(); ++i)
    push_value(tuple->val[i]); // push component
}

@ Wrapping a tuple is a simple matter of allocating a |tuple_value| of the
proper size, and then filling it from back to front with shared values popped
from the stack. This was the first place in the Atlas software where a
descending loop using an unsigned loop variable was used; this requires a
post-decrement operation in the test so as to stop \emph{after} handling the
value~$0$. Similar loops can now be found all over the place, and they could
be even more ubiquitous if would have chosen to sacrifice readability for speed
by preferring decreasing loops to increasing ones whenever there is a choice:
the decreasing variant can be slightly more efficient then its increasing
equivalent, because a test against~$0$ can be more efficient than a test
against another value, especially one not known at compile time.

The concrete loop below will however only be invoked for tuples constructed
according to the user program. For uses where an explicit value of $n$ is
known, as happens when called internally from a wrapper function, we provide a
templated version of |wrap_tuple| in section@#templated wrap_tuple section@>,
and that version does not involve a loop at all.

@< Function definitions @>=
void wrap_tuple(size_t n)
{ std::shared_ptr<tuple_value> result = std::make_shared<tuple_value>(n);
  while (n-->0) // standard idiom; not |(--n>=0)|, since |n| is unsigned!
    result->val[n] =pop_value();
  push_value(std::move(result));
}

@*1 Representation of an evaluation context.
%
While evaluating user programs, values will be given to local identifiers such
as arguments of functions being called. The identification of identifiers is
determined during type analysis (static binding); it results for local
identifiers in a method to locate the associated value in the evaluation
context, which is formed of a stack of frames, each holding a vector of
values.

One of the methods for our frame type will produce a back insert iterator, so
we need the following include in order to declare it.

@< Includes needed... @>=
#include <iterator>

@~Frames are actually allocated on the heap, and their lifetimes do not follow
a stack regime unless a very limited use is made of user-defined functions
(never passing such a function as value out of the expression in which it was
defined), so it is better to just say they are linked lists of frames. A
singly linked list suffices, and by using shared pointers as links,
destruction of frames once inaccessible is automatic.

@s back_insert_iterator vector

@< Type definitions @>=
typedef std::shared_ptr<class evaluation_context> shared_context;
class evaluation_context
{ shared_context next;
  std::vector<shared_value> frame;
  evaluation_context@[(const evaluation_context&) = delete@];
  // never copy contexts
public:
  evaluation_context (const shared_context& next)
@/: next(next), frame() @+{}
  void reserve (size_t n) @+{@; frame.reserve(n); }
  shared_value& elem(size_t i,size_t j);
  std::back_insert_iterator<std::vector<shared_value> > back_inserter ()
  {@; return std::back_inserter(frame); }
  const shared_context& tail() const @+{@; return next; }
  std::vector<shared_value>::const_iterator begin() const
    @+{@; return frame.begin(); }
  std::vector<shared_value>::const_iterator end() const
    @+{@; return frame.end(); }
};

@ The method |evaluation_context::elem| descends the stack and then selects a
value from the proper frame.

@< Function def... @>=
shared_value& evaluation_context::elem(size_t i, size_t j)
{
  evaluation_context* p=this;
  while (i-->0 and (p=p->next.get())!=nullptr) {}
  assert(p!=nullptr and j<p->frame.size());
@/return p->frame[j];
}

@* Values representing type-checked expressions.
%
The parser is a \Cpp-program that upon success returns a value of type |expr|
representing the parse tree. While analysing this expression for
type-correctness, it will be convenient to transform it into a value that can
be efficiently evaluated. This value will be a pointer to an object of one of
a number of classes derived, most often directly, from an empty base class (in
the same way as runtime values are pointers to an object of a class derived
from |value_base|), which classes have a virtual method |evaluate| that
performs the operation described by the expression. We shall now define the
base class.

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
  expression_base@[(const expression_base&) = delete@]; // they are never copied
  expression_base@[(expression_base&&) = delete@]; // nor moved
  expression_base& operator=@[(const expression_base&)=delete@]; // nor assigned
  expression_base& operator=@[(expression_base&&)=delete@]; // nor move-assigned
  virtual ~expression_base() @+ {}
@)// other virtual methods
  virtual void evaluate(level l) const =0;
  virtual void print(std::ostream& out) const =0;
@)// non-virtuals that call the virtual |evaluate|
  void void_eval() const @+{@; evaluate(no_value); }
  void eval() const @+{@; evaluate(single_value); }
  void multi_eval() const @+{@; evaluate(multi_value); }
};
@)
typedef expression_base* expression;
typedef std::unique_ptr<expression_base> expression_ptr;
typedef std::shared_ptr<expression_base> shared_expression;

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
whence we use |shared_value| smart pointers in the stack.

This choice will have consequences in many places in the evaluator, since once
a value is referred to by such a smart pointer, its ownership cannot be
transferred to any other regime; when strict ownership should be needed, the
only option would be to make a copy by calling |clone|. However, it turns out
to be convenient to \emph{always} use shared pointers for runtime values, and
to make the distinction concerning whether one knows this pointer to be unique
by having its type be pointer-to-non-const in that case. After construction or
duplication one can start out with such a pointer, use it to store the proper
value, then convert it pointer-to-const (i.e., |shared_value|), for handing on
the stack and passing around in general; if destructive access is needed one
may reconvert to pointer-to-non-const after having checked unique ownership
(or else having duplicated the value pointed to).

@< Declarations of global variables @>=
extern std::vector<shared_value> execution_stack;

@~We define the stack as a static variable of this compilation unit; it is
initially empty. All usable built-in functions will be provided with a small
wrapper function that takes it values from the stack and places its results
there again. Parameters are placed on the stack in order, and should therefore
be popped from the stack in reverse order.

@< Global variable definitions @>=
std::vector<shared_value> execution_stack;

@ Sometimes we may need to expand a value into tuple components separately
pushed onto the stack, but only if the |level l@;| so indicates and the value
is indeed of tuple type; the function |push_expanded| will help doing this.

@< Declarations of exported functions @>=
void push_expanded(expression_base::level l, const shared_value& v);

@~Type information is not retained in compiled expression values, so
|push_expanded| cannot know which type had been found for |v| (moreover,
|push_expanded| is typically called for arguments for which the type is not
determined at the time that \.{realex} is compiled). But it can use a dynamic
cast do determine whether |v| actually is a tuple value or not.

@< Function definitions @>=
void push_expanded(expression_base::level l, const shared_value& v)
{ if(l==expression_base::single_value)
    push_value(v);
  else if (l==expression_base::multi_value)
  { shared_tuple p = std::dynamic_pointer_cast<const tuple_value>(v);
    if (p==nullptr)
      push_value(v);
    else
      for (size_t i=0; i<p->length(); ++i)
        push_value(p->val[i]); // push components
  }
}

@* Some useful function templates.
%
We now define some inline functions to facilitate manipulating the stack. The
function |push_value| does what its name suggests. For exception safety it
takes a shared pointer as argument. The former form used to take an |auto_ptr|
argument by value, which allowed both to transfer ownership from an lvalue
(i.e., a variable) of the same type, and to bind to an rvalue (for instance
the result of a function). With the change to a representation as |unique_ptr|
instance, the lvalue argument case would no longer bind as-is, and an
invocation of |std::move| had to be inserted into the code in more than~$60$
places for this reason (the rvalue case does not need modification). The
argument passing was also changed to modifiable rvalue reference, with the
same syntactic obligations for the caller; this avoids one transfer of
ownership, doing so only when the pointer is converted to a |shared_ptr| in
the code below. Finally it was realised that there is no advantage to first
creating a unique pointer, so we now always create a shared pointer for values
that will be pushed onto the stack; this could be a pointer-to-non-const to a
type derived from |value_base|, which upon passing to |push_value| will be
converted to a |shared_value| by the appropriate constructor of the
|shared_ptr| template (very conveniently, this constructor is not marked as
|explicit|). The shared pointer version of |push_value| can take its argument
as a constant lvalue reference (since it does not need to modify the pointer,
just the reference count), or as rvalue reference.

@: Push execution stack @>

@< Template and inline function definitions @>=
inline void push_value(const shared_value& v)
  @+{@; execution_stack.push_back(v); }

inline void push_value(shared_value&& v)
  @+{@; execution_stack.push_back(std::move(v)); }

@ There is a counterpart |pop_value| to |push_value|. By move-constructing
from the stack top just before it is popped, we avoid incrementing and then
immediately decrementing the |use_count| value.

@: Pop execution stack @>

@< Template and inline function definitions @>=

inline shared_value pop_value()
{@; shared_value arg(std::move(execution_stack.back()));
  execution_stack.pop_back();
  return arg;
}

@ Most often the result of calling |pop_value| must be dynamically cast to the
type it is known to have (by the type check that was passed); should such a
cast fail, this reveals a flaw of our type system, so we throw a
|std::logic_error|. The function template |get| with explicitly provided type
serves for this purpose; it returns a shared pointer (because values on the
stack are shared pointers) of the proper kind.

Sometimes defeating copy-on-write is desired (to allow changes like filling
internal tables that will \emph{benefit} other shareholders), and
|non_const_get| will const-cast the result of |get| to allow that.

@< Template and inline function definitions @>=

template <typename D> // |D| is a type derived from |value_base|
 inline std::shared_ptr<const D> get() throw(std::logic_error)
{ std::shared_ptr<const D> p=std::dynamic_pointer_cast<const D>(pop_value());
  if (p.get()==nullptr)
    throw std::logic_error(std::string("Argument is no ")+D::name());
  return p;
}
@.Argument is no ...@>

@)
template <typename D> // |D| is a type derived from |value_base|
  inline std::shared_ptr<D> non_const_get() throw(std::logic_error)
{@; return std::const_pointer_cast<D>(get<D>()); }

@ Here is a function template similar to |get|, that applies in situations
where the value whose type is known does not reside on the stack. As for |get|
we convert using a |dynamic_cast|, and to throw a |logic_error| in case our
type prediction was wrong. This function is defined at the level of ordinary
pointers, and it is not intended for use where the caller assumes ownership of
the result; the original pointer is assumed to retain ownership as long as the
result of this call survives, and in particular that pointer should probably
not be obtained by calling the |get| method for a smart pointer temporary, nor
should the result of |force| converted to a smart pointer, lest double
deletion would ensue.

We provide two versions, where overloading will choose one or the other
depending on the const-ness of the argument. Since calling |get| for a
|shared_value| pointer returns a pointer to constant, it will often be the
second one that is selected.

@< Template and inline function definitions @>=
template <typename D> // |D| is a type derived from |value_base|
 D* force(value v) throw(std::logic_error)
{ D* p=dynamic_cast<D*>(v);
  if (p==nullptr) throw
    std::logic_error(std::string("forced value is no ")+D::name());
  return p;
}
@)
template <typename D> // |D| is a type derived from |value_base|
const D* force(const value_base* v) throw(std::logic_error)
{ const D* p=dynamic_cast<const D*>(v);
  if (p==nullptr) throw
    std::logic_error(std::string("forced value is no ")+D::name());
  return p;
}

@ In some cases a wrapper function will want to get unique access to an object
(so that we can modify it to produce its result), which requires duplicating
the value if the value on the stack has |use_count()>1|. In fact if we just
computed the value on the stack from a function call, it is virtually
guaranteed to be unshared. Similarly the component assignment operation must
ensure that the name of the aggregate that is being assigned to is made to
hold a unique (non-shared) instance of its value which can then be modified in
place (this was the original motivation for this functionality). The operation
|uniquify| implements this, and calling it makes clear our destructive
intentions. We make it take a modifiable lvalue argument into which a new
shared (but currently unique) pointer is stored in case duplication was
necessary, and return a |value| raw pointer-to-non-const version of the
possibly modified value of that pointer, which can then be used to make the
change to the unique copy (the original pointer cannot, since it is of course
still a pointer-to-const).

For the more common case of arguments on the stack, we provide a variant
function template |get_own| of |get|. It has the same prototype as
|non_const_get|, but like |uniquify| respects copy-on-write by making a copy
first in case there are other shareholders. Since these functions return
pointers that are guaranteed to be unique, one might wonder why no use of
|std::unique_ptr| is made. The answer is this is not possible, since there is
no way to persuade a |shared_ptr| to release its ownership (as in the
|release| method of unique pointers), even if it happens to be (or is known to
be) the unique owner.

@< Template and inline function def... @>=
template <typename D> // |D| is a type derived from |value_base|
  std::shared_ptr<D> get_own() throw(std::logic_error)
{ std::shared_ptr<const D> p=get<D>();
  if (p.unique())
    return std::const_pointer_cast<D>(p);
  return std::shared_ptr<D>(p->clone());
}
@)
inline value uniquify(shared_value& v)
{ if (not v.unique())
     v=shared_value(v->clone());
  return const_cast<value>(v.get());
}

@ The argument~$n$ to |wrap_tuple| most often is a compile time constant, so
we give a templated version that can be completely unrolled by the compiler.

@:templated wrap_tuple section@>

@< Template and inline... @>=
template<unsigned int> void wrap_tuple();
template<unsigned int>
  void do_wrap(std::vector<shared_value>::iterator it);
template<>
   inline void do_wrap<0u>(std::vector<shared_value>::iterator it) @+{}
template<unsigned int n>
   inline void do_wrap(std::vector<shared_value>::iterator it)
   {@; *--it = pop_value(); do_wrap@[<n-1>@](it); }
template<unsigned int n>
   inline void wrap_tuple()
   { std::shared_ptr<tuple_value> result = std::make_shared<tuple_value>(n);
     do_wrap<n>(result->val.end());
     push_value(std::move(result));
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
integers in positions where rational numbers are required.

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
expression_ptr conform_types
  (const type_expr& found, type_expr& required
  , expression_ptr&& d, const expr& e);
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
  expression_ptr exp;
public:
  conversion(const conversion_info& t,expression_ptr e)
   :type(t),exp(std::move(e)) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ The |evaluate| method for conversions dispatches to the |convert| member,
after evaluating |exp|. The |level| argument is not passed to the |convert|
function, which will always replace a value on the stack by a single other
value. There is a language design decision implicit in this implementation:
there are no implicit conversions that return a tuple type. In fact we tried
some such conversions, for instance from a rational number (and later from a
split integer) to a pair of integers, and this was unsatisfactory even when
correctly implemented. One reason is that it disturbs operator and function
overloading: one can no longer define operators for the converted-from type if
the operator or function already exists for the tuple type converted to; for
instance one could not define unary minus for rational numbers because binary
minus for integers was already defined. Another reason is that using
decomposition of tuples in a let-expression to disassemble the converted-from
type will not work without a cast to a specific tuple type, since in this
context the mere desire to have some unspecified tuple does not suffice to
activate the implicit conversion. For these reasons it is preferable to always
have the conversion to a tuple be an explicit function.

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
and to. Since these types occur in a small fixed collection of types, all of
which are statically created in sections @# first types section @> and~%
@# second types section @>, the type expressions are not stored in the table
itself, which instead stores non-owning pointers to them.

@< Type definitions @>=
struct conversion_record : public conversion_info
{ const type_expr* from, * to;
  // non-owned pointers, themselves could be |const| in gcc 4.8
  conversion_record @| (const type_expr& from_type,
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

@ The function |coercion| simplifies filling the coercion table; it is
externally callable. Its action is simply extending |coerce_table| with a new
|conversion_record|.
@< Function def... @>=
void coercion(const type_expr& from,
              const type_expr& to,
              const char* s, conversion_info::conv_f f)
{@; coerce_table.emplace_back(from,to,s,f); }

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
{ expression_ptr exp;
public:
  voiding(expression_ptr e) : exp(e.release()) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ The |evaluate| method should not ignore its |level| argument completely:
when |l==single_value| an actual empty tuple should be produced, which
|wrap_tuple<0>()| does.

@< Function definitions @>=
void voiding::evaluate(level l) const
{@; exp->void_eval();
  if (l==single_value)
    wrap_tuple<0>();
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
{ for (auto it=coerce_table.begin(); it!=coerce_table.end(); ++it)
    if (from_type==*it->from and to_type==*it->to)
    @/{@; e.reset(new conversion(*it,std::move(e)));
      return true;
    }
  if (to_type==void_type)
  {@; e.reset(new voiding(std::move(e)));
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
A call to |conform_types| will beinto an invariably followed (upon success) by
returning the expression now held in the argument~|d| from the calling
function; we can avoid having to repeat that argument in a return statement by
returning the value in question already from |conform_types|. To indicate that
the expression~|d| is incorporated into to return value, we choose to get |d|
passed by rvalue reference, even though the argument will usually be held in a
variable.

If both attempts to conform the types fail we throw a |type_error|; doing so
we must take a copy of |found| (since it a qualified |const|), but we can move
from |required|, whose owner will be destructed before the error is caught.

@< Function def... @>=
expression_ptr conform_types
(const type_expr& found, type_expr& required, expression_ptr&& d, const expr& e)
{ if (not required.specialise(found) and not coerce(found,required,d))
    throw type_error(e,found.copy(),std::move(required));
  return std::move(d);
}

@ List displays and loops produce a row of values of arbitrary (but identical)
type; when they occur in a context requiring a non-row type, we may be able to
find a coercion that reconciles the requirements. The following function finds
whether this is possible, and if so sets |components_type| to the type
required for the components; if not it will return~|nullptr|. The
implementation is simply to look in |coerce_table| for conversions from some
row type, then try to specialise |component_type| to the component type of
that row type. By using this mechanism, the components of a row display
themselves obtain a context that may generate further conversions to obtain
this type.

Currently all calls to this function have |component_type| initially
undetermined, so the call to of the |specialise| method will always succeed,
but we test the result nonetheless. In case there exist multiple row types
that could convert to |final_type|, the first one in the table is chosen; this
currently happens when |final_type| is \.{mat}, in which case this function
will return |component_type| equal to \.{vec} rather than to \.{[int]}.

@< Function def... @>=
const conversion_record* row_coercion(const type_expr& final_type,
                                            type_expr& component_type)
{ for (auto it=coerce_table.begin(); it!=coerce_table.end(); ++it)
    if (final_type==*it->to and it->from->kind==row_type)
      return component_type.specialise(*it->from->component_type)
        ? &*it
        : nullptr;
  return nullptr;
}

@*1 Proximity of types.
%
The above implicit conversions of types pose a limitation to the possibilities
of overloading operator symbols and function identifiers. If a symbol should
be overloaded for too closely related operand types, situations could occur in
which given operand expressions can be converted to either one of the operand
types. This would either produce unpredictable behaviour, or necessitate a
complicated set of rules to determine which of the definitions of the symbol
is to be used (overloading resolution in \Cpp\ is a good example of such
complications).

There are two ways to avoid the occurrence of complications by restricting the
rules of the language: either forbid type conversions in arguments of
overloaded symbols, or forbid simultaneous definitions of such symbols for too
closely related types. Forbidding all automatic type conversions in case of
overloading would defeat to a large extent the purpose of overloading, namely
as a convenience to the user. Originally we therefore chose the latter
solution of allowing all coercions in arguments, while forbidding overloading
in certain cases. However this both had practical implementation problems (we
needed to try to convert operands with every possible operand type as required
type, which took too much time) and led to severe mutual exclusions between
overloaded types, as often some expression might be simultaneously acceptable
to two different types; as extreme case we needed do exclude \.{void} as
overloaded operand type altogether, since \emph{any} valid expression can be
voided to \.{void}. Therefore we settled for a mixture of both kinds of
restrictions, allowing only certain conversions in operand types and
forbidding overloading between types related by them. This is somewhat along
the model of the language Algol~68, where the ``firm'' context for operands
allows a subset of coercions. Our approach involves analysing operands twice:
first in isolation to determine their \foreign{a priori} type, and then
possibly a second time with a selected overload in the context of the required
operand type (if different from the \foreign{a priori} type), with the
occasion to insert coercions as needed.

Our rules for coercions and overloading will be governed by a single relation
|is_close| between pairs of (argument) types: its resulting value (a small
integer) will tell both whether for a given \foreign{a priori} type another
(operand) type provides a viable candidate, and whether two types can coexist
as operand types for a same overloaded operator or function. (Multiple
operands or arguments are considered as one argument with the tuple type
formed from their individual types.)

@< Declarations of exported functions @>=
unsigned int is_close (const type_expr& x, const type_expr& y);

@ We do allow simultaneous overloading between closely related types in some
cases, namely if they can be ordered so that if $t_1$ precedes $t_2$ then some
expression of (\foreign{a priori}) type~$t_1$ can be converted to type~$t_2$
but no expression of type~$t_2$ can be converted to type~$t_1$; in such cases
reasonable behaviour can be obtained by trying a match for~$t_1$ before trying
one for~$t_2$, and this allows for instance arithmetic operators to be defined
for type \.{(int,int)} as well as for type \.{(rat,rat}). Note that the
conversion does not necessarily apply at the outer level, since the
expressions could be tuple or row displays, with coercions being applied to
individual component expressions; for instance there exists nested displays of
type \.{([vec],[int])} that can be converted to type \.{([[int]],[rat])}
or \.{(mat,ratvec)}. So our relation will be a partial order, and compatible
with tuple and row formation: $x_i\leq y_i$ for all~$i$ implies
$(x_1,\ldots,x_n)\leq(y_1,\ldots,y_n)$ as well as $[x_1]\leq[y_1]$.

Therefore |is_close| returns a value composed of 3~bits: one indicating
whether the types are close at all, and two others for indicating the
conversion partial order in one direction or the other. Thus |is_close|
returns an integer rather than a boolean value, with the following
interpretations: |0x0| means the types are unrelated, |0x4| means the types
are mutually exclusive but neither can be converted to the other (as for
instance \.{(int,rat)} and \.{(rat,int)}), |0x5| means the types are close and
(only) the first can be converted to the second (example, \.{int}
and \.{rat}), |0x6| is the opposite relation, and |0x7| means both types can
be converted to each other (like \.{vec} and \.{[int]}, or any case of equal
types).

For types $t_1$ and~$t_2$ which do not admit a relative priority, we want to
disallow simultaneous overloading with arguments types $t_1$ and~$t_2$ if any
expression given as argument could be converted to either of them. Deciding
the existence of such an expression would require study of all available
language constructs, but the situation is somewhat simplified by the fact
that, for efficiency reasons, overloading resolution is not done using the
complete argument expression, but only using its type. In fact matching will
be done using calls to the very function |is_close| we are discussing here,
testing the bit for conversion towards the required argument type; this
provides us with an opportunity to adjust rules for possible type conversions
of arguments at the same time as defining the exclusion rules.

Empty row displays, or more precisely arguments of type~\.{[*]}, pose a
difficulty: they would be valid in any context requiring a specific row type,
so if we stipulated that one may write \.{[]} to designate an empty row
operand of any row type, then |is_close| would have to consider all row types
close to each other (and therefore mutually exclusive for overloading). This
used to be the convention adopted, but it was found to be rather restrictive
in use, so the rules were change to state that an un-cast expression \.{[]}
will not match overload instances of specific row types; it might match an
parameter of specified type \.{[*]} (once we allow that as type expression),
and such a parameter could \emph{only} take an empty list corresponding
argument (which is very limiting of course, but it allows being explicit about
which overloaded instance should be selected by an argument \.{[]}).

These considerations are not limited to empty lists (although it is the most
common case): whenever an expression has an \foreign{a priori} type
containing \.*, that expression will not select any overload with a concrete
type in its place (overloading does not perform type specialisation).


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
{ expression_ptr dummy(nullptr);
  if (x==y)
    return 0x7;
  if (x.kind==undetermined_type or y.kind==undetermined_type)
    return 0x0;
      // undetermined types do not specialise (or coerce), and are not close
  if (x.kind==primitive_type or y.kind==primitive_type)
  { unsigned int flags=0x0;
    if (coerce(x,y,dummy)) flags |= 0x1;
    if (coerce(y,x,dummy)) flags |= 0x2;
    return flags==0 ? flags : flags|0x4;
  }
  if (x.kind!=y.kind)
    return 0x0;
  if (x.kind==row_type)
    return is_close(*x.component_type,*y.component_type);
  if (x.kind!=tuple_type)
    return 0x0; // non-aggregate types are only close if equal
  auto it0=x.tupple.begin(), it1=y.tupple.begin();
  unsigned int flags=0x7;
  while (not x.tupple.at_end(it0) and not y.tupple.at_end(it1) @|
         and (flags&=is_close(*it0,*it1))!=0)
  @/{@; ++it0; ++it1; }
  return x.tupple.at_end(it0) and y.tupple.at_end(it1) ? flags : 0x0;
}

@* Error values.
%
Before we describe evaluation of expressions we must realise that evaluation
can cause runtime errors. The evaluator may throw exceptions due to
inconsistency of our (rather than the user's) program, which are classified as
|std::logic_error|. It may also throw exceptions due to errors not caught by
the type checker in the user input, such as size mismatch in matrix
operations; in such cases it will throw |std::runtime_error|. These standard
exception types will be used without any type derivation.

@< Includes needed in \.{types.h} @>=
#include <stdexcept>

@ We first derive a general exception class |program_error| from
|std::exception|, used for all errors other than runtime errors. It
encompasses all kind of error of the input detected before evaluation starts
by static analysis; for instance the use of undefined variables falls in this
category. The derived class just stores an error message string.

@< Type definitions @>=
class program_error : public std::exception
{ std::string message;
public:
  explicit program_error(const std::string& s) : message(s) @+{}
  ~program_error () noexcept @+{} // backward compatibility for gcc 4.6
  const char* what() const throw() @+{@; return message.c_str(); }
};

@ We derive from |program_error| an exception type |expr_error| that stores in
addition to the error message a reference to an expression to which the
message applies. Placing a reference in an error object may seem hazardous,
because the error might terminate the lifetime of the object referred to, but
in fact it is safe: all |expr| objects are constructed in dynamic memory
during parsing, and destructed at the disposal of the now translated
expression at the end of the main interpreter loop; all throwing of
|expr_error| (or derived types) happens after the parser has finished, and the
corresponding |catch| happens in the main loop before disposal of the
expression, so the reference certainly survives the lifetime of the
|expr_error| object.

The error type is declared a |struct|, so that the |catch| clause may access
the |offender| expression for use in an error message. At the point where this
error value is constructed (and thrown), we just provide a general error
message |s| for storage in the |program_error| base class, and the offending
expression. There is no need to override either of the virtual methods here.

@< Type definitions @>=
struct expr; // predeclare
struct expr_error : public program_error
{ const expr& offender; // the subexpression causing a problem
@)
  expr_error (const expr& e,const std::string& s) throw()
    : program_error(s),offender(e) @+{}
};

@ We derive from |expr_error| an even more specific exception type
|type_error| that we shall throw when an expression fails to have the right
type. In addition to the offending subexpression it stores the
types that failed to match. For the latter, the error value owns the types
pointed to, so the caller should relinquish ownership of, or copy (in case it
did not own), the types passed when throwing the exception.

@< Type definitions @>=
struct type_error : public expr_error
{ type_expr actual; @+
  type_expr required; // the types that conflicted
@)
  type_error (const expr& e, type_expr&& a, type_expr&& r) noexcept @/
    : expr_error(e,"Type error") @|
      ,actual(std::move(a)),required(std::move(r)) @+{}
  type_error@[(type_error&& e) = default@];
};

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
Cartan_class_type, KGB_element_type, block_type,
module_parameter_type, split_integer_type, virtual_module_type, @[@]


@~The following list must match that of the previous module, for proper
functioning of I/O.

@< Other primitive type names @>=
"vec", "mat", "ratvec",
"LieType","RootDatum", "InnerClass", "RealForm",
"CartanClass", "KGBElt", "Block", "Param", "Split", "ParamPol", @[@]

@* Index.

% Local IspellDict: british
