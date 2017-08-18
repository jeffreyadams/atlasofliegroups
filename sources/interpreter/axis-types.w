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

@( axis-types.h @>=

#ifndef AXIS_TYPES_H
#define AXIS_TYPES_H

@< Includes needed in \.{axis-types.h} @>@;
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

@h "axis-types.h"
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


@< Includes needed in \.{axis-types.h} @>=
#include <memory> // for |std::unique_ptr|, |std::shared_ptr|
#include "../Atlas.h" // for utilities (like |sl_list|); must come first
#include "buffer.h" // for |id_type|


@* Types to represent types.
%
We make a fundamental choice to check types before attempting to execute any
expression entered by the user; thus type errors can be signalled earlier and
expressed more understandably than if one would wait until an accident is
about to happen at runtime. This means that types are to be represented
independently of runtime values. Types will have an effect on the
transformation of expressions into executable internal form (for example in
selecting the instance of an overloaded function to use), but once that is
done, types will no longer be represented at all: the internal code executes
without being ``aware'' of the types of the values it manipulates. This is
only possible because those values are then accessed via generic pointers
(i.e., pointers that, as far as \Cpp\ is concerned, could point to any kind of
value). Currently runtime values also contain an indication of their type,
which will allow us to do some dynamic testing to trap possible errors in the
logic of our interpreter, but the user should never notice this; in an
optimised version of the \.{atlas} program, this runtime type information and
corresponding tests could be dropped altogether.

Types are represented by small tree structures, with nodes of types
|type_expr| that will be detailed later. For now we provide no user-defined
encapsulation (as \Cpp~classes do), so a type expression basically describes
the structure of the corresponding values, and allows all operations
compatible with that structure to be performed. For instance the type printed
as \.{(int->[(bool,mat)])} specifies a function mapping integers to lists of
pairs of a Boolean value and a matrix, and allow operations compatible with
that such as calling it with an integer argument and then selecting the second
component of the final element of the result. However, classes defined in the
Atlas software library itself will be presented to the user as (opaque)
primitive types; only built-in operations declared for that primitive type can
be applied to such values, so we do provide abstraction from internal
representation there. One kind of type is a function type, which specifies the
types of the individual function arguments and of the values it returns.

Trees representing different type expressions will not share any sub-trees
among each other, so they have strict (unshared) ownership of their parts,
which simplifies memory management. This choice means some recursive copying
of tree structures is sometimes required, but (although runtime cost of type
handling is negligible with respect to other factors in the interpreter) we
avoid doing so more than absolutely necessary.

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
also uses ordinary pointers of type~|type_p| when ownership is managed by the
containing structure.

@< Type definitions @>=
class type_expr;
typedef type_expr* type_p;
typedef const type_expr* const_type_p;
typedef std::unique_ptr<type_expr> type_ptr;

@*2 Type lists.
%
We also need type lists as building block for types (for instance for the
arguments of a function). These used to be defined in a similar manner to
types themselves, but since an STL-compatible singly-linked list container
class templates |simple_list| and |sl_list| was added to the Atlas utilities
library, these were used to replace the implementation of type lists.

@< Includes needed in \.{axis-types.h} @>=
#include "sl_list.h" // for internals of the |sl_list| class template

@ The class template |simple_list| is built into the structure itself, while
|sl_list|, which is more flexible but requires more space at the head of the
list (the class |sl_list<type_expr>| itself), is occasionally used for
temporary variables. This provides a good test for the usability of the new
container type; so far it has passed them gracefully, though occasionally
after enriching the repertoire of methods of the container type.

@< Type definitions @>=
typedef containers::simple_list<type_expr> type_list;
typedef atlas::containers::sl_node<type_expr>* raw_type_list;
typedef containers::sl_list<type_expr> dressed_type_list;
typedef containers::weak_sl_list_const_iterator<type_expr> wtl_const_iterator;
typedef containers::weak_sl_list_iterator<type_expr> wtl_iterator;
  // wtl = weak type list

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
for the list itself. Nonetheless, the latter is passed as a modifiable lvalue
reference |dst|, as it will continue to hold the extended list.

@< Declarations of exported functions @>=
type_list empty_tuple();
type_list& prefix(type_expr&& t, type_list& dst);
dressed_type_list& prefix(type_expr&& t, dressed_type_list& dst);

@ Move semantics, introduced in \Cpp11, has solved any difficulty with
ownership management that the functions below previously had to deal with.

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
Furthermore one has function types that map one type to another, types that
are ``row of'' some other type, tuples (Cartesian products) of some given
sequence of types, and disjoint unions (co-products in the category of sets)
of some sequence of types. We also allow for an undetermined type, which can
serve as a wild-card to specify type patterns. Before we can define
|type_expr|, we need to enumerate its variants and the possibilities for
primitive types. Here are enumerations of tags for the basic kinds of types,
and for the sub-cases of |primitive_type| (some of them will be introduced
later).

@< Type definitions @>=
enum type_tag
 { undetermined_type, primitive_type,
   function_type, row_type, tuple_type, union_type };

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
Type expressions are defined by a tagged union. Originally clients could
access the tag and the corresponding variant directly, since coherent variant
usage (accessing only the variant indicated by the tag) is not something that
accessor methods could enforce anyway. However we have none nonetheless made
the fields private with public accessor methods, in order to be able to hide
in our implementation possibilities not visible to clients. Type expressions
are implemented as unshared linked structures (so they are owned by the top
level object) that used to always spell out the type structure in full, but we
want to also be able to represent types with a recursive structure, where
spelling out in full would lead to infinite structures so is not possible.
Anticipating this, we allow for our class to implement an abbreviation
mechanism in a way transparent to clients, which is possible by exposing the
tag and type components only through accessor methods. These methods do return
the stored fields unchanged: non-owning raw pointers into the type expression

The field |tuple_variant| used to be (after inclusion of the
|containers::simple_list| class template in Atlas, which coincided with the
migration to the \Cpp11 standard) of type |type_list|. Having variant members
of a |union| with nontrivial special member functions is allowed in \Cpp11,
although it remains the programmer's responsibility to explicitly call
constructors and destructors as those variants come and go. Currently this is
changed to |raw_type_list|; the main difference is that the variant can now be
initiated by simple assignment rather than placement |new|, and at termination
requires a call of |delete| for the raw pointer rather than an explicit
destructor call for the |type_list| object. This step backwards to raw
pointers is therefore mostly a simplification, but there is one real drawback
of not having a true |type_list| object: it is now impossible to create a
|type_list::iterator| (or its |const| relative) to iterate over the type lists
in tuple types. It turns out however that most of the time weak iterators
(which do not allow for insertion of deletion of nodes) are sufficient, which
is why the |raw_type_list| solution is now chosen.

There is one restriction on types that is not visible in the definition below,
namely that the list of types referred to by the |tuple_variant| field cannot
have length~$1$ (but length~$0$ is allowed, the $0$-tuple being the |void|
type; we might also take the $0$-union to be the uninhabited type \.* that is
for instance the ``return type'' of the |error| function, but this is not
currently the case). This is because anything that would suggest a $1$-tuple
or $1$-union (for instance a parenthesised expression) is identified with its
unique component.

@< Definition of |type_expr| @>=
class type_expr
{ type_tag tag;
  union
  { primitive_tag prim_variant; // when |tag==primitive|
    func_type* func_variant; // when |kind==function_type|
    type_p compon_variant; // when |kind==row_type|
    raw_type_list tuple_variant; // when |kind==tuple_type| or |kind==union_type|
  };
@)
public:
  type_tag kind () const @+{@; return tag; }
  func_type * func() const @+{@; return func_variant; }
  type_p component_type () const @+{@; return compon_variant; }
  primitive_tag prim() const @+{@; return prim_variant; }
  raw_type_list tuple () const @+{@; return tuple_variant; }
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
require, acquire, or lose ownership of the pattern for/by this method.

The constructors for the row, tuple, and union types receive pointers that
will be directly inserted into our |struct| and become owned by it; the
ownership management is simplest when such pointers are ``passed by check'',
i.e., as rvalue references to smart pointers.

@< Methods of the |type_expr| structure @>=

type_expr() noexcept : tag(undetermined_type) @+{}
explicit type_expr(primitive_tag p) noexcept
  : tag(primitive_type), prim_variant(p) @+{}
inline type_expr(type_expr&& arg, type_expr&& result) noexcept;
 // for function types
explicit type_expr(type_ptr&& c) noexcept
  : tag(row_type), compon_variant(c.release()) @+{}
explicit type_expr(type_list&& l,bool is_union=false) noexcept;
  // tuple and union types
@)
#ifdef incompletecpp11
type_expr(const type_expr& t) = @[delete@];
type_expr& operator=(const type_expr& t) = @[delete@];
#endif
type_expr(type_expr&& t) noexcept; // move constructor
type_expr& operator=(type_expr&& t) noexcept; // do move assignment only
void clear() noexcept;
  // resets to undefined state, cleaning up owned pointers
~type_expr() noexcept @+{@; clear(); }
  // that is all that is needed for destruction
void swap(type_expr& t) noexcept;
@)
type_expr copy() const; // in lieu of deep copy constructor
void set_from(type_expr&& p) noexcept; // shallow copy
bool specialise(const type_expr& pattern);
  // try to match pattern, possibly modifying |*this|
bool can_specialise(const type_expr& pattern) const;
  // tell whether |specialise| would succeed

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

@ The constructor for the |type_list| variant with |tuple_variant| field is a
move constructor.

@< Function definitions @>=
type_expr::type_expr(type_list&& l,bool is_union) noexcept
  : tag(is_union ? union_type: tuple_type)
  , tuple_variant(l.release())
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
  switch (result.tag=tag)
  { case undetermined_type: break;
    case primitive_type: result.prim_variant=prim_variant; break;
    case row_type:
      result.compon_variant=new type_expr(compon_variant->copy());
    break;
    case tuple_type: case union_type:
      @< Assign a deep copy of |tuple_variant| to |result.tuple_variant| @>
    break;
    case function_type: result.func_variant=new
    func_type(func_variant->copy());
    break;
  }
  return result;
}

@ We create a duplicate list by applying the |copy| method for all members of
the list accessed from |tuple_variant|. The code below originally
placement-constructed that list into the |result.tuple_variant| field,
initially default constructing the list, and then adding the components using
|insert| through an output iterator |oit| pointing at the end of the list.
This had to be modified when |tuple_variant| was made to by a raw pointer,
from which we cannot create an |type_list::iterator| (this is one of the rare
places in the code where that matters). So now we instead create an initially
empty |type_list dst@;| and make an iterator for it; then after creating the
nodes of the list, the |type_list| is demoted to |raw_type_list| by calling
its |release| method, which raw pointer is assigned to |result.tuple_variant|.

@< Assign a deep copy of |tuple_variant| to |result.tuple_variant| @>=
{
  wtl_const_iterator it(tuple_variant);
  type_list dst;
  for (type_list::iterator oit = dst.begin(); not it.at_end(); ++it,++oit)
    dst.insert(oit,it->copy());
  result.tuple_variant = dst.release(); // incorporate and transfer ownership
}

@ The method |clear|, doing the work for the destructor, must similarly clean
up afterwards, with the recursion again being implicit. We set |tag =
undetermined_type| to indicate that no variant is active; this condition is
tested in the |set_from| method (and also in by the destructor, though if
|clear| is \emph{called by} the destructor, that test is already behind us).

There is no reason to use a smart pointers as variant of a union, since one
must still explicitly propagate destruction to the active variant, as happens
here for the |tuple_variant| field. If it were done here, then all that would
change is that the code below would have to call their destructors rather than
call |delete| directly.

@< Function definitions @>=
void type_expr::clear() noexcept
{ switch (tag)
  { case undetermined_type: case primitive_type: break;
    case row_type: delete compon_variant; break;
    case tuple_type: case union_type: delete tuple_variant; break;
    case function_type: delete func_variant; break;
  }
  tag = undetermined_type;
}

@ The method |set_from| makes a shallow copy of the structure; it implements
move semantics. Correspondingly it takes as argument another |type_expr| by
modifiable rvalue reference. The contents of its top-level structure will be
moved to the current |type_expr| using move assignment where applicable, and
then set to a empty |undetermined_type| value, which effectively detaches any
possible descendants from it. This operation requires that |type_expr|
previously had |tag==undetermined_type|. We test this condition using an
|assert| statement (rather than throwing |logic_error|) to honour the
|noexcept| specification.

The moving copy constructor is nearly identical to the |set_from| method.

@h <stdexcept>
@< Function definitions @>=
void type_expr::set_from(type_expr&& p) noexcept
{ assert (tag==undetermined_type);
   // logic should ensure this; we promised not to |throw|
  switch(tag=p.tag) // copy top node
  { case undetermined_type: break;
    case primitive_type: prim_variant=p.prim_variant; break;
    case function_type: func_variant=p.func_variant; break;
    case row_type: compon_variant=p.compon_variant; break;
    case tuple_type: case union_type: tuple_variant = p.tuple_variant; break;
  }
  p.tag=undetermined_type;
  // detach descendants, so |p.clear()| will destroy top-level only
}
@)
type_expr::type_expr(type_expr&& x) noexcept // move constructor
: tag(x.tag)
{ switch(tag) // move top node
  { case undetermined_type: break;
    case primitive_type: prim_variant=x.prim_variant; break;
    case function_type: func_variant=x.func_variant; break;
    case row_type: compon_variant=x.compon_variant; break;
    case tuple_type: case union_type: tuple_variant = x.tuple_variant; break;
  }
  x.tag=undetermined_type;
  // detach descendants, so destructor of |x| will do nothing
}

@ For move assignment we reuse the |set_from| method, and for |swap| we do a
move construction and two |set_from| calls, unless the |tag| fields match, in
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
  if (tag==other.tag) // an easy case
    switch(tag)
    { case undetermined_type: break; // no need to swap |nothing| fields
      case primitive_type: std::swap(prim_variant,other.prim_variant); break;
      case function_type: std::swap(func_variant,other.func_variant); break;
      case row_type: std::swap(compon_variant,other.compon_variant); break;
      case tuple_type: case union_type: std::swap(tuple_variant,other.tuple_variant); break;
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
{ if (pattern.tag==undetermined_type)
    return true; // specialisation to \.* trivially succeeds.
  if (tag==undetermined_type) // specialising \.* also always succeeds,
    {@; set_from(pattern.copy()); return true; }
     // by setting |*this| to |pattern|
  if (pattern.tag!=tag) return false; // impossible to refine
  switch(tag)
  { case primitive_type: return prim_variant==pattern.prim_variant;
    case function_type:
      return func_variant->arg_type.specialise(pattern.func_variant->arg_type) @|
         and func_variant->result_type.specialise
                                          (pattern.func_variant->result_type);
    case row_type: return compon_variant->specialise(*pattern.compon_variant);
    case tuple_type: case union_type:
     @< Try to specialise types in |tuple_variant| to those in
        |pattern.tuple_variant|,
        and |return| whether this succeeded @>
    default: return true; // to keep the compiler happy, cannot be reached
  }
}

@ For tuples and unions, specialisation is done component by component.

@< Try to specialise types in |tuple_variant| to those in
   |pattern.tuple_variant|... @>=
{
  wtl_iterator it0(tuple_variant);
  wtl_const_iterator it1(pattern.tuple_variant);
  while (not it0.at_end() and not it1.at_end() and it0->specialise(*it1))
    @/{@; ++it0; ++it1; }
  return it0.at_end() and it1.at_end();
  // whether both lists terminated
}

@ Here is the definition of the accessor |can_specialise|, which is quite
similar.

@< Function definitions @>=
bool type_expr::can_specialise(const type_expr& pattern) const
{ if (pattern.tag==undetermined_type)
    return true; // specialisation to \.* trivially succeeds.
  if (tag==undetermined_type)
    return true;
  if (pattern.tag!=tag) return false; // impossible to refine
  switch(tag)
  { case primitive_type: return prim_variant==pattern.prim_variant;
    case function_type:
      return func_variant->arg_type.can_specialise
                                    (pattern.func_variant->arg_type) @|
         and func_variant->result_type.can_specialise
                                    (pattern.func_variant->result_type);
    case row_type:
      return compon_variant->can_specialise(*pattern.compon_variant);
    case tuple_type: case union_type:
      @< Find out and |return| whether we can specialise the types in
         |tuple_variant|
         to those in |pattern.tuple_variant| @>
    default: return true; // to keep the compiler happy, cannot be reached
  }
}

@ For tuples and unions, the test for possible specialisation is done
component by component.

@< Find out and |return| whether we can specialise the types in
   |tuple_variant| to those in |pattern.tuple_variant| @>=
{
  wtl_const_iterator it0(tuple_variant);
  wtl_const_iterator it1(pattern.tuple_variant);
  while (not it0.at_end() and not it1.at_end()
         and it0->can_specialise(*it1))
    @/{@; ++it0; ++it1; }
  return it0.at_end() and it1.at_end();
  // whether both lists terminated
}

@ The constructor for function types cannot be defined inside the structure
definition, since |func_type| is not (and cannot be) a complete type there.
The constructor for |func_type| used will be defined below.

@< Template and inline function definitions @>=
type_expr::type_expr(type_expr&& arg, type_expr&& result) noexcept
@/ : tag(function_type)
   , func_variant(new func_type(std::move(arg),std::move(result)))
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
constructor for~|type_expr|. Like for |type_expr| we do not provide an
ordinary copy constructor, but a value-returning |copy| accessor method.

@s result_type normal

@< Type definitions @>=
struct func_type
{ type_expr arg_type, result_type;
@)
  func_type(type_expr&& a, type_expr&& r)
@/ : arg_type(std::move(a)), result_type(std::move(r)) @+{}
#ifdef incompletecpp11
  func_type(const func_type& f) = @[delete@];
  func_type& operator=(const func_type& f) = @[delete@];
  func_type(func_type&& f)
@/: arg_type(std::move(f.arg_type)), result_type(std::move(f.result_type))
  @+{}
  func_type& operator=(func_type&& f)
  {@; arg_type = std::move(f.arg_type); result_type = std::move(f.result_type);
    return *this;
  }
#else
  func_type(func_type&& f) = @[default@]; // move constructor
  func_type& operator=(func_type&& f) = @[default@]; // move assignment
#endif
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
tuple or union types; defining a separate function for a |type_list|
facilitates our task a bit.

@< Function definitions @>=

void print(std::ostream& out, const raw_type_list& l,char sep)
{ wtl_const_iterator it(l);
  if (not it.at_end())
    while (out << *it, not (++it).at_end())
      out << sep;
}

std::ostream& operator<<(std::ostream& out, const func_type& f)
{
  out << '(';
  if (f.arg_type.kind()==tuple_type or
      f.arg_type.kind()==union_type and f.arg_type.tuple()!=nullptr)
     print(out,f.arg_type.tuple(),f.arg_type.kind()==tuple_type?',':'|');
     // naked tuple or union
  else out << f.arg_type; // other component type
  out << "->";
  if (f.result_type.kind()==tuple_type or
      f.result_type.kind()==union_type and f.result_type.tuple()!=nullptr)
     print(out,f.result_type.tuple(),f.result_type.kind()==tuple_type?',':'|');
     // naked tuple or union
  else out << f.result_type; // other component type
  out << ')';
  return out;
}
@)
std::ostream& operator<<(std::ostream& out, const type_expr& t)
{ switch(t.kind())
  { case undetermined_type: out << '*'; break;
    case primitive_type: out << prim_names[t.prim()]; break;
    case function_type: out << *t.func(); break;
    case row_type: out << '[' << *t.component_type() << ']'; break;
    case tuple_type:
      if (t.tuple()==nullptr)
        out << "void";
      else
      {@;
         print(out << '(', t.tuple(),',');
         out << ')';
      }
    break;
    case union_type:
      if (t.tuple()==nullptr)
        out << "(*)";
      else
      {@;
         print(out << '(', t.tuple(),'|');
         out << ')';
      }
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
{ if (x.kind()!=y.kind()) return false;
  switch (x.kind())
  { @+ default:
// all cases are listed below, but compilers want a |default| to |return|
  @\case undetermined_type: return true;
    case primitive_type: return x.prim()==y.prim();
    case function_type:
      return x.func()->arg_type==y.func()->arg_type
	 and x.func()->result_type==y.func()->result_type;
    case row_type: return *x.component_type()==*y.component_type();
    case tuple_type: case union_type:
       @< Find out and |return| whether all types in |x.tuple()|
          are equal to those in |y.tuple()| @>
  }
}

@ This module has a familiar structure.
@< Find out and |return| whether all types in |x.tuple()| are equal to those in
  |y.tuple()| @>=
{
  wtl_const_iterator it0(x.tuple());
  wtl_const_iterator it1(y.tuple());
  while (not it0.at_end() and not it1.at_end()
         and *it0==*it1)
    @/{@; ++it0; ++it1; }
  return it0.at_end() and it1.at_end();
  // whether both lists terminated
}

@*2 Type constructing functions.
%
Instead of using the constructors directly, we often use the constructing
functions below. The versions with |mk| return |type_ptr| values owning the
constructed expression. They also take such smart pointers, or |type_list|
values, as argument whenever the underlying pointer is to be directly inserted
into the structure built. However for |mk_function_type| this is not the case,
so instead of insisting that the caller hold a unique-pointer to the argument
types it suffices to hold a |type_expr| whose contents can be moved into the
type to be constructed. If |t| is a |type_ptr| held by a client, it can pass
|std::move(*t)| as argument to |mk_function_type|, which will move the
contents of the node |t| points to into the function type, and after return
the destructor of |t| will eventually delete the now empty node.

The functions with |make| are wrappers that do the same thing, but undress the
smart pointers in their interface to raw pointers; these functions
should be used exclusively in the parser. In other words these raw pointers
are owning, as if they were smart pointers, and in fact the functions convert
their arguments into unique pointers right away, and produce their results by
calling the |release| method of a unique pointer. This travesty is done just
so that the parser only sees POD, types which it can put into a |union|.

@< Declarations of exported functions @>=
type_ptr mk_prim_type(primitive_tag p);
type_ptr mk_function_type(type_expr&& a, type_expr&& r);
type_ptr mk_row_type(type_ptr&& c);
type_ptr mk_tuple_type (type_list&& l);
type_ptr mk_union_type (type_list&& l);
@)
type_p make_prim_type(unsigned int p);
type_p make_function_type(type_p a,type_p r);
type_p make_row_type(type_p c);
type_p make_tuple_type(raw_type_list l);
type_p make_union_type(raw_type_list l);
@)
raw_type_list make_type_singleton(type_p raw);
raw_type_list make_type_list(raw_type_list l,type_p t);

@ The functions like |mk_prim| below simply call the constructor in the
context of the operator |new|, and then capture of the resulting (raw) pointer
into a |type_ptr| result. The functions like |make_prim| wrap them for the
parser interface, undressing the pointers to raw again. Using this setup it is
ensured that all pointers are considered owning their target, implicitly so
while being manipulated by the parser (which guarantees that every pointer
placed on the parsing stack will be argument of an interface function exactly
once, possibly some |destroy| function in case it pops symbols during error
recovery).

The first group of functions makes some provisions to facilitate special
relations between types, so that the parser does not have to deal with them.
Thus |mk_prim_type| will return an empty tuple type when called with the type
name for |"void"|, although this is not a primitive type (indeed |"void"|
should probably better be handled just like user-defined type abbreviations).
The function |mk_union_type| will not encapsulate a list of length one into an
invalid |type_expr|, but rather return its unique list element, unpacked. The
reason for this is that having a syntactic category for a list of at least one
component separated by vertical bars simplifies the grammar considerably; when
the list has one component (and no bars) a union of one variant is temporarily
created, but upon incorporation into an encompassing type, the union is then
removed again. For tuples a similar provision is not necessary, as a somewhat
more involved set of grammar rules is used that avoids making type lists of
length~$1$.

The functions |make_tuple_type| and |make_union_type| reverse the list of
types they handle, to compensate for the fact that the left-recursive grammar
rules (easier for the parser generator) construct their type lists in reverse
order.

@< Function definitions @>=
type_ptr mk_prim_type(primitive_tag p)
{ return p<nr_of_primitive_types ?
    type_ptr(new type_expr(p)) :
    type_ptr(mk_tuple_type(empty_tuple()));
}

type_ptr mk_function_type (type_expr&& a, type_expr&& r)
{@; return type_ptr(new type_expr(std::move(a),std::move(r))); }

type_ptr mk_row_type(type_ptr&& c)
{@; return type_ptr (new type_expr(std::move(c))); }

type_ptr mk_tuple_type (type_list&& l)
{@; return type_ptr(new type_expr(std::move(l),false)); }

type_ptr mk_union_type (type_list&& l)
{ return type_ptr(l.singleton() ? new type_expr(std::move(l.front()))
                  : new type_expr(std::move(l),true));
}

@)
type_p make_prim_type(unsigned int p)
{@; return mk_prim_type(static_cast<primitive_tag>(p)).release(); }

type_p make_function_type(type_p a,type_p r)
{@; return
    mk_function_type(std::move(*type_ptr(a)),std::move(*type_ptr(r))).release();
}

type_p make_row_type(type_p c) {@; return mk_row_type(type_ptr(c)).release(); }

type_p make_tuple_type(raw_type_list l)
{@; type_list result(l); result.reverse();
    return mk_tuple_type(std::move(result)).release();
}

type_p make_union_type(raw_type_list l)
{@; type_list result(l); result.reverse();
    return mk_union_type(std::move(result)).release();
}
@)

raw_type_list make_type_singleton(type_p t)
{@; type_list result;
  result.push_front(std::move(*type_ptr(t)));
  return result.release();
}

raw_type_list make_type_list(raw_type_list l,type_p t)
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
denotations in the source code (mostly in calls installing built-in functions
into \.{axis}) rather than from user input, and we are not going to write
incorrect strings (we hope). Therefore we don't care if the error handling is
crude here. The simplest way of parsing ``by hand'' is recursive descent, so
that is what we shall use. By passing a character pointer by reference, we
allow the recursive calls to advance the index within the string read.

The function |scan_type| does the real parsing, |mk_type_expr| calls it,
providing a local modifiable pointer to bind to its reference parameter (which
is important because |scan_type| cannot directly accept a \Cee-string constant
as argument) while also doing error reporting, and |mk_type| is just a wrapper
around |mk_type_expr| that converts the result from a |type_expr| to a smart
pointer to (a freshly allocated instance of) such. Currently the function
|mk_type_expr| is called only during the start-up phase of \.{atlas}, and if an
error is encountered (of type |logic_error|, since this must be an error in
the \.{atlas} program itself), printing of the error message will be followed
by termination of the program.

@< Function definitions @>=
type_expr scan_in_parens(const char*& s);
type_expr scan_union_list(const char*& s);
type_expr scan_tuple_list(const char*& s);
@)
type_expr scan_type(const char*& s)
{ if (*s=='(')
  { type_expr result=scan_in_parens(++s);
    if (*s++!=')')
      throw logic_error("Missing ')' in type");
    return result;
  }
  else if (*s=='[')
  {
    type_ptr p(new type_expr(scan_in_parens(++s)));
    if (*s++!=']')
      throw logic_error("Missing ']' in type");
    return type_expr(std::move(p));
  }
  else if (*s=='*') return ++s,@[type_expr()@]; // undetermined type
  else @< Scan and |return| a primitive type, or |throw| a |logic_error| @>
}
@)
type_expr mk_type_expr(const char* s)
{ const char* orig=s;
  try
  {@; return scan_type(s); }
  catch (logic_error e)
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


@ For primitive types we use the same strings as for printing them. We test as
many characters as the type name has, and the fact that no alphanumeric
character follows, so that a longer type name will not match a prefix of it.

In this module we use the fact that the order in the list |prim_names| matches
that in the enumeration type |primitive_tag|, by casting the integer index
into the former list to an element of that enumeration.

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
  throw logic_error("Type unrecognised");
}

@ The following code demonstrates how simple recursive descent parsing can be.

The function |scan_in_parens| scans a maximal portion of a type string that is
enclosed in parentheses or brackets. It only deals directly with an arrow that
might separate argument and return type in a function type (but might of
course be absent), and relies on |scan_union_list| to recognise the first and
possibly the second part. The function recognises one or more variants
separated by vertical bars (using a |while| loop but which uses a comma
operator to scan a variant before the loop condition), relying on its turn on
|scan_tuple_list| to deal with individual variants. It must take care though
to unpack the list in case it has length one (no union is involved in this
case). Finally |scan_tuple_list| similarly combines components in a list,
recursively calling |scan_type| for individual components. Unlike unions, we
allow the number of components to be zero, which can be recognised by having a
terminating character at the read position right from the start; when the
happens the loop is no entered at all, an the empty list of components will
give an empty tuple.

@< Function definitions @>=
type_expr scan_in_parens(const char*& s)
{ type_expr a=scan_union_list(s);
  if (*s!='-' or s[1]!='>')
    return a;
  return type_expr(std::move(a),scan_union_list(s+=2));
    // construct function type
}
@)
type_expr scan_union_list(const char*& s)
{ dressed_type_list variants;
  while (variants.emplace_back(scan_tuple_list(s)),*s=='|')
    ++s;
  return variants.size()==1
    ? std::move(variants.front())
    : type_expr(variants.undress(),true);
}
@)
type_expr scan_tuple_list(const char*& s)
{ static const std::string term("|-)]");
  dressed_type_list members;
  if (term.find(*s)==std::string::npos)
    // only act on non-terminating characters
    while (members.emplace_back(scan_type(s)),*s==',')
      ++s;
  return members.size()==1
    ? std::move(members.front())
    : type_expr(members.undress(),false);
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
extern const type_expr row_of_ratvec_type; // \.{[ratvec]}
extern const type_expr row_row_of_int_type; // \.{[[int]]}
extern const type_expr row_row_of_rat_type; // \.{[[rat]]}
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
const type_expr row_of_ratvec_type(mk_type_expr("[ratvec]"));
const type_expr row_row_of_int_type(mk_type_expr("[[int]]"));
const type_expr row_row_of_rat_type(mk_type_expr("[[rat]]"));
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

@* Run-time values.
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

@< Includes needed in \.{axis-types.h} @>=
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
in specific contexts (namely the function templates |uniquify|, |get_own| and
|force_own| defined below) where the actual type will in fact be known, so
some derived types may choose to never use this and not implement |clone| at
all, which is why it is not defined pure virtual. In fact this state of
affairs suggests one could do duplication with a non-virtual (but
systematically named) method instead.

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
#ifndef incompletecpp11
  value_base(const value_base& x) = @[default@];
#else
  value_base(const value_base& x) {}
#endif
public:
  value_base& operator=(const value_base& x) = @[delete@];
};
inline value_base::~value_base() @+{} // necessary but empty implementation
@)
typedef const value_base* value;
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

@*1 Row-of values.
%
Here we define a first type derived from |value_base|, namely the type for
``row of'' types. They are implemented using vectors from the standard
template library.

@< Includes needed in \.{axis-types.h} @>=
#include <vector>
#include <cassert>

@~Since the actual values accessed will be of types derived from |value_base|,
we must pass through a level of indirection, so we have a vector of pointers.
We define these pointers to be |shared_value| pointers, so that the row takes
(shared) ownership of its components without needing an explicit destructor.
This has the additional advantage over explicit ownership management that the
copy constructor, needed for the |clone| method, can safely just
copy-construct the vector of pointers: a possible exception thrown during the
copy is guaranteed to clean up any pointers present in the vector (resetting
their reference counts to their original values). Note also that
default-constructed shared pointers are set to null pointers, so the first
constructor below, which just reserves space for |n| shared pointers, has
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
  template <typename I> row_value(I begin, I end) : val(begin,end)
    @+{} // set from iterator range
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


@*1 Tuple values.
%
Since we use dynamically typed values internally, we can collect the
components of a tuple in a vector without problem. In fact we could reuse the
type |row_value| to hold the components of a tuple, if it weren't for the fact
that it would then print with brackets. Therefore we trivially derive a new
class from |row_value|.

@< Type definitions @>=
struct tuple_value : public row_value
{ tuple_value(size_t n) : row_value(n) @+{}
  template <typename I> tuple_value(I begin, I end) : row_value(begin,end) @+{}
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

@*1 Union values.
%
Formally discriminated unions are in a sense dual to Cartesian products, but
the concrete implementation is fairly different. We must just store an
indication (a tag) of the variant of the union that is actually taken in the
value at hand, and then the component value itself. Since it usually takes
only a few bits to represent the tag, and pointers (especially $64$-bit ones)
have some room to spare, it is tempting to somehow cram the tag bits into a
pointer value and avoid an extra level of dereference. However if we want
orthogonality of the language (the component can be any value, including a
union) there is just no way to elegantly do this, so we store a structure with
a tag and a pointer. We take advantage of padding space this would probably
leave in the structure, and also add the name of the injector function that
was used to produce this value, which will be used for printing purposes only.

@< Type definitions @>=
class union_value : public value_base
{ shared_value comp;
  unsigned short tag;
  id_type injector_name;
public:
  union_value(unsigned short tag,shared_value&& v,id_type name) :
     comp(std::move(v)),tag(tag),injector_name(name) @+{}
  unsigned int variant() const @+{@; return tag; }
  const shared_value& contents() const @+{@; return comp; }
  union_value* clone() const @+{@; return new union_value(*this); }
  void print(std::ostream& out) const;
  static const char* name() @+{@; return "union value"; }
private:
  union_value@[(const union_value& v) = default@];
 // copy constructor; used by |clone|
};
@)
typedef std::unique_ptr<union_value> union_ptr;
typedef std::shared_ptr<const union_value> shared_union;
typedef std::shared_ptr<union_value> own_union;

@ The |print| method for unions will print the value, followed by a dot and
the injector name.

@h "lexer.h" // for |main_hash_table|
@< Function definitions @>=
void union_value::print(std::ostream& out) const
{@; out << *comp << '.' << main_hash_table->name_of(injector_name); }


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
typedef std::unique_ptr<const expression_base> expression_ptr;
typedef std::shared_ptr<const expression_base> shared_expression;

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
determined at \.{atlas} compile time). But it can use a dynamic
cast do determine whether |v| actually is a tuple value or not.

@< Function definitions @>=
void push_expanded(expression_base::level l, const shared_value& v)
{ if (l==expression_base::single_value)
    push_value(v);
  else if (l==expression_base::multi_value)
  { shared_tuple p = std::dynamic_pointer_cast<const tuple_value>(v);
    if (p==nullptr)
      push_value(v);
    else
      for (size_t i=0; i<p->length(); ++i)
        push_value(p->val[i]); // push components
  }
} // if |l==expression_base::no_value| then do nothing

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
that will be pushed onto the stack using the |std::make_shared| template
function. This can either be in the argument expression of |push_value|, in
cases where the object pushed can be constructed in place to its definite
value, or earlier (the result of |std::make_shared| being held by shared
pointer to non~|const|, that is, convertible to |own_value|) if the
constructed object needs modification before being pushed. In both cases the
rvalue reference version of |push_value| will be used (in the latter case by
wrapping the |own_value| in |std::move|), which avoids any manipulation of
reference counts; the constant lvalue reference case is used only in the rare
cases (as in |push_tuple_components| above) where a pre-existing
|shared_value| not held in a local variable is being pushed.

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
|logic_error|. The function template |get| with explicitly provided type
serves for this purpose; it returns a shared pointer (because values on the
stack are shared pointers) of the proper kind.

Sometimes defeating copy-on-write is desired (to allow changes like filling
internal tables that will \emph{benefit} other shareholders), and
|non_const_get| will const-cast the result of |get| to allow that.

@< Template and inline function definitions @>=

template <typename D> // |D| is a type derived from |value_base|
 inline std::shared_ptr<const D> get() throw(logic_error)
{ std::shared_ptr<const D> p=std::dynamic_pointer_cast<const D>(pop_value());
  if (p.get()==nullptr)
    throw logic_error(std::string("Argument is no ")+D::name());
  return p;
}
@.Argument is no ...@>

@)
template <typename D> // |D| is a type derived from |value_base|
  inline std::shared_ptr<D> non_const_get() throw(logic_error)
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
|shared_value| pointer returns a |value| (which is pointer to constant), it
will often be the second one that is selected.

@< Template and inline function definitions @>=
template <typename D> // |D| is a type derived from |value_base|
  D* force (value_base* v)
{ D* p=dynamic_cast<D*>(v);
  if (p==nullptr) throw
    logic_error(std::string("forced value is no ")+D::name());
  return p;
}
@)
template <typename D> // |D| is a type derived from |value_base|
  const D* force (value v)
{ const D* p=dynamic_cast<const D*>(v);
  if (p==nullptr) throw
    logic_error(std::string("forced value is no ")+D::name());
  return p;
}

@ The \.{axis} language allows assignment operations to components of
aggregates (such as rows, matrices, strings, tuples), which in fact assign a
new value to the name bound to the aggregate. Since we implement
copy-on-write, this should in principle make a copy of the aggregate before
modifying the component, but the copy can be avoided if the aggregate value is
unshared. The operation |uniquify| implements making a copy if necessary, and
calling it makes clear our destructive intentions. We make it take a
modifiable lvalue argument into which a new shared (but currently unique)
pointer will be stored in case duplication was necessary, and it returns a raw
pointer-to-non-const version of the new value of that shared pointer, which
can then be used to make the change to the unique copy. The original pointer
is of course still a pointer-to-const, so cannot be used for the modification;
it does however retain ownership. Not surprisingly the implementation of
|uniquify| uses a |const_cast| operation.

Similarly some wrapper functions will want to get unique access to an argument
object, so that they can return a modified version of it as result. This
requires duplicating the argument value on the stack only if it has
|use_count()>1|, which is a rare circumstance since most functions that place
a value on the stack will do so with an unshared pointer. For this we provide
a variant function template of |get| called |get_own|. It has the same
prototype as |non_const_get|, but like |uniquify| respects copy-on-write by
making a copy first in case there are other shareholders. The implementation
somewhat differs from |uniquify| in that it uses a |std::const_pointer_cast|,
and only in the case no copy is made (in which case that cast is basically all
it does); in the case a copy is made, it returns a shared pointer to the copy
(modifiable), while leaving the original shared pointer unchanged (except that
its reference count will be decremented).

Since these functions return pointers that are guaranteed to be unique, one
might wonder why no use of |std::unique_ptr| is made. The answer is this is
not possible, since there is no way to persuade a |shared_ptr| to release its
ownership (as in the |release| method of unique pointers), even if it happens
to be (or is known to be) the unique owner.

@< Template and inline function def... @>=
inline value_base* uniquify(shared_value& v)
{ if (not v.unique())
     v=shared_value(v->clone());
  return const_cast<value_base*>(v.get());
}
@)
template <typename D> // |D| is a type derived from |value_base|
  std::shared_ptr<D> get_own()
{ std::shared_ptr<const D> p=get<D>();
  if (p.unique())
    return std::const_pointer_cast<D>(p);
  return std::shared_ptr<D>(p->clone());
}

@ Finally there is |force_own| that intends to be to |force| what |get_own| is
to |get|, both requiring its argument to currently hold a pointer to a
specific type derived from |value_base|, and returning a shared pointer to
non-|const| version of it, with duplication being applied only if the original
pointer was not |unique|. In order for this to work, the original copy of the
pointer must be cleared at the time |unique| is called, so that the call has
some chance of returning |true|; this is obtained by taking the argument as
rvalue reference, and making sure a temporary is move-constructed from it
inside the body of |force_own|.

@< Template and inline function def... @>=
template <typename D> // |D| is a type derived from |value_base|
  std::shared_ptr<D> force_own(shared_value&& q)
{ std::shared_ptr<const D> p=
     std::dynamic_pointer_cast<const D>(shared_value(std::move(q)));
  if (p==nullptr) throw
    logic_error(std::string("forced value is no ")+D::name());
  if (p.unique())
    return std::const_pointer_cast<D>(p);
  return std::shared_ptr<D>(p->clone());
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
expression~|e| is thrown.

The function |row_coercion| specialises if possible |compon_variant| in such a
way that the corresponding row type can be coerced to |final_type|, and
returns a pointer to the |conversion_record| for the coercion in question. The
function |coercion| serves for filling the coercion table.

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
conversions. Which conversion is to be applied is determined by
|conversion_type|, a reference to a |conversion_info| structure containing a
function pointer |convert| that provides the actual conversion routine; this
is easier than deriving a plethora of classes differing only in their virtual
method |evaluate|, although it is slightly less efficient. The type |conv_f|
of conversion function pointer specifies no argument or return value, as we
know beforehand that the |shared_value| objects serving for this are to be
found and left on the runtime stack.

@< Type definitions @>=
struct conversion_info
{ typedef void (*conv_f)();
  conv_f convert;
  std::string name;
  conversion_info(const char* s,conv_f c) : convert(c),name(s) @+{}
};
@)
class conversion : public expression_base
{ const conversion_info& conversion_type;
  expression_ptr exp;
public:
  conversion(const conversion_info& t,expression_ptr e)
@/: conversion_type(t),exp(std::move(e)) @+{}
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
{ exp->eval();
  (*conversion_type.convert)();
  if (l==no_value)
    execution_stack.pop_back();
}
@)
void conversion::print(std::ostream& out) const
@+{@; out << conversion_type.name << ':' << *exp; }

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
operate on any input type: the ``voiding'' coercion. It can be applied during
type analysis to for instance to allow a conditional expression in a void
context to have branches that evaluate to different types (including the
possibility of absent branches which will be taken to deliver an empty tuple):
those branches which do not already have void type will get their resulting
value voided, so that in the end all branches share the void type of the
entire conditional.

At runtime, voiding is mostly taken care of by the |level| argument~|l| passed
around in the evaluation mechanism. Any subexpressions with imposed void type,
like the expression before the semicolon in a sequence expression, will get
their |evaluate| method called with |l==no_value|, which suppresses the
production of any value on the |execution_stack|. This setting will be
inherited down to for instance the branches, if the subexpression was a
conditional, so nothing special needs to be done to ensure that branches which
originally had non-void type get voided. However, there are some rare cases
where the void type does not derive from the syntactic nature of the context,
such as in the right hand side of an assignment to a variable that happens to
have |void| type (a quite useless but valid operation). In these cases the
type analysis will have to explicitly insert a mechanism to avoid a value to
be produced (and in the example, assigned) where none was intended.

The |voiding| expression type serves that purpose. It distinguishes itself
from instances of |conversion|, in that |voiding::evaluate| calls |evaluate|
for the contained expression with |l==no_value|, rather than with
|l==single_value| as |conversion::evaluate| does. If fact that is about all
that |voiding| is about.

@< Type definitions @>=
class voiding : public expression_base
{ expression_ptr exp;
public:
  voiding(expression_ptr e) : exp(e.release()) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ The |voiding::evaluate| method should not ignore its own |level| argument
completely: when called with |l==single_value|, an actual empty tuple should
be produced on the stack, which |wrap_tuple<0>()| does.

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
appropriate entry, and wraps |e| into a corresponding |conversion| if it finds
one. Ownership of the expression pointed to by |e| is handled implicitly: it
is released during the construction of the |conversion|, and immediately
afterwards |reset| reclaims ownership of the pointer to that |conversion|.
Note that the |conversion| constructor uses |*it| only for its
|conversion_info| base type; the types |from| and |to| are not retained at run
time.

When |to_type==void_type|, the conversion always succeeds, as the syntactic
voiding coercion is allowed in all places where |coerce| is called. We used to
do |e.reset(new voiding(std::move(e)));| in that case as well, which is a safe
way to ensure that voiding will take place dynamically whenever present
syntactically. However as argued above, in most cases no runtime action is
necessary, since the context will have already set |l==no_value| for
subexpressions with void type. In removing that statement from the code below,
we have moved responsibility to type analysis, which must now ensure that any
subexpression with void type will have |l==no_value| when evaluated. This
means that in some cases a |voiding| has to be constructed explicitly. This
crops up in many places (for instance a component in a tuple display just
might happen to have void type), but almost all such cases are far-fetched; so
we trade some economy of code for efficiency in execution.

@< Function definitions @>=
bool coerce(const type_expr& from_type, const type_expr& to_type,
	    expression_ptr& e)
{ if (to_type==void_type)
  {@;
     return true;
  } // syntactically voided here, |e| is unchanged
  for (auto it=coerce_table.begin(); it!=coerce_table.end(); ++it)
    if (from_type==*it->from and to_type==*it->to)
    @/{@; e.reset(new conversion(*it,std::move(e)));
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
  return std::move(d); // invoking |std::move| is necessary here
}

@ List displays and loops produce a row of values of arbitrary (but identical)
type; when they occur in a context requiring a non-row type, we may be able to
find a coercion that reconciles the requirements. The following function finds
whether this is possible, and if so sets |component_type| to the type required
for the components; if not it will return~|nullptr|. The implementation is
simply to look in |coerce_table| for conversions from some row type, then try
to specialise |component_type| to the component type of that row type. By
using this mechanism, the components of a row display themselves obtain a
context that may generate further conversions to obtain this type.

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
    if (final_type==*it->to and it->from->kind()==row_type)
      return component_type.specialise(*it->from->component_type())
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
\foreign{a priori} type \.{([vec],[int])} that can be converted either to
type \.{([[int]],[rat])} or to \.{(mat,ratvec)} if the context so requires .
So our relation will be a partial order, and compatible with tuple and row
formation: $x_i\leq y_i$ for all~$i$ implies
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

For types $t_1$ and~$t_2$ which do not admit a relative priority of one over
the other to be established, we shall disallow simultaneous overloading with
arguments types $t_1$ and~$t_2$ if any expression given as argument could be
converted to either of them. Deciding the existence of such an expression
would require study of all available language constructs, but the situation is
somewhat simplified by the fact that, for efficiency reasons, overloading
resolution is not done using the complete argument expression, but only using
its type. In fact matching will be done using calls to the very function
|is_close| we are discussing here, testing the bit for conversion towards the
required argument type; this provides us with an opportunity to adjust rules
for possible type conversions of arguments at the same time as defining the
exclusion rules.

Empty row displays, or more precisely arguments of type~\.{[*]}, pose a
difficulty: they would be valid in any context requiring a specific row type,
so if we stipulated that one may write \.{[]} to designate an empty row
operand of any row type, then |is_close| would have to consider all row types
close to each other (and therefore mutually exclusive for overloading). This
used to be the convention adopted, but it was found to be rather restrictive
in use, so the rules were changed to state that an un-cast expression \.{[]}
will not match overload instances of specific row types; it might match an
parameter of specified type \.{[*]} (we allow that as type specification for
function parameters), and such a parameter can \emph{only} take an empty
list corresponding argument (this is very limiting of course, but it allows
being explicit about which overloaded instance should be selected by an
argument \.{[]}, by defining an overload for \.{[*]} that explicitly calls
the one for say \.{[vec]}).

These considerations are not limited to empty lists (although it is the most
common case): whenever an expression has an \foreign{a priori} type
containing \.*, that expression will not select any overload with a concrete
type in its place (overloading does not perform type specialisation).


@ So here is the (recursive) definition of the relation |is_close|. Equal
types are always close, while undetermined types are not convertible to any
other type. A primitive type is in the relation |is_close| to another type
only if it is identical or if there is a direct conversion between the types,
as decided by~|coerce|. Two row types are in the relation |is_close| if their
component types are, and two tuple types are so if they have the same number
of component types, and if each pair of corresponding component types is
(recursively) in the relation |is_close|.

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
  if (x.kind()==undetermined_type or y.kind()==undetermined_type)
    return 0x0;
      // undetermined types do not specialise (or coerce), and are not close
  if (x==void_type or y==void_type)
    return 0x0; // |void| does not allow coercion for overload, and is not close
  if (x.kind()==primitive_type or y.kind()==primitive_type)
  { unsigned int flags=0x0;
    if (coerce(x,y,dummy)) flags |= 0x1;
    if (coerce(y,x,dummy)) flags |= 0x2;
    return flags==0 ? flags : flags|0x4;
  }
  if (x.kind()!=y.kind())
    return 0x0;
  if (x.kind()==row_type)
    return is_close(*x.component_type(),*y.component_type());
  if (x.kind()!=tuple_type)
    return 0x0; // non-aggregate types are only close if equal
  unsigned int flags=0x7; // now we have two tuple types; compare components
  wtl_const_iterator it0(x.tuple());
  wtl_const_iterator it1(y.tuple());
  while (not it0.at_end() and not it1.at_end() @|
         and (flags&=is_close(*it0,*it1))!=0)
  @/{@; ++it0; ++it1; }
  return it0.at_end() and it1.at_end() ? flags : 0x0;
}

@ For balancing we need a related but slightly different partial ordering on
types. The function |broader_eq| tells whether an expression of \foreign{a
priori} type |b| might also valid (with possible coercions inserted) in the
context of type~|a|; if so we call type |a| broader than~|b|.

@< Declarations of exported functions @>=
bool broader_eq (const type_expr& a, const type_expr& b);

@~The relation |greater_eq(a,b)| does not imply that values of type~|b| can be
converted to type~|a|; what is indicated is a possible conversion of
expressions might depend on the form of the expressions, and only be achieved
by coercions applied to more or less deeply nested subexpressions (for
instance a list display is required if $a$ is \.{[rat]} and |b| is \.{[int]}).
It is not the responsibility of |greater_eq| to decide whether this succeeds
or not (an error will be thrown later if it turns out to fail), but it defines
a partial ordering to guide which types are candidate types will be
considered.

Since these are types deduced for expressions rather than required type
patterns, \.* stands for tho most narrow rather than a very broad type, indeed
one for an expression like \.{die} that can never return a value, and the
type \.{[*]} is narrower than any other row type (the only value that an
expression with this type can return is an empty list). Nonetheless
``broader'' does not imply more possible values, and the type \.{void} is the
broadest of all since all values can be converted to it.

The implementation is by structural recursion, like |is_close|, but some
details are different. We do take |void_type| into consideration here, as well
as the (rare) type |unknown_type|. Since we want to define a partial ordering,
we must forbid one direction of all two-way coercions; since those always
involve exactly one primitive types, we do this by only allowing the
conversion in the direction of the primitive type. For the rest we just do the
recursion in the usual way with just one twist: function types can only be
comparable if they have equal argument types, but there might be a recursive
|broader_eq| relation between the result types (because a coercion might
``creep into'' the body of a lambda expression; this is not likely, but we can
cater for it here).

@< Function definitions @>=
bool broader_eq (const type_expr& a, const type_expr& b)
{ if (a==void_type or b==unknown_type)
    return true;
  if (a==unknown_type or b==void_type)
    return false;
  if (a.kind()==primitive_type)
    return (is_close(a,b)&0x2)!=0; // whether |b| can be converted to |a|
  if (a.kind()!=b.kind())
    // includes remaining cases where |b.kind()==primitive_type|
    return false; // no broader between different kinds on non-primitive types
  if (a.kind()==row_type)
    return broader_eq(*a.component_type(),*b.component_type());
  if (a.kind()==function_type)
    return a.func()->arg_type==b.func()->arg_type and @|
    broader_eq(a.func()->result_type,b.func()->result_type);
  wtl_const_iterator itb(b.tuple());
  for (wtl_const_iterator ita(a.tuple());
       not ita.at_end(); ++ita,++itb)
    if (itb.at_end() or not broader_eq(*ita,*itb)) return false;
  return itb.at_end(); // if list of |a| has ended, that of |b| must as well
}

@* Error values.
%
Before we describe evaluation of expressions we must realise that evaluation
can cause runtime errors.

@< Includes needed in \.{axis-types.h} @>=
#include <stdexcept>

@ We use classes derived from |std::exception| and similar standard ones like
|std::runtime_error|, but we define our own local hierarchy, with
|atlas::interpreter::error_base| as base class. The main reason to do this is
to have a centralised error message to which exception handlers have write
access, so that it is possible to extend the error message and then re-throw
the same error object. The simplest way to allow this is to give public access
to that string member, so we make this a |struct| rather than a |class|.

@< Type definitions @>=
struct error_base : public std::exception
{ std::string message;
  explicit error_base(const std::string& s) : message(s) @+{}
#ifdef incompletecpp11
  ~error_base () throw() @+{} // backward compatibility for gcc 4.4
#endif
  const char* what() const throw() @+{@; return message.c_str(); }
};

@ We classify errors into three classes: those due to inconsistency of our
(rather than the user's) program are classified |logic_error|, those arising
during the analysis of the user program are are classified |program_error|,
and those not caught by analysis but during evaluation are classified
|runtime_error|. The  first and last are similar to exceptions of the same
name in the |std| namespace, but they are not derived from those exception
classes.

@< Type definitions @>=
struct logic_error : public error_base
{ explicit logic_error(const std::string& s) : error_base(s) @+{}
#ifdef incompletecpp11
  ~logic_error () throw() @+{} // backward compatibility for gcc 4.4
#endif
};
@)
struct program_error : public error_base
{ explicit program_error(const std::string& s) : error_base(s) @+{}
#ifdef incompletecpp11
  ~program_error () throw() @+{} // backward compatibility for gcc 4.4
#endif
};
@)
struct runtime_error : public error_base
{ explicit runtime_error(const std::string& s) : error_base(s) @+{}
#ifdef incompletecpp11
  ~runtime_error () throw() @+{} // backward compatibility for gcc 4.4
#endif
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
#ifdef incompletecpp11
  ~expr_error() throw() @+{}
#endif
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
#ifdef incompletecpp11
  type_error(type_error&& e)
  : expr_error(std::move(e))
  , actual(std::move(e.actual)), required(std::move(e.required)) @+{}
  ~type_error () throw() @+{}
#else
  type_error@[(type_error&& e) = default@];
#endif
};

@ For type balancing, we shall use controlled throwing and catching of errors
inside the type checker, which is facilitated by using a dedicated error type.
If balancing ultimately fails, this error will be thrown uncaught by the
balancing code, so |catch| blocks around type checking functions must be
prepared to report the types that are stored in the error value.

@< Type definitions @>=
struct balance_error : public expr_error
{ containers::sl_list<type_expr> variants;
  balance_error(const expr& e)
  : expr_error(e,"No common type found"), variants() @+{}
#ifdef incompletecpp11
  balance_error(const balance_error&)=@[delete@];
  balance_error(balance_error&& o) @|
  : expr_error(std::move(o)),variants(std::move(o.variants)) @+{}
  ~balance_error() throw() @+{}
#endif
};

@ Here is another special purpose error type, throwing of which does not
always signal an error. In fact it is meant to always be caught silently,
since its purpose is to implement a |break| statement that will allow
terminating the execution of loops other than through the regular termination
condition. During analysis we check that any use of |break| occurs within some
loop, and without being inside an intermediate $\lambda$-expression (which
would allow smuggling it out of the loop), so not catching a thrown
|loop_break| should never happen. Therefore we derive this type from
|logic_error|, so that should this ever happen, it will be reported as such.

@< Type definitions @>=
struct loop_break : private logic_error
{ unsigned depth;
  loop_break(int n)
  : logic_error("Uncaught break from loop"), depth(n) @+{}
#ifdef incompletecpp11
  ~loop_break () throw() @+{}
#endif
} ;


@ Similarly to the above, there is a |function_return| object that can be
thrown to implement a |return| expression. Again analysis will have ensured
that the object will always be caught, so this type as well will be derived
from |logic_error|.

@< Type definitions @>=
struct function_return : private logic_error
{ shared_value val;
  function_return(shared_value&& val)
  : logic_error("Uncaught return from function"), val(val) @+{}
#ifdef incompletecpp11
  ~function_return () throw() @+{}
#endif
} ;


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
