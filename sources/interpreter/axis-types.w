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
\def\axis.{\.{axis}}

@* Outline.
%
This file describes a central part of the interpreter for Axis, the (new)
command language of the Atlas of Lie Groups and Representation software. This
part describes the fundamental type declarations used throughout the
interpreter, as well as some functions that operate on types, for instance to
see if one can be converted into another. This compilation unit is
called \.{axis-types}, which refers both to user types and to meta-types of the
interpreter used to represent user types and other fundamental notions.

@( axis-types.h @>=

#ifndef AXIS_TYPES_H
#define AXIS_TYPES_H

#include "axis-types-fwd.h"

@< Includes needed in \.{axis-types.h} @>@;
namespace atlas { namespace interpreter {
@< Type definitions @>@;
@< Declarations of global variables @>@;
@< Declarations of exported functions @>@;
@< Template and inline function definitions @>@;
}}
#endif

@ Some header files need to know about classes and class templates defined here,
but do not want to include \.{axis-types.h} to avoid mutual dependence problems.
For those uses we write a separate header file \.{axis-types-fwd.h} which
introduces the relevant names as declared but not defined (incomplete types).

@( axis-types-fwd.h @>=
#ifndef AXIS_TYPES_FWD_H
#define AXIS_TYPES_FWD_H

@< Includes needed in \.{axis-types-fwd.h} @>@;
namespace atlas { namespace interpreter {
@< Type declarations @>@;
}}
#endif

@ Each compilation unit follows a similar global pattern. While some modules are
absent in each case, the order in the implementation files is local type
definitions, global variables, local variables, local functions, global
functions (the last point being the main goal of the implementation unit). While
the local variable and function definitions are in an anonymous namespace to
prevent the linker from exporting their names, the local type definitions must
not be so enclosed, because they may concern types declared as members of a
global class, to be made complete by a local definition.

@h "axis-types.h"
@h <cstdlib>

@c

namespace atlas { namespace interpreter {
@< Global variable definitions @>@;
@< Local type definitions @>@;
namespace {@;
  @< Local function definitions @>@;
}
@< Function definitions @>@;
}}


@ The parser produces a parse tree, in the form of a value of type |expr|
defined in the unit \.{parsetree}. The task of the evaluator (defined largely in
the \.{axis.w} compilation unit) is to take such an expression, analyse it and
then evaluate it to obtain a value. At various points errors can occur, which we
shall have to handle gracefully: during the analysis static ``type'' errors may
prevent us from undertaking any meaningful action, and in absence of such errors
there still can be dynamic ``runtime'' errors.

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
happens in the container classes of the Standard Template Library. For this
reason the code using these types is to some extent exposed to the linked
nature of the data, and some concern for proper memory management will be
necessary. We shall however provide functions to facilitate safe and
consistent handling of storage of these data, while avoiding excessive
duplication during manipulation. Altogether this part of the program will be
of a quite different nature than that of the main mathematical library.


@< Includes needed in \.{axis-types-fwd.h} @>=
#include <memory> // for |std::unique_ptr|, |std::shared_ptr|
#include "Atlas.h"
  // for utilities (like |sl_list|); this include must come first
#include "buffer.h" // for |id_type|


@* Types to represent types.
%
We make a fundamental choice to check types before attempting to execute any
expression entered by the user; thus type errors can be signalled earlier and
expressed more understandably than if one would wait until an accident is about
to happen at runtime. Therefore types will be represented independently of any
runtime value representation. Types will have an effect on the transformation of
expressions into executable internal form (for example in selecting the instance
of an overloaded function to use), but once that is done, types will no longer
be represented at all: the internal code executes without being ``aware'' of the
types of the values it manipulates. This is only possible because those values
are then accessed via generic pointers (i.e., pointers that, as far as \Cpp\ is
concerned, could point to any kind of value). Run time values have different
representations according to their top level structure (for instance an integer
value and a row value have different structures), which are realised by using
classes derived from a common base class for run time values; by using dynamic
casts the program can check whether the dynamic value of expressions match their
static type, and trap errors in the logic of our interpreter. The user should
never notice these tests, and in the optimised version of the \.{atlas} program,
no run time checking is done at all (dynamic casts are replaced by static ones).

In the simplest case, types are represented by small tree structures, with nodes
of types |type_expr| that will be detailed later. Unlike \Cpp~classes, we never
have types that gives some visibility of the structure of values but which limit
access to handling that structure to privileged methods: whenever we do
encapsulate representation structure, like for primitive types, the type of
expressions for such values is completely opaque (only functions coming with the
type can be applied to them). So a type expression basically describes the
accessible structure of the corresponding values, and allows all operations
compatible with that structure to be performed. For instance the type printed
as \.{(int->[(bool,mat)])} specifies a function mapping integers to lists of
pairs of a Boolean value and a matrix; it allows operations compatible with
that, for instance calling it with an integer argument, then subscripting it
with an integer index, and finally selecting the second component of the result
to obtain a matrix. As mentioned, classes defined in the Atlas software library
itself will be presented to the user as (opaque) primitive types; only built-in
operations declared for that primitive type can be applied to such values. In a
development beyond the original type structure, some expressions can have an
abstract type meaning the corresponding value can be of any kind whatsoever, but
that type will be represented by a type variable that is treated as a distinct
primitive type, so that no usage of the expression that makes any assumption
about the value structure will be allowed. More generally types can be ``second
order types'', namely type expressions that can contain type variables. Then
there are the main kinds of non-primitive non-abstract types: function types,
row-, tuple-, and discriminated union types.

Trees representing different type expressions will not share any sub-trees
among each other, so they have strict (unshared) ownership of their parts,
which simplifies memory management. This choice means some recursive copying
of tree structures is sometimes required, but (although runtime cost of type
handling is negligible with respect to other factors in the interpreter) we
avoid doing so more than absolutely necessary.

Usually types are built from the bottom up (from leaves to the root), although
during type checking the reverse also occurs; also type variables can be
substituted for in certain circumstances which also leads to top-down growth of
type expressions. During bottom-up construction, new nodes are held in local
variables before being moved to dynamically managed storage when becoming
subexpressions of newer nodes. This means move semantics is needed for type
expressions (transferring ownership of descendants when moving a node), which
was originally realised using auto-pointers to descendants. With the advent
of \Cpp11, special syntactic support for move semantics was added, so that one
can distinguish calls that should be realised with shallow copies and ownership
transfer, from the occasionally needed deep copy. At the same time auto-pointers
were replaced by unique-pointers, with essentially the same functionality, but
better adapted to the syntactic facilities. The structure also uses ordinary
pointers of type~|type_p| in places where ownership is managed by some
containing structure rather than by the pointer itself.

@< Type declarations @>=
class type_expr;
@/
using type_p = type_expr*;
using const_type_p = const type_expr*;
using type_ptr = std::unique_ptr<type_expr>;
using const_type_ptr = std::unique_ptr<const type_expr>;

@*2 Type lists.
%
We also need type lists as building block for types (for instance for the
arguments of a function). These used to be defined in a similar manner to
types themselves, but since an STL-compatible singly-linked list container
class templates |simple_list| and |sl_list| was added to the Atlas utilities
library, these were used to replace the implementation of type lists.

When incorporating type lists into the |type_expr| structure, the class template
|simple_list| will be used; the more flexible but less compact |sl_list|
template will be occasionally used for temporary variables, whose type is then
|dressed_type_list|. This usage provides a good practical test for the usability
of these new container types; so far they have passed them gracefully, though
occasionally after enriching the repertoire of methods of the container type.
Also, it turns out to be useful to sometimes not use the provided container
types directly; for instance, for a component of a \Cpp\ |union| type, there is
little advantage of using a smart pointer (it still needs to be managed manually
by the containing union type), so we use a raw pointer to a node
(|raw_type_list| or |const_raw_type_list|) in such occasions. We shall also have
occasions where instead of normal iterators over type lists we use weak
iterators (ones that cannot be used to insert or delete nodes), and the types
|wtl_iterator| and |wtl_const_iterator| are defined for that purpose.

@< Type declarations @>=
using type_list = containers::simple_list<type_expr>;
using dressed_type_list = containers::sl_list<type_expr>;
@)
using raw_type_list = atlas::containers::sl_node<type_expr>*;
@/using const_raw_type_list = atlas::containers::sl_node<type_expr>const *;
using wtl_iterator = containers::weak_sl_list_iterator<type_expr>;
  // wtl = weak type list
using wtl_const_iterator = containers::weak_sl_list_const_iterator<type_expr>;

@ Since types and type lists own their trees, copying them means making a deep
copy. The class |type_expr| will provide no copy constructor but instead a more
visible |copy| method to do a deep copy; for |type_expr| we cannot add any
methods and copying must be done by hand, but it is hardly ever necessary. On
some occasions we need a copy of |type_expr| at the pointer level: we have
access to a pointer to some |type_expr| (possibly by simply taking its address),
but need an owning |type_ptr| instead. The function |acquire| will achieve this.

@< Declarations of exported functions @>=
type_ptr acquire(const_type_p t);

@~The function |acquire| simply creates a unique-pointer from the call of |new|
(using the |std::make_unique| function template), whose constructing expressions
invokes the |copy| method of |type_expr| pointed to. That will recursively
perform a deep copy, as detailed in section @#type expression copy@>.

@< Function definitions @>=
type_ptr acquire(const_type_p t) @+
{@; return std::make_unique<type_expr>(t->copy()); }


@ Type lists are usually built by starting with a default-constructed
|type_list| (or equivalently initialising it with |empty_tuple()| defined below,
which can avoid a Most Vexing Parse situation if |type_list()| were used
instead), and repeatedly calling |prefix| to add nodes in front. The function
|prefix| is efficient both in handling the node and the list itself, due to its
use of move-semantics. Nonetheless, the list is passed as a modifiable lvalue
reference |dst| that will be made to hold the extended list, which is also the
result of |prefix| call, for convenience. A dressed variation of |prefix| is
also defined below.

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

@*2 Kinds of type; primitive types.
%
We have a simple but flexible type model. There is a finite number of
``primitive'' types, many of which are abstractions of complicated classes
defined in the Atlas library, such as root data or reductive groups. We also
have type variables, which as mentioned are treated similarly to primitive
types. Then one has types for functions that map one type to another, types that
are ``row of'' some other type, tuples (Cartesian products) of some given
sequence of types, and disjoint unions (co-products in the category of sets) of
some sequence of types.

In addition to these, we allow for two more possibilities, that do not
correspond to the way in which values can be built up. The final possibility
|tabled| provides a way to reference a type indirectly (its details are to be
found after looking up the reference), which is essential for being able to
specify recursive types. We also allow for an undetermined type, which can
serve as a wild-card when specifying type patterns rather than complete types;
this possibility will be used in the default constructor so as to not commit to
any of the concrete variants.

Before we can define |type_expr|, we need to enumerate its variants and the
possibilities for primitive types. Here are enumerations of tags for the basic
kinds of types, and for the sub-cases of |primitive_type| (some of them will be
introduced later).

@< Type definitions @>=
enum type_tag
 { undetermined_type, primitive_type, variable_type,
   function_type, row_type, tuple_type, union_type,
   tabled
 };

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
Type expressions are defined by a tagged union, where the tags indicate the kind
of type denoted (primitive, function, row, etc.). We always only access the
variant corresponding to the current tag value, but this is not something that
can be statically ascertained in \Cpp. Nonetheless the fields are private, with
public accessor methods.

That variant |tabled|, used for applied user-defined type (constructor)
identifiers, makes use of a structure |tabled_type_cons| (which is larger than
other variants, thus enlarging |type_expr|, but this is no big deal). Its |nr|
field changes interpretation between parsing and type analysis: the parser first
stores the type identifier there, but this will be looked up in the
|global_id_table|, after which either the |type_expr| is replaced by the
definition stored there, or the type remains |tabled| but with |nr| replaced by
an index into the |type_map|, where most ``heavy weight'' user type definitions
are recorded, notably those of recursive types and type constructors, that can
never be full expanded. Functions like the equality operation on |type_expr|
values take care to not indefinitely expand references into |type_map|; they
avoid specifically to do this when both types (or constructors) are recursive,
in which case they just test equivalence of any type arguments present.

Although \Cpp11 allows variant members of a |union| with nontrivial special
member functions, so that they could be smart pointer types, it leaves it to the
programmer's responsibility to explicitly call constructors and destructors as
those variants come and go; this effectively ruins most of the advantages that
smart pointers would give. For this reason we use raw pointers here instead, and
in particular |raw_type_list| rather than |type_list| for the |tuple_variant|
and the |type_args| field in |tabled_type_cons|. One drawback of that is that we
will not be able to create a |type_list::iterator| for traversal of the list,
but in practice weak iterators (specifically |wtl_iterator|), which one can
construct from a raw pointer, will always suffice.

There is one restriction on types that is not visible in the definition below,
namely that the list of types referred to by the |tuple_variant| field cannot
have length~$1$. This is because syntactically anything that would suggest a
$1$-tuple or $1$-union (for instance a parenthesised expression) in fact
designates what would be the unique component of such a tuple or union. We do
allow a list of length~$0$ when |tag==tuple_type|, which defines the |void| type
as we saw above. We do not use an empty list with |tag==union_type|, even though
a $0$-union would naturally represent an uninhabited type. We do have use for
such a type, to assign to expressions whose evaluation never returns normally at
all, like calls of |error|; however we have chosen to instead use the impossible
polymorphic type (for all \.{A}) \.{A} (that is, ``all types at once'') for that
purpose.

@< Type definitions @>=

using type_nr_type = id_type; // interpretation flips after parsing
struct tabled_type_cons {@; type_nr_type nr; raw_type_list type_args; };
class type_expr
{ type_tag tag;
  union
  { primitive_tag prim_variant; // when |tag==primitive|
    unsigned int typevar_variant; // when |tag==variable_type|
    func_type* func_variant; // when |kind==function_type|
    type_p row_variant; // when |kind==row_type|
    raw_type_list tuple_variant; // when |kind==tuple_type| or |kind==union_type|
    tabled_type_cons tabled_variant; // when |kind==tabled|
  };
@)
  class defined_type_mapping;
  static defined_type_mapping type_map;
@)
public:
  @< Ordinary methods of the |type_expr| class @>@;
private:
  @< Private methods of the |type_expr| class @>@;
@)
public:
  @< Static methods of |type_expr| that will access |type_map| @>@;
};

@ We start with a number of public methods that allow testing the actual variant
of a |type_expr|, and after that the data specific to that variant. The method
|raw_kind| should in general be used first, in order to determine which variant
of the union is active. Then according to the value found, one can call one of
the methods |prim|, |typevar_count|, |func|, |component_type|, |tuple| (the last
one serving both for |tuple_type| and |union_type|), or one of |tabled_nr|,
|tabled_args|, and some others mentioned below, when |raw_kind()==tabled|. In
the final case one often is interested in what the definition of the tabled type
expands to, and this information is given by the |top_level|, which returns a
reference to the relevant |type_map| entry. For most purposes, this is not a
practical value to use however, since any arguments of a |tabled| type
constructor are represented by type variables, and the appropriate elements of
|tabled_args()| still need to be substituted for them; the exception to this is
when one just needs to know the |tag| in the type (constructor) definition, and
the method |top_kind| can be used to get this tag.

Before user defined type constructors were introduced it used to be the case
that |top_kind| (then simply called |kind|) was used where we now use
|raw_kind|, and methods like |func|, |component_type| or |tuple| would similarly
call |top_level| implicitly whenever |tag==tabled|. In cases where we want to
treat the tabled case transparently by expanding the definition, we now instead
use the |expanded| method. This method actually performs the substitution of
type constructor arguments, and therefore returns |type_expr| by value rather
than by reference; it will usually require introducing a local variable to hold
the result. Sometimes it is convenient to actually replace a |type_expr| by this
expansion when |tag==tabled|, and the non-|const| method |expand| provides for
this.

The table |type_map| provides some additional information about its types, such
as names of their fields in case of tuple or union types, and whether the
definition is actually recursive. The methods |tabled_arity|, |is_recursive| and
|type_name| access such information for a |type_expr| with |tag==tabled|.

@< Ordinary methods of the |type_expr| class @>=
type_tag raw_kind () const @+{@; return tag; } // don't translate |tabled|
const type_expr& top_level () const; // what |tabled_variant| is equated to
type_tag top_kind () const @+
{@; return raw_kind()==tabled ? top_level().raw_kind() : raw_kind(); }
@)
primitive_tag prim () const @+
    {@; assert(tag==primitive_type); return prim_variant; }
unsigned int typevar_count () const @+
    {@; assert(tag==variable_type); return typevar_variant; }
const func_type* func() const  @+
    {@; assert(tag==function_type); return func_variant; }
  func_type* func() @+
    {@; assert(tag==function_type); return func_variant; }
const type_expr& component_type () const @+
    {@; assert(tag==row_type); return *row_variant; }
 type_expr& component_type () @+
    {@; assert(tag==row_type); return *row_variant; }
const_raw_type_list tuple () const @+
    {@; assert(tag==tuple_type or tag==union_type); return tuple_variant; }
  raw_type_list tuple () @+
    {@; assert(tag==tuple_type or tag==union_type); return tuple_variant; }
@)
type_nr_type tabled_nr () const @+
    {@; assert(tag==tabled); return tabled_variant.nr; }
const raw_type_list tabled_args() const @+
    {@; assert(tag==tabled); return tabled_variant.type_args; }
unsigned short tabled_arity() const; // number of arguments taken
bool is_recursive() const; // whether |tabled| type is recursive
id_type type_name () const; // identifier corresponding to |tabled_variant|
type_expr expanded () const; // top level expansion of |tabled_variant|
type_expr& expand ()
{@; if (tag==tabled)
      *this = expanded();
    return *this;
}

@ There is a default constructor that sets |kind==undetermined_type| and
constructs no variant at all (in \Cpp\ \emph{at most} one variant of a union is
active). For every variant of the union there is a factory method, which first
default-constructs a |type_expr result@;|, and then changes its variant while
assigning to the corresponding field. This is possible because all variants of
the |union| are POD types (such as raw pointers): no assignment will involve
(incorrectly) destructing the ``previous value'' of the field.

The various factory methods all take rvalue reference arguments to the
information to be incorporated, which can be moved into the appropriate field
after setting |tag| appropriately. For the tuple, union, and (tabled) used
defined types, the pilfering is done by calling the |simple_list::release|
method which cracks open the list object and returns its contents as raw
pointer, while leaving a null pointer behind. For function types there are two
rvalue references to |type_expr| fields, for argument and result types, whose
contents will be move assigned inside the factory function. It will prove
convenient to handle tuple and union constructors together in a |tuple_or_union|
factor function, to which we provide the actual tag wanted as an argument. There
are two factory functions for tabled types: |user_type| to explicitly provide
the parts to assemble, and |local_ref| that constructs a default argument list
(as required when mutually recursive type constructors refer to each other) on
the fly, given just the arity |n_args|.

@< Ordinary methods of the |type_expr| class @>=

type_expr() noexcept : tag(undetermined_type) @+{}
static type_expr primitive(primitive_tag p)
@/{@; type_expr result;
  result.tag=primitive_type,
    result.prim_variant=p;
  return result;
}
static type_expr variable(unsigned int nr)
@/{@; type_expr result;
  result.tag=variable_type,
    result.typevar_variant=nr;
  return result;
}
static type_expr function(type_expr&& arg, type_expr&& result);
static type_expr row(type_expr&& t)
@/{ type_expr result;
    result.tag=row_type,
      result.row_variant=new type_expr(std::move(t));
    return result;
  }
static type_expr tuple(type_list&& l)
@/{@; type_expr result;
    result.tag=tuple_type,
      result.tuple_variant=l.release();
    return result;
  }
static type_expr tuple_or_union(type_tag tag,type_list&& l)
@/{@; type_expr result;
  result.tag=tag,
    result.tuple_variant=l.release();
  return result;
}

static type_expr user_type(type_nr_type type_nr,type_list&& l)
{ type_expr result;
  result.tag=tabled,
    result.tabled_variant.nr=type_nr,
    result.tabled_variant.type_args=l.release();
  return result;
}

static type_expr local_ref(type_nr_type type_nr, unsigned int n_args)
  // local reference within definition group
{ type_list args;
  for (unsigned int i=n_args; i-->0; )
    args.push_front(type_expr::variable(i));
  return user_type(type_nr,std::move(args));
}


@ A move constructor for |type_expr| is provided, but no copy constructor;
instead of the latter we provide a |copy| method explicitly creating a deep copy
value, to which one may then apply move-construction. We similarly provide a
move assignment operator, which implicitly makes the copy assignment operator
deleted (rather then implicitly provided), and a |swap| method (as
proof-of-concept, it is not currently used). The destructor for |type_expr| does
(by calling |clear|) take care to destruct the active field, explicitly calling
|delete| if that field is a pointer.

 Two special cases of (partial) move semantics are provided, for replacing an
undefined (or partially undefined) type by a more specific type (note that in
these cases nothing disappears, so no clean-up is necessary). The first case,
the |set_from| method, makes a shallow copy of its argument into the object for
which is was called, which is required to be undetermined initially. This method
that was present long before \Cpp11 allowed proper move semantics is in fact
used to implement the move constructor and move assignment operator. The second
case, the |specialise| method, is used during type analysis, to see if our type
matches a given pattern, or in case it was (partially) undefined whether it can
be made to match the pattern by if necessary replacing some undetermined
descendants by more specific ones. The call returns a value indicating whether
this was possible, and if so makes the necessary specialisations to our type.
That is done by copying, so the caller does not require, acquire, or lose
ownership of the pattern for/by this method.

@< Ordinary methods of the |type_expr| class @>=

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

@ Here are declarations of some more public methods of |type_expr|.

@< Ordinary methods of the |type_expr| class @>=
  bool operator== (const type_expr& y) const;
  bool operator!= (const type_expr& y) const @+{@; return not((*this)==y); }
  bool is_unstable() const;
    // whether |undetermined| components elements are present
@)
void print(std::ostream& out) const;

@ For the definition of |class type_expr@;| to be processed properly, we must
pay some attention to ordering of type definitions, because of the recursions
present. The structure |func_type| will contain instances of |type_expr|, so its
definition must follow, but since we used a pointer to such a structure, it must
be declared as a structure before the definition of |type_expr| is seen. For
type lists this kind of difficulty is taken care of because
|simple_list<type_expr>| could be used in a |typedef| before |type_expr| was a
complete type.

@< Type declarations @>=

struct func_type;

@ The variant for function types needs both an argument type and a result type,
which is obtained by having |type_expr| contain a raw pointer to a |func_type|
structure. There are occasions where outside code needs to refer to |func_type|,
so we cannot make this a local definition to the |type_expr| class. The
definition itself is simple, with all methods inlined. Like for |type_expr| we
do not provide an ordinary copy constructor, but a value-returning |copy|
accessor method.

@< Type definitions @>=
struct func_type
{ type_expr arg_type, result_type;
@)
  func_type(type_expr&& a, type_expr&& r)
@/: arg_type(std::move(a)), result_type(std::move(r)) @+{}
@/func_type(func_type&& f) = default; // move constructor
  func_type& operator=(func_type&& f) = default; // move assignment
  func_type copy() const // in lieu of a copy constructor
  {@; return func_type(arg_type.copy(),result_type.copy()); }
};

@ The factory method for function types could not be defined inside the
definition of |type_expr|, since |func_type| is not (and cannot be) a complete
type there. The code below is put in the header file after the type definitions,
so the |func_type| definition we just saw, and in particular its constructor
used here, will be fully known at that point.

@< Template and inline function definitions @>=
inline type_expr type_expr::function(type_expr&& arg_tp, type_expr&& res_tp)
{ type_expr result;
  result.tag=function_type;
  result.func_variant=new func_type(std::move(arg_tp),std::move(res_tp));
  return result;
}

@ Because a |type_expr| owns all of its descendants, copying one means making a
deep copy, and in order to avoid inadvertently invoking that operation, it is
provided by an explicitly called |copy| method rather than by a regular copy
constructor. If the copy needs to be held in an existing variable, the result of
|copy| can be move-assigned; this is cheap, and probably return value
optimisation from |copy| can avoid all overhead. Using a named function here
also results in the recursion being more visible that it would be for a regular
copy constructor. Since this is not a constructor, the pointers returned from
|new| are immediately owned when stored in a component of |result|, so we need
not worry about that. (We can also note that, even if this had been a
constructor, there would be no need for action to ensure exception safety, since
no exception can happen between storing the pointer and function termination.)

This is a straightforward structural recursion. We do not replace any tabled
types by their definitions, so there is no possibility of infinite recursion; we
do however continue to deep copy to the argument types, if any, of a tabled type
constructor. Note that any owned pointers of |result| are assigned just before
exit from the |switch|, and therefore just before |result.tag| gets assigned so
as to signal ownership of those pointers; no exception (which would cause a
memory leak) can happen between the two assignments. And if an exception should
occur before those pointers are installed, no attempt will be made to |delete|
any pointer fields while destructing |result| because one still has
|result.tag==undetermined_type| (this would be different if we instead had
written |switch (result.tag=tag)|, as we erroneously used to).

@:type expression copy@>

@< Function definitions @>=
type_expr type_expr::copy() const
{ type_expr result;
  switch (tag)
  { case undetermined_type: break;
    case primitive_type: result.prim_variant=prim_variant; break;
    case variable_type: result.typevar_variant=typevar_variant; break;
    case function_type: result.func_variant=new
      func_type(func_variant->copy());
    break;
    case row_type:
      result.row_variant=new type_expr(row_variant->copy());
    break;
    case tuple_type: case union_type:
      @< Assign a deep copy of |tuple_variant| to |result.tuple_variant| @>
    break;
    case tabled: result.tabled_variant.nr=tabled_nr();
      @< Assign a deep copy of |tabled_variant.type_args| to
         |result.tabled_variant.type_args| @>
  }
  result.tag=tag; // henceforth |result| owns added things
  return result;
}

@ The code below exemplifies looping over a list accessed by a raw
pointer-to-node like |raw_type_list|. Since we do not need to alter the
structure (or indeed anything) of the list looped over, we can do with a weak
iterator; we can create on from the raw pointer which is a matter of
repackaging, and then use the usual |simple_list| iteration syntax. For building
the copy we use an actual list container for exception safety; we choose
|dressed_type_list| (which is |sl_list<type_expr>|) to facilitate appending at
the end. Eventually we must |undress| to |simple_list|, and then |release| from
it the raw pointer to be assigned into |result|.

@< Assign a deep copy of |tuple_variant| to |result.tuple_variant| @>=
{
  dressed_type_list dst;
  for (wtl_const_iterator it(tuple_variant); not it.at_end(); ++it)
    dst.push_back(it->copy());
  result.tuple_variant = dst.undress().release();
  // incorporate and transfer ownership
}

@ Just as illustration that we can, we here do the same but using a |simple_list|
rather than |sl_list|.

@< Assign a deep copy of |tabled_variant.type_args| to
  |result.tabled_variant.type_args| @>=
{
  type_list dst; auto dit=dst.begin();
  for (wtl_const_iterator it(tabled_variant.type_args); not it.at_end(); ++it)
    dit = dst.insert(dit,it->copy());
  result.tabled_variant.type_args = dst.release();
  // incorporate and transfer ownership
}

@ The method |clear| does the work for the |type_expr| destructor, cleaning up
the nodes accessible from here before setting |tag=undetermined_type|. The
reason for separating this out from the destructor is that it allows variables
to be made ready for reuse. The recursion of this method is implicit, as the
|delete| calls may run destructors that call |clear| again.

If some variants of the |union| had been smart pointers, then we would have
needed to explicitly call their destructors, rather than |delete|. This would
not really have made the code simpler or more transparent.

@< Function definitions @>=
void type_expr::clear() noexcept
{ switch (tag)
  { case undetermined_type: case primitive_type: case variable_type: break;
    case function_type: delete func_variant; break;
    case row_type: delete row_variant; break;
    case tuple_type: case union_type: delete tuple_variant; break;
    case tabled: delete tabled_variant.type_args;
    // free our share of constructor arguments
  }
  tag = undetermined_type;
}

@ The method |set_from| makes a shallow copy of the structure; it implements
move semantics. Correspondingly it takes as argument another |type_expr| by
modifiable rvalue reference. The contents of its top-level structure is moved to
the current |type_expr|, and then set to a empty |undetermined_type| value,
which effectively detaches any possible descendants from it. This operation
requires that |type_expr| previously had |tag==undetermined_type|. We test this
condition using an |assert| statement (rather than throwing |logic_error|) to
honour the |noexcept| specification.

Using assignments (rather than placement-|new|) in |set_from| is possible since
all variants are POD types. Also we do not use move-assignments here (although
morally they are), because the data consists of either raw pointers or of basic
(integer) data, for which moving is not anything other than copying; the raw
pointers that are left in |p| will not get deleted because |p.tag| is reset to
|undetermined_type|.

The moving copy constructor is nearly identical to |set_from| method, so we call
the default constructor first, and then |set_from| to move the top-level value.

@h <stdexcept>
@< Function definitions @>=
void type_expr::set_from(type_expr&& p) noexcept
{ assert (tag==undetermined_type);
   // logic should ensure this; we promised not to |throw|
  switch(tag=p.tag) // copy top node
  { case undetermined_type: break;
    case primitive_type: prim_variant=p.prim_variant; break;
    case variable_type: typevar_variant=p.typevar_variant; break;
    case function_type: func_variant=p.func_variant; break;
    case row_type: row_variant=p.row_variant; break;
    case tuple_type: case union_type: tuple_variant = p.tuple_variant; break;
    case tabled: tabled_variant=p.tabled_variant;
  }
  p.tag=undetermined_type;
  // detach descendants, so |p.clear()| becomes a no-op
}
@)
type_expr::type_expr(type_expr&& x) noexcept // move constructor
: type_expr()
{@; set_from(std::move(x)); }

@ For move assignment we reuse the |set_from| method, and for |swap| we do a
move construction and two |set_from| calls, unless the |tag| fields match, in
which case we can call |std::swap| directly on the matching variant fields. The
|swap| method is not actually used (so we exclude its implementation from
compilation); it serves just as illustration of how it should be done if ever
needed.

@< Function definitions @>=
type_expr& type_expr::operator=(type_expr&& x) noexcept // move assignment
{ if (this!=&x)
  { clear(); // detach anything previously linked
    set_from(std::move(x)); // move top level structure
  }
  return *this;
}
@)
#if 0
void type_expr::swap(type_expr& other) noexcept
{
  if (tag==other.tag) // an easy case
    switch(tag)
    { case undetermined_type: break; // no need to swap |nothing| fields
      case primitive_type: std::swap(prim_variant,other.prim_variant); break;
      case variable_type:
        std::swap(typevar_variant,other.typevar_variant); break;
      case function_type: std::swap(func_variant,other.func_variant); break;
      case row_type: std::swap(row_variant,other.row_variant); break;
      case tuple_type: case union_type:
        std::swap(tuple_variant,other.tuple_variant); break;
      case tabled: std::swap(tabled_variant,other.tabled_variant); break;
    }
  else
  @/{@;
    type_expr t(std::move(other));
    other.set_from(std::move(*this));
    set_from(std::move(t));
  }
}
#endif

@ The |specialise| method is mostly used to either set (if initially
undetermined) the |type_expr| it is called for to a given |pattern| (which might
be a complete type, for instance when called from |conform_types| below), or to
test if it already matches that pattern. This is done so that the same code can
cater both for strong type contexts (where a value of a specific type must be
produced) and for weak ones (where during type analysis any type found will do),
as well as occasionally for intermediate requirements (like: we need a $3$-tuple
type with |int| as second component). If possible, the initial value of our
|type_expr| is replaced by a refinement (meaning that some undetermined parts
are substituted for) that is also a refinement of |pattern|. This foreshadows
the notion of unification in the context of second order types, but the purpose
of the |specialise| method is mainly this convenience of treating different
types of context at once.

In the case of an |undetermined_type|, |specialise| uses the |set_from| method
to make |*this| a copy of |pattern|. In the other cases we only continue if
the top levels of both type declarers match, in which case we try to
recursively specialise all descendants. We do not guarantee
commit-or-roll-back, in other words, when the specialisation fails, some
modifications to our type may still have been made. This is no problem in most
situations, since failure to specialise $t_1$ to $t_2$ will usually be
followed by an attempt to coerce $t_2$ to $t_1$, or by throwing of an error;
here any specialisation that brings $t_1$ closer to $t_2$ cannot be harmful
(it probably makes no difference at all). And should commit-or-roll-back be
important, we provide an accessor method |can_specialise| that can be tested
before calling |specialise| to avoid unwanted side-effects.

Both of types here can have undetermined parts, but we do not deal with
polymorphism here: a context never requires a polymorphic type, and if |pattern|
is polymorphic it must match an undetermined part of our type, which will then
be set from it. Thus all type variables are treated as primitive types.

@< Function definitions @>=
bool type_expr::specialise(const type_expr& pattern)
{ if (pattern.tag==undetermined_type)
    return true; // specialisation to \.* trivially succeeds.
  if (tag==undetermined_type) // specialising \.* also always succeeds,
    {@; set_from(pattern.copy()); return true; }
     // by setting |*this| to a copy of  |pattern|
  @< Handle cases of |specialise| where at least one type is tabled,
     but fall through if both are tabled with equal |tabled_nr()| @>
@)
  if (pattern.tag!=tag) return false;
    // now it is impossible to refine if tags mismatch
  switch(tag)
  { case primitive_type: return prim_variant==pattern.prim_variant;
    case variable_type: return typevar_variant==pattern.typevar_variant;
      // (fixed) type variable is like primitive
    case function_type:
      return func_variant->arg_type.specialise(pattern.func_variant->arg_type) @|
         and func_variant->result_type.specialise
                                          (pattern.func_variant->result_type);
    case row_type: return row_variant->specialise(*pattern.row_variant);
    case tuple_type: case union_type:
     @< Try to specialise types in |tuple_variant| to those in
        |pattern.tuple_variant|,
        and |return| whether this succeeded @>
    case tabled: assert(tabled_nr()==pattern.tabled_nr());
      // since we fell through above
     @< Try to specialise types in |tabled_args()| to those in
        |pattern.tabled_args()|,
        and |return| whether this succeeded @>
    default: assert(false);
      return true; // to keep the compiler happy, cannot be reached
  }
}

@ This method must pay some nontrivial attention to types with |tag==tabled|,
whose (possibly recursive) meaning is stored in |type_expr::type_map|. We fall
through this code only in the case mentioned in the module name Types or type
constructors that are entered into that table never have undetermined
components, so if |tag==tabled| holds for our own |type_expr|, then no
substitution for undetermined parts needs to be made. However, in order to
ensure for instance that, after a successful |specialise| against a ``row of''
pattern, our |component_type()| can be relied upon, we instead replace, when
|tag==tabled|, ourselves by an expanded version before trying again recursively.
In the unlikely case that our type is neither undetermined nor tabled, but we do
have |pattern.tag==tabled| (apparently the context had imposed some partial
structure, but the type check came up with a tabled type instead), we expand the
definition of |pattern| in a recursive call of specialise, just to make sure we
don't accidentally miss a case where the tabled type does match the context
requirement.

@< Handle cases of |specialise| where at least one type is tabled,
   but fall through if both are tabled with equal |tabled_nr()| @>=
{
  if (tag==tabled)
  { if (pattern.tag==tabled)
    { if (tabled_nr()!=pattern.tabled_nr())
      { if (is_recursive() and pattern.is_recursive()
            or top_kind()!=pattern.top_kind())
          return false;
        return expand().specialise(pattern.expanded());
      }
      else {} // equal tabled types: fall through
    }
    else // only |*this| is tabled
    {
      if (top_kind()!=pattern.tag)
        return false;
      return expand().specialise(pattern);
    }
  }
  else if (pattern.tag==tabled)
    return specialise(pattern.expanded());
}

@ For tuples and unions, specialisation is done component by component. As
before we use weak iterators to loop over the tuple components. The termination
condition is one of the lists running out, or a component specialisation
failing.

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

@ This is quite similar to the previous section, but argument lists must have
equal lengths.

@< Try to specialise types in |tabled_args()| to those in
|pattern.tabled_args()|, and |return| whether this succeeded @>=

{
  wtl_iterator it0(tabled_args());
  wtl_const_iterator it1(pattern.tabled_args());
  while (not it0.at_end() and not it1.at_end())
    if (not it0->specialise(*it1))
      return false;
    else @/{@; ++it0; ++it1; }
  assert (it0.at_end() and it1.at_end());
  return true;
}


@ The accessor |can_specialise| returns the same value as |specialise|, but
without any side effect. It turns out to be perfectly symmetric in |*this| and
|pattern|. If both are tabled, there can be no undetermined parts left, and we
forward to the equality test. For the rest this is just like |specialise|.

@< Function definitions @>=
bool type_expr::can_specialise(const type_expr& pattern) const
{ if (pattern.tag==undetermined_type or tag==undetermined_type)
    return true;
  if (pattern.tag==tabled)
  { if (tag==tabled) return *this==pattern;
    // let equality test deal with this
    else return can_specialise(pattern.expanded());
  }
  else if (tag==tabled)
    return expanded().can_specialise(pattern);
@)
  if (pattern.tag!=tag) return false; // impossible to refine
  switch(tag)
  { case primitive_type: return prim_variant==pattern.prim_variant;
    case variable_type: return typevar_variant==pattern.typevar_variant;
      // like primitive
    case function_type:
      return func_variant->arg_type.can_specialise
                                    (pattern.func_variant->arg_type) @|
         and func_variant->result_type.can_specialise
                                    (pattern.func_variant->result_type);
    case row_type:
      return row_variant->can_specialise(*pattern.row_variant);
    case tuple_type: case union_type:
      @< Find out and |return| whether we can specialise the types in
         |tuple_variant|
         to those in |pattern.tuple_variant| @>
    default: assert(false);
      return true; // to keep the compiler happy, cannot be reached
  }
}

@ This part too is like the one for |specialise|, but we can use |const|
iterators throughout.

@< Find out and |return| whether we can specialise the types in
   |tuple_variant| to those in |pattern.tuple_variant| @>=
{
  wtl_const_iterator it0(tuple_variant), it1(pattern.tuple_variant);
  while (not it0.at_end() and not it1.at_end()
         and it0->can_specialise(*it1))
    @/{@; ++it0; ++it1; }
  return it0.at_end() and it1.at_end();
  // whether both lists terminated
}

@ The equality test differs little from |can_specialise|. We no longer care
specially about undetermined parts (while not forbidden, they should simply
match up identically), and on the other hand we do some extra effort for the
case where both types are tabled, since that was forwarded here. In doing so we
treat user type constructors like any other type constructor (both the
constructor and all argument types must match), and simple tabled types are
treated like type constructors with $0$ arguments, meaning they only match when
identical. That this is correct, and avoids non-termination by repeated
expansion of tabled types, is ensured by the way we process type definitions
when entering them into the tables. Notably, accidental type equalities between
tabled types are detected and squashed at that time, and we ensure that all
intermediate types (even unnamed ones) along any recursive self-relation are
tabled. Thus we avoid the possibility of a loop in which alternately one and
then the other type is a (recursive) tabled one expanding, while the other is
non-tabled at that point.

@< Function definitions @>=
bool type_expr::operator== (const type_expr& y) const
{ if (tag!=y.tag)
  { if (tag!=tabled and y.tag!=tabled) return false; // different structures
    return tag==tabled ? expanded()==y  : (*this)==y.expanded();
  }
  switch (tag) // we know that tags are equal, branch on which they are
  { case undetermined_type: return true;
    case primitive_type: return prim_variant==y.prim_variant;
    case variable_type: return typevar_variant==y.typevar_variant;
    case function_type:
      return func_variant->arg_type==y.func_variant->arg_type @|
	 and func_variant->result_type==y.func_variant->result_type;
    case row_type: return *row_variant==*y.row_variant;
    case tuple_type: case union_type:
       @< Find out and |return| whether all types in |tuple_variant|
          are equal to those in |y.tuple_variant| @>
    case tabled:
       @< Find out and |return| whether type constructor call in
          |tabled_variant| is identical to that in |y.tabled_variant| @>
  }
  assert(false); return true; // cannot be reached, but compilers don't trust it
}

@ This module has a familiar structure, using weak iterators to loop over
|tuple_variant|.

@< Find out and |return| whether all types in |tuple_variant| are equal to
those in |y.tuple_variant| @>=
{
  wtl_const_iterator it0(tuple_variant), it1(y.tuple_variant);
  while (not it0.at_end() and not it1.at_end()
         and *it0==*it1)
    @/{@; ++it0; ++it1; }
  return it0.at_end() and it1.at_end();
  // whether both lists terminated
}

@ For tabled types or type constructor calls, we distinguish between recursive
ones and the rest. For recursive tabled type constructors, instances are
considered equal only if the constructors are identical, and the argument lists
are equal (in length and recursively for each component position, like for tuple
types). For non recursive tabled types or type constructors we use the same rule
when constructors are identical, but otherwise kick the can down the road by
expanding both sides and then recursively testing equivalence of the expansions
(similarly to what was done when only one type was tabled). In the absence of
recursion, the latter does not risk non-termination. The reason for not always
doing the same as in the recursive case is that different tabled type
constructors (even with different arity) can have common instances, the equality
of which would be missed if we insisted on having
|tabled_variant.nr==y.tabled_variant.nr|. We could detect the potential for this
when entering tabled types by testing if their polymorphic defining expressions
can be unified, but besides being a lot of effort, the information it provides
would not simplify the equality test. For recursively defined type constructors
there will be no such non obvious cases of equality, because by decree different
such constructors produce different types.

@< Find out and |return| whether type constructor call in |tabled_variant| is
   identical to that in |y.tabled_variant| @>=
{ if (tabled_variant.nr!=y.tabled_variant.nr)
    return not (is_recursive() and y.is_recursive())
           and expanded()==y.expanded();
  wtl_const_iterator it0(tabled_variant.type_args),
                   it1(y.tabled_variant.type_args);
  while (not it0.at_end() and not it1.at_end())
    if (*it0!=*it1)
      return false;
    else
      {@; ++it0; ++it1; }
  assert (it0.at_end() and it1.at_end());
  // both lengths match arity of the type constructor
  return true;
}

@ The predicate |is_unstable| checks for the presence of undefined components
(\.*) in a type (such types are ``unstable'' since they might get modified by
|specialise|, and in certain places their occurrence could undermine the
consistency of the type system).

@< Function definitions @>=
bool type_expr::is_unstable() const
{ switch (tag)
  {
    case undetermined_type: return true;
    case tabled:  return false;
      // syntax excludes \.* in tabled types and their arguments
    case primitive_type: return false;
    case variable_type: return false;
    case function_type:
      return func_variant->arg_type.is_unstable()
          or func_variant->result_type.is_unstable();
    case row_type: return(row_variant->is_unstable());
    case tuple_type: case union_type:
    {
      for (wtl_const_iterator it(tuple_variant); not it.at_end(); ++it)
        if (it->is_unstable()) return true;
      return false;
    }
  }
  assert(false); return true; // cannot be reached, but compilers don't trust it
}

@ We shall employ various forms of substitution for type variables, that are
variations of the method |type_expr::copy|, but which do something special when
encountering certain type variables. Most will be treated later, when we deal
with second order types, but |simple_subst| is used (amongst others) to
implement application of tabled type constructors, so we define it right away.
It takes a |type_expr| (the right hand side of a type constructor definition)
and a vector of argument expressions, and produces a copy of the former in which
the latter are substituted for type variables numbered from~$0$ upwards. There
is no restriction on the argument expressions, which are not inspected at all
after they have been substituted, and the only a restriction on~|tp| is that one
must not have |is_unstable(tp)|.

@< Declarations of exported functions @>=
type_expr simple_subst
  (const type_expr& tp, const std::vector<type_expr>& assign);

@ While |simple_subst| is a variation of |type_expr::copy|, it is not a method
of |type_expr|, so it cannot directly assign to private members. Instead, we use
the factory methods like |type_expr::row| to build our result on the way back in
a recursive traversal of the type, moving from the result of recursive calls (no
|std::move| is needed, since these results are not stored in variables). Any
type variables must be parameters of the type constructor, since such
constructors cannot be defined in the scope of any bound type variables; in
other words parameters are numbered from~$0$. We simply replace them by copies
the corresponding |type_expr| from |assign|, which are required to exist. The
|tabled| case should not arise when |simple_subst| is used to implement
expansion of a tabled type, but it is allowed when it serves to implement the
operation of a user defined type constructor. When this case arises,
we are therefore translating a user type constructor call from the identifier
table into a tabled constructor call, which is a form that type analysis
functions can deal with later if needed. All that then needs to be done here is
substitute for any type parameters that might occur in the argument list of the
tabled type, which is quite similar to the substitutions made into the
components of a tuple or union type.

@< Function definitions @>=

type_expr simple_subst
  (const type_expr& tp, const std::vector<type_expr>& assign)
{ switch (tp.raw_kind())
  { case primitive_type: return type_expr::primitive(tp.prim());
    case function_type: return
      type_expr::function(simple_subst(tp.func()->arg_type,assign),
                          simple_subst(tp.func()->result_type,assign));
    case row_type: return
      type_expr::row(simple_subst(tp.component_type(),assign));
    case tuple_type:
    case union_type:
    { dressed_type_list aux;
      for (wtl_const_iterator it(tp.tuple()); not it.at_end(); ++it)
        aux.push_back(simple_subst(*it,assign));
      return type_expr::tuple_or_union(tp.raw_kind(),aux.undress());
    }
    case tabled:
    { dressed_type_list aux;
      for (wtl_const_iterator it(tp.tabled_args()); not it.at_end(); ++it)
        aux.push_back(simple_subst(*it,assign));
      return type_expr::user_type(tp.tabled_nr(),aux.undress());
    }
    case variable_type:
  @/{@; auto c = tp.typevar_count();
      assert(c<assign.size());
      return assign[c].copy();
    }
    default: assert(false); // there should be no undetermined type components
  }
  return type_expr(); // cannot be reached
}


@*2 User defined, possibly recursive, types and type constructors.
%
We come to a new part of the |type_expr| type, a static member |type_map| that
allows for storage of ``tabled'' types and type constructors. These can
represent recursive (effectively infinite) type expressions, like a row-of type
whose component type is the type itself, or algebraic types like binary trees
with node labels of some type (the latter would be a recursive type
constructor). The mechanism is separate from the one used to associate types
with user defined type identifiers, of which it can be considered an
internalised form, made accessible to |type_expr| methods. Indeed entries of
|type_map| derive from user type definitions, and the identifier table will
equate such type identifiers to certain tabled types. We need to pre-declare
some types used in the declaration of |type_expr| methods.

@< Type declarations @>=
struct type_binding;
class type;

@ The |type_map| member is of a sub-class |type_expr::defined_type_mapping|,
which is basically a vector of |type_binding|. The variant |tabled_variant| of
|type_expr| will record an index into this vector (we ignore here another use
made of this variant in the parser, which is dealt with in the \.{global.w}
module). This class does not hide its data (though it does have one private
method |dissect_to|), but the unique object of this class is a private static
member of |type_expr|, so access is mostly controlled by static methods of that
class.

The way tabled type are used has evolved a bit since the introduction of user
type constructors. Before, the main importance was being able to do structural
type equivalence testing for all defined types at the point of definition, so
that during type analysis one is sure that distinct tabled types are not
equivalent, which makes it possible to have a finite equality test in the
presence of recursively defined types. We now isolate mutually recursive cliques
in type (constructor) definitions and simply decree the absence of equivalences
among different such cliques, so they function as primitive types or type
constructors; they are then exempt from expansion during testing for type
equivalence, which thus becomes a finite process. We do keep non recursive types
and type constructors in |type_expr::type_map| as well, making sure distinct
entries are not equivalent to speed up testing; it would be more work to
completely separate the recursive and non recursive parts of type definitions,
so we limit our efforts to flagging recursive entries and treating them slightly
differently in equality testing. Entries of |type_map| type also can have an
associated a field list, and provide for look up of union types, which is needed
in discrimination clauses.

@< Type definitions @>=
struct type_binding
{ static constexpr id_type no_id = -1;
  id_type name;
  unsigned short arity;
  bool recursive;
  type_expr tp;
  std::vector<id_type> fields;
  type_binding(type_expr&& t, unsigned short arity)
  : name(no_id), arity(arity), recursive(false), tp(std::move(t)), fields() @+{}
};
@)
class type_expr::defined_type_mapping : public std::vector<type_binding>
{ public:
  defined_type_mapping () : std::vector<type_binding>() @+{}
  const type_expr& definiens(type_nr_type i) const @+
    {@; assert (i<size()); return (*this)[i].tp; }
  unsigned short arity(type_nr_type i) const @+
    {@; assert (i<size()); return (*this)[i].arity; }
};

@~We need to define that declared static class member; it starts out empty.
@< Global variable definitions @>=
type_expr::defined_type_mapping type_expr::type_map;

@ A number of additional |static| methods of |type_expr| serve to regulate
access to the static class member |type_map|. While most of them simply serve as
a hatch (dutch: ``doorgeefluik'', no good English equivalent) to pass on
information, the method |matching_bindings| used to find fields to be used for a
given (supposedly tabled) tuple or union type is not entirely trivial, and the
method |add_typedefs| used to enter a list of newly defined (potentially
recursive) types into |type_map| is quite elaborate. Its argument is a list
|defs| of pairings of a type identifier to a type expression, the latter passed
by non-owning pointer. The potentially recursive nature of these definitions
lies in that they can not only refer, using the |tabled_variant|, to types
already defined in the mapping, but also to the types they define themselves.
For this purpose, those recursive type numbers start to count from
|type_expr::table_size()| as it is before |add_typedefs| method is called. The
return value is a list of the same length giving their type numbers after
applying type equivalencing~; usually these will be the same numbers used
initially, but some may have mapped to equivalent previously known types.

@< Static methods of |type_expr| that will access |type_map| @>=
static type_nr_type table_size();
static void reset_table_size(type_nr_type old_size);
static const std::vector<id_type>& fields(type_nr_type type_number);
static void set_fields (id_type type_number, std::vector<id_type>&& fields);
static sl_list<const type_binding*> matching_bindings (const type& tp);
static std::vector<type_nr_type> add_typedefs
 (const std::vector<std::pair<id_type,const_type_p> >& defs,
  unsigned int n_args);

@ Here are the easy ones among those methods: |table_size| just returns the
current |size| of |type_map| while |reset_table_size| shrinks the table back to
a previous size; |fields| and |set_fields| provide access to the list of field
names that can be associated to tuple and union types. (Note
that being a |static| member, we cannot |const| qualify |fields|, or any other
of these methods, to indicate that they leave |type_map| unchanged.) There
should be no method to remove the name of a tabled type, as doing so might lead
to non-termination of printing recursive types, but |reset_table_size| provides
a way to undo extensions of the table size; it must however be used only in
situations where no type expression outside the table can have captured the
names that are removed.

@< Function definitions @>=
type_nr_type type_expr::table_size() @+{@; return type_map.size(); }
void type_expr::reset_table_size(type_nr_type old_size)
{@; type_map.erase(std::next(type_map.begin(),old_size),type_map.end()); }
@)
const std::vector<id_type>& type_expr::fields(type_nr_type type_number)
{@; assert(type_number<table_size());
  return type_map[type_number].fields;
}
void type_expr::set_fields(id_type type_number, std::vector<id_type>&& fields)
{@; assert(type_number<table_size());
   type_map[type_number].fields=fields;
}

@ The method |matching_bindings| is used to find |fields| associated to a given
tuple or union type, in order to interpret a field assignment or a
discrimination clause, respectively; it returns a list non-owning pointers to
|type_binding| whose |tp| member can unify with |tp|. Since all type variables
in entries of |type_map| are polymorphic, but |tp| can have some fixed type
variables, the former need to be shifted if |tp.floor()>0|, so that all
polymorphic type variables start at the same level.


@< Function definitions @>=
sl_list<const type_binding*> type_expr::matching_bindings (const type& tp)
{ sl_list<const type_binding*> result;
  if (tp.floor()==0)
  { for (auto it=type_map.begin(); it!=type_map.end(); ++it)
      if (not it->fields.empty() and tp.has_unifier(it->tp))
        result.push_back(&*it);
  }
  else // we must shift any polymorphic variables in |it->tp|
  { for (auto it=type_map.begin(); it!=type_map.end(); ++it)
      if (not it->fields.empty() and tp.has_unifier(shift(it->tp,0,tp.floor())))
        result.push_back(&*it);
  }
  return result;
}

@ Here are the accessor methods for |type_expr| values that have
|raw_kind()==tabled|.

@< Function definitions @>=
id_type type_expr::type_name() const @+
{@; return type_map[tabled_variant.nr].name; }

unsigned short type_expr::tabled_arity() const
{ assert(tag==tabled and tabled_variant.nr<table_size());
@/return type_map[tabled_variant.nr].arity;
}

bool type_expr::is_recursive() const
{ assert(tag==tabled and tabled_variant.nr<table_size());
@/return type_map[tabled_variant.nr].recursive;
}

const type_expr& type_expr::top_level() const
{@; return type_map.definiens(tabled_variant.nr); }

@ The |expanded| methods requires more work, but a call to |simple_subst| does
the essential part. We just need to convert the type arguments from
|raw_type_list| to a |std::vector| to prepare for the call. In all cases we must
make a copy, so trying to save work when |arity==0| is not really worth it.

@< Function definitions @>=
type_expr type_expr::expanded () const
  // top level expansion of |tabled_variant|
{ if (tag!=tabled)
    return copy();
  const auto arity = tabled_arity();
  assert(length(tabled_args())==arity);
  std::vector<type_expr> assign;
  assign.reserve(arity);
  for (wtl_const_iterator it(tabled_args()); not it.at_end(); ++it)
    assign.push_back(it->copy());
  return simple_subst(type_map.definiens(tabled_nr()),assign);
}

@ The definition of |type_expr::add_typedefs| is subtle and requires quite a bit
of work. This is due to our requirements about tabled types~: for recursive
types, all types in their recursive cluster (mutual descendants) must be tabled
and marked as recursive, and among remaining tabled types and type constructors
there are no two equivalent ones. (The latter is no a strict necessity, since
tabled types not marked as recursive will be expanded if needed in the equality
test, so equivalence among them would ultimately be detected anyway.) We trade
off getting a simpler and more rapid equivalence test during type analysis for
the additional time spent in processing type definitions.

Type equivalence is defined in a conceptually simple way: two type expressions
are equivalent if they can be made identical by finitely many expansions of user
type and type constructor definitions. (This is somewhat stricter than the pure
structural equivalence that we used to use before the introduction of user type
constructors: (applications of) distinct recursive type (constructor)
definitions may give rise to types for which no difference can be exhibited at
any depth, but which are not equivalent under our current definition. Cases
where this would make a difference seem contrived though, so we do not
expect the change of definition to affect users; the new definition seems
easier to adapt to type constructors.)

Tabled type are introduced by calls of |add_typedefs|, each of which introduces
a group of type constructors of the same arity, namely the number of type
variables abstracted in the context of the grouped type definition (these type
variables denoting type parameters). There right hand sides can mutually refer
to types being defined, as well as using built-in or previously defined types
and type constructors of course. Since the type names being defined are not yet
treated as type constructors, mutual references do not take type arguments; when
later used, the arguments passed will be passed on unchanged in such mutual
references.

Grouped type definitions define an oriented graph, where each vertex is a type
subexpression of one of the right hand sides of the definitions. Edges are for
direct descendant relations, or for references within the definition group,
represented by |tabled| types with |tabled_nr()>=table_size()|. (The |tabled|
case with |tabled_nr()<table_size()| can also occur, as an application of an
already tabled type constructor; it is treated like built-in type constructors
with any argument types being type subexpressions.) This graph can serve to
detect actual recursion patterns as cliques of mutually two-way reachability
among vertices.

To process a group of type definitions, we first decompose the right hand side
type expressions into ``nodes'' for which all type subexpressions will (also)
become tabled references. Then we determine the recursive cliques among these
nodes, and mark all their members with the |recursive| flag. For the remaining
(not directly recursive) nodes we can use |type_expr::operator==| to test for
equality to weed out any equality among nodes or with previously tabled
constructors of the same arity. Finally we remove redundant tabled entries and
readjust references to entries that as a consequence have been shifted in the
table.

@h "preorder.h"

@< Function definitions @>=
std::vector<type_nr_type> type_expr::add_typedefs
  (const std::vector<std::pair<id_type,const_type_p> >& defs,
   unsigned int n_args)
{
  const type_nr_type n_defs=defs.size(), old_table_size=table_size();
  std::vector<type_data> type_array;

@/@< Copy right hand sides of |defs| type |type_array| and add entries for
     any of their component types that have themselves any descendants @>

  preorder::Preorder graph(type_array.begin(),type_array.end());
  auto cliques = graph.closure().cliques();

@/@< Add entries to |type_map| according to the entries of |type_array|,
     while setting the |recursive| flag for members of |cliques| @>

  std::vector<type_nr_type> relocate(type_array.size());
@/@< For new members of |type_map| that are not |recursive|, test whether they
     are equivalent any earlier type and if so set the corresponding entry of
     |relocate| @>
  @< Remove redundant types from |type_map|, and adjust |relocate|
     correspondingly @>
  @< Renumber the tabled entries according to |relocate| @>
  @< Copy type names from |defs| and truncate |relocate| to |n_defs| entries @>
  return relocate;
}

@ A first concern is which component types of the right hand sides should get a
separate slot in |type_array|; the minimal requirement is any intermediate type
in a recursion, i.e., one that both descends from some type being defined, and
that has the same type as descendant. As we are going to use |type_array| to
find out such recursive relations it is not practical to limit ourselves to such
types, so we shall make a slot for any type expression that has itself any
descendants: row, function, tuple and union types, and instances of already
defined tabled type constructors with |tabled_arity()>0| (whose argument types
we count as descendants).In order to easily see which direct descendants of a
node produce a link in |graph| we define the predicate |has_descendants|.

We do not accept the applications of previously tabled type constructors as
entries of |type_map|, to avoid having |type_expr::expanded| on one tabled type
(constructor) returning another one. When inserting into |type_map|, they will
be effectively be |expanded()| to a non tabled type (obtained by substitution
from the |definiens| of the previously tabled constructor). The unlikely
possibility that the expansion turns out to have no descendants will not break
anything, so we ignore it. (Actually, if a parameter type to the old constructor
would turn out to be unused, certain instances of this scenario would lead us to
believe there is type recursion via the argument list, when in fact there is
none; this would be very devious, and maybe we should exclude the possibility by
insisting that type constructors use all their type arguments at least once.)

@< Local function definitions @>=
bool has_descendants (const type_expr& tp)
{
  switch(tp.raw_kind())
  { default: return false;
    // |undetermined_type|, |primitive_type|, |variable_type|
  case row_type: case function_type: case union_type:
    return true;
  case tuple_type: return tp.tuple()!=nullptr;
    case tabled:
  return tp.tabled_nr()<type_expr::table_size() and tp.tabled_arity()>0;
  }
}

@ Now we can detail how we fill |type_array| and what we shall store there.
Since each entry corresponds type subexpression contained in |defs|, it is
practical to maintain a pointer to that |type_expr|, rather than to reconstruct
that information in some form. And since the main purpose of building this array
is to represent cross references between the entries as a graph, we shall
maintain a list of integers for the outgoing edges of this graph. Since the
edges correspond to those directly descendant types~|tp| of the node for which
|has_descendants(tp)| holds, this information suffices to produce an entry
of |type_map| corresponding to the node, with those descendants that are
also present in |type_map| replaced by |tabled| references.

@< Local type definitions @>=

struct type_expr::type_data
{ const_type_p tp;
  sl_list<unsigned short> out;
@)
  type_data() : tp(nullptr), out() @+{} // default creates empty slot
  type_data(const_type_p tp, sl_list<unsigned short>&& out)
  : tp(tp), out(std::move(out)) @+{}
  sl_list<unsigned short>::const_iterator begin() const
  @+{@; return out.begin(); }
  sl_list<unsigned short>::const_iterator end() const
  @+{@; return out.end(); }
};

@ We shall use some private methods to fill |type_array| and then |type_map|,
the first two of which are mutually recursive.

@< Private methods of the |type_expr| class @>=
  struct type_data;
  type_nr_type dissect_to(std::vector<type_data>& type_array) const;
  void record
    (std::vector<type_data>& type_array, sl_list<unsigned short>& out) const;
  type_expr rewrite
    (type_nr_type sz, unsigned int n_args,
     sl_list<unsigned short>::const_iterator& it) const;

@ The |type_expr| method |dissect_to| calls itself recursively on all its
components with descendants to ensure their presence in |type_array|, while
collecting the indices into that array that they return; finally it pushes the
|type_expr| itself with the accumulated list there, and returns its own index.
The list also collects the references to types currently being defined, in the
form of components that are |tabled| with |tabled_nr()| above the current table
size. Since all of this is done for components in all variants it gets a bit
repetitive, so we define another auxiliary method |record| to do the testing,
the recursive calling, and the accumulation of indices. The fact that recursive
calls are only done for actual type subexpressions of |*this|, guarantees
termination of the recursion.

@< Function definitions @>=
type_nr_type type_expr::dissect_to (std::vector<type_data>& type_array) const
{
  assert(has_descendants(*this)); // so we need to handle only those cases below
  sl_list<unsigned short>  out;
  switch(tag)
  { default: assert(false); // kinds that always fail |has_descendants()|
  break; case row_type: row_variant->record(type_array,out);
  break; case function_type:
    func_variant->arg_type.record(type_array,out);
    func_variant->result_type.record(type_array,out);
  break; case tuple_type: case union_type:
    for (wtl_const_iterator it(tuple_variant); not it.at_end(); ++it)
      it->record(type_array,out);
  break; case tabled:
    for (wtl_const_iterator it(tabled_args()); not it.at_end(); ++it)
      it->record(type_array,out);
  }
  const type_nr_type result = type_array.size();
  type_array.emplace_back(this,std::move(out));
  return result;
}
@)
void type_expr::record
  (std::vector<type_data>& type_array, sl_list<unsigned short>& out) const
{ if (has_descendants(*this))
    out.push_back(dissect_to(type_array));
  else if (raw_kind()==tabled and tabled_nr()>=table_size())
    out.push_back(tabled_nr()-table_size());
}

@ In order for references to types currently being defined to point correctly,
their entries must be at the first |n_defs| positions in |type_array|. We
therefore start placing that many empty slots at the start of |type_array|, so
that calls of |dissect_to| will correctly compute later indices. The we make
those calls for each right hand side in |defs|, which takes care of all
descendant types, but places the |type_data| for the right hand side at the end
of |type_array|; we then move it from there to the slot reserved for it at the
beginning. We have to be careful that |dissect_to| cannot be applied for types
without descendants, and although such types seem to be of limited utility as
right hand sides of type definitions, we must be prepared to handle them;
fortunately this is quite easy.

@< Copy right hand sides of |defs| type |type_array| and add entries for
   any of their component types that have themselves any descendants @>=
{
  type_array.resize(n_defs); // slots to be filled later
  for (type_nr_type i=0; i<n_defs; ++i)
    if (has_descendants(*defs[i].second))
    { defs[i].second->dissect_to(type_array);
    @/type_array[i] = std::move(type_array.back());
      type_array.pop_back();
    }
    else
      type_array[i] = type_data(defs[i].second,sl_list<unsigned short>());
}

@ When we install an element into |type_map| according to a |type_array|
entry~|e|, we build a |type_expr| modelled on |e.tp|, but with only direct
descendants. Such a descendant is a copy of the corresponding component of
|e.tp| whenever that component neither satisfies |has_descendants| nor is a
reference to a type currently being defined; for the remaining cases we build a
|tabled| reference to a type currently being added to |type_map|, with the
|tabled_nr()| determined by an entry of |e.out|, and |tabled_args()| set to a
standard list of |n_args| successive type variables (so that the type arguments
will later be passed around unchanged within the definition group). Such local
references arise in two ways: either form a component for which
|has_descendants| holds (here the destination must be obtained from |e.out|), or
from a component that already was such a local reference. We resist the
temptation to simply copy the component in the latter case, since we must attach
an argument list, and advance over an entry from~|e.out| even if we know what
number it is.

In order to organise this efficiently, we use an auxiliary method that produces
such a component type, which needs to know the arity |n_args| of the current
definition group, and uses an iterator |it| to fetch entries from |e.out|. Since
the |sl_list| iterators do not have post-increment operators, their handling
must be done in two statements.

@< Function definitions @>=
type_expr type_expr::rewrite
 (type_nr_type sz, unsigned int n_args,
  sl_list<unsigned short>::const_iterator& it) const
{
  if (has_descendants(*this))
@/{@; type_nr_type nr = sz+*it;
    ++it;
    return local_ref(nr,n_args);
  }
  if (raw_kind()==tabled and tabled_nr()>=sz)
  // this case has a link too
  { type_nr_type nr = sz+*it;
    ++it;
    // use and translate link that was recorded
    assert(nr == tabled_nr()); // this is what it was recorded from
    return local_ref(nr,n_args);
  }
  return copy();
  // remaining cases, including old tabled types with |arity==0|
}

@ Now the transfer of elements of |type_array| to |type_map| is fairly
straightforward: all descendants must be passed through |rewrite|, and then a
new type of the same |kind()| constructed from the results to be moved into
|type_map|. Note that these elements almost always have descendants, unless they
are copied as complete right hand side of a type definition, but since the
grammar forbids a single type identifier as such a right hand side (in order to
rule out circular type definitions), that can only be a primitive type. The only
complication in our code is the case where the type referenced from |type_array|
has |raw_kind()==tabled|, which must be a previously tabled type constructor
application, which application we must expand to prevent having a |type_map|
entry with |tag==tabled|.

@< Add entries to |type_map| according to the entries of |type_array|,
   while setting the |recursive| flag for members of |cliques| @>=
{ BitMap rec(type_array.size());
  for (const auto& clique : cliques)
    rec |= BitMap(type_array.size(),clique.wcbegin(),clique.wcend());
@)
  for (type_nr_type i=0; i<type_array.size(); ++i)
  {
    const auto& data = type_array[i];
    sl_list<unsigned short>::const_iterator oit=data.out.cbegin();
    type_expr tp;
    switch(data.tp->raw_kind())
    { default: tp = data.tp->copy();
    // rare, but |defs| may equate to a type without descendants
    break; case row_type:
      tp = type_expr::row
        (data.tp->component_type().rewrite(old_table_size,n_args,oit));
    break; case function_type:
      { auto tmp = data.tp->func()->arg_type.rewrite(old_table_size,n_args,oit);
        // sequence the |rewrite| calls
        tp = type_expr::function
           (std::move(tmp),
            data.tp->func()->result_type.rewrite(old_table_size,n_args,oit));
      }
    break;case tuple_type: case union_type:
      { dressed_type_list aux;
        for (wtl_const_iterator it(data.tp->tuple()); not it.at_end(); ++it)
          aux.push_back(it->rewrite(old_table_size,n_args,oit));
        tp = type_expr::tuple_or_union(data.tp->raw_kind(),aux.undress());
      }
    break; case tabled:
      @< Apply |rewrite| to every type in |data.tp->tabled_args()|,
         then assign to |tp| the expansion of the tabled type number
         |data.tp->tabled_nr()| to that list of argument types @>
     }
    type_map.emplace_back(std::move(tp),n_args);
    if (rec.isMember(i))
      type_map.back().recursive=true;
  }
}

@ The case where an entry of |type_array| has |tag==tabled| when a previously
tabled type constructor with positive arity is applied. The case had to be taken
into account because the result of the application can have one or more of the
types currently being defined as descendants, if they are passed via type
arguments to the constructor, and thus must be aware of such relations when
identifying recursive cliques. The more complicated way that |tp| is obtained
here can result in a |type_map| entry that is somewhat more complicated than one
would get directly, for instance if we apply a previous recursive type
constructor, we necessarily get components that are a type constructor with a
non-standard argument list: one containing other tabled references rather than
just type variables. It does not seem like this can cause problems.

@< Apply |rewrite| to every type in |data.tp->tabled_args()|... @>=
{ dressed_type_list aux;
  for (wtl_const_iterator it(data.tp->tabled_args()); not it.at_end(); ++it)
    aux.push_back(it->rewrite(old_table_size,n_args,oit));
  tp = type_expr::user_type(data.tp->tabled_nr(),aux.undress()).expanded();
}

@ Once the type definitions are installed in |type_map| with the necessary
|recursive| attributes, we can safely apply equality testing among types to them
(because that test stops at tabled types flagged |recursive|). This makes it
fairly easy to eliminate any duplicates.

@< For new members of |type_map| that are not |recursive|, test whether they
   are equivalent any earlier type and if so set the corresponding entry of
   |relocate| @>=
{
  for (type_nr_type nr = 0; nr<type_array.size(); ++nr)
  { const auto cur = nr+old_table_size;
    relocate[nr] = cur; // by default refer to ourselves
    if (not type_map[cur].recursive)
      for (type_nr_type k=0; k<cur; ++k)
        if (type_map.arity(k)==n_args and
            type_map.definiens(k)==type_map.definiens(cur))
        @/{@; relocate[nr]=k;
          break;
        }
  }
}

@ We traverse |relocate|, and for each entry not pointing to itself, we remove
the corresponding entry from |type_map|. We keep the list |relocate| of current
locations (in |type_map|) of entries that were constructed from those of
|type_array| up to date as follows. Entries are maintained if their initial
|relocate| value points to themselves; however, if the number |removed| of
earlier entries not maintained is positive, we move it down that many places,
and also decrease it address by |removed| to signal the displacement.
If an entry was found equal to a previous type, as witnessed by its initial
|relocate| value, is not maintained, and |removed| is increased; also, if that
|relocate| value points to an entry from the current definition group, it might
have been shifted down by in the current loop, so we replace our |relocate|
value by the current one for the entry pointed to.

@< Remove redundant types from |type_map|, and adjust |relocate|
   correspondingly @>=
{
  assert(type_map.size()==old_table_size+type_array.size());
  auto it = type_map.begin()+old_table_size; type_nr_type removed = 0;
  for (type_nr_type nr = 0; nr<type_array.size(); ++nr)
  { if (relocate[nr]!=old_table_size+nr)
    {
      ++removed; // forget |type_map[old_table_size+nr]|
      if (relocate[nr]>=old_table_size) // then maybe our twin has moved since
        relocate[nr] = relocate[relocate[nr]-old_table_size];
        // so forward to its current home
    }
    else if (removed==0)
      ++it; // include, but no copying is needed yet
    else
    {@; *it++ = std::move(type_map[old_table_size+nr]);
      relocate[nr] -= removed;
    }
  }
  type_map.erase(it,type_map.end());
}

@ Renumbering the tabled entries is not algorithmically hard, but requires a bit
of preparation to do conveniently. We want to replace certain type numbers~|nr|
by |relocate[nr-old_table_size]| inside tabled types, which is best done by a
method of |type_expr| (so as to have non-const access to the |tabled_variant|)
that we shall call |fix|. It needs to access to some variables local to the
method |add_typedefs| that we are (still) defining; we pack them into a
structure |fix_data| to be passed as argument. (A \Cpp11 lambda expression could
more easily capture local variables, but would lack the privilege to access a
|tabled_variant| field.)

@< Private methods of the |type_expr| class @>=
struct fix_data
{@; const std::vector<type_nr_type>& rel;
  type_nr_type sz;
};
void fix(const fix_data& f);

@~So here is how the replacement is done. Since |relocate| only applies to new
additions to |type_map|, indexing is shifted by |sz| (which will hold
|old_table_size|, and we only do anything for tabled types whose |tabled_nr()|
is at least |sz|. When the condition is met we replace the |nr| field by
|rel[nr-sz]|. There is a subtlety for previously tabled types: even though their
|type_nr()| should not change, they can have instances of types currently being
added to |type_map| in their |tabled_args()|, and these instances need to be
renumbered.

@< Function definitions @>=
void type_expr::fix(const fix_data& f)
{ if (tag!=tabled)
    return;
  if (tabled_variant.nr>=f.sz)
    tabled_variant.nr=f.rel[tabled_variant.nr-f.sz];
  else
    for (wtl_iterator it(tabled_args()); not it.at_end(); ++it)
      it->fix(f);
}

@ With that in our baggage, let us set to doing the renumbering. Since old
|type_map| entries cannot refer to new ones, we only run over the latter. After
preparing our |fix_data| as expected from |relocate| and |old_table_size|, we
apply a |fix| to every descendent |type_expr| of the |type_map| entry. The
|tabled| case should never arise as entry of |type_map|, and this means that we
never need to care about argument type lists of tabled constructors.

@< Renumber the tabled entries according to |relocate| @>=
{ fix_data r{relocate,old_table_size};
  for (auto it = type_map.begin()+old_table_size; it!=type_map.end(); ++it)
  { auto& tp = it->tp;
    switch(tp.tag)
    { default: // leave types without descendants unchanged
    break; case row_type: tp.row_variant->fix(r);
    break; case function_type:
      tp.func_variant->arg_type.fix(r);
      tp.func_variant->result_type.fix(r);
    break; case tuple_type: case union_type:
      for (wtl_iterator it(tp.tuple_variant); not it.at_end(); ++it)
        it->fix(r);
    break; case tabled: assert(false); // a |definiens| never should be |tabled|
    }
  }
}

@ And here are the final steps of |add_typedefs|. We do not set the |fields|
field of any defined type; that is handled by our caller, in \.{global.w}, using
the |set_fields| method.

@< Copy type names from |defs| and truncate |relocate| to |n_defs| entries @>=
{
  for (unsigned int i=0; i<n_defs; ++i)
  { type_nr_type nr= relocate[i];
    if (type_map[nr].name==type_binding::no_id
        // don't overwrite an existing type name
       @+and type_map[nr].tp.tag!=primitive_type) // nor name a primitive type
      type_map[nr].name=defs[i].first;
      // but otherwise insert type name into |type_map|
  }
  relocate.erase(relocate.begin()+n_defs,relocate.end()); // truncate
}


@*2 Printing types.
%
For printing types, we pass |type_expr| values to the operator~`|<<|' by
constant reference; passing by pointer could be defined, but we avoid this in
general, as would override a language-provided definition that prints a
hexadecimal value, which (without being very useful) would prevent the compiler
from signalling a possibly forgotten overload definition. Since we often hold
types in |type_ptr| values, this does mean the we must then dereference
explicitly when calling the operators below.

@< Declarations of exported functions @>=
std::ostream& operator<<(std::ostream& out, const type_expr& t);
std::ostream& operator<<(std::ostream& out, const func_type& f);

@~Printing types is relegated to the |type_expr::print| method. Before we define
it, here are two auxiliary function for printing sequences of types (inside
tuple and union types) and function types. For the latter we suppress
additional parentheses around argument and result types in case these are
tuple or union types.

@h "sl_list.h"
@< Function definitions @>=

std::ostream& operator<<(std::ostream& out, const type_expr& t)
{@; t.print(out); return out; }
@)
void print(std::ostream& out, const_raw_type_list l,char sep)
{ wtl_const_iterator it(l);
  if (not it.at_end())
    while (out << *it, not (++it).at_end())
      out << sep;
}

std::ostream& operator<<(std::ostream& out, const func_type& f)
{
  out << '(';
  if (f.arg_type.raw_kind()==tuple_type or
      f.arg_type.raw_kind()==union_type and f.arg_type.tuple()!=nullptr)
     print(out,f.arg_type.tuple(),@|f.arg_type.raw_kind()==tuple_type?',':'|');
     // naked tuple or union
  else out << f.arg_type; // other component type
  out << "->";
  if (f.result_type.raw_kind()==tuple_type or @|
      f.result_type.raw_kind()==union_type and f.result_type.tuple()!=nullptr)
     print(out,f.result_type.tuple(),
           f.result_type.raw_kind()==tuple_type?',':'|');
     // tuple, union
  else out << f.result_type; // other component type
  return out << ')';
}

@ And here is the |type_expr::print| method. Note that when |tag==tabled| we
print its |type_name()| if one is provided; if not then we expand this (unnamed,
yet present in |type_map|) type. Any recursion must pass through a named type,
so infinite recursion should not be possible.

@< Function definitions @>=

void type_expr::print(std::ostream& out) const
{ switch(tag)
  { case undetermined_type: out << '*'; break;
    case primitive_type: out << prim_names[prim()]; break;
    case variable_type: out << static_cast<char>('A' + typevar_variant); break;
    case function_type: out << *func_variant; break;
    case row_type: out << '[' << *row_variant << ']'; break;
    case tuple_type:
      if (tuple_variant==nullptr)
        out << "void";
      else
      {@;
         interpreter::print(out << '(', tuple_variant,',');
         out << ')';
      }
    break;
    case union_type:
      if (tuple_variant==nullptr)
        out << "(*)"; // this should not really occur
      else
      {@;
         interpreter::print(out << '(', tuple_variant,'|');
         out << ')';
      }
    break;
    case tabled:
      if (type_map[tabled_variant.nr].name!=type_binding::no_id)
      {
        out << main_hash_table->name_of(type_name());
        if (tabled_variant.type_args!=nullptr)
        {@;
           interpreter::print(out << '<', tabled_variant.type_args,',');
           out << '>';
        }
      }
      else expanded().print(out); // expand out when no identifier is attached
    break;
  }
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
type_ptr mk_type_variable(unsigned int i);
type_ptr mk_function_type(type_expr&& a, type_expr&& r);
type_ptr mk_row_type(type_expr&& c);
type_ptr mk_tuple_type (type_list&& l);
type_ptr mk_union_type (type_list&& l);
@)
type_p make_prim_type(unsigned int p);
type_p make_type_variable(unsigned int i);
type_p make_function_type(type_p a,type_p r);
type_p make_row_type(type_p c);
type_p make_tuple_type(raw_type_list l);
type_p make_union_type(raw_type_list l);
type_p make_tabled_type(id_type nr,raw_type_list l);
@)
raw_type_list make_type_singleton(type_p raw);
raw_type_list make_type_list(raw_type_list l,type_p t);
type_p unmake_type_singleton(raw_type_list l);

@ The functions like |mk_prim| below simply call the corresponding factory
method in the context of the function |std::make_unique|, which invokes |new|
and captures the resulting pointer, transforming it into a |type_ptr| result.
The functions like |make_prim| wrap them for the parser interface, undressing
the pointers to raw again. Using this setup it is ensured that all pointers are
considered owning their target, implicitly so while being manipulated by the
parser (which guarantees that every pointer placed on the parsing stack will be
argument of an interface function exactly once, possibly some |destroy| function
in case it pops symbols during error recovery).

These functions make some provisions to facilitate special relations between
types, so that the parser does not have to deal with them. Thus |mk_prim_type|
will return an empty tuple type when called with the type name for |"void"|,
although this is not a primitive type (indeed |"void"| should probably better be
handled just like user-defined type abbreviations). The function |mk_union_type|
will not encapsulate a list of length one into an invalid |type_expr|, but
rather return its unique list element, unpacked. The reason for this is that
having a syntactic category for a list of at least one component separated by
vertical bars simplifies the grammar considerably; when the list has one
component (and no bars) a union of one variant is temporarily created, but upon
incorporation into an encompassing type, the union is then removed again. For
tuples a similar provision is not necessary, as a somewhat more involved set of
grammar rules is used that avoids making type lists of length~$1$.

The functions |make_tuple_type| and |make_union_type| reverse the list of
types they handle, to compensate for the fact that the left-recursive grammar
rules (easier for the parser generator) construct their type lists in reverse
order.

@< Function definitions @>=
type_ptr mk_prim_type(primitive_tag p)
{ return p<nr_of_primitive_types ?
    std::make_unique<type_expr>(type_expr::primitive(p)) :
    type_ptr(mk_tuple_type(empty_tuple()));
}

type_ptr mk_type_variable(unsigned int i)
{@; return std::make_unique<type_expr>(type_expr::variable(i)); }

type_ptr mk_function_type (type_expr&& a, type_expr&& r)
{@; return std::make_unique<type_expr>
        (type_expr::function(std::move(a),std::move(r))); }

type_ptr mk_row_type(type_expr&& t)
{@; return std::make_unique<type_expr>(type_expr::row(std::move(t))); }

type_ptr mk_tuple_type (type_list&& l)
{@; return std::make_unique<type_expr>(type_expr::tuple(std::move(l))); }

type_ptr mk_union_type (type_list&& l)
{ if (l.singleton())
    return std::make_unique<type_expr>(std::move(l.front()));
  return std::make_unique<type_expr>
    (type_expr::tuple_or_union(union_type,std::move(l)));
}

type_ptr mk_tabled_type(type_nr_type nr)
{@; return std::make_unique<type_expr>(type_expr::user_type(nr,type_list())); }

@ A second group of functions similarly constructs type bottom-up, but since they
are intended for use by the parser, they return raw pointers (which the parser
treats as owning their referent). Because of this intended use, they use the raw
|type_p| rather than the smart |type_ptr|, both in results and in arguments.

@< Function definitions @>=

type_p make_prim_type(unsigned int p)
{@; return mk_prim_type(static_cast<primitive_tag>(p)).release(); }

type_p make_type_variable(unsigned int i)
{@; return mk_type_variable(i).release(); }

type_p make_function_type(type_p a,type_p r)
{@; return
    mk_function_type(std::move(*type_ptr(a)),std::move(*type_ptr(r))).release();
}

type_p make_row_type(type_p c)
{@; return mk_row_type(std::move(*type_ptr(c))).release(); }

type_p make_tuple_type(raw_type_list l)
{@; type_list result(l); result.reverse();
    return mk_tuple_type(std::move(result)).release();
}

type_p make_union_type(raw_type_list l)
{@; type_list result(l); result.reverse();
    return mk_union_type(std::move(result)).release();
}
@)
type_p make_tabled_type(id_type id,raw_type_list l)
{ type_list args(l); args.reverse();
  return new type_expr(type_expr::user_type(id,std::move(args)));
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

@ We define a function |unmake_type_singleton| to undo |make_type_singleton|
above. It goes in tandem with |unmake_pattern_singleton| defined in the
\.{parsetree.w} module, and the technical reason for needing them is explained
there.

@< Function definitions @>=
type_p unmake_type_singleton(raw_type_list p)
{@; type_list l(p); assert(l.singleton());
  return new type_expr(std::move(l.front()));
}


@*1 Second order types and type unification.
%
In this section we consider questions relating to universally quantified types
and type constructors. Both notions are related to the possibility of describing
families of types using a type with type variables, which can then be subject to
substitution of concrete types for them. This is to be distinguished from the
use of type patterns to represent incomplete type requirements (like ``a
$3$-tuple is needed here'') during type checking, where undetermined positions
are simply holes that can be arbitrarily and independently filled; type
variables can, and typically do, occur more than once in a type and have to be
systematically substituted for. We can write for instance $([[T]]\to[T])$ for
arbitrary type $T$, as type of the row-of-rows concatenation function. Type
constructors can also be represented by types using type variables taken from an
ordered list of one or more ``parameter types'', and when a corresponding list
of ``actual type arguments'' is given, substitution defines the resulting
concrete type.

Having functions with polymorphic types generalises what was already done for
certain general built-in functions, like the mentioned concatenation operation,
in pre-2.0 versions of the \axis. language. Their type checking, notably in
(overloaded) function calls, was handled by \foreign{ad hoc} code that is no
longer necessary. There is also some polymorphism involved in certain
non-function values like the empty list; this used to be handled in an
unsatisfactory way, but now we can handle values of polymorphic types like $[T]$
and $(\&{void}\mid T)$ normally, although we disallow assignable variables of
such types.

Our new type system much resembles in flavour (and is inspired by) the
Hindley-Milner type system used for instance in the Haskell language, but with
important differences, stemming mostly from the fact that we don't care much
about deducing types of \emph{all} expressions without any help from the user;
we have always required function parameters to be given explicit types. Indeed
such deduction might not be possible in the presence of certain features; it
notably does not blend well with``ad hoc'' function overloading (where a same
symbol or identifier can have multiple meanings to choose from) which was in use
long before we had second order types. We also allow declaring recursive types,
which Hindley-Milner typing does not allow (it gives a singular role to
injection maps into a disjoint union type, called constructors of an algebraic
data type, which provides a way around this restriction in cases that are
important in practice).

We already introduced (numbered) type variables into the definition of
|type_expr|. In type of expressions these arise either through explicit
introduction by the user (abstraction clauses), or from using names with a
polymorphic type. In the former case the type variables start out as
``fixed by the context'' within their scope (nothing can be assumed about them
and no substitution made) but when emerging from such a scope, the type
variables in its type become polymorphic (implicitly universally quantified).
We have chosen to use just one kind of type variables, with a frontier between
fixed and polymorphic types determined by the lexical context. Polymorphic type
variables have no distinct identity, and can be freely renumbered, for instance
to avoid a clash with other type variables when moving to a different context.

Polymorphic types require various forms of substitution operations, and a new
process called unification, which in its simplest form answers the question
whether two given polymorphic types admit any types obtainable from both by
substitution. If so, one usually also wants to know the simplest (most general)
such substitution. For instance, with $f,g,h$ denoting distinct type
constructors of one argument and with type variables $S,T$, the expressions
$f(S)$ and $f(g(T))$ can be made equal by the substitution $S:=g(T)$, and this
is the most general such substitution (with less general solutions being those
which in addition substitute some specific type for $T$). On the other hand, if
the second expression had instead been just $g(T)$, the substitution problem
would have had no solutions, as we cannot make the difference between the top
level $f$ and $g$ go away. The initial problem could arise in the context of
trying to match a function with type (for all~$S$) $(f(S)\to{h(S)})$ when called
with an argument expression of type (for all~$T$) $f(g(T))$; here the match
succeeds, and the function application so formed then has type (for all~$T$)
$h(g(T))$. Here the type of the function call reuses the type variable $T$ from
the type of the argument; in general it could involve type variables coming both
from the function and the argument types.

Although we allow recursive user-defined types that represent infinitely
repetitive type expressions, we do not allow unification to invent such types in
order to produce a solution when none exists otherwise. For instance, applying a
function of type $((S\to[S])\to([S]\to[S]))$ to the identity function, of type
$(T\to{T})$, will fail because we cannot make the argument type of $T$ equal
both to $S$ and to $[S]$; inventing a recursive type $L$ with $L=[L]$, and then
substituting $T:=(L\to{L})$, is rejected. Allowing it would complicate
unification immensely for no good purpose; if really needed, a user can
explicitly demand such a type substitution.

@ We come to type transformation and substitution operations. The argument~|tp|
to these functions, while specified as |type_expr|, should not contain any
|undetermined_type| subexpressions; to achieve this, a caller may need to first
replace each such occurrence by a fresh (only used at that place) type variable;
the factory function |type::wrap| given later does such a replacement.

The utility function |shift| renumbers type variables by adding |amount|, but
(in order to respect variables that are fixed in the context) leaves those
numbered less than |fix| unchanged.

The other two functions declared here relate to the unification process of
polymorphic types. They use an auxiliary structure |type_assignment|, defined
below, which specifies the type assignments to perform; it will be filled during
the unification by |can_unify|, and can then later be used in |substitution| to
obtain the unified type. It can however also be used with related type
expressions. For instance in a function application we apply unification to the
argument part of a function type and the type of a concrete argument, and once
the match is made, the |type_assignment| can then be substituted into the result
part of the function type to get the type of the call.

Due to the nature of |can_unify|, the type it assigns in |assign| to a type
variable may itself be polymorphic and subject to later substitutions;
therefore, and in contrast to |simple_subst|, the |substitution| function
continues to traverse any substituted type subexpression for further
substitutions (with measures taken to ensure termination of the process). The
final parameter |shift| of |substitution| indicates an amount by which all type
variables in |tp| are shifted up before being looked up in |assign|. This is in
order to accommodate situations in which a shift was applied to a type
expression before computing the |type_assignment| (to avoid clashes with type
variables in a type expression being unified with), in which case one can later
pass the |amount| of shift that was applied to |substitution|, which avoids
having to shift that (or a related) type expression again when using the
assignment.

@< Declarations of exported functions @>=
type_expr shift (const type_expr& tp, unsigned int fix, unsigned int amount);

@)
type_expr substitution
  (const type_expr& tp, const type_assignment& assign, unsigned int shift=0);
bool can_unify(const type_expr& P, const type_expr& Q, type_assignment& assign);


@ When we need to renumber the polymorphic type variables in one type
to avoid collision with variables bound in another type, we can use the
auxiliary recursive function |shift|, which is similar to |simple_subst|.
Type variables below the threshold |fix| are unchanged.

@< Function definitions @>=
type_expr shift
  (const type_expr& t, unsigned int fix, unsigned int amount)
{ type_ptr result;
  switch (t.raw_kind())
  { case primitive_type: return type_expr::primitive(t.prim());
    case function_type: result =
      mk_function_type(shift(t.func()->arg_type,fix,amount),
                       shift(t.func()->result_type,fix,amount));
    break;
    case row_type:
      result = mk_row_type(shift(t.component_type(),fix,amount));
    break;
    case tuple_type:
    case union_type:
    { dressed_type_list aux;
      for (wtl_const_iterator it(t.tuple()); not it.at_end(); ++it)
        aux.push_back(shift(*it,fix,amount));
      return type_expr::tuple_or_union(t.raw_kind(),aux.undress());
    }
    case tabled:
    { dressed_type_list aux;
      for (wtl_const_iterator it(t.tabled_args()); not it.at_end(); ++it)
        aux.push_back(shift(*it,fix,amount));
      return type_expr::user_type(t.tabled_nr(),aux.undress());
    }
    case variable_type:
  @/{@; auto c = t.typevar_count();
      return type_expr::variable(c<fix ? c : c+amount);
    }
    default: assert(false);
  }
  return std::move(*result);
}

@ Unification produces a set of substitutions for type variables, whose
structure is a bit more complicated than the simple list of |type_expr| values
that was used in |simple_subst|. We define a class |type_assignment| to hold the
relevant information. In fact it still holds a list of |const_type_ptr| values,
but we add a threshold |var_start| where the numbering of the corresponding type
variables starts. Also, the default value for these pointers is |nullptr|,
indicating that no substitution is to be made for the type variable. When a
substitution is called for, the type expression substituted may contain other
type variables, either below or above the threshold; in the latter case these
type variables may need further substitution. The absence of cycles that would
prevent termination of the substitution process is a class invariant.

Some methods are provided to easily assign a type expression to a type variable,
test whether this has been done, and as help in performing substitutions, to
renumber type variables taking into account that those that have been assigned
to will disappear from the type by the substitution.

@< Type definitions @>=
class type_assignment
{ unsigned int threshold; // first |variable_type| number that may vary
  std::vector<const_type_ptr> equiv; // variable substitutions go here
public:
  type_assignment(unsigned int fix_nr, unsigned int var_nr)
  : threshold(fix_nr)
  , equiv(var_nr)
    // default initialises |equiv| entries, we cannot say ``to |nullptr|''
  {}
@)
  type_assignment copy() const;
  unsigned int var_start() const @+{@; return threshold; }
  unsigned int size() const @+{@; return equiv.size(); }
  bool empty() const
  { auto null = @[ [](const const_type_ptr& p){@; return p==nullptr; } @];
    return all_of(equiv.begin(),equiv.end(),null);
  }
  void grow(unsigned int n) @+{@; equiv.resize(size()+n); }
@)
  const_type_p equivalent (unsigned int i) const
  {@; return i<threshold ? nullptr :
      (assert(i<threshold+size()),equiv[i-threshold].get()); }
  bool set_equivalent(unsigned int i, type_ptr&& p);
  bool set_equivalent(unsigned int i, const_type_p p);
@)
  unsigned int renumber(unsigned int i) const
  // reflect removal of variables assigned here
  { for (unsigned int j=i; j-->threshold; )
    // backwards to avoid using |i| after initialisation
      if (equiv[j-threshold]!=nullptr)
        --i; // take into account removal of |j|
    return i;
  }
@)
private:
  bool is_free_in(const type_expr& tp, unsigned int nr) const;
};

@ Sometimes we want to make a copy of a |type_assignment|, and since the |equiv|
filed holds unique pointers, we need to clone the type they point to.

@< Function definitions @>=
type_assignment type_assignment::copy() const
{ type_assignment result(threshold,size());
  result.equiv.reserve(size());
  for (const auto& tp : equiv)
    result.equiv.emplace_back(tp==nullptr ? nullptr : new type_expr(tp->copy()));
  return result;
}

@  The test for non-containment of a type variable is done by a private
recursive method |is_free_in|. In case a type variable is encountered for which
an equivalent is present in |equiv|, the recursion continues into the type
expression it is equated to; since we avoid situations where that expression
directly or indirectly references the same type variable, termination of this
process is ensured.

@< Function definitions @>=
bool type_assignment::is_free_in(const type_expr& tp, unsigned int nr) const
{ switch(tp.raw_kind())
  {
  case variable_type:
    if (@[auto* p=equivalent(tp.typevar_count())@;@])
      return is_free_in(*p,nr);
      // unpack assigned type variable, and retry
    return tp.typevar_count()==nr;
  case row_type: return is_free_in(tp.component_type(),nr);
  case function_type: return
    is_free_in(tp.func()->arg_type,nr) or
    is_free_in(tp.func()->result_type,nr);
  case tuple_type: case union_type:
    for(wtl_const_iterator it(tp.tuple()); not it.at_end(); ++it)
      if (is_free_in(*it,nr))
        return true;
    return false;
  default: // |undetermined| (shouldn't happen), |primitive_type|, |tabled|
    return false;
  }
}

@ There are two versions of |set_equivalent|, the first of which can move from
its type argument, while the second makes a copy of it. Both use |is_free_in| to
test for direct or indirect self-reference that would be introduced by the
assignment; if it is found, the assignment is not done, so as to preserve the
class invariant. These methods return a Boolean value telling whether setting
the type variable succeeded.

@< Function definitions @>=
bool type_assignment::set_equivalent(unsigned int i, type_ptr&& p)
{ assert (i>=threshold and i<threshold+size());
@/if (is_free_in(*p,i))
    return false;
  equiv[i-threshold]=std::move(p);
  return true;
}
bool type_assignment::set_equivalent(unsigned int i, const_type_p p)
{ assert (i>=threshold and i<threshold+size());
@/if (is_free_in(*p,i))
    return false;
  equiv[i-threshold]=std::make_unique<type_expr>(p->copy());
  return true;
}


@ We next give the |substitution| function. It essentially copies |tp| with
substitution according to |assign|, recursively propagated into the substituted
type expressions. The final argument |shift_amount| effectively renumbers all
the type variables in |tp| by adding |shift_amount| before looking up their
equivalents in |assign| (we interpret all type variables of |tp| as polymorphic,
since in calls with |shift_amount>0| we always have that |tp| is a subexpression
of a type in the overload table, in which no fixed type variables can occur).

For |tabled| types we just return a copy of the type, without expansion, which
avoids the non-termination of the recursion. This is a valid possibility,
because tabled types cannot currently involve type variables. It does however
point to difficulty that should be dealt with in the future: one often would
like recursive types (like search trees) to be seen as instances of user defined
type constructors, whose application to a concrete type would result in such a
recursive type; when that becomes possible, one needs to find a solution to
dealing with that here (ideally it would suffice to keep the recursive
constructor and continue the substitution into its arguments, like for other
kinds of types).

The case of a |variable_type| is the only one where something more interesting
happens; we shall deal with this in a separate section.

@< Function definitions @>=

type_expr substitution
  (const type_expr& tp, const type_assignment& assign, unsigned int shift_amount)
{ type_ptr result;
  switch (tp.raw_kind())
  { case primitive_type: return type_expr::primitive(tp.prim());
    case function_type: result =
      mk_function_type(substitution(tp.func()->arg_type,assign,shift_amount),
                       substitution(tp.func()->result_type,assign,shift_amount));
    break;
    case row_type:
      result = mk_row_type(substitution(tp.component_type(),assign,shift_amount));
    break;
    case tuple_type:
    case union_type:
    { dressed_type_list aux;
      for (wtl_const_iterator it(tp.tuple()); not it.at_end(); ++it)
        aux.push_back(substitution(*it,assign,shift_amount));
      return type_expr::tuple_or_union(tp.raw_kind(),aux.undress());
    }
    case tabled:
    { dressed_type_list aux;
      for (wtl_const_iterator it(tp.tabled_args()); not it.at_end(); ++it)
        aux.push_back(substitution(*it,assign,shift_amount));
      return type_expr::user_type(tp.tabled_nr(),aux.undress());
    }
    case variable_type:
      @< If a type is associated in |assign| to type variable
         |tp.typevar_count()|, return a copy of that type into which
         substitution was recursively applied; otherwise return a copy of the
         type variable itself @>
    default: assert(false); // there should be no undetermined type components
  }
  return std::move(*result);
}

@ The |equivalent| method looks up any substitution recorded for a type
variable, with a null pointer signalling a negative result. This makes
performing the substitution, if called for, quite easy. If a substitution is
done, we must go on applying other substitutions into the type expression found
by |equivalent|, but since we are no longer dealing with a subexpression
of~|tp|, the |shift_amount| argument is set to~|0| here. Polymorphic type
variables that are not being substituted for remain type variables, but we apply
|assign.renumber| to fill in the gaps left by type variables that are being
substituted for. This is possible since no context ever ascribes a meaning to
particular polymorphic variables (unlike fixed type variables, which should
never be renumbered), but it is probably not necessary: in places where
polymorphic type variables need to fill a contiguous range, notably in the
|type| class below, we shall apply a compacting renumbering of them anyway.

There is a potential for non-terminating recursion here, but that possibility is
avoided by ensuring |assign| has no direct or indirect self-references. To that
end we shall refuse, when setting a value for type variable in any
|type_assignment|, values (type expressions) that directly or indirectly refer
back to the variable in question itself.

@< If a type is associated in |assign|... @>=
{ auto c = tp.typevar_count()+shift_amount;
  auto p = assign.equivalent(c);
  return p==nullptr
    ? type_expr::variable(assign.renumber(c)) : substitution(*p,assign,0);
}

@ Next we define the unification procedure, which is called |can_unify| since
returns a Boolean value indicating whether unification succeeded. In case of
success the substitutions made to achieve unification are recorded in |assign|;
in case of failure, the caller should ignore (or wipe out) any substitutions
that might be recorded in |assign|. In case of success, the actual unified type
can then be obtained if desired by a subsequent call |substitute(P_orig,assign)|
or |substitute(Q_orig,assign)|, but |assign| can also be used in different
manners.

The treatment of type variables and tabled types is detailed in other sections.
Apart from that, we just perform a recursive traversal of both types, failing as
soon as different tags are found, and succeeding if the traversal completes
without that happening.

@< Function definitions @>=

bool can_unify
  (const type_expr& P_orig, const type_expr& Q_orig, type_assignment& assign)
{ const_type_p P=&P_orig, Q=&Q_orig; // we need assignable pointers internally
  auto P_kind = P->raw_kind(), Q_kind = Q->raw_kind();
@/@< If |P| or |Q| is a type variable with |typevar_count()>=assign.var_start()|,
     expand any existing substitution in |assign| for them, or when none is
     present record a new one and return whether that was possible;
     in expanded cases, |P|, |Q|, |P_kind| and |Q_kind| are updated @>
  assert(P_kind!=undetermined_type and Q_kind!=undetermined_type);
  // no patterns here
  if (P_kind==tabled or Q_kind==tabled)
    @< If both types are tabled, recursive, and different, |return false|;
       otherwise decide and |return| using recursive calls of |can_unify|:
       if both types are tabled and equal, call it for their type arguments,
       or else for their expansions @>
  if (P_kind!=Q_kind)
    return false;
    // with name expansions behind us, the type tags must match to succeed
  switch(P_kind)
  {
  case primitive_type: return P->prim()==Q->prim();
  case variable_type: // by above preparation, both are fixed in the context
    return P->typevar_count()==Q->typevar_count();
    // they are treated as primitive types
  case function_type: return
    can_unify(P->func()->arg_type,Q->func()->arg_type,assign) @|
    and
    can_unify(P->func()->result_type,Q->func()->result_type,assign);
  case row_type: return
    can_unify(P->component_type(),Q->component_type(),assign);
  case tuple_type: case union_type:
    {
      wtl_const_iterator it0(P->tuple()), it1(Q->tuple());
      while (not it0.at_end() and not it1.at_end()
             and can_unify(*it0,*it1,assign))
         ++it0,++it1;
      return it0.at_end() and it1.at_end(); // whether equal length lists
    }
  default: assert(false); // other cases were eliminated before the |switch|
    return false; // keep compiler happy
  }
}

@ The presence of recursive types creates the risk of an infinite recursion if
we simply expand type definitions as we usually do. To avoid that scenario, we
perform a test near the entry of the recursive function |can_unify|, but after
assigned type variables have been substituted. (The order is important here: an
assignment can equate a type variable to a tabled type, but a tabled type cannot
be defined as a type variable.) Types descended from tabled types are also
tabled, so it suffices to catch the case where both types are now simultaneously
tabled and recursive: any non-termination would have to run via such a case.

The code here is similar to the corresponding part of the equality test among
|type_expr| values. Since for non recursive tabled types we only know they are
not identical to another such type, but substitution instances of different
tabled type constructors (possibly of different arity) could be equal, the case
of finding distinct |tabled_nr()| values is only short-circuited if both types
are recursive; in the remaining cases we expand and try again. The recursive
call is necessary, as reassignment of |P| and |Q| (as we do in the case of type
variables) would run into lifetime problems for the expanded types. When we do
find equal |tabled_nr()| values, we proceed similarly to the handling of tuple
and union types we just saw, with argument types in the place of component
types; the only difference is that here we expect (and insist) that the argument
type lists (for the same tabled type constructor) are of equal length.

When we come here we know that at least one type is tabled; if indeed just one
is, we expand that one (whether recursive or not) and call |can_unify|
recursively.
@:avoiding infinite recursion@>

@< If both types are tabled, recursive, and different... @>=
{ if (P_kind==tabled and Q_kind==tabled)
  { if (P->tabled_nr()!=Q->tabled_nr())
    @/return not (P->is_recursive() and Q->is_recursive()) @|
        and can_unify(P->expanded(),Q->expanded(),assign);
    wtl_const_iterator it0(P->tabled_args()), it1(Q->tabled_args());
    while (not it0.at_end() and not it1.at_end())
      if (not can_unify(*it0,*it1,assign))
        return false;
      else
        {@; ++it0; ++it1; }
    assert (it0.at_end() and it1.at_end());
     // both length should match tabled arity
    return true;
  }
  if (P_kind==tabled)
    return can_unify(P->expanded(),*Q,assign);
  else
    return can_unify(*P,Q->expanded(),assign);
}

@ The code below is basically symmetric in |P| and |Q|, but we need to write out
the two cases anyway. When a type variable is encountered that already has an
equivalent in |assign|, we proceed with the substituted type expression in its
place, much like what is done with tabled type names but using their
|assign.equivalent| instead of they |definiens|. A difference is that while a
tabled type cannot refer directly to another one, |assign| can make one type
variable equal to another (through an assignment in one direction, to avoid
circularity); therefore we use |while| loops to iterate replacement until
something other than a type variable is found. Type variables numbered below
|assign.threshold| do not have an |equiv| entry in |assign|, so they will fall
through these loops, just like polymorphic variables without assignment to them.
When one of those polymorphic variables is present and the two are not already
equal, we call |assign.set_equivalent| to try to assign the other type
expression to it, and return whether that succeeds. The way |set_equivalent|
handles the case of a type expression that cannot be assigned because it refers
back to the same type variable, namely by returning |false|, means that this
just makes the unification fail; this is precisely what we want for this case.

@< If |P| or |Q| is a type variable... @>=
{ const_type_p p; // we first substitute already assigned type variables
  while (P_kind==variable_type and
         (p=assign.equivalent(P->typevar_count()))!=nullptr)
    P_kind=(P=p)->raw_kind(); // replace |P| by type previously assigned to it
  while (Q_kind==variable_type and
         (p=assign.equivalent(Q->typevar_count()))!=nullptr)
    Q_kind=(Q=p)->raw_kind(); // replace |Q| by type previously assigned to it
@)
  if (P_kind==variable_type and P->typevar_count()>=assign.var_start())
  { auto c = P->typevar_count();
    return (Q_kind==variable_type and Q->typevar_count()==c)
      or assign.set_equivalent(c,Q);
  }
  if (Q_kind==variable_type and Q->typevar_count()>=assign.var_start())
    return assign.set_equivalent(Q->typevar_count(),P);
}

@*2 Wrapped up polymorphic types.
%
For a long time, the recursive class |type_expr| was used both to represent a
pattern that the context expects for the type of an expression (for
instance \.{(*->*)} for the function part in a call expression), to be filled in
during type analysis, and the type ascribed to an expression after type
checking, or to a variable. With the advent of second order types, we have made
a clearer distinction between patterns and (possibly polymorphic) types, which
is already reflected in definition of |type_expr|. While type patterns typically
contain undetermined parts, types are completely determined but may contain
polymorphic type variables. Whenever a type found gets stored as the type of a
variable or function, we shall wrap the |type_expr| in an object of class |type|
defined below. This class notably provides data and methods to help administrate
polymorphic types. There is however no clear separation between the use of
|type_expr| and |type|, because when methods for polymorphic types are defined
by structural recursion, like the already introduced |substitution|, this is
easier to do using |type_expr|. Indeed, even the main type checking function
|convert_expr| returns the type it deduced in the form of a |type_expr|, which
the caller then often converts to a |type|. During the call of such recursive
functions operating upon |type_expr|, the extra information needed to build a
|type| from them is held elsewhere.

Polymorphic types arise in type checking when identifiers and functions with a
recorded polymorphic type are used. When combined in for instance function
calls, unification may also produce new polymorphic types. However the main way
in which fresh polymorphic types arise is by the use of ``type abstraction''
clauses, in which the user explicitly introduces some type variables using
the \&{any\_type} keyword. Within the scope of the clause introducing it, such a
type variable represents a fixed but unknown type; it is handled like a
primitive type, but with no related operations given for it. These type
variables may end up in the (usually function) type that is derived for the
clause, and when they thus emerge from the scope in which they were introduced,
they become polymorphic type variables. The latter are implicitly universally
quantified, and they can be substituted for during unification. Thus an
important quantity will be the threshold between type variables still in scope
(numbered below the threshold), and polymorphic type variables (numbered from
the threshold upwards). While manipulating |type_expr| values, this threshold is
stored under names like |fix_count| or |fc|, and it gets recorded with the
|type_expr| when wrapping it into a |type|. Our decision to make polymorphic
type variables ``float on top'' of the fixed ones often creates an annoying
requirement to renumber polymorphic type variables. But this probably would be
needed even if we had a separate pool of polymorphic variables, since
unification requires sets of polymorphic type variables from separate types to
be (made) disjoint. Therefore we stick to our decision here.

Polymorphic types are most often function types, and the main use of
polymorphism and unification occurs when we resolve function overloading.
However, unification (in a one-way sense) also occurs whenever a previously
determined polymorphic type arises in a context where a specific (monomorphic)
type is required, in which case the former must successfully unify to the
latter. Some expressions can be given a non-function polymorphic type, of which
the most common example is the empty row display, which in a context that does
not expect any particular component type is given the polymorphic type \.{[A]}
(similar examples arise when a polymorphic injector function for a union type
constructor is applied, and the variants of the union other than that of the
injector remain polymorphic).

We arrange things so that no context ever \emph{requires} a polymorphic type.
One cannot write a polymorphic type in a cast, and whenever an identifier is
given a polymorphic type at its initialisation, it is implicitly given the
constant attribute; this excludes the possibility later writing an assignment to
such a variable, in which the right hand side would be a context requiring value
of (sufficiently) polymorphic type to replace the initial value. As there are
usually not many different values with a given polymorphic type anyway, this is
not expected to be a severe restriction for user.

@ The main information added by a |type| to the |type_expr| it contains, is the
range of type variable numbers that are considered to be polymorphic. This is
contained in a |type_assignment| field, which provides its |var_start()| as the
lower bound of the range, and its |size()| (namely that of its |equiv| table)
for the size of the range. At the same time, that field allows (optionally and
temporarily) recording assignments to the polymorphic type variables, which
possibility is used in the unification process (when this has happened, the type
variables in question are considered to be neither fixed nor polymorphic). When
we don't need the type assignments or are done with them, we have the
possibility to clear them, and revert to the original (more) polymorphic type.

The class |type| has a large number of methods, many of them invoking some form
of unification or substitution of type variables. Here we give those related to
creating instances of |type|. Like for |type_expr| we have a move constructor
and move-assignment operator, while constructing a copy requires explicitly
calling |copy|, but the main ways of constructing |type| are by factory
functions (static methods). The principal one is |wrap| that converts a
|type_expr| to a |type|. It needs an indication |fix_nr| of the threshold above
which to consider variables polymorphic, and it will renumber those polymorphic
values in order of appearance in the type expression to a consecutive range of
numbers. An optional argument |gap| can be supplied which makes renumbering
start at |fix_count+gap| to ensure that |gap| numbers remain unused, presumably
to avoid clashes with some type variables elsewhere. We can alternatively use
the |constructor| method to get a |type| used to represent a type constructor.
In that case all type variables will be taken to be parameter types (since type
constructors cannot be introduced in the scope of a type abstraction), and the
caller passes the desired |degree| explicitly (there is no renumbering here, and
not all type parameters need to actually occur in the |type_expr|). Finally
|type::bottom| produces a neutral |type| consisting of just a type variable, the
pendant of an undetermined |type_expr|; it is used for instance as the starting
value in balancing of types, and thereby provides the component type for the
type of an empty row display.

@< Type definitions @>=
class type
{
  type_expr te;
  type_assignment a;
@)
  type() : te(), a(0,0) @+{}
  type(unsigned int fix_nr,unsigned int var_nr) : te(), a(fix_nr,var_nr) @+{}
public:
  static type wrap(const type_expr& te,
		   unsigned int fix_count, unsigned int gap=0);
  static type constructor(type_expr&& te, unsigned int degree);
  static type bottom(unsigned int fix_count); // a polymorphic type variable
  type(type&& tp) = default;
  type& operator=(type&& tp) = default;
  type copy() const @+
    {@; type result;
      result.te=te.copy();
      result.a=a.copy();
      return result;
    }
@)
  @< Methods of |type| to access component types @>@;
@)
  @< Utility methods of |type| @>@;
};

@ We access the component |type_expr| elements of a |type| using methods of the
same name of those of |type_expr|; we do not attempt to recreate |type| values
for them, so these methods have the same return types as their |type_expr|
counterparts. We include here also a few methods that access the
|type_assignment| field, and some simple tests.

@< Methods of |type| to access component types @>=
type_tag kind () const @+{@; return te.raw_kind(); }
type_tag top_kind () const @+{@; return te.top_kind(); }
primitive_tag prim () const     @+{@; return te.prim(); }
unsigned int typevar_count () const @+{@; return te.typevar_count(); }
const func_type* func() const  @+{@; return te.func(); }
      func_type* func()        @+{@; return te.func(); }
const type_expr& component_type () const @+{@; return te.component_type(); }
      type_expr& component_type ()       @+{@; return te.component_type(); }
const_raw_type_list tuple () const @+{@; return te.tuple(); }
      raw_type_list tuple ()       @+{@; return te.tuple(); }
type_nr_type tabled_nr () const @+{@; return te.tabled_nr(); }
@)
const type_assignment& assign () const @+{@; return a; }
unsigned int floor () const @+{@; return a.var_start(); }
unsigned int degree() const @+{@; return a.size(); }
unsigned int ceil() const @+{@; return floor()+degree(); }
  // start disjoint type variables here
bool is_polymorphic() const @+{@; return degree()>0; }
bool is_clean() const; // absence of pending type assignments
const type_expr& unwrap() const @+{@; return te ; }
bool is_void() const @+
{@; return te.raw_kind()==tuple_type and length(te.tuple())==0; }

@ The method |unwrap| above gives access to the stored |type_expr|, but ignores
any type assignments that were made. If one does want to take into account type
assignments, calling |expunge| will do so and then remove those type variables.
And if in addition one needs a |type_expr| rather than a |type|, one can instead
call |bake| or, if this is the final use of our |type| value, |bake_off|. To
forget any pending type assignments and revert to a previous state one can call
|clear|, possibly passing a degree which needs to be restored.

Unification with another type is done by the |unify| method. This can both
decrease the degree or increase it (by capturing type variables from the other
type); there is no fixed relationship between type variables before and after.
The |const| method |has_unifier| tests whether our type can unify to what the
|type_expr| expects. The method |unify_specialise| preforms one-sided
unification of our |type| (in |pattern| no type variables are changed, but
undetermined parts may be filled in)), recording the substitution required in
our |type_assignment|. The |matches| method is similar, but specific for use in
overload resolution, where our type is that found for the argument, and |formal|
is the argument type specification of one overloaded instance. In case of
success, our |assign()| can then be used to perform substitutions to the type of
that overloaded instance, while in case of failure the caller can call |clear|
before attempting another match. The |degree| parameter and output parameter
|shift_amount| are needed for polymorphic overloads: the former informs about
the degree of polymorphism of the overload, and the latter reports back the
amount by which the polymorphic variables needed to be shifted to avoid
collision with our type variables. Calls where no overloading is involved are
simpler: one just needs to attempt unification for the parameter part of the
function and the argument type; it can be done by calling |matches_argument| on
the the complete function type, and passing it the argument type.

@< Utility methods of |type| @>=
type& expunge(); // eliminate assigned type variables, by substitution
type& expand() {@; te.expand(); return *this; }
type_expr bake() const; // extract |type_expr| after substitution
type_expr bake_off(); // extract |type_expr|, sacrificing self if needed
type& clear(unsigned int d); // remove any type assignments, reserve |d| new ones
type& clear() @+{@; return clear(degree()); }
@)
bool unify(const type& other);
bool has_unifier(const type_expr& t) const;
bool unify_specialise(type_expr& pattern)
  // adapt to pattern while specialising pattern
  {@; return unify_specialise(te,pattern); }
  // recursive helper method does the work
bool unify_specialise(type_expr& pattern, unsigned int fix_count) const;
bool matches
  (const type_expr& formal, unsigned int poly_degree,
   unsigned int& shift_amount);
bool matches_argument(const type& arg_type); // for non-overloaded calls
@)
type_expr skeleton (const type_expr& sub_tp) const;
  // extract a type pattern from polymorphic |sub_tp|
void wrap_row () @+{@; te.set_from(type_expr::row(std::move(te))); }
  // put on a ``row-of''
private:
bool unify_specialise(const type_expr& sub_tp, type_expr& pattern);

@ First some simple methods. One uses |expunge| to incorporate any pending type
assignments, and |clear| to forget any and reset the number of type variables
to~|d| (presumably the number it had before some unification method might have
extended the set). Finally |bake| and |bake_off| provide the type expression
taking into account pending type assignments, the latter leaving the type itself
in an unusable state (since it was possibly moved from).

@< Function definitions @>=
bool type::is_clean() const @+{@; return a.empty(); }
type& type::expunge()
{
  if (is_clean())
    return *this; // nothing to expunge
  return *this = wrap(substitution(te,a),floor());
  // apply |a|, then renumber remaining type variables
}
@)
type& type::clear(unsigned int d)
@+{@;
  a = type_assignment(floor(),d);
  return *this;
}
@)
type_expr type::bake() const
@+{@; return substitution(te,a); }
type_expr type::bake_off() // variant available if |*this| is no longer needed
{
  if (a.empty())
    return std::move(te); // save some work here
  return substitution(te,a);
}

@ When turning a |type_expr| into a |type|, we build a transformed copy using a
local recursive function |pack|. While doing so it also eliminates
undetermined components by replacing them by fresh type variables; this probably
should never be necessary, but will ensure that in a |type| the |te| field never
|is_unstable()|. The implementation of |pack| is similar to that of
|substitution|, but returns a |type_expr| rather than a (smart) pointer to it.
The |translate| argument is a lookup list, mapping new (position) to old (value)
type variable numbers, which initially just maps frozen type variables to
themselves, and is updated with entries for type variables encountered or
created during the recursive traversal.

@< Local function definitions @>=
type_expr pack(const type_expr& te, sl_list<unsigned int>& translate)
{ switch (te.raw_kind())
  { case primitive_type: return type_expr::primitive(te.prim());
    case function_type:
    { auto arg_fix = pack(te.func()->arg_type,translate); // this one first
      auto res_fix = pack(te.func()->result_type,translate); // then this one
      return type_expr::function(std::move(arg_fix),std::move(res_fix));
    }
    case row_type: return type_expr::row(pack(te.component_type(),translate));
    case tuple_type:
    case union_type:
    { dressed_type_list aux;
      for (wtl_const_iterator it(te.tuple()); not it.at_end(); ++it)
        aux.push_back(pack(*it,translate));
      return  type_expr::tuple_or_union(te.raw_kind(),aux.undress());
    }
    case tabled:
    { dressed_type_list aux;
      for (wtl_const_iterator it(te.tabled_args()); not it.at_end(); ++it)
        aux.push_back(pack(*it,translate));
      return type_expr::user_type(te.tabled_nr(),aux.undress());
    }
    case undetermined_type:
    { unsigned int k=translate.size(); translate.push_back(-1);
      return type_expr::variable(k);
    }
    case variable_type:
    { unsigned int k=0;
      for (auto it=translate.begin(); not translate.at_end(it); ++it, ++k)
        if (*it==te.typevar_count()) // then we found a known type variable
          return type_expr::variable(k);
      translate.push_back(te.typevar_count()); // new variable; record it
      return type_expr::variable(k); // return new number for variable
    }
   default: assert(false); return type_expr();
  }
}


@ The call |type::wrap(t,fc,gap)| that converts a |type_expr t@;| to a |type|,
renumbering the type variables from $fc$ upwards into a consecutive range
starting at |fc+gap|. Using |pack| this is straightforward: we prepare a local
list |translate| to map type variables that should remain fixed themselves and
pad it out with |gap| dummy entries, and then let |pack| do the work. The
polymorphic degree of the result is determined by how much |pack| makes the
|translate| table grow, so we need to store the |type_expr| returned by |pack|
in a local variable |te| before constructing our |type|, and finally move |te|
into that type.

@< Function definitions @>=
type type::wrap (const type_expr& t, unsigned int fix_count, unsigned int gap)
{
  sl_list<unsigned int> translate;
  for (unsigned int k=0; k<fix_count; ++k)
    translate.push_back(k);
     // up to |fix_count|, type variables are unchanged,
  for (unsigned int k=0; k<gap; ++k)
    translate.push_back(-1); // reserve |gap| values as ``to remain unused''
@)
  type_expr te=pack(t,translate);
  fix_count += gap; // the new starting value
  type result(fix_count,translate.size()-fix_count);
  result.te = std::move(te);
  return result;
}

@ We can use a |type| to serve as body of a type constructor, in which case the
type variables it contains represent arguments of that constructor, identified
by position. We therefore do not want to do any renumbering, so the factory
function |type::constructor| does not call |pack|, sets the |fix_count| field
to~|0|, and sets the degree as requested by the caller, no questions asked.

The method |bottom| constructs a type whose body is a single polymorphic type
variable, so that it can unify to anything. It similarly avoids the
complications of |type::wrap|, using only the provided |fix_count| to choose a
variable number beyond any type abstractions in scope in the context of the
call.

@< Function definitions @>=
type type::constructor(type_expr&& te, unsigned int degree)
{ type result(0,degree);
  result.te = std::move(te);
  return result;
}
@)
type type::bottom(unsigned int fix_count)
{ type result(fix_count,1);
  result.te = type_expr::variable(fix_count);
  return result;
}

@ The methods |type::unify| and |type::has_unifier| are easily implemented using
|can_unify|.

@< Function definitions @>=
bool type::unify(const type& other)
{ assert(floor()==other.floor());
  const auto d = degree();
  if (other.is_polymorphic())
  {
    a.grow(other.degree()); // make place for other type variables
    if (can_unify(te,shift(other.te,floor(),d),a))
    {@;
      expunge();
      return true;
    }
    else
    {@;
      clear(d);
      return false;
    }
  }
  if (can_unify(te, other.te, a))
  {@;
    expunge();
    return true;
  }
  else
  {@;
    clear(d);
    return false;
  }
}
@)
bool type::has_unifier(const type_expr& t) const
{
  type tp = type::wrap(t,floor(),degree()); // renumber to avoid clashes
  tp.a = type_assignment(a.var_start(),degree()+tp.degree());
  // (ab)use |tp.a| as local variable
  return can_unify(te,tp.te,tp.a);
}

@ The method |type::unify_specialise| is like the function |can_unify|, but on
the side of the pattern the only changes are specialisations of undefined
subexpressions. Any type variables present in |pattern| are not substituted for:
unless they match up with the same type variable on our side, they cause
unification to fail. The necessary substitutions on our |type| side are recorded
in our |a| field, which is convenient for our implementation: any occurrence of
a type variable after the first will get the value that was substituted for it
the first time. Once the unification succeeds, the caller can decide whether to
preserve these type assignments for further unification, or use them to perform
substitutions, or forget them by calling |clear|. If the unification fails,
any type assignments should be forgotten by the caller.

@< Function definitions @>=

bool type::unify_specialise(const type_expr& sub_tp, type_expr& pattern)
{ auto P_kind = sub_tp.raw_kind(), Q_kind = pattern.raw_kind();
  if (P_kind==tabled or Q_kind==tabled)
    @< Decide |unify_specialise| in the presence of tabled types,
       avoiding any recursive calls if both type are tabled and recursive @>
  if (Q_kind==undetermined_type)
    {@; pattern.set_from(sub_tp.copy()); return true; }
  if (P_kind!=Q_kind and P_kind!=variable_type)
    return false;
  switch(P_kind)
  {
  case variable_type:
    { auto c = sub_tp.typevar_count();
      if (c<floor()) // fixed type; |Q| must match
        return Q_kind==variable_type and pattern.typevar_count()==c;
      auto eq=a.equivalent(c);
      if (eq!=nullptr)
        return unify_specialise(*eq,pattern);
      @< If |pattern| has |undetermined| entries, throw an error @>
      a.set_equivalent(c,std::make_unique<type_expr>(pattern.copy()));
      return true;
    }
  case primitive_type: return sub_tp.prim()==pattern.prim();
  case function_type: return
    unify_specialise(sub_tp.func()->arg_type,pattern.func()->arg_type) and @|
    unify_specialise(sub_tp.func()->result_type,pattern.func()->result_type);
  case row_type: return
    unify_specialise(sub_tp.component_type(),pattern.component_type());
  case tuple_type: case union_type:
    { wtl_const_iterator p_it(sub_tp.tuple()); // need two different types here
      for(wtl_iterator q_it(pattern.tuple());
          not (p_it.at_end() and q_it.at_end()); ++p_it,++q_it)
      { if (p_it.at_end() or q_it.at_end() or not unify_specialise(*p_it,*q_it))
          return false; // unequal lengths or some subtype fails unification
      }
      return true;
    }
  default: assert(false);
  // |tabled| impossible, and |undetermined_type| should not happen
  }
  return false; // keep compiler happy
}

@ This code is a bit subtle, and follows the pattern laid out in
section@#avoiding infinite recursion@>. Even though we are in a method of
|type| here, |*this| is only there to provide its |a| field, so recursive calls
can expand |sub_tp| and/or |pattern| as needed, and any assignments to |a| they
make will be picked up as they should. Since the |pattern| argument, expanding
it means assigning the expanded value to |pattern| before passing it down into
the recursive call.

@< Decide |unify_specialise| in the presence of tabled types,
   avoiding any recursive calls if both type are tabled and recursive @>=
{ if (P_kind==tabled and Q_kind==tabled)
  { if (sub_tp.tabled_nr()!=pattern.tabled_nr())
    @/return not (sub_tp.is_recursive() and pattern.is_recursive()) @|
        and unify_specialise(sub_tp.expanded(),pattern=pattern.expanded());
    wtl_const_iterator it0(sub_tp.tabled_args());
    wtl_iterator it1(pattern.tabled_args());
    while (not it0.at_end() and not it1.at_end())
      if (not unify_specialise(*it0,(*it1)=it1->expanded()))
        return false;
      else
        {@; ++it0; ++it1; }
    assert (it0.at_end() and it1.at_end());
     // both length should match tabled arity
    return true;
  }
  if (P_kind==tabled)
    return unify_specialise(sub_tp.expanded(),pattern);
  else
    return unify_specialise(sub_tp,pattern=pattern.expanded());
}

@ The scenario in which we are asked to unify a type variable with a pattern
that has undetermined parts without being completely undetermined is extremely
unlikely if at all possible. Rather than finding the ``right thing'' to do here,
like inventing new type variables to plug the undetermined parts, we prefer to
signal an error in such cases. It is awkward that we need to throw an error from
a method of |type|, with no expression to attach this type to, but a caller may
catch and repackage the error to provide more detail.

@< If |pattern| has |undetermined| entries, throw an error @>=
if (pattern.is_unstable())
{ std::ostringstream o;
  o << "Cannot unify a type variable and an incomplete type " << pattern;
  throw program_error(o.str());
}

@ The second method |type::unify_specialise| is basically the same as the first,
but is a |const| method so that in particular the |type_assignment| field |a| is
unchanged. This is achieved by copying our type before calling
|unify_specialise| on the copy. We use the occasion of copying to also renumber
our bound variables starting from a new level |lvl|, which is useful to make
place for variables in |pattern| that should be considered fixed there, and
disjoint from our bound variables.

@< Function definitions @>=

bool type::unify_specialise(type_expr& pattern, unsigned int lvl) const
{ assert(lvl>=floor()); // we cannot start numbering below our |floor()|
  return wrap(te,floor(),lvl-floor()).unify_specialise(pattern);
}

@ The method |matches| is typically called with as our type the type of an
(argument) expression, and as |f_par_tp| the parameter part of a function type
from the overload table. Both our (actual argument) type and the function type
can be polymorphic. Our type stores its own polymorphic |floor()| and
|degree()|, while for the function type the degree is passed as a separate
argument |f_deg|; coming from a global table, its type variables are all
polymorphic, starting from number~|0|. The task of this method is similar to
that of |f_par_tp.specialise|, but instead of filling undetermined slots, we are
deducing assignments to the free type variables in |f_par_tp|, which are then
(opportunistically) stored in the |type_assignment| field of |*this|. We assume
the caller has cleared all our previous type assignments, so we have a clean
slate of |degree()| type variables.

Since we are calling |can_unify|, we must first make the sets of type variables
disjoint, which we do by shifting any type variables of |f_par_tp| to start at
|ceil()|. The caller of |matches| should be aware that, when afterwards using
the |type_assignment| of the |type| object, the same shift should be applied to
any type expression related to |f_par_tp| substituted into; the optional final
argument of |substitution| can be used for this. We store the amount that our
method shifted by in its output parameter |shift_amount| for the convenience of
the caller (the value of |ceil()| from which it was copied will have been raised
after the call). Typically the above is used by the caller for performing
substitution into the result type of the function type that |f_par_tp| was taken
from.

@< Function definitions @>=

bool type::matches
  (const type_expr& f_par_tp, unsigned int f_deg, unsigned int& shift_amount)
{
  shift_amount = ceil(); // record where our type assignments used to end
  a.grow(f_deg); // create space for new type variables
  if (shift_amount==0 or f_deg==0) // then no need to renumber |f_par_tp|
    return can_unify(f_par_tp,te,a);
  return can_unify(shift(f_par_tp,0,shift_amount),te,a);
}

@ The method |matches_argument| is a variation on |matches| to be used in
function calls in which the function does not come from the overload table. Here
the preconditions are different: both function and argument have a |type| (which
was not the case for a function type in the overload table), which can be
polymorphic and in addition have fixed type variables; their |floor()| values
that separate the two regimes are the same. We choose |matches_argument| to be a
method called for the full function type, passing the actual argument type to it
as (constant) argument (this is the opposite order from what |matches| does); in
case of success, the substitution is recorded in the |type_assignment| of the
function type. As in the case of |matches|, we need to make the polymorphic
variable sets disjoint, so if both sets are non empty, we shift the argument
polymorphic variables to follow those of the function type. In this manner the
type assignment can be applied, after a successful match, to the result type
without any shift.

@< Function definitions @>=

bool type::matches_argument(const type& actual_arg_type)
{
  assert(floor()==actual_arg_type.floor());
  unsigned int fd=degree(), ad=actual_arg_type.degree();
  a.grow(ad);
  const type_expr& arg_type = expand().func()->arg_type;
  if (fd==0 or ad==0) // then no need to renumber |actual_arg_type|
    return can_unify(arg_type,actual_arg_type.te,a);
  return can_unify(arg_type,shift(actual_arg_type.te,floor(),fd),a);
}

@ Before converting the argument for a polymorphic function, we can extract from
its polymorphic type a type pattern that may aid the conversion, namely by
simply replacing all type variables bound in the polymorphic type by
undetermined types; the type variables up to |floor()| are unchanged. The method
|skeleton| achieves this; it could have been implemented by calling
|substitution| with an assignment of undetermined types, but it can be done just
as easily by a direct recursion.

@< Function definitions @>=
type_expr type::skeleton (const type_expr& sub_t) const
{ type_ptr result;
  switch (sub_t.raw_kind())
  { case primitive_type: return type_expr::primitive(sub_t.prim());
    case function_type: result =
      mk_function_type(skeleton(sub_t.func()->arg_type),
                       skeleton(sub_t.func()->result_type));
    break;
    case row_type:
      result = mk_row_type(skeleton(sub_t.component_type()));
    break;
    case tuple_type:
    case union_type:
    { dressed_type_list aux;
      for (wtl_const_iterator it(sub_t.tuple()); not it.at_end(); ++it)
        aux.push_back(skeleton(*it));
      return type_expr::tuple_or_union(sub_t.raw_kind(),aux.undress());
    }
    case tabled:
    { dressed_type_list aux;
      for (wtl_const_iterator it(sub_t.tabled_args()); not it.at_end(); ++it)
        aux.push_back(skeleton(*it));
      return type_expr::user_type(sub_t.tabled_nr(),aux.undress());
    }
    case variable_type:
    { auto c = sub_t.typevar_count();
      return c>=a.var_start() ? type_expr() : type_expr::variable(c);
    }
    default: assert(false);
  }
  return std::move(*result);
}


@ We also provide an overload for printing |type| values without having to call
|expr| explicitly each time.

@< Function definitions @>=
std::ostream& operator<<(std::ostream& strm, const type& t)
{@; return strm << t.bake(); }

@ We need to export the declaration of that last output operator instance.

@< Declarations of exported functions @>=
std::ostream& operator<<(std::ostream& strm, const type& t);

@*1 Specifying types by strings.
%
In practice we shall rarely call functions like |mk_prim_type| and |mk_row_type|
directly to make explicit types, since this is rather laborious. Instead, such
explicit types will be constructed by the function |mk_type_expr| that parses a
(\Cee~type) string |s|, and correspondingly calls the appropriate type
constructing functions to build its return value. Now that we extended this
function, allowing it to scan and return polymorphic types, we add an output
parameter |var_count| that signals how many distinct type variables were found.

In some cases, as in section @# first types section @>, we are interested in
type patterns that allow `\.*' to be used to designate an |undetermined_type|,
but which are not intended to be in other ways polymorphic. For those occasions
we use |mk_type_pattern|, which lacks the |var_count| parameter. It can also be
used when an explicit non polymorphic type string is given, to avoid the need to
supply a |var_count| variable that serves no purpose. Finally |mk_type| scans a
type string and produces a (wrapped) |type| value.

@< Declarations of exported functions @>=
type_expr mk_type_expr(const char* s,unsigned int& var_count);
type_expr mk_type_pattern(const char* s);
type mk_type(const char* s);

@ The task of converting a properly formatted string into a type is one of
parsing a simple kind of expressions. The strings used here come from string
denotations in the source code (mostly in calls installing built-in functions
into \axis.) rather than from user input, and we are not going to write
incorrect strings (we hope). Therefore we don't care if the error handling is
crude here. The simplest way of parsing ``by hand'' is recursive descent, so
that is what we shall use. By passing a character pointer by reference, we
allow the recursive calls to advance the index within the string read.

@< Local function definitions @>=
type_expr scan_type(const char*& s, std::string& vars);
type_expr scan_in_parens(const char*& s, std::string& vars);
type_expr scan_union_list(const char*& s, std::string& vars);
type_expr scan_tuple_list(const char*& s, std::string& vars);

@ The function |scan_type| does the real parsing, |mk_type_expr| calls it,
providing a local modifiable pointer to bind to its reference parameter (which
is important because |scan_type| cannot directly accept a \Cee-string constant
as argument) while also doing error reporting. The function |mk_type_expr| is
called only during the start-up phase of \.{atlas} (but after the table
|prim_names| of primitive type names is installed), and if an error is
encountered (of type |logic_error|, since this must be an error in the \.{atlas}
program itself), printing of the error message will be followed by termination
of the program. Although the function |mk_type_pattern| should not produce
polymorphic types and |mk_type_expr| should produce a type without undefined
type subexpressions, it not really worth the hassle to ensure this by writing
two separate implementations, so we implement |mk_type_pattern| by calling
|mk_type_expr| and ignoring the number of type variables it finds.

@< Function definitions @>=
type_expr mk_type_expr(const char* s,unsigned int& var_count)
{ const char* orig=s; std::string vars;
  try
  { auto result=scan_type(s,vars);
    var_count = vars.length();
    return result;
  }
  catch (logic_error& e)
  { std::cerr << e.what() << "; original string: '" << orig @|
              << "' text remaining: '" << s << "'\n";
    throw;
  // make the error hard to ignore; if thrown probably aborts the program
  }
}
@)
type_expr mk_type_pattern(const char* s)
{@; unsigned int dummy; return mk_type_expr(s,dummy);
}
@)
type mk_type(const char* s)
{ unsigned int dummy;
  auto result = type::wrap(mk_type_expr(s,dummy),0);
  assert (result.degree()==dummy);
  return result;
}

@ The part of |scan_type| dealing with enclosed or atomic type expressions is
quite simple.

@< Local function definitions @>=
type_expr scan_type(const char*& s, std::string& vars)
{ if (*s=='(')
  { type_expr result=scan_in_parens(++s,vars);
    if (*s++!=')')
      throw logic_error("Missing ')' in type");
    return result;
  }
  else if (*s=='[')
  {
    type_expr t = scan_in_parens(++s,vars);
    if (*s++!=']')
      throw logic_error("Missing ']' in type");
    return type_expr::row(std::move(t));
  }
  else if (*s=='*') return ++s,type_expr(); // undetermined type
  else @< Scan and |return| a primitive type, or |throw| a |logic_error| @>
}
@)


@ For primitive types we use the same strings as for printing them. Since none
of them consists of a single character, and single letter names will suffice as
type variables for built-in functions; we only retainin |vars| the collection of
type variables seen within the current expression and the position of each type
variable in it.

In this module we use the fact that the order in the list |prim_names| matches
that in the enumeration type |primitive_tag|, by casting the integer index
into the former list to an element of that enumeration.

@h <cctype> // |isalpha|
@< Scan and |return| a primitive type, or |throw| a |logic_error| @>=
{ std::string str;
  while (isalpha(*s))
    str.push_back(*s++);
  auto l = str.length();
  if (l==1)
  { auto pos = vars.find(str[0]);
    if (pos==vars.npos)
    {@; pos=vars.length();
      vars+=str[0];
    }
    return type_expr::variable(pos);
  }
  if (l>1)
  { for (size_t i=0; i<nr_of_primitive_types; ++i)
      if (str==prim_names[i])
        return type_expr::primitive(static_cast<primitive_tag>(i));
    std::cerr << str << ": ";
    throw logic_error("Primitive type unrecognised");
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

@h <string>

@< Local function definitions @>=
type_expr scan_in_parens(const char*& s, std::string& vars)
{ type_expr a=scan_union_list(s,vars);
  if (*s!='-' or s[1]!='>')
    return a;
  return type_expr::function(std::move(a),scan_union_list(s+=2,vars));
}
@)
type_expr scan_union_list(const char*& s, std::string& vars)
{ dressed_type_list variants;
  while (variants.emplace_back(scan_tuple_list(s,vars)),*s=='|')
    ++s;
  return variants.size()==1
    ? std::move(variants.front())
    : type_expr::tuple_or_union(union_type,variants.undress());
}
@)
type_expr scan_tuple_list(const char*& s, std::string& vars)
{ static const std::string term("|-)]");
  dressed_type_list members;
  if (term.find(*s)==std::string::npos)
    // only act on non-terminating characters
    while (members.emplace_back(scan_type(s,vars)),*s==',')
      ++s;
  return members.size()==1
    ? std::move(members.front())
    : type_expr::tuple(members.undress());
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
const type_expr void_type = type_expr::tuple(empty_tuple());
const type_expr int_type = type_expr::primitive(integral_type);
const type_expr bool_type = type_expr::primitive(boolean_type);
const type_expr row_of_type(mk_type_pattern("[*]"));
const type_expr gen_func_type(mk_type_pattern("(*->*)"));

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
extern const type_expr KType_type; // \.{KType}
extern const type_expr KTypePol_type; // \.{KTypePol}
extern const type_expr param_type; // \.{Param}
extern const type_expr param_pol_type; // \.{ParamPol}

@ The definitions below have all become self-contained, due to the use of
|mk_type_expr| for non-primitive types. Indeed, since we cannot have sharing
between type (sub-)expressions, the economy of using the |copy| method for
previously constructed type constants would be truly marginal. So in their
current form, some of these definitions could now (again) be moved to other
compilation units, where they might even be just local constants.

@: second types section @>

@< Global variable definitions @>=
const type_expr rat_type = type_expr::primitive(rational_type);
const type_expr str_type = type_expr::primitive(string_type);
const type_expr vec_type = type_expr::primitive(vector_type);
const type_expr ratvec_type = type_expr::primitive(rational_vector_type);
const type_expr mat_type = type_expr::primitive(matrix_type);
const type_expr row_of_int_type(mk_type_pattern("[int]"));
const type_expr row_of_rat_type(mk_type_pattern("[rat]"));
const type_expr row_of_vec_type(mk_type_pattern("[vec]"));
const type_expr row_of_ratvec_type(mk_type_pattern("[ratvec]"));
const type_expr row_row_of_int_type(mk_type_pattern("[[int]]"));
const type_expr row_row_of_rat_type(mk_type_pattern("[[rat]]"));
const type_expr pair_type(mk_type_pattern("(*,*)"));
const type_expr int_int_type(mk_type_pattern("(int,int)"));
const type_expr Lie_type_type = type_expr::primitive(complex_lie_type_type);
const type_expr rd_type = type_expr::primitive(root_datum_type);
const type_expr ic_type = type_expr::primitive(inner_class_type);
const type_expr rf_type = type_expr::primitive(real_form_type);
const type_expr split_type = type_expr::primitive(split_integer_type);
const type_expr KType_type = type_expr::primitive(K_type_type);
const type_expr KTypePol_type = type_expr::primitive(K_type_pol_type);
const type_expr param_type = type_expr::primitive(module_parameter_type);
const type_expr param_pol_type = type_expr::primitive(virtual_module_type);

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
  return type_expr::tuple(std::move(tl));
}

@* Run-time values.
%
Now we shall consider run-time values. As we mentioned before, the interpreter
must access values via generic pointers in order to be able to manipulate them
regardless of their types, which could be arbitrarily complicated. We could
either use void pointers to represent generic values and cast them when
necessary, or use inheritance and the dynamic cast feature of \Cpp. We choose
the second option, which is quite convenient to use, although this means that in
reality we have dynamic type information stored within the values, even though
that information had already been determined during type analysis. We shall in
fact use this information to double-check our type analysis at run time.

@< Type declarations @>=
struct value_base;
using value = const value_base*;
using shared_value = std::shared_ptr<const value_base>;
using  own_value = std::shared_ptr<value_base>;

@~We start with a base class for values. For it to be an abstract class, there
must be at least one pure virtual function in the class; the destructor having
to be virtual anyway, we make it pure virtual (this does not mean it is
unimplemented, in fact it must be implemented as it will always be called,
after the destructor for a derived class; it just means derived
classes \emph{must} override the default). The printing function does not have
a useful default, so we make it pure virtual as well, without providing an
implementation (in the base class). This |print| method will demonstrate the
ease of using dynamic typing via inheritance; it will not do any dynamic
casting, but other operations on values will.

@< Includes needed in \.{axis-types.h} @>=
#include <iostream> // needed for specification of |print| method below

@~Apart from virtual methods, we define other methods that will be redefined in
all or some derived classes, and which are selected based on the static type at
hand in the calling code rather than on the dynamic type, as the former will be
known to match the latter due to our type system. The method |name| is used in
reporting logic errors from function templates, notably the failure of a value
to be of the predicted type, where |name| names that predicted type. Since
callers are functions that know the type (via a template argument) but need not
have any object of that type at hand, we define |name| to be a |static|. For
certain derived types there will be operations implemented by making changes to
an existing value; we use copy-on-write when our reference to the initial value
is shared, and these derived classes provide a (usually default) copy
constructor, which will be invoked by the |get_own| function template below. It
used to call a virtual |clone| method, but that forces \emph{all} derived
classes to implement copying, while in the current situation those that do not
use |get_own| may simply not implement copying by deleting the copy constructor.

We disable assignment of |value_base| objects, since they
should always be handled by reference; the base class is abstract anyway, but
this ensures us that for no derived class an implicitly defined assignment
operator is accidentally invoked. Though derived copy constructors will be
|private|, they need to copy the base object, so we provide a |protected| copy
constructor.

Values are always handled via pointers. The raw pointer type is |value|, and a
shared smart pointer-to-constant is |shared_value|. The const-ness of the
latter reflects a copy-on-write policy: we rarely need to modify values
in-place, but when we do, we ensure our shared pointer is actually unique, and
then |const_cast| it to a derived version of |own_value| for modification (this
will be hidden in a function template defined later).

@< Type definitions @>=
struct value_base
{ value_base() @+ {}
@/virtual ~value_base() = 0;
  virtual void print(std::ostream& out) const =0;
@)
// |static const char* name();| just a model; defined in derived classes
@)
  value_base& operator=(const value_base& x) = delete;
};
inline value_base::~value_base() @+{} // necessary but empty implementation

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
we must pass through a level of indirection, so we have a vector of pointers. We
define these pointers to be |shared_value| pointers, so that the row takes
(shared) ownership of its components without needing an explicit destructor.
This has the additional advantage over explicit ownership management that the
copy constructor can safely just copy-construct the vector of pointers: a
possible exception thrown during the copy is guaranteed to clean up any pointers
present in the vector (resetting their reference counts to their original
values). Note also that default-constructed shared pointers are set to null
pointers, so the first constructor below, which just reserves space for |n|
shared pointers, has set them to exception-safe values while waiting for the
slots to be filled.

Of course ownership of pointers to |row_value| objects also needs to be managed.
The type |own_row| will be used after constructing the row while filling in the
contents, or after ensuring unique ownership of the row in order to perform
destructive operations (so although a |shared_ptr<value_base>|, it is known to
actually be unique); at all other times the pointer converted to |shared_row|
(and possibly down-cast to |shared_value|) will be used.

@< Type definitions @>=
struct row_value : public value_base
{ std::vector<shared_value> val;
@)
  explicit row_value(size_t n) : val(n) @+{} // start with |n| null pointers
  template <typename I> row_value(I begin, I end) : val(begin,end)
    @+{} // set from iterator range
  void print(std::ostream& out) const;
  size_t length() const @+{@; return val.size(); }
  static const char* name() @+{@; return "row value"; }
  row_value @[(const row_value& ) = default@];
    // we use |get_own<row_value>|
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
  void print(std::ostream& out) const;
  static const char* name() @+{@; return "tuple value"; }
@/tuple_value (const tuple_value& ) = default;
    // we use |uniquify<tuple_value>|
};
@)
using tuple_ptr = std::unique_ptr<tuple_value>;
using shared_tuple = std::shared_ptr<const tuple_value>;
using own_tuple = std::shared_ptr<tuple_value>;

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
  void print(std::ostream& out) const;
  static const char* name() @+{@; return "union value"; }
  union_value (const union_value& v) = delete;
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

@< Includes needed in \.{axis-types.h} @>=
#include <iterator>

@~Frames are actually allocated on the heap, and their lifetimes do not follow
a stack regime unless a very limited use is made of user-defined functions
(never passing such a function as value out of the expression in which it was
defined), so it is better to just say they are linked lists of frames. A
singly linked list suffices, and by using shared pointers as links,
destruction of frames once inaccessible is automatic.

@s back_insert_iterator vector

@< Type definitions @>=
using shared_context = std::shared_ptr<class evaluation_context>;
class evaluation_context
{ shared_context next;
  std::vector<shared_value> frame;
  evaluation_context(const evaluation_context&) = delete;
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

Most of the time, executable subexpressions will not be shared in any way, so
they will be passed around and linked together via unique-pointer values. We
choose to make |expression_ptr| a pointer-to-constant type, because evaluating
an expression (in some provided context) will never change that expression.
During the building of the expression tree, it will sometimes happen that we
want to modify it after the fact to achieve some kind of optimisation, and in
such cases it will happen that we need to const-cast away the |const| that is
introduced here (usually after also having applied a |dynamic_cast| to a derived
type). When handling user defined functions, we shall have values that refer to
(derived from) |expression| objects, and in doing so share them. So in those
cases, |shared_expression| values will be used.

@< Type declarations @>=
struct expr; // abstract syntax tree representation, see \.{parsetree.w}
struct expression_base; // executable expression
enum class eval_level : unsigned;
using expression = expression_base*;
using expression_ptr = std::unique_ptr<const expression_base>;
using shared_expression = std::shared_ptr<const expression_base>;

@ A fundamental choice is whether to make the result type of the |evaluate| type
equal to |value|. Although this would seem the natural choice, we prefer
to make its result |void|, and to handle all value-passing via an execution
stack; thus we hope to be able to avoid packing and unpacking of tuple values
in most cases when calling functions. In order to do that effectively, we
dynamically pass a parameter to the |evaluate| method telling whether the
result is expected to be ``expanded'' on the runtime stack in case it is of a
tuple type.

@< Type definitions @>=
enum class eval_level : unsigned @+{ no_value, single_value, multi_value };
struct expression_base
{ using level = eval_level;
@)
  expression_base() @+ {}
  expression_base(const expression_base&) = delete; // they are never copied
  expression_base(expression_base&&) = delete; // nor moved
  expression_base& operator=(const expression_base&)=delete; // nor assigned
  expression_base& operator=(expression_base&&)=delete; // nor move-assigned
  virtual ~expression_base() @+ {}
@)// other virtual methods
  virtual void evaluate(level l) const =0;
  virtual void print(std::ostream& out) const =0;
@)// non-virtuals that call the virtual |evaluate|
  void void_eval() const @+{@; evaluate(eval_level::no_value); }
  void eval() const @+{@; evaluate(eval_level::single_value); }
  void multi_eval() const @+{@; evaluate(eval_level::multi_value); }
};

@ Like for values, we can assure right away that printing converted
expressions will work.

@< Declarations of exported functions @>=
inline std::ostream& operator<< (std::ostream& out, const expression_base& e)
{@; e.print(out); return out; }

@*1 The execution stack. Our |evaluate| methods will put values on the execution
stack, which we shall now declare. We decide that values on the execution stack
can be shared with other values (for instance when the user subscripts a row,
vector or matrix bound to an identifier, it would be wasteful to duplicate that
entire structure just so that it can briefly reside on the execution stack),
whence we use |shared_value| smart pointers in the stack.

This choice will have consequences in many places in the evaluator, since once a
value is referred to by such a smart pointer, its ownership cannot be
transferred to any other regime; when strict ownership should be needed, the
only option would be to make a copy. However, it turns out to be convenient
to \emph{always} use shared pointers for runtime values, and to make the
distinction concerning whether one knows this pointer to be unique by having its
type be pointer-to-non-const in that case. After construction or duplication one
can start out with such a pointer, use it to store the proper value, then
convert it pointer-to-const (i.e., |shared_value|), for handing on the stack and
passing around in general; if destructive access is needed one may reconvert to
pointer-to-non-const after having checked unique ownership (or else having
duplicated the value pointed to).

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
pushed onto the stack, but only if the |level l@;| so indicates and the value is
indeed of tuple type; the function |push_expanded| will help doing this. Since
the argument might well be a shared pointer that was just created by
|std::make_shared|, we provide an rvalue version that will avoid changing any
reference count in such cases; it also provides the caller with an opportunity
to explicitly let it give up its ``share'' of an existing shared pointer in
passing it to |push_expanded|, by invoking |std::move| on the argument.

@< Declarations of exported functions @>=
void push_expanded(eval_level l, const shared_value& v);
void push_expanded(eval_level l, shared_value&& v);

@~Type information is not retained in compiled expression values, so
|push_expanded| cannot know which type had been found for |v| (moreover,
|push_expanded| is typically called for arguments for which the type is not
determined at \.{atlas} compile time). But it can use a dynamic
cast do determine whether |v| actually is a tuple value or not.

@< Function definitions @>=
void push_expanded(eval_level l, const shared_value& v)
{ if (l==eval_level::single_value)
    push_value(v);
  else if (l==eval_level::multi_value)
  { shared_tuple p = std::dynamic_pointer_cast<const tuple_value>(v);
    if (p==nullptr)
      push_value(v);
    else
      for (size_t i=0; i<p->length(); ++i)
        push_value(p->val[i]); // push components, copying shared pointers
  }
} // if |l==eval_level::no_value| then do nothing
@)
void push_expanded(eval_level l, shared_value&& v)
{ if (l==eval_level::single_value)
    push_value(std::move(v));
  else if (l==eval_level::multi_value)
  { shared_tuple p = std::dynamic_pointer_cast<const tuple_value>(v);
    if (p==nullptr)
      push_value(std::move(v));
    else if (v=nullptr,p.unique())
      // if caller held unique copy of pointer, we may dismember the tuple
      for (size_t i=0; i<p->length(); ++i)
        push_value(std::move(p->val[i]));
          // push components, moving shared pointers
    else // others than caller might hold a copy
      for (size_t i=0; i<p->length(); ++i)
        push_value(p->val[i]); // push components, copying shared pointers
  }
} // if |l==eval_level::no_value| then do nothing

@*1 Some useful function templates.
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
type it is known to have (by the type check that was passed); should such a cast
fail, this reveals a flaw of our type system, so we throw a |logic_error|. The
function template |get| with explicitly provided type serves for this purpose;
it returns a shared pointer to constant (because values on the stack are such
pointers) of the proper kind. Since the dynamic cast should never fail in
practice (and there is very extensive experience that it does not, although this
of course not exclude the possibility that it will in very rare circumstances,
or after further development of the program), we disable the test when compiling
without debugging, by doing a static cast rather than a dynamic cast; this just
pops and returns the pointer without inspecting it in any way.

Sometimes defeating copy-on-write is desired (to allow changes like filling
internal tables that will \emph{benefit} other shareholders), and
|non_const_get| will const-cast the result of |get| to allow that.

@< Template and inline function definitions @>=

template <typename D> // |D| is a type derived from |value_base|
 inline std::shared_ptr<const D> get()
{
#ifdef NDEBUG
  return std::static_pointer_cast<const D>(pop_value());
#else
  std::shared_ptr<const D> p=std::dynamic_pointer_cast<const D>(pop_value());
  if (p.get()==nullptr)
  { std::ostringstream o; o << "Argument is no " << D::name();
    throw logic_error(o.str());
  }
  return p;
#endif
}
@.Argument is no ...@>

@)
template <typename D> // |D| is a type derived from |value_base|
  inline std::shared_ptr<D> non_const_get()
{@; return std::const_pointer_cast<D>(get<D>()); }

@ Here is a function template similar to |get|, that applies in situations where
the value whose type is known does not reside on the stack. As for |get| we
convert using a dynamic case, and to throw a |logic_error| in case our type
prediction was wrong. This function is defined at the level of ordinary
pointers, and it is not intended for use where the caller assumes ownership of
the result; the original pointer is assumed to retain ownership as long as the
result of this call survives, and in particular that pointer should probably not
be obtained by calling the |get| method for a smart pointer temporary, nor
should the result of |force| converted to a smart pointer, lest double deletion
would ensue. As a consequence, we here use the basic |dynamic_cast| of a raw
pointer rather than a |dynamic_pointer_cast| of a |std::shared_ptr|. Like in the
case of |get| we provide a version using a static cast (omitting any check) when
no debugging is enabled, since the validity of the downcast should be ensured by
haveing passed the type check.

We provide two versions, where overloading will choose one or the other
depending on the const-ness of the argument. Since calling |get| for a
|shared_value| pointer returns a |value| (which is pointer to constant), it
will often be the second one that is selected.

@< Template and inline function definitions @>=
template <typename D> // |D| is a type derived from |value_base|
  D* force (value_base* v)
{
#ifdef NDEBUG
  return static_cast<D*>(v);
#else
  D* p=dynamic_cast<D*>(v);
  if (p==nullptr)
  { std::ostringstream o; o << "forced value is no " << D::name();
    throw logic_error(o.str());
  }
  return p;
#endif
}
@)
template <typename D> // |D| is a type derived from |value_base|
  const D* force (value v)
{
#ifdef NDEBUG
  return static_cast<const D*>(v);
#else
  const D* p=dynamic_cast<const D*>(v);
  if (p==nullptr)
  { std::ostringstream o; o << "forced value is no " << D::name();
    throw logic_error(o.str());
  }
  return p;
#endif
}

@ The \axis. language allows assignment operations to components of aggregates
(such as rows, matrices, strings, tuples), which in fact assign a new value to
the name bound to the aggregate. Since we implement copy-on-write, this should
in principle make a copy of the aggregate before modifying the component, but
the copy can be avoided if the aggregate value is unshared. The operation
|uniquify| implements making a copy if necessary, and calling it makes clear our
destructive intentions. If copying is necessary we invoke the copy constructor
via |std::make_shared|, and |uniquify| is templated by the actual derived type
to specify which copy constructor to use. This also allows us to export a
pointer to that derived type. We give |uniquify| a modifiable lvalue
|shared_value| (shared pointer-to-const) argument, to which a new shared (but
currently unique) pointer will be assigned in case duplication was necessary,
and it returns a raw pointer-to-non-const, which can then be used to make the
change to the unique copy. The argument |v| retains ownership. Not surprisingly
the implementation of |uniquify| uses a |const_cast| operation when no
duplication takes place.

The method |std::shared_ptr::unique| used here was removed from recent versions
of the \Cpp-standard, because it does not play well in multi-threaded
environments where some other thread might either still be working on a method
invoked using a since destroyed copy of the pointer, or resuscitate a copy form
a still existing |std::weak_ptr|. Our interpreter is not (yet) capable of
running simultaneously in multiple threads (the only multi-threading currently
used occurs withing a single built-in function), so this is no concern to us.
Should we for some other reason need to move to a recent version of \Cpp, then
we must make a home grown variant of |std::shared_ptr| that provides an
attribute that makes |unique| possible again, which is set immediately after
|make_shared| and irreversibly cleared as soon as any copy is made.

@< Template and inline function def... @>=
template <typename D> // |D| is a type derived from |value_base|
  D* uniquify(shared_value& v)
{ const D* p = force<D>(v.get());
  if (v.unique())
    return const_cast<D*>(p); // we can now safely write to |*p|
  auto result = std::make_shared<D>(*p);
    // invokes copy constructor; assumes it exists for |D|
  v = result; // upcast and constify (temporarily duplicates shared pointer)
  return result.get();
    // now |v| retains unique shared pointer, but we can modify its target
}

@ Similarly some wrapper functions will want to get unique access to an argument
object, so that they can return a modified version of it as result. This
requires taking an argument pointer from the stack, and duplicating the argument
value if the pointer has |use_count()>1| (a rare circumstance, since most
functions that place a value on the stack will do so with fresh pointer from
|std::make_shared|). We provide a variant function template of |get| called
|get_own|, that operates like |uniquify| except that ownership is not retained
by an external object, but by the result from |get_own| which therefore is a
shared pointer (but guaranteed to be unique at return). It uses a
|std::const_pointer_cast| in the case where no copy is made.

Finally there is |force_own|, which is similar to |get_own| but takes its value
not from the stack but from some other place that will not retain ownership;
like for |get_own| ownership passes to the result, which is a shared (but
unique) pointer to non-|const| for the derived type. In order for this to work,
the original copy of the pointer must be cleared at the time |unique| is called,
so that this call has some chance of returning |true|. To the end we take the
argument as rvalue reference, and make sure a temporary is move-constructed from
it inside the body of |force_own|; the temporary is constructed in the argument
to |std::dynamic_pointer_cast| (which, until \Cpp20, has no overloads that
directly bind to, and upon success move from, and rvalue argument; our
work-around moves from the rvalue even if the dynamic cast fails, but then we
throw a |logic_error| anyway).

Since these functions return pointers that are guaranteed to be unique, one
might wonder why no use of |std::unique_ptr| is made. The answer is this is
not possible, since there is no way to persuade a |shared_ptr| to release its
ownership (as in the |release| method of unique pointers), even if it happens
to be (or is known to be) the unique owner.

@< Template and inline function def... @>=
template <typename D> // |D| is a type derived from |value_base|
  std::shared_ptr<D> get_own()
{ std::shared_ptr<const D> p=get<D>();
  if (p.unique())
    return std::const_pointer_cast<D>(p);
  return std::make_shared<D>(*p);
    // invokes copy constructor; assumes it exists for |D|
}
@)
template <typename D> // |D| is a type derived from |value_base|
  std::shared_ptr<D> force_own(shared_value&& q)
{
  std::shared_ptr<const D> p=
#ifdef NDEBUG
     std::static_pointer_cast<const D>(shared_value(std::move(q)));
#else
     std::dynamic_pointer_cast<const D>(shared_value(std::move(q)));
  if (p==nullptr)
  { std::ostringstream o; o << "forced value is no " << D::name();
    throw logic_error(o.str());
  }
#endif
  if (p.unique())
    return std::const_pointer_cast<D>(p);
  return std::make_shared<D>(*p); // invokes copy constructor; assumes it exists
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
@*1 Built-in functions.
%
Ultimately, the evaluation process leads to the evaluation of built-in
operations or functions, which will be accessed though function pointers called
wrapper functions; these will call library functions after checking
preconditions that the latter assume without checking them on each call. The
keeps the library independent of the interpreter, while providing the user of
the interpreter with a level of safety from accidental crashes that is not (for
efficiency reasons) provided systematically in the library.

@< Type declarations @>=
using wrapper_function = void (* )(eval_level);
struct function_base; // a type derived from |value_base|, defined in \.{axis.w}
using shared_function = std::shared_ptr<const function_base>;
// specialises |shared_value|
template <bool variadic> struct builtin_value;
// derived from |function_base|, defined in \.{axis.w}
using builtin = builtin_value<false>;
using shared_builtin = std::shared_ptr<const builtin>;
using shared_variadic_builtin = std::shared_ptr<const builtin_value<true> >;
struct special_builtin; // derived from |builtin|, defined in \.{axis.w}

@* Implicit conversion of values between types.
%
When interfacing this generic interpreter with a concrete library such as that
of the Atlas of Lie Groups and Representations, a mechanism must be provided to
convert data in the interpreter (represented essentially as nested lists) into
the internal format of the library. For this mechanism to be transparent to the
user, we have chosen to provide the conversion through implicit operations that
are accompanied by type changes; thus when the user enters a list of lists of
integers in a position where an integral matrix is required, the necessary
conversions are automatically inserted during type analysis. In fact we install
a general mechanism of automatic type conversions. This will for instance also
provide the inverse conversions of those just mentioned, and in some instances
merely provide convenience to the user; the latter is the case for instance when
we allow integers in positions where rational numbers are required.

To this end, the function |coerce| requires two fully determined types
|from_type| and |to_type|, and its final argument~|e| is a reference to the
previously converted expression. If a conversion of value of |from_type| to
|to_type| is available, then |coerce| will modify |e| by insertion of a
conversion around it; the return value of |coerce| indicates whether an
applicable conversion was found. The function |conform_types|, available in two
forms, first tries to specialise the type |required| by the context to the one
|found| for the expression itself; if this fails it then tries to
coerce |found| to |required|. If the latter is the case, it wraps the translated
expression |d| in a call of the conversion function found. Should both attempts
fail an, then it trows an error mentioning the expression~|e|.

The function |row_coercion| specialises, if possible, |component_type| in such a
way that the type ``row-of |component_type|'' can be coerced to |final_type|,
and returns a pointer to the |conversion_record| for the coercion in question.
The function |coercion| serves for filling the coercion table.

@< Declarations of exported functions @>=

struct source_location; // defined in \.{parsetree.w}, remains incomplete here

bool coerce(const type_expr& from_type, const type_expr& to_type,
            expression_ptr& e, const source_location& loc);
bool coercible(const type_expr& from_type, const type_expr& to_type);

expression_ptr conform_types
  ( const type_expr& found, type_expr& required
  , expression_ptr&& d, const expr& e);
expression_ptr conform_types
  ( const type& found, type_expr& required, unsigned int fix_count
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

At runtime an implicit conversion is essentially a function call, so we catch an
re-throw errors after adding a line of back-tracing information, similarly to
what we shall do for function calls (as implemented in the \.{axis.w} module).
This only applies to errors that occur during attempted conversion, not those
that might occur during the evaluation of |exp|, so the |try| block below
limited to the conversion call itself.

@< Function def...@>=
void conversion::evaluate(level l) const
{ exp->eval();
  try {@; (*conversion_type.convert)(); }
  @< Catch block for failing implicit conversion @>
  if (l==level::no_value)
    execution_stack.pop_back();
}
@)
void conversion::print(std::ostream& out) const
@+{@; out << conversion_type.name << ':' << *exp; }

@ The errors that might be thrown and temporarily caught here can mostly arise
vector/matrix related conversions (problems while converting big integer or
rational values to bounded size internal representation, or shape problems for
matrices), though the Atlas specific implicit conversion from string to Lie type
can also fail; in all cases however, the error is detected outside the library,
so there is no need to catch |std::exception| (which is what library functions
would throw) here. Contrary to function calls, an implicit conversion stores no
information about the source location for the place where it is invoked, so in
case of an error we can only report the kind of conversion that was attempted.

@< Catch block for failing implicit conversion @>=
catch (error_base& e)
{
  std::ostringstream o;
  o << "In implicit conversion of kind " << conversion_type.name << '.';
  e.trace(o.str());
  throw;
}

@*1 Coercion of types.
%
An important aspect of automatic conversions is that they can be applied in
situations where the result type and only part of the source type is known: for
instance one can be inserted when a list display occurs in a context requiring a
vector, because no row type equals the (primitive) vector type. While this use
requires particular consideration according to the syntactic form of the
expression (here a list display), there is also a simpler form of automatic
conversion that can be applied to a large variety of expressions (identifiers,
function calls, \dots) whenever they are found to have a different type from
what the context requires. The function |coerce| will try to insert an automatic
conversion in such situations, if this can resolve the type conflict. We present
this mechanism first, since the table it employs can then be re-used to handle
the more subtle cases of automatic conversions.

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
type analysis, for instance to allow a conditional expression in a void context
to have branches that evaluate to different types (including the possibility of
absent branches which will be taken to deliver an empty tuple): those branches
which do not already have void type will get their resulting value voided, so
that in the end all branches share the void type of the entire conditional.

At runtime, voiding is mostly taken care of by the |level| argument~|l| passed
around in the evaluation mechanism. Any subexpressions with imposed void type,
like the expression before the semicolon in a sequence expression, will get
their |evaluate| method called with |l==level::no_value|, which suppresses the
production of any value on the |execution_stack|. This setting will be
inherited down to for instance the branches, if the subexpression was a
conditional, so nothing special needs to be done to ensure that branches which
originally had non-void type get voided. However, there are some rare cases
where the void type does not derive from the syntactic nature of the context,
such as in the right hand side of an assignment to a variable that happens to
have |void| type (a quite useless possibility, but valid). In these cases the
type analysis will have to explicitly insert a mechanism to avoid a value to
be produced (and in the example, assigned) where none was intended.

The |voiding| expression type serves that purpose. It distinguishes itself
from instances of |conversion|, in that |voiding::evaluate| calls |evaluate|
for the contained expression with |l==level::no_value|, rather than with
|l==single_value| as |conversion::evaluate| does. In fact, that is about all
that |voiding| is about.

@< Type definitions @>=
class voiding : public expression_base
{ expression_ptr exp;
public:
  voiding(expression_ptr e) : exp(e.release()) @+{}
  virtual void evaluate(level l) const;
  virtual void print(std::ostream& out) const;
};

@ As mentioned, the main point of the |voiding::evaluate| method is that is
calls |void_eval|, which passes |level::no_value| regardless of the |level|
argument it was itself called with. Nonetheless it should not ignore that
argument completely: when called with |l==single_value|, an actual empty tuple
should be produced on the stack, which |wrap_tuple<0>()| does. There is no need
to contribute back-tracing information about the voiding in case of an error at
runtime, since the voiding conversion itself cannot fail.


@< Function definitions @>=
void voiding::evaluate(level l) const
{@; exp->void_eval();
  if (l==level::single_value)
    wrap_tuple<0>();
}
@)
void voiding::print(std::ostream& out) const
@+{@; out << "voided:" << *exp; }


@ The function |coerce| simply traverses the |coerce_table| looking for an
appropriate entry, and wraps |e| into a corresponding |conversion| if it finds
one by calling the |do_conversion| function. Relegating this to a function
defined in the compilation unit~\.{axis.w} (which deals with actual compilation
and execution) is justified by the fact that when applied to a constant argument
value, the value conversion will actually be performed on the value rather than
wrapped in a |conversion| structure (constant folding); therefore this
functionality is not as expression-unaware as the |coerce| function used to be.
We also define a variant function |coercible| that just evaluates the condition,
omitting any action.

When |to_type==void_type|, the conversion always succeeds, as the syntactic
voiding coercion is allowed in all places where |coerce| is called. We now
perform no action in that case, although it would seem safer to, ans we used to
do, wrap the result in a |voiding| in this case to ensure no value will be left
on the runtime stack during evaluation. However as argued above, in most such
cases no runtime action is actually needed, since the context will have already
set |l==level::no_value| for subexpressions with void type. In not providing any
|voiding| here, we have moved responsibility to type analysis, which must now
ensure that any subexpression with void type will have |l==level::no_value| when
evaluated. This means that in some cases a |voiding| will have to be inserted
explicitly. This crops up in many places (for instance a component in a tuple
display just might happen to have void type), but almost all such cases are
far-fetched. So we gain in efficiency of execution by only applying a |voiding|
in such rare cases, but the price to pay is quite a bit of code in our
interpreter to ensure we do in fact find and properly handle all those
exceptional cases.

@h "parse_types.h" // for |source_location|
@h "axis.h" // for |do_conversion|

@< Function definitions @>=
bool coerce(const type_expr& from_type, const type_expr& to_type,
	    expression_ptr& e, const source_location& loc)
{ if (to_type==void_type)
  {@;
     return true;
  } // syntactically voided here, |e| is unchanged
  for (const auto& entry : coerce_table)
    if (from_type==*entry.from and to_type==*entry.to)
    @/{@;
      do_conversion(e,entry,loc);
      return true;
    }
  return false;
}
@)
bool coercible(const type_expr& from_type, const type_expr& to_type)
{ auto match = @[ [&from_type,&to_type]
    @/(const conversion_record& cr)
     {@; return from_type==*cr.from and to_type==*cr.to; } @];
  return to_type==void_type or
    std::any_of(coerce_table.begin(), coerce_table.end(), match);
}

@ Often we first try to specialise a required type to the available type of a
subexpression, or else (if the first fails) coerce the available type to the one
required. The function |conform_types| will facilitate this. The argument |d| is
a possibly already partially converted expression, which should be further
wrapped in a conversion call if appropriate, while |e| is the original
expression that should be mentioned in an error message if both attempts fail. A
call to |conform_types| will be invariably followed (upon success) by returning
the expression that was passed as argument~|d|, possibly modified by
|conform_types|, from that calling function. By returning the modified~|d| from
|conform_types|, we can allow the caller to combine the two by writing |return
conform_types(...)| (knowing that the action of returning might be aborted by an
error thrown from |conform_types|). To indicate that the expression~|d| is
incorporated into to return value, we choose to get |d| passed by rvalue
reference, even though the argument will usually be held in a variable.

If both attempts to conform the types fail we throw a |type_error|; doing so
we must take a copy of |found| (since it a qualified |const|), but we can move
from |required|, whose owner will be destructed before the error is caught.

A second form of |conform_types| is used when a wrapped up, possibly
polymorphic, |found| type needs to be compared against a |required| type. It is
used for applied identifiers, function calls, and on other occasions where
|found| is a type determined independently of the current context and could be
polymorphic. (All forms of assignment expressions use the first form of
|conform_types| for returning there value, because we have forbidden assigning
to identifiers with polymorphic type.) Extra care is needed here
since polymorphic type variables from the stored type |found| might conflict
with those present in |required| as abstract (fixed) type variables currently in
scope. For this reason the current limit |fix_count| is passed as argument, and
our polymorphic type variables will be renumbered to start at that limit. The
actual renumbering takes place in the method |type::unify_specialise| that is
called here.

@< Function def... @>=
expression_ptr conform_types
(const type_expr& found, type_expr& required, expression_ptr&& d, const expr& e)
{ if (not required.specialise(found) and not coerce(found,required,d,e.loc))
    throw type_error(e,found.copy(),required.copy());
  return std::move(d); // invoking |std::move| is necessary here
}
expression_ptr conform_types
  (const type& found, type_expr& required, unsigned int fix_count,
   expression_ptr&& d, const expr& e)
{ if (not found.unify_specialise(required,fix_count) and @|
      not coerce(found.unwrap(),required,d,e.loc))
    throw type_error(e,found.bake(),required.copy());
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
{ assert(component_type==type_expr());
  for (const auto& entry : coerce_table)
    if (final_type==*entry.to and entry.from->raw_kind()==row_type)
      return component_type.specialise(entry.from->component_type())
        ? &entry
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

Our rules for coercions and overloading will be governed by a function
|is_close| computing several relations between pairs of (argument) types; its
resulting value, a small integer, encodes these in separate bits. These tell
whether for a given \foreign{a priori} type another (operand) type provides a
viable candidate, and whether two types can coexist as operand types for a same
overloaded operator or function. (Multiple operands or arguments are considered
as one argument with the tuple type formed from their individual types.)

It may be noted that this function was not changed by the introduction of second
order types into the \axis. language. Making a substitution for a free type
variable, thus producing a less general type, is similar to coercion in that it
changes the type of an expression, and is a conversion that, amongst other uses,
can be applied to the operand/argument expression in a formula or call. It does
not involve any run-time action though, and the reason that this kind of
conversion is not taken into account in |is_close| is that in operator
overloading they are treated as producing exact matches, unlike actual type
coercions. The rules for exact matches are that in concrete instances there must
be one exact match among all overloaded variants, which means that while
introducing those variants we do not have to worry about the \emph{potential}
for such ambiguities, where we do for coercions.

There are only two uses of |is_close| from outside this \.{axis-types.w} module:
in the interpreter it is used to identify overloaded variants that might provide
an inexact match if no exact matches were found, and when adding overloads to
the table it is used to order overloaded definitions so that the more
restrictive ones are tested before ones that would accept a superset of
potential arguments; in the latter use it also helps to detect (and reject)
conflicting overloads that could otherwise allow ambiguous inexact call
expressions. However, |is_close| is also used below in the definition of the
|accepts| relation.

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

We used to worry here about arguments written as empty row display \.{[]}, as it
would be convertible to any row type and therefore make any two concrete row
types incompatible as argument types of same overloaded function, because their
presence allows an ambiguous call. However, now that we use polymorphic types
the empty row display gets the polymorphic type~\.{[A]}, which (like any
polymorphic type) never plays a role in type coercion; it can therefore only
ever give an exact type match, and if in a concrete call it gives more than one
match, that call will be forbidden by the rule that exact matches must be
unique.

In short, there is nothing involving polymorphic types for |is_close| to
worry about. In fact in all cases the arguments types will both be extracted
from |type| values (which means |is_unstable| cannot hold for them: one does
not need to worry about undetermined components either), which types
moreover have been tested to be monomorphic.

@ So here is the (recursive) definition of the relation |is_close|. Equal types
are always close, while undetermined types are not convertible to any other
type. A primitive type is in the relation |is_close| to another type only if it
is identical or if there is a direct conversion between the types, as decided
by~|coerce|. Two row types are in the relation |is_close| if their component
types are, and two tuple types are so if they have the same number of component
types, and if each pair of corresponding component types is (recursively) in the
relation |is_close|. Tabled types are expanded, except if they are recursive in
which case (the case of equality having already been considered) we just say no.

The above applies to the most significant of the three bits used in the result
of the function |is_close|. The two other bits indicate whether, by applying
conversions from~|coerce| to corresponding (nested) component types of row and
tuple types, one of the types can be converted to the other. If the |is_close|
relation is false all bits will be zero, but in the contrary case any of the 4
possible combinations of the remaining bits could be set.

@< Function definitions @>=
unsigned int is_close (const type_expr& x, const type_expr& y)
{ if (x==y)
    return 0x7; // this also makes recursive types equal to themselves
  auto xk=x.raw_kind(), yk=y.raw_kind();
  if (xk==undetermined_type or yk==undetermined_type)
    return 0x0;
      // undetermined types do not specialise (or coerce), and are not close
  if (x==void_type or y==void_type)
    return 0x0; // |void| does not allow coercion for overload, and is not close
  if (xk==tabled)
    return x.is_recursive() ? 0x0 : is_close(x.expanded(),y);
  if (yk==tabled)
    return y.is_recursive() ? 0x0 : is_close(x,y.expanded());
@)
  if (xk==primitive_type or yk==primitive_type)
  { unsigned int flags=0x0;
    if (coercible(x,y)) flags |= 0x1;
    if (coercible(y,x)) flags |= 0x2;
    return flags==0 ? flags : flags|0x4;
  }
  if (xk!=yk)
    return 0x0;
  if (xk==row_type)
    return is_close(x.component_type(),y.component_type());
  if (xk!=tuple_type)
    return 0x0; // non-aggregate types are only close if equal
  unsigned int flags=0x7; // now we have two tuple types; compare components
  wtl_const_iterator it0(x.tuple());
  wtl_const_iterator it1(y.tuple());
  while (not it0.at_end() and not it1.at_end() @|
         and (flags&=is_close(*it0,*it1))!=0)
  @/{@; ++it0; ++it1; }
  return it0.at_end() and it1.at_end() ? flags : 0x0;
}

@ For balancing we need a related but slightly different operation on types. For
|type| values |a,b|, the call |join_to(a,b)| will change |a| if possible to a
type that also accepts expressions of type |b| (possibly no change at all is
necessary), and return whether it could do so. Unlike |is_close|, this operation
does take polymorphic types into account, where |a| accepts a polymorphic type
|b| from which |a| can be obtained by substitution for polymorphic type
variables. But for monomorphic types ``|a| accepts |b|'' holds when there is a
coercion from type~|b| to type~|a|. The separation between polymorphic and
monomorphic behaviour of |join_to| is awkward, but related to a similar
separation in the behaviour of |balance|. The type |type::bottom(c)|, which
unifies to any at all, is a neutral element for |join_to|, and will be used as
initial value in balancing; it is this type that provides the component type of
an empty row display.

We pass |b| by rvalue reference so that we can move this type into |a| in the
case where |b| is the more accepting (monomorphic) type.

@< Declarations of exported functions @>=
bool join_to (type& a, type&& b);

@~The implementation of |join_to| mainly dispatches between unification methods
in the case of polymorphic types, and a local function |accepts| that handles
detection of possible coercions in the monomorphic case.

@< Function definitions @>=
bool join_to (type& a, type&& b)
{
  assert (a.is_clean() and b.is_clean()); // caller should ensure this
  assert(a.floor()==b.floor()); // caller ensures both match the current scope
  if (a.is_polymorphic())
    return a.unify(b);
  if (b.is_polymorphic())
    return b.has_unifier(a.unwrap());
  if (accepts(a.unwrap(),b.unwrap()))
    return true; // nothing to do
  if (accepts(b.unwrap(),a.unwrap()))
    {@; a = std::move(b); return true; }
  return false;
}

@ For locating possible coercions, we define a partial order |accepts| on
monomorphic types (this function is called only from |join_to| above in cases
where polymorphic types have been singled out; there might be fixed type
variables which will be treated like primitive types without associated
coercions). This relation holds whenever there is a coercion from |b| to |a|,
but also when certain expressions of type |b| can be converted into one of type
|a| by inserting coercions inside the expression. For instance
$\\{accepts}(\.{[Split]},\.{[int]})$ holds because a list display of \.{int}
expressions becomes of type \.{[Split]} after inserting coercions of each entry
individually to \.{Split}. If |accepts(a,b)| holds, then an expression that gave
type~|b| can be reconsidered in a context requiring type~|a|, and this second
tentative might succeed (with inserted conversions) or still fail, with an error
being reported at that time.

The implementation of |accepts| is by structural recursion; for this reason it
uses |type_expr| arguments rather than |type|. It is like |is_close|, but some
details are different. We do take |void_type| into consideration here, as type
that will accept any type. Being a partial ordering, we forbid one direction of
all two-way coercions; such cases always involve exactly one primitive type,
and we make that one accept the other type, which does not in return accept the
primitive one. For the rest we just do structural recursion (accepting only if
top-level structure matches and all components accept), but with one twist: for
a pair function types, the parameter types must be equal rather than just
accepting. The reason is that the is no way any expression with \foreign{a
priori} type a function type can be modified by inserting coercions to one with
a different argument type, whereas a change in return type is possible, even
though it is a fairly hypothetical possibility.

@< Local function definitions @>=
bool accepts (const type_expr& a, const type_expr& b)
{
  auto ak=a.raw_kind(), bk=b.raw_kind();
  if (a==void_type or bk==undetermined_type)
    return true; // |void| accepts every type, everything accepts \.*
  if (b==void_type)
    return false; // nothing else accepts |void|
  if (a==b)
    return true;
  if (ak==tabled)
    return not a.is_recursive() and accepts(a.expanded(),b);
  if (bk==tabled)
    return not b.is_recursive() and accepts(a,b.expanded());
  if (ak!=primitive_type and ak!=bk)
    return false;
    // different kinds on non-primitive types are incomparable
@)
  switch(ak)
  {
  case undetermined_type: return false; // should not happen
  case primitive_type:
    return (is_close(a,b)&0x2)!=0; // whether |b| can be converted to |a|
  case variable_type: // assume these fixed; polymorphic types do not get here
    return a.typevar_count()==b.typevar_count();
  case row_type:
    return accepts(a.component_type(),b.component_type());
  case function_type:
    return a.func()->arg_type==b.func()->arg_type and @|
    accepts(a.func()->result_type,b.func()->result_type);
  case tuple_type: case union_type:
  { wtl_const_iterator itb(b.tuple());
    for (wtl_const_iterator ita(a.tuple());
         not ita.at_end(); ++ita,++itb)
      if (itb.at_end() or not accepts(*ita,*itb))
        return false;
    return itb.at_end(); // if list of |a| has ended, that of |b| must as well
  }
  default: assert(false); // should not happen: tabled types were expanded
    return false; // keep compiler happy
  }
}

@* Error values.
%
Before we describe evaluation of expressions we must realise that evaluation
can cause runtime errors.

@< Includes needed in \.{axis-types.h} @>=
#include <stdexcept>
#include <sstream>

@ We use classes derived from |std::exception| and similar standard ones like
|std::runtime_error|, but we define our own local hierarchy, with
|atlas::interpreter::error_base| as base class. The main reason to do this is to
have a centralised error |message| to which exception handlers have write
access, so that it is possible to extend the error message and then re-throw the
same error object. The simplest way to allow this is to give public access to
that string member, so we make this a |struct| rather than a |class|. For errors
happening during evaluation, we want to have the error object thus collect,
during stack unwinding, a trace of the dynamic context (functions being called
and so forth) in which the error occurred, but we want to keep it apart from the
root error |message|, to avoid clobbering the display with possibly massive
output. Therefore a second field |back_trace| with a layered structure is added,
and a method |trace| is provided for extending it. Although this feature is
mostly used with the class |program_error| to be derived below, we provide the
functionality in this base class, as question whether a back trace might be
present in fact depends in principle more on the place where the error is caught
(for instance if type checking produces an error, there can be no back trace)
than on the type of the error itself, even if the two often go together.

We used to provide a method here to extend the message, passing the old message
to first to a temporary |std::ostringstream| object, writing to it, and then
exporting the result as an extended message again. It might seem preferable to
store an |std::ostringstream| object permanently, but that type is really ill
suited to be part of an error value, notably because it is impossible to
implement the |what| method so as to present its contents: the |str| method can
return the contents as a |std::string|, but unless that string is then stored
separately in the error value, its lifetime will be to short to produce a valid
|(char*)| pointer to be returned from |what|. So we finally decided the only
reasonable way to proceed is to not provide a string extension method here, and
instead require the caller to construct (just before throwing) an error string
using a temporary |std::ostringstream| and call its |str| method while throwing
(which string is then copied into the error object); this pattern will occur
repeatedly.

@< Type definitions @>=
struct error_base : public std::exception
{ std::string message;
  simple_list<std::string> back_trace;
  explicit error_base(const std::string& s) : message(s),back_trace() @+{}
  explicit error_base(std::string&& s) : message(std::move(s)),back_trace() @+{}
  error_base (error_base&& other) = default;
  void trace (std::string&& line) @+{@; back_trace.push_front(std::move(line)); }
  const char* what() const throw() @+{@; return message.c_str(); }
};

@ We classify errors into three classes: those due to inconsistency of our
(rather than the user's) program are classified |logic_error|, those arising
during the analysis of the user program are are classified |program_error|,
and those not caught by analysis but during evaluation are classified
|runtime_error|. The first and last are similar to exceptions of the same
name in the |std| namespace, but they are not derived from those exception
classes.

@< Type definitions @>=
struct logic_error : public error_base
{ explicit logic_error(const std::string& s) : error_base(s) @+{}
  explicit logic_error(std::string&& s) : error_base(std::move(s)) @+{}
  logic_error (logic_error&& other) = default;
};
@)
struct program_error : public error_base
{ explicit program_error(const std::string& s) : error_base(s) @+{}
  explicit program_error(std::string&& s) : error_base(std::move(s)) @+{}
  program_error (program_error&& other) = default;
};
@)
struct runtime_error : public error_base
{ explicit runtime_error(const std::string& s) : error_base(s) @+{}
  explicit runtime_error(std::string&& s) : error_base(std::move(s)) @+{}
  runtime_error (runtime_error&& other) = default;
};

@ We derive from |program_error| an exception type |expr_error| that stores in
addition to the error message a reference to an expression to which the
message applies. Placing a reference in an error object may seem hazardous,
because the error might terminate the lifetime of the object referred to, but
in fact in the way we use it, it is safe. Nearly all |expr| objects are
constructed in dynamic memory during parsing, and destructed at the disposal of
the now translated expression at the end of the main interpreter loop; all
throwing of |expr_error| (or derived types) happens after the parser has
finished, and the corresponding |catch| happens in the main loop before disposal
of the expression, so the reference certainly survives the lifetime of the
|expr_error| object. There are a couple of exceptions to this, where a temporary
|expr| value is built during type analysis (the recursive function
|convert_expr| defined in \.{axis.w}), which we arrange to be held in a |static|
variable (with some effort to make this safe in the recursion) at the point
where an error might be thrown, to protect it from being destructed before the
|expr_error| is caught.

The error type is declared a |struct|, so that the |catch| clause may access
the |offender| expression for use in an error message. At the point where this
error value is constructed (and thrown), we just provide a general error
message |s| for storage in the |program_error| base class, and the offending
expression. There is no need to override either of the virtual methods here.

@< Type definitions @>=
struct expr_error : public program_error
{ const expr& offender; // the subexpression causing a problem
@)
  expr_error (const expr& e,const std::string& s) noexcept
    : program_error(s),offender(e) @+{}
  expr_error (const expr& e,std::string&& s) noexcept
    : program_error(std::move(s)),offender(e) @+{}
  expr_error (expr_error&& other) = default;
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
  type_error (type_error&& e) = default;
};

@ For type balancing, we shall use controlled throwing and catching of errors
inside the type checker, which is facilitated by using a dedicated error type.
If balancing ultimately fails, this error will be thrown uncaught by the
balancing code, so |catch| blocks around type checking functions must be
prepared to report the types that are stored in the error value.

When a |balance_error| object is constructed, a descriptive name for the
items being balanced (branches or components of some type of clause) is passed,
which is recorded in the error message. The list of variants is left empty at
construction, but will be filled before actually throwing the |balance_error|.

@< Type definitions @>=
struct balance_error : public expr_error
{ containers::sl_list<type> variants;
  balance_error(const expr& e, const char* items_name)
  : expr_error(e,"No common type found between "), variants()
  @/{@; message+=items_name; }
  balance_error (balance_error&& other) = default;
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
complex_lie_type_type , root_datum_type, Weyl_element_type,
inner_class_type, real_form_type,
Cartan_class_type, KGB_element_type, block_type,
split_integer_type, K_type_type, K_type_pol_type,
module_parameter_type, virtual_module_type, @[@]


@~The following list must match that of the previous module, for proper
functioning of I/O.

@< Other primitive type names @>=
"vec", "mat", "ratvec",
"LieType","RootDatum", "WeylElt", "InnerClass", "RealForm",
"CartanClass", "KGBElt", "Block",
"Split", "KType", "KTypePol", "Param", "ParamPol", @[@]

@* Index.

% Local IspellDict: british
