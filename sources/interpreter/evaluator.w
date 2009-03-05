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
This file describes the unit \.{evaluator}, which implements an evaluator of
expressions, that were produced by the parser with the help of the types and
functions defined in the unit \.{parsetree}.

@h "evaluator.h"
@h <cstdlib>
@c
namespace atlas { namespace interpreter {
namespace {
@< Local type definitions @>@;
@< Declarations of local functions @>@;
}@; // |namespace|
@< Global variable definitions @>@;
namespace {
@< Local variable definitions @>@;
@< Local function definitions @>@;
}@; // |namespace|
@< Function definitions @>@;
}@; }@;

@ As usual the external interface is written to the header file associated to
this file.


@( evaluator.h @>=

#ifndef EVALUATOR_H
#define EVALUATOR_H

@< Includes needed in the header file @>@;
namespace atlas { namespace interpreter {
@< Type definitions @>@;
@< Declarations of global variables @>@;
@< Declarations of exported functions @>@;
@< Template and inline function definitions @>@;
}@; }@;
#endif

@ We start by defining a small template function to help giving sensible error
messages. A return type of |std::string| is useful in throwing for instance
|std::runtime_error|.

@< Includes needed in the header file @>=
#include <string>
#include <sstream>

@~The template function |num| just returns an integer (or any printable value)
represented as a string (in decimal). The reason for using a template is that
without them it is hard to do a decent job for both signed and unsigned types.
Using string streams, the definition of~|num| is trivial.

@< Template and inline function definitions @>=
template <typename T>
  std::string num(T n) @+{@; std::ostringstream s; s<<n; return s.str(); }

@* Types used for evaluation.
%
The parser produces a parse tree, in the form of a value of type |expr|
defined in the unit \.{parsetree}. Our task is to take such an expression,
analyse it and then evaluate it to obtain a value. At various points errors
can occur, which we shall have to handle gracefully: during the analysis
static ``type'' errors may prevent us from undertaking any meaningful action,
and in absence of such errors there still can be dynamic ``runtime'' errors.

We first define some data types used by the evaluator. First of all we shall
consider ``type'' values, represented by the structure |type_declarator|, in
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


@< Include... @>=
#include <memory> // for |std::pair|, |std::auto_ptr|

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
|type_declarator| that will be detailed later. That name indicates that they
expose the structure of the represented types rather than attempting to hide
it in a closed abstraction (as \Cpp~classes do). Note however that classes
defined in the Atlas software library itself will be presented to the user as
primitive types, so we do provide abstraction there. A simple but important
example is a function type that specifies the types of the individual function
arguments and returned values.

Different types will not share any sub-trees among each other, which simplifies
memory management. This choice means some recursive copying of tree structures
is sometimes required, but we avoid doing so more than absolutely necessary.

Usually types are built from the bottom up (from leaves to the root), although
during type checking the reverse also occurs. In bottom-up handling of types,
they should not reside in local |type_declarator| variables, since these would
require (recursive) copying in order to become a descendant type of another
type. Instead they will be referred to by local pointer values, and in order
to ensure exception-safety, these should be auto-pointers; for this reason we
define |type_ptr| to be an auto-pointer type. For the links between a node and
its descendant types we use ordinary types however. We also define a
|type_list| type that will be used inside tuple types, and a corresponding
auto-pointer type |type_list_ptr|.

@< Type definitions @>=
struct type_declarator;
typedef std::auto_ptr<type_declarator> type_ptr;
typedef struct type_node* type_list;
typedef std::auto_ptr<type_node> type_list_ptr;

@ Since types and type lists own their trees, their copy constructors must
make a deep copy. Most invocations of the copy constructor for types, from
outside the implementation of the |type_declarator| class itself, will take
place by an explicit call of a function |copy| for improved visibility. This
function converts a reference to a |type_declarator| into an auto-pointer to a
copied instance of that type declarator.

@< Declarations of exported functions @>=
type_ptr copy(const type_declarator& t);

@~The function |copy| simply creates an auto-pointer from the call of |new|
for |type_declarator(t)|, where the latter invokes the copy constructor of
|type_declarator|. The latter, has a recursive definition (given in section
@#type declarator copy@>) that performs a deep copy, but the details don't
concern us here.

@< Function definitions @>=
type_ptr copy(const type_declarator& t) @+
{@; return type_ptr(new type_declarator(t)); }

@*2 Type lists.
%
The auxiliary type |type_list| gives a preview of matters related to
ownership, that will also apply, in a more complex setting, to
|type_declarator|. It is a simply linked list, whose nodes contain a
|type_declarator| as a sub-object. We might instead have used a pointer there,
but that would be less compact and require more memory (de)allocations. On the
other hand it forced us to introduce the |type_declarator::set_from| method
(see below) to make a shallow copy, avoiding the deep copy that would result
from initialising |t(*head)|. Note how the pointer on the new node is passed
as an auto-pointer that is immediately released when inserted into the node;
this is the unique way to handle this that assures there is exactly one owner
at all times (if an ordinary pointer had been passed, the pointed-to structure
might fail to be deleted if the allocation of this node should throw an
exception).

Because a |type_node| contains a
|type_declarator| we must make sure its definition given here \emph{follows}
that of |type_declarator|, at least in the tangled code.

@< Definition of |struct type_node@;| @>=
struct type_node
{ type_declarator t; @+ type_list next;
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

We pass the auto-pointers by value, rather than by (necessarily non-constant)
reference. This implies a small overhead due to multiple ownership transfer,
but allows (due to a clever definition of auto-pointers) passing anonymous
values obtained from a constructor or construction function as argument, which
would be forbidden if using non-constant reference parameters.


@< Local function def... @>=

type_list_ptr make_type_singleton(type_ptr t)
{@; return type_list_ptr(new type_node(t,type_list_ptr(NULL))); }
@)
type_list_ptr make_type_list(type_ptr t,type_list_ptr l)
{@; return type_list_ptr(new type_node(t,l)); }

@ Here is one more function that is convenient to have around.

@< Local function definitions @>=
size_t length(type_list l)
@+{@; size_t len=0;
  for (; l!=NULL; l=l->next) ++len;
  return len;
}

@*2 Primitive types.
%
Types can be formed in a variety of ways (primitive types, tuple types,
function types, etc.), and correspondingly will involve a |union| of various
alternatives. So before we can define |type_declarator|, we need to enumerate
its variants and the possibilities for primitive types. Some of them will be
introduced later.

@< Type definitions @>=
enum type_tag
{ undetermined_type, primitive_type, row_type, tuple_type, function_type };

enum primitive_tag
{ integral_type, string_type, boolean_type
  , @< Other primitive type tags @>@;@;
  @+nr_of_primitive_types };

@ For printing types (and later for parsing them) we shall need names for the
primitive ones.

@< Declarations of global variables @>=

extern const char* prim_names[];

@~Here are the names of the primitive types already enumerated.

@< Global variable definitions @>=
const char* prim_names[]=
{"int","string","bool",@< Other primitive type names@>@;@; };

@*2 Type declarators.
%
We have a simple but flexible type model. There is a finite number of
``primitive'' types, many of which are abstractions of complicated classes
defined in the Atlas library, such as root data or reductive groups.
Furthermore one has types that are ``row of'' some other type, tuples
(Cartesian products) of some given sequence of types, and function types. We
also allow for an undetermined type, which can serve as a wild-card to specify
type patterns. Here are enumerations of tags for the basic kinds of types, and
for the sub-cases of |primitive_type|.

Type declarators are defined by a tagged union. Code accessing the union will
in general test the tag and access the corresponding variant directly, so we
give public access to the data fields (it is not obvious how coherent variant
selection could be enforced by private data and public methods without greatly
complicating usage; the \Cpp~model of abstraction does not seem particularly
suited for tagged unions). By using an anonymous union, the field selectors
like |prim| of the variants in the union can be used directly on the level of
the |type_declarator|, thus avoiding an additional level of selection to
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

@< Type definitions @>=

struct func_type; // must be predeclared for |type_declarator|
@)
struct type_declarator
{ type_tag kind;
  union
  { primitive_tag prim; // when |kind==primitive|
    type_declarator* component_type; // when |kind==row_type|
    type_list tuple; // when |kind==tuple_type|
    func_type* func; // when |kind==function_type|
  };
@)
  type_declarator() : kind(undetermined_type) @+{}
  explicit type_declarator(primitive_tag p)
    : kind(primitive_type) @+{@; prim=p; }
  explicit type_declarator(type_ptr c)
    : kind(row_type) @+{@; component_type=c.release(); }
  explicit type_declarator(type_list_ptr l) : kind(tuple_type)
  @+{@; tuple=l.release(); }
  inline type_declarator(type_ptr arg, type_ptr result); // for function types
@)
  type_declarator(const type_declarator& t); // copy constructor
  ~type_declarator();
@)
  void set_from(type_ptr p);
  bool specialise(const type_declarator& pattern);
private:
  type_declarator& operator=(const type_declarator& t);
   // assignment is forbidden
};

@< Definition of |struct type_node@;| @>@; // this must \emph{follow}

@ The constructor for function types cannot be defined inside the structure
definition, since |func_type| is not (and cannot be) a complete type there.
The constructor for |func_type| used will be defined below.

@< Function definitions @>=
type_declarator::type_declarator(type_ptr arg, type_ptr result)
   : kind(function_type) {@; func=new func_type(arg,result); }

@ The variant for function types needs both an argument type and a result
type. We cannot define the appropriate structure directly within the |union|
of |type_declarator| where it is needed, since the scope of that definition
would then be too limited to perform the appropriate |new| in the constructor.
Therefore we must declare and name the structure type before using it in
|type_declarator|. Afterwards we define the structure, and while we are doing
that, we also define a constructor to fill the structure (it was used above)
and a copy constructor.

The |func_type| structure has two sub-objects of type |type_declarator|,
rather than pointers to them. This is more compact and practical for most
purposes, but we should take care that the basic constructor only makes a
shallow copy of the topmost nodes pointed to. This is achieved by using the
|set_from| method to fill the initially undetermined slots with a shallow copy
of the types, in the same way as was done in the constructor for~|type_node|.
The copy constructor needs no special care, as it makes a deep copy.

@< Type definitions @>=
struct func_type
{ type_declarator arg_type, result_type;
@)
  func_type(type_ptr a, type_ptr r)
   : arg_type(), result_type()
  {@; arg_type.set_from(a); result_type.set_from(r); }
  func_type(const func_type& f)
   : arg_type(f.arg_type),result_type(f.result_type) @+{}
};


@ As we remarked above, the copy constructor for a |type_declarator|
recursively copies the descendant types; this is necessary because these
descendants are owned by the object. As usual the recursion is implicit in the
use of a copy constructor within |new|. In all cases there is only one object
directly pointed to, so there is no need for intermediate auto-pointers: no
exception can be thrown between the moment that the call to |new| returns (if
it does), and the completion of our constructor, which incorporates the
returned pointer.
@:type declarator copy@>

@< Function definitions @>=
type_declarator::type_declarator(const type_declarator& t) : kind(t.kind)
{ switch (kind)
  { case undetermined_type: break;
    case primitive_type: prim=t.prim; break;
    case row_type:
      component_type=new type_declarator(*t.component_type);
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
type_declarator::~type_declarator()
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
copied to the current |type_declarator|, but without invoking the copy
constructor; afterwards the original node will be destroyed after detaching
any possible descendants from it. This operation is only safe if the
|type_declarator| previously had no descendants, and in fact we insist that it
had |kind==undetermined_type|; if this condition fails we signal a
|std::logic_error|. In a sense this is like a |swap| method, but only defined
if the type was undetermined to begin with; in modern parlance, it implements
move semantics.

@h <stdexcept>
@< Function definitions @>=
void type_declarator::set_from(type_ptr p)
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
pattern.

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

@< Function definitions @>=
bool type_declarator::specialise(const type_declarator& pattern)
{ if (pattern.kind==undetermined_type) return true; // no need to refine
  if (kind==undetermined_type)
    {@;new(this) type_declarator(pattern); return true; }
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

@ For printing types, we shall pass |type_declarator| values to the
operator~`|<<|' by constant reference, which seems more decent than doing
so by pointer (which would override the definition that simply prints the
hexadecimal address of a pointer); we shall not define instances of~`|<<|' for
other pointer types either. Since we often hold types in |type_ptr| values,
this does mean the we must dereference explicitly in printing.

@< Declarations of exported functions @>=
std::ostream& operator<<(std::ostream& out, const type_declarator& t);

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
std::ostream& operator<<(std::ostream& out, const type_declarator& t)
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

@<  Declarations of exported functions @>=
bool operator== (const type_declarator& x,const type_declarator& y);
inline bool operator!= (const type_declarator& x,const type_declarator& y)
{@; return !(x==y); }

@~This code is quite similar to the |specialise| method; in fact one could
often use that method instead of the equality operator, but here we want both
operands to be |const|.

@< Function definitions @>=
bool operator== (const type_declarator& x,const type_declarator& y)
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
which are auto-pointers. The reasons for, and the way of, using them are the
same as explained above for |make_type_list|.

@< Local function definitions @>=
type_ptr make_undetermined_type()
@+{@; return type_ptr(new type_declarator); }
@)
type_ptr make_prim_type(primitive_tag p)
@+{@; return type_ptr(new type_declarator(p)); }
@)
type_ptr make_row_type(type_ptr c)
{@; return type_ptr(new type_declarator(c)); }
@)
type_ptr make_tuple_type (type_list_ptr l)
{@; return type_ptr(new type_declarator(l)); }
@)
type_ptr make_function_type (type_ptr a, type_ptr r)
{@; return type_ptr(new type_declarator(a,r));
}


@*2 Specifying types by strings.
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

@ For primitive types we use the same strings as for printing them. We use the
fact that no name is a prefix of another one, so the first match is decisive.

@< Scan and |return| a primitive type, or |throw| a |logic_error| @>=
{ for (size_t i=0; i<nr_of_primitive_types; ++i)
  { std::string name=prim_names[i];
    if (name.compare(0,name.length(),s,name.length())==0)
    @/{@; s+=name.length();
      return make_prim_type(static_cast<primitive_tag>(i));
    }
  }
  throw std::logic_error("Type unrecognised");
}

@*2 Predefined type declarators.
%
We shall often need to refer to certain types for comparison. Instead of
generating them on the fly each time using |make_type|, we define constant
values that can be used everywhere. We include a few types that will be
introduced later.

@< Declarations of global variables @>=
extern const type_declarator unknown_type; // \.{*}
extern const type_declarator int_type; // \.{int}
extern const type_declarator str_type; // \.{string}
extern const type_declarator bool_type; // \.{bool}
extern const type_declarator vec_type; // \.{vec}
extern const type_declarator mat_type; // \.{mat}
extern const type_declarator row_of_type; // \.{[*]}
extern const type_declarator row_of_int_type; // \.{[int]}
extern const type_declarator row_of_vec_type; // \.{[vec]}
extern const type_declarator row_row_of_int_type; // \.{[[int]]}
extern const type_declarator int_int_type; // \.{(int,int)}

@ The definition of the variables uses the constructors we have seen above,
rather than functions like |make_primitive_type| and |make_row_type|, so that
no dynamic allocation is required for the top level structure. For ``row of''
types we construct the |type_declarator| from an auto-pointer to another one
produce by |copy| applied to a previous |type_declarator|.

@< Global variable definitions @>=

const type_declarator unknown_type; // uses default constructor
const type_declarator int_type(integral_type);
const type_declarator str_type(string_type);
const type_declarator bool_type(boolean_type);
const type_declarator vec_type(vector_type);
const type_declarator mat_type(matrix_type);
const type_declarator row_of_type(copy(unknown_type));
const type_declarator row_of_int_type(copy(int_type));
const type_declarator row_of_vec_type(copy(vec_type));
const type_declarator row_row_of_int_type(copy(row_of_int_type));
const type_declarator int_int_type(*make_type("(int,int)"));

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

@< Includes needed in the header file @>=
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

@~The operator~`|<<|' calls the (virtual) |print| method of the object pointed
to, and returns the reference to the output stream object.

@< Function definitions @>=
std::ostream& operator<< (std::ostream& out, const value_base& v)
{@; v.print(out); return out; }

@ Now we derive the first ``primitive'' value types. For each type we define
a corresponding auto-pointer type (whose name ends with \&{\_ptr}), since we
shall often need to hold such values by pointers, and the risk of exceptions
is ever present.

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

struct string_value : public value_base
{ std::string val;
@)
  explicit string_value(const char* t) : val(t) @+ {}
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

@ Here we define one more type derived from |value_base|, namely the type for
``row of'' types.

@< Includes needed in the header file @>=
#include <vector>

@~Row values are implemented using vectors from the standard template library.
Since the actual values accessed will be of types derived from |value_base|,
we must pass through a level of indirection, so we have a vector of pointers.
We define these pointers to be |shared_value| pointers, so that the row takes
(shared) ownership of its components without needing a explicit destructor.
This as the additional advantage over explicit ownership management that the
copy constructor, needed for the |clone| method, can safely just
copy-construct the vector of pointers: a possible exception thrown during the
copy is guaranteed to clean up any pointers present in the vector.

Of course ownership of pointers to |row_value| objects also needs to be
managed, which could be either by an auto-pointer |row_ptr| (if the pointer is
known to be unshared) or by a shared pointer |shared_row|.

@< Type definitions @>=
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

@< Template and inline function definitions @>=
template <typename D> // |D| is a type derived from |value_base|
 D* force(value v) throw(std::logic_error)
{ D* p=dynamic_cast<D*>(v);
  if (p==NULL) throw
    std::logic_error(std::string("forced value is no ")+D::name());
  return p;
}

@*2 Representation of an evaluation context.
%
When we evaluate a user program, values will be given to local identifiers
such as arguments of functions being called. The identification of such
identifiers is static, and determined during type analysis; at execution time
we only need to record the values bound to those identifiers (although having
the identifier names associated to them may be helpful in reporting runtime
errors). The dynamic evaluation context is formed by a stack of blocks, each
holding a number of values.

@< Type definitions @>=
class context;
typedef std::tr1::shared_ptr<context> context_ptr;
class context
{ const context_ptr next;
  const std::vector<shared_value> frame;
  context(const context&); // contexts should not be copied, just shared
public:
  context(const context_ptr& n,
          const std::vector<shared_value>& f) : next(n), frame(f) @+{}
  shared_value elem(size_t i,size_t j) const;
};

@ The method |context::elem| descends the stack and then selects a value from
the proper frame.

@< Function def... @>=
shared_value context::elem(size_t i, size_t j) const
{
  const context* p=this;
  while (i-->0)
  @/{@; assert(p->next.use_count()>0 and p->next.get()!=NULL);
    p=p->next.get();
  }
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
in most cases when calling functions.

@< Type definitions @>=
struct expression_base
{ expression_base() @+ {}
  virtual ~expression_base() @+ {}
  virtual void evaluate() const =0;
  virtual void print(std::ostream& out) const =0;
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
converted into the latter, in which case the |use_count| will becme~$1$. For
convenience we make these template functions that accept a smart pointer to
any type derived from |value_base| (since a conversion of such pointers from
derived to base is not possible without a cast in a function argument
position). For even more convenience we also provide a variant taking an
ordinary pointer, so that expressions using |new| can be written without cast
in the argument of |push_value|. Since |push_value| has only one argument,
such use of does not compromise exception safety: nothing can throw between
the return of |new| and the conversions of its result into an auto-pointer.

@< Template and inline function definitions @>=
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
shared pointer (because the value on the stack might be shared). Finally we
provide a variant template function |get_own| that returns a privately owned
copy of the value from the stack, to which modifications can be made without
danger of altering shared instances. Since there is no way to persuade a
|shared_ptr| to release its ownership even if ith happens to be the unique
owner, we must unconditionally call |clone| to achieve this.

@< Template and inline function definitions @>=

inline shared_value pop_value()
{@; shared_value arg=execution_stack.back();
  execution_stack.pop_back();
  return arg;
}
@)
template <typename D> // |D| is a type derived from |value_base|
 std::tr1::shared_ptr<D> get() throw(std::logic_error)
{ std::tr1::shared_ptr<D> p=std::tr1::dynamic_pointer_cast<D>(pop_value());
  if (p.get()==NULL)
    throw std::logic_error(std::string("Argument is no ")+D::name());
  return p;
}
@.Argument is no ...@>
@)
template <typename D> // |D| is a type derived from |value_base|
 std::auto_ptr<D> get_own() throw(std::logic_error)
{ D* p=dynamic_cast<D*>(execution_stack.back().get());
  if (p==NULL)
    throw std::logic_error(std::string("Argument is no ")+D::name());
  std::auto_ptr<D> result(p->clone());
  execution_stack.pop_back(); // only now can we release the shared pointer
  return result;
}
@.Argument is no ...@>

@ Now let us define classes derived from |expression_base|. The simplest is
|denotation|, which simply stores a |value|, which it returns upon evaluation.
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
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const
  @+{@; denoted_value->print(out); }
};

@ The following function is not defined inside the class definition, because
that definition precedes the one of the |inline| function |push_value| in the
header file. It is not even made |inline| itself, since there is little point
in doing so for virtual methods (calls via the vtable cannot be inlined).

@< Function def... @>=
void denotation::evaluate() const
@+{@; push_value(denoted_value); }


@ Let us finish with some functions related to the execution stack. The stack
owns the values it contains, but there was no reason to wrap it into a class
with a destructor, since we never intend to destroy the stack: if our program
exits either peacefully or by an uncaught exception we don't care about some
values that are not destroyed. We must remember however to |delete| the values
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


@*1 Error values.
Before we describe evaluation of expressions we must realise that evaluation
can cause runtime errors. The evaluator may throw exceptions due to
inconsistency of our (rather than the user's) program, which are classified as
|std::logic_error|. It may also throw exceptions due to errors not caught by
the type checker in the user input, such as size mismatch in matrix
operations; in such cases it will throw |std::runtime_error|.

@< Include... @>=
#include "parsetree.h"
#include <stdexcept>

@ We first define a general exception class |program_error| derived from
|exception|, which represents any kind of error of the user input determined
by static analysis (for instance use of undefined variables).

Although it does nothing explicitly (the string will be destructed anyway), we
must explicitly define a destructor for |program_error|: the automatically
generated one would lack the |throw()| specifier.

@< Type definitions @>=
class program_error : public std::exception
{ std::string message;
public:
  explicit program_error(std::string s) : message(s) @+{}
  ~program_error() throw() @+{@;} // obligatory definition
  const char* what() const throw() @+{@; return message.c_str(); }
};

@ We derive from |program_error| a more specific exception type |type_error|
that we shall throw when an expression fails to have the right type. It stores
the offending subexpression and pointers to the types that failed to match.
For the latter, the error value owns the types pointed to, so the caller
should relinquish ownership of, or copy (in case it did not own), the types
passed when throwing the exception.

@< Type definitions @>=
struct type_error : public program_error
{ expr offender; // the subexpression with a wrong type
  const type_declarator* actual; @+
  const type_declarator* required; // the types that conflicted
@)
  type_error (const expr& e, type_ptr a, type_ptr r) throw() @/
    : program_error("Type error"),offender(e) @|
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
never be used. This is not easy either, since two copies of type declarators
must be made, which can throw exceptions (let us hope the hypothetical call to
the copy constructor is made \emph{before} the |throw| of |type_error| is
executed, or else we will be dead due to interwoven exceptions\dots), and
during the second copy the pointer to the first copy must be held in an
auto-pointer to avoid a memory leak. This excludes any attempt to initialise
|actual| and |required| directly to their proper values, sigh! In the code
below, we could have said instead |actual=new type_declarator(*e.actual)|,
thus avoiding the use of the anonymous auto-pointer produced by the second
call of |copy|, at the price of some loss of symmetry and readability.

@< Function definitions @>=
type_error::type_error(const type_error& e)
 : program_error(e), offender(e.offender)
{@; type_ptr p=copy(*e.required);
  actual=copy(*e.actual).release(); required=p.release();
}


@* The evaluator.
%
The expression returned by the parser is type-checked before execution starts.
During type checking, there are subsecpressions for which a definite type is
required, and others where none is (like the complete expression). The
difference is important, as in the former case conversions can be inserted to
make types match, for instance between a list of integers and a vector value;
this is in fact the only way the user can produce the latter kind of values.
However, both cases are handled by a single function |convert_expr|, which in
addition builds (upon success) an |expression| value, i.e., one in terms of
which evaluation is defined. Apart from the |expr e@;| produced by the parser,
|convert_expr| takes a type as argument, in the form of a non-constant
reference |type_declarator& type@;|. If |type| is undefined initially, then it
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
expression convert_expr(const expr& e, type_declarator& type)
  throw(std::bad_alloc,program_error);

@ In the function |convert_expr| we shall need a type for storing bindings
between identifiers and types. It will be a vector, and since vectors cannot
hold auto-pointers, we need to define a descructor to clean up the types.
Although we shall in fact need to stack different such binding vectors for
nested scopes, we prefer not to use any STL container to implement the stack,
since that would necessitate defining a copy constructor, which would a painful
operation given that (1)~it must call |copy| on the types held in the
bindings, in order to avoid double destruction, (2)~these calls to|copy|
could throw an exception, and (3)~the constructor won't be complete, and the
destructor therefore not activated, until the final entry is copied, so we
would need a |try|\dots|catch| in the constructor to avoid a memory leak.

So instead we chain the different bindings into a linked list, and forbid any
copy or assignement. In fact the nesting of the various bindings will embed in
the nested calls of |convert_expr|, so we can afford to allocate each
|bindings| as local variable to (some instace of) |convert_expr|, and there is
no ownership of them down the linked list. We do provide methods |push| and
|pop| to prepend and detach a |bindings| from a list, arrive at the following
type definition.

@< Type def... @>=
struct bindings
: public std::vector<std::pair<Hash_table::id_type,type_declarator*> >
{ bindings* next; // non-owned pointer

  typedef std::vector<std::pair<Hash_table::id_type,type_declarator*> >
    base;
  bindings(size_t n=0) : @[ base() @], next(NULL) @+{@; base::reserve(n); }
  ~bindings () @+
  {@; for (size_t i=0; i<base::size(); ++i)
      delete operator[](i).second;
  }
  void add(Hash_table::id_type id,type_ptr t)
  {@; push_back(std::make_pair(id,(type_declarator*)NULL));
      back().second=t.release();
  }
  void push (bindings*& sp) @+{@; next=sp; sp=this; }
  void pop (bindings*& sp) const @+{@; sp=next; }
private:
  bindings(const bindings&); // copy constructor
  void operator= (const bindings&); // assignement operator
};

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
we don't return an |expression_ptr| is entirely pragmatic.


@< Function definitions @>=
expression convert_expr(const expr& e, type_declarator& type)
  throw(std::bad_alloc,program_error)
{

  switch(e.kind)
  { case integer_denotation:
      if(type.specialise(int_type)) // change undefined type to integral
        return new denotation
          (shared_value(new int_value(e.e.int_denotation_variant)));
      else throw type_error(e,copy(int_type),copy(type));
   case string_denotation:
      if (type.specialise(str_type)) // change undefined type to string
        return new denotation
          (shared_value(new string_value(e.e.str_denotation_variant)));
      else throw type_error(e,copy(str_type),copy(type));
   case boolean_denotation:
      if (type.specialise(bool_type)) // change undefined type to boolean
        return new denotation
          (shared_value(new bool_value(e.e.int_denotation_variant)));
      else throw type_error(e,copy(bool_type),copy(type));
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
  virtual void evaluate() const;
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

@ The evaluation of a |list_expression| evaluates the components in a simple
loop and wraps the result in a |row_value|. Since evaluation pushes the
resulting value onto the execution stack, we pop it off into the result
afterwards. We take care to hold the partial result via an auto-pointer
|result|, so that in case of a runtime error during the evaluation of one of
the component expressions the values already computed are cleaned up.

@< Function def... @>=
void list_expression::evaluate() const
{ row_ptr result(new row_value(component.size()));
  for (size_t i=0; i<component.size(); ++i)
    component[i]->evaluate(),result->val[i]=pop_value();
  push_value(result); // result will be shared from here on
}

@*1 Conversion to rigid vectors and matrices.
%
The following sections deal particularly with interfacing to the Atlas
library, but address a problem that would be present for other libraries as
well: how can a user of the interpreter construct data in the internal format
of the library, rather than that (essentially nested lists) of the
interpreter? The answer is via conversion routines, calls of which are
implicitly inserted when the type analysis detects their necessity. Currently
the conversions from the type ``row of integer'' to and from the type
``vector', as well as those from the types ``row of row of integer'' and ``row
of vector'' to and from the type ``matrix'' are implemented. Here ``vector''
really means |latticetypes::Weight|, which stands for |matrix::Vector<int>|,
while ``matrix'' means |latticetypes::LatticeMatrix|, which stands for
|matrix::Matrix<int>|. These declarations are given
in~\.{structure/latticetypes\_fwd.h}, while the template class
|matrix::Matrix| is defined in~\.{utilities/matrix.h}.


@< Includes... @>=
#include "latticetypes_fwd.h"
#include "matrix.h" // this makes |latticetypes::LatticeMatrix| a complete type

@*2 Inserting conversion functions.
%
The conversion to vector or matrix will be represented by a special kind of
expression, which behaves like a function call of a fixed function. For each
function a different structure derived from |expression_base| is used; they
only differ among each other by their virtual methods |evaluate| and |print|.

@< Type definitions @>=
struct vector_conversion : public expression_base
{ expression exp;
@)
  explicit vector_conversion(expression_ptr e) : exp(e.release()) @+{}
  virtual ~vector_conversion()@;{@; delete exp; }
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};
@)
struct matrix_conversion : public vector_conversion
{ explicit matrix_conversion(expression_ptr e) : vector_conversion(e) @+{}
@)
  virtual ~matrix_conversion() @+{}
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};

@ When printed, these expressions provide a brief indication of the
conversion.

@< Function definitions @>=
void vector_conversion::print(std::ostream& out) const
@+{@; out << "V:" << *exp; }
void matrix_conversion::print(std::ostream& out) const
@+{@; out << "M:" << *exp; }

@ When in |convert_expr| we encounter a list display when a non-row is
expected, we single out the cases that the required type is vector or matrix,
in which case we continue to convert the component expressions with expected
type integer respectively vector. Here the insertion of the conversion
function just means wrapping the list display into a call to |new
vector_conversion| of |new matrix_conversion|.

@< If |type| can be converted from some row-of type... @>=
  { bool is_vector; type_declarator component_type;
    if (type==vec_type)
    {@; component_type.specialise(int_type); is_vector=true; }
    else if (type==mat_type)
    {@; component_type.specialise(vec_type); is_vector=false; }
    else throw type_error(e,copy(row_of_type),copy(type));
    std::auto_ptr<list_expression> display@|
      (new list_expression@|(length(e.e.sublist)));
  @/size_t i=0;
    for (expr_list l=e.e.sublist; l!=NULL; l=l->next)
      display->component[i++]=convert_expr(l->e,component_type);
    return
      is_vector ? new vector_conversion(expression_ptr(display))
		: new matrix_conversion(expression_ptr(display));
  }

@*2 Coercion of types.
%
The most important instances where conversions are inserted are the cases
above, where a list display occurs in a context requiring a vector or a
matrix, since they provide a way for the user to construct vectors and
matrices. However, it seems logical to provide the same conversions also when
the expression is not a list display, for instance an applied identifier or a
function call. In those cases the need for a conversion will usually arise
from a difference between the required type and an already fully determined
type of the subexpression. Since multiple syntactic classes are involved, it
is not attractive to handle such conversions separately for the various
branches of the big switch on the expression~|kind| in the type analysis, so
we put in place a general mechanism that will insert a conversion depending on
the pair of the expression type and required type, if it exists.

The mechanism is realised by a function |coerce|. Its final argument~|e| is a
reference to the previously converted expression having |from_type|, and if
values of this type can be converted to |to_type| then the expression will be
modified by insertion of a conversion expression; the return value indicates
whether a conversion was found.

@< Declarations of exported functions @>=

bool coerce( const type_declarator& from_type, const type_declarator& to_type
	   , expression_ptr& e) throw(std::bad_alloc);

@ The implementation of |coerce| will be determined by a simple table lookup.
Our convention will be that conversion functions take an auto-pointer as
argument, release it upon insertion of the new node, and return an ordinary
|expression| pointer, which the caller should then assume ownership of (the
latter is just for convenience, to avoid the necessity to convert to an
auto-pointer in each separate conversion function).

@< Local type definitions @>=
struct conversion_info
{ typedef expression (*converter)(expression_ptr);
  type_declarator from,to; converter conv;
    // function transforming the |expression| by inserting the coercion
  conversion_info (type_ptr from_type, type_ptr to_type, converter c)
   : from(),to(), conv(c)
  @/{@; from.set_from(from_type); to.set_from(to_type); }
};
@)
typedef std::vector<conversion_info> coerce_table_type;

@ Here is the lookup table.

@< Global variable definitions @>=

conversion_info coerce_table[]= { @<Initialiser for |coerce_table| @>@;@;};

const size_t nr_coercions=sizeof(coerce_table)/sizeof(conversion_info);

@ We need |converter| functions to instantiate |conversion_info::conv| from.
We need to declare these functions before we come to the global variable
definitions.

@< Declarations of local functions @>=
expression vector_converter(expression_ptr e);
expression matrix_converter(expression_ptr e);

@ The conversion-inserting functions are simple packaging routines.

@< Local function definitions @>=
expression vector_converter(expression_ptr e) @+
{@; return new vector_conversion(e); }
expression matrix_converter(expression_ptr e) @+
{@; return new matrix_conversion(e); }


@ All that remains is to initialise the |coercion table|.
@<Initialiser for |coerce_table| @>=
conversion_info(copy(row_of_int_type), copy(vec_type), vector_converter), @/
conversion_info(copy(row_of_vec_type), copy(mat_type), matrix_converter), @/
@[@]

@ The function |coerce| simply traverses the |coerce_table| looking for an
appropriate entry, and applies the |converter| if it finds one. As indicated
above, it is the task of |coerce| to keep ownership while calling the
conversion function, and switch to owning the full expression afterwards.

@< Function definitions @>=
bool coerce( const type_declarator& from_type, const type_declarator& to_type
	   , expression_ptr& e) throw(std::bad_alloc)
{ for (conversion_info*
       it=&coerce_table[0]; it<&coerce_table[nr_coercions]; ++it)
    if (from_type==it->from and to_type==it->to)
    @/{@; e.reset(it->conv(e));
      return true;
    }
  return false;
}

@ We extend our range of conversions with other possibilities. For each we
need to define a type derived from |expression_base|; in fact we can derive
from |vector_conversion|, like for |matrix_conversion|.

@<Type definitions@>=

struct matrix2_conversion : public vector_conversion
  // for \.{[[int]]}$\to$\.{mat}
{ explicit matrix2_conversion(expression_ptr e) : vector_conversion(e) @+{}
@)
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};

struct int_list_conversion : public vector_conversion
  // for \.{vec}$\to$\.{[int]}
{ explicit int_list_conversion(expression_ptr e) : vector_conversion(e) @+{}
@)
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};
@)
struct vec_list_conversion : public vector_conversion
  // for \.{mat}$\to$\.{[vec]}
{ explicit vec_list_conversion(expression_ptr e) : vector_conversion(e) @+{}
@)
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};
@)
struct int_list_list_conversion : public vector_conversion
  // for \.{mat}$\to$\.{[[int]]}
{ explicit int_list_list_conversion(expression_ptr e)
   : vector_conversion(e) @+{}
@)
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};

@ The print functions are very simple. As before, these definitions were
nevertheless not in-lined in the structure definition, because in the header
file type definitions come before the point where operators are declared,
notably the instance of the operator~`|<<|' used to print the
|expression_base| reference `|*exp|'.

@< Function definitions @>=
void matrix2_conversion::print(std::ostream& out) const
@+{@; out << "MV:" << *exp; }
void int_list_conversion::print(std::ostream& out) const
@+{@; out << "[I]:" << *exp; }
void vec_list_conversion::print(std::ostream& out) const
@+{@; out << "[V]:" << *exp; }
void int_list_list_conversion::print(std::ostream& out) const
@+{@; out << "[[I]]:" << *exp; }

@ Here are the declarations of their conversion-inserting functions.

@< Declarations of local functions @>=
expression matrix2_converter(expression_ptr e);
expression int_list_converter(expression_ptr e);
expression vec_list_converter(expression_ptr e);
expression int_list_list_converter(expression_ptr e);

@~Their definitions are hardly surprising.

@< Local function definitions @>=
expression matrix2_converter(expression_ptr e)
@+{@; return new matrix2_conversion(e); }
expression int_list_converter(expression_ptr e)
@+{@; return new int_list_conversion(e); }
expression vec_list_converter(expression_ptr e)
@+{@; return new vec_list_conversion(e); }
expression int_list_list_converter(expression_ptr e)
@+{@; return new int_list_list_conversion(e); }


@ All that remains is to initialise the |coercion table|.
@<Initialiser for |coerce_table| @>=
conversion_info(copy(row_row_of_int_type), copy(mat_type)
  , matrix2_converter), @/
conversion_info(copy(vec_type), copy(row_of_int_type)
  , int_list_converter), @/
conversion_info(copy(mat_type), copy(row_of_vec_type)
  , vec_list_converter), @/
conversion_info(copy(mat_type), copy(row_row_of_int_type)
  , int_list_list_converter), @/
@[@]


@*2 Primitive types for vectors and matrices.
%
When packed into a a type derived from |value_base| (defined later), these
types will be considered as primitive at the level of the evaluator.

@< Other primitive type tags @>=
vector_type, matrix_type, @[@]

@~In choosing the short names below, we choose to hide the specific
mathematical meaning that was implied in the names (like |Weight|) these types
have in the Atlas software. We believe that a more extensive name might be
more confusing that helpful to users; besides, the interpretation of the
values is not entirely fixed (vectors are used for coweights and (co)roots as
well as for weights, and matrices could denote either a basis or an
automorphism of a lattice.

@< Other primitive type names @>=
"vec", "mat", @[@]

@ Here are the corresponding types derived from |value_base|. In the
constructor for |vector_value|, the argument is a reference to
|latticetypes::CoeffList|, which stands for
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
@)


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

@ To make a small but visible difference in printing between vectors and lists
of integers, weights will be printed in equal width fields one longer than the
minimum necessary.

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

@ For matrices we align columns, and print vertical bars along the sides.
@< Function def... @>=
void matrix_value::print(std::ostream& out) const
{ size_t k=val.numRows(),l=val.numColumns();
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

@*2 Implementing the conversion functions.
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

For |vector_conversion|, the |evaluate| method is easy, and follows a pattern
that we shall see a lot more: evaluate the subexpression, pop the result off
the |execution_stack| while converting it to an appropriate shared pointer,
then produce and push the resulting |value|. The conversion to the
appropriate kind of pointer is safe since we type-checked the expression, and
the shared pointer will clean up the used argument.

@< Function definition... @>=
latticetypes::Weight row_to_weight(const row_value& r)
{ latticetypes::Weight result(r.val.size());
  for(size_t i=0; i<r.val.size(); ++i)
    result[i]=force<int_value>(r.val[i].get())->val;
  return result;
}
@)
void vector_conversion::evaluate() const
{ exp->evaluate(); @+ shared_row r(get<row_value>());
  push_value(new vector_value(row_to_weight(*r)));
}

@ For |matrix_conversion|, the |evaluate| method is longer, but still
straightforward. Since the type of the argument was checked to be \.{[vec]},
we can safely cast the argument pointer to |(row_value*)|, and each of its
component pointers to |(vector_value*)|. Any ragged columns are silently
extended will null entries to make a rectangular shape for the matrix.

@< Function definitions @>=
void matrix_conversion::evaluate() const
{ exp->evaluate(); @+ shared_row r(get<row_value>());
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

@ The remaining ``internalising'' conversion function, from row of row of
integer to matrix, is very similar. Only this times the argument was checked
to be of type \.{[[int]]}, so we cast the raw component pointers to
|(row_value*)|, and then use |row_to_weight| to build column vectors to be
used by the matrix constructor.

@< Function definitions @>=
void matrix2_conversion::evaluate() const
{ exp->evaluate(); @+ shared_row r(get<row_value>());
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

@ The other three conversion functions perform the inverse transformations of
those given above. Again it will be handy to have a basic function
|weight_to_row| that performs more or less the inverse transformation of
|row_to_weight|, but rather than returning a |row_value| it returns a
|row_ptr| pointing to it.

@< Function definitions @>=
row_ptr weight_to_row(const latticetypes::Weight& v)
{ row_ptr result (new row_value(v.size()));
  for(size_t i=0; i<v.size(); ++i)
    result->val[i]=shared_value(new int_value(v[i]));
  return result;
}
@)
void int_list_conversion::evaluate() const
{ exp->evaluate(); @+ shared_vector v(get<vector_value>());
  push_value(weight_to_row(v->val));
}
@)
void vec_list_conversion::evaluate() const
{ exp->evaluate(); @+ shared_matrix m(get<matrix_value>());
  row_ptr result(new row_value(m->val.numColumns()));
  for(size_t i=0; i<m->val.numColumns(); ++i)
    result->val[i]=shared_value(new vector_value(m->val.column(i)));
  push_value(result);
}
@)
void int_list_list_conversion::evaluate() const
{ exp->evaluate(); @+ shared_matrix m(get<matrix_value>());
  row_ptr result(new row_value(m->val.numColumns()));
  for(size_t i=0; i<m->val.numColumns(); ++i)
    result->val[i]=shared_value(weight_to_row(m->val.column(i)).release());

  push_value(result);
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
  virtual void evaluate() const;
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


@ For type checking, we shall need a tuple pattern with any number of unknown
components.

@< Declarations of local functions @>=
type_ptr unknown_tuple(size_t n);

@~This pattern is built up in a simple loop.

@< Local function definitions @>=
type_ptr unknown_tuple(size_t n)
{ type_list_ptr tl(NULL);
  while (n-->0) tl=make_type_list(copy(unknown_type),tl);
  return make_tuple_type(tl);
}


@ When converting a tuple expression, we first try to specialise an unknown
type to a tuple with the right number of components; unless the type was
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
Evaluating a tuple display is much like evaluating a list display. Since we
use dynamically typed values internally, we can collect the components of a
tuple in a vector without problem. In fact we could reuse the type |row_value|
to hold the components of a tuple, if it weren't for the fact that it would
then print with brackets. Therefore we trivially derive a new class from
|row_value|.

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

@ Here are two functions that pack and unpack tuples from values on the stack;
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
{ shared_tuple tuple(get<tuple_value>());
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

@ The evaluation of a tuple display evaluates the components in a simple loop
and wraps the result in a |tuple_value|, which is accomplished by
|wrap_tuple|.

@< Function def... @>=
void tuple_expression::evaluate() const
{ for (size_t i=0; i<component.size(); ++i) component[i]->evaluate();
  wrap_tuple(component.size());
}

@*1 Function calls.
We shall now extend the evaluator as described until now to allow function
calls. This will turn out to involve several new notions as well.

We start with introducing a type for representing function calls after type
checking.

@< Type definitions @>=
struct fixed_call_expression : public expression_base
{ Hash_table::id_type function;
  expression argument;
@)
  explicit fixed_call_expression(Hash_table::id_type f,expression_ptr a)
   : function(f),argument(a.release()) @+{}
  virtual ~fixed_call_expression() @+ {@; delete argument; }
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};

@ To print a function call we print the function name and the argument, the
latter enclosed in parentheses unless it is a tuple expression. The latter
condition must be tested by a dynamic cast.

@< Function definitions @>=
void fixed_call_expression::print(std::ostream& out) const
{ out << main_hash_table->name_of(function);
  if (dynamic_cast<tuple_expression*>(argument)!=NULL) out << *argument;
  else out << '(' << *argument << ')';
}

@*2 Identifier tables.
%
As we said above, we need an identifier table to record the types of known
functions (and later other types of identifiers). For the moment we use a
single flat table; there will doubtlessly be need for handling scopes later.
We shall actually store both a type and a value in the table. Each such pair
form an |id_data| structure; it is intended only to reside inside a table. The
entry has shared ownership of the value, and the containing table will have
strict ownership of the type. Giving ownership of the type directly to
|id_data| would complicate its duplication, and therefore insertion into the
table. As a consequence we only allow construction in an empty state; the
pointers should be set only after insertion into the table, which then assumes
their ownership for the type.

@< Type definitions @>=

struct id_data
{ shared_value val; @+ type_declarator* type;
  id_data() : val(),type(NULL)@+ {}
};

@ We cannot store auto pointers in a table, so upon entering into the table we
will convert to ordinary pointers, transferring ownership to the table; this
holds for the values stored as well as for their types. The pointers returned
from a table lookup are also ordinary; destruction of the type declarators
referred to only takes place when the table itself is destructed.

@< Includes... @>=
#include <map>
#include "lexer.h"

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
  type_declarator* type_of(Hash_table::id_type id) const; // lookup
  shared_value value_of(Hash_table::id_type id) const; // lookup
@)
  size_t size() const @+{@; return table.size(); }
  void print(std::ostream&) const;
};

@ As was indicated above, the table has strict ownership of the contained
types, so the destructor must explicitly delete them.

@< Function def... @>=
Id_table::~Id_table()
{ for (map_type::const_iterator p=table.begin(); p!=table.end(); ++p)
  @/{@; delete p->second.type; }
}

@ The method |add| tries to insert the new mapping from the key |id| to a new
value-type pair. Doing so, it must distinguish the case that the key |id| is
already present. In that case the |insert| method will not overwrite the old
value, making no change instead, which is fortunate since it gives us the
occasion to clean up the old values; it then does return both a pointer
(iterator) to the obstructing (key,data) pair, and a boolean failure status.
When detecting that status, we perform the clean-up, and finally (in both
cases) we overwrite pointers to |val| and |type| into the node in the table.

@< Function... @>=
void Id_table::add(Hash_table::id_type id, shared_value val, type_ptr type)
{ std::pair<map_type::iterator,bool> trial
     =table.insert(std::make_pair(id,id_data()));
  if (not trial.second) // then key was present; destroy associated type
  {@; delete trial.first->second.type; }

  trial.first->second.val= shared_value(val);
  trial.first->second.type=type.release();
}

@ In order to have |const| lookup methods, we must refrain from inserting into
the table if the key is not found; we return a null pointer in that case. The
pointer returned by |type_of| remains owned by the table.

@< Function... @>=
type_declarator* Id_table::type_of(Hash_table::id_type id) const
{@; map_type::const_iterator p=table.find(id);
  return p==table.end() ? NULL : p->second.type;
}
shared_value Id_table::value_of(Hash_table::id_type id) const
{@; map_type::const_iterator p=table.find(id);
  return p==table.end() ? shared_value(value(NULL)) : p->second.val;
}

@ We provide a |print| member that shows the contents of the entire table.
@< Function... @>=

void Id_table::print(std::ostream& out) const
{ for (map_type::const_iterator p=table.begin(); p!=table.end(); ++p)
    out << main_hash_table->name_of(p->first) << ": " @|
        << *p->second.type << ": " << *p->second.val << std::endl;
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

@*2 Type-checking function calls.
%
The function type stored in the table prescribes the required type for the
argument (which could be a tuple type representing multiple arguments), and
also gives the final type the expression will have. Note that since we do
not own the types coming from the table, we must copy them if we want to
export them, either as type of the call or in the error value thrown in case
of a type mismatch.

In the function |convert_expr| we first get the function from the identifier
table, test if it is known and of function type, then type-check and convert
the argument expression, and build a converted function call~|call|. Finally
we test if the required type matches the return type (in which case we simply
return~|call|), or if the return type can be coerced to it (in which case we
return |call| as transformed by |coerce|); if neither is possible we throw
a~|type_error|.

@< Other cases for type-checking and converting... @>=
case function_call:
{ type_declarator* f_type=global_id_table->type_of(e.e.call_variant->fun);
  if (f_type==NULL || f_type->kind!=function_type)
    throw program_error("Call of "
                        + std::string(f_type==NULL ? "unknown " : "non-" )
			+"function '"
			+ main_hash_table->name_of(e.e.call_variant->fun)
		        +"'");
@.Call of unknown function@>
@.Call of non-function@>
  expression_ptr call (new fixed_call_expression
     (e.e.call_variant->fun,
      expression_ptr(convert_expr(e.e.call_variant->arg,
                                  f_type->func->arg_type))));
  if (type.specialise(f_type->func->result_type) or
      coerce(f_type->func->result_type,type,call))
    return call.release();
  else throw type_error(e,copy(f_type->func->result_type),copy(type));
}

@*2 Evaluating function calls.
%
The evaluation of the call of a built-in function executes is a ``wrapper
function'', that usually consists of a call to a library function sandwiched
between unpacking and repacking statements; in some simple cases a wrapper
function may decide to do the entire job itself.

The arguments and results of wrapper functions will be transferred from and to
stack as a |value|, so a wrapper function has neither arguments nor a result
type. Thus variables that refer to a wrapper function have the type
|wrapper_function| defined below. We shall need to bind values of this type to
identifiers representing built-in functions, so we derive an associated
``primitive type'' from |value_base|.

@< Type definitions @>=
typedef void (* wrapper_function)();
@)
struct builtin_value : public value_base
{ wrapper_function val;
@)
  std::string print_name;
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


@ To evaluate a |fixed_call_expression| object we evaluate the arguments, get
and check the wrapper function, and call it. Since values are handled via the
execution stack, we don't see them at all in this code.

@< Function definitions @>=
void fixed_call_expression::evaluate() const
{ argument->evaluate(); // push evaluated argument on stack
  value f_val=global_id_table->value_of(function).get();
    // ownership assured by table
  if (f_val==NULL) throw std::logic_error("built-in function absent");
@.built-in function absent@>
  force<builtin_value>(f_val)->val();
  // call the wrapper function, leaving result on the stack
}

@*1 Invoking the type checker.
%
Let us recapitulate what will happen. The parser will read what the user
types, and returns an |expr| value. Then we shall call |convert_expr| for this
value, and if this does not throw any exceptions, we are ready to call
|evaluate| (for the |expr| value that may have been modified by the insertion
of conversion functions) of the |evaluate| method of the |expression| value
returned by |convert_expr|; after this main program will print the result. The
call to |convert_expr| is done via |analyse_types|.

@< Declarations of exported functions @>=
type_ptr analyse_types(const expr& e,expression_ptr& p)
   throw(std::bad_alloc,std::runtime_error);

@~An important reason for having this function is that it will give us an
occasion to catch any thrown |type_error| and |program_error| exceptions,
something we did not want to do inside the (implicitly) recursive function
|convert_expr|; since we cannot return normally from |analyse_types| in the
presence of these errors, we map these errors to |std::runtime_error|, an
exception for which the code that calls us will have to provide a handler
anyway.

@< Function definitions @>=
type_ptr analyse_types(const expr& e,expression_ptr& p)
  throw(std::bad_alloc,std::runtime_error)
{ try
  {@; type_ptr type=copy(unknown_type);
    p.reset(convert_expr(e,*type));
    return type;
  }
  catch (type_error& err)
  { std::cerr << err.what() << std::endl <<
    "Subexpression " << err.offender << @|
    " has wrong type: found " << *err.actual << @|
    " while " << *err.required << " was needed.\n";
@.Subexpression has wrong type@>
  }
  catch (program_error& err)
  {@; std::cerr << err.what() <<
          " in expression '" << e << "'\n";
  }
  throw std::runtime_error("Type check failed");
@.Type check failed@>
}

@*1 Wrapper functions.
%
We have not defined any wrapper functions yet, and therefore have nothing in
the |global_id_table|. Our first wrapper function is a trivial one, but which
will be put into the |global_id_table| under different names and signatures,
each forcing a particular argument type.

@< Function definitions @>=

void id_wrapper () @+{} // nothing to do, value stays on the stack

@< Definition of other wrapper functions @>@;

@ All that remains to do for the moment is to install these functions into the
|global_id_table|.

@< Declarations of exported functions @>=
void initialise_evaluator();

@ The following function greatly facilitates this repetitive task.

@< Template and inline... @>=
inline void install_function
 (wrapper_function f,const char*name, const char* type)
{ shared_value val(new builtin_value(f,name));
  global_id_table->add
    (main_hash_table->match_literal(name),val,make_type(type));
}

@ Here we install the trivial wrapper function under two different identity
signatures.

@< Function definitions @>=
void initialise_evaluator()
{ execution_stack.reserve(16); // avoid some early reallocations
@)
  install_function(id_wrapper,"vec","(vec->vec)");
  install_function(id_wrapper,"mat","(mat->mat)");
@/@< Installation of other built-in functions @>
}

@ We shall now introduce some real built-in functions, starting with integer
arithmetic. Arithmetic operators are implemented by wrapper functions with two
integer arguments (or one in the case of unary minus). Note that the values
are pulled from the stack in reverse order, which is important for the
non-commutative operations like `|-|', `|/|' and~`|%|'. Since values are
shared, we must allocate new value objects for the results.

@< Definition of other wrapper functions @>=
void plus_wrapper ()
{ push_tuple_components();
  shared_int j(get<int_value>()); shared_int i(get<int_value>());
  push_value(new int_value(i->val+j->val));
}
@)
void minus_wrapper ()
{ push_tuple_components();
  shared_int j(get<int_value>()); shared_int i(get<int_value>());
  push_value(new int_value(i->val-j->val));
}
@)
void times_wrapper ()
{ push_tuple_components();
  shared_int j(get<int_value>()); shared_int i(get<int_value>());
  push_value(new int_value(i->val*j->val));
}
@)
void divide_wrapper ()
{ push_tuple_components();
  shared_int j(get<int_value>()); shared_int i(get<int_value>());
  if (j->val==0) throw std::runtime_error("Division by zero");
  push_value(new int_value(i->val/j->val));
}
@)
void modulo_wrapper ()
{ push_tuple_components();
  shared_int  j(get<int_value>()); shared_int i(get<int_value>());
  if (j->val==0) throw std::runtime_error("Modulo zero");
  push_value(new int_value(i->val%j->val));
}
@)
void unary_minus_wrapper ()
{@; shared_int i=get<int_value>();  push_value(new int_value(-i->val)); }
@)
void divmod_wrapper ()
{ push_tuple_components();
  shared_int j(get<int_value>()); shared_int i(get<int_value>());
  if (j->val==0) throw std::runtime_error("DivMod by zero");
  push_value(new int_value(i->val%j->val));
  push_value(new int_value(i->val/j->val));
  wrap_tuple(2);
}

@ Here is a simple function that outputs a string without its quotes, but with
a terminating newline. This is the first place in this file where we produce
output to a file. In general, rather than writing directly to |std::cout|, we
shall pass via a pointer whose |output_stream| value is maintained in the main
program, so that redirecting output to a different stream can be easily
implemented. Since this is a wrapper function (hence without arguments) there
is no other way to convey the output stream to be used than via a dedicated
global variable.

@< Definition of other wrapper functions @>=
void print_wrapper ()
{ shared_string s(get<string_value>());
  *output_stream << s->val << std::endl;
  wrap_tuple(0); // don't forget to return a value
}


@ We now define a few functions, to really exercise something, even if it is
modest, from the Atlas library. The remainder of this chapter is not really to
be considered part of the interpreter, but a first step to its interface with
the Atlas library, which is developed in much more detail in the compilation
unit \.{built-in-types}. But it was important to have some experience with
calling function from the library, whence these functions. First of all the
identity matrix and matrix transposition.

@< Declarations of exported functions @>=
void id_mat_wrapper ();
void transpose_mat_wrapper ();

@ Since in general built-in functions may throw exceptions (even |transpose|!)
we hold the pointers to the values created to hold their results in
auto-pointers. The reason that in |id_mat_wrapper| we create a |matrix_value|
around an empty |LatticeMatrix| rather than build a filled matrix object
first, is that the constructor for |matrix_value| would than have to copy that
matrix object (which implies copying its contents). Actually this is a case
where the compiler might avoid such a copy if the identity matrix were
produced directly in the argument to |new matrix_value| (rather than in a
variable), but there is (currently) no constructor or function that produces
an identity matrix as result. We also define |diagonal_wrapper|, a slight
generalisation of |id_mat_wrapper| that produces a diagonal matrix from a
vector.

@< Definition of other wrapper functions @>=
void id_mat_wrapper ()
{ shared_int i(get<int_value>());
  matrix_ptr m
     (new matrix_value(latticetypes::LatticeMatrix()));
  matrix::identityMatrix(m->val,std::abs(i->val)); push_value(m);
}
@)
void transpose_mat_wrapper ()
{ shared_matrix m(get<matrix_value>());
  push_value(new matrix_value(m->val.transposed()));
}
@)
void diagonal_wrapper ()
{ shared_vector d(get<vector_value>());
  size_t n=d->val.size();
  matrix_ptr m
     (new matrix_value(latticetypes::LatticeMatrix(n,n,0)));
  for (size_t i=0; i<n; ++i) m->val(i,i)=d->val[i];
  push_value(m);
}

@ Now the product of a matrix and a vector or matrix. The first of these was
in fact our first function with more than one argument (arithmetic on integer
constants was done inside the parser at that time). We make them callable from
other compilation units.

@< Declarations of exported functions @>=
void mv_prod_wrapper ();
void mm_prod_wrapper ();

@ In |mv_prod_wrapper| we use the matrix method |apply| which requires its
output vector to be already of the proper size, but which nevertheless could
throw an exception (since it copies its input argument, for in case the output
argument should coincide with it). Therefore we use an auto-pointer for the
output vector~|w|. In |mm_prod_wrapper|, the method |operator*=| does its own
resizing (and may also throw an exception).

For wrapper functions with multiple arguments, we must always remember that
they are to be popped from the stack in reverse order, but in the case of
|mm_prod_wrapper| this is particularly important.

@< Definition of other wrapper functions @>=
void mv_prod_wrapper ()
{ push_tuple_components();
  shared_vector v(get<vector_value>());
  shared_matrix m(get<matrix_value>());
  if (m->val.numColumns()!=v->val.size())
    throw std::runtime_error(std::string("Size mismatch ")@|
     + num(m->val.numColumns()) + ":" + num(v->val.size()) + " in mv_prod");
  push_value(new vector_value(m->val.apply(v->val)));
}
@)
void mm_prod_wrapper ()
{ push_tuple_components();
  shared_matrix r(get<matrix_value>()); // right factor
  shared_matrix l(get<matrix_value>()); // left factor
  if (l->val.numColumns()!=r->val.numRows())
  { std::ostringstream s;
    s<< "Size mismatch " << l->val.numColumns() << ":" << r->val.numRows()
    @| << " in mm_prod";
    throw std::runtime_error(s.str());
  }
  push_value(new matrix_value(l->val*r->val));
}

@ As a last example, here is the Smith normal form algorithm. We provide both
the invariant factors and the rewritten basis on which the normal for is
assumed, as separate functions, and to illustrate the possibilities of tuples,
the two combined into a single function.

@h "smithnormal.h"

@< Definition of other wrapper functions @>=
void invfact_wrapper ()
{ shared_matrix m(get<matrix_value>());
  size_t nr=m->val.numRows();
  latticetypes::WeightList b; @+ matrix::initBasis(b,nr);
  vector_ptr inv_factors
     @| (new vector_value(latticetypes::Weight(0)));
  smithnormal::smithNormal(inv_factors->val,b.begin(),m->val);
  push_value(inv_factors);
}
@)
void Smith_basis_wrapper ()
{ shared_matrix m(get<matrix_value>());
  size_t nr=m->val.numRows();
  latticetypes::WeightList b; @+ matrix::initBasis(b,nr);
  latticetypes::Weight inv_factors(0);
  smithnormal::smithNormal(inv_factors,b.begin(),m->val);
  push_value(new matrix_value(latticetypes::LatticeMatrix(b)));
}
@)

void Smith_wrapper ()
{ shared_matrix m(get<matrix_value>());
  size_t nr=m->val.numRows();
  latticetypes::WeightList b; @+ matrix::initBasis(b,nr);
  vector_ptr inv_factors
     @| (new vector_value(latticetypes::Weight(0)));
  smithnormal::smithNormal(inv_factors->val,b.begin(),m->val);
@/push_value(new matrix_value(latticetypes::LatticeMatrix(b)));
  push_value(inv_factors);
  wrap_tuple(2);
}

@ Here is one more wrapper function that uses the Smith normal form algorithm,
but behind the scenes, namely to invert a matrix. Since this cannot be done in
general over the integers, we return an integral matrix and a common
denominator to be applied to all coefficients.
@< Definition of other wrapper functions @>=
void invert_wrapper ()
{ shared_matrix m(get<matrix_value>());
  if (m->val.numRows()!=m->val.numColumns())
  { std::ostringstream s;
    s<< "Cannot invert a " @|
     << m->val.numRows() << "x" << m->val.numColumns() << " matrix";
    throw std::runtime_error(s.str());
  }
  int_ptr denom(new int_value(0));
@/push_value(new matrix_value(m->val.inverse(denom->val)));
  push_value(denom);
  wrap_tuple(2);
}

@ We must not forget to install what we have defined. The names of the
arithmetic operators correspond to the ones used in the parser definition
file \.{parser.y}.

@< Installation of other built-in functions @>=
install_function(plus_wrapper,"+","(int,int->int)");
install_function(minus_wrapper,"-","(int,int->int)");
install_function(times_wrapper,"*","(int,int->int)");
install_function(divide_wrapper,"/","(int,int->int)");
install_function(modulo_wrapper,"%","(int,int->int)");
install_function(unary_minus_wrapper,"-u","(int->int)");
install_function(divmod_wrapper,"/%","(int,int->int,int)");
install_function(print_wrapper,"print","(string->)");
install_function(id_mat_wrapper,"id_mat","(int->mat)");
install_function(transpose_mat_wrapper,"transpose_mat","(mat->mat)");
install_function(diagonal_wrapper,"diagonal_mat","(vec->mat)");
install_function(mv_prod_wrapper,"mv_prod","(mat,vec->vec)");
install_function(mm_prod_wrapper,"mm_prod","(mat,mat->mat)");
install_function(invfact_wrapper,"inv_fact","(mat->vec)");
install_function(Smith_basis_wrapper,"Smith_basis","(mat->mat)");
install_function(Smith_wrapper,"Smith","(mat->mat,vec)");
install_function(invert_wrapper,"invert","(mat->mat,int)");

@*1 Applied identifiers.
%
We can use identifiers, which come in two kinds: local and global ones. The
two cases are treated in fairly different way, and during type analysis the
two are converted into different kinds of |expression|. For the former, we
shall use a stack of variable bindings, which is independent of the
|execution_stack| (that is used for anonymous components of expressions being
evaluated); it is handled as a linked list, and accessed through a (smart)
pointer that is local to the evaluator.

@< Local var... @>=
context_ptr execution_context;

@ Whan an exception is caught at runtime, the execution context must be reset;
being a smart pointer, this takes care of adjusting reference counts and maybe
deleting values that are part of the discarded context.

@< Actions... @>=
execution_context.reset();

@ We start considering
global identifiers, which will be converted into a |global_identifier| object.

@< Type definitions @>=
struct global_identifier : public expression_base
{ Hash_table::id_type code;
@)
  explicit global_identifier(Hash_table::id_type id) : code(id) @+{}
  virtual ~global_identifier() @+ {}
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};

@ To print a global identifier, we get its name from the main hash table.

@< Function definitions @>=
void global_identifier::print(std::ostream& out) const
{@; out<< main_hash_table->name_of(code); }

@ The method |global_identifier::evaluate| just picks the value from the
identifier table, which creates a pointer sharing the value which remains in
the table. Sharing avoids calling the |clone| method here, which was used in a
first implementation (and was introduced for that purpose), but we retain that
method for future use elsewhere. The error thrown here is a |logic_error|,
since a missing identifier should have been caught during the type check.

@< Function definitions @>=
void global_identifier::evaluate() const
{ shared_value p=global_id_table->value_of(code);
  if (p==NULL) throw std::logic_error
  @|   (std::string("Identifier without value:")
        +main_hash_table->name_of(code));
@.Identifier without value@>
  push_value(p);
}

@ Local identifiers are different in that they need to look up their value in
the current execution context.

@< Type definitions @>=
class local_identifier : public global_identifier
{ size_t depth, offset;
public:
  explicit local_identifier(Hash_table::id_type id, size_t i, size_t j)
     : global_identifier(id), depth(i), offset(j) @+{}
  virtual void evaluate() const; // only this method is redefined
};

@ The method |local_identifier::evaluate| just looks up a value in the
|execution_context|.

@< Function definitions @>=
void local_identifier::evaluate() const @+
{@; push_value(execution_context->elem(depth,offset));}

@ For an applied identifier, we first look in |id_context| for a binding of
the identifier, and if found it will be a local identifier, and otherwise we
look in |global_id_table|. If found in either way, the associated type must
equal the expected type (if any), or be convertible to it using |coerce|.

@< Other cases for type-checking and converting... @>=
case applied_identifier:
{ type_declarator* it=NULL; expression_ptr id; size_t i=0;
  for (bindings* p=id_context; p!=NULL; p=p->next,++i)
    for (size_t j=0; j<p->size(); ++j)
      if ((*p)[j].first==e.e.identifier_variant)
      { it=(*p)[j].second;
        id.reset(new local_identifier(e.e.identifier_variant,i,j));
        goto common_ending;
      }

  it = global_id_table->type_of(e.e.identifier_variant);
  if (it==NULL) throw program_error
  @|   (std::string("Undefined identifier ")
	+main_hash_table->name_of(e.e.identifier_variant));
@.Undefined identifier@>
  id.reset(new global_identifier(e.e.identifier_variant));

common_ending:
  if (type.specialise(*it) or coerce(*it,type,id))
    return id.release();
  else throw type_error(e,copy(*it),copy(type));
}

@*1 Let expressions.
%
We shall now consider a simple type of expression in which local variables
occur, the let-expression. A let expression is equivalent to an anonymous
function ($\lambda$-expression) applied to the expression(s) in the
let-declarations. There is however a simplification, in that no types have to
be declared for the function parameters, since the types of the corresponding
expressions can be used for this. Nevertheless, we shall in converting
expressions to internal form forget the syntactic origin of the expression,
and translate to an application of an anonymous function, giving us an
occasion to introduce such functions as new form of |expression|, and
expressions in which such functions rather than just built-in ones can be
called.

@< Type def... @>=
struct call_expression : public expression_base
{ expression function, argument;
@)
  explicit call_expression(expression_ptr f,expression_ptr a)
   : function(f.release()),argument(a.release()) @+{}
  virtual ~call_expression() @+ {@; delete function; delete argument; }
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};
@)
struct lambda_expression : public expression_base
{ std::vector<Hash_table::id_type> param;
  shared_expression body;
@)
  lambda_expression(const std::vector<Hash_table::id_type>& p, expression_ptr b)
  : param(p) , body(b.release()) @+{}
  virtual ~lambda_expression() @+{}
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};


@ To print a function call we print the function expression in parentheses,
and the argument, latter enclosed in parentheses unless it is a tuple
expression (which already has parentheses), which condition is tested by a
dynamic cast. To print an anonymous function, we print the parameter list
enclosed in parentheses, followed by a colon and by the function body. The
parameter list cannot include types with the current setup, as they are not
explicitly stored after type analysis.

@< Function definitions @>=
void call_expression::print(std::ostream& out) const
{ out << '(' << *function << ')';
  if (dynamic_cast<tuple_expression*>(argument)!=NULL) out << *argument;
  else out << '(' << *argument << ')';
}
@)
void lambda_expression::print(std::ostream& out) const
{ out << '(' << main_hash_table->name_of(param[0]);
  for (size_t i=1; i<param.size(); ++i)
    out << ',' << main_hash_table->name_of(param[i]);
  out << "): " << *body;
}

@ In converting let-expressions we must distinguish those with a simgle
declaration and those with multiple declarations, since the latter need to
construct a tuple from the various right hand sides.

@< Other cases for type-checking and converting... @>=
case let_expr:
{ let lexp=e.e.let_variant;
  if (lexp->first.next==NULL)
    @< Convert let-expression with a single declaration @>
  else
     @< Convert let-expression with multiple declarations @>

}

@ We proceed by first deducing the type of the declared identifier from the
right hand side of its declaration, then setting up a new binding for the
identifier with the type found, and converting the body to the required type
in the extended context. Note that the single-variable |bindings| remains a
local variable, but it is temporarily preprended to the context by calling
|push| and |pop| with the pointer variable |id_context|. The |expression|
obtained from converting the body is first turned into a |lambda_expression|,
and then an application of that expression to the argument is produced and
returned.

@< Convert let-expression with a single declaration @>=
{
  type_ptr decl_type(new type_declarator);
  expression_ptr arg(convert_expr(lexp->first.val,*decl_type));
@/bindings new_binding(1);
  new_binding.add(lexp->first.id,decl_type);
  new_binding.push(id_context);
  expression_ptr body(convert_expr(lexp->body,type));
  new_binding.pop(id_context);
  std::vector<Hash_table::id_type> formal(1,lexp->first.id);
  expression_ptr func(new lambda_expression(formal,body));
  return new call_expression(func,arg);
}

@ In case of multiple declarations we first count how many there are. Then we
proceed mostly like befor, but repeat the action for each declaration, and we
assemble a list of formal parameters, a list of converted arguments, and a
list of bindings of identifiers to types wihich is to be pushed onto the
identifier context. The remainder is almost identical, except that a cast is
needed to change the argument tuple auto-pointer into an |expression_ptr|.

@< Convert let-expression with multiple declarations @>=
{ size_t l=1; // start counting declarations
  for (struct let_node* p=lexp->first.next; p!=NULL; p=p->next)
    ++l;
@)
  bindings new_bindings(l);
  std::auto_ptr<tuple_expression> arg(new tuple_expression(l));
  std::vector<Hash_table::id_type> formal; formal.reserve(l);
@/struct let_node* p=&lexp->first; size_t i=0;
  do
  {
    type_ptr decl_type(new type_declarator);
    arg->component[i]=convert_expr(p->val,*decl_type);
  @/formal.push_back(p->id);
    new_bindings.add(p->id,decl_type);
  }
  while (++i,(p=p->next)!=NULL);
@)
  new_bindings.push(id_context); // modify context only now!
  expression_ptr body(convert_expr(lexp->body,type));
  new_bindings.pop(id_context);
  expression_ptr func(new lambda_expression(formal,body));
  return new call_expression(func,expression_ptr(arg));
}

@ Before we can dicuss the evaluation of let-expressions, we need to introduce
a type for the intermediate value produced by the anonymous function, before
it is actually called. Such a value is traditionally called a closure, and it
contains (a reference to) the expression body, as well as the evaluation
context current at the point the anonymous function is produced. (For a
let-expression, the closure will be immediately called the same context, but
we set up things so that they will continue to work in a more general
setting.)

@< Type def... @>=
struct closure_value : public value_base
{ context_ptr cont;
  std::vector<Hash_table::id_type> param;
   // used in evaluation only to count arguments
  shared_expression body;
@)
  closure_value@|(context_ptr c,
                  const std::vector<Hash_table::id_type>& p,
                  shared_expression b) : cont(c), param(p), body(b) @+{}
  void print(std::ostream& out) const;
  closure_value* clone() const @+
  {@; return new closure_value(cont,param,body); }
  static const char* name() @+{@; return "closure"; }
};
typedef std::auto_ptr<closure_value> closure_ptr;
typedef std::tr1::shared_ptr<closure_value> shared_closure;

@ For now a closure prints just the body, preceded by formal parameter list.
One could imagine printing after this body ``where'' followed by the bindings
held in the |context| field. Even better only the bindings for relevant
(because referenced) identifiers could be printed.

@< Function def... @>=
void closure_value::print(std::ostream& out) const
{ for (size_t i=0; i<param.size(); ++i)
    out << (i==0 ? '(' : ',') << param[i];
  out << "): " << *body;
}

@ The following code defines in essence a call-by-value $\lambda$ calculus
evaluator. Later one should extend the method |call_expression::evaluate| to
allow for the possibility the the function part evaluates to a built-in
function, but currently that cannot happen, since |call_expression| values are
only produced from let-expressions.

@< Function def... @>=
void lambda_expression::evaluate() const
{@; closure_ptr result(new closure_value(execution_context,param,body));
  push_value(result);
}
@)
void call_expression::evaluate() const
{ function->evaluate(); @+
  shared_closure f(get<closure_value>());
@)argument->evaluate();
  context_ptr saved_context(execution_context);
  if (f->param.size()==1)
  {
    shared_value arg= pop_value();
    execution_context.reset(new context
      (f->cont,std::vector<shared_value>(1,arg)));
  }
  else
  {
    shared_tuple args(get<tuple_value>());
    execution_context.reset(new context(f->cont,args->val));
  }
  f->body->evaluate();
  execution_context = saved_context;
}

@*1 Array subscription.
%
While we have seen expressions to build list, and vectors and matrices out of
them, we so far are not able to access their components once they are
constructed. To that end we shall now introduce operations to index such
values. We allow subscription of rows, but also of vectors and matrices. Since
after type analysis we know which of the cases applies, we define several
classes derived from |expression_base|. These types differ mostly by their
|evaluate| method, so after the first one we derive the others from it.

@< Type definitions @>=
struct row_subscription : expression_base
{ expression array, index;
@)
  row_subscription(expression_ptr a, expression_ptr i)
  : array(a.release()),index(i.release()) @+{}
  virtual ~row_subscription() @+ {@; delete array; delete index; }
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};
@)
struct vector_subscription : row_subscription
{ vector_subscription(expression_ptr a, expression_ptr i)
  : row_subscription(a,i) @+{}
  virtual void evaluate() const;
};
@)
struct matrix_subscription : row_subscription
{ matrix_subscription(expression_ptr a, expression_ptr ij)
  : row_subscription(a,ij) @+{}
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};

@ These subscriptions are printed in the usual subscription syntax. For matrix
subscriptions, where the index type is \.{(int,int)}, the index expression is
quite likely to be a tuple display, in which case we suppress parentheses.
Since we have passed the type check here, we know that any tuple display is
necessarily a pair.

@< Function definitions @>=
void row_subscription::print(std::ostream& out) const
{@; out << *array << '[' << *index << ']';
}
@)
void matrix_subscription::print(std::ostream& out) const
{ tuple_expression* p=dynamic_cast<tuple_expression*>(index);
  if (p==NULL) out << *array << '[' << *index << ']';
  else
    out << *array << '[' << *p->component[0] << ',' << *p->component[1] << ']';
}

@ When encountering a subscription in |convert_expr|, we determine the types
of array and of the indexing expression separately, ignoring so far any type
required by the context. Then we look if the types agree with any of the three
types of subscription expressions that we can convert to, throwing an error if
none does. Finally we check is the a priori type |subscr_type| of the
subscripted expression equals of specialises the required |type|, or can be
converted to it by |coerce|, again throwing an error if nothing works. For the
indexing expression only equality of types is admitted, since nothing can be
coerced to \.{int} or to \.{(int,int)} anyway, and there is little point in
catering for (indexing) expressions having completely undetermined type (which
can only happen for expressions that cannot be evaluated without error).

@< Other cases for type-checking and converting... @>=
case subscription:
{ type_declarator array_type, index_type, subscr_type;
    // initialised to |undetermined_type|
  expression_ptr array
    (convert_expr(e.e.subscription_variant->array,array_type));
  expression_ptr index
    (convert_expr(e.e.subscription_variant->index,index_type));
  expression_ptr subscr;
  if (array_type.kind==row_type) // a row subscription
    if (index_type!=int_type) throw type_error
      (e.e.subscription_variant->index,copy(index_type),copy(int_type));
    else
    { subscr_type.specialise(*array_type.component_type); // type after indexing
      subscr.reset(new row_subscription(array,index));
    }
  else if (array_type==vec_type)
    if (index_type!=int_type) throw type_error
      (e.e.subscription_variant->index,copy(index_type),copy(int_type));
    else
    @/{@; subscr_type.specialise(int_type);
      subscr.reset(new vector_subscription(array,index));
    }
  else if (array_type==mat_type)
    if (index_type!=int_int_type) throw type_error
      (e.e.subscription_variant->index,copy(index_type),copy(int_int_type));
    else
    @/{@; subscr_type.specialise(int_type);
      subscr.reset(new matrix_subscription(array,index));
    }
  else throw type_error // array expression is not of any aggregate type
      (e.e.subscription_variant->array,copy(array_type),copy(row_of_type));
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
void row_subscription::evaluate() const
{ shared_int i((index->evaluate(),get<int_value>()));
  shared_row r((array->evaluate(),get<row_value>()));
  if (static_cast<unsigned int>(i->val)>=r->val.size())
    throw std::runtime_error("index "+num(i->val)+" out of range");
  push_value(r->val[i->val]);
}
@)
void vector_subscription::evaluate() const
{ shared_int i((index->evaluate(),get<int_value>()));
  shared_vector v((array->evaluate(),get<vector_value>()));
  if (static_cast<unsigned int>(i->val)>=v->val.size())
    throw std::runtime_error("index "+num(i->val)+" out of range");
  push_value(new int_value(v->val[i->val]));
}
@)
void matrix_subscription::evaluate() const
{ index->evaluate(); push_tuple_components();
  shared_int j(get<int_value>());
  shared_int i(get<int_value>());
  shared_matrix m((array->evaluate(),get<matrix_value>()));
  if (static_cast<unsigned int>(i->val)>=m->val.numRows())
    throw std::runtime_error("initial index "+num(i->val)+" out of range");
  if (static_cast<unsigned int>(j->val)>=m->val.numColumns())
    throw std::runtime_error("final index "+num(j->val)+" out of range");
  push_value(new int_value(m->val(i->val,j->val)));
}


@* Operations other than evaluation of expressions.
This section will be devoted to some other interactions between user and
program that do not consist just of evaluating expressions. What will
presented is not particularly related to the evaluator, and is present here
for somewhat opportunistic purposes.

Several kinds of things will be defined: there will be a tiny bit of global
state that can be set from the main program and inspected by any module that
cares to (and that reads \.{evaluator.h}); there are the declarations of
built-in types that are not dealt with directly by the evaluator, but which
the evaluator must know about if only to store the type names at the
appropriate places, and finally there are several functions that can be called
directly by the parser (those are not declared in our header file, but in
\Cee-style in \.{parsetree.h}, which the parser does read).

@< Declarations of global variables @>=
extern int verbosity;

@~By raising the value of |verbosity|, some trace of internal operations can
be activated.

@< Global variable definitions @>=
int verbosity=0;

@ In order for this compilation unit to function properly, it must know of the
existence and names for other built-in types. We could scoop up these names
using clever \&{\#include} directives, but that is not really worth the hassle
(when adding such types you have to recompile this unit anyway; it is not so
much worse to actually extend the lines below as well).

@< Other primitive type names @>=
"LieType","RootDatum", "InnerClass", "RealForm", "DualRealForm",
"CartanClass", @[@]

@~The enumeration values below are never directly used, but they must be
present so that the evaluator be aware of the number of primitive types (via
the final enumeration value |nr_of_primitive_types|).

@< Other primitive type tags @>=
complex_lie_type_type , root_datum_type, inner_class_type, real_form_type,
dual_real_form_type, Cartan_class_type, @[@]

@ It is useful to print type information, either for a single expression or for
all identifiers in the table. We declare the pointer that was already used in
|print_wrapper|.

@< Declarations of global variables @>=
extern std::ostream* output_stream;

@ The |output_stream| will normally point to |std::cout|.

@< Global variable definitions @>=
std::ostream* output_stream= &std::cout;

@ For the moment, applied identifiers can only get their value through the
function |global_set_identifier| that was declared in \.{parsetree.h}. We
define it here since it uses the services of the evaluator. It will be called
by the parser for global identifier definitions (simple or multiple);
therefore it has \Cee-linkage.

Recall that the parser guarantees that |ids| represents a list of identifiers.
What has to be done here is straightforward. We type-check the expression~|e|,
storing its result, then evaluate it, and store the (type,value) pair(s) into
in |global_id_table|. To provide some feedback to the user we report the type
assigned, but not the value since this might result in more output than is
desirable in case of an assignment.

A |std::runtime_error| may be thrown either during the initial analysis of
the assignment statement or during type check or evaluation, and we catch all
those cases here. Note however that no exception can be thrown between the
successful return from |evaluate| and a non-multiple assignment, so we can use
an ordinary pointer to hold the pointer returned. Whether or not an error is
caught, the identifier list |ids| and the expression |rhs| should not be
destroyed here, since the parser which aborts after calling this function
should do that while clearing its parsing stack.

@< Function definitions @>=
extern "C"
void global_set_identifier(expr_list ids, expr rhs)
{ using namespace atlas::interpreter; using namespace std;
  try
  { expression_ptr e;
    type_ptr t=analyse_types(rhs,e);
    if (ids->next!=NULL)
      @< Check that identifiers are distinct and that |t| is an appropriate
         tuple type; if not, |throw| a |runtime_error| @>
    e->evaluate();
    shared_value v=pop_value();
    if (ids->next==NULL)
    { cout << "Identifier " << ids->e << ": " << *t << std::endl;
      global_id_table->add(ids->e.e.identifier_variant,v,t);
    }
    else @< Perform a multiple assignment @>
  }
  catch (runtime_error& err)
  { cerr << err.what() << ", identifier" << (ids->next!=NULL ? "s " :" ");
    for (expr_list l=ids; l!=NULL; l=l->next)
      cerr << main_hash_table->name_of(l->e.e.identifier_variant)
           << (l->next!=NULL?",":"");
    cerr << " not assigned to.\n";
    reset_evaluator(); main_input_buffer->close_includes();
  }
  catch (logic_error& err)
  { cerr << "Unexpected error: " << err.what() << ", evaluation aborted.\n";
    reset_evaluator(); main_input_buffer->close_includes();
  }
  catch (exception& err)
  { cerr << err.what() << ", evaluation aborted.\n";
    reset_evaluator(); main_input_buffer->close_includes();
  }
}

@ Whenever we have an assignment to an expression list, the identifiers must
be distinct and the right hand side must have a tuple type with the proper
number of components; otherwise there are no restrictions. To test whether
identifiers are distinct, we put them into a |set| and check membership.

@h <set>
@< Check that identifiers are distinct and that |t| is an appropriate tuple
   type... @>=
{ if (t->kind!=tuple_type)
    throw runtime_error ("Multiple assignment requires tuple value");
@.Multiple assignment requires tuple@>
  set<Hash_table::id_type> seen;
  type_list tl=t->tuple;
  for (expr_list l=ids; l!=NULL||tl!=NULL; l=l->next,tl=tl->next)
  { if (l==NULL || tl==NULL)
      throw runtime_error
        ("Right hand side has too "+string(l==NULL ? "many":"few")
        +" components");
@.Right hand side has too many...@>
@.Right hand side has too few...@>
    if (!seen.insert(l->e.e.identifier_variant).second)
        // then identifier was already present
      throw runtime_error
        (std::string("Repeated identifier ")
        +main_hash_table->name_of(l->e.e.identifier_variant)
        +" in multiple assignment");
  }
}

@ In the multiple assignment case, the type check above should guarantee that
the result of evaluation is an appropriate tuple value. The multiple
assignment is performed by simultaneously traversing the components of the
identifier list, the tuple value and the tuple type; after each value is
stored into |global_id_table| the pointer to it is removed from the tuple
value immediately to ensure that it will not be deleted upon destruction of
the tuple value. We report the type of each variable assigned separately.


@< Perform a multiple assignment @>=
{ tuple_value* tv = force<tuple_value>(v.get());
  cout << "Identifiers ";
  size_t i=0; type_list tl=t->tuple;
  for (expr_list l=ids; l!=NULL; l=l->next,++i,tl=tl->next)
  { cout << l->e << ": " << tl->t << ( l->next!=NULL ? ", " : ".\n");
    shared_value val = tv->val[i];
    global_id_table->add(l->e.e.identifier_variant,val,copy(tl->t));
  }
}

@ The function |type_of_expr| prints the type of a single expression, without
evaluating it; since the expression is not limited to being a single
identifier, we must cater for the possibility that the type analysis fails, in
which case |analyse types| after catching it re-throws a |std::runtime_error|
(however neither a |std::logic_error| nor any more general |exception| can
be thrown while merely type checking).

@< Function definitions @>=
extern "C"
void type_of_expr(expr e)
{ try
  {@; expression_ptr p;
    *output_stream << "type: " << *analyse_types(e,p) << std::endl;
  }
  catch (std::runtime_error& err) {@; std::cerr<<err.what()<<std::endl; }
}

@ The function |show_ids| prints a table of all known identifiers and their
types.

@< Function definitions @>=
extern "C"
void show_ids()
{@; *output_stream << *global_id_table;
}


@* Index.

% Local IspellDict: british
