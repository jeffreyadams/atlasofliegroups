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
@c
namespace atlas { namespace interpreter {
@< Local type definitions @>@;
@< Declarations of local functions @>@;
@< Global variable definitions @>@;
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
consider ``type'' values, represented by the structure |type_declarer|, in
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

@*1 Types to represent types.
%
We make a fundamental choice to check types before attempting to execute any
expression entered by the use; thus type errors can be signalled earlier and
expressed more understandably than if would wait until an accident is about to
happen at runtime. This means that types are to be represented independently
of runtime values. Of course the interpreter executes the same code regardless
of the types of the values is manipulates; those values are then accessed via
generic pointers. The runtime values do also contain an indication of their
actual type, which will allow us to do some dynamic testing to trap possible
errors in the logic of our interpreter, but the user should never notice this.

@< Include... @>=
#include <memory>

@~Types are represented by small tree structures, with nodes of types
|type_declarator| that will be detailed later. That name indicates that they
describe the structure of values as accessible to the user, without attempting
to present to the user a closed abstraction (as \Cpp~classes do). Note however
that classes defined in the Atlas software will be presented to the user as
primitive types, so we do provide abstraction there.

Different types cannot share any subtrees among each other, so a root of a
type expression will own its tree. Usually types are built from the bottom up
(from leaves to the root), but during type checking a top-down approach where
types grow by substitution into type patters also occurs. In bottom-up
handling of types, they should not reside in local |type_declarator|
variables, since these would require (recursive) copying in order to become a
descendant type of another type. Instead they will be referred to by local
pointer values, and in order to ensure exception-safety, these should be
auto-pointers; for this reason we define |type_ptr| to be an auto-pointer
type. For the links between a node and its descendant types we use ordinary
types however. We also define a |type_list| type that will be used inside
tuple types, and a corresponding auto-pointer type |type_list_ptr|.

@< Type definitions @>=
struct type_declarator;
typedef std::auto_ptr<type_declarator> type_ptr;
typedef struct type_node* type_list;
typedef std::auto_ptr<type_node> type_list_ptr;

@ Since types and type lists own their trees, their copy constructors must
make a deep copy. Most applications of the copy constructor for types will be
via the function |copy|, which converts a pointer to a |type_declarator|, into
an auto-pointer to a copied instance of that type declarator. The change from
pointer to auto-pointer should remind us that we have gained ownership of the
copy.

@< Declarations of exported functions @>=
type_ptr copy(const type_declarator* t);

@~The deep copy is obtained by allocating a fresh |type_declarator|,
copy-constructing its contents from the structure pointed to by the argument
of~|copy|. The recursive copying will be implied by the definition if the copy
constructor that will be defined later.

@< Function definitions @>=
type_ptr copy(const type_declarator* t) @+
{@; return type_ptr(new type_declarator(*t)); }

@*2 Type lists.
%
The auxiliary type |type_list| gives a preview of matters related to
ownership, that will also apply, in a more complex setting, to
|type_declarator|. It is a simply linked list with component type~|t| given by
a sub-object; compared with using a pointer this it is more compact and
requires less (de)allocations. On the other hand it forced us to introduce the
|type_declarator::set_from| (see below) to make the basic constructor below
work efficiently; if we had initialised |t(*head)|, a recursive copy would
have been made. Also, we must make sure this definition
\emph{follows} that of |type_declarator| in the tangled code.

@< Definition of |struct type_node@;| @>=
struct type_node
{ type_declarator t; @+ type_list next;
@)
  type_node(type_ptr head, type_list tail) : t(),next(tail)
  @+{@; t.set_from(head); }
  type_node(const type_node& n);
  ~type_node();
};

@ Copy-constructing a type list is easy, since the node contains only one
pointer (another advantage of not using a pointer for the |t| field):
the recursive copy of the tail list is copy-constructed into the fresh
node.

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

@ Type lists are usually built by calling |make_type_list| or
|make_type_singleton|, rather than using a constructor. These functions are
efficient, due to the use of auto-pointers. A caller holding auto-pointers to
the descendants transfers ownership to the new node at the call, for which it
gets ownership upon returning. Using ordinary pointers one would run into
ownership conflicts that can only be solved by calling the (recursive) copy
constructor, which would also be necessary when using local |type_node|
variables.

We pass the auto-pointers by value, rather than by (necessarily non-constant)
reference. This implies a small overhead due to multiple ownership transfer,
but allows passing anonymous values obtained from a constructor or
construction function as argument (due to a clever definition of
auto-pointers); this is forbidden for non-constant reference parameters.

Note that we are prudent to release the auto-pointers only \emph{after} the
allocation for |new| is completed; it seems that \Cpp~standard does not
specify the order of events for |new type_node(t.release(),l.release())|
precisely enough to consider that safe. The fact that between the completion
of |new| and that of |return| nobody owns the new node (since we use an
ordinary pointer) should be harmless, as there is nothing that could throw an
exception there.


@< Function def... @>=

type_list_ptr make_type_list(type_ptr t,type_list_ptr l)
{@; type_node* p=new type_node(t,l.get()); l.release();
  return type_list_ptr(p);
}
@)
type_list_ptr make_type_singleton(type_ptr t)
{@; return type_list_ptr(new type_node(t,NULL)); }

@ Here is one more function that is convenient to have around.

@< Function definitions @>=
size_t length(type_list l)
@+{@; size_t len=0;
  for (; l!=NULL; l=l->next) ++len;
  return len;
}

@*2 Primitive types.
%
Before we can define |type_declarator|, we need to enumerate its variants and
the possibilities for primitive types. Some of them will be introduced later.

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
defined in the Atlas library, such as inner classes of real forms. Furthermore
one has types that are ``row of'' some other type, tuples (Cartesian products)
of some given sequence of types, and function types. We also allow for an
undetermined type, which can serve as a wild-card to specify type patterns.
Here are enumerations of tags for the basic kinds of types, and for the
sub-cases of |primitive_type|.

Type declarators are defined by a tagged union. Users will in general test the
tag and access the corresponding variant directly, so we give public access to
the data fields (it is not obvious how coherent variant selection could be
enforced by private data and public methods without greatly complicating
usage; the \Cpp~model of abstraction does not seem particularly suited for
tagged unions). By using an anonymous union, the field selectors like |prim|
of the variants in the union can be used directly on the level of the
|type_declarator|, thus avoiding an additional level of selection to access
them.

There is one restriction on types that is not visible in the definition below,
namely that the list of types referred to by the |tuple| field cannot have
length~$1$ (but length~$0$ is allowed). This is because anything that would
suggest a $1$-tuple (for instance a parenthesised expression) is identified
with its unique component.

We forbid the assignment operator, since the default definition would not do
the right thing, the proper definition would be somewhat complicated, and
without much use since types do not change dynamically. Two special cases are
provided for however, for replacing an undefined (or partially undefined) type
by a more specific type (this avoids the problem of cleaning of parts that
disappear). The |set_from| method makes a shallow copy of its argument into
the object for which is was called, which is required to be undetermined
initially; it is intended for use in constructors to avoid the deep copy of
the copy constructor. Since the shallow copy incorporates any descendants
without copying, the caller of |set_from| should have and relinquish ownership
of the argument type passed; therefore the argument is passed as an
auto-pointer by value. The |specialise| method is used during type analysis,
to see if our type matches a given pattern, or in case it was (partially)
undefined whether it can be made to match the pattern by if necessary
replacing some undetermined descendants by more specific ones. The call
returns a value indicating whether this was possible, and if so makes the
necessary specialisations to our type. That is done by copying, so the caller
does not require or lose ownership of the pattern for this method.

@< Type definitions @>=

struct func_type;
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
  explicit type_declarator(type_declarator* c)
    : kind(row_type) @+{@; component_type=c; }
  type_declarator(type_list component_types) : kind(tuple_type)
    @+{@; tuple=component_types; }
  type_declarator(type_ptr arg, type_ptr result); // for function types
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

@< Definition of |struct type_node@;| @>@;

@ This constructor cannot be defined inside the structure definition, since
|func_type| is not (and cannot be) a complete type there. The constructor for
|func_type| used will be defined below.

@< Function definitions @>=
inline type_declarator::type_declarator(type_ptr arg, type_ptr result)
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
purposes, but creates a slight problem for defining the basic constructor, if
we want to avoid using the copy constructor of |type_declarator| to
recursively copy argument and result types from the place they were previously
created into the two slots of one same structure. This problem is solved by
using the |get_from| method to fill the initially undetermined slots with a
shallow copy of the types, in the same way as was done in the constructor
for~|type_node|.

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

@ The method |set_from| is much like |specialise|, but it avoids making a deep
copy. Its argument is an auto-pointer passed by value, to remind users that
this argument must be |new|-allocated and that they are giving up ownership of
this argument. The contents of its top-level structure will be copied to the
current |type_declarator|, but without invoking the copy constructor;
afterwards the node will be destroyed after detaching any possible descendants
from it. This operation is only safe if the |type_declarator| previously had
no descendants, and in fact we insist that it had |kind==undetermined_type|;
if this condition fails we signal a |std::logic_error|. In a sense this is a
|swap| method, but only defined if the type was undetermined to begin with.

@h <stdexcept>
@< Function definitions @>=
void type_declarator::set_from(type_ptr p)
{ if (kind!=undetermined_type) throw std::logic_error("Illegal set_from");
  kind=p->kind;
  switch(kind)
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
we only continue is the top levels of both type declarers match, in which case
we try to recursively specialise all descendants. We do not guarantee
commit-or-roll-back, in other words, when the specialisation fails, some
modifications to our type may still have been made. This is not very serious,
since types are not shared and the initial call of specialisation always
starts with a fresh type, so at worst it might lead to a slightly misleading
type appearing in an error message (but even this is not easy to produce).

@< Function definitions @>=
bool type_declarator::specialise(const type_declarator& pattern)
{ if (kind==undetermined_type)
    {@;new(this) type_declarator(pattern); return true; }
  if (pattern.kind==undetermined_type) return true; // no need to refine
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
operator~`|<<|' by constant reference, which is seems more decent than doing
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

@< Function definitions @>=
type_ptr make_undetermined_type()
@+{@; return type_ptr(new type_declarator); }
@)
type_ptr make_prim_type(primitive_tag p)
@+{@; return type_ptr(new type_declarator(p)); }
@)
type_ptr make_row_type(type_ptr c)
{@; type_declarator* p=new type_declarator(c.get()); c.release();
    return type_ptr(p);
}
@)
type_ptr make_tuple_type (type_list_ptr l)
{@; type_declarator* p=new type_declarator(l.get()); l.release();
  return type_ptr(p);
}
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
    l0.get()!=NULL and l0->next==NULL ? copy(&l0->t) : make_tuple_type(l0);
  if (is_tuple) return a;
  type_ptr r =
    l1.get()!=NULL and l1->next==NULL ? copy(&l1->t) : make_tuple_type(l1);
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
no dynamic allocation is required for the top level structure. For
|row_of_type| and its relatives we use a copy of another type as descendant
(using the type itself would cause double deletion upon program termination),
and we should not forget to release the auto-pointer returned by |copy| lest
the copy be destroyed before we have a chance to use it. Recall that a
|type_declarator| constructed from a \emph{pointer to} the same constructs a
``row of'' type.

@< Global variable definitions @>=

const type_declarator unknown_type; // uses default constructor
const type_declarator int_type(integral_type);
const type_declarator str_type(string_type);
const type_declarator bool_type(boolean_type);
const type_declarator vec_type(vector_type);
const type_declarator mat_type(matrix_type);
const type_declarator row_of_type(copy(&unknown_type).release());
const type_declarator row_of_int_type(copy(&int_type).release());
const type_declarator row_of_vec_type(copy(&vec_type).release());
const type_declarator row_row_of_int_type(copy(&row_of_int_type).release());
const type_declarator int_int_type(*make_type("(int,int)"));

@*1 Dynamically typed values.
%
Now we shall consider runtime values. We could either use void pointers to
represent generic values and cast them when necessary, or use inheritance and
the dynamic cast feature of \Cpp. We choose the second option, which is quite
convenient to use, although this means that in reality we have dynamic type
information stored in the values as well as an external |type_declarator|
value describing their type. We shall use this information to double-check our
type analysis at runtime.

@< Includes needed in the header file @>=
#include <iostream>

@~We start with a base class for values. There must be at least one virtual
function in the class, which can conveniently be a function for printing. This
allows the base class to be defined abstract (one cannot declare a destructor
purely virtual since it will always be called, after the destructor for a
derived class). The printing function will demonstrate the ease of using
dynamic typing via inheritance. It does not even require any dynamic casting,
but other operations on values will. When identifiers were introduced, we also
needed to clone objects of types derived from |value_base|, so we added
another (purely) virtual method, |clone|. The method |name| is useful in
reporting logic errors from template functions; since these usually do not
have any object at hand to call this method from, it is defined |static|
rather than |virtual|. We forbid copying and assignment for the moment (the
class is abstract anyway).

@< Type definitions @>=
struct value_base
{ value_base() @+ {};
  virtual ~value_base() @+ {};
  virtual void print(std::ostream& out) const =0;
  virtual value_base* clone() const =0;
  static const char* name(); // just a model; this instance remains undefined
private: //copying and assignment forbidden
  value_base(const value_base& x);
  value_base& operator=(const value_base& x);
};
@)
typedef value_base* value;
typedef std::auto_ptr<value_base> owned_value;

@ We can already make sure that the operator~`|<<|' will do the right thing
for any of our values.

@< Declarations of exported functions @>=
std::ostream& operator<< (std::ostream& out, const value_base& v);

@~The operator~`|<<|' calls the (virtual) |print| method of the object pointed
to, and returns the reference to the output stream object. The copy
constructor and assignment constructor for |value_base| are currently
forbidden, but we provide loud versions to easily enable tracking them if this
should be needed.

@< Function definitions @>=
std::ostream& operator<< (std::ostream& out, const value_base& v)
{@; v.print(out); return out; }
@)
value_base::value_base(const value_base& x)
{@; std::cerr << "Copying " << x << std::endl; }
value_base& value_base::operator=(const value_base& x)
{@; std::cerr << "Assigning " << *this << "<-" << x << std::endl;
  return *this;
}

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

@ Here we define one more type derived from |value_base|, namely the type for
``row of'' types.

@< Includes needed in the header file @>=
#include <vector>

@~Using vectors from the standard template library, the realisation of row
values is quite easy. Since row values hold pointers that they own, it is
useful to define an auto-pointer type for them as well.

@< Type definitions @>=
struct row_value : public value_base
{ std::vector<value> val;
@)
  explicit row_value(const std::vector<value>& v) : val(v) @+{}
  ~row_value();
  void print(std::ostream& out) const;
  size_t length() const @+{@; return val.size(); }
  row_value* clone() const @+{@; return new row_value(*this); }
  static const char* name() @+{@; return "row value"; }
protected:
  row_value(const row_value& v);
private:
  row_value& operator=(const row_value& v);
};
@)
typedef std::auto_ptr<row_value> row_ptr;

@ Before we forget it, let us define the copy constructor (needed for the
|clone| method) and the destructor for |row_value| objects. Since the |val|
field contains a vector of pointers, we must explicitly clone respectively
delete the pointed-to objects. We use here the fact (without which life would
be much harder) that |delete| will do the right thing even if it is called
with a pointer to a base class of the class for which the pointer was created.

@< Function definitions @>=
row_value::row_value(const row_value& v) : val(v.val)
{ for (std::vector<value>::iterator p=val.begin() ; p!=val.end(); ++p)
    *p=(*p)->clone();
}

row_value::~row_value()
{@; for (std::vector<value>::iterator p=val.begin(); p!=val.end(); ++p)
    delete *p;
}

@ So here is the first occasion where we shall use virtual functions. For the
moment the output routine performs an immediate recursion; later we shall try
to make this more elegant by computing the width needed to output component
values, and adapt the formatting to that.

@< Function definitions @>=
void row_value::print(std::ostream& out) const
{ if (val.empty()) out << "[]";
  else
  { out << '[';
    std::vector<value>::const_iterator p=val.begin();
    do {@; (*p)->print(out); ++p; out << (p==val.end() ? ']' : ','); }
    while (p!=val.end());
  }
}

@ Often we know what variant a |value| object takes, based on the type
analysis. We can convert to that type using a |dynamic_cast|, but at such
moments we wish to throw a |logic_error| in case our type prediction was
wrong. To avoid having such casts and |throw| statements all over the place,
we define a template function to do the casting and throwing.

@< Template and inline function definitions @>=
template <typename D> // |D| is a type derived from |value_base|
 D* force(value v) throw(std::logic_error)
{ D* p=dynamic_cast<D*>(v);
  if (p==NULL) throw
    std::logic_error(std::string("forced value is no ")+D::name());
  return p;
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

@ Like for values, we can assure right away that printing converted
expressions will work.

@< Declarations of exported functions @>=
inline std::ostream& operator<< (std::ostream& out, const expression_base& e)
{@; e.print(out); return out; }


@ Since our |evaluate| methods will put values on the execution stack, let us
declare it right away.

@< Declarations of global variables @>=
extern std::vector<value> execution_stack;

@~We define the stack as a static variable of this compilation unit; it is
initially empty. All usable built-in functions will be provided with a small
wrapper function that takes it values from the stack and places its results
there again. Parameters are placed on the stack in order, and should therefore
be popped from the stack in reverse order.

@< Global variable definitions @>=
std::vector<value> execution_stack;

@ Here are two functions to facilitate manipulating the stack. Being inline,
we define them right away. There is a remote chance that the |push_back|
method runs out of stack space and throws an exception; to prevent a memory
leak we therefore temporarily place the pointer in an |owned_value|. If the
value already resides in an auto-pointer (possibly of a derived type), we also
allow passing that pointer by value (which might be more efficient than
releasing it to be captured just afterwards), by defining a second variant of
|push_value|; this has to be a template function since a conversion from
auto-pointer-of-derived-type to |value_base| is not possible without a cast in
a function argument position. For |pop_value| we don't return an auto-pointer,
since in most cases we have to perform a |dynamic_cast| to the result anyway.

@< Template and inline function definitions @>=
inline void push_value(value v)
{@; owned_value safe(v); execution_stack.push_back(v); safe.release(); }
template<typename D> // |D| is a type derived from |value_base|
  inline void push_value(std::auto_ptr<D> v)
  @+{@; execution_stack.push_back(v.get()); v.release(); }
@)
inline value pop_value()
{@; value arg=execution_stack.back(); execution_stack.pop_back();
    return arg;
}

@ Now let us define classes derived from |expression_base|. The simplest is
|denotation|, which simply stores a |value|, which it returns upon
evaluation. A denotation owns the value it stores, which it must therefore
clone upon evaluation and delete upon destruction. This is not optimally
efficient for denotations that are evaluated just once, but in those
situations (evaluating expressions entered by the user directly) efficiency is
not a major concern.


@< Type definitions @>=
struct denotation : public expression_base
{ value denoted_value;
@)
  explicit denotation(value v) : denoted_value(v) @+{}
  virtual ~denotation() @+ {@; delete denoted_value; }
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
@+{@; push_value(denoted_value->clone()); }


@ Let us finish with some functions related to the execution stack. The stack
owns the values it contains, but there was no reason to wrap it into a class
with a destructor, since we never intend to destroy the stack: if our program
exits either peacefully or by an uncaught exception we don't care about some
values that are not destroyed. We must remember however to |delete| the values
whenever we empty the stack after catching a runtime error, so we provide a
function to clear the stack.

@< Declarations of exported functions @>=
void clear_execution_stack ();

@~We shall monitor disappearing values when clearing the stack. This provides
some output that may not be easy to interpret for an unsuspecting user, so we
limit this to the case that the |verbosity| parameter (defined below) is
nonzero.

@< Function definitions @>=
void clear_execution_stack ()
{ if (!execution_stack.empty())
  { if (verbosity!=0)
      std::cerr << "Discarding from execution stack:" << std::endl;
    do
    { value v=execution_stack.back();
      if (verbosity!=0) std::cerr << *v << std::endl;
      delete v; execution_stack.pop_back();
    }
    while (!execution_stack.empty());
  }
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
thus avoiding the use of one of the two auto-pointers, at the price of some
loss of symmetry and readability.

@< Function definitions @>=
type_error::type_error(const type_error& e)
 : program_error(e), offender(e.offender)
{@; type_ptr p=copy(e.required);
  actual=copy(e.actual).release(); required=p.release();
}


@* The evaluator.
%
With the translation of expressions into |expression| values, there is no need
for an evaluator in the strict sense, since we shall just call the |evaluate|
method. Nevertheless we leave in places the scaffolding of the direct
evaluator that was previously used, until the new system is fully functional
(it actually is, but it is currently only used when the \.{verbose} command
has been given to the interpreter). This evaluator is realised as a function
that maps expressions produced by the parser directly to values of a type
derived from |value_type|

@< Declarations of exported functions @>=
value evaluate(expr e)
  throw(std::bad_alloc,std::logic_error,std::runtime_error);

@~ Here we shall define just the outer level of |evaluate| and the simplest
cases of its main |switch|, namely those for denotations, leaving the cases
for other syntactic categories until later. The functions involved with type
checking will be similarly presented, so that most code will be presented
grouped together based on the syntactic category involved, rather than on the
particular function the code belongs to.

Evaluating a denotation is trivial, and yields the constant value stored in
the denotation. We assign the value to be returned to a variable rather than
calling |return| in the various branches, in order to make the conversion to
|value| less implicit, and to avoid compiler warnings about ending a non-void
function without~|return|. (It is fairly curious that such a warning should be
given when there is no warning for missing cases in the switch, and all cases
end with |return|; other functions defined below show alternative methods to
keep the compiler from complaining.)

@h<cstdlib>
@< Function definitions @>=
value evaluate(expr e)
  throw(std::bad_alloc,std::logic_error,std::runtime_error)
{ value result=NULL;
  switch(e.kind)
  { case integer_denotation:
      result= new int_value(e.e.int_denotation_variant);
      break;
    case string_denotation:
      result= new string_value(e.e.str_denotation_variant);
      break;
    case boolean_denotation:
      result= new bool_value(e.e.int_denotation_variant);
      break;
    @\@< Cases for evaluating other kinds of expressions, which store the
    value obtained in |result|  @>
  }
  return result;
}

@*1 Finding the type of an expression.
%
Although the above code can directly evaluate expressions as returned by the
parser, we wish to proceed with more precaution, and check types first. One
reason is that it is best to signal type errors early and at a level where
they can be understood by the user, rather than discovering a wrong type in
the middle of the computation where the context of the actual  error has
disappeared. Type checking also allows us to embark upon the construction of
rigidly typed values needed  by the Atlas software, such as
|std::vector<int>|, with confidence that no values of the wrong type will pop
up in the process.

Type checking can appear in several forms. The simplest idea is to take an
expression and find its type, which is what the function |find_type| does. In
this approach it may happen that |find_type| is called for a subexpression for
which it is known beforehand that only one type is acceptable; this happens
notably for arguments of (non-overloaded) functions. In this case it is more
useful to proceed with knowledge of that type; among other things, this allows
conversion functions to be inserted in a more flexible way. The function
|check_type| that will be given later does this. Finally, There is an approach
that will eventually replace both |find_type| and |check_type|, which behaves
like |check_type| but with a given type pattern that is partially
undetermined, and which is specialised as a result of the checking process.
This functionality is implemented in |convert_expr|, which in the mean time
also converts the expression from the form returned by the parser into a
\Cpp-object of a class derived from~|expression_base|, whose |evaluate| method
can be used instead of calling the function |evaluate|. As for |evaluate|, we
begin with an outline of this functions that only handle the denotation case.

@< Declarations of exported functions @>=
type_ptr find_type (expr e) throw(std::bad_alloc,program_error);
void check_type (const type_declarator& t,expr& e)
   throw(std::bad_alloc,program_error);
expression convert_expr(const expr& e, type_declarator& type)
  throw(std::bad_alloc,program_error);

@~For all of these functions we shall proceed by a straightforward traversal
of the parse tree; |find_type| is the simplest one, with type information
flowing from the bottom (leaves) up uniquely.

@< Function definitions @>=
type_ptr find_type (expr e) throw(std::bad_alloc,program_error)
{ switch(e.kind)
  { case integer_denotation: return make_prim_type(integral_type);
    case string_denotation: return make_prim_type(string_type);
    case boolean_denotation: return make_prim_type(boolean_type);
  @\@< Cases for finding the type of other kinds of expressions, all of which
       either |return| the type found, or |throw| a |type_error| @>
  }
  return type_ptr(NULL); // keep the compiler happy, never reached
}

@*1 Checking if an expression has a given type.
A refinement of the notion of type checking will occur in contexts that we
have not encountered yet, such as function arguments, where in descending the
parse tree we already know which type is required, and we must check if it
matches the actual type of the expression. The distinction is important
because of implicit type conversions: if the type is already known, such
conversions can be inserted (and this is in fact the only way that we shall be
able to construct such built-in types as matrices), and by taking into account
the required type during the descent, the conversions can be more flexible.
For instance, we have insisted in |find_type| that types in a list display be
exactly equal, but if the (row) type of the lists display is known beforehand,
it will be sufficient that each component expression can individually be
converted into the required component type.

@~In spite of the fact that a type is specified, our first concern is analyse
the expression, as for |find_type|. We adopt a somewhat peculiar convention
that an error is signalled by setting the type |actual| to a non-null value;
this avoids having to place |throw| expressions in many places. When we do
throw a |type_error| after finding |actual| to be set, we need to copy the
required type~|t| for inclusion in the error object.

@< Function definitions @>=
void check_type (const type_declarator& t,expr& e)
   throw(std::bad_alloc,program_error)
{ type_ptr actual(NULL);
  switch(e.kind)
  { case integer_denotation: if (t!=int_type) actual= copy(&int_type);
    break;
    case string_denotation:  if (t!=str_type) actual= copy(&str_type);
    break;
    case boolean_denotation: if (t!=bool_type) actual= copy(&bool_type);
    break;
  @\@< Other cases for testing whether the type of |e| matches |t|,
       which in case of failure set |actual| to the erroneous type @>
  }
  if (actual.get()!=NULL) throw type_error(e,actual,copy(&t));
@.Type error@>
}

@*1 Type checking with conversion to executable values.
%
The function |convert_expr| type-checks an expression, given as |expr| value,
and upon success converts it to an |expression| value. Apart from the
expression, it takes a type as argument in the form of a non-constant
reference to a |type_indication|; this may be a partially defined type (such
as \.{(int,*)} for ``pair on an integer and something'') initially, and upon
success it will have become a completely defined type for the expression
(unless the expression can actually have multiple types, such the empty list
which gets type~\.{[*]}). The type should be owned by the caller, and we never
destroy any of it in this function; all changes are affected by calling the
|specialise| method for |type| or for its descendants, so in particular any
new parts to the type are freshly copied.

@ The function |convert_expr| descends the expression tree in the same way as
|evaluate| does. For the cases of denotations we first try to specialise to
the required type and test whether that succeeded; thus the cases where |type|
is initially undetermined or equal to the correct type are handled together.

@< Function definitions @>=
expression convert_expr(const expr& e, type_declarator& type)
  throw(std::bad_alloc,program_error)
{ switch(e.kind)
  { case integer_denotation:
      if(type.specialise(int_type)) // change undefined type to integral
        return new denotation(new int_value(e.e.int_denotation_variant));
      else throw type_error(e,copy(&int_type),copy(&type));
   case string_denotation:
      if (type.specialise(str_type)) // change undefined type to string
        return new denotation(new string_value(e.e.str_denotation_variant));
      else throw type_error(e,copy(&str_type),copy(&type));
   case boolean_denotation:
      if (type.specialise(bool_type)) // change undefined type to boolean
        return new denotation(new bool_value(e.e.int_denotation_variant));
      else throw type_error(e,copy(&bool_type),copy(&type));
   @\@< Other cases for type-checking and converting expressions, all of which
      either |return| or |throw| a |type_error| @>
 }
 return NULL; // keep compiler happy
}

@*1 List displays.
%
Now we shall consider how these various functions handle list displays,
expressions that build a list of values by constructing each component by an
explicit expression. We start by defining the structure used to represent list
displays after conversion in the type check. The logical name |list_display|
for this structure is already taken, so we call it a |list_expression|
instead.

@< Type definitions @>=
struct list_expression : public expression_base
{ std::vector<expression> component;
@)
  explicit list_expression(std::vector<expression> l) : component(l) @+{}
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


@ If a list display has are multiple components, they must all have the same
type, and the result type will be ``row of'' that type. There are two possible
complications: (1)~when there are no components at all, the component type
could be anything, and (2)~if some non-row type is required (such as vector or
matrix) a type error may be avoided by inserting a conversion to that type,
provided it exists and the components have an appropriate type (the simplest
case is converting ``row of integer'' to ``vector''). Complication (1) only
occurs if no final type is required (as in |find_type|), while (2) only occurs
in the opposite case (as in |check_type|).

In |find_type| we ``solve'' the problem of finding the type of an empty row
display by arbitrarily giving it type ``row of integer''. In other cases the
first type found will be required in all other positions; this is obtained by
calling |check_type| instead of recursively calling |find_type|.

@< Cases for finding the type of other kinds of expression... @>=
case list_display:
{ if (e.e.sublist==NULL) return copy(&row_of_int_type);
  type_ptr component_type=find_type(e.e.sublist->e);
  for (expr_list l=e.e.sublist->next; l!=NULL; l=l->next)
    check_type(*component_type,l->e);
  return make_row_type(component_type);
}

@ As we mentioned, |check_type| has no difficulty with empty list displays,
but it has to do something if the required type is not a row type. In most
cases it will then have to throw a |type_error|, but there are some cases that
can be salvaged by inserting a conversion function; the details for this will
be given later. Note that the possibility of automatically inserting these
conversions was one of the most compelling reasons that we chose to have type
checking even in the very earliest versions of this interpreter.

If the type required for a list display is a row type, we can simply test that
all component expressions have the required component type. In the case where
no known type conversion can solve the conflict between the required type and
a list display, we want to signal an error without further checking the list
display (since this could come across other less relevant type errors, and an
error so thrown would preempt our attempt to signal the error), so we signal
the type found by the pattern \.{[*]}, an unspecified row type.

@< Other cases for testing whether the type of |e| matches |t|... @>=
case list_display:
  if (t.kind==row_type)
    for (expr_list l=e.e.sublist; l!=NULL; l=l->next)
      check_type(*t.component_type,l->e);
  else
  { @< If the type |t| can be converted from a list type, insert the
       conversion into |e|, and check the types of the component expressions
       in |e.e.sublist|, then |break| @>
    actual=copy(&row_of_type);
    // otherwise signal a ``row of'' type where none can be handled
  }
break;


@ The function |convert_expr| handles the cases dealt with above in an
original manner. If the required type was completely undetermined, we
specialise it to a row type with undetermined component type before proceeding
recursively with the component expressions; thus we build up our type from the
top (root) down. As a side effect of this strategy, the type will remain
partly undefined in the case of an empty list, and it is even possible that
for a list display whose first component has such a partly defined type, the
further components will progressively refine this. This approach also allows
us to conveniently merge the cases corresponding to |find_type| and
|check_type| above. Of course we do have the additional obligation here to
build a converted |list_expression| representing the list display; we know its
size beforehand, and upon allocation we fill it with null pointers that will
progressively be replaced, in order to ensure proper destruction in case a
|type_error| is thrown by a component expression. When a type other than ``row
of'' is expected, we must of course take explicit action to see whether some
type conversion can resolve the conflict; as for |check_type| we postpone that
code for now.

@< Other cases for type-checking and converting expressions... @>=
case list_display:
  if (type.specialise(row_of_type))
  { std::auto_ptr<list_expression> result@|
      (new list_expression@|
       (std::vector<expression>(length(e.e.sublist),NULL)));
  @/size_t i=0;
    for (expr_list l=e.e.sublist; l!=NULL; l=l->next)
      result->component[i++]=convert_expr(l->e,*type.component_type);
    return result.release();
  }
  else
  @< If |type| can be converted from some row-of type, check the components of
     |e.e.sublist| against the required type, and apply a conversion function
     to the converted expression; otherwise |throw| a |type_error| @>

@ Evaluating a list display is quite straightforward. We first calculate the
length the resulting value will have, then construct a |row_value| object of
the proper size, pointed to by an auto-pointer and filled with null pointers.
These null pointers are then successively replaced by the results of
evaluating the expressions in the list display. Proceeding in this order
guarantees that if any evaluation should throw an exception, then the results
of previous evaluations will be destroyed.

@< Cases for evaluating other kinds of expressions... @>=
case list_display:
{ row_ptr v(new row_value(std::vector<value>(length(e.e.sublist),NULL)));
  size_t i=0;
  for (expr_list l=e.e.sublist;l!=NULL; l=l->next)
    v->val[i++]=evaluate(l->e);
  result=v.release();
}
break;

@ The evaluation of a |list_expression| evaluates the components in a simple
loop and wraps the result in a |row_value|. Since evaluation pushes the
resulting value onto the execution stack, we pop it off into the result
afterwards. We take care to hold the partial result via an auto-pointer
|result|, so that in case of a runtime error during the evaluation of one of
the component expressions the values already computed are cleaned up.

@< Function def... @>=
void list_expression::evaluate() const
{ row_ptr result(new row_value(std::vector<value> (component.size(),NULL)));
  for (size_t i=0; i<component.size(); ++i)
    component[i]->evaluate(),result->val[i]=pop_value();
  push_value(result);
}

@*1 Conversion to rigid vectors and matrices.
%
The following sections deal with interfacing to the Atlas library rather than
with the evaluator proper: they deal with converting data from the evaluator's
``native'' format of nested lists to vectors and matrices that can be used by
the library. Currently the conversions from the type ``row of integer'' to the
type ``vector'' and from the types ``row of row of integer'' and ``row of
vector'' to the type ``matrix'' are implemented. Here ``vector'' really means
|latticetypes::Weight|, which stands for |std::vector<int>|, while ``matrix''
means |latticetypes::LatticeMatrix|, which stands for |matrix::Matrix<int>|.
These declarations are given in~\.{structure/latticetypes\_fwd.h}, while the
template class |matrix::Matrix| is defined in~\.{utilities/matrix.h}.


@< Includes... @>=
#include "latticetypes_fwd.h"
#include "matrix.h" // this makes |latticetypes::LatticeMatrix| a complete type

@ These conversion functions to vector and matrix will be defined below.

@< Declarations of local functions @>=
latticetypes::Weight cast_intlist_to_weight(const value);
latticetypes::LatticeMatrix cast_intlistlist_to_matrix(const value);

@*2 Inserting conversion functions.
%
In |check_type| we come to the code below when a list display is found in a
position where some non-row type is required. If that required type is vector
or matrix, we proceed as if it were row of integer respectively row of vector,
while inserting a call to the appropriate conversion function into the
expression tree. Since we are already past the point where the case of a
required ``row of'' type is treated, we have to repeat its code, which is that
of recursively calling |check_type| in a loop over the component expressions.
We can however merge the cases of required types vector and matrix, using a
variable~|comp| for the required component type.

@< If the type |t| can be converted from a list type... @>=
{ expr_list l=e.e.sublist; // the list of component expressions
  type_ptr comp(NULL); // the type that they should have
  if (t==vec_type)
  {@; @< Insert a vector conversion at |e| @>
    comp=copy(&int_type);
  }
  else if (t==mat_type)
  {@; @< Insert a matrix conversion at |e| @>
    comp=copy(&vec_type);
  }
  if (comp.get()!=NULL) // a substitute component type was set
  { for (; l!=NULL; l=l->next) check_type(*comp,l->e);
    break; // from enclosing |switch|
  }
}

@ A function call is created by calling |make_application_node| with as
arguments the code for a function identifier (currently the only possibility)
and a pointer to the argument list. Here there is a single argument, namely
the list display, so we build a singleton list from it. The function
identifier is obtained from the main hash table, and since the conversion
functions are not intended to be entered explicitly by the user, they are
entered into the table under names `\.{>vec<[int]:}' and `\.{>mat<[vec]:}'
that cannot result from user input. This module and the following one will be
reused in other cases where a conversion function needs to be inserted.

@< Insert a vector conversion at |e| @>=
{ expr_list arg=make_exprlist_node(e,null_expr_list);
  e=make_application_node
       (main_hash_table->match_literal(">vec<[int]:"),arg);
}

@~The only difference between vector and matrix conversion is the name of the
conversion function.

@< Insert a matrix conversion at |e| @>=
{ expr_list arg=make_exprlist_node(e,null_expr_list);
  e=make_application_node
       (main_hash_table->match_literal(">mat<[vec]:"),arg);
}

@ In the setting of |convert_expr|, the conversion to vector or matrix will be
represented by a special kind of converted expression, which behaves like a
function call of a fixed function. For each function a different structure
derived from |expression_base| is used; they only differ among each other by
their virtual methods |evaluate| and |print|.

@< Type definitions @>=
struct vector_conversion : public expression_base
{ expression exp;
@)
  explicit vector_conversion(expression e) : exp(e) @+{}
  virtual ~vector_conversion()@;{@; delete exp; }
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};
@)
struct matrix_conversion : public vector_conversion
{ explicit matrix_conversion(expression e) : vector_conversion(e) @+{}
@)
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
    else throw type_error(e,copy(&row_of_type),copy(&type));
    std::auto_ptr<list_expression> display@|
      (new list_expression@|
        (std::vector<expression>(length(e.e.sublist),NULL)));
  @/size_t i=0;
    for (expr_list l=e.e.sublist; l!=NULL; l=l->next)
      display->component[i++]=convert_expr(l->e,component_type);
    expression result=
      is_vector ? new vector_conversion(display.get())
		: new matrix_conversion(display.get());
    display.release();
    return result;
  }

@*2 Coercion of types.
%
The cases above, where a list display in a context requiring a vector or a
matrix is transformed by following it by a conversion function, and the most
important application of these conversions, since they provide a way for the
user to enter a value of type vector or matrix that would not otherwise be
available. However, it seems logical to provide the same conversions also when
the expression is not a list display, for instance an applied identifier or a
function call. In those cases the need for a conversion will usually arise
from a difference between the required type and an already fully determined
type of the subexpression (this is clear for an applied identifier, and for a
(non-overloaded) function call the result type is also available even before
analysing the argument expression). Since multiple syntactic classes are
involved, it is not attractive to handle such conversions separately for the
various branches of the big switch on the expression~|kind| in the type
analysis, so we put in place a general mechanism that will insert a conversion
depending on the pair of the expression type and required type, if it exists.
As this mechanism is fairly new, it is not incorporated in the obsolescent
functions |find_type| and |get_type|, but only in |convert_expr|.

The mechanism is realised by a function |coerce|. Its final argument~|e| is a
reference to the previously converted expression having |from_type|, and if
values of this type can be converted to |to_type| then the expression will be
modified by insertion of a conversion expression; the return value indicates
whether a conversion was found.

@< Declarations of exported functions @>=

bool coerce( const type_declarator& from_type, const type_declarator& to_type
	   , expression& e) throw(std::bad_alloc);

@ The implementation of |coerce| will be determined by a simple table lookup.

@< Local type definitions @>=
struct conversion_info
{ typedef expression (*converter)(expression);
  type_declarator from,to; converter conv;
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

@ The function |coerce| simply traverses the |coerce_table| looking for an
appropriate entry, and applies the |converter| if it finds one.

@< Function definitions @>=
bool coerce( const type_declarator& from_type, const type_declarator& to_type
	   , expression& e) throw(std::bad_alloc)
{ for (size_t i=0; i<nr_coercions; ++i)
    if (from_type==coerce_table[i].from and to_type==coerce_table[i].to)
    @/{@; e=coerce_table[i].conv(e); return true; }
  return false;
}

@ We need functions to instantiate |conversion_info::conv| from. We need to
declare these functions before we come to the global variable definitions.

@< Declarations of local functions @>=
expression vector_converter(expression e);
expression matrix_converter(expression e);

@ These conversion-inserting functions are in general quite simple packaging
routines; here are the ones for |vector_conversion| and |matrix_conversion|.

@< Function definitions @>=
expression vector_converter(expression e) @+
{@; return new vector_conversion(e); }
expression matrix_converter(expression e) @+
{@; return new matrix_conversion(e); }


@ All that remains is to initialise the |coercion table|.
@<Initialiser for |coerce_table| @>=
conversion_info(copy(&row_of_int_type), copy(&vec_type), vector_converter), @/
conversion_info(copy(&row_of_vec_type), copy(&mat_type), matrix_converter), @/
@[@]

@ We shall now extend our range of conversions with other possibilities. For
each we need to define a type derived from |expression_base|; in fact we can
derive from |vector_conversion|, like |matrix_conversion| did.

@<Type definitions@>=

struct matrix2_conversion : public vector_conversion
  // for \.{[[int]]}$\to$\.{mat}
{ explicit matrix2_conversion(expression e) : vector_conversion(e) @+{}
@)
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};

struct int_list_conversion : public vector_conversion
  // for \.{vec}$\to$\.{[int]}
{ explicit int_list_conversion(expression e) : vector_conversion(e) @+{}
@)
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};
@)
struct vec_list_conversion : public vector_conversion
  // for \.{mat}$\to$\.{[vec]}
{ explicit vec_list_conversion(expression e) : vector_conversion(e) @+{}
@)
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};
@)
struct int_list_list_conversion : public vector_conversion
  // for \.{mat}$\to$\.{[[int]]}
{ explicit int_list_list_conversion(expression e) : vector_conversion(e) @+{}
@)
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};

@ The print functions are simple. They are not inside the structure
definition, since their second `|<<|' cannot be matched there.

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
expression matrix2_converter(expression e);
expression int_list_converter(expression e);
expression vec_list_converter(expression e);
expression int_list_list_converter(expression e);

@~Their definitions are hardly surprising

@< Function definitions @>=
expression matrix2_converter(expression e)
{@; return new matrix2_conversion(e); }
expression int_list_converter(expression e)
{@; return new int_list_conversion(e); }
expression vec_list_converter(expression e)
{@; return new vec_list_conversion(e); }
expression int_list_list_converter(expression e)
{@; return new int_list_list_conversion(e); }


@ All that remains is to initialise the |coercion table|.
@<Initialiser for |coerce_table| @>=
conversion_info(copy(&row_row_of_int_type), copy(&mat_type)
  , matrix2_converter), @/
conversion_info(copy(&vec_type), copy(&row_of_int_type)
  , int_list_converter), @/
conversion_info(copy(&mat_type), copy(&row_of_vec_type)
  , vec_list_converter), @/
conversion_info(copy(&mat_type), copy(&row_row_of_int_type)
  , int_list_list_converter), @/
@[@]


@*2 Primitive types for vectors and matrices.
%
When packed into a a type derived from |value_base| (defined later), these
types will be considered as primitive at the level of the evaluator.

@< Other primitive type tags @>=
vector_type, matrix_type, @[@]

@~In choosing the short names below, we choose to hide the specific
mathematical meaning that was implied in the names these types have in the
Atlas software. We believe that a more extensive name might be more confusing
that helpful to users; besides, the interpretation of the values is not
entirely fixed (vectors are used for coweights and (co)roots as well as for
weights, and matrices could denote either a basis or an automorphism of a
lattice.

@< Other primitive type names @>=
"vec", "mat", @[@]

@ Here are the corresponding types derived from |value_base|.

@< Type definitions @>=

struct vector_value : public value_base
{ latticetypes::Weight val;
@)
  explicit vector_value(const latticetypes::Weight& v) : val(v) @+ {}
  ~vector_value()@+ {}
  virtual void print(std::ostream& out) const;
  vector_value* clone() const @+{@; return new vector_value(*this); }
  static const char* name() @+{@; return "vector"; }
private:
  vector_value(const vector_value& v) : val(v.val) @+{}
};
@)
typedef std::auto_ptr<vector_value> vector_ptr;
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
@)

@ To make a small but visible difference, weights will be printed in equal
width fields one longer than the minimum necessary.
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
We now define the functions that perform the required conversions at run
time: these are the functions |cast_intlist_to_weight| and
|cast_intlistlist_to_matrix|, both of which use an auxiliary function
|row_to_weight|. The former two use up and destroy their arguments, while the
latter does not. Since exceptions may be thrown in many places, the easiest
way to ensure proper destruction is placement of the passed pointer into an
auto-pointer of type |row_ptr|. In retrospect it might have been better to
define the former two directly as wrapper functions, which currently are left
to be defined later.

Some people might complain that |row_to_weight| makes a copy of the vector
|result| (and its entries) at the |return| statement just before destroying
the local variable. However since there is only one |return| statement, the
code generated by a decent compiler like \.{g++} creates the vector object at
its destination in this case, and does not call the copy constructor (more
generally this is avoided if all |return| statements use the same variable).
So we shall not hesitate to use this idiom henceforth.

@< Function definition... @>=
latticetypes::Weight row_to_weight(const row_value& r)
{ latticetypes::Weight result(r.val.size());
  for(size_t i=0; i<r.val.size(); ++i)
    result[i]=force<int_value>(r.val[i])->val;
  return result;
}
@)
latticetypes::Weight cast_intlist_to_weight(const value v)
{ row_ptr r(force<row_value>(v));
  return row_to_weight(*r); // destroys |*r| on returning
}

@ Converting to a lattice matrix is slightly longer but not really more
complicated. When the call was inserted by the code above, the components of
the list display will have type ``vector'' (possibly after conversion from
``row of integer''), but this function accepts a ``row of row of integer''
value as well (which could be the value of an applied identifier). The
component vectors or lists will become the columns of the matrix, although
this fact is not evident from the code below; it is due to the way in which
the constructor of a |latticetypes::LatticeMatrix| from a
|latticetypes::WeightList| is defined.

@< Function definition... @>=
latticetypes::LatticeMatrix cast_intlistlist_to_latmat(const value v)
{ row_ptr rr(force<row_value>(v)); // we must destroy |*rr| at end
  latticetypes::WeightList res_vec(rr->val.size(),latticetypes::Weight());
  size_t depth=0; // maximal length of vectors
  for(size_t i=0; i<rr->val.size(); ++i)
  { row_value* r=dynamic_cast<row_value*>(rr->val[i]);
    if (r!=NULL)
      {@; latticetypes::Weight v=row_to_weight(*r); res_vec[i].swap(v);}
    else
      res_vec[i].swap(force<vector_value>(rr->val[i])->val);
    if (res_vec[i].size()>depth) depth=res_vec[i].size();
  }
  for(size_t i=0; i<res_vec.size(); ++i)
    if (res_vec[i].size()<depth)
    { size_t j=res_vec[i].size();
      res_vec[i].resize(depth);
      for (;j<depth; ++j) res_vec[i][j]=0; // extend weights if necessary
    }
  return latticetypes::LatticeMatrix(res_vec);
}

@ For |vector_conversion| and |matrix_conversion|, the |evaluate| methods are
relatively easy; we know that when a |matrix_conversion| was generated, the
expression |exp| returns a value of type \.{[vec]} (a row of vectors).

@< Function definitions @>=
void vector_conversion::evaluate() const
{ exp->evaluate(); @+ row_ptr r(get<row_value>());
  push_value(new vector_value(row_to_weight(*r)));
}
@)
void matrix_conversion::evaluate() const
{ exp->evaluate(); @+ row_ptr r(get<row_value>());
  latticetypes::WeightList res_vec(r->val.size(),latticetypes::Weight());
  size_t depth=0; // maximal length of vectors
  for(size_t i=0; i<r->val.size(); ++i)
  { res_vec[i].swap(force<vector_value>(r->val[i])->val);
    if (res_vec[i].size()>depth) depth=res_vec[i].size();
  }
  for(size_t i=0; i<res_vec.size(); ++i)
    if (res_vec[i].size()<depth)
    { size_t j=res_vec[i].size();
      res_vec[i].resize(depth);
      for (;j<depth; ++j) res_vec[i][j]=0; // extend weights if necessary
    }
  push_value(new matrix_value(latticetypes::LatticeMatrix(res_vec)));
}

@ Here is the remaining ``internalising'' conversion function, from row of row
of integer to matrix.

@< Function definitions @>=
void matrix2_conversion::evaluate() const
{ exp->evaluate(); @+ row_ptr r(get<row_value>());
  latticetypes::WeightList res_vec(r->val.size(),latticetypes::Weight());
  size_t depth=0; // maximal length of vectors
  for(size_t i=0; i<r->val.size(); ++i)
  { latticetypes::Weight v=row_to_weight(*force<row_value>(r->val[i]));
    res_vec[i].swap(v);
    if (res_vec[i].size()>depth) depth=res_vec[i].size();
  }
  for(size_t i=0; i<res_vec.size(); ++i)
    if (res_vec[i].size()<depth)
    { size_t j=res_vec[i].size();
      res_vec[i].resize(depth);
      for (;j<depth; ++j) res_vec[i][j]=0; // extend weights if necessary
    }
  push_value(new matrix_value(latticetypes::LatticeMatrix(res_vec)));

}

@ The other three conversion functions perform the inverse transformations of
those given above. Again it will be handy to have a basic function
|weight_to_row| that performs more or less the inverse transformation of
|row_to_weight|, but rather than returning a |row_value| it returns a |value|
pointing to it.

@< Function definitions @>=
value weight_to_row(const latticetypes::Weight v)
{ row_ptr result
    (new row_value(std::vector<value>(v.size(),NULL)));
  for(size_t i=0; i<v.size(); ++i)
    result->val[i]=new int_value(v[i]);
  return result.release();
}
@)
void int_list_conversion::evaluate() const
{ exp->evaluate(); @+ vector_ptr v(get<vector_value>());
  push_value(weight_to_row(v->val));
}
@)
void vec_list_conversion::evaluate() const
{ exp->evaluate(); @+ matrix_ptr m(get<matrix_value>());
  row_ptr result
    (new row_value(std::vector<value>(m->val.numColumns(),NULL)));
  for(size_t i=0; i<m->val.numColumns(); ++i)
  { vector_value*column=new vector_value(latticetypes::Weight());
    result->val[i]=column; // now |column | is owned
    m->val.column(column->val,i); // insert column into result vector
  }
  push_value(result.release());
}
@)
void int_list_list_conversion::evaluate() const
{ exp->evaluate(); @+ matrix_ptr m(get<matrix_value>());
  row_ptr result
    (new row_value(std::vector<value>(m->val.numColumns(),NULL)));
  for(size_t i=0; i<m->val.numColumns(); ++i)
  { latticetypes::Weight column; m->val.column(column,i);
    result->val[i]=weight_to_row(column);
  }
  push_value(result.release());
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
{ explicit tuple_expression(std::vector<expression> l) : list_expression(l)@+{}
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};

@ When we print a tuple display, we just print the component expressions,
enclosed in brackets and separated by commas, to match their input syntax.


@< Function def... @>=
void tuple_expression::print(std::ostream& out) const
{ out << '(';
  if (component.size()==0) out << ')';
  else
    for (size_t i=0; i<component.size(); ++i)
    out << *component[i] << (i<component.size()-1 ? ',' : ')');
}


@ We must now add the relevant cases for the variant |tuple_display| to the
functions |find_type|, |check_type|, and~|convert_expr|. Calling |find_type|
for a tuple display, we must collect the types in a tuple type, which is again
a linked list. This recursive activity most easily done by a recursive
function.

@< Declarations of local functions @>=
type_list_ptr find_type_list(expr_list);

@~Given this function, all we need to do call it and wrap up.

@< Cases for finding the type of other kinds of expressions... @>=
case tuple_display:  return make_tuple_type(find_type_list(e.e.sublist));

@ The definition of this recursive function is straightforward.
@< Function def...@>=
type_list_ptr find_type_list(expr_list l)
{ if (l==NULL) return type_list_ptr(NULL);
  return make_type_list(find_type(l->e),find_type_list(l->next));
}

@ For type checking, we shall need a tuple pattern with any number of unknown
components.

@< Declarations of local functions @>=
type_ptr unknown_tuple(size_t n);

@~This pattern is built up in a simple loop.

@< Function definitions @>=
type_ptr unknown_tuple(size_t n)
{ type_list_ptr tl(NULL);
  while (n-->0) tl=make_type_list(copy(&unknown_type),tl);
  return make_tuple_type(tl);
}


@ If in |check_type| the type required for a tuple display is a tuple type, we
can simply test that all component expressions have the required number and
component types. If not, we call a type error with an instance of
|unknown_tuple|, since we did not in fact type-check the component
expressions.

@< Other cases for testing whether the type of |e| matches |t|... @>=
case tuple_display:
  if (t.kind==tuple_type)
    @< Check that the components of |e.e.sublist| has the number and types
       specified by |t.tuple| @>
  else actual=unknown_tuple(length(e.e.sublist));
break;


@ If we did find a tuple display where a tuple type was required, but the
number of components does not match, we throw a |program_error| rather than a
|type_error| to indicate the problem, since the format of a |type_error| value
is not well suited to this case (but see the next section for a way to use it
anyway). The solution adopted is quite coarse, and it does not single out the
offending tuple. Note that we do signal a |type_error| that is found in the
initial components, even if the number of components should disagree; we only
complain about the number of components if we run out either of component
types or of component expressions to check, without finding any erroneous
types in the preceding components.

@< Check that the components of |e.e.sublist| has the number and types
       specified by |t.tuple| @>=
{ expr_list a=e.e.sublist; type_list l=t.tuple;
  for (; a!=NULL and l!=NULL; a=a->next,l=l->next)
    check_type(l->t,a->e);
  if (a!=NULL or l!=NULL) throw @|
      program_error("Too "+ std::string(a==NULL ? "few" : "many")
		   +" components in tuple");
@.Too few components in tuple@>
@.Too many components in tuple@>
}

@ When converting a tuple expression, we first try to specialise an unknown
type to a tuple with the right number of components; unless the type was
completely undetermined, this just amounts to a test that it is a tuple type
with the right number of components. Here we report a wrong number of
components via a type pattern, which is probably as clear as mentioning too
few or too many components.

@< Other cases for type-checking and converting expressions... @>=
case tuple_display:
{ type_ptr tup=unknown_tuple(length(e.e.sublist));
  if (type.specialise(*tup))
  { std::auto_ptr<tuple_expression> result@|
      (new tuple_expression@|
       (std::vector<expression>(length(e.e.sublist),NULL)));
  @/type_list tl=type.tuple;
    size_t i=0;
    for (expr_list el=e.e.sublist; el!=NULL; el=el->next,tl=tl->next)
      result->component[i++]=convert_expr(el->e,tl->t);
    return result.release();
  }
  else throw type_error(e,tup,copy(&type));
}
break;

@*2 Evaluating tuple displays.
%
Evaluating a tuple display is much like evaluating a list display. Thanks to
dynamically typed values, we can collect the components of a tuple in a vector
without problem. In fact we could reuse the type |row_value| to hold the
components of a tuple, if it weren't for the fact that it would print with
brackets, as a list. Therefore let us trivially derive a new class from
|row_value|.

@< Type definitions @>=
struct tuple_value : public row_value
{ tuple_value(const std::vector<value>& v) : row_value(v) @+{}
  tuple_value* clone() const @+{@; return new tuple_value(*this); }
  void print(std::ostream& out) const;
  static const char* name() @+{@; return "tuple value"; }
private:
  tuple_value(const row_value& v); // copy constructor; used by |clone|
  tuple_value& operator=(const row_value& v); // assignment forbidden
};
@)
typedef std::auto_ptr<tuple_value> tuple_ptr;

@ We just need to redefine the |print| method.
@< Function definitions @>=
void tuple_value::print(std::ostream& out) const
{ if (val.empty()) out << "()";
  else
  { out << '(';
    std::vector<value>::const_iterator p=val.begin();
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
|get<tuple_value>| that will be defined later. We also take care to unlink
each component before handing it to |push_value|, to prevent double ownership
at the moment the stack expands.

@< Function definitions @>=
void push_tuple_components()
{ tuple_ptr tuple(get<tuple_value>());
  for (size_t i=0; i<tuple->length(); ++i)
  { value v=tuple->val[i];
    tuple->val[i]=NULL; // remove component so it remains unshared
    push_value(v); // push component
  }
}

@ We need no auto-pointer in |wrap_tuple|, as shrinking the stack will not
throw any exceptions.

@< Function definitions @>=
void wrap_tuple(size_t n)
{ tuple_value* result=new tuple_value(std::vector<value>(n,NULL));
  while (n-->0) // not |(--n>=0)|, since |n| is unsigned!
    result->val[n]=pop_value();
  push_value(result);
}

@ We could use the same code as for list displays to evaluate a tuple display.
However, since tuples are usually small and we have |wrap_tuple| available, we
can do it more concisely as follows.

@< Cases for evaluating other kinds of expressions... @>=
case tuple_display:
{ for (expr_list l=e.e.sublist; l!=NULL; l=l->next)
    push_value(evaluate(l->e));
  wrap_tuple(length(e.e.sublist));
  result=pop_value();
}
break;

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
struct call_expression : public expression_base
{ Hash_table::id_type function;
  expression argument;
@)
  explicit call_expression(Hash_table::id_type f,expression a)
   : function(f),argument(a) @+{}
  virtual ~call_expression() @+ {@; delete argument; }
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};

@ To print a function call we print the function name and the argument, the
latter enclosed in parentheses unless it is a tuple expression. The latter
condition must be tested by a dynamic cast.

@< Function definitions @>=
void call_expression::print(std::ostream& out) const
{ out << main_hash_table->name_of(function);
  if (dynamic_cast<tuple_expression*>(argument)!=NULL) out << *argument;
  else out << '(' << *argument << ')';
}

@*2 Identifier tables.
As we said above, we need an identifier table to record the types of known
functions (and later other types of identifiers). For the moment we use a
single flat table; there will doubtlessly be need for handling scopes later.
We shall actually store both a type and a value in the table.

@< Type definitions @>=

struct id_data
{ value val; @+ type_declarator* type;
  id_data(value v,type_declarator* t) : val(v),type(t)@+ {}
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
  ~Id_table(); // destructor
@)
  void add(Hash_table::id_type id, value v, type_ptr t); // insertion
  type_declarator* type_of(Hash_table::id_type id) const; // lookup
  value value_of(Hash_table::id_type id) const; // lookup
@)
  size_t size() const @+{@; return table.size(); }
  void print(std::ostream&) const;
};

@ Since we have stored pointers, the destructor must explicitly delete them,
even though we have no immediate plans for destroying identifier tables.

@< Function def... @>=
Id_table::~Id_table()
{ for (map_type::const_iterator p=table.begin(); p!=table.end(); ++p)
  @/{@; delete p->second.val; delete p->second.type; }
}

@ For insertion we must distinguish the case that the key is already present,
since the |insert| method for maps prefers not to overwrite the old value if
there is one, but rather to make no change. It does return in the case  both a
pointer (iterator) to the (key,data) pair that obstructed the insertion and a
boolean failure status, so that we can easily overwrite ``manually'' the old
value.

@< Function... @>=
void Id_table::add(Hash_table::id_type id, value v, type_ptr t)
{ id_data data(v,t.get()); owned_value safe(v);
  std::pair<map_type::iterator,bool> trial
     =table.insert(std::make_pair(id,data));
  safe.release(),t.release(); // no more exception protection needed once here
  if (!trial.second) // then key was present; destroy and replace its data
@/{@; delete trial.first->second.val; delete trial.first->second.type;
    trial.first->second=data;
  }
}

@ In order to have |const| lookup methods, we must refrain from inserting into
the table if the key is not found; we return a null pointer in that case.

@< Function... @>=
type_declarator* Id_table::type_of(Hash_table::id_type id) const
{@; map_type::const_iterator p=table.find(id);
  return p==table.end() ? NULL : p->second.type;
}
value Id_table::value_of(Hash_table::id_type id) const
{@; map_type::const_iterator p=table.find(id);
  return p==table.end() ? NULL : p->second.val;
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

@< Cases for finding the type of other kinds of expressions... @>=
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
  check_type(f_type->func->arg_type,e.e.call_variant->arg);
  return copy(&f_type->func->result_type);
}

@ When a function call appears where a fixed type is expected, we first test
that the function returns this type, and then go on to check the types of the
arguments.

@< Other cases for testing whether the type of |e| matches |t|... @>=
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
  if (f_type->func->result_type!=t)
    actual=copy(&f_type->func->result_type); // signal an error
  else check_type(f_type->func->arg_type,e.e.call_variant->arg);
}
break;

@ In the function |convert_expr| we first get the function from the identifier
table, test if it is known and of function type, then type-check and convert
the argument expression, and build a converted function call~|call|. Finally
we test if the required type matches the return type (in which case we simple
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
  expression call=new call_expression @|
     (e.e.call_variant->fun, convert_expr(e.e.call_variant->arg,
                                          f_type->func->arg_type));
  expression_ptr p(call); // temporarily owns converted call
  if (type.specialise(f_type->func->result_type))
    return p.release();
  else if (coerce(f_type->func->result_type,type,call))
    {@; p.release(); return call; }
  else throw type_error(e,copy(&f_type->func->result_type),copy(&type));
}

@*2 Evaluating function calls.
%
The evaluation of a function call is usually demanded explicitly by the user,
but |check_type| may have inserted calls of conversion functions that will
also be executed by the part of |evaluate| given below (the result of
|conv_expr| however does not encode conversions by ordinary function calls).
In either case what is really executed is a ``wrapper function'', that usually
consists of a call to a library function sandwiched between unpacking and
repacking statements; in some simple cases a wrapper function may decide to do
the entire job itself.

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
  ~builtin_value()@+ {} // don't try to destroy the function pointed to!
  virtual void print(std::ostream& out) const
  @+{@; out << ':' << print_name << ':'; }
  builtin_value* clone() const @+{@; return new builtin_value(*this); }
  static const char* name() @+{@; return "built-in function"; }
private:
  builtin_value(const builtin_value& v)
  : val(v.val), print_name(v.print_name)
  {} // copy constructor, used by |clone|
};
@)
typedef std::auto_ptr<builtin_value> builtin_ptr;


@ Finally we can say what to do when a function call is requested. Recall that
for the moment the function expression is always an identifier, so the
associated (wrapper function) value can be found in the |global_id_table|.

@< Cases for evaluating other kinds of expressions... @>=
case function_call:
{ push_value(evaluate(e.e.call_variant->arg));
    // evaluate and push argument
  std::string name=main_hash_table->name_of(e.e.call_variant->fun);
    // for error messages
  value f_val=global_id_table->value_of(e.e.call_variant->fun);
  if (f_val==NULL) throw
    std::logic_error("built-in function absent: "+name);
@.built-in function absent@>
  force<builtin_value>(f_val)->val();
  // call the wrapper function, leaving result on the stack
  result=pop_value(); // get the result to return it from |evaluate|
}
break;

@ To evaluate a |call_expression| object we evaluate the arguments, get and
check the wrapper function, and call it. Since values are handled via the
execution stack, we don't see them at all in this code.

@< Function definitions @>=
void call_expression::evaluate() const
{ argument->evaluate(); // push evaluated argument on stack
  value f_val=global_id_table->value_of(function);
  if (f_val==NULL) throw std::logic_error("built-in function absent");
@.built-in function absent@>
  force<builtin_value>(f_val)->val();
  // call the wrapper function, leaving result on the stack
}

@*1 Wrapper functions.
%
We have not defined any wrapper functions yet, and therefore have nothing in
the |global_id_table|. Let us start with some preparations for the former.
Wrapper functions will routinely have to do some unpacking of values, which
will involve dynamically casting the |value| to the type they are known to
have because we passed the type checker; should the cast fail we shall throw a
|std::logic_error|. We define a template function to do the unpacking for
every type derived from |value_base|. It is very much like the template
function |force|, but we pull the value off the stack, and therefore take care
to clean it up in case we have to signal a |logic_error|

@< Template and inline function definitions @>=
template <typename D> // |D| is a type derived from |value_base|
 D* get() throw(std::logic_error)
{ value v=pop_value();
  D* p=dynamic_cast<D*>(v);
  if (p==NULL)
  @/{@; delete v;
      throw std::logic_error(std::string("Argument is no ")+D::name());
    }
  return p;
}
@.Argument is no ...@>

@ Here are our first wrapper functions. The function |id_wrapper| is a trivial
function, but it will be put into the |main_id_table| under different names
and signatures, each forcing a particular argument type.

@< Function definitions @>=
void intlist_to_weight_wrapper ()
{@; push_value(new vector_value(cast_intlist_to_weight(pop_value())));
}

void intlistlist_to_latmat_wrapper ()
{@; push_value(new matrix_value(cast_intlistlist_to_latmat(pop_value())));
}

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
{ global_id_table->add(main_hash_table->match_literal(name)
                      ,new builtin_value(f,name),make_type(type));
}

@~Although they would not get type-checked since they are
inserted automatically, we do put types in place for the wrapper functions
|intlist_to_weight_wrapper| and |intlistlist_to_latmat_wrapper| defined above.

@< Function definitions @>=
void initialise_evaluator()
{ execution_stack.reserve(16); // avoid some early reallocations
@)
  install_function(intlist_to_weight_wrapper,">vec<[int]:","([int]->vec)");
  install_function(intlistlist_to_latmat_wrapper
                  ,">mat<[vec]:","([vec]->mat)");
  install_function(id_wrapper,"vec","(vec->vec)");
  install_function(id_wrapper,"mat","(mat->mat)");
@/@< Installation of other built-in functions @>
}

@*1 Invoking the type checker.
%
Let us recapitulate what will happen. The parser will read what the user
types, and returns an |expr| value. Then we shall call |find_type| or
|convert_expr| for this value, and if this does not throw any exceptions, we
are ready to call |evaluate| (for the |expr| value that may have been modified
by the insertion of conversion functions) of the |evaluate| method of the
|expression| value returned by |convert_expr|; after this main program will
print the result. The invocation of the type checker is done by a call to
|analyse_types|. Its signature allows it to be based either on |find_type| or
|convert_expr|; in the former case the pointer~|p| will be made null.

@< Declarations of exported functions @>=
type_ptr analyse_types(const expr& e,expression& p)
   throw(std::bad_alloc,std::runtime_error);

@~The main reason for having this function this function will give us an
occasion to catch any thrown |type_error| and |program_error| exceptions,
something we did not want to do inside the recursive functions |find_type|,
|check_type| and |convert_expr|; since we cannot return normally from
|analyse_types| in the presence of these errors, we map these errors to
|std::runtime_error|, an exception for which the code that calls us will have
to provide a handler anyway.

The choice between |find_type| and |convert_expr| is currently based on
|verbosity|, so that the user may switch between two forms of evaluation
dynamically; this will however go away soon (only |convert_expr| will remain).

@< Function definitions @>=
type_ptr analyse_types(const expr& e,expression& p)
  throw(std::bad_alloc,std::runtime_error)
{ try
  { if (verbosity==0) {@; p=NULL; return find_type(e); }
    type_ptr type=copy(&unknown_type);
    p=convert_expr(e,*type);
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
  { std::cerr << err.what() <<
          " in expression '" << e << "'\n";
  }
  throw std::runtime_error("Type check failed");
@.Type check failed@>
}

@*1 Some built-in functions.
%
We shall now introduce some real built-in functions, starting with integer
arithmetic. Arithmetic operators are implemented by wrapper functions with two
integer arguments (or one in the case of unary minus). Note that the values
are pulled from the stack in reverse order, which is important for the
non-commutative operations like `|-|', `|/|' and~`|%|'. Since values are not
shared, we reuse one the first value object and destroy the second.

@< Definition of other wrapper functions @>=
void plus_wrapper ()
{ push_tuple_components();
  int_ptr j(get<int_value>()); int_ptr i(get<int_value>());
  i->val+=j->val;
  push_value(i);
}
@)
void minus_wrapper ()
{ push_tuple_components();
  int_ptr j(get<int_value>()); int_ptr i(get<int_value>());
  i->val-=j->val;
  push_value(i);
}
@)
void times_wrapper ()
{ push_tuple_components();
  int_ptr j(get<int_value>()); int_ptr i(get<int_value>());
  i->val*=j->val;
  push_value(i);
}
@)
void divide_wrapper ()
{ push_tuple_components();
  int_ptr j(get<int_value>()); int_ptr i(get<int_value>());
  if (j->val==0) throw std::runtime_error("Division by zero");
  i->val/=j->val;
  push_value(i);
}
@)
void modulo_wrapper ()
{ push_tuple_components();
  int_ptr  j(get<int_value>()); int_ptr i(get<int_value>());
  if (j->val==0) throw std::runtime_error("Modulo zero");
  i->val%=j->val;
  push_value(i);
}
@)
void unary_minus_wrapper ()
@+{@; int_value* i=get<int_value>(); i->val =-i->val; push_value(i); }
@)
void divmod_wrapper ()
{ push_tuple_components();
  int_ptr j(get<int_value>()); int_ptr i(get<int_value>());
  if (j->val==0) throw std::runtime_error("DivMod by zero");
  int mod=i->val%j->val;
  i->val/=j->val; j->val=mod;
@/push_value(i); push_value(j); wrap_tuple(2);
}

@ Here is a simple function that outputs a string without its quotes, but with
a terminating newline. This is the first place in this file where we produce
output to a file. In general, rather than writing directly to |std::cout|, we
shall pass via a pointer whose |output_stream| value is maintained in the main
program, so that redirecting output to a different stream can be easily
implemented.

@< Definition of other wrapper functions @>=
void print_wrapper ()
{ string_ptr s(get<string_value>());
  *output_stream << s->val << std::endl;
  wrap_tuple(0); // don't forget to return a value
}

@ We install the function above. The names of  the arithmetic operators
correspond to the ones used in the parser definition file \.{parser.y}.

@< Installation of other built-in functions @>=
install_function(plus_wrapper,"+","(int,int->int)");
install_function(minus_wrapper,"-","(int,int->int)");
install_function(times_wrapper,"*","(int,int->int)");
install_function(divide_wrapper,"/","(int,int->int)");
install_function(modulo_wrapper,"%","(int,int->int)");
install_function(unary_minus_wrapper,"-u","(int->int)");
install_function(divmod_wrapper,"/%","(int,int->int,int)");
install_function(print_wrapper,"print","(string->)");


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
matrix object (which implies copying its contents). In |transpose_mat_wrapper|
we can in fact return the same |matrix_value| that held the argument, but this
does not absolve us from using an auto-pointer while |transpose| is active. We
also define |diagonal_wrapper|, a slight generalisation of |id_mat_wrapper|
that produces a diagonal matrix from a vector.

@< Definition of other wrapper functions @>=
void id_mat_wrapper ()
{ int_ptr i(get<int_value>());
  matrix_ptr m
     (new matrix_value(latticetypes::LatticeMatrix()));
  identityMatrix(m->val,std::abs(i->val)); push_value(m);
}
@)
void transpose_mat_wrapper ()
{@; matrix_ptr m(get<matrix_value>());
  m->val.transpose(); push_value(m);
}
@)
void diagonal_wrapper ()
{ vector_ptr d(get<vector_value>());
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
  vector_ptr v(get<vector_value>());
  matrix_ptr m(get<matrix_value>());
  if (m->val.numColumns()!=v->val.size())
    throw std::runtime_error(std::string("Size mismatch ")@|
     + num(m->val.numColumns()) + ":" + num(v->val.size()) + " in mv_prod");
  vector_ptr w@|
    (new vector_value(latticetypes::Weight(m->val.numRows())));
  m->val.apply(w->val,v->val);
  push_value(w);
}
@)
void mm_prod_wrapper ()
{ push_tuple_components();
  matrix_ptr r(get<matrix_value>()); // right factor
  matrix_ptr l(get<matrix_value>()); // left factor
  if (l->val.numColumns()!=r->val.numRows())
  { std::ostringstream s;
    s<< "Size mismatch " << l->val.numColumns() << ":" << r->val.numRows()
    @| << " in mm_prod";
    throw std::runtime_error(s.str());
  }
  l->val*=r->val;
  push_value(l);
}

@ Here is finally the Smith normal form algorithm. We provide both the
invariant factors and the rewritten basis on which the normal for is assumed,
as separate functions, and to illustrate the possibilities of tuples, the two
combined into a single function.

@h "smithnormal.h"

@< Definition of other wrapper functions @>=
void invfact_wrapper ()
{ matrix_ptr m(get<matrix_value>());
  size_t nr=m->val.numRows();
  latticetypes::WeightList b; @+ matrix::initBasis(b,nr);
  vector_ptr inv_factors
     @| (new vector_value(latticetypes::Weight(0)));
  smithnormal::smithNormal(inv_factors->val,b.begin(),m->val);
  push_value(inv_factors);
}
@)
void Smith_basis_wrapper ()
{ matrix_ptr m(get<matrix_value>());
  size_t nr=m->val.numRows();
  latticetypes::WeightList b; @+ matrix::initBasis(b,nr);
  latticetypes::Weight inv_factors(0);
  smithnormal::smithNormal(inv_factors,b.begin(),m->val);
  latticetypes::LatticeMatrix new_basis(b); // convert basis into matrix
  m->val.swap(new_basis);
  push_value(m);
}
@)

void Smith_wrapper ()
{ matrix_ptr m(get<matrix_value>());
  size_t nr=m->val.numRows();
  latticetypes::WeightList b; @+ matrix::initBasis(b,nr);
  vector_ptr inv_factors
     @| (new vector_value(latticetypes::Weight(0)));
  smithnormal::smithNormal(inv_factors->val,b.begin(),m->val);
  latticetypes::LatticeMatrix new_basis(b); // convert basis into matrix
  m->val.swap(new_basis);
@/push_value(m); push_value(inv_factors); wrap_tuple(2);
}

@ Here is one more wrapper function that uses the Smith normal form algorithm,
but behind the scenes, namely to invert a matrix. Since this cannot be done in
general over the integers, we return an integral matrix and a common
denominator to be applied to all coefficients.
@< Definition of other wrapper functions @>=
void invert_wrapper ()
{ matrix_ptr m(get<matrix_value>());
  if (m->val.numRows()!=m->val.numColumns())
  { std::ostringstream s;
    s<< "Cannot invert a " @|
     << m->val.numRows() << "x" << m->val.numColumns() << " matrix";
    throw std::runtime_error(s.str());
  }
  int_ptr denom(new int_value(0));
  m->val.invert(denom->val);
@/push_value(m); push_value(denom); wrap_tuple(2);
}

@ We must not forget to install what we have defined
@< Installation of other built-in functions @>=
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
We can use identifiers, though currently the definition of identifiers is only
possible at the global level, and is handled separately from expression
evaluation (executing definition statement does call the evaluator, but such a
statement is not a case within the evaluator). So as far as this compilation
unit is concerned, the only point to consider is evaluating applied
identifiers. Here is the type into which applied identifiers are converted.

@< Type definitions @>=
struct identifier_expression : public expression_base
{ Hash_table::id_type code;
@)
  explicit identifier_expression(Hash_table::id_type id) : code(id) @+{}
  virtual ~identifier_expression() @+ {}
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};

@ To print an applied identifier, we get its name from the main hash table.

@< Function definitions @>=
void identifier_expression::print(std::ostream& out) const
{@; out<< main_hash_table->name_of(code); }

@ To evaluate an applied identifier, we just pick the value from the
identifier table, not forgetting to duplicate it, since values are destroyed
after being used in evaluation. Actually implementing this the first time was
more work than it would seem, since we had to introduce the |clone| virtual
method of |value_base| (and types derived from it) in order to do it.

@< Cases for evaluating other kinds of expressions... @>=
case applied_identifier:
{ value p=global_id_table->value_of(e.e.identifier_variant);
  if (p==NULL) throw std::logic_error
  @|   ("Identifier without value:"
	+main_hash_table->name_of(e.e.identifier_variant));
@.Identifier without value@>
  result=p->clone();
}
break;

@ The method |identifier_expression::evaluate| follows the same scheme.

@< Function definitions @>=
void identifier_expression::evaluate() const
{ value p=global_id_table->value_of(code);
  if (p==NULL) throw std::logic_error
  @|   ("Identifier without value:"+main_hash_table->name_of(code));
@.Identifier without value@>
  push_value(p->clone());
}

@ For finding the type, matters are quite the same: we just copy a type
plucked from the table.

@< Cases for finding the type of other kinds of expressions... @>=
case applied_identifier:
{ type_declarator* it=global_id_table->type_of(e.e.identifier_variant);
  if (it==NULL) throw program_error
  @|   ("Undefined identifier "
	+main_hash_table->name_of(e.e.identifier_variant));
@.Undefined identifier@>
  return copy(it);
}

@ For an applied identifier there is not much difference in type checking for
the case where a definite type is required: the type found in the table must
now equal the specified one. Only we must cater for the fact that the types
only match after conversion. In this case we cannot push down the conversion
into the expression as we did for list displays so we must handle simple and
composite conversions.

@< Other cases for testing whether the type of |e| matches |t|... @>=
case applied_identifier:
{ type_declarator& it=*global_id_table->type_of(e.e.identifier_variant);
  if (&it==NULL) throw program_error
  @|   ("Undefined identifier "
	+main_hash_table->name_of(e.e.identifier_variant));
@.Undefined identifier@>
  if (it!=t)
  { if (t==vec_type && it==row_of_int_type)
      @< Insert a vector conversion at |e| @>
    else if (t==mat_type and
           (it==row_row_of_int_type or it==row_of_vec_type))
      @< Insert a matrix conversion at |e| @>
    else actual=copy(&it); // if no conversion works, signal an error
  }
}
break;

@ For |convert_expr| the logic is similar, but we use the general |coerce|
routine.

@< Other cases for type-checking and converting expressions... @>=
case applied_identifier:
{ type_declarator& it=*global_id_table->type_of(e.e.identifier_variant);
  if (&it==NULL) throw program_error
  @|   ("Undefined identifier "
	+main_hash_table->name_of(e.e.identifier_variant));
@.Undefined identifier@>
  expression id=new identifier_expression(e.e.identifier_variant);
  expression_ptr p(id); // temporarily owns |id|
  if (type.specialise(it)) return p.release();
  else if (coerce(it,type,id)) {@; p.release(); return id; }
  else throw type_error(e,copy(&it),copy(&type));
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
  row_subscription(expression a, expression i) : array(a),index(i) @+{}
  virtual ~row_subscription() @+ { delete array; delete index; }
  virtual void evaluate() const;
  virtual void print(std::ostream& out) const;
};
@)
struct vector_subscription : row_subscription
{ vector_subscription(expression a, expression i)
  : row_subscription(a,i) @+{}
  virtual void evaluate() const;
};
@)
struct matrix_subscription : row_subscription
{ matrix_subscription(expression a, expression ij)
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

@ Here is the case for subscriptions in |find_type|; we only treat the case of
subscription from a row, since this implementation is of temporary use anyway.
Upon success we avoid copying the component type from the row type found,
since |array_type| is going to destruct its component just after; rather we
copy the pointer into an auto-pointer and then reset the original pointer
to~|NULL|.

@< Cases for finding the type of other kinds of expressions... @>=
case subscription:
{ type_ptr array_type=find_type(e.e.subscription_variant->array);
  type_ptr index_type=find_type(e.e.subscription_variant->index);
  if (*index_type!=int_type)  throw type_error
      (e.e.subscription_variant->index,index_type,copy(&int_type));
  else if (array_type->kind!=row_type) throw type_error
      (e.e.subscription_variant->array,array_type,copy(&row_of_type));
  else
  { type_ptr result_type(array_type->component_type);
    array_type->component_type=NULL;
    return result_type;
  }
}

@ When a subscription appears where a fixed type is expected, one might
imagine testing the index first and then imposing the corresponding row-of
type on the array expression. However, this is not a good idea in general,
since subscription from different aggregate types, such as \.{[int]} and
\.{vec} might produce the same type as result, so we cannot reliably predict
the type of the array expression. Therefore this code proceeds along the same
lines as that of the previous section.

@< Other cases for testing whether the type of |e| matches |t|... @>=
case subscription:
{ type_ptr array_type=find_type(e.e.subscription_variant->array);
  type_ptr index_type=find_type(e.e.subscription_variant->index);
  if (*index_type!=int_type)  throw type_error
      (e.e.subscription_variant->index,index_type,copy(&int_type));
  else if (array_type->kind!=row_type) throw type_error
      (e.e.subscription_variant->array,array_type,copy(&row_of_type));
  else if (*array_type->component_type!=t)
   // then signal type error for entire subscription
  @/{@; actual=type_ptr(array_type->component_type);
    array_type->component_type=NULL;
  }
}
break;

@ When encountering a subscription in |convert_expr|, we determine the types
of array and index expression separately. Then we look if the types agree with
any of the four types of subscription expressions that we can convert to.
If none match,
 or to signal an error if none
applies.

@< Other cases for type-checking and converting... @>=
case subscription:
{ type_declarator array_type, index_type, subscr_type;
    // initialised to |undetermined_type|
  expression array=
    convert_expr(e.e.subscription_variant->array,array_type);
  expression index=
    convert_expr(e.e.subscription_variant->index,index_type);
  expression subscr=NULL;
  if (array_type.kind==row_type) // a row subscription
    if (index_type!=int_type) throw type_error
      (e.e.subscription_variant->index,copy(&index_type),copy(&int_type));
    else
    { type_ptr comp(array_type.component_type);
      array_type.component_type=NULL; // avoid double ownership
      subscr_type.set_from(comp); // incorporate component type
      subscr=new row_subscription(array,index);
    }
  else if (array_type==vec_type)
    if (index_type!=int_type) throw type_error
      (e.e.subscription_variant->index,copy(&index_type),copy(&int_type));
    else
    @/{@; subscr_type.specialise(int_type);
      subscr= new vector_subscription(array,index);
    }
  else if (array_type==mat_type)
    if (index_type!=int_int_type) throw type_error
      (e.e.subscription_variant->index,copy(&index_type),copy(&int_int_type));
    else
    @/{@; subscr_type.specialise(int_type);
      subscr=new matrix_subscription(array,index);
    }
  else throw type_error // array expression is not of any aggregate type
      (e.e.subscription_variant->array,copy(&array_type),copy(&row_of_type));
@)
  expression_ptr p(subscr);  // temporarily owns |subscr|
  if (type.specialise(subscr_type)) return p.release();
  else if (coerce(subscr_type,type,subscr))
   {@; p.release(); return subscr; }
  else throw type_error(e,copy(&subscr_type),copy(&type));
}


@ Evaluating subscriptions is currently very inefficient, since it requires a
copy of the array object to be evaluated, of which all but one component will
be thrown away. Later we should optimise at least the case where the array
object is given by an applied identifier, to avoid the copy.

@< Cases for evaluating other kinds of expressions... @>=
case subscription:
{ row_ptr array(force<row_value>(evaluate(e.e.subscription_variant->array)));
  int_ptr index(force<int_value>(evaluate(e.e.subscription_variant->index)));
  if (static_cast<unsigned int>(index->val)>=array->val.size())
    throw std::runtime_error("Index "+num(index->val)+" out of range");
  result=array->val[index->val]->clone();
}
break;

@ Here are the |evaluate| methods for the various subscription expressions.
For |matrix_subscription|, note that |push_tuple_components()| takes care of
giving access to the individual index values without ownership conflict (by
the time |j| and |i| are accessed, the tuple is already destroyed, but its
components survive).

@< Function definitions @>=
void row_subscription::evaluate() const
{ int_ptr i((index->evaluate(),get<int_value>()));
  row_ptr r((array->evaluate(),get<row_value>()));
  if (static_cast<unsigned int>(i->val)>=r->val.size())
    throw std::runtime_error("index "+num(i->val)+" out of range");
  push_value(r->val[i->val]->clone());
}
@)
void vector_subscription::evaluate() const
{ int_ptr i((index->evaluate(),get<int_value>()));
  vector_ptr v((array->evaluate(),get<vector_value>()));
  if (static_cast<unsigned int>(i->val)>=v->val.size())
    throw std::runtime_error("index "+num(i->val)+" out of range");
  push_value(new int_value(v->val[i->val]));
}
@)
void matrix_subscription::evaluate() const
{ index->evaluate(); push_tuple_components();
  int_ptr j(get<int_value>());
  int_ptr i(get<int_value>());
  matrix_ptr m((array->evaluate(),get<matrix_value>()));
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

@*1 Built-in types defined elsewhere.
In order for this compilation unit to function properly, it must know of the
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


@*1 Making global definitions.
For the moment, applied identifiers can only get their value through the
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
caught, the identifier list |ids| and the expression |e| should not be
destroyed here, since the parser which aborts after calling this function
should do that while clearing its parsing stack.

@< Function definitions @>=
extern "C"
void global_set_identifier(expr_list ids, expr rhs)
{ using namespace atlas::interpreter; using namespace std;
  try
  { expression e;
    type_ptr t=analyse_types(rhs,e);
    if (ids->next!=NULL)
      @< Check that identifiers are distinct and that |t| is an appropriate
         tuple type; if not, |throw| a |runtime_error| @>
    value v= e==NULL ? evaluate(rhs) :(e->evaluate(),pop_value());
    if (ids->next==NULL)
    { cout << "Identifier " << ids->e << ": " << *t << std::endl;
      global_id_table->add(ids->e.e.identifier_variant,v,t); // releases |t|
    }
    else @< Perform a multiple assignment @>
  }
  catch (runtime_error& err)
  { cerr << err.what() << ", identifier" << (ids->next!=NULL ? "s " :" ");
    for (expr_list l=ids; l!=NULL; l=l->next)
      cerr << main_hash_table->name_of(l->e.e.identifier_variant)
           << (l->next!=NULL?",":"");
    cerr << " not assigned to.\n";
    clear_execution_stack();
  }
  catch (logic_error& err)
  { cerr << "Unexpected error: " << err.what() << ", evaluation aborted.\n";
    clear_execution_stack(); main_input_buffer->close_includes();
  }
  catch (exception& err)
  { cerr << err.what() << ", evaluation aborted.\n";
    clear_execution_stack(); main_input_buffer->close_includes();
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
        ("Repeated identifier "
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
{ tuple_ptr tv(force<tuple_value>(v));
  cout << "Identifiers ";
  size_t i=0; type_list tl=t->tuple;
  for (expr_list l=ids; l!=NULL; l=l->next,++i,tl=tl->next)
  { cout << l->e << ": " << tl->t << ( l->next!=NULL ? ", " : ".\n");
    global_id_table->
      add(l->e.e.identifier_variant,tv->val[i],copy(&tl->t));
    tv->val[i]=NULL; // ensure value in table is unshared
  }
}

@*1 Printing type information.
%
It is useful to print type information, either for a single expression or for
all identifiers in the table. We declare the pointer that was already used in
|print_wrapper|.

@< Declarations of global variables @>=
extern std::ostream* output_stream;

@ The |output_stream| will normally point to |std::cout|.

@< Global variable definitions @>=
std::ostream* output_stream= &std::cout;

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
  { expression p;
    *output_stream << "type: " << *analyse_types(e,p) << std::endl;
  }
  catch (std::runtime_error& err) { std::cerr<<err.what()<<std::endl; }
}

@ The function |show_ids| prints a table of all known identifiers and their
types. It does this by copying the identifiers from the global identifier
table into a |map| object, and then traversing that object to print their
types.

@h <map>
@< Function definitions @>=
extern "C"
void show_ids()
{@; *output_stream << *global_id_table;
}


@* Index.

% Local IspellDict: british
