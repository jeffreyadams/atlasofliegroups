\def\emph#1{{\it#1\/}}

@* Evaluating expressions.
This file describes the unit \.{evaluator}, which implements an evaluator of
expressions, that were produced by the parser with the help of the types and
functions defined in the unit \.{parsetree}.

@h "evaluator.h"
@c
namespace atlas { namespace interpreter {
@< Global variable definitions @>@;
@< Declarations of local functions @>@;
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
}@; }@;
#endif

@*1 Types to represent types.
Before we start writing the evaluator, we must consider the question of
types. Certainly evaluating will involve values of different types, and we do
not want to rewrite every part of the evaluator separately for every type of
value handled, which will in fact not be limited to a fixed finite set of
types. Therefore there is hardly any choice but to handle values as generic
pointers, and to accompany them by some (structured) typing information that
will allow these pointers to be converted into more specific ones when it
comes to actually accessing the data representing the values.

At the very least we want to handle sequences of values, for instance for
vectors and matrices. Often these will be uniform sequences, all of the same
type, but we also want to allow tuples of mixed type. In the former case we
can describe the type in a more compact way than the value itself, for
instance ``row of integers'', while in the latter case the type has as many
components as the value itself. A fundamental question is whether we shall
store values and their type information in the same structure or in distinct
ones. The former may seem more natural, certainly in the case of tuples, but
the latter seems more flexible when types are predictable such as in handling
uniform sequences, and it allows for the possibility that type information can
be computed before actual evaluation starts; we shall adapt the latter
approach.

@< Include... @>=
#include <memory>

@~ So let us now consider how to represent types, which will be small tree
structures. For now we do not attempt sharing of subtrees between types, so a
root of a type expression will own its tree. To avoid useless copying, every
root must be held in an auto-pointer, so that proper destruction can take
place at exceptions, but internally types use ordinary pointers. We also
provide a |type_list| type to be used in future variants.

@< Type definitions @>=
struct type_declarator;
typedef std::auto_ptr<type_declarator> type_ptr;
typedef struct type_node* type_list;
typedef std::auto_ptr<type_node> type_list_ptr;
struct type_node
{ type_declarator* t; @+ type_list next;
  type_node(type_declarator* t, type_list l) : t(t),next(l) @+{}
  type_node(const type_node& n);
  ~type_node();
};
@< Other |struct| definitions @>@;

@ The |type_list| type gives a preview of matters that will apply to types in
general, since types and type lists own their trees. The copy constructor must
make a full copy of the list, which happens recursively by allocation fresh
objects for the ones directly pointed to, while copy-constructing their
contents. We first define a function |copy| that will be used in most cases
where a call of the copy constructor for |type_declarator| is desired. It
converts a pointer to a |type_declarator| (which is assumed not to be owned),
into an auto-pointer to a copied instance of that type declarator.

@< Function definitions @>=
type_ptr copy(const type_declarator* t) @+
{@; return type_ptr(new type_declarator(*t)); }

@ We copy-constructing a type list, we must take care that after a first
allocation has succeeded, the second allocation or the construction of the
object in it may throw an error. In order to ensure the proper destruction of
the result of the initial allocation in this case we hold the result in a
temporary auto-pointer, provided by~|copy|. For the second (final) allocation
this precaution would be useless: either no problem arises and a returned
auto-pointer would be released immediately after creation, or a problem arises
before the return from |new| and no pointer is ever returned; in this case the
cleanup is done by the calls that are aborted (|new| will deallocate if its
allocation succeeds but the construction---here the recursive
copy-construction implicit in |*n.next|--- throws an exception; any nodes of a
half-constructed list will therefore be deallocated). Note that it is
necessary in the code below that there are no calls of |new| after that of any
|release|, which is why we the |next| field is filled below before the |t|
field is.

@< Function definitions @>=
type_node::type_node(const type_node& n)
{ type_ptr head=copy(n.t);
  next=n.next==NULL ? NULL : new type_node(*n.next);
  t=head.release();
}

@ The destructor for |type_node| will possibly be called implicitly by the
destruction of a |type_list_ptr|. The destructor for variants containing a
|type_list| should call |delete| for that pointer, which will also call the
destructor below, which will recursively clean up the whole list.

@< Function def... @>=
type_node::~type_node()
@+{@; delete t; delete next; }

@ The constructors for a type list are quite simple, thanks to the fact that
their arguments are passed as auto-pointers. Types and type lists are built up
from the leaves to the root, and any object holding a type (list) must have a
destructor that destroys the entire tree. When installing a descendant into a
fresh node for the parent by the calling the constructor, the temporary
holding the descendent would destroy its subtree just afterwards, forcing a
full copy to made by the copy constructor to be used for incorporation into
the parent node, a clear waste (which would apply repeatedly during the
construction process). By using auto-pointers we can pass the responsibility
for destroying the descendent first to the constructor, which will pass the
responsibility to the parent once the descendent is inserted; the destruction
of the temporary in the function calling the constructor will do nothing.

We in fact pass the auto-pointers by value rather than by (non-constant)
reference to functions such as |make_type_list|, which could give a small
overhead of using the copy constructor for auto-pointers at the call,
transferring ownership, followed just afterwards by a transfer of ownership by
the |release| method. This declaration however allows |make_list_type| to be
called directly with the result of a constructor or construction function (and
in these cases no copy constructor is called), where otherwise intermediate
storage in a named variable would be required to provide an ``lvalue'' for a
non-constant reference in the call.

In these construction functions we are prudent to release the auto-pointer
only \emph{after} the allocation for |new| is completed; the order would be
undefined if we would write |new type_node(t.release(),l.release())|. If the
standard had specified that arguments to a constructor in |new| are only to be
evaluated after the allocation is completed, we could have used that more
fluent idiom below and in similar situations; unfortunately we are forced to
do it differently in such cases.

@< Function def... @>=

type_list_ptr make_type_list(type_ptr t,type_list_ptr l)
{ type_node* p=new type_node(t.get(),l.get());
  t.release(),l.release(); return type_list_ptr(p);
}
type_list_ptr make_type_singleton(type_ptr t)
{@; type_node* p=new type_node(t.get(),NULL); t.release();
  return type_list_ptr(p);
}


@ Now we must define type declarators themselves. We start with a very simple
type model, but provide hooks for extension. For now, types are either
primitive ones, of which there are a finite number and which can therefore be
represented by an enumeration value, or they are ``row of'' some other type
(later we shall add alternatives to this list). We shall therefore use the
type |type_declarator| defined here to represent types (the name indicates
that it expresses the structure rather than an abbreviated name for the type).
By using an anonymous union, the field selectors like |prim| of the variants
in the union can be used directly on the level of the structure. The variant
|prim| will apply if |kind==primitive_type|, and the variant |comp_type| will
apply if |kind==row_type|.


@< Type definitions @>=
enum primitive @+{ integral_type, string_type, boolean_type
  , @< Other primitive types @>@;@;
@/ nr_of_primitive_types };
enum type_type @+{ primitive_type, row_type
                 , @< Tags for other kinds of types @>@;@; };
struct type_declarator
{ type_type kind;
  union
  { primitive prim; type_declarator* comp_type;
    @< Variants of the |union| in |type_declarator|
       for other kinds of types@>@;
  };
@)
  explicit type_declarator(primitive p)
    : kind(primitive_type) @+{@; prim=p; }
  explicit type_declarator(type_declarator* c)
    : kind(row_type) @+{@; comp_type=c; }
  @< Other constructors for |type_declarator| @>@;
  type_declarator(const type_declarator& t); // copy constructor
  ~type_declarator();
};

@ Since a |type_declarator| possesses all its subtypes, the copy constructor
must in the case of a |row_type| recursively copy of the component type, and
the destructor must clean up afterwards. Since only one object is pointed to,
there is no need for intermediate auto-pointers here.

@< Function definitions @>=
type_declarator::type_declarator(const type_declarator& t) : kind(t.kind)
{ switch (kind)
  { case primitive_type: prim=t.prim; break;
    case row_type:
      comp_type=new type_declarator(*t.comp_type);
    break;
    @\@<Other cases in the copy constructor for |type_declarator| @>
  }
}

type_declarator::~type_declarator()
{ switch (kind)
  { case primitive_type: break;
    case row_type: delete comp_type; break;
    @\@<Other cases in the destructor for |type_declarator| @>
  }
}


@ The constructing functions below are similar to those for type lists;
similar remarks as above apply.

@< Function definitions @>=
type_ptr make_prim_type(primitive p)
@+{@; return type_ptr(new type_declarator(p)); }

type_ptr make_row_type(type_ptr p)
{@; type_declarator* t=new type_declarator(p.get()); p.release();
    return type_ptr(t);
}

@ In practice we shall rarely call functions like |make_prim_type| and
|make_row_type| directly to make explicit types, since this is rather
laborious. Instead, such explicit types will be constructed by the function
|make_type| that parses a string, and correspondingly calls the appropriate
type constructing functions. Since that requires knowledge of the complete set
of possible types, we defer its definition to the end of this file.

@< Declarations of exported functions @>=
type_ptr make_type(const char* s);

@ For printing types (and later for parsing them) we shall need names for the
primitive ones.

@< Declarations of global variables @>=

extern const char* prim_names[];

@~Here are the ones we already know.
@< Global variable definitions @>=
const char* prim_names[]=
{"int","string","bool",@< Other primitive type names@>@;@; };

@ We shall pass |type_declarator| values to the operator~`|<<|' by constant
reference, which is seems more decent than doing so by pointer (overriding the
definition that simply prints the hexadecimal address); for other types as
well we shall not define an instance of~`|<<|' for their pointers. Since we
often hold types in |type_ptr| values, this does mean the we must dereference
explicitly in printing.

@< Declarations of exported functions @>=
std::ostream& operator<<(std::ostream& out, const type_declarator& t);

@ The cases for printing the types seen so far are very simple.

@< Function definitions @>=

std::ostream& operator<<(std::ostream& out, const type_declarator& t)
{ switch(t.kind)
  { case primitive_type: out << prim_names[t.prim]; break;
    case row_type: out << "[" << *t.comp_type << "]"; break;
  @\@< Other cases for printing types @>
  }
  return out;
}

@ Finally we need a comparison for structural (in)equality of type
declarators.

@< Function definitions @>=
bool operator!= (const type_declarator& x,const type_declarator& y)
{ if (x.kind!=y.kind) return true;
  switch (x.kind)
  { case primitive_type: return x.prim!=y.prim;
    case row_type:
      return *x.comp_type!=*y.comp_type;
    @\@<Other cases for |operator!=| for |type_declarator|,
        all performing |return| @>
  }
  return false; // to keep the compiler from complaining, never reached
}

inline bool operator== (const type_declarator& x,const type_declarator& y)
{@; return !(x!=y); }


@*1 Dynamically typed values.
Now we shall consider the values. We could either use void pointers to
represent generic values and cast them when necessary, or use inheritance and
the dynamic cast feature of \Cpp. We shall try the second option since it
seems more sophisticated, and see how it works out, even though this means
that in reality we have dynamic type information stored in the values as well
as an external |type_declarator| value describing their type.

@< Includes needed in the header file @>=
#include <iostream>
@~We start with a base class for values. There must be at least one virtual
function in the class, which could be just the destructor, but we add a
function for printing. This allows the base class to be defined abstract (one
cannot declare a destructor purely virtual since it will always be called,
after the destructor for a derived class). The printing function will
demonstrate the ease of using dynamic typing via inheritance. It does not even
require any dynamic casting, but other operations on values will. Once
identifiers are introduced, we will need to clone objects of types derived
from |value_base| so we introduce a (purely) virtual |clone| method. We forbid
copying and assigning for the moment (the class is abstract anyway); we'll see
if we can get away with that.

@< Type definitions @>=
struct value_base
{ value_base() @+ {};
  virtual ~value_base() @+ {};
  virtual void print(std::ostream& out) const =0;
  virtual value_base* clone() const =0;
private:
  value_base(const value_base& x); //copying and assigning forbidden
  value_base& operator=(const value_base& x);
};
typedef value_base* value_ptr;

@ We can already make sure that the operator~`|<<|' will do the right thing
for any of our values.

@< Declarations of exported functions @>=
std::ostream& operator<< (std::ostream& out, const value_base& v);

@~The operator~`|<<|' calls the (virtual) |print| method of the object pointed
to, and (trivially but importantly; we forgot to do this initially) returns
the reference to the output stream object. The copy constructor and assignment
constructor for |value_base| are currently forbidden, but we provide loud
versions to easily enable tracking them if this should be needed.

@< Function definitions @>=
std::ostream& operator<< (std::ostream& out, const value_base& v)
@+{@; v.print(out); return out; }
@)
value_base::value_base(const value_base& x)
{@; std::cerr << "Copying " << x << std::endl; }
value_base& value_base::operator=(const value_base& x)
{@; std::cerr << "Assigning " << *this << "<-" << x << std::endl;
  return *this;
}

@ Now we derive the first ``primitive'' value types.

@< Type definitions @>=

struct int_value : public value_base
{ int value;
  explicit int_value(int v) : value(v) @+ {}
  ~int_value()@+ {}
  void print(std::ostream& out) const @+{@; out << value; }
  int_value* clone() const @+{@; return new int_value(*this); }
private:
  int_value(const int_value& v) : value(v.value) @+{}
};

struct string_value : public value_base
{ std::string value;
  explicit string_value(char* t) : value(t) @+ {}
  ~string_value()@+ {}
  void print(std::ostream& out) const @+{@; out << '"' << value << '"'; }
  string_value* clone() const @+{@; return new string_value(*this); }
private:
  string_value(const string_value& v) : value(v.value) @+{}
};

struct bool_value : public value_base
{ bool value;
  explicit bool_value(bool v) : value(v) @+ {}
  ~bool_value()@+ {}
  void print(std::ostream& out) const @+{@; out << std::boolalpha << value; }
  bool_value* clone() const @+{@; return new bool_value(*this); }
private:
  bool_value(const bool_value& v) : value(v.value) @+{}
};

@ For the moment there is one more type derived from |value_base|, namely the
type for ``row of'' types.

@< Includes needed in the header file @>=
#include <vector>

@~Using vectors from the standard template library, the realisation of row
values is quite easy. Since row values hold pointers they own, it is useful to
define a auto-pointer type for them as well.

@< Type definitions @>=
struct row_value : public value_base
{ std::vector<value_ptr> value;
  explicit row_value(const std::vector<value_ptr>& v) : value(v) @+{}
  ~row_value();
  void print(std::ostream& out) const;
  size_t length() const @+{@; return value.size(); }
  row_value* clone() const @+{@; return new row_value(*this); }
protected:
  row_value(const row_value& v);
private:
  row_value& operator=(const row_value& v);
};
typedef std::auto_ptr<row_value> row_ptr;

@ Before we forget it, let us define the copy constructor (needed for the
|clone| method) and the destructor for |row_value| objects. Since the |value|
field contains a vector of pointers, we must explicitly clone respectively
delete the pointed-to objects. We use here the fact (without which life would
be much harder) that |delete| will do the right thing even if it is called
with a pointer to a base class of the class for which the pointer was created.

@< Function definitions @>=
row_value::row_value(const row_value& v) : value(v.value)
{ for (std::vector<value_ptr>::iterator p=value.begin() ; p!=value.end(); ++p)
    *p=(*p)->clone();
}

row_value::~row_value()
{ for (std::vector<value_ptr>::iterator p=value.begin() ; p!=value.end(); ++p)
    delete *p;
}

@ So here is the first occasion where we shall use virtual functions. For the
moment the output routine performs an immediate recursion; later we shall try
to make this more elegant by computing the width needed to output component
values, and adapt the formatting to that.

@< Function definitions @>=
void row_value::print(std::ostream& out) const
{ if (value.empty()) out << "[]";
  else
  { out << '[';
    std::vector<value_ptr>::const_iterator p=value.begin();
    do {@; (*p)->print(out); ++p; out << (p==value.end() ? ']' : ','); }
    while (p!=value.end());
  }
}

@*2 The evaluator.
Before we describe evaluation of expressions we must realise that evaluation
can cause runtime errors. The evaluator may throw exceptions due to
inconsistency of our (rather than the user's) program, which are classified as
|std::logic_error|. It may also throw exceptions due to errors not caught by
the type checker in the user input, such as size mismatch in matrix
operations; in such cases it will throw |std::runtime_error|.

@< Include... @>=
#include "parsetree.h"
#include <stdexcept>

@ The evaluator will be realised as a function that maps expressions produced
by the parser to values of a type derived from |value_type|.

@< Declarations of exported functions @>=
value_ptr evaluate(expr e)
  throw(std::bad_alloc,std::logic_error,std::runtime_error);

@~Evaluating a denotation is trivial, and yields the constant value stored in
the denotation. For other types of expressions we shall give the evaluation in
separate sections. We assign the value to be returned to a variable rather
than calling |return| in the various branches, in order to make the conversion
to |value_ptr| less implicit, and to avoid compiler warnings about ending a
non-void function without |return|.

@h<cstdlib>
@< Function definitions @>=
value_ptr evaluate(expr e)
  throw(std::bad_alloc,std::logic_error,std::runtime_error)
{ value_ptr result=NULL;
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
    case list_display:
      @< Evaluate a list display and set |result| to a pointer to the
	 result @>
      break;
    @\@< Cases for evaluating other kinds of expressions @>
  }
  return result;
}

@ Evaluating a list display is quite straightforward. We first calculate the
length the resulting value will have, then construct a |row_value| object of
the proper size filled with null pointers, that are successively replaced by
the results of evaluating the expressions in the list display. Proceeding in
this order guarantees that if any evaluation should throw an exception, then
the results of previous evaluations will be destroyed.

@< Evaluate a list display... @>=
{ size_t length=0;
  for (expr_list l=e.e.sublist; l!=NULL; l=l->next) ++length;
  row_ptr v@|(new
    row_value(std::vector<value_ptr>(length,static_cast<value_ptr>(NULL))));
  expr_list l=e.e.sublist;
  for (size_t i=0; i<length; ++i,l=l->next)
    v->value[i]=evaluate(l->e);
  result=v.release();
}

@*1 Making more rigid types.
That was remarkably easy (apart from the incomprehensible error messages the
compiler threw at us for clerical errors). However, one must realise that this
evaluator only returns a |value_ptr|, which points to an object that is
effectively dynamically typed, and which cannot be used directly by any
library function. In order to produce for instance an object of type
|std::vector<int>|, it is necessary to assure first of all that all
expressions in a list display do indeed have integral type; once this is done
one can evaluate the expression and confidently convert the object pointed to
by the resulting |value_ptr| into a vector. So we are led to define our
type-checking function now.

Before we do that however, we must realise that type analysis can fail: if for
instance a list display has components with different types, then there is no
way in which we can return a type for the display. We first define a general
exception class |program_error| derived from |exception|, which represents any
kind of error due to user input (for instance a wrong number of arguments,
when later we define function calls).

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
  type_error
    (const expr& e, const type_declarator & a, const type_declarator& r)
    throw() @/
    : program_error("Type error"),offender(e),actual(&a),required(&r) @+{}
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
below we, could have said instead |actual=new type_declarator(*e.actual)|, thus
avoiding the use of one of the two auto-pointers, at the price of some loss of
symmetry and readability.

@< Function definitions @>=
type_error::type_error(const type_error& e)
 : program_error(e), offender(e.offender)
{@; type_ptr p=copy(e.required);
  actual=copy(e.actual).release(); required=p.release();
}


@*2 Finding the type of an expression.
For the moment type analysis is restricted to the function |find_type|, which
examines an expression and returns its type.

@< Declarations of local functions @>=
type_ptr find_type (expr e) throw(std::bad_alloc,program_error);

@~We proceed by a straightforward traversal of the parse tree, with type
information flowing up.
@< Function definitions @>=
type_ptr find_type (expr e) throw(std::bad_alloc,program_error)
{ switch(e.kind)
  { case integer_denotation: return make_prim_type(integral_type);
    case string_denotation: return make_prim_type(string_type);
    case boolean_denotation: return make_prim_type(boolean_type);
    case list_display:
     @< Find the component type~|c| and |return| ``row of''~|c|;
        |throw| a |type_error| if there are unequal component types @>
  @\@< Cases for finding the type of other kinds of expressions @>
  }
  return type_ptr(NULL); // keep the compiler happy, never reached
}

@ If there are no components in a list display, we must still decide on a
component type; we choose it to be integral in this case. In other cases the
first type found will be required in all other positions (this must later be
replaced by a more flexible regime, inserting conversion functions if
required).

@< Find the component type... @>=
{ if (e.e.sublist==NULL) return make_type("[int]");
  type_ptr comp_type=find_type(e.e.sublist->e);
  for (expr_list l=e.e.sublist->next; l!=NULL; l=l->next)
  { type_ptr c2=find_type(l->e);
    if (*c2!=*comp_type)
      throw type_error (l->e,*c2.release(),*comp_type.release());
@.Type error@>
  }
  return make_row_type(comp_type);
}

@*2 Checking if an expression has a given type.
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

The function |check_type| does not return any value, but it may modify its
argument~|e|; if it returns at all without throwing an exception, things are
in order and any necessary conversions are inserted into~|e|. If it throws a
|type_error|, it must copy or construct the types reported, since it does not
own its first argument.

@< Declarations of local functions @>=
void check_type (const type_declarator& t,expr& e)
   throw(std::bad_alloc,program_error);

@~In spite of the fact that a type is specified, our first concern is analyse
the expression, as for |find_type|. We adopt a somewhat peculiar convention
that an error is signalled by setting the type |actual| to a non-null value;
this avoids having to place |throw| expressions in many places. When we do
throw a |type_error| after finding |actual| to be set, we first copy the
required type~|t| for inclusion in the error object. After this the
construction of the |type_error| involves no copying, since the
|type_declarator| objects are passed by constant reference, and a pointer to
them is included in the |type_error|. We do make sure that the copy is
complete before constructing the |type_error|, since a possible |bad_alloc|
during the copy should not risk finding the pointer from |actual| already
released.

@< Function definitions @>=
void check_type (const type_declarator& t,expr& e)
   throw(std::bad_alloc,program_error)
{ static type_declarator vect(vector_type);
  static type_declarator matr(matrix_type);
   // precompute these ``primitive'' types
  type_ptr actual(NULL);
  switch(e.kind)
  { case integer_denotation:
      if (t.kind!=primitive_type || t.prim!=integral_type)
        actual= make_prim_type(integral_type);
    break;
    case string_denotation:
      if (t.kind!=primitive_type || t.prim!=string_type)
        actual= make_prim_type(string_type);
    break;
    case boolean_denotation:
      if (t.kind!=primitive_type || t.prim!=boolean_type)
        actual= make_prim_type(boolean_type);
    break;
  @/@< Other cases for testing whether the type of |e| matches |t| @>
  }
  if (actual.get()!=NULL)
  { type_ptr p=copy(&t); // ensure copy is made before throwing
    throw type_error(e,*actual.release(),*p.release());
@.Type error@>
  }
}

@ In the case of a list display we encounter some new kinds of considerations.
Firstly it could happen that the required type is such that a list display
can never yield it, nor any type that can be converted to it. In that was we
can signal a type error, but we do not know the offending type yet. We could
figure this out by calling |find_type|, but this is somewhat wasted effort,
and worse it could provoke another type error which would prevent us from
reporting the one we found first. So we shall prefer to report the type as
``row of something'' without searching further. We therefore introduce a new
``primitive type'' that will not correspond to any real values, but
represents ``some undetermined type'' in error messages.
@< Other primitive types @>=
undetermined_type, @[@]

@~The type will be printed as a |"*"|; a short name will prove to be
practical.
@< Other primitive type names @>=
"*", @[@]

@ If the type required for a list display is a row type, we can simply test
that all component expressions have the required component type. If not, we
encounter the second new consideration, namely that some ``primitive'' types
such as \.{vec} and \.{mat} may be obtained after conversion from a list type.
If this does not work, we shall call a type error with actual type \.{[*]}.

@< Other cases for testing whether the type of |e| matches |t| @>=
case list_display:
  if (t.kind==row_type)
    for (expr_list l=e.e.sublist; l!=NULL; l=l->next)
      check_type(*t.comp_type,l->e);
  else
  { @< If the type |t| can be converted from a list type, insert the
       conversion into |e|, and check the types of the component expressions
       in |e.e.sublist|, then |break| @>
    actual=make_type("[*]");
  }
break;

@ Here is where our implicit conversions from list displays take place. We
provide conversion from a list of integers to a vector, and from a list of
vectors to a matrix. The conversion is realised by inserting a function call
into the expression; later a more complete solution of rebuilding a
new type of expression tree during type analysis will be realised. The
function part is a special name that will be inserted into the identifier
table with an appropriate function value elsewhere.

@< If the type |t| can be converted from a list type... @>=
{ expr_list l=e.e.sublist; // the list of component expressions
  type_ptr comp(NULL); // the type that they should have
  if (t==vect)
  {@; @< Insert a vector conversion at |e| @>
    comp=make_type("int");
  }
  else if (t==matr)
  {@; @< Insert a matrix conversion at |e| @>
    comp=make_type("vec");
  }
  if (comp.get()!=NULL) // a substitute component type was set
  { for (; l!=NULL; l=l->next) check_type(*comp,l->e);
    break; // from enclosing |switch|
  }
}

@ The conversion functions called `\.{>vec<[int]:}' and `\.{>mat<[vec]:}' (the
names are deliberately different from any one that the user could type), which
functions will be defined later. This module and the following one will be
reused in other cases where a conversion function needs to be inserted.

@< Insert a vector conversion at |e| @>=
{ expr_list arg=make_exprlist_node(e,null_expr_list);
  e=make_application_node
       (main_hash_table->match_literal(">vec<[int]:"),arg);
}

@ The only difference here is the name of the conversion function.

@< Insert a matrix conversion at |e| @>=
{ expr_list arg=make_exprlist_node(e,null_expr_list);
  e=make_application_node
       (main_hash_table->match_literal(">mat<[vec]:"),arg);
}


@*2 Conversion to rigid values.
Now that we have a type analyser, we can write a function that converts a list
of integers into an object of type |latticetypes::Weight|, which stands for
|std::vector<int>|. The order of events will be globally the following. First
the type analysis shows that an expression has type ``row of integer'', while
it needs to be the argument of a function the requires a
|latticetypes::Weight|; it will insert a function call to a conversion routine
into the expression tree. Then the function |evaluate| above is called to
evaluate the modified expression tree. This will first evaluate the list
display into an object referred to by a |value_ptr|. Then the conversion
function |cast_intlist_to_weight| defined below will come along to convert it
into an object of the given type.

@< Includes... @>=
#include "latticetypes_fwd.h"
#include "matrix.h" // this makes |latticetypes::LatticeMatrix| a complete type

@~We shall define a conversion |cast_intlistlist_to_matrix| to matrices as
well.

@< Declarations of local functions @>=
latticetypes::Weight cast_intlist_to_weight(const value_ptr);
latticetypes::LatticeMatrix cast_intlistlist_to_matrix(const value_ptr);

@~The functions |cast_intlist_to_weight| and |cast_intlistlist_to_matrix| will
use up their arguments, and therefore take charge of destroying them. On the
other hand |row_to_weight| is just an auxiliary function and leaves the
destroying to its callers. Since exceptions may be thrown in many places, the
easiest way to ensure proper destruction is placement of the passed pointer
into an auto-pointer of type |row_ptr|.

Some people might complain that |row_to_weight| makes a copy of the vector
|result| (and its entries) at the |return| statement just before destroying
the local variable. However since there is only one |return| statement, the
code generated by a decent compiler like \.{g++} creates the vector object at
its destination in this case, and does not call the copy constructor (more
generally this is avoided if all |return| statements use the same variable).
So we shall not hesitate to use this idiom henceforth.

@< Function definition... @>=
latticetypes::Weight row_to_weight(const row_value* r)
{ latticetypes::Weight result(r->value.size());
  for(size_t i=0; i<r->value.size(); ++i)
  { int_value* n=dynamic_cast<int_value*>(r->value[i]);
    if (n==NULL)
       throw std::logic_error("Row display entry failed to be integral");
@.Row display entry failed to...@>
    result[i]=n->value;
  }
  return result;
}

latticetypes::Weight cast_intlist_to_weight(const value_ptr v)
{ row_ptr r(dynamic_cast<row_value*>(v));
  if (r.get()==NULL)
    throw std::logic_error("Row display failed to return a row");
@.Row display entry failed to...@>
  return row_to_weight(r.get()); // destroys |*r| on returning
}

@ Converting to a lattice matrix is slightly longer but not really more
complicated. We accept list entries that return either a dynamically typed
list of integers or a |Weight| value (via the ``primitive type'' |weight_value|
defined below). Short vectors will be extended with zeros if necessary. We use
the same auxiliary function |row_to_weight| as for vectors in the former case.
In either case the objects |*rr->value[i]| will get destroyed when |*rr|~is.

@< Function definition... @>=
latticetypes::LatticeMatrix cast_intlistlist_to_latmat(const value_ptr v)
{ row_ptr rr(dynamic_cast<row_value*>(v));
    // will destroy |*rr| at end
  if (rr.get()==NULL)
    throw std::logic_error("Matrix display failed to return row");
@.Matrix display failed to return row@>
 latticetypes::WeightList res_vec(rr->value.size());
  size_t depth=0; // maximal length of vectors
  for(size_t i=0; i<rr->value.size(); ++i)
  { row_value* r=dynamic_cast<row_value*>(rr->value[i]);
    if (r!=NULL) res_vec[i]=row_to_weight(r);
    else
    { weight_value* lambda=dynamic_cast<weight_value*>(rr->value[i]);
      if (lambda==NULL)
        throw std::logic_error("Matrix column failed to be a weight");
@.Matrix column failed to be a weight@>
      res_vec[i]=lambda->value;
    }
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

@ Now that we can convert to them, we can treat (converted) vectors and
integers as primitive types.

@< Other primitive types @>=
vector_type, matrix_type, @[@]

@~We use simple names, hoping that it is clear they are primitive to the
programming language offered.

@< Other primitive type names @>=
"vec", "mat", @[@]

@ Here are the corresponding types derived from |value_base|.

@< Type definitions @>=

struct weight_value : public value_base
{ latticetypes::Weight value;
  weight_value(const latticetypes::Weight& v) : value(v) @+ {}
  ~weight_value()@+ {}
  virtual void print(std::ostream& out) const;
  weight_value* clone() const @+{@; return new weight_value(*this); }
private:
  weight_value(const weight_value& v) : value(v.value) @+{}
};


struct latmat_value : public value_base
{ latticetypes::LatticeMatrix value;
  latmat_value(const latticetypes::LatticeMatrix& v) : value(v) @+ {}
  ~latmat_value()@+ {}
  virtual void print(std::ostream& out) const;
  latmat_value* clone() const @+{@; return new latmat_value(*this); }
private:
  latmat_value(const latmat_value& v) : value(v.value) @+{}
};

@ To make a small but visible difference, weights will be printed in equal
width fields one longer than the minimum necessary.
@h<sstream>
@h<iomanip>
@< Function def... @>=
void weight_value::print(std::ostream& out) const
{ size_t l=value.size(),w=0; std::vector<std::string> tmp(l);
  for (size_t i=0; i<l; ++i)
  { std::ostringstream s; s<<value[i]; tmp[i]=s.str();
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
void latmat_value::print(std::ostream& out) const
{ size_t k=value.numRows(),l=value.numColumns();
  std::vector<size_t> w(l,0);
  for (size_t i=0; i<k; ++i)
    for (size_t j=0; j<l; ++j)
    { std::ostringstream s; s<<value(i,j); size_t len=s.str().length();
      if (len>w[j]) w[j]=len;
    }
  out << std::endl << std::right;
  for (size_t i=0; i<k; ++i)
  { out << '|';
    for (size_t j=0; j<l; ++j)
      out << std::setw(w[j]+1) << value(i,j) << (j<l-1 ? ',' : ' ');
    out << '|' << std::endl;
  }
}

@*1 Other kinds of types.
At this point the order of this file starts diverging radically from the order
in which it was originally written. We shall now extend our type system in two
ways, first introducing tuple types and then function types. Originally there
were just function types that could take multiple arguments; now there are
tuple types that allow functions to have both multiple arguments and multiple
results, while function types have just one argument and one result type (but
both could be tuple types). This does mean that type checking gets a bit more
involved, for instance when a function requires a tuple argument we can expect
any expression that produces such a type, but the most common case of a tuple
display will require the component types to be distributed over its
components.

@*2 Tuple types.
Tuple types are very simple in conception: they just specify a list of
component types, and there values will consist of a sequence of such values.
Compared with array types they are much more flexible, since no relation among
the types is required; on the other hand one should not imagine performing
loops over the components of a tuple. The list of types can have any length
except~$1$: anything that would suggest a $1$-tuple (for instance a
parenthesised expression) is identified with its unique component. For the
moment we decide that there are no other identifications among tuple types,
for instance a pair whose second component is a pair is not identified with a
triplet.

@< Tags for other kinds of types @>=
tuple_type, @[@]

@ We do not need any structure for tuple types, since the type |type_list|
will do fine to represent a tuple type. However we do have to declare this a
legal variant of |type_declarator|.

@< Variants of the |union|... @>=
type_list tuple;

@ The constructor for to this variant is even simpler that the other ones,
since |new| is not called.

@< Other constructors for |type_declarator| @>=
type_declarator(type_list component_types)
  : kind(tuple_type)
{@; tuple=component_types; }

@~The copy constructor for |type_declarator| makes a deep copy, so we must
duplicate the type lists in this case. The recursion to do this is implicit,
via the copy constructor for |type_node|.

@<Other cases in the copy constructor for |type_declarator| @>=
case tuple_type:
  tuple=new type_node(*t.tuple);
break;

@ The destructor destroys the entire list, again by implicit recursion.

@<Other cases in the destructor for |type_declarator| @>=
case tuple_type:
  delete tuple;
break;

@ The function for actually building tuple types is |make_tuple_type|.

@< Declarations of local functions @>=
type_ptr make_tuple_type (type_list_ptr l);

@~As usual we must release the auto-pointer after allocating, for prudence's
sake.

@< Function def... @>=
type_ptr make_tuple_type (type_list_ptr l)
{@; type_declarator* p=new type_declarator(l.get()); l.release();
  return type_ptr(p);
}

@ Here is a straightforward structural comparison of function types for
inequality. It might seem that the increment operation of the |for| loop
below lacks caution, as the termination condition requires both pointers to be
null, but in fact the loop body quits if only one of the pointers is null.

@<Other cases for |operator!=| for |type_declarator|... @>=
case tuple_type:
  for(type_list l0=x.tuple,
                l1=y.tuple;
      l0!=NULL || l1!=NULL; l0=l0->next, l1=l1->next)
    if (l0==NULL || l1==NULL || *l0->t!=*l1->t) return true;
  return false;

@ And finally the case for printing function types.

@< Other cases for printing types @>=
case tuple_type:
  out << '(';
  for (type_list l=t.tuple; l!=NULL; l=l->next)
    out << *l->t << ( l->next!=NULL ? "," : "" );
  out << ')' ;
break;

@ We must now add the relevant cases for the variant |tuple_display|. Let us
first consider type checking. Calling |find_type| for a tuple display, we must
collect the types in a tuple display, which is a recursive activity most
easily done by a recursive function.

@< Declarations of local functions @>=
type_list_ptr find_type_list(expr_list);

@~Now we can just call and wrap up.
@< Cases for finding the type of other kinds of expressions @>=
case tuple_display:  return make_tuple_type(find_type_list(e.e.sublist));

@ The definition of this recursive function is straightforward.
@< Function def...@>=
type_list_ptr find_type_list(expr_list l)
{ if (l==NULL) return type_list_ptr(NULL);
  type_ptr t=find_type(l->e);
  type_list_ptr tail=find_type_list(l->next);
  type_node* result=new type_node(t.get(),tail.get());
  t.release(),tail.release();
  return type_list_ptr(result);
}


@ Evaluating a tuple display is much like evaluating a list display. Thanks to
dynamically typed values, we can collect the components of a tuple in a vector
without problem. In fact we could reuse the type |row_value| to hold the
components of a tuple, if it weren't for the fact that it would print with
brackets, as a list. Therefore let us trivially derive a new class from
|row_value|.

@< Type definitions @>=
struct tuple_value : public row_value
{ void print(std::ostream& out) const;
  tuple_value(const std::vector<value_ptr>& v) : row_value(v) @+{}
  tuple_value* clone() const @+{@; return new tuple_value(*this); }
private:
  tuple_value(const row_value& v);
  tuple_value& operator=(const row_value& v);
};

@ We just need to redefine the |print| method.
@< Function definitions @>=
void tuple_value::print(std::ostream& out) const
{ if (value.empty()) out << "()";
  else
  { out << '(';
    std::vector<value_ptr>::const_iterator p=value.begin();
    do {@; (*p)->print(out); ++p; out << (p==value.end() ? ')' : ','); }
    while (p!=value.end());
  }
}

@ Now we can use the same code as for list displays to evaluate it. Almost.

@< Cases for evaluating other kinds of expressions @>=
case tuple_display:
{ size_t length=0;
  for (expr_list l=e.e.sublist; l!=NULL; l=l->next) ++length;
  row_ptr v@|(new
    tuple_value(std::vector<value_ptr>(length,static_cast<value_ptr>(NULL))));
  expr_list l=e.e.sublist;
  for (size_t i=0; i<length; ++i,l=l->next)
    v->value[i]=evaluate(l->e);
  result=v.release();
}
break;

@ If the type required for a tuple display is a tuple type, we can simply test
that all component expressions have the required number and component types.
If not, we shall either call a type error with a tuple of unknown component
types.

@< Other cases for testing whether the type of |e| matches |t| @>=
case tuple_display:
  if (t.kind==tuple_type)
    @< Check that the components of |e.e.sublist| has the number and types
       specified by |t.tuple| @>
  else
  { type_list_ptr tl(NULL);
    static type_declarator unknown(undetermined_type);
    for (expr_list l=e.e.sublist; l!=NULL; l=l->next)
      tl=make_type_list(copy(&unknown),tl);
    actual=make_tuple_type(tl);
  }
break;


@ If we did find a tuple display where a tuple type was required, but the
number of components does not match, we throw a |program_error| rather than a
|type_error| to indicate the problem (this is quite coarse, it should be
improved).

@< Check that the components of |e.e.sublist| has the number and types
       specified by |t.tuple| @>=
{ type_list l=t.tuple;
  for (expr_list a=e.e.sublist; a!=NULL || l!=NULL; a=a->next,l=l->next)
    if (a==NULL || l==NULL)
      throw program_error("Too "+ std::string(a==NULL ? "few" : "many")
			 +" components in tuple");
@.Too few components in tuple@>
@.Too many components in tuple@>
    else check_type(*l->t,a->e);
}

@*2 Function types.
Before we can handle function calls, we must extend our system of types to
include function types, and store them into an identifier table so that
function calls can be type-checked. Extending our type system with function
types means we must go over a large number of small details.

@< Tags for other kinds of types @>=
function_type, @[@]

@ The variant for function types needs both an argument type and a result
type. We cannot however define the appropriate structure directly where it is
needed, since the scope of that definition would then be too limited to
perform the appropriate |new| in the constructor. Therefore we must
pre-declare the structure type. While we are doing that, we might as well
define a constructor to fill the structure (it will be used twice below). We
resist the temptation however to define a destructor: this is just an
auxiliary type, and we leave the destruction to the destructor of
|type_declarator| that incorporates (a pointer to) this type.

The |func_type| structure contains pointers to the relevant types rather, than
the |type_declarator| objects themselves, for a somewhat silly reason:
assuming we have constructed both |type_declarator| objects separately, it
would be impossible to place these objects side-by-side in a single node
without invoking the copy constructor, recreating recursively both objects. A
possible remedy would be to define an efficient |swap| method for
|type_declarator| and a dummy variant that could be produced on node creation,
after which the contents could be swapped with the formerly created descendant
objects; we have however not done so.

@< Other |struct| definitions @>=
struct func_type
{ type_declarator* arg_type,* result_type;
@)
  func_type(type_declarator* a, type_declarator* r)
   : arg_type(a), result_type(r) @+{}
};

@~Now we can use a pointer to this structure as a variant for
|type_declarator|.

@< Variants of the |union|... @>=
func_type* func;

@ The constructor corresponding to this variant is as straightforward as the
other ones.

@< Other constructors for |type_declarator| @>=
type_declarator(type_declarator* arg, type_declarator* result)
  : kind(function_type)
{@; func=new func_type(arg,result);
}

@~The copy constructor for |type_declarator| makes a deep copy. Note that this
takes some effort for this case, since we did not bother to define a (deep)
copy constructor for |func_type|.

@<Other cases in the copy constructor for |type_declarator| @>=
case function_type:
{ type_ptr a=copy(t.func->arg_type)
          ,r=copy(t.func->result_type);
  func=new func_type(a.get(),r.get()); a.release(),r.release();
}
break;

@ This code could have been reduced to |delete func| if we had defined a
suitable destructor for the structure |func_type|. But we chose not to do
that. Note that the initial |delete func->arg_types| will start by destroying
the list of argument types itself via the call to the destructor for
|type_node|.

@<Other cases in the destructor for |type_declarator| @>=
case function_type:
  delete func->arg_type; delete func->result_type;
  delete func;
break;

@ The function for actually building function types is |make_function_type|
which uses auto-pointers.

@< Function def... @>=
type_ptr make_function_type (type_ptr a, type_ptr r)
{ type_declarator* p=new type_declarator(a.get(),r.get());
  a.release(),r.release();
  return type_ptr(p);
}

@ Here is a straightforward structural comparison of function types for
inequality. It might seem that the increment operation of the |for| loop
below lacks caution, as the termination condition requires both pointers to be
null, but in fact the loop body quits if only one of the pointers is null.

@<Other cases for |operator!=| for |type_declarator|... @>=
case function_type:
  return  *x.func->arg_type!=*x.func->arg_type
     ||   *x.func->result_type!=*x.func->result_type;

@ And finally the case for printing function types. Here we suppress
additional parentheses in case the argument or result type is a tuple type.

@< Other cases for printing types @>=
case function_type:
  out << '(';
  if (t.func->arg_type->kind!=tuple_type)
     out << *t.func->arg_type << "->";
  else
    for (type_list l=t.func->arg_type->tuple;
         l!=NULL; l=l->next)
    out << *l->t << ( l->next!=NULL ? "," : "->" );
  if (t.func->result_type->kind!=tuple_type)
     out << *t.func->result_type << ")";
  else
    for (type_list l=t.func->result_type->tuple;
         l!=NULL; l=l->next)
    out << *l->t << ( l->next!=NULL ? "," : ")" );
break;

@*1 Function calls.
We shall now extend the evaluator as described until now to allow function
calls. This will turn out to involve several new notions as well.

@*2 Identifier tables.
As we said above, we need an identifier table to record the types of known
functions (and later other types of identifiers). For the moment we use a
single flat table; there will doubtlessly be need for handling scopes later.
We shall actually store both a type and a value in the table.

@< Type definitions @>=

struct id_data
{ value_ptr value; @+ type_declarator* type;
  id_data(value_ptr v,type_declarator* t) : value(v),type(t)@+ {}
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
  void add(Hash_table::id_type id, value_ptr v, type_ptr t); // insertion
  type_declarator* type_of(Hash_table::id_type id) const; // lookup
  value_ptr value_of(Hash_table::id_type id) const; // lookup
@)
  size_t size() const @+{@; return table.size(); }
  void print(std::ostream&) const;
};

@ Since we have stored pointers, the destructor must explicitly delete them,
even though we have no immediate plans for destroying identifier tables.

@< Function def... @>=
Id_table::~Id_table()
{ for (map_type::const_iterator p=table.begin(); p!=table.end(); ++p)
  @/{@; delete p->second.value; delete p->second.type; }
}

@ For insertion we must distinguish the case that the key is already present,
since the |insert| method for maps prefers not to overwrite the old value if
there is one, but rather to make no change. It does return in the case  both a
pointer (iterator) to the (key,data) pair that obstructed the insertion and a
boolean failure status, so that we can easily overwrite ``manually'' the old
value.

@< Function... @>=
void Id_table::add(Hash_table::id_type id, value_ptr v, type_ptr t)
{ id_data data(v,t.get()); std::auto_ptr<value_base> safe(v);
  std::pair<map_type::iterator,bool> trial
     =table.insert(std::make_pair(id,data));
  safe.release(),t.release(); // no more exception protection needed once here
  if (!trial.second) // then key was present; destroy and replace its data
@/{@; delete trial.first->second.value; delete trial.first->second.type;
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
value_ptr Id_table::value_of(Hash_table::id_type id) const
{@; map_type::const_iterator p=table.find(id);
  return p==table.end() ? NULL : p->second.value;
}

@ We provide a |print| member the shows the contents of the entire table.
@< Function... @>=

void Id_table::print(std::ostream& out) const
{ for (map_type::const_iterator p=table.begin(); p!=table.end(); ++p)
    out << main_hash_table->name_of(p->first) << ": " @|
        << *p->second.type << ": " << *p->second.value << std::endl;
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
The function type stored in the table prescribes the required type for the
argument (which could be a tuple type representing multiple arguments), and
also gives the final type the expression will have. Note that since we do
not own the types coming from the table, we must copy them if we want to
export them, either as type of the call or in the error value thrown in case
of a type mismatch.

@< Cases for finding the type of other kinds of expressions @>=
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
  check_type(*f_type->func->arg_type,e.e.call_variant->arg);
  return copy(f_type->func->result_type);
}

@ When a function call appears where a fixed type is expected, we first test
that the function returns this type, and then go on to check the types of the
arguments.

@< Other cases for testing whether the type of |e| matches |t| @>=
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
  if (*f_type->func->result_type!=t)
    actual=copy(f_type->func->result_type); // signal an error
  else check_type(*f_type->func->arg_type,e.e.call_variant->arg);
}
break;

@*2 Evaluating function calls.
The evaluation of a function call may either be demanded explicitly by the
user, or by the insertion performed by the code above. In order to pass a
number of arguments unknown at compile time, we must pass them indirectly, and
the most natural way to do that is via a stack. In fact we could pass a single
tuple-value, but using a stack we can avoid constructing the tuple. For the
moment however we do construct the tuple an then unpack it to the stack.

All usable built-in functions will be provided with a small wrapper function
that takes and unpacks the values from the stack; the last parameter is popped
from the stack first. There should be no error in unpacking, since types have
already been tested. After unpacking, the wrapper function then calls the
built-in function. For the stack we shall use a vector valued variable.
@< Declarations of global variables @>=
extern std::vector<value_ptr> execution_stack;

@~
@< Global variable definitions @>=
std::vector<value_ptr> execution_stack;

@ The stack owns the values it contains, but there is no reason to wrap it
into a class with a destructor, since we never intend to destroy the stack: if
our program exits either peacefully or by an uncaught exception we don't care
about some values that are not destroyed. We must remember however to |delete|
the values whenever we empty the stack after catching a runtime error.
Therefore we provide a function to clear the stack.

@< Declarations of exported functions @>=
void clear_execution_stack ();

@~We monitor disappearing values when clearing the stack. Since currently the
stack only holds values very briefly during entry and exit of function calls,
it is virtually impossible to actually make this function produce output from
user expressions, but some wrapper functions simulate the evaluation of nested
expressions, and in case of errors the values already computed but not
involved in the error may show up.

@< Function definitions @>=
void clear_execution_stack ()
{ if (!execution_stack.empty())
  { std::cerr << "Discarding from execution stack:" << std::endl;
    do
    { value_ptr v=execution_stack.back();
      std::cerr << *v << std::endl;
      delete v; execution_stack.pop_back();
    }
    while (!execution_stack.empty());
  }
}

@ Here are two functions to facilitate manipulating the stack. Being inline,
we define them right away.

@< Declarations of exported functions @>=
inline void push_value(value_ptr p)
@+{@; execution_stack.push_back(p); }
@)
inline value_ptr pop_arg()
{@; value_ptr arg=execution_stack.back(); execution_stack.pop_back();
    return arg;
}

@ The result of wrapper functions will be pushed on the stack as a
|value_ptr|, so a wrapper function has neither arguments nor a result type.
Thus variables that refer to a wrapper function have the type
|wrapper_function| defined below. We shall need to bind values of this type to
identifiers representing built-in functions, so we derive an associated
``primitive type'' from |value_base|.

@< Type definitions @>=
typedef void (* wrapper_function)();

struct builtin_value : public value_base
{ wrapper_function value;
  std::string print_name;
  builtin_value(wrapper_function v,const char* n)
  : value(v), print_name(n) @+ {}
  ~builtin_value()@+ {} // don't try to destroy the function pointed to!
  virtual void print(std::ostream& out) const
  @+{@; out << ':' << print_name << ':'; }
  builtin_value* clone() const @+{@; return new builtin_value(*this); }
private:
  builtin_value(const builtin_value& v)
  : value(v.value), print_name(v.print_name) @+{} // copy constructor

};

@ Finally we can say what to do when a function call is requested. Recall that
for the moment the function expression is always an identifier, so the
associated (wrapper function) value can be found in the |global_id_table|.

@< Cases for evaluating other kinds of expressions @>=
case function_call:
{ push_value(evaluate(e.e.call_variant->arg));
    // evaluate and push argument
  std::string name=main_hash_table->name_of(e.e.call_variant->fun);
    // for error messages
  value_ptr f_val=global_id_table->value_of(e.e.call_variant->fun);
  if (f_val==NULL) throw
    std::logic_error("Built-in function absent: "+name);
@.Built-in function absent@>
  builtin_value* b=dynamic_cast<builtin_value*>(f_val);
  if (b==NULL) throw
    std::logic_error("Built-in not a function: "+name);
@.Built-in not a function@>
  b->value(); // call the wrapper function, leaving result on the stack
  result=pop_arg(); // get the result to return it from |evaluate|
}
break;

@*1 Wrapper functions.
Wow, we should be ready to roll! Except that we have not defined any wrapper
functions yet, and (therefore) have nothing in the |global_id_table|. Let us
start with some preparations for the former. Wrapper functions will routinely
have to do some unpacking  of values, which will involve dynamically casting
the |value_ptr| to the type they are known to have because we passed the type
checker; should the cast fail we shall throw a |std::logic_error|. To avoid
having |throw| statements all over the place, we define auxiliary
functions to do the unpacking, one for every type.

@< Declarations of exported functions @>=
int_value* get_int() throw(std::logic_error);
string_value* get_string() throw(std::logic_error);
bool_value* get_bool() throw(std::logic_error);
weight_value* get_vec() throw(std::logic_error);
latmat_value* get_mat() throw(std::logic_error);
row_value* get_row() throw(std::logic_error);
tuple_value* get_tuple() throw(std::logic_error);

@ The implementation of these functions is so similar that we almost defined a
template for them.

@< Definition of other wrapper functions @>=
int_value* get_int() throw(std::logic_error)
{ value_ptr p=pop_arg();
  int_value* result=dynamic_cast<int_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not an integer"); }
  return result;
}

string_value* get_string() throw(std::logic_error)
{ value_ptr p=pop_arg();
  string_value* result=dynamic_cast<string_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not a string"); }
  return result;
}

bool_value* get_bool() throw(std::logic_error)
{ value_ptr p=pop_arg();
  bool_value* result=dynamic_cast<bool_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not a boolean"); }
  return result;
}


weight_value* get_vec() throw(std::logic_error)
{ value_ptr p=pop_arg();
  weight_value* result=dynamic_cast<weight_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not a vector"); }
  return result;
}

latmat_value* get_mat() throw(std::logic_error)
{ value_ptr p=pop_arg();
  latmat_value* result=dynamic_cast<latmat_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not a matrix"); }
  return result;
}

row_value* get_row() throw(std::logic_error)
{ value_ptr p=pop_arg();
  row_value* result=dynamic_cast<row_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not a row"); }
  return result;
}

tuple_value* get_tuple() throw(std::logic_error)
{ value_ptr p=pop_arg();
  tuple_value* result=dynamic_cast<tuple_value*>(p);
  if (result==NULL)
    {@; delete p; throw std::logic_error("Argument is not a tuple"); }
  return result;
}

@.Argument is not ...@>

@ The function |push_tuple_components| will be called by wrapper functions
that need the tuple components on the stack; the function call |wrap_tuple(n)|
inversely builds a tuple from $n$ components on the stack.

@< Declarations of exported functions @>=
void push_tuple_components();
void wrap_tuple(size_t n);

@ These functions use the same convention for stack order: the last tuple
component is on top of the stack.

@< Function definitions @>=
void push_tuple_components()
{ std::auto_ptr<tuple_value> tuple(get_tuple());
  for (size_t i=0; i<tuple->length(); ++i)
  { push_value(tuple->value[i]); // push component
    tuple->value[i]=NULL; // remove component so it remains unshared
  }
}
@)
void wrap_tuple(size_t n)
{ tuple_value* result=new tuple_value(std::vector<value_ptr>(n));
  while (n-->0) // not |--n>=0| since |n| is unsigned!
    result->value[n]=pop_arg();
  push_value(result);
}

@ Here are our first wrapper functions. The function |id_wrapper| is a trivial
function, but it will be put into the |main_id_table| under different names
and signatures, each forcing a particular argument type.

@< Function definitions @>=
void intlist_to_weight_wrapper ()
{@; push_value(new weight_value(cast_intlist_to_weight(pop_arg())));
}

void intlistlist_to_latmat_wrapper ()
{@; push_value(new latmat_value(cast_intlistlist_to_latmat(pop_arg())));
}

void id_wrapper () @+{} // nothing to do, value stays on the stack

@< Definition of other wrapper functions @>@;

@ All that remains to do for the moment is to install these functions into the
|global_id_table|.

@< Declarations of exported functions @>=
void initialise_evaluator();
@)
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
Let us recapitulate what will happen. The parser will read what the user
types, and returns an |expr| value. If the user typed applications of |"vec"|
or |"mat"|, then the necessity to transform their argument will be detected
during type checking, and appropriate calls to conversion functions will be
inserted (this is why the functions for |"vec"| and |"mat"| have nothing to do
themselves). Then everything will be executed by calling |evaluate|, and the
main program will print the result. Ahh\dots! Nobody ever bothered to call the
type checker yet. We had better do that right now. Since rewriting is
involved, we declare the function with a non-constant reference parameter,
although this does not really make a difference (since the rewriting takes
place below the topmost level of the expression).

@< Declarations of exported functions @>=
type_ptr analyse_types(expr& e) throw(std::bad_alloc,std::runtime_error);

@~In fact it suffices to call |find_type| on the expression given, and return
the result. However this function will give us an occasion to catch any
thrown |type_error| and |program_error| exceptions, something we did not want
to do inside the recursive function |find_type|.

@< Function definitions @>=
type_ptr analyse_types(expr& e) throw(std::bad_alloc,std::runtime_error)
{ try {@; return find_type(e); }
  catch (type_error& e)
  { std::cerr << e.what() << std::endl <<
    "Subexpression " << e.offender << @| " has wrong type: found "
         << *e.actual << " while " << *e.required << " was needed.\n";
@.Subexpression has wrong type@>
  }
  catch (program_error& err)
  { std::cerr << err.what() <<
          " in expression '" << e << "'\n";
  }
  throw std::runtime_error("Type check failed");
@.Type check failed@>
}

@*1 Built-in functions.
We can hardly believe it, but that worked almost straight out of the box! (The
only problems were some glitches in the printing of vectors and matrices. But
it should be said that that was before we put any of the cleanup code in
place; doing the did provoke quite a few core dumps.) Encouraged by that, we
shall now introduce some real built-in functions, starting with integer
arithmetic. Arithmetic operators are implemented by wrapper functions with two
integer arguments (or one in the case of unary minus). Note that the values
are pulled from the stack in reverse order, which is important for the
non-commutative operations like `|-|', `|/|' and~`|%|'. Since values are not
shared, we reuse one the first value object and destroy the second.

@< Definition of other wrapper functions @>=
void plus_wrapper ()
{ push_tuple_components();
  std::auto_ptr<int_value>j(get_int()); std::auto_ptr<int_value>i(get_int());
  i->value+=j->value;
  push_value(i.release());
}
void minus_wrapper ()
{ push_tuple_components();
  std::auto_ptr<int_value>j(get_int()); std::auto_ptr<int_value>i(get_int());
  i->value-=j->value;
  push_value(i.release());
}

void times_wrapper ()
{ push_tuple_components();
  std::auto_ptr<int_value>j(get_int()); std::auto_ptr<int_value>i(get_int());
  i->value*=j->value;
  push_value(i.release());
}

void divide_wrapper ()
{ push_tuple_components();
  std::auto_ptr<int_value>j(get_int()); std::auto_ptr<int_value>i(get_int());
  if (j->value==0) throw std::runtime_error("Division by zero");
  i->value/=j->value;
  push_value(i.release());
}

void modulo_wrapper ()
{ push_tuple_components();
  std::auto_ptr<int_value> j(get_int()); std::auto_ptr<int_value> i(get_int());
  if (j->value==0) throw std::runtime_error("Modulo zero");
  i->value%=j->value;
  push_value(i.release());
}

void unary_minus_wrapper ()
{@; int_value* i=get_int(); i->value =-i->value; push_value(i); }

void divmod_wrapper ()
{ push_tuple_components();
  std::auto_ptr<int_value> j(get_int()); std::auto_ptr<int_value> i(get_int());
  if (j->value==0) throw std::runtime_error("DivMod by zero");
  int mod=i->value%j->value;
  i->value/=j->value; j->value=mod;
@/push_value(i.release()); push_value(j.release()); wrap_tuple(2);
}

@ We install the arithmetic operators. Their names correspond to the ones used
in the parser definition file \.{parser.y}.

@< Installation of other built-in functions @>=
install_function(plus_wrapper,"+","(int,int->int)");
install_function(minus_wrapper,"-","(int,int->int)");
install_function(times_wrapper,"*","(int,int->int)");
install_function(divide_wrapper,"/","(int,int->int)");
install_function(modulo_wrapper,"%","(int,int->int)");
install_function(unary_minus_wrapper,"-u","(int->int)");
install_function(divmod_wrapper,"/%","(int,int->int,int)");


@ We now define a few functions, to really exercise something, first of all
the identity matrix and matrix transposition.

@< Declarations of exported functions @>=
void id_mat_wrapper ();
void transpose_mat_wrapper ();

@ Since in general built-in functions may throw exceptions (even |transpose|!)
we hold the pointers to the values created to hold their results in
auto-pointers. The reason that in |id_mat_wrapper| we create a |latmat_value|
around an empty |LatticeMatrix| rather than build a filled matrix object
first, is that the constructor for |latmat_value| would than have to copy that
matrix object (which implies copying its contents). In |transpose_mat_wrapper|
we can in fact return the same |latmat_value| that held the argument, but this
does not absolve us from using an auto-pointer while |transpose| is active. We
also define |diagonal_wrapper|, a slight generalisation of |id_mat_wrapper|
that produces a diagonal matrix from a vector.

@< Definition of other wrapper functions @>=
void id_mat_wrapper ()
{ std::auto_ptr<int_value> i(get_int());
  std::auto_ptr<latmat_value> m
     (new latmat_value(latticetypes::LatticeMatrix()));
  identityMatrix(m->value,std::abs(i->value)); push_value(m.release());
}

void transpose_mat_wrapper ()
{@; std::auto_ptr<latmat_value>m(get_mat());
  m->value.transpose(); push_value(m.release());
}

void diagonal_wrapper ()
{ std::auto_ptr<weight_value> d(get_vec());
  size_t n=d->value.size();
  std::auto_ptr<latmat_value> m
     (new latmat_value(latticetypes::LatticeMatrix(n,n,0)));
  for (size_t i=0; i<n; ++i) m->value(i,i)=d->value[i];
  push_value(m.release());
}

@ Now the product of a matrix and a vector or matrix. The first of these was
in fact our first function with more than one argument (arithmetic on integer
constants was done inside the parser at that time). We make them callable from
other compilation units.

@< Declarations of exported functions @>=
void mv_prod_wrapper ();
void mm_prod_wrapper ();

@ In |mv_prod_wrapper| we use the matrix method |apply| which requires its
output vector to be already of the proper size, but nevertheless it copies its
constant second argument, for in case it should coincide with the first; it
may therefore throw an exception, and we use an auto-pointer for the
result~|w|. In |mm_prod_wrapper|, the method |operator*=| does its own
resizing (and may also throw an exception).

@< Definition of other wrapper functions @>=
void mv_prod_wrapper ()
{ push_tuple_components();
  std::auto_ptr<weight_value> v(get_vec());
  std::auto_ptr<latmat_value> m(get_mat());
  if (m->value.numColumns()!=v->value.size())
  { std::ostringstream s;
    s<< "Size mismatch " << m->value.numColumns() << ":" << v->value.size()
     << " in mv_prod";
    throw std::runtime_error(s.str());
  }
  std::auto_ptr<weight_value> w@|
    (new weight_value(latticetypes::Weight(m->value.numRows())));
  m->value.apply(w->value,v->value);
  push_value(w.release());
}
@)
void mm_prod_wrapper ()
{ push_tuple_components();
  std::auto_ptr<latmat_value> r(get_mat()); // right factor
  std::auto_ptr<latmat_value> l(get_mat()); // left factor
  if (l->value.numColumns()!=r->value.numRows())
  { std::ostringstream s;
    s<< "Size mismatch " << l->value.numColumns() << ":" << r->value.numRows()
    @| << " in mm_prod";
    throw std::runtime_error(s.str());
  }
  l->value*=r->value;
  push_value(l.release());
}

@ Here is finally the Smith normal form algorithm. We provide both the
invariant factors and the rewritten basis on which the normal for is assumed,
as separate functions, and to illustrate the possibilities of tuples, the two
combined into a single function.

@h "smithnormal.h"
@h "smithnormal_def.h"

@< Definition of other wrapper functions @>=
void invfact_wrapper ()
{ std::auto_ptr<latmat_value> m(get_mat());
  size_t nr=m->value.numRows();
  latticetypes::WeightList b; @+ matrix::initBasis(b,nr);
  std::auto_ptr<weight_value> inv_factors
     @| (new weight_value(latticetypes::Weight(0)));
  smithnormal::smithNormal(inv_factors->value,b.begin(),m->value);
  push_value(inv_factors.release());
}
void Smith_basis_wrapper ()
{ std::auto_ptr<latmat_value> m(get_mat());
  size_t nr=m->value.numRows();
  latticetypes::WeightList b; @+ matrix::initBasis(b,nr);
  latticetypes::Weight inv_factors(0);
  smithnormal::smithNormal(inv_factors,b.begin(),m->value);
  latticetypes::LatticeMatrix new_basis(b); // convert basis into matrix
  m->value.swap(new_basis);
  push_value(m.release());
}

void Smith_wrapper ()
{ std::auto_ptr<latmat_value> m(get_mat());
  size_t nr=m->value.numRows();
  latticetypes::WeightList b; @+ matrix::initBasis(b,nr);
  std::auto_ptr<weight_value> inv_factors
     @| (new weight_value(latticetypes::Weight(0)));
  smithnormal::smithNormal(inv_factors->value,b.begin(),m->value);
  latticetypes::LatticeMatrix new_basis(b); // convert basis into matrix
  m->value.swap(new_basis);
@/push_value(m.release()); push_value(inv_factors.release()); wrap_tuple(2);
}

@ Here is one more wrapper function that uses the Smith normal form algorithm,
but behind the scenes, namely to invert a matrix. Since this cannot be done in
general over the integers, we return an integral matrix and a common
denominator to be applied to all coefficients.
@< Definition of other wrapper functions @>=
void invert_wrapper ()
{ std::auto_ptr<latmat_value> m(get_mat());
  if (m->value.numRows()!=m->value.numColumns())
  { std::ostringstream s;
    s<< "Cannot invert a " @|
     << m->value.numRows() << "x" << m->value.numColumns() << " matrix";
    throw std::runtime_error(s.str());
  }
  std::auto_ptr<int_value> denom(new int_value(0));
  m->value.invert(denom->value);
@/push_value(m.release()); push_value(denom.release()); wrap_tuple(2);
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
After function applications were working, we added identifiers. The definition
of identifiers is currently only at the global level, handled separately from
expression evaluation (it does call the evaluator, but is not a case within
the evaluator). So as far as this compilation unit is concerned, the only
point to consider is evaluating applied identifiers. We just pick the value
from the identifier table, not forgetting to duplicate it, since values are
destroyed after being used in evaluation. Actually implementing this was more
work than it would seem, since we had to introduce the |clone| method to do
it.

@< Cases for evaluating other kinds of expressions @>=
case applied_identifier:
{ value_ptr p=global_id_table->value_of(e.e.identifier_variant);
  if (p==NULL) throw std::logic_error
  @|   ("Identifier without value:"
	+main_hash_table->name_of(e.e.identifier_variant));
@.Identifier without value@>
  result=p->clone();
}
break;

@ For type checking matters are quite the same, just copy a type plucked
from the table.

@< Cases for finding the type of other kinds of expressions @>=
case applied_identifier:
{ type_declarator* t=global_id_table->type_of(e.e.identifier_variant);
  if (t==NULL) throw program_error
  @|   ("Undefined identifier "
	+main_hash_table->name_of(e.e.identifier_variant));
@.Undefined identifier@>
  return copy(t);
}

@ For an applied identifier there is not much difference in type checking for
the case where a definite type is required: the type found in the table must
now equal the specified one. Only we must cater for the fact that the types
only match after conversion. In this case we cannot push down the conversion
into the expression as we did for list displays so we must handle simple and
composite conversions.

@< Other cases for testing whether the type of |e| matches |t| @>=
case applied_identifier:
{ type_declarator* it=global_id_table->type_of(e.e.identifier_variant);
  if (it==NULL) throw program_error
  @|   ("Undefined identifier "
	+main_hash_table->name_of(e.e.identifier_variant));
@.Undefined identifier@>
  if (*it!=t)
  { static type_declarator row_int=*make_type("[int]").release();
    static type_declarator row_row_int=*make_type("[[int]]").release();
    static type_declarator row_vec=*make_type("[vec]").release();
    if (t==vect && *it==row_int)
      @< Insert a vector conversion at |e| @>
    else if (t==matr && (*it==row_row_int || *it==row_vec))
      @< Insert a matrix conversion at |e| @>
    else actual=copy(it); // if no conversion works, signal an error
  }
}
break;

@* Operations other than evaluation of expressions.
This file also defines some operations that can be invoked by the parser that
do not consist of evaluation an expression. The declarations of these
functions are given in \.{parsetree.h} so that the parser can see them.

@*1 Making global definitions.
For the moment, applied identifiers can only get their value through the
function |global_set_identifier| that was declared in \.{parsetree.h}. We
define it here since it uses the services of the evaluator. It will be called
by the parser for global identifier definitions (simple or multiple);
therefore it has \Cee-linkage.

Recall that the parser guarantees that |ids| represents a list of identifiers.
What has to be done here is straightforward. We type-check the expression~|e|,
storing its result; then we test if, in case the left hand side has more than
one identifier, it is an appropriate tuple type. If there is no error, then we
evaluate the expression, and if everything has gone well we store the
(type,value) pair(s) into the global identifier table. Note that although this
function may itself throw a |runtime_error|, the catch clause is also there to
catch errors produced during the call to |evaluate|.

@< Function definitions @>=
extern "C"
void global_set_identifier(expr_list ids, expr e)
{ using namespace atlas::interpreter; using namespace std;
  try
  { type_ptr t=analyse_types(e);
    if (ids->next!=NULL)
      @< Check that identifiers are distinct and that |t| is an appropriate
         tuple type; if not, |throw| a |runtime_error| @>
    value_ptr v=evaluate(e);
    if (ids->next==NULL)
    { std::cout << "Identifier " << ids->e << ": " << *t << std::endl;
      global_id_table->add(ids->e.e.identifier_variant,v,t); // releases |t|
    }
    else
    { auto_ptr<tuple_value> tv(dynamic_cast<tuple_value*>(v));
      if (tv.get()==NULL) throw logic_error("Non-tuple value assigned");
@.Non-tuple value assigned@>
      std::cout << "Identifiers ";
      size_t i=0; type_list tl=t->tuple;
      for (expr_list l=ids; l!=NULL; l=l->next,++i,tl=tl->next)
      { std::cout << l->e << ": " << *tl->t << ( l->next!=NULL ? ", " : ".\n");
        global_id_table->
          add(l->e.e.identifier_variant,tv->value[i],copy(tl->t));
        tv->value[i]=NULL; // ensure value in table is unshared
      }
    }
  }
  catch (runtime_error& err)
  { cerr << err.what() << ", identifier" << (ids->next!=NULL ? "s " :" ");
    for (expr_list l=ids; l!=NULL; l=l->next)
      cerr << main_hash_table->name_of(l->e.e.identifier_variant)
           << (l->next!=NULL?",":"");
    cerr << " not defined.\n";
    destroy_exprlist(ids); destroy_expr(e);
    clear_execution_stack();
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

@*1 Printing type information.
It is useful to print type information, either for a single expression or for
all identifiers in the table.

@< Function definitions @>=
extern "C"
void type_of_expr(expr e)
{ try { std::cout << "type: " << *analyse_types(e) << std::endl; }
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
{ std::cout << *global_id_table;
}


@* Specifying types by strings.
The task of converting a properly formatted string into a type is one of
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
{ try {@; return scan_type(s); }   // provide an lvalue
  catch (std::logic_error e)
  { std::cerr << e.what() << "; text remaining: " << s << std::endl;;
    return make_prim_type(string_type);
  }
}

type_ptr scan_type(const char*& s)
{ if (*s=='[')
    @< Scan and |return| a row type, or |throw| a |logic_error| @>
  if (*s=='(')
    @< Scan and |return| a tuple or function type,
       or |throw| a |logic_error| @>
  @< Scan and |return| a primitive type, or |throw| a |logic_error| @>
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
variable to be able to share the code for constructing the tuple type.

The only complication is that single parenthesised types, and single argument
or return types should not be converted into tuple types with one component,
but just into the constituent type. This is done by assigning the extracted or
constructed types to the |type_ptr| variables |a| and |r|, that were
initialised to~|NULL|; in the case of single types that type is then unlinked
to prevent being destroyed upon the destruction of the type node (the latter
part would be hard to do if we wanted to initialise |a| and~|r| directly to
the right value). Note that this is one of the few places where we really use
the ownership-tracking semantics of auto-pointers, in the sense that their
destruction behaviour at a certain point is variable: the node pointed to by
|l0| and |l1| will only be deleted if the auto-pointer was not passed on in a
call to |make_tuple_type|.

@< Scan and |return| a tuple or function type, or |throw| a |logic_error| @>=
{ type_list_ptr l0=scan_type_list(++s), l1(NULL);
  bool is_tuple=*s==')';
  if (*s=='-' && *++s=='>') l1=scan_type_list(++s);
  if (*s++!=')') throw std::logic_error("Missing ')' in type");
  type_ptr a(NULL);
  if (l0.get()!=NULL && l0->next==NULL) {@; a=type_ptr(l0->t); l0->t=NULL; }
  else a=make_tuple_type(l0);
  if (is_tuple) return a;
  type_ptr r(NULL);
  if (l1.get()!=NULL && l1->next==NULL) {@; r=type_ptr(l1->t); l1->t=NULL; }
  else r=make_tuple_type(l1);
  return make_function_type(a,r);
}

@ A comma-separated list of types is handled by a straightforward recursion.
@< Function definitions @>=
@)
type_list_ptr scan_type_list(const char*& s)
{ type_ptr head=scan_type(s);
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
      return make_prim_type(static_cast<primitive>(i));
    }
  }
  throw std::logic_error("Type unrecognised");
}

@* Built-in types defined elsewhere.
In order for this compilation unit to function properly, it must know of the
existence and names for other built-in types. We could scoop up these names
using clever \&{\#include} directives, but that is not really worth the hassle
(when adding such types you have to recompile this unit anyway; it is not so
much worse to actually extend the lines below as well).

@< Other primitive types @>=
complex_lie_type_type , root_datum_type, complexgroup_type, @[@]

@~@< Other primitive type names @>=
"LieType","RootDatum", "ComplexGroup", @[@]

@* Index.

% Local IspellDict: default
