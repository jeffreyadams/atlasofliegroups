@* Building the parse tree.
This is the program unit \.{parsetree} which produces the implementation file
\.{parsetree.cpp} and the header file \.{parsetree.h}; the latter is read in
by \&{\#include} from the \Cee-file generated from \.{parser.y}, so it should
contain declarations in \Cee-format. Some other functions defined here are not
used in the parser, but can be used by other modules written in~\Cpp; their
declaration is skipped when the header file is read in from a
\Cee-file.

@( parsetree.h @>=
#ifndef PARSETREE_H
#define PARSETREE_H
#if __cplusplus
#include <iostream>
extern "C" {
#endif
@< Declarations in \Cee-style for the parser @>@;
#if __cplusplus
}@;

namespace atlas
{@; namespace interpreter
  {@;
@< Declaration of \Cpp\ functions @>@;
  }@;
}@;

#endif
#endif

@ The main file \.{parsetree.cpp} contains the implementations of the
functions that are needed to build the parse tree. Since these are to be
called from the parser, we declare them all to be callable from~\Cee.

@h "parsetree.h"
@c
extern "C" {
@< Definitions of functions in \Cee-style for the parser @>@;
}@;
namespace atlas
{@; namespace interpreter
  {@;
@< Definitions of \Cpp\ functions @>@;
  }@;
}@;

@ For a large part the declarations for the parser consist of the recursive
definition of the type |expr|.

@< Declarations in \Cee-style for the parser @>=
@< Typedefs that are required in |union expru@;| @>@;
union expru {@; @< Variants of |union expru @;| @>@; };

typedef enum @+
{ @< Enumeration tags for |expr_kind| @> @;@; } expr_kind;
typedef struct {@; union expru e; expr_kind kind; } expr;

@< Structure and typedef declarations for types built upon |expr| @>@;

@< Declaration of functions in \Cee-style for the parser @>@;

@ While we are defining functions to parse expressions, we shall also define
a function to print the expressions once parsed; this provides a useful test
to see if what we have read in corresponds to what was typed, and this
functionality will also be used in producing error messages.

@< Declaration of \Cpp\ functions @>=
std::ostream& operator<< (std::ostream& out, expr e);

@~The definitions of this instance of the operator~`|<<|' are distributed
among the different variants of the |union expru@;| that we shall define.

@< Definitions of \Cpp\ functions @>=
std::ostream& operator<< (std::ostream& out, expr e)
{@; switch (e.kind)
  {@; @< Cases for printing an expression |e| @>
  }
  return out;
}


@ In parallel, we also define a function |destroy_expr| to clean up the memory
occupied by an expression. It is classified as a parsing function since it is
called amongst others by the parser when popping off tokens at syntax errors.

@< Declaration of functions in \Cee-style for the parser @>=
void destroy_expr(expr e);

@~The definition of |destroy_expr| is also distributed among the different
variants of the |union expru@;|.

@< Definitions of functions in \Cee-style for the parser @>=
void destroy_expr(expr e)
{@; switch (e.kind)
  {@; @< Cases for destroying an expression |e| @>
  }
}


@*1 Atomic expressions.
The simplest expressions are atomic constants, which we shall call
denotations. There are recognised by the scanner, and either the scanner or
the parser will build an appropriate node for them, which just stores the
constant value denoted. We need no structures here, since the value itself
will fit comfortably inside the |union expru@;|. The fact that strings are
stored as a character pointer, which should be produced using |new[]| by the
caller of the functions described here, is a consequence of the restriction
that only types representable in \Cee\ can be used.

@< Variants... @>=

int int_denotation_variant;
char* str_denotation_variant;

@~But each type of denotation has a tag identifying it.

@< Enumeration tags for |expr_kind| @>=
integer_denotation, string_denotation, boolean_denotation, @[@]

@ To print an integer denotation we just print its |v| field; for string
denotations we do the same but enclosed in quotes, while for Boolean
denotations we reproduce the keyword that gives the denotation.

@< Cases for printing... @>=
case integer_denotation: out << e.e.int_denotation_variant; break;
case string_denotation:
  out << '"' << e.e.str_denotation_variant << '"'; break;
case boolean_denotation:
  out << (e.e.int_denotation_variant ? "true" : "false"); break;

@~When a string denotation is destroyed, we free the string that was created
by the lexical scanner.

@< Cases for destroying... @>=
case integer_denotation: case boolean_denotation: break;
case string_denotation: delete[] e.e.str_denotation_variant; break;

@ To build the node for denotations, we provide the functions below.

@< Declaration of functions in \Cee-style for the parser @>=
expr make_int_denotation (int val);
expr make_string_denotation(char* val);
expr make_bool_denotation(int val);

@~The definition of these functions is quite trivial, as will be typical for
node-building functions.

@h <cstdio>

@< Definitions of functions in \Cee... @>=
expr make_int_denotation (int val)
{ expr result; result.kind=integer_denotation;
@/result.e.int_denotation_variant=val; return result;
}
expr make_string_denotation(char* val)
{ expr result; result.kind=string_denotation;
@/result.e.str_denotation_variant=val; return result;
}
expr make_bool_denotation(int val)
{ expr result; result.kind=boolean_denotation;
@/result.e.int_denotation_variant=val; return result;
}

@ Another atomic expression is an applied identifier. They use the type
|id_type| that should ideally be |Hash_table::id_type|, but since we are
restricted to using \Cee-syntax for defining |expr|, we see no better way than
to redefine a global typedef.

@< Typedefs... @>=
typedef short id_type;

@ For identifiers we just store their code.
@< Variants... @>=
id_type identifier_variant;

@~Their tag is |applied_identifier|.

@< Enumeration tags for |expr_kind| @>=
applied_identifier, @[@]

@ To print an applied identifier, we look it up in the main hash table.

@< Cases for printing... @>=
case applied_identifier:
  out << main_hash_table->name_of(e.e.identifier_variant);
break;

@~For destroying an applied identifier, there is nothing to do.

@< Cases for destroying... @>=
case applied_identifier: break;

@ To build the node for denotations we provide the function
|make_applied_identifier|.

@< Declaration of functions in \Cee-style for the parser @>=
expr make_applied_identifier (id_type id);

@~The definition of |make_applied_identifier| is entirely trivial.

@< Definitions of functions in \Cee... @>=
expr make_applied_identifier (id_type id)
{@; expr result; result.kind=applied_identifier;
  result.e.identifier_variant=id; return result;
}

@*1 Expression lists, and list and tuple displays.
A first recursive type of expression is the expression list, which will be
used for various purposes.

@< Typedefs... @>=
typedef struct exprlist_node* expr_list;

@~The type is implemented as a simply linked list.
@< Structure and typedef declarations for types built upon |expr| @>=
struct exprlist_node {@; expr e; expr_list next; };

@ Any syntactic category whose parsing value is a list of expressions will use
the variant |sublist|.

@< Variants... @>=
expr_list sublist;

@ An evident category with an |expr_list| as parsing value is a bracketed list
for designating a vector; we shall call such an expression a list display.

@< Enumeration tags for |expr_kind| @>= list_display, @[@]

@ To print a list display, we call the operator~`|<<|' recursively to print
subexpressions.

@< Cases for printing... @>=
case list_display:
  { expr_list l=e.e.sublist;
    if (l==NULL) out << "[]";
    else
    { out << '[';
      do {@; out << l->e; l=l->next; out << (l==NULL ? ']' : ',');}
      while (l!=NULL);
    }
  }
  break;

@ Destroying lists of expressions will be done in a function callable from the
parser, as it may need to discard tokens holding such lists.

@< Declaration of functions in \Cee-style for the parser @>=
void destroy_exprlist(expr_list l);

@~This function recursively destroys subexpressions, and clear up the nodes of
the list themselves when we are done with them. Note that | l=l->next| cannot
be the final statement in the loop body below.

@< Definitions of functions in \Cee-style for the parser @>=
void destroy_exprlist(expr_list l)
{@; while (l!=NULL)
  {@; destroy_expr(l->e); expr_list this_node=l; l=l->next; delete this_node; }
}

@~Destroying a list display is now easily defined.

@< Cases for destroying... @>=
case list_display: destroy_exprlist(e.e.sublist);
break;


@ To build an |exprlist_node|, we just combine an expression with an
|expr_list| inside a freshly created node. To start off the construction, we
shall use |null_expr_list| for the empty list. Often it will be practical to
use right recursive grammar rule that build lists backwards, so we provide a
reversal function to get the proper ordering once the end of the list is
reached. Finally we provide the wrapping function |wrap_expr_list| for list
displays.

@< Declaration of functions in \Cee-style for the parser @>=
extern const expr_list null_expr_list;
expr_list make_exprlist_node(expr e, expr_list l);
expr_list reverse_expr_list(expr_list l);
expr wrap_list_display(expr_list l);


@~The definition of |make_expr_list| is quite trivial, and |reverse_expr_list|
is easy as well. Of the two possible and equivalent list reversal paradigms,
we use the ``hold the head'' style which starts with |t=e|; the other option
is ``hold the tail'', which would start with |t=e->next|. Either one would do
just as well. We do not call |reverse_expr_list| from |wrap_expr_list| although
the two are usually combined, since whether or not the list should be reversed
can only be understood when the grammar rules are given.

@< Definitions of functions in \Cee... @>=
const expr_list null_expr_list=NULL;
expr_list make_exprlist_node(expr e, expr_list l)
{@; expr_list n=new exprlist_node; n->e=e; n->next=l; return n; }
expr_list reverse_expr_list(expr_list l)
{ expr_list r=NULL;
  while (l!=NULL) {@; expr_list t=l; l=t->next; t->next=r; r=t; }
  return r;
}
expr wrap_list_display(expr_list l)
{@; expr result; result.kind=list_display; result.e.sublist=l; return result;
}

@ Besides list displays in which all types must agree, we have tuple displays.
These are rather different from a type-checking point of view, but for
constructing the parse tree, there is no difference, we just need a different
tag to mark the distinction.

@< Enumeration tags for |expr_kind| @>= tuple_display, @[@]

@ We must add cases the appropriate switches to handle tuple displays.
Printing a tuple display is similar to printing a list display, but using
parentheses instead of brackets.

@< Cases for printing... @>=
case tuple_display:
  { expr_list l=e.e.sublist;
    if (l==NULL) out << "()";
    else
    { out << '(';
      do {@; out << l->e; l=l->next; out << (l==NULL ? ')' : ',');}
      while (l!=NULL);
    }
  }
  break;

@~Destroying a tuple display is the same as destroying a list display.

@< Cases for destroying... @>=
case tuple_display: destroy_exprlist(e.e.sublist);
break;

@ To make tuple displays, we use a function similar to that for list displays.
@< Declaration of functions in \Cee-style for the parser @>=
expr wrap_tuple_display(expr_list l);

@~In fact the only difference is the tag inserted.
@< Definitions of functions in \Cee... @>=
expr wrap_tuple_display(expr_list l)
{@; expr result; result.kind=tuple_display; result.e.sublist=l; return result;
}

@*1 Function applications.
Another recursive type of expression is the function application.

@< Typedefs... @>=
typedef struct application_node* app;

@~Now that we have tuples, a function application just takes one argument, so
an |application_node| contains an identifier tag and an argument expression.
Later the identifier tag will be replaced by an arbitrary expression, but for
the moment only named functions are allowed.

@< Structure and typedef declarations for types built upon |expr| @>=
struct application_node {@; id_type fun; expr arg; };

@ The tag used for expressions that will invoke a built-in function is
|function_call|. It is used for operator invocations as well.

@< Enumeration tags for |expr_kind| @>= function_call, @[@]

@~An evident (and unique) category of |expr| values with an |app| as parsing
value is a function call.

@< Variants... @>=
app call_variant;

@ To print a function call, we look up the name of the function in the global
hash table, and either print a tuple display or a single expression enclosed
in parentheses for which we call the operator~`|<<|' recursively. We do not
attempt to reconstruct infix formulae.

@h "lexer.h"
@< Cases for printing... @>=
case function_call:
  { app a=e.e.call_variant;
    out << main_hash_table->name_of(a->fun);
    expr arg=a->arg;
    if (arg.kind==tuple_display) out << arg;
    else out << '(' << arg << ')';
  }
  break;

@~Here we clean up the argument expression, and then the node for the function
call itself.

@< Cases for destroying... @>=
case function_call:
  destroy_expr(e.e.call_variant->arg);
  delete e.e.call_variant;
  break;


@ To build an |application_node|, we combine the function identifier with an
|expr_list| for the argument. For the moment we do the packing or unpacking
(in case of a single argument) of the argument list here, rather than via the
syntax; the latter option would allow avoiding to pack singleton lists. But
the current method should be compatible with providing multiple arguments as a
single tuple value.

@< Declaration of functions in \Cee-style for the parser @>=
expr make_application_node(id_type f, expr_list args);

@~Here for once there is some work to do. If a singleton argument list is
provided, the argument expression must be picked from it, but in all other
cases the argument list must be made into a tuple display. Note that it is
convenient here that |wrap_tuple_display| does not reverse the list, since
this is already done by the parser before calling |make_application_node|.

@< Definitions of functions in \Cee... @>=
expr make_application_node(id_type f, expr_list args)
{ app a=new application_node; a->fun=f;
  if (args!=NULL && args->next==NULL) // a single argument
  { a->arg=args->e; delete args; }
  else a->arg=wrap_tuple_display(args);
  expr result; result.kind=function_call; result.e.call_variant=a;
  return result;
}

@* Other functions callable from the parser.
Here is one function that is not so much a parsing function but just a
wrapper function enabling the parser to interrogate the identifier table.

@< Declaration of functions in \Cee-style for the parser @>=
short lookup_identifier(const char*);

@~The parser will only call this with string constants, so we can use the
|match_literal| method.

@< Definitions of functions in \Cee-style for the parser @>=
id_type lookup_identifier(const char* name)
{@; return atlas::interpreter::main_hash_table->match_literal(name); }

@ The next functions are declared here, because the parser needs to see these
declarations in \Cee-style, but they are defined in in the file
\.{evaluator.w}, since that is where the functionality is available, and we do
not want to make the current compilation unit depend on \.{evaluator.h}.

The function |global_set_identifier| is a temporary feature to handle defining
identifiers while there is not yet a general structure for identifier
definitions within expressions (so that this cannot yet be handled as a part
of evaluation expressions). The left hand side~|ids| refers to a list of
identifiers (this is guaranteed by the parser), which are only represented as
an |expr_list| because this avoids defining a separate type.

@< Declaration of functions in \Cee-style for the parser @>=
void global_set_identifier(expr_list ids, expr e);
void show_ids();
void type_of_expr(expr e);


@* Index.

% Local IspellDict: default