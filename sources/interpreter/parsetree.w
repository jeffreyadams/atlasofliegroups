\def\emph#1{{\it#1\/}}

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
@< Declarations of \Cpp\ functions @>@;
  }@;
}@;

#endif
#endif

@ The main file \.{parsetree.cpp} contains the implementations of the
functions that are needed to build the parse tree. Since these are to be
called from the parser, we declare them all to be callable from~\Cee.

@h "parsetree.h"
@c
@< Definitions of constants... @>@;
namespace atlas
{ namespace interpreter
  {
extern "C" {@;
@< Definitions of functions in \Cee-style for the parser @>@;
}@;
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

@< Declarations of functions in \Cee-style for the parser @>@;

@ While we are defining functions to parse expressions, we shall also define
a function to print the expressions once parsed; this provides a useful test
to see if what we have read in corresponds to what was typed, and this
functionality will also be used in producing error messages.

@< Declarations of \Cpp\ functions @>=
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

@< Declarations of functions in \Cee-style for the parser @>=
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

@ To print an integer denotation we just print its variant field; for string
denotations we do the same but enclosed in quotes, while for Boolean
denotations we reproduce the keyword that gives the denotation.

@< Cases for printing... @>=
case integer_denotation: out << e.e.int_denotation_variant; break;
case string_denotation:
  out << '"' << e.e.str_denotation_variant << '"'; break;
case boolean_denotation:
  out << (e.e.int_denotation_variant!=0 ? "true" : "false"); break;

@~When a string denotation is destroyed, we free the string that was created
by the lexical scanner.

@< Cases for destroying... @>=
case integer_denotation: case boolean_denotation: break;
case string_denotation: delete[] e.e.str_denotation_variant; break;

@ To build the node for denotations, we provide the functions below.

@< Declarations of functions in \Cee-style for the parser @>=
expr make_int_denotation (int val);
expr make_string_denotation(char* val);
expr make_bool_denotation(int val);

@~The definition of these functions is quite trivial, as will be typical for
node-building functions.

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

@ To build the node for applied identifiers we provide the function
|make_applied_identifier|.

@< Declarations of functions in \Cee-style for the parser @>=
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

@ Before we go on to use this type, let us define a simple function for
calculating the length of the list; it will actually be used by the evaluator
rather than by the parser, so we declare it for \Cpp~use.

@< Declarations of \Cpp\ functions @>=
size_t length(expr_list l);

@~The definition is no surprise.

@< Definitions of \Cpp\ functions @>=
size_t length(expr_list l) @+
{@; size_t n=0; while(l!=NULL) {@; ++n; l=l->next;}
  return n;
}

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

@< Declarations of functions in \Cee-style for the parser @>=
void destroy_exprlist(expr_list l);

@~This function recursively destroys subexpressions, and cleans up the nodes of
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

@< Declarations of functions in \Cee-style for the parser @>=
extern const expr_list null_expr_list;
expr_list make_exprlist_node(expr e, expr_list l);
expr_list reverse_expr_list(expr_list l);
expr wrap_list_display(expr_list l);


@~The definition of |make_expr_list| is quite trivial.

@< Definitions of constants visible to the parser @>=
const expr_list null_expr_list=NULL;

@~The function |reverse_expr_list| is easy as well. Of the two possible and
equivalent list reversal paradigms, we use the ``hold the head'' style which
starts with |t=l|; the other option is ``hold the tail'', which would start
with |t=l->next|. Either one would do just as well. We do not call
|reverse_expr_list| from |wrap_expr_list| although the two are usually
combined, since whether or not the list should be reversed can only be
understood when the grammar rules are given.

@< Definitions of functions in \Cee... @>=
expr_list make_exprlist_node(expr e, expr_list l)
{@; expr_list n=new exprlist_node; n->e=e; n->next=l; return n; }
expr_list reverse_expr_list(expr_list l)
{ expr_list r=null_expr_list;
  while (l!=null_expr_list) {@; expr_list t=l; l=t->next; t->next=r; r=t; }
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
@< Declarations of functions in \Cee-style for the parser @>=
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
Function and argument can be arbitrary expressions.

@< Structure and typedef declarations for types built upon |expr| @>=
struct application_node {@; expr fun; expr arg; };

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

@< Cases for printing... @>=
case function_call:
{ app a=e.e.call_variant;
  expr fun=a->fun,arg=a->arg;
  if (fun.kind==applied_identifier) out << fun;
  else out << '(' << fun << ')';
  if (arg.kind==tuple_display) out << arg;
  else out << '(' << arg << ')';
}
break;

@~Here we clean up function and argument, and then the node for the function
call itself.

@< Cases for destroying... @>=
case function_call:
  destroy_expr(e.e.call_variant->fun);
  destroy_expr(e.e.call_variant->arg);
  delete e.e.call_variant;
break;


@ To build an |application_node|, we combine the function identifier with an
|expr_list| for the argument. For the moment we do the packing or unpacking
(in case of a single argument) of the argument list here, rather than via the
syntax; the latter option would allow avoiding to pack singleton lists. But
the current method should be compatible with providing multiple arguments as a
single tuple value.

@< Declarations of functions in \Cee-style for the parser @>=
expr make_application_node(expr f, expr_list args);

@~Here for once there is some work to do. If a singleton argument list is
provided, the argument expression must be picked from it, but in all other
cases the argument list must be made into a tuple display. Note that it is
convenient here that |wrap_tuple_display| does not reverse the list, since
this is already done by the parser before calling |make_application_node|.

@< Definitions of functions in \Cee... @>=
expr make_application_node(expr f, expr_list args)
{ app a=new application_node; a->fun=f;
  if (args!=NULL && args->next==NULL) // a single argument
  {@; a->arg=args->e; delete args; }
  else a->arg=wrap_tuple_display(args);
  expr result; result.kind=function_call; result.e.call_variant=a;
  return result;
}

@*1 Identifier patterns.
%
For each case where local identifiers will be introduced (like let-expressions
of function headings) we shall in fact allow more general patterns. A defining
occurrence of an identifier \\{ident} may be replaced by some tuple, say
$(x,y,z)$ (which assumes the value to be bound will be a $3$-tuple), or it
could be both, as in \\{ident}:$(x,y,z)$; for the nested identifiers the same
options apply recursively, and in addition the identifier may be suppressed
altogether, to allow such partial tagging of components as \\{ident}:$(,,z)$.
To accommodate such possibilities we introduce the following recursive types.
While not visible from these definitions, out grammar will ensure that
whenever a sublist is present, it has at least 2 nodes.

@< Typedefs... @>=
typedef struct pattern_node* patlist;
struct id_pat
{ patlist sublist;
  id_type name;
  unsigned char kind; /* bit 0: has name, bit 1: has sublist */
};
struct pattern_node
{@; patlist next;
  struct id_pat body;
};

@ These types do not themselves represent a variant of |expru|, but will be
used inside such variants. We can already provide a printing function.

@< Declarations of \Cpp... @>=
std::ostream& operator<< (std::ostream& out, const id_pat& p);

@~Only parts whose presence is indicated in |kind| are printed. Although the
grammar for patterns only allows sublists of length at least~2, we also visit
this code when printing user-defined functions, so we take care that even if
|sublist| is marked as present, it might be null. However, the list cannot be
of length~1, so having put aside that case we can print a first comma
unconditionally.

@< Definitions of \Cpp...@>=
std::ostream& operator<< (std::ostream& out, const id_pat& p)
{ if ((p.kind & 0x1)!=0)
    out << main_hash_table->name_of(p.name);
  if (p.kind==0x3) // both parts present
    out << ':';
  if ((p.kind & 0x2)!=0)
    if (p.sublist==NULL)
      out << "()";
    else
    { out << '(' << p.sublist->body << ',';
      for (patlist l=p.sublist->next; l!=NULL; l=l->next)
      out << l->body << (l->next!=NULL ? ',' : ')');
    }
  return out;
}

@ The function to build a node takes a pointer to a structure with the
contents of the current node; in practice this is the address of a local
variable in the parser. Patterns also need cleaning up, which
|destroy_pattern| and |destroy_id_pat| will handle, and reversal as handled by
|reverse_patlist|.

@< Declarations in \Cee... @>=
patlist make_pattern_node(patlist next,struct id_pat* body);
void destroy_pattern(patlist p);
void destroy_id_pat(struct id_pat* p);
patlist reverse_patlist(patlist p);

@ The function just assembles the pieces. In practice the |next| pointer will
point to previously parsed nodes, so (as usual) reversal will be necessary.
Being bored, we add a variation on list reversal.

@< Definitions of functions in \Cee... @>=
patlist make_pattern_node(patlist next,struct id_pat* body)
{@; patlist l=new pattern_node; l->next=next; l->body=*body; return l; }
@)
void destroy_pattern(patlist p)
{@; while (p!=NULL)
  {@; patlist q=p; p=p->next; destroy_id_pat(&q->body); delete q; }
}
@)
void destroy_id_pat(struct id_pat* p)
{@; if ((p->kind & 0x2)!=0)
    destroy_pattern(p->sublist);
}
@)
patlist reverse_patlist(patlist p)
{@; patlist q=NULL;
  while (p!=NULL)
  {@; patlist t=q; q=p; p=q->next; q->next=t; }
  return q;
}

@*1 Let expressions.
We introduce let-expressions that introduce and bind local identifiers, as a
first step towards having use defined functions. Indeed let-expressions will
be implemented (initially) as a user-defined function that is immediately
called with the values bound in the let-expression. The reason to do this
before considereing user-defined functions in general is that this avoids for
now having to specify in the user program the types of the parameters of the
function, these type being determined by the values provided.

@< Typedefs... @>=
typedef struct let_expr_node* let;

@~After parsing, let-expression will have a single let-bindings followed by a
body giving the value to be returned. During parsing however, we may form a
list of let-bindings that will be converted into one with a tuple as left hand
side. We therefore define a node for the list of bindings, and a structure for
the whole let-expression that contains one binding and a body.

@< Structure and typedef declarations for types built upon |expr| @>=
typedef struct let_node* let_list;
struct let_node {@; struct id_pat pattern; expr val; let_list next; };
struct let_expr_node {@; struct id_pat pattern; expr val; expr body; };

@ The tag used for let-expressions is |let_expr|.

@< Enumeration tags for |expr_kind| @>= let_expr, @[@]

@~Without surprise there is a class of |expr| values with a |let| as parsing
value.

@< Variants... @>=
let let_variant;

@ To print a let-expression we do the obvious things, wihtout worrying about
parentheses; this should be fixed (for all printing routines).

@< Cases for printing... @>=
case let_expr:
{ let lexp=e.e.let_variant;
  out << "let " << lexp->pattern << '=' << lexp->val <<  " in " << lexp->body;
}
break;

@ Destroying lists of declarations will be done in a function callable from the
parser, like |destroy_exprlist|.

@< Declarations of functions in \Cee-style for the parser @>=
void destroy_letlist(let_list l);

@~Like |destroy_exprlist|, this function recursively destroys nodes, and the
expressions they contain.

@< Definitions of functions in \Cee-style for the parser @>=
void destroy_letlist(let_list l)
{ while (l!=NULL)
    @/{@; let_list p=l; l=l->next;
      destroy_id_pat(&p->pattern);
      destroy_expr(p->val);
      delete p;
    }
}

@~Here we clean up the declaration, and then the body of the let-expression.

@< Cases for destroying... @>=
case let_expr:
{ let lexp=e.e.let_variant;
  destroy_id_pat(&lexp->pattern);
  destroy_expr(lexp->val);
  destroy_expr(lexp->body);
  delete lexp;
}
break;

@ For building let-expressions, two functions will be defined. The function
|add_let_node| adds one declaration to a list (it is called with |prev==NULL|
for the first clause), while |make_let_expr_node| wraps up the let-expression.

@< Declarations of functions in \Cee-style for the parser @>=
let_list add_let_node(let_list prev, struct id_pat pattern, expr val);
expr make_let_expr_node(let_list decls, expr body);

@~In |make_let_expr_node| we may need to convert multiple declarations to one,
in which case we must take care to reverse the order, since the last one added
will be at the head of the list. Fortunately it is actually easier to build a
merged list in reverse order. We provide local exception safety here, although
we realise that the parser function being written in \Cee, it will not be able
to do any cleaning up of values referred to in its stack, in case of a thrown
exception (here, or in calls to make-functions elsewhere). In fact the problem
is not just one of programming language: even a parser generated as \Cpp~code
would not help without special provisions (such as an exception handler in the
parser function), since the |union| used for values on the parsing stack
cannot have objects that handle their own cleaning up as members. Rather, this
situation should probably be corrected by using a different allocation
strategy while building the parse tree (avoiding explicit calls to |new| but
allocating from local storage pools that can explicitly be emptied), in which
case all cleaning up in the code below should also be removed.

@< Definitions of functions in \Cee... @>=
let_list add_let_node(let_list prev,struct id_pat pattern, expr val)
{@; let_list l=new let_node;
  l->pattern=pattern; l->val=val; l->next=prev;
  return l;
}
@)
expr make_let_expr_node(let_list decls, expr body)
{ let l=new let_expr_node; l->body=body;
  expr result; result.kind=let_expr; result.e.let_variant=l;
  try
  { if (decls->next==NULL) // single declaration
    @/{@; l->pattern=decls->pattern; l->val=decls->val;
      delete decls;
    }
    else
    { l->pattern.kind=0x2; l->pattern.sublist=NULL;
      l->val=wrap_tuple_display(null_expr_list);
      while (decls!=NULL)
      { l->pattern.sublist =
          make_pattern_node(l->pattern.sublist,&decls->pattern);
        l->val.e.sublist=make_exprlist_node(decls->val,l->val.e.sublist);
        let_list p=decls; decls=p->next; delete p;
      }
    }
    return result;
  }
  catch(...)
  {@; destroy_expr(result);
    throw;
  }
}

@*1 Types and user-defined functions.
%
One reason let-expressions were introduced before user-defined functions, is
that it avoids to problem of having to specify types in the user program.
Inevitably we have to deal with that though, since for a function definition
there is no way to know the type of the arguments with certainty, unless the
user specifies them (at least this is true if we want to allow such things as
type coercion and function overloading to be possible, so that types cannot be
deduced by analysis of the \emph{usage} of the parameters in the function
only). The syntax for types is easy enough, but we have to decide what kind of
object the parser will use to represent types specified by the user. The
evaluator has an internal type |struct type_declarator@;| to represent them,
but the \Cpp\ definition of that type cannot be understood by the parser. So
we tried to manage with pointers to incomplete types. That does not work
either, because |struct type_declarator@;| is actually defined within a
|namespace|, so we cannot give the correct name and be understood in \Cee,
while if we lie about the name of the |struct| then \Cpp\ will complain about
converting between pointers to different structures. So, grudgingly, we shall
cast to and from pointers to |void|.

@< Typedefs that are required... @>=
typedef void* ptr;

@~The functions declared below provide an interface to routines defined in
the module \.{evaluator.w}.

@< Declarations in \Cee... @>=
ptr mk_type_singleton(ptr t);
ptr mk_type_list(ptr t,ptr l);
ptr mk_prim_type(int p);
ptr mk_row_type(ptr c);
ptr mk_tuple_type(ptr l);
ptr mk_function_type(ptr a,ptr r);
@)
void destroy_type(ptr t);
void destroy_type_list(ptr t);

@ Some other functions that must also convert pointers from |ptr|, and which
are therefore defined in \.{evaluator.w}, are not used directly by the parser;
their declarations can use \Cpp\ types.

@< Declarations of \Cpp... @>=
ptr first_type(ptr typel);
std::ostream& print_type(std::ostream& out, ptr type);

@ For user-defined functions we shall use a structure |lambda_node|.
@< Typedefs that are required... @>=
typedef struct lambda_node* lambda;

@~It contains a pattern for the formal parameter(s), its type (a void pointer
that actually points to a |type_declarator| structure defined
in \.{evaluator.w}), and an  expression (the body of the function).

@< Structure and typedef... @>=
struct lambda_node
{@; struct id_pat pattern; ptr arg_type; expr body; };

@ The tag used for user-defined functions is |lambda_expr|.
@< Enumeration tags... @>=
lambda_expr,@[@]

@ We introduce the variant of |union expru@;| as usual.
@< Variants of |union... @>=
lambda lambda_variant;

@ We must take care of printing lambda expressions; we avoid a double set of
parentheses.

@< Cases for printing... @>=
case lambda_expr:
{ lambda fun=e.e.lambda_variant;
  if ((fun->pattern.kind&0x1)!=0)
    out << '(' << fun->pattern << ')';
  else
    out << fun->pattern;
  out << ':' << fun->body;
}
break;

@ And we must of course take care of destroying lambda expressions, which just
call handler functions.

@< Cases for destroying an expression |e| @>=
case lambda_expr:
{ lambda fun=e.e.lambda_variant;
  destroy_id_pat(&fun->pattern);
  destroy_type(fun->arg_type);
  destroy_expr(fun->body);
}
break;

@ Finally there is as usual a function for constructing a node, to be called
by the parser.

@< Declarations of functions... @>=
expr make_lambda_node(patlist patl, ptr typel, expr body);

@~There is a twist in building a lambda node, in that it is passed lists of
patterns and types rather than single ones. We must distinguish the case of a
singleton, where the head node must be unpacked, and the multiple case, where
a tuple pattern and type must be wrapped up from the lists.

@< Definitions of functions... @>=
expr make_lambda_node(patlist patl, ptr typel, expr body)
{ lambda fun=new lambda_node; fun->body=body;
  if (patl!=NULL and patl->next==NULL)
  { fun->pattern=patl->body; delete patl; // clean up node
    fun->arg_type = first_type(typel);
  }
  else
  @/{@; fun->pattern.kind=0x2; fun->pattern.sublist=patl;
    fun->arg_type=mk_tuple_type(typel);
  }
  expr result; result.kind=lambda_expr; result.e.lambda_variant=fun;
  return result;
}

@*1 Control structures.

@*2 Conditional expressions.
Of course we need if-then-else expressions.

@< Typedefs... @>=
typedef struct conditional_node* cond;

@~The parser handles \&{elif} constructions, so we only need to handle the
basic two-branch case.

@< Structure and typedef declarations for types built upon |expr| @>=
struct conditional_node
 {@; expr condition; expr then_branch; expr else_branch; };

@ The tag used for these expressions is |conditional_expr|.

@< Enumeration tags for |expr_kind| @>= conditional_expr, @[@]

@~The variant of |expr| values with an |cond| as parsing value is tagged
|if_variant|.

@< Variants... @>=
cond if_variant;

@ To print a conditional expression at parser level, we shall not
reconstruct \&{elif} constructions.

@< Cases for printing... @>=
case conditional_expr:
{ cond c=e.e.if_variant; out << " if " << c->condition
  << " then " << c->then_branch << " else " << c->else_branch << " fi ";
}
break;

@~Here we clean up constituent expressions and then the node for the conditional
call itself.

@< Cases for destroying... @>=
case conditional_expr:
  destroy_expr(e.e.if_variant->condition);
  destroy_expr(e.e.if_variant->then_branch);
  destroy_expr(e.e.if_variant->else_branch);
  delete e.e.if_variant;
break;


@ To build an |conditional_node|, we define a function as usual.
@< Declarations of functions in \Cee-style for the parser @>=
expr make_conditional_node(expr c, expr t, expr e);

@~It is entirely straightforward.

@< Definitions of functions in \Cee... @>=
expr make_conditional_node(expr c, expr t, expr e)
{ cond n=new conditional_node; n->condition=c;
  n->then_branch=t; n->else_branch=e;
@/expr result; result.kind=conditional_expr; result.e.if_variant=n;
  return result;
}

@*1 Array subscriptions.
%
We want to be able to select components from array structures (lists, vectors,
matrices), so we define a subscription expression. If there are multiple
indices, these can be realised as a subscription by a tuple expression, so we
define only one type so subscription expression.

@< Typedefs... @>=
typedef struct subscription_node* sub;

@~In a subscription the array and the index(es) can syntactically be arbitrary
expressions (although the latter should have as type integer, or a tuple of
integers).

@< Structure and typedef declarations for types built upon |expr| @>=
struct subscription_node {@; expr array; expr index; };

@ Here is the tag used for subscriptions.

@< Enumeration tags for |expr_kind| @>= subscription, @[@]

@~And here is the corresponding variant of the |union expru@;|.
value is a function call.

@< Variants... @>=
sub subscription_variant;

@ To print a subscription, we just print the expression of the array, followed
by the expression for the index in brackets. As an exception, the case of a
tuple display as index is handled separately, in order to avoid having
parentheses directly inside the brackets.

@h "lexer.h"
@< Cases for printing... @>=
case subscription:
{ sub s=e.e.subscription_variant; out << s->array << '[';
  expr i=s->index;
  if (i.kind!=tuple_display) out << i;
  else
  { expr_list l=i.e.sublist;
    if (l!=NULL)
    {@; out << l->e;
      while ((l=l->next)!=NULL) out << ',' << l->e;
    }
  }
  out << ']';
}
break;

@~Here we recursively destroy both subexpressions, and then the node for the
subscription call itself.

@< Cases for destroying... @>=
case subscription:
  destroy_expr(e.e.subscription_variant->array);
  destroy_expr(e.e.subscription_variant->index);
  delete e.e.subscription_variant;
break;


@ To build an |subscription_node|, we simply combine the array and the index
part.

@< Declarations of functions in \Cee-style for the parser @>=
expr make_subscription_node(expr a, expr i);

@~This is straightforward, as usual.

@< Definitions of functions in \Cee... @>=
expr make_subscription_node(expr a, expr i)
{ sub s=new subscription_node; s->array=a; s->index=i;
  expr result; result.kind=subscription; result.e.subscription_variant=s;
  return result;
}

@*1 Cast expressions.
%S
These are very simple expressions consisting of a type and an expression,
which is forced to be of that type.

@< Typedefs... @>=
typedef struct cast_node* cast;

@~The type contained in the cast is represented by a void pointer, for the
known reasons.

@< Structure and typedef declarations for types built upon |expr| @>=
struct cast_node {@; ptr type; expr exp; };

@ The tag used for casts is |cast_expr|.

@< Enumeration tags for |expr_kind| @>=
cast_expr, @[@]

@ And there is of course a variant of |expr_union| for casts.
@< Variants of ... @>=
cast cast_variant;

@ Printing cast expressions follows their input syntax.

@< Cases for printing... @>=
case cast_expr:
{@; cast c = e.e.cast_variant;
  print_type(out,c->type) << ':' << c->exp ;
}
break;

@ Casts are built by |make_cast|.

@< Declarations of functions in \Cee... @>=
expr make_cast(ptr type, expr exp);

@~No surprises here.

@< Definitions of functions in \Cee...@>=
expr make_cast(ptr type, expr exp)
{ cast c=new cast_node; c->type=type; c->exp=exp;
@/ expr result; result.kind=cast_expr; result.e.cast_variant=c;
   return result;
}

@ Eventually we want to rid ourselves from the cast.

@< Cases for destr... @>=
case cast_expr:
  destroy_expr(e.e.cast_variant->exp); delete e.e.cast_variant;
break;

@*1 Assignment statements.
%
Simple assignment statements are quite simple as expressions.

@< Typedefs... @>=
typedef struct assignment_node* assignment;

@~In a simple assignment the left hand side is just an identifier.

@< Structure and typedef declarations for types built upon |expr| @>=
struct assignment_node {@; id_type lhs; expr rhs; };

@ The tag used for assignment statements is |ass_stat|.

@< Enumeration tags for |expr_kind| @>=
ass_stat, @[@]

@ And there is of course a variant of |expr_union| for assignments.
@< Variants of ... @>=
assignment assign_variant;

@ Printing assignment statements is absolutely straightforward.

@< Cases for printing... @>=
case ass_stat:
{@; assignment ass = e.e.assign_variant;
  out << main_hash_table->name_of(ass->lhs) << ":=" << ass->rhs ;
}
break;

@ Assignment statements are built by |make_assignment|.

@< Declarations of functions in \Cee... @>=
expr make_assignment(id_type lhs, expr rhs);

@~It does what one would expect it to (except for those who expect their
homework assignment made).

@< Definitions of functions in \Cee...@>=
expr make_assignment(id_type lhs, expr rhs)
{ assignment a=new assignment_node; a->lhs=lhs; a->rhs=rhs;
@/ expr result; result.kind=ass_stat; result.e.assign_variant=a;
   return result;
}

@ What is made must eventually be unmade (even assignments).

@< Cases for destr... @>=
case ass_stat:
  destroy_expr(e.e.assign_variant->rhs); delete e.e.assign_variant;
break;

@*2 Component assignments.
%
We have special expressions for assignments to a component.

@< Typedefs... @>=
typedef struct comp_assignment_node* comp_assignment;

@~In a component assignment has for the left hand side an identifier and an
index.

@< Structure and typedef declarations for types built upon |expr| @>=
struct comp_assignment_node {@; id_type aggr; expr index; expr rhs; };

@ The tag used for assignment statements is |comp_ass_stat|.

@< Enumeration tags for |expr_kind| @>=
comp_ass_stat, @[@]

@ And there is of course a variant of |expr_union| for assignments.
@< Variants of ... @>=
comp_assignment comp_assign_variant;

@ Printing component assignment statements follow the input syntax.

@< Cases for printing... @>=
case comp_ass_stat:
{@; comp_assignment ass = e.e.comp_assign_variant;
  out << main_hash_table->name_of(ass->aggr) << '[' << ass->index << "]:="
      << ass->rhs ;
}
break;

@ Assignment statements are built by |make_assignment|, which for once does
not simply combine the expression components, because for reason of parser
generation the array and index will have already been combined before this
function can be called.

@< Declarations of functions in \Cee... @>=
expr make_comp_ass(expr lhs, expr rhs);

@~Here we have to take the left hand side apart a bit, and clean up its node.

@< Definitions of functions in \Cee...@>=
expr make_comp_ass(expr lhs, expr rhs)
{ comp_assignment a=new comp_assignment_node;
@/a->aggr=lhs.e.subscription_variant->array.e.identifier_variant;
  a->index=lhs.e.subscription_variant->index;
  delete lhs.e.subscription_variant;
  a->rhs=rhs;
@/expr result; result.kind=comp_ass_stat; result.e.comp_assign_variant=a;
  return result;
}

@~Destruction one the other hand is as straightforward as usual.

@< Cases for destr... @>=
case comp_ass_stat:
  destroy_expr(e.e.comp_assign_variant->index);
  destroy_expr(e.e.comp_assign_variant->rhs);
  delete e.e.comp_assign_variant;
break;

@*1 Sequence statements.
%
Having assignments statements, it is logical to be able to build a sequence of
expressions (statements) as well, retaining the value only of the final one.

@< Typedefs... @>=
typedef struct sequence_node* sequence;

@~Since control structures and let-expressions tend to break up long chains,
we do not expect their average length to be very great. So we build up
sequences by chaining pairs, instead of storing a vector of expressions: at
three pointers overhead per vector, a chained representation is more compact
as long as the average length of a sequence is less than~$5$ expressions.

@< Structure and typedef declarations for types built upon |expr| @>=
struct sequence_node {@; expr first; expr last; };

@ The tag used for sequence statements is |seq_expr|.

@< Enumeration tags for |expr_kind| @>=
seq_expr, @[@]

@ And there is of course a variant of |expr_union| for sequences.
@< Variants of ... @>=
sequence sequence_variant;

@ Printing sequences is absolutely straightforward.

@< Cases for printing... @>=
case seq_expr:
{@; sequence seq = e.e.sequence_variant;
  out << seq->first << ";" << seq->last ;
}
break;

@ Sequences are built by |make_sequence|.

@< Declarations of functions in \Cee... @>=
expr make_sequence(expr first, expr last);

@~It does what one would expect it to.

@< Definitions of functions in \Cee...@>=
expr make_sequence(expr first, expr last)
{ sequence s=new sequence_node; s->first=first; s->last=last;
@/ expr result; result.kind=seq_expr; result.e.sequence_variant=s;
   return result;
}

@ Finally sequence nodes need destruction, like everything else.

@< Cases for destr... @>=
case seq_expr:
  destroy_expr(e.e.sequence_variant->first);
  destroy_expr(e.e.sequence_variant->last);
  delete e.e.sequence_variant;
break;

@* Other functions callable from the parser.
Here are some functions that are not so much a parsing functions as just
wrapper functions enabling the parser to call \Cpp~functions.

@< Declarations of functions in \Cee-style for the parser @>=
short lookup_identifier(const char*);
void include_file();

@~The parser will only call this with string constants, so we can use the
|match_literal| method.

@< Definitions of functions in \Cee-style for the parser @>=
id_type lookup_identifier(const char* name)
{@; return main_hash_table->match_literal(name); }

@~To include a file, we call the |push_file| method from the input buffer,
providing a file name that was remembered by the lexical analyser.

@< Definitions of functions in \Cee-style for the parser @>=
void include_file()
{@; main_input_buffer->push_file(lex->scanned_file_name()); }

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

@< Declarations of functions in \Cee-style for the parser @>=
void global_set_identifier(expr_list ids, expr e);
void show_ids();
void type_of_expr(expr e);


@* Index.

% Local IspellDict: british
