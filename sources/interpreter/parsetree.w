\def\emph#1{{\it#1\/}}
\chardef\pow = `\^

@* Building the parse tree.
This is the program unit \.{parsetree} which produces the implementation file
\.{parsetree.cpp} and the header file \.{parsetree.h}; the latter is read in
by \&{\#include} from the \Cee-file \.{parser.tab.c} generated
from \.{parser.y}. It used to be the case that this file was compiled by
a \Cee\ compiler, and therefore produced functions with \Cee-language linking.
As a consequence some jumping through hoops was necessary to cleanly integrate
them into the \Cpp\ program. However it turns out to be possible to compile
the file \.{parser.tab.c} by a \Cpp-compiler, which makes all linkage to be
for \Cpp, and removes any need for |extern "C"| declarations. For the moment
these declarations are simply removed, and the structure of this file still
carries a legacy of the original design.

This file also defines some other functions defined here are not
used in the parser, but can be used by other modules written in~\Cpp; their
declaration is separated fro historic reasons only.

@( parsetree.h @>=
#ifndef PARSETREE_H
#define PARSETREE_H

#include <iostream>
#include <memory>
#include <string>
#include "buffer.h" // for |Hash_table|
#include "sl_list.h"
namespace atlas
{
  namespace interpreter
  {
@< Declarations for the parser @>@;

@< Declaration of functions not for the parser @>@;
  }@;
}@;

#endif


@ The main file \.{parsetree.cpp} contains the implementations of the
functions that are needed to build the parse tree. Since these are to be
called from the parser, we used to declare them all to be callable from~\Cee;
however now the parser is compiled as \Cpp\ code, that is no longer an issue,
and we define everything in a \Cpp\ namespace as usual.

@h "parsetree.h"
@c
namespace atlas
{ namespace interpreter
  {
@< Definitions of functions for the parser @>@;
@< Definitions of functions not for the parser @>@;
  }@;
}@;

@ For a large part the declarations for the parser consist of the recursive
definition of the type |expr|. While that used to be a POD type used directly
on the parser stack, the inflexibility of this solution motivated changing
this to placing row pointers

@< Declarations for the parser @>=
@< Type declarations needed in definition of |struct expr@;| @>@;

enum expr_kind @+
 { @< Enumeration tags for |expr_kind| @>@;@; @+no_expr };
struct expr {
  expr_kind kind;
  union {@; @< Variants of |union expru @;| @>@; };
@)
  @< Methods of |expr| @>@;
};
typedef expr* expr_p; // raw pointer type for use on parser stack
typedef std::unique_ptr<expr> expr_ptr;
@)
@< Structure and typedef declarations for types built upon |expr| @>@;
@< Declarations of functions for the parser @>@;

@ We start right away declaring |expr| as a |struct|, avoiding complaints that
it is not declared.

@< Type declarations needed in definition of |struct expr@;| @>=
struct expr;


@ When default-constructing an |expr| we set its |kind| to |no_expr|.
@< Methods of |expr| @>=
expr();
~expr(); // defined below using a large |switch| statement

@ While we are defining functions to parse expressions, we shall also define a
function |destroy_expr| to clean up the memory occupied by an expression. It
is classified as a parsing function since it is called amongst others by the
parser when popping off tokens at syntax errors.

@< Declarations of functions for the parser @>=
void destroy_expr(expr_p e);

@~The definition of |destroy_expr| just calls |delete|, which will invoke the
destructor of the |expr| pointed to.

@< Definitions of functions for the parser @>=
void destroy_expr(expr_p p) @+ {@; delete p; }
void destroy_expr_body(const expr& e) @+ {@; e.~expr(); }

@ The actual definition of the destructor is distributed among the different
variants of |expr|.

@< Definitions of functions not for the parser @>=
expr::expr() : kind(no_expr)
@+{}

expr::~expr()
{
  if (kind!=no_expr)
  {
    expr& e = *this;
    switch (kind)
    {@; @< Cases for destroying an expression |e| @>
      @+ case no_expr: {}
    }
    kind = no_expr;
  }
}

@ We define a move constructor and move assignment; our first purpose is to
eliminate the need for copying as much as possible by detecting now forbidden
copy constructions hidden in the code.

@< Methods of |expr| @>=
void set_from (expr& other);
expr (expr&& other);
void operator= (expr&& other);

@ For now we use that |expr| is trivially copyable; this should be changed.

@< Definitions of functions not for... @>=
void expr::set_from (expr& other)
{ assert(kind==no_expr);
  switch (kind=other.kind)
    {@; @< Cases for copying an expression from |other| @>
       @+ case no_expr: {}
    }
  other.kind = no_expr;
}

expr::expr (expr&& other) : @[ expr () @]
{@; set_from(other); }

void expr::operator= (expr&& other)
{
  if (this!=&other)
  {@; this->~expr();
    set_from(other);
  }
}

@ In parallel, we also define a function to print the expressions once parsed;
this provides a useful test to see if what we have read in corresponds to what
was typed, and this functionality will also be used in producing error
messages.

@< Declaration of functions not for the parser @>=
std::ostream& operator<< (std::ostream& out, const expr& e);

@~The definitions of this instance of the operator~`|<<|' are also distributed
among the different variants of |expr| that we shall define.

@< Definitions of functions not for the parser @>=
std::ostream& operator<< (std::ostream& out, const expr& e)
{@; switch (e.kind)
  {@; @< Cases for printing an expression |e| @>
    @+ case no_expr: {}
  }
  return out;
}


@*1 Denotations. The simplest expressions are atomic constants, which we shall
call denotations. There are recognised by the scanner, and either the scanner
or the parser will build an appropriate node for them, which just stores the
constant value denoted. For integer and Boolean denotations, the value itself
will fit comfortably inside the |struct expr@;|, and moreover we can share a
variant because the tag will tell whether to interpret it is integer or
Boolean. For strings we store a |std::string|, which is more complicated than
the |char*@[@]@;| that we used to store, but it is instructive for how the
special member functions of |expr| should handle non-POD variants, which can
neither be assigned to non-initialised memory  without calling a constructor,
nor be left in memory that will be reclaimed without calling a destructor.

@< Variants... @>=

int int_denotation_variant;
std::string str_denotation_variant;

@~Each of the three types of denotation has a tag identifying it.

@< Enumeration tags for |expr_kind| @>=
integer_denotation, string_denotation, boolean_denotation, @[@]

@ For each of these variants there is a corresponding constructor that
placement-constructs the constant value into that variant. For the integer and
Boolean case we might alternatively have assigned to the field, bit not for
the string variant.

@< Methods of |expr| @>=
  expr(int n) : kind(integer_denotation), int_denotation_variant(n) @+{}
  expr(bool b)
   : kind(boolean_denotation), int_denotation_variant(b ? 1 : 0) @+{}
  expr(std::string&& s)
   : kind(string_denotation)
   , str_denotation_variant(std::move(s)) @+{}

@~For more explicit construction of these variants in a dynamically allocated
|expr| object, we provide the functions below.

@< Declarations of functions for the parser @>=
expr_p make_int_denotation (int val);
expr_p make_bool_denotation(bool val);
expr_p make_string_denotation(std::string&& val);

@~The definition of these functions is quite trivial easy, as will be typical
for node-building functions.

@< Definitions of functions for the parser @>=
expr_p make_int_denotation (int val)
  @+{@; return new expr(val); }

expr_p make_bool_denotation(bool val)
  @+{@; return new expr(val); }

expr_p make_string_denotation(std::string&& val)
 @+{@; return new expr(std::move(val)); }

@~For integer and Boolean denotations there is nothing to destroy. For string
denotations however we must destroy the |std::string| object. It was quite a
puzzle to find the right syntax for that, because of what is essentially a bug
in \.{gcc}, namely with |string| in place of |basic_string| the look-up of the
destructor fails.

@s basic_string string

@< Cases for destroying... @>=
case integer_denotation: case boolean_denotation: break;
case string_denotation:
  e.str_denotation_variant.~basic_string<char>(); break;

@ In the |expr::set_from| method we change variants both in |*this| (which was
|no_expr|) and in |other| (which was |no_expr|). This means that for non-POD
type we must combine construction into a variant of |*this| and a destruction
of that variant of |other|. We use move construction for efficiency, and it
will probably leave an empty shell to be destructed, but this does not mean we
can omit the destruction.

@< Cases for copying... @>=
  case integer_denotation:
  case boolean_denotation:
    int_denotation_variant = other.int_denotation_variant; break;
  case string_denotation:
    new (&str_denotation_variant)
    std::string(std::move(other.str_denotation_variant));
    other.str_denotation_variant.~basic_string<char>();
  break;

@ To print an integer denotation we just print its variant field; for Boolean
denotations we reproduce the keyword that gives the denotation, while for string
denotations we print the stored string enclosed in quotes.

@< Cases for printing... @>=
case integer_denotation: out << e.int_denotation_variant; break;
case boolean_denotation:
  out << (e.int_denotation_variant!=0 ? "true" : "false"); break;
case string_denotation:
  out << '"' << e.str_denotation_variant << '"'; break;

@ Another atomic expression is an applied identifier. They use the type
|Hash_table::id_type| (a small integer type) of indices into the table of
identifier names, which we lift out of that class by using a |typedef|.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef Hash_table::id_type id_type;

@~Their tag is |applied_identifier|. An expression that behaves somewhat
similarly is `\.\$', which stands for the last value computed.

@< Enumeration tags for |expr_kind| @>=
applied_identifier,
last_value_computed, @[@]

@ For identifiers we just store their code; for |last_value_computed| nothing
at all.
@< Variants... @>=
id_type identifier_variant;

@ Rather than having separate constructors for these variants, which would be
ambiguous with constructors already defined, we build these variants directly
in the following functions.

@< Declarations of functions for the parser @>=
expr_p make_applied_identifier (id_type id);
expr_p make_dollar();

@~In spite of the absence of dedicated constructors, these function have
rather simple definitions. We split off a function |mk_applied_identifier|
that is like a constructor, because it will be serve again below.

@< Definitions of functions for the parser @>=
expr mk_applied_identifier (id_type id)
{@; expr result; result.kind=applied_identifier;
  result.identifier_variant=id; return result;
}
expr_p make_applied_identifier (id_type id)
 {@; return new expr(mk_applied_identifier(id)); }

expr_p make_dollar ()
{@; expr_p result = new expr;
  result->kind=last_value_computed;
  return result;
}

@~Like for integer and boolean denotations, there is nothing to destroy here.

@< Cases for destroying... @>=
case applied_identifier:
case last_value_computed: break;

@ Having a POD type variant, copying an applied identifier can be done by
assignment.

@< Cases for copying... @>=
case applied_identifier: identifier_variant=other.identifier_variant; break;
case last_value_computed: break;

@~To print an applied identifier, we look it up in the main hash table. We
print \.\$ as the user wrote it.

@< Cases for printing... @>=
case applied_identifier:
  out << main_hash_table->name_of(e.identifier_variant);
break;
case last_value_computed: out << '$'; @q$@> break;

@*1 Expression lists.
%
We shall need expression lists for various purposes.
After long refactoring of the code (notably unmaking |expr| a POD-exclusive
tagged union), we can use the atlas class template |containers::simple_list|
to implement it, instead of having to define a custom list type.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef containers::simple_list<expr> expr_list;
typedef containers::sl_node<expr>* raw_expr_list;
  // raw counterpart for use by parser

@ Any syntactic category whose parsing value is a list of expressions will use
the variant |sublist|.

@< Variants... @>=
expr_list sublist;

@~There are two such syntactic categories: tuple displays and list displays,
which are written with parentheses respectively with square brackets.
These are rather different from a type-checking point of view (in a list
display the types of all components must be the same), but for
constructing the parse tree, there is no difference, we just need a different
tag to mark the distinction.

@< Enumeration tags for |expr_kind| @>= tuple_display, list_display, @[@]

@~We provide a constructor that serves both types of lists.

@< Methods of |expr| @>=
expr(expr_list&& nodes, bool is_tuple)
 : kind(is_tuple ? tuple_display : list_display)
 , sublist(std::move(nodes))
 @+{}

@~Destroying a tuple display or list display is easily defined.

@< Cases for destroying... @>=
  case tuple_display:
  case list_display: e.sublist.~expr_list();
break;

@ Destroying lists of expressions will be done in a function callable from the
parser, as it may need to discard tokens holding such lists.

@< Declarations of functions for the parser @>=
void destroy_exprlist(raw_expr_list l);

@~Its definition is easy; since only the head pointer of the list is raw,
simply deleting it will recursively clear up the whole list.

@< Definitions of functions for the parser @>=
void destroy_exprlist(raw_expr_list l)
@+{@; delete l; }

@ Copying can be obtained by move construction followed by destruction of
(what remains of) the original value. For raw pointers simply assigning
|sublist=other.sublist| would also work, but for smart pointers this is the
only proper way to proceed, even though destructing one after
move-constructing out of it is most probably a no-op.

@< Cases for copying... @>=
  case tuple_display:
  case list_display:
    new (&sublist) expr_list (std::move(other.sublist));
    other.sublist.~expr_list();
  break;

@ To build an |exprlist_node|, we provide a function |make_exprlist_node| to
combine an expression with a |raw_expr_list|. To start off the construction,
one may use |raw_expr_list()| for the empty list. Often it will be practical
to use right recursive grammar rule that build lists backwards, so we provide
a reversal function to get the proper ordering once the end of the list is
reached. Finally we provide the wrapping function |wrap_expr_list| for list
displays.

@< Declarations of functions for the parser @>=
raw_expr_list make_exprlist_node(expr_p e, raw_expr_list l);
raw_expr_list reverse_expr_list(raw_expr_list l);
expr_p wrap_tuple_display(raw_expr_list l);
expr_p wrap_list_display(raw_expr_list l);

@~The function |reverse_expr_list| is easy as well. Of the two possible and
equivalent list reversal paradigms, we use the ``hold the head'' style which
starts with |t=l|; the other option is ``hold the tail'', which would start
with |t=l->next|. Either one would do just as well. We do not call
|reverse_expr_list| from |wrap_expr_list| although the two are usually
combined, since whether or not the list should be reversed can only be
understood when the grammar rules are given.

@< Definitions of functions for the parser @>=
raw_expr_list make_exprlist_node(expr_p e, raw_expr_list raw)
  { expr_ptr saf(e); expr_list l(raw);
    l.push_front(std::move(*e));
    return l.release();
  }

raw_expr_list reverse_expr_list(raw_expr_list raw)
{@; expr_list l(raw); l.reverse(); return l.release(); }
@)
expr_p wrap_tuple_display(raw_expr_list l)
{@; return new expr(expr_list(l),true); }
expr_p wrap_list_display(raw_expr_list l)
{@; return new expr(expr_list(l),false); }

@ Printing tuple displays and list displays is entirely similar, using
parentheses in the former case and brackets in the latter..

@< Cases for printing... @>=
case tuple_display:
{ const expr_list& l=e.sublist;
  if (l.empty()) out << "()";
  else
  { for (auto it=l.begin(); not l.at_end(it); ++it)
      out << (it==l.begin() ? '(' : ',') << *it;
    out << ')';
  }
}
break;
case list_display:
{ const expr_list& l=e.sublist;
  if (l.empty()) out << "[]";
  else
  { for (auto it=l.begin(); not l.at_end(it); ++it)
      out << (it==l.begin() ? '[' : ',') << *it;
    out << ']';
  }
}
break;

@ Length 0 tuple displays are sometimes used as a substitute expression where
nothing useful is provided, for instance for a missing else-branch. They are
easy to build, but for recognising them later, it is useful to have a function
at hand.

@< Declaration of functions not for the parser @>=
bool is_empty(const expr& e);

@~The implementation is of course straightforward.

@< Definitions of functions not for the parser @>=
bool is_empty(const expr& e)
@+{@; return e.kind==tuple_display and e.sublist.empty(); }

@*1 Function applications.
Another recursive type of expression is the function application.

@< Type declarations needed in definition of |struct expr@;| @>=
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

@
@< Cases for copying... @>=
  case function_call: call_variant = other.call_variant; break;

@ To print a function call, we look up the name of the function in the global
hash table, and either print a tuple display or a single expression enclosed
in parentheses for which we call the operator~`|<<|' recursively. We do not
attempt to reconstruct infix formulae.

@< Cases for printing... @>=
case function_call:
{ app a=e.call_variant;
  const expr& fun=a->fun; const expr& arg=a->arg;
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
  destroy_expr_body(e.call_variant->fun);
  destroy_expr_body(e.call_variant->arg);
  delete e.call_variant;
break;


@ To build an |application_node|, we combine the function identifier with an
|expr_list| for the argument. For the moment we do the packing or unpacking
(in case of a single argument) of the argument list here, rather than via the
syntax; the latter option would allow avoiding to pack singleton lists. But
the current method should be compatible with providing multiple arguments as a
single tuple value.

@< Declarations of functions for the parser @>=
expr_p make_application_node(expr_p f, raw_expr_list args);

@~Here for once there is some work to do. If a singleton argument list is
provided, the argument expression must be picked from it, but in all other
cases the argument list must be made into a tuple display. Note that it is
convenient here that |wrap_tuple_display| does not reverse the list, since
this is already done by the parser before calling |make_application_node|.

@< Definitions of functions for the parser @>=
expr mk_application_node(expr&& f, expr_list&& args)
{ app a=new application_node; a->fun=std::move(f);
  if (not args.empty() and args.at_end(++args.begin())) // a single argument
    a->arg=std::move(args.front());
  else
    a->arg=expr(std::move(args),true); // make into an argument tuple
  expr result; result.kind=function_call; result.call_variant=a;
  return result;
}
expr_p make_application_node(expr_p f, raw_expr_list args)
 { expr_ptr saf(f);
   return new expr(mk_application_node(std::move(*f),expr_list(args))); }

@ We shall frequently need to form a function application where the function
is accessed by an applied identifier and the argument is either a single
expression or a tuple of two expressions. We provide two functions to
facilitate those constructions.

@< Declarations of functions for the parser @>=
expr_p make_unary_call(id_type name, expr_p arg);
expr_p make_binary_call(id_type name, expr_p x, expr_p y);

@~In the unary case we avoid calling |make_application| with a singleton
list that will be immediately destroyed,

@< Definitions of functions for the parser @>=
expr mk_unary_call(id_type name, expr& arg)
{ app a=new application_node;
  a->fun=mk_applied_identifier(name);
  a->arg=std::move(arg);
  expr result; result.kind=function_call; result.call_variant=a;
  return result;
}
expr_p make_unary_call(id_type name, expr_p arg)
 { expr_ptr saf(arg); return new expr(mk_unary_call(name,*arg)); }
@)
expr mk_binary_call(id_type name, expr& x, expr& y)
{ expr_list args;
  args.push_front(std::move(y));
  args.push_front(std::move(x));
  return mk_application_node(mk_applied_identifier(name),std::move(args));
}
expr_p make_binary_call(id_type name, expr_p x, expr_p y)
 { expr_ptr xx(x), yy(y);
   return new expr(mk_binary_call(name,*x,*y));
 }

@*1 Operators and priority.
%
Applications of operators are transformed during parsing in function calls,
and the priorities of operators, which only serve to define the structure of
the resulting call tree, are a matter of grammar only. Nevertheless, handling
operator priorities purely in the syntax description has certain
disadvantages. It means that, unlike for identifiers which form just one token
category, the parser must distinguish separate tokens for each operator, or at
least for those of each priority level. Consequently it must also specify
separate syntax rules for each operator. Thus both the amount of tokens and
the amount of rules increase linearly with then number of operator symbols (or
priority levels) distinguished.
%; this might well imply a quadratic growth of the size of the parser tables

Moreover, although operator priority can perfectly well be described by syntax
rules alone, the repetitive nature of such a specification makes that parser
generators like \.{bison} propose as convenience the possibility of
supplementing an ambiguous set of syntax rules with a priority-based
disambiguation system. Since disambiguation is done at parser generation time,
the effect on the size of the parser tables is probably about the same as that
of specifying priorities explicitly in syntax rules (although that
necessitates additional non-terminal symbols for formulae at each priority
level, these do not directly enlarge the tables), while understanding that the
resulting parser always produces the desired parse trees requires analysing
the way disambiguation affects the given grammar.

We shall therefore opt for a solution in which handling of operator priorities
is performed neither via syntax rules nor via the disambiguation mechanism of
the parser generator, but rather via an explicit algorithm in the parser
actions. In other words we shall present a grammar whose reductions do not
follow the precise structure of parse tree that we want to construct, but
instead define parsing actions that restructure the tree dynamically to the
desired structure. Thus the grammar describes the language in a simplified way
(with all operators of equal priority and left associative), thus reducing the
amount of terminal symbols and parser states, while the parsing actions
explicitly perform the priority comparisons that would otherwise be encoded
implicitly in the transition tables of the parsing automaton.

Let us consider first the case of infix operators only. During lexical
analysis we shall associate an integral priority level to each operator
symbol; we encode in this level also the desired associativity of the operator
by the convention that associativity is to the left at even levels, and to the
right at odd levels. Thus there is never a conflict of different directions of
associativity at the same priority level.

In the absence of unary operators, a formula is an alternation of operands and
operators, where the former comprises any expression bound more tightly than
any formula, for instance one enclosed in parentheses. For every operator that
comes along we can already determine its left subtree by comparing priorities
with any previous ones. As a consequence of these comparisons, the right
subtree of some previous operators~$\omega$ may also be completed (if~$\omega$
turned out to bind more strongly than the current operator); in that case a
formula with root~$\omega$ is constructed and henceforth becomes a single
operand. Therefore at any point just after seeing an operator, there will be a
list of pending operators of increasing priorities (weakly at odd priority
levels), each with a complete left operand. If the formula would terminate
after one more operand, then each pending operator would get as right operand
the formula recursively constructed from everything that follows it. Now when
a new operand followed by an operator~$\omega$ comes along, the pending
operators are considered from right to left; while the operator has a higher
priority than~$\omega$, or the same even priority, it receives the new(est)
operand as its right operand, and the resulting formula becomes the newest
operand. Once an operator is encountered whose priority is too low to capture
the operand, the newest operand becomes the left subtree of~$\omega$, which
now becomes the leftmost pending operator.

Unary operators complicate the picture. Any priority one would like to
associate to a prefix operator can only be relevant with respect to operators
that follow, not those that precede. For instance one might want unary `\.-'
to have lower priority than `\.\pow' in order to parse \.{-x\pow2+y\pow2} as
$-(x^2)+(y^2)$, but \.{x\pow-2} can still only parse as $x^{(-2)}$. Giving
unary `\.-' a lower priority than `\.*' gives a real problem, since then
\.{-2*y} parses as $-(2*y)$, while certainly \.{x\pow2*y} parses as $(x^2)*y$,
so should \.{x\pow-2*y} parse as $(x^{-2})*y$ or as $x^{-(2*y)}$?

Two reasonable solutions exists for defining a general mechanism: the simple
solution is to give unary operators maximum priority (so that they can be
handled immediately), the other is to give that priority whenever they are
immediately preceded by another operator. We choose the latter option, since
it allows interpreting \.{-x\pow2} less surprisingly as $-(x^2)$; somewhat
more surprisingly the parentheses in $x+(-1)^n*y$ become superfluous. The
example \.{x\pow-2*y} will then parse as $(x^{-2})*y$, which also seems
reasonable. Note however that with this rule one does have the surprise
that \.{x\pow-y\pow2} parses as $x^{(-y)^2}$.

@ The data type necessary to store these intermediate data during priority
resolutions is a dynamic list of triples subtree-operator-priority. We use a
linked list, which can be used without difficulty from the parser; the link
points to an operator further to the left in the formula, whence it is called
|prev| rather than the more usual~|next|. To implement the above solution for
unary operators, we allow for the very first pending operator to not have any
left subtree; the expression is left of type |no_expr|, which can be tested to
detect the end of the list.

Postfix operators are quite rare in mathematics (the factorial exclamation
mark is the clearest example, though certain exponential notations like the
derivative prime could be considered as postfix operators as well) and more
importantly seem to invariably have infinite priority, so the can be handled
in the parser without dynamic priority comparisons. And even if such
comparisons were needed, they could be handled by a new function operating in
the list of partial formulae, and need not be taken into account in the data
structure of that list itself. So here is that structure:

@< Structure and typedef declarations... @>=
typedef struct partial_formula* form_stack;
struct partial_formula
{@; expr left_subtree; id_type op; int prio;
  form_stack prev;
};

@ We define the following functions operating on partial formulae: two to
start them out with a binary or unary operator, the principal one to extend
with a new operand and binary operator, one to finish off the formula with
a final operand, and of course one to clean up.

@< Declarations for the parser @>=
form_stack start_formula (expr_p e, id_type op, int prio);
form_stack start_unary_formula (id_type op, int prio);
form_stack extend_formula (form_stack pre, expr_p e,id_type op, int prio);
expr_p end_formula (form_stack pre, expr_p e);
void destroy_formula(form_stack s);

@ Starting a binary formula simply creates an initial node.
@< Definitions of functions for the parser @>=
form_stack start_formula (expr_p e, Hash_table::id_type op, int prio)
{ expr_ptr saf(e);
  form_stack result = new partial_formula;
  result->left_subtree=std::move(*e); result->op=op; result->prio=prio;
  result->prev=NULL;
  return result;
}
@)
form_stack start_unary_formula (id_type op, int prio)
{ form_stack result = new partial_formula;
  result->op=op; result->prio=prio; // leave |result->left_subtree| empty
  result->prev=NULL;
  return result;
}

@ Extending a formula involves the priority comparisons and manipulations
indicated above. It turns out |start_formula| could be replaced by a call to
|extend_formula| with |pre==NULL|.

@< Definitions of functions for the parser @>=
form_stack extend_formula (form_stack pre, expr_p ee,id_type op, int prio)
{ form_stack result=NULL; expr_ptr saf(ee); expr& e = *ee;
  while (pre!=NULL and (pre->prio>prio or pre->prio==prio and prio%2==0))
  { delete result; result=pre; // clean up, but holding on to one node
    @< Put |e| as right subtree into rightmost element of |pre|, and replace
    |e| by the result, popping it from |pre| @>
  }

  if (result==NULL)
    result=new partial_formula; // allocate if no nodes were combined
  result->left_subtree=std::move(e); result->op=op; result->prio=prio;
  result->prev=pre; return result;
}

@ Here we make either a unary or a binary operator call.

@< Put |e| as right subtree into rightmost element of |pre|...@>=
{ if (pre->left_subtree.kind==no_expr)
    e = mk_unary_call(pre->op,e); // apply initial unary operator
  else
    e = mk_binary_call(pre->op,pre->left_subtree,e);
  pre=pre->prev;
}

@ Wrapping up a formula is similar to the initial part of |extend_formula|,
but with an infinitely low value for the ``current priority'' |prio|. So we
can reuse the main part of the loop of |extend_formula|, with just a minor
modification to make sure all nodes get cleaned up after use.

@< Definitions of functions for the parser @>=
expr nd_formula (form_stack pre, expr& e)
{ while (pre!=NULL)
  { form_stack t=pre;
    @< Put |e| as right subtree into rightmost element of |pre|...@>
    delete t;
  }
  return std::move(e);
}
expr_p end_formula (form_stack pre, expr_p e)
{@; expr_ptr ee(e);
  return new expr(nd_formula(pre,*e));
}

@ Destroying a formula stack is straightforward.
@< Definitions of functions for the parser @>=
void destroy_formula(form_stack s)
{@; while (s!=NULL)
  {@; form_stack t=s; s=s->prev; delete t;
  }
}

@*1 Identifier patterns.
%
For each case where local identifiers will be introduced (like let-expressions
of function headings) we shall in fact allow more general patterns. A defining
occurrence of an identifier \\{ident} may be replaced by some tuple, say
$(x,y,z)$ (which assumes the value to be bound will be a $3$-tuple), or it
could be both, as in $(x,y,z)$:\\{ident} for the nested identifiers the same
options apply recursively, and in addition the identifier may be suppressed
altogether, to allow such partial tagging of components as \\{ident}:$(,,z)$.
To accommodate such possibilities we introduce the following recursive types.
Our grammar allows for a pattern $()$ (but not $()$:\\{ident}), in which case
|sublist.empty()|; it turned out that the possibility to bind no identifiers
at all (while providing a void value) has its uses. However patterns of the
form $(x)$, which would give a sublist of length~$1$, will be forbidden: they
would be confusing since $1$-tuples do not exist.

Since patterns pointed to by |patlist| will need to be duplicated (in the
evaluator), we provide a (default) constructor for |pattern_node| that will
ensure that constructed nodes will immediately be safe for possible
destruction by |destroy_id_pat| (no undefined pointers; the |sublist| pointer
of |id_pat| is ignored when |kind==0|). The structure |id_pat| itself cannot
have a constructor, since it figures in a |union|, where this is not allowed.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct pattern_node* patlist;
struct id_pat
{ patlist sublist;
  id_type name;
  unsigned char kind; /* bit 0: has name, bit 1: has sublist */
};
struct pattern_node
{ patlist next;
  struct id_pat body;
  pattern_node () : next(NULL) @+{@; body.kind=0; }
};

@ These types do not themselves represent a variant of |expr|, but will be
used inside such variants. We can already provide a printing function.

@< Declaration of functions not for the parser @>=
std::ostream& operator<< (std::ostream& out, const id_pat& p);

@~Only parts whose presence is indicated in |kind| are printed. We take care
that even if |sublist| is marked as present, it might be null. However, the
list cannot be of length~1, so having put aside the case of an empty sublist,
we can print a first comma unconditionally.

@< Definitions of functions not for the parser @>=
std::ostream& operator<< (std::ostream& out, const id_pat& p)
{ if ((p.kind & 0x2)!=0)
    if (p.sublist==NULL)
      out << "()";
    else
    { out << '(' << p.sublist->body << ',';
      for (patlist l=p.sublist->next; l!=NULL; l=l->next)
      out << l->body << (l->next!=NULL ? ',' : ')');
    }
  if (p.kind==0x3) // both parts present
    out << ':';
  if ((p.kind & 0x1)!=0)
    out << main_hash_table->name_of(p.name);
  return out;
}

@ The function to build a node takes a pointer to a structure with the
contents of the current node; in practice this is the address of a local
variable in the parser. Patterns also need cleaning up, which
|destroy_pattern| and |destroy_id_pat| will handle, and reversal as handled by
|reverse_patlist|.

@< Declarations for the parser @>=
patlist make_pattern_node(patlist next,struct id_pat* body);
void destroy_pattern(patlist p);
void destroy_id_pat(struct id_pat* p);
patlist reverse_patlist(patlist p);

@ The function just assembles the pieces. In practice the |next| pointer will
point to previously parsed nodes, so (as usual) reversal will be necessary.
Being bored, we add a variation on list reversal.

@< Definitions of functions for the parser @>=
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
first step towards having user defined functions. Indeed let-expressions will
be implemented (initially) as a user-defined functions that is immediately
called with the values bound in the let-expression. The reason to do this
before considereing user-defined functions in general is that this avoids for
now having to specify in the user program the types of the parameters of the
function, these type being determined by the values provided.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct let_expr_node* let;

@~After parsing, let-expression will have a single let-binding followed by a
body giving the value to be returned. During parsing however, we may form
intermediate values containing list of let-bindings, that will later be
converted into a single one with a tuple as left hand side. We therefore
define a node type |let_node| for a list of bindings, and a structure
|let_expr_node| for a complete let-expression, containing only one binding,
but in addition a body.

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

@
@< Cases for copying... @>=
  case let_expr: let_variant = other.let_variant; break;

@ To print a let-expression we do the obvious things, wihtout worrying about
parentheses; this should be fixed (for all printing routines).

@< Cases for printing... @>=
case let_expr:
{ let lexp=e.let_variant;
  out << "let " << lexp->pattern << '=' << lexp->val <<	 " in " << lexp->body;
}
break;

@ Destroying lists of declarations will be done in a function callable from the
parser, like |destroy_exprlist|.

@< Declarations of functions for the parser @>=
void destroy_letlist(let_list l);

@~Like |destroy_exprlist|, this function recursively destroys nodes, and the
expressions they contain.

@< Definitions of functions for the parser @>=
void destroy_letlist(let_list l)
{ while (l!=NULL)
    @/{@; let_list p=l; l=l->next;
      destroy_id_pat(&p->pattern);
      destroy_expr_body(p->val);
      delete p;
    }
}

@~Here we clean up the declaration, and then the body of the let-expression.

@< Cases for destroying... @>=
case let_expr:
{ let lexp=e.let_variant;
  destroy_id_pat(&lexp->pattern);
  destroy_expr_body(lexp->val);
  destroy_expr_body(lexp->body);
  delete lexp;
}
break;

@ For building let-expressions, three functions will be defined. The function
|make_let_node| makes a list of one declaration, while |append_let_node|
appends such a list |cur| (assured to be of length~$1$) to a previously
constructed list |prev| of declaration; finally |make_let_expr_node| wraps up
an entire let-expression.

@< Declarations of functions for the parser @>=
let_list make_let_node(struct id_pat pattern, expr_p val);
let_list append_let_node(let_list prev, let_list cur);
expr_p make_let_expr_node(let_list decls, expr_p body);

@~The functions |make_let_node| and |append_let_node| build a list in reverse
order, which makes the latter function a particularly simple one. In
|make_let_expr_node| combines multiple declarations (that were separated by
commas) to one, taking care to reverse the order at the same time so that the
tuple of patterns being declared comes out in the same order as it was
specified in the program. Fortunately it is actually easier to build a merged
list in reverse order.

We provide local exception safety here, in the sense that if an exception is
thrown we clean up the locally created nodes. As a side remark, we do realise
that, as the parser function is written in \Cee, it will not be able to do any
cleaning up of values referred to in its stack, in case of an exception thrown
during its execution (here, or in calls to make-functions elsewhere, which it
should be noted is fairly unlikely to happen). In fact the problem is not just
one of programming language: even a parser generated as \Cpp~code would not
help without special provisions (such as an exception handler in the parser
function), since the |union| used for values on the parsing stack cannot have
as members objects that handle their own cleaning up. Rather, this situation
should probably be corrected by using a different allocation strategy while
building the parse tree (avoiding explicit calls to |new| but allocating from
local storage pools that can explicitly be emptied), in which case all
cleaning up in the code below should also be removed.

@< Definitions of functions for the parser @>=
let_list make_let_node(struct id_pat pattern, expr_p val)
{ expr_ptr saf(val);
  let_list l=new let_node;
  l->pattern=pattern; l->val=std::move(*val); l->next=NULL;
  return l;
}
@)
let_list append_let_node(let_list prev, let_list cur)
{@; cur->next=prev; return cur; }
@)
expr mk_let_expr_node(let_list decls, expr& body)
{ let l=new let_expr_node; l->body=std::move(body);
  expr result; result.kind=let_expr; result.let_variant=l;
  try
  { if (decls->next==NULL) // single declaration
    @/{@; l->pattern=decls->pattern; l->val=std::move(decls->val);
      delete decls;
    }
    else
    { l->pattern.kind=0x2; l->pattern.sublist=NULL;
      l->val=expr(expr_list(),true);
      while (decls!=NULL)
      { l->pattern.sublist =
	  make_pattern_node(l->pattern.sublist,&decls->pattern);
	l->val.sublist.push_front(std::move(decls->val));
	let_list p=decls; decls=p->next; delete p;
      }
    }
    return result;
  }
  catch(...)
  {@; destroy_expr_body(result);
    throw;
  }
}
expr_p make_let_expr_node(let_list decls, expr_p body)
 { expr_ptr saf(body); return new expr(mk_let_expr_node(decls,*body)); }

@*1 Types and user-defined functions.
%
One reason let-expressions were introduced before user-defined functions, is
that it avoids to problem of having to specify types in the user program.
Inevitably we have to deal with that though, since for a function definition
there is no way to know the type of the arguments with certainty, unless the
user specifies them (at least this is true if we want to allow such things as
type coercion and function overloading to be possible, so that types cannot be
deduced by analysis of the \emph{usage} of the parameters in the function
only). Types are also an essential ingredient of casts.

When the parser was compiled as \Cee~code, we were forced to use void pointers
in the parser, to masquerade for the actual pointers to \Cpp~types; numerous
static casts were then used to get the proper pointer types from them. Now
that this is no longer necessary, everything has been reformulated in terms of
the actual pointer types. Something that remains (for now) is the avoidance of
smart pointers for types in the parser, since other pointers for expressions
that it handles are not smart pointers either.

We avoid including \.{types.h} into our header file \.{parsetree.h}, since the
reverse inclusion is present and necessary. But we do need to know about the
following type names, where it fortunately suffices to know they are pointers
to unspecified structures.

@< Structure and typedef... @>=
struct type_expr;
typedef class atlas::containers::sl_node<type_expr>* raw_type_list; // predeclare;
typedef class atlas::containers::simple_list<type_expr> type_list; // predeclare;
typedef type_expr* type_p;
typedef std::unique_ptr<type_expr> type_ptr;

@ These functions provide an interface to routines defined in the
module \.{types.w}, stripping off the smart pointers.

@< Declarations of functions for the parser @>=
raw_type_list mk_type_singleton(type_p raw);
raw_type_list mk_type_list(type_p t,raw_type_list l);
type_p mk_prim_type(int p);
type_p mk_row_type(type_p c);
type_p mk_tuple_type(raw_type_list l);
type_p mk_function_type(type_p a,type_p r);
@)
void destroy_type(type_p t);
void destroy_type_list(raw_type_list t);

@ The following function is not used in the parser, and never had \Cee~linkage.

@< Declaration of functions not for the parser @>=
std::ostream& print_type(std::ostream& out, type_p type);

@ All that is needed are conversions from ordinary pointer to auto-pointer and
back (or from integer to enumeration type). Nothing can throw during these
conversions, so passing bare pointers is exception-safe.

@< Definitions of functions for the parser @>=

raw_type_list mk_type_singleton(type_p raw)
{ type_ptr t(raw); // ensures node is cleaned up
  type_list result;
  result.push_front(std::move(*t));
  return result.release();
}

raw_type_list mk_type_list(type_p raw,raw_type_list l)
{ type_ptr t(raw); // ensure clean-up
  type_list tmp(l); // since |prefix| needs second argument an lvalue reference
  return prefix(std::move(*t),tmp).release();
}
@)
type_p mk_prim_type(int p)
{@; return make_prim_type(static_cast<primitive_tag>(p)).release(); }

type_p mk_row_type(type_p c)
{@; return make_row_type(type_ptr(c)).release(); }

type_p mk_tuple_type(raw_type_list l)
{@; return make_tuple_type(type_list(l)).release(); }

type_p mk_function_type(type_p a,type_p r)
{ type_ptr pa(a), pr(r); // ensure cleaning up at return
  return make_function_type(std::move(*a),std::move(*r)).release(); }
@)

void destroy_type(type_p t)@+ {@; delete t; }
void destroy_type_list(raw_type_list t)@+ {@; delete t; } // recursive destruction

std::ostream& print_type(std::ostream& out, type_p type) @+
{@; return out << *type; }

@ For user-defined functions we shall use a structure |lambda_node|.
@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct lambda_node* lambda;

@~It contains a pattern for the formal parameter(s), its type (a void pointer
that actually points to a |type_declarator| structure defined
in \.{evaluator.w}), and an  expression (the body of the function).

@< Structure and typedef... @>=
struct lambda_node
{@; struct id_pat pattern; type_p arg_type; expr body; };

@ The tag used for user-defined functions is |lambda_expr|.
@< Enumeration tags... @>=
lambda_expr,@[@]

@ We introduce the variant of |expr| as usual.
@< Variants of |union... @>=
lambda lambda_variant;

@
@< Cases for copying... @>=
  case lambda_expr: lambda_variant = other.lambda_variant; break;

@ We must take care of printing lambda expressions; we avoid a double set of
parentheses.

@< Cases for printing... @>=
case lambda_expr:
{ lambda fun=e.lambda_variant;
  if ((fun->pattern.kind&0x1)!=0)
    out << '(' << fun->pattern << ')';
  else
    out << fun->pattern;
  out << ':' << fun->body;
}
break;

@ And we must of course take care of destroying lambda expressions, which just
call handler functions.

@h "types.h" // so that |type_p| will point to a complete type in |delete|

@< Cases for destroying an expression |e| @>=
case lambda_expr:
{ lambda fun=e.lambda_variant;
  destroy_id_pat(&fun->pattern);
  delete(fun->arg_type);
  destroy_expr_body(fun->body);
}
break;

@ Finally there is as usual a function for constructing a node, to be called
by the parser.

@< Declarations of functions for the parser @>=
expr_p make_lambda_node(patlist pat_l, raw_type_list type_l, expr_p body);

@~There is a twist in building a lambda node, in that for syntactic reasons
the parser passes lists of patterns and types rather than single ones. We must
distinguish the case of a singleton, in which case the head node must be
unpacked, and the multiple case, where a tuple pattern and type must be
wrapped up from the lists. In the former case, |fun->arg_type| wants to have a
pointer to an isolated |type_expr|, but the head of |type_l| is a |type_node|
that contains a |type_expr| as its |t| field; making a (shallow) copy of that
field is the easiest way to obtain an isolated |type_expr|. After the
copy, destruction of |type_l| deletes the original |type_node|.

@< Definitions of functions for the parser @>=
expr mk_lambda_node(patlist pat_l, raw_type_list raw, expr& body)
{ type_list type_l(raw); // smart pointer might clean up a node
  lambda fun=new lambda_node; fun->body=std::move(body);
  if (not type_l.empty() and type_l.at_end(++type_l.begin()))
  { fun->pattern=pat_l->body; delete pat_l; // clean up node
    fun->arg_type = new type_expr(std::move(type_l.front()));
     // make a shallow copy
  }
  else
  { fun->pattern.kind=0x2; fun->pattern.sublist=pat_l;
    fun->arg_type=make_tuple_type(std::move(type_l)).release();
  }
  expr result; result.kind=lambda_expr; result.lambda_variant=fun;
  return result;
}
expr_p make_lambda_node(patlist pat_l, raw_type_list raw, expr_p body)
 { expr_ptr saf(body); return new expr(mk_lambda_node(pat_l,raw,*body)); }

@*1 Control structures.

@*2 Conditional expressions.
Of course we need if-then-else expressions.

@< Type declarations needed in definition of |struct expr@;| @>=
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

@
@< Cases for copying... @>=
  case conditional_expr: if_variant = other.if_variant; break;

@ To print a conditional expression at parser level, we shall not
reconstruct \&{elif} constructions.

@< Cases for printing... @>=
case conditional_expr:
{ cond c=e.if_variant; out << " if " << c->condition @|
  << " then " << c->then_branch << " else " << c->else_branch << " fi ";
}
break;

@~Here we clean up constituent expressions and then the node for the conditional
call itself.

@< Cases for destroying... @>=
case conditional_expr:
  destroy_expr_body(e.if_variant->condition);
  destroy_expr_body(e.if_variant->then_branch);
  destroy_expr_body(e.if_variant->else_branch);
  delete e.if_variant;
break;


@ To build an |conditional_node|, we define a function as usual.
@< Declarations of functions for the parser @>=
expr_p make_conditional_node(expr_p c, expr_p t, expr_p e);

@~It is entirely straightforward.

@< Definitions of functions for the parser @>=
expr mk_conditional_node(expr& c, expr& t, expr& e)
{ cond n=new conditional_node; n->condition=std::move(c);
  n->then_branch=std::move(t); n->else_branch=std::move(e);
@/expr result; result.kind=conditional_expr; result.if_variant=n;
  return result;
}
expr_p make_conditional_node(expr_p c, expr_p t, expr_p e)
 { expr_ptr saf0(c), saf1(t), saf2(e);
   return new expr(mk_conditional_node(*c,*t,*e));
 }

@*2 Loops.
%
Loops are a cornerstone of any form of non-recursive programming. The two main
flavours provided are traditional: |while| and |for| loops; the former provide
open-ended iteration, while the latter fix the range of iteration at entry.
However we provide a relative innovation by having both types deliver a value:
this will be a row value with one entry for each iteration performed.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct while_node* w_loop;
typedef struct for_node* f_loop;
typedef struct cfor_node* c_loop;

@~A |while| loop has two elements: a condition (which determines whether an
iteration will be undertaken), and a body (which contributes an entry to the
result). A |for| loop has three parts, a pattern introducing variables, an
expression iterated over (the in-part) and the loop body. A counted |for| loop
(the simple version of |for| loop) has four parts (an identifier, a count, a
lower bound and a body) and two variants (increasing and decreasing) which can
be distinguished by a boolean.

@< Structure and typedef declarations for types built upon |expr| @>=
struct while_node
 {@; expr condition; expr body; };
struct for_node
 {@; struct id_pat id; expr in_part; expr body; };
struct cfor_node
 {@; id_type id; expr count; expr bound; short up; expr body; };

@ The tags used for these expressions are |while_expr|, |for_expr| and
|cfor_expr|.

@< Enumeration tags for |expr_kind| @>= while_expr, for_expr, cfor_expr, @[@]

@~The variant of |expr| values with an |w_loop| as parsing value is tagged
|while_variant|.

@< Variants... @>=
w_loop while_variant;
f_loop for_variant;
c_loop cfor_variant;

@
@< Cases for copying... @>=
  case while_expr: while_variant = other.while_variant; break;
  case for_expr: for_variant = other.for_variant; break;
  case cfor_expr: cfor_variant = other.cfor_variant; break;

@ To print a |while| or |for| expression at parser level, we reproduce the
input syntax.

@< Cases for printing... @>=
case while_expr:
{ w_loop w=e.while_variant;
  out << " while " << w->condition << " do " << w->body << " od ";
}
break;
case for_expr:
{ f_loop f=e.for_variant;
  out << " for " << f->id.sublist->next->body;
  if (f->id.sublist->body.kind==0x1)
    out << '@@' << f->id.sublist->body;
  out << " in " << f->in_part << " do " << f->body << " od ";
}
break;
case cfor_expr:
{ c_loop c=e.cfor_variant;
  out << " for " << main_hash_table->name_of(c->id) << ": " << c->count;
  if (c->bound.kind!=tuple_display or c->bound.sublist.empty())
    out << (c->up!=0 ? " from " : " downto ") << c->bound;
  out << " do " << c->body << " od ";
}
break;

@~Cleaning up is exactly like that for conditional expressions.

@< Cases for destroying... @>=
case while_expr:
  destroy_expr_body(e.while_variant->condition);
  destroy_expr_body(e.while_variant->body);
  delete e.while_variant;
break;
case for_expr:
  destroy_id_pat(&e.for_variant->id);
  destroy_expr_body(e.for_variant->in_part);
  destroy_expr_body(e.for_variant->body);
  delete e.for_variant;
break;
case cfor_expr:
  destroy_expr_body(e.cfor_variant->count);
  destroy_expr_body(e.cfor_variant->bound);
  destroy_expr_body(e.cfor_variant->body);
  delete e.cfor_variant;
break;


@ To build a |while_node|, |for_node| or |cfor_node|, here are yet three
more \\{make}-functions.

@< Declarations of functions for the parser @>=
expr_p make_while_node(expr_p c, expr_p b);
expr_p make_for_node(struct id_pat id, expr_p ip, expr_p b);
expr_p make_cfor_node(id_type id, expr_p count, expr_p bound, short up, expr_p b);

@~They are quite straightforward, as usual.

@< Definitions of functions for the parser @>=
expr mk_while_node(expr& c, expr& b)
{ w_loop w=new while_node; w->condition=std::move(c); w->body=std::move(b);
@/expr result; result.kind=while_expr; result.while_variant=w;
  return result;
}
expr_p make_while_node(expr_p c, expr_p b)
 { expr_ptr saf0(c),saf1(b); return new expr(mk_while_node(*c,*b)); }

expr mk_for_node(struct id_pat id, expr& ip, expr& b)
{ f_loop f=new for_node;
  f->id=id; f->in_part=std::move(ip); f->body=std::move(b);
@/expr result; result.kind=for_expr; result.for_variant=f;
  return result;
}
expr_p make_for_node(struct id_pat id, expr_p ip, expr_p b)
 { expr_ptr saf0(ip),saf1(b); return new expr(mk_for_node(id,*ip,*b)); }

expr mk_cfor_node(id_type id, expr& count, expr& bound, short up, expr& b)
{ c_loop c=new cfor_node;
  c->id=id; c->count=std::move(count); c->bound=std::move(bound);
  c->up=up; c->body=std::move(b);
@/expr result; result.kind=cfor_expr; result.cfor_variant=c;
  return result;
}
expr_p make_cfor_node(id_type id, expr_p count, expr_p bound, short up, expr_p b)
 { expr_ptr saf0(count),saf1(bound),saf2(b);
   return new expr(mk_cfor_node(id,*count,*bound,up,*b));
 }

@*1 Array subscriptions.
%
We want to be able to select components from array structures (lists, vectors,
matrices), so we define a subscription expression. If there are multiple
indices, these can be realised as a subscription by a tuple expression, so we
define only one type so subscription expression.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct subscription_node* sub;

@~In a subscription the array and the index(es) can syntactically be arbitrary
expressions (although the latter should have as type integer, or a tuple of
integers).

@< Structure and typedef declarations for types built upon |expr| @>=
struct subscription_node {@; expr array; expr index; };

@ Here is the tag used for subscriptions.

@< Enumeration tags for |expr_kind| @>= subscription, @[@]

@~And here is the corresponding variant of |expr|.

@< Variants... @>=
sub subscription_variant;

@
@< Cases for copying... @>=
  case subscription: subscription_variant = other.subscription_variant; break;

@ To print a subscription, we just print the expression of the array, followed
by the expression for the index in brackets. As an exception, the case of a
tuple display as index is handled separately, in order to avoid having
parentheses directly inside the brackets.

@h "lexer.h"
@< Cases for printing... @>=
case subscription:
{ sub s=e.subscription_variant; out << s->array << '[';
  const expr& i=s->index;
  if (i.kind!=tuple_display) out << i;
  else
    for (auto it=i.sublist.begin(); not i.sublist.at_end(it); ++it)
      out << (it==i.sublist.begin() ? "" : ",") << *it;
  out << ']';
}
break;

@~Here we recursively destroy both subexpressions, and then the node for the
subscription call itself.

@< Cases for destroying... @>=
case subscription:
  destroy_expr_body(e.subscription_variant->array);
  destroy_expr_body(e.subscription_variant->index);
  delete e.subscription_variant;
break;


@ To build an |subscription_node|, we simply combine the array and the index
part.

@< Declarations of functions for the parser @>=
expr_p make_subscription_node(expr_p a, expr_p i);

@~This is straightforward, as usual.

@< Definitions of functions for the parser @>=
expr mk_subscription_node(expr& a, expr& i)
{ sub s=new subscription_node; s->array=std::move(a); s->index=std::move(i);
  expr result; result.kind=subscription; result.subscription_variant=s;
  return result;
}
expr_p make_subscription_node(expr_p a, expr_p i)
 { expr_ptr saf0(a),saf1(i); return new expr(mk_subscription_node(*a,*i)); }

@*1 Cast expressions.
%S
These are very simple expressions consisting of a type and an expression,
which is forced to be of that type.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct cast_node* cast;

@~The type contained in the cast is represented by a void pointer, for the
known reasons.

@< Structure and typedef declarations for types built upon |expr| @>=
struct cast_node {@; type_p type; expr exp; };

@ The tag used for casts is |cast_expr|.

@< Enumeration tags for |expr_kind| @>=
cast_expr, @[@]

@ And there is of course a variant of |expr_union| for casts.
@< Variants of ... @>=
cast cast_variant;

@
@< Cases for copying... @>=
  case cast_expr: cast_variant = other.cast_variant; break;

@ Printing cast expressions follows their input syntax.

@< Cases for printing... @>=
case cast_expr:
{@; cast c = e.cast_variant;
  print_type(out,c->type) << ':' << c->exp ;
}
break;

@ Casts are built by |make_cast|.

@< Declarations of functions for the parser @>=
expr_p make_cast(type_p type, expr_p exp);

@~No surprises here.

@< Definitions of functions for the parser@>=
expr mk_cast(type_p type, expr& exp)
{ cast c=new cast_node; c->type=type; c->exp=std::move(exp);
@/ expr result; result.kind=cast_expr; result.cast_variant=c;
   return result;
}
expr_p make_cast(type_p type, expr_p exp)
 { expr_ptr saf(exp); return new expr(mk_cast(type,*exp)); }

@ Eventually we want to rid ourselves from the cast.

@< Cases for destr... @>=
case cast_expr:
  destroy_expr_body(e.cast_variant->exp); delete e.cast_variant;
break;

@ A different kind of cast serves to obtain the current value of an overloaded
operator symbol.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct op_cast_node* op_cast;

@~We store an (operator) identifier and a type, as before represented by a
void pointer.

@< Structure and typedef declarations for types built upon |expr| @>=
struct op_cast_node {@; id_type oper; type_p type; };

@ The tag used for casts is |op_cast_expr|.

@< Enumeration tags for |expr_kind| @>=
op_cast_expr, @[@]

@ And there is of course a variant of |expr_union| for casts.
@< Variants of ... @>=
op_cast op_cast_variant;

@
@< Cases for copying... @>=
  case op_cast_expr: op_cast_variant = other.op_cast_variant; break;

@ Printing operator cast expressions follows their input syntax.

@< Cases for printing... @>=
case op_cast_expr:
{ op_cast c = e.op_cast_variant;
  print_type(out << main_hash_table->name_of(c->oper) << '@@',c->type);
}
break;

@ Casts are built by |make_cast|.

@< Declarations of functions for the parser @>=
expr_p make_op_cast(id_type name,type_p type);

@~No surprises here either.

@< Definitions of functions for the parser@>=
expr mk_op_cast(id_type name,type_p type)
{ op_cast c=new op_cast_node; c->oper=name; c->type=type;
@/ expr result; result.kind=op_cast_expr; result.op_cast_variant=c;
   return result;
}
expr_p make_op_cast(id_type name,type_p type)
 { return new expr(mk_op_cast(name,type)); }

@ Eventually we want to rid ourselves from the operator cast.

@< Cases for destr... @>=
case op_cast_expr:
  delete e.op_cast_variant;
break;

@*1 Assignment statements.
%
Simple assignment statements are quite simple as expressions.

@< Type declarations needed in definition of |struct expr@;| @>=
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

@
@< Cases for copying... @>=
  case ass_stat: assign_variant = other.assign_variant; break;

@ Printing assignment statements is absolutely straightforward.

@< Cases for printing... @>=
case ass_stat:
{@; assignment ass = e.assign_variant;
  out << main_hash_table->name_of(ass->lhs) << ":=" << ass->rhs ;
}
break;

@ Assignment statements are built by |make_assignment|.

@< Declarations of functions for the parser @>=
expr_p make_assignment(id_type lhs, expr_p rhs);

@~It does what one would expect it to (except for those who expect their
homework assignment made).

@< Definitions of functions for the parser@>=
expr mk_assignment(id_type lhs, expr& rhs)
{ assignment a=new assignment_node; a->lhs=lhs; a->rhs=std::move(rhs);
@/ expr result; result.kind=ass_stat; result.assign_variant=a;
   return result;
}
expr_p make_assignment(id_type lhs, expr_p rhs)
 { expr_ptr saf(rhs); return new expr(mk_assignment(lhs,*rhs)); }

@ What is made must eventually be unmade (even assignments).

@< Cases for destr... @>=
case ass_stat:
  destroy_expr_body(e.assign_variant->rhs); delete e.assign_variant;
break;

@*2 Component assignments.
%
We have special expressions for assignments to a component.

@< Type declarations needed in definition of |struct expr@;| @>=
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

@
@< Cases for copying... @>=
  case comp_ass_stat: comp_assign_variant = other.comp_assign_variant; break;

@ Printing component assignment statements follow the input syntax.

@< Cases for printing... @>=
case comp_ass_stat:
{@; comp_assignment ass = e.comp_assign_variant;
  out << main_hash_table->name_of(ass->aggr) << '[' << ass->index << "]:="
      << ass->rhs ;
}
break;

@ Assignment statements are built by |make_assignment|, which for once does
not simply combine the expression components, because for reason of parser
generation the array and index will have already been combined before this
function can be called.

@< Declarations of functions for the parser @>=
expr_p make_comp_ass(expr_p lhs, expr_p rhs);

@~Here we have to take the left hand side apart a bit, and clean up its node.

@< Definitions of functions for the parser@>=
expr mk_comp_ass(expr& lhs, expr& rhs)
{ comp_assignment a=new comp_assignment_node;
@/a->aggr=lhs.subscription_variant->array.identifier_variant;
  a->index=std::move(lhs.subscription_variant->index);
  a->rhs=std::move(rhs);
@/expr result; result.kind=comp_ass_stat; result.comp_assign_variant=a;
  return result;
}
expr_p make_comp_ass(expr_p lhs, expr_p rhs)
 { expr_ptr saf0(lhs),saf1(rhs); return new expr(mk_comp_ass(*lhs,*rhs)); }

@~Destruction one the other hand is as straightforward as usual.

@< Cases for destr... @>=
case comp_ass_stat:
  destroy_expr_body(e.comp_assign_variant->index);
  destroy_expr_body(e.comp_assign_variant->rhs);
  delete e.comp_assign_variant;
break;

@*1 Sequence statements.
%
Having assignments statements, it is logical to be able to build a sequence of
expressions (statements) as well, retaining the value only of the final one.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct sequence_node* sequence;

@~Since control structures and let-expressions tend to break up long chains,
we do not expect their average length to be very great. So we build up
sequences by chaining pairs, instead of storing a vector of expressions: at
three pointers overhead per vector, a chained representation is more compact
as long as the average length of a sequence is less than~$5$ expressions.

In fact we used this expression also for a variant of sequence expressions, in
which the value of the \emph{first} expression is retained as final value,
while the second expression is then evaluated without using its value; the
|forward| field indicates whether the first form was used.

@< Structure and typedef declarations for types built upon |expr| @>=
struct sequence_node {@; expr first; expr last; int forward; };

@ The tag used for sequence statements is |seq_expr|.

@< Enumeration tags for |expr_kind| @>=
seq_expr, @[@]

@ And there is of course a variant of |expr_union| for sequences.
@< Variants of ... @>=
sequence sequence_variant;

@
@< Cases for copying... @>=
  case seq_expr: sequence_variant = other.sequence_variant; break;

@ Printing sequences is absolutely straightforward.

@< Cases for printing... @>=
case seq_expr:
{@; sequence seq = e.sequence_variant;
  out << seq->first << (seq->forward ? ";" : " next ") << seq->last ;
}
break;

@ Sequences are built by |make_sequence|.

@< Declarations of functions for the parser @>=
expr_p make_sequence(expr_p first, expr_p last, int forward);

@~It does what one would expect it to.

@< Definitions of functions for the parser @>=
expr mk_sequence(expr& first, expr& last, int forward)
{ sequence s=new sequence_node;
  s->first=std::move(first); s->last=std::move(last); s->forward=forward;
@/ expr result; result.kind=seq_expr; result.sequence_variant=s;
   return result;
}
expr_p make_sequence(expr_p first, expr_p last, int forward)
 { expr_ptr saf0(first),saf1(last);
   return new expr(mk_sequence(*first,*last,forward));
 }

@ Finally sequence nodes need destruction, like everything else.

@< Cases for destr... @>=
case seq_expr:
  destroy_expr_body(e.sequence_variant->first);
  destroy_expr_body(e.sequence_variant->last);
  delete e.sequence_variant;
break;

@* Other functions callable from the parser.
Here are some functions that are not so much a parsing functions as just
wrapper functions enabling the parser to call \Cpp~functions.

@< Declarations of functions for the parser @>=
id_type lookup_identifier(const char*);
void include_file(int skip_seen);

@~The parser will only call this with string constants, so we can use the
|match_literal| method.

@< Definitions of functions for the parser @>=
id_type lookup_identifier(const char* name)
{@; return main_hash_table->match_literal(name); }

@~To include a file, we call the |push_file| method from the input buffer,
providing a file name that was remembered by the lexical analyser. If this
fails, then we abort all includes, as there is not much point in continuing to
read a file when another on which it depends cannot be found.

@< Definitions of functions for the parser @>=
void include_file(int skip_seen)
{ if (not main_input_buffer->push_file
          (lex->scanned_file_name(),skip_seen!=0))
    main_input_buffer->close_includes();
     // nested include failure aborts all includes
}

@* Index.

% Local IspellDict: british
