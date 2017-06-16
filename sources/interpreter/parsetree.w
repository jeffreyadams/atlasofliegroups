% Copyright (C) 2006-2015 Marc van Leeuwen
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
\chardef\pow = `\^

@* Building the parse tree.
This is the program unit \.{parsetree} which produces the implementation file
\.{parsetree.cpp} and two header files \.{parse\_types.h} and \.{parsetree.h}.
The file \.{parse\_types.h} is read, using \&{\#include}, early in the
file \.{parser.tab.c} generated from \.{parser.y}, so that the type |YYSTYPE|
used on the parsing stack can be properly defined, while the
file \.{parsetree.h} containing declarations of the functions defined here is
read later. This separation is done mainly so that we can without circularity
include \.{parser.tab.h} into \.{parsetree.h}, which requires
including \.{parse\_types.h} first. We do need to pre-declare the type
|YYLTYPE| that will be defined in \.{parser.tab.h}, but which is already
referred to in the interfaces of classes defined in \.{parse\_types.h};
moreover this has to be done outside any namespaces.

@( parse_types.h @>=
#ifndef PARSE_TYPES_H
#define PARSE_TYPES_H
@< Includes needed in \.{parse\_types.h} @>@;
struct YYLTYPE; // pre-declare

namespace atlas
{
  namespace interpreter
  {
@< Type declarations for the parser @>@;
  }@;
}@;
#endif


@ It used to be the case that this file was compiled by a \Cee\ compiler, and
therefore produced functions with \Cee-language linking. As a consequence some
jumping through hoops was necessary to cleanly integrate them into the \Cpp\
program. However it turns out to be possible to compile the
file \.{parser.tab.c} by a \Cpp-compiler, which makes all linkage to be
for \Cpp, and removes any need for |extern "C"| declarations. For the moment
these declarations are simply removed, and the structure of this file still
carries a legacy of the original design.

This file also defines some other functions defined here are not
used in the parser, but can be used by other modules written in~\Cpp; their
declaration is separated for historic reasons only.

@( parsetree.h @>=
#ifndef PARSETREE_H
#define PARSETREE_H

#include <iostream>
#include <memory>
#include <string>

#include "parse_types.h" // this must precede the next include
#include "parser.tab.h" // because \.{parser.tab.h} does not look after itself
namespace atlas
{
  namespace interpreter
  {
@< Declarations of functions for the parser @>@;

@< Declarations of functions not for the parser @>@;
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

@*1 Outline of the type \ {\bf expr}.
%
For a large part the declarations for the parser consist of the recursive
definition of the type |expr|. While that used to be a POD type used directly
on the parser stack, this solution was very inflexible; it was therefore
replaced by one where |expr| is not so constrained, and raw pointers |expr_p|
to it are what is placed on the parser stack. Since the parser rarely needs to
take apart parsing values, which would require dereferencing the pointer, it
might be that just declaring |typedef struct expr* expr_p;| here would have
sufficed, but we leave the detailed type definitions visible to the parser
anyway.

@< Type declarations for the parser @>=
@< Type declarations needed in definition of |struct expr@;| @>@;

enum expr_kind @+
 { @< Enumeration tags for |expr_kind| @>@;@; @+no_expr };
typedef struct expr* expr_p; // raw pointer type for use on parser stack
struct expr {
  expr_kind kind;
  union {@; @< Variants of the anonymous |union| in |expr| @>@; };
  source_location loc;
@)
  @< Methods of |expr| @>@;
};
typedef std::unique_ptr<expr> expr_ptr;
@)
@< Structure and typedef declarations for types built upon |expr| @>@;

@ We start right away declaring |expr| as a |struct|, avoiding complaints that
it is not declared.

@< Type declarations needed in definition of |struct expr@;| @>=
struct expr;

@ To represent identifiers efficiently, and also file names, we shall use the
type |Hash_table::id_type| (a small integer type) of indices into the table of
identifier or of files names, which we lift out of that class by using a
|typedef|.

@< Includes needed... @>=
#include "buffer.h" // for |Hash_table|

@ In order to be able to track in which file a given user-defined function was
defined, we shall include a record including the location of definition into
the runtime values for such functions.

@< Type declarations needed... @>=
struct source_location
{ unsigned int start_line;
  unsigned short extent, first_col, last_col;
  id_type file;
  source_location(const YYLTYPE& loc); // construct from parser-provided data
  source_location() : start_line(~0u) @+{}
    // sometimes we have no location
  source_location@|
    (const source_location& left, const source_location& right);
    // compute span
  bool undefined() const @+{@; return start_line==~0u; }
};

@ Ultimately we shall want to print |source_location| values in human readable
format.
@< Declarations of functions not for the parser @>=
std::ostream& operator<<(std::ostream& out, const source_location& sl);

@ We try to produce a somewhat compact format in case of items contained in one
source line.

@< Definitions of functions not for the parser @>=
std::ostream& operator<<(std::ostream& out, const source_location& sl)
{
  out << "at " << main_input_buffer->name_of(sl.file) @|
      << ':' << sl.start_line << ':' << sl.first_col << '-';
  if (sl.extent>0)
    out << '-' << sl.start_line+sl.extent << ':';
  return out << sl.last_col;
}

@ The following constructor computes a |source_location| structure, given the
|YYLTYPE| structure that the parser computes. In addition to the information
it provides, we need to get the file name from the |main_input_buffer|.

@< Definitions of functions for the parser @>=
source_location::source_location(const YYLTYPE& loc)
: start_line(loc.first_line) // ``narrow'' to |unsigned int|
, extent(loc.last_line-loc.first_line) // narrow to |unsigned short|
, first_col(loc.first_column)
, last_col(loc.last_column) // these are narrowed too
, file(main_input_buffer->current_file())
{}

@ In some cases we need to compute a text range delimited by two
subexpressions as the source location of an expression. In general these are
virtual expressions not explicitly present in the source, but obtained by
syntactic de-sugaring in the parser, like rewriting an infix operator as an
application to a tuples expression formed of the operands.

@< Definitions of functions for the parser @>=
source_location::source_location
  (const source_location& left, const source_location& right)
@/: start_line(left.start_line)
, extent(right.start_line-left.start_line+right.extent)
@/, first_col(left.first_col), last_col(right.last_col)
, file(left.file)
{ // correct the above computations when one location was undefined
  if (right.undefined())
    *this = left;
  else if (left.undefined())
    *this = right;
}

@ When default-constructing an |expr| (which happens rarely, but is
occasionally used to have a local variable whose actual value will be set
differently according to different branches of execution) we set its |kind| to
|no_expr|, which will ensure proper behaviour (no action) when it gets
assigned to, or destroyed.

We don't have a copy constructor or copy assignment operator, because we
always want to use move semantics for expressions. In full \Cpp11 this
circumstance would be implicitly deduced from the explicit declaration of a
move constructor and move assignment operator, but older compiler versions
would implicitly declare a copy constructor or copy assignment operator if we
didn't declare them deleted. Since those implicitly defined functions do not
do the right thing, this would lead to insidious errors where for
instance a copy constructor gets inserted in a |return| statement (which
moreover one would have expected to be elided) and the corresponding
destruction of the local variable leaves the actual (copied) return value
with dangling pointers. Whence the explicit |delete| declarations below, which
do no harm even if they were implicitly assumed.

@< Methods of |expr| @>=
expr() : kind(no_expr), loc() @+{}
@/expr(const expr& x) = @[delete@];
expr& operator=(const expr& x) = @[delete@];
~expr(); // defined below using a large |switch| statement

@ While we are defining functions to parse expressions, we shall also define a
function |destroy_expr| to clean up the memory occupied by an expression. It
is classified as a parsing function since it is called amongst others by the
parser when popping off tokens at syntax errors.

@< Declarations of functions for the parser @>=
void destroy_expr(expr_p e);

@~The function |destroy_expr| just calls |delete|, which will invoke the
destructor of the |expr| pointed to.

@< Definitions of functions for the parser @>=
void destroy_expr(expr_p p) @+ {@; delete p; }

@ The actual definition of the destructor is distributed among the different
variants of |expr|.

@< Definitions of functions not for the parser @>=
expr::~expr()
{
  if (kind!=no_expr)
  {
    switch (kind)
    {@; @< Cases for destroying an expression @>
      @+ case no_expr: {}
    }
    kind = no_expr;
  }
}

@ We define a move constructor and move assignment, a (rarely sued) |swap|
method. Our first purpose in defining these was to eliminate the need for
copying as much as possible by detecting now forbidden copy constructions and
assignments hidden in the code.

@< Methods of |expr| @>=
void set_from (expr& other);
expr (expr&& other);
void operator= (expr&& other);
void swap (expr& other);

@ We do not assume that |expr| is trivially copyable, even though most
variants (and all in case the macro |incompletecpp1| is defined) are raw
pointers. Therefore copying requires a large case distinction (which will turn
out to be very boring). The move constructor and assignment and |swap| method
are easy to define in terms of |set_from|. However move assignment does
involve a subtle point: it could be that |other| is a subexpression of the old
value of |*this|, and it is necessary to move-construct it to a temporary
variable in order to prevent if from destruction when calling |this->~expr|.

@< Definitions of functions not for... @>=
void expr::set_from (expr& other)
{ assert(kind==no_expr);
  switch (kind=other.kind)
    {@; @< Cases for copying an expression from |other| @>
       @+ case no_expr: {}
    }
  other.kind = no_expr;
  loc = other.loc;
}
@)
expr::expr (expr&& other) : kind(no_expr) @+
{@; set_from(other); }

void expr::operator= (expr&& other)
{@;
  if (this!=&other)
  { expr tmp(std::move(other)); // move |other| out of destruction's way
    this->~expr(); // now destruct old value
    set_from(tmp);
  }
}
@)
void expr::swap (expr& other)
{@; expr tmp(std::move(*this));
  set_from(other);
  other.set_from(tmp);
}

@ In parallel, we also define a function to print the expressions once parsed;
this provides a useful test to see if what we have read in corresponds to what
was typed, and this functionality will also be used in producing error
messages.

@< Declarations of functions not for the parser @>=
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


@*1 Atomic expression. The simplest expressions are the ones consisting of
just one symbol, namely the denotations (constants), applied identifiers, and
the symbol \.\$ referring to the last value computed.

@*2 Denotations.
%
We call simple constant expressions ``denotations''. There are recognised by
the scanner, and the parser will build an appropriate node for them, which
just stores the constant value denoted, or in case of integers the string of
characters. The latter will be converted below to a |arithmetic::big_int| to be
stored in the |expr| value.

@< Includes needed... @>=
#include "bigint.h" // defines |arithmetic::big_int|

@ To keep the size of the |union| in |expr| small, we only use variants that
are no bigger than a pointer, storing a pointer when more than that is needed.
For Boolean denotations, the value itself will fit comfortably inside the
|struct expr@;|. For integers and strings we store a |std::string| value since
it turns out to have the size of a single pointer; it will an example of how
the special member functions of |expr| should handle non-POD variants (they
can neither be assigned to non-initialised memory without calling a
constructor, nor be left in memory that will be reclaimed without calling a
destructor). If \Cpp11 support is incomplete, we store a raw character pointer
|char*@[@]@;| instead. In the integer case we leave the conversion from
character string to numeric value until after type checking.

@< Variants of ... @>=

bool bool_denotation_variant;
#ifdef incompletecpp11
const char* str_denotation_variant;
#else
std::string str_denotation_variant;
#endif

@~Each of the three types of denotation has a tag identifying it.

@< Enumeration tags for |expr_kind| @>=
integer_denotation, string_denotation, boolean_denotation, @[@]

@ For most of the variants there is a corresponding constructor that
placement-constructs the constant value into that variant. For the Boolean
case we might alternatively have assigned to the field, but not for the
integer and string variants. Constructors whose arguments other than |loc| have
easily convertible argument types like |int| or |bool| are given an additional
tag argument to avoid accidentally invoking an unintended constructor. This
also removes in particular the danger of accidentally invoking the Boolean
constructor by accidentally passing some unrelated type by pointer (which
would implicitly convert to |bool|).

@< Methods of |expr| @>=
  struct int_tag @+{}; @+
  struct bool_tag @+{}; @+
  struct string_tag @+{};
  expr (bool b, const YYLTYPE& loc, bool_tag)
@/: kind(boolean_denotation), bool_denotation_variant(b), loc(loc) @+{}
#ifdef incompletecpp11
  expr(std::string&& s, const YYLTYPE& loc, int_tag)
@/: kind(integer_denotation)
  , str_denotation_variant(std::strcpy(new char[s.size()+1],s.c_str()))
  , loc(loc) @+{}
  expr(std::string&& s, const YYLTYPE& loc, string_tag)
@/ : kind(string_denotation)
   , str_denotation_variant(std::strcpy(new char[s.size()+1],s.c_str()))
   , loc(loc) @+{}
#else
  expr(std::string&& s, const YYLTYPE& loc, int_tag)
@/: kind(integer_denotation)
  , str_denotation_variant(std::move(s))
  , loc(loc) @+{}
  expr(std::string&& s, const YYLTYPE& loc, string_tag)
@/ : kind(string_denotation)
   , str_denotation_variant(std::move(s))
   , loc(loc) @+{}
#endif

@~For more explicit construction of these variants in a dynamically allocated
|expr| object, we provide the functions below. Note that |make_int_denotation|
and |make_string_denotation| untypically take a pointer to a |std::string| as
argument; this is because the value of a string token is maintained on the
parser stack as a variant of a |union| and so cannot be a type with a
nontrivial destructor (destruction being handled by explicit calls).

@< Declarations of functions for the parser @>=
expr_p make_int_denotation (std::string* val_p, const YYLTYPE& loc);
expr_p make_bool_denotation(bool val, const YYLTYPE& loc);
expr_p make_string_denotation(std::string* val_p, const YYLTYPE& loc);

@~The definition of these functions is quite easy, as will be typical for
node-building functions. Since the parser stack only contains plain types, the
argument to |make_int_denotation| and |make_string_denotation| are raw
pointers to |std::string|, for which they will take care to call |delete|, so
that the temporarily allocated variable gets cleaned up (as this happens
whenever a string token becomes an expression, the parser only has to perform
clean-up during error recovery, when such a token gets popped without becoming
an expression). We can move from the pointed-to string before it gets
destroyed, avoiding to copy the string.

As a special service to the parser, which for some rules will generate
integer denotations with the value $0$ not coming from the program text (for
instance for certain omitted arguments), we allow calling
|make_int_denotation| with |nullptr| as first argument, that will be replaced
by a string specifying $0$ here. It happens that the empty string converts to
the value $0$, so we provide that string rather than |"0"| here.

@< Definitions of functions for the parser @>=
expr_p make_int_denotation (std::string* val_p, const YYLTYPE& loc)
{ std::unique_ptr<std::string>p(val_p); // this ensures clean-up
  if (val_p==nullptr)
    p.reset(new std::string); // empty string will convert to $0$
  return new expr(std::move(*p),loc,expr::int_tag());
}

expr_p make_bool_denotation(bool val, const YYLTYPE& loc)
{@; return new expr (val,loc,expr::bool_tag()); }

expr_p make_string_denotation(std::string* val_p, const YYLTYPE& loc)
{@; expr_p result=new expr(std::move(*val_p),loc,expr::string_tag());
  delete val_p;
  return result;
}

@ For Boolean denotations there is nothing to destroy. For integer and string
denotations however we must destroy the object held through the variant. Since
there is no level of pointer, we should call the destructor for the variant
itself explicitly (in other variants we shall see raw pointers as variants,
for which calling |delete| is the proper action). It was quite a puzzle to
find the right syntax for that, because of what is essentially a bug
in \.{gcc}, namely with |string| in place of |basic_string<char>| below, the
look-up of the destructor fails.

@s basic_string string

@< Cases for destroying... @>=
case boolean_denotation: break;
case integer_denotation:
case string_denotation:
#ifdef incompletecpp11
  delete[](str_denotation_variant); break;
#else
  str_denotation_variant.~basic_string<char>(); break;
#endif

@ In the |expr::set_from| method we changed variants both in |*this| (which
was |no_expr|) and in |other| (which will become |no_expr|). This means that
for non-POD type we must combine construction into a variant of |*this| and a
destruction of that variant of |other|. We use move construction for
efficiency, and it will probably leave an empty shell to be destructed, but
this does not mean we can omit the destruction.

@< Cases for copying... @>=
  case boolean_denotation:
    bool_denotation_variant = other.bool_denotation_variant; break;
  case integer_denotation:
  case string_denotation:
#ifdef incompletecpp11
    str_denotation_variant = other.str_denotation_variant;
    other.str_denotation_variant=nullptr;
#else
    new (&str_denotation_variant)
    std::string(std::move(other.str_denotation_variant));
    other.str_denotation_variant.~basic_string<char>();
#endif
  break;

@ To print an integer or Boolean denotation we just print its variant field;
for Boolean denotations this requires making sure that the stream has its
|boolalpha| status set, which we do on the fly here. For string denotations we
print the stored string enclosed in quotes.

@< Cases for printing... @>=
case integer_denotation: out << e.str_denotation_variant; break; // no quotes
case boolean_denotation:
   out << std::boolalpha << e.bool_denotation_variant; break;
case string_denotation:
  out << '"' << e.str_denotation_variant << '"'; break;

@*2 Applied identifiers, the last value computed, break, return, die.
%
For representing applied identifiers, we use the integer type |id_type|
defined above. Their tag is |applied_identifier|. An expression that behaves
somewhat similarly is `\.\$', which stands for the last value computed.
Finally we have expressions |break|, |return e| and \&{die} that are
type-less, and whose evaluation breaks from the most current loop respectively
aborts evaluation.

@< Enumeration tags for |expr_kind| @>=
applied_identifier,
last_value_computed,break_expr,return_expr,
die_expr, @[@]

@ For |applied_identifier|s we just store their code, for |break_expr| their
depth, for |return_expr| a pointer to another |expr|; for
|last_value_computed| and |die_expr| nothing at all.
@< Variants of ... @>=
id_type identifier_variant;
unsigned break_variant;
expr_p return_variant;

@ We need new tags here to define new constructors, which for the rest are
straightforward.

@< Methods of |expr| @>=
  struct identifier_tag @+{}; @+
  struct dollar_tag @+{};
  struct break_tag @+{@; unsigned depth; };
  struct return_tag @+{};
  struct die_tag @+{};
@)
  expr(id_type id, const YYLTYPE& loc, identifier_tag)
@/: kind(applied_identifier), identifier_variant(id), loc(loc) @+{}
  expr (const YYLTYPE& loc, dollar_tag)
  : kind(last_value_computed), loc(loc) @+{}
  expr (const YYLTYPE& loc, break_tag t)
  : kind(break_expr), break_variant(t.depth), loc(loc) @+{}
  expr (expr_p exp, const YYLTYPE& loc, return_tag t)
  : kind(return_expr), return_variant(exp), loc(loc) @+{}
  expr (const YYLTYPE& loc, die_tag)
  : kind(die_expr), loc(loc) @+{}

@ As usual there are interface function to the parser.

@< Declarations of functions for the parser @>=
expr_p make_applied_identifier (id_type id, const YYLTYPE& loc);
expr_p make_dollar(const YYLTYPE& loc);
expr_p make_break(unsigned n,const YYLTYPE& loc);
expr_p make_die(const YYLTYPE& loc);
expr_p make_return(expr_p exp,const YYLTYPE& loc);

@~These function have rather simple definitions.

@< Definitions of functions for the parser @>=
expr_p make_applied_identifier (id_type id, const YYLTYPE& loc)
 {@; return new expr(id,loc,expr::identifier_tag()); }

expr_p make_dollar (const YYLTYPE& loc)
@+{@; return new expr(loc,expr::dollar_tag()); }
expr_p make_break (unsigned n,const YYLTYPE& loc)
{@; return new expr(loc,@[expr::break_tag{n}@]); }
expr_p make_return (expr_p exp,const YYLTYPE& loc)
{@; return new expr(exp,loc,@[expr::return_tag{}@]); }
expr_p make_die (const YYLTYPE& loc)
@+{@; return new expr(loc,expr::die_tag()); }

@~Like for integer and boolean denotations, there is nothing to destroy here,
except for a |return_expr|.

@< Cases for destroying... @>=
case applied_identifier:
case last_value_computed:
case break_expr:
case die_expr: break;
case return_expr: delete return_variant; break;

@ Having a POD type variant, copying an applied identifier or break can be
done by assignment; the other two cases have no data at all.

@< Cases for copying... @>=
case applied_identifier: identifier_variant=other.identifier_variant; break;
case break_expr: break_variant=other.break_variant;
case return_expr: return_variant=other.return_variant;
case last_value_computed: case die_expr: break;

@~To print an applied identifier, we look it up in the main hash table. We
also print the other cases as the user wrote them.

@< Cases for printing... @>=
case applied_identifier:
  out << main_hash_table->name_of(e.identifier_variant);
break;
case last_value_computed: out << '$'; @q$@> break;
case break_expr:
{@;
  out << " break ";
  if (e.break_variant>0)
    out << e.break_variant << ' ';
}
  break;
case return_expr: @+{@; out << " return " << *e.return_variant; }
  break;
case die_expr: out << " die "; break;

@*1 Expression lists.
%
We use expression lists in various contexts, so they will occur as variant
of~|expr|.

@< Includes needed... @>=
#include "sl_list.h" // lists are used in parsing types

@~Historically implementing these lists using the Atlas class template
|containers::simple_list| was the first time a non-POD variant of |expr| was
introduced, and it required passage to \Cpp11 as well as substantial
refactoring of the code to make this possible and safe. However, it ultimately
turned out to be simpler to actually store a raw pointer version of the list
(obtained by calling its |release| method) in |expr|; especially so if
alternative compilation with compilers that (still) forbid non-POD variants of
a union is to be supported, but even independently of that there is little
benefit from storing smart pointers. In the mean time this illustrates the
usefulness of having a container type with a lightweight conversion to raw
pointer and back.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef containers::simple_list<expr> expr_list;
typedef containers::sl_node<expr>* raw_expr_list;
  // raw counterpart for use by parser
@)
typedef containers::weak_sl_list_const_iterator<expr> wel_const_iterator;
typedef containers::weak_sl_list_iterator<expr> wel_iterator;
  // wel = weak expression list

@ Any syntactic category whose parsing value is a list of expressions will use
the variant |sublist|.

@< Variants of ... @>=
raw_expr_list sublist;

@~There are two such syntactic categories: tuple displays and list displays,
which are written with parentheses respectively with square brackets.
These are rather different from a type-checking point of view (in a list
display the types of all components must be the same), but for
constructing the parse tree, there is no difference, we just need a different
tag to mark the distinction.

@< Enumeration tags for |expr_kind| @>= tuple_display, list_display, @[@]

@~We provide constructors for both types of lists. Most constructors expect
the parser to directly or indirectly (though an intermediately constructed
|source_location| structure) provide source location information, but the
final constructor (without |loc| argument) extracts the location information
for the component expressions (this is intended for situations where the
expression list does not correspond to an actual enclosed expression, for
instance combined right hand sides of parallel \&{let}-declarations).

@< Methods of |expr| @>=
struct tuple_display_tag @+ {};@+ struct list_display_tag @+{};
expr(expr_list&& nodes, tuple_display_tag, const YYLTYPE& loc)
@/: kind(tuple_display)
  , sublist(nodes.release())
  , loc(loc)
@+{}
expr(expr_list&& nodes, list_display_tag, const YYLTYPE& loc)
@/: kind(list_display)
  , sublist(nodes.release())
  , loc(loc)
@+{}
expr(expr_list&& nodes, tuple_display_tag, const source_location& loc)
@/: kind(tuple_display)
  , sublist(nodes.release())
  , loc(loc)
@+{}
expr(expr_list&& nodes, tuple_display_tag);

@ To determine the extent of a list of expressions, we try to find the first
and the last ones, if any.

@< Definitions of functions for the parser @>=
expr::expr(expr_list&& nodes, tuple_display_tag)
@/: kind(tuple_display)
  , sublist(nodes.release())
  , loc()
{
  if (sublist==nullptr)
    return; // cannot extract something from nothing
  const source_location& left(sublist->contents.loc);
  wel_const_iterator it(sublist);
  while (not std::next(it).at_end())
    ++it;
  loc = source_location(left,it->loc);
}


@ To build an |exprlist_node|, we provide a function |make_exprlist_node| to
combine an expression with a |raw_expr_list|. To start off the construction,
one may use |raw_expr_list(nullptr)| for the empty list. Often it will be
practical to use right recursive grammar rule that build lists backwards
(which means the source expression corresponding to~|e| will actually be to
the right of the one corresponding to~|l|), so we provide a reversal function
to get the proper ordering once the end of the list is reached. Finally we
provide the wrapping function |wrap_expr_list| for list displays.

@< Declarations of functions for the parser @>=
raw_expr_list make_exprlist_node(expr_p e, raw_expr_list l);
raw_expr_list reverse_expr_list(raw_expr_list l);
expr_p wrap_tuple_display(raw_expr_list l, const YYLTYPE& loc);
expr_p wrap_list_display(raw_expr_list l, const YYLTYPE& loc);

@~The function |reverse_expr_list| is easy as well. Of the two possible and
equivalent list reversal paradigms, we use the ``hold the head'' style which
starts with |t=l|; the other option is ``hold the tail'', which would start
with |t=l->next|. Either one would do just as well. We do not call
|reverse_expr_list| from |wrap_expr_list| although the two are usually
combined, since whether or not the list should be reversed can only be
understood when the grammar rules are given.

@< Definitions of functions for the parser @>=
raw_expr_list make_exprlist_node(expr_p e, raw_expr_list raw)
  {@; expr_ptr saf(e); expr_list l(raw);
    l.push_front(std::move(*e));
    return l.release();
  }

raw_expr_list reverse_expr_list(raw_expr_list raw)
{@; expr_list l(raw); l.reverse(); return l.release(); }
@)
expr_p wrap_tuple_display(raw_expr_list l, const YYLTYPE& loc)
{@; return new expr(expr_list(l),expr::tuple_display_tag(),loc); }
expr_p wrap_list_display(raw_expr_list l, const YYLTYPE& loc)
{@; return new expr(expr_list(l),expr::list_display_tag(),loc); }

@ Destroying a tuple display or list display is easily defined.

@< Cases for destroying... @>=
  case tuple_display:
  case list_display: delete sublist; break;

@ Destroying lists of expressions will also be done via a function callable
from the parser, as it may need to discard tokens holding such lists; these
lists are at that point represented by raw pointers.

@< Declarations of functions for the parser @>=
void destroy_exprlist(raw_expr_list l);

@~Its definition is easy; we convert the raw pointer from the parser into an
|expr_list| then simply let it die, which cleans up the list. This is
equivalent to calling |delete l| directly, but depends less on knowing
implementation details (though the fact that |expr_list| is just a smart
pointer version of |raw_expr_list| is not really a secret). We enclose
the case in additional parentheses to avoid interpreting it as a declaration.

@< Definitions of functions for the parser @>=
void destroy_exprlist(raw_expr_list l)
@+{@; (expr_list(l)); }

@ In general copying can be obtained by move construction followed by
destruction of (what remains of) the original value. However for raw pointers
simply assigning |sublist=other.sublist| works. We do not bother to assign
|nullptr| to |other.sublist|, since immediately after executing the code below
we shall set |other.kind=no_expr| which effectively forgets the remaining raw
pointer.

@< Cases for copying... @>=
  case tuple_display:
  case list_display: sublist = other.sublist; other.sublist=nullptr; break;


@ Printing tuple displays and list displays is entirely similar, using
parentheses in the former case and brackets in the latter..

@< Cases for printing... @>=
case tuple_display:
{
  const raw_expr_list l=e.sublist;
  wel_const_iterator it(l);
  if (it.at_end())
    out << "()";
  else
  { out << '(' << *it;
    while (not (++it).at_end())
      out << ',' << *it;
    out << ')';
  }
}
break;
case list_display:
{
  const raw_expr_list l=e.sublist;
  wel_const_iterator it(l);
  if (it.at_end())
    out << "[]";
  else
  { out << '[' << *it;
    while (not (++it).at_end())
      out << ',' << *it;
    out << ']';
  }
}
break;

@ Length 0 tuple displays are sometimes used as a substitute expression where
nothing useful is provided, for instance for a missing else-branch. They are
easy to build, but for recognising them later, it is useful to have a function
at hand.

@< Declarations of functions not for the parser @>=
bool is_empty(const expr& e);

@~The implementation is of course straightforward.

@< Definitions of functions not for the parser @>=
bool is_empty(const expr& e)
  {@; return e.kind==tuple_display and e.sublist==nullptr; }

@*1 Function applications.
%
Another recursive type of expression is the function application. Since it
will store two sub-expressions by value, we must use a pointer when including
it as variant of |expr|. We use a raw pointer, as there is not much point in
having a smart pointer as a variant field in a |union|: the variant does not
benefit from automatic destruction, so the destructor |expr::~expr| would
still have to explicitly call the destructor for the variant; with a raw
pointer it will just directly called |delete| for the pointer.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct application_node* app;

@~Since \.{axis} has tuples, we convene that every function call takes just
one argument. So we define an |application_node| to contain a function and an
argument expression; both are allowed to be arbitrary expressions, though the
function is often an applied identifier, and the argument is often a tuple
display.

The moving constructor does what the braced initialiser-list syntax would do
by default; it is present only for backward compatibility \.{gcc}~4.6.

@< Structure and typedef declarations for types built upon |expr| @>=
struct application_node
{ expr fun; @+ expr arg;
@)
  application_node(expr&& fun, expr&& arg)
  : fun(std::move(fun)), arg(std::move(arg)) @+{}
};

@ The tag used for expressions that will invoke a built-in function is
|function_call|. It is used as internal representation for operator
applications (entered as formulae) as well.

@< Enumeration tags for |expr_kind| @>= function_call, @[@]

@~An evident (and unique) category of |expr| values with an |app| as parsing
value is a function call.

@< Variants of ... @>=
app call_variant;

@ Since function applications are generated both directly by the parser and
indirectly inside parser action functions, where in the latter case an
appropriate |source_location| is computed there, we define two constructors
for building this variant, differing by the type through which they take the
source location.

@< Methods of |expr| @>=
expr(app&& fx, const YYLTYPE& loc)
@/: kind(function_call)
 , call_variant(std::move(fx))
 , loc(loc)
 @+{}
expr(app&& fx, const source_location& loc)
@/: kind(function_call)
 , call_variant(std::move(fx))
 , loc(loc)
 @+{}

@ Building an |application_node| usually combines the function expression with
an |expr_list| for the argument. The argument list will either be packed into a
tuple, or if it has length$~1$ unpacked into a single expression. For tracing
the location of the argument tuple, two additional locations will be provided
by the parser, those of the left and right parentheses that delimit the
argument list (and which are present even for an empty argument list).

Another version is provided if the argument is already grouped as a single
expression, for use with the (generalised) field selector syntax, where the
(projection) function follows the argument. Here no separate location
for the argument is needed, as this information has already been included.

@< Declarations of functions for the parser @>=
expr_p make_application_node(expr_p f, raw_expr_list args,
 const YYLTYPE& loc, const YYLTYPE& left, const YYLTYPE& right);
expr_p make_application_node(expr_p f, expr_p arg, const YYLTYPE& loc);

@~Here for once there is some work to do. Since there are two cases to deal
with, the argument expression is initially default-constructed, after which
either it is set from the unique element of the argument list, or to a tuple
display for the argument list. In the latter case we use a move assignment,
since binding a freshly constructed |expr| to the modifiable lvalue that
|set_from| wants would require introducing a dummy name.

@< Definitions of functions for the parser @>=
expr_p make_application_node(expr_p f, raw_expr_list r_args,
 const YYLTYPE& loc, const YYLTYPE& left, const YYLTYPE& right)
{ expr_ptr ff(f); expr_list args(r_args);
  app a(new application_node @[{ std::move(*ff), expr() }@]);
  if (args.singleton())
    a->arg.set_from(args.front());
  else
    a->arg=expr(std::move(args),expr::tuple_display_tag(),@|
      source_location(source_location(left),source_location(right)));
  return new expr(std::move(a),loc); // move construct application expression
}
@)
expr_p make_application_node(expr_p f, expr_p arg, const YYLTYPE& loc)
{ expr_ptr ff(f); expr_ptr args(arg);
  app a(new application_node { std::move(*ff), std::move(*args) } );
  return new expr(std::move(a),loc); // move construct application expression
}
@ Destroying a raw pointer field just means calling |delete| on it.

@< Cases for destroying... @>=
case function_call: delete call_variant; break;

@ As before, copying a raw pointer is simple.

@< Cases for copying... @>=
   case function_call:
 @/ call_variant=other.call_variant; break;

@ To print a function call, we print the function and a tuple display for the
 arguments, even if there is only one. The function part will also be enclosed
 in parentheses, except in the (very common) case where the function is given
 by an identifier (possibly an operator name). The recursive call of |<<|
 handles all details, such as looking up the function name. We do not yet
 attempt to (re)construct infix formulae, but this is the place where is could
 be done.

@< Cases for printing... @>=
case function_call:
{ const app& a=e.call_variant;
  const expr& fun=a->fun; const expr& arg=a->arg;
  if (fun.kind==applied_identifier) out << fun;
  else out << '(' << fun << ')';
  if (arg.kind==tuple_display) out << arg;
  else out << '(' << arg << ')';
}
break;

@ We shall need to form a function application where the function is accessed
by an applied identifier and the argument is either a single expression or a
tuple of two expressions. However formulae will mostly build the operator
applications they involve by other means; they only call one of the functions
defined here in case of a unary operator that follows a binary operator. The
expression built contains the operator, for which no location information has
yet been processed, so we provide a second location argument |op_loc| in which
the parser can pass this information.

@< Declarations of functions for the parser @>=
expr_p make_unary_call(id_type name, expr_p arg,
   const YYLTYPE& loc, const YYLTYPE& op_loc);
expr_p make_binary_call(id_type name, expr_p x, expr_p y,
   const YYLTYPE& loc, const YYLTYPE& op_loc);

@~In the binary case, we must construct an argument list of two expressions.
We were tempted to initialise |args| below using a two-element initialiser
list, but initialiser lists are incompatible with move semantics. So instead
we push the arguments in reverse order onto an expression list, and then make
a function call using that argument list. The same reasons for manually
constructing arguments lists will apply more often below.

@< Definitions of functions for the parser @>=

expr_p make_binary_call(id_type name, expr_p x, expr_p y,
   const YYLTYPE& loc, const YYLTYPE& op_loc)
{
  expr_ptr xx(x), yy(y); // wrap in smart pointers for exception safety
  const source_location range(xx->loc,yy->loc); // extent of operands
  expr_list args; // start with empty argument list
  args.push_front(std::move(*yy));
  args.push_front(std::move(*xx)); // build up lest back-to-front
  expr arg_pack(std::move(args),expr::tuple_display_tag(),range);
  app a(new application_node @|
    @[{ expr(name,op_loc,expr::identifier_tag()), std::move(arg_pack) }@]);
  return new expr(std::move(a),loc); // move construct application expression
}

@~In the unary case we avoid making an argument list, constructing the
function call directly.

@< Definitions of functions for the parser @>=
expr_p make_unary_call(id_type name, expr_p a,
   const YYLTYPE& loc, const YYLTYPE& op_loc)
{
  expr_ptr aa(a); // wrap in smart pointer for exception safety
  return new expr(new @| application_node(
    expr(name,op_loc,expr::identifier_tag()),std::move(*aa)),loc);
}

@*2 Boolean negation.

It used to be the case that a Boolean negation ``|not E|'' was converted inside
the parser to the equivalent conditional expression ``|if E| \&{then}~|false|
|else|~|true|~\&{fi}'', similarly to the way the Boolean operations |and| and
|or| still are translated (which ensures their ``lazy'' evaluation). For |not|
it is however more efficient to translate into the call of a simple built-in
negation function. We therefore define an expression type for Boolean
negation.

Since there is just an |expr| to store in a negation node, we won't define a
node type, but just store a raw pointer to |expr|. The tag used for these
nodes is |negation_expr|.

@< Enumeration tags for |expr_kind| @>=
negation_expr, @[@]

@ And there is of course a variant of |expr_union| for negations.
@< Variants of ... @>=
expr* negation_variant;

@ There is a constructor for building the new variant.
@< Methods of |expr| @>=
struct negate_tag @+{}; @+
expr(expr* p, const YYLTYPE& loc,negate_tag)
 : kind(negation_expr)
 , negation_variant(p)
 , loc(loc)
@+{}

@ The parser builds negations by calling |make_negation|.

@< Declarations of functions for the parser @>=
expr_p make_negation(expr_p exp, const YYLTYPE& loc);

@~For once, the fact that the parser gives us a raw pointer to the argument
expression is directly useful: we store the pointer right away into a new
|expr| value, with no need to copy and clean up the |expr| pointed to. The day
that we shall be able to modify the parser to pass around the |expr| values it
constructs by value rather than by pointer, this function will become more
complicated where all the others become simpler.

@< Definitions of functions for the parser@>=
expr_p make_negation(expr_p e, const YYLTYPE& loc)
@/{@;
  return new expr(e,loc,expr::negate_tag());
}

@ Here too things are simplified: we can just copy the raw pointer.
@< Cases for copying... @>=
case negation_expr: negation_variant=other.negation_variant;
break;

@ Eventually we want to rid ourselves from the cast.

@< Cases for destroying... @>=
case negation_expr: delete negation_variant; break;

@ Printing cast expressions follows their input syntax.

@< Cases for printing... @>=
case negation_expr: out << " not " << *e.negation_variant;
break;

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
operators, where the former comprises any expression types bound more tightly
than formulae, for instance expressions enclosed in parentheses. For every
operator that comes along we can already determine its left subtree by
comparing priorities with any previous ones. As a consequence of these
comparisons, the right subtree of some previous operators~$\omega$ may also be
completed (if~$\omega$ turned out to bind more strongly than the current
operator); in that case a formula with root~$\omega$ is constructed and
henceforth becomes a single operand. Therefore at any point just after seeing
an operator, there will be a list of pending operators of increasing
priorities (weakly at odd priority levels), each with a complete left operand.
If the formula would terminate after one more operand, then each pending
operator would get as right operand the formula recursively constructed from
everything that follows it. Now when a new operand followed by an
operator~$\omega$ comes along, the pending operators are considered from right
to left; while the operator has a higher priority than~$\omega$, or the same
even priority, it receives the new(est) operand as its right operand, and the
resulting formula becomes the newest operand. Once an operator is encountered
whose priority is too low to capture the operand, the newest operand becomes
the left subtree of~$\omega$, which now becomes the leftmost pending operator.

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
handled immediately), the other is to give them maximal priority only when
they are immediately preceded by another operator. We choose the latter
option, since it allows interpreting \.{-x\pow2} less surprisingly as
$-(x^2)$; somewhat more surprisingly the parentheses in $x+(-1)^n*y$ become
superfluous. The example \.{x\pow-2*y} will then parse as $(x^{-2})*y$, which
also seems reasonable. Note however that with this rule one does have the
surprise that \.{x\pow-y\pow2} parses as $x^{(-y)^2}$.

@ The data type necessary to store these intermediate data during priority
resolutions is a dynamic list of triples subtree-operator-priority. We use a
|mirrored_simple_list|, which is a list re-looked so that its head appears at
|back|, like for |std::stack|. That head represents the rightmost part of the
formula seen do far. To implement the above solution for unary operators, we
allow for the very first pending operator (at the tail of the list) to not
have any left subtree; the expression is left of type |no_expr|, which can be
tested to detect the end of the list.

All nodes will be directly constructed in place, through |emplace_back|, so we
declare a simple moving constructor. For compilers with incomplete \Cpp11
support we explicitly delete the copy constructor assignment to prevent bad
surprises (though the current code would not be affected by their presence).
Even a moving constructor is not needed.

Postfix operators are quite rare in mathematics (the factorial exclamation
mark is the clearest example, though certain exponential notations like the
derivative prime could be considered as postfix operators as well) and more
importantly seem to invariably have infinite priority, so they could be
handled in the parser without dynamic priority comparisons. So here is that
structure:

@< Structure and typedef declarations... @>=
struct formula_node @/{@; expr left_subtree; expr op_exp; int prio;
formula_node(expr&& l,expr&& o, int prio)
 : left_subtree(std::move(l)), op_exp(std::move(o)), prio(prio) @+{}
#ifdef incompletecpp11
@/formula_node(const formula_node& x) = @[delete@];
@/formula_node operator=(const formula_node& x) = @[delete@];
#endif
};

typedef containers::mirrored_simple_list<formula_node> form_stack;
typedef containers::sl_node<formula_node>* raw_form_stack;

@ We define the following functions operating on partial formulae: two to
start them out with a binary or unary operator, the principal one to extend
with a new operand and binary operator, one to finish off the formula with
a final operand, and of course one to clean up.

@< Declarations of functions for the parser @>=
raw_form_stack start_formula
  (expr_p e, id_type op, int prio, const YYLTYPE& op_loc);
raw_form_stack start_unary_formula (id_type op, int prio, const YYLTYPE& op_loc);
raw_form_stack extend_formula
  (raw_form_stack pre, expr_p e,id_type op, int prio, const YYLTYPE& op_loc);
expr_p end_formula (raw_form_stack pre, expr_p e, const YYLTYPE& loc);
void destroy_formula(raw_form_stack s);

@ Starting a formula simply creates an initial node, with a left operand in
the case of a binary formula.
@< Definitions of functions for the parser @>=
raw_form_stack start_formula
   (expr_p e, id_type op, int prio, const YYLTYPE& op_loc)
{ form_stack result;
  result.emplace_back (
   std::move(*expr_ptr(e)), expr(op,op_loc,expr::identifier_tag()), prio );
  return result.release();
}
@)
raw_form_stack start_unary_formula (id_type op, int prio, const YYLTYPE& op_loc)
{ form_stack result;    // leave |left_subtree| empty
  result.emplace_back (expr(), expr(op,op_loc,expr::identifier_tag()), prio );
@/  return result.release();
}

@ Extending a formula involves the priority comparisons and manipulations
indicated above. It turns out |start_formula| could have been replaced by a
call to |extend_formula| with |pre==nullptr|. The second part of the condition
in the while loop is short for |s.front().prio>prio or (s.front().prio==prio
and prio%2==0)|.

@< Definitions of functions for the parser @>=

raw_form_stack extend_formula
  (raw_form_stack pre, expr_p ep,id_type op, int prio, const YYLTYPE& op_loc)
{ expr e(std::move(*expr_ptr (ep))); form_stack s(pre);
  while (not s.empty() and s.front().prio>=prio+prio%2)
    @< Replace |e| by |oper(left_subtree,e)| where |oper| and |left_subtree|
       come from popped |s.front()| @>
  s.emplace_back(std::move(e),expr(op,op_loc,expr::identifier_tag()),prio);
  return s.release();
}

@ Here we make either a unary or a binary operator call. The unary case only
applies for the last node of the stack (since only |start_unary_formula| can
create a node without left operand), but that fact is not used here. Only once
the node is emptied do we pop it with |s.pop_back()|.

@< Replace |e| by |oper(left_subtree,e)|...@>=
{ expr& lt = s.front().left_subtree;
  if (lt.kind==no_expr) // apply initial unary operator
  {
    expr& oper = s.front().op_exp;
    const source_location range(oper.loc,e.loc); // extent of operands
    e = expr(new application_node(std::move(oper),std::move(e)),range);
  }
  else
  {
    expr_list args; // start with empty argument list
    const source_location range(lt.loc,e.loc); // extent of operands
    args.push_front(std::move(e));
    args.push_front(std::move(lt));
    expr arg_pack(std::move(args),expr::tuple_display_tag(),range);
    app a(new application_node @|
      ( std::move(s.front().op_exp), std::move(arg_pack) ));
    e= expr(std::move(a),range); // move construct application expression
  }
  s.pop_back();
}

@ Wrapping up a formula is similar to the initial part of |extend_formula|,
but with an infinitely low value for the ``current priority'' |prio|. So we
can reuse the main part of the loop of |extend_formula|, with just a minor
modification to make sure all nodes get cleaned up after use.

@< Definitions of functions for the parser @>=
expr_p end_formula (raw_form_stack pre, expr_p ep, const YYLTYPE& loc)
{ expr e(std::move(*expr_ptr (ep))); form_stack s(pre);
  while (not s.empty())
    @< Replace |e| by |oper(left_subtree,e)|...@>
  return new expr(std::move(e));
}

@ Destroying a formula stack is straightforward.
@< Definitions of functions for the parser @>=
void destroy_formula(raw_form_stack s)
@+{@; (form_stack(s)); }

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
|sublist.empty()| holds, even though |kind ==0x2|; it turned out that the
possibility to bind no identifiers at all (while providing a void value) has
its uses. However patterns of the form $(x)$, which would give a sublist of
length~$1$, will be forbidden: they would be confusing since $1$-tuples do not
exist.

The structure |raw_id_pat| cannot have a constructor, since it figures in a
bare |union|, where this is not allowed. However the value from them will be
transformed into |id_pat| which is fully equipped with (move) constructors and
destructors. Its default constructor will ensure that constructed nodes will
immediately be safe for possible destruction by |destroy_id_pat| (no undefined
pointers; the |sublist| pointer of |id_pat| is ignored when |kind==0|). Other
constructors are from a |raw_id_pat| reference and from individual components.

@h "axis-types.h"
   // so complete type definitions will be known in \.{parsetree.cpp}

@< Type declarations needed in definition of |struct expr@;| @>=
typedef containers::simple_list<struct id_pat> patlist;
typedef containers::sl_node<struct id_pat>* raw_patlist;
@)
struct raw_id_pat
{ id_type name;
  unsigned char kind;
  // bits 0,1,2: whether pattern has a name, has sublist, is const
  raw_patlist sublist;
};
struct id_pat
{ id_type name;
  unsigned char kind;
  // bits 0,1,2: whether pattern has a name, has sublist, is const
  patlist sublist;
@)
  id_pat() :name(), kind(0x0), sublist() @+{}
  explicit id_pat(const raw_id_pat& x)
  @/ : name(x.name), kind(x.kind)
  , sublist((kind & 0x2)==0 ? nullptr : x.sublist)
  @+{}
  id_pat (id_type n) : name(n), kind(0x1), sublist() @+{}
    // name only,  no constness
  id_pat(patlist&& l): name(), kind(0x2),  sublist(std::move(l)) @+{}
    // no name, no constness
  id_pat (id_type n, unsigned char k, patlist&& l)
  : name(n), kind(k), sublist(std::move(l)) @+{}
@)
  id_pat (const id_pat& x) = @[ delete @];
@/id_pat& operator=(const id_pat& x) = @[ delete @];
#ifdef incompletecpp11
  id_pat (id_pat&& x)
  : name(x.name), kind(x.kind), sublist(std::move(x.sublist)) @+{}
@/id_pat& operator=(id_pat&& x)
  @/{@; name = x.name; kind = x.kind; sublist = std::move(x.sublist);
    return *this;
  }
#else
@/id_pat (id_pat&& x) = @[ default @];
@/id_pat& operator=(id_pat&& x) = @[ default @];
#endif
  raw_id_pat release()
  @+{@; return { name, kind, sublist.release() }; }
};

@ The function |make_pattern_node| to build a node takes a modifiable
reference to a structure |pattern| with the contents of the current node; in
practice this is a local variable in the parser. The argument names and order
reflects the fact that patterns are often recognised by left-recursive rules,
so that a new pattern is tacked onto an existing pattern list (the resulting
list will need reversal when further integrated).

@< Declarations of functions for the parser @>=
raw_patlist make_pattern_node(raw_patlist prev,raw_id_pat& pattern);

@ With the mentioned proviso about order, the function implementation just
assembles the pieces.

@< Definitions of functions for the parser @>=
raw_patlist make_pattern_node(raw_patlist prev,raw_id_pat& pattern)
{@; patlist l(prev); l.push_front(id_pat(pattern)); return l.release(); }

@ Patterns also need cleaning up, which is what |destroy_pattern| and
|destroy_id_pat| will handle, and reversal as handled by |reverse_patlist|.

@< Declarations of functions for the parser @>=
void destroy_pattern(raw_patlist p);
void destroy_id_pat(const raw_id_pat& p);
raw_patlist reverse_patlist(raw_patlist p);

@ The function |destroy_pattern| just converts to a smart pointer, and
|destroy_id_pat| calls it whenever there is a sub-list in a pattern.

@< Definitions of functions for the parser @>=
void destroy_pattern(raw_patlist p)
@+{@; (patlist(p)); }
@)
void destroy_id_pat(const raw_id_pat& p)
@+{@; if ((p.kind & 0x2)!=0)
    destroy_pattern(p.sublist);
}
@)
raw_patlist reverse_patlist(raw_patlist raw)
@+{@; patlist p(raw);
  p.reverse();
  return p.release();
}

@ We can provide a printing function, that can be used from within the one for
|expr|.

@< Declarations of functions not for the parser @>=
std::ostream& operator<< (std::ostream& out, const id_pat& p);

@~Only parts whose presence is indicated in |kind| are printed. We take care
that even if |sublist| is marked as present, it might be null. The
list cannot be of length~$1$ but this fact is not used below.

@< Definitions of functions not for the parser @>=
std::ostream& operator<< (std::ostream& out, const id_pat& p)
{ if ((p.kind & 0x2)!=0)
  { out << '(';
    if (not p.sublist.empty())
    { auto it = p.sublist.begin();
      out << *it;
      while (not (++it).at_end())
        out << ',' << *it;
    }
    out << ')';
  }
  if ((p.kind&0x3)==0x3) // both parts present
    out << ':';
  if ((p.kind&0x4)!=0) // identifier declared constant
    out << '!';
  if ((p.kind & 0x1)!=0)
    out << main_hash_table->name_of(p.name);
  return out;
}

@*1 Let expressions.
%
We now consider let-expressions, which introduce and bind local identifiers,
which historically were a first step towards having user defined functions.
Indeed a let-expression used to be implemented as a user-defined function that
were immediately called with the (tuple of) values bound in the
let-expression (the implementation of let-expressions has been somewhat
optimised since). The reason to have had let expressions before user-defined
functions was that it could be done before user-specified types were
introduced in the language; in a let expression types of variables are deduced
from the values provided, whereas function parameters need to have explicitly
specified types.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct let_expr_node* let;

@~After parsing, let-expression will have a single let-binding followed by a
body giving the value to be returned. During parsing however, we may form
intermediate values containing lists of let-bindings, that will later be
converted into a single one with a tuple as left hand side. We therefore
define a list type |let_list| for a list of bindings, and a structure
|let_expr_node| for a complete let-expression, containing (the components of)
only one binding, and containing in addition a body.

The moving constructor does what the braced initialiser-list syntax would do
by default; it is present only for backward compatibility \.{gcc}~4.6.

@< Structure and typedef declarations for types built upon |expr| @>=
struct let_pair { id_pat pattern; expr val;
#ifdef incompletecpp11
  let_pair (const let_pair&) = @[delete@];
  let_pair& operator=(const let_pair&) = @[delete@];
  let_pair (let_pair&& x)
  : pattern(std::move(x.pattern)), val(std::move(x.val)) @+{}
  let_pair (id_pat&& p, expr&& v)
  : pattern(std::move(p)), val(std::move(v)) @+{}
#endif
};
typedef containers::simple_list<let_pair> let_list;
typedef containers::sl_node<let_pair>* raw_let_list;
@)
struct let_expr_node
{ id_pat pattern; expr val; expr body;
@)
  let_expr_node(id_pat&& pattern, expr&& val, expr&& body)
@/: pattern(std::move(pattern)), val(std::move(val)), body(std::move(body))@+{}
};

@ The tag used for let-expressions is |let_expr|.

@< Enumeration tags for |expr_kind| @>= let_expr, @[@]

@~Without surprise there is a class of |expr| values with a |let| as parsing
value.

@< Variants of ... @>=
let let_variant;

@ There is a constructor for building this variant.
@< Methods of |expr| @>=
expr(let&& declaration, const YYLTYPE& loc)
 : kind(let_expr)
 , let_variant(std::move(declaration))
 , loc(loc)
 @+{}

@ For building let-expressions, four functions will be defined. The function
|make_let_node| makes a list of one declaration, while |append_let_node|
appends such a list |cur| (assured to be of length~$1$) to a previously
constructed list |prev| of declarations; finally |make_let_expr_node| wraps up
an entire let-expression. The latter will uses an auxiliary |zip_decls| to do
the main work, which function can then also be used when processing
(global) \&{set} definitions, which do not have any |body|.

@< Declarations of functions for the parser @>=
raw_let_list make_let_node(raw_id_pat& pattern, expr_p val);
raw_let_list append_let_node(raw_let_list prev, raw_let_list cur);
expr_p make_let_expr_node(raw_let_list decls, expr_p body, const YYLTYPE& loc);
std::pair<id_pat,expr> zip_decls(raw_let_list decls);

@ The functions |make_let_node| and |append_let_node| build a list in reverse
order, which makes the latter function a particularly simple one. The purpose
of |make_let_expr_node| is to combine multiple declarations (that were
separated by commas) to one, taking care to reverse the order at the same time
so that the tuple of patterns being declared comes out in the same order as it
was specified in the program. Fortunately it is actually easier to build a
merged list in reverse order.

The argument |cur| for |append_let_node| is always the result of an
application of |make_let_node|, so it is a list of length$~1$. Using this is
just a trick to avoid having to deal with yet another type (corresponding to a
|let_pair|, but demoted to POD) in the parser.

@< Definitions of functions for the parser @>=
raw_let_list make_let_node(raw_id_pat& pattern, expr_p val)
{ let_list result;
  result.push_front ( let_pair { id_pat(pattern), std::move(*expr_ptr(val)) } );
  return result.release();
}

raw_let_list append_let_node(raw_let_list prev, raw_let_list cur)
{@; let_list result(prev);
  result.push_front(std::move(let_list(cur).front()));
  return result.release();
 }
@)
std::pair<id_pat,expr> zip_decls(raw_let_list d)
{ let_list decls(d);
  std::pair<id_pat,expr> result;
  id_pat& pattern=result.first; expr& val=result.second;
  if (decls.singleton()) // single declaration
  @/{@; pattern=std::move(decls.front().pattern);
    val=std::move(decls.front().val);
  }
  else
  { patlist patl; expr_list expl;
    for (auto it=decls.begin(); not decls.at_end(it); ++it)
      // zip open |decls|, reversing
    {@; patl.push_front(std::move(it->pattern));
      expl.push_front(std::move(it->val));
    }
    pattern = id_pat(std::move(patl));
    val=expr(std::move(expl),expr::tuple_display_tag());
      // make a tuple expression
  }
  return result;
}
expr_p make_let_expr_node(raw_let_list d, expr_p b, const YYLTYPE& loc)
{
  expr_ptr bb(b);
  expr& body=*bb; // ensure exception safety
  auto pat_expr = zip_decls(d);
@/return new expr (new  @| let_expr_node
  { std::move(pat_expr.first), std::move(pat_expr.second),std::move(body) }
  ,loc);
}

@ For the raw pointer |let|, copying is done by assignment, as before.

@< Cases for copying... @>=
 case let_expr:
   let_variant=other.let_variant;
 break;

@~While |let_list| is only handled at the outermost level during parsing, it
is |let| smart pointers that get built into |expr| values. As for all
variants of |expr|, calling the destructor must be done explicitly.

@< Cases for destroying... @>=
case let_expr: delete let_variant; break;

@ Destroying lists of declarations will be done in a function callable from the
parser, like |destroy_exprlist|.

@< Declarations of functions for the parser @>=
void destroy_letlist(raw_let_list l);

@~Like in |destroy_exprlist|, merely reconstructing the non-raw list and
letting that be destructed suffices to recursively destroy all nodes of the
list, and anything accessed from them.

@< Definitions of functions for the parser @>=
void destroy_letlist(raw_let_list l)
@+{@; (let_list(l)); }

@ To print a let-expression we just do the obvious things, for now without
worrying about parentheses.

@< Cases for printing... @>=
case let_expr:
{ const let& lexp=e.let_variant;
  out << "let " << lexp->pattern << '=' << lexp->val <<	 " in " << lexp->body;
}
break;

@*1 Types and user-defined functions.
%
As was indicated above, let-expressions were introduced before user-defined
functions because it avoids to problem of having to specify types in the user
program. When defining function one does have to specify parameter types,
since in this case there is no way to know those types with certainty: in an
interactive program we cannot wait for \emph{usage} of the function to deduce
the types of its parameters, and such things as type coercion and function
overloading would make such type deduction doubtful even if it could be done).
Types are also an essential ingredient of casts. Types have an elaborate
(though straightforward) internal structure, and the necessary constructing
functions like |make_tuple_type| are defined in the module \.{axis-types.w}
rather than here.

When the parser was compiled as \Cee~code, we were forced to use void pointers
in the parser, to masquerade for the actual pointers to \Cpp~types; numerous
static casts were then used to get the proper pointer types from them. Now
that this is no longer necessary, everything has been reformulated in terms of
the actual pointer types. Something that remains (for now) is the avoidance of
smart pointers for types while being handled in the parser, since other
pointers for expressions that it handles are not smart pointers either.

@< Includes needed... @>=
#include "axis-types.h" // parsing types need |type_p| and such

@ Most functionality for these types is given in \.{axis-types.w}; we just
need to say how to destroy them.

@< Declarations of functions for the parser @>=
void destroy_type(type_p t);
void destroy_type_list(raw_type_list t);

@ Since |type_ptr| and |type_list| are fully equipped types, it suffices to
convert to them and forget.

@< Definitions of functions for the parser @>=

void destroy_type(type_p t)@+ {@; (type_ptr(t)); }
void destroy_type_list(raw_type_list t)@+ {@; (type_list(t)); }
  // recursive destruction

@ For user-defined functions we shall use a structure |lambda_node|.
@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct lambda_node* lambda;

@~It contains a pattern for the formal parameter(s), its type (a smart pointer
defined in \.{axis-types.w}), an expression (the body of the function), and
finally the source location of the function.


@< Structure and typedef... @>=
struct lambda_node
{ id_pat pattern; type_expr parameter_type; expr body;
@)
  lambda_node(id_pat&& pattern, type_expr&& type, expr&& body)
  @/: pattern(std::move(pattern))
  , parameter_type(std::move(type))
  , body(std::move(body))
@+{}
};

@ The tag used for user-defined functions is |lambda_expr|.
@< Enumeration tags... @>=
lambda_expr,@[@]

@ We introduce the variant of |expr| as usual.
@< Variants of ... @>=
lambda lambda_variant;

@ There is a constructor for building lambda expressions.
@< Methods of |expr| @>=
expr(lambda&& fun, const YYLTYPE& loc)
 : kind(lambda_expr)
 , lambda_variant(std::move(fun))
 , loc(loc)
@+{}

@ There is as usual a function for constructing a node, to be called
by the parser.

@< Declarations of functions for the parser @>=
expr_p make_lambda_node(raw_patlist pat_l, raw_type_list type_l, expr_p body,
 const YYLTYPE& loc);

@ There is a twist in building a lambda node, similar to what we saw for
building let-expressions, in that for syntactic reasons the parser passes
lists of patterns and their types, rather than passing single ones. We must
distinguish the case of a singleton, in which case the head node must be
unpacked, and the multiple case, where a tuple pattern and type must be
wrapped up from the lists. In the former case, |fun->parameter_type| wants to
have a pointer to an isolated |type_expr|, but the head of |type_l| is a
|type_node| that contains a |type_expr| as its |t| field; making a (shallow)
copy of that field is the easiest way to obtain an isolated |type_expr|. After
the copy, destruction of |type_l| deletes the original |type_node|. In the
latter case we apply list reversal here to both pattern list and type list.

@< Definitions of functions for the parser @>=
expr_p make_lambda_node(raw_patlist p, raw_type_list tl, expr_p b,
 const YYLTYPE& loc)
{
  patlist pat_l(p);
  type_list type_l(tl);
  expr_ptr body_p(b);
  expr& body=*body_p; // safety
  id_pat pattern; type_expr parameter_type;
  if (type_l.singleton())
@/{@; pattern=std::move(pat_l.front());
    parameter_type = std::move(type_l.front());
  }
  else
@/{ pat_l.reverse(); pattern=id_pat(std::move(pat_l));
  @/type_l.reverse(); parameter_type=type_expr(std::move(type_l));
  // make tuple type
  }
  return new expr(lambda(new@| lambda_node
      (std::move(pattern),std::move(parameter_type),std::move(body))),loc);
}

@ Since |lambda| is a raw pointer, we can just assign.
@< Cases for copying... @>=
case lambda_expr: lambda_variant=other.lambda_variant;
break;

@ And we must of course take care of destroying lambda expressions, which is
done correctly by the implicit destructions provoked by calling |delete|.

@< Cases for destroying... @>=
case lambda_expr: delete lambda_variant; break;

@ Because of the above transformations, lambda expressions are printed with
all parameter types grouped into one tuple (unless there was exactly one
parameter). In case of no parameters at all, we do not print |"(() ())"| but
rather |"@@"| (which though unrelated of other uses of .\\'@@ is easy to get
used to).

@< Cases for printing... @>=
case lambda_expr:
{ const lambda& fun=e.lambda_variant;
  if (fun->parameter_type==void_type)
    out << '@@';
  else
    out << '(' << fun->parameter_type << ' ' << fun->pattern << ')';
  out << ':' << fun->body;
}
break;


@*1 Control structures.

@*2 Conditional expressions.
Of course we need if-then-else expressions.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct conditional_node* cond;

@~The parser handles \&{elif} constructions, so we only need to handle the
basic two-branch case. Nonetheless we prefer to have the two branches as a
linked list, so as to dispose of a |raw_expr_list| value that can be passed to
the same balancing function as for instance the components of a list display
(it would be quite tedious to reorganise the contents of two explicit |expr|
fields to make them occur in a list without temporarily removing the
originals, given that copying |expr| values is expensive). On the other hand,
given that the list always has length~$2$, we can represent the list by its
initial node.

@< Structure and typedef declarations for types built upon |expr| @>=
struct conditional_node
{ expr condition; containers::sl_node<expr> branches;
@)
  conditional_node(expr&& condition, expr&& then_branch, expr&& else_branch)
@/: condition(std::move(condition))
  , branches(std::move(then_branch))
  {@; branches.next.reset(new
      containers::sl_node<expr>(std::move(else_branch)));
  }
  @< Other constructor for |conditional_node| @>@;
};

@ The tag used for these expressions is |conditional_expr|.

@< Enumeration tags for |expr_kind| @>= conditional_expr, @[@]

@~The variant of |expr| values with |cond| as parsing value is tagged
|if_variant|.

@< Variants of ... @>=
cond if_variant;

@ There is a constructor for building conditional expressions. As explained
below, the variant |if_variant| will be reused (for case expressions), so
exceptionally we pass the |expr_kind| tag explicitly.
@< Methods of |expr| @>=
expr(expr_kind which, cond&& conditional, const YYLTYPE& loc)
 : kind(which)
 , if_variant(std::move(conditional))
 , loc(loc)
@+{}

@ To build an |conditional_node|, we define a function as usual.
@< Declarations of functions for the parser @>=
expr_p make_conditional_node(expr_p c, expr_p t, expr_p e, const YYLTYPE& loc);

@~Here we pass the |conditional_expr| tag explicitly to the |expr| constructor
above.

@< Definitions of functions for the parser @>=
expr_p make_conditional_node(expr_p c, expr_p t, expr_p e, const YYLTYPE& loc)
{
  expr_ptr cc(c), tt(t), ee(e);
  expr& cnd=*cc; expr& thn=*tt; expr& els=*ee;
@/return new @| expr(conditional_expr, new @|
      conditional_node { std::move(cnd), std::move(thn), std::move(els) },loc);
}

@ We follow the usual coding pattern for copying raw pointers.  We include
some other case labels for expression kinds to be introduced later that have
isomorphic nodes, and therefore use |if_variant|.

@< Cases for copying... @>=
case conditional_expr:
case int_case_expr:
case union_case_expr:
  if_variant=other.if_variant;
break;

@~Again we just need to activate the variant destructor, the rest is
automatic.

@< Cases for destroying... @>=
case conditional_expr:
case int_case_expr:
case union_case_expr:
   delete if_variant; break;

@ To print a conditional expression at parser level, we shall not
reconstruct \&{elif} constructions.

@< Cases for printing... @>=
case conditional_expr:
{ const auto& c=*e.if_variant; out << " if " << c.condition @|
  << " then " << c.branches.contents @|
  << " else " << c.branches.next->contents << " fi ";
}
break;

@*2 Integer case statements.
%
For a long time the language had no multi-way branching expression at all.
Currently we implement a \&{case} expression selecting based on an integer
value, explicit branches being given for a finite segment of the possible
values. With the later introduction of distinguished union type there will be
a richer set of \&{case} expressions.

As a by-effect of representing the branches of a conditional expression as a
non-empty list, the integer case expression becomes structurally
\emph{equivalent} to the conditional expression (this time the
branches can form a non-empty list of any length). So for the parse we just
need a tag value to discriminate the two, but not a new variant in the union.

@< Enumeration tags for |expr_kind| @>= int_case_expr, @[@]

@~On the other hand, we do need a different constructor for
|conditional_node| to handle the more general list. It is in fact more
straightforward than the previous one.

@< Other constructor for |conditional_node| @>=
conditional_node(expr&& condition, containers::sl_node<expr>&& branches)
@/: condition(std::move(condition))
  , branches(std::move(branches))
  @+{}

@
We define a function for constructing the expression as usual.

@< Declarations of functions for the parser @>=
expr_p make_int_case_node(expr_p selector, raw_expr_list ins, const YYLTYPE& loc);

@~This implementation used to be not straightforward, because these case
expressions were converted using a subscription from a list of functions
without arguments. With a dedicated expression type, the handling has now
become a lot simpler. We separate out an auxiliary function |make_case_node|
that does the main work, and that will be reused later for construction of an
isomorphic kind of expression.

There is still an unusual point here in that we move from the complete first
node of |in_list| into the new |conditional_node|, something that is only
possible because we still hold the raw pointer |i| to that node, passed to us
as argument (the type |expr_list| does not allow this).

@< Definitions of functions for the parser @>=
expr_p make_case_node
  (expr_kind kind,expr_p s, raw_expr_list i, const YYLTYPE& loc)
{
  expr_ptr sel(s); expr_list in_list(i);
  assert(not in_list.empty());
  expr& selector=*sel;
  containers::sl_node<expr>& ins = *i; // list of \&{in}-expressions
  return new expr(kind, new @|
                  conditional_node(std::move(selector),std::move(ins)),loc);
}

expr_p make_int_case_node(expr_p s, raw_expr_list i, const YYLTYPE& loc)
{@; return make_case_node(int_case_expr,s,i,loc); }

@ To print an integer case expression at parser level, we do the obvious.

@< Cases for printing... @>=
case int_case_expr:
{ const auto& c=*e.if_variant;
  wel_const_iterator it(&c.branches);
  out << " case " << c.condition << " in " << *it;
  for (++it; not it.at_end(); ++it)
    out  << ", " << *it;
  out << " esac ";
}
break;


@*2 Discrimination case statements.
%
Union values can be tested for their actual variant using a multi-way switch
similar to the case statement. However here branches are not simply
expressions to be evaluated, they behave like functions that will be applied
to the body of the union value being discriminated upon, if it has the
appropriate variant, with argument type matching the type of that variant.

In fact we define two forms of discrimination case statements. A first, more
primitive form, provides exactly what we just described, a separate function
for each branch; it has the same syntactic form as an integer case statement,
except that we write vertical bars instead of commas to make the distinction.
A second form provides more flexibility to the user, but can only be used if
the union type has occurred in a type definition where the different variants
have been given names (tags); these can then be used to identify branches, and
frees the user from the obligation of writing the branches in the same order
as the corresponding variants are listed in the union, and also of explicitly
specifying types (since they can be deduced from the context) for the
identifiers (function arguments) that are introduced at the beginning of each
branch.

The tag used for the first form is |union_case_expr|, that for the second form
|discrimination_expr|.

@< Enumeration tags for |expr_kind| @>=
union_case_expr, discrimination_expr, @[@]

@~The first case can reuse the |conditional_node| structure, just like the
integral case expressions, so all we need is a trivial variation of
|make_int_case_node|.

@< Definitions of functions for the parser @>=

expr_p make_union_case_node(expr_p s, raw_expr_list i, const YYLTYPE& loc)
{@; return make_case_node(union_case_expr,s,i,loc); }

@ Printing an union case expression at parser level is a very slight variation
of the same for an integer case expression.

@< Cases for printing... @>=
case union_case_expr:
{ const auto& c=*e.if_variant;
  wel_const_iterator it(&c.branches);
  out << " case " << c.condition << " in " << *it;
  for (++it; not it.at_end(); ++it)
    out  << " | " << *it;
  out << " esac ";
}
break;

@ But for the second kind we need a more elaborate
parsing structure.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct discrimination_node* disc;

@ Concretely branches have an identifier that indicates the
case, and then a pattern and a body as in a \&{let} statement. In a
|discrimination_node|, the expression being discriminated is called its
|subject|, which is followed by a non-empty list of |branches|.

@< Structure and typedef declarations for types built upon |expr| @>=
struct case_variant
{ id_type label;
  id_pat pattern;
  expr branch;
#ifdef incompletecpp11
  case_variant (const case_variant&) = @[delete@];
  case_variant& operator=(const case_variant&) = @[delete@];
  case_variant (case_variant&& x)
  : label(x.label)
  , pattern(std::move(x.pattern))
  , branch(std::move(x.branch)) @+{}
  case_variant (id_type lab, id_pat&& p, expr&& v)
  : label(lab), pattern(std::move(p)), branch(std::move(v)) @+{}
#endif
};
typedef containers::simple_list<case_variant> case_list;
typedef containers::sl_node<case_variant>* raw_case_list;
@)
struct discrimination_node
{ expr subject; containers::sl_node<case_variant> branches;
@)
  discrimination_node(expr&& subject, raw_case_list branches);
};

@ The |discrimination_node| constructor is less straightforward than usual,
since we want to move from the head node of the |branches| argument. Doing so
requires passing a raw pointer to the constructor (whereas we usually pass
smart pointers at such occasions), which can be used for this purpose. But we
also want to clean up the old node using a smart pointer; as the lifetime of a
temporary inside the initialiser expression equals (the containing
full-expression which is) that initialiser expression, the code
below achieves the desired effect using a comma operator (destroying the smart
pointer of its left hand side takes place only after the construction of
|subject| has moved the data out of the node pointed to).

@< Definitions of functions for the parser @>=
discrimination_node::discrimination_node(expr&& subject, raw_case_list branches)
@/: subject(std::move(subject))
  , branches((case_list(branches),std::move(*branches)))
  @+{}

@ The variant of |expr| values with |disc| as parsing value is tagged
|disc_variant|.

@< Variants of ... @>=
disc disc_variant;

@ There is of course an |expr| constructor for building discrimination
expressions.

@< Methods of |expr| @>=
expr(disc&& discrimination, const YYLTYPE& loc)
 : kind(discrimination_expr)
 , disc_variant(std::move(discrimination))
 , loc(loc)
@+{}

@ To build a |discrimination_node|, we need several functions. The first two
are similar to what we did for \&{let} lists: create or extend a
|raw_case_list| with the data for a single node.

@< Declarations of functions for the parser @>=
raw_case_list make_case_node(id_type sel, raw_id_pat& pattern, expr_p val);
raw_case_list append_case_node
  (raw_case_list prev, id_type sel, raw_id_pat& pattern, expr_p val);
expr_p make_union_case_node
  (expr_p selector, raw_expr_list ins, const YYLTYPE& loc);
expr_p make_discrimination_node
  (expr_p s, containers::sl_node<case_variant>* branches, const YYLTYPE& loc);

@~The implementation is basic list handling. The first two functions build
individual nodes, which the prepend to the list being constructed using
|push_front|. In |make_discrimination_node| we then reverse the list and then
release the pointer, since that is what the constructor for
|discrimination_node| needs (and it takes possession of the node), for reasons
that were explained just above.

@< Definitions of functions for the parser @>=
raw_case_list append_case_node
  (raw_case_list prev, id_type sel, raw_id_pat& pattern, expr_p val)
{ case_list result(prev);
  result.push_front (
    case_variant { sel, id_pat(pattern), std::move(*expr_ptr(val)) } );
@/return result.release();
}

raw_case_list make_case_node(id_type sel, raw_id_pat& pattern, expr_p val)
{@; return append_case_node(nullptr,sel,pattern,val); }
// just append to an empty list

@)
expr_p make_discrimination_node
  (expr_p s, containers::sl_node<case_variant>* b, const YYLTYPE& loc)
{
  expr_ptr subj(s); case_list branches(b);
  branches.reverse(); // undo effect of left-recursive grammar rules
@/return new @| expr(new @|
      discrimination_node { std::move(*subj), branches.release() },loc);
}

@ We follow the usual coding pattern for copying raw pointers.
@< Cases for copying... @>=
case discrimination_expr:
  disc_variant=other.disc_variant;
break;

@~Again we just need to activate the variant destructor, the rest is
automatic.

@< Cases for destroying... @>=
case discrimination_expr: delete disc_variant; break;

@ Destroying lists of variants  will be done in a function callable from the
parser, like |destroy_letlist|.

@< Declarations of functions for the parser @>=
void destroy_case_list(raw_case_list l);

@~Like in |destroy_let_list|, merely reconstructing the non-raw list and
letting that be destructed suffices to recursively destroy all nodes of the
list, and anything accessed from them.

@< Definitions of functions for the parser @>=
void destroy_case_list(raw_case_list l)
@+{@; (case_list(l)); }


@ To print a discrimination expression at parser level, we follow the input
syntax.

@< Cases for printing... @>=
case discrimination_expr:
{ auto& d=*e.disc_variant; raw_case_list l=&d.branches;
  out << " case " << d.subject;

  for (containers::weak_sl_list_const_iterator<case_variant> it(l);
    @| not it.at_end(); ++it)
    if (it->label==type_table::no_id)
      out << " | else " << it->branch;
    else
      out << " |" << main_hash_table->name_of(it->label)
       @| << '(' << it->pattern << "):" << it->branch;
  out << " esac ";
}
break;


@*2 Loops.
%
Loops are a cornerstone of any form of non-recursive programming. The three
flavours are |while| loops, and |for| loops iterating either over a row value,
or over an integer range (counted |for|-loops).

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct while_node* w_loop;
typedef struct for_node* f_loop;
typedef struct cfor_node* c_loop;

@ We shall use a small |BitSet| value to record reversal attributes of for
loops and of slice expressions.

@< Includes needed... @>=
#include "bitset.h"

@~A |while| loop has two parts: a condition (which determines whether an
iteration will be undertaken), and a body (which contributes an entry to the
result). A |for| loop has three parts: a pattern introducing variables, an
expression iterated over (the in-part) and the loop body. A counted |for| loop
(the simple version of |for| loop) has four parts: an identifier, a count, a
lower bound and a body. Apart from this there are some flags to record static
Boolean attributes of loops: for a |while| loop $2$ bits telling whether the
condition is negated and whether the resulting row should be reversed; for a
|for| loop $2$ buts telling whether object should be traversed in reverse and
whether the resulting row should be reversed, and for counted for loops $3$
bits, the first tow having the same meaning as for |for| loops (the ``object
traversed'' being a sequence of consecutive integers) and a third bit telling
whether the loop index is absent.

@< Structure and typedef declarations for types built upon |expr| @>=
struct while_node
{ expr body; // contains a |do|-expression producing condition and body
  BitSet<2> flags;
@)
  while_node(expr&& body, unsigned flags)
@/: body(std::move(body))
  , flags(flags)@+{}
  // backward compatibility for gcc 4.6
};
struct for_node
{ id_pat id; @+ expr in_part; @+ expr body;
  BitSet<2> flags;
@)
  for_node(id_pat&& id, expr&& in_part, expr&& body, unsigned flags)
@/: id(std::move(id))
  , in_part(std::move(in_part))
  , body(std::move(body))
  , flags(flags)@+{}
};
struct cfor_node
{ id_type id; @+ expr count; @+ expr bound; @+ expr body;
  BitSet<3> flags;
@)
  cfor_node
     (id_type id, expr&& count, expr&& bound, expr&& body, unsigned flags)
@/: id(id)
  , count(std::move(count))
  , bound(std::move(bound))
  , body(std::move(body))
  , flags(flags)@+{}
};

@ The tags used for these expressions are |while_expr|, |for_expr| and
|cfor_expr|.

@< Enumeration tags for |expr_kind| @>= while_expr, for_expr, cfor_expr, @[@]

@~The variant of |expr| values with an |w_loop| as parsing value is tagged
|while_variant|.

@< Variants of ... @>=
w_loop while_variant;
f_loop for_variant;
c_loop cfor_variant;

@ There is a constructor for building each type of loop expression.
@< Methods of |expr| @>=
expr(w_loop&& loop, const YYLTYPE& loc)
 : kind(while_expr)
 , while_variant(std::move(loop))
 , loc(loc)
@+{}
expr(f_loop&& loop, const YYLTYPE& loc)
 : kind(for_expr)
 , for_variant(std::move(loop))
 , loc(loc)
@+{}
expr(c_loop&& loop, const YYLTYPE& loc)
 : kind(cfor_expr)
 , cfor_variant(std::move(loop))
 , loc(loc)
@+{}

@ To build a |while_node|, |for_node| or |cfor_node|, here are yet three
more \\{make}-functions.

@< Declarations of functions for the parser @>=
expr_p make_while_node(expr_p b, unsigned flags, const YYLTYPE& loc);
expr_p make_for_node
  (raw_id_pat& id, expr_p ip, expr_p b, unsigned flags, const YYLTYPE& loc);
expr_p make_cfor_node
 (id_type id, expr_p count, expr_p bound, expr_p b, unsigned flags,
  const YYLTYPE& loc);

@ They are quite straightforward, as usual.

@< Definitions of functions for the parser @>=
expr_p make_while_node(expr_p b, unsigned flags, const YYLTYPE& loc)
{
  expr_ptr bb(b); expr& body=*bb;
  return new expr(new while_node { std::move(body), flags},loc) ;
}
@)
expr_p make_for_node
  (raw_id_pat& id, expr_p ip, expr_p b, unsigned flags, const YYLTYPE& loc)
{
  id_pat ind(id); expr_ptr iip(ip), bb(b);
  expr& in=*iip;  expr& body=*bb;
  return new expr(new @|
    for_node { std::move(ind), std::move(in), std::move(body), flags },loc);
}
@)
expr_p make_cfor_node
  (id_type id, expr_p c, expr_p l, expr_p b, unsigned flags,
   const YYLTYPE& loc)
{
  expr_ptr cc(c), ll(l), bb(b);
  expr& cnt=*cc; expr& lim=*ll; expr& body=*bb;
  return new expr (new @|
    cfor_node { id, std::move(cnt),std::move(lim),std::move(body),flags },loc);
}

@ Again we apply assignment for raw pointer variants.
@< Cases for copying... @>=
  case while_expr: while_variant=other.while_variant; break;
  case for_expr:   for_variant=other.for_variant;     break;
  case cfor_expr:  cfor_variant=other.cfor_variant;   break;

@~Cleaning up is exactly like that for conditional expressions, or indeed most
other kinds.

@< Cases for destroying... @>=
case while_expr: delete while_variant; break;
case for_expr:   delete for_variant;   break;
case cfor_expr:  delete cfor_variant;  break;

@ To print a |while| or |for| expression at parser level, we reproduce the
input syntax.

@< Cases for printing... @>=
case while_expr:
{ const w_loop& w=e.while_variant;
@/ out << " while " << w->body << " od ";
}
break;
case for_expr:
{ const f_loop& f=e.for_variant;
  const patlist& pl = f->id.sublist;
@/const id_pat& index = pl.front();
  const id_pat& entry = *++pl.begin(); // |pl.size()==2|
@/out << " for " << entry;
  if ((index.kind&0x1)!=0)
    out << '@@' << index;
  out << " in " << f->in_part
      << (f->flags[0] ? " ~do " : " do ") @| << f->body
      << (f->flags[1] ? " ~od " : " od ");
}
break;
case cfor_expr:
{ const cfor_node& c=*e.cfor_variant;
@/out << " for ";
  if (c.id!=id_type(-1))
    out << main_hash_table->name_of(c.id);
  out << ": " << c.count;
  if (not is_empty(c.bound))
    out << " from " << c.bound;
  out << (c.flags[0] ? " ~do " : " do ") << c.body
      << (c.flags[1] ? " ~od " : " od ");
}
break;

@*1 Array subscriptions and slices.
%
We want to be able to select components from array structures (lists, vectors,
matrices), so we define a subscription expression. If there are multiple
indices, as in selecting a matrix entry, these can be realised as a
subscription by a tuple expression, so we define only one type of subscription
expression.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct subscription_node* sub;
typedef struct slice_node* slc;

@~In a subscription the array and the index(es) can syntactically be arbitrary
expressions (although the latter should have as type integer, or a tuple of
integers).

@< Structure and typedef declarations for types built upon |expr| @>=
struct subscription_node
{ expr array; expr index; bool reversed;
@)
  subscription_node(expr&& array, expr&& index, bool reversed)
@/: array(std::move(array))
  , index(std::move(index))
  , reversed(reversed)@+{}
};
struct slice_node
{ expr array; expr lower,upper; BitSet<3> flags;
@)
  slice_node(expr&& array, expr&& lower, expr&& upper, unsigned flags)
@/: array(std::move(array))
  , lower(std::move(lower))
  , upper(std::move(upper))
  , flags(flags)
  @+{}
};

@ Here is the tag used for subscriptions and slices.

@< Enumeration tags for |expr_kind| @>= subscription, slice, @[@]

@~And here is the corresponding variant of |expr|.

@< Variants of ... @>=
sub subscription_variant;
slc slice_variant;

@ There are constructors for building the new variants. For the furst one a
variation will be useful too.
@< Methods of |expr| @>=
expr(sub&& s, const YYLTYPE& loc)
 : kind(subscription)
 , subscription_variant(std::move(s))
 , loc(loc)
 @+{}
expr(sub&& s, const source_location& loc)
 : kind(subscription)
 , subscription_variant(std::move(s))
 , loc(loc)
 @+{}
expr(slc&& s, const YYLTYPE& loc)
 : kind(slice)
 , slice_variant(std::move(s))
 , loc(loc)
 @+{}

@ To build an |subscription_node| or |slice_node|, we simply combine the array
and the index part.

@< Declarations of functions for the parser @>=
expr_p make_subscription_node
  (expr_p a, expr_p i, bool reversed, const YYLTYPE& loc);
expr_p make_slice_node
  (expr_p a, expr_p lower, expr_p upper, unsigned flags, const YYLTYPE& loc);

@~This is straightforward, as usual.

@< Definitions of functions for the parser @>=
expr_p make_subscription_node
  (expr_p a, expr_p i, bool reversed, const YYLTYPE& loc)
{ expr_ptr aa(a); expr_ptr ii(i);
  expr& arr=*aa; expr& ind=*ii;
  return new expr(new
     subscription_node {std::move(arr), std::move(ind), reversed },loc);
}
expr_p make_slice_node
  (expr_p a, expr_p lower, expr_p upper, unsigned flags, const YYLTYPE& loc)
{ expr_ptr aa(a); expr_ptr ll(lower); expr_ptr uu(upper);
@/expr& arr=*aa; expr& lo=*ll; expr& up=*uu;
  return new expr(new @|
     slice_node {std::move(arr), std::move(lo), std::move(up), flags },loc);
}

@ Another boring but inevitable section.
@< Cases for copying... @>=
case subscription: subscription_variant=other.subscription_variant; break;
case slice:        slice_variant       =other.slice_variant;        break;

@~Here we recursively destroy the subexpressions, and then the node for the
subscription or slice itself.

@< Cases for destroying... @>=
case subscription: delete subscription_variant; break;
case slice: delete slice_variant; break;

@ To print a subscription, we just print the expression of the array, followed
by the expression for the index in brackets. As an exception, the case of a
tuple display as index is handled separately, in order to avoid having
parentheses directly inside the brackets.

For slices, we print the expression with optional occurrences of \.{\char\~},
governed by |flags|.
@h "lexer.h"
@< Cases for printing... @>=
case subscription:
{ const sub& s=e.subscription_variant;
  const expr& i=s->index;
  out << s->array << '[';
  if (i.kind!=tuple_display) out << i;
  else
  {
    wel_const_iterator it (i.sublist);
    out << *it;
    while (not (++it).at_end())
      out << ',' << *it;
 }
  out << ']';
}
break;
case slice:
{ const slc& s=e.slice_variant;
@/out << s->array << (s->flags[0] ? "~[" : "[") @|
      << s->lower << (s->flags[1] ? "~:" : ":") @|
      << s->upper << (s->flags[2] ? "~]" : "]");
}
break;

@*1 Cast expressions.
%S
These are a way to force an expression to get a specified type (provided that
type analysis succeeds with that target type).

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct cast_node* cast;

@~The corresponding node type stores both type and expression.

@< Structure and typedef declarations for types built upon |expr| @>=
struct cast_node
{ type_expr type; expr exp;
@)
  cast_node(type_expr&& type, expr&& exp)
@/: type(std::move(type))
  , exp(std::move(exp))@+{}
  // backward compatibility for gcc 4.6
};

@ The tag used for casts is |cast_expr|.

@< Enumeration tags for |expr_kind| @>=
cast_expr, @[@]

@ And there is of course a variant of |expr_union| for casts.
@< Variants of ... @>=
cast cast_variant;

@ There is a constructor for building the new variant.
@< Methods of |expr| @>=
expr(cast&& c, const YYLTYPE& loc)
 : kind(cast_expr)
 , cast_variant(std::move(c))
 , loc(loc)
@+{}

@ Casts are built by |make_cast|.

@< Declarations of functions for the parser @>=
expr_p make_cast(type_p type, expr_p exp, const YYLTYPE& loc);

@~No surprises here.

@< Definitions of functions for the parser@>=
expr_p make_cast(type_p t, expr_p e, const YYLTYPE& loc)
{
  type_ptr tt(t); expr_ptr ee(e);
  type_expr& type=*tt; expr& exp=*ee;
  return new expr(new cast_node { std::move(type), std::move(exp) },loc);
}

@ Nor here.
@< Cases for copying... @>=
case cast_expr: cast_variant=other.cast_variant;
break;

@ Eventually we want to rid ourselves from the cast.

@< Cases for destroying... @>=
case cast_expr: delete cast_variant; break;
@ Printing cast expressions follows their input syntax.

@< Cases for printing... @>=
case cast_expr:
{@; const cast& c = e.cast_variant;
  out << c->type << ':' << c->exp ;
}
break;

@ A different kind of cast serves to obtain the current value of an overloaded
operator symbol.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct op_cast_node* op_cast;

@~We store an (operator) identifier and a type, as before represented by a
void pointer.

@< Structure and typedef declarations for types built upon |expr| @>=
struct op_cast_node
{ id_type oper; type_expr type;
@)
  op_cast_node(id_type oper,type_expr&& type)
@/: oper(oper)
  , type(std::move(type))@+{}
  // backward compatibility for gcc 4.6
};

@ The tag used for casts is |op_cast_expr|.

@< Enumeration tags for |expr_kind| @>=
op_cast_expr, @[@]

@ And there is of course a variant of |expr_union| for casts.
@< Variants of ... @>=
op_cast op_cast_variant;

@ There is a constructor for building the new variant.
@< Methods of |expr| @>=
expr(op_cast&& c, const YYLTYPE& loc)
 : kind(op_cast_expr)
 , op_cast_variant(std::move(c))
 , loc(loc)
@+{}

@ Operator casts are built by |make_op_cast|.

@< Declarations of functions for the parser @>=
expr_p make_op_cast(id_type name,type_p type, const YYLTYPE& loc);

@~No surprises here either.

@< Definitions of functions for the parser@>=
expr_p make_op_cast(id_type name,type_p t, const YYLTYPE& loc)
{
  type_ptr tt(t); type_expr& type=*tt;
  return new expr(new op_cast_node { name, std::move(type) },loc);
}

@ Nor here.
@< Cases for copying... @>=
case op_cast_expr: op_cast_variant=other.op_cast_variant; break;

@ Eventually we want to rid ourselves from the operator cast.

@< Cases for destroying... @>=
case op_cast_expr: delete op_cast_variant; break;

@ Printing operator cast expressions follows their input syntax.

@< Cases for printing... @>=
case op_cast_expr:
{ const op_cast& c = e.op_cast_variant;
  out << main_hash_table->name_of(c->oper) << '@@' << c->type;
}
break;

@*1 Assignment statements.
%
Simple assignment statements are quite simple as expressions.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct assignment_node* assignment;

@~In a simple assignment the left hand side is just an identifier. However, we
now also cater for general identifier patterns as left hand side, but keep
having a special constructor for the single identifier case.

@< Structure and typedef declarations for types built upon |expr| @>=
struct assignment_node
{ id_pat lhs; expr rhs;
@)
  assignment_node(id_type lhs, expr&& rhs)
@/: lhs(lhs)
  , rhs(std::move(rhs))@+{}
  assignment_node(id_pat&& lhs, expr&& rhs)
@/: lhs(std::move(lhs))
  , rhs(std::move(rhs))@+{}
#ifdef incompletecpp11
  assignment_node(const assignment_node& x) = @[delete@];
  assignment_node(assignment_node&& x)
@/: lhs(std::move(x.lhs))
  , rhs(std::move(x.rhs))
  @+{}
#endif
};

@ The tag used for assignment statements is |ass_stat|.

@< Enumeration tags for |expr_kind| @>=
ass_stat, @[@]

@ And there is of course a variant of |expr_union| for assignments.
@< Variants of ... @>=
assignment assign_variant;

@ As always there is a constructor for building the new variant.
@< Methods of |expr| @>=
expr(assignment&& a, const YYLTYPE& loc)
 : kind(ass_stat)
 , assign_variant(std::move(a))
 , loc(loc)
@+{}

@ Assignment statements are built by |make_assignment| (when simply assigning
to an identifier) or by |make_multi_assignment| (for a more general pattern as
left hand side, which requires using a different syntax).

@< Declarations of functions for the parser @>=
expr_p make_assignment(id_type lhs, expr_p rhs, const YYLTYPE& loc);
expr_p make_multi_assignment(raw_id_pat& lhs, expr_p rhs, const YYLTYPE& loc);

@~These functions do what one would expect them to (except for those who
expect them to make their homework assignments).

@< Definitions of functions for the parser@>=

expr_p make_assignment(id_type id, expr_p r, const YYLTYPE& loc)
{
  expr_ptr rr(r); expr& rhs=*rr;
  return new expr(new @| assignment_node { id, std::move(rhs) },loc);
}
@)
expr_p make_multi_assignment(raw_id_pat& lhs, expr_p r, const YYLTYPE& loc)
{
  expr_ptr rr(r); expr& rhs=*rr;
  return new expr(new assignment_node { id_pat(lhs), std::move(rhs) },loc);
}

@ Copying is done by assignment of raw pointers, not really a novelty.

@< Cases for copying... @>=
case ass_stat: assign_variant=other.assign_variant; break;

@ What is made must eventually be unmade (even assignments).

@< Cases for destroying... @>=
case ass_stat: delete assign_variant; break;

@ Printing assignment statements, we include |"set"| only if this cannot be a
simple assignment

@< Cases for printing... @>=
case ass_stat:
{ const assignment& ass = e.assign_variant;
  if ((ass->lhs.kind&0x3)!=0x1)
    out << "set ";
  out << ass->lhs << ":=" << ass->rhs ;
}
break;

@*2 Component and field assignments.
%
We have special expressions for assignments to an indexed component of an
indexable value (like a row, vector, or matrix) that is bound to an
identifier, and for assignments to a named field of a tuple that is bound to
an identifier. The |typedef| names are shortened here to avoid a name conflict
with types defined in \.{axis.w}.

@< Type declarations needed in definition of |struct expr@;| @>=
typedef struct comp_assignment_node* comp_assignment;
typedef struct field_assignment_node* fld_assignment;

@~In a component assignment, the left hand side records an identifier (holding
the whole value) and an index. In addition the user may have specified reverse
indexing as in $a\sim[0]:=x$, and whether this is the case will be recorded in
the |reversed| field.

In a field assignment, the the left hand side records two identifiers, one
for the identifier holding the tuple, and one for the selection function
indicating the field to modify.

@< Structure and typedef declarations for types built upon |expr| @>=
struct comp_assignment_node
{ id_type aggr; expr index; expr rhs; bool reversed;
@)
  comp_assignment_node(id_type aggr, expr&& index, expr&& rhs, bool reversed)
@/: aggr(aggr)
  , index(std::move(index))
  , rhs(std::move(rhs))
  , reversed(reversed)
  @+{}
#ifdef incompletecpp11
  comp_assignment_node(const comp_assignment_node& x) = @[delete@];
  comp_assignment_node(comp_assignment_node&& x)
@/: aggr(x.aggr)
  , index(std::move(x.index))
  , rhs(std::move(x.rhs))
  , reversed(x.reversed)
  @+{}
#endif
};
@)
struct field_assignment_node
{ id_type aggr, selector; expr rhs;
@)
  field_assignment_node(id_type aggr, id_type selector, expr&& rhs)
@/: aggr(aggr)
  , selector(selector)
  , rhs(std::move(rhs))
  @+{}
#ifdef incompletecpp11
  field_assignment_node(const field_assignment_node& x) = @[delete@];
  field_assignment_node(field_assignment_node&& x)
@/: aggr(x.aggr)
  , selector(selector)
  , rhs(std::move(x.rhs))
  @+{}
#endif
};

@ The tags used for component assignment statements are |comp_ass_stat| for
the indexed case, and |select_ass_stat| for the field-selected case.

@< Enumeration tags for |expr_kind| @>=
comp_ass_stat, field_ass_stat, @[@]

@ And there are of course variants of |expr_union| for these component
assignments.
@< Variants of ... @>=
comp_assignment comp_assign_variant;
fld_assignment field_assign_variant;

@ As always there are constructors for building the new variants.
@< Methods of |expr| @>=
expr(comp_assignment ca, const YYLTYPE& loc)
 : kind(comp_ass_stat)
 , comp_assign_variant(ca)
 , loc(loc)
@+{}
@)
expr(fld_assignment&& fa, const YYLTYPE& loc)
 : kind(field_ass_stat)
 , field_assign_variant(std::move(fa))
 , loc(loc)
@+{}

@ Compound assignment statements are built by |make_comp_ass|, which for once
does not simply combine the expression components, because for reason of
parser generation the array and index will have already been combined before
this function can be called. There is also |make_comp_upd_ass| for the
operate-assign form of component assignments, as in $v[i]\mathrel+:=a$.
Finally, similar functions are also given for field assignments. The plain
assignment just needs the identifiers, but the operate-assign forms wraps them
into |applied_identifier| expressions (which attaches the proper location
information) since these are necessary in the right hand side to be
constructed.

@< Declarations of functions for the parser @>=
expr_p make_comp_ass(expr_p lhs, expr_p rhs, const YYLTYPE& loc);
expr_p make_comp_upd_ass(expr_p lhs, id_type op, expr_p rhs,
   const YYLTYPE& loc,  const YYLTYPE& op_loc);
@)
expr_p make_field_ass(id_type tup, id_type sel, expr_p rhs, const YYLTYPE& loc);
expr_p make_field_upd_ass(expr_p tupex, expr_p selex, id_type op, expr_p rhs,
   const YYLTYPE& loc,  const YYLTYPE& op_loc);

@ Here we have to take the left hand side apart a bit; the grammar ensures
that |lhs| presents the necessary structure for this code to work.

@< Definitions of functions for the parser@>=
expr_p make_comp_ass(expr_p l, expr_p r, const YYLTYPE& loc)
{
  expr_ptr ll(l), rr(r);
  assert(ll->kind==subscription); // grammar ensures this
  subscription_node& s = *ll->subscription_variant;
  assert(s.array.kind==applied_identifier); // grammar ensures this
  return new expr(new comp_assignment_node @|
  { s.array.identifier_variant
@|, std::move(s.index)
  , std::move(*rr)
  , s.reversed
  },loc);
}

@ Here is the even simpler code for field assignments.

@< Definitions of functions for the parser@>=
expr_p make_field_ass(id_type tup, id_type sel, expr_p r, const YYLTYPE& loc)
{@;
  expr_ptr rr(r);
  return new expr(new field_assignment_node { tup , sel, std::move(*rr) },loc);
}


@ This variation is somewhat complicated both at input (for the same reason as
|make_comp_ass|) and the output side, because $v[I]\mathrel\star:= E$ is
translated into essentially ``\&{let}~$\$=I$~\&{in}~$v[\$]:=v[\$]\star E$'',
where $\$$ is a local variable (to avoid evaluating $I$ twice) that
can't conflict with $v$ or any names used in the expression~$E$.

In doing everything by hand (except for some constructions done by
|make_binary_call|), this code gives an idea of the amount of work that the
parser actions go through in handling even relatively simple expressions.

@< Definitions of functions for the parser@>=
expr_p make_comp_upd_ass(expr_p l, id_type op, expr_p r,
   const YYLTYPE& loc, const YYLTYPE& op_loc)
{
  expr_ptr ll(l), rr(r);
  assert(ll->kind==subscription); // grammar ensures this
  auto& s = *ll->subscription_variant;
  assert(s.array.kind==applied_identifier); // likewise
  id_type array_name=s.array.identifier_variant;
      // save before move from |s.array|
  id_type hidden=lookup_identifier("$");@q$@>
  id_pat v(hidden);
  expr::identifier_tag id_tag;
  expr_ptr subs(new expr (new subscription_node@|
    (std::move(s.array)
    ,expr(hidden,loc,id_tag)
    ,s.reversed
    ),ll->loc));
  expr_ptr formula(make_binary_call(op,subs.release(),rr.release(),loc,op_loc));
  expr ca (new comp_assignment_node @|
     {array_name
     ,expr(hidden,loc,id_tag)
     ,std::move(*formula)
     ,s.reversed},loc);
  return new expr(new let_expr_node @|
    (std::move(v),std::move(s.index), std::move(ca)),loc);
}

@ Here is the corresponding code for field assignments.

@< Definitions of functions for the parser@>=
expr_p make_field_upd_ass(expr_p tupex, expr_p selex, id_type op, expr_p r,
   const YYLTYPE& loc,  const YYLTYPE& op_loc)
{
  expr_ptr rr(r);
  expr_ptr tup_exp(tupex);
  expr_ptr sel_exp(selex);
  assert(tup_exp->kind==applied_identifier and
         sel_exp->kind==applied_identifier); // grammar
  id_type tup=tup_exp->identifier_variant, sel=sel_exp->identifier_variant;
   // save, then move
  expr_ptr selec(new expr (new application_node@|
    (std::move(*sel_exp),std::move(*tup_exp)),loc));
  expr_ptr formula
    (make_binary_call(op,selec.release(),rr.release(),loc,op_loc));
  return new expr(new field_assignment_node(tup,sel,std::move(*formula)),loc);
}

@ This is getting boring; fortunately we are almost done with the syntax.

@< Cases for copying... @>=
case comp_ass_stat: comp_assign_variant=other.comp_assign_variant; break;
case field_ass_stat: field_assign_variant=other.field_assign_variant; break;

@~Destruction one the other hand is as straightforward as usual.

@< Cases for destroying... @>=
case comp_ass_stat: delete comp_assign_variant; break;
case field_ass_stat: delete field_assign_variant; break;

@ Printing component assignment statements follow the input syntax.

@< Cases for printing... @>=
case comp_ass_stat:
{ const comp_assignment& ass = e.comp_assign_variant;
  out << main_hash_table->name_of(ass->aggr)
      << (ass->reversed ? "~[" : "[") << ass->index << @| "]:="
      << ass->rhs ;
}
break;
case field_ass_stat:
{ const fld_assignment& ass = e.field_assign_variant;
  out << main_hash_table->name_of(ass->aggr)
      << '.' << main_hash_table->name_of(ass->selector) @|
      << ":=" << ass->rhs ;
}
break;

@*1 Recursive function expressions.
%
The basic \.{axis} language is not friendly for recursion, since identifiers
defined in a local or global definition only come into scope after the body of
the definition. The fact that recursive functions can nonetheless be defined
is due to the possibility to call functions from a variable, where a runtime
assignment to the variable ensures that the by the time it gets called, it
refers to the very function (body) that contains the call. We provide
syntactic sugar to make this easier for the user; this being so, no new kind
of |expr| is needed to implement it.

@< Declarations of functions for the parser @>=
expr_p make_recfun(id_type f, expr_p d,
   const YYLTYPE& loc, const YYLTYPE& f_loc);

@ The \&{rec\_fun} syntax requires a function name $f$ and a definition |def|
which, apart from specifying its argument type(s)~$A$ as usual, also specifies
an explicit result type~$R$; this is represented as having a body that is a
cast to~$R$. We shall translate this into ``\&{let} $f=(A~.)
R:$~\&{die} \&{in} $(f:=\\{def})$'' where the dot designates an empty
identifier pattern (the \&{die} must be wrapped in a function to prevent the
evaluation of \&{let} from blowing up).

@< Definitions of functions for the parser@>=
expr_p make_recfun(id_type f, expr_p d,
   const YYLTYPE& loc, const YYLTYPE& f_loc)
{
  expr_ptr dd(d); expr& definition=*dd;
  assert(definition.kind==lambda_expr); // grammar ensures this
  lambda_node& lam = *definition.lambda_variant;
  assert(lam.body.kind==cast_expr); // grammar ensures this
  cast_node& body = *lam.body.cast_variant;
  expr die_expr(f_loc,expr::die_tag());
  expr dummy_body (new cast_node(body.type.copy(),std::move(die_expr)),f_loc);
  expr dummy_f
   (new lambda_node@|(id_pat(),lam.parameter_type.copy(),std::move(dummy_body))
    ,loc);
  expr rec_assign (new assignment_node(f,std::move(definition)),loc);
  return new expr(new let_expr_node @|
    (id_pat(f),std::move(dummy_f),std::move(rec_assign)),loc);
}


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

@< Structure and typedef declarations for types built upon |expr| @>=
struct sequence_node
{ expr first; expr last;
@)
  sequence_node(expr&& first, expr&& last)
@/: first(std::move(first))
  , last(std::move(last))@+{}
  // backward compatibility for gcc 4.6
};

@ In fact we use this expression also for two variants variant of sequence
expressions, one in which the value of the \emph{first} expression is retained
as final value and the second expression is then evaluated without using its
value (the \&{next} keyword), and one where the result of the left operand
determines whether the right operand will be evaluated or not (the |do| in
|while| loops).

@< Enumeration tags for |expr_kind| @>=
seq_expr, next_expr, do_expr, @[@]

@ And there is of course a variant of |expr_union| for sequences.
@< Variants of ... @>=
sequence sequence_variant;

@ As always there is a constructor for building the new variant.
@< Methods of |expr| @>=
expr(sequence&& s, unsigned which, const YYLTYPE& loc)
 : kind(which==0 ? seq_expr : which==1 ? next_expr : do_expr)
 , sequence_variant(std::move(s))
 , loc(loc)
@+{}

@ Sequences are built by |make_sequence|.

@< Declarations of functions for the parser @>=
expr_p make_sequence
  (expr_p first, expr_p last, unsigned which, const YYLTYPE& loc);

@~It does what one would expect it to.

@< Definitions of functions for the parser @>=
expr_p make_sequence
  (expr_p f, expr_p l, unsigned which, const YYLTYPE& loc)
{
  expr_ptr ff(f), ll(l); expr& first=*ff; expr& last=*ll;
  return new expr(new @|
    sequence_node { std::move(first), std::move(last) }, which ,loc);
}

@ Are these the final cases? For now, they are!

@< Cases for copying... @>=
case seq_expr:
case next_expr:
case do_expr:
   sequence_variant=other.sequence_variant; break;

@ Finally sequence nodes need destruction, like everything else.

@< Cases for destroying... @>=
case seq_expr:
case next_expr:
case do_expr:
   delete sequence_variant; break;

@ Printing sequences is absolutely straightforward.

@< Cases for printing... @>=
case seq_expr:
{@; const sequence& seq = e.sequence_variant;
  out << seq->first << "; " << seq->last ;
}
break;
case next_expr:
{@; const sequence& seq = e.sequence_variant;
  out << seq->first << " next " << seq->last ;
}
break;
case do_expr:
{@; const sequence& seq = e.sequence_variant;
  out << seq->first << " do " << seq->last ;
}
break;

@* Other functions callable from the parser.
Here are some functions that are not so much a parsing functions as just
wrapper functions that used to enable the parser to call \Cpp~functions.

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
