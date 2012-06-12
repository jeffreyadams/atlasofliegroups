@* Parsing the program code. @:title@>
The most intricate part of \.{\me.} is its mechanism for converting
\Cee-like code into \TeX\ code, and we shall consider this aspect of the
program now. This parsing mechanism must be able to deal with fragmentary
constructions whose overall ``part of speech'' is not known (e.g., those
within~`\pb'), which implies that recognition must in principle be able to
proceed without any context information. Therefore a ``bottom-up'' approach
is used, that collects successive tokens, and only takes action once it has
seen enough to recognise a complete match with a syntax rule. Bottom-up
parsing is a very powerful technique (it can handle grammars that are not
tractable by many other parsing methods); it is one of the many significant
contributions to computer science that the designer of the original \.{WEB}
system, D.~E. Knuth, has made. Even so, the technique used here is less
powerful than that of traditional bottom-up parsers produced by parser
generators, since those can use information derived from the full
left-context (and a small amount of right-context) to help decide difficult
situations. In practice such context is not often required to unambiguously
parse program fragments, and most situations where there is a possibility of
incorrect recognition can be avoided be careful specification of the syntax
rules to be used; this does however make formulation of the grammar a
slightly subtle matter.

As we have already seen, the input is represented as a sequence of
{\sl scraps}, each of which specifies a {\sl category\/} and a
{\sl translation}. The category defines a syntactic class, and the
translation is a token list that represents \TeX\ code; the latter is
effectively just a single string of characters, although it is stored in a
different form for efficiency reasons. Rules of syntax and semantics tell us
how to combine adjacent scraps into larger ones, and if we are lucky an
entire \Cee\ text that starts out as hundreds of small scraps will join
together into one gigantic scrap whose translation is the desired
\TeX\ code. If we are unlucky (i.e., if the input is not properly formed
for recognition by our parser), we will be left with several scraps that
don't combine; their translations will simply be output, one by one.

@*1  Introduction to the grammar.
The parsing process is governed by a set of reduction rules, each of which
specifies how a sequence of consecutive scraps with certain given categories
can be combined into a single new scrap, and it also tells what category the
new scrap has and how its translation is obtained from those of the original
scraps.  The set of reduction rules is given by a table, which is actually
an initialiser for a static array, starting in section@#rules@>.
In many cases a rule is simply given by listing the category codes
forming its pattern (its left hand side) and the category code of the
resulting scrap (the right hand side). For example, one of the reduction
rules could be represented as
$$
\hbox{|expression| |binop| |expression|}\Rightarrow\hbox{|expression|}.
$$
For such simple rules the translations of the original scraps are simply
concatenated to form the translation of the resulting scrap. The line in the
table corresponding to this rule actually looks like this:
$$
\hbox{|{2,{{expression,binop,expression}},{expression, NULL}}|}.
$$
The `2' is the number used to identify the rule, which is useful when
tracing the activities of the parser, and the |NULL| signifies that the
default translation applies; the braces reflect the structure used to
represent syntax rules, which has some additional components that are not
specified explicitly here, and that are therefore initialised to~0.

For other rules specific additional items are to be inserted in between the
translations of the scraps that match the pattern. For instance, here is
another rule present in the table
$$
\hbox{|{ 6,{{expression, comma, expression}},{expression,"__p1_"}}|}.
$$
Here the translation is specified by the ``format string'' |"__p1_"|, where
the underscores stand for the translations of the three participating
scraps, the `\.{p1}' for an |opt| token with argument~`\.1' (eventually
producing `\.{\\penalty10}') that comes after the translation of the second
scrap, with category |comma|. That translation will almost always be |","|,
since there are no rules with |comma| as their right hand side, but if a
comment follows the comma, it will have been incorporated into its scrap
while preparing for the reduction, and in this case the optional break will
follow the comment (and have no effect, because there is already a forced
break after comments). The interpretation of all characters that can appear
in the format strings can be found in section@#format@>.

To handle all cases that can arise in programs, in combination with a large
number of formatting options that might be selected, our mechanism for rules
of a somewhat more general form than those presented above, but the
extensions are used only for a few rules. Consider the following rule:
$$
\hbox{|{182,{{lpar,statement,statement},1},
			{statement,"_B_"},forced_statements}|}.
$$
Here, the `1' indicates that the first scrap from the left in the pattern,
|lpar|, serves as context only, and does not take part in the actual
reduction. The rule therefore reduces a sequence of two consecutive
statements into a single one, inserting a |break_space| (from the `\.B'
format character) between their translations, but only if the statements are
preceded by a left parenthesis (category~|lpar|); this restriction means
that the rule will only apply to the first two semicolon-separated
expressions following `|for|~|(|'. The entry |forced_statements| means
moreover that the rule will not be used at all unless a `\.{+f}'~or~`\.{+a}'
command line option was specified, which forces line breaks between
consecutive statements; this rule then avoids such line breaks between the
controlling expressions of a |for|-loop, since it takes precedence (see
below) over the rule that normally combines statements.

@ The rules are applied from left to right, as follows. Suppose we are
currently working on the sequence of scraps with categories $c_1\,c_2\ldots
c_n$. We try first to find a rule whose pattern matches an initial substring
$c_1\,c_2\ldots\,$, and if no such rule exists, we try to find a rule
applicable to the next substring $c_2\,c_3\ldots\,$; if that fails, we try
to match $c_3\,c_4\ldots\,$, etc. When several patterns match, starting at
the same scrap, the longest one is selected (no two rules can have identical
patterns). For instance, there is a rule that reduces `|struct_like|
|expression| |lbrace|' to `|struct_head|', and another one that reduces
`|struct_like| |expression|' to `|int_like|'; the latter is only applied
when no |lbrace| scrap follows. This is a sensible rule, since the longer
pattern is the more specific one, and without such a rule it could never
match. One should be aware however that this only works because the scrap
with category |lbrace| represents a single token that requires no reduction
to create it, for otherwise the two-scrap reduction would be applied before
the three-scrap reduction would have a chance to match.

One might say that we use a left-to-right eager strategy for choosing
reductions; this strategy is chosen on heuristic grounds, and there is no
guarantee that it will find a successful sequence of reductions if one
exists. In other words, if we interchange left and right hand sides of the
rules and view them as {\sl production\/} rules of a context free grammar,
then the language generated by this grammar can be much larger than the one
recognised by our parser. For instance, the context free grammar would, from
the rules given above, generate expressions with multiple operators where
implicit parentheses are placed arbitrarily (and therefore it would be
ambiguous); our parsing strategy however will always place implicit
parentheses to the left. (We don't care if implicit parentheses are placed
incorrectly, because it does not influence typesetting.) It is an
interesting theoretic problem to find an algorithm that will transform a
set of reductions rules into a context free grammar generating exactly the
language recognised; that would help to verify that the rules will work as
intended. We must be cautious in formulating the rules, not only to give
rules that are sufficient to reduce all desired programs, but also not to
specify rules that could in certain circumstances match with precedence over
the ones we intended for the situation.

Some of the formatting of the output is controlled by the \TeX~format in the
file \.{cwebxmac.tex}, rather than directly by \.{\me.}. For instance, the
`\.{\\penalty10}' mentioned above is actually written as `\.{\\31}'; the
macro `\.{\\3}' will in fact raise the penalty for any enclosing pair of
braces or parentheses; thus breaks in argument lists or initialiser lists
will be avoided in favour of ones outside those lists. There are other
effects not visible in the grammar, like optional breaks automatically
associated with certain operators. In the whole it is a delicate interplay
between the set of reduction rules, the reduction mechanism, the macro
definitions, and the (mathematical) typesetting rules of \TeX\ that
determines how programs will be typeset. This division of labour is quite
convenient (e.g., it would tremendously complicate the task of \.{\me.} if
it had to decide which optional breaks should actually be taken), but it has
the unfortunate consequence that it is not easy to gain a comprehensive
understanding of the translation process. Indeed there have been quite a few
surprises during the construction of the system, and we cannot be sure that
there are no more to come; yet we believe the situation is more transparent
than in the original |WEAVE|, since we have avoided rules that are oriented
towards the target language (\TeX) rather than towards the source language
(\Cee).

@*2 Categories of scraps.
Here is a list of the category codes that scraps can have. The category codes
that apply to reserved words (e.g., |while_like|, but also |declaration| for
|va_dcl|) as well as |expression| (that is used for |type_defined| identifiers
in their typedef declaration) are also used as their |ilk| value, and they are
sufficiently high that they can be distinguished from |ilk| values that are
not category codes, like |type_defined|, |TeX_like|, |NULL_like|,
|const_like|, |typedef_like| and |namespace_like|; thus reserved words will
never be mistaken for identifiers with such an ilk, and we can safely
substitute certain such |ilk| values by a different category code as happens
in section@#ilk change@>. A number of categories can occur in scraps, but do
not occur in any of the reduction rules, since they are handled by other
means; they have values $\geq{}$|min_nonrule_cat|. The macro |valid_cat|
checks whether |c| is a category that might match in a rule; it uses its
argument twice, so its argument should not cause a side effect. If this
section is changed, section@#category output@> should be changed
correspondingly.

@d max_category end_expr /* largest scrap category */
@d min_nonrule_cat lproc /* first category forbidden in rules */
@d valid_cat(c) ((c)>0 && (c)<min_nonrule_cat)

@<Typedef and enum...@>=
enum @/ @:categories@>
{ unop = 1, /* a unary operator like `|!|' */
  binop, /* a binary operator like `|<|' */
  unorbinop, /* an operator that can be either, like `|-|' */
  select, /* structure selection: `|.|' or `|->|' */
  question, /* a question mark operator */
  lbrace, rbrace, lpar, rpar, lbrack, rbrack,
    /* `|{|', \ `|}|', `(', \ `)', `[', \ `]' */
  comma, semi, colon, colcol, magic, /* `,', `;', `:', `$\CC$', \:; */
  subscript, /* an array subscript, like `|[]|' or `|[i++]|' */
  struct_head, /* the beginning of a struct specifier, like `|struct s{|' */
  short_lbrace, short_struct_head, /* from `\.{\{@@;}', for one-liners */
  compound_statement, /* a complete compound statement */
  statement, /* a complete statement, possibly compound */
  function, /* a complete function definition */
  function_head, /* a function identifier followed by formal parameters */
  parameters,
  /* parameters in function declaration, or casting operator like `|(int)|' */
  label, /* a statement label */
  if_head, /* `|if|' followed by a (parenthesised) expression */
  if_else_head, /* \lq|if @t\dots@>@; else|', \lq|while(@t\dots@>)|'
                or \lq|switch(@t\dots@>)|' */
  do_head, /* `|do @t\dots@>@; while|' */
  mod_scrap, /* module name */
  declarator, /* abstract declarator, like `|(*)(int,char*[])|' */
  declaration, /* a complete declaration */
  expression, /* an expression, possibly a single identifier */
  while_like, /* `|for|', `|while|', `|switch|' */
  do_like,
  if_like,
  else_like,
  extern_like,
  throw_like, try_like, catch_like, /* single tokens */
  int_like, /* `|int|', `|char|', \dots  */
  case_like, /* `|case|', `|default|' */
  sizeof_like, /* `|sizeof|', `|const|', `\&{new}', `\&{delete} */
  struct_like, /* `|struct|', `|union|', `|enum|',
                  `\&{class}', `\&{typename}' */
  return_like, /* `|return|', `|break|', `|continue|', `|goto|' */
  template_like, langle, rangle,
              /* `\&{template}', `$\langle$', `$\rangle$', for \Cpp */
  templ_params, /* template parameters */
  lproc, /* `\&\#' and following identifier starting preprocessor directive */
  rproc, /* end of a preprocessor directive */
  insert, /* comment or other syntactically inert item */
  begin_expr, end_expr /* \:[ and \:] */
};

@ As we have already seen, tokens are converted to elementary scraps by
the function |C_read|; these scraps form the `terminal symbols' of our
grammar. The translation of tokens to scraps is largely governed by the
static array |trans_ini|, whose initialisation values we shall now give.

@d yes_math 1 /* should be in math mode */
@d no_math 2 /* should be in horizontal mode */
@d maybe_math 0 /* works in either horizontal or math mode */

@< Initialiser for |trans_ini| @>=
{ '!', unop,	yes_math, "\\R" }, @/ @.\\R@>
{ '~', unop,	yes_math, "\\CM" }, @/ @.\\CM@>
{ '/', binop,	yes_math, "/" }, @/
{ '<', binop,	yes_math, "<" }, @/
{ '>', binop,	yes_math, ">" }, @/
{ '.', select,	yes_math, "." }, @/
{ '=', binop,	yes_math, "\\K" }, @/ @.\\K@>
{ '|', binop,	yes_math, "\\OR" }, @.\\OR@>
   { or, binop,	yes_math, "\\OR" }, @/
{ '^', binop,	yes_math, "\\XOR" }, @/ @.\\XOR@>
{ '%', binop,	yes_math, "\\MOD" }, @/ @.\\MOD@>
{ '+', unorbinop,	yes_math, "+" }, @/
{ '-', unorbinop,	yes_math, "-" }, @/
{ '*', unorbinop,	yes_math, "*" }, @/
{ '&', unorbinop,	yes_math, "\\AND" }, @/ @.\\AND@>
{ '?', question,	yes_math, "\\?" }, @/
{ '(', lpar,	yes_math, "(" }, @/
{ ')', rpar,	yes_math, ")" }, @/
{ '[', lbrack,	maybe_math, "[" }, @/
{ ']', rbrack,	maybe_math, "]" }, @/
{ '{', lbrace,	yes_math, "\\{" }, @/
{ '}', rbrace,	yes_math, "\\}" }, @/
{ ',', comma,	yes_math, "," }, @/
{ ';', semi,	yes_math, ";" }, @/
{ ':', colon,	maybe_math, ":" }, @/
{ '#', insert,	maybe_math, "\\#" },
    /* this should occur only in macro definitions */ @/ @.\\\#@>
{ at_sign_image, insert, maybe_math, "@@" },
    /* this should not occur in legal \Cee~text */
@)
{ not_eq,	binop,	yes_math, "\\I" }, @/ @.\\I@>
{ lt_eq,	binop,	yes_math, "\\Z" }, @/ @.\\Z@>
{ gt_eq,	binop,	yes_math, "\\G" }, @/ @.\\G@>
{ eq_eq,	binop,	yes_math, "\\E" }, @/ @.\\E@>
{ and_and,	binop,	yes_math, "\\W" }, @/ @.\\W@>
{ or_or,	binop,	yes_math, "\\V" }, @/ @.\\V@>
{ plus_plus,	unop,	yes_math, "\\PP" }, @/ @.\\PP@>
{ minus_minus,	unop, 	yes_math, "\\MM" }, @/ @.\\MM@>
{ minus_gt,	select,	yes_math, "\\MG" }, @/ @.\\MG@>
{ gt_gt,	binop,	yes_math, "\\GG" }, @/ @.\\GG@>
{ lt_lt,	binop,	yes_math, "\\LL" }, @/ @.\\LL@>
{ mul_assign,	binop,	yes_math, "\\KK*" }, @/ @.\\KK@>
{ div_assign,	binop,	yes_math, "\\KK/" }, @/
{ mod_assign,	binop,	yes_math, "\\KK\\MOD" }, @/ @.\\MOD@>
{ plus_assign,	binop,	yes_math, "\\KK+" }, @/
{ minus_assign,	binop,	yes_math, "\\KK-" }, @/
{ left_assign,	binop,	yes_math, "\\KK\\LL" }, @/
{ right_assign,	binop,	yes_math, "\\KK\\GG" }, @/
{ and_assign,	binop,	yes_math, "\\KK\\AND" }, @/ @.\\AND@>
{ xor_assign,	binop,	yes_math, "\\KK\\XOR" }, @/ @.\\XOR@>
{ or_assign,	binop,	yes_math, "\\KK\\OR" }, @/ @.\\OR@>
{ thin_space,	insert, yes_math, "\\," }, @/ @.\\,@>
{ pseudo_semi,	magic,	maybe_math, "" }, @/
{ force_expr_open,  begin_expr,	  maybe_math, "" }, @/
{ force_expr_close, end_expr,	  maybe_math, "" }, @/
{ join,		insert,	no_math, "\\J" }, @/ @.\\J@>
{ ellipsis,   int_like,	yes_math, "\\ldots" }, @/
{ sh_sh,	binop,	yes_math, "\\SS" }, @/ @.\\SS@>
{ colon_colon,	colcol, yes_math, "\\CC" }

@ For \Cpp\ we give `\.<' and `\.>' separate categories, which is necessary to
be able to recognise template arguments.

@< Fix the categories of angle brackets for \Cpp @>=
{@; token_trans['<'].cat=langle; token_trans['>'].cat=rangle; }

@ Certain tokens that lead to fixed scraps are not included in the
|trans_ini| array because their translations involve non-character tokens.
Since there are only a few of them the easiest solution is to install each
one explicitly into the |token_trans| array.
@d start_scrap(s,c,m) p=&token_trans[s],p->cat=c, p->mathness=5*(m)
@d end_scrap p->trans=text_ptr, freeze_text();

@< Install the translations of tokens involving line breaks @>=
{ scrap* p;
  start_scrap(math_break,insert,maybe_math); /* \:\v */
  app(opt), app('0'); end_scrap;
@/start_scrap(line_break,insert,no_math); /* \:/ */
  app(force); end_scrap;
@/start_scrap(end_preproc,rproc,no_math); /* end of preprocessor directive */
  app(force); end_scrap;
@/start_scrap(' ',insert,no_math); /* space within preprocessor directive */
  app(break_space); end_scrap;
@/start_scrap(big_line_break,insert,no_math); /* \:) */
  app(big_force); end_scrap;
@/start_scrap(backup_line,insert,no_math); /* \:\\ */
  app(backup); end_scrap;
@/start_scrap(no_line_break,insert,no_math); /* \:+ */
  app(cancel),app(relax),app(break_space),app(relax),app(cancel); end_scrap;
@/start_scrap(include_preproc,insert,yes_math); /* \:p */
  app(force),app_str("\\ATP"),app(force); end_scrap; @.\\ATP@>
}

@ When \.{\me.} is compiled with the |DEBUG| switch, it can display
its parsing steps. The order of strings in |cat_name| must match that
in the |enum| declaration in section@#categories@>.

@c
@:category output@>
#ifdef DEBUG
void print_cat (int c) /* symbolic printout of a category */
{ static char* cat_name[]=
  { "unop", "binop", "op", "select"
  , "?", "{", "}", "(", ")", "[", "]", ",", ";", ":", "::", "@@;"
  , "subscr", "struct_head", "short_{", "short_struct_head"
  , "cmp_stmt", "stmt"
  , "function", "function_head", "params", "label"
  , "if_head", "if_else_head", "do_head"
  , "mod_name", "declarator", "decl", "exp"
  , "for", "do", "if", "else", "extern"
  , "throw", "try", "catch,"
  , "int", "case", "sizeof", "struct", "return"
  , "template", "<", ">", "templ_params"
  , "#{", "#}", "insert", "@@[", "@@]"
  };
  if (c<=max_category && c>0) printf("%s",cat_name[c-1]);
  else printf ("IMPOSSIBLE");
}
#endif /* |DEBUG| */

@ Another major class of terminal symbols is formed by the reserved words.
If a name exists in the hash table with an |ilk| specifying a reserved word,
then |id_lookup| will return the reserved word when called with that name,
and |C_read| will use the |ilk| to set the category of the resulting scrap.
So all that has to be done is to get all the reserved words into the hash
table with the right ilks initially. The simplest way to do this is to call
|id_lookup| for all reserved words with the proper |ilk| at the beginning of
each run of \.{\me.}.  Fortunately there are not too many reserved words.
This code below uses the fact that instead of using pointers to beginning
and end of a string for |id_lookup|, one may also pass a single pointer to a
null-terminated string provided the other pointer is null.

@^reserved words@>

@<Store all the reserved words@>=
{ int i; static char* int_likes[]=
    { "auto","char","double","float","int","long","register"
    , "short","signed","static","unsigned","void" };
  static char* defined_types[] =
    { "FILE", "size_t", "ptrdiff_t", "wchar_t"
    , "jmp_buf", "sig_atomic_t", "fpos_t", "div_t", "ldiv_t"
    , "clock_t","time_t"
    , "va_list"
    };
  static char* return_likes[]=
    {"break","continue","goto","return"};
  int int_like_nr=array_size(int_likes),
      defined_type_nr=array_size(defined_types),
      return_like_nr=array_size(return_likes);

  for (i=0; i<int_like_nr; ++i) id_lookup(int_likes[i],NULL,int_like);
  for (i=0; i<defined_type_nr; ++i)
    id_lookup(defined_types[i],NULL,type_defined);
  for (i=0; i<return_like_nr; ++i) id_lookup(return_likes[i],NULL,return_like);

  id_lookup("case", NULL, case_like);
  id_lookup("const", NULL, const_like);
  id_lookup("default", NULL, case_like);
  id_lookup("do", NULL, do_like);
  id_lookup("else", NULL, else_like);
  id_lookup("enum", NULL, struct_like);
  id_lookup("extern", NULL, C_plus_plus ? extern_like : int_like);
  id_lookup("for", NULL, while_like);
  id_lookup("if", NULL, if_like);
  id_lookup("sizeof", NULL, sizeof_like);
  id_lookup("struct", NULL, struct_like);
  id_lookup("switch", NULL, while_like);
  id_lookup("typedef", NULL, typedef_like);
  id_lookup("union", NULL, struct_like);
  id_lookup("va_dcl",NULL, declaration);
  id_lookup("volatile", NULL, const_like);
  id_lookup("while", NULL, while_like);
  id_lookup("NULL", NULL, NULL_like);
  id_lookup("TeX", NULL, TeX_like);
  if (C_plus_plus) @<Store reserved words for \Cpp@>
}

@ One main difference between \Cee\ and \Cpp, as far as \.{\me.} is
concerned, is that the latter has a number of additional reserved words. Most
of them are sufficiently like some \Cee-reserved word (or category) that we
can simply make it behave like that \Cee~symbol, without changing the syntax.
For `\&{new}', `\&{delete}', and `\&{operator}', some additional syntax rules
will be needed however; nevertheless we do not need to extend the set of
syntactic categories. For `\&{operator}' we abuse the category |case_like|,
since its proper use is rather restricted (`|case|' it is always followed by
an expression, while `|default|', `\&{private}' and its relatives are always
followed by a colon), so there will be no confusion with `\&{operator}', which
is always followed by an operator symbol. Similarly, the different cast
operators are always followed by an operator (`$<$') which |sizeof| is never,
so they should not cause confusion.

@<Store reserved words...@>=
{ int i;
  static char* cpp_types[] =
  { "exception", "bad_exception", "bad_cast", "bad_typeid", "logic_error",
    "domain_error", "invalid_argument", "length_error", "out_of_range",
    "bad_alloc", "runtime_error", "range_error", "overflow_error",
    "underflow_error",
    "string",
    "iterator", "const_iterator", "reverse_iterator", "const_reverse_iterator",
    "size_type", "value_type",
    "ios_base","ios", "istream", "ostream", "iostream",
    "istringstream", "ifstream", "ostringstream", "ofstream",
    "stringstream", "fstream",
    "streambuf", "stringbuf", "filebuf",
    "streamoff", "streampos",
    "input_iterator_tag", "output_iterator_tag", "forward_iterator_tag",
    "bidirectional_iterator_tag", "random_access_iterator_tag",
    "pair", "auto_ptr", "allocator", "raw_storage_iterator",
    "vector", "list", "deque",
    "set", "multiset", "map", "multimap",
    "stack", "queue", "priority_queue", "bitset",
    "shared_ptr", "weak_ptr", "unique_ptr"
  };
  for (i=0; i<array_size(cpp_types); ++i)
    id_lookup(cpp_types[i],NULL,type_defined);
  id_lookup("asm", NULL, int_like);
  id_lookup("and", NULL, and_like);
  id_lookup("bool", NULL, int_like);
  id_lookup("catch", NULL, catch_like);
  id_lookup("class", NULL, struct_like);
  id_lookup("delete", NULL, sizeof_like);
  id_lookup("explicit", NULL, int_like);
  id_lookup("false", NULL, expression);
  id_lookup("friend", NULL, int_like);
  id_lookup("inline", NULL, int_like);
  id_lookup("namespace", NULL, namespace_like);
  id_lookup("new", NULL, sizeof_like);
  id_lookup("not", NULL, not_like);
  id_lookup("operator", NULL, case_like);
  id_lookup("or", NULL, and_like);
  id_lookup("private", NULL, case_like);
  id_lookup("protected", NULL, case_like);
  id_lookup("public", NULL, case_like);
  id_lookup("template", NULL, template_like);
  id_lookup("this", NULL, expression);
  id_lookup("throw", NULL, throw_like);
  id_lookup("true", NULL, expression);
  id_lookup("try", NULL, try_like);
  id_lookup("typeid", NULL, sizeof_like);
  id_lookup("typename", NULL, typename_like);
  id_lookup("using", NULL, int_like);
  id_lookup("virtual", NULL, int_like);
  id_lookup("xor", NULL, and_like);
  id_lookup("const_cast", NULL, sizeof_like);
  id_lookup("static_cast", NULL, sizeof_like);
  id_lookup("dynamic_cast", NULL, sizeof_like);
  id_lookup("reinterpret_cast", NULL, sizeof_like);
}

@ There are a few more kinds of elementary scraps that the functions we have
given before can produce, which we mention here for completeness. Ordinary
identifiers get category |expression|, and their names will be expanded on
output as argument to a control sequence that provides the proper formatting.
For strings, constants, verbatim constructions, and \TeX~strings, the
applicable control sequences and the constituent characters (escaped with
backslashes where necessary) are written explicitly into token memory;
their scraps also have category |expression|. Comments are converted to
scraps of category |insert|, and their contents are also stored literally;
in case of `\pb' fragments a reference to a text marked with
|inner_text_flag| is stored, for the production of which the parsing
mechanism has in fact already been invoked. By contrast module names are
stored by reference to the name table, just like identifiers, and their
scraps have category |mod_scrap|; in their case the parsing mechanism may
be called during the {\sl output\/} process, if any `\pb' constructions
occur.

The final \TeX\ output produced for elementary scraps will often be marked
with special control sequences.  Ordinary multi-character identifiers are
enclosed in `\.{\\\\\{}$\,\ldots\,$\.\}' (single character identifiers are
merely preceded by a space; they will be set in math italic), identifiers
whose |ilk| is |TeX_like| will become control sequences that are
also enclosed in `\.{\\\\\{}$\,\ldots\,$\.\}' (to establish italic type),
reserved words are enclosed in `\.{\\\&\{}$\,\ldots\,$\.\}', strings and
all-caps identifiers in `\.{\\.\{}$\,\ldots\,$\.\}', constants in
`\.{\\T\{}$\,\ldots\,$\.\}', and verbatim constructions in
`\.{\\vb\{}$\,\ldots\,$\.\}'. Comments are enclosed in
`\.{\\C\{}$\,\ldots\,$\.\}' and usually followed by `\.{\\6}' (a forced
break), and module names take the form `\.{\\X$n$:}$\,\ldots\,$\.{\\X}'
where |n| is the section number (since module names have
|mathness==yes_math|, there is no danger that the final `\.{\\X}' will
disable a following space when coming from a `\hbox{\.{\v@@< ... @@>\v}}'
construction).

@*1 Internal representation of the grammar. @:title@>
We will now consider the how the reduction rules themselves are
represented and used. As we have seen, a rule must define a sequence of
categories for its left hand side, and for its right hand side a category
and a prescription for constructing its translation. In addition, some
categories of the left hand side may be marked as context, so that they will
not take part in the reduction, and there is a way to specify conditional
loading of rules. A few more pieces of information are included for
convenience and efficiency.

Individual reduction rules are stored in a structure called |reduction|.
It is organised in a way that allows for semi-static initialisation, i.e.,
the essential parts of information are stored near the beginning of the
structure or of one of its sub-structures, so that they can be defined by
an initialiser expression, while some further fields are computed from them
and assigned at startup time. Within the fields that are statically
initialised some fields that usually are~0 are put at the end, so that in
the default case they can be omitted from the initialiser.

The field |id| holds an identification number for the rule, which is used in
debugging. Then follows the left hand side information, consisting of an
array of at most |max_lhs_length| categories (which include those of a
possible context; if less than the maximal number of categories are present
they are padded with zeros), followed by integers |context| and |length|.
The field |context| specifies which categories, if any, form the context:
this can be a sequence of one or more categories at either end of the left
hand side of the rule, but not at both ends. If |context==0| (as is the
case if no explicit initialiser is specified), there are no context
categories; when |context>0|, the first |context| categories form the left
context, and when |context<0|, the last |abs(context)| categories form the
right context. (In practice it is wise to use only token categories (ones
that do not require reduction to be formed) for a right context, unless one
can be quite sure that no unintended reduction will affect the categories
taking part in the intended reduction before the right context has been
reduced.) There must be at least one category that is not part of the
context, lest a ``reduction'' would increase the number of scraps. The right
hand side of a rule specifies a category and a string used as a format to
build up the translation; for the common case that the translation is formed
by concatenating the translations of all scraps of the left hand side (not
including those those of the context), a null pointer may be given instead
of a format string. The following field |mask| can be used to specify
selective loading of rule at startup time: any bit set in it will suppress
loading under some condition dependent on the setting of option flags in the
call of \.{\me.}. The field |displacement| is computed at startup time to
record the number of positions (usually negative) by which the parsing
pointer~|pp| should be changed after application of the rule to have any
chance of subsequently matching another (or the same) rule.

@d max_lhs_length 4

@< Typedef and enumeration declarations @>=

typedef struct
{ short id; /* for debugging */
  struct
    {@; eight_bits category[max_lhs_length]; signed char context,length; } lhs;
  struct {@; eight_bits category; char *translation; } rhs;
  sixteen_bits mask;
  short displacement;
} reduction;

@ We shall organise the rules in a ``trie'' structure for fast matching. Such
a structure is a set of nodes, where each node represents a sequence of
categories that occurs as an initial part of at least one rule; each node
contains links to other nodes that may be reached by adding one more category
to the sequence, and possibly to a rule if its sequence of categories is a
complete left hand side. The very first node is the root, and represents the
empty sequence of categories; each node can be reached by an unique sequence
of links from this root, so the graph structure defined by the links is that
of a directed tree. This kind of structure is suitable for rapid searching,
since at each point where some categories match, the question whether any of
the rules can extend that match with the next category found can be answered
by testing the presence of a single link. If the match fails one could imagine
having some pointer to the node that would be reached by the largest final
sub-sequence of the one already matched that could be the starting point for
the next attempt at a match, which would speed up matching even more. However
we do not include such pointers in our structure, so after a failed match we
must start again at the root of our trie after advancing the place where we
will try to match the first category.

Here are the details of our trie structure. If |q| points to a trie node
reached after matching some sequence of categories, and that sequence
corresponds to the left hand side (including context) of some rule (which
should be unique), then |q->rule| points to that rule, otherwise
|q->rule==NULL|. If that sequence of categories is a proper prefix of the left
hand side of a rule (which may happen whether or not |q->rule==NULL|), and $c$
is the next category in that left hand side (which implies |valid_cat(c)|, and
in particular $c\neq0$), then |q->next[c-1]| is the index of the trie node
reached after a further match of~$c$. The entries of |q->next| that do not
correspond to any such successor node are set to~0, which is unambiguous
because the root of the trie does not figure as a successor of any node. We do
not attempt a sparse representation (which would avoid storage of such 0's),
but we do use a relatively compact |sixteen_bits| representation for the
entries of |q->next|, which saves a considerable amount of space, since there
is a total of |(min_nonrule_cat-1)*max_no_of_nodes| such entries? (The entries
used to be |eight_bits| number, but for a decent grammar for \Cpp\ a limit of
256 trie nodes was too restrictive. One could extend the possibilities of such
small numbers a but by using the fact that the successor of a trie node always
has a strictly greater index that the node itself, so one could store the
difference and change the logic for computing successors below, hoping that
even with more than 256 nodes, one would not need strides larger than~255.
However anno~2006 when memory is measured in gigabytes, so much parsimony of
memory seems ridiculous.)

@<Typedef...@>=

typedef struct {@; reduction* rule; sixteen_bits next[min_nonrule_cat-1]; }
  trie_node;

@ Trie nodes are allocated from an array |trie_node_mem|, with the root of
the trie at |trie_node_mem[0]|. The root of the trie is special in that is
allocated initially (namely |node_no| is initialised to~|1|), but it is not
initialised. It corresponds to the state in which no scraps have been matched
yet, and only its |next| field is used; its entries will be set as a
by-product of installing reduction rules. In general a rule will be determined
by the |rule| field of its final node, and appropriate entries of the |next|
fields of certain nodes, which determine the unique path to that final
node.

We introduce some macros that will help us find our way around the trie. The
address the successor of the trie node pointed to by |q| for category |c| can
be written as |successor(q,c)|. The absence of such a successor can be found
by testing |no_successor(q,c)|. If |x| is the address of any node in the tree
(except the root), then we can make that node the successor of |q| for
category |c| by invoking |set_successor(q,c,x)|.

@d trie_root (&trie_node_mem[0])
@d successor(q,c)
   (&trie_node_mem[(q)->next[(c)-1]]) /* category codes start at 1 */
@d no_successor(q,c) ((q)->next[(c)-1]==0)
@d set_successor(q,c,x) ((q)->next[(c)-1]=(sixteen_bits)((x)-trie_node_mem))

@<Global...@>=
trie_node trie_node_mem[max_no_of_nodes];
int node_no = 1;  /* number of trie nodes allocated */
#ifdef DEBUG
boolean install_failed=false;
#endif

@ Trie nodes are allocated in a straightforward sequential way. We don't
trust that uninitialised statically allocated pointers will be |NULL|
(although they should), especially not on machines where |NULL| is not
represented as ``all bits cleared'', so we do a bit of extra work here.  For
the entries of the array |next| we expect no problems however (since
|eight_bits| is an integral type), so we do not explicitly initialise them.

@c trie_node *get_new_trie_node(void)
{ if (node_no>=max_no_of_nodes)
    overflow("trie node"); @.trie node capacity exceeded@>
  trie_node_mem[node_no].rule=NULL;
  return &trie_node_mem[node_no++];
}

@ The function |install_rule| installs a reduction, and if \.{\me.} was
compiled with |DEBUG| set, it also performs some checks on the validity of
the rule; if any check fails the variable |install_failed| is set to |true|.
Since |print| has a variable number of arguments, the macro |rule_error|
does not incorporate them but just prepends |print| to the argument list;
the replacement text of |rule_error| can therefore not be parenthesised, and
the macro should be used with some care. Note however that a call to
|rule_error| (with argument list) can safely by used as the then-branch of an
|if|-statement.

@d rule_error install_failed=true,print
   /* report a problematic syntax rule */

@ A global variable |rule_mask| is set at startup time according to the
relevant option flags; any rules~|r| for which |rule_mask & r->mask!=0| are
suppressed. For each set of mutually exclusive settings, a number of bits in
|rule_mask| is reserved equal to the size of the set; the current setting
will have its corresponding bit set while the others are cleared. Therefore
the bits set in |r->mask| specify the option settings for which the rule is
disabled, and the default state |r->rule==0| means that the rule always
applies, regardless of any optional settings.

@< Global... @>= sixteen_bits rule_mask;

@ The function |install_rule| enters a rule into the trie structure, and if
|DEBUG| is defined, performs some sanity checks on the rule.

@< Prototypes @>=
void install_rule (reduction* rule);

@~The length of rules can be found because category~0 is not used for scraps.
@c
void install_rule(reduction *rule)
{ if ((rule_mask & rule->mask)==0)
  { eight_bits* p=rule->lhs.category, i=0;
    while (i<max_lhs_length && p[i]!=0) ++i;
    rule->lhs.length=i;
#ifdef DEBUG
    @< Check left-hand side of |rule| @>
    @< Check right-hand side of |rule| @>
#endif
    @< Install |rule| in the trie structure @>
    @< Compute |rule->displacement|, and if context is specified modify
       |rule->lhs.length| and |rule->context| @>
  }
}

@ The left-hand side should not be only context, and all categories should
be legal ones.

@< Check left... @>=
{ if (rule->lhs.length<=abs(rule->lhs.context))
    rule_error("\nNo scraps to replace in rule %d.\n", rule->id);
		@.No scraps to replace...@>
  for(i=0; i<rule->lhs.length; ++i)
    if (!valid_cat(p[i]))
      rule_error("\nUnknown category %d in LHS of rule %d.\n", p[i], rule->id);
		  @.Unknown category...@>
}

@ The right-hand side should have a valid category, and unless its
translation is |NULL| (or |""| which we treat as if it were |NULL|), it
should contain as many times |'_'| as there are non-context categories in
the left-hand side.

@< Check right... @>=
{ int c=rule->rhs.category; char* s=rule->rhs.translation;
  if (!valid_cat(c))
      rule_error("\nUnknown category %d in RHS of rule %d.\n", c, rule->id);
		  @.Unknown category...@>
  if (s!=NULL)
  { if (*s=='\0') s=rule->rhs.translation=NULL; /* replace empty string */
    else
    { i=0;
      do
	if (*s!='p') i+= *s++=='_'; /* count underscores */
	else if (++s,isdigit((eight_bits)*s)) ++s; /* skip digit and advance */
	else rule_error("\nDigit should follow 'p' in format of rule %d.\n"
			 @.Digit should follow 'p'...@>	, rule->id);
      while (*s!='\0');
      if (i!=rule->lhs.length-abs(rule->lhs.context))
	rule_error("\nCount of '_' not equal to length LHS in rule %d.\n"
		    @.Count of '\_' ...@> , rule->id);
    }
  }
}

@ Since trie nodes are not represented sparsely, insertion is easy. Just
before moving to a successor node we make sure it exists, by adding one to the
trie if necessary. It does not matter if the final nodes already existed, but
there should not already be a rule attached to it, which would mean that a
previous rule has an identical left hand side (here context is taken into
account to determine equality).

@<Install |rule|...@>=
{ trie_node* q=trie_root;
  for (i=0; i<rule->lhs.length; ++i)
  { if (no_successor(q,p[i])) set_successor(q,p[i],get_new_trie_node());
    q=successor(q,p[i]);
  }
#ifdef DEBUG
  if (q->rule!=NULL)
    rule_error("\nIdentical left-hand sides in rules %d and %d.\n"
	        @.Identical left-hand sides...@>   , q->rule->id, rule->id);
#endif
  q->rule=rule;
}

@ We compute |displacement| conservatively, based on local considerations;
alternatively we might also consider the whole set of rules to find larger
(less negative) values that would make parsing go a bit faster. A rule can
have a left hand side of length |max_lhs_length|. This means that it is safe
to move |pp| so that it will afterwards be |max_lhs_length-1| positions to the
left of the first scrap that may have obtained a new category. That scrap is
the one corresponding to first non-context category in the left hand side of
the rule, unless that category is the same as the category of the right hand
side, in which case it is the next scrap (this is for reduction rules whose
effect on the category sequence is the ``absorption'' of some categories
further to the right; then the mentioned next scrap did not participate in the
reduction and we cannot predict its category (unless there was right context,
but we do not attempt to push that exceptional case any further)).

After a rule has been installed, there is no need to record the full length
of the left hand side, including context, since this is implicit from the
place in the trie where the pointer to this rule is located; rather we store
the number of scraps that will be replaced. Similarly it is more useful to
know the offset of the first scrap to be replaced (which is~0 in case of a
right context) rather than the value of |context| as stored at
initialisation. This means that if no context is specified nothing changes,
and if there are $k$ categories of context then $k$ is subtracted from
|rule->lhs.length| and the offset will be~$k$ in case of a left context,
and~$0$ otherwise. The code below calculates the offset in~|k|, and the
non-positive displacement is obtained by subtracting |max_lhs_length-1| from
|k|, or from |k+1| in the absorption case mentioned above. We check that an
absorption rule actually absorbs at least one scrap; a rule that violates this
requirement would lead to a fruitless reduction, and what is worse to a value
|displacement>0| that would invalidate the logic of the reduction code below.

@< Compute |rule->displacement|... @>=
{ int k=rule->lhs.context,d;
  if (k<0) {@; rule->lhs.length+=k; k=0; } @+else rule->lhs.length-=k;
  d=k-(max_lhs_length-1); /* this cannot be positive */
  if (rule->lhs.category[k]==rule->rhs.category) /* no category change */
  { ++d;
#ifdef DEBUG
    if (rule->lhs.length==1)
      rule_error("\nNo categories change in rule %d.\n", rule->id);
		  @.No categories change...@>
#endif
  }
  rule->lhs.context=k;
  rule->displacement=d; /* if positive, an error was reported */
}

@ Once the trie has been installed, the function |match| tests whether the
category pointed to by~|p|, and its successors, match the categories in the
trie structure, starting at the root of the trie, and up to a node that
contains a rule. Only valid categories are taken into account, so in
particular the category~0 used to mark the end of the scrap sequence will stop
the matching loop. Multiple matches are possible, in which case the longest
one takes precedence; this means we use the last non-null |rule| field
encountered. This need not be the one in the last trie node seen, since we may
have successfully matched a category that only occurs in rules that require
some further scrap of a category that is not present (it seems unlikely that
this actually occurs in the given set of rules, but we do not count on it);
therefore we maintain a variable |rule| to store last non-null |rule| field
seen. In any case it may happen that no rule at all has completely matched,
possibly even that no single scrap has matched; in such cases |match|
returns~|NULL|. Note the we avoid using a side effect in the argument of
|valid_cat|, as we must.

@c
reduction *match (scrap_pointer p)
{ trie_node* q=trie_root; reduction* rule=NULL; int c;
  while (c=p++->cat,valid_cat(c) && !no_successor(q,c))
    if ((q=successor(q,c))->rule!=NULL) rule=q->rule;
  return rule;
}

@*1 The parsing mechanism. @:title@>
Conceptually the working sequence of scraps is like a deck of cards, in
which we repeatedly replace a sequence of consecutive cards by a single new
card. Since such replacements never increase the number of cards, we can use
sequential allocation for the current sequence of scraps, and our only
difficulty will be how to conveniently fill the holes that might be left after
each reduction step. Now reduction usually takes place near the beginning of
the scrap sequence (assuming that the scrap sequence makes syntactical sense)
because that is where we are looking first, and we want to avoid shifting
down the whole remainder of the scrap sequence each time. Therefore the
sequence of scraps, which initially occupies the positions from |scrap_base|
to |scrap_ptr|, is allowed to have a hole in its middle, the low and high
end of which are pointed to by variables |lo_ptr| and~|hi_ptr|. There is
also a variable that points to the place where reductions are currently
taking place, which is the parsing pointer~|pp|. It will always point into
the area below the hole, and when it approaches the hole so closely that a
potential reduction might involve scraps from above, the situation is
remedied by sliding down scraps to the lower region, effectively raising the
hole. Therefore the scraps in the higher region are those that have never
been considered for a reduction yet. Eventually all scraps have been moved
down (i.e., we have |hi_ptr==scrap_ptr|), and after that has happened a
scrap with category~0 (which is not otherwise used) is copied down to signal
the imminent end of the reduction process. When finally no more rules
match the scraps in the lower region, the parsing stops.

@< Global variables @>=
scrap_pointer pp; /* current position for reducing scraps */
scrap_pointer lo_ptr; /* end of sequence of scraps that have been examined */
scrap_pointer hi_ptr; /* first scrap that has not been examined */

@ The |mathness| is an attribute of scraps that says whether their
translation is to be processed in a math mode context or not. Since the
translation can be concatenated from a large number of other scraps, there
can be switches in and out of math inside the translation, and we need to
specify the mathness at each of the boundaries. For some scraps it either
makes no difference whether their translation is processed in math mode or
not, or the required mathness is to be determined by the grammatical context
rather than by the scrap itself. Such scraps have mathness |maybe_math| at
both ends; otherwise a definite mathness is specified at either end. The
least significant pair of bits of the |mathness| field of a scrap controls
the right boundary, and the pair of bits to its left controls the left
boundary.

@d left_math(a) (a->mathness>>2)
@d right_math(a) (a->mathness&0x3)

@ If we combine two scraps neither of which has mathness |maybe_math| at its
boundaries, then a `\.\$' is inserted in between if and only if the
mathnesses at the common boundary do not agree; if a scrap with |maybe_math|
joins one with a definite mathness, that mathness is propagated across the
former scrap to its other boundary. In order to implement this, we maintain
two mathness values while building up a text: |init_mathness| and
|cur_mathness| which represent the values at the left and right boundaries
of the part contributed so far; these are local variables of whatever
function is concatenating translations, but they should be called by these
names since they are addressed by the macros below. As a consequence of the
left-to-right order of combining translations, a |maybe_math| scrap that is
combined with scraps with definite mathnesses, will actually be set in the
mode inherited from its left (unless it appears as the leftmost scrap in a
reduction); this can be used to make certain symbols, such as colons, behave
in two slightly different ways depending on their syntactic function. (This
method is not infallible however, as a comment following the symbol will
always force it to be processed in horizontal mode; this happens because
|insert| scraps are tacked onto the scrap before them before any ordinary
reduction can affect it.)

Before scraps requiring some definite mathness are contributed, we invoke
either |set_mode(yes_math)| or |set_mode(no_math)| as appropriate; the first
time this happens will determine the value of |init_mathness|.

@d set_mode(x)
if (cur_mathness==maybe_math) cur_mathness=init_mathness=x;
else if (cur_mathness!=x) {@; app('$'); cur_mathness=x; }
@q $ emacs-cookie @>
else @; /* expect semicolon */

@ The macro |app_trans| is invoked with a |scrap_pointer| as argument,
and appends its translation as a single token; |add_trans| will in addition
to this administrate |init_mathness| and |cur_mathness|, and interpolate
any necessary math shifts.

@d app_trans(a) app_tok(text_flag+text_index((a)->trans))
@d add_trans(a)
{ scrap_pointer scr=a; /* so that evaluating |a| may have side effects */
  if (left_math(scr)!=maybe_math)
  { if (cur_mathness==maybe_math) init_mathness=left_math(scr);
    else if (cur_mathness!=left_math(scr)) app('$');
@q $ emacs-cookie @>
    cur_mathness=right_math(scr);
  }
  app_trans(scr);
}

@ The function call |fuse(s,n)| will concatenate the translations of |n|
scraps starting with |*s|, taking care of the mathnesses, and install the
resulting text into |s->trans|.

@c
void fuse (scrap_pointer s, int n)
{ int cur_mathness=maybe_math, init_mathness=maybe_math; scrap_pointer p=s;
  check_toks(n); check_text();
  do add_trans(p++)@; while (--n>0); /* gather all the translations */
  s->trans=text_ptr; freeze_text();
  s->mathness=(init_mathness<<2)+cur_mathness;
}

@ When a matching rule has been found, the function |reduce| is called to
perform the corresponding actions. At that point |pp| points to the first
scrap involved in the match, and the argument |rule| to |reduce| points to
the matching rule.

If a rule has a left hand side of length~1 (not counting context) and also
the default translation (plain concatenation), then all that is to be done
is to change the category of a scrap, and part of the processing can be
skipped.

@c
@< Utility functions needed during reduction @>@;
void reduce (reduction* rule)
{ int k=rule->lhs.context, l=rule->lhs.length;
  scrap_pointer s = pp+k, p=s;/* position of the new scrap */
  char* f=rule->rhs.translation; /* format string for translation */

  s->cat=rule->rhs.category;

  if (l>1 || f!=NULL) /* otherwise ready already */
  { if (f==NULL) fuse(s,l), p+=l; /* default translation */
    else @<Generate token list according to format string |f| @>
    if (l>1) @<Fill vacant scrap positions@>
  }

  @<Print a snapshot of the scrap list if debugging@>
  @<Change |pp| to |max(pp+rule->displacement,scrap_base)| @>
}

@ When we have applied a reduction to the sequence of scraps, we usually
remove scraps (we never create more scraps than we remove), thereby creating
a small ``hole'' in the sequence. We fix that hole by sliding scraps
downward across it, thereby moving the hole upwards, until it reaches the
``official'' hole at |lo_ptr|; then |lo_ptr| is adjusted so that the small
hole is incorporated in the official hole. During the translation process
the pointer~|p| was moved across all scraps that took part in the reduction,
so the scraps to move are at positions~|i| with |p<=i<lo_ptr|.

@<Fill vacant...@>=
{ scrap_pointer q=s+1; /* position after the newly formed scrap */
  while (p<lo_ptr) *q++=*p++;
  lo_ptr=q;
}

@ Here we are careful not to create any pointers below |scrap_base|,
since they are not guaranteed to exist (but mentioning them in the module
name will do no harm). The code below uses the fact that
|rule->displacement<=0|.

@<Change |pp| to...@>=
if (pp<scrap_base-rule->displacement) pp=scrap_base; @+
else pp+=rule->displacement;

@ We need no extensive coding mechanism for describing translations, since
they all follow a similar pattern. In all cases all the translations of the
scraps in the left hand side (not including the context scraps) are used
exactly once, in left to right order (violation of these principles would
result in very strange effects indeed for the printed output). The only
things that need to be added are formatting controls like |indent|, |force|
or |break_space|, math shifts (but these are already taken care of by
|add_trans|), and white space. We may also specify calls of
|make_underlined| or |make_unreserved| for certain scraps in the left hand
side. Although the items inserted are of a modest variety, one should
realise that their presence is the only reason we need to parse at all;
without them the translations could have been computed by purely lexical
means.

In format strings an underscore indicates the translation of the next scrap
of the left hand side. Other characters each encode a formatting control; the
character `\.p' encodes |opt| and is followed by a digit that becomes its
argument. The characters `\.!' and `\.\$' respectively cause
|make_underlined| and |make_nonreserved| to be called for the next scrap;
`\.\ ' and~`\.\~' produce a space in the translation, where the latter is
non-breakable and the former forces horizontal mode. To just force
horizontal or math mode there are `\.h'~and~`\.m'; the latter avoids the
possibility of a completely empty formula by adding a space inside math
mode. The precise meaning of these and other formatting characters is easily
read off from the code below. The `\.o' and `\.r' format characters (the
latter is used only in compatibility mode) affect the math category used by
\TeX\ for the next symbol (character or control sequence); the syntax rules
using them do not put braces around that symbol, since these could also
capture a following comment, causing a \TeX\ error.

The number of free tokens we require to be available by calling |check_toks|
is a conservative estimate, based on a hypothetic ``worst case'' reduction
rule with a left hand side of length~4 and with 8 additional items in its
translation (counting format codes with multi-character translations with
multiplicity), whose mathnesses alternate, so that a maximal number of math
shifts is required.

@< Generate token list... @>=
{ int cur_mathness=maybe_math, init_mathness=maybe_math;
  check_toks(23); check_text();
  do
    switch (*f++) @:format@>
    { case '+': app(indent); break;
      case '-': app(outdent); break;
      case 'p': app(opt); app(*f++); break; /* penalty with numeric argument */
      case 'f': set_mode(no_math); app(force); break;
      case 'F': set_mode(no_math); app(big_force); break;
      case 'b': set_mode(no_math); app(backup); break;
      case 'B': set_mode(no_math); app(break_space); break;
      case 't': set_mode(yes_math); app_str("\\a"); break; @.\\a@>
					/* next item in tab space */
      case ',': set_mode(yes_math); app_str("\\,"); break; @.\\,@>
					/* thin space */
      case 'h': set_mode(no_math); break; /* force horizontal mode */
      case 'm': set_mode(yes_math); app(' '); break;
				/* force math mode, avoid `\.{\$\$}' */
      case 'o': set_mode(yes_math); app_str("\\m"); break; @.\\m@>
					/* make ``mathord'' */
      case 'r': set_mode(yes_math); app_str("\\MRL"); break; @.\\MRL@>
					/* make ``mathrel'' */
      case 'a': set_mode(yes_math); app_str("\\ang"); break; @.\\ang@>
		         /* change `\.<' or `\.>' to `$\ang<$' or~$\ang>$'  */
      case '!': make_underlined(p); break;
      case '$': make_nonreserved(p); break; @q $ emacs-cookie @>
      case ' ': set_mode(no_math); app(' '); break;
      case '~': case '{': case '}': app(f[-1]); break;
               /* insert character literally */
      default: printf("%c: ",f[-1]);
        confusion("illegal character in format string");
      case '_':	add_trans(p++);
    }
  while (*f!='\0');
  s->trans=text_ptr; freeze_text();
  s->mathness=(init_mathness<<2)+cur_mathness;
}

@ An |int_like| identifier following a |struct_like| token, a selection
operator (`|.|' or `|->|'), or the special code \:;, is to be typeset as an
ordinary identifier. The function |make_nonreserved| alters the flag of the
token representing such an identifier occurrence. The scrap representing the
|int_like| identifier should not be formed by any reduction, but come directly
from |C_read|, so in principle we expect the translation of our scrap to be an
unnested text consisting of a single token, and this text is only used as
translation of our scrap, so we can modify its token in place. However, a
comment (or other |insert|) directly following the |int_like| identifier may
complicate this picture slightly, because |insert| scraps are tacked onto the
previous scrap before it gets the chance to take part in any reduction; this
means our token may be buried inside one or more levels of text nesting, but
still is the very first token of the translation.  Any such levels of nesting
are soaked off by the |while| loop below; after this process a single-token
text should remain, possibly containing a reserved word. If this is the case
we replace |res_flag| by |id_flag| to make it print as an ordinary identifier.

@< Utility functions needed during reduction @>=
void make_nonreserved (scrap_pointer p)
{ text_pointer q=p->trans; token t;
  while (text_flag<=(t=text_begin(q)[0]) && t<inner_text_flag)
    q=text_at(t-text_flag);
  if (text_end(q)==text_begin(q)+1 && res_flag<=t && t<mod_flag)
    text_begin(q)[0] = t-res_flag+id_flag;
}

@ We wish to mark defining occurrences of identifiers, triggered by the
matching of certain reduction rules. Thanks to the unsurpassed
declaration syntax of~\Cee, it is a non-trivial task to locate the
identifier in question, since it may already have been wrapped up in a more
or less complicated declarator or function heading. Fortunately the
identifier is always the first one occurring in that expression
\unskip\footnote{${}^\dagger$}
{Due to the possible occurrence of the keywords |const| and |volatile| in
 strange places within declarations, this is not entirely true. We will
 arrange things however in such a way that unless we are forced otherwise by
 enclosing parentheses, these keywords will remain outside the text that is
 passed to |first_ident|. If we ignore the possibility of specifying
 redundant parentheses, this only causes problems when one wishes to
 declare certain peculiar things like constant function pointers.
}, so a left-to-right scan through the, possibly nested, text forming the
translation of the expression should reveal it. We accept either an ordinary
identifier or a reserved word; the latter possibility is needed for typedef
declarations. @^recursion@>

@< Utility functions needed during reduction @>=
id_pointer first_ident(text_pointer p)
{ token_pointer q; token t;
  if (p>=text_ptr) confusion("first_ident"); @.first\_ident@>
  for (q=text_begin(p); q<text_end(p); ++q)
    if (id_flag<=(t=*q) && t<mod_flag) return id_at(t%id_flag);
    else if (text_flag<=t) /* text or inner text */
    {@; id_pointer r=first_ident(text_at(t%id_flag));
      if (r!=NULL) return r;
    }
  return NULL;
}

@ In the following situations we want to mark the occurrence of an identifier
as a definition: in a declaration of a variable, as |r| in `|xref_pointer
*r@;|', when a label is declared, as |found| in `|found:|', in the declaration
of a function, as |f| in `|char* f(int n){@t\dots@>@;}|', after a
|struct_like| token (declaring a structure-, union-, or enumeration tag) and
for all identifiers in an enumeration list. This is accomplished by the
invocation of |make_underlined| at appropriate times, which modifies the cross
reference list for this identifier.

@< Utility functions needed during reduction @>=
void make_underlined (scrap_pointer p)
  /* underline entry for first identifier in |p->trans| */
{ id_pointer name=first_ident(p->trans); /* name of first identifier */
  if (name==NULL) return;
    /* this happens for unsyntactical things like `|int 3;|' */
  { sixteen_bits head=xref_index(name->xref),* r=&head; int n;
    while ((n=xnum(*r)&num_mask)!=0 && n<section_count) r=&xlink(*r);
    if (n==section_count) xnum(*r)|=def_flag;
    else /* this may happen for one-letter identifiers */
    { make_xref(section_count+def_flag,*r);
      if (r==&head) name->xref=xref_ptr; @+ else *r=xref_index(xref_ptr);
    }
  }
}

@ We have now seen how a match is made, and what is done once a matching
rule is found; we still have to consider how everything is set up properly,
and how rules are repeatedly applied until no more reduction is possible,
responding properly to successful and failing matches. All this is
performed by the function |translate|; as we have seen it is called by
|do_C| and |finish_C| after scraps have been stored from |scrap_base|
to |scrap_ptr| of the |scrap_info| array, and it returns a pointer to the
text representing the result of parsing all those scraps.

We start with appending a dummy scrap if either no scraps are present at
all, or some tokens remain that have not been packed into a scrap yet; the
latter can only be due to `\.{@@t...@@>}' as a final item in compatibility
mode. Then we set |lo_ptr| and |hi_ptr| appropriately, and begin to apply
rules as long as possible. When this is done and more than one scrap
remains, their translations are wrapped rather bluntly together to a single
text.
@:translate definition@>

@c
text_pointer translate (void) /* converts a sequence of scraps */
{ pp=lo_ptr=hi_ptr=scrap_base;
  if (scrap_ptr==pp || dangling_tokens()) /* then append dummy scrap */
  {@; check_scrap(); pack_scrap(insert,no_math); }
  @< If tracing, print an indication of where we are @>
  @< Reduce the scraps using the rules until no more rules apply @>
  @< Combine the irreducible scraps that remain,
     and |return| the resulting text @>
}

@ Before applying |match|, we must make sure it has good input (at least
|max_lhs_length| scraps). If a match at |pp| exists, |reduce| will perform
the required processing and updating of |pp| (in this case |pp| is
certainly not increased), if not, we move to the right and try again.


@< Reduce the scraps... @>=
do
{ reduction *rule;
  @< Make sure the entries |pp| through |pp+max_lhs_length-1|
     of |scrap_info| are defined, or that |lo_ptr->cat==0| @>
  if ((rule=match(pp))!=NULL) reduce(rule);
  else
  {@; ++pp;
    @< Take special action if |pp->cat| is |end_expr| or |rproc| @>
  }
}
while (pp<lo_ptr);

@ If we get to the end of the scrap list, we set |lo_ptr->cat=0|, which will
prevent any matches to scraps beyond those that are actually present.  All
scraps of category |insert| pass across the hole between |lo_ptr| and
|hi_ptr|, and we take the opportunity to remove them, tacking their
translation to the scrap below the hole. We even make sure that
|hi_ptr->cat!=insert|, so that the scraps with which the |insert| scrap is
combined will not have undergone any ordinary reduction yet. The only
possible remaining |insert| scrap is one at the very start of the list; it
will be handled at the end of all reductions.

@< Make sure the entries... @>=
{ scrap_pointer lo_min = pp+max_lhs_length;
  while (lo_ptr<lo_min && lo_ptr->cat!=0)
    if (hi_ptr>=scrap_ptr) lo_ptr->cat=0;
    else
    { *lo_ptr++ = *hi_ptr++;
      while (hi_ptr<scrap_ptr && hi_ptr->cat==insert)
	{@; *lo_ptr = *hi_ptr++; fuse(lo_ptr-1,2); }
    }
}

@ The category pairs |lproc|--|rproc| and |begin_expr|--|end_expr| are
special in that the don't occur in any rules, but rather serve only as
markers. When the material in between will not reduce any further the whole
construction will be wrapped up, their translations concatenated, and the
result treated as an |insert| respectively an |expression|. We do not
actually form an |insert| scrap if it is not at the start of the scrap
sequence, but rather combine everything directly with the scrap immediately
before the |lproc| scrap, leaving the category of that scrap as it is; this
is necessary, because the code that combines ordinary |insert| scraps with
their predecessor only looks at scraps when they cross the hole from |hi_ptr|
to~|lo_ptr|. We use here that |rproc-lproc==end_expr-begin_expr==1|. The
variable names |s| and |p| are chosen with the same meaning as in |reduce|,
so that we could reuse a module.

@< Take special action... @>=
if (pp->cat==end_expr || pp->cat==rproc)
{ int start=pp->cat-1; /* the opening category matching |pp->cat| */
  scrap_pointer s=pp, p=pp+1;
  while ((--s)->cat!=start && s>scrap_base) {}
  if (s->cat==start) /* if opening symbol is missing, take no action */
  { if (start==begin_expr) s->cat=expression;
    else if (s==scrap_base) s->cat=insert; @+
    else --s; /* position of new scrap */
    fuse(s,(int)(p-s));
    @< Fill vacant scrap positions @> /* using values of |p| and |s| */
    pp= s-scrap_base<max_lhs_length ? scrap_base : s+1-max_lhs_length;
  }
}

@ If the initial sequence of scraps does not reduce to a single scrap,
we concatenate the translations of all remaining scraps, separated by
blank spaces, with dollar signs surrounding the translations of scraps
where appropriate. Because of this last action we must build a new text
even if only one scrap remains.

@< Combine the irreducible... @>=
{ scrap_pointer j;
  if (scrap_base->cat==insert && lo_ptr>scrap_base+1)
  { fuse(scrap_base,2); /* merge initial |insert| into the next scrap */
    j=scrap_base; j->cat=j[1].cat; --lo_ptr; while (++j<lo_ptr) *j=j[1];
  }
  @<If semi-tracing, show the irreducible scraps@>
  check_toks(1);
  for (j=scrap_base; j<lo_ptr; j++)
  { if (j!=scrap_base) app_char_tok(' ');
    if (left_math(j)==yes_math) app('$');
    app_trans(j);
    if (right_math(j)==yes_math) app('$');
  }
  freeze_text();
  return text_ptr-1;
}

@ If \.{\me.} is compiled with |DEBUG| defined, it can be put into
debugging mode by setting |tracing| to a positive value by means of \:1,
\:2, or \:3. When |tracing| is set to level~1 any sequence of two or more
irreducible scraps remaining at the end of a call to |translate| will be
printed out. When |tracing| is set to level~2, the rule numbers and current
scrap categories will be printed out after each reduction that takes place,
and if it is set to level~3, then mathnesses at the boundaries of scraps
will also be indicated.

@<Global...@>=
#ifdef DEBUG
int tracing; /* how much parsing details to show */
#endif

@ When parsing fails to reduce everything to a single scrap, pleasing results
will probably not be obtained; it is therefore advisable to run
\.{\me.} with |tracing==trace1| before a final version of a |CWEB|
program is fixed. In order to allow this without changing the source file
itself, we initialise |tracing| to |trace1| if the flag `\.{+d}' is supplied
to \.{\me.}.

@< Set initial values @>=
#ifdef DEBUG
tracing = flags['d'] ? trace1 : trace0;
#endif

@ The following code is activated by the `\.{+d}' flag or the \:1 control
code.

@<If semi-tracing, show the irreducible scraps@>=
#ifdef DEBUG
{ if (tracing==trace1 && lo_ptr>=scrap_base+2)
  { print("\nIrreducible scrap sequence at line %d in section %d:\n"
	   @.Irreducible scrap sequence...@>  ,cur_line, section_count);
    mark_harmless();
    for (j=scrap_base; j<lo_ptr-1; j++)
      print_cat(j->cat), putchar(' ');
    print_cat(j->cat); new_line(); /* |term_line_empty| is still valid */
  }
}
#endif

@ When full tracing is enabled the following message indicates which piece
of \Cee~text is being parsed (but for section bodies it will generally show
the first line of the next section, since that has already been fetched).

@<If tracing,...@>=
#ifdef DEBUG
{ if (tracing>=trace2)
  { print("\nTracing after l.%d:\n", cur_line); @.Tracing after...@>
    if (loc>buffer+50) /* shorten long lines to keep within margins */
    {@; printf("..."); term_write (loc-50,50); }
    else term_write(buffer,loc-buffer);
    new_line(); /* |term_line_empty| is still valid */
  }
}
#endif

@ After each reduction, full tracing will print a line starting with the
rule number, followed by a display of all the categories of scraps which
have been considered until now, i.e., those at positions below |lo_ptr|.
The scrap that was produced by this reduction, which is pointed to by~|s|,
has its category highlighted by enclosing it in inverted angle brackets.
If tracing is set to~3, extra-full tracing is active, and mathnesses at the
boundaries of scraps are indicated.

@d math_char(x) ((x)==yes_math ? '+' : (x)==no_math ? '-' : '?')

@<Print a snapsh...@>=
#ifdef DEBUG
{ scrap_pointer k; /* pointer into |scrap_info| */
  if (tracing>=trace2)
  { print("\n%3d:", rule->id);
    for (k=scrap_base; k<lo_ptr; k++)
    { putchar (' ');
      if (tracing==trace3) putchar(math_char(left_math(k)));
      if (k==s) putchar('>'), print_cat(k->cat), putchar('<');
      else print_cat(k->cat);
      if (tracing==trace3) putchar(math_char(right_math(k)));
    }
    print("%s\n", hi_ptr<scrap_ptr ? " ..." : ".");
  }
}
#endif
