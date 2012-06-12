@*1 Rules of the grammar. @:title@>
\def\:#1{`\.{@@#1}'}% for in case this file is processed in isolation
We first arrange the proper setting of |rule_mask|, which will control the
selection of rules actually used. Recall that any bits set in the mask of a
rule prescribe its {\it suppression\/} when the same bit is set in
|rule_mask|; therefore for instance the bit characterising \Cpp\ is called
|no_plus_plus|, so that rules specifying it will not be loaded for \Cpp. In
some cases two masks will be combined using the bitwise-or operator `|@v|',
this means (somewhat counterintuitively) that the rule will only be selected
if the conditions represented by the two masks are {\it both\/} satisfied.
The use of the bitwise-and operator `|&|' is even more exceptional: it is
only meaningful if its two operands both select one setting of the same
three-way switch; the rule will then be selected if that switch is in either
of the two indicated positions. The |merged_decls| flag is special in that
setting `\.{+m}' only enables an extra rule, but does no disable any rules;
therefore only one bit is used for this option, and raising this bit in
|rule_mask| suppresses the rule marked with |merged_decls|.

@d cwebx		0x0001 /* use normally */
@d compatibility	0x0002 /* use in compatibility mode */
@d only_plus_plus	0x0004 /* use in \Cpp\ */
@d no_plus_plus		0x0008 /* use in ordinary \Cee\ only */
@d unaligned_braces	0x0050 /* use if `\.{+u}' flag was set */
@d aligned_braces	0x0020 /* use unless `\.{+u}' flag was set */
@d wide_braces		0x0030 /* use if `\.{+w}' set */
@d standard_braces	0x0060 /* use unless `\.{+u}' or `\.{+w}' set */
@d merged_decls		0x0080 /* use if `\.{+m}' set */
@d forced_statements	0x0100 /* use if `\.{+a}' or `\.{+f}' set */
@d no_forced_statements 0x0600 /* use unless `\.{+a}' or `\.{+f}' set */
@d all_stats_forced	0x0300 /* use if `\.{+a}' set */
@d not_all_stats_forced	0x0400 /* use unless `\.{+a}' set */

@< Set initial values @>=
rule_mask= (compatibility_mode ? 0x0001 : 0x0002)
	 | (C_plus_plus ? 0x0008 : 0x0004)
	 | (flags['w'] ? 0x0040 : flags['u'] ? 0x0020 : 0x0010)
	 | (flags['m'] ? 0x0000 : 0x0080)
	 | (flags['a'] ? 0x0400 : flags['f'] ? 0x0200 : 0x0100)
	 ;
{ static reduction rule[] = { @< Rules @>@;@; };
@/int i=array_size(rule); @+ do install_rule(&rule[--i]); while (i>0);
#ifdef DEBUG
  if (install_failed) fatal("inconsistent grammar",0);
#endif
}

@*2 Expressions. @:rules@>
These rules should be obvious. Variants of rule~2 are added for \Cpp\ to
handle the cases where `\.<' or `\.>' has nothing to do with templates. Rule~5
allows typedef identifiers to be used as field selectors in structures; rules
7~and~8 attach a parameter list in a function call. In rule~14 we prefix a
potentially binary operator such as `|*|' that is used in a unary way by a
`\.{\\mathord}' command to make sure that \TeX\ will not mistake it for a
binary operator. In simple cases such as |*p| this is redundant, but if such
operators are repeated more than one level deep, as in |**p|, \TeX\ would
otherwise treat the first operator as the left operand of the second, and
insert the wrong spacing. Moreover, typical \Cee~constructions as a cast
|(void*) &x| or a declaration |char *p@;| would confuse \TeX\ even more. In
rule~13 we need not insert `\.{\\mathord}', since operators of category |unop|
are already treated as ordinary symbols by~\TeX. Rule~15 allows destructors,
which consist of a tilde followed by a class name, to be recognised in \Cpp.

@< Rules @>=
{ 1, {{expression, unop}},			{expression, NULL}},	@/
{ 2, {{expression, binop, expression}},		{expression, NULL}},	@/
{ 2, {{expression, langle, expression}},
				{expression, NULL},only_plus_plus},	@/
{ 2, {{expression, rangle, expression}},
				{expression, NULL},only_plus_plus},	@/
{ 3, {{expression, unorbinop, expression}},	{expression, NULL}},	@/
{ 4, {{expression, select, expression}},	{expression, NULL}},	@/
{ 5, {{expression, select, int_like}},		{expression, "__$_"}},	@/
@q $ emacs-cookie @>
{ 6, {{expression, comma, expression}},		{expression, "__p3_"}},	@/
{ 7, {{expression, expression}},		{expression, NULL}},	@/
{ 8, {{expression, lpar, rpar}},		{expression, "__,_"}},	@/
{ 9, {{expression, subscript}},			{expression, NULL}},	@/
{10, {{lpar, expression, rpar}},		{expression, NULL}},	@/
{11, {{lbrack, expression, rbrack}},		{subscript,  NULL}},	@/
{12, {{lbrack, rbrack}},			{subscript,  "_,_"}},	@/
{13, {{unop, expression}},			{expression, NULL}},	@/
{14, {{unorbinop, expression}},			{expression, "o__"}},	@/
{15, {{unop, int_like}},	      {int_like, NULL},only_plus_plus}, @[@]

@ Here are some less common kinds of formulae. Processing the colon belonging
to the question mark operator in math mode will give it the proper spacing,
which is different from that of a colon following a label. Rule~21 processes
casts, since the category |parameters|, which represents parenthesised lists
specifying function argument types, encompasses the case of a single
parenthesised type specification. The argument of |sizeof| may be a type
specification rather than an expression; in \Cee\ it then must be
parenthesised (rule 22); in \Cpp\ this is not necessary, (and |sizeof_like|
might be `|new|' (rule~24). In \Cpp\ we also have a functional form of the
cast (rule~25), which also doubles as an explicit call to a constructor. The
angle-bracketed cast (using a keyword like \&{static\_cast}) will employ that
same rule after applying rule~26. A drawback of rule~25 is that
it will spuriously trigger in declaration of function pointers, since the
declared identifier will be wrapped in parentheses and preceded by a type
specifier. It seems rather hard to preempt application of rule~25
automatically for declarations involving function pointers: even though the
case can be recognised by the following left parenthesis, there is no way to
force the reduction of rule~10 with priority over rule~25 without giving both
left and right context, which we cannot. Since constructors probably occur
more frequently than function pointers, we leave it as it is; for \Cpp\ users
of \.{\me.} it would be advisable to always enclose the parentheses around a
function pointer name by `\.{@@[}' and `\.{@@]}'.

@s new sizeof

@< Rules @>=
{20, {{question, expression, colon}},		{binop, "__m_"}},	@/
{21, {{parameters, expression}},		{expression, "_,_"}},	@/
{22, {{sizeof_like, parameters}},		{expression, NULL}},	@/
{23, {{sizeof_like, expression}},		{expression,"_~_"}},	@/
{24, {{sizeof_like, int_like}},	                {expression,"_~_"}
						      ,only_plus_plus},	@/
{25, {{int_like, lpar, expression, rpar}},	{expression, NULL}
						      ,only_plus_plus},	@/
{26, {{sizeof_like, templ_params}}, {int_like,NULL}   ,only_plus_plus}, @[@]

@*2 Declarations. In a declaration in \Cee, the identifier being declared is
wrapped up in a declarator, which looks like an expression of a restricted
kind: only prefix asterisk, postfix subscript and formal parameters, and
parentheses are used. In a bottom-up parser of the kind we are using, it is
natural, and hardly avoidable, that declarators are parsed as expressions.
Therefore we start recognising a declaration when we see a type specifier
followed by the first declarator; at that point we have a succession
`|int_like| |expression| |semi|' or `|int_like| |expression| |comma|' (rules
31~and~33). It is also possible that there are no declarators at all, namely
when a |struct|, |union|, or~|enum| specifier is introduced without declaring
any variables; in that case we have `|int_like| |semi|' (rule~32). Because the
type specifier might be composite, like |unsigned long int|, and there might
moreover be storage class specifiers and type modifiers (like `|const|'), we
first contract any sequence of |int_like| items to a single one (rule~30).
This rule must be applied with precedence over rule~24 above to properly parse
for instance |new|~|unsigned long int|, so we add a variant of rule~30 that
makes this happen (rule~35). In case a declarator is followed by a comma we
reduce to |int_like|, so that the next declarator can be matched, otherwise we
reduce to |declaration|.

It is not quite true that declarators always look like expressions, since the
type modifiers `|const|' and `|volatile|' may penetrate into declarators. When
they do they will almost always be preceded by an asterisk, and rule~34 will
treat such cases. The choice for |int_like| as the result category is not
completely obvious, since it makes the modifier and the preceding asterisk
part of the type specifier rather than of the declarator, which strictly
speaking is not correct; the choice for |unop| or |unorbinop| might therefore
seem a more logical one. One reason for not doing that is that a space would
have to be inserted into the translation after the modifier scrap, which would
not look right in abstract declarators for contrived cases like \hbox{|int
f(char *const)@;|}; more importantly, if the modifier would become part of the
declarator, it would be a (reserved) identifier that precedes the identifier
actually being declared, and when the declarator then receives a call from
|make_underlined| by rule 31~or~33, it would mislead |first_ident|. The
current solution has a small flaw as well, since it cannot handle the
situation where the modifier is separated from the type specifier by a
parenthesis, as in $\&{void}~(\m*\&{const}~\m*f)~(\&{int})$; apart from the
|make_underlined| problem, such rare cases are hard to handle without causing
trouble in other situations, so we do not attempt to handle them.

@< Rules @>=
{30, {{int_like, int_like}},		   {int_like, "_~_"}},	@/
{31, {{int_like, expression, semi}},	   {declaration, "_~!__"}},	@/
{32, {{int_like, semi}},		   {declaration, NULL}},	@/
{33, {{int_like, expression, comma}},	   {int_like, "_~!__p8"}},	@/
{34, {{unorbinop, int_like}},		   {int_like, "o__"}},		@/
{35, {{sizeof_like,int_like, int_like},1},
				 {int_like, "_~_"},only_plus_plus},	@[@]

@ If a typedef identifier is simultaneously used as a field selector in a
|struct| or |union| declaration, it must be made to parse as expression and
be printed in italic type; this can be achieved by placing the magic wand
\:; before the identifier, by rule~40. The reason that we place \:; at the
beginning rather than at the end of the construction here, is to prevent the
|int_like| identifier from combining with something before it first.
Rule~40 only applies if the \:; does not match by any rule with what comes
before it.

Rule~41 parses new-style (\caps{ANSI/ISO}) headings of function definitions
(while K|&|R style headings are parsed as expressions); the resulting
|function_head| will not be incorporated into a |declaration| (unless a comma
or semicolon follows) but rather into a |function|. If the parameter
specifications include identifiers (as in the case of function headings), the
arguments look like declarations without the final semicolon; rule~42 (with
aid of rule~33) constructs such parameter lists, and assures that if it is
necessary to break such a function heading across lines, the continuation
lines have three notches of (additional) hanging indentation. Parameter
specifications using abstract declarators (without identifiers) will be
treated below. In |struct| declarations we may encounter bit-field
specifications consisting of a colon followed by a number; we contract that
into the following semicolon in rule~43, after which de declaration of the
field can be handled by rule 31~or~32.

@< Rules @>=
{40, {{magic, int_like}},		   {expression, "_$_"}},	@/
@q $ emacs-cookie @>
{41, {{expression, parameters}},	   {function_head, "_B_"}},	@/
{42, {{lpar, int_like, expression, rpar}}, {parameters, "_+++_~!_---_"}},@/
{43, {{colon, expression, semi}},          {semi, "m___"}},		@[@]

@ Abstract declarators are used after type specifiers in casts and for
specifying the argument types in declarations of functions or function
pointers. They are like ordinary declarators, except that the defined
identifier has been ``abstracted''; an example is `|**(* )(int)|' in `|void
g(char**(* )(int))@;|', which tells that |g| takes as argument a pointer to
a function with |int| parameter yielding a pointer to pointer to |char|. A
difficulty with abstract declarators is that they are built up around the
vacuum left by abstracting the identifier, and since for more than one
reason we cannot allow rules with empty left hand side, we have to find an
alternative way to get them started.

The natural solution to this problem is to look for sequences that can only
occur if an identifier has been abstracted from between them, for instance
`\.{*)}' (in categories: |unorbinop| |rpar|). The most compelling reason why
in |C_read| we had to laboriously change the category of a |type_defined|
identifier to |expression| instead of |int_like| inside its defining typedef
declaration, is that it allows us to ensure that any remaining |int_like|
scrap that is followed by a |subscript| is a sure sign of an abstract
declarator.

Here are the cases that start off abstract declarators (these are the first
examples of rules that need context categories in their left hand side). As a
visual hint to the reader we leave a little bit of white space on the spot
where the identifier has vanished. Rules 50~and~51 handle declarators for
pointer arguments, where the vanished identifier is preceded by an asterisk,
which either stands at the end of the declarator, or is parenthesised (for
function pointer arguments). In these rules there is no need to prefix the
asterisk with `\.{\\mathord}', since the right context makes an interpretation
as binary operator impossible. Rule~52 is added for \Cpp, in order to allow
declarators inside angle brackets of a template argument. Rules 53~and~54
treat declarators for arrays, possibly of pointers; there are no corresponding
rules with |parameters| instead of |subscript| since abstract declarators
never specify functions themselves, only function pointers. In fact the
``function analogue'' of rule~54 would incorrectly match a cast following an
operator like `|*|' or `|-|'. Rule~55 treats an abstract declarator consisting
of subscripts only, which are redundantly parenthesised; here too the
corresponding pattern with |parameters| is not only never needed, it would
also spuriously trigger on parenthesised expressions that start with a cast.

@< Rules @>=
{50, {{unorbinop, rpar}, -1},			{declarator, "_,"}},	@/
{51, {{unorbinop, comma},-1},			{declarator, "_,"}},	@/
{52, {{unorbinop, rangle},-1},	{declarator, "_,"},only_plus_plus},	@/
{53, {{int_like, subscript},1},			{declarator, ",_"}},	@/
{54, {{unorbinop, subscript},1},		{declarator, ",_"}},	@/
{55, {{lpar, subscript},1},			{declarator, ",_"}},	@[@]

@~ Abstract declarators may grow just like ordinary declarators, to include
prefixed asterisks, as well as postfixed subscripts and parameters, and
grouping parentheses.

@< Rules @>=
{60, {{unorbinop, declarator}},  {declarator, "o__"}},		@/
{61, {{declarator, subscript}},  {declarator, NULL}},		@/
{62, {{declarator, parameters}}, {declarator, NULL}},		@/
{63, {{lpar, declarator, rpar}}, {declarator, NULL}},		@[@]

@ Here is how abstract declarators are assembled into |parameters|, keeping in
mind that the ``abstract declarator'' might be completely empty (i.e., absent)
as after `|int|' in `|void f(int);|' (rules 71~and~73). We increase hanging
indentation as in rule~42 that dealt with function headings with identifiers.
We put no space after the type specifier here, since it is followed either by
an abstract declarator, a right parenthesis or comma, so certainly not by an
identifier; therefore a space is neither necessary, nor would it improve
readability. The \caps{ANSI/ISO} syntax allows empty parentheses as a
parameter specification in abstract declarators, although this is an old-style
form; rule~74 has been included to handle this case. It also serves in \Cpp\
as the proper way to specify an empty parameter list. Fortunately a
parenthesised list of identifiers (which would parse as |expression|) is not
allowed as parameter specification. On the other hand empty parentheses are
allowed in a call to a parameterless function, which done by rule~8 which
applies with precedence over rule~74.

For \Cpp\ one has a similar list of possibilities for parsing one or more type
arguments in angle brackets for a template instantiation; here an empty
argument list or an expression argument is also permitted. We shall reduce
this to |templ_params|, rules 75--79, and the last of these gets an extra
incarnation as rule~80 to preempt the recognition of a comparison expression
using the operator~'|<|'. The result will then combine with the preceding
template class identifier by rule~81, or with the preceding ordinary
identifier by rule~82. Rule~78 handles template arguments that are ordinary
expressions.

@< Rules @>=
{70, {{lpar, int_like, declarator, comma}}, {lpar, "____p5"}},		@/
{71, {{lpar, int_like, comma}},		    {lpar, "___p5"}},		@/
{72, {{lpar, int_like, declarator, rpar}}, {parameters, "_+++__---_"}},	@/
{73, {{lpar, int_like, rpar}},		   {parameters, "_+++_---_"}},	@/
{74, {{lpar, rpar}},			    {parameters, "_,_"}},	@/

{75, {{langle, int_like, declarator, comma}},
				{langle, "____p5"},only_plus_plus},	@/
{76, {{langle, int_like, comma}},{langle, "___p5"},only_plus_plus},	@/
{77, {{langle, int_like, declarator, rangle}},
			{templ_params, "a___a_"},only_plus_plus},	@/
{77, {{langle, int_like, rangle}},
			{templ_params, "a__a_"},only_plus_plus},	@/
{78, {{langle, rangle}}, {templ_params, "a_a_"},only_plus_plus},	@/
{79, {{langle, expression, rangle}},
			{templ_params, "a__a_"},only_plus_plus},	@/
{80, {{expression, langle, expression, rangle},1},
			{templ_params, "a__a_"},only_plus_plus},	@/
{81, {{int_like,templ_params}},{int_like, NULL},only_plus_plus},	@/
{82, {{expression,templ_params}},{expression, NULL},only_plus_plus},	@[@]

@*2 Structure, union, and enumeration specifiers. It is permissible to use
typedef identifiers as structure, union, or enumeration tags as well, so we
include cases where an |int_like| follows a |struct_like| token. Rule~94
allows using |struct s| in place of a type identifier, but in \Cpp\ this is
never useful, and such a rule would preempt rule~92 in definitions of derived
classes (where the left brace is not immediately in view), so we exclude the
rule there. However, we may expect to find something like `\&{class C};' as a
pre-declaration in \Cpp, so for that case we replace rule~94 by one which will
recognise such declarations. Another difference in \Cpp\ is that one can find
|function| in a class specifier (rule~96) and also things like `\&{private}:';
the latter are parsed just like `|default:|', i.e., as a |label| (rule~98).
Rule~97 is added to allow completely empty structures or classes, which
appears to amuse mostly the users of \Cpp\ (they can now overload functions
with dummy arguments that have a distinguished type), but they are allowed in
\Cee\ as well.

@< Rules @>=
{90, {{struct_like, lbrace}}, {struct_head, "_ft_"},standard_braces},	@/
{90, {{struct_like, lbrace}}, {struct_head, "_~_"},unaligned_braces},	@/
{90, {{struct_like, lbrace}}, {struct_head, "_f_"},wide_braces},	@/
{91, {{struct_like, expression, lbrace}},
		{struct_head, "_~!_ft_"},standard_braces},		@/
{91, {{struct_like, expression, lbrace}},
		{struct_head, "_~!_~_"},unaligned_braces},		@/
{91, {{struct_like, expression, lbrace}},
		{struct_head, "_~!_f_"},wide_braces},			@/
{92, {{struct_like, int_like, lbrace}},
	{struct_head, "_~!$_ft_"},standard_braces|no_plus_plus},	@/
@q $ emacs-cookie @>
{92, {{struct_like, int_like, lbrace}},
	{struct_head, "_~!_ft_"},standard_braces|only_plus_plus},	@/
{92, {{struct_like, int_like, lbrace}},
	{struct_head, "_~!$_~_"},unaligned_braces|no_plus_plus},	@/
@q $ emacs-cookie @>
{92, {{struct_like, int_like, lbrace}},
	{struct_head, "_~!_~_"},unaligned_braces|only_plus_plus},	@/
{92, {{struct_like, int_like, lbrace}},
	{struct_head, "_~!$_f_"},wide_braces|no_plus_plus},		@/
@q $ emacs-cookie @>
{92, {{struct_like, int_like, lbrace}},
	{struct_head, "_~!_f_"},wide_braces|only_plus_plus},		@/
{93, {{struct_like, expression}}, {int_like, "_~_"}},			@/
{94, {{struct_like, int_like}}, {int_like, "_~$_"},no_plus_plus},	@/
@q $ emacs-cookie @>
{94, {{struct_like, int_like, semi}},
				{declaration, "_~!__"},only_plus_plus},	@/
{95, {{struct_head, declaration, rbrace}},
		{int_like, "_+_-f_"},standard_braces},			@/
{95, {{struct_head, declaration, rbrace}},
		{int_like, "_+f_-f_"},unaligned_braces & wide_braces},	@/
{96, {{struct_head, function, rbrace}},
		{int_like, "_+_-f_"},standard_braces|only_plus_plus},	@/
{96, {{struct_head, function, rbrace}},
   {int_like, "_+f_-f_"},(unaligned_braces&wide_braces)|only_plus_plus},@/
{97, {{struct_head, rbrace}}, {int_like,"_B_"}},			@/
{98, {{label, declaration}}, {declaration, "b_f_"},only_plus_plus},     @[@]

@ In \Cpp\ we treat |extern| in the context of ``|extern|~|"C"|~|{|'' as if it
were |struct|, otherwise |extern|~|"C"| behaves as |int|, and finally |extern|
itself can be like |static|, as in~\Cee\ (but the rule in only active for
\Cpp\ since for \Cee\ we directly produced an |int_like| category for
|extern|).

@< Rules @>=
{100, {{extern_like, expression, lbrace},-2},
				{struct_like, NULL},only_plus_plus},	@/
{101, {{extern_like, expression}},
				{int_like, "_~_"},only_plus_plus},	@/
{102, {{extern_like}},		{int_like, NULL},only_plus_plus},	@[@]

@ Rules 105--109 are for enumerations. Rules 105~and~106 preempt rule~90
above, avoiding forced line breaks. Rule 107 contracts successive enumeration
constants, and makes sure that a |break_space| and a call to |make_underlined|
follow each comma. After such contractions either a single |expression| or one
followed by a final comma remain between the braces. Rules 108~and~109 deal
with these cases, assure that the enumeration list is indented if it should be
broken across lines, and finally apply |make_underlined| to the first
enumeration constant. We want to be sure that indentation is balanced, so we
introduce |indent| and |outdent| in the same rule; it is for this reason that
rule~107 cannot use the |struct_head| in the reduction itself.

@< Rules @>=
{105, {{struct_like, lbrace, expression},-1}, {struct_head, "_B_"}},	@/
{106, {{struct_like, expression, lbrace, expression},-1},
		{struct_head, "_~_B_"}},				@/
{107, {{struct_head, expression, comma, expression},1},
		{expression, "__B!_"}},					@/
{108, {{struct_head, expression, rbrace}}, {int_like, "_~+!_-B_"}},	@/
{109, {{struct_head, expression, comma, rbrace}},
					{int_like, "_~+!__-B_"}},	@[@]

@ The following rules are added to allow short structure and union
specifiers to be kept on one line without having to repeatedly specify \:+.
The idea is to place \:; after the left brace; this will cause the rules
below to be invoked instead of those above, which avoids introducing forced
line breaks.

@< Rules @>=
{110, {{struct_like, lbrace, magic}}, {short_struct_head, "_B__+"}},	@/
{111, {{struct_like, expression, lbrace, magic}},
		{short_struct_head, "_~!_B__+"}},			@/
{112, {{struct_like, int_like, lbrace, magic}},
		{short_struct_head, "_~!$_B__+"}, no_plus_plus},	@/
@q $ emacs-cookie @>
{112, {{struct_like, int_like, lbrace, magic}},
		{short_struct_head, "_~!_B__+"}, only_plus_plus},	@/
{113, {{short_struct_head, declaration}}, {short_struct_head, "_B_"}},	@/
{114, {{short_struct_head, rbrace}}, {int_like, "_-B_"}},		@[@]


@*2 Statements.
Rule~120 gives the usual way statements are formed, while rule~121 handles the
anomalous case of an empty statement. Its use can always be avoided by using
an empty pair of braces instead, which much more visibly indicates the
absence of a statement (e.g., an empty loop body); when the empty statement
is used however, it will either be preceded by a space or start a new line
(like any other statement), so there is always some distinction between a
|while| loop with empty body and the |while| that ends a |do|~statement. A
rule like this with left hand side of length~1 makes the corresponding
category (viz.~|semi|) ``unstable'', and can only be useful for categories
that usually are scooped up (mostly from the left) by a longer rule.  Rules
122--124 make labels (ordinary, case and default), and rules 125~and~126 attach
the labels to statements. Rule~127 makes \:; behave like an invisible
semicolon when it does not match any of the rules designed for it, for
instance if it follows an expression.

@< Rules @>=
{120, {{expression, semi}},		{statement, NULL}},	@/
{121, {{semi}},				{statement, NULL}},	@/
{122, {{expression, colon}},		{label, "!_h_"}},	@/
{123, {{case_like, expression, colon}},	{label, "_ _h_"}},	@/
{124, {{case_like, colon}},		{label, "_h_"}},	@/
{125, {{label, label}},			{label, "_B_"}},	@/
{126, {{label, statement}}, {statement, "b_B_"},not_all_stats_forced},	@/
{126, {{label, statement}}, {statement, "b_f_"},all_stats_forced},	@/
{127, {{magic}},				{semi, NULL}},		@[@]

@ The following rules format compound statements and aggregate initialisers.
Rules 130--134 combine declarations and statements within compound statements.
A newline is forced between declarations by rule~130, unless the declarations
are local (preceded by a left brace) and `\.{+m}' was specified (rule~131);
this rule does not apply to structure specifiers, because the left brace
will already have been captured in a |struct_head| before the rule can match.
If `\.{+f}'~or~`\.{+a}' was specified, then a newline is forced between
statements as well (rule~133). Between the declarations and statements some
extra white space appears in ordinary \Cee\ (rule~132), but not in \Cpp, where
declarations and statements may be arbitrarily mixed (rule~134).  Rules
135--137 then build compound statements, where the last case is the unusual
one where a compound statement ends with a declaration; empty compound
statements are made into simple statements so that they will look better
when used in a loop statement or after a label. If compound statements are
not engulfed by a conditional or loop statement (see below) then they decay
to ordinary statements by rule~138. Rules 139~and~140 reduce aggregate
initialiser expressions, where the reduction of comma-separated lists of
expressions is already handled by the expression syntax.

@< Rules @>=
{130, {{declaration, declaration}}, {declaration, "_f_"}},		@/
{131, {{lbrace, declaration, declaration},1},
				  {declaration, "_B_"},merged_decls},	@/
{132, {{declaration, statement}},  {statement, "_F_"},no_plus_plus},	@/
{132, {{declaration, statement}},
         {statement, "_f_"},only_plus_plus|forced_statements},  	@/
{132, {{declaration, statement}},
         {statement, "_B_"},only_plus_plus|no_forced_statements},	@/
{133, {{statement, statement}},{statement, "_f_"},forced_statements},	@/
{133, {{statement, statement}},{statement, "_B_"},no_forced_statements},@/
{134, {{statement, declaration}}, {declaration, "_f_"},only_plus_plus},	@/
{135, {{lbrace, rbrace}},	  {statement, "_,_"}},			@/
{136, {{lbrace, statement, rbrace}},
	{compound_statement, "ft_+_-f_"},standard_braces},		@/
{136, {{lbrace, statement, rbrace}},
	{compound_statement, "_+f_-f_"},unaligned_braces},		@/
{136, {{lbrace, statement, rbrace}},
	{compound_statement, "f_+f_-f_"},wide_braces},			@/
{137, {{lbrace, declaration, rbrace}},
	{compound_statement, "ft_+_-f_"},standard_braces},		@/
{137, {{lbrace, declaration, rbrace}},
	{compound_statement, "_+f_-f_"},unaligned_braces},		@/
{137, {{lbrace, declaration, rbrace}},
	{compound_statement, "f_+f_-f_"},wide_braces},			@/
{138, {{compound_statement}},			{statement, "f_f"}},	@/
{139, {{lbrace, expression, comma, rbrace}},	{expression, "_,__,_"}},@/
{140, {{lbrace, expression, rbrace}},		{expression, "_,_,_"}},	@[@]

@ Like for structure and union specifiers, we allow compound statements to
be kept on one line by inserting \:; after the left brace. Such statements
will reduce to |statement| rather that to |compound_statement|, so that they
will be treated as if they were simple statements.

@< Rules @>=
{150, {{lbrace, magic}},		{short_lbrace, "__+"}},	@/
{151, {{short_lbrace, declaration}},	{short_lbrace, "_B_"}},	@/
{152, {{short_lbrace, statement}},	{short_lbrace, "_B_"}},	@/
{153, {{short_lbrace, rbrace}},		{statement, "_-B_"}},	@[@]

@*2 Selection, iteration and jump statements.
There are three intermediate categories involved in the recognition of
conditional statements. The category |if_like| stands for `|if|' or an
initial segment of a repeated if-clause, up to and including `|else|~|if|'.
An |if_head| is an |if_like| followed by its (parenthesised) condition
(rules 160~and~161). If the statement following the condition is followed by
`|else|~|if|', the whole construct reduces to |if_like| (so that the
indentation will not increase after the second condition, rules 162~and~163),
otherwise, if only `|else|' follows, reduction is to an |if_else_head|
(rules 164~and~165), and finally, if no |else| follows at all, we reduce with
only the if-branch to |statement| (rules 166~and~167). The reduction rules for
|if_else_head| differ from those for |if_head| in that it will not combine
with an |else|, even if it is present; the formatting is identical to that
of an |else|-less |if_head| (rules 168~and~169).  (It might be tempting to
replace rules 166~and~167 by a reduction from |if_head| to |if_else_head| to
be applied if no matching `|else|' is found, but that would require some
subtle measures to prevent this decay at times when the right context is
insufficiently reduced to decide whether an `|else|' is present or not.) The
formatting of the if and else branches depends on whether they are compound
statements or some other kind of statement (possibly another conditional
statement), and on the flags for statement forcing and brace alignment.

@< Rules @>=
{160, {{if_like, expression}},	      {if_head, "f_~_"}},		@/
{161, {{lbrace,if_like,expression},1}, {if_head, "_~_"},standard_braces},@/
{162, {{if_head, compound_statement, else_like, if_like}},
			{if_like, "__f_~_"},aligned_braces},		@/
{162, {{if_head, compound_statement, else_like, if_like}},
			{if_like, "_~_~_~_"},unaligned_braces},		@/
{163, {{if_head, statement, else_like, if_like}},
			{if_like, "_+B_-f_~_"},not_all_stats_forced},	@/
{163, {{if_head, statement, else_like, if_like}},
			{if_like, "_+f_-f_~_"},all_stats_forced},	@/
{164, {{if_head, compound_statement, else_like}},
			{if_else_head, "__f_"},aligned_braces},		@/
{164, {{if_head, compound_statement, else_like}},
			{if_else_head, "_~_~_"},unaligned_braces},	@/
{165, {{if_head, statement, else_like}},
			{if_else_head, "_+B_-f_"},not_all_stats_forced},@/
{165, {{if_head, statement, else_like}},
			{if_else_head, "_+f_-f_"},all_stats_forced},	@/
{166, {{if_head, compound_statement}},
			{statement, "__f"},aligned_braces},		@/
{166, {{if_head, compound_statement}},
			{statement, "_~_f"},unaligned_braces},		@/
{167, {{if_head, statement}},
			{statement, "_+B_-f"},not_all_stats_forced},	@/
{167, {{if_head, statement}},
			{statement, "_+f_-f"},all_stats_forced},	@/
{168, {{if_else_head, compound_statement}},
			{statement, "__f"},aligned_braces},		@/
{168, {{if_else_head, compound_statement}},
			{statement, "_~_f"},unaligned_braces},		@/
{169, {{if_else_head, statement}},
			{statement, "_+B_-f"},not_all_stats_forced},	@/
{169, {{if_else_head, statement}},
			{statement, "_+f_-f"},all_stats_forced},	@[@]

@ The following rules prevent forced line breaks from conditional statements
that occur within a one-line compound statement.

@< Rules @>=
{170, {{short_lbrace, if_like, expression},1},
		{if_head, "_~_"}},				@/
{171, {{short_lbrace, if_head, statement, else_like}},
		{short_lbrace, "_B_B_B_"}},			@/
{172, {{short_lbrace, if_head, statement}},
		{short_lbrace, "_B_B_"}},			@[@]

@ Switch and loop statements make use of the syntax for conditionals by
reducing to |if_else_head| which will take one further statement and indent it
(rules 180, and rule~181 which avoids a break before a loop at the beginning
of a compound statement in the |standard_braces| style). Recall that `|for|'
and `|switch|' are both |while_like|; the parenthesised object following
`|for|' looks like nothing we have seen before, however, so we need extra
rules to come to terms with it (rules 182--184). Rule~182 is needed to avoid a
line break when these are normally inserted between statements, and rule~184
is needed in case the third expression is empty. For \Cpp\ one needs rule~185
that is similar to rule~182, but this time to preempt the line break that
would otherwise be inserted by rule~132 in case of a |for| loop declaring the
loop variable.

The |do|-|while| loops have to be treated separately. Because we want to
distinguish the case of a loop body that is a |compound_statement| from other
kinds of statements, we cannot wait until the |while| combines with the loop
control condition to an |if_else_head|, since by then a |compound_statement|
will have decayed to |statement|. Hence we pick up the unreduced `|while|'
token and form a new category |do_head| (rules 186~and~187); in case of a
compound statement the `|while|' will be on the same line as the closing
brace. Rules 188~and~189 then combine this with the condition and the
ridiculous mandatory semicolon at the end to form a |statement|.

@< Rules @>=
{180, {{while_like, expression}}, {if_else_head, "f_~_"}},		@/
{181, {{lbrace, while_like, expression},1},
			{if_else_head, "_~_"},standard_braces},		@/
{182, {{lpar, statement, statement}, 1},
			{statement, "_B_"}, forced_statements},		@/
{183, {{lpar, statement, expression, rpar}},	{expression, "__B__"}},	@/
{184, {{lpar, statement, rpar}},	{expression, NULL}},		@/
{185, {{lpar, declaration, statement}, 1},
			{statement, "_B_"}, only_plus_plus},		@/
{186, {{do_like, compound_statement, while_like}},
			{do_head, "__~_"},standard_braces},		@/
{186, {{do_like, compound_statement, while_like}},
			{do_head, "_~_~_"},unaligned_braces},		@/
{186, {{do_like, compound_statement, while_like}},
			{do_head, "__f_"},wide_braces},			@/
{187, {{do_like, statement, while_like}},
			{do_head, "_+B_-B_"},not_all_stats_forced},	@/
{187, {{do_like, statement, while_like}},
			{do_head, "_+f_-f_"},all_stats_forced},		@/
{188, {{do_head, expression, semi}}, {statement, "f_~__f"}},		@/
{189, {{lbrace, do_head, expression, semi},1}, {statement, "_~__f"}},	@[@]

@ The following rules prevent forced line breaks from loop statements
that occur within a one-line compound statement. Since no special layout is
required between the heading of a |while| loop and its body, rule~200
incorporates the heading as if it were a separate statement. For a
|do|-|while| loop we must take a bit more effort to get the spacing
following the |while| correct.

@< Rules @>=
{200, {{short_lbrace, while_like, expression}},
					{short_lbrace, "_B_~_"}},	@/
{201, {{short_lbrace, do_like, statement, while_like},1},
					{do_head, "_B_B_"}},		@/
{202, {{short_lbrace, do_head, expression, semi}},
					{short_lbrace, "_B_~__"}},	@[@]

@ The tokens `|goto|', `|continue|', `|break|', and `|return|' are all
|return_like|; although what may follow them is not the same in all cases,
the following two rules cover all legal uses. Note that rule~211 does not
wait for a semicolon to come along; this may lead to a premature match as in
`|return a+b;|', but this does not affect formatting, while the rule allows
saying things like `|return home|' in a module name (or elsewhere) without
risking irreducible scraps.

@< Rules @>=
{210, {{return_like, semi}},	  {statement, NULL}},		@/
{211, {{return_like, expression}}, {expression, "_~_"}},	@[@]

@*2 Function definitions and external declarations. Apart from the initial
specification of the result type (which is optional, defaulting to |int|), a
new-style function heading will parse as an |function_head| (see the
declaration syntax above), while an old-style function heading is an
|expression| possibly followed by a |declaration| (specifying the function
parameters). Rules 220--222 parse these two kinds of function headings
together with the function body, yielding category |function|. Rule~223
attaches the optional result type specifier. Although the \Cee~syntax requires
that the function body is a compound statement, we allow it to be a
|statement| (to which |compound_statement| will decay), for in case a very
short function body is specified using `\.{\{@@;}'. In \Cpp\ the result type
is obligatory and old-style declarations are forbidden, but on the other hand
there are new syntactic forms such as constructors and destructors; therefore
rules 221--223 are replaced by ones that will be given later. However the case
of a definition of a function without parameters in \Cpp\ looks like  an
old-style function heading with result type, so a rule for this case is used
in place of rule~221 (the corresponding function declaration will be parsed as
an ordinary declaration, which doesn't bother us).

At the outer level declarations and functions can be mixed; when they do a bit
of white space surrounds the functions (rules 224--226), except in \Cpp, where
mixing declarations and definitions is normal in class definitions. Since the
distinction is not made in \Cpp, we might as well let |function| decay to
|declaration| there (rule~224), which will properly handle all combinations,
including those where something like \&{private} or nothing at all follows in
a class declaration. The combination of several declarations is already taken
care of by the syntax for compound statements; no extra white space is
involved there. Rules 227--230 take care of function declarations that are not
definitions (i.e., there is no function body); if followed by a semicolon,
comma, assignment operator or a right parenthesis, the |function_head| decays
to an |expression|, and the rest of the syntax will take care of recognising a
|declaration| or |parameters|. It will also take care of explicit calls of a
default constructor used as an expression, which is why we add the cases of an
|unorbinop|, |langle| or |rangle| as right context. Rules 223~and~227 will be
replaced in~\Cpp, for reasons explained below (incidentally, this is the
reason the category |function_head| was introduced; it used to be simply
|expression|).

@< Rules @>=
{220, {{function_head, statement}}, {function, "!_f_"}},		@/
{221, {{expression, statement}}, {function, "!_f_"},no_plus_plus},	@/
{221, {{int_like,expression, statement}},
                                 {function, "!_ _f_"},only_plus_plus},	@/
{222, {{expression, declaration, statement}},
			       {function, "!_++f_--f_"},no_plus_plus},	@/
{223, {{int_like, function}},	 {function, "_ _"},no_plus_plus},	@/
{224, {{declaration, function}}, {function, "_F_"},no_plus_plus},	@/
{225, {{function, declaration}}, {declaration, "_F_"},no_plus_plus},	@/
{226, {{function, function}},	 {function, "_F_"},no_plus_plus},	@/
{224, {{function}},		 {declaration, NULL},only_plus_plus},	@/
{227, {{function_head, semi},-1},  {expression, NULL},no_plus_plus},	@/
{228, {{function_head, comma},-1}, {expression, NULL}},			@/
{229, {{function_head, binop},-1}, {expression, NULL}},			@/
{229, {{function_head, unorbinop},-1}, {expression, NULL}},		@/
{229, {{function_head, langle},-1}, {expression, NULL}},		@/
{229, {{function_head, rangle},-1}, {expression, NULL}},		@/
{230, {{function_head, rpar},-1},  {expression, NULL}},			@[@]

@*2 Module names.
Although module names nearly always stand for statements, they can be made
to stand for a declaration by appending \:;, or for an expression by
appending `\.{@@;@@;}'. The latter possibility is most likely to be useful
if the module stands for (part of) an initialiser list. A module name can
also be made into an expression by enclosing it in \:[ and~\:], but in that
case rule~230 will apply first, placing a forced break after the module
name. Rules 241, 244,~and~245 prevent a module name from generating forced
breaks if it occurs on a one-line compound statement or structure or union
specifier, while rules 247~and~248 serve to prevent rules 243~and~244 from
matching with priority over rule~246. The rules given here will be replaced by
other ones in compatibility mode.

@< Rules @>=
{241, {{mod_scrap}},			{statement, "_f"},cwebx},	@/
{242, {{short_lbrace, mod_scrap},1},	{statement, NULL},cwebx},	@/
{243, {{mod_scrap, magic}},		{declaration, "f__f"},cwebx},	@/
{244, {{lbrace, mod_scrap, magic},1},
			{declaration, "__f"},cwebx|standard_braces},	@/
{245, {{short_lbrace, mod_scrap, magic},1}, {declaration, NULL},cwebx},	@/
{246, {{short_struct_head, mod_scrap, magic},1},
					    {declaration,NULL},cwebx},	@/
{247, {{mod_scrap, magic, magic}},	{expression, NULL},cwebx},	@/
{248, {{lbrace, mod_scrap, magic, magic},1},
			{expression, NULL},cwebx|standard_braces},	@/
{249, {{short_lbrace, mod_scrap, magic, magic},1},
			{expression, NULL},cwebx},			@[@]

@*2 Additional rules for compatibilty mode.
@^Levy/Knuth \.{CWEB}@>
Although our grammar differs completely from the one used in \LKC., we use
most of it also in compatibility mode (the exception is formed by the rules
concerning module names). We do add a few rules in compatibility mode, mostly
do deal with circumstances that are different for some reason or other.

We start with module names, which behave in a completely different way. In
compatibility mode, as in \LKC., a module name normally stands for an
expression (rule~244) and in practice is almost always followed by a visible
or invisible (|magic|) semicolon. Rules 240~and~241 treat these cases
explicitly, in order to insert a forced break after the semicolon; rule~241
for the case of an invisible semicolon is needed because if we would wait for
the |magic| semicolon to decay to an ordinary one, it might instead combine
with an |int_like| token following it. Rules 242~and~243 are provided to allow
the short form of compound statements even in compatibility mode (even though
it is not present in \LKC.): they preempt rules 240~and~241, avoiding the
forced break. Since in compatibility mode one has no means of indicating that
a module name stands for a set of declarations, we add rule~245 to allow them
nevertheless to be used before a function definition.

Rules 250~and~251 compensate for the fact that compound assignment operators
like `|+=|' are scanned as two tokens in compatibility mode (see
section@#truly stupid@> for an explanation why this is done).
Rule~252 allows types to be used in the argument lists of macros, without
enclosing them between \:[~and~\:], in compatibility mode; this is done
frequently in the Stanford GraphBase. @^Stanford GraphBase@> It is sufficient
to remove expressions from the beginning of the argument list, since types,
and more generally types followed by declarators, are already removed by the
standard rules for |parameters|. As a result the argument list will either
reduce to an |expression| or to |parameters|, depending on whether the final
item was an expression. In both cases it will combine with the macro name to
an |expression|, although the spacing will be a bit too wide in the
|parameters| case. But then, one ought to use \:[~and~\:] anyway, which avoids
this problem.

@< Rules @>=
{241, {{mod_scrap, semi}}, 	{statement, "__f"},compatibility},	@/
{242, {{mod_scrap, magic}}, 	{statement, "__f"},compatibility},	@/
{243, {{short_lbrace, mod_scrap, semi},1},
				{statement, NULL},compatibility},	@/
{244, {{short_lbrace, mod_scrap, magic},1},
				{statement, NULL},compatibility},	@/
{245, {{mod_scrap}},        	{expression, NULL},compatibility},	@/
{246, {{statement, function}},	 {function, "_F_"},compatibility},	@/
@)
{250, {{binop, binop}},			{binop,"r__"},compatibility},	@/
{251, {{unorbinop, binop}},		{binop,"r__"},compatibility},	@/
{252, {{lpar, expression, comma}}, {lpar, "___p3"}, compatibility},	@[@]
@[@]

@*2 Additional rules for \Cpp.
Up to this point we have included some specific rules for \Cpp, in places
where a slight deviation from the \Cee~syntax was required. There are however
a large number of syntactic possibilities of \Cpp\ that are not even remotely
similar to those of~\Cee, so it is most convenient to collect them in a
separate section.

We start with rules for `\&{operator}', which are simple: it should combine
with a following operator symbol of any type, or |subscript|, or |lpar, rpar|,
to form an expression (rules 260--262). We force the operators to be
interpreted as ordinary symbols in math mode for proper spacing (except for
|unop|, operators which already get proper spacing); the argument
for the ``mathord'' control is enclosed in braces in case of a |binop| token,
since this might be a compound operator like |+=| (this could cause problems
is a comment follows the operator token; don't do that). Then rules 263--266
take care of the `::'~operator: either a class name or an expression (a
namespace identifier) or nothing is expected at the left, and either an
ordinary or a class identifier at the right; the resulting category is that of
the right hand side. However, we must make sure that the contraction of
|int_like| scraps does not absorb class name used before `::' into a possibly
preceding type name; therefore rule~266 doubles rule~265 in case of a
preceding |int_like| scrap.

@< Rules @>=
{260, {{case_like, binop}},	{expression, "_o{_}"},only_plus_plus},	@/
{260, {{case_like, unorbinop}},	{expression, "_o_"},only_plus_plus},	@/
{260, {{case_like, unop}},	{expression, NULL},only_plus_plus},	@/
{260, {{case_like, langle}},	{expression, "_o_"},only_plus_plus},	@/
{260, {{case_like, rangle}},	{expression, "_o_"},only_plus_plus},	@/
{261, {{case_like, subscript}},	{expression, NULL},only_plus_plus},	@/
{262, {{case_like, lpar, rpar}}, {expression, NULL},only_plus_plus},	@/
{263, {{colcol, expression}},	{expression, NULL},only_plus_plus},	@/
{263, {{colcol, int_like}},	{int_like, NULL},only_plus_plus},	@/
{264, {{expression, colcol, expression}},
				{expression, NULL},only_plus_plus},	@/
{264, {{expression, colcol, int_like}},
				{int_like, NULL},only_plus_plus},	@/
{265, {{int_like, colcol, expression}},
				{expression, NULL},only_plus_plus},	@/
{265, {{int_like, colcol, int_like}},
				{int_like, NULL},only_plus_plus},	@/
{266, {{int_like, int_like, colcol, expression},1},
				{expression, NULL},only_plus_plus},	@/
{266, {{int_like, int_like, colcol, int_like},1},
				{int_like, NULL},only_plus_plus},	@[@]

@ Type identifiers may appear as the left hand side of an assignment within a
list of formal parameters, indicating a default argument; in this case the
whole assignment should behave as a type identifier (rule~270). Rules 271
and~272 serve to absorb the base class in a class definition into the name of
the class being defined. Rules 273--274 deal with template definitions. What
comes inside the angle brackets of such a definition looks much like what is
inside such brackets in an instantiation of a template class, so rules 75--78
may reduce the whole, angles included, to |templ_params|. We can however also
expect parameters specified as `\&{class X}' or `\&{typename X}'; since
rule~94 requires the presence of a semicolon for \Cpp, which is not present in
the case considered here, we must provide additional rules that will fit such
parameters into the logic of rules 75--78, which is what rules 273~and~274 do.
Also non-type template argument look like ordinary function arguments, and
rules 275~and~276 include these.

While we are at it, we deal with the other use of `\&{typename}', namely to
make clear that the following expression is a type. To this end we reduce that
keyword and the following |expression| or |int_like| scrap, probably produced
by a reduction involving |colcol|, to an |int_like| scrap (rules 277~and~278),
but only if an |expression| follows since otherwise this rule would match too
soon, when some |colcol| still has to be reduced. The recognition of the
template definition heading is then simple: the keyword `\&{template}' should
be followed by such an angled group, and then a valid declaration; we reduce
to |declaration| forcing a line break before the declaration and indenting it
by one notch (rule~279).

@< Rules @>=
{270, {{int_like, binop, expression}},
				{int_like, NULL},only_plus_plus},	@/
{271, {{int_like, colon, case_like, int_like}},
                                {int_like,"_m__~_"},only_plus_plus },	@/
{272, {{int_like,colon, int_like}}, {int_like,"_m__"},only_plus_plus },	@/
{273, {{struct_like, int_like, comma},-1},
				{int_like,"_~_"},only_plus_plus},	@/
{274, {{struct_like, int_like, rangle},-1},
				{int_like,"_~_"},only_plus_plus},	@/
{275, {{langle, int_like, expression, comma}},
				{langle,"__~!__p5"},only_plus_plus},	@/
{276, {{langle, int_like, expression, rangle}},
			{templ_params,"a__~!_a_"},only_plus_plus},	@/
{277, {{struct_like,int_like,expression},-1},
				{int_like,"_~_"},only_plus_plus},	@/
{278, {{struct_like,expression,expression},-1},
				{int_like,"_~_"},only_plus_plus},	@/
{279, {{template_like,templ_params,declaration}},
				{declaration,"__f+_-"},only_plus_plus},	@[@]

@ Next we give rules catering with constructor declarations in class
definitions. First of all we must recognise the fact that the class name is
being used as a function name here; the simplest solution is to recognise the
combination of an |int_like| followed by a parameter list (rule~280). We
cannot let a |function_head| (possibly created by rule~280) decay to an
|expression| when followed be a semicolon, as we do for~\Cee, since
declarations of constructor and destructor members of a class lack an initial
type specification, so the |expression| would fail to become part of a
|declaration|. Therefore, special measures are necessary: the simplest
solution is to simply absorb (rules 281~and~282) any preceding type specifier
into the |function_head| (thereby removing the distinction between its
presence or absence), and construct de |declaration| explicitly from the
|function_head| and the following semicolon (rule~283). This absorption is
only allowed when a semicolon does follow (or a statement, in case of a
function definition rather than declaration), since when followed by a comma
or right parenthesis it would perturb the recognition of parameter lists
containing function (pointer) parameters). We must also absorb |int_like| from
the other side into a |function_head| (this time independently of what
follows), since there may be a |const| after the parameters in declarations of
member functions (rule~284); we need a special rule here to handle
parameterless function, since the function name followed by empty parentheses
will have been parsed as an |expression| rather than as a |function_head| in
this case (rule~285).

These rules do create another problem: while inside the declaration of a
\&{class}~\&{example} the sequence ``\&{example}$(\,);$'' is a declaration,
indicating the presence of a default constructor in the class, the same
sequence is an expression standing for an explicit call of that default
constructor in other contexts, such as ``|example x=@[example()@];|''. We
already explained that we do not want to reduce the |function_head| to
|expression| based on the presence of the following semicolon, so we shall try
to demand such a reduction based on the left context. One such context is a
binary operator as in the declaration example given or in an actual
assignment, which is taken care of by rule~286. Usage in the arguments of
functions or followed by operators will already be handled by rules 228--230.
One use that has to be considered explicitly is in a |return| statement, which
is what rule~287 below does.

Rule~288 allows member object initialisers to be specified between the
function head of a constructor and its body; we insert a penalty before the
colon lower than the one that is emitted after commas, so the list of
initialisers will have tendency to stay together on the line following the
function heading (unless it will fit on the same line. The increase of
indentation level given by this rule will usually lead to {\it less\/} hanging
indentation, due to the way \TeX\ and our macros handle this: at every
|indent| the indentation level is increased and set as the current hanging
indentation, but at an |outdent| the level is decreased but not set as hanging
indentation; as a result the hanging indentation level, which is only taken
into account by \TeX\ at the next forced break, is usually the maximum level
reached since the preceding forced break. But here we have descended three
levels since our maximum inside the |function_head|, and doing an |indent|
will reset the effective hanging indentation to two levels below that.
Rule~289 caters for calls of
\&{delete}[].

@< Rules @>=
{280, {{int_like, parameters}},	{function_head, "!__"},only_plus_plus},	@/
{281, {{int_like, function_head, semi},-1},
				{function_head, "_~_"},only_plus_plus},	@/
{282, {{int_like, function_head, statement},-1},
				{function_head, "_~_"},only_plus_plus},	@/
{283, {{function_head, semi}},	{declaration, NULL},only_plus_plus},	@/
{284, {{function_head, int_like}},
				{function_head, "_ _"},only_plus_plus},	@/
{285, {{int_like,expression,int_like},1},
                                {function_head,"_ _"},only_plus_plus},	@/
{286, {{expression, binop, function_head}},
				{expression,NULL},only_plus_plus},	@/
{287, {{return_like, function_head}},
				{expression,"_~_"},only_plus_plus},	@/
@)
{288, {{function_head,colon,expression,lbrace},-1},
			{function_head,"_+p1m__-"},only_plus_plus},	@/
{289, {{sizeof_like,subscript}},{sizeof_like,NULL},only_plus_plus},	@[@]

@ Now we get to rules for exception handling.
Rules 290 and~291 treat a simple \&{throw} call just like a |return|
statement, in particular, we do not wait for a semicolon for the same reasons
as mentioned for rule~211. Rule~292 handles the use of \&{throw} at the end of
a |function_head| to indicate the list of exceptions that can be thrown when
the function is called; the resulting |int_like| will treated just like a
|const| by rule~284. Rule~@93 handles an isolated \&{throw} that re-throws the
last exception caught. We treat \&{try} just like for instance
``|while(busy)|'', in other words, we let the rules for an |if_else_head| do
the real work. The same goes for \&{catch} followed by the exception to catch,
the latter will parse as |parameters|.

@< Rules @>=
{290, {{throw_like,expression}},  {expression,"_~_"},only_plus_plus},	@/
{291, {{throw_like,function_head}}, {expression,"_~_"},only_plus_plus},	@/
{292, {{throw_like,parameters}},  {int_like,NULL},only_plus_plus},	@/
{293, {{throw_like,semi}},   	  {statement,NULL},only_plus_plus},	@/
{295, {{try_like}},		  {if_else_head,"f_"},only_plus_plus},	@/
{296, {{catch_like,parameters}},  {if_else_head,"f__"},only_plus_plus},	@[@]