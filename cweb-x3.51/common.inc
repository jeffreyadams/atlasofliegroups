@ This is is the beginning of material included (by means of `\.{@@i}') from
the file \.{common.inc}, so that it appears in an identical way in the
programs for \.{CTANGLE}, \.{CWEAVE}, and in the code shared between them
(from the file \.{common.w}). There will therefore not be any code producing
material, just like in a header file. However, the major part of the
declarative information that is shared between all the compilation units is
contained in an actual header file \.{common.h}, which is written as an
auxiliary output file when \.{CTANGLE} processes \.{common.w}; in particular
that file contains all the typedefs, |extern| variable declarations and
function prototype declarations relating to code shared by \.{CTANGLE} and
\.{CWEAVE}. Therefore this file is almost exclusively devoted to some macro
definitions of general utility and a few function prototype declarations for
functions that are {\it used\/} by the shared code but defined separately in
\.{CTANGLE} and \.{CWEAVE}.

@h <stdlib.h>
@h <stdio.h>
@h <string.h>
@h <ctype.h>
@h <limits.h>
@h "common.h"

@ We start with some dimensioning parameters that are used both in
\.{CTANGLE} and \.{CWEAVE}. Some of them have been decreased with respect to
their earlier values which were sufficient in the original \.{WEB} to handle
\TeX. Nevertheless they should be sufficient for most applications of
\.{CWEB}, since in \Cee\ there is no need to generate large programs
entirely as a single file (and \TeX\ is not written in \.{CWEB}!).  The
smaller values allow \.{CWEB} to run on rather moderate size computers.

@d max_bytes 50000L
   /* the number of bytes in identifier and module names */
@d max_modules 1000 /* maximal number of module names */
@d max_idents 5000 /* maximal number of identifiers */
@d max_sections 4000
   /* greater than the total number of sections, and less than $10240$ */
@d hash_size 353 /* hash modulus, preferably odd */
@d buf_size 200 /* maximum length of input line, plus two */
@d longest_name 1000
   /* module names and strings shouldn't be longer than this */
@d long_buf_size (buf_size+longest_name) /* for \.{CWEAVE} */

@ Here are some macros of general utility. We use |local| in place of
|static| for functions and variables when we want to stress the property
of file-scope rather than static initialisation and permanence of value.
The macro |array_size| can be used to compute the number of elements in a
statically initialised array. For truth values we use the type |boolean|
and the values |true| and |false|. We also make the general convention that
whenever we mention a pointer to the ``end'' of a subsequence of a linear
array, we mean the address of first entry beyond the subsequence itself.
This is in keeping with general practice in~\Cee, and the language definition
in fact guarantees the existence of the address even if the subsequence
should contain the last entry of the array; this convention avoids repeated
use of turgid phrases such as ``pointing one place beyond the last entry
of~\dots''.

@d local static
@f local static
@d array_size(a) ((int)(sizeof(a)/sizeof(a[0])))
@d false (boolean) 0
@d true (boolean) 1
@d ctangle 0
@d cweave 1

@ Although \.{CWEB} uses the conventions of \Cee\ programs found in
the standard \.{<ctype.h>} header file, it does assume that the character set
is the \caps{ASCII} code. This dependency is mild however, and limited to
the assumption that certain character codes below |040| are not occupied by
ordinary characters. To be able to use such vacant spots in a character
code independent way is possible (as is done in \TeX) by mapping all
characters on input to their \caps{ASCII} positions (thereby ensuring that
certain positions remain unused) and unmapping them on output; such an
approach was deemed too tedious for the \.{CWEB} system however. Rather,
the few places that have to be modified for non-\caps{ASCII} codes can be
located by looking up the entry named ``\caps{ASCII} code dependencies''
in the index of each of the programs.

@ A few character pairs are encoded internally as single characters, using
the definitions below. These definitions are consistent with an extension of
\caps{ASCII} code originally developed at \caps{MIT} and explained in
Appendix~C of {\sl The \TeX book\/}. Thus, users who have such a character
set can type things like `\.{\char'32}' and `\.{\char'4}' instead of
`\.{!=}' and `\.{\&\&}'; however, their files will not be too portable until
more people adopt the extended code. To be precise, when moved to other
installations their files can still be processed by \.{CTANGLE} and
\.{CWEAVE}, producing the required output, but users will have problems
reading the source files. Therefore it is advised to use the two-character
form in the \.{CWEB} files, which will be converted to single characters on
input.

@^ASCII code dependencies@> @^system dependencies@>

@d and_and 04 /* `\.{\&\&}'; this corresponds to \caps{MIT}'s \.{\char'4} */
@d lt_lt 020 /* `\.{<<}'; this corresponds to \caps{MIT}'s \.{\char'20} */
@d gt_gt 021 /* `\.{>>}'; this corresponds to \caps{MIT}'s \.{\char'21} */
@d plus_plus 013 /* `\.{++}'; this corresponds to \caps{MIT}'s \.{\char'13} */
@d minus_minus 01 /* `\.{--}'; this corresponds to \caps{MIT}'s \.{\char'1} */
@d minus_gt 031 /* `\.{->}'; this corresponds to \caps{MIT}'s \.{\char'31} */
@d not_eq 032 /* `\.{!=}'; this corresponds to \caps{MIT}'s \.{\char'32} */
@d lt_eq 034 /* `\.{<=}'; this corresponds to \caps{MIT}'s \.{\char'34} */
@d gt_eq 035 /* `\.{>=}'; this corresponds to \caps{MIT}'s \.{\char'35} */
@d eq_eq 036 /* `\.{==}'; this corresponds to \caps{MIT}'s \.{\char'36} */
@d or_or 037 /* `\.{\v\v}'; this corresponds to \caps{MIT}'s \.{\char'37} */

@ Invoking the macro |find_char| will fetch a new line if the current one
has been completely read, and it will either return |false| to indicate that
there is no input left, or otherwise ensure that it is safe to inspect
|*loc| (although it might be the line-ending space).

@d find_char() (loc<=limit || get_line())

@ The following declarations give the interface to the searching algorithms.
They are defined separately in \.{CTANGLE} and \.{CWEAVE} and used by the
common code, rather than the other way around. The small number of such
functions does not seem to justify creating a separate header file for them,
and moreover, it would be difficult to decide whether such a header file
should be produced by the source file for \.{CTANGLE} or by the source file
for \.{CWEAVE}.

@< Function prototypes used but not defined in the shared code @>=

boolean names_match (id_pointer,char*,int,int);
void init_id_name (id_pointer,int);
void init_module_name (mod_pointer);

@ Here are some macros associated with the storage of names. The first four
of them convert pointers to structures representing identifiers or module
names to smaller |sixteen_bits| indices into the arrays they are stored in,
and back again. The macros |name_begin|, |name_end| and |length| provide
the basic attributes of the actual string referred to by a name pointer.
The macro |complete_name| tells whether the complete name for a module has
been encountered so far, and |print_mod| prints a module name on the user's
terminal for error reporting, providing a colon, angle brackets and an
ellipsis if necessary.

@d id_index(p) ((sixteen_bits)((p)-id_table))
@d id_at(i)    (&id_table[i])
@d mod_index(p) ((sixteen_bits)((p)-mod_table))
@d mod_at(i)    (&mod_table[i])
@d name_begin(p) ((p)->byte_start)
@d length(p) ((int)(strlen(name_begin(p))))
@d name_end(p) (name_begin(p)+length(p))
@d complete_name(p) ((p)->byte_start[-1]=='\0')
@d print_mod(p)
   printf(": <%s%s>",name_begin(p), complete_name(p) ? "" : "..." )

@ Three levels of error severity are distinguished: informative, serious and
fatal. An overflow stop occurs if \.{CWEB}'s tables aren't large enough.
Sometimes the program will detect a violation of one of its supposed
invariants (at least this could happen if the program still contains some
errors), and \.{CWEB} then prints an error message that is really for the
\.{CWEB} maintenance person, not the user. In such cases the program says
|confusion("indication of where we are")|.

@d spotless 0 /* |history| value for normal jobs */
@d harmless_message 1 /* |history| value when non-serious info was printed */
@d error_message 2 /* |history| value when an error was noted */
@d fatal_message 3 /* |history| value when we had to stop prematurely */
@d mark_harmless()
   @+ if (history==spotless) history=harmless_message; @+ else @;
@d mark_error() (history=error_message)
@d overflow(t) fatal("\n! Sorry, %s capacity exceeded",t)
			@.Sorry,... capacity exceeded@>
@d confusion(s) fatal("\n! This can't happen: %s",s) @.This can't happen@>

@ Command line flag settings are stored in an array |flags|. For some of
them we use symbolic names.

@d show_banner flags['b'] /* should the banner line be printed? */
@d show_happiness flags['h'] /* should lack of errors be announced? */
@d show_progress flags['p'] /* should progress reports be printed? */
@d show_stats flags['s'] /* should statistics be printed at end of run? */
@d C_plus_plus flags['+'] /* is the language `\Cpp' rather than `\Cee'? */
@d compatibility_mode flags['c'] /* emulate \LKC.? */

@ Here are some macros for terminal output. The macro |update_terminal| is
invoked when we want to make sure that everything we have output to the
terminal so far has actually left the computer's internal buffers and been
sent. Note that the cast on the |leng| parameter of |term_write| is
necessary because a pointer difference need not be of type |int|.

@d update_terminal() fflush(stdout) /* empty the terminal output buffer */
@d new_line() putchar('\n')
@d term_write(string,leng) printf("%.*s",(int)(leng),string) /* write on the standard output */
