% This file is part of CWEBx.
% This program by Marc van Leeuwen based on earlier versions by
% D. E. Knuth., Silvio Levy and Frank Jensen.
% It is distributed WITHOUT ANY WARRANTY, express or implied.
% CWEB (Revision: 2.0) % Don Knuth, July 1990
% Version 3.x, Marc van Leeuwen, December 1993
% CWEBx 2+1.0, Marc van Leeuwen, August 1994
% CWEBx 3.0, Marc van Leeuwen, Januari 1995
% CWEBx 3.02, Marc van Leeuwen, April 1996
% CWEBx 3.03, Marc van Leeuwen, January 1998
% CWEBx 3.05, Marc van Leeuwen, May 2000
% CWEBx 3.06, Marc van Leeuwen, December 2005
% CWEBx 3.1, Marc van Leeuwen, May 2006
% CWEBx 3.5, Marc van Leeuwen, June 2009

% Copyright (C) 1987,1990 Silvio Levy and Donald E. Knuth
% Copyright 1994-2009 Marc A. A. van Leeuwen

% Permission is granted to make and distribute verbatim copies of this
% document provided that the copyright notice and this permission notice
% are preserved on all copies.

% Permission is granted to copy and distribute modified versions of this
% document under the conditions for verbatim copying, provided that the
% entire resulting derived work is distributed under the terms of a
% permission notice identical to this one.

\def\me.{CWEAVE} \def\myroots{\.{WEAVE}}

% Here is TeX material that gets inserted after \input cwebxmac

\font\bfit=cmbxti10
\def\cwebxsectionnumber#1{{\ifcase \gdepth \bf \or \bfit \or \it \fi#1.}\quad}

\def\xr.{cross-reference}
\def\Cpp{\Cee\PP}
\def\TeXxstring{\TeX\_string}
\def\TeXxlike{\TeX\_like}
\def\skipxTeX{skip\_\TeX}
\def\copyxTeX{copy\_\TeX}

@i intro.inc  % Here is some source text that matches the start of CTANGLE

@d banner "This is CWEAVE (Version "@+version_string@+")"
@q the C macro |version_string| was defined at the end of intro.inc@>

@ The following parameters are specific to \.{\me.}; those which are
common to \.{CTANGLE} and \.{CWEAVE} are defined in the file \.{common.inc}
and appear below. Some of these values have been decreased with respect to
their earlier values which were sufficient in the original |WEB| to handle
\TeX; a motivation is given at the common declarations.

@d max_refs 10000 /* number of \xr.s; must be less than 65536 */
@d max_toks 10000 /* number of symbols in \Cee\ texts being parsed;
  must be less than |65536| */
@d max_texts 2500 /* number of phrases in \Cee\ texts being parsed;
  must be less than |10240| */
@d max_scraps 4000 /* number of tokens in \Cee\ texts being parsed */
@d max_no_of_nodes 300 /* number of nodes in search trie for grammar rules,
  must be at most |65536| */
@d line_length 80 /* maximal line length for \TeX\ output;
  should be less than |256| */
@d stack_size 400 /* number of simultaneous output levels */
@d sort_stack_size 500 /* number of identifier lists during sorting */

@ The program is built from two compilation units, one with source file
\.{common.w}, which contains a collection of routines and data shared
between \.{CTANGLE} and \.{CWEAVE}, and a second with master source file
\.{cweave.w} containing all code specific to \.{\me.}, and whose
typeset version you are now reading. All compilation units of the \.{CWEB}
system incorporate the file \.{common.inc} containing common declarations.

\.{\me.} has the following outline. It operates in three phases: first
it reads the source file once for collecting \xr. data, then it
reads the source file for a second time, meanwhile producing the bulk of
the \TeX\ output file, and finally it outputs the information in its tables
from which \TeX\ will produce the list of module names, index and table
of contents. Syntax errors would often be reported identically in the first
two phases, and problems in printing module names similarly in the second
and third phase; to avoid this we check after each phase if any serious
errors were found, and if so, we stop prematurely. \.{\me.} can optionally
be compiled with preprocessor symbols |DEBUG| and/or |STAT| defined; the
former is useful for the user who requires detailed information about the
parsing process, the latter if one wishes to keep track of how much of the
static resources were actually used.

@c
@< Function prototypes used but not defined in the shared code @>@;
@< Typedef and enumeration declarations @>@;
@< Prototypes @>@;
@< Global variables @>@;

int main (int argc,char** argv)
{ program=cweave;
  make_xrefs=true;
  common_init(argc,argv,banner);
  @< Set initial values @>
  if (show_banner)
    print("%s, in %s mode.\n",banner,C_plus_plus ? "C++" : "C");
    /* print a ``banner line'' */
  @< Store all the reserved words @>
  phase_one (); /* read all the user's text and store the \xr.s */
  if (history>harmless_message) wrap_up(); /* stop in case of trouble */
  open_output_file();
  phase_two (); /* read all the text again and translate it to \TeX\ form */
  if (history>harmless_message) wrap_up(); /* stop in case of trouble */
  phase_three (); /* output the \xr. index */
  wrap_up (); /* and exit gracefully */
  return 0; /* for completeness---not reached */
}

@ The macro |variant| needs to be defined before we include the file
\.{common.h}, in order to affect the definition of |id_info| and |mod_info|
(see \.{common.w} for an explanation); it refers to the |struct xref_info|
that will be declared later.

@d variant @;xref_info
@< Prototypes @>=
void phase_one (void);
  /* read all the user's text and store the \xr.s */
void phase_two (void);
  /* read all the text again and translate it to \TeX\ form */
void phase_three (void); /* output the \xr. index */

@i common.inc


@*1 Data structures exclusive to {\tt \me.}. @:title@>
%
Although the eventual typeset output produced after running \.{\me.} is meant
to closely resemble the \.{CWEB} source from which it was created (or maybe it
is the other way around), the task of \.{\me.} is more complicated than that
of \.{CTANGLE}, because it involves a much more detailed ``understanding'' of
the program fragments present in the source. Therefore the bulk of the program
text is concerned with the parsing process during Phase~II. Detailed
information about these matters shall be given at more appropriate places
below; here we only discuss how individual tokens are tagged on input, and
how \xr. infromation is stored.

The |ilk| field of an |id_info| structure is used to distinguish between
various types of identifiers, as follows:

\yskip\hang |normal| identifiers are part of the \Cee\ program and
will usually appear in italic type.

\yskip\hang |roman| identifiers are index entries that appear after \:\^ in
the \.{CWEB} file.

\yskip\hang |wildcard| identifiers are index entries that appear after \:?
in the \.{CWEB} file.

\yskip\hang |typewriter| identifiers are index entries that appear after \:.
in the \.{CWEB} file.

@:section number example@>
\yskip\hang |reference| identifiers are used for explicit \xr.s by the user
within pieces of \TeX~text, allowing us for instance to remark that this is
section@#section number example@> out of a total of@#index@> in this file,
and that the sections having titles are number@#title@>. These \xr.s marks
are set by~\:: and referred to by~\:\#.

\yskip\hang |header_file_name| is used for strings that have been seen as
name of a scanned header file, so multiple inclusions of the same file can be
avoided

\yskip\hang |type_defined| identifiers have occurred in a |typedef|
declaration, or in a \&{class}, |struct|, |union|, or |enum| declaration
in~\Cpp. They will parse as types like |int|, but \xr.s to
them will be collected.

\yskip\hang |TeX_like| identifiers like |TeX| are typeset as user-definable
control sequences.

\yskip\hang |NULL_like| identifiers like `|NULL|' (\.{NULL}) are also
typeset as (\TeX) control sequences, but in addition are treated like
reserved words for indexing purposes, i.e., only underlined references are
generated.

\yskip\hang |const_like| and similar identifiers are \Cee\ reserved words that
will be typeset in boldface; their |ilk| values all exceed |NULL_like|. Most
of these values are also the category codes for these tokens (defining their
parsing behaviour), and defined in section@#categories@> together with the
other category codes. Some special cases are given here, for tokens that need
a distinguished |ilk| only during lexical analysis: |const_like| and
|typedef_like| will become |int_like| in parsing (and the same holds for
|type_defined| identifiers). The remaining values are used only for \Cpp:
|and_like| and |not_like| become |binop| and |unop| respectively, but also
trigger |reserved|; |namespace_like| and |typename_like| will become
|struct_like|, but remain different during phase~I to (conditionally for
|typename_like|) avoid making the identifier following it |type_defined|.

@f TeX_like TeX

@d reserved(a) (a->ilk>=type_defined)
  /* whether identifier automatically gets this |ilk| */
@d unindexed(a) (a->ilk>=NULL_like)
  /* whether cross-referencing is suppressed */

@< Typedef and enum... @>=
enum @/
{ normal, /* |ilk| of ordinary identifiers */
  roman, /* |ilk| of roman type index entries */
  wildcard, /* |ilk| of user-formatted index entries */
  typewriter, /* |ilk| of typewriter type entries */
  reference, /* |ilk| of identifiers used for explicit \xr.s */
  header_file_name, /* |ilk| of file names seen as included header file */
  type_defined, /* |ilk| of identifiers that are defined by |typedef| */
  TeX_like, NULL_like,
    /* |ilk| of identifiers with user-given control sequences */
  const_like, typedef_like, and_like, not_like, namespace_like, typename_like
     /* special reserved words */
};

@ Besides the memory used for storing names, another large memory area is
used in \.{\me.} for keeping the \xr. data. All uses of the name~|p| are
recorded in a linked list beginning at |p->xref|, which is an |xref_pointer|
pointing into the |xmem| array (here we use the pointer field for additional
information in |id_info| and |mod_info| structures, which was called |equiv|
in~\.{CTANGLE}). The elements of |xmem| are structures consisting of an
integer, |num|, and the index |next| of another element of |xmem|. If
|x==p->xref|, the value of |x->num| is a section number where |p| occurs,
plus a multiple of |cite_flag| that indicates the nature of the occurrence,
or it is |file_flag| marking the fact that |p| is a module name which is the
name of an auxiliary output file (for \.{CTANGLE}). The next such \xr.
for~|p|, if any, is |xmem[x->next]|.

Since the entries of |xmem| are small but numerous, and the \xr. lists are
not traversed very frequently, a significant amount of space is saved by
storing a |sixteen_bits| index as a link to the next node rather than an
|xref_pointer|, and at very little cost in performance. The main price paid
is in terms of elegance, as somewhat different methods of traversal of the
lists are appropriate depending on whether the traversal may lead to the
insertion of a new node into the list or not. If no such insertion is
needed, then we can work with |xref_pointer| values as if we were dealing
with an ordinary linked list, provided only that the pointer to a successor
node is obtained by invoking a macro |next_xref| instead of selecting the
link field. If however we may need to insert a new node into the list, the
the most convenient method is to use a pointer to a link field, which then
must have type |(sixteen_bits*)|; for that case the macros |xnum| and
|xlink|, which select the |num| and |next| fields respectively of a node
specified by its index, and which yield expressions that one can assign to
or take the address of, are more useful than |next_xref|. Since the
beginning of the list is indicated by a |xref_pointer| rather than by its
index, we need to convert that pointer into an index (for which the macro
|xref_index| is supplied) when using the second method, and store this
index in a local variable, not forgetting to update the original
|xref_pointer| in case an insertion occurs at the front of the list.

During collection of \xr.s, the lists are basically in decreasing order by
section number, as a result of the fact that new nodes are added at the
beginning of the list (for module names the situation is actually a bit more
complicated, as described below), so that in that phase the |next| field
really refers to the previous \xr. for the same name; at the end of Phase~I,
however, the lists are reversed.

The global variable |xref_switch| is set either to |def_flag| or to zero,
depending on whether the next \xr. to an identifier is to be underlined or
not in the index. This switch is set to |def_flag| when \:! or \:d or \:f is
scanned, and it is cleared to zero when the next identifier or index entry
\xr. has been made. Similarly, the global variable |mod_xref_switch| is
either |def_flag|, |cite_flag|, or zero, depending on whether a module name
is being defined, cited, or used. During Phase~II a number of \xr.s for
identifiers will be changed into underlined ones, when the identifier is
found to occur in a declaration or function definition, since these
occurrences cannot be reliably recognised in Phase~I.

@d xref equiv_or_xref
@d next_xref(x) (&xmem[(x)->next])
@d xnum(i) (xmem[i].num)
@d xlink(i) (xmem[i].next)
@d xref_index(p) ((sixteen_bits)((p)-xmem))
@d cite_flag 0x4000
       /* a |sixteen_bits| power of 2 that is at least |max_sections+2| */
@d def_flag 0x8000 /* twice that */
@d num_mask (cite_flag-1) /* a bit-mask for working modulo |cite_flag| */

@<Typedef...@>=
typedef struct xref_info
{ sixteen_bits num; /* section number plus a multiple of |cite_flag| */
  sixteen_bits next; /* index of the next \xr. in the list */
} xref_info, *xref_pointer;

@ We use an allocation pointer~|xref_ptr| that indicates how much of |xmem|
is already in use. Unlike in the other cases of sequential allocation, where
the allocation pointer points to the first free position, |xref_ptr| points
to the last occupied position in |xmem|; the first node |xmem[0]| is already
in use when \.{\me.} starts. All allocation of |xref_info| nodes is
performed by calling |make_xref(n,i)|, where |n| is the |num| field of the
new node, and |i| is the index of its successor in the list; after
|make_xref| has been invoked, |xref_ptr| points to the new node.

@d make_xref(n,i)
    /* create \xr. node with |num==n| and successor |xmem[i]| */
  if (++xref_ptr >= &xmem[max_refs])
    overflow ("cross-reference"); @.cross-reference capacity exceeded@>
  else xref_ptr->num=n, xref_ptr->next=i @;

@<Global...@>=
xref_info xmem[max_refs]; /* contains \xr. information */
xref_pointer xref_ptr = &xmem[0]; /* the last used position in |xmem| */
sixteen_bits xref_switch = 0, mod_xref_switch = 0;
	/* either zero or |def_flag| */

@ A first node is initially allocated, with its |num==0|; the node is used
to initialise \xr. lists (see |init_id_name| and |init_module_name| below),
and serves as a sentinel at the end of those lists.

@<Set init...@>=
xnum(0)=0;  /* sentinel node terminating \xr. lists */

@ A new \xr. for an identifier~|p| is formed by calling |new_id_xref(p)|
with |section_count| and |xref_switch| set to the appropriate value.
Multiple references to the same section are merged (with underlining
activated if any reference requires so) while non-underlined references to
one-letter identifiers, reserved words or identifiers with |ilk==NULL_like|
(like `\.{NULL}') are ignored. Explicit \xr.s set by~\:: are never
underlined, so |xref_switch| is ignored in this case. If the user has set
the |no_xref| flag (the \.{-x} option of the command line), |new_id_xref| is
active only for explicit \xr.s; when we are reading a header file included
by~\:h, it is inactive altogether (but since \:: is disabled in such header
files, this case is treated together with the |no_xref|~case).

@d no_xref (!flags['x'])
@d make_xrefs flags['x'] /* should cross references be output? */

@c
void new_id_xref (id_pointer p)
{ sixteen_bits f=xref_switch; xref_switch=0;
  if (p->ilk==reference) f=0;
  else if (f==0 && (unindexed(p) || length(p)==1)
        || no_xref || including_header_file) return;
  if ((p->xref->num&num_mask)==section_count) p->xref->num|=f;
  else
  {@; make_xref(section_count|f,xref_index(p->xref)); p->xref=xref_ptr; }
}

@ The \xr. lists for module names are slightly different. Suppose that a
module name is defined in sections $d_1$,~\dots,~$d_k$, cited in sections
$c_1$,~\dots,~$c_l$, and used in sections $u_1$,~\dots,~$u_m$, where the
sequences of $d$'s, $c$'s and $u$'s each are in increasing order. Then its
list will contain |@t$d_k$@>+def_flag|, \dots, |d1+def_flag|,
|@t$c_l$@>+cite_flag|, \dots, |c1+cite_flag|, $u_m$, \dots, $u_1$,~$0$, in
this order; the final~$0$ is the sentinel node, which allows the loops below
to be a little faster. If the module name specifies an output file (with the
\:( feature) then a node containing the special value |file_flag| is
prepended to this list. The special ordering described here only serves for
efficiency of insertion, and after Phase~II the order will be adjusted to
the more natural sequence |d1+def_flag|, \dots, |@t$d_k$@>+def_flag|,
|c1+cite_flag|, \dots, |@t$c_l$@>+cite_flag|, $u_1$, \dots,~$u_m$,~0.

There can be multiple applied or cited occurrences of some module name
within one section, but only one defining occurrence. Therefore, in the
former cases we perform a test to avoid duplicate references.

@d file_flag (cite_flag-1)
		/* a distinguished value, not divisible by |cite_flag| */

@c
void new_mod_xref (mod_pointer p)
{ sixteen_bits head, *q, m=section_count+mod_xref_switch;
  if (p->xref->num==file_flag) q=&p->xref->next; /* skip |file_flag| */
  else head=xref_index(p->xref),q=&head;
  if (mod_xref_switch!=def_flag)
  { while (xnum(*q)>m)
      q=&xlink(*q); /* skip the $d_i$'s and possibly $c_i$'s */
    if (xnum(*q)==m) return; /* don't duplicate */
  }
  make_xref(m,*q); mod_xref_switch=0;
  if (q==&head) p->xref=xref_ptr; @+ else *q=xref_index(xref_ptr);
}

@ When a module name starts with \:(, we will call |set_file_flag|.

@c
void set_file_flag (mod_pointer p)
{@; if (p->xref->num!=file_flag)
    {@; make_xref(file_flag,xref_index(p->xref)); p->xref=xref_ptr; }
}

@ A third large area of memory is used for sixteen-bit `tokens', which
appear in short lists similar to the strings of characters in |byte_mem|.
Token lists are used to contain the result of \Cee\ code translated into
\TeX\ form; further details about them will be explained later. Sequences of
tokens which have been delimited as such are called texts; in fact the extent
of a text is implicitly defined by the start of the next text. To this end an
array~|text_mem| of pointers to the starting tokens of texts, in incresing
order, is maintained (a final pointer to beyond the last text is also
inclduded). To access both the pointers to the beginning and the end of a
text, one needs a pointer into~|text_mem|, which will have |text_pointer| as
type.

@<Typedef...@>=
typedef sixteen_bits token, * token_pointer, ** text_pointer;

@~The first position of |tok_mem| that is unoccupied by replacement text is
called |tok_ptr|, and the first unused location of |text_mem| is called
|text_ptr|. Since we already know that the next text to be stored will start
at |tok_ptr|, we make sure that |*text_ptr==tok_ptr| whenever we are not in
the process of appending tokens. In this way we can also always find the end
of any stored text as the beginning of the next one, which is computed by
the macro |text_end|.

@d tok_mem_end  (&tok_mem[max_toks]) /* end of |tok_mem| */
@d text_mem_end  (&text_mem[max_texts]) /* end of |text_mem| */
@d text_index(p) ((sixteen_bits)((p)-text_mem))
                 /* index of a |text_pointer p@;| in |text_mem| */
@d text_at(i) (&text_mem[i]) /* inverse operation of |text_index| */
@d text_begin(p) (*(p))
              /* get pointer to start of text from a |text_pointer p@;| */
@d text_end(p) (*(p+1))
              /* get pointer to end of text from a |text_pointer p@;| */

@<Global...@>=
token tok_mem[max_toks]; /* tokens */
token_pointer text_mem[max_texts]; /* directory into |tok_mem| */
token_pointer tok_ptr = tok_mem; /* first unused position in |tok_mem| */
text_pointer text_ptr = text_mem; /* first unused position in |text_mem| */
#ifdef STAT
token_pointer max_tok_ptr = tok_mem; /* largest value of |tok_ptr| */
text_pointer max_text_ptr = text_mem; /* largest value of |text_ptr| */
#endif

@ We initialise our invariant.

@<Set init...@>=
*text_ptr=tok_ptr;

@ Here are the three functions needed to complete |id_lookup|. For a
name~|x| stored in the table to match a string of length~|l| starting at~|q|
and of requested ilk~|ilk|, it is required that the names match exactly, and
either |x->ilk==ilk|, or |ilk==normal| and |x->ilk| specifies some reserved
word (where the ilks |TeX_like| and |NULL_like| are considered as reserved).
This rule means that if we look up an ``identifier'' whose name matches
reserved word stored in the table, the reserved word rather than an
identifier is returned, so that recognition of reserved words is automatic
when the tables are properly initialised. The other two functions install
the sentinel node |xmem[0]| at the end of each \xr. list as it is created.

@c
boolean names_match (id_pointer x, char* q, int l, int ilk)
{ char* p=name_begin(x);
  if ((x->ilk==ilk || ilk==normal && reserved(x)))
  @/{@; while (--l>=0) if (*p++!=*q++) return false; return *p=='\0'; }
  else return false;
}

void init_id_name (id_pointer p, int t)
{@; p->ilk = t; p->xref = &xmem[0]; }

void init_module_name (mod_pointer p) {@; p->xref=&xmem[0]; }


@*1 Skipping and copying \TeX\ material. @:title@>
The source file consists roughly speaking of two sorts of input, namely
\TeX\ material and pieces of \Cee~text. The limbo part is entirely the
realm of \TeX, and the sections are divided (possibly quite unevenly)
into a \TeX~part, and the remainder which is~\Cee, although the distinction
is not quite as pure here due to `\pb' interjections and comments.
The meddling of \.{\me.} (and indeed of all of \.{CWEB}) in the \TeX\
parts is extremely superficial: it is restricted to skipping or copying it,
only replacing \:@@ by~`\.{@@}', and the main interest in these parts
of the source is to find out where they end. By contrast the \Cee\ fragments
have to be broken up into meaningful parts by \.{\me.} in order to perform
the required operations of collecting \xr.s and formatting for
pretty-printing output. Before we study that more complicated process of
lexical scanning, let us study the easier question of what happens with the
\TeX~parts.

The three tasks of passing over limbo material, over the \TeX~part of a
section, and over a comment embedded in the \Cee~part are sufficiently
different that they merit separate routines. Moreover, the entire source is
scanned twice by \.{\me.}, with different purposes, and so the scanning
routines have two flavours. We could choose either to use distinct scanning
routines on both passes, or to have a single all-purpose scanning routine.
The first option is more attractive if the functions are simple and small,
and can perform tasks specific to their pass on-the-fly, while the second
option is to be favoured if the scanning task becomes complicated, to avoid
near-duplication of substantial pieces of code. We have chosen the first
option for the routines that pass over ordinary \TeX\ material, both in
limbo and in sections, and the second option for the function passing over
comments; for the lexical scanning of \Cee~text we shall also write a single
function for both phases. Where distinct routines are used we must of
course make sure that they are sufficiently similar that they will make the
same decisions about where these parts of the source text end.

@ Although the functions to be described here have a quite trivial task,
they will eventually run into a control code that terminates their action,
and so we have to be aware of tokens that do not belong in the \TeX~part.
The relevant tokens form only a small subset of the tokens that can be
recognised by the function |get_next| that scans \Cee~text, so we defer
an enumeration of possible tokens to a later point. At this point it
is sufficient to know there is an array |ccode| which translates characters
that may follow `\.{@@}' in a control code into a numeric value greater
than any (unsigned) character, which values have symbolic names and are
arranged in order of increasing significance. In particular any value
greater than or equal to |format| terminates the \TeX~part of a section,
and the largest value of all is |new_section|. By using the same encoding,
the \TeX~text scanning functions can return a value directly usable for
continuing the scan into \Cee~territory. Tokens consisting of a single
character are represented as that character itself by |get_next|, and at
this point the relevant case is the character~`\.\v' which interrupts
\TeX~text.

Like in |CTANGLE|, we always access |ccode| by means of the macro |code_of|.

@d code_of(c) ccode[(unsigned char)(c)]

@<Global...@>=
int ccode[UCHAR_MAX + 1]; /* meaning of a character following `\.{@@}' */

@ This section performs the simplest of all scanning operations, namely to
skip through portions of the input that are not in any sections, i.e., that
precede the first section, on the first pass. It uses the fact that
|get_line| places a~|' '| at~|*limit|, so that placing the sentinel cannot
inadvertently create a token \.{@@@@}. Although a few control codes like \:q
are allowed in limbo, they can be ignored here, since they do not affect the
recognition of the end of the limbo part. An exception is the format code
\:s, which must be obeyed in this first phase; the code doing this will be
given together with the processing of other format codes.  After the code
below is executed, the value of |input_has_ended| will tell whether or not a
section has actually been found.

@< Skip the limbo part @>=
while (find_char())
{ limit[1]='@@'; /* place a sentinel */
  while (*loc++!='@@') {}
  if (loc<=limit) /* note that |loc!=limit+1| since |*limit==' '| */
  { int c=code_of(*loc++);
    if (c==new_section) break;
    if (c==format) @< Process a format code in limbo @>
  }
}

@ In Phase~II, the corresponding task is slightly less trivial, as the limbo
material must also be copied to the output; it is performed by the
|copy_limbo| function. Output is generated by calling |out| for ordinary
characters, and |finish_line| when a completed line is to be sent out.  No
spaces or tab marks are copied by |copy_limbo| into the beginning of a line
(indicated by the macro |output_line_empty|), nor by the function |copy_TeX|
below. As a consequence a line with only spaces and removed items (like \:q,
\:s, or for |copy_TeX| \xr. entries for the index) will produce a completely
empty line of output; such lines will only be actually written out by
|finish_line| if they come from a completely blank input line. Any pair \:@@
is replaced by `\.{@@}', and apart from this only the control codes \:q,
\:s, and \:l are allowed to be present before the first section is
encountered. Note that we have chosen to detect violations of this rule in
Phase~II; on other occasions however, errors were checked in Phase~I, and we
can then assume that errors are absent in Phase~II, as we would not even get
there otherwise.

@c
void copy_limbo (void) /* copy \TeX\ code until the next section begins */
{ while (loc<=limit || (finish_line(),get_line()))
  { eight_bits c;
    limit[1]='@@'; /* place a sentinel */
    while ((c=*loc++)!='@@')
      if (!(output_line_empty() && isspace(c))) out(c);
    if (loc<=limit) /* then we have hit a control code */
      switch(code_of(*loc++))
      {
      case new_section: return; /* the only exit, unless no sections exist */
      case ignored_text: get_control_text(); break;
      case format: get_next(); get_next(); break; /* skip two identifiers */
      case char_trans: out_str("\\ATL "); break; @.\\ATL@>
      default: err_print("! Double @@ required in limbo part");
			  @.Double @@ required...@>
        /* fall through */
      case at_sign_image: out('@@');
      }
  }
}

@ The function |skip_TeX| is used on the first pass to skip through the
\TeX\ code at the beginning of a section. It returns the next control code
or `\.\v' found in the input. A |new_section| is assumed to exist at the
very end of the file. Any comment character `\.{\%}' that is not escaped
(with a backslash) will disable recognition of `\.\v' and cross-referencing
control codes for the remainder of the line, so that no spurious \xr.s
will be created to commented-out text. Recognition of control codes that
terminate the \TeX~part of the section will still be enabled however, to
maintain synchronisation with the activities of |CTANGLE|.

@f skip_TeX TeX

@c
int skip_TeX (void) /* skip past pure \TeX\ code */
{ char c;
  while (find_char())
  { limit[1]='@@';
    while ((c=*loc++)!='@@' && c!='%')
      if (c=='|') return c;
      else if (c=='\\' && *loc!='@@') ++loc;
	/* ignore `\.{\\\%}'  and `\.{\\\v}' */
    if (loc<=limit)
      if (c=='@@') return code_of(*loc++);
      else /* ignore remainder of line unless a major control code occurs */
	do @+
	  if ((c=*loc++)=='@@' && code_of(*loc++)>=format)
	    return code_of(loc[-1]);
	while (loc<limit);
  }
  return new_section;
}

@ During Phase~II, the function |copy_TeX| processes the \TeX\ code at the
beginning of a section; for example, the words you are now reading were
copied in this way. Like |skip_TeX|, it returns the next control code or
`\.\v' found in the input, and ignores anything that is commented out,
except control codes that terminate the \TeX~part; the characters after the
`\.{\%}' are not even copied to the output. (Note that in the limbo part
comments {\sl are\/} copied, allowing any commented-out header information
such as copyright notices to be passed on to the \TeX~file.)

@f copy_TeX TeX

@c
int copy_TeX (void) /* copy pure \TeX\ material */
{ eight_bits c; /* current character being copied */
  while (loc<=limit || (finish_line(),get_line()))
  { limit[1]='@@';
    while((c=*loc++)!='@@')
    { if (c=='|') return '|';
      if (!(output_line_empty() && isspace(c))) out(c);
      if (c=='%') break;
      if (c=='\\' && *loc!='@@') out(*loc++);
	/* copy `\.{\\\%}' and `\.{\\\v}' */
    }
    if (loc<=limit)
      if (c=='@@') return code_of(*loc++);
      else /* ignore remainder of line unless a major control code occurs */
	do
	  if ((c=*loc++)=='@@' && code_of(*loc++)>=format)
	    return finish_line(),code_of(loc[-1]);
	while(loc<limit);
  }
  return new_section;
}

@ The final \TeX~scanning function is |scan_comment| which scans the
\TeX~text in comments. It is used both in Phase~I and in Phase~II, but
stores the characters it sees only if |phase==2|. The function returns the
token that terminated its action, which is one of |end_comment|, `\.\v' or
(in erroneous cases) |new_module|.

@< Prototypes @>=
int scan_comment(int* bal, boolean one_liner);
   /* skip or copy \TeX\ text in comments */

@~The function |scan_comment| counts the braces it encounters to see if they
are balanced; the parameter |bal| points to an integer variable keeping
track of the brace level (it can be positive initially if |scan_comment| is
called after `\pb' occurring within braces in a comment).  This feature is
mainly a remnant from \.{CTANGLE}'s Pascal origins, since there comments are
closed by `\.\}', and counting is necessary to establish which braces are
for \TeX\ and which are not. Although there is no such need in \Cee, the
feature is retained. The parameter |one_liner| tells whether we are dealing
with a \Cpp\ one-line comment.

@c
int scan_comment (int* bal, boolean one_liner)
{ char c; boolean forced_out=false; /* prematurely terminated? */
  while (one_liner ? loc<limit
		   : find_char() && (*loc!='*' || loc[1]!='/' ))
  @< Handle next character; if |'|'|, |return| it,
     if a new section starts, set |forced_out| and |goto done| @>
  if (input_has_ended)
    forced_out=true,err_print("! Input ended in mid-comment");
			       @.Input ended in mid-comment@>
  else if (!one_liner) loc+=2; /* move past `\.{*{}/}' */
done:
  if (*bal>0) err_print("! Too few closing braces in comment");
			 @.Too few closing braces...@>
  return forced_out ? new_section : end_comment;
}

@ Like elsewhere, `\.@@' should be doubled in comments; |scan_comment|
replaces them by a single `\.@@'. In Phase~II, instead of copying the \TeX\
material into the output buffer like |copy_TeX|, |scan_comment| copies it
into the token memory, since comments will be output together with the rest
of the formatted \Cee~code. To this end it calls the macro |app_char_tok(c)|
rather than |out(c)|. When including header files comments should be skipped
with no processing at all, so this module is almost completely disabled in
that case.

@< Handle next character... @>=
if (including_header_file) ++loc; /* don't process characters here */
else
{ switch(c=*loc++)
  {
  case '|': return '|'; /* beginning of `\pb' inside comment */
  case '@@':
    if (*loc++!='@@')
      if (code_of(loc[-1])!=new_section)
	err_print("! Double @@ required in comment");
		   @.Double @@ required...@>
      else
      {@; err_print("! Section ended in mid-comment");
		     @.Section ended in mid-comment@>
	forced_out=true; goto done;
      }
    break;
  case '\\':
    if (*loc!='@@') {@; if (phase==2) app_char_tok(c); c=*loc++; }
    break;
  case '{': ++*bal; break;
  case '}':
    @+if (*bal>0) --*bal;
    @+else err_print("! Extra } in comment");
		     @.Extra \} in comment@>
    @+break;
  case '/': @+ if (*loc=='*') err_print("! Nested comment");
					 @.Nested comment@>

  }
  if (phase==2) app_char_tok(c);
}


@*1 Getting the next token. @:title@>
We now come to the most important lexical scanning function, |get_next|,
which locates and classifies the next token of \Cee~text. Before we give the
function itself, we shall first specify which kind of values it can return.
The result value is an integer, which can either be at most |UCHAR_MAX|, in
which case the character code represents that character, or possibly a
compressed multi-character symbol as defined by the macros
in~\.{common.inc}, or it is one of the values greater than |UCHAR_MAX|
defined in the enumeration below, which indicates a particular token or
class of tokens or a control code. (Actually, |underline| and
|trace0|,~\dots,~|trace3| are never returned by |get_next| because they are
treated within the scanner; they are included however because they are
distinguished control codes.)  The ordering of this enumeration is designed
to simplify \.{\me.}'s logic; for example, from |format| on, larger
numbers are given to the control codes that denote more significant
milestones, and the code of |new_section| is the largest of all.  Also, the
sequence |identifier|,~\dots,~|xref_mark| runs parallel to the first five
|ilk|~codes. For efficiency all codes that produce a fixed scrap are placed
at the beginning; these are all codes preceding |ignore|.  The code \:> is
treated as ignored because it should not occur in places where we are not
searching for it.

@f TeX_string TeX

@<Typedef and enum...@>=
enum @/
{ at_sign_image = UCHAR_MAX+1, /* quoted `\.{@@}' */
  or, /* \:v */
  mul_assign, div_assign, mod_assign, plus_assign, minus_assign,
  left_assign, right_assign, and_assign, xor_assign, or_assign,
  sh_sh, ellipsis, colon_colon, @/
  start_preproc, end_preproc, /* begin and end of a preprocessor directive */
  join, /* \:\& */
  thin_space, /* \:, */
  math_break, /* \:\v */
  line_break, /* \:/ */
  big_line_break, /* \:) */
  no_line_break, /* \:+ */
  backup_line, /* \:\\ */
  pseudo_semi, /* \:; */
  force_expr_open, force_expr_close, /* \:[, \ \:] */
  include_preproc, /* \:p */
@)
  ignore, /* control code of no interest to \.{\me.} */
  constant, string, /* the next five codes should remain in this order */
  identifier, /* any (possibly reserved) word found in \Cee\ text */
  xref_roman, xref_wildcard, xref_typewriter, xref_mark,
    /* \:\^, \ \:?, \ \:., \ \:: */
  refer, /* \:\# */
  TeX_string, /* \:t */
  verbatim, /* \:= */
  ignored_text, /* \:q */
  char_trans, /* \:l */
  ASCII_code, /* \:' */
  begin_comment, end_comment,
  underline, /* \:! */
#ifdef DEBUG
  trace0, trace1, trace2, trace3, /* \:0, \dots, \:3 */
#endif
  format, /* \:f */
  definition, /* \:d */
  header, /* \:h */
  begin_C, /* \:c */
  module_name, /* \:< and \:( */
  new_section /* \:\ , \:\~ and \:* */
};

@ Here we initialise the |ccode| table in accordance with the comments given
in the enumeration above.

@<Set ini...@>=
{ unsigned char c=0;
  do ccode[c] = isspace(c) ? new_section : ignore; while(c++!=UCHAR_MAX);
  ccode['@@'] = at_sign_image;
@/ccode['v'] = ccode['V'] = or;
@/ccode['!'] = underline; /* set definition flag */
@/ccode['^'] = xref_roman; /* index entry to be typeset normally */
@/ccode['?'] = xref_wildcard; /* index entry to be in user format */
@/ccode['.'] = xref_typewriter; /* index entry to be in typewriter type */
@/ccode[':'] = xref_mark; ccode['#']=refer; /* explicit \xr.s */
@/ccode['t'] = ccode['T'] = TeX_string; /* \TeX\ box within \Cee\ text */
@/ccode['='] = verbatim;
@/ccode['q'] = ccode['Q'] = ignored_text;
@/ccode['l'] = ccode['L'] = char_trans;
@/ccode['\''] = ASCII_code;
@/ccode['&'] = join; /* concatenate two tokens */
@/ccode[','] = thin_space;
@/ccode['|'] = math_break;
@/ccode['/'] = line_break;
@/ccode[')'] = big_line_break;
@/ccode['\\']= backup_line;
@/ccode['+'] = no_line_break;
@/ccode[';'] = pseudo_semi;
@/ccode['['] = force_expr_open; ccode[']'] = force_expr_close;
ccode['p'] = ccode['P'] = include_preproc;
#ifdef DEBUG
  ccode['0'] = trace0; ccode['1'] = trace1;
  ccode['2'] = trace2; ccode['3'] = trace3;
#endif
@/ccode['f'] = ccode['F'] = ccode['s'] = ccode['S'] = format;
@/ccode['d'] = ccode['D'] = definition;
@/ccode['h'] = ccode['H'] = header;
@/ccode['c'] = ccode['C'] = begin_C;
      /* \Cee\ text in unnamed module */
@/ccode['<'] = ccode['('] = module_name; /* beginning of a module name */
@/ccode['~'] = ccode['*'] = new_section; /* beginning of a new section */
  if (compatibility_mode)
    @< Reset some control codes to match \LKC. @>
}

@ In \.{CWEBx} there are a few control codes that also exist in Levy/Knuth
\.{CWEB} but have a different meaning. In compatibility mode we reassign the
meaning of these codes to that of \LKC., making their usual function
inaccessible, since it is not intended that hybrid programs should be
written using the codes of \LKC. together with features particular to
\.{CWEBx}.
@^Levy/Knuth \.{CWEB}@>

@< Reset some control codes... @>=
{ ccode['h']=ccode['H']=include_preproc; /* \:h means \:p */
  ccode['p']=ccode['P']=begin_C; /* \:p means \:c */
  ccode['#']=big_line_break; /* \:\# means \:) */
  ccode[':']=xref_wildcard; /* \:: means \:? */
}

@ We come now to the definition of |get_next| itself. When returning certain
values it will have performed some additional actions, as follows.

\yskip\hang |constant|, |string|, |TeX_string|, |verbatim|: The token is
    copied into |mod_text|, with slight modifications; the global variables
    |id_first| and |id_loc| are set to the beginning and ending-plus-one
    locations in |mod_text|.

\yskip\hang |identifier|, |xref_roman|, |xref_wildcard|, |xref_typewriter|,
    |xref_mark|, |refer|, |module_name|: The global variable |cur_id| or
    |cur_mod| will point to the identifier, control text or module name that
    has just been scanned (for the |xref_roman|, \dots, |xref_mark| this is
    only true if |phase==1|).

\yskip\hang |underline|: this value is not even returned. If |get_next| sees
    \:!, it sets |xref_switch| to |def_flag| and goes on to the next token.
\yskip

Preprocessing directives complicate scanning in two ways: first, their
lexical structure is different from that of ordinary text, and second, a
preprocessor directive can occur at any place in a \Cee~text, so
syntactically it should be treated like a comment. The first issue is
resolved by maintaining a static variable |preprocessing| which is~$0$ in
ordinary \Cee~text, $2$ in \&{\#include} directives, and~$1$ in other
preprocessor directives. The second issue must be dealt with during parsing,
and in order to be able to do so, |get_next| emits special tokens
|start_preproc| and |end_preproc| when it has sensed the boundaries of a
preprocessor directive.

@<Global...@>=
id_pointer cur_id; /* identifier or index entry just scanned */
mod_pointer cur_mod; /* module name just scanned */
int preprocessing=0;

@ As one might expect, |get_next| consists mostly of a big switch that
branches to the various cases that can arise. Any character |c>=0x80| that
is not contained in a string or comment is assumed to belong to an
identifier.

@< Prototypes @>= int get_next (void);
@~@c
int get_next (void) /* produces the next input token */
{ eight_bits c; /* the current character */
restart:
  if (!find_char()) {@; preprocessing=0; return new_section; }
  @< If a preprocessor line has ended, handle it and |return end_preproc| @>
  if ((c=*loc++)=='@@')
    @< Get control code and possibly module name, and either |return| it,
       or |goto restart| if ignored or handled within |get_next| @>
  if (isspace(c))
    if (preprocessing>0) return ' '; /* keep spaces in preprocessor lines */
    else goto restart; /* ignore other white space */
  if (c=='L' && (*loc=='\'' || *loc=='"'))
  {@; get_string(); return string; }
  if (isalpha(c) || c=='_' || c>=0x80)
    {@; @< Get an identifier @> return identifier; }
  if (isdigit(c) || c=='.' && isdigit((eight_bits)*loc))
    @/{@; @< Get a numeric constant @> return constant; }
  if (c=='\'' || c=='"' || (c=='<' && preprocessing==2))
    {@; get_string(); return string; }
  if (c=='#' && loc==&buffer[1])
  @/{@; @< Handle start of preprocessor directive; maybe |goto restart| @>
     return start_preproc; }
  if (c=='\\' && preprocessing>0 && loc==limit)
    { ++loc; /* move past |limit|, so |get_line| will be called */
      goto restart;
    }
  @< Compress multi-character tokens @>
  return c;
}

@ We will start scanning a header file when \:h is encountered in the first
phase, and also sometimes when a `\.{\#include}' directive is seen. To avoid
code duplication, we define a function to be called in both cases.

@< Prototypes @>=
boolean push_header_file(boolean suspend);

@~When a `\.\#' is seen as the first character of a line, |get_next| returns
a special code |start_preproc| and sets |preprocessing| to a non-zero value.
Because of the freakish use of `\.<' and `\.>' to delimit a file name in
lines that start with `\.{\#include}', those lines get an extra-special
treatment, and |preprocessing| is set to~$2$ rather than to~$1$. If however
we encounter a `\.{\#include}' directive when we are already busy reading a
header file due to \:h in Phase~I, then we actually execute the directive,
and start reading the nested header file; this action will be transparent to
the stream of tokens produced, so we |goto restart| to fetch the first token
from the newly opened file.

@<Handle start of prep...@>=
{ while (loc<limit && isspace((eight_bits)*loc)) ++loc;
    /* allow spaces after `\.\#' */
  if (limit-loc>=7 && strncmp(loc,"include",7)==0) /* `\.{\#include}' line */
    if (including_header_file) /* start nested header file */
    @/{@; loc+=7; push_header_file(false); goto restart; }
    else preprocessing=2;
  else preprocessing=1;
}

@ When about to open a header file, we check whether it has already been
scanned before. This simulates the often used mechanism (based on preprocessor
symbols) for avoiding multiple inclusions; even if the compiler should
actually see the same file multiple times, there is no point for us to scan in
multiple times. In any case we must avoid inclusions loops  (and we cannot use
the preprocessor), and even in the absence of those, multiple inclusion of the
same file would give an important performance penalty. We use the |xref| field
of an entry specifically created for the header file; it will point to |xmem[0]|
upon first encounter of the name, and is set to a different value below once
the name has been seen.

@c
boolean push_header_file(boolean suspend)
{ id_pointer p;
  if (locate_file_name() &&
      (p=id_lookup(id_first,id_loc,header_file_name))->xref==&xmem[0])
  { p->xref=&xmem[1]; /* mark file as seen, to avoid multiple inclusion */
    return push_input_file(true,suspend);
  }
  return false;
}

@ When we get to the end of a preprocessor line, we lower the flag and
send a code |end_preproc| (unless the last character was a `\.\\', but that
case has already been taken out and never comes here).

@<If a prep...@>=
if (preprocessing>0 && loc==limit)
{@; preprocessing=0; return end_preproc; }

@ The following code assigns values to the compound operators `\.{++}',
`\.{--}', `\.{->}', `\.{>=}', `\.{<=}', `\.{==}', `\.{<<}', `\.{>>}',
`\.{!=}', `\.{\v\v}', and `\.{\&\&}', to the special symbols `\.{/*}',
`\.{*/}', `\.{//}', `\.{...}', `\.{::}', and `\.{\#\#}', and moreover, if not
in compatibility mode, to the assignment operators `\.{*=}', `\.{/=}',
`\.{+=}', `\.{-=}', `\.{>>=}', `\.{<<=}', `\.{\&=}', `\.{\^=}', `\.{\v=}'.
Although the comment ending token `\.{*/}' should never occur when we are
scanning \Cee~text, we must recognise it in order to detect unclosed `\pb'
constructions within comments (fortunately no legal combination of operators
causes an adjacent sequence `\.{*/}'; only `\.*' immediately followed by a
comment could cause it, and the user should simply not do this). The reason
that in compatibility mode we are forced to follow suit with \LKC.
@^Levy/Knuth \.{CWEB}@> in not combining assignment operators is a truly
stupid one: @:truly stupid@> the Stanford GraphBase @^Stanford GraphBase@>
contains fragments like `\.{\$\v n1\v=n\_1\$}' which will cause trouble if
`\.{\v=}' is parsed as a single symbol (of course the fragments should have
been written as `\.{\v n1\v\$\{\}=n\_1\$}', but we have set the goal to handle
the GraphBase as it is, and moreover its files cannot be changed). For
splitting up the assignment operators here, we shall have to pay the price of
including syntax rules (in compatibility mode) for recombining them.

The macros defined below are for strictly local use, so we don't mind that
they could mess things up if used in the |if|-part of an |if|-|else|
statement. Also, there is no need to test whether |*loc| or |loc[1]| lie
beyond |limit|, since if they do, they are not evaluated anyway because we
know |*limit==' '|. Note that the three-symbol assignment operators
`\.{>>=}' and~`\.{<<=}' must be tested before their non-assignment
counterparts `\.{>>}' and~`\.{<<}'.

@d compress2(char2,code) if (*loc==char2) return ++loc, code @;
@d compress3(char2,char3,code)
  if (*loc==char2 && loc[1]==char3) return loc+=2, code @;
@d comp_ass_op2(code)
  if (*loc=='=' && !compatibility_mode) return ++loc, code @;
@d comp_ass_op3(char2,code)
  if (*loc==char2 && loc[1]=='=' && !compatibility_mode) return loc+=2,code @;

@<Compress multi...@>=
switch (c) {
case '/': compress2('*',begin_comment); @+
  if (C_plus_plus) compress2('/',begin_comment);
  comp_ass_op2(div_assign); break;
case '*': compress2('/',end_comment);	comp_ass_op2(mul_assign);  break;
case '%': comp_ass_op2(mod_assign);	break;
case '+': compress2('+',plus_plus);	comp_ass_op2(plus_assign); break;
case '-': compress2('-',minus_minus);	compress2 ('>', minus_gt);
	  comp_ass_op2(minus_assign);	break;
case '=': compress2('=',eq_eq);		break;
case '>': compress2('=',gt_eq);		comp_ass_op3('>',right_assign);
	  compress2 ('>',gt_gt);	break;
case '<': compress2('=', lt_eq);	comp_ass_op3('<',left_assign);
	  compress2 ('<', lt_lt);	break;
case '&': compress2('&',and_and);	comp_ass_op2(and_assign);   break;
case '^': comp_ass_op2(xor_assign);	break;
case '|': compress2('|',or_or);		comp_ass_op2(or_assign);    break;
case '!': compress2('=',not_eq);	break;
case '.': compress3('.','.', ellipsis); break;
case '#': compress2 ('#', sh_sh);	break;
case ':': @+ if (C_plus_plus)		compress2 (':',colon_colon);
}

@ The code below is almost identical to the corresponding module in
|CTANGLE|; the difference is that we accept characters |c>=0x80| without
testing whether a translation is defined for them (since we have not
recorded such information) and that we always look up identifiers, even if
they have length~$1$.

@< Get an identifier @>=
{  id_first=--loc; /* mark beginning of identifier */
   do c=*++loc; while (isalnum(c) || c=='_' || c>=0x80);
   cur_id= id_lookup(id_first,loc,normal);
}

@ In \Cee~text, numeric constants are specified in the ordinary \Cee~manner:
octals start with `\.0', hexadecimals with `\.{0x}', and anything starting
with a non-zero digit is a decimal constant; however for octal and
hexadecimal constants \.{\me.} will produce output using italics or
typewriter font, respectively, and introduced by a raised circle or hash
mark. Forgivably in contradiction with the definition of~\Cee, we treat the
ubiquitous constant~`0' as decimal rather than as octal. When the kind of
constant we are dealing with has been recognised, we represent this
information internally by special marker characters, which replace the marks
used in the \Cee~source (like `\.{0x}' for hexadecimal or `\.E' for the
exponent of a floating point constant). These markers are characters like
`\.\^' that will get a backslash prepended by the output routine for
constants, so that formatting of the constants can be controlled by defining
the corresponding \TeX\ control words (like `\.{\\\^}' in the case
mentioned) appropriately.

@d shift_and_store(ch) (*id_loc++=ch,c=*++loc)

@<Get a numeric constant@>=
{ id_first=id_loc=&mod_text[1];

  if (c=='0' && (isdigit(c=*loc) || tolower(c)=='x')) /* octal or hex */
  { if (isdigit(c)) /* octal constant with at least two digits */
    { *id_loc++ = '~'; /* store `\.\~' in place of leading `\.0' */
      do shift_and_store(c); while (isdigit(c));
	/* copy second and following digits */
    }
    else /* hex constant */
    { shift_and_store('^'); /* replace `\.{0x}' by `\.\^' */
      while (isxdigit(c)) shift_and_store(c);
    }
  }
  else /* decimal constant */
  { c=*--loc; /* recover first digit or decimal point */
    while (isdigit(c)) shift_and_store(c);
    if (c=='.') @+ do shift_and_store(c); while (isdigit(c));
    if (tolower(c)== 'e') /* floating point constant with exponent */
    { shift_and_store('_'); /* replace `\.e' by `\.\_' */
      if (c=='+' || c=='-') {@; *id_loc++ = c; c=*++loc; }
      while (isdigit(c)) shift_and_store(c); /* exponent */
    }
  }
  if (isalpha(c)) /* `\.{U}', `\.{L}', and/or `\.{F}' suffix */
  @/{@; *id_loc++ = '$'; @q $ emacs-cookie @>
    do shift_and_store(c); while (isalpha(c)); }
}

@ After an `\.{@@}' sign has been scanned, the next character tells us
whether there is more work to do.  This code uses the fact that our internal
code numbers |xref_roman|, |xref_wildcard|, |xref_typewriter|, and
|xref_mark| are consecutive in the stated order, as are the |ilk| codes
|roman|, |wildcard|, |typewriter|, and |reference|. We silently eliminate the
possibility of indexing the empty string, since it would cause anomalous
situations in hashing and sorting of the index, and it would look rather
silly anyway.

@<Get control code and possibly module name...@>=
if (including_header_file) goto restart; /* ignore `\.@@' in header files */
else
{ int cc=code_of(*loc++);
  switch (cc)
  { case ignore: goto restart;
    case underline: xref_switch=def_flag; goto restart;
#ifdef DEBUG
    case trace0: case trace1: case trace2: case trace3: @+
      if (phase==2) tracing=cc; @+ goto restart;
#endif
    case char_trans:
      err_print("! `@@l' only allowed in limbo"); goto restart;
		 @.`@@l' only allowed in limbo@>
    case ASCII_code: @< Scan an \caps{ASCII} constant @> @+ return string;
    case module_name:
      @< Scan the module name and make |cur_mod| point to it @> @+ break;
    case ignored_text: get_control_text(); goto restart;
    case verbatim: case TeX_string: get_control_text(); break;
    case xref_roman: case xref_wildcard: case xref_typewriter:
    case xref_mark: case refer:
      if (get_control_text()) goto restart; /* don't index empty strings */
      if (cc==refer) cur_id=id_lookup(id_first,id_loc,reference);
      else if (phase==1)
	cur_id=id_lookup(id_first,id_loc,cc-xref_roman+roman);
  }
  return cc;
}

@ There is no reason why we should allow a newline within an \caps{ASCII}
constant, even if it is escaped.

@< Scan an \caps{ASCII} constant @>=
{ id_first=&mod_text[1]; strncpy(id_first,"@@'",2); id_loc=&id_first[2];
  while ((*id_loc++=c=*loc++)!='\'')
  { if (c=='\\')
      *id_loc++=*loc++; /* copy any character following backslash */
    else if (c=='@@' && *loc++!='@@')
    {@; err_print("! Double @@ required in strings"); --loc; }
		   @.Double @@ required...@>
    if (loc>=limit) {@; err_print("! ASCII constant didn't end"); break; }
				   @.ASCII constant didn't end@>
  }
}

@ Here |get_module_name| does nearly all the work; we only need to recognise
when \:( is used rather than~\:<, and if so insert |file_flag| in the
appropriate \xr. list.

@< Scan the module name... @>=
{ boolean file_module=loc[-1]=='(';
  cur_mod=get_module_name();
  if (file_module && phase==1 && cur_mod!=NULL) set_file_flag(cur_mod);
}


@* Phase I processing. @:title@>
We now have accumulated enough functions to make it possible to carry out
\.{\me.}'s first pass over the source file.  If everything works right,
both Phase~I and Phase~II of \.{\me.} will assign the same numbers to
sections, and these numbers will agree with what \.{CTANGLE} does.

We keep track of the current section number in |section_count|, which
is the total number of sections that have started.  Sections which have
been altered by a change file entry have their |changed_section| flag
turned on during the first phase. Meanwhile we also keep track using
|change_exists| of whether any change was made at all, which will tell us
whether the index has changed. The global variable |next_control| often
contains the most recent output of |get_next|; in interesting cases, this
will be the control code that ended a section or part of a section.

@d shift() (next_control=get_next())

@<Global...@>=
boolean change_exists=false; /* has any section changed? */
int next_control; /* control code waiting to be acted upon */

@ The overall processing strategy in Phase~I has the following
straightforward outline.

@c
void phase_one (void)
  /* read all the user's text and store the \xr.s */
{ phase=1; reset_input(); section_count=0;
  @< Skip the limbo part @>
  while (!input_has_ended)
    @< Store cross-reference data for the current section @>
  if (change_exists) mark_section_as_changed(section_count);
    /* the index changes if anything does */
  @< Print error messages about unused or undefined module names @>
  @< Reverse the \xr. lists for identifiers @>
}

@ The outline for each section is equally straightforward.

@< Store cross-reference data... @>=
{ if (++section_count==max_sections)
    overflow("section number"); @.section number capacity exceeded@>
  if (loc[-1]=='*') print_section_progress ();
  @< Store cross-references in the \TeX~part of a section @>
  @< Store cross-references in the definition part of a section @>
  @< Store cross-references in the \Cee~part of a section @>
  if (section_changed(section_count)) change_exists=true;
 }

@ We interrupt our refinement of |phase_one| temporarily for some auxiliary
functions that are used in its various parts.

@< Prototypes @>=
void C_xref (boolean); /* make \xr.s within in straight \Cee~text */
void outer_xref (void); /* make \xr.s in \Cee~text with comments */
void mod_check (mod_pointer); /* check \xr.s for module names */

@ The function |C_xref| stores references to identifiers in \Cee~text,
either for an entire fragment enclosed in `\pb', or for a portion of a macro
or section body delimited by comments or module names. The boolean parameter
|inner| tells whether the former is the case; if so |C_xref| should stop
with |next_control=='|'|, and otherwise it should stop with either
|next_control==begin_comment| or |next_control>=format|. In fact, setting
|inner| will make |C_xref| stop at |'|'| but proceed past module names,
while comment delimiters and major control codes will make |C_xref| stop
regardless of |inner|. If |next_control>=format| when |C_xref| is called,
nothing will happen, but it is safe to call the function when
|next_control=='|'| or |next_control==end_comment|, which will be stepped
over rather than considered as termination condition. Thus we can avoid
saying |shift()| immediately before calling |C_xref| on several occasions.
After a `\.\#' that starts a preprocessor directive, an identifier must
follow, which must not be \xr.d; this is achieved by performing an extra
|shift()|. If no identifier follows, we report an error, but do not perform
the extra |shift()|.  The code below uses the fact that our internal code
numbers |identifier|, |xref_roman|, |xref_wildcard|, |xref_typewriter|, and
|xref_mark| are consecutive.

The other task of Phase~I is to collect all unusual |ilk| assignments,
and the subtler part of this is the processing of |typedef| declarations.
The details of this process will be explained later, but this is where
it is hooked into the other actions, since all relevant tokens (which do
not include module names or comments) pass here one by one; we exempt tokens
from preprocessor lines from consideration.

@c
void C_xref (boolean inner)
{ while (next_control<format || next_control==module_name && inner)
  { if (preprocessing==0)
      @< Keep track of tokens relevant to |typedef| declarations @>
    if (next_control>=identifier && next_control<=xref_mark)
      new_id_xref(cur_id);
    else if (next_control==module_name && cur_mod!=NULL)
      mod_xref_switch=cite_flag,new_mod_xref(cur_mod);
    if (next_control==start_preproc && shift()!=end_preproc
      &&next_control!=identifier)
      err_print("! Identifier should follow `#'");
		 @.Identifier should follow `\#'@>
    else shift();
    if (next_control=='|' && inner
     || next_control==begin_comment || next_control==end_comment)
      return;
  }
}

@ The function |outer_xref| is like |C_xref|, but is used to scan an entire
macro body, or a portion of a section body delimited by module names; it
handles \Cee~text with embedded comments. It is called after \:d, a format
definition, \:c, or a module name (either defining or applied) has been
scanned, and after \:h has fired up its header file; in all cases
|next_control| is already processed, and we start with |shift()|. (There is
also one call that can occur during Phase~II, namely if illegal items
following a \:s format definition were found; its purpose then is merely to
ensure |next_control>=format| without producing any output.)

While a comment is being scanned, tokens that pass |C_xref| should not be
considered as part of a possible |typedef| that is in progress; this is
achieved by invoking the macro |typedef_tracking| with the proper boolean
value at the beginning and end of the comment.

@c
void outer_xref (void) /* extension of |C_xref| */
{ shift(); /* move past previously processed token */
  while (next_control<format)
    if (next_control!=begin_comment) C_xref(false);
    else
    { boolean one_liner=loc[-1]=='/'; int bal=0; /* brace level in comment */
      typedef_tracking(false);
      while ((next_control=scan_comment(&bal,one_liner))=='|')
	@/{@; C_xref(true); if (next_control!='|') break; }
      typedef_tracking(true);
    }
}

@ In the \TeX~part of a section, \xr. entries are made only for the
identifiers in \Cee\ texts enclosed in `\pb', or for control texts
introduced by \:\^, \:., \:?, or \::.

@< Store cross-references in the \TeX... @>=
do
  switch (next_control=skip_TeX())
  { case underline: xref_switch=def_flag; break;
    case '|': C_xref(true); break;
    case module_name: case refer: loc-=2; get_next(); break;
    case ignored_text: get_control_text(); break;
    case char_trans: err_print("! `@@l' only allowed in limbo"); break;
				@.`@@l' only allowed in limbo@>
    case xref_roman: case xref_wildcard: case xref_typewriter:
    case xref_mark: loc-=2; get_next(); new_id_xref(cur_id);
  }
while (next_control<format);

@ During the definition and \Cee~parts of a section, \xr.s are made for all
identifiers except reserved words and the right hand sides of format
definitions. The \TeX\ code in comments is, of course, ignored, except for
\Cee\ portions enclosed in `\pb'; the text of a module name is skipped
entirely, even if it contains `\pb' constructions.

@ When we get to the following code we have |next_control>=format|.

@< Store cross-references in the def... @>=
while (next_control<begin_C) /* |format|, |definition| or |header| */
  if(next_control!=header)
  { xref_switch=def_flag; /* implied \:! for first identifier */
    if (next_control==format) @< Process a format code in a section @>
    outer_xref(); /* macro definition or comment after format definition */
  }
  else @< Read a header file, scanning it for |typedef| declarations @>

@ Before we handle format codes occurring in a section, let us consider
their treatment in limbo. Here the code must be a non-printing \:s rather
than \:f, since we are not prepared to emit formatted output in limbo. The
syntax and semantics are simple: two identifiers must follow, and the |ilk|
of the latter is is assigned to the former.

@< Process a format code in limbo @>=
if (tolower((eight_bits)loc[-1])=='f')
  err_print("! Double @@ required in limbo part");
	     @.Double @@ required...@>
else
{ id_pointer lhs;
  if (shift()==identifier && (lhs=cur_id,shift()==identifier))
    lhs->ilk=cur_id->ilk;
  else err_print("! Improper format definition");
		  @.Improper format definition@>
}

@ In a section \:s is processed in the same way; for \:f we additionally
produce a defining \xr. for the left hand side. We do not call |shift| since
this will be done by |outer_xref|, which is called after this code to
process any comments that follow the format definition.

@< Process a format code in a section @>=
{ boolean f= tolower((eight_bits)loc[-1])=='f';
  id_pointer lhs;
  if (shift()==identifier && (lhs=cur_id,shift()==identifier))
  { if (f) new_id_xref(lhs); @+ else xref_switch=0;
    lhs->ilk=cur_id->ilk;
  }
  else err_print("! Improper format definition");
		  @.Improper format definition@>
}

@ After \:h a file name follows, enclosed in white space, double quotes, or
angle brackets. We open the header file by the same routine used to open \:i
files, but here we do suspend reading from the change file (via the |true|
argument to |push_header_file|). If all is well the whole file will contain
pure \Cee~text without any control codes, so |outer_xref| will come to a halt
shortly after returning from that file.

@< Read a header file... @>=
{ if (push_header_file(true)) /* prepare for reading header file */
    including_header_file=true; /* will be reset on closing the file */
  typedef_tracking(true); /* this is what we are doing it for */
  outer_xref();
    /* |shift()| and  collect typedefs until |next_control>=format| */
  typedef_tracking(false);
}


@ Finally, when the \TeX\ and definition parts have been treated, we have
|next_control>=begin_C|. The loop repeatedly marks a module name
\xr. and then calls |outer_xref| to scan everything up to the
next module name; if |next_control==module_name| initially, we raise
|mod_xref_switch| which will cause that first module name to be marked as
defining.

@< Store cross-references in the \Cee... @>=
{ if (next_control<new_section) /* |begin_C| or |module_name| */
  { typedef_tracking(true);
    mod_xref_switch= next_control==module_name ? def_flag : 0;
    do
    { if (next_control==module_name && cur_mod!=NULL)
	new_mod_xref(cur_mod);
      outer_xref();
    } while (next_control<new_section);
    typedef_tracking(false);
  }
}

@ After we have seen everything, we want to check that each module name was
both defined and used. The following recursive function walks through the
tree of module names alphabetically and prints out anomalies. It also
reverses the sublists of \xr.s to sections where a module is being defined,
cited, and used. In this section we use the fact that \xr. lists end with a
sentinel with |num==0|.
@^recursion@>

@c
void mod_check (mod_pointer p) /* print anomalies in subtree |p| */
{ if (p != NULL)
  { mod_check (p->llink); /* traverse left subtree */
    { boolean file_module = p->xref->num==file_flag;
      sixteen_bits head, *q, threshold;
        /* lower limit of |num| values of current interest */
      if (file_module) q=&p->xref->next; @+
      else head=xref_index(p->xref),q=&head;
      if (!complete_name(p))
      @/{@; print("\n! Never completed"); print_mod(p); mark_harmless(); }
		     @.Never completed: <module name>@>
      if (xnum(*q)<=(threshold=def_flag))
      @/{@; print("\n! Never defined"); print_mod(p); mark_harmless(); }
		     @.Never defined: <module name>@>
      else
      @< Reverse sublist after |*q| with entries |num>threshold|;
	 make |q| point to final |next| field @>
      if (xnum(*q)>(threshold=cite_flag)) @/@< Reverse sublist... @>
      if (xnum(*q)==(threshold=0))
      @/{@; if(!file_module)
	  {@; print("\n! Never used"); print_mod(p); mark_harmless(); }
		       @.Never used: <module name>@>
      }
      else @< Reverse sublist... @>
      if (!file_module) p->xref=&xmem[head];
	/* set pointer to possibly modified value */
    }
    mod_check (p->rlink); /* traverse right subtree */
  }
}

@~@<Print error messages about un...@>=
mod_check(root);

@ We now come to the reversal of the \xr. lists, which is necessary because
by repeatedly prepending elements to these lists, or in the case of module
names to one of three sublists, these (sub)lists have obtained reverse
ordering. The method of traversal of the set of all identifiers is different
from that for the set of all module names, whence these tasks are linked
into the program at different points, but the reversal routines themselves
are quite similar. As we already traverse the module names in |check_root|,
the reversal code for module \xr.s was simply inserted at the proper place
in that function; for traversal of the set of identifiers we use the |hash|
table, by following all non-empty hash lists.

@< Reverse the \xr. lists... @>=
{ id_pointer name; id_pointer *h; /* pointer into |hash| */
  for (h=hash; h<hash_end; h++)
    for (name=*h; name!=NULL; name=name->hash_link)
      if (name->ilk!=header_file_name)
      /* traverse hash lists, except for header file names */
    @< Reverse the list |name->xref| @>
}

@ As Knuth keeps reminding us, list reversal can be thought of as a process of
repeatedly popping values off one list~|x| and pushing them onto the reversed
list~|y| (or you may read ``stack'' for ``list'' if you like). It can also be
useful to remember that the basic action can be performed by a four-stroke
engine, where the left hand side of each assignment equals the right hand side
of the previous one. The basic cycle can actually take different forms, each
using an auxiliary variable~|t|. One way is to use~|t| to hold the entry moved
(``hold the head''), repeating |{ t=x; x=t->next; t->next=y; y=t;}| until
|x==NULL|; another way is to use~|t| to hold the remainder of the list to be
reversed (``hold the tail''), repeating |{ t=x->next; x->next=y; y=x; x=t;}|,
again until |x==NULL|. For reversing the \xr. lists of identifiers, we use
hold-the-head.

@< Reverse the list |name->xref| @>=
{
  sixteen_bits x=xref_index(name->xref),t,y=0;
    /* index of the sentinel node */
  while (xnum(x)!=0)
   {@; t=x; x=xlink(t); xlink(t)=y; y=t; }
  name->xref=&xmem[y]; /* don't forget to link in the reversed list */
}

@ The reversal of sublists of the \xr. list attached to module names is only
slightly more complicated. At the three places where the code below is used,
things have been set up so that |q| points to the location of the link
pointing to the start of the sublist (since this link is going to be changed,
we need a pointer to it) and the end of the list is implicitly indicated by
|threshold|: the sublist ends before the first entry with |num<=threshold|
(which always exists because of the sentinel with |num==0|). It has also been
ensured that the code is only invoked when the indicated sublist is not empty,
so that we can use a |do|-|while| loop; we have chosen the hold-the-tail
paradigm, mainly to demonstrate the alternative possibility, although it also
allows the code to be slightly shorter here.

After the sublist has been reversed, some links must be redirected to install
it in its proper place. The link |*q| must be pointed to the head of the
reversed list, which is in |y|, while the link at the end of the sublist must
be pointed to the unaffected remainder of the list (this remainder should
ideally have been assigned to |y| initially, but it is only located once we
have arrived at the first entry with |num<=threshold|). Fortunately the final
node of the sublist is not only pointed to by its predecessor, but also by
|*q| (before it is changed) since it used to be the first node of the sublist;
therefore a small sequence of carefully ordered assignments will do the trick.
It is instructive to check that if the sublist to be reversed has length~$1$,
then all variables will eventually return to their original state, except that
|q| is advanced to point to the |xlink| field of the single node. The
initialisation of~|y| is only present to keep certain compilers from
complaining that its value is used before it is first assigned to; the initial
value is irrelevant since it will be overwritten.

@< Reverse sublist... @>=
{ sixteen_bits x=*q,y=0,t;
  do {@; t=xlink(x); xlink(x)=y; y=x; } while (xnum(x=t)>threshold);
  xlink(t=*q)=x; *q=y; q=&xlink(t);
}


@*1 The problem of typedef declarations. @:title@>
%
We now consider how |typedef| declarations are processed. It is a slightly
problematic matter, since it involves the syntactic structure of the program,
but we are doing only lexical analysis during Phase~I. Nevertheless, we don't
want to delay recognition of |typedef| declarations to Phase~II, since this
would restrict the user's freedom of ordering the sections so that all typedef
declarations remain before any of their uses; and even if the user accepts
such a discipline, it would make it cumbersome to mention a typedef identifier
in the commentary directly before its definition as it is quite natural to do.
Our approach will be to specify simple rules that pinpoint the identifiers
which are subject to a |typedef| definition in any syntactically correct
program; we don't care too much about strange behaviour in the presence of
syntax errors, and also assume some basic decency on the part of the user,
e.g., the identifier should occur in the same section as the |typedef| token.
On the other hand this code is active while scanning header files whose
contents the user cannot always control, so it is required to be quite robust
when confronted with plain (non-literate) \Cee\ or \Cpp~code. The main
non-negotiable point is avoiding the use of the same identifier as ordinary
identifier in some places and as type identifier elsewhere; to keep those
apart, \.{\me.} would need to understand all scope and name-space rules of the
language, which would mean integrating functionality of both \.{CTANGLE} and of
a \Cpp\ compiler, and it obviously is not up to that task.

There is little to care about until a |typedef| token (recognisable by its
|ilk|) comes along, but then it becomes tricky. Ordinarily the identifier
being defined is the first non-type identifier that follows, but |struct| (and
|union| and |enum|) tokens complicate the situation; also there may be more
than one identifier subject to the |typedef|, separated by commas from each
other. Because of |struct| we must keep track of the nesting of braces, and
because the mentioned commas should be distinguished from those occurring in
the parameter specifications of functions, we should also keep track of
parenthesis nesting. A semicolon at the proper level of brace nesting signals
the end of a |typedef| declaration.

In \Cpp\ there is an extra complication due to template parameters. Between
the angle brackets `\.<' and `\.>' we should suspend our search for the
identifier being defined, just like we do between braces, and such brackets
can be nested. We shall use the same counter as for braces to handle angle
brackets for \Cpp.

Fortunately |typedef| declarations cannot be nested inside each other,
so we can use global variables to keep the proper counts. Three integer
counters are used, two for the nesting levels of braces and parentheses,
which are set to~0 whenever a |typedef| is scanned and properly maintained
thereafter, and a master counter which determines if we are paying
attention at all.

@d typedef_tracking(b) (typedef_master += b ? 5 : -5)

@< Global... @>=
local int typedef_master=-5; /* tracking disabled outside \Cee~parts */
local int brace_level, par_level;

@ The master counter is ordinarily equal to~0 when tracking is enabled, and
negative if it is disabled. When it is~0 and a |typedef| is seen, it is raised
to~2, after which an |int_like| or |type_defined| identifier will further
raise it to~4, indicating that any |normal| identifier coming along at the
same brace level will be made |type_defined|. When that happens the master
counter drops to~1, indicating that it can still be rekindled by a comma at
the proper parenthesis level. We must be prepared to see more than one typedef
for the same identifier, as we may have scanned a header file before seeing
the source that produced that file, so we also let the master counter drop
to~1 if it was~4 and an identifier that is already |type_defined| is seen. If
however following the |typedef|, when the master counter is~2, a |struct_like|
identifier is seen, the master counter is raised only to~3, so that a
following identifier will not be made |type_defined|, but rather pass the
honour on by setting the master counter to~4; this also happens when instead
of an identifier a left brace is seen. Finally, a semicolon at the right brace
level will return a positive value of the master counter to~$0$; if the master
counter is not~1 at this point no identifier has been marked, and the user is
warned.

In \Cpp\ we must track angle brackets, as said above, but there is one more
complication. In a typedef definition like $\&{typedef}\
\&{Some\_class}{::}\&{type\_member}\ \&{new\_type};$ we must not make the
identifier \&{type\_member} following |colon_colon| type-defined (which it
already is), although it is encountered at level~$4$ (the class name
\&{Some\_class} has raised it to that level), but wait for \&{new\_type}. To
this end we undo the effect of seeing the class name by setting
|typedef_master| back to the level~$2$ when the |colon_colon| is seen.

@< Keep track of tokens relevant to |typedef| declarations @>=
{ if (typedef_master==0 &&
      next_control==identifier && cur_id->ilk==typedef_like)
  @/{@; typedef_master=2; brace_level=par_level=0; }
  else if (typedef_master>0) switch(next_control)
  { case identifier:
      if (brace_level==0)
	if (typedef_master==2)
	{ if (cur_id->ilk==int_like || cur_id->ilk==type_defined)
	    typedef_master=4;
	  else if (cur_id->ilk==struct_like) typedef_master=3;
	}
	else if (typedef_master==4)
	{ if(cur_id->ilk==normal||cur_id->ilk==type_defined) /* this is it */
	    cur_id->ilk=type_defined, typedef_master=1;
	}
	else if (typedef_master==3) typedef_master=4;
      break;
    case '{': @+
      if (brace_level++==0 && typedef_master==3) typedef_master=4; @+ break;
    case '}': --brace_level; break;
    case '<': @+ if (C_plus_plus) brace_level++; @+ break;
    case '>': @+ if (C_plus_plus) --brace_level; @+ break;
    case ',': @+
      if (typedef_master==1 && par_level==0) typedef_master=4; @+ break;
    case '(': ++par_level; break;
    case ')': --par_level; break;
    case ';': @+
      if (brace_level==0)
      { if (typedef_master>=2)
          @< Issue a warning about an unrecognised typedef @>
        typedef_master=0;
      } @+
      break;
    case colon_colon:
       @+ if (C_plus_plus && brace_level==0 && typedef_master==4)
       typedef_master=2; @+ break;
  }
  if (C_plus_plus)
  @< Take action to mark identifiers following \&{class} as |type_defined| @>
}

@ When we see a semicolon that terminates a |typedef| definition is seen, it
should normally be the case that |typedef_master==1|, marking that fact that
the (final) identifier being type-defined has in fact been seen. If this is
not the case, then we issue a warning, since a failure to recognise the type
definition may cause a lot of syntactical trouble later on with the identifier
involved. Since this situation can easily arise while scanning a header file,
we make the error message more understandable by reporting the file name if it
is not the main file.

@< Issue a warning about an unrecognised typedef @>=
{ if (including_header_file)
    print("In file %s:\n",cur_file_name);
  print("\nUnrecognised typedef at line %d in section %d:\n"
	 @.Unrecognised typedef...@>  ,cur_line, section_count);
  mark_harmless();
}

@ In \Cpp, any identifier that follows \&{class} (or |struct|, |union|, or
|enum|) is considered as a typedef identifier, i.e., every time some
\&{class~x} is encountered, an implicit `|typedef| \&{class~x~x}' is assumed.
Since (except in contrived examples) the identifier is the one immediately
following the |class| keyword, we can use a static two-state variable to
signal that an implicit type definition is pending.

There is a problem with the keyword \&{typename} though. If it occurs to
introduce a template parameter, then that template parameter should be
interpreted as type so that the template body can be properly parsed; this is
what would normally happen, just like when one would write \&{class} instead
of \&{typename}. If however it occurs to force the following expression to be
interpreted as a type name, then that expression will probably involve some
namespace resolution, and we do not want to do anything to the {\it first\/}
identifier in that expression (if it is template parameter or a class name
then it has category |int_like| already, but it could be a namespace
identifier as in \&{typename}~$\\{std}\CC\&{vector}\ang<\&T\ang>
\CC\&{iterator}~p$, and we do not want to make \\{std} a type name!
Unfortunately the only way to tell these two uses of \&{typename} apart during
phase~I is to look ahead one more token, which means we must add two more
states to our state variable. Upon seeing \&{class} we move to state~$1$ which
will simply make the next token |type_defined| if it is an identifier. Upon
seeing \&{typename} we move to state~$2$, and a successive identifier advances
to state~$3$ while storing its identity in another static variable |this_id|,
then if the following token (in state~$3$) is {\it not\/} the namespace
resolution operation |colon_colon|, then the save identifier will be made
|type_defined|.

Should the logic below fail and make an identifier after \&{typename}
|type_defined| that shouldn't, the literate programmer can as a last resort
interject a \:\v or \:+ control between the \&{typename} keyword and the
identifier, which brings the identifier out of reach (this used to be the only
mechanism to disable \&{typename}). This is not possible in a header file
being scanned, so as an additional precaution we do nothing when \&{typename}
seen in header files. This is not so bad, because template arguments serve a
local purpose in the header file anyway.

@< Take action to mark identifiers following \&{class}... @>=
{ static int class_seen=0; static id_pointer this_id;
  switch (class_seen)
  { case 0:
    if (next_control==identifier)
      if (cur_id->ilk==struct_like) class_seen=1;
      else if(cur_id->ilk==typename_like)
        class_seen=2;
  break;
    case 1:
    if (next_control==identifier && cur_id->ilk==normal)
      cur_id->ilk=type_defined;
    class_seen=0;
  break;
    case 2:
    if (next_control==identifier && cur_id->ilk==normal)
      {@; this_id=cur_id; class_seen=3; }
    else class_seen=0;
  break;
    case 3:
    if (next_control!=colon_colon) this_id->ilk=type_defined;
    class_seen=0;
  }
}

@ We now consider |typedef| declarations in Phase~II. In Phase~I we have set
the |ilk| of the defined identifiers to~|type_defined|, which will make them
behave as |int_like|; although this works fine everywhere else, it thwarts a
correct parse of their |typedef| declaration itself during Phase~II. In
fact, according to the \caps{ANSI}/\caps{ISO}~\Cee\ syntax, the identifier
being declared in a typedef declaration should not be considered to be a
\\{typedef-name}, since that would make the declaration unsyntactical; this
is justified by the fact that the typedef identifier only comes into scope
after the declaration is completed. This is not a problem for single-pass
compilers, but it is for us: we have recorded the typedef declaration on the
first pass, and it will be active during the entire second pass, since we
cannot delimit the declaration to its proper range (we are not processing
the code in the same order that the \Cee~compiler will, and for the
commentary parts of sections, the concept of range does not even apply). To
attempt to write the grammar in such a way that it will accept typedef
declarations in which the defined identifier is |int_like| would be very
difficult, since without the help of some rather remote context, a
declarator of this kind can not always be distinguished form an abstract
declarator; compare the declaration `\hbox{|typedef char *(example[4]);|}',
which declares |example| to specify the type ``array of 4 pointers to
character'', with a valid declaration `|void f(example[3]);|' that might
follow it, declaring |f| as a function that takes as an argument an array of
3 such |example| objects (in fact the argument will be passed as a pointer to
|example|, of course).

So, rather than solving the problem in a syntactic way, we stoop down to
emulating a one-pass system by setting the |category| of the defining
occurrence of an identifier in a typedef declaration explicitly to
|expression|, despite the fact that its |ilk| is |type_defined|. The defining
occurrence is located by the same lexical means used in Phase~I, in fact, by
using the same intricate succession of states of |typedef_master| that was
used there. A difference is that at the {\sl moment supr\`eme\/} the
identifier found is now |type_defined| rather than |normal|. The fact that
this identifier is not necessarily the only or first |type_defined|
identifier in the declaration, and that is has to be recognised in a
left-to-right pass, may explain some of the details of our code, for
instance why |const_like| had to be distinguished from both |int_like| and
|type_defined|; one may compare `\hbox{|typedef unsigned long int * const
example;|}' with `\hbox{|typedef const example * examp2;|}'.

With respect to the code for Phase~I, there is also a slight difference in
the way the code is hooked into the program, since the function |C_read|
already contains a |switch| on the value of |next_control|, which has a
special case for identifiers. So let us first consider the identifier cases
(including reserved words).

@f example scrap /* pretend |example| is |type_defined| */
@f examp2 example /* and |examp2| is another such identifier */

@< Track identifiers... @>=
{ if (typedef_master==0 && cat==typedef_like)
    typedef_master=2, brace_level=par_level=0;
  else if (typedef_master>0 && brace_level==0)
    if (typedef_master==2)
    { if (cat==int_like || cat==type_defined) typedef_master=4;
      else if (cat==struct_like) typedef_master=3;
    }
    else if (typedef_master==4 && cat==type_defined) /* this is it */
      cat=expression, typedef_master=1;
    else if (typedef_master==3) typedef_master=4;
}

@ And here are the relevant cases of non-identifiers.

@< Cases necessary for proper parsing of typedefs... @>=

case '{': @+
  if (typedef_master>0 && brace_level++==0 && typedef_master==3)
    typedef_master=4;
  @+break;
case '}':@+ if (typedef_master>0) --brace_level; @+break;
case '<':@+ if (C_plus_plus && typedef_master>0) ++brace_level; @+break;
case '>':@+ if (C_plus_plus && typedef_master>0) --brace_level; @+break;
case ',':@+ if (typedef_master==1 && par_level==0) typedef_master=4; @+break;
case '(':@+ if (typedef_master>0) ++par_level; @+break;
case ')':@+ if (typedef_master>0) --par_level; @+break;
case ';':@+ if (typedef_master>0 && brace_level==0) typedef_master=0; @+break;
case colon_colon:
       @+ if (C_plus_plus && typedef_master>0 && brace_level==0)
          typedef_master=2 ; @+ break;

@* Phase II processing. @:title@>
With the description of Phase~I still fresh in our memory, let us look at
the general outline of Phase~II, which is analogous, although it is more
complicated. The extra complication is due to the fact that much more has to
be done during Phase~II, notably the \Cee~texts have to be parsed in order
to determine their proper formatting, and the resulting token lists have to
be written to a file in a form that \TeX\ will be able to process. Most of
the actual work however is localised in a few powerful functions that will
be defined in detail later on, so that the definition of |phase_two| can be
given here without much problems.

There are three stages in the processing of a piece of \Cee~text during
Phase~II: first it is scanned lexically (using |get_next|), and the
resulting tokens are collected in the form of `scraps' that form the input
for the second stage, the parsing algorithm, which transforms the scraps into
a recursively nested token list, which is converted in the third stage to
textual output.
For the small parts of \Cee~text enclosed in `\pb', a function~|do_C| is
available which handles all three stages of processing: it is called when an
opening `\.\v' is seen, and when it is completed one has arrived at the
closing `\.\v', and the required output is written on the \TeX~file.
For the larger parts of \Cee~text, the function |outer_read| will read in,
and convert to scraps, chunks of \Cee~text delimited by control codes
|c>=format| (so like |outer_xref| it handles comments, but it will not
incorporate module names), and when enough of these have been accumulated,
a call on |finish_C| will invoke the parsing algorithm and send the
resulting tokens to the output file.

A number of simple functions for producing output are also called explicitly
at certain points. We have already seen |out| and~|finish_line| for
character-based output; there is also |out_str| for writing a string,
|out_sec_nr| for a section number, |list_refs| for generating the text for
\:\#, and |footnote| for producing the module \xr. information at the end of
sections. Invoking the macro |tex_new_line| immediately after |finish_line|
was called will produce an empty line on the output. Furthermore certain
small pieces of code which have been scanned directly rather than via |do_C|
or |outer_read| (for instance after format or macro definitions) are
converted into scraps: first a number of tokens are appended by means of
|app| or~|app_str| and then the whole sequence is converted to a scrap by
calling |pack_scrap|; |app| and~|pack_scrap| are macros.

@< Prototypes @>=
void do_C (void); /* handle \Cee~text enclosed in `\pb' */
void outer_read (void); /* transform input into scraps */
void finish_C (void); /* finishes a definition or a \Cee~part */
void finish_line(void); /* send out a line of output */
void out_str (char*); /* write multiple characters */
void out_sec_nr (int); /* output a section number */
xref_pointer list_refs (xref_pointer,sixteen_bits);
  /* output module \xr.s */
void footnote(xref_pointer*,sixteen_bits); /* same with heading text */
void app_str(char*); /* append a sequence of character tokens */

@ Like in |phase_one|, we loop over the sections after passing over the limbo
part.

@c
void phase_two (void)
   /* read all the text again and translate it to \TeX\ form */
{ phase=2; reset_input ();
  print_progress("\nWriting the output file...");
		  @.Writing the output file...@>
  section_count=0; copy_limbo(); finish_line();
  tex_new_line(); /* insert a blank line, it looks nice */
  while (!input_has_ended) @<Translate the current section@>
}

@ The output file will contain the control sequence `\.{\\Y}' before a
non-empty definition portion of a section, and before a non-empty
\Cee~portion (the \TeX~portion is always considered to be non-empty, since
it contains at least the section number). This puts a little white space
between adjacent portions when they are printed.

@d emit_space() out_str ("\\Y"); @.\\Y@>

@< Translate the current section @>=
{ section_count++;
  @< Output the code for the beginning of a new section @>
  @< Translate the \TeX~part of the current section @>
  if (next_control<begin_C)
  @/{@; emit_space();
    @< Translate the definition part of the current section @>
  }
  if (next_control<new_section)
  { mod_pointer this_module=NULL; /* the current module name */
    emit_space(); @< Translate the \Cee~part of the current section @>
    @< Show cross-references to this section @>
  }
  @< Output the code for the end of a section @>
}

@ Sections beginning with the \.{CWEB} control sequence \:{\ } start in the
output with the \TeX\ control sequence `\.{\\M}', followed by the section
number. Similarly, \:* sections lead to the control sequence `\.{\\N}', and
\:\~ sections to `\.{\\n}'. If this is a changed section, we put `\.*' just
before the section number.

@< Output the code for the beginning... @>=
{ out('\\'); out(loc[-1]=='*' ? 'N' : loc[-1]=='~' ? 'n' : 'M' );
  @.\\N@> @.\\n@> @.\\M@>
  if (loc[-1]=='*')
  {@; print_section_progress(); @< Handle title level @>@+ }
  out_sec_nr(section_count); out_str(". ");
}

@ Between \:* and the title that follows, a level can be specified in the
form of another `\.*', or a decimal number; the absence of a number will be
interpreted as level~0. The level will be written out after `\.{\\N}' as a
first argument, delimited by a space, after which the second argument is the
section number and the third specifies the title.

@< Handle title level @>=
{ if (*loc=='*') ++loc,out_str("-1");
  else if (!isdigit((eight_bits)*loc)) out('0');
  else do out(*loc++); while (isdigit((eight_bits)*loc));
  out(' '); /* terminate level by a space */
}

@ In the \TeX~part of a section, we simply copy the source text, except that
index entries are not copied and \Cee\ text within `\pb' is translated;
during this translation we track typedef definitions so that any complete
typedef declaration within `\pb' will be parsed correctly, as will be
explained below.

@< Translate the \TeX... @>=
do
  switch (next_control=copy_TeX())
  { case '|': typedef_master=0; do_C(); break;
    case at_sign_image: out('@@'); break;
    case thin_space: case math_break: case ASCII_code: case line_break:
    case big_line_break: case no_line_break: case join: case pseudo_semi:
    case force_expr_open: case force_expr_close:
      err_print("! You can't do that in TeX text");
		 @.You can't do that...@>
      break;
#ifdef DEBUG
    case trace0: case trace1: case trace2: case trace3: tracing=next_control;
      break;
#endif
    case module_name: loc-=2; get_next(); break; /* get module name */
    case refer: loc-=2; get_next(); /* get name referred to */
      if (cur_id->xref->num==0) err_print("! Undefined reference");
      else list_refs(cur_id->xref,0);
      break;
    case TeX_string: err_print("! TeX string should be in C text only");
				@.TeX string should be...@>
    /* fall through */
    case xref_roman: case xref_wildcard: case xref_typewriter:
    case xref_mark: case ignored_text:
      get_control_text(); /* skip to \:> */
  }
while (next_control<format);

@ When we get to the following code we have |format<=next_control<begin_C|.
The first conditional statement below processes the first few tokens of a
preprocessor directive, which do not follow the ordinary syntax rules; the
remainder is then treated by |outer_read| using the regular parsing
mechanism, although in the case of a format definition this should only
involve possible comments and formatting controls, and in case of \:s even
these should be absent.

@< Translate the def... @>=
{ typedef_tracking(false);
  do
  { boolean suppressed=false; /* whether output suppressed by \:s */
    if (next_control==format) @< Start a format definition @>
    else if (next_control==definition) @< Start a macro definition @>
    else @< Start a header file inclusion @>
    if (!suppressed) outer_read(), finish_C();
    else if (next_control<format)
    { err_print("! Improper stuff after `@@s' format definition");
		 @.Improper stuff after `@@s'...@>
      if (next_control==begin_comment) loc-=2; /* try to get back in phase */
      outer_xref(); /* skip illegal stuff */
    }
  } while (next_control<begin_C); /* |format|, |definition|, or |header| */
}

@ The syntax of a format definition has already been checked. It suffices to
build scraps that will produce the desired output. The trickiest point is
formatting the identifiers in a way that is least confusing to readers.  We
try to keep analogy to lines for macro definitions in case these define one
identifier to stand for another: the left hand side is set in italics
regardless of its |ilk|, while the right hand side has its normal
appearance, prescribed by its |ilk|. If the left hand side is a |TeX_like| or
|NULL_like| identifier, we stick to the representation in italics, even though
this costs some extra work, and is not what would happen for a macro
definition; in addition we append the usual formatted form of the identifier
in parentheses, to show the correspondence of name and printed symbol.

@< Start a format... @>=
if (tolower((eight_bits)loc[-1])=='s')
@/{@; suppressed=true; shift(); shift(); shift(); }
  /* skip format definition */
else
{ int saved_code=0,saved_mathness;
  app_str("\\F"); shift(); /* this will produce `\&{format}' */ @.\\F@>
  if (cur_id->ilk!=TeX_like && cur_id->ilk!=NULL_like)
    app(id_flag+id_index(cur_id));
  else @< Expand identifier and set |saved_code| and |saved_mathness| @>
  app('~'); pack_scrap(insert,yes_math); shift();
  app((cur_id->ilk==normal || cur_id->ilk==TeX_like || cur_id->ilk==NULL_like
      ? id_flag : res_flag
      )+id_index(cur_id));
@/check_scrap();
  pack_scrap(insert,cur_id->ilk==TeX_like ? no_math : yes_math);
  shift();
  if (saved_code!=0)
  { app_str("\\quad("); app(saved_code); app(')');
  @/check_scrap(); pack_scrap(insert,saved_mathness);
  }
}

@ Since conversion of an identifier into a \TeX\ control sequence is
performed by the output routine, we need to circumvent this to force italic
type; this is done by expanding the name into characters directly rather
than leaving this to the output routine.

@< Expand identifier... @>=
{ char* p=name_begin(cur_id);
  saved_mathness=cur_id->ilk==TeX_like ? no_math : yes_math;
  saved_code=id_flag+id_index(cur_id);/* save to print afterwards */
  app_str("\\\\{"); @.\\\\@>
  do {@; if (*p=='_') app('\\'); app_tok(*p); } while (*++p!='\0');
  app('}'); check_toks(10);
}

@ Keeping in line with the conventions of the \Cee\ preprocessor (and
otherwise contrary to the rules of \.{CWEB}) we distinguish here between the
cases that a `\.(' immediately follows the identifier being defined, and the
case that anything else (possibly a space) does. In the latter case, the
replacement text starts immediately after the identifier, in the former
case, it starts after we scan the matching `\.)', which must be simply the
first `\.)' that follows.

@<Start a macro...@>=
{ if (shift()!=identifier)
    err_print("! Improper macro definition");
	       @.Improper macro definition@>
  else
  { app_str("\\D$"); @q $ emacs-cookie @>
          /* this will produce \&{\#define} */ @.\\D@>
    app(id_flag+id_index(cur_id));
    if (*loc=='(')
    { shift();
      do
      { app_char_tok(next_control);
	if (shift()!=identifier) break;
	app(id_flag+id_index(cur_id));
      } while(shift()==',');
      check_toks(2);
      if (next_control==')') {@; app(')'); shift(); }
      else err_print("! Improper macro definition");
    }
    else shift();
    app('$');@q $ emacs-cookie @>
    app(break_space); pack_scrap(insert,no_math);
  }
}

@ For scanning the token following \:h we temporarily set |preprocessing=2|,
so that angle brackets will be recognised as string quotes.

@<Start a header file...@>=
{ app_str("\\h"); /* this will produce \&{\#include} */ @.\\h@>
  pack_scrap(insert,no_math);
  { int save=preprocessing; preprocessing=2; /* emulate `\.{\#include}' */
    while (shift()==' ') {} /* skip spaces and read file name as string */
    preprocessing=save;
  }
}

@ Finally, when the \TeX\ and definition parts have been treated, we have
|next_control>=begin_C|. If the section defines a module name, we assign the
name to the variable |this_module|, so that the proper \xr. information can
be listed at the end of the section. Like in Phase~I it is necessary to pay
special attention to |typedef| declarations, and this time tracking is
enabled both within `\pb' in \TeX~text and in the \Cee~part of a section;
since the former may involve incomplete pieces of syntax like a sole
`|typedef|', we reset the master counter to its neutral state at the
beginning of a \Cee~part. After the heading of the \Cee~text has been
processed we alternatively read ordinary pieces of \Cee~text and module
names until the module has ended; we start with calling |outer_read| before
testing termination, in order to ensure that the overflow tests contained in
|outer_read| will be executed even in case of a section with a \Cee~part
consisting only of a heading.

@<Translate the \Cee...@>=
{ typedef_master=0;
  if (next_control==begin_C) shift();
  else
  { this_module=cur_mod; /* register the name for this module */
    @< Check that `\.{=}' or `\.{==}' follows this module name, and
       emit the scraps to start the module definition @>
  }
  do
  { outer_read();
    if (next_control==new_section) break;
    if (next_control==module_name) @< Append a module name scrap @>
    else err_print("! You can't do that in C text");
		    @.You can't do that...@>
      /* |format|, |definition| or |begin_C| */
    shift();
  } while (true);
  finish_C();
}

@ Despite the name of this module, we allow `\.+' to precede the `\.{=}' or
`\.{==}', just as |CTANGLE| does. Note however that, unlike in |CTANGLE|,
the `\.+' will be scanned as part of an `\.{+=}' compound operator (whence
in fact `\.{+= =}' is allowed here whereas `\.{+ ==}' is not; we hope
nobody minds this). Note also that if for whatever reason the current
section number should fail to appear in the \xr. list for the
module name (e.g., if section numbers have inadvertently got out of
synchronisation with respect to Phase~I), then listing the
\xr. information at the end of the section is suppressed by setting
|this_module=NULL|.

@< Check that `\.{=}' ... @>=
{ if (shift()=='=' || next_control==eq_eq || next_control==plus_assign)
  @/{@; if (next_control!=plus_assign || shift()=='=') shift(); }
      /* accept `\.=', `\.{==}', `\.{+=}' or `\.{+==}' */
  else err_print("! You need an = sign after the module name");
		  @.You need an = sign...@>
  if (this_module!=NULL) /* i.e., unless module name was bad */
  { xref_pointer x=this_module->xref;
    if (x->num==file_flag) x=next_xref(x);
    app_str("\\4$"); @q $ emacs-cookie @>
                     /* module name will be flush left */ @.\\4@>
    app(mod_flag+mod_index(this_module));
    if (x->num != section_count+def_flag)
    { app_str("\\PE"); /* module has also been defined before */ @.\\PE@>
      this_module = NULL; /* so we won't give \xr. info here */
    }
    else app_str("\\EQ"); /* output a module definition sign */ @.\\EQ@>
    app_str("{}$"); @q $ emacs-cookie @>
    app(force); pack_scrap(insert,no_math);
      /* this forces a line break unless \:+ follows */
  }
}

@ Cross references relating to a named module are given after its first
defining section ends (for further defining sections of this name we will
have put |this_module=NULL|).

@< Show cross-references... @>=
{ if (this_module != NULL)
  { xref_pointer foot_ref=this_module->xref;
    if (foot_ref->num==file_flag) foot_ref=next_xref(foot_ref);
    foot_ref=next_xref(foot_ref); /* don't \xr. to yourself */
    footnote(&foot_ref,def_flag);
      /* display further defining sections; advance |foot_ref| */
    footnote(&foot_ref,cite_flag); /* display any citations */
    footnote(&foot_ref,0); /* display uses */
  }
}

@ The `\.{\\fi}' closes a \TeX~conditional that was initiated by the macro
that started off the section; this allows printing of only the changed
sections in a simple way.

@<Output the code for the end of a section@>=
{@; out_str ("\\fi"); finish_line ();  tex_new_line(); }
@.\\fi@> /* insert a blank line, it looks nice */


@*1 Auxiliary functions used in Phase~II. @:title@>
We now define the functions that do the actual processing of \Cee~code
during Phase~II, but without going into the details of parsing and output.
We explain the functions |do_C|, |outer_read|, |finish_C|, |list_refs|,
and |footnote| used above, and also two further auxiliaries |C_read|
and~|C_translate|.

@< Prototypes @>=
text_pointer translate(void);
  /* build formatted text from collected scraps */
void make_output(text_pointer,mode);
  /* output text in |inner| or |outer| mode */

@ Before we discuss these functions, we must first discuss what scraps are.
Scraps are the objects manipulated during parsing, and they have two main
attributes: a syntactic category |cat|, that determines the way they will be
treated by the parser, and a translation |trans|, which is a pointer
into~|text_mem| denoting a (possibly recursively nested) sequence of tokens,
that determines the representation of the scrap upon output. Since some
parts of the output are to be processed in \TeX's math mode, and other parts
in horizontal mode, an additional field |mathness| tells which mode is
required at each end of the translation of the scrap.

@<Typedef...@>=
typedef struct
{ eight_bits cat; /* category code */
  eight_bits mathness; /* whether in math mode at left and right boundary */
  text_pointer trans; /* translation text */
} scrap, *scrap_pointer;

@ When \Cee\ text is converted into scraps for parsing, the resulting scraps
are placed in an array |scrap_info|, between the locations pointed to by
|scrap_base| and~|scrap_ptr|. Actually, |scrap_info| is one field of a
|union|, since the same memory is used for a different purpose during
Phase~III.

Basic scraps are created by invoking |app|, |app_tok| or |app_char_tok| a
number of times creating the constituent tokens, and then consolidating the
text by means of |freeze_text|; the resulting text is accessible as
|text_ptr| before, and as |text_ptr-1| after the call of |freeze_text|. In
the common case that the text forms the translation of a new scrap that is
to be added to the scrap sequence, |pack_scrap| can be used in place of
|freeze_text|; a category and `mathness' should be supplied in this case.
The latter can take one of three values as explained later, and by
multiplying it by~5 (binary~$0101$) it is duplicated into the two least
significant pairs of bits, because for elementary scraps the value is the
same at its left and right boundaries. Note that none of |app|, |freeze_text|
and |pack_scrap| do bound checks, since it is assumed that these have been
done beforehand; |app_tok| and |app_char_tok| however can be called at more
uncertain times.

@d scrap_info scrap_union.scrap_field
@d scrap_info_end  (&scrap_info[max_scraps]) /* end of |scrap_info| */
@)
@d app(a) (*tok_ptr++ = a)
@d app_tok(a) @+
  if (tok_ptr>tok_mem_end-2) overflow("token"); @.token capacity exceeded@>
  @+ else app(a) @;
@d app_char_tok(c) app_tok((unsigned char)(c))
@d freeze_text() (*++text_ptr = tok_ptr)
@d pack_scrap(c,m)
 ( scrap_ptr->cat = c, scrap_ptr->trans = text_ptr, freeze_text(),
  (scrap_ptr++)->mathness = 5*(m) )

@<Global...@>=
union
{ scrap scrap_field[max_scraps]; /* memory array for scraps */
  @< Alternative use of |scrap_union| @>@;
} scrap_union;
scrap_pointer scrap_base=scrap_info;
	/* beginning of the current scrap sequence */
scrap_pointer scrap_ptr = scrap_info;
	/* points to end of the current scrap sequence */
#ifdef STAT
scrap_pointer max_scr_ptr = scrap_info;
	/* largest value assumed by |scrap_ptr| */
#endif

@ Token lists are stored in |tok_mem| and represent text to be output to the
\TeX~file. Because during parsing token lists will often be formed by
concatenation of existing ones, a representation is chosen where this can be
done easily; in particular a token can be a reference to a token list stored
elsewhere. Also identifiers, reserved words and module names are represented
by a reference to the name table rather than by their constituent
characters. All other items are stored as list of characters, which have
been widened to fill a 16-bit |token|, and special layout codes that will be
explained below. More precisely, a |token t@;| is interpreted as follows.
\yskip

\item{$\bullet$} |t<=UCHAR_MAX|: the character~|t|, which possibly is
     a compressed operator like `\.{\&\&}';
\item{$\bullet$} |UCHAR_MAX<t<id_flag|: a special layout feature such as
     |indent|,
\item{$\bullet$}|id_flag<=t<res_flag|: the identifier with name
     |id_at(t-id_flag)|;
\item{$\bullet$}|res_flag<=t<mod_flag|: the reserved word
     |id_at(t-res_flag)|;
\item{$\bullet$}|mod_flag<=t<text_flag|: the module named
     |mod_at(t-mod_flag)|;
\item{$\bullet$}|text_flag<=t<inner_text_flag|: a reference to the token list
     |text_at(t-text_flag)|;
\item{$\bullet$}|inner_text_flag<=t|: a reference to the token list
     |text_at(t-inner_text_flag)|, which is to be translated without
     line-break controls.

@d id_flag 10240U /* signifies an identifier */
@d res_flag (2*id_flag) /* signifies a reserved word */
@d mod_flag (3*id_flag) /* signifies a module name */
@d text_flag (4*id_flag) /* signifies a token list */
@d inner_text_flag (5*id_flag) /* signifies a token list in `\pb' */

@~The special layout tokens are the following:
\yskip

\item{$\bullet$} |indent| causes future lines to be indented one more unit;
\item{$\bullet$} |outdent| causes future lines to be indented one less unit;
\item{$\bullet$} |opt| denotes an optional line break, it is followed by an
     integer~|n|, and the break will normally occur with penalty~$10n$ (in
     fact \TeX\ raises some penalties, but this remains invisible to
     \.{\me.});
\item{$\bullet$} |flush_left| denotes that the line will be printed flush
     left, regardless of the current indentation level;
\item{$\bullet$} |break_space| denotes an optional line break or an ``en''
     space;
\item{$\bullet$} |force| denotes a forced line break;
\item{$\bullet$} |big_force| denotes a forced line break with additional
     vertical space;
\item{$\bullet$} |backup| denotes a forced line break followed by a
     backspace of one indentation unit;
\item{$\bullet$} |big_backup| denotes a forced line break with additional
     vertical space, followed by a backspace of one indentation unit;
\item{$\bullet$} |cancel| obliterates any space, |break_space|, |force|,
     |big_force|, |backup|, or |big_backup| tokens that immediately precede
     or follow it and also cancels any |opt| tokens that follow it;
\item{$\bullet$} |relax| does nothing, but serves as a ``stopper''.

\yskip\noindent The character tokens |' '|~and~|'~'|, together with the
sequence of tokens from |break_space| to |big_backup| form a hierarchy, in
the sense that on output any consecutive sequence of such tokens is
equivalent to their maximum, except that |big_force| and |backup| combine to
|big_backup|.  To facilitate the computation we added |space|~and~|tilde| to
the enumeration below, that overlay |opt|~and~|flush_left| respectively
(this causes no problems, since |space|~and~|tilde| are not used as token
values). Formally, we define a partial ordering on the set |{' ', '~',
break_space, force, big_force, backup, big_backup }| that differs from the
linear order in which they were listed only by the fact that |big_force| and
|backup| are incomparable; then any consecutive list of such tokens is
replaced by their least upper bound.

@< Typedef and enum... @>=
enum @/
{ cancel=UCHAR_MAX+1,/* the following 9 items should remain in this order */
  indent, outdent, opt, flush_left, break_space, force, big_force,
  backup, big_backup, @/
  relax, @/
  space=opt, tilde=flush_left @/
};

@ The memory management for tokens follows a simple block regime.
The source file is divided into blocks, and whenever a block is entered,
the current states of the |text_mem| and~|tok_mem| are marked, after
which they will gradually get filled up further; at the end of the block
all memory used during the block is released by resetting the pointers
into these arrays to the values they had on block entry. There are three
kinds of blocks: the global block, which is filled during initialisation and
is never released, the section blocks, which correspond to each individual
section with a non-empty \Cee~part (or to a macro or format definition),
and the inner blocks, which correspond to each `\pb' contained in the
\TeX~part of a section or in a module name (but not in a comment). Because
nesting is at most three blocks deep (for `\pb' inside module names) a
two-element stack will suffice to hold the saved markers (nothing needs to
be saved for the global block), and we can address these elements directly
without a stack pointer since we know at which level we are. In fact the
section blocks butt together, and nothing is added to the global block
except at initialisation time, so that each section block starts in the same
state. Therefore, if we call |enter_block(0)| after initialisation to record
this state, then there is no need to call it any more, and it suffices to
call |leave_block(0)| each time upon leaving a section.

Scraps follow a different regime than texts and tokens, since they must be
assembled into a single contiguous sequence before each translation, and can
be discarded when the translation is over. In particular the scraps used to
parse `\pb' within a comment can and must be released before reading in the
\Cee~text following the comment, but the tokens and texts formed while
translating the comment must remain until they are output.

@d enter_block(i) save[i].txt=text_ptr, save[i].tok=tok_ptr;
@d leave_block(i) text_ptr=save[i].txt, tok_ptr=save[i].tok;

@< Global variables @>=
struct {@; text_pointer txt; token_pointer tok; } save[2];

@ The conversion of input tokens as obtained by |get_next| into scraps
that can be processed by the parser is mainly handled by the function
|C_read|, which is analogous to the |C_xref| routine used during Phase~I.
Like |C_xref|, the function |C_read| takes a boolean argument telling
whether it is processing `\pb'; it starts with the current value of
|next_control| and it uses the operation |shift| repeatedly to read
\Cee~text until encountering a terminating token. Also like |C_ref|, the
initial conditions |next_control=='|'| or |next_control==end_comment| will
not lead to immediate termination; in these cases the first time through the
|while| loop will have no effect.

@c
void C_read (boolean inner) /* creates scraps from \Cee\ tokens */
{ while (next_control<format || next_control==module_name && inner)
  { @<Append the scrap appropriate to |next_control|@>
    if (shift()=='|' && inner
      || next_control==begin_comment || next_control==end_comment) return;
  }
}

@ Many input tokens are completely determined by the value returned from
|get_next|; we have arranged it that such a value is always less than
|ignore|. For those tokens we will install corresponding scrap at initialisation
time in an array |token_trans|, so that when such a token comes along, we
can simply copy the scrap from |token_trans| into scrap memory.

@< Global variables @>=
scrap token_trans[ignore];

@~Since the scraps in |token_trans| contain pointers into |tok_mem|, that
array cannot be initialised statically; rather we do the initialisation
dynamically based on information stored statically. For each token three
kinds of information are supplied: the category of its translation, a
string of characters giving the translation itself, and an indication for
its the |mathness| (whether the translation must or must not occur in math
mode, or whether both are allowed). The actual initialisation values will
be given later.

@< Set initial values @>=
{ static struct {@; short tok; eight_bits cat, mathness; char* tr; }
    trans_ini [] = @/{ @< Initialiser for |trans_ini| @>@;@; };@/
  int i,n=array_size(trans_ini);
  scrap o={ insert, 5*maybe_math, NULL }; /* completely inert scrap */

  for (i=0; i<n; ++i)
  { scrap* p=&token_trans[trans_ini[i].tok];
    p->cat=trans_ini[i].cat; p->mathness=5*trans_ini[i].mathness;
    app_str(trans_ini[i].tr); p->trans=text_ptr; freeze_text();
  }
  if (C_plus_plus) @< Fix the categories of angle brackets for \Cpp @>
  @< Install the translations of tokens involving line breaks @>

  o.trans=text_ptr; freeze_text(); /* empty translation */
  for (i=0; i<=UCHAR_MAX; ++i) /* clear all remaining tokens */
    if (token_trans[i].cat==0) token_trans[i]=o;
  enter_block(0); /* fix tokens; will be restored after each section */
}

@ Having installed |token_trans|, the number of cases in the
switch statement below is greatly reduced.

@< Append the scr... @>=
{ check_scrap(); check_toks(6); /* `\.{\\hbox\{}' */
  switch (next_control)
  { case string: case constant: case verbatim:
      @< Append a string or constant @> @+ goto done;
    case TeX_string: @< Append a \TeX\ string scrap @> @+ goto done;
    case identifier: @< Append an identifier scrap @> @+ goto done;
    case module_name: @< Append a module name scrap @> @+ goto done;
    case start_preproc:
      @< Append a scrap starting a preprocessing directive @> @+ goto done;
    case refer:
      err_print("! You can't use `@@#' in C text"); /*fall through */
    case ignore: case begin_comment: case end_comment:
    case xref_roman: case xref_wildcard: case xref_typewriter:
    case xref_mark: goto done;
  @\@< Cases necessary for proper parsing of typedefs,
       each followed by |break| @>
    case '|': @+ if (inner) goto done; /* skip initial `\.\v' of `\pb' */
  }
  *scrap_ptr=token_trans[next_control]; /* fixed scrap for this input token */
  @< Possibly scoop up some dangling output tokens in compatibility mode @>
  ++scrap_ptr; /* incorporate the scrap */
  done: {}
}

@ When we need to be sure there is enough space to store |n|~more tokens,
we say |check_toks(n)|; when a scrap has to be appended we invoke
|check_scrap|, and when a text is to be formed otherwise than by
|pack_scrap|, we invoke |check_text|. The key points in the program where
we make such checks is when reading in \Cee~text in the functions |C_read|
and |outer_read|, and when translating the text in |reduce|. At these points
we make sure that there is enough room to spare so that in cases where
explicit tokens occasionally need to be appended (e.g., by |app_str|) no
test is necessary. When we cannot be sure about this, for instance because
an indefinite number of tokens is appended, we use |app_tok| or
|app_char_tok| instead of |app|; after this it will in fact be safe to do
one additional |app|.

@d check_toks(n) @+ if (tok_ptr>tok_mem_end-n)
   overflow("token"); @.token capacity exceeded@> @+ else @;
@d check_text()  @+ if (text_ptr>=text_mem_end-1)
   overflow("text"); @.text capacity exceeded@> @+ else @;
@d check_scrap() @+ if (scrap_ptr>=scrap_info_end)
   overflow("scrap"); @.scrap capacity exceeded@> @+ else check_text() @;

@ As was just explained, we can use |app| rather than |app_char_tok| here.

@c
void app_str(char* s) @+{@; while(*s!='\0') app(*s++); }

@ In long strings we insert a discretionary break every 20~characters, so
that \TeX\ is less likely to run into problems, especially if the string
occurs in the \TeX~part of a section; if we are directly after a backslash
however, we postpone the break since otherwise the quote escaping the line
break would appear to be escaped itself. Many of the special characters in
a string must be prefixed by `\.\\' so that \TeX\ will print them properly.
In the case of constants however, this `\.\\' converts marker characters
that were inserted during scanning into control sequences used in formatting
the constant.
@^special string characters@>

@< Append a string or... @>=
{ int count = -1; /* characters remaining before string break */

  if (next_control==constant) app_str("\\T{"); @.\\T@>
  else if (next_control==string) {@; count=20; app_str("\\.{"); } @.\\.@>
  else app_str("\\vb{"); @.\\vb@>

  while (id_first<id_loc)
  { if (count--==0) /* insert a discretionary break in a long string */
      if (id_first[-1]=='\\') count=0; /* no break after backslash */
      else {@; check_toks(2); app_str("\\)"); count = 20; } @.\\)@>
    if (strchr(" \\#%$^{}~&_",*id_first)!=NULL) app('\\');
@.\\\ @> @.\\\\@> @.\\\#@> @.\\\%@> @.\\\$@> @.\\\^@>
@.\\\{@> @.\\\}@> @.\\\~@> @.\\\~@> @.\\\&@> @.\\\_@>
    app_char_tok(*id_first++);
  }
  app('}');
  if (next_control==verbatim) pack_scrap(insert,maybe_math);
  else pack_scrap(expression,yes_math);
}

@ A \TeX~string is boxed and copied without further ado; undoubling of any
\:@@ has already been done by |get_control_text|. In compatibility mode we
do not however produce a scrap, but rather leave the output tokens produced to
be picked up by the next scrap appended; this makes |dangling_tokens()| hold.
In order to make sure that this next scrap exists, we will append (in
|translate|, section @#translate definition@>) a dummy scrap if necessary.
@:TeX string@>

@d dangling_tokens() (compatibility_mode && tok_ptr>*text_ptr)

@< Append a \TeX\ string scrap @>=
{ app_str("\\hbox{");
  while (id_first<id_loc) app_char_tok(*id_first++);
  app('}');
  if (!compatibility_mode) pack_scrap(expression,maybe_math);
}

@ Although dangling tokens left behind by the code in section@#TeX string@> in
compatibility mode are usually picked up automatically while constructing the
next scrap, this will not be the case if the next input token has a
prefabricated scrap, since |freeze_text| is not called for such tokens. So we
must take some action to avoid that \TeX~strings will only appear in the
output when a token comes along that does call~|freeze_text|. When the code
below is encountered, the fixed scrap has been copied to |*scrap_ptr|. Its
category and mathness are correct, but if |dangling_tokens()| holds, we need
to prepend the dangling tokens to its (otherwise fixed) translation. This is
done by appending that translation as a token to the dangling tokens (this is
what |app_trans|, defined below, does), then wrapping them up together, and
replacing the original translation by a pointer to the combined result
(actually we insert the pointer before consolidating the text it points to).

@< Possibly scoop up... @>=

{@; if (dangling_tokens())
  {@; app_trans(scrap_ptr); scrap_ptr->trans=text_ptr; freeze_text(); }
}

@ It is during the conversion of identifiers to scraps that their |ilk| plays
a crucial r\^ole. Ordinary identifiers, and those with |ilk| equal to
|TeX_like| or |NULL_like|, will get the category |expression|, and are set in
math mode except in the case of |TeX_like| identifiers. If the |ilk| specifies
some reserved word on the other hand, that |ilk| becomes the category of the
identifier, determining its behaviour during parsing. The translation of these
reserved words is done using |res_flag| rather than |id_flag|, as a result of
which they will be printed in boldface; they get mathness |maybe_math|,
indicating that they can be set equally well inside and outside math mode.
Identifiers whose |ilk| is |type_defined|, |const_like|, or |typedef_like|
will become reserved words with category |int_like|, and |namespace_like| of
|typename_like| becomes |struct_like|.
@:ilk change@>
These special |ilk| values have served their purpose during the scanning of
typedef declarations; this involves subtle manoeuvres that will be explained
later.

@< Append an identifier scrap @>=
{ id_pointer p=cur_id; int cat=p->ilk;
  @< Track identifiers relevant to typedef;
     maybe change |cat| from |int_like| to |expression| @>
  if (cat==normal || cat==TeX_like || cat==NULL_like)
  { app(id_flag+id_index(p));
    pack_scrap(expression
              , cat==TeX_like && !compatibility_mode ? no_math : yes_math);
  }
  else
  { if (cat==and_like || cat==not_like) /* provide text operators with space */
    { if (cat==and_like) app('~');
      app(res_flag+id_index(p)); app('~');
    }
    else app(res_flag+id_index(p)); /* append reserved word */
@)
    if (cat==type_defined || cat==const_like || cat==typedef_like)
      cat=int_like;
    else if (cat==and_like) cat=binop;
    else if (cat==not_like) cat=unop;
    else if (cat==namespace_like || cat==typename_like) cat=struct_like;
    pack_scrap(cat,maybe_math);
  }
}

@ For bad module names (e.g., an ambiguous prefix) an error has already
been reported, and they are silently suppressed from the output.

@< Append a module name scrap @>=
{@; if (cur_mod!=NULL)
    app(mod_flag+mod_index(cur_mod)), pack_scrap(mod_scrap,yes_math);
}

@ We tested in Phase~I that `\.\#' is followed by an identifier or by a
newline (which will cause |get_next| to return |end_preproc|). In the former
case we incorporate the identifier into the scrap for the preprocessor
directive; in the latter case the |lproc| scrap will just contain the
`\.\#', but since we already scanned the following |end_preproc|, we append
the corresponding |rproc| scrap as well.

@< Append a scrap starting a preprocessing directive @>=
{ app(force); app(flush_left); app_str("\\&\\#");
  if (shift()==identifier)
  {@; app(res_flag+id_index(cur_id)); pack_scrap(lproc,no_math); }
  else if (next_control==end_preproc)
  @/{@; pack_scrap(lproc,no_math);
    check_scrap(); *scrap_ptr++=token_trans[end_preproc];
  }
  else confusion("no identifier after `#'");
		@.no identifier after `\#'@>

}

@ When the `\.\v' that introduces \Cee\ text is sensed, a call on
|C_translate| will return a pointer to the \TeX\ translation of that
text.  If scraps exist in |scrap_info|, they are unaffected by this
translation process, which is useful since we convert comments to single
scraps with help of |C_translate| while building the scrap sequence for the
surrounding piece of \Cee~text.

@c
text_pointer C_translate(void)
{ text_pointer p; scrap_pointer save_base=scrap_base;
  scrap_base=scrap_ptr;
  C_read(true); /* get the scraps together */
  if (next_control != '|') err_print("! Missing `|' after C text");
				      @.Missing `|'...@>
  p=translate(); /* make the translation */
#ifdef STAT
  if (scrap_ptr>max_scr_ptr) max_scr_ptr=scrap_ptr;
#endif
  scrap_ptr=scrap_base; scrap_base=save_base; return p;
}

@ The function |outer_read| is to |C_read| as |outer_xref| is to |C_xref|:
it constructs a sequence of scraps for \Cee~text until
|next_control>=format|, taking care of embedded comments. It is called
between each occurrence of \:d, \:f, \:h, or a module name, which makes it a
convenient place to test whether the memory arrays have enough spare room to
cater for the stuff that could be needed for processing the next such token;
the most demanding requirements are for a module name heading the \Cee~part
of a section. These tests must be made outside the main loop of |outer_read|.

We use that if |next_control==end_comment| when |C_read| is called,
this value is effectively ignored.

@c
void outer_read (void) /* makes scraps from \Cee\ tokens and comments */
{ while (next_control<format)
    if (next_control!=begin_comment) C_read(false);
    else @< Read a comment, and convert it into a scrap @>
  check_scrap(); check_toks(11); /* `\.{\$\\4$m$\\PE\{\}\$$f$}' */
}

@ Since the call on |C_translate| used to process `\pb' inside comments will
itself create tokens, the token sequence for the comment under construction
must be wrapped up each time before such a call is made; the resulting token
referring to that initial segment will be contributed as first new item after
the translation is made.  Tests on the availability of two more tokens and a
text must be made inside the loop that incorporates successive `\pb'
fragments; the space for tokens outside these fragments is tested within
|scan_comment|.

The code below used to start a comment by emitting |cancel|, making each
comment go on the same line as the statement (or other code) that precedes it.
This is usually what is wanted, even if in the source file the comment starts
on a new line, but this situation made it all but impossible to make the
comment go on a line of its own in the rare cases that this is desired; for
instance putting \:/ or \:) before the comment made no difference, and also
two successive comments would force themselves onto the same line. The logic
was therefore modified to omit the |cancel| if the preceding scrap is an
|insert|, which takes care of the mentioned cases. However, we then realised
that |insert| scraps are the {\it only\/} way in which spacing can occur
before a comment, since |insert| scraps are fused with the preceding
non-|insert| scrap before the parser gets a chance at inserting formatting
controls; such controls therefore always go {\it after\/} comments. Therefore
the modified code removed the |cancel| in the only case where it could
possibly have made any difference, so we might as well omit the attempt
altogether. We have left the code but instructed the preprocessor to exclude
it, just to document this piece of the history of our \.{\me.} program.

@< Read a comment, and convert it into a scrap @>=
{ boolean one_liner=loc[-1]=='/'; int bal=0; /* brace level in comment */
  typedef_tracking(false);
  check_scrap(); check_toks(4);
#if 0
  if (scrap_ptr==scrap_base || scrap_ptr[-1].cat!=insert)
    app(cancel);
#endif
  app_str(one_liner ? "\\SHC{" : "\\C{"); @.\\C@> @.\\SHC@>
  while ((next_control=scan_comment(&bal,one_liner))=='|')
  { text_pointer p=text_ptr, q=(freeze_text(), C_translate());
    check_toks(7);
    app_tok(text_flag+text_index(p)); /* initial text */
    if (compatibility_mode) app_str("\\PB{"); @.\\PB@>
    app(inner_text_flag+text_index(q)); /* text from `\pb' */
    if (compatibility_mode) app('}');
    check_text();
  }
  app_char_tok('}'); app(force); pack_scrap(insert, no_math);
	    /* the full comment becomes a scrap */
  typedef_tracking(true);
}

@ The function |do_C| does the scanning, translation, and output of
\Cee~text within `\pb' brackets. It is called during the scanning of the
\TeX~part of a section and during the output of module names. As we have
seen, this function is not called when processing comments, where
|C_translate| is used instead, because no direct output should be produced
at such times.

@c
void do_C (void) /* read, translate, and output \Cee~text in `\pb' */
{ enter_block(1);
  if (compatibility_mode) out_str("\\PB{"); @.\\PB@>
  make_output(C_translate(),inner); /* output the list */
  if (compatibility_mode) out('}');
#ifdef STAT
  if (text_ptr>max_text_ptr) max_text_ptr = text_ptr;
  if (tok_ptr>max_tok_ptr) max_tok_ptr = tok_ptr;
#endif
  leave_block(1); /* forget the tokens */
}

@ The function |finish_C| outputs the translation of the current scraps,
preceded by the control sequence `\.{\\B}' and followed by the control
sequence `\.{\\par}'. It also restores the token and scrap memories to their
state as immediately after initialisation.

@c
void finish_C (void)
{ out_str ("\\B"); @.\\B@>
  make_output(translate(),outer);
  out_str("\\par"); finish_line();
#ifdef STAT
  if (text_ptr>max_text_ptr) max_text_ptr=text_ptr;
  if (tok_ptr>max_tok_ptr) max_tok_ptr=tok_ptr;
  if (scrap_ptr>max_scr_ptr) max_scr_ptr=scrap_ptr;
#endif
  leave_block(0); scrap_ptr=scrap_info;
	/* forget the tokens and the scraps */
}

@ The function |footnote| gives \xr. information about
further definitions of a module name (if |flag==def_flag|),
about citations of a module name (if |flag==cite_flag|),
or about the uses of a module name (if |flag==0|).
It assumes that |*p| points to the first \xr. entry
of interest, and it leaves |*p| pointing to the first element
not printed (possibly the sentinel with |num==0|).
Typical outputs are `\hbox{\.{\\Q 2001.}}',
`\hbox{\.{\\Us 370\\ET1009.}}'\ETs `\hbox{\.{\\As 8, 27\\*, 51\\ETs64.}}'.

@c
void footnote (xref_pointer* p,sixteen_bits flag)
{ if ((*p)->num<=flag) return;
  finish_line(); out('\\');
  out(flag==0 ? 'U' : flag==cite_flag ? 'Q' : 'A'); @.\\A@> @.\\Q@> @.\\U@>
  *p=list_refs(*p,flag);
  out('.');
}

@~The function |list_refs|, which does the main work for |footnote| and is
also used to produce the text replacing \:\#, distinguishes three cases,
according as the number of relevant \xr.s is one, two, or more than two. The
function always produces at least one \xr.: it is never called with
|x->num<=flag|. The value of |x| after traversing the references is
returned, for the benefit of |footnote|.

@c
xref_pointer list_refs (xref_pointer x,sixteen_bits flag)
{ xref_pointer q=next_xref(x); /* second element in \xr. list */
  if (q->num>flag) out('s');
     /* use `\.{\\As}', `\.{\\Qs}' or `\.{\\Us}' */
     @.\\As@> @.\\Qs@> @.\\Us@>
  out(' ');
  while (out_sec_nr(x->num&num_mask),x=next_xref(x),x->num>flag)
    if (next_xref(x)->num>flag) out_str(", "); /* |x| is not the last */
    else
    { out_str("\\ET"); /* next number printed will be the last */
      if (x!=q) out('s'); /* `\.{\\ETs}' for the last of more than two */
    } @.\\ET@> @.\\ETs@>
  return x;
}



@i parser.w
@i rules.w


@* Output of tokens. @:title@>
Now that we have treated the highest level of processing by \.{\me.}, we
shall have to descend again to the level of character strings, which
eventually have to written to the output file. Our first concern is to
linearise the multi-layered token lists into a sequence of output tokens.
The output of special layout tokens is affected by whether or not they were
generated from within `\pb', so there are two modes of output: during output
of `\pb' we are in |inner| mode, and otherwise in |outer| mode.
A switch from |outer| to |inner| mode can occur in the middle of a text,
namely for `\pb' fragments inside comments; this is indicated by the fact
that the reference to the translation of the fragment is tagged with
|inner_text_flag|. Apart from the fact that the mode is set to |inner|
during the output of such subtrees, no traces of the tree structure of the
internal representation of texts are left after linearisation.

Linearising texts is therefore a straightforward process, that is easy to
implement using a stack. In the linearised stream of tokens certain
sequences of tokens, particularly layout tokens, have to be considered
together, e.g., |cancel| will remove any adjacent line-breaking tokens. No
line-breaking tokens will occur at either end of the output of a complete
text either; if a line break is required there, the caller of the output
routines will supply it explicitly. For identifiers, reserved words and
module names, the transformation of a single token into a sequence of
characters is handled by separate functions to be discussed later. For
module names the transformation may result in a recursive call of the output
routines via the processing of `\pb' fragments by~|do_C|, but this recursion
in never more than one level deep, since such fragments cannot contain
module names.
@^recursion@>

@< Prototypes @>=
void out_identifier (id_pointer);
void out_keyword (id_pointer);
xref_pointer out_module_name(mod_pointer);

@ The stack that is used to keep track of token lists at different levels
of output is similar to the one used in |CTANGLE|. Entries have three parts:
for the token list at each level |tok_field| and |end_field| record where we
are respectively where we should stop, and |mode_field| records the output
mode. The current values of these quantities are referred to quite
frequently, so they are stored in a separate place, and are called
|cur_tok|, |cur_end|, and |cur_mode|.

@d cur_tok cur_state.tok_field /* location of next output token in |tok_mem| */
@d cur_end cur_state.end_field /* current ending location in |tok_mem| */
@d cur_mode cur_state.mode_field /* current mode of interpretation */

@<Typedef...@>=
typedef enum { inner, outer } mode;
typedef struct
{ token_pointer tok_field; /* present location within token list */
  token_pointer end_field; /* ending location of token list */
  mode mode_field; /* interpretation of control tokens */
} output_stack_element, *stack_pointer;

@ The stack grows upwards with the global variable |stack_ptr| pointing to
the first vacant location above the top of the stack. The entry |stack[0]|
at the bottom of the stack is not used; when |stack_ptr| points to it then
|cur_state| itself has become invalid and the output process is completed.
Therefore we can use |stack[0]| to store |cur_state| in.

@d cur_state stack[0] /* the currently active state variables */
@d stack_end (&stack[stack_size]) /* end of |stack| */

@<Global...@>=
output_stack_element stack[stack_size]; /* info for non-current levels */
stack_pointer stack_ptr=&stack[0];
  /* first unused location in the output state stack */
#ifdef STAT
stack_pointer max_stack_ptr = stack; /* largest value assumed by |stack_ptr| */
#endif

@ To insert token list |p| into the output, the function |push_level| is
called; it saves the old level of output and gets a new one going.
Conversely, the macro |pop_level| restores the conditions that were in force
when the current level was begun. The value of |cur_mode| is not changed by
|push_level|, but it might be changed explicitly directly after the call; if
so this setting will remain in effect until the matching invocation of
|pop_level|. At the beginning of |make_output|, a call to |push_level| is
made to put the root text into the output stream. If this is not a recursive
call to |make_output|, the old (undefined) value of |cur_state| will be
``saved'' into |stack[0]|, which is |cur_state| itself; although this is
redundant, no harm is done.

@d pop_level()  cur_state = *--stack_ptr
@c
void push_level (text_pointer p) /* suspends the current level */
{ if (stack_ptr==stack_end) overflow("stack"); @.stack capacity exceeded@>
  *stack_ptr++=cur_state;
#ifdef STAT
  if (stack_ptr>max_stack_ptr) max_stack_ptr=stack_ptr;
#endif
  cur_tok=text_begin(p); cur_end=text_end(p);
}

@ The function |make_output| traverses the nested token lists using the
stack, and producing the corresponding output. Since it can be called
recursively, the initial value of |stack_ptr| is saved in |stack_bot|; we
return from |make_output| when |stack_ptr| drops back again to this level.
When we encounter tokens marked with |id_flag|, |res_flag| and |mod_flag|, we
respectively call |out_identifier|, |out_keyword| and |out_module_name|.
As mentioned above, the last of these may result in a recursive call to
|make_output| via |do_C|; because calls of |do_C| also involve
lexical scanning, the values of |next_control|, |cur_id| are saved, and
restored when |make_output| is complete.
@^recursion@>

@c
void make_output(text_pointer t,mode m) /* output a complete text */
{ int save_next_control=next_control; id_pointer save_cur_id=cur_id;
  stack_pointer stack_bot=stack_ptr; token state=cancel;
  push_level(t); cur_mode=m;
  do
    if (cur_tok==cur_end) pop_level();
    else
    { token a= *cur_tok % id_flag;
      switch (*cur_tok++/id_flag)
      {
      @\@< Cases 1, 2, 3: identifiers, reserved words and module names @>
        case 4: push_level(text_at(a)); break;
	case 5: push_level(text_at(a)); cur_mode=inner; break;
	case 0: @< Output the character or format control |a| @>
      }
    }
  while(stack_ptr>stack_bot);
  @< Complete any unfinished work at the end of the translation @>
  next_control=save_next_control; cur_id=save_cur_id;
}

@ To handle the proper interaction of adjacent output tokens, we keep track
of an integer |state| variable to record that a possibly unfinished sequence
is being processed. The state is mainly affected by the output of format
control tokens like |break_space|, so for convenience we often set the state
equal to the value of such a token. If |state==0| the most recent token was
an ordinary one, and no special action is called for. If |state| has one of
the values |space|, |tilde|, |break_space|, |force|, |big_force|, |backup|,
or |big_backup|, it is the upper bound in the sense mentioned before of a
sequence of such tokens that has recently passed (where |space|~and~|tilde|
represent the tokens |' '|~and~|'~'|, respectively). We call such states
`white-space states'; when an ordinary character comes along in a
white-space state, the token recorded in |state| is first output.

@< Output white space if specified by |state| @>=
{ if (state>=space)
    if (state<break_space) out(state==space ? ' ' : '~');
    else
    { out('\\');
      if (state<backup) out(state-break_space+'5');
	/* `\.{\\5}', `\.{\\6}', or `\.{\\7}' */ @.\\5@> @.\\6@> @.\\7@>
      else out(state-backup+'6'),out_str("\\4");
        /* `\.{\\6\\4}' or `\.{\\7\\4}' */ @.\\6@> @.\\7@> @.\\4@>
      finish_line();
    }
  state=0;
}

@ The following three cases now become easy.
@< Cases 1, 2, 3... @>=
case 1: @< Output white space... @>@+ out_identifier(id_at(a)); break;
case 2: @< Output white space... @>@+ out_keyword(id_at(a)); break;
case 3: @< Output white space... @>@+ out_module_name(mod_at(a)); break;

@ Whether or not a forced break is required at the end of the translation is
determined by the function calling |make_output|, not by the code being
translated; therefore we do not invoke |@< Output white space... @>| at the
end, and if a white-space state was set, it will sink silently into
oblivion. An exception is made for |big_force| or |big_backup|, for which we
output `\.{\\Y}' so that the extra vertical white space (which must have been
inserted explicitly by a \:) control code) will not vanish, for instance
when it occurs between macro definitions.

@< Complete any unfinished work... @>=
{@; if (cur_mode==outer && (state==big_force || state==big_backup))
    out_str("\\Y");
}

@ White-space format controls operate by modifying |state| when appropriate,
as described above. The same holds basically for |' '| and |'~'|, but these
characters can also occur in ordinary text inside comments, where they must
not be contracted, so we take some care that in such a position (where
|cur_mode==outer| holds), they are output as ordinary characters.  In
white-space states |opt| and the digit following it are ignored. Like the
white-space format controls, the token |cancel| also sets |state| to its own
value, but it ignores the value it previously had. In this state any spaces
and format control tokens except |indent|, |outdent| and |flush_left| are
ignored. We also have set |state==cancel| at the beginning of |make_output|,
so that there will be no forced break or white space at the beginning of the
translation.

The tokens |indent| and |outdent| are output directly without altering
|state|; effectively this means that they will be moved in front of any
white-space format controls that preceded them. This is important because
the format controls that are explicitly inserted by codes like \:/ will
stick to the token to their left, so that any indentation changes generated
by the syntax rules at the same point will necessarily come {\it after\/}
them, while in fact these indentation changes must be allowed to affect the
indentation at explicitly inserted line breaks. With the exceptions noted
above, |opt| and |flush_left| are transmitted like normal tokens.

When |cur_mode==inner| things are a bit different. The tokens |indent|,
|outdent|, and |flush_left| are completely ignored, and all white-space
tokens except |'~'| are treated like |' '|, i.e., they set |state=space|.
Optional breaks are produced using `\.{\\0}' rather than `\.{\\3}' to avoid
a ragged right margin within a paragraph.

@< Output the character or format control |a| @>=
{ switch (a)
  { case relax: @< Output white space... @>@+ break;
    case cancel: state=cancel; break;
    case indent: case outdent:
      if (cur_mode==outer) {@; out('\\'); out(a-indent+'1'); }
      break;   @.\\1@> @.\\2@>
    case opt:
    { int digit=*cur_tok++;
      if (state==0)
      {@; out('\\'); out(cur_mode==outer ? '3' : '0'); out(digit); }
					@.\\3@> @.\\0@>
      break;
    }
    case flush_left:
      if (cur_mode==outer)
      {@; @< Output white space... @> out_str("\\8"); }  @.\\8@>
      break;
    case big_force: case backup:
      if (a+state==big_force+backup) a=big_backup; /* fall through */
    case break_space: case force: case big_backup:
      if (cur_mode==inner) a=space;
    up_state:
      if (state!=cancel && state<a) state=a;
      break;
    case ' ': case '~':
      if (cur_mode==inner) {@; a= a==' ' ? space : tilde; goto up_state; }
      if (state==cancel || state>=break_space) break;
        /* else fall through */
    default: @< Output white space... @>@+ out(a);
  }
}

@*1 Low-level output routines. @:title@>
We now come to the functions that write individual tokens and characters.
We start with a simple function to output a section number in decimal
notation. The number to be converted by |out_sec_nr| is known to be less
than |def_flag|, so it cannot have more than five decimal digits. If the
section was changed, we output `\.{\\*}' just after the number.

@c
void out_sec_nr (int n) /* output a section number */
{ char s[6];
  sprintf(s,"%d",n); out_str(s);
  if (section_changed(n)) out_str ("\\*"); @.\\*@>
}

@ Since we want to typeset identifiers such as \.{catch22} as |catch22|, which
is achieved by writing the string `\.{\$\\\\\{catch\}\_\{22\}\$}' to the
\TeX~file, it is useful to have an auxiliary function |out_id_part| that
will output a substring of an identifier (or reserved word), enclosing it in
braces if it contains more than one character. The substring is specified by
a pointer and a length count. A call to |out_id_part| for the full name
referred to by an |id_pointer p@;| is achieved by invoking |out_id_full(p)|.


@d out_id_full(p) out_id_part(name_begin(p),length(p))

@c
void out_id_part (char* s, int l)
{ boolean b=l!=1;
  if (b) out ('{');
  while (--l>=0) {@; if (*s=='_') out ('\\'); out (*s++); }
  if (b) out('}');
}

@ User specified index entries, those produced by `\.{@@\^}', `\.{@@.}', or
`\.{@@?}', are different from identifiers since they are not subject to the
restrictions of the \Cee\ programming language; if users should decide to use
any characters with a special meaning to \TeX\ here, it is their
responsibility to properly escape them. Therefore there is no reason to give
underscores a special treatment, and we shall define a function |out_index| to
use in place of |out_id_full| for those cases. Index entries of length~0 are
discarded on input, so they never come here; for index entries of length~1 we
omit braces.

In \LKC., there is no function like |out_index|, and the equivalent of
|out_id_part| is used instead, which has as a consequence that underscores are
automatically escaped in index entries, but other special characters for \TeX\
are not. While there is no particular rhyme or reason to such a behaviour,
users have learned to adapt to it, and so, in an attempt at bug-to-bug
compatibility, we transfer control to |out_id_part| in compatibility mode.

@c
void out_index (id_pointer p)
{ char* s=name_begin(p); boolean b= s[1]!='\0';
  if (compatibility_mode) {@; out_id_full(p); return; }
  if (b) out ('{');
  out_str(s);
  if (b) out('}');
}

@ The function |out_keyword| just prefixes `\.{\\\&}' to the output
of |out_id_full|.

@c
void out_keyword(id_pointer p)
{@; out_str("\\&"); out_id_full(p); }
@.\\\&@> @.\\templ@>

@ The function |out_identifier| formats ordinary identifiers (those that were
stored with |id_flag|). Yet there is some variation in the way these are
formatted. Identifiers whose |ilk| is |TeX_like| or |NULL_like| are converted
into control sequences. Otherwise the variable |kind| will be used to classify
the identifier. If only one alphabetic character is present, it is set in math
(rather than text) italic (|kind==single| or |kind==indexed_single|),
otherwise if all characters are either upper case letters or digits, the
identifier is set in typewriter type (|kind==caps|), and in the remaining,
most common, case, it is set in text italic (|kind==ord| or |kind==indexed|).
If math or text italic are used and we are not in compatibility mode, then a
trailing sequence of digits in the name, if present, is set as a subscript to
the rest of the identifier (|kind==indexed| or |kind==indexed_single|).

@c
void out_identifier (id_pointer p)
{ int k=length(p); eight_bits* ch=(eight_bits*)name_begin(p);
@/enum { ord, indexed, caps, single, indexed_single } kind;
  if (p->ilk==TeX_like || p->ilk==NULL_like) @.\\NULL@> @.\\TeX@>
  @/{@;
    @< Output the identifier |p| in the form of a \TeX\ control sequence @>
    return;
  }
  if (!compatibility_mode) /* then search for possibly trailing digits */
  { do --k; while (isdigit((eight_bits)ch[k]));
    /* terminates because |!isdigit(ch[0])| */
    ++k; /* point to end of identifier without its index (if any) */
  }
  @< Determine the |kind| of |p|, and output some leading characters,
     distinguishing the case |kind==caps| @>
  if (kind==indexed || kind==indexed_single)
  { out_id_part(name_begin(p),k); /* main part */
    out ('_'); out_id_part(name_begin(p)+k,length(p)-k); /* subscript */
  }
  else out_id_full (p);
}

@ When an identifier is output in the form of a control sequence,
underscores are replaced by `\.x' so that they will become part of the
control sequence. If the |ilk| of the identifier was |TeX_like|, the control
sequence will be processed in horizontal mode with italic type selected; if
the |ilk| was |NULL_like|, it will be processed in math mode.

@< Output the identifier |p| in the form of a \TeX\ control sequence @>=
{ if (p->ilk==TeX_like) out_str("\\\\{");
  out('\\'); @+ do out(*ch=='_' ? 'x' : *ch); while (*++ch!='\0');
  if (p->ilk==TeX_like) out('}');
}

@ For identifiers that consist of a single character, not counting any
trailing digits, that character is written in unadorned form, since it will
be processed in math mode, so that the math italic font will be used. We do
output a space before the character however, to eliminate the possibility
that it is captured by a control word at the end of the previous output.
In other (multi-letter) cases, we output an initial `\.{\\\\}', or `\.{\\.}'
for all-caps identifiers.

@< Determine the |kind|... @>=
if (k==1) {@; out(' '); kind= length(p)==1 ? single : indexed_single; }
else
{ int i=k; @+
  while (--i>=0) @+ if (!isupper(ch[i])&&!isdigit(ch[i])&&ch[i]!='_') break;
  kind= i<0 ? caps : k<length(p) ? indexed : ord;
@/out('\\'); out(kind==caps ? '.' : '\\');	@.\\\\@>  @.\\.@>
}

@ The last function for output of tokens, |out_module_name|, represents the
most complicated aspect of such output, since up to here we have not paid
any attention to `\pb' constructions in module names. When these are present
they will need a complete treatment by the scanning, parsing and general
output routines described above. Now |out_module_name| does not involve
itself directly in such actions---they are contained in a call to
|do_C|---but we do have to be aware that we are in the midst of scanning and
output operations going on at a different level, and be careful not to
disturb them. The text between `\pb' will be placed at the end of the active
input buffer (this is why the buffer has size |long_buf_size| rather than
|buf_size|) and the translation process uses the end of the active |tok_mem|
area.

On the other hand |out_module_name| is also used during Phase~III while
writing the list of module names (which implies that scanning and parsing is
not over at the end of Phase~II). We want |out_module_name| to behave
slightly differently in Phase~III however, by printing a full list of
section numbers in which the module is defined rather than just the first
one. Also, immediately after calling |out_module_name| we will need a
pointer to the sublist of \xr.s to sections in which the module is used, and
since this sublist starts at the end of the list of defining sections, we
might as well return it as result from |out_module_name|.

@c
xref_pointer out_module_name(mod_pointer name)
{ xref_pointer x=name->xref; boolean file_module= x->num==file_flag;
  if (file_module) x=next_xref(x);
  out_str ("\\X"); @.\\X@>
  if (x->num>=def_flag)
  { out_sec_nr(x->num-def_flag); /* output the defining section number */
    if (phase==3) /* all of them in Phase III */
      while (x=next_xref(x), x->num>=def_flag)
	out_str (", "), out_sec_nr(x->num-def_flag);
  }
  else out ('0'); /* section number `0' means `nowhere defined' */
  out (':'); if (file_module) out_str("\\.{"); @.\\.@>
  @<Output the text of the module name@>
  if (file_module) out_str ("}");
  out_str ("\\X");
  return x;
}

@ Copying \TeX\ text from a module name to the output is fairly simple, but
we shouldn't forget to escape special characters if this is actually a file
name (which is typeset verbatim), even though they are probably quite rare.

@< Output the text... @>=
{ char* k=name_begin(name),c;
  while ((c=*k++)!='\0')
  { if (file_module) @+{@; if (strchr(" \\#%$^{}~&_",c)!=NULL) out ('\\'); }
    if (c=='@@' && *k++!='@@') @< Report illegal control code in module name @>
    if (file_module || c!='|') out(c);
    else
    { char* save_loc=loc, *save_limit=limit;
      @< Copy the \Cee\ text into the |buffer| array @>
      do_C(); loc=save_loc; *(limit=save_limit)=' ';
    }
  }
}

@ We haven't checked for illegal control codes in module names yet, so
we should report an error if we encounter one.

@< Report illegal control... @>=
{@; print("\n! Illegal control code in module name"); print_mod(name);
	     @.Illegal control code...@>
  mark_error();
}

@ Within `\pb' we should be aware of character and string constants, since
any `\.\v' occurring there is not the closing one. The variable |delimiter|
is zero outside such constants, otherwise it equals the delimiter that began
the constant.  We copy the opening and closing `\.\v' into the buffer, so
that an error message that displays the whole buffer will look a little bit
sensible. We also add a space at the end (just like |get_line| does) so that
the closing `\.\v' cannot by accident be parsed as part of an operator such
as `\.{\v=}'. Putting |next_control='|'| completes the proper initial
conditions for calling |do_C|. We need not test for overflow of |buffer|
here, since we start filling it at position |limit| which is at most
|&buffer[buf_size-2]|, and add at most |longest_name| characters from the
module name, so that if worst comes to worst, the final space stored at the
end is at position |&buffer[long_buf_size-2]| leaving one more place after
|limit|, which is occasionally needed by |get_next|.

@< Copy the \Cee\ text into... @>=
{ char delimiter='\0'; /* |'"'| or |'\''|, or |'\0'| outside of strings */
  next_control=*limit++='|'; loc=limit;
  do
    if ((c=*k++)=='\0') @< Report a runaway \Cee~text @>
    else
    { *limit++=c;
    @/@< In special cases copy the character after |c|,
         or change |delimiter| @>
    }
  while (c!='|' || delimiter!='\0');
  *limit=' ';
}

@ Here cases like `\.{1+@@\v2}', `\.{@@'\v'}' and `\.{'\\''}' within `\pb'
in module names are handled properly, as well as the more ordinary cases of
character constants and strings. In case of a control code, the character
following the `\.@@' cannot be the null character that terminates the module
name, because it is impossible to enter a module name that ends with a
single `\.@@'. If the user is being really devious by saying `\.{@@< Sudden
\v"death\\ @@>}', then the final `\.\\' is silently removed to avoid
disaster to \.{\me.}.

@< In special cases... @>=
{ if (c=='@@') /* control code; now |*k!='\0'| */
  { if ((*limit++=*k++)=='\'' && delimiter=='\0')	/* copy code, test for \:' */
      delimiter='\''; /* which behaves like `\.'' */
  }
  else if (c=='\\' && delimiter!='\0')
  { char d=*limit++=*k++; /* escaped character, possibly |delimiter| */
    if (d=='\0') --k,limit-=2; /* remove backslash, error is issued anyway */
  }
  else if (c=='\'' || c=='"')
    if (delimiter=='\0') delimiter=c;
    else if (delimiter==c) delimiter='\0';
}

@ If the module name ends before the \Cee~text does, we add an |'|'|, and if
a string or character constant was left unclosed we prepend the closing
delimiter to it: we must prevent that |do_C| will ever call |get_line|, even
in erroneous situations. In these cases there is a slight chance that we
might overflow |buffer|, in which case we quit since the user is just trying
to break \.{\me.} anyway.

@< Report a runaway \Cee~text @>=
{ print("\n! C text in module name didn't end"); print_mod(name);
	   @.C text...didn't end@>
  mark_error();
  if (delimiter!='\0') *limit++=delimiter;
  *limit++='|';
  if (limit>&buffer[long_buf_size-2])
    fatal("fix that first, you sneaky devil");
         @.fix that first...@>
  break;
}

@ Finally we come to the point where the characters produced are sent off to
the \TeX~file. The \TeX\ output is supposed to appear in lines at most
|line_length| characters long, so we place it into an output buffer
|out_buf|. The first character of the output buffer is used as a sentinel,
and is not actually written to the output file, so we usually refer to the
output buffer via |out_line|. The pointer |out_ptr| indicates the next
position to be used in the output buffer. During the output process,
|out_line_nr| will hold the current line number of the line about to be
output; it is only used for (rather unlikely) diagnostic messages.

@d out_line (&out_buf[1]) /* start of actual output line */
@d out_buf_end  (&out_line[line_length]) /* end of |out_buf| */

@<Global...@>=
char out_buf[line_length+1]; /* assembled characters */
char* out_ptr; /* first unused position in |out_buf| */
int out_line_nr=1; /* number of next line to be output */

@ The auxiliary function |flush_buffer| empties the buffer up to a given
breakpoint~|b|, and moves any remaining characters to the beginning of the
buffer, so that they will appear on the next line. If the |percent|
parameter is |true| a |'%'| is appended to the line that is being output; in
this case the breakpoint |b| should be strictly less than |out_buf_end|. If
the |percent| parameter is |false|, trailing blanks are suppressed.  The
characters emptied from the buffer form a new line of output.

@d tex_putc(c) putc (c, tex_file)
@d tex_new_line() (putc('\n', tex_file),++out_line_nr)
@d tex_printf(format) fprintf(tex_file, format)

@c
void flush_buffer(char* b, boolean percent)
   /* output from |out_line| to |b|, where |b<=out_ptr| */
{ int j=(int)(b-out_line); /* number of characters to be output */
  if (!percent) @+ while (j>0 && out_line[j-1]==' ') --j;
    /* remove trailing blanks */
  fprintf(tex_file, "%.*s",j,out_line);
  if (percent) tex_putc('%');
  tex_new_line();
  { char* p=out_line;
    while (b<out_ptr) *p++=*b++; /* shift back remainder of line */
    out_ptr=p; /* adjust to end of (possibly empty) shifted part */
  }
}

@ The usual way to send a completed line to the output is to call
|finish_line|. When we are copying \TeX\ source material, we retain line
breaks that occur in the input, except that an empty line is not output when
the \TeX\ source line was not entirely blank. For example, a line of the
\TeX\ file that contains only an index \xr. entry will not be copied. Since
|get_line| has removed trailing blanks, the input line is blank if and only
if it is empty; similarly the emptiness of the output line is tested by a
simple comparison because |copy_limbo| and |copy_TeX| have not copied spaces
into an empty output line. The |finish_line| routine is called just before
|get_line| inputs a new line, and just after a line break token has been
emitted during the output of translated \Cee~text.

@d output_line_empty() (out_ptr==out_line)

@c
void finish_line(void) /* do this at the end of a line */
{ if (!output_line_empty()) flush_buffer(out_ptr, false);
  else if (limit==buffer) tex_new_line(); /* copy blank input line */
}

@ In particular, the function |finish_line| is called near the very
beginning of Phase~II.  We initialize the output buffer in such a way that
the first line of the output file will be `\.{\\input cwebxmac}', or
`\.{\\input cwebcmac}' in compatibility mode. At position~0 of the output
buffer, which is never really output, we put a space that will help make the
function |break_out| a little faster.

@< Set initial... @>=
{ char* line1=" \\input cwebxmac";
  out_ptr=&out_buf[0]; @+ do *out_ptr++=*line1++; while (*line1!='\0');
  if (compatibility_mode) out_ptr[-4]='c'; /* change to \.{cwebcmac} */
}

@ When we wish to append one character |c| to the output buffer, we write
`|out(c)|'; this will cause the buffer to be broken at a sensible place and
flushed, if it was already full.  If we want to append more than one
character at once, we say |out_str(s)|, where |s| is a string containing the
characters.

@d out(c)
  *(out_ptr>=out_buf_end ? (break_out(),out_ptr++) : out_ptr++)=c

@c
void out_str (char* s) @+{@; while (*s!='\0') out (*s++); }

@ The function |break_out| will determine a good break point in the output
buffer when it is about to overflow.

@< Prototypes @>=
void break_out (void);

@~For speed we search from right to left for a proper break point, although
in this direction we cannot know whether we are inside a control sequence.
Nevertheless any blank space is a safe break point, as is the position just
before a backslash that isn't preceded by another backslash. If the break is
not at a space, a |'%'| is output at the break.

@c
void break_out (void) /* finds a way to break the output line */
{ char* k=out_ptr,c; int count=0; /* number of backslashes seen */
  do @+ if ((c=*--k)==' ') goto found; @+ while (c!='\\');
  do ++count; while ((c=*--k)=='\\');
found:
  if (++k>out_line) flush_buffer(k,c!=' '); @+ else @< Break peculiar line @>
}

@ We get to this section only in the unusual case that the entire output line
consists of a string of backslashes followed by a string of characters that
contains no spaces or backslashes. Depending on the number |count| of
backslashes we handle this case as well as possible. If the |count>=2| we
split off a maximal even number of initial backslashes, and if |count==0| we
break the line by putting a |'%'| at the position of the last character. In
the latter case we do not place the |'%'| after the last character, since
that would make the output line one character longer than we promised;
similarly we are cautious in the former case not to place a |'%'| after a
buffer completely filled with backslashes. If |count==1| the line is
probably one enormous control word, and in this unlikely case we simply
output the whole line as is stands, after issuing a warning message.

@< Break peculiar line @>=
if (count==0) flush_buffer(out_ptr-1,true);
else if (count>=2) flush_buffer
  (&out_line[count&=~1]==out_buf_end ? out_buf_end-2 : &out_line[count],true);
else
{ print("\n! Line had to be broken (output l.%d):\n",out_line_nr);
	   @.Line had to be broken@>
  term_write(out_line,out_ptr-out_line); new_line(); mark_harmless();
  flush_buffer(out_ptr,false);
}


@* Phase III processing. @:title@>
After Phase~II is completed, \.{\me.}'s only remaining task is writing out
the index, after sorting the identifiers and index entries, and a list of
module names.  If the user has set the |no_xref| flag (the \.{-x} option on
the command line), we just finish off the page, omitting the index, module
name list, and table of contents.

@d triple_file_output flags['t']
@d even_out_pages flags['e']

@c
void phase_three (void) /* output the \xr. index */
{ finish_line();
  if (no_xref) out_str("\\endcodemode\\vfill\\end");
  else
  { phase=3; print_progress("\nWriting the index...");
			     @.Writing the index...@>
    typedef_tracking(false); /* during parse of `\pb' in module names */
    if (change_exists)
    {@; @<Tell about changed sections @> finish_line(); }
    if (triple_file_output)
      @< Finish the \TeX~file and open the index file @>
    else {@; out_str("\\inx"); finish_line(); }  @.\\inx@> /* index */
    @<Sort and output the index@>
    if (triple_file_output)
      @< Finish the index file and open the module listing file @>
    else {@; out_str("\\fin"), finish_line(); } @.\\fin@> /* end of index */
    @<Output all the module names@>
    if (!triple_file_output)
    { out_str("\\con"); @.\\con@> /* table of contents */
      if (even_out_pages) out_str("even"); @.\\coneven@>
    }
  }
  finish_line(); fclose(tex_file);
  print_progress("\nDone.\n");
  check_complete(); /* was all of the change file used? */
}

@ When output is to be distributed over three files, we can already finish
off the main output file, since it ends in a standard way: with `\.{\\input}'
commands for the two other files (that still have to be written) and either
`\.{\\con}' or `\.{\\coneven}' to produce the table of contents. Since the
three output files are written one after the other, we can reuse the file
pointer |tex_file| for all of them.

@< Finish the \TeX~file and open the index file @>=
{ out_str("\\inx \\input \\jobname.idx\n" @+
	  "\\fin \\input \\jobname.scn\n" @+
	  "\\con");
  if (even_out_pages) out_str("even");
  finish_line(); fclose(tex_file);
  if ((tex_file=fopen(idx_file_name,"w"))==NULL)
    fatal("! Cannot open \"%s\" as output file",idx_file_name);
}

@~The switch from the index file to the module listing file is similar, but
simpler.

@< Finish the index file and open the module listing file @>=
{ finish_line(); fclose(tex_file);
  if ((tex_file=fopen(scn_file_name,"w"))==NULL)
    fatal("! Cannot open \"%s\" as output file",scn_file_name);
}

@ Just before the index comes a list of all the changed sections, including
the index section itself. The outer |while| loop will terminate without
overstepping the bounds of |changed_section| in the inner loop since
|section_changed(section_count)| holds.

@<Tell about changed sections@>=
{ int k=0; boolean first=true;
  out_str("\\ch "); @.\\ch@> /* changes */
  while (k<section_count)
  { do ++k; while (!section_changed(k));
    if (first) first=false; @+ else out_str(", ");
    out_sec_nr(k);
  }
  out('.');
}

@ A left-to-right radix sorting method is used, since this makes it easy to
adjust the collating sequence, and the running time will be at worst
proportional to the total length of all entries in the index.
During the sorting phase the partially sorted index entries that have not
yet been output are held on a stack of identifier lists called |sort_info|.
Each of these lists is unsorted, but between elements of different lists the
ordering is completely known: for any pair of lists, all elements of the
list nearer to the top of the stack precede all elements of the other list
(this ordering is chosen because elements are output at the top of the
stack). The entries of |sort_info| are structures of type |sort_node|; they
have a |head| field that points to the first name of the list, and a
|depth| field that tells how many initial characters are known to be equal
among the entries of the list.

@<Typedef...@>=
typedef struct {@; id_pointer head; int depth; } sort_node, * sort_pointer;

@~In fact |sort_info| shares its memory with the |scrap_info| array that
serves no purpose at this moment (although it will be used again during
the output of module names).

@d sort_info scrap_union.id_list_field
@d sort_info_end (&sort_info[sort_stack_size])
@< Alternative use... @>=
sort_node id_list_field[sort_stack_size];

@ The variable |sort_ptr| points to the first unused slot above the top of
the stack in~|sort_info|.

@<Global...@>=
sort_pointer sort_ptr=sort_info;
#ifdef STAT
sort_pointer max_sort_ptr=sort_info; /* largest value of |sort_ptr| */
#endif

@ The desired alphabetic order is specified by the |collate| array; namely,
|collate[0] < collate[1] < @t\dots@> < collate[end_collate-1]|. Upper case
letters are treated like the corresponding lower case letters, since we want
to have `|t<TeX<@[token@]|'. Therefore the collation array should have length
|UCHAR_MAX+1-26|. We are a bit extra cautious however, in case there are
implementations in which the verdict |isalpha(c)| is given (incorrectly) to
certain values |c>=128|, which would cause the part of |collate| that is
actually used to be shorter than that; therefore we compute the actual
length in |end_collate|.

@<Global...@>=
eight_bits collate[UCHAR_MAX-25]; /* collation order */
int end_collate; /* length of the part of |collate| actually used */

@ We use the order |'\0' < ' ' < @tother characters@> < '_' < 'A'=='a' <
@t\dots@> < 'Z'=='z' < '0' <@t\dots@> < '9'|. If there should be any
characters~|c| for which |tolower(c)| does not occur in |collate| (this can
only happen if |isalnum(c)| incorrectly holds), then entries containing the
character~|c| will not appear in the index. By computing |end_collate|
rather than using a constant |UCHAR_MAX-25|, we avoid however that the mere
existence of such characters would disrupt the proper order of the
index entries that do appear, due to spurious bytes |'\0'| at the end
of~|collate|.

@<Set init...@>=
{ char *p="_abcdefghijklmnopqrstuvwxyz0123456789";
  int c='\1', k=2;
  collate[0]='\0'; collate[1]=' ';
  do @+ if (!isalnum(c) && c!='_' && c!=' ') collate[k++]=c; @+
  while (++c<=UCHAR_MAX);
  while ((c=*p++)!='\0') collate[k++]=c;
  end_collate=k; /* record the length */
}

@ The lists of identifiers used for sorting the index cannot be linked
together using |hash_link|, since we still need to be able to look up
identifiers when writing the list of module names. Therefore we declare a
separate array |index_link| that will provide the necessary links, and a
macro |ilink| used to access it: |ilink(p)| points to the successor of~|p|
(or is~|NULL|) for any |id_pointer p@;|.

At each sorting step we partition a list of identifiers into at most
|UCHAR_MAX-25| sublists, based on the first character position where the
entries of the list are not known to be equal. After the partitioning step
the sublist for character~|c| is pointed to by |bucket[tolower(c)]|.

@d ilink(p) index_link[id_index(p)]

@<Global...@>=
id_pointer index_link[max_idents]; /* links identifiers during sorting */
id_pointer bucket[UCHAR_MAX+1];

@ The basic sorting step is to pop the list from the top of the stack, and
either output it, if needs no further refinement, or to split it up into the
buckets, after which the function |unbucket| collects the non-empty buckets
and pushes them back on the stack. At the very first step however, the
identifiers do not come from the stack but directly from the identifier
table. In the remaining cases, the |depth| field of the element popped off
the stack is used to determine the character position used for splitting,
and is passed to |unbucket| so that it can set the proper depth when it
pushes lists back on the stack. When |depth| has risen to~|255| (or more
likely, has been set explicitly to that value by a previous |unbucket|), no
attempt is made to split up a list, even if it is not reduced to a
singleton.

@d infinity 255 /* $\infty$ (approximately) */

@< Sort and output... @>=
{ @<Do the first pass of sorting@>
   /* the first time, entries do not come from |sort_info| */
  unbucket(1); /* pick up first-order bucketed lists */
  while (sort_ptr>sort_info) /* i.e., the stack is not empty */
  { eight_bits depth=(--sort_ptr)->depth;
    id_pointer name=sort_ptr->head;
    if (ilink(name)==NULL || depth==infinity)
        /* singleton or set of look-alikes */
      @< Output index entries for the list starting at |name| @>
    else
    {@; @< Split the list starting at |name| into further lists @>
      unbucket(depth+1);
    }
  }
}

@ To begin the sorting, we go through all the hash lists and put each entry
having a non-empty \xr. list into the proper bucket. The buckets are emptied
initially; they will also be emptied by each call to~|unbucket|. The entries
of the |index_link| array are initialised when their corresponding
identifiers are prepended to the list of some bucket.

@<Do the first pass...@>=
{ id_pointer name;
  eight_bits c=UCHAR_MAX;
  id_pointer *h; /* pointer into |hash| */

  do bucket[c]=NULL; while (c--!=0);
  for (h=hash; h<hash_end; h++)
    for (name=*h; name!=NULL; name=name->hash_link)
	/* traverse all hash lists */
      if (name->ilk!=reference && name->ilk!=header_file_name @|
         && name->xref->num!=0)
      	/* leave out non-identifiers and unreferenced names */
      {@; c=name_begin(name)[0]; c=tolower(c);
	ilink(name)=bucket[c]; bucket[c]=name;
      }
}

@ Here we split the list at the top of the stack into buckets. Entries of
length |depth| will bet sent to |bucket['\0']| due to the null character
stored at the end of identifier names; shorter entries do not come here
since they will already have been output. Since we are changing the |ilink|
of the front node of the list when putting it into its bucket, we must be a
bit more cautious then usual in traversing the list.

@< Split the list... @>=
do
{ eight_bits c= tolower((eight_bits)name_begin(name)[depth]);
  id_pointer next_name=ilink(name); /* save link */
  ilink(name)=bucket[c]; bucket[c]=name; name=next_name;
    /* put into bucket */
}
while (name!=NULL);

@ The function |unbucket| goes through the buckets in the reverse order of
the collating sequence specified in the |collate| array, and adds non-empty
lists to the stack. The parameter |d| to |unbucket| tells the current depth
in the buckets; it will be recorded in all lists pushed on the stack during
this call, except that the contents of the bucket containing the identifiers
of length~|d-1| will be given depth |infinity| so that they will be output
directly afterwards. Any two sequences that agree up to the case of their
characters (with any characters beyond position~255 being ignored) are
regarded as identical and will be output in a random order. It is only
because of this case-merging that there could possibly be more than one
identifier of length~|d-1|, which makes setting the depth to |infinity|
necessary; in other cases the special attention is superfluous, since
singleton lists will be output directly anyway.

@< Prototypes @>=
void unbucket (eight_bits);

@~@c
void unbucket (eight_bits d) /* empties buckets having depth |d| */
{ int i=end_collate; /* index into |collate| */
  while(--i>=0)	@+ if (bucket[collate[i]]!=NULL)
  { if (sort_ptr>=sort_info_end)
      overflow("sorting"); @.sorting capacity exceeded@>
    sort_ptr->depth= i==0 ? infinity : d;
      /* |infinity| means there is nothing left to compare */
    sort_ptr++->head=bucket[collate[i]];
    bucket[collate[i]]=NULL; /* push and empty bucket */
#ifdef STAT
    if (sort_ptr>max_sort_ptr) max_sort_ptr=sort_ptr;
#endif
  }
}

@ Each line in the index file starts with the macro `\.{\\@@}', which will
be properly defined when it is encountered.

@< Output index... @>=
do
{ out_str("\\@@");   @.\\@@@>
  @<Output the index entry |name|@>
  @<Output the \xr.s of |name|@>
} while ((name=ilink(name))!= NULL);

@ In this section the distinction between |xref_roman|, |xref_wildcard|, and
|xref_typewriter| finally becomes effective. As first character (immediately
after `\.{\\@@}') we output a character `\.h' or `\.m',indicating whether the
index entry should be processed in horizontal or in math mode.
@.\\NULL@> @.\\TeX@>

@<Output the index entry...@>=
switch (name->ilk)
{ case normal: case NULL_like: out('m'); out_identifier(name); break;
  case TeX_like: out('h'); out_identifier(name); break;
  case roman: out('h'); out_index(name); break;
  case wildcard: out_str("h\\9"); out_index(name); break;   @.\\9@>
  case typewriter: out_str("h\\."); out_index(name); break;   @.\\.@>
  default: out('h'); out_keyword(name);
}

@ Section numbers that are to be underlined are enclosed in
`\.{\\[}\dots\.]'. The first `\.{,\ }' will be scooped up by the macro
`\.{\\@@}'.

@<Output the \xr.s...@>=
{ xref_pointer x=name->xref;
  do
  { sixteen_bits n=x->num;
    out_str(", ");
    if (n<def_flag) out_sec_nr(n);
    else {@; out_str("\\["); out_sec_nr(n-def_flag); out(']'); }   @.\\[@>
  } while ((x=next_xref(x))->num!=0);
  out('.'); finish_line();
}

@ The following recursive function traverses the tree of module names and
prints them. We use the order already present in the tree, which means that
we do not use the collation sequence that was employed for sorting the
index. If the user starts all module names with capital letters however, the
difference should hardly be noticeable.
@^recursion@>

@< Prototypes @>=
void list_modules (mod_pointer);

@~The macro `\.{\\@@}' is redefined to serve for lines in the list of module
names as well. Since module names are to be processed in math mode, we
enclose them in dollar signs.

@c
void list_modules (mod_pointer p) /* print all module names in subtree |p| */
{ if (p != NULL)
  { list_modules(p->llink);
@/  out_str("\\@@$");   @.\\@@@>
    leave_block(0); scrap_ptr=scrap_info; /* get ready for parsing */
    {@; xref_pointer x=out_module_name(p); out('$');
      footnote(&x,cite_flag); footnote(&x,0);
    }@/
    finish_line();
@/  list_modules(p->rlink);
  }
}

@~Initially |list_modules| is called for the root of the tree, of course.

@<Output all the module names@>= list_modules(root);

@ At the end of the run, if |STAT| was defined and the `\.{+s}' flag
present, we report how much of all the arrays was actually needed.

@d report(k,c,m)
  printf("\t%lu %ss (out of %lu)\n",(unsigned long)(c),k,(unsigned long)(m))
@c
#ifdef STAT
void print_stats()
{ print("\nMemory usage statistics:\n"); @.Memory usage statistics@>
@/report("identifier",		id_index(id_ptr),	max_idents);
@/report("module name",		mod_index(mod_ptr),	max_modules);
@/report("byte",		byte_ptr-byte_mem,	max_bytes);
@/report("cross-reference",	xref_ptr-xmem,		max_refs-1);
@)printf("Parsing:\n");
@/report("scrap",		max_scr_ptr-scrap_info,	max_scraps);
@/report("text",		max_text_ptr-text_mem,	max_texts);
@/report("token",		max_tok_ptr-tok_mem,	max_toks);
@/report("trie node",		node_no,		max_no_of_nodes);
@/report("level",		max_stack_ptr-stack,	stack_size);
@)printf("Sorting:\n");
@/report("level",		max_sort_ptr-sort_info,	sort_stack_size);
}
#endif


@* Index. @:index@> @:title@>
If you have read and understood the code for Phase~III above,
you know what is in this index and how it got here.  All sections in
which an identifier is used are listed with that identifier, except
that reserved words are indexed only when they appear in format
definitions, and the appearances of identifiers in module names are
not indexed. Underlined entries correspond to where the identifier
was declared. Error messages, control sequences put into the output,
and a few other things like ``recursion'' are indexed here too.
