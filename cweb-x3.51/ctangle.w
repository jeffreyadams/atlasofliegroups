% This file is part of CWEBx.
% This program by Marc van Leeuwen based on earlier versions by
% D. E. Knuth., Silvio Levy and Frank Jensen.
% It is distributed WITHOUT ANY WARRANTY, express or implied.
% CWEB (Revision: 2.0) % Don Knuth, July 1990
% Version 3.x, Marc van Leeuwen, December 1993
% CWEBx 2+1.0, Marc van Leeuwen, August 1994
% CWEBx 3.0, Marc van Leeuwen, Januari 1995
% CWEBx 3.02, Marc van Leeuwen, April 1996
% CWEBx 3.02a, Marc van Leeuwen, September 1996
% CWEBx 3.03, Marc van Leeuwen, January 1998 (ctangle unchanged)
% CWEBx 3.1, Marc van Leeuwen, May 2006
% CWEBx 3.5, Marc van Leeuwen, June 2009

% Copyright (C) 1987,1990 Silvio Levy and Donald E. Knuth
% Copyright 1994--2006 Marc A. A. van Leeuwen

% Permission is granted to make and distribute verbatim copies of this
% document provided that the copyright notice and this permission notice
% are preserved on all copies.

% Permission is granted to copy and distribute modified versions of this
% document under the conditions for verbatim copying, provided that the
% entire resulting derived work is distributed under the terms of a
% permission notice identical to this one.

\def\me.{CTANGLE} \def\myroots{\.{TANGLE}}

% Here is TeX material that gets inserted after \input cwebxmac
\def\ASCII.{\caps{ASCII}}

@i intro.inc  % Here is some source text that matches the start of CWEAVE

@d banner "This is CTANGLE (version "@+version_string@+")"
@q the C macro |version_string| was defined at the end of intro.inc@>

@ The following parameters are specific to \.{CTANGLE}; those which are
common to \.{CTANGLE} and \.{CWEAVE} are defined in the file \.{common.inc}
and appear below. Some of these values have been decreased with respect to
their earlier values which were sufficient in the original \.{WEB} to handle
\TeX; a motivation is given at the common declarations. The macro |variant|
is used in the include file \.{common.h}, and should therefore be defined
early; it refers to the |struct text| that will be declared later.

@d max_toks 150000L /* number of bytes in compressed \Cee~code */
@d max_texts 2500 /* number of replacement texts, must be less than 10240 */
@d max_files 50 /* number of auxiliary output files */
@d stack_size_max 50 /* module nesting level during output */
@d max_indent 1000 /* size of buffer for saving indentation characters */
@)
@d variant @;text

@ The program is built from two compilation units, one with source file
\.{common.w}, which contains a collection of routines and data shared
between \.{CTANGLE} and \.{CWEAVE}, and a second with source file
\.{ctangle.w} containing all code specific to \.{\me.}, and whose
typeset version you are now reading. All compilation units of the \.{CWEB}
system incorporate the file \.{common.inc} containing common declarations.

\.{\me.} has a fairly straightforward outline.  It operates in
two phases: first it reads the source file, saving the \Cee~code in
compressed form; then outputs the code, after shuffling it around.
It can optionally be compiled with the preprocessor symbol |STAT| defined,
in which case it will keep track of how much of \.{\me.}'s resources
were actually used. One command line argument is specifically designated to
modify the behaviour of \.{\me.}. Normally \.{\me.} will produce
machine-oriented output that supplies \&{\#line} directives (to allow
compilers and debuggers to locate the original source lines of statements),
while preserving line breaks but otherwise ignoring the lay-out and comments
in the source text. When the flag `\.{-l}' is supplied on the command line
however, the output cannot be traced back easily to the source lines, but
will be easier for humans to read: no \&{\#line} directives are produced,
but original lay-out and comments are preserved as well as possible. The
flag controlling this behaviour is called |line_output|; in the standard
mode of operation we have |line_output==true|.

@d line_output flags['l']
@c
@< Function prototypes used but not defined in the shared code @>@;
@< Typedef and enumeration declarations @>@;
@< Prototypes @>@;
@< Global variables @>@;

int main (int argc, char** argv)
{ program=ctangle;
  line_output=true;
  common_init(argc,argv,banner);
  @<Set initial values@>
  if (show_banner)
    print("%s, in %s mode.\n",banner,C_plus_plus ? "C++" : "C");
    /* print a ``banner line'' */
  phase_one(); /* read all the user's text and compress it into |tok_mem| */
  phase_two(); /* output the contents of the compressed tables */
  wrap_up(); /* and exit gracefully */
  return 0;
}

@i common.inc


@* Data structures exclusive to {\tt \me.}.
The basic mission of \.{\me.} is a relatively simple one, namely to
separate the \Cee~code from the commentary, to forget about the latter, and
to output the former almost verbatim, but plugging in the program text
corresponding to a module each time a reference is made. To this end the
bodies of all the sections are collected during Phase~I, and the proper
links are laid, after which the code is linearised for output in a
straightforward manner during Phase~II.

The \Cee~code itself will be stored as a sequence of bytes (more precisely,
of |eight_bits|) in a single large array |tok_mem|, much in the same way as
the characters for identifiers and module names are stored in |byte_mem|;
the encoding used will be described below. The individual sections are
represented by structures of type |text|, which contain pointers to the
|tok_mem| array as well as links by which they may be connected into lists;
the same structures are also used to represent the replacement texts of
macros. All |text| structures reside in an array |text_table|, which is
comparable to |id_table| and |mod_table|.

@<Typedef...@>=
typedef struct text
{ eight_bits* tok_start; /* pointer into |tok_mem| */
  struct text* text_link; /* relates replacement texts */
} text,* text_pointer;

@ Since all allocation is during Phase~I and all deallocation during
Phase~II, there is no reason for reclaiming memory when it is no longer
needed. Hence allocation is strictly linear, and a single pointer into
|tok_mem| for each |text| suffices, since a texts ends where the next one
starts. The first position of |tok_mem| that is unoccupied by a replacement
text is called |tok_ptr|, and the first unused location of |text_table| is
called |text_ptr|. Whenever we are not in the process of adding tokens to
|tok_mem|, we have the identity |tok_begin(text_ptr)==tok_ptr|.

@d tok_begin(p) (p)->tok_start
@d tok_end(p) ((p)+1)->tok_start
@d text_table_end (&text_table[max_texts])
@d tok_mem_end (&tok_mem[max_toks])

@<Glob...@>=
text text_table[max_texts];
text_pointer text_ptr=&text_table[0]; /* first unused position in |text_table| */
eight_bits tok_mem[max_toks];
eight_bits *tok_ptr=&tok_mem[0]; /* first unused position in |tok_mem| */

@ Invariants must be initialised.

@<Set init...@>=
tok_begin(text_ptr)=tok_ptr;

@ The following macro and function are used to enter one- and two-byte
tokens into |tok_mem| when a replacement text is being generated.

@d store_byte(c) @+
if (tok_ptr==tok_mem_end) overflow("token"); @+ else *tok_ptr++=c@;

@c
void store_two_bytes (sixteen_bits x)
{ if (tok_ptr+2>tok_mem_end) overflow("token");
  *tok_ptr++ = x >> 8; /* store high byte */
  *tok_ptr++ = x & 0xFF; /* store low byte */
}

@ When by successive calls of |store_byte| and |store_two_bytes| a complete
replacement text has been stored, the following code makes it into an
official text object by advancing |text_ptr| and storing the starting
location of the next text to be generated into the |text_table| array.

@< Wrap up the accumulated bytes into a completed text,
   pointed to by |cur_text| @>=
{ cur_text=text_ptr++; /* consolidate the replacement text */
  if (text_ptr>=text_table_end) overflow("text");
  tok_begin(text_ptr)=tok_ptr; /* mark end of replacement text */
}

@ The |text_link| field links a section to a possible continuation of it (the
next section with the same name, or for unnamed sections the next unnamed
section). The set of all sections linked together in such a way will be termed
a `module'. Apart from the value |NULL| used to indicate the end of such a
list, we need two other distinguished values: |macro_flag|, which indicates
that the text is a macro replacement text rather than a section body, and
|header_flag| similarly indicating that the text is the (quoted) file name
of a header file that is included. By a fortunate coincidence, there are two
unused pointer values near the end of |text_table|: since the macro |tok_end|
should work for all pointers to a replacement text, the last entry of
|text_table| is never used to store a replacement text; its address can be
used as one reserved pointer value, and the definition of~\Cee\ generously
offers us the address of the non-existent entry following it as another
valid pointer value. The pointer to the first section of the unnamed module
appears in |text_root|, and the address of the final link of this module is
recorded in |last_unnamed|. The named modules are accessed by the extra
pointer in the |mod_info| structure for a module name, which is a
|text_pointer| pointing to the body of the first section of that module;
in~\.{\me.} this field will be called |equiv|.

@d macro_flag (text_table_end-1) /* address of unused entry */
@d header_flag text_table_end /* address just beyond |text_table| */
@d next_sec(m) ((m)->text_link) /* next section of the same module, if any */
@d equiv equiv_or_xref /* info corresponding to names */

@<Glob...@>=
text_pointer text_root=NULL;
text_pointer* last_unnamed=&text_root; /* where to link new unnamed sections */

@ Here are the functions specific to \.{\me.}, which are used by the
common lookup routines. The function |names_match| decides whether a name of
length~|l| starting at position |q| equals the identifier represented
by~|x|. The parameter |dummy| is present only for compatibility with
\.{CWEAVE}, where it contains the |ilk| code.  The common lookup routines
refer to separate routines |init_module_name| and |init_id_name| when the
data structure grows. Actually |init_id_name| is called only when
|program==cweave|, but we need to declare a dummy version so that the linker
won't complain about its absence.

@c
boolean names_match (id_pointer x,char* q,int l,int dummy)
{@; char* p=name_begin(x); while (--l>=0) if (*p++!=*q++) return false;
  return *p=='\0';
}

void init_module_name(mod_pointer node)
@+{@; node->equiv=NULL; }

void init_id_name (id_pointer dummy,int ilk) @+ {}


@* Tokens.
Clearly, \.{\me.} must hold all the code for a complete program file at
the end of Phase~I (not counting the possible multiple use of some modules,
which will probably not make a great difference), and so it is worth while
to use a compact representation for the internal storage of the \Cee~code.
Using full-fledged data compression techniques would probably be overdoing
things a bit, but a few simple methods will make a great difference.  There
are likely to be many identical occurrences of identifiers and keywords of
the language, so it is efficient to collect them into a table and replace
them in the internal format by references to the table; module names have
to be entered in a table anyway for the purpose of associating occurrences
of the same name with each other. Furthermore white space and comments can
be removed upon storing the texts; any necessary white space can be
reinserted during the output process (but note that newlines must be
represented internally, in order to be able to match output lines with
input lines).

@ These compressed texts representing \Cee~code, while physically
represented as a sequence of eight-bit bytes in |tok_mem|, logically consist
of streams of `tokens', some of which occupy two or more consecutive byte
positions, while others take just one byte. If the first byte is |a| and,
in case |a>=0x80|, the next byte is~|b|, then the interpretation of the
token is as follows.

\Y
\item{$\bullet$} |0<=a<0x80|: the token represents the character~|a|;
\item{$\bullet$} |0x80<=a<0xA8|: the two-byte token represents the
  identifier with name |id_at((a-0x80)*@t$2^8$@>+b)|;
\item{$\bullet$} |0xA8<=a<0xD0|: the two-byte token represents the module
  with name |mod_at((a-0xA8)*@t$2^8$@>+b)|;
\item{$\bullet$} |0xD0<=a<0xF8|: this two-byte token marks the beginning of
  (a part of) the repacement text for the current module, defined in section
  number |(a-0xD0)*@t$2^8$@>+b|;
\item{$\bullet$} |a==0xF8|: the two-byte token represents the (8-bit)
  character~|b|.
\item{$\bullet$} |a==0xF9|: the token consists of 5 bytes, say |a|,~|b|,
  |c|, |d| and~|e|, and represents a \&{\#line} directive for line number
  |b*@t$2^8$@>+c| in the file |id_at(d*@t$2^8$@>+e)|.
\item{$\bullet$} |a==0xFA|: the token consists of a single byte, and
  represents the location, indicated by the \:p control codes in the source
  file, where the preprocessor directives produced by \:d and \:h will be
  included in the output; this location will be marked by the presence of
  the virtual module called `|@p|' in the typeset output.
\Y

Any \Cee~token that consists of a single 7-bit character is represented by
that character itself; in particular, a single-character identifier like
`|x|' will be a one-byte token, while all longer identifiers will occupy two
bytes. Some of the 7-bit codes will not be present however, since they are
control codes and they do not represent a symbol useful for \Cee\ in the
\caps{MIT} extension of \ASCII., so we can use them for special purposes.
The following symbolic names are used:

\yskip
\hang |join| denotes the concatenation of adjacent items with no space or
  line breaks allowed between them (the \:\& operation of \.{CWEB}).

\hang |verb_quote| denotes the beginning or end of a sequence of bytes that
  is to be copied verbatim to the output, i.e., no spaces should be inserted.
  Such sequences are used to transmit numerical constants and strings (in
  which case the ordinary string quotes appear within the |verb_quote|
  bytes), and to implement the \:= operation.

@^ASCII code dependencies@>

@d verb_quote 0x2 /* takes the place of extended \ASCII. \.{\char2} */
@d join 0x3 /* takes the place of extended \ASCII. \.{\char3} */

@< Global var... @>=
boolean atp_seen=false; /* whether \:p occurs anywhere */


@* Phase II processing.
We will start to explain the structure of the tokens lists built up in
memory by considering how they are output during Phase~II; this is the
simpler part of the program. We will proceed in top-down fashion for this
part.

@<Prototypes@>=
void phase_two (void); /* output the contents of the compressed tables */
void output (text_pointer); /* recursively write out modules */
void output_preproc_directives(void); /* write out all \:d and \:h stuff */
void C_newline (void); /* send a newline to the \Cee~file */
void out_char (eight_bits); /* write a single token */

@ Here is the general control routine.	After Phase~I, the quantity
|output_file_count| will be equal to the number of extra output files to be
written. There should at least be either an unnamed section or one that
produces an additional output file, for otherwise there will be no output at
all; if the is no unnamed section, \:d~and~\:h commands will only have
effect if \:p occurs somewhere.

@c void phase_two (void)
{ phase=2;
  if (text_root==NULL && output_file_count==0)
  @/{@; print("\n! No program text was specified."); mark_harmless(); }
		 @.No program text...@>
  else
  { if (show_progress)
    { print("\nWriting the output file%s"   @.Writing the output...@>
	   ,(text_root!=NULL)+output_file_count>1 ? "s" : "");
      if (text_root!=NULL) printf(" (%s):",C_file_name);
      update_terminal();
    }
    if (text_root==NULL) C_file=NULL;
    else
    { open_output_file(); cur_line=1;
      if (!atp_seen) output_preproc_directives();
      output(text_root);
    }
    @<Write all the named output files@>
    print_progress("\nDone.\n");
  }
}

@ \.{\me.} can write output on multiple files. If a module name is
introduced in at least one place by \:( instead of \:<, it is treated as the
name of a file. All these special module names are saved in the array
|output_files|.

@<Glob...@>=
mod_pointer output_file[max_files];
int output_file_count=0;

@ To write the named output files, we proceed as for the unnamed module, the
only subtlety is that we have to open each one after closing the previous
one. However, if there is no main file, there is no file to close the first
time, and this situation can be recognised because we have set |C_file=NULL|
in that case. The length of a module name can be at most |longest_name|
(which is rather more generous than allowed elsewhere for file names), so we
need no test for overflowing |output_file_name|.

@<Write all the named output files@>=
{ int i;
  char output_file_name[longest_name+1]; /* name of the file */
  for (i=0; i<output_file_count; i++)
  { mod_pointer output_module=output_file[i];
    @< Copy name of |output_module| into |output_file_name| @>
    if (C_file!=NULL) fclose(C_file);
    if (output_module->equiv==NULL)
    @/{@; print("\n! Module not present");
		   @.Module not present@>
      print_mod(output_module); err_print("");
    }
    else if ((C_file=fopen(output_file_name,"w"))==NULL)
    @/{@; print("\n! Cannot open \"%s\" as output file",output_file_name);
		   @.Cannot open output file@>
      err_print("");
    }
    else
    { if (show_progress) print("\n(%s):",output_file_name);
      cur_line=1;
      output(output_module->equiv);
    }
  }
}

@ If the module name defining the output file contains any occurrences of
`\.{@@@@}', they are undoubled in forming the actual file name; apart from
this should be no control codes in the module name.

@< Copy name of |output_module| into |output_file_name| @>=
{ char* p=output_file_name, *q=name_begin(output_module);
  while (*q!='\0')
    if ((*p++=*q++)=='@@')
      if (*q++!='@@')
      { print("\n! Illegal control code in file name");
		 @.Illegal control code in file name@>
	print_mod(output_module); err_print("");
      }
  *p='\0';
}

@ Here is how a file name gets added to the list of output files.

@<Adjoin |cur_mod| to the set of output files@>=
{ int i=0;
  while(i<output_file_count)
    if (output_file[i]==cur_mod) break; @+ else ++i;
  if (i==output_file_count) /* if not present, add |cur_mod| */
    if (output_file_count==max_files) overflow ("output files");
    else output_file[output_file_count++]=cur_mod;
}


@* Writing replacement texts.
We now come to the functions |output| and |output_preproc_directives|; the
task of the former is a recursive one, and it may also call the latter if
\:p occurs. The function |output| is not recursive, however, but it uses a
stack to keep track of what is going on at different ``levels'' as the
sections are being written out. Entries on this stack have five parts:

\yskip

\hang |repl_field| points to replacement text of the active section;

\hang |byte_field| is the |tok_mem| location from which the next token
      on a particular level will be read;

\hang |end_field| is the |tok_mem| location where the replacement
      text of a particular level will end;

\hang |sec_nr_field| is the section number.

\hang |indent_field| determines the amount of indentation, for in case
      |line_output==false|

\yskip\noindent
The current values of these quantities are referred to quite frequently, so
they are stored in a separate place instead of in the |stack| array. We call
the current values |cur_byte|, |cur_end|, |cur_repl|, |cur_sec| and
|cur_ind|.

@d cur_repl cur_state.repl_field
@d cur_byte cur_state.byte_field
@d cur_end cur_state.end_field
@d cur_sec cur_state.sec_nr_field
@d cur_ind cur_state.indent_field

@<Typedef...@>=
typedef struct
{ text_pointer repl_field; /* replacement text of active section */
  eight_bits *byte_field; /* present location within replacement text */
  eight_bits *end_field; /* ending location of replacement text */
  sixteen_bits sec_nr_field; /* section number */
  sixteen_bits indent_field; /* amount of indentation */
} output_state, * stack_pointer;

@ The global variable |stack_ptr| tells how many levels of output are
currently in progress; the stack grows upwards with |stack_ptr| pointing to
the first vacant location above the top of the stack. The entry |stack[0]|
at the bottom of the stack is not used; when |stack_ptr| points to it then
|cur_state| itself has become invalid and the output process is completed;
this condition can be tested as |stack_empty()|. This somewhat strange use of
the stack is forced upon us because the \Cee\ language does not guarantee
that a pointer value at offset~$-1$ of an allocated array exists. Since the
location |stack[0]| is vacant we might as well use it to store |cur_state|
in.

@d cur_state stack[0]
@d stack_end (&stack[stack_size_max])
@d stack_empty() (stack_ptr==&stack[0])

@<Global...@>=
output_state stack[stack_size_max]; /* info for non-current levels */
stack_pointer stack_ptr; /* points above top of output state stack */

@ When |line_output==false| we keep track of the current indentation level
in an array of white space characters |indent_buffer|. Every time a newline
is output, a sequence of |cur_ind| characters from |indent_buffer| are also
output by invoking |put_indent|. When relevant non-newline characters are
output, a white space counterpart of the character is appended to
|indent_buffer| by invoking the macro |append_white|, whose argument should
not involve any side effects. The current position of writing the characters
is maintained in |ind_i|, and reset to |cur_ind| on output of a newline.

@d C_printf(format,x) fprintf(C_file,format,x)
@d C_putc(c) putc(c,C_file)
@d put_indent()
   (indent_buffer[ind_i=cur_ind]='\0',C_printf("%s",indent_buffer))
@d append_white(c)
  if (ind_i>=max_indent) overflow("indent buffer");
  else indent_buffer[ind_i++]= isspace(c) ? c : ' ' @;
@< Global var... @>=
char indent_buffer[max_indent];
   /* white space characters to produce indentation */
sixteen_bits ind_i;
   /* number of those characters needed to reach current position */

@ When the replacement text for name |p| is to be inserted into the output,
the following function is called to save the old level of output and get
the new one going. The value of |cur_sec| is not set since its value is
recorded in the first token to be read after this function is called. It is
by calling |push_level| that the current value of |ind_i| becomes
significant, by being copied to |cur_ind|.

@c
void push_level (mod_pointer p) /* suspends the current level */
{ if (stack_ptr==stack_end) overflow("output stack");
  *stack_ptr++=cur_state;
@/cur_repl=p->equiv;
  cur_byte=tok_begin(cur_repl); cur_end=tok_end(cur_repl);
@/cur_ind=ind_i;
}

@ When we come to the end of a replacement text, the function
|continue_or_pop_level| does what its name suggests: it either moves to the
continuation of this replacement text (i.e., to the next section of the same
module) or returns the state to the most recently stacked level when a
module is completed.  In the latter case all the fields of |cur_state| are
restored; in the former case, like for |push_level|, the value of |cur_sec|
is not yet set, but will be adjusted immediately afterwards.

@c
void continue_or_pop_level (void)
{ if (cur_repl->text_link != NULL)
    /* then link to a continuation, staying on the same level */
  {@; cur_repl=next_sec(cur_repl);
    cur_byte=tok_begin(cur_repl); cur_end=tok_end(cur_repl);
  }
  else if (--stack_ptr > &stack[0]) cur_state=*stack_ptr;
     /* otherwise |stack_empty()| holds */
  ind_i=cur_ind; /* reset indentation to base of current level */
}

@ The function |output| handles output of tokens by sending them to a lower
level function |out_char|, until the condition |stack_empty()| holds. The
task of |output| is to isolate and decode tokens, to handle stacking and
unstacking as necessary, and to issue special calls to |out_char| to mark
the beginning and end of replacements texts and to produce appropriate
\&{\#line} directives. As its name suggests, |out_char| usually handles the
output of single characters, but some special codes will make it perform a
few other operations. Some of these, like |verb_quote| and compressed
operator codes, occur as bytes with values below |0x80| inside the
replacements texts, and these get no special treatment by |output|. A few
other ones are generated explicitly by |output|, and these have values from
|0x80| upwards, defined in the enumeration below. Information about the
precise value of two-byte tokens is passed to |out_char| via the global
variable |cur_val|.

@<Global...@>=
int cur_val; /* additional information corresponding to output token */

@~The code |identifier| is used for identifiers of length two or more, in
which case |cur_val| indicates the identifier name. The codes
|section_start| or |section_end| are sent at the beginning and end of parts
of the replacement text defined in a one section, in which case |cur_val| is
the section number; this number will be recorded in a comment to aid human
readers of the \Cee~file produced. When a \&{\#line} directive is to be
output, |line_mark| is sent to |out_char|, which will fetch the following
4~bytes (specifying the line number and file name) from the replacement text
itself; this is the only occasion where |output| does not read a complete
token before calling |out_char|. Any 8-bit characters that were escaped by a
byte `|0xF8|' are sent directly to |out_char|; they do not interfere with
the codes below since the input routine guarantees that they only occur
between |verb_quote| tokens, at which times |out_char| does not honour any
special codes.

@< Typedef and enumeration... @>=
enum {identifier=0x80, section_start, section_end, line_mark };

@ Due to the output stack operations, the function |output| can be written
iteratively rather than recursively. A user who is not afraid of recursion
may however easily rewrite this code into a recursive form, and dispose of
the output stack and its routines altogether; the recursion depth will be
equal to the level of nesting of modules, which is not likely to cause any
problems in practice. We have retained the iterative version out of
reverence for the author of the original |WEB| system; it has the
additional advantage that, should there be a cycle in the directed graph of
module references, a reasonable overflow message will be produced by
|push_level|, whereas otherwise system stack overflow would occur, and it
would depend on the \Cee~runtime support whether it is properly detected.

@c
void output (text_pointer repl) /* sends tokens to |out_char| */
{ stack_ptr=&stack[1];
@/cur_repl=repl; cur_byte=tok_begin(cur_repl); cur_end=tok_end(cur_repl);
  cur_ind=ind_i=0;
  do
    if (cur_byte==cur_end)
    { cur_val=cur_sec;
      continue_or_pop_level();
      out_char(section_end); /* output ending section number comment */
    }
    else
    @< Output the token starting at |*cur_byte|, advancing |cur_byte|
       correspondingly @>
  while(!stack_empty());
  C_newline();
}

@ Most tokens lead to calling |out_char| with an appropriate value, after
possibly setting |cur_val| and |cur_sec|. Exceptions are module names, which
will push their replacement text on the stack but not produce any direct
output, and tokens |0xFA| which will produce the collected preprocessor
directives by calling |output_preproc_directives|, after ensuring that this
will start on a fresh line of output.

@< Output the token starting at |*cur_byte|... @>=
{ int a = *cur_byte++;
  if (a<0x80) out_char(a); /* single byte token */
  else if (a>=0xF8)
    if (a<0xFA) out_char(a==0xF8 ? *cur_byte++ : line_mark);
    else {@; C_newline(); output_preproc_directives(); } /* |a==0xFA| */
  else
  { cur_val=(((a-=0x80)%0x28)<<8)+*cur_byte++;
    switch (a/0x28)
    { case 0: out_char(identifier); break;
      case 1: @<Expand module name |cur_val| @> @+ break;
      case 2: cur_sec=cur_val; out_char(section_start);
	  /* set the correct section number and output comment */
    }
  }
}

@ If any defining occurrence of a module name was encountered during the
Phase~I, its replacement text will have been be linked into to |equiv| field
of that name; otherwise that field will contain a null pointer and we must
report an error.

@<Expand module name...@>=
{ mod_pointer mod_name=mod_at(cur_val);
  if (mod_name->equiv!=NULL) push_level(mod_name);
  else
  {@; print("\n! Module not present"); print_mod(mod_name); err_print("");
	       @.Module not present@>
  }
}

@ Output of macro definitions and inclusions of header files is sufficiently
different from ordinary output that we use a special function
|output_preproc_directives| for it. We bypass |output| and call |out_char|
directly, since in these cases there can be no nesting of replacement texts.
We go through the list of all replacement texts and copy the ones that refer
to macros or header files, preceded respectively by \&{\#define} and
\&{\#include}. During the output of a single directive any line breaks
present in the source file will be protected by backslashes. We use local
variables |p| and |end| instead of |cur_byte| and |cur_end|, since the
current level of the stack is already in use if this function is called from
|output|.

@c
void output_preproc_directives (void)
{ text_pointer repl,l;
  eight_bits* p, *end;
  protect=true; /* newlines should be preceded by |'\\'| */
  for (repl=&text_table[0]; repl<text_ptr; repl++)
    if ((l=repl->text_link)==macro_flag || l==header_flag)
    { p=tok_begin(repl); end=tok_end(repl);
      C_printf ("#%se ",l==macro_flag ? "defin" : "includ");
      out_state=no_space;
      while (p<end)
      @< Output preprocessor token starting at |p|, advancing |p|
         correspondingly @>
      C_newline(); /* this newline is not escaped */
    }
  protect=false;
}

@ The situation here is simpler than in |output|, since module names and \:p
control codes are forbidden in the replacement texts handled here, nor will
codes for \.{\#line} directives or section number indications have been
included in them during input. Should any such token nevertheless turn up,
then one of the two calls of |confusion| below is invoked; the diagnostic is
not quite exhaustive as to the possible cause in either case, but we do not
wish to elaborate on code that will never be executed.

@< Output preprocessor token... @>=
{ int a=*p++;
  if (a<0x80) out_char(a); /* single byte token */
  else if (a>=0xF8)
    if (a==0xF8) out_char(*p++); @+
    else confusion("`@@p' within macro"); @.`@@p' within macro@>
  else
  { cur_val=(((a-=0x80)%0x28)<<8)+*p++;
    if (a<0x28) out_char(identifier); @+
    else confusion("module within macro"); @.module within macro@>
  }
}


@* Writing characters. The |output| routine above handles the global
structure of output generation; we now present the routines that transform
the lexical items produced by |output| into characters.
First we give a function that is called whenever we want to finish off a
line sent to the \Cee~file.
It keeps |cur_line| equal to the number of the next line to be output,
and displays a progress report every 100 lines.

@c
void C_newline (void) /* writes one line to output file */
{ C_putc('\n');
  if (!line_output) put_indent();
  if (cur_line%100==0 && show_progress)
  { if (cur_line%500!=0) print("."); /* progress report */
    else print(cur_line%2500==0 ? "%u\n" : "%u",cur_line);
    update_terminal();
  }
  ++cur_line;
}

@ The function |out_char| must make sure that the output has the proper
``surface structure''. When |line_output==false| everything should
essentially be written as it is stored, but when |line_output==true| only
the essential tokens have been stored, and white space has to be inserted
where appropriate. Spaces should not occur at certain places (e.g., not in
the middle of a string or a constant or an identifier, not at a \:\&
position where quantities are being joined together), while in other places
they are required (e.g., between identifiers, which for \.{\me.} includes
reserved words, and between certain operators that might otherwise be
considered as a single token). Such surface structure can very nicely be
obtained by attaching a small finite state machine to the output generator,
recording information about the most recent token (as a matter of fact, even
\TeX's fabulous formatting of math formulae is largely based on such a simple
device). In our present case, the state of the output process is recorded in
the global variable |out_state|. Furthermore there is a variable |protect|
that is explicitly managed by the function calling |out_char| (notably by
|output_preproc_directives|), being set to |true| whenever newlines are to
be escaped by a backslash.

@<Global...@>=
eight_bits out_state; /* current status of partial output */
boolean protect; /* should newline characters be quoted? */

@ The output process can be in one of following states:

\yskip\hang |no_space| means that no space will precede the following item.
This state is set by \:\&, and after the output of anything that could not
possibly combine into a larger token, such as punctuation.

\yskip\hang |num_or_id| means that the last item in the buffer is a number or
identifier, hence a blank space or line break must be inserted if the next
item is also a number or identifier.

\yskip\hang |operator| means that the last item in the buffer is an operator,
hence a blank space or line break must be inserted if the next item is also
an operator.

\yskip\hang |literal| means we're copying only character tokens, and
that they are to be output exactly as stored.  This is the case during
strings, verbatim constructions and numerical constants; the |verb_quote|
character serves as a quoting mechanism which switches this state
on~and~off.

@< Typedef and enum... @>=
enum { no_space, num_or_id, operator, literal};

@ The function |out_char| is a many-way switch on the kind of character or
token transmitted to it. The |verb_quote| mechanism is used amongst others
to encapsulate numeric constants, so we make the whole construction behave
as an identifier with respect to spacing (the user can still defeat any
spaces around it by means of \:\&). The code |line_mark| will only be stored
(and output) if |line_output==true|. We do not anticipate the need of
indentation across literals, operators and such, and assume that module
names are preceded on the same line only by white space, identifiers
(keywords) and things like braces; apart from identifiers all of these come
through the |default| case below. Therefore we only include an invocation of
|append_white| in that case and for identifiers; to handle all possible
cases the assiduous user could add such invocations to all other calls of
|C_putc| and |C_printf|.

@c
void out_char (eight_bits c)
{ if (out_state==literal)
    if (c==verb_quote) out_state=num_or_id; /* close literal mode */
    else C_putc(c); /* write literal character, possibly 8-bit */
  else if (isalnum(c)&&c<0x80 || c=='_')
    /* single character identifier or number */
  { if (out_state==num_or_id && line_output) C_putc(' ');
    C_putc(c); out_state=num_or_id;
    if (!line_output) append_white(c);
  }
  else switch (c)
  { case verb_quote:
      if (out_state==num_or_id && line_output) C_putc(' ');
      out_state=literal; break;
    case join: out_state=no_space; break;
    case '\n': @< Output a newline @> @+ break;
    case identifier: @< Output an identifier @> @+ break;
    case section_start:
      if (line_output) C_printf("/*%d:*/",cur_val); @+ else C_newline();
      out_state=no_space; break;
    case section_end:
      if (line_output) C_printf("/*:%d*/",cur_val); @+ else C_newline();
      out_state=no_space; break;
    case line_mark: @< Output a \.{\#line} directive @> @+ break;
    @\@<Cases of operator symbols@>
    @\@<Cases of compressed operators@>
    default: C_putc(c); out_state=no_space;
      if (!line_output) append_white(c);
  }
}

@ Newlines are escaped whenever |protect| holds; a space is prepended so
that the newline is effectively replaced by a space. This also means that
the user who cannot break the habit of escaping newlines in multi-line
preprocessor directives, even though this is unnecessary after \:d and \:h,
is not punished for this by the creation of the sequence `\.{\\\\}' at the
end of the lines, but rather gets `\.{\\\ \\}', which is harmless.

@< Output a newline @>=
{ if (protect) {@; C_putc(' '); C_putc('\\'); }
  C_newline();
  if (out_state!=literal) out_state=no_space;
}

@ Apart from alphanumeric characters and underscores, we allow identifiers
to contain characters in the range |0x80<=c<=UCHAR_MAX|, if a translation
into a string of legal characters is known to \.{\me.}; such translations
can be specified by \:l directives in limbo. These translations are stored
in an array |c_trans|.

@d trans_limit 9 /* maximal length of a translation string */
@d trans_of(c) c_trans[(unsigned char)(c)-0x80]
@d translation_exists(c) (trans_of(c)[0]!='\0') /* non-empty translation */

@< Global var... @>=
char c_trans[UCHAR_MAX+1-0x80][trans_limit+1];

@~In compatibility mode, default translations exist of the form `\.X\\{NN}',
where \\{NN} is the (two digit) upper case hexadecimal representation of the
character in question. Otherwise the static initialisation will have made all
translation strings empty.

@< Set init... @>=
if (compatibility_mode)
@/{@; unsigned char c=UCHAR_MAX;
  do sprintf(trans_of(c),"X%X",c); while (--c>=0x80);
}

@ In case of identifiers the name indexed by |cur_val| is written out. The
translation of characters in the range from |0x80| upwards is performed
here, rather than on input of such characters, since this is slightly more
convenient. The check that |translation_exists| was done on input, however.

@< Output an identifier @>=
{ char* p=name_begin(id_at(cur_val)); int l=0;
  if (out_state==num_or_id && line_output) C_putc(' ');
  do
    if ((unsigned char)(*p)<0x80) {@; C_putc(*p);++l; }
    else
    {@; char* q=trans_of(*p); do {@; C_putc(*q); ++l; } while (*++q!='\0'); }
  while (*++p!='\0');
  out_state=num_or_id;
  if (!line_output) do append_white(' '); while (--l>0);
}

@ As indicated above, a call |out_char(line_mark)| causes four further bytes
to be fetched directly from the token memory, increasing |cur_byte| so that
|output| will resume at the correct byte.

@< Output a \.{\#line} directive @>=
{ sixteen_bits a;
  a=(*cur_byte++)<<8; a+=*cur_byte++; /* get the line number */
  C_newline(); C_printf("#line %u \"",a); @/
  a=(*cur_byte++)<<8; a+=*cur_byte++; /* get the file name index */
  C_printf("%s\"",name_begin(id_at(a)));
  C_newline(); out_state=no_space;
}

@ In some cases when two operator symbols are adjacent, a space is required
in between to avoid interpretation as a compound symbol. The language
definition is not explicit about when such a space is mandatory, but
examples are the somewhat unusual constructions | x /@, *p |, | a++ + b != a
+ ++b| and |sum = a - -b|; certain backward \Cee~compilers even think |x=-x|
is ambiguous without a space before the minus sign. Rather than trying to
detect precisely the problematic cases we take the conservative approach of
always putting in a space, unless the second operator is~`\.='; this
exception is made because we store operators like `\.{+=}' as two separate
tokens, but they form a single token according to the \caps{ANSI/ISO}
standard. Note that other multi-character operators like `\.{\&\&}' stay
intact because they will be compressed on input, and therefore do not involve
the present section.

@<Cases of operator symbols@>=
case '+': case '-': case '*': case '/': case '%': case'?':
case '<': case '>': case '&': case '|':
  if (out_state==operator && line_output) C_putc(' '); /* fall through */
case '=': C_putc(c); out_state=operator; break;

@ Compilers don't mind doing repetitive work, so we use a macro to produce
several almost identical cases. Had the cases been much longer or more
numerous however, then combining the cases and using an array of strings
would have been preferable.

@d comp_op(op)
  (C_printf(out_state==operator && line_output ? " %s" : "%s",op)
  ,out_state=operator)

@<Cases of compressed operators@>=
case plus_plus: comp_op("++"); break;
case minus_minus: comp_op("--"); break;
case minus_gt: comp_op("->"); break;
case gt_gt: comp_op(">>"); break;
case eq_eq: comp_op("=="); break;
case lt_lt: comp_op("<<"); break;
case gt_eq: comp_op(">="); break;
case lt_eq: comp_op("<="); break;
case not_eq: comp_op("!="); break;
case and_and: comp_op("&&"); break;
case or_or: comp_op("||"); break;


@* Introduction to the phase I.  We have now seen that \.{\me.} will be
able to output the full \Cee\ program, if we can only get that program into
the byte memory in the proper format. The input process is something like
the output process in reverse, since we compress the text as we read it in
and we expand it as we write it out (but of course, as any programmer knows,
input requires a bit more work than the corresponding output). To preserve
the symmetry we shall proceed in a bottom-up fashion in presenting this part
of the program.

The basic character input routines are defined in the code shared with
\.{CWEAVE}. At the next higher level there are three main input routines.
The most interesting is the one that gets the next token of a \Cee\ text;
the other two are used to scan rapidly past \TeX\ text in the \.{CWEB}
source code. One of the latter routines will jump to the next token that
starts with `\.{@@}', and the other skips to the end of a \Cee~comment.

@ Control codes in \.{CWEB} begin with `\.{@@}', and the next character
identifies the code. Some of these are of interest only to \.{CWEAVE}, so
\.{\me.} ignores them; the others are converted by \.{\me.} into
internal code numbers. The code numbers have been chosen such that they are
distinct from ordinary characters (so that a function can return either a
character or a control code), and their ordering is such as to simplify the
program logic; larger numbers are given to the control codes that denote
more significant milestones. This order is only relevant starting from
|format| however. The code \:> is treated as ignored because it should not
occur in places where we are not searching for it. The following enumeration
lists the non-character values that may be returned from~|get_next|;
therefore the two values |id_code| and |constant| are included although they
do not belong to any control code.

@<Typedef and enum...@>=
enum @/
{ ignore=0x80, /* control code of no interest to \.{\me.}, or \:> */
  id_code, constant, /* token codes but not control codes */
  verbatim, /* control code for \:= */
  at_sign_image, /* control code for \:@@ */
  join_code, /* control code for \:\& */
  ord, /* control code for \:' */
  control_text, /* control code for \:t, \:\^, etc. */
  include_preproc, /* control code for \:p */
  char_trans, /* control code for \:l */
  format, /* control code for \:f */
  definition, /* control code for \:d */
  header, /* control code for \:h */
  begin_C, /* control code for \:c */
  module_name, /* control code for \:< or \:( */
  new_section /* control code for \:\ , \:\~, or \:* */
};

@ The conversion from the character following `\.@@' to the corresponding
numeric code is performed by means of the table~|ccode|. Since we don't know
whether characters are signed or not, we always access |ccode| via the macro
|code_of|.

@d code_of(c) ccode[(unsigned char)(c)]

@<Global...@>=
eight_bits ccode[UCHAR_MAX+1]; /* meaning of a char following `\.@@' */

@~Here we initialise |ccode|. The code code for~\:v is the only case where a
value stored in the |ccode| table is an ordinary character; this is done so
that it may be used in the \Cee~part of a section to stand for~`|@v|', even
though there is no need for using it there.

@<Set ini...@>=
{ unsigned char c=0;
  do ccode[c] = isspace (c) ? new_section : ignore; while(c++!=UCHAR_MAX);
@/ccode['v']=ccode['V']='|';
@/ccode['=']=verbatim;
@/ccode['@@']=at_sign_image;
@/ccode['&']=join_code;
@/ccode['\'']=ord;
@/ccode['^']=ccode['?']=ccode['.']=ccode[':']=ccode['#']=
    ccode['t']=ccode['T']=ccode['q']=ccode['Q']=control_text;
@/ccode['p']=ccode['P']=include_preproc;
@/ccode['l']=ccode['L']=char_trans;
@/ccode['f']=ccode['F']=ccode['s']=ccode['S']=format;
@/ccode['d']=ccode['D']=definition;
@/ccode['h']=ccode['H']=header;
@/ccode['c']=ccode['C']=begin_C;
@/ccode['<']=ccode['(']=module_name;
@/ccode['~']=ccode['*']=new_section;
if (compatibility_mode)
  @< Reset some control codes to match Levy/Knuth \.{CWEB} @>
}

@ In \.{CWEBx} there are a few control codes that also exist in \LKC. but
have a different meaning. In compatibility mode we reassign the meaning of
these codes to that of \LKC., making their usual function inaccessible,
since it is not intended that hybrid programs should be written using the
codes of \LKC. together with features particular to \.{CWEBx}. For
\.{\me.} the change of meaning of \:: (which is used instead of \:? in
\LKC.) has no effect since it introduces a control text anyway.

@< Reset some control codes... @>=
{ ccode['h']=ccode['H']=include_preproc; /* \:h means \:p */
  ccode['p']=ccode['P']=begin_C; /* \:p means \:c */
  ccode['#']=ignore; /* \:\# means \:) */
}

@ The function |skip_ahead| reads through the input at fairly high speed
until finding the next non-ignorable control code, which it returns. It is
used to skip over \TeX~text, and text in limbo. It uses the fact that
|get_line| places a~|' '| at~|*limit|, so that placing the sentinel cannot
inadvertently create a token \:@@.

@c
eight_bits skip_ahead (void) /* skip to next control code */
{ eight_bits c; /* control code found */
  while (find_char())
  { limit[1]='@@'; /* place a sentinel */
    while (*loc++!='@@') {}
    if (loc<=limit && (c=code_of(*loc++))!=ignore) return c;
  }
  return new_section;
}

@ The function |skip_comment| reads through the input until finding the
end-comment token `\.{*/}' or a newline, or in case of a \Cpp\ one-line
comment (starting with `$/\!/$') just until finding a newline. The |one_liner|
parameter tells whether the latter applies, and a boolean result is returned
telling whether a newline in the middle of a comment was scanned (which can
only happen if |one_liner==false|). Returning such a newline is necessary so
that the each newline in the \Cee\ part of a section may be copied to the
output; otherwise the \&{\#line} commands inserted into the \Cee\ file by
the output routines become useless. When |line_output==false|, all characters
encountered are also stored, so that they can be output along with the
\Cee~code.

The function gives an error message if it encounters a sequence `\.{/*}',
since this almost certainly means the user has forgotten to close the
comment. In this way we avoid the most pernicious consequence of \Cee's rule
that nested comments are forbidden, namely that such an omission often leads
to a syntactically correct, but unintended program. This means that one
cannot write either `\.{/*}' or `\.{*/}' in a comment; this difficulty can
usually be circumvented by inserting something like `\.{\{\}}' between the
two characters.

As a safety measure, |skip_comment| comes to an end if it runs into the next
section, and prints an error message; it is slightly non-trivial that this
test does not interfere with searching for the end of the comment. If a
module name occurs in a comment (presumably in `\pb'), it is scanned by the
normal routines that will be given later: this could be the only place where
the name is spelled out in full, and/or using \:(.

@c
boolean skip_comment (boolean one_liner) /* skips over comments */
{ char c; /* current character */
  do
  { if (loc>=limit)
      if (one_liner) return false;
      else if (get_line ()) return true;
      else
      {@; err_print("! Input ended in mid-comment"); return false; }
		     @.Input ended in mid-comment@>
    if ((c=*loc++)=='/' && *loc=='*')
      err_print("! `/*' inside comment, did you forget `*/' before? ");
		 @.`/*' inside comment...@>
    if (c=='@@')
    { eight_bits cc=code_of(*loc++);
      if (cc==new_section) /* watch out for \:\ , \:\~, and \:* */
      @/{@; err_print("! Section ended in mid-comment"); loc-=2; return false; }
		       @.Section ended in mid-comment@>
      if (cc==module_name) @/{@; @< Scan the module name... @> continue; }
    }
    if (!line_output)
    {@; store_byte(c);
      if (c=='@@' && (c=loc[-1])!='@@') store_byte(c);
    }
  } while (c!='*' || *loc!='/' || one_liner);
  ++loc; @+ if (!line_output) store_byte('/');
  return false;
}


@* Getting the next token.
We now come to the function |get_next|, which does more extensive scanning.
It will move the input pointer over a single token of input, and if that
token consists of a single character that is not an identifier or a number,
it will return that character; certain two-character operators are also
returned as a single character representing the operator. In the case of
longer tokens, it will return one of the codes of the enumeration given
earlier, after possibly setting some global variables in such a way that the
complete token can be retrieved. For numeric constants |id_first| and
|id_loc| will be made to point to the beginning and end of the token
within~|buffer|. For string constants, which need not reside in a single
line of input, the characters are copied into the separate buffer |mod_text|
and this is where |id_first| and |id_loc| will then point. Identifiers and
module names are even looked up, and |cur_id| or |cur_mod| is made to point
to the name found. In the case of single-character identifiers or bad module
names, which are not entered into the name tables, |cur_id| respectively
|cur_mod| is set to~|NULL|; in the former case |id_first| is made to point
to the character that forms the identifier.

@<Global...@>=
id_pointer cur_id; /* identifier just scanned */
mod_pointer cur_mod; /* module name just scanned */

@ As one might expect, |get_next| consists mostly of a big switch that
branches to the various special cases that can arise. The static flag
|preprocessing| is raised while we are scanning an explicit preprocessor
directive (i.e., not one introduced by \:d or \:h). The static flag
|comment_continues| is raised when |get_next| must return a newline while
scanning a comment, so that it will remember what it was doing the next time
it is called. Although we allow 8-bit character input, any characters
|c>=0x80| should be contained in a comment or a character or string constant,
or in case |translation_exists(c)| an identifier, since they cannot
otherwise be part of a valid \Cee~symbol.

@c
eight_bits get_next (void) /* produces the next input token */
{ static boolean preprocessing=false; /* did this line start with `\.\#'? */
  static boolean comment_continues=false; /* were we scanning a comment? */
  eight_bits c; /* the current character */
restart:
  @< Handle cases involving newlines and comments;
     if a token is found |return| it, or |goto restart| if it is ignored;
     otherwise read the next character into |c| @>
  if (c=='L' && (*loc=='\'' || *loc=='\"'))
  {@; get_string(); return constant; }
  if (c<0x80 ? isalpha(c) || c=='_' : translation_exists(c))
  @/{@; @< Get an identifier @> return id_code; }
  if (c>=0x80) {@; err_print("! Illegal 8-bit character"); goto restart; }
			      @.Illegal 8-bit character@>
  if (isdigit(c) || c=='.' && isdigit((eight_bits)*loc))
  @/{@; @< Get a numeric constant @> return constant; }
  switch(c)
  { case '\'': case '"': get_string(); return constant;
    case '@@': @< Get a control code and possibly a complete module name,
		  and either |return| it,
		  or |goto restart| if control code is |ignore| @>
  @\@< In applicable cases |return| compressed two-symbol operator @>
  }
  return c;
}

@ The flag |preprocessing| raised when the first character of a line
is~`\.\#', and is lowered at a non-escaped newline; while it is raised spaces
are not ignored. The test whether a newline is escaped is not entirely
sound, since a final `\.{\\}' might be part of a \:\\ control code or of a
\Cpp\ one-line comment (after `$/\!/$'); a completely correct check would be
more difficult, and leaving |preprocessing| raised a bit too long does no
real harm. Contrary to |preprocessing|, the flag |comment_continues| will
always be given a new value on the first call to |get_next| after it is set,
and therefore can remain set only if |skip_comment| finds another newline.

@< Handle cases involving newlines... @>=
{ if (loc>=limit)
  { if (preprocessing && limit>buffer && limit[-1]!='\\')
      preprocessing=false;
    return get_line() ? '\n' : new_section;
  }
  if (comment_continues
   ||(c=*loc++)=='/' && (*loc=='*' || C_plus_plus && *loc=='/'))
  @< Scan to the end of the comment or of the current line,
     setting |comment_continues| to indicate which case applies;
     in the former case |goto restart|, in the latter |return '\n'| @>
  if (isspace(c))
    if (line_output)
      if (preprocessing) return ' '; @+ else goto restart;
      /* ignore other spaces */
    else return c;
      /* when |line_output==false| preserve white space faithfully */
  if (c=='#' && loc==buffer+1) preprocessing=true;
}

@ When we come to this code |loc| either points to the beginning of a
multi-line comment or to the second character of `\.{//}' or `\.{/*}'; in
the latter case we must advance |loc| before looking for `\.{*/}', lest
`\.{/*/}' would be scanned as a complete comment. Another subtlety arises
when both |preprocessing| and |comment_continues| are raised, which means
that a preprocessor line has ended inside a comment; although in \Cee\ this
means that the preprocessor line continues, \.{\me.} will terminate it
because the comment is removed and a non-escaped newline remains. The
warning message below can be ignored by the user unless something follows
the comment, but the proper way to resolve this is to start the comment on
the line after the preprocessor line, which avoids all confusion.

@< Scan to the end of the comment... @>=
{ boolean one_liner=false;
  if (!comment_continues)
  { if (!line_output) {@; store_byte('/'); store_byte(*loc); }
    one_liner=*loc++=='/';
      /*/ record kind of comment, and advance to comment proper {/}*/
  }
  else if (preprocessing)
  { print("\nWarning: Multi-line comment in preprocessor line");
		    @.Multi-line comment...@>
    preprocessing=false; mark_harmless();
  }
  if (comment_continues=skip_comment(one_liner))
    /* scan to end of comment or newline */
    return '\n'; /* comment contains a newline; |get_line| has been called */
  else goto restart; /* comment complete, get next token */
}

@ The following code assigns values to the combinations `\.{++}', `\.{--}',
`\.{->}', `\.{>=}', `\.{<=}', `\.{==}', `\.{<<}', `\.{>>}', `\.{!=}',
`\.{\v\v}' and `\.{\&\&}'.  The compound assignment operators (e.g., `\.{+=}')
are separate tokens, according to {\sl The C Reference Manual\/} (so there may
be white space and even comments in between), but according to the
\caps{ANSI/ISO} standard they are single tokens. We have not allocated
single-byte codes for them, so we scan them as two separate symbols; these
will be adjacent on output. Therefore, as far as \.{\me.} is concerned, one
may write a separation between the two characters even if the \Cee-compiler
does not accept this; however, this practice should still be avoided because
it will confuse the parser of |CWEAVE|, leading to ill-formatted output.

@d compress(char2,code) @+
  if (*loc==char2) return ++loc,code @;

@< In applicable cases... @>=
case '+': compress('+',plus_plus); break;
case '-': compress('-',minus_minus); compress('>',minus_gt); break;
case '=': compress('=',eq_eq); break;
case '>': compress('=',gt_eq); compress('>',gt_gt); break;
case '<': compress('=',lt_eq); compress('<',lt_lt); break;
case '&': compress('&',and_and); break;
case '|': compress('|',or_or); break;
case '!': compress('=',not_eq); break;

@ Identifiers are looked up unless they consist of a single character
|c<0x80|, in which case that character can be found as |*id_first|; there is
no need to set |id_loc| in either case. In case any 8-bit characters appear in
identifiers, they are only stored in the name table, and cannot interfere with
any of the special codes we use.

@< Get an identifier @>=
{ id_first=--loc; /* mark beginning of identifier */
  do c=*++loc;
  while (c<0x80 ? isalnum(c) || c=='_' : translation_exists(c));
  cur_id=
    loc==id_first+1 && (eight_bits)(*id_first)<0x80
      ? NULL : id_lookup(id_first,loc,0);
}

@ Scanning numeric constants is straightforward. The only subtle point is
that we must not call functions like |isdigit| with an argument of type
|char|, like |*loc|, since that may fail for 8-bits characters if |char| is
a signed type; therefore we copy characters into |c|, which has type
|unsigned char|, before applying such functions. We don't mind \.{8}'s and
\.{9}'s appearing in octal constants (although the \Cee~compiler will), so
there is no need here to distinguish between octal and decimal constants.

@< Get a numeric constant @>=
{ if (*(id_first=loc-1)=='0' && tolower((eight_bits)*loc)=='x')
    /* hex constant */
    do c=*++loc; while (isxdigit(c));
  else /* octal, decimal or float constant */
  { while (isdigit(c)) c=*loc++;
    if (c=='.') @+ do c=*loc++; while (isdigit(c));
    if (tolower(c)=='e') /* floating point constant with exponent */
    { if ((c=*loc)=='+' || c=='-') c=*++loc;
      while (isdigit(c)) c=*++loc;
    }
    else --loc; /* back up to first character after constant */
  }
  while (isalpha(c)) c=*++loc;
    /* incorporate any `\.{U}', `\.{L}', or `\.{F}' suffixes */
  id_loc=loc;
}

@ After an `\.@@' sign has been scanned, the next character tells us whether
there is more work to do. Although control codes like \:\^ and \:t are
ignored by \.{\me.}, their control text has to be skipped, which is
performed by |get_control_text|. Verbatim constructions also use
|get_control_text|, but then return the |verbatim| control code, unless their
control text was empty, in which case they are ignored.

@< Get a control code... @>=
{ eight_bits cc=code_of(*loc++);
  switch(cc)
  { case ignore: goto restart;
    case control_text: get_control_text(); goto restart;
    case verbatim: if (get_control_text()) goto restart; @+ else break;
    case ord: @<Scan an \ASCII. constant@> @+ break;
    case module_name:
    @<Scan the module name and make |cur_mod| point to it@>
  }
  return cc;
}

@ Before we scan a module name, we record whether it started with \:(, so
that we may record it afterwards as an auxiliary file name if necessary.

@<Scan the module name...@>=
{ boolean file_module= loc[-1]=='(';
  /* does this module define an output file? */
  cur_mod=get_module_name();
  if (file_module && cur_mod!=NULL)
     @<Adjoin |cur_mod| to the set of output files@>
}

@ There is no reason why we should allow a newline within an \ASCII.
constant, even if it is escaped.

@<Scan an \ASCII. constant@>=
id_first=loc; /* first character after opening quote */
while (*loc!='\'')
{ if (*loc++=='\\') loc++; /* accept any character following backslash */
  if (loc>=limit) {@; err_print("! ASCII constant didn't end"); break; }
				 @.ASCII constant didn't end@>
}
id_loc=loc++; /* move past closing quote */


@* Scanning a piece of \Cee\ text.
Having acquired the skills of skipping over text and isolating tokens that
really matter, we are now ready to tackle complete pieces of \Cee~text. The
rules for generating the replacement texts for preprocessor directives
following \:d and \:h and for bodies of sections are almost identical; the
only differences are that

\yskip \item{a)}
  Module names are not possible in preprocessor directives; indeed, the
  appearance of a module name terminates such directives and starts the
  \Cee~part of the section.

\item{b)}
  The codes \:d, \:h, \:f, and \:c are not allowed in a \Cee~part, while
  they terminate preprocessor directives.

\item{c)}
  The code \:p is not allowed inside preprocessor directives.
\yskip

\noindent Therefore there is a single function |scan_repl| with a parameter
|context| indicating which of the two kinds of replacement texts is being
scanned; its value is either |preproc_directive| or |section_body|. After
|scan_repl| has acted, |cur_text| will point to the replacement text just
generated, and |next_control| will contain the control code that terminated
the activity.

@d preproc_directive  0
@d section_body 1

@<Global...@>=
text_pointer cur_text; /* replacement text formed by |scan_repl| */
eight_bits next_control; /* control code which has already been scanned */

@ We will try to reduce storage requirements by discarding any trailing
newlines both after preprocessor directives and section bodies. This does not
affect the line sequencing indicated by \&{\#line} directives in either case,
since in the former case such directives are not issued in the first place,
and in the latter case a new \&{\#line} directive will be issued before any
further code is added. Note that we cannot reliably recognise newlines in a
right-to-left motion when the end of a replacement text is reached, since
there is a mixture of one-byte, two-byte and five-byte tokens. Therefore,
instead of tracing backwards, we mark the position of the last non-discardable
token every time one is appended, and back up to the last position recorded
when reaching the end of a replacement text. Given this mechanism, it is easy
to discard other things than newlines, such as \&{\#line} directives that
apply to no lines at all, and to back up at other occasions where the line
sequencing is interrupted as well, like before applied module names.

This is also a good place to add a simple but useful check, namely that in
each section and macro replacement text the parentheses and braces should be
matched (failure in this respect can lead to syntax errors that are very
hard to locate in the \.{CWEB} source, since the real problem may be buried
inside a module name or macro).

@c
void scan_repl (eight_bits context) /* creates a replacement text */
{ eight_bits a; /* the current token */
  eight_bits* keep=tok_ptr; /* non-discardable stuff up to this point */
  int brace_level=0, par_level=0;
  if (context==section_body)
    @< Insert the line number into |tok_mem| @>
  do @< Scan and store a token, setting |keep=tok_ptr| afterwards unless it
	is discardable; when the replacement text has ended |goto done| @>
  while (true);
done: tok_ptr=keep; /* backup over trailing newlines and \&{\#line} marks */
  next_control=a; /* this control code must be reconsidered */
  if (par_level>0 || brace_level>0) @< Report unmatched opening symbols @>
  @< Wrap up the accumulated bytes into a completed text,
     pointed to by |cur_text| @>
}

@ When the function |scan_repl| is called we have either just scanned \:h or
the identifier following \:d or a defining occurrence of a module name (with
the following `\.='); in all cases we start by calling |get_next|.

@< Scan and store a token... @>=
{ switch (a=get_next())
  {
  @\@< Cases where |a| is a special token (|id_code|,
      |module_name|, etc.): either store the corresponding bytes and
      |continue|, or |goto done| if |a| signals the end of this
      replacement text@>
  case format: case definition: case header: case begin_C:
    if (context==preproc_directive) goto done;
    err_print
      ("! `@@f', `@@d', `@@h', and `@@c' are ignored in section body");
		     @.`@@f', `@@d', ... are ignored...@>
    continue;
  case new_section: goto done;
  @\@< Cases that keep track of |par_level| and |brace_level| @>
  case '\n': store_byte('\n');
    if (context==section_body && print_where) /* input file was switched */
    { tok_ptr=keep; /* back up over discardable items */
      @< Insert the line number into |tok_mem| @>
    }
    continue;
  case join_code: a=join;
  }
  store_byte(a); keep=tok_ptr; /* mark as non-discardable */
}

@ We test matching parentheses and braces separately because it is easier,
and in the unlikely case that a mismatch should go undetected the compiler
will surely report an error that is easy to locate. All parentheses and
braces that are returned from |get_next| are real ones, i.e., not part of a
string or comment.

@< Cases that keep track of |par_level| and |brace_level| @>=
case '(': ++par_level; break;
case ')':
  if (par_level<=0)
  {@; err_print("! Unmatched closing parenthesis"); continue; }
		 @.Unmatched closing...@>
  --par_level; break;
case '{': ++brace_level; break;
case '}':
  if (brace_level<=0)
  {@; err_print("! Unmatched closing brace"); continue; }
  --brace_level; break;

@ The most difficult part of matching parentheses and braces is reporting any
errors at the end of the macro or section in proper English. We supply the
missing closing symbols so that there is a slightly larger chance that the
produced \Cee~code will compile. Note that there is no need to update |keep|
here, since no backing up will follow.

@< Report unmatched opening symbols @>=
{ char *p, *s; int l;
  if (par_level>0) l=par_level, s="parenthes", p= l>1 ? "es" : "is";
  else l=brace_level, s="brace", p=l>1 ? "s" : "";
  print("\n! There %s %d unclosed %s%s"
       ,par_level+brace_level>1 ? "are" : "is", l, s, p );
  @.There are unclosed...@>
  if (par_level>0 && brace_level>0)
    print(" and %d unclosed brace%s"
         , brace_level, brace_level>1 ? "s" : "");
  print(" in the previous ");
  err_print(context==preproc_directive ? "macro" : "section");
  while (--par_level>=0) store_byte(')');
  while (--brace_level>=0) store_byte('}');
}

@ When a module reference is inserted the line sequencing is interrupted,
so we must restore it afterwards; on the other hand newlines before the
insertion can be discarded.

@< Cases where |a| is...@>=
case id_code:
  @< Store the identifier just scanned @> @+ keep=tok_ptr; continue;
case module_name: if (context==preproc_directive) goto done;
  if (cur_mod!=NULL) /* don't record bad module name */
  { sixteen_bits n = mod_index(cur_mod); /* index of module name */
    @< If this looks like a defining occurrence, report a runaway section @>
    if (line_output) tok_ptr=keep; /* back up */
    store_two_bytes(0xA800+n); keep=tok_ptr;
      /* store reference to module name */
    @< Insert the line number into |tok_mem| @>
      /* to get in phase after module insertion */
  }
  continue;
case constant: case verbatim:
  @/@< Copy a constant or verbatim construction @> @+ keep=tok_ptr;
  continue;
case ord:
  @<Translate an \ASCII. constant@> @+ keep=tok_ptr; continue;
case include_preproc:
  if (context==preproc_directive)
    err_print("! `@@p' is forbidden in preprocessor directive");
	       @.`@@p' is forbidden...@>
  else
  { if (line_output) tok_ptr=keep;
    store_byte(0xFA); atp_seen=true; keep=tok_ptr;
    @< Insert the line number... @>
  }
  continue;
case char_trans:
  err_print("! `@@l' is only allowed in limbo");
	     @.`@@l' is only allowed in limbo@>
  continue;

@ Identifiers of length~$1$ have no name pointer, but |id_first| still
points at the relevant character.

@< Store the identifier... @>=
if (cur_id==NULL) store_byte(*id_first); /* single-character identifier */
else store_two_bytes(0x8000+id_index(cur_id));

@ Here we look to see whether a module name is followed by~`\.{=}'
or~`\.{+=}' (whether another `\.{=}' follows is irrelevant at this point),
in which case we assume that introduces the \Cee-part of a new section, of
which we somehow missed the opening \:\ . Note that this need not limit the
use of module names standing for single expressions, even if followed by the
`\.{=}' or the `\.{+=}' operator, since such names should be followed by
`\.{@@;@@;}' for proper treatment by \.{CWEAVE}.

@< If this looks like a defining occurrence, report... @>=
{ char *p=loc;
  while (*p==' ' && p<limit) ++p;
  if (*p=='+') ++p; /* an optional `\.+' is allowed */
  if (*p=='=')
    err_print
      ("! Illegal defining occurrence of module name; did you forget `@@ '?");
	@.Illegal defining occurrence...@>
}

@ Constants of length one can only be single digit numbers, so they can
be stored in a single (unquoted) byte like single letter identifiers.

Inside strings (and verbatim constructions) {\it any\/} character may in
principle show up, including the character |verb_quote| used to delimit the
verbatim storage of the string. Even if this is hardly a sign of good
programming, the effect on output would be so disastrous that some counter
measure is called for: we replace |verb_quote| by its three-digit octal
escape code, which the user should have written to begin with. Because a \:i
line or a change file match can occur inside a string, we must test
|print_where| at the end of a string.

@< Copy a constant... @>=
if (id_loc==id_first+1) store_byte(*id_first); /* single digit number */
else
{ store_byte(verb_quote); /* opening */
  do
    if (*id_first==verb_quote)
    { ++id_first; store_byte('\\'); store_byte('0'+(verb_quote>>6));
      store_byte('0'+((verb_quote>>3)&7)); store_byte('0'+(verb_quote&7));
    }
    else
    { if ((eight_bits)(*id_first)>=0x80) store_byte(0xF8);
	/* quote 8-bit characters */
      store_byte(*id_first++);
    }
  while (id_first<id_loc);
  store_byte(verb_quote); /* closing */
  if (context==section_body && print_where)
    @< Insert the line number into |tok_mem| @>
}

@ This section should be rewritten on machines that don't use \ASCII.
code internally. It is probably best to use a table indexed by characters
and containing \ASCII. codes for this purpose (like the |xord| array in
\TeX). Incidentally, it would be tempting to write statements like
`|case'\\': c=@'\\';|' here, which would in principle be correct in any
other place, but which would cause circularity in a bootstrapped system like
this.
@^ASCII code dependencies@>

@< Translate an \ASCII. constant @>=
{ int c=*id_first++;
  if (c=='@@' && *id_first++ !='@@')
  @/{@; err_print("! Double `@@' should be used in ASCII constant");
		   @.Double `@@' should be used...@>
      --id_first;
  }
  if (c=='\\')
  { c=*id_first++;
    switch (c)
    { case 't':c='\011';break;
      case 'n':c='\012';break;
      case 'b':c='\010';break;
      case 'f':c='\014';break;
      case 'v':c='\013';break;
      case 'r':c='\015';break;
      case '0':c='\0';break;
      case '\\':c='\134';break;
      case '\'':c='\047'; break;
      case '\"':c='\042'; break;
      default: err_print("! Unrecognised escape sequence");
			  @.Unrecognised escape sequence@>
    }
  }
  else
  { /* at this point |c| should be converted to its \ASCII. code number */
  @+}
  if (id_first!=id_loc)
    if (id_loc>id_first)
      err_print("! ASCII constant should be single character");
		 @.ASCII constant should be...@>
    else {@; err_print("! Empty ASCII constant"); c=0; }
			@.Empty ASCII constant@>
  store_byte(verb_quote);
  if (c>=100) store_byte('0'+c/100); @+ if (c>=10) store_byte('0'+(c/10)%10);
  store_byte('0'+c%10); @/
  store_byte(verb_quote);
}

@ Line number marks, which will lead to \&{\#line} directives on output, are
inserted in several places, so we shall call a function to do the necessary
work. For reasons explained below this function scans ahead a bit through
the input, and in doing so it may hit the end of the input; when this
happens it must be reported back to |scan_repl|. Therefore |mark_line|
returns a boolean value which tells whether |input_has_ended|.

@< Prototypes @>= boolean mark_line(void);

@~When |mark_line| returns |true|, we should terminate |scan_repl| with
|next_control==new_section|; this is achieved by jumping to |done|, where
some exitialisations are done, and where |next_control| is assigned from the
variable~|a|.

@< Insert the line... @>=
{@; if (mark_line()) {@; a=new_section; goto done; } }

@~If |line_output==true|, the following items are stored for recording the
line number: first a byte equal to |0xF9|, then the numeric line number,
and finally the |id_index| of the file name, which is looked up as if it were
an identifier (except that we use the feature of |id_lookup| to accept
null-terminated strings as first argument if the second argument is~|NULL|).
Note that |keep| is not increased; if nothing essential follows this line
number indication then it need not be recorded. We improve compactness of
storage and output a bit by discarding any newlines that are immediately
ahead; the inserted line number will automatically be increased properly.
While reading we might find there is no input left, in which case no line
number is recorded and |true| is returned.

@c
boolean mark_line(void)
{ while (loc>=limit)
    @+ if (!get_line()) return true; /* skip over newlines */
  print_where=false; /* request is being serviced */
  if (line_output)
  { store_byte(0xF9);@/
    id_first=changing ? change.name : cur_file_name;
    store_two_bytes(changing ?  change_line : cur_line);
    store_two_bytes(id_index(id_lookup(id_first,NULL,0)));
  }
  return false;
}


@* Scanning a section.
The function |scan_section| starts when \:\ , \:\~, or \:* has been sensed
in the input, and it proceeds until the end of that section.  It uses
|section_count| to keep track of the current section number; hopefully
\.{CWEAVE} and \.{\me.} will both assign the same numbers to sections.

@c
void scan_section (void)
{ ++section_count;
  if (loc[-1]=='*') print_section_progress();
  @<Scan the \TeX\ and definition parts of the current section@>
  @<Scan the \Cee\ part of the current section@>
}

@ For \.{\me.} there is no real difference between the \TeX\ and
definition parts of a section; after all, the latter can contain \&{format}
commands that are ignored just like \TeX~text is. So we set up a loop that
stops when |next_control>=begin_C|, and handles the cases |definition| and
|header|, while ignoring anything else. The price we pay for our haste in
skipping the \TeX~part of a section (and \:f commands) is that we may have
to back up after finding \:< in order to scan the whole module name and set
|cur_mod|, which requires calling~|get_next|. In case \:< follows a
|definition| or |header| code, it is scanned by |scan_repl|, and no backing
up is required; therefore the code for backing up and rescanning cannot be
moved outside the loop.

@<Scan the \TeX\ and definition part...@>=
{ next_control=ignore; /* clear |new_section| */
  do
  { if (next_control<definition)
    { if ((next_control=skip_ahead())==module_name)
      @/{@; loc -= 2; get_next(); }
	/* back up and scan the module name itself */
    }
    else if (next_control==definition)
    { @< Scan and store the identifier after \:d, if not found |continue| @>
      scan_repl(preproc_directive);
      cur_text->text_link=macro_flag; /* characterise as macro */
    }
    else /* |next_control==header| */
    {@; scan_repl(preproc_directive); cur_text->text_link=header_flag; }
    if (next_control==module_name)
    @< If this is not a defining occurrence, set |next_control=ignore| @>
  } while (next_control<begin_C);
}

@ We retain the strange lexical rule of the \Cee\ preprocessor that a space
after the identifier defined as a macro is quite significant, especially
when followed by `\.{(}'. Together with the spaces in preprocessor lines
contained in a \Cee~part, these are the only kind of unquoted spaces that
are ever stored as tokens, when |line_output| holds.

@< Scan and store the id... @>=
{ do next_control=get_next(); /* allow white space before identifier */
  while (line_output ? next_control=='\n'
        : next_control<0x80 && isspace(next_control));
  if (next_control!=id_code)
  @/{@; err_print("! Macro definition ignored, must start with identifier");
		   @.Macro definition ignored...@>
    continue;
  }
  @< Store the id... @>
  if (isspace((eight_bits)*loc)) store_byte(' ');
      /* separate the replacement text from a parameterless macro */
}

@ Due to the fact that module names are allowed enclosed in `\pb', it is not
certain that the occurrence of a module name signals the start of the
\Cee~part of the section. The decisive factor will be what follows the
module name: if this is `\.=' or `\.{==}', possibly preceded by `\.+', then
this is really the name of the module defined in this section. If not, then
we require that the next token is a `\.\v' (the one closing `\pb'), unless
we are in compatibility mode, in which case anything goes. We use the fact
that scanning the `\.+' and `\.=' symbols will not destroy the value in
|cur_mod|.

@< If this is not a defining occurrence, set |next_control=ignore| @>=
{ eight_bits t=get_next();
  if (t=='+') t=get_next(); /* skip an optional `\.+' */
  if (t!='=' && t!=eq_eq) /* then apparently a cited occurrence, ignore it */
  { next_control=ignore;
    if (t!='|' && !compatibility_mode)
      err_print("! `=' sign missing, module name ignored");
		 @.`=' sign missing, module name ignored@>
  }
}

@ If we come to this module with |next_control==module_name|, then the code
of the previous section has just removed the `\.=' or `\.{==}' and optional
`\.+' following the module name.

@<Scan the \Cee...@>=
{ mod_pointer p; /* module name for the current section */
  switch (next_control)
  { default: return; /* no \Cee\ part present */
    case begin_C: p=NULL; break; /* section contributing to unnamed module */
    case module_name: p=cur_mod; /* section contributing to named module */
  }
  store_two_bytes(0xD000+section_count); /* record section number */
  scan_repl(section_body);
      /* now |cur_text| points to the replacement text */
  @< Link the section body at |cur_text| to the module named |p|,
     or to the unnamed module if |p==NULL| @>
}

@ Sections whose body starts with a bad name (e.g., an ambiguous prefix) will
be treated as if they were unnamed, so the user had better not ignore the
error message. Unnamed sections get a slightly special treatment, which makes
accumulating many of them a bit more efficient than building up a named module
that is defined in many sections.

@<Link the section...@>=
{ if (p == NULL) /* unnamed section, or bad module name */
  {@; *last_unnamed = cur_text; last_unnamed = &cur_text->text_link; }
  else
  { text_pointer* q=&p->equiv; /* text for the current module */
    while(*q!=NULL) q=&(*q)->text_link; /* find end of list */
    *q=cur_text; /* add section to module */
  }
  cur_text->text_link=NULL;
  /* end list, also marking replacement text as a non-macro */
}


@* Phase I processing.
Finally we can wrap everything up and define the global structure
of~Phase~I.

@<Proto...@>=
void phase_one (void);
  /* read all the user's text and compress it into |tok_mem| */

@~The only thing left to do is picking up character translations in the limbo
part.

@c
void phase_one (void)
{ phase=1; section_count=0; reset_input();
  while ((next_control=skip_ahead())!=new_section)
    if (next_control==char_trans) @< Store a character translation @>
  while (!input_has_ended) scan_section(); /* load all sections */
  check_complete(); /* verify that change file hasn't got out of sync */
}

@ The code \:l should be followed (after optional space) by a two-digit hex
number~|c| with |c>=0x80|, white space, and a translation string of at most
|trans_limit| characters, either alphanumeric or underscores, terminated by
another space.

@< Store a character translation @>=
{ int c;
  while (loc<limit && isspace((eight_bits)*loc)) ++loc; /* skip space */
  if (!(isxdigit((eight_bits)loc[0])&&
	isxdigit((eight_bits)loc[1])&&
	isspace ((eight_bits)loc[2])))
    err_print("! Two-digit hex number and space should follow `@@l'");
	       @.Two-digit hex number...@>
  else if (sscanf(loc,"%x",&c),c<0x80)
    err_print("! You cannot translate characters < 0x80");
	       @.You cannot translate...@>
  else
  { char* p=trans_of(c); int i=0;
    loc+=3; @+
    while (find_char() && isspace((eight_bits)*loc)) ++loc; /* skip space */
    if (!input_has_ended)
      while (isalnum((eight_bits)*loc) || *loc=='_')
	if (++i<=trans_limit) *p++=*loc++; @+ else break;
    if (i>0) *p='\0'; /* terminate translation unless none was given */
    @< Report any problems with the translation string @>
  }
}

@ The main contribution of the \:l feature to \.{CWEB} is new syntax,
leading to a variety of new error messages. It is possible to give an empty
translation for a character by terminating the translation string by a
non-space character that is not allowed in identifiers; this is forbidden
however, since an empty translation for~|c| would just make
|translation_exists(c)| fail, making~|c| useless anyway. Even if the
translation string is not empty, terminating it by a non-space character is
considered an error, since this would make the \TeX\ macro used for
formatting the corresponding line of the printed document fail.

@< Report any problems... @>=
if (i==0) err_print("! Translation string absent after `@@l'");
		     @.Translation string...@>
else if (i>trans_limit) err_print("! Translation string too long");
else if (!isspace((eight_bits)*loc))
  err_print("! Translation string not terminated by space");


@ At the end of the run, if |STAT| was defined and the `\.{+s}' flag
present, we report how much of all the arrays was actually needed.

@d report(k,c,m)
  printf("%lu %ss (out of %lu)\n",(unsigned long)(c),k,(unsigned long)(m))

@c
#ifdef STAT
void print_stats (void)
{ print("\nMemory usage statistics:\n");
@/report("identifier", id_index(id_ptr), max_idents);
@/report("module name", mod_index(mod_ptr), max_modules);
@/report("byte", byte_ptr-byte_mem, max_bytes);
@/report("replacement text", text_ptr-text_table, max_texts);
@/report("token", tok_ptr-tok_mem, max_toks);
}
#endif


@* Index.  Here is a cross-reference table for the \.{\me.}
processor.  All sections in which an identifier is used are listed with
that identifier, except that reserved words are indexed only when they
appear in format definitions, and the appearances of identifiers in
module names are not indexed. Underlined entries correspond to where
the identifier was declared. Error messages and a few other things
like ``\ASCII. code dependencies'' are indexed here too.
