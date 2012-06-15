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

\def\currentCWEBxversion{x3.51}

\def\pb{$\.|\ldots\.|$} % C brackets (|...|)
\def\LKC.{Levy/Knuth \.{CWEB}}
\def\:#1{`\.{@@#1}'}

\def\title{Common code for CTANGLE and CWEAVE (Version \currentCWEBxversion)}
\def\topofcontents
{\topglue 0pt plus .5 fill
 \centerline{\titlefont Common code for {\ttitlefont CTANGLE} and
   {\ttitlefont CWEAVE}}
 \vskip 15pt
 \centerline{(\.{CWEB} version \currentCWEBxversion)}
}
\def\botofcontents
{\vfill\noindent
 Copyright \copyright\ 1987,\thinspace1990 Silvio Levy and Donald E. Knuth
 \par\noindent
 Copyright 1994--2009 Marc A. A. van Leeuwen
 \bigskip\noindent
 Permission is granted to make and distribute verbatim copies of this
 document provided that the copyright notice and this permission notice
 are preserved on all copies.

 \smallskip\noindent
 Permission is granted to copy and distribute modified versions of this
 document under the conditions for verbatim copying, provided that the
 entire resulting derived work is distributed under the terms of a
 permission notice identical to this one.
}

@* Introduction. This file contains code common to both the \.{CTANGLE}
and \.{CWEAVE} programs, contained in a separate compilation unit. It roughly
concerns the following problems: input routines, name table handling, error
handling and handling of the command line.

@i common.inc


@* Generalities.
That completes the contents of \.{common.inc}.
In the texts below we will sometimes use \.{CWEB} to refer to either
of the two component programs, if no confusion can arise.

Here is the overall appearance of this file, except for the function
definitions for (local and public) functions, which will follow in the
remaining unnamed sections.

@c
@< Function prototypes used but not defined in the shared code @>@;
@< Definitions of variables common to \.{CTANGLE} and \.{CWEAVE} @>@;
@< Prototypes of local functions @>@;

@ For all functions and variables defined here that are accessible to
\.{CTANGLE} and \.{CWEAVE}, prototype respectively |extern| declarations are
placed on the file \.{common.h} that is included by all three source files.
Typedef declarations that are publicly visible also appear in this file.

@( common.h @>=
@< Public typedef declarations @>@;
@< Declarations of public variables and function prototypes @>@;

@ Since a number of large arrays are used to store various kinds of data, we
wish to have some control over the number of bytes occupied by the basic
items. In such cases we use one of the following types rather than |int|.

@< Public typedef declarations @>=
typedef char boolean;
typedef unsigned char eight_bits;
typedef unsigned short sixteen_bits;

@ In certain cases \.{CTANGLE} and \.{CWEAVE} should do almost, but not
quite, the same thing.	In these cases we've written common code for both,
differentiating between the two by means of the global variable |program|.
Furthermore, |CTANGLE| operates in two phases (input and output), and
similarly |CWEAVE| operates in three phases (cross-reference collection,
translation of the source, and output of the index); the global variable
|phase| tells which phase we are in.

@< Declarations... @>=
extern int program, phase;

@~@< Definitions... @>=
int program, phase;

@ There's an initialisation function that gets both \.{CTANGLE} and
\.{CWEAVE} off to a good start.

@< Declarations...@>=
void common_init (int argc,char** argv,char* version);

@~We will fill in the details of this function later.
@c
void common_init (int argc,char** argv,char* version)
{ @< If called with argument \.{--version}, print |version| and quit @>
  @< Initialise variables @>
  @< Set the default options common to \.{CTANGLE} and \.{CWEAVE} @>
  scan_args(argc,argv);
}

@ In case either \.{CTANGLE} or \.{CWEAVE} is called with a single argument
equal to \.{--version}, we just print the |version| (which includes the version
number) and call |exit|.

@< If called with argument \.{--version}, print |version| and quit @>=
if (argc==2 && strcmp(argv[1],"--version")==0)
  {@; print("%s\n",version); exit(0); }

@* Input routines.
The lowest level of input to the \.{CWEB} programs is performed by
|input_ln|, which must be told which file to read from. The return value
of |input_ln| is |true| if the read is successful and |false| if not
(i.e., if file has ended). The conventions of \TeX\ are followed; i.e.,
the characters of the next line of the file are copied into the |buffer|
array, and the global variable |limit| will point to the first unoccupied
position; trailing white space is ignored.

@< Declarations... @>=
extern char buffer[], *loc, *limit;

@~The value of |limit| must be less than or equal to |buffer_end|, so that
|*buffer_end| is never filled by |input_ln|. The characters |*limit| and
|limit[1]| are reserved for placing a sentinel at the end of the line. For
the convenience of |CWEAVE|, the buffer is extended by |longest_name|
characters, so that there is enough space to place any module name after
the input.

@d buffer_end  (&buffer[buf_size-2]) /* practical end of |buffer| */

@< Definitions... @>=
char buffer[long_buf_size]; /* where each line of input goes */
char *loc=buffer;
    /* points to the next character to be read from the buffer */
char *limit=buffer; /* points to the end of the input line */

@ If a non-empty line follows the last newline in a file, we return it with
a success status when |EOF| is read; in this case the next call will read
|EOF| once again, and that time return failure. We are careful in case
|isspace| is a macro using its argument more than once. As a service to the
include and change file handling functions, |get_line| replaces any initial
\:I, \:X, \:Y, or \:Z by its lower-case counterpart. The value of |loc| is
usually irrelevant during the lower level input functions, and will later be
set properly by |get_line|; however, when error messages are given, the
value of |loc| will determine the way that the current line of input is
displayed, so we give it an appropriate value in those cases.

@c
local boolean input_ln (FILE *f)
  /* copies a line into |buffer| or returns |false| */
{ register int c; /* the character read */
  register char* k=limit=buffer;  /* where next character goes */
  while ((c=getc(f))!='\n' && c!=EOF)
    if (k<=buffer_end) {@; *k++=c; @+ if (!isspace(c)) limit=k; }
  if (k>buffer_end)
  { loc=&buffer[0]; /* now |err_print| will display unbroken input line */
    err_print ("! Input line too long");   @.Input line too long@>
    if (limit>buffer_end) limit=buffer_end; /* truncate line */
  }
  if (buffer[0]=='@@' && limit>&buffer[1] && strchr("IXYZ",buffer[1])!=NULL)
    buffer[1]=tolower(buffer[1]);
  return c!=EOF || limit>buffer; /* whether anything new has been found */
}

@ Now comes the problem of deciding which file to read from next.  Recall
that the actual text that \.{CWEB} should process comes from two streams: a
|web_file|, which can contain possibly nested include commands \:i,
and a |change_file|.
For each file we store its name and line number for error reporting and for
the production of \&{\#line} directives by |CTANGLE|. Information for the
|web_file| together with the currently open include files is kept on a stack
|file| with stack pointer |include_depth|, while the change file has its own
record. The boolean |changing| tells whether or not we're reading from the
|change_file|. Whenever we switch from the |cur_file| to the |change_file|
or vice versa, or if the |cur_file| has changed, we tell |CTANGLE| to print
this information by means of a \&{\#line} directive in the \Cee\ file, by
raising the |print_where| flag. This flag is handled like an interrupt
request, i.e., it remains raised until it is serviced at an appropriate time
(when a complete token has been scanned), whereupon it is cleared by the
service routine. In Phase~I of \.{CWEAVE} header files following \:h
are swiftly scanned for typedef declarations, which creates one or more
extra levels of input, for which we use the same stack as for~\:i.
Certain functions operate differently at such times, so we must be able
to distinguish these inclusions; this is achieved by maintaining a boolean
variable |including_header_file|.

@< Declarations...@>=
#define max_file_name_length 256
extern struct f
{ FILE *file; char name[max_file_name_length]; sixteen_bits line; }
file[], change;
extern int include_depth;
extern boolean input_has_ended, changing, web_file_open, print_where
	     , including_header_file;
@)
boolean locate_file_name();
boolean push_input_file(boolean,boolean); /* start a new level of input */
boolean get_line (void); /* get the next line of merged input */

#define cur_file file[include_depth].file /* current file */
#define cur_file_name file[include_depth].name /* current file name */
#define cur_line file[include_depth].line
  /* number of current line in current file */
#define web_file file[0].file
#define change_file change.file
#define change_line change.line

@ We also keep an array |at_h_path| of alternative search paths for locating
\:h files, and a path |at_i_path| for locating \:i files.

@d max_include_depth 32
  /* maximum nesting depth of source files, not counting the change file */
@d max_include_paths 16
  /* maximum number of additional search paths for \:h files */
@d max_path_length 200 /* maximal length of a search path */

@<Definitions...@>=
struct f file[max_include_depth]; /* stack of non-change files */
struct f change; /* change file */
local char web_file_name[max_file_name_length]
         , change_file_name[max_file_name_length]
         , alt_web_file_name[max_file_name_length];
int include_depth; /* current level of nesting */
boolean input_has_ended; /* whether there is no more input */
boolean changing; /* whether the current line is from |change_file| */
boolean web_file_open=false; /* whether the web file is being read */
boolean print_where=false; /* should |CTANGLE| print line and file info? */
local struct { char* name; int length; }
  at_h_path[max_include_paths],at_i_path;
  /* alternative search paths for \:h and \:i */
boolean including_header_file=false; /* are we processing \:h? */

@ Before we consider merging the main input with the change file, we first
handle pushing and popping files in the main stream. We first define a
function to delimit the name in the input buffer, leaving the caller the
possibility to inspect it before proceeding to |push_input_file| below. The
function |locate_file_name| points |id_first| at the first character and
|id_loc| beyond the last character of the file name, signaling a possible
absence of the file name altogether; it also sets |loc| to the end of the line
since we shall not read anything more than the file name from the current
line, regardless of whether a file will actually be opened. When the include
file name is delimited by spaces, any white-space character will end the name;
in any case we accept the end of the line as end of the file name.

@c
boolean locate_file_name()
{ char delim=' '; /* the character being used to delimit the file name */
  while (loc<limit && (isspace((eight_bits)*loc))) ++loc;
  if (*loc=='"') delim=*loc++; /* file name in quotes */
  else if (*loc=='<') delim='>',++loc; /* file name in angle brackets */
  if (loc>=limit)
  {@; err_print("! Include file name not given");
	       @.Include file name not given@>
    return false;
  }
  id_first = loc;
  while (loc<limit &&(delim==' ' ? !isspace((eight_bits)*loc) : *loc!=delim))
    ++loc;
  id_loc=loc;
  loc=&limit[1]; /* force |input_ln| before next character is read */
  return true;
}

@ The function |push_input_file| opens a new level of input, e.g., when a \:i
line is found. It is used both from within the common code, namely to open
files included by~\:i, and by a direct call from |CWEAVE|, to open files
included by~\:h (and nested include files) during its first phase. When it is
called, the file name is supposed to start at the first non-blank character
from |loc|; it is delimited by either blank space, double quotes or angle
brackets. Once the file name is located and the file is opened, any further
input from the current line will be discarded. This function has two boolean
parameters, the first telling whether the file to be included is a header file
(rather than a \:i file), the second telling whether changes should be
suspended during the file inclusion; the boolean result tells whether a file
was actually opened.

@c

boolean push_input_file(boolean header,boolean suspend)
{ boolean success=false; /* whether a file has been opened */
  boolean non_system = *id_loc!='>';
    /* true unless a system header is asked for */
 if (++include_depth>=max_include_depth)
  { --include_depth;
    err_print("! Too many nested includes");
@.Too many nested includes@>
    print(" (%d)",max_include_depth);
  }
  else
  { @< Read file name into |cur_file_name| @>
    if (non_system && (cur_file=fopen(cur_file_name,"r"))!=NULL)
      success=true;
    else @< Try to open file in alternative directory, and set |success| @>
    if (success)
    {@; cur_line=0; print_where=true;
      @< If necessary deactivate the change file @>
    }
    else
    { --include_depth;
      if (non_system) /* don't complain about system header files */
	err_print("! Cannot open include file");
		   @.Cannot open include file@>
    }
  }
  return success;
}

@ Here we just do a string copy, bounded by the length of the destination.

@< Read file name into |cur_file_name| @>=
{ char* k=cur_file_name;
  while (id_first<id_loc)
    if (k==&cur_file_name[max_file_name_length-1])
    @/{@; err_print("! Include file name truncated"); break; }
		     @.Include file name truncated@>
    else *k++=*id_first++;
  *k='\0';
}

@ At initialisation time, paths may have been stored in |at_i_path| and
|at_h_path|, that will be prefixed to given file name in an attempt to find
include files. For a file opened in this manner |cur_file_name| will not
contain the full path name but just the final component; this affects error
messages and \&{\#line} directives produced while reading the file. The only
problem that this might cause is that a debugger could be unable to locate
a source file for code included from a \:i file; however, it is not likely
that such a file outside the current directory will contain any code
producing material, and good debuggers have their own means to specify a
search path for source files.

@< Try to open file in alternative directory... @>=
{ char name_buf[max_path_length+max_file_name_length]; int i;
  if (header) /* \:h include file, or subsidiary \.{\#include} */
    for (i=0; i<max_include_paths; ++i)
      if (at_h_path[i].name==NULL) break;
      else
      { strcpy(name_buf,at_h_path[i].name);
	strcpy(&name_buf[at_h_path[i].length],cur_file_name);
        if ((cur_file=fopen(name_buf,"r"))!=NULL) {@; success=true; break; }
      }
  else if (at_i_path.name!=NULL) /* \:i include file */
  { strcpy(name_buf,at_i_path.name);
    strcpy(&name_buf[at_i_path.length],cur_file_name);
    success= (cur_file=fopen(name_buf,"r"))!=NULL;
  }
}

@ We allow values for |at_h_path[0]| and |at_i_path| to be built into the
program by defining the preprocessor symbols |CWEBHEADERS| respectively
|CWEBINPUTS| to be the appropriate strings; these should be complete
prefixes that can be attached to file names, so they probably have to end
with a pathname separator. For |CWEBHEADERS| one could take the system
include file area, although this should usually not be necessary since the
typedefs contained in system header files are already built into |CWEAVE|.
Defining |CWEBINPUTS| only makes sense if \.{CWEB} is compiled specifically
for one particular project, as long as there are no general purpose
\:i~files.  Following \LKC. we also allow |CWEBINPUTS|
to be overridden by an environment variable. For specifying further values of
|at_h_path| we use the more practical method of command line arguments (just
like most \Cee~compilers use for specifying alternative include
directories), which are handled by |scan_args|.

@< Init... @>=
{ char* cwebinputs=getenv("CWEBINPUTS");
  at_h_path[0].name=at_i_path.name=NULL; /* defaults */
#ifdef CWEBHEADERS
  at_h_path[0].name=CWEBHEADERS;
  at_h_path[0].length=(int)strlen(CWEBHEADERS);
#endif
  if (cwebinputs!=NULL)
  { at_i_path.length=(int)strlen(cwebinputs);
    at_i_path.name=strcpy(byte_ptr,cwebinputs);
    byte_ptr+=at_i_path.length+1;
  }
  else
  {
#ifdef CWEBINPUTS
    at_i_path.name=CWEBINPUTS; at_i_path.length=(int)strlen(CWEBINPUTS);
#endif
  }
}

@ In some cases we wish to suspend changes during the inclusion of a file,
and any nested inclusions. This happens when during the first pass of
\.{CWEAVE} a header file specified after \:h is being read in, or when
in compatibility mode a \:i inclusion is issued under control of the
change file. In such cases we suspend the change file by setting
|changing=false| and |change_limit=change_buffer| as if the change file had
ended, after having saved the old values of |changing| and |change_limit|.
The change file will be reactivated when the included file ends, which
requires saving |include_depth| as well.

@<Definitions...@>=
local boolean saved_changing; /* were we changing before it was suspended? */
local char* saved_change_limit; /* end of suspended change line */
local int saved_include_depth=0; /* depth after opening \:h file */

@~Although this could be easily changed, the current code relies on the fact
that suspension cannot be nested (as would be the case if \:h were
activated in a file included under control of the change file), since
compatibility mode does not support \:h file inclusion.

@< If necessary deactivate the change file @>=
if (suspend)
{ saved_changing=changing; changing=false;
@/saved_change_limit=change_limit; change_limit=change_buffer;
@/saved_include_depth=include_depth;
}

@ The function |get_web_line| fetches the next line from the main input
stream, taking care of the interpretation on \:i and of restoring input to
the parent file when an include file has ended. Like |input_ln| it returns a
boolean value telling whether a line could be found. When this is not the
case, it means that either the main input stream has dried up, or we have
come to the end of a header file included by an \:h code issued from the
change file, in which case |changing==true| afterwards, and we should not
read from the main input stream after all.

In compatibility mode |get_web_line| operates differently, since in
\LKC. include files are expanded after matching against the change file
rather than before, although an expanded include file will be reconsidered
for matching against (the same line of) the change file, unless it was
included under control of the change file itself. This means that in
compatibility mode no include file should be opened while preparing the
main input stream for a match. On the other hand, if an include file ends,
there is no other sensible action but to close the file and read on at the
previous level.

@c
local boolean get_web_line(void)
{ do
    if (++cur_line,input_ln(cur_file)) /* then a line has been found */
      if (!compatibility_mode
       && limit>&buffer[1] && buffer[0]=='@@' && buffer[1]=='i')
      { loc=&buffer[2]; print_where=true;
        if (locate_file_name()) /* expand \:i */
         push_input_file(false,false);
      }
      else return true; /* return the line without further action */
    else if (include_depth==0) /* then end of input has been reached */
    @/{@; input_has_ended=true; web_file_open=false; return false; }
    else
    { fclose(cur_file); print_where=true;
      if (include_depth--==saved_include_depth) /* then restore |changing| */
      { changing=saved_changing; change_limit=saved_change_limit;
	saved_include_depth=0; including_header_file=false;
	if (changing) return false; /* fall back into change file */
      }
    }
  while (true);
}

@ Now we come to merging the main input stream with the change file.
When |changing| is false, the first non-empty line of |change_file| after
the next \:x to be matched is kept in |change_buffer|, for purposes of
comparison with the next line of |cur_file|. The test |lines_match| is
used for equality between the two lines; it will never return |true| when
|change_limit==change_buffer| because it will not be invoked when
|limit==buffer|.

@d lines_match()
  (change_limit-change_buffer==limit-buffer
  && strncmp(buffer, change_buffer, limit-buffer)==0)

@<Definitions...@>=
local char change_buffer[buf_size]; /* next line of |change_file| */
local char *change_limit; /* points to the effective end of |change_buffer| */

@ The function |prime_the_change_buffer| sets |change_buffer| in preparation
for the next matching operation. After the change file has been completely
input, we set |change_limit=change_buffer|, so that no further matches will
be made; since blank lines in the change file are not used for matching, we
have |(change_limit==change_buffer && !changing)| if and only if the change
file is exhausted (or suspended). This function is called only when
|changing| is true; hence error messages will be reported correctly.

@c
local void prime_the_change_buffer (void)
{ change_limit=change_buffer;
    /* this value is used if the change file ends */
  @<Skip over comment lines in the change file; |return| if end of file@>
  @<Skip to the next non-blank line; |return| if end of file@>
  @<Move |buffer| and |limit| to |change_buffer| and |change_limit|@>
}

@ While looking for a line that begins with \:x in the change file, we allow
lines that begin with `\.{@@}', as long as they don't begin with \:y or \:z
(which would probably indicate that the change file is fouled up).

@<Skip over comment lines in the change file...@>=
do
{ if (++change_line,!input_ln(change_file)) return;
  if (limit>&buffer[1] && buffer[0]=='@@')
    if (buffer[1]=='x') break;
    else if (buffer[1]=='y' || buffer[1]=='z')
    { loc=&buffer[2]; /* point out error after \:y or \:z */
      err_print ("! Where is the matching @@x?"); @.Where is the match...@>
    }
    else @< Check for erroneous \:i @>
} while (true);

@ When not in compatibility mode, \:i lines are expanded by |get_web_line|,
so it makes no sense to place such a code between \:x and~\:y or between \:y
and~\:z; to allow them outside of the changes they would only cause
confusion by suggesting the inclusion of a subsidiary change file. Therefore
we normally do not allow \:i at the beginning of any line of the change file;
in compatibility mode however, \:i can be used in both sides of a change,
and outside of the changes the code is allowed but ignored.

@< Check for erron... @>=
{ if (buffer[1]=='i' && !compatibility_mode)
  @/{@; loc=&buffer[2]; err_print ("! No includes allowed in change file"); }
}				    @.No includes allowed...@>

@ After a \:x has been found, we ignore the rest of the line, as well
as any blank lines that follow it; since |input_ln| removes trailing blanks,
we can simply test for empty lines.

@< Skip to the next non-blank line... @>=
do
  if (++change_line,!input_ln(change_file))
  @/{@; loc=&buffer[0]; err_print("! Change file ended after @@x"); return; }
while (limit==buffer); @.Change file ended...@>

@ @<Move |buffer| and |limit| to |change_buffer| and |change_limit|@>=
{@; int n=(int)(limit-buffer); change_limit=change_buffer+n;
  strncpy(change_buffer,buffer,n);
}

@ The function |check_change| is used to see if the next change entry should
go into effect; it is called only when |changing| is false. The idea is to
test whether or not the current contents of |buffer| matches the current
contents of |change_buffer|. If not, there's nothing more to do, but if so,
a change is called for. When this happens, all of the text down to the \:y
is supposed to match, and an error message is issued if any discrepancy is
found; after finding \:y we have |changing==true|, so that subsequent lines
will be read from the change file. Since |check_change| is called only when
|change_limit>change_buffer|, i.e., when the change file is active, we don't
have to consider the case here that |get_web_line| returns |false| after
reactivating the suspended change file.

@c
local void check_change (void)
  /* switches to |change_file| if the buffers match */
{ int n=0; /* the number of discrepancies found */
  if (!lines_match()) return;
  print_where=true; /* indicate interrupted line sequencing */
  do
  { changing=true;
    @< Read a line from the change file into the change buffer;
       if \:y is found, |break|; if the change file ends, |return| @>
    changing=false;
    if (!get_web_line())
    @/{@; loc=&buffer[0];
      err_print("! CWEB file ended during a change"); return;
		 @.CWEB file ended...@>
    }
    if (!lines_match()) ++n;
  } while (true);
  if (n>0)
  { loc=&buffer[2];
    print("\n! Hmm... %d of the preceding lines failed to match",n);
	    @.Hmm... $n$ of the preceding...@>
    err_print("");
  }
}

@ Since we read a line from the change file before reading the line from the
main input stream that should match it, we can use |input_ln| to read the
line into |buffer| first, and move it to |change_buffer| afterwards.
When expecting \:y, we signal and ignore any \:x or \:z.

@< Read a line from the change file into the change buffer... @>=
{ if (++change_line,!input_ln(change_file))
  { loc=&buffer[0]; err_print("! Change file ended before @@y");
			       @.Change file ended...@>
    change_limit=change_buffer; changing=false; return;
  }
  if (limit>&buffer[1] && buffer[0]=='@@')
    if (buffer[1]=='y') break;
    else if (buffer[1]=='x' || buffer[1]=='z')
    @/{@; loc=&buffer[2]; err_print("! Where is the matching @@y?"); }
				     @.Where is the match...@>
    else @< Check for erron... @>
  @< Move |buffer| and |limit|... @>
}

@ The function |reset_input|, which gets \.{CWEB} ready to read the \.{CWEB}
source file(s), is used at the beginning of Phase~I of |CTANGLE|, and
of Phases I~and~II of |CWEAVE|.

@< Declarations... @>=
void reset_input (void);

@~Although |reset_input| will not read anything from |web_file| after
opening it, it will move up to the first change line in |change_file|.

@c
void reset_input (void) /* initialise to read the web file and change file */
{ boolean use_change_file= change_file_name[0]!='\0';
  @<Open input files@>
  cur_line=0; change_line=0; include_depth=0;
  if (use_change_file) {@; changing=true; prime_the_change_buffer(); }
    /* prepare change file */
  else change_limit=change_buffer;
    /* emulate that change file that has ended */
  limit=buffer; loc=&buffer[1]; /* now |find_char()| will read a line */
  changing=false; input_has_ended=false;
}

@ The following code opens the input files. We complain about a missing
change file only if it was explicitly mentioned as a command line argument.

@<Open input files@>=
{ if ((web_file=fopen(web_file_name,"r"))!=NULL)
    strcpy(file[0].name,web_file_name);
  else if ((web_file=fopen(alt_web_file_name,"r"))!=NULL)
    strcpy(file[0].name,alt_web_file_name);
  else fatal("! Cannot open \"%s\" as input file", web_file_name);
	      @.Cannot open input file@>
  web_file_open=true;
  if (use_change_file)
    if ((change_file=fopen(change_file_name,"r"))!=NULL)
      strcpy(change.name,change_file_name);
    else if (!change_file_explicit)
      use_change_file=false; /* forget about the change file */
    else fatal("! Cannot open \"%s\" as change file", change_file_name);
		@.Cannot open change file@>
}

@ Here are some more variables relevant to the reading of input files.
Every time an input line is read coming from the change file, and on
returning to reading from |cur_file| (possibly after one or more lines of
the main input stream have been removed by the change file), we mark the
current section as having changed by setting a bit in the bitmap
|changed_section|.

@< Declarations... @>=
extern sixteen_bits section_count;
extern eight_bits changed_section[];
#define mark_section_as_changed(n) (changed_section[(n)>>3]|=1<<((n)&7))
#define section_changed(n) ((changed_section[(n)>>3]&(1<<((n)&7)))!=0)

@~
@<Defin...@>=
sixteen_bits section_count; /* the current section number */
eight_bits changed_section[(max_sections+7)/8]; /* is the section changed? */

@ The function |get_line| puts the next line of merged input into the
buffer and updates the variables appropriately. A space is placed at the
right end of the line (i.e., at |*limit|), serving many purposes, like
ensuring that a final `\.@@' on a line will be interpreted as \:\ ,
and at other times allowing us to test for interesting characters without
first testing for line end. The function returns |!input_has_ended| because
we often want to check the value of that variable after calling the
function. Usually |get_line| is called after the space at |*limit| has been
processed; therefore a call often takes the form of the macro |find_char|,
which calls |get_line| if necessary, and returns whether it has succeeded
in making |loc| point to a valid character (possibly a line-ending space).

The logic of marking sections as changed, which is implemented below, is a
bit subtle. The status of |changing| is noted for every line read, but since
the test is made at the start of |get_line|, it actually happens just before
the next line is read in. This means that if the first replacement line of a
change involves the start of a new section, then the new section is marked
as changed rather than the one before it, which is the right choice,
assuming that sections are started at the beginning of a line. It is
possible that |changing| is switched on and then right back off again (if
the line after \:y starts with \:z); in that case we mark the
current section as changed and restart |get_line|, since we still have found
no actual line.

@c
boolean get_line (void) /* inputs the next line */
{
restart:
  if (changing) mark_section_as_changed(section_count);
  else @<Read from |cur_file| and maybe turn on |changing|@>
  if (changing)
  { @<Read from |change_file| and maybe turn off |changing|@>
    if (!changing)
    {@; mark_section_as_changed(section_count); goto restart; }
  }
  loc=&buffer[0]; *limit= ' '; /* place sentinel space */
  if (compatibility_mode && buffer[0]=='@@' && buffer[1]=='i')
  { loc+=2; print_where=true;
    if (locate_file_name())
      push_input_file(false,changing);
    goto restart;
  }
  if (limit-buffer>5
   && strncmp(buffer,"#line",5)==0 && isspace((eight_bits)buffer[5]))
    @< Set file name and line number according to \&{\#line}
       directive and |goto restart| @>
  return !input_has_ended;
}

@ A \&{\#line} directive should have the form `\.{\#line 85 "common.w"}'.
The line number and file name simply override |cur_line| and |cur_file_name|.

@< Set file name and line number... @>=
{ sixteen_bits line=0;
  print_where=true; /* output a \&{\#line} directive soon */
  loc=&buffer[6]; @+ while (loc<limit && isspace((eight_bits)*loc)) ++loc;
  if (isdigit((eight_bits)*loc))
  { do line=10*line + *loc++ -'0'; while (isdigit((eight_bits)*loc));
    while (loc<limit && isspace((eight_bits)*loc)) ++loc;
    if (*loc++=='"')
    { int i=0; @+ while (&loc[i]<limit && loc[i]!='"') ++i;
      if (loc[i]=='"' && i<max_file_name_length)
      { struct f* cur_f= changing ? &change : &file[include_depth];
        cur_f->line=line-1; /* directive applies to next line, not this one */
        strncpy(cur_f->name,loc,i); cur_f->name[i]='\0';
	goto restart;
      }
    }
  }
  err_print("! Improper #line directive"); goto restart;
	    @.Improper \#line directive@>
}

@ After checking that a line could be obtained from the main input stream,
some quick tests are made that will avoid calling |check_change| in most
cases. The switch to |changing| mentioned in the module name is
usually brought about by |check_change|, but may also happen when the
suspension of the change file during inclusion of a header file is ended by
|get_web_line|. In either case further input lines should be taken from
the change file.

@< Read from |cur_file|... @>=
{ if (get_web_line()
   && change_limit>change_buffer
   && limit-buffer==change_limit-change_buffer
   && buffer[0]==change_buffer[0]
     ) check_change();
}

@ Here we get a line of input from the change file, unless it starts
with~\:z. The statements `|loc=buffer; *limit=' ';|' were performed in
|get_web_line| for lines from the main input stream, but must be issued
explicitly for lines from the change file.

@< Read from |change_file|... @>=
{ if (++change_line,!input_ln (change_file))
  { err_print("! Change file ended without @@z"); @.Change file ended...@>
  @/buffer[0]='@@'; buffer[1]='z'; limit=&buffer[2];
  }
  if (limit>&buffer[1] && buffer[0]=='@@') /* check if the change has ended */
    if (buffer[1]=='z')
    @/{@; prime_the_change_buffer(); changing=false; print_where=true; }
    else if (buffer[1]=='x' || buffer[1]=='y')
    @/{@; loc=&buffer[2]; err_print("! Where is the matching @@z?"); }
				     @.Where is the match...@>
    else @< Check for erron... @>
}

@ The function |check_complete| will be called at the end of \.{CTANGLE}
and \.{CWEAVE} to check for an unfinished state of the change file.

@< Declarations...@>=
extern void check_complete (void);

@~When |check_complete| is called we have |input_has_ended|, which implies
|!changing|. The only thing to test for is that there is no change line still
waiting for a match. In order to get a decent display of the non-matching
line, we copy it from the |change_buffer| to the |buffer|, and set some other
variable appropriately. There is no need to restore any variables, since
|check_complete| is called at the very end of the run.

@c
void check_complete (void) /* checks that all changes were picked up */
{ if (change_limit!=change_buffer)
  { int l=(int)(change_limit-change_buffer);
    strncpy(buffer,change_buffer,l); limit=&buffer[l];
    changing=true; loc=buffer; web_file_open=true;
      /* prepare unmatched line for display */
    err_print("! Change file entry did not match");
	       @.Change file entry did not match@>
  }
}


@* Storage of names and strings.
Both \.{CTANGLE} and \.{CWEAVE} store the strings representing identifiers,
module names and (in case of \.{CWEAVE}) index entries in a large array of
characters, called |byte_mem|. These strings are not accessed directly,
but via structures that collect further information about these objects.
These structures come in two kinds, depending on whether they correspond
to identifiers (or index entries) or to module names; these structures are
called |id_info| and |mod_info| respectively.

@< Public typedef...@>=

typedef struct id_info
{ char *byte_start; /* beginning of the name in |byte_mem| */
  @<More fields of |struct id_info|@>@;
} id_info, *id_pointer;
@)
typedef struct mod_info
{ char *byte_start; /* beginning of the name in |byte_mem| */
  @<More fields of |struct mod_info|@>@;
} mod_info, *mod_pointer;

@ All |id_info| and |mod_info| structures are stored in one of two arrays,
called |id_table| and |mod_table| respectively, and hence all |id_pointer|
and |mod_pointer| values point into these arrays. Therefore we can freely
convert between such a pointer and an index into the appropriate array;
the macros |id_index|, |id_at|, |mod_index| and |mod_at| defined above
perform these conversions.

@< Declarations... @>=
extern char byte_mem[], *byte_ptr;
extern id_info id_table[], *id_ptr;
extern mod_info mod_table[], *mod_ptr;

@~The first unused position in |byte_mem| is kept in |byte_ptr|, and the
first unused positions in |id_table| and |mod_table| are similarly kept in
|id_ptr| and |mod_ptr|, respectively. We want to keep
|byte_ptr<=byte_mem_end|, |id_ptr<=id_table_end| and
|mod_ptr<=mod_table_end|.

@d byte_mem_end  (&byte_mem[max_bytes]) /* end of |byte_mem| */
@d id_table_end  (&id_table[max_idents]) /* end of |id_table| */
@d mod_table_end (&mod_table[max_modules]) /* end of |mod_table| */

@<Definitions...@>=
char byte_mem[max_bytes]; /* characters of names */
char *byte_ptr=&byte_mem[0]; /* first unused position in |byte_mem| */
id_info id_table[max_idents]; /* information about identifiers */
id_pointer id_ptr=&id_table[0]; /* first unused position in |id_table| */
mod_info mod_table[max_modules]; /* information about module names */
mod_pointer mod_ptr=&mod_table[0]; /* first unused position in |mod_table| */

@ Here is a simple function that copies a (not necessarily null-terminated)
string~|s| of length~|l| into |byte_mem| and returns a pointer to the copied
string.

@c
char* store_string(char* s, int l)
{ char* dest=byte_ptr;
  if (byte_mem_end-byte_ptr<=l) overflow ("byte memory");
  byte_ptr+=l; *byte_ptr++='\0'; return strncpy(dest,s,l);
}

@ A component is present in both |id_info| and |mod_info| structures,
whose function is different for \.{CTANGLE} and \.{CWEAVE}. In |CTANGLE|,
it is used only in |mod_info| structures, and is a pointer to a replacement
text for the module name. In |CWEAVE| it is a pointer to a list of
cross-references for the identifier or module name. The precise nature of
these pointers is of no interest to the common code, and at this point we do
not have the types available to express this field as a |union|. However,
since in both cases it will be a pointer to a structure, we can use a little
trick by declaring it as a pointer to a |struct variant| where |variant| is
not further specified here. Then in \.{CTANGLE} and \.{CWEAVE} we define
|variant| as a macro for the appropriate |struct| specifier, while in the
common code we leave it undefined (the \:d defining |variant| is given
very early, so that it will precede the `\.{@@h "common.h"}' that will read
in the definition of |struct id_info| and |struct mod_info|).

It is not quite proper that in the common code these structures contain a
pointer to a named but undefined structure, while in the other compilation
unit of \.{CTANGLE} and \.{CWEAVE} respectively, they contain a pointer to
a known structure with a {\it different\/} specifier (due to the macro
definition); in \Cee\ structures with different specifiers can never be the
same type. However, since the linker does not type-check, and the compiled
code is not likely to depend on the name of a structure specifier, this
should cause no problems, except possibly some mild confusion to a debugger.
An alternative solution of defining a structure of the same name to denote
entirely different types in \.{CTANGLE} and \.{CWEAVE} would be clearer to
the computer, but more confusing to humans. Of course we could also have
used a |(void*)| field instead of a pointer to a structure, but this would
require a lot of explicit or implicit casts, allowing possible type errors
to go undetected by the compiler.

@< More fields of |struct id_info| @>=
struct variant* equiv_or_xref; /* extra information about this identifier */
@~
@< More fields of |struct mod_info| @>=
struct variant* equiv_or_xref; /* extra information about this module */

@ The |id_info| structures are linked together into (hopefully) short lists
sharing the same hash key by means of the |hash_link| field. In |CWEAVE|,
an additional field is used to store information about the syntactic r\^ole
of the identifier, its so-called |ilk|.

@< More fields of |struct id_info| @>=
struct id_info *hash_link; /* links identifiers with same hash code */
int ilk; /* syntactic information, used in |CWEAVE| only */

@ The pointers to the heads of the hash lists are stored in the
array~|hash|. Identifiers can be found by computing a hash code~|h| and then
looking at the identifiers represented by the |id_pointer|s |hash[h]|,
|hash[h]->hash_link|, |hash[h]->hash_link->hash_link|,~\dots, until either
finding the desired name or encountering the null pointer. Thus the hash
table consists of entries of type |id_pointer|; it is maintained by the
function |id_lookup|, which finds a given identifier and returns the
appropriate |id_pointer|. The basic matching is done by the function
|names_match|, which is slightly different in \.{CTANGLE} and \.{CWEAVE}. If
there is no match for the identifier, it is inserted into the table.

@< Declarations... @>=
extern id_pointer hash[];
#define hash_end  (&hash[hash_size]) /* end of |hash| */
id_pointer id_lookup(char*,char*,int);

@~The |hash| table has size~|hash_size|.

@<Definitions...@>=
id_pointer hash[hash_size]; /* heads of hash lists */

@ Initially all the hash lists are empty. Although strictly speaking the
array |hash| should be initialised to null pointers by the compiler,
we don't count too much on this in case such pointers are not represented
by null bit patterns.

@<Init...@>=
{@; int i=hash_size; do hash[--i]=NULL; while(i>0); }

@ Here is the main function for finding identifiers (and index entries).
The parameters |first| and |last| point to the beginning and end of the
string. The parameter |ilk| is used by |CWEAVE| only in the course of the
matching process; within |id_lookup| it is just passed on to the functions
|names_match| and~|init_id_name|, which are defined in \.{CTANGLE} and
\.{CWEAVE} separately. We facilitate the initialisation of the reserved
words in |CWEAVE| by allowing the use of null-terminated strings rather
than using a pointer to the end of the string: in this case |first| should
point to such a string, and |last==NULL|.

@c
id_pointer id_lookup (char* first,char* last,int ilk)
  /* look up an identifier */
{ int l,h; /* length and hash code of the given identifier */
  if (last==NULL) last=first+(l=(int)strlen(first)); /* null-terminated string */
  else l=(int)(last-first); /* compute the length */
  @<Compute the hash code |h| of the string from |first| to |last| @>
  @<Compute and |return| the name location,
    possibly by entering a new identifier into the tables@>
}

@ A simple hash code is used: If the sequence of character codes is
$c_1c_2\ldots c_n$, its hash value will be
$$ (2^{n-1}c_1+2^{n-2}c_2+\cdots+c_n)\bmod\hbox{|hash_size|}. $$
The purist might worry about the empty string, which, although of course
impossible as identifier, might be flagged as an index entry; such an index
entry is however explicitly ruled out in |CWEAVE|.

@<Compute the hash...@>=
{ char* p=first;
  h=*p; while (++p<last) h=((h<<1)+*p)%hash_size;
}

@ Here we either find an existing node for the identifier, or create one,
inserting it into the name table and the hash list. The maximal average
length |max_idents/hash_prime| of hash lists should be sufficiently small
that linearly searching them is quite fast.

@<Compute and |return|...@>=
{ id_pointer p=hash[h]; /* the head of the hash list */
  while (p!=NULL && !names_match(p,first,l,ilk)) p=p->hash_link;
  if (p==NULL) /* we haven't seen this identifier before */
  @< Make |p| point to a fresh |id_node| that refers to a string copied from
     |first|, and prepend the node to hash list |h| @>
  return p;
}

@ The information associated with a new identifier must be initialised
in |CWEAVE| in a way that does not apply to~\.{CTANGLE}; hence the
function |init_id_name|.

@< Make |p| point to a fresh |id_node|... @>=
{ p=id_ptr; /* this is where the new name entry will be created */
  if (id_ptr++>=id_table_end) overflow ("identifier");
  name_begin(p)=store_string(first,l);
  if (program==cweave) init_id_name(p,ilk);
  p->hash_link=hash[h]; hash[h]=p; /* insert |p| at beginning of hash list */
}

@ The names of modules are stored in |byte_mem| together with the identifier
names, but a hash table is not used for them because \.{CWEB} needs to be
able to recognise a module name when given a prefix of that name. To this
end the |mod_info| structures are linked into a binary search tree. The only
unconventional thing about this search tree is that the search key is not
necessarily equal to the full string stored but that a prefix is decisive
for matching purposes: once that prefix matches, not matching the remainder
of the string is considered to be an error. The reason for this is that the
module name may have been specified first by a short prefix and is later
extended; we need to be aware of the ambiguity when a second incompatible
extension of the original prefix is attempted. Therefore a field
|key_length| is included in |struct mod_info| that indicates the length of
the prefix that should determine a unique match; its value is the length of
the shortest partial specification of this module name that has been
encountered so far.

@< More fields of |struct mod_info| @>=
struct mod_info *llink,*rlink;
  /* left and right links in binary search tree */
int key_length; /* number of characters decisive for match */

@ There is one more attribute that must be recorded for a module name,
namely whether the string accessed by the |byte_start| field is the full
module name, or just the longest prefix specified so far (i.e., all
occurrences so far ended with `\.{...}'). Rather than extending |struct
mod_info| with another (boolean) field or encoding this information in one
of the other fields we use a small trick to record this information. All
full module names will be stored in |byte_mem|, and since that is filled
with null-terminated strings, we will have for any |mod_pointer p@;|
representing a complete module name that |p->byte_start[-1]=='\0'|. The
macro |complete_name| defined above tests this condition, and we shall
store incomplete module names in such a way that the condition fails for
them. This method costs just one byte of temporary storage for each module
name that is incompletely specified on its first occurrence, except that a
one-time sacrifice of another byte is needed to guarantee that even the
first name is preceded by a null byte.

@< Initialise... @>=
*byte_ptr++='\0'; /* prefix a null byte to the first string */

@~Like most trees, the search tree for module names starts at its |root|.

@< Declarations... @>=
extern mod_pointer root;

@~The search tree starts out empty.
@< Definitions... @>=
mod_pointer root=NULL;
 /* the root of the binary search tree for module names */

@ The set of all possible strings can be conceptually arranged in an
infinite tree with the empty string at its root, where the direct descendents
of a string are those strings obtained by appending a single letter to it,
listed in alphabetical order. The position of one string relative to another
in this tree is one of five possibilities: to the left of, to the right of,
ancestor of, descendent of, or equal to. The type |enum mod_comparison|
encodes these five possibilities.

@< Prototypes of local functions @>=
enum mod_comparison @/
{ less, /* the first name is lexicographically less than,
	   but no prefix of the second */
  equal, /* the first name is equal to the second */
  greater, /* the first name is lexicographically greater than,
	      but no extension of the second */
  prefix, /* the first name is a proper prefix of the second */
  extension /* the first name is a proper extension of the second */
};

@~The function |mod_name_cmp| will determine which one of the relations
described above holds for a given pair of strings, each represented by a
pointer to its beginning and its length.

@c
local enum mod_comparison mod_name_cmp
    (char* p, int l1, char* q, int l2)
{ int l= l1<l2 ? l1 : l2;
  while (--l>=0) if (*p++!=*q++) return *--p<*--q ? less : greater;
  return l1<l2 ? prefix : l1>l2 ? extension : equal;
}

@ When a new module name is found and approved, we first store the string
representing it in an appropriate place, and then call |new_mod_node| to
install a node for the name in the |mod_table| array and prepare it for
inclusion in the search tree. The information associated with the new node
must be initialised in a slightly different way in |CWEAVE| than in
|CTANGLE|; hence the call of |init_module_name|.

@c
local mod_pointer make_mod_node (char* name)
{ mod_pointer node=mod_ptr; /* allocate new node */
  if (mod_ptr++>=mod_table_end) overflow ("module name");
  name_begin(node)=name;
  node->llink=NULL; node->rlink=NULL;
  init_module_name(node); /* initialise new node */
  return node;
}

@ The function |mod_name_lookup| is used to look up a complete module name
in the search tree, inserting it if necessary, and returns a pointer to
where it was found. Its parameters give the beginning and length of the
module name. The name might be illegal, for instance if it is a prefix of an
existing name, in which case |NULL| is returned rather than a pointer into
|mod_table|, after printing an error message.
@c
local mod_pointer mod_name_lookup (char* name, int l)
  /* finds complete module name */
{ mod_pointer p; /* current node of the search tree */
  mod_pointer* loc=&root; /* |p| will come from this location */
  while ((p=*loc)!=NULL)
  { int l0=p->key_length; char* key=name_begin(p);
    switch (mod_name_cmp(name,l,key,l0))
    { case less: loc=&p->llink; break;
      case greater: loc=&p->rlink; break;
      case equal: case extension:
	@< Check that |name| matches the remainder of |key|;
	   if so, complete |p| if necessary and |return p|,
	   if |name| is a prefix of~|key| fall through,
	   and otherwise report an error and |return NULL| @>
      case prefix:
	err_print("! Incompatible module name"); @.Incompatible module name@>
	print("\nName is a prefix of <%s%s>.\n" @.Name is a prefix of <...>@>
	     ,key, complete_name(p) ? "" : "...");
      return NULL; /* dummy module name */
    }
  }
  @< Copy |name| into |byte_mem|, install a new node in the tree at |*loc|,
     and |return| it @>
}

@ When the name is not found in the tree, a new module has been specified
whose name is complete from the start, and can therefore be installed in
|byte_mem|. Its |key_length| field is set to the length~|l| of the name as
specified here, but it might be decreased later.

@< Copy |name| into |byte_mem|... @>=
{ (p=make_mod_node(store_string(name,l)))->key_length=l; /* prepare new node */
  return *loc=p; /* install new node into tree */
}

@ When the module name |name| being looked up matches the string |key| in a
node~|p| of the search tree for the first |p->key_length| characters, then
this is guaranteed to be the only such match in the tree. It is required
that the rest of |key| matches the remainder of~|name| as well, and if
|complete_name(p)| holds then the match must exhaust~|name|.

@< Check that |name| matches the remainder of |key|... @>=
{ enum mod_comparison cmp=
    mod_name_cmp(name+l0,l-l0,key+l0,(int)strlen(key+l0));
  switch(cmp)
  { case less: case greater:
      err_print("! Incompatible module name"); @.Incompatible module name@>
      print("\nName inconsistently extends <%.*s...>.\n",l0,key);
	     @.Name inconsistently extends <...>@>
      return NULL;
    case extension: case equal:
      if (complete_name(p))
	if (cmp==equal) return p;
	else
	{ err_print("! Incompatible module name"); @.Incompatible module name@>
	  print("\nPrefix exists: <%s>.\n",key); return NULL;
		 @.Prefix exists@>
	}
      name_begin(p)=store_string(name,l);
        /* install |name| in place of |key| */
      @< Deallocate memory for |key| @>
      return p;
  }
}

@ The function |prefix_lookup| is similar to |mod_name_lookup|, but it is
called when the specification of a module name ends with `\.{...}'. Here,
unlike in |mod_name_lookup|, there is the possibility that there is more
than one match. The function decides whether there are $0$,~$1$, or~$>1$
matches present in the tree, but in the final case it need not find all
matches.  There is a match of |name| with node~|p| if the first
|p->key_length| characters at |name_begin(p)| are equal, an extension or a
prefix of~|name| (in the last case we should not have |completed(p)|). It is
clear that if there is a match with two distinct nodes, then there is also a
match with their closest common ancestor in the search tree, and since we do
not allow inclusions among the significant parts of the keys in the tree,
this situation can only occur if |name| is a proper prefix for all its
matches.  Hence, once a first match is found, any further match must be
among its descendents in the tree.

The value of |loc| is maintained as in |mod_name_lookup| to allow insertion
of a new node if no match is found; however, once we have |match!=NULL|, the
value of~|loc| becomes irrelevant.

@c
local mod_pointer prefix_lookup (char* name,int l)
  /* finds module name given a prefix */
{ mod_pointer p=root,* loc=&root; /* current node and where it comes from */
  mod_pointer match=NULL; /* the first matching node, if any */
  mod_pointer saved=NULL; /* another subtree that might have matches */
  while (p!=NULL)
  { int l0=p->key_length; char* key=name_begin(p);
    switch (mod_name_cmp(name,l,key,l0))
    { case less: p=*(loc=&p->llink); break;
      case greater: p=*(loc=&p->rlink); break;
      case equal: return p; /* a match, and no other matches are possible */
      case extension:
	@< Check that |name| matches the remainder of~|key|;
	   if so, extend its string if necessary and |return p|,
	   otherwise report an error and |return NULL| @>
      case prefix:
	if (match!=NULL)
	{@; err_print("! Ambiguous prefix"); return NULL; }
		       @.Ambiguous prefix@>
	match=p; saved=p->rlink; p=p->llink; /* |loc| is irrelevant now */
   }
    if (p==NULL && match!=NULL)
      p=saved, saved=NULL; /* search other subtree */
  }
  if (match==NULL)
    @< Copy the incomplete |name| in a temporary place, install a new node
       in the tree at |*loc|, and |return| it @>
  match->key_length=l; /* |name| is a shorter prefix than used before */
  return match;
}

@ Incomplete nodes are not stored in |byte_mem|, but temporarily set aside
in dynamic memory, which will be freed when the name is completed. A
non-zero byte is installed at the beginning to make the |complete_name|
macro function properly. Like for complete names, the initial |key_length|
is the full length of the string specified.

@< Copy the incomplete |name|... @>=
{ char* key=(char*)malloc(l+2);
  if (key==NULL) fatal("Out of dynamic memory!");
  *key++='\1'; /* ensure that |complete_name(p)| is false afterwards */
  strncpy(key,name,l); key[l]='\0'; /* store the incomplete name */
  (p=make_mod_node(key))->key_length=l; /* prepare new node */
  return *loc=p; /* install new node into tree */
}

@~When freeing the dynamic memory, we mustn't forget to subtract~1.

@< Deallocate memory for |key| @>=  @+ free(key-1);

@ When we come here, we have found that |name| extends the first
|p->key_length| characters of the name stored in~|p|, and therefore cannot
similarly match any other nodes.

@< Check that |name| matches the remainder of~|key|;
   if so, extend its string if necessary and |return p|,
   otherwise report an error and |return NULL| @>=
{ enum mod_comparison cmp=
    mod_name_cmp(name+l0,l-l0,key+l0,(int)strlen(key+l0));
  switch(cmp)
  { case less: case greater:
      err_print("! Incompatible module name"); @.Incompatible module name@>
      print("\nName inconsistently extends <%.*s...>.\n",l0,key);
	     @.Name inconsistently extends <...>@>
      return NULL;
    case prefix: case equal: return p;
    case extension:
      if (complete_name(p))
      { err_print("! Incompatible module name"); @.Incompatible module name@>
	print("\nPrefix exists: <%s>.\n",key); return NULL; @.Prefix exists@>
      }
      @< Replace name stored in |p| by the larger prefix |name| @>
      return p;
  }
}

@ We come here in the rather unusual case that a second larger prefix of the
same module name is given before the full name is specified. We discard the
old memory for |key| and replace it by a fresh copy of |name|, which is
easier and probably more efficient than to try to reuse the old |key| by
using |realloc|. As for the initial dynamic allocation, we must not forget
the byte at~|key[-1]| here.

@< Replace name stored in |p| by the larger prefix |name| @>=
{ @< Deallocate memory for |key| @>
  if ((key=(char*)malloc(l+2))==NULL) fatal("Out of dynamic memory!");
  *key++='\1'; /* ensure that |complete_name(p)| is false afterwards */
  strncpy(key,name,l); key[l]='\0'; /* store the incomplete name */
  name_begin(p)=key; /* install new name in node |p| */
}

@* Lexical routines.
Much of the code of \.{CWEB} deals with lexical matters such as recognising
tokens, control texts, etc. This inevitably involves some complications
since the relevant set of lexical rules varies from place to place in the
source file, being sometimes that of \TeX, sometimes of \Cee, and sometimes
rules defined by \.{CWEB} itself.  Furthermore the source file is read in on
three different occasions: once by \.{CTANGLE} and twice by \.{CWEAVE}; the
level of detail in which the source is inspected differs between these
passes. We can alleviate the task of lexical recognition by collecting in
the common code a few functions that will scan specific lexical items; we do
this for various classes of ``large'' lexical items: module names, control
texts and strings.

In all cases we copy the characters of the lexical item to a separate place,
making minor modifications as necessary (for instance replacing \:@@
by~`\.@@'). To this end there is a buffer |mod_text| used to copy the
characters into; its name is derived from its most important purpose of
storing module names, for which a separate buffer is definitely needed since
these names can be spread over several lines of input. There are also two
pointers |id_first| and~|id_loc| that are used to point to the beginning
and end of lexical entities.

@< Declarations... @>=
extern char mod_text[], *id_first, *id_loc;
#define mod_text_end (&mod_text[longest_name+1]) /* end of |mod_text| */
mod_pointer get_module_name (void);
boolean get_control_text(void);
void get_string(void);

@~The character |mod_text[0]| will not be used to store any actual
characters, whence we make |mod_text| one slot longer than |longest_name|.

@<Definitions...@>=
char mod_text[longest_name+1]; /* name being sought for */
char *id_first; /* where the current identifier begins in the buffer */
char *id_loc; /* just after the current identifier in the buffer */

@ The function for scanning and looking up module names is
|get_module_name|; it will read the module name from the input, and decide
how to search depending on the presence or absence of an ellipsis. This
function is called when \:< has just been scanned; it will scan the
module name and look it up by |mod_name_lookup| or |prefix_lookup| as
appropriate; it returns a |mod_pointer| to the name found, or |NULL| if the
name was found to be erroneous. Consequent to the limit, \.{CWEB} will treat
`\.{@@<...@@>}' as a valid module name if and only if if there exactly one
module name in the entire program.

@c
mod_pointer get_module_name (void)
{ @< Put module name into |mod_text|,
     set |id_first| and |id_loc| appropriately @>
  { int l=(int)(id_loc-id_first);
    return l>=3 && strncmp(id_loc-3,"...",3)==0
    @/	? prefix_lookup(id_first,l-3) : mod_name_lookup(id_first,l);
  }
}

@ Module names are placed into the |mod_text| array with consecutive spaces,
tabs, and carriage-returns replaced by single spaces. There will be no
spaces at the beginning or the end. We set |mod_text[0]=' '| to facilitate
this, and put |id_first=&mod_text[1]| for the calls to the lexical routines
that use |mod_text|.

@<Init...@>=
mod_text[0]=' ';

@~A module name of exactly |longest_name| characters will cause an error
message, even though all its characters are stored; this is because we
don't take the effort to distinguish a `full' state of the buffer from an
`overflowed' one.

@< Put module name... @>=
{ eight_bits c; char* k=mod_text; /* points to last recorded character */
  do
  { if (!find_char())
    {@; err_print("! Input ended in module name"); break; }
		   @.Input ended in module name@>
    c=*loc++;
    @<Handle |c=='@@'|; if end of name, |break|@>
    if (isspace(c)) c=' '; /* convert tabs, newlines etc. */
    if (k<mod_text_end-1 && !(c==' ' && *k==' ') ) *++k=c;
  } while(true);
  id_first=&mod_text[1];
  if (k>=mod_text_end-1)
@/{@; print("\n! Module name too long: "); @.Module name too long@>
    term_write(id_first,25); err_print("..");
  }
  id_loc= *k==' ' && k>mod_text ? k : k+1;
    /* point after last non-space character */
}

@ Except for the beginning of new sections we are not fussy about unexpected
control sequences within module names, since we might be within `\pb';
any problems will be eventually be detected during the output of module
names by |CWEAVE|. So we are only looking for \:> (incidentally this
means that control texts like those introduced by \:t are forbidden
in module names).

@< Handle |c=='@@'|... @>=
if (c=='@@')
{ if ((c=*loc++)=='>') break;
  if (isspace(c) || c=='*' || c=='~')
  @/{@; err_print("! Module name didn't end"); loc-=2; break; }
		   @.Module name didn't end@>
  if (k<mod_text_end-1) *++k='@@';
    /* record the `\.{@@}', now |c==loc[-1]| again */
}

@ The function |get_control_text| handles control texts. Like
|get_module_name| it copies text until \:> into the |mod_text| array,
and sets |id_first| and |id_loc| appropriately. The function returns a
boolean value telling whether the control text was in fact empty, since that
case sometimes needs special treatment. Always using |get_control_text|
guarantees uniformity in the recognition of various types of control texts
(not counting module names).

@c
boolean get_control_text(void)
{ char c,* k=id_first=&mod_text[1]; /* points after last recorded character */
  do
    if ((*k++=*loc++)=='@@')
      if ((c=*loc++)!='@@')
      { if (c!='>')
	  err_print("! Control codes are forbidden in control text");
		     @.Control codes are forbidden...@>
	return (id_loc=k-1)==id_first;
      }
  while(loc<=limit);
  err_print("! Control text didn't end"); @.Control text didn't end@>
  return (id_loc=k)==id_first;
}

@ And here is a similar function |get_string| that is used to scan strings
and character constants. These can contain newlines or instances of their own
delimiters if they are protected by a backslash. We follow this convention,
but do not allow the string to be longer than |longest_name|. The macro
|copy_char| copies a character to |mod_text| without risking overflow; the
|else| branch of its definition serves only to ensure that any side effect
of its argument is always performed.

@d copy_char(c) @+ if (id_loc<mod_text_end) *id_loc++=c; @+ else (void)(c) @;

@c
void get_string(void)
{ char c, delim = loc[-1]; /* what started the string */
  id_loc=id_first = &mod_text[1]; copy_char(delim);
  if (delim=='L')
    *id_loc++=delim=*loc++; /* after `\.L' comes the real delimiter */
  else if (delim=='<') delim='>'; /* for file names in \&{\#include} lines */
  do
  { if (loc>=limit)
    {@; err_print("! String didn't end"); loc=limit; break; }
		   @.String didn't end@>
    copy_char(c=*loc++);
    if (c=='\\')
      if (loc<limit) copy_char(*loc++);
        /* record `\.\\|c|' with |c!='\n'| literally */
      else if (get_line())
        if (program==cweave) --id_loc; /* |CWEAVE| erases backslash */
	else copy_char('\n'); /* but |CTANGLE| copies escaped newline */
      else
      {@; loc=buffer;
        err_print("! Input ended in middle of string");
		   @.Input ended in middle of string@>
        break;
      }
    else if (!including_header_file && c=='@@')
      if (*loc=='@@') ++loc; /* undouble \:@@ */
      else err_print("! Double @@ required in strings");
		      @.Double @@ required...@>
  }
  while (c!=delim);
  if (id_loc>=mod_text_end)
  @/{@; print("\n! String too long: "); @.String too long@>
    term_write(mod_text+1,25); err_print("..");
  }
}


@* Reporting errors to the user.
A global variable called |history| will contain one of four values at the
end of every run: |spotless| means that no unusual messages were printed;
|harmless_message| means that a message of possible interest was printed
but no serious errors were detected; |error_message| means that at least
one error was found; |fatal_message| means that the program terminated
abnormally. The value of |history| does not influence the behaviour of the
program; it is simply computed for the convenience of systems that might
want to use such information.

@< Declarations... @>=
extern history; /* indicates how bad this run was */
extern void err_print (char*), wrap_up (void), print_stats (void),
       fatal (char*,...);

@~Despite the name |history|, each run starts in an unblemished state.
@< Definitions... @>=
int history=spotless; /* indicates how bad this run was */

@ The call |err_print("! Error message")| will report an error to the user,
printing the error message at the beginning of a new line and then giving an
indication of where the error was spotted in the source file. Note that the
string passed to |err_print| does not end with a period, since one will be
automatically supplied, as will an initial newline if the string begins
with~|"!"|.

@c
void err_print (char *s) /* prints `\..' and location of error message */
{ print(*s=='!' ? "\n%s." : "%s.",s);
  if (web_file_open) @<Print error location based on input buffer@>
  update_terminal(); mark_error();
}

@ The error locations can be indicated by using the global variables
|loc|, |cur_line|, |cur_file_name| and |changing|, which tell
respectively the first unlooked-at position in |buffer|, the current
line number, the current file, and whether the current line is from
|change_file| or |cur_file|.

@< Print error location based on input buffer @>=
{ char *k, *l=(loc<limit) ? loc : limit; /* pointers into |buffer| */
  if (changing) printf(" (l. %d of change file)\n", change_line);
  else if (include_depth==0) printf(" (l. %d)\n", cur_line);
  else printf(" (l. %d of include file %s)\n", cur_line, cur_file_name);
  if (l>buffer)
  { for (k=buffer; k<l; k++) putchar(*k=='\t' ? ' ': *k);
    new_line();
    for (k=buffer; k<l; k++) putchar(' '); /* space out the next line */
  }
  for (k=l; k<limit; k++) putchar(*k); /* print the part not yet read */
}

@ When a fatal error is encountered, or when the program comes to its normal
termination, the function |wrap_up| is called, which will not return to the
caller. Some implementations may wish to pass the |history| value to the
operating system so that it can be used to govern whether or not other
programs are started. Here, for instance, we pass the operating system
a status of~0 if and only if at worst only harmless messages were printed.
@^system dependencies@>

@c
void wrap_up (void)
{
#ifdef STAT
  if (show_stats) print_stats(); /* print statistics about memory usage */
#endif
  @< Print the job |history| @>
  exit(history>harmless_message);
}

@ A spotless status is reported if |show_happiness| is true, any other
exit status is always reported.

@< Print the job |history| @>=
{ static char* mess[]=
  { "No errors were found.",
    "Did you see the warning message above?",
    "Pardon me, but I think I spotted something wrong",
    "That was a fatal error, my friend."
  };
  if (show_happiness || history>0) print("\n(%s)\n",mess[history]);
}

@ When an error or overflow condition is detected for which no recovery is
possible, we call |fatal|, which aborts the program after pointing out the
source of trouble and wrapping up as graciously as possible.

@h <stdarg.h> /* needed to define variadic functions */

@c
void fatal(char* s,...)
{ va_list p; va_start(p,s);
  vprintf(s,p); va_end(p); err_print("");
    /* print reason and location of fatal stop */
  history=fatal_message; wrap_up();
}


@* Command line arguments.
The user calls \.{CWEAVE} and \.{CTANGLE} with one or more arguments on the
command line. These are either file names or sets of flags to be turned on
(beginning with |"+"|) or off (beginning with |"-"|); in case the special
flag |'i'| occurs, the remainder of that argument is interpreted as a string
rather than as a set of flags.

The following globals are for communicating the user's desires to the rest
of the program. The various file name variables contain strings with the
names of those files. Most of the flags are undefined but available for
future extensions.

@< Declarations... @>=
extern boolean flags[];
extern char C_file_name[],idx_file_name[],scn_file_name[];

@~@< Definitions... @>=
boolean flags[UCHAR_MAX+1]; /* an option for each character code */
char C_file_name[max_file_name_length]; /* name of |C_file| */
local char tex_file_name[max_file_name_length]; /* name of |tex_file| */
char idx_file_name[max_file_name_length]; /* name of index file */
char scn_file_name[max_file_name_length]; /* name of module names file */
local boolean change_file_explicit=false;
  /* was a change file argument specified? */

@ We now must look at the command line arguments and set the file names
accordingly.

@< Prototypes of local functions @>=
local void scan_args (int argc,char** argv);

@~At least one file name must be present: the \.{CWEB} file.  It may have an
extension, or it may omit the extension to get |".w"| or |".web"| added.
The \TeX\ output file name is formed by replacing the \.{CWEB} file name
extension by |".tex"|, and the \Cee\ file name by replacing the extension
by~|".c"|.

If a second file name is given, it is the change file, again either with an
extension or without one to get |".ch"|. An omitted change file argument
means that the change file name is that of the \.{CWEB} file with the
extension replaced by |".ch"| if that file exists, or none at all otherwise.
If present, a third file name replaces the default output file name,
possibly including the extension.

@c
local void scan_args (int argc,char** argv)
{ char *dot_pos; /* position of rightmost |'.'| in the argument */
  int files_found=0, paths_found= at_h_path[0].name==NULL ? 0 : 1;
  while (--argc>0) /* first ``argument'' (program name) is irrelevant */
    if (((*++argv)[0]=='+' || (*argv)[0]=='-') && (*argv)[1]!='\0')
      @<Handle flag argument@>
    else
    { if (strlen(*argv)+5>max_file_name_length)
	/* we need room to add things like `\.{.web}' */
	fatal("! Filename too long:\n%s", *argv);
      dot_pos=strrchr(*argv,'.');
      switch (++files_found)
      { case 1:
	@< Make |web_file_name|, and defaults for other file names,
	   from~|*argv| @> @+
      break; case 2: @< Make |change_file_name| from |*argv| @> @+
      break; case 3: @< Make output file names from |*argv| @> @+
      break; default: @< Print usage error message and quit @>
      }
    }
  if (files_found==0) @< Print usage error message and quit @>
  if (paths_found<max_include_paths)
    at_h_path[paths_found].name=NULL; /* mark end of list */
}

@ We use |*argv| for the |web_file_name| if there is a |'.'| in it,
otherwise we add |".w"|. We prepare an alternative name |alt_web_file_name|
by adding |"web"| after the dot. The other file names come from adding other
things after the dot. Since there is no general agreement about the proper
extension for \Cpp\ file names, we use a macro |CPPEXT| for it, which
defaults to |"C"| but can be overridden by predefining this preprocessor
symbol when \.{common.c} is compiled.

@< Make |web_file_name|... @>=
#ifndef CPPEXT
#define CPPEXT "cpp"
  /* extension for \Cpp\ file names; should not exceed 3 characters */
#endif
{ if (dot_pos==NULL) sprintf(web_file_name,"%s.w",*argv);
  else
  { sprintf(web_file_name,"%s",*argv); /* use file name and extension */
    *dot_pos='\0'; /* truncate the name before the dot */
  }
  sprintf(alt_web_file_name,"%s.web",*argv);
  sprintf(change_file_name,"%s.ch",*argv);
  if (program==ctangle)
    sprintf(C_file_name,"%s.%s",*argv, C_plus_plus ? CPPEXT : "c");
  else
  { sprintf(tex_file_name,"%s.tex",*argv);
    sprintf(idx_file_name,"%s.idx",*argv);
    sprintf(scn_file_name,"%s.scn",*argv);
  }
}

@ If in place of the change file name a `\.-' is specified, we clear the
|change_file_name|, and if a `\.+' is specified we retain the default
(but |files_found| is increased so that the next file name will be the
output file); otherwise, |change_file_explicit| is raised and the default
change file name is replaced by the given one.

@< Make |change_file_name|... @>=
if ((*argv)[0]=='-') change_file_name[0]='\0';
else if ((*argv)[0]!='+')
@/{@; change_file_explicit=true;
  sprintf(change_file_name,dot_pos==NULL ? "%s.ch" : "%s", *argv);
}

@ If an output file name is given when calling \.{CWEAVE}, its base name is
also used for |idx_file_name| and |scn_file_name|.

@< Make output file names... @>=
if (program==ctangle)
  if (dot_pos!=NULL) sprintf(C_file_name, "%s", *argv);
  else sprintf(C_file_name,"%s.%s", *argv, C_plus_plus ? CPPEXT : "c");
else
{ if (dot_pos!=NULL)
    {@; sprintf(tex_file_name, "%s", *argv); *dot_pos='\0'; }
  else sprintf(tex_file_name,"%s.tex", *argv);
  sprintf(idx_file_name,"%s.idx",*argv);
  sprintf(scn_file_name,"%s.scn",*argv);
}

@ The |flags| are turned off initially; any flags that are on by default are
set before calling |common_init|.

@< Set the default options common to \.{CTANGLE} and \.{CWEAVE} @>=
show_banner=show_happiness=show_progress=true;

@~Flags are made case-insensitive, since some operating systems do not
allow passing of both lower and upper case arguments. The |'i'| flag gets
a special treatment.

@< Handle flag... @>=
{ boolean flag_change=(**argv == '+');
  char* p=&(*argv)[1]; unsigned char c;
  while ((c=*p++)!='\0')
    if ((c=tolower(c))!='i') flags[c]=flag_change;
    else @< Store string |p| as new include path and |break| @>
}

@ We copy include paths into |byte_mem|; since we are at the very beginning
of the run, there is no chance that we overflow |byte_mem| already, so that
we need not include a fourth error message below.

@< Store string |p| as new include path and |break| @>=
{ size_t l=strlen(p);
  if (l==0) err_print("! Empty include path");
		       @.Empty include path@>
  else if (l>max_path_length) err_print("! Include path too long");
  					 @.Include path too long@>
  else if (paths_found>=max_include_paths)
    err_print("! Too many include paths");
	       @.Too many include paths@>
  else
  { at_h_path[paths_found].length=(int)l;
    at_h_path[paths_found++].name=strcpy(byte_ptr,p);
    byte_ptr+=l+1;
  }
  break;
}

@ The usage message gives only a sketch of the proper way to call the
\.{CWEB}. programs; more details can be found in the manual.

@< Print usage error message and quit @>=
fatal("! Usage:\n" @+
"c%se [(+|-)options] cwebfile[.w] [(changefile[.ch]|+|-) [outputfile[.%s]]]"
      , program==ctangle ? "tangl" : "weav"
      , program==ctangle ? "c" : "tex");

@* Output. The only thing done in the common code for the output files is
declaring the variables through which they are accessed, and opening them.

@< Declarations... @>=
extern FILE *C_file, *tex_file;
void open_output_file(void);

@~Having separate variables for the output files of \.{CTANGLE} and
\.{CWEAVE} only serves to give them more meaningful names.

@<Definitions...@>=
FILE *C_file; /* where output of \.{CTANGLE} goes */
FILE *tex_file; /* where output of \.{CWEAVE} goes */

@~@c
void open_output_file(void)
{ char* name; FILE** file;
  if (program==ctangle) {@; name=C_file_name; file=&C_file; }
  else  {@; name=tex_file_name; file=&tex_file; }
  if ((*file=fopen(name,"w"))==NULL)
    fatal("! Cannot open \"%s\" as output file",name);
	   @.Cannot open output file@>
}

@ All regular terminal output passes through the function |print|, which
takes the place of |printf|; it suppresses initial newlines if the line is
already empty. The function |print_progress| and |print_section_progress|
are used to display progress reports on the screen.

@< Declarations...@>=
void print(char*,...), print_progress(char*), print_section_progress(void);

@~A variable |term_line_empty| keeps track of whether we are at the
beginning of an output line, so that we can force certain terminal output
to start on a new line without producing empty lines.

@<Definitions...@>=
local boolean term_line_empty=true;
	/* has anything been written to the current line? */

@~Any initial or final newline produced by~|print| is caused by the format
string, not by any of the arguments following it; this makes it easy to
maintain |term_line_empty|.

If a function produces terminal output otherwise than through one of the
three functions below, it should take care of |term_line_empty| itself;
for instance |err_print| produces the error context using |printf| and
|putchar|, knowing that |term_line_empty==false| at that time, which is
also the correct value afterwards.

@c
void print(char* s,...)
{ va_list p; va_start(p,s);
  if (term_line_empty && *s=='\n') ++s; /* avoid printing empty line */
  vprintf(s,p); va_end(p); /* print formatted value */
  term_line_empty= s[strlen(s)-1]=='\n'; update_terminal();
}

void print_progress (char* s) @+{@; if (show_progress) print(s); }

void print_section_progress (void)
@+{@; if (show_progress) print("*%u",section_count); }


@* Index. % CWEAVE will add the index to this section.
