% wc: An example of CWEB by Silvio Levy and Donald E. Knuth
% Adapted to CWEB version 3.0 by Marc van Leeuwen -- CWI Amsterdam

\nocon % omit table of contents
\datethis % print date on listing

@* An example of \.{CWEB}.  This example presents the ``word count''
program from \caps{UNIX}, rewritten in \.{CWEB} to demonstrate
literate programming in \Cee.  The level of detail is intentionally
high, for didactic purposes; many of the things spelled out here don't
need to be explained in other programs.

The purpose of \.{wc} is to count lines, words, and/or characters in a list
of files. The number of lines is the number of newline characters in the
file. The number of characters is the file length in bytes.  A ``word'' is a
maximal sequence of consecutive characters other than white space,
containing at least one visible character.

@ Most \.{CWEB} programs share a common structure.  It's probably a good
idea to have one module that states this structure explicitly, even though
the elements could all be introduced in sections contributing to of the
unnamed module if they don't need to appear in any special order.

@c

@< Global variables @>@;
@< Functions @>@;
@< The main program @>@;

@ We must include the standard I/O definitions, since we want to send
formatted output to |stdout| and |stderr|. We also use the character
classification macro |isgraph| to detect visible characters.

@h <stdio.h>
@h <ctype.h>
@ The |status| variable tells the operating system if the run was
successful or not, and |prog_name| is used in case there's an error message to
be printed.

@d OK 0 /* |status| code for successful run */
@d usage_error 1 /* |status| code for improper syntax */
@d cannot_open_file 2 /* |status| code for file access error */

@< Global variables @>=
int status=OK; /* exit status of command, initially |OK| */
char *prog_name; /* who we are */

@ Now we come to the general layout of the |main| function. 

@< The main... @>=
main (int argc,char** argv)
{
  @< Variables local to |main| @>@;
  prog_name=*argv++; --argc; /* process program name */
  @< Set up option selection @>
  @< Process all the files @>
  @< Print the grand totals if there were multiple files @>
  exit(status);
}

@ If the first argument begins with a `\.-' the user is choosing
the desired counts and specifying the order in which they should be
displayed.  Each selection is given by the initial character
(lines, words, or characters).  For example, `\.{-cl}' would cause
just the number of characters and the number of lines to be printed.

We do not process this string now.  It will be used to control the
formatting at output time.

@< Var... @>=
int file_count; /* how many files there are */
char *which; /* which counts to print */

@ @< Set up o... @>=
which="lwc"; /* if no option is given, print all three values */
if (argc>0 && (*argv)[0] == '-') {@; which=&(*argv++)[1]; --argc; }
file_count=argc;

@ Now we scan the remaining arguments and try to open a file, if
possible.  The file is processed and its statistics are given.
We use a |do|~\dots~|while| loop because we should read from the
standard input if no file name is given.

@< Process... @>=
do {
  @< If a file is given try to open |*argv|; |continue| if unsuccesful @>
  @< Initialize pointers and counters @>
  @< Scan file @>
  @< Write statistics for file @>
  @< Close file @>
  @< Update grand totals @> /* even if there is only one file */
} while (++argv,--argc>0);

@ Here's the code to open the file. We use the low-level functions |open|,
|read|, and |close| that operate work file descriptors rather than with
|FILE|s. A special trick allows us to handle input from |stdin| when no name
is given.  Recall that the file descriptor to |stdin| is 0; that's what we
initialize our file descriptor to.

@< Variabl... @>=
int fd=0; /* file descriptor, initialized to |stdin| */

@~@d READ_ONLY 0 /* read access code for system |open| routine */

@< If a fi... @>=
if (file_count>0 && (fd=open(*argv,READ_ONLY))<0) {
  fprintf (stderr, "%s: cannot open file %s\n", prog_name, *argv);
@.cannot open file@>
  status|=cannot_open_file;
  --file_count;
  continue;
}

@ @< Close file @>=
close(fd);

@ We will do some homemade buffering in order to speed things up: Characters
will be read into the |buffer| array before we process them.
To do this we set up appropriate pointers and counters.

@d buf_size BUFSIZ /* \.{stdio.h}'s |BUFSIZ| is chosen for efficiency*/

@< Var... @>=
char buffer[buf_size]; /* we read the input into this array */
register char *ptr; /* the first unprocessed character in |buffer| */
register char *buf_end; /* the first unused position in |buffer| */
register int c; /* current character, or number of characters just read */
int in_word; /* are we within a word? */
long word_count, line_count, char_count; /* number of words, lines, 
    and characters found in the file so far */

@ @< Init... @>=
ptr=buf_end=buffer; line_count=word_count=char_count=0; in_word=0;

@ The grand totals must be initialized to zero at the beginning of the
program. If we made these variables local to |main|, we would have to
do this initialization explicitly; however, \Cee's globals are automatically
zeroed. (Or rather, ``statically zeroed.'') (Get it?)
@^Joke@>

@< Global var... @>=
long tot_word_count, tot_line_count, tot_char_count; /* total number of words, lines and chars */

@ The present module, which does the counting that is \.{wc}'s {\it raison
d'\^etre}, was actually one of the simplest to write. We look at each
character and change state if it begins or ends a word.

@< Scan... @>=
while (1) {
  @< Fill |buffer| if it is empty; |break| at end of file @>
  c=*ptr++;
  if (isgraph(c)) /* visible character */
  {@; if (!in_word) ++word_count, in_word=1; }
  else if (isspace(c))
  { in_word=0; /* |c| white space */
    if (c=='\n') ++line_count;
  }
}

@ Buffered I/O allows us to count the number of characters almost for free.

@< Fill |buff... @>=
if (ptr>=buf_end) {
  ptr=buffer; c=read(fd,ptr,buf_size);
  if (c<=0) break;
  char_count+=c; buf_end=buffer+c;
}

@ It's convenient to output the statistics by defining a new function
|wc_print|; then the same function can be used for the totals.
Additionally we must decide here if we know the name of the file
we have processed or if it was just |stdin|.

@< Write... @>=
wc_print(which, char_count, word_count, line_count);
if (file_count) printf (" %s\n", *argv); /* not |stdin| */
else printf ("\n"); /* |stdin| */

@ @< Upda... @>=
tot_line_count+=line_count;
tot_word_count+=word_count;
tot_char_count+=char_count;

@ We might as well improve a bit on \caps{UNIX}'s \.{wc} by counting the
files too.

@< Print the... @>=
if (file_count>1) {
  wc_print(which, tot_char_count, tot_word_count, tot_line_count);
  printf(" total in %d files\n",file_count);
}

@ Here now is the function that prints the values according to the
specified options.  The calling routine is supposed to supply a
newline. If an invalid option character is found we inform
the user about proper usage of the command. Counts are printed in
10-digit fields so that they will line up in columns.

@d print_count(n) printf("%10ld",n)

@< Fun... @>=
wc_print(char* which, long char_count, long word_count, long line_count)
{
  while (*which) 
    switch (*which++) {
    case 'l': print_count(line_count); break;
    case 'w': print_count(word_count); break;
    case 'c': print_count(char_count); break;
    default: if ((status & usage_error)==0) {
        fprintf (stderr, "\nUsage: %s [-lwc] [filename ...]\n", prog_name);
@.Usage: ...@>
        status|=usage_error;
      }
    }
}

@ Incidentally, a test of this program against the system \.{wc} command
on a SPARCstation showed that the ``official'' \.{wc} was slower. Furthermore,
although that \.{wc} gave an appropriate error message for the options
`\.{-abc}', it made no complaints about the options `\.{-labc}'!
Perhaps the system routine would have been better if its programmer had
been more literate?

@* Index.
Here is a list of the identifiers used, and where they appear. Underlined
entries indicate the place of definition. Error messages are also shown.
