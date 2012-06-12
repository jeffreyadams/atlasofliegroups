% Adapted to CWEB version 3.0 by Marc van Leeuwen -- CWI Amsterdam
\noinx

@* Introduction. This file contains the program |wmerge|, which takes two
or more files and merges them according to the conventions of \.{CWEB}. We
use the routines of \.{CWEB} itself. The function |common_init| takes care of
processing command line arguments. Since the result of the merge will be
produced on the standard output, we prevent distraction as much as possible
by clearing flags |'h'| and |'p'| (for |'b'| it is not required since no
banner is produced anyway.

@h <stdio.h>
@h "../common.h" /* the header file for \.{CWEB}'s \.{common.w} */
@c
@< Prototype @>@;
main (int argc,char** argv)
{ common_init(argc,argv);
  flags['h']=flags['p']=0;
  reset_input();
  while (get_line())
    put_line();
  wrap_up();
}

@ This file should be linked together with the object file produced from
|"common.w"|, which is also used in both |CTANGLE| and |CWEAVE|.
That file defines the functions |common_init|, |reset_input|, |get_line|, and
|wrap_up|. There are however a number of functions that are required by that
compilation unit although they are not actually used; we define them with
trivial function bodies. Since the linker doesn't check types anyway we
don't specify any here either.

@c
void print_stats() @+ {}
void names_match () @+ {}
void init_module_name() @+ {}
void init_id_name () @+ {}


@ All that remains is to define |put_line| which is trivial. The external
variable |buffer| holds the characters read by |get_line|, up to |limit|,
and |loc| points to the next character to be read, i.e., after calling
|get_line| it points to |buffer[0]|.

@< Prototype @>= void put_line(void);
@~@c

void put_line(void)
{
  while (loc<limit) putchar(*loc++);
  putchar('\n');
}
