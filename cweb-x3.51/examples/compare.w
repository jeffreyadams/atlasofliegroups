@* Comparing text files.
This is an entirely trivial program, that tests whether two text files are
equal, and if not so, points out the first point of difference.

@h <stdio.h>
@h <stdlib.h>
@c typedef char bool;

@ The outline of the program is simple. We read characters from both input
files into |c1| and~|c2| until the comparison is complete. Line and column
counts are maintained in |line| and~|col|.

@d left_margin 1 /* leftmost column number; change to 0 if you prefer */
@c
@<Functions@>@;
int main(int n, char** arg)
{ FILE *f1,*f2; /* the two input files */
  int c1,c2,col=left_margin; long line=1;
  @< Open the files |f1| and~|f2|,
     taking their names from the command line or from the terminal;
     in case of an error for which no recovery is possible,
     call |exit(1)| @>
  @< Search for first difference,
     leaving |c1!=c2| if and only if a difference was found @>
  @< Report the outcome of the comparison @>
  return 0;  /* successful completion */
}

@ The heart of the program is this simple loop.  When we reach the end of
one of the files, the files match if and only if the other file has also
reached its end. For this reason the test |c1==c2|, which requires
characters to be read from both files, must precede the test for file end;
when only one file ends, it is the former test which breaks the loop.

@< Search for first difference... @>=
while ((c1=getc(f1))==(c2=getc(f2)) && c1!=EOF)
  if (c1=='\n') {@; ++line; col=left_margin; } @+ else ++col;

@ When the first difference occurs at the end of one of the files, or at the
end of a line, we give a message indicating this fact.

@< Report... @>=
if (c1==c2) printf("Files match.\n");
else 
{ printf("Files differ.\n");
  if (c1==EOF || c2==EOF)
  @/{@;
    the_file(c1==EOF);
    printf("is contained in the other as initial segment.\n");
  }
  else if (c1=='\n' || c2=='\n')
  @/{@;
    the_file(c1=='\n');
    printf("has a shorter line number %ld than the other.\n",line);
  }
  else printf("First difference at line %ld, column %d.\n",line,col);
}

@ The function |the_file| starts a sentence about the first or second file,
depending on its boolean argument.

@<Functions@>=
void the_file(bool is_first)
@+{@; printf("The %s file ", is_first ? "first" : "second" ); }

@ There can be be zero, one or two command line arguments. If there are none,
the user is prompted to supply them, and if there are two these are taken as
the file names, prompting the user only in case a file could not be opened.
In case just one argument is present, the first file is assumed to be the
standard input, which does not have to be opened; in this case however we
will not read a file name from terminal in case the second file cannot be
opened.

@d read_mode "r"
@< Open... @>=
--n; ++arg; /* ignore ``argument'' 0, which is the program name */
if (n==0)
@/{@;
  open_file(&f1,"First file to compare", NULL);
  open_file(&f2,"Second file to compare", NULL);
}
else if (n==1)
{ f1=stdin;
  if ((f2=fopen(*arg,read_mode))==NULL)
    {@;  printf("Could not open file %s.\n",*arg); exit(1); }
}
else if (n==2)
{ open_file(&f1,"Give another first file", *arg++);
  open_file(&f2,"Give another second file", *arg);
}
else
{@; printf("No more than two command line arguments are allowed.\n");
  exit(1);
}

@ The function |open_file| will try to open the file |name| for reading, and
if this fails it will prompt for another file name until it has success. If
called with |name==NULL|, the function starts with prompting right away.

@<Functions@>=
void open_file(FILE** f,char* prompt,char* name)
{ char buf[80];
  if (name==NULL || (*f=fopen(name,read_mode))==NULL)
    do {@; printf("%s: ",prompt); fflush(stdout); scanf("%79s",buf); }@/
    while ((*f=fopen(buf,read_mode))==NULL);
}

@* Index.
