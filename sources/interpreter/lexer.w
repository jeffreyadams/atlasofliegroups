% Copyright (C) 2006 Marc van Leeuwen
% This file is part of the Atlas of Reductive Lie Groups software (the Atlas)

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

@* Introduction.
This file describes a lexical scanner that provides a layer situated between
the class |BufferedInput| that provides lines of input, and the parser that is
to receive a sequence of tokens. The functionality provided is the
recognition of tokens, in particular keywords and constants, and the
management of a table of identifiers, so that in the remainder of the program
identifiers can be represented by a small identification number rather than by
a string (with a conversion to string for output still being possible).

@h "buffer.h"


@ As usual the external interface is written to the header file associated to
this file.

@( lexer.h @>=
#ifndef LEXER_H
#define LEXER_H
@p@;
namespace atlas { namespace interpreter {

@< Class declarations @>@;
@< Declarations of exported functions @>@;
@< Declarations of static variables @>@;

}@; }@;
#endif

@ Here we include only our own header file \.{lexer.h}; all other include
directives will go to that file, as indicated in the previous section (this
might be a bit wasteful, but only system header files are concerned).

@c
#include "lexer.h"
using namespace std; @/
namespace atlas
{ namespace interpreter
   {
@< Definitions of class members @>@;
@< Function definitions @>@;
@< Definitions of static variables @>@;
   }@;
}@;



@* Tables recording identifiers. All distinct tokens recognised by the scanner
are stored as null-terminated strings in dynamically allocated memory; these
include all keywords, operators, and identifiers used, but not literal
constants or strings (since there is no point in remembering those). To
provide storage for such strings, we provide a simple class |String_pool|
which basically allows to copy any string to a safe place. When the string
pool is destructed, all character pointers it has given out become invalid.
Only one string pool is envisaged, and it its lifetime will only end at
termination of the program; still the following class is adapted to manage
strings in various pools with different lifetimes. A |String_pool| will be
implemented as a list of blocks containing the individual strings, from newest
to oldest. The blocks are usually of the size specified when the pool is
constructed, but if in rare cases a string larger than the block size is
requested, a single block of that size will be created (but |block_size| will
not change). To know how much remains in the current block, we keep a separate
record, in |chars_left|.


@< Class declarations @>=

class String_pool
{ public:
    static const size_t default_block_size=512;
    String_pool(size_t b=default_block_size); // create pool and set block size
    char* store(const char* s, size_t length); // store string of given length
    ~String_pool();
  private:
    String_pool(const String_pool&); // copying forbidden
    String_pool& operator=(const String_pool&); // assignment forbidden
@)
    size_t block_size; // block size of the pool
    char* start; // start of block in current usage
    char* point; // first character available
    size_t chars_left; // number of characters left in current block
    String_pool* prev; // start of list of previous blocks, for deallocation
};

@ The constructor for a |String_pool| sets the block size, but does not yet
allocate a block.

@< Definitions of class members @>=
String_pool::String_pool(size_t b)
: block_size(b),start(NULL),point(NULL),chars_left(0),prev(NULL)
{}

@ Whenever a identifier |str| of length |n| has to be stored in |pool|, we
shall call |pool.store(str,n)|. It will allocate |n+1| bytes and copy |str|
including the terminating null byte into it, returning a pointer to the
allocated string.

@< Definitions of class members @>=
char* String_pool::store(const char* s,size_t len)
{ size_t n= len>=block_size ? len+1 : block_size;
            // size of new block if needed
  if (len>=chars_left) // then the current block cannot store |s|, so allocate
  { if (start!=NULL) // initial block needs no backing up
    { String_pool* p=new String_pool(block_size);
      p->prev=prev; p->start=start; // only these fields need backing up
      prev=p; // link to old blocks, needed for destruction
    }
    start=point=new char[n]; // now we can copy to where |point| points
    chars_left=n;
  }
  strncpy(point,s,len); @+ char* result=point;
  chars_left-=len+1;
  point+=len;  *point++='\0';
  return result;
}

@ The destructor for |String_pool| recursively frees all blocks allocated.
This happens because |delete| calls the destructor for the |String_pool|
object pointed to by its argument |prev| before freeing its memory.

@< Definitions of class members @>=
String_pool::~String_pool(void)
@+{@; delete prev; delete[] start; }

@*1 The hash table. Our hash table will use an instance of |String_pool| for
storing its keys. The functionality of the hash table is quite simple, and it
hides the actual hashing. Upon construction the hash table is empty, and its
stores every distinct string that is passed to it, never removing one. It
matches the string against all strings it has seen before, and returns a
sequence number: either the number that was given to the same string when it
was first encountered, or the next unused number if the string was never seen
before. Methods are also provided for the inverse conversion, from a sequence
number to the corresponding string, and for determining the current number of
entries stored. This will allow working everywhere with sequence numbers to
represent identifiers.

@h<vector>
@< Class declarations @>=

class Hash_table
{ public: typedef short id_type; // type of value representing identifiers
private: // data members
    String_pool pool;
    id_type mod;  // hash modulus, the number of slots present
    std::vector<id_type> hash_tab;
    std::vector<const char*> name_tab;
public: // interface
    Hash_table(size_t init=initial_hash_mod
              ,size_t block_size=String_pool::default_block_size);
    id_type match(const char* s, size_t l) @+{@; return do_match(s,l,true); }
    id_type match_literal(const char* s) /* literal null-terminated string */
      {@; return do_match(s,strlen(s),false); }
    const char* name_of(id_type nr) const @+{@; return name_tab[nr]; }
    id_type nr_entries() const @+{@; return name_tab.size(); }
private: // auxiliary functions
    static const id_type initial_hash_mod=97;
    id_type hash(const char* s, size_t l) const; // auxiliary hash function
    id_type do_match(const char* s, size_t l, bool copy_string);
    size_t max_fill() const {@; return static_cast<id_type>(mod*3L/4); }
                    // maximum of entries in use
};

@ Here is the constructor for a hash table. By nature |hash_tab| is larger
than the number of entries currently stored; its size is always equal to~|mod|
(which is stored separately to speed up hash calculations). The initialiser
expression for |mod| takes care to ensure that |max_fill()<mod| always holds,
by rejecting ridiculously small values. Although |name_tab| will grow
incrementally, we prefer to reserve enough space for it to avoid reallocation
until rehashing is necessary. Therefore the allocated size of |name_tab| will
always have size |max_fill| which grows with (but lags behind) |mod|.

@< Definitions of class members @>=

Hash_table::Hash_table(size_t init,size_t block_size)
: pool(block_size)
, mod(init<2?2:init), hash_tab(mod), name_tab(0)
{@; name_tab.reserve(max_fill()); for(int i=0; i<mod; ++i) hash_tab[i]=-1;
}

@ The hash code of an identifier is determined by interpreting its characters
as digits of a number in base~256, and then taking the number modulo |mod|.
While the final result is a short value, intermediate results may be almost as
large as |mod<<8|, so we use a |long| value here.

The code below shows a bad property of the rigid type system of \Cpp: it is
virtually impossible to transform a pointer of type |(char*)| to one of type
|(unsigned char*)|, short of using a |reinterpret_cast| or (even worse)
calling the placement-|new| operator (as in |unsigned char p=new(s) unsigned
char@;|), even though this is a fairly innocent thing to do. Here we just want
to assure that we are doing arithmetic on non-negative values, even if some
characters would classify as negative |char| values (failing to do this could
cause the program the program to crash). And we cannot help having pointers to
(possibly) signed characters in the first place, because that is for instance
what string denotations give us (and what many standard functions require).
The simplest solution would be to convert |s| to |(const unsigned char*)|, and
everything would be safe. We dare not do such a conversion for fear of
losing our job, so instead we convert every character accessed via~|s|
laboriously into an unsigned value, which can be done by a standard
conversion; this accounts for two of the three |static_cast|s below
(the third is for documentation only).

@< Definitions of class members @>=

Hash_table::id_type Hash_table::hash
        (const char* s, size_t l) const // |s| a string of length |l>0|
{ long r=static_cast<unsigned char>(*s++)%mod;
  for ( --l; l>0 ; --l) r=(static_cast<unsigned char>(*s++)+(r<<8))%mod;
  return static_cast<id_type>(r);
}


@ When an identifier has been isolated in the input stream, a call to |match|
returns its index in |name_tab|, possibly after adding the identifier to the
tables. This small number is henceforth used to characterise the identifier.
The private function |do_match| does the real work; it is also called by
|match_literal| in case of keywords stored at initialisation time, for which a
null-terminated string is supplied that does not need to be copied to the
string pool. The parameter |copy_string| distinguishes the two uses.

Searching the hash table starts at the hash key and proceeds linearly through
the table until either a match or an empty slot is found; in the latter case
no matching entry is present, and we add one in the place of the empty slot.
This method is valid because we never remove an entry from the hash table. We
use here the fact that |max_fill()<mod| holds, so that there will be at least
one empty slot present, and it is not necessary to test for complete wrap
around the table. When a new entry would make the table too full, we extend
the size of the table and rehash. After this we could return
|do_match(str,l,copy_string)| recursively, but since we assume that the
too-full condition has been resolved, and we know that |str| is new, we can
avoid recursion and perform a simplified computation of the new hash code in
this case.

@h <cstring>
@< Definitions of class members @>=
Hash_table::id_type Hash_table::do_match
	(const char* str, size_t l, bool copy_string)
{ id_type i,h=hash(str,l);
  while ((i=hash_tab[h])>=0)
    if (std::strncmp(name_tab[i],str,l)==0 && name_tab[i][l]=='\0') return i;
          /* identifier found in table */
    else if (++h==mod) h=0; /* move past occupied slot */
  if (name_tab.size()>=max_fill())
  { @< Extend hash table @>
    h=hash(str,l); @+
    while (hash_tab[h]>=0) @+
      if (++h==mod) h=0; /* find free slot */
  }
  hash_tab[h]= name_tab.size(); /* this number represents the new identifier */
  name_tab.push_back(copy_string ? pool.store(str,l) : str);
  return hash_tab[h];
}

@ When the hash table becomes too full, we can throw away the |hash_tab|
because the |name_tab| contains all the information required; therefore we can
assign a larger vector to |hash_tab| rather than using |hash_tab.resize()|
which would copy the old entries. By performing the assignment via the |swap|
method (necessarily taken from the newly constructed value) we can avoid any
risk of copying of the new (empty) slots that an assignment might imply. The
code for inserting the indices of rehashed strings from |name_tab| is simpler
than in the basic matching loop, since we know that all those strings are
distinct; it suffices to find for each an empty slot. The sequence number
associated to each string does not change. We finally reserve enough space for
|name_tab| to last until the next time rehashing is necessary.

@< Extend hash table @>=
{ vector<id_type>(mod=2*mod-1).swap(hash_tab); // keep it odd
  for (int h=0; h<mod; ++h) hash_tab[h]=-1;
  for (size_t i=0; i<name_tab.size(); ++i)
  { int h=hash(name_tab[i],strlen(name_tab[i])); // rehash the string
    while (hash_tab[h]>=0) @+ if (++h==mod) h=0; /* find empty slot */
    hash_tab[h]=i;
  }
  name_tab.reserve(max_fill());
}

@ Since there is one central hash table, and other parts of the program must
have access to it, for instance to look up names of identifiers, we define a
static variable with a pointer to this hash table.

@< Declarations of static variables @>=
extern Hash_table* main_hash_table;

@~We initialise this variable to the null pointer; the main program will make
it point to the main hash table once it is allocated.

@< Definitions of static variables @>=
Hash_table* main_hash_table=NULL;

@*1 Identifier completion.
We define a completion function |id_completion_func| that will be used by the
\.{readline} library. This is done here since it deals with the hash table,
but strictly speaking this has nothing to do with lexical analysis. The
completion function will be installed in the main program, and it will be
called from the |readline| function (which is probably called by
|BufferedInput::getline|) if the user asks for it. The function prototype is
dictated by the \.{readline} library.

@< Declarations of exported functions @>=
extern "C"
char* id_completion_func(const char* text, int state);

@~The completion function will find and identifiers currently present in the
main hash table. This includes keywords, built-in functions and identifiers
introduced by the user, the latter possibly even by error (which is possibly
annoying for those who make many typos). It has some strange characteristics
that are dictated by the \.{readline} library. It must perform a loop that is
actually started \emph{outside} the function body, so the only possible way to
keep track of the loop state is using static variables. The |state| parameter
signals (by being~|0|) when a new loop starts, so in that case it is time to
(re-)initialise the static variables. A part of the loop can be picked up
inside the function body, namely the search for the next match to the supplied
prefix |text|. A tricky point is that once a partial match is found, we must
increment the static iterator before returning, since that |return| jumps out
of our local loop. For this reason we increment |i| right away while picking
its identifier from the hash table. Otherwise there are no other complications,
except having to produce a string allocated by |malloc|. If our local
loop terminates normally there are no more matches and we return |NULL| to
indicate that.

@h <cstdlib>
@< Function definitions @>=
extern "C"
char* id_completion_func(const char* text, int state)
{ static size_t l; static Hash_table::id_type i,n;
  if (state==0)
  { i=0; n=main_hash_table->nr_entries();
    l=std::strlen(text); // fix length for during search
  }
  while (i<n)
    // |i| is initialised above when |state==0|, and incremented below
  { const char* s=main_hash_table->name_of(i++);
      // get stored identifier and increment loop
    if (std::strncmp(text,s,l) == 0) // is |text| a prefix of |s|?
    { char* res=static_cast<char*>(std::malloc(std::strlen(s)+1));
        // we need a |malloc|ed string
      if (res==NULL)
	return NULL; // if memory gets full here, that's just too bad
      std::strcpy(res,s);
      return res;
    }
  }
  return NULL; /* if loop terminates, report failure */
}

@* The lexical analyser class.
%
We now come to the lexical analyser proper. Although only one lexical analyser
is envisaged, we shall define a class for it. The file \.{parser.tab.h}
contains definitions of constants defined by the parser and used in the code
below, but on its turn it uses types defined in \.{parsetree.h} which is
therefore loaded before it (we would like to have put an \&{\#include} of the
file \.{parsetree.h} into \.{parser.tab.h}, but do not know if or how this
could be arranged).

@h "parsetree.h"
@h "parser.tab.h"

@< Class declarations @>=
class Lexical_analyser
{ enum states @+ { initial, normal, ended };
@)BufferedInput& input;
  Hash_table& id_table;
  Hash_table::id_type keyword_limit; // first non-keyword identifier
  int nesting; // number of pending opening symbols
  char prevent_termination; // either |'\0'| or character requiring more input
  int comment_start, comment_end; // characters that start/end a comment
  states state; // to trigger special behaviour
  std::string file_name; // stores a file name for I/O redirection
public:
  Lexical_analyser(BufferedInput&, Hash_table&, const char**);
  int get_token(YYSTYPE *valp, YYLTYPE* locp);
  bool reset(); // get clean slate, return |false| if |getline()| fails
  void set_comment_delims (char c, char d) @+
          {@; comment_start=c; comment_end=d; }
  const char* scanned_file_name() const {@; return file_name.c_str(); }
private:
  void skip_space();
  char* scan_quoted_string();
};

@ Since there is one lexical analyser object, and other parts of the program
must have access to it, we define a static variable with a pointer to it.

@< Declarations of static variables @>=
extern Lexical_analyser* lex;

@~We initialise this variable to the null pointer; the main program will make
it point to the main hash table once it is allocated.

@< Definitions of static variables @>=
Lexical_analyser* lex=NULL;


@ Here is the constructor for the lexical analyser, which assumes that a
buffered input object and an empty hash table object have been previously
constructed, and are passed by reference. Currently it is called with a list
of keyword strings as final parameter, which keywords will be installed into
the hash table and determine the value of |keyword_limit|. There will probably
be a need to further parametrise the lexical analyser, if we do not want to
hard-code all lexical details into it (there is nothing wrong with that as
long as there is only one object of this class, but the class concept invites
us to envision some more flexible use). One such parametrisation is via the
|comment_start| and |command_end| characters, that if set will automatically
skip text enclosed between them (they may or may not be equal). In the unset
state they are set to an integer that cannot match any |char| value; we would
have like to use |EOF| defined in \.{ctype.h} here, but it is only guaranteed
to be non-|(unsigned char)|, and since using the type |(unsigned char*)| is
unwieldy, we use another value.

@< Definitions of class members @>=
Lexical_analyser::Lexical_analyser
  (BufferedInput& source, Hash_table& hash, const char** keywords)
: input(source),id_table(hash),nesting(0)
 ,prevent_termination('\0'),state(initial)
{ @< Install |keywords| into |id_table| @>
  comment_start=comment_end=256; // a non-|char| value
}

@ Keywords are identified by sequence number in the order by which they are
entered into |id_table|; therefore the order in the list |keyword| passed to
the constructor should match the numeric \.{\%token} values define in
\.{parser.y}. The actual code transmitted for keywords will be obtained by
adding the constant |FIRST_KEYWORD_CODE| to the value returned from the hash
table look-up.

@< Install |keywords| into |id_table| @>=
{ for (size_t i=0; keywords[i]!=0; ++i)
    id_table.match_literal(keywords[i]);
  keyword_limit=id_table.nr_entries();
}

@ The member function |reset| can be called to reset the lexical analyser,
discarding any remaining input on the current line and clearing the |nesting|
level. It is safe to call |reset| after a successfully executed command, which
will fetch a new line without getting any tokens, but |reset| should not be
called twice in succession, as this will discard the newly fetched line. The
result returned tells whether a fresh line was successfully obtained.

@< Definitions of class members @>=
bool Lexical_analyser::reset()
{@; nesting=0; state=initial; input.reset(); return input.getline(); }

@ Skipping spaces is a rather common activity during scanning; it is performed
by |skip_space|. When it is called, the first potential space character has
been shifted in, and when it returns a non-space character is left in
shifted-in position. Processing of newline-escaping backslashes and skipping
comments is also performed by |skip_space|. If |comment_start| is set and
|comment_end!='\n'|, then newlines will be skipped inside comments, but if
|comment_end=='\n'| then a closing newline will be considered for the
possibility of ending input. In the former case we change the input prompt by
pushing the comment character to warn the user that no action has been
performed. In case |prevent_termination| is set, that character is also pushed
input the prompt for the duration of skipping spaces.

In case end of input occurs one obtains |shift()=='\0'|. If this happens in
the main loop then we break from it like for any non-space character, and
although the following |input.unshift()| does nothing, the value |'\0'| should
reappear at the next call of |shift|. If end of input occurs during a comment,
then the comment is terminated due to the explicit test, and again the
persistence of |shift()=='\0'| will guarantee that the end of input eventually
turns up in non-skipping context.

@h<cctype>

@< Definitions of class members @>=
void Lexical_analyser::skip_space(void)
{ if (prevent_termination!='\0') input.push_prompt(prevent_termination);
  do
  { char c=input.shift();
    if (isspace(c))
    { if (c=='\n' && prevent_termination=='\0' && nesting==0) break;
         // newline is not skipped if some command could end here
      continue; // all other space is skipped
    }
    if (c==comment_start)
    { input.push_prompt(comment_start);
      do c=input.shift(); while (c!=comment_end && c!='\0');
      input.pop_prompt();
      if (c=='\n') input.unshift(); // reconsider newline character
    }
    else break; // non-space and non-comment character
  } while(true);
  if (prevent_termination!='\0') input.pop_prompt();
  input.unshift(); // prepare to re-read character that ended space
}

@ Another auxiliary method is used for scanning a quoted string; it should be
called when an initial double-quote character has been recognised, and after
scanning the string returns a pointer to the designated string (with quotes
and escapes removed), null terminated and allocated by |new[]|. We currently
use a simple model for strings. They should be contained in a single line, and
the only escapes used are the doubling of double-quote characters. We copy the
string while reducing doubled double-quote characters to single ones.


@< Definitions of class members @>=
char* Lexical_analyser::scan_quoted_string()
{ const char* start=input.point(),*end;
  int nr_quotes=0; // number of escaped quotes
  do
  { char c;
    do c=input.shift(); while (c!='"' && c!='\n' && c!='\0');
    if (c!='"')
    { input.unshift(); end=input.point();
      int l0,c0,l1,c1;
      input.locate(start,l0,c0); input.locate(end,l1,c1);
      input.show_range(cerr,l0,c0,l1,c1);
      cerr << "Closing string denotation.\n";
      break;
    }
    else if ((c=input.shift())!='"')
    {@; input.unshift(); end=input.point()-1; break; }
    else ++nr_quotes; // for doubled quotes, continue
  } while (true);
  size_t len=end-start-nr_quotes;
  char* s=new char[len+1];
  while (start<end) // copy characters, undoubling doubled quotes
    if ((*s++=*start++)=='"') ++start;
  *s='\0'; return s-len;
}


@*1 The main scanning routine.
%
Now we come to the function |get_token| that will actually be called by the
parser and return a token to it. The token is an integral value that is either
a character code in the case of single character tokens, or a token value
defined in \.{parser.tab.h}. The null token value has the special meaning of
signalling the end of the input to be recognised by the parser, which it wants
to see before returning successfully. Since our parser is written to recognise
individual commands rather than the full input source, we take measures to
send a null token after each command, even though the input buffer does
not signal end of input.

If it should happen that the input buffer does signal end of input by making
|input.shift()=='\0'|, then this will be transmitted to the parser as the end
of a command by the code below, by the general mechanism that characters that
are not specifically recognised are transmitted with a token value equal to
the character value. But this is really a marginal case, since the input
buffer guarantees that end of input cannot occur in the middle of a line;
normally the previous line will therefore be completely processed, and the
end-of-input condition is signalled by a failure of the call of the |reset|
method of the lexical analyser rather than by |input.shift()=='\0'| occurring
in its |getline| method. The only way the latter can happen is if the
preceding newline character was ignored by |skip_space|, due to
|prevent_termination| or |nesting|.

The way in which we arrange to signal the end of a command from |get_token| is
by sending \emph{two} successive tokens, a |'\n'| followed by a null token.
This strange setup is a work-around of a nasty property of the parsers
generated by \.{bison}, namely that on one hand the required terminal null
token is not representable in syntax rules (it should follow input matching
the target symbol, but is not part of its production), while on the other hand
the parser will perform any unambiguous reduction of the input it has seen,
even if the lookahead token is wrong for that reduction (so a syntax error
will be reported \emph{after} the reduction is performed). This means that if
we would define the syntax of commands normally, and if such a command would
occur correctly in the input but followed by extra tokens, then the reduction
recognising and executing the command would be performed \emph{before}
flagging a syntax error for the spurious tokens, which is undesirable.
Therefore we include into the production for any command the |'\n'| token that
is guaranteed by our lexical analyser to precede the terminating null token,
thus ensuring that the command will only be executed if it constitutes the
complete user input.

(One may imagine that premature reduction could cause problems in other
positions as well, but this is by far the most irritating case, and one where
no grammar rewriting can help.)

The implementation of this two-token termination reporting is by setting
|state=ended| when the token |'\n'| is returned, and testing this condition as
the first action in |get_token|, sending the following null token if it holds.
The |state| variable is also used to allow the scanner to behave differently
while scanning the very first token of a command, when |state==initial| will
hold; once a token is scanned, |state| is set to |normal| unless scanning had
set it to |ended|.

Besides returning a token code, each token defines a value for
|prevent_termination|. Since the previous value is only used in |skip_space|,
we set this variable to its most common value |'\0'| after the call to that
function, so that only for those tokens that cannot end a command we have to
set that value explicitly.

The code below also deals with setting the fields of the |locp| to delimit the
token scanned.

@< Definitions of class members @>=
int Lexical_analyser::get_token(YYSTYPE *valp, YYLTYPE* locp)
{ if (state==ended) {@; state=initial; return 0; } // send end of input
  skip_space(); prevent_termination='\0';
  input.locate(input.point(),locp->first_line,locp->first_column);
  int code; char c=input.shift();
  if (isalpha(c)) @< Scan an identifier or a keyword @>
  else if (isdigit(c)) @< Scan a number @>
  else @< Scan a token starting with a non alpha-numeric character @>
  input.locate(input.point(),locp->last_line,locp->last_column);
@/if (state==initial) state=normal; return code;
}

@ Everything that looks like an identifier is either that or a keyword. In any
case we start with looking it up in the table, and then the numeric value of
the code returned will allow us to discriminate the possibilities. In this
scanner we cannot yet handle multiple distinct keywords that should scan as
the same category because of a similar syntactic role, although this is
probably desirable.

@< Scan an identifier or a keyword @>=
{ const char* p=input.point@[()@]-1; // start of token
  do c=input.shift(); while(isalpha(c) || isdigit(c) || c=='_');
  input.unshift();
  Hash_table::id_type id_code=id_table.match(p,input.point()-p);
  if (id_code>=keyword_limit) {@; valp->id_code=id_code; code=IDENT; }
  else
  { code=QUIT+id_code;
    switch(code)
    {
      case LET: ++nesting; input.push_prompt('L'); break;
      case IN: --nesting; input.pop_prompt(); prevent_termination='I'; break;
      case WHATTYPE: prevent_termination='W'; break;
    }
  }
}

@ Numbers do not include a sign, and are always scanned as |unsigned long|
values.

@< Scan a number @>=
{ const char* p=input.point@[()@]-1; // start of token
  do c=input.shift(); while(isdigit(c));
  input.unshift();
  const char* end=input.point();
  unsigned long val=*p++-'0'; @+  while (p<end) val=10*val+(*p++-'0');
  valp->val=val; code=INT;
}

@ After splitting off the alphanumeric characters, the scanner plunges into a
large switch on the value of the look-ahead character~|c|. We increase and
decrease nesting on obvious grouping characters, and for certain characters we
store them into |prevent_termination| to prevent a newline to be interpreted
as command termination. If in the |initial| state we find a |'>'| character,
then we prepare for output redirection.

@< Scan a token... @>=
{ switch(c)
  {      case '"': @< Scan a string denotation @>
  break; case '(':
         case '{':
         case '[': ++nesting; input.push_prompt(c); code=c;
  break; case ')':
         case '}':
         case ']': --nesting; input.pop_prompt(); code=c;
  break; case '<':
         case '>':
         if (state==initial)
         { code= c=='<' ? FROMFILE :
                 input.shift()=='>' ? ADDTOFILE :
                 (input.unshift(),TOFILE) ;
	   @< Read in |file_name| @> break;
         }
    @/// |else| {\bf fall through}
         case '=':
         case '+':
         case '-':
         case '*':
         case '%':
         case '^':
         case ':': prevent_termination=c; code=c;
  break; case '/': prevent_termination=c;
    code= input.shift()=='%' ? DIVMOD : (input.unshift(),'/');
  break; case '\n': state=ended; // and {\bf fall through}.
         default: code=c;
  }
}

@ We hand to the parser a string denotation expression in |valp|, rather than
the |(char *)| value returned by |scan_quoted_string()|, thus performing some
work that is usually left to the parser. The reason for this is that it avoids
having to extend union of possible token values with a variant for |(char *)|
and to define a corresponding destructor for in case the parser should decide
to drop the token because of a syntax error. In the current situation the
token value is a (string denotation) expression, and its destruction is
handled by |destroy_expr|.

@< Scan a string... @>=
{@; valp->expression=make_string_denotation(scan_quoted_string());
  code=STRING;
}

@ Since file names have a different lexical structure than identifiers, they
are treated separately in the scanner; moreover since at most one file name
can occur per command, we store it in a field |file_name| of the lexical
scanner reserved for that purpose. We allow a file name to be either a
sequence of non-space characters not starting with a quote, or a quoted
string.

@< Read in |file_name| @>=
if ((skip_space(),c=input.shift())=='"')
@/{@; char* s=scan_quoted_string(); file_name=s; delete[] s; }
else
@/{@; file_name="";
    while (!isspace(c)){@; file_name+=c; c=input.shift(); }
    input.unshift();
}

@* Index.

% Local IspellDict: british
