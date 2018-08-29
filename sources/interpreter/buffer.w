% Copyright (C) 2006-2016 Marc van Leeuwen
% This file is part of the Atlas of Lie Groups and Representations (the Atlas)

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


\def\emph#1{{\it #1\/}}

@* Buffered input.
This file defines some classes used at the input level below that of the
lexical analyser. First it provides classes |String_pool| and |Hash_table| for
the management of collections of identifiers, so that in the remainder of the
program identifiers can be represented by a small identification number rather
than by a string (with a conversion to string for output still being possible).

Next it provides a simple class |BufferedInput|, which provides an interface
to input streams that can be used by the lexical analyser. Its purpose is to
deal with details that should not influence the lexical analysis, such as the
question whether input comes from a file or from standard input, in the latter
case whether the GNU \.{readline} library is called when a new line of input
is required, and the handling of escaped newlines. Since the input is copied
to a line buffer in any case, we can also provide character pointers into that
line buffer; thus one can read identifiers and constants that are contained in
a single line by first scanning characters until the end of the token, and then
extracting the string between the initial and final position. This avoids
having to copy the characters to a separate string during the initial scan
(but admittedly this will only be practical if the token characters can be used
``as is'', in other words there are no escaped or encoded characters). The
buffer also contains a line counter, so that the lexical analyser need not
worry about that (and indeed it could not, since escaped newlines are
transparently handled, and never seen by the lexical analyser).

We have included functionality to open a subsidiary input file upon an
existing buffered input object, thus creating a stack of files from which
input will be read; exhausted files will be transparently popped, returning
control to the file active at the point they were opened.

@( buffer.h @>=

#ifndef BUFFER_H
#define BUFFER_H

@< Includes needed in the header file @>@;


namespace atlas
{ namespace interpreter
  {
    @< Class declarations @>@;
    @< Inline function definitions @>@;
    @< Declarations of static variables @>@;
    @< Declarations functions used but defined elsewhere @>@;
  }@;
}@;
#endif

@ Our main file will provide the class implementations.

@h "buffer.h"
@c

namespace atlas
{ namespace interpreter
  {
@< Definitions of class members @>@;
@< Definitions of static variables @>@;
  }@;
}@;



@* Tables for storing identifiers.
%
All distinct tokens recognised by the scanner will be stored as
null-terminated strings in dynamically allocated memory; these include all
keywords, operators, and identifiers used, but not literal constants or
strings (since there is no point in remembering those). To provide storage for
such strings, we provide a simple class |String_pool| which basically allows
to copy any string to a safe place and use pointers to those copies henceforth
(we probably would have used a different approach if we had been more
acquainted with the standard library when we wrote this code). Only when the
string pool is destructed do the character pointers it has given out become
invalid. A string pool will be incorporated into the |Hash_table| class, and
their number and lifetimes will be determined by that usage; in particular the
|String_pool| class will be adapted to manage strings in various pools with
different lifetimes. A |String_pool| basically holds a pointer |start| to a
block of characters it currently distributes from (which if non null is
obtained from a call of |new[]|), a pointer |point| to the first among these
characters that is free for future distribution, and a count |chars_left| of
such characters. When a request cannot be served from the current block,
|start| will be set to fresh memory obtained from |new[]| but only after
saving the previous value to |prev|, effectively a linked list of old blocks
that the destructor of the |String_pool| will all free, as well as the current
block pointed to be |start|.


@< Class declarations @>=

class String_pool
{ public:
    static const size_t default_block_size=512;
    String_pool(size_t b=default_block_size); // create pool and set block size
    char* store(const char* s, size_t length); // store string of given length
    ~String_pool() @+{@; delete start; delete prev; }
  private:
    String_pool(const String_pool&); // copying forbidden
    String_pool& operator=(const String_pool&); // assignment forbidden
@)
    const size_t block_size; // block size of the pool
    char* start; // start of block in current usage
    char* point; // first character available
    size_t chars_left; // number of characters left in current block
    struct backup
    { const char* const block_start; @+
      const backup* const link;
    public:
      backup (char* start, backup* prev) : block_start(start), link(prev) @+{}
      ~backup ();
    };
    backup* prev; // list of previous blocks, for deallocation
};

@ The constructor for a |String_pool| sets the block size, but does not yet
allocate a block.

@< Definitions of class members @>=
String_pool::String_pool(size_t b)
: block_size(b),start(nullptr),point(nullptr),chars_left(0),prev(nullptr)
{}

@ Whenever a identifier |str| of length |n| has to be stored in |pool|, we
shall call |pool.store(str,n)|. It will allocate reserve |n+1| bytes, either
previously allocated or by allocating at this point, and copy |str|
including a terminating null byte into those bytes; it then returns a
pointer to those bytes.

@< Definitions of class members @>=
char* String_pool::store(const char* s,size_t len)
{ size_t n= std::max(len+1,block_size); // size of new block if needed
  if (len>=chars_left) // then the current block cannot store |s|, so allocate
  { if (start!=nullptr) // initial block needs no backing up
      prev=new backup(start,prev); // push pointers to chain, for destruction
    start=point=new char[n]; // now we can copy to where |point| points
    chars_left=n;
  }
  char* result=strncpy(point,s,len);
  point+=len;  *point++='\0'; // write |'\0'| whether or not present in |s|
  chars_left-=len+1;
  return result;
}

@ The destructor for |String_pool::backup| recursively frees all blocks of
characters it accesses via |block_start| pointers, as well as all the records
for the |backup| structures themselves. The recursive nature comes from the
fact that |delete link| calls the destructor |~backup| (an indirect recursion)
for the object pointed to by |link| before freeing the memory of that object.

@< Definitions of class members @>=
String_pool::backup::~backup() @+{@; delete[] block_start; delete link; }

@*1 Hash tables.
%
Our hash tables will use an instance of |String_pool| for storing its keys.
The functionality of hash tables is quite simple, and hides the actual
hashing. Upon construction the hash table is empty, and its stores every
distinct string that is passed to it, never removing any of them. It matches
the string against all strings it has seen before, and returns a sequence
number: either the number that was given to the same string when it was first
encountered, or the next unused number if the string was never seen before.
Methods are also provided for the inverse conversion, from a sequence number
to the corresponding string, and for finding the current number of entries
stored. This will allow working everywhere with sequence numbers to represent
identifiers.

@h<vector>
@< Class declarations @>=

typedef unsigned short id_type; // type of value representing identifiers

class Hash_table
{
public:
  static const id_type empty; // value (|~0u|) reserved for empty slots
private: // data members
  String_pool pool;
  id_type mod;  // hash modulus, the number of slots present
  std::vector<id_type> hash_tab;
  std::vector<const char*> name_tab;
@)
public: // interface
  Hash_table(size_t init=initial_hash_mod
            ,size_t block_size=String_pool::default_block_size);
@)
  id_type nr_entries() const @+{@; return name_tab.size(); }
  const char* name_of(id_type nr) const @+{@; return name_tab[nr]; }
  bool knows(const char* name) const;
@)
  id_type match(const char* s, size_t l) @+{@; return do_match(s,l,true); }
  id_type match_literal(const char* s) /* literal null-terminated string */
    {@; return do_match(s,std::strlen(s),false); }
private: // auxiliary functions
  static const id_type initial_hash_mod=97;
  id_type hash(const char* s, size_t l) const; // auxiliary hash function
  id_type do_match(const char* s, size_t l, bool copy_string);
  size_t max_fill() const {@; return static_cast<id_type>(mod*3L/4); }
                  // maximum of entries in use
};

@ Since we used the |vector| template and the |strlen| function, we must make
sure its header is included at the proper place.

@< Includes needed in the header file @>=
#include <vector>
#include <cstring>

@ Let's get that static class constant defined before we forget it, since a
temporary may be need to be created when it is passed as (constant) reference
to |id_type| argument, so the value supplied at definition does not suffice.
It is amazing how \Cpp\ makes it a pain to define a simple compile time
constant. Note that the common |enum| trick won't work here, since
|Hash_table::id_type| is a |short| type, so we definitely don't want |empty|
to have type |int| in comparisons. The |static_cast| is necessary to ensure we
use the |unsigned short| overload of the bitwise-complement operator |~|,
because there is no such thing as a short integer denotation; thus we avoid a
warning message about truncating a large unsigned constant to a shorter value,
even though that is well-defined behaviour. (The perplexed reader might be
reassured, or not, to learn that writing just |-1| as right hand side would
also get the desired result, and without provoking any warning message.)

@< Definitions of class members @>=

const id_type Hash_table::empty = ~static_cast<id_type>(0);

@ Here is the constructor for a hash table. By nature |hash_tab| is larger
than the number of entries currently stored; its size is always equal to~|mod|
(which is stored separately to speed up hash calculations). The initialiser
expression for |mod| takes care to ensure that |max_fill()<mod| always holds,
by rejecting ridiculously small values. Although |name_tab| will grow
incrementally, we prefer to reserve enough space for it to avoid reallocation
until rehashing is necessary. Therefore the allocated size of |name_tab| will
always have size |max_fill()|, which grows with (but lags behind) |mod|.

@< Definitions of class members @>=

Hash_table::Hash_table(size_t init,size_t block_size)
: pool(block_size)
, mod(init<2?2:init), hash_tab(mod,empty), name_tab(0)
{@; name_tab.reserve(max_fill());
}

@ The hash code of an identifier is determined by interpreting its characters
as digits of a number in base~256, and then taking the number modulo |mod|.
While the final result is a short value, intermediate results may be almost as
large as |mod<<8|, so we use an |unsigned long| value here.

The code below shows a bad property of the rigid type system of \Cpp: even
though this is a fairly innocent thing to do, it is virtually impossible to
transform a pointer of type |(char*)| to one of type |(unsigned char*)|, short
of using a |reinterpret_cast| or (even worse) calling the placement-|new|
operator (as in |unsigned char* p=new(s) unsigned char@;|). Here we just want
to assure that we are doing arithmetic on non-negative values, even if some
characters would classify as negative |char| values (failing to do this could
cause the program to crash). And we cannot help having pointers to (possibly)
signed characters in the first place, because that is for instance what string
denotations give us (and what many standard functions require). The simplest
solution would be to cast |s| to |(const unsigned char*)|, and everything
would be safe. We dare not do such a conversion for fear of losing our job, so
instead we convert every character accessed via~|s| laboriously into an
unsigned value, which can be done by a standard conversion; this accounts for
two of the three |static_cast|s below (the third is for documentation only).

The code below is safe if either |l>0| (as will always be the case for
identifiers) or if |l==0| and |*s=='\0'| in which case it will return |0| as
hash code (this may occur when the user rather foolishly attempts
to open a file with an empty name).

@< Definitions of class members @>=

id_type Hash_table::hash (const char* s, size_t l) const
  // |s| a string of length |l>0|
{ unsigned long r=static_cast<unsigned char>(*s++)%mod;
  while (l-->1) // one less, because we already read one character
    r=(static_cast<unsigned char>(*s++)+(r<<8))%mod;
  return static_cast<id_type>(r);
}

@ The following accessor method is a simplified version of |do_match| defined
below; it looks up whether a string is present, but does not attempt to insert
it if it is absent. It should not be used if the intention is to insert the
string anyway; in that case one should first record the old size of the table,
match the string and test whether the returned value is less that the old size
(which indicates that the identifier was already known). Indeed the code below
is only used with file names, in an optimisation that avoids opening a file if
the name corresponds to one already successfully opened before (in that
application it is too soon to know whether the name should be added, in case
it is new).

Given the fact that the hash table can only grow, we can be sure upon finding
an empty slot that the name searched is absent from the table. In case of
finding a non-empty slot, an exact name match will result in a positive
answer, while a failure to match leads to trying the next slot. Since it is
guaranteed that the table has at least one empty slot at all times, the loop
will certainly terminate.

@h <cstring>
@< Definitions of class members @>=
bool Hash_table::knows (const char* name) const
{ size_t l = std::strlen(name);
  id_type i,h=hash(name,l);
  while ((i=hash_tab[h])!=empty) // search breaks off at first empty slot
    if (std::strncmp(name_tab[i],name,l)==0 && name_tab[i][l]=='\0')
      return true; // identifier found in table
    else
      if (++h==mod) // move past occupied slot
        h=0; // possibly wrapping around
  return false;
}

@ When an identifier is passed to |match|, this method (which calls |do_match|
with |copy_string==true|) returns its index in |name_tab|, possibly after
adding the identifier to the tables. This small number is henceforth used to
characterise the identifier. The private function |do_match| can also be
called by |match_literal| in case of keywords stored at initialisation time,
for which a null-terminated string is supplied that does not need to be copied
to the string pool. The parameter |copy_string| distinguishes the two uses.

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

@< Definitions of class members @>=
id_type Hash_table::do_match (const char* str, size_t l, bool copy_string)
{ id_type i,h=hash(str,l);
  while ((i=hash_tab[h])!=empty)
    if (std::strncmp(name_tab[i],str,l)==0 && name_tab[i][l]=='\0')
      return i; /* identifier found in table */
    else if (++h==mod)
      h=0; /* move past occupied slot */
  if (name_tab.size()>=max_fill())
  { @< Extend hash table @>
    h=hash(str,l); @+
    while (hash_tab[h]!=empty) @+
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
|name_tab| to last until the next time rehashing is necessary. This is really
unnecessary, since |name_tab| can very well determine for itself when it needs
reallocation; the current code ensures that |name_tab| is reallocated exactly
when |hash_tab| grows.

@< Extend hash table @>=
{ std::vector<id_type>(mod=2*mod-1,empty).swap(hash_tab); // make length odd
  for (size_t i=0; i<name_tab.size(); ++i)
  { int h=hash(name_tab[i],std::strlen(name_tab[i])); // rehash the string
    while (hash_tab[h]!=empty) @+
      if (++h==mod)
        h=0; /* find empty slot */
    hash_tab[h]=i;
  }
  name_tab.reserve(max_fill());
}

@* The class for buffered input.
%
We come now to the definition of the |BufferedInput| class. Since its
declaration requires some types to have been defined, we must first make sure
that is done.

Input will be collected from one or more streams (a terminal input stream and
one or more files), and these will be accessed through |std::istream| values.
The copy constructor for that class is private, so there is no question of
actually containing data members of that type; they will always be handled by
reference or by pointer. As a consequence it should suffice in the header file
to know that |std::istream| is a class, which is done by including \.{iosfwd},
and include \.{iostream} in the implementation. The actual pointer to input
streams will be contained in an |input_record| structure that we shall discuss
in more detail later.

@h <iostream>

@< Includes needed... @>=
#include <iosfwd>

@~Here is the main part of the class definition, giving the principal user
interface. The methods |shift|, |unshift|, |eol| and |getline| give the basic
functionality of a line-buffered input stream. The method |push_file| opens a
new level of input, taking control immediately, and until it is exhausted
(although it can itself be interrupted by further calls of |push_file|). Such
input files should be opened after processing a complete line from the current
file, since calling |push_file| causes the current line of input to be
abandoned and forgotten. Another method |close_includes| is provided to allow
aborting all includes in case of errors. The private method |pop_file| is
called when an additional input stream is exhausted.

@< Class declarations @>=
class BufferedInput
{ typedef char* (*rl_type)(const char* );
     // |typedef| avoids inclusion of \.{readline} headers
  typedef void (*add_hist_type) (const char* );  // this one too
  @< Define |struct input_record@;| @> // a local class
  @< Define |struct file_stack@;| @> // another one, using the previous

@)
  @< Data members of |BufferedInput| @>@;
@)
 public:
    BufferedInput (std::istream& s);
        // associate line buffer to raw input stream
    BufferedInput (const char* prompt
                  , rl_type rl=nullptr, add_hist_type=nullptr @|
		  , const char* prompt2="> "
                  , const char* def_ext=".at"@|);
        // use |stdin|, maybe with readline
@)
    char shift (); // inspect a new character
    void unshift (); // back up so last character will be reconsidered
    bool eol () const; // end of line: |true| if |shift| would return a newline
    bool getline ();
         // fetch a new line to replace current one; |true| if successful
    bool push_file (const char* file_name,bool skip_seen);
      // |true| if successful
    void close_includes();
    @< Other methods of |BufferedInput| @>@;
private:
    void pop_file();
};

@ There will be one main input buffer, and other parts of the program must
have access to it, for instance to call the |push_file| method; therefore we
define a static variable with a pointer to this input buffer.

@< Declarations of static variables @>=
extern BufferedInput* main_input_buffer;

@~We define |main_input_buffer| here (it must be defined somewhere, and doing
so in \.{main.w} where it is mostly used would be awkward since there is no
header file~\.{main.h}). We initialise this variable to the null pointer; the
main program will make it point to the main input buffer once it is allocated.

@< Definitions of static variables @>=
BufferedInput* main_input_buffer=nullptr;

@ In our implementation we use the classes |std::string| and |atlas::BitMap|.

@< Includes needed in the header file @>=
#include <string>
#include "../Atlas.h" // for common |using| declarations
#include "bitmap.h" // so that the type |BitMap| is complete
#include "sl_list.h" // used for stack of input files

@~At construction a |BufferedInput| object will fix a reference |base_stream| to
the input stream to be used when no additional input files are open. Furthermore
it stores the last line read in its |line_buffer| field, and the pointer~|p|
points to the next character to be produced by |shift()|. In case this buffer is
used for terminal input, a prompt string is stored and will be printed every
time the program requests more input; this is the |prompt| pointer so the
lifetime of the value provided at construction should contain that of the
|BufferedInput| object created. We also store a secondary prompt |prompt2| and
provide a dynamic string |temp_prompt| that together can be used to temporarily
change the prompt to indicate that some input needs completion. In case of
terminal input, the input may be pre-processed by the function |readline|, which
is typically going to be the function of that name from the GNU~\.{readline}
library; then one may also want to use another function stored in |add_hist| to
store lines for interactive retrieval, typically the |add_history| function from
the history library. The member |line_no| holds the number of the current line,
or more precisely of the first of |cur_lines| lines that are currently read into
|line_buffer| (if more than one, they were joined by escaped newlines). We
record the length of the last prompt printed in |prompt_length| for the purpose
of pointing to tokens. Finally there is a stack |input_stack| of active input
files together with their file names, tables |input_files_seen| and
|input_files_completed| that together record the status (unread, partially read,
fully read) of all file names, and a pointer |stream| that tracks the current
input stream.

@< Data members... @>=
std::istream& base_stream;
std::string line_buffer; // current line
const char* p; // buffer pointer
const char *const prompt, *const prompt2; // primary and secondary prompt
const char* def_ext; // default file extension
std::string temp_prompt; // varying addition to secondary prompt
const rl_type readline; // readline function
const add_hist_type add_hist; // history function
unsigned long line_no; // current line number
int cur_lines,prompt_length; // local variables
file_stack input_stack;
 // active input streams
Hash_table input_files_seen; // files successfully opened at least once
BitMap input_files_completed; // marks those files that were read without error
std::istream* stream; // points to the current input stream

@ A line is ended if |p| is either null (indicating an absence of any line
buffered, as is initially the case), or points to a null character (the one
terminating the line buffer, which always follows a newline character).

@< Inline fun... @>=
inline bool BufferedInput::eol() const @+{@; return p==nullptr or *p=='\0'; }

@ There are two constructors, one for associating an input buffer to some
(raw) |istream| object (which may represent a disk file or pipe), another for
associating it to probably interactive input from |stdin|. A prompt and
readline function only apply to the second case and are set to |nullptr| in the
first case, which will disable certain interactions when fetching new lines; a
null pointer may also be passed explicitly in the second case to obtain this
disabling. Constructing the class does not yet fetch a line.

@< Definitions of class members @>=
BufferedInput::BufferedInput (std::istream& s)
@/:base_stream(s)
@/,line_buffer()
@/,p(nullptr)
@/,prompt(nullptr)
@/,prompt2(nullptr)
@/,def_ext(nullptr)
@/,temp_prompt(nullptr)
@/,@|readline(nullptr)
@/,add_hist(nullptr)
@/,line_no(1)
@/,cur_lines(0)
@/,input_stack()
@/,input_files_seen()
@/,input_files_completed()
@/,stream(&base_stream)
@/{}
@)
BufferedInput::BufferedInput
(const char* pr, rl_type rl, add_hist_type ah,const char* pr2,const char* de)
@/:base_stream(std::cin)
@/,line_buffer()
@/,p(nullptr)
@/,prompt(pr)
@/,prompt2(pr2)
@/,def_ext(de)
@/,temp_prompt("")
@/,@|readline(rl)
@/,add_hist(ah)
@/,line_no(1)
@/,cur_lines(0)
@/,input_stack()
@/,input_files_seen()
@/,input_files_completed()
@/,stream(&base_stream)
@/{}


@*1 Input from auxiliary files.
%
We shall define the class of the records making up the |input_stack|. They
will contain an |std::ifstream| field, so we need to include the corresponding
header file.
@< Includes needed... @>=

#include <fstream>

@~For every currently open \emph{auxiliary} input stream (necessarily a file
stream), we store a pointer to the |ifstream| object, the file name (for error
reporting) and the line number of \emph{the stream that was interrupted} by
opening this auxiliary file (the line number in the file itself will be held
in the input buffer, as long as it is not interrupted).

It used to be the case that |input_record| stored a pointer to |ifstream|
rather than such a structure itself, which was necessary because the
|input_stack| was implemented as a |std::vector|, and copying of |ifstream|
values is not possible (and this was before move semantics was made possible).
Now that the stack is a list, and records can be constructed in place, copying
is no issue, and they might as well contain an |ifstream| value; this takes
care of closing and cleaning up whenever stack records are popped or
destructed (in exception handling).

@< Define |struct input_record@;| @>=
struct input_record
{ std::ifstream f_stream; // the actual stream record
  id_type name; // identifies the file name for |stream|
  unsigned long line_no; // this refers to the older input stream!
@)
  input_record(BufferedInput&, const char* file_name);
  input_record@[(const input_record& rec) = delete@]; // cannot copy
  input_record& operator=@[(const input_record& rec) = delete@]; // or assign
};

@ Our file stack is essentially a |std::stack| of |input_record| using
|containers::mirrored_sl_list| as container (which allows for the element type
|input_record| the is not copy-assignable). Using |containers::mirrored_sl_list|
rather than |containers::mirrored_simple_list| allows us to directly call |size|
for our file stack, which will end up calling |containers::sl_list::size|.
However in one place we also need to iterate over the input stack (to check for
files currently being read), for which we add |ctop| and |cbottom| methods to
produce (int that order) a range to traverse for visiting all active
|input_record|s.

@< Define |struct file_stack@;| @>=
struct file_stack
: public std::stack<input_record,containers::mirrored_sl_list<input_record> >
{ typedef containers::simple_list<input_record>::weak_const_iterator
  const_iterator;
@)
  file_stack () :
  @[ std::stack<input_record,containers::mirrored_sl_list<input_record> >() @]
  @+{}
  const_iterator ctop() const @+{@; return c.wcbegin(); }
  const_iterator cbottom() const @+{@; return c.wcend(); }
};

@ When an input file is exhausted, the stored line number is restored, the
|stream| pointer deleted (which also closes the file) and the record popped
from the stack; then the |stream| is set to the previous input stream. When
|close_includes| is called all auxiliary input files are closed.

@h "global.h" // for |output_stream|

@< Definitions of class members @>=
void BufferedInput::pop_file()
{ line_no = input_stack.top().line_no;
  *output_stream << "Completely read file '" << cur_fname() << "'."
                 << std::endl;
  input_files_completed.insert (input_stack.top().name); // reading succeeded
  input_stack.pop(); // also closes the current file
  stream= input_stack.empty() ? &base_stream : &input_stack.top().f_stream;
}
@)
void BufferedInput::close_includes()
{ while (not input_stack.empty())
  { std::cerr << "Abandoning reading of file '" << cur_fname() @|
              << "' at line " << line_no << std::endl;
  @/line_no = input_stack.top().line_no;
    input_stack.pop();
  }
  stream= &base_stream;
}

@ Often we need the number or name of the topmost file on the |input_stack|.

@< Other methods of |BufferedInput| @>=
id_type current_file() const
{@; return input_stack.empty() ? Hash_table::empty : input_stack.top().name; }
const char* name_of(id_type f) const
{@; return
   f==Hash_table::empty ? "<standard input>" : input_files_seen.name_of(f); }
const char* cur_fname() const
  @+{@; return name_of(current_file()); }

@ When opening files, we shall look them up in a list of places, this search
path being specified by the user. For each of these places in turn, we shall
prefix the given file name with a string. How this string is stored and
recovered does not concern us here; it is defined elsewhere (in fact
in \.{main.w}). We shall just use the following functions that tell how many
path components there are, and give access to each component.

@< Declarations functions used but defined elsewhere @>=
unsigned int input_path_size();
const std::string& input_path_component(unsigned int i);

@ Pushing a new input file requires constructing the new |input_record| on the
stack, which will be done in place (the record is ``emplaced''). Opening the
file is done during member initialisation by the constructor of that record.
We also set |line_no| to the line number at which we shall resume reading the
interrupted file, since switching back to that file will happen when we pop
the record created here from the |input_stack|; this happens \emph{after} an
attempt to read from the associated file failed, at which point incrementing
the line number is already done, and won't be repeated.

If opening the file fails the constructor still succeeds; we leave it to the
calling function to detect this condition and destroy the created record.
However, in case the file is not open after initialising |stream|, presumably
because the file was not found, we do a second attempt after extending the
file name with |def_ext| (provided it was non-null, and not already a suffix
of |name|).

@< Definitions of class members @>=
BufferedInput::input_record::input_record
(BufferedInput& parent, const char* file_name)
@/: f_stream() // this tries to open the file
  , name(~0)
  , line_no(parent.line_no+parent.cur_lines) // record where reading will resume

{ unsigned int path_size = input_path_size();
  for (unsigned int i=0; i<=path_size; ++i)
  { std::string pathname(i<path_size ? input_path_component(i) : std::string());
    pathname += file_name;
    f_stream.open(pathname);
    if (not f_stream.is_open() and parent.def_ext!=nullptr)
    { long inx=pathname.size()-std::strlen(parent.def_ext);
      if (inx<0 or pathname.substr(inx)!=parent.def_ext)
        // only add extension if absent
        f_stream.open(pathname+= parent.def_ext);
          // try to reopen |stream| with extended name
    }
    if (f_stream.is_open())
       // once opening succeeds, record name in |input_files_seen|
    {@;name = parent.input_files_seen.match(pathname.c_str(),pathname.size());
      return;
    }
  }
}

@ If |push_file| is called with |skip_seen==true|, this means it should not
re-read a file that was already successfully read before. However, we want
this to apply to the actual file name, including the default extension
|def_ext| that might be added by the |input_record| constructor, and we don't
want to store file names that were unsuccessfully opened, so we postpone this
test until the |input_record| is constructed. To limit the drawback that this
strategy would lead us to open and then immediately close files in case they
are skipped, we add a short-cut that avoids any work in case |name| exactly
matches a file name previously read; this should be an incentive to write
an explicit extension when including a file from another file.

If the file just pushed on the stack was not successfully opened, we destroy
the stack record, print an error message, and report failure. Otherwise we may
still decide not to read from it (because this was previously done) and
immediately pop the record, but in such cases |push_file| reports success to
its caller nonetheless.

@< Definitions of class members @>=
bool BufferedInput::push_file(const char* name, bool skip_seen)
{
  if (skip_seen and input_files_seen.knows(name) @|
      and input_files_completed.isMember(input_files_seen.match_literal(name)))
    return true;
  input_stack.emplace(*this,name); // this tries to open the file
@)
  if (input_stack.top().f_stream.is_open())
  { @< In cases where reading from this file should be avoided,
       pop its record from the |input_stack|;
       otherwise prepare for reading from the new file @>
    return true; // succeed whether or not a file was actually pushed
  }
  else // opening file failed
  { input_stack.pop();
    // don't call |pop_file|: |stream| and |line_no| were not updated
    std::cerr << "failed to open input file '" << name << "'." << std::endl;
    return false;
  }
}

@ The conditions that lead to silently popping a just successfully opened file
are the presence of a file with the same name below it on the |input_stack| (a
test necessary to avoid input loops), or both the |skip_seen| argument is set,
and the file is marked as completely read in the |input_files_completed|
bitmap. Although the conditions are tested here in the opposite order, |skip|
is set to the logical ``or'' of both conditions in all cases.

When the file is really going to be read from, this is reported to
|*output_stream|, and the member variable |stream| made to point to
|input_stack.top().f_stream| so that the next reading operation will be from
the new file. We also take care to get the line numbering for the new file off
to a good start.

@< In cases where reading from this file should be avoided,... @>=
{ bool skip=false;
  const id_type file_nr = input_stack.top().name;
  if (file_nr==input_files_completed.capacity()) // it wasn't seen before
    input_files_completed.extend_capacity(false);
       // add new empty slot to bitmap; may set it later
  else // old name; need to do some checks
  { skip = skip_seen and input_files_completed.isMember(file_nr);
    auto it=input_stack.ctop();
    while (not skip and ++it!=input_stack.cbottom())
      if (file_nr==it->name)
        skip=true; // avoid recursive inclusion of active file
  }
  if (skip)
    input_stack.pop(); // silently skip including the file
  else
  { *output_stream << "Starting to read from file '" << cur_fname()
                   << "'." << std::endl;
  @/stream= &input_stack.top().f_stream;
    line_no=1; // prepare to read from pushed file
    cur_lines=0;
      // so we won't advance |line_no| when getting first line of new file
  }
}

@*1 Getting a line from the buffer.
%
The member function |getline| is usually called at times when |eol()| would
return |true| (otherwise some input will be discarded). This should only
happen when a new request for input arrives \emph{after} the end-of-line
condition has been transmitted at an earlier request (otherwise we would wait
for new input before executing a command), which is why we always return a
newline character from |shift()| before making |eol()| true. Thus our caller
does not have to distinguish between a fresh end-of-line condition and one
that has already been acted upon. For us it means that we should actually
store a newline character at the end of line buffer, in spite of the fact that
|std::getline| and |readline| strip that character. This also means that there
is no objection to calling our |getline| automatically in case |shift| is
called when |eol()| is true (any end-of-line action has been done at such a
point); we shall do so, and in fact this removes the burden of monitoring
|eol()| from our callers (one could even consider making that function
private, but its publicity does no real harm). We return success if either no
error flags were set on our |istream| parent, or else if at least some
characters were read (probably followed by an end-of-file condition in absence
of a final newline); in the latter case the next call will return false.

As we said in the introduction, we shall treat escaped newlines here, so that
the lexical analyser does not need to worry about them (and this will allow
escaping newlines in the middle of tokens, for instance string constants, to
be handled painlessly). A price to pay is that we cannot analyse the left
context of a possibly escaping character to see if it is not actually part of
a token. So we simply assert that a final backslash character on a line is
always one that escapes the newline; for tokens whose representation is or
ends with a backslash (in \.{axis} this is the case for the integer division
operator) the user should simply refrain from using them at the very end of a
line (one may add a comment after it to avoid the problem).

Since it is easy and efficient to do so, we also skip trailing spaces here;
this should make little difference since the lexical analyser will probably
ignore spaces anyway, except to separate tokens, and a non-escaped line end
will do so perfectly well without trailing spaces. However, we think it is
reassuring for users to know that invisible characters are eliminated at a low
level of input processing, so that they should really never make a difference.

@h <cctype>

@< Definitions of class members @>=
bool BufferedInput::getline()
{ line_buffer="";
  line_no+=cur_lines; cur_lines=0;
  std::string pr=@< Prompt for the next line@>@;@;;
  bool go_on, popped=false;
  do
  { std::string line;
  @/@< Get |line| without newline from |stream| if there is one; if none can be
       obtained then |break| if |cur_lines>0|, otherwise pop |input_stack|,
       set |popped=true| and |break|; if nothing works, |return false| @>
    ++cur_lines;
    std::string::size_type l=line.length();
    while (l>0 and std::isspace(line[l-1]))
      --l;
    if (l<line.length())
      line.erase(l); // remove trailing space
    pr= "\\ "; // continuation prompt
    if (go_on=(l>0 and line[l-1]=='\\'))
      // keep going only if final backslash present
      line.erase(l-1); // in which case it is chopped off
    line_buffer+=line; // append |line| to buffer in all cases
  }
  while(go_on);
  line_buffer.push_back(popped ? '\f' :'\n');
     // add suppressed newline, or end-of-file indication
  p=line_buffer.data();
@/return true;
  // delay reporting end of input if anything was read at all
}

@ When a temporary prompt is non-empty (the lexical analyser will typically
set it to indicate that for some reason the previous line is incomplete) then
the initial prompt (one not following an escaped newline) will be that
temporary prompt followed by a space and the secondary prompt.

@< Prompt for the next line@>=
(temp_prompt.empty() ? std::string(prompt==nullptr?"":prompt)
                     : temp_prompt+" "+prompt2)

@ Reading a line might fail, in which case an error condition will be set on
the input stream |*stream|. If the error condition is set, we don't even try
to get a line, but if it is cleared initially it could still be set after
trying to get a line; in the latter case we only act on the condition if
moreover we read nothing at all during the current call of this function.
In case we have read nothing this time, but |cur_lines>0| indicates that some
earlier input processed during this call ended with a backslash, we protest
against this abuse, but still return as if the final backslash had not been
present. If however |cur_lines==0| indicates that this is the first time
around, then failure to get anything from |stream| leads to popping of the
current input file, or failure of the current call if there is nothing on the
stack. Note that in case the main input stream ends immediately after a
command, then that command will be returned with a newline added to it the
first time (and presumably the command gets executed), and the following call
to |InputBuffer::getline| will fail, probably leading to program termination.

@< Get |line| without newline from |stream| ... @>=
{ if (stream->good())
    {@; @< Read a line into |line|, prompting with |pr| if appropriate @> }
  if (line.empty() and not stream->good()) // nothing found and file ended
  { if (cur_lines>0)
    { std::cerr << "end of file after backslash ";
      if (input_stack.empty())
        std::cerr << (stream==&std::cin ? "on standard input" : "in main file")
                  << std::endl;
      else
      { std::cerr << "in file '" << cur_fname() << '\'' << std::endl;
@/      pop_file(); popped=true;
      }
    }
    else if (input_stack.empty())
      return false; // all hope of returning something has gone
    else {@; pop_file(); popped=true; }
    break;
// in these cases terminate line-continuation; return the part that was found
  }
  else {} // something (at least a bona fide empty line) was contributed here
}

@ Actually getting a line is taken care of by the |readline| function if
present, or by |std::getline|. Both chop the newline from the line. A
particularity of the GNU |readline| function is that upon typing \.{\^D} at
the start of a line it will signal a file-ended condition by returning a null
pointer, but it will not actually set the end-of-file condition on the
standard input stream; to obtain uniform behaviour between the various input
streams we therefore manually set the end-of-file on |std::cin| when this
happens; we also echo \.{\^D} to give the user feedback that the end-of-file
signal was understood.

The GNU |readline| function is written in \Cee, so it cannot help returning
its result in a |malloc|-ed \Cee-style string. We copy it to the |string|
object |line| and then free the string returned by |readline|, unless
|add_hist| is also set and the line is non-empty: then the string is
presumable stored away without copying (the history library documentation
suggests this without stating it clearly), so calling |free| would be an
error.

@h <cstdlib>

@< Read a line... @>=
if (input_stack.empty() and prompt!=nullptr)
  // do only at top level, and only if prompt enabled
{ prompt_length= pr.length();
  if (readline!=nullptr) // skip calling 'readline' if no function is supplied
  { char* l=readline(pr.c_str());
    if (l==nullptr) // then |readline| failed, flag end of file
    {@; line=""; stream->setstate(std::ios_base::eofbit);
      std::cout << "^D\n";
    }
    else
    @/{@; line=l;
      if (add_hist!=nullptr and *l!='\0')
        add_hist(l);
      else
        std::free(l);
    }
  }
  else
  {@; std::cout << pr;
      std::getline(std::cin,line,'\n'); }
}
else
   std::getline(*stream,line,'\n');


@ We have seen that |shift| should call |getline| when necessary. It does so
at most once, since there is at least a newline character to return from the
line fetched. If |getline| is called and returns failure, then we return a
null character to signal end of input, while setting |p| to a null pointer to
ensure that the |eol()| condition holds, and that subsequent calls to |shift|
will continue to fail. Using a null character to signal file end makes this
buffer unsuited for reading binary files, but that is not our purpose here,
and it is more practical than returning a non-|char| value (which would
require a different signature). We want behaviour to be well defined when
|unshift| is called after |shift| returned a null character or in other
strange ways (before calling |shift| for instance), so we perform
simple test before decrementing |p|.

@< Inline fun... @>=
inline char BufferedInput::shift()
{ if (eol())
    if (not getline())
      {@; p=nullptr; return '\0'; } // signals file end
  return *p++;
}
@)
inline void BufferedInput::unshift()
{@; if (p!=nullptr and p>line_buffer.data())
      --p;
}


@*1 Secondary buffer methods.
%
The |BufferedInput| class provides some methods that are not directly related
to scanning tokens, but provide handles to manage information that is directly
related to the input process. First of all, we provide a method
|include_depth| to find out the number of currently open additional input
files. Then we provide the method |point| to obtain a pointer into the current
line buffer, which will be useful to isolate a token without reassembling it
character by character. We ensure that it will not return |nullptr|,
redirecting it |line_buffer.c_str()| in that case, which will point to a null
character. The method |set_line_no| can be used to set the internal line
counter. The user method |locate| works in the opposite direction as |point|:
it provides the line and column number of a pointer~|p| into |line_buffer|.

@< Other methods of |BufferedInput| @>=
unsigned int include_depth() const @+{@; return input_stack.size(); }
const char* point() const @+ {@;return p==nullptr ? line_buffer.c_str(): p;}
void set_line_no (unsigned long l) @+
{@; line_no=l; }
void locate (const char* p, int& line, int& column) const;

@ We deliver the current line and the offset of the pointer from
|line_buffer.data()|. This describes points that lie in lines continued by
escaped newlines somewhat cryptically (the line number is that of the first of
those lines, and the column number is the accumulated offset), but this could
be improved in the future by maintaining a list of indices at which escaped
newlines were removed.

@< Definitions of class members @>=
void BufferedInput::locate (const char* p, int& line, int& column) const
{@; line=line_no; column=p-line_buffer.data(); }

@ For error reporting, the parser may provide a range from one (line,column)
pair to another. The method |show_range| will do its best to translate this
into a comprehensible indication of the location.

@< Other methods of |BufferedInput| @>=
void show_range
 (std::ostream& out,unsigned long l0, int c0, unsigned long l1, int c1)
 const;


@ When we show a range, most of the time it will be entirely within
the current line. In that case we shall just print a line with spaces until
the start of the token, then carets for the duration of the token.
@< Definitions of class members @>=
void BufferedInput::show_range
  (std::ostream& out,unsigned long l0, int c0, unsigned long l1, int c1)
  const
{ if (l1==line_no)
  { int pl=prompt_length; // offset of last line on the screen
    if (stream!=&std::cin or cur_lines>1)
    { if (not input_stack.empty())
        out << "In input file '" << cur_fname() << "', ";
      if (stream!=&std::cin)
        out << "line " << l0 << ":\n";
      out<<line_buffer; pl=0; // echo line in these cases
    }
    if (l0<l1)
      c0=0; // forget start column of multi-line range
    for (int i=pl+c0; i>0; --i)
      out<<' ';
    for (int i=c1; i>c0; --i)
      out<<'^';
    out << std::endl;
    if (l0<l1)
      if (l1-l0==1)
        out << "Range started in previous line\n";
      else
        out << "Range started " << (l1-l0) << "lines above\n";
  }
  else // range ended before the current line
    out << "Range from line " << l0 << " column " << c0 @|
	<< " to line " << l1 << " column " << c1 << ".\n";
}


@ To change the prompt, the user may push characters into the temporary prompt
and pop them off.

@< Other methods of |BufferedInput| @>=
void push_prompt(char c);
char top_prompt() const;
void pop_prompt();
void reset();

@ When popping we silently ignore the case of an empty
prompt,since it is up to the parser to issue an error message, which it cannot
do if we already throw an exception here. One can clear the temporary prompt
by calling the |reset| method.

@h <stdexcept>
@< Definitions of class members @>=
void BufferedInput::push_prompt(char c) @+
{@; temp_prompt.push_back(c); }
@)
char BufferedInput::top_prompt() const
{@; size_t l=temp_prompt.length();
  return l==0 ? '\0' : temp_prompt[l-1];
}
@)
void BufferedInput::pop_prompt()
{@; size_t l=temp_prompt.length();
  if (l>0)
    temp_prompt.erase(l-1);
}
@)
void BufferedInput::reset() @+
{@; temp_prompt=""; }


@* Index.

% Local IspellDict: british

%%  LocalWords:  BufferedInput
