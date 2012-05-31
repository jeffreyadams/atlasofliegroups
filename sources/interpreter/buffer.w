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
different lifetimes. A |String_pool| will be implemented as a list of blocks
containing the individual strings, from newest to oldest. The blocks are
usually of the size specified when the pool is constructed, but if in rare
cases a string larger than the block size is requested, a single block of that
size will be created (but |block_size| will not change). To know how much
remains in the current block, we keep a separate record, in |chars_left|.


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
    const size_t block_size; // block size of the pool
    char* start; // start of block in current usage
    char* point; // first character available
    size_t chars_left; // number of characters left in current block
    const String_pool* prev;
      // start of list of previous blocks, for deallocation
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

@*1 Hash tables.
%
Our hash tables will use an instance of |String_pool| for storing its keys.
Their functionality is quite simple, and hides the actual hashing. Upon
construction the hash table is empty, and its stores every distinct string
that is passed to it, never removing one. It matches the string against all
strings it has seen before, and returns a sequence number: either the number
that was given to the same string when it was first encountered, or the next
unused number if the string was never seen before. Methods are also provided
for the inverse conversion, from a sequence number to the corresponding
string, and for finding the current number of entries stored. This will
allow working everywhere with sequence numbers to represent identifiers.

@h<vector>
@< Class declarations @>=

class Hash_table
{
public:
  typedef unsigned short id_type; // type of value representing identifiers
  static const id_type empty = ~0; // value reserved for empty slots
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
cause the program the program to crash). And we cannot help having pointers to
(possibly) signed characters in the first place, because that is for instance
what string denotations give us (and what many standard functions require).
The simplest solution would be to cast |s| to |(const unsigned char*)|, and
everything would be safe. We dare not do such a conversion for fear of losing
our job, so instead we convert every character accessed via~|s| laboriously
into an unsigned value, which can be done by a standard conversion; this
accounts for two of the three |static_cast|s below (the third is for
documentation only).

The code below is safe if either |l>0| (as will always be the case for
identifiers) or if |l==0| and |*s=='\0'| in which case it will return |0| as
hash code (this may occur when the user rather foolishly attempts
to open a file with an empty name).

@< Definitions of class members @>=

Hash_table::id_type Hash_table::hash
        (const char* s, size_t l) const // |s| a string of length |l>0|
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
Hash_table::id_type Hash_table::do_match
	(const char* str, size_t l, bool copy_string)
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
to know that |std::istream| is a class in the header file, which is done by
including \.{iosfwd}, and include \.{iostream} in the implementation. The
actual pointer to input streams will be contained in an |input_record|
structure that we shall discuss in more detail later.

@h <iostream>

@< Class declarations @>=
#include <iosfwd>
@< Define |struct input_record@;| @>

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
@)
  @< Data members of |BufferedInput| @>@;
@)
 public:
    BufferedInput (std::istream& s);
        // associate line buffer to raw input stream
    BufferedInput (const char* prompt
                  , rl_type rl=NULL, add_hist_type=NULL @|
		  , const char* prompt2="> "
                  , const char* def_ext=".rx"@|);
        // use |stdin|, maybe with readline
    ~BufferedInput();
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
BufferedInput* main_input_buffer=NULL;

@ In our implementation we use the class |std::string| and the template class
|std::stack|.

@< Includes needed in the header file @>=
#include <string>
#include <stack>

@~At construction a |BufferedInput| object will fix a reference |base_stream|
to the input stream to be used when no additional input files are open.
Furthermore it stores the last line read in its |line_buffer| field, and the
pointer~|p| points to the next character to be produced by |shift()|. In case
this buffer is used for terminal input, a prompt string is stored and will be
printed every time the program requests more input; this is the |prompt|
pointer so the lifetime of the value provided at construction should contain
that of the |BufferedInput| object created. We also store a secondary prompt
|prompt2| and provide a dynamic string |temp_prompt| that together can be used
to temporarily change the prompt to indicate that some input needs completion.
In case of terminal input, the input may be pre-processed by the function
|readline|, which is typically going to be the function of that name from the
GNU~\.{readline} library; then one may also want to use another function
stored in |add_hist| to store lines for interactive retrieval, typically the
|add_history| function from the history library. The member |line_no| holds
the number of the current line, or more precisely of the first of |cur_lines|
lines that are currently read into |line_buffer| (if more than one, they were
joined by escaped newlines). We record the length of the last prompt printed
in |prompt_length| for the purpose of pointing to tokens. Finally there is a
stack |input_stack| of active input files (which are accessed by pointers
since |std::ifstream| objects cannot be put in a stack themselves) together
with their file names, and a pointer |stream| that always points to the
current input stream.

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
std::stack<input_record> input_stack; // active input streams
Hash_table input_files_seen;
std::istream* stream; // points to the current input stream

@ There are two constructors, one for associating an input buffer to some
(raw) |istream| object (which may represent a disk file or pipe), another for
associating it to interactive input from |stdin|. A prompt and readline
function only apply to the second case (and |rl| may be a null pointer to
request no input editing). Constructing the class does not yet fetch a line.

@< Definitions of class members @>=
BufferedInput::BufferedInput (std::istream& s)@/
:base_stream(s)
,line_buffer()
,p(NULL)
,prompt("")
,prompt2("")
,def_ext(NULL)
,temp_prompt("")
,@|readline(NULL)
,add_hist(NULL)
,line_no(1)
,cur_lines(0)
,input_stack()
,input_files_seen()
,stream(&base_stream)
@+{}

BufferedInput::BufferedInput
(const char* pr, rl_type rl, add_hist_type ah,const char* pr2,const char*
de)@/
:base_stream(std::cin)
,line_buffer()
,p(NULL)
,prompt(pr)
,prompt2(pr2)
,def_ext(de)
,temp_prompt("")
,@|readline(rl)
,add_hist(ah)
,line_no(1)
,cur_lines(0)
,input_stack()
,stream(&base_stream)
@+{}


@*1 Input from auxiliary files.
%
For every currently open \emph{auxiliary} input stream (necessarily a file
stream), we store a pointer to the |ifstream| object, the file name (for error
reporting) and the line number of \emph{the stream that was interrupted} by
opening this auxiliary file (the line number in the file itself will be held
in the input buffer, as long as it is not interrupted). All this goes into an
|input_record| structure that has a constructor as unique method.

@< Define |struct input_record@;| @>=
struct input_record
{ std::string name; // the name of the file on which |stream| is opened
  std::ifstream* stream; // pointer owned by parent object
  unsigned long line_no; // this refers to the older input stream!
@)
  input_record(const char* file_name, const char* def_ext, unsigned long line);
};

@ It would have been convenient if each |input_record| owned its own
|std::ifstream| pointer, so that its destructor could take care of deleting
it, which would also close the associated file. This would however pose a
problem for the copy constructor (obligatory for objects in a |std::stack|),
since as it would have to duplicate the |std::ifstream| object pointed to, so
as to avoid double destruction of the same instance, while no copy constructor
for |std::ifstream| is accessible. The only viable solution would then be
using a reference-counted shared pointer, which seems an unnecessarily heavy
measure to take just in order to have ownership by individual |input_record|
values. So instead we have decided to give the ownership of the mentioned
pointers to the common |BufferedInput| object whose |input_stack| holds the
stack record.

Thus the destructor for |BufferedInput| takes care of deleting the |stream|
pointers before popping the containing record; this still closes the
associated file automatically. It is however exceptional that an instance of
|BufferedInput| would be destructed while it still has input files on its
stack; in normal operation the files will be closed when their record is
popped off the stack.

@h <fstream>

@< Definitions of class members @>=
BufferedInput::~BufferedInput()
{@; while (not input_stack.empty())
  {@; delete input_stack.top().stream;
    input_stack.pop();
  }
}

@ When an input file is exhausted, the stored line number is restored, the
|stream| pointer deleted (which also closes the file) and the record popped
from the stack; then the |stream| is set to the previous input stream. When
|close_includes| is called all auxiliary input files are closed. Although the
initial loop of that method coincides with the action of the destructor
function, it would seem inappropriate to call that function explicitly, since
the object continues to exist.

@< Definitions of class members @>=
void BufferedInput::pop_file()
{ line_no = input_stack.top().line_no;
  std::cout << "Completely read file '" << input_stack.top().name
            << "'." << std::endl;
  delete input_stack.top().stream;
  input_stack.pop();
  stream= input_stack.empty() ? &base_stream : input_stack.top().stream;
}
@)
void BufferedInput::close_includes()
{ while (not input_stack.empty())
  { std::cerr << "Abandoning reading of file '" << input_stack.top().name
              << "' at line " << line_no << std::endl;
    line_no = input_stack.top().line_no;
    delete input_stack.top().stream;
    input_stack.pop();
  }
  stream= &base_stream;
}

@ Pushing a new input file requires constructing the new |input_record| on the
stack. Opening the file will be done by the constructor of that record. If
opening the file fails the constructor still succeeds; we leave it to the
calling function to detect this condition and destroy the created record.

However, in case the file is not open after constructing the |fstream|,
presumably because the file was not found, we do a second attempt after
extending the file name with |def_ext| (provided it was non-null). Since
extending the string could throw an exception, we arrange for deleting the
|ifstream| object in this case through an explicit |catch| clause. Because
that data member |stream| follows |name|, the same problem cannot occur during
initialisation of data members: if the |stream| field gets initialised at all,
execution will certainly reach the body of the constructor, and so no memory
leak is possible here.

@< Definitions of class members @>=
input_record::input_record
(const char* file_name, const char* def_ext,unsigned long line) :
 name(file_name), stream (new std::ifstream(file_name)), line_no(line)
{ try
  {
    if (def_ext!=NULL and not stream->is_open())
    @/{@; name += def_ext;
      stream->open(name.c_str());
    }
  }
  catch (...)
  {@; delete stream; throw; }
}

@ If |push_file| is called with |skip_seen==true|, this means it should not
re-read a file that was already opened before. However, we want this to apply
to the actual file name, with addition of the default extension from |def_ext|
that might be added by the |input_record| constructor if necessary, and we
don't want to store file names that were unsuccessfully opened, so we postpone
this test until the |input_record| is constructed. To limit the drawback that
this strategy would lead us to open and then immediately close files in case
they are skipped, we add an optimisation clause that avoids any work in case
|name| exactly matches a file name previously opened; this should be an
incentive to include the extension when including a file from another file. It
can have a noticeable effect, but only in the weird case that a user first
successfully reads a file without extension, then deletes it and again tries
to load it without mentioning an extension; without the optimisation clause
our implementation would in that case try to load the file name with
extension, but in fact we won't do that.

The above constructor for |input_record| is then called with |line| equal to
the line number at which we shall resume reading the interrupted file, since
switching back to that file will happen after a future failing attempt to read
from the file opened here; at that point advancing the line number has already
taken place, and will not be repeated after popping the |input_stack|.

Only if the new file stream is in |good| state after opening do we set to
buffer |stream| member to it, and reset the line number to~1. In the contrary
case we print an error message, destroy the stack record, and report failure.
We eave it up to the caller to decide whether or not to abandon the currently
open input file(s): when including a second file from an included file fails,
it may not be very sensible to continue reading from the first included file,
since its commands probably depended on those of the second file. However such
considerations are not really appropriate at the buffer level, and should be
made by the caller.

@< Definitions of class members @>=
bool BufferedInput::push_file(const char* name, bool skip_seen)
{
  if (input_files_seen.knows(name))
    return true;
  input_stack.push(input_record(name,def_ext,line_no+cur_lines));
  // reading will resume there
  if (input_stack.top().stream->good())
  { const std::string& name=input_stack.top().name;
    size_t old_size = input_files_seen.nr_entries();
      // the value of |match| for a new name
    if (skip_seen and
        input_files_seen.match(name.c_str(),name.size())<old_size)
      // seen before
    @/{@; delete input_stack.top().stream;
      input_stack.pop();
    } // so just close file and pop record
    else
    { std::cout << "Starting to read from file '" << input_stack.top().name
                << "'." << std::endl;
      stream= input_stack.top().stream;
      line_no=1; // prepare to read from pushed file
      cur_lines=0;
        // so we won't advance |line_no| when getting first line of new file
    }
    return true; // succeed whether or not a file was actually pushed
  }
  else
  { std::cerr << "failed to open input file '" << input_stack.top().name
              << "'." << std::endl;
    delete input_stack.top().stream;
    input_stack.pop();
// no need to call |pop_file|: |stream|, |line_no| and |cur_lines| are unchanged
    return false;
  }
}

@ A line is ended if |p| is either null (indicating an absence of any line
buffered, as is initially the case), or points to a null character (the one
terminating the line buffer, which always follows a newline character).

@< Inline fun... @>=
inline bool BufferedInput::eol() const @+{@; return p==0 or *p=='\0'; }


@ The member function |getline| is usually called at times when |eol()| would
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
always one that escapes the newline; for tokens whose representation can ends
with a backslash (like the integer division operator in \.{realex}) the user
should simply refrain from using then at the end of a line (on can add a
comment after it if necessary).

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
  const char* pr=@< Prompt for the next line@>@;@;;
  bool go_on, popped=false;
  do
  { std::string line;
  @/@< Get |line| without newline from |stream| if there is one; if none can be
       obtained then |break| if |cur_lines>0|, otherwise pop |input_stack|,
       set |popped=true| and |break|; if nothing works return |false| @>
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
     // add suppressed newline, or form feed if file was popped
  p=line_buffer.data();
@/return true;
  // delay reporting end of input if anything was read at all
}

@ When a temporary prompt is non-empty (the lexical analyser will typically
set it to indicate that for some reason the previous line is incomplete) then
the initial prompt (one not following an escaped newline) will be that
temporary prompt followed by a space and the secondary prompt.

@< Prompt for the next line@>=
(temp_prompt.empty() ? prompt : (temp_prompt+" "+prompt2).c_str())

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
      { std::cerr << "in file '" << input_stack.top().name << '\'' << std::endl;
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
if (stream==&std::cin) // which implies |input_stack.empty()|
{ prompt_length= std::strlen(pr);
  if (readline!=NULL)
  { char* l=readline(pr);
    if (l==NULL) // then |readline| failed, flag end of file
    {@; line=""; stream->setstate(std::ios_base::eofbit);
      std::cout << "^D\n";
    }
    else
    @/{@; line=l;
      if (add_hist!=NULL and *l!='\0')
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
      {@; p=NULL; return '\0'; } // signals file end
  return *p++;
}
@)
inline void BufferedInput::unshift()
{@; if (p!=NULL and p>line_buffer.data())
      --p;
}


@* Secondary buffer methods.
%
The |BufferedInput| class provides some methods
that are not directly related to scanning tokens, but provide handles to
manage information that is directly related to the input process. First of
all, we provide a method |include_depth| to find out the number of currently
open additional input files. Then we provide the method |point| to obtain a
pointer into the current line buffer, which will be useful to isolate a token
without reassembling it character by character. The method |set_line_no| can
be used to set the internal line counter. The user method |locate| works in
the opposite direction as |point|: it provides the line and column number of a
pointer~|p| into |line_buffer|.

@< Other methods of |BufferedInput| @>=
unsigned int include_depth() const @+{@; return input_stack.size(); }
const char* point() const @+
{@;return p;} // points at next char in line buffer
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
        out << "In input file '" << input_stack.top().name << "', ";
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
char BufferedInput::top_prompt() const
{ size_t l=temp_prompt.length();
  return l==0 ? '\0' : temp_prompt[l-1];
}
void BufferedInput::pop_prompt()
{ size_t l=temp_prompt.length();
  if (l>0)
    temp_prompt.erase(l-1);
}
void BufferedInput::reset() @+
{@; temp_prompt=""; }


@* Index.

% Local IspellDict: british
