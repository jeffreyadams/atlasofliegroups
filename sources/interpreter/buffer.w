\def\emph#1{{\it #1\/}}

@* The buffered input class.
This file defines a simple class |BufferedInput|, which provides an interface
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
worry about that (and indeed it could not, if newlines are transparently
escaped).

Currently this is all the functionality implemented, but in the future one
could easily add methods to recursively open a subsidiary input file upon an
existing buffered input object, thus creating a stack of files from which
input will be read, possibly with transparent popping of exhausted files; if
one similarly enables temporarily redirecting input to a given string, one
could even implement a simple macro mechanism this way. But none of that is
done yet.


@( buffer.h @>=

#ifndef BUFFER_H
#define BUFFER_H
#include<iostream>
#include<string>
namespace atlas
{ namespace interpreter
   {@;
@< Class declaration @>@;
@< Inline function definitions @>@;
   }@;
}@;
#endif

@ Our main file will provide the class implementation
@h<iostream>
@h "buffer.h"
@c
using namespace std;
namespace atlas
{@; namespace interpreter
   {@;
@< Definitions of class members @>@;
   }@;
}@;


@ The class will be based on an |istream|, which seems the proper point in the
I/O class hierarchy, since we want to hide the distinction between input from
a file or from a terminal. However, we want to be able to build upon an
existing |istream|, notably |cin|, and since their copy constructor is
private, we must use a reference or pointer to the |istream|, which excludes
using heritage. This is no problem, since we do not want to allow calling such
|istream| member functions as |>>|, |getline|, |get| or |unget| directly
anyway, since this would bypass the line buffer (these calls would apply to
input beyond the current line, disregarding the unprocessed part of that
line).

@< Class declaration @>=
class BufferedInput
{ typedef char* (*rl_type)(const char* );
     // |rl_type| is type of pointer to readline function
  typedef void (*add_hist_type) (const char* );
     // |add_hist_type| is type of pointer to |add_history| function
  @< Data members of |BufferedInput| @>@;
  public:
    BufferedInput (std::istream& s);
        // associate line buffer to raw input stream
    BufferedInput (const char* prompt
                  , rl_type rl=0, add_hist_type=0
		  , const char* prompt2="> ");
        // use |stdin|, maybe with readline
@)
    char shift (); // inspect a new character
    void unshift (); // back up so last character will be reconsidered
    bool eol () const; // end of line: |true| if |shift| would return a newline
    bool getline ();
         // fetch a new line to replace current one; |true| if successful
    @< Other methods of |BufferedInput| @>@;
};


@ Our object stores the last line read in its |line_buffer| field, and the
pointer~|p| points to the next character to be produced by |shift()|. In case
this buffer is used for terminal input, a prompt string is stored and will be
printed every time the program requests more input; this is the |prompt|
pointer so the lifetime of the value provided at construction should contain
that of the |BufferedInput| object created. We also store a secondary prompt
|prompt2| and provide a dynamic string |temp_prompt| that together can be used
to temporarily change the prompt to indicate that some input needs completion.
In case of terminal input, the input may be pre-processed by the function
|readline|, which is typically going to be the function of that name from the
\.{readline} library; then one may also want to use another function stored in
|add_hist| to store lines for interactive retrieval, typically the
|add_history| function from the history library. The member |line_no| holds
the number of the current line, or more precisely of the first of |cur_lines|
lines that are currently read into |line_buffer| (if more than one, they were
joined by escaped newlines). We record the length of the last prompt printed
in |prompt_length| for the purpose of pointing to tokens.

@< Data members... @>=
std::istream& stream;
std::string line_buffer;
const char* p, *prompt, *prompt2;
std::string temp_prompt;
rl_type readline;
add_hist_type add_hist;
unsigned long line_no;
int cur_lines,prompt_length;

@ There are two constructors, one for associating an input buffer to some
(raw) |istream| object (which may represent a disk file or pipe), another for
associating it to interactive input from |stdin|. A prompt and readline
function only apply to the second case (and |rl| may be a null pointer to
request no input editing). Constructing the class does not yet fetch a line.

@< Definitions of class members @>=
BufferedInput::BufferedInput (istream& s)
@/:
stream(s),line_buffer(),p(NULL),prompt(""),prompt2(""),temp_prompt(""),@|
readline(NULL),add_hist(NULL),line_no(1),cur_lines(0)
{}

BufferedInput::BufferedInput
(const char* pr, rl_type rl, add_hist_type ah,const char* pr2)
@/:
stream(cin),line_buffer(),p(NULL),prompt(pr),prompt2(pr2),temp_prompt(""),@|
readline(rl),add_hist(ah),line_no(1),cur_lines(0)
{}

@ A line is ended if |p| is either null (indicating an absence of any line
buffered, as is initially the case), or points to a null character (the one
terminating the line buffer, which always follows a newline character).

@< Inline fun... @>=
inline bool BufferedInput::eol() const @+{@; return p==0 || *p=='\0'; }


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
private, but its publicity does no harm). We return success if either no
error flags were set on our |istream| parent, or else if at least some
characters were read (probably followed by an end-of-file condition in absence
of a final newline); in the latter case the next call will return false.

As we said in the introduction, we shall treat escaped newlines here, so that
the lexical analyser does not need to worry about them (and this will allow
escaping newlines in the middle of tokens, for instance string constants, to
be handled painlessly). A price to pay is that we cannot analyse the left
context of a possibly escaping character to see if it is not actually part of
a token. So we simply assert that a final backslash character on a line is
always one that escapes the newline; one should not use tokens whose
representation can end with a backslash (we hope this is not too restrictive).
Since it is easy and efficient to do so, we also skip trailing spaces here;
this should make little difference since the lexical analyser will probably
ignore spaces anyway, except to separate tokens (and a non-escaped line end
will do so perfectly well without trailing spaces), but we think it can be
useful to know that invisible characters are eliminated very early, so that
they should really never make a difference.

@h <cctype>

@< Definitions of class members @>=
bool BufferedInput::getline()
{ string line; line_buffer="";
  line_no+=cur_lines; cur_lines=0;
  const char* pr=@< Prompt for the next line@>@;@;;
  do
  { @< Read a line into |line| omitting the final newline,
       and increment |cur_lines|@>
    string::size_type l=line.length();
    while (l>0 && isspace(line[l-1])) --l;
    if (l<line.length()) line.erase(l); // remove trailing space
    if (l==0 || line[l-1]!='\\') {@; line_buffer+=line; break; }
    line.erase(l-1); line_buffer+=line; pr= "\\ "; // continuation prompt
  } while (true);
  line_buffer.push_back('\n'); // add newline that was suppressed on reading
  p=line_buffer.data();
  return stream.good() || line_buffer.length()>1;
}

@ When a temporary prompt is set (it is a non-empty string) then the prompt
for the next line will be that string followed by a space and the secondary
prompt.

@< Prompt for the next line@>=
(temp_prompt.empty() ? prompt : (temp_prompt+" "+prompt2).c_str())

@ Actually getting a line is taken care of by the |readline| function if
present, or by |std::getline|. Both chop the newline from the line, which is
the reason our module is called the way it is. The GNU |readline| function is
written in \Cee, so it cannot help returning its result in a |malloc|-ed
\Cee-style string. We copy it to the |string| object |line| and then free the
string returned by |readline|, unless |add_hist| is also set and the line is
non-empty: then the string is presumable stored away without copying (the
history library documentation suggests this without stating it clearly), so
calling |free| would be an error.

@h <cstring>
@< Read a line... @>=
{ prompt_length= &stream==&cin ? strlen(pr) : 0;
  if (readline!=0)
  { char* l=readline(pr);
    if (l==0) line="";
    else
    { line=l;
      if (add_hist!=0 && *l!='\0') add_hist(l);
      else free(l);
    }
  }
  else
  {@; if (&stream==&cin) cout << pr;
    std::getline(stream,line,'\n');
  }
  ++cur_lines;
}


@ We have seen that |shift| should call |getline| when necessary. It does so
at most once, since there is at least a newline character to return from the
line fetched. If |getline| is called and returns failure, then we return a
null character to signal end of input (this should probably not happen on
|stdin|). Using a null character to signal file end makes this buffer unsuited
for reading binary files, but that is not our purpose here, and it is more
practical than returning a non-|char| value (which would require a different
signature).

@< Inline fun... @>=
inline char BufferedInput::shift()
{ if (eol())
    if (!getline()) return '\0'; // signals file end
  return *p++;
}

inline void BufferedInput::unshift() {@; --p; }


@* Secondary buffer methods.
The |BufferedInput| class provides some methods that are not directly related
to scanning tokens, but provide handles to manage information that is directly
related to the input process. First of all, we provide a method |point| to
obtain a pointer into the current line buffer, which will be useful to isolate
a token without reassembling it character by character. The method
|set_line_no| can be used to set the internal line counter. The user method
|locate| works in the opposite direction as |point|: it  provides the line and
column number of a pointer~|p| into |line_buffer|.

@< Other methods of |BufferedInput| @>=
const char* point() const @+{@;return p;}
   // points at next char in line buffer
void set_line_no (unsigned long l) @+{@; line_no=l; }
void locate (const char* p, int& line, int& column) const;

@ We deliver the current line and the offset of the pointer from
|line_buffer.data()|. This describes points that lie in lines continued by
escaped newlines somewhat cryptically (the line number is that of the first of
those lines, and the column number is the accumulated offset), but this could
be improved later by maintaining a list of indices at which escaped newlines
were removed.

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


@ When we have to show a range, most of the time it will be entirely within
the current line. In that case we shall just print a line with spaces until
the start of the token, then carets for the duration of the token.
@< Definitions of class members @>=
void BufferedInput::show_range
  (std::ostream& out,unsigned long l0, int c0, unsigned long l1, int c1)
  const
{ if (l1==line_no)
  { int pl=prompt_length; // offset of last line on the screen
    if (&stream!=&cin || cur_lines>1)
      {@; out<<line_buffer; pl=0;} // echo line in these cases
    if (l0<l1) c0=0; // forget start column of multi-line range
    for (int i=pl+c0; i>0; --i) out<<' ';
    for (int i=c1; i>c0; --i) out<<'^';
    out << endl;
    if (l0<l1)
      if (l1-l0==1) out << "Range started in previous line\n";
      else out << "Range started " << (l1-l0) << "lines above\n";
  }
  else // range ended before the current line
  { out << "Range from line " << l0 << " column " << c0 @|
	<< " to line " << l1 << " column " << c1 << ".\n";
  }
}


@ To change the prompt, the user may push characters into the temporary prompt
and pop them off.

@< Other methods of |BufferedInput| @>=
void push_prompt(char c);
void pop_prompt();
void reset();

@ When popping we silently ignore the case of an empty
prompt,since it is up to the parser to issue an error message, which it cannot
do if we already throw an exception here. One can clear the temporary prompt
by calling the |reset| method.

@h <stdexcept>
@< Definitions of class members @>=
void BufferedInput::push_prompt(char c) @+{@; temp_prompt.push_back(c); }
void BufferedInput::pop_prompt()
{ size_t l=temp_prompt.length();
  if (l>0) temp_prompt.erase(l-1);
}
void BufferedInput::reset()  @+{@; temp_prompt=""; }


@* Index.

% Local IspellDict: british