/*
  This is input_readline.c

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For copyright and license information see the LICENSE file
*/

/*
  This file is one of the two variants of input.cpp. It is called .c so
  that no object file is created by the makefile.
*/

#include "input.h"

#include <iostream>
#include <stdio.h>
#include <cstdlib>
// #include <cstdio> // didn't work for gcc 4.4

#include <readline/history.h>
#include <readline/readline.h>

#include "commands.h"

namespace atlas {

namespace {

const char* readLine(const char* prompt = "", bool toHistory = true);
char* completionGenerator(const char*, int);
void displayCompletions(char**, int, int);

}

/*****************************************************************************

        Chapter I -- The InputBuffer class

  An InputBuffer object can hold a line of input, and can be used as in
  |std::istringstream| object (since it is publicly derived from that class),
  in particular one can read from it using the 'source >> variable' syntax.
  In addition, it provides a |getline| method that will fill the string that
  it holds from standard input; this calls readline unless NREADLINE was set.

******************************************************************************/

namespace input {

/*!
  Forget string in input buffer, read a new line, and position at start
*/
void InputBuffer::getline(const char* prompt, bool toHistory)
{
  const char* line = readLine(prompt,toHistory); // non-owned pointer

  str(line==NULL ? std::cout << "qq\n","qq" : line); // 'qq' at end of input
  reset(); // clear flags and start reading at beginning
}

void InputBuffer::reset() // clear flags and start reading at beginning
{
  clear();
  seekg(0,std::ios_base::beg);
}


/*
  Synopsis: resets the stream to position pos.

  Clears the flags as well; the idea is to undo a peek-forward operation.
*/
void InputBuffer::reset(std::streampos pos)
{
  clear();
  seekg(pos);
}

} // |namespace input|

/*****************************************************************************

        Chapter II -- The HistoryBuffer class -- assumes NREADLINE undefined

  A HistoryBuffer is an input buffer with its own history. Each time a line is
  read, the history pointer is temporarily made to point to our own history
  record.

  The implementation below cannot be understood without looking into the
  sources of the history library. Rather than packing its static history
  variables into a static |HISTORY_STATE| structure, every call to
  |history_get_history_state| |malloc|s a fresh such structure, but pointers
  to such structures are never freed in the library (in particular not at the
  end of a |history_set_history_state| call). It it therefore our sorry duty
  as user to clean up the mess left by the history library. On th upside, we
  can call |history_set_history_state| with a pointer to our |state| record.

******************************************************************************/

namespace input {

/*
  Synopsis: the default constructor.

  Creates a clean slate history record, and empty input buffer.
*/
HistoryBuffer::HistoryBuffer()
  : InputBuffer(), state()
{
  state.entries=NULL;
  state.offset=0;
  state.length=0;
  state.size=0;
  state.flags=0;
}

/*
  Synopsis: constructor for a HistoryBuffer whose input buffer is initialized
  with str, while its history record is a clean slate.
*/

HistoryBuffer::HistoryBuffer(const std::string& str)
  : InputBuffer(str), state()
{
  state.entries=NULL;
  state.offset=0;
  state.length=0;
  state.size=0;
  state.flags=0;
}

/*
  Synopsis: destructor

  Frees memory occupied by the history record about to be forgotten.
*/
HistoryBuffer::~HistoryBuffer()

{
  HISTORY_STATE* global_history=history_get_history_state();
  history_set_history_state(&state); // substitue our history record

  clear_history(); // free all memory occupied by our history entries
  std::free(history_list()); // also free the array of history entries

  history_set_history_state(global_history); // restore history record
  std::free(global_history); // and free temporary that was used to hold it
}

// redefine the virtual |getline| method to use our own history record
void HistoryBuffer::getline(const char* prompt, bool toHistory)
{
  HISTORY_STATE* global_history=history_get_history_state();
  history_set_history_state(&state); // substitue our history record

  InputBuffer::getline(prompt,toHistory); // now call the base method

  HISTORY_STATE* our_history=history_get_history_state();
  state=*our_history; // copy history back
  std::free (our_history); // and free temporary that was used to hold it

  history_set_history_state(global_history); // restore history record
  std::free(global_history); // and free temporary that was used to hold it
}

} // namespace input

/*****************************************************************************

        Chapter III -- Functions declared in input.h
                    -- version _with_ readline

******************************************************************************/

namespace input {

bool hasQuestionMark(InputBuffer& buf)

/*
  Synopsis: looks if the next character assignment from buf would return '?'
*/

{
  std::streampos pos = buf.tellg();
  char x = 0;
  buf >> x;
  buf.reset(pos);

  return x == '?';
}

void initReadLine()

{
  using namespace commands;

  // allow setting of atlas-specific user preferences
  rl_readline_name = "Atlas";

  rl_completion_entry_function = completionGenerator;
  rl_completion_display_matches_hook = displayCompletions;
}

}

/*****************************************************************************

        Chapter III -- Local functions
                    -- version _with_ readline

******************************************************************************/

namespace {

char* completionGenerator(const char* text, int state)

/*
  Synopsis: function passed to the readline library for command completion.

  Precondition: char* is the word to be completed; state is a variable which
  tells if this is the first attempt to complete.

  If I understand correctly, it should return the next valid completion of
  text, 0 if there is none such. The passed char* should be written in a string
  correctly allocated by malloc, as it will be deallocated by readline!
*/

{
  using namespace commands;

  static std::vector<const char*> e;
  static std::vector<const char*>::iterator prev;

  const CommandTree* mode = currentMode();

  if (state == 0) { // compute list of completions
    mode->extensions(e,text);
    prev = e.begin();
  }

  if (prev == e.end())
    return 0;
  else {
    char* val = (char*)std::malloc(strlen(*prev)+1);
    std::strcpy(val,*prev);
    ++prev;
    return val;
  }
}

void displayCompletions(char** matches, int num, int)

/*
  Synopsis: function passed to the readline library for the display of
  completion lists
*/

{
  rl_crlf();
  fprintf(rl_outstream,"completions are: ");

  for (int j = 1; j <= num; ++j) {
    fprintf(rl_outstream,"%s",matches[j]);
    if (j < num)
      fprintf(rl_outstream,",");
  }

  rl_crlf();
  rl_reset_line_state ();

  return;
}


/*
  Synopsis: gets a line of input using the GNU readline library.

  NOTE: this code is more or less taken from the readline manual.
*/
const char* readLine (const char* prompt, bool toHistory)
{
  /* since |readline| allocates, and we shall |free| in a subsequent call, we
     need to hold the buffer pointer in a |static| variable
  */
  static char *line_read = NULL;

  if (line_read!=NULL) // deallocate first on every call except the first
    std::free(line_read); // no need to clear |line_read|, is overwritten next

  /* Get a line from the user. */
  line_read = readline(prompt); // this returns a freshly |malloc|ed buffer

  if (toHistory and line_read!=NULL and line_read[0]!='\0') // skip empty lines
    add_history(line_read); // add to the global (static) history record

  return line_read;
}

} // namespace

} // namespace atlas
