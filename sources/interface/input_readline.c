/*
  This is input_readline.c
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

/*
  This file is one of the two variants of input.cpp. It is called .c so
  that no object file is created by the makefile.
*/

#include "input.h"

#include <iostream>

namespace { // localize to avoid namespace pollution

#include <readline/history.h>
#include <readline/readline.h>

}

#include "commands.h"

namespace atlas {

namespace {

char* completionGenerator(const char*, int);
void displayCompletions(char**, int, int);
char* readLine(const char* prompt = "", bool toHistory = true);

}

/*****************************************************************************

        Chapter I -- The InputBuffer class

  ... explain here when it is stable ...

******************************************************************************/

namespace input {

std::istream& InputBuffer::getline(std::istream& is, const char* prompt,
				   bool toHistory)

/*
  Synopsis: reads 
*/

{
  std::string line = readLine(prompt,toHistory);

  str(line);
  reset();

  return is;
}

void InputBuffer::reset()

/*
  Synopsis: rewinds the stream for re-reading.
*/

{
  clear();
  seekg(0,std::ios_base::beg);

  return;
}

void InputBuffer::reset(std::streampos pos)

/*
  Synopsis: resets the stream to position pos.

  Clears the flags as well; the idea is to undo a peek-forward operation.
*/

{  
  clear();
  seekg(pos);

  return;
}

}
 	

/*****************************************************************************

        Chapter II -- The HistoryBuffer class -- assumes NREADLINE undefined

  A HistoryBuffer is an input buffer with its own history. On exit, the
  original history is restored.

******************************************************************************/

namespace input {

HistoryBuffer::HistoryBuffer()
  :InputBuffer()

/*
  Synopsis: the default constructor.

  Saves the current state of history, and starts a new one.
*/

{
  d_history = history_get_history_state();

  HISTORY_STATE* hs = (HISTORY_STATE*)malloc(sizeof(HISTORY_STATE));
  memset(hs,0,sizeof(HISTORY_STATE));

  history_set_history_state(hs);
}

HistoryBuffer::HistoryBuffer(const std::string& str)
  :InputBuffer(str)

/*
  Synopsis: constructor for a HistoryBuffer initialized with str.

  Saves the current state of history, and starts a new one.
*/

{
  d_history = history_get_history_state();

  HISTORY_STATE* hs = (HISTORY_STATE*)malloc(sizeof(HISTORY_STATE));
  memset(hs,0,sizeof(HISTORY_STATE));

  history_set_history_state(hs);
}

HistoryBuffer::~HistoryBuffer()

/*
  Synopsis: destructor

  Resets the history to its original state
*/

{  
  clear_history();

  HISTORY_STATE* hs = history_get_history_state();
  free(hs);

  history_set_history_state((HISTORY_STATE*)d_history);
}
  
}

/*****************************************************************************

        Chapter III -- Functions declared in input.h 
                    -- version _with_ readline

  ... explain here when it is stable ...

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

  return;
}

}

/*****************************************************************************

        Chapter III -- Local functions
                    -- version _with_ readline

  ... explain here when it is stable ...

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

  const CommandMode* mode = currentMode();

  if (state == 0) { // compute list of completions
    mode->extensions(e,text);
    prev = e.begin();
  }

  if (prev == e.end())
    return 0;
  else {
    char* val = (char*)malloc(strlen(*prev)+1);
    strcpy(val,*prev);
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

char * readLine (const char* prompt, bool toHistory)

/*
  Synopsis: gets a line of input using the GNU readline library.

  NOTE: this is lifted from the readline manual.

  NOTE: this is sure to use C-style i/o, and in particular stdin instead of
  cin. So we cannot forego synchronizing standard streams, unless and until
  a C++-version makes its appearance!
*/

{
  static char *line_read = 0;

  /* If the buffer has already been allocated,
     return the memory to the free pool. */
  if (line_read) {
      free (line_read);
      line_read = 0;
    }

  /* Get a line from the user. */
  line_read = readline (prompt);

  /* If the line has any text in it, and toHistory is true,
     save it on the history. */
  if (toHistory)
    if (line_read && *line_read)
      add_history (line_read);

  return (line_read);
}

}

}
