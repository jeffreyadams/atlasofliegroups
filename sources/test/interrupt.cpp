/* This file, producing an object file that gets linked into the Fokko program
   only, serves to provide a trivial target for the functions defined in the
   interpreter that are called from the Atlas library, when those library files
   are linked into the Fokko program which does not have the interpreter.
   Such interpreter calls are extremely rare (in fact absent for the most part
   of Atlas develment) but in order for long running library functions to be
   interruptible by the user (and the interrupt being caught in the interpreter
   loop) it is necessary to compile in a call to |check_interrupt()|.
*/

// The Fokko program just ignores calls of |interpreter::check_interrupt|
namespace atlas { namespace interpreter {
    void check_interrupt() {}
}}
