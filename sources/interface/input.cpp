/*
  This is input.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#ifndef NREADLINE

// use the readline version
#include "input_readline.c"

#else

// use the non-readline version
#include "input_noreadline.c"

#endif
