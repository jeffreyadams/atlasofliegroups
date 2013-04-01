/*
  This is main.cpp

  Copyright (c) 2004 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations project

  For license information see the LICENSE file

*/

#include <cstdio>   // not obviously used, but appears helpful for Windows
#include <iostream>
#include <stdexcept>
#include <cstdlib> // for |exit|
#include <cstring> // for |strcmp|

#include "commands.h"
#include "emptymode.h"
#include "error.h"
#include "input.h"
#include "version.h"
#include "io.h" // so modify message directory

// This one function should not reside in the |atlas| (or any other) namespace
int main(int argc, char* argv[])
{
  ++argv, --argc; //we don't need to know the name we were called by; skip it
  if (argc>0 and std::strcmp(argv[argc-1],"--version")==0)
  {
    std::cout << atlas::version::NAME << ' ' << atlas::version::VERSION
	      << std::endl;
    std::exit(0);
  }
  if (argc>0) // then use first argument to override message directory
    atlas::io::MESSAGE_DIR=*argv++, --argc; // use and skip this argument

  try

  {

    atlas::input::initReadLine();
    atlas::commands::empty_mode.run();

    std::exit(0);
  }

  catch (atlas::error::NumericOverflow& e) {

    std::cerr << "error: uncaught NumericOverflow" << std::endl;

  }

  catch (atlas::error::NumericUnderflow& e) {

    std::cerr << "error: uncaught NumericUnderflow" << std::endl;

  }

  catch (std::exception& e) {

    std::cerr << "error: uncaught standard exceptions on exit"
	      << std::endl << e.what()
	      << std::endl;
  }


  catch (...) {

    std::cerr << std::endl << "error: uncaught exceptions on exit"
	      << std::endl;

  }

  std::exit(1);
}
