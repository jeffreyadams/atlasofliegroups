/*
  This is emptymode.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#include "emptymode.h"

#include <iostream>

#include "helpmode.h"
#include "io.h"
#include "interactive.h"
#include "wgraph.h"
#include "wgraph_io.h"

#include "mainmode.h"
#include "test.h"
#include "version.h"

/****************************************************************************

  This file contains the definitions for the "empty" command mode, which
  is the startup mode for the program. It defines a number of commands that
  will be inherited by all modes; notably the "", "help" and "qq" commands.

*****************************************************************************/

namespace atlas {

namespace {

  using namespace emptymode;

  void printVersion();

  // functions for the predefined commands

  void help_f();
  void q_h();
  void type_f();
  void extract_graph_f();
  void extract_cells_f();

}

/****************************************************************************

        Chapter I -- The empty mode |CommandMode|

  One instance of |CommandMode| for the empty mode is created at the
  first call of |emptyMode()|; further calls just return a reference to it.

*****************************************************************************/

namespace emptymode {


/*
  Synopsis: returns a |CommandMode| object that is constructed on first call.
*/
commands::CommandMode& emptyMode()
{
  static commands::CommandMode
    empty_mode("empty: ",printVersion,commands::relax_f);

  if (empty_mode.empty()) // amazingly only true on first call!
  {
    empty_mode.add("",commands::relax_f); // so typing nothing is OK
    empty_mode.add("help",help_f);
    empty_mode.add("q",q_h); // in empty mode "q" just gives help information
    empty_mode.add("qq",commands::exitInteractive);
    // the type command needs to be recognized in the empty mode, or else
    // it will trigger activation of the main mode and _then_ execute, which
    // leads to setting the type twice!
    empty_mode.add("type",type_f);
    empty_mode.add("extract-graph",extract_graph_f);
    empty_mode.add("extract-cells",extract_cells_f);

    test::addTestCommands(empty_mode,EmptymodeTag());
  }
  return empty_mode;
}

} // namespace emptymode

/*****************************************************************************

        Chapter II --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/

namespace {

void help_f()

{
  helpmode::intro_h();
  activate(helpmode::helpMode());
}

void q_h()

{
  io::printFile(std::cerr,"q.help",io::MESSAGE_DIR);
}

void type_f()

{
  try {
    activate(mainmode::mainMode());
  }
  catch (commands::EntryError) { // silently ignore failure to enter
  }
}


void extract_graph_f()
{
  ioutils::InputFile block_file("block information");
  ioutils::InputFile matrix_file("matrix information");
  ioutils::InputFile polynomial_file("polynomial information");
  ioutils::OutputFile file;

  wgraph::WGraph wg=wgraph::wGraph(block_file,matrix_file,polynomial_file);
  wgraph_io::printWGraph(file,wg);
}

void extract_cells_f()
{
  ioutils::InputFile block_file("block information");
  ioutils::InputFile matrix_file("matrix information");
  ioutils::InputFile polynomial_file("polynomial information");
  ioutils::OutputFile file;

  wgraph::WGraph wg=wgraph::wGraph(block_file,matrix_file,polynomial_file);
  wgraph::DecomposedWGraph dg(wg);
  wgraph_io::printWDecomposition(file,dg);
}


}

/****************************************************************************

        Chapter III --- Utilities.

  This section contains some utility functions used in this module :

    - printVersion() : prints version info on startup;

*****************************************************************************/

namespace {

void printVersion()

/*
  Prints an opening message and the version number.
*/

{
  std::cout << "This is " << version::NAME << " version " << version::VERSION
	    << "." << std::endl;
  std::cout << "Build date: " << version::COMPILEDATE << "." << std::endl;
  std::cout <<
    "Enter \"help\" if you need assistance."
	    << std::endl << std::endl;
}

} // namespace

} // namespace atlas
