/*
  This is mainmode.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "mainmode.h"

#include "basic_io.h"
#include "complexredgp.h"
#include "complexredgp_io.h"
#include "emptymode.h"
#include "error.h"
#include "helpmode.h"
#include "interactive.h"
#include "io.h"
#include "ioutils.h"
#include "realform_io.h"
#include "realmode.h"
#include "rootdata.h"
#include "special.h"
#include "test.h"

namespace atlas {

namespace mainmode {

  using namespace commands;

}

/****************************************************************************

  This file contains the definitions for the "main" command mode, which
  is the central mode for the program.

  Basically, entering the main mode means setting the group type. In an
  interactive session, there is always a "current group"; this may be
  exported through the function "currentGroup". The current group may be
  changed with the "type" command, which amounts to exiting and re-entering
  the main mode.

  NOTE : it could be useful to have several active groups simultaneously, say
  for comparison purposes. This should be easy to do (an easy scheme would
  be for instance to have "first", "second" switch between two groups; but
  more sophisticated things are possible.)

*****************************************************************************/

namespace {

  using namespace mainmode;

  class ThisMode:public CommandMode {
  public:
  // constructors and destructors
    ThisMode();
    virtual ~ThisMode();
  // manipulators
    virtual const std::vector<const commands::CommandMode*>& next() const;
  };

  void this_entry() throw(EntryError);
  void this_exit();

  // functions for the predefined commands

  void blocksizes_f();
  void coroots_f();
  void help_f();
  void poscoroots_f();
  void posroots_f();
  void q_h();
  void roots_f();
  void realform_f();
  void showdualforms_f();
  void showrealforms_f();
  void simplecoroots_f();
  void simpleroots_f();
  void type_f();

  // local variables

  complexredgp::ComplexReductiveGroup G_C;
  complexredgp_io::Interface G_I(G_C);

 }

/****************************************************************************

        Chapter I -- The MainMode class

  Only one instance of this class will be constructed, on the first call
  of the function mainMode().

*****************************************************************************/

namespace {

ThisMode::ThisMode()
  :CommandMode("main: ",this_entry,this_exit)

{
  using namespace emptymode;

  // set parent mode
  d_prev = &emptyMode();

  // add the commands from the previous mode
  commands::addCommands(*this,prev());

  // add the commands from the current mode
  add("",relax_f);
  add("blocksizes",blocksizes_f);
  add("coroots",coroots_f);
  add("help",help_f);
  add("q",q_h);
  add("poscoroots",poscoroots_f);
  add("posroots",posroots_f);
  add("realform",realform_f);
  add("roots",roots_f);
  add("type",type_f);
  add("showdualforms",showdualforms_f);
  add("showrealforms",showrealforms_f);
  add("simplecoroots",simplecoroots_f);
  add("simpleroots",simpleroots_f);

  // add special commands

  special::addSpecialCommands(*this,MainmodeTag());

  // add test commands

  test::addTestCommands(*this,MainmodeTag());
}

ThisMode::~ThisMode()

{}

const std::vector<const commands::CommandMode*>& ThisMode::next() const

/*
  Synopsis: returns the list of direct descendants of this mode.

  The list is constructed on first call.

  Used for look-ahead in command completion.
*/

{
  using namespace commands;
  using namespace realmode;

  static std::vector<const CommandMode*> nextList;

  if (nextList.size() == 0) {
    const CommandMode* modePtr = &realMode();
    nextList.push_back(modePtr);
  }

  return nextList;
}

void this_entry() throw(EntryError)

/*
  Synopsis: entry function to the main program mode. 

  It attempts to set the group type interactively. Throws an EntryError on 
  failure.
*/

{
  using namespace interactive;
  using namespace realredgp;

  try {
    getInteractive(G_I);
  }
  catch(error::InputError& e) {
    e("complex group not set");
    throw EntryError();
  }

  return;
}

void this_exit()

/*
  Synopsis: resets G_C and G_I to their (meaningless) default values. This
  ensures the destruction of the previous ones.
*/

{
  using namespace complexredgp;
  using namespace complexredgp_io;

  ComplexReductiveGroup G_def;
  G_C.swap(G_def);
  Interface I_def(G_C);
  G_I.swap(I_def);

  return;
}

}


namespace mainmode {

const CommandMode& mainMode()

/*
  Synopsis: returns the ThisMode object. 

  It is constructed on first call.
*/

{
  static ThisMode mode;
  return mode;
}

}

/*****************************************************************************

        Chapter II --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/

namespace {

void blocksizes_f()

/*
  Synopsis: prints the matrix of blocksizes.
*/

{
  using namespace complexredgp;
  using namespace complexredgp_io;

  try {
    G_C.fillCartan();
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  printBlockSizes(std::cout,G_I);

  return;
}

void coroots_f()

/*
  Synopsis: prints the coroots in the lattice basis.
*/

{
  using namespace basic_io;
  using namespace ioutils;
  using namespace latticetypes;
  using namespace rootdata;

  ioutils::OutputFile file;

  const RootDatum& rd = G_C.rootDatum();

  WeightList::const_iterator first = rd.beginCoroot();
  WeightList::const_iterator last = rd.endCoroot();
  seqPrint(file,first,last,"\n") << std::endl;

  return;
}

void help_f()

{
  activate(helpmode::helpMode());
  return;
}

void poscoroots_f()

/*
  Synopsis: prints the positive coroots in the lattice basis.
*/

{
  using namespace basic_io;
  using namespace ioutils;
  using namespace rootdata;

  ioutils::OutputFile file;

  const RootDatum& rd = G_C.rootDatum();

  WRootIterator first = rd.beginPosCoroot();
  WRootIterator last = rd.endPosCoroot();
  seqPrint(file,first,last,"\n") << std::endl;

  return;
}

void posroots_f()

/*
  Synopsis: prints the positive roots in the lattice basis.
*/

{
  using namespace basic_io;
  using namespace ioutils;
  using namespace rootdata;

  ioutils::OutputFile file;

  const RootDatum& rd = G_C.rootDatum();

  WRootIterator first = rd.beginPosRoot();
  WRootIterator last = rd.endPosRoot();
  seqPrint(file,first,last,"\n") << std::endl;

  return;
}

void q_h()

{  
  io::printFile(std::cout,"q.help",io::MESSAGE_DIR);
  return;
}

void realform_f()

/*
  Synopsis: activates real mode.
*/

{
  using namespace commands;

  try {
    activate(realmode::realMode());
  }
  catch (EntryError) {
    return;
  }

  return;
}

void roots_f()

/*
  Prints the roots in the natural coordinates.
*/

{
  using namespace basic_io;
  using namespace ioutils;
  using namespace latticetypes;
  using namespace rootdata;

  ioutils::OutputFile file;

  const RootDatum& rd = G_C.rootDatum();

  WeightList::const_iterator first = rd.beginRoot();
  WeightList::const_iterator last = rd.endRoot();
  seqPrint(file,first,last,"\n") << std::endl;

  return;
}

void showdualforms_f()

{  
  using namespace realform_io;

  const realform_io::Interface rfi = G_I.dualRealFormInterface();

  std::cout << "(weak) dual real forms are:" << std::endl;
  printRealForms(std::cout,rfi);

  return;
}

void showrealforms_f()

{  
  using namespace realform_io;

  const realform_io::Interface rfi = G_I.realFormInterface();

  std::cout << "(weak) real forms are:" << std::endl;
  printRealForms(std::cout,rfi);

  return;
}

void simplecoroots_f()

/*
  Prints the simple coroots in the lattice coordinates.
*/

{
  using namespace basic_io;
  using namespace ioutils;
  using namespace rootdata;

  const RootDatum& rd = G_C.rootDatum();

  WRootIterator first = rd.beginSimpleCoroot();
  WRootIterator last = rd.endSimpleCoroot();
  seqPrint(std::cout,first,last,"\n") << std::endl;

  return;
}

void simpleroots_f()

/*
  Prints the simple roots in the lattice coordinates.
*/

{
  using namespace basic_io;
  using namespace ioutils;
  using namespace rootdata;

  const RootDatum& rd = G_C.rootDatum();

  WRootIterator first = rd.beginSimpleRoot();
  WRootIterator last = rd.endSimpleRoot();
  seqPrint(std::cout,first,last,"\n") << std::endl;

  return;
}

void type_f()

/*
  Resets the type, effectively restarting the program. If the construction
  of the new type fails, the current type remains in force.
*/

{
  using namespace interactive;
  using namespace realredgp;

  try {
    getInteractive(G_I);
  }
  catch(error::InputError& e) {
    e("complex group not changed");
    return;
  }

  return;
}

}

/*****************************************************************************

        Chapter III -- Functions declared in mainmode.h

  ... explain here when it is stable ...
  
******************************************************************************/

namespace mainmode {

complexredgp::ComplexReductiveGroup& currentComplexGroup()

{
  return G_C;
}

complexredgp_io::Interface& currentComplexInterface()

{
  return G_I;
}

}

}
