/*
  This is mainmode.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For copyright and license information see the LICENSE file
*/

#include "mainmode.h"

#include "basic_io.h"
#include "complexredgp.h"
#include "complexredgp_io.h"
#include "realredgp.h"
#include "emptymode.h"
#include "error.h"
#include "helpmode.h"
#include "interactive.h"
#include "io.h"
#include "ioutils.h"
#include "realform_io.h"
#include "realmode.h"
#include "rootdata.h"
#include "kgb.h"
#include "kgb_io.h"
#include "testprint.h"
#include "prettyprint.h"
#include "special.h"
#include "test.h"

namespace atlas {


/****************************************************************************

  This file contains the definitions for the "main" command mode, which
  is the central mode for the program.

  Basically, entering the main mode means setting the group type. In an
  interactive session, there is always a "current group"; this may be exported
  through the function "currentComplexGroup". The current group may be changed
  with the "type" command, which amounts to exiting and re-entering the main
  mode.

  NOTE : it could be useful to have several active groups simultaneously, say
  for comparison purposes. This should be easy to do (an easy scheme would
  be for instance to have "first", "second" switch between two groups; but
  more sophisticated things are possible.)

*****************************************************************************/

namespace {

  using namespace mainmode;

  void main_mode_entry() throw(commands::EntryError);
  void main_mode_exit();

  // functions for the predefined commands

  void type_f();
  void cmatrix_f();
  void rootdatum_f();
  void roots_f();
  void coroots_f();
  void simpleroots_f();
  void simplecoroots_f();
  void posroots_f();
  void poscoroots_f();
  void realform_f();
  void showdualforms_f();
  void showrealforms_f();
  void blocksizes_f();
  void dual_kgb_f();
  void help_f();

  // local variables
  // these have been changed to pointers to avoid swapping of G_C

  complexredgp::ComplexReductiveGroup* G_C_pointer=NULL;
  complexredgp_io::Interface* G_I_pointer=NULL;

 }


/*****************************************************************************

        Chapter I -- Functions declared in mainmode.h

******************************************************************************/

namespace mainmode {

complexredgp::ComplexReductiveGroup& currentComplexGroup()

{
  return *G_C_pointer;
}

complexredgp_io::Interface& currentComplexInterface()
{
  return *G_I_pointer;
}

void replaceComplexGroup(complexredgp::ComplexReductiveGroup* G
			,complexredgp_io::Interface* I)
{
  delete G_C_pointer;
  delete G_I_pointer;
  G_C_pointer=G;
  G_I_pointer=I;

}

} // namespace mainmode

/****************************************************************************

        Chapter II -- The main mode |CommandMode|

  One instance of |CommandMode| for the main mode is created at the
  first call of |mainMode()|; further calls just return a reference to it.

*****************************************************************************/

namespace {

/*
  Synopsis: entry function to the main program mode.

  It attempts to set the group type interactively. Throws an EntryError on
  failure.
*/
void main_mode_entry() throw(commands::EntryError)
{
  try {
    interactive::getInteractive(G_C_pointer,G_I_pointer);
  }
  catch(error::InputError& e) {
    e("complex group not set");
    throw commands::EntryError();
  }
}

// function only called from |commands::exitMode|
void main_mode_exit()
{
  replaceComplexGroup(NULL,NULL);
}

} // namespace


namespace mainmode {

/*
  Synopsis: returns a |CommandMode| object that is constructed on first call.
*/
commands::CommandMode& mainMode()
{
  static commands::CommandMode main_mode
    ("main: ",main_mode_entry,main_mode_exit);
  if (main_mode.empty()) // true on first call
  {
    // add the commands from the empty mode
    commands::addCommands(main_mode,emptymode::emptyMode());

    // add the commands from the current mode
    main_mode.add("type",type_f);
    main_mode.add("cmatrix",cmatrix_f);
    main_mode.add("rootdatum",rootdatum_f);
    main_mode.add("roots",roots_f);
    main_mode.add("coroots",coroots_f);
    main_mode.add("simpleroots",simpleroots_f);
    main_mode.add("simplecoroots",simplecoroots_f);
    main_mode.add("posroots",posroots_f);
    main_mode.add("poscoroots",poscoroots_f);
    main_mode.add("realform",realform_f);
    main_mode.add("showdualforms",showdualforms_f);
    main_mode.add("showrealforms",showrealforms_f);
    main_mode.add("blocksizes",blocksizes_f);
    main_mode.add("dualkgb",dual_kgb_f); // here, since no real form needed
    main_mode.add("help",help_f); // override
    main_mode.add("q",commands::exitMode);

    // add special commands

    special::addSpecialCommands(main_mode,MainmodeTag());

    // add test commands

    test::addTestCommands(main_mode,MainmodeTag());
  }
  return main_mode;
}

} // namespace mainmode

/*****************************************************************************

        Chapter III --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/

namespace {

// Print the Cartan matrix on stdout.
void cmatrix_f()
{
  latticetypes::LatticeMatrix q;

  rootdata::cartanMatrix(q,currentComplexGroup().rootDatum());
  prettyprint::printMatrix(std::cout,q);

}


// Print information about the root datum (see testprint.cpp for details).
void rootdatum_f()
{
  ioutils::OutputFile file;
  testprint::print(file,currentComplexGroup().rootDatum());
}


// Print the roots in the lattice basis.
void roots_f()
{
  ioutils::OutputFile file;

  const rootdata::RootDatum& rd = currentComplexGroup().rootDatum();

  latticetypes::WeightList::const_iterator first = rd.beginRoot();
  latticetypes::WeightList::const_iterator last = rd.endRoot();
  basic_io::seqPrint(file,first,last,"\n") << std::endl;
}

// Print the coroots in the lattice basis.
void coroots_f()
{
  ioutils::OutputFile file;

  const rootdata::RootDatum& rd = currentComplexGroup().rootDatum();

  latticetypes::WeightList::const_iterator first = rd.beginCoroot();
  latticetypes::WeightList::const_iterator last = rd.endCoroot();
  basic_io::seqPrint(file,first,last,"\n") << std::endl;
}


// Print the simple roots in the lattice coordinates.
void simpleroots_f()
{
  const rootdata::RootDatum& rd = currentComplexGroup().rootDatum();

  rootdata::WRootIterator first = rd.beginSimpleRoot();
  rootdata::WRootIterator last = rd.endSimpleRoot();
  basic_io::seqPrint(std::cout,first,last,"\n") << std::endl;
}

// Print the simple coroots in the lattice coordinates.
void simplecoroots_f()
{
  const rootdata::RootDatum& rd = currentComplexGroup().rootDatum();

  rootdata::WRootIterator first = rd.beginSimpleCoroot();
  rootdata::WRootIterator last = rd.endSimpleCoroot();
  basic_io::seqPrint(std::cout,first,last,"\n") << std::endl;
}

// Print the positive roots in the lattice basis.
void posroots_f()
{
  ioutils::OutputFile file;

  const rootdata::RootDatum& rd = currentComplexGroup().rootDatum();

  rootdata::WRootIterator first = rd.beginPosRoot();
  rootdata::WRootIterator last = rd.endPosRoot();
  basic_io::seqPrint(file,first,last,"\n") << std::endl;
}

// Print the positive coroots in the lattice basis.
void poscoroots_f()
{
  ioutils::OutputFile file;

  const rootdata::RootDatum& rd = currentComplexGroup().rootDatum();

  rootdata::WRootIterator first = rd.beginPosCoroot();
  rootdata::WRootIterator last = rd.endPosCoroot();
  basic_io::seqPrint(file,first,last,"\n") << std::endl;
}

void help_f() // override more extensive help of empty mode by simple help
{
  activate(helpmode::helpMode());
}

// Print the matrix of blocksizes.
void blocksizes_f()
{
  currentComplexGroup().fillCartan();
  complexredgp_io::printBlockSizes(std::cout,currentComplexInterface());
}

// Activates real mode (user will select real form)
void realform_f()
{
  commands::activate(realmode::realMode());
}


void showrealforms_f()
{
  const realform_io::Interface& rfi =
    currentComplexInterface().realFormInterface();

  std::cout << "(weak) real forms are:" << std::endl;
  realform_io::printRealForms(std::cout,rfi);
}

void showdualforms_f()
{
  const realform_io::Interface& rfi =
    currentComplexInterface().dualRealFormInterface();

  std::cout << "(weak) dual real forms are:" << std::endl;
  realform_io::printRealForms(std::cout,rfi);
}


// Print a kgb table for a dual real form.
void dual_kgb_f()
{
  complexredgp::ComplexReductiveGroup& G_C = currentComplexGroup();
  G_C.fillCartan(); // must generate all Cartans: no real form chosen

  const complexredgp_io::Interface& G_I = currentComplexInterface();
  const realform::RealFormList rfl = // get list of all dual real forms
    G_C.dualRealFormLabels(G_C.mostSplit(G_C.quasisplit()));

  realform::RealForm drf;

  interactive::getInteractive(drf,G_I,rfl,tags::DualTag());

  // the complex group must be in a variable: is non-const for real group
  complexredgp::ComplexReductiveGroup dG_C(G_C,tags::DualTag());
  realredgp::RealReductiveGroup dG(dG_C,drf);
  dG.fillCartan();

  std::cout << "dual kgbsize: " << dG.kgbSize() << std::endl;
  ioutils::OutputFile file;

  kgb::KGB kgb(dG);
  kgb_io::printKGB(file,kgb);
}

/*
  Reset the type, effectively reentering the main mode. If the construction
  of the new type fails, the current type remains in force.
*/
void type_f()
{
  try
  {
    complexredgp::ComplexReductiveGroup* G;
    complexredgp_io::Interface* I;
    interactive::getInteractive(G,I);
    replaceComplexGroup(G,I);
  }
  catch(error::InputError& e) {
    e("complex group not changed");
  }

}

} // namespace

} // namespace atlas
