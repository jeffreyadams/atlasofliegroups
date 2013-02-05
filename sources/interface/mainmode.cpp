/*
  This is mainmode.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright 2013 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#include "mainmode.h"
#include "commands.h"
#include "test.h"     // to absorb test commands

#include "basic_io.h"
#include "complexredgp.h"
#include "complexredgp_io.h"
#include "realredgp.h"
#include "realredgp_io.h"
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
#include "tags.h"

namespace atlas {

namespace commands {

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

  void main_mode_entry() throw(EntryError);
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
  void showrealforms_f();
  void showdualforms_f();
  void blocksizes_f();
  void gradings_f();
  void strongreal_f();
  void dual_kgb_f();
  void help_f();

  // local variables
  // these have been changed to pointers to avoid swapping of G_C

  ComplexReductiveGroup* G_C_pointer=NULL;
  ComplexReductiveGroup* dual_G_C_pointer=NULL;
  complexredgp_io::Interface* G_I_pointer=NULL;

} // |namespace|


/*****************************************************************************

        Chapter I -- Functions declared in mainmode.h

******************************************************************************/

ComplexReductiveGroup& currentComplexGroup()

{
  return *G_C_pointer;
}

ComplexReductiveGroup& current_dual_group()
{
  if (dual_G_C_pointer==NULL)
    dual_G_C_pointer = new ComplexReductiveGroup
      (currentComplexGroup(), tags::DualTag());
  return *dual_G_C_pointer;
}

complexredgp_io::Interface& currentComplexInterface()
{
  return *G_I_pointer;
}

void replaceComplexGroup(ComplexReductiveGroup* G
			,complexredgp_io::Interface* I)
{
  delete G_C_pointer;
  delete dual_G_C_pointer;
  delete G_I_pointer;
  G_C_pointer=G;
  dual_G_C_pointer=NULL;
  G_I_pointer=I;
}


/****************************************************************************

        Chapter II -- The main mode |CommandNode|

  One instance of |CommandNode| for the main mode is created at the
  first call of |mainNode()|; further calls just return a reference to it.

*****************************************************************************/

namespace {

/*
  Synopsis: entry function to the main program mode.

  It attempts to set the group type interactively. Throws an EntryError on
  failure.
*/
void main_mode_entry() throw(EntryError)
{
  try {
    interactive::getInteractive(G_C_pointer,G_I_pointer);
  }
  catch(error::InputError& e) {
    e("complex group not set");
    throw EntryError();
  }
}

// function only called from |exitMode|
void main_mode_exit()
{
  replaceComplexGroup(NULL,NULL);
}

} // |namespace|


/*
  Synopsis: returns a |CommandNode| object that is constructed on first call.
*/
CommandNode mainNode()
{
  CommandNode result ("main: ",main_mode_entry,main_mode_exit);
  result.add("type",type_f);
  result.add("cmatrix",cmatrix_f);
  result.add("rootdatum",rootdatum_f);
  result.add("roots",roots_f);
  result.add("coroots",coroots_f);
  result.add("simpleroots",simpleroots_f);
  result.add("simplecoroots",simplecoroots_f);
  result.add("posroots",posroots_f);
  result.add("poscoroots",poscoroots_f);
  result.add("realform",realform_f);
  result.add("showrealforms",showrealforms_f);
  result.add("showdualforms",showdualforms_f);
  result.add("blocksizes",blocksizes_f);
  result.add("gradings",gradings_f);
  result.add("strongreal",strongreal_f);
  result.add("dualkgb",dual_kgb_f); // here, since no real form needed
  result.add("help",help_f); // override
  result.add("q",exitMode);

  // add test commands

  test::addTestCommands(result,MainmodeTag());
  return result;
}

/*****************************************************************************

        Chapter III --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/

namespace {

// Print the Cartan matrix on stdout.
void cmatrix_f()
{
  prettyprint::printMatrix
    (std::cout,currentComplexGroup().rootDatum().cartanMatrix());

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

  const RootDatum& rd = currentComplexGroup().rootDatum();

  WeightList::const_iterator first = rd.beginRoot();
  WeightList::const_iterator last = rd.endRoot();
  basic_io::seqPrint(file,first,last,"\n") << std::endl;
}

// Print the coroots in the lattice basis.
void coroots_f()
{
  ioutils::OutputFile file;

  const RootDatum& rd = currentComplexGroup().rootDatum();

  CoweightList::const_iterator first = rd.beginCoroot();
  CoweightList::const_iterator last = rd.endCoroot();
  basic_io::seqPrint(file,first,last,"\n") << std::endl;
}


// Print the simple roots in the lattice coordinates.
void simpleroots_f()
{
  const RootDatum& rd = currentComplexGroup().rootDatum();

  rootdata::WRootIterator first = rd.beginSimpleRoot();
  rootdata::WRootIterator last = rd.endSimpleRoot();
  basic_io::seqPrint(std::cout,first,last,"\n") << std::endl;
}

// Print the simple coroots in the lattice coordinates.
void simplecoroots_f()
{
  const RootDatum& rd = currentComplexGroup().rootDatum();

  rootdata::WRootIterator first = rd.beginSimpleCoroot();
  rootdata::WRootIterator last = rd.endSimpleCoroot();
  basic_io::seqPrint(std::cout,first,last,"\n") << std::endl;
}

// Print the positive roots in the lattice basis.
void posroots_f()
{
  ioutils::OutputFile file;

  const RootDatum& rd = currentComplexGroup().rootDatum();

  rootdata::WRootIterator first = rd.beginPosRoot();
  rootdata::WRootIterator last = rd.endPosRoot();
  basic_io::seqPrint(file,first,last,"\n") << std::endl;
}

// Print the positive coroots in the lattice basis.
void poscoroots_f()
{
  ioutils::OutputFile file;

  const RootDatum& rd = currentComplexGroup().rootDatum();

  rootdata::WRootIterator first = rd.beginPosCoroot();
  rootdata::WRootIterator last = rd.endPosCoroot();
  basic_io::seqPrint(file,first,last,"\n") << std::endl;
}

void help_f() // override more extensive help of empty mode by simple help
{
  help_mode.activate();
}

// Print the matrix of blocksizes.
void blocksizes_f()
{
  complexredgp_io::printBlockSizes(std::cout,currentComplexInterface());
}

// Activates real mode (user will select real form)
void realform_f()
{
  real_mode.activate();
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


// Print the gradings associated to the weak real forms.
void gradings_f()
{
  ComplexReductiveGroup& G_C = currentComplexGroup();

  // get Cartan class; abort if unvalid
  size_t cn=interactive::get_Cartan_class(G_C.Cartan_set(G_C.quasisplit()));

  ioutils::OutputFile file;

  static_cast<std::ostream&>(file) << std::endl;
  complexredgp_io::printGradings(file,cn,currentComplexInterface())
      << std::endl;

}

// Print information about strong real forms.
void strongreal_f()
{
  ComplexReductiveGroup& G_C = currentComplexGroup();

  // get Cartan class; abort if unvalid
  size_t cn=interactive::get_Cartan_class(G_C.Cartan_set(G_C.quasisplit()));

  ioutils::OutputFile file;
  file << "\n";
  realredgp_io::printStrongReal
    (file,
     currentComplexGroup(),
     currentComplexInterface().realFormInterface(),
     cn);
}

// Print a kgb table for a dual real form.
void dual_kgb_f()
{
  ComplexReductiveGroup& G_C = currentComplexGroup();

  const complexredgp_io::Interface& G_I = currentComplexInterface();
  const RealFormNbrList rfl = // get list of all dual real forms
    G_C.dualRealFormLabels(G_C.mostSplit(G_C.quasisplit()));

  RealFormNbr drf;

  interactive::getInteractive(drf,G_I,rfl,tags::DualTag());

  // the complex group must be in a variable: is non-const for real group
  ComplexReductiveGroup dG_C(G_C,tags::DualTag());
  RealReductiveGroup dG(dG_C,drf);

  std::cout << "dual kgbsize: " << dG.KGB_size() << std::endl;
  ioutils::OutputFile file;

  KGB kgb(dG,dG.Cartan_set());
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
    ComplexReductiveGroup* G;
    complexredgp_io::Interface* I;
    interactive::getInteractive(G,I);
    replaceComplexGroup(G,I);
    drop_to(main_mode); // drop invalidated descendant modes if called from them
  }
  catch(error::InputError& e) {
    e("complex group not changed");
  }

}

} // |namespace|

} // |namespace|

} // namespace atlas
