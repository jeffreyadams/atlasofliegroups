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

#include "lietype.h"
#include "basic_io.h"
#include "innerclass.h"
#include "output.h"
#include "realredgp.h"
#include "emptymode.h"
#include "error.h"
#include "helpmode.h"
#include "interactive.h"
#include "io.h"
#include "ioutils.h"
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
  through the function "current_inner_class". The current group may be changed
  with the "type" command, which amounts to exiting and re-entering the main
  mode.

  NOTE (by Fokko): it could be useful to have several active groups
  simultaneously, say for comparison purposes. This should be easy to do (an
  easy scheme would be for instance to have "first", "second" switch between
  two groups; but more sophisticated things are possible.)

*****************************************************************************/

namespace {

  void main_mode_entry();
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
  void dualkgb_f();
  void help_f();

  // local variables
  // most of these are pointers to avoid swapping of inner classes

  lietype::Layout layout;
  LatticeMatrix lattice_basis; // allows mapping Lie type info to complex group
  InnerClass* ic_pointer=nullptr;
  InnerClass* dual_ic_pointer=nullptr;
  output::Interface* G_I_pointer=nullptr;

} // |namespace|


/*****************************************************************************

        Chapter I -- Functions declared in mainmode.h

******************************************************************************/

InnerClass& current_inner_class()

{
  return *ic_pointer;
}

InnerClass& current_dual_inner_class()
{
  if (dual_ic_pointer==nullptr)
    dual_ic_pointer = new InnerClass(current_inner_class(),tags::DualTag());
  return *dual_ic_pointer;
}


const lietype::Layout& current_layout() { return layout; }
const LatticeMatrix& current_lattice_basis() { return lattice_basis; }

output::Interface& currentComplexInterface()
{
  return *G_I_pointer;
}

void replace_inner_class(InnerClass* G,output::Interface* I)
{
  delete ic_pointer;
  delete dual_ic_pointer;
  delete G_I_pointer;
  ic_pointer=G;
  dual_ic_pointer=nullptr;
  G_I_pointer=I;
}


/****************************************************************************

        Chapter II -- The main mode |CommandNode|

*****************************************************************************/

namespace {

/*
  Synopsis: entry function to the main program mode.

  It attempts to set the group type interactively. Throws an EntryError on
  failure.
*/
void main_mode_entry()
{
  try {
    interactive::get_group_type(ic_pointer,G_I_pointer,layout,lattice_basis);
    dual_ic_pointer=nullptr; // this one is lazily initialised
  }
  catch(error::InputError& e) {
    e("complex group not set");
    throw EntryError();
  }
}

// function only called from |exitMode|
void main_mode_exit()
{
  replace_inner_class(nullptr,nullptr); lattice_basis.clear();
}

} // |namespace|


/*
  Synopsis: returns a |CommandNode| object that is constructed during the call.
*/
CommandNode mainNode()
{
  CommandNode result ("main: ",main_mode_entry,main_mode_exit);
  result.add("type",type_f,"override"); // tag is a dummy
  result.add("cmatrix",cmatrix_f,"prints the Cartan matrix",std_help);
  result.add("rootdatum",rootdatum_f,"outputs the root datum",std_help);
  result.add("roots",roots_f,
	     "outputs the roots in the lattice basis",std_help);
  result.add("coroots",coroots_f,
	     "outputs the coroots in the lattice basis",std_help);
  result.add("simpleroots",simpleroots_f,
	     "outputs the simple roots in the lattice basis",std_help);
  result.add("simplecoroots",simplecoroots_f,
	     "outputs the simple coroots in the lattice basis",std_help);
  result.add("posroots",posroots_f,
	     "outputs the positive roots in the lattice basis",std_help);
  result.add("poscoroots",poscoroots_f,
	     "outputs the positive coroots in the lattice basis",std_help);
  result.add("realform",realform_f,
	     "sets the real form for the group",std_help);
  result.add("showrealforms",showrealforms_f,
	     "outputs the weak real forms for this complex group",std_help);
  result.add("showdualforms",showdualforms_f,
	     "outputs the weak real forms for the dual group",std_help);
  result.add("blocksizes",blocksizes_f,
	     "outputs the matrix of blocksizes",std_help);
  result.add("gradings",gradings_f,
	     "prints gradings of imaginary roots for real forms",std_help);
  result.add("strongreal",strongreal_f,
	     "outputs information about strong real forms",std_help);
  // the next function is in main mode since no real form needed
  result.add("dualkgb",dualkgb_f,
	     "prints the KGB data for a dual real form",std_help);
  result.add("help",help_f,"override");
  result.add("q",exitMode,"override"); // q defined but inactive in empty mode

  // add test commands

  test::addTestCommands<MainmodeTag>(result);
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
    (std::cout,current_inner_class().rootDatum().cartanMatrix());

}


// Print information about the root datum (see testprint.cpp for details).
void rootdatum_f()
{
  ioutils::OutputFile file;
  testprint::print(file,current_inner_class().rootDatum());
}


// Print the roots in the lattice basis.
void roots_f()
{
  ioutils::OutputFile file;

  const RootDatum& rd = current_inner_class().rootDatum();

  WeightList::const_iterator first = rd.beginRoot();
  WeightList::const_iterator last = rd.endRoot();
  basic_io::seqPrint(file,first,last,"\n") << std::endl;
}

// Print the coroots in the lattice basis.
void coroots_f()
{
  ioutils::OutputFile file;

  const RootDatum& rd = current_inner_class().rootDatum();

  CoweightList::const_iterator first = rd.beginCoroot();
  CoweightList::const_iterator last = rd.endCoroot();
  basic_io::seqPrint(file,first,last,"\n") << std::endl;
}


// Print the simple roots in the lattice coordinates.
void simpleroots_f()
{
  const RootDatum& rd = current_inner_class().rootDatum();

  rootdata::WRootIterator first = rd.beginSimpleRoot();
  rootdata::WRootIterator last = rd.endSimpleRoot();
  basic_io::seqPrint(std::cout,first,last,"\n") << std::endl;
}

// Print the simple coroots in the lattice coordinates.
void simplecoroots_f()
{
  const RootDatum& rd = current_inner_class().rootDatum();

  rootdata::WRootIterator first = rd.beginSimpleCoroot();
  rootdata::WRootIterator last = rd.endSimpleCoroot();
  basic_io::seqPrint(std::cout,first,last,"\n") << std::endl;
}

// Print the positive roots in the lattice basis.
void posroots_f()
{
  ioutils::OutputFile file;

  const RootDatum& rd = current_inner_class().rootDatum();

  rootdata::WRootIterator first = rd.beginPosRoot();
  rootdata::WRootIterator last = rd.endPosRoot();
  basic_io::seqPrint(file,first,last,"\n") << std::endl;
}

// Print the positive coroots in the lattice basis.
void poscoroots_f()
{
  ioutils::OutputFile file;

  const RootDatum& rd = current_inner_class().rootDatum();

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
  output::printBlockSizes(std::cout,
			  current_inner_class(),currentComplexInterface());
}

// Activates real mode (user will select real form)
void realform_f()
{
  real_mode.activate();
}


void showrealforms_f()
{
  const output::FormNumberMap& rfi =
    currentComplexInterface().realFormInterface();

  std::cout << "(weak) real forms are:" << std::endl;
  output::printRealForms(std::cout,rfi);
}

void showdualforms_f()
{
  const output::FormNumberMap& rfi =
    currentComplexInterface().dualRealFormInterface();

  std::cout << "(weak) dual real forms are:" << std::endl;
  output::printRealForms(std::cout,rfi);
}


// Print the gradings associated to the weak real forms.
void gradings_f()
{
  InnerClass& ic = current_inner_class();

  // get Cartan class; abort if unvalid
  size_t cn=interactive::get_Cartan_class(ic.Cartan_set(ic.quasisplit()));

  ioutils::OutputFile file;

  static_cast<std::ostream&>(file) << std::endl;
  output::printGradings(file,ic,cn,currentComplexInterface())
      << std::endl;

}

// Print information about strong real forms.
void strongreal_f()
{
  InnerClass& ic = current_inner_class();

  // get Cartan class; abort if unvalid
  size_t cn=interactive::get_Cartan_class(ic.Cartan_set(ic.quasisplit()));

  ioutils::OutputFile file;
  file << "\n";
  output::printStrongReal
    (file,
     current_inner_class(),
     currentComplexInterface().realFormInterface(),
     cn);
}

// Print a kgb table for a dual real form.
void dualkgb_f()
{
  InnerClass& ic = current_inner_class();
  output::Interface& G_I = currentComplexInterface();

  RealFormNbr drf = interactive::get_dual_real_form(G_I,ic,ic.numRealForms());

  // the complex group must be in a variable: is non-const for real group
  InnerClass dic(ic,tags::DualTag());
  RealReductiveGroup dG(dic,drf);

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
    InnerClass* G;
    output::Interface* I;
    interactive::get_group_type(G,I,layout,lattice_basis);
    replace_inner_class(G,I);
    drop_to(main_mode); // drop invalidated descendant modes if called from them
  }
  catch(error::InputError& e) {
    e("complex group not changed");
  }

}


} // |namespace|

} // |namespace commands|

} // |namespace atlas|
