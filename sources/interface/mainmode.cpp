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
  void dualkgb_f();
  void help_f();

  // help commands
  void type_h();
  void cmatrix_h();
  void rootdatum_h();
  void roots_h();
  void coroots_h();
  void simpleroots_h();
  void simplecoroots_h();
  void posroots_h();
  void poscoroots_h();
  void realform_h();
  void showdualforms_h();
  void showrealforms_h();
  void blocksizes_h();
  void gradings_h();
  void strongreal_h();
  void dualkgb_h();

  // command tags for the help facility
  const char* cmatrix_tag = "prints the Cartan matrix";
  const char* rootdatum_tag = "outputs the root datum";
  const char* roots_tag = "outputs the roots in the lattice basis";
  const char* coroots_tag = "outputs the coroots in the lattice basis";
  const char* simpleroots_tag =
    "outputs the simple roots in the lattice basis";
  const char* simplecoroots_tag =
    "outputs the simple coroots in the lattice basis";
  const char* posroots_tag = "outputs the positive roots in the lattice basis";
  const char* poscoroots_tag =
    "outputs the positive coroots in the lattice basis";
  const char* realform_tag = "sets the real form for the group";
  const char* showdualforms_tag =
    "outputs the weak real forms for the dual group";
  const char* showrealforms_tag =
    "outputs the weak real forms for this complex group";
  const char* blocksizes_tag = "outputs the matrix of blocksizes";
  const char* gradings_tag="prints gradings of imaginary roots for real forms";
  const char* strongreal_tag = "outputs information about strong real forms";
  const char* dualkgb_tag = "prints the KGB data for a dual real form";

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
  Synopsis: returns a |CommandNode| object that is constructed during the call.
*/
CommandNode mainNode()
{
  CommandNode result ("main: ",main_mode_entry,main_mode_exit);
  result.add("type",type_f,"override"); // tag is a dummy
  result.add("cmatrix",cmatrix_f,cmatrix_tag,cmatrix_h);
  result.add("rootdatum",rootdatum_f,rootdatum_tag,rootdatum_h);
  result.add("roots",roots_f,roots_tag,roots_h);
  result.add("coroots",coroots_f,coroots_tag,coroots_h);
  result.add("simpleroots",simpleroots_f,simpleroots_tag,simpleroots_h);
  result.add("simplecoroots",simplecoroots_f,simplecoroots_tag,simplecoroots_h);
  result.add("posroots",posroots_f,posroots_tag,posroots_h);
  result.add("poscoroots",poscoroots_f,poscoroots_tag,poscoroots_h);
  result.add("realform",realform_f,realform_tag,realform_h);
  result.add("showrealforms",showrealforms_f,showrealforms_tag,showrealforms_h);
  result.add("showdualforms",showdualforms_f,showdualforms_tag,showdualforms_h);
  result.add("blocksizes",blocksizes_f,blocksizes_tag,blocksizes_h);
  result.add("gradings",gradings_f,gradings_tag,gradings_h);
  result.add("strongreal",strongreal_f,strongreal_tag,strongreal_h);
  // the next function is in main mode since no real form needed
  result.add("dualkgb",dualkgb_f,dualkgb_tag,dualkgb_h);
  result.add("help",help_f,"override");
  result.add("q",exitMode,"override"); // q is inactive in empty mode

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
void dualkgb_f()
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


void cmatrix_h()
{
  io::printFile(std::cerr,"cmatrix.help",io::MESSAGE_DIR);
}

void rootdatum_h()
{
  io::printFile(std::cerr,"rootdatum.help",io::MESSAGE_DIR);
}

void roots_h()
{
  io::printFile(std::cerr,"roots.help",io::MESSAGE_DIR);
}

void coroots_h()
{
  io::printFile(std::cerr,"coroots.help",io::MESSAGE_DIR);
}

void simpleroots_h()
{
  io::printFile(std::cerr,"simpleroots.help",io::MESSAGE_DIR);
}

void simplecoroots_h()
{
  io::printFile(std::cerr,"simplecoroots.help",io::MESSAGE_DIR);
}

void posroots_h()
{
  io::printFile(std::cerr,"posroots.help",io::MESSAGE_DIR);
}

void poscoroots_h()
{
  io::printFile(std::cerr,"poscoroots.help",io::MESSAGE_DIR);
}

void realform_h()
{
  io::printFile(std::cerr,"realform.help",io::MESSAGE_DIR);
}

void showrealforms_h()
{
  io::printFile(std::cerr,"showrealforms.help",io::MESSAGE_DIR);
}

void showdualforms_h()
{
  io::printFile(std::cerr,"showdualforms.help",io::MESSAGE_DIR);
}

void blocksizes_h()
{
  io::printFile(std::cerr,"blocksizes.help",io::MESSAGE_DIR);
}

void gradings_h()
{
  io::printFile(std::cerr,"gradings.help",io::MESSAGE_DIR);
}

void strongreal_h()
{
  io::printFile(std::cerr,"strongreal.help",io::MESSAGE_DIR);
}

void dualkgb_h()
{
  io::printFile(std::cerr,"dualkgb.help",io::MESSAGE_DIR);
}

} // |namespace|

} // |namespace commands|

} // |namespace atlas|
