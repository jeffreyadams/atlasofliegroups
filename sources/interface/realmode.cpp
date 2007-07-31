/*
  This is realmode.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  See file main.cpp for full copyright notice
*/

#include "realmode.h"

#include "cartan_io.h"
#include "complexredgp.h"
#include "complexredgp_io.h"
#include "error.h"
#include "helpmode.h"
#include "interactive.h"
#include "io.h"
#include "ioutils.h"
#include "mainmode.h"
#include "realredgp.h"
#include "realredgp_io.h"
#include "special.h"
#include "test.h"

/****************************************************************************

  This file contains the commands defined in the "real" mode of the program.
  This means that a real form has been chosen.

  NOTE: just as in the main mode, it would be nice to be capable of having
  several active real forms ...

*****************************************************************************/

namespace atlas {

namespace {

  using namespace realmode;

  class ThisMode:public commands::CommandMode {
  public:
    ThisMode();
    ~ThisMode();
  };

  void this_entry() throw(commands::EntryError);
  void this_exit();

  // functions for the predefined commands

  void cartan_f();
  void gradings_f();
  void help_f();
  void q_h();
  void realform_f();
  void realweyl_f();
  void strongreal_f();
  void type_f();

  // local variables

  realredgp::RealReductiveGroup G_R;
  realredgp_io::Interface G_RI;

}

/****************************************************************************

        Chapter I -- The ThisMode class

  Only one instance of this class will be constructed, on the first call
  of the function realMode().

*****************************************************************************/

namespace {

ThisMode::ThisMode()
  :commands::CommandMode("real: ",this_entry,this_exit)

{
  // set parent mode
  d_prev = &mainmode::mainMode();

  // add the commands from the main mode
  commands::addCommands(*this,prev());

  // add commands for this mode
  add("cartan",cartan_f);
  add("gradings",gradings_f);
  add("help",help_f);
  add("q",commands::exitMode);
  add("realform",realform_f);
  add("realweyl",realweyl_f);
  add("strongreal",strongreal_f);
  // the "type" command should be redefined here because it needs to exit
  // the real mode
  add("type",type_f);

  // add special commands
  special::addSpecialCommands(*this,RealmodeTag());

  // add test commands
  test::addTestCommands(*this,RealmodeTag());
}

ThisMode::~ThisMode()

{}

void this_entry() throw(commands::EntryError)

/*
  Synopsis: attempts to set a real form interactively. In case of failure,
  throws an InputError and returns.
*/

{
  using namespace interactive;
  using namespace mainmode;

  try {
    getInteractive(G_R,currentComplexInterface());
  }
  catch(error::InputError& e) {
    e("no real form set");
    throw commands::EntryError();
  }

  realredgp_io::Interface RI(G_R,currentComplexInterface());
  G_RI.swap(RI);

  return;
}

void this_exit()

/*
  Synopsis: resets the real form to the default (meaningless!) real form,
  and destructs the current one.
*/

{
  using namespace realredgp;

  RealReductiveGroup newG;
  G_R.swap(newG);

  return;
}

}

/*****************************************************************************

        Chapter II --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/

namespace {

void cartan_f()

/*
  Synopsis: prints the cartan classes of G_R.
*/

{
  using namespace mainmode;
  using namespace realredgp_io;

  try {
    G_R.fillCartan();
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  ioutils::OutputFile file;

  (std::ostream&)(file) << std::endl;
  printCartanClasses(file,G_RI) << std::endl;

  return;
}

void gradings_f()

/*
  Synopsis: prints the gradings associated to the weak real forms.
*/

{
  using namespace complexredgp_io;
  using namespace commands;
  using namespace input;
  using namespace interactive;
  using namespace mainmode;

  try {
    G_R.fillCartan();
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  size_t cn;

  // get Cartan class; abort if unvalid
  try {
    getCartanClass(cn,G_R.cartanSet(),currentLine());
  }
  catch (error::InputError& e) {
    e("aborted");
    return;
  }

  ioutils::OutputFile file;

  (std::ostream&)(file) << std::endl;
  printGradings(file,cn,G_RI.complexInterface()) << std::endl;

  return;
}

void help_f()

{
  activate(helpmode::helpMode());
  return;
}

void q_h()

/*
  Synopsis: exits the current mode.
*/

{
  commands::exitMode();
  return;
}

void realform_f()

/*
  Synopsis: resets the type, effectively re-entering the real mode. If the
  construction of the new type fails, the current type remains in force.
*/

{
  using namespace interactive;
  using namespace mainmode;

  try {
    getInteractive(G_R,currentComplexInterface());
  }
  catch (error::InputError& e) {
    e("real form not changed");
    return;
  }
  catch (error::InnerClassError& e) {
    e("that real form is not in the current inner class\n\
real form not changed");
    return;
  }

  realredgp_io::Interface RI(G_R,currentComplexInterface());
  G_RI.swap(RI);

  return;
}

void realweyl_f()

/*
  Synopsis: prints out the structure of the real weyl group.
*/

{
  using namespace commands;
  using namespace input;
  using namespace interactive;
  using namespace ioutils;
  using namespace realmode;
  using namespace realredgp;
  using namespace realredgp_io;

  try {
    G_R.fillCartan();
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  size_t cn;

  // get Cartan class; abort if unvalid
  try {
    getCartanClass(cn,G_R.cartanSet(),currentLine());
  }
  catch (error::InputError& e) {
    e("aborted");
    return;
  }

  OutputFile file;
  file << "\n";
  printRealWeyl(file,G_R,cn);

  return;
}

void strongreal_f()

/*
  Synopsis: outputs information about strong real forms.
*/

{

  using namespace commands;
  using namespace input;
  using namespace interactive;
  using namespace ioutils;
  using namespace realmode;
  using namespace realredgp;
  using namespace realredgp_io;

  try {
    G_R.fillCartan();
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  size_t cn;

  // get Cartan class; abort if unvalid
  try {
    getCartanClass(cn,G_R.cartanSet(),currentLine());
  }
  catch (error::InputError& e) {
    e("aborted");
    return;
  }

  OutputFile file;
  file << "\n";
  printStrongReal(file,G_R,G_RI.realFormInterface(),cn);

  return;
}

void type_f()

/*
  Synopsis: resets the type of the complex group.

  In case of success, the real forms are invalidated, and therefore we
  should exit real mode; in case of failure, we don't need to.
*/

{
  using namespace complexredgp_io;
  using namespace commands;
  using namespace interactive;
  using namespace mainmode;

  try {
    complexredgp::ComplexReductiveGroup* G;
    complexredgp_io::Interface* I;

    interactive::getInteractive(G,I);
    mainmode::replaceComplexGroup(G,I);
    exitMode();
  }
  catch (error::InputError& e) {
    e("complex group and real form not changed");
  }
}

}

/*****************************************************************************

        Chapter III -- Functions declared in realmode.h

  The following functions are defined here :

    - currentRealGroup() : exports G_R to the special and test modules;

******************************************************************************/

namespace realmode {

realredgp::RealReductiveGroup& currentRealGroup()

{
  return G_R;
}

realredgp_io::Interface& currentRealInterface()

{
  return G_RI;
}

commands::CommandMode& realMode()

/*
  Synopsis: returns the ThisMode object.

  It is constructed on first call.
*/

{
  static ThisMode mode;
  return mode;
}

}

}
