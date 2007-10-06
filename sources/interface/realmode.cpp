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

/*
  Synopsis: prints the cartan classes of G_R.
*/
void cartan_f()
{
  try {
    G_R.fillCartan();

    ioutils::OutputFile file;

    static_cast<std::ostream&>(file) << std::endl;
    realredgp_io::printCartanClasses(file,G_RI) << std::endl;
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("abort");
  }
}

/*
  Synopsis: prints the gradings associated to the weak real forms.
*/
void gradings_f()
{
  try {
    G_R.fillCartan();

    size_t cn;

    // get Cartan class; abort if unvalid
    interactive::getCartanClass(cn,G_R.cartanSet(),commands::currentLine());

    ioutils::OutputFile file;

    static_cast<std::ostream&>(file) << std::endl;
    complexredgp_io::printGradings(file,cn,G_RI.complexInterface())
      << std::endl;
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

void help_f()
{
  activate(helpmode::helpMode());
}

/*
  Synopsis: exits the current mode.
*/
void q_h()
{
  commands::exitMode();
}

/*
  Synopsis: resets the type, effectively re-entering the real mode. If the
  construction of the new type fails, the current type remains in force.
*/
void realform_f()
{
  try {
    interactive::getInteractive(G_R,mainmode::currentComplexInterface());

    realredgp_io::Interface RI(G_R,mainmode::currentComplexInterface());
    G_RI.swap(RI);

  }
  catch (error::InputError& e) {
    e("real form not changed");
  }
  catch (error::InnerClassError& e) {
    e("that real form is not in the current inner class\n"
      "real form not changed");
  }
}

/*
  Synopsis: prints out the structure of the real weyl group.
*/
void realweyl_f()
{
  try {
    G_R.fillCartan();

    size_t cn;

    // get Cartan class; abort if unvalid
    interactive::getCartanClass(cn,G_R.cartanSet(),commands::currentLine());

    ioutils::OutputFile file;
    file << "\n";
    realredgp_io::printRealWeyl(file,G_R,cn);

  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }
}


/*
  Synopsis: outputs information about strong real forms.
*/
void strongreal_f()
{
  try {
    G_R.fillCartan();

    size_t cn;

    // get Cartan class; abort if unvalid
    interactive::getCartanClass(cn,G_R.cartanSet(),commands::currentLine());

    ioutils::OutputFile file;
    file << "\n";
    realredgp_io::printStrongReal(file,G_R,G_RI.realFormInterface(),cn);
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }
}


/*
  Synopsis: resets the type of the complex group.

  In case of success, the real forms are invalidated, and therefore we
  should exit real mode; in case of failure, we don't need to.
*/
void type_f()
{
  try {
    complexredgp::ComplexReductiveGroup* G;
    complexredgp_io::Interface* I;

    interactive::getInteractive(G,I);
    mainmode::replaceComplexGroup(G,I);
    commands::exitMode();
  }
  catch (error::InputError& e) {
    e("complex group and real form not changed");
  }
}

}

/*****************************************************************************

        Chapter III -- Functions declared in realmode.h

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

/*
  Synopsis: returns the ThisMode object.

  It is constructed on first call.
*/
commands::CommandMode& realMode()
{
  static ThisMode mode;
  return mode;
}

}

}
