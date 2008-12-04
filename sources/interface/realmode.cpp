/*
  This is realmode.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For copyright and license information see the LICENSE file
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
#include "kgb.h"
#include "kgb_io.h"
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

  void real_mode_entry() throw(commands::EntryError);
  void real_mode_exit();

  // functions for the predefined commands

  void components_f();
  void cartan_f();
  void corder_f();
  void gradings_f();
  void realform_f();
  void realweyl_f();
  void strongreal_f();
  void kgb_f();
  void kgborder_f();
  void type_f();


  // local variables

  realredgp::RealReductiveGroup* G_R_pointer=NULL;
  realredgp_io::Interface* RI_pointer=NULL;

}

/*****************************************************************************

        Chapter I -- Functions declared in realmode.h

******************************************************************************/

namespace realmode {

// Return a |CommandMode| object that is constructed on first call.
commands::CommandMode& realMode()
{
  static commands::CommandMode real_mode
    ("real: ",real_mode_entry,real_mode_exit);
  if (real_mode.empty()) // true upon first call
  {
    // add the commands from the main mode
    commands::addCommands(real_mode,mainmode::mainMode());

    // add commands for this mode
    real_mode.add("components",components_f);
    real_mode.add("cartan",cartan_f);
    real_mode.add("corder",corder_f);
    real_mode.add("gradings",gradings_f);
    real_mode.add("realform",realform_f);
    real_mode.add("realweyl",realweyl_f);
    real_mode.add("strongreal",strongreal_f);
    real_mode.add("kgb",kgb_f);
    real_mode.add("kgborder",kgborder_f);
    // the "type" command should be redefined here because it needs to exit
    // the real mode
    real_mode.add("type",type_f); // override

    // add special commands
    special::addSpecialCommands(real_mode,RealmodeTag());

    // add test commands
    test::addTestCommands(real_mode,RealmodeTag());
  }
  return real_mode;
}

realredgp::RealReductiveGroup& currentRawRealGroup()
{
  return *G_R_pointer;
}

realredgp::RealReductiveGroup& currentRealGroup()
{
  mainmode::currentComplexGroup().fillCartan(); // generate ALL Cartan classes
  return *G_R_pointer;
}

realform::RealForm currentRawRealForm()
{
  return G_R_pointer->realForm();
}

realform::RealForm currentRealForm()
{
  mainmode::currentComplexGroup().fillCartan(); // generate ALL Cartan classes
  return G_R_pointer->realForm();
}

realredgp_io::Interface& currentRealInterface()
{
  return *RI_pointer;
}


} // namespace realmode


/****************************************************************************

        Chapter II -- The real mode |CommandMode|

  One instance of |CommandMode| for the real mode is created at the
  first call of |realMode()|; further calls just return a reference to it.

*****************************************************************************/

namespace {

/*
  Synopsis: attempts to set a real form interactively. In case of failure,
  throws an InputError and returns.
*/
void real_mode_entry() throw(commands::EntryError)
{
  try {
    G_R_pointer=new realredgp::RealReductiveGroup;
    interactive::getInteractive(currentRawRealGroup(),
				mainmode::currentComplexInterface());

    RI_pointer=new realredgp_io::Interface
      (currentRawRealGroup(),mainmode::currentComplexInterface());
  }
  catch(error::InputError& e) {
    real_mode_exit(); // clean up, even if unnecessary for |RI_pointer|
    e("no real form set");
    throw commands::EntryError();
  }
}


/*
  Synopsis: destroys the real group and its interface, resoring NULL pointers
*/
void real_mode_exit()
{
  delete G_R_pointer; delete RI_pointer; G_R_pointer=NULL; RI_pointer=NULL;
}


} // namespace

/*****************************************************************************

        Chapter III --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/

namespace {

/*
  Print the (dual) component group of the current group. We print it out
  in terms of the canonical basis of $T(2)^\vee$
*/
void components_f()
{
  const realredgp::RealReductiveGroup& G = currentRawRealGroup();
  const latticetypes::ComponentList& c = G.dualComponentReps();

  if (c.size() > 0)
    std::cout << "component group is (Z/2)^" << c.size() << std::endl;
  else
    std::cout << "group is connected" << std::endl;

}


// Print the cartan classes of G_R.
void cartan_f()
{
  mainmode::currentComplexGroup().fillCartan(); // generate ALL Cartan classes

  ioutils::OutputFile file;

  static_cast<std::ostream&>(file) << std::endl;
  realredgp_io::printCartanClasses(file,currentRealInterface()) << std::endl;
}

// Print the Hasse diagram of the ordering of Cartan classes.
void corder_f()
{
  realredgp::RealReductiveGroup& G_R = currentRealGroup();

  std::cout << "Hasse diagram of Cartan class ordering:" << std::endl;
  realredgp_io::printCartanOrder(std::cout,G_R);
}


// Print the gradings associated to the weak real forms.
void gradings_f()
{
  // Cartan classes will be generated by call to |currentRealGroup|

  // get Cartan class; abort if unvalid
  size_t cn=interactive::get_Cartan_class(currentRealGroup().Cartan_set());

  ioutils::OutputFile file;

  static_cast<std::ostream&>(file) << std::endl;
  complexredgp_io::printGradings(file,cn,
				 currentRealInterface().complexInterface())
      << std::endl;

}

/*
  Reset the type, effectively re-entering the real mode. If the
  construction of the new type fails, the current real form remains in force.
*/
void realform_f()
{
  try {
    interactive::getInteractive
      (currentRawRealGroup(),mainmode::currentComplexInterface());

    realredgp_io::Interface RI
      (currentRawRealGroup(),mainmode::currentComplexInterface());
    currentRealInterface().swap(RI);

  }
  catch (error::InputError& e) {
    e("real form not changed");
  }
  catch (error::InnerClassError& e) {
    e("that real form is not in the current inner class\n"
      "real form not changed");
  }
}

// Show the structure of the real weyl group.
void realweyl_f()
{
  // Cartan classes will be generated by call to |currentRealGroup|

  // get Cartan class; abort if unvalid
  size_t cn=interactive::get_Cartan_class(currentRealGroup().Cartan_set());

  ioutils::OutputFile file;
  file << "\n";
  realredgp_io::printRealWeyl(file,currentRawRealGroup(),cn);
}


// Print information about strong real forms.
void strongreal_f()
{
  // Cartan classes will be generated by call to |currentRealGroup|

  // get Cartan class; abort if unvalid
  size_t cn=interactive::get_Cartan_class(currentRealGroup().Cartan_set());

  ioutils::OutputFile file;
  file << "\n";
  realredgp_io::printStrongReal(file,
				currentRawRealGroup(),
				currentRealInterface().realFormInterface(),
				cn);
}


// Print the kgb table
void kgb_f()
{
  realredgp::RealReductiveGroup& G_R = currentRealGroup();
  std::cout << "kgbsize: " << G_R.KGB_size() << std::endl;
  ioutils::OutputFile file;
  kgb::KGB kgb(G_R);
  kgb_io::printKGB(file,kgb);
}

// Print the Hasse diagram of the ordering of K orbits on G/B.
void kgborder_f()
{
  realredgp::RealReductiveGroup& G = currentRealGroup();

  std::cout << "kgbsize: " << G.KGB_size() << std::endl;
  ioutils::OutputFile file;

  kgb::KGB kgb(G);
  kgb_io::printBruhatOrder(file,kgb.bruhatOrder());
}

/* Reset the type of the complex group.

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
    commands::exitMode(); // upon success pop real mode, destroying real group
  }
  catch (error::InputError& e) {
    e("complex group and real form not changed");
  }
}


} // namespace

} // namespace atlas
