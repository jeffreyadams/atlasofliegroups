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
  void realweyl_f();
  void kgb_f();
  void KGB_f();
  void kgborder_f();

  void type_f();
  void realform_f();


  // local variables

  RealReductiveGroup* G_R_pointer=NULL;

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
    real_mode.add("realform",realform_f);
    real_mode.add("realweyl",realweyl_f);
    real_mode.add("kgb",kgb_f);
    real_mode.add("KGB",KGB_f);
    real_mode.add("kgborder",kgborder_f);
    // the "type" command should be redefined here because it needs to exit
    // the real mode
    real_mode.add("type",type_f); // override

    // add test commands
    test::addTestCommands(real_mode,RealmodeTag());
  }
  return real_mode;
}

RealReductiveGroup& currentRealGroup()
{
  return *G_R_pointer;
}

RealFormNbr currentRealForm()
{
  return G_R_pointer->realForm();
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
  try
  {
    G_R_pointer=new RealReductiveGroup
      (interactive::getRealGroup(mainmode::currentComplexInterface()));
  }
  catch(error::InputError& e)
  {
    real_mode_exit(); // clean up
    e("no real form set");
    throw commands::EntryError();
  }
}


/*
  Synopsis: destroys the real group and its interface, resoring NULL pointers
*/
void real_mode_exit()
{
  delete G_R_pointer; G_R_pointer=NULL;
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
  const RealReductiveGroup& G = currentRealGroup();
  size_t r = G.component_rank();

  if (r>0)
    std::cout << "component group is (Z/2)^" << r << std::endl;
  else
    std::cout << "real group is connected" << std::endl;

}


// Print the cartan classes of G_R.
void cartan_f()
{
  ioutils::OutputFile file;

  static_cast<std::ostream&>(file) << std::endl;
  realredgp_io::printCartanClasses(file,currentRealForm(),
				   mainmode::currentComplexInterface())
    << std::endl;
}

// Print the Hasse diagram of the ordering of Cartan classes.
void corder_f()
{
  RealReductiveGroup& G_R = currentRealGroup();

  std::cout << "Hasse diagram of Cartan class ordering:" << std::endl;
  realredgp_io::printCartanOrder(std::cout,G_R);
}


/*
  Reset the type, effectively re-entering the real mode. If the
  construction of the new type fails, the current real form remains in force.
*/
void realform_f()
{
  try
  { // we can call the swap method for rvalues, but not with and rvalue arg
    interactive::getRealGroup(mainmode::currentComplexInterface()).swap
      (realmode::currentRealGroup());
  }
  catch (error::InputError& e)
  {
    e("real form not changed");
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
  realredgp_io::printRealWeyl(file,currentRealGroup(),cn);
}


// Print the kgb table
void kgb_f()
{
  RealReductiveGroup& G_R = currentRealGroup();
  std::cout << "kgbsize: " << G_R.KGB_size() << std::endl;
  ioutils::OutputFile file;
  KGB kgb(G_R);
  kgb_io::printKGB(file,kgb); // generate traditionally, and forget
}

// Same, but with more info and non-traditional generation
void KGB_f()
{
  RealReductiveGroup& G_R = realmode::currentRealGroup();

  ioutils::OutputFile f;

  kgb_io::var_print_KGB(f,mainmode::currentComplexGroup(),G_R.kgb());
}

// Print the Hasse diagram of the ordering of K orbits on G/B.
void kgborder_f()
{
  RealReductiveGroup& G = currentRealGroup();

  std::cout << "kgbsize: " << G.KGB_size() << std::endl;
  ioutils::OutputFile file;

  kgb_io::printBruhatOrder(file,G.Bruhat_KGB());
}

/* Reset the type of the complex group.

  In case of success, the real forms are invalidated, and therefore we
  should exit real mode; in case of failure, we don't need to.
*/
void type_f()
{
  try {
    ComplexReductiveGroup* G;
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
