/*
  This is realmode.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#include "realmode.h"

#include <cstdio>   // not obviously used, but appears helpful for Windows

#include "innerclass.h"
#include "output.h"
#include "error.h"
#include "interactive.h"
#include "io.h"
#include "ioutils.h"
#include "realredgp.h"
#include "kgb.h"
#include "kgb_io.h"
#include "test.h"
#include "kgp.h"
#include "repr.h"

#include "mainmode.h"
#include "blockmode.h"
#include "reprmode.h"
#include "helpmode.h"

/****************************************************************************

  This file contains the commands defined in the "real" mode of the program.
  This means that a real form has been chosen.

  NOTE: just as in the main mode, it would be nice to be capable of having
  several active real forms ...

*****************************************************************************/

namespace atlas {

namespace commands {

namespace {

  void real_mode_entry() throw(EntryError);
  void real_mode_exit();

  // functions for the predefined commands

  void components_f();
  void cartan_f();
  void corder_f();
  void realweyl_f();
  void kgb_f();
  void KGB_f();
  void kgborder_f();
  void kgbtwist_f();
  void kgbgraph_f();

  void kgp_f();
  void kgporder_f();
  void kgpgraph_f();

  void realform_f();
  void dualrealform_f();
  void repr_f();

  // help commands
  void KGB_h();

  // command tags for the help facility

  // local variables

  RealReductiveGroup* G_R_pointer=NULL;
  Rep_table* rt=NULL;
} // |namespace|

/*****************************************************************************

        Chapter I -- Functions declared in realmode.h

******************************************************************************/

// Return a |CommandNode| object that is constructed on first call.
CommandNode realNode()
{
  CommandNode result("real: ",real_mode_entry,real_mode_exit);

  result.add("realform",realform_f,"override");
  result.add("components",components_f,
	     "describes component group of the real group",std_help);
  result.add("cartan",cartan_f,
	     "prints the conjugacy classes of Cartan subgroups",std_help);
  result.add("corder",corder_f,
	     "shows Hasse diagram of ordering of Cartan classes",std_help);
  result.add("realweyl",realweyl_f,
	     "outputs the structure of the real Weyl group",std_help);
  result.add("kgb",kgb_f,"prints the orbits of K on G/B",std_help);
  result.add("KGB",KGB_f,
	     "computes KGB data (more information than the kgb command)",KGB_h);
  result.add("kgborder",kgborder_f,
	     "prints the Bruhat ordering on K\\G/B",std_help);
  result.add("kgbtwist",kgbtwist_f,"shows twist orbits on K\\G/B");
  result.add("kgbgraph",kgbgraph_f,
	     "makes a 'dot' file for the Bruhat ordering on K\\G/B",std_help);
  result.add("kgp", kgp_f,"prints the orbits of K on G/P",std_help);
  result.add("kgporder", kgporder_f,
	     "prints the Bruhat ordering on K\\G/P",std_help);
  result.add("kgpgraph", kgpgraph_f,
	     "makes a 'dot' file for the Bruhat ordering on K\\G/P",std_help);
  result.add("dualrealform",dualrealform_f,"sets the dual real form",use_tag);
  result.add("repr",repr_f,"sets the parameter for a representation",use_tag);

  // add test commands
  test::addTestCommands<RealmodeTag>(result);
  return result;
}

RealReductiveGroup& currentRealGroup()
{
  return *G_R_pointer;
}

RealFormNbr currentRealForm()
{
  return G_R_pointer->realForm();
}

const Rep_context& currentRepContext() { return *rt; }
Rep_table& currentRepTable() { return *rt; }


/****************************************************************************

        Chapter II -- The real mode |CommandNode|

  One instance of |CommandNode| for the real mode is created at the
  first call of |realMode()|; further calls just return a reference to it.

*****************************************************************************/

namespace {

/*
  Synopsis: attempts to set a real form interactively. In case of failure,
  throws an InputError and returns.
*/
void real_mode_entry() throw(EntryError)
{
  try
  {
    RealFormNbr rf = interactive::get_real_form(currentComplexInterface());
    G_R_pointer=new RealReductiveGroup(current_inner_class(),rf);
    rt = new Rep_table(currentRealGroup());
  }
  catch(error::InputError& e)
  {
    real_mode_exit(); // clean up
    e("no real form set");
    throw EntryError();
  }
}

/*
  Reset the real form, effectively re-entering the real mode. If the choice
  of a new real form fails, the current real form remains in force.
*/
void realform_f()
{
  try
  { // we can call the swap method for rvalues, but not with and rvalue arg
    RealFormNbr rf = interactive::get_real_form(currentComplexInterface());
  RealReductiveGroup(current_inner_class(),rf).swap(currentRealGroup());
    delete rt; rt = new Rep_table(currentRealGroup());
    drop_to(real_mode); // drop invalidated descendant modes if called from them
  }
  catch (error::InputError& e)
  {
    e("real form not changed");
  }
}



/*
  Synopsis: destroys the real group and its interface, resoring NULL pointers
*/
void real_mode_exit()
{
  delete G_R_pointer; G_R_pointer=NULL;
  delete rt; rt=NULL;
}



/*****************************************************************************

        Chapter III --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/

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
  output::printCartanClasses(file,currentRealGroup(),
				   currentComplexInterface())
    << std::endl;
}

// Print the Hasse diagram of the ordering of Cartan classes.
void corder_f()
{
  RealReductiveGroup& G_R = currentRealGroup();

  std::cout << "Hasse diagram of Cartan class ordering:" << std::endl;
  output::printCartanOrder(std::cout,G_R);
}


// enter block mode
void dualrealform_f()
{
  drop_to(real_mode); // ensure we have left any descendent mode
  block_mode.activate();
}

// enter repr mode
void repr_f()
{
  drop_to(real_mode); // ensure we have left any descendent mode
  repr_mode.activate();
}

// Show the structure of the real weyl group.
void realweyl_f()
{
  // Cartan classes will be generated by call to |currentRealGroup|

  // get Cartan class; abort if unvalid
  size_t cn=interactive::get_Cartan_class(currentRealGroup().Cartan_set());

  ioutils::OutputFile file;
  file << "\n";
  output::printRealWeyl(file,currentRealGroup(),cn);
}


// Print the kgb table
void kgb_f()
{
  RealReductiveGroup& G_R = currentRealGroup();
  std::cout << "kgbsize: " << G_R.KGB_size() << std::endl;
  ioutils::OutputFile file;
  KGB kgb(G_R,G_R.Cartan_set());
  kgb_io::printKGB(file,kgb); // generate traditionally, and forget
}

// Same, but with more info and non-traditional generation
void KGB_f()
{
  RealReductiveGroup& G_R = currentRealGroup();
  std::cout << "kgbsize: " << G_R.KGB_size() << std::endl;

  ioutils::OutputFile f;

  InnerClass& G=G_R.innerClass();
  kgb_io::var_print_KGB(f,G,G_R.kgb());
}


// Print the Hasse diagram of the ordering of K orbits on G/B.
void kgborder_f()
{
  RealReductiveGroup& G = currentRealGroup();

  std::cout << "kgbsize: " << G.KGB_size() << std::endl;
  ioutils::OutputFile file;

  kgb_io::printBruhatOrder(file,G.Bruhat_KGB());
}

void kgbtwist_f()
{
  RealReductiveGroup& G = currentRealGroup();

  ioutils::OutputFile file;

  kgb_io::print_twist(file,G.kgb());
}


void kgbgraph_f()
{
  RealReductiveGroup& G_R = currentRealGroup();
  std::cout << "kgbsize: " << G_R.KGB_size() << std::endl;

  ioutils::OutputFile file;
  kgb_io::makeDotFile(file,G_R.kgb(),G_R.Bruhat_KGB());
}

void kgp_f()
{
  // build KGP
  RealReductiveGroup& G_R = currentRealGroup();
  Parabolic psg;
  interactive::getInteractive(psg, G_R.rank());
  kgb::KGP kgp(G_R,psg);

  // print the size and simple roots
  std::cout << "kgp size for roots {";
  bool first = true;
  for (Parabolic::iterator it=psg.begin(); it(); ++it)
  {
    if (first)
      first=false;
    else
      std::cout << ",";
    std::cout << *it+1;
  }
  std::cout << "}: " << kgp.size() << std::endl;

  ioutils::OutputFile file;
  kgp.print(file);
}

void kgporder_f()
{
  // build KGP
  RealReductiveGroup& G_R = currentRealGroup();
  atlas::Parabolic psg;
  interactive::getInteractive(psg, G_R.rank());
  kgb::KGP kgp(G_R,psg);

  // compute closure
  kgp.fillClosure();

  // print the size and simple roots
  std::cout << "kgp size for roots {";
  bool first = true;
  for (Parabolic::iterator it=psg.begin(); it(); ++it)
  {
    if (first)
      first=false;
    else
      std::cout << ",";
    std::cout << *it+1;
  }
  std::cout << "}: " << kgp.size() << std::endl;

  ioutils::OutputFile file;
  kgp.printClosure(file);
}

void kgpgraph_f()
{
  // build KGP and fill closure
  RealReductiveGroup& G_R = currentRealGroup();
  atlas::Parabolic psg;
  interactive::getInteractive(psg, G_R.rank());
  kgb::KGP kgp(G_R,psg);
  kgp.fillClosure();

  // print the size and simple roots
  std::cout << "kgp size for roots {";
  bool first = true;
  for (Parabolic::iterator it=psg.begin(); it(); ++it)
  {
    if (first)
      first=false;
    else
      std::cout << ",";
    std::cout << *it+1;
  }
  std::cout << "}: " << kgp.size() << std::endl;

  ioutils::OutputFile file;
  if (file.is_std_cout()) // make sure the user entered an actual filename
    throw error::InputError(); // as standard output makes no sense here

  // make the dot file
  kgp.makeDotFile(file);
}



void KGB_h()
{
  io::printFile(std::cerr,"KGB_.help",io::MESSAGE_DIR);
}


} // |namespace|

} // |namespace commands|

} // |namespace atlas|
