/*
  This is test.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "test.h"

#include <cassert>
#include <iostream>
#include <fstream>

#include "helpmode.h"
#include "emptymode.h"
#include "mainmode.h"
#include "realmode.h"
#include "blockmode.h"

#include "commands.h"
#include "rootdata.h"
#include "complexredgp.h"
#include "realredgp.h"
#include "interactive.h"
#include "ioutils.h"
#include "prettyprint.h"
#include "blocks.h"
#include "kgb.h"
#include "kgb_io.h"
#include "klsupport.h"
#include "kl.h"
#include "kltest.h"
#include "standardrepk.h"
#include "free_abelian.h"

/*****************************************************************************

  This module contains some commands for testing the program.

******************************************************************************/

namespace atlas {

namespace {
  using namespace test;

  // functions for the test commands

  void coroots_rootbasis_f();
  void poscoroots_rootbasis_f();
  void posroots_rootbasis_f();
  void roots_rootbasis_f();
  void KGB_f();
  void sub_KGB_f();
  void trivial_f();
  void test_f();

  // help functions


  // tags

  const char* test_tag = "(test command)";

/*
  For convenience, the "test" command is added to the mode that is flagged by
  the testMode constant defined here; therefore "test" appears (conditionally)
  in every template instance of |addTestCommands| and of |addTestHelp| below.
  Set this constant according to the requirements of the |test_f| function.
*/
  enum TestMode {EmptyMode, MainMode, RealMode, BlockMode, numTestMode};
  const TestMode testMode = RealMode; // currently does a KL matrix test

  // utilities
  const rootdata::RootDatum& currentRootDatum();

}

/*****************************************************************************

        Chapter I -- Functions declared by test.h

  This section defines the functions declared in test.h :

    - addTestCommands() : adds the test commands to the main command
      tree;
    - addTestHelp() : adds help functionality;

******************************************************************************/

namespace test {

/* |addTestCommands| is a template only so that its declaration is shorter;
   its defining instances are defined just as overloaded functions would,
   since each of them needs to test a specific value of |testMode|.
*/


// Add to the empty mode the test commands that require that mode.
template<>
void addTestCommands<emptymode::EmptymodeTag>
  (commands::CommandMode& mode, emptymode::EmptymodeTag)
{
  if (testMode == EmptyMode)
    mode.add("test",test_f);

}


// Add to the main mode the test commands that require that mode.
template<>
void addTestCommands<mainmode::MainmodeTag>
  (commands::CommandMode& mode, mainmode::MainmodeTag)
{
  if (testMode == MainMode)
    mode.add("test",test_f);

  // add additional commands here :

  mode.add("roots_rootbasis",roots_rootbasis_f);
  mode.add("posroots_rootbasis",posroots_rootbasis_f);
  mode.add("coroots_rootbasis",coroots_rootbasis_f);
  mode.add("poscoroots_rootbasis",poscoroots_rootbasis_f);
}


// Add to the real mode the test commands that require that mode.
template<>
void addTestCommands<realmode::RealmodeTag>
  (commands::CommandMode& mode, realmode::RealmodeTag)
{
  if (testMode == RealMode)
    mode.add("test",test_f);

  // add additional commands here :

  mode.add("KGB",KGB_f);
  mode.add("sub_KGB",sub_KGB_f);
  mode.add("trivial",trivial_f);

}

// Add to the block mode the test commands that require that mode.
template<>
void addTestCommands<blockmode::BlockmodeTag>
  (commands::CommandMode& mode, blockmode::BlockmodeTag)
{
  if (testMode == BlockMode)
    mode.add("test",test_f);

}


// Add to the help mode the test commands that require that mode.
template<> void addTestHelp<emptymode::EmptymodeTag>
             (commands::CommandMode& mode, commands::TagDict& t,
	      emptymode::EmptymodeTag)
{
  using namespace commands;
  using namespace helpmode;

  if (testMode == EmptyMode) {
    mode.add("test",helpmode::nohelp_h);
    commands::insertTag(t,"test",test_tag);
  }


  // add additional help commands here:

  // add additional command tags here:

}


// Add to the main mode the help commands for test commands with that mode
template<> void addTestHelp<mainmode::MainmodeTag>
             (commands::CommandMode& mode, commands::TagDict& t,
	      mainmode::MainmodeTag)
{
  using namespace commands;
  using namespace helpmode;

  if (testMode == MainMode) {
    mode.add("test",helpmode::nohelp_h);
    insertTag(t,"test",test_tag);
  }

  // add additional help commands here:

  // add additional command tags here:

  mode.add("rootdatum",helpmode::nohelp_h);
  mode.add("roots_rootbasis",helpmode::nohelp_h);
  mode.add("posroots_rootbasis",helpmode::nohelp_h);
  mode.add("coroots_rootbasis",helpmode::nohelp_h);
  mode.add("poscoroots_rootbasis",helpmode::nohelp_h);


  // add additional command tags here :

  insertTag(t,"coroots_rootbasis",test_tag);
  insertTag(t,"poscoroots_rootbasis",test_tag);
  insertTag(t,"posroots_rootbasis",test_tag);
  insertTag(t,"roots_rootbasis",test_tag);
  insertTag(t,"rootdatum",test_tag);

}


// Add to the real mode the help commands for test commands with that mode
template<> void addTestHelp<realmode::RealmodeTag>
             (commands::CommandMode& mode, commands::TagDict& t,
	      realmode::RealmodeTag)
{
  using namespace commands;
  using namespace helpmode;

  if (testMode == RealMode) {
    mode.add("test",helpmode::nohelp_h);
    insertTag(t,"test",test_tag);
  }

  mode.add("KGB",helpmode::nohelp_h);
  mode.add("sub_KGB",helpmode::nohelp_h);
  mode.add("trivial",helpmode::nohelp_h);


  // add additional command tags here :

  insertTag(t,"KGB",test_tag);
  insertTag(t,"sub_KGB",test_tag);
  insertTag(t,"trivial",test_tag);

}

// Add to the block mode the help commands for test commands with that mode
template<> void addTestHelp<blockmode::BlockmodeTag>
             (commands::CommandMode& mode, commands::TagDict& t,
	      blockmode::BlockmodeTag)
{
  using namespace commands;
  using namespace helpmode;

  if (testMode == BlockMode) {
    mode.add("test",helpmode::nohelp_h);
    insertTag(t,"test",test_tag);
  }

  // add additional help commands here:



  // add additional command tags here:

}

} // namespace test

namespace {






}


/*****************************************************************************

        Chapter II -- Functions for the test commands

******************************************************************************/

namespace {

  // Main mode functions


// Print the roots in the simple root coordinates.
void roots_rootbasis_f()
{
  try {
    const rootdata::RootDatum& rd=mainmode::currentComplexGroup().rootDatum();
    ioutils::OutputFile file;

    for (rootdata::RootNbr i=0; i<rd.numRoots(); ++i)
      prettyprint::printInRootBasis(file,i,rd) << std::endl;
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

// Print the positive roots in the simple root coordinates.
void posroots_rootbasis_f()

{
  try {
    const rootdata::RootDatum& rd=mainmode::currentComplexGroup().rootDatum();
    ioutils::OutputFile file;

    prettyprint::printInRootBasis(file,rd.posRootSet(),rd);
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

// Print the coroots in the simple coroot coordinates.
void coroots_rootbasis_f()
{
  try {
    const rootdata::RootDatum rd (mainmode::currentComplexGroup().rootDatum(),
				  tags::DualTag());
    ioutils::OutputFile file;

    for (rootdata::RootNbr i=0; i<rd.numRoots(); ++i)
      prettyprint::printInRootBasis(file,i,rd) << std::endl;
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

// Print the positive coroots in the simple coroot coordinates.
void poscoroots_rootbasis_f()
{
  try {
    const rootdata::RootDatum rd (mainmode::currentComplexGroup().rootDatum(),
				  tags::DualTag());
    ioutils::OutputFile file;

    prettyprint::printInRootBasis(file,rd.posRootSet(),rd);
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

// Real mode functions

void KGB_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();
  G_R.fillCartan(); // must not forget this!
  kgb::KGB kgb(G_R,G_R.cartanSet());
  kgb_io::var_print_KGB(std::cout,mainmode::currentComplexGroup(),kgb);
}

void sub_KGB_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();
  G_R.fillCartan(); // must not forget this!
  kgb::KGB kgb(G_R,G_R.cartanSet());
  unsigned long cn;
  interactive::getInteractive(cn,"Cartan class: ",G_R.cartanSet());
  standardrepk::KHatComputations khc(G_R,kgb);

  weyl::WeylWord ww;
  standardrepk::PSalgebra q=
    khc.theta_stable_parabolic(ww,G_R.complexGroup().cartanClasses(),cn);
  kgb::KGBEltList sub=khc.sub_KGB(q);

  std::cout << "Conjugating word [" << ww << "]\n";
  kgb_io::print_sub_KGB(std::cout,kgb,sub);
}

void trivial_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();
  G_R.fillCartan(); // must not forget this!
  const rootdata::RootDatum& rd=G_R.rootDatum();

  kgb::KGB kgb(G_R,G_R.cartanSet());
  standardrepk::KHatComputations khc(G_R,kgb);

  standardrepk::SR_rewrites::combination sum;

  size_t max_l=kgb.length(kgb.size()-1);

  for (kgb::KGBElt x=0; x<kgb.size(); ++x)
  {
    standardrepk::StandardRepK sr=khc.std_rep(rd.twoRho(),kgb.titsElt(x));
    standardrepk::SR_rewrites::combination c=khc.standardize(sr);
    if ((max_l-kgb.length(x))%2 == 0)
      sum += c;
    else
      sum-=c;
  }

  std::cout << "Standard normal final limit representations:\n";
  for (standardrepk::SR_rewrites::seq_no i=0; i<khc.nr_reps(); ++i)
  {
    const standardrepk::StandardRepK& sr=khc.rep_no(i);
    khc.print(std::cout << 'R' << i << ": ",sr) << std::endl;
  }

  std::cout << "Standardized expression:\n";
  {
    std::ostringstream s;
    for (standardrepk::SR_rewrites::combination::const_iterator
	   it=sum.begin(); it!=sum.end(); ++it)
    {
      s << (it->second>0 ? it==sum.begin() ? "" : " + " : " - ");
      long int ac=intutils::abs<long int>(it->second);
      if (ac!=1)
	s << ac << '*';
      s << 'R' << it->first;
    }
    ioutils::foldLine(std::cout,s.str(),"+-","",1) << std::endl;
  }

}



// Block mode functions



// Empty mode functions


/*
  Function invoked by the "test" command.
*/
void test_f()
{
  // put your code here, and define testMode at top of file appropriately

  try
  {
    realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();

    G_R.fillCartan(); // must not forget this!
    kgb::KGB kgb(G_R,G_R.cartanSet());

//     kgb_io::var_print_KGB(std::cout,mainmode::currentComplexGroup(),kgb);

    unsigned long x;
    latticetypes::Weight lambda;
    interactive::getInteractive(x,"Choose KGB element: ",kgb.size());
    prettyprint::printVector(std::cout<<"2rho = ",G_R.rootDatum().twoRho())
      << std::endl;
    interactive::getInteractive(lambda,"Give lambda-rho: ",G_R.rank());
    standardrepk::KHatComputations khc(G_R,kgb);
    khc.go(x,lambda);
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
  }
  catch (error::InputError& e) {
    e("aborted");
  }
}

} // namespace

} // namespace atlas
