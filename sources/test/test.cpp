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
#include "smithnormal.h"

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
  void Ktypeform_f();
  void Ktypemat_f();
  void mod_lattice_f();
  void branch_f();
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
  mode.add("Ktypeform",Ktypeform_f);
  mode.add("Ktypemat",Ktypemat_f);
  mode.add("mod_lattice",mod_lattice_f);
  mode.add("branch",branch_f);

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
  mode.add("Ktypeform",helpmode::nohelp_h);
  mode.add("Ktypemat",helpmode::nohelp_h);
  mode.add("mod_lattice",helpmode::nohelp_h);
  mode.add("branch",helpmode::nohelp_h);


  // add additional command tags here :

  insertTag(t,"KGB",test_tag);
  insertTag(t,"sub_KGB",test_tag);
  insertTag(t,"trivial",test_tag);
  insertTag(t,"Ktypeform",test_tag);
  insertTag(t,"Ktypemat",test_tag);
  insertTag(t,"mod_lattice",test_tag);
  insertTag(t,"branch",test_tag);

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

  ioutils::OutputFile f;

  kgb::KGB kgb(G_R,G_R.cartanSet());
  kgb_io::var_print_KGB(f,mainmode::currentComplexGroup(),kgb);
}

void sub_KGB_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();
  G_R.fillCartan(); // must not forget this!
  kgb::KGB kgb(G_R,G_R.cartanSet());
  standardrepk::KhatContext khc(G_R,kgb);

  standardrepk::StandardRepK sr=interactive::get_standardrep(khc);

  weyl::WeylWord ww;
  standardrepk::PSalgebra q= khc.theta_stable_parabolic(sr,ww);
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
  standardrepk::KhatContext khc(G_R,kgb);

  kgb::KGBElt last=kgb.size()-1;

  weyl::WeylWord ww;
  standardrepk::PSalgebra q=
    khc.theta_stable_parabolic(khc.KGB_elt_rep(last),ww);

  kgb::KGBEltList subset=khc.sub_KGB(q);
  size_t max_l=kgb.length(subset.back());

  standardrepk::combination sum(khc.height_order());
  for (size_t i=0; i<subset.size(); ++i)
  {
    kgb::KGBElt x=subset[i];
    standardrepk::StandardRepK sr=khc.std_rep(rd.twoRho(),kgb.titsElt(x));
    standardrepk::combination c=khc.standardize(sr);
    if ((max_l-kgb.length(x))%2 == 0)
      sum += c;
    else
      sum-=c;
  }


  {
    std::ostringstream s; khc.print(s,sum);
    ioutils::foldLine(std::cout,s.str(),"+\n- ","",1) << std::endl;
  }

} // trivial

void Ktypeform_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();
  G_R.fillCartan(); // must not forget this!

  kgb::KGB kgb(G_R,G_R.cartanSet());
  standardrepk::KhatContext khc(G_R,kgb);

  standardrepk::StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G_R.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  ioutils::OutputFile f;

  standardrepk::CharForm kf= khc.K_type_formula(sr);

  khc.print(f << "K-type formula for mu(",kf.first) << "):\n";
  {
    std::ostringstream s; khc.print(s,kf.second);
    ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
  }

  for (size_t i=0; i<khc.nr_reps(); ++i)
    khc.print(f << 'R' << i << ": ",khc.rep_no(i))
      << ", height: " << khc.height(i) << std::endl;


  standardrepk::combination sum(khc.height_order());

  for (standardrepk::Char::const_iterator
	 it=kf.second.begin(); it!=kf.second.end(); ++it)
  {
#ifdef VERBOSE
    khc.print(f,it->first) << " has height " << khc.height(it->first)
				   << std::endl;
    size_t old_size=khc.nr_reps();
#endif
    standardrepk::combination st=khc.standardize(it->first);
#ifdef VERBOSE
    for (size_t i=old_size; i<khc.nr_reps(); ++i)
      khc.print(f << 'R' << i << ": ",khc.rep_no(i))
        << ", height: " << khc.height(i) << std::endl;

    std::ostringstream s; khc.print(s,it->first) << " = ";
    khc.print(s,st,true);
    ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
#endif
    sum.add_multiple(st,it->second);
  }

  f << "Converted to Standard normal final limit form:\n";
  {
    std::ostringstream s; khc.print(s,sum);
    ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
  }
} // |Kypeform_f|

void Ktypemat_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();
  G_R.fillCartan(); // must not forget this!

  kgb::KGB kgb(G_R,G_R.cartanSet());
  standardrepk::KhatContext khc(G_R,kgb);

  standardrepk::StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isStandard(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not standard, as witnessed by coroot "
	<< G_R.rootDatum().coroot(khc.fiber(sr).simpleImaginary(witness))
	<< ".\n";
      return;
    }
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G_R.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  standardrepk::combination c=khc.standardize(sr);

  if (c.empty())
  {
    khc.print(std::cout << "Representation ",sr) << " is zero.\n";
    return;
  }

  assert(c.size()==1 and khc.rep_no(c.begin()->first)==sr);

  khc.print(std::cout << "Height of representation ",sr) << " is "
    << khc.height(c.begin()->first) << ".\n";
  unsigned long bound=
    interactive::get_bounded_int(interactive::common_input(),
				 "Give height bound: ",
				 9999);

  ioutils::OutputFile f;

  std::set<standardrepk::equation> singleton;
  singleton.insert(khc.mu_equation(c.begin()->first,bound));

  {
    standardrepk::equation init=*singleton.begin();
    khc.print(f << "Initial formula: mu(",khc.rep_no(init.first))
      << ") =\n";
    {
      std::ostringstream s; khc.print(s,init.second);
      ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
    }
  }

  std::vector<standardrepk::equation> system =
    khc.saturate(singleton,bound);

  std::vector<standardrepk::seq_no> new_order;
  matrix::Matrix<standardrepk::CharCoeff> m =
    standardrepk::triangularize(system,new_order);

  f << "Ordering of representations/K-types:\n";
  for (std::vector<standardrepk::seq_no>::const_iterator
	 it=new_order.begin(); it!=new_order.end(); ++it)
    khc.print(f,khc.rep_no(*it)) << ", height " << khc.height(*it)
      << std::endl;

#ifdef VERBOSE
  prettyprint::printMatrix(f<<"Triangular system:\n",m,3);
#endif

  prettyprint::printMatrix(f<<"Matrix of K-type multiplicites:\n",
			   standardrepk::inverse_lower_triangular(m),
			   3);

} // |Ktypemat_f|

void mod_lattice_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();
  G_R.fillCartan(); // must not forget this!

  unsigned long cn=interactive::get_Cartan_class(G_R.cartanSet());

  latticetypes::LatticeMatrix q = G_R.cartan(cn).involution();
  for (size_t j = 0; j<q.numRows(); ++j)
    q(j,j) -= 1;

  latticetypes::WeightList bs; matrix::initBasis(bs,q.numRows());
  latticetypes::CoeffList invf; smithnormal::smithNormal(invf,bs.begin(),q);

   size_t l = invf.size();
   size_t f=0; while (f<l and invf[f]==1) ++f; // skip invariant factors 1

   std::cout << "At Cartan class " << cn;
   if (l==0)
     std::cout << " weights are used unchanged";
   else
   {
     std::cout << " weights are modulo";
     if (f>0)
     {
       std::cout << " multiples of ";
       basic_io::seqPrint(std::cout,&bs[0],&bs[f],", ", "", "");
       if (f<l)
	 std::cout << " and";
     }
     if (f<l)
     {
       std::cout << " even multiples of ";
       basic_io::seqPrint(std::cout,&bs[f],&bs[l],", ", "", "");
     }
   }
   std::cout << ".\n";

} // |mod_lattice_f|

void branch_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();
  G_R.fillCartan(); // must not forget this!

  kgb::KGB kgb(G_R,G_R.cartanSet());
  standardrepk::KhatContext khc(G_R,kgb);

  standardrepk::StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isStandard(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not standard, as witnessed by coroot "
	<< G_R.rootDatum().coroot(khc.fiber(sr).simpleImaginary(witness))
	<< ".\n";
      return;
    }
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G_R.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  standardrepk::combination c=khc.standardize(sr);

  if (c.empty())
  {
    khc.print(std::cout << "Representation ",sr) << " is zero.\n";
    return;
  }

  assert(c.size()==1 and khc.rep_no(c.begin()->first)==sr);

  khc.print(std::cout << "Height of representation ",sr) << " is "
    << khc.height(c.begin()->first) << ".\n";
  unsigned long bound=
    interactive::get_bounded_int(interactive::common_input(),
				 "Give height bound: ",
				 9999);

  ioutils::OutputFile f;

  standardrepk::combination result=khc.branch(c.begin()->first,bound);

  {
    std::ostringstream s; khc.print(s,result);
    ioutils::foldLine(f,s.str(),"+","",1) << std::endl;
  }

} // |branch_f|

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
    const complexredgp::ComplexReductiveGroup& G=G_R.complexGroup();

    G_R.fillCartan(); // must not forget this!
    kgb::KGB kgb(G_R,G_R.cartanSet());

    unsigned long x=interactive::get_bounded_int
      (interactive::sr_input(),"Choose KGB element: ",kgb.size());

    prettyprint::printVector(std::cout<<"2rho = ",G_R.rootDatum().twoRho())
      << std::endl;

    latticetypes::Weight lambda=
      interactive::get_weight(interactive::sr_input(),
			      "Give lambda-rho: ",
			      G_R.rank());
    standardrepk::KhatContext khc(G_R,kgb);

    standardrepk::StandardRepK sr=khc.std_rep_rho_plus(lambda,kgb.titsElt(x));

    (lambda *= 2) += G.rootDatum().twoRho();
    prettyprint::printVector(std::cout << "Weight (1/2)",lambda);
    prettyprint::printVector(std::cout << " converted to (1/2)",khc.lift(sr));

    const weyl::TwistedInvolution& canonical=G.twistedInvolution(sr.Cartan());
    if (kgb.involution(x)!=canonical)
      prettyprint::printWeylElt(std::cout << " at involution ",
				canonical, G.weylGroup());
    std::cout << "\nHeight is " << khc.height(sr) << std::endl;

    khc.go(sr);
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
