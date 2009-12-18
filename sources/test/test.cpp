/*
  This is test.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include "io.h"    //needed for help commands

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
#include "prerootdata.h"
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
#include "repr.h"
#include "free_abelian.h"
#include "matreduc.h"
#include "testrun.h"
#include "basic_io.h"

/*****************************************************************************

  This module contains some commands for testing the program.

******************************************************************************/

namespace atlas {

namespace {
  using namespace test;

  // functions for the test commands

  void test_f();

  void roots_rootbasis_f();
  void coroots_rootbasis_f();
  void posroots_rootbasis_f();
  void poscoroots_rootbasis_f();
  void checkbasept_f();
  void sub_KGB_f();
  void trivial_f();
  void Ktypeform_f();
  void qKtypeform_f();
  void Ktypemat_f();
  void qKtypemat_f();
  void mod_lattice_f();
  void branch_f();
  void qbranch_f();
  void srtest_f();
  void testrun_f();
  void exam_f();

  void X_f();


  // help functions

  // no |test_h|, test command remains development-reserved; uses |nohelp_h|

  void roots_rootbasis_h();
  void coroots_rootbasis_h();
  void posroots_rootbasis_h();
  void poscoroots_rootbasis_h();
  void mod_lattice_h();
  void checkbasept_h();
  void sub_KGB_h();
  void trivial_h();
  void Ktypeform_h();
  void Ktypemat_h();
  void branch_h();
  void srtest_h();
  void testrun_h();
  void exam_h();

  // tags
  const char* test_tag = "test command (for development only)";

  const char* roots_rootbasis_tag = "outputs the roots in the simple rootbasis";
  const char* coroots_rootbasis_tag =
    "outputs the coroots in the simple coroot basis";
  const char* posroots_rootbasis_tag =
     "outputs the positive roots in the simple root basis";
  const char* poscoroots_rootbasis_tag =
      "outputs the positive coroots in the simple coroot basis";
  const char* sub_KGB_tag =
    "computes subset of KGB data used in Ktypeform";
  const char* trivial_tag =
    "tries to compute Ktypeform for trivial representation";
  const char* Ktypeform_tag = "computes formula for a K-type";
  const char* Ktypemat_tag =
     "computes matrix relating K-types and standard modules";
  const char* mod_lattice_tag =
    "gives a basis of quotient of character lattice";
  const char* branch_tag = "computes restriction of representation to K";
  const char* srtest_tag = "gives information about a representation";
  const char* testrun_tag =
                "iterates over root data of given rank, calling examine";
  const char* examine_tag =
    "tests if kgb and KGB commands generate identical numberings";

/*
  For convenience, the "test" command is added to the mode that is flagged by
  the testMode constant defined here; therefore "test" appears (conditionally)
  in every template instance of |addTestCommands| and of |addTestHelp| below.
  Set this constant according to the requirements of the |test_f| function.
*/
  enum TestMode {EmptyMode, MainMode, RealMode, BlockMode, numTestMode};
  const TestMode testMode = RealMode; // currently does representation test

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
  mode.add("testrun",testrun_f);
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

  mode.add("X",X_f);
}


// Add to the real mode the test commands that require that mode.
template<>
void addTestCommands<realmode::RealmodeTag>
  (commands::CommandMode& mode, realmode::RealmodeTag)
{
  if (testMode == RealMode)
    mode.add("test",test_f);

  // add additional commands here :

  mode.add("checkbasept",checkbasept_f);
  mode.add("sub_KGB",sub_KGB_f);
  mode.add("trivial",trivial_f);
  mode.add("Ktypeform",Ktypeform_f);
  mode.add("qKtypeform",qKtypeform_f);
  mode.add("Ktypemat",Ktypemat_f);
  mode.add("qKtypemat",qKtypemat_f);
  mode.add("mod_lattice",mod_lattice_f);
  mode.add("branch",branch_f);
  mode.add("qbranch",qbranch_f);
  mode.add("srtest",srtest_f);

  mode.add("examine",exam_f);
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
    mode.add("test",nohelp_h);
    commands::insertTag(t,"test",test_tag);
  }


  // add additional help commands here:
  mode.add("testrun",helpmode::nohelp_h);

  // add additional command tags here:
  commands::insertTag(t,"testrun",testrun_tag);

}


// Add to the main mode the help commands for test commands with that mode
template<> void addTestHelp<mainmode::MainmodeTag>
             (commands::CommandMode& mode, commands::TagDict& t,
	      mainmode::MainmodeTag)
{
  using namespace commands;
  using namespace helpmode;

  if (testMode == MainMode) {
    mode.add("test",nohelp_h);
    insertTag(t,"test",test_tag);
  }

  // add additional help commands here:

  mode.add("roots_rootbasis",roots_rootbasis_h);
  mode.add("posroots_rootbasis",posroots_rootbasis_h);
  mode.add("coroots_rootbasis",coroots_rootbasis_h);
  mode.add("poscoroots_rootbasis",poscoroots_rootbasis_h);

  // add additional command tags here :

  insertTag(t,"roots_rootbasis",roots_rootbasis_tag);
  insertTag(t,"posroots_rootbasis",posroots_rootbasis_tag);
  insertTag(t,"coroots_rootbasis",coroots_rootbasis_tag);
  insertTag(t,"poscoroots_rootbasis",poscoroots_rootbasis_tag);

}


// Add to the real mode the help commands for test commands with that mode
template<> void addTestHelp<realmode::RealmodeTag>
             (commands::CommandMode& mode, commands::TagDict& t,
	      realmode::RealmodeTag)
{
  using namespace commands;
  using namespace helpmode;

  if (testMode == RealMode) {
    mode.add("test",nohelp_h);
    insertTag(t,"test",test_tag);
  }

  // add additional help commands here:

  mode.add("checkbasept",checkbasept_h);
  mode.add("sub_KGB",helpmode::nohelp_h);
  mode.add("trivial",helpmode::nohelp_h);
  mode.add("Ktypeform",Ktypeform_h);
  mode.add("Ktypemat",Ktypemat_h);
  mode.add("mod_lattice",mod_lattice_h);
  mode.add("branch",branch_h);
  mode.add("srtest",srtest_h);
  mode.add("examine",helpmode::nohelp_h);


  // add additional command tags here :

  insertTag(t,"sub_KGB",sub_KGB_tag);
  insertTag(t,"trivial",trivial_tag);
  insertTag(t,"Ktypeform",Ktypeform_tag);
  insertTag(t,"Ktypemat",Ktypemat_tag);
  insertTag(t,"mod_lattice",mod_lattice_tag);
  insertTag(t,"branch",branch_tag);
  insertTag(t,"srtest",srtest_tag);
  insertTag(t,"examine",examine_tag);

}

// Add to the block mode the help commands for test commands with that mode
template<> void addTestHelp<blockmode::BlockmodeTag>
             (commands::CommandMode& mode, commands::TagDict& t,
	      blockmode::BlockmodeTag)
{
  using namespace commands;
  using namespace helpmode;

  if (testMode == BlockMode) {
    mode.add("test",nohelp_h);
    insertTag(t,"test",test_tag);
  }

  // add additional help commands here:



  // add additional command tags here:

}

} // namespace test


/*****************************************************************************

        Chapter II -- Functions for the test commands

******************************************************************************/

namespace {

  // Main mode functions

// Print the roots in the simple root coordinates.
void roots_rootbasis_f()
{
  try {
    const rootdata::RootSystem& rs =
      mainmode::currentComplexGroup().rootSystem();
    ioutils::OutputFile file;

    for (rootdata::RootNbr i=0; i<rs.numRoots(); ++i)
      prettyprint::printInRootBasis(file,i,rs) << std::endl;
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

// Print the positive roots in the simple root coordinates.
void posroots_rootbasis_f()

{
  try {
    const rootdata::RootSystem& rs =
      mainmode::currentComplexGroup().rootSystem();
    ioutils::OutputFile file;

    prettyprint::printInRootBasis(file,rs.posRootSet(),rs);
  }
  catch (error::InputError& e) {
    e("aborted");
  }

}

// Print the coroots in the simple coroot coordinates.
void coroots_rootbasis_f()
{
  try
  {
    const rootdata::RootSystem rs
      (mainmode::currentComplexGroup().dualRootSystem());
    ioutils::OutputFile file;

    for (rootdata::RootNbr i=0; i<rs.numRoots(); ++i)
      prettyprint::printInRootBasis(file,i,rs) << std::endl;
  }
  catch (error::InputError& e)
  {
    e("aborted");
  }

}

// Print the positive coroots in the simple coroot coordinates.
void poscoroots_rootbasis_f()
{
  try
  {
    const rootdata::RootSystem rs
      (mainmode::currentComplexGroup().dualRootSystem());
    ioutils::OutputFile file;

    prettyprint::printInRootBasis(file,rs.posRootSet(),rs);
  }
  catch (error::InputError& e)
  {
    e("aborted");
  }

}

// Real mode functions

void checkbasept_f()
{
  realredgp::RealReductiveGroup& G_R = realmode::currentRealGroup();

  kgb::KGB kgb(G_R,G_R.Cartan_set());
  kltest::checkBasePoint(kgb);
}

void sub_KGB_f()
{
  realredgp::RealReductiveGroup& G = realmode::currentRealGroup();
  standardrepk::KhatContext khc(G);

  standardrepk::StandardRepK sr=interactive::get_standardrep(khc);

  weyl::WeylWord ww;
  standardrepk::PSalgebra p= khc.theta_stable_parabolic(sr,ww);
  kgb::KGBEltList sub=khc.sub_KGB(p);

  std::cout << "Conjugating word [" << ww << "]\n";
  kgb_io::print_sub_KGB(std::cout,G.kgb(),sub);
}

void trivial_f()
{
  realredgp::RealReductiveGroup& G = realmode::currentRealGroup();
  const rootdata::RootDatum& rd=G.rootDatum();
  const kgb::KGB& kgb = G.kgb();

  standardrepk::KhatContext khc(G);

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
  realredgp::RealReductiveGroup& G = realmode::currentRealGroup();

  standardrepk::KhatContext khc(G);

  standardrepk::StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
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

void qKtypeform_f()
{
  realredgp::RealReductiveGroup& G = realmode::currentRealGroup();

  standardrepk::qKhatContext khc(G);

  standardrepk::StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  ioutils::OutputFile f;

  standardrepk::q_CharForm kf= khc.q_K_type_formula(sr);

  khc.print(f << "q-K-type formula for mu(",kf.first) << "):\n";
  {
    std::ostringstream s; khc.print(s,kf.second);
    ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
  }

  standardrepk::q_combin sum(khc.height_order());

  for (standardrepk::q_Char::const_iterator
	 it=kf.second.begin(); it!=kf.second.end(); ++it)
  {
#ifdef VERBOSE
    khc.print(f,it->first) << " has height " << khc.height(it->first)
				   << std::endl;
    size_t old_size=khc.nr_reps();
#endif
    standardrepk::q_combin st=khc.standardize(it->first);
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
} // |qKypeform_f|

void Ktypemat_f()
{
  realredgp::RealReductiveGroup& G = realmode::currentRealGroup();

  standardrepk::KhatContext khc(G);

  standardrepk::StandardRepK sr=interactive::get_standardrep(khc);
  khc.normalize(sr);

  {
    size_t witness;
    if (not khc.isStandard(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not standard, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleImaginary(witness))
	<< ".\n";
      return;
    }
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  standardrepk::combination c=khc.standardize(sr);

  if (c.empty())
  {
    khc.print(std::cout << "Representation ",sr) << " is zero.\n";
    return;
  }

  assert(c.size()==1);
  assert(khc.rep_no(c.begin()->first)==sr);

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
  std::vector<standardrepk::seq_no> new_order;

#ifdef VERBOSE
  matrix::Matrix_base<standardrepk::CharCoeff> m;
  matrix::Matrix_base<standardrepk::CharCoeff> ktypemat =
    khc.K_type_matrix(singleton,bound,new_order,&m);
#else
  matrix::Matrix_base<standardrepk::CharCoeff> ktypemat =
    khc.K_type_matrix(singleton,bound,new_order,NULL);
#endif

  std::cout << "Ordering of representations/K-types:\n";
  for (std::vector<standardrepk::seq_no>::const_iterator
	 it=new_order.begin(); it!=new_order.end(); ++it)
    khc.print(std::cout,khc.rep_no(*it)) << ", height " << khc.height(*it)
       << std::endl;

#ifdef VERBOSE
  prettyprint::printMatrix(std::cout<<"Triangular system:\n",m,3);
#endif

  prettyprint::printMatrix(f<<"Matrix of K-type multiplicites:\n",ktypemat,3);

} // |Ktypemat_f|

void qKtypemat_f()
{
  realredgp::RealReductiveGroup& G = realmode::currentRealGroup();

  standardrepk::qKhatContext khc(G);

  standardrepk::StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isNormal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not normal, as witnessed by coroot sum "
	<< khc.info(sr.Cartan()).coroot_sum(witness)
	<< ".\n";
      return;
    }
    if (not khc.isStandard(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not standard, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleImaginary(witness))
	<< ".\n";
      return;
    }
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  standardrepk::q_combin c=khc.standardize(sr);

  if (c.empty())
  {
    khc.print(std::cout << "Representation ",sr) << " is zero.\n";
    return;
  }

  assert(c.size()==1);
  assert(khc.rep_no(c.begin()->first)==sr);

  khc.print(std::cout << "Height of representation ",sr) << " is "
    << khc.height(c.begin()->first) << ".\n";
  unsigned long bound=
    interactive::get_bounded_int(interactive::common_input(),
				 "Give height bound: ",
				 9999);

  ioutils::OutputFile f;

  std::set<standardrepk::q_equation> singleton;
  singleton.insert(khc.mu_equation(c.begin()->first,bound));

  {
    standardrepk::q_equation init=*singleton.begin();
    khc.print(f << "Initial formula: mu(",khc.rep_no(init.first))
      << ") =\n";
    {
      std::ostringstream s; khc.print(s,init.second);
      ioutils::foldLine(f,s.str(),"+\n- ","",1) << std::endl;
    }
  }

#ifdef VERBOSE

  std::vector<standardrepk::q_equation> system =
    khc.saturate(singleton,bound);
  std::cout << "System of equations:\n";
  for (size_t i=0; i<system.size(); ++i)
  {
    const standardrepk::q_equation& si=system[i];
    khc.print(std::cout<< si.first << ' ',khc.rep_no(si.first))
      << " [" << khc.height(si.first) << "]\n     ";
    for (standardrepk::q_combin::const_iterator
	   it=si.second.begin(); it!=si.second.end(); ++it)
      std::cout << '+' << it->second << "*I(" << it->first << ')';
    std::cout << std::endl;
  }
#endif

  std::vector<standardrepk::seq_no> new_order;

#ifdef VERBOSE
  matrix::Matrix_base<standardrepk::q_CharCoeff> m;
  matrix::Matrix_base<standardrepk::q_CharCoeff> ktypemat =
    khc.K_type_matrix(singleton,bound,new_order,&m);
#else
  matrix::Matrix_base<standardrepk::q_CharCoeff> ktypemat =
    khc.K_type_matrix(singleton,bound,new_order,NULL);
#endif

  f << "Ordering of representations/K-types:\n";
  for (std::vector<standardrepk::seq_no>::const_iterator
	 it=new_order.begin(); it!=new_order.end(); ++it)
    khc.print(f,khc.rep_no(*it)) << ", height " << khc.height(*it)
      << std::endl;

#ifdef VERBOSE
  prettyprint::printMatrix(std::cout<<"Triangular system:\n",m,3);
#endif

  prettyprint::printMatrix(f<<"Matrix of K-type multiplicites:\n",
			   standardrepk::inverse_lower_triangular(m),
			   3);

} // |qKtypemat_f|

void mod_lattice_f()
{
  realredgp::RealReductiveGroup& G = realmode::currentRealGroup();

  unsigned long cn=interactive::get_Cartan_class(G.Cartan_set());

  latticetypes::LatticeMatrix q = G.cartan(cn).involution();
  for (size_t j = 0; j<q.numRows(); ++j)
    q(j,j) -= 1;

  latticetypes::CoeffList factor;
  latticetypes::LatticeMatrix b = matreduc::adapted_basis(q,factor);

  bitset::RankFlags units, doubles;
  unsigned n1=0,n2=0;

  for (size_t i=0; i<factor.size(); ++i)
    if (factor[i]==1)
      units.set(i),++n1;
    else if (factor[i]==2)
      doubles.set(i),++n2;

   std::cout << "At Cartan class " << cn;
   if (n1+n2==0)
     std::cout << " weights are used unchanged";
   else
   {
     std::cout << " weights are modulo";
     if (n1>0)
     {
       std::cout << " multiples of ";
       for (bitset::RankFlags::iterator it=units.begin(); it(); ++it,--n1)
	 std::cout << b.column(*it) << (n1>2 ? ", " : n1>1 ? ", and " : "");
       if (n2>0)
	 std::cout << " and";
     }
     if (n2>0)
     {
       std::cout << " even multiples of ";
       for (bitset::RankFlags::iterator it=doubles.begin(); it(); ++it,--n2)
	 std::cout << b.column(*it) << (n2>2 ? ", " : n2>1 ? ", and " : "");
     }
   }
   std::cout << ".\n";

} // |mod_lattice_f|

void branch_f()
{
  realredgp::RealReductiveGroup& G = realmode::currentRealGroup();

  standardrepk::KhatContext khc(G);

  standardrepk::StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isStandard(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not standard, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleImaginary(witness))
	<< ".\n";
      return;
    }
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
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

void qbranch_f()
{
  realredgp::RealReductiveGroup& G = realmode::currentRealGroup();

  standardrepk::qKhatContext khc(G);

  standardrepk::StandardRepK sr=interactive::get_standardrep(khc);

  {
    size_t witness;
    if (not khc.isStandard(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not standard, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleImaginary(witness))
	<< ".\n";
      return;
    }
    if (not khc.isFinal(sr,witness))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.fiber(sr).simpleReal(witness)) << ".\n";
      return;
    }
  }

  standardrepk::q_combin c=khc.standardize(sr);

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

  standardrepk::q_combin result=khc.branch(c.begin()->first,bound);

  {
    std::ostringstream s; khc.print(s,result);
    ioutils::foldLine(f,s.str(),"+","",1) << std::endl;
  }

} // |qbranch_f|


// Block mode functions



// Empty mode functions


/*
  Function invoked by the "srtest" command.
*/
void srtest_f()
{
  // put your code here, and define testMode at top of file appropriately

  try
  {
    realredgp::RealReductiveGroup& G = realmode::currentRealGroup();
    const kgb::KGB& kgb = G.kgb();

    unsigned long x=interactive::get_bounded_int
      (interactive::sr_input(),"Choose KGB element: ",G.kgb().size());

    prettyprint::printVector(std::cout<<"2rho = ",G.rootDatum().twoRho())
      << std::endl;

    latticetypes::Weight lambda=
      interactive::get_weight(interactive::sr_input(),
			      "Give lambda-rho: ",
			      G.rank());
    standardrepk::KhatContext khc(G);

    standardrepk::StandardRepK sr=khc.std_rep_rho_plus(lambda,kgb.titsElt(x));

    (lambda *= 2) += G.rootDatum().twoRho();
    prettyprint::printVector(std::cout << "Weight (1/2)",lambda);
    prettyprint::printVector(std::cout << " converted to (1/2)",khc.lift(sr));

    const weyl::TwistedInvolution& canonical =
      G.complexGroup().twistedInvolution(sr.Cartan());
    if (kgb.involution(x)!=canonical)
      prettyprint::printWeylElt(std::cout << " at involution ",
				canonical, G.weylGroup());
    std::cout << "\nHeight is " << khc.height(sr) << std::endl;

    khc.go(sr);
  }
  catch (error::MemoryOverflow& e)
  {
    e("error: memory overflow");
  }
  catch (error::InputError& e)
  {
    e("aborted");
  }
}

bool examine(realredgp::RealReductiveGroup& G)
{
  kgb::KGB kgb1(G);
  kgb::KGB kgb2(G,G.Cartan_set());
  if (kgb1.size()!=kgb2.size()) return false;
  for (size_t i=0; i<kgb1.size(); ++i)
  {
    for (size_t s=0; s<G.semisimpleRank(); ++s)
    {
      if (kgb1.cross(s,i)!=kgb2.cross(s,i)) return false;
      if (kgb1.cayley(s,i)!=kgb2.cayley(s,i)) return false;
    }
    if(kgb1.status(i)!=kgb2.status(i)) return false;
  }
  return true;
}

void testrun_f()
{
  unsigned long rank=interactive::get_bounded_int
    (interactive::common_input(),"rank: ",constants::RANK_MAX+1);
  std::cout << "Testing agreement of kgb and KGB:\n";
  for (testrun::LieTypeIterator it(testrun::Semisimple,rank); it(); ++it)
  {
    std::cout<< *it << std::endl;
    size_t count=0;
    for (testrun::CoveringIterator cit(*it); cit(); ++cit)
    {
      if (count>0) std::cout << ',' << std::flush;
      rootdata::RootDatum rd(*cit);
      latticetypes::LatticeMatrix id(rd.rank()); // identity
      complexredgp::ComplexReductiveGroup G(rd,id);
      for (realform::RealForm rf=0; rf<G.numRealForms(); ++rf)
      {
	realredgp::RealReductiveGroup G_R(G,rf);
	if (not examine(G_R))
	  std::cout << "Failure at real form " << rf << std::endl;
      }
      std::cout << ++count;
    }
    std::cout << '.' << std::endl;
  }

}

void exam_f()
{
  std::cout << "kgb and KGB "
            << (examine(realmode::currentRealGroup()) ? "agree" : "differ")
	    << std::endl;
}

void X_f()
{
  complexredgp::ComplexReductiveGroup& G=mainmode::currentComplexGroup();
  tits::GlobalTitsGroup Tg(G);
  kgb::global_KGB kgb(G,Tg);
  ioutils::OutputFile f;
  kgb_io::print_X(f,kgb);
}

void test_f()
{
  try
  {
    realredgp::RealReductiveGroup& G = realmode::currentRealGroup();

    standardrepk::KhatContext khc(G);
    standardrepk::StandardRepK srk=interactive::get_standardrep(khc);

    latticetypes::RatWeight nu =
      interactive::get_ratweight(interactive::sr_input(),"nu: ",G.rank());

    repr::Rep_context rc(G);
    repr::StandardRepr sr = rc.sr(srk,khc,nu);
    tits::GlobalTitsElement y=rc.y(sr);
    complexredgp::ComplexReductiveGroup& dual_G=mainmode::current_dual_group();
    kgb::global_KGB dual_KGB (dual_G,tits::GlobalTitsGroup(dual_G),y);
    ioutils::OutputFile f;
    f << "y is element " << dual_KGB.lookup(y) << " in dual KGB set:\n";
    kgb_io::print_X(f,dual_KGB);
  }
  catch (error::MemoryOverflow& e)
  {
    e("error: memory overflow");
  }
  catch (error::InputError& e)
  {
    e("aborted");
  }
}

//Help commands

void roots_rootbasis_h()
{
  io::printFile(std::cerr,"roots_rootbasis.help",io::MESSAGE_DIR);
}

void coroots_rootbasis_h()
{
  io::printFile(std::cerr,"coroots_rootbasis.help",io::MESSAGE_DIR);
}

void posroots_rootbasis_h()
{
  io::printFile(std::cerr,"posroots_rootbasis.help",io::MESSAGE_DIR);
}

void poscoroots_rootbasis_h()
{
  io::printFile(std::cerr,"poscoroots_rootbasis.help",io::MESSAGE_DIR);
}

void checkbasept_h()
{
  io::printFile(std::cerr,"checkbasept.help",io::MESSAGE_DIR);
}

void sub_KGB_h()
{
  io::printFile(std::cerr,"sub_KGB.help",io::MESSAGE_DIR);
}

void trivial_h()
{
  io::printFile(std::cerr,"trivial.help",io::MESSAGE_DIR);
}

void Ktypemat_h()
{
  io::printFile(std::cerr,"Ktypemat.help",io::MESSAGE_DIR);
}

void Ktypeform_h()
{
  io::printFile(std::cerr,"Ktypeform.help",io::MESSAGE_DIR);
}

void mod_lattice_h()
{
  io::printFile(std::cerr,"mod_lattice.help",io::MESSAGE_DIR);
}

void branch_h()
{
  io::printFile(std::cerr,"branch.help",io::MESSAGE_DIR);
}

void srtest_h()
{
  io::printFile(std::cerr,"srtest.help",io::MESSAGE_DIR);
}

} // |namespace|

} // |namespace atlas|
