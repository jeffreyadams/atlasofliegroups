/*
  This is test.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "io.h"    //needed for help commands

#include "test.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <map>

#include "atlas_types.h" // here to preempt double inclusion of _fwd files

#include "free_abelian.h"
#include "permutations.h"
#include "matreduc.h"

#include "dynkin.h"
#include "prerootdata.h"
#include "rootdata.h"
#include "cartanclass.h"
#include "complexredgp.h"
#include "realredgp.h"

#include "kgb.h"
#include "blocks.h"
#include "klsupport.h"
#include "kl.h"
#include "standardrepk.h"
#include "repr.h"
#include "ext_block.h"

#include "ioutils.h"
#include "basic_io.h"
#include "prettyprint.h"
#include "interactive.h"
#include "realform_io.h"
#include "kgb_io.h"
#include "block_io.h"

#include "commands.h"
#include "helpmode.h"
#include "emptymode.h"
#include "mainmode.h"
#include "realmode.h"
#include "blockmode.h"
#include "reprmode.h"

#include "testrun.h"
#include "kltest.h"


namespace atlas {

namespace test {

/*****************************************************************************

  This module contains some commands for testing the program.

******************************************************************************/

namespace {
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
  void embedding_f();

/*
  For convenience, the "test" command is added to the mode that is flagged by
  the testMode constant defined here; therefore "test" appears (conditionally)
  in every template instance of |addTestCommands| and of |addTestHelp| below.
  Set this constant according to the requirements of the |test_f| function.
*/
  enum TestMode {EmptyMode, MainMode, RealMode, BlockMode, ReprMode,
		 numTestMode};
  const TestMode testMode = ReprMode; // currently does test of extended block

  // utilities
  const RootDatum& currentRootDatum();

} // |namespace|

/*****************************************************************************

        Chapter I -- Functions declared by test.h

  This section defines the functions declared in test.h :

    - addTestCommands() : adds the test commands to the main command
      tree;
    - addTestHelp() : adds help functionality;

******************************************************************************/

  const char* test_tag = "test command (for development only)";


/* |addTestCommands| is a template only so that its declaration is shorter;
   its defining instances are defined just as overloaded functions would,
   since each of them needs to test a specific value of |testMode|.
*/


// Add to the empty mode node the test commands that require that mode.
template<>
void addTestCommands<commands::EmptymodeTag> (commands::CommandNode& mode)
{
  mode.add("testrun",testrun_f,
	   "iterates over root data of given rank, calling examine",
	   commands::use_tag);
  if (testMode == EmptyMode)
    mode.add("test",test_f,"test command (for development only)");
}


// Add to the main mode node the test commands that require that mode.
template<>
void addTestCommands<commands::MainmodeTag> (commands::CommandNode& mode)
{
  mode.add("roots_rootbasis",roots_rootbasis_f,
	   "outputs the roots in the simple root basis",commands::std_help);
  mode.add("posroots_rootbasis",posroots_rootbasis_f,
	   "outputs the positive roots in the simple root basis",
	   commands::std_help);
  mode.add("coroots_rootbasis",coroots_rootbasis_f,
	   "outputs the coroots in the simple coroot basis",
	   commands::std_help);
  mode.add("poscoroots_rootbasis",poscoroots_rootbasis_f,
	   "outputs the positive coroots in the simple coroot basis",
	   commands::std_help);

  mode.add("X",X_f,"prints union of K\\G/B for real forms in inner class");

  if (testMode == MainMode)
    mode.add("test",test_f,test_tag);

}

// Add to the real mode node the test commands that require that mode.
template<>
void addTestCommands<commands::RealmodeTag> (commands::CommandNode& mode)
{
  mode.add("checkbasept",checkbasept_f,
	   "checks basepoint conjecture",commands::std_help);
  mode.add("sub_KGB",sub_KGB_f,
	   "computes subset of KGB data used in Ktypeform",
	   commands::use_tag);
  mode.add("trivial",trivial_f,
	   "tries to compute Ktypeform for trivial representation",
	   commands::use_tag);
  mode.add("Ktypeform",Ktypeform_f,
	   "computes formula for a K-type",commands::std_help);
  mode.add("qKtypeform",qKtypeform_f,
	   "q version of Ktypeform",commands::use_tag);
  mode.add("Ktypemat",Ktypemat_f,
	   "computes matrix relating K-types and standard modules",
	   commands::std_help);
  mode.add("qKtypemat",qKtypemat_f,"q version of Ktypemat",commands::use_tag);
  mode.add("mod_lattice",mod_lattice_f,
	   "gives a basis of quotient of character lattice",commands::std_help);
  mode.add("branch",branch_f,
	   "computes restriction of representation to K",commands::std_help);
  mode.add("qbranch",qbranch_f,"q version of branch");
  mode.add("srtest",srtest_f,
	   "gives information about a representation",commands::std_help);

  mode.add("examine",exam_f,
	   "tests whether block is ordered by Weyl length",commands::use_tag);
  mode.add("embedding",embedding_f,
	   "give information about embedding of KGB of subsystem",
	   commands::std_help);

  if (testMode == RealMode)
    mode.add("test",test_f,test_tag);

}


// Add to the block mode the test commands that require that mode.
template<>
void addTestCommands<commands::BlockmodeTag> (commands::CommandNode& mode)
{
  if (testMode == BlockMode)
    mode.add("test",test_f,test_tag);

}

// Add to the repr mode the test commands that require that mode.
template<>
void addTestCommands<commands::ReprmodeTag> (commands::CommandNode& mode)
{
  if (testMode == ReprMode)
    mode.add("test",test_f,test_tag);

  // add additional commands here :
}


/*****************************************************************************

        Chapter II -- Functions for the test commands

******************************************************************************/

namespace {

  // Main mode functions

// Print the roots in the simple root coordinates.
void roots_rootbasis_f()
{
  const RootSystem& rs =  commands::currentComplexGroup().rootSystem();
  ioutils::OutputFile file;

  for (RootNbr i=0; i<rs.numRoots(); ++i)
    prettyprint::printInRootBasis(file,i,rs) << std::endl;
}

// Print the positive roots in the simple root coordinates.
void posroots_rootbasis_f()

{
  const RootSystem& rs = commands::currentComplexGroup().rootSystem();

  ioutils::OutputFile file;
  prettyprint::printInRootBasis(file,rs.posRootSet(),rs);
}

// Print the coroots in the simple coroot coordinates.
void coroots_rootbasis_f()
{
  const RootSystem rs (commands::currentComplexGroup().dualRootSystem());

  ioutils::OutputFile file;
  for (RootNbr i=0; i<rs.numRoots(); ++i)
    prettyprint::printInRootBasis(file,i,rs) << std::endl;

}

// Print the positive coroots in the simple coroot coordinates.
void poscoroots_rootbasis_f()
{
  const RootSystem rs (commands::currentComplexGroup().dualRootSystem());

  ioutils::OutputFile file;
  prettyprint::printInRootBasis(file,rs.posRootSet(),rs);
}

// Real mode functions

void checkbasept_f()
{
  RealReductiveGroup& G_R = commands::currentRealGroup();

  KGB kgb(G_R,G_R.Cartan_set());
  kltest::checkBasePoint(kgb);
}

void sub_KGB_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();
  standardrepk::KhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);

  WeylWord ww;
  standardrepk::PSalgebra p= khc.theta_stable_parabolic(sr,ww);
  KGBEltList sub=khc.sub_KGB(p);

  std::cout << "Conjugating word [" << ww << "]\n";
  kgb_io::print_sub_KGB(std::cout,G.kgb(),sub);
}

void trivial_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();
  const RootDatum& rd=G.rootDatum();
  const KGB& kgb = G.kgb();

  standardrepk::KhatContext khc(G);

  KGBElt last=kgb.size()-1;

  WeylWord ww;
  standardrepk::PSalgebra q=
    khc.theta_stable_parabolic(khc.KGB_elt_rep(last),ww);

  KGBEltList subset=khc.sub_KGB(q);
  size_t max_l=kgb.length(subset.back());

  standardrepk::combination sum(khc.height_order());
  for (size_t i=0; i<subset.size(); ++i)
  {
    KGBElt x=subset[i];
    StandardRepK sr=khc.std_rep(rd.twoRho(),kgb.titsElt(x));
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
  RealReductiveGroup& G = commands::currentRealGroup();

  standardrepk::KhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);

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
  RealReductiveGroup& G = commands::currentRealGroup();

  standardrepk::qKhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);

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
  RealReductiveGroup& G = commands::currentRealGroup();

  standardrepk::KhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);
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
  RealReductiveGroup& G = commands::currentRealGroup();

  standardrepk::qKhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);

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

  prettyprint::printMatrix(f<<"Matrix of K-type multiplicites:\n",ktypemat,3);
} // |qKtypemat_f|

void mod_lattice_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();

  unsigned long cn=interactive::get_Cartan_class(G.Cartan_set());

  WeightInvolution q = G.cartan(cn).involution();
  for (size_t j = 0; j<q.numRows(); ++j)
    q(j,j) -= 1;

  CoeffList factor;
  int_Matrix b = matreduc::adapted_basis(q,factor);

  RankFlags units, doubles;
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
       for (RankFlags::iterator it=units.begin(); it(); ++it,--n1)
	 std::cout << b.column(*it) << (n1>2 ? ", " : n1>1 ? ", and " : "");
       if (n2>0)
	 std::cout << " and";
     }
     if (n2>0)
     {
       std::cout << " even multiples of ";
       for (RankFlags::iterator it=doubles.begin(); it(); ++it,--n2)
	 std::cout << b.column(*it) << (n2>2 ? ", " : n2>1 ? ", and " : "");
     }
   }
   std::cout << ".\n";

} // |mod_lattice_f|

void branch_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();

  standardrepk::KhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);

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
  RealReductiveGroup& G = commands::currentRealGroup();

  standardrepk::qKhatContext khc(G);

  StandardRepK sr=interactive::get_standardrep(khc);

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
  RealReductiveGroup& G = commands::currentRealGroup();
  const KGB& kgb = G.kgb();

  unsigned long x=interactive::get_bounded_int
    (interactive::sr_input(),"Choose KGB element: ",G.kgb().size());

  prettyprint::printVector(std::cout<<"2rho = ",G.rootDatum().twoRho())
    << std::endl;

  Weight lambda=
    interactive::get_weight(interactive::sr_input(),
			    "Give lambda-rho: ",
			    G.rank());
  standardrepk::KhatContext khc(G);

  StandardRepK sr=khc.std_rep_rho_plus(lambda,kgb.titsElt(x));

  (lambda *= 2) += G.rootDatum().twoRho();
  prettyprint::printVector(std::cout << "Weight (1/2)",lambda);
  prettyprint::printVector(std::cout << " converted to (1/2)",khc.lift(sr));

  const TwistedInvolution& canonical =
    G.complexGroup().twistedInvolution(sr.Cartan());
  if (kgb.involution(x)!=canonical)
    prettyprint::printWeylElt(std::cout << " at involution ",
			      canonical, G.weylGroup());
  std::cout << "\nHeight is " << khc.height(sr) << std::endl;

  khc.go(sr);
}

bool examine(RealReductiveGroup& G)
{
  const WeylGroup& W = G.weylGroup();
  const KGB& kgb=G.kgb();
  size_t l = W.length(kgb.involution(0)),t;
  for (size_t i=1; i<kgb.size(); ++i)
    if ((t=W.length(kgb.involution(i)))<l)
      return false;
    else
      l=t;
  return true;
}

void testrun_f()
{
  unsigned long rank=interactive::get_bounded_int
    (interactive::common_input(),"rank: ",constants::RANK_MAX+1);
  std::cout << "Testing W-length monotonicity.\n";
  for (testrun::LieTypeIterator it(testrun::Semisimple,rank); it(); ++it)
  {
    std::cout<< *it << std::endl;
    size_t count=0;
    for (testrun::CoveringIterator cit(*it); cit(); ++cit)
    {
      if (count>0) std::cout << ',';
      std::cout << ++count;
      PreRootDatum prd = *cit;
      WeightInvolution id(prd.rank()); // identity
      ComplexReductiveGroup G(prd,id);
      for (RealFormNbr rf=0; rf<G.numRealForms(); ++rf)
      {
	RealReductiveGroup G_R(G,rf);
	if (not examine(G_R))
	{
	  lietype::InnerClassType ict;
	  for (size_t i=0; i<it->size(); ++i)
	    ict.push_back('e');
	  lietype::Layout lay(*it,ict);
	  realform_io::Interface itf(G,lay);
	  std::cout << " Failure at real form " << itf.out(rf) << std::endl;
	}
	std::cout << std::flush;
      }
    }
    std::cout << '.' << std::endl;
  }

}

void exam_f()
{
  std::cout << "W-length monotine in KGB? "
            << (examine(commands::currentRealGroup())
		? "yes" : "no")
	    << std::endl;
}

void X_f()
{
  ComplexReductiveGroup& G=commands::currentComplexGroup();
  kgb::global_KGB kgb(G); // build global Tits group, "all" square classes
  ioutils::OutputFile f;
  kgb_io::print_X(f,kgb);
}



TorusElement torus_part
  (const RootDatum& rd,
   const WeightInvolution& theta,
   const RatWeight& lambda, // discrete parameter
   const RatWeight& gamma // infinitesimal char
  )
{
  InvolutionData id(rd,theta);
  Weight cumul(rd.rank(),0);
  arithmetic::Numer_t n=gamma.denominator();
  const Ratvec_Numer_t& v=gamma.numerator();
  const RootNbrSet pos_real = id.real_roots() & rd.posRootSet();
  for (RootNbrSet::iterator it=pos_real.begin(); it(); ++it)
    if (rd.coroot(*it).dot(v) %n !=0) // nonintegral
      cumul+=rd.root(*it);
  // now |cumul| is $2\rho_\Re(G)-2\rho_\Re(G(\gamma))$

  return y_values::exp_pi(gamma-lambda+RatWeight(cumul,2));
}

void embedding_f()
{
  RealReductiveGroup& GR = commands::currentRealGroup();
  ComplexReductiveGroup& G = GR.complexGroup();
  const RootDatum& rd = G.rootDatum();

  Weight lambda_rho =
    interactive::get_weight(interactive::sr_input(),"Give lambda-rho: ",
			    G.rank());

  RatWeight lambda(lambda_rho *2 + rd.twoRho(),2);

  RatWeight nu=
    interactive::get_ratweight
    (interactive::sr_input(),"rational parameter nu: ",rd.rank());

  const KGB& kgb = GR.kgb();
  unsigned long x=interactive::get_bounded_int
    (interactive::common_input(),"KGB element: ",kgb.size());

  WeightInvolution theta = G.involutionMatrix(kgb.involution(x));

  nu = RatWeight // make |nu| fixed by $-\theta$
    (nu.numerator()- theta*nu.numerator(),2*nu.denominator());

  RatWeight gamma = nu + RatWeight
    (lambda.numerator()+theta*lambda.numerator(),2*lambda.denominator());

  ioutils::OutputFile f;
  f << "Infinitesimal character is " << gamma << std::endl;

  SubSystemWithGroup sub = SubSystemWithGroup::integral(rd,gamma);

  WeylWord ww;
  const tits::SubTitsGroup sub_gTg
    (G,sub,G.involutionMatrix(kgb.involution(x)),ww);

  Permutation pi;

  f << "Subsystem on dual side is ";
  if (sub.rank()==0)
    f << "empty.\n";
  else
  {
    f << "of type " << dynkin::Lie_type(sub.cartanMatrix(),true,false,pi)
      << ", with roots ";
    for (weyl::Generator s=0; s<sub.rank(); ++s)
      f << sub.parent_nr_simple(pi[s]) << (s<sub.rank()-1 ? "," : ".\n");
  }
  f << "Twisted involution in subsystem: " << ww << ".\n";

  weyl::TI_Entry::Pooltype pool;
  HashTable<weyl::TI_Entry,unsigned int> hash_table(pool);
  std::vector<unsigned int> stops;
  {
    std::vector<TwistedInvolution> queue(1,TwistedInvolution());
    size_t qi = 0; // |queue| inspected up to |qi|
    const TwistedWeylGroup& tW = sub_gTg; // base object, subgroup

    while (qi<queue.size())
    {
      weyl::TI_Entry tw=queue[qi++];
      if (hash_table.find(tw)==hash_table.empty) // Cartan class not seen yet
      { stops.push_back(hash_table.size());
	// generate like kgb::FiberData::complete_class|
	for (unsigned int i=hash_table.match(tw); i<hash_table.size(); ++i)
	{
	  tw=pool[i];
	  for (weyl::Generator s=0; s<tW.rank(); ++s)
	    if (tW.hasTwistedCommutation(s,tw)) // |s| real or imaginary
	    {
	      if (not tW.hasDescent(s,tw)) // |s| imaginary
		queue.push_back(tW.prod(s,tw));
	    }
	    else // |s| complex
	      hash_table.match(tW.twistedConjugated(tw,s));
	} // |for(i)|
      } // |if| new Cartan class
    } // |while| queue not completely inspected
  } // forget |queue|
  stops.push_back(~0); // sentinel
  size_t si=0; // index into |stops|
  for (unsigned int i=0; i<pool.size(); ++i)
  {
    if (i==stops[si])
    { (std::ostream&)f << std::endl; ++si; } // separate classes
    TwistedInvolution tw=pool[i];
    f << sub_gTg.base_point_offset(tw).log_2pi()
      << "  \t@  " << sub_gTg.word(tw) <<std::endl;
  }
} // |embedding_f|

// help functions for |test_f|

// find lift |h| of |t| such that |theta_t*h = -h|; |t| is |theta_t|-stable
Coweight minus_stable_lift(TorusPart t, const CoweightInvolution& theta_t)
{
  Coweight h(t.size(),0);
  for (TorusPart::base_set::iterator it = t.data().begin(); it(); ++it)
    h[*it]=1; // lift torus part to a coweight
  CoweightInvolution theta1 = theta_t;
  for (unsigned int i=0; i<theta1.numRows(); ++i)
    theta1(i,i)+=1; // add identity
  Coweight deviation = theta1*h; // failure of |h| to be |-theta| fixed
  deviation /= 2; // may, but should not, |throw std::runtime_error|
  Coweight correction = matreduc::find_solution(theta1,deviation); // idem
  correction *= 2; // the correction must be by an element of $2X_*$
  h -= correction;
  assert ( theta_t*h == -h );
  assert(TorusPart(h) == t); // the correction kept |h| a lift of |t|
  return h;
}


struct z_data { TorusPart t; Coweight h; Weight lamnum; int z; };

void z_choice(TorusPart t, const ComplexReductiveGroup& G,
	      const CoweightInvolution& theta_t, const Weight& lambda,
	      z_data& out)
{
  const TitsGroup& Tg= G.titsGroup();
  out.t = t; out.lamnum=lambda;
  t -= Tg.twisted(t);
  Coweight h = minus_stable_lift(t,theta_t);
  out.h=h;
  int z =
    arithmetic::remainder(lambda.dot(h+G.distinguished().transposed()*h),4);
  out.z=z; // reinterpreted modulo 8, thus choosing the "positive" square root
}

void test_f()
{
  ioutils::OutputFile file;
  const ComplexReductiveGroup& G = commands::currentComplexGroup();
  const RootDatum& rd = G.rootDatum();
  const TitsCoset& bTg = commands::currentRealGroup().basedTitsGroup();
  const KGB& kgb = commands::currentRealGroup().kgb();
  param_block& block = commands::current_param_block();
  const RatWeight& gamma = block.gamma();
  const Weight gam_num (gamma.numerator().begin(),gamma.numerator().end());
  ext_block::extended_block eblock(block,G.twistedWeylGroup());
  std::vector<z_data> data(eblock.size());
  for (BlockElt n=0; n<eblock.size(); ++n)
  {
    BlockElt fix = eblock.z(n);
    const CoweightInvolution theta_t =
      G.involution_table().matrix(block.involution(fix)).transposed();
    const RatWeight& lambda = block.lambda(fix);
    const Weight lamnum (lambda.numerator().begin(),lambda.numerator().end());

    // just transmitting the torus part $t$ suffices, because one has
    // $\sigma_w\delta((\sigma_w)^{-1})=1$ as fix implies $\delta(w)=w$
    z_choice(kgb.torus_part(block.parent_x(fix)),G,theta_t,lamnum,data[n]);
    file << fix << ": " << data[n].z << '/' << 4*lambda.denominator()
	 << std::endl;
  }
  // now check how cross links relate to the choices made

  CoweightInvolution deltr1 = G.distinguished().transposed();
  for (unsigned int i=0; i<deltr1.numRows(); ++i)
    deltr1(i,i)+=1; // add identity
  for (BlockElt n=0; n<eblock.size(); ++n)
  {
    const z_data& here=data[n];
    BlockElt fix = eblock.z(n);
    const CoweightInvolution theta_t =
      G.involution_table().matrix(block.involution(fix)).transposed();
    std::cout << eblock.z(n) << ": ";
    for (weyl::Generator s=0; s<eblock.rank(); ++s)
      if (eblock.cross(s,n)!=UndefBlock)
      {
	BlockElt sn = eblock.cross(s,n);
	const z_data& there=data[sn];
	TitsElt a = kgb.titsElt(block.parent_x(eblock.z(sn)));
	ext_gen g = block.orbit(s);
	for (unsigned int i=0; i<g.w_tau.size(); ++i)
	  bTg.strict_based_twisted_conjugate(a,g.w_tau[i]); // pull |a| here

	TorusPart shift = here.t+bTg.titsGroup().left_torus_part(a);
	Coweight lift = minus_stable_lift(shift,theta_t);
	Coweight wh = there.h;
	for (unsigned int i=0; i<g.w_tau.size(); ++i)
	  rd.simpleCoreflect(wh,g.w_tau[i]); // SHOULD use subsystem....
	wh += lift; // multiply (under $exp(\pi\ii/2 * )$) by |lift|
	wh -= G.distinguished().transposed()*lift; // and by $lift^{-\delta}$
	// now |wh| is compatible with the torus part choice down |here|

	Coweight offset = here.h-wh;
	(offset - theta_t*offset)/=4; // |assert((1-theta_t)*offset % 4 == 0)|;
	int Lambda_h_over_h = here.lamnum.dot(offset); // significant modulo 4

	// now compute $(w.\gamma-\gamma)(h)$, by which $z$ gets multiplied
	Weight gn =gam_num;
	for (unsigned int i=0; i<g.w_tau.size(); ++i)
	  rd.simpleReflect(gn,g.w_tau[i]);
	gn -= gam_num;
	gn /= gamma.denominator(); // the result should lie in root lattice
	int zz = gn.dot(there.h); // offset, significant modulo 4
	zz = there.z + 2*zz ; // transport here, result is mod 8
	assert ( (here.lamnum.dot(deltr1*wh)-zz)%4==0 );
	int compare = (zz - here.z) - 2*Lambda_h_over_h;
	compare = arithmetic::remainder(compare,8);
	assert ( compare%4 == 0);
	if (compare!=0)
	  std::cout << (int)s << '(' << eblock.z(sn) << ") ";
      }
    std::cout << std::endl;
  }
} // |test_f|


} // |namespace|


} // |namespace test|

} // |namespace atlas|
