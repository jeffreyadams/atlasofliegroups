/*
  This is test.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/


#include "block_io.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <cassert>
#include <map>

#include "polynomials.h"
#include "kgb.h"     // |kgb.size()|
#include "complexredgp.h" // |twoRho| in |nu_block::print|
#include "blocks.h"
#include "ext_block.h"
#include "kl.h"
#include "repr.h"

#include "basic_io.h"	// operator |<<|
#include "ioutils.h"	// |digits|
#include "prettyprint.h" // |printWeylElt|
//extra

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
#include "tits.h"
#include "weyl.h"
#include "involutions.h" // for |InvolutionTable|
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

#include "ext_kl.h"

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
  void braid_f();
  void go_f();
  void mytest_f();
  void fix_braid_f();

  void roots_rootbasis_f();
  void coroots_rootbasis_f();
  void posroots_rootbasis_f();
  void poscoroots_rootbasis_f();
  void Ktypeform_f();
  void qKtypeform_f();
  void Ktypemat_f();
  void qKtypemat_f();
  void mod_lattice_f();
  void branch_f();
  void qbranch_f();
  void srtest_f();

  void X_f();

/*
  For convenience, the "test" command is added to the mode that is flagged by
  the testMode constant defined here; therefore "test" appears (conditionally)
  in every template instance of |addTestCommands| and of |addTestHelp| below.
  Set this constant according to the requirements of the |test_f| function.
*/
  enum TestMode {EmptyMode, MainMode, RealMode, BlockMode, ReprMode,
		 numTestMode};
  const TestMode testMode = BlockMode; // currently does test of extended block

  // utilities
  const RootDatum& currentRootDatum();

} // |namespace|

/*****************************************************************************

        Chapter I -- Functions declared by test.h

  This section defines the functions declared in test.h :

    - addTestCommands() : adds the test commands to the main command tree;

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
  mode.add("go",go_f,
	   "generates difficult SO(5,5) extended block, and runs 'braid'",
	   commands::use_tag);
  if (testMode == EmptyMode)
    mode.add("test",test_f,test_tag);
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

  mode.add("X",X_f,"prints union of K\\G/B for real forms in inner class",
	   commands::std_help);

  if (testMode == MainMode)
    mode.add("test",test_f,test_tag);

}

// Add to the real mode node the test commands that require that mode.
template<>
void addTestCommands<commands::RealmodeTag> (commands::CommandNode& mode)
{
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

  if (testMode == RealMode)
    mode.add("test",test_f,test_tag);

}


// Add to the block mode the test commands that require that mode.
template<>
void addTestCommands<commands::BlockmodeTag> (commands::CommandNode& mode)
{
  mode.add("braid",braid_f,
	   "tests braid relations on an extended block",commands::use_tag);
  if (testMode == BlockMode)
    mode.add("test",test_f,test_tag);

  mode.add("mytest",mytest_f,
	   "my test command",commands::use_tag);
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

  std::ostream& print(std::ostream& strm);



// Empty mode functions


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


void X_f()
{
  ComplexReductiveGroup& G=commands::currentComplexGroup();
  kgb::global_KGB kgb(G); // build global Tits group, "all" square classes
  ioutils::OutputFile f;
  kgb_io::print_X(f,kgb);
}


// Real mode functions

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

void test_f() // trial of twisted KLV computation
{

  ext_block::extended_block
    eblock(commands::currentBlock(),
	   commands::currentComplexGroup().twistedWeylGroup());

  BlockElt last; input::InputBuffer& cl= commands::currentLine();
  cl >> last; // maybe get threshold for filling
  if (cl.fail())
    last=eblock.size();
  else // convert to block number
    last=eblock.element(last);

  std::vector<ext_kl::Pol> pool;
  ext_kl::KL_table twisted_KLV(eblock,pool);
  twisted_KLV.fill_columns(last);

  ioutils::OutputFile f;
  for (BlockElt y=0; y<last; ++y)
    for (BlockElt x=y+1; x-->0; )
      if (not twisted_KLV.P(x,y).isZero())
      {
	f << "P(" << eblock.z(x) << ',' << eblock.z(y) << ")=";
	f << twisted_KLV.P(x,y) << std::endl;
      }

}

bool isDirectRecursion(ext_block::DescValue v)
{
  return is_descent(v) and is_unique_image(v);
}


void mytest(KGB kgb){
  //  prettyprint::printInvolution(strm,kgb.involution(1),kgb.twistedWeylGroup())<< std::endl;

  }
void mytest_f(){
  // RealReductiveGroup& G = commands::currentRealGroup();
  // ComplexReductiveGroup& GC=commands::currentComplexGroup();
  // const global_KGB& global_kgb = global_KGB(GC);
  // const KGB& kgb = G.kgb();
  // const TitsGroup& Tg = G.titsGroup();
  // const GlobalTitsGroup& global_Tg = GlobalTitsGroup(GC);
  // ioutils::OutputFile f;

  // RatWeight zero=global_kgb.torus_part(0).log_2pi();
  // WeightInvolution delta = G.distinguished();
  // std::cout<<  "delta:" << std::endl;
  // for (unsigned i=0; i<delta.numRows(); ++i)
  // 	{   std::cout << " " << std::endl;
  // 	  for (unsigned j=0; j<delta.numColumns(); ++j)
  // 	    std::cout << delta(i,j);}


  // for (size_t i=0; i<kgb.size(); ++i){
  //   KGBElt dual_x=kgb.Hermitian_dual(i);
  //   if (dual_x==i){
  //     f << "" << std::endl << "x=" << i << " w=";
  //     prettyprint::printInvolution(f,kgb.involution(i),kgb.twistedWeylGroup());      

  //     //    WeightInvolution matrix=kgb.involution_matrix(i);
  //     //    RatWeight shift=y_shift-matrix*y_shift;
      
  //   //    TorusPart t=kgb.torus_part(i);
  //   //TorusElement global_t=global_kgb.torus_part(i);
  //     //    RootDatum rd=GC.rootDatum();

  //     TorusElement s=global_kgb.torus_part(i);
  //     WeightInvolution theta=kgb.involution_matrix(i);
  //     RatWeight log_s=s.log_2pi();

  //     RatWeight gamma=(theta*log_s-log_s);
  //     RatWeight log_h=(gamma-delta*gamma);
  //     log_h=RatWeight(log_h.numerator(),4*log_h.denominator()).normalize();

  //     RatWeight log_s_new=(log_s+theta*log_s);
  //     log_s_new=RatWeight(log_s_new.numerator(),2*log_s_new.denominator()).normalize();
  //     bool equal=(log_s==log_s_new);
  //     //      std::cout << "  s=" << log_t << "  t_new: " << t_new;
  //     std::cout << "  s=" << log_s << "  ";
  //     if (!equal) std::cout << " s'=" << log_s_new;
  //     if (log_h != zero)      std::cout << " " <<  "  h=" << log_h;
  //     RatWeight log_t=log_h+delta*log_h;
  //     log_t=RatWeight(log_t.numerator(),log_t.denominator()).normalize();
  //     if (log_t != zero) std::cout << " t=" << log_t;
      


      

  //   // WeightInvolution matrix=kgb.involution_matrix(i);
  //   // if (matrix*rw==rw) std::cout << "yes equal"; else std::cout << "NOT not equal";
  //   // std::cout << "rw=" << rw << "theta(rw)=" << matrix*rw << std::endl;



  //   //    TorusElement t=kgb.torus_part_global(rd,i);
  //   //    TwistedInvolution ti=kgb.involution(i);
  //   //    WeylElt w=WeylElt(ti);
  //   //    GlobalTitsElement xi=GlobalTitsElement(t,w);

  //   //    TorusPart tp=kgb.torus_part(i);
  //   //    TitsElt xi_tits=TitsElt(Tg,tp,w);
  //   //    WeylElt wxi=xi_tits.w();



  //   }

  //   //  f << "" << std::endl;
  // }
}


void mytest_f_old(){


  RealReductiveGroup& G = commands::currentRealGroup();
  ComplexReductiveGroup& GC=commands::currentComplexGroup();
  const global_KGB& global_kgb = global_KGB(GC);
  const KGB& kgb = G.kgb();
  
  const TitsGroup& Tg = G.titsGroup();
  const GlobalTitsGroup& global_Tg = GlobalTitsGroup(GC);
  ioutils::OutputFile f;

  TorusElement yrho =y_values::exp_2pi(global_kgb.globalTitsGroup().torus_part_offset());
  //  const GlobalTitsElement& hdelta=GlobalTitsElement(yrho);
  //  GlobalTitsElement minus_hdelta=GlobalTitsElement(minus_yrho);
  RatWeight y_shift=yrho.log_2pi();
  f << "m_rho=exp(i\\pi(" << y_shift << "), terms with h=m_rho labelled *" << std::endl;
  static TorusElement id =global_kgb.torus_part(0);

  for (size_t i=0; i<kgb.size(); ++i){
    KGBElt dual_x=kgb.Hermitian_dual(i);
        if (dual_x==i){
	  f << "" << std::endl << "x=" << i << " w=";
      prettyprint::printInvolution(f,kgb.involution(i),kgb.twistedWeylGroup());      

      //    WeightInvolution matrix=kgb.involution_matrix(i);
      //    RatWeight shift=y_shift-matrix*y_shift;
      
    //    TorusPart t=kgb.torus_part(i);
    //TorusElement global_t=global_kgb.torus_part(i);
    RootDatum rd=GC.rootDatum();




    TorusElement t=kgb.torus_part_global(rd,i);
    TwistedInvolution ti=kgb.involution(i);
    WeylElt w=WeylElt(ti);
    GlobalTitsElement xi=GlobalTitsElement(t,w);

    //    TorusPart tp=kgb.torus_part(i);
    //    TitsElt xi_tits=TitsElt(Tg,tp,w);
    //    WeylElt wxi=xi_tits.w();

    

    RatWeight log_torus_part_xi=t.log_2pi()+y_shift;
    f << " h(xi): ";
    if (t==id){
      f << "*";
    }



  //    eblock.toggle_edge(153,213); // 4, 2Ci/2Cr
  //    eblock.toggle_edge(254,196); // 4, 2Ci/2Cr
  //    eblock.toggle_edge(240,284); // 4, 2Ci/2Cr




  //    eblock.toggle_edge(153,213); // 4, 2Ci/2Cr
  //    eblock.toggle_edge(254,196); // 4, 2Ci/2Cr
  //    eblock.toggle_edge(240,284); // 4, 2Ci/2Cr




  //    eblock.toggle_edge(153,213); // 4, 2Ci/2Cr
  //    eblock.toggle_edge(254,196); // 4, 2Ci/2Cr
  //    eblock.toggle_edge(240,284); // 4, 2Ci/2Cr




  //    eblock.toggle_edge(153,213); // 4, 2Ci/2Cr
  //    eblock.toggle_edge(254,196); // 4, 2Ci/2Cr
  //    eblock.toggle_edge(240,284); // 4, 2Ci/2Cr

    f << log_torus_part_xi << " ";

     GlobalTitsElement delta_xi=global_Tg.twisted(xi);
     TorusElement torus_part_delta_xi=delta_xi.torus_part();
     RatWeight log_torus_part_delta_xi=torus_part_delta_xi.log_2pi()+y_shift;
     f << "h(delta(xi)):" << log_torus_part_delta_xi << " ";
     RatWeight difference=log_torus_part_delta_xi-log_torus_part_xi;
     
     if (xi==delta_xi)
       f << "true";
     else
       f <<  "false  " << difference;

     const WeylGroup& W = Tg.weylGroup();
     // test Cayleys                         
     //     weyl::Generator s;
     unsigned s;
     for (s=0; s<G.semisimpleRank(); ++s)
       {if (kgb.status(s,i)==gradings::Status::Real)
	   {  f << "" << std::endl;
	     f << " s: " << s << " " ;
	     WeylElt sw=w;
	     W.leftMult(sw,s);
	     GlobalTitsElement my_cayley=GlobalTitsElement(t,sw);
	     GlobalTitsElement cayley=xi;
	     global_Tg.do_inverse_Cayley(s,cayley);

	     TwistedInvolution xi_w=xi.tw();
	     TwistedInvolution cayley_wtemp=cayley.tw();
	     f << "first ";
	     prettyprint::printInvolution(f,xi_w,kgb.twistedWeylGroup());      
	     f << "second ";
	     prettyprint::printInvolution(f,cayley_wtemp,kgb.twistedWeylGroup());      

	     TorusElement torus_part_my_cayley=my_cayley.torus_part();
	     RatWeight log_torus_part_my_cayley=torus_part_my_cayley.log_2pi()+y_shift;
	     TwistedInvolution my_cayley_w=my_cayley.tw();
	     f << "  my_cayley:" << log_torus_part_my_cayley << " ";
	     prettyprint::printInvolution(f,my_cayley_w,kgb.twistedWeylGroup());      
	     f << "  ";

	     TorusElement torus_part_cayley=cayley.torus_part();
	     RatWeight log_torus_part_cayley=torus_part_cayley.log_2pi()+y_shift;
	     TwistedInvolution cayley_w=cayley.tw();
	     f << "cayley: " << log_torus_part_cayley << " ";
	     prettyprint::printInvolution(f,cayley_w,kgb.twistedWeylGroup());      
	     f << "  ";
	     
	     if (my_cayley==cayley)
	       f << "cayleys are equal ";
	     else
	       f << "cayleys not equal ";
	     

	     //	     WeylElt w = WeylElt();  // the identity
	     //W.leftMult(w,s);   // s as an element of W
	     


	     //	     TitsElt muin=TitsElt(Tg,w);
	     //	     f << "  muin: ";
	     //	     prettyprint::printTitsElt(f,muin,Tg);

	     //	     GlobalTitsElement sigma=GlobalTitsElement(w,rd.rank());
	     //	     TitsElt sigma_tits=TitsElt(Tg,w);

	     //	     GlobalTitsElement my_Cayley=global_Tg.prod(sigma,xi);
			      //	     GlobalTitsElement cayley=xi;
			      //	     global_Tg.do_inverse_Cayley(s,cayley);

			      //	     if (my_Cayley==cayley)
			      //	       f << " cayleys are equal ";
	     

	     //	     TitsElt mu=TitsElt(Tg,w);
	     //	     f << "  mu: ";
	     //	     prettyprint::printTitsElt(f,mu,Tg);

	   }



	 //	 GlobalTitsElement test=global_Tg.prod(axi,xi);
	 
	 //std::ostream& printTitsElt(std::ostream& strm, const TitsElt& a,const TitsGroup& Tg)

      //      W.leftMult(xi.w,s);
      //GlobalTitsGroup::do_inverse_Cayley(weyl::Generator s,GlobalTitsElement& a)
	 //      global_Tg.do_inverse_Cayley(s,xi);
      //      WeylElt test=prod(WeylElt(),s);
      //      GlobalTitsElement sigma=GlobalTitsElement(s,rd.rank());
      //(const WeylElt& we,size_t rank) : t(rank),w(we) {}

       //      global_Tg.mult_sigma(xi,s);
      //      mult_sigma(TitsElt&, weyl::Generator);
      
    //	  assert(kgb.status(sub.simple(s),conj_x)==gradings::Status::Real);
	}




  // RatWeight twisted_weight=twisted_te.log_2pi();
  // TorusElement modified=y_values::exp_2pi(shift+twisted_weight);
  // GlobalTitsElement modified_tits_element=GlobalTitsElement(modified,w);

  // f << "  modified: " << shift+twisted_weight+y_shift << " ";
  
  // if (modified_tits_element==global_te)
  //   f << "  modified is true"; 
  // else
  //   f << "  modified is false";


      }//  if (dual_x==i){




  }
  //    TorusElement torus_twist=y_values::exp_2pi(shift);
    
  




  f << "" << std::endl;
}

  


// Check for nasty endgame cases in block
void test_braid(ext_block::extended_block my_eblock)
{
  ext_block::extended_block eblock=my_eblock;
  std::cout << "testing braids" << std::endl;
  bool OK=true; int count=0;
  for (weyl::Generator t=1; t<eblock.rank(); ++t)
    for (weyl::Generator s=0; s<t; ++s)
    {
      BitMap seen(eblock.size());
      for (BlockElt x=0; x<eblock.size(); ++x)
	if (not seen.isMember(x))
	  {
	    BitMap cluster(eblock.size());
	    if (ext_block::check_braid(eblock,s,t,x,cluster))
	      {
		++count;
		// if ((eblock.z(x)==59 or eblock.z(x)==97 or eblock.z(x)==237 or eblock.z(x)==32) and (s+1==3 or s+1==2) and t+1==4){
		//    std::cout << std::endl << "Braid relation SUCCESS: " << eblock.z(x)
		//  	    << ", s=" << s+1 << ", t=" << t+1;
		//    for (BitMap::iterator it=cluster.begin(); it(); ++it)
		//      std::cout << (it==cluster.begin() ? " (" : ",")
		//  	      << eblock.z(*it) ;
		//    std::cout << ')' << std::endl << std::endl << std::endl;
		//    seen |= cluster; // don't do elements of same cluster again}
		// }
	      }
	    else
	    {
	      OK = false;
	      std::cout << std::endl << "Braid relation failure: " << eblock.z(x)
			<< ", s=" << s+1 << ", t=" << t+1;
	      for (BitMap::iterator it=cluster.begin(); it(); ++it){
		std::cout << (it==cluster.begin() ? " (" : ",")
			  << eblock.z(*it) ;
		//		const ext_block::DescValue type = eblock.descent_type(t,eblock.z(*it));
		//		std::cout << " " << descent_code(type);
	      }
	      std::cout << ')';
	      // for (BitMap::iterator iter=cluster.begin(); iter(); ++iter){
	      // 	BlockElt x=eblock.z(*iter);
	      // 	const ext_block::DescValue type = eblock.descent_type(t,x);
	      // 	if (type==atlas::ext_block::two_semi_imaginary){
	      // 	  BlockElt y=eblock.Cayley(t,*iter);
	      // 	  std::cout << "y=" << eblock.Cayley(t,x) << "  " << eblock.z(y) << std::endl;
	      // 	  eblock.toggle_edge(x,eblock.z(y));

	      // 	  break;
	      // 	}
	      // }
	      //	      std::cout << ')' << std::endl;

	      seen |= cluster; // don't do elements of same cluster again
	    }
	  }
    }
  if (OK)
    std::cout << "All " << count << " relations hold!\n";
  std::cout << std::endl;
} // |braid_f|



ext_block::extended_block fix_braid(ext_block::extended_block eblock)
{
  std::cout << std::endl << "Fixing braids" << std::endl;
  bool OK=true; int count=0;
  for (weyl::Generator t=1; t<eblock.rank(); ++t)
    for (weyl::Generator s=0; s<t; ++s)
    {
      BitMap seen(eblock.size());
      for (BlockElt x=0; x<eblock.size(); ++x)
	if (not seen.isMember(x))
	  {
	    BitMap cluster(eblock.size());
	    if (ext_block::check_braid(eblock,s,t,x,cluster))
		++count;
	    else
	    {
	      OK = false;
	      std::cout << std::endl << "Braid relation failure: " << eblock.z(x)
			<< ", s=" << s+1 << ", t=" << t+1;
	      for (BitMap::iterator it=cluster.begin(); it(); ++it){
		std::cout << (it==cluster.begin() ? " (" : ",")
			  << eblock.z(*it);
		std::cout << "[" << descent_code(eblock.descent_type(t,*it)) << "]";
	      }
	      std::cout << ")";
	      //	      for (BitMap::iterator iter=cluster.begin(); iter(); ++iter){
	      int found_2Ci=0;
	      for (BitMap::iterator iter=cluster.begin(); iter(); ++iter){
		BlockElt x=eblock.z(*iter);
		const ext_block::DescValue type = eblock.descent_type(t,*iter);
		//		std::cout << " " << descent_code(type) << " ";
		if (type==atlas::ext_block::two_semi_imaginary){
		  if (found_2Ci==0){
		    ++found_2Ci;
		  }else{
		    BlockElt y=eblock.Cayley(t,*iter);
		    eblock.toggle_edge(x,eblock.z(y));
		    std::cout << std::endl;
		  break;
		  }
		}
	      }
	      //	      std::cout << ')' << std::endl;
	      //	      std::cout << ')';
	      seen |= cluster; // don't do elements of same cluster again
	    }
	  }
    }
  if (OK)
    std::cout << "All " << count << " relations hold!\n";
  std::cout << std::endl;
  std::cout << "returning eblock" << std::endl;
  return eblock;
} // |braid_f|



void braid_f()
{
  ext_block::extended_block
    eblock(commands::currentBlock(),
	   commands::currentComplexGroup().twistedWeylGroup());
  test_braid(eblock);
}


void toggle(int v[]){
  int w[4];
  int i;
  for (i=0; i<4; ++i){
    w[i]=v[i];
  }
  if (v[4]==1){
    w[0]=v[0];w[1]=v[3];w[2]=v[1];w[3]=v[3];
  } 
  else if (v[4]==1){
    w[0]=v[1];w[1]=v[2];w[2]=v[2];w[3]=v[3];
  }
  else if (v[4]==2){
    w[0]=v[0];w[1]=v[2];w[2]=v[1];w[3]=v[3];
  }
  else if (v[4]==3){
    w[0]=v[0];w[1]=v[2];w[2]=v[1];w[3]=v[2];
  }
  else if (v[4]==4){
    w[0]=v[0];w[1]=v[2];w[2]=v[1];w[3]=v[3];
  }
  else if (v[4]==5){
    w[0]=v[0];w[1]=v[3];w[2]=v[1];w[3]=v[2];
  }
  else if (v[4]==6){
    w[0]=v[0];w[1]=v[3];w[2]=v[1];w[3]=v[3];
  }
  else if (v[4]==7){
    w[0]=v[0];w[1]=v[3];w[2]=v[1];w[3]=v[3];
  }
  v[0]=w[0];v[1]=w[1];v[2]=w[2];v[3]=w[3];
}

void go_f()
{
  drop_to(commands::empty_mode); // make sure all modes will need re-entering
  //commands::currentLine().str("D6 sc u 2 2"); // type-ahead in input buffer
  //  commands::currentLine().str("D5 sc s 2 3"); // type-ahead in input buffer
  // commands::currentLine().reset(); // and reset to start at beginning
  commands::main_mode.activate();
  commands::real_mode.activate();
  commands::block_mode.activate();

  ext_block::extended_block eblock(commands::currentBlock(),
	   commands::currentComplexGroup().twistedWeylGroup());
  //std::cout <<  "default braids"<< std::endl;
  //WORKS
  //  eblock.toggle_edge(153,213); // 4, 2Ci/2Cr
  //  eblock.toggle_edge(256,198); // 4, 2Ci/2Cr
  //  eblock.toggle_edge(240,284); // 4, 2Ci/2Cr
  //WORKS
  test_braid(eblock);
  // eblock.toggle_edge(153,213); // 4, 2Ci/2Cr
  // eblock.toggle_edge(198,256); // 4, 2Ci/2Cr
  // eblock.toggle_edge(240,284); // 4, 2Ci/2Cr
  // test_braid(eblock);exit(0);

  ext_block::extended_block fixed_eblock=fix_braid(eblock);
  test_braid(fixed_eblock);

  bool polynomials=false;
  if (polynomials){
  std::vector<ext_kl::Pol> pool;
  ext_kl::KL_table twisted_KLV(eblock,pool);
  twisted_KLV.fill_columns();

  ioutils::OutputFile f;
  for (BlockElt y=0; y<twisted_KLV.size(); ++y)
    for (BlockElt x=y+1; x-->0; )
      if (not twisted_KLV.P(x,y).isZero())
      {
	f << "P(" << eblock.z(x) << ',' << eblock.z(y) << ")=";
	f << twisted_KLV.P(x,y) << std::endl;
      }
  }
}


// Block mode functions


} // |namespace|


} // |namespace test|

} // |namespace atlas|
