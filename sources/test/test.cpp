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
#include "innerclass.h" // |twoRho| in |nu_block::print|
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

#include "../Atlas.h" // here to preempt double inclusion of _fwd files

#include "free_abelian.h"
#include "permutations.h"
#include "matreduc.h"

#include "dynkin.h"
#include "prerootdata.h"
#include "rootdata.h"
#include "cartanclass.h"
#include "innerclass.h"
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
#include "output.h"
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
  //  void block_braid_f();
  void repr_braid_f();

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
  void testrun_f();
  void exam_f();

  void X_f();

/*
  For convenience, the "test" command is added to the mode that is flagged by
  the testMode constant defined here; therefore "test" appears (conditionally)
  in every template instance of |addTestCommands| and of |addTestHelp| below.
  Set this constant according to the requirements of the |test_f| function.
*/
  enum TestMode {EmptyMode, MainMode, RealMode, BlockMode, ReprMode,
		 numTestMode};
  const TestMode testMode = EmptyMode; // currently does nothing, so empty mode

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
  mode.add("testrun",testrun_f,
	   "iterates over root data of given rank, calling examine",
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

  mode.add("examine",exam_f,
	   "tests whether x0 will change",commands::use_tag);


  if (testMode == RealMode)
    mode.add("test",test_f,test_tag);

}


// Add to the block mode the test commands that require that mode.
template<>
void addTestCommands<commands::BlockmodeTag> (commands::CommandNode& mode)
{
#if 0 // function cannot be defined without sign-fixing for such blocks
  mode.add("bbraid",block_braid_f,
	   "tests braid relations on an extended block",commands::use_tag);
#endif

  if (testMode == BlockMode)
    mode.add("test",test_f,test_tag);

}

// Add to the repr mode the test commands that require that mode.
template<>
void addTestCommands<commands::ReprmodeTag> (commands::CommandNode& mode)
{
  mode.add("braid",repr_braid_f,
	   "tests braid relations on an extended block",commands::use_tag);
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

bool examine(RealReductiveGroup& G);
void testrun_f()
{
  unsigned long rank=interactive::get_bounded_int
    (interactive::common_input(),"rank: ",constants::RANK_MAX+1);
  std::cout << "Testing x0 torus bits.\n"; // adapt to |examine|
  for (testrun::LieTypeIterator it(testrun::Semisimple,rank); it(); ++it)
  {
    std::cout<< *it << std::endl;
    for (testrun::CoveringIterator cit(*it); cit(); ++cit)
    {
      PreRootDatum prd = *cit;
      WeightInvolution id(prd.rank()); // identity
      InnerClass G(prd,id);
      for (RealFormNbr rf=0; rf<G.numRealForms(); ++rf)
      {
	RealReductiveGroup G_R(G,rf);
	if (not examine(G_R))
	{
	  WeightList subattice_basis;
	  cit.makeBasis(subattice_basis);
	  basic_io::seqPrint(std::cout,
			     subattice_basis.begin(),subattice_basis.end(),
			     ", ","\n Sublattice basis: ","\n");
	  lietype::InnerClassType ict; // need layout to convert form number
	  for (size_t i=0; i<it->size(); ++i)
	    ict.push_back('e');
	  lietype::Layout lay(*it,ict);
	  output::FormNumberMap itf(G,lay);
	  std::cout << " Failure at real form " << itf.out(rf) << std::endl;
	}
	std::cout << std::flush;
      }
    }
    std::cout << '.' << std::endl;
  }

}


// Main mode functions

// Print the roots in the simple root coordinates.
void roots_rootbasis_f()
{
  const RootSystem& rs =  commands::current_inner_class().rootSystem();
  ioutils::OutputFile file;

  for (RootNbr i=0; i<rs.numRoots(); ++i)
    prettyprint::printInRootBasis(file,i,rs) << std::endl;
}

// Print the positive roots in the simple root coordinates.
void posroots_rootbasis_f()

{
  const RootSystem& rs = commands::current_inner_class().rootSystem();

  ioutils::OutputFile file;
  prettyprint::printInRootBasis(file,rs.posRootSet(),rs);
}

// Print the coroots in the simple coroot coordinates.
void coroots_rootbasis_f()
{
  const RootSystem rs (commands::current_inner_class().dualRootSystem());

  ioutils::OutputFile file;
  for (RootNbr i=0; i<rs.numRoots(); ++i)
    prettyprint::printInRootBasis(file,i,rs) << std::endl;

}

// Print the positive coroots in the simple coroot coordinates.
void poscoroots_rootbasis_f()
{
  const RootSystem rs (commands::current_inner_class().dualRootSystem());

  ioutils::OutputFile file;
  prettyprint::printInRootBasis(file,rs.posRootSet(),rs);
}


void X_f()
{
  InnerClass& G=commands::current_inner_class();
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
    if (not khc.isFinal(sr))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.witness()) << ".\n";
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
    if (not khc.isFinal(sr))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.witness()) << ".\n";
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
    if (not khc.isStandard(sr))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not standard, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.witness())
	<< ".\n";
      return;
    }
    if (not khc.isFinal(sr))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.witness()) << ".\n";
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
    khc.K_type_matrix(singleton,bound,new_order,nullptr);
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
    if (not khc.isNormal(sr))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not normal, as witnessed by coroot sum "
	<< khc.info(sr.Cartan()).coroot_sum(khc.witness())
	<< ".\n";
      return;
    }
    if (not khc.isStandard(sr))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not standard, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.witness())
	<< ".\n";
      return;
    }
    if (not khc.isFinal(sr))
    {
      khc.print(std::cout << "Representation ",sr)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.witness()) << ".\n";
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
    khc.K_type_matrix(singleton,bound,new_order,nullptr);
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

  CoeffList factor;
  int_Matrix b = matreduc::adapted_basis(G.cartan(cn).involution()-1,factor);

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

  StandardRepK srk=interactive::get_standardrep(khc);

  {
    if (not khc.isStandard(srk))
    {
      khc.print(std::cout << "Representation ",srk)
        << " is not standard, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.witness())
	<< ".\n";
      return;
    }
    if (not khc.isFinal(srk))
    {
      khc.print(std::cout << "Representation ",srk)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.witness()) << ".\n";
      return;
    }
  }

  khc.normalize(srk); // needed so |K_type_formula| gives itself as leading term
  standardrepk::combination c=khc.standardize(srk);

  if (c.empty())
  {
    khc.print(std::cout << "Representation ",srk) << " is zero.\n";
    return;
  }

  assert(c.size()==1 and khc.rep_no(c.begin()->first)==srk);

  khc.print(std::cout << "Height of representation ",srk) << " is "
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

  StandardRepK srk=interactive::get_standardrep(khc);

  {
    if (not khc.isStandard(srk))
    {
      khc.print(std::cout << "Representation ",srk)
        << " is not standard, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.witness())
	<< ".\n";
      return;
    }
    if (not khc.isFinal(srk))
    {
      khc.print(std::cout << "Representation ",srk)
        << " is not final, as witnessed by coroot "
	<< G.rootDatum().coroot(khc.witness()) << ".\n";
      return;
    }
  }

  khc.normalize(srk); // needed so |K_type_formula| gives itself as leading term
  standardrepk::q_combin c=khc.standardize(srk);

  if (c.empty())
  {
    khc.print(std::cout << "Representation ",srk) << " is zero.\n";
    return;
  }

  assert(c.size()==1 and khc.rep_no(c.begin()->first)==srk);

  khc.print(std::cout << "Height of representation ",srk) << " is "
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
    G.innerClass().involution_of_Cartan(sr.Cartan());
  if (kgb.involution(x)!=canonical)
    prettyprint::printWeylElt(std::cout << " at involution ",
			      canonical, G.weylGroup());
  std::cout << "\nHeight is " << khc.height(sr) << std::endl;

  khc.go(sr);
}

bool examine(RealReductiveGroup& G)
{
  const KGB& kgb=G.kgb();
  TorusPart t0 = kgb.torus_part(0);
  TorusPart t1 = G.innerClass().x0_torus_part(G.realForm());
  return t0==t1;
}

void exam_f()
{
  RealReductiveGroup& G = commands::currentRealGroup();
  if (examine(G))
    std::cout << "x0 torus bits constistent with traditional ones";
  else
    std::cout << "x0 torus bits changed from " << G.kgb().torus_part(0)
	      << " to " << G.innerClass().x0_torus_part(G.realForm());
  std::cout << std::endl;
}

void test_f() // trial of twisted KLV computation
{
  std::cout << "No test function in this version of  Fokko" << std::endl;
  return;
  // this function is disabled as it fails for the following reason:
  // currently we cannot correct signs in block mode extended block
  ext_block::ext_block
    eblock(commands::current_inner_class(),
	   commands::currentBlock(),
	   commands::currentRealGroup().kgb(),
	   commands::currentDualRealGroup().kgb(),
	   commands::current_inner_class().distinguished()
	   );

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


// Check braid relations extended block
int test_braid(const ext_block::ext_block& eblock)
{
  int failed=0;
  std::cout << "testing quadratic relations" << std::endl;
  for (BlockElt x=0; x<eblock.size(); ++x)
    for (weyl::Generator s=0; s<eblock.rank(); ++s)
      if (not ext_block::check_quadratic(eblock,s,x))
      {
	std::cout <<  "Quadratic relation failure: " << eblock.z(x)
		  << ", s=" << (unsigned)s << std::endl;
	++failed;
      }
  if (failed>0)
    return failed; // avoid testing braid relations when quadratics fail

  std::cout << "testing braids" << std::endl;
  int count=0;
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
	      ++failed;
	      std::cout <<  "Braid relation failure: " << eblock.z(x)
			<< ", s=" << (unsigned)s << ", t=" << (unsigned)t;
	      for (BitMap::iterator it=cluster.begin(); it(); ++it)
		std::cout << (it==cluster.begin() ? " (" : ",")
			  << eblock.z(*it) ;

	      std::cout << ')' << std::endl;
	      seen |= cluster; // don't do elements of failed cluster again
	    }
	  }
    }
  if (failed==0 and count>0)
  { if (count==1)
      std::cout << "The braid relation holds!\n";
    else
      std::cout << "All " << count << " braid relations hold!\n";
    std::cout << std::endl;
  }
  return failed;
} // |test_braid|

#if 0 // call of |check| below needs implementing
void block_braid_f() // in block mode
{
  WeightInvolution delta = interactive::get_commuting_involution
    (commands::current_layout(), commands::current_lattice_basis());

  auto& block = commands::currentBlock();

  ext_block::ext_block eblock(commands::current_inner_class(),block,
			      commands::currentRealGroup().kgb(),
			      commands::currentDualRealGroup().kgb(),
			      delta);
  if (check(eblock,block,true))
    test_braid(eblock);
}
#endif

void repr_braid_f()
{
  commands::ensure_full_block();
  WeightInvolution delta = interactive::get_commuting_involution
    (commands::current_layout(), commands::current_lattice_basis());

  auto& block = commands::current_param_block();
  if (not ((delta-1)*block.gamma().numerator()).isZero())
  {
    std::cout << "Chosen delta does not fix gamma=" << block.gamma()
	      << " for the current block." << std::endl;
    return;
  }
  ext_block::ext_block eblock(commands::current_inner_class(),block,delta,true);
  test_braid(eblock);
}



// Block mode functions

} // |namespace|


} // |namespace test|

} // |namespace atlas|
