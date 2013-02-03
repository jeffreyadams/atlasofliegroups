/*
  This is reprmode.cpp

  Copyright (C) 2013 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "reprmode.h"
#include "realmode.h"
#include "mainmode.h"
#include "blockmode.h" // for re-use of |currentKL| and such

#include "complexredgp.h"
#include "complexredgp_io.h"
#include "error.h"
#include "helpmode.h"
#include "interactive.h"
#include "io.h"
#include "ioutils.h"
#include "basic_io.h"
#include "filekl.h"

#include "dynkin.h"
#include "lietype.h"
#include "realredgp.h"
#include "realredgp_io.h"
#include "kgb.h"
#include "kgb_io.h"
#include "blocks.h"
#include "repr.h"
#include "block_io.h"
#include "kl.h"
#include "kl_io.h"
#include "wgraph.h"
#include "wgraph_io.h"
#include "test.h"

/****************************************************************************

  This file contains the commands defined in the "repr" mode of the program.
  This means that a real form and a representation parameter have been chosen

*****************************************************************************/

namespace atlas {

namespace commands {

  void repr_mode_entry() throw(EntryError);
  void repr_mode_exit();

  // functions for the predefined commands

  void small_kgb_f(); // not yet implemented
  void small_dual_kgb_f(); // not yet implemented
  void iblock_f();
  void nblock_f();
  void partial_block_f();
  void block_f();
  void blockorder_f();
  void blockwrite_f(); // not yet implemented
  void blockstabilizer_f(); // not yet implemented
  void blocktwist_f();
  void deform_f();
  void kl_f();
  void klbasis_f();
  void kllist_f();
  void primkl_f();
  void klwrite_f();
  void wgraph_f();
  void wcells_f();

  void type_f();
  void realform_f();
  void repr_f();

  // mode-local variables
  block_type state=noblock;
  BlockElt entry_z = blocks::UndefBlock;
  SubSystemWithGroup* sub=NULL;
  StandardRepr* sr=NULL;
  param_block* block_pointer=NULL; // block contains |KLContext| pointer
  wgraph::WGraph* WGr_pointer=NULL;


/*****************************************************************************

        Chapter I -- Functions declared in reprmode.h

******************************************************************************/


// Returns a |CommandNode| object that is constructed on first call.
CommandNode reprNode()
{
  CommandNode result("repr: ",repr_mode_entry,repr_mode_exit);
  result.add("type",type_f); // override
  result.add("realform",realform_f); // this one too
  // result.add("smallkgb",small_kgb_f);
  // result.add("smalldualkgb",small_dual_kgb_f);
  result.add("iblock",iblock_f);
  result.add("nblock",nblock_f);
  result.add("partial_block",partial_block_f);
  result.add("block",block_f);
  result.add("blockorder",blockorder_f);
  // result.add("blockwrite",blockwrite_f);
  // result.add("blockstabilizer",blockstabilizer_f);
  result.add("blocktwist",blocktwist_f);
  result.add("deform",deform_f);
  result.add("kl",kl_f);
  result.add("klbasis",klbasis_f);
  result.add("kllist",kllist_f);
  result.add("primkl",primkl_f);
  result.add("klwrite",klwrite_f);
  result.add("wcells",wcells_f);
  result.add("wgraph",wgraph_f);

  // add test commands
  test::addTestCommands(result,ReprmodeTag());

  return result;
}

param_block& current_param_block()
{
  if (state==noblock) // we have entered reprmode without setting block
  {
    block_pointer = // partial block default
      new non_integral_block(currentRepTable(),*sr);
    state=partial_block;
    entry_z = block_pointer->size()-1;
  }
  return *block_pointer;
}

const SubSystemWithGroup& currentSubSystem() { return *sub; }

const StandardRepr& currentStandardRepr() { return *sr; }


/****************************************************************************

        Chapter II -- The repr mode |CommandNode|

  One instance of |CommandNode| for the repr mode is created at the
  first call of |reprMode()|; further calls just return a reference to it.

*****************************************************************************/

/*
  Synopsis: attempts to set a real form and dual real form interactively.
  In case of failure, throws an InputError and returns.
*/
void repr_mode_entry() throw(EntryError)
{
  try
  {
    RealReductiveGroup& GR = currentRealGroup();

    Weight lambda_rho;
    RatWeight gamma(0);
    KGBElt x;

    sub = new SubSystemWithGroup
      (interactive::get_parameter(GR,x,lambda_rho,gamma));

    Permutation pi;

    std::cout << "Subsystem on dual side is ";
    if (sub->rank()==0)
      std::cout << "empty.\n";
    else
    {
      std::cout << "of type "
		<< dynkin::Lie_type(sub->cartanMatrix(),true,false,pi)
		<< ", with roots ";
      for (weyl::Generator s=0; s<sub->rank(); ++s)
	std::cout << sub->parent_nr_simple(pi[s])
		  << (s<sub->rank()-1 ? "," : ".\n");
    }

    sr = new
      StandardRepr(currentRepContext().sr(x,lambda_rho,gamma));
  }
  catch(error::InputError& e)
  {
    repr_mode_exit(); // clean up
    e("no parameter was set");
    throw EntryError();
  }
}

/*
  Reset the parameter, effectively re-entering repr mode. If the choice
  of a new parameter fails, the current parameter remains in force.
*/
void repr_f()
{
  try
  {
    RealReductiveGroup& GR = currentRealGroup();

    Weight lambda_rho;
    RatWeight gamma(0);
    KGBElt x;

    sub = new SubSystemWithGroup
      (interactive::get_parameter(GR,x,lambda_rho,gamma));

    Permutation pi;

    std::cout << "Subsystem on dual side is ";
    if (sub->rank()==0)
      std::cout << "empty.\n";
    else
    {
      std::cout << "of type "
		<< dynkin::Lie_type(sub->cartanMatrix(),true,false,pi)
		<< ", with roots ";
      for (weyl::Generator s=0; s<sub->rank(); ++s)
	std::cout << sub->parent_nr_simple(pi[s])
		  << (s<sub->rank()-1 ? "," : ".\n");
    }
    delete sr;
    sr = new
      StandardRepr(currentRepContext().sr(x,lambda_rho,gamma));
    delete block_pointer; block_pointer=NULL;
    delete WGr_pointer; WGr_pointer=NULL;
  }
  catch (error::InputError& e)
  {
    e("parameter not changed");
  }
}
/*
  Synopsis: destroys any local data, resoring NULL pointers
*/
void repr_mode_exit()
{
  state=noblock;
  delete sr; sr=NULL;
  delete block_pointer; block_pointer=NULL;
  delete WGr_pointer; WGr_pointer=NULL;
}


/*****************************************************************************

        Chapter III --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/

/*
  Synopsis: resets the type of the complex group.

  In case of success, the real forms are invalidated, and therefore we
  should exit real mode; in case of failure, we don't need to.
*/
void type_f()
{
  try {
    ComplexReductiveGroup* G;
    complexredgp_io::Interface* I;

    interactive::getInteractive(G,I);
    replaceComplexGroup(G,I);
    exitMode(); // upon success pop block mode, destroying dual group
    exitMode(); // and pop real mode, destroying real group
  }
  catch (error::InputError& e) {
    e("complex group and real form not changed");
  }
}


/*
  Synopsis: resets the type, effectively re-entering the real mode. If the
  construction of the new type fails, the current block remains in force.
*/
void realform_f()
{
  try
  { // we can call the swap method for rvalues, but not with and rvalue arg
    interactive::getRealGroup(currentComplexInterface()).swap
      (currentRealGroup());

    exitMode(); // upon success pop repr mode, destroying data
  }
  catch (error::InputError& e) {
    e("real form not changed");
  }
}

void iblock_f()
{
  if (state!=iblock)
  {
    delete WGr_pointer; WGr_pointer=NULL;
    delete block_pointer; // destroy any installed block first
    block_pointer =
      new blocks::gamma_block(currentRepContext(),
			      currentSubSystem(),
			      currentStandardRepr(),
			      entry_z);
    state=iblock;
  }
  block_f();
} // |iblock_f|

void nblock_f()
{
  if (state!=nblock)
  {
    delete WGr_pointer; WGr_pointer=NULL;
    delete block_pointer; // destroy installed block first
    block_pointer =
      new non_integral_block(currentRepContext(),
			     currentStandardRepr(),
			     entry_z);
    state=nblock;
  }
  block_f();
} // |nblock_f|

void partial_block_f()
{
  if (state!=partial_block)
  {
    delete WGr_pointer; WGr_pointer=NULL;
    delete block_pointer; // destroy installed block first
    block_pointer =
      new non_integral_block(currentRepContext(),
			     currentStandardRepr());
    state=partial_block;
    entry_z = current_param_block().size()-1;
  }
  block_f();
} // |partial_block_f|

// Print the current block
void block_f()
{
  ioutils::OutputFile file;
  current_param_block().print_to(file,false);
  file << "Input parameters define element " << entry_z
       << " of this block." << std::endl;
}

// Print the Hasse diagram for the Bruhat order on the current block
void blockorder_f()
{
  param_block& block = current_param_block();
  std::cout << "block size: " << block.size() << std::endl;
  ioutils::OutputFile file;
  kgb_io::printBruhatOrder(file,block.bruhatOrder());
}

void blocktwist_f()
{
  ioutils::OutputFile file;
  block_io::print_twist(file,current_param_block());
}

void kl_f()
{
  ioutils::OutputFile file;
  param_block& block=current_param_block(); // now |entry_z| is defined
  block_io::print_KL(file,block,entry_z);
}


void deform_f()
{

  Rep_table& rt = currentRepTable();
  param_block& block = current_param_block();
  repr::SR_poly terms = rt.deformation_terms(block,entry_z);

  std::vector<StandardRepr> pool;
  HashTable<StandardRepr,unsigned long> hash(pool);

  ioutils::OutputFile f;

  f << "Orientation numbers:\n";
  bool first=true;
  for (BlockElt x=0; x<=entry_z; ++x)
    if (block.survives(x))
    {
      hash.match(rt.sr(block,x));
      if (first) first=false;
      else f<< ", ";
      StandardRepr r = rt.sr(block,x);
      f << x << ": " <<  rt.orientation_number(r);
    }
  f << ".\n";

  if (block.survives(entry_z))
  {
    f << "Deformation terms for I(" << entry_z << ")_c: (1-s) times\n";
    std::ostringstream os;
    for (repr::SR_poly::const_iterator it=terms.begin(); it!=terms.end(); ++it)
    {
      int eval=it->second.e();
      os << ' ';
      if (eval==1 or eval==-1)
	os << (eval==1 ? '+' : '-'); // sign of evaluation
      else
	os << std::setiosflags(std::ios_base::showpos) << eval;
      os <<"I(" << hash.find(it->first) << ")_c";
    }
    ioutils::foldLine(f,os.str()) << std::endl;

  }
} // |deform_f|


/* For each element $y$ in the block, outputs the list of non-zero K-L
   polynomials $P_{x,y}$.

   This is what is required to write down the K-L basis element $c_y$.
*/
void klbasis_f()
{
  const kl::KLContext& klc = currentKL();

  ioutils::OutputFile file;
  file << "Full list of non-zero Kazhdan-Lusztig-Vogan polynomials:"
       << std::endl << std::endl;
  kl_io::printAllKL(file,klc,current_param_block());
}


// Print the list of all distinct Kazhdan-Lusztig-Vogan polynomials
void kllist_f()
{
  const kl::KLContext& klc = currentKL();

  ioutils::OutputFile file;
  kl_io::printKLList(file,klc);
}

/*
  Print out the list of all K-L polynomials for primitive pairs.

  Explanation: x is primitive w.r.t. y, if any descent for y is also a
  descent for x, or a type II imaginary ascent. Ths means that none of
  the easy recursion formulas applies to P_{x,y}.
*/

void primkl_f()
{
  const kl::KLContext& klc = currentKL();

  ioutils::OutputFile file;
  file << "Kazhdan-Lusztig-Vogan polynomials for primitive pairs:"
       << std::endl << std::endl;
  kl_io::printPrimitiveKL(file,klc,current_param_block());
}

// Write the results of the KL computations to a pair of binary files
void klwrite_f()
{
  std::ofstream matrix_out, coefficient_out; // binary output files
  interactive::open_binary_file(matrix_out,"File name for matrix output: ");
  interactive::open_binary_file
    (coefficient_out,"File name for polynomial output: ");

  const kl::KLContext& klc = currentKL();

  if (matrix_out.is_open())
  {
    std::cout << "Writing matrix entries... " << std::flush;
    filekl::write_matrix_file(klc,matrix_out);
    std::cout << "Done." << std::endl;
  }
  if (coefficient_out.is_open())
  {
    std::cout << "Writing polynomial coefficients... " << std::flush;
    filekl::write_KL_store(klc.polStore(),coefficient_out);
    std::cout << "Done." << std::endl;
  }
}

// Print the W-graph corresponding to a block.
void wgraph_f()
{
  const wgraph::WGraph& wg = currentWGraph();
  ioutils::OutputFile file; wgraph_io::printWGraph(file,wg);
}

// Print the cells of the W-graph of the block.
void wcells_f()
{
  const wgraph::WGraph& wg = currentWGraph();
  wgraph::DecomposedWGraph dg(wg);

  ioutils::OutputFile file; wgraph_io::printWDecomposition(file,dg);
}



//      Chapter IV ---    H E L P    F U N C T I O N S


// Install help functions for block functions into help mode

namespace {

  const char* small_kgb_tag =
    "prints part of the KGB data pertinent to the current block";
  const char* small_dual_kgb_tag =
    "prints part of the dual KGB data pertinent to the current block";
  const char* block_tag = "prints all the representations in a block";
  const char* blockorder_tag =
   "shows Hasse diagram of the Bruhat order on the blocks";
  const char* blockwrite_tag = "writes the block information to disk";
  const char* blockstabilizer_tag = "print the real Weyl group for the block";
  const char* klbasis_tag = "prints the KL basis for the Hecke module";
  const char* kllist_tag = "prints the list of distinct KL polynomials";
  const char* klprim_tag = "prints the KL polynomials for primitive pairs";
  const char* klwrite_tag = "writes the KL polynomials to disk";
  const char* wgraph_tag = "prints the W-graph for the block";
  const char* wcells_tag = "prints the Kazhdan-Lusztig cells for the block";

void small_kgb_h()
{
  io::printFile(std::cerr,"smallkgb.help",io::MESSAGE_DIR);
}

void small_dual_kgb_h()
{
  io::printFile(std::cerr,"smalldualkgb.help",io::MESSAGE_DIR);
}

void block_h()
{
  io::printFile(std::cerr,"block.help",io::MESSAGE_DIR);
}

void blockorder_h()
{
  io::printFile(std::cerr,"blockorder.help",io::MESSAGE_DIR);
}

void blockwrite_h()
{
  io::printFile(std::cerr,"blockwrite.help",io::MESSAGE_DIR);
}

void block_stabilizer_h()
{
  io::printFile(std::cerr,"blockstabilizer.help",io::MESSAGE_DIR);
}

void klbasis_h()
{
  io::printFile(std::cerr,"klbasis.help",io::MESSAGE_DIR);
}

void kllist_h()
{
  io::printFile(std::cerr,"kllist.help",io::MESSAGE_DIR);
}

void primkl_h()
{
  io::printFile(std::cerr,"primkl.help",io::MESSAGE_DIR);
}

void klwrite_h()
{
  io::printFile(std::cerr,"klwrite.help",io::MESSAGE_DIR);
}

void wcells_h()
{
  io::printFile(std::cerr,"wcells.help",io::MESSAGE_DIR);
}

void wgraph_h()
{
  io::printFile(std::cerr,"wgraph.help",io::MESSAGE_DIR);
}



void addBlockHelp(CommandNode& mode, TagDict& tagDict)

{
  mode.add("smallkgb",small_kgb_h);
  mode.add("smalldualkgb",small_dual_kgb_h);
  mode.add("block",block_h);
  mode.add("blockorder",blockorder_h);
  mode.add("klwrite",klwrite_h);
  mode.add("kllist",kllist_h);
  mode.add("primkl",primkl_h);
  mode.add("klbasis",klbasis_h);
  mode.add("wcells",wcells_h);
  mode.add("wgraph",wgraph_h);

  // using commands::insertTag; // argument-dependent lookup will find it

  // insertTag(tagDict,"components",components_tag);
  insertTag(tagDict,"smallkgb",small_kgb_tag);
  insertTag(tagDict,"smalldualkgb",small_dual_kgb_tag);
  insertTag(tagDict,"block",block_tag);
  insertTag(tagDict,"blockorder",blockorder_tag);
  insertTag(tagDict,"blockwrite",blockwrite_tag);
  insertTag(tagDict,"blockstabilizer",blockstabilizer_tag);
  insertTag(tagDict,"klbasis",klbasis_tag);
  insertTag(tagDict,"kllist",kllist_tag);
  insertTag(tagDict,"primkl",klprim_tag);
  insertTag(tagDict,"klwrite",klwrite_tag);
  insertTag(tagDict,"wcells",wcells_tag);
  insertTag(tagDict,"wgraph",wgraph_tag);
}


} // |namespace|

} // |namespace commands|

} // |namespace atlas|
