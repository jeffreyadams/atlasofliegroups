/*
  This is reprmode.cpp

  Copyright (C) 2013 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#include <iostream>
#include <fstream>

#include "reprmode.h"
#include "realmode.h"
#include "mainmode.h"

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

namespace reprmode {

  void repr_mode_entry() throw(commands::EntryError);
  void repr_mode_exit();

  // functions for the predefined commands

  void small_kgb_f(); // not yet implemented
  void small_dual_kgb_f(); // not yet implemented
  void block_f();
  void blockorder_f();
  void blockwrite_f(); // not yet implemented
  void blockstabilizer_f(); // not yet implemented
  void blocktwist_f();
  void kl_f();
  void klbasis_f();
  void kllist_f();
  void primkl_f();
  void klwrite_f();
  void wgraph_f();
  void wcells_f();

  void type_f();
  void realform_f();

  // mode-local variables
  block_type state=noblock;
  BlockElt entry_z = blocks::UndefBlock;
  Rep_context* rc=NULL;
  SubSystemWithGroup* sub=NULL;
  StandardRepr* sr=NULL;
  param_block* block_pointer=NULL; // block contains |KLContext| pointer
  wgraph::WGraph* WGr_pointer=NULL;


/*****************************************************************************

        Chapter I -- Functions declared in reprmode.h

******************************************************************************/


// Returns a |CommandMode| object that is constructed on first call.
commands::CommandMode& reprMode()
{
  static commands::CommandMode repr_mode
    ("repr: ",repr_mode_entry,repr_mode_exit);
  if (repr_mode.empty()) // true upon first call
  {
    // add the commands from the real mode
    commands::addCommands(repr_mode,realmode::realMode());

    // add commands for this mode
    // the "type" command should be redefined here because it needs to exit
    // the block and real modes
    repr_mode.add("type",type_f); // override
    repr_mode.add("realform",realform_f); // this one too
    // repr_mode.add("smallkgb",small_kgb_f);
    // repr_mode.add("smalldualkgb",small_dual_kgb_f);
    // repr_mode.add("block",block_f);
    repr_mode.add("blockorder",blockorder_f);
    // repr_mode.add("blockwrite",blockwrite_f);
    // repr_mode.add("blockstabilizer",blockstabilizer_f);
    repr_mode.add("blocktwist",blocktwist_f);
    repr_mode.add("kl",kl_f);
    repr_mode.add("klbasis",klbasis_f);
    repr_mode.add("kllist",kllist_f);
    repr_mode.add("primkl",primkl_f);
    repr_mode.add("klwrite",klwrite_f);
    repr_mode.add("wcells",wcells_f);
    repr_mode.add("wgraph",wgraph_f);

    // add test commands
    test::addTestCommands(repr_mode,ReprmodeTag());
  }
  return repr_mode;
}

param_block& currentBlock()
{
  if (state==noblock) // we have entered reprmode without setting block
  {
    block_pointer =  new non_integral_block(*rc,*sr); // partial block default
    state=partial_block;
    entry_z = block_pointer->size()-1;
  }
  return *block_pointer;
}

const Rep_context& currentRepContext() { return *rc; }

const SubSystemWithGroup& currentSubSystem() { return *sub; }

const StandardRepr& currentStandardRepr() { return *sr; }

kl::KLContext& currentKL()
{
  return currentBlock().klc(currentBlock().size()-1,true);
}

const wgraph::WGraph& currentWGraph()
{
  if (WGr_pointer==NULL)
  {
    const kl::KLContext& c=currentKL();
    WGr_pointer=new wgraph::WGraph(c.rank());
    kl::wGraph(*WGr_pointer,c);
  }
  return *WGr_pointer;
}


/****************************************************************************

        Chapter II -- The repr mode |CommandMode|

  One instance of |CommandMode| for the repr mode is created at the
  first call of |reprMode()|; further calls just return a reference to it.

*****************************************************************************/

/*
  Synopsis: attempts to set a real form and dual real form interactively.
  In case of failure, throws an InputError and returns.
*/
void repr_mode_entry() throw(commands::EntryError)
{
  try
  {
    RealReductiveGroup& GR = realmode::currentRealGroup();

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

    rc = new Rep_context(GR);
    sr = new StandardRepr(rc->sr(x,lambda_rho,gamma));
  }
  catch(error::InputError& e)
  {
    repr_mode_exit(); // clean up
    e("no parameter was set");
    throw commands::EntryError();
  }
}


/*
  Synopsis: destroys any local data, resoring NULL pointers
*/
void repr_mode_exit()
{
  state=noblock;
  delete rc; rc=NULL;
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
    mainmode::replaceComplexGroup(G,I);
    commands::exitMode(); // upon success pop block mode, destroying dual group
    commands::exitMode(); // and pop real mode, destroying real group
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
    interactive::getRealGroup(mainmode::currentComplexInterface()).swap
      (realmode::currentRealGroup());

    commands::exitMode(); // upon success pop repr mode, destroying data
  }
  catch (error::InputError& e) {
    e("real form not changed");
  }
}


// Print the current block
void block_f()
{
  ioutils::OutputFile file;
  currentBlock().print_to(file,false);
}

// Print the Hasse diagram for the Bruhat order on the current block
void blockorder_f()
{
  param_block& block = currentBlock();
  std::cout << "block size: " << block.size() << std::endl;
  ioutils::OutputFile file;
  kgb_io::printBruhatOrder(file,block.bruhatOrder());
}

void blocktwist_f()
{
  ioutils::OutputFile file;
  block_io::print_twist(file,currentBlock());
}

void kl_f()
{
  ioutils::OutputFile file;
  block_io::print_KL(file,currentBlock(),entry_z);
}

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
  kl_io::printAllKL(file,klc,currentBlock());
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
  kl_io::printPrimitiveKL(file,klc,currentBlock());
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




} // namespace




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



void addBlockHelp(commands::CommandMode& mode, commands::TagDict& tagDict)

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

} // namespace reprmode

} // namespace atlas
