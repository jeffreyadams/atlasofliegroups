/*
  This is blockmode.cpp

  Copyright (C) 2007 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For copyright and license information see the LICENSE file
*/

#include <iostream>
#include <fstream>

#include "blockmode.h"
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
#include "realredgp.h"
#include "realredgp_io.h"
#include "kgb.h"
#include "kgb_io.h"
#include "blocks.h"
#include "block_io.h"
#include "kl.h"
#include "kl_io.h"
#include "wgraph.h"
#include "wgraph_io.h"
#include "test.h"

/****************************************************************************

  This file contains the commands defined in the "block" mode of the program.
  This means that a real form and a dual real form have been chosen.

*****************************************************************************/

namespace atlas {

namespace {

  using namespace blockmode;

  void block_mode_entry() throw(commands::EntryError);
  void block_mode_exit();

  // functions for the predefined commands

  void small_kgb_f();
  void small_dual_kgb_f();
  void block_f();
  void smallblock_f();
  void dual_block_f();
  void small_dual_block_f();
  void dual_map_f();
  void blockd_f();
  void blocku_f();
  void blockorder_f();
  void blockwrite_f();
  void blockstabilizer_f();
  void blocktwist_f();
  void klbasis_f();
  void kllist_f();
  void primkl_f();
  void klwrite_f();
  void wgraph_f();
  void wcells_f();

  void type_f();
  void realform_f();

  // local variables

  ComplexReductiveGroup* dual_G_C_pointer=NULL;
  RealReductiveGroup* dual_G_R_pointer=NULL;
  Block* block_pointer=NULL;
  wgraph::WGraph* WGr_pointer=NULL;
}

/*****************************************************************************

        Chapter I -- Functions declared in blockmode.h

******************************************************************************/

namespace blockmode {

// Returns a |CommandMode| object that is constructed on first call.
commands::CommandMode& blockMode()
{
  static commands::CommandMode block_mode
    ("block: ",block_mode_entry,block_mode_exit);
  if (block_mode.empty()) // true upon first call
  {
    // add the commands from the real mode
    commands::addCommands(block_mode,realmode::realMode());

    // add commands for this mode
    // the "type" command should be redefined here because it needs to exit
    // the block and real modes
    block_mode.add("type",type_f); // override
    block_mode.add("realform",realform_f); // this one too
    block_mode.add("smallkgb",small_kgb_f);
    block_mode.add("smalldualkgb",small_dual_kgb_f);
    block_mode.add("block",block_f);
    block_mode.add("smallblock",smallblock_f);
    block_mode.add("dualblock",dual_block_f);
    block_mode.add("smalldualblock",small_dual_block_f);
    block_mode.add("dualmap",dual_map_f);
    block_mode.add("blockd",blockd_f);
    block_mode.add("blocku",blocku_f);
    block_mode.add("blockorder",blockorder_f);
    block_mode.add("blockwrite",blockwrite_f);
    block_mode.add("blockstabilizer",blockstabilizer_f);
    block_mode.add("blocktwist",blocktwist_f);
    block_mode.add("klbasis",klbasis_f);
    block_mode.add("kllist",kllist_f);
    block_mode.add("primkl",primkl_f);
    block_mode.add("klwrite",klwrite_f);
    block_mode.add("wcells",wcells_f);
    block_mode.add("wgraph",wgraph_f);

    // add test commands
    test::addTestCommands(block_mode,BlockmodeTag());
  }
  return block_mode;
}

ComplexReductiveGroup& currentDualComplexGroup()
{
  return *dual_G_C_pointer;
}

RealReductiveGroup& currentDualRealGroup()
{
  return *dual_G_R_pointer;
}

RealFormNbr currentDualRealForm()
{
  return dual_G_R_pointer->realForm();
}

Block& currentBlock()
{
  if (block_pointer==NULL)
  {
    block_pointer =
      new Block(Block::build
			 (realmode::currentRealGroup(),
			  currentDualRealGroup()));
  }
  return *block_pointer;
}

kl::KLContext& currentKL()
{
  return currentBlock().klc(currentBlock().size()-1,true);
}

wgraph::WGraph& currentWGraph()
{
  if (WGr_pointer==NULL)
  {
    const kl::KLContext& c=currentKL();
    WGr_pointer=new wgraph::WGraph(c.rank());
    kl::wGraph(*WGr_pointer,c);
  }
  return *WGr_pointer;
}

} // namespace blockmode

/****************************************************************************

        Chapter II -- The block mode |CommandMode|

  One instance of |CommandMode| for the block mode is created at the
  first call of |blockMode()|; further calls just return a reference to it.

*****************************************************************************/

namespace {

/*
  Synopsis: attempts to set a real form and dual real form interactively.
  In case of failure, throws an InputError and returns.
*/
void block_mode_entry() throw(commands::EntryError)
{
  try
  {
    RealReductiveGroup& G_R = realmode::currentRealGroup();

    ComplexReductiveGroup& G_C = G_R.complexGroup();
    const complexredgp_io::Interface& G_I = mainmode::currentComplexInterface();

    // get dual real form
    RealFormNbr drf;

    interactive::getInteractive
      (drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),tags::DualTag());

    dual_G_C_pointer=new
      ComplexReductiveGroup(G_C,tags::DualTag());
    dual_G_R_pointer=new RealReductiveGroup(*dual_G_C_pointer,drf);
  }
  catch(error::InputError& e)
  {
    block_mode_exit(); // clean up
    e("no dual real form set");
    throw commands::EntryError();
  }
}


/*
  Synopsis: destroys any local data, resoring NULL pointers
*/
void block_mode_exit()
{
  delete dual_G_C_pointer; dual_G_C_pointer=NULL;
  delete dual_G_R_pointer; dual_G_R_pointer=NULL;
  delete block_pointer; block_pointer=NULL;
  delete WGr_pointer; WGr_pointer=NULL;
}

} // namespace

/*****************************************************************************

        Chapter III --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/

namespace {

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

    commands::exitMode(); // upon success pop block mode, destroying dual group
  }
  catch (error::InputError& e) {
    e("real form not changed");
  }
}

// Print the kgb table, only the necessary part for one block
void small_kgb_f()
{
  RealReductiveGroup& G_R = realmode::currentRealGroup();
  RealReductiveGroup& dGR = currentDualRealGroup();

  BitMap common=blocks::common_Cartans(G_R,dGR);

  std::cout << "relevant Cartan classes: ";
  basic_io::seqPrint(std::cout,common.begin(),common.end(),",","{","}\n");

  std::cout
    << "partial kgb size: "
    << mainmode::currentComplexGroup().KGB_size
         (realmode::currentRealForm(),common)
    << std::endl;

  ioutils::OutputFile file;
  KGB kgb(G_R,common);
  kgb_io::printKGB(file,kgb);
}

void small_dual_kgb_f()
{
  RealReductiveGroup& G_R = realmode::currentRealGroup();
  RealReductiveGroup& dGR = currentDualRealGroup();
  ComplexReductiveGroup& dGC = currentDualComplexGroup();

  BitMap common=blocks::common_Cartans(dGR,G_R);

  std::cout << "relevant Cartan classes for dual group: ";
  basic_io::seqPrint(std::cout,common.begin(),common.end(),",","{","}\n");

  std::cout << "partial kgb size: " <<
    dGC.KGB_size(currentDualRealForm(),common) << std::endl;
  ioutils::OutputFile file;

  KGB kgb(dGR,common);
  kgb_io::printKGB(file,kgb);
}

// Print the current block
void block_f()
{
  ioutils::OutputFile file;
  currentBlock().print_to(file,false);
}

void smallblock_f()
{
  ioutils::OutputFile file;
  // must unfortunatly regenerate the block here
  Block::build(mainmode::currentComplexGroup(),
	       realmode::currentRealForm(),
	       currentDualRealForm()).print_to(file,false);
}

// Print the dual block of the current block
void dual_block_f()
{
  Block block =
    Block::build(currentDualRealGroup(),realmode::currentRealGroup());

  ioutils::OutputFile file;
  block.print_to(file,false);
}

void small_dual_block_f()
{
  ComplexReductiveGroup& dG = currentDualComplexGroup();

  Block block =
    Block::build(dG,currentDualRealForm(),realmode::currentRealForm());

  ioutils::OutputFile file;
  block.print_to(file,false);
}

// Print the correspondence of the current block with its dual block
void dual_map_f()
{
  const Block& block = currentBlock();
  Block dual_block =
    Block::build(currentDualRealGroup(),realmode::currentRealGroup());

  const std::vector<BlockElt> v=blocks::dual_map(block,dual_block);

  std::ostringstream s("");
  basic_io::seqPrint(s,v.begin(),v.end(),", ","[","]\n");
  ioutils::OutputFile file;
  foldLine(file,s.str()," ");
}

// Print the current block with involutions in involution-reduced form
void blockd_f()
{
  ioutils::OutputFile file;
  currentBlock().print_to(file,true);
}

// Print the unitary elements of the block.
void blocku_f()
{
  ioutils::OutputFile file;
  block_io::printBlockU(file,currentBlock());
}


// Print the Hasse diagram for the Bruhat order on the current block
void blockorder_f()
{
  Block& block = currentBlock();
  std::cout << "block size: " << block.size() << std::endl;
  ioutils::OutputFile file;
  kgb_io::printBruhatOrder(file,block.bruhatOrder());
}

// Writes a binary file containing descent sets and ascent sets for block
void blockwrite_f()
{
  std::ofstream block_out; // binary output files
  if (interactive::open_binary_file
      (block_out,"File name for block output: "))
  {
    filekl::write_block_file(currentBlock(),block_out);
    std::cout<< "Binary file written.\n" << std::endl;
  }
  else
    std::cout << "No file written.\n";
}

/*
  Synopsis: prints out information about the stabilizer of a representation
  under the cross action
*/
void blockstabilizer_f()
{
  RealReductiveGroup& G_R = realmode::currentRealGroup();
  RealReductiveGroup& dGR = currentDualRealGroup();

  // get Cartan class; abort if unvalid
  size_t cn=interactive::get_Cartan_class(blocks::common_Cartans(G_R,dGR));

  ioutils::OutputFile file;
  realredgp_io::printBlockStabilizer
    (file,realmode::currentRealGroup(),cn,currentDualRealForm());
}

void blocktwist_f()
{
  ioutils::OutputFile file;

  block_io::print_twist(file,currentBlock());
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
  wgraph::WGraph& wg = currentWGraph();
  ioutils::OutputFile file; wgraph_io::printWGraph(file,wg);
}

// Print the cells of the W-graph of the block.
void wcells_f()
{
  wgraph::WGraph& wg = currentWGraph();
  wgraph::DecomposedWGraph dg(wg);

  ioutils::OutputFile file; wgraph_io::printWDecomposition(file,dg);
}




} // namespace




//      Chapter IV ---    H E L P    F U N C T I O N S


// Install help functions for block functions into help mode

namespace {

  const char* small_kgb_tag =
    "prints part of the KGB data pertinent to one block";
  const char* small_dual_kgb_tag =
    "prints part of the dual KGB data pertinent to one block";
  const char* block_tag = "prints all the representations in a block";
  const char* small_block_tag =
    "generates block using partial KGB and dual KGB data";
  const char* dual_block_tag = "prints a block for the dual group";
  const char* small_dual_block_tag =
    "generates dual block using partial KGB and dual KGB data";
  const char* dual_map_tag =
    "prints the bijection from block to its dual block";
  const char* blockd_tag =
   "prints all representations in the block, alternative format";
  const char* blocku_tag =
   "prints the unitary representations in the block at rho";
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

void small_block_h()
{
  io::printFile(std::cerr,"smallblock.help",io::MESSAGE_DIR);
}

void dualblock_h()
{
  io::printFile(std::cerr,"dualblock.help",io::MESSAGE_DIR);
}

void small_dual_block_h()
{
  io::printFile(std::cerr,"smalldualblock.help",io::MESSAGE_DIR);
}

void dualmap_h()
{
  io::printFile(std::cerr,"dualmap.help",io::MESSAGE_DIR);
}

void blockd_h()
{
  io::printFile(std::cerr,"blockd.help",io::MESSAGE_DIR);
}

void blocku_h()
{
  io::printFile(std::cerr,"blocku.help",io::MESSAGE_DIR);
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

} // namespace

namespace blockmode {

void addBlockHelp(commands::CommandMode& mode, commands::TagDict& tagDict)

{
  mode.add("smallkgb",small_kgb_h);
  mode.add("smalldualkgb",small_dual_kgb_h);
  mode.add("block",block_h);
  mode.add("smallblock",small_block_h);
  mode.add("dualblock",dualblock_h);
  mode.add("smalldualblock",small_dual_block_h);
  mode.add("dualmap",dualmap_h);
  mode.add("blockd",blockd_h);
  mode.add("blocku",blocku_h);
  mode.add("blockorder",blockorder_h);
  mode.add("blockwrite",blockwrite_h);
  mode.add("blockstabilizer",block_stabilizer_h);
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
  insertTag(tagDict,"smallblock",small_block_tag);
  insertTag(tagDict,"dualblock",dual_block_tag);
  insertTag(tagDict,"smalldualblock",small_dual_block_tag);
  insertTag(tagDict,"blockd",blockd_tag);
  insertTag(tagDict,"blocku",blocku_tag);
  insertTag(tagDict,"blockorder",blockorder_tag);
  insertTag(tagDict,"blockwrite",blockwrite_tag);
  insertTag(tagDict,"blockstabilizer",blockstabilizer_tag);
  insertTag(tagDict,"dualblock",dual_block_tag);
  insertTag(tagDict,"dualmap",dual_map_tag);
  insertTag(tagDict,"klbasis",klbasis_tag);
  insertTag(tagDict,"kllist",kllist_tag);
  insertTag(tagDict,"primkl",klprim_tag);
  insertTag(tagDict,"klwrite",klwrite_tag);
  insertTag(tagDict,"wcells",wcells_tag);
  insertTag(tagDict,"wgraph",wgraph_tag);
}

} // namespace blockmode

} // namespace atlas
