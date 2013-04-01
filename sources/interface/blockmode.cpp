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

#include "dynkin.h"
#include "lietype.h"
#include "realredgp.h"
#include "realredgp_io.h"
#include "kgb.h"
#include "kgb_io.h"
#include "blocks.h"
#include "ext_block.h"
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

namespace commands {

namespace {

  void block_mode_entry() throw(EntryError);
  void block_mode_exit();

  // functions for the predefined commands

  void smallkgb_f();
  void smalldualkgb_f();
  void block_f();
  void smallblock_f();
  void dualblock_f();
  void smalldualblock_f();
  void dualmap_f();
  void blockd_f();
  void blocku_f();
  void blockorder_f();
  void blockwrite_f();
  void blockstabilizer_f();
  void blocktwist_f();
  void extblock_f();
  void klbasis_f();
  void kllist_f();
  void primkl_f();
  void klwrite_f();
  void wgraph_f();
  void wcells_f();

  void dualrealform_f();

  void smallkgb_h();
  void smalldualkgb_h();
  void block_h();
  void smallblock_h();
  void dualblock_h();
  void smalldualblock_h();
  void dualmap_h();
  void blockd_h();
  void blocku_h();
  void blockorder_h();
  void blockwrite_h();
  void blockstabilizer_h();
  void blocktwist_h();
  void extblock_h();
  void klbasis_h();
  void kllist_h();
  void primkl_h();
  void klwrite_h();
  void wgraph_h();
  void wcells_h();


  const char* smallkgb_tag =
    "prints part of the KGB data pertinent to one block";
  const char* smalldualkgb_tag =
    "prints part of the dual KGB data pertinent to one block";
  const char* block_tag = "prints all the representations in a block";
  const char* smallblock_tag =
    "generates block using partial KGB and dual KGB data";
  const char* dualblock_tag = "prints a block for the dual group";
  const char* smalldualblock_tag =
    "generates dual block using partial KGB and dual KGB data";
  const char* dualmap_tag =
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
  const char* primkl_tag = "prints the KL polynomials for primitive pairs";
  const char* klwrite_tag = "writes the KL polynomials to disk";
  const char* wgraph_tag = "prints the W-graph for the block";
  const char* wcells_tag = "prints the Kazhdan-Lusztig cells for the block";

  // local variables

  ComplexReductiveGroup* dual_G_C_pointer=NULL;
  RealReductiveGroup* dual_G_R_pointer=NULL;
  Block* block_pointer=NULL;
  wgraph::WGraph* WGr_pointer=NULL;
} // |namespace|

/*****************************************************************************

        Chapter I -- Functions declared in blockmode.h

******************************************************************************/

// Returns a |CommandNode| object that is constructed on first call.
CommandNode blockNode()
{
  CommandNode result("block: ",block_mode_entry,block_mode_exit);

  result.add("dualrealform",dualrealform_f,"override");
  result.add("smallkgb",smallkgb_f,smallkgb_tag,smallkgb_h);
  result.add("smalldualkgb",smalldualkgb_f,smalldualkgb_tag,smalldualkgb_h);
  result.add("block",block_f,block_tag,block_h);
  result.add("smallblock",smallblock_f,smallblock_tag,smallblock_h);
  result.add("dualblock",dualblock_f,dualblock_tag,dualblock_h);
  result.add("smalldualblock",smalldualblock_f,
	     smalldualblock_tag,smalldualblock_h);
  result.add("dualmap",dualmap_f,dualmap_tag,dualmap_h);
  result.add("blockd",blockd_f,blockd_tag,blockd_h);
  result.add("blocku",blocku_f,blocku_tag,blocku_h);
  result.add("blockorder",blockorder_f,blockorder_tag,blockorder_h);
  result.add("blockwrite",blockwrite_f,blockwrite_tag,blockwrite_h);
  result.add("blockstabilizer",blockstabilizer_f,
	     blockstabilizer_tag,blockstabilizer_h);
  result.add("blocktwist",blocktwist_f,"shows twist orbits on block");
  result.add("extblock",extblock_f,"prints block for extended group");
  result.add("klbasis",klbasis_f,klbasis_tag,klbasis_h);
  result.add("kllist",kllist_f,kllist_tag,kllist_h);
  result.add("primkl",primkl_f,primkl_tag,primkl_h);
  result.add("klwrite",klwrite_f,klwrite_tag,klwrite_h);
  result.add("wcells",wcells_f,wcells_tag,wcells_h);
  result.add("wgraph",wgraph_f,wgraph_tag,wgraph_h);

  // add test commands
  test::addTestCommands<BlockmodeTag>(result);

  return result;
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
    block_pointer =
      new Block(Block::build(currentRealGroup(),currentDualRealGroup()));

  return *block_pointer;
}

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

        Chapter II -- The block mode |CommandNode|

  One instance of |CommandNode| for the block mode is created at the
  first call of |blockMode()|; further calls just return a reference to it.

*****************************************************************************/

namespace {

/*
  Synopsis: attempts to set a real form and dual real form interactively.
  In case of failure, throws an InputError and returns.
*/
void block_mode_entry() throw(EntryError)
{
  try
  {
    RealReductiveGroup& G_R = currentRealGroup();

    ComplexReductiveGroup& G_C = G_R.complexGroup();
    const complexredgp_io::Interface& G_I = currentComplexInterface();

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
    throw EntryError();
  }
}

/*
  Reset the dual real form, effectively re-entering block mode. If the choice
  of a new dual real form fails, the current dual real form remains in force.
*/
void dualrealform_f()
{
  try
  {
    RealReductiveGroup& G_R = currentRealGroup();
    ComplexReductiveGroup& G_C = G_R.complexGroup();
    const complexredgp_io::Interface& G_I = currentComplexInterface();

    // get dual real form
    RealFormNbr drf;

    interactive::getInteractive
      (drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),tags::DualTag());

    // we can call the swap method for rvalues, but not with and rvalue arg
    RealReductiveGroup(*dual_G_C_pointer,drf).swap(*dual_G_R_pointer);

    delete block_pointer; block_pointer=NULL;
    delete WGr_pointer; WGr_pointer=NULL;
    drop_to(block_mode); // exit from (hypothetical) descendant modes
  }
  catch (error::InputError& e) {
    e("dual real form not changed");
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

/*****************************************************************************

        Chapter III --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/


// Print the kgb table, only the necessary part for one block
void smallkgb_f()
{
  RealReductiveGroup& G_R = currentRealGroup();
  RealReductiveGroup& dGR = currentDualRealGroup();

  BitMap common=blocks::common_Cartans(G_R,dGR);

  std::cout << "relevant Cartan classes: ";
  basic_io::seqPrint(std::cout,common.begin(),common.end(),",","{","}\n");

  std::cout
    << "partial kgb size: "
    << currentComplexGroup().KGB_size
         (currentRealForm(),common)
    << std::endl;

  ioutils::OutputFile file;
  KGB kgb(G_R,common);
  kgb_io::printKGB(file,kgb);
}

void smalldualkgb_f()
{
  RealReductiveGroup& G_R = currentRealGroup();
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
  Block::build(currentComplexGroup(),
	       currentRealForm(),
	       currentDualRealForm()).print_to(file,false);
}

// Print the dual block of the current block
void dualblock_f()
{
  Block block =
    Block::build(currentDualRealGroup(),currentRealGroup());

  ioutils::OutputFile file;
  block.print_to(file,false);
}

void smalldualblock_f()
{
  ComplexReductiveGroup& dG = currentDualComplexGroup();

  Block block =
    Block::build(dG,currentDualRealForm(),currentRealForm());

  ioutils::OutputFile file;
  block.print_to(file,false);
}

// Print the correspondence of the current block with its dual block
void dualmap_f()
{
  const Block& block = currentBlock();
  Block dual_block =
    Block::build(currentDualRealGroup(),currentRealGroup());

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
  RealReductiveGroup& G_R = currentRealGroup();
  RealReductiveGroup& dGR = currentDualRealGroup();

  // get Cartan class; abort if unvalid
  size_t cn=interactive::get_Cartan_class(blocks::common_Cartans(G_R,dGR));

  ioutils::OutputFile file;
  realredgp_io::printBlockStabilizer
    (file,currentRealGroup(),cn,currentDualRealForm());
}

void blocktwist_f()
{
  ioutils::OutputFile file;
  block_io::print_twist(file,currentBlock());
}

void extblock_f()
{
  ext_block::extended_block eblock(currentBlock(),
				   currentComplexGroup().twistedWeylGroup());
  ioutils::OutputFile file;
  eblock.print_to(file);
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
  const wgraph::DecomposedWGraph dg(wg);

  ioutils::OutputFile file; wgraph_io::printWDecomposition(file,dg);
}





//      Chapter IV ---    H E L P    F U N C T I O N S



void smallkgb_h()
{
  io::printFile(std::cerr,"smallkgb.help",io::MESSAGE_DIR);
}

void smalldualkgb_h()
{
  io::printFile(std::cerr,"smalldualkgb.help",io::MESSAGE_DIR);
}

void block_h()
{
  io::printFile(std::cerr,"block.help",io::MESSAGE_DIR);
}

void smallblock_h()
{
  io::printFile(std::cerr,"smallblock.help",io::MESSAGE_DIR);
}

void dualblock_h()
{
  io::printFile(std::cerr,"dualblock.help",io::MESSAGE_DIR);
}

void smalldualblock_h()
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

void blockstabilizer_h()
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

} // |namespace commands|

} // |namespace atlas|
