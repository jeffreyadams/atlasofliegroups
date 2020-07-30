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

#include "innerclass.h"
#include "output.h"
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
#include "kgb.h"
#include "kgb_io.h"
#include "blocks.h"
#include "ext_block.h"
#include "subsystem.h"
#include "repr.h"
#include "block_io.h"
#include "kl.h"
#include "kl_io.h"
#include "wgraph.h"
#include "wgraph_io.h"
#include "test.h"
#include "ext_kl.h"

/****************************************************************************

  This file contains the commands defined in the "repr" mode of the program.
  This means that a real form and a representation parameter have been chosen

*****************************************************************************/

namespace atlas {

namespace commands {

  void repr_mode_entry();
  void repr_mode_exit();

  // functions for the predefined commands

  void full_block_f();
  void block_f();
  void blockorder_f();
  void extblock_f();
  void gextblock_f();
  void deform_f();
  void kl_f();
  void extkl_f();
  void klbasis_f();
  void kllist_f();
  void primkl_f();
  void klwrite_f();
  void wgraph_f();
  void wcells_f();

  void repr_f(); // assign a new parameter while staying in this mode

  // mode-local variables
  block_type state=noblock;
  BlockElt entry_z = UndefBlock;
  SubSystemWithGroup* sub=nullptr;
  StandardRepr* sr=nullptr;
  blocks::common_block* common_block_pointer=nullptr; // accesses a |KL_table|
  wgraph::WGraph* WGr_pointer=nullptr;


/*****************************************************************************

        Chapter I -- Functions declared in reprmode.h

******************************************************************************/


// Returns a |CommandNode| object that is constructed on first call.
CommandNode reprNode()
{
  CommandNode result("repr: ",repr_mode_entry,repr_mode_exit);

  result.add("repr",repr_f,"override");
  result.add("full_block",full_block_f,"computes a full non-integral block",
	     std_help);
  result.add("block",block_f,"second"); // block mode sets tag
  result.add("blockorder",blockorder_f,"second");
  result.add("extblock",extblock_f,"second");
  result.add("gextblock",gextblock_f,"second");
  result.add("kl",kl_f,
	     "computes KL polynomials in character formula for this parameter",
	     std_help);
  result.add("klbasis",klbasis_f,"second");
  result.add("kllist",kllist_f,"second");
  result.add("primkl",primkl_f,"second");
  result.add("klwrite",klwrite_f,"second");
  result.add("wcells",wcells_f,"second");
  result.add("wgraph",wgraph_f,"second");
  result.add("extkl",extkl_f,"computes the KL polynomials for extended block");

  // add test commands
  test::addTestCommands<ReprmodeTag>(result);

  return result;
}

blocks::common_block& current_common_block()
{
  if (state!=full_block or common_block_pointer==nullptr)
  { // we have entered reprmode without setting block
    const Rep_context rc(currentRealGroup());
    auto srm = repr::StandardReprMod::mod_reduce(rc,*sr);
    common_block_pointer = // generate full block and set |entry_z|
      new blocks::common_block(currentRepTable(),srm,entry_z);
    state=full_block;
    entry_z = common_block_pointer->size()-1;
  }
  return *common_block_pointer;
}

const SubSystemWithGroup& currentSubSystem() { return *sub; }

const StandardRepr& currentStandardRepr() { return *sr; }

kl::KL_table& current_KL()
{
  return current_common_block().kl_tab(nullptr,0,true);
}

void ensure_full_block()
{
  if (state!=full_block)
    current_common_block(); // install (full) common block
}

const wgraph::WGraph& current_WGraph()
{ if (WGr_pointer==nullptr)
  { ensure_full_block();
    const kl::KL_table& c=current_KL();
    WGr_pointer=new wgraph::WGraph(kl::wGraph(c));
  }
  return *WGr_pointer;
}

/****************************************************************************

        Chapter II -- The repr mode |CommandNode|

  One instance of |CommandNode| for the repr mode is created at the
  first call of |reprMode()|; further calls just return a reference to it.

*****************************************************************************/

/*
  Attempt to set a real form and dual real form interactively.
  In case of failure catches an |InputError|, signals the user,
  and rethrows |EntryError|.
*/
void repr_mode_entry()
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
    state = noblock;
    delete common_block_pointer; common_block_pointer=nullptr;
    delete WGr_pointer; WGr_pointer=nullptr;
    drop_to(repr_mode); // exit from (hypothetical) descendant modes
  }
  catch (error::InputError& e)
  {
    e("parameter not changed");
  }
}

// Destroy any local data, restoring |nullptr| pointers
void repr_mode_exit()
{
  state=noblock;
  delete sr; sr=nullptr;
  delete common_block_pointer; common_block_pointer=nullptr;
  delete WGr_pointer;WGr_pointer=nullptr;
}


/*****************************************************************************

        Chapter III --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/

void full_block_f()
{
  ensure_full_block();
  block_f();
} // |full_block_f|

// Print the current block
void block_f()
{
  ioutils::OutputFile file;
  current_common_block().print_to(file,false);
  file << "Input parameters define element " << entry_z
       << " of this block." << std::endl;
}

// Print the Hasse diagram for the Bruhat order on the current block
void blockorder_f()
{
  auto& block = current_common_block();
  std::cout << "block size: " << block.size() << std::endl;
  ioutils::OutputFile file;
  kgb_io::printBruhatOrder(file,block.bruhatOrder());
}

void extblock_f()
{
  const auto& delta=current_inner_class().distinguished(); // implicit here
  if (not ((delta-1)*sr->gamma().numerator()).isZero())
  {
    std::cout << "Infinitesimal character " << sr->gamma()
	      <<" not fixed by distinguished involution." << std::endl;
    return;
  }
  ensure_full_block();
  ext_block::ext_block eblock
    (current_common_block(),
     current_inner_class().distinguished(),
     nullptr);
  ioutils::OutputFile file;
  eblock.print_to(file);
}

void gextblock_f()
{
  WeightInvolution delta = interactive::get_commuting_involution
    (commands::current_layout(), commands::current_lattice_basis());
  if (not ((delta-1)*sr->gamma().numerator()).isZero())
  {
    std::cout << "Chosen delta does not fix gamma=" << sr->gamma()
	      << " for the current block." << std::endl;
    return;
  }

  ensure_full_block();
  auto& block = current_common_block();
  ext_block::ext_block eblock(block,delta,nullptr);
  std::cout << "Extended block structure checked successfully." << std::endl;
  ioutils::OutputFile file;
  eblock.print_to(file);
}

void extkl_f()
{
  WeightInvolution delta = interactive::get_commuting_involution
    (commands::current_layout(), commands::current_lattice_basis());
  if (not ((delta-1)*sr->gamma().numerator()).isZero())
  {
    std::cout << "Chosen delta does not fix gamma=" << sr->gamma()
	      << " for the current block." << std::endl;
    return;
  }

  ensure_full_block();
  auto& block = current_common_block();
  ext_block::ext_block eblock(block,delta,nullptr);
  ext_kl::KL_table twisted_KLV(eblock,nullptr);
  twisted_KLV.fill_columns();

  ioutils::OutputFile f;
  for (BlockElt y=0; y<eblock.size(); ++y)
    for (BlockElt x=y+1; x-->0; )
      if (not twisted_KLV.P(x,y).isZero())
      {
	f << "P(" << eblock.z(x) << ',' << eblock.z(y) << ")=";
	f << twisted_KLV.P(x,y) << std::endl;
      }
}


void kl_f()
{
  ioutils::OutputFile file;
  auto& block=current_common_block(); // now |entry_z| is defined
  auto singular = block.singular(currentStandardRepr().gamma());
  block_io::print_KL(file,block,entry_z,singular);
}



/* For each element $y$ in the block, outputs the list of non-zero K-L
   polynomials $P_{x,y}$.

   This is what is required to write down the K-L basis element $c_y$.
*/
void klbasis_f()
{
  const kl::KL_table& kl_tab = current_KL();

  ioutils::OutputFile file;
  file << "Full list of non-zero Kazhdan-Lusztig-Vogan polynomials:"
       << std::endl << std::endl;
  kl_io::printAllKL(file,kl_tab,current_common_block());
}


// Print the list of all distinct Kazhdan-Lusztig-Vogan polynomials
void kllist_f()
{
  const kl::KL_table& kl_tab = current_KL();

  ioutils::OutputFile file;
  kl_io::printKLList(file,kl_tab);
}

/*
  Print out the list of all K-L polynomials for primitive pairs.

  Explanation: x is primitive w.r.t. y, if any descent for y is also a
  descent for x, or a type II imaginary ascent. Ths means that none of
  the easy recursion formulas applies to P_{x,y}.
*/

void primkl_f()
{
  const kl::KL_table& kl_tab = current_KL();

  ioutils::OutputFile file;
  file << "Kazhdan-Lusztig-Vogan polynomials for primitive pairs:"
       << std::endl << std::endl;
  kl_io::printPrimitiveKL(file,kl_tab,current_common_block());
}

// Write the results of the KL computations to a pair of binary files
void klwrite_f()
{
  std::ofstream matrix_out, coefficient_out; // binary output files
  interactive::open_binary_file(matrix_out,"File name for matrix output: ");
  interactive::open_binary_file
    (coefficient_out,"File name for polynomial output: ");

  const kl::KL_table& kl_tab = current_KL();

  if (matrix_out.is_open())
  {
    std::cout << "Writing matrix entries... " << std::flush;
    filekl::write_matrix_file(kl_tab,matrix_out);
    std::cout << "Done." << std::endl;
  }
  if (coefficient_out.is_open())
  {
    std::cout << "Writing polynomial coefficients... " << std::flush;
    filekl::write_KL_store(kl_tab.pol_store(),coefficient_out);
    std::cout << "Done." << std::endl;
  }
}

// Print the W-graph corresponding to a block.
void wgraph_f()
{
  const wgraph::WGraph& wg = current_WGraph();
  ioutils::OutputFile file; wgraph_io::printWGraph(file,wg);
}

// Print the cells of the W-graph of the block.
void wcells_f()
{
  const wgraph::WGraph& wg = current_WGraph();
  wgraph::DecomposedWGraph dg(wg);

  ioutils::OutputFile file; wgraph_io::printWDecomposition(file,dg);
}


} // |namespace commands|

} // |namespace atlas|
