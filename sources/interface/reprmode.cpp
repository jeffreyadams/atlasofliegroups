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

  void nblock_f();
  void partial_block_f();
  void block_f();
  void blockorder_f();
  void blocktwist_f();
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
  param_block* param_block_pointer=nullptr; // block gives access to |KLContext|
  wgraph::WGraph* param_WGr_pointer=nullptr;


/*****************************************************************************

        Chapter I -- Functions declared in reprmode.h

******************************************************************************/


// Returns a |CommandNode| object that is constructed on first call.
CommandNode reprNode()
{
  CommandNode result("repr: ",repr_mode_entry,repr_mode_exit);

  result.add("repr",repr_f,"override");
  result.add("nblock",nblock_f,"computes a non-integral block",std_help);
  result.add("partial_block",partial_block_f,
	     "computes the part of a non-integral block below given parameter",
	     use_tag);
  result.add("block",block_f,"second"); // block mode sets tag
  result.add("blockorder",blockorder_f,"second");
  result.add("blocktwist",blocktwist_f,"second");
  result.add("extblock",extblock_f,"second");
  result.add("gextblock",gextblock_f,"second");
  result.add("deform",deform_f,"computes deformation terms",std_help);
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

param_block& current_param_block()
{
  if (state==noblock) // we have entered reprmode without setting block
  {
    param_block_pointer = // partial block default
      new param_block(currentRepTable(),*sr);
    state=partial_block;
    entry_z = param_block_pointer->size()-1;
  }
  return *param_block_pointer;
}

const SubSystemWithGroup& currentSubSystem() { return *sub; }

const StandardRepr& currentStandardRepr() { return *sr; }

kl::KLContext& current_param_KL()
{
  return current_param_block().klc(current_param_block().size()-1,true);
}

void ensure_full_block()
{
  if (state!=nblock)
  {
    delete param_block_pointer; // destroy installed block first
    param_block_pointer =
      new param_block(currentRepContext(),currentStandardRepr(),entry_z);
    state=nblock;
  }
}

const wgraph::WGraph& current_param_WGraph()
{ if (param_WGr_pointer==nullptr)
  { ensure_full_block();
    const kl::KLContext& c=current_param_KL();
    param_WGr_pointer=new wgraph::WGraph(kl::wGraph(c));
  }
  return *param_WGr_pointer;
}

/****************************************************************************

        Chapter II -- The repr mode |CommandNode|

  One instance of |CommandNode| for the repr mode is created at the
  first call of |reprMode()|; further calls just return a reference to it.

*****************************************************************************/

/*
  Synopsis: attempts to set a real form and dual real form interactively.
  In case of failure, throws an InputError and returns.
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
    delete param_block_pointer; param_block_pointer=nullptr;
    delete param_WGr_pointer; param_WGr_pointer=nullptr;
    drop_to(repr_mode); // exit from (hypothetical) descendant modes
  }
  catch (error::InputError& e)
  {
    e("parameter not changed");
  }
}
/*
  Synopsis: destroys any local data, resoring nullptr pointers
*/
void repr_mode_exit()
{
  state=noblock;
  delete sr; sr=nullptr;
  delete param_block_pointer; param_block_pointer=nullptr;
  delete param_WGr_pointer; param_WGr_pointer=nullptr;
}


/*****************************************************************************

        Chapter III --- Functions for the predefined commands

  This section contains the definitions of the functions associated to the
  various commands defined in this mode.

******************************************************************************/

void nblock_f()
{
  ensure_full_block();
  block_f();
} // |nblock_f|

void partial_block_f()
{
  if (state!=partial_block)
  {
    delete param_WGr_pointer; param_WGr_pointer=nullptr;
    delete param_block_pointer; // destroy installed block first
    param_block_pointer =
      new param_block(currentRepContext(),currentStandardRepr());
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
    (current_param_block(),
     current_inner_class().distinguished());
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
  auto& block = current_param_block();
  ext_block::ext_block eblock(block,delta,true);
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
  auto& block = current_param_block();
  ext_block::ext_block eblock(block,delta,true);
  std::vector<ext_kl::Pol> pool;
  ext_kl::KL_table twisted_KLV(eblock,pool);
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
      hash.match(block.sr(x));
      if (first) first=false;
      else f<< ", ";
      StandardRepr r = block.sr(x);
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
  const kl::KLContext& klc = current_param_KL();

  ioutils::OutputFile file;
  file << "Full list of non-zero Kazhdan-Lusztig-Vogan polynomials:"
       << std::endl << std::endl;
  kl_io::printAllKL(file,klc,current_param_block());
}


// Print the list of all distinct Kazhdan-Lusztig-Vogan polynomials
void kllist_f()
{
  const kl::KLContext& klc = current_param_KL();

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
  const kl::KLContext& klc = current_param_KL();

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

  const kl::KLContext& klc = current_param_KL();

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
  const wgraph::WGraph& wg = current_param_WGraph();
  ioutils::OutputFile file; wgraph_io::printWGraph(file,wg);
}

// Print the cells of the W-graph of the block.
void wcells_f()
{
  const wgraph::WGraph& wg = current_param_WGraph();
  wgraph::DecomposedWGraph dg(wg);

  ioutils::OutputFile file; wgraph_io::printWDecomposition(file,dg);
}


} // |namespace commands|

} // |namespace atlas|
