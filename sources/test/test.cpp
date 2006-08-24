/*
  This is test.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "test.h"

#include <cassert>
#include <iostream>
#include <iterator>
#include <set>
#include <sstream>

#include "basic_io.h"
#include "bitmap.h"
#include "block_io.h"
#include "blocks.h"
#include "bruhat.h"
#include "cartan_io.h"
#include "commands.h"
#include "complexredgp.h"
#include "complexredgp_io.h"
#include "dynkin.h"
#include "emptymode.h"
#include "error.h"
#include "gradings.h"
#include "helpmode.h"
#include "input.h"
#include "interactive.h"
#include "involutions.h"
#include "io.h"
#include "ioutils.h"
#include "kl.h"
#include "kl_io.h"
#include "klsupport.h"
#include "kltest.h"
#include "kgb.h"
#include "kgb_io.h"
#include "lattice.h"
#include "mainmode.h"
#include "poset_io.h"
#include "prettyprint.h"
#include "realform.h"
#include "realform_io.h"
#include "realmode.h"
#include "realredgp_io.h"
#include "realweyl.h"
#include "realweyl_io.h"
#include "rootdata.h"
#include "size.h"
#include "smithnormal.h"
#include "tits.h"
#include "tori.h"
#include "weyl.h"
#include "wgraph.h"
#include "wgraph_io.h"

#include "testprint.h"
#include "testrun.h"

/*****************************************************************************

  This module contains some commands for testing the program.

******************************************************************************/

namespace atlas {

namespace {
  using namespace test;

  // functions for the test commands

  void block_f();
  void blockd_f();
  void blocku_f();
  void blockstabilizer_f();
  void checkbasept_f();
  void cmatrix_f();
  void corder_f();
  void components_f();
  void coroots_rootbasis_f();
  void extrkl_f();
  void kgb_f();
  void klbasis_f();
  void kllist_f();
  void poscoroots_rootbasis_f();
  void posroots_rootbasis_f();
  void wcells_f();
  void wgraph_f();
  void roots_rootbasis_f();
  void rootdatum_f();
  void test_f();

  // help functions

  void block_h();
  void blocku_h();
  void cmatrix_h();
  void extrkl_h();
  void kgb_h();
  void klbasis_h();
  void kllist_h();
  void wcells_h();
  void wgraph_h();

  // tags

  const char* test_tag = "(test command)";
  const char* block_tag = "prints all the representations in a block";
  const char* blocku_tag = "prints the unitary representations in the block at rho";
  const char* cmatrix_tag = "prints the Cartan matrix";
  const char* kgb_tag = "prints the orbits of K on G/B";
  const char* klbasis_tag = "prints the KL basis for the Hecke module";
  const char* kllist_tag = "prints the list of distinct KL polynomials";
  const char* wcells_tag = "prints the Kazhdan-Lusztig cells for the block";
  const char* wgraph_tag = "prints the W-graph for the block";
  
  enum TestMode {EmptyMode, MainMode, RealMode, numTestMode};
  const TestMode testMode = MainMode;

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

template<>
void addTestCommands<emptymode::EmptymodeTag> 
  (commands::CommandMode& mode, emptymode::EmptymodeTag)

/*
  Synopsis: adds to the empty mode the test commands that require that mode.

  NOTE: the mode still needs to be passed as an argument, because this is
  called during the construction of the mode!

  NOTE: for convenience, the "test" command is added to the mode that is
  flagged by the testMode constant; this is EmptyMode by default, but should
  be redefined to the correct value.
*/

{
  if (testMode == EmptyMode)
    mode.add("test",test_f);

  return;
}

template<> 
void addTestCommands<mainmode::MainmodeTag> 
  (commands::CommandMode& mode, mainmode::MainmodeTag)

/*
  Synopsis: adds to the main mode the test commands that require that mode.

  NOTE: the mode still needs to be passed as an argument, because this is
  called during the construction of the mode!

  NOTE: for convenience, the "test" command is added to the mode that is
  flagged by the testMode constant; this is EmptyMode by default, but should
  be redefined to the correct value.
*/

{
  if (testMode == MainMode)
    mode.add("test",test_f);

  // add additional commands here :

  mode.add("cmatrix",cmatrix_f);
  mode.add("coroots_rootbasis",coroots_rootbasis_f);
  mode.add("poscoroots_rootbasis",poscoroots_rootbasis_f);
  mode.add("posroots_rootbasis",posroots_rootbasis_f);
  mode.add("roots_rootbasis",roots_rootbasis_f);
  mode.add("rootdatum",rootdatum_f);

  return;
}

template<> 
void addTestCommands<realmode::RealmodeTag> 
  (commands::CommandMode& mode, realmode::RealmodeTag)

/*
  Synopsis: adds to the real mode the test commands that require that mode.

  NOTE: the mode still needs to be passed as an argument, because this is
  called during the construction of the mode!

  NOTE: for convenience, the "test" command is added to the mode that is
  flagged by the testMode constant; this is EmptyMode by default, but should
  be redefined to the correct value.
*/

{  
  if (testMode == RealMode)
    mode.add("test",test_f);

  mode.add("block",block_f);
  mode.add("blockd",blockd_f);
  mode.add("blocku",blocku_f);
  mode.add("blockstabilizer",blockstabilizer_f);
  mode.add("checkbasept",checkbasept_f);
  mode.add("components",components_f);
  mode.add("corder",corder_f);
  mode.add("extrkl",extrkl_f);
  mode.add("kgb",kgb_f);
  mode.add("klbasis",klbasis_f);
  mode.add("kllist",kllist_f);
  mode.add("wcells",wcells_f);
  mode.add("wgraph",wgraph_f);

  return;
}

template<> void addTestHelp<emptymode::EmptymodeTag>
             (commands::CommandMode& mode, commands::TagDict& t, 
	      emptymode::EmptymodeTag)

/*
  Synopsis: adds to the empty help mode the help commands for the test
  commands that require that mode.

  NOTE: the mode still needs to be passed as an argument, because this is
  called during the construction of the mode!


  NOTE: for convenience, the "test" command is added to the mode that is
  flagged by the testMode constant; this is EmptyMode by default, but should
  be redefined to the correct value.
*/

{  
  using namespace commands;
  using namespace helpmode;

  if (testMode == EmptyMode) {
    mode.add("test",nohelp_h);
    insertTag(t,"test",test_tag);
  }

  // add additional help commands here:

  /******** tags ************************************************************/

  // add additional command tags here:

}

template<> void addTestHelp<mainmode::MainmodeTag>
             (commands::CommandMode& mode, commands::TagDict& t,
	      mainmode::MainmodeTag)

/*
  Synopsis: adds to the main help mode the help commands for the test
  commands that require that mode.

  NOTE: the mode still needs to be passed as an argument, because this is
  called during the construction of the mode!

  NOTE: for convenience, the "test" command is added to the mode that is
  flagged by the testMode constant; this is EmptyMode by default, but should
  be redefined to the correct value.
*/

{
  using namespace commands;
  using namespace helpmode;

  if (testMode == MainMode) {
    mode.add("test",nohelp_h);
    insertTag(t,"test",test_tag);
  }

  // add additional help commands here:

  /******** tags ************************************************************/

  // add additional command tags here:

  mode.add("cmatrix",cmatrix_h);
  mode.add("coroots_rootbasis",nohelp_h);
  mode.add("gradings",nohelp_h);
  mode.add("poscoroots_rootbasis",nohelp_h);
  mode.add("posroots_rootbasis",nohelp_h);
  mode.add("roots_rootbasis",nohelp_h);
  mode.add("rootdatum",nohelp_h);

  // add additional command tags here :

  insertTag(t,"cmatrix",cmatrix_tag);
  insertTag(t,"coroots_rootbasis",test_tag);
  insertTag(t,"gradings",test_tag);
  insertTag(t,"poscoroots_rootbasis",test_tag);
  insertTag(t,"posroots_rootbasis",test_tag);
  insertTag(t,"roots_rootbasis",test_tag);
  insertTag(t,"rootdatum",test_tag);
 
  return;
}

template<> void addTestHelp<realmode::RealmodeTag>
             (commands::CommandMode& mode, commands::TagDict& t, 
	      realmode::RealmodeTag)

/*
  Synopsis: adds to the real help mode the help commands for the test
  commands that require that mode.

  NOTE: the mode still needs to be passed as an argument, because this is
  called during the construction of the mode!

  NOTE: for convenience, the "test" command is added to the mode that is
  flagged by the testMode constant; this is EmptyMode by default, but should
  be redefined to the correct value.
*/

{
  using namespace commands;
  using namespace helpmode;

  if (testMode == RealMode) {
    insertTag(t,"test",test_tag);
  }

  // add additional help commands here:

  mode.add("block",block_h);
  mode.add("blockd",nohelp_h);
  mode.add("blocku",blocku_h);
  mode.add("blockstabilizer",nohelp_h);
  mode.add("checkbasept",nohelp_h);
  mode.add("components",nohelp_h);
  mode.add("corder",nohelp_h);
  mode.add("extrkl",extrkl_h);
  mode.add("involution",nohelp_h);
  mode.add("kgb",kgb_h);
  mode.add("klbasis",klbasis_h);
  mode.add("kllist",kllist_h);
  mode.add("wcells",wcells_h);
  mode.add("wgraph",wgraph_h);

  /******** tags ************************************************************/

  // add additional command tags here:
  insertTag(t,"block",block_tag);
  insertTag(t,"blockd",test_tag);
  insertTag(t,"blocku",blocku_tag);
  insertTag(t,"blockstabilizer",test_tag);
  insertTag(t,"checkbasept",test_tag);
  insertTag(t,"components",test_tag);
  insertTag(t,"corder",test_tag);
  insertTag(t,"extrkl",test_tag);
  insertTag(t,"involution",test_tag);
  insertTag(t,"kgb",kgb_tag);
  insertTag(t,"klbasis",klbasis_tag);
  insertTag(t,"kllist",kllist_tag);
  insertTag(t,"wcells",wcells_tag);
  insertTag(t,"wgraph",wgraph_tag);

  return;
}

}

namespace {

void block_h()

{
  io::printFile(std::cerr,"block.help",io::MESSAGE_DIR);
  return;
}

void blocku_h()

{
  io::printFile(std::cerr,"blocku.help",io::MESSAGE_DIR);
  return;
}

void cmatrix_h()

{
  io::printFile(std::cerr,"cmatrix.help",io::MESSAGE_DIR);
  return;
}

void extrkl_h()

{
  io::printFile(std::cerr,"extrkl.help",io::MESSAGE_DIR);
  return;
}

void kgb_h()

{
  io::printFile(std::cerr,"kgb.help",io::MESSAGE_DIR);
  return;
}

void klbasis_h()

{
  io::printFile(std::cerr,"klbasis.help",io::MESSAGE_DIR);
  return;
}

void kllist_h()

{
  io::printFile(std::cerr,"kllist.help",io::MESSAGE_DIR);
  return;
}

void wcells_h()

{
  io::printFile(std::cerr,"wcells.help",io::MESSAGE_DIR);
  return;
}

void wgraph_h()

{
  io::printFile(std::cerr,"wgraph.help",io::MESSAGE_DIR);
  return;
}

}

/*****************************************************************************

        Chapter II -- Utility functions

******************************************************************************/

namespace {

const rootdata::RootDatum& currentRootDatum()

{
  return mainmode::currentComplexGroup().rootDatum();
}

}

/*****************************************************************************

        Chapter III -- Functions for the test commands

******************************************************************************/

namespace {

void block_f()

/*
  Synopsis: constructs the block of the category of Harish-Chandra modules
  corresponding to a given real form and dual real form.
*/

{
  using namespace block_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();
  }
  catch (MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
  const realredgp_io::Interface& G_RI = currentRealInterface();
  const complexredgp_io::Interface& G_I = G_RI.complexInterface();

  // get dual real form
  RealForm drf;

  try {
    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());
  }
  catch (InputError& e) {
    e("aborted");
    return;
  }

  Block block(G_C,G_R.realForm(),drf);

  OutputFile file;
  printBlock(file,block);

  return;
}

void blockd_f()

/*
  Synopsis: constructs the block of the category of Harish-Chandra modules
  corresponding to a given real form and dual real form.

  NOTE: the difference with the ordinary blocks is that it outputs involutions
  in reduced-involution form.
*/

{
  using namespace block_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();
  }
  catch (MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
  const realredgp_io::Interface& G_RI = currentRealInterface();
  const complexredgp_io::Interface& G_I = G_RI.complexInterface();

  // get dual real form
  RealForm drf;

  try {
    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());
  }
  catch (InputError& e) {
    e("aborted");
    return;
  }

  Block block(G_C,G_R.realForm(),drf);

  OutputFile file;
  printBlockD(file,block);

  return;
}

void blocku_f()

/*
  Synopsis: outputs the _unitary_ elements of the block.
*/

{
  using namespace block_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();
  }
  catch (MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
  const realredgp_io::Interface& G_RI = currentRealInterface();
  const complexredgp_io::Interface& G_I = G_RI.complexInterface();

  // get dual real form
  RealForm drf;

  try {
    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());
  }
  catch (InputError& e) {
    e("aborted");
    return;
  }

  Block block(G_C,G_R.realForm(),drf);

  OutputFile file;
  printBlockU(file,block);

  return;
}

void blockstabilizer_f()

/*
  Synopsis: prints out information about the stabilizer of a representation 
  under the cross action
*/

{
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();
  }
  catch (MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  size_t cn;

  // get Cartan class; abort if unvalid
  try {
    getCartanClass(cn,G_R.cartanSet(),currentLine());
  }
  catch (InputError& e) {
    e("aborted");
    return;
  }

  const complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
  const realredgp_io::Interface& G_RI = currentRealInterface();
  const complexredgp_io::Interface& G_I = G_RI.complexInterface();

  // get dual real form
  RealForm drf;

  try {
    getInteractive(drf,G_I,G_C.dualRealFormLabels(cn),DualTag());
  }
  catch (InputError& e) {
    e("aborted");
    return;
  }

  OutputFile file;
  printBlockStabilizer(file,G_RI,cn,drf);

  return;
}

void checkbasept_f()

/*
  Synopsis: checks if the conjectural basepoint construction is ok.
*/

{
  using namespace kgb;
  using namespace kltest;
  using namespace realmode;
  using namespace realredgp;

  RealReductiveGroup& G = currentRealGroup();
  G.fillCartan();

  KGB kgb(G);
  if (checkBasePoint(kgb))
    std::cerr << "true" << std::endl;
  else
    std::cerr << "false" << std::endl;

  return;
}

void cmatrix_f()

/*
  Prints the Cartan matrix on stdout.
*/

{  
  using namespace latticetypes;
  using namespace prettyprint;
  using namespace rootdata;
  using namespace testprint;

  LatticeMatrix q;

  cartanMatrix(q,currentRootDatum());
  printMatrix(std::cout,q);

  return;
}

void components_f()

/*
  Prints the (dual) component group of the current group. We print it out
  in terms of the canonical basis of T(2)^v
*/

{
  using namespace realredgp;
  using namespace latticetypes;
  using namespace ioutils;
  using namespace testprint;

  const RealReductiveGroup& G = realmode::currentRealGroup();
  const ComponentList& c = G.componentReps();
  
  if (c.size() > 0)
    std::cout << "component group is (Z/2)^" << c.size() << std::endl;
  else
    std::cout << "group is connected" << std::endl;

  return;
}

void coroots_rootbasis_f()

/*
  Prints the coroots in the simple coroot coordinates.
*/

{
  using namespace basic_io;
  using namespace lattice;
  using namespace latticetypes;

  const rootdata::RootDatum& rd = currentRootDatum();
  ioutils::OutputFile file;

  std::vector<Weight> v;
  baseChange(rd.beginCoroot(),rd.endCoroot(),back_inserter(v),
	     rd.beginSimpleCoroot(),rd.endSimpleCoroot());
  seqPrint(file,v.begin(),v.end(),"\n","","") << std::endl;

  return;
}

void corder_f()

/*
  Synopsis: prints the Hasse diagram of the ordering of Cartan classes.
*/

{
  using namespace realmode;
  using namespace realredgp;
  using namespace realredgp_io;
  using namespace ioutils;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();
  }
  catch (error::MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  std::cout << "hasse diagram of Cartan ordering:" << std::endl;
  printCartanOrder(std::cout,G_R);

  return;
}

void extrkl_f()

/*!
  \brief Prints out the list of all non-zero k-l polynomials for primitive 
  pairs.

  Explanation: x is primitive w.r.t. y, if any descent for y is also a
  descent for x, or a type II imaginary ascent. Ths means that none of
  the easy recursion formulas applies to P_{x,y}.
*/

{
  using namespace basic_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace kl;
  using namespace kl_io;
  using namespace klsupport;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();
  }
  catch (MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
  const realredgp_io::Interface& G_RI = currentRealInterface();
  const complexredgp_io::Interface& G_I = G_RI.complexInterface();

  // get dual real form
  RealForm drf;

  try {
    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());
  }
  catch (InputError& e) {
    e("aborted");
    return;
  }

  Block block(G_C,G_R.realForm(),drf);

  KLSupport kls(block);
  kls.fill();

  KLContext klc(kls);
  klc.fill();

  OutputFile file;
  file << "Non-zero Kazhdan-Lusztig-Vogan polynomials for primitive pairs:"
       << std::endl << std::endl;
  printExtremalKL(file,klc);

  return;
}

void kgb_f()

/*
  Outputs the kgb table.
*/

{
  using namespace basic_io;
  using namespace kgb;
  using namespace kgb_io;
  using namespace realmode;
  using namespace realredgp;

  RealReductiveGroup& G = currentRealGroup();
  G.fillCartan();

  std::cout << "kgbsize: " << G.kgbSize() << std::endl;
  ioutils::OutputFile file;

  KGB kgb(G);
  printKGB(file,kgb);

  return;
}

void klbasis_f()

/*
  Synopsis: for each element y in the block, outputs the list of non-zero
  k-l polynomials P_{x,y}.

  This is what is required to write down the k-l basis element c_y.
*/

{
  using namespace basic_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace kl;
  using namespace kl_io;
  using namespace klsupport;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();
  }
  catch (MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
  const realredgp_io::Interface& G_RI = currentRealInterface();
  const complexredgp_io::Interface& G_I = G_RI.complexInterface();

  // get dual real form
  RealForm drf;

  try {
    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());
  }
  catch (InputError& e) {
    e("aborted");
    return;
  }

  Block block(G_C,G_R.realForm(),drf);

  KLSupport kls(block);
  kls.fill();

  KLContext klc(kls);
  klc.fill();

  OutputFile file;
  file << "Full list of non-zero Kazhdan-Lusztig-Vogan polynomials:"
       << std::endl << std::endl;
  printAllKL(file,klc);

  return;
}

void kllist_f()

/*
  Synopsis: outputs the list of all distinct Kazhdan-Lusztig-Vogan polynomials
*/

{
  using namespace basic_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace kl;
  using namespace kl_io;
  using namespace klsupport;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();
  }
  catch (MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
  const realredgp_io::Interface& G_RI = currentRealInterface();
  const complexredgp_io::Interface& G_I = G_RI.complexInterface();

  // get dual real form
  RealForm drf;

  try {
    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());
  }
  catch (InputError& e) {
    e("aborted");
    return;
  }

  Block block(G_C,G_R.realForm(),drf);

  KLSupport kls(block);
  kls.fill();

  KLContext klc(kls);
  klc.fill();

  OutputFile file;
  printKLList(file,klc);

  return;
}
  
void poscoroots_rootbasis_f()

/*
  Prints the positive coroots in the simple coroot coordinates.
*/

{
  using namespace basic_io;
  using namespace lattice;
  using namespace latticetypes;

  const rootdata::RootDatum& rd = currentRootDatum();
  ioutils::OutputFile file;

  std::vector<Weight> v;
  baseChange(rd.beginPosCoroot(),rd.endPosCoroot(),back_inserter(v),
	     rd.beginSimpleCoroot(),rd.endSimpleCoroot());
  seqPrint(file,v.begin(),v.end(),"\n","","") << std::endl;

  return;
}
  
void posroots_rootbasis_f()

/*
  Prints the positive roots in the simple root coordinates.
*/

{
  using namespace basic_io;
  using namespace lattice;
  using namespace latticetypes;

  const rootdata::RootDatum& rd = currentRootDatum();
  ioutils::OutputFile file;

  std::vector<Weight> v;
  baseChange(rd.beginPosRoot(),rd.endPosRoot(),back_inserter(v),
	     rd.beginSimpleRoot(),rd.endSimpleRoot());
  seqPrint(file,v.begin(),v.end(),"\n","","") << std::endl;

  return;
}
  
void roots_rootbasis_f()

/*
  Prints the roots in the simple root coordinates.
*/

{
  using namespace basic_io;
  using namespace lattice;
  using namespace latticetypes;

  const rootdata::RootDatum& rd = currentRootDatum();
  ioutils::OutputFile file;

  std::vector<Weight> v;
  baseChange(rd.beginRoot(),rd.endRoot(),back_inserter(v),
	     rd.beginSimpleRoot(),rd.endSimpleRoot());
  seqPrint(file,v.begin(),v.end(),"\n","","") << std::endl;

  return;
}

void rootdatum_f()

/*
  Prints information about the root datum (see testprint.cpp for details).
*/

{
  const rootdata::RootDatum& rd = currentRootDatum();
  ioutils::OutputFile file;

  testprint::print(file,rd);

  return;
}

void wcells_f()

/*
  Synopsis: outputs the cells of the W-graph of the block.
*/

{
  using namespace basic_io;
  using namespace blocks;
  using namespace bruhat;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace kl;
  using namespace klsupport;
  using namespace partition;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;
  using namespace wgraph;
  using namespace wgraph_io;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();
  }
  catch (MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
  const realredgp_io::Interface& G_RI = currentRealInterface();
  const complexredgp_io::Interface& G_I = G_RI.complexInterface();

  // get dual real form
  RealForm drf;

  try {
    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());
  }
  catch (InputError& e) {
    e("aborted");
    return;
  }

  Block block(G_C,G_R.realForm(),drf);

  KLSupport kls(block);
  kls.fill();

  KLContext klc(kls);
  klc.fill();

  WGraph wg(klc.rank());
  kl::wGraph(wg,klc);
  
  OutputFile file;
  printCells(file,wg);

  return;
}

void wgraph_f()

/*
  Synopsis: outputs the W-graph corresponding to a block.
*/

{
  using namespace basic_io;
  using namespace blocks;
  using namespace commands;
  using namespace error;
  using namespace interactive;
  using namespace ioutils;
  using namespace kl;
  using namespace kl_io;
  using namespace klsupport;
  using namespace realform;
  using namespace realmode;
  using namespace realredgp;
  using namespace tags;
  using namespace wgraph;
  using namespace wgraph_io;

  RealReductiveGroup& G_R = currentRealGroup();

  try {
    G_R.fillCartan();
  }
  catch (MemoryOverflow& e) {
    e("error: memory overflow");
    return;
  }

  complexredgp::ComplexReductiveGroup& G_C = G_R.complexGroup();
  const realredgp_io::Interface& G_RI = currentRealInterface();
  const complexredgp_io::Interface& G_I = G_RI.complexInterface();

  // get dual real form
  RealForm drf;

  try {
    getInteractive(drf,G_I,G_C.dualRealFormLabels(G_R.mostSplit()),DualTag());
  }
  catch (InputError& e) {
    e("aborted");
    return;
  }

  Block block(G_C,G_R.realForm(),drf);

  KLSupport kls(block);
  kls.fill();

  KLContext klc(kls);
  klc.fill();

  WGraph wg(klc.rank());
  kl::wGraph(wg,klc);
  
  OutputFile file;
  printWGraph(file,wg);

  return;
}

void test_f()

/*
  Response to the "test" command.
*/

{
  using namespace involutions;
  using namespace mainmode;
  using namespace complexredgp;

  ComplexReductiveGroup& G = currentComplexGroup();

  InvolutionSet inv(G);

  return;
}

}

}
