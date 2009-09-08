/*
  This is mainhelp.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For copyright and license information see the LICENSE file
*/

#include "mainhelp.h"

#include <iostream>

#include "io.h"
#include "mainmode.h"
#include "special.h"
#include "test.h"

/****************************************************************************

  This file contains the help mode for the "main" command mode, which
  is the startup mode for the program.

*****************************************************************************/

namespace atlas {

namespace {

  // help commands
  void type_h();
  void cmatrix_h();
  void rootdatum_h();
  void roots_h();
  void coroots_h();
  void simpleroots_h();
  void simplecoroots_h();
  void posroots_h();
  void poscoroots_h();
  void realform_h();
  void showdualforms_h();
  void showrealforms_h();
  void blocksizes_h();
  void gradings_h();
  void strongreal_h();
  void dualkgb_h();

  // command tags for the help facility
  const char* type_tag = "resets the group type";
  const char* cmatrix_tag = "prints the Cartan matrix";
  const char* rootdatum_tag = "outputs the root datum";
  const char* roots_tag = "outputs the roots in the lattice basis";
  const char* coroots_tag = "outputs the coroots in the lattice basis";
  const char* simpleroots_tag =
    "outputs the simple roots in the lattice basis";
  const char* simplecoroots_tag =
    "outputs the simple coroots in the lattice basis";
  const char* posroots_tag = "outputs the positive roots in the lattice basis";
  const char* poscoroots_tag =
    "outputs the positive coroots in the lattice basis";
  const char* realform_tag = "sets the real form for the group";
  const char* showdualforms_tag =
    "outputs the weak real forms for the dual group";
  const char* showrealforms_tag =
    "outputs the weak real forms for this complex group";
  const char* blocksizes_tag = "outputs the matrix of blocksizes";
  const char* gradings_tag="prints gradings of imaginary roots for real forms";
  const char* strongreal_tag = "outputs information about strong real forms";
  const char* dual_kgb_tag = "prints the KGB data for a dual real form";
}

/*****************************************************************************

        Chapter I --- Functions declared in mainhelp.h

******************************************************************************/

namespace mainhelp {

void addMainHelp(commands::CommandMode& mode, commands::TagDict& tagDict)

{
  using namespace commands;

  mode.add("type",type_h);
  mode.add("cmatrix",cmatrix_h);
  mode.add("rootdatum",rootdatum_h);
  mode.add("roots",roots_h);
  mode.add("coroots",coroots_h);
  mode.add("simpleroots",simpleroots_h);
  mode.add("simplecoroots",simplecoroots_h);
  mode.add("posroots",posroots_h);
  mode.add("poscoroots",poscoroots_h);
  mode.add("realform",realform_h);
  mode.add("showrealforms",showrealforms_h);
  mode.add("showdualforms",showdualforms_h);
  mode.add("blocksizes",blocksizes_h);
  mode.add("gradings",gradings_h);
  mode.add("strongreal",strongreal_h);
  mode.add("dualkgb",dualkgb_h);

  insertTag(tagDict,"type",type_tag);
  insertTag(tagDict,"cmatrix",cmatrix_tag);
  insertTag(tagDict,"rootdatum",rootdatum_tag);
  insertTag(tagDict,"roots",roots_tag);
  insertTag(tagDict,"coroots",coroots_tag);
  insertTag(tagDict,"simpleroots",simpleroots_tag);
  insertTag(tagDict,"simplecoroots",simplecoroots_tag);
  insertTag(tagDict,"poscoroots",poscoroots_tag);
  insertTag(tagDict,"posroots",posroots_tag);
  insertTag(tagDict,"realform",realform_tag);
  insertTag(tagDict,"showdualforms",showdualforms_tag);
  insertTag(tagDict,"showrealforms",showrealforms_tag);
  insertTag(tagDict,"blocksizes",blocksizes_tag);
  insertTag(tagDict,"gradings",gradings_tag);
  insertTag(tagDict,"strongreal",strongreal_tag);
  insertTag(tagDict,"dualkgb",dual_kgb_tag);
}

}

/*****************************************************************************

        Chapter II --- Functions for the help commands

  This section contains the definitions of the help functions associated to
  the various commands defined in this mode.

******************************************************************************/

namespace {

void type_h()
{
  io::printFile(std::cerr,"type.help",io::MESSAGE_DIR);
  return;
}

void cmatrix_h()
{
  io::printFile(std::cerr,"cmatrix.help",io::MESSAGE_DIR);
}

void rootdatum_h()
{
  io::printFile(std::cerr,"rootdatum.help",io::MESSAGE_DIR);
}

void roots_h()
{
  io::printFile(std::cerr,"roots.help",io::MESSAGE_DIR);
}

void coroots_h()
{
  io::printFile(std::cerr,"coroots.help",io::MESSAGE_DIR);
}

void simpleroots_h()
{
  io::printFile(std::cerr,"simpleroots.help",io::MESSAGE_DIR);
}

void simplecoroots_h()
{
  io::printFile(std::cerr,"simplecoroots.help",io::MESSAGE_DIR);
}

void posroots_h()
{
  io::printFile(std::cerr,"posroots.help",io::MESSAGE_DIR);
}

void poscoroots_h()
{
  io::printFile(std::cerr,"poscoroots.help",io::MESSAGE_DIR);
}

void realform_h()
{
  io::printFile(std::cerr,"realform.help",io::MESSAGE_DIR);
}

void showrealforms_h()
{
  io::printFile(std::cerr,"showrealforms.help",io::MESSAGE_DIR);
}

void showdualforms_h()
{
  io::printFile(std::cerr,"showdualforms.help",io::MESSAGE_DIR);
}

void blocksizes_h()
{
  io::printFile(std::cerr,"blocksizes.help",io::MESSAGE_DIR);
}

void gradings_h()
{
  io::printFile(std::cerr,"gradings.help",io::MESSAGE_DIR);
}

void strongreal_h()
{
  io::printFile(std::cerr,"strongreal.help",io::MESSAGE_DIR);
}

void dualkgb_h()
{
  io::printFile(std::cerr,"dualkgb.help",io::MESSAGE_DIR);
}

} // |namespace|

} // |namespace atlas|
