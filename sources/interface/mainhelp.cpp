/*
  This is mainhelp.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
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
  void blocksizes_h();
  void coroots_h();
  void poscoroots_h();
  void posroots_h();
  void roots_h();
  void showdualforms_h();
  void showrealforms_h();
  void simplecoroots_h();
  void simpleroots_h();
  void type_h();

  // command tags for the help facility
  const char* blocksizes_tag = "outputs the matrix of blocksizes";
  const char* coroots_tag = "outputs the coroots in the lattice basis";
  const char* poscoroots_tag = 
    "outputs the positive coroots in the lattice basis";
  const char* posroots_tag = "outputs the positive roots in the lattice basis";
  const char* roots_tag = "outputs the roots in the lattice basis";
  const char* showdualforms_tag = 
    "outputs the weak real forms for the dual group";
  const char* showrealforms_tag = 
    "outputs the weak real forms for this complex group";
  const char* simplecoroots_tag = 
    "outputs the simple coroots in the lattice basis";
  const char* simpleroots_tag = 
    "outputs the simple roots in the lattice basis";
  const char* type_tag = "resets the group type";

}

/*****************************************************************************

        Chapter I --- Functions declared in mainhelp.h

******************************************************************************/

namespace mainhelp {

void addMainHelp(commands::CommandMode& mode, commands::TagDict& tagDict)

{  
  using namespace commands;

  mode.add("blocksizes",blocksizes_h);
  mode.add("coroots",coroots_h);
  mode.add("poscoroots",poscoroots_h);
  mode.add("posroots",posroots_h);
  mode.add("roots",roots_h);
  mode.add("showdualforms",showdualforms_h);
  mode.add("showrealforms",showrealforms_h);
  mode.add("simpleroots",simpleroots_h);
  mode.add("simplecoroots",simplecoroots_h);
  mode.add("type",type_h);

  TagDict::value_type blocksizes_entry = 
    std::make_pair("blocksizes",blocksizes_tag);
  tagDict.insert(blocksizes_entry);

  TagDict::value_type coroots_entry = std::make_pair("coroots",coroots_tag);
  tagDict.insert(coroots_entry);

  TagDict::value_type poscoroots_entry = 
    std::make_pair("poscoroots",poscoroots_tag);
  tagDict.insert(poscoroots_entry);

  TagDict::value_type posroots_entry = std::make_pair("posroots",posroots_tag);
  tagDict.insert(posroots_entry);

  TagDict::value_type roots_entry = std::make_pair("roots",roots_tag);
  tagDict.insert(roots_entry);

  TagDict::value_type showdualforms_entry 
    = std::make_pair("showdualforms",showdualforms_tag);
  tagDict.insert(showdualforms_entry);

  TagDict::value_type showrealforms_entry 
    = std::make_pair("showrealforms",showrealforms_tag);
  tagDict.insert(showrealforms_entry);

  TagDict::value_type simpleroots_entry 
    = std::make_pair("simpleroots",simpleroots_tag);
  tagDict.insert(simpleroots_entry);

  TagDict::value_type simplecoroots_entry 
    = std::make_pair("simplecoroots",simplecoroots_tag);
  tagDict.insert(simplecoroots_entry);

  TagDict::value_type type_entry = std::make_pair("type",type_tag);
  tagDict.insert(type_entry);

  return;
}

}

/*****************************************************************************

        Chapter II --- Functions for the help commands

  This section contains the definitions of the help functions associated to 
  the various commands defined in this mode.

******************************************************************************/

namespace {

void blocksizes_h()

{  
  io::printFile(std::cerr,"blocksizes.help",io::MESSAGE_DIR);
  return;
}

void coroots_h()

{  
  io::printFile(std::cerr,"coroots.help",io::MESSAGE_DIR);
  return;
}

void posroots_h()

{  
  io::printFile(std::cerr,"posroots.help",io::MESSAGE_DIR);
  return;
}

void poscoroots_h()

{  
  io::printFile(std::cerr,"poscoroots.help",io::MESSAGE_DIR);
  return;
}
   
void roots_h()

{  
  io::printFile(std::cerr,"roots.help",io::MESSAGE_DIR);
  return;
}

void showdualforms_h()
  
{  
  io::printFile(std::cerr,"showdualforms.help",io::MESSAGE_DIR);
  return;
}

void showrealforms_h()
  
{  
  io::printFile(std::cerr,"showrealforms.help",io::MESSAGE_DIR);
  return;
}

void simpleroots_h()

{  
  io::printFile(std::cerr,"simpleroots.help",io::MESSAGE_DIR);
  return;
}
  
void simplecoroots_h()

{  
  io::printFile(std::cerr,"simplecoroots.help",io::MESSAGE_DIR);
  return;
}

void type_h()

{
  io::printFile(std::cerr,"type.help",io::MESSAGE_DIR);
  return;
}

}

}
