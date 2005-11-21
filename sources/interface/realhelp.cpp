/*
  This is realhelp.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups version 0.2.4 

  See file main.cpp for full copyright notice
*/

#include "realhelp.h"

#include <iostream>

#include "io.h"
#include "realmode.h"
#include "special.h"
#include "test.h"

/****************************************************************************

  This file contains the help mode for the "real" command mode.

*****************************************************************************/

namespace atlas {

namespace {

  // help commands
  void cartan_h();
  void realform_h();
  void realweyl_h();

  // command tags for the help facility
  const char* cartan_tag = "prints the conjugacy classes of Cartan subgroups";
  const char* realform_tag = "sets the real form for the group";
  const char* realweyl_tag = "outputs the structure of the real Weyl group";

}

/****************************************************************************

        Chapter I --- Functions declared in realhelp.h

*****************************************************************************/

namespace realhelp {

void addRealHelp(commands::CommandMode& mode, commands::TagDict& tagDict)

{  
  using namespace commands;

  mode.add("cartan",&cartan_h);
  mode.add("realform",realform_h);
  mode.add("realweyl",realweyl_h);

  TagDict::value_type cartan_entry = std::make_pair("cartan",cartan_tag);
  tagDict.insert(cartan_entry);

  TagDict::value_type realform_entry = std::make_pair("realform",realform_tag);
  tagDict.insert(realform_entry);

  TagDict::value_type realweyl_entry = std::make_pair("realweyl",realweyl_tag);
  tagDict.insert(realweyl_entry);
  return;
}

}

/*****************************************************************************

        Chapter II --- Functions for the help commands

  This section contains the definitions of the help functions associated to 
  the various commands defined in this mode.

******************************************************************************/

namespace {

void cartan_h()

{
  io::printFile(std::cerr,"cartan.help",io::MESSAGE_DIR);
  return;
}

void realform_h()

{
  io::printFile(std::cerr,"realform.help",io::MESSAGE_DIR);
  return;
}

void realweyl_h()

{
  io::printFile(std::cerr,"realweyl.help",io::MESSAGE_DIR);
  return;
}

}

}
