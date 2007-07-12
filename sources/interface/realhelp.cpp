/*
  This is realhelp.cpp
  
  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups 

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
  void strongreal_h();

  // command tags for the help facility
  const char* cartan_tag = "prints the conjugacy classes of Cartan subgroups";
  const char* realform_tag = "sets the real form for the group";
  const char* realweyl_tag = "outputs the structure of the real Weyl group";
  const char* strongreal_tag = "outputs information about strong real forms";

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
  mode.add("strongreal",strongreal_h);

  insertTag(tagDict,"cartan",cartan_tag);
  insertTag(tagDict,"realform",realform_tag);
  insertTag(tagDict,"realweyl",realweyl_tag);
  insertTag(tagDict,"strongreal",strongreal_tag);

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


void strongreal_h()

{
  io::printFile(std::cerr,"strongreal.help",io::MESSAGE_DIR);
  return;
}

}

}
