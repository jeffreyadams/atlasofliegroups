/*
  This is kgb_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "kgb_io.h"

#include <iomanip>
#include <iostream>
#include <set>

#include "complexredgp.h"
#include "kgb.h"

#include "ioutils.h"	// |digits|
#include "basic_io.h"	// |operator<<|
#include "prettyprint.h"// |printStatus|

#include "bruhat.h"	// |BruhatOrder|

/*****************************************************************************

  Input/output functions for the kgb data structure, defined in
  sources/kl/kgb.h

******************************************************************************/

namespace atlas {

/*****************************************************************************

        Chapter I -- Functions declared in kgb_io.h

******************************************************************************/

namespace kgb_io {


/*
  Basic  print of |KGB_base| object |kgb| to |strm|.

  Explanation: for each parameter, we output the length, the Cartan class,
  root types, cross-actions and Cayley transforms for each generator, and the
  underlying root datum involution (or rather, the corresponding Weyl group
  element). We use a '*' for undefined Cayley transforms.

  NOTE: this will print reasonably on 80 columns only for groups that are
  not too large (up to rank 4 or so). We haven't tried to go over to more
  sophisticated formatting for larger groups.
*/
std::ostream& print(std::ostream& strm,
		    const KGB_base& kgb,
		    bool traditional,
		    const ComplexReductiveGroup* G,
		    const KGBEltList* which)
{
  bool subset= which!=NULL;

  // compute maximal width of entry
  int width = ioutils::digits(kgb.size()-1,10ul);
  int cwidth = ioutils::digits(kgb.Cartan_class(kgb.size()-1),10ul);
  int lwidth = ioutils::digits(kgb.length(kgb.size()-1),10ul);
  const int pad = 2;

  size_t size= subset ? which->size() : kgb.size();
  for (size_t i = 0; i<size; ++i)
  {
    size_t j = subset ? (*which)[i] : i;
    strm << std::setw(width) << j << ":  ";

    // print length
    strm << std::setw(lwidth) << kgb.length(j) << std::setw(pad) << "";

    if (traditional)  // print Cartan class here in traditional mode
    {
      strm << std::setw(cwidth) << kgb.Cartan_class(j)
	   << std::setw(pad) << "";
    }

    // print status
    prettyprint::printStatus(strm,kgb.status(j),kgb.rank()) << ' ';

    // print cross actions
    for (size_t s = 0; s < kgb.rank(); ++s) {
      strm << std::setw(width+pad) << kgb.cross(s,j);
    }
    strm << std::setw(pad) << "";

    // print Cayley transforms
    for (size_t s = 0; s < kgb.rank(); ++s)
    {
      KGBElt z = kgb.cayley(s,j);
      if (z != UndefKGB)
	strm << std::setw(width+pad) << z;
      else
	strm << std::setw(width+pad) << '*';
    }

    strm << std::setw(pad) << "";
    if (not traditional)
    {
      kgb.print(strm,j) // virtual print method
	<< (G!=NULL and
	    kgb.involution(j)==G->twistedInvolution(kgb.Cartan_class(j))
	    ? '#' : ' ');
      unsigned cc= kgb.Cartan_class(j);
      if (cc!=~0u)
	strm << std::setw(cwidth) << cc << ' ';
    }

   // print root datum involution
    prettyprint::printInvolution(strm,kgb.involution(j),kgb.twistedWeylGroup())
      << std::endl;
  }

  return strm;
}


std::ostream& printKGB(std::ostream& strm, const KGB& kgb)
{
  return print(strm,kgb,true,NULL,NULL);
}

std::ostream& print_sub_KGB(std::ostream& strm,
			    const KGB& kgb,
			    const KGBEltList& which)
{
  prettyprint::prettyPrint(strm << "Base grading: [",
			   kgb.base_grading(),kgb.rank()) << "].\n";
  return print(strm,kgb,false,NULL,&which);
}

std::ostream& var_print_KGB(std::ostream& strm,
			    const ComplexReductiveGroup& G,
			    const KGB& kgb)
{
  prettyprint::prettyPrint(strm << "Base grading: [",
			   kgb.base_grading(),kgb.rank()) << "].\n";
  return print(strm,kgb,false,&G,NULL);
}


std::ostream& print_X(std::ostream& strm, const kgb::global_KGB& kgb)
{
  {
    TorusElement yrho =
      y_values::exp_2pi(kgb.globalTitsGroup().torus_part_offset());

    strm << "\\exp(i\\pi\\check\\rho) = \\exp(2i\\pi("
	 << yrho.log_2pi() << "))" << std::endl;
  }
  return print(strm,kgb,false,NULL,NULL);
}


// Print the Hasse diagram of the Bruhat ordering |bruhat| to |strm|.
std::ostream&
printBruhatOrder(std::ostream& strm, const BruhatOrder& bruhat)
{
  size_t size = bruhat.size();
  strm << "0:" << std::endl;
  for (size_t j = 1; j < size; ++j) {
    const set::EltList& e = bruhat.hasse(j);
    strm << j << ": ";
    basic_io::seqPrint(strm,e.begin(),e.end()) << std::endl;
  }

  unsigned long nc=bruhat.n_comparable();
  strm << "Number of comparable pairs = " << nc << std::endl;

  return strm;
} //printBruhatOrder

// make a '.dot' file that can be processed by the 'dot' program
// see www.graphviz.org for more info
void makeDotFile(std::ostream& strm, const KGB& kgb, const BruhatOrder& bruhat) {
  // local data
  size_t size = kgb.size();
  size_t rank = kgb.rank();

  // edge colors - feel free to change these
  // but remember to change them in the help file as well - spc
  std::string colorca("black"); // cross action 
  std::string colorctI("blue"); // cayley transform type I
  std::string colorctII("green"); // cayley transform type II
  std::string colorcl("gray"); // additional R/S closure edge

  // write header
  strm << "digraph G {" << std::endl << "ratio=\"1.5\"" << std::endl << "size=\"7.5,10.0\"" << std::endl;

  // create the vertices
  for (size_t i=0; i<size; i++) {
    strm << "v" << i << std::endl;    
  }

  // vector of sets to track which closure edges come from
  // cross actions and cayley transforms
  std::vector<std::set<size_t> > edges(size);

  // build the c/a and c/t graph
  for (size_t i=0; i<size; i++) {
    for (size_t j=0; j<rank; j++) {
		  // depending on the type of root, add an edge if appropriate
      const gradings::Status& type = kgb.status(i);
      
      // CASE: complex cross action
      if (type[j] == gradings::Status::Complex) {
        // get the cross action
        KGBElt ca = kgb.cross(j,i);

        // its only a closure edge if it inceases length
        if (ca > i) {
          // add an edge in the graph
          strm << "v" << ca << " -> v" << i << "[color=" << colorca << "] [arrowhead=none] [style=bold]" << std::endl;
          edges[ca].insert(i);
        }
      }

      // CASE: cayley transform
      else if (type[j] == gradings::Status::ImaginaryNoncompact) {
        // get the cross action and cayley transform
        KGBElt ca = kgb.cross(j,i);
        KGBElt ct = kgb.cayley(j,i);

        // CASE: type I
        if (ca != i) {
          // add an edge in the graph
          strm << "v" << ct << " -> v" << i << "[color=" << colorctI << "] [arrowhead=none] [style=bold]" << std::endl;
          edges[ct].insert(i);
        }

        // CASE: type II
        else {
          // add an edge in the graph
          strm << "v" << ct << " -> v" << i << "[color=" << colorctII << "] [arrowhead=none] [style=bold]" << std::endl;
          edges[ct].insert(i);
        }
      }
    }
  }

  // finally, add the closure edges
  for (size_t i=0; i<size; i++) {
    // get the list
    const set::EltList& clist = bruhat.hasse(i);
    size_t clsize = clist.size();

    // add edges for the ones that arent already there
    for (size_t j=0; j<clsize; j++) {
      set::Elt e = clist[j];
      if (edges[i].count(e) == 0) {
          // add an edge in the graph
          strm << "v" << i << " -> v" << e << "[color=" << colorcl << "] [arrowhead=none]" << std::endl;        
      }
    }
  }

  // write footer
  strm << "}" << std::endl;
}

} // namespace kgb_io

} // namespace atlas
