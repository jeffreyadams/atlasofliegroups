/*
  This is kgb_io.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  part of the Atlas of Reductive Lie Groups

  For license information see the LICENSE file
*/

#include <iomanip>
#include <iostream>

#include "bruhat.h"
#include "kgb_io.h"

#include "ioutils.h"
#include "complexredgp.h"
#include "kgb.h"
#include "prettyprint.h"
#include "set.h"

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
  Print the data from |kgb| to |strm|.

  Explanation: for each parameter, we output the length, the Cartan class,
  root types, cross-actions and Cayley transforms for each generator, and the
  underlying root datum involution (or rather, the corresponding Weyl group
  element). We use a '*' for undefined Cayley transforms.

  NOTE: this will print reasonably on 80 columns only for groups that are
  not too large (up to rank 4 or so). We haven't tried to go over to more
  sophisticated formatting for larger groups.
*/
std::ostream& print(std::ostream& strm,
		    const kgb::KGB& kgb,
		    const complexredgp::ComplexReductiveGroup* G,
		    const kgb::KGBEltList* which)
{
  bool extra= G!=NULL;
  bool subset= which!=NULL;
  bool traditional = not (extra or subset);

  if (not traditional)
    prettyprint::prettyPrint(strm << "Base grading: [",
			     kgb.base_grading(),kgb.rank()) << "].\n";
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
    strm << std::setw(lwidth) << kgb.length(j);
    strm << std::setw(pad) << "";

    if (traditional)  // print Cartan class (in traditional mode)
    {
      strm << std::setw(cwidth) << kgb.Cartan_class(j);
      strm << std::setw(pad) << "";
    }

    // print status
    prettyprint::printStatus(strm,kgb.status(j),kgb.rank());
    strm << ' ';

    // print cross actions
    for (size_t s = 0; s < kgb.rank(); ++s) {
      strm << std::setw(width+pad) << kgb.cross(s,j);
    }
    strm << std::setw(pad) << "";

    // print Cayley transforms
    for (size_t s = 0; s < kgb.rank(); ++s) {
      kgb::KGBElt z = kgb.cayley(s,j);
      if (z != kgb::UndefKGB)
	strm << std::setw(width+pad) << z;
      else
	strm << std::setw(width+pad) << '*';
    }
    strm << std::setw(pad) << "";

    if (not traditional)
    {
      tits::TitsElt a=kgb.titsElt(j);
      const tits::TitsGroup& Tg=kgb.titsGroup();
      Tg.mult(a,kgb.basedTitsGroup().twisted(a));
      assert(a==tits::TitsElt(Tg));

    // print torus part
      prettyprint::prettyPrint(strm,kgb.torus_part(j)) <<
	(extra and kgb.involution(j)==G->twistedInvolution(kgb.Cartan_class(j))
	? '#' : ' ') <<
	std::setw(cwidth) << kgb.Cartan_class(j) << std::setw(pad) << "";
    }
    // print root datum involution
    prettyprint::printWeylElt(strm,kgb.involution(j),kgb.weylGroup());

    strm << std::endl;
  }

  return strm;
}

std::ostream& printKGB(std::ostream& strm, const kgb::KGB& kgb)
{
  return print(strm,kgb,NULL,NULL);
}

std::ostream& print_sub_KGB(std::ostream& strm,
			    const kgb::KGB& kgb,
			    const kgb::KGBEltList& which)
{
  return print(strm,kgb,NULL,&which);
}

std::ostream& var_print_KGB(std::ostream& strm,
			    const complexredgp::ComplexReductiveGroup& G,
			    const kgb::KGB& kgb)
{
  return print(strm,kgb,&G,NULL);
}


std::ostream& print_X(std::ostream& strm, const kgb::global_KGB& kgb)
{
  {
    tits::TorusElement
      yrho(latticetypes::RatWeight(kgb.rootDatum().dual_twoRho(),4));

    strm << "\\exp(i\\pi\\check\\rho) = \\exp(2i\\pi(" << yrho.as_rational()
	 << "))" << std::endl;
  }

  size_t size = kgb.size();
  // compute maximal width of entry
  int width = ioutils::digits(size-1,10ul);
  int cwidth = ioutils::digits(kgb.Cartan_class(size-1),10ul);
  int lwidth = ioutils::digits(kgb.length(size-1),10ul);
  const int pad = 2;

  for (size_t i = 0; i<size; ++i)
  {
    strm << std::setw(width) << i << ":  ";

    // print length
    strm << std::setw(lwidth) << kgb.length(i);
    strm << std::setw(pad) << "";

    // print status
    prettyprint::printStatus(strm,kgb.status(i),kgb.rank());
    strm << ' ';

    // print cross actions
    for (size_t s = 0; s < kgb.rank(); ++s)
      strm << std::setw(width+pad) << kgb.cross(s,i);
    strm << std::setw(pad) << "";

    // print Cayley transforms
    for (size_t s = 0; s < kgb.rank(); ++s)
    {
      kgb::KGBElt C = kgb.cayley(s,i), iC1=kgb.inverseCayley(s,i).first;
      strm << std::setw(width+pad);
      if (C != kgb::UndefKGB) strm << C;
      else if (iC1 != kgb::UndefKGB) strm << iC1;
      else strm << '*';
    }
    strm << std::setw(pad) << "";

    // print torus part
    tits::TorusElement a=kgb.torus_part(i);
    strm << a.as_rational() << ' ' <<
      std::setw(cwidth) << kgb.Cartan_class(i) << std::setw(pad) << "";

    // print root datum involution
    prettyprint::printWeylElt(strm,kgb.involution(i),kgb.weylGroup());

    strm << std::endl;
  }

  return strm;
}


// Print the Hasse diagram of the Bruhat ordering |bruhat| to |strm|.
std::ostream&
printBruhatOrder(std::ostream& strm, const bruhat::BruhatOrder& bruhat)
{
  size_t size = bruhat.size();
  strm << "0:" << std::endl;
  for (size_t j = 1; j < size; ++j) {
    const set::SetEltList& e = bruhat.hasse(j);
    strm << j << ": ";
    basic_io::seqPrint(strm,e.begin(),e.end()) << std::endl;
  }

  unsigned long nc=bruhat.n_comparable();
  strm << "Number of comparable pairs = " << nc << std::endl;

  return strm;
} //printBruhatOrder

} // namespace kgb_io

} // namespace atlas
