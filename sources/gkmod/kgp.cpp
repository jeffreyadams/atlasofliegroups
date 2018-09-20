/*
  This is kgp.cpp

  Copyright (C) 2011 Scott Crofts
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "kgp.h"

#include "constants.h"
#include "gradings.h"
#include "realredgp.h"
#include "ioutils.h"

#include <vector>
#include <bitset>
#include <algorithm>
#include <iomanip>

namespace atlas {
  namespace kgb {

    // Methods of |KGP_orbit|

std::ostream& KGP_orbit::print(std::ostream& strm) const
{
  size_t size = members.size();

  // print the orbit elements
  strm << "{";
  for (size_t i=0; i<size-1; i++)
    strm << members[i] << ",";
  strm << members[size-1] << "}";

  return strm;
}


    // Methods of |KGP|


    // Constructor
KGP::KGP(realredgp::RealReductiveGroup& G_R, RankFlags generators)
: kgb(G_R.kgb())
, kgborder(G_R.Bruhat_KGB())
, kgptable(kgb.size(),KGPElt(~0))
, bruhat(nullptr)
, msize(0)
{
  // local data
  size_t kgbsize = kgb.size();
  KGPElt unassigned = ~0;

  // determine the KGP orbit for each KGB orbit
  size_t count = 0;
  for (KGBElt i=0; i<kgbsize; i++)
  { // see if we found a new orbit
    KGP_queue q;
    if (kgptable[i] == unassigned)
    { // new orbit
      kgptable[i] = count++;
      q.push(i);
    }

    // for each element of the queue, explore each of its edges
    while (not q.empty())
    { // get the orbit
      KGBElt kgbelt = q.front();

      // for each simple root, determine the other elements
      // look for other elements in the orbit
      for (RankFlags::iterator it=generators.begin(); it(); ++it)
      {
	weyl::Generator j=*it;
	// the root is in the parabolic, check the cross action
	KGBElt ca = kgb.cross(j,kgbelt);

	// see if we found a new element
	if (kgptable[ca] == unassigned)
	{ // new element - add it to the queue
	  kgptable[ca] = kgptable[kgbelt];
	  q.push(ca);
	}

	// if the root is noncompact, also check the Cayley transform
	gradings::Status::Value rt = kgb.status(j,kgbelt);
	if (rt == gradings::Status::ImaginaryNoncompact)
	{ // see if the cayley transform gives a new element
	  KGBElt ct = kgb.cayley(j,kgbelt);
	  if (kgptable[ct] == unassigned)
	  { // new element - add it to the queue
	    kgptable[ct] = kgptable[kgbelt];
	    q.push(ct);
	  }
	}
      } // |for(it)|

      // at this point we have processed the orbit
      // so remove it from the queue
      q.pop();
    }  // |while (not q.empty())|
  } // |for(i)|

  // allocate enough memory to hold the orbits
  data.resize(count);

  // fill the orbits
  for (KGBElt i=0; i<kgbsize; i++)
  {
    data[kgptable[i]].insert(i);
    if (data[kgptable[i]].size() > msize)
      msize = data[kgptable[i]].size();
  }

  // sort by dimension - i.e. length of the open orbit
  std::sort(data.begin(), data.end());

  // the above sort undoubtedly messed up the mapping, so we need to fix it
  for (KGPElt i=0; i<count; i++)
    for (auto elt : data[i])
      kgptable[elt] = i;

}  // |KGP::KGP|


    // fill closure function
void KGP::fillClosure()
{
  if (bruhat!= nullptr)
    return;

  size_t kgp_size = data.size();

  // build the Hasse diagram
  // use a bit vector to keep track of closure relations
  std::vector<Poset::EltList> hasse(kgp_size);
  std::vector<bool> closure(kgp_size,0);
  for (KGPElt i=0; i<kgp_size; ++i)
  { // for each KGB orbit in this KGP orbit, examine closure edges of degree one
    KGP_orbit kgp_orbit = data[i];
    size_t min_elt = kgp_size;
    for (auto elt : kgp_orbit)
    { // get closure edges
      const Poset::EltList& clist = kgborder.hasse(elt);
      size_t lsize = clist.size();

      // fill the KGP closure list
      for (size_t k=0; k<lsize; ++k)
      {
	KGPElt cur_orbit = kgptable[clist[k]];
	closure[cur_orbit]=1;
	if (cur_orbit < min_elt)
	  min_elt = cur_orbit;
      }
    }

    // at this point, closure contains a generating set of edges for
    // the closure relation. We now reduce this set to a minimal
    // generating set
    KGP_queue q;
    for (KGPElt j=i; j-->min_elt;)
      if (closure[j]==1)
      {
	hasse[i].push_back(j);
	q.push(j);
	reduce(q, closure, hasse, min_elt);
      }


    closure[i]=0;
    std::sort(hasse[i].begin(), hasse[i].end());
  }

  // store the hasse diagram
  bruhat = new bruhat::BruhatOrder(hasse);
} // |KGP::fillClosure|

// helper function - removes redundant edges from a closure relation
void KGP::reduce(KGP_queue& q,
		 std::vector<bool>& closure,
		 std::vector<Poset::EltList>& hasse,
		 KGPElt min_elt)
{
  // while the queue is not empty, recursively remove edges
  while(!q.empty())
  {
    // get the next element
    KGPElt cur_elt = q.front();

    // remove it from the list
    closure[cur_elt]=0;

    // walk the list of lower edges
    Poset::EltList& clist = hasse[cur_elt];
    size_t lsize = clist.size();
    for (size_t i=0; i<lsize; i++)
      if (clist[i] >= min_elt)
	q.push(clist[i]);

    // remove the element
    q.pop();
  }
} // |KGP::reduce|

// print functions
std::ostream& KGP::print(std::ostream& strm) const
{
  size_t kgb_size = kgb.size();
  size_t kgp_size = data.size();
  size_t kgbwidth = ioutils::digits(kgb_size-1,10);
  size_t kgpwidth = ioutils::digits(kgp_size-1,10);
  size_t lwidth = ioutils::digits(kgb.length(kgb_size-1)-1,10);
  size_t cwidth = ioutils::digits(msize,10);

  // print orbits
  for (size_t i=0; i<kgp_size; i++)
  {
    strm << std::setw(kgpwidth) << i << ":[" << std::setw(kgbwidth)
	 << data[i].open() << "*] ";
    strm << std::setw(lwidth) << kgb.length(data[i].open()) << " "
	 << std::setw(cwidth) << data[i].size() << " ";
    data[i].print(strm);
    strm << std::endl;
  }

  return strm;
} // |KGP::print|

std::ostream& KGP::printClosure(std::ostream& strm) const
{
  size_t kgp_size = data.size();
  size_t kgpwidth = ioutils::digits(kgp_size-1,10);

  for (KGPElt i=0; i<kgp_size; i++)
  {
    // print the orbit
    strm << std::setw(kgpwidth) << i << ": ";

    // print the list
    const Poset::EltList& clist = bruhat->hasse(i);
    size_t lsize = clist.size();

    for (size_t j=0; j<lsize; j++)
    {
      strm << clist[j];
      if (j!=lsize-1) strm << ",";
    }
    strm << std::endl;
  }

  return strm;
}

// make a '.dot' file that can be processed by the 'dot' program
// see www.graphviz.org for more info
void KGP::makeDotFile(std::ostream& strm)
{
  // make sure the closure order has been computed
  fillClosure();

  size_t kgp_size = data.size();

  // write header
  strm << "digraph G {" << std::endl << "ratio=\"1.5\"" << std::endl
       << "size=\"7.5,10.0\"" << std::endl;

  // create vertices
  for (size_t i=0; i<kgp_size; i++)
    strm << "v" << i << std::endl; // create the vertex

  // add edges
  for (size_t i=0; i<kgp_size; i++)
  {
    const Poset::EltList& clist = bruhat->hasse(i);
    size_t clsize = clist.size();
    for (size_t j=0; j<clsize; j++) // add an edge in the graph
      strm << "v" << i << " -> v" << clist[j]
	   << "[color=gray] [arrowhead=none] [style=bold]" << std::endl;
  }

  // write footer
  strm << "}" << std::endl;
} // |KGP::makeDotFile|

  } // |namespace kgb|
} // |namespace atlas|
