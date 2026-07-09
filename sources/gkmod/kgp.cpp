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

#include <cassert>
#include <vector>
#include <bitset>
#include <algorithm>
#include <iomanip>

namespace atlas {
  namespace kgb {

    typedef containers::queue<KGPElt> KGP_queue;

    // Methods of |KGP_orbit|

std::ostream& KGP_orbit::print(std::ostream& strm) const
{
  // print the orbit elements
  for (auto it=members.begin(); it(); ++it)
    strm << (it==members.begin() ? '{' : ',') << *it;
  strm << '}';

  return strm;
}


    // Methods of |KGP|


    // Constructor
KGP::KGP(realredgp::RealReductiveGroup& G_R, RankFlags generators)
: kgb(G_R.kgb())
, kgborder(G_R.Bruhat_KGB())
, kgptable(kgb.size(),KGPElt(~0))
, data()
, bruhat(nullptr)
, msize(0)
{
  // local data
  size_t kgbsize = kgb.size();
  KGPElt unassigned = ~0;

  // determine the KGP orbit for each KGB orbit
  size_t count = 0;
  for (KGBElt i=0; i<kgbsize; i++) // increasing order, necessitated by Cayleys
  { // see if we found a new orbit
    KGP_queue q;
    if (kgptable[i] != unassigned)
      continue; // skip elements previously attributed to an orbit

    auto cur_orbit =  count++; // start a new orbit
    kgptable[i] = cur_orbit;
    q.push(i); // all of whose elements are queued for further processing

    // for each element of the queue, explore each of its edges
    while (not q.empty())
    { // get the orbit
      KGBElt kgbelt = q.front();
      q.pop();

      // for each simple root, determine the other elements
      // look for other elements in the orbit
      for (RankFlags::iterator it=generators.begin(); it(); ++it)
      {
	weyl::Generator j=*it;
	// the root is in the parabolic, check the cross action
	KGBElt ca = kgb.cross(j,kgbelt);

	// see if we found a new element
	if (kgptable[ca] == unassigned)
	{
	  kgptable[ca] = cur_orbit;
	  q.push(ca);
	}
	else
	  assert(kgptable[ca]==cur_orbit); // we don't expect fusion with old orbit

	// if the root is noncompact, also check the Cayley transform
	gradings::Status::Value rt = kgb.status(j,kgbelt);
	if (rt == gradings::Status::ImaginaryNoncompact)
	{ // see if the cayley transform gives a new element
	  KGBElt ct = kgb.cayley(j,kgbelt);
	  if (kgptable[ct] == unassigned)
	  {
	    kgptable[ct] = cur_orbit;
	    q.push(ct);
	  }
	  else
	    assert(kgptable[ca]==cur_orbit); // don't expect fusion with old orbit
	}
      } // |for(it)|

    }  // |while (not q.empty())|

  } // |for(i)|

  // allocate enough memory to hold the orbits
  data.reserve(count);
  for (unsigned i=0; i<count; ++i)
    data.emplace_back(kgbsize);

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
  // use a bitmap to keep track of closure relations
  std::vector<Poset::EltList> hasse(kgp_size);
  for (KGPElt i=0; i<kgp_size; ++i)
  { // for each KGB orbit in this KGP orbit, examine closure edges of degree one
    const KGP_orbit& kgp_orbit = data[i];
    bitmap::BitMap closure(kgp_size);
    for (auto rep : kgp_orbit)
    { // get KGB elements covered by our current |kgp_orbit| representative |rep|
      const Poset::EltList& covered_list = kgborder.hasse(rep);

      // fill the KGP closure list
      for (auto covered : covered_list)
	closure.insert(kgptable[covered]); // add orbit containing covered element

    } // |for (auto rep : kgp_orbit)|

    // filter |closure|, removing redundancies from earlier |hasse|
    hasse[i] = reduce(i, std::move(closure), hasse);

  } // |for (KGPElt i)|

  // store the hasse diagram
  bruhat = new bruhat::BruhatOrder(hasse);
} // |KGP::fillClosure|

// helper function - removes redundant edges from a closure relation
Poset::EltList KGP::reduce
  (KGPElt i, // element for which |closure| was computed
   bitmap::BitMap closure, // elements in topological closure of |i|
   std::vector<Poset::EltList>& hasse) // earlier covered lists
{
  if (closure.empty())
    return Poset::EltList();
  const KGPElt min_covered = *closure.begin(); // lower bound for result
  containers::simple_list<KGPElt> result; // truly covered elements, decreasing

  unsigned long j=i; // the |unsigned long| type is imposed by |BitMap::back_up|
  while (closure.back_up(j)) // closure is filtered while we are descending it
  {
    result.push_front(j);
    KGP_queue q { KGPElt(j) };

    // find and remove from |closure| elements reachable from |j| via |hasse|
    while(not q.empty())
    {
      // get the next element, and remove it from the queue, and from |closure|
      KGPElt cur_elt = q.front();
      q.pop();
      closure.remove(cur_elt); // filter what we encounter (need not be present)

      // propagate work to consider elements reachable downwards from |cur_elt|
      Poset::EltList& covered_list = hasse[cur_elt];
      for (auto covered : covered_list)
	if (covered >= min_covered) // don't bother pushing too small elements
	  q.push(covered);

    } // |while(not q.empty())|
  } // |while closure.back_up(j)|

  return Poset::EltList(result.wbegin(),result.wend()); // convert to vector

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
    const Poset::EltList& covered_list = bruhat->hasse(i);
    size_t lsize = covered_list.size();

    for (size_t j=0; j<lsize; j++)
    {
      strm << covered_list[j];
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
    const Poset::EltList& covered_list = bruhat->hasse(i);
    size_t clsize = covered_list.size();
    for (size_t j=0; j<clsize; j++) // add an edge in the graph
      strm << "v" << i << " -> v" << covered_list[j]
	   << "[color=gray] [arrowhead=none] [style=bold]" << std::endl;
  }

  // write footer
  strm << "}" << std::endl;
} // |KGP::makeDotFile|

  } // |namespace kgb|
} // |namespace atlas|
