/*
  This is kgp.h

  Copyright (C) 2011 Scott Crofts
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#ifndef KGP_H
#define KGP_H

#include "kgb.h"
#include "bruhat.h"
#include "bitmap.h"
#include "sl_list.h"

namespace atlas {
namespace kgb {

class KGP_orbit
{
  std::vector<KGBElt> members;

public:
  // Constructor
  KGP_orbit() {};

  // print function
  std::ostream& print(std::ostream& strm) const;

  // comparison based on dimension (i.e. length of open orbit)
  bool operator< (const KGP_orbit &elt) const  { return (open() < elt.open()); }

  // accessors
  KGBElt open() const { return members.back(); }
  size_t size() const { return members.size(); }

  std::vector<KGBElt>::const_iterator begin() const { return members.begin(); }
  std::vector<KGBElt>::const_iterator end() const   { return members.end(); }

  // manipulators
  void insert (KGBElt x) { members.push_back(x); }

};

class KGP
{
  // link to KGB graph
  const KGB& kgb;

  // link to the bruhat order on kgb
  const bruhat::BruhatOrder& kgborder;

  // table to hold KGP orbit numbers for each KGB orbit
  std::vector<KGPElt> kgptable;

  // list of KGP orbits
  std::vector<KGP_orbit> data;

  // hasse diagram for closure relations
  // (pointer since there is no default constructor)
  bruhat::BruhatOrder *bruhat;

  // maximum orbit size - needed for print formatting
  size_t msize;

public:
  // Constructor
  KGP(realredgp::RealReductiveGroup& G_R, RankFlags generators);

  // Destructor (free memory if necessary)
  ~KGP() { if (bruhat != NULL) delete bruhat; }

  // compute (strong) closure order on the KGP orbits
  void fillClosure();

  // print functions
  std::ostream& print(std::ostream& strm) const;
  std::ostream& printClosure(std::ostream& strm) const;

  // make a '.dot' file that can be processed by the 'dot' program
  // see www.graphviz.org for more info
  // (not const since could cause the closure order to be filled)
  void makeDotFile(std::ostream& strm);

  // accessors
  size_t size() const { return data.size(); }

private:
  // helper function - removes redundant edges from a closure relation
  typedef containers::queue<KGPElt> KGP_queue;
  void reduce(KGP_queue& q, std::vector<bool>& closure,
	      std::vector<Poset::EltList>& hasse, KGPElt minelt);

}; // |class KGP|


} // |namespace kgb|
} // |namespace atlas|

#endif



