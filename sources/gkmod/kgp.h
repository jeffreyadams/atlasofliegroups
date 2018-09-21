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
  bitmap::BitMap members;

public:
  // Constructor
  KGP_orbit(KGBElt size) : members(size) {};

  // print function
  std::ostream& print(std::ostream& strm) const;

  // comparison based on dimension (i.e. length of open orbit)
  bool operator< (const KGP_orbit &elt) const  { return (open() < elt.open()); }

  // accessors
  KGBElt open() const
  { unsigned long last=members.capacity();
    bool success = members.back_up(last); // set |last| to final member
    assert(success); ndebug_use(success);
    return static_cast<KGBElt>(last); // return result as |KGBElt|
  }

  size_t size() const { return members.size(); }

  bitmap::BitMap::iterator begin() const { return members.begin(); }
  bitmap::BitMap::iterator end() const   { return members.end(); }

  // manipulators
  void insert (KGBElt x) { members.insert(x); }

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
  Poset::EltList reduce
  (KGPElt i, bitmap::BitMap closure, std::vector<Poset::EltList>& hasse);

}; // |class KGP|


} // |namespace kgb|
} // |namespace atlas|

#endif



